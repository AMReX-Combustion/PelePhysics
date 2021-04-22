#include <reactor.h> 
#include <AMReX_ParmParse.H>
#include "mechanism.H"
#include <PelePhysics.H>
#include <AMReX_Gpu.H>
#include <AMReX_SUNMemory.H>
#include <AMREX_misc.H>

using namespace amrex;
using namespace amrex::sundials;

int eint_rho = 1; // in/out = rhoE/rhoY
int enth_rho = 2; // in/out = rhoH/rhoY 
int use_erkstep=0;
Real relTol    = 1.0e-6;
Real absTol    = 1.0e-10;
Array<Real,NUM_SPECIES+1> typVals = {-1};

/******************************************************************************************/
/* Initialization routine, called once at the begining of the problem */
int reactor_info(int reactor_type, int Ncells)
{
    ParmParse pp("ode");
    pp.query("use_erkstep", use_erkstep);
    pp.query("rtol",relTol);
    pp.query("atol",absTol);
    if(use_erkstep==1)
    {
        Print()<<"Using ERK Step on GPU\n";
    }
    else
    {
        Print()<<"Using ARK Step on GPU\n";
    }
    Print() << "Setting ARK/ERKODE tolerances rtol = " << relTol 
                << " atol = " << absTol << " in PelePhysics \n";
    return(0);
}
/******************************************************************************************/
void SetTypValsODE(const std::vector<double>& ExtTypVals) 
{
    Vector<std::string> kname;
    pele::physics::eos::speciesNames(kname);

    Print() << "Set the typVals in PelePhysics: \n  ";
    int size_ETV = ExtTypVals.size();
    AMREX_ASSERT(size_ETV == typVals.size());
    for (int i=0; i<size_ETV-1; i++) {
        typVals[i] = ExtTypVals[i];
        Print() << kname[i] << ":" << typVals[i] << "  ";
    }
    typVals[size_ETV-1] = ExtTypVals[size_ETV-1];
    Print() << "Temp:"<< typVals[size_ETV-1] <<  " \n";
}
/******************************************************************************************/
/*void SetTolFactODE(double relative_tol,double absolute_tol) 
{
    relTol = relative_tol;
    absTol = absolute_tol;
    Print() << "Set RTOL, ATOL = "<<relTol<< " "<<absTol<<  " in PelePhysics\n";
}*/
/******************************************************************************************/
/* react call with array4 of data */
int react(const amrex::Box& box,
          amrex::Array4<amrex::Real> const& rY_in,
          amrex::Array4<amrex::Real> const& rY_src_in, 
          amrex::Array4<amrex::Real> const& T_in,
          amrex::Array4<amrex::Real> const& rEner_in,  
          amrex::Array4<amrex::Real> const& rEner_src_in,
          amrex::Array4<amrex::Real> const& FC_in,
          amrex::Array4<int> const& mask,
          amrex::Real &dt_react,
          amrex::Real &time,
          const int &reactor_type,
          amrex::gpuStream_t stream) 
{
    int NCELLS, NEQ, neq_tot,flag;
    realtype time_init, time_out;

    void *arkode_mem    = NULL;
    N_Vector y          = NULL;

    NEQ            = NUM_SPECIES+1;
    NCELLS         = box.numPts();
    neq_tot        = NEQ * NCELLS;
    AMREX_ASSERT(NCELLS < std::numeric_limits<int>::max());

    /* User data */
    UserData user_data;
    BL_PROFILE_VAR("AllocsInARKODE", AllocsARKODE);
    user_data = (ARKODEUserData *) The_Arena()->alloc(sizeof(struct ARKODEUserData));  
    BL_PROFILE_VAR_STOP(AllocsARKODE);
    user_data->ncells_d             = NCELLS;
    user_data->neqs_per_cell        = NEQ;
    user_data->ireactor_type        = reactor_type; 
    user_data->iverbose             = 1;
    user_data->stream               = stream;
    user_data->nbBlocks             = std::max(1,NCELLS/32);
    user_data->nbThreads            = 32;

    /* Definition of main vector */
#if defined(AMREX_USE_CUDA)
    y = N_VNewWithMemHelp_Cuda(neq_tot, /*use_managed_mem=*/true, *The_SUNMemory_Helper());
    if(check_flag((void*)y, "N_VNewWithMemHelp_Cuda", 0)) return(1);
#elif defined(AMREX_USE_HIP)
    y = N_VNewWithMemHelp_Hip(neq_tot, /*use_managed_mem=*/true, *The_SUNMemory_Helper());
    if(check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0)) return(1);
#endif

    /* Use a non-default stream for kernel execution */
#if defined(AMREX_USE_CUDA)
    SUNCudaExecPolicy* stream_exec_policy = new SUNCudaThreadDirectExecPolicy(256, stream);
    SUNCudaExecPolicy* reduce_exec_policy = new SUNCudaBlockReduceExecPolicy(256, 0, stream);
    N_VSetKernelExecPolicy_Cuda(y, stream_exec_policy, reduce_exec_policy);
#elif defined(AMREX_USE_HIP)
    SUNHipExecPolicy* stream_exec_policy = new SUNHipThreadDirectExecPolicy(512, stream);
    SUNHipExecPolicy* reduce_exec_policy = new SUNHipBlockReduceExecPolicy(512, 0, stream);
    N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
#endif

    /* Define vectors to be used later in creact */
    BL_PROFILE_VAR_START(AllocsARKODE);
    user_data->rhoe_init_d   = (double*) The_Device_Arena()->alloc(NCELLS* sizeof(double));
    user_data->rhoesrc_ext_d = (double*) The_Device_Arena()->alloc(NCELLS* sizeof(double));
    user_data->rYsrc_d       = (double*) The_Device_Arena()->alloc(NCELLS*NUM_SPECIES*sizeof(double));
    BL_PROFILE_VAR_STOP(AllocsARKODE);

    /* Get Device pointer of solution vector */
#if defined(AMREX_USE_CUDA)
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y);
#elif defined(AMREX_USE_HIP)
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Hip(y);
#else
    Abort("No device arrary pointer");
#endif

    BL_PROFILE_VAR("reactor::FlatStuff", FlatStuff);
    /* Fill the full box_ncells length vectors from input Array4*/
    const auto len        = amrex::length(box);
    const auto lo         = amrex::lbound(box);
    amrex::ParallelFor(box,
    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        int icell = (k-lo.z)*len.x*len.y + (j-lo.y)*len.x + (i-lo.x);
        box_flatten(icell, i, j, k, user_data->ireactor_type,
                    rY_in, rY_src_in, T_in,
                    rEner_in, rEner_src_in,
                    yvec_d,
                    user_data->rYsrc_d,
                    user_data->rhoe_init_d,
                    user_data->rhoesrc_ext_d);
    });
    BL_PROFILE_VAR_STOP(FlatStuff);

    /* Initial time and time to reach after integration */
    time_init = time;
    time_out  = time + dt_react;

    if(use_erkstep == 0)
    {
        arkode_mem = ARKStepCreate(cF_RHS, NULL, time, y);
        flag = ARKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
        flag = ARKStepSStolerances(arkode_mem, relTol, absTol);
        flag = ARKStepResStolerance(arkode_mem, absTol);
        /* call integrator */
        flag = ARKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
    }
    else
    {
        arkode_mem = ERKStepCreate(cF_RHS, time, y);
        flag = ERKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
        flag = ERKStepSStolerances(arkode_mem, relTol, absTol);
        /* call integrator */
        flag = ERKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
    }

#ifdef MOD_REACTOR
    /* If reactor mode is activated, update time */
    dt_react = time_init - time;
    time  = time_init;
#endif

    /* Get estimate of how hard the integration process was */
    long int nfe,nfi;
    if(use_erkstep==0)
    {
        flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
    }
    else
    {
        flag = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
    }
    
    /* Update the input/output Array4 rY_in and rEner_in*/
    BL_PROFILE_VAR_START(FlatStuff);
    amrex::ParallelFor(box, 
    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        int icell = (k-lo.z)*len.x*len.y + (j-lo.y)*len.x + (i-lo.x);
        box_unflatten(icell, i, j, k, user_data->ireactor_type,
               rY_in, T_in, rEner_in, rEner_src_in, FC_in,
               yvec_d, user_data->rhoe_init_d, nfe, dt_react);
    });
    BL_PROFILE_VAR_STOP(FlatStuff);

    //cleanup
    N_VDestroy(y);          /* Free the y vector */
    if(use_erkstep==0)
    {
        ARKStepFree(&arkode_mem);
    }
    else
    {
        ERKStepFree(&arkode_mem);
    }

    The_Device_Arena()->free(user_data->rhoe_init_d);
    The_Device_Arena()->free(user_data->rhoesrc_ext_d);
    The_Device_Arena()->free(user_data->rYsrc_d);

    The_Arena()->free(user_data);

    return nfe;
}
/******************************************************************************************/
/* react call with 1d array of data */
int react(realtype *rY_in, realtype *rY_src_in,
        realtype *rX_in, realtype *rX_src_in,
        realtype &dt_react, realtype &time,
        int reactor_type, int Ncells, amrex::gpuStream_t stream)
{

    int NCELLS, NEQ, neq_tot,flag;
    realtype time_init, time_out;
    void *arkode_mem    = NULL;
    N_Vector y         = NULL;

    NEQ            = NUM_SPECIES+1;
    NCELLS         = Ncells;
    neq_tot        = NEQ * NCELLS;

    /* User data */
    UserData user_data;
    BL_PROFILE_VAR("AllocsInARKODE", AllocsARKODE);
    user_data = (ARKODEUserData *) The_Arena()->alloc(sizeof(struct ARKODEUserData));
    BL_PROFILE_VAR_STOP(AllocsARKODE);
    user_data->ncells_d             = NCELLS;
    user_data->neqs_per_cell        = NEQ;
    user_data->ireactor_type        = reactor_type; 
    user_data->iverbose             = 1;
    user_data->stream               = stream;
    user_data->nbBlocks             = std::max(1,NCELLS/32);
    user_data->nbThreads            = 32;

    /* Definition of main vector */
#if defined(AMREX_USE_CUDA)
    y = N_VNewWithMemHelp_Cuda(neq_tot, /*use_managed_mem=*/true, *The_SUNMemory_Helper());
    if(check_flag((void*)y, "N_VNewWithMemHelp_Cuda", 0)) return(1);
#elif defined(AMREX_USE_HIP)
    y = N_VNewWithMemHelp_Hip(neq_tot, /*use_managed_mem=*/true, *The_SUNMemory_Helper());
    if(check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0)) return(1);
#endif

    /* Use a non-default stream for kernel execution */
#if defined(AMREX_USE_CUDA)
    SUNCudaExecPolicy* stream_exec_policy = new SUNCudaThreadDirectExecPolicy(256, stream);
    SUNCudaExecPolicy* reduce_exec_policy = new SUNCudaBlockReduceExecPolicy(256, 0, stream);
    N_VSetKernelExecPolicy_Cuda(y, stream_exec_policy, reduce_exec_policy);
#elif defined(AMREX_USE_HIP)
    SUNHipExecPolicy* stream_exec_policy = new SUNHipThreadDirectExecPolicy(256, stream);
    SUNHipExecPolicy* reduce_exec_policy = new SUNHipBlockReduceExecPolicy(256, 0, stream);
    N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
#endif

    /* Define vectors to be used later in creact */
    BL_PROFILE_VAR_START(AllocsARKODE);
    user_data->rhoe_init_d   = (double*) The_Device_Arena()->alloc(NCELLS* sizeof(double));
    user_data->rhoesrc_ext_d = (double*) The_Device_Arena()->alloc(NCELLS* sizeof(double));
    user_data->rYsrc_d       = (double*) The_Device_Arena()->alloc(NCELLS*NUM_SPECIES*sizeof(double));
    BL_PROFILE_VAR_STOP(AllocsARKODE);

    /* Get Device pointer of solution vector */
#if defined(AMREX_USE_CUDA)
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y);
#elif defined(AMREX_USE_HIP)
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Hip(y);
#else
    Abort("No device arrary pointer");
#endif

    BL_PROFILE_VAR("AsyncCpy", AsyncCpy);
    // rhoY,T
    Gpu::htod_memcpy_async(yvec_d, rY_in, sizeof(realtype) * (NEQ*NCELLS));
    // rhoY_src_ext
    Gpu::htod_memcpy_async(user_data->rYsrc_d, rY_src_in, (NUM_SPECIES*NCELLS) *sizeof(realtype) );
    // rhoE/rhoH
    Gpu::htod_memcpy_async(user_data->rhoe_init_d,rX_in, sizeof(realtype) * NCELLS);
    Gpu::htod_memcpy_async(user_data->rhoesrc_ext_d, rX_src_in, sizeof(realtype) * NCELLS);
    BL_PROFILE_VAR_STOP(AsyncCpy);

    /* Initial time and time to reach after integration */
    time_init = time;
    time_out  = time + dt_react;
    
    if(use_erkstep == 0)
    {
        arkode_mem = ARKStepCreate(cF_RHS, NULL, time, y);
        flag = ARKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
        flag = ARKStepSStolerances(arkode_mem, relTol, absTol);
        flag = ARKStepResStolerance(arkode_mem, absTol);
        flag = ARKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);      /* call integrator */
    }
    else
    {
        arkode_mem = ERKStepCreate(cF_RHS, time, y);
        flag = ERKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
        flag = ERKStepSStolerances(arkode_mem, relTol, absTol);
        flag = ERKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);      /* call integrator */
    }


#ifdef MOD_REACTOR
    /* If reactor mode is activated, update time */
    dt_react = time_init - time;
    time  = time_init;
#endif

    /* Pack data to return in main routine external */
    BL_PROFILE_VAR_START(AsyncCpy);
    Gpu::dtoh_memcpy_async(rY_in, yvec_d, (NEQ*NCELLS)*sizeof(realtype));
    for  (int i = 0; i < NCELLS; i++)
    {
        rX_in[i] = rX_in[i] + dt_react * rX_src_in[i];
    }
    BL_PROFILE_VAR_STOP(AsyncCpy);


    /* Get estimate of how hard the integration process was */
    long int nfe,nfi;
    if(use_erkstep==0)
    {
        flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
    }
    else
    {
        flag = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
    }

    //cleanup
    N_VDestroy(y);          /* Free the y vector */
    if(use_erkstep==0)
    {
        ARKStepFree(&arkode_mem);
    }
    else
    {
        ERKStepFree(&arkode_mem);
    }

    The_Device_Arena()->free(user_data->rhoe_init_d);
    The_Device_Arena()->free(user_data->rhoesrc_ext_d);
    The_Device_Arena()->free(user_data->rYsrc_d);

    The_Arena()->free(user_data);

    return nfe;
}
/******************************************************************************************/
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in,
        void *user_data)
{

    BL_PROFILE_VAR("fKernelSpec()", fKernelSpec);

    /* Get Device pointers for Kernel call */
#if defined(AMREX_USE_CUDA)
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y_in);
    realtype *ydot_d      = N_VGetDeviceArrayPointer_Cuda(ydot_in);
#elif defined(AMREX_USE_HIP)
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Hip(y_in);
    realtype *ydot_d      = N_VGetDeviceArrayPointer_Hip(ydot_in);
#endif

    // allocate working space 
    UserData udata = static_cast<ARKODEUserData*>(user_data);
    udata->dt_save = t;

    const auto ec = Gpu::ExecutionConfig(udata->ncells_d);
    //launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, udata->stream>>>(
    launch_global<<<udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>
        ([=] AMREX_GPU_DEVICE () noexcept 
    {
         for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 icell < udata->ncells_d; icell += stride)
         {
             fKernelSpec(icell, udata->dt_save, udata->ireactor_type, yvec_d, ydot_d, udata->rhoe_init_d,
                 udata->rhoesrc_ext_d, udata->rYsrc_d);
         }

    });

    Gpu::Device::streamSynchronize();

    BL_PROFILE_VAR_STOP(fKernelSpec);

    return(0);
}
/******************************************************************************************/
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void 
fKernelSpec(int icell, double dt_save, int reactor_type,
        realtype *yvec_d, realtype *ydot_d,  
        double *rhoe_init, double *rhoesrc_ext, double *rYs)
{
    GpuArray<Real,NUM_SPECIES> mw;
    GpuArray<Real,NUM_SPECIES> massfrac;
    GpuArray<Real,NUM_SPECIES> ei_pt;
    GpuArray<Real,NUM_SPECIES> cdots_pt;
    Real Cv_pt = 0.0;
    Real rho_pt = 0.0;
    Real temp_pt = 0.0;
    Real nrg_pt = 0.0;

    int offset = icell * (NUM_SPECIES + 1);

    /* MW CGS */
    get_mw(mw.arr);

    /* rho */ 
    rho_pt = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) 
    {
        rho_pt = rho_pt + yvec_d[offset + n];
    }

    /* Yks, C CGS*/
    for (int i = 0; i < NUM_SPECIES; i++)
    {
        massfrac[i] = yvec_d[offset + i] / rho_pt;
    }

    /* NRG CGS */
    nrg_pt = (rhoe_init[icell] + rhoesrc_ext[icell]*dt_save) /rho_pt;

    /* temp */
    temp_pt = yvec_d[offset + NUM_SPECIES];

    /* Additional var needed */
    auto eos = pele::physics::PhysicsType::eos();
    if (reactor_type == 1)
    {
        /* UV REACTOR */
        eos.EY2T(nrg_pt, massfrac.arr, temp_pt);
        eos.T2Ei(temp_pt, ei_pt.arr);
        eos.TY2Cv(temp_pt, massfrac.arr, Cv_pt);
    }
    else 
    {
        /* HP REACTOR */
        eos.HY2T(nrg_pt, massfrac.arr, temp_pt);
        eos.TY2Cp(temp_pt, massfrac.arr, Cv_pt);
        eos.T2Hi(temp_pt, ei_pt.arr);
    }

    eos.RTY2WDOT(rho_pt, temp_pt, massfrac.arr, cdots_pt.arr);

    /* Fill ydot vect */
    ydot_d[offset + NUM_SPECIES] = rhoesrc_ext[icell];
    for (int i = 0; i < NUM_SPECIES; i++)
    {
        ydot_d[offset + i]           = cdots_pt[i] + rYs[icell * NUM_SPECIES + i];
        ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES]
            - ydot_d[offset + i] * ei_pt[i];
    }
    ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES] /(rho_pt * Cv_pt);
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
      if (ParallelDescriptor::IOProcessor()) {
          fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                  funcname);
          Abort("abort");
      }
      return(1);
  }
  /* Check if flag < 0 */
  else if (opt == 1) {
      errflag = (int *) flagvalue;
      if (*errflag < 0) {
          if (ParallelDescriptor::IOProcessor()) {
              fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                      funcname, *errflag);
              Abort("abort");
          }
          return(1);
      }
  }
  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
      if (ParallelDescriptor::IOProcessor()) {
          fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                  funcname);
          Abort("abort");
      }
      return(1);
  }

  return(0);
}
