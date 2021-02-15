#include <reactor.h> 
#include <AMReX_ParmParse.H>
#include <chemistry_file.H>
#include "mechanism.h"
#include <EOS.H>
#include <AMReX_Gpu.H>
#include <AMReX_SUNMemory.H>
#include <AMREX_misc.H>

using namespace amrex;

AMREX_GPU_DEVICE_MANAGED  int eint_rho = 1; // in/out = rhoE/rhoY
AMREX_GPU_DEVICE_MANAGED  int enth_rho = 2; // in/out = rhoH/rhoY 
AMREX_GPU_DEVICE_MANAGED int use_erkstep=0;
AMREX_GPU_DEVICE_MANAGED Real relTol    = 1.0e-6;
AMREX_GPU_DEVICE_MANAGED Real absTol    = 1.0e-10;
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
    EOS::speciesNames(kname);

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
          cudaStream_t stream) 
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
    user_data = (ARKODEUserData *) The_Managed_Arena()->alloc(sizeof(struct ARKODEUserData));  
    BL_PROFILE_VAR_STOP(AllocsARKODE);
    user_data->ncells_d[0]             = NCELLS;
    user_data->neqs_per_cell[0]        = NEQ;
    user_data->ireactor_type           = reactor_type; 
    user_data->iverbose                = 1;
    user_data->stream                  = stream;
    user_data->nbBlocks                = std::max(1,NCELLS/32);
    user_data->nbThreads               = 32;

    y = N_VMakeWithManagedAllocator_Cuda(neq_tot, sunalloc, sunfree);
    N_VSetCudaStream_Cuda(y, &stream);


    BL_PROFILE_VAR_START(AllocsARKODE);
    /* Define vectors to be used later in creact */
    user_data->rhoe_init   = (double *) sunalloc(NCELLS* sizeof(double));
    user_data->rhoesrc_ext = (double *) sunalloc(NCELLS* sizeof(double));
    user_data->rYsrc       = (double *) sunalloc(NCELLS*NUM_SPECIES*sizeof(double));
    BL_PROFILE_VAR_STOP(AllocsARKODE);

    /* Get Device MemCpy of in arrays */
    /* Get Device pointer of solution vector */
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y);
    
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
                    yvec_d, user_data->rYsrc, 
                    user_data->rhoe_init, user_data->rhoesrc_ext);
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
    
    BL_PROFILE_VAR_START(FlatStuff);
    /* Update the input/output Array4 rY_in and rEner_in*/
    amrex::ParallelFor(box, 
    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
    {
        int icell = (k-lo.z)*len.x*len.y + (j-lo.y)*len.x + (i-lo.x);
        box_unflatten(icell, i, j, k, user_data->ireactor_type,
               rY_in, T_in, rEner_in, rEner_src_in, FC_in,
               yvec_d, user_data->rhoe_init, nfe, dt_react);
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

    sunfree(user_data->rhoe_init);
    sunfree(user_data->rhoesrc_ext);
    sunfree(user_data->rYsrc);

    The_Managed_Arena()->free(user_data);

    return nfe;
}
/******************************************************************************************/
/* react call with 1d array of data */
int react(realtype *rY_in, realtype *rY_src_in, 
        realtype *rX_in, realtype *rX_src_in,
        realtype &dt_react, realtype &time,
        int reactor_type, int Ncells, cudaStream_t stream)
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
    cudaMallocManaged(&user_data, sizeof(struct ARKODEUserData));
    BL_PROFILE_VAR_STOP(AllocsARKODE);
    user_data->ncells_d[0]             = NCELLS;
    user_data->neqs_per_cell[0]        = NEQ;
    user_data->ireactor_type           = reactor_type; 
    user_data->iverbose                = 1;
    user_data->stream                  = stream;
    user_data->nbBlocks                = std::max(1,NCELLS/32);
    user_data->nbThreads               = 32;

    //y = N_VNewManaged_Cuda(neq_tot);
    y = N_VNewWithMemHelp_Cuda(neq_tot, /*use_managed_mem=*/true, *amrex::sundials::The_SUNMemory_Helper());
    N_VSetCudaStream_Cuda(y, &stream);


    BL_PROFILE_VAR_START(AllocsARKODE);
    /* Define vectors to be used later in creact */
    cudaMalloc(&(user_data->rhoe_init),   NCELLS*sizeof(double));
    cudaMalloc(&(user_data->rhoesrc_ext), NCELLS*sizeof(double));
    cudaMalloc(&(user_data->rYsrc), (NCELLS*NUM_SPECIES)*sizeof(double));
    BL_PROFILE_VAR_STOP(AllocsARKODE);

    /* Get Device MemCpy of in arrays */
    /* Get Device pointer of solution vector */
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y);
    BL_PROFILE_VAR("AsyncCpy", AsyncCpy);
    // rhoY,T
    cudaMemcpyAsync(yvec_d, rY_in, sizeof(realtype) * (NEQ*NCELLS), cudaMemcpyHostToDevice);
    // rhoY_src_ext
    cudaMemcpyAsync(user_data->rYsrc, rY_src_in, (NUM_SPECIES*NCELLS) *sizeof(double), cudaMemcpyHostToDevice);
    // rhoE/rhoH
    cudaMemcpyAsync(user_data->rhoe_init, rX_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
    cudaMemcpyAsync(user_data->rhoesrc_ext, rX_src_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
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
    cudaMemcpyAsync(rY_in, yvec_d, (NEQ*NCELLS)*sizeof(realtype), cudaMemcpyDeviceToHost);

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

    cudaFree(user_data->rhoe_init);
    cudaFree(user_data->rhoesrc_ext);
    cudaFree(user_data->rYsrc);
    cudaFree(user_data);

    return nfe;
}
/******************************************************************************************/
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, 
        void *user_data)
{

    BL_PROFILE_VAR("fKernelSpec()", fKernelSpec);

    cudaError_t cuda_status = cudaSuccess;

    /* Get Device pointers for Kernel call */
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y_in);
    realtype *ydot_d      = N_VGetDeviceArrayPointer_Cuda(ydot_in);

    // allocate working space 
    UserData udata = static_cast<ARKODEUserData*>(user_data);
    udata->dt_save = t;

    const auto ec = Gpu::ExecutionConfig(udata->ncells_d[0]);   
    //launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, udata->stream>>>(
    launch_global<<<udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>
        ([=] AMREX_GPU_DEVICE () noexcept 
    {
         for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 icell < udata->ncells_d[0]; icell += stride) 
         {
             fKernelSpec(icell, user_data, yvec_d, ydot_d, udata->rhoe_init, 
                 udata->rhoesrc_ext, udata->rYsrc);    
         }

    }); 

    cuda_status = cudaStreamSynchronize(udata->stream);  
    assert(cuda_status == cudaSuccess);

    BL_PROFILE_VAR_STOP(fKernelSpec);

    return(0);
}
/******************************************************************************************/
    AMREX_GPU_DEVICE inline void 
fKernelSpec(int icell, void *user_data, 
        realtype *yvec_d, realtype *ydot_d,  
        double *rhoe_init, double *rhoesrc_ext, double *rYs)
{
    UserData udata = static_cast<ARKODEUserData*>(user_data);

    GpuArray<Real,NUM_SPECIES> mw;
    GpuArray<Real,NUM_SPECIES> massfrac;
    GpuArray<Real,NUM_SPECIES> ei_pt;
    GpuArray<Real,NUM_SPECIES> cdots_pt;
    Real Cv_pt, rho_pt, temp_pt, nrg_pt;

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
    nrg_pt = (rhoe_init[icell] + rhoesrc_ext[icell]*(udata->dt_save)) /rho_pt;

    /* temp */
    temp_pt = yvec_d[offset + NUM_SPECIES];

    /* Additional var needed */
    if (udata->ireactor_type == 1)
    {
        /* UV REACTOR */
        EOS::EY2T(nrg_pt, massfrac.arr, temp_pt);
        EOS::T2Ei(temp_pt, ei_pt.arr);
        EOS::TY2Cv(temp_pt, massfrac.arr, Cv_pt);
    }
    else 
    {
        /* HP REACTOR */
        EOS::HY2T(nrg_pt, massfrac.arr, temp_pt);
        EOS::TY2Cp(temp_pt, massfrac.arr, Cv_pt);
        EOS::T2Hi(temp_pt, ei_pt.arr);
    }

    EOS::RTY2WDOT(rho_pt, temp_pt, massfrac.arr, cdots_pt.arr);

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
