#include "reactor.H"

int eint_rho = 1;
int enth_rho = 2;
int use_erkstep = 0;
amrex::Real relTol = 1.0e-6;
amrex::Real absTol = 1.0e-10;
amrex::Array<amrex::Real, NUM_SPECIES + 1> typVals = {-1};

int
reactor_init(int reactor_type, int Ncells)
{
  BL_PROFILE("Pele::reactor_init()");
  amrex::ParmParse pp("ode");
  pp.query("use_erkstep", use_erkstep);
  pp.query("rtol", relTol);
  pp.query("atol", absTol);
  if (use_erkstep == 1) {
    amrex::Print() << "Using ERK Step\n";
  } else {
    amrex::Print() << "Using ARK Step\n";
  }
  amrex::Print() << "Setting ARK/ERKODE tolerances rtol = " << relTol
                 << " atol = " << absTol << " in PelePhysics \n";
  return (0);
}

// React with array4
int
react(
  const amrex::Box& box,
  amrex::Array4<amrex::Real> const& rY_in,
  amrex::Array4<amrex::Real> const& rY_src_in,
  amrex::Array4<amrex::Real> const& T_in,
  amrex::Array4<amrex::Real> const& rEner_in,
  amrex::Array4<amrex::Real> const& rEner_src_in,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Array4<int> const& mask,
  amrex::Real& dt_react,
  amrex::Real& time,
  const int& reactor_type
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  BL_PROFILE("Pele::react()");
  int NCELLS, NEQ, neq_tot;
  realtype time_init, time_out;

  void* arkode_mem = NULL;
  N_Vector y = NULL;

  NEQ = NUM_SPECIES + 1;
  NCELLS = box.numPts();
  neq_tot = NEQ * NCELLS;
  AMREX_ASSERT(NCELLS < std::numeric_limits<int>::max());

  UserData user_data;
  user_data =
    (ARKODEUserData*)amrex::The_Arena()->alloc(sizeof(struct ARKODEUserData));
  user_data->ncells_d = NCELLS;
  user_data->neqs_per_cell = NEQ;
  user_data->ireactor_type = reactor_type;
  user_data->iverbose = 1;
#ifdef AMREX_USE_GPU
  user_data->stream = stream;
#endif
  user_data->nbBlocks = std::max(1, NCELLS / 32);
  user_data->nbThreads = 32;

#if defined(AMREX_USE_CUDA)
  y = N_VNewWithMemHelp_Cuda(
    neq_tot, /*use_managed_mem=*/true,
    *amrex::sundials::The_SUNMemory_Helper());
  if (check_flag((void*)y, "N_VNewWithMemHelp_Cuda", 0))
    return (1);
  SUNCudaExecPolicy* stream_exec_policy =
    new SUNCudaThreadDirectExecPolicy(256, stream);
  SUNCudaExecPolicy* reduce_exec_policy =
    new SUNCudaBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Cuda(y, stream_exec_policy, reduce_exec_policy);
#elif defined(AMREX_USE_HIP)
  y = N_VNewWithMemHelp_Hip(
    neq_tot, /*use_managed_mem=*/true,
    *amrex::sundials::The_SUNMemory_Helper());
  if (check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0))
    return (1);
  SUNHipExecPolicy* stream_exec_policy =
    new SUNHipThreadDirectExecPolicy(512, stream);
  SUNHipExecPolicy* reduce_exec_policy =
    new SUNHipBlockReduceExecPolicy(512, 0, stream);
  N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
#else
  y = N_VNew_Serial(neq_tot);
  if (check_flag((void*)y, "N_VNew_Serial", 0))
    return (1);
#endif

  user_data->rhoe_init_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
    NCELLS * sizeof(amrex::Real));
  user_data->rhoesrc_ext_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
    NCELLS * sizeof(amrex::Real));
  user_data->rYsrc_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
    NCELLS * NUM_SPECIES * sizeof(amrex::Real));

#if defined(AMREX_USE_CUDA)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Cuda(y);
#elif defined(AMREX_USE_HIP)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Hip(y);
#else
  realtype* yvec_d = N_VGetArrayPointer(y);
#endif

  const auto len = amrex::length(box);
  const auto lo = amrex::lbound(box);
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    int icell = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
    box_flatten(
      icell, i, j, k, user_data->ireactor_type, rY_in, rY_src_in, T_in,
      rEner_in, rEner_src_in, yvec_d, user_data->rYsrc_d,
      user_data->rhoe_init_d, user_data->rhoesrc_ext_d);
  });

  time_init = time;
  time_out = time + dt_react;

  if (use_erkstep == 0) {
    arkode_mem = ARKStepCreate(cF_RHS, NULL, time, y);
    ARKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
    ARKStepSStolerances(arkode_mem, relTol, absTol);
    ARKStepResStolerance(arkode_mem, absTol);
    ARKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
  } else {
    arkode_mem = ERKStepCreate(cF_RHS, time, y);
    ERKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
    ERKStepSStolerances(arkode_mem, relTol, absTol);
    ERKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
  }

#ifdef MOD_REACTOR
  dt_react = time_init - time;
  time = time_init;
#endif

  long int nfe, nfi;
  if (use_erkstep == 0) {
    ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  } else {
    ERKStepGetNumRhsEvals(arkode_mem, &nfe);
  }

  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    int icell = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
    box_unflatten(
      icell, i, j, k, user_data->ireactor_type, rY_in, T_in, rEner_in,
      rEner_src_in, FC_in, yvec_d, user_data->rhoe_init_d, nfe, dt_react);
  });

  N_VDestroy(y);
  if (use_erkstep == 0) {
    ARKStepFree(&arkode_mem);
  } else {
    ERKStepFree(&arkode_mem);
  }

  amrex::The_Device_Arena()->free(user_data->rhoe_init_d);
  amrex::The_Device_Arena()->free(user_data->rhoesrc_ext_d);
  amrex::The_Device_Arena()->free(user_data->rYsrc_d);

  amrex::The_Arena()->free(user_data);

  return nfe;
}

// React for 1d array
int
react(
  realtype* rY_in,
  realtype* rY_src_in,
  realtype* rX_in,
  realtype* rX_src_in,
  realtype& dt_react,
  realtype& time,
  int reactor_type,
  int Ncells
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  BL_PROFILE("Pele::react()");
  int NCELLS, NEQ, neq_tot;
  realtype time_init, time_out;
  void* arkode_mem = NULL;
  N_Vector y = NULL;
  NEQ = NUM_SPECIES + 1;
  NCELLS = Ncells;
  neq_tot = NEQ * NCELLS;
  UserData user_data;
  user_data =
    (ARKODEUserData*)amrex::The_Arena()->alloc(sizeof(struct ARKODEUserData));
  user_data->ncells_d = NCELLS;
  user_data->neqs_per_cell = NEQ;
  user_data->ireactor_type = reactor_type;
  user_data->iverbose = 1;
#ifdef AMREX_USE_GPU
  user_data->stream = stream;
#endif
  user_data->nbBlocks = std::max(1, NCELLS / 32);
  user_data->nbThreads = 32;

#if defined(AMREX_USE_CUDA)
  y = N_VNewWithMemHelp_Cuda(
    neq_tot, /*use_managed_mem=*/true,
    *amrex::sundials::The_SUNMemory_Helper());
  if (check_flag((void*)y, "N_VNewWithMemHelp_Cuda", 0))
    return (1);
  SUNCudaExecPolicy* stream_exec_policy =
    new SUNCudaThreadDirectExecPolicy(256, stream);
  SUNCudaExecPolicy* reduce_exec_policy =
    new SUNCudaBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Cuda(y, stream_exec_policy, reduce_exec_policy);
#elif defined(AMREX_USE_HIP)
  y = N_VNewWithMemHelp_Hip(
    neq_tot, /*use_managed_mem=*/true,
    *amrex::sundials::The_SUNMemory_Helper());
  if (check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0))
    return (1);
  SUNHipExecPolicy* stream_exec_policy =
    new SUNHipThreadDirectExecPolicy(256, stream);
  SUNHipExecPolicy* reduce_exec_policy =
    new SUNHipBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
#else
  y = N_VNew_Serial(neq_tot);
  if (check_flag((void*)y, "N_VNew_Serial", 0))
    return (1);
#endif

  user_data->rhoe_init_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
    NCELLS * sizeof(amrex::Real));
  user_data->rhoesrc_ext_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
    NCELLS * sizeof(amrex::Real));
  user_data->rYsrc_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
    NCELLS * NUM_SPECIES * sizeof(amrex::Real));

#if defined(AMREX_USE_CUDA)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Cuda(y);
#elif defined(AMREX_USE_HIP)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Hip(y);
#else
  realtype* yvec_d = N_VGetArrayPointer(y);
#endif

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy_async(
    yvec_d, rY_in, sizeof(realtype) * (NEQ * NCELLS));
  amrex::Gpu::htod_memcpy_async(
    user_data->rYsrc_d, rY_src_in, (NUM_SPECIES * NCELLS) * sizeof(realtype));
  amrex::Gpu::htod_memcpy_async(
    user_data->rhoe_init_d, rX_in, sizeof(realtype) * NCELLS);
  amrex::Gpu::htod_memcpy_async(
    user_data->rhoesrc_ext_d, rX_src_in, sizeof(realtype) * NCELLS);
#else
  std::memcpy(yvec_d, rY_in, sizeof(realtype) * (NEQ * NCELLS));
  std::memcpy(
    user_data->rYsrc_d, rY_src_in, (NUM_SPECIES * NCELLS) * sizeof(realtype));
  std::memcpy(user_data->rhoe_init_d, rX_in, sizeof(realtype) * NCELLS);
  std::memcpy(user_data->rhoesrc_ext_d, rX_src_in, sizeof(realtype) * NCELLS);
#endif

  time_init = time;
  time_out = time + dt_react;

  if (use_erkstep == 0) {
    arkode_mem = ARKStepCreate(cF_RHS, NULL, time, y);
    ARKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
    ARKStepSStolerances(arkode_mem, relTol, absTol);
    ARKStepResStolerance(arkode_mem, absTol);
    ARKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
  } else {
    arkode_mem = ERKStepCreate(cF_RHS, time, y);
    ERKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
    ERKStepSStolerances(arkode_mem, relTol, absTol);
    ERKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
  }
#ifdef MOD_REACTOR
  dt_react = time_init - time;
  time = time_init;
#endif

#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy_async(
#else
  std::memcpy(
#endif
    rY_in, yvec_d, (NEQ * NCELLS) * sizeof(realtype));
  for (int i = 0; i < NCELLS; i++) {
    rX_in[i] = rX_in[i] + dt_react * rX_src_in[i];
  }

  long int nfe, nfi;
  if (use_erkstep == 0) {
    ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  } else {
    ERKStepGetNumRhsEvals(arkode_mem, &nfe);
  }

  N_VDestroy(y);
  if (use_erkstep == 0) {
    ARKStepFree(&arkode_mem);
  } else {
    ERKStepFree(&arkode_mem);
  }

  amrex::The_Device_Arena()->free(user_data->rhoe_init_d);
  amrex::The_Device_Arena()->free(user_data->rhoesrc_ext_d);
  amrex::The_Device_Arena()->free(user_data->rYsrc_d);

  amrex::The_Arena()->free(user_data);

  return nfe;
}

int
cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, void* user_data)
{
  BL_PROFILE("Pele::cF_RHS()");
#if defined(AMREX_USE_CUDA)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Cuda(y_in);
  realtype* ydot_d = N_VGetDeviceArrayPointer_Cuda(ydot_in);
#elif defined(AMREX_USE_HIP)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Hip(y_in);
  realtype* ydot_d = N_VGetDeviceArrayPointer_Hip(ydot_in);
#else
  realtype* yvec_d = N_VGetArrayPointer(y_in);
  realtype* ydot_d = N_VGetArrayPointer(ydot_in);
#endif

  UserData udata = static_cast<ARKODEUserData*>(user_data);
  udata->dt_save = t;

#ifdef AMREX_USE_GPU
  const auto ec = amrex::Gpu::ExecutionConfig(udata->ncells_d);
  amrex::launch_global<<<
    udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>(
    [=] AMREX_GPU_DEVICE() noexcept {
      for (int icell = blockDim.x * blockIdx.x + threadIdx.x,
               stride = blockDim.x * gridDim.x;
           icell < udata->ncells_d; icell += stride) {
        fKernelSpec(
          icell, udata->dt_save, udata->ireactor_type, yvec_d, ydot_d,
          udata->rhoe_init_d, udata->rhoesrc_ext_d, udata->rYsrc_d);
      }
    });
#else
  for (int icell = 0; icell < udata->ncells_d; icell++) {
    fKernelSpec(
      icell, udata->dt_save, udata->ireactor_type, yvec_d, ydot_d,
      udata->rhoe_init_d, udata->rhoesrc_ext_d, udata->rYsrc_d);
  }
#endif

  amrex::Gpu::Device::streamSynchronize();

  return (0);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fKernelSpec(
  int icell,
  amrex::Real dt_save,
  int reactor_type,
  realtype* yvec_d,
  realtype* ydot_d,
  amrex::Real* rhoe_init,
  amrex::Real* rhoesrc_ext,
  amrex::Real* rYs)
{
  amrex::GpuArray<amrex::Real, NUM_SPECIES> mw = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> ei_pt = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> cdots_pt = {0.0};
  amrex::Real Cv_pt = 0.0;
  amrex::Real rho_pt = 0.0;
  amrex::Real temp_pt = 0.0;
  amrex::Real nrg_pt = 0.0;

  int offset = icell * (NUM_SPECIES + 1);

  get_mw(mw.arr);

  rho_pt = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    rho_pt = rho_pt + yvec_d[offset + n];
  }

  for (int i = 0; i < NUM_SPECIES; i++) {
    massfrac[i] = yvec_d[offset + i] / rho_pt;
  }

  nrg_pt = (rhoe_init[icell] + rhoesrc_ext[icell] * dt_save) / rho_pt;

  temp_pt = yvec_d[offset + NUM_SPECIES];

  auto eos = pele::physics::PhysicsType::eos();
  if (reactor_type == 1) {
    eos.EY2T(nrg_pt, massfrac.arr, temp_pt);
    eos.T2Ei(temp_pt, ei_pt.arr);
    eos.TY2Cv(temp_pt, massfrac.arr, Cv_pt);
  } else {
    eos.HY2T(nrg_pt, massfrac.arr, temp_pt);
    eos.TY2Cp(temp_pt, massfrac.arr, Cv_pt);
    eos.T2Hi(temp_pt, ei_pt.arr);
  }

  eos.RTY2WDOT(rho_pt, temp_pt, massfrac.arr, cdots_pt.arr);

  ydot_d[offset + NUM_SPECIES] = rhoesrc_ext[icell];
  for (int i = 0; i < NUM_SPECIES; i++) {
    ydot_d[offset + i] = cdots_pt[i] + rYs[icell * NUM_SPECIES + i];
    ydot_d[offset + NUM_SPECIES] =
      ydot_d[offset + NUM_SPECIES] - ydot_d[offset + i] * ei_pt[i];
  }
  ydot_d[offset + NUM_SPECIES] =
    ydot_d[offset + NUM_SPECIES] / (rho_pt * Cv_pt);
}

// Check function return value...
//     opt == 0 means SUNDIALS function allocates memory so check if
//              returned NULL pointer
//     opt == 1 means SUNDIALS function returns a flag so check if
//              flag >= 0
//     opt == 2 means function allocates memory so check if returned
//              NULL pointer
int
check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  if (opt == 0 && flagvalue == NULL) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      fprintf(
        stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      amrex::Abort("abort");
    }
    return (1);
  } else if (opt == 1) {
    errflag = (int*)flagvalue;
    if (*errflag < 0) {
      if (amrex::ParallelDescriptor::IOProcessor()) {
        fprintf(
          stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname,
          *errflag);
        amrex::Abort("abort");
      }
      return (1);
    }
  } else if (opt == 2 && flagvalue == NULL) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      fprintf(
        stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      amrex::Abort("abort");
    }
    return (1);
  }

  return (0);
}
