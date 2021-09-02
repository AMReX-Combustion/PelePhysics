#include "ReactorArkode.H"

namespace pele {
namespace physics {
namespace reactions {

int
ReactorArkode::init(int reactor_type, int /*Ncells*/)
{
  BL_PROFILE("Pele::ReactorArkode::init()");
  m_reactor_type = reactor_type;
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

  pp.query("rk_method", rk_method);
  pp.query("rk_controller", rk_controller);
  switch (rk_method) {
  case 20:
    rk_method = HEUN_EULER_2_1_2;
    amrex::Print() << "Using HEUN_EULER_2_1_2 method\n";
    break;
  case 30:
    rk_method = BOGACKI_SHAMPINE_4_2_3;
    amrex::Print() << "Using BOGACKI_SHAMPINE_4_2_3 method\n";
    break;
  case 31:
    rk_method = ARK324L2SA_ERK_4_2_3;
    amrex::Print() << "Using ARK324L2SA_ERK_4_2_3 method\n";
    break;
  case 40:
    rk_method = ZONNEVELD_5_3_4;
    amrex::Print() << "Using ZONNEVELD_5_3_4 method\n";
    break;
  case 41:
    rk_method = ARK436L2SA_ERK_6_3_4;
    amrex::Print() << "Using ARK436L2SA_ERK_6_3_4 method\n";
    break;
  case 42:
    rk_method = SAYFY_ABURUB_6_3_4;
    amrex::Print() << "Using SAYFY_ABURUB_6_3_4 method\n";
    break;
  case 43:
    rk_method = ARK437L2SA_ERK_7_3_4;
    amrex::Print() << "Using ARK437L2SA_ERK_7_3_4 method\n";
    break;
  case 50:
    rk_method = CASH_KARP_6_4_5;
    amrex::Print() << "Using CASH_KARP_6_4_5 method\n";
    break;
  case 51:
    rk_method = FEHLBERG_6_4_5;
    amrex::Print() << "Using FEHLBERG_6_4_5 method\n";
    break;
  case 52:
    rk_method = DORMAND_PRINCE_7_4_5;
    amrex::Print() << "Using DORMAND_PRINCE_7_4_5 method\n";
    break;
  case 53:
    rk_method = ARK548L2SA_ERK_8_4_5;
    amrex::Print() << "Using ARK548L2SA_ERK_8_4_5 method\n";
    break;
  case 54:
    rk_method = ARK548L2SAb_ERK_8_4_5;
    amrex::Print() << "Using ARK548L2SAb_ERK_8_4_5 method\n";
    break;
  case 60:
    rk_method = VERNER_8_5_6;
    amrex::Print() << "Using VERNER_8_5_6 method\n";
    break;
  case 80:
    rk_method = FEHLBERG_13_7_8;
    amrex::Print() << "Using FEHLBERG_13_7_8 method\n";
    break;
  default:
    rk_method = ZONNEVELD_5_3_4;
    amrex::Print() << "Using ZONNEVELD_5_3_4 method\n";
    break;
  }

  switch (rk_controller) {
  case 0:
    rk_controller = ARK_ADAPT_PID;
    amrex::Print() << "Using the PID controller\n";
    break;
  case 1:
    rk_controller = ARK_ADAPT_PI;
    amrex::Print() << "Using the PI controller\n";
    break;
  case 2:
    rk_controller = ARK_ADAPT_I;
    amrex::Print() << "Using the I controller\n";
    break;
  case 3:
    rk_controller = ARK_ADAPT_EXP_GUS;
    amrex::Print() << "Using the explicit Gustafsson controller\n";
    break;
  default:
    rk_controller = ARK_ADAPT_PID;
    amrex::Print() << "Using the PID controller\n";
    break;
  }

  return (0);
}

int
ReactorArkode::react(
  const amrex::Box& box,
  amrex::Array4<amrex::Real> const& rY_in,
  amrex::Array4<amrex::Real> const& rY_src_in,
  amrex::Array4<amrex::Real> const& T_in,
  amrex::Array4<amrex::Real> const& rEner_in,
  amrex::Array4<amrex::Real> const& rEner_src_in,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Array4<int> const& /*mask*/,
  amrex::Real& dt_react,
  amrex::Real& time
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  BL_PROFILE("Pele::ReactorArkode::react()");

  void* arkode_mem = NULL;
  N_Vector y = NULL;

  int NEQ = NUM_SPECIES + 1;
  int NCELLS = box.numPts();
  int neq_tot = NEQ * NCELLS;
  AMREX_ASSERT(NCELLS < std::numeric_limits<int>::max());

  ARKODEUserData* user_data;
  user_data =
    (ARKODEUserData*)amrex::The_Arena()->alloc(sizeof(struct ARKODEUserData));
  user_data->ncells_d = NCELLS;
  user_data->neqs_per_cell = NEQ;
  user_data->ireactor_type = m_reactor_type;
  user_data->iverbose = 1;
  //#ifdef AMREX_USE_GPU
  //  user_data->stream = stream;
  //#endif
  //  user_data->nbThreads = 32;
  //  user_data->nbBlocks = std::max(1, NCELLS / user_data->nbThreads);

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

  flatten(
    box, NCELLS, rY_in, rY_src_in, T_in, rEner_in, rEner_src_in, yvec_d,
    user_data->rYsrc_d, user_data->rhoe_init_d, user_data->rhoesrc_ext_d);

  realtype time_init = time;
  realtype time_out = time + dt_react;

  if (use_erkstep == 0) {
    arkode_mem = ARKStepCreate(cF_RHS, NULL, time, y);
    ARKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
    ARKStepSStolerances(arkode_mem, relTol, absTol);
    ARKStepResStolerance(arkode_mem, absTol);
    ARKStepSetTableNum(arkode_mem, -1, rk_method);
    ARKStepSetAdaptivityMethod(arkode_mem, rk_controller, 1, 0, NULL);
    ARKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
  } else {
    arkode_mem = ERKStepCreate(cF_RHS, time, y);
    ERKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
    ERKStepSStolerances(arkode_mem, relTol, absTol);
    ERKStepSetTableNum(arkode_mem, rk_method);
    ERKStepSetAdaptivityMethod(arkode_mem, rk_controller, 1, 0, NULL);
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

  unflatten(
    box, NCELLS, rY_in, T_in, rEner_in, rEner_src_in, FC_in, yvec_d,
    user_data->rhoe_init_d, &nfe, dt_react);

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
ReactorArkode::react(
  realtype* rY_in,
  realtype* rY_src_in,
  realtype* rX_in,
  realtype* rX_src_in,
  realtype& dt_react,
  realtype& time,
  int Ncells
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  BL_PROFILE("Pele::ReactorArkode::react()");
  void* arkode_mem = NULL;
  N_Vector y = NULL;
  int NEQ = NUM_SPECIES + 1;
  int NCELLS = Ncells;
  int neq_tot = NEQ * NCELLS;
  ARKODEUserData* user_data;
  user_data =
    (ARKODEUserData*)amrex::The_Arena()->alloc(sizeof(struct ARKODEUserData));
  user_data->ncells_d = NCELLS;
  user_data->neqs_per_cell = NEQ;
  user_data->ireactor_type = m_reactor_type;
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

  realtype time_init = time;
  realtype time_out = time + dt_react;

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
ReactorArkode::cF_RHS(
  realtype t, N_Vector y_in, N_Vector ydot_in, void* user_data)
{
  BL_PROFILE("Pele::ReactorArkode::cF_RHS()");
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

  ARKODEUserData* udata = static_cast<ARKODEUserData*>(user_data);
  udata->dt_save = t;

  auto ncells = udata->ncells_d;
  auto dt_save = udata->dt_save;
  auto reactor_type = udata->ireactor_type;
  auto rhoe_init = udata->rhoe_init_d;
  auto rhoesrc_ext = udata->rhoesrc_ext_d;
  auto rYsrc = udata->rYsrc_d;
  amrex::ParallelFor(udata->ncells_d, [=] AMREX_GPU_DEVICE(int icell) noexcept {
    fKernelSpec<Ordering>(
      icell, ncells, dt_save, reactor_type, yvec_d, ydot_d, rhoe_init,
      rhoesrc_ext, rYsrc);
  });

  amrex::Gpu::Device::streamSynchronize();

  return (0);
}

void
ReactorArkode::SetTypValsODE(const std::vector<amrex::Real>& /*ExtTypVals*/)
{
  amrex::Print() << "WARNING: ignoring TypVals for this reactor." << std::endl;
}
} // namespace reactions
} // namespace physics
} // namespace pele
