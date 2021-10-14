#include "ReactorArkode.H"

namespace pele {
namespace physics {
namespace reactions {

int
ReactorArkode::init(int reactor_type, int /*ncells*/)
{
  BL_PROFILE("Pele::ReactorArkode::init()");
  m_reactor_type = reactor_type;
  ReactorTypes::check_reactor_type(m_reactor_type);
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
  const int ncells = box.numPts();
  AMREX_ASSERT(ncells < std::numeric_limits<int>::max());

  const int neq = NUM_SPECIES + 1;
  const int neq_tot = neq * ncells;
#if defined(AMREX_USE_CUDA)
  N_Vector y = N_VNewWithMemHelp_Cuda(
    neq_tot, false, *amrex::sundials::The_SUNMemory_Helper());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Cuda", 0))
    return (1);
  SUNCudaExecPolicy* stream_exec_policy =
    new SUNCudaThreadDirectExecPolicy(256, stream);
  SUNCudaExecPolicy* reduce_exec_policy =
    new SUNCudaBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Cuda(y, stream_exec_policy, reduce_exec_policy);
  realtype* yvec_d = N_VGetDeviceArrayPointer_Cuda(y);
#elif defined(AMREX_USE_HIP)
  N_Vector y = N_VNewWithMemHelp_Hip(
    neq_tot, false, *amrex::sundials::The_SUNMemory_Helper());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0))
    return (1);
  SUNHipExecPolicy* stream_exec_policy =
    new SUNHipThreadDirectExecPolicy(512, stream);
  SUNHipExecPolicy* reduce_exec_policy =
    new SUNHipBlockReduceExecPolicy(512, 0, stream);
  N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
  realtype* yvec_d = N_VGetDeviceArrayPointer_Hip(y);
#else
  N_Vector y = N_VNew_Serial(neq_tot);
  if (utils::check_flag((void*)y, "N_VNew_Serial", 0)) {
    return (1);
  }
  realtype* yvec_d = N_VGetArrayPointer(y);
#endif

  const int verbose = 1;
  const auto captured_reactor_type = m_reactor_type;
  ARKODEUserData* user_data =
    (ARKODEUserData*)amrex::The_Arena()->alloc(sizeof(struct ARKODEUserData));
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, &ncells, &ncells + 1, &(user_data->ncells));
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, &neq, &neq + 1, &(user_data->neq));
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, &captured_reactor_type,
    &captured_reactor_type + 1, &(user_data->reactor_type));
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, &verbose, &verbose + 1, &(user_data->verbose));
  amrex::Gpu::DeviceVector<amrex::Real> v_rhoe_init(ncells, 0);
  amrex::Gpu::DeviceVector<amrex::Real> v_rhoesrc_ext(ncells, 0);
  amrex::Gpu::DeviceVector<amrex::Real> v_rYsrc(ncells * NUM_SPECIES, 0);
  const auto p_rhoe_init = v_rhoe_init.begin();
  const auto p_rhoesrc_ext = v_rhoesrc_ext.begin();
  const auto p_rYsrc = v_rYsrc.begin();
  amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int) {
    user_data->rhoe_init = p_rhoe_init;
    user_data->rhoesrc_ext = p_rhoesrc_ext;
    user_data->rYsrc = p_rYsrc;
  });

  flatten(
    box, ncells, rY_in, rY_src_in, T_in, rEner_in, rEner_src_in, yvec_d,
    user_data->rYsrc, user_data->rhoe_init, user_data->rhoesrc_ext);

  realtype time_init = time;
  realtype time_out = time + dt_react;

  void* arkode_mem = nullptr;
  if (use_erkstep == 0) {
    arkode_mem = ARKStepCreate(cF_RHS, nullptr, time, y);
    ARKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
    ARKStepSStolerances(arkode_mem, relTol, absTol);
    ARKStepResStolerance(arkode_mem, absTol);
    ARKStepSetTableNum(arkode_mem, -1, rk_method);
    ARKStepSetAdaptivityMethod(arkode_mem, rk_controller, 1, 0, nullptr);
    ARKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
  } else {
    arkode_mem = ERKStepCreate(cF_RHS, time, y);
    ERKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
    ERKStepSStolerances(arkode_mem, relTol, absTol);
    ERKStepSetTableNum(arkode_mem, rk_method);
    ERKStepSetAdaptivityMethod(arkode_mem, rk_controller, 1, 0, nullptr);
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

  amrex::Gpu::DeviceVector<long int> v_nfe(ncells, nfe);
  long int* d_nfe = v_nfe.data();
  unflatten(
    box, ncells, rY_in, T_in, rEner_in, rEner_src_in, FC_in, yvec_d,
    user_data->rhoe_init, d_nfe, dt_react);

  N_VDestroy(y);
  if (use_erkstep == 0) {
    ARKStepFree(&arkode_mem);
  } else {
    ERKStepFree(&arkode_mem);
  }

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
  int ncells
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  BL_PROFILE("Pele::ReactorArkode::react()");
  AMREX_ASSERT(ncells < std::numeric_limits<int>::max());

  int neq = NUM_SPECIES + 1;
  int neq_tot = neq * ncells;
#if defined(AMREX_USE_CUDA)
  N_Vector y = N_VNewWithMemHelp_Cuda(
    neq_tot, /*use_managed_mem=*/true,
    *amrex::sundials::The_SUNMemory_Helper());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Cuda", 0))
    return (1);
  SUNCudaExecPolicy* stream_exec_policy =
    new SUNCudaThreadDirectExecPolicy(256, stream);
  SUNCudaExecPolicy* reduce_exec_policy =
    new SUNCudaBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Cuda(y, stream_exec_policy, reduce_exec_policy);
  realtype* yvec_d = N_VGetDeviceArrayPointer_Cuda(y);
#elif defined(AMREX_USE_HIP)
  N_Vector y = N_VNewWithMemHelp_Hip(
    neq_tot, /*use_managed_mem=*/true,
    *amrex::sundials::The_SUNMemory_Helper());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0))
    return (1);
  SUNHipExecPolicy* stream_exec_policy =
    new SUNHipThreadDirectExecPolicy(256, stream);
  SUNHipExecPolicy* reduce_exec_policy =
    new SUNHipBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
  realtype* yvec_d = N_VGetDeviceArrayPointer_Hip(y);
#else
  N_Vector y = N_VNew_Serial(neq_tot);
  if (utils::check_flag((void*)y, "N_VNew_Serial", 0)) {
    return (1);
  }
  realtype* yvec_d = N_VGetArrayPointer(y);
#endif

  const int verbose = 1;
  const auto captured_reactor_type = m_reactor_type;
  ARKODEUserData* user_data =
    (ARKODEUserData*)amrex::The_Arena()->alloc(sizeof(struct ARKODEUserData));
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, &ncells, &ncells + 1, &(user_data->ncells));
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, &neq, &neq + 1, &(user_data->neq));
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, &captured_reactor_type,
    &captured_reactor_type + 1, &(user_data->reactor_type));
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, &verbose, &verbose + 1, &(user_data->verbose));
  amrex::Gpu::DeviceVector<amrex::Real> v_rhoe_init(ncells, 0);
  amrex::Gpu::DeviceVector<amrex::Real> v_rhoesrc_ext(ncells, 0);
  amrex::Gpu::DeviceVector<amrex::Real> v_rYsrc(ncells * NUM_SPECIES, 0);
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, v_rhoe_init.data(), v_rhoe_init.data() + 1,
    &(user_data->rhoe_init));
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, v_rhoesrc_ext.data(), v_rhoesrc_ext.data() + 1,
    &(user_data->rhoesrc_ext));
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, v_rYsrc.data(), v_rYsrc.data() + 1,
    &(user_data->rYsrc));

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy_async(
    yvec_d, rY_in, sizeof(realtype) * (neq * ncells));
  amrex::Gpu::htod_memcpy_async(
    user_data->rYsrc, rY_src_in, (NUM_SPECIES * ncells) * sizeof(realtype));
  amrex::Gpu::htod_memcpy_async(
    user_data->rhoe_init, rX_in, sizeof(realtype) * ncells);
  amrex::Gpu::htod_memcpy_async(
    user_data->rhoesrc_ext, rX_src_in, sizeof(realtype) * ncells);
#else
  std::memcpy(yvec_d, rY_in, sizeof(realtype) * (neq * ncells));
  std::memcpy(
    user_data->rYsrc, rY_src_in, (NUM_SPECIES * ncells) * sizeof(realtype));
  std::memcpy(user_data->rhoe_init, rX_in, sizeof(realtype) * ncells);
  std::memcpy(user_data->rhoesrc_ext, rX_src_in, sizeof(realtype) * ncells);
#endif

  realtype time_init = time;
  realtype time_out = time + dt_react;

  void* arkode_mem = nullptr;
  if (use_erkstep == 0) {
    arkode_mem = ARKStepCreate(cF_RHS, nullptr, time, y);
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
    rY_in, yvec_d, (neq * ncells) * sizeof(realtype));
  for (int i = 0; i < ncells; i++) {
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

  auto* udata = static_cast<ARKODEUserData*>(user_data);
  udata->dt_save = t;

  const auto ncells = udata->ncells;
  const auto dt_save = udata->dt_save;
  const auto reactor_type = udata->reactor_type;
  auto* rhoe_init = udata->rhoe_init;
  auto* rhoesrc_ext = udata->rhoesrc_ext;
  auto* rYsrc = udata->rYsrc;
  amrex::ParallelFor(udata->ncells, [=] AMREX_GPU_DEVICE(int icell) noexcept {
    utils::fKernelSpec<Ordering>(
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
