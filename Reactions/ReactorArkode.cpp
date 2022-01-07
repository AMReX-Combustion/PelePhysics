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
  pp.query("rk_method", rk_method);
  pp.query("rk_controller", rk_controller);
  std::string method_string = "ZONNEVELD_5_3_4";
  std::string controller_string = "PID";

  switch (rk_method) {
  case 20:
    rk_method = HEUN_EULER_2_1_2;
    method_string = "HEUN_EULER_2_1_2";
    break;
  case 30:
    rk_method = BOGACKI_SHAMPINE_4_2_3;
    method_string = "BOGACKI_SHAMPINE_4_2_3";
    break;
  case 31:
    rk_method = ARK324L2SA_ERK_4_2_3;
    method_string = "ARK324L2SA_ERK_4_2_3";
    break;
  case 40:
    rk_method = ZONNEVELD_5_3_4;
    method_string = "ZONNEVELD_5_3_4";
    break;
  case 41:
    rk_method = ARK436L2SA_ERK_6_3_4;
    method_string = "ARK436L2SA_ERK_6_3_4";
    break;
  case 42:
    rk_method = SAYFY_ABURUB_6_3_4;
    method_string = "SAYFY_ABURUB_6_3_4";
    break;
  case 43:
    rk_method = ARK437L2SA_ERK_7_3_4;
    method_string = "ARK437L2SA_ERK_7_3_4";
    break;
  case 50:
    rk_method = CASH_KARP_6_4_5;
    method_string = "CASH_KARP_6_4_5";
    break;
  case 51:
    rk_method = FEHLBERG_6_4_5;
    method_string = "FEHLBERG_6_4_5";
    break;
  case 52:
    rk_method = DORMAND_PRINCE_7_4_5;
    method_string = "DORMAND_PRINCE_7_4_5";
    break;
  case 53:
    rk_method = ARK548L2SA_ERK_8_4_5;
    method_string = "ARK548L2SA_ERK_8_4_5";
    break;
  case 54:
    rk_method = ARK548L2SAb_ERK_8_4_5;
    method_string = "ARK548L2SAb_ERK_8_4_5";
    break;
  case 60:
    rk_method = VERNER_8_5_6;
    method_string = "VERNER_8_5_6";
    break;
  case 80:
    rk_method = FEHLBERG_13_7_8;
    method_string = "FEHLBERG_13_7_8";
    break;
  default:
    rk_method = ZONNEVELD_5_3_4;
    method_string = "ZONNEVELD_5_3_4";
    break;
  }

  switch (rk_controller) {
  case 0:
    rk_controller = ARK_ADAPT_PID;
    controller_string = "PID";
    break;
  case 1:
    rk_controller = ARK_ADAPT_PI;
    controller_string = "PI";
    break;
  case 2:
    rk_controller = ARK_ADAPT_I;
    controller_string = "I";
    break;
  case 3:
    rk_controller = ARK_ADAPT_EXP_GUS;
    controller_string = "explicit Gustafsson";
    break;
  default:
    rk_controller = ARK_ADAPT_PID;
    controller_string = "PID";
    break;
  }

  if (use_erkstep == 1) {
    amrex::Print() << "ERK Step:" << std::endl;
  } else {
    amrex::Print() << "ARK Step:" << std::endl;
  }
  amrex::Print() << "  Using " << method_string << " method" << std::endl;
  amrex::Print() << "  Using the " << controller_string << " controller"
                 << std::endl;

  return (0);
}

int
ReactorArkode::react(
  const amrex::Box& box,
  amrex::Array4<amrex::Real> const& rY_in,
  amrex::Array4<amrex::Real> const& rYsrc_in,
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
    neq_tot, false, *amrex::sundials::The_SUNMemory_Helper(),
    *amrex::sundials::The_Sundials_Context());
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
    neq_tot, false, *amrex::sundials::The_SUNMemory_Helper(),
    *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0))
    return (1);
  SUNHipExecPolicy* stream_exec_policy =
    new SUNHipThreadDirectExecPolicy(256, stream);
  SUNHipExecPolicy* reduce_exec_policy =
    new SUNHipBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
  realtype* yvec_d = N_VGetDeviceArrayPointer_Hip(y);
#elif defined(AMREX_USE_DPCPP)
  N_Vector y = N_VNewWithMemHelp_Sycl(
    neq_tot, false, *amrex::sundials::The_SUNMemory_Helper(),
    &amrex::Gpu::Device::streamQueue(),
    *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Sycl", 0))
    return (1);
  SUNSyclExecPolicy* stream_exec_policy =
    new SUNSyclThreadDirectExecPolicy(256);
  SUNSyclExecPolicy* reduce_exec_policy =
    new SUNSyclBlockReduceExecPolicy(256, 0);
  N_VSetKernelExecPolicy_Sycl(y, stream_exec_policy, reduce_exec_policy);
  realtype* yvec_d = N_VGetDeviceArrayPointer_Sycl(y);
#else
  N_Vector y = N_VNew_Serial(neq_tot, *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNew_Serial", 0)) {
    return (1);
  }
  realtype* yvec_d = N_VGetArrayPointer(y);
#endif

  amrex::ParmParse pp("ode");
  int verbose = 0;
  pp.query("verbose", verbose);
  const auto captured_reactor_type = m_reactor_type;
  auto* user_data = new ARKODEUserData{};
  amrex::Gpu::DeviceVector<amrex::Real> v_rhoe_init(ncells, 0);
  amrex::Gpu::DeviceVector<amrex::Real> v_rhoesrc_ext(ncells, 0);
  amrex::Gpu::DeviceVector<amrex::Real> v_rYsrc_ext(ncells * NUM_SPECIES, 0);
  user_data->ncells = ncells;
  user_data->neq = neq;
  user_data->reactor_type = captured_reactor_type;
  user_data->verbose = verbose;
  user_data->rhoe_init = v_rhoe_init.begin();
  user_data->rhoesrc_ext = v_rhoesrc_ext.begin();
  user_data->rYsrc_ext = v_rYsrc_ext.begin();

  flatten(
    box, ncells, rY_in, rYsrc_in, T_in, rEner_in, rEner_src_in, yvec_d,
    user_data->rYsrc_ext, user_data->rhoe_init, user_data->rhoesrc_ext);

  realtype time_init = time;
  realtype time_out = time + dt_react;

  void* arkode_mem = nullptr;
  if (use_erkstep == 0) {
    arkode_mem = ARKStepCreate(
      cF_RHS, nullptr, time, y, *amrex::sundials::The_Sundials_Context());
    ARKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
    set_sundials_solver_tols(
      *amrex::sundials::The_Sundials_Context(), arkode_mem, user_data->ncells,
      user_data->verbose, relTol, absTol, "arkstep");
    ARKStepSetTableNum(
      arkode_mem, ARKODE_DIRK_NONE, static_cast<ARKODE_ERKTableID>(rk_method));
    ARKStepSetAdaptivityMethod(arkode_mem, rk_controller, 1, 0, nullptr);
    ARKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
  } else {
    arkode_mem =
      ERKStepCreate(cF_RHS, time, y, *amrex::sundials::The_Sundials_Context());
    ERKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
    set_sundials_solver_tols(
      *amrex::sundials::The_Sundials_Context(), arkode_mem, user_data->ncells,
      user_data->verbose, relTol, absTol, "erkstep");
    ERKStepSetTableNum(arkode_mem, static_cast<ARKODE_ERKTableID>(rk_method));
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

  if (user_data->verbose > 1) {
    print_final_stats(arkode_mem);
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

  delete user_data;

  return nfe;
}

// React for 1d array
int
ReactorArkode::react(
  realtype* rY_in,
  realtype* rYsrc_in,
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
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(),
    *amrex::sundials::The_Sundials_Context());
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
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(),
    *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0))
    return (1);
  SUNHipExecPolicy* stream_exec_policy =
    new SUNHipThreadDirectExecPolicy(256, stream);
  SUNHipExecPolicy* reduce_exec_policy =
    new SUNHipBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
  realtype* yvec_d = N_VGetDeviceArrayPointer_Hip(y);
#elif defined(AMREX_USE_DPCPP)
  N_Vector y = N_VNewWithMemHelp_Sycl(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(),
    &amrex::Gpu::Device::streamQueue(),
    *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Sycl", 0))
    return (1);
  SUNSyclExecPolicy* stream_exec_policy =
    new SUNSyclThreadDirectExecPolicy(256);
  SUNSyclExecPolicy* reduce_exec_policy =
    new SUNSyclBlockReduceExecPolicy(256, 0);
  N_VSetKernelExecPolicy_Sycl(y, stream_exec_policy, reduce_exec_policy);
  realtype* yvec_d = N_VGetDeviceArrayPointer_Sycl(y);
#else
  N_Vector y = N_VNew_Serial(neq_tot, *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNew_Serial", 0)) {
    return (1);
  }
  realtype* yvec_d = N_VGetArrayPointer(y);
#endif

  amrex::ParmParse pp("ode");
  int verbose = 0;
  pp.query("verbose", verbose);
  const auto captured_reactor_type = m_reactor_type;
  auto* user_data = new ARKODEUserData{};
  amrex::Gpu::DeviceVector<amrex::Real> v_rhoe_init(ncells, 0);
  amrex::Gpu::DeviceVector<amrex::Real> v_rhoesrc_ext(ncells, 0);
  amrex::Gpu::DeviceVector<amrex::Real> v_rYsrc_ext(ncells * NUM_SPECIES, 0);
  user_data->ncells = ncells;
  user_data->neq = neq;
  user_data->reactor_type = captured_reactor_type;
  user_data->verbose = verbose;
  user_data->rhoe_init = v_rhoe_init.begin();
  user_data->rhoesrc_ext = v_rhoesrc_ext.begin();
  user_data->rYsrc_ext = v_rYsrc_ext.begin();

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy_async(
    yvec_d, rY_in, sizeof(realtype) * (neq * ncells));
  amrex::Gpu::htod_memcpy_async(
    user_data->rYsrc_ext, rYsrc_in, (NUM_SPECIES * ncells) * sizeof(realtype));
  amrex::Gpu::htod_memcpy_async(
    user_data->rhoe_init, rX_in, sizeof(realtype) * ncells);
  amrex::Gpu::htod_memcpy_async(
    user_data->rhoesrc_ext, rX_src_in, sizeof(realtype) * ncells);
#else
  std::memcpy(yvec_d, rY_in, sizeof(realtype) * (neq * ncells));
  std::memcpy(
    user_data->rYsrc_ext, rYsrc_in, (NUM_SPECIES * ncells) * sizeof(realtype));
  std::memcpy(user_data->rhoe_init, rX_in, sizeof(realtype) * ncells);
  std::memcpy(user_data->rhoesrc_ext, rX_src_in, sizeof(realtype) * ncells);
#endif

  realtype time_init = time;
  realtype time_out = time + dt_react;

  void* arkode_mem = nullptr;
  if (use_erkstep == 0) {
    arkode_mem = ARKStepCreate(
      cF_RHS, nullptr, time, y, *amrex::sundials::The_Sundials_Context());
    ARKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
    set_sundials_solver_tols(
      *amrex::sundials::The_Sundials_Context(), arkode_mem, user_data->ncells,
      user_data->verbose, relTol, absTol, "arkstep");
    ARKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
  } else {
    arkode_mem =
      ERKStepCreate(cF_RHS, time, y, *amrex::sundials::The_Sundials_Context());
    ERKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
    set_sundials_solver_tols(
      *amrex::sundials::The_Sundials_Context(), arkode_mem, user_data->ncells,
      user_data->verbose, relTol, absTol, "erkstep");
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

  if (user_data->verbose > 1) {
    print_final_stats(arkode_mem);
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
#elif defined(AMREX_USE_DPCPP)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Sycl(y_in);
  realtype* ydot_d = N_VGetDeviceArrayPointer_Sycl(ydot_in);
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
  auto* rYsrc_ext = udata->rYsrc_ext;
  amrex::ParallelFor(udata->ncells, [=] AMREX_GPU_DEVICE(int icell) noexcept {
    utils::fKernelSpec<Ordering>(
      icell, ncells, dt_save, reactor_type, yvec_d, ydot_d, rhoe_init,
      rhoesrc_ext, rYsrc_ext);
  });

  amrex::Gpu::Device::streamSynchronize();

  return (0);
}

void
ReactorArkode::print_final_stats(void* arkode_mem)
{
  long int nst, nst_a, nfe, nfi;
  long lenrw, leniw;
  int flag;

  if (use_erkstep) {
    flag = ERKStepGetWorkSpace(arkode_mem, &lenrw, &leniw);
    utils::check_flag(&flag, "ERKStepGetWorkSpace", 1);
    flag = ERKStepGetNumSteps(arkode_mem, &nst);
    utils::check_flag(&flag, "ERKStepGetNumSteps", 1);
    flag = ERKStepGetNumStepAttempts(arkode_mem, &nst_a);
    utils::check_flag(&flag, "ERKStepGetNumStepAttempts", 1);
    flag = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
    utils::check_flag(&flag, "ERKStepGetNumRhsEvals", 1);

  } else {
    flag = ARKStepGetWorkSpace(arkode_mem, &lenrw, &leniw);
    utils::check_flag(&flag, "ARKStepGetWorkSpace", 1);
    flag = ARKStepGetNumSteps(arkode_mem, &nst);
    utils::check_flag(&flag, "ARKStepGetNumSteps", 1);
    flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
    utils::check_flag(&flag, "ARKStepGetNumStepAttempts", 1);
    flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
    utils::check_flag(&flag, "ARKStepGetNumRhsEvals", 1);
  }

#ifdef AMREX_USE_OMP
  amrex::Print() << "\nFinal Statistics: "
                 << "(thread:" << omp_get_thread_num() << ", ";
  amrex::Print() << "arkodeMem:" << arkode_mem << ")\n";
#else
  amrex::Print() << "\nFinal Statistics:\n";
#endif

  amrex::Print() << "   Internal solver steps = " << nst
                 << " (attempted = " << nst_a << ")\n";
  amrex::Print() << "   Total RHS evals:  Fe = " << nfe << "\n";
  amrex::Print() << "lenrw      = " << lenrw << "    leniw         = " << leniw
                 << "\n";
}
} // namespace reactions
} // namespace physics
} // namespace pele
