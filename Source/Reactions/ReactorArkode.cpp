#include "ReactorArkode.H"

namespace pele::physics::reactions {

int
ReactorArkode::init(int reactor_type, int /*ncells*/)
{
  BL_PROFILE("Pele::ReactorArkode::init()");
  m_reactor_type = reactor_type;
  ReactorTypes::check_reactor_type(m_reactor_type);
  amrex::ParmParse pp("ode");
  pp.query("verbose", verbose);
  pp.query("use_erkstep", use_erkstep);
  pp.query("rtol", relTol);
  pp.query("atol", absTol);
  pp.query("atomic_reductions", atomic_reductions);
  pp.query("rk_method", rk_method);
  pp.query("rk_controller", rk_controller);
  pp.query("clean_init_massfrac", m_clean_init_massfrac);
  std::string method_string = "ARKODE_ZONNEVELD_5_3_4";
  std::string controller_string = "PID";

  amrex::Print() << "Initializing ARKODE:\n";

  switch (rk_method) {
  case 20:
    rk_method = ARKODE_HEUN_EULER_2_1_2;
    method_string = "ARKODE_HEUN_EULER_2_1_2";
    break;
  case 30:
    rk_method = ARKODE_BOGACKI_SHAMPINE_4_2_3;
    method_string = "ARKODE_BOGACKI_SHAMPINE_4_2_3";
    break;
  case 31:
    rk_method = ARKODE_ARK324L2SA_ERK_4_2_3;
    method_string = "ARKODE_ARK324L2SA_ERK_4_2_3";
    break;
  case 40:
    rk_method = ARKODE_ZONNEVELD_5_3_4;
    method_string = "ARKODE_ZONNEVELD_5_3_4";
    break;
  case 41:
    rk_method = ARKODE_ARK436L2SA_ERK_6_3_4;
    method_string = "ARKODE_ARK436L2SA_ERK_6_3_4";
    break;
  case 42:
    rk_method = ARKODE_SAYFY_ABURUB_6_3_4;
    method_string = "ARKODE_SAYFY_ABURUB_6_3_4";
    break;
  case 43:
    rk_method = ARKODE_ARK437L2SA_ERK_7_3_4;
    method_string = "ARKODE_ARK437L2SA_ERK_7_3_4";
    break;
  case 50:
    rk_method = ARKODE_CASH_KARP_6_4_5;
    method_string = "ARKODE_CASH_KARP_6_4_5";
    break;
  case 51:
    rk_method = ARKODE_FEHLBERG_6_4_5;
    method_string = "ARKODE_FEHLBERG_6_4_5";
    break;
  case 52:
    rk_method = ARKODE_DORMAND_PRINCE_7_4_5;
    method_string = "ARKODE_DORMAND_PRINCE_7_4_5";
    break;
  case 53:
    rk_method = ARKODE_ARK548L2SA_ERK_8_4_5;
    method_string = "ARKODE_ARK548L2SA_ERK_8_4_5";
    break;
  case 54:
    rk_method = ARKODE_ARK548L2SAb_ERK_8_4_5;
    method_string = "ARKODE_ARK548L2SAb_ERK_8_4_5";
    break;
  case 60:
    rk_method = ARKODE_VERNER_8_5_6;
    method_string = "ARKODE_VERNER_8_5_6";
    break;
  case 80:
    rk_method = ARKODE_FEHLBERG_13_7_8;
    method_string = "ARKODE_FEHLBERG_13_7_8";
    break;
  default:
    rk_method = ARKODE_ZONNEVELD_5_3_4;
    method_string = "ARKODE_ZONNEVELD_5_3_4";
    break;
  }

  switch (rk_controller) {
  case 0:
    controller_string = "PID";
    sun_controller =
      SUNAdaptController_PID(*amrex::sundials::The_Sundials_Context());
    break;
  case 1:
    controller_string = "PI";
    sun_controller =
      SUNAdaptController_PI(*amrex::sundials::The_Sundials_Context());
    break;
  case 2:
    controller_string = "I";
    sun_controller =
      SUNAdaptController_I(*amrex::sundials::The_Sundials_Context());
    break;
  case 3:
    controller_string = "explicit Gustafsson";
    sun_controller =
      SUNAdaptController_ExpGus(*amrex::sundials::The_Sundials_Context());
    break;
  default:
    controller_string = "PID";
    sun_controller =
      SUNAdaptController_PID(*amrex::sundials::The_Sundials_Context());
    break;
  }

  if (use_erkstep == 1) {
    amrex::Print() << "  Using ERKStep" << std::endl;
  } else {
    amrex::Print() << "  Using ARKStep" << std::endl;
  }
  amrex::Print() << "  Using " << method_string << " method" << std::endl;
  amrex::Print() << "  Using the " << controller_string << " controller"
                 << std::endl;

  if (atomic_reductions != 0) {
    amrex::Print() << "  Using atomic reductions\n";
  } else {
    amrex::Print() << "  Using LDS reductions\n";
  }

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

  const int ncells = static_cast<int>(box.numPts());
  AMREX_ASSERT(ncells < std::numeric_limits<int>::max());

  const int neq = NUM_SPECIES + 1;
  const int neq_tot = neq * ncells;

#ifdef SUNDIALS_BUILD_WITH_PROFILING
  SUNProfiler sun_profiler = nullptr;
  SUNContext_GetProfiler(
    *amrex::sundials::The_Sundials_Context(), &sun_profiler);
#endif

  // Solution vector and execution policy
#ifdef AMREX_USE_GPU
  auto y = utils::setNVectorGPU(neq_tot, atomic_reductions, stream);
  sunrealtype* yvec_d = N_VGetDeviceArrayPointer(y);
#else
  N_Vector y = N_VNew_Serial(neq_tot, *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag(static_cast<void*>(y), "N_VNew_Serial", 0) != 0) {
    return (1);
  }
  sunrealtype* yvec_d = N_VGetArrayPointer(y);
#endif

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

  sunrealtype time_init = time;
  sunrealtype time_out = time + dt_react;

  void* arkode_mem = nullptr;
  if (use_erkstep == 0) {
    arkode_mem = ARKStepCreate(
      cF_RHS, nullptr, time, y, *amrex::sundials::The_Sundials_Context());
    ARKodeSetUserData(arkode_mem, static_cast<void*>(user_data));
    utils::set_sundials_solver_tols<Ordering>(
      *amrex::sundials::The_Sundials_Context(), arkode_mem, user_data->ncells,
      relTol, absTol, m_typ_vals, "arkstep", verbose);
    ARKStepSetTableNum(
      arkode_mem, ARKODE_DIRK_NONE, static_cast<ARKODE_ERKTableID>(rk_method));
    int flag = ARKodeSetAdaptController(arkode_mem, sun_controller);
    utils::check_flag(&flag, "ARKodeSetAdaptController", 1);
    BL_PROFILE_VAR(
      "Pele::ReactorArkode::react():ARKStepEvolve", AroundARKEvolve);
    ARKodeEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
    BL_PROFILE_VAR_STOP(AroundARKEvolve);
  } else {
    arkode_mem =
      ERKStepCreate(cF_RHS, time, y, *amrex::sundials::The_Sundials_Context());
    ARKodeSetUserData(arkode_mem, static_cast<void*>(user_data));
    utils::set_sundials_solver_tols<Ordering>(
      *amrex::sundials::The_Sundials_Context(), arkode_mem, user_data->ncells,
      relTol, absTol, m_typ_vals, "erkstep", verbose);
    ERKStepSetTableNum(arkode_mem, static_cast<ARKODE_ERKTableID>(rk_method));
    int flag = ARKodeSetAdaptController(arkode_mem, sun_controller);
    utils::check_flag(&flag, "ARKodeSetAdaptController", 1);
    BL_PROFILE_VAR(
      "Pele::ReactorArkode::react():ERKStepEvolve", AroundERKEvolve);
    ARKodeEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
    BL_PROFILE_VAR_STOP(AroundERKEvolve);
  }

#ifdef MOD_REACTOR
  dt_react = time_init - time;
  time = time_init;
#endif

  long int nfe;
  ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);

  if (user_data->verbose > 1) {
    print_final_stats(arkode_mem);
  }

  amrex::Gpu::DeviceVector<long int> v_nfe(ncells, nfe);
  long int* d_nfe = v_nfe.data();
  unflatten(
    box, ncells, rY_in, T_in, rEner_in, rEner_src_in, FC_in, yvec_d,
    user_data->rhoe_init, d_nfe, dt_react);

  N_VDestroy(y);
  ARKodeFree(&arkode_mem);

  delete user_data;

#ifdef SUNDIALS_BUILD_WITH_PROFILING
  SUNProfiler_Print(sun_profiler, stdout);
#endif

  return static_cast<int>(nfe);
}

// React for 1d array
int
ReactorArkode::react(
  sunrealtype* rY_in,
  sunrealtype* rYsrc_in,
  sunrealtype* rX_in,
  sunrealtype* rX_src_in,
  sunrealtype& dt_react,
  sunrealtype& time,
  int ncells
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  BL_PROFILE("Pele::ReactorArkode::react()");
  AMREX_ASSERT(ncells < std::numeric_limits<int>::max());

  std::cout << "Reacting (flattened)\n";

  int neq = NUM_SPECIES + 1;
  int neq_tot = neq * ncells;

#ifdef AMREX_USE_GPU
  auto y = utils::setNVectorGPU(neq_tot, atomic_reductions, stream);
  sunrealtype* yvec_d = N_VGetDeviceArrayPointer(y);
#else
  N_Vector y = N_VNew_Serial(neq_tot, *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag(static_cast<void*>(y), "N_VNew_Serial", 0) != 0) {
    return (1);
  }
  sunrealtype* yvec_d = N_VGetArrayPointer(y);
#endif

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
    yvec_d, rY_in, sizeof(sunrealtype) * (neq * ncells));
  amrex::Gpu::htod_memcpy_async(
    user_data->rYsrc_ext, rYsrc_in,
    (NUM_SPECIES * ncells) * sizeof(sunrealtype));
  amrex::Gpu::htod_memcpy_async(
    user_data->rhoe_init, rX_in, sizeof(sunrealtype) * ncells);
  amrex::Gpu::htod_memcpy_async(
    user_data->rhoesrc_ext, rX_src_in, sizeof(sunrealtype) * ncells);
#else
  std::memcpy(yvec_d, rY_in, sizeof(sunrealtype) * (neq * ncells));
  std::memcpy(
    user_data->rYsrc_ext, rYsrc_in,
    (NUM_SPECIES * ncells) * sizeof(sunrealtype));
  std::memcpy(user_data->rhoe_init, rX_in, sizeof(sunrealtype) * ncells);
  std::memcpy(user_data->rhoesrc_ext, rX_src_in, sizeof(sunrealtype) * ncells);
#endif

  sunrealtype time_init = time;
  sunrealtype time_out = time + dt_react;

  void* arkode_mem = nullptr;
  if (use_erkstep == 0) {
    arkode_mem = ARKStepCreate(
      cF_RHS, nullptr, time, y, *amrex::sundials::The_Sundials_Context());
    ARKodeSetUserData(arkode_mem, static_cast<void*>(user_data));
    utils::set_sundials_solver_tols<Ordering>(
      *amrex::sundials::The_Sundials_Context(), arkode_mem, user_data->ncells,
      relTol, absTol, m_typ_vals, "arkstep", verbose);
    BL_PROFILE_VAR(
      "Pele::ReactorArkode::react():ARKStepEvolve", AroundARKEvolve);
    ARKodeEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
    BL_PROFILE_VAR_STOP(AroundARKEvolve);
  } else {
    arkode_mem =
      ERKStepCreate(cF_RHS, time, y, *amrex::sundials::The_Sundials_Context());
    ARKodeSetUserData(arkode_mem, static_cast<void*>(user_data));
    utils::set_sundials_solver_tols<Ordering>(
      *amrex::sundials::The_Sundials_Context(), arkode_mem, user_data->ncells,
      relTol, absTol, m_typ_vals, "erkstep", verbose);
    BL_PROFILE_VAR(
      "Pele::ReactorArkode::react():ERKStepEvolve", AroundERKEvolve);
    ARKodeEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);
    BL_PROFILE_VAR_STOP(AroundERKEvolve);
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
    rY_in, yvec_d, (neq * ncells) * sizeof(sunrealtype));
  for (int i = 0; i < ncells; i++) {
    rX_in[i] = rX_in[i] + dt_react * rX_src_in[i];
  }

  long int nfe;
  ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);

  if (user_data->verbose > 1) {
    print_final_stats(arkode_mem);
  }

  N_VDestroy(y);
  ARKodeFree(&arkode_mem);

  delete user_data;

  return static_cast<int>(nfe);
}

int
ReactorArkode::cF_RHS(
  sunrealtype t, N_Vector y_in, N_Vector ydot_in, void* user_data)
{
  BL_PROFILE("Pele::ReactorArkode::cF_RHS()");
#ifdef AMREX_USE_GPU
  sunrealtype* yvec_d = N_VGetDeviceArrayPointer(y_in);
  sunrealtype* ydot_d = N_VGetDeviceArrayPointer(ydot_in);
#else
  sunrealtype* yvec_d = N_VGetArrayPointer(y_in);
  sunrealtype* ydot_d = N_VGetArrayPointer(ydot_in);
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
  long int nst, nst_a, netf, nfe;
  int flag;

  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  utils::check_flag(&flag, "ARKodeGetNumSteps", 1);
  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  utils::check_flag(&flag, "ARKodeGetNumStepAttempts", 1);
  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  utils::check_flag(&flag, "ARKodeGetNumErrTestFails", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  utils::check_flag(&flag, "ARKodeGetNumRhsEvals", 1);

#ifdef AMREX_USE_OMP
  amrex::Print() << "\nFinal Statistics: "
                 << "(thread:" << omp_get_thread_num() << ", ";
  amrex::Print() << "arkodeMem:" << arkode_mem << ")\n";
#else
  amrex::Print() << "\nFinal Statistics:\n";
#endif

  amrex::Print() << "   Internal steps   = " << nst << "\n";
  amrex::Print() << "   Attempted steps  = " << nst_a << "\n";
  amrex::Print() << "   Error test fails = " << netf << "\n";
  amrex::Print() << "   Total RHS evals  = " << nfe << "\n";
}
} // namespace pele::physics::reactions
