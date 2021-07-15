#include "reactor.H"

#ifdef AMREX_USE_GPU

void
SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals)
{
}

void
SetTolFactODE(amrex::Real relative_tol, amrex::Real absolute_tol)
{
}

int
reactor_info(int /*reactor_type*/, int /*ode_ncells*/)
{
  return (0);
}

// Array4 version
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
  const int& reactor_type,
  amrex::gpuStream_t stream)
{

  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  constexpr int eint_rho = 1;
  constexpr int enth_rho = 2;

  amrex::Real time_init = time;

  ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real renergy_loc =
      rEner_in(i, j, k, 0) + rEner_src_in(i, j, k, 0) * dt_react;
    rEner_in(i, j, k, 0) = renergy_loc;
    amrex::Real rY_loc[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      rY_loc[n] = rY_in(i, j, k, n) + rY_src_in(i, j, k, n) * dt_react;
      rY_in(i, j, k, n) = rY_loc[n];
    }
    amrex::Real rho_loc = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      rho_loc += rY_loc[n];
    }
    amrex::Real Y_loc[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      Y_loc[n] = rY_loc[n] / rho_loc;
    }
    amrex::Real energy_loc = renergy_loc / rho_loc;
    amrex::Real T_loc = T_in(i, j, k, 0);
    auto eos = pele::physics::PhysicsType::eos();
    if (reactor_type == eint_rho) {
      eos.EY2T(energy_loc, Y_loc, T_loc);
    } else {
      eos.HY2T(energy_loc, Y_loc, T_loc);
    }
    T_in(i, j, k, 0) = T_loc;
    FC_in(i, j, k, 0) = 0.0;
  });

#ifdef MOD_REACTOR
  time = time_init + dt_react;
#endif
  return 0;
}

// 1D version
int
react(
  amrex::Real* rY_in,
  amrex::Real* rY_src_in,
  amrex::Real* rX_in,
  amrex::Real* rX_src_in,
  amrex::Real& dt_react,
  amrex::Real& time,
  int cvode_iE,
  int Ncells,
  gpuStream_t stream)
{
  amrex::Real time_init = time;
#ifdef MOD_REACTOR
  time = time_init + dt_react;
#endif
  return 0;
}

void
reactor_close()
{
}

#else

int ireactor_type;

void
SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals)
{
}

void
SetTolFactODE(amrex::Real relative_tol, amrex::Real absolute_tol)
{
}

void
ReSetTolODE()
{
}

int
reactor_init(int reactor_type, int /*ode_ncells*/)
{
  ireactor_type = reactor_type;
  return (0);
}

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
  const int& reactor_type)
{
  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  constexpr int eint_rho = 1;
  constexpr int enth_rho = 2;

  amrex::Real time_init = time;

  ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real renergy_loc =
      rEner_in(i, j, k, 0) + rEner_src_in(i, j, k, 0) * dt_react;
    rEner_in(i, j, k, 0) = renergy_loc;
    amrex::Real rY_loc[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      rY_loc[n] = rY_in(i, j, k, n) + rY_src_in(i, j, k, n) * dt_react;
      rY_in(i, j, k, n) = rY_loc[n];
    }
    amrex::Real rho_loc = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      rho_loc += rY_loc[n];
    }
    amrex::Real Y_loc[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      Y_loc[n] = rY_loc[n] / rho_loc;
    }
    amrex::Real energy_loc = renergy_loc / rho_loc;
    amrex::Real T_loc = T_in(i, j, k, 0);
    auto eos = pele::physics::PhysicsType::eos();
    if (reactor_type == eint_rho) {
      eos.EY2T(energy_loc, Y_loc, T_loc);
    } else {
      eos.HY2T(energy_loc, Y_loc, T_loc);
    }
    T_in(i, j, k, 0) = T_loc;
    FC_in(i, j, k, 0) = 0.0;
  });

#ifdef MOD_REACTOR
  time = time_init + dt_react;
#endif
  return 0;
}

int
react(
  amrex::Real* rY_in,
  amrex::Real* rY_src_in,
  amrex::Real* rX_in,
  amrex::Real* rX_src_in,
  amrex::Real& dt_react,
  amrex::Real& time)
{
  amrex::Real time_init = time;
#ifdef MOD_REACTOR
  time = time_init + dt_react;
#endif
  return 0;
}

void
reactor_close()
{
}

#endif
