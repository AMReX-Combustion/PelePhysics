#include "ReactorNull.H"

namespace pele::physics::reactions {

int
ReactorNull::init(int reactor_type, int /*ncells*/)
{
  BL_PROFILE("Pele::ReactorNull::init()");
  m_reactor_type = reactor_type;
  ReactorTypes::check_reactor_type(m_reactor_type);
  amrex::ParmParse pp("ode");
  pp.query("verbose", verbose);
  return (0);
}

// Array4 version
int
ReactorNull::react(
  const amrex::Box& box,
  amrex::Array4<amrex::Real> const& rY_in,
  amrex::Array4<amrex::Real> const& rYsrc_in,
  amrex::Array4<amrex::Real> const& T_in,
  amrex::Array4<amrex::Real> const& rEner_in,
  amrex::Array4<amrex::Real> const& rEner_src_in,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Array4<int> const& /*mask*/,
  amrex::Real& dt_react,
  amrex::Real&
#ifdef MOD_REACTOR
    time
#endif
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t /*stream*/
#endif
)
{
#ifdef MOD_REACTOR
  amrex::Real time_init = time;
#endif

  int captured_reactor_type = m_reactor_type;
  const auto* leosparm = m_d_eosparm;

  ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real renergy_loc =
      rEner_in(i, j, k, 0) + rEner_src_in(i, j, k, 0) * dt_react;
    rEner_in(i, j, k, 0) = renergy_loc;
    amrex::Real rY_loc[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      rY_loc[n] = rY_in(i, j, k, n) + rYsrc_in(i, j, k, n) * dt_react;
      rY_in(i, j, k, n) = rY_loc[n];
    }
    amrex::Real rho_loc = 0.0;
    for (amrex::Real n : rY_loc) {
      rho_loc += n;
    }
    amrex::Real Y_loc[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      Y_loc[n] = rY_loc[n] / rho_loc;
    }
    amrex::Real energy_loc = renergy_loc / rho_loc;
    amrex::Real T_loc = T_in(i, j, k, 0);
    auto eos = pele::physics::PhysicsType::eos(leosparm);
    if (captured_reactor_type == ReactorTypes::e_reactor_type) {
      eos.REY2T(rho_loc, energy_loc, Y_loc, T_loc);
    } else if (captured_reactor_type == ReactorTypes::h_reactor_type) {
      eos.RHY2T(rho_loc, energy_loc, Y_loc, T_loc);
    } else {
      amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
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
ReactorNull::react(
  amrex::Real* /*rY_in*/,
  amrex::Real* /*rYsrc_in*/,
  amrex::Real* /*rX_in*/,
  amrex::Real* /*rX_src_in*/,
  amrex::Real&
#ifdef MOD_REACTOR
    dt_react
#endif
  /*dt_react*/,
  amrex::Real&
#ifdef MOD_REACTOR
    time
#endif
  /*time*/,
  int /*ncells*/
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t /*stream*/
#endif
)
{
#ifdef MOD_REACTOR
  amrex::Real time_init = time;
  time = time_init + dt_react;
#endif
  return 0;
}
} // namespace pele::physics::reactions
