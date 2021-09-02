#include "ReactorRK64.H"

namespace pele {
namespace physics {
namespace reactions {

int
ReactorRK64::init(int reactor_type, int /*Ncells*/)
{
  BL_PROFILE("Pele::ReactorRK64::init()");
  m_reactor_type = reactor_type;
  amrex::ParmParse pp("ode");
  pp.query("atol", absTol);
  pp.query("rk64_nsubsteps_guess", rk64_nsubsteps_guess);
  pp.query("rk64_nsubsteps_min", rk64_nsubsteps_min);
  pp.query("rk64_nsubsteps_max", rk64_nsubsteps_max);
  return (0);
}

int
ReactorRK64::react(
  amrex::Real* rY_in,
  amrex::Real* rY_src_in,
  amrex::Real* rX_in,
  amrex::Real* rX_src_in,
  amrex::Real& dt_react,
  amrex::Real& time,
  int Ncells
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  BL_PROFILE("Pele::ReactorRK64::react()");

  amrex::Real time_init = time;
  amrex::Real time_out = time + dt_react;
  const amrex::Real tinyval = 1e-50;

  // capture reactor type
  const int captured_reactor_type = m_reactor_type;
  const int captured_nsubsteps_guess = rk64_nsubsteps_guess;
  const int captured_nsubsteps_min = rk64_nsubsteps_min;
  const int captured_nsubsteps_max = rk64_nsubsteps_max;
  const amrex::Real captured_abstol = absTol;
  RK64Params rkp;

  int* nstepsvec;
  nstepsvec = new int[Ncells]();

  amrex::ParallelFor(Ncells, [=] AMREX_GPU_DEVICE(int icell) noexcept {
    amrex::Real soln_reg[NUM_SPECIES + 1] = {0.0};
    amrex::Real carryover_reg[NUM_SPECIES + 1] = {0.0};
    amrex::Real error_reg[NUM_SPECIES + 1] = {0.0};
    amrex::Real rhs[NUM_SPECIES + 1] = {0.0};
    amrex::Real rYsrc[NUM_SPECIES] = {0.0};
    amrex::Real current_time = time_init;
    const int neq = (NUM_SPECIES + 1);

    for (int sp = 0; sp < neq; sp++) {
      soln_reg[sp] = rY_in[icell * neq + sp];
      carryover_reg[sp] = soln_reg[sp];
    }

    amrex::Real dt_rk = dt_react / amrex::Real(captured_nsubsteps_guess);
    amrex::Real dt_rk_min = dt_react / amrex::Real(captured_nsubsteps_max);
    amrex::Real dt_rk_max = dt_react / amrex::Real(captured_nsubsteps_min);

    amrex::Real rhoe_init[] = {rX_in[icell]};
    amrex::Real rhoesrc_ext[] = {rX_src_in[icell]};

    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      rYsrc[sp] = rY_src_in[icell * NUM_SPECIES + sp];
    }

    int nsteps = 0;
    amrex::Real change_factor;
    while (current_time < time_out) {
      for (int sp = 0; sp < neq; sp++) {
        error_reg[sp] = 0.0;
      }
      for (int stage = 0; stage < rkp.nstages_rk64; stage++) {
        fKernelSpec<Ordering>(
          0, 1, current_time - time_init, captured_reactor_type, soln_reg, rhs,
          rhoe_init, rhoesrc_ext, rYsrc);

        for (int sp = 0; sp < neq; sp++) {
          error_reg[sp] += rkp.err_rk64[stage] * dt_rk * rhs[sp];
          soln_reg[sp] =
            carryover_reg[sp] + rkp.alpha_rk64[stage] * dt_rk * rhs[sp];
          carryover_reg[sp] =
            soln_reg[sp] + rkp.beta_rk64[stage] * dt_rk * rhs[sp];
        }
      }

      current_time += dt_rk;
      nsteps++;

      amrex::Real max_err = tinyval;
      for (int sp = 0; sp < neq; sp++) {
        if (fabs(error_reg[sp]) > max_err) {
          max_err = fabs(error_reg[sp]);
        }
      }

      if (max_err < captured_abstol) {
        change_factor =
          rkp.betaerr_rk64 * pow((captured_abstol / max_err), rkp.exp1_rk64);
        dt_rk = amrex::min<amrex::Real>(dt_rk_max, dt_rk * change_factor);
      } else {
        change_factor =
          rkp.betaerr_rk64 * pow((captured_abstol / max_err), rkp.exp2_rk64);
        dt_rk = amrex::max<amrex::Real>(dt_rk_min, dt_rk * change_factor);
      }
    }
    nstepsvec[icell] = nsteps;
    // copy data back
    for (int sp = 0; sp < neq; sp++) {
      rY_in[icell * neq + sp] = soln_reg[sp];
    }
    rX_in[icell] = rhoe_init[0] + dt_react * rhoesrc_ext[0];
  });

#ifdef MOD_REACTOR
  time = time_out;
#endif

  amrex::Real avgsteps = 0.0;
  for (int i = 0; i < Ncells; i++) {
    avgsteps += nstepsvec[i];
  }

  return (int(avgsteps / amrex::Real(Ncells)));
}

int
ReactorRK64::react(
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
  BL_PROFILE("Pele::ReactorRK64::react()");

  amrex::Real time_init = time;
  amrex::Real time_out = time + dt_react;
  const amrex::Real tinyval = 1e-50;

  // capture reactor type
  const int captured_reactor_type = m_reactor_type;
  const int captured_nsubsteps_guess = rk64_nsubsteps_guess;
  const int captured_nsubsteps_min = rk64_nsubsteps_min;
  const int captured_nsubsteps_max = rk64_nsubsteps_max;
  const amrex::Real captured_abstol = absTol;
  RK64Params rkp;

  int* nstepsvec;
  int Ncells = box.numPts();
  const auto len = amrex::length(box);
  const auto lo = amrex::lbound(box);
  nstepsvec = new int[Ncells]();

  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real soln_reg[NUM_SPECIES + 1] = {0.0};
    amrex::Real carryover_reg[NUM_SPECIES + 1] = {0.0};
    amrex::Real error_reg[NUM_SPECIES + 1] = {0.0};
    amrex::Real rhs[NUM_SPECIES + 1] = {0.0};
    amrex::Real rYsrc[NUM_SPECIES] = {0.0};
    amrex::Real current_time = time_init;
    const int neq = (NUM_SPECIES + 1);

    amrex::Real rho = 0.0;
    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      soln_reg[sp] = rY_in(i, j, k, sp);
      carryover_reg[sp] = soln_reg[sp];
      rho += rY_in(i, j, k, sp);
    }
    amrex::Real rho_inv = 1.0 / rho;
    amrex::Real mass_frac[NUM_SPECIES] = {0.0};
    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      mass_frac[sp] = rY_in(i, j, k, sp) * rho_inv;
    }
    amrex::Real temp = T_in(i, j, k, 0);

    amrex::Real Enrg_loc = rEner_in(i, j, k, 0) * rho_inv;
    auto eos = pele::physics::PhysicsType::eos();
    if (captured_reactor_type == 1) {
      eos.REY2T(rho, Enrg_loc, mass_frac, temp);
    } else {
      eos.RHY2T(rho, Enrg_loc, mass_frac, temp);
    }
    soln_reg[NUM_SPECIES] = temp;
    carryover_reg[NUM_SPECIES] = soln_reg[NUM_SPECIES];

    amrex::Real dt_rk = dt_react / amrex::Real(captured_nsubsteps_guess);
    amrex::Real dt_rk_min = dt_react / amrex::Real(captured_nsubsteps_max);
    amrex::Real dt_rk_max = dt_react / amrex::Real(captured_nsubsteps_min);

    amrex::Real rhoe_init[] = {rEner_in(i, j, k, 0)};
    amrex::Real rhoesrc_ext[] = {rEner_src_in(i, j, k, 0)};

    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      rYsrc[sp] = rY_src_in(i, j, k, sp);
    }

    int nsteps = 0;
    amrex::Real change_factor;
    while (current_time < time_out) {
      for (int sp = 0; sp < neq; sp++) {
        error_reg[sp] = 0.0;
      }
      for (int stage = 0; stage < rkp.nstages_rk64; stage++) {
        fKernelSpec<Ordering>(
          0, 1, current_time - time_init, captured_reactor_type, soln_reg, rhs,
          rhoe_init, rhoesrc_ext, rYsrc);

        for (int sp = 0; sp < neq; sp++) {
          error_reg[sp] += rkp.err_rk64[stage] * dt_rk * rhs[sp];
          soln_reg[sp] =
            carryover_reg[sp] + rkp.alpha_rk64[stage] * dt_rk * rhs[sp];
          carryover_reg[sp] =
            soln_reg[sp] + rkp.beta_rk64[stage] * dt_rk * rhs[sp];
        }
      }

      current_time += dt_rk;
      nsteps++;

      amrex::Real max_err = tinyval;
      for (int sp = 0; sp < neq; sp++) {
        if (fabs(error_reg[sp]) > max_err) {
          max_err = fabs(error_reg[sp]);
        }
      }

      if (max_err < captured_abstol) {
        change_factor =
          rkp.betaerr_rk64 * pow((captured_abstol / max_err), rkp.exp1_rk64);
        dt_rk = amrex::min<amrex::Real>(dt_rk_max, dt_rk * change_factor);
      } else {
        change_factor =
          rkp.betaerr_rk64 * pow((captured_abstol / max_err), rkp.exp2_rk64);
        dt_rk = amrex::max<amrex::Real>(dt_rk_min, dt_rk * change_factor);
      }
    }

    // copy data back
    int icell = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
    nstepsvec[icell] = nsteps;
    rho = 0.0;
    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      rY_in(i, j, k, sp) = soln_reg[sp];
      rho += rY_in(i, j, k, sp);
    }
    rho_inv = 1.0 / rho;
    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      mass_frac[sp] = rY_in(i, j, k, sp) * rho_inv;
    }
    temp = soln_reg[NUM_SPECIES];
    rEner_in(i, j, k, 0) = rhoe_init[0] + dt_react * rhoesrc_ext[0];
    Enrg_loc = rEner_in(i, j, k, 0) * rho_inv;

    if (captured_reactor_type == 1) {
      eos.REY2T(rho, Enrg_loc, mass_frac, temp);
    } else {
      eos.RHY2T(rho, Enrg_loc, mass_frac, temp);
    }
    T_in(i, j, k, 0) = temp;
    FC_in(i, j, k, 0) = nsteps;
  });

#ifdef MOD_REACTOR
  time = time_out;
#endif

  amrex::Real avgsteps = 0.0;
  for (int i = 0; i < Ncells; i++) {
    avgsteps += nstepsvec[i];
  }

  return (int(avgsteps / amrex::Real(Ncells)));
}

void
ReactorRK64::SetTypValsODE(const std::vector<amrex::Real>& /*ExtTypVals*/)
{
  amrex::Print() << "WARNING: ignoring TypVals for this reactor." << std::endl;
}
} // namespace reactions
} // namespace physics
} // namespace pele
