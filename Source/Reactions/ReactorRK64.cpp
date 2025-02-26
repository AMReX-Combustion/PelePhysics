#include "AMReX_Reduce.H"
#include "ReactorRK64.H"

namespace pele::physics::reactions {

int
ReactorRK64::init(int reactor_type, int /*ncells*/)
{
  BL_PROFILE("Pele::ReactorRK64::init()");
  m_reactor_type = reactor_type;
  ReactorTypes::check_reactor_type(m_reactor_type);
  amrex::ParmParse pp("ode");
  pp.query("verbose", verbose);
  pp.query("atol", absTol);
  pp.query("rk64_nsubsteps_guess", rk64_nsubsteps_guess);
  pp.query("rk64_nsubsteps_min", rk64_nsubsteps_min);
  pp.query("rk64_nsubsteps_max", rk64_nsubsteps_max);
  pp.query("clean_init_massfrac", m_clean_init_massfrac);
  return (0);
}

int
ReactorRK64::react(
  amrex::Real* rY_in,
  amrex::Real* rYsrc_in,
  amrex::Real* rX_in,
  amrex::Real* rX_src_in,
  amrex::Real& dt_react,
  amrex::Real& time,
  int ncells
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t /*stream*/
#endif
)
{
  BL_PROFILE("Pele::ReactorRK64::react()");

  amrex::Real time_init = time;
  amrex::Real time_out = time + dt_react;
  const amrex::Real tinyval = 1e-50;

  // Copy to device
  amrex::Gpu::DeviceVector<amrex::Real> rY(ncells * (NUM_SPECIES + 1), 0);
  amrex::Gpu::DeviceVector<amrex::Real> rYsrc(ncells * NUM_SPECIES, 0);
  amrex::Gpu::DeviceVector<amrex::Real> rX(ncells, 0);
  amrex::Gpu::DeviceVector<amrex::Real> rX_src(ncells, 0);
  amrex::Real* d_rY = rY.data();
  amrex::Real* d_rYsrc = rYsrc.data();
  amrex::Real* d_rX = rX.data();
  amrex::Real* d_rX_src = rX_src.data();
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, rY_in, rY_in + ncells * (NUM_SPECIES + 1), d_rY);
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, rYsrc_in, rYsrc_in + ncells * NUM_SPECIES,
    d_rYsrc);
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, rX_in, rX_in + ncells, d_rX);
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, rX_src_in, rX_src_in + ncells, d_rX_src);

  // capture reactor type
  const int captured_reactor_type = m_reactor_type;
  const int captured_nsubsteps_guess = rk64_nsubsteps_guess;
  const int captured_nsubsteps_min = rk64_nsubsteps_min;
  const int captured_nsubsteps_max = rk64_nsubsteps_max;
  const amrex::Real captured_abstol = absTol;
  const auto* leosparm = m_d_eosparm;
  RK64Params rkp;

  amrex::Gpu::DeviceVector<int> v_nsteps(ncells, 0);
  int* d_nsteps = v_nsteps.data();

  amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(int icell) noexcept {
    amrex::Real soln_reg[NUM_SPECIES + 1] = {0.0};
    amrex::Real carryover_reg[NUM_SPECIES + 1] = {0.0};
    amrex::Real error_reg[NUM_SPECIES + 1] = {0.0};
    amrex::Real ydot[NUM_SPECIES + 1] = {0.0};
    amrex::Real rYsrc_ext[NUM_SPECIES] = {0.0};
    amrex::Real current_time = time_init;
    const int neq = (NUM_SPECIES + 1);

    for (int sp = 0; sp < neq; sp++) {
      soln_reg[sp] = d_rY[icell * neq + sp];
      carryover_reg[sp] = soln_reg[sp];
    }

    amrex::Real dt_rk = dt_react / amrex::Real(captured_nsubsteps_guess);
    amrex::Real dt_rk_min = dt_react / amrex::Real(captured_nsubsteps_max);
    amrex::Real dt_rk_max = dt_react / amrex::Real(captured_nsubsteps_min);

    amrex::Real rhoe_init[] = {d_rX[icell]};
    amrex::Real rhoesrc_ext[] = {d_rX_src[icell]};

    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      rYsrc_ext[sp] = d_rYsrc[icell * NUM_SPECIES + sp];
    }

    int nsteps = 0;
    amrex::Real change_factor;
    while (current_time < time_out) {
      for (amrex::Real& sp : error_reg) {
        sp = 0.0;
      }
      for (int stage = 0; stage < rkp.nstages_rk64; stage++) {
        utils::fKernelSpec<Ordering>(
          0, 1, current_time - time_init, captured_reactor_type, soln_reg, ydot,
          rhoe_init, rhoesrc_ext, rYsrc_ext, leosparm);

        for (int sp = 0; sp < neq; sp++) {
          error_reg[sp] += rkp.err_rk64[stage] * dt_rk * ydot[sp];
          soln_reg[sp] =
            carryover_reg[sp] + rkp.alpha_rk64[stage] * dt_rk * ydot[sp];
          carryover_reg[sp] =
            soln_reg[sp] + rkp.beta_rk64[stage] * dt_rk * ydot[sp];
        }
      }

      current_time += dt_rk;
      nsteps++;

      amrex::Real max_err = tinyval;
      for (amrex::Real sp : error_reg) {
        max_err = fabs(sp) > max_err ? fabs(sp) : max_err;
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
      // Don't overstep the integration time
      dt_rk = amrex::min<amrex::Real>(dt_rk, time_out - current_time);
    }
    d_nsteps[icell] = nsteps;
    // copy data back
    for (int sp = 0; sp < neq; sp++) {
      d_rY[icell * neq + sp] = soln_reg[sp];
    }
    d_rX[icell] = rhoe_init[0] + dt_react * rhoesrc_ext[0];
  });

#ifdef MOD_REACTOR
  time = time_out;
#endif

  const int avgsteps = amrex::Reduce::Sum<int>(
    ncells, [=] AMREX_GPU_DEVICE(int i) noexcept -> int { return d_nsteps[i]; },
    0);

  amrex::Gpu::copy(
    amrex::Gpu::deviceToHost, d_rY, d_rY + ncells * (NUM_SPECIES + 1), rY_in);
  amrex::Gpu::copy(
    amrex::Gpu::deviceToHost, d_rYsrc, d_rYsrc + ncells * NUM_SPECIES,
    rYsrc_in);
  amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_rX, d_rX + ncells, rX_in);
  amrex::Gpu::copy(
    amrex::Gpu::deviceToHost, d_rX_src, d_rX_src + ncells, rX_src_in);

  return (int(avgsteps / amrex::Real(ncells)));
}

int
ReactorRK64::react(
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
  amrex::gpuStream_t /*stream*/
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
  const auto* leosparm = m_d_eosparm;
  RK64Params rkp;

  int ncells = static_cast<int>(box.numPts());
  const auto len = amrex::length(box);
  const auto lo = amrex::lbound(box);

  amrex::Gpu::DeviceVector<int> v_nsteps(ncells, 0);
  int* d_nsteps = v_nsteps.data();

  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real soln_reg[NUM_SPECIES + 1] = {0.0};
    amrex::Real carryover_reg[NUM_SPECIES + 1] = {0.0};
    amrex::Real error_reg[NUM_SPECIES + 1] = {0.0};
    amrex::Real ydot[NUM_SPECIES + 1] = {0.0};
    amrex::Real rYsrc_ext[NUM_SPECIES] = {0.0};
    amrex::Real current_time = time_init;
    const int neq = (NUM_SPECIES + 1);

    auto eos = pele::physics::PhysicsType::eos(leosparm);
    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      soln_reg[sp] = rY_in(i, j, k, sp);
      carryover_reg[sp] = soln_reg[sp];
    }
    amrex::Real rho = 0.0, rho_inv = 0.0;
    amrex::Real mass_frac[NUM_SPECIES] = {0.0};
    eos.RY2RRinvY(soln_reg, rho, rho_inv, mass_frac);

    amrex::Real temp = T_in(i, j, k, 0);

    amrex::Real Enrg_loc = rEner_in(i, j, k, 0) * rho_inv;
    if (captured_reactor_type == ReactorTypes::e_reactor_type) {
      eos.REY2T(rho, Enrg_loc, mass_frac, temp);
    } else if (captured_reactor_type == ReactorTypes::h_reactor_type) {
      eos.RHY2T(rho, Enrg_loc, mass_frac, temp);
    } else {
      amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
    }
    soln_reg[NUM_SPECIES] = temp;
    carryover_reg[NUM_SPECIES] = soln_reg[NUM_SPECIES];

    amrex::Real dt_rk = dt_react / amrex::Real(captured_nsubsteps_guess);
    amrex::Real dt_rk_min = dt_react / amrex::Real(captured_nsubsteps_max);
    amrex::Real dt_rk_max = dt_react / amrex::Real(captured_nsubsteps_min);

    amrex::Real rhoe_init[] = {rEner_in(i, j, k, 0)};
    amrex::Real rhoesrc_ext[] = {rEner_src_in(i, j, k, 0)};

    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      rYsrc_ext[sp] = rYsrc_in(i, j, k, sp);
    }

    int nsteps = 0;
    amrex::Real change_factor;
    while (current_time < time_out) {
      for (amrex::Real& sp : error_reg) {
        sp = 0.0;
      }
      for (int stage = 0; stage < rkp.nstages_rk64; stage++) {
        utils::fKernelSpec<Ordering>(
          0, 1, current_time - time_init, captured_reactor_type, soln_reg, ydot,
          rhoe_init, rhoesrc_ext, rYsrc_ext, leosparm);

        for (int sp = 0; sp < neq; sp++) {
          error_reg[sp] += rkp.err_rk64[stage] * dt_rk * ydot[sp];
          soln_reg[sp] =
            carryover_reg[sp] + rkp.alpha_rk64[stage] * dt_rk * ydot[sp];
          carryover_reg[sp] =
            soln_reg[sp] + rkp.beta_rk64[stage] * dt_rk * ydot[sp];
        }
      }

      current_time += dt_rk;
      nsteps++;

      amrex::Real max_err = tinyval;
      for (amrex::Real sp : error_reg) {
        max_err = fabs(sp) > max_err ? fabs(sp) : max_err;
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
      // Don't overstep the integration time
      dt_rk = amrex::min<amrex::Real>(dt_rk, time_out - current_time);
    }

    // copy data back
    int icell = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
    d_nsteps[icell] = nsteps;
    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      rY_in(i, j, k, sp) = soln_reg[sp];
    }
    eos.RY2RRinvY(soln_reg, rho, rho_inv, mass_frac);

    temp = soln_reg[NUM_SPECIES];
    rEner_in(i, j, k, 0) = rhoe_init[0] + dt_react * rhoesrc_ext[0];
    Enrg_loc = rEner_in(i, j, k, 0) * rho_inv;

    if (captured_reactor_type == ReactorTypes::e_reactor_type) {
      eos.REY2T(rho, Enrg_loc, mass_frac, temp);
    } else if (captured_reactor_type == ReactorTypes::h_reactor_type) {
      eos.RHY2T(rho, Enrg_loc, mass_frac, temp);
    } else {
      amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
    }
    T_in(i, j, k, 0) = temp;
    FC_in(i, j, k, 0) = nsteps;
  });

#ifdef MOD_REACTOR
  time = time_out;
#endif

  const int avgsteps = amrex::Reduce::Sum<int>(
    ncells, [=] AMREX_GPU_DEVICE(int i) noexcept -> int { return d_nsteps[i]; },
    0);
  return (int(avgsteps / amrex::Real(ncells)));
}

} // namespace pele::physics::reactions
