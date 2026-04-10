/**
 * @file ReactorQSS.cpp
 * @brief QSS chemistry reactor implementation for PeleLMeX.
 *
 * Based on CHEMEQ2 (Mott, Oran & van Leer, JCP 164(2), 2000).
 * CPU-only: QssOde virtual dispatch is not GPU-callable.
 */

#include "ReactorQSS.H"
#include "ReactorTypes.H"
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <cmath>

namespace pele::physics::reactions {

// ============================================================
// init
// ============================================================

int
ReactorQSS::init(int reactor_type, int /*ncells*/)
{
  BL_PROFILE("Pele::ReactorQSS::init()");

  m_reactor_type = reactor_type;
  ReactorTypes::check_reactor_type(m_reactor_type);

  // Common ODE parameters (shared prefix with other reactors)
  {
    amrex::ParmParse pp("ode");
    pp.query("verbose", verbose);
    pp.query("clean_init_massfrac", m_clean_init_massfrac);
  }

  // QSS-specific parameters
  {
    amrex::ParmParse pp("qss");
    pp.query("epsmax",          m_epsmax);
    pp.query("epsmin",          m_epsmin);
    pp.query("dtmin",           m_dtmin);
    pp.query("dtmax",           m_dtmax);
    pp.query("itermax",         m_itermax);
    pp.query("tfd",             m_tfd);
    pp.query("abstol",          m_abstol);
    pp.query("Tmax",            m_Tmax);
    int stab_check = m_stability_check ? 1 : 0;
    pp.query("stability_check", stab_check);
    m_stability_check = (stab_check != 0);
  }

  if (verbose > 0) {
    amrex::Print() << "ReactorQSS::init()\n"
                   << "  reactor_type  = " << m_reactor_type << "\n"
                   << "  epsmax        = " << m_epsmax << "\n"
                   << "  epsmin        = " << m_epsmin << "\n"
                   << "  dtmin         = " << m_dtmin << "\n"
                   << "  dtmax         = " << m_dtmax << "\n"
                   << "  itermax       = " << m_itermax << "\n"
                   << "  abstol        = " << m_abstol << "\n"
                   << "  stability_check = " << m_stability_check << "\n";
  }

  return 0;
}

// ============================================================
// CellOde::setup
// ============================================================

void
ReactorQSS::CellOde::setup(
  int reactor_type,
  const amrex::Real rYsrc[NUM_SPECIES],
  amrex::Real rhoesrc_ext,
  amrex::Real tstart,
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    eosparm) noexcept
{
  m_reactor_type = reactor_type;
  for (int n = 0; n < NUM_SPECIES; n++) {
    m_rYsrc[n] = rYsrc[n];
  }
  m_rhoesrc_ext = rhoesrc_ext;
  m_tstart      = tstart;
  m_eosparm     = eosparm;
  m_nfeval      = 0;
}

// ============================================================
// CellOde::odefun — sign-split chemistry RHS for QSS
//
// State: y[0..NUM_SPECIES-1] = rhoY_n [g/cm^3]
//        y[NUM_SPECIES]       = T      [K]
//
// Returns q >= 0 (production) and d >= 0 (destruction) such that
//   dy/dt = q - d
//
// Species split: q_i = max(0, wdot_i + rYsrc_i)
//                d_i = max(0, -(wdot_i + rYsrc_i))
//
// Temperature:   Tdot = (rhoesrc_ext - sum_i((wdot_i+rYsrc_i)*ei)) / (rho*CvCp)
//                q_T  = max(0, Tdot)
//                d_T  = max(0, -Tdot)
// ============================================================

void
ReactorQSS::CellOde::odefun(
  amrex::Real /*t*/,
  const amrex::Vector<amrex::Real>& y,
  amrex::Vector<amrex::Real>& q,
  amrex::Vector<amrex::Real>& d,
  bool /*corrector*/)
{
  AMREX_ASSERT(static_cast<int>(y.size()) == NUM_SPECIES + 1);
  ++m_nfeval;

  // ------------------------------------------------------------------
  // Unpack state: rhoY_n and T
  // ------------------------------------------------------------------
  amrex::Real massdens[NUM_SPECIES];
  amrex::Real rho = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    massdens[n] = y[n];
    rho += massdens[n];
  }
  const amrex::Real rho_inv = 1.0 / rho;

  amrex::Real massfrac[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = massdens[n] * rho_inv;
  }

  // T is bounded to [ymin, ymax] = [250, 4000] K by the QssIntegrator before
  // every odefun call, so it is always in a valid range for the EOS.
  const amrex::Real T = y[NUM_SPECIES];

  // ------------------------------------------------------------------
  // Chemistry: net production rates wdot [g/cm^3/s]
  // ------------------------------------------------------------------
  auto eos = pele::physics::PhysicsType::eos(m_eosparm);

  amrex::Real wdot[NUM_SPECIES] = {0.0};
  eos.RTY2WDOT(rho, T, massfrac, wdot);

  // ------------------------------------------------------------------
  // Thermal properties for energy equation
  // ------------------------------------------------------------------
  amrex::Real ei[NUM_SPECIES] = {0.0}; // specific enthalpy or internal energy
  amrex::Real CvCp = 1.0;             // Cp (h reactor) or Cv (e reactor)

  if (m_reactor_type == ReactorTypes::e_reactor_type) {
    eos.RTY2Ei(rho, T, massfrac, ei);
    eos.RTY2Cv(rho, T, massfrac, CvCp);
  } else { // h_reactor_type
    eos.RTY2Hi(rho, T, massfrac, ei);
    eos.RTY2Cp(rho, T, massfrac, CvCp);
  }

  // ------------------------------------------------------------------
  // Species RHS: net chemistry + external source, then sign-split
  // Also accumulate the energy budget.
  // ------------------------------------------------------------------
  amrex::Real rhoesrc = m_rhoesrc_ext;
  for (int n = 0; n < NUM_SPECIES; n++) {
    const amrex::Real rhs_n = wdot[n] + m_rYsrc[n];
    q[n] = (rhs_n > 0.0) ? rhs_n : 0.0;
    d[n] = (rhs_n < 0.0) ? -rhs_n : 0.0;
    // Subtract the enthalpy/energy carried by this species source
    rhoesrc -= rhs_n * ei[n];
  }

  // ------------------------------------------------------------------
  // Temperature RHS: energy residual → dT/dt, then sign-split
  //
  // rhoesrc [erg/cm^3/s] is what's left after all species carry away
  // their enthalpy. Divided by rho*Cp (or Cv) gives K/s.
  // ------------------------------------------------------------------
  const amrex::Real Tdot = rhoesrc * rho_inv / CvCp;
  q[NUM_SPECIES] = (Tdot > 0.0) ? Tdot : 0.0;
  d[NUM_SPECIES] = (Tdot < 0.0) ? -Tdot : 0.0;
}

// ============================================================
// integrate_cells — CPU loop, one QssIntegrator per call
// ============================================================

void
ReactorQSS::integrate_cells(
  amrex::Real* y_vect,
  const amrex::Real* src_vect,
  const amrex::Real* src_vect_energy,
  amrex::Real dt_react,
  amrex::Real time,
  int ncells,
  long int* FCunt)
{
  const int neq = NUM_SPECIES + 1;

  // Build a fresh QssIntegrator and CellOde for this call (thread-safe).
  QssIntegrator integrator;
  integrator.initialize(neq);
  integrator.epsmax         = m_epsmax;
  integrator.epsmin         = m_epsmin;
  integrator.dtmin          = m_dtmin;
  integrator.dtmax          = m_dtmax;
  integrator.itermax        = m_itermax;
  integrator.tfd            = m_tfd;
  integrator.abstol         = m_abstol;
  integrator.stabilityCheck = m_stability_check;

  // Temperature bounds: keep T inside the valid range of the thermodynamic
  // polynomials at every predictor/corrector step.
  //   ymin  — prevents driving T to 1e-20 K (the species default) which would
  //           make RTY2WDOT produce NaN.
  //   ymax  — prevents the predictor from overshooting into ranges where EOS
  //           polynomials extrapolate wildly (NASA 7-term polynomials are only
  //           fitted to ~5000 K).  Without this bound the predictor can jump
  //           to 1e4+ K in one dt_max step, making eps explode and collapsing
  //           the step to below dtmin in the very first corrector iteration.
  // enforce_ymin stays at 1 (default) so CHEMEQ2 computes rtau for T and
  // damps temperature oscillations during ignition.
  integrator.ymin[NUM_SPECIES]         = 250.0;   // [K]
  integrator.ymax[NUM_SPECIES]         = m_Tmax;  // [K] set via qss.Tmax
  integrator.enforce_ymax[NUM_SPECIES] = 1.0;

  CellOde cell_ode;
  integrator.setOde(&cell_ode);

  amrex::Vector<amrex::Real> y_cell(neq);
  amrex::Vector<amrex::Real> y_out(neq);

  for (int icell = 0; icell < ncells; icell++) {
    // ----------------------------------------------------------------
    // Unpack flat state for this cell into y_cell
    // ----------------------------------------------------------------
    for (int n = 0; n < neq; n++) {
      y_cell[n] = y_vect[icell * neq + n];
    }

    // ----------------------------------------------------------------
    // Build per-cell external sources for CellOde
    // ----------------------------------------------------------------
    const amrex::Real* rYsrc_cell = src_vect + icell * NUM_SPECIES;
    const amrex::Real  rhoesrc    = src_vect_energy[icell];

    cell_ode.setup(m_reactor_type, rYsrc_cell, rhoesrc, time, m_h_eosparm);

    // ----------------------------------------------------------------
    // Integrate from time to time + dt_react
    // ----------------------------------------------------------------
    integrator.setState(y_cell, time);
    const int ret = integrator.integrateToTime(time + dt_react);

    if (ret != 0) {
      // Integration failed (step collapsed below dtmin).  Freeze this cell
      // at its initial state rather than accepting a potentially-NaN result.
      // The cell will be re-attempted next SDC iteration with updated flow
      // conditions.
      if (verbose > 0) {
        amrex::Print() << "ReactorQSS: integration failed for cell " << icell
                       << " (ret=" << ret << "). Freezing cell state.\n";
      }
      // y_vect already contains y_cell values; no write needed.
      FCunt[icell] = static_cast<long int>(cell_ode.nfeval());
      continue;
    }

    integrator.getState(y_out);

    // ----------------------------------------------------------------
    // Write integrated state back to flat vector
    // ----------------------------------------------------------------
    for (int n = 0; n < neq; n++) {
      y_vect[icell * neq + n] = y_out[n];
    }

    FCunt[icell] = static_cast<long int>(cell_ode.nfeval());
  }
}

// ============================================================
// react — Array4 version (primary PeleLMeX path)
// ============================================================

int
ReactorQSS::react(
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
    /*time*/
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t /*stream*/
#endif
)
{
  BL_PROFILE("Pele::ReactorQSS::react()");

#ifdef AMREX_USE_GPU
  // The QssOde interface uses virtual dispatch (host code) and amrex::Vector
  // (host memory). GPU execution is not supported.
  amrex::Abort(
    "ReactorQSS is CPU-only. GPU builds must use ReactorCvode instead.");
  return -1;
#else
  const int ncells = static_cast<int>(box.numPts());
  const int neq    = NUM_SPECIES + 1;

#ifdef MOD_REACTOR
  const amrex::Real time_init = time;
#else
  // Without MOD_REACTOR, the time argument is unused. We still need a
  // local start time for the QSS integrator state.  Use 0.0 as the
  // offset — only the interval length dt_react matters for chemistry.
  const amrex::Real time_init = 0.0;
#endif

  // ------------------------------------------------------------------
  // Allocate flat host arrays.
  // On CPU (no GPU), amrex::Gpu::DeviceVector<T> is amrex::Vector<T>,
  // so no actual device allocation happens here.
  // ------------------------------------------------------------------
  amrex::Gpu::DeviceVector<amrex::Real> d_y(ncells * neq);
  amrex::Gpu::DeviceVector<amrex::Real> d_src(ncells * NUM_SPECIES);
  amrex::Gpu::DeviceVector<amrex::Real> d_energy(ncells);
  amrex::Gpu::DeviceVector<amrex::Real> d_energy_src(ncells);
  amrex::Gpu::DeviceVector<long int>    d_fc(ncells, 0L);

  // ------------------------------------------------------------------
  // Flatten: Array4 → flat 1D vectors.
  // Recomputes T from energy+species for consistency (see box_flatten).
  // ------------------------------------------------------------------
  flatten_ops.flatten(
    box, ncells, m_reactor_type, m_clean_init_massfrac,
    rY_in, rYsrc_in, T_in, rEner_in, rEner_src_in,
    d_y.data(), d_src.data(), d_energy.data(), d_energy_src.data());

  // ------------------------------------------------------------------
  // QSS integration on the CPU.
  // d_y, d_src, d_energy_src are host memory on CPU builds.
  // ------------------------------------------------------------------
  integrate_cells(
    d_y.data(), d_src.data(), d_energy_src.data(),
    dt_react, time_init, ncells, d_fc.data());

  // ------------------------------------------------------------------
  // Unflatten: flat 1D → Array4.
  // Updates rEner_in = initial + dt * external_src, then recomputes T
  // from updated energy and integrated species for thermodynamic
  // consistency. The T from QSS is used only as the initial guess for
  // the Newton solve in REY2T / RHY2T.
  // ------------------------------------------------------------------
  flatten_ops.unflatten(
    box, ncells, m_reactor_type, m_clean_init_massfrac,
    rY_in, T_in, rEner_in, rEner_src_in, FC_in,
    d_y.data(), d_energy.data(), d_fc.data(), dt_react);

#ifdef MOD_REACTOR
  time = time_init + dt_react;
#endif

  return 0;
#endif // AMREX_USE_GPU
}

// ============================================================
// react — 1D raw-pointer version (testing / non-GPU fallback)
// ============================================================

int
ReactorQSS::react(
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
  BL_PROFILE("Pele::ReactorQSS::react()");

#ifdef AMREX_USE_GPU
  amrex::Abort(
    "ReactorQSS is CPU-only. GPU builds must use ReactorCvode instead.");
  return -1;
#else
  const amrex::Real time_init = time;
  const int neq = NUM_SPECIES + 1;

  // rY_in is already in YCOrder layout: [rhoY_0..rhoY_{N-1}, T]_{icell}
  // rX_src_in is the per-cell energy source.
  amrex::Vector<long int> FCunt(ncells, 0L);

  integrate_cells(
    rY_in, rYsrc_in, rX_src_in, dt_react, time_init, ncells, FCunt.data());

  // Update energy with external source (chemistry does not change total energy;
  // only external forcing does, matching what unflatten does in the Array4 path).
  const auto* leosparm = m_h_eosparm;
  const int captured_reactor_type = m_reactor_type;
  for (int icell = 0; icell < ncells; icell++) {
    rX_in[icell] += dt_react * rX_src_in[icell];

    // Recompute T from updated energy and species for consistency.
    amrex::Real massdens[NUM_SPECIES];
    amrex::Real rho = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      massdens[n] = rY_in[icell * neq + n];
      rho += massdens[n];
    }
    const amrex::Real rho_inv = 1.0 / rho;
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = massdens[n] * rho_inv;
    }
    amrex::Real Enrg = rX_in[icell] * rho_inv;
    amrex::Real T    = rY_in[icell * neq + NUM_SPECIES]; // QSS T as initial guess
    auto eos = pele::physics::PhysicsType::eos(leosparm);
    if (captured_reactor_type == ReactorTypes::e_reactor_type) {
      eos.REY2T(rho, Enrg, massfrac, T);
    } else {
      eos.RHY2T(rho, Enrg, massfrac, T);
    }
    rY_in[icell * neq + NUM_SPECIES] = T;
  }

#ifdef MOD_REACTOR
  time = time_init + dt_react;
#endif
  return 0;
#endif // AMREX_USE_GPU
}

} // namespace pele::physics::reactions
