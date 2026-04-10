#include "ReactorAdaptive.H"
#include "ReactorTypes.H"
#include "ReactorUtils.H"
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Loop.H>

namespace pele::physics::reactions {

// ============================================================
// init
// ============================================================

int
ReactorAdaptive::init(int reactor_type, int ncells)
{
  BL_PROFILE("Pele::ReactorAdaptive::init()");

  m_reactor_type = reactor_type;
  ReactorTypes::check_reactor_type(m_reactor_type);

  amrex::ParmParse pp("adaptive");
  pp.query("primary_integrator",  m_primary_name);
  pp.query("fallback_integrator", m_fallback_name);
  pp.query("verbose",             m_adaptive_verbose);
  pp.query("T_threshold",         m_T_threshold);
  pp.query("species_name",        m_spec_name);
  pp.query("Y_threshold",         m_Y_threshold);

  std::string cond_type_str  = "temperature";
  std::string cond_logic_str = "AND";
  pp.query("condition_type",  cond_type_str);
  pp.query("condition_logic", cond_logic_str);

  // Parse condition type
  if (cond_type_str == "temperature") {
    m_cond_type = ConditionType::Temperature;
  } else if (cond_type_str == "species") {
    m_cond_type = ConditionType::Species;
  } else if (cond_type_str == "combined") {
    m_cond_type = (cond_logic_str == "OR") ? ConditionType::CombinedOr
                                            : ConditionType::CombinedAnd;
  } else {
    amrex::Abort(
      "ReactorAdaptive: unknown condition_type '" + cond_type_str +
      "'. Choose: temperature | species | combined");
  }

  // Create and initialise sub-reactors.
  // Note: set_eos_parm() is called after init() by PeleLMeX; species index
  // resolution is deferred there.
  m_primary  = ReactorBase::create(m_primary_name);
  m_fallback = ReactorBase::create(m_fallback_name);

  m_primary->init(reactor_type, ncells);
  m_fallback->init(reactor_type, ncells);

  amrex::Print() << "Initializing ReactorAdaptive:\n"
                 << "  Primary  integrator : " << m_primary_name << "\n"
                 << "  Fallback integrator : " << m_fallback_name << "\n"
                 << "  Condition type      : " << cond_type_str;
  if (m_cond_type == ConditionType::CombinedAnd ||
      m_cond_type == ConditionType::CombinedOr) {
    amrex::Print() << " (" << cond_logic_str << ")";
  }
  amrex::Print() << "\n";
  if (m_cond_type != ConditionType::Species) {
    amrex::Print() << "  T_threshold         : " << m_T_threshold << " K\n";
  }
  if (m_cond_type != ConditionType::Temperature) {
    amrex::Print() << "  Species             : " << m_spec_name << "\n"
                   << "  Y_threshold         : " << m_Y_threshold << "\n";
  }
  amrex::Print()
    << "  (condition TRUE → primary; FALSE → fallback)\n";

  return 0;
}

// ============================================================
// set_eos_parm — propagate to sub-reactors and resolve species index
// ============================================================

void
ReactorAdaptive::set_eos_parm(
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    h_eosparm,
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    d_eosparm)
{
  // Store locally (via base class)
  ReactorBase::set_eos_parm(h_eosparm, d_eosparm);

  // Forward to both sub-reactors
  m_primary->set_eos_parm(h_eosparm, d_eosparm);
  m_fallback->set_eos_parm(h_eosparm, d_eosparm);

  // Resolve species index now that we have the EOS params
  if (m_cond_type != ConditionType::Temperature) {
    amrex::Vector<std::string> kname;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
      kname, m_h_eosparm);
    m_spec_idx = -1;
    for (int i = 0; i < NUM_SPECIES; i++) {
      if (kname[i] == m_spec_name) {
        m_spec_idx = i;
        break;
      }
    }
    if (m_spec_idx < 0) {
      amrex::Abort(
        "ReactorAdaptive: species '" + m_spec_name +
        "' not found in mechanism. Check adaptive.species_name.");
    }
    amrex::Print() << "ReactorAdaptive: species '" << m_spec_name
                   << "' resolved to index " << m_spec_idx << "\n";
  }
}

// ============================================================
// set_typ_vals_ode — propagate to both sub-reactors
// ============================================================

void
ReactorAdaptive::set_typ_vals_ode(const std::vector<amrex::Real>& ExtTypVals)
{
  ReactorBase::set_typ_vals_ode(ExtTypVals); // stores in m_typ_vals
  m_primary->set_typ_vals_ode(ExtTypVals);
  m_fallback->set_typ_vals_ode(ExtTypVals);
}

// ============================================================
// use_primary — per-cell condition evaluation
// ============================================================

bool
ReactorAdaptive::use_primary(
  amrex::Real T,
  amrex::Array4<amrex::Real> const& rY_in,
  amrex::Real rho,
  int i,
  int j,
  int k) const noexcept
{
  const bool cond_T = (T < m_T_threshold);

  bool cond_Y = true;
  if (m_cond_type != ConditionType::Temperature && m_spec_idx >= 0) {
    const amrex::Real Y_spec = rY_in(i, j, k, m_spec_idx) / rho;
    cond_Y = (Y_spec < m_Y_threshold);
  }

  switch (m_cond_type) {
  case ConditionType::Temperature: return cond_T;
  case ConditionType::Species:     return cond_Y;
  case ConditionType::CombinedAnd: return cond_T && cond_Y;
  case ConditionType::CombinedOr:  return cond_T || cond_Y;
  default:                         return cond_T;
  }
}

// ============================================================
// react — box-level (primary PeleLMeX path)
// ============================================================

int
ReactorAdaptive::react(
  const amrex::Box& box,
  amrex::Array4<amrex::Real> const& rY_in,
  amrex::Array4<amrex::Real> const& rYsrc_in,
  amrex::Array4<amrex::Real> const& T_in,
  amrex::Array4<amrex::Real> const& rEner_in,
  amrex::Array4<amrex::Real> const& rEner_src_in,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Array4<int> const& mask,
  amrex::Real& dt_react,
  amrex::Real& time
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t /*stream*/
#endif
)
{
  BL_PROFILE("Pele::ReactorAdaptive::react()");

#ifdef AMREX_USE_GPU
  amrex::Abort(
    "ReactorAdaptive (Phase 1) is CPU-only. "
    "GPU builds are not yet supported.");
  return -1;
#else

  const int neq = NUM_SPECIES + 1;
  int total_nfe = 0;
  int n_primary = 0, n_fallback = 0;

  // Iterate over all cells in the box sequentially on CPU.
  // Each cell is integrated independently with a 1-cell buffer (ncells=1).
  // With ncells=1, YCOrder and CYOrder are identical, so the same buffer
  // layout works for both ReactorQSS (YCOrder) and ReactorCvode (CYOrder).

  amrex::LoopOnCpu(box, [&](int i, int j, int k) noexcept {

    // Skip EB-covered cells
    if (mask(i, j, k) == -1) {
      FC_in(i, j, k, 0) = 0.0;
      return;
    }

    // ── Pack: one-cell flat buffers ───────────────────────────────────────
    // Layout (YCOrder, ncells=1):
    //   rY_cell[s]           = rhoY_s  for s in [0, NUM_SPECIES)
    //   rY_cell[NUM_SPECIES] = T       (recomputed from energy for consistency)
    //   rYsrc_cell[s]        = rhoYsrc_s
    //   rX_cell[0]           = rhoH (or rhoE)
    //   rXsrc_cell[0]        = energy source

    amrex::Real rY_cell[NUM_SPECIES + 1];
    amrex::Real rYsrc_cell[NUM_SPECIES];
    amrex::Real rX_cell[1];
    amrex::Real rXsrc_cell[1];

    // Use box_flatten to pack + recompute T from energy (same as other reactors)
    utils::box_flatten<utils::YCOrder>(
      /*icell=*/0, i, j, k, /*ncells=*/1,
      m_reactor_type, /*clean_init_massfrac=*/false,
      rY_in, rYsrc_in, T_in, rEner_in, rEner_src_in,
      rY_cell, rYsrc_cell, rX_cell, rXsrc_cell);

    // ── Condition: evaluate on recomputed T and pre-react species ─────────
    // T after box_flatten is consistent with the energy
    const amrex::Real T_cell = rY_cell[NUM_SPECIES];
    amrex::Real rho = 0.0;
    for (int s = 0; s < NUM_SPECIES; s++) {
      rho += rY_cell[s];
    }

    ReactorBase* reactor =
      use_primary(T_cell, rY_in, rho, i, j, k) ? m_primary.get()
                                                : m_fallback.get();

    if (reactor == m_primary.get()) {
      ++n_primary;
    } else {
      ++n_fallback;
    }

    // ── Integrate: single-cell call ───────────────────────────────────────
    // After react(), the buffers hold the updated state:
    //   rY_cell[s]           = updated rhoY_s
    //   rY_cell[NUM_SPECIES] = updated T (recomputed internally by sub-reactor)
    //   rX_cell[0]           = updated rhoH = old_rhoH + dt * rhoHsrc
    amrex::Real dt_cell   = dt_react;
    amrex::Real time_cell = time;
    const int nfe =
      reactor->react(rY_cell, rYsrc_cell, rX_cell, rXsrc_cell,
                     dt_cell, time_cell, /*ncells=*/1);
    total_nfe += nfe;

    // ── Unpack: direct copy back to Array4 ───────────────────────────────
    // The 1D react already applied the energy source contribution and
    // recomputed T. box_unflatten is intentionally NOT used here because it
    // would re-apply the energy source term a second time.
    for (int s = 0; s < NUM_SPECIES; s++) {
      rY_in(i, j, k, s) = rY_cell[s];
    }
    T_in(i, j, k, 0)    = rY_cell[NUM_SPECIES];
    rEner_in(i, j, k, 0) = rX_cell[0];
    FC_in(i, j, k, 0)    = static_cast<amrex::Real>(nfe);
  });

  if (m_adaptive_verbose > 0) {
    const int ncells = static_cast<int>(box.numPts());
    amrex::Print() << "[ReactorAdaptive] box=" << box
                   << "  primary=" << n_primary
                   << "  fallback=" << n_fallback
                   << "  skipped=" << (ncells - n_primary - n_fallback)
                   << "  total_nfe=" << total_nfe << "\n";
  }

  return total_nfe;
#endif // AMREX_USE_GPU
}

// ============================================================
// react — 1D pointer interface (delegates to primary)
// ============================================================

int
ReactorAdaptive::react(
  amrex::Real* rY_in,
  amrex::Real* rYsrc_in,
  amrex::Real* rX_in,
  amrex::Real* rX_src_in,
  amrex::Real& dt_react,
  amrex::Real& time,
  int ncells
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  BL_PROFILE("Pele::ReactorAdaptive::react()");
  // The 1D interface is not called by PeleLMeX in the normal path.
  // Delegate entirely to the primary sub-reactor.
  return m_primary->react(
    rY_in, rYsrc_in, rX_in, rX_src_in, dt_react, time, ncells
#ifdef AMREX_USE_GPU
    ,
    stream
#endif
  );
}

// ============================================================
// flatten / unflatten — delegate to primary
// ============================================================

void
ReactorAdaptive::flatten(
  const amrex::Box& box,
  int ncells,
  amrex::Array4<const amrex::Real> const& rhoY,
  amrex::Array4<const amrex::Real> const& frcExt,
  amrex::Array4<const amrex::Real> const& temperature,
  amrex::Array4<const amrex::Real> const& rhoE,
  amrex::Array4<const amrex::Real> const& frcEExt,
  amrex::Real* y_vect,
  amrex::Real* src_vect,
  amrex::Real* vect_energy,
  amrex::Real* src_vect_energy)
{
  m_primary->flatten(
    box, ncells, rhoY, frcExt, temperature, rhoE, frcEExt, y_vect, src_vect,
    vect_energy, src_vect_energy);
}

void
ReactorAdaptive::unflatten(
  const amrex::Box& box,
  int ncells,
  amrex::Array4<amrex::Real> const& rhoY,
  amrex::Array4<amrex::Real> const& temperature,
  amrex::Array4<amrex::Real> const& rhoE,
  amrex::Array4<amrex::Real> const& frcEExt,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Real* y_vect,
  amrex::Real* vect_energy,
  long int* FCunt,
  amrex::Real dt)
{
  m_primary->unflatten(
    box, ncells, rhoY, temperature, rhoE, frcEExt, FC_in, y_vect, vect_energy,
    FCunt, dt);
}

// ============================================================
// close / print_final_stats
// ============================================================

void
ReactorAdaptive::close()
{
  if (m_primary)  m_primary->close();
  if (m_fallback) m_fallback->close();
}

void
ReactorAdaptive::print_final_stats(void* sundials_mem)
{
  if (m_primary)  m_primary->print_final_stats(sundials_mem);
  if (m_fallback) m_fallback->print_final_stats(sundials_mem);
}

} // namespace pele::physics::reactions
