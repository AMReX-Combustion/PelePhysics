#ifndef REACTORADAPTIVE_H
#define REACTORADAPTIVE_H

#include "ReactorBase.H"
#include <memory>
#include <string>
#include <vector>

namespace pele::physics::reactions {

/**
 * @brief Adaptive chemistry reactor that routes each cell to one of two
 *        sub-reactors based on a per-cell thermophysical condition.
 *
 * Motivation: explicit integrators (QSS, ARKODE) succeed in low-stiffness
 * (cold / non-reacting) cells and are fast; an implicit solver (CVODE) is
 * needed in stiff (flame / ignition) cells but is expensive.  This reactor
 * evaluates a user-specified condition on (T, Y_species) per cell each
 * timestep and dispatches accordingly.
 *
 * Condition semantics:  TRUE  → primary reactor (default: ReactorQSS)
 *                       FALSE → fallback reactor (default: ReactorCvode)
 *
 * Condition types (adaptive.condition_type):
 *   "temperature" : T_cell < T_threshold
 *   "species"     : Y_species_cell < Y_threshold
 *   "combined"    : temperature AND/OR species  (adaptive.condition_logic)
 *
 * Implementation notes (Phase 1, CPU-only):
 *   - Uses amrex::LoopOnCpu over the box; each cell is a separate 1-cell
 *     react() call to the selected sub-reactor.
 *   - ncells=1 ensures YCOrder == CYOrder, so the same flat buffer works for
 *     both ReactorQSS (YCOrder) and ReactorCvode (CYOrder).
 *   - clean_init_massfrac is NOT supported; keep ode.clean_init_massfrac = 0.
 *
 * ParmParse parameters (prefix "adaptive"):
 *   adaptive.primary_integrator   (string, default "ReactorQSS")
 *   adaptive.fallback_integrator  (string, default "ReactorCvode")
 *   adaptive.condition_type       (string, default "temperature")
 *   adaptive.condition_logic      (string, default "AND")   -- combined only
 *   adaptive.T_threshold          (Real,   default 1200.0)  -- K
 *   adaptive.species_name         (string, default "OH")
 *   adaptive.Y_threshold          (Real,   default 1.0e-6)  -- mass fraction
 *   adaptive.verbose              (int,    default 0)
 */
class ReactorAdaptive : public ReactorBase::Register<ReactorAdaptive>
{
public:
  static std::string identifier() { return "ReactorAdaptive"; }

  int  init(int reactor_type, int ncells) override;
  void close() override;
  void print_final_stats(void* sundials_mem) override;

  // Primary PeleLMeX path: box-level react, CPU-only
  int react(
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
    amrex::gpuStream_t stream
#endif
  ) override;

  // 1D pointer interface: delegates entirely to primary sub-reactor
  int react(
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
  ) override;

  // Flatten / unflatten: delegate to primary (not called by PeleLMeX for us)
  void flatten(
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
    amrex::Real* src_vect_energy) override;

  void unflatten(
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
    amrex::Real dt) override;

  // Forward EOS params to both sub-reactors; also resolves species index
  void set_eos_parm(
    const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
      h_eosparm,
    const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
      d_eosparm) override;

  // Forward typical values to both sub-reactors
  void set_typ_vals_ode(const std::vector<amrex::Real>& ExtTypVals) override;

private:
  enum class ConditionType { Temperature, Species, CombinedAnd, CombinedOr };

  // Returns true if this cell should use the primary (fast) reactor
  bool use_primary(
    amrex::Real T,
    amrex::Array4<amrex::Real> const& rY_in,
    amrex::Real rho,
    int i,
    int j,
    int k) const noexcept;

  std::unique_ptr<ReactorBase> m_primary;
  std::unique_ptr<ReactorBase> m_fallback;

  std::string m_primary_name  = "ReactorQSS";
  std::string m_fallback_name = "ReactorCvode";
  std::string m_spec_name     = "OH";

  ConditionType m_cond_type   = ConditionType::Temperature;
  amrex::Real   m_T_threshold = 1200.0;   // K
  int           m_spec_idx    = -1;        // resolved in set_eos_parm
  amrex::Real   m_Y_threshold = 1.0e-6;   // mass fraction

  int m_reactor_type     = 2; // h_reactor_type (matches PeleLMeX default)
  int m_adaptive_verbose = 0;
};

} // namespace pele::physics::reactions
#endif
