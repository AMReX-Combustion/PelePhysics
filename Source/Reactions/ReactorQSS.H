#ifndef REACTORQSS_H
#define REACTORQSS_H

#include "ReactorBase.H"
#include "qss_int.h"

namespace pele::physics::reactions {

/**
 * @brief QSS (Quasi-Steady-State) chemistry reactor for PeleLMeX.
 *
 * Implements the CHEMEQ2 algorithm (Mott, Oran & van Leer, JCP 164(2), 2000)
 * using the QssIntegrator in qss_int.h. The ODE is split as dy/dt = q(y) - d(y)
 * where q >= 0 (production) and d >= 0 (destruction). The production/destruction
 * split uses a sign-split of the net chemistry RHS: q = max(0, wdot + src),
 * d = max(0, -(wdot + src)).
 *
 * State vector per cell (YCOrder): [rhoY_0, ..., rhoY_{N-1}, T]
 *
 * CPU-ONLY: QssOde uses virtual dispatch which is not GPU-callable. On GPU-enabled
 * builds the react() methods abort immediately. Use ReactorCvode for GPU runs.
 *
 * ParmParse parameters:
 *   ode.verbose            (int,  default 0)     - verbosity level
 *   ode.clean_init_massfrac (int, default 0)     - clip negative mass fractions
 *   qss.epsmax             (Real, default 20.0)  - max step error multiplier
 *   qss.epsmin             (Real, default 1e-2)  - min step error multiplier
 *   qss.dtmin              (Real, default 1e-15) - min internal substep [s]
 *   qss.dtmax              (Real, default 1e-6)  - max internal substep [s]
 *   qss.itermax            (int,  default 2)     - max corrector iterations
 *   qss.tfd                (Real, default 1.000008) - finish-time tolerance factor
 *   qss.abstol             (Real, default 1e-8)  - absolute tolerance for error
 *   qss.stability_check    (int,  default 1)     - enable stability limiter
 */
class ReactorQSS : public ReactorBase::Register<ReactorQSS>
{
public:
  static std::string identifier() { return "ReactorQSS"; }

  // YCOrder: [rhoY_0..rhoY_{N-1}, T]_{cell0}, [...]_{cell1}, ...
  using Ordering = utils::YCOrder;

  int init(int reactor_type, int ncells) override;

  void close() override {}

  void print_final_stats(void* /*mem*/) override {}

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

  void flatten(
    const amrex::Box& box,
    const int ncells,
    amrex::Array4<const amrex::Real> const& rhoY,
    amrex::Array4<const amrex::Real> const& frcExt,
    amrex::Array4<const amrex::Real> const& temperature,
    amrex::Array4<const amrex::Real> const& rhoE,
    amrex::Array4<const amrex::Real> const& frcEExt,
    amrex::Real* y_vect,
    amrex::Real* src_vect,
    amrex::Real* vect_energy,
    amrex::Real* src_vect_energy) override
  {
    flatten_ops.flatten(
      box, ncells, m_reactor_type, m_clean_init_massfrac, rhoY, frcExt,
      temperature, rhoE, frcEExt, y_vect, src_vect, vect_energy,
      src_vect_energy);
  }

  void unflatten(
    const amrex::Box& box,
    const int ncells,
    amrex::Array4<amrex::Real> const& rhoY,
    amrex::Array4<amrex::Real> const& temperature,
    amrex::Array4<amrex::Real> const& rhoE,
    amrex::Array4<amrex::Real> const& frcEExt,
    amrex::Array4<amrex::Real> const& FC_in,
    amrex::Real* y_vect,
    amrex::Real* vect_energy,
    long int* FCunt,
    amrex::Real dt) override
  {
    flatten_ops.unflatten(
      box, ncells, m_reactor_type, m_clean_init_massfrac, rhoY, temperature,
      rhoE, frcEExt, FC_in, y_vect, vect_energy, FCunt, dt);
  }

  void set_eos_parm(
    const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
      h_eosparm,
    const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
      d_eosparm) override
  {
    m_h_eosparm = h_eosparm;
    m_d_eosparm = d_eosparm;
  }

private:
  /**
   * @brief Per-cell chemistry ODE for QSS integration.
   *
   * State vector y has NUM_SPECIES + 1 entries:
   *   y[0..NUM_SPECIES-1] = rhoY_n   (mass densities, g/cm^3)
   *   y[NUM_SPECIES]       = T        (temperature, K)
   *
   * odefun provides the production (q) and destruction (d) split
   * such that dy/dt = q - d, where q >= 0 and d >= 0.
   * The split is by sign of the net chemistry RHS.
   *
   * Not thread-safe individually — one CellOde is created locally per
   * integrate_cells call so concurrent react() calls are safe.
   */
  class CellOde : public QssOde
  {
  public:
    void setup(
      int reactor_type,
      const amrex::Real rYsrc[NUM_SPECIES],
      amrex::Real rhoesrc_ext,
      amrex::Real tstart,
      const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
        eosparm) noexcept;

    void odefun(
      amrex::Real t,
      const amrex::Vector<amrex::Real>& y,
      amrex::Vector<amrex::Real>& q,
      amrex::Vector<amrex::Real>& d,
      bool corrector = false) override;

    int nfeval() const noexcept { return m_nfeval; }
    void resetCount() noexcept { m_nfeval = 0; }

  private:
    int m_reactor_type{2}; // ReactorTypes::h_reactor_type (enthalpy)
    amrex::Real m_rYsrc[NUM_SPECIES] = {0.0};
    amrex::Real m_rhoesrc_ext{0.0};
    amrex::Real m_tstart{0.0};
    const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
      m_eosparm{nullptr};
    int m_nfeval{0};
  };

  /**
   * @brief Integrate all cells in y_vect from time to time+dt_react on the CPU.
   *
   * Thread-safe: creates local QssIntegrator and CellOde instances.
   *
   * @param y_vect       [in/out] flat state: [rhoY_0..rhoY_{N-1}, T] per cell
   * @param src_vect     [in]     flat species external sources per cell
   * @param src_vect_energy [in]  per-cell energy external source (rhoesrc_ext)
   * @param dt_react     [in]     integration time interval
   * @param time         [in]     current simulation time (absolute)
   * @param ncells       [in]     number of cells to integrate
   * @param FCunt        [out]    per-cell function evaluation count
   */
  void integrate_cells(
    amrex::Real* y_vect,
    const amrex::Real* src_vect,
    const amrex::Real* src_vect_energy,
    amrex::Real dt_react,
    amrex::Real time,
    int ncells,
    long int* FCunt);

  utils::FlattenOps<Ordering> flatten_ops;
  int m_reactor_type{0};
  bool m_clean_init_massfrac{false};

  // QSS integrator parameters (set in init(), configurable via ParmParse)
  amrex::Real m_epsmax{20.0};
  amrex::Real m_epsmin{1.0e-2};
  amrex::Real m_dtmin{1.0e-15};
  amrex::Real m_dtmax{1.0e-6};
  int m_itermax{2};
  amrex::Real m_tfd{1.000008};
  amrex::Real m_abstol{1.0e-8};
  amrex::Real m_Tmax{4000.0}; // upper T bound for integrator state [K]
  bool m_stability_check{true};
};

} // namespace pele::physics::reactions
#endif
