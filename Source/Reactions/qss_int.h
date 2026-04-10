/**
 * @file qss_int.h
 * @brief Quasi-Steady-State (QSS) integrator for stiff ODEs (PelePhysics reactor format)
 *
 * Based on CHEMEQ2 / Mott, Oran, van Leer (2000) JCP 164(2), 407-428.
 * Splits ODEs as dy/dt = q(y) - d(y) (production minus destruction).
 */

 #ifndef QSS_INT_H
 #define QSS_INT_H
 
 #include <AMReX_REAL.H>
 #include <AMReX_Vector.H>
 #include <cstddef>
 
 namespace pele::physics::reactions {
 
 /** Abstract ODE interface for QSS: provides production q and destruction d. */
 class QssOde
 {
 public:
   virtual ~QssOde() = default;
   /** Compute q and d at (t, y). corrector=true for corrector step. */
   virtual void odefun(
     amrex::Real t,
     const amrex::Vector<amrex::Real>& y,
     amrex::Vector<amrex::Real>& q,
     amrex::Vector<amrex::Real>& d,
     bool corrector = false) = 0;
 };
 
 /** QSS integrator: advances state y from tstart to tf using q-d split. */
 class QssIntegrator
 {
 public:
   QssIntegrator();
   void setOde(QssOde* ode);
   /** Allocate for neq equations (e.g. NUM_SPECIES + 1). */
   void initialize(size_t neq);
   /** Set initial state and start time; resets internal step state. */
   void setState(const amrex::Vector<amrex::Real>& yIn, amrex::Real tstart_);
   /** Integrate from tstart to tf. Returns 0 on success, non-zero on failure. */
   int integrateToTime(amrex::Real tf);
   /** Get solution after integrateToTime: copy internal y into yOut. */
   void getState(amrex::Vector<amrex::Real>& yOut) const;
 
   // Tolerances and limits (ParmParse in reactor init)
   amrex::Real epsmax{20.0};
   amrex::Real epsmin{1e-2};
   amrex::Real dtmin{1e-15};
   amrex::Real dtmax{1e-6};
   int itermax{2};
   amrex::Real tfd{1.000008};
   amrex::Real abstol{1e-8};
   bool stabilityCheck{true};
   amrex::Vector<amrex::Real> ymin;
   amrex::Vector<amrex::Real> enforce_ymin;
   // Optional per-variable upper bound.  When enforce_ymax[i] != 0, y[i] is
   // clamped to ymax[i] after every update.  Use for temperature to stay
   // inside the valid range of the thermodynamic polynomials and prevent the
   // predictor from overshooting into regions where EOS returns garbage.
   amrex::Vector<amrex::Real> ymax;
   amrex::Vector<amrex::Real> enforce_ymax;

 private:
   void getInitialStepSize(amrex::Real tf);
   int integrateOneStep(amrex::Real tf);
 
   QssOde* ode_{nullptr};
   size_t N{0};
   bool firstStep{true};
   amrex::Real dt{0};
   amrex::Real ts{0};
   amrex::Real tstart{0};
   amrex::Real tn{0};
   int gcount{0};
   int rcount{0};
 
   amrex::Vector<amrex::Real> y, q, d, rtaus, y1, ys, rtau, qs;
   amrex::Vector<amrex::Real> ym1, ym2, scratch;
 };
 
 } // namespace pele::physics::reactions
 
 #endif
