/**
 * @file qss_int.cpp
 * @brief Implementation of the Quasi-Steady-State (QSS) integrator for stiff ODEs
 *
 * Based on CHEMEQ2 / Mott, Oran, van Leer (2000) JCP 164(2), 407-428.
 * dy/dt = q(y) - d(y).
 */

 #include "qss_int.h"
 #include <AMReX_Print.H>
 #include <algorithm>
 #include <cmath>
 #include <iostream>
 
 namespace pele::physics::reactions {
 
 namespace {
 inline int
 sign(amrex::Real x)
 {
   return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
 }
 inline bool
 notnan(const amrex::Vector<amrex::Real>& v)
 {
   for (amrex::Real val : v) {
     if (std::isnan(val)) return false;
   }
   return true;
 }
 } // namespace
 
 QssIntegrator::QssIntegrator() = default;
 
 void
 QssIntegrator::setOde(QssOde* ode)
 {
   ode_ = ode;
 }
 
 void
 QssIntegrator::initialize(size_t neq)
 {
   N = neq;
   y.resize(N);
   q.resize(N);
   d.resize(N);
   rtaus.resize(N);
   y1.resize(N);
   ys.resize(N);
   rtau.resize(N);
   qs.resize(N);
   ym1.resize(N);
   ym2.resize(N);
   scratch.resize(N);
   ymin.assign(N, 1e-20);
   enforce_ymin.assign(N, 1.0);
   ymax.assign(N, 1e30);      // no upper bound by default
   enforce_ymax.assign(N, 0.0);
 }
 
 void
 QssIntegrator::setState(const amrex::Vector<amrex::Real>& yIn, amrex::Real tstart_)
 {
   AMREX_ASSERT(yIn.size() == N);
   AMREX_ASSERT(notnan(yIn));
 
   for (size_t i = 0; i < N; i++) {
     y[i] = enforce_ymin[i] ? std::max(yIn[i], ymin[i]) : yIn[i];
   }
   gcount = 0;
   rcount = 0;
   tstart = tstart_;
   tn = 0.0;
   firstStep = true;
 }
 
 void
 QssIntegrator::getState(amrex::Vector<amrex::Real>& yOut) const
 {
   yOut.resize(N);
   for (size_t i = 0; i < N; i++) {
     yOut[i] = y[i];
   }
 }
 
 void
 QssIntegrator::getInitialStepSize(amrex::Real tf)
 {
   firstStep = false;
   amrex::Real scratch_value = 1.0e-25;
 
   for (size_t i = 0; i < N; i++) {
     if (std::abs(y[i]) > abstol) {
       const amrex::Real absq = std::abs(q[i]);
       const amrex::Real scr2 =
         std::abs(1.0 / y[i]) * static_cast<amrex::Real>(sign(0.1 * epsmin * absq - d[i]));
       const amrex::Real scr1 = scr2 * d[i];
       scratch_value =
         std::max(std::max(scr1, -std::abs(absq - d[i]) * scr2), scratch_value);
     }
   }
 
   const amrex::Real sqreps = 0.5;
   dt = std::min(sqreps / scratch_value, tf);
   dt = std::min(dt, dtmax);
 }
 
 int
 QssIntegrator::integrateToTime(amrex::Real tf)
 {
   while (tfd * tn < tf) {
     int ret = integrateOneStep(tf);
     if (ret != 0) return ret;
   }
   return 0;
 }
 
 int
 QssIntegrator::integrateOneStep(amrex::Real tf)
 {
   AMREX_ASSERT(notnan(y));
   ode_->odefun(tn + tstart, y, q, d);
   AMREX_ASSERT(notnan(q));
   AMREX_ASSERT(notnan(d));
   gcount += 1;
 
   if (firstStep) {
     getInitialStepSize(tf);
   }
 
   ts = tn;
   for (size_t i = 0; i < N; i++) {
     rtau[i] = enforce_ymin[i] ? dt * d[i] / y[i] : 0.0;
   }
   qs = q;
   ys = y;
   rtaus = rtau;
 
   while (true) {
     for (size_t i = 0; i < N; i++) {
       amrex::Real denom =
         1.0 +
         rtau[i] * (180 + rtau[i] * (60 + rtau[i] * (11 + rtau[i]))) /
           (360 + rtau[i] * (60 + rtau[i] * (12 + rtau[i])));
       scratch[i] = (q[i] - d[i]) / denom;
     }
 
     amrex::Real eps = 1e-10;
     for (int iter = 0; iter < itermax; iter++) {
       if (stabilityCheck) {
         ym2 = ym1;
         ym1 = y;
       }
 
       for (size_t i = 0; i < N; i++) {
         amrex::Real new_val = ys[i] + dt * scratch[i];
         new_val = std::max(new_val, ymin[i]);
         if (enforce_ymax[i]) new_val = std::min(new_val, ymax[i]);
         y[i] = new_val;
       }

       if (iter == 0) {
         tn = ts + dt;
         y1 = y;
       }
 
       ode_->odefun(tn + tstart, y, q, d, true);
       AMREX_ASSERT(notnan(q));
       AMREX_ASSERT(notnan(d));
       gcount += 1;
 
       amrex::Vector<amrex::Real> rtaub(N);
       for (size_t i = 0; i < N; i++) {
         rtaub[i] =
           enforce_ymin[i] ? 0.5 * (rtaus[i] + dt * d[i] / y[i]) : 0.0;
       }
 
       amrex::Vector<amrex::Real> alpha(N);
       for (size_t i = 0; i < N; i++) {
         alpha[i] =
           (180. + rtaub[i] * (60. + rtaub[i] * (11. + rtaub[i]))) /
           (360. + rtaub[i] * (60. + rtaub[i] * (12. + rtaub[i])));
       }
 
       for (size_t i = 0; i < N; i++) {
         scratch[i] =
           (qs[i] * (1.0 - alpha[i]) + q[i] * alpha[i] - ys[i] * rtaub[i] / dt) /
           (1.0 + alpha[i] * rtaub[i]);
       }
     }
 
     eps = 0.0;
     for (size_t i = 0; i < N; i++) {
       amrex::Real new_y = ys[i] + dt * scratch[i];
       if (enforce_ymin[i]) new_y = std::max(new_y, ymin[i]);
       if (enforce_ymax[i]) new_y = std::min(new_y, ymax[i]);
       amrex::Real error = std::abs(new_y - y1[i]);
       new_y = std::max(new_y, ymin[i]);
       if (enforce_ymax[i]) new_y = std::min(new_y, ymax[i]);
       y[i] = new_y;
 
       if (std::abs(y[i]) > abstol && 0.25 * (ys[i] + y[i]) > ymin[i]) {
         error /= y[i];
         eps = std::max(
           .5 * (error + std::min(std::abs(q[i] - d[i]) / (q[i] + d[i] + 1e-30), error)),
           eps);
       }
     }
     AMREX_ASSERT(notnan(y));
 
     if (stabilityCheck) {
       ym2 = ym1;
       ym1 = y;
     }
 
     eps /= epsmin;
 
     if (dt <= dtmin + 1e-16 * tn) {
       amrex::Print() << "QssIntegrator failed: timestep too small: dt = " << dt
                      << ", tn = " << tn << ", dtmin = " << dtmin << "\n";
       return -1;
     }
 
     amrex::Real stab = 0;
     if (stabilityCheck && itermax >= 3) {
       stab = 0.01;
       for (size_t i = 0; i < N; i++) {
         if (std::abs(y[i]) > abstol) {
           stab = std::max(
             stab,
             std::abs(y[i] - ym1[i]) / (std::abs(ym1[i] - ym2[i]) + 1e-20 * y[i]));
         }
       }
     }
 
     if (eps <= epsmax && stab <= 1.0) {
       if (tf <= tn * tfd) return 0;
     } else {
       tn = ts;
     }
 
     amrex::Real rteps = 0.5 * (eps + 1.0);
     rteps = 0.5 * (rteps + eps / rteps);
     rteps = 0.5 * (rteps + eps / rteps);
 
     amrex::Real dto = dt;
     dt = std::min(dt * (1.0 / rteps + 0.005), tfd * (tf - tn));
     dt = std::min(dt, dtmax);
     if (stabilityCheck) {
       dt = std::min(dt, dto / (stab + 0.001));
     }
 
     if (eps > epsmax || stab > 1.0) {
       rcount += 1;
       for (size_t i = 0; i < N; i++) {
         rtaus[i] *= dt / dto;
       }
     } else {
       return 0;
     }
   }
 }
 
 } // namespace pele::physics::reactions
 