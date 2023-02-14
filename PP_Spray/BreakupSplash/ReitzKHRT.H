
#ifndef REITZKHRT_H
#define REITZKHRT_H

// This is the implementation of the wave model detailed by Reitz (1987)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
updateBreakupKHRT(
  const int pid,
  SprayParticleContainer::ParticleType& p,
  const amrex::Real& Reyn_d,
  const amrex::Real& dt,
  const amrex::Real& avg_inject_d3,
  const amrex::Real B0,
  const amrex::Real B1,
  const amrex::Real C3,
  const GasPhaseVals& gpv,
  const SprayData& fdat,
  splash_breakup* N_SB,
  const SBPtrs& rf)
{
  // Percent of mass lost before new parcel is formed
  const amrex::Real massperc = 0.03;
  const amrex::Real We_crit = 6.;
  // RT model constants
  const amrex::Real A1 = 0.188;
  // Retreive particle data
  amrex::RealVect vel_part(AMREX_D_DECL(
    p.rdata(SprayComps::pstateVel), p.rdata(SprayComps::pstateVel + 1),
    p.rdata(SprayComps::pstateVel + 2)));
  amrex::Real T_part = p.rdata(SprayComps::pstateT);
  amrex::Real dia_part = p.rdata(SprayComps::pstateDia);
  amrex::Real num_dens = p.rdata(SprayComps::pstateNumDens);
  amrex::Real rho_part = 0.;
  amrex::Real mu_part = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    amrex::Real Y_part = p.rdata(SprayComps::pstateY + spf);
    rho_part += Y_part / fdat.rhoL(T_part, spf);
    mu_part += Y_part * fdat.muL(T_part, spf);
  }
  rho_part = 1. / rho_part;
  const amrex::Real rjet3 = avg_inject_d3 / 8.;
  const amrex::Real rjet = std::cbrt(r03);
  const amrex::Real Njet = fdat.inj_num_ppp;
  amrex::Real d0 = p.rdata(SprayComps::pstateBphi1);
  amrex::Real r0 = 0.5 * d0;
  amrex::Real N0 = num_dens * std::pow(dia_part / d0, 3);
  amrex::Real rt_time = p.rdata(SprayComps::pstateBphi2);
  amrex::Real sigma = fdat.sigma;
  amrex::Real rad_part = 0.5 * dia_part;
  amrex::RealVect diff_vel = gpv.vel_fluid - vel_part;
  amrex::Real vel_mag = diff_vel.vectorLength();
  amrex::Real vel_sq = vel_mag * vel_mag;
  amrex::Real Reyn = Reyn_d;
  amrex::Real C_D = 0.;
  if (Reyn > 1000.) {
    C_D = 0.424;
  } else if (Reyn > 1.) {
    C_D = 24. / Reyn * (1. + std::cbrt(Reyn * Reyn) / 6.);
  } else if (Reyn > 0.) {
    C_D = 24. / Reyn;
  }
  amrex::Real Re_L = rho_part * vel_mag * rad_part / mu_part;
  amrex::Real We_L = rho_part * vel_sq * rad_part / sigma;
  amrex::Real We_g = gpv.rho_fluid * vel_sq * rad_part / sigma;
  amrex::Real Z = mu_part / std::sqrt(rho_part * rad_part * sigma);
  amrex::Real T = Z * std::sqrt(We_g);
  bool breakupRT = false;
  // Check of RT breakup occurs
  if (rt_time >= 0.) {
    // Add gravity here when applicable
    amrex::Real a =
      3. / 8. * C_D * gpv.rho_fluid * vel_sq / (rho_part * rad_part);
    amrex::Real momdiff = a * (rho_part - gpv.rho_fluid);
    amrex::Real omega_RT = std::sqrt(
      2. / (3. * std::sqrt(3. * sigma)) * momdiff * std::sqrt(momdiff) /
      (rho_part + gpv.rho_fluid));
    amrex::Real KRT = std::sqrt(momdiff / (3. * sigma));
    amrex::Real lam_RT = C3 * M_PI / KRT;
    if (rad_part < lam_RT) {
      rt_time = 0.;
    } else {
      rt_time += dt;
    }
    amrex::Real tau_RT = 1. / omega_RT;
    // RT breakup occurs
    if (rt_time > tau_RT && lam_RT < rad_part) {
      rt_time = 0.;
      num_dens *= std::pow(rad_part / lam_RT, 3);
      rad_part = lam_RT;
      breakupRT = true;
      if (num_dens > 1.) {
        N_SB[pid] = splash_breakup::breakup_RT;
        amrex::Real Utan = A1 * omega_RT * lam_RT;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          rf.loc[AMREX_SPACEDIM * pid + dir] = p.pos(dir);
          rf.vel[AMREX_SPACEDIM * pid + dir] = vel_part[dir];
          rf.norm[AMREX_SPACEDIM * pid + dir] = vel_part[dir] / vel_mag;
        }
        rf.numDens[pid] = num_dens;
        rf.ref_dia[pid] = 2. * rad_part;
        rf.dtpp[pid] = 0.;
        rf.phi1[pid] = Utan;
        rf.phi2[pid] = 0.; // unused
        rf.T0[pid] = T_part;
#if SPRAY_FUEL_NUM > 1
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          rf.Y0[SPRAY_FUEL_NUM * pid + spf] =
            p.rdata(SprayComps::pstateY + spf);
        }
#endif
        p.id() = -1;
        return;
      }
    }
  }
  if (!breakupRT) {
    // Check if KH breakup occurs
    amrex::Real lam_KH = rad_part * 9.02 * (1. + 0.45 * std::sqrt(Z)) *
                         (1. + 0.4 * std::pow(T, 0.7)) /
                         std::pow(1. + 0.865 * std::pow(We_g, 1.67), 0.6);
    amrex::Real omega_KH =
      (0.34 + 0.385 * We_g * std::sqrt(We_g)) /
      ((1. + Z) * (1. + 1.4 * std::pow(T, 0.6))) *
      std::sqrt(sigma / (rho_part * std::pow(rad_part, 3)));
    amrex::Real tau_KH = 3.726 * B1 * rad_part / (lam_KH * omega_KH);
    amrex::Real rs = B0 * lam_KH;
    // KH breakup occurs
    if (rs < rad_part && We_g > We_crit) {
      amrex::Real rp = (dt / tau_KH * rs + rad_part) / (1. + dt / tau_KH);
      amrex::Real np = num_dens * std::pow(rad_part / rp, 3);
      // KH breakup occurs
      if (rs < rad_part && We_g > We_crit) {
        rt_time = 0.;
        amrex::Real rp2 = rp * rp;
        amrex::Real rs2 = rs * rs;
        amrex::Real shed_mass = rp2 * (np - N0);
        amrex::Real ns = shed_mass / (rs2 * rs);
        amrex::Real removed_perc =
          N0 * (std::pow(r0, 3) - rp2 * rp) / (Njet * rjet3);
        amrex::Real Utan = A1 * omega_KH * lam_KH;
        // If removed mass is greater than 3%, break particle up
        if (ns > num_dens && removed_perc > massperc) {
          N_SB[pid] = splash_breakup::breakup_KH;
          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            rf.loc[AMREX_SPACEDIM * pid + dir] = p.pos(dir);
            rf.vel[AMREX_SPACEDIM * pid + dir] = vel_part[dir];
            rf.norm[AMREX_SPACEDIM * pid + dir] = vel_part[dir] / vel_mag;
          }
          rf.numDens[pid] = ns;
          rf.ref_dia[pid] = 2. * rs;
          rf.dtpp[pid] = 0.;
          rf.T0[pid] = T_part;
          rf.phi1[pid] = Utan;
          rf.phi2[pid] = 0.; // Unused
#if SPRAY_FUEL_NUM > 1
          for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
            rf.Y0[SPRAY_FUEL_NUM * pid + spf] =
              p.rdata(SprayComps::pstateY + spf);
          }
#endif
          rad_part = rp;
          num_dens = (np * rp2 * rp - ns * rs2 * rs) / (rp2 * rp);
          d0 = 2. * rp; // Reset original diameter
        } else {
          rad_part = rp;
          num_dens = np;
        }
      }
    }
  }
  p.rdata(SprayComps::pstateDia) = 2. * rad_part;
  p.rdata(SprayComps::pstateNumDens) = num_dens;
  p.rdata(SprayComps::pstateBphi1) = d0;
  p.rdata(SprayComps::pstateBphi2) = rt_time;
  return;
}

#endif
