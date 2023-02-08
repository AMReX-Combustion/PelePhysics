
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
  const GasPhaseVals& gpv,
  const SprayData& fdat,
  splash_breakup* N_SB,
  const SBPtrs& rf)
{
  // Percent of mass lost before new parcel is formed
  const amrex::Real massperc = 0.03;
  const amrex::Real We_crit = 6.;
  // Wave model constants
  const amrex::Real B0 = 0.61; // Make this an input
  const amrex::Real B1 = 7.;   // Make this an input
  // RT model constants
  const amrex::Real A1 = 0.188;
  // Retreive particle data
  amrex::RealVect vel_part(AMREX_D_DECL(
    p.rdata(SprayComps::pstateVel), p.rdata(SprayComps::pstateVel + 1),
    p.rdata(SprayComps::pstateVel + 2)));
  amrex::Real T_part = p.rdata(SprayComps::pstateT);
  amrex::Real dia_part = p.rdata(SprayComps::pstateDia);
  amrex::Real num_dens = p.rdata(SprayComps::pstateNumDens);
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_part; // Liquid mass fractions
  amrex::Real rho_part = 0.;
  amrex::Real mu_part = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_part[spf] = p.rdata(SprayComps::pstateY + spf);
    rho_part += Y_part[spf] / fdat.rhoL(T_part, spf);
    mu_part += Y_part[spf] * fdat.muL(T_part, spf);
  }
  rho_part = 1. / rho_part;
  const amrex::Real r03 = avg_inject_d3 / 8.;
  const amrex::Real r0 = std::cbrt(r03);
  const amrex::Real L_break = B1 * r0 * std::sqrt(rho_part / gpv.rho_fluid);
  amrex::Real dis = p.rdata(SprayComps::pstateBphi1);
  amrex::Real rt_time = p.rdata(SprayComps::pstateBphi2);
  amrex::Real CTau = 9.;
  if (dis > L_break) {
    CTau = 1.;
  }
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
    amrex::Real lam_RT = M_PI / KRT;
    if (rt_time > 0. || lam_RT < rad_part) {
      rt_time += dt;
    }
    amrex::Real tau_RT = CTau / omega_RT;
    // RT breakup occurs
    if (rt_time > tau_RT && lam_RT < rad_part) {
      rt_time = 0.;
      num_dens *= std::pow(rad_part / lam_RT, 3);
      rad_part = lam_RT;
      breakupRT = true;
      if (num_dens > fdat.max_num_ppp) {
        amrex::Real Utan = A1 * omega_RT * lam_RT;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          rf.loc[AMREX_SPACEDIM * pid + dir] = p.pos(dir);
          rf.vel[AMREX_SPACEDIM * pid + dir] = vel_part[dir];
          rf.norm[AMREX_SPACEDIM * pid + dir] = vel_part[dir] / vel_mag;
        }
        rf.numDens[pid] = p.rdata(SprayComps::pstateNumDens);
        rf.d0[pid] = p.rdata(SprayComps::pstateDia);
        rf.dtpp[pid] = 0.;
        rf.phi1[pid] = Utan;
        rf.phi2[pid] = dis;
        rf.T0[pid] = p.rdata(SprayComps::pstateT);
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
      amrex::Real rad_mod = (dt / tau_KH * rs + rad_part) / (1. + dt / tau_KH);
      amrex::Real num_mod = num_dens * std::pow(rad_part / rad_mod, 3);
      // KH breakup occurs
      if (rs < rad_part && We_g > We_crit) {
        rt_time = 0.;
        amrex::Real rs2 = rs * rs;
        amrex::Real rp = rad_part;
        amrex::Real rp2 = rp * rp;
        amrex::Real np = num_dens;
        // Determine if new parcels should be formed to represent stripped
        // droplets. This is the method outlined by Patterson and Reitz
        // (1998). Solve a cubic equation by transforming it into a depressed
        // cubic
        amrex::Real nb = fdat.inj_num_ppp;
        amrex::Real pv = -rs2 / 3.;
        amrex::Real qv = -2. * rs2 * rs / 27. + np * rp2 * (rs - rp) / nb;
        amrex::Real R1 = std::sqrt(qv * qv / 4. + std::pow(pv, 3) / 27.);
        amrex::Real yt = std::cbrt(-qv / 2. + R1) + std::cbrt(-qv / 2. - R1);
        // Radius of larger droplet after smaller droplets are stripped
        amrex::Real rb = yt + rs / 3.;
        // Number density of smaller stripped droplets
        amrex::Real ns = (np * rp2 * rp - nb * std::pow(rb, 3)) / (rs2 * rs);
        amrex::Real removed_perc = (r03 - std::pow(rad_mod, 3)) / r03;
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
          rf.d0[pid] = 2. * rs;
          rf.dtpp[pid] = 0.;
          rf.T0[pid] = T_part;
          rf.phi1[pid] = Utan;
          rf.phi2[pid] = dis;
#if SPRAY_FUEL_NUM > 1
          for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
            rf.Y0[SPRAY_FUEL_NUM * pid + spf] =
              p.rdata(SprayComps::pstateY + spf);
          }
#endif
          rad_part = rb;
          num_dens = nb;
        } else {
          rad_part = rad_mod;
          num_dens = num_mod;
        }
      }
    }
  }
  p.rdata(SprayComps::pstateDia) = 2. * rad_part;
  p.rdata(SprayComps::pstateNumDens) = num_dens;
  p.rdata(SprayComps::pstateBphi2) = rt_time;
  return;
}

#endif
