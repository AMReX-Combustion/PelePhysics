
#ifndef TABBREAKUP_H
#define TABBREAKUP_H

// This is the implementation of the TAB breakup model by O'Rourke and Amsden
// (1987)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
updateBreakup(
  const amrex::Real C_D,
  const amrex::Real rem_dt,
  const amrex::Real dt,
  int pid,
  const GasPhaseVals& gpv,
  const SprayComps& SPI,
  const SprayData& fdat,
  SprayParticleContainer::ParticleType& p,
  splash_breakup* N_SB,
  const SBPtrs& rf)
{
  // Model constants
  const amrex::Real C_k = 8.;
  // O'Rourke and Amsden list this value as 5 but Marek in "The double-mass
  // model of drop..." from 2013 says this value should actually be 10
  const amrex::Real C_d = 10.;
  const amrex::Real C_b = 0.5;
  const amrex::Real C_F = 1. / 3.;
  // K is used to determine the SMD of the resulting droplets
  const amrex::Real K = 10. / 3.;

  // Retreive particle data
  amrex::RealVect vel_part(AMREX_D_DECL(
    p.rdata(SPI.pstateVel), p.rdata(SPI.pstateVel + 1),
    p.rdata(SPI.pstateVel + 2)));
  amrex::Real T_part = p.rdata(SPI.pstateT);
  amrex::Real dia_part = p.rdata(SPI.pstateDia);
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_part; // Liquid mass fractions
  amrex::Real rho_part = 0.;
  amrex::Real mu_part = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_part[spf] = p.rdata(SPI.pstateY + spf);
    rho_part += Y_part[spf] / fdat.rhoL(T_part, spf);
    mu_part += Y_part[spf] * fdat.mu[spf];
  }
  rho_part = 1. / rho_part;
  amrex::Real sigma = fdat.sigma;
  amrex::Real rad_part = 0.5 * dia_part;
  amrex::RealVect diff_vel = gpv.vel_fluid - vel_part;
  amrex::Real We = gpv.rho_fluid * diff_vel.radSquared() * rad_part / sigma;
  amrex::Real We_crit = C_k * C_b / C_F;
  amrex::Real Wer = We / We_crit;
  amrex::Real denom = rho_part * rad_part * rad_part;
  amrex::Real td = 2. * denom / (C_d * mu_part);
  amrex::Real omega2 = C_k * sigma / (denom * rad_part) - 1. / (td * td);
  if (omega2 <= 0.) {
    p.rdata(SPI.pstateTABY) = 0.;
    p.rdata(SPI.pstateTABYdot) = 0.;
    return;
  }
  amrex::Real yn = p.rdata(SPI.pstateTABY);
  amrex::Real ydotn = p.rdata(SPI.pstateTABYdot);
  amrex::Real A2 = std::pow(yn - Wer, 2) + ydotn * ydotn / omega2;
  amrex::Real A = std::sqrt(A2);
  amrex::Real omega = std::sqrt(omega2);
  amrex::Real ynp =
    Wer + std::exp(-dt / td) *
            ((yn - Wer) * std::cos(omega * dt) +
             1. / omega * (ydotn + (yn - Wer) / td) * std::sin(omega * dt));
  amrex::Real ydotnp =
    (Wer - ynp) / td +
    omega * std::exp(-dt / td) *
      (1. / omega * (ydotn + (yn - Wer) / td) * std::cos(omega * dt) -
       (yn - Wer) * std::sin(omega * dt));
  if (Wer + A <= 1.) {
    p.rdata(SPI.pstateTABY) = ynp;
    p.rdata(SPI.pstateTABYdot) = ydotnp;
    return;
  }
  amrex::Real pv1 = (yn - Wer) / A;
  amrex::Real pv2 = -yn / (A * omega);
  amrex::Real phi1 = std::acos(pv1);
  amrex::Real phi2 = std::asin(pv2);
  amrex::Real mnv = std::acos((1. - Wer) / A);
  amrex::Real tb1 = (mnv - phi1) / omega;
  amrex::Real tb2 = (mnv - phi2) / omega;
  amrex::Real tbv1 = (tb1 > 0.) ? tb1 : 2. * dt;
  amrex::Real tbv2 = (tb2 > 0.) ? tb2 : 2. * dt;
  amrex::Real tbv = amrex::min(tbv1, tbv2);
  // Check if breakup occurs
  if (tbv < dt) {
    // Evaluate ydot at tb
    amrex::Real ytb =
      Wer + std::exp(-tbv / td) *
              ((yn - Wer) * std::cos(omega * tbv) +
               1. / omega * (ydotn + (yn - Wer) / td) * std::sin(omega * tbv));
    amrex::Real ydottb =
      (Wer - ytb) / td +
      omega * std::exp(-tbv / td) *
        (1. / omega * (ydotn + (yn - Wer) / td) * std::cos(omega * tbv) -
         (yn - Wer) * std::sin(omega * tbv));
    amrex::Real Utan, r32;
    // Original TAB model
    if (fdat.use_ETAB) {
      amrex::Real Wet = 80.;
      amrex::Real k1 = 2. / 9.;
      amrex::Real k2 = 2. / 9.;
      k1 *= (k2 / k1 * (std::sqrt(Wet) - 1.)) * std::pow(We / Wet, 4) + 1.;
      amrex::Real Kbr = k1 * omega;
      if (We > Wet) {
        Kbr = k2 * omega * std::sqrt(We);
      }
      r32 = rad_part * std::exp(-Kbr * tbv);
      amrex::Real rsmr = std::sqrt(std::pow(rad_part, 3) / (r32 * r32));
      amrex::Real AE2 = 3. * (1. - rad_part / rsmr + 5. * C_D * We / 72.) *
                        omega2 / (ydottb * ydottb);
      amrex::Real AE = std::sqrt(AE2);
      Utan = AE * C_b * rad_part * ydottb;
    } else {
      r32 = 3. * rad_part /
            (1. + 8. * K / 20. +
             denom * rad_part / sigma * ydottb * ydottb * (6. * K - 5.) / 120.);
      Utan = C_b * rad_part * ydottb;
    }
    N_SB[pid] = splash_breakup::breakup;
    amrex::Real velMag = vel_part.vectorLength();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      // Save location of droplet at tb, which is between tb < dt
      rf.loc[AMREX_SPACEDIM * pid + dir] =
        p.pos(dir) + (tbv - dt) * vel_part[dir];
      rf.vel[AMREX_SPACEDIM * pid + dir] = vel_part[dir];
      rf.norm[AMREX_SPACEDIM * pid + dir] = vel_part[dir] / velMag;
    }
    rf.d0[pid] = p.rdata(SPI.pstateDia);
    rf.dtpp[pid] = dt - tbv + rem_dt;
    rf.phi1[pid] = r32;
    rf.phi2[pid] = Utan;
    rf.T0[pid] = p.rdata(SPI.pstateT);
#if SPRAY_FUEL_NUM > 1
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      rf.Y0[SPRAY_FUEL_NUM * pid + spf] = p.rdata(SPI.pstateY + spf);
    }
#endif
    p.id() = -1;
  } else {
    p.rdata(SPI.pstateTABY) = ynp;
    p.rdata(SPI.pstateTABYdot) = ydotnp;
    return;
  }
}
#endif
