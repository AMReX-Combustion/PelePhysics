
#ifndef TABBREAKUP_H
#define TABBREAKUP_H

// This is the implementation of the TAB breakup model by O'Rourke and Amsden
// (1987) and the ETAB model by Tanner (1997)
// This function subcycles based on a dt_breakup from large We
// Continues looping until dt is reached in case multiple breakups occur within
// a timestep Returns the total tangential velocity from all breakups occuring
// over dt

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
updateBreakup(
  const amrex::Real& C_D,
  const amrex::Real& cur_time,
  const amrex::Real& dt,
  const GasPhaseVals& gpv,
  const SprayData& fdat,
  SprayParticleContainer::ParticleType& p,
  amrex::Real& breakup_time)
{
  // Model constants
  const amrex::Real C_k = 8.;
  // The original paper misprints this value as 5
  const amrex::Real C_d = 10.;
  const amrex::Real C_b = 0.5;
  const amrex::Real C_F = 1. / 3.;

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
  amrex::Real Utan_total = 0.;
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
    p.rdata(SprayComps::pstateTABY) = 0.;
    p.rdata(SprayComps::pstateTABYdot) = 0.;
    return 0.;
  }
  amrex::Real yn = p.rdata(SprayComps::pstateTABY);
  amrex::Real ydotn = p.rdata(SprayComps::pstateTABYdot);
  // Estimate breakup time assuming large We
  amrex::Real tbconst =
    2. * std::sqrt(3. * rho_part / gpv.rho_fluid) / diff_vel.vectorLength();
  amrex::Real tbest = rad_part * tbconst;
  amrex::Real subdt = dt / (0.1 * tbest);
  amrex::Real curt = 0.;
  while (curt < dt) {
    amrex::Real A2 = std::pow(yn - Wer, 2) + ydotn * ydotn / omega2;
    amrex::Real A = std::sqrt(A2);
    amrex::Real omega = std::sqrt(omega2);
    amrex::Real ynp =
      Wer + std::exp(-subdt / td) * ((yn - Wer) * std::cos(omega * subdt) +
                                     1. / omega * (ydotn + (yn - Wer) / td) *
                                       std::sin(omega * subdt));
    amrex::Real ydotnp =
      (Wer - ynp) / td +
      omega * std::exp(-subdt / td) *
        (1. / omega * (ydotn + (yn - Wer) / td) * std::cos(omega * subdt) -
         (yn - Wer) * std::sin(omega * subdt));
    if (Wer + A <= 1.) {
      yn = ynp;
      ydotn = ydotnp;
    } else {
      amrex::Real tbv = 0.;
      if (std::abs(yn) < 1.) {
        // This minimum value of tb comes from a maximum phi value
        amrex::Real pv1 = (yn - Wer);
        amrex::Real pv2 = -ydotn / omega;
        amrex::Real phi = std::atan2(pv2, pv1);
        if (phi < 0.) {
          phi = 2. * M_PI + phi;
        }
        amrex::Real theta = std::acos((1. - Wer) / A);
        // Minimum tb occurs at a minimum theta that is still above phi
        if (theta < phi) {
          if (2. * M_PI - theta >= phi) {
            theta = -theta;
          }
          theta += 2. * M_PI;
        }
        tbv = (theta - phi) / omega;
      }
      // Check if breakup occurs
      if (tbv < subdt) {
        const amrex::Real Wet = 80.;
        const amrex::Real k2 = 2. / 9.;
        const amrex::Real k1 =
          k2 * ((std::sqrt(Wet) - 1.) * std::pow(We / Wet, 4) + 1.);
        amrex::Real Kbr = k1 * omega;
        if (We >= Wet) {
          Kbr = k2 * omega * std::sqrt(We);
        }
        amrex::Real etheta = amrex::max(-1., amrex::min(1., 1. - 1. / Wer));
        amrex::Real etb = std::acos(etheta) / omega;
        amrex::Real r32 = rad_part * std::exp(-Kbr * etb);
        amrex::Real rsmr = std::sqrt(std::pow(rad_part, 3) / (r32 * r32));
        amrex::Real AE2 =
          3. * (1. - rad_part / rsmr + 5. * C_D * We / 72.) * omega2;
        amrex::Real AE = std::sqrt(AE2);
        amrex::Real Utan = AE * C_b * rad_part;
        Utan_total += Utan;
        num_dens *= std::pow(rad_part / r32, 3);
        rad_part = r32;
        yn = 0.;
        ydotn = 0.;
        if (curt + tbv < breakup_time) {
          breakup_time = curt + tbv;
        }
        tbest = rad_part * tbconst;
        subdt = dt / (0.1 * tbest);
      } else {
        yn = ynp;
        ydotn = ydotnp;
      }
    }
    curt += subdt;
  }
  p.rdata(SprayComps::pstateDia) = 2. * rad_part;
  p.rdata(SprayComps::pstateNumDens) = num_dens;
  p.rdata(SprayComps::pstateTABY) = yn;
  p.rdata(SprayComps::pstateTABYdot) = ydotn;
  return Utan_total;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
splitDroplet(
  const int pid,
  SprayParticleContainer::ParticleType& p,
  const SprayData& fdat,
  splash_breakup* N_SB,
  const SBPtrs& rf,
  const amrex::Real rem_dt,
  const amrex::Real Utan_total)
{
  amrex::Real num_dens = p.rdata(SprayComps::pstateNumDens);
  if (num_dens > fdat.max_num_ppp) {
    N_SB[pid] = splash_breakup::breakup;
    amrex::RealVect vel_part(AMREX_D_DECL(
      p.rdata(SprayComps::pstateVel), p.rdata(SprayComps::pstateVel + 1),
      p.rdata(SprayComps::pstateVel + 2)));
    amrex::Real velMag = vel_part.vectorLength();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      rf.loc[AMREX_SPACEDIM * pid + dir] = p.pos(dir);
      rf.vel[AMREX_SPACEDIM * pid + dir] = vel_part[dir];
      rf.norm[AMREX_SPACEDIM * pid + dir] = vel_part[dir] / velMag;
    }
    rf.numDens[pid] = p.rdata(SprayComps::pstateNumDens);
    rf.d0[pid] = p.rdata(SprayComps::pstateDia);
    rf.dtpp[pid] = rem_dt;
    rf.phi1[pid] = Utan_total;
    rf.phi2[pid] = p.rdata(SprayComps::pstateTABY);
    rf.T0[pid] = p.rdata(SprayComps::pstateT);
#if SPRAY_FUEL_NUM > 1
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      rf.Y0[SPRAY_FUEL_NUM * pid + spf] = p.rdata(SprayComps::pstateY + spf);
    }
#endif
    p.id() = -1;
  }
}
#endif
