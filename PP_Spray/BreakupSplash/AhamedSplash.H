
#ifndef AHAMEDSPLASH_H
#define AHAMEDSPLASH_H

#include "SprayFuelData.H"
#include "SBData.H"

/****************************************************************
 Functions for imposing wall boundaries including adding
 splashed droplets, modifying particles to be wall film,
 and rebounding particles. These occur on the host
 ***************************************************************/

enum splash_type { rebound = 0, deposit, splash, thermal_breakup };

// Determine the type of particle impact at the wall
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
splash_type
splash_criteria(
  const amrex::Real Kv,
  const amrex::Real Tstar,
  const amrex::Real alpha,
  const amrex::Real Kcrit)
{
  splash_type crit = splash_type::thermal_breakup;
  if (Tstar < 1.1) {
    if (Kv < Kcrit) {
      crit = splash_type::deposit;
    } else {
      crit = splash_type::splash;
    }
  } else {
    if (Kv < Kcrit) {
      crit = splash_type::rebound;
    }
  }
  return crit;
}

// Determine if droplet is splashing, fill SB data if it is
// TODO: Add wet wall droplet splashing functions
AMREX_GPU_DEVICE
AMREX_INLINE
void
droplet_splashing(
  SprayParticleContainer::ParticleType& p,
  int pid,
  const amrex::RealVect& dx,
  const amrex::RealVect& plo,
  const SprayData& fdat,
  const amrex::Real p_gas,
  amrex::RealVect& vel_part,
  const amrex::Real par_dot,
  const amrex::RealVect normal,
  const amrex::Real* cBoilT,
  splash_breakup* N_SB,
  const SBPtrs& rf,
  const bool dry_wall,
  const amrex::Real rem_dt)
{
  amrex::Real Nw_Vp = normal.dotProduct(vel_part);
  amrex::Real T_part = p.rdata(SprayComps::pstateT);
  amrex::Real mu_part = 0.;
  amrex::Real rho_part = 0.;
  // TODO: Determine correct method for handling multi-component liquids
  amrex::Real Tboil = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int ysp = SprayComps::pstateY + spf;
    mu_part += p.rdata(ysp) * fdat.mu[spf];
    rho_part += p.rdata(ysp) / fdat.rhoL(T_part, spf);
    Tboil += p.rdata(ysp) * cBoilT[spf];
  }
  rho_part = 1. / rho_part;
  amrex::Real Tstar = fdat.wall_T / Tboil;
  const amrex::Real dia_part = p.rdata(SprayComps::pstateDia);
  const amrex::Real pmass = M_PI / 6. * rho_part * std::pow(dia_part, 3);
  // Weber number
  const amrex::Real We = rho_part * dia_part * Nw_Vp * Nw_Vp / fdat.sigma;
  const amrex::Real Re_L =
    amrex::Math::abs(Nw_Vp) * dia_part * rho_part / mu_part;
  const amrex::Real Kv = std::sqrt(We * std::sqrt(Re_L));
  amrex::Real umag = vel_part.vectorLength();
  // Inclination of on-coming particle to surface in radians
  amrex::Real alpha =
    amrex::max(M_PI / 6., std::asin(amrex::Math::abs(Nw_Vp) / umag));
#ifdef PELELM_USE_SPRAY
  const amrex::Real p_bar = p_gas * 1.E-5;
#else
  const amrex::Real p_bar = p_gas * 1.E-6;
#endif
  const amrex::Real Kcrit = (26.788 * std::log(p_bar) + 33.546) *
                              std::exp(-6.E-3 * (fdat.wall_T - 260.)) +
                            36.;
  splash_type splash_flag = splash_criteria(Kv, Tstar, alpha, Kcrit);
  // Amount of mass that gets deposited
  amrex::Real dm = 0.;
  if (splash_flag == splash_type::rebound) {
    amrex::Real rbound = 0.823 * std::exp(-1.9835E-3 * std::pow(Kv, 1.6));
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      p.pos(dir) -= (1. + rbound) * par_dot * vel_part[dir] / Nw_Vp;
      vel_part[dir] -= Nw_Vp * (1. + rbound) * normal[dir];
    }
  } else if (splash_flag == splash_type::deposit) {
    dm = pmass;
  } else {
    amrex::Real ms;
    // If splashing, turn current particle into wall film
    if (splash_flag == splash_type::splash) {
      N_SB[pid] = splash_breakup::splash_splash;
      ms = pmass * 0.75 * (1. - std::exp(-3.E-3 * (Kv - Kcrit)));
      dm = pmass - ms;
    } else {
      // For thermal breakup, no deposition occurs and current particle is
      // deleted
      N_SB[pid] = splash_breakup::splash_thermal_breakup;
      ms = pmass * 0.75 * (1. - std::exp(-0.0175 * (Kv - Kcrit)));
      p.id() = -1;
    }
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      // Save location where droplet hit wall
      rf.loc[AMREX_SPACEDIM * pid + dir] =
        (p.pos(dir) - par_dot * vel_part[dir] / Nw_Vp) * dx[dir] + plo[dir];
      rf.norm[AMREX_SPACEDIM * pid + dir] = normal[dir];
      rf.vel[AMREX_SPACEDIM * pid + dir] = vel_part[dir];
    }
    rf.d0[pid] = p.rdata(SprayComps::pstateDia);
    rf.dtpp[pid] = par_dot * dx[0] / Nw_Vp + rem_dt;
    rf.phi1[pid] = Kv;
    rf.phi2[pid] = ms;
    rf.T0[pid] = p.rdata(SprayComps::pstateT);
#if SPRAY_FUEL_NUM > 1
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      rf.Y0[SPRAY_FUEL_NUM * pid + spf] = p.rdata(SprayComps::pstateY + spf);
    }
#endif
  }
  if (dm > 0.) {
    // Convert parcel to wall film
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      p.pos(dir) -= 1.1 * par_dot * vel_part[dir] / Nw_Vp;
    }
    vel_part = amrex::RealVect::TheZeroVector();
    p.rdata(SprayComps::pstateFilmVol) = dm / rho_part;
  }
}

#endif
