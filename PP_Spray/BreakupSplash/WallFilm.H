
#ifndef WALLFILM_H
#define WALLFILM_H

#include "Drag.H"

namespace WFIndx {
static const int Y = 0;
static const int Vol = Y + SPRAY_FUEL_NUM;
static const int T = Vol + 1;
static const int Id = T + 1;
static const int wf_num = Id + 1;
} // namespace WFIndx

// Fill wall film fab with sum of wall film parcel values
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fillFilmFab(
  amrex::Array4<amrex::Real> const& wall_film,
  const SprayData& fdat,
  SprayParticleContainer::ParticleType& p,
  amrex::RealVect plo,
  amrex::RealVect dx)
{
  amrex::RealVect lxc = (p.pos() - plo) / dx;
  amrex::IntVect ijkc = lxc.floor(); // Cell with particle
  amrex::Real T_part = p.rdata(SprayComps::pstateT);
  amrex::Real rho_film = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    rho_film += p.rdata(SprayComps::pstateY + spf) / fdat.rhoL(T_part, spf);
  }
  rho_film = 1. / rho_film;
  amrex::Real film_vol = p.rdata(SprayComps::pstateFilmVol);
  amrex::Real film_mass = film_vol * rho_film;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    amrex::Gpu::Atomic::Add(
      &wall_film(ijkc, WFIndx::Y + spf),
      p.rdata(SprayComps::pstateY + spf) * film_mass);
  }
  amrex::Gpu::Atomic::Add(&wall_film(ijkc, WFIndx::Vol), film_vol);
  amrex::Gpu::Atomic::Add(&wall_film(ijkc, WFIndx::T), film_mass * T_part);
  // Will be set by only one particle, that's the one we will reference in the
  // update
  wall_film(ijkc, WFIndx::Id) = static_cast<amrex::Real>(p.id());
}

// Fill one parcel with total wall film values at that cell, delete the other
// wall film parcels
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
getWallFilm(
  const amrex::IntVect ijkc,
  amrex::Array4<const amrex::Real> const& wall_film,
  SprayParticleContainer::ParticleType& p,
  amrex::RealVect dx)
{
  if (static_cast<amrex::Long>(wall_film(ijkc, WFIndx::Id)) == p.id()) {
    amrex::Real vol = wall_film(ijkc, WFIndx::Vol);
    // TODO: Allow this to change in case EB is used
    amrex::Real area = AMREX_D_TERM(1., *dx[0], *dx[0]);
    amrex::Real height = vol / area;
    amrex::Real wall_mass = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      wall_mass += wall_film(ijkc, WFIndx::Y + spf);
    }
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      p.rdata(SprayComps::pstateY + spf) =
        wall_film(ijkc, WFIndx::Y + spf) / wall_mass;
    }
    p.rdata(SprayComps::pstateFilmVol) = wall_film(ijkc, WFIndx::Vol);
    p.rdata(SprayComps::pstateDia) = height;
    p.rdata(SprayComps::pstateT) = wall_film(ijkc, WFIndx::T) / wall_mass;
  } else {
    p.id() = -1;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
calculateFilmSource(
  const amrex::Real flow_dt,
  const amrex::RealVect dx,
  GasPhaseVals& gpv,
  SprayData fdat,
  SprayParticleContainer::ParticleType& p,
  amrex::Real* cBoilT,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* trans_parm)
{
  auto eos = pele::physics::PhysicsType::eos();
  const amrex::Real dt = flow_dt * fdat.dtmod;
  SprayUnits SPU;
  const amrex::Real rule = 1. / 2.;
  amrex::Real C_eps = 1.E-15;
  amrex::Real min_height = 1.E-5 * SPU.len_conv;
  bool get_xi = false;
  bool get_Ddiag = true;
  bool get_lambda = true;
  bool get_mu = true;
  bool get_chi = false;
  if (!fdat.mass_trans) {
    get_Ddiag = false;
    get_lambda = false;
  }
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_skin;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> h_film;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> cp_n;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Ddiag;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> L_fuel = {{0.0}};
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> mi_dot = {{0.0}};
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_film; // Liquid mass fractions
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Psat = {{0.0}};
  amrex::Real T_film = p.rdata(SprayComps::pstateT);
  amrex::Real rho_film = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_film[spf] = p.rdata(SprayComps::pstateY + spf);
    rho_film += Y_film[spf] / fdat.rhoL(amrex::min(T_film, cBoilT[spf]), spf);
  }
  rho_film = 1. / rho_film;
  amrex::Real film_vol = p.rdata(SprayComps::pstateFilmVol);
  amrex::Real film_mass = film_vol * rho_film;
  amrex::Real start_mass = film_mass;
  amrex::Real film_height = p.rdata(SprayComps::pstateDia);
  // TODO: Allow this to change in case EB is used
  amrex::Real film_area = AMREX_D_TERM(1., *dx[0], *dx[0]);
  int isub = 1;
  int nsub = 1;
  while (isub <= nsub) {
    amrex::Real Tcrit = 0.;
    amrex::Real cp_film = 0.;
    amrex::Real lambda_film = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      cp_film += Y_film[spf] * fdat.cp[spf];
      lambda_film += Y_film[spf] * fdat.lambda[spf];
      Tcrit += Y_film[spf] * fdat.critT[spf];
    }
    T_film = amrex::min(0.999 * Tcrit, T_film);
    // Model the fuel vapor using the one-third rule
    amrex::Real delT = amrex::max(gpv.T_fluid - T_film, 0.);
    amrex::Real T_skin = gpv.T_fluid;
    if (fdat.mass_trans) {
      T_skin = T_film + rule * delT;
    }

    // Calculate the C_p at the skin temperature for each species
    eos.T2Cpi(T_skin, cp_n.data());
    eos.T2Hi(T_film, h_film.data());
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y_skin[n] = 0.;
      h_film[n] *= SPU.eng_conv;
      cp_n[n] *= SPU.eng_conv;
    }
    amrex::Real cp_skin = 0.; // Average C_p in modeled skin phase
    amrex::Real mw_skin = 0.; // Average molar mass of skin phase
    amrex::Real B_M = 0.;     // Mass Spalding number
    amrex::Real sumXVap = 0.; // Sum of Y_L Psat_f / mw_f
    calcVaporState(
      fdat, gpv, rule, T_film, C_eps, Y_film.data(), h_film.data(), cp_n.data(),
      cBoilT, Y_skin.data(), Psat.data(), L_fuel.data(), B_M, sumXVap, cp_skin,
      mw_skin);
    amrex::Real lambda_skin = 0.;
    amrex::Real mu_skin = 0.;
    amrex::Real xi_skin = 0.;
    amrex::Real rho_skin =
      mw_skin * gpv.p_fluid /
      (pele::physics::Constants::RU * T_skin * SPU.ru_conv);
    amrex::Real rho_cgs = rho_skin / SPU.rho_conv;
    auto trans = pele::physics::PhysicsType::transport();
    trans.transport(
      get_xi, get_mu, get_lambda, get_Ddiag, get_chi, T_skin, rho_cgs,
      Y_skin.data(), Ddiag.data(), mu_skin, xi_skin, lambda_skin, trans_parm);
    mu_skin *= SPU.mu_conv;
    lambda_skin *= SPU.lambda_conv;
    amrex::Real Pr_skin = mu_skin * cp_skin / lambda_skin;
    // Local Reynolds number
    amrex::Real vel_mag = gpv.vel_fluid.vectorLength();
    amrex::Real Reyn = rho_skin * vel_mag * film_height / mu_skin;
    amrex::Real Nu =
      amrex::max(1., 0.0287 * std::pow(Pr_skin, 0.6) * std::pow(Reyn, 0.8));
    // If gas phase is not saturated
    if (sumXVap > 0.) {
      // TODO: Determine a better way to estimate convective heat transfer
      amrex::Real refL = 0.5 * film_height;
      amrex::Real h_conv = Nu * lambda_skin / refL;
      amrex::Real rhoDtotal = 0.;
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        const int fspec = fdat.indx[spf];
        // Convert mass diffusion coefficient from mixture average
        // to binary for fuel only, not concerned with other species
        Ddiag[fspec] *= mw_skin * gpv.invmw[fspec] * SPU.rhod_conv;
        // Get normalized vapor molar fraction
        amrex::Real X_norm =
          Y_film[spf] * gpv.invmw[fspec] * Psat[spf] / sumXVap;
        Ddiag[fspec] *= X_norm;
        rhoDtotal += Ddiag[fspec];
      }
      amrex::Real Sc_skin = mu_skin / rhoDtotal;
      amrex::Real q_conv = film_area * h_conv * (gpv.T_fluid - T_film);
      amrex::Real q_cond =
        -film_area * lambda_film * 2. * (T_film - fdat.wall_T) / film_height;
      amrex::Real gstar_coeff = film_area * h_conv / lambda_skin *
                                std::pow(Sc_skin / Pr_skin, 0.6) *
                                std::log1p(B_M);
      amrex::Real sumL = 0.;
      amrex::Real m_dot = 0.;
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        // Species index
        const int fspec = fdat.indx[spf];
        // Gas species index, might be the same
        const int fdspec = fdat.dep_indx[spf];
        if (
          Y_film[spf] > 1.E-12 &&
          Psat[spf] > gpv.X_fluid[fdspec] * gpv.p_fluid) {
          mi_dot[spf] = -amrex::max(gstar_coeff * Ddiag[fspec], 0.);
          m_dot += mi_dot[spf];
          sumL += mi_dot[spf] * L_fuel[spf];
        }
      }
      if (isub == 1) {
        amrex::Real inv_tau_m = -m_dot / (0.2 * film_mass);
        nsub = amrex::max(nsub, static_cast<int>(dt * inv_tau_m));
      }
      amrex::Real part_temp_src =
        (q_conv + q_cond + sumL) / (film_mass * cp_film);
      gpv.fluid_eng_src += q_conv;
      auto film_dt = dt / static_cast<amrex::Real>(nsub);
      T_film += film_dt * part_temp_src;
      amrex::Real new_mass = film_mass + film_dt * m_dot;
      film_height = new_mass / (rho_film * film_area);
      if (film_height > min_height) {
        if (SPRAY_FUEL_NUM > 1) {
          rho_film = 0.;
          for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
            Y_film[spf] =
              (Y_film[spf] * film_mass + mi_dot[spf] * film_dt) / new_mass;
            rho_film += Y_film[spf] / fdat.rhoL(T_film, spf);
          }
          rho_film = 1. / rho_film;
        } else {
          rho_film = fdat.rhoL(T_film, 0);
        }
      } else {
        new_mass = 0.;
        p.id() = -1;
        nsub = isub;
        isub = nsub + 1;
      }
      film_mass = new_mass;
      isub++;
    }
    if (nsub > 1) {
      gpv.fluid_eng_src /= static_cast<amrex::Real>(nsub);
    }
    p.rdata(SprayComps::pstateT) = T_film;
    amrex::Real mdot_total = (film_mass - start_mass) / dt;
    gpv.fluid_mass_src = mdot_total;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      const int fdspec = fdat.dep_indx[spf];
      amrex::Real oldY = p.rdata(SprayComps::pstateY + spf);
      amrex::Real newY = Y_film[spf];
      amrex::Real midot = (newY * film_mass - oldY * start_mass) / dt;
      gpv.fluid_Y_dot[spf] = midot;
      p.rdata(SprayComps::pstateY + spf) = newY;
      gpv.fluid_eng_src += midot * h_film[fdspec];
    }
    p.rdata(SprayComps::pstateDia) = film_height;
    p.rdata(SprayComps::pstateFilmVol) = film_mass / rho_film;
  }
}
#endif
