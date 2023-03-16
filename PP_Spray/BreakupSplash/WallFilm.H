
#ifndef WALLFILM_H
#define WALLFILM_H

#include "Drag.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
calculateFilmSource(
  const amrex::Real flow_dt,
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
  amrex::Real theta_c = fdat.theta_c;
  amrex::Real tanthc = std::tan(theta_c);
  amrex::Real film_height = p.rdata(SprayComps::pstateFilmHght);
  amrex::Real film_dia = p.rdata(SprayComps::pstateDia);
  amrex::Real rad2 = std::pow(0.5 * film_dia, 2);
  amrex::Real hght2 = std::pow(film_height, 2);
  amrex::Real film_mass =
    M_PI / 6. * rho_film * film_height * (3. * rad2 + hght2);
  // Film surface area exposed to gas
  amrex::Real top_area = M_PI * (rad2 + hght2);
  // Film area exposed to the wall
  amrex::Real bottom_area = M_PI * rad2;
  amrex::Real start_mass = film_mass;
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
  amrex::Real rho_skin = mw_skin * gpv.p_fluid /
                         (pele::physics::Constants::RU * T_skin * SPU.ru_conv);
  amrex::Real rho_cgs = rho_skin / SPU.rho_conv;
  auto trans = pele::physics::PhysicsType::transport();
  trans.transport(
    get_xi, get_mu, get_lambda, get_Ddiag, get_chi, T_skin, rho_cgs,
    Y_skin.data(), Ddiag.data(), nullptr, mu_skin, xi_skin, lambda_skin,
    trans_parm);
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
      amrex::Real X_norm = Y_film[spf] * gpv.invmw[fspec] * Psat[spf] / sumXVap;
      Ddiag[fspec] *= X_norm;
      rhoDtotal += Ddiag[fspec];
    }
    amrex::Real Sc_skin = mu_skin / rhoDtotal;
    amrex::Real q_conv = top_area * h_conv * (gpv.T_fluid - T_film);
    amrex::Real q_cond =
      -bottom_area * lambda_film * 2. * (T_film - fdat.wall_T) / film_height;
    amrex::Real gstar_coeff = top_area * h_conv / lambda_skin *
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
        Y_film[spf] > 1.E-12 && Psat[spf] > gpv.X_fluid[fdspec] * gpv.p_fluid) {
        mi_dot[spf] = -amrex::max(gstar_coeff * Ddiag[fspec], 0.);
        m_dot += mi_dot[spf];
        sumL += mi_dot[spf] * L_fuel[spf];
      }
    }
    amrex::Real part_temp_src =
      (q_conv + q_cond + sumL) / (film_mass * cp_film);
    gpv.fluid_eng_src += q_conv;
    T_film += dt * part_temp_src;
    amrex::Real new_mass = film_mass + dt * m_dot;
    // Compute the new height assuming H = 0.5 * dia * tan(theta_c)
    film_height = std::cbrt(
      6. * new_mass / (rho_film * M_PI * (3. / (tanthc * tanthc) + 1.)));
    film_dia = 2. * film_height / tanthc;
    if (film_height > min_height) {
      if (SPRAY_FUEL_NUM > 1) {
        rho_film = 0.;
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          Y_film[spf] = (Y_film[spf] * film_mass + mi_dot[spf] * dt) / new_mass;
          rho_film += Y_film[spf] / fdat.rhoL(T_film, spf);
        }
        rho_film = 1. / rho_film;
      } else {
        rho_film = fdat.rhoL(T_film, 0);
      }
    } else {
      new_mass = 0.;
      p.id() = -1;
    }
    film_mass = new_mass;
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
  p.rdata(SprayComps::pstateDia) = film_dia;
  p.rdata(SprayComps::pstateFilmHght) = film_height;
}
#endif
