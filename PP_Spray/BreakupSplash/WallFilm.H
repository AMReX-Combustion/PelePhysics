
#ifndef WALLFILM_H
#define WALLFILM_H

#include "Drag.H"

// Determine the wall film thickness within the cell by summing film volume
// divided by the face area
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fillFilmFab(
  amrex::Array4<amrex::Real> const& wf_arr,
  SprayParticleContainer::ParticleyType& p,
  const amrex::Real& face_area,
  const amrex::RealVect& plo,
  const amrex::RealVect& dx)
{
  amrex::RealVect lxc = (p.pos() - plo) / dx;
  amrex::IntVect ijkc = lxc.floor(); // Cell with particle
  amrex::Real film_hght = p.rdata(SprayComps::pstateFilmHght);
  amrex::Real film_dia = p.rdata(SprayComps::pstateDia);
  amrex::Real film_vol =
    M_PI * film_hght * (0.75 * film_dia * film_dia + film_hght * film_hght);
  amrex::Gpu::Atomic::Add(&wf_arr(ijkc, 0), film_vol / face_area);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
bool
interpolateFilm(
  const amrex::IntVect& ijkc,
  const amrex::IntVect& bflags,
  amrex::RealVect& film_normal,
#ifdef AMREX_USE_EB
  bool eb_in_box,
  amrex::Array4<const amrex::EBCellFlag> const& flags,
  amrex::Array4<const amrex::Real> const& bnorm,
#endif
  amrex::IntVect* indx_array,
  amrex::Real* weights)
{
  bool do_fe_interp = false;
#ifdef AMREX_USE_EB
  if (eb_in_box) {
    if (flags(ijkc).isSingleValued()) {
      do_fe_interp = true;
      film_normal = {AMREX_D_DECL(bnorm(ijkc, 0), bnorm(ijkc, 1), bnorm(ijkc, 2))};
    }
  } else
#endif
  {
    amrex::Real magsum = 0.;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      if (bflags[dir] < 0) {
        film_normal[dir] = 1.;
        magsum += 1.;
      } else if (bflags[dir] > 0) {
        film_normal[dir] = -1.;
        magsum += 1.;
      }
    }
    film_normal /= std::sqrt(magsum);
  }
  // TODO: Consider providing more indices in case in a small vol cell
  for (int aindx = 0; aindx < AMREX_D_PICK(2, 4, 8); ++aindx) {
    indx_array[aindx] = ijkc;
    weights[aindx] = 0.;
  }
  weights[0] = 1.;
  return do_fe_interp;
}

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
  bool get_mu = false;
  bool get_chi = false;
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
  amrex::Real Tcrit = 0.;
  amrex::Real cp_film = 0.;
  amrex::Real lambda_film = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    amrex::Real minT = amrex::min(T_film, cBoilT[spf]);
    Y_film[spf] = p.rdata(SprayComps::pstateY + spf);
    rho_film += Y_film[spf] / fdat.rhoL(minT, spf);
    cp_film += Y_film[spf] * fdat.cp[spf];
    lambda_film += Y_film[spf] * fdat.lambdaL(minT, spf);
    Tcrit += Y_film[spf] * fdat.critT[spf];
  }
  rho_film = 1. / rho_film;
  T_film = amrex::min(0.999 * Tcrit, T_film);
  amrex::Real film_height = p.rdata(SprayComps::pstateFilmHght);
  amrex::Real film_dia = p.rdata(SprayComps::pstateDia);
  amrex::Real rad2 = std::pow(0.5 * film_dia, 2);
  // Surface area assuming film is a cylinder
  amrex::Real film_area = M_PI * rad2;
  amrex::Real film_mass = rho_film * film_height * film_area;
  amrex::Real start_mass = film_mass;
  // Model the fuel vapor using the one-third rule
  amrex::Real delT = amrex::max(gpv.T_fluid - T_film, 0.);
  amrex::Real T_skin = T_film + rule * delT;

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
  lambda_skin *= SPU.lambda_conv;
  // If gas phase is not saturated
  if (sumXVap > 0.) {
    // Estimate convective heat transfer using O'Rourke and Amsden 1996 assuming
    // y+ < y+c
    // TODO: This is very rudimentary and should be improved to account for gas
    // phase velocity
    amrex::Real h_heat = lambda_skin / film_height;
    amrex::Real q_conv = film_area * h_heat * (gpv.T_fluid - T_film);
    amrex::Real q_cond =
      -film_area * lambda_film * 2. * (T_film - fdat.wall_T) / film_height;
    // Use Chilton-Colburn analogy to find the convective mass transfer
    // coefficient; like droplet evaporation, this uses mixture averaged
    // values for lambda, cp, and B_M but uses normalized species mass
    // diffusivity
    amrex::Real h_mass = film_area * h_heat * std::cbrt(lambda_skin / cp_skin) *
                         std::log1p(B_M) / lambda_skin;
    amrex::Real sumL = 0.;
    amrex::Real m_dot = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      // Species index
      const int fspec = fdat.indx[spf];
      // Gas species index, might be the same
      const int fdspec = fdat.dep_indx[spf];
      // Convert mass diffusion coefficient from mixture average
      // to binary for fuel only, not concerned with other species
      Ddiag[fspec] *= mw_skin / gpv.mw[fspec] * SPU.rhod_conv;
      if (X_vapor[spf] > 0.) {
        // Normalize mass diffusivity by fuel vapor molar fraction
        amrex::Real cur_rhoD =
          Ddiag[fspec] * Y_film[spf] * gpv.invmw[fspec] * Psat[spf] / sumXVap;
        mi_dot[spf] = -amrex::max(h_mass * std::pow(cur_rhoD, 2. / 3.), 0.);
        m_dot += mi_dot[spf];
        sumL += mi_dot[spf] * L_fuel[spf];
      }
    }
    amrex::Real part_temp_src =
      (q_conv + q_cond + sumL) / (film_mass * cp_film);
    gpv.fluid_eng_src += q_conv;
    T_film += dt * part_temp_src;
    amrex::Real new_mass = film_mass + dt * m_dot;
    // Assumes mass is only lost in the wall normal direction and the diameter
    // remains constant; this is not physically correct
    film_height = new_mass / (rho_film * film_area);
    if (film_height > min_height) {
      if (SPRAY_FUEL_NUM > 1) {
        rho_film = 0.;
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          Y_film[spf] = amrex::min(
            1., amrex::max(
                  (Y_film[spf] * film_mass + mi_dot[spf] * dt) / new_mass, 0.));
          if (Y_film[spf] < 1.E-12) {
            Y_film[spf] = 0.;
          }
          rho_film += Y_film[spf] / fdat.rhoL(T_film, spf);
        }
        rho_film = 1. / rho_film;
      } else {
        rho_film = fdat.rhoL(T_film, 0);
      }
    } else {
      film_height = 0.;
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
