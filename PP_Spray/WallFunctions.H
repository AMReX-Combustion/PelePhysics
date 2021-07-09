#ifndef _WALLFUNCTIONS_H_
#define _WALLFUNCTIONS_H_

#include "Drag.H"
#include "SprayFuelData.H"
#include "SprayInterpolation.H"

using namespace amrex;

// Check if tile is adjacent to non-periodic boundaries
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
bool
tile_at_bndry(
  const Box& in_box,
  const IntVect& bndry_lo,
  const IntVect& bndry_hi,
  const Box& domain)
{
  Box testBox(in_box);
  testBox.grow(1);
  if (domain.contains(testBox))
    return false;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    Box boxlo(in_box);
    boxlo.growLo(dir, 1);
    Box boxhi(in_box);
    boxhi.growHi(dir, 1);
    if (!domain.contains(boxlo) && bndry_lo[dir] != 0) {
      return true;
    } else if (!domain.contains(boxhi) && bndry_hi[dir] != 0) {
      return true;
    }
  }
  return false;
}

/****************************************************************
 Functions computing wall film source terms
 ***************************************************************/

AMREX_GPU_DEVICE AMREX_INLINE void
calculateWallFilmSource(
  const Real flow_dt,
  GasPhaseVals& gpv,
  const SprayComps& SPI,
  const SprayData& fdat,
  SprayParticleContainer::ParticleType& p,
#ifdef USE_SPRAY_SOA
  const std::array<SprayParticleContainer::SoA::RealVector, NAR_SPR>& attribs,
  const int pid,
#endif
  const Real T_wall,
  const Real face_area,
  const Real diff_cent,
  pele::physics::transport::TransParm const* trans_parm)
{
  auto eos = pele::physics::PhysicsType::eos();
  SprayUnits SPU;
  const Real C_eps = 1.E-15;
  const Real B_eps = 1.E-7;
  const Real ht_tol = 2.E-6;
  const Real part_dt = 0.5 * flow_dt;
  bool get_xi = false;
  bool get_Ddiag = true;
  bool get_lambda = true;
  bool get_mu = true;
  GpuArray<Real, SPRAY_FUEL_NUM> Y_film;
  GpuArray<Real, SPRAY_FUEL_NUM> Y_vapor;
  GpuArray<Real, SPRAY_FUEL_NUM> L_fuel;
  GpuArray<Real, SPRAY_FUEL_NUM> mi_dot;
  GpuArray<Real, SPRAY_FUEL_NUM> cBoilT;
  GpuArray<Real, NUM_SPECIES> Y_skin;
  GpuArray<Real, NUM_SPECIES> h_part;
  GpuArray<Real, NUM_SPECIES> Ddiag;
  for (int sp = 0; sp < NUM_SPECIES; ++sp)
    Y_skin[sp] = 0.;
    // Get particle variables
    // Note: the diameter is set assuming a sphere
    // of the same volume as the film
    // and temperature is set to a negative value
#ifdef USE_SPRAY_SOA
  Real T_film = -attribs[SPI.pstateT].data()[pid];
  Real vol = attribs[SPI.pstateVol].data()[pid];
  Real ht_film = attribs[SPI.pstateHt].data()[pid];
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
    Y_film[spf] = attribs[SPI.pstateY + spf].data()[pid];
#else
  Real T_film = -p.rdata(SPI.pstateT);
  Real vol = p.rdata(SPI.pstateVol);
  Real ht_film = p.rdata(SPI.pstateHt);
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
    Y_film[spf] = p.rdata(SPI.pstateY + spf);
#endif
  calcBoilT(fdat, gpv, cBoilT.data());
  // TODO: Add model for when T > T_boil
  Real T_i = gpv.T_fluid;
  Real area_film = vol / ht_film;
  Real dia_film = std::sqrt(4. * area_film / M_PI);
  eos.T2Hi(T_film, h_skin.data());
  for (int n = 0; n < NUM_SPECIES; ++n) {
    h_part[n] *= SPU.eng_conv;
  }
  calcVaporY(
    fdat, gpv, T_film, C_eps, Y_film.data(), h_part.data(), cBoilT.data(),
    Y_vapor.data(), L_fuel.data());
  Real T_vapor = 0.5 * (T_film + T_i);
  Real sumYSkin = 0.;
  Real sumYFuel = 0.;
  Real rho_film = 0.;
  Real lambda_film = 0.;
  Real cp_film = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fspec = fdat.indx[spf];
    rho_film += Y_film[spf] / fdat.rho[spf];
    lambda_film += Y_film[spf] * fdat.lambda[spf];
    cp_film += Y_film[spf] * fdat.cp[spf];
    Y_skin[fspec] = 0.5 * (Y_vapor[spf] + gpv.Y_fluid[fspec]);
    sumYSkin += Y_skin[fspec];
    sumYFuel += gpv.Y_fluid[fspec];
  }
  rho_film = 1. / rho_film;
  const Real restYSkin = 1. - sumYSkin;
  const Real restYfluid = 1. - sumYFuel;
  Real renorm = restYSkin / restYfluid;
  Real mw_vap = 0.;
  for (int sp = 0; sp < NUM_SPECIES; ++sp) {
    if (Y_skin[sp] == 0.) {
      Y_skin[sp] = gpv.Y_fluid[sp] * renorm;
    }
    mw_vap += Y_skin[sp] * gpv.invmw[sp];
  }
  mw_vap = 1. / mw_vap;
  Real pmass = vol * rho_film;
  Real lambda_skin = 0.;
  Real mu_skin = 0.;
  Real xi_skin = 0.;
  auto trans = pele::physics::PhysicsType::transport();
  trans.transport(
    get_xi, get_mu, get_lambda, get_Ddiag, T_vapor, gpv.rho_fluid,
    Y_skin.data(), Ddiag.data(), mu_skin, xi_skin, lambda_skin, trans_parm);
  mu_skin *= SPU.mu_conv;
  lambda_skin *= SPU.lambda_conv;
  Real dy_i = diff_cent - ht_film; // Distance from film surface to cell center
  // Determine the mass evaporation values
  Real m_dot = 0.;
  Real qvap = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fspec = fdat.indx[spf];
    Ddiag[fspec] *= mw_vap * gpv.invmw[fspec] * SPU.rhod_conv;
    // TODO: Ensure condensed mass does not exceed gas phase mass
    mi_dot[spf] = -Ddiag[fspec] / (1. - Y_vapor[spf]) *
                  (gpv.Y_fluid[fspec] - Y_vapor[spf]) / dy_i;
    m_dot += mi_dot[spf];
    qvap += mi_dot[spf] * L_fuel[spf];
    gpv.fluid_Y_dot[spf] += mi_dot[spf];
    gpv.fluid_eng_src += mi_dot[spf] * h_skin[fspec];
  }
  gpv.fluid_mass_src = m_dot;
  // Determine the temperature at wall film surface using energy balance
  Real fs1 = lambda_skin * ht_film;
  Real fs2 = lambda_film * dy_i;
  Real T_s = (fs1 * T_i + fs2 * T_wall - ht_film * dy_i * qvap) / (fs1 + fs2);
  Real qconv = lambda_skin * (T_i - T_s) / dy_i;
  gpv.fluid_eng_src += qconv;
  // Update film height
  Real new_mass = pmass + part_dt * m_dot;
  Real new_rho = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Real newY = (Y_film[spf] * pmass + mi_dot[spf] * part_dt) / new_mass;
    p.rdata(SPI.pstateY + spf) = newY;
    new_rho += newY / fdat.rho[spf];
  }
  new_rho = 1. / new_rho;
  Real new_vol = new_mass / new_rho;
  p.rdata(SPI.pstateVol) = new_vol;
  // Adjust the height by assuming radius remains unchanged
  Real new_ht = 4. * new_vol / (M_PI * dia_film * dia_film);
  p.rdata(SPI.pstateHt) = new_ht;
  p.rdata(SPI.pstateT) = -0.5 * (T_s + T_wall);
  if (new_ht < ht_tol)
    p.id() = -1;
}

/****************************************************************
 Functions for imposing wall boundaries including adding
 splashed droplets, modifying particles to be wall film,
 and rebounding particles. These occur on the host
 ***************************************************************/

enum splash_type {
  rebound = 0,
  deposit,
  splash,
  thermal_breakup,
  no_impact,
  wall_film
};

// Determine the type of particle impact at the wall
// Return 0 - rebound, 1 - deposit, 2 - splash, 3 - thermal breakup
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
splash_type
splash_criteria(const Real Kv, const Real Tstar, const Real alpha)
{
  Real Kcrit = 20. + 2. * alpha / M_PI * 20.;
  if (Tstar < 1.1) {
    Kcrit = 130.;
    if (Tstar < 1.)
      Kcrit = 54. + 76. * std::exp(13. * (Tstar - 1.));
    if (Kv < Kcrit)
      return splash_type::deposit;
    else
      return splash_type::splash;
  }
  if (Kv < Kcrit)
    return splash_type::rebound;
  return splash_type::thermal_breakup;
}

// Find tangents along surface
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
void
find_tangents(
  const RealVect& testvec,
  const RealVect& norm,
  RealVect& tanBeta,
  RealVect& tanPsi)
{
#if AMREX_SPACEDIM == 3
  tanPsi = testvec.crossProduct(norm);
  tanBeta = tanPsi.crossProduct(norm);
  tanPsi /= tanPsi.vectorLength();
  tanBeta /= tanBeta.vectorLength();
#else
  tanBeta[0] = -norm[1];
  tanBeta[1] = norm[0];
#endif
}

// Compute the angles of reflecting droplets
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
void
compute_angles(
  const Real& alpha,
  const Real& Tstar,
  const Real& We,
  const bool dry,
  SprayRefl& SPRF)
{
  // Inclination angle in degrees
  const Real alphad = alpha * 180. / M_PI;
  Real omega = 0.;
  if (alphad <= 80. && AMREX_SPACEDIM == 3)
    omega = std::sqrt(
      (1. + 8.872 * std::cos(1.152 * alpha)) / (1. - std::cos(alpha)));
  SPRF.omega = omega;
  SPRF.expomega = 1. - std::exp(-omega);
  Real beta_mean;
  if (dry)
    beta_mean = 9.3 + 0.22 * alphad;
  else if (Tstar > 1.1)
    beta_mean = 0.225 * alphad * std::exp(std::pow(0.017 * alphad - 0.937, 2));
  else
    beta_mean = 0.96 * alphad * std::exp(-4.5E-3 * We);
  Real stdev = 4.; // 4 degrees
  // Now convert mean and stdev to log terms
  Real mean2 = beta_mean * beta_mean;
  Real st2 = stdev * stdev;
  Real term1 = std::log(beta_mean);
  Real term2 = std::log(mean2 + st2);
  SPRF.beta_mean = 2. * term1 - 0.5 * term2;
  SPRF.beta_stdv = std::sqrt(amrex::max(-2. * term1 + term2, 0.));
}

AMREX_GPU_HOST_DEVICE
AMREX_INLINE
bool
check_wall(
  SprayParticleContainer::ParticleType& p,
  const IntVect bflags,
  const IntVect ijkc,
  const IntVect ijkc_prev,
#ifdef AMREX_USE_EB
  const bool use_EB,
  Array4<EBCellFlag const> const& flags,
  Array4<Real const> const& bcent,
  Array4<Real const> const& bnorm,
#endif
  IntVect& bloc,
  RealVect& normal,
  RealVect& bcentv)
{
  bool wall_check = false; // Check if particle is reflected
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (bflags[dir] == -1 || bflags[dir] == 1) {
      normal[dir] = -bflags[dir];
      bcentv[dir] = -0.5 * bflags[dir];
      wall_check = true;
    }
  }
#ifdef AMREX_USE_EB
  if (use_EB && !wall_check) {
    // If particle has moved into a cut-cell
    if (flags(ijkc).isSingleValued()) {
      // Use the normal and boundary centroid from the new cell
      normal = {-bnorm(ijkc, 0), -bnorm(ijkc, 1), -bnorm(ijkc, 2)};
      bcentv = {bcent(ijkc, 0), bcent(ijkc, 1), bcent(ijkc, 2)};
      wall_check = true;
      // If the particle has moved into a covered cell
    } else if (flags(ijkc).isCovered()) {
      // If the cell the particle came from is a cut-cell
      if (flags(ijkc_prev).isSingleValued()) {
        // Use the normal and boundary centroid from the previous cell
        bloc = ijkc_prev;
        normal = {
          -bnorm(ijkc_prev, 0), -bnorm(ijkc_prev, 1), -bnorm(ijkc_prev, 2)};
        bcentv = {
          bcent(ijkc_prev, 0), bcent(ijkc_prev, 1), bcent(ijkc_prev, 2)};
        // otherwise an error has occurred
      } else {
        Abort("WallFunctions::check_wall: Particle is outside EB");
      }
      wall_check = true;
    }
  }
#endif
  return wall_check;
}

AMREX_GPU_HOST_DEVICE
AMREX_INLINE
splash_type
impose_wall(
  SprayParticleContainer::ParticleType& p,
  SprayComps SPI,
  SprayData fdat,
  const RealVect& dx,
  const RealVect& plo,
  const RealVect& phi,
  const Real T_wall,
  const IntVect bloc,
  const RealVect normal,
  const RealVect bcentv,
  SprayRefl& SPRF,
  bool isActive,
  const bool dry_wall,
  amrex::RandomEngine const& engine)
{
  const Real tolerance = std::numeric_limits<Real>::epsilon();
  Real sigma = fdat.sigma;
  splash_type splash_flag = splash_type::no_impact;
  RealVect pvel = {AMREX_D_DECL(
    p.rdata(SPI.pstateVel), p.rdata(SPI.pstateVel + 1),
    p.rdata(SPI.pstateVel + 2))};
  // Projection of vector pointing from EB centroid to particle onto EB normal
  const Real par_dot_EB = AMREX_D_TERM(
    (p.pos(0) - (bloc[0] + 0.5 + bcentv[0]) * dx[0] - plo[0]) * normal[0],
    +(p.pos(1) - (bloc[1] + 0.5 + bcentv[1]) * dx[1] - plo[1]) * normal[1],
    +(p.pos(2) - (bloc[2] + 0.5 + bcentv[2]) * dx[2] - plo[2]) * normal[2]);
  // Check if particle is on the wrong side of the EB
  if (par_dot_EB < tolerance) {
    splash_flag = splash_type::rebound;
    // Plane ref point
    Real Nw_Vp = 0.;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
      Nw_Vp += normal[dir] * pvel[dir];
    // Parcel normal velocity
    RealVect Vpn = {
      AMREX_D_DECL(Nw_Vp * normal[0], Nw_Vp * normal[1], Nw_Vp * normal[2])};
    // Parcel tangential velocity
    RealVect Vpt = pvel - Vpn;
    // If using a splash model
    if (sigma > 0. && isActive) {
      Real mu_part = 0.;
      Real rho_part = 0.;
      // TODO: Determine correct method for handling multi-component liquids
      Real Tstar = 0.; // Average sum Y_i T_wall/T_boil,i
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        const int ysp = SPI.pstateY + spf;
        mu_part += p.rdata(ysp) * fdat.mu[spf];
        rho_part += p.rdata(ysp) / fdat.rho[spf];
        Tstar += T_wall * p.rdata(ysp) / fdat.boilT[spf];
      }
      rho_part = 1. / rho_part;
      const Real dia_part = p.rdata(SPI.pstateDia);
      const Real d3 = std::pow(dia_part, 3);
      const Real pmass = M_PI / 6. * rho_part * d3;
      // Weber number
      const Real We = rho_part * dia_part * Nw_Vp * Nw_Vp / sigma;
      const Real Re_L = std::abs(Nw_Vp) * dia_part * rho_part / mu_part;
      const Real sqReyn = std::sqrt(Re_L);
      const Real Kv = std::sqrt(We) * std::pow(Re_L, 0.25);
      Real umag = pvel.vectorLength();
      // Inclination of on-coming particle to surface in radians
      const Real alpha = std::asin(std::abs(Nw_Vp) / umag);
      splash_flag = splash_criteria(Kv, Tstar, alpha);
      // Boundary layer thickness of particle estimate provided by
      // Pasandideh-Fard et al. 1996
      const Real delta = 2. * dia_part / sqReyn;
      // splash_flag == 0, then particle rebounds, no breakup, no deposit
      if (splash_flag != splash_type::rebound) {
        // Energy loss due to viscious dissipation during impact
        Real E_diss = 4. * M_PI / 3. * rho_part * Nw_Vp * Nw_Vp * d3 / sqReyn;
        // If splash_flag == 1, particle is completely deposited onto the wall
        if (splash_flag == splash_type::deposit) {
          // TODO: Determine better way to estimate wall film temperature,
          // diameter, etc
          Real depot_vol = pmass / rho_part;
          // Assume deposited film is a dome shape of height delta, determine
          // the diameter
          Real depot_dia =
            2. * std::sqrt(
                   (6. * depot_vol / M_PI - std::pow(delta, 3)) / (3. * delta));
          // Now estimate the height to be a cylinder of diameter depot_dia
          Real depot_height = 4. * depot_vol / (M_PI * depot_dia * depot_dia);
          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
            p.rdata(SPI.pstateVel + dir) = 0.;
          p.rdata(SPI.pstateVol) = depot_vol;
          p.rdata(SPI.pstateHt) = depot_height;
        } else {
          // Determine diameter of secondary droplets
          Real expon = 3.6 * (alpha / M_PI) * (alpha / M_PI);
          if (dry_wall)
            SPRF.dia_refl =
              dia_part * 3.3 * std::exp(expon) * std::pow(We, -0.65);
          else
            SPRF.dia_refl =
              dia_part * 2.2 * std::exp(expon) * std::pow(We, -0.36);
          // Compute angles of reflecting droplets
          compute_angles(alpha, Tstar, We, dry_wall, SPRF);
          // If thermal breakup, all mass is converted to secondary droplets
          Real splash_mass = pmass;
          // If splash, some mass rebounds and some deposits
          if (splash_flag == splash_type::splash) {
            // Determine the fraction of the droplet mass
            // that forms secondary drops from Kuhnke 2004
            // TODO: This is incomplete by always assuming a dry wall
            const Real B = 0.2 + 0.6 * amrex::Random(engine);
            splash_mass *= amrex::min(1., (Tstar - 0.8) / 0.3 * (1. - B) + B);
            // Mass deposited into wall film
            Real depot_mass = pmass - splash_mass;
            Real depot_vol = depot_mass / rho_part;
            // Assume deposited film is a dome shape of height delta,
            // determine the diameter
            Real depot_dia =
              2. *
              std::sqrt(
                (6. * depot_vol / M_PI - std::pow(delta, 3)) / (3. * delta));
            // Now estimate the height to be a cylinder of diameter depot_dia
            Real depot_height = 4. * depot_vol / (M_PI * depot_dia * depot_dia);
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
              p.rdata(SPI.pstateVel + dir) = 0.;
            p.rdata(SPI.pstateVol) = depot_vol;
            p.rdata(SPI.pstateHt) = depot_height;
          } else {
            p.id() = -1;
          }
          Real mass_refl = M_PI / 6. * rho_part * std::pow(SPRF.dia_refl, 3);
          Real Ns = splash_mass / mass_refl;
          SPRF.Ns_refl = int(Ns);
          // TODO: Not sure what this variable is but it ranges from 1.45-2
          Real nu32 = 1.45;
          Real We_out =
            SPRF.dia_refl / dia_part *
              (We * (1. - 0.85 * std::pow(std::sin(alpha * M_PI / 180.), 2)) +
               12.) -
            12. / nu32;
          SPRF.Unorm = std::sqrt(sigma * We_out / (rho_part * SPRF.dia_refl));
          // This is the dt to move the secondary droplets
          SPRF.dt_pp = par_dot_EB / Nw_Vp;
          // Test vector for finding tangents
          RealVect testvec = -pvel;
          // Find two tangent vectors by taking the cross product with the
          // velocity vector
          find_tangents(testvec, normal, SPRF.tanBeta, SPRF.tanPsi);
        }
      } // if (splash_flag != splash_type::rebound)
    }   // if (sigma > 0.)
    // Reflect particle
    if (sigma < 0. || splash_flag == splash_type::rebound) {
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        p.pos(dir) -= 2. * par_dot_EB * normal[dir];
        p.rdata(SPI.pstateVel + dir) = -Vpn[dir] + Vpt[dir];
      }
    } else {
      // Or place particle at wall
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        SPRF.norm[dir] = normal[dir];
        SPRF.pos_refl[dir] -= par_dot_EB * normal[dir];
        p.pos(dir) -= (1. + 1.E-4) * par_dot_EB * normal[dir];
      }
    }
  } // if (par_dot_EB < tolerance)
  return splash_flag;
}

void
create_splash_droplet(
  SprayParticleContainer::ParticleType& p,
  SprayComps SPI,
  SprayRefl SPRF,
  SprayUnits SPU)
{
  Real rand1 = amrex::Random();
  Real mean = SPRF.beta_mean;
  Real stdev = SPRF.beta_stdv;
  Real beta = amrex::RandomNormal(mean, stdev);
  beta = std::exp(beta) * M_PI / 180.;
  Real omega = SPRF.omega;
  Real expb = SPRF.expomega;
  // If incoming angle is greater than 80 degrees, azimuthal
  // angle of outgoing drop can be between 0 to 2*pi
  Real psi = AMREX_D_PICK(0., 0., rand1 * 2. * M_PI);
  // Otherwise, follow we modify the recommendation by Naber and Reitz 1988,
  // so the azimuthal angle distribution favors the pre-splash drop path
  // as the inclination angle decreases
  if (omega > 0.) {
    Real rand2 = std::copysign(1., 0.5 - amrex::Random());
    psi = -rand2 / omega * std::log(1. - rand1 * expb) * M_PI;
  }
  Real costhetad = std::cos(beta);
  Real sinthetad = std::sin(beta);
  AMREX_D_TERM(Real un = SPRF.Unorm * sinthetad;
               , Real utBeta = SPRF.Unorm * costhetad * std::cos(psi);
               , Real utPsi = SPRF.Unorm * costhetad * std::sin(psi););
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    Real pvel = AMREX_D_TERM(
      un * SPRF.norm[dir], +utBeta * SPRF.tanBeta[dir],
      +utPsi * SPRF.tanPsi[dir]);
    Gpu::Atomic::Add(&p.pos(dir), SPRF.dt_pp * pvel);
    Gpu::Atomic::Add(&p.rdata(SPI.pstateVel + dir), pvel);
  }
}

#endif // _WALLFUNCTIONS_H_
