#ifndef _WALLFUNCTIONS_H_
#define _WALLFUNCTIONS_H_

#include "Drag.H"
#include "SprayFuelData.H"
#include "SprayInterpolation.H"

// Check if tile is adjacent to non-periodic boundaries
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
bool
tile_at_bndry(
  const amrex::Box& in_box,
  const amrex::IntVect& bndry_lo,
  const amrex::IntVect& bndry_hi,
  const amrex::Box& domain)
{
  amrex::Box testBox(in_box);
  testBox.grow(1);
  if (domain.contains(testBox)) {
    return false;
  }
  bool at_bndry = false;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    amrex::Box boxlo(in_box);
    boxlo.growLo(dir, 1);
    amrex::Box boxhi(in_box);
    boxhi.growHi(dir, 1);
    if (
      (!domain.contains(boxlo) && bndry_lo[dir] != 0) ||
      (!domain.contains(boxhi) && bndry_hi[dir] != 0)) {
      at_bndry = true;
    }
  }
  return at_bndry;
}

/****************************************************************
 Functions computing wall film source terms
 ***************************************************************/
#if 0
AMREX_GPU_DEVICE AMREX_INLINE void
calculateWallFilmSource(
  const amrex::Real flow_dt,
  GasPhaseVals& gpv,
  const SprayComps& SPI,
  const SprayData& fdat,
  SprayParticleContainer::ParticleType& p,
  const amrex::Real T_wall,
  const amrex::Real /*face_area*/,
  const amrex::Real diff_cent,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* trans_parm)
{
  auto eos = pele::physics::PhysicsType::eos();
  SprayUnits SPU;
  const amrex::Real C_eps = 1.E-15;
  const amrex::Real ht_tol = 2.E-6;
  const amrex::Real part_dt = 0.5 * flow_dt;
  bool get_xi = false;
  bool get_Ddiag = true;
  bool get_lambda = true;
  bool get_mu = true;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_film;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_vapor;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> L_fuel;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> mi_dot;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cBoilT;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Psat;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_skin;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> h_skin;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Ddiag;
  for (int sp = 0; sp < NUM_SPECIES; ++sp) {
    Y_skin[sp] = 0.;
  }
  // Get particle variables
  // Note: the diameter is set assuming a sphere
  // of the same volume as the film
  // and temperature is set to a negative value
  amrex::Real T_film = -p.rdata(SPI.pstateT);
  amrex::Real vol = p.rdata(SPI.pstateVol);
  amrex::Real ht_film = p.rdata(SPI.pstateHt);
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_film[spf] = p.rdata(SPI.pstateY + spf);
  }
  calcBoilT(fdat, gpv, cBoilT.data());
  // TODO: Add model for when T > T_boil
  amrex::Real T_i = gpv.T_fluid;
  amrex::Real area_film = vol / ht_film;
  amrex::Real dia_film = std::sqrt(4. * area_film / M_PI);
  eos.T2Hi(T_film, h_skin.data());
  for (int n = 0; n < NUM_SPECIES; ++n) {
    h_skin[n] *= SPU.eng_conv;
  }
  calcVaporY(
    fdat, gpv, T_film, C_eps, Y_film.data(), h_skin.data(), cBoilT.data(),
    Y_vapor.data(), Psat.data(), L_fuel.data());
  amrex::Real T_vapor = 0.5 * (T_film + T_i);
  amrex::Real sumYSkin = 0.;
  amrex::Real sumYFuel = 0.;
  amrex::Real rho_film = 0.;
  amrex::Real lambda_film = 0.;
  amrex::Real cp_film = 0.;
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
  const amrex::Real restYSkin = 1. - sumYSkin;
  const amrex::Real restYfluid = 1. - sumYFuel;
  amrex::Real renorm = restYSkin / restYfluid;
  amrex::Real mw_vap = 0.;
  for (int sp = 0; sp < NUM_SPECIES; ++sp) {
    if (Y_skin[sp] == 0.) {
      Y_skin[sp] = gpv.Y_fluid[sp] * renorm;
    }
    mw_vap += Y_skin[sp] * gpv.invmw[sp];
  }
  mw_vap = 1. / mw_vap;
  amrex::Real pmass = vol * rho_film;
  amrex::Real lambda_skin = 0.;
  amrex::Real mu_skin = 0.;
  amrex::Real xi_skin = 0.;
  auto trans = pele::physics::PhysicsType::transport();
  trans.transport(
    get_xi, get_mu, get_lambda, get_Ddiag, T_vapor, gpv.rho_fluid,
    Y_skin.data(), Ddiag.data(), mu_skin, xi_skin, lambda_skin, trans_parm);
  mu_skin *= SPU.mu_conv;
  lambda_skin *= SPU.lambda_conv;
  amrex::Real dy_i = diff_cent - ht_film; // Distance from film surface to cell center
  // Determine the mass evaporation values
  amrex::Real m_dot = 0.;
  amrex::Real qvap = 0.;
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
  amrex::Real fs1 = lambda_skin * ht_film;
  amrex::Real fs2 = lambda_film * dy_i;
  amrex::Real T_s = (fs1 * T_i + fs2 * T_wall - ht_film * dy_i * qvap) / (fs1 + fs2);
  amrex::Real qconv = lambda_skin * (T_i - T_s) / dy_i;
  gpv.fluid_eng_src += qconv;
  // Update film height
  amrex::Real new_mass = pmass + part_dt * m_dot;
  amrex::Real new_rho = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    amrex::Real newY = (Y_film[spf] * pmass + mi_dot[spf] * part_dt) / new_mass;
    p.rdata(SPI.pstateY + spf) = newY;
    new_rho += newY / fdat.rho[spf];
  }
  new_rho = 1. / new_rho;
  amrex::Real new_vol = new_mass / new_rho;
  p.rdata(SPI.pstateVol) = new_vol;
  // Adjust the height by assuming radius remains unchanged
  amrex::Real new_ht = 4. * new_vol / (M_PI * dia_film * dia_film);
  p.rdata(SPI.pstateHt) = new_ht;
  p.rdata(SPI.pstateT) = -0.5 * (T_s + T_wall);
  if (new_ht < ht_tol) {
    p.id() = -1;
  }
}
#endif
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
splash_criteria(
  const amrex::Real Kv, const amrex::Real Tstar, const amrex::Real alpha)
{
  splash_type crit = splash_type::thermal_breakup;
  amrex::Real Kcrit = 20. + 2. * alpha / M_PI * 20.;
  if (Tstar < 1.1) {
    Kcrit = 130.;
    if (Tstar < 1.) {
      Kcrit = 54. + 76. * std::exp(13. * (Tstar - 1.));
    }
    if (Kv < Kcrit) {
      crit = splash_type::deposit;
    } else {
      crit = splash_type::splash;
    }
  }
  if (Kv < Kcrit) {
    crit = splash_type::rebound;
  }
  return crit;
}

// Find tangents along surface
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
void
find_tangents(
#if AMREX_SPACEDIM == 3
  const amrex::RealVect& testvec,
  amrex::RealVect& tanPsi,
#else
  const amrex::RealVect& /*testvec*/,
  amrex::RealVect& /*tanPsi*/,
#endif
  const amrex::RealVect& norm,
  amrex::RealVect& tanBeta)
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
  const amrex::Real& alpha,
  const amrex::Real& Tstar,
  const amrex::Real& We,
  const bool dry,
  SprayRefl& SPRF)
{
  // Inclination angle in degrees
  const amrex::Real alphad = alpha * 180. / M_PI;
  amrex::Real omega = 0.;
  if (alphad <= 80. && AMREX_SPACEDIM == 3) {
    omega = std::sqrt(
      (1. + 8.872 * std::cos(1.152 * alpha)) / (1. - std::cos(alpha)));
  }
  SPRF.omega = omega;
  SPRF.expomega = 1. - std::exp(-omega);
  amrex::Real beta_mean;
  if (dry) {
    beta_mean = 9.3 + 0.22 * alphad;
  } else if (Tstar > 1.1) {
    beta_mean = 0.225 * alphad * std::exp(std::pow(0.017 * alphad - 0.937, 2));
  } else {
    beta_mean = 0.96 * alphad * std::exp(-4.5E-3 * We);
  }
  amrex::Real stdev = 4.; // 4 degrees
  // Now convert mean and stdev to log terms
  amrex::Real mean2 = beta_mean * beta_mean;
  amrex::Real st2 = stdev * stdev;
  amrex::Real term1 = std::log(beta_mean);
  amrex::Real term2 = std::log(mean2 + st2);
  SPRF.beta_mean = 2. * term1 - 0.5 * term2;
  SPRF.beta_stdv = std::sqrt(amrex::max(-2. * term1 + term2, 0.));
}

AMREX_GPU_HOST_DEVICE
AMREX_INLINE
bool
check_wall(
  const amrex::IntVect bflags,
#ifdef AMREX_USE_EB
  const amrex::IntVect ijkc,
  const amrex::IntVect ijkc_prev,
  const bool use_EB,
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<amrex::Real const> const& bcent,
  amrex::Array4<amrex::Real const> const& bnorm,
  amrex::IntVect& bloc,
#endif
  amrex::RealVect& normal,
  amrex::RealVect& bcentv)
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
      normal = {
        AMREX_D_DECL(-bnorm(ijkc, 0), -bnorm(ijkc, 1), -bnorm(ijkc, 2))};
      bcentv = {AMREX_D_DECL(bcent(ijkc, 0), bcent(ijkc, 1), bcent(ijkc, 2))};
      wall_check = true;
      // If the particle has moved into a covered cell
    } else if (flags(ijkc).isCovered()) {
      // If the cell the particle came from is a cut-cell
      if (flags(ijkc_prev).isSingleValued()) {
        // Use the normal and boundary centroid from the previous cell
        bloc = ijkc_prev;
        normal = {AMREX_D_DECL(
          -bnorm(ijkc_prev, 0), -bnorm(ijkc_prev, 1), -bnorm(ijkc_prev, 2))};
        bcentv = {AMREX_D_DECL(
          bcent(ijkc_prev, 0), bcent(ijkc_prev, 1), bcent(ijkc_prev, 2))};
        // otherwise an error has occurred
      } else {
        amrex::Abort("WallFunctions::check_wall: Particle is outside EB");
      }
      wall_check = true;
    }
  }
#endif
  return wall_check;
}

AMREX_GPU_DEVICE
AMREX_INLINE
void
droplet_splashing(
  splash_type& splash_flag,
  SprayParticleContainer::ParticleType& p,
  const SprayComps& SPI,
  const SprayData& fdat,
  const amrex::RealVect pvel,
  const amrex::Real par_dot_EB,
  const amrex::Real Nw_Vp,
  const amrex::Real T_wall,
  const amrex::RealVect normal,
  SprayRefl& SPRF,
  const bool dry_wall)
{
  amrex::Real mu_part = 0.;
  amrex::Real rho_part = 0.;
  // TODO: Determine correct method for handling multi-component liquids
  amrex::Real Tstar = 0.; // Average sum Y_i T_wall/T_boil,i
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int ysp = SPI.pstateY + spf;
    mu_part += p.rdata(ysp) * fdat.mu[spf];
    rho_part += p.rdata(ysp) / fdat.rho[spf];
    Tstar += T_wall * p.rdata(ysp) / fdat.boilT[spf];
  }
  rho_part = 1. / rho_part;
  const amrex::Real dia_part = p.rdata(SPI.pstateDia);
  const amrex::Real d3 = std::pow(dia_part, 3);
  const amrex::Real pmass = M_PI / 6. * rho_part * d3;
  // Weber number
  const amrex::Real We = rho_part * dia_part * Nw_Vp * Nw_Vp / fdat.sigma;
  const amrex::Real Re_L = std::abs(Nw_Vp) * dia_part * rho_part / mu_part;
  const amrex::Real sqReyn = std::sqrt(Re_L);
  const amrex::Real Kv = std::sqrt(We) * std::pow(Re_L, 0.25);
  amrex::Real umag = pvel.vectorLength();
  // Inclination of on-coming particle to surface in radians
  const amrex::Real alpha = std::asin(std::abs(Nw_Vp) / umag);
  splash_flag = splash_criteria(Kv, Tstar, alpha);
  // Boundary layer thickness of particle estimate provided by
  // Pasandideh-Fard et al. 1996
  const amrex::Real delta = 2. * dia_part / sqReyn;
  // splash_flag == 0, then particle rebounds, no breakup, no deposit
  if (splash_flag != splash_type::rebound) {
    // If splash_flag == 1, particle is completely deposited onto the wall
    if (splash_flag == splash_type::deposit) {
      // TODO: Determine better way to estimate wall film temperature,
      // diameter, etc
      amrex::Real depot_vol = pmass / rho_part;
      // Assume deposited film is a dome shape of height delta, determine
      // the diameter
      amrex::Real depot_dia =
        2. *
        std::sqrt((6. * depot_vol / M_PI - std::pow(delta, 3)) / (3. * delta));
      // Now estimate the height to be a cylinder of diameter depot_dia
      amrex::Real depot_height =
        4. * depot_vol / (M_PI * depot_dia * depot_dia);
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        p.rdata(SPI.pstateVel + dir) = 0.;
      }
      p.rdata(SPI.pstateVol) = depot_vol;
      p.rdata(SPI.pstateHt) = depot_height;
    } else {
      // Determine diameter of secondary droplets
      amrex::Real expon = 3.6 * (alpha / M_PI) * (alpha / M_PI);
      if (dry_wall) {
        SPRF.dia_refl = dia_part * 3.3 * std::exp(expon) * std::pow(We, -0.65);
      } else {
        SPRF.dia_refl = dia_part * 2.2 * std::exp(expon) * std::pow(We, -0.36);
      }
      // Compute angles of reflecting droplets
      compute_angles(alpha, Tstar, We, dry_wall, SPRF);
      // If thermal breakup, all mass is converted to secondary droplets
      amrex::Real splash_mass = pmass;
      // If splash, some mass rebounds and some deposits
      if (splash_flag == splash_type::splash) {
        // Determine the fraction of the droplet mass
        // that forms secondary drops from Kuhnke 2004
        // TODO: This is incomplete by always assuming a dry wall
        // TODO: Random(engine) degrades performance too much, must add
        // random values here somehow
        const amrex::Real B = 0.2 + 0.6; //* amrex::Random(engine);
        splash_mass *= amrex::min(1., (Tstar - 0.8) / 0.3 * (1. - B) + B);
        // Mass deposited into wall film
        amrex::Real depot_mass = pmass - splash_mass;
        amrex::Real depot_vol = depot_mass / rho_part;
        // Assume deposited film is a dome shape of height delta,
        // determine the diameter
        amrex::Real depot_dia =
          2. * std::sqrt(
                 (6. * depot_vol / M_PI - std::pow(delta, 3)) / (3. * delta));
        // Now estimate the height to be a cylinder of diameter depot_dia
        amrex::Real depot_height =
          4. * depot_vol / (M_PI * depot_dia * depot_dia);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          p.rdata(SPI.pstateVel + dir) = 0.;
        }
        p.rdata(SPI.pstateVol) = depot_vol;
        p.rdata(SPI.pstateHt) = depot_height;
      } else {
        p.id() = -1;
      }
      amrex::Real mass_refl = M_PI / 6. * rho_part * std::pow(SPRF.dia_refl, 3);
      amrex::Real Ns = splash_mass / mass_refl;
      SPRF.Ns_refl = int(Ns);
      // TODO: Not sure what this variable is but it ranges from 1.45-2
      amrex::Real nu32 = 1.45;
      amrex::Real We_out =
        SPRF.dia_refl / dia_part *
          (We * (1. - 0.85 * std::pow(std::sin(alpha * M_PI / 180.), 2)) +
           12.) -
        12. / nu32;
      SPRF.Unorm = std::sqrt(fdat.sigma * We_out / (rho_part * SPRF.dia_refl));
      // This is the dt to move the secondary droplets
      SPRF.dt_pp = par_dot_EB / Nw_Vp;
      // Test vector for finding tangents
      amrex::RealVect testvec = -pvel;
      // Find two tangent vectors by taking the cross product with the
      // velocity vector
      find_tangents(testvec, SPRF.tanPsi, normal, SPRF.tanBeta);
    }
  } // if (splash_flag != splash_type::rebound)
}

AMREX_GPU_DEVICE
AMREX_INLINE
void
impose_wall(
  SprayParticleContainer::ParticleType& p,
  const SprayComps& SPI,
  const SprayData& fdat,
  const amrex::RealVect& dx,
  const amrex::RealVect& plo,
  const amrex::Real T_wall,
  const amrex::IntVect bloc,
  const amrex::RealVect normal,
  const amrex::RealVect bcentv,
  SprayRefl& SPRF,
  bool isActive,
  const bool dry_wall)
//  amrex::RandomEngine const& engine) // Degrades performance significantly
{
  const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();
  splash_type splash_flag;
  amrex::RealVect pvel = {AMREX_D_DECL(
    p.rdata(SPI.pstateVel), p.rdata(SPI.pstateVel + 1),
    p.rdata(SPI.pstateVel + 2))};
  // Projection of vector pointing from EB centroid to particle onto EB normal
  const amrex::Real par_dot_EB = AMREX_D_TERM(
    (p.pos(0) - (bloc[0] + 0.5 + bcentv[0]) * dx[0] - plo[0]) * normal[0],
    +(p.pos(1) - (bloc[1] + 0.5 + bcentv[1]) * dx[1] - plo[1]) * normal[1],
    +(p.pos(2) - (bloc[2] + 0.5 + bcentv[2]) * dx[2] - plo[2]) * normal[2]);
  // Check if particle is on the wrong side of the EB
  if (par_dot_EB < tolerance) {
    splash_flag = splash_type::rebound;
    // Plane ref point
    amrex::Real Nw_Vp = 0.;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      Nw_Vp += normal[dir] * pvel[dir];
    }
    // Parcel normal velocity
    amrex::RealVect Vpn = {
      AMREX_D_DECL(Nw_Vp * normal[0], Nw_Vp * normal[1], Nw_Vp * normal[2])};
    // Parcel tangential velocity
    amrex::RealVect Vpt = pvel - Vpn;
    // If using a splash model
    // TODO: Implement splash model
    if (fdat.sigma > 0. && isActive) {
      droplet_splashing(
        splash_flag, p, SPI, fdat, pvel, par_dot_EB, Nw_Vp, T_wall, normal,
        SPRF, dry_wall);
    } // if (sigma > 0.)
    // Reflect particle
    if (fdat.sigma < 0. || splash_flag == splash_type::rebound) {
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
}

void
create_splash_droplet(
  SprayParticleContainer::ParticleType& p, SprayComps SPI, SprayRefl SPRF)
{
  amrex::Real rand1 = amrex::Random();
  amrex::Real mean = SPRF.beta_mean;
  amrex::Real stdev = SPRF.beta_stdv;
  amrex::Real beta = amrex::RandomNormal(mean, stdev);
  beta = std::exp(beta) * M_PI / 180.;
  amrex::Real omega = SPRF.omega;
  amrex::Real expb = SPRF.expomega;
  // If incoming angle is greater than 80 degrees, azimuthal
  // angle of outgoing drop can be between 0 to 2*pi
  amrex::Real psi = AMREX_D_PICK(0., 0., rand1 * 2. * M_PI);
  // Otherwise, follow we modify the recommendation by Naber and Reitz 1988,
  // so the azimuthal angle distribution favors the pre-splash drop path
  // as the inclination angle decreases
  if (omega > 0.) {
    amrex::Real rand2 = std::copysign(1., 0.5 - amrex::Random());
    psi = -rand2 / omega * std::log(1. - rand1 * expb) * M_PI;
  }
  amrex::Real costhetad = std::cos(beta);
  amrex::Real sinthetad = std::sin(beta);
  AMREX_D_TERM(amrex::Real un = SPRF.Unorm * sinthetad;
               , amrex::Real utBeta = SPRF.Unorm * costhetad * std::cos(psi);
               , amrex::Real utPsi = SPRF.Unorm * costhetad * std::sin(psi););
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    amrex::Real pvel = AMREX_D_TERM(
      un * SPRF.norm[dir], +utBeta * SPRF.tanBeta[dir],
      +utPsi * SPRF.tanPsi[dir]);
    amrex::Gpu::Atomic::Add(&p.pos(dir), SPRF.dt_pp * pvel);
    amrex::Gpu::Atomic::Add(&p.rdata(SPI.pstateVel + dir), pvel);
  }
}

#endif // _WALLFUNCTIONS_H_
