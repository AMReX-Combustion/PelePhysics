
#ifndef WALLFUNCTIONS_H_
#define WALLFUNCTIONS_H_

#include "Drag.H"
#include "SprayFuelData.H"

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

/**
Determines if particle must be reflected off BC or EB
@param[in] pos Position of the particle after advecting, normalized by dx
@param[in] bflags Flags if particle is outside reflective boundaries
@param[in] ijkc Grid cell index containing particle
@param[in] ijkc_prev Previous grid cell index containing particle
@param[in] vel_part Velocity of the particle
@param[in] use_EB Flag if EB is used in the current box
@param[in] flags Array of flags denoting if a cell has EB in it
@param[in] bcent Array of EB centroids for each cell
@param[in] bnorm Array of EB normal vectors for each cell
@param[in] vfrac Array of EB volume fractions for each cell
@param[in] min_eb_vfrac Minimum volume fraction to allow particles in
@param[out] par_dot Magnitude of vector from wall plane to particle, negative
means outside the domain
@param[out] normal Boundary normal vector
@return Flag if particle must be reflected
*/
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
bool
check_wall(
  const amrex::RealVect pos,
  const amrex::IntVect bflags,
  const amrex::IntVect ijkc,
#ifdef AMREX_USE_EB
  const amrex::IntVect ijkc_prev,
  const amrex::RealVect vel_part,
  const bool use_EB,
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<amrex::Real const> const& bcent,
  amrex::Array4<amrex::Real const> const& bnorm,
  amrex::Array4<amrex::Real const> const& vfrac,
  const amrex::Real min_eb_vfrac,
#endif
  amrex::Real& par_dot,
  amrex::RealVect& normal)
{
  bool wall_check = false; // Check if particle is reflected
  amrex::RealVect bcentv;  // Absolute location of point on boundary face
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    // -1 - Outside lower reflective boundary
    // 1 - Outside upper reflective boundary
    // Only accounts for 1  Cartesian direction
    if ((bflags[dir] == -1 || bflags[dir] == 1) && !wall_check) {
      normal[dir] = -amrex::Real(bflags[dir]);
      bcentv[dir] = amrex::Real(ijkc[dir]) + 0.5 * amrex::Real(1 - bflags[dir]);
      wall_check = true;
    } else {
      normal[dir] = 0.;
      bcentv[dir] = 0.;
    }
  }
#ifdef AMREX_USE_EB
  if (use_EB && !wall_check) {
    // If particle has moved into a cut-cell
    if (flags(ijkc).isSingleValued()) {
      wall_check = true;
      // Use the normal and boundary centroid from the new cell
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        normal[dir] = -bnorm(ijkc, dir);
        bcentv[dir] = amrex::Real(ijkc[dir]) + 0.5 + bcent(ijkc, dir);
      }
      // If volume fraction is too small, treat cell face as the wall to prevent
      // particle from residing in cell
      if (vfrac(ijkc) < min_eb_vfrac) {
        // Direction of face to treat as wall
        int fdir = normal.maxDir(true);
        amrex::Real face_sgn = (bcent(ijkc, fdir) < 0.) ? 1. : -1.;
        normal = amrex::RealVect::TheZeroVector();
        normal[fdir] = -face_sgn;
        bcentv[fdir] = amrex::Real(ijkc[fdir]) + 0.5 * (1. - face_sgn);
      }
      // If the particle has moved into a covered cell
    } else if (flags(ijkc).isCovered()) {
      wall_check = true;
      // If the particle came from a cut-cell, check if it is behind that cells
      // EB
      amrex::Real test_par = 2.;
      if (flags(ijkc_prev).isSingleValued()) {
        test_par = 0.;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          normal[dir] = -bnorm(ijkc_prev, dir);
          bcentv[dir] =
            amrex::Real(ijkc_prev[dir]) + 0.5 + bcent(ijkc_prev, dir);
          test_par += (pos[dir] - bcentv[dir]) * normal[dir];
        }
      }
      // If it didn't come from a cut-cell or isn't behind the EB of the
      // previous cell, more tests are required
      if (test_par > 0.) {
        amrex::Real max_diff = 0.;
        amrex::RealVect cur_bcent;
        amrex::RealVect cur_norm;
        // Determine the distance from the particle to the nearest wall
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          if (ijkc_prev[dir] != ijkc[dir]) {
            amrex::Real cur_diff = 0.;
            amrex::IntVect tsti = ijkc_prev;
            tsti[dir] = ijkc[dir];
            // If droplet passes through adjacent cell with EB, find distance to
            // EB
            if (flags(tsti).isSingleValued()) {
              amrex::Real diff1 = 0.;
              amrex::Real diff2 = 0.;
              for (int okdir = 0; okdir < AMREX_SPACEDIM; ++okdir) {
                cur_norm[okdir] = -bnorm(tsti, okdir);
                cur_bcent[okdir] =
                  amrex::Real(tsti[okdir]) + 0.5 + bcent(tsti, okdir);
                diff1 += (pos[okdir] - cur_bcent[okdir]) * cur_norm[okdir];
                diff2 += vel_part[okdir] * cur_norm[okdir];
              }
              if (diff2 != 0.) {
                cur_diff = diff1 / diff2;
              }
              // If droplet passes through another covered cell, find distance
              // to normal cell face
            } else if (flags(tsti).isCovered()) {
              for (int okdir = 0; okdir < AMREX_SPACEDIM; ++okdir) {
                cur_bcent[okdir] = amrex::Real(tsti[okdir]) + 0.5;
                cur_norm[okdir] = 0.;
              }
              amrex::Real sgn = (ijkc_prev[dir] > ijkc[dir]) ? 1. : -1.;
              cur_bcent[dir] += sgn * 0.5;
              cur_norm[dir] = sgn;
              cur_diff = (pos[dir] - cur_bcent[dir]) / vel_part[dir];
            }
            if (cur_diff > max_diff) {
              max_diff = cur_diff;
              normal = cur_norm;
              bcentv = cur_bcent;
            }
          }
        } // for (int dir = 0;...
      }   // if (test_par > 0.)
    }
  }
#endif
  bool refl_check = false;
  if (wall_check) {
    const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();
    // Projection of vector pointing from EB centroid to particle onto EB normal
    par_dot = AMREX_D_TERM(
      (pos[0] - bcentv[0]) * normal[0], +(pos[1] - bcentv[1]) * normal[1],
      +(pos[2] - bcentv[2]) * normal[2]);
    if (par_dot < tolerance) {
      refl_check = true;
    }
  }
  return refl_check;
}

/*
Checks and reflects particles off BC or EB. If particle reflects too many times,
it places the particle near the wall sets the normal velocity component to zero.
*/
AMREX_GPU_DEVICE
AMREX_INLINE
void
impose_wall(
  int pid,
  SprayParticleContainer::ParticleType& p,
  const SprayComps& SPI,
  const amrex::RealVect& dx,
  const amrex::RealVect& plo,
  const amrex::IntVect bflags,
  const bool use_EB,
#ifdef AMREX_USE_EB
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<amrex::Real const> const& bcent,
  amrex::Array4<amrex::Real const> const& bnorm,
  amrex::Array4<amrex::Real const> const& vfrac,
  const amrex::Real min_eb_vfrac,
#endif
  amrex::IntVect& ijkc,
  amrex::IntVect& ijkc_prev,
  int* Nrefl_ptr,
  amrex::Real* norm_ptr,
  amrex::Real* reflvel_ptr,
  amrex::Real* reflloc_ptr)
{
  amrex::Real par_dot = 2.;
  amrex::RealVect normal;
  amrex::RealVect vel_part(AMREX_D_DECL(
    p.rdata(SPI.pstateVel), p.rdata(SPI.pstateVel + 1),
    p.rdata(SPI.pstateVel + 2)));
  // Normalize pos
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    p.pos(dir) = (p.pos(dir) - plo[dir]) / dx[dir];
  }
  // Check if particle is behind a wall
  bool do_refl = check_wall(
    p.pos(), bflags, ijkc,
#ifdef AMREX_USE_EB
    ijkc_prev, vel_part, use_EB, flags, bcent, bnorm, vfrac, min_eb_vfrac,
#endif
    par_dot, normal);
  // Keep track of number of times particle bounces
  int bounce_counter = 0;
  // Maximum allowable bounces
  int max_bounce = 5;
  while (do_refl) {
    // Plane ref point
    amrex::Real Nw_Vp = normal.dotProduct(vel_part);
    if (!use_EB) {
      ijkc_prev = ijkc;
    }
#ifdef AMREX_USE_EB
    else {
      if (!flags(ijkc).isCovered()) {
        ijkc_prev = ijkc;
      }
    }
#endif
    // Reflect particle position and velocity around the boundary
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      p.pos(dir) -= 2. * par_dot * normal[dir];
      p.rdata(SPI.pstateVel + dir) -= 2. * Nw_Vp * normal[dir];
      vel_part[dir] = p.rdata(SPI.pstateVel + dir);
      ijkc[dir] = static_cast<int>(amrex::Math::floor(p.pos(dir)));
    }
    do_refl = check_wall(
      p.pos(), bflags, ijkc,
#ifdef AMREX_USE_EB
      ijkc_prev, vel_part, use_EB, flags, bcent, bnorm, vfrac, min_eb_vfrac,
#endif
      par_dot, normal);
    bounce_counter++;
    // If particle exceeds max bounces, place near wall and set normal velocity
    // to zero
    if (bounce_counter > max_bounce && do_refl) {
      Nw_Vp = normal.dotProduct(vel_part);
      Nrefl_ptr[pid] = 10;
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        norm_ptr[3*pid + dir] = normal[dir];
        reflvel_ptr[3*pid+dir] = p.rdata(SPI.pstateVel+dir);
        p.pos(dir) +=
          (0.5 * p.rdata(SPI.pstateDia) / dx[dir] - par_dot) * normal[dir];
        reflloc_ptr[3*pid+dir] = p.pos(dir);
        p.rdata(SPI.pstateVel + dir) -= Nw_Vp * normal[dir];
        ijkc[dir] = static_cast<int>(amrex::Math::floor(p.pos(dir)));
      }
      do_refl = false;
    }
  }
  // Undo normalization
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    p.pos(dir) = p.pos(dir) * dx[dir] + plo[dir];
  }
}

/****************************************************************
 Functions computing wall film source terms
 WIP
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
  const amrex::Real Re_L = amrex::Math::abs(Nw_Vp) * dia_part * rho_part / mu_part;
  const amrex::Real sqReyn = std::sqrt(Re_L);
  const amrex::Real Kv = std::sqrt(We) * std::pow(Re_L, 0.25);
  amrex::Real umag = pvel.vectorLength();
  // Inclination of on-coming particle to surface in radians
  const amrex::Real alpha = std::asin(amrex::Math::abs(Nw_Vp) / umag);
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
        amrex::Real depot_dia = std::cbrt(6. / M_PI * depot_mass / rho_part);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          p.rdata(SPI.pstateVel + dir) = 0.;
        }
        p.rdata(SPI.pstateDia) = depot_dia;
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
#endif // #if 0
#endif // _WALLFUNCTIONS_H_
