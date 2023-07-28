
#ifndef WALLFUNCTIONS_H
#define WALLFUNCTIONS_H

#include "AhamedSplash.H"

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

// This is the same as the check_bounds function but pos and phi have been
// normalized by plo and dx
AMREX_GPU_DEVICE AMREX_INLINE bool
check_bounds_norm(
  const amrex::RealVect pos,
  const amrex::RealVect phi,
  const amrex::IntVect bndry_lo,
  const amrex::IntVect bndry_hi,
  amrex::IntVect& bflags)
{
  for (int hilo = 0; hilo < 2; ++hilo) {
    int fact = (hilo == 0) ? 1 : -1;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      amrex::Real dom_loc = (hilo == 0) ? 0. : phi[dir];
      int bflag = (hilo == 0) ? bndry_lo[dir] : bndry_hi[dir];
      if (bflag != 0) { // Not periodic
        amrex::Real diff =
          static_cast<amrex::Real>(fact) * (pos[dir] - dom_loc);
        if (diff < 0.) {
          // Outside domain
          if (bflag == 1) { // Outside reflective BC
            bflags[dir] = -fact;
          } else { // Outside non-reflective BC
            return true;
          }
          // Adjacent to non-periodic boundary
        } else if (diff < 0.5) {
          bflags[dir] = -2 * fact;
        }
      }
    }
  }
  return false;
}

/**
Determines if particle must be reflected off BC or EB
@param[in] pos Position of the particle after advecting, normalized by dx
@param[in] bflags Flags if particle is outside reflective boundaries
@param[in] dx Grid spacing
@param[in] ijkc Grid cell index containing particle
@param[in] ijkc_prev Previous grid cell index containing particle
@param[in] vel_part Velocity of the particle
@param[in] use_EB Flag if EB is used in the current box
@param[in] flags Array of flags denoting if a cell has EB in it
@param[in] bcent Array of EB centroids for each cell
@param[in] bnorm Array of EB normal vectors for each cell
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
  const amrex::RealVect dx,
  const amrex::IntVect ijkc,
#ifdef AMREX_USE_EB
  const amrex::IntVect& ijkc_prev,
  const amrex::RealVect& vel_part,
  const bool use_EB,
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<amrex::Real const> const& bcent,
  amrex::Array4<amrex::Real const> const& bnorm,
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
      normal[dir] = -static_cast<amrex::Real>(bflags[dir]);
      bcentv[dir] = static_cast<amrex::Real>(ijkc[dir]) +
                    0.5 * static_cast<amrex::Real>(1 - bflags[dir]);
      wall_check = true;
    } else {
      normal[dir] = 0.;
      bcentv[dir] = 0.;
    }
  }
#ifdef AMREX_USE_EB
  if (use_EB && !wall_check) {
    if (flags(ijkc).isSingleValued()) {
      wall_check = true;
      normal = {
        AMREX_D_DECL(-bnorm(ijkc, 0), -bnorm(ijkc, 1), -bnorm(ijkc, 2))};
      bcentv = {AMREX_D_DECL(
        ijkc[0] + 0.5 + bcent(ijkc, 0), ijkc[1] + 0.5 + bcent(ijkc, 1),
        ijkc[2] + 0.5 + bcent(ijkc, 2))};
    } else if (flags(ijkc).isCovered()) {
      wall_check = true;
      if (flags(ijkc_prev).isSingleValued()) {
        normal = {AMREX_D_DECL(
          -bnorm(ijkc_prev, 0), -bnorm(ijkc_prev, 1), -bnorm(ijkc_prev, 2))};
        bcentv = {AMREX_D_DECL(
          ijkc_prev[0] + 0.5 + bcent(ijkc_prev, 0),
          ijkc_prev[1] + 0.5 + bcent(ijkc_prev, 1),
          ijkc_prev[2] + 0.5 + bcent(ijkc_prev, 2))};
      } else {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          if (ijkc[dir] != ijkc_prev[dir]) {
            amrex::Real face_sgn = (ijkc[dir] > ijkc_prev[dir]) ? -1. : 1.;
            normal[dir] = face_sgn;
            bcentv[dir] = static_cast<amrex::Real>(ijkc[dir]) + 0.5 * face_sgn;
          } else {
            normal[dir] = 0.;
            bcentv[dir] = static_cast<amrex::Real>(ijkc[dir]);
          }
        }
        normal /= normal.vectorLength();
      }
    }
  }
#endif
  bool refl_check = false;
  if (wall_check) {
    const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();
    // Projection of vector pointing from EB centroid to particle onto EB normal
    par_dot = AMREX_D_TERM(
      (pos[0] - bcentv[0]) * dx[0] * normal[0], +(pos[1] - bcentv[1]) * dx[1] * normal[1],
      +(pos[2] - bcentv[2]) * dx[2] * normal[2]);
    if (par_dot < tolerance) {
      refl_check = true;
    }
  }
  return refl_check;
}

/*
Imposes BC or EB surfaces. Splash model is applied if it is turned on.
If the particle reflects multiple times, splash model is not applied to avoid
creating particles outside of boundary
*/
AMREX_GPU_DEVICE
AMREX_INLINE
void
impose_wall(
  bool do_splash,
  int pid,
  SprayParticleContainer::ParticleType& p,
  const amrex::RealVect& dx,
  const amrex::RealVect& plo,
  const amrex::RealVect& phi,
  amrex::IntVect& bflags,
  const amrex::Real* cBoilT,
  const amrex::Real p_fluid,
  const bool use_EB,
#ifdef AMREX_USE_EB
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<amrex::Real const> const& bcent,
  amrex::Array4<amrex::Real const> const& bnorm,
  amrex::Array4<amrex::Real const> const& vfrac,
  const amrex::Real min_eb_vfrac,
  const amrex::IntVect ijkc_prev,
#endif
  const amrex::IntVect ijkc)
{
  amrex::Real par_dot = 2.;
  amrex::RealVect normal;
  amrex::RealVect vel_part(AMREX_D_DECL(
    p.rdata(SprayComps::pstateVel), p.rdata(SprayComps::pstateVel + 1),
    p.rdata(SprayComps::pstateVel + 2)));
  amrex::RealVect normphi;
  amrex::RealVect normpos;
  // Normalize pos
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    normpos[dir] = (p.pos(dir) - plo[dir]) / dx[dir];
    normphi[dir] = (phi[dir] - plo[dir]) / dx[dir];
  }
  // Check if particle is behind a wall
  bool do_refl = check_wall(
    normpos, bflags, dx, ijkc,
#ifdef AMREX_USE_EB
    ijkc_prev, vel_part, use_EB, flags, bcent, bnorm,
#endif
    par_dot, normal);
  if (do_refl) {
    amrex::Real Nw_Vp = normal.dotProduct(vel_part);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      p.rdata(SprayComps::pstateVel + dir) -= 2. * Nw_Vp * normal[dir];
      p.pos(dir) -= 2. * par_dot * normal[dir];
    }
  }
}
#endif // WALLFUNCTIONS_H
