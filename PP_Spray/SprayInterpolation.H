#ifndef _SPRAYINTERPOLATION_H_
#define _SPRAYINTERPOLATION_H_

#include <AMReX_RealVect.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
InterpolateFaceVelocity(
  const amrex::RealVect& len,
  const amrex::IntVect& dom_lo,
  const amrex::IntVect& dom_hi,
  amrex::GpuArray<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> const& umac,
  amrex::RealVect& vel_fluid)
{
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    amrex::RealVect face_len(AMREX_D_DECL(
      len[0] + (dir == 0) * 0.5, len[1] + (dir == 1) * 0.5,
      len[2] + (dir == 2) * 0.5));
    // Compute distances for trilinear interpolation
    amrex::IntVect indx(AMREX_D_DECL(
      amrex::Math::floor(face_len[0]), amrex::Math::floor(face_len[1]),
      amrex::Math::floor(face_len[2])));
    amrex::RealVect delL(AMREX_D_DECL(
      face_len[0] - indx[0], face_len[1] - indx[1], face_len[2] - indx[2]));
    int iin[2] = {
      amrex::max(indx[0] - 1, dom_lo[0]), amrex::min(indx[0], dom_hi[0])};
    amrex::Real ssx[2] = {1. - delL[0], delL[0]};
    int zindx = 0;
    int yindx = 0;
    int xindx = 0;
#if AMREX_SPACEDIM > 1
    int jjn[2] = {
      amrex::max(indx[1] - 1, dom_lo[1]), amrex::min(indx[1], dom_hi[1])};
    amrex::Real ssy[2] = {1. - delL[1], delL[1]};
#if AMREX_SPACEDIM == 3
    int kkn[2] = {
      amrex::max(indx[2] - 1, dom_lo[2]), amrex::min(indx[2], dom_hi[2])};
    amrex::Real ssz[2] = {1. - delL[2], delL[2]};
    for (int kk = 0; kk < 2; ++kk) {
      zindx = kkn[kk];
#endif
      for (int jj = 0; jj < 2; ++jj) {
        yindx = jjn[jj];
#endif
        for (int ii = 0; ii < 2; ++ii) {
          xindx = iin[ii];
          vel_fluid[dir] += (umac[dir])(xindx, yindx, zindx, 0) *
                            AMREX_D_TERM(ssx[ii], *ssy[jj], *ssz[kk]);
        }
#if AMREX_SPACEDIM > 1
      }
#if AMREX_SPACEDIM == 3
    }
#endif
#endif
  }
}

// Check if outside of bounds of Cartesian grid
// Return IntVect that where the flags are per direction
// 0 - Not outside domain or at periodic boundary
// -1 - Outside lower reflective domain
// 1 - Outside upper reflective domain
// -2 - Adjacent to lower reflective domain boundary
// 2 - Adjacent to upper reflective domain boundary
// And true is returned if outside non-reflective domain
// If near a domain, ijk is shifted inward for interpolation purposes
AMREX_GPU_HOST_DEVICE AMREX_INLINE bool
check_bounds(
  const amrex::RealVect pos,
  const amrex::RealVect plo,
  const amrex::RealVect phi,
  const amrex::RealVect dx,
  const amrex::IntVect bndry_lo,
  const amrex::IntVect bndry_hi,
  amrex::IntVect& ijk,
  amrex::IntVect& dflags)
{
  for (int hilo = 0; hilo < 2; ++hilo) {
    amrex::Real fact = (hilo == 0) ? 1. : -1.;
    amrex::RealVect dom_ex = (hilo == 0) ? plo : phi;
    amrex::IntVect bndry_flag = (hilo == 0) ? bndry_lo : bndry_hi;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      amrex::Real dom_loc = dom_ex[dir];
      int bflag = bndry_flag[dir];
      if (bflag != 0) { // Not periodic
        amrex::Real diff = fact * (pos[dir] - dom_loc);
        if (diff < 0.) {
          if (bflag == 1)
            dflags[dir] = -int(fact);
          else
            return true;
          // Adjacent to non-periodic boundary, shift interpolation cells
        } else if (diff < 0.5 * dx[dir]) {
          dflags[dir] = -2 * int(fact);
          ijk[dir] += fact;
        }
      }
    }
  }
  return false;
}

// Slightly modified from MFIX code

/****************************************************************
 Functions for interpolation on non-EB mesh
 ***************************************************************/

AMREX_GPU_HOST_DEVICE AMREX_INLINE void
trilinear_interp(
  const amrex::IntVect& ijk,
  const amrex::RealVect& lx,
  amrex::IntVect* indx_array,
  amrex::Real* weights,
  const amrex::IntVect& bflags)
{
  // Note: if near a reflective boundary, ijk has been shifted in check_bounds
  AMREX_D_TERM(int i = ijk[0];, int j = ijk[1];, int k = ijk[2];);

  const amrex::RealVect sx_hi = lx - ijk;
  const amrex::RealVect sx_lo = 1. - sx_hi;
  amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>, 2> ssv;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    // Check if adjacent to upper boundary
    // or lower boundary and make it zeroth-order interpolation
    if (bflags[dir] > 1) {
      ssv[0][dir] = 0.;
      ssv[1][dir] = 1.;
    } else if (bflags[dir] < -1) {
      ssv[0][dir] = 1.;
      ssv[1][dir] = 0.;
    } else {
      ssv[0][dir] = sx_lo[dir];
      ssv[1][dir] = sx_hi[dir];
    }
  }
  int cc = 0;
#if AMREX_SPACEDIM == 3
  for (int kk(-1); kk < 1; kk++) {
#endif
#if AMREX_SPACEDIM > 1
    for (int jj(-1); jj < 1; jj++) {
#endif
      for (int ii(-1); ii < 1; ii++) {
        AMREX_D_TERM(indx_array[cc][0] = i + ii;, indx_array[cc][1] = j + jj;
                     , indx_array[cc][2] = k + kk;);
        weights[cc] =
          AMREX_D_TERM(ssv[ii + 1][0], *ssv[jj + 1][1], *ssv[kk + 1][2]);
        cc++;
      }
#if AMREX_SPACEDIM > 1
    }
#endif
#if AMREX_SPACEDIM == 3
  }
#endif
}

#ifdef AMREX_USE_EB

/****************************************************************
 Functions for the Newtons solver
 ***************************************************************/

AMREX_GPU_HOST_DEVICE AMREX_INLINE amrex::Real
f(const int dir,
  const amrex::RealVect& pos,
  const amrex::GpuArray<amrex::GpuArray<amrex::Real, 3>, 8>& nodes,
  const amrex::Real& xi,
  const amrex::Real& eta,
  const amrex::Real& zeta)
{
  amrex::Real a0 = (nodes[0][dir] - pos[dir]);
  amrex::Real a1 = (nodes[1][dir] - nodes[0][dir]);
  amrex::Real a2 = (nodes[2][dir] - nodes[0][dir]);
  amrex::Real a3 = (nodes[4][dir] - nodes[0][dir]);
  amrex::Real a4 =
    (nodes[0][dir] - nodes[1][dir] + nodes[3][dir] - nodes[2][dir]);
  amrex::Real a5 =
    (nodes[0][dir] - nodes[1][dir] - nodes[4][dir] + nodes[5][dir]);
  amrex::Real a6 =
    (nodes[0][dir] - nodes[2][dir] - nodes[4][dir] + nodes[6][dir]);
  amrex::Real a7 =
    (nodes[1][dir] - nodes[3][dir] + nodes[2][dir] + nodes[4][dir] -
     nodes[5][dir] + nodes[7][dir] - nodes[6][dir] - nodes[0][dir]);

  return a0 + a1 * xi + a2 * eta + a3 * zeta + a4 * xi * eta + a5 * xi * zeta +
         a6 * eta * zeta + a7 * xi * eta * zeta;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE amrex::Real
dfdxi(
  const int dir,
  const amrex::GpuArray<amrex::GpuArray<amrex::Real, 3>, 8>& nodes,
  const amrex::Real& xi,
  const amrex::Real& eta,
  const amrex::Real& zeta)
{
  amrex::Real a1 = (nodes[1][dir] - nodes[0][dir]);
  amrex::Real a4 =
    (nodes[0][dir] - nodes[1][dir] + nodes[3][dir] - nodes[2][dir]);
  amrex::Real a5 =
    (nodes[0][dir] - nodes[1][dir] - nodes[4][dir] + nodes[5][dir]);
  amrex::Real a7 =
    (nodes[1][dir] - nodes[3][dir] + nodes[2][dir] + nodes[4][dir] -
     nodes[5][dir] + nodes[7][dir] - nodes[6][dir] - nodes[0][dir]);

  return a1 + a4 * eta + a5 * zeta + a7 * eta * zeta;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE amrex::Real
dfdeta(
  const int dir,
  const amrex::GpuArray<amrex::GpuArray<amrex::Real, 3>, 8>& nodes,
  const amrex::Real& xi,
  const amrex::Real& eta,
  const amrex::Real& zeta)
{

  amrex::Real a2 = (nodes[2][dir] - nodes[0][dir]);
  amrex::Real a4 =
    (nodes[0][dir] - nodes[1][dir] + nodes[3][dir] - nodes[2][dir]);
  amrex::Real a6 =
    (nodes[0][dir] - nodes[2][dir] - nodes[4][dir] + nodes[6][dir]);
  amrex::Real a7 =
    (nodes[1][dir] - nodes[3][dir] + nodes[2][dir] + nodes[4][dir] -
     nodes[5][dir] + nodes[7][dir] - nodes[6][dir] - nodes[0][dir]);

  return a2 + a4 * xi + a6 * zeta + a7 * xi * zeta;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE amrex::Real
dfdzeta(
  const int dir,
  const amrex::GpuArray<amrex::GpuArray<amrex::Real, 3>, 8>& nodes,
  const amrex::Real& xi,
  const amrex::Real& eta,
  const amrex::Real& zeta)
{

  amrex::Real a3 = (nodes[4][dir] - nodes[0][dir]);
  amrex::Real a5 =
    (nodes[0][dir] - nodes[1][dir] - nodes[4][dir] + nodes[5][dir]);
  amrex::Real a6 =
    (nodes[0][dir] - nodes[2][dir] - nodes[4][dir] + nodes[6][dir]);
  amrex::Real a7 =
    (nodes[1][dir] - nodes[3][dir] + nodes[2][dir] + nodes[4][dir] -
     nodes[5][dir] + nodes[7][dir] - nodes[6][dir] - nodes[0][dir]);

  return a3 + a5 * xi + a6 * eta + a7 * xi * eta;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE void
get_interp_mapping(
  const amrex::RealVect& pos,
  const amrex::GpuArray<amrex::GpuArray<amrex::Real, 3>, 8>& nodes,
  amrex::Real& xi,
  amrex::Real& eta,
  amrex::Real& zeta)
{
  int lc(0);
  amrex::Real err(1.0);

  while (err > 1.0e-3 && lc < 10) {

    amrex::Real f0 = f(0, pos, nodes, xi, eta, zeta);
    amrex::Real f1 = f(1, pos, nodes, xi, eta, zeta);
    amrex::Real f2 = f(2, pos, nodes, xi, eta, zeta);

    amrex::Real df0dxi = dfdxi(0, nodes, xi, eta, zeta);
    amrex::Real df0deta = dfdeta(0, nodes, xi, eta, zeta);
    amrex::Real df0dzeta = dfdzeta(0, nodes, xi, eta, zeta);

    amrex::Real df1dxi = dfdxi(1, nodes, xi, eta, zeta);
    amrex::Real df1deta = dfdeta(1, nodes, xi, eta, zeta);
    amrex::Real df1dzeta = dfdzeta(1, nodes, xi, eta, zeta);

    amrex::Real df2dxi = dfdxi(2, nodes, xi, eta, zeta);
    amrex::Real df2deta = dfdeta(2, nodes, xi, eta, zeta);
    amrex::Real df2dzeta = dfdzeta(2, nodes, xi, eta, zeta);

    amrex::Real detJ = df0dxi * (df1deta * df2dzeta - df1dzeta * df2deta) -
                       df0deta * (df1dxi * df2dzeta - df1dzeta * df2dxi) +
                       df0dzeta * (df1dxi * df2deta - df1deta * df2dxi);

    amrex::Real detJ_xi = f0 * (df1deta * df2dzeta - df1dzeta * df2deta) -
                          df0deta * (f1 * df2dzeta - df1dzeta * f2) +
                          df0dzeta * (f1 * df2deta - df1deta * f2);

    amrex::Real detJ_eta = df0dxi * (f1 * df2dzeta - df1dzeta * f2) -
                           f0 * (df1dxi * df2dzeta - df1dzeta * df2dxi) +
                           df0dzeta * (df1dxi * f2 - f1 * df2dxi);

    amrex::Real detJ_zeta = df0dxi * (df1deta * f2 - f1 * df2deta) -
                            df0deta * (df1dxi * f2 - f1 * df2dxi) +
                            f0 * (df1dxi * df2deta - df1deta * df2dxi);
    amrex::Real new_xi = xi - detJ_xi / detJ;
    amrex::Real new_eta = eta - detJ_eta / detJ;
    amrex::Real new_zeta = zeta - detJ_zeta / detJ;

    lc += 1;
    err = amrex::max(
      amrex::Math::abs(xi - new_xi), amrex::Math::abs(eta - new_eta),
      amrex::Math::abs(zeta - new_zeta));

    xi = new_xi;
    eta = new_eta;
    zeta = new_zeta;
  }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE void
fe_interp(
  const amrex::RealVect& pos,
  const int ip,
  const int jp,
  const int kp,
  const amrex::RealVect& dx,
  const amrex::RealVect& dxi,
  const amrex::RealVect& plo,
  amrex::Array4<const amrex::EBCellFlag> const& flags,
  amrex::Array4<const amrex::Real> const& ccent,
  amrex::Array4<const amrex::Real> const& bcent,
  amrex::Array4<const amrex::Real> const& bnorm,
  amrex::Array4<const amrex::Real> const& vfrac,
  amrex::IntVect* indx_array,
  amrex::Real* weights)
{

  const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

  const amrex::Real cdist_x =
    pos[0] - (ip + 0.5 + ccent(ip, jp, kp, 0)) * dx[0] - plo[0];
  const amrex::Real cdist_y =
    pos[1] - (jp + 0.5 + ccent(ip, jp, kp, 1)) * dx[1] - plo[1];
  const amrex::Real cdist_z =
    pos[2] - (kp + 0.5 + ccent(ip, jp, kp, 2)) * dx[2] - plo[2];
  const amrex::Real vtol = 0.05;

  // Distance between particle and cell centoid.
  const amrex::Real cdist =
    std::sqrt(cdist_x * cdist_x + cdist_y * cdist_y + cdist_z * cdist_z);

  // Before doing anything fancy, just check that the particle isn't overlapping
  // the cell center. If it is, use the cell value.
  if (cdist < tolerance) {
    indx_array[0] = {ip, jp, kp};
    weights[0] = 1.;
    for (int aindx = 1; aindx < 8; ++aindx) {
      indx_array[aindx] = {ip, jp, kp};
      weights[aindx] = 0.;
    }
    return;
  }

  amrex::Real par_dot_EB(2.);
  amrex::Real cent_dot_EB(1.);

  // If the particle is inside a cut-cell, verify that it is on the correct
  // side of the EB before trying to interpolate.
  if (flags(ip, jp, kp).isSingleValued()) {
    const amrex::RealVect normal = {
      -bnorm(ip, jp, kp, 0), -bnorm(ip, jp, kp, 1), -bnorm(ip, jp, kp, 2)};

    // Projection of vector pointing from EB centroid to particle onto EB normal
    par_dot_EB =
      (pos[0] - (ip + 0.5 + bcent(ip, jp, kp, 0)) * dx[0] - plo[0]) *
        normal[0] +
      (pos[1] - (jp + 0.5 + bcent(ip, jp, kp, 1)) * dx[1] - plo[1]) *
        normal[1] +
      (pos[2] - (kp + 0.5 + bcent(ip, jp, kp, 2)) * dx[2] - plo[2]) * normal[2];

    // Projection of vector pointing from EB centroid to cell centroid onto EB
    // normal
    cent_dot_EB =
      (ccent(ip, jp, kp, 0) - bcent(ip, jp, kp, 0)) * dx[0] * normal[0] +
      (ccent(ip, jp, kp, 1) - bcent(ip, jp, kp, 1)) * dx[1] * normal[1] +
      (ccent(ip, jp, kp, 2) - bcent(ip, jp, kp, 2)) * dx[2] * normal[2];
    // Temporary sanity check
    AMREX_ASSERT_WITH_MESSAGE(
      cent_dot_EB > tolerance,
      "cent_dot_EB < tolerance ... this makes no sense!");
  }
  // Particle position relative to cell center [-0.5, 0.5]
  amrex::Real gx = (pos[0] - plo[0]) * dxi[0] - (ip + 0.5);
  amrex::Real gy = (pos[1] - plo[1]) * dxi[1] - (jp + 0.5);
  amrex::Real gz = (pos[2] - plo[2]) * dxi[2] - (kp + 0.5);

  // Use the centroid location of the cell containing the particle
  // to guess the interpolation stencil.
  const int i = (gx < ccent(ip, jp, kp, 0)) ? ip : ip + 1;
  const int j = (gy < ccent(ip, jp, kp, 1)) ? jp : jp + 1;
  const int k = (gz < ccent(ip, jp, kp, 2)) ? kp : kp + 1;

  const int di = i - ip; // 0 or 1
  const int dj = j - jp; // 0 or 1
  const int dk = k - kp; // 0 or 1

  // Count the number of non-connected cells in the stencil
  int covered = 0;
  for (int kk(-1); kk < 1; kk++) {
    for (int jj(-1); jj < 1; jj++) {
      for (int ii(-1); ii < 1; ii++) {
        if (not flags(ip, jp, kp).isConnected(di + ii, dj + jj, dk + kk))
          covered += 1;
      }
    }
  }
  // A negative value implies that the particle is 'behind' the EB (and
  // therefore inside the wall). Although this shouldn't occur often, it could
  // for fast moving particles.  Set the fluid velocity to zero and use the
  // cell value for remaining items.
  if (
    (flags(ip, jp, kp).isSingleValued() and (par_dot_EB <= tolerance)) ||
    flags(ip, jp, kp).isCovered()) {
    Abort("Particle has penetrated EB boundary");
    // The particle is near the EB. It is either
    // 1) between the cell centroid and the EB, or
    // 2) one or more cells in the 2^3 stencil is covered.
    // Either way, interpolating to the particle is not straight forward so
    // (for now) we do a 1D interpolation to the EB.

  } else if ((covered > 0) or (par_dot_EB < cent_dot_EB)) {
    for (int aindx = 0; aindx < 8; ++aindx) {
      indx_array[aindx] = {ip, jp, kp};
      weights[aindx] = 0.;
    }
    weights[0] = 1.;
  } else {
    amrex::GpuArray<amrex::GpuArray<amrex::Real, 3>, 8> nodes;

    // Here we pack the local 2x2x2 stencil into a local array. Note that
    // the node ordering is not consistent with the previous version.
    // Nodes 2 and 3 have been swapped as well as nodes 6 and 7. This was
    // to allow for a more compact for-loop filling.

    int lc(0);

    for (int kk(-1); kk < 1; kk++) {
      for (int jj(-1); jj < 1; jj++) {
        for (int ii(-1); ii < 1; ii++) {
          nodes[lc][0] =
            (i + ii + 0.5 + ccent(i + ii, j + jj, k + kk, 0)) * dx[0];
          nodes[lc][1] =
            (j + jj + 0.5 + ccent(i + ii, j + jj, k + kk, 1)) * dx[1];
          nodes[lc][2] =
            (k + kk + 0.5 + ccent(i + ii, j + jj, k + kk, 2)) * dx[2];
          lc += 1;
        }
      }
    }

    amrex::Real xi = (pos[0] - nodes[0][0] - plo[0]) * dxi[0];
    amrex::Real eta = (pos[1] - nodes[0][1] - plo[1]) * dxi[1];
    amrex::Real zeta = (pos[2] - nodes[0][2] - plo[2]) * dxi[2];

    get_interp_mapping(pos - plo, nodes, xi, eta, zeta);
    indx_array[0] = {i - 1, j - 1, k - 1};
    indx_array[1] = {i, j - 1, k - 1};
    indx_array[2] = {i, j, k - 1};
    indx_array[3] = {i - 1, j, k - 1};
    indx_array[4] = {i - 1, j - 1, k};
    indx_array[5] = {i, j - 1, k};
    indx_array[6] = {i, j, k};
    indx_array[7] = {i - 1, j, k};
    weights[0] =
      1 - xi - eta - zeta + xi * eta + xi * zeta + eta * zeta - xi * eta * zeta;
    weights[1] = xi * (1 - eta - zeta + eta * zeta);
    weights[2] = xi * eta - xi * eta * zeta;
    weights[3] = eta * (1 - xi - zeta + xi * zeta);
    weights[4] = zeta * (1 - xi - eta + xi * eta);
    weights[5] = xi * zeta - xi * eta * zeta;
    weights[6] = xi * eta * zeta;
    weights[7] = eta * zeta - xi * eta * zeta;
    // Rudimentary attempt to prevent source terms from being added
    // to small volume cells to avoid having to redistribute later
    // Will probably have to add a redistribute function
    amrex::Real rw = 0.;
    for (int aindx = 0; aindx < 8; ++aindx) {
      if (vfrac(indx_array[aindx]) < vtol) {
        weights[aindx] = 0.;
      }
      rw += weights[aindx];
    }
    rw = 1. / rw;
    for (int aindx = 0; aindx < 8; ++aindx) {
      weights[aindx] *= rw;
    }
  }
  return;
}

#endif // AMREX_USE_EB
#endif // _SPRAYINTERPOLATION_H_
