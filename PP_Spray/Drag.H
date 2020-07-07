#ifndef DRAG_H_
#define DRAG_H_

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include <SprayParticles.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void AdjIndexWeights(const RealVect& len,
                     IntVect*        indx_array,
                     Real*           coef,
                     const IntVect&  dom_lo,
                     const IntVect&  dom_hi)
{
  // Compute distances for trilinear interpolation
  IntVect indx(AMREX_D_DECL(amrex::Math::floor(len[0]),
                            amrex::Math::floor(len[1]),
                            amrex::Math::floor(len[2])));
  RealVect delL(AMREX_D_DECL(len[0] - indx[0], len[1] - indx[1], len[2] - indx[2]));
  int iin[2] = {amrex::max(indx[0] - 1, dom_lo[0]),
                amrex::min(indx[0], dom_hi[0])};
  Real ssx[2] = {1. - delL[0], delL[0]};
  int aindx = 0;
#if AMREX_SPACEDIM > 1
  int jjn[2] = {amrex::max(indx[1] - 1, dom_lo[1]),
                amrex::min(indx[1], dom_hi[1])};
  Real ssy[2] = {1. - delL[1], delL[1]};
#if AMREX_SPACEDIM == 3
  int kkn[2] = {amrex::max(indx[2] - 1, dom_lo[2]),
                amrex::min(indx[2], dom_hi[2])};
  Real ssz[2] = {1. - delL[2], delL[2]};
  for (int kk = 0; kk != 2; ++kk) {
    int zindx = kkn[kk];
#endif
    for (int jj = 0; jj != 2; ++jj) {
      int yindx = jjn[jj];
#endif
      for (int ii = 0; ii != 2; ++ii) {
        int xindx = iin[ii];
        AMREX_D_TERM(indx_array[aindx][0] = xindx;,
                     indx_array[aindx][1] = yindx;,
                     indx_array[aindx][2] = zindx;);
        coef[aindx] = AMREX_D_TERM(ssx[ii],*ssy[jj],*ssz[kk]);
        ++aindx;
      }
#if AMREX_SPACEDIM > 1
    }
#if AMREX_SPACEDIM == 3
  }
#endif
#endif
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void InterpolateFaceVelocity(const RealVect& len,
                             const IntVect&  dom_lo,
                             const IntVect&  dom_hi,
                             GpuArray<Array4<Real const>, AMREX_SPACEDIM> const& umac,
                             RealVect&       vel_fluid)
{
  for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
    RealVect face_len(AMREX_D_DECL(len[0] + (dir == 0)*0.5,
                                   len[1] + (dir == 1)*0.5,
                                   len[2] + (dir == 2)*0.5));
    // Compute distances for trilinear interpolation
    IntVect indx(AMREX_D_DECL(amrex::Math::floor(face_len[0]),
                              amrex::Math::floor(face_len[1]),
                              amrex::Math::floor(face_len[2])));
    RealVect delL(AMREX_D_DECL(face_len[0] - indx[0],
                               face_len[1] - indx[1],
                               face_len[2] - indx[2]));
    int iin[2] = {amrex::max(indx[0] - 1, dom_lo[0]),
                  amrex::min(indx[0], dom_hi[0])};
    Real ssx[2] = {1. - delL[0], delL[0]};
    int zindx = 0;
    int yindx = 0;
    int xindx = 0;
#if AMREX_SPACEDIM > 1
    int jjn[2] = {amrex::max(indx[1] - 1, dom_lo[1]),
                  amrex::min(indx[1], dom_hi[1])};
    Real ssy[2] = {1. - delL[1], delL[1]};
#if AMREX_SPACEDIM == 3
    int kkn[2] = {amrex::max(indx[2] - 1, dom_lo[2]),
                  amrex::min(indx[2], dom_hi[2])};
    Real ssz[2] = {1. - delL[2], delL[2]};
    for (int kk = 0; kk != 2; ++kk) {
      zindx = kkn[kk];
#endif
      for (int jj = 0; jj != 2; ++jj) {
        yindx = jjn[jj];
#endif
        for (int ii = 0; ii != 2; ++ii) {
          xindx = iin[ii];
          vel_fluid[dir] += (umac[dir])(xindx, yindx, zindx, 0)*
            AMREX_D_TERM(ssx[ii],*ssy[jj],*ssz[kk]);
        }
#if AMREX_SPACEDIM > 1
      }
#if AMREX_SPACEDIM == 3
    }
#endif
#endif
  }
}

// Compute the heat transfer coefficient using the
// corrected Nusselt number and B_T value
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real calcHeatCoeff(const Real& ratio,
                   const Real& B_M,
                   const Real& B_eps,
                   const Real& C_eps,
                   const Real& Nu_0)
{
  if (B_M <= C_eps) return 0.;
  const int maxIter = 100;
  const Real NU2 = Nu_0 - 2.;
  const Real BM1 = 1. + B_M;
  Real phi = ratio/Nu_0;
  Real B_T_old = std::pow(BM1, phi) - 1.;
  Real BT1 = 1. + B_T_old;
  Real logB = std::log(BT1);
  Real invFT = B_T_old/(logB*std::pow(BT1, 0.7));
  Real Nu_num = 2. + NU2*invFT;
  phi = ratio/Nu_num;
  Real B_T = std::pow(BM1, phi) - 1.;
  Real error = std::abs(B_T - B_T_old);
  int k = 0;
  while (k < maxIter && error > B_eps) {
    B_T_old = B_T;
    BT1 = 1. + B_T;
    logB = std::log(BT1);
    invFT = B_T/(logB*std::pow(BT1, 0.7));
    Nu_num = 2. + NU2*invFT;
    phi = ratio/Nu_num;
    B_T = std::pow(BM1, phi) - 1.;
    error = std::abs(B_T - B_T_old);
    k++;
  }
  logB = std::log(1. + B_T);
  invFT = B_T/(logB*std::pow(1. + B_T, 0.7));
  Nu_num = 2. + NU2*invFT;
  return Nu_num*logB/B_T;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
bool checkWall(Real&      pos,
               Real&      vel,
               const Real phi,
               const Real plo,
               const int  hi_bound,
               const int  lo_bound)
{
  if (pos > phi) {
    switch (hi_bound) {
    case 1 :
      pos *= -1.;
      pos += 2.*phi;
      vel *= -1.;
      break;
    case -1 :
      return true;
      break;
    default :
      break;
    }
  }
  if (pos < plo) {
    switch (lo_bound) {
    case 1 :
      pos *= -1.;
      pos += 2.*plo;
      vel *= -1.;
      break;
    case -1 :
      return true;
      break;
    default :
      break;
    }
  }
  return false;
}

#endif
