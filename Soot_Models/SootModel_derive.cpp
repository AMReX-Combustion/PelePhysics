
#include <AMReX_CONSTANTS.H>
#include <PeleC.H>
#include <SootModel.H>
#include "SootModel_derive.H"

using namespace amrex;

void soot_largeparticledata (const Box&       bx,
                             FArrayBox&       slfab,
                             const int        dcomp,
                             const int        ncomp,
                             const FArrayBox& datafab,
                             const Geometry&  geomdata,
                             const Real&      time,
                             const int*       bcrec,
                             const int        level)
{
  auto const dat = datafab.array();
  auto       sl  = slfab.array();
  const auto lo = lbound(bx);
  const auto hi = ubound(bx);
  for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
        Real M00 = dat(i,j,k,PeleC::FirstSootVar);
        Real M01 = dat(i,j,k,PeleC::FirstSootVar + 2)/M00;
        Real M10 = dat(i,j,k,PeleC::FirstSootVar + 1)/M00;
        Real N0 = dat(i,j,k,PeleC::FirstSootVar + PeleC::NumSootVars - 1);
        // N_L
        sl(i,j,k,dcomp) = M00 - N0;
        // V_L
        sl(i,j,k,dcomp+1) = (M10 - N0*SootModel::m_V0)/(M00 - N0 + 1.E-30);
        // S_L
        sl(i,j,k,dcomp+2) = (M01 - N0*SootModel::m_S0)/(M00 - N0 + 1.E-30);
      }
    }
  }
}

void soot_genvars (const Box&       bx,
                   FArrayBox&       slfab,
                   const int        dcomp,
                   const int        ncomp,
                   const FArrayBox& datafab,
                   const Geometry&  geomdata,
                   const Real&      time,
                   const int*       bcrec,
                   const int        level)
{
  auto const dat = datafab.array();
  auto       sl  = slfab.array();
  const auto lo = lbound(bx);
  const auto hi = ubound(bx);
  for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
        // Density
        Real rho = dat(i,j,k,PeleC::Density);
        Real M00 = dat(i,j,k,PeleC::FirstSootVar);
        Real fv = dat(i,j,k,PeleC::FirstSootVar + 1);
        Real Sf = dat(i,j,k,PeleC::FirstSootVar + 2);
        Real Vs = fv/M00;
        Real Ss = Sf/M00;
        // Soot mass is volume * density
        Real soot_mass = fv*SootModel::m_SootDensity;
        // Soot concentration
        sl(i,j,k,dcomp) = soot_mass;
        // dp = 6*V/S
        sl(i,j,k,dcomp+1) = 6.*Vs/Ss;
        // np = S^3/(36 pi V^2)
        sl(i,j,k,dcomp+2) = std::pow(Ss, 3)/(36.*M_PI*Vs*Vs);
      }
    }
  }
}
