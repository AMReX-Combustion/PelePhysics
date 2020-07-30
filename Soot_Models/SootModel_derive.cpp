
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
                             const Real       time,
                             const int*       bcrec,
                             const int        level)
{
  const Real V0 = SootConst::V0;
  const Real S0 = SootConst::S0;
  auto const dat = datafab.array();
  auto       sl  = slfab.array();
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      Real M00 = dat(i,j,k,UFSOOT);
      Real M01 = dat(i,j,k,UFSOOT + 2)/M00;
      Real M10 = dat(i,j,k,UFSOOT + 1)/M00;
      Real N0 = dat(i,j,k,UFSOOT + NUM_SOOT_VARS - 1);
      // N_L
      sl(i,j,k,dcomp) = M00 - N0;
      // V_L
      sl(i,j,k,dcomp+1) = (M10 - N0*V0)/(M00 - N0 + 1.E-30);
      // S_L
      sl(i,j,k,dcomp+2) = (M01 - N0*S0)/(M00 - N0 + 1.E-30);
    });
}

void soot_genvars (const Box&       bx,
                   FArrayBox&       slfab,
                   const int        dcomp,
                   const int        ncomp,
                   const FArrayBox& datafab,
                   const Geometry&  geomdata,
                   const Real       time,
                   const int*       bcrec,
                   const int        level)
{
  const Real sootRho = SootConst::SootDensity;
  auto const dat = datafab.array();
  auto       sl  = slfab.array();
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      Real rho = dat(i,j,k,URHO);
      Real M00 = dat(i,j,k,UFSOOT);
      Real fv = dat(i,j,k,UFSOOT + 1);
      Real Sf = dat(i,j,k,UFSOOT + 2);
      Real Vs = fv/M00;
      Real Ss = Sf/M00;
      // Soot mass is volume * density
      Real soot_mass = fv*sootRho;
      // Soot concentration
      sl(i,j,k,dcomp) = soot_mass;
      // dp = 6*V/S
      sl(i,j,k,dcomp+1) = 6.*Vs/Ss;
      // np = S^3/(36 pi V^2)
      sl(i,j,k,dcomp+2) = std::pow(Ss, 3)/(36.*M_PI*Vs*Vs);
    });
}
