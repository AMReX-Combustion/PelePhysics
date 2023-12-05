
#include "SootModel.H"
#include "SootModel_derive.H"

using namespace amrex;

void
soot_largeparticledata(
  const Box& bx,
  FArrayBox& slfab,
  const int dcomp,
  const int /*ncomp*/,
  const FArrayBox& datafab,
  const Geometry& /*geomdata*/,
  const Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  SootConst sc;
  const Real V0 = sc.V0 / std::pow(sc.len_conv, 3);
  const Real S0 = sc.S0 / std::pow(sc.len_conv, 2);
  auto const dat = datafab.array();
  auto sl = slfab.array(dcomp);
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
    Real M00 = dat(i, j, k, 0);
    Real M01 = dat(i, j, k, 2);
    Real M10 = dat(i, j, k, 1);
    Real N0 = dat(i, j, k, NUM_SOOT_MOMENTS);
    Real N_L = M00 - N0;
    // Number density for the second mode (large particles)
    sl(i, j, k, 0) = N_L;
    if (N_L < 1.) {
      sl(i, j, k, 1) = 0.;
      sl(i, j, k, 2) = 0.;
    } else {
      // Mean volume of the second mode
      sl(i, j, k, 1) = (M10 - N0 * V0) / N_L;
      // Mean surface area of the second mode
      sl(i, j, k, 2) = (M01 - N0 * S0) / N_L;
    }
  });
}

void
soot_genvars(
  const Box& bx,
  FArrayBox& slfab,
  const int dcomp,
  const int /*ncomp*/,
  const FArrayBox& datafab,
  const Geometry& /*geomdata*/,
  const Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  SootConst sc;
  const Real sootRho = sc.SootDensity / sc.rho_conv;
  auto const dat = datafab.array();
  auto sl = slfab.array(dcomp);
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    Real rho = dat(i, j, k, 0);
    Real fv = dat(i, j, k, 2);
    // Soot mass is volume * density
    Real soot_mass = fv * sootRho;
    // Soot concentration
    sl(i, j, k, 0) = soot_mass;
    // Total mass
    sl(i, j, k, 1) = rho + soot_mass;
  });
}
