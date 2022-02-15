#ifndef _SOOTREACTIONS_H_
#define _SOOTREACTIONS_H_

#include "Constants_Soot.H"

struct SootReaction
{
  Real SootDensityC;
  Real SootChi;
  // Vectors of the forward and backward pre-exponential rates,
  // temperature exponents, and activation energies/R
  GpuArray<Real, NUM_SOOT_REACT> A_f = {{0.}};
  GpuArray<Real, NUM_SOOT_REACT> n_f = {{0.}};
  GpuArray<Real, NUM_SOOT_REACT> ER_f = {{0.}};
  GpuArray<Real, NUM_SOOT_REACT> A_b = {{0.}};
  GpuArray<Real, NUM_SOOT_REACT> n_b = {{0.}};
  GpuArray<Real, NUM_SOOT_REACT> ER_b = {{0.}};
  // Vector of number of reactants and products
  GpuArray<int, NUM_SOOT_REACT> rNum = {{0}};
  GpuArray<int, NUM_SOOT_REACT> pNum = {{0}};
  // Vector of species indices
  // Maximum 3 species per reaction side
  GpuArray<int, 3 * NUM_SOOT_REACT> nIndx_f = {{0}};
  GpuArray<int, 3 * NUM_SOOT_REACT> nIndx_b = {{0}};
  // Soot surface reaction reference indexing
  // Forward and backward reaction soot indices
  // TODO: Assumes only 1 mole of soot is on each side
  GpuArray<int, NUM_SOOT_REACT> sIndx_f = {{0}};
  GpuArray<int, NUM_SOOT_REACT> sIndx_b = {{0}};
  // Vector of stoichiometric coefficients
  GpuArray<Real, 3 * NUM_SOOT_REACT> nu_f = {{0.}};
  GpuArray<Real, 3 * NUM_SOOT_REACT> nu_b = {{0.}};

  // Compute the dimerization rate
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real
  dimerRate(const Real& T, const Real& xi_PAH) const
  {
    // The effective rate of dimerization
    // This should coincide with the last reaction in the reaction lists
    const int fr = NUM_SOOT_REACT - 1;
    return amrex::max(
      0., A_f[fr] * std::pow(T, n_f[fr]) * std::exp(-ER_f[fr] / T) * xi_PAH *
            xi_PAH);
  }

  // Compute the surface and gas phase chemistry rates
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void chemicalSrc(
    const Real& T,
    const Real& surf,
    const Real xi_n[],
    const Real moments[],
    Real& k_sg,
    Real& k_ox,
    Real& k_o2,
    Real omega_src[]) const
  {
    // Number of surface reactions
    const int nsr = NUM_SOOT_REACT;
    GpuArray<Real, NUM_SOOT_REACT> k_fwd;
    GpuArray<Real, NUM_SOOT_REACT> k_bkwd;
    GpuArray<Real, NUM_SOOT_REACT> w_fwd;
    GpuArray<Real, NUM_SOOT_REACT> w_bkwd;
    const Real invT = 1. / T;
    // Loop over reactions
    for (int i = 0; i < nsr; ++i) {
      k_fwd[i] = A_f[i] * std::pow(T, n_f[i]) * std::exp(-ER_f[i] * invT);
      k_bkwd[i] = A_b[i] * std::pow(T, n_b[i]) * std::exp(-ER_b[i] * invT);
      Real fwdM = 1.;
      for (int j = 0; j < rNum[i]; ++j) {
        // Reactant gas species index
        const int rIndx = nIndx_f[3 * i + j];
        fwdM *= std::pow(xi_n[rIndx], nu_f[3 * i + j]);
      }
      w_fwd[i] = k_fwd[i] * fwdM;
      Real bkwdM = 1.;
      for (int j = 0; j < pNum[i]; ++j) {
        // Product gas species index
        const int pIndx = nIndx_b[3 * i + j];
        bkwdM *= std::pow(xi_n[pIndx], nu_b[3 * i + j]);
      }
      w_bkwd[i] = k_bkwd[i] * bkwdM;
    }
    // TODO: This will depend on the surface reactions, currently hardcoded
    Real fSootStar = computeRadSiteConc(xi_n, k_fwd.data(), k_bkwd.data());
    computeSurfRates(w_fwd.data(), w_bkwd.data(), fSootStar, k_sg, k_ox, k_o2);
    // Determine the concentration of hydrogenated and radical soot surface
    // sites Quasi-steady state for surface radical sites on soot
    GpuArray<Real, SootIndx::numSootSpecs> C_Soot;
    C_Soot[SootIndx::indxSootS] = fSootStar * SootDensityC * moments[2];
    C_Soot[SootIndx::indxSootH] = (1. - fSootStar) * SootDensityC * moments[2];
    // TODO: The last two reaction (6 and 7) are special
    // they are not treated the same in this loop
    // TODO: Currently assumes soot is involved in every reaction on both sides
    for (int i = 0; i < nsr - 2; ++i) {
      int sootIndx = sIndx_f[i];
      w_fwd[i] *= C_Soot[sootIndx];
      sootIndx = sIndx_b[i];
      w_bkwd[i] *= C_Soot[sootIndx];
    }
    // TODO: This is specific to the hardcoded 6th reaction, not universal
    w_fwd[5] *= surf;
    // Loop over gas species and solve for omegas
    for (int i = 0; i < nsr; ++i) {
      // Creation/destruction rates for the i-th reaction
      Real rate = w_fwd[i] - w_bkwd[i];
      // Loop over reactant species and subtract nu'_(i,r)*rate
      for (int r = 0; r < rNum[i]; ++r) {
        const int rIndx = nIndx_f[3 * i + r];
        Real nup = nu_f[3 * i + r];
        omega_src[rIndx] -= nup * rate;
      }
      // Loop over product species and add nu"_(i,r)*rate
      for (int p = 0; p < pNum[i]; ++p) {
        const int pIndx = nIndx_b[3 * i + p];
        Real nupp = nu_b[3 * i + p];
        omega_src[pIndx] += nupp * rate;
      }
    }
  }

  // Return fSootStar, fraction of hydrogenated sites
  // that are radical sites
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real computeRadSiteConc(
    const Real xi_n[], const Real k_fwd[], const Real k_bkwd[]) const
  {
    Real C_OH = xi_n[SootGasSpecIndx::indxOH];
    Real C_H = xi_n[SootGasSpecIndx::indxH];
    Real C_H2 = xi_n[SootGasSpecIndx::indxH2];
    Real C_H2O = xi_n[SootGasSpecIndx::indxH2O];
    Real C_C2H2 = xi_n[SootGasSpecIndx::indxC2H2];
    // Factor r for the quasi-steady state concentration of radical sites, r1/r2
    Real r1 = (k_fwd[0] * C_OH + k_fwd[1] * C_H + k_fwd[2]);
    Real r2 =
      (k_bkwd[0] * C_H2O + k_bkwd[1] * C_H2 + k_bkwd[2] * C_H +
       k_fwd[3] * C_C2H2);
    return r1 / (r2 + r1);
  }

  // Compute the surface chemistry rates (1/s)
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void computeSurfRates(
    const Real w_fwd[],
    const Real* /*w_bkwd*/,
    const Real& fSootStar,
    Real& k_sg,
    Real& k_ox,
    Real& k_o2) const
  {
    k_sg = w_fwd[3] * fSootStar;
    k_ox = w_fwd[4] * fSootStar + w_fwd[5] * 0.5 / SootChi;
    k_o2 = w_fwd[4] * fSootStar;
  }
};

#endif
