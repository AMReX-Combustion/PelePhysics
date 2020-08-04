
#ifndef _SOOTREACTIONS_H_
#define _SOOTREACTIONS_H_

#include "EOS.H"
#include "Constants_Soot.H"

using namespace amrex;

struct SootReactions
{
  using DVReal = Gpu::DeviceVector<Real>;
  using DVInt = Gpu::DeviceVector<int>;

  SootReactions (const DVReal& A_fin,
                 const DVReal& n_fin,
                 const DVReal& ER_fin,
                 const DVReal& A_bin,
                 const DVReal& n_bin,
                 const DVReal& ER_bin,
                 const DVInt&  rNumin,
                 const DVInt&  pNumin,
                 const DVInt&  sIndx_fin,
                 const DVInt&  sIndx_bin,
                 const DVInt&  nIndx_fin,
                 const DVInt&  nIndx_bin,
                 const DVReal& nu_fin,
                 const DVReal& nu_bin)
  : A_f(A_fin.dataPtr()),
    n_f(n_fin.dataPtr()),
    ER_f(ER_fin.dataPtr()),
    A_b(A_bin.dataPtr()),
    n_b(n_bin.dataPtr()),
    ER_b(ER_bin.dataPtr()),
    rNum(rNumin.dataPtr()),
    pNum(pNumin.dataPtr()),
    sIndx_f(sIndx_fin.dataPtr()),
    sIndx_b(sIndx_bin.dataPtr()),
    nIndx_f(nIndx_fin.dataPtr()),
    nIndx_b(nIndx_bin.dataPtr()),
    nu_f(nu_fin.dataPtr()),
    nu_b(nu_bin.dataPtr())
  {}

  // Compute the dimerization rate
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real dimerRate (const Real& T,
                  const Real& xi_PAH) const
  {
    // The effective rate of dimerization
    // This should coincide with the last reaction in the reaction lists
    const int fr = NUM_SOOT_REACT - 1;
    return amrex::max(0., A_f[fr]*std::pow(T, n_f[fr])*std::exp(-ER_f[fr]/T)*xi_PAH*xi_PAH);
  }

  // Compute the surface and gas phase chemistry rates
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void chemicalSrc (const Real& T,
                    const Real& surf,
                    const Real  xi_n[],
                    const Real  moments[],
                    Real        k_fwd[],
                    Real        k_bkwd[],
                    Real        w_fwd[],
                    Real        w_bkwd[],
                    Real&       k_sg,
                    Real&       k_ox,
                    Real&       k_o2,
                    Real        omega_src[]) const
  {
    // Number of surface reactions
    const int nsr = NUM_SOOT_REACT;
    const Real invT = 1./T;
    // Loop over reactions
    for (int i = 0; i != nsr; ++i) {
      k_fwd[i] = A_f[i]*std::pow(T, n_f[i])*std::exp(-ER_f[i]*invT);
      k_bkwd[i] = A_b[i]*std::pow(T, n_b[i])*std::exp(-ER_b[i]*invT);
      Real fwdM = 1.;
      for (int j = 0; j != rNum[i]; ++j) {
        // Reactant gas species index
        const int rIndx = nIndx_f[3*i+j];
        fwdM *= std::pow(xi_n[rIndx], nu_f[3*i+j]);
      }
      w_fwd[i] = k_fwd[i]*fwdM;
      Real bkwdM = 1.;
      for (int j = 0; j != pNum[i]; ++j) {
        // Product gas species index
        const int pIndx = nIndx_b[3*i+j];
        bkwdM *= std::pow(xi_n[pIndx], nu_b[3*i+j]);
      }
      w_bkwd[i] = k_bkwd[i]*bkwdM;
    }
    // TODO: This will depend on the surface reactions, currently hardcoded
    Real fSootStar = computeRadSiteConc(xi_n, k_fwd, k_bkwd);
    computeSurfRates(w_fwd, w_bkwd, fSootStar, k_sg, k_ox, k_o2);
    // Determine the concentration of hydrogenated and radical soot surface sites
    // Quasi-steady state for surface radical sites on soot
    Real C_Soot[SootConst::SootIndx::numSootSpecs];
    C_Soot[SootConst::SootIndx::indxSootS] =
      fSootStar*SootConst::SootDensityC*moments[2];
    C_Soot[SootConst::SootIndx::indxSootH] =
      (1. - fSootStar)*SootConst::SootDensityC*moments[2];
    // TODO: The last two reaction (6 and 7) are special
    // they are not treated the same in this loop
    // TODO: Currently assumes soot is involved in every reaction on both sides
    for (int i = 0; i != nsr - 2; ++i) {
      int sootIndx = sIndx_f[i];
      w_fwd[i] *= C_Soot[sootIndx];
      sootIndx = sIndx_b[i];
      w_bkwd[i] *= C_Soot[sootIndx];
    }
    // TODO: This is specific to the hardcoded 6th reaction, not universal
    w_fwd[5] *= surf;
    // Loop over gas species and solve for omegas
    for (int i = 0; i != nsr; ++i) {
      // Creation/destruction rates for the i-th reaction
      Real rate = w_fwd[i] - w_bkwd[i];
      // Loop over reactant species and subtract nu'_(i,r)*rate
      for (int r = 0; r != rNum[i]; ++r) {
        const int rIndx = nIndx_f[3*i+r];
        Real nup = nu_f[3*i+r];
        omega_src[rIndx] -= nup*rate;
      }
      // Loop over product species and add nu"_(i,r)*rate
      for (int p = 0; p != pNum[i]; ++p) {
        const int pIndx = nIndx_b[3*i+p];
        Real nupp = nu_b[3*i+p];
        omega_src[pIndx] += nupp*rate;
      }
    }
  }

  // Return fSootStar, fraction of hydrogenated sites
  // that are radical sites
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real computeRadSiteConc (const Real xi_n[],
                           const Real k_fwd[],
                           const Real k_bkwd[]) const
  {
    Real C_OH = xi_n[SootConst::GasSpecIndx::indxOH];
    Real C_H = xi_n[SootConst::GasSpecIndx::indxH];
    Real C_H2 = xi_n[SootConst::GasSpecIndx::indxH2];
    Real C_H2O = xi_n[SootConst::GasSpecIndx::indxH2O];
    Real C_C2H2 = xi_n[SootConst::GasSpecIndx::indxC2H2];
    // Factor r for the quasi-steady state concentration of radical sites, r1/r2
    Real r1 = (k_fwd[0]*C_OH + k_fwd[1]*C_H + k_fwd[2]);
    Real r2 = (k_bkwd[0]*C_H2O + k_bkwd[1]*C_H2 + k_bkwd[2]*C_H
               + k_fwd[3]*C_C2H2);
    return r1/(r2 + r1);
  }

  // Compute the surface chemistry rates (1/s)
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void computeSurfRates (const Real  w_fwd[],
                         const Real  w_bkwd[],
                         const Real& fSootStar,
                         Real&       k_sg,
                         Real&       k_ox,
                         Real&       k_o2) const
  {
    k_sg = w_fwd[3]*fSootStar;
    k_ox = w_fwd[4]*fSootStar + w_fwd[5]*0.5/SootConst::SootChi;
    k_o2 = w_fwd[4]*fSootStar;
  }

private:
  const Real * A_f;
  const Real * n_f;
  const Real * ER_f;
  const Real * A_b;
  const Real * n_b;
  const Real * ER_b;
  const int * rNum;
  const int * pNum;
  const int * sIndx_f;
  const int * sIndx_b;
  const int * nIndx_f;
  const int * nIndx_b;
  const Real * nu_f;
  const Real * nu_b;
};

class SootReactionContainer
{

  using DVReal = Gpu::DeviceVector<Real>;
  using DVInt = Gpu::DeviceVector<int>;
  using HVReal = Gpu::HostVector<Real>;
  using HVInt = Gpu::HostVector<int>;

public:

  SootReactionContainer () {}

  ~SootReactionContainer () {}

  void build (const HVReal& A_fin,
              const HVReal& n_fin,
              const HVReal& ER_fin,
              const HVReal& A_bin,
              const HVReal& n_bin,
              const HVReal& ER_bin,
              const HVInt&  rNumin,
              const HVInt&  pNumin,
              const HVInt&  sIndx_fin,
              const HVInt&  sIndx_bin,
              const HVInt&  nIndx_fin,
              const HVInt&  nIndx_bin,
              const HVReal& nu_fin,
              const HVReal& nu_bin)
  {
    const int nsr = NUM_SOOT_REACT;
    const int nsr_sp = nsr*3;
    A_f.resize(nsr);
    n_f.resize(nsr);
    ER_f.resize(nsr);
    A_b.resize(nsr);
    n_b.resize(nsr);
    ER_b.resize(nsr);
    rNum.resize(nsr);
    pNum.resize(nsr);
    sIndx_f.resize(nsr);
    sIndx_b.resize(nsr);
    nIndx_f.resize(nsr_sp);
    nIndx_b.resize(nsr_sp);
    nu_f.resize(nsr_sp);
    nu_b.resize(nsr_sp);
    Gpu::copy(Gpu::hostToDevice, A_fin.begin(), A_fin.end(),
              A_f.begin());
    Gpu::copy(Gpu::hostToDevice, A_bin.begin(), A_bin.end(),
              A_b.begin());
    Gpu::copy(Gpu::hostToDevice, n_fin.begin(), n_fin.end(),
              n_f.begin());
    Gpu::copy(Gpu::hostToDevice, n_bin.begin(), n_bin.end(),
              n_b.begin());
    Gpu::copy(Gpu::hostToDevice, ER_fin.begin(), ER_fin.end(),
              ER_f.begin());
    Gpu::copy(Gpu::hostToDevice, ER_bin.begin(), ER_bin.end(),
              ER_b.begin());
    Gpu::copy(Gpu::hostToDevice, rNumin.begin(), rNumin.end(),
              rNum.begin());
    Gpu::copy(Gpu::hostToDevice, pNumin.begin(), pNumin.end(),
              pNum.begin());
    Gpu::copy(Gpu::hostToDevice, sIndx_fin.begin(), sIndx_fin.end(),
              sIndx_f.begin());
    Gpu::copy(Gpu::hostToDevice, sIndx_bin.begin(), sIndx_bin.end(),
              sIndx_b.begin());
    Gpu::copy(Gpu::hostToDevice, nIndx_fin.begin(), nIndx_fin.end(),
              nIndx_f.begin());
    Gpu::copy(Gpu::hostToDevice, nIndx_bin.begin(), nIndx_bin.end(),
              nIndx_b.begin());
    Gpu::copy(Gpu::hostToDevice, nu_fin.begin(), nu_fin.end(),
              nu_f.begin());
    Gpu::copy(Gpu::hostToDevice, nu_bin.begin(), nu_bin.end(),
              nu_b.begin());
  }

  SootReactions getSootReactions () const
  {
    return SootReactions(A_f, n_f, ER_f, A_b, n_b, ER_b, rNum, pNum,
                         sIndx_f, sIndx_b, nIndx_f, nIndx_b, nu_f, nu_b);
  }

private:
  DVReal A_f;
  DVReal n_f;
  DVReal ER_f;
  DVReal A_b;
  DVReal n_b;
  DVReal ER_b;
  DVInt  rNum;
  DVInt  pNum;
  DVInt  sIndx_f;
  DVInt  sIndx_b;
  DVInt  nIndx_f;
  DVInt  nIndx_b;
  DVReal nu_f;
  DVReal nu_b;
};

#endif
