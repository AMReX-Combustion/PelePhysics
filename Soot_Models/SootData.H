
#ifndef _SOOTDATA_H_
#define _SOOTDATA_H_

#include "EOS.H"
#include "Constants_Soot.H"

using namespace amrex;

struct SootData
{

  using DVReal = Gpu::DeviceVector<Real>;
  SootData (const DVReal& dime6,
            const DVReal& nve3,
            const DVReal& nve6,
            const DVReal& momFact,
            const DVReal& ssfmCF,
            const DVReal& sscnCF,
            const DVReal& smallOF,
            const DVReal& fragFact,
            const Real&   condFact,
            const Real&   lambdaCoagFact)
    : m_dime6(dime6.dataPtr()),
      m_nve3(nve3.dataPtr()),
      m_nve6(nve6.dataPtr()),
      m_momFact(momFact.dataPtr()),
      m_ssfmCF(ssfmCF.dataPtr()),
      m_sscnCF(sscnCF.dataPtr()),
      m_smallOF(smallOF.dataPtr()),
      m_fragFact(fragFact.dataPtr()),
      m_condFact(condFact),
      m_lambdaCF(lambdaCoagFact)
  {}

  // Given an exponent numerator, return dimer_vol^(exp/6)
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real getDimerExp6 (const int exp) const
  {
    int indx = (exp + 3)/2;
    return m_dime6[indx];
  }

  // Given an exponent numerator, return (nucl_vol)^(exp/3)
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real getNuclExp3 (const int exp) const
  {
    int indx = exp + 3;
    return m_nve3[indx];
  }

  // Given an exponent numerator, return (nucl_vol)^(exp/6)
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real getNuclExp6 (const int exp) const
  {
    int indx = exp + 3;
    return m_nve6[indx];
  }

/********************************************************************
  Moment source terms
********************************************************************/

  // Nucleation source term
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void nucleationMomSrc (const Real& betaNucl,
                         const Real& dimerConc,
                         Real        mom_src[]) const
  {
    const Real dimerConc2 = dimerConc*dimerConc;
    for (int i = 0; i != NUM_SOOT_MOMENTS; ++i)
      mom_src[i] += 0.5*betaNucl*dimerConc2*m_momFact[i];
    mom_src[NUM_SOOT_MOMENTS] += 0.5*betaNucl*dimerConc2;
  }

  // Condensation source term
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void condensationMomSrc (const Real& colConst,
                           const Real& dimerConc,
                           const Real  momFV[],
                           Real        mom_src[]) const
  {
    /** Compute condensation source values
        @param colConst Constant for free molecular collisions
        @param dimerConc Concentration of dimer
        @param momFV Vector of factors used in moment interpolation
        @param mom_src Moment source values
    */
    Real weightDelta = momFV[NUM_SOOT_MOMENTS];
    for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
      const Real momV = SootConst::MomOrderV[i];
      const Real momS = SootConst::MomOrderS[i];
      const Real vv1 = momV + 2.*SootConst::SootAv;
      const Real vs1 = momS + 2.*SootConst::SootAs;
      const Real vv2 = momV + SootConst::SootAv;
      const Real vs2 = momS + SootConst::SootAs;
      Real volTerm =
              fracMom(vv1 -  1.,  vs1,  momFV)*getDimerExp6(3)
        +  2.*fracMom(vv2 -  1.,  vs2,  momFV)*getDimerExp6(5)
        +     fracMom(momV - 1., momS,  momFV)*getDimerExp6(7)
        + 0.5*fracMom(vv1 -  2.,  vs1,  momFV)*getDimerExp6(9)
        +     fracMom(vv2 -  2.,  vs2,  momFV)*getDimerExp6(11)
        + 0.5*fracMom(momV - 2., momS,  momFV)*getDimerExp6(13);
      const Real ss3 = momS + 3.*SootConst::SootFitE;
      const Real sv3 = momV - 2.*SootConst::SootFitE;
      const Real ss2 = ss3 + SootConst::SootAs;
      const Real sv2 = sv3 + SootConst::SootAv;
      const Real ss1 = ss3 + 2.*SootConst::SootAs;
      const Real sv1 = sv3 + 2.*SootConst::SootAv;
      const Real surfTerm = fracMom(sv1 - 1., ss1, momFV)*getDimerExp6(3)
        +  2.*fracMom(sv2 - 1., ss2, momFV)*getDimerExp6(5)
        +     fracMom(sv3 - 1., ss3, momFV)*getDimerExp6(7)
        + 0.5*fracMom(sv1 - 2., ss1, momFV)*getDimerExp6(9)
        +     fracMom(sv2 - 2., ss2, momFV)*getDimerExp6(11)
        + 0.5*fracMom(sv3 - 2., ss3, momFV)*getDimerExp6(13);
      mom_src[i] += colConst*(momV*volTerm +
                              SootConst::SootFitC*momS*surfTerm)*dimerConc;
    }
    // Source for the weight of the delta function
    mom_src[NUM_SOOT_MOMENTS] -=
      m_condFact*colConst*dimerConc*weightDelta;
  }

  // Surface growth source term
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void surfaceGrowthMomSrc (const Real& k_sg,
                            const Real  momFV[],
                            Real        mom_src[]) const
  {
    // Index of the weight of the delta function
    const int dwIndx = NUM_SOOT_MOMENTS;
    const Real weightDelta = momFV[dwIndx];
    const Real factor = SootConst::SootDensityC*SootConst::dVol*k_sg;
    for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
      Real fact1 = fracMom(SootConst::MomOrderV[i] - 1.,
                           SootConst::MomOrderS[i] + 1., momFV);
      Real fact2 = fracMom(SootConst::MomOrderV[i] - 1. - 2.*SootConst::SootFitE,
                           SootConst::MomOrderS[i] + 1. + 3.*SootConst::SootFitE, momFV);
      mom_src[i] += (SootConst::MomOrderV[i]*fact1 +
                     SootConst::MomOrderS[i]*SootConst::SootFitC*fact2)*factor;
    }
    // Weight of the delta function
    mom_src[dwIndx] -= SootConst::nuclSurf*k_sg*SootConst::SootDensityC*weightDelta;
  }

  // Oxidation and fragmentation source terms
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void oxidFragMomSrc (const Real& k_ox,
                       const Real& k_o2,
                       const Real  momFV[],
                       Real        mom_src[]) const
  {
    // Index of the weight of the delta function
    const int dwIndx = NUM_SOOT_MOMENTS;
    const Real weightDelta = momFV[dwIndx];
    const Real factOx = k_ox*SootConst::dVol*SootConst::SootDensityC;
    const Real factO2 = 2.*k_o2*SootConst::dVol*SootConst::SootDensityC;
    for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
      // Oxidation of the small particles
      Real small = -factOx*m_smallOF[i]*weightDelta;
      // Oxidation of the larger particles
      Real factLarge = fracMomLarge(SootConst::MomOrderV[i] - 1.,
                                    SootConst::MomOrderS[i] + 1., momFV);
      Real large = -factOx*(SootConst::MomOrderV[i] + 2./3.*SootConst::MomOrderS[i])*factLarge;
      // Add oxidation source
      mom_src[i] += large + small;
      // Add fragmentation source
      mom_src[i] += m_fragFact[i]*factO2*factLarge;
    }
    Real fracLarge = fracMomLarge(-1., 1., momFV);
    Real small = -factOx*m_smallOF[dwIndx]*weightDelta;
    Real inter = SootConst::nuclVol/(fracMomLarge(1., 0., momFV)/
                          fracMomLarge(0., 0., momFV));
    Real large = factOx*inter*fracLarge;
    // Add oxidation source for weight of delta function
    mom_src[dwIndx] += (small + large);
    // Add fragmentation source for weight of delta function
    mom_src[dwIndx] += inter*factO2*fracLarge;
  }

  // Return the dimer concentration
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real dimerization (const Real& convT,
                     const Real& betaNucl,
                     const Real& dimerRate,
                     const Real  momFV[]) const
  {
    // Collision coefficient for condensation
    const Real betaCond = getBetaCond(convT, momFV);
    // Using the following quadratic equation:
    // betaNucl*[DIMER]^2 + betaCond*[DIMER] - dimerRate = 0
    // compute the [DIMER] using the quadratic formula
    // x = -b + sqrt(b^2 - 4ac)/(2a)
    const Real delta = betaCond*betaCond + 4.*betaNucl*dimerRate;
    return (std::sqrt(delta) - betaCond)/(2.*betaNucl);
  }

  // Clip moment values
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void clipMoments (Real moments[]) const
  {
    Real weightDelta = moments[NUM_SOOT_MOMENTS];
    const Real tolV = SootConst::smallWeight*SootConst::nuclVol;
    const Real tolS = SootConst::smallWeight*SootConst::nuclSurf;
    // Check for globally small moments
    if (moments[0] < SootConst::smallWeight || moments[1] < tolV
        || moments[2] < tolS) {
      moments[1] = amrex::max(moments[1], tolV);
      moments[0] = moments[1]/SootConst::nuclVol;
      moments[2] = SootConst::nuclSurf*moments[0];
      weightDelta = moments[0];
    }
    // Check for size of second mode
    if (moments[1] < SootConst::nuclVol*moments[0] ||
        moments[2] < SootConst::nuclSurf*moments[0]) {
      moments[0] = moments[1]/SootConst::nuclVol;
      moments[2] = moments[0]*SootConst::nuclSurf;
    }
#if NUM_SOOT_MOMENTS == 6
    // Check for (co)variance of second mode
    moments[3] = amrex::max(moments[3], moments[1]*moments[1]/moments[0]);
    moments[4] = amrex::max(moments[4], moments[1]*moments[2]/moments[0]);
    moments[5] = amrex::max(moments[5], moments[2]*moments[2]/moments[0]);
#endif
    if (weightDelta < SootConst::smallWeight) {
      for (int i = 0; i != NUM_SOOT_MOMENTS; ++i)
        moments[i] += (SootConst::smallWeight - weightDelta)*
          std::pow(SootConst::nuclVol, SootConst::MomOrderV[i])*
          std::pow(SootConst::nuclSurf, SootConst::MomOrderS[i]);
      weightDelta = SootConst::smallWeight;
    }
    if (weightDelta > moments[0])
      weightDelta = moments[0];
    moments[NUM_SOOT_MOMENTS] = weightDelta;
  }

  // Convert from CGS to mol, clip moments, and convert back
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void momConvClipConv (Real mom[]) const
  {
    SootConst::convertCGStoMol(mom);
    clipMoments(mom);
    SootConst::convertMoltoCGS(mom);
  }

  // Compute the coagulation source term
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void coagulationMomSrc (const Real& colConst,
                          const Real& T,
                          const Real& mu,
                          const Real& rho,
                          const Real& molMass,
                          const Real  momFV[],
                          Real        mom_src[]) const
  {
    // Index of the weight of the delta function
    const int dwIndx = NUM_SOOT_MOMENTS;
    // Free molecular collision coefficient with van der Waals enhancements
    const Real C_fm = 2.2*colConst;
    // Continuum collision coefficient
    const Real C_cn = 8.*EOS::RU*T/(3.*mu);
    // Mean free path for finite nudsen number correction in continuum regimes
    const Real lambda = 3.*mu/rho*std::sqrt(M_PI*molMass/(8.*EOS::RU*T))*m_lambdaCF;
    Real weightDelta2 = std::pow(momFV[dwIndx], 2);
    for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
      // Collisions between two first mode particles
      // Collision model: pure coalescence
      // S_(0+0) = (144*pi)^(1/3)*V0^(2/3)
      // Free molecular regime
      Real ss_fm = C_fm*m_ssfmCF[i]*weightDelta2;
      // Continuum regime
      Real ss_cn = 4.*C_cn*(1. + 1.257*lambda*getNuclExp3(-1))*m_sscnCF[i]*weightDelta2;
      Real prodss = ss_fm*ss_cn;
      // Harmonic mean for transitional regime
      Real ss = (prodss == 0.) ? 0. : prodss/(ss_fm + ss_cn);

      // Collision between a particle in each mode
      // Collision model: "Splashing"
      // S_(i+0) = S_i + delta S
      // delta S = S*delta V / V *2/3*n_p^(-0.2043)
      // delta V = 2*W_C/rho_soot
      // Free molecular regime
      Real sl_fm = C_fm*FMCoagSL(i, momFV);
      // Continuum regime
      Real sl_cn = C_cn*CNCoagSL(i, lambda, momFV);
      Real prodsl = sl_fm*sl_cn;
      // Harmonic mean for transitional regime
      Real sl = (prodsl == 0.) ? 0. : prodsl/(sl_fm + sl_cn);
      // Collision between two second mode particles
      // Collision model: Pure aggregation
      // S_(i+j) = S_i + S_j
      // Free molecular regime
      Real ll_fm = C_fm*FMCoagLL(i, momFV);
      // Continuum regime
      Real ll_cn = C_cn*CNCoagLL(i, lambda, momFV);
      Real prodll = ll_fm*ll_cn;
      // Harmonic mean for transitional regime
      Real ll = (prodll == 0.) ? 0. : prodll/(ll_fm + ll_cn);
      mom_src[i] += (ss + sl + ll);
    }
    // Free molecular regime
    Real ss_fm = -C_fm*weightDelta2*m_ssfmCF[dwIndx];
    // Continuum regime
    Real ss_cn = -4.*C_cn*(1. + 1.257*lambda*getNuclExp3(-1))*weightDelta2;
    // Harmonic mean for transitional regime
    Real prodss = ss_fm*ss_cn;
    Real ss = (prodss == 0.) ? 0. : prodss/(ss_fm + ss_cn);
    // Free molecular regime
    Real sl_fm = C_fm*FMCoagSL(dwIndx, momFV);
    // Continuum regime
    Real sl_cn = C_cn*CNCoagSL(dwIndx, lambda, momFV);
    // Harmonic mean for transitional regime
    Real prodsl = sl_fm*sl_cn;
    Real sl = (prodsl == 0.) ? 0. : prodsl/(sl_fm + sl_cn);
    mom_src[dwIndx] += (ss + sl);
  }

/*********************************************************************
  Moment interpolation functions
*********************************************************************/

  // Compute the moment interpolation array
  /*
    momFV contains factors for interpolating the moments
    It is ordered as the following
    momFV[0-NUM_SOOT_MOMENTS-1] - Corresponding factor for moment interpolation
    momFV[NUM_SOOT_MOMENTS] - Weight of the delta function
    momFV[NUM_SOOT_MOMENTS+1] - modeCoef
    modeCoef signifies the number of modes to be used
    If the moments are effectively zero, modeCoef = 0. and only 1 mode is used
    Otherwise, modeCoef = 1. and both modes are used
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void computeFracMomVect (const Real moments[],
                           Real       momFV[]) const
  {
    // See above for description of modeCoef
    Real modeCoef;
    // Copy over the weight of the delta function
    momFV[NUM_SOOT_MOMENTS] = moments[NUM_SOOT_MOMENTS];
#if NUM_SOOT_MOMENTS == 3
    const Real M00 = moments[0] - m_momFact[0]*moments[3];
    const Real M10 = moments[1] - m_momFact[1]*moments[3];
    const Real M01 = moments[2] - m_momFact[2]*moments[3];
    // If moments are effectively zero, only use one mode
    if (M00 < 1.E-36 || M10 < 1.E-36 || M01 < 1.E-36) {
      // Contribution from only one mode
      momFV[0] = moments[0];
      momFV[1] = moments[1];
      momFV[2] = moments[2];
      modeCoef = 0.;
    } else {
      // Contribution from both modes
      momFV[0] = M00;
      momFV[1] = M10;
      momFV[2] = M01;
      modeCoef = 1.;
    }
#elif NUM_SOOT_MOMENTS == 6
    const Real M00 = moments[0] - m_momFact[0]*moments[6];
    const Real M10 = moments[1] - m_momFact[1]*moments[6];
    const Real M01 = moments[2] - m_momFact[2]*moments[6];
    const Real M20 = moments[3] - m_momFact[3]*moments[6];
    const Real M11 = moments[4] - m_momFact[4]*moments[6];
    const Real M02 = moments[5] - m_momFact[5]*moments[6];
    Real minMom = amrex::min(M00, amrex::min(M10, M01));
    minMom = amrex::min(minMom, amrex::min(M20, amrex::min(M11, M02)));
    // If moments are effectively zero, only use one mode
    if (minMom < 1.E-36) {
      const Real c1 = std::pow(moments[0], -1.5);
      const Real c2 = std::pow(moments[0], 0.5);
      momFV[0] = moments[0];
      momFV[1] = std::pow(moments[1], 2.)*c1*std::pow(moments[3], -0.5);
      momFV[2] = std::pow(moments[2], 2.)*c1*std::pow(moments[5], -0.5);
      momFV[3] = std::pow(moments[3], 0.5)*c2*std::pow(moments[1], -1.);
      momFV[4] = moments[4]*moments[0]/(moments[1]*moments[2]);
      momFV[5] = std::pow(moments[5], 0.5)*c2*std::pow(moments[2], -1.);
      modeCoef = 0.;
    } else {
      const Real c1 = std::pow(M00, -1.5);
      const Real c2 = std::pow(M00, 0.5);
      momFV[0] = M00;
      momFV[1] = std::pow(M10, 2.)*c1*std::pow(M20, -0.5);
      momFV[2] = std::pow(M01, 2.)*c1*std::pow(M02, -0.5);
      momFV[3] = std::pow(M20, 0.5)*c2*std::pow(M10, -1.);
      momFV[4] = M11*M00/(M10*M01);
      momFV[5] = std::pow(M02, 0.5)*c2*std::pow(M01, -1.);
      modeCoef = 1.;
    }
#endif
    momFV[NUM_SOOT_MOMENTS + 1] = modeCoef;
  }

  // Moment interpolation
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real fracMomLarge (const Real volOrd,
                     const Real surfOrd,
                     const Real momFV[]) const
  {
    // Weight of the delta function
    Real dwVal = momFV[NUM_SOOT_MOMENTS];
    Real factor = std::pow(SootConst::nuclVol, volOrd)*
      std::pow(SootConst::nuclSurf, surfOrd);
    // Remove the contribution from the first mode
    Real outMom = fracMom(volOrd, surfOrd, momFV) - dwVal*factor;
    // If the moment is negative, return a small (consistent) value
    if (outMom <= 0. || outMom != outMom)
      return factor*1.E-66;
    return outMom;
  }

  // Moment interpolation
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real fracMom (const Real volOrd,
                const Real surfOrd,
                const Real momFV[]) const
  {
    // If modeCoef = 0.; only first mode is used
    // If modeCoef = 1.; both modes are used
    const Real modeCoef = momFV[NUM_SOOT_MOMENTS + 1];
    Real bothPFact = 0.;
    if (modeCoef > 0.)
      bothPFact = momFV[NUM_SOOT_MOMENTS]*
        std::pow(SootConst::nuclVol, volOrd)*
        std::pow(SootConst::nuclSurf, surfOrd);
#if NUM_SOOT_MOMENTS == 3
    Real peak = std::pow(momFV[0], 1. - volOrd - surfOrd)*
      std::pow(momFV[1], volOrd)*std::pow(momFV[2], surfOrd);
    return bothPFact + peak;
#elif NUM_SOOT_MOMENTS == 6
    Real prod = momFV[0];
    prod *= std::pow(momFV[1], volOrd);
    prod *= std::pow(momFV[2], surfOrd);
    prod *= std::pow(momFV[3], volOrd*volOrd);
    prod *= std::pow(momFV[4], volOrd*surfOrd);
    prod *= std::pow(momFV[5], surfOrd*surfOrd);
    return bothPFact + prod;
#endif
  }

  // Interpolation for the reduced mass term (square root of sum) in the
  // collision kernel for collision between a particle in each mode
  // Only two grid functions used for all moments
  // Limited sensitivity to increasing the number of grid functions
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real psiSL (const Real x,
              const Real y,
              const Real a,
              const Real b,
              const Real momFV[]) const
  {
    const Real weightDelta = momFV[NUM_SOOT_MOMENTS];
    const Real factor = weightDelta*std::pow(SootConst::nuclVol, a + 2./3.*b);
    Real VF[3] = {2.*SootConst::SootAv + x, SootConst::SootAv + x, x};
    Real SF[3] = {2.*SootConst::SootAs + y, SootConst::SootAs + y, y};
    const Real FML_1 = fracMomLarge(VF[0] - 0.5, SF[0], momFV);
    const Real FML_2 = fracMomLarge(VF[1] - 0.5, SF[1], momFV);
    const Real FML_3 = fracMomLarge(VF[2] - 0.5, SF[2], momFV);
    // nuclVolExp6[i] = nuclVol^(2*i - 3)/6
    Real psi1 = factor*(getNuclExp6(-3)*FML_1 + 2.*getNuclExp6(-1)*FML_2 +
                        getNuclExp6(1)*FML_3);
    const Real FPL_1 = fracMomLarge(VF[0] + 0.5, SF[0], momFV);
    const Real FPL_2 = fracMomLarge(VF[1] + 0.5, SF[1], momFV);
    const Real FPL_3 = fracMomLarge(VF[2] + 0.5, SF[2], momFV);
    Real psi2_1 = factor*(getNuclExp6(-3)*FPL_1 + 2.*getNuclExp6(-1)*FPL_2 +
                          getNuclExp6(1)*FPL_3);
    Real psi2_2 = factor*(getNuclExp6(3)*FML_1 + 2.*getNuclExp6(5)*FML_2 +
                          getNuclExp6(7)*FML_3);
    return std::sqrt(psi1*(psi2_1 + psi2_2));
  }

  // Interpolation for the reduced mass term (square root of sum) in the
  // collision kernel for collision between two particles in the second mode
  // Only two grid functions used for all moments
  // Limited sensitivity to increasing the number of grid functions
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real psiLL (const Real x,
              const Real y,
              const Real a,
              const Real b,
              const Real momFV[]) const
  {
    Real VF_xy[3] = {2.*SootConst::SootAv + x, SootConst::SootAv + x, x};
    Real SF_xy[3] = {2.*SootConst::SootAs + y, SootConst::SootAs + y, y};
    Real VF_ab[3] = {a, SootConst::SootAv + a, 2.*SootConst::SootAv + a};
    Real SF_ab[3] = {b, SootConst::SootAs + b, 2.*SootConst::SootAs + b};
    Real xy_M[3] = {fracMomLarge(VF_xy[0] - 0.5, SF_xy[0], momFV),
                    fracMomLarge(VF_xy[1] - 0.5, SF_xy[1], momFV),
                    fracMomLarge(VF_xy[2] - 0.5, SF_xy[2], momFV)};
    Real xy_P[3] = {fracMomLarge(VF_xy[0] + 0.5, SF_xy[0], momFV),
                    fracMomLarge(VF_xy[1] + 0.5, SF_xy[1], momFV),
                    fracMomLarge(VF_xy[2] + 0.5, SF_xy[2], momFV)};
    Real ab_M[3] = {fracMomLarge(VF_ab[0] - 0.5, SF_ab[0], momFV),
                    fracMomLarge(VF_ab[1] - 0.5, SF_ab[1], momFV),
                    fracMomLarge(VF_ab[2] - 0.5, SF_ab[2], momFV)};
    Real ab_P[3] = {fracMomLarge(VF_ab[0] + 0.5, SF_ab[0], momFV),
                    fracMomLarge(VF_ab[1] + 0.5, SF_ab[1], momFV),
                    fracMomLarge(VF_ab[2] + 0.5, SF_ab[2], momFV)};
    Real psi1 = xy_M[0]*ab_M[0] + 2.*xy_M[1]*ab_M[1] + xy_M[2]*ab_M[2];
    Real psi2_1 = xy_P[0]*ab_M[0] + 2.*xy_P[1]*ab_M[1] + xy_P[2]*ab_M[2];
    Real psi2_2 = xy_M[0]*ab_P[0] + 2.*xy_M[1]*ab_P[1] + xy_M[2]*ab_P[2];
    return std::sqrt(psi1*(psi2_1 + psi2_2));
  }

  // Free molecular coagulation source term
  // Small-Large: "Splashing"
  // -Generalized grid function follows terms
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real FMCoagSL (const int  i,
                 const Real momFV[]) const
  {
    // Weight of delta function N0 and M00
    if (i == NUM_SOOT_MOMENTS || i == 0) {
      return -psiSL(0., 0., 0., 0., momFV);
    }
    const Real fact1 = -2.*SootConst::SootFitE;
    const Real fact2 = 3.*SootConst::SootFitE;
    switch (i) {
    case 1: // M10
      return 0.;
    case 2: // M01
      return SootConst::SootFitC*psiSL(fact1 - 1., fact2 + 1., 1., 0., momFV)
        - psiSL(0., 0., 0., 1., momFV);
    case 3: // M20
      return 2.*psiSL(1., 0., 1., 0., momFV);
    case 4: // M11
      return SootConst::SootFitC*psiSL(fact1, fact2 + 1., 1., 0., momFV)
        + psiSL(0., 1., 1., 0., momFV)
        + SootConst::SootFitC*psiSL(fact1 - 1., fact2 + 1., 2., 0., momFV)
        - psiSL(0., 0., 1., 1., momFV);
    case 5: // M02
      return 2.*SootConst::SootFitC*psiSL(fact1 - 1., fact2 + 2., 1., 0., momFV)
        + SootConst::SootFitC*SootConst::SootFitC*psiSL(2.*fact1 - 2., -3.*fact1 + 2.,
                                                        2., 0., momFV)
        - psiSL(0., 0., 0., 2., momFV);
    default:
      Abort("SootModel::FMCoagSL: Moment not contained in number of moments!");
    }
    return 0.;
  }

  // Free molecular coagulation source term
  // Large-Large: Pure aggregation
  // -Generalized grid function follows terms
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real FMCoagLL (const int  i,
                 const Real momFV[]) const
  {
    switch (i) {
    case 0: // M00
      return -0.5*psiLL(0., 0., 0., 0., momFV);
    case 1: // M10
      return 0.;
    case 2: // M01
      return 0.;
    case 3: // M20
      return psiLL(1., 0., 1., 0., momFV);
    case 4: // M11
      return psiLL(1., 0., 0., 1., momFV);
    case 5: // M02
      return psiLL(0., 1., 0., 1., momFV);
    default:
      Abort("SootModel::FMCoagLL: Moment not contained in number of moments!");
    }
    return 0.;
  }

  // Continuum coagulation source terms
  // Small-Large: "Splashing"
  // Large-Large: Pure aggregation
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real CNCoagSL (const int   i,
                 const Real& lambda,
                 const Real  momFV[]) const
  {
    const Real weightDelta = momFV[NUM_SOOT_MOMENTS];
    // Mean free path for finite Knudsen number correction in continuum regime
    if (i == NUM_SOOT_MOMENTS || i == 0) { // N0 or M00
      int n[] = {0, 1, -1, -2};
      Real x = 0.;
      Real y = 0.;
      return -weightDelta*CNCoagSLFunc(n, x, y, lambda, momFV);
    }
    switch (i) {
    case 1: // M10
      return 0.;
    case 2: // M01
      {
        Real p1, p2;
        {
          int n[] = {3, 4, 2, 1};
          Real x = -2.*SootConst::SootFitE - 1.;
          Real y = 3.*SootConst::SootFitE + 1.;
          p1 = SootConst::SootFitC*CNCoagSLFunc(n, x, y, lambda, momFV);
        }
        {
          int n[] = {2, 3, 1, 0};
          p2 = -CNCoagSLFunc(n, 0., 0., lambda, momFV);
        }
        return weightDelta*(p1 + p2);
      }
    case 3: // M20
      {
        int n[] = {3, 4, 2, 1};
        return 2.*weightDelta*CNCoagSLFunc(n, 1., 0., lambda, momFV);
      }
    case 4: // M11
      {
        Real p1, p2, p3, p4;
        {
          int n[] = {3, 4, 2, 1};
          Real x = -2.*SootConst::SootFitE;
          Real y = 3.*SootConst::SootFitE + 1.;
          p1 = SootConst::SootFitC*CNCoagSLFunc(n, x, y, lambda, momFV);
        }
        {
          int n[] = {3, 4, 2, 1};
          p2 = CNCoagSLFunc(n, 0., 1., lambda, momFV);
        }
        {
          int n[] = {6, 7, 5, 4};
          Real x = -2.*SootConst::SootFitE - 1.;
          Real y = 3.*SootConst::SootFitE + 1.;
          p3 = SootConst::SootFitC*CNCoagSLFunc(n, x, y, lambda, momFV);
        }
        {
          int n[] = {5, 6, 4, 3};
          p4 = -CNCoagSLFunc(n, 0., 0., lambda, momFV);
        }
        return weightDelta*(p1 + p2 + p3 + p4);
      }
    case 5: // M02
      {
        Real p1, p2, p3;
        {
          int n[] = {3, 4, 2, 1};
          Real x = -2.*SootConst::SootFitE - 1.;
          Real y = 3.*SootConst::SootFitE + 2.;
          p1 = 2.*SootConst::SootFitC*CNCoagSLFunc(n, x, y, lambda, momFV);
        }
        {
          int n[] = {6, 7, 5, 4};
          Real x = -4.*SootConst::SootFitE - 2.;
          Real y = 6.*SootConst::SootFitE + 2.;
          p2 = SootConst::SootFitC*SootConst::SootFitC*CNCoagSLFunc(n, x, y, lambda, momFV);
        }
        {
          int n[] = {4, 5, 3, 2};
          p3 = -CNCoagSLFunc(n, 0., 0., lambda, momFV);
        }
        return 2.*weightDelta*(p1 + p2 + p3);
      }
    }
    return 0.;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real CNCoagSLFunc (int         n[4],
                     const Real  x,
                     const Real  y,
                     const Real& lambda,
                     const Real  momFV[]) const
  {
    Real xy_1 = fracMomLarge(x, y, momFV);
    Real xy_2 = fracMomLarge(x - SootConst::SootAv, y - SootConst::SootAs, momFV);
    Real xy_3 = fracMomLarge(x + SootConst::SootAv, y + SootConst::SootAs, momFV);
    Real xy_4 = fracMomLarge(x - 2.*SootConst::SootAv, y - 2.*SootConst::SootAs, momFV);
    Real n_1 = getNuclExp3(n[0]);
    Real n_2 = getNuclExp3(n[1]);
    Real n_3 = getNuclExp3(n[2]);
    Real n_4 = getNuclExp3(n[3]);
    return 2.*xy_1*n_1 + xy_2*n_2 + xy_3*n_3 +
      1.257*lambda*(xy_1*n_3 + xy_2*n_1 + xy_3*n_4 + xy_4*n_2);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real CNCoagLL (const int   i,
                 const Real& lambda,
                 const Real  momFV[]) const
  {
    switch (i) {
    case 0: // M00
      return -0.5*CNCoagLLFunc(0., 0., lambda, momFV);
    case 1: // M10
      return 0.;
    case 2: // M01
      return 0.;
    case 3: // M20
      return CNCoagLLFunc(1., 0., lambda, momFV);
    case 4: // M11
      return CNCoagLLFunc(1., 0., 0., 1., lambda, momFV);
    case 5: // M02
      return CNCoagLLFunc(0., 1., lambda, momFV);
    }
    return 0.;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real CNCoagLLFunc (const Real  x,
                     const Real  y,
                     const Real& lambda,
                     const Real  momFV[]) const
  {
    const Real stav = SootConst::SootAv;
    const Real stas = SootConst::SootAs;
    Real xy_1 = fracMomLarge(x, y, momFV);
    Real xy_2 = fracMomLarge(x -    stav, y -    stas, momFV);
    Real xy_3 = fracMomLarge(x +    stav, y +    stas, momFV);
    Real xy_4 = fracMomLarge(x - 2.*stav, y - 2.*stas, momFV);
    return 2.*xy_1*xy_1 + xy_2*xy_3 + xy_3*xy_2
      + 1.257*lambda*(xy_1*xy_2 + xy_2*xy_1 + xy_3*xy_4 + xy_4*xy_3);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real CNCoagLLFunc(const Real  x,
                    const Real  y,
                    const Real  a,
                    const Real  b,
                    const Real& lambda,
                    const Real  momFV[]) const
  {
    const Real stav = SootConst::SootAv;
    const Real stas = SootConst::SootAs;
    Real xy_1 = fracMomLarge(x,           y,           momFV);
    Real xy_2 = fracMomLarge(x -    stav, y -    stas, momFV);
    Real xy_3 = fracMomLarge(x +    stav, y +    stas, momFV);
    Real xy_4 = fracMomLarge(x - 2.*stav, y - 2.*stas, momFV);
    Real ab_1 = fracMomLarge(a, b, momFV);
    Real ab_2 = fracMomLarge(a -    stav, b -    stas, momFV);
    Real ab_3 = fracMomLarge(a +    stav, b +    stas, momFV);
    Real ab_4 = fracMomLarge(a - 2.*stav, b - 2.*stas, momFV);
    return  2.*ab_1*xy_1 + ab_2*xy_3 + ab_3*xy_2
      + 1.257*lambda*(ab_1*xy_2 + ab_2*xy_1 + ab_3*xy_4 + ab_4*xy_3);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real getBetaCond (const Real& convT,
                    const Real  momFV[]) const
  {
    // Collision frequency between two dimer in the free
    // molecular regime WITHOUT van der Waals enhancement
    // Units: 1/s
    const Real stav = SootConst::SootAv;
    const Real stas = SootConst::SootAs;
    const Real Cfm = SootConst::colFactPi23*convT*SootConst::colFact16*SootConst::avogadros;
    const Real SN = fracMom(  2.*stav,  2.*stas, momFV)*getDimerExp6(-3)
      +          2.*fracMom(     stav,     stas, momFV)*getDimerExp6(-1)
      +             fracMom(       0.,       0., momFV)*getDimerExp6(1)
      +         0.5*fracMom(2.*stav-1., 2.*stas, momFV)*getDimerExp6(3)
      +             fracMom(   stav-1.,    stas, momFV)*getDimerExp6(5)
      +         0.5*fracMom(       -1.,      0., momFV)*getDimerExp6(7);
    return Cfm*SN;
  }

private:
  const Real * m_dime6;
  const Real * m_nve3;
  const Real * m_nve6;
  const Real * m_momFact;
  const Real * m_ssfmCF;
  const Real * m_sscnCF;
  const Real * m_smallOF;
  const Real * m_fragFact;
  const Real m_condFact;
  const Real m_lambdaCF;
};

class SootDataContainer
{
  using DVReal = Gpu::DeviceVector<Real>;
  using HVReal = Gpu::HostVector<Real>;

public:

  SootDataContainer () {}

  ~SootDataContainer () {}

  void build (const HVReal& dimerExp6,
              const HVReal& nuclVolExp3,
              const HVReal& nuclVolExp6,
              const HVReal& momFact,
              const HVReal& ssfmCoagFact,
              const HVReal& sscnCoagFact,
              const HVReal& smallOxidFact,
              const HVReal& fragFact,
              const Real&   condFact,
              const Real&   lambdaCoagFact)
  {
    m_dime6.resize(9);
    Gpu::copy(Gpu::hostToDevice, dimerExp6.begin(), dimerExp6.end(),
              m_dime6.begin());
    m_nve3.resize(11);
    Gpu::copy(Gpu::hostToDevice, nuclVolExp3.begin(), nuclVolExp3.end(),
              m_nve3.begin());
    m_nve6.resize(11);
    Gpu::copy(Gpu::hostToDevice, nuclVolExp6.begin(), nuclVolExp6.end(),
              m_nve6.begin());
    m_momFact.resize(NUM_SOOT_MOMENTS);
    Gpu::copy(Gpu::hostToDevice, momFact.begin(), momFact.end(),
              m_momFact.begin());
    m_ssfmCF.resize(NUM_SOOT_VARS);
    Gpu::copy(Gpu::hostToDevice, ssfmCoagFact.begin(), ssfmCoagFact.end(),
              m_ssfmCF.begin());
    m_sscnCF.resize(NUM_SOOT_MOMENTS);
    Gpu::copy(Gpu::hostToDevice, sscnCoagFact.begin(), sscnCoagFact.end(),
              m_sscnCF.begin());
    m_smallOF.resize(NUM_SOOT_VARS);
    Gpu::copy(Gpu::hostToDevice, smallOxidFact.begin(), smallOxidFact.end(),
              m_smallOF.begin());
    m_fragFact.resize(NUM_SOOT_MOMENTS);
    Gpu::copy(Gpu::hostToDevice, fragFact.begin(), fragFact.end(),
              m_fragFact.begin());
    m_condFact = condFact;
    m_lambdaCoagFact = lambdaCoagFact;
  }

  SootData getSootData () const
  {
    return SootData(m_dime6, m_nve3, m_nve6, m_momFact, m_ssfmCF,
                    m_sscnCF, m_smallOF, m_fragFact, m_condFact,
                    m_lambdaCoagFact);
  }

private:
  DVReal m_dime6;
  DVReal m_nve3;
  DVReal m_nve6;
  DVReal m_momFact;
  DVReal m_ssfmCF;
  DVReal m_sscnCF;
  DVReal m_smallOF;
  DVReal m_fragFact;
  Real m_condFact;
  Real m_lambdaCoagFact;
};

#endif
