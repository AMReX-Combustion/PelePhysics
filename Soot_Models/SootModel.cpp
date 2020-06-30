
// Standard library includes
#include <string>
#include <map>

// AMReX include statements
#include <AMReX_Utility.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_FArrayBox.H>

// PeleC include statements
#include <PeleC.H>
#include <PeleC_F.H>

// PeleC-Physics include statements
#include <Fuego_EOS.H>

// PeleC-MP include statements
#include <SootModel.H>

using namespace amrex;

Real SootModel::m_V0 = 0.;
Real SootModel::m_S0 = 0.;
Real SootModel::m_dimerVol = 0.;
Real SootModel::m_nuclVol = 0.;
#if NUM_SOOT_MOMENTS == 3
Real SootModel::MomUnitConv[4] = {0., 0., 0., 0.};
#elif NUM_SOOT_MOMENTS == 6
Real SootModel::MomUnitConv[7] = {0., 0., 0., 0., 0., 0., 0.};
#endif

// Default constructor
SootModel::SootModel()
  :
  m_sootVerbosity(0),
  m_readSootParams(false),
  m_memberDataDefined(false),
  m_conserveMass(true),
  m_SootAv(1. - (2./m_SootDf)),
  m_SootAs(3./m_SootDf - 1.),
  m_PAHindx(-1),
  m_sootVarName(NUM_SOOT_MOMENTS + 1, ""),
  local_soot_dt(1.E18),
  m_maxDtRate(0.2),
  m_reactDataFilled(false),
  m_gasSpecNames(numGasSpecs, ""),
  m_gasSpecRefs(numGasSpecs, -1),
  m_numSurfReacts(-1)
{
  m_sootVarName[NUM_SOOT_MOMENTS] = "soot_N0";
  m_sootVarName[0] = "soot_N";
  m_sootVarName[1] = "soot_fv";
  m_sootVarName[2] = "soot_S";
#if NUM_SOOT_MOMENTS == 6
  m_sootVarName[3] = "soot_V_V";
  m_sootVarName[4] = "soot_V_S";
  m_sootVarName[5] = "soot_S_S";
#endif
}

// Fill in moment source and reaction member data
void
SootModel::define()
{
  // This should be called after readSootParams()
  AMREX_ASSERT(m_readSootParams);

  // From SootModel_react.cpp
  // Initialize reaction and species member data
  initializeReactData();

  // From SootModel_setup.cpp
  // Assign moment factor member data
  defineMemberData();

  // From SootModel_react.cpp
  // Fill surface reaction data
  fillReactionData();
  if (m_sootVerbosity) {
    Print() << "SootModel::define(): Soot model successfully defined"
            << std::endl;
  }
}

// Clip moment values
void
SootModel::clipMoments(Vector<Real>& moments) const
{
  Real weightDelta = moments[NUM_SOOT_MOMENTS];;
  const Real dimer_nbrc = m_dimerVol; // Set dimer size
  // Nucleated particle size
  const Real nucl_nbrc = 2.*dimer_nbrc;
  const Real nucl_surf = std::pow(nucl_nbrc, m_two3rd);
  const Real tolV = m_smallWeight*nucl_nbrc;
  const Real tolS = m_smallWeight*nucl_surf;
  // Check for globally small moments
  if (moments[0] < m_smallWeight || moments[1] < tolV
      || moments[2] < tolS) {
    moments[1] = amrex::max(moments[1], tolV);
    moments[0] = moments[1]/nucl_nbrc;
    moments[2] = nucl_surf*moments[0];
    weightDelta = moments[0];
  }

  // Check for size of second mode
  if (moments[1] < nucl_nbrc*moments[0] ||
      moments[2] < nucl_surf*moments[0]) {
    moments[0] = moments[1]/nucl_nbrc;
    moments[2] = moments[0]*nucl_surf;
  }
#if NUM_SOOT_MOMENTS == 6
  // Check for (co)variance of second mode
  moments[3] = amrex::max(moments[3], moments[1]*moments[1]/moments[0]);
  moments[4] = amrex::max(moments[4], moments[1]*moments[2]/moments[0]);
  moments[5] = amrex::max(moments[5], moments[2]*moments[2]/moments[0]);
#endif
  if (weightDelta < m_smallWeight) {
    for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
      moments[i] += (m_smallWeight - weightDelta)*
        std::pow(nucl_nbrc, MomOrderV[i])*std::pow(nucl_surf, MomOrderS[i]);
    }
    weightDelta = m_smallWeight;
  }
  if (weightDelta > moments[0])
    weightDelta = moments[0];
  moments[NUM_SOOT_MOMENTS] = weightDelta;
}

// Add soot source term
void
SootModel::addSootSourceTerm(const Box&       vbox,
                             const FArrayBox& Ufab,
                             const FArrayBox& Qfab,
                             const FArrayBox& coeff_cc,
                             FArrayBox&       Sfab,
                             const Real       time,
                             const Real       dt) const
{
  AMREX_ASSERT(m_memberDataDefined);
  BL_PROFILE("SootModel::addSootSourceTerm");
  if (m_sootVerbosity) {
    Print() << "SootModel::addSootSourceTerm(): Adding soot source term to "
            << vbox << std::endl;
  }
  const int nCompTr = dComp_lambda + 1; // Total number of components
  // Primitive components
  const int qRhoIndx = QRHO;
  const int qTempIndx = QTEMP;
  const int qSpecIndx = QFS;
  const int qSootIndx = QFSOOT;
  const int numSootMom = NUM_SOOT_MOMENTS;
  const int numSootVar = NUM_SOOT_VARS;
  const int rhoIndx = PeleC::Density;
  const int engIndx = PeleC::Eden;
  const int specIndx = PeleC::FirstSpec;
  // Number of relevant gas species
  const int numGS = GasSpecIndx::numGasSpecs;
  // Molar concentrations (mol/cm^3)
  // Ordered based on GasSpecIndx, not PeleC species
  Vector<Real> xi_n(numGS);
  // Vector of moment values M_xy (cm^(3(x + 2/3y))cm^(-3))
  // M00, M10, M01,..., N0
  Real moments[NUM_SOOT_VARS];
  /*
    These are the values inside the terms in fracMom
    momFV[NUM_SOOT_MOMENTS] - Weight of the delta function
    momFV[NUM_SOOT_MOMENTS+1] - modeCoef
    where modeCoef signifies the number of modes to be used
    If the moments are effectively zero, modeCoef = 0 and only 1 mode is used
    Otherwise, modeCoef = 1 and both modes are used
    The rest of the momFV values are used in fracMom fact1 = momFV[0],
    fact2 = momFV[1]^volOrd, fact2 = momFV[2]^surfOrd, etc.
  */
  Real momFV[NUM_SOOT_VARS + 1];
  // Vector of source terms for moment equations
  Real mom_src[NUM_SOOT_VARS];
  // Vector of reaction source terms for species
  Vector<Real> omega_src(numGS);
  const int nsr = m_numSurfReacts;
  // Forward and backward reaction rates
  Vector<Real> k_fwd(nsr);
  Vector<Real> k_bkwd(nsr);
  // Creation and destruction rates
  Vector<Real> w_fwd(nsr);
  Vector<Real> w_bkwd(nsr);

  // H2 absorbs the error from surface reactions
  const int absorbIndx = GasSpecIndx::indxH2;
  const int absorbIndxP = specIndx + local2GlobalIndx(absorbIndx);

  const auto lo = lbound(vbox);
  const auto hi = ubound(vbox);
  const auto Qstate = Qfab.array();
  const auto Ustate = Ufab.array();
  const auto coeff_state = coeff_cc.array();
  auto Sstate = Sfab.array();
  for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j ) {
      for (int i = lo.x; i <= hi.x; ++i) {
        Real Hi[NUM_SPECIES];
        const Real rho = Qstate(i,j,k,qRhoIndx);
        const Real T = Qstate(i,j,k,qTempIndx);
        // Dynamic viscosity
        const Real mu = coeff_state(i,j,k,dComp_mu);
        // Compute species enthalpy
        EOS::T2Hi(T, Hi);
        // Extract mass fractions for gas phases corresponding to GasSpecIndx
        // Compute the average molar mass (g/mol)
        Real molarMass = 0.;
        for (int sp = 0; sp != numGS; ++sp) {
          const int specConvIndx = local2GlobalIndx(sp);
          const int peleIndx = qSpecIndx + specConvIndx;
          Real cn = Qstate(i,j,k,peleIndx);
          enth_n[sp] = Hi[specConvIndx];
          Real conv = cn/m_gasMW[sp];
          xi_n[sp] = rho*conv;
          molarMass += conv;
          // Reset the reaction source term
          omega_src[sp] = 0.;
        }
        molarMass = 1./molarMass;
        // Molar concentration of the PAH inception species
        Real xi_PAH = xi_n[GasSpecIndx::indxPAH];
        // Extract moment values
        for (int mom = 0; mom != numSootVar; ++mom) {
          const int peleIndx = qSootIndx + mom;
          moments[mom] = Qstate(i,j,k,peleIndx);
          // Reset moment source
          mom_src[mom] = 0.;
        }
        // Convert moments from CGS to mol of C
        convertCGStoMol(moments);
        // Clip moments
        clipMoments(moments);
        // Compute constant values used throughout
        // (R*T*Pi/(2*A*rho_soot))^(1/2)
        const Real convT = std::sqrt(m_colFact*T);
        // Constant for free molecular collisions
        const Real colConst = convT*m_colFactPi23*m_colFact16*avogadros;
        // Collision frequency between two dimer in the free
        // molecular regime with van der Waals enhancement
        // Units: cm^3/mol-s
        const Real betaNucl = convT*m_betaNuclFact;
        const Real betaDimer = convT*m_betaDimerFact;
        // Compute the vector of factors used for moment interpolation
        computeFracMomVect(moments, momFV);
        // Estimate [DIMER]
        Real dimerConc = dimerization(T, convT, betaNucl, betaDimer,
                                      xi_PAH, momFV);
        // Add the nucleation source term to mom_src
        nucleationMomSrc(betaNucl, dimerConc, mom_src);
        // Add the condensation source term to mom_src
        condensationMomSrc(colConst, dimerConc, momFV, mom_src);
        // Add the coagulation source term to mom_src
        coagulationMomSrc(colConst, T, mu, rho, molarMass, momFV, mom_src);
        // Reaction rates for surface growth (k_sg), oxidation (k_ox),
        // and fragmentation (k_o2)
        Real k_sg = 0.;
        Real k_ox = 0.;
        Real k_o2 = 0.;
        // Compute the species reaction source terms into omega_src
        chemicalSrc(T, xi_n, moments, momFV, k_fwd, k_bkwd, w_fwd, w_bkwd,
                    k_sg, k_ox, k_o2, omega_src);
        // Compute the continuity and energy source term
        Real rho_src = 0.;
        Real eng_src = 0.;
        Real local_rate = 0.;
        for (int sp = 0; sp != numGS; ++sp) {
          // Convert to proper units
          omega_src[sp] *= m_gasMW[sp];
          // Convert from local gas species index to global gas species index
          const int specConvIndx = local2GlobalIndx(sp);
          const int peleIndx = specIndx + specConvIndx;
          Real spec_rate = std::abs(omega_src[sp])/
            (Ustate(i,j,k,peleIndx) + 1.E-12);
          Sstate(i,j,k,peleIndx) += omega_src[sp];
          local_rate = amrex::max(local_rate, spec_rate);
          rho_src += omega_src[sp];
          eng_src += omega_src[sp]*Hi[specConvIndx];
        }
        // Add the surface growth source to mom_src
        surfaceGrowthMomSrc(k_sg, momFV, mom_src);
        oxidFragMomSrc(k_ox, k_o2, momFV, mom_src);
        // Convert moment source terms back to CGS units
        convertMoltoCGS(mom_src, moments);
        if (m_conserveMass) {
          // Difference between mass lost from fluid and mass gained to soot
          Real del_rho_dot = rho_src + mom_src[1]*m_SootDensity;
          // Add that mass to H2
          Sstate(i,j,k,absorbIndxP) -= del_rho_dot;
          rho_src -= del_rho_dot;
          eng_src -= enth_n[absorbIndx]*del_rho_dot;
        }
        // Add density source term
        Real rho_rate = std::abs(rho_src)/rho;
        Sstate(i,j,k,rhoIndx) += rho_src;
        Real eng_rate = std::abs(eng_src/Ustate(i,j,k,engIndx));
        Sstate(i,j,k,engIndx) -= eng_src;
        local_rate = amrex::max(local_rate, amrex::max(eng_rate, rho_rate));
        // Add moment source terms
        for (int mom = 0; mom != numSootVar; ++mom) {
          const int peleIndx = sootIndx + mom;
          Sstate(i,j,k,peleIndx) += mom_src[mom];
          // Only limit based on positive moment sources
          local_rate = std::max(local_rate, -mom_src[mom]/moments[mom]);
        }
        if (dt*local_rate > m_maxDtRate) {
          local_soot_dt = std::min(local_soot_dt, m_maxDtRate/local_rate);
        }
      }
    }
  }
}

/********************************************************************
  Moment source terms
********************************************************************/

// Nucleation source term
void
SootModel::nucleationMomSrc(const Real&   betaNucl,
                            const Real&   dimerConc,
                            Vector<Real>& mom_src) const
{
  const Real dimerConc2 = dimerConc*dimerConc;
  for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
    mom_src[i] += 0.5*betaNucl*dimerConc2*m_momFact[i];
  }
  mom_src[NUM_SOOT_MOMENTS] += 0.5*betaNucl*dimerConc2;
}

// Condensation source term
void
SootModel::condensationMomSrc(const Real&         colConst,
                              const Real&         dimerConc,
                              const Vector<Real>& momFV,
                              Vector<Real>&       mom_src) const
{  
/** Compute condensation source values
    @param colConst Constant for free molecular collisions
    @param dimerConc Concentration of dimer
    @param momFV Vector of factors used in moment interpolation
    @param mom_src Moment source values
*/
  if (m_sootVerbosity >= 2) {
    Print() << "SootModel::condensationMomSrc(): Adding condensation"
            << std::endl;
  }
  Real weightDelta = momFV[NUM_SOOT_MOMENTS];
  for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
    const Real vv1 = MomOrderV[i] + 2.*m_SootAv;
    const Real vs1 = MomOrderS[i] + 2.*m_SootAs;
    const Real vv2 = MomOrderV[i] + m_SootAv;
    const Real vs2 = MomOrderS[i] + m_SootAs;
    Real volTerm = fracMom(vv1 - 1.,            vs1, momFV)*getDimerExp6(3)
      +  2.*fracMom(vv2 - 1.,                   vs2, momFV)*getDimerExp6(5)
      +     fracMom(MomOrderV[i] - 1., MomOrderS[i], momFV)*getDimerExp6(7)
      + 0.5*fracMom(vv1 - 2.,                   vs1, momFV)*getDimerExp6(9)
      +     fracMom(vv2 - 2.,                   vs2, momFV)*getDimerExp6(11)
      + 0.5*fracMom(MomOrderV[i] - 2., MomOrderS[i], momFV)*getDimerExp6(13);
    const Real ss3 = MomOrderS[i] + 3.*m_SootFitE;
    const Real sv3 = MomOrderV[i] - 2.*m_SootFitE;
    const Real ss2 = ss3 + m_SootAs;
    const Real sv2 = sv3 + m_SootAv;
    const Real ss1 = ss3 + 2.*m_SootAs;
    const Real sv1 = sv3 + 2.*m_SootAv;
    const Real surfTerm = fracMom(sv1 - 1., ss1, momFV)*getDimerExp6(3)
      +  2.*fracMom(sv2 - 1., ss2, momFV)*getDimerExp6(5)
      +     fracMom(sv3 - 1., ss3, momFV)*getDimerExp6(7)
      + 0.5*fracMom(sv1 - 2., ss1, momFV)*getDimerExp6(9)
      +     fracMom(sv2 - 2., ss2, momFV)*getDimerExp6(11)
      + 0.5*fracMom(sv3 - 2., ss3, momFV)*getDimerExp6(13);
    mom_src[i] += colConst*(MomOrderV[i]*volTerm + 
                            m_SootFitC*MomOrderS[i]*surfTerm)*dimerConc;
  }
  // Source for the weight of the delta function
  mom_src[NUM_SOOT_MOMENTS] -= 
    m_condFact*colConst*dimerConc*weightDelta;
}

// Compute the coagulation source term
void
SootModel::coagulationMomSrc(const Real&         colConst,
                             const Real&         T,
                             const Real&         mu,
                             const Real&         rho,
                             const Real&         molMass,
                             const Vector<Real>& momFV,
                             Vector<Real>&       mom_src) const
{
  if (m_sootVerbosity >= 2) {
    Print() << "SootModel::coagulationMomSrc(): Adding coagulation"
            << std::endl;
  }
  // Index of the weight of the delta function
  const int dwIndx = NUM_SOOT_MOMENTS;
  // Free molecular collision coefficient with van der Waals enhancements
  const Real C_fm = 2.2*colConst;
  // Continuum collision coefficient
  const Real C_cn = 8.*EOS::RU*T/(3.*mu);
  // Mean free path for finite nudsen number correction in continuum regimes
  const Real lambda = 3.*mu/rho*std::sqrt(M_PI*molMass/(8.*EOS::RU*T))*m_lambdaCoagFact;
  Real weightDelta2 = std::pow(momFV[dwIndx], 2);
  for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
    // Collisions between two first mode particles
    // Collision model: pure coalescence
    // S_(0+0) = (144*pi)^(1/3)*V0^(2/3)
    // Free molecular regime
    Real ss_fm = C_fm*m_ssfmCoagFact[i]*weightDelta2;
    // Continuum regime
    Real ss_cn = 4.*C_cn*(1. + 1.257*lambda*getNuclExp3(-1))*m_sscnCoagFact[i]*weightDelta2;
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
  Real ss_fm = -C_fm*weightDelta2*m_ssfmCoagFact[dwIndx];
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

// Surface growth source term
void
SootModel::surfaceGrowthMomSrc(const Real&         k_sg,
                               const Vector<Real>& momFV,
                               Vector<Real>&       mom_src) const
{
  if (m_sootVerbosity >= 2) {
    Print() << "SootModel::surfaceGrowthMomSrc(): Adding surface growth"
            << std::endl;
  }
  // Index of the weight of the delta function
  const int dwIndx = NUM_SOOT_MOMENTS;
  const Real weightDelta = momFV[dwIndx];
  const Real factor = m_SootDensityC*m_dVol*k_sg;
  for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
    Real fact1 = fracMom(MomOrderV[i] - 1., MomOrderS[i] + 1., momFV);
    Real fact2 = fracMom(MomOrderV[i] - 1. - 2.*m_SootFitE, 
                         MomOrderS[i] + 1. + 3.*m_SootFitE, momFV);
    mom_src[i] += factor*(MomOrderV[i]*fact1 + MomOrderS[i]*m_SootFitC*fact2);
  }
  // Weight of the delta function
  mom_src[dwIndx] -= m_nuclSurf*k_sg*m_SootDensityC*weightDelta;
}

// Oxidation and fragmentation source terms
void
SootModel::oxidFragMomSrc(const Real&         k_ox,
                          const Real&         k_o2,
                          const Vector<Real>& momFV,
                          Vector<Real>&       mom_src) const
{
  if (m_sootVerbosity >= 2) {
    Print() << "SootModel::oxidFragMomSrc(): Oxidation and fragmentation"
            << std::endl;
  }
  // Index of the weight of the delta function
  const int dwIndx = NUM_SOOT_MOMENTS;
  const Real weightDelta = momFV[dwIndx];
  const Real factOx = k_ox*m_dVol*m_SootDensityC;
  const Real factO2 = 2.*k_o2*m_dVol*m_SootDensityC;
  for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
    // Oxidation of the small particles
    Real small = -factOx*m_smallOxidFact[i]*weightDelta;
    // Oxidation of the larger particles
    Real factLarge = fracMomLarge(MomOrderV[i] - 1., MomOrderS[i] + 1., 
                                  momFV);
    Real large = -factOx*(MomOrderV[i] + m_two3rd*MomOrderS[i])*factLarge;
    // Add oxidation source
    mom_src[i] += large + small;
    // Add fragmentation source
    mom_src[i] += m_fragFact[i]*factO2*factLarge;
  }
  Real fracLarge = fracMomLarge(-1., 1., momFV);
  Real small = -factOx*m_smallOxidFact[dwIndx]*weightDelta;
  Real inter = m_nuclVol/(fracMomLarge(1., 0., momFV)/
                          fracMomLarge(0., 0., momFV));
  Real large = factOx*inter*fracLarge;
  // Add oxidation source for weight of delta function
  mom_src[dwIndx] += (small + large);
  // Add fragmentation source for weight of delta function
  mom_src[dwIndx] += inter*factO2*fracLarge;
}

// Return the dimer concentration
Real
SootModel::dimerization(const Real&         T,
                        const Real&         convT,
                        const Real&         betaNucl,
                        const Real&         betaDimer,
                        const Real&         xi_PAH,
                        const Vector<Real>& momFV) const
{
  // Collision coefficient for condensation
  const Real betaCond = getBetaCond(convT, momFV);
  // The effective rate of dimerization
  // This should coincide with the last reaction in the reaction lists
  const int fr = m_numSurfReacts - 1;
  const Real dimerRate =
    A_f[fr]*std::pow(T, n_f[fr])*std::exp(-ER_f[fr]/T)*xi_PAH*xi_PAH;

  // Using the following quadratic equation:
  // betaNucl*[DIMER]^2 + betaCond*[DIMER] - dimerRate = 0
  // compute the [DIMER] using the quadratic formula
  // x = -b + sqrt(b^2 - 4ac)/(2a)
  const Real delta = betaCond*betaCond + 4.*betaNucl*dimerRate;
  return (std::sqrt(delta) - betaCond)/(2.*betaNucl);
}

/*********************************************************************
  Moment interpolation functions
*********************************************************************/

// Compute the moment interpolation vector
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
void
SootModel::computeFracMomVect(const Vector<Real>& moments,
                              Vector<Real>&       momFV) const
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

  const Real minMom = std::min({M00, M10, M01, M20, M11, M02});
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
Real
SootModel::fracMomLarge(const Real          volOrd,
                        const Real          surfOrd,
                        const Vector<Real>& momFV) const
{
  // Weight of the delta function
  Real dwVal = momFV[NUM_SOOT_MOMENTS];
  Real factor = std::pow(m_nuclVol, volOrd)*std::pow(m_nuclSurf, surfOrd);
  // Remove the contribution from the first mode
  Real outMom = fracMom(volOrd, surfOrd, momFV) - dwVal*factor;
  
  // If the moment is negative, return a small (consistent) value
  if (outMom <= 0. || outMom != outMom) {
    return factor*1.E-66;
  }
  return outMom;
}

// Moment interpolation
Real
SootModel::fracMom(const Real          volOrd,
                   const Real          surfOrd,
                   const Vector<Real>& momFV) const
{
  // If modeCoef = 0.; only first mode is used
  // If modeCoef = 1.; both modes are used
  const Real modeCoef = momFV[NUM_SOOT_MOMENTS + 1];
  Real bothPFact = modeCoef*std::pow(m_nuclVol, volOrd)*
    std::pow(m_nuclSurf, surfOrd);
#if NUM_SOOT_MOMENTS == 3
  Real peak = std::pow(momFV[0], 1. - volOrd - surfOrd)*
    std::pow(momFV[1], volOrd)*std::pow(momFV[2], surfOrd);

  return momFV[3]*bothPFact + peak;
#elif NUM_SOOT_MOMENTS == 6
  Real prod = momFV[0];
  prod *= std::pow(momFV[1], volOrd);
  prod *= std::pow(momFV[2], surfOrd);
  prod *= std::pow(momFV[3], volOrd*volOrd);
  prod *= std::pow(momFV[4], volOrd*surfOrd);
  prod *= std::pow(momFV[5], surfOrd*surfOrd);

  return momFV[6]*bothPFact + prod;
#endif
}

// Interpolation for the reduced mass term (square root of sum) in the
// collision kernel for collision between a particle in each mode
// Only two grid functions used for all moments
// Limited sensitivity to increasing the number of grid functions
Real
SootModel::psiSL(const Real          x,
                 const Real          y,
                 const Real          a,
                 const Real          b,
                 const Vector<Real>& momFV) const
{
  const Real weightDelta = momFV[NUM_SOOT_MOMENTS];
  const Real factor = weightDelta*std::pow(m_nuclVol, a + m_two3rd*b);
  Real VF[3] = {2.*m_SootAv + x, m_SootAv + x, x};
  Real SF[3] = {2.*m_SootAs + y, m_SootAs + y, y};

  const Real FML_1 = fracMomLarge(VF[0] - 0.5, SF[0], momFV);
  const Real FML_2 = fracMomLarge(VF[1] - 0.5, SF[1], momFV);
  const Real FML_3 = fracMomLarge(VF[2] - 0.5, SF[2], momFV);
  // m_nuclVolExp6[i] = m_nuclVol^(2*i - 3)/6
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
Real
SootModel::psiLL(const Real          x,
                 const Real          y,
                 const Real          a,
                 const Real          b,
                 const Vector<Real>& momFV) const
{
  Real VF_xy[3] = {2.*m_SootAv + x, m_SootAv + x, x};
  Real SF_xy[3] = {2.*m_SootAs + y, m_SootAs + y, y};

  Real VF_ab[3] = {a, m_SootAv + a, 2.*m_SootAv + a};
  Real SF_ab[3] = {b, m_SootAs + b, 2.*m_SootAs + b};

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
Real
SootModel::FMCoagSL(const int           i,
                    const Vector<Real>& momFV) const
{
  // Weight of delta function N0 and M00
  if (i == NUM_SOOT_MOMENTS || i == 0) {
    return -psiSL(0., 0., 0., 0., momFV);
  }
  const Real fact1 = -2.*m_SootFitE;
  const Real fact2 = 3.*m_SootFitE;
  switch (i) {
  case 1: // M10
    return 0.;
  case 2: // M01
    return m_SootFitC*psiSL(fact1 - 1., fact2 + 1., 1., 0., momFV)
      - psiSL(0., 0., 0., 1., momFV);
  case 3: // M20
    return 2.*psiSL(1., 0., 1., 0., momFV);
  case 4: // M11
    return m_SootFitC*psiSL(fact1, fact2 + 1., 1., 0., momFV)
      + psiSL(0., 1., 1., 0., momFV)
      + m_SootFitC*psiSL(fact1 - 1., fact2 + 1., 2., 0., momFV)
      - psiSL(0., 0., 1., 1., momFV);
  case 5: // M02
    return 2.*m_SootFitC*psiSL(fact1 - 1., fact2 + 2., 1., 0., momFV)
      + m_SootFitC*m_SootFitC*psiSL(2.*fact1 - 2., -3.*fact1 + 2., 
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
Real
SootModel::FMCoagLL(const int           i,
                    const Vector<Real>& momFV) const
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
Real
SootModel::CNCoagSL(const int           i,
                    const Real&         lambda,
                    const Vector<Real>& momFV) const
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
        Real x = -2.*m_SootFitE - 1.;
        Real y = 3.*m_SootFitE + 1.;
        p1 = m_SootFitC*CNCoagSLFunc(n, x, y, lambda, momFV);
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
        Real x = -2.*m_SootFitE;
        Real y = 3.*m_SootFitE + 1.;
        p1 = m_SootFitC*CNCoagSLFunc(n, x, y, lambda, momFV);
      }
      {
        int n[] = {3, 4, 2, 1};
        p2 = CNCoagSLFunc(n, 0., 1., lambda, momFV);
      }
      {
        int n[] = {6, 7, 5, 4};
        Real x = -2.*m_SootFitE - 1.;
        Real y = 3.*m_SootFitE + 1.;
        p3 = m_SootFitC*CNCoagSLFunc(n, x, y, lambda, momFV);
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
        Real x = -2.*m_SootFitE - 1.;
        Real y = 3.*m_SootFitE + 2.;
        p1 = 2.*m_SootFitC*CNCoagSLFunc(n, x, y, lambda, momFV);
      }
      {
        int n[] = {6, 7, 5, 4};
        Real x = -4.*m_SootFitE - 2.;
        Real y = 6.*m_SootFitE + 2.;
        p2 = m_SootFitC*m_SootFitC*CNCoagSLFunc(n, x, y, lambda, momFV);
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

Real
SootModel::CNCoagSLFunc(int                 n[4],
                        const Real          x,
                        const Real          y,
                        const Real&         lambda,
                        const Vector<Real>& momFV) const
{
  Real xy_1 = fracMomLarge(x, y, momFV);
  Real xy_2 = fracMomLarge(x - m_SootAv, y - m_SootAs, momFV);
  Real xy_3 = fracMomLarge(x + m_SootAv, y + m_SootAs, momFV);
  Real xy_4 = fracMomLarge(x - 2.*m_SootAv, y - 2.*m_SootAs, momFV);
  Real n_1 = getNuclExp3(n[0]);
  Real n_2 = getNuclExp3(n[1]);
  Real n_3 = getNuclExp3(n[2]);
  Real n_4 = getNuclExp3(n[3]);
  return 2.*xy_1*n_1 + xy_2*n_2 + xy_3*n_3 + 
    1.257*lambda*(xy_1*n_3 + xy_2*n_1 + xy_3*n_4 + xy_4*n_2);
}

Real
SootModel::CNCoagLL(const int           i,
                    const Real&         lambda,
                    const Vector<Real>& momFV) const
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

Real
SootModel::CNCoagLLFunc(const Real          x,
                        const Real          y,
                        const Real&         lambda,
                        const Vector<Real>& momFV) const
{
  Real xy_1 = fracMomLarge(x, y, momFV);
  Real xy_2 = fracMomLarge(x - m_SootAv, y - m_SootAs, momFV);
  Real xy_3 = fracMomLarge(x + m_SootAv, y + m_SootAs, momFV);
  Real xy_4 = fracMomLarge(x - 2.*m_SootAv, y - 2.*m_SootAs, momFV);
  return 2.*xy_1*xy_1 + xy_2*xy_3 + xy_3*xy_2 
    + 1.257*lambda*(xy_1*xy_2 + xy_2*xy_1 + xy_3*xy_4 + xy_4*xy_3);
}

Real
SootModel::CNCoagLLFunc(const Real          x,
                        const Real          y,
                        const Real          a,
                        const Real          b,
                        const Real&         lambda,
                        const Vector<Real>& momFV) const
{
  Real xy_1 = fracMomLarge(x, y, momFV);
  Real xy_2 = fracMomLarge(x - m_SootAv, y - m_SootAs, momFV);
  Real xy_3 = fracMomLarge(x + m_SootAv, y + m_SootAs, momFV);
  Real xy_4 = fracMomLarge(x - 2.*m_SootAv, y - 2.*m_SootAs, momFV);
  Real ab_1 = fracMomLarge(a, b, momFV);
  Real ab_2 = fracMomLarge(a - m_SootAv, b - m_SootAs, momFV);
  Real ab_3 = fracMomLarge(a + m_SootAv, b + m_SootAs, momFV);
  Real ab_4 = fracMomLarge(a - 2.*m_SootAv, b - 2.*m_SootAs, momFV);
  return  2.*ab_1*xy_1 + ab_2*xy_3 + ab_3*xy_2 
    + 1.257*lambda*(ab_1*xy_2 + ab_2*xy_1 + ab_3*xy_4 + ab_4*xy_3);
}

Real
SootModel::getBetaCond(const Real&         convT,
                       const Vector<Real>& momFV) const
{
  // Collision frequency between two dimer in the free
  // molecular regime WITHOUT van der Waals enhancement
  // Units: 1/s
  const Real Cfm = m_colFactPi23*convT*m_colFact16*avogadros;
  const Real SN = fracMom(2.*m_SootAv,    2.*m_SootAs, momFV)*getDimerExp6(-3)
    +          2.*fracMom(   m_SootAv,       m_SootAs, momFV)*getDimerExp6(-1)
    +             fracMom(         0.,             0., momFV)*getDimerExp6(1)
    +         0.5*fracMom(2.*m_SootAv-1., 2.*m_SootAs, momFV)*getDimerExp6(3)
    +             fracMom(   m_SootAv-1.,    m_SootAs, momFV)*getDimerExp6(5)
    +         0.5*fracMom(           -1.,          0., momFV)*getDimerExp6(7);
  return Cfm*SN;
}
