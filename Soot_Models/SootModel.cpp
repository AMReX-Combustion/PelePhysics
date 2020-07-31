
// Standard library includes
#include <string>
#include <map>

// AMReX include statements
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_FArrayBox.H>

// PeleC include statements
#include "PeleC.H"

// PeleC-MP include statements
#include "SootModel.H"
#include "SootModel_derive.H"

using namespace amrex;

namespace SootConst
{
  AMREX_GPU_DEVICE_MANAGED Real V0 = 0.;
  AMREX_GPU_DEVICE_MANAGED Real S0 = 0.;
  AMREX_GPU_DEVICE_MANAGED Real SootDensityC = 0.;
  AMREX_GPU_DEVICE_MANAGED Real colFact = 0.;
  AMREX_GPU_DEVICE_MANAGED Real colFact23 = 0.;
  AMREX_GPU_DEVICE_MANAGED Real colFact16 = 0.;
  AMREX_GPU_DEVICE_MANAGED Real colFactPi23 = 0.;
  AMREX_GPU_DEVICE_MANAGED Real dimerVol = 0.;
  AMREX_GPU_DEVICE_MANAGED Real nuclVol = 0.;
  AMREX_GPU_DEVICE_MANAGED Real nuclSurf = 0.;
#if NUM_SOOT_MOMENTS == 3
  AMREX_GPU_DEVICE_MANAGED Real MomOrderV[3] = {0., 1., 0.};
  AMREX_GPU_DEVICE_MANAGED Real MomOrderS[3] = {0., 0., 1.};
  AMREX_GPU_DEVICE_MANAGED Real unitConv[4] = {0., 0., 0., 0.};
#elif NUM_SOOT_MOMENTS == 6
  AMREX_GPU_DEVICE_MANAGED Real MomOrderV[6] = {0., 1., 0., 2., 1., 0.};
  AMREX_GPU_DEVICE_MANAGED Real MomOrderS[6] = {0., 0., 1., 0., 1., 2.};
  AMREX_GPU_DEVICE_MANAGED Real unitConv[7] = {0., 0., 0., 0., 0., 0., 0.};
#endif
  AMREX_GPU_DEVICE_MANAGED int refIndx[NUM_SOOT_GS] = {-1};
};

// Default constructor
SootModel::SootModel()
  :
  m_sootVerbosity(0),
  m_readSootParams(false),
  m_memberDataDefined(false),
  m_conserveMass(true),
  m_PAHindx(-1),
  m_sootVarName(NUM_SOOT_MOMENTS + 1, ""),
  local_soot_dt(1.E18),
  m_maxDtRate(0.2),
  m_reactDataFilled(false),
  m_gasSpecNames(NUM_SOOT_GS, ""),
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
  const int ngs = NUM_SOOT_GS;
  // This should be called after readSootParams()
  AMREX_ASSERT(m_readSootParams);
  // Relevant species names for the surface chemistry model
  // TODO: Currently must correspond to GasSpecIndx enum in SootModel.H
  m_gasSpecNames = {"H2", "H", "OH", "H2O", "CO", "C2H2", "O2", m_PAHname};
  // Number of accounted for PAH particles (hard-coded)
  const int numPAH = 3;
  // Names of PAH species
  // TODO: Currently can only handle naphthalene(C10H8), phenathrene(C14H10),
  // and pyrene(C16H10) and can only handle 1 PAH inception species
  std::string PAH_names[] = {"A2", "A3", "A4"};
  // Corresponding sticking coefficients
  const Real PAH_gammas[] = {0.002, 0.015, 0.025};
  // Average number of C atoms per PAH species
  const Real PAH_numC[] = {10., 14., 16.};

  // Determine which PAH inception species is being used
  bool corrPAH = false;
  for (int cpah = 0; cpah != numPAH; ++cpah) {
    if (PAH_names[cpah] == m_PAHname) {
      m_gammaStick = PAH_gammas[cpah];
      SootConst::dimerVol = 2.*PAH_numC[cpah];
      corrPAH = true;
    }
  }
  SootConst::nuclVol = 2.*SootConst::dimerVol;
  SootConst::nuclSurf = std::pow(SootConst::nuclVol, 2./3.);
  if (!corrPAH) {
    Abort(m_PAHname + " not recognized as PAH inception species. Must be" +
          PAH_names[0] + ", " + PAH_names[1] + ", or " + PAH_names[2]);
  }
  Vector<std::string> spec_names(NUM_SPECIES);
  EOS::speciesNames(spec_names);
  // Loop over all species
  for (int i = 0; i != NUM_SPECIES; ++i) {
    // Check if species is the PAH inceptor
    if (spec_names[i] == m_PAHname) m_PAHindx = i;
    // Check if species matches others for surface reactions
    for (int sootSpec = 0; sootSpec != ngs; ++sootSpec) {
      if (spec_names[i] == m_gasSpecNames[sootSpec])
        SootConst::refIndx[sootSpec] = i;
    }
  }
  // Return error if no PAH inception species is specified in mechanism
  if (m_PAHindx == -1) {
    Abort("PAH inception species was not found in PelePhysics mechanism");
  }
  // Return error if not all soot species are present
  for (int sootSpec = 0; sootSpec != ngs; ++sootSpec) {
    if (SootConst::refIndx[sootSpec] == -1)
      Abort("Species " + m_gasSpecNames[sootSpec] +
            " must be present in PelePhysics mechanism");
  }

  // Assign moment factor member data
  defineMemberData();
  // From SootModel_react.cpp
  // Initialize reaction and species member data
  initializeReactData();

  if (m_sootVerbosity && ParallelDescriptor::IOProcessor()) {
    Print() << "SootModel::define(): Soot model successfully defined"
            << std::endl;
  }
}

// Read soot related inputs
void
SootModel::readSootParams()
{
  ParmParse pp("soot");
  pp.get("incept_pah", m_PAHname);
  m_sootVerbosity = 0;
  pp.query("v", m_sootVerbosity);
  if (m_sootVerbosity) {
    Print() << "SootModel::readSootParams(): Reading input parameters"
            << std::endl;
  }
  // Set the maximum allowable change for variables during soot source terms
  m_maxDtRate = 0.2;
  pp.query("max_dt_rate", m_maxDtRate);
  // Determines if mass is conserved by adding lost mass to H2
  m_conserveMass = false;
  pp.query("conserve_mass", m_conserveMass);
  m_readSootParams = true;
}

// Define all member data
void
SootModel::defineMemberData()
{
  if (m_sootVerbosity) {
    Print() << "SootModel::defineMemberData(): Defining member data"
            << std::endl;
  }

  SootConst::V0 = SootConst::SootMolarMass/(SootConst::avogadros*SootConst::SootDensity); // cm^3
  SootConst::S0 = std::pow(36.*M_PI, 1./3.)*std::pow(SootConst::V0, 2./3.); // cm^2
  SootConst::SootDensityC = SootConst::SootChi*SootConst::S0;
  SootConst::colFact = M_PI*EOS::RU/(2.*SootConst::avogadros*SootConst::SootDensity);
  SootConst::colFact23 = std::pow(SootConst::V0, 2./3.);
  SootConst::colFact16 = std::pow(SootConst::V0, 1./6.); // Units: cm^0.5
  SootConst::colFactPi23 = std::pow(6./M_PI, 2./3.);
  Gpu::HostVector<Real> dimerExp6(9);
  Gpu::HostVector<Real> nuclVolExp3(11);
  Gpu::HostVector<Real> nuclVolExp6(11);
  Gpu::HostVector<Real> momFact(NUM_SOOT_MOMENTS);
  Gpu::HostVector<Real> ssfmCoagFact(NUM_SOOT_VARS);
  Gpu::HostVector<Real> sscnCoagFact(NUM_SOOT_MOMENTS);
  Gpu::HostVector<Real> smallOxidFact(NUM_SOOT_VARS);
  Gpu::HostVector<Real> fragFact(NUM_SOOT_MOMENTS);
  // Compute V_nucl and V_dimer to fractional powers
  for (int i = 0; i != 9; ++i) {
    Real exponent = 2.*(Real)i - 3.;
    dimerExp6[i] = std::pow(SootConst::dimerVol, exponent/6.);
  }
  for (int i = 0; i != 11; ++i) {
    Real exponent = (Real)i - 3.;
    nuclVolExp3[i] = std::pow(SootConst::nuclVol, exponent/3.);
    nuclVolExp6[i] = std::pow(SootConst::nuclVol, exponent/6.);
  }

  for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
    // Used to convert moments to mol of C
    SootConst::unitConv[i] = std::pow(SootConst::V0, SootConst::MomOrderV[i])*
      std::pow(SootConst::S0, SootConst::MomOrderS[i])*SootConst::avogadros;
    // Used for computing nucleation source term
    momFact[i] = std::pow(SootConst::nuclVol, SootConst::MomOrderV[i])*
      std::pow(SootConst::nuclSurf, SootConst::MomOrderS[i]);
  }
  // and to convert the weight of the delta function
  SootConst::unitConv[NUM_SOOT_MOMENTS] = SootConst::avogadros;

  // Coagulation, oxidation, and fragmentation factors
  m_lambdaCoagFact =
    1./(std::pow(6.*SootConst::SootMolarMass/(M_PI*SootConst::SootDensity*SootConst::avogadros), 1./3.));
  for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
    const Real expFact = SootConst::MomOrderV[i] + 2./3.*SootConst::MomOrderS[i];
    const Real factor = (std::pow(2., expFact - 1.) - 1.);
    ssfmCoagFact[i] =
      std::pow(2., 2.5)*factor*std::pow(SootConst::nuclVol, expFact + 1./6.);
    sscnCoagFact[i] = factor*std::pow(SootConst::nuclVol, expFact);
    smallOxidFact[i] = std::pow(SootConst::nuclVol, expFact - 1./3.);
    fragFact[i] = std::pow(2., 1. - SootConst::MomOrderV[i] - SootConst::MomOrderS[i]) - 1.;
  }
  const Real ne32 = nuclVolExp3[2 + 3];
  const Real ne3m = nuclVolExp3[-1 + 3];
  ssfmCoagFact[NUM_SOOT_MOMENTS] = std::pow(2., 2.5)*
    std::pow(1./SootConst::nuclVol, 0.5)*ne32;
  smallOxidFact[NUM_SOOT_MOMENTS] = ne3m;
  // Build the vectors in the data container
  m_SootDataContainer.build_vectors(dimerExp6, nuclVolExp3, nuclVolExp6,
                                    momFact, ssfmCoagFact, sscnCoagFact,
                                    smallOxidFact, fragFact);
  // Beta and dimer factors
  // cm^0.5 mol^-1
  const Real dnfact = 4.*std::sqrt(2.)*SootConst::colFact16*
    SootConst::colFactPi23*SootConst::avogadros;
  m_betaDimerFact = dnfact*std::pow(0.5*SootConst::dimerVol, 1./6.);
  m_betaNuclFact = 2.2*dnfact*dimerExp6[2];

  /// Condensation factor
  m_condFact = std::sqrt((1./SootConst::nuclVol) + (1./SootConst::dimerVol))*
    std::pow((std::pow(SootConst::nuclVol, 1./3.) + std::pow(SootConst::dimerVol, 1./3.)), 2.);
  m_SootDataContainer.build_scalars(m_betaDimerFact, m_betaNuclFact, m_condFact, m_lambdaCoagFact);
  m_memberDataDefined = true;
}

// Add derive plot variables
void
SootModel::addSootDerivePlotVars(DeriveList&           derive_lst,
                                 const DescriptorList& desc_lst)
{
  // Add in soot variables
  Vector<std::string> sootNames = {"rho_soot", "soot_dp", "soot_np"};
  derive_lst.add("soot_vars",IndexType::TheCellType(),sootNames.size(),
                 sootNames,soot_genvars, the_same_box);
  derive_lst.addComponent("soot_vars",desc_lst,0,PeleC::Density, NVAR);

  // Variables associated with the second mode (large particles)
  Vector<std::string> large_part_names = {"NL", "VL", "SL"};
  derive_lst.add("soot_large_particles",IndexType::TheCellType(),
                 large_part_names.size(),large_part_names,
                 soot_largeparticledata, the_same_box);
  derive_lst.addComponent("soot_large_particles",desc_lst,0,PeleC::Density, NVAR);
}

// Add soot source term
void
SootModel::addSootSourceTerm(const Box&                vbox,
                             Array4<const Real> const& Ustate,
                             Array4<const Real> const& Qstate,
                             Array4<const Real> const& coeff_state,
                             Array4<Real> const&       soot_state,
                             const Real                time,
                             const Real                dt) const
{
  AMREX_ASSERT(m_memberDataDefined);
  BL_PROFILE("SootModel::addSootSourceTerm");
  if (m_sootVerbosity && ParallelDescriptor::IOProcessor()) {
    Print() << "SootModel::addSootSourceTerm(): Adding soot source term to "
            << vbox << std::endl;
  }
  // Primitive components
  const int qRhoIndx = QRHO;
  const int qTempIndx = QTEMP;
  const int qSpecIndx = QFS;
  const int qSootIndx = QFSOOT;
  const int numSootVar = NUM_SOOT_VARS;
  const int rhoIndx = URHO;
  const int engIndx = UEDEN;
  const int specIndx = UFS;
  const int sootIndx = UFSOOT;
  const Real betaDF = m_betaDimerFact;
  const Real betaNF = m_betaNuclFact;

  const bool conserveMass = m_conserveMass;
  Real mw_fluid[NUM_SPECIES];
  EOS::molecular_weight(mw_fluid);

  // H2 absorbs the error from surface reactions
  const int absorbIndx = SootConst::GasSpecIndx::indxH2;
  const int absorbIndxP = specIndx + SootConst::refIndx[absorbIndx];

  const Real maxDtRate = m_maxDtRate;
  SootData sd = m_SootDataContainer.getSootData();
  SootReactions sr = m_SootReactionContainer.getSootReactions();
  amrex::ParallelFor(vbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      Real Hi[NUM_SPECIES];
      Real enth_n[NUM_SOOT_GS];
      Real omega_src[NUM_SOOT_GS];
      // Molar concentrations (mol/cm^3)
      Real xi_n[NUM_SOOT_GS];
      // Array of moment values M_xy (cm^(3(x + 2/3y))cm^(-3))
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
      Real momFV[NUM_SOOT_VARS+1];
      // Array of source terms for moment equations
      Real mom_src[NUM_SOOT_VARS];
      // Forward and backward reaction rates
      Real k_fwd[NUM_SOOT_REACT];
      Real k_bkwd[NUM_SOOT_REACT];
      Real w_fwd[NUM_SOOT_REACT];
      Real w_bkwd[NUM_SOOT_REACT];
      const Real rho = Qstate(i,j,k,qRhoIndx);
      const Real T = Qstate(i,j,k,qTempIndx);
      AMREX_ASSERT(std::abs(T) < 1.E20);
      // Dynamic viscosity
      const Real mu = coeff_state(i,j,k,dComp_mu);
      // Compute species enthalpy
      EOS::T2Hi(T, Hi);
      // Extract mass fractions for gas phases corresponding to GasSpecIndx
      // Compute the average molar mass (g/mol)
      Real molarMass = 0.;
      for (int sp = 0; sp != NUM_SOOT_GS; ++sp) {
        const int specConvIndx = SootConst::refIndx[sp];
        const int peleIndx = qSpecIndx + specConvIndx;
        Real cn = Qstate(i,j,k,peleIndx);
        enth_n[sp] = Hi[specConvIndx];
        Real conv = cn/mw_fluid[specConvIndx];
        xi_n[sp] = rho*conv;
        molarMass += conv;
        // Reset the reaction source term
        omega_src[sp] = 0.;
      }
      molarMass = 1./molarMass;
      // Molar concentration of the PAH inception species
      Real xi_PAH = xi_n[SootConst::GasSpecIndx::indxPAH];
      // Extract moment values
      for (int mom = 0; mom != numSootVar; ++mom) {
        const int peleIndx = qSootIndx + mom;
        moments[mom] = Qstate(i,j,k,peleIndx);
        // Reset moment source
        mom_src[mom] = 0.;
      }
      // Convert moments from CGS to mol of C
      SootConst::convertCGStoMol(moments);
      // Clip moments
      sd.clipMoments(moments);
      // Compute constant values used throughout
      // (R*T*Pi/(2*A*rho_soot))^(1/2)
      const Real convT = std::sqrt(SootConst::colFact*T);
      // Constant for free molecular collisions
      const Real colConst = convT*SootConst::colFactPi23*
        SootConst::colFact16*SootConst::avogadros;
      // Collision frequency between two dimer in the free
      // molecular regime with van der Waals enhancement
      // Units: cm^3/mol-s
      const Real betaNucl = convT*betaNF;
      const Real betaDimer = convT*betaDF;
      // Compute the vector of factors used for moment interpolation
      sd.computeFracMomVect(moments, momFV);
      // Compute the dimerization rate
      const Real dimerRate = sr.dimerRate(T, xi_PAH);
      // Estimate [DIMER]
      Real dimerConc = sd.dimerization(T, convT, betaNucl, betaDimer,
                                       xi_PAH, dimerRate, momFV);
      // Add the nucleation source term to mom_src
      sd.nucleationMomSrc(betaNucl, dimerConc, mom_src);
      // Add the condensation source term to mom_src
      sd.condensationMomSrc(colConst, dimerConc, momFV, mom_src);
      // Add the coagulation source term to mom_src
      sd.coagulationMomSrc(colConst, T, mu, rho, molarMass, momFV, mom_src);
      Real surf = SootConst::S0*sd.fracMom(0., 1., momFV);
      // Reaction rates for surface growth (k_sg), oxidation (k_ox),
      // and fragmentation (k_o2)
      Real k_sg = 0.;
      Real k_ox = 0.;
      Real k_o2 = 0.;
      // Compute the species reaction source terms into omega_src
      sr.chemicalSrc(T, surf, xi_n, moments, momFV, k_fwd, k_bkwd, w_fwd, w_bkwd,
                     k_sg, k_ox, k_o2, omega_src);
      // Compute the continuity and energy source term
      Real rho_src = 0.;
      Real eng_src = 0.;
      for (int sp = 0; sp != NUM_SOOT_GS; ++sp) {
        // Convert from local gas species index to global gas species index
        const int specConvIndx = SootConst::refIndx[sp];
        // Convert to proper units
        omega_src[sp] *= mw_fluid[specConvIndx];
        const int peleIndx = specIndx + specConvIndx;
        soot_state(i,j,k,peleIndx) += omega_src[sp];
        rho_src += omega_src[sp];
        eng_src += omega_src[sp]*Hi[specConvIndx];
      }
      // Add the surface growth source to mom_src
      sd.surfaceGrowthMomSrc(k_sg, momFV, mom_src);
      sd.oxidFragMomSrc(k_ox, k_o2, momFV, mom_src);
      // Convert moment source terms back to CGS units
      SootConst::convertMoltoCGS(mom_src, moments);
      if (conserveMass) {
        // Difference between mass lost from fluid and mass gained to soot
        Real del_rho_dot = rho_src + mom_src[1]*SootConst::SootDensity;
        // Add that mass to H2
        soot_state(i,j,k,absorbIndxP) -= del_rho_dot;
        rho_src -= del_rho_dot;
        eng_src -= enth_n[absorbIndx]*del_rho_dot;
      }
      // Add density source term
      soot_state(i,j,k,rhoIndx) += rho_src;
      soot_state(i,j,k,engIndx) -= eng_src;
      // Add moment source terms
      for (int mom = 0; mom != numSootVar; ++mom) {
        const int peleIndx = sootIndx + mom;
        soot_state(i,j,k,peleIndx) += mom_src[mom];
      }
    });
}
