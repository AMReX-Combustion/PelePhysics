
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

// Pele include statements
#ifdef SOOT_PELE_LM
#include "PeleLM.H"
#else
#include "PeleC.H"
#endif

// PeleC-MP include statements
#include "SootModel.H"
#include "SootModel_derive.H"

static amrex::Box
the_same_box(const amrex::Box& b)
{
  return b;
}

// Default constructor
SootModel::SootModel()
  : m_sootVerbosity(0),
    m_sootIndx(),
    m_setIndx(false),
    m_readSootParams(false),
    m_memberDataDefined(false),
    m_conserveMass(true),
    m_PAHindx(-1),
    m_sootVarName(NUM_SOOT_MOMENTS + 1, ""),
    local_soot_dt(1.E18),
    m_maxDtRate(-1.),
    m_reactDataFilled(false),
    m_gasSpecNames(NUM_SOOT_GS, "")
{
  m_sootData = std::make_unique<SootData>();
  m_sootReact = std::make_unique<SootReaction>();
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
  // Double check indices are set
  m_setIndx = m_sootIndx.checkIndices();
  if (!m_setIndx)
    Abort("SootModel::define(): Must set indices before defining");
  // Relevant species names for the surface chemistry model
  // TODO: Currently must correspond to GasSpecIndx enum in Constants_Soot.H
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
  Real dimerVol;
  for (int cpah = 0; cpah < numPAH; ++cpah) {
    if (PAH_names[cpah] == m_PAHname) {
      m_gammaStick = PAH_gammas[cpah];
      dimerVol = 2. * PAH_numC[cpah];
      corrPAH = true;
    }
  }

  if (!corrPAH) {
    Abort(
      m_PAHname + " not recognized as PAH inception species. Must be " +
      PAH_names[0] + ", " + PAH_names[1] + ", or " + PAH_names[2]);
  }
  Vector<std::string> spec_names(NUM_SPECIES);
  pele::physics::eos::speciesNames(spec_names);
  // Loop over all species
  for (int i = 0; i < NUM_SPECIES; ++i) {
    // Check if species is the PAH inceptor
    if (spec_names[i] == m_PAHname)
      m_PAHindx = i;
    // Check if species matches others for surface reactions
    for (int sootSpec = 0; sootSpec < ngs; ++sootSpec) {
      if (spec_names[i] == m_gasSpecNames[sootSpec]) {
        m_sootData->refIndx[sootSpec] = i;
      }
    }
  }
  // Return error if no PAH inception species is specified in mechanism
  if (m_PAHindx == -1) {
    Abort("PAH inception species was not found in PelePhysics mechanism");
  }
  // Return error if not all soot species are present
  for (int sootSpec = 0; sootSpec < ngs; ++sootSpec) {
    if (m_sootData->refIndx[sootSpec] == -1) {
      Abort(
        "Species " + m_gasSpecNames[sootSpec] +
        " must be present in PelePhysics mechanism");
    }
  }

  // Assign moment factor member data
  defineMemberData(dimerVol);
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
  m_maxDtRate = 0.05;
  pp.query("max_dt_rate", m_maxDtRate);
  // Determines if mass is conserved by adding lost mass to H2
  m_conserveMass = false;
  pp.query("conserve_mass", m_conserveMass);
  m_readSootParams = true;
}

// Define all member data
void
SootModel::defineMemberData(const Real dimerVol)
{
  if (m_sootVerbosity) {
    Print() << "SootModel::defineMemberData(): Defining member data"
            << std::endl;
  }
  SootConst sc;
  Real nuclVol = 2. * dimerVol;
  Real nuclSurf = std::pow(nuclVol, 2. / 3.);
  m_sootData->nuclVol = nuclVol;
  m_sootData->nuclSurf = nuclSurf;
  // Compute V_nucl and V_dimer to fractional powers
  for (int i = 0; i < 9; ++i) {
    Real exponent = 2. * (Real)i - 3.;
    m_sootData->dime6[i] = std::pow(dimerVol, exponent / 6.);
  }
  for (int i = 0; i < 11; ++i) {
    Real exponent = (Real)i - 3.;
    m_sootData->nve3[i] = std::pow(nuclVol, exponent / 3.);
    m_sootData->nve6[i] = std::pow(nuclVol, exponent / 6.);
  }

  for (int i = 0; i < NUM_SOOT_MOMENTS; ++i) {
    // Used to convert moments to mol of C
    m_sootData->unitConv[i] = std::pow(sc.V0, sc.MomOrderV[i]) *
                              std::pow(sc.S0, sc.MomOrderS[i]) *
                              pele::physics::Constants::Avna;
    // First convert moments from SI to CGS, only necessary for PeleLM
    const Real expFact = 3. * sc.MomOrderV[i] + 2. * sc.MomOrderS[i];
    m_sootData->unitConv[i] *= std::pow(sc.len_conv, 3. - expFact);
    // Used for computing nucleation source term
    m_sootData->momFact[i] =
      std::pow(nuclVol, sc.MomOrderV[i]) * std::pow(nuclSurf, sc.MomOrderS[i]);
  }
  // and to convert the weight of the delta function
  m_sootData->unitConv[NUM_SOOT_MOMENTS] =
    pele::physics::Constants::Avna * std::pow(sc.len_conv, 3);

  // Coagulation, oxidation, and fragmentation factors
  m_sootData->lambdaCF =
    1. / (std::pow(
           6. * sc.SootMolarMass /
             (M_PI * sc.SootDensity * pele::physics::Constants::Avna),
           1. / 3.));
  for (int i = 0; i < NUM_SOOT_MOMENTS; ++i) {
    const Real expFact = sc.MomOrderV[i] + 2. / 3. * sc.MomOrderS[i];
    const Real factor = (std::pow(2., expFact - 1.) - 1.);
    m_sootData->ssfmCF[i] =
      std::pow(2., 2.5) * factor * std::pow(nuclVol, expFact + 1. / 6.);
    m_sootData->sscnCF[i] = factor * std::pow(nuclVol, expFact);
    m_sootData->smallOF[i] = std::pow(nuclVol, expFact - 1. / 3.);
    m_sootData->fragFact[i] =
      std::pow(2., 1. - sc.MomOrderV[i] - sc.MomOrderS[i]) - 1.;
  }
  const Real ne32 = m_sootData->nve3[2 + 3];
  const Real ne3m = m_sootData->nve3[-1 + 3];
  m_sootData->ssfmCF[NUM_SOOT_MOMENTS] =
    std::pow(2., 2.5) * std::pow(1. / nuclVol, 0.5) * ne32;
  m_sootData->smallOF[NUM_SOOT_MOMENTS] = ne3m;

  // Beta and dimer factors
  // cm^0.5 mol^-1
  const Real dnfact = 4. * std::sqrt(2.) * sc.colFact16 * sc.colFactPi23 *
                      pele::physics::Constants::Avna;
  m_betaDimerFact = dnfact * std::pow(0.5 * dimerVol, 1. / 6.);
  m_betaNuclFact = 2.2 * dnfact * m_sootData->dime6[2];

  /// Condensation factor
  m_sootData->condFact =
    std::sqrt((1. / nuclVol) + (1. / dimerVol)) *
    std::pow((std::pow(nuclVol, 1. / 3.) + std::pow(dimerVol, 1. / 3.)), 2.);
  m_memberDataDefined = true;
}

// Add derive plot variables
void
SootModel::addSootDerivePlotVars(
  DeriveList& derive_lst,
  const DescriptorList& desc_lst,
  const int rhoIndx,
  const int sootIndx)
{
  // Double check indices are set
  m_setIndx = m_sootIndx.checkIndices();
  if (!m_setIndx)
    Abort("SootModel::addSootDerivePlotVars(): Must set indices first");
  // Add in soot variables
  Vector<std::string> sootNames = {"rho_soot", "sum_rho_soot"};
  derive_lst.add(
    "soot_vars", IndexType::TheCellType(), sootNames.size(), sootNames,
    soot_genvars, the_same_box);
  derive_lst.addComponent("soot_vars", desc_lst, State_Type, rhoIndx, 1);
  derive_lst.addComponent(
    "soot_vars", desc_lst, State_Type, sootIndx, NUM_SOOT_MOMENTS + 1);

  // Variables associated with the second mode (large particles)
  Vector<std::string> large_part_names = {"NL", "rho_soot_L", "soot_S_L"};
  derive_lst.add(
    "soot_large_particles", IndexType::TheCellType(), large_part_names.size(),
    large_part_names, soot_largeparticledata, the_same_box);
  derive_lst.addComponent(
    "soot_large_particles", desc_lst, State_Type, sootIndx,
    NUM_SOOT_MOMENTS + 1);
}

// Add soot source term
void
SootModel::addSootSourceTerm(
  const Box& vbox,
  Array4<const Real> const& Qstate,
  Array4<const Real> const& coeff_mu,
  Array4<Real> const& soot_state,
  const Real time,
  const Real dt) const
{
  AMREX_ASSERT(m_memberDataDefined);
  AMREX_ASSERT(m_setIndx);
  BL_PROFILE("SootModel::addSootSourceTerm");
  if (m_sootVerbosity && ParallelDescriptor::IOProcessor()) {
    Print() << "SootModel::addSootSourceTerm(): Adding soot source term to "
            << vbox << std::endl;
  }
  // Primitive components
  const int qRhoIndx = m_sootIndx.qRhoIndx;
  const int qTempIndx = m_sootIndx.qTempIndx;
  const int qSpecIndx = m_sootIndx.qSpecIndx;
  const int qSootIndx = m_sootIndx.qSootIndx;
  const int rhoIndx = m_sootIndx.rhoIndx;
  const int engIndx = m_sootIndx.engIndx;
  const int specIndx = m_sootIndx.specIndx;
  const int sootIndx = m_sootIndx.sootIndx;
  const Real betaNF = m_betaNuclFact;

  const bool conserveMass = m_conserveMass;

  // H2 absorbs the error from surface reactions
  const int absorbIndx = SootGasSpecIndx::indxH2;
  const int absorbIndxN = m_sootData->refIndx[absorbIndx];
  const int absorbIndxP = specIndx + absorbIndxN;

  SootData const* sd = m_sootData.get();
  SootReaction const* sr = m_sootReact.get();
  amrex::ParallelFor(vbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    auto eos = pele::physics::PhysicsType::eos();
    SootConst sc;
    GpuArray<Real, NUM_SPECIES> mw_fluid;
    eos.molecular_weight(mw_fluid.data());
    GpuArray<Real, NUM_SPECIES> Hi;
    GpuArray<Real, NUM_SOOT_GS> omega_src;
    // Molar concentrations (mol/cm^3)
    GpuArray<Real, NUM_SOOT_GS> xi_n;
    // Array of moment values M_xy (cm^(3(x + 2/3y))cm^(-3))
    // M00, M10, M01,..., N0
    GpuArray<Real, NUM_SOOT_MOMENTS + 1> moments;
    Real* momentsPtr = moments.data();
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
    GpuArray<Real, NUM_SOOT_MOMENTS + 2> mom_fv;
    Real* mom_fvPtr = mom_fv.data();
    // Array of source terms for moment equations
    GpuArray<Real, NUM_SOOT_MOMENTS + 1> mom_src;
    Real* mom_srcPtr = mom_src.data();
    const Real rho = Qstate(i, j, k, qRhoIndx) * sc.rho_conv;
    const Real T = Qstate(i, j, k, qTempIndx);
    // Dynamic viscosity
    const Real mu = coeff_mu(i, j, k) * sc.mu_conv;
    // Compute species enthalpy
    eos.T2Hi(T, Hi.data());
    // Extract mass fractions for gas phases corresponding to GasSpecIndx
    // Compute the average molar mass (g/mol)
    Real molarMass = 0.;
    for (int sp = 0; sp < NUM_SOOT_GS; ++sp) {
      const int specConvIndx = sd->refIndx[sp];
      const int peleIndx = qSpecIndx + specConvIndx;
      // State provided by PeleLM is the concentration, rhoY
#ifdef SOOT_PELE_LM
      Real cn = Qstate(i, j, k, peleIndx) / Qstate(i, j, k, qRhoIndx);
#else
      Real cn = Qstate(i, j, k, peleIndx);
#endif
      if (cn < 2.E-10)
        cn = 0.;
      Real conv = cn / mw_fluid[specConvIndx];
      xi_n[sp] = rho * conv;
      molarMass += conv;
      // Reset the reaction source term
      omega_src[sp] = 0.;
    }
    molarMass = 1. / molarMass;
    // Molar concentration of the PAH inception species
    Real xi_PAH = xi_n[SootGasSpecIndx::indxPAH];
    // Extract moment values
    for (int mom = 0; mom < NUM_SOOT_MOMENTS + 1; ++mom) {
      const int peleIndx = qSootIndx + mom;
      moments[mom] = Qstate(i, j, k, peleIndx);
      // Reset moment source
      mom_src[mom] = 0.;
    }
    // Convert moments from CGS to mol of C
    sd->convertToMol(moments.data());
    // Clip moments
    sd->clipMoments(moments.data());
    // Compute constant values used throughout
    // (R*T*Pi/(2*A*rho_soot))^(1/2)
    const Real convT = std::sqrt(sc.colFact * T);
    // Constant for free molecular collisions
    const Real colConst =
      convT * sc.colFactPi23 * sc.colFact16 * pele::physics::Constants::Avna;
    // Collision frequency between two dimer in the free
    // molecular regime with van der Waals enhancement
    // Units: cm^3/mol-s
    const Real betaNucl = convT * betaNF;
    // Compute the vector of factors used for moment interpolation
    sd->computeFracMomVect(moments.data(), mom_fvPtr);
    // Compute the dimerization rate
    const Real dimerRate = sr->dimerRate(T, xi_PAH);
    // Estimate [DIMER]
    Real dimerConc = sd->dimerization(convT, betaNucl, dimerRate, mom_fvPtr);
    // Add the nucleation source term to mom_src
    sd->nucleationMomSrc(betaNucl, dimerConc, mom_srcPtr);
    // Add the condensation source term to mom_src
    sd->condensationMomSrc(colConst, dimerConc, mom_fvPtr, mom_srcPtr);
    // Add the coagulation source term to mom_src
    sd->coagulationMomSrc(
      colConst, T, mu, rho, molarMass, mom_fvPtr, mom_srcPtr);
    Real surf = sc.S0 * sd->fracMom(0., 1., mom_fvPtr);
    // Reaction rates for surface growth (k_sg), oxidation (k_ox),
    // and fragmentation (k_o2)
    Real k_sg = 0.;
    Real k_ox = 0.;
    Real k_o2 = 0.;
    // Compute the species reaction source terms into omega_src
    sr->chemicalSrc(
      T, surf, xi_n.data(), moments.data(), k_sg, k_ox, k_o2, omega_src.data());
    // Compute the continuity and energy source term
    Real rho_src = 0.;
    Real eng_src = 0.;
    for (int sp = 0; sp < NUM_SOOT_GS; ++sp) {
      // Convert from local gas species index to global gas species index
      const int specConvIndx = sd->refIndx[sp];
      // Convert to proper units
      omega_src[sp] *= mw_fluid[specConvIndx];
      const int peleIndx = specIndx + specConvIndx;
      soot_state(i, j, k, peleIndx) += omega_src[sp] * sc.mass_src_conv;
      rho_src += omega_src[sp];
      eng_src += omega_src[sp] * Hi[specConvIndx];
    }
    // Add the surface growth source to mom_src
    sd->surfaceGrowthMomSrc(k_sg, mom_fvPtr, mom_srcPtr);
    sd->oxidFragMomSrc(k_ox, k_o2, mom_fvPtr, mom_srcPtr);
    // Convert moment source terms back to CGS/SI units
    sd->convertFromMol(mom_srcPtr, momentsPtr);
    if (conserveMass) {
      // Difference between mass lost from fluid and mass gained to soot
      Real del_rho_dot = rho_src + mom_src[1] * sc.SootDensity;
      // Add that mass to H2
      soot_state(i, j, k, absorbIndxP) -= del_rho_dot * sc.mass_src_conv;
      rho_src -= del_rho_dot;
      eng_src -= Hi[absorbIndxN] * del_rho_dot;
    }
    // Add density source term
    soot_state(i, j, k, rhoIndx) += rho_src * sc.mass_src_conv;
    soot_state(i, j, k, engIndx) -= eng_src * sc.eng_src_conv;
    // Add moment source terms
    for (int mom = 0; mom < NUM_SOOT_MOMENTS + 1; ++mom) {
      const int peleIndx = sootIndx + mom;
      soot_state(i, j, k, peleIndx) += mom_src[mom];
    }
  });
}
