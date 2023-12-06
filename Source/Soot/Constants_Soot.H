#ifndef CONSTANTS_SOOT_H
#define CONSTANTS_SOOT_H

#include "PhysicsConstants.H"

// Note: These will change if a different soot surface chemistry model is
// implemented Define the number of species involved in soot model
#define NUM_SOOT_GS 8
#define NUM_SOOT_REACT 7

struct SootComps
{
  int qRhoIndx = -1;
  int qTempIndx = -1;
  int qSpecIndx = -1;
  int qSootIndx = -1;
  int rhoIndx = -1;
  int engIndx = -1;
  int specIndx = -1;
  int sootIndx = -1;
  bool checkIndices() const
  {
    bool check = true;
    if (
      qRhoIndx < 0 || qTempIndx < 0 || qSpecIndx < 0 || qSootIndx < 0 ||
      rhoIndx < 0 || engIndx < 0 || specIndx < 0 || sootIndx < 0) {
      check = false;
    }
    return check;
  }
};

// Local gas species index, not the same as what comes from PelePhysics
// Used as input to m_sootData->getPPIndx()
enum SootGasSpecIndx {
  indxH2 = 0,
  indxH,
  indxOH,
  indxH2O,
  indxCO,
  indxC2H2,
  indxO2,
  indxPAH, // PAH species
  numGasSpecs
};

enum SootIndx {
  indxSootS = 0, // Soot-*
  indxSootH,     // Soot-H
  numSootSpecs
};

struct SootConst
{
#ifdef PELELM_USE_SOOT
  amrex::Real len_conv = 100.;  // Turn m to cm
  amrex::Real rho_conv = 0.001; // Turn kg/m^3 to g/cm^3
  amrex::Real mu_conv = 10.;
  amrex::Real mass_src_conv = 1000.;
  amrex::Real eng_src_conv = 0.1;
#else
  amrex::Real len_conv = 1.;
  amrex::Real rho_conv = 1.;
  amrex::Real mu_conv = 1.;
  amrex::Real mass_src_conv = 1.;
  amrex::Real eng_src_conv = 1.;
#endif
  /// Soot density (g/cm^3)
  amrex::Real SootDensity = 1.8;
  /// Soot molar mass (g/mol)
  amrex::Real SootMolarMass = 12.011150; // Carbon
  /// Surface density of hydrogenate sites (1/cm^2)
  amrex::Real SootChi = 1.7E15;
  // Small weight used for initialization and clipping
  amrex::Real smallWeight = 1.E-26;
  /// Soot fractal dimension
  amrex::Real SootDf = 1.8;
  /// Coefficients for fit to small surface area change
  amrex::Real SootFitC = 2. / 3.;
  amrex::Real SootFitE = -0.2043;
  /// Free molecular size exponents
  amrex::Real SootAv = 1. - (2. / SootDf);
  amrex::Real SootAs = 3. / SootDf - 1.;
  // Increase in number of carbon atoms for surface reactions
  amrex::Real dVol = 2.;
#if NUM_SOOT_MOMENTS == 3
  /// Volume moment orders
  amrex::GpuArray<amrex::Real, 3> MomOrderV = {{0., 1., 0.}};
  /// Surface areas moment orders
  amrex::GpuArray<amrex::Real, 3> MomOrderS = {{0., 0., 1.}};
#elif NUM_SOOT_MOMENTS == 6
  amrex::GpuArray<amrex::Real, 6> MomOrderV = {{0., 1., 0., 2., 1., 0.}};
  amrex::GpuArray<amrex::Real, 6> MomOrderS = {{0., 0., 1., 0., 1., 2.}};
#endif
  // Volume of smallest soot particles, units L^3
  amrex::Real V0 =
    SootMolarMass / (pele::physics::Constants::Avna * SootDensity);
  // Surface area of smallest soot particles, units L^2
  amrex::Real S0 = std::pow(36. * M_PI, 1. / 3.) * std::pow(V0, 2. / 3.);
  // Surface density (mol of C)
  amrex::Real SootDensityC = SootChi * S0;
  // Pi*R/(2*A*rho_soot)
  amrex::Real colFact = M_PI * pele::physics::Constants::RU /
                        (2. * pele::physics::Constants::Avna * SootDensity);
  // (M_soot/(A*rho_soot))^2/3
  amrex::Real colFact23 = std::pow(V0, 2. / 3.);
  // (M_soot/(A*rho_soot))^1/6, units sqrt(L)
  amrex::Real colFact16 = std::pow(V0, 1. / 6.);
  // (6/Pi)^2/3
  amrex::Real colFactPi23 = std::pow(6. / M_PI, 2. / 3.);
};

#endif
