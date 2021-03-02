#ifndef _CONSTANTS_SOOT_H_
#define _CONSTANTS_SOOT_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

// TODO: This should be specified in makefile when
// different soot formation reactions are implemented
// Define the number of species involved in soot reactions
#define NUM_SOOT_GS 8
#define NUM_SOOT_REACT 7

using namespace amrex;

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
  bool checkIndices()
  {
    if (
      qRhoIndx < 0 || qTempIndx < 0 || qSpecIndx < 0 || qSootIndx < 0 ||
      rhoIndx < 0 || engIndx < 0 || specIndx < 0 || specIndx < 0 ||
      sootIndx < 0) {
      return false;
    }
    return true;
  }
};

// Local gas species index, not the same as what comes from PelePhysics
// Use gasIndxtoPeleP to convert between
enum SootGasSpecIndx {
  // Hardcoded indices for species in m_specSootRefs
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
#ifdef SOOT_PELE_LM
  Real len_conv = 100.; // Turn m to cm
  Real rho_conv = 0.001; // Turn kg/m^3 to g/cm^3
  Real mu_conv = 10.;
  Real mass_src_conv = 1000.;
  Real eng_src_conv = 1.E-1;
#else
  Real len_conv = 1.;
  Real rho_conv = 1.;
  Real mu_conv = 1.;
  Real mass_src_conv = 1.;
  Real eng_src_conv = 1.;
#endif
  /// Soot density (g/cm^3)
  Real SootDensity = 1.8;
  /// Soot molar mass (g/mol)
  Real SootMolarMass = 12.011150; // Carbon
  /// Surface density of hydrogenate sites (1/cm^2)
  Real SootChi = 1.7E15;
  // Small weight used for initialization and clipping
  Real smallWeight = 1.E-26;
  /// Soot fractal dimension
  Real SootDf = 1.8;
  /// Coefficients for fit to small surface area change
  Real SootFitC = 2. / 3.;
  Real SootFitE = -0.2043;
  /// Free molecular size exponents
  Real SootAv = 1. - (2. / SootDf);
  Real SootAs = 3. / SootDf - 1.;
  // Increase in number of carbon atoms for surface reactions
  Real dVol = 2.;
  /// Avogadro's number in 1/mol
  // TODO: Should be pulling this value from PeleMP/constants
  Real avogadros = 6.02214076E23;
#if NUM_SOOT_MOMENTS == 3
  /// Volume moment orders
  GpuArray<Real, 3> MomOrderV = {{0., 1., 0.}};
  /// Surface areas moment orders
  GpuArray<Real, 3> MomOrderS = {{0., 0., 1.}};
#elif NUM_SOOT_MOMENTS == 6
  GpuArray<Real, 6> MomOrderV = {{0., 1., 0., 2., 1., 0.}};
  GpuArray<Real, 6> MomOrderS = {{0., 0., 1., 0., 1., 2.}};
#endif
  // Volume of smallest soot particles, units L^3
  Real V0 = SootMolarMass / (avogadros * SootDensity);
  // Surface area of smallest soot particles, units L^2
  Real S0 = std::pow(36. * M_PI, 1. / 3.) * std::pow(V0, 2. / 3.);
  // Surface density (mol of C)
  Real SootDensityC = SootChi * S0;
  // Pi*R/(2*A*rho_soot)
  Real colFact = M_PI * EOS::RU / (2. * avogadros * SootDensity);
  // (M_soot/(A*rho_soot))^2/3
  Real colFact23 = std::pow(V0, 2. / 3.);
  // (M_soot/(A*rho_soot))^1/6, units sqrt(L)
  Real colFact16 = std::pow(V0, 1. / 6.);
  // (6/Pi)^2/3
  Real colFactPi23 = std::pow(6. / M_PI, 2. / 3.);
};

#endif
