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

namespace SootConst
{
  // Local gas species index, not the same as what comes from PelePhysics
  // Use gasIndxtoPeleP to convert between
  enum GasSpecIndx {
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
  /// Soot density (g/cm^3)
  constexpr Real SootDensity = 1.8;
  /// Soot molar mass (g/mol)
  constexpr Real SootMolarMass = 12.011150; // Carbon
  /// Soot fractal dimension
  constexpr Real SootDf = 1.8;
  /// Surface density of hydrogenate sites (1/cm^2)
  constexpr Real SootChi = 1.7E15;
  /// Coefficients for fit to small surface area change
  constexpr Real SootFitC = 2./3.;
  constexpr Real SootFitE = -0.2043;
  /// Free molecular size exponents
  constexpr Real SootAv = 1. - (2./SootDf);
  constexpr Real SootAs = 3./SootDf - 1.;
  // Increase in number of carbon atoms for surface reactions
  constexpr Real dVol = 2.;
  // Small weight used for initialization and clipping
  constexpr Real smallWeight = 1.E-26;

  /// Avogadro's number in 1/mol
  // TODO: Should be pulling this value from PeleC/constants
  constexpr Real avogadros = 6.02214076E23;

  /// Volume moment orders
  extern AMREX_GPU_DEVICE_MANAGED Real MomOrderV[NUM_SOOT_MOMENTS];
  /// Surface areas moment orders
  extern AMREX_GPU_DEVICE_MANAGED Real MomOrderS[NUM_SOOT_MOMENTS];

  extern AMREX_GPU_DEVICE_MANAGED Real V0; // Volume of smallest soot particles
  extern AMREX_GPU_DEVICE_MANAGED Real S0; // Surface area of smallest soot particles
  extern AMREX_GPU_DEVICE_MANAGED Real SootDensityC; // Surface density (mol of C)
  extern AMREX_GPU_DEVICE_MANAGED Real colFact; // Pi*R/(2*A*rho_soot)
  extern AMREX_GPU_DEVICE_MANAGED Real colFact23; // (M_soot/(A*rho_soot))^2/3
  extern AMREX_GPU_DEVICE_MANAGED Real colFact16; // (M_soot/(A*rho_soot))^1/6
  extern AMREX_GPU_DEVICE_MANAGED Real colFactPi23; // (6/Pi)^2/3
  extern AMREX_GPU_DEVICE_MANAGED Real dimerVol;
  extern AMREX_GPU_DEVICE_MANAGED Real nuclVol;
  extern AMREX_GPU_DEVICE_MANAGED Real nuclSurf;
  extern AMREX_GPU_DEVICE_MANAGED Real unitConv[NUM_SOOT_MOMENTS+1];
  extern AMREX_GPU_DEVICE_MANAGED int refIndx[NUM_SOOT_GS];

  // Convert moments from CGS to mol of C
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void convertCGStoMol (Real moments[])
  {
    for (int i = 0; i != NUM_SOOT_VARS; ++i)
      moments[i] /= unitConv[i];
  }

  // Convert moments from mol of C to CGS
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void convertMoltoCGS (Real mom_src[])
  {
    for (int i = 0; i != NUM_SOOT_VARS; ++i)
      mom_src[i] *= unitConv[i];
  }

  // Convert moments from mol of C to CGS
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void convertMoltoCGS (Real mom_src[],
                        Real moments[])
  {
    for (int i = 0; i != NUM_SOOT_VARS; ++i) {
      mom_src[i] *= unitConv[i];
      moments[i] *= unitConv[i];
    }
  }

  // Initialize moment values to very small
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void initialSmallMomVals(Real moments[])
  {
    for (int i = 0; i != NUM_SOOT_MOMENTS; ++i) {
      moments[i] = smallWeight*
        std::pow(nuclVol, 3.*(MomOrderV[i] + 2./3.*MomOrderS[i]));
    }
    moments[NUM_SOOT_MOMENTS] = smallWeight;
    convertMoltoCGS(moments);
  }

};

#endif
