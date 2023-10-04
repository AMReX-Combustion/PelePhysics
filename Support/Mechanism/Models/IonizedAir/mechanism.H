#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  O
1  N
2  E
*/

// Species
#define O2_ID 0
#define N2_ID 1
#define E_ID 2
#define On_ID 3
#define O2n_ID 4

#define NUM_ELEMENTS 3
#define NUM_SPECIES 5
#define NUM_IONS 3
#define NUM_REACTIONS 0

#define NUM_FIT 4

//  ALWAYS on CPU stuff -- can have different def depending on if we are CPU or
//  GPU based. Defined in mechanism.cpp
void atomicWeight(amrex::Real* awt);
//  MISC
void CKAWT(amrex::Real* awt);
void CKNCF(int* ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
void GET_RMAP(int* _rmap);
void CKINU(const int i, int& nspec, int* ki, int* nu);
void CKKFKR(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real* x,
  amrex::Real* q_f,
  amrex::Real* q_r);
void progressRateFR(
  amrex::Real* q_f, amrex::Real* q_r, amrex::Real* sc, amrex::Real T);
//  SPARSE INFORMATION
void SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS);
void SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP);
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS);
void SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base);

// A few mechanism parameters
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKINDX(int& mm, int& kk, int& ii, int& nfit)
{
  mm = 3;
  kk = 5;
  ii = 0;
  nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[5] = {
  0.0312519532470779,    // O2
  0.0356964374955379,    // N2
  1822.8884868472639482, // E
  0.0625017634094834,    // On
  0.0312514174667246,    // O2n
};
#endif
const amrex::Real h_global_imw[5] = {
  0.0312519532470779,    // O2
  0.0356964374955379,    // N2
  1822.8884868472639482, // E
  0.0625017634094834,    // On
  0.0312514174667246,    // O2n
};

//  molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[5] = {
  31.998000, // O2
  28.014000, // N2
  0.000549,  // E
  15.999549, // On
  31.998549, // O2n
};
#endif
const amrex::Real h_global_mw[5] = {
  31.998000, // O2
  28.014000, // N2
  0.000549,  // E
  15.999549, // On
  31.998549, // O2n
};

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_imw(amrex::Real* imw_new)
{
  imw_new[0] = 0.0312519532470779;    // O2
  imw_new[1] = 0.0356964374955379;    // N2
  imw_new[2] = 1822.8884868472639482; // E
  imw_new[3] = 0.0625017634094834;    // On
  imw_new[4] = 0.0312514174667246;    // O2n
}

//  inverse molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
imw(const int n)
{
#if AMREX_DEVICE_COMPILE
  return global_imw[n];
#else
  return h_global_imw[n];
#endif
}
//  molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_mw(amrex::Real* mw_new)
{
  mw_new[0] = 31.998000; // O2
  mw_new[1] = 28.014000; // N2
  mw_new[2] = 0.000549;  // E
  mw_new[3] = 15.999549; // On
  mw_new[4] = 31.998549; // O2n
}

//  molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
mw(const int n)
{
#if AMREX_DEVICE_COMPILE
  return global_mw[n];
#else
  return h_global_mw[n];
#endif
}

//  Returns R, Rc, Patm
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRP(amrex::Real& ru, amrex::Real& ruc, amrex::Real& pa)
{
  ru = 8.31446261815324e+07;
  ruc = 1.98721558317399615845;
  pa = 1.01325e+06;
}

// compute Cv/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cv_R(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: O2
    species[0] = +2.21293640e+00 + 1.12748635e-03 * tc[1] -
                 5.75615047e-07 * tc[2] + 1.31387723e-09 * tc[3] -
                 8.76855392e-13 * tc[4];
    // species 1: N2
    species[1] = +2.29867700e+00 + 1.40824000e-03 * tc[1] -
                 3.96322200e-06 * tc[2] + 5.64151500e-09 * tc[3] -
                 2.44485500e-12 * tc[4];
    // species 3: On
    species[3] = +1.90805921e+00 - 1.69804907e-03 * tc[1] +
                 2.98069955e-06 * tc[2] - 2.43835127e-09 * tc[3] +
                 7.61229311e-13 * tc[4];
    // species 4: O2n
    species[4] = +2.66442522e+00 - 9.28741138e-04 * tc[1] +
                 6.45477082e-06 * tc[2] - 7.74703380e-09 * tc[3] +
                 2.93332662e-12 * tc[4];
  } else {
    // species 0: O2
    species[0] = +2.69757819e+00 + 6.13519689e-04 * tc[1] -
                 1.25884199e-07 * tc[2] + 1.77528148e-11 * tc[3] -
                 1.13643531e-15 * tc[4];
    // species 1: N2
    species[1] = +1.92664000e+00 + 1.48797700e-03 * tc[1] -
                 5.68476100e-07 * tc[2] + 1.00970400e-10 * tc[3] -
                 6.75335100e-15 * tc[4];
    // species 3: On
    species[3] = +1.54474869e+00 - 4.66695513e-05 * tc[1] +
                 1.84912357e-08 * tc[2] - 3.18159223e-12 * tc[3] +
                 1.98962956e-16 * tc[4];
    // species 4: O2n
    species[4] = +2.95666294e+00 + 5.98141823e-04 * tc[1] -
                 2.12133905e-07 * tc[2] + 3.63267581e-11 * tc[3] -
                 2.24989228e-15 * tc[4];
  }

  // species with no change at a midpoint T
  // species 2: E
  species[2] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4];
}

// compute Cp/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cp_R(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: O2
    species[0] = +3.21293640e+00 + 1.12748635e-03 * tc[1] -
                 5.75615047e-07 * tc[2] + 1.31387723e-09 * tc[3] -
                 8.76855392e-13 * tc[4];
    // species 1: N2
    species[1] = +3.29867700e+00 + 1.40824000e-03 * tc[1] -
                 3.96322200e-06 * tc[2] + 5.64151500e-09 * tc[3] -
                 2.44485500e-12 * tc[4];
    // species 3: On
    species[3] = +2.90805921e+00 - 1.69804907e-03 * tc[1] +
                 2.98069955e-06 * tc[2] - 2.43835127e-09 * tc[3] +
                 7.61229311e-13 * tc[4];
    // species 4: O2n
    species[4] = +3.66442522e+00 - 9.28741138e-04 * tc[1] +
                 6.45477082e-06 * tc[2] - 7.74703380e-09 * tc[3] +
                 2.93332662e-12 * tc[4];
  } else {
    // species 0: O2
    species[0] = +3.69757819e+00 + 6.13519689e-04 * tc[1] -
                 1.25884199e-07 * tc[2] + 1.77528148e-11 * tc[3] -
                 1.13643531e-15 * tc[4];
    // species 1: N2
    species[1] = +2.92664000e+00 + 1.48797700e-03 * tc[1] -
                 5.68476100e-07 * tc[2] + 1.00970400e-10 * tc[3] -
                 6.75335100e-15 * tc[4];
    // species 3: On
    species[3] = +2.54474869e+00 - 4.66695513e-05 * tc[1] +
                 1.84912357e-08 * tc[2] - 3.18159223e-12 * tc[3] +
                 1.98962956e-16 * tc[4];
    // species 4: O2n
    species[4] = +3.95666294e+00 + 5.98141823e-04 * tc[1] -
                 2.12133905e-07 * tc[2] + 3.63267581e-11 * tc[3] -
                 2.24989228e-15 * tc[4];
  }

  // species with no change at a midpoint T
  // species 2: E
  species[2] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4];
}

// compute the g/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
gibbs(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: O2
    species[0] = -1.005249020000000e+03 * invT - 2.821801190000000e+00 -
                 3.212936400000000e+00 * tc[0] - 5.637431750000000e-04 * tc[1] +
                 9.593584116666666e-08 * tc[2] - 1.094897691666667e-10 * tc[3] +
                 4.384276960000000e-14 * tc[4];
    // species 1: N2
    species[1] = -1.020900000000000e+03 * invT - 6.516950000000001e-01 -
                 3.298677000000000e+00 * tc[0] - 7.041200000000000e-04 * tc[1] +
                 6.605369999999999e-07 * tc[2] - 4.701262500000001e-10 * tc[3] +
                 1.222427500000000e-13 * tc[4];
    // species 3: On
    species[3] = +1.143577170000000e+04 * invT + 1.046682399999996e-01 -
                 2.908059210000000e+00 * tc[0] + 8.490245350000000e-04 * tc[1] -
                 4.967832583333333e-07 * tc[2] + 2.031959391666667e-10 * tc[3] -
                 3.806146555000000e-14 * tc[4];
    // species 4: O2n
    species[4] = -6.870769830000000e+03 * invT - 6.869815900000003e-01 -
                 3.664425220000000e+00 * tc[0] + 4.643705690000000e-04 * tc[1] -
                 1.075795136666667e-06 * tc[2] + 6.455861500000001e-10 * tc[3] -
                 1.466663310000000e-13 * tc[4];
  } else {
    // species 0: O2
    species[0] = -1.233930180000000e+03 * invT + 5.084126000000002e-01 -
                 3.697578190000000e+00 * tc[0] - 3.067598445000000e-04 * tc[1] +
                 2.098069983333333e-08 * tc[2] - 1.479401233333333e-12 * tc[3] +
                 5.682176550000000e-17 * tc[4];
    // species 1: N2
    species[1] = -9.227977000000000e+02 * invT - 3.053888000000000e+00 -
                 2.926640000000000e+00 * tc[0] - 7.439885000000000e-04 * tc[1] +
                 9.474601666666666e-08 * tc[2] - 8.414199999999999e-12 * tc[3] +
                 3.376675500000000e-16 * tc[4];
    // species 3: On
    species[3] = +1.150420890000000e+04 * invT - 1.976561460000000e+00 -
                 2.544748690000000e+00 * tc[0] + 2.333477565000000e-05 * tc[1] -
                 3.081872616666667e-09 * tc[2] + 2.651326858333333e-13 * tc[3] -
                 9.948147799999999e-18 * tc[4];
    // species 4: O2n
    species[4] = -7.062872290000000e+03 * invT + 1.677952770000000e+00 -
                 3.956662940000000e+00 * tc[0] - 2.990709115000000e-04 * tc[1] +
                 3.535565083333334e-08 * tc[2] - 3.027229841666666e-12 * tc[3] +
                 1.124946140000000e-16 * tc[4];
  }

  // species with no change at a midpoint T
  // species 2: E
  species[2] = -7.453750000000000e+02 * invT + 1.422081220000000e+01 -
               2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
               0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
               0.000000000000000e+00 * tc[4];
}

// compute the a/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
helmholtz(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: O2
    species[0] = -1.00524902e+03 * invT - 3.82180119e+00 -
                 3.21293640e+00 * tc[0] - 5.63743175e-04 * tc[1] +
                 9.59358412e-08 * tc[2] - 1.09489769e-10 * tc[3] +
                 4.38427696e-14 * tc[4];
    // species 1: N2
    species[1] = -1.02090000e+03 * invT - 1.65169500e+00 -
                 3.29867700e+00 * tc[0] - 7.04120000e-04 * tc[1] +
                 6.60537000e-07 * tc[2] - 4.70126250e-10 * tc[3] +
                 1.22242750e-13 * tc[4];
    // species 3: On
    species[3] = +1.14357717e+04 * invT - 8.95331760e-01 -
                 2.90805921e+00 * tc[0] + 8.49024535e-04 * tc[1] -
                 4.96783258e-07 * tc[2] + 2.03195939e-10 * tc[3] -
                 3.80614656e-14 * tc[4];
    // species 4: O2n
    species[4] = -6.87076983e+03 * invT - 1.68698159e+00 -
                 3.66442522e+00 * tc[0] + 4.64370569e-04 * tc[1] -
                 1.07579514e-06 * tc[2] + 6.45586150e-10 * tc[3] -
                 1.46666331e-13 * tc[4];
  } else {
    // species 0: O2
    species[0] = -1.23393018e+03 * invT - 4.91587400e-01 -
                 3.69757819e+00 * tc[0] - 3.06759845e-04 * tc[1] +
                 2.09806998e-08 * tc[2] - 1.47940123e-12 * tc[3] +
                 5.68217655e-17 * tc[4];
    // species 1: N2
    species[1] = -9.22797700e+02 * invT - 4.05388800e+00 -
                 2.92664000e+00 * tc[0] - 7.43988500e-04 * tc[1] +
                 9.47460167e-08 * tc[2] - 8.41420000e-12 * tc[3] +
                 3.37667550e-16 * tc[4];
    // species 3: On
    species[3] = +1.15042089e+04 * invT - 2.97656146e+00 -
                 2.54474869e+00 * tc[0] + 2.33347756e-05 * tc[1] -
                 3.08187262e-09 * tc[2] + 2.65132686e-13 * tc[3] -
                 9.94814780e-18 * tc[4];
    // species 4: O2n
    species[4] = -7.06287229e+03 * invT + 6.77952770e-01 -
                 3.95666294e+00 * tc[0] - 2.99070912e-04 * tc[1] +
                 3.53556508e-08 * tc[2] - 3.02722984e-12 * tc[3] +
                 1.12494614e-16 * tc[4];
  }

  // species with no change at a midpoint T
  // species 2: E
  species[2] = -7.45375000e+02 * invT + 1.32208122e+01 -
               2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
               0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
               0.00000000e+00 * tc[4];
}

// compute the e/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesInternalEnergy(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: O2
    species[0] = +2.21293640e+00 + 5.63743175e-04 * tc[1] -
                 1.91871682e-07 * tc[2] + 3.28469308e-10 * tc[3] -
                 1.75371078e-13 * tc[4] - 1.00524902e+03 * invT;
    // species 1: N2
    species[1] = +2.29867700e+00 + 7.04120000e-04 * tc[1] -
                 1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
                 4.88971000e-13 * tc[4] - 1.02090000e+03 * invT;
    // species 3: On
    species[3] = +1.90805921e+00 - 8.49024535e-04 * tc[1] +
                 9.93566517e-07 * tc[2] - 6.09587818e-10 * tc[3] +
                 1.52245862e-13 * tc[4] + 1.14357717e+04 * invT;
    // species 4: O2n
    species[4] = +2.66442522e+00 - 4.64370569e-04 * tc[1] +
                 2.15159027e-06 * tc[2] - 1.93675845e-09 * tc[3] +
                 5.86665324e-13 * tc[4] - 6.87076983e+03 * invT;
  } else {
    // species 0: O2
    species[0] = +2.69757819e+00 + 3.06759845e-04 * tc[1] -
                 4.19613997e-08 * tc[2] + 4.43820370e-12 * tc[3] -
                 2.27287062e-16 * tc[4] - 1.23393018e+03 * invT;
    // species 1: N2
    species[1] = +1.92664000e+00 + 7.43988500e-04 * tc[1] -
                 1.89492033e-07 * tc[2] + 2.52426000e-11 * tc[3] -
                 1.35067020e-15 * tc[4] - 9.22797700e+02 * invT;
    // species 3: On
    species[3] = +1.54474869e+00 - 2.33347756e-05 * tc[1] +
                 6.16374523e-09 * tc[2] - 7.95398057e-13 * tc[3] +
                 3.97925912e-17 * tc[4] + 1.15042089e+04 * invT;
    // species 4: O2n
    species[4] = +2.95666294e+00 + 2.99070912e-04 * tc[1] -
                 7.07113017e-08 * tc[2] + 9.08168952e-12 * tc[3] -
                 4.49978456e-16 * tc[4] - 7.06287229e+03 * invT;
  }

  // species with no change at a midpoint T
  // species 2: E
  species[2] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
}

// compute the h/(RT) at the given temperature (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEnthalpy(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: O2
    species[0] = +3.21293640e+00 + 5.63743175e-04 * tc[1] -
                 1.91871682e-07 * tc[2] + 3.28469308e-10 * tc[3] -
                 1.75371078e-13 * tc[4] - 1.00524902e+03 * invT;
    // species 1: N2
    species[1] = +3.29867700e+00 + 7.04120000e-04 * tc[1] -
                 1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
                 4.88971000e-13 * tc[4] - 1.02090000e+03 * invT;
    // species 3: On
    species[3] = +2.90805921e+00 - 8.49024535e-04 * tc[1] +
                 9.93566517e-07 * tc[2] - 6.09587818e-10 * tc[3] +
                 1.52245862e-13 * tc[4] + 1.14357717e+04 * invT;
    // species 4: O2n
    species[4] = +3.66442522e+00 - 4.64370569e-04 * tc[1] +
                 2.15159027e-06 * tc[2] - 1.93675845e-09 * tc[3] +
                 5.86665324e-13 * tc[4] - 6.87076983e+03 * invT;
  } else {
    // species 0: O2
    species[0] = +3.69757819e+00 + 3.06759845e-04 * tc[1] -
                 4.19613997e-08 * tc[2] + 4.43820370e-12 * tc[3] -
                 2.27287062e-16 * tc[4] - 1.23393018e+03 * invT;
    // species 1: N2
    species[1] = +2.92664000e+00 + 7.43988500e-04 * tc[1] -
                 1.89492033e-07 * tc[2] + 2.52426000e-11 * tc[3] -
                 1.35067020e-15 * tc[4] - 9.22797700e+02 * invT;
    // species 3: On
    species[3] = +2.54474869e+00 - 2.33347756e-05 * tc[1] +
                 6.16374523e-09 * tc[2] - 7.95398057e-13 * tc[3] +
                 3.97925912e-17 * tc[4] + 1.15042089e+04 * invT;
    // species 4: O2n
    species[4] = +3.95666294e+00 + 2.99070912e-04 * tc[1] -
                 7.07113017e-08 * tc[2] + 9.08168952e-12 * tc[3] -
                 4.49978456e-16 * tc[4] - 7.06287229e+03 * invT;
  }

  // species with no change at a midpoint T
  // species 2: E
  species[2] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
}

// compute the S/R at the given temperature (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEntropy(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: O2
    species[0] = +3.21293640e+00 * tc[0] + 1.12748635e-03 * tc[1] -
                 2.87807523e-07 * tc[2] + 4.37959077e-10 * tc[3] -
                 2.19213848e-13 * tc[4] + 6.03473759e+00;
    // species 1: N2
    species[1] = +3.29867700e+00 * tc[0] + 1.40824000e-03 * tc[1] -
                 1.98161100e-06 * tc[2] + 1.88050500e-09 * tc[3] -
                 6.11213750e-13 * tc[4] + 3.95037200e+00;
    // species 3: On
    species[3] = +2.90805921e+00 * tc[0] - 1.69804907e-03 * tc[1] +
                 1.49034977e-06 * tc[2] - 8.12783757e-10 * tc[3] +
                 1.90307328e-13 * tc[4] + 2.80339097e+00;
    // species 4: O2n
    species[4] = +3.66442522e+00 * tc[0] - 9.28741138e-04 * tc[1] +
                 3.22738541e-06 * tc[2] - 2.58234460e-09 * tc[3] +
                 7.33331655e-13 * tc[4] + 4.35140681e+00;
  } else {
    // species 0: O2
    species[0] = +3.69757819e+00 * tc[0] + 6.13519689e-04 * tc[1] -
                 6.29420995e-08 * tc[2] + 5.91760493e-12 * tc[3] -
                 2.84108828e-16 * tc[4] + 3.18916559e+00;
    // species 1: N2
    species[1] = +2.92664000e+00 * tc[0] + 1.48797700e-03 * tc[1] -
                 2.84238050e-07 * tc[2] + 3.36568000e-11 * tc[3] -
                 1.68833775e-15 * tc[4] + 5.98052800e+00;
    // species 3: On
    species[3] = +2.54474869e+00 * tc[0] - 4.66695513e-05 * tc[1] +
                 9.24561785e-09 * tc[2] - 1.06053074e-12 * tc[3] +
                 4.97407390e-17 * tc[4] + 4.52131015e+00;
    // species 4: O2n
    species[4] = +3.95666294e+00 * tc[0] + 5.98141823e-04 * tc[1] -
                 1.06066953e-07 * tc[2] + 1.21089194e-11 * tc[3] -
                 5.62473070e-16 * tc[4] + 2.27871017e+00;
  }

  // species with no change at a midpoint T
  // species 2: E
  species[2] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 1.17208122e+01;
}

// Returns the mean specific heat at CP (Eq. 33)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cpbl)
{
  amrex::Real result = 0;
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real cpor[5];                                //  temporary storage
  cp_R(cpor, tc);

  // perform dot product
  for (int id = 0; id < 5; ++id) {
    result += x[id] * cpor[id];
  }

  cpbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CP (Eq. 34)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBS(const amrex::Real T, const amrex::Real y[], amrex::Real& cpbs)
{
  amrex::Real result = 0.0;
  const amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache

  // compute Cp/R at the given temperature

  // species with midpoint at T=1000 kelvin
  if (tT < 1000) {
    // species 0: O2
    result +=
      y[0] *
      (+3.21293640e+00 + 1.12748635e-03 * tc[1] - 5.75615047e-07 * tc[2] +
       1.31387723e-09 * tc[3] - 8.76855392e-13 * tc[4]) *
      0.0312519532470779;
    // species 1: N2
    result +=
      y[1] *
      (+3.29867700e+00 + 1.40824000e-03 * tc[1] - 3.96322200e-06 * tc[2] +
       5.64151500e-09 * tc[3] - 2.44485500e-12 * tc[4]) *
      0.0356964374955379;
    // species 3: On
    result +=
      y[3] *
      (+2.90805921e+00 - 1.69804907e-03 * tc[1] + 2.98069955e-06 * tc[2] -
       2.43835127e-09 * tc[3] + 7.61229311e-13 * tc[4]) *
      0.0625017634094834;
    // species 4: O2n
    result +=
      y[4] *
      (+3.66442522e+00 - 9.28741138e-04 * tc[1] + 6.45477082e-06 * tc[2] -
       7.74703380e-09 * tc[3] + 2.93332662e-12 * tc[4]) *
      0.0312514174667246;
  } else {
    // species 0: O2
    result +=
      y[0] *
      (+3.69757819e+00 + 6.13519689e-04 * tc[1] - 1.25884199e-07 * tc[2] +
       1.77528148e-11 * tc[3] - 1.13643531e-15 * tc[4]) *
      0.0312519532470779;
    // species 1: N2
    result +=
      y[1] *
      (+2.92664000e+00 + 1.48797700e-03 * tc[1] - 5.68476100e-07 * tc[2] +
       1.00970400e-10 * tc[3] - 6.75335100e-15 * tc[4]) *
      0.0356964374955379;
    // species 3: On
    result +=
      y[3] *
      (+2.54474869e+00 - 4.66695513e-05 * tc[1] + 1.84912357e-08 * tc[2] -
       3.18159223e-12 * tc[3] + 1.98962956e-16 * tc[4]) *
      0.0625017634094834;
    // species 4: O2n
    result +=
      y[4] *
      (+3.95666294e+00 + 5.98141823e-04 * tc[1] - 2.12133905e-07 * tc[2] +
       3.63267581e-11 * tc[3] - 2.24989228e-15 * tc[4]) *
      0.0312514174667246;
  }

  // species with no change at a midpoint T
  // species 2: E
  result += y[2] *
            (+2.50000000e+00 + 0.00000000e+00 * tc[1] + 0.00000000e+00 * tc[2] +
             0.00000000e+00 * tc[3] + 0.00000000e+00 * tc[4]) *
            1822.8884868472639482;

  cpbs = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 35)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cvbl)
{
  amrex::Real result = 0;
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real cvor[5];                                //  temporary storage
  cv_R(cvor, tc);

  // perform dot product
  for (int id = 0; id < 5; ++id) {
    result += x[id] * cvor[id];
  }

  cvbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 36)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBS(const amrex::Real T, const amrex::Real y[], amrex::Real& cvbs)
{
  amrex::Real result = 0.0;
  const amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  // compute Cv/R at the given temperature

  // species with midpoint at T=1000 kelvin
  if (tT < 1000) {
    // species 0: O2
    result +=
      y[0] *
      (+2.21293640e+00 + 1.12748635e-03 * tc[1] - 5.75615047e-07 * tc[2] +
       1.31387723e-09 * tc[3] - 8.76855392e-13 * tc[4]) *
      0.0312519532470779;
    // species 1: N2
    result +=
      y[1] *
      (+2.29867700e+00 + 1.40824000e-03 * tc[1] - 3.96322200e-06 * tc[2] +
       5.64151500e-09 * tc[3] - 2.44485500e-12 * tc[4]) *
      0.0356964374955379;
    // species 3: On
    result +=
      y[3] *
      (+1.90805921e+00 - 1.69804907e-03 * tc[1] + 2.98069955e-06 * tc[2] -
       2.43835127e-09 * tc[3] + 7.61229311e-13 * tc[4]) *
      0.0625017634094834;
    // species 4: O2n
    result +=
      y[4] *
      (+2.66442522e+00 - 9.28741138e-04 * tc[1] + 6.45477082e-06 * tc[2] -
       7.74703380e-09 * tc[3] + 2.93332662e-12 * tc[4]) *
      0.0312514174667246;
  } else {
    // species 0: O2
    result +=
      y[0] *
      (+2.69757819e+00 + 6.13519689e-04 * tc[1] - 1.25884199e-07 * tc[2] +
       1.77528148e-11 * tc[3] - 1.13643531e-15 * tc[4]) *
      0.0312519532470779;
    // species 1: N2
    result +=
      y[1] *
      (+1.92664000e+00 + 1.48797700e-03 * tc[1] - 5.68476100e-07 * tc[2] +
       1.00970400e-10 * tc[3] - 6.75335100e-15 * tc[4]) *
      0.0356964374955379;
    // species 3: On
    result +=
      y[3] *
      (+1.54474869e+00 - 4.66695513e-05 * tc[1] + 1.84912357e-08 * tc[2] -
       3.18159223e-12 * tc[3] + 1.98962956e-16 * tc[4]) *
      0.0625017634094834;
    // species 4: O2n
    result +=
      y[4] *
      (+2.95666294e+00 + 5.98141823e-04 * tc[1] - 2.12133905e-07 * tc[2] +
       3.63267581e-11 * tc[3] - 2.24989228e-15 * tc[4]) *
      0.0312514174667246;
  }

  // species with no change at a midpoint T
  // species 2: E
  result += y[2] *
            (+1.50000000e+00 + 0.00000000e+00 * tc[1] + 0.00000000e+00 * tc[2] +
             0.00000000e+00 * tc[3] + 0.00000000e+00 * tc[4]) *
            1822.8884868472639482;

  cvbs = result * 8.31446261815324e+07;
}

// Returns the mean enthalpy of the mixture in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBML(const amrex::Real T, const amrex::Real x[], amrex::Real& hbml)
{
  amrex::Real result = 0;
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real hml[5];                                 //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesEnthalpy(hml, tc);

  // perform dot product
  for (int id = 0; id < 5; ++id) {
    result += x[id] * hml[id];
  }

  hbml = result * RT;
}

// Returns mean enthalpy of mixture in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBMS(const amrex::Real T, const amrex::Real y[], amrex::Real& hbms)
{
  amrex::Real result = 0.0;
  const amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (tT < 1000) {
    // species 0: O2
    result += y[0] *
              (+3.21293640e+00 + 5.63743175e-04 * tc[1] -
               1.91871682e-07 * tc[2] + 3.28469308e-10 * tc[3] -
               1.75371078e-13 * tc[4] - 1.00524902e+03 * invT) *
              0.0312519532470779;
    // species 1: N2
    result += y[1] *
              (+3.29867700e+00 + 7.04120000e-04 * tc[1] -
               1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
               4.88971000e-13 * tc[4] - 1.02090000e+03 * invT) *
              0.0356964374955379;
    // species 3: On
    result += y[3] *
              (+2.90805921e+00 - 8.49024535e-04 * tc[1] +
               9.93566517e-07 * tc[2] - 6.09587818e-10 * tc[3] +
               1.52245862e-13 * tc[4] + 1.14357717e+04 * invT) *
              0.0625017634094834;
    // species 4: O2n
    result += y[4] *
              (+3.66442522e+00 - 4.64370569e-04 * tc[1] +
               2.15159027e-06 * tc[2] - 1.93675845e-09 * tc[3] +
               5.86665324e-13 * tc[4] - 6.87076983e+03 * invT) *
              0.0312514174667246;
  } else {
    // species 0: O2
    result += y[0] *
              (+3.69757819e+00 + 3.06759845e-04 * tc[1] -
               4.19613997e-08 * tc[2] + 4.43820370e-12 * tc[3] -
               2.27287062e-16 * tc[4] - 1.23393018e+03 * invT) *
              0.0312519532470779;
    // species 1: N2
    result += y[1] *
              (+2.92664000e+00 + 7.43988500e-04 * tc[1] -
               1.89492033e-07 * tc[2] + 2.52426000e-11 * tc[3] -
               1.35067020e-15 * tc[4] - 9.22797700e+02 * invT) *
              0.0356964374955379;
    // species 3: On
    result += y[3] *
              (+2.54474869e+00 - 2.33347756e-05 * tc[1] +
               6.16374523e-09 * tc[2] - 7.95398057e-13 * tc[3] +
               3.97925912e-17 * tc[4] + 1.15042089e+04 * invT) *
              0.0625017634094834;
    // species 4: O2n
    result += y[4] *
              (+3.95666294e+00 + 2.99070912e-04 * tc[1] -
               7.07113017e-08 * tc[2] + 9.08168952e-12 * tc[3] -
               4.49978456e-16 * tc[4] - 7.06287229e+03 * invT) *
              0.0312514174667246;
  }

  // species with no change at a midpoint T
  // species 2: E
  result +=
    y[2] *
    (+2.50000000e+00 + 0.00000000e+00 * tc[1] + 0.00000000e+00 * tc[2] +
     0.00000000e+00 * tc[3] + 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT) *
    1822.8884868472639482;

  const amrex::Real RT = 8.31446261815324e+07 * tT; // R*T

  hbms = result * RT;
}

// get mean internal energy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBML(const amrex::Real T, const amrex::Real x[], amrex::Real& ubml)
{
  amrex::Real result = 0;
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real uml[5];                                 //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesInternalEnergy(uml, tc);

  // perform dot product
  for (int id = 0; id < 5; ++id) {
    result += x[id] * uml[id];
  }

  ubml = result * RT;
}

// get mean internal energy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBMS(const amrex::Real T, const amrex::Real y[], amrex::Real& ubms)
{
  amrex::Real result = 0.0;
  const amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (tT < 1000) {
    // species 0: O2
    result += y[0] *
              (+2.21293640e+00 + 5.63743175e-04 * tc[1] -
               1.91871682e-07 * tc[2] + 3.28469308e-10 * tc[3] -
               1.75371078e-13 * tc[4] - 1.00524902e+03 * invT) *
              0.0312519532470779;
    // species 1: N2
    result += y[1] *
              (+2.29867700e+00 + 7.04120000e-04 * tc[1] -
               1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
               4.88971000e-13 * tc[4] - 1.02090000e+03 * invT) *
              0.0356964374955379;
    // species 3: On
    result += y[3] *
              (+1.90805921e+00 - 8.49024535e-04 * tc[1] +
               9.93566517e-07 * tc[2] - 6.09587818e-10 * tc[3] +
               1.52245862e-13 * tc[4] + 1.14357717e+04 * invT) *
              0.0625017634094834;
    // species 4: O2n
    result += y[4] *
              (+2.66442522e+00 - 4.64370569e-04 * tc[1] +
               2.15159027e-06 * tc[2] - 1.93675845e-09 * tc[3] +
               5.86665324e-13 * tc[4] - 6.87076983e+03 * invT) *
              0.0312514174667246;
  } else {
    // species 0: O2
    result += y[0] *
              (+2.69757819e+00 + 3.06759845e-04 * tc[1] -
               4.19613997e-08 * tc[2] + 4.43820370e-12 * tc[3] -
               2.27287062e-16 * tc[4] - 1.23393018e+03 * invT) *
              0.0312519532470779;
    // species 1: N2
    result += y[1] *
              (+1.92664000e+00 + 7.43988500e-04 * tc[1] -
               1.89492033e-07 * tc[2] + 2.52426000e-11 * tc[3] -
               1.35067020e-15 * tc[4] - 9.22797700e+02 * invT) *
              0.0356964374955379;
    // species 3: On
    result += y[3] *
              (+1.54474869e+00 - 2.33347756e-05 * tc[1] +
               6.16374523e-09 * tc[2] - 7.95398057e-13 * tc[3] +
               3.97925912e-17 * tc[4] + 1.15042089e+04 * invT) *
              0.0625017634094834;
    // species 4: O2n
    result += y[4] *
              (+2.95666294e+00 + 2.99070912e-04 * tc[1] -
               7.07113017e-08 * tc[2] + 9.08168952e-12 * tc[3] -
               4.49978456e-16 * tc[4] - 7.06287229e+03 * invT) *
              0.0312514174667246;
  }

  // species with no change at a midpoint T
  // species 2: E
  result +=
    y[2] *
    (+1.50000000e+00 + 0.00000000e+00 * tc[1] + 0.00000000e+00 * tc[2] +
     0.00000000e+00 * tc[3] + 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT) *
    1822.8884868472639482;

  const amrex::Real RT = 8.31446261815324e+07 * tT; // R*T

  ubms = result * RT;
}

// get mixture entropy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSBML(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real& sbml)
{
  amrex::Real result = 0;
  // Log of normalized pressure in cgs units dynes/cm^2 by Patm
  amrex::Real logPratio = log(P / 1013250.0);
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real sor[5]; //  temporary storage
  speciesEntropy(sor, tc);

  // Compute Eq 42
  for (int id = 0; id < 5; ++id) {
    result += x[id] * (sor[id] - log((x[id] + 1e-100)) - logPratio);
  }

  sbml = result * 8.31446261815324e+07;
}

// get mixture entropy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSBMS(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real& sbms)
{
  amrex::Real result = 0;
  // Log of normalized pressure in cgs units dynes/cm^2 by Patm
  amrex::Real logPratio = log(P / 1013250.0);
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real sor[5];  //  temporary storage
  amrex::Real x[5];    //  need a ytx conversion
  amrex::Real YOW = 0; // See Eq 4, 6 in CK Manual

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 5; i++) {
    YOW += y[i] * imw(i);
  }
  // Now compute y to x conversion
  x[0] = y[0] / (31.998000 * YOW);
  x[1] = y[1] / (28.014000 * YOW);
  x[2] = y[2] / (0.000549 * YOW);
  x[3] = y[3] / (15.999549 * YOW);
  x[4] = y[4] / (31.998549 * YOW);
  speciesEntropy(sor, tc);
  // Perform computation in Eq 42 and 43
  for (int i = 0; i < 5; i++) {
    result += x[i] * (sor[i] - log((x[i] + 1e-100)) - logPratio);
  }
  // Scale by R/W
  sbms = result * 8.31446261815324e+07 * YOW;
}

//  get temperature given internal energy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_EY(
  const amrex::Real e, const amrex::Real y[], amrex::Real& t, int& ierr)
{
#ifdef CONVERGENCE
  const int maxiter = 5000;
  const amrex::Real tol = 1.e-12;
#else
  const int maxiter = 200;
  const amrex::Real tol = 1.e-6;
#endif
  amrex::Real tmin = 90;   // max lower bound for thermo def
  amrex::Real tmax = 4000; // min upper bound for thermo def
  amrex::Real e1, emin, emax, cv, t1, dt;
  CKUBMS(tmin, y, emin);
  CKUBMS(tmax, y, emax);
  if (e < emin) {
    // Linear Extrapolation below tmin
    CKCVBS(tmin, y, cv);
    t = tmin - (emin - e) / cv;
    ierr = 1;
    return;
  }
  if (e > emax) {
    // Linear Extrapolation above tmax
    CKCVBS(tmax, y, cv);
    t = tmax - (emax - e) / cv;
    ierr = 1;
    return;
  }
  t1 = t;
  if (t1 < tmin || t1 > tmax) {
    t1 = tmin + (tmax - tmin) / (emax - emin) * (e - emin);
  }
  for (int i = 0; i < maxiter; ++i) {
    CKUBMS(t1, y, e1);
    CKCVBS(t1, y, cv);
    dt = (e - e1) / cv;
    if (dt > 100.) {
      dt = 100.;
    } else if (dt < -100.) {
      dt = -100.;
    } else if (fabs(dt) < tol) {
      break;
    }
    t1 += dt;
  }
  t = t1;
  ierr = 0;
}

//  get temperature given enthalpy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_HY(
  const amrex::Real h, const amrex::Real y[], amrex::Real& t, int& ierr)
{
#ifdef CONVERGENCE
  const int maxiter = 5000;
  const amrex::Real tol = 1.e-12;
#else
  const int maxiter = 200;
  const amrex::Real tol = 1.e-6;
#endif
  amrex::Real tmin = 90;   // max lower bound for thermo def
  amrex::Real tmax = 4000; // min upper bound for thermo def
  amrex::Real h1, hmin, hmax, cp, t1, dt;
  CKHBMS(tmin, y, hmin);
  CKHBMS(tmax, y, hmax);
  if (h < hmin) {
    // Linear Extrapolation below tmin
    CKCPBS(tmin, y, cp);
    t = tmin - (hmin - h) / cp;
    ierr = 1;
    return;
  }
  if (h > hmax) {
    // Linear Extrapolation above tmax
    CKCPBS(tmax, y, cp);
    t = tmax - (hmax - h) / cp;
    ierr = 1;
    return;
  }
  t1 = t;
  if (t1 < tmin || t1 > tmax) {
    t1 = tmin + (tmax - tmin) / (hmax - hmin) * (h - hmin);
  }
  for (int i = 0; i < maxiter; ++i) {
    CKHBMS(t1, y, h1);
    CKCPBS(t1, y, cp);
    dt = (h - h1) / cp;
    if (dt > 100.) {
      dt = 100.;
    } else if (dt < -100.) {
      dt = -100.;
    } else if (fabs(dt) < tol) {
      break;
    }
    t1 += dt;
  }
  t = t1;
  ierr = 0;
}

// Compute P = rhoRT/W(x)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPX(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real& P)
{
  amrex::Real XW = 0;                      //  To hold mean molecular wt
  XW += x[0] * 31.998000;                  // O2
  XW += x[1] * 28.014000;                  // N2
  XW += x[2] * 0.000549;                   // E
  XW += x[3] * 15.999549;                  // On
  XW += x[4] * 31.998549;                  // O2n
  P = rho * 8.31446261815324e+07 * T / XW; // P = rho*R*T/W
}

// Compute P = rhoRT/W(y)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPY(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real& P)
{
  amrex::Real YOW = 0; //  for computing mean MW

  for (int i = 0; i < 5; i++) {
    YOW += y[i] * imw(i);
  }
  P = rho * 8.31446261815324e+07 * T * YOW; // P = rho*R*T/W
}

// Compute P = rhoRT/W(c)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPC(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real c[],
  amrex::Real& P)
{
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 31.998000; // O2
  W += c[1] * 28.014000; // N2
  W += c[2] * 0.000549;  // E
  W += c[3] * 15.999549; // On
  W += c[4] * 31.998549; // O2n

  for (int id = 0; id < 5; ++id) {
    sumC += c[id];
  }
  P = rho * 8.31446261815324e+07 * T * sumC / W; // P = rho*R*T/W
}

// Compute rho = PW(x)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOX(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real& rho)
{
  amrex::Real XW = 0;                        //  To hold mean molecular wt
  XW += x[0] * 31.998000;                    // O2
  XW += x[1] * 28.014000;                    // N2
  XW += x[2] * 0.000549;                     // E
  XW += x[3] * 15.999549;                    // On
  XW += x[4] * 31.998549;                    // O2n
  rho = P * XW / (8.31446261815324e+07 * T); // rho = P*W/(R*T)
}

// Compute rho = P*W(y)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOY(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real& rho)
{
  amrex::Real YOW = 0;

  for (int i = 0; i < 5; i++) {
    YOW += y[i] * imw(i);
  }

  rho = P / (8.31446261815324e+07 * T * YOW); // rho = P*W/(R*T)
}

// Compute rho = P*W(c)/(R*T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOC(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real c[],
  amrex::Real& rho)
{
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 31.998000; // O2
  W += c[1] * 28.014000; // N2
  W += c[2] * 0.000549;  // E
  W += c[3] * 15.999549; // On
  W += c[4] * 31.998549; // O2n

  for (int id = 0; id < 5; ++id) {
    sumC += c[id];
  }
  rho = P * W / (sumC * T * 8.31446261815324e+07); // rho = PW/(R*T)
}

// get molecular weight for all species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWT(amrex::Real wt[])
{
  get_mw(wt);
}

// given y[species]: mass fractions
// s mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWY(const amrex::Real y[], amrex::Real& wtm)
{
  amrex::Real YOW = 0;

  for (int i = 0; i < 5; i++) {
    YOW += y[i] * imw(i);
  }

  wtm = 1.0 / YOW;
}

// given x[species]: mole fractions
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWX(const amrex::Real x[], amrex::Real& wtm)
{
  amrex::Real XW = 0;     //  see Eq 4 in CK Manual
  XW += x[0] * 31.998000; // O2
  XW += x[1] * 28.014000; // N2
  XW += x[2] * 0.000549;  // E
  XW += x[3] * 15.999549; // On
  XW += x[4] * 31.998549; // O2n
  wtm = XW;
}

// given c[species]: molar concentration
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWC(const amrex::Real c[], amrex::Real& wtm)
{
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 31.998000; // O2
  W += c[1] * 28.014000; // N2
  W += c[2] * 0.000549;  // E
  W += c[3] * 15.999549; // On
  W += c[4] * 31.998549; // O2n

  for (int id = 0; id < 5; ++id) {
    sumC += c[id];
  }
  //  CK provides no guard against division by zero
  wtm = W / sumC;
}

// get Cp/R as a function of T
// for all species (Eq 19)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPOR(const amrex::Real T, amrex::Real cpor[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cp_R(cpor, tc);
}

// get H/RT as a function of T
// for all species (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHORT(const amrex::Real T, amrex::Real hort[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  speciesEnthalpy(hort, tc);
}

// get S/R as a function of T
// for all species (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSOR(const amrex::Real T, amrex::Real sor[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  speciesEntropy(sor, tc);
}

// convert y[species] (mass fracs) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTX(const amrex::Real y[], amrex::Real x[])
{
  amrex::Real YOW = 0;

  for (int i = 0; i < 5; i++) {
    YOW += y[i] * imw(i);
  }

  amrex::Real YOWINV = 1.0 / YOW;

  for (int i = 0; i < 5; i++) {
    x[i] = y[i] * imw(i) * YOWINV;
  }
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTCP(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real c[])
{
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 5; i++) {
    c[i] = y[i] * imw(i);
  }
  for (int i = 0; i < 5; i++) {
    YOW += c[i];
  }

  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // Now compute conversion

  for (int i = 0; i < 5; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTCR(
  const amrex::Real rho,
  amrex::Real /*T*/,
  const amrex::Real y[],
  amrex::Real c[])
{

  for (int i = 0; i < 5; i++) {
    c[i] = rho * y[i] * imw(i);
  }
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTY(const amrex::Real x[], amrex::Real y[])
{
  amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
  // Compute mean molecular wt first
  XW += x[0] * 31.998000; // O2
  XW += x[1] * 28.014000; // N2
  XW += x[2] * 0.000549;  // E
  XW += x[3] * 15.999549; // On
  XW += x[4] * 31.998549; // O2n
  // Now compute conversion
  amrex::Real XWinv = 1.0 / XW;
  y[0] = x[0] * 31.998000 * XWinv;
  y[1] = x[1] * 28.014000 * XWinv;
  y[2] = x[2] * 0.000549 * XWinv;
  y[3] = x[3] * 15.999549 * XWinv;
  y[4] = x[4] * 31.998549 * XWinv;
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTCP(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real c[])
{
  amrex::Real PORT = P / (8.31446261815324e+07 * T); // P/RT

  // Compute conversion, see Eq 10
  for (int id = 0; id < 5; ++id) {
    c[id] = x[id] * PORT;
  }
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTCR(
  const amrex::Real rho,
  const amrex::Real /*T*/,
  const amrex::Real x[],
  amrex::Real c[])
{
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 31.998000; // O2
  XW += x[1] * 28.014000; // N2
  XW += x[2] * 0.000549;  // E
  XW += x[3] * 15.999549; // On
  XW += x[4] * 31.998549; // O2n
  ROW = rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 5; ++id) {
    c[id] = x[id] * ROW;
  }
}

// convert c[species] (molar conc) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTX(const amrex::Real c[], amrex::Real x[])
{
  amrex::Real sumC = 0;

  // compute sum of c
  for (int id = 0; id < 5; ++id) {
    sumC += c[id];
  }

  //  See Eq 13
  amrex::Real sumCinv = 1.0 / sumC;
  for (int id = 0; id < 5; ++id) {
    x[id] = c[id] * sumCinv;
  }
}

// convert c[species] (molar conc) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTY(const amrex::Real c[], amrex::Real y[])
{
  amrex::Real CW = 0; // See Eq 12 in CK Manual
  // compute denominator in eq 12 first
  CW += c[0] * 31.998000; // O2
  CW += c[1] * 28.014000; // N2
  CW += c[2] * 0.000549;  // E
  CW += c[3] * 15.999549; // On
  CW += c[4] * 31.998549; // O2n
  // Now compute conversion
  amrex::Real CWinv = 1.0 / CW;
  y[0] = c[0] * 31.998000 * CWinv;
  y[1] = c[1] * 28.014000 * CWinv;
  y[2] = c[2] * 0.000549 * CWinv;
  y[3] = c[3] * 15.999549 * CWinv;
  y[4] = c[4] * 31.998549 * CWinv;
}

// get specific heat at constant volume as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVML(const amrex::Real T, amrex::Real cvml[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cv_R(cvml, tc);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    cvml[id] *= 8.31446261815324e+07;
  }
}

// get specific heat at constant pressure as a
// function of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPML(const amrex::Real T, amrex::Real cpml[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cp_R(cpml, tc);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    cpml[id] *= 8.31446261815324e+07;
  }
}

// get internal energy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUML(const amrex::Real T, amrex::Real uml[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesInternalEnergy(uml, tc);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    uml[id] *= RT;
  }
}

// get enthalpy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHML(const amrex::Real T, amrex::Real hml[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesEnthalpy(hml, tc);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    hml[id] *= RT;
  }
}

// Returns the standard-state entropies in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSML(const amrex::Real T, amrex::Real sml[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  speciesEntropy(sml, tc);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    sml[id] *= 8.31446261815324e+07;
  }
}

// Returns the specific heats at constant volume
// in mass units (Eq. 29)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVMS(const amrex::Real T, amrex::Real cvms[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cv_R(cvms, tc);
  // multiply by R/molecularweight
  cvms[0] *= 2.598431970171023e+06; // O2
  cvms[1] *= 2.967966951578939e+06; // N2
  cvms[2] *= 1.515633818095350e+11; // E
  cvms[3] *= 5.196685754368074e+06; // On
  cvms[4] *= 2.598387422913829e+06; // O2n
}

// Returns the specific heats at constant pressure
// in mass units (Eq. 26)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPMS(const amrex::Real T, amrex::Real cpms[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cp_R(cpms, tc);
  // multiply by R/molecularweight
  cpms[0] *= 2.598431970171023e+06; // O2
  cpms[1] *= 2.967966951578939e+06; // N2
  cpms[2] *= 1.515633818095350e+11; // E
  cpms[3] *= 5.196685754368074e+06; // On
  cpms[4] *= 2.598387422913829e+06; // O2n
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUMS(const amrex::Real T, amrex::Real ums[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T

  speciesInternalEnergy(ums, tc);

  for (int i = 0; i < 5; i++) {
    ums[i] *= RT * imw(i);
  }
}

// Returns enthalpy in mass units (Eq 27.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHMS(const amrex::Real T, amrex::Real hms[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T

  speciesEnthalpy(hms, tc);

  for (int i = 0; i < 5; i++) {
    hms[i] *= RT * imw(i);
  }
}

// Returns the entropies in mass units (Eq 28.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSMS(const amrex::Real T, amrex::Real sms[])
{
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  speciesEntropy(sms, tc);
  // multiply by R/molecularweight
  sms[0] *= 2.598431970171023e+06; // O2
  sms[1] *= 2.967966951578939e+06; // N2
  sms[2] *= 1.515633818095350e+11; // E
  sms[3] *= 5.196685754368074e+06; // On
  sms[4] *= 2.598387422913829e+06; // O2n
}

// GPU version of productionRate: no more use of thermo namespace vectors
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_qfqr(
  amrex::Real* /*qf*/,
  amrex::Real* /*qr*/,
  const amrex::Real* /*sc*/,
  const amrex::Real* /*sc_qss*/,
  const amrex::Real* /*tc*/,
  const amrex::Real /*invT*/)
{
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
productionRate(
  amrex::Real* wdot, const amrex::Real* /*sc*/, const amrex::Real /*T*/)
{

  for (int i = 0; i < 5; ++i) {
    wdot[i] = 0.0;
  }
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[])
{

  // convert to SI
  for (int id = 0; id < 5; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    C[id] *= 1.0e-6;
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWYP(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real wdot[])
{
  amrex::Real c[5]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 5; i++) {
    YOW += y[i] * imw(i);
  }
  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // multiply by 1e6 so c goes to SI
  PWORT *= 1e6;
  // Now compute conversion (and go to SI)
  for (int i = 0; i < 5; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWXP(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real wdot[])
{
  amrex::Real c[5]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 5; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWYR(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real wdot[])
{
  amrex::Real c[5]; // temporary storage

  // See Eq 8 with an extra 1e6 so c goes to SI
  for (int i = 0; i < 5; i++) {
    c[i] = 1e6 * rho * y[i] * imw(i);
  }

  // call productionRate
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWXR(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real wdot[])
{
  amrex::Real c[5];   // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 31.998000; // O2
  XW += x[1] * 28.014000; // N2
  XW += x[2] * 0.000549;  // E
  XW += x[3] * 15.999549; // On
  XW += x[4] * 31.998549; // O2n
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 5; ++id) {
    c[id] = x[id] * ROW;
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int kcharge[])
{
  kcharge[0] = 0;  // O2
  kcharge[1] = 0;  // N2
  kcharge[2] = -1; // E
  kcharge[3] = -1; // On
  kcharge[4] = -1; // O2n
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real zk[])
{

  int kchrg[5];
  CKCHRG(kchrg);

  for (int id = 0; id < 5; ++id) {
    zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
  }
}

// compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dcvpRdT(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: O2
    species[0] = +1.12748635e-03 - 1.15123009e-06 * tc[1] +
                 3.94163169e-09 * tc[2] - 3.50742157e-12 * tc[3];
    // species 1: N2
    species[1] = +1.40824000e-03 - 7.92644400e-06 * tc[1] +
                 1.69245450e-08 * tc[2] - 9.77942000e-12 * tc[3];
    // species 3: On
    species[3] = -1.69804907e-03 + 5.96139910e-06 * tc[1] -
                 7.31505381e-09 * tc[2] + 3.04491724e-12 * tc[3];
    // species 4: O2n
    species[4] = -9.28741138e-04 + 1.29095416e-05 * tc[1] -
                 2.32411014e-08 * tc[2] + 1.17333065e-11 * tc[3];
  } else {
    // species 0: O2
    species[0] = +6.13519689e-04 - 2.51768398e-07 * tc[1] +
                 5.32584444e-11 * tc[2] - 4.54574124e-15 * tc[3];
    // species 1: N2
    species[1] = +1.48797700e-03 - 1.13695220e-06 * tc[1] +
                 3.02911200e-10 * tc[2] - 2.70134040e-14 * tc[3];
    // species 3: On
    species[3] = -4.66695513e-05 + 3.69824714e-08 * tc[1] -
                 9.54477669e-12 * tc[2] + 7.95851824e-16 * tc[3];
    // species 4: O2n
    species[4] = +5.98141823e-04 - 4.24267810e-07 * tc[1] +
                 1.08980274e-10 * tc[2] - 8.99956912e-15 * tc[3];
  }

  // species with no change at a midpoint T
  // species 2: E
  species[2] = +0.00000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3];
}

// compute an approx to the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian_precond(
  amrex::Real* J,
  const amrex::Real* /*sc*/,
  const amrex::Real /*T*/,
  const int /*HP*/)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  for (int i = 0; i < 36; i++) {
    J[i] = 0.0;
  }
#else
  amrex::Abort();
#endif
}

// compute an approx to the reaction Jacobian (for preconditioning)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT_SIMPLIFIED(
  amrex::Real* J, const amrex::Real* sc, const amrex::Real* Tp, const int* HP)
{
  amrex::Real c[5];

  for (int k = 0; k < 5; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian_precond(J, c, *Tp, *HP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 5; k++) {
    J[30 + k] *= 1.e-6;
    J[k * 6 + 5] *= 1.e6;
  }
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian(
  amrex::Real* J,
  const amrex::Real* /*sc*/,
  const amrex::Real /*T*/,
  const int /*consP*/)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  for (int i = 0; i < 36; i++) {
    J[i] = 0.0;
  }
#else
  amrex::Abort();
#endif
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT(
  amrex::Real* J,
  const amrex::Real* sc,
  const amrex::Real* Tp,
  const int* consP)
{
  amrex::Real c[5];

  for (int k = 0; k < 5; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian(J, c, *Tp, *consP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 5; k++) {
    J[30 + k] *= 1.e-6;
    J[k * 6 + 5] *= 1.e6;
  }
}

// Transport function declarations

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENIMC(int* LENIMC)
{
  *LENIMC = 21;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENRMC(int* LENRMC)
{
  *LENRMC = 630;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNO(int* NO)
{
  *NO = 4;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKK(int* KK)
{
  *KK = 5;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLITE(int* NLITE)
{
  *NLITE = 1;
}

// Patm in ergs/cm3
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPATM(amrex::Real* PATM)
{
  *PATM = 0.1013250000000000E+07;
}

// the molecular weights in g/mol
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetWT(amrex::Real* WT)
{
  WT[0] = 3.19980000E+01;
  WT[1] = 2.80140000E+01;
  WT[2] = 5.48579909E-04;
  WT[3] = 1.59995486E+01;
  WT[4] = 3.19985486E+01;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetEPS(amrex::Real* EPS)
{
  EPS[0] = 1.07400000E+02;
  EPS[1] = 9.75300000E+01;
  EPS[2] = 8.50000000E+02;
  EPS[3] = 8.00000000E+01;
  EPS[4] = 1.07400000E+02;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetSIG(amrex::Real* SIG)
{
  SIG[0] = 3.45800000E+00;
  SIG[1] = 3.62100000E+00;
  SIG[2] = 4.25000000E+02;
  SIG[3] = 2.75000000E+00;
  SIG[4] = 3.45800000E+00;
}

// the dipole moment in Debye
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetDIP(amrex::Real* DIP)
{
  DIP[0] = 0.00000000E+00;
  DIP[1] = 0.00000000E+00;
  DIP[2] = 0.00000000E+00;
  DIP[3] = 0.00000000E+00;
  DIP[4] = 0.00000000E+00;
}

// the polarizability in cubic Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPOL(amrex::Real* POL)
{
  POL[0] = 1.60000000E+00;
  POL[1] = 1.76000000E+00;
  POL[2] = 0.00000000E+00;
  POL[3] = 0.00000000E+00;
  POL[4] = 1.60000000E+00;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetZROT(amrex::Real* ZROT)
{
  ZROT[0] = 3.80000000E+00;
  ZROT[1] = 4.00000000E+00;
  ZROT[2] = 1.00000000E+00;
  ZROT[3] = 0.00000000E+00;
  ZROT[4] = 3.80000000E+00;
}

// 0: monoatomic, 1: linear, 2: nonlinear
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLIN(int* NLIN)
{
  NLIN[0] = 1;
  NLIN[1] = 1;
  NLIN[2] = 0;
  NLIN[3] = 0;
  NLIN[4] = 1;
}

// Poly fits for the viscosities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFETA(amrex::Real* COFETA)
{
  COFETA[0] = -1.60066453E+01;
  COFETA[1] = 2.16753735E+00;
  COFETA[2] = -1.97226850E-01;
  COFETA[3] = 8.50065468E-03;
  COFETA[4] = -1.55270224E+01;
  COFETA[5] = 1.92766908E+00;
  COFETA[6] = -1.66518287E-01;
  COFETA[7] = 7.19100649E-03;
  COFETA[8] = -2.69531463E+01;
  COFETA[9] = -7.88142336E-01;
  COFETA[10] = 3.18038402E-01;
  COFETA[11] = -1.85151697E-02;
  COFETA[12] = -1.40007552E+01;
  COFETA[13] = 1.44311220E+00;
  COFETA[14] = -1.02767741E-01;
  COFETA[15] = 4.39662983E-03;
  COFETA[16] = -1.60066368E+01;
  COFETA[17] = 2.16753735E+00;
  COFETA[18] = -1.97226850E-01;
  COFETA[19] = 8.50065468E-03;
}

// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFLAM(amrex::Real* COFLAM)
{
  COFLAM[0] = -2.11868696E+00;
  COFLAM[1] = 2.98568648E+00;
  COFLAM[2] = -2.86879126E-01;
  COFLAM[3] = 1.23850874E-02;
  COFLAM[4] = 7.60995874E+00;
  COFLAM[5] = -1.18418483E+00;
  COFLAM[6] = 3.03558391E-01;
  COFLAM[7] = -1.54159450E-02;
  COFLAM[8] = 1.12880184E-01;
  COFLAM[9] = -7.88142336E-01;
  COFLAM[10] = 3.18038402E-01;
  COFLAM[11] = -1.85151697E-02;
  COFLAM[12] = 2.78453320E+00;
  COFLAM[13] = 1.44311220E+00;
  COFLAM[14] = -1.02767741E-01;
  COFLAM[15] = 4.39662983E-03;
  COFLAM[16] = -4.17000350E+00;
  COFLAM[17] = 3.85409511E+00;
  COFLAM[18] = -3.99415041E-01;
  COFLAM[19] = 1.68318738E-02;
}

// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFD(amrex::Real* COFD)
{
  COFD[0] = -1.47079535E+01;
  COFD[1] = 3.10657376E+00;
  COFD[2] = -1.85922460E-01;
  COFD[3] = 7.92680827E-03;
  COFD[4] = -1.44285962E+01;
  COFD[5] = 2.99858376E+00;
  COFD[6] = -1.72232643E-01;
  COFD[7] = 7.34804765E-03;
  COFD[8] = -2.43418414E+01;
  COFD[9] = 5.33526425E+00;
  COFD[10] = -4.46630688E-01;
  COFD[11] = 1.81086658E-02;
  COFD[12] = -1.34772625E+01;
  COFD[13] = 2.80526779E+00;
  COFD[14] = -1.47524446E-01;
  COFD[15] = 6.29376912E-03;
  COFD[16] = -1.47079578E+01;
  COFD[17] = 3.10657376E+00;
  COFD[18] = -1.85922460E-01;
  COFD[19] = 7.92680827E-03;
  COFD[20] = -1.44285962E+01;
  COFD[21] = 2.99858376E+00;
  COFD[22] = -1.72232643E-01;
  COFD[23] = 7.34804765E-03;
  COFD[24] = -1.42056777E+01;
  COFD[25] = 2.91297621E+00;
  COFD[26] = -1.61544771E-01;
  COFD[27] = 6.90271324E-03;
  COFD[28] = -2.41515325E+01;
  COFD[29] = 5.28689054E+00;
  COFD[30] = -4.42932076E-01;
  COFD[31] = 1.80445515E-02;
  COFD[32] = -1.32254691E+01;
  COFD[33] = 2.69833584E+00;
  COFD[34] = -1.33541002E-01;
  COFD[35] = 5.68346517E-03;
  COFD[36] = -1.44286002E+01;
  COFD[37] = 2.99858376E+00;
  COFD[38] = -1.72232643E-01;
  COFD[39] = 7.34804765E-03;
  COFD[40] = -2.43418414E+01;
  COFD[41] = 5.33526425E+00;
  COFD[42] = -4.46630688E-01;
  COFD[43] = 1.81086658E-02;
  COFD[44] = -2.41515325E+01;
  COFD[45] = 5.28689054E+00;
  COFD[46] = -4.42932076E-01;
  COFD[47] = 1.80445515E-02;
  COFD[48] = -1.82216979E+01;
  COFD[49] = 1.68992430E+00;
  COFD[50] = 1.09701574E-01;
  COFD[51] = -8.81888811E-03;
  COFD[52] = -2.36348036E+01;
  COFD[53] = 5.13665150E+00;
  COFD[54] = -4.28328574E-01;
  COFD[55] = 1.75969700E-02;
  COFD[56] = -2.43418414E+01;
  COFD[57] = 5.33526425E+00;
  COFD[58] = -4.46630688E-01;
  COFD[59] = 1.81086658E-02;
  COFD[60] = -1.34772625E+01;
  COFD[61] = 2.80526779E+00;
  COFD[62] = -1.47524446E-01;
  COFD[63] = 6.29376912E-03;
  COFD[64] = -1.32254691E+01;
  COFD[65] = 2.69833584E+00;
  COFD[66] = -1.33541002E-01;
  COFD[67] = 5.68346517E-03;
  COFD[68] = -2.36348036E+01;
  COFD[69] = 5.13665150E+00;
  COFD[70] = -4.28328574E-01;
  COFD[71] = 1.75969700E-02;
  COFD[72] = -1.22676560E+01;
  COFD[73] = 2.49197560E+00;
  COFD[74] = -1.06603364E-01;
  COFD[75] = 4.50971238E-03;
  COFD[76] = -1.34772654E+01;
  COFD[77] = 2.80526779E+00;
  COFD[78] = -1.47524446E-01;
  COFD[79] = 6.29376912E-03;
  COFD[80] = -1.47079578E+01;
  COFD[81] = 3.10657376E+00;
  COFD[82] = -1.85922460E-01;
  COFD[83] = 7.92680827E-03;
  COFD[84] = -1.44286002E+01;
  COFD[85] = 2.99858376E+00;
  COFD[86] = -1.72232643E-01;
  COFD[87] = 7.34804765E-03;
  COFD[88] = -2.43418414E+01;
  COFD[89] = 5.33526425E+00;
  COFD[90] = -4.46630688E-01;
  COFD[91] = 1.81086658E-02;
  COFD[92] = -1.34772654E+01;
  COFD[93] = 2.80526779E+00;
  COFD[94] = -1.47524446E-01;
  COFD[95] = 6.29376912E-03;
  COFD[96] = -1.47079620E+01;
  COFD[97] = 3.10657376E+00;
  COFD[98] = -1.85922460E-01;
  COFD[99] = 7.92680827E-03;
}

// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKTDIF(int* KTDIF)
{
  KTDIF[0] = 2;
}

// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFTD(amrex::Real* COFTD)
{
  COFTD[0] = -8.18816845E-02;
  COFTD[1] = 6.25102118E-04;
  COFTD[2] = -1.93882483E-07;
  COFTD[3] = 1.90693525E-11;
  COFTD[4] = -6.27142421E-02;
  COFTD[5] = 6.17357859E-04;
  COFTD[6] = -1.93698208E-07;
  COFTD[7] = 1.91874330E-11;
  COFTD[8] = 0.00000000E+00;
  COFTD[9] = 0.00000000E+00;
  COFTD[10] = 0.00000000E+00;
  COFTD[11] = 0.00000000E+00;
  COFTD[12] = -2.08077509E-02;
  COFTD[13] = 5.95806873E-04;
  COFTD[14] = -1.91160904E-07;
  COFTD[15] = 1.92016568E-11;
  COFTD[16] = -8.18816846E-02;
  COFTD[17] = 6.25102119E-04;
  COFTD[18] = -1.93882484E-07;
  COFTD[19] = 1.90693525E-11;
}

// compute the critical parameters for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS(
  amrex::Real* Tci, amrex::Real* ai, amrex::Real* bi, amrex::Real* acentric_i)
{

  amrex::Real EPS[5];
  amrex::Real SIG[5];
  amrex::Real wt[5];
  amrex::Real Rcst = 83.144598; // in bar [CGS] !
  amrex::Real avogadro = 6.02214199e23;
  amrex::Real boltzmann = 1.3806503e-16; // we work in CGS

  egtransetEPS(EPS);
  egtransetSIG(SIG);
  get_mw(wt);

  // species 0: O2
  // Imported from NIST
  Tci[0] = 154.581000;
  ai[0] = 1e6 * 0.42748 * Rcst * Rcst * Tci[0] * Tci[0] /
          (31.998800 * 31.998800 * 50.430466);
  bi[0] = 0.08664 * Rcst * Tci[0] / (31.998800 * 50.430466);
  acentric_i[0] = 0.022200;

  // species 1: N2
  // Imported from NIST
  Tci[1] = 126.192000;
  ai[1] = 1e6 * 0.42748 * Rcst * Rcst * Tci[1] * Tci[1] /
          (28.013400 * 28.013400 * 33.958000);
  bi[1] = 0.08664 * Rcst * Tci[1] / (28.013400 * 33.958000);
  acentric_i[1] = 0.037200;

  // species 2: E
  Tci[2] = 1.316 * EPS[2];
  ai[2] = (5.55 * avogadro * avogadro * EPS[2] * boltzmann * 1e-24 * SIG[2] *
           SIG[2] * SIG[2]) /
          (wt[2] * wt[2]);
  bi[2] = 0.855 * avogadro * 1e-24 * SIG[2] * SIG[2] * SIG[2] / (wt[2]);
  acentric_i[2] = 0.0;

  // species 3: On
  Tci[3] = 1.316 * EPS[3];
  ai[3] = (5.55 * avogadro * avogadro * EPS[3] * boltzmann * 1e-24 * SIG[3] *
           SIG[3] * SIG[3]) /
          (wt[3] * wt[3]);
  bi[3] = 0.855 * avogadro * 1e-24 * SIG[3] * SIG[3] * SIG[3] / (wt[3]);
  acentric_i[3] = 0.0;

  // species 4: O2n
  Tci[4] = 1.316 * EPS[4];
  ai[4] = (5.55 * avogadro * avogadro * EPS[4] * boltzmann * 1e-24 * SIG[4] *
           SIG[4] * SIG[4]) /
          (wt[4] * wt[4]);
  bi[4] = 0.855 * avogadro * 1e-24 * SIG[4] * SIG[4] * SIG[4] / (wt[4]);
  acentric_i[4] = 0.0;
}

// compute the critical parameter quantities for each species for SRK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS_SRK(
  amrex::Real* sqrtOneOverTc,
  amrex::Real* sqrtAsti,
  amrex::Real* Bi,
  amrex::Real* Fomega)
{

  // species 0: O2
  // Imported from NIST
  sqrtOneOverTc[0] = 8.0430717653241e-02;
  sqrtAsti[0] = 3.6980081924740e+04;
  Bi[0] = 6.9005216444999e-01;
  Fomega[0] = 5.1945301904908e-01;

  // species 1: N2
  // Imported from NIST
  sqrtOneOverTc[1] = 8.9019282240563e-02;
  sqrtAsti[1] = 4.2022983929562e+04;
  Bi[1] = 9.5560052221563e-01;
  Fomega[1] = 5.4259343186608e-01;

  // species 2: E
  // Computed from Lennard-Jones
  sqrtOneOverTc[2] = 2.9899408206870e-02;
  sqrtAsti[2] = 7.7623333318521e+12;
  Bi[2] = 7.2051661156636e+10;
  Fomega[2] = 4.8508000000000e-01;

  // species 3: On
  // Computed from Lennard-Jones
  sqrtOneOverTc[3] = 9.7460146948181e-02;
  sqrtAsti[3] = 4.2498698800541e+04;
  Bi[3] = 6.6927940000755e-01;
  Fomega[3] = 4.8508000000000e-01;

  // species 4: O2n
  // Computed from Lennard-Jones
  sqrtOneOverTc[4] = 8.4114315417892e-02;
  sqrtAsti[4] = 3.4717518285005e+04;
  Bi[4] = 6.6536782042076e-01;
  Fomega[4] = 4.8508000000000e-01;
}

// gauss-jordan solver external routine
// Replace this routine with the one generated by the Gauss Jordan solver of DW
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
sgjsolve(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/)
{
  amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

// Replace this routine with the one generated by the Gauss Jordan solver of DW
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
sgjsolve_simplified(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/)
{
  amrex::Abort(
    "sgjsolve_simplified not implemented, choose a different solver ");
}

#endif
