#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  C
1  O
*/

// Species
#define CO_ID 0
#define CO2_ID 1

#define NUM_ELEMENTS 2
#define NUM_SPECIES 2
#define NUM_IONS 0
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
  mm = 2;
  kk = 2;
  ii = 0;
  nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[2] = {
  0.0357015351660121, // CO
  0.0227226249176305, // CO2
};
#endif
const amrex::Real h_global_imw[2] = {
  0.0357015351660121, // CO
  0.0227226249176305, // CO2
};

//  molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[2] = {
  28.010000, // CO
  44.009000, // CO2
};
#endif
const amrex::Real h_global_mw[2] = {
  28.010000, // CO
  44.009000, // CO2
};

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_imw(amrex::Real* imw_new)
{
  imw_new[0] = 0.0357015351660121; // CO
  imw_new[1] = 0.0227226249176305; // CO2
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
  mw_new[0] = 28.010000; // CO
  mw_new[1] = 44.009000; // CO2
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
cv_R(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;

  // species with inflection points at T = (1000.0, 6000.0) kelvin
  if (T < 1000) {
    // species 0: CO
    species[0] = +1.48904533e+04 * invT2 - 2.92228594e+02 * invT +
                 4.72452717e+00 - 8.17623503e-03 * T + 1.45690347e-05 * T2 -
                 1.08774630e-08 * T3 + 3.02794183e-12 * T4;
    // species 1: CO2
    species[1] = +4.94365054e+04 * invT2 - 6.26411601e+02 * invT +
                 4.30172524e+00 + 2.50381382e-03 * T - 2.12730873e-07 * T2 -
                 7.68998878e-10 * T3 + 2.84967780e-13 * T4;
  } else if ((1000 <= T) && (T < 6000)) {
    // species 0: CO
    species[0] = +4.61919725e+05 * invT2 - 1.94470486e+03 * invT +
                 4.91671418e+00 - 5.66428283e-04 * T + 1.39881454e-07 * T2 -
                 1.78768036e-11 * T3 + 9.62093557e-16 * T4;
    // species 1: CO2
    species[1] = +1.17696242e+05 * invT2 - 1.78879148e+03 * invT +
                 7.29152319e+00 - 9.22315678e-05 * T + 4.86367688e-09 * T2 -
                 1.89105331e-12 * T3 + 6.33003659e-16 * T4;
  } else {
    // species 0: CO
    species[0] = +8.86866296e+08 * invT2 - 7.50037784e+05 * invT +
                 2.48547498e+02 - 3.95635110e-02 * T + 3.29777208e-06 * T2 -
                 1.31840993e-10 * T3 + 1.99893795e-15 * T4;
    // species 1: CO2
    species[1] = -1.54442329e+09 * invT2 + 1.01684706e+06 * invT -
                 2.57140523e+02 + 3.36940108e-02 * T - 2.18118434e-06 * T2 +
                 6.99142084e-11 * T3 - 8.84235150e-16 * T4;
  }
}

// compute Cp/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cp_R(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;

  // species with inflection points at T = (1000.0, 6000.0) kelvin
  if (T < 1000) {
    // species 0: CO
    species[0] = +1.48904533e+04 * invT2 - 2.92228594e+02 * invT +
                 5.72452717e+00 - 8.17623503e-03 * T + 1.45690347e-05 * T2 -
                 1.08774630e-08 * T3 + 3.02794183e-12 * T4;
    // species 1: CO2
    species[1] = +4.94365054e+04 * invT2 - 6.26411601e+02 * invT +
                 5.30172524e+00 + 2.50381382e-03 * T - 2.12730873e-07 * T2 -
                 7.68998878e-10 * T3 + 2.84967780e-13 * T4;
  } else if ((1000 <= T) && (T < 6000)) {
    // species 0: CO
    species[0] = +4.61919725e+05 * invT2 - 1.94470486e+03 * invT +
                 5.91671418e+00 - 5.66428283e-04 * T + 1.39881454e-07 * T2 -
                 1.78768036e-11 * T3 + 9.62093557e-16 * T4;
    // species 1: CO2
    species[1] = +1.17696242e+05 * invT2 - 1.78879148e+03 * invT +
                 8.29152319e+00 - 9.22315678e-05 * T + 4.86367688e-09 * T2 -
                 1.89105331e-12 * T3 + 6.33003659e-16 * T4;
  } else {
    // species 0: CO
    species[0] = +8.86866296e+08 * invT2 - 7.50037784e+05 * invT +
                 2.49547498e+02 - 3.95635110e-02 * T + 3.29777208e-06 * T2 -
                 1.31840993e-10 * T3 + 1.99893795e-15 * T4;
    // species 1: CO2
    species[1] = -1.54442329e+09 * invT2 + 1.01684706e+06 * invT -
                 2.56140523e+02 + 3.36940108e-02 * T - 2.18118434e-06 * T2 +
                 6.99142084e-11 * T3 - 8.84235150e-16 * T4;
  }
}

// compute the g/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
gibbs(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;
  const amrex::Real logT = log(T);

  // species with inflection points at T = (1000.0, 6000.0) kelvin
  if (T < 1000) {
    // species 0: CO
    species[0] = -7.445226630000000e+03 * invT2 - 1.332354737390000e+04 * invT -
                 2.922285939000000e+02 * logT * invT -
                 5.724527170000000e+00 * logT + 1.358376852000000e+01 +
                 4.088117515000000e-03 * T - 2.428172448333333e-06 * T2 +
                 9.064552516666667e-10 * T3 - 1.513970913500000e-13 * T4;
    // species 1: CO2
    species[1] = -2.471825270000000e+04 * invT2 - 4.590839620100000e+04 * invT -
                 6.264116010000000e+02 * logT * invT -
                 5.301725240000000e+00 * logT + 1.235000468000000e+01 -
                 1.251906908000000e-03 * T + 3.545514546666667e-08 * T2 +
                 6.408323983333334e-11 * T3 - 1.424838900500000e-14 * T4;
  } else if ((1000 <= T) && (T < 6000)) {
    // species 0: CO
    species[0] = -2.309598625000000e+05 * invT2 - 4.410965947000001e+03 * invT -
                 1.944704863000000e+03 * logT * invT -
                 5.916714180000000e+00 * logT + 1.979084526000000e+01 +
                 2.832141415000000e-04 * T - 2.331357566666667e-08 * T2 +
                 1.489733634166667e-12 * T3 - 4.810467785000000e-17 * T4;
    // species 1: CO2
    species[1] = -5.884812095000000e+04 * invT2 - 4.087229737700000e+04 * invT -
                 1.788791477000000e+03 * logT * invT -
                 8.291523189999999e+00 * logT + 3.481821600000000e+01 +
                 4.611578390000000e-05 * T - 8.106128133333333e-10 * T2 +
                 1.575877760000000e-13 * T3 - 3.165018295000000e-17 * T4;
  } else {
    // species 0: CO
    species[0] = -4.434331480000000e+08 * invT2 + 4.951383346000000e+06 * invT -
                 7.500377840000000e+05 * logT * invT -
                 2.495474979000000e+02 * logT + 2.310252283900000e+03 +
                 1.978175550000000e-02 * T - 5.496286800000000e-07 * T2 +
                 1.098674944166667e-11 * T3 - 9.994689739999999e-17 * T4;
    // species 1: CO2
    species[1] = +7.722116435000000e+08 * invT2 - 7.026367454000000e+06 * invT +
                 1.016847056000000e+06 * logT * invT +
                 2.561405230000000e+02 * logT - 2.510318016000000e+03 -
                 1.684700540000000e-02 * T + 3.635307228333333e-07 * T2 -
                 5.826184033333333e-12 * T3 + 4.421175750000000e-17 * T4;
  }
}

// compute the a/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
helmholtz(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;
  const amrex::Real logT = log(T);

  // species with inflection points at T = (1000.0, 6000.0) kelvin
  if (T < 1000) {
    // species 0: CO
    species[0] = -7.44522663e+03 * invT2 - 1.33235474e+04 * invT -
                 2.92228594e+02 * logT * invT - 5.72452717e+00 * logT +
                 1.25837685e+01 + 4.08811751e-03 * T - 2.42817245e-06 * T2 +
                 9.06455252e-10 * T3 - 1.51397091e-13 * T4;
    // species 1: CO2
    species[1] = -2.47182527e+04 * invT2 - 4.59083962e+04 * invT -
                 6.26411601e+02 * logT * invT - 5.30172524e+00 * logT +
                 1.13500047e+01 - 1.25190691e-03 * T + 3.54551455e-08 * T2 +
                 6.40832398e-11 * T3 - 1.42483890e-14 * T4;
  } else if ((1000 <= T) && (T < 6000)) {
    // species 0: CO
    species[0] = -2.30959862e+05 * invT2 - 4.41096595e+03 * invT -
                 1.94470486e+03 * logT * invT - 5.91671418e+00 * logT +
                 1.87908453e+01 + 2.83214141e-04 * T - 2.33135757e-08 * T2 +
                 1.48973363e-12 * T3 - 4.81046778e-17 * T4;
    // species 1: CO2
    species[1] = -5.88481209e+04 * invT2 - 4.08722974e+04 * invT -
                 1.78879148e+03 * logT * invT - 8.29152319e+00 * logT +
                 3.38182160e+01 + 4.61157839e-05 * T - 8.10612813e-10 * T2 +
                 1.57587776e-13 * T3 - 3.16501830e-17 * T4;
  } else {
    // species 0: CO
    species[0] = -4.43433148e+08 * invT2 + 4.95138335e+06 * invT -
                 7.50037784e+05 * logT * invT - 2.49547498e+02 * logT +
                 2.30925228e+03 + 1.97817555e-02 * T - 5.49628680e-07 * T2 +
                 1.09867494e-11 * T3 - 9.99468974e-17 * T4;
    // species 1: CO2
    species[1] = +7.72211644e+08 * invT2 - 7.02636745e+06 * invT +
                 1.01684706e+06 * logT * invT + 2.56140523e+02 * logT -
                 2.51131802e+03 - 1.68470054e-02 * T + 3.63530723e-07 * T2 -
                 5.82618403e-12 * T3 + 4.42117575e-17 * T4;
  }
}

// compute the e/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesInternalEnergy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;
  const amrex::Real logT = log(T);

  // species with inflection points at T = (1000.0, 6000.0) kelvin
  if (T < 1000) {
    // species 0: CO
    species[0] = -1.48904533e+04 * invT2 - 2.92228594e+02 * logT * invT +
                 4.72452717e+00 - 4.08811751e-03 * T + 4.85634490e-06 * T2 -
                 2.71936576e-09 * T3 + 6.05588365e-13 * T4 -
                 1.30313188e+04 * invT;
    // species 1: CO2
    species[1] = -4.94365054e+04 * invT2 - 6.26411601e+02 * logT * invT +
                 4.30172524e+00 + 1.25190691e-03 * T - 7.09102909e-08 * T2 -
                 1.92249720e-10 * T3 + 5.69935560e-14 * T4 -
                 4.52819846e+04 * invT;
  } else if ((1000 <= T) && (T < 6000)) {
    // species 0: CO
    species[0] = -4.61919725e+05 * invT2 - 1.94470486e+03 * logT * invT +
                 4.91671418e+00 - 2.83214141e-04 * T + 4.66271513e-08 * T2 -
                 4.46920090e-12 * T3 + 1.92418711e-16 * T4 -
                 2.46626108e+03 * invT;
    // species 1: CO2
    species[1] = -1.17696242e+05 * invT2 - 1.78879148e+03 * logT * invT +
                 7.29152319e+00 - 4.61157839e-05 * T + 1.62122563e-09 * T2 -
                 4.72763328e-13 * T3 + 1.26600732e-16 * T4 -
                 3.90835059e+04 * invT;
  } else {
    // species 0: CO
    species[0] = -8.86866296e+08 * invT2 - 7.50037784e+05 * logT * invT +
                 2.48547498e+02 - 1.97817555e-02 * T + 1.09925736e-06 * T2 -
                 3.29602483e-11 * T3 + 3.99787590e-16 * T4 +
                 5.70142113e+06 * invT;
    // species 1: CO2
    species[1] = +1.54442329e+09 * invT2 + 1.01684706e+06 * logT * invT -
                 2.57140523e+02 + 1.68470054e-02 * T - 7.27061446e-07 * T2 +
                 1.74785521e-11 * T3 - 1.76847030e-16 * T4 -
                 8.04321451e+06 * invT;
  }
}

// compute the h/(RT) at the given temperature (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEnthalpy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;
  const amrex::Real logT = log(T);

  // species with inflection points at T = (1000.0, 6000.0) kelvin
  if (T < 1000) {
    // species 0: CO
    species[0] = -1.48904533e+04 * invT2 - 2.92228594e+02 * logT * invT +
                 5.72452717e+00 - 4.08811751e-03 * T + 4.85634490e-06 * T2 -
                 2.71936576e-09 * T3 + 6.05588365e-13 * T4 -
                 1.30313188e+04 * invT;
    // species 1: CO2
    species[1] = -4.94365054e+04 * invT2 - 6.26411601e+02 * logT * invT +
                 5.30172524e+00 + 1.25190691e-03 * T - 7.09102909e-08 * T2 -
                 1.92249720e-10 * T3 + 5.69935560e-14 * T4 -
                 4.52819846e+04 * invT;
  } else if ((1000 <= T) && (T < 6000)) {
    // species 0: CO
    species[0] = -4.61919725e+05 * invT2 - 1.94470486e+03 * logT * invT +
                 5.91671418e+00 - 2.83214141e-04 * T + 4.66271513e-08 * T2 -
                 4.46920090e-12 * T3 + 1.92418711e-16 * T4 -
                 2.46626108e+03 * invT;
    // species 1: CO2
    species[1] = -1.17696242e+05 * invT2 - 1.78879148e+03 * logT * invT +
                 8.29152319e+00 - 4.61157839e-05 * T + 1.62122563e-09 * T2 -
                 4.72763328e-13 * T3 + 1.26600732e-16 * T4 -
                 3.90835059e+04 * invT;
  } else {
    // species 0: CO
    species[0] = -8.86866296e+08 * invT2 - 7.50037784e+05 * logT * invT +
                 2.49547498e+02 - 1.97817555e-02 * T + 1.09925736e-06 * T2 -
                 3.29602483e-11 * T3 + 3.99787590e-16 * T4 +
                 5.70142113e+06 * invT;
    // species 1: CO2
    species[1] = +1.54442329e+09 * invT2 + 1.01684706e+06 * logT * invT -
                 2.56140523e+02 + 1.68470054e-02 * T - 7.27061446e-07 * T2 +
                 1.74785521e-11 * T3 - 1.76847030e-16 * T4 -
                 8.04321451e+06 * invT;
  }
}

// compute the S/R at the given temperature (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEntropy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;
  const amrex::Real logT = log(T);

  // species with inflection points at T = (1000.0, 6000.0) kelvin
  if (T < 1000) {
    // species 0: CO
    species[0] = -7.44522663e+03 * invT2 + 2.92228594e+02 * invT +
                 5.72452717e+00 * logT - 8.17623503e-03 * T +
                 7.28451735e-06 * T2 - 3.62582101e-09 * T3 +
                 7.56985457e-13 * T4 - 7.85924135e+00;
    // species 1: CO2
    species[1] = -2.47182527e+04 * invT2 + 6.26411601e+02 * invT +
                 5.30172524e+00 * logT + 2.50381382e-03 * T -
                 1.06365436e-07 * T2 - 2.56332959e-10 * T3 +
                 7.12419450e-14 * T4 - 7.04827944e+00;
  } else if ((1000 <= T) && (T < 6000)) {
    // species 0: CO
    species[0] = -2.30959862e+05 * invT2 + 1.94470486e+03 * invT +
                 5.91671418e+00 * logT - 5.66428283e-04 * T +
                 6.99407270e-08 * T2 - 5.95893454e-12 * T3 +
                 2.40523389e-16 * T4 - 1.38741311e+01;
    // species 1: CO2
    species[1] = -5.88481209e+04 * invT2 + 1.78879148e+03 * invT +
                 8.29152319e+00 * logT - 9.22315678e-05 * T +
                 2.43183844e-09 * T2 - 6.30351104e-13 * T3 +
                 1.58250915e-16 * T4 - 2.65266928e+01;
  } else {
    // species 0: CO
    species[0] = -4.43433148e+08 * invT2 + 7.50037784e+05 * invT +
                 2.49547498e+02 * logT - 3.95635110e-02 * T +
                 1.64888604e-06 * T2 - 4.39469978e-11 * T3 +
                 4.99734487e-16 * T4 - 2.06070479e+03;
    // species 1: CO2
    species[1] = +7.72211644e+08 * invT2 - 1.01684706e+06 * invT -
                 2.56140523e+02 * logT + 3.36940108e-02 * T -
                 1.09059217e-06 * T2 + 2.33047361e-11 * T3 -
                 2.21058788e-16 * T4 + 2.25417749e+03;
  }
}

// compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dcvpRdT(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;
  const amrex::Real invT3 = invT * invT * invT;

  // species with inflection points at T = (1000.0, 6000.0) kelvin
  if (T < 1000) {
    // species 0: CO
    species[0] = -2.97809065e+04 * invT3 + 2.92228594e+02 * invT2 -
                 8.17623503e-03 + 2.91380694e-05 * T - 3.26323891e-08 * T2 +
                 1.21117673e-11 * T3;
    // species 1: CO2
    species[1] = -9.88730108e+04 * invT3 + 6.26411601e+02 * invT2 +
                 2.50381382e-03 - 4.25461746e-07 * T - 2.30699663e-09 * T2 +
                 1.13987112e-12 * T3;
  } else if ((1000 <= T) && (T < 6000)) {
    // species 0: CO
    species[0] = -9.23839450e+05 * invT3 + 1.94470486e+03 * invT2 -
                 5.66428283e-04 + 2.79762908e-07 * T - 5.36304108e-11 * T2 +
                 3.84837423e-15 * T3;
    // species 1: CO2
    species[1] = -2.35392484e+05 * invT3 + 1.78879148e+03 * invT2 -
                 9.22315678e-05 + 9.72735376e-09 * T - 5.67315994e-12 * T2 +
                 2.53201464e-15 * T3;
  } else {
    // species 0: CO
    species[0] = -1.77373259e+09 * invT3 + 7.50037784e+05 * invT2 -
                 3.95635110e-02 + 6.59554416e-06 * T - 3.95522980e-10 * T2 +
                 7.99575179e-15 * T3;
    // species 1: CO2
    species[1] = +3.08884657e+09 * invT3 - 1.01684706e+06 * invT2 +
                 3.36940108e-02 - 4.36236867e-06 * T + 2.09742625e-10 * T2 -
                 3.53694060e-15 * T3;
  }
}

// Returns the mean specific heat at CP (Eq. 33)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cpbl)
{
  amrex::Real result = 0;
  amrex::Real cpor[2]; //  temporary storage
  cp_R(cpor, T);

  // perform dot product
  for (int id = 0; id < 2; ++id) {
    result += x[id] * cpor[id];
  }

  cpbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CP (Eq. 34)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBS(const amrex::Real T, const amrex::Real y[], amrex::Real& cpbs)
{
  amrex::Real result = 0.0;

  // compute Cp/R at the given temperature
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;

  // species with inflection points at T = (1000.0, 6000.0) kelvin
  if (T < 1000) {
    // species 0: CO
    result += y[0] *
              (+1.48904533e+04 * invT2 - 2.92228594e+02 * invT +
               5.72452717e+00 - 8.17623503e-03 * T + 1.45690347e-05 * T2 -
               1.08774630e-08 * T3 + 3.02794183e-12 * T4) *
              0.0357015351660121;
    // species 1: CO2
    result += y[1] *
              (+4.94365054e+04 * invT2 - 6.26411601e+02 * invT +
               5.30172524e+00 + 2.50381382e-03 * T - 2.12730873e-07 * T2 -
               7.68998878e-10 * T3 + 2.84967780e-13 * T4) *
              0.0227226249176305;
  } else if ((1000 <= T) && (T < 6000)) {
    // species 0: CO
    result += y[0] *
              (+4.61919725e+05 * invT2 - 1.94470486e+03 * invT +
               5.91671418e+00 - 5.66428283e-04 * T + 1.39881454e-07 * T2 -
               1.78768036e-11 * T3 + 9.62093557e-16 * T4) *
              0.0357015351660121;
    // species 1: CO2
    result += y[1] *
              (+1.17696242e+05 * invT2 - 1.78879148e+03 * invT +
               8.29152319e+00 - 9.22315678e-05 * T + 4.86367688e-09 * T2 -
               1.89105331e-12 * T3 + 6.33003659e-16 * T4) *
              0.0227226249176305;
  } else {
    // species 0: CO
    result += y[0] *
              (+8.86866296e+08 * invT2 - 7.50037784e+05 * invT +
               2.49547498e+02 - 3.95635110e-02 * T + 3.29777208e-06 * T2 -
               1.31840993e-10 * T3 + 1.99893795e-15 * T4) *
              0.0357015351660121;
    // species 1: CO2
    result += y[1] *
              (-1.54442329e+09 * invT2 + 1.01684706e+06 * invT -
               2.56140523e+02 + 3.36940108e-02 * T - 2.18118434e-06 * T2 +
               6.99142084e-11 * T3 - 8.84235150e-16 * T4) *
              0.0227226249176305;
  }

  cpbs = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 35)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cvbl)
{
  amrex::Real result = 0;
  amrex::Real cvor[2]; //  temporary storage
  cv_R(cvor, T);

  // perform dot product
  for (int id = 0; id < 2; ++id) {
    result += x[id] * cvor[id];
  }

  cvbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 36)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBS(const amrex::Real T, const amrex::Real y[], amrex::Real& cvbs)
{
  amrex::Real result = 0.0;
  // compute Cv/R at the given temperature
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;

  // species with inflection points at T = (1000.0, 6000.0) kelvin
  if (T < 1000) {
    // species 0: CO
    result += y[0] *
              (+1.48904533e+04 * invT2 - 2.92228594e+02 * invT +
               4.72452717e+00 - 8.17623503e-03 * T + 1.45690347e-05 * T2 -
               1.08774630e-08 * T3 + 3.02794183e-12 * T4) *
              0.0357015351660121;
    // species 1: CO2
    result += y[1] *
              (+4.94365054e+04 * invT2 - 6.26411601e+02 * invT +
               4.30172524e+00 + 2.50381382e-03 * T - 2.12730873e-07 * T2 -
               7.68998878e-10 * T3 + 2.84967780e-13 * T4) *
              0.0227226249176305;
  } else if ((1000 <= T) && (T < 6000)) {
    // species 0: CO
    result += y[0] *
              (+4.61919725e+05 * invT2 - 1.94470486e+03 * invT +
               4.91671418e+00 - 5.66428283e-04 * T + 1.39881454e-07 * T2 -
               1.78768036e-11 * T3 + 9.62093557e-16 * T4) *
              0.0357015351660121;
    // species 1: CO2
    result += y[1] *
              (+1.17696242e+05 * invT2 - 1.78879148e+03 * invT +
               7.29152319e+00 - 9.22315678e-05 * T + 4.86367688e-09 * T2 -
               1.89105331e-12 * T3 + 6.33003659e-16 * T4) *
              0.0227226249176305;
  } else {
    // species 0: CO
    result += y[0] *
              (+8.86866296e+08 * invT2 - 7.50037784e+05 * invT +
               2.48547498e+02 - 3.95635110e-02 * T + 3.29777208e-06 * T2 -
               1.31840993e-10 * T3 + 1.99893795e-15 * T4) *
              0.0357015351660121;
    // species 1: CO2
    result += y[1] *
              (-1.54442329e+09 * invT2 + 1.01684706e+06 * invT -
               2.57140523e+02 + 3.36940108e-02 * T - 2.18118434e-06 * T2 +
               6.99142084e-11 * T3 - 8.84235150e-16 * T4) *
              0.0227226249176305;
  }

  cvbs = result * 8.31446261815324e+07;
}

// Returns the mean enthalpy of the mixture in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBML(const amrex::Real T, const amrex::Real x[], amrex::Real& hbml)
{
  amrex::Real result = 0;
  amrex::Real hml[2];                        //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesEnthalpy(hml, T);

  // perform dot product
  for (int id = 0; id < 2; ++id) {
    result += x[id] * hml[id];
  }

  hbml = result * RT;
}

// Returns mean enthalpy of mixture in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBMS(const amrex::Real T, const amrex::Real y[], amrex::Real& hbms)
{
  amrex::Real result = 0.0;
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;
  const amrex::Real logT = log(T);

  // species with inflection points at T = (1000.0, 6000.0) kelvin
  if (T < 1000) {
    // species 0: CO
    result +=
      y[0] *
      (-1.48904533e+04 * invT2 - 2.92228594e+02 * logT * invT + 5.72452717e+00 -
       4.08811751e-03 * T + 4.85634490e-06 * T2 - 2.71936576e-09 * T3 +
       6.05588365e-13 * T4 - 1.30313188e+04 * invT) *
      0.0357015351660121;
    // species 1: CO2
    result +=
      y[1] *
      (-4.94365054e+04 * invT2 - 6.26411601e+02 * logT * invT + 5.30172524e+00 +
       1.25190691e-03 * T - 7.09102909e-08 * T2 - 1.92249720e-10 * T3 +
       5.69935560e-14 * T4 - 4.52819846e+04 * invT) *
      0.0227226249176305;
  } else if ((1000 <= T) && (T < 6000)) {
    // species 0: CO
    result +=
      y[0] *
      (-4.61919725e+05 * invT2 - 1.94470486e+03 * logT * invT + 5.91671418e+00 -
       2.83214141e-04 * T + 4.66271513e-08 * T2 - 4.46920090e-12 * T3 +
       1.92418711e-16 * T4 - 2.46626108e+03 * invT) *
      0.0357015351660121;
    // species 1: CO2
    result +=
      y[1] *
      (-1.17696242e+05 * invT2 - 1.78879148e+03 * logT * invT + 8.29152319e+00 -
       4.61157839e-05 * T + 1.62122563e-09 * T2 - 4.72763328e-13 * T3 +
       1.26600732e-16 * T4 - 3.90835059e+04 * invT) *
      0.0227226249176305;
  } else {
    // species 0: CO
    result +=
      y[0] *
      (-8.86866296e+08 * invT2 - 7.50037784e+05 * logT * invT + 2.49547498e+02 -
       1.97817555e-02 * T + 1.09925736e-06 * T2 - 3.29602483e-11 * T3 +
       3.99787590e-16 * T4 + 5.70142113e+06 * invT) *
      0.0357015351660121;
    // species 1: CO2
    result +=
      y[1] *
      (+1.54442329e+09 * invT2 + 1.01684706e+06 * logT * invT - 2.56140523e+02 +
       1.68470054e-02 * T - 7.27061446e-07 * T2 + 1.74785521e-11 * T3 -
       1.76847030e-16 * T4 - 8.04321451e+06 * invT) *
      0.0227226249176305;
  }

  const amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  hbms = result * RT;
}

// get mean internal energy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBML(const amrex::Real T, const amrex::Real x[], amrex::Real& ubml)
{
  amrex::Real result = 0;
  amrex::Real uml[2];                        //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesInternalEnergy(uml, T);

  // perform dot product
  for (int id = 0; id < 2; ++id) {
    result += x[id] * uml[id];
  }

  ubml = result * RT;
}

// get mean internal energy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBMS(const amrex::Real T, const amrex::Real y[], amrex::Real& ubms)
{
  amrex::Real result = 0.0;
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;
  const amrex::Real logT = log(T);

  // species with inflection points at T = (1000.0, 6000.0) kelvin
  if (T < 1000) {
    // species 0: CO
    result +=
      y[0] *
      (-1.48904533e+04 * invT2 - 2.92228594e+02 * logT * invT + 4.72452717e+00 -
       4.08811751e-03 * T + 4.85634490e-06 * T2 - 2.71936576e-09 * T3 +
       6.05588365e-13 * T4 - 1.30313188e+04 * invT) *
      0.0357015351660121;
    // species 1: CO2
    result +=
      y[1] *
      (-4.94365054e+04 * invT2 - 6.26411601e+02 * logT * invT + 4.30172524e+00 +
       1.25190691e-03 * T - 7.09102909e-08 * T2 - 1.92249720e-10 * T3 +
       5.69935560e-14 * T4 - 4.52819846e+04 * invT) *
      0.0227226249176305;
  } else if ((1000 <= T) && (T < 6000)) {
    // species 0: CO
    result +=
      y[0] *
      (-4.61919725e+05 * invT2 - 1.94470486e+03 * logT * invT + 4.91671418e+00 -
       2.83214141e-04 * T + 4.66271513e-08 * T2 - 4.46920090e-12 * T3 +
       1.92418711e-16 * T4 - 2.46626108e+03 * invT) *
      0.0357015351660121;
    // species 1: CO2
    result +=
      y[1] *
      (-1.17696242e+05 * invT2 - 1.78879148e+03 * logT * invT + 7.29152319e+00 -
       4.61157839e-05 * T + 1.62122563e-09 * T2 - 4.72763328e-13 * T3 +
       1.26600732e-16 * T4 - 3.90835059e+04 * invT) *
      0.0227226249176305;
  } else {
    // species 0: CO
    result +=
      y[0] *
      (-8.86866296e+08 * invT2 - 7.50037784e+05 * logT * invT + 2.48547498e+02 -
       1.97817555e-02 * T + 1.09925736e-06 * T2 - 3.29602483e-11 * T3 +
       3.99787590e-16 * T4 + 5.70142113e+06 * invT) *
      0.0357015351660121;
    // species 1: CO2
    result +=
      y[1] *
      (+1.54442329e+09 * invT2 + 1.01684706e+06 * logT * invT - 2.57140523e+02 +
       1.68470054e-02 * T - 7.27061446e-07 * T2 + 1.74785521e-11 * T3 -
       1.76847030e-16 * T4 - 8.04321451e+06 * invT) *
      0.0227226249176305;
  }

  const amrex::Real RT = 8.31446261815324e+07 * T; // R*T

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
  amrex::Real sor[2]; //  temporary storage
  speciesEntropy(sor, T);

  // Compute Eq 42
  for (int id = 0; id < 2; ++id) {
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
  amrex::Real sor[2];  //  temporary storage
  amrex::Real x[2];    //  need a ytx conversion
  amrex::Real YOW = 0; // See Eq 4, 6 in CK Manual

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 2; i++) {
    YOW += y[i] * imw(i);
  }
  // Now compute y to x conversion
  x[0] = y[0] / (28.010000 * YOW);
  x[1] = y[1] / (44.009000 * YOW);
  speciesEntropy(sor, T);
  // Perform computation in Eq 42 and 43
  for (int i = 0; i < 2; i++) {
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
  XW += x[0] * 28.010000;                  // CO
  XW += x[1] * 44.009000;                  // CO2
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

  for (int i = 0; i < 2; i++) {
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
  W += c[0] * 28.010000; // CO
  W += c[1] * 44.009000; // CO2

  for (int id = 0; id < 2; ++id) {
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
  XW += x[0] * 28.010000;                    // CO
  XW += x[1] * 44.009000;                    // CO2
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

  for (int i = 0; i < 2; i++) {
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
  W += c[0] * 28.010000; // CO
  W += c[1] * 44.009000; // CO2

  for (int id = 0; id < 2; ++id) {
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

  for (int i = 0; i < 2; i++) {
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
  XW += x[0] * 28.010000; // CO
  XW += x[1] * 44.009000; // CO2
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
  W += c[0] * 28.010000; // CO
  W += c[1] * 44.009000; // CO2

  for (int id = 0; id < 2; ++id) {
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
  cp_R(cpor, T);
}

// get H/RT as a function of T
// for all species (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHORT(const amrex::Real T, amrex::Real hort[])
{
  speciesEnthalpy(hort, T);
}

// get S/R as a function of T
// for all species (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSOR(const amrex::Real T, amrex::Real sor[])
{
  speciesEntropy(sor, T);
}

// convert y[species] (mass fracs) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTX(const amrex::Real y[], amrex::Real x[])
{
  amrex::Real YOW = 0;

  for (int i = 0; i < 2; i++) {
    YOW += y[i] * imw(i);
  }

  amrex::Real YOWINV = 1.0 / YOW;

  for (int i = 0; i < 2; i++) {
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
  for (int i = 0; i < 2; i++) {
    c[i] = y[i] * imw(i);
  }
  for (int i = 0; i < 2; i++) {
    YOW += c[i];
  }

  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // Now compute conversion

  for (int i = 0; i < 2; i++) {
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

  for (int i = 0; i < 2; i++) {
    c[i] = rho * y[i] * imw(i);
  }
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTY(const amrex::Real x[], amrex::Real y[])
{
  amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
  // Compute mean molecular wt first
  XW += x[0] * 28.010000; // CO
  XW += x[1] * 44.009000; // CO2
  // Now compute conversion
  amrex::Real XWinv = 1.0 / XW;
  y[0] = x[0] * 28.010000 * XWinv;
  y[1] = x[1] * 44.009000 * XWinv;
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
  for (int id = 0; id < 2; ++id) {
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
  XW += x[0] * 28.010000; // CO
  XW += x[1] * 44.009000; // CO2
  ROW = rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 2; ++id) {
    c[id] = x[id] * ROW;
  }
}

// convert c[species] (molar conc) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTX(const amrex::Real c[], amrex::Real x[])
{
  amrex::Real sumC = 0;

  // compute sum of c
  for (int id = 0; id < 2; ++id) {
    sumC += c[id];
  }

  //  See Eq 13
  amrex::Real sumCinv = 1.0 / sumC;
  for (int id = 0; id < 2; ++id) {
    x[id] = c[id] * sumCinv;
  }
}

// convert c[species] (molar conc) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTY(const amrex::Real c[], amrex::Real y[])
{
  amrex::Real CW = 0; // See Eq 12 in CK Manual
  // compute denominator in eq 12 first
  CW += c[0] * 28.010000; // CO
  CW += c[1] * 44.009000; // CO2
  // Now compute conversion
  amrex::Real CWinv = 1.0 / CW;
  y[0] = c[0] * 28.010000 * CWinv;
  y[1] = c[1] * 44.009000 * CWinv;
}

// get specific heat at constant volume as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVML(const amrex::Real T, amrex::Real cvml[])
{
  cv_R(cvml, T);

  // convert to chemkin units
  for (int id = 0; id < 2; ++id) {
    cvml[id] *= 8.31446261815324e+07;
  }
}

// get specific heat at constant pressure as a
// function of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPML(const amrex::Real T, amrex::Real cpml[])
{
  cp_R(cpml, T);

  // convert to chemkin units
  for (int id = 0; id < 2; ++id) {
    cpml[id] *= 8.31446261815324e+07;
  }
}

// get internal energy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUML(const amrex::Real T, amrex::Real uml[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesInternalEnergy(uml, T);

  // convert to chemkin units
  for (int id = 0; id < 2; ++id) {
    uml[id] *= RT;
  }
}

// get enthalpy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHML(const amrex::Real T, amrex::Real hml[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesEnthalpy(hml, T);

  // convert to chemkin units
  for (int id = 0; id < 2; ++id) {
    hml[id] *= RT;
  }
}

// Returns the standard-state entropies in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSML(const amrex::Real T, amrex::Real sml[])
{
  speciesEntropy(sml, T);

  // convert to chemkin units
  for (int id = 0; id < 2; ++id) {
    sml[id] *= 8.31446261815324e+07;
  }
}

// Returns the specific heats at constant volume
// in mass units (Eq. 29)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVMS(const amrex::Real T, amrex::Real cvms[])
{
  cv_R(cvms, T);
  // multiply by R/molecularweight
  cvms[0] *= 2.968390795484913e+06; // CO
  cvms[1] *= 1.889264154639560e+06; // CO2
}

// Returns the specific heats at constant pressure
// in mass units (Eq. 26)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPMS(const amrex::Real T, amrex::Real cpms[])
{
  cp_R(cpms, T);
  // multiply by R/molecularweight
  cpms[0] *= 2.968390795484913e+06; // CO
  cpms[1] *= 1.889264154639560e+06; // CO2
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUMS(const amrex::Real T, amrex::Real ums[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  speciesInternalEnergy(ums, T);

  for (int i = 0; i < 2; i++) {
    ums[i] *= RT * imw(i);
  }
}

// Returns enthalpy in mass units (Eq 27.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHMS(const amrex::Real T, amrex::Real hms[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  speciesEnthalpy(hms, T);

  for (int i = 0; i < 2; i++) {
    hms[i] *= RT * imw(i);
  }
}

// Returns the entropies in mass units (Eq 28.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSMS(const amrex::Real T, amrex::Real sms[])
{
  speciesEntropy(sms, T);
  // multiply by R/molecularweight
  sms[0] *= 2.968390795484913e+06; // CO
  sms[1] *= 1.889264154639560e+06; // CO2
}

// GPU version of productionRate: no more use of thermo namespace vectors
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_qfqr(
  amrex::Real* /*qf*/,
  amrex::Real* /*qr*/,
  const amrex::Real* /*sc*/,
  const amrex::Real* /*sc_qss*/,
  const amrex::Real /*T*/,
  const amrex::Real /*invT*/,
  const amrex::Real /*logT*/)
{
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
productionRate(
  amrex::Real* wdot, const amrex::Real* /*sc*/, const amrex::Real /*T*/)
{

  for (int i = 0; i < 2; ++i) {
    wdot[i] = 0.0;
  }
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[])
{

  // convert to SI
  for (int id = 0; id < 2; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, T);

  // convert to chemkin units
  for (int id = 0; id < 2; ++id) {
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
  amrex::Real c[2]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 2; i++) {
    YOW += y[i] * imw(i);
  }
  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // multiply by 1e6 so c goes to SI
  PWORT *= 1e6;
  // Now compute conversion (and go to SI)
  for (int i = 0; i < 2; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 2; ++id) {
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
  amrex::Real c[2]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 2; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 2; ++id) {
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
  amrex::Real c[2]; // temporary storage

  // See Eq 8 with an extra 1e6 so c goes to SI
  for (int i = 0; i < 2; i++) {
    c[i] = 1e6 * rho * y[i] * imw(i);
  }

  // call productionRate
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 2; ++id) {
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
  amrex::Real c[2];   // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 28.010000; // CO
  XW += x[1] * 44.009000; // CO2
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 2; ++id) {
    c[id] = x[id] * ROW;
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 2; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int kcharge[])
{
  kcharge[0] = 0; // CO
  kcharge[1] = 0; // CO2
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real zk[])
{

  int kchrg[2];
  CKCHRG(kchrg);

  for (int id = 0; id < 2; ++id) {
    zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
  }
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
  for (int i = 0; i < 9; i++) {
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
  amrex::Real c[2];

  for (int k = 0; k < 2; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian_precond(J, c, *Tp, *HP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 2; k++) {
    J[6 + k] *= 1.e-6;
    J[k * 3 + 2] *= 1.e6;
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
  for (int i = 0; i < 9; i++) {
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
  amrex::Real c[2];

  for (int k = 0; k < 2; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian(J, c, *Tp, *consP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 2; k++) {
    J[6 + k] *= 1.e-6;
    J[k * 3 + 2] *= 1.e6;
  }
}

// Transport function declarations

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENIMC(int* LENIMC)
{
  *LENIMC = 8;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENRMC(int* LENRMC)
{
  *LENRMC = 130;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNO(int* NO)
{
  *NO = 4;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKK(int* KK)
{
  *KK = 2;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLITE(int* NLITE)
{
  *NLITE = 0;
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
  WT[0] = 2.80100000E+01;
  WT[1] = 4.40090000E+01;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetEPS(amrex::Real* EPS)
{
  EPS[0] = 9.81000000E+01;
  EPS[1] = 2.44000000E+02;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetSIG(amrex::Real* SIG)
{
  SIG[0] = 3.65000000E+00;
  SIG[1] = 3.76300000E+00;
}

// the dipole moment in Debye
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetDIP(amrex::Real* DIP)
{
  DIP[0] = 0.00000000E+00;
  DIP[1] = 0.00000000E+00;
}

// the polarizability in cubic Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPOL(amrex::Real* POL)
{
  POL[0] = 1.95000000E+00;
  POL[1] = 2.65000000E+00;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetZROT(amrex::Real* ZROT)
{
  ZROT[0] = 1.80000000E+00;
  ZROT[1] = 2.10000000E+00;
}

// 0: monoatomic, 1: linear, 2: nonlinear
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLIN(int* NLIN)
{
  NLIN[0] = 1;
  NLIN[1] = 1;
}

// Poly fits for the viscosities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFETA(amrex::Real* COFETA)
{
  COFETA[0] = -1.51727145E+01;
  COFETA[1] = 1.71554871E+00;
  COFETA[2] = -1.30197682E-01;
  COFETA[3] = 5.24650173E-03;
  COFETA[4] = -1.95762007E+01;
  COFETA[5] = 3.21383321E+00;
  COFETA[6] = -2.97666118E-01;
  COFETA[7] = 1.14407820E-02;
}

// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFLAM(amrex::Real* COFLAM)
{
  COFLAM[0] = -6.59930001E+00;
  COFLAM[1] = 5.15042377E+00;
  COFLAM[2] = -6.24566626E-01;
  COFLAM[3] = 2.92450231E-02;
  COFLAM[4] = -7.43758586E+00;
  COFLAM[5] = 4.56044161E+00;
  COFLAM[6] = -4.38897355E-01;
  COFLAM[7] = 1.71585898E-02;
}

// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFD(amrex::Real* COFD)
{
  COFD[0] = -1.36864543E+01;
  COFD[1] = 2.64212959E+00;
  COFD[2] = -1.18476768E-01;
  COFD[3] = 4.71039964E-03;
  COFD[4] = -1.55661602E+01;
  COFD[5] = 3.20553742E+00;
  COFD[6] = -1.80905498E-01;
  COFD[7] = 6.99284071E-03;
  COFD[8] = -1.55661602E+01;
  COFD[9] = 3.20553742E+00;
  COFD[10] = -1.80905498E-01;
  COFD[11] = 6.99284071E-03;
  COFD[12] = -1.80696015E+01;
  COFD[13] = 3.97303883E+00;
  COFD[14] = -2.66996069E-01;
  COFD[15] = 1.02065369E-02;
}

// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKTDIF(int* /*KTDIF*/)
{
}

// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFTD(amrex::Real* /*COFTD*/)
{
}

// compute the critical parameters for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS(
  amrex::Real* Tci, amrex::Real* ai, amrex::Real* bi, amrex::Real* acentric_i)
{

  amrex::Real EPS[2];
  amrex::Real SIG[2];
  amrex::Real wt[2];
  amrex::Real Rcst = 83.144598; // in bar [CGS] !

  egtransetEPS(EPS);
  egtransetSIG(SIG);
  get_mw(wt);

  // species 0: CO
  // Imported from NIST
  Tci[0] = 132.850000;
  ai[0] = 1e6 * 0.42748 * Rcst * Rcst * Tci[0] * Tci[0] /
          (28.010000 * 28.010000 * 34.940000);
  bi[0] = 0.08664 * Rcst * Tci[0] / (28.010000 * 34.940000);
  acentric_i[0] = 0.045000;

  // species 1: CO2
  // Imported from NIST
  Tci[1] = 304.120000;
  ai[1] = 1e6 * 0.42748 * Rcst * Rcst * Tci[1] * Tci[1] /
          (44.009950 * 44.009950 * 73.740000);
  bi[1] = 0.08664 * Rcst * Tci[1] / (44.009950 * 73.740000);
  acentric_i[1] = 0.225000;
}

// compute the critical parameter quantities for each species for SRK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS_SRK(
  amrex::Real* sqrtOneOverTc,
  amrex::Real* sqrtAsti,
  amrex::Real* Bi,
  amrex::Real* Fomega)
{

  // species 0: CO
  // Imported from NIST
  sqrtOneOverTc[0] = 8.6759935530451e-02;
  sqrtAsti[0] = 4.3619324712923e+04;
  Bi[0] = 9.7786303823350e-01;
  Fomega[0] = 5.5459948367500e-01;

  // species 1: CO2
  // Imported from NIST
  sqrtOneOverTc[1] = 5.7342616962522e-02;
  sqrtAsti[1] = 4.3745610142919e+04;
  Bi[1] = 6.7506167841745e-01;
  Fomega[1] = 8.2653709187500e-01;
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
