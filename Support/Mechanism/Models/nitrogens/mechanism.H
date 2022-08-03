#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  N
*/

// Species
#define N2_ID 0
#define N2a_ID 1

#define NUM_ELEMENTS 1
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
void CKINU(int* i, int* nspec, int* ki, int* nu);
void CKKFKR(
  amrex::Real* P,
  amrex::Real* T,
  amrex::Real* x,
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
CKINDX(int* mm, int* kk, int* ii, int* nfit)
{
  *mm = 1;
  *kk = 2;
  *ii = 0;
  *nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_imw(amrex::Real* imw_new)
{
  imw_new[0] = 1.0 / 28.014000; // N2
  imw_new[1] = 1.0 / 28.014000; // N2a
}

//  molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_mw(amrex::Real* mw_new)
{
  mw_new[0] = 28.014000; // N2
  mw_new[1] = 28.014000; // N2a
}

//  Returns R, Rc, Patm
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRP(amrex::Real* ru, amrex::Real* ruc, amrex::Real* pa)
{
  *ru = 8.31446261815324e+07;
  *ruc = 1.98721558317399615845;
  *pa = 1.01325e+06;
}

// compute Cv/R at the given temperature
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cv_R(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: N2
    species[0] = +2.29867700e+00 + 1.40824000e-03 * tc[1] -
                 3.96322200e-06 * tc[2] + 5.64151500e-09 * tc[3] -
                 2.44485500e-12 * tc[4];
    // species 1: N2a
    species[1] = +2.29867700e+00 + 1.40824000e-03 * tc[1] -
                 3.96322200e-06 * tc[2] + 5.64151500e-09 * tc[3] -
                 2.44485500e-12 * tc[4];
  } else {
    // species 0: N2
    species[0] = +1.92664000e+00 + 1.48797700e-03 * tc[1] -
                 5.68476100e-07 * tc[2] + 1.00970400e-10 * tc[3] -
                 6.75335100e-15 * tc[4];
    // species 1: N2a
    species[1] = +1.92664000e+00 + 1.48797700e-03 * tc[1] -
                 5.68476100e-07 * tc[2] + 1.00970400e-10 * tc[3] -
                 6.75335100e-15 * tc[4];
  }
  return;
}

// compute Cp/R at the given temperature
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cp_R(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: N2
    species[0] = +3.29867700e+00 + 1.40824000e-03 * tc[1] -
                 3.96322200e-06 * tc[2] + 5.64151500e-09 * tc[3] -
                 2.44485500e-12 * tc[4];
    // species 1: N2a
    species[1] = +3.29867700e+00 + 1.40824000e-03 * tc[1] -
                 3.96322200e-06 * tc[2] + 5.64151500e-09 * tc[3] -
                 2.44485500e-12 * tc[4];
  } else {
    // species 0: N2
    species[0] = +2.92664000e+00 + 1.48797700e-03 * tc[1] -
                 5.68476100e-07 * tc[2] + 1.00970400e-10 * tc[3] -
                 6.75335100e-15 * tc[4];
    // species 1: N2a
    species[1] = +2.92664000e+00 + 1.48797700e-03 * tc[1] -
                 5.68476100e-07 * tc[2] + 1.00970400e-10 * tc[3] -
                 6.75335100e-15 * tc[4];
  }
  return;
}

// compute the g/(RT) at the given temperature
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
gibbs(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: N2
    species[0] = -1.020900000000000e+03 * invT - 6.516950000000001e-01 -
                 3.298677000000000e+00 * tc[0] - 7.041200000000000e-04 * tc[1] +
                 6.605369999999999e-07 * tc[2] - 4.701262500000001e-10 * tc[3] +
                 1.222427500000000e-13 * tc[4];
    // species 1: N2a
    species[1] = -1.020900000000000e+03 * invT - 6.516950000000001e-01 -
                 3.298677000000000e+00 * tc[0] - 7.041200000000000e-04 * tc[1] +
                 6.605369999999999e-07 * tc[2] - 4.701262500000001e-10 * tc[3] +
                 1.222427500000000e-13 * tc[4];
  } else {
    // species 0: N2
    species[0] = -9.227977000000000e+02 * invT - 3.053888000000000e+00 -
                 2.926640000000000e+00 * tc[0] - 7.439885000000000e-04 * tc[1] +
                 9.474601666666666e-08 * tc[2] - 8.414199999999999e-12 * tc[3] +
                 3.376675500000000e-16 * tc[4];
    // species 1: N2a
    species[1] = -9.227977000000000e+02 * invT - 3.053888000000000e+00 -
                 2.926640000000000e+00 * tc[0] - 7.439885000000000e-04 * tc[1] +
                 9.474601666666666e-08 * tc[2] - 8.414199999999999e-12 * tc[3] +
                 3.376675500000000e-16 * tc[4];
  }
  return;
}

// compute the a/(RT) at the given temperature
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
helmholtz(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: N2
    species[0] = -1.02090000e+03 * invT - 1.65169500e+00 -
                 3.29867700e+00 * tc[0] - 7.04120000e-04 * tc[1] +
                 6.60537000e-07 * tc[2] - 4.70126250e-10 * tc[3] +
                 1.22242750e-13 * tc[4];
    // species 1: N2a
    species[1] = -1.02090000e+03 * invT - 1.65169500e+00 -
                 3.29867700e+00 * tc[0] - 7.04120000e-04 * tc[1] +
                 6.60537000e-07 * tc[2] - 4.70126250e-10 * tc[3] +
                 1.22242750e-13 * tc[4];
  } else {
    // species 0: N2
    species[0] = -9.22797700e+02 * invT - 4.05388800e+00 -
                 2.92664000e+00 * tc[0] - 7.43988500e-04 * tc[1] +
                 9.47460167e-08 * tc[2] - 8.41420000e-12 * tc[3] +
                 3.37667550e-16 * tc[4];
    // species 1: N2a
    species[1] = -9.22797700e+02 * invT - 4.05388800e+00 -
                 2.92664000e+00 * tc[0] - 7.43988500e-04 * tc[1] +
                 9.47460167e-08 * tc[2] - 8.41420000e-12 * tc[3] +
                 3.37667550e-16 * tc[4];
  }
  return;
}

// compute the e/(RT) at the given temperature
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesInternalEnergy(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: N2
    species[0] = +2.29867700e+00 + 7.04120000e-04 * tc[1] -
                 1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
                 4.88971000e-13 * tc[4] - 1.02090000e+03 * invT;
    // species 1: N2a
    species[1] = +2.29867700e+00 + 7.04120000e-04 * tc[1] -
                 1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
                 4.88971000e-13 * tc[4] - 1.02090000e+03 * invT;
  } else {
    // species 0: N2
    species[0] = +1.92664000e+00 + 7.43988500e-04 * tc[1] -
                 1.89492033e-07 * tc[2] + 2.52426000e-11 * tc[3] -
                 1.35067020e-15 * tc[4] - 9.22797700e+02 * invT;
    // species 1: N2a
    species[1] = +1.92664000e+00 + 7.43988500e-04 * tc[1] -
                 1.89492033e-07 * tc[2] + 2.52426000e-11 * tc[3] -
                 1.35067020e-15 * tc[4] - 9.22797700e+02 * invT;
  }
  return;
}

// compute the h/(RT) at the given temperature (Eq 20)
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEnthalpy(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: N2
    species[0] = +3.29867700e+00 + 7.04120000e-04 * tc[1] -
                 1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
                 4.88971000e-13 * tc[4] - 1.02090000e+03 * invT;
    // species 1: N2a
    species[1] = +3.29867700e+00 + 7.04120000e-04 * tc[1] -
                 1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
                 4.88971000e-13 * tc[4] - 1.02090000e+03 * invT;
  } else {
    // species 0: N2
    species[0] = +2.92664000e+00 + 7.43988500e-04 * tc[1] -
                 1.89492033e-07 * tc[2] + 2.52426000e-11 * tc[3] -
                 1.35067020e-15 * tc[4] - 9.22797700e+02 * invT;
    // species 1: N2a
    species[1] = +2.92664000e+00 + 7.43988500e-04 * tc[1] -
                 1.89492033e-07 * tc[2] + 2.52426000e-11 * tc[3] -
                 1.35067020e-15 * tc[4] - 9.22797700e+02 * invT;
  }
  return;
}

// compute the S/R at the given temperature (Eq 21)
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEntropy(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: N2
    species[0] = +3.29867700e+00 * tc[0] + 1.40824000e-03 * tc[1] -
                 1.98161100e-06 * tc[2] + 1.88050500e-09 * tc[3] -
                 6.11213750e-13 * tc[4] + 3.95037200e+00;
    // species 1: N2a
    species[1] = +3.29867700e+00 * tc[0] + 1.40824000e-03 * tc[1] -
                 1.98161100e-06 * tc[2] + 1.88050500e-09 * tc[3] -
                 6.11213750e-13 * tc[4] + 3.95037200e+00;
  } else {
    // species 0: N2
    species[0] = +2.92664000e+00 * tc[0] + 1.48797700e-03 * tc[1] -
                 2.84238050e-07 * tc[2] + 3.36568000e-11 * tc[3] -
                 1.68833775e-15 * tc[4] + 5.98052800e+00;
    // species 1: N2a
    species[1] = +2.92664000e+00 * tc[0] + 1.48797700e-03 * tc[1] -
                 2.84238050e-07 * tc[2] + 3.36568000e-11 * tc[3] -
                 1.68833775e-15 * tc[4] + 5.98052800e+00;
  }
  return;
}

// Returns the mean specific heat at CP (Eq. 33)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBL(amrex::Real* T, amrex::Real* x, amrex::Real* cpbl)
{
  int id; // loop counter
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real cpor[2];                                //  temporary storage
  cp_R(cpor, tc);

  // perform dot product
  for (id = 0; id < 2; ++id) {
    result += x[id] * cpor[id];
  }

  *cpbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CP (Eq. 34)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBS(amrex::Real* T, amrex::Real* y, amrex::Real* cpbs)
{
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real cpor[2], tresult[2];                    //  temporary storage
  amrex::Real imw[2];

  get_imw(imw);

  cp_R(cpor, tc);

  for (int i = 0; i < 2; i++) {
    tresult[i] = cpor[i] * y[i] * imw[i];
  }
  for (int i = 0; i < 2; i++) {
    result += tresult[i];
  }

  *cpbs = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 35)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBL(amrex::Real* T, amrex::Real* x, amrex::Real* cvbl)
{
  int id; // loop counter
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real cvor[2];                                //  temporary storage
  cv_R(cvor, tc);

  // perform dot product
  for (id = 0; id < 2; ++id) {
    result += x[id] * cvor[id];
  }

  *cvbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 36)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBS(amrex::Real* T, amrex::Real* y, amrex::Real* cvbs)
{
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real cvor[2];                                //  temporary storage
  amrex::Real imw[2];

  get_imw(imw);

  cv_R(cvor, tc);

  // multiply by y/molecularweight
  result += cvor[0] * y[0] * imw[0]; // N2
  result += cvor[1] * y[1] * imw[1]; // N2a

  *cvbs = result * 8.31446261815324e+07;
}

// Returns the mean enthalpy of the mixture in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBML(amrex::Real* T, amrex::Real* x, amrex::Real* hbml)
{
  int id; // loop counter
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real hml[2];                                 //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesEnthalpy(hml, tc);

  // perform dot product
  for (id = 0; id < 2; ++id) {
    result += x[id] * hml[id];
  }

  *hbml = result * RT;
}

// Returns mean enthalpy of mixture in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBMS(amrex::Real* T, amrex::Real* y, amrex::Real* hbms)
{
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real hml[2], tmp[2];                         //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[2];

  get_imw(imw);

  speciesEnthalpy(hml, tc);

  int id;
  for (id = 0; id < 2; ++id) {
    tmp[id] = y[id] * hml[id] * imw[id];
  }
  for (id = 0; id < 2; ++id) {
    result += tmp[id];
  }

  *hbms = result * RT;
}

// get mean internal energy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBML(amrex::Real* T, amrex::Real* x, amrex::Real* ubml)
{
  int id; // loop counter
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real uml[2];                                 //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesInternalEnergy(uml, tc);

  // perform dot product
  for (id = 0; id < 2; ++id) {
    result += x[id] * uml[id];
  }

  *ubml = result * RT;
}

// get mean internal energy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBMS(amrex::Real* T, amrex::Real* y, amrex::Real* ubms)
{
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real ums[2];                                 //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[2];

  get_imw(imw);

  speciesInternalEnergy(ums, tc);

  // perform dot product + scaling by wt
  result += y[0] * ums[0] * imw[0]; // N2
  result += y[1] * ums[1] * imw[1]; // N2a

  *ubms = result * RT;
}

// get mixture entropy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSBML(amrex::Real* P, amrex::Real* T, amrex::Real* x, amrex::Real* sbml)
{
  int id; // loop counter
  amrex::Real result = 0;
  // Log of normalized pressure in cgs units dynes/cm^2 by Patm
  amrex::Real logPratio = log(*P / 1013250.0);
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real sor[2]; //  temporary storage
  speciesEntropy(sor, tc);

  // Compute Eq 42
  for (id = 0; id < 2; ++id) {
    result += x[id] * (sor[id] - log((x[id] + 1e-100)) - logPratio);
  }

  *sbml = result * 8.31446261815324e+07;
}

// get mixture entropy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSBMS(amrex::Real* P, amrex::Real* T, amrex::Real* y, amrex::Real* sbms)
{
  amrex::Real result = 0;
  // Log of normalized pressure in cgs units dynes/cm^2 by Patm
  amrex::Real logPratio = log(*P / 1013250.0);
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real sor[2];  //  temporary storage
  amrex::Real x[2];    //  need a ytx conversion
  amrex::Real YOW = 0; // See Eq 4, 6 in CK Manual
  amrex::Real imw[2];

  get_imw(imw);

  // Compute inverse of mean molecular wt first
  YOW += y[0] * imw[0]; // N2
  YOW += y[1] * imw[1]; // N2a
  // Now compute y to x conversion
  x[0] = y[0] / (28.014000 * YOW);
  x[1] = y[1] / (28.014000 * YOW);
  speciesEntropy(sor, tc);
  // Perform computation in Eq 42 and 43
  result += x[0] * (sor[0] - log((x[0] + 1e-100)) - logPratio);
  result += x[1] * (sor[1] - log((x[1] + 1e-100)) - logPratio);
  // Scale by R/W
  *sbms = result * 8.31446261815324e+07 * YOW;
}

//  get temperature given internal energy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_EY(amrex::Real* e, amrex::Real* y, amrex::Real* t, int* ierr)
{
#ifdef CONVERGENCE
  const int maxiter = 5000;
  const amrex::Real tol = 1.e-12;
#else
  const int maxiter = 200;
  const amrex::Real tol = 1.e-6;
#endif
  amrex::Real ein = *e;
  amrex::Real tmin = 90;   // max lower bound for thermo def
  amrex::Real tmax = 4000; // min upper bound for thermo def
  amrex::Real e1, emin, emax, cv, t1, dt;
  int i; //  loop counter
  CKUBMS(&tmin, y, &emin);
  CKUBMS(&tmax, y, &emax);
  if (ein < emin) {
    // Linear Extrapolation below tmin
    CKCVBS(&tmin, y, &cv);
    *t = tmin - (emin - ein) / cv;
    *ierr = 1;
    return;
  }
  if (ein > emax) {
    // Linear Extrapolation above tmax
    CKCVBS(&tmax, y, &cv);
    *t = tmax - (emax - ein) / cv;
    *ierr = 1;
    return;
  }
  t1 = *t;
  if (t1 < tmin || t1 > tmax) {
    t1 = tmin + (tmax - tmin) / (emax - emin) * (ein - emin);
  }
  for (i = 0; i < maxiter; ++i) {
    CKUBMS(&t1, y, &e1);
    CKCVBS(&t1, y, &cv);
    dt = (ein - e1) / cv;
    if (dt > 100.) {
      dt = 100.;
    } else if (dt < -100.) {
      dt = -100.;
    } else if (fabs(dt) < tol)
      break;
    else if (t1 + dt == t1)
      break;
    t1 += dt;
  }
  *t = t1;
  *ierr = 0;
  return;
}

//  get temperature given enthalpy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_HY(amrex::Real* h, amrex::Real* y, amrex::Real* t, int* ierr)
{
#ifdef CONVERGENCE
  const int maxiter = 5000;
  const amrex::Real tol = 1.e-12;
#else
  const int maxiter = 200;
  const amrex::Real tol = 1.e-6;
#endif
  amrex::Real hin = *h;
  amrex::Real tmin = 90;   // max lower bound for thermo def
  amrex::Real tmax = 4000; // min upper bound for thermo def
  amrex::Real h1, hmin, hmax, cp, t1, dt;
  int i; //  loop counter
  CKHBMS(&tmin, y, &hmin);
  CKHBMS(&tmax, y, &hmax);
  if (hin < hmin) {
    // Linear Extrapolation below tmin
    CKCPBS(&tmin, y, &cp);
    *t = tmin - (hmin - hin) / cp;
    *ierr = 1;
    return;
  }
  if (hin > hmax) {
    // Linear Extrapolation above tmax
    CKCPBS(&tmax, y, &cp);
    *t = tmax - (hmax - hin) / cp;
    *ierr = 1;
    return;
  }
  t1 = *t;
  if (t1 < tmin || t1 > tmax) {
    t1 = tmin + (tmax - tmin) / (hmax - hmin) * (hin - hmin);
  }
  for (i = 0; i < maxiter; ++i) {
    CKHBMS(&t1, y, &h1);
    CKCPBS(&t1, y, &cp);
    dt = (hin - h1) / cp;
    if (dt > 100.) {
      dt = 100.;
    } else if (dt < -100.) {
      dt = -100.;
    } else if (fabs(dt) < tol)
      break;
    else if (t1 + dt == t1)
      break;
    t1 += dt;
  }
  *t = t1;
  *ierr = 0;
  return;
}

// Compute P = rhoRT/W(x)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPX(amrex::Real* rho, amrex::Real* T, amrex::Real* x, amrex::Real* P)
{
  amrex::Real XW = 0;                           //  To hold mean molecular wt
  XW += x[0] * 28.014000;                       // N2
  XW += x[1] * 28.014000;                       // N2a
  *P = *rho * 8.31446261815324e+07 * (*T) / XW; // P = rho*R*T/W

  return;
}

// Compute P = rhoRT/W(y)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPY(amrex::Real* rho, amrex::Real* T, amrex::Real* y, amrex::Real* P)
{
  amrex::Real YOW = 0; //  for computing mean MW
  amrex::Real imw[2];

  get_imw(imw);

  YOW += y[0] * imw[0];                          // N2
  YOW += y[1] * imw[1];                          // N2a
  *P = *rho * 8.31446261815324e+07 * (*T) * YOW; // P = rho*R*T/W

  return;
}

// Compute P = rhoRT/W(c)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPC(amrex::Real* rho, amrex::Real* T, amrex::Real* c, amrex::Real* P)
{
  int id; // loop counter
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 28.014000; // N2
  W += c[1] * 28.014000; // N2a

  for (id = 0; id < 2; ++id) {
    sumC += c[id];
  }
  *P = *rho * 8.31446261815324e+07 * (*T) * sumC / W; // P = rho*R*T/W

  return;
}

// Compute rho = PW(x)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOX(amrex::Real* P, amrex::Real* T, amrex::Real* x, amrex::Real* rho)
{
  amrex::Real XW = 0;                             //  To hold mean molecular wt
  XW += x[0] * 28.014000;                         // N2
  XW += x[1] * 28.014000;                         // N2a
  *rho = *P * XW / (8.31446261815324e+07 * (*T)); // rho = P*W/(R*T)

  return;
}

// Compute rho = P*W(y)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOY(amrex::Real* P, amrex::Real* T, amrex::Real* y, amrex::Real* rho)
{
  amrex::Real YOW = 0;
  amrex::Real tmp[2];
  amrex::Real imw[2];

  get_imw(imw);

  for (int i = 0; i < 2; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 2; i++) {
    YOW += tmp[i];
  }

  *rho = *P / (8.31446261815324e+07 * (*T) * YOW); // rho = P*W/(R*T)
  return;
}

// Compute rho = P*W(c)/(R*T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOC(amrex::Real* P, amrex::Real* T, amrex::Real* c, amrex::Real* rho)
{
  int id; // loop counter
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 28.014000; // N2
  W += c[1] * 28.014000; // N2a

  for (id = 0; id < 2; ++id) {
    sumC += c[id];
  }
  *rho = *P * W / (sumC * (*T) * 8.31446261815324e+07); // rho = PW/(R*T)

  return;
}

// get molecular weight for all species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWT(amrex::Real* wt)
{
  get_mw(wt);
}

// given y[species]: mass fractions
// s mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWY(amrex::Real* y, amrex::Real* wtm)
{
  amrex::Real YOW = 0;
  amrex::Real tmp[2];
  amrex::Real imw[2];

  get_imw(imw);

  for (int i = 0; i < 2; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 2; i++) {
    YOW += tmp[i];
  }

  *wtm = 1.0 / YOW;
  return;
}

// given x[species]: mole fractions
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWX(amrex::Real* x, amrex::Real* wtm)
{
  amrex::Real XW = 0;     //  see Eq 4 in CK Manual
  XW += x[0] * 28.014000; // N2
  XW += x[1] * 28.014000; // N2a
  *wtm = XW;

  return;
}

// given c[species]: molar concentration
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWC(amrex::Real* c, amrex::Real* wtm)
{
  int id; // loop counter
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 28.014000; // N2
  W += c[1] * 28.014000; // N2a

  for (id = 0; id < 2; ++id) {
    sumC += c[id];
  }
  //  CK provides no guard against divison by zero
  *wtm = W / sumC;

  return;
}

// get Cp/R as a function of T
// for all species (Eq 19)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPOR(amrex::Real* T, amrex::Real* cpor)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cp_R(cpor, tc);
}

// get H/RT as a function of T
// for all species (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHORT(amrex::Real* T, amrex::Real* hort)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  speciesEnthalpy(hort, tc);
}

// get S/R as a function of T
// for all species (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSOR(amrex::Real* T, amrex::Real* sor)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  speciesEntropy(sor, tc);
}

// convert y[species] (mass fracs) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTX(amrex::Real* y, amrex::Real* x)
{
  amrex::Real YOW = 0;
  amrex::Real tmp[2];
  amrex::Real imw[2];

  get_imw(imw);

  for (int i = 0; i < 2; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 2; i++) {
    YOW += tmp[i];
  }

  amrex::Real YOWINV = 1.0 / YOW;

  for (int i = 0; i < 2; i++) {
    x[i] = y[i] * imw[i] * YOWINV;
  }
  return;
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTCP(amrex::Real* P, amrex::Real* T, amrex::Real* y, amrex::Real* c)
{
  amrex::Real YOW = 0;
  amrex::Real PWORT;
  amrex::Real imw[2];

  get_imw(imw);

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 2; i++) {
    c[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 2; i++) {
    YOW += c[i];
  }

  // PW/RT (see Eq. 7)
  PWORT = (*P) / (YOW * 8.31446261815324e+07 * (*T));
  // Now compute conversion

  for (int i = 0; i < 2; i++) {
    c[i] = PWORT * y[i] * imw[i];
  }
  return;
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTCR(amrex::Real* rho, amrex::Real* /*T*/, amrex::Real* y, amrex::Real* c)
{
  amrex::Real imw[2];

  get_imw(imw);

  for (int i = 0; i < 2; i++) {
    c[i] = (*rho) * y[i] * imw[i];
  }
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTY(amrex::Real* x, amrex::Real* y)
{
  amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
  // Compute mean molecular wt first
  XW += x[0] * 28.014000; // N2
  XW += x[1] * 28.014000; // N2a
  // Now compute conversion
  amrex::Real XWinv = 1.0 / XW;
  y[0] = x[0] * 28.014000 * XWinv;
  y[1] = x[1] * 28.014000 * XWinv;

  return;
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTCP(amrex::Real* P, amrex::Real* T, amrex::Real* x, amrex::Real* c)
{
  int id;                                                  // loop counter
  amrex::Real PORT = (*P) / (8.31446261815324e+07 * (*T)); // P/RT

  // Compute conversion, see Eq 10
  for (id = 0; id < 2; ++id) {
    c[id] = x[id] * PORT;
  }

  return;
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTCR(amrex::Real* rho, amrex::Real* /*T*/, amrex::Real* x, amrex::Real* c)
{
  int id;             // loop counter
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 28.014000; // N2
  XW += x[1] * 28.014000; // N2a
  ROW = (*rho) / XW;

  // Compute conversion, see Eq 11
  for (id = 0; id < 2; ++id) {
    c[id] = x[id] * ROW;
  }

  return;
}

// convert c[species] (molar conc) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTX(amrex::Real* c, amrex::Real* x)
{
  int id; // loop counter
  amrex::Real sumC = 0;

  // compute sum of c
  for (id = 0; id < 2; ++id) {
    sumC += c[id];
  }

  //  See Eq 13
  amrex::Real sumCinv = 1.0 / sumC;
  for (id = 0; id < 2; ++id) {
    x[id] = c[id] * sumCinv;
  }

  return;
}

// convert c[species] (molar conc) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTY(amrex::Real* c, amrex::Real* y)
{
  amrex::Real CW = 0; // See Eq 12 in CK Manual
  // compute denominator in eq 12 first
  CW += c[0] * 28.014000; // N2
  CW += c[1] * 28.014000; // N2a
  // Now compute conversion
  amrex::Real CWinv = 1.0 / CW;
  y[0] = c[0] * 28.014000 * CWinv;
  y[1] = c[1] * 28.014000 * CWinv;

  return;
}

// get specific heat at constant volume as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVML(amrex::Real* T, amrex::Real* cvml)
{
  int id;              // loop counter
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cv_R(cvml, tc);

  // convert to chemkin units
  for (id = 0; id < 2; ++id) {
    cvml[id] *= 8.31446261815324e+07;
  }
}

// get specific heat at constant pressure as a
// function of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPML(amrex::Real* T, amrex::Real* cpml)
{
  int id;              // loop counter
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cp_R(cpml, tc);

  // convert to chemkin units
  for (id = 0; id < 2; ++id) {
    cpml[id] *= 8.31446261815324e+07;
  }
}

// get internal energy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUML(amrex::Real* T, amrex::Real* uml)
{
  int id;              // loop counter
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesInternalEnergy(uml, tc);

  // convert to chemkin units
  for (id = 0; id < 2; ++id) {
    uml[id] *= RT;
  }
}

// get enthalpy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHML(amrex::Real* T, amrex::Real* hml)
{
  int id;              // loop counter
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesEnthalpy(hml, tc);

  // convert to chemkin units
  for (id = 0; id < 2; ++id) {
    hml[id] *= RT;
  }
}

// Returns the standard-state entropies in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSML(amrex::Real* T, amrex::Real* sml)
{
  int id;              // loop counter
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  speciesEntropy(sml, tc);

  // convert to chemkin units
  for (id = 0; id < 2; ++id) {
    sml[id] *= 8.31446261815324e+07;
  }
}

// Returns the specific heats at constant volume
// in mass units (Eq. 29)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVMS(amrex::Real* T, amrex::Real* cvms)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cv_R(cvms, tc);
  // multiply by R/molecularweight
  cvms[0] *= 2.967966951578939e+06; // N2
  cvms[1] *= 2.967966951578939e+06; // N2a
}

// Returns the specific heats at constant pressure
// in mass units (Eq. 26)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPMS(amrex::Real* T, amrex::Real* cpms)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cp_R(cpms, tc);
  // multiply by R/molecularweight
  cpms[0] *= 2.967966951578939e+06; // N2
  cpms[1] *= 2.967966951578939e+06; // N2a
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUMS(amrex::Real* T, amrex::Real* ums)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[2];

  get_imw(imw);

  speciesInternalEnergy(ums, tc);

  for (int i = 0; i < 2; i++) {
    ums[i] *= RT * imw[i];
  }
}

// Returns enthalpy in mass units (Eq 27.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHMS(amrex::Real* T, amrex::Real* hms)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[2];

  get_imw(imw);

  speciesEnthalpy(hms, tc);

  for (int i = 0; i < 2; i++) {
    hms[i] *= RT * imw[i];
  }
}

// Returns the entropies in mass units (Eq 28.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSMS(amrex::Real* T, amrex::Real* sms)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  speciesEntropy(sms, tc);
  // multiply by R/molecularweight
  sms[0] *= 2.967966951578939e+06; // N2
  sms[1] *= 2.967966951578939e+06; // N2a
}

// GPU version of productionRate: no more use of thermo namespace vectors
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_qfqr(
  amrex::Real* qf,
  amrex::Real* qr,
  amrex::Real* sc,
  amrex::Real* sc_qss,
  const amrex::Real* tc,
  amrex::Real invT)
{

  return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
productionRate(amrex::Real* wdot, amrex::Real* sc, amrex::Real T)
{
  const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
  const amrex::Real invT = 1.0 / tc[1];

  for (int i = 0; i < 2; ++i) {
    wdot[i] = 0.0;
  }

  return;
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(amrex::Real* T, amrex::Real* C, amrex::Real* wdot)
{
  int id; // loop counter

  // convert to SI
  for (id = 0; id < 2; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, *T);

  // convert to chemkin units
  for (id = 0; id < 2; ++id) {
    C[id] *= 1.0e-6;
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWYP(amrex::Real* P, amrex::Real* T, amrex::Real* y, amrex::Real* wdot)
{
  int id;           // loop counter
  amrex::Real c[2]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;
  amrex::Real imw[2];

  get_imw(imw);

  // Compute inverse of mean molecular wt first
  YOW += y[0] * imw[0]; // N2
  YOW += y[1] * imw[1]; // N2a
  // PW/RT (see Eq. 7)
  PWORT = (*P) / (YOW * 8.31446261815324e+07 * (*T));
  // multiply by 1e6 so c goes to SI
  PWORT *= 1e6;
  // Now compute conversion (and go to SI)
  c[0] = PWORT * y[0] * imw[0];
  c[1] = PWORT * y[1] * imw[1];

  // convert to chemkin units
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 2; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWXP(amrex::Real* P, amrex::Real* T, amrex::Real* x, amrex::Real* wdot)
{
  int id;           // loop counter
  amrex::Real c[2]; // temporary storage
  amrex::Real PORT =
    1e6 * (*P) /
    (8.31446261815324e+07 * (*T)); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (id = 0; id < 2; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 2; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWYR(amrex::Real* rho, amrex::Real* T, amrex::Real* y, amrex::Real* wdot)
{
  int id;           // loop counter
  amrex::Real c[2]; // temporary storage
  amrex::Real imw[2];

  get_imw(imw);

  // See Eq 8 with an extra 1e6 so c goes to SI
  c[0] = 1e6 * (*rho) * y[0] * imw[0];
  c[1] = 1e6 * (*rho) * y[1] * imw[1];

  // call productionRate
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 2; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWXR(amrex::Real* rho, amrex::Real* T, amrex::Real* x, amrex::Real* wdot)
{
  int id;             // loop counter
  amrex::Real c[2];   // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 28.014000; // N2
  XW += x[1] * 28.014000; // N2a
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * (*rho) / XW;

  // Compute conversion, see Eq 11
  for (id = 0; id < 2; ++id) {
    c[id] = x[id] * ROW;
  }

  // convert to chemkin units
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 2; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int* kcharge)
{
  kcharge[0] = 0; // N2
  kcharge[1] = 0; // N2a
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real* zk)
{
  amrex::Real imw[2];
  get_imw(imw);

  int kchrg[2];
  CKCHRG(kchrg);

  for (int id = 0; id < 2; ++id) {
    zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw[id];
  }
}

// compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dcvpRdT(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: N2
    species[0] = +1.40824000e-03 - 7.92644400e-06 * tc[1] +
                 1.69245450e-08 * tc[2] - 9.77942000e-12 * tc[3];
    // species 1: N2a
    species[1] = +1.40824000e-03 - 7.92644400e-06 * tc[1] +
                 1.69245450e-08 * tc[2] - 9.77942000e-12 * tc[3];
  } else {
    // species 0: N2
    species[0] = +1.48797700e-03 - 1.13695220e-06 * tc[1] +
                 3.02911200e-10 * tc[2] - 2.70134040e-14 * tc[3];
    // species 1: N2a
    species[1] = +1.48797700e-03 - 1.13695220e-06 * tc[1] +
                 3.02911200e-10 * tc[2] - 2.70134040e-14 * tc[3];
  }
  return;
}

// compute an approx to the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian_precond(amrex::Real* J, amrex::Real* sc, amrex::Real T, const int HP)
{

  for (int i = 0; i < 9; i++) {
    J[i] = 0.0;
  }

  amrex::Real wdot[2];
  for (int k = 0; k < 2; k++) {
    wdot[k] = 0.0;
  }

  const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
  amrex::Real invT = 1.0 / tc[1];
  amrex::Real invT2 = invT * invT;

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 / T;
  amrex::Real refCinv = 1.0 / refC;

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int k = 0; k < 2; ++k) {
    mixture += sc[k];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[2];
  gibbs(g_RT, tc);

  // compute the species enthalpy
  amrex::Real h_RT[2];
  speciesEnthalpy(h_RT, tc);

  amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  amrex::Real dqdci, dcdc_fac, dqdc[2];
  amrex::Real Pr, fPr, F, k_0, logPr;
  amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
  amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
  amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const amrex::Real ln10 = log(10.0);
  const amrex::Real log10e = 1.0 / log(10.0);
  amrex::Real c_R[2], dcRdT[2], e_RT[2];
  amrex::Real* eh_RT;
  if (HP) {
    cp_R(c_R, tc);
    dcvpRdT(dcRdT, tc);
    eh_RT = &h_RT[0];
  } else {
    cv_R(c_R, tc);
    dcvpRdT(dcRdT, tc);
    speciesInternalEnergy(e_RT, tc);
    eh_RT = &e_RT[0];
  }

  amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT = 0.0, dehmixdT = 0.0;
  for (int k = 0; k < 2; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[6 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 2; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 2; ++m) {
      dehmixdc += eh_RT[m] * J[k * 3 + m];
    }
    J[k * 3 + 2] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[8] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
}

// compute an approx to the reaction Jacobian (for preconditioning)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT_SIMPLIFIED(
  amrex::Real* J, amrex::Real* sc, amrex::Real* Tp, const int* HP)
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

  return;
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian(amrex::Real* J, amrex::Real* sc, amrex::Real T, const int consP)
{

  for (int i = 0; i < 9; i++) {
    J[i] = 0.0;
  }

  amrex::Real wdot[2];
  for (int k = 0; k < 2; k++) {
    wdot[k] = 0.0;
  }

  const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
  amrex::Real invT = 1.0 / tc[1];
  amrex::Real invT2 = invT * invT;

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 / T;
  amrex::Real refCinv = 1.0 / refC;

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int k = 0; k < 2; ++k) {
    mixture += sc[k];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[2];
  gibbs(g_RT, tc);

  // compute the species enthalpy
  amrex::Real h_RT[2];
  speciesEnthalpy(h_RT, tc);

  amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  amrex::Real dqdci, dcdc_fac, dqdc[2];
  amrex::Real Pr, fPr, F, k_0, logPr;
  amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
  amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
  amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const amrex::Real ln10 = log(10.0);
  const amrex::Real log10e = 1.0 / log(10.0);
  amrex::Real c_R[2], dcRdT[2], e_RT[2];
  amrex::Real* eh_RT;
  if (consP) {
    cp_R(c_R, tc);
    dcvpRdT(dcRdT, tc);
    eh_RT = &h_RT[0];
  } else {
    cv_R(c_R, tc);
    dcvpRdT(dcRdT, tc);
    speciesInternalEnergy(e_RT, tc);
    eh_RT = &e_RT[0];
  }

  amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT = 0.0, dehmixdT = 0.0;
  for (int k = 0; k < 2; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[6 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 2; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 2; ++m) {
      dehmixdc += eh_RT[m] * J[k * 3 + m];
    }
    J[k * 3 + 2] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[8] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT(amrex::Real* J, amrex::Real* sc, amrex::Real* Tp, const int* consP)
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

  return;
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
  WT[0] = 2.80140000E+01;
  WT[1] = 2.80140000E+01;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetEPS(amrex::Real* EPS)
{
  EPS[0] = 9.75300000E+01;
  EPS[1] = 9.75300000E+01;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetSIG(amrex::Real* SIG)
{
  SIG[0] = 3.62100000E+00;
  SIG[1] = 3.62100000E+00;
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
  POL[0] = 1.76000000E+00;
  POL[1] = 1.76000000E+00;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetZROT(amrex::Real* ZROT)
{
  ZROT[0] = 4.00000000E+00;
  ZROT[1] = 4.00000000E+00;
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
  COFETA[0] = -1.55270224E+01;
  COFETA[1] = 1.92766908E+00;
  COFETA[2] = -1.66518287E-01;
  COFETA[3] = 7.19100649E-03;
  COFETA[4] = -1.55270224E+01;
  COFETA[5] = 1.92766908E+00;
  COFETA[6] = -1.66518287E-01;
  COFETA[7] = 7.19100649E-03;
}

// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFLAM(amrex::Real* COFLAM)
{
  COFLAM[0] = 7.60995874E+00;
  COFLAM[1] = -1.18418483E+00;
  COFLAM[2] = 3.03558391E-01;
  COFLAM[3] = -1.54159450E-02;
  COFLAM[4] = 7.60995874E+00;
  COFLAM[5] = -1.18418483E+00;
  COFLAM[6] = 3.03558391E-01;
  COFLAM[7] = -1.54159450E-02;
}

// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFD(amrex::Real* COFD)
{
  COFD[0] = -1.42056777E+01;
  COFD[1] = 2.91297621E+00;
  COFD[2] = -1.61544771E-01;
  COFD[3] = 6.90271324E-03;
  COFD[4] = -1.42056777E+01;
  COFD[5] = 2.91297621E+00;
  COFD[6] = -1.61544771E-01;
  COFD[7] = 6.90271324E-03;
  COFD[8] = -1.42056777E+01;
  COFD[9] = 2.91297621E+00;
  COFD[10] = -1.61544771E-01;
  COFD[11] = 6.90271324E-03;
  COFD[12] = -1.42056777E+01;
  COFD[13] = 2.91297621E+00;
  COFD[14] = -1.61544771E-01;
  COFD[15] = 6.90271324E-03;
}

// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKTDIF(int* KTDIF)
{
}

// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFTD(amrex::Real* COFTD)
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
  amrex::Real avogadro = 6.02214199e23;
  amrex::Real boltzmann = 1.3806503e-16; // we work in CGS
  amrex::Real Rcst = 83.144598;          // in bar [CGS] !

  egtransetEPS(EPS);
  egtransetSIG(SIG);
  get_mw(wt);

  // species 0: N2
  // Imported from NIST
  Tci[0] = 126.192000;
  ai[0] = 1e6 * 0.42748 * Rcst * Rcst * Tci[0] * Tci[0] /
          (28.013400 * 28.013400 * 33.958000);
  bi[0] = 0.08664 * Rcst * Tci[0] / (28.013400 * 33.958000);
  acentric_i[0] = 0.037200;

  // species 1: N2a
  Tci[1] = 1.316 * EPS[1];
  ai[1] = (5.55 * avogadro * avogadro * EPS[1] * boltzmann * 1e-24 * SIG[1] *
           SIG[1] * SIG[1]) /
          (wt[1] * wt[1]);
  bi[1] = 0.855 * avogadro * 1e-24 * SIG[1] * SIG[1] * SIG[1] / (wt[1]);
  acentric_i[1] = 0.0;

  return;
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
