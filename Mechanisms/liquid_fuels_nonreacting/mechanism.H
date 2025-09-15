#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  O
1  H
2  C
3  N
4  Ar
5  He
*/

// Species
#define POSF10264_ID 0
#define POSF10289_ID 1
#define POSF10325_ID 2
#define POSF11498_ID 3
#define NC7H16_ID 4
#define NC10H22_ID 5
#define NC12H26_ID 6
#define IC8H18_ID 7
#define O2_ID 8
#define N2_ID 9
#define H2O_ID 10
#define CO2_ID 11

#define NUM_GAS_ELEMENTS 6  // Elements in the homogeneous phase
#define NUM_GAS_SPECIES 12  // Species in the homogeneous phase
#define NUM_GAS_REACTIONS 0 // Reactions in the homogeneous phase

#define SITE_DENSITY 0.000000E+00 // mol/cm^2

#define NUM_SURFACE_ELEMENTS 0  // Additional elements in heterogeneous phase
#define NUM_SURFACE_SPECIES 0   // Species in the heterogeneous phase
#define NUM_SURFACE_REACTIONS 0 // Reactions in the heterogeneous phase

#define NUM_ELEMENTS (NUM_GAS_ELEMENTS + NUM_SURFACE_ELEMENTS)
#define NUM_SPECIES (NUM_GAS_SPECIES + NUM_SURFACE_SPECIES)
#define NUM_REACTIONS (NUM_GAS_REACTIONS + NUM_SURFACE_REACTIONS)

#define NUM_IONS 0

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
  mm = 6;
  kk = 12;
  ii = 0;
  nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[12] = {
  6.4810074077914686e-03, // POSF10264
  5.9767147194530112e-03, // POSF10289
  6.4810074077914686e-03, // POSF10325
  5.4239641584448413e-03, // POSF11498
  9.9795419390249981e-03, // NC7H16
  7.0280983371519331e-03, // NC10H22
  5.8706117177409886e-03, // NC12H26
  8.7541144337838789e-03, // IC8H18
  3.1251953247077940e-02, // O2
  3.5696437495537948e-02, // N2
  5.5509297807382736e-02, // H2O
  2.2722624917630486e-02, // CO2
};
#endif
const amrex::Real h_global_imw[12] = {
  6.4810074077914686e-03, // POSF10264
  5.9767147194530112e-03, // POSF10289
  6.4810074077914686e-03, // POSF10325
  5.4239641584448413e-03, // POSF11498
  9.9795419390249981e-03, // NC7H16
  7.0280983371519331e-03, // NC10H22
  5.8706117177409886e-03, // NC12H26
  8.7541144337838789e-03, // IC8H18
  3.1251953247077940e-02, // O2
  3.5696437495537948e-02, // N2
  5.5509297807382736e-02, // H2O
  2.2722624917630486e-02, // CO2
};

//  molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[12] = {
  1.5429699999999997e+02, // POSF10264
  1.6731600000000000e+02, // POSF10289
  1.5429699999999997e+02, // POSF10325
  1.8436699999999999e+02, // POSF11498
  1.0020500000000000e+02, // NC7H16
  1.4228600000000000e+02, // NC10H22
  1.7034000000000000e+02, // NC12H26
  1.1423200000000000e+02, // IC8H18
  3.1998000000000001e+01, // O2
  2.8013999999999999e+01, // N2
  1.8015000000000001e+01, // H2O
  4.4009000000000000e+01, // CO2
};
#endif
const amrex::Real h_global_mw[12] = {
  1.5429699999999997e+02, // POSF10264
  1.6731600000000000e+02, // POSF10289
  1.5429699999999997e+02, // POSF10325
  1.8436699999999999e+02, // POSF11498
  1.0020500000000000e+02, // NC7H16
  1.4228600000000000e+02, // NC10H22
  1.7034000000000000e+02, // NC12H26
  1.1423200000000000e+02, // IC8H18
  3.1998000000000001e+01, // O2
  2.8013999999999999e+01, // N2
  1.8015000000000001e+01, // H2O
  4.4009000000000000e+01, // CO2
};

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_imw(amrex::Real* imw_new)
{
  imw_new[0] = 6.4810074077914686e-03;  // POSF10264
  imw_new[1] = 5.9767147194530112e-03;  // POSF10289
  imw_new[2] = 6.4810074077914686e-03;  // POSF10325
  imw_new[3] = 5.4239641584448413e-03;  // POSF11498
  imw_new[4] = 9.9795419390249981e-03;  // NC7H16
  imw_new[5] = 7.0280983371519331e-03;  // NC10H22
  imw_new[6] = 5.8706117177409886e-03;  // NC12H26
  imw_new[7] = 8.7541144337838789e-03;  // IC8H18
  imw_new[8] = 3.1251953247077940e-02;  // O2
  imw_new[9] = 3.5696437495537948e-02;  // N2
  imw_new[10] = 5.5509297807382736e-02; // H2O
  imw_new[11] = 2.2722624917630486e-02; // CO2
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
  mw_new[0] = 1.5429699999999997e+02;  // POSF10264
  mw_new[1] = 1.6731600000000000e+02;  // POSF10289
  mw_new[2] = 1.5429699999999997e+02;  // POSF10325
  mw_new[3] = 1.8436699999999999e+02;  // POSF11498
  mw_new[4] = 1.0020500000000000e+02;  // NC7H16
  mw_new[5] = 1.4228600000000000e+02;  // NC10H22
  mw_new[6] = 1.7034000000000000e+02;  // NC12H26
  mw_new[7] = 1.1423200000000000e+02;  // IC8H18
  mw_new[8] = 3.1998000000000001e+01;  // O2
  mw_new[9] = 2.8013999999999999e+01;  // N2
  mw_new[10] = 1.8015000000000001e+01; // H2O
  mw_new[11] = 4.4009000000000000e+01; // CO2
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10264
    species[0] = +3.70491270e+00 + 5.89113270e-02 * T + 1.05000140e-04 * T2 -
                 1.85010880e-07 * T3 + 8.22384310e-11 * T4;
    // species 1: POSF10289
    species[1] = -9.59735394e-01 + 8.87493420e-02 * T + 4.27242780e-05 * T2 -
                 1.17163540e-07 * T3 + 5.39969040e-11 * T4;
    // species 2: POSF10325
    species[2] = +1.57855180e+00 + 7.26248320e-02 * T + 7.32363700e-05 * T2 -
                 1.51154440e-07 * T3 + 6.89355790e-11 * T4;
    // species 3: POSF11498
    species[3] = +3.21777530e+00 + 9.30842010e-02 * T + 8.07107430e-05 * T2 -
                 1.74530840e-07 * T3 + 7.99705090e-11 * T4;
    // species 8: O2
    species[8] = +2.78245636e+00 - 2.99673416e-03 * T + 9.84730201e-06 * T2 -
                 9.68129509e-09 * T3 + 3.24372837e-12 * T4;
    // species 9: N2
    species[9] = +2.29867700e+00 + 1.40824040e-03 * T - 3.96322200e-06 * T2 +
                 5.64151500e-09 * T3 - 2.44485400e-12 * T4;
    // species 10: H2O
    species[10] = +3.19864056e+00 - 2.03643410e-03 * T + 6.52040211e-06 * T2 -
                  5.48797062e-09 * T3 + 1.77197817e-12 * T4;
    // species 11: CO2
    species[11] = +1.35677352e+00 + 8.98459677e-03 * T - 7.12356269e-06 * T2 +
                  2.45919022e-09 * T3 - 1.43699548e-13 * T4;
  } else {
    // species 0: POSF10264
    species[0] = +2.48974230e+01 + 5.54620920e-02 * T - 1.73377380e-05 * T2 +
                 1.75824520e-09 * T3 + 6.39718990e-14 * T4;
    // species 1: POSF10289
    species[1] = +1.91150530e+01 + 6.91818150e-02 * T - 2.31719970e-05 * T2 +
                 1.95197810e-09 * T3 + 2.70418640e-13 * T4;
    // species 2: POSF10325
    species[2] = +2.03578910e+01 + 6.65890870e-02 * T - 2.59674690e-05 * T2 +
                 4.47750010e-09 * T3 - 2.36071520e-13 * T4;
    // species 3: POSF11498
    species[3] = +3.26430400e+01 + 6.79280910e-02 * T - 2.05460000e-05 * T2 +
                 2.49695640e-09 * T3 - 6.96320510e-14 * T4;
    // species 8: O2
    species[8] = +2.28253784e+00 + 1.48308754e-03 * T - 7.57966669e-07 * T2 +
                 2.09470555e-10 * T3 - 2.16717794e-14 * T4;
    // species 9: N2
    species[9] = +1.92664000e+00 + 1.48797680e-03 * T - 5.68476000e-07 * T2 +
                 1.00970380e-10 * T3 - 6.75335100e-15 * T4;
    // species 10: H2O
    species[10] = +2.03399249e+00 + 2.17691804e-03 * T - 1.64072518e-07 * T2 -
                  9.70419870e-11 * T3 + 1.68200992e-14 * T4;
    // species 11: CO2
    species[11] = +2.85746029e+00 + 4.41437026e-03 * T - 2.21481404e-06 * T2 +
                  5.23490188e-10 * T3 - 4.72084164e-14 * T4;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 4: NC7H16
    species[4] = -2.26836187e+00 + 8.54355820e-02 * T - 5.25346786e-05 * T2 +
                 1.62945721e-08 * T3 - 2.02394925e-12 * T4;
    // species 5: NC10H22
    species[5] = -3.08416969e+00 + 1.22535012e-01 * T - 7.76815739e-05 * T2 +
                 2.49834877e-08 * T3 - 3.23548038e-12 * T4;
    // species 6: NC12H26
    species[6] = -3.62181594e+00 + 1.47237711e-01 * T - 9.43970271e-05 * T2 +
                 3.07441268e-08 * T3 - 4.03602230e-12 * T4;
  } else {
    // species 4: NC7H16
    species[4] = +2.12148969e+01 + 3.47675750e-02 * T - 1.18407129e-05 * T2 +
                 1.83298478e-09 * T3 - 1.06130266e-13 * T4;
    // species 5: NC10H22
    species[5] = +3.09882239e+01 + 4.77244922e-02 * T - 1.62276391e-05 * T2 +
                 2.50963259e-09 * T3 - 1.45215772e-13 * T4;
    // species 6: NC12H26
    species[6] = +3.75095037e+01 + 5.63550048e-02 * T - 1.91493200e-05 * T2 +
                 2.96024862e-09 * T3 - 1.71244150e-13 * T4;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 7: IC8H18
    species[7] = -5.20868893e+00 + 1.11440581e-01 * T - 7.91346582e-05 * T2 +
                 2.92406242e-08 * T3 - 4.43743191e-12 * T4;
  } else {
    // species 7: IC8H18
    species[7] = +2.61373590e+01 + 3.79004890e-02 * T - 1.29437358e-05 * T2 +
                 2.00760372e-09 * T3 - 1.16400580e-13 * T4;
  }
}

// compute Cp/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cp_R(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10264
    species[0] = +4.70491270e+00 + 5.89113270e-02 * T + 1.05000140e-04 * T2 -
                 1.85010880e-07 * T3 + 8.22384310e-11 * T4;
    // species 1: POSF10289
    species[1] = +4.02646060e-02 + 8.87493420e-02 * T + 4.27242780e-05 * T2 -
                 1.17163540e-07 * T3 + 5.39969040e-11 * T4;
    // species 2: POSF10325
    species[2] = +2.57855180e+00 + 7.26248320e-02 * T + 7.32363700e-05 * T2 -
                 1.51154440e-07 * T3 + 6.89355790e-11 * T4;
    // species 3: POSF11498
    species[3] = +4.21777530e+00 + 9.30842010e-02 * T + 8.07107430e-05 * T2 -
                 1.74530840e-07 * T3 + 7.99705090e-11 * T4;
    // species 8: O2
    species[8] = +3.78245636e+00 - 2.99673416e-03 * T + 9.84730201e-06 * T2 -
                 9.68129509e-09 * T3 + 3.24372837e-12 * T4;
    // species 9: N2
    species[9] = +3.29867700e+00 + 1.40824040e-03 * T - 3.96322200e-06 * T2 +
                 5.64151500e-09 * T3 - 2.44485400e-12 * T4;
    // species 10: H2O
    species[10] = +4.19864056e+00 - 2.03643410e-03 * T + 6.52040211e-06 * T2 -
                  5.48797062e-09 * T3 + 1.77197817e-12 * T4;
    // species 11: CO2
    species[11] = +2.35677352e+00 + 8.98459677e-03 * T - 7.12356269e-06 * T2 +
                  2.45919022e-09 * T3 - 1.43699548e-13 * T4;
  } else {
    // species 0: POSF10264
    species[0] = +2.58974230e+01 + 5.54620920e-02 * T - 1.73377380e-05 * T2 +
                 1.75824520e-09 * T3 + 6.39718990e-14 * T4;
    // species 1: POSF10289
    species[1] = +2.01150530e+01 + 6.91818150e-02 * T - 2.31719970e-05 * T2 +
                 1.95197810e-09 * T3 + 2.70418640e-13 * T4;
    // species 2: POSF10325
    species[2] = +2.13578910e+01 + 6.65890870e-02 * T - 2.59674690e-05 * T2 +
                 4.47750010e-09 * T3 - 2.36071520e-13 * T4;
    // species 3: POSF11498
    species[3] = +3.36430400e+01 + 6.79280910e-02 * T - 2.05460000e-05 * T2 +
                 2.49695640e-09 * T3 - 6.96320510e-14 * T4;
    // species 8: O2
    species[8] = +3.28253784e+00 + 1.48308754e-03 * T - 7.57966669e-07 * T2 +
                 2.09470555e-10 * T3 - 2.16717794e-14 * T4;
    // species 9: N2
    species[9] = +2.92664000e+00 + 1.48797680e-03 * T - 5.68476000e-07 * T2 +
                 1.00970380e-10 * T3 - 6.75335100e-15 * T4;
    // species 10: H2O
    species[10] = +3.03399249e+00 + 2.17691804e-03 * T - 1.64072518e-07 * T2 -
                  9.70419870e-11 * T3 + 1.68200992e-14 * T4;
    // species 11: CO2
    species[11] = +3.85746029e+00 + 4.41437026e-03 * T - 2.21481404e-06 * T2 +
                  5.23490188e-10 * T3 - 4.72084164e-14 * T4;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 4: NC7H16
    species[4] = -1.26836187e+00 + 8.54355820e-02 * T - 5.25346786e-05 * T2 +
                 1.62945721e-08 * T3 - 2.02394925e-12 * T4;
    // species 5: NC10H22
    species[5] = -2.08416969e+00 + 1.22535012e-01 * T - 7.76815739e-05 * T2 +
                 2.49834877e-08 * T3 - 3.23548038e-12 * T4;
    // species 6: NC12H26
    species[6] = -2.62181594e+00 + 1.47237711e-01 * T - 9.43970271e-05 * T2 +
                 3.07441268e-08 * T3 - 4.03602230e-12 * T4;
  } else {
    // species 4: NC7H16
    species[4] = +2.22148969e+01 + 3.47675750e-02 * T - 1.18407129e-05 * T2 +
                 1.83298478e-09 * T3 - 1.06130266e-13 * T4;
    // species 5: NC10H22
    species[5] = +3.19882239e+01 + 4.77244922e-02 * T - 1.62276391e-05 * T2 +
                 2.50963259e-09 * T3 - 1.45215772e-13 * T4;
    // species 6: NC12H26
    species[6] = +3.85095037e+01 + 5.63550048e-02 * T - 1.91493200e-05 * T2 +
                 2.96024862e-09 * T3 - 1.71244150e-13 * T4;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 7: IC8H18
    species[7] = -4.20868893e+00 + 1.11440581e-01 * T - 7.91346582e-05 * T2 +
                 2.92406242e-08 * T3 - 4.43743191e-12 * T4;
  } else {
    // species 7: IC8H18
    species[7] = +2.71373590e+01 + 3.79004890e-02 * T - 1.29437358e-05 * T2 +
                 2.00760372e-09 * T3 - 1.16400580e-13 * T4;
  }
}

// compute the g/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
gibbs(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10264
    species[0] = -3.739186300000000e+04 * invT - 8.734536299999998e+00 -
                 4.704912700000000e+00 * logT - 2.945566350000000e-02 * T -
                 1.750002333333333e-05 * T2 + 1.541757333333333e-08 * T3 -
                 4.111921550000000e-12 * T4;
    // species 1: POSF10289
    species[1] = -3.519053500000000e+04 * invT - 3.243793839400000e+01 -
                 4.026460600000000e-02 * logT - 4.437467100000000e-02 * T -
                 7.120713000000000e-06 * T2 + 9.763628333333333e-09 * T3 -
                 2.699845200000000e-12 * T4;
    // species 2: POSF10325
    species[2] = -3.802985500000000e+04 * invT - 1.995518320000000e+01 -
                 2.578551800000000e+00 * logT - 3.631241600000000e-02 * T -
                 1.220606166666667e-05 * T2 + 1.259620333333333e-08 * T3 -
                 3.446778950000000e-12 * T4;
    // species 3: POSF11498
    species[3] = -4.701783600000000e+04 * invT - 1.421897370000000e+01 -
                 4.217775300000000e+00 * logT - 4.654210050000000e-02 * T -
                 1.345179050000000e-05 * T2 + 1.454423666666667e-08 * T3 -
                 3.998525450000000e-12 * T4;
    // species 8: O2
    species[8] = -1.063943560000000e+03 * invT + 1.247806300000001e-01 -
                 3.782456360000000e+00 * logT + 1.498367080000000e-03 * T -
                 1.641217001666667e-06 * T2 + 8.067745908333334e-10 * T3 -
                 1.621864185000000e-13 * T4;
    // species 9: N2
    species[9] = -1.020899900000000e+03 * invT - 6.516950000000001e-01 -
                 3.298677000000000e+00 * logT - 7.041202000000000e-04 * T +
                 6.605369999999999e-07 * T2 - 4.701262500000001e-10 * T3 +
                 1.222427000000000e-13 * T4;
    // species 10: H2O
    species[10] = -3.029372670000000e+04 * invT + 5.047672768000000e+00 -
                  4.198640560000000e+00 * logT + 1.018217050000000e-03 * T -
                  1.086733685000000e-06 * T2 + 4.573308850000000e-10 * T3 -
                  8.859890850000000e-14 * T4;
    // species 11: CO2
    species[11] = -4.837196970000000e+04 * invT - 7.544278700000000e+00 -
                  2.356773520000000e+00 * logT - 4.492298385000000e-03 * T +
                  1.187260448333333e-06 * T2 - 2.049325183333333e-10 * T3 +
                  7.184977399999999e-15 * T4;
  } else {
    // species 0: POSF10264
    species[0] = -4.633780500000000e+04 * invT + 1.359452230000000e+02 -
                 2.589742300000000e+01 * logT - 2.773104600000000e-02 * T +
                 2.889623000000000e-06 * T2 - 1.465204333333333e-10 * T3 -
                 3.198594950000000e-15 * T4;
    // species 1: POSF10289
    species[1] = -4.254968400000000e+04 * invT + 1.000664460000000e+02 -
                 2.011505300000000e+01 * logT - 3.459090750000000e-02 * T +
                 3.861999500000000e-06 * T2 - 1.626648416666667e-10 * T3 -
                 1.352093200000000e-14 * T4;
    // species 2: POSF10325
    species[2] = -4.579699600000000e+04 * invT + 1.074939770000000e+02 -
                 2.135789100000000e+01 * logT - 3.329454350000000e-02 * T +
                 4.327911500000000e-06 * T2 - 3.731250083333333e-10 * T3 +
                 1.180357600000000e-14 * T4;
    // species 3: POSF11498
    species[3] = -5.836167200000000e+04 * invT + 1.816835200000000e+02 -
                 3.364304000000000e+01 * logT - 3.396404550000000e-02 * T +
                 3.424333333333333e-06 * T2 - 2.080797000000000e-10 * T3 +
                 3.481602550000000e-15 * T4;
    // species 8: O2
    species[8] = -1.088457720000000e+03 * invT - 2.170693450000000e+00 -
                 3.282537840000000e+00 * logT - 7.415437700000000e-04 * T +
                 1.263277781666667e-07 * T2 - 1.745587958333333e-11 * T3 +
                 1.083588970000000e-15 * T4;
    // species 9: N2
    species[9] = -9.227977000000000e+02 * invT - 3.053888000000000e+00 -
                 2.926640000000000e+00 * logT - 7.439884000000000e-04 * T +
                 9.474600000000001e-08 * T2 - 8.414198333333333e-12 * T3 +
                 3.376675500000000e-16 * T4;
    // species 10: H2O
    species[10] = -3.000429710000000e+04 * invT - 1.932777610000000e+00 -
                  3.033992490000000e+00 * logT - 1.088459020000000e-03 * T +
                  2.734541966666666e-08 * T2 + 8.086832250000000e-12 * T3 -
                  8.410049600000000e-16 * T4;
    // species 11: CO2
    species[11] = -4.875916600000000e+04 * invT + 1.585822230000000e+00 -
                  3.857460290000000e+00 * logT - 2.207185130000000e-03 * T +
                  3.691356733333334e-07 * T2 - 4.362418233333334e-11 * T3 +
                  2.360420820000000e-15 * T4;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 4: NC7H16
    species[4] = -2.565865650000000e+04 * invT - 3.664165307000000e+01 +
                 1.268361870000000e+00 * logT - 4.271779100000000e-02 * T +
                 8.755779766666667e-06 * T2 - 1.357881008333333e-09 * T3 +
                 1.011974625000000e-13 * T4;
    // species 5: NC10H22
    species[5] = -3.430218630000000e+04 * invT - 4.631018369000000e+01 +
                 2.084169690000000e+00 * logT - 6.126750600000000e-02 * T +
                 1.294692898333333e-05 * T2 - 2.081957308333334e-09 * T3 +
                 1.617740190000000e-13 * T4;
    // species 6: NC12H26
    species[6] = -4.006542530000000e+04 * invT - 5.272127854000000e+01 +
                 2.621815940000000e+00 * logT - 7.361885550000000e-02 * T +
                 1.573283785000000e-05 * T2 - 2.562010566666667e-09 * T3 +
                 2.018011150000000e-13 * T4;
  } else {
    // species 4: NC7H16
    species[4] = -3.427600810000000e+04 * invT + 1.145189165000000e+02 -
                 2.221489690000000e+01 * logT - 1.738378750000000e-02 * T +
                 1.973452150000000e-06 * T2 - 1.527487316666667e-10 * T3 +
                 5.306513300000000e-15 * T4;
    // species 5: NC10H22
    species[5] = -4.663928400000000e+04 * invT + 1.724923449000000e+02 -
                 3.198822390000000e+01 * logT - 2.386224610000000e-02 * T +
                 2.704606516666666e-06 * T2 - 2.091360491666667e-10 * T3 +
                 7.260788599999999e-15 * T4;
    // species 6: NC12H26
    species[6] = -5.488434650000000e+04 * invT + 2.111804257000000e+02 -
                 3.850950370000000e+01 * logT - 2.817750240000000e-02 * T +
                 3.191553333333333e-06 * T2 - 2.466873850000000e-10 * T3 +
                 8.562207500000000e-15 * T4;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 7: IC8H18
    species[7] = -2.994468750000000e+04 * invT - 4.916085903000000e+01 +
                 4.208688930000000e+00 * logT - 5.572029050000000e-02 * T +
                 1.318910970000000e-05 * T2 - 2.436718683333333e-09 * T3 +
                 2.218715955000000e-13 * T4;
  } else {
    // species 7: IC8H18
    species[7] = -4.079581770000000e+04 * invT + 1.504148540000000e+02 -
                 2.713735900000000e+01 * logT - 1.895024450000000e-02 * T +
                 2.157289300000000e-06 * T2 - 1.673003100000000e-10 * T3 +
                 5.820029000000000e-15 * T4;
  }
}

// compute the a/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
helmholtz(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10264
    species[0] = -3.73918630e+04 * invT - 9.73453630e+00 -
                 4.70491270e+00 * logT - 2.94556635e-02 * T -
                 1.75000233e-05 * T2 + 1.54175733e-08 * T3 -
                 4.11192155e-12 * T4;
    // species 1: POSF10289
    species[1] = -3.51905350e+04 * invT - 3.34379384e+01 -
                 4.02646060e-02 * logT - 4.43746710e-02 * T -
                 7.12071300e-06 * T2 + 9.76362833e-09 * T3 -
                 2.69984520e-12 * T4;
    // species 2: POSF10325
    species[2] = -3.80298550e+04 * invT - 2.09551832e+01 -
                 2.57855180e+00 * logT - 3.63124160e-02 * T -
                 1.22060617e-05 * T2 + 1.25962033e-08 * T3 -
                 3.44677895e-12 * T4;
    // species 3: POSF11498
    species[3] = -4.70178360e+04 * invT - 1.52189737e+01 -
                 4.21777530e+00 * logT - 4.65421005e-02 * T -
                 1.34517905e-05 * T2 + 1.45442367e-08 * T3 -
                 3.99852545e-12 * T4;
    // species 8: O2
    species[8] = -1.06394356e+03 * invT - 8.75219370e-01 -
                 3.78245636e+00 * logT + 1.49836708e-03 * T -
                 1.64121700e-06 * T2 + 8.06774591e-10 * T3 -
                 1.62186418e-13 * T4;
    // species 9: N2
    species[9] = -1.02089990e+03 * invT - 1.65169500e+00 -
                 3.29867700e+00 * logT - 7.04120200e-04 * T +
                 6.60537000e-07 * T2 - 4.70126250e-10 * T3 +
                 1.22242700e-13 * T4;
    // species 10: H2O
    species[10] = -3.02937267e+04 * invT + 4.04767277e+00 -
                  4.19864056e+00 * logT + 1.01821705e-03 * T -
                  1.08673369e-06 * T2 + 4.57330885e-10 * T3 -
                  8.85989085e-14 * T4;
    // species 11: CO2
    species[11] = -4.83719697e+04 * invT - 8.54427870e+00 -
                  2.35677352e+00 * logT - 4.49229839e-03 * T +
                  1.18726045e-06 * T2 - 2.04932518e-10 * T3 +
                  7.18497740e-15 * T4;
  } else {
    // species 0: POSF10264
    species[0] = -4.63378050e+04 * invT + 1.34945223e+02 -
                 2.58974230e+01 * logT - 2.77310460e-02 * T +
                 2.88962300e-06 * T2 - 1.46520433e-10 * T3 -
                 3.19859495e-15 * T4;
    // species 1: POSF10289
    species[1] = -4.25496840e+04 * invT + 9.90664460e+01 -
                 2.01150530e+01 * logT - 3.45909075e-02 * T +
                 3.86199950e-06 * T2 - 1.62664842e-10 * T3 -
                 1.35209320e-14 * T4;
    // species 2: POSF10325
    species[2] = -4.57969960e+04 * invT + 1.06493977e+02 -
                 2.13578910e+01 * logT - 3.32945435e-02 * T +
                 4.32791150e-06 * T2 - 3.73125008e-10 * T3 +
                 1.18035760e-14 * T4;
    // species 3: POSF11498
    species[3] = -5.83616720e+04 * invT + 1.80683520e+02 -
                 3.36430400e+01 * logT - 3.39640455e-02 * T +
                 3.42433333e-06 * T2 - 2.08079700e-10 * T3 +
                 3.48160255e-15 * T4;
    // species 8: O2
    species[8] = -1.08845772e+03 * invT - 3.17069345e+00 -
                 3.28253784e+00 * logT - 7.41543770e-04 * T +
                 1.26327778e-07 * T2 - 1.74558796e-11 * T3 +
                 1.08358897e-15 * T4;
    // species 9: N2
    species[9] = -9.22797700e+02 * invT - 4.05388800e+00 -
                 2.92664000e+00 * logT - 7.43988400e-04 * T +
                 9.47460000e-08 * T2 - 8.41419833e-12 * T3 +
                 3.37667550e-16 * T4;
    // species 10: H2O
    species[10] = -3.00042971e+04 * invT - 2.93277761e+00 -
                  3.03399249e+00 * logT - 1.08845902e-03 * T +
                  2.73454197e-08 * T2 + 8.08683225e-12 * T3 -
                  8.41004960e-16 * T4;
    // species 11: CO2
    species[11] = -4.87591660e+04 * invT + 5.85822230e-01 -
                  3.85746029e+00 * logT - 2.20718513e-03 * T +
                  3.69135673e-07 * T2 - 4.36241823e-11 * T3 +
                  2.36042082e-15 * T4;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 4: NC7H16
    species[4] = -2.56586565e+04 * invT - 3.76416531e+01 +
                 1.26836187e+00 * logT - 4.27177910e-02 * T +
                 8.75577977e-06 * T2 - 1.35788101e-09 * T3 +
                 1.01197462e-13 * T4;
    // species 5: NC10H22
    species[5] = -3.43021863e+04 * invT - 4.73101837e+01 +
                 2.08416969e+00 * logT - 6.12675060e-02 * T +
                 1.29469290e-05 * T2 - 2.08195731e-09 * T3 +
                 1.61774019e-13 * T4;
    // species 6: NC12H26
    species[6] = -4.00654253e+04 * invT - 5.37212785e+01 +
                 2.62181594e+00 * logT - 7.36188555e-02 * T +
                 1.57328378e-05 * T2 - 2.56201057e-09 * T3 +
                 2.01801115e-13 * T4;
  } else {
    // species 4: NC7H16
    species[4] = -3.42760081e+04 * invT + 1.13518917e+02 -
                 2.22148969e+01 * logT - 1.73837875e-02 * T +
                 1.97345215e-06 * T2 - 1.52748732e-10 * T3 +
                 5.30651330e-15 * T4;
    // species 5: NC10H22
    species[5] = -4.66392840e+04 * invT + 1.71492345e+02 -
                 3.19882239e+01 * logT - 2.38622461e-02 * T +
                 2.70460652e-06 * T2 - 2.09136049e-10 * T3 +
                 7.26078860e-15 * T4;
    // species 6: NC12H26
    species[6] = -5.48843465e+04 * invT + 2.10180426e+02 -
                 3.85095037e+01 * logT - 2.81775024e-02 * T +
                 3.19155333e-06 * T2 - 2.46687385e-10 * T3 +
                 8.56220750e-15 * T4;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 7: IC8H18
    species[7] = -2.99446875e+04 * invT - 5.01608590e+01 +
                 4.20868893e+00 * logT - 5.57202905e-02 * T +
                 1.31891097e-05 * T2 - 2.43671868e-09 * T3 +
                 2.21871595e-13 * T4;
  } else {
    // species 7: IC8H18
    species[7] = -4.07958177e+04 * invT + 1.49414854e+02 -
                 2.71373590e+01 * logT - 1.89502445e-02 * T +
                 2.15728930e-06 * T2 - 1.67300310e-10 * T3 +
                 5.82002900e-15 * T4;
  }
}

// compute the e/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesInternalEnergy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10264
    species[0] = +3.70491270e+00 + 2.94556635e-02 * T + 3.50000467e-05 * T2 -
                 4.62527200e-08 * T3 + 1.64476862e-11 * T4 -
                 3.73918630e+04 * invT;
    // species 1: POSF10289
    species[1] = -9.59735394e-01 + 4.43746710e-02 * T + 1.42414260e-05 * T2 -
                 2.92908850e-08 * T3 + 1.07993808e-11 * T4 -
                 3.51905350e+04 * invT;
    // species 2: POSF10325
    species[2] = +1.57855180e+00 + 3.63124160e-02 * T + 2.44121233e-05 * T2 -
                 3.77886100e-08 * T3 + 1.37871158e-11 * T4 -
                 3.80298550e+04 * invT;
    // species 3: POSF11498
    species[3] = +3.21777530e+00 + 4.65421005e-02 * T + 2.69035810e-05 * T2 -
                 4.36327100e-08 * T3 + 1.59941018e-11 * T4 -
                 4.70178360e+04 * invT;
    // species 8: O2
    species[8] = +2.78245636e+00 - 1.49836708e-03 * T + 3.28243400e-06 * T2 -
                 2.42032377e-09 * T3 + 6.48745674e-13 * T4 -
                 1.06394356e+03 * invT;
    // species 9: N2
    species[9] = +2.29867700e+00 + 7.04120200e-04 * T - 1.32107400e-06 * T2 +
                 1.41037875e-09 * T3 - 4.88970800e-13 * T4 -
                 1.02089990e+03 * invT;
    // species 10: H2O
    species[10] = +3.19864056e+00 - 1.01821705e-03 * T + 2.17346737e-06 * T2 -
                  1.37199266e-09 * T3 + 3.54395634e-13 * T4 -
                  3.02937267e+04 * invT;
    // species 11: CO2
    species[11] = +1.35677352e+00 + 4.49229839e-03 * T - 2.37452090e-06 * T2 +
                  6.14797555e-10 * T3 - 2.87399096e-14 * T4 -
                  4.83719697e+04 * invT;
  } else {
    // species 0: POSF10264
    species[0] = +2.48974230e+01 + 2.77310460e-02 * T - 5.77924600e-06 * T2 +
                 4.39561300e-10 * T3 + 1.27943798e-14 * T4 -
                 4.63378050e+04 * invT;
    // species 1: POSF10289
    species[1] = +1.91150530e+01 + 3.45909075e-02 * T - 7.72399900e-06 * T2 +
                 4.87994525e-10 * T3 + 5.40837280e-14 * T4 -
                 4.25496840e+04 * invT;
    // species 2: POSF10325
    species[2] = +2.03578910e+01 + 3.32945435e-02 * T - 8.65582300e-06 * T2 +
                 1.11937502e-09 * T3 - 4.72143040e-14 * T4 -
                 4.57969960e+04 * invT;
    // species 3: POSF11498
    species[3] = +3.26430400e+01 + 3.39640455e-02 * T - 6.84866667e-06 * T2 +
                 6.24239100e-10 * T3 - 1.39264102e-14 * T4 -
                 5.83616720e+04 * invT;
    // species 8: O2
    species[8] = +2.28253784e+00 + 7.41543770e-04 * T - 2.52655556e-07 * T2 +
                 5.23676387e-11 * T3 - 4.33435588e-15 * T4 -
                 1.08845772e+03 * invT;
    // species 9: N2
    species[9] = +1.92664000e+00 + 7.43988400e-04 * T - 1.89492000e-07 * T2 +
                 2.52425950e-11 * T3 - 1.35067020e-15 * T4 -
                 9.22797700e+02 * invT;
    // species 10: H2O
    species[10] = +2.03399249e+00 + 1.08845902e-03 * T - 5.46908393e-08 * T2 -
                  2.42604967e-11 * T3 + 3.36401984e-15 * T4 -
                  3.00042971e+04 * invT;
    // species 11: CO2
    species[11] = +2.85746029e+00 + 2.20718513e-03 * T - 7.38271347e-07 * T2 +
                  1.30872547e-10 * T3 - 9.44168328e-15 * T4 -
                  4.87591660e+04 * invT;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 4: NC7H16
    species[4] = -2.26836187e+00 + 4.27177910e-02 * T - 1.75115595e-05 * T2 +
                 4.07364302e-09 * T3 - 4.04789850e-13 * T4 -
                 2.56586565e+04 * invT;
    // species 5: NC10H22
    species[5] = -3.08416969e+00 + 6.12675060e-02 * T - 2.58938580e-05 * T2 +
                 6.24587193e-09 * T3 - 6.47096076e-13 * T4 -
                 3.43021863e+04 * invT;
    // species 6: NC12H26
    species[6] = -3.62181594e+00 + 7.36188555e-02 * T - 3.14656757e-05 * T2 +
                 7.68603170e-09 * T3 - 8.07204460e-13 * T4 -
                 4.00654253e+04 * invT;
  } else {
    // species 4: NC7H16
    species[4] = +2.12148969e+01 + 1.73837875e-02 * T - 3.94690430e-06 * T2 +
                 4.58246195e-10 * T3 - 2.12260532e-14 * T4 -
                 3.42760081e+04 * invT;
    // species 5: NC10H22
    species[5] = +3.09882239e+01 + 2.38622461e-02 * T - 5.40921303e-06 * T2 +
                 6.27408147e-10 * T3 - 2.90431544e-14 * T4 -
                 4.66392840e+04 * invT;
    // species 6: NC12H26
    species[6] = +3.75095037e+01 + 2.81775024e-02 * T - 6.38310667e-06 * T2 +
                 7.40062155e-10 * T3 - 3.42488300e-14 * T4 -
                 5.48843465e+04 * invT;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 7: IC8H18
    species[7] = -5.20868893e+00 + 5.57202905e-02 * T - 2.63782194e-05 * T2 +
                 7.31015605e-09 * T3 - 8.87486382e-13 * T4 -
                 2.99446875e+04 * invT;
  } else {
    // species 7: IC8H18
    species[7] = +2.61373590e+01 + 1.89502445e-02 * T - 4.31457860e-06 * T2 +
                 5.01900930e-10 * T3 - 2.32801160e-14 * T4 -
                 4.07958177e+04 * invT;
  }
}

// compute the h/(RT) at the given temperature (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEnthalpy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10264
    species[0] = +4.70491270e+00 + 2.94556635e-02 * T + 3.50000467e-05 * T2 -
                 4.62527200e-08 * T3 + 1.64476862e-11 * T4 -
                 3.73918630e+04 * invT;
    // species 1: POSF10289
    species[1] = +4.02646060e-02 + 4.43746710e-02 * T + 1.42414260e-05 * T2 -
                 2.92908850e-08 * T3 + 1.07993808e-11 * T4 -
                 3.51905350e+04 * invT;
    // species 2: POSF10325
    species[2] = +2.57855180e+00 + 3.63124160e-02 * T + 2.44121233e-05 * T2 -
                 3.77886100e-08 * T3 + 1.37871158e-11 * T4 -
                 3.80298550e+04 * invT;
    // species 3: POSF11498
    species[3] = +4.21777530e+00 + 4.65421005e-02 * T + 2.69035810e-05 * T2 -
                 4.36327100e-08 * T3 + 1.59941018e-11 * T4 -
                 4.70178360e+04 * invT;
    // species 8: O2
    species[8] = +3.78245636e+00 - 1.49836708e-03 * T + 3.28243400e-06 * T2 -
                 2.42032377e-09 * T3 + 6.48745674e-13 * T4 -
                 1.06394356e+03 * invT;
    // species 9: N2
    species[9] = +3.29867700e+00 + 7.04120200e-04 * T - 1.32107400e-06 * T2 +
                 1.41037875e-09 * T3 - 4.88970800e-13 * T4 -
                 1.02089990e+03 * invT;
    // species 10: H2O
    species[10] = +4.19864056e+00 - 1.01821705e-03 * T + 2.17346737e-06 * T2 -
                  1.37199266e-09 * T3 + 3.54395634e-13 * T4 -
                  3.02937267e+04 * invT;
    // species 11: CO2
    species[11] = +2.35677352e+00 + 4.49229839e-03 * T - 2.37452090e-06 * T2 +
                  6.14797555e-10 * T3 - 2.87399096e-14 * T4 -
                  4.83719697e+04 * invT;
  } else {
    // species 0: POSF10264
    species[0] = +2.58974230e+01 + 2.77310460e-02 * T - 5.77924600e-06 * T2 +
                 4.39561300e-10 * T3 + 1.27943798e-14 * T4 -
                 4.63378050e+04 * invT;
    // species 1: POSF10289
    species[1] = +2.01150530e+01 + 3.45909075e-02 * T - 7.72399900e-06 * T2 +
                 4.87994525e-10 * T3 + 5.40837280e-14 * T4 -
                 4.25496840e+04 * invT;
    // species 2: POSF10325
    species[2] = +2.13578910e+01 + 3.32945435e-02 * T - 8.65582300e-06 * T2 +
                 1.11937502e-09 * T3 - 4.72143040e-14 * T4 -
                 4.57969960e+04 * invT;
    // species 3: POSF11498
    species[3] = +3.36430400e+01 + 3.39640455e-02 * T - 6.84866667e-06 * T2 +
                 6.24239100e-10 * T3 - 1.39264102e-14 * T4 -
                 5.83616720e+04 * invT;
    // species 8: O2
    species[8] = +3.28253784e+00 + 7.41543770e-04 * T - 2.52655556e-07 * T2 +
                 5.23676387e-11 * T3 - 4.33435588e-15 * T4 -
                 1.08845772e+03 * invT;
    // species 9: N2
    species[9] = +2.92664000e+00 + 7.43988400e-04 * T - 1.89492000e-07 * T2 +
                 2.52425950e-11 * T3 - 1.35067020e-15 * T4 -
                 9.22797700e+02 * invT;
    // species 10: H2O
    species[10] = +3.03399249e+00 + 1.08845902e-03 * T - 5.46908393e-08 * T2 -
                  2.42604967e-11 * T3 + 3.36401984e-15 * T4 -
                  3.00042971e+04 * invT;
    // species 11: CO2
    species[11] = +3.85746029e+00 + 2.20718513e-03 * T - 7.38271347e-07 * T2 +
                  1.30872547e-10 * T3 - 9.44168328e-15 * T4 -
                  4.87591660e+04 * invT;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 4: NC7H16
    species[4] = -1.26836187e+00 + 4.27177910e-02 * T - 1.75115595e-05 * T2 +
                 4.07364302e-09 * T3 - 4.04789850e-13 * T4 -
                 2.56586565e+04 * invT;
    // species 5: NC10H22
    species[5] = -2.08416969e+00 + 6.12675060e-02 * T - 2.58938580e-05 * T2 +
                 6.24587193e-09 * T3 - 6.47096076e-13 * T4 -
                 3.43021863e+04 * invT;
    // species 6: NC12H26
    species[6] = -2.62181594e+00 + 7.36188555e-02 * T - 3.14656757e-05 * T2 +
                 7.68603170e-09 * T3 - 8.07204460e-13 * T4 -
                 4.00654253e+04 * invT;
  } else {
    // species 4: NC7H16
    species[4] = +2.22148969e+01 + 1.73837875e-02 * T - 3.94690430e-06 * T2 +
                 4.58246195e-10 * T3 - 2.12260532e-14 * T4 -
                 3.42760081e+04 * invT;
    // species 5: NC10H22
    species[5] = +3.19882239e+01 + 2.38622461e-02 * T - 5.40921303e-06 * T2 +
                 6.27408147e-10 * T3 - 2.90431544e-14 * T4 -
                 4.66392840e+04 * invT;
    // species 6: NC12H26
    species[6] = +3.85095037e+01 + 2.81775024e-02 * T - 6.38310667e-06 * T2 +
                 7.40062155e-10 * T3 - 3.42488300e-14 * T4 -
                 5.48843465e+04 * invT;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 7: IC8H18
    species[7] = -4.20868893e+00 + 5.57202905e-02 * T - 2.63782194e-05 * T2 +
                 7.31015605e-09 * T3 - 8.87486382e-13 * T4 -
                 2.99446875e+04 * invT;
  } else {
    // species 7: IC8H18
    species[7] = +2.71373590e+01 + 1.89502445e-02 * T - 4.31457860e-06 * T2 +
                 5.01900930e-10 * T3 - 2.32801160e-14 * T4 -
                 4.07958177e+04 * invT;
  }
}

// compute the S/R at the given temperature (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEntropy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10264
    species[0] = +4.70491270e+00 * logT + 5.89113270e-02 * T +
                 5.25000700e-05 * T2 - 6.16702933e-08 * T3 +
                 2.05596077e-11 * T4 + 1.34394490e+01;
    // species 1: POSF10289
    species[1] = +4.02646060e-02 * logT + 8.87493420e-02 * T +
                 2.13621390e-05 * T2 - 3.90545133e-08 * T3 +
                 1.34992260e-11 * T4 + 3.24782030e+01;
    // species 2: POSF10325
    species[2] = +2.57855180e+00 * logT + 7.26248320e-02 * T +
                 3.66181850e-05 * T2 - 5.03848133e-08 * T3 +
                 1.72338947e-11 * T4 + 2.25337350e+01;
    // species 3: POSF11498
    species[3] = +4.21777530e+00 * logT + 9.30842010e-02 * T +
                 4.03553715e-05 * T2 - 5.81769467e-08 * T3 +
                 1.99926272e-11 * T4 + 1.84367490e+01;
    // species 8: O2
    species[8] = +3.78245636e+00 * logT - 2.99673416e-03 * T +
                 4.92365101e-06 * T2 - 3.22709836e-09 * T3 +
                 8.10932092e-13 * T4 + 3.65767573e+00;
    // species 9: N2
    species[9] = +3.29867700e+00 * logT + 1.40824040e-03 * T -
                 1.98161100e-06 * T2 + 1.88050500e-09 * T3 -
                 6.11213500e-13 * T4 + 3.95037200e+00;
    // species 10: H2O
    species[10] = +4.19864056e+00 * logT - 2.03643410e-03 * T +
                  3.26020105e-06 * T2 - 1.82932354e-09 * T3 +
                  4.42994543e-13 * T4 - 8.49032208e-01;
    // species 11: CO2
    species[11] = +2.35677352e+00 * logT + 8.98459677e-03 * T -
                  3.56178134e-06 * T2 + 8.19730073e-10 * T3 -
                  3.59248870e-14 * T4 + 9.90105222e+00;
  } else {
    // species 0: POSF10264
    species[0] = +2.58974230e+01 * logT + 5.54620920e-02 * T -
                 8.66886900e-06 * T2 + 5.86081733e-10 * T3 +
                 1.59929747e-14 * T4 - 1.10047800e+02;
    // species 1: POSF10289
    species[1] = +2.01150530e+01 * logT + 6.91818150e-02 * T -
                 1.15859985e-05 * T2 + 6.50659367e-10 * T3 +
                 6.76046600e-14 * T4 - 7.99513930e+01;
    // species 2: POSF10325
    species[2] = +2.13578910e+01 * logT + 6.65890870e-02 * T -
                 1.29837345e-05 * T2 + 1.49250003e-09 * T3 -
                 5.90178800e-14 * T4 - 8.61360860e+01;
    // species 3: POSF11498
    species[3] = +3.36430400e+01 * logT + 6.79280910e-02 * T -
                 1.02730000e-05 * T2 + 8.32318800e-10 * T3 -
                 1.74080128e-14 * T4 - 1.48040480e+02;
    // species 8: O2
    species[8] = +3.28253784e+00 * logT + 1.48308754e-03 * T -
                 3.78983334e-07 * T2 + 6.98235183e-11 * T3 -
                 5.41794485e-15 * T4 + 5.45323129e+00;
    // species 9: N2
    species[9] = +2.92664000e+00 * logT + 1.48797680e-03 * T -
                 2.84238000e-07 * T2 + 3.36567933e-11 * T3 -
                 1.68833775e-15 * T4 + 5.98052800e+00;
    // species 10: H2O
    species[10] = +3.03399249e+00 * logT + 2.17691804e-03 * T -
                  8.20362590e-08 * T2 - 3.23473290e-11 * T3 +
                  4.20502480e-15 * T4 + 4.96677010e+00;
    // species 11: CO2
    species[11] = +3.85746029e+00 * logT + 4.41437026e-03 * T -
                  1.10740702e-06 * T2 + 1.74496729e-10 * T3 -
                  1.18021041e-14 * T4 + 2.27163806e+00;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 4: NC7H16
    species[4] = -1.26836187e+00 * logT + 8.54355820e-02 * T -
                 2.62673393e-05 * T2 + 5.43152403e-09 * T3 -
                 5.05987313e-13 * T4 + 3.53732912e+01;
    // species 5: NC10H22
    species[5] = -2.08416969e+00 * logT + 1.22535012e-01 * T -
                 3.88407869e-05 * T2 + 8.32782923e-09 * T3 -
                 8.08870095e-13 * T4 + 4.42260140e+01;
    // species 6: NC12H26
    species[6] = -2.62181594e+00 * logT + 1.47237711e-01 * T -
                 4.71985136e-05 * T2 + 1.02480423e-08 * T3 -
                 1.00900558e-12 * T4 + 5.00994626e+01;
  } else {
    // species 4: NC7H16
    species[4] = +2.22148969e+01 * logT + 3.47675750e-02 * T -
                 5.92035645e-06 * T2 + 6.10994927e-10 * T3 -
                 2.65325665e-14 * T4 - 9.23040196e+01;
    // species 5: NC10H22
    species[5] = +3.19882239e+01 * logT + 4.77244922e-02 * T -
                 8.11381955e-06 * T2 + 8.36544197e-10 * T3 -
                 3.63039430e-14 * T4 - 1.40504121e+02;
    // species 6: NC12H26
    species[6] = +3.85095037e+01 * logT + 5.63550048e-02 * T -
                 9.57466000e-06 * T2 + 9.86749540e-10 * T3 -
                 4.28110375e-14 * T4 - 1.72670922e+02;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 7: IC8H18
    species[7] = -4.20868893e+00 * logT + 1.11440581e-01 * T -
                 3.95673291e-05 * T2 + 9.74687473e-09 * T3 -
                 1.10935798e-12 * T4 + 4.49521701e+01;
  } else {
    // species 7: IC8H18
    species[7] = +2.71373590e+01 * logT + 3.79004890e-02 * T -
                 6.47186790e-06 * T2 + 6.69201240e-10 * T3 -
                 2.91001450e-14 * T4 - 1.23277495e+02;
  }
}

// compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dcvpRdT(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T2;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10264
    species[0] = +5.89113270e-02 + 2.10000280e-04 * T - 5.55032640e-07 * T2 +
                 3.28953724e-10 * T3;
    // species 1: POSF10289
    species[1] = +8.87493420e-02 + 8.54485560e-05 * T - 3.51490620e-07 * T2 +
                 2.15987616e-10 * T3;
    // species 2: POSF10325
    species[2] = +7.26248320e-02 + 1.46472740e-04 * T - 4.53463320e-07 * T2 +
                 2.75742316e-10 * T3;
    // species 3: POSF11498
    species[3] = +9.30842010e-02 + 1.61421486e-04 * T - 5.23592520e-07 * T2 +
                 3.19882036e-10 * T3;
    // species 8: O2
    species[8] = -2.99673416e-03 + 1.96946040e-05 * T - 2.90438853e-08 * T2 +
                 1.29749135e-11 * T3;
    // species 9: N2
    species[9] = +1.40824040e-03 - 7.92644400e-06 * T + 1.69245450e-08 * T2 -
                 9.77941600e-12 * T3;
    // species 10: H2O
    species[10] = -2.03643410e-03 + 1.30408042e-05 * T - 1.64639119e-08 * T2 +
                  7.08791268e-12 * T3;
    // species 11: CO2
    species[11] = +8.98459677e-03 - 1.42471254e-05 * T + 7.37757066e-09 * T2 -
                  5.74798192e-13 * T3;
  } else {
    // species 0: POSF10264
    species[0] = +5.54620920e-02 - 3.46754760e-05 * T + 5.27473560e-09 * T2 +
                 2.55887596e-13 * T3;
    // species 1: POSF10289
    species[1] = +6.91818150e-02 - 4.63439940e-05 * T + 5.85593430e-09 * T2 +
                 1.08167456e-12 * T3;
    // species 2: POSF10325
    species[2] = +6.65890870e-02 - 5.19349380e-05 * T + 1.34325003e-08 * T2 -
                 9.44286080e-13 * T3;
    // species 3: POSF11498
    species[3] = +6.79280910e-02 - 4.10920000e-05 * T + 7.49086920e-09 * T2 -
                 2.78528204e-13 * T3;
    // species 8: O2
    species[8] = +1.48308754e-03 - 1.51593334e-06 * T + 6.28411665e-10 * T2 -
                 8.66871176e-14 * T3;
    // species 9: N2
    species[9] = +1.48797680e-03 - 1.13695200e-06 * T + 3.02911140e-10 * T2 -
                 2.70134040e-14 * T3;
    // species 10: H2O
    species[10] = +2.17691804e-03 - 3.28145036e-07 * T - 2.91125961e-10 * T2 +
                  6.72803968e-14 * T3;
    // species 11: CO2
    species[11] = +4.41437026e-03 - 4.42962808e-06 * T + 1.57047056e-09 * T2 -
                  1.88833666e-13 * T3;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 4: NC7H16
    species[4] = +8.54355820e-02 - 1.05069357e-04 * T + 4.88837163e-08 * T2 -
                 8.09579700e-12 * T3;
    // species 5: NC10H22
    species[5] = +1.22535012e-01 - 1.55363148e-04 * T + 7.49504631e-08 * T2 -
                 1.29419215e-11 * T3;
    // species 6: NC12H26
    species[6] = +1.47237711e-01 - 1.88794054e-04 * T + 9.22323804e-08 * T2 -
                 1.61440892e-11 * T3;
  } else {
    // species 4: NC7H16
    species[4] = +3.47675750e-02 - 2.36814258e-05 * T + 5.49895434e-09 * T2 -
                 4.24521064e-13 * T3;
    // species 5: NC10H22
    species[5] = +4.77244922e-02 - 3.24552782e-05 * T + 7.52889777e-09 * T2 -
                 5.80863088e-13 * T3;
    // species 6: NC12H26
    species[6] = +5.63550048e-02 - 3.82986400e-05 * T + 8.88074586e-09 * T2 -
                 6.84976600e-13 * T3;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 7: IC8H18
    species[7] = +1.11440581e-01 - 1.58269316e-04 * T + 8.77218726e-08 * T2 -
                 1.77497276e-11 * T3;
  } else {
    // species 7: IC8H18
    species[7] = +3.79004890e-02 - 2.58874716e-05 * T + 6.02281116e-09 * T2 -
                 4.65602320e-13 * T3;
  }
}

// Returns the mean specific heat at CP (Eq. 33)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cpbl)
{
  amrex::Real result = 0;
  amrex::Real cpor[12]; //  temporary storage
  cp_R(cpor, T);

  // perform dot product
  for (int id = 0; id < 12; ++id) {
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10264
    result += y[0] *
              (+4.70491270e+00 + 5.89113270e-02 * T + 1.05000140e-04 * T2 -
               1.85010880e-07 * T3 + 8.22384310e-11 * T4) *
              0.0064810074077915;
    // species 1: POSF10289
    result += y[1] *
              (+4.02646060e-02 + 8.87493420e-02 * T + 4.27242780e-05 * T2 -
               1.17163540e-07 * T3 + 5.39969040e-11 * T4) *
              0.0059767147194530;
    // species 2: POSF10325
    result += y[2] *
              (+2.57855180e+00 + 7.26248320e-02 * T + 7.32363700e-05 * T2 -
               1.51154440e-07 * T3 + 6.89355790e-11 * T4) *
              0.0064810074077915;
    // species 3: POSF11498
    result += y[3] *
              (+4.21777530e+00 + 9.30842010e-02 * T + 8.07107430e-05 * T2 -
               1.74530840e-07 * T3 + 7.99705090e-11 * T4) *
              0.0054239641584448;
    // species 8: O2
    result += y[8] *
              (+3.78245636e+00 - 2.99673416e-03 * T + 9.84730201e-06 * T2 -
               9.68129509e-09 * T3 + 3.24372837e-12 * T4) *
              0.0312519532470779;
    // species 9: N2
    result += y[9] *
              (+3.29867700e+00 + 1.40824040e-03 * T - 3.96322200e-06 * T2 +
               5.64151500e-09 * T3 - 2.44485400e-12 * T4) *
              0.0356964374955379;
    // species 10: H2O
    result += y[10] *
              (+4.19864056e+00 - 2.03643410e-03 * T + 6.52040211e-06 * T2 -
               5.48797062e-09 * T3 + 1.77197817e-12 * T4) *
              0.0555092978073827;
    // species 11: CO2
    result += y[11] *
              (+2.35677352e+00 + 8.98459677e-03 * T - 7.12356269e-06 * T2 +
               2.45919022e-09 * T3 - 1.43699548e-13 * T4) *
              0.0227226249176305;
  } else {
    // species 0: POSF10264
    result += y[0] *
              (+2.58974230e+01 + 5.54620920e-02 * T - 1.73377380e-05 * T2 +
               1.75824520e-09 * T3 + 6.39718990e-14 * T4) *
              0.0064810074077915;
    // species 1: POSF10289
    result += y[1] *
              (+2.01150530e+01 + 6.91818150e-02 * T - 2.31719970e-05 * T2 +
               1.95197810e-09 * T3 + 2.70418640e-13 * T4) *
              0.0059767147194530;
    // species 2: POSF10325
    result += y[2] *
              (+2.13578910e+01 + 6.65890870e-02 * T - 2.59674690e-05 * T2 +
               4.47750010e-09 * T3 - 2.36071520e-13 * T4) *
              0.0064810074077915;
    // species 3: POSF11498
    result += y[3] *
              (+3.36430400e+01 + 6.79280910e-02 * T - 2.05460000e-05 * T2 +
               2.49695640e-09 * T3 - 6.96320510e-14 * T4) *
              0.0054239641584448;
    // species 8: O2
    result += y[8] *
              (+3.28253784e+00 + 1.48308754e-03 * T - 7.57966669e-07 * T2 +
               2.09470555e-10 * T3 - 2.16717794e-14 * T4) *
              0.0312519532470779;
    // species 9: N2
    result += y[9] *
              (+2.92664000e+00 + 1.48797680e-03 * T - 5.68476000e-07 * T2 +
               1.00970380e-10 * T3 - 6.75335100e-15 * T4) *
              0.0356964374955379;
    // species 10: H2O
    result += y[10] *
              (+3.03399249e+00 + 2.17691804e-03 * T - 1.64072518e-07 * T2 -
               9.70419870e-11 * T3 + 1.68200992e-14 * T4) *
              0.0555092978073827;
    // species 11: CO2
    result += y[11] *
              (+3.85746029e+00 + 4.41437026e-03 * T - 2.21481404e-06 * T2 +
               5.23490188e-10 * T3 - 4.72084164e-14 * T4) *
              0.0227226249176305;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 4: NC7H16
    result += y[4] *
              (-1.26836187e+00 + 8.54355820e-02 * T - 5.25346786e-05 * T2 +
               1.62945721e-08 * T3 - 2.02394925e-12 * T4) *
              0.0099795419390250;
    // species 5: NC10H22
    result += y[5] *
              (-2.08416969e+00 + 1.22535012e-01 * T - 7.76815739e-05 * T2 +
               2.49834877e-08 * T3 - 3.23548038e-12 * T4) *
              0.0070280983371519;
    // species 6: NC12H26
    result += y[6] *
              (-2.62181594e+00 + 1.47237711e-01 * T - 9.43970271e-05 * T2 +
               3.07441268e-08 * T3 - 4.03602230e-12 * T4) *
              0.0058706117177410;
  } else {
    // species 4: NC7H16
    result += y[4] *
              (+2.22148969e+01 + 3.47675750e-02 * T - 1.18407129e-05 * T2 +
               1.83298478e-09 * T3 - 1.06130266e-13 * T4) *
              0.0099795419390250;
    // species 5: NC10H22
    result += y[5] *
              (+3.19882239e+01 + 4.77244922e-02 * T - 1.62276391e-05 * T2 +
               2.50963259e-09 * T3 - 1.45215772e-13 * T4) *
              0.0070280983371519;
    // species 6: NC12H26
    result += y[6] *
              (+3.85095037e+01 + 5.63550048e-02 * T - 1.91493200e-05 * T2 +
               2.96024862e-09 * T3 - 1.71244150e-13 * T4) *
              0.0058706117177410;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 7: IC8H18
    result += y[7] *
              (-4.20868893e+00 + 1.11440581e-01 * T - 7.91346582e-05 * T2 +
               2.92406242e-08 * T3 - 4.43743191e-12 * T4) *
              0.0087541144337839;
  } else {
    // species 7: IC8H18
    result += y[7] *
              (+2.71373590e+01 + 3.79004890e-02 * T - 1.29437358e-05 * T2 +
               2.00760372e-09 * T3 - 1.16400580e-13 * T4) *
              0.0087541144337839;
  }

  cpbs = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 35)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cvbl)
{
  amrex::Real result = 0;
  amrex::Real cvor[12]; //  temporary storage
  cv_R(cvor, T);

  // perform dot product
  for (int id = 0; id < 12; ++id) {
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10264
    result += y[0] *
              (+3.70491270e+00 + 5.89113270e-02 * T + 1.05000140e-04 * T2 -
               1.85010880e-07 * T3 + 8.22384310e-11 * T4) *
              0.0064810074077915;
    // species 1: POSF10289
    result += y[1] *
              (-9.59735394e-01 + 8.87493420e-02 * T + 4.27242780e-05 * T2 -
               1.17163540e-07 * T3 + 5.39969040e-11 * T4) *
              0.0059767147194530;
    // species 2: POSF10325
    result += y[2] *
              (+1.57855180e+00 + 7.26248320e-02 * T + 7.32363700e-05 * T2 -
               1.51154440e-07 * T3 + 6.89355790e-11 * T4) *
              0.0064810074077915;
    // species 3: POSF11498
    result += y[3] *
              (+3.21777530e+00 + 9.30842010e-02 * T + 8.07107430e-05 * T2 -
               1.74530840e-07 * T3 + 7.99705090e-11 * T4) *
              0.0054239641584448;
    // species 8: O2
    result += y[8] *
              (+2.78245636e+00 - 2.99673416e-03 * T + 9.84730201e-06 * T2 -
               9.68129509e-09 * T3 + 3.24372837e-12 * T4) *
              0.0312519532470779;
    // species 9: N2
    result += y[9] *
              (+2.29867700e+00 + 1.40824040e-03 * T - 3.96322200e-06 * T2 +
               5.64151500e-09 * T3 - 2.44485400e-12 * T4) *
              0.0356964374955379;
    // species 10: H2O
    result += y[10] *
              (+3.19864056e+00 - 2.03643410e-03 * T + 6.52040211e-06 * T2 -
               5.48797062e-09 * T3 + 1.77197817e-12 * T4) *
              0.0555092978073827;
    // species 11: CO2
    result += y[11] *
              (+1.35677352e+00 + 8.98459677e-03 * T - 7.12356269e-06 * T2 +
               2.45919022e-09 * T3 - 1.43699548e-13 * T4) *
              0.0227226249176305;
  } else {
    // species 0: POSF10264
    result += y[0] *
              (+2.48974230e+01 + 5.54620920e-02 * T - 1.73377380e-05 * T2 +
               1.75824520e-09 * T3 + 6.39718990e-14 * T4) *
              0.0064810074077915;
    // species 1: POSF10289
    result += y[1] *
              (+1.91150530e+01 + 6.91818150e-02 * T - 2.31719970e-05 * T2 +
               1.95197810e-09 * T3 + 2.70418640e-13 * T4) *
              0.0059767147194530;
    // species 2: POSF10325
    result += y[2] *
              (+2.03578910e+01 + 6.65890870e-02 * T - 2.59674690e-05 * T2 +
               4.47750010e-09 * T3 - 2.36071520e-13 * T4) *
              0.0064810074077915;
    // species 3: POSF11498
    result += y[3] *
              (+3.26430400e+01 + 6.79280910e-02 * T - 2.05460000e-05 * T2 +
               2.49695640e-09 * T3 - 6.96320510e-14 * T4) *
              0.0054239641584448;
    // species 8: O2
    result += y[8] *
              (+2.28253784e+00 + 1.48308754e-03 * T - 7.57966669e-07 * T2 +
               2.09470555e-10 * T3 - 2.16717794e-14 * T4) *
              0.0312519532470779;
    // species 9: N2
    result += y[9] *
              (+1.92664000e+00 + 1.48797680e-03 * T - 5.68476000e-07 * T2 +
               1.00970380e-10 * T3 - 6.75335100e-15 * T4) *
              0.0356964374955379;
    // species 10: H2O
    result += y[10] *
              (+2.03399249e+00 + 2.17691804e-03 * T - 1.64072518e-07 * T2 -
               9.70419870e-11 * T3 + 1.68200992e-14 * T4) *
              0.0555092978073827;
    // species 11: CO2
    result += y[11] *
              (+2.85746029e+00 + 4.41437026e-03 * T - 2.21481404e-06 * T2 +
               5.23490188e-10 * T3 - 4.72084164e-14 * T4) *
              0.0227226249176305;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 4: NC7H16
    result += y[4] *
              (-2.26836187e+00 + 8.54355820e-02 * T - 5.25346786e-05 * T2 +
               1.62945721e-08 * T3 - 2.02394925e-12 * T4) *
              0.0099795419390250;
    // species 5: NC10H22
    result += y[5] *
              (-3.08416969e+00 + 1.22535012e-01 * T - 7.76815739e-05 * T2 +
               2.49834877e-08 * T3 - 3.23548038e-12 * T4) *
              0.0070280983371519;
    // species 6: NC12H26
    result += y[6] *
              (-3.62181594e+00 + 1.47237711e-01 * T - 9.43970271e-05 * T2 +
               3.07441268e-08 * T3 - 4.03602230e-12 * T4) *
              0.0058706117177410;
  } else {
    // species 4: NC7H16
    result += y[4] *
              (+2.12148969e+01 + 3.47675750e-02 * T - 1.18407129e-05 * T2 +
               1.83298478e-09 * T3 - 1.06130266e-13 * T4) *
              0.0099795419390250;
    // species 5: NC10H22
    result += y[5] *
              (+3.09882239e+01 + 4.77244922e-02 * T - 1.62276391e-05 * T2 +
               2.50963259e-09 * T3 - 1.45215772e-13 * T4) *
              0.0070280983371519;
    // species 6: NC12H26
    result += y[6] *
              (+3.75095037e+01 + 5.63550048e-02 * T - 1.91493200e-05 * T2 +
               2.96024862e-09 * T3 - 1.71244150e-13 * T4) *
              0.0058706117177410;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 7: IC8H18
    result += y[7] *
              (-5.20868893e+00 + 1.11440581e-01 * T - 7.91346582e-05 * T2 +
               2.92406242e-08 * T3 - 4.43743191e-12 * T4) *
              0.0087541144337839;
  } else {
    // species 7: IC8H18
    result += y[7] *
              (+2.61373590e+01 + 3.79004890e-02 * T - 1.29437358e-05 * T2 +
               2.00760372e-09 * T3 - 1.16400580e-13 * T4) *
              0.0087541144337839;
  }

  cvbs = result * 8.31446261815324e+07;
}

// Returns the mean enthalpy of the mixture in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBML(const amrex::Real T, const amrex::Real x[], amrex::Real& hbml)
{
  amrex::Real result = 0;
  amrex::Real hml[12];                       //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesEnthalpy(hml, T);

  // perform dot product
  for (int id = 0; id < 12; ++id) {
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10264
    result +=
      y[0] *
      (+4.70491270e+00 + 2.94556635e-02 * T + 3.50000467e-05 * T2 -
       4.62527200e-08 * T3 + 1.64476862e-11 * T4 - 3.73918630e+04 * invT) *
      0.0064810074077915;
    // species 1: POSF10289
    result +=
      y[1] *
      (+4.02646060e-02 + 4.43746710e-02 * T + 1.42414260e-05 * T2 -
       2.92908850e-08 * T3 + 1.07993808e-11 * T4 - 3.51905350e+04 * invT) *
      0.0059767147194530;
    // species 2: POSF10325
    result +=
      y[2] *
      (+2.57855180e+00 + 3.63124160e-02 * T + 2.44121233e-05 * T2 -
       3.77886100e-08 * T3 + 1.37871158e-11 * T4 - 3.80298550e+04 * invT) *
      0.0064810074077915;
    // species 3: POSF11498
    result +=
      y[3] *
      (+4.21777530e+00 + 4.65421005e-02 * T + 2.69035810e-05 * T2 -
       4.36327100e-08 * T3 + 1.59941018e-11 * T4 - 4.70178360e+04 * invT) *
      0.0054239641584448;
    // species 8: O2
    result +=
      y[8] *
      (+3.78245636e+00 - 1.49836708e-03 * T + 3.28243400e-06 * T2 -
       2.42032377e-09 * T3 + 6.48745674e-13 * T4 - 1.06394356e+03 * invT) *
      0.0312519532470779;
    // species 9: N2
    result +=
      y[9] *
      (+3.29867700e+00 + 7.04120200e-04 * T - 1.32107400e-06 * T2 +
       1.41037875e-09 * T3 - 4.88970800e-13 * T4 - 1.02089990e+03 * invT) *
      0.0356964374955379;
    // species 10: H2O
    result +=
      y[10] *
      (+4.19864056e+00 - 1.01821705e-03 * T + 2.17346737e-06 * T2 -
       1.37199266e-09 * T3 + 3.54395634e-13 * T4 - 3.02937267e+04 * invT) *
      0.0555092978073827;
    // species 11: CO2
    result +=
      y[11] *
      (+2.35677352e+00 + 4.49229839e-03 * T - 2.37452090e-06 * T2 +
       6.14797555e-10 * T3 - 2.87399096e-14 * T4 - 4.83719697e+04 * invT) *
      0.0227226249176305;
  } else {
    // species 0: POSF10264
    result +=
      y[0] *
      (+2.58974230e+01 + 2.77310460e-02 * T - 5.77924600e-06 * T2 +
       4.39561300e-10 * T3 + 1.27943798e-14 * T4 - 4.63378050e+04 * invT) *
      0.0064810074077915;
    // species 1: POSF10289
    result +=
      y[1] *
      (+2.01150530e+01 + 3.45909075e-02 * T - 7.72399900e-06 * T2 +
       4.87994525e-10 * T3 + 5.40837280e-14 * T4 - 4.25496840e+04 * invT) *
      0.0059767147194530;
    // species 2: POSF10325
    result +=
      y[2] *
      (+2.13578910e+01 + 3.32945435e-02 * T - 8.65582300e-06 * T2 +
       1.11937502e-09 * T3 - 4.72143040e-14 * T4 - 4.57969960e+04 * invT) *
      0.0064810074077915;
    // species 3: POSF11498
    result +=
      y[3] *
      (+3.36430400e+01 + 3.39640455e-02 * T - 6.84866667e-06 * T2 +
       6.24239100e-10 * T3 - 1.39264102e-14 * T4 - 5.83616720e+04 * invT) *
      0.0054239641584448;
    // species 8: O2
    result +=
      y[8] *
      (+3.28253784e+00 + 7.41543770e-04 * T - 2.52655556e-07 * T2 +
       5.23676387e-11 * T3 - 4.33435588e-15 * T4 - 1.08845772e+03 * invT) *
      0.0312519532470779;
    // species 9: N2
    result +=
      y[9] *
      (+2.92664000e+00 + 7.43988400e-04 * T - 1.89492000e-07 * T2 +
       2.52425950e-11 * T3 - 1.35067020e-15 * T4 - 9.22797700e+02 * invT) *
      0.0356964374955379;
    // species 10: H2O
    result +=
      y[10] *
      (+3.03399249e+00 + 1.08845902e-03 * T - 5.46908393e-08 * T2 -
       2.42604967e-11 * T3 + 3.36401984e-15 * T4 - 3.00042971e+04 * invT) *
      0.0555092978073827;
    // species 11: CO2
    result +=
      y[11] *
      (+3.85746029e+00 + 2.20718513e-03 * T - 7.38271347e-07 * T2 +
       1.30872547e-10 * T3 - 9.44168328e-15 * T4 - 4.87591660e+04 * invT) *
      0.0227226249176305;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 4: NC7H16
    result +=
      y[4] *
      (-1.26836187e+00 + 4.27177910e-02 * T - 1.75115595e-05 * T2 +
       4.07364302e-09 * T3 - 4.04789850e-13 * T4 - 2.56586565e+04 * invT) *
      0.0099795419390250;
    // species 5: NC10H22
    result +=
      y[5] *
      (-2.08416969e+00 + 6.12675060e-02 * T - 2.58938580e-05 * T2 +
       6.24587193e-09 * T3 - 6.47096076e-13 * T4 - 3.43021863e+04 * invT) *
      0.0070280983371519;
    // species 6: NC12H26
    result +=
      y[6] *
      (-2.62181594e+00 + 7.36188555e-02 * T - 3.14656757e-05 * T2 +
       7.68603170e-09 * T3 - 8.07204460e-13 * T4 - 4.00654253e+04 * invT) *
      0.0058706117177410;
  } else {
    // species 4: NC7H16
    result +=
      y[4] *
      (+2.22148969e+01 + 1.73837875e-02 * T - 3.94690430e-06 * T2 +
       4.58246195e-10 * T3 - 2.12260532e-14 * T4 - 3.42760081e+04 * invT) *
      0.0099795419390250;
    // species 5: NC10H22
    result +=
      y[5] *
      (+3.19882239e+01 + 2.38622461e-02 * T - 5.40921303e-06 * T2 +
       6.27408147e-10 * T3 - 2.90431544e-14 * T4 - 4.66392840e+04 * invT) *
      0.0070280983371519;
    // species 6: NC12H26
    result +=
      y[6] *
      (+3.85095037e+01 + 2.81775024e-02 * T - 6.38310667e-06 * T2 +
       7.40062155e-10 * T3 - 3.42488300e-14 * T4 - 5.48843465e+04 * invT) *
      0.0058706117177410;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 7: IC8H18
    result +=
      y[7] *
      (-4.20868893e+00 + 5.57202905e-02 * T - 2.63782194e-05 * T2 +
       7.31015605e-09 * T3 - 8.87486382e-13 * T4 - 2.99446875e+04 * invT) *
      0.0087541144337839;
  } else {
    // species 7: IC8H18
    result +=
      y[7] *
      (+2.71373590e+01 + 1.89502445e-02 * T - 4.31457860e-06 * T2 +
       5.01900930e-10 * T3 - 2.32801160e-14 * T4 - 4.07958177e+04 * invT) *
      0.0087541144337839;
  }

  const amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  hbms = result * RT;
}

// get mean internal energy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBML(const amrex::Real T, const amrex::Real x[], amrex::Real& ubml)
{
  amrex::Real result = 0;
  amrex::Real uml[12];                       //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesInternalEnergy(uml, T);

  // perform dot product
  for (int id = 0; id < 12; ++id) {
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10264
    result +=
      y[0] *
      (+3.70491270e+00 + 2.94556635e-02 * T + 3.50000467e-05 * T2 -
       4.62527200e-08 * T3 + 1.64476862e-11 * T4 - 3.73918630e+04 * invT) *
      0.0064810074077915;
    // species 1: POSF10289
    result +=
      y[1] *
      (-9.59735394e-01 + 4.43746710e-02 * T + 1.42414260e-05 * T2 -
       2.92908850e-08 * T3 + 1.07993808e-11 * T4 - 3.51905350e+04 * invT) *
      0.0059767147194530;
    // species 2: POSF10325
    result +=
      y[2] *
      (+1.57855180e+00 + 3.63124160e-02 * T + 2.44121233e-05 * T2 -
       3.77886100e-08 * T3 + 1.37871158e-11 * T4 - 3.80298550e+04 * invT) *
      0.0064810074077915;
    // species 3: POSF11498
    result +=
      y[3] *
      (+3.21777530e+00 + 4.65421005e-02 * T + 2.69035810e-05 * T2 -
       4.36327100e-08 * T3 + 1.59941018e-11 * T4 - 4.70178360e+04 * invT) *
      0.0054239641584448;
    // species 8: O2
    result +=
      y[8] *
      (+2.78245636e+00 - 1.49836708e-03 * T + 3.28243400e-06 * T2 -
       2.42032377e-09 * T3 + 6.48745674e-13 * T4 - 1.06394356e+03 * invT) *
      0.0312519532470779;
    // species 9: N2
    result +=
      y[9] *
      (+2.29867700e+00 + 7.04120200e-04 * T - 1.32107400e-06 * T2 +
       1.41037875e-09 * T3 - 4.88970800e-13 * T4 - 1.02089990e+03 * invT) *
      0.0356964374955379;
    // species 10: H2O
    result +=
      y[10] *
      (+3.19864056e+00 - 1.01821705e-03 * T + 2.17346737e-06 * T2 -
       1.37199266e-09 * T3 + 3.54395634e-13 * T4 - 3.02937267e+04 * invT) *
      0.0555092978073827;
    // species 11: CO2
    result +=
      y[11] *
      (+1.35677352e+00 + 4.49229839e-03 * T - 2.37452090e-06 * T2 +
       6.14797555e-10 * T3 - 2.87399096e-14 * T4 - 4.83719697e+04 * invT) *
      0.0227226249176305;
  } else {
    // species 0: POSF10264
    result +=
      y[0] *
      (+2.48974230e+01 + 2.77310460e-02 * T - 5.77924600e-06 * T2 +
       4.39561300e-10 * T3 + 1.27943798e-14 * T4 - 4.63378050e+04 * invT) *
      0.0064810074077915;
    // species 1: POSF10289
    result +=
      y[1] *
      (+1.91150530e+01 + 3.45909075e-02 * T - 7.72399900e-06 * T2 +
       4.87994525e-10 * T3 + 5.40837280e-14 * T4 - 4.25496840e+04 * invT) *
      0.0059767147194530;
    // species 2: POSF10325
    result +=
      y[2] *
      (+2.03578910e+01 + 3.32945435e-02 * T - 8.65582300e-06 * T2 +
       1.11937502e-09 * T3 - 4.72143040e-14 * T4 - 4.57969960e+04 * invT) *
      0.0064810074077915;
    // species 3: POSF11498
    result +=
      y[3] *
      (+3.26430400e+01 + 3.39640455e-02 * T - 6.84866667e-06 * T2 +
       6.24239100e-10 * T3 - 1.39264102e-14 * T4 - 5.83616720e+04 * invT) *
      0.0054239641584448;
    // species 8: O2
    result +=
      y[8] *
      (+2.28253784e+00 + 7.41543770e-04 * T - 2.52655556e-07 * T2 +
       5.23676387e-11 * T3 - 4.33435588e-15 * T4 - 1.08845772e+03 * invT) *
      0.0312519532470779;
    // species 9: N2
    result +=
      y[9] *
      (+1.92664000e+00 + 7.43988400e-04 * T - 1.89492000e-07 * T2 +
       2.52425950e-11 * T3 - 1.35067020e-15 * T4 - 9.22797700e+02 * invT) *
      0.0356964374955379;
    // species 10: H2O
    result +=
      y[10] *
      (+2.03399249e+00 + 1.08845902e-03 * T - 5.46908393e-08 * T2 -
       2.42604967e-11 * T3 + 3.36401984e-15 * T4 - 3.00042971e+04 * invT) *
      0.0555092978073827;
    // species 11: CO2
    result +=
      y[11] *
      (+2.85746029e+00 + 2.20718513e-03 * T - 7.38271347e-07 * T2 +
       1.30872547e-10 * T3 - 9.44168328e-15 * T4 - 4.87591660e+04 * invT) *
      0.0227226249176305;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 4: NC7H16
    result +=
      y[4] *
      (-2.26836187e+00 + 4.27177910e-02 * T - 1.75115595e-05 * T2 +
       4.07364302e-09 * T3 - 4.04789850e-13 * T4 - 2.56586565e+04 * invT) *
      0.0099795419390250;
    // species 5: NC10H22
    result +=
      y[5] *
      (-3.08416969e+00 + 6.12675060e-02 * T - 2.58938580e-05 * T2 +
       6.24587193e-09 * T3 - 6.47096076e-13 * T4 - 3.43021863e+04 * invT) *
      0.0070280983371519;
    // species 6: NC12H26
    result +=
      y[6] *
      (-3.62181594e+00 + 7.36188555e-02 * T - 3.14656757e-05 * T2 +
       7.68603170e-09 * T3 - 8.07204460e-13 * T4 - 4.00654253e+04 * invT) *
      0.0058706117177410;
  } else {
    // species 4: NC7H16
    result +=
      y[4] *
      (+2.12148969e+01 + 1.73837875e-02 * T - 3.94690430e-06 * T2 +
       4.58246195e-10 * T3 - 2.12260532e-14 * T4 - 3.42760081e+04 * invT) *
      0.0099795419390250;
    // species 5: NC10H22
    result +=
      y[5] *
      (+3.09882239e+01 + 2.38622461e-02 * T - 5.40921303e-06 * T2 +
       6.27408147e-10 * T3 - 2.90431544e-14 * T4 - 4.66392840e+04 * invT) *
      0.0070280983371519;
    // species 6: NC12H26
    result +=
      y[6] *
      (+3.75095037e+01 + 2.81775024e-02 * T - 6.38310667e-06 * T2 +
       7.40062155e-10 * T3 - 3.42488300e-14 * T4 - 5.48843465e+04 * invT) *
      0.0058706117177410;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 7: IC8H18
    result +=
      y[7] *
      (-5.20868893e+00 + 5.57202905e-02 * T - 2.63782194e-05 * T2 +
       7.31015605e-09 * T3 - 8.87486382e-13 * T4 - 2.99446875e+04 * invT) *
      0.0087541144337839;
  } else {
    // species 7: IC8H18
    result +=
      y[7] *
      (+2.61373590e+01 + 1.89502445e-02 * T - 4.31457860e-06 * T2 +
       5.01900930e-10 * T3 - 2.32801160e-14 * T4 - 4.07958177e+04 * invT) *
      0.0087541144337839;
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
  amrex::Real sor[12]; //  temporary storage
  speciesEntropy(sor, T);

  // Compute Eq 42
  for (int id = 0; id < 12; ++id) {
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
  amrex::Real sor[12]; //  temporary storage
  amrex::Real x[12];   //  need a ytx conversion
  amrex::Real YOW = 0; // See Eq 4, 6 in CK Manual

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 12; i++) {
    YOW += y[i] * imw(i);
  }
  // Now compute y to x conversion
  x[0] = y[0] / (1.5429699999999997e+02 * YOW);
  x[1] = y[1] / (1.6731600000000000e+02 * YOW);
  x[2] = y[2] / (1.5429699999999997e+02 * YOW);
  x[3] = y[3] / (1.8436699999999999e+02 * YOW);
  x[4] = y[4] / (1.0020500000000000e+02 * YOW);
  x[5] = y[5] / (1.4228600000000000e+02 * YOW);
  x[6] = y[6] / (1.7034000000000000e+02 * YOW);
  x[7] = y[7] / (1.1423200000000000e+02 * YOW);
  x[8] = y[8] / (3.1998000000000001e+01 * YOW);
  x[9] = y[9] / (2.8013999999999999e+01 * YOW);
  x[10] = y[10] / (1.8015000000000001e+01 * YOW);
  x[11] = y[11] / (4.4009000000000000e+01 * YOW);
  speciesEntropy(sor, T);
  // Perform computation in Eq 42 and 43
  for (int i = 0; i < 12; i++) {
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
  amrex::Real tmin = 300.0;  // max lower bound for thermo def
  amrex::Real tmax = 3000.0; // min upper bound for thermo def
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
  amrex::Real tmin = 300.0;  // max lower bound for thermo def
  amrex::Real tmax = 3000.0; // min upper bound for thermo def
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
  XW += x[0] * 1.5429699999999997e+02;     // POSF10264
  XW += x[1] * 1.6731600000000000e+02;     // POSF10289
  XW += x[2] * 1.5429699999999997e+02;     // POSF10325
  XW += x[3] * 1.8436699999999999e+02;     // POSF11498
  XW += x[4] * 1.0020500000000000e+02;     // NC7H16
  XW += x[5] * 1.4228600000000000e+02;     // NC10H22
  XW += x[6] * 1.7034000000000000e+02;     // NC12H26
  XW += x[7] * 1.1423200000000000e+02;     // IC8H18
  XW += x[8] * 3.1998000000000001e+01;     // O2
  XW += x[9] * 2.8013999999999999e+01;     // N2
  XW += x[10] * 1.8015000000000001e+01;    // H2O
  XW += x[11] * 4.4009000000000000e+01;    // CO2
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

  for (int i = 0; i < 12; i++) {
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
  W += c[0] * 1.5429699999999997e+02;  // POSF10264
  W += c[1] * 1.6731600000000000e+02;  // POSF10289
  W += c[2] * 1.5429699999999997e+02;  // POSF10325
  W += c[3] * 1.8436699999999999e+02;  // POSF11498
  W += c[4] * 1.0020500000000000e+02;  // NC7H16
  W += c[5] * 1.4228600000000000e+02;  // NC10H22
  W += c[6] * 1.7034000000000000e+02;  // NC12H26
  W += c[7] * 1.1423200000000000e+02;  // IC8H18
  W += c[8] * 3.1998000000000001e+01;  // O2
  W += c[9] * 2.8013999999999999e+01;  // N2
  W += c[10] * 1.8015000000000001e+01; // H2O
  W += c[11] * 4.4009000000000000e+01; // CO2

  for (int id = 0; id < 12; ++id) {
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
  XW += x[0] * 1.5429699999999997e+02;       // POSF10264
  XW += x[1] * 1.6731600000000000e+02;       // POSF10289
  XW += x[2] * 1.5429699999999997e+02;       // POSF10325
  XW += x[3] * 1.8436699999999999e+02;       // POSF11498
  XW += x[4] * 1.0020500000000000e+02;       // NC7H16
  XW += x[5] * 1.4228600000000000e+02;       // NC10H22
  XW += x[6] * 1.7034000000000000e+02;       // NC12H26
  XW += x[7] * 1.1423200000000000e+02;       // IC8H18
  XW += x[8] * 3.1998000000000001e+01;       // O2
  XW += x[9] * 2.8013999999999999e+01;       // N2
  XW += x[10] * 1.8015000000000001e+01;      // H2O
  XW += x[11] * 4.4009000000000000e+01;      // CO2
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

  for (int i = 0; i < 12; i++) {
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
  W += c[0] * 1.5429699999999997e+02;  // POSF10264
  W += c[1] * 1.6731600000000000e+02;  // POSF10289
  W += c[2] * 1.5429699999999997e+02;  // POSF10325
  W += c[3] * 1.8436699999999999e+02;  // POSF11498
  W += c[4] * 1.0020500000000000e+02;  // NC7H16
  W += c[5] * 1.4228600000000000e+02;  // NC10H22
  W += c[6] * 1.7034000000000000e+02;  // NC12H26
  W += c[7] * 1.1423200000000000e+02;  // IC8H18
  W += c[8] * 3.1998000000000001e+01;  // O2
  W += c[9] * 2.8013999999999999e+01;  // N2
  W += c[10] * 1.8015000000000001e+01; // H2O
  W += c[11] * 4.4009000000000000e+01; // CO2

  for (int id = 0; id < 12; ++id) {
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

  for (int i = 0; i < 12; i++) {
    YOW += y[i] * imw(i);
  }

  wtm = 1.0 / YOW;
}

// given x[species]: mole fractions
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWX(const amrex::Real x[], amrex::Real& wtm)
{
  amrex::Real XW = 0;                   //  see Eq 4 in CK Manual
  XW += x[0] * 1.5429699999999997e+02;  // POSF10264
  XW += x[1] * 1.6731600000000000e+02;  // POSF10289
  XW += x[2] * 1.5429699999999997e+02;  // POSF10325
  XW += x[3] * 1.8436699999999999e+02;  // POSF11498
  XW += x[4] * 1.0020500000000000e+02;  // NC7H16
  XW += x[5] * 1.4228600000000000e+02;  // NC10H22
  XW += x[6] * 1.7034000000000000e+02;  // NC12H26
  XW += x[7] * 1.1423200000000000e+02;  // IC8H18
  XW += x[8] * 3.1998000000000001e+01;  // O2
  XW += x[9] * 2.8013999999999999e+01;  // N2
  XW += x[10] * 1.8015000000000001e+01; // H2O
  XW += x[11] * 4.4009000000000000e+01; // CO2
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
  W += c[0] * 1.5429699999999997e+02;  // POSF10264
  W += c[1] * 1.6731600000000000e+02;  // POSF10289
  W += c[2] * 1.5429699999999997e+02;  // POSF10325
  W += c[3] * 1.8436699999999999e+02;  // POSF11498
  W += c[4] * 1.0020500000000000e+02;  // NC7H16
  W += c[5] * 1.4228600000000000e+02;  // NC10H22
  W += c[6] * 1.7034000000000000e+02;  // NC12H26
  W += c[7] * 1.1423200000000000e+02;  // IC8H18
  W += c[8] * 3.1998000000000001e+01;  // O2
  W += c[9] * 2.8013999999999999e+01;  // N2
  W += c[10] * 1.8015000000000001e+01; // H2O
  W += c[11] * 4.4009000000000000e+01; // CO2

  for (int id = 0; id < 12; ++id) {
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

  for (int i = 0; i < 12; i++) {
    YOW += y[i] * imw(i);
  }

  amrex::Real YOWINV = 1.0 / YOW;

  for (int i = 0; i < 12; i++) {
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
  for (int i = 0; i < 12; i++) {
    c[i] = y[i] * imw(i);
  }
  for (int i = 0; i < 12; i++) {
    YOW += c[i];
  }

  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // Now compute conversion

  for (int i = 0; i < 12; i++) {
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

  for (int i = 0; i < 12; i++) {
    c[i] = rho * y[i] * imw(i);
  }
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTY(const amrex::Real x[], amrex::Real y[])
{
  amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
  // Compute mean molecular wt first
  XW += x[0] * 1.5429699999999997e+02;  // POSF10264
  XW += x[1] * 1.6731600000000000e+02;  // POSF10289
  XW += x[2] * 1.5429699999999997e+02;  // POSF10325
  XW += x[3] * 1.8436699999999999e+02;  // POSF11498
  XW += x[4] * 1.0020500000000000e+02;  // NC7H16
  XW += x[5] * 1.4228600000000000e+02;  // NC10H22
  XW += x[6] * 1.7034000000000000e+02;  // NC12H26
  XW += x[7] * 1.1423200000000000e+02;  // IC8H18
  XW += x[8] * 3.1998000000000001e+01;  // O2
  XW += x[9] * 2.8013999999999999e+01;  // N2
  XW += x[10] * 1.8015000000000001e+01; // H2O
  XW += x[11] * 4.4009000000000000e+01; // CO2
  // Now compute conversion
  amrex::Real XWinv = 1.0 / XW;
  y[0] = x[0] * 1.5429699999999997e+02 * XWinv;
  y[1] = x[1] * 1.6731600000000000e+02 * XWinv;
  y[2] = x[2] * 1.5429699999999997e+02 * XWinv;
  y[3] = x[3] * 1.8436699999999999e+02 * XWinv;
  y[4] = x[4] * 1.0020500000000000e+02 * XWinv;
  y[5] = x[5] * 1.4228600000000000e+02 * XWinv;
  y[6] = x[6] * 1.7034000000000000e+02 * XWinv;
  y[7] = x[7] * 1.1423200000000000e+02 * XWinv;
  y[8] = x[8] * 3.1998000000000001e+01 * XWinv;
  y[9] = x[9] * 2.8013999999999999e+01 * XWinv;
  y[10] = x[10] * 1.8015000000000001e+01 * XWinv;
  y[11] = x[11] * 4.4009000000000000e+01 * XWinv;
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
  for (int id = 0; id < 12; ++id) {
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
  XW += x[0] * 1.5429699999999997e+02;  // POSF10264
  XW += x[1] * 1.6731600000000000e+02;  // POSF10289
  XW += x[2] * 1.5429699999999997e+02;  // POSF10325
  XW += x[3] * 1.8436699999999999e+02;  // POSF11498
  XW += x[4] * 1.0020500000000000e+02;  // NC7H16
  XW += x[5] * 1.4228600000000000e+02;  // NC10H22
  XW += x[6] * 1.7034000000000000e+02;  // NC12H26
  XW += x[7] * 1.1423200000000000e+02;  // IC8H18
  XW += x[8] * 3.1998000000000001e+01;  // O2
  XW += x[9] * 2.8013999999999999e+01;  // N2
  XW += x[10] * 1.8015000000000001e+01; // H2O
  XW += x[11] * 4.4009000000000000e+01; // CO2
  ROW = rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 12; ++id) {
    c[id] = x[id] * ROW;
  }
}

// convert c[species] (molar conc) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTX(const amrex::Real c[], amrex::Real x[])
{
  amrex::Real sumC = 0;

  // compute sum of c
  for (int id = 0; id < 12; ++id) {
    sumC += c[id];
  }

  //  See Eq 13
  amrex::Real sumCinv = 1.0 / sumC;
  for (int id = 0; id < 12; ++id) {
    x[id] = c[id] * sumCinv;
  }
}

// convert c[species] (molar conc) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTY(const amrex::Real c[], amrex::Real y[])
{
  amrex::Real CW = 0; // See Eq 12 in CK Manual
  // compute denominator in eq 12 first
  CW += c[0] * 1.5429699999999997e+02;  // POSF10264
  CW += c[1] * 1.6731600000000000e+02;  // POSF10289
  CW += c[2] * 1.5429699999999997e+02;  // POSF10325
  CW += c[3] * 1.8436699999999999e+02;  // POSF11498
  CW += c[4] * 1.0020500000000000e+02;  // NC7H16
  CW += c[5] * 1.4228600000000000e+02;  // NC10H22
  CW += c[6] * 1.7034000000000000e+02;  // NC12H26
  CW += c[7] * 1.1423200000000000e+02;  // IC8H18
  CW += c[8] * 3.1998000000000001e+01;  // O2
  CW += c[9] * 2.8013999999999999e+01;  // N2
  CW += c[10] * 1.8015000000000001e+01; // H2O
  CW += c[11] * 4.4009000000000000e+01; // CO2
  // Now compute conversion
  amrex::Real CWinv = 1.0 / CW;
  y[0] = c[0] * 1.5429699999999997e+02 * CWinv;
  y[1] = c[1] * 1.6731600000000000e+02 * CWinv;
  y[2] = c[2] * 1.5429699999999997e+02 * CWinv;
  y[3] = c[3] * 1.8436699999999999e+02 * CWinv;
  y[4] = c[4] * 1.0020500000000000e+02 * CWinv;
  y[5] = c[5] * 1.4228600000000000e+02 * CWinv;
  y[6] = c[6] * 1.7034000000000000e+02 * CWinv;
  y[7] = c[7] * 1.1423200000000000e+02 * CWinv;
  y[8] = c[8] * 3.1998000000000001e+01 * CWinv;
  y[9] = c[9] * 2.8013999999999999e+01 * CWinv;
  y[10] = c[10] * 1.8015000000000001e+01 * CWinv;
  y[11] = c[11] * 4.4009000000000000e+01 * CWinv;
}

// get specific heat at constant volume as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVML(const amrex::Real T, amrex::Real cvml[])
{
  cv_R(cvml, T);

  // convert to chemkin units
  for (int id = 0; id < 12; ++id) {
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
  for (int id = 0; id < 12; ++id) {
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
  for (int id = 0; id < 12; ++id) {
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
  for (int id = 0; id < 12; ++id) {
    hml[id] *= RT;
  }
}

// Returns the standard-state entropies in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSML(const amrex::Real T, amrex::Real sml[])
{
  speciesEntropy(sml, T);

  // convert to chemkin units
  for (int id = 0; id < 12; ++id) {
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
  cvms[0] *= 5.388609382005639e+05;  // POSF10264
  cvms[1] *= 4.969317111425829e+05;  // POSF10289
  cvms[2] *= 5.388609382005639e+05;  // POSF10325
  cvms[3] *= 4.509734723759263e+05;  // POSF11498
  cvms[4] *= 8.297452839831585e+05;  // NC7H16
  cvms[5] *= 5.843486090095469e+05;  // NC10H22
  cvms[6] *= 4.881098167284983e+05;  // NC12H26
  cvms[7] *= 7.278575721473177e+05;  // IC8H18
  cvms[8] *= 2.598431970171023e+06;  // O2
  cvms[9] *= 2.967966951578939e+06;  // N2
  cvms[10] *= 4.615299815794193e+06; // H2O
  cvms[11] *= 1.889264154639560e+06; // CO2
}

// Returns the specific heats at constant pressure
// in mass units (Eq. 26)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPMS(const amrex::Real T, amrex::Real cpms[])
{
  cp_R(cpms, T);
  // multiply by R/molecularweight
  cpms[0] *= 5.388609382005639e+05;  // POSF10264
  cpms[1] *= 4.969317111425829e+05;  // POSF10289
  cpms[2] *= 5.388609382005639e+05;  // POSF10325
  cpms[3] *= 4.509734723759263e+05;  // POSF11498
  cpms[4] *= 8.297452839831585e+05;  // NC7H16
  cpms[5] *= 5.843486090095469e+05;  // NC10H22
  cpms[6] *= 4.881098167284983e+05;  // NC12H26
  cpms[7] *= 7.278575721473177e+05;  // IC8H18
  cpms[8] *= 2.598431970171023e+06;  // O2
  cpms[9] *= 2.967966951578939e+06;  // N2
  cpms[10] *= 4.615299815794193e+06; // H2O
  cpms[11] *= 1.889264154639560e+06; // CO2
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUMS(const amrex::Real T, amrex::Real ums[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  speciesInternalEnergy(ums, T);

  for (int i = 0; i < 12; i++) {
    ums[i] *= RT * imw(i);
  }
}

// Returns enthalpy in mass units (Eq 27.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHMS(const amrex::Real T, amrex::Real hms[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  speciesEnthalpy(hms, T);

  for (int i = 0; i < 12; i++) {
    hms[i] *= RT * imw(i);
  }
}

// Returns the entropies in mass units (Eq 28.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSMS(const amrex::Real T, amrex::Real sms[])
{
  speciesEntropy(sms, T);
  // multiply by R/molecularweight
  sms[0] *= 5.388609382005639e+05;  // POSF10264
  sms[1] *= 4.969317111425829e+05;  // POSF10289
  sms[2] *= 5.388609382005639e+05;  // POSF10325
  sms[3] *= 4.509734723759263e+05;  // POSF11498
  sms[4] *= 8.297452839831585e+05;  // NC7H16
  sms[5] *= 5.843486090095469e+05;  // NC10H22
  sms[6] *= 4.881098167284983e+05;  // NC12H26
  sms[7] *= 7.278575721473177e+05;  // IC8H18
  sms[8] *= 2.598431970171023e+06;  // O2
  sms[9] *= 2.967966951578939e+06;  // N2
  sms[10] *= 4.615299815794193e+06; // H2O
  sms[11] *= 1.889264154639560e+06; // CO2
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

  for (int i = 0; i < 12; ++i) {
    wdot[i] = 0.0;
  }
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[])
{

  // convert to SI
  for (int id = 0; id < 12; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, T);

  // convert to chemkin units
  for (int id = 0; id < 12; ++id) {
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
  amrex::Real c[12]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 12; i++) {
    YOW += y[i] * imw(i);
  }
  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // multiply by 1e6 so c goes to SI
  PWORT *= 1e6;
  // Now compute conversion (and go to SI)
  for (int i = 0; i < 12; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 12; ++id) {
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
  amrex::Real c[12]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 12; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 12; ++id) {
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
  amrex::Real c[12]; // temporary storage

  // See Eq 8 with an extra 1e6 so c goes to SI
  for (int i = 0; i < 12; i++) {
    c[i] = 1e6 * rho * y[i] * imw(i);
  }

  // call productionRate
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 12; ++id) {
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
  amrex::Real c[12];  // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 1.5429699999999997e+02;  // POSF10264
  XW += x[1] * 1.6731600000000000e+02;  // POSF10289
  XW += x[2] * 1.5429699999999997e+02;  // POSF10325
  XW += x[3] * 1.8436699999999999e+02;  // POSF11498
  XW += x[4] * 1.0020500000000000e+02;  // NC7H16
  XW += x[5] * 1.4228600000000000e+02;  // NC10H22
  XW += x[6] * 1.7034000000000000e+02;  // NC12H26
  XW += x[7] * 1.1423200000000000e+02;  // IC8H18
  XW += x[8] * 3.1998000000000001e+01;  // O2
  XW += x[9] * 2.8013999999999999e+01;  // N2
  XW += x[10] * 1.8015000000000001e+01; // H2O
  XW += x[11] * 4.4009000000000000e+01; // CO2
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 12; ++id) {
    c[id] = x[id] * ROW;
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 12; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int kcharge[])
{
  kcharge[0] = 0;  // POSF10264
  kcharge[1] = 0;  // POSF10289
  kcharge[2] = 0;  // POSF10325
  kcharge[3] = 0;  // POSF11498
  kcharge[4] = 0;  // NC7H16
  kcharge[5] = 0;  // NC10H22
  kcharge[6] = 0;  // NC12H26
  kcharge[7] = 0;  // IC8H18
  kcharge[8] = 0;  // O2
  kcharge[9] = 0;  // N2
  kcharge[10] = 0; // H2O
  kcharge[11] = 0; // CO2
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real zk[])
{

  int kchrg[12];
  CKCHRG(kchrg);

  for (int id = 0; id < 12; ++id) {
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
  for (int i = 0; i < 169; i++) {
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
  amrex::Real c[12];

  for (int k = 0; k < 12; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian_precond(J, c, *Tp, *HP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 12; k++) {
    J[156 + k] *= 1.e-6;
    J[k * 13 + 12] *= 1.e6;
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
  for (int i = 0; i < 169; i++) {
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
  amrex::Real c[12];

  for (int k = 0; k < 12; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian(J, c, *Tp, *consP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 12; k++) {
    J[156 + k] *= 1.e-6;
    J[k * 13 + 12] *= 1.e6;
  }
}

// Transport function declarations

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENIMC(int* LENIMC)
{
  *LENIMC = 48;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENRMC(int* LENRMC)
{
  *LENRMC = 3060;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNO(int* NO)
{
  *NO = 4;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKK(int* KK)
{
  *KK = 12;
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
  WT[0] = 1.54297000E+02;
  WT[1] = 1.67316000E+02;
  WT[2] = 1.54297000E+02;
  WT[3] = 1.84367000E+02;
  WT[4] = 1.00205000E+02;
  WT[5] = 1.42286000E+02;
  WT[6] = 1.70340000E+02;
  WT[7] = 1.14232000E+02;
  WT[8] = 3.19980000E+01;
  WT[9] = 2.80140000E+01;
  WT[10] = 1.80150000E+01;
  WT[11] = 4.40090000E+01;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetEPS(amrex::Real* EPS)
{
  EPS[0] = 1.77645000E+02;
  EPS[1] = 1.78548000E+02;
  EPS[2] = 1.77645000E+02;
  EPS[3] = 1.79321000E+02;
  EPS[4] = 4.59600000E+02;
  EPS[5] = 7.04917000E+02;
  EPS[6] = 7.89980000E+02;
  EPS[7] = 4.58500000E+02;
  EPS[8] = 1.07400000E+02;
  EPS[9] = 9.75300000E+01;
  EPS[10] = 5.72400000E+02;
  EPS[11] = 2.44000000E+02;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetSIG(amrex::Real* SIG)
{
  SIG[0] = 9.60500000E+00;
  SIG[1] = 1.01150000E+01;
  SIG[2] = 9.60500000E+00;
  SIG[3] = 1.06070000E+01;
  SIG[4] = 6.25300000E+00;
  SIG[5] = 6.67500000E+00;
  SIG[6] = 7.04700000E+00;
  SIG[7] = 6.41400000E+00;
  SIG[8] = 3.45800000E+00;
  SIG[9] = 3.62100000E+00;
  SIG[10] = 2.60500000E+00;
  SIG[11] = 3.76300000E+00;
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
  DIP[5] = 0.00000000E+00;
  DIP[6] = 0.00000000E+00;
  DIP[7] = 0.00000000E+00;
  DIP[8] = 0.00000000E+00;
  DIP[9] = 0.00000000E+00;
  DIP[10] = 1.84400000E+00;
  DIP[11] = 0.00000000E+00;
}

// the polarizability in cubic Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPOL(amrex::Real* POL)
{
  POL[0] = 0.00000000E+00;
  POL[1] = 0.00000000E+00;
  POL[2] = 0.00000000E+00;
  POL[3] = 0.00000000E+00;
  POL[4] = 0.00000000E+00;
  POL[5] = 0.00000000E+00;
  POL[6] = 0.00000000E+00;
  POL[7] = 0.00000000E+00;
  POL[8] = 1.60000000E+00;
  POL[9] = 1.76000000E+00;
  POL[10] = 0.00000000E+00;
  POL[11] = 2.65000000E+00;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetZROT(amrex::Real* ZROT)
{
  ZROT[0] = 1.00000000E+00;
  ZROT[1] = 1.00000000E+00;
  ZROT[2] = 1.00000000E+00;
  ZROT[3] = 1.00000000E+00;
  ZROT[4] = 1.00000000E+00;
  ZROT[5] = 1.00000000E+00;
  ZROT[6] = 1.00000000E+00;
  ZROT[7] = 1.00000000E+00;
  ZROT[8] = 3.80000000E+00;
  ZROT[9] = 4.00000000E+00;
  ZROT[10] = 4.00000000E+00;
  ZROT[11] = 2.10000000E+00;
}

// 0: monoatomic, 1: linear, 2: nonlinear
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLIN(int* NLIN)
{
  NLIN[0] = 2;
  NLIN[1] = 2;
  NLIN[2] = 2;
  NLIN[3] = 2;
  NLIN[4] = 2;
  NLIN[5] = 2;
  NLIN[6] = 2;
  NLIN[7] = 2;
  NLIN[8] = 1;
  NLIN[9] = 1;
  NLIN[10] = 2;
  NLIN[11] = 1;
}

// Poly fits for the viscosities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFETA(amrex::Real* COFETA)
{
  COFETA[0] = -2.27773842E+01;
  COFETA[1] = 4.31076442E+00;
  COFETA[2] = -4.79725560E-01;
  COFETA[3] = 2.09551504E-02;
  COFETA[4] = -2.28855988E+01;
  COFETA[5] = 4.32699538E+00;
  COFETA[6] = -4.81698313E-01;
  COFETA[7] = 2.10350457E-02;
  COFETA[8] = -2.27773842E+01;
  COFETA[9] = 4.31076442E+00;
  COFETA[10] = -4.79725560E-01;
  COFETA[11] = 2.09551504E-02;
  COFETA[12] = -2.29704633E+01;
  COFETA[13] = 4.34078437E+00;
  COFETA[14] = -4.83377394E-01;
  COFETA[15] = 2.11032300E-02;
  COFETA[16] = -2.19841200E+01;
  COFETA[17] = 3.46341268E+00;
  COFETA[18] = -2.80516687E-01;
  COFETA[19] = 8.70427548E-03;
  COFETA[20] = -1.07789332E+01;
  COFETA[21] = -1.57700426E+00;
  COFETA[22] = 4.51648325E-01;
  COFETA[23] = -2.60205017E-02;
  COFETA[24] = -7.22654243E+00;
  COFETA[25] = -3.12842246E+00;
  COFETA[26] = 6.69477571E-01;
  COFETA[27] = -3.60312934E-02;
  COFETA[28] = -2.20110343E+01;
  COFETA[29] = 3.48324809E+00;
  COFETA[30] = -2.83512675E-01;
  COFETA[31] = 8.85089297E-03;
  COFETA[32] = -1.71618439E+01;
  COFETA[33] = 2.68036374E+00;
  COFETA[34] = -2.72570227E-01;
  COFETA[35] = 1.21650964E-02;
  COFETA[36] = -1.65695492E+01;
  COFETA[37] = 2.39056562E+00;
  COFETA[38] = -2.34558144E-01;
  COFETA[39] = 1.05024037E-02;
  COFETA[40] = -1.05420930E+01;
  COFETA[41] = -1.37777266E+00;
  COFETA[42] = 4.20502561E-01;
  COFETA[43] = -2.40627344E-02;
  COFETA[44] = -2.40015087E+01;
  COFETA[45] = 5.14359547E+00;
  COFETA[46] = -5.74269731E-01;
  COFETA[47] = 2.44937679E-02;
}

// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFLAM(amrex::Real* COFLAM)
{
  COFLAM[0] = -2.14027130E+01;
  COFLAM[1] = 9.10680973E+00;
  COFLAM[2] = -9.14963217E-01;
  COFLAM[3] = 3.28409140E-02;
  COFLAM[4] = -2.61362776E+01;
  COFLAM[5] = 1.10505805E+01;
  COFLAM[6] = -1.18630318E+00;
  COFLAM[7] = 4.54329501E-02;
  COFLAM[8] = -2.28849559E+01;
  COFLAM[9] = 9.72792992E+00;
  COFLAM[10] = -1.00067773E+00;
  COFLAM[11] = 3.67484895E-02;
  COFLAM[12] = -1.96463710E+01;
  COFLAM[13] = 8.38330580E+00;
  COFLAM[14] = -8.20313442E-01;
  COFLAM[15] = 2.88217869E-02;
  COFLAM[16] = -1.79582385E+01;
  COFLAM[17] = 7.27686839E+00;
  COFLAM[18] = -5.88898376E-01;
  COFLAM[19] = 1.49980246E-02;
  COFLAM[20] = -1.14916042E+01;
  COFLAM[21] = 4.27409243E+00;
  COFLAM[22] = -1.46234967E-01;
  COFLAM[23] = -6.16222086E-03;
  COFLAM[24] = -9.35546056E+00;
  COFLAM[25] = 3.31628813E+00;
  COFLAM[26] = -1.07145282E-02;
  COFLAM[27] = -1.24086090E-02;
  COFLAM[28] = -2.24747075E+01;
  COFLAM[29] = 9.23650708E+00;
  COFLAM[30] = -8.67708821E-01;
  COFLAM[31] = 2.80758651E-02;
  COFLAM[32] = -1.93717634E+00;
  COFLAM[33] = 2.89110257E+00;
  COFLAM[34] = -2.71096988E-01;
  COFLAM[35] = 1.15344939E-02;
  COFLAM[36] = 1.29305954E+01;
  COFLAM[37] = -3.52816965E+00;
  COFLAM[38] = 6.45498434E-01;
  COFLAM[39] = -3.19375023E-02;
  COFLAM[40] = 2.33729836E+01;
  COFLAM[41] = -8.96536178E+00;
  COFLAM[42] = 1.52828630E+00;
  COFLAM[43] = -7.58551616E-02;
  COFLAM[44] = -1.13649201E+01;
  COFLAM[45] = 5.88177354E+00;
  COFLAM[46] = -5.68651792E-01;
  COFLAM[47] = 2.03561484E-02;
}

// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFD(amrex::Real* COFD)
{
  COFD[0] = -2.18108044E+01;
  COFD[1] = 4.73185706E+00;
  COFD[2] = -3.97136602E-01;
  COFD[3] = 1.71287819E-02;
  COFD[4] = -2.19011342E+01;
  COFD[5] = 4.73819890E+00;
  COFD[6] = -3.97895603E-01;
  COFD[7] = 1.71591023E-02;
  COFD[8] = -2.19909049E+01;
  COFD[9] = 4.74452326E+00;
  COFD[10] = -3.98652154E-01;
  COFD[11] = 1.71893100E-02;
  COFD[12] = -2.18108044E+01;
  COFD[13] = 4.73185706E+00;
  COFD[14] = -3.97136602E-01;
  COFD[15] = 1.71287819E-02;
  COFD[16] = -2.19011342E+01;
  COFD[17] = 4.73819890E+00;
  COFD[18] = -3.97895603E-01;
  COFD[19] = 1.71591023E-02;
  COFD[20] = -2.18108044E+01;
  COFD[21] = 4.73185706E+00;
  COFD[22] = -3.97136602E-01;
  COFD[23] = 1.71287819E-02;
  COFD[24] = -2.19884881E+01;
  COFD[25] = 4.74358845E+00;
  COFD[26] = -3.98540349E-01;
  COFD[27] = 1.71848467E-02;
  COFD[28] = -2.20779057E+01;
  COFD[29] = 4.74987970E+00;
  COFD[30] = -3.99294082E-01;
  COFD[31] = 1.72150143E-02;
  COFD[32] = -2.19884881E+01;
  COFD[33] = 4.74358845E+00;
  COFD[34] = -3.98540349E-01;
  COFD[35] = 1.71848467E-02;
  COFD[36] = -2.21651184E+01;
  COFD[37] = 4.75530685E+00;
  COFD[38] = -3.99946130E-01;
  COFD[39] = 1.72411958E-02;
  COFD[40] = -2.40906670E+01;
  COFD[41] = 5.61211378E+00;
  COFD[42] = -4.90899453E-01;
  COFD[43] = 2.03798140E-02;
  COFD[44] = -2.41733837E+01;
  COFD[45] = 5.61195916E+00;
  COFD[46] = -4.90702336E-01;
  COFD[47] = 2.03642813E-02;
  COFD[48] = -2.40906670E+01;
  COFD[49] = 5.61211378E+00;
  COFD[50] = -4.90899453E-01;
  COFD[51] = 2.03798140E-02;
  COFD[52] = -2.42535569E+01;
  COFD[53] = 5.61184021E+00;
  COFD[54] = -4.90536060E-01;
  COFD[55] = 2.03511202E-02;
  COFD[56] = -2.17797069E+01;
  COFD[57] = 4.48837319E+00;
  COFD[58] = -2.95423315E-01;
  COFD[59] = 9.93861345E-03;
  COFD[60] = -2.41898976E+01;
  COFD[61] = 5.42163057E+00;
  COFD[62] = -4.47403738E-01;
  COFD[63] = 1.77447667E-02;
  COFD[64] = -2.42637999E+01;
  COFD[65] = 5.41668113E+00;
  COFD[66] = -4.46506166E-01;
  COFD[67] = 1.76957158E-02;
  COFD[68] = -2.41898976E+01;
  COFD[69] = 5.42163057E+00;
  COFD[70] = -4.47403738E-01;
  COFD[71] = 1.77447667E-02;
  COFD[72] = -2.43381182E+01;
  COFD[73] = 5.41271446E+00;
  COFD[74] = -4.45773700E-01;
  COFD[75] = 1.76552987E-02;
  COFD[76] = -1.90524160E+01;
  COFD[77] = 3.10412895E+00;
  COFD[78] = -8.65603859E-02;
  COFD[79] = -2.35997273E-04;
  COFD[80] = -1.51008484E+01;
  COFD[81] = 1.22557031E+00;
  COFD[82] = 1.86816623E-01;
  COFD[83] = -1.31748031E-02;
  COFD[84] = -2.40663867E+01;
  COFD[85] = 5.28427912E+00;
  COFD[86] = -4.23523371E-01;
  COFD[87] = 1.64688184E-02;
  COFD[88] = -2.41362914E+01;
  COFD[89] = 5.27734218E+00;
  COFD[90] = -4.22336184E-01;
  COFD[91] = 1.64058868E-02;
  COFD[92] = -2.40663867E+01;
  COFD[93] = 5.28427912E+00;
  COFD[94] = -4.23523371E-01;
  COFD[95] = 1.64688184E-02;
  COFD[96] = -2.42088391E+01;
  COFD[97] = 5.27218129E+00;
  COFD[98] = -4.21432600E-01;
  COFD[99] = 1.63574348E-02;
  COFD[100] = -1.80172717E+01;
  COFD[101] = 2.58457980E+00;
  COFD[102] = -1.02083759E-02;
  COFD[103] = -3.88694132E-03;
  COFD[104] = -1.39518811E+01;
  COFD[105] = 6.68848217E-01;
  COFD[106] = 2.66350983E-01;
  COFD[107] = -1.68784757E-02;
  COFD[108] = -1.28665148E+01;
  COFD[109] = 1.40407102E-01;
  COFD[110] = 3.41555146E-01;
  COFD[111] = -2.03644864E-02;
  COFD[112] = -2.41476812E+01;
  COFD[113] = 5.61216368E+00;
  COFD[114] = -4.90989108E-01;
  COFD[115] = 2.03869804E-02;
  COFD[116] = -2.42311225E+01;
  COFD[117] = 5.61204488E+00;
  COFD[118] = -4.90797373E-01;
  COFD[119] = 2.03717141E-02;
  COFD[120] = -2.41476812E+01;
  COFD[121] = 5.61216368E+00;
  COFD[122] = -4.90989108E-01;
  COFD[123] = 2.03869804E-02;
  COFD[124] = -2.43121018E+01;
  COFD[125] = 5.61187373E+00;
  COFD[126] = -4.90623569E-01;
  COFD[127] = 2.03581931E-02;
  COFD[128] = -2.18490678E+01;
  COFD[129] = 4.49449726E+00;
  COFD[130] = -2.96375694E-01;
  COFD[131] = 9.98604682E-03;
  COFD[132] = -1.91360670E+01;
  COFD[133] = 3.11407022E+00;
  COFD[134] = -8.80209213E-02;
  COFD[135] = -1.66273706E-04;
  COFD[136] = -1.81001252E+01;
  COFD[137] = 2.59333725E+00;
  COFD[138] = -1.15105084E-02;
  COFD[139] = -3.82426753E-03;
  COFD[140] = -2.19202028E+01;
  COFD[141] = 4.50060260E+00;
  COFD[142] = -2.97325595E-01;
  COFD[143] = 1.00333718E-02;
  COFD[144] = -1.87890298E+01;
  COFD[145] = 4.11991398E+00;
  COFD[146] = -3.22514659E-01;
  COFD[147] = 1.40877517E-02;
  COFD[148] = -1.88865073E+01;
  COFD[149] = 4.12490756E+00;
  COFD[150] = -3.23120293E-01;
  COFD[151] = 1.41123098E-02;
  COFD[152] = -1.87890298E+01;
  COFD[153] = 4.11991398E+00;
  COFD[154] = -3.22514659E-01;
  COFD[155] = 1.40877517E-02;
  COFD[156] = -1.89772705E+01;
  COFD[157] = 4.12916640E+00;
  COFD[158] = -3.23636780E-01;
  COFD[159] = 1.41332501E-02;
  COFD[160] = -2.13968193E+01;
  COFD[161] = 5.25183817E+00;
  COFD[162] = -4.57376333E-01;
  COFD[163] = 1.94504429E-02;
  COFD[164] = -2.26794657E+01;
  COFD[165] = 5.59897074E+00;
  COFD[166] = -4.91655801E-01;
  COFD[167] = 2.05094389E-02;
  COFD[168] = -2.28722333E+01;
  COFD[169] = 5.60324057E+00;
  COFD[170] = -4.88278590E-01;
  COFD[171] = 2.02076422E-02;
  COFD[172] = -2.14347202E+01;
  COFD[173] = 5.24832862E+00;
  COFD[174] = -4.56959708E-01;
  COFD[175] = 1.94339074E-02;
  COFD[176] = -1.55511233E+01;
  COFD[177] = 3.48070094E+00;
  COFD[178] = -2.40859499E-01;
  COFD[179] = 1.05972514E-02;
  COFD[180] = -1.84342672E+01;
  COFD[181] = 4.00051993E+00;
  COFD[182] = -3.07414279E-01;
  COFD[183] = 1.34475779E-02;
  COFD[184] = -1.85349237E+01;
  COFD[185] = 4.00765621E+00;
  COFD[186] = -3.08333282E-01;
  COFD[187] = 1.34872455E-02;
  COFD[188] = -1.84342672E+01;
  COFD[189] = 4.00051993E+00;
  COFD[190] = -3.07414279E-01;
  COFD[191] = 1.34475779E-02;
  COFD[192] = -1.86280811E+01;
  COFD[193] = 4.01370657E+00;
  COFD[194] = -3.09111839E-01;
  COFD[195] = 1.35208245E-02;
  COFD[196] = -2.10685665E+01;
  COFD[197] = 5.15027524E+00;
  COFD[198] = -4.46126111E-01;
  COFD[199] = 1.90401391E-02;
  COFD[200] = -2.24185383E+01;
  COFD[201] = 5.53218678E+00;
  COFD[202] = -4.85657384E-01;
  COFD[203] = 2.03539768E-02;
  COFD[204] = -2.27628265E+01;
  COFD[205] = 5.60192583E+00;
  COFD[206] = -4.91455403E-01;
  COFD[207] = 2.04776556E-02;
  COFD[208] = -2.11096851E+01;
  COFD[209] = 5.14898175E+00;
  COFD[210] = -4.46021078E-01;
  COFD[211] = 1.90381749E-02;
  COFD[212] = -1.52414498E+01;
  COFD[213] = 3.35922578E+00;
  COFD[214] = -2.25181399E-01;
  COFD[215] = 9.92132878E-03;
  COFD[216] = -1.49828551E+01;
  COFD[217] = 3.25781069E+00;
  COFD[218] = -2.12199367E-01;
  COFD[219] = 9.36657283E-03;
  COFD[220] = -2.30218519E+01;
  COFD[221] = 5.57736917E+00;
  COFD[222] = -4.77603538E-01;
  COFD[223] = 1.94525527E-02;
  COFD[224] = -2.31048740E+01;
  COFD[225] = 5.57414350E+00;
  COFD[226] = -4.76952003E-01;
  COFD[227] = 1.94149792E-02;
  COFD[228] = -2.30218519E+01;
  COFD[229] = 5.57736917E+00;
  COFD[230] = -4.77603538E-01;
  COFD[231] = 1.94525527E-02;
  COFD[232] = -2.31821707E+01;
  COFD[233] = 5.57110765E+00;
  COFD[234] = -4.76354539E-01;
  COFD[235] = 1.93809021E-02;
  COFD[236] = -1.91334458E+01;
  COFD[237] = 3.82263611E+00;
  COFD[238] = -1.93983472E-01;
  COFD[239] = 4.95789388E-03;
  COFD[240] = -1.57430646E+01;
  COFD[241] = 2.18116680E+00;
  COFD[242] = 4.90755572E-02;
  COFD[243] = -6.71158203E-03;
  COFD[244] = -1.47122747E+01;
  COFD[245] = 1.67339458E+00;
  COFD[246] = 1.22536705E-01;
  COFD[247] = -1.01692680E-02;
  COFD[248] = -1.91915956E+01;
  COFD[249] = 3.82896318E+00;
  COFD[250] = -1.94955707E-01;
  COFD[251] = 5.00579392E-03;
  COFD[252] = -2.12652444E+01;
  COFD[253] = 5.59961822E+00;
  COFD[254] = -4.91624856E-01;
  COFD[255] = 2.05035546E-02;
  COFD[256] = -2.10643264E+01;
  COFD[257] = 5.53614864E+00;
  COFD[258] = -4.86046753E-01;
  COFD[259] = 2.03659193E-02;
  COFD[260] = -1.19157788E+01;
  COFD[261] = 9.28952533E-01;
  COFD[262] = 2.42107483E-01;
  COFD[263] = -1.59824150E-02;
  COFD[264] = -2.18171057E+01;
  COFD[265] = 5.12625274E+00;
  COFD[266] = -4.43796461E-01;
  COFD[267] = 1.89713627E-02;
  COFD[268] = -2.19142145E+01;
  COFD[269] = 5.13044729E+00;
  COFD[270] = -4.44231042E-01;
  COFD[271] = 1.89857477E-02;
  COFD[272] = -2.18171057E+01;
  COFD[273] = 5.12625274E+00;
  COFD[274] = -4.43796461E-01;
  COFD[275] = 1.89713627E-02;
  COFD[276] = -2.20051534E+01;
  COFD[277] = 5.13402627E+00;
  COFD[278] = -4.44601346E-01;
  COFD[279] = 1.89979689E-02;
  COFD[280] = -2.28671163E+01;
  COFD[281] = 5.50522401E+00;
  COFD[282] = -4.63604304E-01;
  COFD[283] = 1.86600785E-02;
  COFD[284] = -2.21024981E+01;
  COFD[285] = 4.94882195E+00;
  COFD[286] = -3.68164829E-01;
  COFD[287] = 1.36011930E-02;
  COFD[288] = -2.17131539E+01;
  COFD[289] = 4.69815223E+00;
  COFD[290] = -3.28454103E-01;
  COFD[291] = 1.15970636E-02;
  COFD[292] = -2.29181126E+01;
  COFD[293] = 5.50621930E+00;
  COFD[294] = -4.63843664E-01;
  COFD[295] = 1.86748763E-02;
  COFD[296] = -1.84688303E+01;
  COFD[297] = 4.49330851E+00;
  COFD[298] = -3.68208715E-01;
  COFD[299] = 1.59565402E-02;
  COFD[300] = -1.81432499E+01;
  COFD[301] = 4.37565431E+00;
  COFD[302] = -3.53906025E-01;
  COFD[303] = 1.53760786E-02;
  COFD[304] = -2.07653629E+01;
  COFD[305] = 5.01091988E+00;
  COFD[306] = -3.77985581E-01;
  COFD[307] = 1.40968617E-02;
  COFD[308] = -2.13425604E+01;
  COFD[309] = 5.40460130E+00;
  COFD[310] = -4.72718910E-01;
  COFD[311] = 1.99362717E-02;
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

  amrex::Real EPS[12];
  amrex::Real SIG[12];
  amrex::Real wt[12];
  amrex::Real Rcst = 83.144598; // in bar [CGS] !
  amrex::Real avogadro = 6.02214199e23;
  amrex::Real boltzmann = 1.3806503e-16; // we work in CGS

  egtransetEPS(EPS);
  egtransetSIG(SIG);
  get_mw(wt);

  // species 0: POSF10264
  Tci[0] = 1.316 * EPS[0];
  ai[0] = (5.55 * avogadro * avogadro * EPS[0] * boltzmann * 1e-24 * SIG[0] *
           SIG[0] * SIG[0]) /
          (wt[0] * wt[0]);
  bi[0] = 0.855 * avogadro * 1e-24 * SIG[0] * SIG[0] * SIG[0] / (wt[0]);
  acentric_i[0] = 0.0;

  // species 1: POSF10289
  Tci[1] = 1.316 * EPS[1];
  ai[1] = (5.55 * avogadro * avogadro * EPS[1] * boltzmann * 1e-24 * SIG[1] *
           SIG[1] * SIG[1]) /
          (wt[1] * wt[1]);
  bi[1] = 0.855 * avogadro * 1e-24 * SIG[1] * SIG[1] * SIG[1] / (wt[1]);
  acentric_i[1] = 0.0;

  // species 2: POSF10325
  Tci[2] = 1.316 * EPS[2];
  ai[2] = (5.55 * avogadro * avogadro * EPS[2] * boltzmann * 1e-24 * SIG[2] *
           SIG[2] * SIG[2]) /
          (wt[2] * wt[2]);
  bi[2] = 0.855 * avogadro * 1e-24 * SIG[2] * SIG[2] * SIG[2] / (wt[2]);
  acentric_i[2] = 0.0;

  // species 3: POSF11498
  Tci[3] = 1.316 * EPS[3];
  ai[3] = (5.55 * avogadro * avogadro * EPS[3] * boltzmann * 1e-24 * SIG[3] *
           SIG[3] * SIG[3]) /
          (wt[3] * wt[3]);
  bi[3] = 0.855 * avogadro * 1e-24 * SIG[3] * SIG[3] * SIG[3] / (wt[3]);
  acentric_i[3] = 0.0;

  // species 4: NC7H16
  Tci[4] = 1.316 * EPS[4];
  ai[4] = (5.55 * avogadro * avogadro * EPS[4] * boltzmann * 1e-24 * SIG[4] *
           SIG[4] * SIG[4]) /
          (wt[4] * wt[4]);
  bi[4] = 0.855 * avogadro * 1e-24 * SIG[4] * SIG[4] * SIG[4] / (wt[4]);
  acentric_i[4] = 0.0;

  // species 5: NC10H22
  Tci[5] = 1.316 * EPS[5];
  ai[5] = (5.55 * avogadro * avogadro * EPS[5] * boltzmann * 1e-24 * SIG[5] *
           SIG[5] * SIG[5]) /
          (wt[5] * wt[5]);
  bi[5] = 0.855 * avogadro * 1e-24 * SIG[5] * SIG[5] * SIG[5] / (wt[5]);
  acentric_i[5] = 0.0;

  // species 6: NC12H26
  Tci[6] = 1.316 * EPS[6];
  ai[6] = (5.55 * avogadro * avogadro * EPS[6] * boltzmann * 1e-24 * SIG[6] *
           SIG[6] * SIG[6]) /
          (wt[6] * wt[6]);
  bi[6] = 0.855 * avogadro * 1e-24 * SIG[6] * SIG[6] * SIG[6] / (wt[6]);
  acentric_i[6] = 0.0;

  // species 7: IC8H18
  Tci[7] = 1.316 * EPS[7];
  ai[7] = (5.55 * avogadro * avogadro * EPS[7] * boltzmann * 1e-24 * SIG[7] *
           SIG[7] * SIG[7]) /
          (wt[7] * wt[7]);
  bi[7] = 0.855 * avogadro * 1e-24 * SIG[7] * SIG[7] * SIG[7] / (wt[7]);
  acentric_i[7] = 0.0;

  // species 8: O2
  // Imported from NIST
  Tci[8] = 154.581000;
  ai[8] = 1e6 * 0.42748 * Rcst * Rcst * Tci[8] * Tci[8] /
          (31.998800 * 31.998800 * 50.430466);
  bi[8] = 0.08664 * Rcst * Tci[8] / (31.998800 * 50.430466);
  acentric_i[8] = 0.022200;

  // species 9: N2
  // Imported from NIST
  Tci[9] = 126.192000;
  ai[9] = 1e6 * 0.42748 * Rcst * Rcst * Tci[9] * Tci[9] /
          (28.013400 * 28.013400 * 33.958000);
  bi[9] = 0.08664 * Rcst * Tci[9] / (28.013400 * 33.958000);
  acentric_i[9] = 0.037200;

  // species 10: H2O
  // Imported from NIST
  Tci[10] = 647.096000;
  ai[10] = 1e6 * 0.42748 * Rcst * Rcst * Tci[10] * Tci[10] /
           (18.015340 * 18.015340 * 220.640000);
  bi[10] = 0.08664 * Rcst * Tci[10] / (18.015340 * 220.640000);
  acentric_i[10] = 0.344300;

  // species 11: CO2
  // Imported from NIST
  Tci[11] = 304.120000;
  ai[11] = 1e6 * 0.42748 * Rcst * Rcst * Tci[11] * Tci[11] /
           (44.009950 * 44.009950 * 73.740000);
  bi[11] = 0.08664 * Rcst * Tci[11] / (44.009950 * 73.740000);
  acentric_i[11] = 0.225000;
}

// compute the critical parameter quantities for each species for SRK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS_SRK(
  amrex::Real* sqrtOneOverTc,
  amrex::Real* sqrtAsti,
  amrex::Real* Bi,
  amrex::Real* Fomega)
{

  // species 0: POSF10264
  // Computed from Lennard-Jones
  sqrtOneOverTc[0] = 6.5402682478301e-02;
  sqrtAsti[0] = 4.2865154585193e+04;
  Bi[0] = 2.9570027693584e+00;
  Fomega[0] = 4.8508000000000e-01;

  // species 1: POSF10289
  // Computed from Lennard-Jones
  sqrtOneOverTc[1] = 6.5237086998630e-02;
  sqrtAsti[1] = 4.2828038124846e+04;
  Bi[1] = 3.1847642808904e+00;
  Fomega[1] = 4.8508000000000e-01;

  // species 2: POSF10325
  // Computed from Lennard-Jones
  sqrtOneOverTc[2] = 6.5402682478301e-02;
  sqrtAsti[2] = 4.2865154585193e+04;
  Bi[2] = 2.9570027693584e+00;
  Fomega[2] = 4.8508000000000e-01;

  // species 3: POSF11498
  // Computed from Lennard-Jones
  sqrtOneOverTc[3] = 6.5096326209175e-02;
  sqrtAsti[3] = 4.1827371033283e+04;
  Bi[3] = 3.3328181938405e+00;
  Fomega[3] = 4.8508000000000e-01;

  // species 4: NC7H16
  // Computed from Lennard-Jones
  sqrtOneOverTc[4] = 4.0661373572610e-02;
  sqrtAsti[4] = 5.5766397379575e+04;
  Bi[4] = 1.2562989577857e+00;
  Fomega[4] = 4.8508000000000e-01;

  // species 5: NC10H22
  // Computed from Lennard-Jones
  sqrtOneOverTc[5] = 3.2832432738125e-02;
  sqrtAsti[5] = 5.3644234851217e+04;
  Bi[5] = 1.0762390490260e+00;
  Fomega[5] = 4.8508000000000e-01;

  // species 6: NC12H26
  // Computed from Lennard-Jones
  sqrtOneOverTc[6] = 3.1014444518975e-02;
  sqrtAsti[6] = 5.1456132490395e+04;
  Bi[6] = 1.0578236621500e+00;
  Fomega[6] = 4.8508000000000e-01;

  // species 7: IC8H18
  // Computed from Lennard-Jones
  sqrtOneOverTc[7] = 4.0710120264334e-02;
  sqrtAsti[7] = 5.0759186273762e+04;
  Bi[7] = 1.1893678422487e+00;
  Fomega[7] = 4.8508000000000e-01;

  // species 8: O2
  // Imported from NIST
  sqrtOneOverTc[8] = 8.0430717653241e-02;
  sqrtAsti[8] = 3.6980081924740e+04;
  Bi[8] = 6.9005216444999e-01;
  Fomega[8] = 5.1945301904908e-01;

  // species 9: N2
  // Imported from NIST
  sqrtOneOverTc[9] = 8.9019282240563e-02;
  sqrtAsti[9] = 4.2022983929562e+04;
  Bi[9] = 9.5560052221563e-01;
  Fomega[9] = 5.4259343186608e-01;

  // species 10: H2O
  // Imported from NIST
  sqrtOneOverTc[10] = 3.9311140369715e-02;
  sqrtAsti[10] = 1.3145468419592e+05;
  Bi[10] = 1.1727204136223e+00;
  Fomega[10] = 1.0013577274636e+00;

  // species 11: CO2
  // Imported from NIST
  sqrtOneOverTc[11] = 5.7342616962522e-02;
  sqrtAsti[11] = 4.3745610142919e+04;
  Bi[11] = 6.7506167841745e-01;
  Fomega[11] = 8.2653709187500e-01;
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
