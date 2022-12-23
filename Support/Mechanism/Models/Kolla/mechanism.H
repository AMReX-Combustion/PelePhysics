#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  C
1  H
2  O
3  N
*/

// Species
static constexpr int H2_ID = 0;
static constexpr int O2_ID = 1;
static constexpr int O_ID = 2;
static constexpr int OH_ID = 3;
static constexpr int H2O_ID = 4;
static constexpr int H_ID = 5;
static constexpr int HO2_ID = 6;
static constexpr int H2O2_ID = 7;
static constexpr int CO_ID = 8;
static constexpr int CO2_ID = 9;
static constexpr int HCO_ID = 10;
static constexpr int N2_ID = 11;

static constexpr int NUM_ELEMENTS = 4;
static constexpr int NUM_SPECIES = 12;
static int constexpr NUM_IONS = 0;
static constexpr int NUM_REACTIONS = 29;

static constexpr int NUM_FIT = 4;

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
  *mm = 4;
  *kk = 12;
  *ii = 29;
  *nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_imw(amrex::Real* imw_new)
{
  imw_new[0] = 1.0 / 2.016000;   // H2
  imw_new[1] = 1.0 / 31.998000;  // O2
  imw_new[2] = 1.0 / 15.999000;  // O
  imw_new[3] = 1.0 / 17.007000;  // OH
  imw_new[4] = 1.0 / 18.015000;  // H2O
  imw_new[5] = 1.0 / 1.008000;   // H
  imw_new[6] = 1.0 / 33.006000;  // HO2
  imw_new[7] = 1.0 / 34.014000;  // H2O2
  imw_new[8] = 1.0 / 28.010000;  // CO
  imw_new[9] = 1.0 / 44.009000;  // CO2
  imw_new[10] = 1.0 / 29.018000; // HCO
  imw_new[11] = 1.0 / 28.014000; // N2
}

//  molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_mw(amrex::Real* mw_new)
{
  mw_new[0] = 2.016000;   // H2
  mw_new[1] = 31.998000;  // O2
  mw_new[2] = 15.999000;  // O
  mw_new[3] = 17.007000;  // OH
  mw_new[4] = 18.015000;  // H2O
  mw_new[5] = 1.008000;   // H
  mw_new[6] = 33.006000;  // HO2
  mw_new[7] = 34.014000;  // H2O2
  mw_new[8] = 28.010000;  // CO
  mw_new[9] = 44.009000;  // CO2
  mw_new[10] = 29.018000; // HCO
  mw_new[11] = 28.014000; // N2
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
    // species 0: H2
    species[0] = +2.29812400e+00 + 8.24944200e-04 * tc[1] -
                 8.14301500e-07 * tc[2] - 9.47543400e-11 * tc[3] +
                 4.13487200e-13 * tc[4];
    // species 1: O2
    species[1] = +2.21293600e+00 + 1.12748600e-03 * tc[1] -
                 5.75615000e-07 * tc[2] + 1.31387700e-09 * tc[3] -
                 8.76855400e-13 * tc[4];
    // species 2: O
    species[2] = +1.94642900e+00 - 1.63816600e-03 * tc[1] +
                 2.42103200e-06 * tc[2] - 1.60284300e-09 * tc[3] +
                 3.89069600e-13 * tc[4];
    // species 3: OH
    species[3] = +3.12530561e+00 - 3.22544939e-03 * tc[1] +
                 6.52764691e-06 * tc[2] - 5.79853643e-09 * tc[3] +
                 2.06237379e-12 * tc[4];
    // species 4: H2O
    species[4] = +2.38684200e+00 + 3.47498200e-03 * tc[1] -
                 6.35469600e-06 * tc[2] + 6.96858100e-09 * tc[3] -
                 2.50658800e-12 * tc[4];
    // species 6: HO2
    species[6] = +3.30179801e+00 - 4.74912051e-03 * tc[1] +
                 2.11582891e-05 * tc[2] - 2.42763894e-08 * tc[3] +
                 9.29225124e-12 * tc[4];
    // species 7: H2O2
    species[7] = +2.38875400e+00 + 6.56922600e-03 * tc[1] -
                 1.48501300e-07 * tc[2] - 4.62580600e-09 * tc[3] +
                 2.47151500e-12 * tc[4];
    // species 8: CO
    species[8] = +2.26245200e+00 + 1.51194100e-03 * tc[1] -
                 3.88175500e-06 * tc[2] + 5.58194400e-09 * tc[3] -
                 2.47495100e-12 * tc[4];
    // species 9: CO2
    species[9] = +1.27572500e+00 + 9.92207200e-03 * tc[1] -
                 1.04091100e-05 * tc[2] + 6.86668700e-09 * tc[3] -
                 2.11728000e-12 * tc[4];
    // species 10: HCO
    species[10] = +1.89833000e+00 + 6.19914700e-03 * tc[1] -
                  9.62308400e-06 * tc[2] + 1.08982500e-08 * tc[3] -
                  4.57488500e-12 * tc[4];
    // species 11: N2
    species[11] = +2.29867700e+00 + 1.40824000e-03 * tc[1] -
                  3.96322200e-06 * tc[2] + 5.64151500e-09 * tc[3] -
                  2.44485500e-12 * tc[4];
  } else {
    // species 0: H2
    species[0] = +1.99142300e+00 + 7.00064400e-04 * tc[1] -
                 5.63382900e-08 * tc[2] - 9.23157800e-12 * tc[3] +
                 1.58275200e-15 * tc[4];
    // species 1: O2
    species[1] = +2.69757800e+00 + 6.13519700e-04 * tc[1] -
                 1.25884200e-07 * tc[2] + 1.77528100e-11 * tc[3] -
                 1.13643500e-15 * tc[4];
    // species 2: O
    species[2] = +1.54206000e+00 - 2.75506200e-05 * tc[1] -
                 3.10280300e-09 * tc[2] + 4.55106700e-12 * tc[3] -
                 4.36805200e-16 * tc[4];
    // species 3: OH
    species[3] = +1.86472886e+00 + 1.05650448e-03 * tc[1] -
                 2.59082758e-07 * tc[2] + 3.05218674e-11 * tc[3] -
                 1.33195876e-15 * tc[4];
    // species 4: H2O
    species[4] = +1.67214600e+00 + 3.05629300e-03 * tc[1] -
                 8.73026000e-07 * tc[2] + 1.20099600e-10 * tc[3] -
                 6.39161800e-15 * tc[4];
    // species 6: HO2
    species[6] = +3.01721090e+00 + 2.23982013e-03 * tc[1] -
                 6.33658150e-07 * tc[2] + 1.14246370e-10 * tc[3] -
                 1.07908535e-14 * tc[4];
    // species 7: H2O2
    species[7] = +3.57316700e+00 + 4.33613600e-03 * tc[1] -
                 1.47468900e-06 * tc[2] + 2.34890400e-10 * tc[3] -
                 1.43165400e-14 * tc[4];
    // species 8: CO
    species[8] = +2.02507800e+00 + 1.44268900e-03 * tc[1] -
                 5.63082800e-07 * tc[2] + 1.01858100e-10 * tc[3] -
                 6.91095200e-15 * tc[4];
    // species 9: CO2
    species[9] = +3.45362300e+00 + 3.14016900e-03 * tc[1] -
                 1.27841100e-06 * tc[2] + 2.39399700e-10 * tc[3] -
                 1.66903300e-14 * tc[4];
    // species 10: HCO
    species[10] = +2.55727100e+00 + 3.34557300e-03 * tc[1] -
                  1.33500600e-06 * tc[2] + 2.47057300e-10 * tc[3] -
                  1.71385100e-14 * tc[4];
    // species 11: N2
    species[11] = +1.92664000e+00 + 1.48797700e-03 * tc[1] -
                  5.68476100e-07 * tc[2] + 1.00970400e-10 * tc[3] -
                  6.75335100e-15 * tc[4];
  }

  // species with no change at a midpoint T
  // species 5: H
  species[5] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4];
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
    // species 0: H2
    species[0] = +3.29812400e+00 + 8.24944200e-04 * tc[1] -
                 8.14301500e-07 * tc[2] - 9.47543400e-11 * tc[3] +
                 4.13487200e-13 * tc[4];
    // species 1: O2
    species[1] = +3.21293600e+00 + 1.12748600e-03 * tc[1] -
                 5.75615000e-07 * tc[2] + 1.31387700e-09 * tc[3] -
                 8.76855400e-13 * tc[4];
    // species 2: O
    species[2] = +2.94642900e+00 - 1.63816600e-03 * tc[1] +
                 2.42103200e-06 * tc[2] - 1.60284300e-09 * tc[3] +
                 3.89069600e-13 * tc[4];
    // species 3: OH
    species[3] = +4.12530561e+00 - 3.22544939e-03 * tc[1] +
                 6.52764691e-06 * tc[2] - 5.79853643e-09 * tc[3] +
                 2.06237379e-12 * tc[4];
    // species 4: H2O
    species[4] = +3.38684200e+00 + 3.47498200e-03 * tc[1] -
                 6.35469600e-06 * tc[2] + 6.96858100e-09 * tc[3] -
                 2.50658800e-12 * tc[4];
    // species 6: HO2
    species[6] = +4.30179801e+00 - 4.74912051e-03 * tc[1] +
                 2.11582891e-05 * tc[2] - 2.42763894e-08 * tc[3] +
                 9.29225124e-12 * tc[4];
    // species 7: H2O2
    species[7] = +3.38875400e+00 + 6.56922600e-03 * tc[1] -
                 1.48501300e-07 * tc[2] - 4.62580600e-09 * tc[3] +
                 2.47151500e-12 * tc[4];
    // species 8: CO
    species[8] = +3.26245200e+00 + 1.51194100e-03 * tc[1] -
                 3.88175500e-06 * tc[2] + 5.58194400e-09 * tc[3] -
                 2.47495100e-12 * tc[4];
    // species 9: CO2
    species[9] = +2.27572500e+00 + 9.92207200e-03 * tc[1] -
                 1.04091100e-05 * tc[2] + 6.86668700e-09 * tc[3] -
                 2.11728000e-12 * tc[4];
    // species 10: HCO
    species[10] = +2.89833000e+00 + 6.19914700e-03 * tc[1] -
                  9.62308400e-06 * tc[2] + 1.08982500e-08 * tc[3] -
                  4.57488500e-12 * tc[4];
    // species 11: N2
    species[11] = +3.29867700e+00 + 1.40824000e-03 * tc[1] -
                  3.96322200e-06 * tc[2] + 5.64151500e-09 * tc[3] -
                  2.44485500e-12 * tc[4];
  } else {
    // species 0: H2
    species[0] = +2.99142300e+00 + 7.00064400e-04 * tc[1] -
                 5.63382900e-08 * tc[2] - 9.23157800e-12 * tc[3] +
                 1.58275200e-15 * tc[4];
    // species 1: O2
    species[1] = +3.69757800e+00 + 6.13519700e-04 * tc[1] -
                 1.25884200e-07 * tc[2] + 1.77528100e-11 * tc[3] -
                 1.13643500e-15 * tc[4];
    // species 2: O
    species[2] = +2.54206000e+00 - 2.75506200e-05 * tc[1] -
                 3.10280300e-09 * tc[2] + 4.55106700e-12 * tc[3] -
                 4.36805200e-16 * tc[4];
    // species 3: OH
    species[3] = +2.86472886e+00 + 1.05650448e-03 * tc[1] -
                 2.59082758e-07 * tc[2] + 3.05218674e-11 * tc[3] -
                 1.33195876e-15 * tc[4];
    // species 4: H2O
    species[4] = +2.67214600e+00 + 3.05629300e-03 * tc[1] -
                 8.73026000e-07 * tc[2] + 1.20099600e-10 * tc[3] -
                 6.39161800e-15 * tc[4];
    // species 6: HO2
    species[6] = +4.01721090e+00 + 2.23982013e-03 * tc[1] -
                 6.33658150e-07 * tc[2] + 1.14246370e-10 * tc[3] -
                 1.07908535e-14 * tc[4];
    // species 7: H2O2
    species[7] = +4.57316700e+00 + 4.33613600e-03 * tc[1] -
                 1.47468900e-06 * tc[2] + 2.34890400e-10 * tc[3] -
                 1.43165400e-14 * tc[4];
    // species 8: CO
    species[8] = +3.02507800e+00 + 1.44268900e-03 * tc[1] -
                 5.63082800e-07 * tc[2] + 1.01858100e-10 * tc[3] -
                 6.91095200e-15 * tc[4];
    // species 9: CO2
    species[9] = +4.45362300e+00 + 3.14016900e-03 * tc[1] -
                 1.27841100e-06 * tc[2] + 2.39399700e-10 * tc[3] -
                 1.66903300e-14 * tc[4];
    // species 10: HCO
    species[10] = +3.55727100e+00 + 3.34557300e-03 * tc[1] -
                  1.33500600e-06 * tc[2] + 2.47057300e-10 * tc[3] -
                  1.71385100e-14 * tc[4];
    // species 11: N2
    species[11] = +2.92664000e+00 + 1.48797700e-03 * tc[1] -
                  5.68476100e-07 * tc[2] + 1.00970400e-10 * tc[3] -
                  6.75335100e-15 * tc[4];
  }

  // species with no change at a midpoint T
  // species 5: H
  species[5] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4];
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
    // species 0: H2
    species[0] = -1.012521000000000e+03 * invT + 6.592218000000000e+00 -
                 3.298124000000000e+00 * tc[0] - 4.124721000000000e-04 * tc[1] +
                 1.357169166666667e-07 * tc[2] + 7.896194999999999e-12 * tc[3] -
                 2.067436000000000e-14 * tc[4];
    // species 1: O2
    species[1] = -1.005249000000000e+03 * invT - 2.821802000000000e+00 -
                 3.212936000000000e+00 * tc[0] - 5.637430000000000e-04 * tc[1] +
                 9.593583333333333e-08 * tc[2] - 1.094897500000000e-10 * tc[3] +
                 4.384277000000000e-14 * tc[4];
    // species 2: O
    species[2] = +2.914764000000000e+04 * invT - 1.756599999999997e-02 -
                 2.946429000000000e+00 * tc[0] + 8.190830000000000e-04 * tc[1] -
                 4.035053333333333e-07 * tc[2] + 1.335702500000000e-10 * tc[3] -
                 1.945348000000000e-14 * tc[4];
    // species 3: OH
    species[3] = +3.346309130000000e+03 * invT + 4.815738570000000e+00 -
                 4.125305610000000e+00 * tc[0] + 1.612724695000000e-03 * tc[1] -
                 1.087941151666667e-06 * tc[2] + 4.832113691666666e-10 * tc[3] -
                 1.031186895000000e-13 * tc[4];
    // species 4: H2O
    species[4] = -3.020811000000000e+04 * invT + 7.966090000000001e-01 -
                 3.386842000000000e+00 * tc[0] - 1.737491000000000e-03 * tc[1] +
                 1.059116000000000e-06 * tc[2] - 5.807150833333333e-10 * tc[3] +
                 1.253294000000000e-13 * tc[4];
    // species 6: HO2
    species[6] = +2.948080400000000e+02 * invT + 5.851355599999999e-01 -
                 4.301798010000000e+00 * tc[0] + 2.374560255000000e-03 * tc[1] -
                 3.526381516666666e-06 * tc[2] + 2.023032450000000e-09 * tc[3] -
                 4.646125620000001e-13 * tc[4];
    // species 7: H2O2
    species[7] = -1.766315000000000e+04 * invT - 3.396609000000000e+00 -
                 3.388754000000000e+00 * tc[0] - 3.284613000000000e-03 * tc[1] +
                 2.475021666666666e-08 * tc[2] + 3.854838333333333e-10 * tc[3] -
                 1.235757500000000e-13 * tc[4];
    // species 8: CO
    species[8] = -1.431054000000000e+04 * invT - 1.586445000000000e+00 -
                 3.262452000000000e+00 * tc[0] - 7.559705000000000e-04 * tc[1] +
                 6.469591666666667e-07 * tc[2] - 4.651620000000000e-10 * tc[3] +
                 1.237475500000000e-13 * tc[4];
    // species 9: CO2
    species[9] = -4.837314000000000e+04 * invT - 7.912765000000000e+00 -
                 2.275725000000000e+00 * tc[0] - 4.961036000000000e-03 * tc[1] +
                 1.734851666666667e-06 * tc[2] - 5.722239166666667e-10 * tc[3] +
                 1.058640000000000e-13 * tc[4];
    // species 10: HCO
    species[10] = +4.159922000000000e+03 * invT - 6.085284000000000e+00 -
                  2.898330000000000e+00 * tc[0] -
                  3.099573500000000e-03 * tc[1] +
                  1.603847333333333e-06 * tc[2] -
                  9.081875000000000e-10 * tc[3] + 2.287442500000000e-13 * tc[4];
    // species 11: N2
    species[11] = -1.020900000000000e+03 * invT - 6.516950000000001e-01 -
                  3.298677000000000e+00 * tc[0] -
                  7.041200000000000e-04 * tc[1] +
                  6.605369999999999e-07 * tc[2] -
                  4.701262500000001e-10 * tc[3] + 1.222427500000000e-13 * tc[4];
  } else {
    // species 0: H2
    species[0] = -8.350340000000000e+02 * invT + 4.346533000000000e+00 -
                 2.991423000000000e+00 * tc[0] - 3.500322000000000e-04 * tc[1] +
                 9.389715000000000e-09 * tc[2] + 7.692981666666667e-13 * tc[3] -
                 7.913760000000000e-17 * tc[4];
    // species 1: O2
    species[1] = -1.233930000000000e+03 * invT + 5.084119999999999e-01 -
                 3.697578000000000e+00 * tc[0] - 3.067598500000000e-04 * tc[1] +
                 2.098070000000000e-08 * tc[2] - 1.479400833333333e-12 * tc[3] +
                 5.682175000000001e-17 * tc[4];
    // species 2: O
    species[2] = +2.923080000000000e+04 * invT - 2.378248000000000e+00 -
                 2.542060000000000e+00 * tc[0] + 1.377531000000000e-05 * tc[1] +
                 5.171338333333333e-10 * tc[2] - 3.792555833333334e-13 * tc[3] +
                 2.184026000000000e-17 * tc[4];
    // species 3: OH
    species[3] = +3.683628750000000e+03 * invT - 2.836911870000000e+00 -
                 2.864728860000000e+00 * tc[0] - 5.282522400000000e-04 * tc[1] +
                 4.318045966666667e-08 * tc[2] - 2.543488950000000e-12 * tc[3] +
                 6.659793800000000e-17 * tc[4];
    // species 4: H2O
    species[4] = -2.989921000000000e+04 * invT - 4.190671000000000e+00 -
                 2.672146000000000e+00 * tc[0] - 1.528146500000000e-03 * tc[1] +
                 1.455043333333333e-07 * tc[2] - 1.000830000000000e-11 * tc[3] +
                 3.195809000000000e-16 * tc[4];
    // species 6: HO2
    species[6] = +1.118567130000000e+02 * invT + 2.321087500000001e-01 -
                 4.017210900000000e+00 * tc[0] - 1.119910065000000e-03 * tc[1] +
                 1.056096916666667e-07 * tc[2] - 9.520530833333334e-12 * tc[3] +
                 5.395426750000000e-16 * tc[4];
    // species 7: H2O2
    species[7] = -1.800696000000000e+04 * invT + 4.072030000000000e+00 -
                 4.573167000000000e+00 * tc[0] - 2.168068000000000e-03 * tc[1] +
                 2.457815000000000e-07 * tc[2] - 1.957420000000000e-11 * tc[3] +
                 7.158270000000000e-16 * tc[4];
    // species 8: CO
    species[8] = -1.426835000000000e+04 * invT - 3.083140000000000e+00 -
                 3.025078000000000e+00 * tc[0] - 7.213445000000000e-04 * tc[1] +
                 9.384713333333334e-08 * tc[2] - 8.488174999999999e-12 * tc[3] +
                 3.455476000000000e-16 * tc[4];
    // species 9: CO2
    species[9] = -4.896696000000000e+04 * invT + 5.409018900000000e+00 -
                 4.453623000000000e+00 * tc[0] - 1.570084500000000e-03 * tc[1] +
                 2.130685000000000e-07 * tc[2] - 1.994997500000000e-11 * tc[3] +
                 8.345165000000000e-16 * tc[4];
    // species 10: HCO
    species[10] = +3.916324000000000e+03 * invT - 1.995028000000000e+00 -
                  3.557271000000000e+00 * tc[0] -
                  1.672786500000000e-03 * tc[1] +
                  2.225010000000000e-07 * tc[2] -
                  2.058810833333333e-11 * tc[3] + 8.569255000000000e-16 * tc[4];
    // species 11: N2
    species[11] = -9.227977000000000e+02 * invT - 3.053888000000000e+00 -
                  2.926640000000000e+00 * tc[0] -
                  7.439885000000000e-04 * tc[1] +
                  9.474601666666666e-08 * tc[2] -
                  8.414199999999999e-12 * tc[3] + 3.376675500000000e-16 * tc[4];
  }

  // species with no change at a midpoint T
  // species 5: H
  species[5] = +2.547163000000000e+04 * invT + 2.960117600000000e+00 -
               2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
               0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
               0.000000000000000e+00 * tc[4];
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
    // species 0: H2
    species[0] = -1.01252100e+03 * invT + 5.59221800e+00 -
                 3.29812400e+00 * tc[0] - 4.12472100e-04 * tc[1] +
                 1.35716917e-07 * tc[2] + 7.89619500e-12 * tc[3] -
                 2.06743600e-14 * tc[4];
    // species 1: O2
    species[1] = -1.00524900e+03 * invT - 3.82180200e+00 -
                 3.21293600e+00 * tc[0] - 5.63743000e-04 * tc[1] +
                 9.59358333e-08 * tc[2] - 1.09489750e-10 * tc[3] +
                 4.38427700e-14 * tc[4];
    // species 2: O
    species[2] = +2.91476400e+04 * invT - 1.01756600e+00 -
                 2.94642900e+00 * tc[0] + 8.19083000e-04 * tc[1] -
                 4.03505333e-07 * tc[2] + 1.33570250e-10 * tc[3] -
                 1.94534800e-14 * tc[4];
    // species 3: OH
    species[3] = +3.34630913e+03 * invT + 3.81573857e+00 -
                 4.12530561e+00 * tc[0] + 1.61272470e-03 * tc[1] -
                 1.08794115e-06 * tc[2] + 4.83211369e-10 * tc[3] -
                 1.03118689e-13 * tc[4];
    // species 4: H2O
    species[4] = -3.02081100e+04 * invT - 2.03391000e-01 -
                 3.38684200e+00 * tc[0] - 1.73749100e-03 * tc[1] +
                 1.05911600e-06 * tc[2] - 5.80715083e-10 * tc[3] +
                 1.25329400e-13 * tc[4];
    // species 6: HO2
    species[6] = +2.94808040e+02 * invT - 4.14864440e-01 -
                 4.30179801e+00 * tc[0] + 2.37456025e-03 * tc[1] -
                 3.52638152e-06 * tc[2] + 2.02303245e-09 * tc[3] -
                 4.64612562e-13 * tc[4];
    // species 7: H2O2
    species[7] = -1.76631500e+04 * invT - 4.39660900e+00 -
                 3.38875400e+00 * tc[0] - 3.28461300e-03 * tc[1] +
                 2.47502167e-08 * tc[2] + 3.85483833e-10 * tc[3] -
                 1.23575750e-13 * tc[4];
    // species 8: CO
    species[8] = -1.43105400e+04 * invT - 2.58644500e+00 -
                 3.26245200e+00 * tc[0] - 7.55970500e-04 * tc[1] +
                 6.46959167e-07 * tc[2] - 4.65162000e-10 * tc[3] +
                 1.23747550e-13 * tc[4];
    // species 9: CO2
    species[9] = -4.83731400e+04 * invT - 8.91276500e+00 -
                 2.27572500e+00 * tc[0] - 4.96103600e-03 * tc[1] +
                 1.73485167e-06 * tc[2] - 5.72223917e-10 * tc[3] +
                 1.05864000e-13 * tc[4];
    // species 10: HCO
    species[10] = +4.15992200e+03 * invT - 7.08528400e+00 -
                  2.89833000e+00 * tc[0] - 3.09957350e-03 * tc[1] +
                  1.60384733e-06 * tc[2] - 9.08187500e-10 * tc[3] +
                  2.28744250e-13 * tc[4];
    // species 11: N2
    species[11] = -1.02090000e+03 * invT - 1.65169500e+00 -
                  3.29867700e+00 * tc[0] - 7.04120000e-04 * tc[1] +
                  6.60537000e-07 * tc[2] - 4.70126250e-10 * tc[3] +
                  1.22242750e-13 * tc[4];
  } else {
    // species 0: H2
    species[0] = -8.35034000e+02 * invT + 3.34653300e+00 -
                 2.99142300e+00 * tc[0] - 3.50032200e-04 * tc[1] +
                 9.38971500e-09 * tc[2] + 7.69298167e-13 * tc[3] -
                 7.91376000e-17 * tc[4];
    // species 1: O2
    species[1] = -1.23393000e+03 * invT - 4.91588000e-01 -
                 3.69757800e+00 * tc[0] - 3.06759850e-04 * tc[1] +
                 2.09807000e-08 * tc[2] - 1.47940083e-12 * tc[3] +
                 5.68217500e-17 * tc[4];
    // species 2: O
    species[2] = +2.92308000e+04 * invT - 3.37824800e+00 -
                 2.54206000e+00 * tc[0] + 1.37753100e-05 * tc[1] +
                 5.17133833e-10 * tc[2] - 3.79255583e-13 * tc[3] +
                 2.18402600e-17 * tc[4];
    // species 3: OH
    species[3] = +3.68362875e+03 * invT - 3.83691187e+00 -
                 2.86472886e+00 * tc[0] - 5.28252240e-04 * tc[1] +
                 4.31804597e-08 * tc[2] - 2.54348895e-12 * tc[3] +
                 6.65979380e-17 * tc[4];
    // species 4: H2O
    species[4] = -2.98992100e+04 * invT - 5.19067100e+00 -
                 2.67214600e+00 * tc[0] - 1.52814650e-03 * tc[1] +
                 1.45504333e-07 * tc[2] - 1.00083000e-11 * tc[3] +
                 3.19580900e-16 * tc[4];
    // species 6: HO2
    species[6] = +1.11856713e+02 * invT - 7.67891250e-01 -
                 4.01721090e+00 * tc[0] - 1.11991006e-03 * tc[1] +
                 1.05609692e-07 * tc[2] - 9.52053083e-12 * tc[3] +
                 5.39542675e-16 * tc[4];
    // species 7: H2O2
    species[7] = -1.80069600e+04 * invT + 3.07203000e+00 -
                 4.57316700e+00 * tc[0] - 2.16806800e-03 * tc[1] +
                 2.45781500e-07 * tc[2] - 1.95742000e-11 * tc[3] +
                 7.15827000e-16 * tc[4];
    // species 8: CO
    species[8] = -1.42683500e+04 * invT - 4.08314000e+00 -
                 3.02507800e+00 * tc[0] - 7.21344500e-04 * tc[1] +
                 9.38471333e-08 * tc[2] - 8.48817500e-12 * tc[3] +
                 3.45547600e-16 * tc[4];
    // species 9: CO2
    species[9] = -4.89669600e+04 * invT + 4.40901890e+00 -
                 4.45362300e+00 * tc[0] - 1.57008450e-03 * tc[1] +
                 2.13068500e-07 * tc[2] - 1.99499750e-11 * tc[3] +
                 8.34516500e-16 * tc[4];
    // species 10: HCO
    species[10] = +3.91632400e+03 * invT - 2.99502800e+00 -
                  3.55727100e+00 * tc[0] - 1.67278650e-03 * tc[1] +
                  2.22501000e-07 * tc[2] - 2.05881083e-11 * tc[3] +
                  8.56925500e-16 * tc[4];
    // species 11: N2
    species[11] = -9.22797700e+02 * invT - 4.05388800e+00 -
                  2.92664000e+00 * tc[0] - 7.43988500e-04 * tc[1] +
                  9.47460167e-08 * tc[2] - 8.41420000e-12 * tc[3] +
                  3.37667550e-16 * tc[4];
  }

  // species with no change at a midpoint T
  // species 5: H
  species[5] = +2.54716300e+04 * invT + 1.96011760e+00 -
               2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
               0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
               0.00000000e+00 * tc[4];
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
    // species 0: H2
    species[0] = +2.29812400e+00 + 4.12472100e-04 * tc[1] -
                 2.71433833e-07 * tc[2] - 2.36885850e-11 * tc[3] +
                 8.26974400e-14 * tc[4] - 1.01252100e+03 * invT;
    // species 1: O2
    species[1] = +2.21293600e+00 + 5.63743000e-04 * tc[1] -
                 1.91871667e-07 * tc[2] + 3.28469250e-10 * tc[3] -
                 1.75371080e-13 * tc[4] - 1.00524900e+03 * invT;
    // species 2: O
    species[2] = +1.94642900e+00 - 8.19083000e-04 * tc[1] +
                 8.07010667e-07 * tc[2] - 4.00710750e-10 * tc[3] +
                 7.78139200e-14 * tc[4] + 2.91476400e+04 * invT;
    // species 3: OH
    species[3] = +3.12530561e+00 - 1.61272470e-03 * tc[1] +
                 2.17588230e-06 * tc[2] - 1.44963411e-09 * tc[3] +
                 4.12474758e-13 * tc[4] + 3.34630913e+03 * invT;
    // species 4: H2O
    species[4] = +2.38684200e+00 + 1.73749100e-03 * tc[1] -
                 2.11823200e-06 * tc[2] + 1.74214525e-09 * tc[3] -
                 5.01317600e-13 * tc[4] - 3.02081100e+04 * invT;
    // species 6: HO2
    species[6] = +3.30179801e+00 - 2.37456025e-03 * tc[1] +
                 7.05276303e-06 * tc[2] - 6.06909735e-09 * tc[3] +
                 1.85845025e-12 * tc[4] + 2.94808040e+02 * invT;
    // species 7: H2O2
    species[7] = +2.38875400e+00 + 3.28461300e-03 * tc[1] -
                 4.95004333e-08 * tc[2] - 1.15645150e-09 * tc[3] +
                 4.94303000e-13 * tc[4] - 1.76631500e+04 * invT;
    // species 8: CO
    species[8] = +2.26245200e+00 + 7.55970500e-04 * tc[1] -
                 1.29391833e-06 * tc[2] + 1.39548600e-09 * tc[3] -
                 4.94990200e-13 * tc[4] - 1.43105400e+04 * invT;
    // species 9: CO2
    species[9] = +1.27572500e+00 + 4.96103600e-03 * tc[1] -
                 3.46970333e-06 * tc[2] + 1.71667175e-09 * tc[3] -
                 4.23456000e-13 * tc[4] - 4.83731400e+04 * invT;
    // species 10: HCO
    species[10] = +1.89833000e+00 + 3.09957350e-03 * tc[1] -
                  3.20769467e-06 * tc[2] + 2.72456250e-09 * tc[3] -
                  9.14977000e-13 * tc[4] + 4.15992200e+03 * invT;
    // species 11: N2
    species[11] = +2.29867700e+00 + 7.04120000e-04 * tc[1] -
                  1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
                  4.88971000e-13 * tc[4] - 1.02090000e+03 * invT;
  } else {
    // species 0: H2
    species[0] = +1.99142300e+00 + 3.50032200e-04 * tc[1] -
                 1.87794300e-08 * tc[2] - 2.30789450e-12 * tc[3] +
                 3.16550400e-16 * tc[4] - 8.35034000e+02 * invT;
    // species 1: O2
    species[1] = +2.69757800e+00 + 3.06759850e-04 * tc[1] -
                 4.19614000e-08 * tc[2] + 4.43820250e-12 * tc[3] -
                 2.27287000e-16 * tc[4] - 1.23393000e+03 * invT;
    // species 2: O
    species[2] = +1.54206000e+00 - 1.37753100e-05 * tc[1] -
                 1.03426767e-09 * tc[2] + 1.13776675e-12 * tc[3] -
                 8.73610400e-17 * tc[4] + 2.92308000e+04 * invT;
    // species 3: OH
    species[3] = +1.86472886e+00 + 5.28252240e-04 * tc[1] -
                 8.63609193e-08 * tc[2] + 7.63046685e-12 * tc[3] -
                 2.66391752e-16 * tc[4] + 3.68362875e+03 * invT;
    // species 4: H2O
    species[4] = +1.67214600e+00 + 1.52814650e-03 * tc[1] -
                 2.91008667e-07 * tc[2] + 3.00249000e-11 * tc[3] -
                 1.27832360e-15 * tc[4] - 2.98992100e+04 * invT;
    // species 6: HO2
    species[6] = +3.01721090e+00 + 1.11991006e-03 * tc[1] -
                 2.11219383e-07 * tc[2] + 2.85615925e-11 * tc[3] -
                 2.15817070e-15 * tc[4] + 1.11856713e+02 * invT;
    // species 7: H2O2
    species[7] = +3.57316700e+00 + 2.16806800e-03 * tc[1] -
                 4.91563000e-07 * tc[2] + 5.87226000e-11 * tc[3] -
                 2.86330800e-15 * tc[4] - 1.80069600e+04 * invT;
    // species 8: CO
    species[8] = +2.02507800e+00 + 7.21344500e-04 * tc[1] -
                 1.87694267e-07 * tc[2] + 2.54645250e-11 * tc[3] -
                 1.38219040e-15 * tc[4] - 1.42683500e+04 * invT;
    // species 9: CO2
    species[9] = +3.45362300e+00 + 1.57008450e-03 * tc[1] -
                 4.26137000e-07 * tc[2] + 5.98499250e-11 * tc[3] -
                 3.33806600e-15 * tc[4] - 4.89669600e+04 * invT;
    // species 10: HCO
    species[10] = +2.55727100e+00 + 1.67278650e-03 * tc[1] -
                  4.45002000e-07 * tc[2] + 6.17643250e-11 * tc[3] -
                  3.42770200e-15 * tc[4] + 3.91632400e+03 * invT;
    // species 11: N2
    species[11] = +1.92664000e+00 + 7.43988500e-04 * tc[1] -
                  1.89492033e-07 * tc[2] + 2.52426000e-11 * tc[3] -
                  1.35067020e-15 * tc[4] - 9.22797700e+02 * invT;
  }

  // species with no change at a midpoint T
  // species 5: H
  species[5] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] + 2.54716300e+04 * invT;
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
    // species 0: H2
    species[0] = +3.29812400e+00 + 4.12472100e-04 * tc[1] -
                 2.71433833e-07 * tc[2] - 2.36885850e-11 * tc[3] +
                 8.26974400e-14 * tc[4] - 1.01252100e+03 * invT;
    // species 1: O2
    species[1] = +3.21293600e+00 + 5.63743000e-04 * tc[1] -
                 1.91871667e-07 * tc[2] + 3.28469250e-10 * tc[3] -
                 1.75371080e-13 * tc[4] - 1.00524900e+03 * invT;
    // species 2: O
    species[2] = +2.94642900e+00 - 8.19083000e-04 * tc[1] +
                 8.07010667e-07 * tc[2] - 4.00710750e-10 * tc[3] +
                 7.78139200e-14 * tc[4] + 2.91476400e+04 * invT;
    // species 3: OH
    species[3] = +4.12530561e+00 - 1.61272470e-03 * tc[1] +
                 2.17588230e-06 * tc[2] - 1.44963411e-09 * tc[3] +
                 4.12474758e-13 * tc[4] + 3.34630913e+03 * invT;
    // species 4: H2O
    species[4] = +3.38684200e+00 + 1.73749100e-03 * tc[1] -
                 2.11823200e-06 * tc[2] + 1.74214525e-09 * tc[3] -
                 5.01317600e-13 * tc[4] - 3.02081100e+04 * invT;
    // species 6: HO2
    species[6] = +4.30179801e+00 - 2.37456025e-03 * tc[1] +
                 7.05276303e-06 * tc[2] - 6.06909735e-09 * tc[3] +
                 1.85845025e-12 * tc[4] + 2.94808040e+02 * invT;
    // species 7: H2O2
    species[7] = +3.38875400e+00 + 3.28461300e-03 * tc[1] -
                 4.95004333e-08 * tc[2] - 1.15645150e-09 * tc[3] +
                 4.94303000e-13 * tc[4] - 1.76631500e+04 * invT;
    // species 8: CO
    species[8] = +3.26245200e+00 + 7.55970500e-04 * tc[1] -
                 1.29391833e-06 * tc[2] + 1.39548600e-09 * tc[3] -
                 4.94990200e-13 * tc[4] - 1.43105400e+04 * invT;
    // species 9: CO2
    species[9] = +2.27572500e+00 + 4.96103600e-03 * tc[1] -
                 3.46970333e-06 * tc[2] + 1.71667175e-09 * tc[3] -
                 4.23456000e-13 * tc[4] - 4.83731400e+04 * invT;
    // species 10: HCO
    species[10] = +2.89833000e+00 + 3.09957350e-03 * tc[1] -
                  3.20769467e-06 * tc[2] + 2.72456250e-09 * tc[3] -
                  9.14977000e-13 * tc[4] + 4.15992200e+03 * invT;
    // species 11: N2
    species[11] = +3.29867700e+00 + 7.04120000e-04 * tc[1] -
                  1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
                  4.88971000e-13 * tc[4] - 1.02090000e+03 * invT;
  } else {
    // species 0: H2
    species[0] = +2.99142300e+00 + 3.50032200e-04 * tc[1] -
                 1.87794300e-08 * tc[2] - 2.30789450e-12 * tc[3] +
                 3.16550400e-16 * tc[4] - 8.35034000e+02 * invT;
    // species 1: O2
    species[1] = +3.69757800e+00 + 3.06759850e-04 * tc[1] -
                 4.19614000e-08 * tc[2] + 4.43820250e-12 * tc[3] -
                 2.27287000e-16 * tc[4] - 1.23393000e+03 * invT;
    // species 2: O
    species[2] = +2.54206000e+00 - 1.37753100e-05 * tc[1] -
                 1.03426767e-09 * tc[2] + 1.13776675e-12 * tc[3] -
                 8.73610400e-17 * tc[4] + 2.92308000e+04 * invT;
    // species 3: OH
    species[3] = +2.86472886e+00 + 5.28252240e-04 * tc[1] -
                 8.63609193e-08 * tc[2] + 7.63046685e-12 * tc[3] -
                 2.66391752e-16 * tc[4] + 3.68362875e+03 * invT;
    // species 4: H2O
    species[4] = +2.67214600e+00 + 1.52814650e-03 * tc[1] -
                 2.91008667e-07 * tc[2] + 3.00249000e-11 * tc[3] -
                 1.27832360e-15 * tc[4] - 2.98992100e+04 * invT;
    // species 6: HO2
    species[6] = +4.01721090e+00 + 1.11991006e-03 * tc[1] -
                 2.11219383e-07 * tc[2] + 2.85615925e-11 * tc[3] -
                 2.15817070e-15 * tc[4] + 1.11856713e+02 * invT;
    // species 7: H2O2
    species[7] = +4.57316700e+00 + 2.16806800e-03 * tc[1] -
                 4.91563000e-07 * tc[2] + 5.87226000e-11 * tc[3] -
                 2.86330800e-15 * tc[4] - 1.80069600e+04 * invT;
    // species 8: CO
    species[8] = +3.02507800e+00 + 7.21344500e-04 * tc[1] -
                 1.87694267e-07 * tc[2] + 2.54645250e-11 * tc[3] -
                 1.38219040e-15 * tc[4] - 1.42683500e+04 * invT;
    // species 9: CO2
    species[9] = +4.45362300e+00 + 1.57008450e-03 * tc[1] -
                 4.26137000e-07 * tc[2] + 5.98499250e-11 * tc[3] -
                 3.33806600e-15 * tc[4] - 4.89669600e+04 * invT;
    // species 10: HCO
    species[10] = +3.55727100e+00 + 1.67278650e-03 * tc[1] -
                  4.45002000e-07 * tc[2] + 6.17643250e-11 * tc[3] -
                  3.42770200e-15 * tc[4] + 3.91632400e+03 * invT;
    // species 11: N2
    species[11] = +2.92664000e+00 + 7.43988500e-04 * tc[1] -
                  1.89492033e-07 * tc[2] + 2.52426000e-11 * tc[3] -
                  1.35067020e-15 * tc[4] - 9.22797700e+02 * invT;
  }

  // species with no change at a midpoint T
  // species 5: H
  species[5] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] + 2.54716300e+04 * invT;
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
    // species 0: H2
    species[0] = +3.29812400e+00 * tc[0] + 8.24944200e-04 * tc[1] -
                 4.07150750e-07 * tc[2] - 3.15847800e-11 * tc[3] +
                 1.03371800e-13 * tc[4] - 3.29409400e+00;
    // species 1: O2
    species[1] = +3.21293600e+00 * tc[0] + 1.12748600e-03 * tc[1] -
                 2.87807500e-07 * tc[2] + 4.37959000e-10 * tc[3] -
                 2.19213850e-13 * tc[4] + 6.03473800e+00;
    // species 2: O
    species[2] = +2.94642900e+00 * tc[0] - 1.63816600e-03 * tc[1] +
                 1.21051600e-06 * tc[2] - 5.34281000e-10 * tc[3] +
                 9.72674000e-14 * tc[4] + 2.96399500e+00;
    // species 3: OH
    species[3] = +4.12530561e+00 * tc[0] - 3.22544939e-03 * tc[1] +
                 3.26382346e-06 * tc[2] - 1.93284548e-09 * tc[3] +
                 5.15593447e-13 * tc[4] - 6.90432960e-01;
    // species 4: H2O
    species[4] = +3.38684200e+00 * tc[0] + 3.47498200e-03 * tc[1] -
                 3.17734800e-06 * tc[2] + 2.32286033e-09 * tc[3] -
                 6.26647000e-13 * tc[4] + 2.59023300e+00;
    // species 6: HO2
    species[6] = +4.30179801e+00 * tc[0] - 4.74912051e-03 * tc[1] +
                 1.05791445e-05 * tc[2] - 8.09212980e-09 * tc[3] +
                 2.32306281e-12 * tc[4] + 3.71666245e+00;
    // species 7: H2O2
    species[7] = +3.38875400e+00 * tc[0] + 6.56922600e-03 * tc[1] -
                 7.42506500e-08 * tc[2] - 1.54193533e-09 * tc[3] +
                 6.17878750e-13 * tc[4] + 6.78536300e+00;
    // species 8: CO
    species[8] = +3.26245200e+00 * tc[0] + 1.51194100e-03 * tc[1] -
                 1.94087750e-06 * tc[2] + 1.86064800e-09 * tc[3] -
                 6.18737750e-13 * tc[4] + 4.84889700e+00;
    // species 9: CO2
    species[9] = +2.27572500e+00 * tc[0] + 9.92207200e-03 * tc[1] -
                 5.20455500e-06 * tc[2] + 2.28889567e-09 * tc[3] -
                 5.29320000e-13 * tc[4] + 1.01884900e+01;
    // species 10: HCO
    species[10] = +2.89833000e+00 * tc[0] + 6.19914700e-03 * tc[1] -
                  4.81154200e-06 * tc[2] + 3.63275000e-09 * tc[3] -
                  1.14372125e-12 * tc[4] + 8.98361400e+00;
    // species 11: N2
    species[11] = +3.29867700e+00 * tc[0] + 1.40824000e-03 * tc[1] -
                  1.98161100e-06 * tc[2] + 1.88050500e-09 * tc[3] -
                  6.11213750e-13 * tc[4] + 3.95037200e+00;
  } else {
    // species 0: H2
    species[0] = +2.99142300e+00 * tc[0] + 7.00064400e-04 * tc[1] -
                 2.81691450e-08 * tc[2] - 3.07719267e-12 * tc[3] +
                 3.95688000e-16 * tc[4] - 1.35511000e+00;
    // species 1: O2
    species[1] = +3.69757800e+00 * tc[0] + 6.13519700e-04 * tc[1] -
                 6.29421000e-08 * tc[2] + 5.91760333e-12 * tc[3] -
                 2.84108750e-16 * tc[4] + 3.18916600e+00;
    // species 2: O
    species[2] = +2.54206000e+00 * tc[0] - 2.75506200e-05 * tc[1] -
                 1.55140150e-09 * tc[2] + 1.51702233e-12 * tc[3] -
                 1.09201300e-16 * tc[4] + 4.92030800e+00;
    // species 3: OH
    species[3] = +2.86472886e+00 * tc[0] + 1.05650448e-03 * tc[1] -
                 1.29541379e-07 * tc[2] + 1.01739558e-11 * tc[3] -
                 3.32989690e-16 * tc[4] + 5.70164073e+00;
    // species 4: H2O
    species[4] = +2.67214600e+00 * tc[0] + 3.05629300e-03 * tc[1] -
                 4.36513000e-07 * tc[2] + 4.00332000e-11 * tc[3] -
                 1.59790450e-15 * tc[4] + 6.86281700e+00;
    // species 6: HO2
    species[6] = +4.01721090e+00 * tc[0] + 2.23982013e-03 * tc[1] -
                 3.16829075e-07 * tc[2] + 3.80821233e-11 * tc[3] -
                 2.69771337e-15 * tc[4] + 3.78510215e+00;
    // species 7: H2O2
    species[7] = +4.57316700e+00 * tc[0] + 4.33613600e-03 * tc[1] -
                 7.37344500e-07 * tc[2] + 7.82968000e-11 * tc[3] -
                 3.57913500e-15 * tc[4] + 5.01137000e-01;
    // species 8: CO
    species[8] = +3.02507800e+00 * tc[0] + 1.44268900e-03 * tc[1] -
                 2.81541400e-07 * tc[2] + 3.39527000e-11 * tc[3] -
                 1.72773800e-15 * tc[4] + 6.10821800e+00;
    // species 9: CO2
    species[9] = +4.45362300e+00 * tc[0] + 3.14016900e-03 * tc[1] -
                 6.39205500e-07 * tc[2] + 7.97999000e-11 * tc[3] -
                 4.17258250e-15 * tc[4] - 9.55395900e-01;
    // species 10: HCO
    species[10] = +3.55727100e+00 * tc[0] + 3.34557300e-03 * tc[1] -
                  6.67503000e-07 * tc[2] + 8.23524333e-11 * tc[3] -
                  4.28462750e-15 * tc[4] + 5.55229900e+00;
    // species 11: N2
    species[11] = +2.92664000e+00 * tc[0] + 1.48797700e-03 * tc[1] -
                  2.84238050e-07 * tc[2] + 3.36568000e-11 * tc[3] -
                  1.68833775e-15 * tc[4] + 5.98052800e+00;
  }

  // species with no change at a midpoint T
  // species 5: H
  species[5] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 4.60117600e-01;
}

// Returns the mean specific heat at CP (Eq. 33)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBL(const amrex::Real* T, const amrex::Real* x, amrex::Real* cpbl)
{
  int id; // loop counter
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real cpor[12];                               //  temporary storage
  cp_R(cpor, tc);

  // perform dot product
  for (id = 0; id < 12; ++id) {
    result += x[id] * cpor[id];
  }

  *cpbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CP (Eq. 34)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBS(const amrex::Real* T, const amrex::Real* y, amrex::Real* cpbs)
{
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real cpor[12], tresult[12];                  //  temporary storage
  amrex::Real imw[12];

  get_imw(imw);

  cp_R(cpor, tc);

  for (int i = 0; i < 12; i++) {
    tresult[i] = cpor[i] * y[i] * imw[i];
  }
  for (const auto& tr : tresult) {
    result += tr;
  }

  *cpbs = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 35)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBL(const amrex::Real* T, const amrex::Real* x, amrex::Real* cvbl)
{
  int id; // loop counter
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real cvor[12];                               //  temporary storage
  cv_R(cvor, tc);

  // perform dot product
  for (id = 0; id < 12; ++id) {
    result += x[id] * cvor[id];
  }

  *cvbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 36)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBS(const amrex::Real* T, const amrex::Real* y, amrex::Real* cvbs)
{
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real cvor[12];                               //  temporary storage
  amrex::Real imw[12];

  get_imw(imw);

  cv_R(cvor, tc);

  // multiply by y/molecularweight
  result += cvor[0] * y[0] * imw[0];    // H2
  result += cvor[1] * y[1] * imw[1];    // O2
  result += cvor[2] * y[2] * imw[2];    // O
  result += cvor[3] * y[3] * imw[3];    // OH
  result += cvor[4] * y[4] * imw[4];    // H2O
  result += cvor[5] * y[5] * imw[5];    // H
  result += cvor[6] * y[6] * imw[6];    // HO2
  result += cvor[7] * y[7] * imw[7];    // H2O2
  result += cvor[8] * y[8] * imw[8];    // CO
  result += cvor[9] * y[9] * imw[9];    // CO2
  result += cvor[10] * y[10] * imw[10]; // HCO
  result += cvor[11] * y[11] * imw[11]; // N2

  *cvbs = result * 8.31446261815324e+07;
}

// Returns the mean enthalpy of the mixture in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBML(const amrex::Real* T, const amrex::Real* x, amrex::Real* hbml)
{
  int id; // loop counter
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real hml[12];                                //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesEnthalpy(hml, tc);

  // perform dot product
  for (id = 0; id < 12; ++id) {
    result += x[id] * hml[id];
  }

  *hbml = result * RT;
}

// Returns mean enthalpy of mixture in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBMS(const amrex::Real* T, const amrex::Real* y, amrex::Real* hbms)
{
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real hml[12], tmp[12];                       //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[12];

  get_imw(imw);

  speciesEnthalpy(hml, tc);

  int id;
  for (id = 0; id < 12; ++id) {
    tmp[id] = y[id] * hml[id] * imw[id];
  }
  for (id = 0; id < 12; ++id) {
    result += tmp[id];
  }

  *hbms = result * RT;
}

// get mean internal energy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBML(const amrex::Real* T, const amrex::Real* x, amrex::Real* ubml)
{
  int id; // loop counter
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real uml[12];                                //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesInternalEnergy(uml, tc);

  // perform dot product
  for (id = 0; id < 12; ++id) {
    result += x[id] * uml[id];
  }

  *ubml = result * RT;
}

// get mean internal energy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBMS(const amrex::Real* T, const amrex::Real* y, amrex::Real* ubms)
{
  amrex::Real result = 0;
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real ums[12];                                //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[12];

  get_imw(imw);

  speciesInternalEnergy(ums, tc);

  // perform dot product + scaling by wt
  result += y[0] * ums[0] * imw[0];    // H2
  result += y[1] * ums[1] * imw[1];    // O2
  result += y[2] * ums[2] * imw[2];    // O
  result += y[3] * ums[3] * imw[3];    // OH
  result += y[4] * ums[4] * imw[4];    // H2O
  result += y[5] * ums[5] * imw[5];    // H
  result += y[6] * ums[6] * imw[6];    // HO2
  result += y[7] * ums[7] * imw[7];    // H2O2
  result += y[8] * ums[8] * imw[8];    // CO
  result += y[9] * ums[9] * imw[9];    // CO2
  result += y[10] * ums[10] * imw[10]; // HCO
  result += y[11] * ums[11] * imw[11]; // N2

  *ubms = result * RT;
}

// get mixture entropy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSBML(
  const amrex::Real* P,
  const amrex::Real* T,
  const amrex::Real* x,
  amrex::Real* sbml)
{
  int id; // loop counter
  amrex::Real result = 0;
  // Log of normalized pressure in cgs units dynes/cm^2 by Patm
  amrex::Real logPratio = log(*P / 1013250.0);
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real sor[12]; //  temporary storage
  speciesEntropy(sor, tc);

  // Compute Eq 42
  for (id = 0; id < 12; ++id) {
    result += x[id] * (sor[id] - log((x[id] + 1e-100)) - logPratio);
  }

  *sbml = result * 8.31446261815324e+07;
}

// get mixture entropy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSBMS(
  const amrex::Real* P,
  const amrex::Real* T,
  const amrex::Real* y,
  amrex::Real* sbms)
{
  amrex::Real result = 0;
  // Log of normalized pressure in cgs units dynes/cm^2 by Patm
  amrex::Real logPratio = log(*P / 1013250.0);
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real sor[12]; //  temporary storage
  amrex::Real x[12];   //  need a ytx conversion
  amrex::Real YOW = 0; // See Eq 4, 6 in CK Manual
  amrex::Real imw[12];

  get_imw(imw);

  // Compute inverse of mean molecular wt first
  YOW += y[0] * imw[0];   // H2
  YOW += y[1] * imw[1];   // O2
  YOW += y[2] * imw[2];   // O
  YOW += y[3] * imw[3];   // OH
  YOW += y[4] * imw[4];   // H2O
  YOW += y[5] * imw[5];   // H
  YOW += y[6] * imw[6];   // HO2
  YOW += y[7] * imw[7];   // H2O2
  YOW += y[8] * imw[8];   // CO
  YOW += y[9] * imw[9];   // CO2
  YOW += y[10] * imw[10]; // HCO
  YOW += y[11] * imw[11]; // N2
  // Now compute y to x conversion
  x[0] = y[0] / (2.016000 * YOW);
  x[1] = y[1] / (31.998000 * YOW);
  x[2] = y[2] / (15.999000 * YOW);
  x[3] = y[3] / (17.007000 * YOW);
  x[4] = y[4] / (18.015000 * YOW);
  x[5] = y[5] / (1.008000 * YOW);
  x[6] = y[6] / (33.006000 * YOW);
  x[7] = y[7] / (34.014000 * YOW);
  x[8] = y[8] / (28.010000 * YOW);
  x[9] = y[9] / (44.009000 * YOW);
  x[10] = y[10] / (29.018000 * YOW);
  x[11] = y[11] / (28.014000 * YOW);
  speciesEntropy(sor, tc);
  // Perform computation in Eq 42 and 43
  result += x[0] * (sor[0] - log((x[0] + 1e-100)) - logPratio);
  result += x[1] * (sor[1] - log((x[1] + 1e-100)) - logPratio);
  result += x[2] * (sor[2] - log((x[2] + 1e-100)) - logPratio);
  result += x[3] * (sor[3] - log((x[3] + 1e-100)) - logPratio);
  result += x[4] * (sor[4] - log((x[4] + 1e-100)) - logPratio);
  result += x[5] * (sor[5] - log((x[5] + 1e-100)) - logPratio);
  result += x[6] * (sor[6] - log((x[6] + 1e-100)) - logPratio);
  result += x[7] * (sor[7] - log((x[7] + 1e-100)) - logPratio);
  result += x[8] * (sor[8] - log((x[8] + 1e-100)) - logPratio);
  result += x[9] * (sor[9] - log((x[9] + 1e-100)) - logPratio);
  result += x[10] * (sor[10] - log((x[10] + 1e-100)) - logPratio);
  result += x[11] * (sor[11] - log((x[11] + 1e-100)) - logPratio);
  // Scale by R/W
  *sbms = result * 8.31446261815324e+07 * YOW;
}

//  get temperature given internal energy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_EY(
  const amrex::Real* e, const amrex::Real* y, amrex::Real* t, int* ierr)
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
    } else if (fabs(dt) < tol) {
      break;
    }
    t1 += dt;
  }
  *t = t1;
  *ierr = 0;
}

//  get temperature given enthalpy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_HY(
  const amrex::Real* h, const amrex::Real* y, amrex::Real* t, int* ierr)
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
    } else if (fabs(dt) < tol) {
      break;
    }
    t1 += dt;
  }
  *t = t1;
  *ierr = 0;
}

// Compute P = rhoRT/W(x)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPX(
  const amrex::Real* rho,
  const amrex::Real* T,
  const amrex::Real* x,
  amrex::Real* P)
{
  amrex::Real XW = 0;                           //  To hold mean molecular wt
  XW += x[0] * 2.016000;                        // H2
  XW += x[1] * 31.998000;                       // O2
  XW += x[2] * 15.999000;                       // O
  XW += x[3] * 17.007000;                       // OH
  XW += x[4] * 18.015000;                       // H2O
  XW += x[5] * 1.008000;                        // H
  XW += x[6] * 33.006000;                       // HO2
  XW += x[7] * 34.014000;                       // H2O2
  XW += x[8] * 28.010000;                       // CO
  XW += x[9] * 44.009000;                       // CO2
  XW += x[10] * 29.018000;                      // HCO
  XW += x[11] * 28.014000;                      // N2
  *P = *rho * 8.31446261815324e+07 * (*T) / XW; // P = rho*R*T/W
}

// Compute P = rhoRT/W(y)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPY(
  const amrex::Real* rho,
  const amrex::Real* T,
  const amrex::Real* y,
  amrex::Real* P)
{
  amrex::Real YOW = 0; //  for computing mean MW
  amrex::Real imw[12];

  get_imw(imw);

  YOW += y[0] * imw[0];                          // H2
  YOW += y[1] * imw[1];                          // O2
  YOW += y[2] * imw[2];                          // O
  YOW += y[3] * imw[3];                          // OH
  YOW += y[4] * imw[4];                          // H2O
  YOW += y[5] * imw[5];                          // H
  YOW += y[6] * imw[6];                          // HO2
  YOW += y[7] * imw[7];                          // H2O2
  YOW += y[8] * imw[8];                          // CO
  YOW += y[9] * imw[9];                          // CO2
  YOW += y[10] * imw[10];                        // HCO
  YOW += y[11] * imw[11];                        // N2
  *P = *rho * 8.31446261815324e+07 * (*T) * YOW; // P = rho*R*T/W
}

// Compute P = rhoRT/W(c)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPC(
  const amrex::Real* rho,
  const amrex::Real* T,
  const amrex::Real* c,
  amrex::Real* P)
{
  int id; // loop counter
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 2.016000;   // H2
  W += c[1] * 31.998000;  // O2
  W += c[2] * 15.999000;  // O
  W += c[3] * 17.007000;  // OH
  W += c[4] * 18.015000;  // H2O
  W += c[5] * 1.008000;   // H
  W += c[6] * 33.006000;  // HO2
  W += c[7] * 34.014000;  // H2O2
  W += c[8] * 28.010000;  // CO
  W += c[9] * 44.009000;  // CO2
  W += c[10] * 29.018000; // HCO
  W += c[11] * 28.014000; // N2

  for (id = 0; id < 12; ++id) {
    sumC += c[id];
  }
  *P = *rho * 8.31446261815324e+07 * (*T) * sumC / W; // P = rho*R*T/W
}

// Compute rho = PW(x)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOX(
  const amrex::Real* P,
  const amrex::Real* T,
  const amrex::Real* x,
  amrex::Real* rho)
{
  amrex::Real XW = 0;                             //  To hold mean molecular wt
  XW += x[0] * 2.016000;                          // H2
  XW += x[1] * 31.998000;                         // O2
  XW += x[2] * 15.999000;                         // O
  XW += x[3] * 17.007000;                         // OH
  XW += x[4] * 18.015000;                         // H2O
  XW += x[5] * 1.008000;                          // H
  XW += x[6] * 33.006000;                         // HO2
  XW += x[7] * 34.014000;                         // H2O2
  XW += x[8] * 28.010000;                         // CO
  XW += x[9] * 44.009000;                         // CO2
  XW += x[10] * 29.018000;                        // HCO
  XW += x[11] * 28.014000;                        // N2
  *rho = *P * XW / (8.31446261815324e+07 * (*T)); // rho = P*W/(R*T)
}

// Compute rho = P*W(y)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOY(
  const amrex::Real* P,
  const amrex::Real* T,
  const amrex::Real* y,
  amrex::Real* rho)
{
  amrex::Real YOW = 0;
  amrex::Real tmp[12];
  amrex::Real imw[12];

  get_imw(imw);

  for (int i = 0; i < 12; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (const auto& tp : tmp) {
    YOW += tp;
  }

  *rho = *P / (8.31446261815324e+07 * (*T) * YOW); // rho = P*W/(R*T)
}

// Compute rho = P*W(c)/(R*T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOC(
  const amrex::Real* P,
  const amrex::Real* T,
  const amrex::Real* c,
  amrex::Real* rho)
{
  int id; // loop counter
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 2.016000;   // H2
  W += c[1] * 31.998000;  // O2
  W += c[2] * 15.999000;  // O
  W += c[3] * 17.007000;  // OH
  W += c[4] * 18.015000;  // H2O
  W += c[5] * 1.008000;   // H
  W += c[6] * 33.006000;  // HO2
  W += c[7] * 34.014000;  // H2O2
  W += c[8] * 28.010000;  // CO
  W += c[9] * 44.009000;  // CO2
  W += c[10] * 29.018000; // HCO
  W += c[11] * 28.014000; // N2

  for (id = 0; id < 12; ++id) {
    sumC += c[id];
  }
  *rho = *P * W / (sumC * (*T) * 8.31446261815324e+07); // rho = PW/(R*T)
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
CKMMWY(const amrex::Real* y, amrex::Real* wtm)
{
  amrex::Real YOW = 0;
  amrex::Real tmp[12];
  amrex::Real imw[12];

  get_imw(imw);

  for (int i = 0; i < 12; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (const auto& tp : tmp) {
    YOW += tp;
  }

  *wtm = 1.0 / YOW;
}

// given x[species]: mole fractions
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWX(const amrex::Real* x, amrex::Real* wtm)
{
  amrex::Real XW = 0;      //  see Eq 4 in CK Manual
  XW += x[0] * 2.016000;   // H2
  XW += x[1] * 31.998000;  // O2
  XW += x[2] * 15.999000;  // O
  XW += x[3] * 17.007000;  // OH
  XW += x[4] * 18.015000;  // H2O
  XW += x[5] * 1.008000;   // H
  XW += x[6] * 33.006000;  // HO2
  XW += x[7] * 34.014000;  // H2O2
  XW += x[8] * 28.010000;  // CO
  XW += x[9] * 44.009000;  // CO2
  XW += x[10] * 29.018000; // HCO
  XW += x[11] * 28.014000; // N2
  *wtm = XW;
}

// given c[species]: molar concentration
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWC(const amrex::Real* c, amrex::Real* wtm)
{
  int id; // loop counter
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 2.016000;   // H2
  W += c[1] * 31.998000;  // O2
  W += c[2] * 15.999000;  // O
  W += c[3] * 17.007000;  // OH
  W += c[4] * 18.015000;  // H2O
  W += c[5] * 1.008000;   // H
  W += c[6] * 33.006000;  // HO2
  W += c[7] * 34.014000;  // H2O2
  W += c[8] * 28.010000;  // CO
  W += c[9] * 44.009000;  // CO2
  W += c[10] * 29.018000; // HCO
  W += c[11] * 28.014000; // N2

  for (id = 0; id < 12; ++id) {
    sumC += c[id];
  }
  //  CK provides no guard against divison by zero
  *wtm = W / sumC;
}

// get Cp/R as a function of T
// for all species (Eq 19)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPOR(const amrex::Real* T, amrex::Real* cpor)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cp_R(cpor, tc);
}

// get H/RT as a function of T
// for all species (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHORT(const amrex::Real* T, amrex::Real* hort)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  speciesEnthalpy(hort, tc);
}

// get S/R as a function of T
// for all species (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSOR(const amrex::Real* T, amrex::Real* sor)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  speciesEntropy(sor, tc);
}

// convert y[species] (mass fracs) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTX(const amrex::Real* y, amrex::Real* x)
{
  amrex::Real YOW = 0;
  amrex::Real tmp[12];
  amrex::Real imw[12];

  get_imw(imw);

  for (int i = 0; i < 12; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (const auto& tp : tmp) {
    YOW += tp;
  }

  amrex::Real YOWINV = 1.0 / YOW;

  for (int i = 0; i < 12; i++) {
    x[i] = y[i] * imw[i] * YOWINV;
  }
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTCP(
  const amrex::Real* P,
  const amrex::Real* T,
  const amrex::Real* y,
  amrex::Real* c)
{
  amrex::Real YOW = 0;
  amrex::Real PWORT;
  amrex::Real imw[12];

  get_imw(imw);

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 12; i++) {
    c[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 12; i++) {
    YOW += c[i];
  }

  // PW/RT (see Eq. 7)
  PWORT = (*P) / (YOW * 8.31446261815324e+07 * (*T));
  // Now compute conversion

  for (int i = 0; i < 12; i++) {
    c[i] = PWORT * y[i] * imw[i];
  }
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTCR(
  const amrex::Real* rho,
  amrex::Real* /*T*/,
  const amrex::Real* y,
  amrex::Real* c)
{
  amrex::Real imw[12];

  get_imw(imw);

  for (int i = 0; i < 12; i++) {
    c[i] = (*rho) * y[i] * imw[i];
  }
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTY(const amrex::Real* x, amrex::Real* y)
{
  amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
  // Compute mean molecular wt first
  XW += x[0] * 2.016000;   // H2
  XW += x[1] * 31.998000;  // O2
  XW += x[2] * 15.999000;  // O
  XW += x[3] * 17.007000;  // OH
  XW += x[4] * 18.015000;  // H2O
  XW += x[5] * 1.008000;   // H
  XW += x[6] * 33.006000;  // HO2
  XW += x[7] * 34.014000;  // H2O2
  XW += x[8] * 28.010000;  // CO
  XW += x[9] * 44.009000;  // CO2
  XW += x[10] * 29.018000; // HCO
  XW += x[11] * 28.014000; // N2
  // Now compute conversion
  amrex::Real XWinv = 1.0 / XW;
  y[0] = x[0] * 2.016000 * XWinv;
  y[1] = x[1] * 31.998000 * XWinv;
  y[2] = x[2] * 15.999000 * XWinv;
  y[3] = x[3] * 17.007000 * XWinv;
  y[4] = x[4] * 18.015000 * XWinv;
  y[5] = x[5] * 1.008000 * XWinv;
  y[6] = x[6] * 33.006000 * XWinv;
  y[7] = x[7] * 34.014000 * XWinv;
  y[8] = x[8] * 28.010000 * XWinv;
  y[9] = x[9] * 44.009000 * XWinv;
  y[10] = x[10] * 29.018000 * XWinv;
  y[11] = x[11] * 28.014000 * XWinv;
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTCP(
  const amrex::Real* P,
  const amrex::Real* T,
  const amrex::Real* x,
  amrex::Real* c)
{
  int id;                                                  // loop counter
  amrex::Real PORT = (*P) / (8.31446261815324e+07 * (*T)); // P/RT

  // Compute conversion, see Eq 10
  for (id = 0; id < 12; ++id) {
    c[id] = x[id] * PORT;
  }
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTCR(
  const amrex::Real* rho,
  const amrex::Real* /*T*/,
  const amrex::Real* x,
  amrex::Real* c)
{
  int id;             // loop counter
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 2.016000;   // H2
  XW += x[1] * 31.998000;  // O2
  XW += x[2] * 15.999000;  // O
  XW += x[3] * 17.007000;  // OH
  XW += x[4] * 18.015000;  // H2O
  XW += x[5] * 1.008000;   // H
  XW += x[6] * 33.006000;  // HO2
  XW += x[7] * 34.014000;  // H2O2
  XW += x[8] * 28.010000;  // CO
  XW += x[9] * 44.009000;  // CO2
  XW += x[10] * 29.018000; // HCO
  XW += x[11] * 28.014000; // N2
  ROW = (*rho) / XW;

  // Compute conversion, see Eq 11
  for (id = 0; id < 12; ++id) {
    c[id] = x[id] * ROW;
  }
}

// convert c[species] (molar conc) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTX(const amrex::Real* c, amrex::Real* x)
{
  int id; // loop counter
  amrex::Real sumC = 0;

  // compute sum of c
  for (id = 0; id < 12; ++id) {
    sumC += c[id];
  }

  //  See Eq 13
  amrex::Real sumCinv = 1.0 / sumC;
  for (id = 0; id < 12; ++id) {
    x[id] = c[id] * sumCinv;
  }
}

// convert c[species] (molar conc) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTY(const amrex::Real* c, amrex::Real* y)
{
  amrex::Real CW = 0; // See Eq 12 in CK Manual
  // compute denominator in eq 12 first
  CW += c[0] * 2.016000;   // H2
  CW += c[1] * 31.998000;  // O2
  CW += c[2] * 15.999000;  // O
  CW += c[3] * 17.007000;  // OH
  CW += c[4] * 18.015000;  // H2O
  CW += c[5] * 1.008000;   // H
  CW += c[6] * 33.006000;  // HO2
  CW += c[7] * 34.014000;  // H2O2
  CW += c[8] * 28.010000;  // CO
  CW += c[9] * 44.009000;  // CO2
  CW += c[10] * 29.018000; // HCO
  CW += c[11] * 28.014000; // N2
  // Now compute conversion
  amrex::Real CWinv = 1.0 / CW;
  y[0] = c[0] * 2.016000 * CWinv;
  y[1] = c[1] * 31.998000 * CWinv;
  y[2] = c[2] * 15.999000 * CWinv;
  y[3] = c[3] * 17.007000 * CWinv;
  y[4] = c[4] * 18.015000 * CWinv;
  y[5] = c[5] * 1.008000 * CWinv;
  y[6] = c[6] * 33.006000 * CWinv;
  y[7] = c[7] * 34.014000 * CWinv;
  y[8] = c[8] * 28.010000 * CWinv;
  y[9] = c[9] * 44.009000 * CWinv;
  y[10] = c[10] * 29.018000 * CWinv;
  y[11] = c[11] * 28.014000 * CWinv;
}

// get specific heat at constant volume as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVML(const amrex::Real* T, amrex::Real* cvml)
{
  int id;              // loop counter
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cv_R(cvml, tc);

  // convert to chemkin units
  for (id = 0; id < 12; ++id) {
    cvml[id] *= 8.31446261815324e+07;
  }
}

// get specific heat at constant pressure as a
// function of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPML(const amrex::Real* T, amrex::Real* cpml)
{
  int id;              // loop counter
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cp_R(cpml, tc);

  // convert to chemkin units
  for (id = 0; id < 12; ++id) {
    cpml[id] *= 8.31446261815324e+07;
  }
}

// get internal energy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUML(const amrex::Real* T, amrex::Real* uml)
{
  int id;              // loop counter
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesInternalEnergy(uml, tc);

  // convert to chemkin units
  for (id = 0; id < 12; ++id) {
    uml[id] *= RT;
  }
}

// get enthalpy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHML(const amrex::Real* T, amrex::Real* hml)
{
  int id;              // loop counter
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesEnthalpy(hml, tc);

  // convert to chemkin units
  for (id = 0; id < 12; ++id) {
    hml[id] *= RT;
  }
}

// Returns the standard-state entropies in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSML(const amrex::Real* T, amrex::Real* sml)
{
  int id;              // loop counter
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  speciesEntropy(sml, tc);

  // convert to chemkin units
  for (id = 0; id < 12; ++id) {
    sml[id] *= 8.31446261815324e+07;
  }
}

// Returns the specific heats at constant volume
// in mass units (Eq. 29)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVMS(const amrex::Real* T, amrex::Real* cvms)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cv_R(cvms, tc);
  // multiply by R/molecularweight
  cvms[0] *= 4.124237409798234e+07;  // H2
  cvms[1] *= 2.598431970171023e+06;  // O2
  cvms[2] *= 5.196863940342046e+06;  // O
  cvms[3] *= 4.888847308845322e+06;  // OH
  cvms[4] *= 4.615299815794193e+06;  // H2O
  cvms[5] *= 8.248474819596468e+07;  // H
  cvms[6] *= 2.519076112874398e+06;  // HO2
  cvms[7] *= 2.444423654422661e+06;  // H2O2
  cvms[8] *= 2.968390795484913e+06;  // CO
  cvms[9] *= 1.889264154639560e+06;  // CO2
  cvms[10] *= 2.865277627042952e+06; // HCO
  cvms[11] *= 2.967966951578939e+06; // N2
}

// Returns the specific heats at constant pressure
// in mass units (Eq. 26)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPMS(const amrex::Real* T, amrex::Real* cpms)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  cp_R(cpms, tc);
  // multiply by R/molecularweight
  cpms[0] *= 4.124237409798234e+07;  // H2
  cpms[1] *= 2.598431970171023e+06;  // O2
  cpms[2] *= 5.196863940342046e+06;  // O
  cpms[3] *= 4.888847308845322e+06;  // OH
  cpms[4] *= 4.615299815794193e+06;  // H2O
  cpms[5] *= 8.248474819596468e+07;  // H
  cpms[6] *= 2.519076112874398e+06;  // HO2
  cpms[7] *= 2.444423654422661e+06;  // H2O2
  cpms[8] *= 2.968390795484913e+06;  // CO
  cpms[9] *= 1.889264154639560e+06;  // CO2
  cpms[10] *= 2.865277627042952e+06; // HCO
  cpms[11] *= 2.967966951578939e+06; // N2
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUMS(const amrex::Real* T, amrex::Real* ums)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[12];

  get_imw(imw);

  speciesInternalEnergy(ums, tc);

  for (int i = 0; i < 12; i++) {
    ums[i] *= RT * imw[i];
  }
}

// Returns enthalpy in mass units (Eq 27.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHMS(const amrex::Real* T, amrex::Real* hms)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[12];

  get_imw(imw);

  speciesEnthalpy(hms, tc);

  for (int i = 0; i < 12; i++) {
    hms[i] *= RT * imw[i];
  }
}

// Returns the entropies in mass units (Eq 28.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSMS(const amrex::Real* T, amrex::Real* sms)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    log(tT), tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  speciesEntropy(sms, tc);
  // multiply by R/molecularweight
  sms[0] *= 4.124237409798234e+07;  // H2
  sms[1] *= 2.598431970171023e+06;  // O2
  sms[2] *= 5.196863940342046e+06;  // O
  sms[3] *= 4.888847308845322e+06;  // OH
  sms[4] *= 4.615299815794193e+06;  // H2O
  sms[5] *= 8.248474819596468e+07;  // H
  sms[6] *= 2.519076112874398e+06;  // HO2
  sms[7] *= 2.444423654422661e+06;  // H2O2
  sms[8] *= 2.968390795484913e+06;  // CO
  sms[9] *= 1.889264154639560e+06;  // CO2
  sms[10] *= 2.865277627042952e+06; // HCO
  sms[11] *= 2.967966951578939e+06; // N2
}

// GPU version of productionRate: no more use of thermo namespace vectors
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_qfqr(
  amrex::Real* qf,
  amrex::Real* qr,
  const amrex::Real* sc,
  const amrex::Real* /*sc_qss*/,
  const amrex::Real* tc,
  const amrex::Real invT)
{

  // reaction 8: H + O2 (+M) <=> HO2 (+M)
  qf[0] = sc[1] * sc[5];
  qr[0] = sc[6];

  // reaction 15: H2O2 (+M) <=> 2 OH (+M)
  qf[1] = sc[7];
  qr[1] = (sc[3] * sc[3]);

  // reaction 21: CO + O (+M) <=> CO2 (+M)
  qf[2] = sc[2] * sc[8];
  qr[2] = sc[9];

  // reaction 4: H2 + M <=> 2 H + M
  qf[3] = sc[0];
  qr[3] = (sc[5] * sc[5]);

  // reaction 5: 2 O + M <=> O2 + M
  qf[4] = (sc[2] * sc[2]);
  qr[4] = sc[1];

  // reaction 6: H + O + M <=> OH + M
  qf[5] = sc[2] * sc[5];
  qr[5] = sc[3];

  // reaction 7: H + OH + M <=> H2O + M
  qf[6] = sc[3] * sc[5];
  qr[6] = sc[4];

  // reaction 25: HCO + M <=> CO + H + M
  qf[7] = sc[10];
  qr[7] = sc[5] * sc[8];

  // reaction 0: H + O2 <=> O + OH
  qf[8] = sc[1] * sc[5];
  qr[8] = sc[2] * sc[3];

  // reaction 1: H2 + O <=> H + OH
  qf[9] = sc[0] * sc[2];
  qr[9] = sc[3] * sc[5];

  // reaction 2: H2 + OH <=> H + H2O
  qf[10] = sc[0] * sc[3];
  qr[10] = sc[4] * sc[5];

  // reaction 3: H2O + O <=> 2 OH
  qf[11] = sc[2] * sc[4];
  qr[11] = (sc[3] * sc[3]);

  // reaction 9: H + HO2 <=> H2 + O2
  qf[12] = sc[5] * sc[6];
  qr[12] = sc[0] * sc[1];

  // reaction 10: H + HO2 <=> 2 OH
  qf[13] = sc[5] * sc[6];
  qr[13] = (sc[3] * sc[3]);

  // reaction 11: HO2 + O <=> O2 + OH
  qf[14] = sc[2] * sc[6];
  qr[14] = sc[1] * sc[3];

  // reaction 12: HO2 + OH <=> H2O + O2
  qf[15] = sc[3] * sc[6];
  qr[15] = sc[1] * sc[4];

  // reaction 13: 2 HO2 <=> H2O2 + O2
  qf[16] = (sc[6] * sc[6]);
  qr[16] = sc[1] * sc[7];

  // reaction 14: 2 HO2 <=> H2O2 + O2
  qf[17] = (sc[6] * sc[6]);
  qr[17] = sc[1] * sc[7];

  // reaction 16: H + H2O2 <=> H2O + OH
  qf[18] = sc[5] * sc[7];
  qr[18] = sc[3] * sc[4];

  // reaction 17: H + H2O2 <=> H2 + HO2
  qf[19] = sc[5] * sc[7];
  qr[19] = sc[0] * sc[6];

  // reaction 18: H2O2 + O <=> HO2 + OH
  qf[20] = sc[2] * sc[7];
  qr[20] = sc[3] * sc[6];

  // reaction 19: H2O2 + OH <=> H2O + HO2
  qf[21] = sc[3] * sc[7];
  qr[21] = sc[4] * sc[6];

  // reaction 20: H2O2 + OH <=> H2O + HO2
  qf[22] = sc[3] * sc[7];
  qr[22] = sc[4] * sc[6];

  // reaction 22: CO + O2 <=> CO2 + O
  qf[23] = sc[1] * sc[8];
  qr[23] = sc[2] * sc[9];

  // reaction 23: CO + HO2 <=> CO2 + OH
  qf[24] = sc[6] * sc[8];
  qr[24] = sc[3] * sc[9];

  // reaction 24: CO + OH <=> CO2 + H
  qf[25] = sc[3] * sc[8];
  qr[25] = sc[5] * sc[9];

  // reaction 26: HCO + O2 <=> CO + HO2
  qf[26] = sc[1] * sc[10];
  qr[26] = sc[6] * sc[8];

  // reaction 27: H + HCO <=> CO + H2
  qf[27] = sc[5] * sc[10];
  qr[27] = sc[0] * sc[8];

  // reaction 28: HCO + O <=> CO2 + H
  qf[28] = sc[2] * sc[10];
  qr[28] = sc[5] * sc[9];

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 12; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[12];
  gibbs(g_RT, tc);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 * invT;
  amrex::Real refCinv = 1 / refC;

  // Evaluate the kfs
  amrex::Real k_f, Corr;
  amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

  // reaction 0:  H + O2 <=> O + OH
  k_f = 3547000000 * exp((-0.406) * tc[0] - (8352.89343569254) * invT);
  qf[8] *= k_f;
  qr[8] *= k_f * exp(-(g_RT[1] - g_RT[2] - g_RT[3] + g_RT[5]));
  // reaction 1:  H2 + O <=> H + OH
  k_f = 0.0508 * exp((2.67) * tc[0] - (3165.23282791169) * invT);
  qf[9] *= k_f;
  qr[9] *= k_f * exp(-(g_RT[0] + g_RT[2] - g_RT[3] - g_RT[5]));
  // reaction 2:  H2 + OH <=> H + H2O
  k_f = 216 * exp((1.51) * tc[0] - (1726.03316371019) * invT);
  qf[10] *= k_f;
  qr[10] *= k_f * exp(-(g_RT[0] + g_RT[3] - g_RT[4] - g_RT[5]));
  // reaction 3:  H2O + O <=> 2 OH
  k_f = 2.97 * exp((2.02) * tc[0] - (6743.10332178324) * invT);
  qf[11] *= k_f;
  qr[11] *= k_f * exp(-(g_RT[2] - 2.000000 * g_RT[3] + g_RT[4]));
  // reaction 4:  H2 + M <=> 2 H + M
  k_f = 45770000000000 * exp((-1.4) * tc[0] - (52525.7555766966) * invT);
  Corr = mixture + (1.5) * sc[0] + (11) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
  qf[3] *= Corr * k_f;
  qr[3] *= Corr * k_f * exp(-(g_RT[0] - 2.000000 * g_RT[5])) * (refCinv);
  // reaction 5:  2 O + M <=> O2 + M
  k_f = 6165 * exp((-0.5) * tc[0]);
  Corr = mixture + (1.5) * sc[0] + (11) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
  qf[4] *= Corr * k_f;
  qr[4] *= Corr * k_f * exp(-(-g_RT[1] + 2.000000 * g_RT[2])) * (refC);
  // reaction 6:  H + O + M <=> OH + M
  k_f = 4714000 * exp((-1) * tc[0]);
  Corr = mixture + (1.5) * sc[0] + (11) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
  qf[5] *= Corr * k_f;
  qr[5] *= Corr * k_f * exp(-(g_RT[2] - g_RT[3] + g_RT[5])) * (refC);
  // reaction 7:  H + OH + M <=> H2O + M
  k_f = 38000000000 * exp((-2) * tc[0]);
  Corr = mixture + (1.5) * sc[0] + (11) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
  qf[6] *= Corr * k_f;
  qr[6] *= Corr * k_f * exp(-(g_RT[3] - g_RT[4] + g_RT[5])) * (refC);
  // reaction 8:  H + O2 (+M) <=> HO2 (+M)
  k_f = 1475000 * exp((0.6) * tc[0]);
  Corr = mixture + sc[0] + (-0.22) * sc[1] + (10) * sc[4] + (0.9) * sc[8] +
         (2.8) * sc[9];
  redP =
    Corr / k_f * 636600000 * exp(-1.72 * tc[0] - (264.088106214317) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent =
    log10((0.2) * exp(-tc[1] * 1e+30) + 0.8 * exp(-tc[1] * 1e-30) + 0.0);
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[0] *= Corr * k_f;
  qr[0] *= Corr * k_f * exp(-(g_RT[1] + g_RT[5] - g_RT[6])) * (refC);
  // reaction 9:  H + HO2 <=> H2 + O2
  k_f = 16600000 * exp(-(414.147315957284) * invT);
  qf[12] *= k_f;
  qr[12] *= k_f * exp(-(-g_RT[0] - g_RT[1] + g_RT[5] + g_RT[6]));
  // reaction 10:  H + HO2 <=> 2 OH
  k_f = 70790000 * exp(-(148.448916412392) * invT);
  qf[13] *= k_f;
  qr[13] *= k_f * exp(-(-2.000000 * g_RT[3] + g_RT[5] + g_RT[6]));
  // reaction 11:  HO2 + O <=> O2 + OH
  k_f = 32500000;
  qf[14] *= k_f;
  qr[14] *= k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[6]));
  // reaction 12:  HO2 + OH <=> H2O + O2
  k_f = 28900000 * exp(-(-250.098682904946) * invT);
  qf[15] *= k_f;
  qr[15] *= k_f * exp(-(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[6]));
  // reaction 13:  2 HO2 <=> H2O2 + O2
  k_f = 420000000 * exp(-(6029.54208967215) * invT);
  qf[16] *= k_f;
  qr[16] *= k_f * exp(-(-g_RT[1] + 2.000000 * g_RT[6] - g_RT[7]));
  // reaction 14:  2 HO2 <=> H2O2 + O2
  k_f = 130000 * exp(-(-819.89091359563) * invT);
  qf[17] *= k_f;
  qr[17] *= k_f * exp(-(-g_RT[1] + 2.000000 * g_RT[6] - g_RT[7]));
  // reaction 15:  H2O2 (+M) <=> 2 OH (+M)
  k_f = 295100000000000 * exp(-(24370.7831249226) * invT);
  Corr = mixture + (1.5) * sc[0] + (11) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
  redP = Corr / k_f * 120200000000 * exp(-(22896.3582941147) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent =
    log10((0.5) * exp(-tc[1] * 1e+30) + 0.5 * exp(-tc[1] * 1e-30) + 0.0);
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[1] *= Corr * k_f;
  qr[1] *= Corr * k_f * exp(-(-2.000000 * g_RT[3] + g_RT[7])) * (refCinv);
  // reaction 16:  H + H2O2 <=> H2O + OH
  k_f = 24100000 * exp(-(1997.77016324474) * invT);
  qf[18] *= k_f;
  qr[18] *= k_f * exp(-(-g_RT[3] - g_RT[4] + g_RT[5] + g_RT[7]));
  // reaction 17:  H + H2O2 <=> H2 + HO2
  k_f = 48200000 * exp(-(4000.57249314752) * invT);
  qf[19] *= k_f;
  qr[19] *= k_f * exp(-(-g_RT[0] + g_RT[5] - g_RT[6] + g_RT[7]));
  // reaction 18:  H2O2 + O <=> HO2 + OH
  k_f = 9.55 * exp((2) * tc[0] - (1997.77016324474) * invT);
  qf[20] *= k_f;
  qr[20] *= k_f * exp(-(g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7]));
  // reaction 19:  H2O2 + OH <=> H2O + HO2
  k_f = 1000000;
  qf[21] *= k_f;
  qr[21] *= k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]));
  // reaction 20:  H2O2 + OH <=> H2O + HO2
  k_f = 580000000 * exp(-(4809.24167509571) * invT);
  qf[22] *= k_f;
  qr[22] *= k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]));
  // reaction 21:  CO + O (+M) <=> CO2 (+M)
  k_f = 18000 * exp(-(1199.66853127845) * invT);
  Corr = mixture + (1.5) * sc[0] + (11) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
  redP =
    Corr / k_f * 1550000000000 * exp(-2.79 * tc[0] - (2108.98104638758) * invT);
  Corr = redP / (1. + redP);
  qf[2] *= Corr * k_f;
  qr[2] *= Corr * k_f * exp(-(g_RT[2] + g_RT[8] - g_RT[9])) * (refC);
  // reaction 22:  CO + O2 <=> CO2 + O
  k_f = 2530000 * exp(-(24003.4349588851) * invT);
  qf[23] *= k_f;
  qr[23] *= k_f * exp(-(g_RT[1] - g_RT[2] + g_RT[8] - g_RT[9]));
  // reaction 23:  CO + HO2 <=> CO2 + OH
  k_f = 30100000 * exp(-(11573.9833135086) * invT);
  qf[24] *= k_f;
  qr[24] *= k_f * exp(-(-g_RT[3] + g_RT[6] + g_RT[8] - g_RT[9]));
  // reaction 24:  CO + OH <=> CO2 + H
  k_f = 0.2229 * exp((1.89) * tc[0] - (-583.077150667929) * invT);
  qf[25] *= k_f;
  qr[25] *= k_f * exp(-(g_RT[3] - g_RT[5] + g_RT[8] - g_RT[9]));
  // reaction 25:  HCO + M <=> CO + H + M
  k_f = 474850 * exp((0.659) * tc[0] - (7484.8446871794) * invT);
  Corr = mixture + (1.5) * sc[0] + (5) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
  qf[7] *= Corr * k_f;
  qr[7] *= Corr * k_f * exp(-(-g_RT[5] - g_RT[8] + g_RT[10])) * (refCinv);
  // reaction 26:  HCO + O2 <=> CO + HO2
  k_f = 7580000 * exp(-(206.318832979935) * invT);
  qf[26] *= k_f;
  qr[26] *= k_f * exp(-(g_RT[1] - g_RT[6] - g_RT[8] + g_RT[10]));
  // reaction 27:  H + HCO <=> CO + H2
  k_f = 72300000;
  qf[27] *= k_f;
  qr[27] *= k_f * exp(-(-g_RT[0] + g_RT[5] - g_RT[8] + g_RT[10]));
  // reaction 28:  HCO + O <=> CO2 + H
  k_f = 30000000;
  qf[28] *= k_f;
  qr[28] *= k_f * exp(-(g_RT[2] - g_RT[5] - g_RT[9] + g_RT[10]));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
productionRate(amrex::Real* wdot, const amrex::Real* sc, const amrex::Real T)
{
  const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
  const amrex::Real invT = 1.0 / tc[1];

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  const amrex::Real refC = 101325 / 8.31446 * invT;
  const amrex::Real refCinv = 1 / refC;

  for (int i = 0; i < 12; ++i) {
    wdot[i] = 0.0;
  }

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 12; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[12];
  gibbs(g_RT, tc);

  {
    // reaction 8:  H + O2 (+M) <=> HO2 (+M)
    const amrex::Real k_f = 1475000 * exp((0.6) * tc[0]);
    amrex::Real Corr = mixture + sc[0] + (-0.22) * sc[1] + (10) * sc[4] +
                       (0.9) * sc[8] + (2.8) * sc[9];
    const amrex::Real redP =
      Corr / k_f * 636600000 * exp(-1.72 * tc[0] - 264.088106214317 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent =
      log10(0.2 * exp(-tc[1] * 1e+30) + 0.8 * exp(-tc[1] * 1e-30) + 0.0);
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[1] * sc[5]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] + g_RT[5] - g_RT[6])) * (refC) * (sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
  }

  {
    // reaction 15:  H2O2 (+M) <=> 2 OH (+M)
    const amrex::Real k_f = 295100000000000 * exp(-(24370.7831249226) * invT);
    amrex::Real Corr =
      mixture + (1.5) * sc[0] + (11) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
    const amrex::Real redP =
      Corr / k_f * 120200000000 * exp(-22896.3582941147 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent =
      log10(0.5 * exp(-tc[1] * 1e+30) + 0.5 * exp(-tc[1] * 1e-30) + 0.0);
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[7]);
    const amrex::Real qr = Corr * k_f * exp(-(-2.000000 * g_RT[3] + g_RT[7])) *
                           (refCinv) * ((sc[3] * sc[3]));
    const amrex::Real qdot = qf - qr;
    wdot[3] += 2.000000 * qdot;
    wdot[7] -= qdot;
  }

  {
    // reaction 21:  CO + O (+M) <=> CO2 (+M)
    const amrex::Real k_f = 18000 * exp(-(1199.66853127845) * invT);
    amrex::Real Corr =
      mixture + (1.5) * sc[0] + (11) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
    const amrex::Real redP =
      Corr / k_f * 1550000000000 * exp(-2.79 * tc[0] - 2108.98104638758 * invT);
    Corr = redP / (1.0 + redP);
    const amrex::Real qf = Corr * k_f * (sc[2] * sc[8]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[2] + g_RT[8] - g_RT[9])) * (refC) * (sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 4:  H2 + M <=> 2 H + M
    const amrex::Real k_f =
      45770000000000 * exp((-1.4) * tc[0] - (52525.7555766966) * invT);
    const amrex::Real Corr =
      mixture + (1.5) * sc[0] + (11) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
    const amrex::Real qf = Corr * k_f * (sc[0]);
    const amrex::Real qr = Corr * k_f * exp(-(g_RT[0] - 2.000000 * g_RT[5])) *
                           (refCinv) * ((sc[5] * sc[5]));
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[5] += 2.000000 * qdot;
  }

  {
    // reaction 5:  2 O + M <=> O2 + M
    const amrex::Real k_f = 6165 * exp((-0.5) * tc[0]);
    const amrex::Real Corr =
      mixture + (1.5) * sc[0] + (11) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
    const amrex::Real qf = Corr * k_f * ((sc[2] * sc[2]));
    const amrex::Real qr =
      Corr * k_f * exp(-(-g_RT[1] + 2.000000 * g_RT[2])) * (refC) * (sc[1]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= 2.000000 * qdot;
  }

  {
    // reaction 6:  H + O + M <=> OH + M
    const amrex::Real k_f = 4714000 * exp((-1) * tc[0]);
    const amrex::Real Corr =
      mixture + (1.5) * sc[0] + (11) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
    const amrex::Real qf = Corr * k_f * (sc[2] * sc[5]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[2] - g_RT[3] + g_RT[5])) * (refC) * (sc[3]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;
  }

  {
    // reaction 7:  H + OH + M <=> H2O + M
    const amrex::Real k_f = 38000000000 * exp((-2) * tc[0]);
    const amrex::Real Corr =
      mixture + (1.5) * sc[0] + (11) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
    const amrex::Real qf = Corr * k_f * (sc[3] * sc[5]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[3] - g_RT[4] + g_RT[5])) * (refC) * (sc[4]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;
  }

  {
    // reaction 25:  HCO + M <=> CO + H + M
    const amrex::Real k_f =
      474850 * exp((0.659) * tc[0] - (7484.8446871794) * invT);
    const amrex::Real Corr =
      mixture + (1.5) * sc[0] + (5) * sc[4] + (0.9) * sc[8] + (2.8) * sc[9];
    const amrex::Real qf = Corr * k_f * (sc[10]);
    const amrex::Real qr = Corr * k_f * exp(-(-g_RT[5] - g_RT[8] + g_RT[10])) *
                           (refCinv) * (sc[5] * sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;
  }

  {
    // reaction 0:  H + O2 <=> O + OH
    const amrex::Real k_f =
      3547000000 * exp((-0.406) * tc[0] - (8352.89343569254) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[5]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[2] - g_RT[3] + g_RT[5])) * (sc[2] * sc[3]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;
  }

  {
    // reaction 1:  H2 + O <=> H + OH
    const amrex::Real k_f =
      0.0508 * exp((2.67) * tc[0] - (3165.23282791169) * invT);
    const amrex::Real qf = k_f * (sc[0] * sc[2]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] + g_RT[2] - g_RT[3] - g_RT[5])) * (sc[3] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 2:  H2 + OH <=> H + H2O
    const amrex::Real k_f =
      216 * exp((1.51) * tc[0] - (1726.03316371019) * invT);
    const amrex::Real qf = k_f * (sc[0] * sc[3]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] + g_RT[3] - g_RT[4] - g_RT[5])) * (sc[4] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 3:  H2O + O <=> 2 OH
    const amrex::Real k_f =
      2.97 * exp((2.02) * tc[0] - (6743.10332178324) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[4]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - 2.000000 * g_RT[3] + g_RT[4])) * ((sc[3] * sc[3]));
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[3] += 2.000000 * qdot;
    wdot[4] -= qdot;
  }

  {
    // reaction 9:  H + HO2 <=> H2 + O2
    const amrex::Real k_f = 16600000 * exp(-(414.147315957284) * invT);
    const amrex::Real qf = k_f * (sc[5] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] - g_RT[1] + g_RT[5] + g_RT[6])) * (sc[0] * sc[1]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 10:  H + HO2 <=> 2 OH
    const amrex::Real k_f = 70790000 * exp(-(148.448916412392) * invT);
    const amrex::Real qf = k_f * (sc[5] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(-2.000000 * g_RT[3] + g_RT[5] + g_RT[6])) * ((sc[3] * sc[3]));
    const amrex::Real qdot = qf - qr;
    wdot[3] += 2.000000 * qdot;
    wdot[5] -= qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 11:  HO2 + O <=> O2 + OH
    const amrex::Real k_f = 32500000;
    const amrex::Real qf = k_f * (sc[2] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[6])) * (sc[1] * sc[3]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 12:  HO2 + OH <=> H2O + O2
    const amrex::Real k_f = 28900000 * exp(-(-250.098682904946) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[6])) * (sc[1] * sc[4]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 13:  2 HO2 <=> H2O2 + O2
    const amrex::Real k_f = 420000000 * exp(-(6029.54208967215) * invT);
    const amrex::Real qf = k_f * ((sc[6] * sc[6]));
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + 2.000000 * g_RT[6] - g_RT[7])) * (sc[1] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[6] -= 2.000000 * qdot;
    wdot[7] += qdot;
  }

  {
    // reaction 14:  2 HO2 <=> H2O2 + O2
    const amrex::Real k_f = 130000 * exp(-(-819.89091359563) * invT);
    const amrex::Real qf = k_f * ((sc[6] * sc[6]));
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + 2.000000 * g_RT[6] - g_RT[7])) * (sc[1] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[6] -= 2.000000 * qdot;
    wdot[7] += qdot;
  }

  {
    // reaction 16:  H + H2O2 <=> H2O + OH
    const amrex::Real k_f = 24100000 * exp(-(1997.77016324474) * invT);
    const amrex::Real qf = k_f * (sc[5] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[3] - g_RT[4] + g_RT[5] + g_RT[7])) * (sc[3] * sc[4]);
    const amrex::Real qdot = qf - qr;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[7] -= qdot;
  }

  {
    // reaction 17:  H + H2O2 <=> H2 + HO2
    const amrex::Real k_f = 48200000 * exp(-(4000.57249314752) * invT);
    const amrex::Real qf = k_f * (sc[5] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] + g_RT[5] - g_RT[6] + g_RT[7])) * (sc[0] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;
  }

  {
    // reaction 18:  H2O2 + O <=> HO2 + OH
    const amrex::Real k_f = 9.55 * exp((2) * tc[0] - (1997.77016324474) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7])) * (sc[3] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;
  }

  {
    // reaction 19:  H2O2 + OH <=> H2O + HO2
    const amrex::Real k_f = 1000000;
    const amrex::Real qf = k_f * (sc[3] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7])) * (sc[4] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;
  }

  {
    // reaction 20:  H2O2 + OH <=> H2O + HO2
    const amrex::Real k_f = 580000000 * exp(-(4809.24167509571) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7])) * (sc[4] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;
  }

  {
    // reaction 22:  CO + O2 <=> CO2 + O
    const amrex::Real k_f = 2530000 * exp(-(24003.4349588851) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[2] + g_RT[8] - g_RT[9])) * (sc[2] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 23:  CO + HO2 <=> CO2 + OH
    const amrex::Real k_f = 30100000 * exp(-(11573.9833135086) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[3] + g_RT[6] + g_RT[8] - g_RT[9])) * (sc[3] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 24:  CO + OH <=> CO2 + H
    const amrex::Real k_f =
      0.2229 * exp((1.89) * tc[0] - (-583.077150667929) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[5] + g_RT[8] - g_RT[9])) * (sc[5] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 26:  HCO + O2 <=> CO + HO2
    const amrex::Real k_f = 7580000 * exp(-(206.318832979935) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[10]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[6] - g_RT[8] + g_RT[10])) * (sc[6] * sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;
  }

  {
    // reaction 27:  H + HCO <=> CO + H2
    const amrex::Real k_f = 72300000;
    const amrex::Real qf = k_f * (sc[5] * sc[10]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] + g_RT[5] - g_RT[8] + g_RT[10])) * (sc[0] * sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;
  }

  {
    // reaction 28:  HCO + O <=> CO2 + H
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[2] * sc[10]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - g_RT[5] - g_RT[9] + g_RT[10])) * (sc[5] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;
  }
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(const amrex::Real* T, amrex::Real* C, amrex::Real* wdot)
{
  int id; // loop counter

  // convert to SI
  for (id = 0; id < 12; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, *T);

  // convert to chemkin units
  for (id = 0; id < 12; ++id) {
    C[id] *= 1.0e-6;
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWYP(
  const amrex::Real* P,
  const amrex::Real* T,
  const amrex::Real* y,
  amrex::Real* wdot)
{
  int id;            // loop counter
  amrex::Real c[12]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;
  amrex::Real imw[12];

  get_imw(imw);

  // Compute inverse of mean molecular wt first
  YOW += y[0] * imw[0];   // H2
  YOW += y[1] * imw[1];   // O2
  YOW += y[2] * imw[2];   // O
  YOW += y[3] * imw[3];   // OH
  YOW += y[4] * imw[4];   // H2O
  YOW += y[5] * imw[5];   // H
  YOW += y[6] * imw[6];   // HO2
  YOW += y[7] * imw[7];   // H2O2
  YOW += y[8] * imw[8];   // CO
  YOW += y[9] * imw[9];   // CO2
  YOW += y[10] * imw[10]; // HCO
  YOW += y[11] * imw[11]; // N2
  // PW/RT (see Eq. 7)
  PWORT = (*P) / (YOW * 8.31446261815324e+07 * (*T));
  // multiply by 1e6 so c goes to SI
  PWORT *= 1e6;
  // Now compute conversion (and go to SI)
  c[0] = PWORT * y[0] * imw[0];
  c[1] = PWORT * y[1] * imw[1];
  c[2] = PWORT * y[2] * imw[2];
  c[3] = PWORT * y[3] * imw[3];
  c[4] = PWORT * y[4] * imw[4];
  c[5] = PWORT * y[5] * imw[5];
  c[6] = PWORT * y[6] * imw[6];
  c[7] = PWORT * y[7] * imw[7];
  c[8] = PWORT * y[8] * imw[8];
  c[9] = PWORT * y[9] * imw[9];
  c[10] = PWORT * y[10] * imw[10];
  c[11] = PWORT * y[11] * imw[11];

  // convert to chemkin units
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 12; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWXP(
  const amrex::Real* P,
  const amrex::Real* T,
  const amrex::Real* x,
  amrex::Real* wdot)
{
  int id;            // loop counter
  amrex::Real c[12]; // temporary storage
  amrex::Real PORT =
    1e6 * (*P) /
    (8.31446261815324e+07 * (*T)); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (id = 0; id < 12; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 12; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWYR(
  const amrex::Real* rho,
  const amrex::Real* T,
  const amrex::Real* y,
  amrex::Real* wdot)
{
  int id;            // loop counter
  amrex::Real c[12]; // temporary storage
  amrex::Real imw[12];

  get_imw(imw);

  // See Eq 8 with an extra 1e6 so c goes to SI
  c[0] = 1e6 * (*rho) * y[0] * imw[0];
  c[1] = 1e6 * (*rho) * y[1] * imw[1];
  c[2] = 1e6 * (*rho) * y[2] * imw[2];
  c[3] = 1e6 * (*rho) * y[3] * imw[3];
  c[4] = 1e6 * (*rho) * y[4] * imw[4];
  c[5] = 1e6 * (*rho) * y[5] * imw[5];
  c[6] = 1e6 * (*rho) * y[6] * imw[6];
  c[7] = 1e6 * (*rho) * y[7] * imw[7];
  c[8] = 1e6 * (*rho) * y[8] * imw[8];
  c[9] = 1e6 * (*rho) * y[9] * imw[9];
  c[10] = 1e6 * (*rho) * y[10] * imw[10];
  c[11] = 1e6 * (*rho) * y[11] * imw[11];

  // call productionRate
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 12; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWXR(
  const amrex::Real* rho,
  const amrex::Real* T,
  const amrex::Real* x,
  amrex::Real* wdot)
{
  int id;             // loop counter
  amrex::Real c[12];  // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 2.016000;   // H2
  XW += x[1] * 31.998000;  // O2
  XW += x[2] * 15.999000;  // O
  XW += x[3] * 17.007000;  // OH
  XW += x[4] * 18.015000;  // H2O
  XW += x[5] * 1.008000;   // H
  XW += x[6] * 33.006000;  // HO2
  XW += x[7] * 34.014000;  // H2O2
  XW += x[8] * 28.010000;  // CO
  XW += x[9] * 44.009000;  // CO2
  XW += x[10] * 29.018000; // HCO
  XW += x[11] * 28.014000; // N2
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * (*rho) / XW;

  // Compute conversion, see Eq 11
  for (id = 0; id < 12; ++id) {
    c[id] = x[id] * ROW;
  }

  // convert to chemkin units
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 12; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int* kcharge)
{
  kcharge[0] = 0;  // H2
  kcharge[1] = 0;  // O2
  kcharge[2] = 0;  // O
  kcharge[3] = 0;  // OH
  kcharge[4] = 0;  // H2O
  kcharge[5] = 0;  // H
  kcharge[6] = 0;  // HO2
  kcharge[7] = 0;  // H2O2
  kcharge[8] = 0;  // CO
  kcharge[9] = 0;  // CO2
  kcharge[10] = 0; // HCO
  kcharge[11] = 0; // N2
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real* zk)
{
  amrex::Real imw[12];
  get_imw(imw);

  int kchrg[12];
  CKCHRG(kchrg);

  for (int id = 0; id < 12; ++id) {
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
    // species 0: H2
    species[0] = +8.24944200e-04 - 1.62860300e-06 * tc[1] -
                 2.84263020e-10 * tc[2] + 1.65394880e-12 * tc[3];
    // species 1: O2
    species[1] = +1.12748600e-03 - 1.15123000e-06 * tc[1] +
                 3.94163100e-09 * tc[2] - 3.50742160e-12 * tc[3];
    // species 2: O
    species[2] = -1.63816600e-03 + 4.84206400e-06 * tc[1] -
                 4.80852900e-09 * tc[2] + 1.55627840e-12 * tc[3];
    // species 3: OH
    species[3] = -3.22544939e-03 + 1.30552938e-05 * tc[1] -
                 1.73956093e-08 * tc[2] + 8.24949516e-12 * tc[3];
    // species 4: H2O
    species[4] = +3.47498200e-03 - 1.27093920e-05 * tc[1] +
                 2.09057430e-08 * tc[2] - 1.00263520e-11 * tc[3];
    // species 6: HO2
    species[6] = -4.74912051e-03 + 4.23165782e-05 * tc[1] -
                 7.28291682e-08 * tc[2] + 3.71690050e-11 * tc[3];
    // species 7: H2O2
    species[7] = +6.56922600e-03 - 2.97002600e-07 * tc[1] -
                 1.38774180e-08 * tc[2] + 9.88606000e-12 * tc[3];
    // species 8: CO
    species[8] = +1.51194100e-03 - 7.76351000e-06 * tc[1] +
                 1.67458320e-08 * tc[2] - 9.89980400e-12 * tc[3];
    // species 9: CO2
    species[9] = +9.92207200e-03 - 2.08182200e-05 * tc[1] +
                 2.06000610e-08 * tc[2] - 8.46912000e-12 * tc[3];
    // species 10: HCO
    species[10] = +6.19914700e-03 - 1.92461680e-05 * tc[1] +
                  3.26947500e-08 * tc[2] - 1.82995400e-11 * tc[3];
    // species 11: N2
    species[11] = +1.40824000e-03 - 7.92644400e-06 * tc[1] +
                  1.69245450e-08 * tc[2] - 9.77942000e-12 * tc[3];
  } else {
    // species 0: H2
    species[0] = +7.00064400e-04 - 1.12676580e-07 * tc[1] -
                 2.76947340e-11 * tc[2] + 6.33100800e-15 * tc[3];
    // species 1: O2
    species[1] = +6.13519700e-04 - 2.51768400e-07 * tc[1] +
                 5.32584300e-11 * tc[2] - 4.54574000e-15 * tc[3];
    // species 2: O
    species[2] = -2.75506200e-05 - 6.20560600e-09 * tc[1] +
                 1.36532010e-11 * tc[2] - 1.74722080e-15 * tc[3];
    // species 3: OH
    species[3] = +1.05650448e-03 - 5.18165516e-07 * tc[1] +
                 9.15656022e-11 * tc[2] - 5.32783504e-15 * tc[3];
    // species 4: H2O
    species[4] = +3.05629300e-03 - 1.74605200e-06 * tc[1] +
                 3.60298800e-10 * tc[2] - 2.55664720e-14 * tc[3];
    // species 6: HO2
    species[6] = +2.23982013e-03 - 1.26731630e-06 * tc[1] +
                 3.42739110e-10 * tc[2] - 4.31634140e-14 * tc[3];
    // species 7: H2O2
    species[7] = +4.33613600e-03 - 2.94937800e-06 * tc[1] +
                 7.04671200e-10 * tc[2] - 5.72661600e-14 * tc[3];
    // species 8: CO
    species[8] = +1.44268900e-03 - 1.12616560e-06 * tc[1] +
                 3.05574300e-10 * tc[2] - 2.76438080e-14 * tc[3];
    // species 9: CO2
    species[9] = +3.14016900e-03 - 2.55682200e-06 * tc[1] +
                 7.18199100e-10 * tc[2] - 6.67613200e-14 * tc[3];
    // species 10: HCO
    species[10] = +3.34557300e-03 - 2.67001200e-06 * tc[1] +
                  7.41171900e-10 * tc[2] - 6.85540400e-14 * tc[3];
    // species 11: N2
    species[11] = +1.48797700e-03 - 1.13695220e-06 * tc[1] +
                  3.02911200e-10 * tc[2] - 2.70134040e-14 * tc[3];
  }

  // species with no change at a midpoint T
  // species 5: H
  species[5] = +0.00000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3];
}

// compute an approx to the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian_precond(
  amrex::Real* J, const amrex::Real* sc, const amrex::Real T, const int HP)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  for (int i = 0; i < 169; i++) {
    J[i] = 0.0;
  }

  amrex::Real wdot[12];
  for (auto& val : wdot) {
    val = 0.0;
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
  for (int k = 0; k < 12; ++k) {
    mixture += sc[k];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[12];
  gibbs(g_RT, tc);

  // compute the species enthalpy
  amrex::Real h_RT[12];
  speciesEnthalpy(h_RT, tc);

  amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  amrex::Real dqdci, dcdc_fac, dqdc[12];
  amrex::Real Pr, fPr, F, k_0, logPr;
  amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
  amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
  amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const amrex::Real ln10 = log(10.0);
  const amrex::Real log10e = 1.0 / log(10.0);
  // reaction 8: H + O2 (+M) <=> HO2 (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + sc[0] + (0.78 - 1) * sc[1] + (11 - 1) * sc[4] +
          (1.9 - 1) * sc[8] + (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[1] * sc[5];
  k_f = 1475000 * exp(0.6 * tc[0]);
  dlnkfdT = 0.6 * invT;
  // pressure-fall-off
  k_0 = 6.366e+20 * exp(-1.72 * tc[0] - (264.088106214317) * invT);
  Pr = 1e-12 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -1.72 * invT + (264.088106214317) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.8) * exp(-T / 1e-30);
  Fcent2 = 0.8 * exp(-T / 1e+30);
  Fcent3 = 0.;
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT = log10e / Fcent * (-Fcent1 / 1e-30 - Fcent2 / 1e+30);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[6];
  Kc = refCinv * exp(g_RT[1] + g_RT[5] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[6]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[1] -= q; // O2
  wdot[5] -= q; // H
  wdot[6] += q; // HO2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = 2 * dcdc_fac;
  dqdc[1] = 0.78 * dcdc_fac + k_f * sc[5];
  dqdc[2] = dcdc_fac;
  dqdc[3] = dcdc_fac;
  dqdc[4] = 11 * dcdc_fac;
  dqdc[5] = dcdc_fac + k_f * sc[1];
  dqdc[6] = dcdc_fac - k_r;
  dqdc[7] = dcdc_fac;
  dqdc[8] = 1.9 * dcdc_fac;
  dqdc[9] = 3.8 * dcdc_fac;
  dqdc[10] = dcdc_fac;
  dqdc[11] = dcdc_fac;
  for (int k = 0; k < 12; k++) {
    J[13 * k + 1] -= dqdc[k];
    J[13 * k + 5] -= dqdc[k];
    J[13 * k + 6] += dqdc[k];
  }
  J[157] -= dqdT; // dwdot[O2]/dT
  J[161] -= dqdT; // dwdot[H]/dT
  J[162] += dqdT; // dwdot[HO2]/dT

  // reaction 15: H2O2 (+M) <=> 2 OH (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (12 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[7];
  k_f = 295100000000000 * exp(-(24370.7831249226) * invT);
  dlnkfdT = (24370.7831249226) * invT2;
  // pressure-fall-off
  k_0 = 1.202e+17 * exp(-(22896.3582941147) * invT);
  Pr = 1e-6 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = (22896.3582941147) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.5) * exp(-T / 1e-30);
  Fcent2 = 0.5 * exp(-T / 1e+30);
  Fcent3 = 0.;
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT = log10e / Fcent * (-Fcent1 / 1e-30 - Fcent2 / 1e+30);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = (sc[3] * sc[3]);
  Kc = refC * exp(-2.000000 * g_RT[3] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[7]) + (2.000000 * h_RT[3]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[3] += 2 * q; // OH
  wdot[7] -= q;     // H2O2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = 2.5 * dcdc_fac;
  dqdc[1] = dcdc_fac;
  dqdc[2] = dcdc_fac;
  dqdc[3] = dcdc_fac - k_r * 2.000000 * sc[3];
  dqdc[4] = 12 * dcdc_fac;
  dqdc[5] = dcdc_fac;
  dqdc[6] = dcdc_fac;
  dqdc[7] = dcdc_fac + k_f;
  dqdc[8] = 1.9 * dcdc_fac;
  dqdc[9] = 3.8 * dcdc_fac;
  dqdc[10] = dcdc_fac;
  dqdc[11] = dcdc_fac;
  for (int k = 0; k < 12; k++) {
    J[13 * k + 3] += 2 * dqdc[k];
    J[13 * k + 7] -= dqdc[k];
  }
  J[159] += 2 * dqdT; // dwdot[OH]/dT
  J[163] -= dqdT;     // dwdot[H2O2]/dT

  // reaction 21: CO + O (+M) <=> CO2 (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (12 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[2] * sc[8];
  k_f = 18000 * exp(-(1199.66853127845) * invT);
  dlnkfdT = (1199.66853127845) * invT2;
  // pressure-fall-off
  k_0 = 1.55e+24 * exp(-2.79 * tc[0] - (2108.98104638758) * invT);
  Pr = 1e-12 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -2.79 * invT + (2108.98104638758) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Lindemann form
  F = 1.0;
  // dlogFdlogPr is 0.0 and unused
  dlogFdT = 0.0;
  // reverse
  phi_r = sc[9];
  Kc = refCinv * exp(g_RT[2] + g_RT[8] - g_RT[9]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[9]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[2] -= q; // O
  wdot[8] -= q; // CO
  wdot[9] += q; // CO2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = 2.5 * dcdc_fac;
  dqdc[1] = dcdc_fac;
  dqdc[2] = dcdc_fac + k_f * sc[8];
  dqdc[3] = dcdc_fac;
  dqdc[4] = 12 * dcdc_fac;
  dqdc[5] = dcdc_fac;
  dqdc[6] = dcdc_fac;
  dqdc[7] = dcdc_fac;
  dqdc[8] = 1.9 * dcdc_fac + k_f * sc[2];
  dqdc[9] = 3.8 * dcdc_fac - k_r;
  dqdc[10] = dcdc_fac;
  dqdc[11] = dcdc_fac;
  for (int k = 0; k < 12; k++) {
    J[13 * k + 2] -= dqdc[k];
    J[13 * k + 8] -= dqdc[k];
    J[13 * k + 9] += dqdc[k];
  }
  J[158] -= dqdT; // dwdot[O]/dT
  J[164] -= dqdT; // dwdot[CO]/dT
  J[165] += dqdT; // dwdot[CO2]/dT

  // reaction 4: H2 + M <=> 2 H + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (12 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[0];
  k_f = 45770000000000 * exp(-1.4 * tc[0] - (52525.7555766966) * invT);
  dlnkfdT = -1.4 * invT + (52525.7555766966) * invT2;
  // reverse
  phi_r = (sc[5] * sc[5]);
  Kc = refC * exp(g_RT[0] - 2.000000 * g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0]) + (2.000000 * h_RT[5]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] -= q;     // H2
  wdot[5] += 2 * q; // H
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = 2.5 * q_nocor + k_f;
  dqdc[1] = q_nocor;
  dqdc[2] = q_nocor;
  dqdc[3] = q_nocor;
  dqdc[4] = 12 * q_nocor;
  dqdc[5] = q_nocor - k_r * 2.000000 * sc[5];
  dqdc[6] = q_nocor;
  dqdc[7] = q_nocor;
  dqdc[8] = 1.9 * q_nocor;
  dqdc[9] = 3.8 * q_nocor;
  dqdc[10] = q_nocor;
  dqdc[11] = q_nocor;
  for (int k = 0; k < 12; k++) {
    J[13 * k + 0] -= dqdc[k];
    J[13 * k + 5] += 2 * dqdc[k];
  }
  J[156] -= dqdT;     // dwdot[H2]/dT
  J[161] += 2 * dqdT; // dwdot[H]/dT

  // reaction 5: 2 O + M <=> O2 + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (12 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = (sc[2] * sc[2]);
  k_f = 6165 * exp(-0.5 * tc[0]);
  dlnkfdT = -0.5 * invT;
  // reverse
  phi_r = sc[1];
  Kc = refCinv * exp(-g_RT[1] + 2.000000 * g_RT[2]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[2]) + (h_RT[1]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;     // O2
  wdot[2] -= 2 * q; // O
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = 2.5 * q_nocor;
  dqdc[1] = q_nocor - k_r;
  dqdc[2] = q_nocor + k_f * 2.000000 * sc[2];
  dqdc[3] = q_nocor;
  dqdc[4] = 12 * q_nocor;
  dqdc[5] = q_nocor;
  dqdc[6] = q_nocor;
  dqdc[7] = q_nocor;
  dqdc[8] = 1.9 * q_nocor;
  dqdc[9] = 3.8 * q_nocor;
  dqdc[10] = q_nocor;
  dqdc[11] = q_nocor;
  for (int k = 0; k < 12; k++) {
    J[13 * k + 1] += dqdc[k];
    J[13 * k + 2] += -2 * dqdc[k];
  }
  J[157] += dqdT;      // dwdot[O2]/dT
  J[158] += -2 * dqdT; // dwdot[O]/dT

  // reaction 6: H + O + M <=> OH + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (12 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[2] * sc[5];
  k_f = 4714000 * exp(-1 * tc[0]);
  dlnkfdT = -1 * invT;
  // reverse
  phi_r = sc[3];
  Kc = refCinv * exp(g_RT[2] - g_RT[3] + g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[3]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // O
  wdot[3] += q; // OH
  wdot[5] -= q; // H
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = 2.5 * q_nocor;
  dqdc[1] = q_nocor;
  dqdc[2] = q_nocor + k_f * sc[5];
  dqdc[3] = q_nocor - k_r;
  dqdc[4] = 12 * q_nocor;
  dqdc[5] = q_nocor + k_f * sc[2];
  dqdc[6] = q_nocor;
  dqdc[7] = q_nocor;
  dqdc[8] = 1.9 * q_nocor;
  dqdc[9] = 3.8 * q_nocor;
  dqdc[10] = q_nocor;
  dqdc[11] = q_nocor;
  for (int k = 0; k < 12; k++) {
    J[13 * k + 2] -= dqdc[k];
    J[13 * k + 3] += dqdc[k];
    J[13 * k + 5] -= dqdc[k];
  }
  J[158] -= dqdT; // dwdot[O]/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[161] -= dqdT; // dwdot[H]/dT

  // reaction 7: H + OH + M <=> H2O + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (12 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[3] * sc[5];
  k_f = 38000000000 * exp(-2 * tc[0]);
  dlnkfdT = -2 * invT;
  // reverse
  phi_r = sc[4];
  Kc = refCinv * exp(g_RT[3] - g_RT[4] + g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[4]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q; // OH
  wdot[4] += q; // H2O
  wdot[5] -= q; // H
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = 2.5 * q_nocor;
  dqdc[1] = q_nocor;
  dqdc[2] = q_nocor;
  dqdc[3] = q_nocor + k_f * sc[5];
  dqdc[4] = 12 * q_nocor - k_r;
  dqdc[5] = q_nocor + k_f * sc[3];
  dqdc[6] = q_nocor;
  dqdc[7] = q_nocor;
  dqdc[8] = 1.9 * q_nocor;
  dqdc[9] = 3.8 * q_nocor;
  dqdc[10] = q_nocor;
  dqdc[11] = q_nocor;
  for (int k = 0; k < 12; k++) {
    J[13 * k + 3] -= dqdc[k];
    J[13 * k + 4] += dqdc[k];
    J[13 * k + 5] -= dqdc[k];
  }
  J[159] -= dqdT; // dwdot[OH]/dT
  J[160] += dqdT; // dwdot[H2O]/dT
  J[161] -= dqdT; // dwdot[H]/dT

  // reaction 25: HCO + M <=> CO + H + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (6 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[10];
  k_f = 474850 * exp(0.659 * tc[0] - (7484.8446871794) * invT);
  dlnkfdT = 0.659 * invT + (7484.8446871794) * invT2;
  // reverse
  phi_r = sc[5] * sc[8];
  Kc = refC * exp(-g_RT[5] - g_RT[8] + g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[10]) + (h_RT[5] + h_RT[8]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H
  wdot[8] += q;  // CO
  wdot[10] -= q; // HCO
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = 2.5 * q_nocor;
  dqdc[1] = q_nocor;
  dqdc[2] = q_nocor;
  dqdc[3] = q_nocor;
  dqdc[4] = 6 * q_nocor;
  dqdc[5] = q_nocor - k_r * sc[8];
  dqdc[6] = q_nocor;
  dqdc[7] = q_nocor;
  dqdc[8] = 1.9 * q_nocor - k_r * sc[5];
  dqdc[9] = 3.8 * q_nocor;
  dqdc[10] = q_nocor + k_f;
  dqdc[11] = q_nocor;
  for (int k = 0; k < 12; k++) {
    J[13 * k + 5] += dqdc[k];
    J[13 * k + 8] += dqdc[k];
    J[13 * k + 10] -= dqdc[k];
  }
  J[161] += dqdT; // dwdot[H]/dT
  J[164] += dqdT; // dwdot[CO]/dT
  J[166] -= dqdT; // dwdot[HCO]/dT

  // reaction 0: H + O2 <=> O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[5];
  k_f = 3547000000 * exp(-0.406 * tc[0] - (8352.89343569254) * invT);
  dlnkfdT = -0.406 * invT + (8352.89343569254) * invT2;
  // reverse
  phi_r = sc[2] * sc[3];
  Kc = exp(g_RT[1] - g_RT[2] - g_RT[3] + g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[2] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q; // O2
  wdot[2] += q; // O
  wdot[3] += q; // OH
  wdot[5] -= q; // H
  // d()/d[O2]
  dqdci = +k_f * sc[5];
  J[14] -= dqdci; // dwdot[O2]/d[O2]
  J[15] += dqdci; // dwdot[O]/d[O2]
  J[16] += dqdci; // dwdot[OH]/d[O2]
  J[18] -= dqdci; // dwdot[H]/d[O2]
  // d()/d[O]
  dqdci = -k_r * sc[3];
  J[27] -= dqdci; // dwdot[O2]/d[O]
  J[28] += dqdci; // dwdot[O]/d[O]
  J[29] += dqdci; // dwdot[OH]/d[O]
  J[31] -= dqdci; // dwdot[H]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[2];
  J[40] -= dqdci; // dwdot[O2]/d[OH]
  J[41] += dqdci; // dwdot[O]/d[OH]
  J[42] += dqdci; // dwdot[OH]/d[OH]
  J[44] -= dqdci; // dwdot[H]/d[OH]
  // d()/d[H]
  dqdci = +k_f * sc[1];
  J[66] -= dqdci; // dwdot[O2]/d[H]
  J[67] += dqdci; // dwdot[O]/d[H]
  J[68] += dqdci; // dwdot[OH]/d[H]
  J[70] -= dqdci; // dwdot[H]/d[H]
  // d()/dT
  J[157] -= dqdT; // dwdot[O2]/dT
  J[158] += dqdT; // dwdot[O]/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[161] -= dqdT; // dwdot[H]/dT

  // reaction 1: H2 + O <=> H + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[2];
  k_f = 0.0508 * exp(2.67 * tc[0] - (3165.23282791169) * invT);
  dlnkfdT = 2.67 * invT + (3165.23282791169) * invT2;
  // reverse
  phi_r = sc[3] * sc[5];
  Kc = exp(g_RT[0] + g_RT[2] - g_RT[3] - g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] -= q; // H2
  wdot[2] -= q; // O
  wdot[3] += q; // OH
  wdot[5] += q; // H
  // d()/d[H2]
  dqdci = +k_f * sc[2];
  J[0] -= dqdci; // dwdot[H2]/d[H2]
  J[2] -= dqdci; // dwdot[O]/d[H2]
  J[3] += dqdci; // dwdot[OH]/d[H2]
  J[5] += dqdci; // dwdot[H]/d[H2]
  // d()/d[O]
  dqdci = +k_f * sc[0];
  J[26] -= dqdci; // dwdot[H2]/d[O]
  J[28] -= dqdci; // dwdot[O]/d[O]
  J[29] += dqdci; // dwdot[OH]/d[O]
  J[31] += dqdci; // dwdot[H]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[5];
  J[39] -= dqdci; // dwdot[H2]/d[OH]
  J[41] -= dqdci; // dwdot[O]/d[OH]
  J[42] += dqdci; // dwdot[OH]/d[OH]
  J[44] += dqdci; // dwdot[H]/d[OH]
  // d()/d[H]
  dqdci = -k_r * sc[3];
  J[65] -= dqdci; // dwdot[H2]/d[H]
  J[67] -= dqdci; // dwdot[O]/d[H]
  J[68] += dqdci; // dwdot[OH]/d[H]
  J[70] += dqdci; // dwdot[H]/d[H]
  // d()/dT
  J[156] -= dqdT; // dwdot[H2]/dT
  J[158] -= dqdT; // dwdot[O]/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[161] += dqdT; // dwdot[H]/dT

  // reaction 2: H2 + OH <=> H + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[3];
  k_f = 216 * exp(1.51 * tc[0] - (1726.03316371019) * invT);
  dlnkfdT = 1.51 * invT + (1726.03316371019) * invT2;
  // reverse
  phi_r = sc[4] * sc[5];
  Kc = exp(g_RT[0] + g_RT[3] - g_RT[4] - g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[3]) + (h_RT[4] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] -= q; // H2
  wdot[3] -= q; // OH
  wdot[4] += q; // H2O
  wdot[5] += q; // H
  // d()/d[H2]
  dqdci = +k_f * sc[3];
  J[0] -= dqdci; // dwdot[H2]/d[H2]
  J[3] -= dqdci; // dwdot[OH]/d[H2]
  J[4] += dqdci; // dwdot[H2O]/d[H2]
  J[5] += dqdci; // dwdot[H]/d[H2]
  // d()/d[OH]
  dqdci = +k_f * sc[0];
  J[39] -= dqdci; // dwdot[H2]/d[OH]
  J[42] -= dqdci; // dwdot[OH]/d[OH]
  J[43] += dqdci; // dwdot[H2O]/d[OH]
  J[44] += dqdci; // dwdot[H]/d[OH]
  // d()/d[H2O]
  dqdci = -k_r * sc[5];
  J[52] -= dqdci; // dwdot[H2]/d[H2O]
  J[55] -= dqdci; // dwdot[OH]/d[H2O]
  J[56] += dqdci; // dwdot[H2O]/d[H2O]
  J[57] += dqdci; // dwdot[H]/d[H2O]
  // d()/d[H]
  dqdci = -k_r * sc[4];
  J[65] -= dqdci; // dwdot[H2]/d[H]
  J[68] -= dqdci; // dwdot[OH]/d[H]
  J[69] += dqdci; // dwdot[H2O]/d[H]
  J[70] += dqdci; // dwdot[H]/d[H]
  // d()/dT
  J[156] -= dqdT; // dwdot[H2]/dT
  J[159] -= dqdT; // dwdot[OH]/dT
  J[160] += dqdT; // dwdot[H2O]/dT
  J[161] += dqdT; // dwdot[H]/dT

  // reaction 3: H2O + O <=> 2 OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[4];
  k_f = 2.97 * exp(2.02 * tc[0] - (6743.10332178324) * invT);
  dlnkfdT = 2.02 * invT + (6743.10332178324) * invT2;
  // reverse
  phi_r = (sc[3] * sc[3]);
  Kc = exp(g_RT[2] - 2.000000 * g_RT[3] + g_RT[4]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2.000000 * h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q;     // O
  wdot[3] += 2 * q; // OH
  wdot[4] -= q;     // H2O
  // d()/d[O]
  dqdci = +k_f * sc[4];
  J[28] -= dqdci;     // dwdot[O]/d[O]
  J[29] += 2 * dqdci; // dwdot[OH]/d[O]
  J[30] -= dqdci;     // dwdot[H2O]/d[O]
  // d()/d[OH]
  dqdci = -k_r * 2.000000 * sc[3];
  J[41] -= dqdci;     // dwdot[O]/d[OH]
  J[42] += 2 * dqdci; // dwdot[OH]/d[OH]
  J[43] -= dqdci;     // dwdot[H2O]/d[OH]
  // d()/d[H2O]
  dqdci = +k_f * sc[2];
  J[54] -= dqdci;     // dwdot[O]/d[H2O]
  J[55] += 2 * dqdci; // dwdot[OH]/d[H2O]
  J[56] -= dqdci;     // dwdot[H2O]/d[H2O]
  // d()/dT
  J[158] -= dqdT;     // dwdot[O]/dT
  J[159] += 2 * dqdT; // dwdot[OH]/dT
  J[160] -= dqdT;     // dwdot[H2O]/dT

  // reaction 9: H + HO2 <=> H2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[6];
  k_f = 16600000 * exp(-(414.147315957284) * invT);
  dlnkfdT = (414.147315957284) * invT2;
  // reverse
  phi_r = sc[0] * sc[1];
  Kc = exp(-g_RT[0] - g_RT[1] + g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (h_RT[0] + h_RT[1]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] += q; // H2
  wdot[1] += q; // O2
  wdot[5] -= q; // H
  wdot[6] -= q; // HO2
  // d()/d[H2]
  dqdci = -k_r * sc[1];
  J[0] += dqdci; // dwdot[H2]/d[H2]
  J[1] += dqdci; // dwdot[O2]/d[H2]
  J[5] -= dqdci; // dwdot[H]/d[H2]
  J[6] -= dqdci; // dwdot[HO2]/d[H2]
  // d()/d[O2]
  dqdci = -k_r * sc[0];
  J[13] += dqdci; // dwdot[H2]/d[O2]
  J[14] += dqdci; // dwdot[O2]/d[O2]
  J[18] -= dqdci; // dwdot[H]/d[O2]
  J[19] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[H]
  dqdci = +k_f * sc[6];
  J[65] += dqdci; // dwdot[H2]/d[H]
  J[66] += dqdci; // dwdot[O2]/d[H]
  J[70] -= dqdci; // dwdot[H]/d[H]
  J[71] -= dqdci; // dwdot[HO2]/d[H]
  // d()/d[HO2]
  dqdci = +k_f * sc[5];
  J[78] += dqdci; // dwdot[H2]/d[HO2]
  J[79] += dqdci; // dwdot[O2]/d[HO2]
  J[83] -= dqdci; // dwdot[H]/d[HO2]
  J[84] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[156] += dqdT; // dwdot[H2]/dT
  J[157] += dqdT; // dwdot[O2]/dT
  J[161] -= dqdT; // dwdot[H]/dT
  J[162] -= dqdT; // dwdot[HO2]/dT

  // reaction 10: H + HO2 <=> 2 OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[6];
  k_f = 70790000 * exp(-(148.448916412392) * invT);
  dlnkfdT = (148.448916412392) * invT2;
  // reverse
  phi_r = (sc[3] * sc[3]);
  Kc = exp(-2.000000 * g_RT[3] + g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (2.000000 * h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += 2 * q; // OH
  wdot[5] -= q;     // H
  wdot[6] -= q;     // HO2
  // d()/d[OH]
  dqdci = -k_r * 2.000000 * sc[3];
  J[42] += 2 * dqdci; // dwdot[OH]/d[OH]
  J[44] -= dqdci;     // dwdot[H]/d[OH]
  J[45] -= dqdci;     // dwdot[HO2]/d[OH]
  // d()/d[H]
  dqdci = +k_f * sc[6];
  J[68] += 2 * dqdci; // dwdot[OH]/d[H]
  J[70] -= dqdci;     // dwdot[H]/d[H]
  J[71] -= dqdci;     // dwdot[HO2]/d[H]
  // d()/d[HO2]
  dqdci = +k_f * sc[5];
  J[81] += 2 * dqdci; // dwdot[OH]/d[HO2]
  J[83] -= dqdci;     // dwdot[H]/d[HO2]
  J[84] -= dqdci;     // dwdot[HO2]/d[HO2]
  // d()/dT
  J[159] += 2 * dqdT; // dwdot[OH]/dT
  J[161] -= dqdT;     // dwdot[H]/dT
  J[162] -= dqdT;     // dwdot[HO2]/dT

  // reaction 11: HO2 + O <=> O2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[6];
  k_f = 32500000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[3];
  Kc = exp(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[1] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q; // O2
  wdot[2] -= q; // O
  wdot[3] += q; // OH
  wdot[6] -= q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[3];
  J[14] += dqdci; // dwdot[O2]/d[O2]
  J[15] -= dqdci; // dwdot[O]/d[O2]
  J[16] += dqdci; // dwdot[OH]/d[O2]
  J[19] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[O]
  dqdci = +k_f * sc[6];
  J[27] += dqdci; // dwdot[O2]/d[O]
  J[28] -= dqdci; // dwdot[O]/d[O]
  J[29] += dqdci; // dwdot[OH]/d[O]
  J[32] -= dqdci; // dwdot[HO2]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[1];
  J[40] += dqdci; // dwdot[O2]/d[OH]
  J[41] -= dqdci; // dwdot[O]/d[OH]
  J[42] += dqdci; // dwdot[OH]/d[OH]
  J[45] -= dqdci; // dwdot[HO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[2];
  J[79] += dqdci; // dwdot[O2]/d[HO2]
  J[80] -= dqdci; // dwdot[O]/d[HO2]
  J[81] += dqdci; // dwdot[OH]/d[HO2]
  J[84] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[157] += dqdT; // dwdot[O2]/dT
  J[158] -= dqdT; // dwdot[O]/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[162] -= dqdT; // dwdot[HO2]/dT

  // reaction 12: HO2 + OH <=> H2O + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[6];
  k_f = 28900000 * exp(-(-250.098682904946) * invT);
  dlnkfdT = (-250.098682904946) * invT2;
  // reverse
  phi_r = sc[1] * sc[4];
  Kc = exp(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[1] + h_RT[4]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q; // O2
  wdot[3] -= q; // OH
  wdot[4] += q; // H2O
  wdot[6] -= q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[4];
  J[14] += dqdci; // dwdot[O2]/d[O2]
  J[16] -= dqdci; // dwdot[OH]/d[O2]
  J[17] += dqdci; // dwdot[H2O]/d[O2]
  J[19] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[OH]
  dqdci = +k_f * sc[6];
  J[40] += dqdci; // dwdot[O2]/d[OH]
  J[42] -= dqdci; // dwdot[OH]/d[OH]
  J[43] += dqdci; // dwdot[H2O]/d[OH]
  J[45] -= dqdci; // dwdot[HO2]/d[OH]
  // d()/d[H2O]
  dqdci = -k_r * sc[1];
  J[53] += dqdci; // dwdot[O2]/d[H2O]
  J[55] -= dqdci; // dwdot[OH]/d[H2O]
  J[56] += dqdci; // dwdot[H2O]/d[H2O]
  J[58] -= dqdci; // dwdot[HO2]/d[H2O]
  // d()/d[HO2]
  dqdci = +k_f * sc[3];
  J[79] += dqdci; // dwdot[O2]/d[HO2]
  J[81] -= dqdci; // dwdot[OH]/d[HO2]
  J[82] += dqdci; // dwdot[H2O]/d[HO2]
  J[84] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[157] += dqdT; // dwdot[O2]/dT
  J[159] -= dqdT; // dwdot[OH]/dT
  J[160] += dqdT; // dwdot[H2O]/dT
  J[162] -= dqdT; // dwdot[HO2]/dT

  // reaction 13: 2 HO2 <=> H2O2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[6] * sc[6]);
  k_f = 420000000 * exp(-(6029.54208967215) * invT);
  dlnkfdT = (6029.54208967215) * invT2;
  // reverse
  phi_r = sc[1] * sc[7];
  Kc = exp(-g_RT[1] + 2.000000 * g_RT[6] - g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[6]) + (h_RT[1] + h_RT[7]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;     // O2
  wdot[6] -= 2 * q; // HO2
  wdot[7] += q;     // H2O2
  // d()/d[O2]
  dqdci = -k_r * sc[7];
  J[14] += dqdci;      // dwdot[O2]/d[O2]
  J[19] += -2 * dqdci; // dwdot[HO2]/d[O2]
  J[20] += dqdci;      // dwdot[H2O2]/d[O2]
  // d()/d[HO2]
  dqdci = +k_f * 2.000000 * sc[6];
  J[79] += dqdci;      // dwdot[O2]/d[HO2]
  J[84] += -2 * dqdci; // dwdot[HO2]/d[HO2]
  J[85] += dqdci;      // dwdot[H2O2]/d[HO2]
  // d()/d[H2O2]
  dqdci = -k_r * sc[1];
  J[92] += dqdci;      // dwdot[O2]/d[H2O2]
  J[97] += -2 * dqdci; // dwdot[HO2]/d[H2O2]
  J[98] += dqdci;      // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[157] += dqdT;      // dwdot[O2]/dT
  J[162] += -2 * dqdT; // dwdot[HO2]/dT
  J[163] += dqdT;      // dwdot[H2O2]/dT

  // reaction 14: 2 HO2 <=> H2O2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[6] * sc[6]);
  k_f = 130000 * exp(-(-819.89091359563) * invT);
  dlnkfdT = (-819.89091359563) * invT2;
  // reverse
  phi_r = sc[1] * sc[7];
  Kc = exp(-g_RT[1] + 2.000000 * g_RT[6] - g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[6]) + (h_RT[1] + h_RT[7]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;     // O2
  wdot[6] -= 2 * q; // HO2
  wdot[7] += q;     // H2O2
  // d()/d[O2]
  dqdci = -k_r * sc[7];
  J[14] += dqdci;      // dwdot[O2]/d[O2]
  J[19] += -2 * dqdci; // dwdot[HO2]/d[O2]
  J[20] += dqdci;      // dwdot[H2O2]/d[O2]
  // d()/d[HO2]
  dqdci = +k_f * 2.000000 * sc[6];
  J[79] += dqdci;      // dwdot[O2]/d[HO2]
  J[84] += -2 * dqdci; // dwdot[HO2]/d[HO2]
  J[85] += dqdci;      // dwdot[H2O2]/d[HO2]
  // d()/d[H2O2]
  dqdci = -k_r * sc[1];
  J[92] += dqdci;      // dwdot[O2]/d[H2O2]
  J[97] += -2 * dqdci; // dwdot[HO2]/d[H2O2]
  J[98] += dqdci;      // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[157] += dqdT;      // dwdot[O2]/dT
  J[162] += -2 * dqdT; // dwdot[HO2]/dT
  J[163] += dqdT;      // dwdot[H2O2]/dT

  // reaction 16: H + H2O2 <=> H2O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[7];
  k_f = 24100000 * exp(-(1997.77016324474) * invT);
  dlnkfdT = (1997.77016324474) * invT2;
  // reverse
  phi_r = sc[3] * sc[4];
  Kc = exp(-g_RT[3] - g_RT[4] + g_RT[5] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[3] + h_RT[4]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q; // OH
  wdot[4] += q; // H2O
  wdot[5] -= q; // H
  wdot[7] -= q; // H2O2
  // d()/d[OH]
  dqdci = -k_r * sc[4];
  J[42] += dqdci; // dwdot[OH]/d[OH]
  J[43] += dqdci; // dwdot[H2O]/d[OH]
  J[44] -= dqdci; // dwdot[H]/d[OH]
  J[46] -= dqdci; // dwdot[H2O2]/d[OH]
  // d()/d[H2O]
  dqdci = -k_r * sc[3];
  J[55] += dqdci; // dwdot[OH]/d[H2O]
  J[56] += dqdci; // dwdot[H2O]/d[H2O]
  J[57] -= dqdci; // dwdot[H]/d[H2O]
  J[59] -= dqdci; // dwdot[H2O2]/d[H2O]
  // d()/d[H]
  dqdci = +k_f * sc[7];
  J[68] += dqdci; // dwdot[OH]/d[H]
  J[69] += dqdci; // dwdot[H2O]/d[H]
  J[70] -= dqdci; // dwdot[H]/d[H]
  J[72] -= dqdci; // dwdot[H2O2]/d[H]
  // d()/d[H2O2]
  dqdci = +k_f * sc[5];
  J[94] += dqdci; // dwdot[OH]/d[H2O2]
  J[95] += dqdci; // dwdot[H2O]/d[H2O2]
  J[96] -= dqdci; // dwdot[H]/d[H2O2]
  J[98] -= dqdci; // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[160] += dqdT; // dwdot[H2O]/dT
  J[161] -= dqdT; // dwdot[H]/dT
  J[163] -= dqdT; // dwdot[H2O2]/dT

  // reaction 17: H + H2O2 <=> H2 + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[7];
  k_f = 48200000 * exp(-(4000.57249314752) * invT);
  dlnkfdT = (4000.57249314752) * invT2;
  // reverse
  phi_r = sc[0] * sc[6];
  Kc = exp(-g_RT[0] + g_RT[5] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[0] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] += q; // H2
  wdot[5] -= q; // H
  wdot[6] += q; // HO2
  wdot[7] -= q; // H2O2
  // d()/d[H2]
  dqdci = -k_r * sc[6];
  J[0] += dqdci; // dwdot[H2]/d[H2]
  J[5] -= dqdci; // dwdot[H]/d[H2]
  J[6] += dqdci; // dwdot[HO2]/d[H2]
  J[7] -= dqdci; // dwdot[H2O2]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[7];
  J[65] += dqdci; // dwdot[H2]/d[H]
  J[70] -= dqdci; // dwdot[H]/d[H]
  J[71] += dqdci; // dwdot[HO2]/d[H]
  J[72] -= dqdci; // dwdot[H2O2]/d[H]
  // d()/d[HO2]
  dqdci = -k_r * sc[0];
  J[78] += dqdci; // dwdot[H2]/d[HO2]
  J[83] -= dqdci; // dwdot[H]/d[HO2]
  J[84] += dqdci; // dwdot[HO2]/d[HO2]
  J[85] -= dqdci; // dwdot[H2O2]/d[HO2]
  // d()/d[H2O2]
  dqdci = +k_f * sc[5];
  J[91] += dqdci; // dwdot[H2]/d[H2O2]
  J[96] -= dqdci; // dwdot[H]/d[H2O2]
  J[97] += dqdci; // dwdot[HO2]/d[H2O2]
  J[98] -= dqdci; // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[156] += dqdT; // dwdot[H2]/dT
  J[161] -= dqdT; // dwdot[H]/dT
  J[162] += dqdT; // dwdot[HO2]/dT
  J[163] -= dqdT; // dwdot[H2O2]/dT

  // reaction 18: H2O2 + O <=> HO2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[7];
  k_f = 9.55 * exp(2 * tc[0] - (1997.77016324474) * invT);
  dlnkfdT = 2 * invT + (1997.77016324474) * invT2;
  // reverse
  phi_r = sc[3] * sc[6];
  Kc = exp(g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[3] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // O
  wdot[3] += q; // OH
  wdot[6] += q; // HO2
  wdot[7] -= q; // H2O2
  // d()/d[O]
  dqdci = +k_f * sc[7];
  J[28] -= dqdci; // dwdot[O]/d[O]
  J[29] += dqdci; // dwdot[OH]/d[O]
  J[32] += dqdci; // dwdot[HO2]/d[O]
  J[33] -= dqdci; // dwdot[H2O2]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[6];
  J[41] -= dqdci; // dwdot[O]/d[OH]
  J[42] += dqdci; // dwdot[OH]/d[OH]
  J[45] += dqdci; // dwdot[HO2]/d[OH]
  J[46] -= dqdci; // dwdot[H2O2]/d[OH]
  // d()/d[HO2]
  dqdci = -k_r * sc[3];
  J[80] -= dqdci; // dwdot[O]/d[HO2]
  J[81] += dqdci; // dwdot[OH]/d[HO2]
  J[84] += dqdci; // dwdot[HO2]/d[HO2]
  J[85] -= dqdci; // dwdot[H2O2]/d[HO2]
  // d()/d[H2O2]
  dqdci = +k_f * sc[2];
  J[93] -= dqdci; // dwdot[O]/d[H2O2]
  J[94] += dqdci; // dwdot[OH]/d[H2O2]
  J[97] += dqdci; // dwdot[HO2]/d[H2O2]
  J[98] -= dqdci; // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[158] -= dqdT; // dwdot[O]/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[162] += dqdT; // dwdot[HO2]/dT
  J[163] -= dqdT; // dwdot[H2O2]/dT

  // reaction 19: H2O2 + OH <=> H2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[7];
  k_f = 1000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[4] * sc[6];
  Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q; // OH
  wdot[4] += q; // H2O
  wdot[6] += q; // HO2
  wdot[7] -= q; // H2O2
  // d()/d[OH]
  dqdci = +k_f * sc[7];
  J[42] -= dqdci; // dwdot[OH]/d[OH]
  J[43] += dqdci; // dwdot[H2O]/d[OH]
  J[45] += dqdci; // dwdot[HO2]/d[OH]
  J[46] -= dqdci; // dwdot[H2O2]/d[OH]
  // d()/d[H2O]
  dqdci = -k_r * sc[6];
  J[55] -= dqdci; // dwdot[OH]/d[H2O]
  J[56] += dqdci; // dwdot[H2O]/d[H2O]
  J[58] += dqdci; // dwdot[HO2]/d[H2O]
  J[59] -= dqdci; // dwdot[H2O2]/d[H2O]
  // d()/d[HO2]
  dqdci = -k_r * sc[4];
  J[81] -= dqdci; // dwdot[OH]/d[HO2]
  J[82] += dqdci; // dwdot[H2O]/d[HO2]
  J[84] += dqdci; // dwdot[HO2]/d[HO2]
  J[85] -= dqdci; // dwdot[H2O2]/d[HO2]
  // d()/d[H2O2]
  dqdci = +k_f * sc[3];
  J[94] -= dqdci; // dwdot[OH]/d[H2O2]
  J[95] += dqdci; // dwdot[H2O]/d[H2O2]
  J[97] += dqdci; // dwdot[HO2]/d[H2O2]
  J[98] -= dqdci; // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[159] -= dqdT; // dwdot[OH]/dT
  J[160] += dqdT; // dwdot[H2O]/dT
  J[162] += dqdT; // dwdot[HO2]/dT
  J[163] -= dqdT; // dwdot[H2O2]/dT

  // reaction 20: H2O2 + OH <=> H2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[7];
  k_f = 580000000 * exp(-(4809.24167509571) * invT);
  dlnkfdT = (4809.24167509571) * invT2;
  // reverse
  phi_r = sc[4] * sc[6];
  Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q; // OH
  wdot[4] += q; // H2O
  wdot[6] += q; // HO2
  wdot[7] -= q; // H2O2
  // d()/d[OH]
  dqdci = +k_f * sc[7];
  J[42] -= dqdci; // dwdot[OH]/d[OH]
  J[43] += dqdci; // dwdot[H2O]/d[OH]
  J[45] += dqdci; // dwdot[HO2]/d[OH]
  J[46] -= dqdci; // dwdot[H2O2]/d[OH]
  // d()/d[H2O]
  dqdci = -k_r * sc[6];
  J[55] -= dqdci; // dwdot[OH]/d[H2O]
  J[56] += dqdci; // dwdot[H2O]/d[H2O]
  J[58] += dqdci; // dwdot[HO2]/d[H2O]
  J[59] -= dqdci; // dwdot[H2O2]/d[H2O]
  // d()/d[HO2]
  dqdci = -k_r * sc[4];
  J[81] -= dqdci; // dwdot[OH]/d[HO2]
  J[82] += dqdci; // dwdot[H2O]/d[HO2]
  J[84] += dqdci; // dwdot[HO2]/d[HO2]
  J[85] -= dqdci; // dwdot[H2O2]/d[HO2]
  // d()/d[H2O2]
  dqdci = +k_f * sc[3];
  J[94] -= dqdci; // dwdot[OH]/d[H2O2]
  J[95] += dqdci; // dwdot[H2O]/d[H2O2]
  J[97] += dqdci; // dwdot[HO2]/d[H2O2]
  J[98] -= dqdci; // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[159] -= dqdT; // dwdot[OH]/dT
  J[160] += dqdT; // dwdot[H2O]/dT
  J[162] += dqdT; // dwdot[HO2]/dT
  J[163] -= dqdT; // dwdot[H2O2]/dT

  // reaction 22: CO + O2 <=> CO2 + O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[8];
  k_f = 2530000 * exp(-(24003.4349588851) * invT);
  dlnkfdT = (24003.4349588851) * invT2;
  // reverse
  phi_r = sc[2] * sc[9];
  Kc = exp(g_RT[1] - g_RT[2] + g_RT[8] - g_RT[9]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[8]) + (h_RT[2] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q; // O2
  wdot[2] += q; // O
  wdot[8] -= q; // CO
  wdot[9] += q; // CO2
  // d()/d[O2]
  dqdci = +k_f * sc[8];
  J[14] -= dqdci; // dwdot[O2]/d[O2]
  J[15] += dqdci; // dwdot[O]/d[O2]
  J[21] -= dqdci; // dwdot[CO]/d[O2]
  J[22] += dqdci; // dwdot[CO2]/d[O2]
  // d()/d[O]
  dqdci = -k_r * sc[9];
  J[27] -= dqdci; // dwdot[O2]/d[O]
  J[28] += dqdci; // dwdot[O]/d[O]
  J[34] -= dqdci; // dwdot[CO]/d[O]
  J[35] += dqdci; // dwdot[CO2]/d[O]
  // d()/d[CO]
  dqdci = +k_f * sc[1];
  J[105] -= dqdci; // dwdot[O2]/d[CO]
  J[106] += dqdci; // dwdot[O]/d[CO]
  J[112] -= dqdci; // dwdot[CO]/d[CO]
  J[113] += dqdci; // dwdot[CO2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[2];
  J[118] -= dqdci; // dwdot[O2]/d[CO2]
  J[119] += dqdci; // dwdot[O]/d[CO2]
  J[125] -= dqdci; // dwdot[CO]/d[CO2]
  J[126] += dqdci; // dwdot[CO2]/d[CO2]
  // d()/dT
  J[157] -= dqdT; // dwdot[O2]/dT
  J[158] += dqdT; // dwdot[O]/dT
  J[164] -= dqdT; // dwdot[CO]/dT
  J[165] += dqdT; // dwdot[CO2]/dT

  // reaction 23: CO + HO2 <=> CO2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[8];
  k_f = 30100000 * exp(-(11573.9833135086) * invT);
  dlnkfdT = (11573.9833135086) * invT2;
  // reverse
  phi_r = sc[3] * sc[9];
  Kc = exp(-g_RT[3] + g_RT[6] + g_RT[8] - g_RT[9]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[3] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q; // OH
  wdot[6] -= q; // HO2
  wdot[8] -= q; // CO
  wdot[9] += q; // CO2
  // d()/d[OH]
  dqdci = -k_r * sc[9];
  J[42] += dqdci; // dwdot[OH]/d[OH]
  J[45] -= dqdci; // dwdot[HO2]/d[OH]
  J[47] -= dqdci; // dwdot[CO]/d[OH]
  J[48] += dqdci; // dwdot[CO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[8];
  J[81] += dqdci; // dwdot[OH]/d[HO2]
  J[84] -= dqdci; // dwdot[HO2]/d[HO2]
  J[86] -= dqdci; // dwdot[CO]/d[HO2]
  J[87] += dqdci; // dwdot[CO2]/d[HO2]
  // d()/d[CO]
  dqdci = +k_f * sc[6];
  J[107] += dqdci; // dwdot[OH]/d[CO]
  J[110] -= dqdci; // dwdot[HO2]/d[CO]
  J[112] -= dqdci; // dwdot[CO]/d[CO]
  J[113] += dqdci; // dwdot[CO2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[3];
  J[120] += dqdci; // dwdot[OH]/d[CO2]
  J[123] -= dqdci; // dwdot[HO2]/d[CO2]
  J[125] -= dqdci; // dwdot[CO]/d[CO2]
  J[126] += dqdci; // dwdot[CO2]/d[CO2]
  // d()/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[162] -= dqdT; // dwdot[HO2]/dT
  J[164] -= dqdT; // dwdot[CO]/dT
  J[165] += dqdT; // dwdot[CO2]/dT

  // reaction 24: CO + OH <=> CO2 + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[8];
  k_f = 0.2229 * exp(1.89 * tc[0] - (-583.077150667929) * invT);
  dlnkfdT = 1.89 * invT + (-583.077150667929) * invT2;
  // reverse
  phi_r = sc[5] * sc[9];
  Kc = exp(g_RT[3] - g_RT[5] + g_RT[8] - g_RT[9]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[5] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q; // OH
  wdot[5] += q; // H
  wdot[8] -= q; // CO
  wdot[9] += q; // CO2
  // d()/d[OH]
  dqdci = +k_f * sc[8];
  J[42] -= dqdci; // dwdot[OH]/d[OH]
  J[44] += dqdci; // dwdot[H]/d[OH]
  J[47] -= dqdci; // dwdot[CO]/d[OH]
  J[48] += dqdci; // dwdot[CO2]/d[OH]
  // d()/d[H]
  dqdci = -k_r * sc[9];
  J[68] -= dqdci; // dwdot[OH]/d[H]
  J[70] += dqdci; // dwdot[H]/d[H]
  J[73] -= dqdci; // dwdot[CO]/d[H]
  J[74] += dqdci; // dwdot[CO2]/d[H]
  // d()/d[CO]
  dqdci = +k_f * sc[3];
  J[107] -= dqdci; // dwdot[OH]/d[CO]
  J[109] += dqdci; // dwdot[H]/d[CO]
  J[112] -= dqdci; // dwdot[CO]/d[CO]
  J[113] += dqdci; // dwdot[CO2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[5];
  J[120] -= dqdci; // dwdot[OH]/d[CO2]
  J[122] += dqdci; // dwdot[H]/d[CO2]
  J[125] -= dqdci; // dwdot[CO]/d[CO2]
  J[126] += dqdci; // dwdot[CO2]/d[CO2]
  // d()/dT
  J[159] -= dqdT; // dwdot[OH]/dT
  J[161] += dqdT; // dwdot[H]/dT
  J[164] -= dqdT; // dwdot[CO]/dT
  J[165] += dqdT; // dwdot[CO2]/dT

  // reaction 26: HCO + O2 <=> CO + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[10];
  k_f = 7580000 * exp(-(206.318832979935) * invT);
  dlnkfdT = (206.318832979935) * invT2;
  // reverse
  phi_r = sc[6] * sc[8];
  Kc = exp(g_RT[1] - g_RT[6] - g_RT[8] + g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[6] + h_RT[8]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q;  // O2
  wdot[6] += q;  // HO2
  wdot[8] += q;  // CO
  wdot[10] -= q; // HCO
  // d()/d[O2]
  dqdci = +k_f * sc[10];
  J[14] -= dqdci; // dwdot[O2]/d[O2]
  J[19] += dqdci; // dwdot[HO2]/d[O2]
  J[21] += dqdci; // dwdot[CO]/d[O2]
  J[23] -= dqdci; // dwdot[HCO]/d[O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[8];
  J[79] -= dqdci; // dwdot[O2]/d[HO2]
  J[84] += dqdci; // dwdot[HO2]/d[HO2]
  J[86] += dqdci; // dwdot[CO]/d[HO2]
  J[88] -= dqdci; // dwdot[HCO]/d[HO2]
  // d()/d[CO]
  dqdci = -k_r * sc[6];
  J[105] -= dqdci; // dwdot[O2]/d[CO]
  J[110] += dqdci; // dwdot[HO2]/d[CO]
  J[112] += dqdci; // dwdot[CO]/d[CO]
  J[114] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[1];
  J[131] -= dqdci; // dwdot[O2]/d[HCO]
  J[136] += dqdci; // dwdot[HO2]/d[HCO]
  J[138] += dqdci; // dwdot[CO]/d[HCO]
  J[140] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[157] -= dqdT; // dwdot[O2]/dT
  J[162] += dqdT; // dwdot[HO2]/dT
  J[164] += dqdT; // dwdot[CO]/dT
  J[166] -= dqdT; // dwdot[HCO]/dT

  // reaction 27: H + HCO <=> CO + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[10];
  k_f = 72300000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[0] * sc[8];
  Kc = exp(-g_RT[0] + g_RT[5] - g_RT[8] + g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[10]) + (h_RT[0] + h_RT[8]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] += q;  // H2
  wdot[5] -= q;  // H
  wdot[8] += q;  // CO
  wdot[10] -= q; // HCO
  // d()/d[H2]
  dqdci = -k_r * sc[8];
  J[0] += dqdci;  // dwdot[H2]/d[H2]
  J[5] -= dqdci;  // dwdot[H]/d[H2]
  J[8] += dqdci;  // dwdot[CO]/d[H2]
  J[10] -= dqdci; // dwdot[HCO]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[10];
  J[65] += dqdci; // dwdot[H2]/d[H]
  J[70] -= dqdci; // dwdot[H]/d[H]
  J[73] += dqdci; // dwdot[CO]/d[H]
  J[75] -= dqdci; // dwdot[HCO]/d[H]
  // d()/d[CO]
  dqdci = -k_r * sc[0];
  J[104] += dqdci; // dwdot[H2]/d[CO]
  J[109] -= dqdci; // dwdot[H]/d[CO]
  J[112] += dqdci; // dwdot[CO]/d[CO]
  J[114] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[5];
  J[130] += dqdci; // dwdot[H2]/d[HCO]
  J[135] -= dqdci; // dwdot[H]/d[HCO]
  J[138] += dqdci; // dwdot[CO]/d[HCO]
  J[140] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[156] += dqdT; // dwdot[H2]/dT
  J[161] -= dqdT; // dwdot[H]/dT
  J[164] += dqdT; // dwdot[CO]/dT
  J[166] -= dqdT; // dwdot[HCO]/dT

  // reaction 28: HCO + O <=> CO2 + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[10];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[5] * sc[9];
  Kc = exp(g_RT[2] - g_RT[5] - g_RT[9] + g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[10]) + (h_RT[5] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q;  // O
  wdot[5] += q;  // H
  wdot[9] += q;  // CO2
  wdot[10] -= q; // HCO
  // d()/d[O]
  dqdci = +k_f * sc[10];
  J[28] -= dqdci; // dwdot[O]/d[O]
  J[31] += dqdci; // dwdot[H]/d[O]
  J[35] += dqdci; // dwdot[CO2]/d[O]
  J[36] -= dqdci; // dwdot[HCO]/d[O]
  // d()/d[H]
  dqdci = -k_r * sc[9];
  J[67] -= dqdci; // dwdot[O]/d[H]
  J[70] += dqdci; // dwdot[H]/d[H]
  J[74] += dqdci; // dwdot[CO2]/d[H]
  J[75] -= dqdci; // dwdot[HCO]/d[H]
  // d()/d[CO2]
  dqdci = -k_r * sc[5];
  J[119] -= dqdci; // dwdot[O]/d[CO2]
  J[122] += dqdci; // dwdot[H]/d[CO2]
  J[126] += dqdci; // dwdot[CO2]/d[CO2]
  J[127] -= dqdci; // dwdot[HCO]/d[CO2]
  // d()/d[HCO]
  dqdci = +k_f * sc[2];
  J[132] -= dqdci; // dwdot[O]/d[HCO]
  J[135] += dqdci; // dwdot[H]/d[HCO]
  J[139] += dqdci; // dwdot[CO2]/d[HCO]
  J[140] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[158] -= dqdT; // dwdot[O]/dT
  J[161] += dqdT; // dwdot[H]/dT
  J[165] += dqdT; // dwdot[CO2]/dT
  J[166] -= dqdT; // dwdot[HCO]/dT

  amrex::Real c_R[12], dcRdT[12], e_RT[12];
  amrex::Real* eh_RT;
  if (HP == 1) {
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
  for (int k = 0; k < 12; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[156 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 12; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 12; ++m) {
      dehmixdc += eh_RT[m] * J[k * 13 + m];
    }
    J[k * 13 + 12] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[168] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
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
  amrex::Real* J, const amrex::Real* sc, const amrex::Real T, const int consP)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  for (int i = 0; i < 169; i++) {
    J[i] = 0.0;
  }

  amrex::Real wdot[12];
  for (auto& val : wdot) {
    val = 0.0;
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
  for (int k = 0; k < 12; ++k) {
    mixture += sc[k];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[12];
  gibbs(g_RT, tc);

  // compute the species enthalpy
  amrex::Real h_RT[12];
  speciesEnthalpy(h_RT, tc);

  amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  amrex::Real dqdci, dcdc_fac, dqdc[12];
  amrex::Real Pr, fPr, F, k_0, logPr;
  amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
  amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
  amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const amrex::Real ln10 = log(10.0);
  const amrex::Real log10e = 1.0 / log(10.0);
  // reaction 8: H + O2 (+M) <=> HO2 (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + sc[0] + (0.78 - 1) * sc[1] + (11 - 1) * sc[4] +
          (1.9 - 1) * sc[8] + (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[1] * sc[5];
  k_f = 1475000 * exp(0.6 * tc[0]);
  dlnkfdT = 0.6 * invT;
  // pressure-fall-off
  k_0 = 6.366e+20 * exp(-1.72 * tc[0] - (264.088106214317) * invT);
  Pr = 1e-12 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -1.72 * invT + (264.088106214317) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.8) * exp(-T / 1e-30);
  Fcent2 = 0.8 * exp(-T / 1e+30);
  Fcent3 = 0.;
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT = log10e / Fcent * (-Fcent1 / 1e-30 - Fcent2 / 1e+30);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[6];
  Kc = refCinv * exp(g_RT[1] + g_RT[5] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[6]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[1] -= q; // O2
  wdot[5] -= q; // H
  wdot[6] += q; // HO2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2 - 1) * dcdc_fac;
    J[1] -= dqdci; // dwdot[O2]/d[H2]
    J[5] -= dqdci; // dwdot[H]/d[H2]
    J[6] += dqdci; // dwdot[HO2]/d[H2]
    // d()/d[O2]
    dqdci = (0.78 - 1) * dcdc_fac + k_f * sc[5];
    J[14] -= dqdci; // dwdot[O2]/d[O2]
    J[18] -= dqdci; // dwdot[H]/d[O2]
    J[19] += dqdci; // dwdot[HO2]/d[O2]
    // d()/d[H2O]
    dqdci = (11 - 1) * dcdc_fac;
    J[53] -= dqdci; // dwdot[O2]/d[H2O]
    J[57] -= dqdci; // dwdot[H]/d[H2O]
    J[58] += dqdci; // dwdot[HO2]/d[H2O]
    // d()/d[H]
    dqdci = +k_f * sc[1];
    J[66] -= dqdci; // dwdot[O2]/d[H]
    J[70] -= dqdci; // dwdot[H]/d[H]
    J[71] += dqdci; // dwdot[HO2]/d[H]
    // d()/d[HO2]
    dqdci = -k_r;
    J[79] -= dqdci; // dwdot[O2]/d[HO2]
    J[83] -= dqdci; // dwdot[H]/d[HO2]
    J[84] += dqdci; // dwdot[HO2]/d[HO2]
    // d()/d[CO]
    dqdci = (1.9 - 1) * dcdc_fac;
    J[105] -= dqdci; // dwdot[O2]/d[CO]
    J[109] -= dqdci; // dwdot[H]/d[CO]
    J[110] += dqdci; // dwdot[HO2]/d[CO]
    // d()/d[CO2]
    dqdci = (3.8 - 1) * dcdc_fac;
    J[118] -= dqdci; // dwdot[O2]/d[CO2]
    J[122] -= dqdci; // dwdot[H]/d[CO2]
    J[123] += dqdci; // dwdot[HO2]/d[CO2]
  } else {
    dqdc[0] = 2 * dcdc_fac;
    dqdc[1] = 0.78 * dcdc_fac + k_f * sc[5];
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = 11 * dcdc_fac;
    dqdc[5] = dcdc_fac + k_f * sc[1];
    dqdc[6] = dcdc_fac - k_r;
    dqdc[7] = dcdc_fac;
    dqdc[8] = 1.9 * dcdc_fac;
    dqdc[9] = 3.8 * dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    for (int k = 0; k < 12; k++) {
      J[13 * k + 1] -= dqdc[k];
      J[13 * k + 5] -= dqdc[k];
      J[13 * k + 6] += dqdc[k];
    }
  }
  J[157] -= dqdT; // dwdot[O2]/dT
  J[161] -= dqdT; // dwdot[H]/dT
  J[162] += dqdT; // dwdot[HO2]/dT

  // reaction 15: H2O2 (+M) <=> 2 OH (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (12 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[7];
  k_f = 295100000000000 * exp(-(24370.7831249226) * invT);
  dlnkfdT = (24370.7831249226) * invT2;
  // pressure-fall-off
  k_0 = 1.202e+17 * exp(-(22896.3582941147) * invT);
  Pr = 1e-6 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = (22896.3582941147) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.5) * exp(-T / 1e-30);
  Fcent2 = 0.5 * exp(-T / 1e+30);
  Fcent3 = 0.;
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT = log10e / Fcent * (-Fcent1 / 1e-30 - Fcent2 / 1e+30);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = (sc[3] * sc[3]);
  Kc = refC * exp(-2.000000 * g_RT[3] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[7]) + (2.000000 * h_RT[3]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[3] += 2 * q; // OH
  wdot[7] -= q;     // H2O2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2.5 - 1) * dcdc_fac;
    J[3] += 2 * dqdci; // dwdot[OH]/d[H2]
    J[7] -= dqdci;     // dwdot[H2O2]/d[H2]
    // d()/d[OH]
    dqdci = -k_r * 2.000000 * sc[3];
    J[42] += 2 * dqdci; // dwdot[OH]/d[OH]
    J[46] -= dqdci;     // dwdot[H2O2]/d[OH]
    // d()/d[H2O]
    dqdci = (12 - 1) * dcdc_fac;
    J[55] += 2 * dqdci; // dwdot[OH]/d[H2O]
    J[59] -= dqdci;     // dwdot[H2O2]/d[H2O]
    // d()/d[H2O2]
    dqdci = +k_f;
    J[94] += 2 * dqdci; // dwdot[OH]/d[H2O2]
    J[98] -= dqdci;     // dwdot[H2O2]/d[H2O2]
    // d()/d[CO]
    dqdci = (1.9 - 1) * dcdc_fac;
    J[107] += 2 * dqdci; // dwdot[OH]/d[CO]
    J[111] -= dqdci;     // dwdot[H2O2]/d[CO]
    // d()/d[CO2]
    dqdci = (3.8 - 1) * dcdc_fac;
    J[120] += 2 * dqdci; // dwdot[OH]/d[CO2]
    J[124] -= dqdci;     // dwdot[H2O2]/d[CO2]
  } else {
    dqdc[0] = 2.5 * dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac - k_r * 2.000000 * sc[3];
    dqdc[4] = 12 * dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac + k_f;
    dqdc[8] = 1.9 * dcdc_fac;
    dqdc[9] = 3.8 * dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    for (int k = 0; k < 12; k++) {
      J[13 * k + 3] += 2 * dqdc[k];
      J[13 * k + 7] -= dqdc[k];
    }
  }
  J[159] += 2 * dqdT; // dwdot[OH]/dT
  J[163] -= dqdT;     // dwdot[H2O2]/dT

  // reaction 21: CO + O (+M) <=> CO2 (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (12 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[2] * sc[8];
  k_f = 18000 * exp(-(1199.66853127845) * invT);
  dlnkfdT = (1199.66853127845) * invT2;
  // pressure-fall-off
  k_0 = 1.55e+24 * exp(-2.79 * tc[0] - (2108.98104638758) * invT);
  Pr = 1e-12 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -2.79 * invT + (2108.98104638758) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Lindemann form
  F = 1.0;
  dlogFdlogPr = 0.0;
  dlogFdT = 0.0;
  // reverse
  phi_r = sc[9];
  Kc = refCinv * exp(g_RT[2] + g_RT[8] - g_RT[9]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[9]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[2] -= q; // O
  wdot[8] -= q; // CO
  wdot[9] += q; // CO2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2.5 - 1) * dcdc_fac;
    J[2] -= dqdci; // dwdot[O]/d[H2]
    J[8] -= dqdci; // dwdot[CO]/d[H2]
    J[9] += dqdci; // dwdot[CO2]/d[H2]
    // d()/d[O]
    dqdci = +k_f * sc[8];
    J[28] -= dqdci; // dwdot[O]/d[O]
    J[34] -= dqdci; // dwdot[CO]/d[O]
    J[35] += dqdci; // dwdot[CO2]/d[O]
    // d()/d[H2O]
    dqdci = (12 - 1) * dcdc_fac;
    J[54] -= dqdci; // dwdot[O]/d[H2O]
    J[60] -= dqdci; // dwdot[CO]/d[H2O]
    J[61] += dqdci; // dwdot[CO2]/d[H2O]
    // d()/d[CO]
    dqdci = (1.9 - 1) * dcdc_fac + k_f * sc[2];
    J[106] -= dqdci; // dwdot[O]/d[CO]
    J[112] -= dqdci; // dwdot[CO]/d[CO]
    J[113] += dqdci; // dwdot[CO2]/d[CO]
    // d()/d[CO2]
    dqdci = (3.8 - 1) * dcdc_fac - k_r;
    J[119] -= dqdci; // dwdot[O]/d[CO2]
    J[125] -= dqdci; // dwdot[CO]/d[CO2]
    J[126] += dqdci; // dwdot[CO2]/d[CO2]
  } else {
    dqdc[0] = 2.5 * dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac + k_f * sc[8];
    dqdc[3] = dcdc_fac;
    dqdc[4] = 12 * dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = 1.9 * dcdc_fac + k_f * sc[2];
    dqdc[9] = 3.8 * dcdc_fac - k_r;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    for (int k = 0; k < 12; k++) {
      J[13 * k + 2] -= dqdc[k];
      J[13 * k + 8] -= dqdc[k];
      J[13 * k + 9] += dqdc[k];
    }
  }
  J[158] -= dqdT; // dwdot[O]/dT
  J[164] -= dqdT; // dwdot[CO]/dT
  J[165] += dqdT; // dwdot[CO2]/dT

  // reaction 4: H2 + M <=> 2 H + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (12 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[0];
  k_f = 45770000000000 * exp(-1.4 * tc[0] - (52525.7555766966) * invT);
  dlnkfdT = -1.4 * invT + (52525.7555766966) * invT2;
  // reverse
  phi_r = (sc[5] * sc[5]);
  Kc = refC * exp(g_RT[0] - 2.000000 * g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0]) + (2.000000 * h_RT[5]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] -= q;     // H2
  wdot[5] += 2 * q; // H
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2.5 - 1) * q_nocor + k_f;
    J[0] -= dqdci;     // dwdot[H2]/d[H2]
    J[5] += 2 * dqdci; // dwdot[H]/d[H2]
    // d()/d[H2O]
    dqdci = (12 - 1) * q_nocor;
    J[52] -= dqdci;     // dwdot[H2]/d[H2O]
    J[57] += 2 * dqdci; // dwdot[H]/d[H2O]
    // d()/d[H]
    dqdci = -k_r * 2.000000 * sc[5];
    J[65] -= dqdci;     // dwdot[H2]/d[H]
    J[70] += 2 * dqdci; // dwdot[H]/d[H]
    // d()/d[CO]
    dqdci = (1.9 - 1) * q_nocor;
    J[104] -= dqdci;     // dwdot[H2]/d[CO]
    J[109] += 2 * dqdci; // dwdot[H]/d[CO]
    // d()/d[CO2]
    dqdci = (3.8 - 1) * q_nocor;
    J[117] -= dqdci;     // dwdot[H2]/d[CO2]
    J[122] += 2 * dqdci; // dwdot[H]/d[CO2]
  } else {
    dqdc[0] = 2.5 * q_nocor + k_f;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = 12 * q_nocor;
    dqdc[5] = q_nocor - k_r * 2.000000 * sc[5];
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = 1.9 * q_nocor;
    dqdc[9] = 3.8 * q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    for (int k = 0; k < 12; k++) {
      J[13 * k + 0] -= dqdc[k];
      J[13 * k + 5] += 2 * dqdc[k];
    }
  }
  J[156] -= dqdT;     // dwdot[H2]/dT
  J[161] += 2 * dqdT; // dwdot[H]/dT

  // reaction 5: 2 O + M <=> O2 + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (12 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = (sc[2] * sc[2]);
  k_f = 6165 * exp(-0.5 * tc[0]);
  dlnkfdT = -0.5 * invT;
  // reverse
  phi_r = sc[1];
  Kc = refCinv * exp(-g_RT[1] + 2.000000 * g_RT[2]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[2]) + (h_RT[1]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;     // O2
  wdot[2] -= 2 * q; // O
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2.5 - 1) * q_nocor;
    J[1] += dqdci;      // dwdot[O2]/d[H2]
    J[2] += -2 * dqdci; // dwdot[O]/d[H2]
    // d()/d[O2]
    dqdci = -k_r;
    J[14] += dqdci;      // dwdot[O2]/d[O2]
    J[15] += -2 * dqdci; // dwdot[O]/d[O2]
    // d()/d[O]
    dqdci = +k_f * 2.000000 * sc[2];
    J[27] += dqdci;      // dwdot[O2]/d[O]
    J[28] += -2 * dqdci; // dwdot[O]/d[O]
    // d()/d[H2O]
    dqdci = (12 - 1) * q_nocor;
    J[53] += dqdci;      // dwdot[O2]/d[H2O]
    J[54] += -2 * dqdci; // dwdot[O]/d[H2O]
    // d()/d[CO]
    dqdci = (1.9 - 1) * q_nocor;
    J[105] += dqdci;      // dwdot[O2]/d[CO]
    J[106] += -2 * dqdci; // dwdot[O]/d[CO]
    // d()/d[CO2]
    dqdci = (3.8 - 1) * q_nocor;
    J[118] += dqdci;      // dwdot[O2]/d[CO2]
    J[119] += -2 * dqdci; // dwdot[O]/d[CO2]
  } else {
    dqdc[0] = 2.5 * q_nocor;
    dqdc[1] = q_nocor - k_r;
    dqdc[2] = q_nocor + k_f * 2.000000 * sc[2];
    dqdc[3] = q_nocor;
    dqdc[4] = 12 * q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = 1.9 * q_nocor;
    dqdc[9] = 3.8 * q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    for (int k = 0; k < 12; k++) {
      J[13 * k + 1] += dqdc[k];
      J[13 * k + 2] += -2 * dqdc[k];
    }
  }
  J[157] += dqdT;      // dwdot[O2]/dT
  J[158] += -2 * dqdT; // dwdot[O]/dT

  // reaction 6: H + O + M <=> OH + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (12 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[2] * sc[5];
  k_f = 4714000 * exp(-1 * tc[0]);
  dlnkfdT = -1 * invT;
  // reverse
  phi_r = sc[3];
  Kc = refCinv * exp(g_RT[2] - g_RT[3] + g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[3]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // O
  wdot[3] += q; // OH
  wdot[5] -= q; // H
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2.5 - 1) * q_nocor;
    J[2] -= dqdci; // dwdot[O]/d[H2]
    J[3] += dqdci; // dwdot[OH]/d[H2]
    J[5] -= dqdci; // dwdot[H]/d[H2]
    // d()/d[O]
    dqdci = +k_f * sc[5];
    J[28] -= dqdci; // dwdot[O]/d[O]
    J[29] += dqdci; // dwdot[OH]/d[O]
    J[31] -= dqdci; // dwdot[H]/d[O]
    // d()/d[OH]
    dqdci = -k_r;
    J[41] -= dqdci; // dwdot[O]/d[OH]
    J[42] += dqdci; // dwdot[OH]/d[OH]
    J[44] -= dqdci; // dwdot[H]/d[OH]
    // d()/d[H2O]
    dqdci = (12 - 1) * q_nocor;
    J[54] -= dqdci; // dwdot[O]/d[H2O]
    J[55] += dqdci; // dwdot[OH]/d[H2O]
    J[57] -= dqdci; // dwdot[H]/d[H2O]
    // d()/d[H]
    dqdci = +k_f * sc[2];
    J[67] -= dqdci; // dwdot[O]/d[H]
    J[68] += dqdci; // dwdot[OH]/d[H]
    J[70] -= dqdci; // dwdot[H]/d[H]
    // d()/d[CO]
    dqdci = (1.9 - 1) * q_nocor;
    J[106] -= dqdci; // dwdot[O]/d[CO]
    J[107] += dqdci; // dwdot[OH]/d[CO]
    J[109] -= dqdci; // dwdot[H]/d[CO]
    // d()/d[CO2]
    dqdci = (3.8 - 1) * q_nocor;
    J[119] -= dqdci; // dwdot[O]/d[CO2]
    J[120] += dqdci; // dwdot[OH]/d[CO2]
    J[122] -= dqdci; // dwdot[H]/d[CO2]
  } else {
    dqdc[0] = 2.5 * q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor + k_f * sc[5];
    dqdc[3] = q_nocor - k_r;
    dqdc[4] = 12 * q_nocor;
    dqdc[5] = q_nocor + k_f * sc[2];
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = 1.9 * q_nocor;
    dqdc[9] = 3.8 * q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    for (int k = 0; k < 12; k++) {
      J[13 * k + 2] -= dqdc[k];
      J[13 * k + 3] += dqdc[k];
      J[13 * k + 5] -= dqdc[k];
    }
  }
  J[158] -= dqdT; // dwdot[O]/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[161] -= dqdT; // dwdot[H]/dT

  // reaction 7: H + OH + M <=> H2O + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (12 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[3] * sc[5];
  k_f = 38000000000 * exp(-2 * tc[0]);
  dlnkfdT = -2 * invT;
  // reverse
  phi_r = sc[4];
  Kc = refCinv * exp(g_RT[3] - g_RT[4] + g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[4]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q; // OH
  wdot[4] += q; // H2O
  wdot[5] -= q; // H
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2.5 - 1) * q_nocor;
    J[3] -= dqdci; // dwdot[OH]/d[H2]
    J[4] += dqdci; // dwdot[H2O]/d[H2]
    J[5] -= dqdci; // dwdot[H]/d[H2]
    // d()/d[OH]
    dqdci = +k_f * sc[5];
    J[42] -= dqdci; // dwdot[OH]/d[OH]
    J[43] += dqdci; // dwdot[H2O]/d[OH]
    J[44] -= dqdci; // dwdot[H]/d[OH]
    // d()/d[H2O]
    dqdci = (12 - 1) * q_nocor - k_r;
    J[55] -= dqdci; // dwdot[OH]/d[H2O]
    J[56] += dqdci; // dwdot[H2O]/d[H2O]
    J[57] -= dqdci; // dwdot[H]/d[H2O]
    // d()/d[H]
    dqdci = +k_f * sc[3];
    J[68] -= dqdci; // dwdot[OH]/d[H]
    J[69] += dqdci; // dwdot[H2O]/d[H]
    J[70] -= dqdci; // dwdot[H]/d[H]
    // d()/d[CO]
    dqdci = (1.9 - 1) * q_nocor;
    J[107] -= dqdci; // dwdot[OH]/d[CO]
    J[108] += dqdci; // dwdot[H2O]/d[CO]
    J[109] -= dqdci; // dwdot[H]/d[CO]
    // d()/d[CO2]
    dqdci = (3.8 - 1) * q_nocor;
    J[120] -= dqdci; // dwdot[OH]/d[CO2]
    J[121] += dqdci; // dwdot[H2O]/d[CO2]
    J[122] -= dqdci; // dwdot[H]/d[CO2]
  } else {
    dqdc[0] = 2.5 * q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor + k_f * sc[5];
    dqdc[4] = 12 * q_nocor - k_r;
    dqdc[5] = q_nocor + k_f * sc[3];
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = 1.9 * q_nocor;
    dqdc[9] = 3.8 * q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    for (int k = 0; k < 12; k++) {
      J[13 * k + 3] -= dqdc[k];
      J[13 * k + 4] += dqdc[k];
      J[13 * k + 5] -= dqdc[k];
    }
  }
  J[159] -= dqdT; // dwdot[OH]/dT
  J[160] += dqdT; // dwdot[H2O]/dT
  J[161] -= dqdT; // dwdot[H]/dT

  // reaction 25: HCO + M <=> CO + H + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[0] + (6 - 1) * sc[4] + (1.9 - 1) * sc[8] +
          (3.8 - 1) * sc[9];
  // forward
  phi_f = sc[10];
  k_f = 474850 * exp(0.659 * tc[0] - (7484.8446871794) * invT);
  dlnkfdT = 0.659 * invT + (7484.8446871794) * invT2;
  // reverse
  phi_r = sc[5] * sc[8];
  Kc = refC * exp(-g_RT[5] - g_RT[8] + g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[10]) + (h_RT[5] + h_RT[8]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H
  wdot[8] += q;  // CO
  wdot[10] -= q; // HCO
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2.5 - 1) * q_nocor;
    J[5] += dqdci;  // dwdot[H]/d[H2]
    J[8] += dqdci;  // dwdot[CO]/d[H2]
    J[10] -= dqdci; // dwdot[HCO]/d[H2]
    // d()/d[H2O]
    dqdci = (6 - 1) * q_nocor;
    J[57] += dqdci; // dwdot[H]/d[H2O]
    J[60] += dqdci; // dwdot[CO]/d[H2O]
    J[62] -= dqdci; // dwdot[HCO]/d[H2O]
    // d()/d[H]
    dqdci = -k_r * sc[8];
    J[70] += dqdci; // dwdot[H]/d[H]
    J[73] += dqdci; // dwdot[CO]/d[H]
    J[75] -= dqdci; // dwdot[HCO]/d[H]
    // d()/d[CO]
    dqdci = (1.9 - 1) * q_nocor - k_r * sc[5];
    J[109] += dqdci; // dwdot[H]/d[CO]
    J[112] += dqdci; // dwdot[CO]/d[CO]
    J[114] -= dqdci; // dwdot[HCO]/d[CO]
    // d()/d[CO2]
    dqdci = (3.8 - 1) * q_nocor;
    J[122] += dqdci; // dwdot[H]/d[CO2]
    J[125] += dqdci; // dwdot[CO]/d[CO2]
    J[127] -= dqdci; // dwdot[HCO]/d[CO2]
    // d()/d[HCO]
    dqdci = +k_f;
    J[135] += dqdci; // dwdot[H]/d[HCO]
    J[138] += dqdci; // dwdot[CO]/d[HCO]
    J[140] -= dqdci; // dwdot[HCO]/d[HCO]
  } else {
    dqdc[0] = 2.5 * q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = 6 * q_nocor;
    dqdc[5] = q_nocor - k_r * sc[8];
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = 1.9 * q_nocor - k_r * sc[5];
    dqdc[9] = 3.8 * q_nocor;
    dqdc[10] = q_nocor + k_f;
    dqdc[11] = q_nocor;
    for (int k = 0; k < 12; k++) {
      J[13 * k + 5] += dqdc[k];
      J[13 * k + 8] += dqdc[k];
      J[13 * k + 10] -= dqdc[k];
    }
  }
  J[161] += dqdT; // dwdot[H]/dT
  J[164] += dqdT; // dwdot[CO]/dT
  J[166] -= dqdT; // dwdot[HCO]/dT

  // reaction 0: H + O2 <=> O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[5];
  k_f = 3547000000 * exp(-0.406 * tc[0] - (8352.89343569254) * invT);
  dlnkfdT = -0.406 * invT + (8352.89343569254) * invT2;
  // reverse
  phi_r = sc[2] * sc[3];
  Kc = exp(g_RT[1] - g_RT[2] - g_RT[3] + g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[2] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q; // O2
  wdot[2] += q; // O
  wdot[3] += q; // OH
  wdot[5] -= q; // H
  // d()/d[O2]
  dqdci = +k_f * sc[5];
  J[14] -= dqdci; // dwdot[O2]/d[O2]
  J[15] += dqdci; // dwdot[O]/d[O2]
  J[16] += dqdci; // dwdot[OH]/d[O2]
  J[18] -= dqdci; // dwdot[H]/d[O2]
  // d()/d[O]
  dqdci = -k_r * sc[3];
  J[27] -= dqdci; // dwdot[O2]/d[O]
  J[28] += dqdci; // dwdot[O]/d[O]
  J[29] += dqdci; // dwdot[OH]/d[O]
  J[31] -= dqdci; // dwdot[H]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[2];
  J[40] -= dqdci; // dwdot[O2]/d[OH]
  J[41] += dqdci; // dwdot[O]/d[OH]
  J[42] += dqdci; // dwdot[OH]/d[OH]
  J[44] -= dqdci; // dwdot[H]/d[OH]
  // d()/d[H]
  dqdci = +k_f * sc[1];
  J[66] -= dqdci; // dwdot[O2]/d[H]
  J[67] += dqdci; // dwdot[O]/d[H]
  J[68] += dqdci; // dwdot[OH]/d[H]
  J[70] -= dqdci; // dwdot[H]/d[H]
  // d()/dT
  J[157] -= dqdT; // dwdot[O2]/dT
  J[158] += dqdT; // dwdot[O]/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[161] -= dqdT; // dwdot[H]/dT

  // reaction 1: H2 + O <=> H + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[2];
  k_f = 0.0508 * exp(2.67 * tc[0] - (3165.23282791169) * invT);
  dlnkfdT = 2.67 * invT + (3165.23282791169) * invT2;
  // reverse
  phi_r = sc[3] * sc[5];
  Kc = exp(g_RT[0] + g_RT[2] - g_RT[3] - g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] -= q; // H2
  wdot[2] -= q; // O
  wdot[3] += q; // OH
  wdot[5] += q; // H
  // d()/d[H2]
  dqdci = +k_f * sc[2];
  J[0] -= dqdci; // dwdot[H2]/d[H2]
  J[2] -= dqdci; // dwdot[O]/d[H2]
  J[3] += dqdci; // dwdot[OH]/d[H2]
  J[5] += dqdci; // dwdot[H]/d[H2]
  // d()/d[O]
  dqdci = +k_f * sc[0];
  J[26] -= dqdci; // dwdot[H2]/d[O]
  J[28] -= dqdci; // dwdot[O]/d[O]
  J[29] += dqdci; // dwdot[OH]/d[O]
  J[31] += dqdci; // dwdot[H]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[5];
  J[39] -= dqdci; // dwdot[H2]/d[OH]
  J[41] -= dqdci; // dwdot[O]/d[OH]
  J[42] += dqdci; // dwdot[OH]/d[OH]
  J[44] += dqdci; // dwdot[H]/d[OH]
  // d()/d[H]
  dqdci = -k_r * sc[3];
  J[65] -= dqdci; // dwdot[H2]/d[H]
  J[67] -= dqdci; // dwdot[O]/d[H]
  J[68] += dqdci; // dwdot[OH]/d[H]
  J[70] += dqdci; // dwdot[H]/d[H]
  // d()/dT
  J[156] -= dqdT; // dwdot[H2]/dT
  J[158] -= dqdT; // dwdot[O]/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[161] += dqdT; // dwdot[H]/dT

  // reaction 2: H2 + OH <=> H + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[3];
  k_f = 216 * exp(1.51 * tc[0] - (1726.03316371019) * invT);
  dlnkfdT = 1.51 * invT + (1726.03316371019) * invT2;
  // reverse
  phi_r = sc[4] * sc[5];
  Kc = exp(g_RT[0] + g_RT[3] - g_RT[4] - g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[3]) + (h_RT[4] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] -= q; // H2
  wdot[3] -= q; // OH
  wdot[4] += q; // H2O
  wdot[5] += q; // H
  // d()/d[H2]
  dqdci = +k_f * sc[3];
  J[0] -= dqdci; // dwdot[H2]/d[H2]
  J[3] -= dqdci; // dwdot[OH]/d[H2]
  J[4] += dqdci; // dwdot[H2O]/d[H2]
  J[5] += dqdci; // dwdot[H]/d[H2]
  // d()/d[OH]
  dqdci = +k_f * sc[0];
  J[39] -= dqdci; // dwdot[H2]/d[OH]
  J[42] -= dqdci; // dwdot[OH]/d[OH]
  J[43] += dqdci; // dwdot[H2O]/d[OH]
  J[44] += dqdci; // dwdot[H]/d[OH]
  // d()/d[H2O]
  dqdci = -k_r * sc[5];
  J[52] -= dqdci; // dwdot[H2]/d[H2O]
  J[55] -= dqdci; // dwdot[OH]/d[H2O]
  J[56] += dqdci; // dwdot[H2O]/d[H2O]
  J[57] += dqdci; // dwdot[H]/d[H2O]
  // d()/d[H]
  dqdci = -k_r * sc[4];
  J[65] -= dqdci; // dwdot[H2]/d[H]
  J[68] -= dqdci; // dwdot[OH]/d[H]
  J[69] += dqdci; // dwdot[H2O]/d[H]
  J[70] += dqdci; // dwdot[H]/d[H]
  // d()/dT
  J[156] -= dqdT; // dwdot[H2]/dT
  J[159] -= dqdT; // dwdot[OH]/dT
  J[160] += dqdT; // dwdot[H2O]/dT
  J[161] += dqdT; // dwdot[H]/dT

  // reaction 3: H2O + O <=> 2 OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[4];
  k_f = 2.97 * exp(2.02 * tc[0] - (6743.10332178324) * invT);
  dlnkfdT = 2.02 * invT + (6743.10332178324) * invT2;
  // reverse
  phi_r = (sc[3] * sc[3]);
  Kc = exp(g_RT[2] - 2.000000 * g_RT[3] + g_RT[4]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2.000000 * h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q;     // O
  wdot[3] += 2 * q; // OH
  wdot[4] -= q;     // H2O
  // d()/d[O]
  dqdci = +k_f * sc[4];
  J[28] -= dqdci;     // dwdot[O]/d[O]
  J[29] += 2 * dqdci; // dwdot[OH]/d[O]
  J[30] -= dqdci;     // dwdot[H2O]/d[O]
  // d()/d[OH]
  dqdci = -k_r * 2.000000 * sc[3];
  J[41] -= dqdci;     // dwdot[O]/d[OH]
  J[42] += 2 * dqdci; // dwdot[OH]/d[OH]
  J[43] -= dqdci;     // dwdot[H2O]/d[OH]
  // d()/d[H2O]
  dqdci = +k_f * sc[2];
  J[54] -= dqdci;     // dwdot[O]/d[H2O]
  J[55] += 2 * dqdci; // dwdot[OH]/d[H2O]
  J[56] -= dqdci;     // dwdot[H2O]/d[H2O]
  // d()/dT
  J[158] -= dqdT;     // dwdot[O]/dT
  J[159] += 2 * dqdT; // dwdot[OH]/dT
  J[160] -= dqdT;     // dwdot[H2O]/dT

  // reaction 9: H + HO2 <=> H2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[6];
  k_f = 16600000 * exp(-(414.147315957284) * invT);
  dlnkfdT = (414.147315957284) * invT2;
  // reverse
  phi_r = sc[0] * sc[1];
  Kc = exp(-g_RT[0] - g_RT[1] + g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (h_RT[0] + h_RT[1]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] += q; // H2
  wdot[1] += q; // O2
  wdot[5] -= q; // H
  wdot[6] -= q; // HO2
  // d()/d[H2]
  dqdci = -k_r * sc[1];
  J[0] += dqdci; // dwdot[H2]/d[H2]
  J[1] += dqdci; // dwdot[O2]/d[H2]
  J[5] -= dqdci; // dwdot[H]/d[H2]
  J[6] -= dqdci; // dwdot[HO2]/d[H2]
  // d()/d[O2]
  dqdci = -k_r * sc[0];
  J[13] += dqdci; // dwdot[H2]/d[O2]
  J[14] += dqdci; // dwdot[O2]/d[O2]
  J[18] -= dqdci; // dwdot[H]/d[O2]
  J[19] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[H]
  dqdci = +k_f * sc[6];
  J[65] += dqdci; // dwdot[H2]/d[H]
  J[66] += dqdci; // dwdot[O2]/d[H]
  J[70] -= dqdci; // dwdot[H]/d[H]
  J[71] -= dqdci; // dwdot[HO2]/d[H]
  // d()/d[HO2]
  dqdci = +k_f * sc[5];
  J[78] += dqdci; // dwdot[H2]/d[HO2]
  J[79] += dqdci; // dwdot[O2]/d[HO2]
  J[83] -= dqdci; // dwdot[H]/d[HO2]
  J[84] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[156] += dqdT; // dwdot[H2]/dT
  J[157] += dqdT; // dwdot[O2]/dT
  J[161] -= dqdT; // dwdot[H]/dT
  J[162] -= dqdT; // dwdot[HO2]/dT

  // reaction 10: H + HO2 <=> 2 OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[6];
  k_f = 70790000 * exp(-(148.448916412392) * invT);
  dlnkfdT = (148.448916412392) * invT2;
  // reverse
  phi_r = (sc[3] * sc[3]);
  Kc = exp(-2.000000 * g_RT[3] + g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (2.000000 * h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += 2 * q; // OH
  wdot[5] -= q;     // H
  wdot[6] -= q;     // HO2
  // d()/d[OH]
  dqdci = -k_r * 2.000000 * sc[3];
  J[42] += 2 * dqdci; // dwdot[OH]/d[OH]
  J[44] -= dqdci;     // dwdot[H]/d[OH]
  J[45] -= dqdci;     // dwdot[HO2]/d[OH]
  // d()/d[H]
  dqdci = +k_f * sc[6];
  J[68] += 2 * dqdci; // dwdot[OH]/d[H]
  J[70] -= dqdci;     // dwdot[H]/d[H]
  J[71] -= dqdci;     // dwdot[HO2]/d[H]
  // d()/d[HO2]
  dqdci = +k_f * sc[5];
  J[81] += 2 * dqdci; // dwdot[OH]/d[HO2]
  J[83] -= dqdci;     // dwdot[H]/d[HO2]
  J[84] -= dqdci;     // dwdot[HO2]/d[HO2]
  // d()/dT
  J[159] += 2 * dqdT; // dwdot[OH]/dT
  J[161] -= dqdT;     // dwdot[H]/dT
  J[162] -= dqdT;     // dwdot[HO2]/dT

  // reaction 11: HO2 + O <=> O2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[6];
  k_f = 32500000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[3];
  Kc = exp(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[1] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q; // O2
  wdot[2] -= q; // O
  wdot[3] += q; // OH
  wdot[6] -= q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[3];
  J[14] += dqdci; // dwdot[O2]/d[O2]
  J[15] -= dqdci; // dwdot[O]/d[O2]
  J[16] += dqdci; // dwdot[OH]/d[O2]
  J[19] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[O]
  dqdci = +k_f * sc[6];
  J[27] += dqdci; // dwdot[O2]/d[O]
  J[28] -= dqdci; // dwdot[O]/d[O]
  J[29] += dqdci; // dwdot[OH]/d[O]
  J[32] -= dqdci; // dwdot[HO2]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[1];
  J[40] += dqdci; // dwdot[O2]/d[OH]
  J[41] -= dqdci; // dwdot[O]/d[OH]
  J[42] += dqdci; // dwdot[OH]/d[OH]
  J[45] -= dqdci; // dwdot[HO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[2];
  J[79] += dqdci; // dwdot[O2]/d[HO2]
  J[80] -= dqdci; // dwdot[O]/d[HO2]
  J[81] += dqdci; // dwdot[OH]/d[HO2]
  J[84] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[157] += dqdT; // dwdot[O2]/dT
  J[158] -= dqdT; // dwdot[O]/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[162] -= dqdT; // dwdot[HO2]/dT

  // reaction 12: HO2 + OH <=> H2O + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[6];
  k_f = 28900000 * exp(-(-250.098682904946) * invT);
  dlnkfdT = (-250.098682904946) * invT2;
  // reverse
  phi_r = sc[1] * sc[4];
  Kc = exp(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[1] + h_RT[4]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q; // O2
  wdot[3] -= q; // OH
  wdot[4] += q; // H2O
  wdot[6] -= q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[4];
  J[14] += dqdci; // dwdot[O2]/d[O2]
  J[16] -= dqdci; // dwdot[OH]/d[O2]
  J[17] += dqdci; // dwdot[H2O]/d[O2]
  J[19] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[OH]
  dqdci = +k_f * sc[6];
  J[40] += dqdci; // dwdot[O2]/d[OH]
  J[42] -= dqdci; // dwdot[OH]/d[OH]
  J[43] += dqdci; // dwdot[H2O]/d[OH]
  J[45] -= dqdci; // dwdot[HO2]/d[OH]
  // d()/d[H2O]
  dqdci = -k_r * sc[1];
  J[53] += dqdci; // dwdot[O2]/d[H2O]
  J[55] -= dqdci; // dwdot[OH]/d[H2O]
  J[56] += dqdci; // dwdot[H2O]/d[H2O]
  J[58] -= dqdci; // dwdot[HO2]/d[H2O]
  // d()/d[HO2]
  dqdci = +k_f * sc[3];
  J[79] += dqdci; // dwdot[O2]/d[HO2]
  J[81] -= dqdci; // dwdot[OH]/d[HO2]
  J[82] += dqdci; // dwdot[H2O]/d[HO2]
  J[84] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[157] += dqdT; // dwdot[O2]/dT
  J[159] -= dqdT; // dwdot[OH]/dT
  J[160] += dqdT; // dwdot[H2O]/dT
  J[162] -= dqdT; // dwdot[HO2]/dT

  // reaction 13: 2 HO2 <=> H2O2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[6] * sc[6]);
  k_f = 420000000 * exp(-(6029.54208967215) * invT);
  dlnkfdT = (6029.54208967215) * invT2;
  // reverse
  phi_r = sc[1] * sc[7];
  Kc = exp(-g_RT[1] + 2.000000 * g_RT[6] - g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[6]) + (h_RT[1] + h_RT[7]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;     // O2
  wdot[6] -= 2 * q; // HO2
  wdot[7] += q;     // H2O2
  // d()/d[O2]
  dqdci = -k_r * sc[7];
  J[14] += dqdci;      // dwdot[O2]/d[O2]
  J[19] += -2 * dqdci; // dwdot[HO2]/d[O2]
  J[20] += dqdci;      // dwdot[H2O2]/d[O2]
  // d()/d[HO2]
  dqdci = +k_f * 2.000000 * sc[6];
  J[79] += dqdci;      // dwdot[O2]/d[HO2]
  J[84] += -2 * dqdci; // dwdot[HO2]/d[HO2]
  J[85] += dqdci;      // dwdot[H2O2]/d[HO2]
  // d()/d[H2O2]
  dqdci = -k_r * sc[1];
  J[92] += dqdci;      // dwdot[O2]/d[H2O2]
  J[97] += -2 * dqdci; // dwdot[HO2]/d[H2O2]
  J[98] += dqdci;      // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[157] += dqdT;      // dwdot[O2]/dT
  J[162] += -2 * dqdT; // dwdot[HO2]/dT
  J[163] += dqdT;      // dwdot[H2O2]/dT

  // reaction 14: 2 HO2 <=> H2O2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[6] * sc[6]);
  k_f = 130000 * exp(-(-819.89091359563) * invT);
  dlnkfdT = (-819.89091359563) * invT2;
  // reverse
  phi_r = sc[1] * sc[7];
  Kc = exp(-g_RT[1] + 2.000000 * g_RT[6] - g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[6]) + (h_RT[1] + h_RT[7]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;     // O2
  wdot[6] -= 2 * q; // HO2
  wdot[7] += q;     // H2O2
  // d()/d[O2]
  dqdci = -k_r * sc[7];
  J[14] += dqdci;      // dwdot[O2]/d[O2]
  J[19] += -2 * dqdci; // dwdot[HO2]/d[O2]
  J[20] += dqdci;      // dwdot[H2O2]/d[O2]
  // d()/d[HO2]
  dqdci = +k_f * 2.000000 * sc[6];
  J[79] += dqdci;      // dwdot[O2]/d[HO2]
  J[84] += -2 * dqdci; // dwdot[HO2]/d[HO2]
  J[85] += dqdci;      // dwdot[H2O2]/d[HO2]
  // d()/d[H2O2]
  dqdci = -k_r * sc[1];
  J[92] += dqdci;      // dwdot[O2]/d[H2O2]
  J[97] += -2 * dqdci; // dwdot[HO2]/d[H2O2]
  J[98] += dqdci;      // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[157] += dqdT;      // dwdot[O2]/dT
  J[162] += -2 * dqdT; // dwdot[HO2]/dT
  J[163] += dqdT;      // dwdot[H2O2]/dT

  // reaction 16: H + H2O2 <=> H2O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[7];
  k_f = 24100000 * exp(-(1997.77016324474) * invT);
  dlnkfdT = (1997.77016324474) * invT2;
  // reverse
  phi_r = sc[3] * sc[4];
  Kc = exp(-g_RT[3] - g_RT[4] + g_RT[5] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[3] + h_RT[4]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q; // OH
  wdot[4] += q; // H2O
  wdot[5] -= q; // H
  wdot[7] -= q; // H2O2
  // d()/d[OH]
  dqdci = -k_r * sc[4];
  J[42] += dqdci; // dwdot[OH]/d[OH]
  J[43] += dqdci; // dwdot[H2O]/d[OH]
  J[44] -= dqdci; // dwdot[H]/d[OH]
  J[46] -= dqdci; // dwdot[H2O2]/d[OH]
  // d()/d[H2O]
  dqdci = -k_r * sc[3];
  J[55] += dqdci; // dwdot[OH]/d[H2O]
  J[56] += dqdci; // dwdot[H2O]/d[H2O]
  J[57] -= dqdci; // dwdot[H]/d[H2O]
  J[59] -= dqdci; // dwdot[H2O2]/d[H2O]
  // d()/d[H]
  dqdci = +k_f * sc[7];
  J[68] += dqdci; // dwdot[OH]/d[H]
  J[69] += dqdci; // dwdot[H2O]/d[H]
  J[70] -= dqdci; // dwdot[H]/d[H]
  J[72] -= dqdci; // dwdot[H2O2]/d[H]
  // d()/d[H2O2]
  dqdci = +k_f * sc[5];
  J[94] += dqdci; // dwdot[OH]/d[H2O2]
  J[95] += dqdci; // dwdot[H2O]/d[H2O2]
  J[96] -= dqdci; // dwdot[H]/d[H2O2]
  J[98] -= dqdci; // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[160] += dqdT; // dwdot[H2O]/dT
  J[161] -= dqdT; // dwdot[H]/dT
  J[163] -= dqdT; // dwdot[H2O2]/dT

  // reaction 17: H + H2O2 <=> H2 + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[7];
  k_f = 48200000 * exp(-(4000.57249314752) * invT);
  dlnkfdT = (4000.57249314752) * invT2;
  // reverse
  phi_r = sc[0] * sc[6];
  Kc = exp(-g_RT[0] + g_RT[5] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[0] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] += q; // H2
  wdot[5] -= q; // H
  wdot[6] += q; // HO2
  wdot[7] -= q; // H2O2
  // d()/d[H2]
  dqdci = -k_r * sc[6];
  J[0] += dqdci; // dwdot[H2]/d[H2]
  J[5] -= dqdci; // dwdot[H]/d[H2]
  J[6] += dqdci; // dwdot[HO2]/d[H2]
  J[7] -= dqdci; // dwdot[H2O2]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[7];
  J[65] += dqdci; // dwdot[H2]/d[H]
  J[70] -= dqdci; // dwdot[H]/d[H]
  J[71] += dqdci; // dwdot[HO2]/d[H]
  J[72] -= dqdci; // dwdot[H2O2]/d[H]
  // d()/d[HO2]
  dqdci = -k_r * sc[0];
  J[78] += dqdci; // dwdot[H2]/d[HO2]
  J[83] -= dqdci; // dwdot[H]/d[HO2]
  J[84] += dqdci; // dwdot[HO2]/d[HO2]
  J[85] -= dqdci; // dwdot[H2O2]/d[HO2]
  // d()/d[H2O2]
  dqdci = +k_f * sc[5];
  J[91] += dqdci; // dwdot[H2]/d[H2O2]
  J[96] -= dqdci; // dwdot[H]/d[H2O2]
  J[97] += dqdci; // dwdot[HO2]/d[H2O2]
  J[98] -= dqdci; // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[156] += dqdT; // dwdot[H2]/dT
  J[161] -= dqdT; // dwdot[H]/dT
  J[162] += dqdT; // dwdot[HO2]/dT
  J[163] -= dqdT; // dwdot[H2O2]/dT

  // reaction 18: H2O2 + O <=> HO2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[7];
  k_f = 9.55 * exp(2 * tc[0] - (1997.77016324474) * invT);
  dlnkfdT = 2 * invT + (1997.77016324474) * invT2;
  // reverse
  phi_r = sc[3] * sc[6];
  Kc = exp(g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[3] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // O
  wdot[3] += q; // OH
  wdot[6] += q; // HO2
  wdot[7] -= q; // H2O2
  // d()/d[O]
  dqdci = +k_f * sc[7];
  J[28] -= dqdci; // dwdot[O]/d[O]
  J[29] += dqdci; // dwdot[OH]/d[O]
  J[32] += dqdci; // dwdot[HO2]/d[O]
  J[33] -= dqdci; // dwdot[H2O2]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[6];
  J[41] -= dqdci; // dwdot[O]/d[OH]
  J[42] += dqdci; // dwdot[OH]/d[OH]
  J[45] += dqdci; // dwdot[HO2]/d[OH]
  J[46] -= dqdci; // dwdot[H2O2]/d[OH]
  // d()/d[HO2]
  dqdci = -k_r * sc[3];
  J[80] -= dqdci; // dwdot[O]/d[HO2]
  J[81] += dqdci; // dwdot[OH]/d[HO2]
  J[84] += dqdci; // dwdot[HO2]/d[HO2]
  J[85] -= dqdci; // dwdot[H2O2]/d[HO2]
  // d()/d[H2O2]
  dqdci = +k_f * sc[2];
  J[93] -= dqdci; // dwdot[O]/d[H2O2]
  J[94] += dqdci; // dwdot[OH]/d[H2O2]
  J[97] += dqdci; // dwdot[HO2]/d[H2O2]
  J[98] -= dqdci; // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[158] -= dqdT; // dwdot[O]/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[162] += dqdT; // dwdot[HO2]/dT
  J[163] -= dqdT; // dwdot[H2O2]/dT

  // reaction 19: H2O2 + OH <=> H2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[7];
  k_f = 1000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[4] * sc[6];
  Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q; // OH
  wdot[4] += q; // H2O
  wdot[6] += q; // HO2
  wdot[7] -= q; // H2O2
  // d()/d[OH]
  dqdci = +k_f * sc[7];
  J[42] -= dqdci; // dwdot[OH]/d[OH]
  J[43] += dqdci; // dwdot[H2O]/d[OH]
  J[45] += dqdci; // dwdot[HO2]/d[OH]
  J[46] -= dqdci; // dwdot[H2O2]/d[OH]
  // d()/d[H2O]
  dqdci = -k_r * sc[6];
  J[55] -= dqdci; // dwdot[OH]/d[H2O]
  J[56] += dqdci; // dwdot[H2O]/d[H2O]
  J[58] += dqdci; // dwdot[HO2]/d[H2O]
  J[59] -= dqdci; // dwdot[H2O2]/d[H2O]
  // d()/d[HO2]
  dqdci = -k_r * sc[4];
  J[81] -= dqdci; // dwdot[OH]/d[HO2]
  J[82] += dqdci; // dwdot[H2O]/d[HO2]
  J[84] += dqdci; // dwdot[HO2]/d[HO2]
  J[85] -= dqdci; // dwdot[H2O2]/d[HO2]
  // d()/d[H2O2]
  dqdci = +k_f * sc[3];
  J[94] -= dqdci; // dwdot[OH]/d[H2O2]
  J[95] += dqdci; // dwdot[H2O]/d[H2O2]
  J[97] += dqdci; // dwdot[HO2]/d[H2O2]
  J[98] -= dqdci; // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[159] -= dqdT; // dwdot[OH]/dT
  J[160] += dqdT; // dwdot[H2O]/dT
  J[162] += dqdT; // dwdot[HO2]/dT
  J[163] -= dqdT; // dwdot[H2O2]/dT

  // reaction 20: H2O2 + OH <=> H2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[7];
  k_f = 580000000 * exp(-(4809.24167509571) * invT);
  dlnkfdT = (4809.24167509571) * invT2;
  // reverse
  phi_r = sc[4] * sc[6];
  Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q; // OH
  wdot[4] += q; // H2O
  wdot[6] += q; // HO2
  wdot[7] -= q; // H2O2
  // d()/d[OH]
  dqdci = +k_f * sc[7];
  J[42] -= dqdci; // dwdot[OH]/d[OH]
  J[43] += dqdci; // dwdot[H2O]/d[OH]
  J[45] += dqdci; // dwdot[HO2]/d[OH]
  J[46] -= dqdci; // dwdot[H2O2]/d[OH]
  // d()/d[H2O]
  dqdci = -k_r * sc[6];
  J[55] -= dqdci; // dwdot[OH]/d[H2O]
  J[56] += dqdci; // dwdot[H2O]/d[H2O]
  J[58] += dqdci; // dwdot[HO2]/d[H2O]
  J[59] -= dqdci; // dwdot[H2O2]/d[H2O]
  // d()/d[HO2]
  dqdci = -k_r * sc[4];
  J[81] -= dqdci; // dwdot[OH]/d[HO2]
  J[82] += dqdci; // dwdot[H2O]/d[HO2]
  J[84] += dqdci; // dwdot[HO2]/d[HO2]
  J[85] -= dqdci; // dwdot[H2O2]/d[HO2]
  // d()/d[H2O2]
  dqdci = +k_f * sc[3];
  J[94] -= dqdci; // dwdot[OH]/d[H2O2]
  J[95] += dqdci; // dwdot[H2O]/d[H2O2]
  J[97] += dqdci; // dwdot[HO2]/d[H2O2]
  J[98] -= dqdci; // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[159] -= dqdT; // dwdot[OH]/dT
  J[160] += dqdT; // dwdot[H2O]/dT
  J[162] += dqdT; // dwdot[HO2]/dT
  J[163] -= dqdT; // dwdot[H2O2]/dT

  // reaction 22: CO + O2 <=> CO2 + O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[8];
  k_f = 2530000 * exp(-(24003.4349588851) * invT);
  dlnkfdT = (24003.4349588851) * invT2;
  // reverse
  phi_r = sc[2] * sc[9];
  Kc = exp(g_RT[1] - g_RT[2] + g_RT[8] - g_RT[9]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[8]) + (h_RT[2] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q; // O2
  wdot[2] += q; // O
  wdot[8] -= q; // CO
  wdot[9] += q; // CO2
  // d()/d[O2]
  dqdci = +k_f * sc[8];
  J[14] -= dqdci; // dwdot[O2]/d[O2]
  J[15] += dqdci; // dwdot[O]/d[O2]
  J[21] -= dqdci; // dwdot[CO]/d[O2]
  J[22] += dqdci; // dwdot[CO2]/d[O2]
  // d()/d[O]
  dqdci = -k_r * sc[9];
  J[27] -= dqdci; // dwdot[O2]/d[O]
  J[28] += dqdci; // dwdot[O]/d[O]
  J[34] -= dqdci; // dwdot[CO]/d[O]
  J[35] += dqdci; // dwdot[CO2]/d[O]
  // d()/d[CO]
  dqdci = +k_f * sc[1];
  J[105] -= dqdci; // dwdot[O2]/d[CO]
  J[106] += dqdci; // dwdot[O]/d[CO]
  J[112] -= dqdci; // dwdot[CO]/d[CO]
  J[113] += dqdci; // dwdot[CO2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[2];
  J[118] -= dqdci; // dwdot[O2]/d[CO2]
  J[119] += dqdci; // dwdot[O]/d[CO2]
  J[125] -= dqdci; // dwdot[CO]/d[CO2]
  J[126] += dqdci; // dwdot[CO2]/d[CO2]
  // d()/dT
  J[157] -= dqdT; // dwdot[O2]/dT
  J[158] += dqdT; // dwdot[O]/dT
  J[164] -= dqdT; // dwdot[CO]/dT
  J[165] += dqdT; // dwdot[CO2]/dT

  // reaction 23: CO + HO2 <=> CO2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[8];
  k_f = 30100000 * exp(-(11573.9833135086) * invT);
  dlnkfdT = (11573.9833135086) * invT2;
  // reverse
  phi_r = sc[3] * sc[9];
  Kc = exp(-g_RT[3] + g_RT[6] + g_RT[8] - g_RT[9]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[3] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q; // OH
  wdot[6] -= q; // HO2
  wdot[8] -= q; // CO
  wdot[9] += q; // CO2
  // d()/d[OH]
  dqdci = -k_r * sc[9];
  J[42] += dqdci; // dwdot[OH]/d[OH]
  J[45] -= dqdci; // dwdot[HO2]/d[OH]
  J[47] -= dqdci; // dwdot[CO]/d[OH]
  J[48] += dqdci; // dwdot[CO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[8];
  J[81] += dqdci; // dwdot[OH]/d[HO2]
  J[84] -= dqdci; // dwdot[HO2]/d[HO2]
  J[86] -= dqdci; // dwdot[CO]/d[HO2]
  J[87] += dqdci; // dwdot[CO2]/d[HO2]
  // d()/d[CO]
  dqdci = +k_f * sc[6];
  J[107] += dqdci; // dwdot[OH]/d[CO]
  J[110] -= dqdci; // dwdot[HO2]/d[CO]
  J[112] -= dqdci; // dwdot[CO]/d[CO]
  J[113] += dqdci; // dwdot[CO2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[3];
  J[120] += dqdci; // dwdot[OH]/d[CO2]
  J[123] -= dqdci; // dwdot[HO2]/d[CO2]
  J[125] -= dqdci; // dwdot[CO]/d[CO2]
  J[126] += dqdci; // dwdot[CO2]/d[CO2]
  // d()/dT
  J[159] += dqdT; // dwdot[OH]/dT
  J[162] -= dqdT; // dwdot[HO2]/dT
  J[164] -= dqdT; // dwdot[CO]/dT
  J[165] += dqdT; // dwdot[CO2]/dT

  // reaction 24: CO + OH <=> CO2 + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[8];
  k_f = 0.2229 * exp(1.89 * tc[0] - (-583.077150667929) * invT);
  dlnkfdT = 1.89 * invT + (-583.077150667929) * invT2;
  // reverse
  phi_r = sc[5] * sc[9];
  Kc = exp(g_RT[3] - g_RT[5] + g_RT[8] - g_RT[9]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[5] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q; // OH
  wdot[5] += q; // H
  wdot[8] -= q; // CO
  wdot[9] += q; // CO2
  // d()/d[OH]
  dqdci = +k_f * sc[8];
  J[42] -= dqdci; // dwdot[OH]/d[OH]
  J[44] += dqdci; // dwdot[H]/d[OH]
  J[47] -= dqdci; // dwdot[CO]/d[OH]
  J[48] += dqdci; // dwdot[CO2]/d[OH]
  // d()/d[H]
  dqdci = -k_r * sc[9];
  J[68] -= dqdci; // dwdot[OH]/d[H]
  J[70] += dqdci; // dwdot[H]/d[H]
  J[73] -= dqdci; // dwdot[CO]/d[H]
  J[74] += dqdci; // dwdot[CO2]/d[H]
  // d()/d[CO]
  dqdci = +k_f * sc[3];
  J[107] -= dqdci; // dwdot[OH]/d[CO]
  J[109] += dqdci; // dwdot[H]/d[CO]
  J[112] -= dqdci; // dwdot[CO]/d[CO]
  J[113] += dqdci; // dwdot[CO2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[5];
  J[120] -= dqdci; // dwdot[OH]/d[CO2]
  J[122] += dqdci; // dwdot[H]/d[CO2]
  J[125] -= dqdci; // dwdot[CO]/d[CO2]
  J[126] += dqdci; // dwdot[CO2]/d[CO2]
  // d()/dT
  J[159] -= dqdT; // dwdot[OH]/dT
  J[161] += dqdT; // dwdot[H]/dT
  J[164] -= dqdT; // dwdot[CO]/dT
  J[165] += dqdT; // dwdot[CO2]/dT

  // reaction 26: HCO + O2 <=> CO + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[10];
  k_f = 7580000 * exp(-(206.318832979935) * invT);
  dlnkfdT = (206.318832979935) * invT2;
  // reverse
  phi_r = sc[6] * sc[8];
  Kc = exp(g_RT[1] - g_RT[6] - g_RT[8] + g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[6] + h_RT[8]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q;  // O2
  wdot[6] += q;  // HO2
  wdot[8] += q;  // CO
  wdot[10] -= q; // HCO
  // d()/d[O2]
  dqdci = +k_f * sc[10];
  J[14] -= dqdci; // dwdot[O2]/d[O2]
  J[19] += dqdci; // dwdot[HO2]/d[O2]
  J[21] += dqdci; // dwdot[CO]/d[O2]
  J[23] -= dqdci; // dwdot[HCO]/d[O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[8];
  J[79] -= dqdci; // dwdot[O2]/d[HO2]
  J[84] += dqdci; // dwdot[HO2]/d[HO2]
  J[86] += dqdci; // dwdot[CO]/d[HO2]
  J[88] -= dqdci; // dwdot[HCO]/d[HO2]
  // d()/d[CO]
  dqdci = -k_r * sc[6];
  J[105] -= dqdci; // dwdot[O2]/d[CO]
  J[110] += dqdci; // dwdot[HO2]/d[CO]
  J[112] += dqdci; // dwdot[CO]/d[CO]
  J[114] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[1];
  J[131] -= dqdci; // dwdot[O2]/d[HCO]
  J[136] += dqdci; // dwdot[HO2]/d[HCO]
  J[138] += dqdci; // dwdot[CO]/d[HCO]
  J[140] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[157] -= dqdT; // dwdot[O2]/dT
  J[162] += dqdT; // dwdot[HO2]/dT
  J[164] += dqdT; // dwdot[CO]/dT
  J[166] -= dqdT; // dwdot[HCO]/dT

  // reaction 27: H + HCO <=> CO + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[10];
  k_f = 72300000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[0] * sc[8];
  Kc = exp(-g_RT[0] + g_RT[5] - g_RT[8] + g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[10]) + (h_RT[0] + h_RT[8]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] += q;  // H2
  wdot[5] -= q;  // H
  wdot[8] += q;  // CO
  wdot[10] -= q; // HCO
  // d()/d[H2]
  dqdci = -k_r * sc[8];
  J[0] += dqdci;  // dwdot[H2]/d[H2]
  J[5] -= dqdci;  // dwdot[H]/d[H2]
  J[8] += dqdci;  // dwdot[CO]/d[H2]
  J[10] -= dqdci; // dwdot[HCO]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[10];
  J[65] += dqdci; // dwdot[H2]/d[H]
  J[70] -= dqdci; // dwdot[H]/d[H]
  J[73] += dqdci; // dwdot[CO]/d[H]
  J[75] -= dqdci; // dwdot[HCO]/d[H]
  // d()/d[CO]
  dqdci = -k_r * sc[0];
  J[104] += dqdci; // dwdot[H2]/d[CO]
  J[109] -= dqdci; // dwdot[H]/d[CO]
  J[112] += dqdci; // dwdot[CO]/d[CO]
  J[114] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[5];
  J[130] += dqdci; // dwdot[H2]/d[HCO]
  J[135] -= dqdci; // dwdot[H]/d[HCO]
  J[138] += dqdci; // dwdot[CO]/d[HCO]
  J[140] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[156] += dqdT; // dwdot[H2]/dT
  J[161] -= dqdT; // dwdot[H]/dT
  J[164] += dqdT; // dwdot[CO]/dT
  J[166] -= dqdT; // dwdot[HCO]/dT

  // reaction 28: HCO + O <=> CO2 + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[10];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[5] * sc[9];
  Kc = exp(g_RT[2] - g_RT[5] - g_RT[9] + g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[10]) + (h_RT[5] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q;  // O
  wdot[5] += q;  // H
  wdot[9] += q;  // CO2
  wdot[10] -= q; // HCO
  // d()/d[O]
  dqdci = +k_f * sc[10];
  J[28] -= dqdci; // dwdot[O]/d[O]
  J[31] += dqdci; // dwdot[H]/d[O]
  J[35] += dqdci; // dwdot[CO2]/d[O]
  J[36] -= dqdci; // dwdot[HCO]/d[O]
  // d()/d[H]
  dqdci = -k_r * sc[9];
  J[67] -= dqdci; // dwdot[O]/d[H]
  J[70] += dqdci; // dwdot[H]/d[H]
  J[74] += dqdci; // dwdot[CO2]/d[H]
  J[75] -= dqdci; // dwdot[HCO]/d[H]
  // d()/d[CO2]
  dqdci = -k_r * sc[5];
  J[119] -= dqdci; // dwdot[O]/d[CO2]
  J[122] += dqdci; // dwdot[H]/d[CO2]
  J[126] += dqdci; // dwdot[CO2]/d[CO2]
  J[127] -= dqdci; // dwdot[HCO]/d[CO2]
  // d()/d[HCO]
  dqdci = +k_f * sc[2];
  J[132] -= dqdci; // dwdot[O]/d[HCO]
  J[135] += dqdci; // dwdot[H]/d[HCO]
  J[139] += dqdci; // dwdot[CO2]/d[HCO]
  J[140] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[158] -= dqdT; // dwdot[O]/dT
  J[161] += dqdT; // dwdot[H]/dT
  J[165] += dqdT; // dwdot[CO2]/dT
  J[166] -= dqdT; // dwdot[HCO]/dT

  amrex::Real c_R[12], dcRdT[12], e_RT[12];
  amrex::Real* eh_RT;
  if (consP == 1) {
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
  for (int k = 0; k < 12; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[156 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 12; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 12; ++m) {
      dehmixdc += eh_RT[m] * J[k * 13 + m];
    }
    J[k * 13 + 12] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[168] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
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
  *LENIMC = 50;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENRMC(int* LENRMC)
{
  *LENRMC = 3156;
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
  *NLITE = 2;
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
  WT[0] = 2.01600000E+00;
  WT[1] = 3.19980000E+01;
  WT[2] = 1.59990000E+01;
  WT[3] = 1.70070000E+01;
  WT[4] = 1.80150000E+01;
  WT[5] = 1.00800000E+00;
  WT[6] = 3.30060000E+01;
  WT[7] = 3.40140000E+01;
  WT[8] = 2.80100000E+01;
  WT[9] = 4.40090000E+01;
  WT[10] = 2.90180000E+01;
  WT[11] = 2.80140000E+01;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetEPS(amrex::Real* EPS)
{
  EPS[0] = 3.80000000E+01;
  EPS[1] = 1.07400000E+02;
  EPS[2] = 8.00000000E+01;
  EPS[3] = 8.00000000E+01;
  EPS[4] = 5.72400000E+02;
  EPS[5] = 1.45000000E+02;
  EPS[6] = 1.07400000E+02;
  EPS[7] = 1.07400000E+02;
  EPS[8] = 9.81000000E+01;
  EPS[9] = 2.44000000E+02;
  EPS[10] = 4.98000000E+02;
  EPS[11] = 9.75300000E+01;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetSIG(amrex::Real* SIG)
{
  SIG[0] = 2.92000000E+00;
  SIG[1] = 3.45800000E+00;
  SIG[2] = 2.75000000E+00;
  SIG[3] = 2.75000000E+00;
  SIG[4] = 2.60500000E+00;
  SIG[5] = 2.05000000E+00;
  SIG[6] = 3.45800000E+00;
  SIG[7] = 3.45800000E+00;
  SIG[8] = 3.65000000E+00;
  SIG[9] = 3.76300000E+00;
  SIG[10] = 3.59000000E+00;
  SIG[11] = 3.62100000E+00;
}

// the dipole moment in Debye
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetDIP(amrex::Real* DIP)
{
  DIP[0] = 0.00000000E+00;
  DIP[1] = 0.00000000E+00;
  DIP[2] = 0.00000000E+00;
  DIP[3] = 0.00000000E+00;
  DIP[4] = 1.84400000E+00;
  DIP[5] = 0.00000000E+00;
  DIP[6] = 0.00000000E+00;
  DIP[7] = 0.00000000E+00;
  DIP[8] = 0.00000000E+00;
  DIP[9] = 0.00000000E+00;
  DIP[10] = 0.00000000E+00;
  DIP[11] = 0.00000000E+00;
}

// the polarizability in cubic Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPOL(amrex::Real* POL)
{
  POL[0] = 7.90000000E-01;
  POL[1] = 1.60000000E+00;
  POL[2] = 0.00000000E+00;
  POL[3] = 0.00000000E+00;
  POL[4] = 0.00000000E+00;
  POL[5] = 0.00000000E+00;
  POL[6] = 0.00000000E+00;
  POL[7] = 0.00000000E+00;
  POL[8] = 1.95000000E+00;
  POL[9] = 2.65000000E+00;
  POL[10] = 0.00000000E+00;
  POL[11] = 1.76000000E+00;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetZROT(amrex::Real* ZROT)
{
  ZROT[0] = 2.80000000E+02;
  ZROT[1] = 3.80000000E+00;
  ZROT[2] = 0.00000000E+00;
  ZROT[3] = 0.00000000E+00;
  ZROT[4] = 4.00000000E+00;
  ZROT[5] = 0.00000000E+00;
  ZROT[6] = 1.00000000E+00;
  ZROT[7] = 3.80000000E+00;
  ZROT[8] = 1.80000000E+00;
  ZROT[9] = 2.10000000E+00;
  ZROT[10] = 0.00000000E+00;
  ZROT[11] = 4.00000000E+00;
}

// 0: monoatomic, 1: linear, 2: nonlinear
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLIN(int* NLIN)
{
  NLIN[0] = 1;
  NLIN[1] = 1;
  NLIN[2] = 0;
  NLIN[3] = 1;
  NLIN[4] = 2;
  NLIN[5] = 0;
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
  COFETA[0] = -1.37549291E+01;
  COFETA[1] = 9.65530587E-01;
  COFETA[2] = -4.45720114E-02;
  COFETA[3] = 2.05871810E-03;
  COFETA[4] = -1.68118998E+01;
  COFETA[5] = 2.52362554E+00;
  COFETA[6] = -2.49309128E-01;
  COFETA[7] = 1.10211025E-02;
  COFETA[8] = -1.48001711E+01;
  COFETA[9] = 1.79491990E+00;
  COFETA[10] = -1.54008440E-01;
  COFETA[11] = 6.86719439E-03;
  COFETA[12] = -1.47696217E+01;
  COFETA[13] = 1.79491990E+00;
  COFETA[14] = -1.54008440E-01;
  COFETA[15] = 6.86719439E-03;
  COFETA[16] = -1.17770995E+01;
  COFETA[17] = -8.26744826E-01;
  COFETA[18] = 3.39009392E-01;
  COFETA[19] = -2.00674472E-02;
  COFETA[20] = -1.98744352E+01;
  COFETA[21] = 3.41660514E+00;
  COFETA[22] = -3.63206306E-01;
  COFETA[23] = 1.58671021E-02;
  COFETA[24] = -1.67963918E+01;
  COFETA[25] = 2.52362554E+00;
  COFETA[26] = -2.49309128E-01;
  COFETA[27] = 1.10211025E-02;
  COFETA[28] = -1.67813504E+01;
  COFETA[29] = 2.52362554E+00;
  COFETA[30] = -2.49309128E-01;
  COFETA[31] = 1.10211025E-02;
  COFETA[32] = -1.63031343E+01;
  COFETA[33] = 2.26143219E+00;
  COFETA[34] = -2.15114671E-01;
  COFETA[35] = 9.53461976E-03;
  COFETA[36] = -2.36749638E+01;
  COFETA[37] = 4.99775518E+00;
  COFETA[38] = -5.52687718E-01;
  COFETA[39] = 2.34353338E-02;
  COFETA[40] = -2.11306886E+01;
  COFETA[41] = 3.26961843E+00;
  COFETA[42] = -2.51355092E-01;
  COFETA[43] = 7.35605058E-03;
  COFETA[44] = -1.62526676E+01;
  COFETA[45] = 2.24839597E+00;
  COFETA[46] = -2.13428438E-01;
  COFETA[47] = 9.46192413E-03;
}

// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFLAM(amrex::Real* COFLAM)
{
  COFLAM[0] = 1.11035309E+01;
  COFLAM[1] = -1.31883723E+00;
  COFLAM[2] = 2.44042222E-01;
  COFLAM[3] = -8.99835317E-03;
  COFLAM[4] = -2.51295522E+00;
  COFLAM[5] = 3.15165699E+00;
  COFLAM[6] = -3.10007687E-01;
  COFLAM[7] = 1.34522334E-02;
  COFLAM[8] = 1.98515155E+00;
  COFLAM[9] = 1.79491990E+00;
  COFLAM[10] = -1.54008440E-01;
  COFLAM[11] = 6.86719439E-03;
  COFLAM[12] = 1.60618776E+01;
  COFLAM[13] = -4.10626624E+00;
  COFLAM[14] = 6.63570999E-01;
  COFLAM[15] = -2.97906170E-02;
  COFLAM[16] = 2.21730553E+01;
  COFLAM[17] = -8.46935470E+00;
  COFLAM[18] = 1.46153792E+00;
  COFLAM[19] = -7.29502313E-02;
  COFLAM[20] = -3.24554543E-01;
  COFLAM[21] = 3.41660514E+00;
  COFLAM[22] = -3.63206306E-01;
  COFLAM[23] = 1.58671021E-02;
  COFLAM[24] = 5.56033781E-01;
  COFLAM[25] = 1.59073629E+00;
  COFLAM[26] = -5.28054537E-02;
  COFLAM[27] = 4.07605156E-04;
  COFLAM[28] = 1.48802076E+00;
  COFLAM[29] = 1.06176261E+00;
  COFLAM[30] = 5.72195029E-02;
  COFLAM[31] = -6.38391114E-03;
  COFLAM[32] = 9.94279641E+00;
  COFLAM[33] = -2.29161583E+00;
  COFLAM[34] = 4.74393158E-01;
  COFLAM[35] = -2.40686330E-02;
  COFLAM[36] = -1.21375517E+01;
  COFLAM[37] = 6.23624372E+00;
  COFLAM[38] = -6.22471355E-01;
  COFLAM[39] = 2.30613270E-02;
  COFLAM[40] = 1.70862641E+00;
  COFLAM[41] = -1.65304969E-01;
  COFLAM[42] = 3.30690188E-01;
  COFLAM[43] = -2.30097301E-02;
  COFLAM[44] = 1.15507226E+01;
  COFLAM[45] = -2.91453572E+00;
  COFLAM[46] = 5.55045264E-01;
  COFLAM[47] = -2.75173246E-02;
}

// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFD(amrex::Real* COFD)
{
  COFD[0] = -1.02395385E+01;
  COFD[1] = 2.15403244E+00;
  COFD[2] = -6.97480266E-02;
  COFD[3] = 3.23666871E-03;
  COFD[4] = -1.15797897E+01;
  COFD[5] = 2.43235504E+00;
  COFD[6] = -1.02890179E-01;
  COFD[7] = 4.52903603E-03;
  COFD[8] = -1.06250314E+01;
  COFD[9] = 2.15849701E+00;
  COFD[10] = -6.53886401E-02;
  COFD[11] = 2.81453370E-03;
  COFD[12] = -1.06283588E+01;
  COFD[13] = 2.15849701E+00;
  COFD[14] = -6.53886401E-02;
  COFD[15] = 2.81453370E-03;
  COFD[16] = -1.68759075E+01;
  COFD[17] = 4.49460341E+00;
  COFD[18] = -3.64766179E-01;
  COFD[19] = 1.56457173E-02;
  COFD[20] = -1.11808845E+01;
  COFD[21] = 2.66936727E+00;
  COFD[22] = -1.34411514E-01;
  COFD[23] = 5.92957488E-03;
  COFD[24] = -1.15806955E+01;
  COFD[25] = 2.43235504E+00;
  COFD[26] = -1.02890179E-01;
  COFD[27] = 4.52903603E-03;
  COFD[28] = -1.15815492E+01;
  COFD[29] = 2.43235504E+00;
  COFD[30] = -1.02890179E-01;
  COFD[31] = 4.52903603E-03;
  COFD[32] = -1.13541221E+01;
  COFD[33] = 2.31999438E+00;
  COFD[34] = -8.75064804E-02;
  COFD[35] = 3.82656365E-03;
  COFD[36] = -1.35545391E+01;
  COFD[37] = 3.13878730E+00;
  COFD[38] = -1.94980335E-01;
  COFD[39] = 8.53744486E-03;
  COFD[40] = -1.57161352E+01;
  COFD[41] = 3.96062263E+00;
  COFD[42] = -2.98964970E-01;
  COFD[43] = 1.29322565E-02;
  COFD[44] = -1.13253618E+01;
  COFD[45] = 2.31195095E+00;
  COFD[46] = -8.63988037E-02;
  COFD[47] = 3.77573452E-03;
  COFD[48] = -1.15797897E+01;
  COFD[49] = 2.43235504E+00;
  COFD[50] = -1.02890179E-01;
  COFD[51] = 4.52903603E-03;
  COFD[52] = -1.53110597E+01;
  COFD[53] = 3.37317428E+00;
  COFD[54] = -2.24900439E-01;
  COFD[55] = 9.81228151E-03;
  COFD[56] = -1.40864784E+01;
  COFD[57] = 3.07458927E+00;
  COFD[58] = -1.86899591E-01;
  COFD[59] = 8.19829781E-03;
  COFD[60] = -1.41066358E+01;
  COFD[61] = 3.07458927E+00;
  COFD[62] = -1.86899591E-01;
  COFD[63] = 8.19829781E-03;
  COFD[64] = -2.10639926E+01;
  COFD[65] = 5.50980703E+00;
  COFD[66] = -4.78335492E-01;
  COFD[67] = 1.98515434E-02;
  COFD[68] = -1.43713019E+01;
  COFD[69] = 3.70920439E+00;
  COFD[70] = -2.67274113E-01;
  COFD[71] = 1.15967481E-02;
  COFD[72] = -1.53187536E+01;
  COFD[73] = 3.37317428E+00;
  COFD[74] = -2.24900439E-01;
  COFD[75] = 9.81228151E-03;
  COFD[76] = -1.53261011E+01;
  COFD[77] = 3.37317428E+00;
  COFD[78] = -2.24900439E-01;
  COFD[79] = 9.81228151E-03;
  COFD[80] = -1.50371688E+01;
  COFD[81] = 3.26249588E+00;
  COFD[82] = -2.10658287E-01;
  COFD[83] = 9.20032462E-03;
  COFD[84] = -1.81197250E+01;
  COFD[85] = 4.33684042E+00;
  COFD[86] = -3.44981265E-01;
  COFD[87] = 1.48142449E-02;
  COFD[88] = -2.05107395E+01;
  COFD[89] = 5.21473296E+00;
  COFD[90] = -4.48646311E-01;
  COFD[91] = 1.89013813E-02;
  COFD[92] = -1.50096252E+01;
  COFD[93] = 3.25515933E+00;
  COFD[94] = -2.09710110E-01;
  COFD[95] = 9.15941830E-03;
  COFD[96] = -1.06250314E+01;
  COFD[97] = 2.15849701E+00;
  COFD[98] = -6.53886401E-02;
  COFD[99] = 2.81453370E-03;
  COFD[100] = -1.40864784E+01;
  COFD[101] = 3.07458927E+00;
  COFD[102] = -1.86899591E-01;
  COFD[103] = 8.19829781E-03;
  COFD[104] = -1.29877255E+01;
  COFD[105] = 2.80841511E+00;
  COFD[106] = -1.52629888E-01;
  COFD[107] = 6.72604927E-03;
  COFD[108] = -1.30027669E+01;
  COFD[109] = 2.80841511E+00;
  COFD[110] = -1.52629888E-01;
  COFD[111] = 6.72604927E-03;
  COFD[112] = -1.91096700E+01;
  COFD[113] = 5.02608697E+00;
  COFD[114] = -4.26959993E-01;
  COFD[115] = 1.80709910E-02;
  COFD[116] = -1.31860264E+01;
  COFD[117] = 3.38003453E+00;
  COFD[118] = -2.25783856E-01;
  COFD[119] = 9.85028660E-03;
  COFD[120] = -1.40915944E+01;
  COFD[121] = 3.07458927E+00;
  COFD[122] = -1.86899591E-01;
  COFD[123] = 8.19829781E-03;
  COFD[124] = -1.40964555E+01;
  COFD[125] = 3.07458927E+00;
  COFD[126] = -1.86899591E-01;
  COFD[127] = 8.19829781E-03;
  COFD[128] = -1.39007309E+01;
  COFD[129] = 2.99164244E+00;
  COFD[130] = -1.76293106E-01;
  COFD[131] = 7.74575100E-03;
  COFD[132] = -1.67115471E+01;
  COFD[133] = 3.98859394E+00;
  COFD[134] = -3.02316219E-01;
  COFD[135] = 1.30661099E-02;
  COFD[136] = -1.91045520E+01;
  COFD[137] = 4.87977047E+00;
  COFD[138] = -4.10448693E-01;
  COFD[139] = 1.74535827E-02;
  COFD[140] = -1.38756380E+01;
  COFD[141] = 2.98558426E+00;
  COFD[142] = -1.75507216E-01;
  COFD[143] = 7.71173691E-03;
  COFD[144] = -1.06283588E+01;
  COFD[145] = 2.15849701E+00;
  COFD[146] = -6.53886401E-02;
  COFD[147] = 2.81453370E-03;
  COFD[148] = -1.41066358E+01;
  COFD[149] = 3.07458927E+00;
  COFD[150] = -1.86899591E-01;
  COFD[151] = 8.19829781E-03;
  COFD[152] = -1.30027669E+01;
  COFD[153] = 2.80841511E+00;
  COFD[154] = -1.52629888E-01;
  COFD[155] = 6.72604927E-03;
  COFD[156] = -1.30182749E+01;
  COFD[157] = 2.80841511E+00;
  COFD[158] = -1.52629888E-01;
  COFD[159] = 6.72604927E-03;
  COFD[160] = -1.91256173E+01;
  COFD[161] = 5.02608697E+00;
  COFD[162] = -4.26959993E-01;
  COFD[163] = 1.80709910E-02;
  COFD[164] = -1.31877859E+01;
  COFD[165] = 3.38003453E+00;
  COFD[166] = -2.25783856E-01;
  COFD[167] = 9.85028660E-03;
  COFD[168] = -1.41119635E+01;
  COFD[169] = 3.07458927E+00;
  COFD[170] = -1.86899591E-01;
  COFD[171] = 8.19829781E-03;
  COFD[172] = -1.41170278E+01;
  COFD[173] = 3.07458927E+00;
  COFD[174] = -1.86899591E-01;
  COFD[175] = 8.19829781E-03;
  COFD[176] = -1.39199573E+01;
  COFD[177] = 2.99164244E+00;
  COFD[178] = -1.76293106E-01;
  COFD[179] = 7.74575100E-03;
  COFD[180] = -1.67337674E+01;
  COFD[181] = 3.98859394E+00;
  COFD[182] = -3.02316219E-01;
  COFD[183] = 1.30661099E-02;
  COFD[184] = -1.91240292E+01;
  COFD[185] = 4.87977047E+00;
  COFD[186] = -4.10448693E-01;
  COFD[187] = 1.74535827E-02;
  COFD[188] = -1.38948654E+01;
  COFD[189] = 2.98558426E+00;
  COFD[190] = -1.75507216E-01;
  COFD[191] = 7.71173691E-03;
  COFD[192] = -1.68759075E+01;
  COFD[193] = 4.49460341E+00;
  COFD[194] = -3.64766179E-01;
  COFD[195] = 1.56457173E-02;
  COFD[196] = -2.10639926E+01;
  COFD[197] = 5.50980703E+00;
  COFD[198] = -4.78335492E-01;
  COFD[199] = 1.98515434E-02;
  COFD[200] = -1.91096700E+01;
  COFD[201] = 5.02608697E+00;
  COFD[202] = -4.26959993E-01;
  COFD[203] = 1.80709910E-02;
  COFD[204] = -1.91256173E+01;
  COFD[205] = 5.02608697E+00;
  COFD[206] = -4.26959993E-01;
  COFD[207] = 1.80709910E-02;
  COFD[208] = -1.31492509E+01;
  COFD[209] = 1.48004045E+00;
  COFD[210] = 1.60499956E-01;
  COFD[211] = -1.19765871E-02;
  COFD[212] = -1.93611201E+01;
  COFD[213] = 5.51579726E+00;
  COFD[214] = -4.76061961E-01;
  COFD[215] = 1.96329391E-02;
  COFD[216] = -2.04177394E+01;
  COFD[217] = 5.31457079E+00;
  COFD[218] = -4.58216496E-01;
  COFD[219] = 1.91825910E-02;
  COFD[220] = -2.04229988E+01;
  COFD[221] = 5.31457079E+00;
  COFD[222] = -4.58216496E-01;
  COFD[223] = 1.91825910E-02;
  COFD[224] = -2.08943723E+01;
  COFD[225] = 5.44718673E+00;
  COFD[226] = -4.72082974E-01;
  COFD[227] = 1.96531328E-02;
  COFD[228] = -2.12021420E+01;
  COFD[229] = 5.20775027E+00;
  COFD[230] = -4.07348285E-01;
  COFD[231] = 1.55473262E-02;
  COFD[232] = -1.87171338E+01;
  COFD[233] = 4.00967621E+00;
  COFD[234] = -2.21153539E-01;
  COFD[235] = 6.31528745E-03;
  COFD[236] = -2.08123331E+01;
  COFD[237] = 5.42470175E+00;
  COFD[238] = -4.69700438E-01;
  COFD[239] = 1.95706912E-02;
  COFD[240] = -1.11808845E+01;
  COFD[241] = 2.66936727E+00;
  COFD[242] = -1.34411514E-01;
  COFD[243] = 5.92957488E-03;
  COFD[244] = -1.43713019E+01;
  COFD[245] = 3.70920439E+00;
  COFD[246] = -2.67274113E-01;
  COFD[247] = 1.15967481E-02;
  COFD[248] = -1.31860264E+01;
  COFD[249] = 3.38003453E+00;
  COFD[250] = -2.25783856E-01;
  COFD[251] = 9.85028660E-03;
  COFD[252] = -1.31877859E+01;
  COFD[253] = 3.38003453E+00;
  COFD[254] = -2.25783856E-01;
  COFD[255] = 9.85028660E-03;
  COFD[256] = -1.93611201E+01;
  COFD[257] = 5.51579726E+00;
  COFD[258] = -4.76061961E-01;
  COFD[259] = 1.96329391E-02;
  COFD[260] = -1.43693219E+01;
  COFD[261] = 4.03992999E+00;
  COFD[262] = -3.08044800E-01;
  COFD[263] = 1.32757775E-02;
  COFD[264] = -1.43717684E+01;
  COFD[265] = 3.70920439E+00;
  COFD[266] = -2.67274113E-01;
  COFD[267] = 1.15967481E-02;
  COFD[268] = -1.43722077E+01;
  COFD[269] = 3.70920439E+00;
  COFD[270] = -2.67274113E-01;
  COFD[271] = 1.15967481E-02;
  COFD[272] = -1.40524219E+01;
  COFD[273] = 3.56261348E+00;
  COFD[274] = -2.48287981E-01;
  COFD[275] = 1.07752947E-02;
  COFD[276] = -1.72994129E+01;
  COFD[277] = 4.71931868E+00;
  COFD[278] = -3.91258152E-01;
  COFD[279] = 1.66866639E-02;
  COFD[280] = -1.95312879E+01;
  COFD[281] = 5.47046983E+00;
  COFD[282] = -4.74577605E-01;
  COFD[283] = 1.97408822E-02;
  COFD[284] = -1.40298991E+01;
  COFD[285] = 3.55837688E+00;
  COFD[286] = -2.47785790E-01;
  COFD[287] = 1.07555332E-02;
  COFD[288] = -1.15806955E+01;
  COFD[289] = 2.43235504E+00;
  COFD[290] = -1.02890179E-01;
  COFD[291] = 4.52903603E-03;
  COFD[292] = -1.53187536E+01;
  COFD[293] = 3.37317428E+00;
  COFD[294] = -2.24900439E-01;
  COFD[295] = 9.81228151E-03;
  COFD[296] = -1.40915944E+01;
  COFD[297] = 3.07458927E+00;
  COFD[298] = -1.86899591E-01;
  COFD[299] = 8.19829781E-03;
  COFD[300] = -1.41119635E+01;
  COFD[301] = 3.07458927E+00;
  COFD[302] = -1.86899591E-01;
  COFD[303] = 8.19829781E-03;
  COFD[304] = -2.04177394E+01;
  COFD[305] = 5.31457079E+00;
  COFD[306] = -4.58216496E-01;
  COFD[307] = 1.91825910E-02;
  COFD[308] = -1.43717684E+01;
  COFD[309] = 3.70920439E+00;
  COFD[310] = -2.67274113E-01;
  COFD[311] = 1.15967481E-02;
  COFD[312] = -1.53265677E+01;
  COFD[313] = 3.37317428E+00;
  COFD[314] = -2.24900439E-01;
  COFD[315] = 9.81228151E-03;
  COFD[316] = -1.53340319E+01;
  COFD[317] = 3.37317428E+00;
  COFD[318] = -2.24900439E-01;
  COFD[319] = 9.81228151E-03;
  COFD[320] = -1.50443476E+01;
  COFD[321] = 3.26249588E+00;
  COFD[322] = -2.10658287E-01;
  COFD[323] = 9.20032462E-03;
  COFD[324] = -1.81286456E+01;
  COFD[325] = 4.33684042E+00;
  COFD[326] = -3.44981265E-01;
  COFD[327] = 1.48142449E-02;
  COFD[328] = -2.05180548E+01;
  COFD[329] = 5.21473296E+00;
  COFD[330] = -4.48646311E-01;
  COFD[331] = 1.89013813E-02;
  COFD[332] = -1.50168047E+01;
  COFD[333] = 3.25515933E+00;
  COFD[334] = -2.09710110E-01;
  COFD[335] = 9.15941830E-03;
  COFD[336] = -1.15815492E+01;
  COFD[337] = 2.43235504E+00;
  COFD[338] = -1.02890179E-01;
  COFD[339] = 4.52903603E-03;
  COFD[340] = -1.53261011E+01;
  COFD[341] = 3.37317428E+00;
  COFD[342] = -2.24900439E-01;
  COFD[343] = 9.81228151E-03;
  COFD[344] = -1.40964555E+01;
  COFD[345] = 3.07458927E+00;
  COFD[346] = -1.86899591E-01;
  COFD[347] = 8.19829781E-03;
  COFD[348] = -1.41170278E+01;
  COFD[349] = 3.07458927E+00;
  COFD[350] = -1.86899591E-01;
  COFD[351] = 8.19829781E-03;
  COFD[352] = -2.04229988E+01;
  COFD[353] = 5.31457079E+00;
  COFD[354] = -4.58216496E-01;
  COFD[355] = 1.91825910E-02;
  COFD[356] = -1.43722077E+01;
  COFD[357] = 3.70920439E+00;
  COFD[358] = -2.67274113E-01;
  COFD[359] = 1.15967481E-02;
  COFD[360] = -1.53340319E+01;
  COFD[361] = 3.37317428E+00;
  COFD[362] = -2.24900439E-01;
  COFD[363] = 9.81228151E-03;
  COFD[364] = -1.53416091E+01;
  COFD[365] = 3.37317428E+00;
  COFD[366] = -2.24900439E-01;
  COFD[367] = 9.81228151E-03;
  COFD[368] = -1.50511964E+01;
  COFD[369] = 3.26249588E+00;
  COFD[370] = -2.10658287E-01;
  COFD[371] = 9.20032462E-03;
  COFD[372] = -1.81371853E+01;
  COFD[373] = 4.33684042E+00;
  COFD[374] = -3.44981265E-01;
  COFD[375] = 1.48142449E-02;
  COFD[376] = -2.05250357E+01;
  COFD[377] = 5.21473296E+00;
  COFD[378] = -4.48646311E-01;
  COFD[379] = 1.89013813E-02;
  COFD[380] = -1.50236540E+01;
  COFD[381] = 3.25515933E+00;
  COFD[382] = -2.09710110E-01;
  COFD[383] = 9.15941830E-03;
  COFD[384] = -1.13541221E+01;
  COFD[385] = 2.31999438E+00;
  COFD[386] = -8.75064804E-02;
  COFD[387] = 3.82656365E-03;
  COFD[388] = -1.50371688E+01;
  COFD[389] = 3.26249588E+00;
  COFD[390] = -2.10658287E-01;
  COFD[391] = 9.20032462E-03;
  COFD[392] = -1.39007309E+01;
  COFD[393] = 2.99164244E+00;
  COFD[394] = -1.76293106E-01;
  COFD[395] = 7.74575100E-03;
  COFD[396] = -1.39199573E+01;
  COFD[397] = 2.99164244E+00;
  COFD[398] = -1.76293106E-01;
  COFD[399] = 7.74575100E-03;
  COFD[400] = -2.08943723E+01;
  COFD[401] = 5.44718673E+00;
  COFD[402] = -4.72082974E-01;
  COFD[403] = 1.96531328E-02;
  COFD[404] = -1.40524219E+01;
  COFD[405] = 3.56261348E+00;
  COFD[406] = -2.48287981E-01;
  COFD[407] = 1.07752947E-02;
  COFD[408] = -1.50443476E+01;
  COFD[409] = 3.26249588E+00;
  COFD[410] = -2.10658287E-01;
  COFD[411] = 9.20032462E-03;
  COFD[412] = -1.50511964E+01;
  COFD[413] = 3.26249588E+00;
  COFD[414] = -2.10658287E-01;
  COFD[415] = 9.20032462E-03;
  COFD[416] = -1.48061406E+01;
  COFD[417] = 3.16912473E+00;
  COFD[418] = -1.98792456E-01;
  COFD[419] = 8.69726395E-03;
  COFD[420] = -1.77672912E+01;
  COFD[421] = 4.20234040E+00;
  COFD[422] = -3.28057658E-01;
  COFD[423] = 1.41006192E-02;
  COFD[424] = -2.02361782E+01;
  COFD[425] = 5.11785645E+00;
  COFD[426] = -4.37867828E-01;
  COFD[427] = 1.85047543E-02;
  COFD[428] = -1.47850505E+01;
  COFD[429] = 3.16433919E+00;
  COFD[430] = -1.98191564E-01;
  COFD[431] = 8.67209742E-03;
  COFD[432] = -1.35545391E+01;
  COFD[433] = 3.13878730E+00;
  COFD[434] = -1.94980335E-01;
  COFD[435] = 8.53744486E-03;
  COFD[436] = -1.81197250E+01;
  COFD[437] = 4.33684042E+00;
  COFD[438] = -3.44981265E-01;
  COFD[439] = 1.48142449E-02;
  COFD[440] = -1.67115471E+01;
  COFD[441] = 3.98859394E+00;
  COFD[442] = -3.02316219E-01;
  COFD[443] = 1.30661099E-02;
  COFD[444] = -1.67337674E+01;
  COFD[445] = 3.98859394E+00;
  COFD[446] = -3.02316219E-01;
  COFD[447] = 1.30661099E-02;
  COFD[448] = -2.12021420E+01;
  COFD[449] = 5.20775027E+00;
  COFD[450] = -4.07348285E-01;
  COFD[451] = 1.55473262E-02;
  COFD[452] = -1.72994129E+01;
  COFD[453] = 4.71931868E+00;
  COFD[454] = -3.91258152E-01;
  COFD[455] = 1.66866639E-02;
  COFD[456] = -1.81286456E+01;
  COFD[457] = 4.33684042E+00;
  COFD[458] = -3.44981265E-01;
  COFD[459] = 1.48142449E-02;
  COFD[460] = -1.81371853E+01;
  COFD[461] = 4.33684042E+00;
  COFD[462] = -3.44981265E-01;
  COFD[463] = 1.48142449E-02;
  COFD[464] = -1.77672912E+01;
  COFD[465] = 4.20234040E+00;
  COFD[466] = -3.28057658E-01;
  COFD[467] = 1.41006192E-02;
  COFD[468] = -2.10907633E+01;
  COFD[469] = 5.29211327E+00;
  COFD[470] = -4.56068366E-01;
  COFD[471] = 1.91195062E-02;
  COFD[472] = -2.20758670E+01;
  COFD[473] = 5.52171573E+00;
  COFD[474] = -4.63284984E-01;
  COFD[475] = 1.85570924E-02;
  COFD[476] = -1.77350630E+01;
  COFD[477] = 4.19328271E+00;
  COFD[478] = -3.26911461E-01;
  COFD[479] = 1.40520357E-02;
  COFD[480] = -1.57161352E+01;
  COFD[481] = 3.96062263E+00;
  COFD[482] = -2.98964970E-01;
  COFD[483] = 1.29322565E-02;
  COFD[484] = -2.05107395E+01;
  COFD[485] = 5.21473296E+00;
  COFD[486] = -4.48646311E-01;
  COFD[487] = 1.89013813E-02;
  COFD[488] = -1.91045520E+01;
  COFD[489] = 4.87977047E+00;
  COFD[490] = -4.10448693E-01;
  COFD[491] = 1.74535827E-02;
  COFD[492] = -1.91240292E+01;
  COFD[493] = 4.87977047E+00;
  COFD[494] = -4.10448693E-01;
  COFD[495] = 1.74535827E-02;
  COFD[496] = -1.87171338E+01;
  COFD[497] = 4.00967621E+00;
  COFD[498] = -2.21153539E-01;
  COFD[499] = 6.31528745E-03;
  COFD[500] = -1.95312879E+01;
  COFD[501] = 5.47046983E+00;
  COFD[502] = -4.74577605E-01;
  COFD[503] = 1.97408822E-02;
  COFD[504] = -2.05180548E+01;
  COFD[505] = 5.21473296E+00;
  COFD[506] = -4.48646311E-01;
  COFD[507] = 1.89013813E-02;
  COFD[508] = -2.05250357E+01;
  COFD[509] = 5.21473296E+00;
  COFD[510] = -4.48646311E-01;
  COFD[511] = 1.89013813E-02;
  COFD[512] = -2.02361782E+01;
  COFD[513] = 5.11785645E+00;
  COFD[514] = -4.37867828E-01;
  COFD[515] = 1.85047543E-02;
  COFD[516] = -2.20758670E+01;
  COFD[517] = 5.52171573E+00;
  COFD[518] = -4.63284984E-01;
  COFD[519] = 1.85570924E-02;
  COFD[520] = -1.98983685E+01;
  COFD[521] = 4.38041133E+00;
  COFD[522] = -2.77538214E-01;
  COFD[523] = 9.06748822E-03;
  COFD[524] = -2.02052920E+01;
  COFD[525] = 5.10993120E+00;
  COFD[526] = -4.36931630E-01;
  COFD[527] = 1.84677592E-02;
  COFD[528] = -1.13253618E+01;
  COFD[529] = 2.31195095E+00;
  COFD[530] = -8.63988037E-02;
  COFD[531] = 3.77573452E-03;
  COFD[532] = -1.50096252E+01;
  COFD[533] = 3.25515933E+00;
  COFD[534] = -2.09710110E-01;
  COFD[535] = 9.15941830E-03;
  COFD[536] = -1.38756380E+01;
  COFD[537] = 2.98558426E+00;
  COFD[538] = -1.75507216E-01;
  COFD[539] = 7.71173691E-03;
  COFD[540] = -1.38948654E+01;
  COFD[541] = 2.98558426E+00;
  COFD[542] = -1.75507216E-01;
  COFD[543] = 7.71173691E-03;
  COFD[544] = -2.08123331E+01;
  COFD[545] = 5.42470175E+00;
  COFD[546] = -4.69700438E-01;
  COFD[547] = 1.95706912E-02;
  COFD[548] = -1.40298991E+01;
  COFD[549] = 3.55837688E+00;
  COFD[550] = -2.47785790E-01;
  COFD[551] = 1.07555332E-02;
  COFD[552] = -1.50168047E+01;
  COFD[553] = 3.25515933E+00;
  COFD[554] = -2.09710110E-01;
  COFD[555] = 9.15941830E-03;
  COFD[556] = -1.50236540E+01;
  COFD[557] = 3.25515933E+00;
  COFD[558] = -2.09710110E-01;
  COFD[559] = 9.15941830E-03;
  COFD[560] = -1.47850505E+01;
  COFD[561] = 3.16433919E+00;
  COFD[562] = -1.98191564E-01;
  COFD[563] = 8.67209742E-03;
  COFD[564] = -1.77350630E+01;
  COFD[565] = 4.19328271E+00;
  COFD[566] = -3.26911461E-01;
  COFD[567] = 1.40520357E-02;
  COFD[568] = -2.02052920E+01;
  COFD[569] = 5.10993120E+00;
  COFD[570] = -4.36931630E-01;
  COFD[571] = 1.84677592E-02;
  COFD[572] = -1.47639411E+01;
  COFD[573] = 3.15955654E+00;
  COFD[574] = -1.97590757E-01;
  COFD[575] = 8.64692156E-03;
}

// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKTDIF(int* KTDIF)
{
  KTDIF[0] = 0;
  KTDIF[1] = 5;
}

// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFTD(amrex::Real* COFTD)
{
  COFTD[0] = 0.00000000E+00;
  COFTD[1] = 0.00000000E+00;
  COFTD[2] = 0.00000000E+00;
  COFTD[3] = 0.00000000E+00;
  COFTD[4] = 4.42736017E-01;
  COFTD[5] = 7.11765887E-05;
  COFTD[6] = -3.84765396E-08;
  COFTD[7] = 6.86318682E-12;
  COFTD[8] = 4.15577509E-01;
  COFTD[9] = 1.09736860E-05;
  COFTD[10] = -3.96016410E-09;
  COFTD[11] = 1.14412838E-12;
  COFTD[12] = 4.21927217E-01;
  COFTD[13] = 1.11413555E-05;
  COFTD[14] = -4.02067239E-09;
  COFTD[15] = 1.16160979E-12;
  COFTD[16] = 6.02021050E-02;
  COFTD[17] = 5.61555733E-04;
  COFTD[18] = -2.55370065E-07;
  COFTD[19] = 3.63385926E-11;
  COFTD[20] = -1.52534742E-01;
  COFTD[21] = -5.46404022E-05;
  COFTD[22] = 2.93412470E-08;
  COFTD[23] = -4.87091914E-12;
  COFTD[24] = 4.44449675E-01;
  COFTD[25] = 7.14520855E-05;
  COFTD[26] = -3.86254673E-08;
  COFTD[27] = 6.88975155E-12;
  COFTD[28] = 4.46067449E-01;
  COFTD[29] = 7.17121674E-05;
  COFTD[30] = -3.87660620E-08;
  COFTD[31] = 6.91482989E-12;
  COFTD[32] = 4.44650439E-01;
  COFTD[33] = 5.06628083E-05;
  COFTD[34] = -2.69818971E-08;
  COFTD[35] = 5.01286176E-12;
  COFTD[36] = 3.25740914E-01;
  COFTD[37] = 3.03631979E-04;
  COFTD[38] = -1.55289598E-07;
  COFTD[39] = 2.41465297E-11;
  COFTD[40] = 1.61391212E-01;
  COFTD[41] = 5.01080793E-04;
  COFTD[42] = -2.38272376E-07;
  COFTD[43] = 3.49342098E-11;
  COFTD[44] = 4.45261429E-01;
  COFTD[45] = 4.94696577E-05;
  COFTD[46] = -2.63023124E-08;
  COFTD[47] = 4.90305626E-12;
  COFTD[48] = 1.52534742E-01;
  COFTD[49] = 5.46404022E-05;
  COFTD[50] = -2.93412470E-08;
  COFTD[51] = 4.87091914E-12;
  COFTD[52] = 2.20482081E-01;
  COFTD[53] = 4.80162630E-04;
  COFTD[54] = -2.32927140E-07;
  COFTD[55] = 3.46469239E-11;
  COFTD[56] = 2.70008280E-01;
  COFTD[57] = 3.61552588E-04;
  COFTD[58] = -1.80743499E-07;
  COFTD[59] = 2.75319341E-11;
  COFTD[60] = 2.72039997E-01;
  COFTD[61] = 3.64273144E-04;
  COFTD[62] = -1.82103531E-07;
  COFTD[63] = 2.77391022E-11;
  COFTD[64] = -1.41882969E-01;
  COFTD[65] = 7.66554625E-04;
  COFTD[66] = -3.06548330E-07;
  COFTD[67] = 4.02957302E-11;
  COFTD[68] = 0.00000000E+00;
  COFTD[69] = 0.00000000E+00;
  COFTD[70] = 0.00000000E+00;
  COFTD[71] = 0.00000000E+00;
  COFTD[72] = 2.20907136E-01;
  COFTD[73] = 4.81088308E-04;
  COFTD[74] = -2.33376187E-07;
  COFTD[75] = 3.47137178E-11;
  COFTD[76] = 2.21307723E-01;
  COFTD[77] = 4.81960701E-04;
  COFTD[78] = -2.33799385E-07;
  COFTD[79] = 3.47766667E-11;
  COFTD[80] = 2.39409087E-01;
  COFTD[81] = 4.47195587E-04;
  COFTD[82] = -2.18950923E-07;
  COFTD[83] = 3.27972343E-11;
  COFTD[84] = 2.44368810E-02;
  COFTD[85] = 7.18240808E-04;
  COFTD[86] = -3.19717751E-07;
  COFTD[87] = 4.48827629E-11;
  COFTD[88] = -1.24647578E-01;
  COFTD[89] = 7.96522972E-04;
  COFTD[90] = -3.24997704E-07;
  COFTD[91] = 4.32515301E-11;
  COFTD[92] = 2.40744276E-01;
  COFTD[93] = 4.45343183E-04;
  COFTD[94] = -2.18173742E-07;
  COFTD[95] = 3.26958310E-11;
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

  // species 0: H2
  // Imported from NIST
  Tci[0] = 33.145000;
  ai[0] = 1e6 * 0.42748 * Rcst * Rcst * Tci[0] * Tci[0] /
          (2.015880 * 2.015880 * 12.964000);
  bi[0] = 0.08664 * Rcst * Tci[0] / (2.015880 * 12.964000);
  acentric_i[0] = -0.219000;

  // species 1: O2
  // Imported from NIST
  Tci[1] = 154.581000;
  ai[1] = 1e6 * 0.42748 * Rcst * Rcst * Tci[1] * Tci[1] /
          (31.998800 * 31.998800 * 50.430466);
  bi[1] = 0.08664 * Rcst * Tci[1] / (31.998800 * 50.430466);
  acentric_i[1] = 0.022200;

  // species 2: O
  Tci[2] = 1.316 * EPS[2];
  ai[2] = (5.55 * avogadro * avogadro * EPS[2] * boltzmann * 1e-24 * SIG[2] *
           SIG[2] * SIG[2]) /
          (wt[2] * wt[2]);
  bi[2] = 0.855 * avogadro * 1e-24 * SIG[2] * SIG[2] * SIG[2] / (wt[2]);
  acentric_i[2] = 0.0;

  // species 3: OH
  Tci[3] = 1.316 * EPS[3];
  ai[3] = (5.55 * avogadro * avogadro * EPS[3] * boltzmann * 1e-24 * SIG[3] *
           SIG[3] * SIG[3]) /
          (wt[3] * wt[3]);
  bi[3] = 0.855 * avogadro * 1e-24 * SIG[3] * SIG[3] * SIG[3] / (wt[3]);
  acentric_i[3] = 0.0;

  // species 4: H2O
  // Imported from NIST
  Tci[4] = 647.096000;
  ai[4] = 1e6 * 0.42748 * Rcst * Rcst * Tci[4] * Tci[4] /
          (18.015340 * 18.015340 * 220.640000);
  bi[4] = 0.08664 * Rcst * Tci[4] / (18.015340 * 220.640000);
  acentric_i[4] = 0.344300;

  // species 5: H
  Tci[5] = 1.316 * EPS[5];
  ai[5] = (5.55 * avogadro * avogadro * EPS[5] * boltzmann * 1e-24 * SIG[5] *
           SIG[5] * SIG[5]) /
          (wt[5] * wt[5]);
  bi[5] = 0.855 * avogadro * 1e-24 * SIG[5] * SIG[5] * SIG[5] / (wt[5]);
  acentric_i[5] = 0.0;

  // species 6: HO2
  Tci[6] = 1.316 * EPS[6];
  ai[6] = (5.55 * avogadro * avogadro * EPS[6] * boltzmann * 1e-24 * SIG[6] *
           SIG[6] * SIG[6]) /
          (wt[6] * wt[6]);
  bi[6] = 0.855 * avogadro * 1e-24 * SIG[6] * SIG[6] * SIG[6] / (wt[6]);
  acentric_i[6] = 0.0;

  // species 7: H2O2
  Tci[7] = 1.316 * EPS[7];
  ai[7] = (5.55 * avogadro * avogadro * EPS[7] * boltzmann * 1e-24 * SIG[7] *
           SIG[7] * SIG[7]) /
          (wt[7] * wt[7]);
  bi[7] = 0.855 * avogadro * 1e-24 * SIG[7] * SIG[7] * SIG[7] / (wt[7]);
  acentric_i[7] = 0.0;

  // species 8: CO
  // Imported from NIST
  Tci[8] = 132.850000;
  ai[8] = 1e6 * 0.42748 * Rcst * Rcst * Tci[8] * Tci[8] /
          (28.010000 * 28.010000 * 34.940000);
  bi[8] = 0.08664 * Rcst * Tci[8] / (28.010000 * 34.940000);
  acentric_i[8] = 0.045000;

  // species 9: CO2
  // Imported from NIST
  Tci[9] = 304.120000;
  ai[9] = 1e6 * 0.42748 * Rcst * Rcst * Tci[9] * Tci[9] /
          (44.009950 * 44.009950 * 73.740000);
  bi[9] = 0.08664 * Rcst * Tci[9] / (44.009950 * 73.740000);
  acentric_i[9] = 0.225000;

  // species 10: HCO
  Tci[10] = 1.316 * EPS[10];
  ai[10] = (5.55 * avogadro * avogadro * EPS[10] * boltzmann * 1e-24 * SIG[10] *
            SIG[10] * SIG[10]) /
           (wt[10] * wt[10]);
  bi[10] = 0.855 * avogadro * 1e-24 * SIG[10] * SIG[10] * SIG[10] / (wt[10]);
  acentric_i[10] = 0.0;

  // species 11: N2
  // Imported from NIST
  Tci[11] = 126.192000;
  ai[11] = 1e6 * 0.42748 * Rcst * Rcst * Tci[11] * Tci[11] /
           (28.013400 * 28.013400 * 33.958000);
  bi[11] = 0.08664 * Rcst * Tci[11] / (28.013400 * 33.958000);
  acentric_i[11] = 0.037200;
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
