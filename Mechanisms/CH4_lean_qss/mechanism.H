#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  H
1  O
2  C
3  N
*/

// Species
#define H2_ID 0
#define H_ID 1
#define O_ID 2
#define O2_ID 3
#define OH_ID 4
#define H2O_ID 5
#define HO2_ID 6
#define CH3_ID 7
#define CH4_ID 8
#define CO_ID 9
#define CO2_ID 10
#define CH2O_ID 11
#define N2_ID 12

#define NUM_ELEMENTS 4
#define NUM_SPECIES 13
#define NUM_IONS 0
#define NUM_REACTIONS 73

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
  mm = 4;
  kk = 13;
  ii = 73;
  nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[13] = {
  0.4960317460317460, // H2
  0.9920634920634921, // H
  0.0625039064941559, // O
  0.0312519532470779, // O2
  0.0587993179279120, // OH
  0.0555092978073827, // H2O
  0.0302975216627280, // HO2
  0.0665114732291320, // CH3
  0.0623324814560868, // CH4
  0.0357015351660121, // CO
  0.0227226249176305, // CO2
  0.0333044694598015, // CH2O
  0.0356964374955379, // N2
};
#endif
const amrex::Real h_global_imw[13] = {
  0.4960317460317460, // H2
  0.9920634920634921, // H
  0.0625039064941559, // O
  0.0312519532470779, // O2
  0.0587993179279120, // OH
  0.0555092978073827, // H2O
  0.0302975216627280, // HO2
  0.0665114732291320, // CH3
  0.0623324814560868, // CH4
  0.0357015351660121, // CO
  0.0227226249176305, // CO2
  0.0333044694598015, // CH2O
  0.0356964374955379, // N2
};

//  molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[13] = {
  2.016000,  // H2
  1.008000,  // H
  15.999000, // O
  31.998000, // O2
  17.007000, // OH
  18.015000, // H2O
  33.006000, // HO2
  15.035000, // CH3
  16.043000, // CH4
  28.010000, // CO
  44.009000, // CO2
  30.026000, // CH2O
  28.014000, // N2
};
#endif
const amrex::Real h_global_mw[13] = {
  2.016000,  // H2
  1.008000,  // H
  15.999000, // O
  31.998000, // O2
  17.007000, // OH
  18.015000, // H2O
  33.006000, // HO2
  15.035000, // CH3
  16.043000, // CH4
  28.010000, // CO
  44.009000, // CO2
  30.026000, // CH2O
  28.014000, // N2
};

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_imw(amrex::Real* imw_new)
{
  imw_new[0] = 0.4960317460317460;  // H2
  imw_new[1] = 0.9920634920634921;  // H
  imw_new[2] = 0.0625039064941559;  // O
  imw_new[3] = 0.0312519532470779;  // O2
  imw_new[4] = 0.0587993179279120;  // OH
  imw_new[5] = 0.0555092978073827;  // H2O
  imw_new[6] = 0.0302975216627280;  // HO2
  imw_new[7] = 0.0665114732291320;  // CH3
  imw_new[8] = 0.0623324814560868;  // CH4
  imw_new[9] = 0.0357015351660121;  // CO
  imw_new[10] = 0.0227226249176305; // CO2
  imw_new[11] = 0.0333044694598015; // CH2O
  imw_new[12] = 0.0356964374955379; // N2
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
  mw_new[0] = 2.016000;   // H2
  mw_new[1] = 1.008000;   // H
  mw_new[2] = 15.999000;  // O
  mw_new[3] = 31.998000;  // O2
  mw_new[4] = 17.007000;  // OH
  mw_new[5] = 18.015000;  // H2O
  mw_new[6] = 33.006000;  // HO2
  mw_new[7] = 15.035000;  // CH3
  mw_new[8] = 16.043000;  // CH4
  mw_new[9] = 28.010000;  // CO
  mw_new[10] = 44.009000; // CO2
  mw_new[11] = 30.026000; // CH2O
  mw_new[12] = 28.014000; // N2
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

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: H2
    species[0] = +1.34433112e+00 + 7.98052075e-03 * T - 1.94781510e-05 * T2 +
                 2.01572094e-08 * T3 - 7.37611761e-12 * T4;
    // species 1: H
    species[1] = +1.50000000e+00 + 7.05332819e-13 * T - 1.99591964e-15 * T2 +
                 2.30081632e-18 * T3 - 9.27732332e-22 * T4;
    // species 2: O
    species[2] = +2.16826710e+00 - 3.27931884e-03 * T + 6.64306396e-06 * T2 -
                 6.12806624e-09 * T3 + 2.11265971e-12 * T4;
    // species 3: O2
    species[3] = +2.78245636e+00 - 2.99673416e-03 * T + 9.84730201e-06 * T2 -
                 9.68129509e-09 * T3 + 3.24372837e-12 * T4;
    // species 4: OH
    species[4] = +2.99201543e+00 - 2.40131752e-03 * T + 4.61793841e-06 * T2 -
                 3.88113333e-09 * T3 + 1.36411470e-12 * T4;
    // species 5: H2O
    species[5] = +3.19864056e+00 - 2.03643410e-03 * T + 6.52040211e-06 * T2 -
                 5.48797062e-09 * T3 + 1.77197817e-12 * T4;
    // species 6: HO2
    species[6] = +3.30179801e+00 - 4.74912051e-03 * T + 2.11582891e-05 * T2 -
                 2.42763894e-08 * T3 + 9.29225124e-12 * T4;
    // species 7: CH3
    species[7] = +2.67359040e+00 + 2.01095175e-03 * T + 5.73021856e-06 * T2 -
                 6.87117425e-09 * T3 + 2.54385734e-12 * T4;
    // species 8: CH4
    species[8] = +4.14987613e+00 - 1.36709788e-02 * T + 4.91800599e-05 * T2 -
                 4.84743026e-08 * T3 + 1.66693956e-11 * T4;
    // species 9: CO
    species[9] = +2.57953347e+00 - 6.10353680e-04 * T + 1.01681433e-06 * T2 +
                 9.07005884e-10 * T3 - 9.04424499e-13 * T4;
    // species 10: CO2
    species[10] = +1.35677352e+00 + 8.98459677e-03 * T - 7.12356269e-06 * T2 +
                  2.45919022e-09 * T3 - 1.43699548e-13 * T4;
    // species 11: CH2O
    species[11] = +3.79372315e+00 - 9.90833369e-03 * T + 3.73220008e-05 * T2 -
                  3.79285261e-08 * T3 + 1.31772652e-11 * T4;
    // species 12: N2
    species[12] = +2.29867700e+00 + 1.40824040e-03 * T - 3.96322200e-06 * T2 +
                  5.64151500e-09 * T3 - 2.44485400e-12 * T4;
  } else {
    // species 0: H2
    species[0] = +2.33727920e+00 - 4.94024731e-05 * T + 4.99456778e-07 * T2 -
                 1.79566394e-10 * T3 + 2.00255376e-14 * T4;
    // species 1: H
    species[1] = +1.50000001e+00 - 2.30842973e-11 * T + 1.61561948e-14 * T2 -
                 4.73515235e-18 * T3 + 4.98197357e-22 * T4;
    // species 2: O
    species[2] = +1.56942078e+00 - 8.59741137e-05 * T + 4.19484589e-08 * T2 -
                 1.00177799e-11 * T3 + 1.22833691e-15 * T4;
    // species 3: O2
    species[3] = +2.28253784e+00 + 1.48308754e-03 * T - 7.57966669e-07 * T2 +
                 2.09470555e-10 * T3 - 2.16717794e-14 * T4;
    // species 4: OH
    species[4] = +2.09288767e+00 + 5.48429716e-04 * T + 1.26505228e-07 * T2 -
                 8.79461556e-11 * T3 + 1.17412376e-14 * T4;
    // species 5: H2O
    species[5] = +2.03399249e+00 + 2.17691804e-03 * T - 1.64072518e-07 * T2 -
                 9.70419870e-11 * T3 + 1.68200992e-14 * T4;
    // species 6: HO2
    species[6] = +3.01721090e+00 + 2.23982013e-03 * T - 6.33658150e-07 * T2 +
                 1.14246370e-10 * T3 - 1.07908535e-14 * T4;
    // species 7: CH3
    species[7] = +1.28571772e+00 + 7.23990037e-03 * T - 2.98714348e-06 * T2 +
                 5.95684644e-10 * T3 - 4.67154394e-14 * T4;
    // species 8: CH4
    species[8] = -9.25148505e-01 + 1.33909467e-02 * T - 5.73285809e-06 * T2 +
                 1.22292535e-09 * T3 - 1.01815230e-13 * T4;
    // species 9: CO
    species[9] = +1.71518561e+00 + 2.06252743e-03 * T - 9.98825771e-07 * T2 +
                 2.30053008e-10 * T3 - 2.03647716e-14 * T4;
    // species 10: CO2
    species[10] = +2.85746029e+00 + 4.41437026e-03 * T - 2.21481404e-06 * T2 +
                  5.23490188e-10 * T3 - 4.72084164e-14 * T4;
    // species 11: CH2O
    species[11] = +7.60690080e-01 + 9.20000082e-03 * T - 4.42258813e-06 * T2 +
                  1.00641212e-09 * T3 - 8.83855640e-14 * T4;
    // species 12: N2
    species[12] = +1.92664000e+00 + 1.48797680e-03 * T - 5.68476000e-07 * T2 +
                  1.00970380e-10 * T3 - 6.75335100e-15 * T4;
  }
}

// compute Cp/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cp_R(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: H2
    species[0] = +2.34433112e+00 + 7.98052075e-03 * T - 1.94781510e-05 * T2 +
                 2.01572094e-08 * T3 - 7.37611761e-12 * T4;
    // species 1: H
    species[1] = +2.50000000e+00 + 7.05332819e-13 * T - 1.99591964e-15 * T2 +
                 2.30081632e-18 * T3 - 9.27732332e-22 * T4;
    // species 2: O
    species[2] = +3.16826710e+00 - 3.27931884e-03 * T + 6.64306396e-06 * T2 -
                 6.12806624e-09 * T3 + 2.11265971e-12 * T4;
    // species 3: O2
    species[3] = +3.78245636e+00 - 2.99673416e-03 * T + 9.84730201e-06 * T2 -
                 9.68129509e-09 * T3 + 3.24372837e-12 * T4;
    // species 4: OH
    species[4] = +3.99201543e+00 - 2.40131752e-03 * T + 4.61793841e-06 * T2 -
                 3.88113333e-09 * T3 + 1.36411470e-12 * T4;
    // species 5: H2O
    species[5] = +4.19864056e+00 - 2.03643410e-03 * T + 6.52040211e-06 * T2 -
                 5.48797062e-09 * T3 + 1.77197817e-12 * T4;
    // species 6: HO2
    species[6] = +4.30179801e+00 - 4.74912051e-03 * T + 2.11582891e-05 * T2 -
                 2.42763894e-08 * T3 + 9.29225124e-12 * T4;
    // species 7: CH3
    species[7] = +3.67359040e+00 + 2.01095175e-03 * T + 5.73021856e-06 * T2 -
                 6.87117425e-09 * T3 + 2.54385734e-12 * T4;
    // species 8: CH4
    species[8] = +5.14987613e+00 - 1.36709788e-02 * T + 4.91800599e-05 * T2 -
                 4.84743026e-08 * T3 + 1.66693956e-11 * T4;
    // species 9: CO
    species[9] = +3.57953347e+00 - 6.10353680e-04 * T + 1.01681433e-06 * T2 +
                 9.07005884e-10 * T3 - 9.04424499e-13 * T4;
    // species 10: CO2
    species[10] = +2.35677352e+00 + 8.98459677e-03 * T - 7.12356269e-06 * T2 +
                  2.45919022e-09 * T3 - 1.43699548e-13 * T4;
    // species 11: CH2O
    species[11] = +4.79372315e+00 - 9.90833369e-03 * T + 3.73220008e-05 * T2 -
                  3.79285261e-08 * T3 + 1.31772652e-11 * T4;
    // species 12: N2
    species[12] = +3.29867700e+00 + 1.40824040e-03 * T - 3.96322200e-06 * T2 +
                  5.64151500e-09 * T3 - 2.44485400e-12 * T4;
  } else {
    // species 0: H2
    species[0] = +3.33727920e+00 - 4.94024731e-05 * T + 4.99456778e-07 * T2 -
                 1.79566394e-10 * T3 + 2.00255376e-14 * T4;
    // species 1: H
    species[1] = +2.50000001e+00 - 2.30842973e-11 * T + 1.61561948e-14 * T2 -
                 4.73515235e-18 * T3 + 4.98197357e-22 * T4;
    // species 2: O
    species[2] = +2.56942078e+00 - 8.59741137e-05 * T + 4.19484589e-08 * T2 -
                 1.00177799e-11 * T3 + 1.22833691e-15 * T4;
    // species 3: O2
    species[3] = +3.28253784e+00 + 1.48308754e-03 * T - 7.57966669e-07 * T2 +
                 2.09470555e-10 * T3 - 2.16717794e-14 * T4;
    // species 4: OH
    species[4] = +3.09288767e+00 + 5.48429716e-04 * T + 1.26505228e-07 * T2 -
                 8.79461556e-11 * T3 + 1.17412376e-14 * T4;
    // species 5: H2O
    species[5] = +3.03399249e+00 + 2.17691804e-03 * T - 1.64072518e-07 * T2 -
                 9.70419870e-11 * T3 + 1.68200992e-14 * T4;
    // species 6: HO2
    species[6] = +4.01721090e+00 + 2.23982013e-03 * T - 6.33658150e-07 * T2 +
                 1.14246370e-10 * T3 - 1.07908535e-14 * T4;
    // species 7: CH3
    species[7] = +2.28571772e+00 + 7.23990037e-03 * T - 2.98714348e-06 * T2 +
                 5.95684644e-10 * T3 - 4.67154394e-14 * T4;
    // species 8: CH4
    species[8] = +7.48514950e-02 + 1.33909467e-02 * T - 5.73285809e-06 * T2 +
                 1.22292535e-09 * T3 - 1.01815230e-13 * T4;
    // species 9: CO
    species[9] = +2.71518561e+00 + 2.06252743e-03 * T - 9.98825771e-07 * T2 +
                 2.30053008e-10 * T3 - 2.03647716e-14 * T4;
    // species 10: CO2
    species[10] = +3.85746029e+00 + 4.41437026e-03 * T - 2.21481404e-06 * T2 +
                  5.23490188e-10 * T3 - 4.72084164e-14 * T4;
    // species 11: CH2O
    species[11] = +1.76069008e+00 + 9.20000082e-03 * T - 4.42258813e-06 * T2 +
                  1.00641212e-09 * T3 - 8.83855640e-14 * T4;
    // species 12: N2
    species[12] = +2.92664000e+00 + 1.48797680e-03 * T - 5.68476000e-07 * T2 +
                  1.00970380e-10 * T3 - 6.75335100e-15 * T4;
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
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: H2
    species[0] = -9.179351730000000e+02 * invT + 1.661320882000000e+00 -
                 2.344331120000000e+00 * logT - 3.990260375000000e-03 * T +
                 3.246358500000000e-06 * T2 - 1.679767450000000e-09 * T3 +
                 3.688058805000000e-13 * T4;
    // species 1: H
    species[1] = +2.547365990000000e+04 * invT + 2.946682853000000e+00 -
                 2.500000000000000e+00 * logT - 3.526664095000000e-13 * T +
                 3.326532733333333e-16 * T2 - 1.917346933333333e-19 * T3 +
                 4.638661660000000e-23 * T4;
    // species 2: O
    species[2] = +2.912225920000000e+04 * invT + 1.116333640000000e+00 -
                 3.168267100000000e+00 * logT + 1.639659420000000e-03 * T -
                 1.107177326666667e-06 * T2 + 5.106721866666666e-10 * T3 -
                 1.056329855000000e-13 * T4;
    // species 3: O2
    species[3] = -1.063943560000000e+03 * invT + 1.247806300000001e-01 -
                 3.782456360000000e+00 * logT + 1.498367080000000e-03 * T -
                 1.641217001666667e-06 * T2 + 8.067745908333334e-10 * T3 -
                 1.621864185000000e-13 * T4;
    // species 4: OH
    species[4] = +3.615080560000000e+03 * invT + 4.095940888000000e+00 -
                 3.992015430000000e+00 * logT + 1.200658760000000e-03 * T -
                 7.696564016666666e-07 * T2 + 3.234277775000000e-10 * T3 -
                 6.820573500000000e-14 * T4;
    // species 5: H2O
    species[5] = -3.029372670000000e+04 * invT + 5.047672768000000e+00 -
                 4.198640560000000e+00 * logT + 1.018217050000000e-03 * T -
                 1.086733685000000e-06 * T2 + 4.573308850000000e-10 * T3 -
                 8.859890850000000e-14 * T4;
    // species 6: HO2
    species[6] = +2.948080400000000e+02 * invT + 5.851355599999999e-01 -
                 4.301798010000000e+00 * logT + 2.374560255000000e-03 * T -
                 3.526381516666666e-06 * T2 + 2.023032450000000e-09 * T3 -
                 4.646125620000001e-13 * T4;
    // species 7: CH3
    species[7] = +1.644499880000000e+04 * invT + 2.069026070000000e+00 -
                 3.673590400000000e+00 * logT - 1.005475875000000e-03 * T -
                 9.550364266666668e-07 * T2 + 5.725978541666666e-10 * T3 -
                 1.271928670000000e-13 * T4;
    // species 8: CH4
    species[8] = -1.024664760000000e+04 * invT + 9.791179889999999e+00 -
                 5.149876130000000e+00 * logT + 6.835489400000000e-03 * T -
                 8.196676650000000e-06 * T2 + 4.039525216666667e-09 * T3 -
                 8.334697800000000e-13 * T4;
    // species 9: CO
    species[9] = -1.434408600000000e+04 * invT + 7.112418999999992e-02 -
                 3.579533470000000e+00 * logT + 3.051768400000000e-04 * T -
                 1.694690550000000e-07 * T2 - 7.558382366666667e-11 * T3 +
                 4.522122495000000e-14 * T4;
    // species 10: CO2
    species[10] = -4.837196970000000e+04 * invT - 7.544278700000000e+00 -
                  2.356773520000000e+00 * logT - 4.492298385000000e-03 * T +
                  1.187260448333333e-06 * T2 - 2.049325183333333e-10 * T3 +
                  7.184977399999999e-15 * T4;
    // species 11: CH2O
    species[11] = -1.430895670000000e+04 * invT + 4.190910250000000e+00 -
                  4.793723150000000e+00 * logT + 4.954166845000000e-03 * T -
                  6.220333466666666e-06 * T2 + 3.160710508333333e-09 * T3 -
                  6.588632600000000e-13 * T4;
    // species 12: N2
    species[12] = -1.020899900000000e+03 * invT - 6.516950000000001e-01 -
                  3.298677000000000e+00 * logT - 7.041202000000000e-04 * T +
                  6.605369999999999e-07 * T2 - 4.701262500000001e-10 * T3 +
                  1.222427000000000e-13 * T4;
  } else {
    // species 0: H2
    species[0] = -9.501589220000000e+02 * invT + 6.542302510000000e+00 -
                 3.337279200000000e+00 * logT + 2.470123655000000e-05 * T -
                 8.324279633333333e-08 * T2 + 1.496386616666667e-11 * T3 -
                 1.001276880000000e-15 * T4;
    // species 1: H
    species[1] = +2.547365990000000e+04 * invT + 2.946682924000000e+00 -
                 2.500000010000000e+00 * logT + 1.154214865000000e-11 * T -
                 2.692699133333334e-15 * T2 + 3.945960291666667e-19 * T3 -
                 2.490986785000000e-23 * T4;
    // species 2: O
    species[2] = +2.921757910000000e+04 * invT - 2.214917859999999e+00 -
                 2.569420780000000e+00 * logT + 4.298705685000000e-05 * T -
                 6.991409816666667e-09 * T2 + 8.348149916666666e-13 * T3 -
                 6.141684549999999e-17 * T4;
    // species 3: O2
    species[3] = -1.088457720000000e+03 * invT - 2.170693450000000e+00 -
                 3.282537840000000e+00 * logT - 7.415437700000000e-04 * T +
                 1.263277781666667e-07 * T2 - 1.745587958333333e-11 * T3 +
                 1.083588970000000e-15 * T4;
    // species 4: OH
    species[4] = +3.858657000000000e+03 * invT - 1.383808430000000e+00 -
                 3.092887670000000e+00 * logT - 2.742148580000000e-04 * T -
                 2.108420466666667e-08 * T2 + 7.328846300000000e-12 * T3 -
                 5.870618800000000e-16 * T4;
    // species 5: H2O
    species[5] = -3.000429710000000e+04 * invT - 1.932777610000000e+00 -
                 3.033992490000000e+00 * logT - 1.088459020000000e-03 * T +
                 2.734541966666666e-08 * T2 + 8.086832250000000e-12 * T3 -
                 8.410049600000000e-16 * T4;
    // species 6: HO2
    species[6] = +1.118567130000000e+02 * invT + 2.321087500000001e-01 -
                 4.017210900000000e+00 * logT - 1.119910065000000e-03 * T +
                 1.056096916666667e-07 * T2 - 9.520530833333334e-12 * T3 +
                 5.395426750000000e-16 * T4;
    // species 7: CH3
    species[7] = +1.677558430000000e+04 * invT - 6.194354070000000e+00 -
                 2.285717720000000e+00 * logT - 3.619950185000000e-03 * T +
                 4.978572466666667e-07 * T2 - 4.964038700000000e-11 * T3 +
                 2.335771970000000e-15 * T4;
    // species 8: CH4
    species[8] = -9.468344590000001e+03 * invT - 1.836246650500000e+01 -
                 7.485149500000000e-02 * logT - 6.695473350000000e-03 * T +
                 9.554763483333333e-07 * T2 - 1.019104458333333e-10 * T3 +
                 5.090761500000000e-15 * T4;
    // species 9: CO
    species[9] = -1.415187240000000e+04 * invT - 5.103502110000000e+00 -
                 2.715185610000000e+00 * logT - 1.031263715000000e-03 * T +
                 1.664709618333334e-07 * T2 - 1.917108400000000e-11 * T3 +
                 1.018238580000000e-15 * T4;
    // species 10: CO2
    species[10] = -4.875916600000000e+04 * invT + 1.585822230000000e+00 -
                  3.857460290000000e+00 * logT - 2.207185130000000e-03 * T +
                  3.691356733333334e-07 * T2 - 4.362418233333334e-11 * T3 +
                  2.360420820000000e-15 * T4;
    // species 11: CH2O
    species[11] = -1.399583230000000e+04 * invT - 1.189563292000000e+01 -
                  1.760690080000000e+00 * logT - 4.600000410000000e-03 * T +
                  7.370980216666666e-07 * T2 - 8.386767666666666e-11 * T3 +
                  4.419278200000001e-15 * T4;
    // species 12: N2
    species[12] = -9.227977000000000e+02 * invT - 3.053888000000000e+00 -
                  2.926640000000000e+00 * logT - 7.439884000000000e-04 * T +
                  9.474600000000001e-08 * T2 - 8.414198333333333e-12 * T3 +
                  3.376675500000000e-16 * T4;
  }
}

// compute the g/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
gibbs_qss(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: CH2
    species[0] = +4.600404010000000e+04 * invT + 2.200146820000000e+00 -
                 3.762678670000000e+00 * logT - 4.844360715000000e-04 * T -
                 4.658164016666667e-07 * T2 + 3.209092941666667e-10 * T3 -
                 8.437085950000000e-14 * T4;
    // species 1: CH2(S)
    species[1] = +5.049681630000000e+04 * invT + 4.967723077000000e+00 -
                 4.198604110000000e+00 * logT + 1.183307095000000e-03 * T -
                 1.372160366666667e-06 * T2 + 5.573466508333334e-10 * T3 -
                 9.715736850000000e-14 * T4;
    // species 2: HCO
    species[2] = +3.839564960000000e+03 * invT + 8.268134100000002e-01 -
                 4.221185840000000e+00 * logT + 1.621962660000000e-03 * T -
                 2.296657433333333e-06 * T2 + 1.109534108333333e-09 * T3 -
                 2.168844325000000e-13 * T4;
    // species 3: CH2OH
    species[3] = -3.193913670000000e+03 * invT - 1.609133250000000e+00 -
                 3.863889180000000e+00 * logT - 2.798361520000000e-03 * T -
                 9.887863183333334e-07 * T2 + 8.711001000000001e-10 * T3 -
                 2.184836390000000e-13 * T4;
  } else {
    // species 0: CH2
    species[0] = +4.626360400000000e+04 * invT - 3.297092110000000e+00 -
                 2.874101130000000e+00 * logT - 1.828196460000000e-03 * T +
                 2.348243283333333e-07 * T2 - 2.168162908333333e-11 * T3 +
                 9.386378350000000e-16 * T4;
    // species 1: CH2(S)
    species[1] = +5.092599970000000e+04 * invT - 6.334463270000000e+00 -
                 2.292038420000000e+00 * logT - 2.327943185000000e-03 * T +
                 3.353199116666667e-07 * T2 - 3.482550000000000e-11 * T3 +
                 1.698581825000000e-15 * T4;
    // species 2: HCO
    species[2] = +4.011918150000000e+03 * invT - 7.026170540000000e+00 -
                 2.772174380000000e+00 * logT - 2.478477630000000e-03 * T +
                 4.140760216666667e-07 * T2 - 4.909681483333334e-11 * T3 +
                 2.667543555000000e-15 * T4;
    // species 3: CH2OH
    species[3] = -3.242506270000000e+03 * invT - 2.117766459999999e+00 -
                 3.692665690000000e+00 * logT - 4.322883985000000e-03 * T +
                 6.251685333333334e-07 * T2 - 6.560288633333334e-11 * T3 +
                 3.242771005000000e-15 * T4;
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
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: H2
    species[0] = -9.17935173e+02 * invT + 6.61320882e-01 -
                 2.34433112e+00 * logT - 3.99026037e-03 * T +
                 3.24635850e-06 * T2 - 1.67976745e-09 * T3 +
                 3.68805881e-13 * T4;
    // species 1: H
    species[1] = +2.54736599e+04 * invT + 1.94668285e+00 -
                 2.50000000e+00 * logT - 3.52666409e-13 * T +
                 3.32653273e-16 * T2 - 1.91734693e-19 * T3 +
                 4.63866166e-23 * T4;
    // species 2: O
    species[2] = +2.91222592e+04 * invT + 1.16333640e-01 -
                 3.16826710e+00 * logT + 1.63965942e-03 * T -
                 1.10717733e-06 * T2 + 5.10672187e-10 * T3 -
                 1.05632985e-13 * T4;
    // species 3: O2
    species[3] = -1.06394356e+03 * invT - 8.75219370e-01 -
                 3.78245636e+00 * logT + 1.49836708e-03 * T -
                 1.64121700e-06 * T2 + 8.06774591e-10 * T3 -
                 1.62186418e-13 * T4;
    // species 4: OH
    species[4] = +3.61508056e+03 * invT + 3.09594089e+00 -
                 3.99201543e+00 * logT + 1.20065876e-03 * T -
                 7.69656402e-07 * T2 + 3.23427778e-10 * T3 -
                 6.82057350e-14 * T4;
    // species 5: H2O
    species[5] = -3.02937267e+04 * invT + 4.04767277e+00 -
                 4.19864056e+00 * logT + 1.01821705e-03 * T -
                 1.08673369e-06 * T2 + 4.57330885e-10 * T3 -
                 8.85989085e-14 * T4;
    // species 6: HO2
    species[6] = +2.94808040e+02 * invT - 4.14864440e-01 -
                 4.30179801e+00 * logT + 2.37456025e-03 * T -
                 3.52638152e-06 * T2 + 2.02303245e-09 * T3 -
                 4.64612562e-13 * T4;
    // species 7: CH3
    species[7] = +1.64449988e+04 * invT + 1.06902607e+00 -
                 3.67359040e+00 * logT - 1.00547588e-03 * T -
                 9.55036427e-07 * T2 + 5.72597854e-10 * T3 -
                 1.27192867e-13 * T4;
    // species 8: CH4
    species[8] = -1.02466476e+04 * invT + 8.79117989e+00 -
                 5.14987613e+00 * logT + 6.83548940e-03 * T -
                 8.19667665e-06 * T2 + 4.03952522e-09 * T3 -
                 8.33469780e-13 * T4;
    // species 9: CO
    species[9] = -1.43440860e+04 * invT - 9.28875810e-01 -
                 3.57953347e+00 * logT + 3.05176840e-04 * T -
                 1.69469055e-07 * T2 - 7.55838237e-11 * T3 +
                 4.52212249e-14 * T4;
    // species 10: CO2
    species[10] = -4.83719697e+04 * invT - 8.54427870e+00 -
                  2.35677352e+00 * logT - 4.49229839e-03 * T +
                  1.18726045e-06 * T2 - 2.04932518e-10 * T3 +
                  7.18497740e-15 * T4;
    // species 11: CH2O
    species[11] = -1.43089567e+04 * invT + 3.19091025e+00 -
                  4.79372315e+00 * logT + 4.95416684e-03 * T -
                  6.22033347e-06 * T2 + 3.16071051e-09 * T3 -
                  6.58863260e-13 * T4;
    // species 12: N2
    species[12] = -1.02089990e+03 * invT - 1.65169500e+00 -
                  3.29867700e+00 * logT - 7.04120200e-04 * T +
                  6.60537000e-07 * T2 - 4.70126250e-10 * T3 +
                  1.22242700e-13 * T4;
  } else {
    // species 0: H2
    species[0] = -9.50158922e+02 * invT + 5.54230251e+00 -
                 3.33727920e+00 * logT + 2.47012365e-05 * T -
                 8.32427963e-08 * T2 + 1.49638662e-11 * T3 -
                 1.00127688e-15 * T4;
    // species 1: H
    species[1] = +2.54736599e+04 * invT + 1.94668292e+00 -
                 2.50000001e+00 * logT + 1.15421486e-11 * T -
                 2.69269913e-15 * T2 + 3.94596029e-19 * T3 -
                 2.49098679e-23 * T4;
    // species 2: O
    species[2] = +2.92175791e+04 * invT - 3.21491786e+00 -
                 2.56942078e+00 * logT + 4.29870569e-05 * T -
                 6.99140982e-09 * T2 + 8.34814992e-13 * T3 -
                 6.14168455e-17 * T4;
    // species 3: O2
    species[3] = -1.08845772e+03 * invT - 3.17069345e+00 -
                 3.28253784e+00 * logT - 7.41543770e-04 * T +
                 1.26327778e-07 * T2 - 1.74558796e-11 * T3 +
                 1.08358897e-15 * T4;
    // species 4: OH
    species[4] = +3.85865700e+03 * invT - 2.38380843e+00 -
                 3.09288767e+00 * logT - 2.74214858e-04 * T -
                 2.10842047e-08 * T2 + 7.32884630e-12 * T3 -
                 5.87061880e-16 * T4;
    // species 5: H2O
    species[5] = -3.00042971e+04 * invT - 2.93277761e+00 -
                 3.03399249e+00 * logT - 1.08845902e-03 * T +
                 2.73454197e-08 * T2 + 8.08683225e-12 * T3 -
                 8.41004960e-16 * T4;
    // species 6: HO2
    species[6] = +1.11856713e+02 * invT - 7.67891250e-01 -
                 4.01721090e+00 * logT - 1.11991006e-03 * T +
                 1.05609692e-07 * T2 - 9.52053083e-12 * T3 +
                 5.39542675e-16 * T4;
    // species 7: CH3
    species[7] = +1.67755843e+04 * invT - 7.19435407e+00 -
                 2.28571772e+00 * logT - 3.61995018e-03 * T +
                 4.97857247e-07 * T2 - 4.96403870e-11 * T3 +
                 2.33577197e-15 * T4;
    // species 8: CH4
    species[8] = -9.46834459e+03 * invT - 1.93624665e+01 -
                 7.48514950e-02 * logT - 6.69547335e-03 * T +
                 9.55476348e-07 * T2 - 1.01910446e-10 * T3 +
                 5.09076150e-15 * T4;
    // species 9: CO
    species[9] = -1.41518724e+04 * invT - 6.10350211e+00 -
                 2.71518561e+00 * logT - 1.03126372e-03 * T +
                 1.66470962e-07 * T2 - 1.91710840e-11 * T3 +
                 1.01823858e-15 * T4;
    // species 10: CO2
    species[10] = -4.87591660e+04 * invT + 5.85822230e-01 -
                  3.85746029e+00 * logT - 2.20718513e-03 * T +
                  3.69135673e-07 * T2 - 4.36241823e-11 * T3 +
                  2.36042082e-15 * T4;
    // species 11: CH2O
    species[11] = -1.39958323e+04 * invT - 1.28956329e+01 -
                  1.76069008e+00 * logT - 4.60000041e-03 * T +
                  7.37098022e-07 * T2 - 8.38676767e-11 * T3 +
                  4.41927820e-15 * T4;
    // species 12: N2
    species[12] = -9.22797700e+02 * invT - 4.05388800e+00 -
                  2.92664000e+00 * logT - 7.43988400e-04 * T +
                  9.47460000e-08 * T2 - 8.41419833e-12 * T3 +
                  3.37667550e-16 * T4;
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

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: H2
    species[0] = +1.34433112e+00 + 3.99026037e-03 * T - 6.49271700e-06 * T2 +
                 5.03930235e-09 * T3 - 1.47522352e-12 * T4 -
                 9.17935173e+02 * invT;
    // species 1: H
    species[1] = +1.50000000e+00 + 3.52666409e-13 * T - 6.65306547e-16 * T2 +
                 5.75204080e-19 * T3 - 1.85546466e-22 * T4 +
                 2.54736599e+04 * invT;
    // species 2: O
    species[2] = +2.16826710e+00 - 1.63965942e-03 * T + 2.21435465e-06 * T2 -
                 1.53201656e-09 * T3 + 4.22531942e-13 * T4 +
                 2.91222592e+04 * invT;
    // species 3: O2
    species[3] = +2.78245636e+00 - 1.49836708e-03 * T + 3.28243400e-06 * T2 -
                 2.42032377e-09 * T3 + 6.48745674e-13 * T4 -
                 1.06394356e+03 * invT;
    // species 4: OH
    species[4] = +2.99201543e+00 - 1.20065876e-03 * T + 1.53931280e-06 * T2 -
                 9.70283332e-10 * T3 + 2.72822940e-13 * T4 +
                 3.61508056e+03 * invT;
    // species 5: H2O
    species[5] = +3.19864056e+00 - 1.01821705e-03 * T + 2.17346737e-06 * T2 -
                 1.37199266e-09 * T3 + 3.54395634e-13 * T4 -
                 3.02937267e+04 * invT;
    // species 6: HO2
    species[6] = +3.30179801e+00 - 2.37456025e-03 * T + 7.05276303e-06 * T2 -
                 6.06909735e-09 * T3 + 1.85845025e-12 * T4 +
                 2.94808040e+02 * invT;
    // species 7: CH3
    species[7] = +2.67359040e+00 + 1.00547588e-03 * T + 1.91007285e-06 * T2 -
                 1.71779356e-09 * T3 + 5.08771468e-13 * T4 +
                 1.64449988e+04 * invT;
    // species 8: CH4
    species[8] = +4.14987613e+00 - 6.83548940e-03 * T + 1.63933533e-05 * T2 -
                 1.21185757e-08 * T3 + 3.33387912e-12 * T4 -
                 1.02466476e+04 * invT;
    // species 9: CO
    species[9] = +2.57953347e+00 - 3.05176840e-04 * T + 3.38938110e-07 * T2 +
                 2.26751471e-10 * T3 - 1.80884900e-13 * T4 -
                 1.43440860e+04 * invT;
    // species 10: CO2
    species[10] = +1.35677352e+00 + 4.49229839e-03 * T - 2.37452090e-06 * T2 +
                  6.14797555e-10 * T3 - 2.87399096e-14 * T4 -
                  4.83719697e+04 * invT;
    // species 11: CH2O
    species[11] = +3.79372315e+00 - 4.95416684e-03 * T + 1.24406669e-05 * T2 -
                  9.48213152e-09 * T3 + 2.63545304e-12 * T4 -
                  1.43089567e+04 * invT;
    // species 12: N2
    species[12] = +2.29867700e+00 + 7.04120200e-04 * T - 1.32107400e-06 * T2 +
                  1.41037875e-09 * T3 - 4.88970800e-13 * T4 -
                  1.02089990e+03 * invT;
  } else {
    // species 0: H2
    species[0] = +2.33727920e+00 - 2.47012365e-05 * T + 1.66485593e-07 * T2 -
                 4.48915985e-11 * T3 + 4.00510752e-15 * T4 -
                 9.50158922e+02 * invT;
    // species 1: H
    species[1] = +1.50000001e+00 - 1.15421486e-11 * T + 5.38539827e-15 * T2 -
                 1.18378809e-18 * T3 + 9.96394714e-23 * T4 +
                 2.54736599e+04 * invT;
    // species 2: O
    species[2] = +1.56942078e+00 - 4.29870569e-05 * T + 1.39828196e-08 * T2 -
                 2.50444497e-12 * T3 + 2.45667382e-16 * T4 +
                 2.92175791e+04 * invT;
    // species 3: O2
    species[3] = +2.28253784e+00 + 7.41543770e-04 * T - 2.52655556e-07 * T2 +
                 5.23676387e-11 * T3 - 4.33435588e-15 * T4 -
                 1.08845772e+03 * invT;
    // species 4: OH
    species[4] = +2.09288767e+00 + 2.74214858e-04 * T + 4.21684093e-08 * T2 -
                 2.19865389e-11 * T3 + 2.34824752e-15 * T4 +
                 3.85865700e+03 * invT;
    // species 5: H2O
    species[5] = +2.03399249e+00 + 1.08845902e-03 * T - 5.46908393e-08 * T2 -
                 2.42604967e-11 * T3 + 3.36401984e-15 * T4 -
                 3.00042971e+04 * invT;
    // species 6: HO2
    species[6] = +3.01721090e+00 + 1.11991006e-03 * T - 2.11219383e-07 * T2 +
                 2.85615925e-11 * T3 - 2.15817070e-15 * T4 +
                 1.11856713e+02 * invT;
    // species 7: CH3
    species[7] = +1.28571772e+00 + 3.61995018e-03 * T - 9.95714493e-07 * T2 +
                 1.48921161e-10 * T3 - 9.34308788e-15 * T4 +
                 1.67755843e+04 * invT;
    // species 8: CH4
    species[8] = -9.25148505e-01 + 6.69547335e-03 * T - 1.91095270e-06 * T2 +
                 3.05731338e-10 * T3 - 2.03630460e-14 * T4 -
                 9.46834459e+03 * invT;
    // species 9: CO
    species[9] = +1.71518561e+00 + 1.03126372e-03 * T - 3.32941924e-07 * T2 +
                 5.75132520e-11 * T3 - 4.07295432e-15 * T4 -
                 1.41518724e+04 * invT;
    // species 10: CO2
    species[10] = +2.85746029e+00 + 2.20718513e-03 * T - 7.38271347e-07 * T2 +
                  1.30872547e-10 * T3 - 9.44168328e-15 * T4 -
                  4.87591660e+04 * invT;
    // species 11: CH2O
    species[11] = +7.60690080e-01 + 4.60000041e-03 * T - 1.47419604e-06 * T2 +
                  2.51603030e-10 * T3 - 1.76771128e-14 * T4 -
                  1.39958323e+04 * invT;
    // species 12: N2
    species[12] = +1.92664000e+00 + 7.43988400e-04 * T - 1.89492000e-07 * T2 +
                  2.52425950e-11 * T3 - 1.35067020e-15 * T4 -
                  9.22797700e+02 * invT;
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

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: H2
    species[0] = +2.34433112e+00 + 3.99026037e-03 * T - 6.49271700e-06 * T2 +
                 5.03930235e-09 * T3 - 1.47522352e-12 * T4 -
                 9.17935173e+02 * invT;
    // species 1: H
    species[1] = +2.50000000e+00 + 3.52666409e-13 * T - 6.65306547e-16 * T2 +
                 5.75204080e-19 * T3 - 1.85546466e-22 * T4 +
                 2.54736599e+04 * invT;
    // species 2: O
    species[2] = +3.16826710e+00 - 1.63965942e-03 * T + 2.21435465e-06 * T2 -
                 1.53201656e-09 * T3 + 4.22531942e-13 * T4 +
                 2.91222592e+04 * invT;
    // species 3: O2
    species[3] = +3.78245636e+00 - 1.49836708e-03 * T + 3.28243400e-06 * T2 -
                 2.42032377e-09 * T3 + 6.48745674e-13 * T4 -
                 1.06394356e+03 * invT;
    // species 4: OH
    species[4] = +3.99201543e+00 - 1.20065876e-03 * T + 1.53931280e-06 * T2 -
                 9.70283332e-10 * T3 + 2.72822940e-13 * T4 +
                 3.61508056e+03 * invT;
    // species 5: H2O
    species[5] = +4.19864056e+00 - 1.01821705e-03 * T + 2.17346737e-06 * T2 -
                 1.37199266e-09 * T3 + 3.54395634e-13 * T4 -
                 3.02937267e+04 * invT;
    // species 6: HO2
    species[6] = +4.30179801e+00 - 2.37456025e-03 * T + 7.05276303e-06 * T2 -
                 6.06909735e-09 * T3 + 1.85845025e-12 * T4 +
                 2.94808040e+02 * invT;
    // species 7: CH3
    species[7] = +3.67359040e+00 + 1.00547588e-03 * T + 1.91007285e-06 * T2 -
                 1.71779356e-09 * T3 + 5.08771468e-13 * T4 +
                 1.64449988e+04 * invT;
    // species 8: CH4
    species[8] = +5.14987613e+00 - 6.83548940e-03 * T + 1.63933533e-05 * T2 -
                 1.21185757e-08 * T3 + 3.33387912e-12 * T4 -
                 1.02466476e+04 * invT;
    // species 9: CO
    species[9] = +3.57953347e+00 - 3.05176840e-04 * T + 3.38938110e-07 * T2 +
                 2.26751471e-10 * T3 - 1.80884900e-13 * T4 -
                 1.43440860e+04 * invT;
    // species 10: CO2
    species[10] = +2.35677352e+00 + 4.49229839e-03 * T - 2.37452090e-06 * T2 +
                  6.14797555e-10 * T3 - 2.87399096e-14 * T4 -
                  4.83719697e+04 * invT;
    // species 11: CH2O
    species[11] = +4.79372315e+00 - 4.95416684e-03 * T + 1.24406669e-05 * T2 -
                  9.48213152e-09 * T3 + 2.63545304e-12 * T4 -
                  1.43089567e+04 * invT;
    // species 12: N2
    species[12] = +3.29867700e+00 + 7.04120200e-04 * T - 1.32107400e-06 * T2 +
                  1.41037875e-09 * T3 - 4.88970800e-13 * T4 -
                  1.02089990e+03 * invT;
  } else {
    // species 0: H2
    species[0] = +3.33727920e+00 - 2.47012365e-05 * T + 1.66485593e-07 * T2 -
                 4.48915985e-11 * T3 + 4.00510752e-15 * T4 -
                 9.50158922e+02 * invT;
    // species 1: H
    species[1] = +2.50000001e+00 - 1.15421486e-11 * T + 5.38539827e-15 * T2 -
                 1.18378809e-18 * T3 + 9.96394714e-23 * T4 +
                 2.54736599e+04 * invT;
    // species 2: O
    species[2] = +2.56942078e+00 - 4.29870569e-05 * T + 1.39828196e-08 * T2 -
                 2.50444497e-12 * T3 + 2.45667382e-16 * T4 +
                 2.92175791e+04 * invT;
    // species 3: O2
    species[3] = +3.28253784e+00 + 7.41543770e-04 * T - 2.52655556e-07 * T2 +
                 5.23676387e-11 * T3 - 4.33435588e-15 * T4 -
                 1.08845772e+03 * invT;
    // species 4: OH
    species[4] = +3.09288767e+00 + 2.74214858e-04 * T + 4.21684093e-08 * T2 -
                 2.19865389e-11 * T3 + 2.34824752e-15 * T4 +
                 3.85865700e+03 * invT;
    // species 5: H2O
    species[5] = +3.03399249e+00 + 1.08845902e-03 * T - 5.46908393e-08 * T2 -
                 2.42604967e-11 * T3 + 3.36401984e-15 * T4 -
                 3.00042971e+04 * invT;
    // species 6: HO2
    species[6] = +4.01721090e+00 + 1.11991006e-03 * T - 2.11219383e-07 * T2 +
                 2.85615925e-11 * T3 - 2.15817070e-15 * T4 +
                 1.11856713e+02 * invT;
    // species 7: CH3
    species[7] = +2.28571772e+00 + 3.61995018e-03 * T - 9.95714493e-07 * T2 +
                 1.48921161e-10 * T3 - 9.34308788e-15 * T4 +
                 1.67755843e+04 * invT;
    // species 8: CH4
    species[8] = +7.48514950e-02 + 6.69547335e-03 * T - 1.91095270e-06 * T2 +
                 3.05731338e-10 * T3 - 2.03630460e-14 * T4 -
                 9.46834459e+03 * invT;
    // species 9: CO
    species[9] = +2.71518561e+00 + 1.03126372e-03 * T - 3.32941924e-07 * T2 +
                 5.75132520e-11 * T3 - 4.07295432e-15 * T4 -
                 1.41518724e+04 * invT;
    // species 10: CO2
    species[10] = +3.85746029e+00 + 2.20718513e-03 * T - 7.38271347e-07 * T2 +
                  1.30872547e-10 * T3 - 9.44168328e-15 * T4 -
                  4.87591660e+04 * invT;
    // species 11: CH2O
    species[11] = +1.76069008e+00 + 4.60000041e-03 * T - 1.47419604e-06 * T2 +
                  2.51603030e-10 * T3 - 1.76771128e-14 * T4 -
                  1.39958323e+04 * invT;
    // species 12: N2
    species[12] = +2.92664000e+00 + 7.43988400e-04 * T - 1.89492000e-07 * T2 +
                  2.52425950e-11 * T3 - 1.35067020e-15 * T4 -
                  9.22797700e+02 * invT;
  }
}

// compute the h/(RT) at the given temperature (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEnthalpy_qss(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: CH2
    species[0] = +3.76267867e+00 + 4.84436072e-04 * T + 9.31632803e-07 * T2 -
                 9.62727883e-10 * T3 + 3.37483438e-13 * T4 +
                 4.60040401e+04 * invT;
    // species 1: CH2(S)
    species[1] = +4.19860411e+00 - 1.18330710e-03 * T + 2.74432073e-06 * T2 -
                 1.67203995e-09 * T3 + 3.88629474e-13 * T4 +
                 5.04968163e+04 * invT;
    // species 2: HCO
    species[2] = +4.22118584e+00 - 1.62196266e-03 * T + 4.59331487e-06 * T2 -
                 3.32860233e-09 * T3 + 8.67537730e-13 * T4 +
                 3.83956496e+03 * invT;
    // species 3: CH2OH
    species[3] = +3.86388918e+00 + 2.79836152e-03 * T + 1.97757264e-06 * T2 -
                 2.61330030e-09 * T3 + 8.73934556e-13 * T4 -
                 3.19391367e+03 * invT;
  } else {
    // species 0: CH2
    species[0] = +2.87410113e+00 + 1.82819646e-03 * T - 4.69648657e-07 * T2 +
                 6.50448872e-11 * T3 - 3.75455134e-15 * T4 +
                 4.62636040e+04 * invT;
    // species 1: CH2(S)
    species[1] = +2.29203842e+00 + 2.32794318e-03 * T - 6.70639823e-07 * T2 +
                 1.04476500e-10 * T3 - 6.79432730e-15 * T4 +
                 5.09259997e+04 * invT;
    // species 2: HCO
    species[2] = +2.77217438e+00 + 2.47847763e-03 * T - 8.28152043e-07 * T2 +
                 1.47290445e-10 * T3 - 1.06701742e-14 * T4 +
                 4.01191815e+03 * invT;
    // species 3: CH2OH
    species[3] = +3.69266569e+00 + 4.32288399e-03 * T - 1.25033707e-06 * T2 +
                 1.96808659e-10 * T3 - 1.29710840e-14 * T4 -
                 3.24250627e+03 * invT;
  }
}

// compute the S/R at the given temperature (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEntropy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: H2
    species[0] = +2.34433112e+00 * logT + 7.98052075e-03 * T -
                 9.73907550e-06 * T2 + 6.71906980e-09 * T3 -
                 1.84402940e-12 * T4 + 6.83010238e-01;
    // species 1: H
    species[1] = +2.50000000e+00 * logT + 7.05332819e-13 * T -
                 9.97959820e-16 * T2 + 7.66938773e-19 * T3 -
                 2.31933083e-22 * T4 - 4.46682853e-01;
    // species 2: O
    species[2] = +3.16826710e+00 * logT - 3.27931884e-03 * T +
                 3.32153198e-06 * T2 - 2.04268875e-09 * T3 +
                 5.28164927e-13 * T4 + 2.05193346e+00;
    // species 3: O2
    species[3] = +3.78245636e+00 * logT - 2.99673416e-03 * T +
                 4.92365101e-06 * T2 - 3.22709836e-09 * T3 +
                 8.10932092e-13 * T4 + 3.65767573e+00;
    // species 4: OH
    species[4] = +3.99201543e+00 * logT - 2.40131752e-03 * T +
                 2.30896920e-06 * T2 - 1.29371111e-09 * T3 +
                 3.41028675e-13 * T4 - 1.03925458e-01;
    // species 5: H2O
    species[5] = +4.19864056e+00 * logT - 2.03643410e-03 * T +
                 3.26020105e-06 * T2 - 1.82932354e-09 * T3 +
                 4.42994543e-13 * T4 - 8.49032208e-01;
    // species 6: HO2
    species[6] = +4.30179801e+00 * logT - 4.74912051e-03 * T +
                 1.05791445e-05 * T2 - 8.09212980e-09 * T3 +
                 2.32306281e-12 * T4 + 3.71666245e+00;
    // species 7: CH3
    species[7] = +3.67359040e+00 * logT + 2.01095175e-03 * T +
                 2.86510928e-06 * T2 - 2.29039142e-09 * T3 +
                 6.35964335e-13 * T4 + 1.60456433e+00;
    // species 8: CH4
    species[8] = +5.14987613e+00 * logT - 1.36709788e-02 * T +
                 2.45900299e-05 * T2 - 1.61581009e-08 * T3 +
                 4.16734890e-12 * T4 - 4.64130376e+00;
    // species 9: CO
    species[9] = +3.57953347e+00 * logT - 6.10353680e-04 * T +
                 5.08407165e-07 * T2 + 3.02335295e-10 * T3 -
                 2.26106125e-13 * T4 + 3.50840928e+00;
    // species 10: CO2
    species[10] = +2.35677352e+00 * logT + 8.98459677e-03 * T -
                  3.56178134e-06 * T2 + 8.19730073e-10 * T3 -
                  3.59248870e-14 * T4 + 9.90105222e+00;
    // species 11: CH2O
    species[11] = +4.79372315e+00 * logT - 9.90833369e-03 * T +
                  1.86610004e-05 * T2 - 1.26428420e-08 * T3 +
                  3.29431630e-12 * T4 + 6.02812900e-01;
    // species 12: N2
    species[12] = +3.29867700e+00 * logT + 1.40824040e-03 * T -
                  1.98161100e-06 * T2 + 1.88050500e-09 * T3 -
                  6.11213500e-13 * T4 + 3.95037200e+00;
  } else {
    // species 0: H2
    species[0] = +3.33727920e+00 * logT - 4.94024731e-05 * T +
                 2.49728389e-07 * T2 - 5.98554647e-11 * T3 +
                 5.00638440e-15 * T4 - 3.20502331e+00;
    // species 1: H
    species[1] = +2.50000001e+00 * logT - 2.30842973e-11 * T +
                 8.07809740e-15 * T2 - 1.57838412e-18 * T3 +
                 1.24549339e-22 * T4 - 4.46682914e-01;
    // species 2: O
    species[2] = +2.56942078e+00 * logT - 8.59741137e-05 * T +
                 2.09742295e-08 * T2 - 3.33925997e-12 * T3 +
                 3.07084227e-16 * T4 + 4.78433864e+00;
    // species 3: O2
    species[3] = +3.28253784e+00 * logT + 1.48308754e-03 * T -
                 3.78983334e-07 * T2 + 6.98235183e-11 * T3 -
                 5.41794485e-15 * T4 + 5.45323129e+00;
    // species 4: OH
    species[4] = +3.09288767e+00 * logT + 5.48429716e-04 * T +
                 6.32526140e-08 * T2 - 2.93153852e-11 * T3 +
                 2.93530940e-15 * T4 + 4.47669610e+00;
    // species 5: H2O
    species[5] = +3.03399249e+00 * logT + 2.17691804e-03 * T -
                 8.20362590e-08 * T2 - 3.23473290e-11 * T3 +
                 4.20502480e-15 * T4 + 4.96677010e+00;
    // species 6: HO2
    species[6] = +4.01721090e+00 * logT + 2.23982013e-03 * T -
                 3.16829075e-07 * T2 + 3.80821233e-11 * T3 -
                 2.69771337e-15 * T4 + 3.78510215e+00;
    // species 7: CH3
    species[7] = +2.28571772e+00 * logT + 7.23990037e-03 * T -
                 1.49357174e-06 * T2 + 1.98561548e-10 * T3 -
                 1.16788599e-14 * T4 + 8.48007179e+00;
    // species 8: CH4
    species[8] = +7.48514950e-02 * logT + 1.33909467e-02 * T -
                 2.86642905e-06 * T2 + 4.07641783e-10 * T3 -
                 2.54538075e-14 * T4 + 1.84373180e+01;
    // species 9: CO
    species[9] = +2.71518561e+00 * logT + 2.06252743e-03 * T -
                 4.99412886e-07 * T2 + 7.66843360e-11 * T3 -
                 5.09119290e-15 * T4 + 7.81868772e+00;
    // species 10: CO2
    species[10] = +3.85746029e+00 * logT + 4.41437026e-03 * T -
                  1.10740702e-06 * T2 + 1.74496729e-10 * T3 -
                  1.18021041e-14 * T4 + 2.27163806e+00;
    // species 11: CH2O
    species[11] = +1.76069008e+00 * logT + 9.20000082e-03 * T -
                  2.21129406e-06 * T2 + 3.35470707e-10 * T3 -
                  2.20963910e-14 * T4 + 1.36563230e+01;
    // species 12: N2
    species[12] = +2.92664000e+00 * logT + 1.48797680e-03 * T -
                  2.84238000e-07 * T2 + 3.36567933e-11 * T3 -
                  1.68833775e-15 * T4 + 5.98052800e+00;
  }
}

// compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dcvpRdT(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: H2
    species[0] = +7.98052075e-03 - 3.89563020e-05 * T + 6.04716282e-08 * T2 -
                 2.95044704e-11 * T3;
    // species 1: H
    species[1] = +7.05332819e-13 - 3.99183928e-15 * T + 6.90244896e-18 * T2 -
                 3.71092933e-21 * T3;
    // species 2: O
    species[2] = -3.27931884e-03 + 1.32861279e-05 * T - 1.83841987e-08 * T2 +
                 8.45063884e-12 * T3;
    // species 3: O2
    species[3] = -2.99673416e-03 + 1.96946040e-05 * T - 2.90438853e-08 * T2 +
                 1.29749135e-11 * T3;
    // species 4: OH
    species[4] = -2.40131752e-03 + 9.23587682e-06 * T - 1.16434000e-08 * T2 +
                 5.45645880e-12 * T3;
    // species 5: H2O
    species[5] = -2.03643410e-03 + 1.30408042e-05 * T - 1.64639119e-08 * T2 +
                 7.08791268e-12 * T3;
    // species 6: HO2
    species[6] = -4.74912051e-03 + 4.23165782e-05 * T - 7.28291682e-08 * T2 +
                 3.71690050e-11 * T3;
    // species 7: CH3
    species[7] = +2.01095175e-03 + 1.14604371e-05 * T - 2.06135228e-08 * T2 +
                 1.01754294e-11 * T3;
    // species 8: CH4
    species[8] = -1.36709788e-02 + 9.83601198e-05 * T - 1.45422908e-07 * T2 +
                 6.66775824e-11 * T3;
    // species 9: CO
    species[9] = -6.10353680e-04 + 2.03362866e-06 * T + 2.72101765e-09 * T2 -
                 3.61769800e-12 * T3;
    // species 10: CO2
    species[10] = +8.98459677e-03 - 1.42471254e-05 * T + 7.37757066e-09 * T2 -
                  5.74798192e-13 * T3;
    // species 11: CH2O
    species[11] = -9.90833369e-03 + 7.46440016e-05 * T - 1.13785578e-07 * T2 +
                  5.27090608e-11 * T3;
    // species 12: N2
    species[12] = +1.40824040e-03 - 7.92644400e-06 * T + 1.69245450e-08 * T2 -
                  9.77941600e-12 * T3;
  } else {
    // species 0: H2
    species[0] = -4.94024731e-05 + 9.98913556e-07 * T - 5.38699182e-10 * T2 +
                 8.01021504e-14 * T3;
    // species 1: H
    species[1] = -2.30842973e-11 + 3.23123896e-14 * T - 1.42054571e-17 * T2 +
                 1.99278943e-21 * T3;
    // species 2: O
    species[2] = -8.59741137e-05 + 8.38969178e-08 * T - 3.00533397e-11 * T2 +
                 4.91334764e-15 * T3;
    // species 3: O2
    species[3] = +1.48308754e-03 - 1.51593334e-06 * T + 6.28411665e-10 * T2 -
                 8.66871176e-14 * T3;
    // species 4: OH
    species[4] = +5.48429716e-04 + 2.53010456e-07 * T - 2.63838467e-10 * T2 +
                 4.69649504e-14 * T3;
    // species 5: H2O
    species[5] = +2.17691804e-03 - 3.28145036e-07 * T - 2.91125961e-10 * T2 +
                 6.72803968e-14 * T3;
    // species 6: HO2
    species[6] = +2.23982013e-03 - 1.26731630e-06 * T + 3.42739110e-10 * T2 -
                 4.31634140e-14 * T3;
    // species 7: CH3
    species[7] = +7.23990037e-03 - 5.97428696e-06 * T + 1.78705393e-09 * T2 -
                 1.86861758e-13 * T3;
    // species 8: CH4
    species[8] = +1.33909467e-02 - 1.14657162e-05 * T + 3.66877605e-09 * T2 -
                 4.07260920e-13 * T3;
    // species 9: CO
    species[9] = +2.06252743e-03 - 1.99765154e-06 * T + 6.90159024e-10 * T2 -
                 8.14590864e-14 * T3;
    // species 10: CO2
    species[10] = +4.41437026e-03 - 4.42962808e-06 * T + 1.57047056e-09 * T2 -
                  1.88833666e-13 * T3;
    // species 11: CH2O
    species[11] = +9.20000082e-03 - 8.84517626e-06 * T + 3.01923636e-09 * T2 -
                  3.53542256e-13 * T3;
    // species 12: N2
    species[12] = +1.48797680e-03 - 1.13695200e-06 * T + 3.02911140e-10 * T2 -
                  2.70134040e-14 * T3;
  }
}

// Returns the mean specific heat at CP (Eq. 33)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cpbl)
{
  amrex::Real result = 0;
  amrex::Real cpor[13]; //  temporary storage
  cp_R(cpor, T);

  // perform dot product
  for (int id = 0; id < 13; ++id) {
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

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: H2
    result += y[0] *
              (+2.34433112e+00 + 7.98052075e-03 * T - 1.94781510e-05 * T2 +
               2.01572094e-08 * T3 - 7.37611761e-12 * T4) *
              0.4960317460317460;
    // species 1: H
    result += y[1] *
              (+2.50000000e+00 + 7.05332819e-13 * T - 1.99591964e-15 * T2 +
               2.30081632e-18 * T3 - 9.27732332e-22 * T4) *
              0.9920634920634921;
    // species 2: O
    result += y[2] *
              (+3.16826710e+00 - 3.27931884e-03 * T + 6.64306396e-06 * T2 -
               6.12806624e-09 * T3 + 2.11265971e-12 * T4) *
              0.0625039064941559;
    // species 3: O2
    result += y[3] *
              (+3.78245636e+00 - 2.99673416e-03 * T + 9.84730201e-06 * T2 -
               9.68129509e-09 * T3 + 3.24372837e-12 * T4) *
              0.0312519532470779;
    // species 4: OH
    result += y[4] *
              (+3.99201543e+00 - 2.40131752e-03 * T + 4.61793841e-06 * T2 -
               3.88113333e-09 * T3 + 1.36411470e-12 * T4) *
              0.0587993179279120;
    // species 5: H2O
    result += y[5] *
              (+4.19864056e+00 - 2.03643410e-03 * T + 6.52040211e-06 * T2 -
               5.48797062e-09 * T3 + 1.77197817e-12 * T4) *
              0.0555092978073827;
    // species 6: HO2
    result += y[6] *
              (+4.30179801e+00 - 4.74912051e-03 * T + 2.11582891e-05 * T2 -
               2.42763894e-08 * T3 + 9.29225124e-12 * T4) *
              0.0302975216627280;
    // species 7: CH3
    result += y[7] *
              (+3.67359040e+00 + 2.01095175e-03 * T + 5.73021856e-06 * T2 -
               6.87117425e-09 * T3 + 2.54385734e-12 * T4) *
              0.0665114732291320;
    // species 8: CH4
    result += y[8] *
              (+5.14987613e+00 - 1.36709788e-02 * T + 4.91800599e-05 * T2 -
               4.84743026e-08 * T3 + 1.66693956e-11 * T4) *
              0.0623324814560868;
    // species 9: CO
    result += y[9] *
              (+3.57953347e+00 - 6.10353680e-04 * T + 1.01681433e-06 * T2 +
               9.07005884e-10 * T3 - 9.04424499e-13 * T4) *
              0.0357015351660121;
    // species 10: CO2
    result += y[10] *
              (+2.35677352e+00 + 8.98459677e-03 * T - 7.12356269e-06 * T2 +
               2.45919022e-09 * T3 - 1.43699548e-13 * T4) *
              0.0227226249176305;
    // species 11: CH2O
    result += y[11] *
              (+4.79372315e+00 - 9.90833369e-03 * T + 3.73220008e-05 * T2 -
               3.79285261e-08 * T3 + 1.31772652e-11 * T4) *
              0.0333044694598015;
    // species 12: N2
    result += y[12] *
              (+3.29867700e+00 + 1.40824040e-03 * T - 3.96322200e-06 * T2 +
               5.64151500e-09 * T3 - 2.44485400e-12 * T4) *
              0.0356964374955379;
  } else {
    // species 0: H2
    result += y[0] *
              (+3.33727920e+00 - 4.94024731e-05 * T + 4.99456778e-07 * T2 -
               1.79566394e-10 * T3 + 2.00255376e-14 * T4) *
              0.4960317460317460;
    // species 1: H
    result += y[1] *
              (+2.50000001e+00 - 2.30842973e-11 * T + 1.61561948e-14 * T2 -
               4.73515235e-18 * T3 + 4.98197357e-22 * T4) *
              0.9920634920634921;
    // species 2: O
    result += y[2] *
              (+2.56942078e+00 - 8.59741137e-05 * T + 4.19484589e-08 * T2 -
               1.00177799e-11 * T3 + 1.22833691e-15 * T4) *
              0.0625039064941559;
    // species 3: O2
    result += y[3] *
              (+3.28253784e+00 + 1.48308754e-03 * T - 7.57966669e-07 * T2 +
               2.09470555e-10 * T3 - 2.16717794e-14 * T4) *
              0.0312519532470779;
    // species 4: OH
    result += y[4] *
              (+3.09288767e+00 + 5.48429716e-04 * T + 1.26505228e-07 * T2 -
               8.79461556e-11 * T3 + 1.17412376e-14 * T4) *
              0.0587993179279120;
    // species 5: H2O
    result += y[5] *
              (+3.03399249e+00 + 2.17691804e-03 * T - 1.64072518e-07 * T2 -
               9.70419870e-11 * T3 + 1.68200992e-14 * T4) *
              0.0555092978073827;
    // species 6: HO2
    result += y[6] *
              (+4.01721090e+00 + 2.23982013e-03 * T - 6.33658150e-07 * T2 +
               1.14246370e-10 * T3 - 1.07908535e-14 * T4) *
              0.0302975216627280;
    // species 7: CH3
    result += y[7] *
              (+2.28571772e+00 + 7.23990037e-03 * T - 2.98714348e-06 * T2 +
               5.95684644e-10 * T3 - 4.67154394e-14 * T4) *
              0.0665114732291320;
    // species 8: CH4
    result += y[8] *
              (+7.48514950e-02 + 1.33909467e-02 * T - 5.73285809e-06 * T2 +
               1.22292535e-09 * T3 - 1.01815230e-13 * T4) *
              0.0623324814560868;
    // species 9: CO
    result += y[9] *
              (+2.71518561e+00 + 2.06252743e-03 * T - 9.98825771e-07 * T2 +
               2.30053008e-10 * T3 - 2.03647716e-14 * T4) *
              0.0357015351660121;
    // species 10: CO2
    result += y[10] *
              (+3.85746029e+00 + 4.41437026e-03 * T - 2.21481404e-06 * T2 +
               5.23490188e-10 * T3 - 4.72084164e-14 * T4) *
              0.0227226249176305;
    // species 11: CH2O
    result += y[11] *
              (+1.76069008e+00 + 9.20000082e-03 * T - 4.42258813e-06 * T2 +
               1.00641212e-09 * T3 - 8.83855640e-14 * T4) *
              0.0333044694598015;
    // species 12: N2
    result += y[12] *
              (+2.92664000e+00 + 1.48797680e-03 * T - 5.68476000e-07 * T2 +
               1.00970380e-10 * T3 - 6.75335100e-15 * T4) *
              0.0356964374955379;
  }

  cpbs = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 35)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cvbl)
{
  amrex::Real result = 0;
  amrex::Real cvor[13]; //  temporary storage
  cv_R(cvor, T);

  // perform dot product
  for (int id = 0; id < 13; ++id) {
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

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: H2
    result += y[0] *
              (+1.34433112e+00 + 7.98052075e-03 * T - 1.94781510e-05 * T2 +
               2.01572094e-08 * T3 - 7.37611761e-12 * T4) *
              0.4960317460317460;
    // species 1: H
    result += y[1] *
              (+1.50000000e+00 + 7.05332819e-13 * T - 1.99591964e-15 * T2 +
               2.30081632e-18 * T3 - 9.27732332e-22 * T4) *
              0.9920634920634921;
    // species 2: O
    result += y[2] *
              (+2.16826710e+00 - 3.27931884e-03 * T + 6.64306396e-06 * T2 -
               6.12806624e-09 * T3 + 2.11265971e-12 * T4) *
              0.0625039064941559;
    // species 3: O2
    result += y[3] *
              (+2.78245636e+00 - 2.99673416e-03 * T + 9.84730201e-06 * T2 -
               9.68129509e-09 * T3 + 3.24372837e-12 * T4) *
              0.0312519532470779;
    // species 4: OH
    result += y[4] *
              (+2.99201543e+00 - 2.40131752e-03 * T + 4.61793841e-06 * T2 -
               3.88113333e-09 * T3 + 1.36411470e-12 * T4) *
              0.0587993179279120;
    // species 5: H2O
    result += y[5] *
              (+3.19864056e+00 - 2.03643410e-03 * T + 6.52040211e-06 * T2 -
               5.48797062e-09 * T3 + 1.77197817e-12 * T4) *
              0.0555092978073827;
    // species 6: HO2
    result += y[6] *
              (+3.30179801e+00 - 4.74912051e-03 * T + 2.11582891e-05 * T2 -
               2.42763894e-08 * T3 + 9.29225124e-12 * T4) *
              0.0302975216627280;
    // species 7: CH3
    result += y[7] *
              (+2.67359040e+00 + 2.01095175e-03 * T + 5.73021856e-06 * T2 -
               6.87117425e-09 * T3 + 2.54385734e-12 * T4) *
              0.0665114732291320;
    // species 8: CH4
    result += y[8] *
              (+4.14987613e+00 - 1.36709788e-02 * T + 4.91800599e-05 * T2 -
               4.84743026e-08 * T3 + 1.66693956e-11 * T4) *
              0.0623324814560868;
    // species 9: CO
    result += y[9] *
              (+2.57953347e+00 - 6.10353680e-04 * T + 1.01681433e-06 * T2 +
               9.07005884e-10 * T3 - 9.04424499e-13 * T4) *
              0.0357015351660121;
    // species 10: CO2
    result += y[10] *
              (+1.35677352e+00 + 8.98459677e-03 * T - 7.12356269e-06 * T2 +
               2.45919022e-09 * T3 - 1.43699548e-13 * T4) *
              0.0227226249176305;
    // species 11: CH2O
    result += y[11] *
              (+3.79372315e+00 - 9.90833369e-03 * T + 3.73220008e-05 * T2 -
               3.79285261e-08 * T3 + 1.31772652e-11 * T4) *
              0.0333044694598015;
    // species 12: N2
    result += y[12] *
              (+2.29867700e+00 + 1.40824040e-03 * T - 3.96322200e-06 * T2 +
               5.64151500e-09 * T3 - 2.44485400e-12 * T4) *
              0.0356964374955379;
  } else {
    // species 0: H2
    result += y[0] *
              (+2.33727920e+00 - 4.94024731e-05 * T + 4.99456778e-07 * T2 -
               1.79566394e-10 * T3 + 2.00255376e-14 * T4) *
              0.4960317460317460;
    // species 1: H
    result += y[1] *
              (+1.50000001e+00 - 2.30842973e-11 * T + 1.61561948e-14 * T2 -
               4.73515235e-18 * T3 + 4.98197357e-22 * T4) *
              0.9920634920634921;
    // species 2: O
    result += y[2] *
              (+1.56942078e+00 - 8.59741137e-05 * T + 4.19484589e-08 * T2 -
               1.00177799e-11 * T3 + 1.22833691e-15 * T4) *
              0.0625039064941559;
    // species 3: O2
    result += y[3] *
              (+2.28253784e+00 + 1.48308754e-03 * T - 7.57966669e-07 * T2 +
               2.09470555e-10 * T3 - 2.16717794e-14 * T4) *
              0.0312519532470779;
    // species 4: OH
    result += y[4] *
              (+2.09288767e+00 + 5.48429716e-04 * T + 1.26505228e-07 * T2 -
               8.79461556e-11 * T3 + 1.17412376e-14 * T4) *
              0.0587993179279120;
    // species 5: H2O
    result += y[5] *
              (+2.03399249e+00 + 2.17691804e-03 * T - 1.64072518e-07 * T2 -
               9.70419870e-11 * T3 + 1.68200992e-14 * T4) *
              0.0555092978073827;
    // species 6: HO2
    result += y[6] *
              (+3.01721090e+00 + 2.23982013e-03 * T - 6.33658150e-07 * T2 +
               1.14246370e-10 * T3 - 1.07908535e-14 * T4) *
              0.0302975216627280;
    // species 7: CH3
    result += y[7] *
              (+1.28571772e+00 + 7.23990037e-03 * T - 2.98714348e-06 * T2 +
               5.95684644e-10 * T3 - 4.67154394e-14 * T4) *
              0.0665114732291320;
    // species 8: CH4
    result += y[8] *
              (-9.25148505e-01 + 1.33909467e-02 * T - 5.73285809e-06 * T2 +
               1.22292535e-09 * T3 - 1.01815230e-13 * T4) *
              0.0623324814560868;
    // species 9: CO
    result += y[9] *
              (+1.71518561e+00 + 2.06252743e-03 * T - 9.98825771e-07 * T2 +
               2.30053008e-10 * T3 - 2.03647716e-14 * T4) *
              0.0357015351660121;
    // species 10: CO2
    result += y[10] *
              (+2.85746029e+00 + 4.41437026e-03 * T - 2.21481404e-06 * T2 +
               5.23490188e-10 * T3 - 4.72084164e-14 * T4) *
              0.0227226249176305;
    // species 11: CH2O
    result += y[11] *
              (+7.60690080e-01 + 9.20000082e-03 * T - 4.42258813e-06 * T2 +
               1.00641212e-09 * T3 - 8.83855640e-14 * T4) *
              0.0333044694598015;
    // species 12: N2
    result += y[12] *
              (+1.92664000e+00 + 1.48797680e-03 * T - 5.68476000e-07 * T2 +
               1.00970380e-10 * T3 - 6.75335100e-15 * T4) *
              0.0356964374955379;
  }

  cvbs = result * 8.31446261815324e+07;
}

// Returns the mean enthalpy of the mixture in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBML(const amrex::Real T, const amrex::Real x[], amrex::Real& hbml)
{
  amrex::Real result = 0;
  amrex::Real hml[13];                       //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesEnthalpy(hml, T);

  // perform dot product
  for (int id = 0; id < 13; ++id) {
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

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: H2
    result +=
      y[0] *
      (+2.34433112e+00 + 3.99026037e-03 * T - 6.49271700e-06 * T2 +
       5.03930235e-09 * T3 - 1.47522352e-12 * T4 - 9.17935173e+02 * invT) *
      0.4960317460317460;
    // species 1: H
    result +=
      y[1] *
      (+2.50000000e+00 + 3.52666409e-13 * T - 6.65306547e-16 * T2 +
       5.75204080e-19 * T3 - 1.85546466e-22 * T4 + 2.54736599e+04 * invT) *
      0.9920634920634921;
    // species 2: O
    result +=
      y[2] *
      (+3.16826710e+00 - 1.63965942e-03 * T + 2.21435465e-06 * T2 -
       1.53201656e-09 * T3 + 4.22531942e-13 * T4 + 2.91222592e+04 * invT) *
      0.0625039064941559;
    // species 3: O2
    result +=
      y[3] *
      (+3.78245636e+00 - 1.49836708e-03 * T + 3.28243400e-06 * T2 -
       2.42032377e-09 * T3 + 6.48745674e-13 * T4 - 1.06394356e+03 * invT) *
      0.0312519532470779;
    // species 4: OH
    result +=
      y[4] *
      (+3.99201543e+00 - 1.20065876e-03 * T + 1.53931280e-06 * T2 -
       9.70283332e-10 * T3 + 2.72822940e-13 * T4 + 3.61508056e+03 * invT) *
      0.0587993179279120;
    // species 5: H2O
    result +=
      y[5] *
      (+4.19864056e+00 - 1.01821705e-03 * T + 2.17346737e-06 * T2 -
       1.37199266e-09 * T3 + 3.54395634e-13 * T4 - 3.02937267e+04 * invT) *
      0.0555092978073827;
    // species 6: HO2
    result +=
      y[6] *
      (+4.30179801e+00 - 2.37456025e-03 * T + 7.05276303e-06 * T2 -
       6.06909735e-09 * T3 + 1.85845025e-12 * T4 + 2.94808040e+02 * invT) *
      0.0302975216627280;
    // species 7: CH3
    result +=
      y[7] *
      (+3.67359040e+00 + 1.00547588e-03 * T + 1.91007285e-06 * T2 -
       1.71779356e-09 * T3 + 5.08771468e-13 * T4 + 1.64449988e+04 * invT) *
      0.0665114732291320;
    // species 8: CH4
    result +=
      y[8] *
      (+5.14987613e+00 - 6.83548940e-03 * T + 1.63933533e-05 * T2 -
       1.21185757e-08 * T3 + 3.33387912e-12 * T4 - 1.02466476e+04 * invT) *
      0.0623324814560868;
    // species 9: CO
    result +=
      y[9] *
      (+3.57953347e+00 - 3.05176840e-04 * T + 3.38938110e-07 * T2 +
       2.26751471e-10 * T3 - 1.80884900e-13 * T4 - 1.43440860e+04 * invT) *
      0.0357015351660121;
    // species 10: CO2
    result +=
      y[10] *
      (+2.35677352e+00 + 4.49229839e-03 * T - 2.37452090e-06 * T2 +
       6.14797555e-10 * T3 - 2.87399096e-14 * T4 - 4.83719697e+04 * invT) *
      0.0227226249176305;
    // species 11: CH2O
    result +=
      y[11] *
      (+4.79372315e+00 - 4.95416684e-03 * T + 1.24406669e-05 * T2 -
       9.48213152e-09 * T3 + 2.63545304e-12 * T4 - 1.43089567e+04 * invT) *
      0.0333044694598015;
    // species 12: N2
    result +=
      y[12] *
      (+3.29867700e+00 + 7.04120200e-04 * T - 1.32107400e-06 * T2 +
       1.41037875e-09 * T3 - 4.88970800e-13 * T4 - 1.02089990e+03 * invT) *
      0.0356964374955379;
  } else {
    // species 0: H2
    result +=
      y[0] *
      (+3.33727920e+00 - 2.47012365e-05 * T + 1.66485593e-07 * T2 -
       4.48915985e-11 * T3 + 4.00510752e-15 * T4 - 9.50158922e+02 * invT) *
      0.4960317460317460;
    // species 1: H
    result +=
      y[1] *
      (+2.50000001e+00 - 1.15421486e-11 * T + 5.38539827e-15 * T2 -
       1.18378809e-18 * T3 + 9.96394714e-23 * T4 + 2.54736599e+04 * invT) *
      0.9920634920634921;
    // species 2: O
    result +=
      y[2] *
      (+2.56942078e+00 - 4.29870569e-05 * T + 1.39828196e-08 * T2 -
       2.50444497e-12 * T3 + 2.45667382e-16 * T4 + 2.92175791e+04 * invT) *
      0.0625039064941559;
    // species 3: O2
    result +=
      y[3] *
      (+3.28253784e+00 + 7.41543770e-04 * T - 2.52655556e-07 * T2 +
       5.23676387e-11 * T3 - 4.33435588e-15 * T4 - 1.08845772e+03 * invT) *
      0.0312519532470779;
    // species 4: OH
    result +=
      y[4] *
      (+3.09288767e+00 + 2.74214858e-04 * T + 4.21684093e-08 * T2 -
       2.19865389e-11 * T3 + 2.34824752e-15 * T4 + 3.85865700e+03 * invT) *
      0.0587993179279120;
    // species 5: H2O
    result +=
      y[5] *
      (+3.03399249e+00 + 1.08845902e-03 * T - 5.46908393e-08 * T2 -
       2.42604967e-11 * T3 + 3.36401984e-15 * T4 - 3.00042971e+04 * invT) *
      0.0555092978073827;
    // species 6: HO2
    result +=
      y[6] *
      (+4.01721090e+00 + 1.11991006e-03 * T - 2.11219383e-07 * T2 +
       2.85615925e-11 * T3 - 2.15817070e-15 * T4 + 1.11856713e+02 * invT) *
      0.0302975216627280;
    // species 7: CH3
    result +=
      y[7] *
      (+2.28571772e+00 + 3.61995018e-03 * T - 9.95714493e-07 * T2 +
       1.48921161e-10 * T3 - 9.34308788e-15 * T4 + 1.67755843e+04 * invT) *
      0.0665114732291320;
    // species 8: CH4
    result +=
      y[8] *
      (+7.48514950e-02 + 6.69547335e-03 * T - 1.91095270e-06 * T2 +
       3.05731338e-10 * T3 - 2.03630460e-14 * T4 - 9.46834459e+03 * invT) *
      0.0623324814560868;
    // species 9: CO
    result +=
      y[9] *
      (+2.71518561e+00 + 1.03126372e-03 * T - 3.32941924e-07 * T2 +
       5.75132520e-11 * T3 - 4.07295432e-15 * T4 - 1.41518724e+04 * invT) *
      0.0357015351660121;
    // species 10: CO2
    result +=
      y[10] *
      (+3.85746029e+00 + 2.20718513e-03 * T - 7.38271347e-07 * T2 +
       1.30872547e-10 * T3 - 9.44168328e-15 * T4 - 4.87591660e+04 * invT) *
      0.0227226249176305;
    // species 11: CH2O
    result +=
      y[11] *
      (+1.76069008e+00 + 4.60000041e-03 * T - 1.47419604e-06 * T2 +
       2.51603030e-10 * T3 - 1.76771128e-14 * T4 - 1.39958323e+04 * invT) *
      0.0333044694598015;
    // species 12: N2
    result +=
      y[12] *
      (+2.92664000e+00 + 7.43988400e-04 * T - 1.89492000e-07 * T2 +
       2.52425950e-11 * T3 - 1.35067020e-15 * T4 - 9.22797700e+02 * invT) *
      0.0356964374955379;
  }

  const amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  hbms = result * RT;
}

// get mean internal energy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBML(const amrex::Real T, const amrex::Real x[], amrex::Real& ubml)
{
  amrex::Real result = 0;
  amrex::Real uml[13];                       //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesInternalEnergy(uml, T);

  // perform dot product
  for (int id = 0; id < 13; ++id) {
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

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: H2
    result +=
      y[0] *
      (+1.34433112e+00 + 3.99026037e-03 * T - 6.49271700e-06 * T2 +
       5.03930235e-09 * T3 - 1.47522352e-12 * T4 - 9.17935173e+02 * invT) *
      0.4960317460317460;
    // species 1: H
    result +=
      y[1] *
      (+1.50000000e+00 + 3.52666409e-13 * T - 6.65306547e-16 * T2 +
       5.75204080e-19 * T3 - 1.85546466e-22 * T4 + 2.54736599e+04 * invT) *
      0.9920634920634921;
    // species 2: O
    result +=
      y[2] *
      (+2.16826710e+00 - 1.63965942e-03 * T + 2.21435465e-06 * T2 -
       1.53201656e-09 * T3 + 4.22531942e-13 * T4 + 2.91222592e+04 * invT) *
      0.0625039064941559;
    // species 3: O2
    result +=
      y[3] *
      (+2.78245636e+00 - 1.49836708e-03 * T + 3.28243400e-06 * T2 -
       2.42032377e-09 * T3 + 6.48745674e-13 * T4 - 1.06394356e+03 * invT) *
      0.0312519532470779;
    // species 4: OH
    result +=
      y[4] *
      (+2.99201543e+00 - 1.20065876e-03 * T + 1.53931280e-06 * T2 -
       9.70283332e-10 * T3 + 2.72822940e-13 * T4 + 3.61508056e+03 * invT) *
      0.0587993179279120;
    // species 5: H2O
    result +=
      y[5] *
      (+3.19864056e+00 - 1.01821705e-03 * T + 2.17346737e-06 * T2 -
       1.37199266e-09 * T3 + 3.54395634e-13 * T4 - 3.02937267e+04 * invT) *
      0.0555092978073827;
    // species 6: HO2
    result +=
      y[6] *
      (+3.30179801e+00 - 2.37456025e-03 * T + 7.05276303e-06 * T2 -
       6.06909735e-09 * T3 + 1.85845025e-12 * T4 + 2.94808040e+02 * invT) *
      0.0302975216627280;
    // species 7: CH3
    result +=
      y[7] *
      (+2.67359040e+00 + 1.00547588e-03 * T + 1.91007285e-06 * T2 -
       1.71779356e-09 * T3 + 5.08771468e-13 * T4 + 1.64449988e+04 * invT) *
      0.0665114732291320;
    // species 8: CH4
    result +=
      y[8] *
      (+4.14987613e+00 - 6.83548940e-03 * T + 1.63933533e-05 * T2 -
       1.21185757e-08 * T3 + 3.33387912e-12 * T4 - 1.02466476e+04 * invT) *
      0.0623324814560868;
    // species 9: CO
    result +=
      y[9] *
      (+2.57953347e+00 - 3.05176840e-04 * T + 3.38938110e-07 * T2 +
       2.26751471e-10 * T3 - 1.80884900e-13 * T4 - 1.43440860e+04 * invT) *
      0.0357015351660121;
    // species 10: CO2
    result +=
      y[10] *
      (+1.35677352e+00 + 4.49229839e-03 * T - 2.37452090e-06 * T2 +
       6.14797555e-10 * T3 - 2.87399096e-14 * T4 - 4.83719697e+04 * invT) *
      0.0227226249176305;
    // species 11: CH2O
    result +=
      y[11] *
      (+3.79372315e+00 - 4.95416684e-03 * T + 1.24406669e-05 * T2 -
       9.48213152e-09 * T3 + 2.63545304e-12 * T4 - 1.43089567e+04 * invT) *
      0.0333044694598015;
    // species 12: N2
    result +=
      y[12] *
      (+2.29867700e+00 + 7.04120200e-04 * T - 1.32107400e-06 * T2 +
       1.41037875e-09 * T3 - 4.88970800e-13 * T4 - 1.02089990e+03 * invT) *
      0.0356964374955379;
  } else {
    // species 0: H2
    result +=
      y[0] *
      (+2.33727920e+00 - 2.47012365e-05 * T + 1.66485593e-07 * T2 -
       4.48915985e-11 * T3 + 4.00510752e-15 * T4 - 9.50158922e+02 * invT) *
      0.4960317460317460;
    // species 1: H
    result +=
      y[1] *
      (+1.50000001e+00 - 1.15421486e-11 * T + 5.38539827e-15 * T2 -
       1.18378809e-18 * T3 + 9.96394714e-23 * T4 + 2.54736599e+04 * invT) *
      0.9920634920634921;
    // species 2: O
    result +=
      y[2] *
      (+1.56942078e+00 - 4.29870569e-05 * T + 1.39828196e-08 * T2 -
       2.50444497e-12 * T3 + 2.45667382e-16 * T4 + 2.92175791e+04 * invT) *
      0.0625039064941559;
    // species 3: O2
    result +=
      y[3] *
      (+2.28253784e+00 + 7.41543770e-04 * T - 2.52655556e-07 * T2 +
       5.23676387e-11 * T3 - 4.33435588e-15 * T4 - 1.08845772e+03 * invT) *
      0.0312519532470779;
    // species 4: OH
    result +=
      y[4] *
      (+2.09288767e+00 + 2.74214858e-04 * T + 4.21684093e-08 * T2 -
       2.19865389e-11 * T3 + 2.34824752e-15 * T4 + 3.85865700e+03 * invT) *
      0.0587993179279120;
    // species 5: H2O
    result +=
      y[5] *
      (+2.03399249e+00 + 1.08845902e-03 * T - 5.46908393e-08 * T2 -
       2.42604967e-11 * T3 + 3.36401984e-15 * T4 - 3.00042971e+04 * invT) *
      0.0555092978073827;
    // species 6: HO2
    result +=
      y[6] *
      (+3.01721090e+00 + 1.11991006e-03 * T - 2.11219383e-07 * T2 +
       2.85615925e-11 * T3 - 2.15817070e-15 * T4 + 1.11856713e+02 * invT) *
      0.0302975216627280;
    // species 7: CH3
    result +=
      y[7] *
      (+1.28571772e+00 + 3.61995018e-03 * T - 9.95714493e-07 * T2 +
       1.48921161e-10 * T3 - 9.34308788e-15 * T4 + 1.67755843e+04 * invT) *
      0.0665114732291320;
    // species 8: CH4
    result +=
      y[8] *
      (-9.25148505e-01 + 6.69547335e-03 * T - 1.91095270e-06 * T2 +
       3.05731338e-10 * T3 - 2.03630460e-14 * T4 - 9.46834459e+03 * invT) *
      0.0623324814560868;
    // species 9: CO
    result +=
      y[9] *
      (+1.71518561e+00 + 1.03126372e-03 * T - 3.32941924e-07 * T2 +
       5.75132520e-11 * T3 - 4.07295432e-15 * T4 - 1.41518724e+04 * invT) *
      0.0357015351660121;
    // species 10: CO2
    result +=
      y[10] *
      (+2.85746029e+00 + 2.20718513e-03 * T - 7.38271347e-07 * T2 +
       1.30872547e-10 * T3 - 9.44168328e-15 * T4 - 4.87591660e+04 * invT) *
      0.0227226249176305;
    // species 11: CH2O
    result +=
      y[11] *
      (+7.60690080e-01 + 4.60000041e-03 * T - 1.47419604e-06 * T2 +
       2.51603030e-10 * T3 - 1.76771128e-14 * T4 - 1.39958323e+04 * invT) *
      0.0333044694598015;
    // species 12: N2
    result +=
      y[12] *
      (+1.92664000e+00 + 7.43988400e-04 * T - 1.89492000e-07 * T2 +
       2.52425950e-11 * T3 - 1.35067020e-15 * T4 - 9.22797700e+02 * invT) *
      0.0356964374955379;
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
  amrex::Real sor[13]; //  temporary storage
  speciesEntropy(sor, T);

  // Compute Eq 42
  for (int id = 0; id < 13; ++id) {
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
  amrex::Real sor[13]; //  temporary storage
  amrex::Real x[13];   //  need a ytx conversion
  amrex::Real YOW = 0; // See Eq 4, 6 in CK Manual

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 13; i++) {
    YOW += y[i] * imw(i);
  }
  // Now compute y to x conversion
  x[0] = y[0] / (2.016000 * YOW);
  x[1] = y[1] / (1.008000 * YOW);
  x[2] = y[2] / (15.999000 * YOW);
  x[3] = y[3] / (31.998000 * YOW);
  x[4] = y[4] / (17.007000 * YOW);
  x[5] = y[5] / (18.015000 * YOW);
  x[6] = y[6] / (33.006000 * YOW);
  x[7] = y[7] / (15.035000 * YOW);
  x[8] = y[8] / (16.043000 * YOW);
  x[9] = y[9] / (28.010000 * YOW);
  x[10] = y[10] / (44.009000 * YOW);
  x[11] = y[11] / (30.026000 * YOW);
  x[12] = y[12] / (28.014000 * YOW);
  speciesEntropy(sor, T);
  // Perform computation in Eq 42 and 43
  for (int i = 0; i < 13; i++) {
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
  XW += x[0] * 2.016000;                   // H2
  XW += x[1] * 1.008000;                   // H
  XW += x[2] * 15.999000;                  // O
  XW += x[3] * 31.998000;                  // O2
  XW += x[4] * 17.007000;                  // OH
  XW += x[5] * 18.015000;                  // H2O
  XW += x[6] * 33.006000;                  // HO2
  XW += x[7] * 15.035000;                  // CH3
  XW += x[8] * 16.043000;                  // CH4
  XW += x[9] * 28.010000;                  // CO
  XW += x[10] * 44.009000;                 // CO2
  XW += x[11] * 30.026000;                 // CH2O
  XW += x[12] * 28.014000;                 // N2
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

  for (int i = 0; i < 13; i++) {
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
  W += c[0] * 2.016000;   // H2
  W += c[1] * 1.008000;   // H
  W += c[2] * 15.999000;  // O
  W += c[3] * 31.998000;  // O2
  W += c[4] * 17.007000;  // OH
  W += c[5] * 18.015000;  // H2O
  W += c[6] * 33.006000;  // HO2
  W += c[7] * 15.035000;  // CH3
  W += c[8] * 16.043000;  // CH4
  W += c[9] * 28.010000;  // CO
  W += c[10] * 44.009000; // CO2
  W += c[11] * 30.026000; // CH2O
  W += c[12] * 28.014000; // N2

  for (int id = 0; id < 13; ++id) {
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
  XW += x[0] * 2.016000;                     // H2
  XW += x[1] * 1.008000;                     // H
  XW += x[2] * 15.999000;                    // O
  XW += x[3] * 31.998000;                    // O2
  XW += x[4] * 17.007000;                    // OH
  XW += x[5] * 18.015000;                    // H2O
  XW += x[6] * 33.006000;                    // HO2
  XW += x[7] * 15.035000;                    // CH3
  XW += x[8] * 16.043000;                    // CH4
  XW += x[9] * 28.010000;                    // CO
  XW += x[10] * 44.009000;                   // CO2
  XW += x[11] * 30.026000;                   // CH2O
  XW += x[12] * 28.014000;                   // N2
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

  for (int i = 0; i < 13; i++) {
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
  W += c[0] * 2.016000;   // H2
  W += c[1] * 1.008000;   // H
  W += c[2] * 15.999000;  // O
  W += c[3] * 31.998000;  // O2
  W += c[4] * 17.007000;  // OH
  W += c[5] * 18.015000;  // H2O
  W += c[6] * 33.006000;  // HO2
  W += c[7] * 15.035000;  // CH3
  W += c[8] * 16.043000;  // CH4
  W += c[9] * 28.010000;  // CO
  W += c[10] * 44.009000; // CO2
  W += c[11] * 30.026000; // CH2O
  W += c[12] * 28.014000; // N2

  for (int id = 0; id < 13; ++id) {
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

  for (int i = 0; i < 13; i++) {
    YOW += y[i] * imw(i);
  }

  wtm = 1.0 / YOW;
}

// given x[species]: mole fractions
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWX(const amrex::Real x[], amrex::Real& wtm)
{
  amrex::Real XW = 0;      //  see Eq 4 in CK Manual
  XW += x[0] * 2.016000;   // H2
  XW += x[1] * 1.008000;   // H
  XW += x[2] * 15.999000;  // O
  XW += x[3] * 31.998000;  // O2
  XW += x[4] * 17.007000;  // OH
  XW += x[5] * 18.015000;  // H2O
  XW += x[6] * 33.006000;  // HO2
  XW += x[7] * 15.035000;  // CH3
  XW += x[8] * 16.043000;  // CH4
  XW += x[9] * 28.010000;  // CO
  XW += x[10] * 44.009000; // CO2
  XW += x[11] * 30.026000; // CH2O
  XW += x[12] * 28.014000; // N2
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
  W += c[0] * 2.016000;   // H2
  W += c[1] * 1.008000;   // H
  W += c[2] * 15.999000;  // O
  W += c[3] * 31.998000;  // O2
  W += c[4] * 17.007000;  // OH
  W += c[5] * 18.015000;  // H2O
  W += c[6] * 33.006000;  // HO2
  W += c[7] * 15.035000;  // CH3
  W += c[8] * 16.043000;  // CH4
  W += c[9] * 28.010000;  // CO
  W += c[10] * 44.009000; // CO2
  W += c[11] * 30.026000; // CH2O
  W += c[12] * 28.014000; // N2

  for (int id = 0; id < 13; ++id) {
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

  for (int i = 0; i < 13; i++) {
    YOW += y[i] * imw(i);
  }

  amrex::Real YOWINV = 1.0 / YOW;

  for (int i = 0; i < 13; i++) {
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
  for (int i = 0; i < 13; i++) {
    c[i] = y[i] * imw(i);
  }
  for (int i = 0; i < 13; i++) {
    YOW += c[i];
  }

  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // Now compute conversion

  for (int i = 0; i < 13; i++) {
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

  for (int i = 0; i < 13; i++) {
    c[i] = rho * y[i] * imw(i);
  }
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTY(const amrex::Real x[], amrex::Real y[])
{
  amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
  // Compute mean molecular wt first
  XW += x[0] * 2.016000;   // H2
  XW += x[1] * 1.008000;   // H
  XW += x[2] * 15.999000;  // O
  XW += x[3] * 31.998000;  // O2
  XW += x[4] * 17.007000;  // OH
  XW += x[5] * 18.015000;  // H2O
  XW += x[6] * 33.006000;  // HO2
  XW += x[7] * 15.035000;  // CH3
  XW += x[8] * 16.043000;  // CH4
  XW += x[9] * 28.010000;  // CO
  XW += x[10] * 44.009000; // CO2
  XW += x[11] * 30.026000; // CH2O
  XW += x[12] * 28.014000; // N2
  // Now compute conversion
  amrex::Real XWinv = 1.0 / XW;
  y[0] = x[0] * 2.016000 * XWinv;
  y[1] = x[1] * 1.008000 * XWinv;
  y[2] = x[2] * 15.999000 * XWinv;
  y[3] = x[3] * 31.998000 * XWinv;
  y[4] = x[4] * 17.007000 * XWinv;
  y[5] = x[5] * 18.015000 * XWinv;
  y[6] = x[6] * 33.006000 * XWinv;
  y[7] = x[7] * 15.035000 * XWinv;
  y[8] = x[8] * 16.043000 * XWinv;
  y[9] = x[9] * 28.010000 * XWinv;
  y[10] = x[10] * 44.009000 * XWinv;
  y[11] = x[11] * 30.026000 * XWinv;
  y[12] = x[12] * 28.014000 * XWinv;
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
  for (int id = 0; id < 13; ++id) {
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
  XW += x[0] * 2.016000;   // H2
  XW += x[1] * 1.008000;   // H
  XW += x[2] * 15.999000;  // O
  XW += x[3] * 31.998000;  // O2
  XW += x[4] * 17.007000;  // OH
  XW += x[5] * 18.015000;  // H2O
  XW += x[6] * 33.006000;  // HO2
  XW += x[7] * 15.035000;  // CH3
  XW += x[8] * 16.043000;  // CH4
  XW += x[9] * 28.010000;  // CO
  XW += x[10] * 44.009000; // CO2
  XW += x[11] * 30.026000; // CH2O
  XW += x[12] * 28.014000; // N2
  ROW = rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 13; ++id) {
    c[id] = x[id] * ROW;
  }
}

// convert c[species] (molar conc) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTX(const amrex::Real c[], amrex::Real x[])
{
  amrex::Real sumC = 0;

  // compute sum of c
  for (int id = 0; id < 13; ++id) {
    sumC += c[id];
  }

  //  See Eq 13
  amrex::Real sumCinv = 1.0 / sumC;
  for (int id = 0; id < 13; ++id) {
    x[id] = c[id] * sumCinv;
  }
}

// convert c[species] (molar conc) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTY(const amrex::Real c[], amrex::Real y[])
{
  amrex::Real CW = 0; // See Eq 12 in CK Manual
  // compute denominator in eq 12 first
  CW += c[0] * 2.016000;   // H2
  CW += c[1] * 1.008000;   // H
  CW += c[2] * 15.999000;  // O
  CW += c[3] * 31.998000;  // O2
  CW += c[4] * 17.007000;  // OH
  CW += c[5] * 18.015000;  // H2O
  CW += c[6] * 33.006000;  // HO2
  CW += c[7] * 15.035000;  // CH3
  CW += c[8] * 16.043000;  // CH4
  CW += c[9] * 28.010000;  // CO
  CW += c[10] * 44.009000; // CO2
  CW += c[11] * 30.026000; // CH2O
  CW += c[12] * 28.014000; // N2
  // Now compute conversion
  amrex::Real CWinv = 1.0 / CW;
  y[0] = c[0] * 2.016000 * CWinv;
  y[1] = c[1] * 1.008000 * CWinv;
  y[2] = c[2] * 15.999000 * CWinv;
  y[3] = c[3] * 31.998000 * CWinv;
  y[4] = c[4] * 17.007000 * CWinv;
  y[5] = c[5] * 18.015000 * CWinv;
  y[6] = c[6] * 33.006000 * CWinv;
  y[7] = c[7] * 15.035000 * CWinv;
  y[8] = c[8] * 16.043000 * CWinv;
  y[9] = c[9] * 28.010000 * CWinv;
  y[10] = c[10] * 44.009000 * CWinv;
  y[11] = c[11] * 30.026000 * CWinv;
  y[12] = c[12] * 28.014000 * CWinv;
}

// get specific heat at constant volume as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVML(const amrex::Real T, amrex::Real cvml[])
{
  cv_R(cvml, T);

  // convert to chemkin units
  for (int id = 0; id < 13; ++id) {
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
  for (int id = 0; id < 13; ++id) {
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
  for (int id = 0; id < 13; ++id) {
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
  for (int id = 0; id < 13; ++id) {
    hml[id] *= RT;
  }
}

// Returns the standard-state entropies in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSML(const amrex::Real T, amrex::Real sml[])
{
  speciesEntropy(sml, T);

  // convert to chemkin units
  for (int id = 0; id < 13; ++id) {
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
  cvms[0] *= 4.124237409798234e+07;  // H2
  cvms[1] *= 8.248474819596468e+07;  // H
  cvms[2] *= 5.196863940342046e+06;  // O
  cvms[3] *= 2.598431970171023e+06;  // O2
  cvms[4] *= 4.888847308845322e+06;  // OH
  cvms[5] *= 4.615299815794193e+06;  // H2O
  cvms[6] *= 2.519076112874398e+06;  // HO2
  cvms[7] *= 5.530071578419182e+06;  // CH3
  cvms[8] *= 5.182610869633635e+06;  // CH4
  cvms[9] *= 2.968390795484913e+06;  // CO
  cvms[10] *= 1.889264154639560e+06; // CO2
  cvms[11] *= 2.769087663409458e+06; // CH2O
  cvms[12] *= 2.967966951578939e+06; // N2
}

// Returns the specific heats at constant pressure
// in mass units (Eq. 26)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPMS(const amrex::Real T, amrex::Real cpms[])
{
  cp_R(cpms, T);
  // multiply by R/molecularweight
  cpms[0] *= 4.124237409798234e+07;  // H2
  cpms[1] *= 8.248474819596468e+07;  // H
  cpms[2] *= 5.196863940342046e+06;  // O
  cpms[3] *= 2.598431970171023e+06;  // O2
  cpms[4] *= 4.888847308845322e+06;  // OH
  cpms[5] *= 4.615299815794193e+06;  // H2O
  cpms[6] *= 2.519076112874398e+06;  // HO2
  cpms[7] *= 5.530071578419182e+06;  // CH3
  cpms[8] *= 5.182610869633635e+06;  // CH4
  cpms[9] *= 2.968390795484913e+06;  // CO
  cpms[10] *= 1.889264154639560e+06; // CO2
  cpms[11] *= 2.769087663409458e+06; // CH2O
  cpms[12] *= 2.967966951578939e+06; // N2
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUMS(const amrex::Real T, amrex::Real ums[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  speciesInternalEnergy(ums, T);

  for (int i = 0; i < 13; i++) {
    ums[i] *= RT * imw(i);
  }
}

// Returns enthalpy in mass units (Eq 27.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHMS(const amrex::Real T, amrex::Real hms[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  speciesEnthalpy(hms, T);

  for (int i = 0; i < 13; i++) {
    hms[i] *= RT * imw(i);
  }
}

// Returns the entropies in mass units (Eq 28.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSMS(const amrex::Real T, amrex::Real sms[])
{
  speciesEntropy(sms, T);
  // multiply by R/molecularweight
  sms[0] *= 4.124237409798234e+07;  // H2
  sms[1] *= 8.248474819596468e+07;  // H
  sms[2] *= 5.196863940342046e+06;  // O
  sms[3] *= 2.598431970171023e+06;  // O2
  sms[4] *= 4.888847308845322e+06;  // OH
  sms[5] *= 4.615299815794193e+06;  // H2O
  sms[6] *= 2.519076112874398e+06;  // HO2
  sms[7] *= 5.530071578419182e+06;  // CH3
  sms[8] *= 5.182610869633635e+06;  // CH4
  sms[9] *= 2.968390795484913e+06;  // CO
  sms[10] *= 1.889264154639560e+06; // CO2
  sms[11] *= 2.769087663409458e+06; // CH2O
  sms[12] *= 2.967966951578939e+06; // N2
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_qss_coeff(
  amrex::Real* k_f,
  amrex::Real* qf,
  amrex::Real* qr,
  const amrex::Real* sc,
  const amrex::Real T,
  amrex::Real* g_RT,
  amrex::Real* g_RT_qss)
{
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  const amrex::Real refC = 101325 / 8.31446 * invT;
  const amrex::Real refCinv = 1. / refC;
  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 13; ++i) {
    mixture += sc[i];
  }

  {
    // reaction 0: CH2 + H (+M) <=> CH3 (+M)
    amrex::Real Corr =
      mixture + sc[0] + (6 - 1) * sc[5] + sc[8] + (1.5 - 1) * sc[9] + sc[10];
    const amrex::Real redP =
      Corr / k_f[0] * 3.2e+15 * exp(-3.14 * logT - 618.956498939805 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.32 * exp(-T * 0.0128205128205128) +
      0.68 * exp(-T * 0.00050125313283208) + exp(-5590 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    qf[0] = Corr * k_f[0] * (sc[1]);
    qr[0] = Corr * k_f[0] * exp(-(g_RT[1] - g_RT[7] + g_RT_qss[0])) * (refC) *
            (sc[7]);
  }

  {
    // reaction 2: H + HCO (+M) <=> CH2O (+M)
    amrex::Real Corr =
      mixture + sc[0] + (6 - 1) * sc[5] + sc[8] + (1.5 - 1) * sc[9] + sc[10];
    const amrex::Real redP = Corr / k_f[1] * 1350000000000 *
                             exp(-2.57 * logT - 717.083748771726 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.2176 * exp(-T * 0.003690036900369) +
      0.7824 * exp(-T * 0.000362976406533575) + exp(-6570 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    qf[1] = Corr * k_f[1] * (sc[1]);
    qr[1] = Corr * k_f[1] * exp(-(g_RT[1] - g_RT[11] + g_RT_qss[2])) * (refC) *
            (sc[11]);
  }

  {
    // reaction 3: CH2O + H (+M) <=> CH2OH (+M)
    amrex::Real Corr =
      mixture + sc[0] + (6 - 1) * sc[5] + sc[8] + (1.5 - 1) * sc[9] + sc[10];
    const amrex::Real redP =
      Corr / k_f[2] * 1.27e+20 * exp(-4.82 * logT - 3286.00482770482 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.2813 * exp(-T * 0.00970873786407767) +
      0.7187 * exp(-T * 0.000774593338497289) + exp(-4160 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    qf[2] = Corr * k_f[2] * (sc[1] * sc[11]);
    qr[2] =
      Corr * k_f[2] * exp(-(g_RT[1] + g_RT[11] - g_RT_qss[3])) * (refC) * (1.0);
  }

  {
    // reaction 17: HCO + H2O <=> CO + H + H2O
    qf[3] = k_f[3] * (sc[5]);
    qr[3] = k_f[3] * exp(-(-g_RT[1] - g_RT[9] + g_RT_qss[2])) * (refCinv) *
            (sc[1] * sc[5] * sc[9]);
  }

  {
    // reaction 18: HCO + M <=> CO + H + M
    const amrex::Real Corr =
      mixture + sc[0] + (0 - 1) * sc[5] + sc[8] + (1.5 - 1) * sc[9] + sc[10];
    qf[4] = Corr * k_f[4] * (1.0);
    qr[4] = Corr * k_f[4] * exp(-(-g_RT[1] - g_RT[9] + g_RT_qss[2])) *
            (refCinv) * (sc[1] * sc[9]);
  }

  {
    // reaction 21: CH2 + O <=> H + HCO
    qf[5] = k_f[5] * (sc[2]);
    qr[5] =
      k_f[5] * exp(-(-g_RT[1] + g_RT[2] + g_RT_qss[0] - g_RT_qss[2])) * (sc[1]);
  }

  {
    // reaction 22: CH2(S) + O <=> CO + H2
    qf[6] = k_f[6] * (sc[2]);
    qr[6] = k_f[6] * exp(-(-g_RT[0] + g_RT[2] - g_RT[9] + g_RT_qss[1])) *
            (sc[0] * sc[9]);
  }

  {
    // reaction 25: HCO + O <=> CO + OH
    qf[7] = k_f[7] * (sc[2]);
    qr[7] = k_f[7] * exp(-(g_RT[2] - g_RT[4] - g_RT[9] + g_RT_qss[2])) *
            (sc[4] * sc[9]);
  }

  {
    // reaction 26: HCO + O <=> CO2 + H
    qf[8] = k_f[8] * (sc[2]);
    qr[8] = k_f[8] * exp(-(-g_RT[1] + g_RT[2] - g_RT[10] + g_RT_qss[2])) *
            (sc[1] * sc[10]);
  }

  {
    // reaction 27: CH2O + O <=> HCO + OH
    qf[9] = k_f[9] * (sc[2] * sc[11]);
    qr[9] =
      k_f[9] * exp(-(g_RT[2] - g_RT[4] + g_RT[11] - g_RT_qss[2])) * (sc[4]);
  }

  {
    // reaction 28: CH2OH + O <=> CH2O + OH
    qf[10] = k_f[10] * (sc[2]);
    qr[10] = k_f[10] * exp(-(g_RT[2] - g_RT[4] - g_RT[11] + g_RT_qss[3])) *
             (sc[4] * sc[11]);
  }

  {
    // reaction 30: CH2O + O2 <=> HCO + HO2
    qf[11] = k_f[11] * (sc[3] * sc[11]);
    qr[11] =
      k_f[11] * exp(-(g_RT[3] - g_RT[6] + g_RT[11] - g_RT_qss[2])) * (sc[6]);
  }

  {
    // reaction 36: H + HCO <=> CO + H2
    qf[12] = k_f[12] * (sc[1]);
    qr[12] = k_f[12] * exp(-(-g_RT[0] + g_RT[1] - g_RT[9] + g_RT_qss[2])) *
             (sc[0] * sc[9]);
  }

  {
    // reaction 37: CH2O + H <=> H2 + HCO
    qf[13] = k_f[13] * (sc[1] * sc[11]);
    qr[13] =
      k_f[13] * exp(-(-g_RT[0] + g_RT[1] + g_RT[11] - g_RT_qss[2])) * (sc[0]);
  }

  {
    // reaction 38: CH2OH + H <=> CH2O + H2
    qf[14] = k_f[14] * (sc[1]);
    qr[14] = k_f[14] * exp(-(-g_RT[0] + g_RT[1] - g_RT[11] + g_RT_qss[3])) *
             (sc[0] * sc[11]);
  }

  {
    // reaction 39: CH2OH + H <=> CH3 + OH
    qf[15] = k_f[15] * (sc[1]);
    qr[15] = k_f[15] * exp(-(g_RT[1] - g_RT[4] - g_RT[7] + g_RT_qss[3])) *
             (sc[4] * sc[7]);
  }

  {
    // reaction 40: CH2OH + H <=> CH2(S) + H2O
    qf[16] = k_f[16] * (sc[1]);
    qr[16] =
      k_f[16] * exp(-(g_RT[1] - g_RT[5] - g_RT_qss[1] + g_RT_qss[3])) * (sc[5]);
  }

  {
    // reaction 44: CH2 + OH <=> CH2O + H
    qf[17] = k_f[17] * (sc[4]);
    qr[17] = k_f[17] * exp(-(-g_RT[1] + g_RT[4] - g_RT[11] + g_RT_qss[0])) *
             (sc[1] * sc[11]);
  }

  {
    // reaction 45: CH2(S) + OH <=> CH2O + H
    qf[18] = k_f[18] * (sc[4]);
    qr[18] = k_f[18] * exp(-(-g_RT[1] + g_RT[4] - g_RT[11] + g_RT_qss[1])) *
             (sc[1] * sc[11]);
  }

  {
    // reaction 46: CH3 + OH <=> CH2 + H2O
    qf[19] = k_f[19] * (sc[4] * sc[7]);
    qr[19] =
      k_f[19] * exp(-(g_RT[4] - g_RT[5] + g_RT[7] - g_RT_qss[0])) * (sc[5]);
  }

  {
    // reaction 47: CH3 + OH <=> CH2(S) + H2O
    qf[20] = k_f[20] * (sc[4] * sc[7]);
    qr[20] =
      k_f[20] * exp(-(g_RT[4] - g_RT[5] + g_RT[7] - g_RT_qss[1])) * (sc[5]);
  }

  {
    // reaction 50: HCO + OH <=> CO + H2O
    qf[21] = k_f[21] * (sc[4]);
    qr[21] = k_f[21] * exp(-(g_RT[4] - g_RT[5] - g_RT[9] + g_RT_qss[2])) *
             (sc[5] * sc[9]);
  }

  {
    // reaction 51: CH2O + OH <=> H2O + HCO
    qf[22] = k_f[22] * (sc[4] * sc[11]);
    qr[22] =
      k_f[22] * exp(-(g_RT[4] - g_RT[5] + g_RT[11] - g_RT_qss[2])) * (sc[5]);
  }

  {
    // reaction 52: CH2OH + OH <=> CH2O + H2O
    qf[23] = k_f[23] * (sc[4]);
    qr[23] = k_f[23] * exp(-(g_RT[4] - g_RT[5] - g_RT[11] + g_RT_qss[3])) *
             (sc[5] * sc[11]);
  }

  {
    // reaction 53: CH2 + HO2 <=> CH2O + OH
    qf[24] = k_f[24] * (sc[6]);
    qr[24] = k_f[24] * exp(-(-g_RT[4] + g_RT[6] - g_RT[11] + g_RT_qss[0])) *
             (sc[4] * sc[11]);
  }

  {
    // reaction 56: CH2 + O2 <=> HCO + OH
    qf[25] = k_f[25] * (sc[3]);
    qr[25] =
      k_f[25] * exp(-(g_RT[3] - g_RT[4] + g_RT_qss[0] - g_RT_qss[2])) * (sc[4]);
  }

  {
    // reaction 57: CH2 + H2 <=> CH3 + H
    qf[26] = k_f[26] * (sc[0]);
    qr[26] = k_f[26] * exp(-(g_RT[0] - g_RT[1] - g_RT[7] + g_RT_qss[0])) *
             (sc[1] * sc[7]);
  }

  {
    // reaction 58: CH2 + CH4 <=> 2 CH3
    qf[27] = k_f[27] * (sc[8]);
    qr[27] = k_f[27] * exp(-(-2.000000 * g_RT[7] + g_RT[8] + g_RT_qss[0])) *
             (sc[7] * sc[7]);
  }

  {
    // reaction 59: CH2(S) + N2 <=> CH2 + N2
    qf[28] = k_f[28] * (sc[12]);
    qr[28] = k_f[28] * exp(-(g_RT[12] - g_RT[12] - g_RT_qss[0] + g_RT_qss[1])) *
             (sc[12]);
  }

  {
    // reaction 60: CH2(S) + O2 <=> CO + H + OH
    qf[29] = k_f[29] * (sc[3]);
    qr[29] = k_f[29] *
             exp(-(-g_RT[1] + g_RT[3] - g_RT[4] - g_RT[9] + g_RT_qss[1])) *
             (refCinv) * (sc[1] * sc[4] * sc[9]);
  }

  {
    // reaction 61: CH2(S) + O2 <=> CO + H2O
    qf[30] = k_f[30] * (sc[3]);
    qr[30] = k_f[30] * exp(-(g_RT[3] - g_RT[5] - g_RT[9] + g_RT_qss[1])) *
             (sc[5] * sc[9]);
  }

  {
    // reaction 62: CH2(S) + H2 <=> CH3 + H
    qf[31] = k_f[31] * (sc[0]);
    qr[31] = k_f[31] * exp(-(g_RT[0] - g_RT[1] - g_RT[7] + g_RT_qss[1])) *
             (sc[1] * sc[7]);
  }

  {
    // reaction 63: CH2(S) + H2O <=> CH2 + H2O
    qf[32] = k_f[32] * (sc[5]);
    qr[32] =
      k_f[32] * exp(-(g_RT[5] - g_RT[5] - g_RT_qss[0] + g_RT_qss[1])) * (sc[5]);
  }

  {
    // reaction 64: CH2(S) + CH4 <=> 2 CH3
    qf[33] = k_f[33] * (sc[8]);
    qr[33] = k_f[33] * exp(-(-2.000000 * g_RT[7] + g_RT[8] + g_RT_qss[1])) *
             (sc[7] * sc[7]);
  }

  {
    // reaction 65: CH2(S) + CO <=> CH2 + CO
    qf[34] = k_f[34] * (sc[9]);
    qr[34] =
      k_f[34] * exp(-(g_RT[9] - g_RT[9] - g_RT_qss[0] + g_RT_qss[1])) * (sc[9]);
  }

  {
    // reaction 66: CH2(S) + CO2 <=> CH2 + CO2
    qf[35] = k_f[35] * (sc[10]);
    qr[35] = k_f[35] * exp(-(g_RT[10] - g_RT[10] - g_RT_qss[0] + g_RT_qss[1])) *
             (sc[10]);
  }

  {
    // reaction 67: CH2(S) + CO2 <=> CH2O + CO
    qf[36] = k_f[36] * (sc[10]);
    qr[36] = k_f[36] * exp(-(-g_RT[9] + g_RT[10] - g_RT[11] + g_RT_qss[1])) *
             (sc[9] * sc[11]);
  }

  {
    // reaction 69: CH3 + HCO <=> CH4 + CO
    qf[37] = k_f[37] * (sc[7]);
    qr[37] = k_f[37] * exp(-(g_RT[7] - g_RT[8] - g_RT[9] + g_RT_qss[2])) *
             (sc[8] * sc[9]);
  }

  {
    // reaction 70: CH2O + CH3 <=> CH4 + HCO
    qf[38] = k_f[38] * (sc[7] * sc[11]);
    qr[38] =
      k_f[38] * exp(-(g_RT[7] - g_RT[8] + g_RT[11] - g_RT_qss[2])) * (sc[8]);
  }

  {
    // reaction 71: HCO + O2 <=> CO + HO2
    qf[39] = k_f[39] * (sc[3]);
    qr[39] = k_f[39] * exp(-(g_RT[3] - g_RT[6] - g_RT[9] + g_RT_qss[2])) *
             (sc[6] * sc[9]);
  }

  {
    // reaction 72: CH2OH + O2 <=> CH2O + HO2
    qf[40] = k_f[40] * (sc[3]);
    qr[40] = k_f[40] * exp(-(g_RT[3] - g_RT[6] - g_RT[11] + g_RT_qss[3])) *
             (sc[6] * sc[11]);
  }

  return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_k_f_qss(
  const amrex::Real T,
  const amrex::Real invT,
  const amrex::Real logT,
  amrex::Real* k_f)
{
  // reaction 0: CH2 + H (+M) <=> CH3 (+M)
  k_f[0] = 25000000000 * exp((-0.8) * logT);
  // reaction 2: H + HCO (+M) <=> CH2O (+M)
  k_f[1] = 1090000 * exp((0.48) * logT - (-130.836333109227) * invT);
  // reaction 3: CH2O + H (+M) <=> CH2OH (+M)
  k_f[2] = 540000 * exp((0.454) * logT - (1811.57999689699) * invT);
  // reaction 17: HCO + H2O <=> CO + H + H2O
  k_f[3] = 2244000000000 * exp((-1) * logT - (8554.68331868023) * invT);
  // reaction 18: HCO + M <=> CO + H + M
  k_f[4] = 187000000000 * exp((-1) * logT - (8554.68331868023) * invT);
  // reaction 21: CH2 + O <=> H + HCO
  k_f[5] = 80000000;
  // reaction 22: CH2(S) + O <=> CO + H2
  k_f[6] = 15000000;
  // reaction 25: HCO + O <=> CO + OH
  k_f[7] = 30000000;
  // reaction 26: HCO + O <=> CO2 + H
  k_f[8] = 30000000;
  // reaction 27: CH2O + O <=> HCO + OH
  k_f[9] = 39000000 * exp(-(1781.38699694871) * invT);
  // reaction 28: CH2OH + O <=> CH2O + OH
  k_f[10] = 10000000;
  // reaction 30: CH2O + O2 <=> HCO + HO2
  k_f[11] = 100000000 * exp(-(20128.6666321888) * invT);
  // reaction 36: H + HCO <=> CO + H2
  k_f[12] = 73400000;
  // reaction 37: CH2O + H <=> H2 + HCO
  k_f[13] = 23000 * exp((1.05) * logT - (1648.03458051046) * invT);
  // reaction 38: CH2OH + H <=> CH2O + H2
  k_f[14] = 20000000;
  // reaction 39: CH2OH + H <=> CH3 + OH
  k_f[15] = 12000000;
  // reaction 40: CH2OH + H <=> CH2(S) + H2O
  k_f[16] = 6000000;
  // reaction 44: CH2 + OH <=> CH2O + H
  k_f[17] = 20000000;
  // reaction 45: CH2(S) + OH <=> CH2O + H
  k_f[18] = 30000000;
  // reaction 46: CH3 + OH <=> CH2 + H2O
  k_f[19] = 56 * exp((1.6) * logT - (2727.43432866158) * invT);
  // reaction 47: CH3 + OH <=> CH2(S) + H2O
  k_f[20] = 25010000;
  // reaction 50: HCO + OH <=> CO + H2O
  k_f[21] = 50000000;
  // reaction 51: CH2O + OH <=> H2O + HCO
  k_f[22] = 3430 * exp((1.18) * logT - (-224.93784961471) * invT);
  // reaction 52: CH2OH + OH <=> CH2O + H2O
  k_f[23] = 5000000;
  // reaction 53: CH2 + HO2 <=> CH2O + OH
  k_f[24] = 20000000;
  // reaction 56: CH2 + O2 <=> HCO + OH
  k_f[25] = 13200000 * exp(-(754.82499870708) * invT);
  // reaction 57: CH2 + H2 <=> CH3 + H
  k_f[26] = 0.5 * exp((2) * logT - (3638.25649376812) * invT);
  // reaction 58: CH2 + CH4 <=> 2 CH3
  k_f[27] = 2.46 * exp((2) * logT - (4161.60182620503) * invT);
  // reaction 59: CH2(S) + N2 <=> CH2 + N2
  k_f[28] = 15000000 * exp(-(301.929999482832) * invT);
  // reaction 60: CH2(S) + O2 <=> CO + H + OH
  k_f[29] = 28000000;
  // reaction 61: CH2(S) + O2 <=> CO + H2O
  k_f[30] = 12000000;
  // reaction 62: CH2(S) + H2 <=> CH3 + H
  k_f[31] = 70000000;
  // reaction 63: CH2(S) + H2O <=> CH2 + H2O
  k_f[32] = 30000000;
  // reaction 64: CH2(S) + CH4 <=> 2 CH3
  k_f[33] = 16000000 * exp(-(-286.83349950869) * invT);
  // reaction 65: CH2(S) + CO <=> CH2 + CO
  k_f[34] = 9000000;
  // reaction 66: CH2(S) + CO2 <=> CH2 + CO2
  k_f[35] = 7000000;
  // reaction 67: CH2(S) + CO2 <=> CH2O + CO
  k_f[36] = 14000000;
  // reaction 69: CH3 + HCO <=> CH4 + CO
  k_f[37] = 26480000;
  // reaction 70: CH2O + CH3 <=> CH4 + HCO
  k_f[38] = 0.00332 * exp((2.81) * logT - (2948.84966161566) * invT);
  // reaction 71: HCO + O2 <=> CO + HO2
  k_f[39] = 7600000 * exp(-(201.286666321888) * invT);
  // reaction 72: CH2OH + O2 <=> CH2O + HO2
  k_f[40] = 18000000 * exp(-(452.894999224248) * invT);

  return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_sc_qss(amrex::Real* sc_qss, amrex::Real* qf_co, amrex::Real* qr_co)
{
  amrex::Real epsilon = 1e-12;

  // QSSA coupling between HCO  CH2OH  CH2(S)  CH2
  // QSSA species 2: HCO

  amrex::Real HCO_num = epsilon + 1.0 * qr_co[1] + 1.0 * qr_co[3] +
                        1.0 * qr_co[4] + 1.0 * qr_co[7] + 1.0 * qr_co[8] +
                        1.0 * qf_co[9] + 1.0 * qf_co[11] + 1.0 * qr_co[12] +
                        1.0 * qf_co[13] + 1.0 * qr_co[21] + 1.0 * qf_co[22] +
                        1.0 * qr_co[37] + 1.0 * qf_co[38] + 1.0 * qr_co[39];
  amrex::Real HCO_denom = epsilon - qf_co[1] - qf_co[3] - qf_co[4] - qr_co[5] -
                          qf_co[7] - qf_co[8] - qr_co[9] - qr_co[11] -
                          qf_co[12] - qr_co[13] - qf_co[21] - qr_co[22] -
                          qr_co[25] - qf_co[37] - qr_co[38] - qf_co[39];
  amrex::Real HCO_rhs = -HCO_num / HCO_denom;

  amrex::Real HCO_CH2 =
    (epsilon + 1.0 * qf_co[5] + 1.0 * qf_co[25]) / HCO_denom;

  // QSSA species 3: CH2OH

  amrex::Real CH2OH_num = epsilon + 1.0 * qf_co[2] + 1.0 * qr_co[10] +
                          1.0 * qr_co[14] + 1.0 * qr_co[15] + 1.0 * qr_co[23] +
                          1.0 * qr_co[40];
  amrex::Real CH2OH_denom = epsilon - qr_co[2] - qf_co[10] - qf_co[14] -
                            qf_co[15] - qf_co[16] - qf_co[23] - qf_co[40];
  amrex::Real CH2OH_rhs = -CH2OH_num / CH2OH_denom;

  amrex::Real CH2OH_CH2S = (epsilon + 1.0 * qr_co[16]) / CH2OH_denom;

  // QSSA species 1: CH2(S)

  amrex::Real CH2S_num = epsilon + 1.0 * qr_co[6] + 1.0 * qr_co[18] +
                         1.0 * qf_co[20] + 1.0 * qr_co[29] + 1.0 * qr_co[30] +
                         1.0 * qr_co[31] + 1.0 * qr_co[33] + 1.0 * qr_co[36];
  amrex::Real CH2S_denom = epsilon - qf_co[6] - qr_co[16] - qf_co[18] -
                           qr_co[20] - qf_co[28] - qf_co[29] - qf_co[30] -
                           qf_co[31] - qf_co[32] - qf_co[33] - qf_co[34] -
                           qf_co[35] - qf_co[36];
  amrex::Real CH2S_rhs = -CH2S_num / CH2S_denom;

  amrex::Real CH2S_CH2OH = (epsilon + 1.0 * qf_co[16]) / CH2S_denom;
  amrex::Real CH2S_CH2 = (epsilon + 1.0 * qr_co[28] + 1.0 * qr_co[32] +
                          1.0 * qr_co[34] + 1.0 * qr_co[35]) /
                         CH2S_denom;

  // QSSA species 0: CH2

  amrex::Real CH2_num = epsilon + 1.0 * qr_co[0] + 1.0 * qr_co[17] +
                        1.0 * qf_co[19] + 1.0 * qr_co[24] + 1.0 * qr_co[26] +
                        1.0 * qr_co[27];
  amrex::Real CH2_denom = epsilon - qf_co[0] - qf_co[5] - qf_co[17] -
                          qr_co[19] - qf_co[24] - qf_co[25] - qf_co[26] -
                          qf_co[27] - qr_co[28] - qr_co[32] - qr_co[34] -
                          qr_co[35];
  amrex::Real CH2_rhs = -CH2_num / CH2_denom;

  amrex::Real CH2_HCO =
    (epsilon + 1.0 * qr_co[5] + 1.0 * qr_co[25]) / CH2_denom;
  amrex::Real CH2_CH2S = (epsilon + 1.0 * qf_co[28] + 1.0 * qf_co[32] +
                          1.0 * qf_co[34] + 1.0 * qf_co[35]) /
                         CH2_denom;

  // Putting it all together
  amrex::Real H_0 = CH2_HCO;
  amrex::Real H_1 = CH2S_CH2OH;
  amrex::Real H_2 = CH2_CH2S / (1 - CH2OH_CH2S * H_1);
  sc_qss[0] = ((CH2_rhs - HCO_rhs * H_0) - (CH2S_rhs - CH2OH_rhs * H_1) * H_2) /
              ((1 - HCO_CH2 * H_0) - CH2S_CH2 * H_2);

  sc_qss[1] = ((CH2S_rhs - CH2OH_rhs * H_1) - (CH2S_CH2 * sc_qss[0])) /
              (1 - CH2OH_CH2S * H_1);

  sc_qss[3] = CH2OH_rhs - (+CH2OH_CH2S * sc_qss[1]);

  sc_qss[2] = HCO_rhs - (HCO_CH2 * sc_qss[0]);

  return;
}

// GPU version of productionRate: no more use of thermo namespace vectors
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_qfqr(
  amrex::Real* qf,
  amrex::Real* qr,
  const amrex::Real* sc,
  const amrex::Real* sc_qss,
  const amrex::Real T,
  const amrex::Real invT,
  const amrex::Real logT)
{

  // reaction 0: CH2 + H (+M) <=> CH3 (+M)
  qf[0] = sc[1] * sc_qss[0];
  qr[0] = sc[7];

  // reaction 1: CH3 + H (+M) <=> CH4 (+M)
  qf[1] = sc[1] * sc[7];
  qr[1] = sc[8];

  // reaction 2: H + HCO (+M) <=> CH2O (+M)
  qf[2] = sc[1] * sc_qss[2];
  qr[2] = sc[11];

  // reaction 3: CH2O + H (+M) <=> CH2OH (+M)
  qf[3] = sc[1] * sc[11];
  qr[3] = sc_qss[3];

  // reaction 4: CO + H2 (+M) <=> CH2O (+M)
  qf[4] = sc[0] * sc[9];
  qr[4] = sc[11];

  // reaction 5: 2 O + M <=> O2 + M
  qf[5] = (sc[2] * sc[2]);
  qr[5] = sc[3];

  // reaction 6: H + O + M <=> OH + M
  qf[6] = sc[1] * sc[2];
  qr[6] = sc[4];

  // reaction 7: CO + O + M <=> CO2 + M
  qf[7] = sc[2] * sc[9];
  qr[7] = sc[10];

  // reaction 8: H + O2 + M <=> HO2 + M
  qf[8] = sc[1] * sc[3];
  qr[8] = sc[6];

  // reaction 9: H + O2 + O2 <=> HO2 + O2
  qf[9] = sc[1] * (sc[3] * sc[3]);
  qr[9] = sc[3] * sc[6];

  // reaction 10: H + O2 + H2O <=> HO2 + H2O
  qf[10] = sc[1] * sc[3] * sc[5];
  qr[10] = sc[5] * sc[6];

  // reaction 11: H + O2 + N2 <=> HO2 + N2
  qf[11] = sc[1] * sc[3] * sc[12];
  qr[11] = sc[6] * sc[12];

  // reaction 12: 2 H + M <=> H2 + M
  qf[12] = (sc[1] * sc[1]);
  qr[12] = sc[0];

  // reaction 13: 2 H + H2 <=> H2 + H2
  qf[13] = sc[0] * (sc[1] * sc[1]);
  qr[13] = (sc[0] * sc[0]);

  // reaction 14: 2 H + H2O <=> H2 + H2O
  qf[14] = (sc[1] * sc[1]) * sc[5];
  qr[14] = sc[0] * sc[5];

  // reaction 15: 2 H + CO2 <=> H2 + CO2
  qf[15] = (sc[1] * sc[1]) * sc[10];
  qr[15] = sc[0] * sc[10];

  // reaction 16: H + OH + M <=> H2O + M
  qf[16] = sc[1] * sc[4];
  qr[16] = sc[5];

  // reaction 17: HCO + H2O <=> CO + H + H2O
  qf[17] = sc[5] * sc_qss[2];
  qr[17] = sc[1] * sc[5] * sc[9];

  // reaction 18: HCO + M <=> CO + H + M
  qf[18] = sc_qss[2];
  qr[18] = sc[1] * sc[9];

  // reaction 19: H2 + O <=> H + OH
  qf[19] = sc[0] * sc[2];
  qr[19] = sc[1] * sc[4];

  // reaction 20: HO2 + O <=> O2 + OH
  qf[20] = sc[2] * sc[6];
  qr[20] = sc[3] * sc[4];

  // reaction 21: CH2 + O <=> H + HCO
  qf[21] = sc[2] * sc_qss[0];
  qr[21] = sc[1] * sc_qss[2];

  // reaction 22: CH2(S) + O <=> CO + H2
  qf[22] = sc[2] * sc_qss[1];
  qr[22] = sc[0] * sc[9];

  // reaction 23: CH3 + O <=> CH2O + H
  qf[23] = sc[2] * sc[7];
  qr[23] = sc[1] * sc[11];

  // reaction 24: CH4 + O <=> CH3 + OH
  qf[24] = sc[2] * sc[8];
  qr[24] = sc[4] * sc[7];

  // reaction 25: HCO + O <=> CO + OH
  qf[25] = sc[2] * sc_qss[2];
  qr[25] = sc[4] * sc[9];

  // reaction 26: HCO + O <=> CO2 + H
  qf[26] = sc[2] * sc_qss[2];
  qr[26] = sc[1] * sc[10];

  // reaction 27: CH2O + O <=> HCO + OH
  qf[27] = sc[2] * sc[11];
  qr[27] = sc[4] * sc_qss[2];

  // reaction 28: CH2OH + O <=> CH2O + OH
  qf[28] = sc[2] * sc_qss[3];
  qr[28] = sc[4] * sc[11];

  // reaction 29: CO + O2 <=> CO2 + O
  qf[29] = sc[3] * sc[9];
  qr[29] = sc[2] * sc[10];

  // reaction 30: CH2O + O2 <=> HCO + HO2
  qf[30] = sc[3] * sc[11];
  qr[30] = sc[6] * sc_qss[2];

  // reaction 31: H + O2 <=> O + OH
  qf[31] = sc[1] * sc[3];
  qr[31] = sc[2] * sc[4];

  // reaction 32: H + HO2 <=> H2O + O
  qf[32] = sc[1] * sc[6];
  qr[32] = sc[2] * sc[5];

  // reaction 33: H + HO2 <=> H2 + O2
  qf[33] = sc[1] * sc[6];
  qr[33] = sc[0] * sc[3];

  // reaction 34: H + HO2 <=> 2 OH
  qf[34] = sc[1] * sc[6];
  qr[34] = (sc[4] * sc[4]);

  // reaction 35: CH4 + H <=> CH3 + H2
  qf[35] = sc[1] * sc[8];
  qr[35] = sc[0] * sc[7];

  // reaction 36: H + HCO <=> CO + H2
  qf[36] = sc[1] * sc_qss[2];
  qr[36] = sc[0] * sc[9];

  // reaction 37: CH2O + H <=> H2 + HCO
  qf[37] = sc[1] * sc[11];
  qr[37] = sc[0] * sc_qss[2];

  // reaction 38: CH2OH + H <=> CH2O + H2
  qf[38] = sc[1] * sc_qss[3];
  qr[38] = sc[0] * sc[11];

  // reaction 39: CH2OH + H <=> CH3 + OH
  qf[39] = sc[1] * sc_qss[3];
  qr[39] = sc[4] * sc[7];

  // reaction 40: CH2OH + H <=> CH2(S) + H2O
  qf[40] = sc[1] * sc_qss[3];
  qr[40] = sc[5] * sc_qss[1];

  // reaction 41: H2 + OH <=> H + H2O
  qf[41] = sc[0] * sc[4];
  qr[41] = sc[1] * sc[5];

  // reaction 42: 2 OH <=> H2O + O
  qf[42] = (sc[4] * sc[4]);
  qr[42] = sc[2] * sc[5];

  // reaction 43: HO2 + OH <=> H2O + O2
  qf[43] = sc[4] * sc[6];
  qr[43] = sc[3] * sc[5];

  // reaction 44: CH2 + OH <=> CH2O + H
  qf[44] = sc[4] * sc_qss[0];
  qr[44] = sc[1] * sc[11];

  // reaction 45: CH2(S) + OH <=> CH2O + H
  qf[45] = sc[4] * sc_qss[1];
  qr[45] = sc[1] * sc[11];

  // reaction 46: CH3 + OH <=> CH2 + H2O
  qf[46] = sc[4] * sc[7];
  qr[46] = sc[5] * sc_qss[0];

  // reaction 47: CH3 + OH <=> CH2(S) + H2O
  qf[47] = sc[4] * sc[7];
  qr[47] = sc[5] * sc_qss[1];

  // reaction 48: CH4 + OH <=> CH3 + H2O
  qf[48] = sc[4] * sc[8];
  qr[48] = sc[5] * sc[7];

  // reaction 49: CO + OH <=> CO2 + H
  qf[49] = sc[4] * sc[9];
  qr[49] = sc[1] * sc[10];

  // reaction 50: HCO + OH <=> CO + H2O
  qf[50] = sc[4] * sc_qss[2];
  qr[50] = sc[5] * sc[9];

  // reaction 51: CH2O + OH <=> H2O + HCO
  qf[51] = sc[4] * sc[11];
  qr[51] = sc[5] * sc_qss[2];

  // reaction 52: CH2OH + OH <=> CH2O + H2O
  qf[52] = sc[4] * sc_qss[3];
  qr[52] = sc[5] * sc[11];

  // reaction 53: CH2 + HO2 <=> CH2O + OH
  qf[53] = sc[6] * sc_qss[0];
  qr[53] = sc[4] * sc[11];

  // reaction 54: CH3 + HO2 <=> CH4 + O2
  qf[54] = sc[6] * sc[7];
  qr[54] = sc[3] * sc[8];

  // reaction 55: CO + HO2 <=> CO2 + OH
  qf[55] = sc[6] * sc[9];
  qr[55] = sc[4] * sc[10];

  // reaction 56: CH2 + O2 <=> HCO + OH
  qf[56] = sc[3] * sc_qss[0];
  qr[56] = sc[4] * sc_qss[2];

  // reaction 57: CH2 + H2 <=> CH3 + H
  qf[57] = sc[0] * sc_qss[0];
  qr[57] = sc[1] * sc[7];

  // reaction 58: CH2 + CH4 <=> 2 CH3
  qf[58] = sc[8] * sc_qss[0];
  qr[58] = (sc[7] * sc[7]);

  // reaction 59: CH2(S) + N2 <=> CH2 + N2
  qf[59] = sc[12] * sc_qss[1];
  qr[59] = sc[12] * sc_qss[0];

  // reaction 60: CH2(S) + O2 <=> CO + H + OH
  qf[60] = sc[3] * sc_qss[1];
  qr[60] = sc[1] * sc[4] * sc[9];

  // reaction 61: CH2(S) + O2 <=> CO + H2O
  qf[61] = sc[3] * sc_qss[1];
  qr[61] = sc[5] * sc[9];

  // reaction 62: CH2(S) + H2 <=> CH3 + H
  qf[62] = sc[0] * sc_qss[1];
  qr[62] = sc[1] * sc[7];

  // reaction 63: CH2(S) + H2O <=> CH2 + H2O
  qf[63] = sc[5] * sc_qss[1];
  qr[63] = sc[5] * sc_qss[0];

  // reaction 64: CH2(S) + CH4 <=> 2 CH3
  qf[64] = sc[8] * sc_qss[1];
  qr[64] = (sc[7] * sc[7]);

  // reaction 65: CH2(S) + CO <=> CH2 + CO
  qf[65] = sc[9] * sc_qss[1];
  qr[65] = sc[9] * sc_qss[0];

  // reaction 66: CH2(S) + CO2 <=> CH2 + CO2
  qf[66] = sc[10] * sc_qss[1];
  qr[66] = sc[10] * sc_qss[0];

  // reaction 67: CH2(S) + CO2 <=> CH2O + CO
  qf[67] = sc[10] * sc_qss[1];
  qr[67] = sc[9] * sc[11];

  // reaction 68: CH3 + O2 <=> CH2O + OH
  qf[68] = sc[3] * sc[7];
  qr[68] = sc[4] * sc[11];

  // reaction 69: CH3 + HCO <=> CH4 + CO
  qf[69] = sc[7] * sc_qss[2];
  qr[69] = sc[8] * sc[9];

  // reaction 70: CH2O + CH3 <=> CH4 + HCO
  qf[70] = sc[7] * sc[11];
  qr[70] = sc[8] * sc_qss[2];

  // reaction 71: HCO + O2 <=> CO + HO2
  qf[71] = sc[3] * sc_qss[2];
  qr[71] = sc[6] * sc[9];

  // reaction 72: CH2OH + O2 <=> CH2O + HO2
  qf[72] = sc[3] * sc_qss[3];
  qr[72] = sc[6] * sc[11];

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 13; ++i) {
    mixture += sc[i];
  }

  for (int i = 0; i < 4; ++i) {
    mixture += sc_qss[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[13];
  gibbs(g_RT, T);
  amrex::Real g_RT_qss[4];
  gibbs_qss(g_RT_qss, T);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 * invT;
  amrex::Real refCinv = 1 / refC;

  // Evaluate the kfs
  amrex::Real k_f, Corr;
  amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

  // reaction 0:  CH2 + H (+M) <=> CH3 (+M)
  k_f = 25000000000 * exp((-0.8) * logT);
  Corr = mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
  redP = Corr / k_f * 3.2e+15 * exp(-3.14 * logT - (618.956498939805) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (0.32) * exp(-T * 0.0128205128205128) +
    0.68 * exp(-T * 0.00050125313283208) + exp(-5590 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[0] *= Corr * k_f;
  qr[0] *= Corr * k_f * exp(-(g_RT[1] - g_RT[7] + g_RT_qss[0])) * (refC);
  // reaction 1:  CH3 + H (+M) <=> CH4 (+M)
  k_f = 12700000000 * exp((-0.63) * logT - (192.731983003208) * invT);
  Corr = mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
  redP = Corr / k_f * 2.477e+21 * exp(-4.76 * logT - (1227.84866456352) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (0.217) * exp(-T * 0.0135135135135135) +
    0.783 * exp(-T * 0.000340020401224073) + exp(-6964 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[1] *= Corr * k_f;
  qr[1] *= Corr * k_f * exp(-(g_RT[1] + g_RT[7] - g_RT[8])) * (refC);
  // reaction 2:  H + HCO (+M) <=> CH2O (+M)
  k_f = 1090000 * exp((0.48) * logT - (-130.836333109227) * invT);
  Corr = mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
  redP =
    Corr / k_f * 1350000000000 * exp(-2.57 * logT - (717.083748771726) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (0.2176) * exp(-T * 0.003690036900369) +
    0.7824 * exp(-T * 0.000362976406533575) + exp(-6570 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[2] *= Corr * k_f;
  qr[2] *= Corr * k_f * exp(-(g_RT[1] - g_RT[11] + g_RT_qss[2])) * (refC);
  // reaction 3:  CH2O + H (+M) <=> CH2OH (+M)
  k_f = 540000 * exp((0.454) * logT - (1811.57999689699) * invT);
  Corr = mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
  redP = Corr / k_f * 1.27e+20 * exp(-4.82 * logT - (3286.00482770482) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (0.2813) * exp(-T * 0.00970873786407767) +
    0.7187 * exp(-T * 0.000774593338497289) + exp(-4160 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[3] *= Corr * k_f;
  qr[3] *= Corr * k_f * exp(-(g_RT[1] + g_RT[11] - g_RT_qss[3])) * (refC);
  // reaction 4:  CO + H2 (+M) <=> CH2O (+M)
  k_f = 43 * exp((1.5) * logT - (40056.0465980557) * invT);
  Corr = mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
  redP = Corr / k_f * 5.07e+15 * exp(-3.42 * logT - (42446.3257606281) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (0.0679999999999999) * exp(-T * 0.0050761421319797) +
    0.932 * exp(-T * 0.000649350649350649) + exp(-10300 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[4] *= Corr * k_f;
  qr[4] *= Corr * k_f * exp(-(g_RT[0] + g_RT[9] - g_RT[11])) * (refC);
  // reaction 5:  2 O + M <=> O2 + M
  k_f = 120000 * exp((-1) * logT);
  Corr = mixture + (1.4) * sc[0] + (14.4) * sc[5] + sc[8] + (0.75) * sc[9] +
         (2.6) * sc[10];
  qf[5] *= Corr * k_f;
  qr[5] *= Corr * k_f * exp(-(2.000000 * g_RT[2] - g_RT[3])) * (refC);
  // reaction 6:  H + O + M <=> OH + M
  k_f = 500000 * exp((-1) * logT);
  Corr = mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
  qf[6] *= Corr * k_f;
  qr[6] *= Corr * k_f * exp(-(g_RT[1] + g_RT[2] - g_RT[4])) * (refC);
  // reaction 7:  CO + O + M <=> CO2 + M
  k_f = 602 * exp(-(1509.64999741416) * invT);
  Corr = mixture + sc[0] + (5) * sc[3] + (5) * sc[5] + sc[8] + (0.5) * sc[9] +
         (2.5) * sc[10];
  qf[7] *= Corr * k_f;
  qr[7] *= Corr * k_f * exp(-(g_RT[2] + g_RT[9] - g_RT[10])) * (refC);
  // reaction 8:  H + O2 + M <=> HO2 + M
  k_f = 2800000 * exp((-0.86) * logT);
  Corr = mixture + (-1) * sc[3] + (-1) * sc[5] + (-0.25) * sc[9] +
         (0.5) * sc[10] + (-1) * sc[12];
  qf[8] *= Corr * k_f;
  qr[8] *= Corr * k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC);
  // reaction 9:  H + O2 + O2 <=> HO2 + O2
  k_f = 300000000 * exp((-1.72) * logT);
  qf[9] *= k_f;
  qr[9] *= k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC);
  // reaction 10:  H + O2 + H2O <=> HO2 + H2O
  k_f = 9380000 * exp((-0.76) * logT);
  qf[10] *= k_f;
  qr[10] *= k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC);
  // reaction 11:  H + O2 + N2 <=> HO2 + N2
  k_f = 375000000 * exp((-1.72) * logT);
  qf[11] *= k_f;
  qr[11] *= k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC);
  // reaction 12:  2 H + M <=> H2 + M
  k_f = 1000000 * exp((-1) * logT);
  Corr = mixture + (-1) * sc[0] + (-1) * sc[5] + sc[8] + (-1) * sc[10];
  qf[12] *= Corr * k_f;
  qr[12] *= Corr * k_f * exp(-(-g_RT[0] + 2.000000 * g_RT[1])) * (refC);
  // reaction 13:  2 H + H2 <=> H2 + H2
  k_f = 90000 * exp((-0.6) * logT);
  qf[13] *= k_f;
  qr[13] *= k_f * exp(-(-g_RT[0] + 2.000000 * g_RT[1])) * (refC);
  // reaction 14:  2 H + H2O <=> H2 + H2O
  k_f = 60000000 * exp((-1.25) * logT);
  qf[14] *= k_f;
  qr[14] *= k_f * exp(-(-g_RT[0] + 2.000000 * g_RT[1])) * (refC);
  // reaction 15:  2 H + CO2 <=> H2 + CO2
  k_f = 550000000 * exp((-2) * logT);
  qf[15] *= k_f;
  qr[15] *= k_f * exp(-(-g_RT[0] + 2.000000 * g_RT[1])) * (refC);
  // reaction 16:  H + OH + M <=> H2O + M
  k_f = 22000000000 * exp((-2) * logT);
  Corr = mixture + (-0.27) * sc[0] + (2.65) * sc[5] + sc[8];
  qf[16] *= Corr * k_f;
  qr[16] *= Corr * k_f * exp(-(g_RT[1] + g_RT[4] - g_RT[5])) * (refC);
  // reaction 17:  HCO + H2O <=> CO + H + H2O
  k_f = 2244000000000 * exp((-1) * logT - (8554.68331868023) * invT);
  qf[17] *= k_f;
  qr[17] *= k_f * exp(-(-g_RT[1] - g_RT[9] + g_RT_qss[2])) * (refCinv);
  // reaction 18:  HCO + M <=> CO + H + M
  k_f = 187000000000 * exp((-1) * logT - (8554.68331868023) * invT);
  Corr = mixture + sc[0] + (-1) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
  qf[18] *= Corr * k_f;
  qr[18] *= Corr * k_f * exp(-(-g_RT[1] - g_RT[9] + g_RT_qss[2])) * (refCinv);
  // reaction 19:  H2 + O <=> H + OH
  k_f = 0.05 * exp((2.67) * logT - (3165.23282791169) * invT);
  qf[19] *= k_f;
  qr[19] *= k_f * exp(-(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4]));
  // reaction 20:  HO2 + O <=> O2 + OH
  k_f = 20000000;
  qf[20] *= k_f;
  qr[20] *= k_f * exp(-(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]));
  // reaction 21:  CH2 + O <=> H + HCO
  k_f = 80000000;
  qf[21] *= k_f;
  qr[21] *= k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT_qss[0] - g_RT_qss[2]));
  // reaction 22:  CH2(S) + O <=> CO + H2
  k_f = 15000000;
  qf[22] *= k_f;
  qr[22] *= k_f * exp(-(-g_RT[0] + g_RT[2] - g_RT[9] + g_RT_qss[1]));
  // reaction 23:  CH3 + O <=> CH2O + H
  k_f = 84300000;
  qf[23] *= k_f;
  qr[23] *= k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[11]));
  // reaction 24:  CH4 + O <=> CH3 + OH
  k_f = 1020 * exp((1.5) * logT - (4327.66332592059) * invT);
  qf[24] *= k_f;
  qr[24] *= k_f * exp(-(g_RT[2] - g_RT[4] - g_RT[7] + g_RT[8]));
  // reaction 25:  HCO + O <=> CO + OH
  k_f = 30000000;
  qf[25] *= k_f;
  qr[25] *= k_f * exp(-(g_RT[2] - g_RT[4] - g_RT[9] + g_RT_qss[2]));
  // reaction 26:  HCO + O <=> CO2 + H
  k_f = 30000000;
  qf[26] *= k_f;
  qr[26] *= k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[10] + g_RT_qss[2]));
  // reaction 27:  CH2O + O <=> HCO + OH
  k_f = 39000000 * exp(-(1781.38699694871) * invT);
  qf[27] *= k_f;
  qr[27] *= k_f * exp(-(g_RT[2] - g_RT[4] + g_RT[11] - g_RT_qss[2]));
  // reaction 28:  CH2OH + O <=> CH2O + OH
  k_f = 10000000;
  qf[28] *= k_f;
  qr[28] *= k_f * exp(-(g_RT[2] - g_RT[4] - g_RT[11] + g_RT_qss[3]));
  // reaction 29:  CO + O2 <=> CO2 + O
  k_f = 2500000 * exp(-(24053.7566254656) * invT);
  qf[29] *= k_f;
  qr[29] *= k_f * exp(-(-g_RT[2] + g_RT[3] + g_RT[9] - g_RT[10]));
  // reaction 30:  CH2O + O2 <=> HCO + HO2
  k_f = 100000000 * exp(-(20128.6666321888) * invT);
  qf[30] *= k_f;
  qr[30] *= k_f * exp(-(g_RT[3] - g_RT[6] + g_RT[11] - g_RT_qss[2]));
  // reaction 31:  H + O2 <=> O + OH
  k_f = 83000000 * exp(-(7252.86180424342) * invT);
  qf[31] *= k_f;
  qr[31] *= k_f * exp(-(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4]));
  // reaction 32:  H + HO2 <=> H2O + O
  k_f = 3970000 * exp(-(337.658382754967) * invT);
  qf[32] *= k_f;
  qr[32] *= k_f * exp(-(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6]));
  // reaction 33:  H + HO2 <=> H2 + O2
  k_f = 28000000 * exp(-(537.435399079441) * invT);
  qf[33] *= k_f;
  qr[33] *= k_f * exp(-(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6]));
  // reaction 34:  H + HO2 <=> 2 OH
  k_f = 134000000 * exp(-(319.542582785997) * invT);
  qf[34] *= k_f;
  qr[34] *= k_f * exp(-(g_RT[1] - 2.000000 * g_RT[4] + g_RT[6]));
  // reaction 35:  CH4 + H <=> CH3 + H2
  k_f = 660 * exp((1.62) * logT - (5454.86865732316) * invT);
  qf[35] *= k_f;
  qr[35] *= k_f * exp(-(-g_RT[0] + g_RT[1] - g_RT[7] + g_RT[8]));
  // reaction 36:  H + HCO <=> CO + H2
  k_f = 73400000;
  qf[36] *= k_f;
  qr[36] *= k_f * exp(-(-g_RT[0] + g_RT[1] - g_RT[9] + g_RT_qss[2]));
  // reaction 37:  CH2O + H <=> H2 + HCO
  k_f = 23000 * exp((1.05) * logT - (1648.03458051046) * invT);
  qf[37] *= k_f;
  qr[37] *= k_f * exp(-(-g_RT[0] + g_RT[1] + g_RT[11] - g_RT_qss[2]));
  // reaction 38:  CH2OH + H <=> CH2O + H2
  k_f = 20000000;
  qf[38] *= k_f;
  qr[38] *= k_f * exp(-(-g_RT[0] + g_RT[1] - g_RT[11] + g_RT_qss[3]));
  // reaction 39:  CH2OH + H <=> CH3 + OH
  k_f = 12000000;
  qf[39] *= k_f;
  qr[39] *= k_f * exp(-(g_RT[1] - g_RT[4] - g_RT[7] + g_RT_qss[3]));
  // reaction 40:  CH2OH + H <=> CH2(S) + H2O
  k_f = 6000000;
  qf[40] *= k_f;
  qr[40] *= k_f * exp(-(g_RT[1] - g_RT[5] - g_RT_qss[1] + g_RT_qss[3]));
  // reaction 41:  H2 + OH <=> H + H2O
  k_f = 216 * exp((1.51) * logT - (1726.03316371019) * invT);
  qf[41] *= k_f;
  qr[41] *= k_f * exp(-(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5]));
  // reaction 42:  2 OH <=> H2O + O
  k_f = 0.0357 * exp((2.4) * logT - (-1061.78716484796) * invT);
  qf[42] *= k_f;
  qr[42] *= k_f * exp(-(-g_RT[2] + 2.000000 * g_RT[4] - g_RT[5]));
  // reaction 43:  HO2 + OH <=> H2O + O2
  k_f = 29000000 * exp(-(-251.60833290236) * invT);
  qf[43] *= k_f;
  qr[43] *= k_f * exp(-(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6]));
  // reaction 44:  CH2 + OH <=> CH2O + H
  k_f = 20000000;
  qf[44] *= k_f;
  qr[44] *= k_f * exp(-(-g_RT[1] + g_RT[4] - g_RT[11] + g_RT_qss[0]));
  // reaction 45:  CH2(S) + OH <=> CH2O + H
  k_f = 30000000;
  qf[45] *= k_f;
  qr[45] *= k_f * exp(-(-g_RT[1] + g_RT[4] - g_RT[11] + g_RT_qss[1]));
  // reaction 46:  CH3 + OH <=> CH2 + H2O
  k_f = 56 * exp((1.6) * logT - (2727.43432866158) * invT);
  qf[46] *= k_f;
  qr[46] *= k_f * exp(-(g_RT[4] - g_RT[5] + g_RT[7] - g_RT_qss[0]));
  // reaction 47:  CH3 + OH <=> CH2(S) + H2O
  k_f = 25010000;
  qf[47] *= k_f;
  qr[47] *= k_f * exp(-(g_RT[4] - g_RT[5] + g_RT[7] - g_RT_qss[1]));
  // reaction 48:  CH4 + OH <=> CH3 + H2O
  k_f = 100 * exp((1.6) * logT - (1570.03599731073) * invT);
  qf[48] *= k_f;
  qr[48] *= k_f * exp(-(g_RT[4] - g_RT[5] - g_RT[7] + g_RT[8]));
  // reaction 49:  CO + OH <=> CO2 + H
  k_f = 47.6 * exp((1.228) * logT - (35.2251666063304) * invT);
  qf[49] *= k_f;
  qr[49] *= k_f * exp(-(-g_RT[1] + g_RT[4] + g_RT[9] - g_RT[10]));
  // reaction 50:  HCO + OH <=> CO + H2O
  k_f = 50000000;
  qf[50] *= k_f;
  qr[50] *= k_f * exp(-(g_RT[4] - g_RT[5] - g_RT[9] + g_RT_qss[2]));
  // reaction 51:  CH2O + OH <=> H2O + HCO
  k_f = 3430 * exp((1.18) * logT - (-224.93784961471) * invT);
  qf[51] *= k_f;
  qr[51] *= k_f * exp(-(g_RT[4] - g_RT[5] + g_RT[11] - g_RT_qss[2]));
  // reaction 52:  CH2OH + OH <=> CH2O + H2O
  k_f = 5000000;
  qf[52] *= k_f;
  qr[52] *= k_f * exp(-(g_RT[4] - g_RT[5] - g_RT[11] + g_RT_qss[3]));
  // reaction 53:  CH2 + HO2 <=> CH2O + OH
  k_f = 20000000;
  qf[53] *= k_f;
  qr[53] *= k_f * exp(-(-g_RT[4] + g_RT[6] - g_RT[11] + g_RT_qss[0]));
  // reaction 54:  CH3 + HO2 <=> CH4 + O2
  k_f = 1000000;
  qf[54] *= k_f;
  qr[54] *= k_f * exp(-(-g_RT[3] + g_RT[6] + g_RT[7] - g_RT[8]));
  // reaction 55:  CO + HO2 <=> CO2 + OH
  k_f = 150000000 * exp(-(11875.9133129914) * invT);
  qf[55] *= k_f;
  qr[55] *= k_f * exp(-(-g_RT[4] + g_RT[6] + g_RT[9] - g_RT[10]));
  // reaction 56:  CH2 + O2 <=> HCO + OH
  k_f = 13200000 * exp(-(754.82499870708) * invT);
  qf[56] *= k_f;
  qr[56] *= k_f * exp(-(g_RT[3] - g_RT[4] + g_RT_qss[0] - g_RT_qss[2]));
  // reaction 57:  CH2 + H2 <=> CH3 + H
  k_f = 0.5 * exp((2) * logT - (3638.25649376812) * invT);
  qf[57] *= k_f;
  qr[57] *= k_f * exp(-(g_RT[0] - g_RT[1] - g_RT[7] + g_RT_qss[0]));
  // reaction 58:  CH2 + CH4 <=> 2 CH3
  k_f = 2.46 * exp((2) * logT - (4161.60182620503) * invT);
  qf[58] *= k_f;
  qr[58] *= k_f * exp(-(-2.000000 * g_RT[7] + g_RT[8] + g_RT_qss[0]));
  // reaction 59:  CH2(S) + N2 <=> CH2 + N2
  k_f = 15000000 * exp(-(301.929999482832) * invT);
  qf[59] *= k_f;
  qr[59] *= k_f * exp(-(g_RT[12] - g_RT[12] - g_RT_qss[0] + g_RT_qss[1]));
  // reaction 60:  CH2(S) + O2 <=> CO + H + OH
  k_f = 28000000;
  qf[60] *= k_f;
  qr[60] *= k_f * exp(-(-g_RT[1] + g_RT[3] - g_RT[4] - g_RT[9] + g_RT_qss[1])) *
            (refCinv);
  // reaction 61:  CH2(S) + O2 <=> CO + H2O
  k_f = 12000000;
  qf[61] *= k_f;
  qr[61] *= k_f * exp(-(g_RT[3] - g_RT[5] - g_RT[9] + g_RT_qss[1]));
  // reaction 62:  CH2(S) + H2 <=> CH3 + H
  k_f = 70000000;
  qf[62] *= k_f;
  qr[62] *= k_f * exp(-(g_RT[0] - g_RT[1] - g_RT[7] + g_RT_qss[1]));
  // reaction 63:  CH2(S) + H2O <=> CH2 + H2O
  k_f = 30000000;
  qf[63] *= k_f;
  qr[63] *= k_f * exp(-(g_RT[5] - g_RT[5] - g_RT_qss[0] + g_RT_qss[1]));
  // reaction 64:  CH2(S) + CH4 <=> 2 CH3
  k_f = 16000000 * exp(-(-286.83349950869) * invT);
  qf[64] *= k_f;
  qr[64] *= k_f * exp(-(-2.000000 * g_RT[7] + g_RT[8] + g_RT_qss[1]));
  // reaction 65:  CH2(S) + CO <=> CH2 + CO
  k_f = 9000000;
  qf[65] *= k_f;
  qr[65] *= k_f * exp(-(g_RT[9] - g_RT[9] - g_RT_qss[0] + g_RT_qss[1]));
  // reaction 66:  CH2(S) + CO2 <=> CH2 + CO2
  k_f = 7000000;
  qf[66] *= k_f;
  qr[66] *= k_f * exp(-(g_RT[10] - g_RT[10] - g_RT_qss[0] + g_RT_qss[1]));
  // reaction 67:  CH2(S) + CO2 <=> CH2O + CO
  k_f = 14000000;
  qf[67] *= k_f;
  qr[67] *= k_f * exp(-(-g_RT[9] + g_RT[10] - g_RT[11] + g_RT_qss[1]));
  // reaction 68:  CH3 + O2 <=> CH2O + OH
  k_f = 36000 * exp(-(4498.75699229419) * invT);
  qf[68] *= k_f;
  qr[68] *= k_f * exp(-(g_RT[3] - g_RT[4] + g_RT[7] - g_RT[11]));
  // reaction 69:  CH3 + HCO <=> CH4 + CO
  k_f = 26480000;
  qf[69] *= k_f;
  qr[69] *= k_f * exp(-(g_RT[7] - g_RT[8] - g_RT[9] + g_RT_qss[2]));
  // reaction 70:  CH2O + CH3 <=> CH4 + HCO
  k_f = 0.00332 * exp((2.81) * logT - (2948.84966161566) * invT);
  qf[70] *= k_f;
  qr[70] *= k_f * exp(-(g_RT[7] - g_RT[8] + g_RT[11] - g_RT_qss[2]));
  // reaction 71:  HCO + O2 <=> CO + HO2
  k_f = 7600000 * exp(-(201.286666321888) * invT);
  qf[71] *= k_f;
  qr[71] *= k_f * exp(-(g_RT[3] - g_RT[6] - g_RT[9] + g_RT_qss[2]));
  // reaction 72:  CH2OH + O2 <=> CH2O + HO2
  k_f = 18000000 * exp(-(452.894999224248) * invT);
  qf[72] *= k_f;
  qr[72] *= k_f * exp(-(g_RT[3] - g_RT[6] - g_RT[11] + g_RT_qss[3]));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
productionRate(amrex::Real* wdot, const amrex::Real* sc, const amrex::Real T)
{
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  const amrex::Real refC = 101325 / 8.31446 * invT;
  const amrex::Real refCinv = 1 / refC;

  for (int i = 0; i < 13; ++i) {
    wdot[i] = 0.0;
  }

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 13; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[13];
  gibbs(g_RT, T);
  amrex::Real g_RT_qss[4];
  gibbs_qss(g_RT_qss, T);

  amrex::Real sc_qss[4];
  amrex::Real kf_qss[41], qf_qss[41], qr_qss[41];
  // Fill sc_qss here
  comp_k_f_qss(T, invT, logT, kf_qss);
  comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, T, g_RT, g_RT_qss);
  comp_sc_qss(sc_qss, qf_qss, qr_qss);

  {
    // reaction 0:  CH2 + H (+M) <=> CH3 (+M)
    const amrex::Real k_f = 25000000000 * exp((-0.8) * logT);
    amrex::Real Corr =
      mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real redP =
      Corr / k_f * 3.2e+15 * exp(-3.14 * logT - 618.956498939805 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.32 * exp(-T * 0.0128205128205128) +
      0.68 * exp(-T * 0.00050125313283208) + exp(-5590 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[1] * sc_qss[0]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] - g_RT[7] + g_RT_qss[0])) * (refC) * (sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[7] += qdot;
  }

  {
    // reaction 1:  CH3 + H (+M) <=> CH4 (+M)
    const amrex::Real k_f =
      12700000000 * exp((-0.63) * logT - (192.731983003208) * invT);
    amrex::Real Corr =
      mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real redP =
      Corr / k_f * 2.477e+21 * exp(-4.76 * logT - 1227.84866456352 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.217 * exp(-T * 0.0135135135135135) +
      0.783 * exp(-T * 0.000340020401224073) + exp(-6964 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[1] * sc[7]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] + g_RT[7] - g_RT[8])) * (refC) * (sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;
  }

  {
    // reaction 2:  H + HCO (+M) <=> CH2O (+M)
    const amrex::Real k_f =
      1090000 * exp((0.48) * logT - (-130.836333109227) * invT);
    amrex::Real Corr =
      mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real redP =
      Corr / k_f * 1350000000000 * exp(-2.57 * logT - 717.083748771726 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.2176 * exp(-T * 0.003690036900369) +
      0.7824 * exp(-T * 0.000362976406533575) + exp(-6570 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[1] * sc_qss[2]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] - g_RT[11] + g_RT_qss[2])) * (refC) * (sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 3:  CH2O + H (+M) <=> CH2OH (+M)
    const amrex::Real k_f =
      540000 * exp((0.454) * logT - (1811.57999689699) * invT);
    amrex::Real Corr =
      mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real redP =
      Corr / k_f * 1.27e+20 * exp(-4.82 * logT - 3286.00482770482 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.2813 * exp(-T * 0.00970873786407767) +
      0.7187 * exp(-T * 0.000774593338497289) + exp(-4160 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[1] * sc[11]);
    const amrex::Real qr = Corr * k_f *
                           exp(-(g_RT[1] + g_RT[11] - g_RT_qss[3])) * (refC) *
                           (sc_qss[3]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[11] -= qdot;
  }

  {
    // reaction 4:  CO + H2 (+M) <=> CH2O (+M)
    const amrex::Real k_f = 43 * exp((1.5) * logT - (40056.0465980557) * invT);
    amrex::Real Corr =
      mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real redP =
      Corr / k_f * 5.07e+15 * exp(-3.42 * logT - 42446.3257606281 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.0679999999999999 * exp(-T * 0.0050761421319797) +
      0.932 * exp(-T * 0.000649350649350649) + exp(-10300 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[0] * sc[9]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[0] + g_RT[9] - g_RT[11])) * (refC) * (sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[9] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 5:  2 O + M <=> O2 + M
    const amrex::Real k_f = 120000 * exp((-1) * logT);
    const amrex::Real Corr = mixture + (1.4) * sc[0] + (14.4) * sc[5] + sc[8] +
                             (0.75) * sc[9] + (2.6) * sc[10];
    const amrex::Real qf = Corr * k_f * ((sc[2] * sc[2]));
    const amrex::Real qr =
      Corr * k_f * exp(-(2.000000 * g_RT[2] - g_RT[3])) * (refC) * (sc[3]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= 2.000000 * qdot;
    wdot[3] += qdot;
  }

  {
    // reaction 6:  H + O + M <=> OH + M
    const amrex::Real k_f = 500000 * exp((-1) * logT);
    const amrex::Real Corr =
      mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real qf = Corr * k_f * (sc[1] * sc[2]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] + g_RT[2] - g_RT[4])) * (refC) * (sc[4]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;
  }

  {
    // reaction 7:  CO + O + M <=> CO2 + M
    const amrex::Real k_f = 602 * exp(-(1509.64999741416) * invT);
    const amrex::Real Corr = mixture + sc[0] + (5) * sc[3] + (5) * sc[5] +
                             sc[8] + (0.5) * sc[9] + (2.5) * sc[10];
    const amrex::Real qf = Corr * k_f * (sc[2] * sc[9]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[2] + g_RT[9] - g_RT[10])) * (refC) * (sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 8:  H + O2 + M <=> HO2 + M
    const amrex::Real k_f = 2800000 * exp((-0.86) * logT);
    const amrex::Real Corr = mixture + (-1) * sc[3] + (-1) * sc[5] +
                             (-0.25) * sc[9] + (0.5) * sc[10] + (-1) * sc[12];
    const amrex::Real qf = Corr * k_f * (sc[1] * sc[3]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC) * (sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
  }

  {
    // reaction 9:  H + O2 + O2 <=> HO2 + O2
    const amrex::Real k_f = 300000000 * exp((-1.72) * logT);
    const amrex::Real qf = k_f * (sc[1] * (sc[3] * sc[3]));
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC) * (sc[3] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
  }

  {
    // reaction 10:  H + O2 + H2O <=> HO2 + H2O
    const amrex::Real k_f = 9380000 * exp((-0.76) * logT);
    const amrex::Real qf = k_f * (sc[1] * sc[3] * sc[5]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC) * (sc[5] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
  }

  {
    // reaction 11:  H + O2 + N2 <=> HO2 + N2
    const amrex::Real k_f = 375000000 * exp((-1.72) * logT);
    const amrex::Real qf = k_f * (sc[1] * sc[3] * sc[12]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC) * (sc[6] * sc[12]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
  }

  {
    // reaction 12:  2 H + M <=> H2 + M
    const amrex::Real k_f = 1000000 * exp((-1) * logT);
    const amrex::Real Corr =
      mixture + (-1) * sc[0] + (-1) * sc[5] + sc[8] + (-1) * sc[10];
    const amrex::Real qf = Corr * k_f * ((sc[1] * sc[1]));
    const amrex::Real qr =
      Corr * k_f * exp(-(-g_RT[0] + 2.000000 * g_RT[1])) * (refC) * (sc[0]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
  }

  {
    // reaction 13:  2 H + H2 <=> H2 + H2
    const amrex::Real k_f = 90000 * exp((-0.6) * logT);
    const amrex::Real qf = k_f * (sc[0] * (sc[1] * sc[1]));
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] + 2.000000 * g_RT[1])) * (refC) * ((sc[0] * sc[0]));
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
  }

  {
    // reaction 14:  2 H + H2O <=> H2 + H2O
    const amrex::Real k_f = 60000000 * exp((-1.25) * logT);
    const amrex::Real qf = k_f * ((sc[1] * sc[1]) * sc[5]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] + 2.000000 * g_RT[1])) * (refC) * (sc[0] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
  }

  {
    // reaction 15:  2 H + CO2 <=> H2 + CO2
    const amrex::Real k_f = 550000000 * exp((-2) * logT);
    const amrex::Real qf = k_f * ((sc[1] * sc[1]) * sc[10]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] + 2.000000 * g_RT[1])) * (refC) * (sc[0] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
  }

  {
    // reaction 16:  H + OH + M <=> H2O + M
    const amrex::Real k_f = 22000000000 * exp((-2) * logT);
    const amrex::Real Corr = mixture + (-0.27) * sc[0] + (2.65) * sc[5] + sc[8];
    const amrex::Real qf = Corr * k_f * (sc[1] * sc[4]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] + g_RT[4] - g_RT[5])) * (refC) * (sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 17:  HCO + H2O <=> CO + H + H2O
    const amrex::Real k_f =
      2244000000000 * exp((-1) * logT - (8554.68331868023) * invT);
    const amrex::Real qf = k_f * (sc[5] * sc_qss[2]);
    const amrex::Real qr = k_f * exp(-(-g_RT[1] - g_RT[9] + g_RT_qss[2])) *
                           (refCinv) * (sc[1] * sc[5] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 18:  HCO + M <=> CO + H + M
    const amrex::Real k_f =
      187000000000 * exp((-1) * logT - (8554.68331868023) * invT);
    const amrex::Real Corr =
      mixture + sc[0] + (-1) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real qf = Corr * k_f * (sc_qss[2]);
    const amrex::Real qr = Corr * k_f *
                           exp(-(-g_RT[1] - g_RT[9] + g_RT_qss[2])) *
                           (refCinv) * (sc[1] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 19:  H2 + O <=> H + OH
    const amrex::Real k_f =
      0.05 * exp((2.67) * logT - (3165.23282791169) * invT);
    const amrex::Real qf = k_f * (sc[0] * sc[2]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4])) * (sc[1] * sc[4]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;
  }

  {
    // reaction 20:  HO2 + O <=> O2 + OH
    const amrex::Real k_f = 20000000;
    const amrex::Real qf = k_f * (sc[2] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6])) * (sc[3] * sc[4]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 21:  CH2 + O <=> H + HCO
    const amrex::Real k_f = 80000000;
    const amrex::Real qf = k_f * (sc[2] * sc_qss[0]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT_qss[0] - g_RT_qss[2])) *
      (sc[1] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
  }

  {
    // reaction 22:  CH2(S) + O <=> CO + H2
    const amrex::Real k_f = 15000000;
    const amrex::Real qf = k_f * (sc[2] * sc_qss[1]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[0] + g_RT[2] - g_RT[9] + g_RT_qss[1])) *
                           (sc[0] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 23:  CH3 + O <=> CH2O + H
    const amrex::Real k_f = 84300000;
    const amrex::Real qf = k_f * (sc[2] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[11])) * (sc[1] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[7] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 24:  CH4 + O <=> CH3 + OH
    const amrex::Real k_f =
      1020 * exp((1.5) * logT - (4327.66332592059) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - g_RT[4] - g_RT[7] + g_RT[8])) * (sc[4] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 25:  HCO + O <=> CO + OH
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[2] * sc_qss[2]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - g_RT[4] - g_RT[9] + g_RT_qss[2])) * (sc[4] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 26:  HCO + O <=> CO2 + H
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[2] * sc_qss[2]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[1] + g_RT[2] - g_RT[10] + g_RT_qss[2])) *
                           (sc[1] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 27:  CH2O + O <=> HCO + OH
    const amrex::Real k_f = 39000000 * exp(-(1781.38699694871) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[11]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[2] - g_RT[4] + g_RT[11] - g_RT_qss[2])) *
                           (sc[4] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[11] -= qdot;
  }

  {
    // reaction 28:  CH2OH + O <=> CH2O + OH
    const amrex::Real k_f = 10000000;
    const amrex::Real qf = k_f * (sc[2] * sc_qss[3]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[2] - g_RT[4] - g_RT[11] + g_RT_qss[3])) *
                           (sc[4] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 29:  CO + O2 <=> CO2 + O
    const amrex::Real k_f = 2500000 * exp(-(24053.7566254656) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[9]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[3] + g_RT[9] - g_RT[10])) * (sc[2] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 30:  CH2O + O2 <=> HCO + HO2
    const amrex::Real k_f = 100000000 * exp(-(20128.6666321888) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[11]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[3] - g_RT[6] + g_RT[11] - g_RT_qss[2])) *
                           (sc[6] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[11] -= qdot;
  }

  {
    // reaction 31:  H + O2 <=> O + OH
    const amrex::Real k_f = 83000000 * exp(-(7252.86180424342) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[3]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4])) * (sc[2] * sc[4]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
  }

  {
    // reaction 32:  H + HO2 <=> H2O + O
    const amrex::Real k_f = 3970000 * exp(-(337.658382754967) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6])) * (sc[2] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 33:  H + HO2 <=> H2 + O2
    const amrex::Real k_f = 28000000 * exp(-(537.435399079441) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6])) * (sc[0] * sc[3]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 34:  H + HO2 <=> 2 OH
    const amrex::Real k_f = 134000000 * exp(-(319.542582785997) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - 2.000000 * g_RT[4] + g_RT[6])) * ((sc[4] * sc[4]));
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[4] += 2.000000 * qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 35:  CH4 + H <=> CH3 + H2
    const amrex::Real k_f =
      660 * exp((1.62) * logT - (5454.86865732316) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] + g_RT[1] - g_RT[7] + g_RT[8])) * (sc[0] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 36:  H + HCO <=> CO + H2
    const amrex::Real k_f = 73400000;
    const amrex::Real qf = k_f * (sc[1] * sc_qss[2]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[0] + g_RT[1] - g_RT[9] + g_RT_qss[2])) *
                           (sc[0] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 37:  CH2O + H <=> H2 + HCO
    const amrex::Real k_f =
      23000 * exp((1.05) * logT - (1648.03458051046) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[11]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[0] + g_RT[1] + g_RT[11] - g_RT_qss[2])) *
                           (sc[0] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[11] -= qdot;
  }

  {
    // reaction 38:  CH2OH + H <=> CH2O + H2
    const amrex::Real k_f = 20000000;
    const amrex::Real qf = k_f * (sc[1] * sc_qss[3]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[0] + g_RT[1] - g_RT[11] + g_RT_qss[3])) *
                           (sc[0] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 39:  CH2OH + H <=> CH3 + OH
    const amrex::Real k_f = 12000000;
    const amrex::Real qf = k_f * (sc[1] * sc_qss[3]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[4] - g_RT[7] + g_RT_qss[3])) * (sc[4] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[7] += qdot;
  }

  {
    // reaction 40:  CH2OH + H <=> CH2(S) + H2O
    const amrex::Real k_f = 6000000;
    const amrex::Real qf = k_f * (sc[1] * sc_qss[3]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[5] - g_RT_qss[1] + g_RT_qss[3])) *
      (sc[5] * sc_qss[1]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 41:  H2 + OH <=> H + H2O
    const amrex::Real k_f =
      216 * exp((1.51) * logT - (1726.03316371019) * invT);
    const amrex::Real qf = k_f * (sc[0] * sc[4]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5])) * (sc[1] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 42:  2 OH <=> H2O + O
    const amrex::Real k_f =
      0.0357 * exp((2.4) * logT - (-1061.78716484796) * invT);
    const amrex::Real qf = k_f * ((sc[4] * sc[4]));
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + 2.000000 * g_RT[4] - g_RT[5])) * (sc[2] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[4] -= 2.000000 * qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 43:  HO2 + OH <=> H2O + O2
    const amrex::Real k_f = 29000000 * exp(-(-251.60833290236) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6])) * (sc[3] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 44:  CH2 + OH <=> CH2O + H
    const amrex::Real k_f = 20000000;
    const amrex::Real qf = k_f * (sc[4] * sc_qss[0]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[1] + g_RT[4] - g_RT[11] + g_RT_qss[0])) *
                           (sc[1] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 45:  CH2(S) + OH <=> CH2O + H
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[4] * sc_qss[1]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[1] + g_RT[4] - g_RT[11] + g_RT_qss[1])) *
                           (sc[1] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 46:  CH3 + OH <=> CH2 + H2O
    const amrex::Real k_f = 56 * exp((1.6) * logT - (2727.43432866158) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[7]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[4] - g_RT[5] + g_RT[7] - g_RT_qss[0])) *
                           (sc[5] * sc_qss[0]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;
  }

  {
    // reaction 47:  CH3 + OH <=> CH2(S) + H2O
    const amrex::Real k_f = 25010000;
    const amrex::Real qf = k_f * (sc[4] * sc[7]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[4] - g_RT[5] + g_RT[7] - g_RT_qss[1])) *
                           (sc[5] * sc_qss[1]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;
  }

  {
    // reaction 48:  CH4 + OH <=> CH3 + H2O
    const amrex::Real k_f = 100 * exp((1.6) * logT - (1570.03599731073) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[4] - g_RT[5] - g_RT[7] + g_RT[8])) * (sc[5] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 49:  CO + OH <=> CO2 + H
    const amrex::Real k_f =
      47.6 * exp((1.228) * logT - (35.2251666063304) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[9]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[4] + g_RT[9] - g_RT[10])) * (sc[1] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 50:  HCO + OH <=> CO + H2O
    const amrex::Real k_f = 50000000;
    const amrex::Real qf = k_f * (sc[4] * sc_qss[2]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[4] - g_RT[5] - g_RT[9] + g_RT_qss[2])) * (sc[5] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 51:  CH2O + OH <=> H2O + HCO
    const amrex::Real k_f =
      3430 * exp((1.18) * logT - (-224.93784961471) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[11]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[4] - g_RT[5] + g_RT[11] - g_RT_qss[2])) *
                           (sc[5] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[11] -= qdot;
  }

  {
    // reaction 52:  CH2OH + OH <=> CH2O + H2O
    const amrex::Real k_f = 5000000;
    const amrex::Real qf = k_f * (sc[4] * sc_qss[3]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[4] - g_RT[5] - g_RT[11] + g_RT_qss[3])) *
                           (sc[5] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 53:  CH2 + HO2 <=> CH2O + OH
    const amrex::Real k_f = 20000000;
    const amrex::Real qf = k_f * (sc[6] * sc_qss[0]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[4] + g_RT[6] - g_RT[11] + g_RT_qss[0])) *
                           (sc[4] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 54:  CH3 + HO2 <=> CH4 + O2
    const amrex::Real k_f = 1000000;
    const amrex::Real qf = k_f * (sc[6] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[3] + g_RT[6] + g_RT[7] - g_RT[8])) * (sc[3] * sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;
  }

  {
    // reaction 55:  CO + HO2 <=> CO2 + OH
    const amrex::Real k_f = 150000000 * exp(-(11875.9133129914) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[9]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[4] + g_RT[6] + g_RT[9] - g_RT[10])) * (sc[4] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 56:  CH2 + O2 <=> HCO + OH
    const amrex::Real k_f = 13200000 * exp(-(754.82499870708) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc_qss[0]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[4] + g_RT_qss[0] - g_RT_qss[2])) *
      (sc[4] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[4] += qdot;
  }

  {
    // reaction 57:  CH2 + H2 <=> CH3 + H
    const amrex::Real k_f = 0.5 * exp((2) * logT - (3638.25649376812) * invT);
    const amrex::Real qf = k_f * (sc[0] * sc_qss[0]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] - g_RT[1] - g_RT[7] + g_RT_qss[0])) * (sc[1] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[7] += qdot;
  }

  {
    // reaction 58:  CH2 + CH4 <=> 2 CH3
    const amrex::Real k_f = 2.46 * exp((2) * logT - (4161.60182620503) * invT);
    const amrex::Real qf = k_f * (sc[8] * sc_qss[0]);
    const amrex::Real qr = k_f *
                           exp(-(-2.000000 * g_RT[7] + g_RT[8] + g_RT_qss[0])) *
                           ((sc[7] * sc[7]));
    const amrex::Real qdot = qf - qr;
    wdot[7] += 2.000000 * qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 59:  CH2(S) + N2 <=> CH2 + N2
    const amrex::Real k_f = 15000000 * exp(-(301.929999482832) * invT);
    const amrex::Real qf = k_f * (sc[12] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[12] - g_RT[12] - g_RT_qss[0] + g_RT_qss[1])) *
      (sc[12] * sc_qss[0]);
    const amrex::Real qdot = qf - qr;
    wdot[12] -= qdot;
    wdot[12] += qdot;
  }

  {
    // reaction 60:  CH2(S) + O2 <=> CO + H + OH
    const amrex::Real k_f = 28000000;
    const amrex::Real qf = k_f * (sc[3] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[3] - g_RT[4] - g_RT[9] + g_RT_qss[1])) *
      (refCinv) * (sc[1] * sc[4] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 61:  CH2(S) + O2 <=> CO + H2O
    const amrex::Real k_f = 12000000;
    const amrex::Real qf = k_f * (sc[3] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[5] - g_RT[9] + g_RT_qss[1])) * (sc[5] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 62:  CH2(S) + H2 <=> CH3 + H
    const amrex::Real k_f = 70000000;
    const amrex::Real qf = k_f * (sc[0] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] - g_RT[1] - g_RT[7] + g_RT_qss[1])) * (sc[1] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[7] += qdot;
  }

  {
    // reaction 63:  CH2(S) + H2O <=> CH2 + H2O
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[5] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[5] - g_RT[5] - g_RT_qss[0] + g_RT_qss[1])) *
      (sc[5] * sc_qss[0]);
    const amrex::Real qdot = qf - qr;
    wdot[5] -= qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 64:  CH2(S) + CH4 <=> 2 CH3
    const amrex::Real k_f = 16000000 * exp(-(-286.83349950869) * invT);
    const amrex::Real qf = k_f * (sc[8] * sc_qss[1]);
    const amrex::Real qr = k_f *
                           exp(-(-2.000000 * g_RT[7] + g_RT[8] + g_RT_qss[1])) *
                           ((sc[7] * sc[7]));
    const amrex::Real qdot = qf - qr;
    wdot[7] += 2.000000 * qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 65:  CH2(S) + CO <=> CH2 + CO
    const amrex::Real k_f = 9000000;
    const amrex::Real qf = k_f * (sc[9] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[9] - g_RT[9] - g_RT_qss[0] + g_RT_qss[1])) *
      (sc[9] * sc_qss[0]);
    const amrex::Real qdot = qf - qr;
    wdot[9] -= qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 66:  CH2(S) + CO2 <=> CH2 + CO2
    const amrex::Real k_f = 7000000;
    const amrex::Real qf = k_f * (sc[10] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[10] - g_RT[10] - g_RT_qss[0] + g_RT_qss[1])) *
      (sc[10] * sc_qss[0]);
    const amrex::Real qdot = qf - qr;
    wdot[10] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 67:  CH2(S) + CO2 <=> CH2O + CO
    const amrex::Real k_f = 14000000;
    const amrex::Real qf = k_f * (sc[10] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[9] + g_RT[10] - g_RT[11] + g_RT_qss[1])) *
      (sc[9] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 68:  CH3 + O2 <=> CH2O + OH
    const amrex::Real k_f = 36000 * exp(-(4498.75699229419) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[4] + g_RT[7] - g_RT[11])) * (sc[4] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[7] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 69:  CH3 + HCO <=> CH4 + CO
    const amrex::Real k_f = 26480000;
    const amrex::Real qf = k_f * (sc[7] * sc_qss[2]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[7] - g_RT[8] - g_RT[9] + g_RT_qss[2])) * (sc[8] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 70:  CH2O + CH3 <=> CH4 + HCO
    const amrex::Real k_f =
      0.00332 * exp((2.81) * logT - (2948.84966161566) * invT);
    const amrex::Real qf = k_f * (sc[7] * sc[11]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[7] - g_RT[8] + g_RT[11] - g_RT_qss[2])) *
                           (sc[8] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;
  }

  {
    // reaction 71:  HCO + O2 <=> CO + HO2
    const amrex::Real k_f = 7600000 * exp(-(201.286666321888) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc_qss[2]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[6] - g_RT[9] + g_RT_qss[2])) * (sc[6] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 72:  CH2OH + O2 <=> CH2O + HO2
    const amrex::Real k_f = 18000000 * exp(-(452.894999224248) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc_qss[3]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[3] - g_RT[6] - g_RT[11] + g_RT_qss[3])) *
                           (sc[6] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[11] += qdot;
  }
}

// GPU version of productionRate: no more use of thermo namespace vectors
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
productionRate_light(
  amrex::Real* wdot,
  amrex::Real* sc,
  amrex::Real* g_RT,
  amrex::Real* g_RT_qss,
  amrex::Real* sc_qss,
  amrex::Real* kf_qss,
  amrex::Real* qf_qss,
  amrex::Real* qr_qss,
  const amrex::Real T,
  const amrex::Real invT,
  const amrex::Real logT)
{
  // reference concentration: P_atm / (RT) in inverse mol/m^3
  const amrex::Real refC = 101325 / 8.31446 * invT;
  const amrex::Real refCinv = 1 / refC;

  for (int i = 0; i < 13; ++i) {
    wdot[i] = 0.0;
  }

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 13; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  gibbs(g_RT, T);
  gibbs_qss(g_RT_qss, T);

  // Fill sc_qss here
  comp_k_f_qss(T, invT, logT, kf_qss);
  comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, T, g_RT, g_RT_qss);
  comp_sc_qss(sc_qss, qf_qss, qr_qss);

  {
    // reaction 0:  CH2 + H (+M) <=> CH3 (+M)
    const amrex::Real k_f = 25000000000 * exp((-0.8) * logT);
    amrex::Real Corr =
      mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real redP =
      Corr / k_f * 3.2e+15 * exp(-3.14 * logT - 618.956498939805 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.32 * exp(-T * 0.0128205128205128) +
      0.68 * exp(-T * 0.00050125313283208) + exp(-5590 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[1] * sc_qss[0]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] - g_RT[7] + g_RT_qss[0])) * (refC) * (sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[7] += qdot;
  }

  {
    // reaction 1:  CH3 + H (+M) <=> CH4 (+M)
    const amrex::Real k_f =
      12700000000 * exp((-0.63) * logT - (192.731983003208) * invT);
    amrex::Real Corr =
      mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real redP =
      Corr / k_f * 2.477e+21 * exp(-4.76 * logT - 1227.84866456352 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.217 * exp(-T * 0.0135135135135135) +
      0.783 * exp(-T * 0.000340020401224073) + exp(-6964 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[1] * sc[7]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] + g_RT[7] - g_RT[8])) * (refC) * (sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;
  }

  {
    // reaction 2:  H + HCO (+M) <=> CH2O (+M)
    const amrex::Real k_f =
      1090000 * exp((0.48) * logT - (-130.836333109227) * invT);
    amrex::Real Corr =
      mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real redP =
      Corr / k_f * 1350000000000 * exp(-2.57 * logT - 717.083748771726 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.2176 * exp(-T * 0.003690036900369) +
      0.7824 * exp(-T * 0.000362976406533575) + exp(-6570 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[1] * sc_qss[2]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] - g_RT[11] + g_RT_qss[2])) * (refC) * (sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 3:  CH2O + H (+M) <=> CH2OH (+M)
    const amrex::Real k_f =
      540000 * exp((0.454) * logT - (1811.57999689699) * invT);
    amrex::Real Corr =
      mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real redP =
      Corr / k_f * 1.27e+20 * exp(-4.82 * logT - 3286.00482770482 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.2813 * exp(-T * 0.00970873786407767) +
      0.7187 * exp(-T * 0.000774593338497289) + exp(-4160 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[1] * sc[11]);
    const amrex::Real qr = Corr * k_f *
                           exp(-(g_RT[1] + g_RT[11] - g_RT_qss[3])) * (refC) *
                           (sc_qss[3]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[11] -= qdot;
  }

  {
    // reaction 4:  CO + H2 (+M) <=> CH2O (+M)
    const amrex::Real k_f = 43 * exp((1.5) * logT - (40056.0465980557) * invT);
    amrex::Real Corr =
      mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real redP =
      Corr / k_f * 5.07e+15 * exp(-3.42 * logT - 42446.3257606281 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.0679999999999999 * exp(-T * 0.0050761421319797) +
      0.932 * exp(-T * 0.000649350649350649) + exp(-10300 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[0] * sc[9]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[0] + g_RT[9] - g_RT[11])) * (refC) * (sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[9] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 5:  2 O + M <=> O2 + M
    const amrex::Real k_f = 120000 * exp((-1) * logT);
    const amrex::Real Corr = mixture + (1.4) * sc[0] + (14.4) * sc[5] + sc[8] +
                             (0.75) * sc[9] + (2.6) * sc[10];
    const amrex::Real qf = Corr * k_f * ((sc[2] * sc[2]));
    const amrex::Real qr =
      Corr * k_f * exp(-(2.000000 * g_RT[2] - g_RT[3])) * (refC) * (sc[3]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= 2.000000 * qdot;
    wdot[3] += qdot;
  }

  {
    // reaction 6:  H + O + M <=> OH + M
    const amrex::Real k_f = 500000 * exp((-1) * logT);
    const amrex::Real Corr =
      mixture + sc[0] + (5) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real qf = Corr * k_f * (sc[1] * sc[2]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] + g_RT[2] - g_RT[4])) * (refC) * (sc[4]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;
  }

  {
    // reaction 7:  CO + O + M <=> CO2 + M
    const amrex::Real k_f = 602 * exp(-(1509.64999741416) * invT);
    const amrex::Real Corr = mixture + sc[0] + (5) * sc[3] + (5) * sc[5] +
                             sc[8] + (0.5) * sc[9] + (2.5) * sc[10];
    const amrex::Real qf = Corr * k_f * (sc[2] * sc[9]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[2] + g_RT[9] - g_RT[10])) * (refC) * (sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 8:  H + O2 + M <=> HO2 + M
    const amrex::Real k_f = 2800000 * exp((-0.86) * logT);
    const amrex::Real Corr = mixture + (-1) * sc[3] + (-1) * sc[5] +
                             (-0.25) * sc[9] + (0.5) * sc[10] + (-1) * sc[12];
    const amrex::Real qf = Corr * k_f * (sc[1] * sc[3]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC) * (sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
  }

  {
    // reaction 9:  H + O2 + O2 <=> HO2 + O2
    const amrex::Real k_f = 300000000 * exp((-1.72) * logT);
    const amrex::Real qf = k_f * (sc[1] * (sc[3] * sc[3]));
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC) * (sc[3] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
  }

  {
    // reaction 10:  H + O2 + H2O <=> HO2 + H2O
    const amrex::Real k_f = 9380000 * exp((-0.76) * logT);
    const amrex::Real qf = k_f * (sc[1] * sc[3] * sc[5]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC) * (sc[5] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
  }

  {
    // reaction 11:  H + O2 + N2 <=> HO2 + N2
    const amrex::Real k_f = 375000000 * exp((-1.72) * logT);
    const amrex::Real qf = k_f * (sc[1] * sc[3] * sc[12]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC) * (sc[6] * sc[12]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
  }

  {
    // reaction 12:  2 H + M <=> H2 + M
    const amrex::Real k_f = 1000000 * exp((-1) * logT);
    const amrex::Real Corr =
      mixture + (-1) * sc[0] + (-1) * sc[5] + sc[8] + (-1) * sc[10];
    const amrex::Real qf = Corr * k_f * ((sc[1] * sc[1]));
    const amrex::Real qr =
      Corr * k_f * exp(-(-g_RT[0] + 2.000000 * g_RT[1])) * (refC) * (sc[0]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
  }

  {
    // reaction 13:  2 H + H2 <=> H2 + H2
    const amrex::Real k_f = 90000 * exp((-0.6) * logT);
    const amrex::Real qf = k_f * (sc[0] * (sc[1] * sc[1]));
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] + 2.000000 * g_RT[1])) * (refC) * ((sc[0] * sc[0]));
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
  }

  {
    // reaction 14:  2 H + H2O <=> H2 + H2O
    const amrex::Real k_f = 60000000 * exp((-1.25) * logT);
    const amrex::Real qf = k_f * ((sc[1] * sc[1]) * sc[5]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] + 2.000000 * g_RT[1])) * (refC) * (sc[0] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
  }

  {
    // reaction 15:  2 H + CO2 <=> H2 + CO2
    const amrex::Real k_f = 550000000 * exp((-2) * logT);
    const amrex::Real qf = k_f * ((sc[1] * sc[1]) * sc[10]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] + 2.000000 * g_RT[1])) * (refC) * (sc[0] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
  }

  {
    // reaction 16:  H + OH + M <=> H2O + M
    const amrex::Real k_f = 22000000000 * exp((-2) * logT);
    const amrex::Real Corr = mixture + (-0.27) * sc[0] + (2.65) * sc[5] + sc[8];
    const amrex::Real qf = Corr * k_f * (sc[1] * sc[4]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[1] + g_RT[4] - g_RT[5])) * (refC) * (sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 17:  HCO + H2O <=> CO + H + H2O
    const amrex::Real k_f =
      2244000000000 * exp((-1) * logT - (8554.68331868023) * invT);
    const amrex::Real qf = k_f * (sc[5] * sc_qss[2]);
    const amrex::Real qr = k_f * exp(-(-g_RT[1] - g_RT[9] + g_RT_qss[2])) *
                           (refCinv) * (sc[1] * sc[5] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 18:  HCO + M <=> CO + H + M
    const amrex::Real k_f =
      187000000000 * exp((-1) * logT - (8554.68331868023) * invT);
    const amrex::Real Corr =
      mixture + sc[0] + (-1) * sc[5] + sc[8] + (0.5) * sc[9] + sc[10];
    const amrex::Real qf = Corr * k_f * (sc_qss[2]);
    const amrex::Real qr = Corr * k_f *
                           exp(-(-g_RT[1] - g_RT[9] + g_RT_qss[2])) *
                           (refCinv) * (sc[1] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 19:  H2 + O <=> H + OH
    const amrex::Real k_f =
      0.05 * exp((2.67) * logT - (3165.23282791169) * invT);
    const amrex::Real qf = k_f * (sc[0] * sc[2]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4])) * (sc[1] * sc[4]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;
  }

  {
    // reaction 20:  HO2 + O <=> O2 + OH
    const amrex::Real k_f = 20000000;
    const amrex::Real qf = k_f * (sc[2] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6])) * (sc[3] * sc[4]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 21:  CH2 + O <=> H + HCO
    const amrex::Real k_f = 80000000;
    const amrex::Real qf = k_f * (sc[2] * sc_qss[0]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT_qss[0] - g_RT_qss[2])) *
      (sc[1] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
  }

  {
    // reaction 22:  CH2(S) + O <=> CO + H2
    const amrex::Real k_f = 15000000;
    const amrex::Real qf = k_f * (sc[2] * sc_qss[1]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[0] + g_RT[2] - g_RT[9] + g_RT_qss[1])) *
                           (sc[0] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 23:  CH3 + O <=> CH2O + H
    const amrex::Real k_f = 84300000;
    const amrex::Real qf = k_f * (sc[2] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[11])) * (sc[1] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[7] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 24:  CH4 + O <=> CH3 + OH
    const amrex::Real k_f =
      1020 * exp((1.5) * logT - (4327.66332592059) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - g_RT[4] - g_RT[7] + g_RT[8])) * (sc[4] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 25:  HCO + O <=> CO + OH
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[2] * sc_qss[2]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - g_RT[4] - g_RT[9] + g_RT_qss[2])) * (sc[4] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 26:  HCO + O <=> CO2 + H
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[2] * sc_qss[2]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[1] + g_RT[2] - g_RT[10] + g_RT_qss[2])) *
                           (sc[1] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 27:  CH2O + O <=> HCO + OH
    const amrex::Real k_f = 39000000 * exp(-(1781.38699694871) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[11]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[2] - g_RT[4] + g_RT[11] - g_RT_qss[2])) *
                           (sc[4] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[11] -= qdot;
  }

  {
    // reaction 28:  CH2OH + O <=> CH2O + OH
    const amrex::Real k_f = 10000000;
    const amrex::Real qf = k_f * (sc[2] * sc_qss[3]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[2] - g_RT[4] - g_RT[11] + g_RT_qss[3])) *
                           (sc[4] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 29:  CO + O2 <=> CO2 + O
    const amrex::Real k_f = 2500000 * exp(-(24053.7566254656) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[9]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[3] + g_RT[9] - g_RT[10])) * (sc[2] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 30:  CH2O + O2 <=> HCO + HO2
    const amrex::Real k_f = 100000000 * exp(-(20128.6666321888) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[11]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[3] - g_RT[6] + g_RT[11] - g_RT_qss[2])) *
                           (sc[6] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[11] -= qdot;
  }

  {
    // reaction 31:  H + O2 <=> O + OH
    const amrex::Real k_f = 83000000 * exp(-(7252.86180424342) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[3]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4])) * (sc[2] * sc[4]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
  }

  {
    // reaction 32:  H + HO2 <=> H2O + O
    const amrex::Real k_f = 3970000 * exp(-(337.658382754967) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6])) * (sc[2] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 33:  H + HO2 <=> H2 + O2
    const amrex::Real k_f = 28000000 * exp(-(537.435399079441) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6])) * (sc[0] * sc[3]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 34:  H + HO2 <=> 2 OH
    const amrex::Real k_f = 134000000 * exp(-(319.542582785997) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - 2.000000 * g_RT[4] + g_RT[6])) * ((sc[4] * sc[4]));
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[4] += 2.000000 * qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 35:  CH4 + H <=> CH3 + H2
    const amrex::Real k_f =
      660 * exp((1.62) * logT - (5454.86865732316) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[0] + g_RT[1] - g_RT[7] + g_RT[8])) * (sc[0] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 36:  H + HCO <=> CO + H2
    const amrex::Real k_f = 73400000;
    const amrex::Real qf = k_f * (sc[1] * sc_qss[2]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[0] + g_RT[1] - g_RT[9] + g_RT_qss[2])) *
                           (sc[0] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 37:  CH2O + H <=> H2 + HCO
    const amrex::Real k_f =
      23000 * exp((1.05) * logT - (1648.03458051046) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[11]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[0] + g_RT[1] + g_RT[11] - g_RT_qss[2])) *
                           (sc[0] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[11] -= qdot;
  }

  {
    // reaction 38:  CH2OH + H <=> CH2O + H2
    const amrex::Real k_f = 20000000;
    const amrex::Real qf = k_f * (sc[1] * sc_qss[3]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[0] + g_RT[1] - g_RT[11] + g_RT_qss[3])) *
                           (sc[0] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 39:  CH2OH + H <=> CH3 + OH
    const amrex::Real k_f = 12000000;
    const amrex::Real qf = k_f * (sc[1] * sc_qss[3]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[4] - g_RT[7] + g_RT_qss[3])) * (sc[4] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[7] += qdot;
  }

  {
    // reaction 40:  CH2OH + H <=> CH2(S) + H2O
    const amrex::Real k_f = 6000000;
    const amrex::Real qf = k_f * (sc[1] * sc_qss[3]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[5] - g_RT_qss[1] + g_RT_qss[3])) *
      (sc[5] * sc_qss[1]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 41:  H2 + OH <=> H + H2O
    const amrex::Real k_f =
      216 * exp((1.51) * logT - (1726.03316371019) * invT);
    const amrex::Real qf = k_f * (sc[0] * sc[4]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5])) * (sc[1] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 42:  2 OH <=> H2O + O
    const amrex::Real k_f =
      0.0357 * exp((2.4) * logT - (-1061.78716484796) * invT);
    const amrex::Real qf = k_f * ((sc[4] * sc[4]));
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + 2.000000 * g_RT[4] - g_RT[5])) * (sc[2] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[4] -= 2.000000 * qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 43:  HO2 + OH <=> H2O + O2
    const amrex::Real k_f = 29000000 * exp(-(-251.60833290236) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[6]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6])) * (sc[3] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 44:  CH2 + OH <=> CH2O + H
    const amrex::Real k_f = 20000000;
    const amrex::Real qf = k_f * (sc[4] * sc_qss[0]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[1] + g_RT[4] - g_RT[11] + g_RT_qss[0])) *
                           (sc[1] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 45:  CH2(S) + OH <=> CH2O + H
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[4] * sc_qss[1]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[1] + g_RT[4] - g_RT[11] + g_RT_qss[1])) *
                           (sc[1] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 46:  CH3 + OH <=> CH2 + H2O
    const amrex::Real k_f = 56 * exp((1.6) * logT - (2727.43432866158) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[7]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[4] - g_RT[5] + g_RT[7] - g_RT_qss[0])) *
                           (sc[5] * sc_qss[0]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;
  }

  {
    // reaction 47:  CH3 + OH <=> CH2(S) + H2O
    const amrex::Real k_f = 25010000;
    const amrex::Real qf = k_f * (sc[4] * sc[7]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[4] - g_RT[5] + g_RT[7] - g_RT_qss[1])) *
                           (sc[5] * sc_qss[1]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;
  }

  {
    // reaction 48:  CH4 + OH <=> CH3 + H2O
    const amrex::Real k_f = 100 * exp((1.6) * logT - (1570.03599731073) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[4] - g_RT[5] - g_RT[7] + g_RT[8])) * (sc[5] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 49:  CO + OH <=> CO2 + H
    const amrex::Real k_f =
      47.6 * exp((1.228) * logT - (35.2251666063304) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[9]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[4] + g_RT[9] - g_RT[10])) * (sc[1] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 50:  HCO + OH <=> CO + H2O
    const amrex::Real k_f = 50000000;
    const amrex::Real qf = k_f * (sc[4] * sc_qss[2]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[4] - g_RT[5] - g_RT[9] + g_RT_qss[2])) * (sc[5] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 51:  CH2O + OH <=> H2O + HCO
    const amrex::Real k_f =
      3430 * exp((1.18) * logT - (-224.93784961471) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[11]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[4] - g_RT[5] + g_RT[11] - g_RT_qss[2])) *
                           (sc[5] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[11] -= qdot;
  }

  {
    // reaction 52:  CH2OH + OH <=> CH2O + H2O
    const amrex::Real k_f = 5000000;
    const amrex::Real qf = k_f * (sc[4] * sc_qss[3]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[4] - g_RT[5] - g_RT[11] + g_RT_qss[3])) *
                           (sc[5] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 53:  CH2 + HO2 <=> CH2O + OH
    const amrex::Real k_f = 20000000;
    const amrex::Real qf = k_f * (sc[6] * sc_qss[0]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[4] + g_RT[6] - g_RT[11] + g_RT_qss[0])) *
                           (sc[4] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 54:  CH3 + HO2 <=> CH4 + O2
    const amrex::Real k_f = 1000000;
    const amrex::Real qf = k_f * (sc[6] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[3] + g_RT[6] + g_RT[7] - g_RT[8])) * (sc[3] * sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;
  }

  {
    // reaction 55:  CO + HO2 <=> CO2 + OH
    const amrex::Real k_f = 150000000 * exp(-(11875.9133129914) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[9]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[4] + g_RT[6] + g_RT[9] - g_RT[10])) * (sc[4] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 56:  CH2 + O2 <=> HCO + OH
    const amrex::Real k_f = 13200000 * exp(-(754.82499870708) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc_qss[0]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[4] + g_RT_qss[0] - g_RT_qss[2])) *
      (sc[4] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[4] += qdot;
  }

  {
    // reaction 57:  CH2 + H2 <=> CH3 + H
    const amrex::Real k_f = 0.5 * exp((2) * logT - (3638.25649376812) * invT);
    const amrex::Real qf = k_f * (sc[0] * sc_qss[0]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] - g_RT[1] - g_RT[7] + g_RT_qss[0])) * (sc[1] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[7] += qdot;
  }

  {
    // reaction 58:  CH2 + CH4 <=> 2 CH3
    const amrex::Real k_f = 2.46 * exp((2) * logT - (4161.60182620503) * invT);
    const amrex::Real qf = k_f * (sc[8] * sc_qss[0]);
    const amrex::Real qr = k_f *
                           exp(-(-2.000000 * g_RT[7] + g_RT[8] + g_RT_qss[0])) *
                           ((sc[7] * sc[7]));
    const amrex::Real qdot = qf - qr;
    wdot[7] += 2.000000 * qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 59:  CH2(S) + N2 <=> CH2 + N2
    const amrex::Real k_f = 15000000 * exp(-(301.929999482832) * invT);
    const amrex::Real qf = k_f * (sc[12] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[12] - g_RT[12] - g_RT_qss[0] + g_RT_qss[1])) *
      (sc[12] * sc_qss[0]);
    const amrex::Real qdot = qf - qr;
    wdot[12] -= qdot;
    wdot[12] += qdot;
  }

  {
    // reaction 60:  CH2(S) + O2 <=> CO + H + OH
    const amrex::Real k_f = 28000000;
    const amrex::Real qf = k_f * (sc[3] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[3] - g_RT[4] - g_RT[9] + g_RT_qss[1])) *
      (refCinv) * (sc[1] * sc[4] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 61:  CH2(S) + O2 <=> CO + H2O
    const amrex::Real k_f = 12000000;
    const amrex::Real qf = k_f * (sc[3] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[5] - g_RT[9] + g_RT_qss[1])) * (sc[5] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 62:  CH2(S) + H2 <=> CH3 + H
    const amrex::Real k_f = 70000000;
    const amrex::Real qf = k_f * (sc[0] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] - g_RT[1] - g_RT[7] + g_RT_qss[1])) * (sc[1] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[7] += qdot;
  }

  {
    // reaction 63:  CH2(S) + H2O <=> CH2 + H2O
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[5] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[5] - g_RT[5] - g_RT_qss[0] + g_RT_qss[1])) *
      (sc[5] * sc_qss[0]);
    const amrex::Real qdot = qf - qr;
    wdot[5] -= qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 64:  CH2(S) + CH4 <=> 2 CH3
    const amrex::Real k_f = 16000000 * exp(-(-286.83349950869) * invT);
    const amrex::Real qf = k_f * (sc[8] * sc_qss[1]);
    const amrex::Real qr = k_f *
                           exp(-(-2.000000 * g_RT[7] + g_RT[8] + g_RT_qss[1])) *
                           ((sc[7] * sc[7]));
    const amrex::Real qdot = qf - qr;
    wdot[7] += 2.000000 * qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 65:  CH2(S) + CO <=> CH2 + CO
    const amrex::Real k_f = 9000000;
    const amrex::Real qf = k_f * (sc[9] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[9] - g_RT[9] - g_RT_qss[0] + g_RT_qss[1])) *
      (sc[9] * sc_qss[0]);
    const amrex::Real qdot = qf - qr;
    wdot[9] -= qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 66:  CH2(S) + CO2 <=> CH2 + CO2
    const amrex::Real k_f = 7000000;
    const amrex::Real qf = k_f * (sc[10] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[10] - g_RT[10] - g_RT_qss[0] + g_RT_qss[1])) *
      (sc[10] * sc_qss[0]);
    const amrex::Real qdot = qf - qr;
    wdot[10] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 67:  CH2(S) + CO2 <=> CH2O + CO
    const amrex::Real k_f = 14000000;
    const amrex::Real qf = k_f * (sc[10] * sc_qss[1]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[9] + g_RT[10] - g_RT[11] + g_RT_qss[1])) *
      (sc[9] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 68:  CH3 + O2 <=> CH2O + OH
    const amrex::Real k_f = 36000 * exp(-(4498.75699229419) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[4] + g_RT[7] - g_RT[11])) * (sc[4] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[7] -= qdot;
    wdot[11] += qdot;
  }

  {
    // reaction 69:  CH3 + HCO <=> CH4 + CO
    const amrex::Real k_f = 26480000;
    const amrex::Real qf = k_f * (sc[7] * sc_qss[2]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[7] - g_RT[8] - g_RT[9] + g_RT_qss[2])) * (sc[8] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 70:  CH2O + CH3 <=> CH4 + HCO
    const amrex::Real k_f =
      0.00332 * exp((2.81) * logT - (2948.84966161566) * invT);
    const amrex::Real qf = k_f * (sc[7] * sc[11]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[7] - g_RT[8] + g_RT[11] - g_RT_qss[2])) *
                           (sc[8] * sc_qss[2]);
    const amrex::Real qdot = qf - qr;
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;
  }

  {
    // reaction 71:  HCO + O2 <=> CO + HO2
    const amrex::Real k_f = 7600000 * exp(-(201.286666321888) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc_qss[2]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[6] - g_RT[9] + g_RT_qss[2])) * (sc[6] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[9] += qdot;
  }

  {
    // reaction 72:  CH2OH + O2 <=> CH2O + HO2
    const amrex::Real k_f = 18000000 * exp(-(452.894999224248) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc_qss[3]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[3] - g_RT[6] - g_RT[11] + g_RT_qss[3])) *
                           (sc[6] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[11] += qdot;
  }

  return;
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[])
{

  // convert to SI
  for (int id = 0; id < 13; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, T);

  // convert to chemkin units
  for (int id = 0; id < 13; ++id) {
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
  amrex::Real c[13]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 13; i++) {
    YOW += y[i] * imw(i);
  }
  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // multiply by 1e6 so c goes to SI
  PWORT *= 1e6;
  // Now compute conversion (and go to SI)
  for (int i = 0; i < 13; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 13; ++id) {
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
  amrex::Real c[13]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 13; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 13; ++id) {
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
  amrex::Real c[13]; // temporary storage

  // See Eq 8 with an extra 1e6 so c goes to SI
  for (int i = 0; i < 13; i++) {
    c[i] = 1e6 * rho * y[i] * imw(i);
  }

  // call productionRate
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 13; ++id) {
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
  amrex::Real c[13];  // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 2.016000;   // H2
  XW += x[1] * 1.008000;   // H
  XW += x[2] * 15.999000;  // O
  XW += x[3] * 31.998000;  // O2
  XW += x[4] * 17.007000;  // OH
  XW += x[5] * 18.015000;  // H2O
  XW += x[6] * 33.006000;  // HO2
  XW += x[7] * 15.035000;  // CH3
  XW += x[8] * 16.043000;  // CH4
  XW += x[9] * 28.010000;  // CO
  XW += x[10] * 44.009000; // CO2
  XW += x[11] * 30.026000; // CH2O
  XW += x[12] * 28.014000; // N2
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 13; ++id) {
    c[id] = x[id] * ROW;
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 13; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int kcharge[])
{
  kcharge[0] = 0;  // H2
  kcharge[1] = 0;  // H
  kcharge[2] = 0;  // O
  kcharge[3] = 0;  // O2
  kcharge[4] = 0;  // OH
  kcharge[5] = 0;  // H2O
  kcharge[6] = 0;  // HO2
  kcharge[7] = 0;  // CH3
  kcharge[8] = 0;  // CH4
  kcharge[9] = 0;  // CO
  kcharge[10] = 0; // CO2
  kcharge[11] = 0; // CH2O
  kcharge[12] = 0; // N2
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real zk[])
{

  int kchrg[13];
  CKCHRG(kchrg);

  for (int id = 0; id < 13; ++id) {
    zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
  }
}

// compute an approx to the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian_precond(
  amrex::Real* J, const amrex::Real* sc, const amrex::Real T, const int HP)
{

  // Analytical Jacobian with QSSA is only supported with symbolic
  // implementation. Re-build in ceptr with -qsj flag.
  amrex::Abort();
}

// compute an approx to the reaction Jacobian (for preconditioning)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT_SIMPLIFIED(
  amrex::Real* J, const amrex::Real* sc, const amrex::Real* Tp, const int* HP)
{
  amrex::Real c[13];

  for (int k = 0; k < 13; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian_precond(J, c, *Tp, *HP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 13; k++) {
    J[182 + k] *= 1.e-6;
    J[k * 14 + 13] *= 1.e6;
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian(amrex::Real* J, amrex::Real* sc, amrex::Real T, const int consP)
{
#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  const amrex::Real refC = 101325 / 8.31446 * invT;
  const amrex::Real refCinv = 1 / refC;
  amrex::Real g_RT[13];
  amrex::Real h_RT[13];
  amrex::Real g_RT_qss[4];
  amrex::Real h_RT_qss[4];
  amrex::Real sc_qss[4];
  amrex::Real kf_qss[41];
  amrex::Real T_pert1, pertT;
  amrex::Real wdot_pert1[13], wdot[13];

  // dwdot/dT by finite difference
  pertT = 1e-2;
  T_pert1 = T + pertT;

  const amrex::Real invT_pert1 = 1.0 / T_pert1;
  const amrex::Real logT_pert1 = log(T_pert1);
  productionRate_light(
    wdot_pert1, sc, g_RT, g_RT_qss, sc_qss, kf_qss, &J[0], &J[41], T_pert1,
    invT_pert1, logT_pert1);
  productionRate_light(
    wdot, sc, g_RT, g_RT_qss, sc_qss, kf_qss, &J[0], &J[41], T, invT, logT);

  // compute the species enthalpy
  speciesEnthalpy(h_RT, T);
  for (int i = 0; i < 196; i++) {
    J[i] = 0.0;
  }
  // Species terms
  const amrex::Real x2 =
    sc[11] + sc[12] + sc[1] + sc[2] + sc[4] + sc[6] + sc[7] + 2 * sc[8];
  const amrex::Real x6 = 2 * sc[0] + 2 * sc[10] + sc[3] + 1.5 * sc[9] + x2;
  const amrex::Real x7 = 6.0 * sc[5] + x6;
  amrex::Real x10 = -g_RT[11] - g_RT[1] + g_RT_qss[3];
  amrex::Real x11 = -3286.0048277048199 * invT - 4.82 * logT;
  amrex::Real x15 = 1.27e+20 * x7 * exp(x11);
  amrex::Real x17 = log(
    exp(-4160.0 * invT) + 0.2813 * exp(-0.0097087378640776708 * T) +
    0.7187 * exp(-0.00077459333849728897 * T));
  amrex::Real x20 = -0.29097730287517898 * x17 +
                    0.43429448190325198 * log(x15 / kf_qss[2]) - 0.4;
  amrex::Real x22 = -0.55155399201712996 * x17 - 0.14 * x20 + 0.75;
  amrex::Real x25 = ((x20) * (x20)) / ((x22) * (x22)) + 1.0;
  amrex::Real x28 =
    exp(M_LN10 * (0.43429448190325198 * x17 / x25)) / (1.0 + x15 / kf_qss[2]);
  amrex::Real x30 = refC * x28 * exp(x10 + x11);
  amrex::Real x33 = -kf_qss[10] * sc[2] - kf_qss[14] * sc[1] -
                    kf_qss[15] * sc[1] - kf_qss[16] * sc[1] -
                    kf_qss[23] * sc[4] - kf_qss[40] * sc[3] -
                    1.27e+20 * x30 * x7 + 10.0e-13;
  const amrex::Real x37 = exp(-g_RT[1] + g_RT[5] + g_RT_qss[1] - g_RT_qss[3]);
  const amrex::Real x41 = (kf_qss[16] * sc[5] * x37 + 10.0e-13) / x33;
  amrex::Real x43 = g_RT[4] - g_RT[5];
  const amrex::Real x45 = exp(-g_RT[7] + g_RT_qss[1] - x43);
  amrex::Real x50 =
    -kf_qss[16] * sc[5] * x37 - kf_qss[18] * sc[4] - kf_qss[20] * sc[5] * x45 -
    kf_qss[28] * sc[12] - kf_qss[29] * sc[3] - kf_qss[30] * sc[3] -
    kf_qss[31] * sc[0] - kf_qss[32] * sc[5] - kf_qss[33] * sc[8] -
    kf_qss[34] * sc[9] - kf_qss[35] * sc[10] - kf_qss[36] * sc[10] -
    kf_qss[6] * sc[2] + 10.0e-13;
  amrex::Real x52 = kf_qss[16] * sc[1] + 10.0e-13;
  amrex::Real x53 = x52 / x50;
  amrex::Real x54 = -x41 * x53 + 1;
  amrex::Real x57 = exp(g_RT_qss[0] - g_RT_qss[1]);
  amrex::Real x66 = kf_qss[28] * sc[12] * x57 + kf_qss[32] * sc[5] * x57 +
                    kf_qss[34] * sc[9] * x57 + kf_qss[35] * sc[10] * x57 +
                    10.0e-13;
  const amrex::Real x68 = x66 / (x50 * x54);
  amrex::Real x70 = -717.08374877172605 * invT - 2.57 * logT;
  amrex::Real x72 = 1350000000000.0 * exp(x70);
  amrex::Real x75 = log(
    exp(-6570.0 * invT) + 0.2176 * exp(-0.0036900369003690001 * T) +
    0.7824 * exp(-0.00036297640653357502 * T));
  amrex::Real x78 = -0.29097730287517898 * x75 +
                    0.43429448190325198 * log(x7 * x72 / kf_qss[1]) - 0.4;
  amrex::Real x81 = -0.55155399201712996 * x75 - 0.14 * x78 + 0.75;
  amrex::Real x84 = ((x78) * (x78)) / ((x81) * (x81)) + 1.0;
  amrex::Real x87 = exp(M_LN10 * (0.43429448190325198 * x75 / x84)) /
                    (1.0 + x7 * x72 / kf_qss[1]);
  amrex::Real x92 = g_RT[2] - g_RT[4];
  amrex::Real x94 = g_RT[11] - g_RT_qss[2];
  amrex::Real x97 = g_RT[7] - g_RT[8];
  amrex::Real x103 =
    kf_qss[25] * exp(-g_RT[3] + g_RT[4] - g_RT_qss[0] + g_RT_qss[2]);
  const amrex::Real x107 = exp(g_RT[1] - g_RT[2] - g_RT_qss[0] + g_RT_qss[2]);
  amrex::Real x110 = g_RT[3] - g_RT[6];
  amrex::Real x117 = -g_RT[0] + g_RT[11] + g_RT[1] - g_RT_qss[2];
  amrex::Real x119 = -kf_qss[11] * sc[6] * exp(-x110 - x94) -
                     kf_qss[12] * sc[1] - kf_qss[13] * sc[0] * exp(-x117) -
                     kf_qss[21] * sc[4] - kf_qss[22] * sc[5] * exp(-x43 - x94) -
                     kf_qss[37] * sc[7] - kf_qss[38] * sc[8] * exp(-x94 - x97) -
                     kf_qss[39] * sc[3] - kf_qss[3] * sc[5] - kf_qss[4] * x6 -
                     kf_qss[5] * sc[1] * x107 - kf_qss[7] * sc[2] -
                     kf_qss[8] * sc[2] - kf_qss[9] * sc[4] * exp(-x92 - x94) -
                     sc[1] * x7 * x72 * x87 - sc[4] * x103 + 10.0e-13;
  const amrex::Real x124 =
    (kf_qss[25] * sc[3] + kf_qss[5] * sc[2] + 10.0e-13) / x119;
  amrex::Real x128 =
    3200000000000000.0 * exp(-618.95649893980499 * invT - 3.14 * logT);
  amrex::Real x131 = log(
    exp(-5590.0 * invT) + 0.32 * exp(-0.012820512820512799 * T) +
    0.68 * exp(-0.00050125313283207998 * T));
  amrex::Real x134 = -0.29097730287517898 * x131 +
                     0.43429448190325198 * log(x128 * x7 / kf_qss[0]) - 0.4;
  amrex::Real x139 = 2.3703703703703702 * ((x134) * (x134)) /
                     ((-0.73540532269 * x131 - 0.18666666666666668 * x134 + 1) *
                      (-0.73540532269 * x131 - 0.18666666666666668 * x134 + 1) *
                      (-0.73540532269 * x131 - 0.18666666666666668 * x134 + 1));
  amrex::Real x140 = 0.43429448190325198 / x7;
  amrex::Real x143 =
    1.7777777777777777 * x134 /
    ((-0.73540532268950665 * x131 - 0.18666666666666668 * x134 + 1) *
     (-0.73540532268950665 * x131 - 0.18666666666666668 * x134 + 1));
  x134 = 1.7777777777777777 * ((x134) * (x134)) /
           ((-0.73540532268950665 * x131 - 0.18666666666666668 * x134 + 1) *
            (-0.73540532268950665 * x131 - 0.18666666666666668 * x134 + 1)) +
         1.0;
  amrex::Real x149 = exp(M_LN10 * (0.43429448190325198 * x131 / x134)) /
                     (1.0 + x128 * x7 / kf_qss[0]);
  amrex::Real x154 = sc[1] * x128 * x131 * x149 * x7 / ((x134) * (x134));
  amrex::Real x155 =
    sc[1] * x149 * exp(-618.95649893980499 * invT - 3.14 * logT);
  amrex::Real x156 = exp(M_LN10 * (0.43429448190325198 * x131 / x134)) * x7 /
                     (kf_qss[0] * ((1.0 + x128 * x7 / kf_qss[0]) *
                                   (1.0 + x128 * x7 / kf_qss[0])));
  amrex::Real x161 =
    2.048e+31 * sc[1] * x156 * exp(-1237.91299787961 * invT - 6.28 * logT) +
    x154 * (0.56 * x139 * x140 + 4.0 * x140 * x143) - 6400000000000000.0 * x155;
  amrex::Real x165 = -kf_qss[17] * sc[4] -
                     kf_qss[19] * sc[5] * exp(-g_RT[7] + g_RT_qss[0] - x43) -
                     kf_qss[24] * sc[6] - kf_qss[25] * sc[3] -
                     kf_qss[26] * sc[0] - kf_qss[27] * sc[8] -
                     kf_qss[28] * sc[12] * x57 - kf_qss[32] * sc[5] * x57 -
                     kf_qss[34] * sc[9] * x57 - kf_qss[35] * sc[10] * x57 -
                     kf_qss[5] * sc[2] - sc[1] * x128 * x149 * x7 + 10.0e-13;
  amrex::Real x167 = kf_qss[5] * sc[1] * x107 + sc[4] * x103 + 10.0e-13;
  amrex::Real x169 = x124 * x167 / ((x165) * (x165));
  amrex::Real x174 = x140 * ((x78) * (x78)) / ((x81) * (x81) * (x81));
  amrex::Real x175 = x78 / ((x81) * (x81));
  amrex::Real x179 = sc[1] * x7 * x72 * x75 * x87 / ((x84) * (x84));
  amrex::Real x180 =
    exp(M_LN10 * (0.43429448190325198 * x75 / x84)) * x7 /
    (kf_qss[1] * ((1.0 + x7 * x72 / kf_qss[1]) * (1.0 + x7 * x72 / kf_qss[1])));
  amrex::Real x183 =
    sc[1] * x180 * exp(-1434.1674975434501 * invT - 5.14 * logT);
  amrex::Real x185 = -2 * kf_qss[4] - 2700000000000.0 * sc[1] * x87 * exp(x70) +
                     x179 * (4.0 * x140 * x175 + 0.56 * x174) +
                     3.645e+24 * x183;
  amrex::Real x186 = -kf_qss[13] * exp(-x117) + x185;
  amrex::Real x188 = x167 / x165;
  amrex::Real x191 = x188 *
                     (kf_qss[25] * sc[3] + kf_qss[5] * sc[2] + 10.0e-13) /
                     ((x119) * (x119));
  amrex::Real x195 = x140 * ((x20) * (x20)) / ((x22) * (x22) * (x22));
  amrex::Real x197 = x140 * x20 / ((x22) * (x22));
  amrex::Real x200 = 1.27e+20 * x17 * x30 * x7 / ((x25) * (x25));
  amrex::Real x202 =
    exp(M_LN10 * (0.43429448190325198 * x17 / x25)) /
    (kf_qss[2] * ((1.0 + x15 / kf_qss[2]) * (1.0 + x15 / kf_qss[2])));
  const amrex::Real x204 = refC * x7;
  amrex::Real x205 = x204 * exp(-6572.0096554096399 * invT - 9.64 * logT + x10);
  amrex::Real x207 = -2.54e+20 * refC * x28 * exp(x10 + x11) +
                     x200 * (0.56 * x195 + 4.0 * x197) +
                     3.2258e+40 * x202 * x205;
  amrex::Real x209 = (kf_qss[16] * sc[5] * x37 + 10.0e-13) / ((x33) * (x33));
  amrex::Real x210 = x209 * x53;
  amrex::Real x214 = x41 * x52 / ((x50) * (x50));
  amrex::Real x217 = kf_qss[28] * sc[12] + kf_qss[32] * sc[5] +
                     kf_qss[34] * sc[9] + kf_qss[35] * sc[10] + 10.0e-13;
  amrex::Real x220 = x217 * x66 / (x165 * x50 * ((x54) * (x54)));
  amrex::Real x222 = x217 * x68 / ((x165) * (x165));
  amrex::Real x223 = x217 / (x165 * x54);
  amrex::Real x225 = x223 * x66 / ((x50) * (x50));
  const amrex::Real x228 = exp(g_RT[11] - g_RT_qss[3] - x43);
  amrex::Real x232 = kf_qss[40] * exp(g_RT[11] - g_RT_qss[3] - x110);
  const amrex::Real x234 = exp(g_RT[0] + g_RT[11] - g_RT[1] - g_RT_qss[3]);
  const amrex::Real x237 = exp(g_RT[11] - g_RT_qss[3] - x92);
  const amrex::Real x245 = exp(-g_RT[1] + g_RT[4] + g_RT[7] - g_RT_qss[3]);
  amrex::Real x248 =
    kf_qss[10] * sc[11] * sc[4] * x237 + kf_qss[14] * sc[0] * sc[11] * x234 +
    kf_qss[15] * sc[4] * sc[7] * x245 + kf_qss[23] * sc[11] * sc[5] * x228 +
    sc[11] * sc[1] * x15 * x28 + sc[11] * sc[6] * x232 + 10.0e-13;
  const amrex::Real x249 = x248 / x33;
  const amrex::Real x253 = exp(-g_RT[3] + g_RT[5] + g_RT[9] - g_RT_qss[1]);
  const amrex::Real x256 = exp(g_RT[0] - g_RT[2] + g_RT[9] - g_RT_qss[1]);
  const amrex::Real x260 = -2.0 * g_RT[7] + g_RT[8];
  const amrex::Real x265 = exp(g_RT[11] + g_RT[1] - g_RT[4] - g_RT_qss[1]);
  const amrex::Real x269 = refCinv * sc[1];
  const amrex::Real x270 = sc[9] * x269;
  const amrex::Real x271 =
    exp(g_RT[1] - g_RT[3] + g_RT[4] + g_RT[9] - g_RT_qss[1]);
  const amrex::Real x274 = exp(-g_RT[10] + g_RT[11] + g_RT[9] - g_RT_qss[1]);
  const amrex::Real x279 = exp(-g_RT[0] + g_RT[1] + g_RT[7] - g_RT_qss[1]);
  amrex::Real x282 =
    kf_qss[18] * sc[11] * sc[1] * x265 + kf_qss[20] * sc[4] * sc[7] +
    kf_qss[29] * sc[4] * x270 * x271 + kf_qss[30] * sc[5] * sc[9] * x253 +
    kf_qss[31] * sc[1] * sc[7] * x279 +
    kf_qss[33] * ((sc[7]) * (sc[7])) * exp(-g_RT_qss[1] - x260) +
    kf_qss[36] * sc[11] * sc[9] * x274 + kf_qss[6] * sc[0] * sc[9] * x256 +
    10.0e-13;
  amrex::Real x283 = x249 * x53 - x282 / x50;
  amrex::Real x286 = -g_RT[9] + g_RT_qss[2];
  const amrex::Real x287 = exp(-x286 - x43);
  const amrex::Real x292 = exp(g_RT[0] - g_RT[1] + g_RT[9] - g_RT_qss[2]);
  const amrex::Real x295 = exp(-x286 - x92);
  amrex::Real x299 = exp(g_RT[1] - x286);
  amrex::Real x304 = kf_qss[39] * exp(-x110 - x286);
  const amrex::Real x308 = exp(g_RT[10] + g_RT[1] - g_RT[2] - g_RT_qss[2]);
  const amrex::Real x311 = exp(-x286 - x97);
  const amrex::Real x314 = refC * sc[11];
  amrex::Real x316 = g_RT[11] - g_RT[1] - g_RT_qss[2];
  amrex::Real x321 =
    kf_qss[11] * sc[11] * sc[3] + kf_qss[12] * sc[0] * sc[9] * x292 +
    kf_qss[13] * sc[11] * sc[1] + kf_qss[21] * sc[5] * sc[9] * x287 +
    kf_qss[22] * sc[11] * sc[4] + kf_qss[37] * sc[8] * sc[9] * x311 +
    kf_qss[38] * sc[11] * sc[7] + kf_qss[3] * sc[5] * x270 * x299 +
    kf_qss[4] * x270 * x299 * x6 + kf_qss[7] * sc[4] * sc[9] * x295 +
    kf_qss[8] * sc[10] * sc[1] * x308 + kf_qss[9] * sc[11] * sc[2] +
    sc[6] * sc[9] * x304 + 1350000000000.0 * x314 * x7 * x87 * exp(x316 + x70) +
    10.0e-13;
  const amrex::Real x322 = x321 / x119;
  const amrex::Real x323 = exp(g_RT[11] + g_RT[1] - g_RT[4] - g_RT_qss[0]);
  amrex::Real x326 = -g_RT[1] + g_RT[7] - g_RT_qss[0];
  amrex::Real x327 =
    refC * exp(-618.95649893980499 * invT - 3.14 * logT + x326);
  const amrex::Real x333 = exp(g_RT[11] + g_RT[4] - g_RT[6] - g_RT_qss[0]);
  amrex::Real x339 =
    kf_qss[26] * exp(-g_RT[0] + g_RT[1] + g_RT[7] - g_RT_qss[0]);
  amrex::Real x342 =
    kf_qss[17] * sc[11] * sc[1] * x323 + kf_qss[19] * sc[4] * sc[7] +
    kf_qss[24] * sc[11] * sc[4] * x333 +
    kf_qss[27] * ((sc[7]) * (sc[7])) * exp(-g_RT_qss[0] - x260) +
    sc[1] * sc[7] * x339 + 3200000000000000.0 * sc[7] * x149 * x327 * x7 +
    10.0e-13;
  amrex::Real x345 = (x188 * x322 - x223 * x283 - x342 / x165) /
                     ((-x124 * x188 + 1 - x217 * x68 / x165) *
                      (-x124 * x188 + 1 - x217 * x68 / x165));
  const amrex::Real x346 =
    x345 *
    (-kf_qss[31] * x225 + x169 * (-kf_qss[26] + x161) + x186 * x191 +
     x220 * (-kf_qss[31] * x214 + x207 * x210) + x222 * (-kf_qss[26] + x161));
  amrex::Real x347 = x342 / ((x165) * (x165));
  x134 =
    3200000000000000.0 * sc[7] * x131 * x149 * x327 * x7 / ((x134) * (x134));
  amrex::Real x351 =
    refC * sc[7] * exp(-1237.91299787961 * invT - 6.28 * logT + x326);
  amrex::Real x352 = -(6400000000000000.0 * sc[7] * x149 * x327 -
                       x134 * (0.56 * x139 * x140 + 4.0 * x140 * x143) -
                       2.048e+31 * x156 * x351) /
                     x165;
  x84 =
    1350000000000.0 * x314 * x7 * x75 * x87 * exp(x316 + x70) / ((x84) * (x84));
  amrex::Real x356 =
    x314 * exp(-1434.1674975434501 * invT - 5.14 * logT + x316);
  amrex::Real x359 = 2 * kf_qss[4] * x270 * x299 - 3.645e+24 * x180 * x356 +
                     2700000000000.0 * x314 * x87 * exp(x316 + x70) -
                     x84 * (4.0 * x140 * x175 + 0.56 * x174);
  amrex::Real x361 = x217 * x283 / (x165 * ((x54) * (x54)));
  x321 = x321 / ((x119) * (x119));
  x167 = x167 * x322 / ((x165) * (x165));
  x282 = x282 / ((x50) * (x50));
  x248 = x248 / ((x33) * (x33));
  x25 = sc[11] * sc[1] * x15 * x17 * x28 / ((x25) * (x25));
  amrex::Real x371 = sc[11] * sc[1] * x28 * exp(x11);
  amrex::Real x373 =
    sc[11] * sc[1] * x7 * exp(-6572.0096554096399 * invT - 9.64 * logT);
  amrex::Real x375 = -3.2258e+40 * x202 * x373 -
                     x25 * (0.56 * x195 + 4.0 * x197) + 2.54e+20 * x371;
  amrex::Real x378 = x249 * x52 / ((x50) * (x50));
  amrex::Real x379 = -kf_qss[31] * x282 + kf_qss[31] * x378 -
                     kf_qss[6] * sc[9] * x256 / x50 - x207 * x248 * x53 +
                     x53 * (kf_qss[14] * sc[11] * x234 + x375) / x33;
  amrex::Real x381 = x217 * x283 / (((x165) * (x165)) * x54);
  amrex::Real x382 = (1.0 / ((-x124 * x188 + 1 - x217 * x68 / x165)));
  const amrex::Real x383 =
    x382 *
    (-x167 * (-kf_qss[26] + x161) - x186 * x188 * x321 - x223 * x379 +
     x347 * (-kf_qss[26] + x161) + x352 +
     x361 * (-kf_qss[31] * x214 + x207 * x210) + x381 * (-kf_qss[26] + x161) +
     x188 * (kf_qss[12] * sc[9] * x292 + x359) / x119);
  amrex::Real x389 = -kf_qss[4] - sc[1] * x72 * x87 +
                     x179 * (2.0 * x140 * x175 + 0.28 * x140 * ((x78) * (x78)) /
                                                   ((x81) * (x81) * (x81))) +
                     1.8225e+24 * x183;
  x72 = -kf_qss[12] - kf_qss[5] * x107 + x389 - x7 * x72 * x87;
  amrex::Real x392 = sc[1] * x156 * exp(-1237.91299787961 * invT - 6.28 * logT);
  amrex::Real x393 = -sc[1] * x128 * x149 +
                     x154 * (0.28 * x139 * x140 + 2.0 * x140 * x143) +
                     1.024e+31 * x392;
  x128 = -x128 * x149 * x7 + x393;
  x20 = 0.28 * x140 * ((x20) * (x20)) / ((x22) * (x22) * (x22)) +
        2.0 * x140 * x20 / ((x22) * (x22));
  x22 = x20 * x200 + 1.6129e+40 * x202 * x205 - 1.27e+20 * x30;
  amrex::Real x397 = -kf_qss[14] - kf_qss[15] - kf_qss[16] + x22;
  amrex::Real x399 = -kf_qss[16] * x41 / x50 + x210 * x397;
  const amrex::Real x400 =
    x345 * (-kf_qss[5] * x107 * x124 / x165 + x128 * x169 + x128 * x222 +
            x191 * x72 + x220 * x399);
  amrex::Real x402 = 3200000000000000.0 * sc[7] * x149 * x327 -
                     x134 * (0.28 * x139 * x140 + 2.0 * x140 * x143) -
                     1.024e+31 * x156 * x351;
  x78 = kf_qss[4] * x270 * x299 - 1.8225e+24 * x180 * x356 +
        1350000000000.0 * x314 * x87 * exp(x316 + x70) -
        x84 * (2.0 * x140 * x175 +
               0.28 * x140 * ((x78) * (x78)) / ((x81) * (x81) * (x81)));
  x81 = (kf_qss[13] * sc[11] + kf_qss[3] * refCinv * sc[5] * sc[9] * x299 +
         kf_qss[4] * refCinv * sc[9] * x299 * x6 + kf_qss[8] * sc[10] * x308 +
         x78) /
        x119;
  x20 = -x20 * x25 - 1.6129e+40 * x202 * x373 + 1.27e+20 * x371;
  const amrex::Real x413 =
    (kf_qss[16] * x249 / x50 - x248 * x397 * x53 -
     (kf_qss[18] * sc[11] * x265 + kf_qss[29] * refCinv * sc[4] * sc[9] * x271 +
      kf_qss[31] * sc[7] * x279) /
       x50 +
     x53 * (sc[11] * x15 * x28 + x20) / x33) /
    x54;
  x128 = x382 * (kf_qss[5] * x107 * x322 / x165 - x128 * x167 + x128 * x347 +
                 x128 * x381 - x188 * x321 * x72 + x188 * x81 + x361 * x399 -
                 x217 * x413 / x165 -
                 (kf_qss[17] * sc[11] * x323 + sc[7] * x339 + x402) / x165);
  amrex::Real x417 = -kf_qss[6] * x214 + x210 * (-kf_qss[10] + x22);
  amrex::Real x419 = -kf_qss[7] - kf_qss[8] + x389;
  const amrex::Real x420 = x345 * (-kf_qss[5] * x188 / x119 - kf_qss[6] * x225 +
                                   x169 * (-kf_qss[5] + x393) + x191 * x419 +
                                   x220 * x417 + x222 * (-kf_qss[5] + x393));
  amrex::Real x427 = -kf_qss[6] * x282 + kf_qss[6] * x378 + x20 * x53 / x33 -
                     x248 * x53 * (-kf_qss[10] + x22);
  const amrex::Real x428 =
    x382 *
    (-x167 * (-kf_qss[5] + x393) - x188 * x321 * x419 - x223 * x427 +
     x347 * (-kf_qss[5] + x393) + x361 * x417 + x381 * (-kf_qss[5] + x393) -
     x402 / x165 + x188 * (kf_qss[9] * sc[11] + x78) / x119);
  amrex::Real x429 = -kf_qss[29] - kf_qss[30];
  amrex::Real x432 = x210 * (-kf_qss[40] + x22) + x214 * x429;
  const amrex::Real x434 =
    x345 * (-kf_qss[25] * x188 / x119 + x169 * (-kf_qss[25] + x393) +
            x191 * (-kf_qss[39] + x389) + x220 * x432 +
            x222 * (-kf_qss[25] + x393) + x225 * x429);
  amrex::Real x438 = x20 * x53 / x33 - x248 * x53 * (-kf_qss[40] + x22) +
                     x282 * x429 - x378 * x429;
  const amrex::Real x439 =
    x382 * (-x167 * (-kf_qss[25] + x393) - x188 * x321 * (-kf_qss[39] + x389) -
            x223 * x438 + x347 * (-kf_qss[25] + x393) + x361 * x432 +
            x381 * (-kf_qss[25] + x393) - x402 / x165 +
            x188 * (kf_qss[11] * sc[11] + x78) / x119);
  amrex::Real x443 = -kf_qss[18] * x214 + x210 * (-kf_qss[23] + x22);
  amrex::Real x444 = -kf_qss[21] - kf_qss[9] * exp(-x92 - x94) - x103 + x389;
  const amrex::Real x445 = x345 * (-kf_qss[18] * x225 - x103 * x124 / x165 +
                                   x169 * (-kf_qss[17] + x393) + x191 * x444 +
                                   x220 * x443 + x222 * (-kf_qss[17] + x393));
  amrex::Real x450 =
    -kf_qss[18] * x282 + kf_qss[18] * x378 - x248 * x53 * (-kf_qss[23] + x22) -
    (kf_qss[20] * sc[7] + kf_qss[29] * x270 * x271) / x50 +
    x53 * (kf_qss[10] * sc[11] * x237 + kf_qss[15] * sc[7] * x245 + x20) / x33;
  x103 = x382 *
         (x103 * x322 / x165 - x167 * (-kf_qss[17] + x393) -
          x188 * x321 * x444 - x223 * x450 + x347 * (-kf_qss[17] + x393) +
          x361 * x443 + x381 * (-kf_qss[17] + x393) -
          (kf_qss[19] * sc[7] + kf_qss[24] * sc[11] * x333 + x402) / x165 +
          x188 * (kf_qss[22] * sc[11] + kf_qss[7] * sc[9] * x295 + x78) / x119);
  x183 = -kf_qss[22] * exp(-x43 - x94) - kf_qss[3] -
         8100000000000.0 * sc[1] * x87 * exp(x70) +
         x179 * (12.0 * x140 * x175 + 1.68 * x174) + 1.0935e+25 * x183;
  amrex::Real x457 = -kf_qss[19] * exp(-g_RT[7] + g_RT_qss[0] - x43) -
                     kf_qss[32] * x57 +
                     x154 * (1.68 * x139 * x140 + 12.0 * x140 * x143) -
                     1.92e+16 * x155 + 6.144e+31 * x392;
  amrex::Real x460 = x200 * (1.68 * x195 + 12.0 * x197) +
                     9.6774e+40 * x202 * x205 - 7.62e+20 * x30;
  amrex::Real x462 = -kf_qss[16] * x37 - kf_qss[20] * x45 - kf_qss[32];
  amrex::Real x463 = -kf_qss[16] * x37 * x53 / x33 + x210 * x460 + x214 * x462;
  const amrex::Real x466 =
    x345 *
    (-kf_qss[32] * x223 * x57 / x50 - kf_qss[32] * x68 / x165 + x169 * x457 +
     x183 * x191 + x220 * x463 + x222 * x457 + x225 * x462);
  amrex::Real x468 =
    (kf_qss[21] * sc[9] * x287 + kf_qss[3] * x270 * x299 -
     1.0935e+25 * x180 * x356 + 8100000000000.0 * x314 * x87 * exp(x316 + x70) -
     x84 * (12.0 * x140 * x175 + 1.68 * x174)) /
    x119;
  amrex::Real x471 = (kf_qss[23] * sc[11] * x228 - 9.6774e+40 * x202 * x373 -
                      x25 * (1.68 * x195 + 12.0 * x197) + 7.62e+20 * x371) /
                     x33;
  amrex::Real x472 = -kf_qss[30] * sc[9] * x253 / x50 - x248 * x460 * x53 +
                     x282 * x462 - x378 * x462 + x471 * x53;
  x457 = x382 *
         (-kf_qss[32] * x283 / (x165 * x54) - x167 * x457 - x183 * x188 * x321 +
          x188 * x468 - x223 * x472 + x347 * x457 + x361 * x463 + x381 * x457 -
          (1.92e+16 * sc[7] * x149 * x327 -
           x134 * (1.68 * x139 * x140 + 12.0 * x140 * x143) -
           6.144e+31 * x156 * x351) /
            x165);
  amrex::Real x475 =
    x209 * x217 * x22 * x52 * x66 / (x165 * ((x50) * (x50)) * ((x54) * (x54)));
  amrex::Real x476 = -kf_qss[11] * exp(-x110 - x94) + x389;
  const amrex::Real x477 = x345 * (x169 * (-kf_qss[24] + x393) + x191 * x476 +
                                   x222 * (-kf_qss[24] + x393) + x475);
  const amrex::Real x483 =
    (-x22 * x248 * x53 + x53 * (sc[11] * x232 + x20) / x33) / x54;
  const amrex::Real x486 =
    x382 *
    (-x167 * (-kf_qss[24] + x393) - x188 * x321 * x476 + x210 * x22 * x361 +
     x347 * (-kf_qss[24] + x393) + x381 * (-kf_qss[24] + x393) -
     x217 * x483 / x165 - x402 / x165 + x188 * (sc[9] * x304 + x78) / x119);
  const amrex::Real x489 =
    x345 * (x169 * x393 + x191 * (-kf_qss[37] + x389) + x222 * x393 + x475);
  const amrex::Real x494 =
    (-x22 * x248 * x53 -
     (kf_qss[20] * sc[4] + kf_qss[31] * sc[1] * x279 +
      2.0 * kf_qss[33] * ((sc[7])) * exp(-g_RT_qss[1] - x260)) /
       x50 +
     x53 * (kf_qss[15] * sc[4] * x245 + x20) / x33) /
    x54;
  x339 =
    x382 * (-x167 * x393 - x188 * x321 * (-kf_qss[37] + x389) +
            x210 * x22 * x361 + x347 * x393 + x381 * x393 - x217 * x494 / x165 -
            (kf_qss[19] * sc[4] +
             2.0 * kf_qss[27] * ((sc[7])) * exp(-g_RT_qss[0] - x260) +
             sc[1] * x339 + 3200000000000000.0 * x149 * x327 * x7 + x402) /
              x165 +
            x188 * (kf_qss[38] * sc[11] + x78) / x119);
  amrex::Real x497 = -kf_qss[38] * exp(-x94 - x97) + x185;
  const amrex::Real x500 =
    x345 *
    (-kf_qss[33] * x225 + x169 * (-kf_qss[27] + x161) + x191 * x497 +
     x220 * (-kf_qss[33] * x214 + x207 * x210) + x222 * (-kf_qss[27] + x161));
  amrex::Real x504 = -x207 * x248 * x53 + x375 * x53 / x33;
  const amrex::Real x506 =
    x382 *
    (-x167 * (-kf_qss[27] + x161) - x188 * x321 * x497 -
     x223 * (-kf_qss[33] * x282 + kf_qss[33] * x378 + x504) +
     x347 * (-kf_qss[27] + x161) + x352 +
     x361 * (-kf_qss[33] * x214 + x207 * x210) + x381 * (-kf_qss[27] + x161) +
     x188 * (kf_qss[37] * sc[9] * x311 + x359) / x119);
  x154 = -kf_qss[34] * x57 + x154 * (0.42 * x139 * x140 + 3.0 * x140 * x143) -
         4800000000000000.0 * x155 + 1.536e+31 * x392;
  x30 = x200 * (0.42 * x195 + 3.0 * x197) + 2.41935e+40 * x202 * x205 -
        1.905e+20 * x30;
  x155 =
    -1.5 * kf_qss[4] +
    2.73375e+24 * sc[1] * x180 * exp(-1434.1674975434501 * invT - 5.14 * logT) -
    2025000000000.0 * sc[1] * x87 * exp(x70) +
    x179 * (3.0 * x140 * x175 + 0.42 * x174);
  x179 = x345 * (-kf_qss[34] * x223 * x57 / x50 - kf_qss[34] * x225 -
                 kf_qss[34] * x68 / x165 + x154 * x169 + x154 * x222 +
                 x155 * x191 + x220 * (-kf_qss[34] * x214 + x210 * x30));
  x174 = (kf_qss[12] * sc[0] * x292 + kf_qss[21] * sc[5] * x287 +
          kf_qss[37] * sc[8] * x311 + kf_qss[3] * sc[5] * x269 * x299 +
          1.5 * kf_qss[4] * sc[9] * x269 * x299 + kf_qss[4] * x269 * x299 * x6 +
          kf_qss[7] * sc[4] * x295 + sc[6] * x304 - 2.73375e+24 * x180 * x356 +
          2025000000000.0 * x314 * x87 * exp(x316 + x70) -
          x84 * (3.0 * x140 * x175 + 0.42 * x174)) /
         x119;
  x175 = (-2.41935e+40 * x202 * x373 - x25 * (0.42 * x195 + 3.0 * x197) +
          1.905e+20 * x371) /
         x33;
  x180 = -kf_qss[34] * x282 + kf_qss[34] * x378 + x175 * x53 -
         x248 * x30 * x53 -
         (kf_qss[29] * sc[4] * x269 * x271 + kf_qss[30] * sc[5] * x253 +
          kf_qss[36] * sc[11] * x274 + kf_qss[6] * sc[0] * x256) /
           x50;
  x139 = x382 * (-kf_qss[34] * x283 / (x165 * x54) - x154 * x167 + x154 * x347 +
                 x154 * x381 - x155 * x188 * x321 + x174 * x188 - x180 * x223 +
                 x361 * (-kf_qss[34] * x214 + x210 * x30) -
                 (4800000000000000.0 * sc[7] * x149 * x327 -
                  x134 * (0.42 * x139 * x140 + 3.0 * x140 * x143) -
                  1.536e+31 * x156 * x351) /
                   x165);
  x143 = -kf_qss[35] - kf_qss[36];
  x149 = x143 * x41 * x52 / ((x50) * (x50)) + x207 * x210;
  x156 = -kf_qss[35] * x57 + x161;
  x161 = x345 * (-kf_qss[35] * x223 * x57 / x50 - kf_qss[35] * x68 / x165 +
                 x143 * x223 * x66 / ((x50) * (x50)) + x149 * x220 +
                 x156 * x169 + x156 * x222 + x185 * x191);
  x195 =
    x382 * (-kf_qss[35] * x283 / (x165 * x54) + x149 * x361 - x156 * x167 +
            x156 * x347 + x156 * x381 - x185 * x188 * x321 -
            x223 * (-x143 * x249 * x52 / ((x50) * (x50)) + x143 * x282 + x504) +
            x352 + x188 * (kf_qss[8] * sc[1] * x308 + x359) / x119);
  x197 = x345 * (x169 * x393 + x191 * x389 + x222 * x393 + x475);
  x87 = (kf_qss[11] * sc[3] + kf_qss[13] * sc[1] + kf_qss[22] * sc[4] +
         kf_qss[38] * sc[7] + kf_qss[9] * sc[2] +
         1350000000000.0 * refC * x7 * x87 * exp(x316 + x70) + x78) /
        x119;
  x200 = (kf_qss[10] * sc[4] * x237 + kf_qss[14] * sc[0] * x234 +
          kf_qss[23] * sc[5] * x228 + sc[1] * x15 * x28 + sc[6] * x232 + x20) /
         x33;
  x202 = (x200 * x53 - x22 * x248 * x53 -
          (kf_qss[18] * sc[1] * x265 + kf_qss[36] * sc[9] * x274) / x50) /
         x54;
  x205 =
    x382 *
    (-x167 * x393 - x188 * x321 * x389 + x188 * x87 + x210 * x22 * x361 +
     x347 * x393 + x381 * x393 - x202 * x217 / x165 -
     (kf_qss[17] * sc[1] * x323 + kf_qss[24] * sc[4] * x333 + x402) / x165);
  x217 = -kf_qss[28] * x57 + x393;
  x225 = -kf_qss[28] * x41 * x52 / ((x50) * (x50)) + x210 * x22;
  x169 = x345 *
         (-kf_qss[28] * x223 * x57 / x50 -
          kf_qss[28] * x223 * x66 / ((x50) * (x50)) - kf_qss[28] * x68 / x165 +
          x169 * x217 + x191 * x389 + x217 * x222 + x220 * x225);
  x53 = kf_qss[28] * x249 * x52 / ((x50) * (x50)) - kf_qss[28] * x282 +
        x20 * x53 / x33 - x22 * x248 * x53;
  x191 = x382 * (-kf_qss[28] * x283 / (x165 * x54) - x167 * x217 -
                 x188 * x321 * x389 + x217 * x347 + x217 * x381 - x223 * x53 +
                 x225 * x361 - x402 / x165 + x188 * x78 / x119);
  x220 = sc_qss[0] * x66 / ((x50) * (x50));
  x222 = (-kf_qss[31] * x220 + x379) / x54;
  x299 = (-sc_qss[0] * x66 / x50 + x283) / ((x54) * (x54));
  x327 = x299 * (-kf_qss[31] * x214 + x207 * x210);
  x345 = x299 * x399;
  x347 = (-kf_qss[6] * x220 + x427) / x54;
  x134 = x299 * x417;
  x351 = (x220 * x429 + x438) / x54;
  x352 = x299 * x432;
  x84 = (-kf_qss[18] * x220 + x450) / x54;
  x356 = x299 * x443;
  x361 = (-kf_qss[32] * sc_qss[0] * x57 / x50 + x220 * x462 + x472) / x54;
  x167 = x299 * x463;
  x25 = x210 * x22 * x299;
  x371 =
    (-kf_qss[33] * x220 - kf_qss[33] * x282 + kf_qss[33] * x378 + x504) / x54;
  x373 = x299 * (-kf_qss[33] * x214 + x207 * x210);
  x378 = (-kf_qss[34] * sc_qss[0] * x57 / x50 - kf_qss[34] * x220 + x180) / x54;
  x210 = x299 * (-kf_qss[34] * x214 + x210 * x30);
  x52 = (-kf_qss[35] * sc_qss[0] * x57 / x50 + x143 * x220 -
         x143 * x249 * x52 / ((x50) * (x50)) + x143 * x282 + x504) /
        x54;
  x214 = x149 * x299;
  x57 = (-kf_qss[28] * sc_qss[0] * x57 / x50 -
         kf_qss[28] * sc_qss[0] * x66 / ((x50) * (x50)) + x53) /
        x54;
  x282 = x225 * x299;
  x379 = sc_qss[0] * (kf_qss[25] * sc[3] + kf_qss[5] * sc[2] + 10.0e-13) /
         ((x119) * (x119));
  x186 = x186 * x321 + x186 * x379 - (kf_qss[12] * sc[9] * x292 + x359) / x119;
  x381 = x321 * x72 + x379 * x72 - x81;
  x72 = -kf_qss[5] * sc_qss[0] / x119 + x321 * x419 + x379 * x419 -
        (kf_qss[9] * sc[11] + x78) / x119;
  x392 = -kf_qss[25] * sc_qss[0] / x119 + x321 * (-kf_qss[39] + x389) +
         x379 * (-kf_qss[39] + x389) - (kf_qss[11] * sc[11] + x78) / x119;
  x393 = x321 * x444 + x379 * x444 -
         (kf_qss[22] * sc[11] + kf_qss[7] * sc[9] * x295 + x78) / x119;
  x399 = x183 * x321 + x183 * x379 - x468;
  x304 = x321 * x476 + x379 * x476 - (sc[9] * x304 + x78) / x119;
  x402 = x321 * (-kf_qss[37] + x389) + x379 * (-kf_qss[37] + x389) -
         (kf_qss[38] * sc[11] + x78) / x119;
  x81 = x321 * x497 + x379 * x497 - (kf_qss[37] * sc[9] * x311 + x359) / x119;
  x417 = x155 * x321 + x155 * x379 - x174;
  x185 = x185 * x321 + x185 * x379 - (kf_qss[8] * sc[1] * x308 + x359) / x119;
  x359 = x321 * x389 + x379 * x389 - x87;
  x119 = x321 * x389 + x379 * x389 - x78 / x119;
  x321 = sc_qss[1] * x207 * x209 + x207 * x248 -
         (kf_qss[14] * sc[11] * x234 + x375) / x33;
  x15 =
    sc_qss[1] * x209 * x397 + x248 * x397 - (sc[11] * x15 * x28 + x20) / x33;
  x28 = sc_qss[1] * x209 * (-kf_qss[10] + x22) - x20 / x33 +
        x248 * (-kf_qss[10] + x22);
  x389 = sc_qss[1] * x209 * (-kf_qss[40] + x22) - x20 / x33 +
         x248 * (-kf_qss[40] + x22);
  x397 = sc_qss[1] * x209 * (-kf_qss[23] + x22) + x248 * (-kf_qss[23] + x22) -
         (kf_qss[10] * sc[11] * x237 + kf_qss[15] * sc[7] * x245 + x20) / x33;
  x78 = -kf_qss[16] * sc_qss[1] * x37 / x33 + sc_qss[1] * x209 * x460 +
        x248 * x460 - x471;
  x22 = sc_qss[1] * x209 * x22 + x22 * x248;
  x232 = x22 - (sc[11] * x232 + x20) / x33;
  x419 = x22 - (kf_qss[15] * sc[4] * x245 + x20) / x33;
  x207 = sc_qss[1] * x207 * x209 + x207 * x248 - x375 / x33;
  x209 = sc_qss[1] * x209 * x30 - x175 + x248 * x30;
  x248 = -x200 + x22;
  x33 = -x20 / x33 + x22;
  x117 = 23000.0 * exp(-1648.03458051046 * invT + 1.05 * logT - x117);
  x375 = sc[0] * x117;
  x20 = 6000000.0 * sc[5] * x37;
  x427 = -717.08374877172594 * invT - 2.57 * logT;
  x70 =
    1238532.11009174 * x7 * exp(-130.83633310922701 * invT - 0.48 * logT + x70);
  x429 = -0.29097730287517898 * x75 + 0.43429448190325198 * log(x70) - 0.4;
  x432 = 1.7777777777777777 * ((x429) * (x429)) /
           ((-0.18666666666666668 * x429 - 0.73540532268950665 * x75 + 1) *
            (-0.18666666666666668 * x429 - 0.73540532268950665 * x75 + 1)) +
         1.0;
  x438 = exp(M_LN10 * (0.43429448190325198 * x75 / x432)) / (x70 + 1.0);
  x443 = 1350000000000.0 * sc[1] * x438 * x7 * exp(x427);
  x444 = 80000000.0 * sc[1] * x107;
  x450 = exp(-8554.6833186802305 * invT - logT);
  x183 = 2244000000000.0 * sc[5] * x450;
  x460 = 187000000000.0 * x450 * x6;
  x462 = -3286.0048277048204 * invT - 4.82 * logT;
  x11 =
    235185185185185.0 * x7 * exp(1811.57999689699 * invT - 0.454 * logT + x11);
  x463 = -0.29097730287517898 * x17 + 0.43429448190325198 * log(x11) - 0.4;
  x468 = 1.7777777777777777 * ((x463) * (x463)) /
           ((-0.73540532268950665 * x17 - 0.18666666666666668 * x463 + 1) *
            (-0.73540532268950665 * x17 - 0.18666666666666668 * x463 + 1)) +
         1.0;
  x471 = exp(M_LN10 * (0.43429448190325198 * x17 / x468)) / (x11 + 1.0);
  x472 = 1.27e+20 * refC * x471 * x7 * exp(x10 + x462);
  x475 = 39000000.0 * exp(-1781.38699694871 * invT - x92 - x94);
  x476 = sc[4] * x475;
  x497 = 7600000.0 * sc[3] * exp(-201.28666632188799 * invT);
  x504 = 13200000.0 * exp(
                        -g_RT[3] + g_RT[4] - g_RT_qss[0] + g_RT_qss[2] -
                        754.82499870708 * invT);
  x154 = sc[4] * x504;
  x30 = 100000000.0 * exp(-20128.666632188801 * invT - x110 - x94);
  x155 = sc[6] * x30;
  x174 = 18000000.0 * sc[3] * exp(-452.89499922424801 * invT);
  x175 = 25010000.0 * sc[5] * x45;
  x180 = 3430.0 * exp(224.93784961470999 * invT + 1.18 * logT - x43 - x94);
  x143 = sc[5] * x180;
  x149 = 32000000.0 * sc[8] * exp(286.83349950869001 * invT);
  x94 = 0.00332 * exp(-2948.8496616156599 * invT + 2.81 * logT - x94 - x97);
  x156 = sc[8] * x94;
  x87 = 16000000.0 * sc[8] * exp(286.83349950869001 * invT);
  x200 = 28000000.0 *
         exp(g_RT[0] - g_RT[1] + g_RT[3] - g_RT[6] - 537.43539907944103 * invT);
  x217 = 0.05 * exp(-3165.23282791169 * invT + 2.67 * logT);
  x225 = 216.0 * sc[4] * exp(-1726.0331637101899 * invT + 1.51 * logT);
  x53 = 0.5 * exp(-3638.25649376812 * invT + 2.0 * logT);
  x220 = sc[0] * x53;
  x299 = sc_qss[0] * x124;
  x379 = -x299 - x322;
  x22 = g_RT[0] - 2.0 * g_RT[1];
  const amrex::Real x775 = refC * exp(-logT + x22);
  x165 = x382 * (x188 * x322 - x223 * x283 - x342 / x165);
  x188 = 660.0 * exp(
                   g_RT[0] - g_RT[1] + g_RT[7] - g_RT[8] -
                   5454.8686573231598 * invT + 1.62 * logT);
  x50 = (-sc_qss[0] * x66 / x50 + x283) / x54;
  x54 = 117906976744186.0 * x7 * exp(-2390.2791625723985 * invT - 4.92 * logT);
  x66 = log(
    exp(-10300.0 * invT) + 0.068 * exp(-0.0050761421319797002 * T) +
    0.932 * exp(-0.00064935064935064902 * T));
  x223 = -0.29097730287517898 * x66 + 0.43429448190325198 * log(x54) - 0.4;
  x283 = -0.14 * x223 - 0.55155399201712996 * x66 + 0.75;
  x342 = exp(
           M_LN10 * (0.43429448190325198 * x66 /
                     (((x223) * (x223)) / ((x283) * (x283)) + 1.0))) /
         (x54 + 1.0);
  x382 = -42446.325760628104 * invT - 3.42 * logT;
  const amrex::Real x815 = exp(-g_RT[0] + g_RT[11] - g_RT[9] + x382);
  x54 = exp(
          M_LN10 * (0.43429448190325198 * x66 /
                    (((x223) * (x223)) / ((x283) * (x283)) + 1.0))) /
        ((x54 + 1.0) * (x54 + 1.0));
  const amrex::Real x822 =
    x314 * x7 *
    exp(
      -g_RT[0] + g_RT[11] - g_RT[9] - 44836.604923200604 * invT - 8.34 * logT);
  amrex::Real x826 = x140 * ((x223) * (x223)) / ((x283) * (x283) * (x283));
  amrex::Real x827 = 4.0 * x140 * x223 / ((x283) * (x283)) + 0.56 * x826;
  x66 = x66 / ((((x223) * (x223)) / ((x283) * (x283)) + 1.0) *
               (((x223) * (x223)) / ((x283) * (x283)) + 1.0));
  const amrex::Real x831 = 5070000000000000.0 * x314 * x342 * x815;
  const amrex::Real x833 = x66 * x7 * x831;
  const amrex::Real x837 = sc[0] * sc[9] * exp(x382);
  const amrex::Real x839 = 5070000000000000.0 * x342 * x7;
  const amrex::Real x841 = x66 * x837 * x839;
  const amrex::Real x843 =
    sc[0] * sc[9] * x7 * exp(-44836.604923200604 * invT - 8.34 * logT);
  const amrex::Real x845 =
    1.014e+16 * x314 * x342 * x815 - 1.014e+16 * x342 * x837 -
    1.19557674418605e+30 * x54 * x822 + 1.19557674418605e+30 * x54 * x843 -
    x827 * x833 + x827 * x841;
  const amrex::Real x853 =
    73400000.0 * sc[1] * x186 + 15000000.0 * sc[2] * x222 -
    15000000.0 * sc[2] * x327 - 15000000.0 * sc[9] * x256 -
    73400000.0 * sc[9] * x292 - sc[9] * x839 * exp(x382);
  const amrex::Real x856 = 28000000.0 * exp(-537.43539907944103 * invT);
  const amrex::Real x858 =
    0.05 *
    exp(-g_RT[0] + g_RT[1] - 3165.23282791169 * invT + 2.67 * logT - x92);
  const amrex::Real x861 = exp(-2.0 * logT);
  const amrex::Real x866 = exp(-logT);
  const amrex::Real x869 = 0.5 * exp(
                                   -g_RT[0] + g_RT[1] + g_RT[7] - g_RT_qss[0] -
                                   3638.25649376812 * invT + 2.0 * logT);
  const amrex::Real x871 =
    216.0 *
    exp(-g_RT[0] + g_RT[1] - 1726.0331637101899 * invT + 1.51 * logT - x43);
  const amrex::Real x877 =
    23000.0 * exp(-1648.03458051046 * invT + 1.05 * logT);
  const amrex::Real x880 =
    660.0 * sc[8] * exp(-5454.8686573231598 * invT + 1.62 * logT);
  amrex::Real x890 = 2.0 * x140 * x223 / ((x283) * (x283)) + 0.28 * x826;
  const amrex::Real x894 =
    -5070000000000000.0 * x342 * x837 - 5.9778837209302301e+29 * x54 * x822 +
    5.9778837209302301e+29 * x54 * x843 + x831 - x833 * x890 + x841 * x890;
  const amrex::Real x895 =
    -1000000.0 * sc[0] * x775 + 1000000.0 * ((sc[1]) * (sc[1])) * x866 + x894;
  const amrex::Real x901 =
    73400000.0 * sc[1] * x381 - 15000000.0 * sc[2] * x345 +
    15000000.0 * sc[2] * x413 - 73400000.0 * x299 - 73400000.0 * x322;
  const amrex::Real x902 = sc_qss[1] * x41;
  const amrex::Real x934 =
    216.0 * sc[0] * exp(-1726.0331637101899 * invT + 1.51 * logT);
  x66 = x66 * (12.0 * x140 * x223 / ((x283) * (x283)) + 1.68 * x826);
  const amrex::Real x963 =
    73400000.0 * sc[1] * x399 - 15000000.0 * sc[2] * x167 +
    15000000.0 * sc[2] * x361 + 3.042e+16 * x314 * x342 * x815 -
    3.042e+16 * x342 * x837 - 3.5867302325581401e+30 * x54 * x822 +
    3.5867302325581401e+30 * x54 * x843 - x66 * x7 * x831 + x66 * x837 * x839;
  const amrex::Real x973 =
    70000000.0 * sc[0] * x25 - 15000000.0 * sc[2] * x25 + x895;
  const amrex::Real x994 =
    660.0 * sc[1] * exp(-5454.8686573231598 * invT + 1.62 * logT);
  x223 = 3.0 * x140 * x223 / ((x283) * (x283)) + 0.42 * x826;
  x283 = -15000000.0 * sc[0] * x256 - 73400000.0 * sc[0] * x292 -
         sc[0] * x839 * exp(x382) + 73400000.0 * sc[1] * x417 -
         15000000.0 * sc[2] * x210 + 15000000.0 * sc[2] * x378 - x223 * x833 +
         x223 * x841 + 7605000000000000.0 * x314 * x342 * x815 -
         7605000000000000.0 * x342 * x837 -
         8.9668255813953501e+29 * x54 * x822 +
         8.9668255813953501e+29 * x54 * x843;
  x826 = x438 * exp(x316 + x427);
  amrex::Real x1062 =
    128000.0 * x7 * exp(-618.95649893980499 * invT - 2.34 * logT);
  amrex::Real x1063 =
    -0.29097730287517898 * x131 + 0.43429448190325198 * log(x1062) - 0.4;
  amrex::Real x1065 = -0.14 * x1063 - 0.55155399201712996 * x131 + 0.75;
  const amrex::Real x1067 =
    ((x1063) * (x1063)) * x140 / ((x1065) * (x1065) * (x1065));
  const amrex::Real x1069 = x1063 * x140 / ((x1065) * (x1065));
  amrex::Real x1071 =
    x131 / ((((x1063) * (x1063)) / ((x1065) * (x1065)) + 1.0) *
            (((x1063) * (x1063)) / ((x1065) * (x1065)) + 1.0));
  amrex::Real x1078 =
    exp(
      M_LN10 * (0.43429448190325198 * x131 /
                (((x1063) * (x1063)) / ((x1065) * (x1065)) + 1.0))) /
    (x1062 + 1.0);
  const amrex::Real x1079 =
    3200000000000000.0 * x1078 * exp(-618.95649893980499 * invT - 3.14 * logT);
  const amrex::Real x1082 =
    sc[1] * x1071 * x1079 * x165 * x7 * (0.56 * x1067 + 4.0 * x1069);
  amrex::Real x1084 =
    1.7777777777777777 * x429 /
    ((-0.18666666666666668 * x429 - 0.73540532268950665 * x75 + 1) *
     (-0.18666666666666668 * x429 - 0.73540532268950665 * x75 + 1));
  x429 = 2.3703703703703702 * ((x429) * (x429)) /
         ((-0.18666666666666668 * x429 - 0.73540532269 * x75 + 1) *
          (-0.18666666666666668 * x429 - 0.73540532269 * x75 + 1) *
          (-0.18666666666666668 * x429 - 0.73540532269 * x75 + 1));
  const amrex::Real x1086 = 4.0 * x1084 * x140 + 0.56 * x140 * x429;
  const amrex::Real x1089 =
    1350000000000.0 * x314 * x7 * x75 * x826 / ((x432) * (x432));
  const amrex::Real x1092 =
    sc[1] * x1078 * x165 * exp(-618.95649893980499 * invT - 3.14 * logT);
  x70 = exp(M_LN10 * (0.43429448190325198 * x75 / x432)) /
        ((x70 + 1.0) * (x70 + 1.0));
  x316 = x314 * x7 * exp(-1565.003830652674 * invT - 5.62 * logT + x316);
  const amrex::Real x1102 =
    refC * sc[7] * x1078 * exp(-618.95649893980499 * invT - 3.14 * logT + x326);
  const amrex::Real x1105 =
    sc[1] * x379 * x7 * exp(-1565.003830652674 * invT - 5.62 * logT);
  x1078 = 3200000000000000.0 * x1078 * x204 *
          exp(-618.95649893980499 * invT - 3.14 * logT + x326);
  const amrex::Real x1109 =
    sc[7] * x1071 * x1078 * (0.56 * x1067 + 4.0 * x1069);
  x131 = exp(
           M_LN10 * (0.43429448190325198 * x131 /
                     (((x1063) * (x1063)) / ((x1065) * (x1065)) + 1.0))) /
         ((x1062 + 1.0) * (x1062 + 1.0));
  x326 = sc[7] * x204 * exp(-1237.91299787961 * invT - 5.48 * logT + x326);
  x1062 = sc[1] * x165 * x7 * exp(-1237.91299787961 * invT - 5.48 * logT);
  x1063 = sc[1] * x379 * x438 * exp(x427);
  x75 = 1350000000000.0 * x1063 * x7 * x75 / ((x432) * (x432));
  x432 = sc[1] * x1079 * x7;
  x1065 = refC * sc[6];
  const amrex::Real x1142 = exp(-g_RT[1] - 0.86 * logT - x110);
  const amrex::Real x1147 = sc[1] * sc[3] * exp(-0.86 * logT);
  const amrex::Real x1149 = 2800000.0 * x1065 * x1142 - 2800000.0 * x1147;
  const amrex::Real x1150 = refC * exp(-g_RT[1] - logT - x92);
  const amrex::Real x1151 = sc[4] * x1150;
  const amrex::Real x1154 = sc[1] * sc[2] * x866;
  const amrex::Real x1156 = 1000000.0 * x1151 - 1000000.0 * x1154;
  amrex::Real x1158 = exp(g_RT[1] - 8554.6833186802305 * invT - logT - x286);
  const amrex::Real x1160 =
    -374000000000.0 * x1158 * x270 + 374000000000.0 * x379 * x450;
  const amrex::Real x1163 = refC * exp(-g_RT[1] - 2.0 * logT - x43);
  const amrex::Real x1166 = sc[1] * sc[4] * x861;
  const amrex::Real x1168 =
    sc[2] * x217 + 28000000.0 * sc[3] * x222 - 28000000.0 * sc[3] * x327 +
    16060000000.0 * sc[5] * x1163 - 16060000000.0 * x1166;
  amrex::Real x1172 =
    195039370078.74 * x7 * exp(-1035.1166815603121 * invT - 4.13 * logT);
  amrex::Real x1173 = log(
    exp(-6964.0 * invT) + 0.217 * exp(-0.0135135135135135 * T) +
    0.783 * exp(-0.00034002040122407302 * T));
  amrex::Real x1175 =
    -0.29097730287517898 * x1173 + 0.43429448190325198 * log(x1172) - 0.4;
  amrex::Real x1176 = -0.55155399201712996 * x1173 - 0.14 * x1175 + 0.75;
  amrex::Real x1179 = ((x1175) * (x1175)) / ((x1176) * (x1176)) + 1.0;
  const amrex::Real x1182 =
    exp(M_LN10 * (0.43429448190325198 * x1173 / x1179)) / (x1172 + 1.0);
  const amrex::Real x1185 = -1227.8486645635121 * invT - 4.76 * logT;
  const amrex::Real x1186 = exp(-g_RT[1] + x1185 - x97);
  x1172 = exp(M_LN10 * (0.43429448190325198 * x1173 / x1179)) /
          ((x1172 + 1.0) * (x1172 + 1.0));
  const amrex::Real x1192 =
    sc[1] * sc[7] * x7 * exp(-2262.965346123824 * invT - 8.89 * logT);
  amrex::Real x1197 =
    ((x1175) * (x1175)) * x140 / ((x1176) * (x1176) * (x1176));
  amrex::Real x1200 = 4.0 * x1175 * x140 / ((x1176) * (x1176)) + 0.56 * x1197;
  const amrex::Real x1201 = refC * sc[8] * x1182 * x1186;
  const amrex::Real x1204 =
    2.477e+21 * x1173 * x1201 * x7 / ((x1179) * (x1179));
  const amrex::Real x1210 = sc[1] * sc[7] * x1182 * exp(x1185);
  x1173 = 2.477e+21 * x1173 * x1210 * x7 / ((x1179) * (x1179));
  x1179 = refC * sc[8] * x7 *
          exp(-g_RT[1] - 2262.965346123824 * invT - 8.89 * logT - x97);
  amrex::Real x1217 = 4.954e+21 * refC * sc[8] * x1182 * x1186 -
                      4.954e+21 * sc[1] * sc[7] * x1182 * exp(x1185) -
                      9.6622503937007901e+32 * x1172 * x1179 +
                      9.6622503937007901e+32 * x1172 * x1192 + x1173 * x1200 -
                      x1200 * x1204;
  const amrex::Real x1218 = 70000000.0 * sc[0] * x222 -
                            70000000.0 * sc[0] * x327 + x1217 + x165 * x53 -
                            x220 * x346 + x220 * x383 + 70000000.0 * x50;
  amrex::Real x1225 =
    1.7777777777777777 * x463 /
    ((-0.73540532268950665 * x17 - 0.18666666666666668 * x463 + 1) *
     (-0.73540532268950665 * x17 - 0.18666666666666668 * x463 + 1));
  x463 = 2.3703703703703702 * ((x463) * (x463)) /
         ((-0.73540532269 * x17 - 0.18666666666666668 * x463 + 1) *
          (-0.73540532269 * x17 - 0.18666666666666668 * x463 + 1) *
          (-0.73540532269 * x17 - 0.18666666666666668 * x463 + 1));
  amrex::Real x1229 = sc[11] * sc[1] * exp(x462);
  amrex::Real x1230 = 1.27e+20 * x17 * x471 * x7 / ((x468) * (x468));
  x17 = exp(M_LN10 * (0.43429448190325198 * x17 / x468)) /
        ((x11 + 1.0) * (x11 + 1.0));
  x11 = -x249 - x902;
  x468 = x11 * x204 * exp(-4760.4296585126494 * invT - 10.094 * logT + x10);
  amrex::Real x1238 = refC * x11 * exp(x10 + x462);
  amrex::Real x1240 =
    sc[11] * sc[1] * x7 * exp(-4760.4296585126494 * invT - 10.094 * logT);
  x10 = 2.54e+20 * refC * x11 * x471 * exp(x10 + x462) +
        x1229 * x1230 * (4.0 * x1225 * x140 + 0.56 * x140 * x463) -
        2.54e+20 * x1229 * x471 -
        x1230 * x1238 * (4.0 * x1225 * x140 + 0.56 * x140 * x463) +
        5.9737037037037004e+34 * x1240 * x17 -
        5.9737037037037004e+34 * x17 * x468;
  const amrex::Real x1243 =
    30000000.0 * sc[4] * x222 - 30000000.0 * sc[4] * x327 -
    20000000.0 * sc[4] * x346 + 20000000.0 * sc[4] * x383 + x10 + x117 * x379 +
    x186 * x375 + x321 * x472;
  const amrex::Real x1250 = 28000000.0 * sc[4] * x271;
  const amrex::Real x1267 = 47.6 * exp(
                                     g_RT[10] + g_RT[1] - g_RT[4] - g_RT[9] -
                                     35.225166606330397 * invT + 1.228 * logT);
  const amrex::Real x1275 = 3970000.0 * exp(-337.65838275496702 * invT);
  const amrex::Real x1277 =
    -sc[6] * x1275 - sc[6] * x856 -
    134000000.0 * sc[6] * exp(-319.54258278599701 * invT);
  const amrex::Real x1278 = sc[0] + 1.5 * sc[10] + sc[11] + sc[1] + sc[2] +
                            sc[4] + sc[6] + sc[7] + sc[8] + 0.75 * sc[9];
  const amrex::Real x1279 = 2800000.0 * x1278 * exp(-0.86 * logT);
  const amrex::Real x1282 = 9380000.0 * exp(-0.76 * logT);
  const amrex::Real x1285 = exp(-1.72 * logT);
  const amrex::Real x1291 = 83000000.0 * exp(-7252.86180424342 * invT);
  const amrex::Real x1293 = -375000000.0 * sc[12] * sc[3] * x1285 -
                            300000000.0 * ((sc[3]) * (sc[3])) * x1285 -
                            sc[3] * sc[5] * x1282 - sc[3] * x1279 -
                            sc[3] * x1291;
  const amrex::Real x1296 = -2244000000000.0 * refCinv * sc[5] * sc[9] * x1158 -
                            187000000000.0 * refCinv * sc[9] * x1158 * x6 +
                            x183 * x381 + x381 * x460;
  const amrex::Real x1300 =
    -70000000.0 * sc[0] * x345 + 70000000.0 * sc[0] * x413 -
    2.477e+21 * sc[7] * x1182 * x7 * exp(x1185) - 70000000.0 * sc[7] * x279 -
    sc[7] * x869 + x128 * x220 - x220 * x400;
  const amrex::Real x1301 = 16060000000.0 * sc[0] + 22000000000.0 * sc[10] +
                            22000000000.0 * sc[3] + 80300000000.0 * sc[5] +
                            22000000000.0 * sc[9] + 22000000000.0 * x2;
  const amrex::Real x1303 = -28000000.0 * sc[3] * x345 +
                            28000000.0 * sc[3] * x413 - sc[4] * x1301 * x861 -
                            sc[4] * x858;
  const amrex::Real x1311 =
    84300000.0 * exp(g_RT[11] + g_RT[1] - g_RT[2] - g_RT[7]);
  const amrex::Real x1315 =
    -sc[11] * x1311 - 30000000.0 * sc[11] * x265 - 20000000.0 * sc[11] * x323 -
    1.27e+20 * sc[11] * x471 * x7 * exp(x462) - sc[11] * x877 +
    20000000.0 * sc[4] * x128 - 30000000.0 * sc[4] * x345 -
    20000000.0 * sc[4] * x400 + 30000000.0 * sc[4] * x413 + x15 * x472 +
    x375 * x381;
  amrex::Real x1316 = 2.0 * x1175 * x140 / ((x1176) * (x1176)) + 0.28 * x1197;
  const amrex::Real x1322 = -4.83112519685039e+32 * x1172 * x1179 +
                            4.83112519685039e+32 * x1172 * x1192 +
                            x1173 * x1316 + 2.477e+21 * x1201 - x1204 * x1316 -
                            2.477e+21 * x1210;
  const amrex::Real x1325 =
    x1229 * x1230 * (2.0 * x1225 * x140 + 0.28 * x140 * x463) -
    1.27e+20 * x1229 * x471 -
    x1230 * x1238 * (2.0 * x1225 * x140 + 0.28 * x140 * x463) +
    1.27e+20 * x1238 * x471 + 2.9868518518518502e+34 * x1240 * x17 -
    2.9868518518518502e+34 * x17 * x468;
  const amrex::Real x1328 = 2.0 * x1084 * x140 + 0.28 * x140 * x429;
  const amrex::Real x1334 = sc[1] * x1071 * x1079 * x165 * x7;
  x1071 = sc[7] * x1071 * x1078;
  const amrex::Real x1344 =
    -187000000000.0 * x1158 * x270 + 187000000000.0 * x379 * x450;
  const amrex::Real x1347 =
    22000000000.0 * sc[5] * x1163 - 22000000000.0 * x1166;
  const amrex::Real x1349 =
    2000000.0 * sc[0] * x775 - 2000000.0 * ((sc[1]) * (sc[1])) * x866 -
    sc[1] * x1079 * x165 + 4.096e+20 * x1062 * x131 - 1350000000000.0 * x1063 -
    x1071 * (0.28 * x1067 + 2.0 * x1069) - x1089 * x1328 +
    3200000000000000.0 * x1102 + 1.67201834862385e+18 * x1105 * x70 +
    500000.0 * x1151 - 500000.0 * x1154 - 4.096e+20 * x131 * x326 +
    x1328 * x75 + x1334 * (0.28 * x1067 + 2.0 * x1069) + x1344 + x1347 +
    1350000000000.0 * x314 * x826 - 1.67201834862385e+18 * x316 * x70;
  const amrex::Real x1354 =
    3970000.0 *
    exp(-g_RT[1] + g_RT[2] + g_RT[5] - g_RT[6] - 337.65838275496702 * invT);
  const amrex::Real x1356 =
    83000000.0 *
    exp(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[4] - 7252.86180424342 * invT);
  const amrex::Real x1369 =
    30000000.0 * sc[2] * x72 - 30000000.0 * x299 - 30000000.0 * x322;
  const amrex::Real x1378 =
    -30000000.0 * sc[4] * x134 + 30000000.0 * sc[4] * x347 -
    20000000.0 * sc[4] * x420 + 20000000.0 * sc[4] * x428 + 84300000.0 * sc[7] +
    x28 * x472 + x375 * x72;
  const amrex::Real x1389 = exp(-g_RT[1] - 1.72 * logT - x110);
  const amrex::Real x1398 = -375000000.0 * sc[12] * sc[1] * x1285 -
                            600000000.0 * sc[1] * sc[3] * x1285 -
                            sc[1] * sc[5] * x1282 - sc[1] * x1279 -
                            sc[1] * x1291 + 300000000.0 * x1065 * x1389;
  const amrex::Real x1405 =
    30000000.0 * sc[4] * x351 - 30000000.0 * sc[4] * x352 -
    20000000.0 * sc[4] * x434 + 20000000.0 * sc[4] * x439 + x375 * x392 +
    x389 * x472;
  const amrex::Real x1407 =
    47.6 * sc[9] * exp(-35.225166606330397 * invT + 1.228 * logT);
  const amrex::Real x1409 =
    sc[4] * exp(-g_RT[1] + 2.0 * g_RT[4] - g_RT[6] - 319.54258278599701 * invT);
  const amrex::Real x1434 =
    20000000.0 * sc[4] * x103 - 30000000.0 * sc[4] * x356 -
    20000000.0 * sc[4] * x445 + 30000000.0 * sc[4] * x84 + 20000000.0 * x165 +
    x375 * x393 + x397 * x472 + 30000000.0 * x50;
  const amrex::Real x1436 = -sc[1] * x1301 * x861 - sc[1] * x858 -
                            28000000.0 * sc[3] * x356 +
                            28000000.0 * sc[3] * x84;
  const amrex::Real x1442 = 12.0 * x1084 * x140 + 1.68 * x140 * x429;
  const amrex::Real x1463 = 9380000.0 * exp(-g_RT[1] - 0.76 * logT - x110);
  const amrex::Real x1472 = 30000000.0 * sc[2] * x399 -
                            2244000000000.0 * x1158 * x270 + x183 * x399 +
                            2244000000000.0 * x379 * x450 + x399 * x460;
  const amrex::Real x1473 =
    12.0 * x1175 * x140 / ((x1176) * (x1176)) + 1.68 * x1197;
  const amrex::Real x1480 = -70000000.0 * sc[0] * x167 +
                            70000000.0 * sc[0] * x361 -
                            2.8986751181102398e+33 * x1172 * x1179 +
                            2.8986751181102398e+33 * x1172 * x1192 +
                            x1173 * x1473 + 1.4862e+22 * x1201 - x1204 * x1473 -
                            1.4862e+22 * x1210 + x220 * x457 - x220 * x466;
  const amrex::Real x1484 =
    -28000000.0 * sc[3] * x167 + 28000000.0 * sc[3] * x361 +
    80300000000.0 * sc[5] * x1163 + x1163 * x1301 - 80300000000.0 * x1166;
  const amrex::Real x1492 =
    -30000000.0 * sc[4] * x167 + 30000000.0 * sc[4] * x361 +
    20000000.0 * sc[4] * x457 - 20000000.0 * sc[4] * x466 +
    x1229 * x1230 * (12.0 * x1225 * x140 + 1.68 * x140 * x463) -
    7.62e+20 * x1229 * x471 -
    x1230 * x1238 * (12.0 * x1225 * x140 + 1.68 * x140 * x463) +
    7.62e+20 * x1238 * x471 + 1.7921111111111101e+35 * x1240 * x17 -
    1.7921111111111101e+35 * x17 * x468 + x375 * x399 + x472 * x78;
  const amrex::Real x1501 =
    -sc[1] * x1275 - sc[1] * x856 -
    134000000.0 * sc[1] * exp(-319.54258278599701 * invT);
  const amrex::Real x1507 =
    375000000.0 * refC * sc[12] * x1389 + 300000000.0 * refC * sc[3] * x1389 +
    refC * sc[5] * x1463 + 2800000.0 * refC * x1142 * x1278;
  const amrex::Real x1511 = -28000000.0 * sc[3] * x25;
  const amrex::Real x1516 =
    -20000000.0 * sc[4] * x477 + 30000000.0 * sc[4] * x483 +
    20000000.0 * sc[4] * x486 + x232 * x472 + x304 * x375;
  const amrex::Real x1521 = -70000000.0 * sc[0] * x25 -
                            30000000.0 * sc[4] * x25 + x1149 + x1322 + x1325 +
                            x1349 - x20 * x25;
  const amrex::Real x1539 =
    84300000.0 * sc[2] + 20000000.0 * sc[4] * x339 - 20000000.0 * sc[4] * x489 +
    30000000.0 * sc[4] * x494 + x375 * x402 + x419 * x472;
  const amrex::Real x1541 =
    70000000.0 * sc[0] * x494 - 2.477e+21 * sc[1] * x1182 * x7 * exp(x1185) -
    70000000.0 * sc[1] * x279 - sc[1] * x869 + x220 * x339 - x220 * x489;
  const amrex::Real x1561 =
    30000000.0 * sc[4] * x371 - 30000000.0 * sc[4] * x373 -
    20000000.0 * sc[4] * x500 + 20000000.0 * sc[4] * x506 + x10 + x207 * x472 +
    x375 * x81;
  const amrex::Real x1563 =
    70000000.0 * sc[0] * x371 - 70000000.0 * sc[0] * x373 +
    2.477e+21 * x1182 * x1186 * x204 + x1217 - x220 * x500 + x220 * x506;
  const amrex::Real x1564 =
    -38000000.0 * sc[1] * x207 + 8.192e+20 * x1062 * x131 -
    2700000000000.0 * x1063 + x1082 - x1086 * x1089 + x1086 * x75 -
    6400000000000000.0 * x1092 + 6400000000000000.0 * x1102 +
    3.3440366972477102e+18 * x1105 * x70 - x1109 + x1156 + x1160 -
    8.192e+20 * x131 * x326 + 2700000000000.0 * x314 * x826 -
    3.3440366972477102e+18 * x316 * x70;
  const amrex::Real x1565 =
    47.6 * sc[4] * exp(-35.225166606330397 * invT + 1.228 * logT);
  x1084 = 3.0 * x1084 * x140 + 0.42 * x140 * x429;
  x429 = 3.0 * x1175 * x140 / ((x1176) * (x1176)) + 0.42 * x1197;
  x1175 = -70000000.0 * sc[0] * x210 + 70000000.0 * sc[0] * x378 -
          7.2466877952755894e+32 * x1172 * x1179 +
          7.2466877952755894e+32 * x1172 * x1192 + x1173 * x429 +
          3.7155e+21 * x1201 - x1204 * x429 - 3.7155e+21 * x1210 + x139 * x220 -
          x179 * x220;
  x1158 = 30000000.0 * sc[2] * x417 - 2244000000000.0 * sc[5] * x1158 * x269 -
          187000000000.0 * x1158 * x269 * x6 - 280500000000.0 * x1158 * x270 +
          x183 * x417 + 280500000000.0 * x379 * x450 + x417 * x460;
  x140 = 20000000.0 * sc[4] * x139 - 20000000.0 * sc[4] * x179 -
         30000000.0 * sc[4] * x210 + 30000000.0 * sc[4] * x378 +
         x1229 * x1230 * (3.0 * x1225 * x140 + 0.42 * x140 * x463) -
         1.905e+20 * x1229 * x471 -
         x1230 * x1238 * (3.0 * x1225 * x140 + 0.42 * x140 * x463) +
         1.905e+20 * x1238 * x471 + 4.4802777777777799e+34 * x1240 * x17 -
         4.4802777777777799e+34 * x17 * x468 + x209 * x472 + x375 * x417;
  x1176 = -70000000.0 * sc[0] * x214 + 70000000.0 * sc[0] * x52 + x1217 -
          x161 * x220 + x195 * x220;
  x1197 = -20000000.0 * sc[4] * x161 + 20000000.0 * sc[4] * x195 -
          30000000.0 * sc[4] * x214 + 30000000.0 * sc[4] * x52 + x10 +
          x185 * x375 + x207 * x472;
  x462 = -sc[1] * x1311 - 30000000.0 * sc[1] * x265 -
         20000000.0 * sc[1] * x323 - 1.27e+20 * sc[1] * x471 * x7 * exp(x462) -
         sc[1] * x877 - 20000000.0 * sc[4] * x197 + 30000000.0 * sc[4] * x202 +
         20000000.0 * sc[4] * x205 + x248 * x472 + x359 * x375;
  x471 = -70000000.0 * sc[0] * x282 + 70000000.0 * sc[0] * x57 + x1322 -
         x169 * x220 + x191 * x220;
  x1217 = -20000000.0 * sc[4] * x169 + 20000000.0 * sc[4] * x191 -
          30000000.0 * sc[4] * x282 + 30000000.0 * sc[4] * x57 + x119 * x375 +
          x1325 + x33 * x472;
  x1225 = ((sc[2]) * (sc[2])) * x866;
  x463 = refC * exp(-2.0 * g_RT[2] + g_RT[3] - logT);
  x1229 = sc[3] * x463;
  x1230 = exp(g_RT[10] - g_RT[2] - g_RT[9] - 1509.64999741416 * invT);
  x17 = refC * sc[10] * x1230;
  x468 = exp(-1509.64999741416 * invT);
  x1238 = sc[2] * sc[9] * x468;
  x1240 = 500000.0 * x1151 - 500000.0 * x1154 - 240000.0 * x1225 +
          240000.0 * x1229 - 602.0 * x1238 + 602.0 * x17;
  x10 = 2500000.0 *
        exp(g_RT[10] + g_RT[2] - g_RT[3] - g_RT[9] - 24053.756625465601 * invT);
  const amrex::Real x1709 =
    1020.0 * exp(-4327.6633259205901 * invT + 1.5 * logT);
  const amrex::Real x1712 = 1204.0 * sc[0] + 2107.0 * sc[10] + 3612.0 * sc[3] +
                            3612.0 * sc[5] + 903.0 * sc[9] + 602.0 * x2;
  const amrex::Real x1716 =
    2.4 * sc[0] + 3.6 * sc[10] + sc[3] + 15.4 * sc[5] + 1.75 * sc[9] + x2;
  const amrex::Real x1722 = exp(
    g_RT[2] - 2.0 * g_RT[4] + g_RT[5] + 1061.78716484796 * invT + 2.4 * logT);
  const amrex::Real x1725 = 39000000.0 * exp(-1781.38699694871 * invT);
  const amrex::Real x1730 = 2500000.0 * sc[9] * exp(-24053.756625465601 * invT);
  const amrex::Real x1732 = 20000000.0 * exp(g_RT[3] - g_RT[6] - x92);
  x92 = 1020.0 *
        exp(g_RT[7] - g_RT[8] - 4327.6633259205901 * invT + 1.5 * logT - x92);
  const amrex::Real x1750 = sc[4] * exp(1061.78716484796 * invT + 2.4 * logT);
  const amrex::Real x1772 = 2500000.0 * sc[3] * exp(-24053.756625465601 * invT);
  const amrex::Real x1795 = 13200000.0 * sc[3] * exp(-754.82499870708 * invT);
  const amrex::Real x1812 = 120000.0 * x1225 - 120000.0 * x1229;
  x97 = 1000000.0 * exp(g_RT[3] - g_RT[6] - x97);
  const amrex::Real x1826 =
    29000000.0 * exp(g_RT[3] - g_RT[6] + 251.60833290235999 * invT - x43);
  const amrex::Real x1834 = 36000.0 * exp(-4498.7569922941902 * invT);
  const amrex::Real x1847 = 100000000.0 * exp(-20128.666632188801 * invT);
  const amrex::Real x1850 =
    36000.0 *
    exp(g_RT[11] - g_RT[3] + g_RT[4] - g_RT[7] - 4498.7569922941902 * invT);
  const amrex::Real x1854 = 29000000.0 * sc[6] * exp(251.60833290235999 * invT);
  const amrex::Real x1877 = 29000000.0 * sc[4] * exp(251.60833290235999 * invT);
  x286 = 7600000.0 * exp(-201.28666632188799 * invT - x110 - x286);
  x110 =
    18000000.0 * exp(g_RT[11] - g_RT_qss[3] - 452.89499922424801 * invT - x110);
  const amrex::Real x1961 =
    56.0 *
    exp(-g_RT[7] + g_RT_qss[0] - 2727.4343286615799 * invT + 1.6 * logT - x43);
  const amrex::Real x1962 = sc[5] * x1961;
  const amrex::Real x1966 = 12000000.0 * sc[1] * x321 + x175 * x222 -
                            x175 * x327 - x1962 * x346 + x1962 * x383;
  const amrex::Real x1970 = -500000.0 * x1151 + 500000.0 * x1154 + x1347;
  const amrex::Real x1976 = 12000000.0 * sc[1] * x15 + x128 * x1962 -
                            x175 * x345 + x175 * x413 - x1962 * x400 -
                            12000000.0 * x249 - 12000000.0 * x902;
  const amrex::Real x1988 = 12000000.0 * sc[1] * x28 + sc[8] * x1709 -
                            x134 * x175 + x175 * x347 - x1962 * x420 +
                            x1962 * x428;
  const amrex::Real x1992 = 10000000.0 * sc[2] * x28 -
                            20000000.0 * sc[6] * x420 +
                            20000000.0 * sc[6] * x428 + x143 * x72 -
                            10000000.0 * x249 - 10000000.0 * x902;
  const amrex::Real x1999 = 12000000.0 * sc[1] * x389 + x175 * x351 -
                            x175 * x352 - x1962 * x434 + x1962 * x439;
  const amrex::Real x2003 =
    10000000.0 * sc[2] * x389 - 20000000.0 * sc[6] * x434 +
    20000000.0 * sc[6] * x439 + sc[7] * x1834 + x143 * x392;
  const amrex::Real x2007 =
    100.0 * sc[8] * exp(-1570.03599731073 * invT + 1.6 * logT);
  const amrex::Real x2017 =
    150000000.0 *
    exp(g_RT[10] + g_RT[4] - g_RT[6] - g_RT[9] - 11875.913312991401 * invT);
  const amrex::Real x2026 = 56.0 * exp(-2727.4343286615799 * invT + 1.6 * logT);
  const amrex::Real x2028 = 12000000.0 * sc[1] * x397 - sc[7] * x2026 -
                            12000000.0 * sc[7] * x245 - sc[7] * x92 -
                            25010000.0 * sc[7] + x103 * x1962 - x175 * x356 +
                            x175 * x84 - x1962 * x445;
  const amrex::Real x2030 =
    3430.0 * exp(224.93784961470999 * invT + 1.18 * logT);
  const amrex::Real x2035 =
    -sc[11] * x1850 - sc[11] * x2030 - 10000000.0 * sc[11] * x237 -
    20000000.0 * sc[11] * x333 + 10000000.0 * sc[2] * x397 +
    20000000.0 * sc[6] * x103 - 20000000.0 * sc[6] * x445 + x143 * x393;
  x43 =
    100.0 * exp(g_RT[7] - g_RT[8] - 1570.03599731073 * invT + 1.6 * logT - x43);
  const amrex::Real x2050 = 12000000.0 * sc[1] * x78 + x165 * x1961 -
                            x167 * x175 + x175 * x361 + x1962 * x457 -
                            x1962 * x466 + 25010000.0 * x45 * x50;
  const amrex::Real x2055 =
    10000000.0 * sc[2] * x78 + 20000000.0 * sc[6] * x457 -
    20000000.0 * sc[6] * x466 + x143 * x399 + x180 * x379;
  const amrex::Real x2059 = 150000000.0 * exp(-11875.913312991401 * invT);
  const amrex::Real x2067 = 12000000.0 * sc[1] * x232 - x175 * x25 +
                            x175 * x483 - x1962 * x477 + x1962 * x486;
  const amrex::Real x2072 =
    10000000.0 * sc[2] * x232 - 20000000.0 * sc[6] * x477 +
    20000000.0 * sc[6] * x486 + x143 * x304 + 20000000.0 * x165;
  const amrex::Real x2080 = 10000000.0 * sc[2] * x419 + sc[3] * x1834 +
                            20000000.0 * sc[6] * x339 -
                            20000000.0 * sc[6] * x489 + x143 * x402;
  const amrex::Real x2086 = 12000000.0 * sc[1] * x419 - sc[4] * x2026 -
                            12000000.0 * sc[4] * x245 - sc[4] * x92 -
                            25010000.0 * sc[4] + x175 * x494 + x1962 * x339 -
                            x1962 * x489;
  const amrex::Real x2087 =
    100.0 * sc[4] * exp(-1570.03599731073 * invT + 1.6 * logT);
  const amrex::Real x2101 =
    sc[2] * x1709 + x175 * x371 - x175 * x373 - x1962 * x500 + x1962 * x506;
  const amrex::Real x2111 = 12000000.0 * sc[1] * x209 + x139 * x1962 -
                            x175 * x210 + x175 * x378 - x179 * x1962;
  const amrex::Real x2141 =
    10000000.0 * sc[2] * x248 - sc[4] * x1850 - sc[4] * x2030 -
    10000000.0 * sc[4] * x237 - 20000000.0 * sc[4] * x333 -
    20000000.0 * sc[6] * x197 + 20000000.0 * sc[6] * x205 + x143 * x359;
  const amrex::Real x2152 = 12000000.0 * sc[1] * x33 - x169 * x1962 -
                            x175 * x282 + x175 * x57 + x191 * x1962;
  const amrex::Real x2153 =
    -22000000000.0 * sc[5] * x1163 + 22000000000.0 * x1166;
  const amrex::Real x2159 =
    -12000000.0 * sc[3] * x25 + x175 * x25 + x20 * x25 + x2153;
  const amrex::Real x2164 = -2800000.0 * x1065 * x1142 + 2800000.0 * x1147;
  const amrex::Real x2176 = -4161.6018262050302 * invT + 2.0 * logT;
  const amrex::Real x2179 = 4.92 * sc[8] * exp(x2176);
  const amrex::Real x2181 =
    -8.192e+20 * x1062 * x131 - x1082 + 6400000000000000.0 * x1092 -
    6400000000000000.0 * x1102 + x1109 + 8.192e+20 * x131 * x326;
  const amrex::Real x2184 =
    sc[1] * x1079 * x165 - 4.096e+20 * x1062 * x131 +
    x1071 * (0.28 * x1067 + 2.0 * x1069) - 3200000000000000.0 * x1102 +
    4.096e+20 * x131 * x326 - x1334 * (0.28 * x1067 + 2.0 * x1069);
  const amrex::Real x2196 =
    -70000000.0 * sc[0] * x25 + x1322 - x149 * x25 + x2184;
  const amrex::Real x2204 =
    0.00332 * exp(-2948.8496616156599 * invT + 2.81 * logT);
  const amrex::Real x2227 = 2.46 * sc[8] * exp(x2176);
  x1200 = -4.954e+21 * refC * sc[8] * x1182 * x1186 +
          4.954e+21 * sc[1] * sc[7] * x1182 * exp(x1185) +
          9.6622503937007901e+32 * x1172 * x1179 -
          9.6622503937007901e+32 * x1172 * x1192 - x1173 * x1200 +
          x1200 * x1204;
  x1316 = 4.83112519685039e+32 * x1172 * x1179 -
          4.83112519685039e+32 * x1172 * x1192 - x1173 * x1316 -
          2.477e+21 * x1201 + x1204 * x1316 + 2.477e+21 * x1210;
  const amrex::Real x2240 = -602.0 * x1238 + x1344 + 602.0 * x17 + x894;
  const amrex::Real x2258 = -14000000.0 * sc[10] * x25;
  const amrex::Real x2260 =
    -15000000.0 * sc[2] * x25 - 40000000.0 * sc[3] * x25 + x2240;
  const amrex::Real x2281 = 602.0 * x1238 - 602.0 * x17;
  x827 = 2700000000000.0 * x1063 + x1086 * x1089 - x1086 * x75 -
         3.3440366972477102e+18 * x1105 * x70 - 1.014e+16 * x314 * x342 * x815 -
         2700000000000.0 * x314 * x826 + 3.3440366972477102e+18 * x316 * x70 +
         1.014e+16 * x342 * x837 + 1.19557674418605e+30 * x54 * x822 -
         1.19557674418605e+30 * x54 * x843 + x827 * x833 - x827 * x841;
  x890 = 1350000000000.0 * x1063 + x1089 * x1328 -
         1.67201834862385e+18 * x1105 * x70 - x1328 * x75 -
         1350000000000.0 * x314 * x826 + 1.67201834862385e+18 * x316 * x70 +
         5070000000000000.0 * x342 * x837 +
         5.9778837209302301e+29 * x54 * x822 -
         5.9778837209302301e+29 * x54 * x843 - x831 + x833 * x890 - x841 * x890;
  // Write the full Jacobian expression...
  J[0] = -180000.0 * refC * sc[0] * exp(-0.6 * logT + x22) -
         550000000.0 * refC * sc[10] * exp(-2.0 * logT + x22) -
         60000000.0 * refC * sc[5] * exp(-1.25 * logT + x22) -
         70000000.0 * sc[0] * x222 + 70000000.0 * sc[0] * x327 -
         20000000.0 * sc[11] * x234 +
         90000.0 * ((sc[1]) * (sc[1])) * exp(-0.6 * logT) +
         20000000.0 * sc[1] * x321 - sc[2] * x217 - sc[3] * x200 -
         sc[7] * x188 - x117 * x379 - x165 * x53 - x186 * x375 + x220 * x346 -
         x220 * x383 - x225 - 70000000.0 * x50 -
         1000000.0 * x775 * (sc[3] + sc[9] + x2) + x845 + x853 +
         (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
          15000000.0 * sc[2] * x68 + x124 * x375) *
           (-x346 + x383) +
         (-20000000.0 * sc[1] * x41) * (x222 - x327 + (-x68) * (-x346 + x383));
  J[14] =
    180000.0 * sc[0] * sc[1] * exp(-0.6 * logT) + 70000000.0 * sc[0] * x345 -
    70000000.0 * sc[0] * x413 + 1100000000.0 * sc[10] * sc[1] * x861 +
    sc[11] * x877 + 120000000.0 * sc[1] * sc[5] * exp(-1.25 * logT) +
    20000000.0 * sc[1] * x15 + 2000000.0 * sc[1] * x866 * (sc[3] + sc[9] + x2) +
    sc[4] * x858 + sc[5] * x871 + sc[6] * x856 + 70000000.0 * sc[7] * x279 +
    sc[7] * x869 - x128 * x220 + x220 * x400 - 20000000.0 * x249 - x375 * x381 +
    x880 + x895 + x901 - 20000000.0 * x902 +
    (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
     15000000.0 * sc[2] * x68 + x124 * x375) *
      (x128 - x400) +
    (-20000000.0 * sc[1] * x41) * (-x345 + x413 + (-x68) * (x128 - x400));
  J[28] =
    70000000.0 * sc[0] * x134 - sc[0] * x217 - 70000000.0 * sc[0] * x347 +
    20000000.0 * sc[1] * x28 + 73400000.0 * sc[1] * x72 -
    15000000.0 * sc[2] * x134 + 15000000.0 * sc[2] * x347 + x220 * x420 -
    x220 * x428 - x375 * x72 + 15000000.0 * x50 + x895 +
    (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
     15000000.0 * sc[2] * x68 + x124 * x375) *
      (-x420 + x428) +
    (-20000000.0 * sc[1] * x41) * (-x134 + x347 + (-x68) * (-x420 + x428));
  J[42] = -sc[0] * x200 - 70000000.0 * sc[0] * x351 +
          70000000.0 * sc[0] * x352 + 20000000.0 * sc[1] * x389 +
          73400000.0 * sc[1] * x392 + 15000000.0 * sc[2] * x351 -
          15000000.0 * sc[2] * x352 + x220 * x434 - x220 * x439 - x375 * x392 +
          x895 +
          (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
           15000000.0 * sc[2] * x68 + x124 * x375) *
            (-x434 + x439) +
          (-20000000.0 * sc[1] * x41) * (x351 - x352 + (-x68) * (-x434 + x439));
  J[56] = 70000000.0 * sc[0] * x356 - 70000000.0 * sc[0] * x84 +
          73400000.0 * sc[1] * x393 + 20000000.0 * sc[1] * x397 + sc[1] * x858 -
          15000000.0 * sc[2] * x356 + 15000000.0 * sc[2] * x84 - x103 * x220 +
          x220 * x445 - x375 * x393 + x895 - x934 +
          (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
           15000000.0 * sc[2] * x68 + x124 * x375) *
            (x103 - x445) +
          (-20000000.0 * sc[1] * x41) * (-x356 + x84 + (-x68) * (x103 - x445));
  J[70] = -60000000.0 * refC * sc[0] * exp(-1.25 * logT + x22) +
          70000000.0 * sc[0] * x167 - 70000000.0 * sc[0] * x361 +
          60000000.0 * ((sc[1]) * (sc[1])) * exp(-1.25 * logT) +
          20000000.0 * sc[1] * x78 + sc[1] * x871 - x220 * x457 + x220 * x466 -
          x375 * x399 + x963 +
          (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
           15000000.0 * sc[2] * x68 + x124 * x375) *
            (x457 - x466) +
          (-20000000.0 * sc[1] * x41) * (-x167 + x361 + (-x68) * (x457 - x466));
  J[84] = -70000000.0 * sc[0] * x483 + 20000000.0 * sc[1] * x232 +
          73400000.0 * sc[1] * x304 + sc[1] * x856 + 15000000.0 * sc[2] * x483 +
          x220 * x477 - x220 * x486 - x304 * x375 + x973 +
          (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
           15000000.0 * sc[2] * x68 + x124 * x375) *
            (-x477 + x486) +
          (-20000000.0 * sc[1] * x41) * (-x25 + x483 + (-x68) * (-x477 + x486));
  J[98] = -sc[0] * x188 - 70000000.0 * sc[0] * x494 +
          70000000.0 * sc[1] * x279 + 73400000.0 * sc[1] * x402 +
          20000000.0 * sc[1] * x419 + sc[1] * x869 + 15000000.0 * sc[2] * x494 -
          x220 * x339 + x220 * x489 - x375 * x402 + x973 +
          (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
           15000000.0 * sc[2] * x68 + x124 * x375) *
            (x339 - x489) +
          (-20000000.0 * sc[1] * x41) * (-x25 + x494 + (-x68) * (x339 - x489));
  J[112] =
    -70000000.0 * sc[0] * x371 + 70000000.0 * sc[0] * x373 -
    2000000.0 * sc[0] * x775 + 2000000.0 * ((sc[1]) * (sc[1])) * x866 +
    20000000.0 * sc[1] * x207 + 73400000.0 * sc[1] * x81 +
    15000000.0 * sc[2] * x371 - 15000000.0 * sc[2] * x373 + x220 * x500 -
    x220 * x506 - x375 * x81 + x845 + x994 +
    (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
     15000000.0 * sc[2] * x68 + x124 * x375) *
      (-x500 + x506) +
    (-20000000.0 * sc[1] * x41) * (x371 - x373 + (-x68) * (-x500 + x506));
  J[126] =
    70000000.0 * sc[0] * x210 - 70000000.0 * sc[0] * x378 -
    1000000.0 * sc[0] * x775 + 1000000.0 * ((sc[1]) * (sc[1])) * x866 +
    20000000.0 * sc[1] * x209 - x139 * x220 + x179 * x220 + x283 - x375 * x417 +
    (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
     15000000.0 * sc[2] * x68 + x124 * x375) *
      (x139 - x179) +
    (-20000000.0 * sc[1] * x41) * (-x210 + x378 + (-x68) * (x139 - x179));
  J[140] =
    -550000000.0 * refC * sc[0] * exp(-2.0 * logT + x22) +
    70000000.0 * sc[0] * x214 - 70000000.0 * sc[0] * x52 +
    550000000.0 * ((sc[1]) * (sc[1])) * x861 + 73400000.0 * sc[1] * x185 +
    20000000.0 * sc[1] * x207 - 15000000.0 * sc[2] * x214 +
    15000000.0 * sc[2] * x52 + x161 * x220 - x185 * x375 - x195 * x220 + x845 +
    (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
     15000000.0 * sc[2] * x68 + x124 * x375) *
      (-x161 + x195) +
    (-20000000.0 * sc[1] * x41) * (-x214 + x52 + (-x68) * (-x161 + x195));
  J[154] = -70000000.0 * sc[0] * x202 - 20000000.0 * sc[0] * x234 +
           20000000.0 * sc[1] * x248 + 73400000.0 * sc[1] * x359 +
           sc[1] * x877 + 15000000.0 * sc[2] * x202 + x197 * x220 +
           5070000000000000.0 * x204 * x342 * x815 - x205 * x220 - x359 * x375 +
           x973 +
           (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
            15000000.0 * sc[2] * x68 + x124 * x375) *
             (-x197 + x205) +
           (-20000000.0 * sc[1] * x41) * (x202 - x25 + (-x68) * (-x197 + x205));
  J[168] =
    70000000.0 * sc[0] * x282 - 70000000.0 * sc[0] * x57 +
    73400000.0 * sc[1] * x119 + 20000000.0 * sc[1] * x33 -
    15000000.0 * sc[2] * x282 + 15000000.0 * sc[2] * x57 - x119 * x375 +
    x169 * x220 - x191 * x220 + x895 +
    (70000000.0 * sc[0] * x68 - 73400000.0 * sc[1] * x124 -
     15000000.0 * sc[2] * x68 + x124 * x375) *
      (-x169 + x191) +
    (-20000000.0 * sc[1] * x41) * (-x282 + x57 + (-x68) * (-x169 + x191));
  J[1] = 360000.0 * refC * sc[0] * exp(-0.6 * logT + x22) +
         1100000000.0 * refC * sc[10] * exp(-2.0 * logT + x22) +
         120000000.0 * refC * sc[5] * exp(-1.25 * logT + x22) +
         20000000.0 * sc[11] * x234 -
         180000.0 * ((sc[1]) * (sc[1])) * exp(-0.6 * logT) -
         73400000.0 * sc[1] * x186 - 38000000.0 * sc[1] * x321 +
         30000000.0 * sc[2] * x186 - 80000000.0 * sc[2] * x346 +
         80000000.0 * sc[2] * x383 + sc[3] * x200 + sc[7] * x188 +
         73400000.0 * sc[9] * x292 + 8.192e+20 * x1062 * x131 -
         2700000000000.0 * x1063 + x1082 - x1086 * x1089 + x1086 * x75 -
         6400000000000000.0 * x1092 + 6400000000000000.0 * x1102 +
         3.3440366972477102e+18 * x1105 * x70 - x1109 + x1149 + x1156 + x1160 +
         x1168 + x1218 + x1243 - 8.192e+20 * x131 * x326 + x183 * x186 -
         x186 * x443 - x186 * x444 + x186 * x460 + x20 * x222 - x20 * x327 +
         x225 + 2700000000000.0 * x314 * x826 -
         3.3440366972477102e+18 * x316 * x70 + x346 * x432 - x383 * x432 +
         2000000.0 * x775 * (sc[3] + sc[9] + x2) +
         (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
          30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
          30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
          x124 * x444 - x124 * x460 - x20 * x68) *
           (-x346 + x383) +
         (38000000.0 * sc[1] * x41 - x41 * x472) *
           (x222 - x327 + (-x68) * (-x346 + x383));
  J[15] =
    -refCinv * sc[9] * x1250 - 360000.0 * sc[0] * sc[1] * exp(-0.6 * logT) -
    2200000000.0 * sc[10] * sc[1] * x861 - sc[10] * x1267 -
    30000000.0 * sc[10] * x308 -
    240000000.0 * sc[1] * sc[5] * exp(-1.25 * logT) - 38000000.0 * sc[1] * x15 -
    73400000.0 * sc[1] * x381 -
    4000000.0 * sc[1] * x866 * (sc[3] + sc[9] + x2) +
    80000000.0 * sc[2] * x128 + 30000000.0 * sc[2] * x381 -
    80000000.0 * sc[2] * x400 - 500000.0 * sc[2] * x7 * x866 - sc[5] * x871 -
    80000000.0 * x107 * x379 - x1079 * x165 * x7 + x1149 + x1277 - x128 * x432 +
    x1293 + x1296 + x1300 + x1303 + x1315 + x1322 + x1325 + x1349 - x20 * x345 +
    x20 * x413 + 38000000.0 * x249 + 73400000.0 * x299 + 73400000.0 * x322 -
    1350000000000.0 * x379 * x438 * x7 * exp(x427) - x381 * x443 - x381 * x444 +
    x400 * x432 - x880 + 38000000.0 * x902 +
    (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
     30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
     x124 * x444 - x124 * x460 - x20 * x68) *
      (x128 - x400) +
    (38000000.0 * sc[1] * x41 - x41 * x472) *
      (-x345 + x413 + (-x68) * (x128 - x400));
  J[29] = -70000000.0 * sc[0] * x134 + sc[0] * x217 +
          70000000.0 * sc[0] * x347 - 38000000.0 * sc[1] * x28 -
          500000.0 * sc[1] * x7 * x866 - 73400000.0 * sc[1] * x72 -
          80000000.0 * sc[2] * x420 + 80000000.0 * sc[2] * x428 -
          28000000.0 * sc[3] * x134 + 28000000.0 * sc[3] * x347 +
          sc[4] * x1356 + sc[5] * x1354 + x1149 + x1322 + x1325 - x134 * x20 +
          x1349 + x1369 + x1378 + 80000000.0 * x165 + x183 * x72 + x20 * x347 -
          x220 * x420 + x220 * x428 + x420 * x432 - x428 * x432 - x443 * x72 -
          x444 * x72 + x460 * x72 +
          (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
           30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
           30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
           x124 * x444 - x124 * x460 - x20 * x68) *
            (-x420 + x428) +
          (38000000.0 * sc[1] * x41 - x41 * x472) *
            (-x134 + x347 + (-x68) * (-x420 + x428));
  J[43] = sc[0] * x200 + 70000000.0 * sc[0] * x351 - 70000000.0 * sc[0] * x352 -
          38000000.0 * sc[1] * x389 - 73400000.0 * sc[1] * x392 +
          30000000.0 * sc[2] * x392 - 80000000.0 * sc[2] * x434 +
          80000000.0 * sc[2] * x439 + 28000000.0 * sc[3] * x351 -
          28000000.0 * sc[3] * x352 + x1322 + x1325 + x1349 + x1398 + x1405 +
          x183 * x392 + x20 * x351 - x20 * x352 - x220 * x434 + x220 * x439 -
          x392 * x443 - x392 * x444 + x392 * x460 + x432 * x434 - x432 * x439 +
          28000000.0 * x50 +
          (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
           30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
           30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
           x124 * x444 - x124 * x460 - x20 * x68) *
            (-x434 + x439) +
          (38000000.0 * sc[1] * x41 - x41 * x472) *
            (x351 - x352 + (-x68) * (-x434 + x439));
  J[57] = -70000000.0 * sc[0] * x356 + 70000000.0 * sc[0] * x84 -
          73400000.0 * sc[1] * x393 - 38000000.0 * sc[1] * x397 +
          80000000.0 * sc[2] * x103 + sc[2] * x1356 +
          30000000.0 * sc[2] * x393 - 80000000.0 * sc[2] * x445 +
          12000000.0 * sc[7] * x245 + x103 * x220 - x103 * x432 + x1149 +
          500000.0 * x1150 * x7 + x1322 + x1325 + x1349 + x1407 +
          268000000.0 * x1409 + x1434 + x1436 + x183 * x393 - x20 * x356 +
          x20 * x84 - x220 * x445 - 28000000.0 * x270 * x271 - x393 * x443 -
          x393 * x444 + x393 * x460 + x432 * x445 + x934 +
          (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
           30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
           30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
           x124 * x444 - x124 * x460 - x20 * x68) *
            (x103 - x445) +
          (38000000.0 * sc[1] * x41 - x41 * x472) *
            (-x356 + x84 + (-x68) * (x103 - x445));
  J[71] =
    120000000.0 * refC * sc[0] * exp(-1.25 * logT + x22) -
    120000000.0 * ((sc[1]) * (sc[1])) * exp(-1.25 * logT) -
    sc[1] * sc[3] * x1282 - 73400000.0 * sc[1] * x399 -
    38000000.0 * sc[1] * x78 - sc[1] * x871 + sc[2] * x1354 +
    80000000.0 * sc[2] * x457 - 80000000.0 * sc[2] * x466 +
    2.4576e+21 * x1062 * x131 - 8100000000000.0 * x1063 + x1065 * x1463 -
    x1071 * (1.68 * x1067 + 12.0 * x1069) - x1089 * x1442 - 1.92e+16 * x1092 +
    1.92e+16 * x1102 + 1.0032110091743101e+19 * x1105 * x70 +
    3000000.0 * x1151 - 3000000.0 * x1154 - 2.4576e+21 * x131 * x326 +
    x1334 * (1.68 * x1067 + 12.0 * x1069) + x1442 * x75 + x1472 + x1480 +
    x1484 + x1492 - x167 * x20 + x20 * x361 + 8100000000000.0 * x314 * x826 -
    1.0032110091743101e+19 * x316 * x70 + 6000000.0 * x37 * x50 - x399 * x443 -
    x399 * x444 - x432 * x457 + x432 * x466 +
    (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
     30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
     x124 * x444 - x124 * x460 - x20 * x68) *
      (x457 - x466) +
    (38000000.0 * sc[1] * x41 - x41 * x472) *
      (-x167 + x361 + (-x68) * (x457 - x466));
  J[85] = 70000000.0 * sc[0] * x483 - 38000000.0 * sc[1] * x232 -
          73400000.0 * sc[1] * x304 + 30000000.0 * sc[2] * x304 -
          80000000.0 * sc[2] * x477 + 80000000.0 * sc[2] * x486 +
          28000000.0 * sc[3] * x483 + x1501 + x1507 + x1511 + x1516 + x1521 +
          x183 * x304 + x20 * x483 - x220 * x477 + x220 * x486 - x304 * x443 -
          x304 * x444 + x304 * x460 + x432 * x477 - x432 * x486 +
          (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
           30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
           30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
           x124 * x444 - x124 * x460 - x20 * x68) *
            (-x477 + x486) +
          (38000000.0 * sc[1] * x41 - x41 * x472) *
            (-x25 + x483 + (-x68) * (-x477 + x486));
  J[99] = sc[0] * x188 - 73400000.0 * sc[1] * x402 - 38000000.0 * sc[1] * x419 +
          80000000.0 * sc[2] * x339 + 30000000.0 * sc[2] * x402 -
          80000000.0 * sc[2] * x489 + 28000000.0 * sc[3] * x494 +
          12000000.0 * sc[4] * x245 + 2244000000000.0 * sc[5] * x402 * x450 +
          x1078 + x1511 + x1521 + x1539 + x1541 + x20 * x494 - x339 * x432 -
          x402 * x443 - x402 * x444 + 187000000000.0 * x402 * x450 * x6 +
          x432 * x489 +
          (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
           30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
           30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
           x124 * x444 - x124 * x460 - x20 * x68) *
            (x339 - x489) +
          (38000000.0 * sc[1] * x41 - x41 * x472) *
            (-x25 + x494 + (-x68) * (x339 - x489));
  J[113] = 4000000.0 * sc[0] * x775 - 4000000.0 * ((sc[1]) * (sc[1])) * x866 -
           73400000.0 * sc[1] * x81 - 80000000.0 * sc[2] * x500 +
           80000000.0 * sc[2] * x506 + 30000000.0 * sc[2] * x81 +
           28000000.0 * sc[3] * x371 - 28000000.0 * sc[3] * x373 +
           44000000000.0 * sc[5] * x1163 + x1149 - 44000000000.0 * x1166 +
           x1561 + x1563 + x1564 + x183 * x81 + x20 * x371 - x20 * x373 +
           x432 * x500 - x432 * x506 - x443 * x81 - x444 * x81 + x460 * x81 -
           x994 +
           (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
            30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
            30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
            x124 * x444 - x124 * x460 - x20 * x68) *
             (-x500 + x506) +
           (38000000.0 * sc[1] * x41 - x41 * x472) *
             (x371 - x373 + (-x68) * (-x500 + x506));
  J[127] = 73400000.0 * sc[0] * x292 + 2000000.0 * sc[0] * x775 -
           2000000.0 * ((sc[1]) * (sc[1])) * x866 - 38000000.0 * sc[1] * x209 -
           73400000.0 * sc[1] * x417 + 80000000.0 * sc[2] * x139 -
           80000000.0 * sc[2] * x179 - 28000000.0 * sc[3] * x210 +
           28000000.0 * sc[3] * x378 + 6.144e+20 * x1062 * x131 -
           2025000000000.0 * x1063 + 2100000.0 * x1065 * x1142 -
           x1071 * (0.42 * x1067 + 3.0 * x1069) - x1084 * x1089 + x1084 * x75 -
           4800000000000000.0 * x1092 + 4800000000000000.0 * x1102 +
           2.5080275229357798e+18 * x1105 * x70 - 2100000.0 * x1147 +
           750000.0 * x1151 - 750000.0 * x1154 + x1158 + x1175 - x1250 * x269 -
           6.144e+20 * x131 * x326 + x1334 * (0.42 * x1067 + 3.0 * x1069) +
           x1347 - x139 * x432 + x140 + x1565 + x179 * x432 - x20 * x210 +
           x20 * x378 + 2025000000000.0 * x314 * x826 -
           2.5080275229357798e+18 * x316 * x70 - x417 * x443 - x417 * x444 +
           (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
            30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
            30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
            x124 * x444 - x124 * x460 - x20 * x68) *
             (x139 - x179) +
           (38000000.0 * sc[1] * x41 - x41 * x472) *
             (-x210 + x378 + (-x68) * (x139 - x179));
  J[141] = 1100000000.0 * refC * sc[0] * exp(-2.0 * logT + x22) -
           1100000000.0 * ((sc[1]) * (sc[1])) * x861 - sc[1] * x1267 -
           73400000.0 * sc[1] * x185 - 30000000.0 * sc[1] * x308 -
           80000000.0 * sc[2] * x161 + 30000000.0 * sc[2] * x185 +
           80000000.0 * sc[2] * x195 - 28000000.0 * sc[3] * x214 +
           28000000.0 * sc[3] * x52 + 4200000.0 * x1065 * x1142 -
           4200000.0 * x1147 + x1176 + x1197 + x1347 + x1564 + x161 * x432 +
           x183 * x185 - x185 * x443 - x185 * x444 + x185 * x460 - x195 * x432 -
           x20 * x214 + x20 * x52 +
           (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
            30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
            30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
            x124 * x444 - x124 * x460 - x20 * x68) *
             (-x161 + x195) +
           (38000000.0 * sc[1] * x41 - x41 * x472) *
             (-x214 + x52 + (-x68) * (-x161 + x195));
  J[155] = 70000000.0 * sc[0] * x202 + 20000000.0 * sc[0] * x234 -
           38000000.0 * sc[1] * x248 - 73400000.0 * sc[1] * x359 -
           80000000.0 * sc[2] * x197 + 80000000.0 * sc[2] * x205 +
           30000000.0 * sc[2] * x359 + 28000000.0 * sc[3] * x202 +
           2244000000000.0 * sc[5] * x359 * x450 + x1511 + x1521 - x197 * x220 +
           x197 * x432 + x20 * x202 + 1350000000000.0 * x204 * x826 +
           x205 * x220 - x205 * x432 - x359 * x443 - x359 * x444 +
           187000000000.0 * x359 * x450 * x6 + x462 +
           (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
            30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
            30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
            x124 * x444 - x124 * x460 - x20 * x68) *
             (-x197 + x205) +
           (38000000.0 * sc[1] * x41 - x41 * x472) *
             (x202 - x25 + (-x68) * (-x197 + x205));
  J[169] = -375000000.0 * sc[1] * sc[3] * x1285 - 73400000.0 * sc[1] * x119 -
           38000000.0 * sc[1] * x33 + 30000000.0 * sc[2] * x119 -
           80000000.0 * sc[2] * x169 + 80000000.0 * sc[2] * x191 -
           28000000.0 * sc[3] * x282 + 28000000.0 * sc[3] * x57 +
           375000000.0 * x1065 * x1389 + x119 * x183 - x119 * x443 -
           x119 * x444 + x119 * x460 + x1217 + x1349 + x169 * x432 -
           x191 * x432 - x20 * x282 + x20 * x57 + x471 +
           (-70000000.0 * sc[0] * x68 + 73400000.0 * sc[1] * x124 -
            30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 -
            30000000.0 * sc[4] * x68 - x124 * x183 - x124 * x375 + x124 * x443 +
            x124 * x444 - x124 * x460 - x20 * x68) *
             (-x169 + x191) +
           (38000000.0 * sc[1] * x41 - x41 * x472) *
             (-x282 + x57 + (-x68) * (-x169 + x191));
  J[2] = -60000000.0 * sc[2] * x186 - sc[2] * x217 - 15000000.0 * sc[2] * x222 -
         10000000.0 * sc[2] * x321 + 15000000.0 * sc[2] * x327 +
         80000000.0 * sc[2] * x346 - 80000000.0 * sc[2] * x383 +
         15000000.0 * sc[9] * x256 + x1156 - 576000.0 * x1225 +
         576000.0 * x1229 - 1204.0 * x1238 + 1204.0 * x17 + x186 * x444 +
         x186 * x476 +
         (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
          x124 * x476) *
           (-x346 + x383) +
         (10000000.0 * sc[2] * x41) * (x222 - x327 + (-x68) * (-x346 + x383));
  J[16] = 30000000.0 * sc[10] * x308 + sc[11] * x1311 -
          80000000.0 * sc[2] * x128 - 10000000.0 * sc[2] * x15 +
          15000000.0 * sc[2] * x345 - 60000000.0 * sc[2] * x381 +
          80000000.0 * sc[2] * x400 - 15000000.0 * sc[2] * x413 -
          500000.0 * sc[2] * x7 * x866 + sc[3] * x1291 + sc[4] * x858 +
          sc[6] * x1275 + 80000000.0 * x107 * x379 + x1240 + x381 * x444 +
          x381 * x476 +
          (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
           x124 * x476) *
            (x128 - x400) +
          (10000000.0 * sc[2] * x41) * (-x345 + x413 + (-x68) * (x128 - x400));
  J[30] = -sc[0] * x217 - sc[10] * x10 - sc[11] * x1725 -
          500000.0 * sc[1] * x7 * x866 + 15000000.0 * sc[2] * x134 -
          480000.0 * sc[2] * x1716 * x866 - 10000000.0 * sc[2] * x28 -
          15000000.0 * sc[2] * x347 + 80000000.0 * sc[2] * x420 -
          80000000.0 * sc[2] * x428 - 60000000.0 * sc[2] * x72 - sc[4] * x1356 -
          sc[5] * x1354 - 0.0357 * sc[5] * x1722 - 20000000.0 * sc[6] -
          84300000.0 * sc[7] - sc[8] * x1709 - sc[9] * x1712 * x468 + x1240 -
          80000000.0 * x165 + 10000000.0 * x249 + 60000000.0 * x299 +
          60000000.0 * x322 + x444 * x72 + x476 * x72 - 15000000.0 * x50 +
          10000000.0 * x902 +
          (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
           x124 * x476) *
            (-x420 + x428) +
          (10000000.0 * sc[2] * x41) * (-x134 + x347 + (-x68) * (-x420 + x428));
  J[44] = sc[1] * x1291 - 15000000.0 * sc[2] * x351 +
          15000000.0 * sc[2] * x352 - 10000000.0 * sc[2] * x389 -
          60000000.0 * sc[2] * x392 + 80000000.0 * sc[2] * x434 -
          80000000.0 * sc[2] * x439 + sc[4] * x1732 + 500000.0 * x1151 -
          500000.0 * x1154 - 240000.0 * x1225 + 240000.0 * x1229 -
          3612.0 * x1238 + 3612.0 * x17 + 240000.0 * x1716 * x463 + x1730 +
          x392 * x444 + x392 * x476 +
          (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
           x124 * x476) *
            (-x434 + x439) +
          (10000000.0 * sc[2] * x41) * (x351 - x352 + (-x68) * (-x434 + x439));
  J[58] = 10000000.0 * sc[11] * x237 + sc[1] * x858 -
          80000000.0 * sc[2] * x103 - sc[2] * x1356 +
          15000000.0 * sc[2] * x356 - 60000000.0 * sc[2] * x393 -
          10000000.0 * sc[2] * x397 + 80000000.0 * sc[2] * x445 -
          15000000.0 * sc[2] * x84 + sc[3] * x1732 + sc[7] * x92 +
          30000000.0 * sc[9] * x295 + 500000.0 * x1150 * x7 + x1240 +
          0.0714 * x1750 + x379 * x475 + x393 * x444 + x393 * x476 +
          (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
           x124 * x476) *
            (x103 - x445) +
          (10000000.0 * sc[2] * x41) * (-x356 + x84 + (-x68) * (x103 - x445));
  J[72] = -sc[2] * x1354 + 15000000.0 * sc[2] * x167 - 0.0357 * sc[2] * x1722 -
          15000000.0 * sc[2] * x361 - 60000000.0 * sc[2] * x399 -
          80000000.0 * sc[2] * x457 + 80000000.0 * sc[2] * x466 -
          10000000.0 * sc[2] * x78 + 3000000.0 * x1151 - 3000000.0 * x1154 -
          3696000.0 * x1225 + 3696000.0 * x1229 - 3612.0 * x1238 +
          3612.0 * x17 + x399 * x444 + x399 * x476 +
          (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
           x124 * x476) *
            (x457 - x466) +
          (10000000.0 * sc[2] * x41) * (-x167 + x361 + (-x68) * (x457 - x466));
  J[86] = sc[1] * x1275 - 10000000.0 * sc[2] * x232 + 15000000.0 * sc[2] * x25 -
          60000000.0 * sc[2] * x304 + 80000000.0 * sc[2] * x477 -
          15000000.0 * sc[2] * x483 - 80000000.0 * sc[2] * x486 -
          20000000.0 * sc[2] + x1240 + x304 * x444 + x304 * x476 +
          (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
           x124 * x476) *
            (-x477 + x486) +
          (10000000.0 * sc[2] * x41) * (-x25 + x483 + (-x68) * (-x477 + x486));
  J[100] = 15000000.0 * sc[2] * x25 - 80000000.0 * sc[2] * x339 -
           60000000.0 * sc[2] * x402 - 10000000.0 * sc[2] * x419 +
           80000000.0 * sc[2] * x489 - 15000000.0 * sc[2] * x494 -
           84300000.0 * sc[2] + sc[4] * x92 + x1240 + x402 * x444 +
           x402 * x476 +
           (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
            x124 * x476) *
             (x339 - x489) +
           (10000000.0 * sc[2] * x41) * (-x25 + x494 + (-x68) * (x339 - x489));
  J[114] = -sc[2] * x1709 - 10000000.0 * sc[2] * x207 -
           15000000.0 * sc[2] * x371 + 15000000.0 * sc[2] * x373 +
           80000000.0 * sc[2] * x500 - 80000000.0 * sc[2] * x506 -
           60000000.0 * sc[2] * x81 + x1156 - 480000.0 * x1225 +
           480000.0 * x1229 - 1204.0 * x1238 + 1204.0 * x17 + x444 * x81 +
           x476 * x81 +
           (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
            x124 * x476) *
             (-x500 + x506) +
           (10000000.0 * sc[2] * x41) * (x371 - x373 + (-x68) * (-x500 + x506));
  J[128] = 15000000.0 * sc[0] * x256 - 80000000.0 * sc[2] * x139 -
           sc[2] * x1712 * x468 + 80000000.0 * sc[2] * x179 -
           10000000.0 * sc[2] * x209 + 15000000.0 * sc[2] * x210 -
           15000000.0 * sc[2] * x378 - 60000000.0 * sc[2] * x417 +
           30000000.0 * sc[4] * x295 + 750000.0 * x1151 - 750000.0 * x1154 -
           420000.0 * x1225 + 420000.0 * x1229 - 903.0 * x1238 + 903.0 * x17 +
           x1772 + x417 * x444 + x417 * x476 +
           (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
            x124 * x476) *
             (x139 - x179) +
           (10000000.0 * sc[2] * x41) * (-x210 + x378 + (-x68) * (x139 - x179));
  J[142] = refC * x1230 * x1712 + 30000000.0 * sc[1] * x308 - sc[2] * x10 +
           80000000.0 * sc[2] * x161 - 60000000.0 * sc[2] * x185 -
           80000000.0 * sc[2] * x195 - 10000000.0 * sc[2] * x207 +
           15000000.0 * sc[2] * x214 - 15000000.0 * sc[2] * x52 + x1156 -
           864000.0 * x1225 + 864000.0 * x1229 - 2107.0 * x1238 + 2107.0 * x17 +
           x185 * x444 + x185 * x476 +
           (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
            x124 * x476) *
             (-x161 + x195) +
           (10000000.0 * sc[2] * x41) * (-x214 + x52 + (-x68) * (-x161 + x195));
  J[156] = sc[1] * x1311 - sc[2] * x1725 + 80000000.0 * sc[2] * x197 -
           15000000.0 * sc[2] * x202 - 80000000.0 * sc[2] * x205 -
           10000000.0 * sc[2] * x248 + 15000000.0 * sc[2] * x25 -
           60000000.0 * sc[2] * x359 + 10000000.0 * sc[4] * x237 + x1240 +
           x359 * x444 + x359 * x476 +
           (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
            x124 * x476) *
             (-x197 + x205) +
           (10000000.0 * sc[2] * x41) * (x202 - x25 + (-x68) * (-x197 + x205));
  J[170] = -60000000.0 * sc[2] * x119 + 80000000.0 * sc[2] * x169 -
           80000000.0 * sc[2] * x191 + 15000000.0 * sc[2] * x282 -
           10000000.0 * sc[2] * x33 - 15000000.0 * sc[2] * x57 + x119 * x444 +
           x119 * x476 + x1240 +
           (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x68 - x124 * x444 -
            x124 * x476) *
             (-x169 + x191) +
           (10000000.0 * sc[2] * x41) * (-x282 + x57 + (-x68) * (-x169 + x191));
  J[3] = -sc[3] * x200 - 40000000.0 * sc[3] * x222 + 40000000.0 * sc[3] * x327 +
         x1149 + 288000.0 * x1225 - 288000.0 * x1229 + x154 * x186 +
         x155 * x186 - x174 * x321 + x1795 * x346 - x1795 * x383 - x186 * x497 +
         (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
           (-x346 + x383) +
         (x174 * x41) * (x222 - x327 + (-x68) * (-x346 + x383));
  J[17] = refCinv * sc[9] * x1250 + 40000000.0 * sc[3] * x345 -
          40000000.0 * sc[3] * x413 + sc[6] * x856 + x1149 - x128 * x1795 +
          x1293 - x15 * x174 + x154 * x381 + x155 * x381 + x1795 * x400 +
          x1812 - x381 * x497 +
          (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
            (x128 - x400) +
          (x174 * x41) * (-x345 + x413 + (-x68) * (x128 - x400));
  J[31] = sc[10] * x10 + 240000.0 * sc[2] * x1716 * x866 +
          40000000.0 * sc[3] * x134 - 40000000.0 * sc[3] * x347 +
          sc[4] * x1356 + 20000000.0 * sc[6] + x1149 + x154 * x72 + x155 * x72 -
          x174 * x28 + x1795 * x420 - x1795 * x428 + x1812 - x497 * x72 +
          (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
            (-x420 + x428) +
          (x174 * x41) * (-x134 + x347 + (-x68) * (-x420 + x428));
  J[45] =
    -sc[0] * x200 - sc[11] * x1847 - 40000000.0 * sc[3] * x351 +
    40000000.0 * sc[3] * x352 - sc[4] * x1732 - sc[5] * x1826 - sc[7] * x1834 -
    sc[8] * x97 - 18000000.0 * x11 * exp(-452.89499922424801 * invT) + x1398 +
    x154 * x392 + x155 * x392 -
    13200000.0 * x165 * exp(-754.82499870708 * invT) - 120000.0 * x1716 * x463 -
    x1730 - x174 * x389 + x1795 * x434 - x1795 * x439 + x1812 -
    7600000.0 * x379 * exp(-201.28666632188799 * invT) - x392 * x497 -
    40000000.0 * x50 +
    (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
      (-x434 + x439) +
    (x174 * x41) * (x351 - x352 + (-x68) * (-x434 + x439));
  J[59] = sc[11] * x1850 + sc[2] * x1356 - sc[3] * x1732 +
          40000000.0 * sc[3] * x356 - 40000000.0 * sc[3] * x84 - x103 * x1795 +
          x1149 + x154 * x393 + x155 * x393 - x174 * x397 + x1795 * x445 +
          x1812 + x1854 + 28000000.0 * x270 * x271 + x379 * x504 - x393 * x497 +
          (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
            (x103 - x445) +
          (x174 * x41) * (-x356 + x84 + (-x68) * (x103 - x445));
  J[73] = -sc[1] * sc[3] * x1282 + 40000000.0 * sc[3] * x167 - sc[3] * x1826 -
          40000000.0 * sc[3] * x361 + 12000000.0 * sc[9] * x253 +
          x1065 * x1463 + 1848000.0 * x1225 - 1848000.0 * x1229 + x154 * x399 +
          x155 * x399 - x174 * x78 - x1795 * x457 + x1795 * x466 - x399 * x497 +
          (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
            (x457 - x466) +
          (x174 * x41) * (-x167 + x361 + (-x68) * (x457 - x466));
  J[87] = sc[11] * x110 + sc[1] * x856 + 20000000.0 * sc[2] +
          40000000.0 * sc[3] * x25 - 40000000.0 * sc[3] * x483 +
          1000000.0 * sc[7] + sc[9] * x286 + x1149 + x1507 + x154 * x304 +
          x155 * x304 - x174 * x232 + x1795 * x477 - x1795 * x486 + x1812 +
          x1877 + x30 * x379 - x304 * x497 +
          (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
            (-x477 + x486) +
          (x174 * x41) * (-x25 + x483 + (-x68) * (-x477 + x486));
  J[101] =
    -sc[3] * x1834 + 40000000.0 * sc[3] * x25 - 40000000.0 * sc[3] * x494 +
    1000000.0 * sc[6] + x1149 + x154 * x402 + x155 * x402 - x174 * x419 -
    x1795 * x339 + x1795 * x489 + x1812 - x402 * x497 +
    (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
      (x339 - x489) +
    (x174 * x41) * (-x25 + x494 + (-x68) * (x339 - x489));
  J[115] =
    -40000000.0 * sc[3] * x371 + 40000000.0 * sc[3] * x373 - sc[3] * x97 +
    x1149 + 240000.0 * x1225 - 240000.0 * x1229 + x154 * x81 + x155 * x81 -
    x174 * x207 + x1795 * x500 - x1795 * x506 - x497 * x81 +
    (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
      (-x500 + x506) +
    (x174 * x41) * (x371 - x373 + (-x68) * (-x500 + x506));
  J[129] =
    40000000.0 * sc[3] * x210 - 40000000.0 * sc[3] * x378 +
    12000000.0 * sc[5] * x253 + sc[6] * x286 + 2100000.0 * x1065 * x1142 -
    2100000.0 * x1147 + 210000.0 * x1225 - 210000.0 * x1229 + x1250 * x269 -
    x139 * x1795 + x154 * x417 + x155 * x417 - x174 * x209 - x1772 +
    x179 * x1795 - x417 * x497 +
    (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
      (x139 - x179) +
    (x174 * x41) * (-x210 + x378 + (-x68) * (x139 - x179));
  J[143] =
    sc[2] * x10 + 40000000.0 * sc[3] * x214 - 40000000.0 * sc[3] * x52 +
    4200000.0 * x1065 * x1142 - 4200000.0 * x1147 + 432000.0 * x1225 -
    432000.0 * x1229 + x154 * x185 + x155 * x185 + x161 * x1795 - x174 * x207 -
    x1795 * x195 - x185 * x497 +
    (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
      (-x161 + x195) +
    (x174 * x41) * (-x214 + x52 + (-x68) * (-x161 + x195));
  J[157] =
    -sc[3] * x1847 - 40000000.0 * sc[3] * x202 + 40000000.0 * sc[3] * x25 +
    sc[4] * x1850 + sc[6] * x110 + x1149 + x154 * x359 + x155 * x359 -
    x174 * x248 + x1795 * x197 - x1795 * x205 + x1812 - x359 * x497 +
    (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
      (-x197 + x205) +
    (x174 * x41) * (x202 - x25 + (-x68) * (-x197 + x205));
  J[171] =
    -375000000.0 * sc[1] * sc[3] * x1285 + 40000000.0 * sc[3] * x282 -
    40000000.0 * sc[3] * x57 + 375000000.0 * x1065 * x1389 + x119 * x154 +
    x119 * x155 - x119 * x497 + x169 * x1795 - x174 * x33 - x1795 * x191 +
    x1812 +
    (40000000.0 * sc[3] * x68 - x124 * x154 - x124 * x155 + x124 * x497) *
      (-x169 + x191) +
    (x174 * x41) * (-x282 + x57 + (-x68) * (-x169 + x191));
  J[4] = 30000000.0 * sc[2] * x186 + 10000000.0 * sc[2] * x321 -
         50000000.0 * sc[4] * x186 - 30000000.0 * sc[4] * x222 -
         5000000.0 * sc[4] * x321 + 30000000.0 * sc[4] * x327 +
         20000000.0 * sc[4] * x346 - 20000000.0 * sc[4] * x383 -
         20000000.0 * sc[6] * x346 + 20000000.0 * sc[6] * x383 -
         1000000.0 * x1151 + 1000000.0 * x1154 + x1168 + x143 * x186 -
         x154 * x186 - x1795 * x346 + x1795 * x383 - x186 * x476 + x1966 -
         x225 +
         (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
          50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
          x124 * x154 + x124 * x476 - x175 * x68) *
           (-x346 + x383) +
         (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
          5000000.0 * sc[4] * x41) *
           (x222 - x327 + (-x68) * (-x346 + x383));
  J[18] = -refCinv * sc[9] * x1250 + sc[10] * x1267 +
          30000000.0 * sc[11] * x265 + 20000000.0 * sc[11] * x323 +
          10000000.0 * sc[2] * x15 + 30000000.0 * sc[2] * x381 +
          500000.0 * sc[2] * x7 * x866 + sc[3] * x1291 -
          20000000.0 * sc[4] * x128 - 5000000.0 * sc[4] * x15 +
          30000000.0 * sc[4] * x345 - 50000000.0 * sc[4] * x381 +
          20000000.0 * sc[4] * x400 - 30000000.0 * sc[4] * x413 + sc[5] * x871 +
          20000000.0 * sc[6] * x128 - 20000000.0 * sc[6] * x400 +
          268000000.0 * sc[6] * exp(-319.54258278599701 * invT) + x128 * x1795 +
          x1303 + x143 * x381 - x154 * x381 - x1795 * x400 + x1970 + x1976 -
          x381 * x476 +
          (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
           50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
           x124 * x154 + x124 * x476 - x175 * x68) *
            (x128 - x400) +
          (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
           5000000.0 * sc[4] * x41) *
            (-x345 + x413 + (-x68) * (x128 - x400));
  J[32] = sc[0] * x217 + sc[11] * x1725 + 500000.0 * sc[1] * x7 * x866 -
          28000000.0 * sc[3] * x134 + 28000000.0 * sc[3] * x347 +
          30000000.0 * sc[4] * x134 - sc[4] * x1356 - 5000000.0 * sc[4] * x28 -
          30000000.0 * sc[4] * x347 + 20000000.0 * sc[4] * x420 -
          20000000.0 * sc[4] * x428 - 50000000.0 * sc[4] * x72 +
          0.0714 * sc[5] * x1722 + 20000000.0 * sc[6] + x1369 - x154 * x72 -
          x1795 * x420 + x1795 * x428 + x1970 + x1988 + x1992 - x476 * x72 +
          (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
           50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
           x124 * x154 + x124 * x476 - x175 * x68) *
            (-x420 + x428) +
          (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
           5000000.0 * sc[4] * x41) *
            (-x134 + x347 + (-x68) * (-x420 + x428));
  J[46] =
    sc[1] * x1291 + 30000000.0 * sc[2] * x392 + 28000000.0 * sc[3] * x351 -
    28000000.0 * sc[3] * x352 - sc[4] * x1732 - 30000000.0 * sc[4] * x351 +
    30000000.0 * sc[4] * x352 - 5000000.0 * sc[4] * x389 -
    50000000.0 * sc[4] * x392 + 20000000.0 * sc[4] * x434 -
    20000000.0 * sc[4] * x439 + sc[5] * x1826 - x154 * x392 +
    13200000.0 * x165 * exp(-754.82499870708 * invT) - x1795 * x434 +
    x1795 * x439 + x1970 + x1999 + x2003 - x392 * x476 + 28000000.0 * x50 +
    (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
     50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
     x124 * x154 + x124 * x476 - x175 * x68) *
      (-x434 + x439) +
    (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
     5000000.0 * sc[4] * x41) *
      (x351 - x352 + (-x68) * (-x434 + x439));
  J[60] = -sc[10] * x2017 - sc[2] * x1356 + 30000000.0 * sc[2] * x393 -
          sc[3] * x1732 - 20000000.0 * sc[4] * x103 +
          30000000.0 * sc[4] * x356 - 50000000.0 * sc[4] * x393 -
          5000000.0 * sc[4] * x397 + 20000000.0 * sc[4] * x445 -
          30000000.0 * sc[4] * x84 - 30000000.0 * sc[9] * x295 + x103 * x1795 -
          500000.0 * x1150 * x7 - x1407 - 536000000.0 * x1409 + x1436 -
          x154 * x393 - 20000000.0 * x165 - 0.1428 * x1750 - x1795 * x445 -
          x1854 + x1970 - x2007 + x2028 + x2035 + 5000000.0 * x249 -
          28000000.0 * x270 * x271 + 50000000.0 * x299 + 50000000.0 * x322 -
          x379 * x475 - x379 * x504 - x393 * x476 - 30000000.0 * x50 +
          5000000.0 * x902 - x934 +
          (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
           50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
           x124 * x154 + x124 * x476 - x175 * x68) *
            (x103 - x445) +
          (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
           5000000.0 * sc[4] * x41) *
            (-x356 + x84 + (-x68) * (x103 - x445));
  J[74] = 5000000.0 * sc[11] * x228 + sc[1] * x871 + 0.0714 * sc[2] * x1722 +
          30000000.0 * sc[2] * x399 + sc[3] * x1826 +
          30000000.0 * sc[4] * x167 - 30000000.0 * sc[4] * x361 -
          50000000.0 * sc[4] * x399 - 20000000.0 * sc[4] * x457 +
          20000000.0 * sc[4] * x466 - 5000000.0 * sc[4] * x78 + sc[7] * x43 +
          50000000.0 * sc[9] * x287 - 3000000.0 * x1151 + 3000000.0 * x1154 +
          x1484 - x154 * x399 + x1795 * x457 - x1795 * x466 + x2050 + x2055 -
          x399 * x476 +
          (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
           50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
           x124 * x154 + x124 * x476 - x175 * x68) *
            (x457 - x466) +
          (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
           5000000.0 * sc[4] * x41) *
            (-x167 + x361 + (-x68) * (x457 - x466));
  J[88] = 268000000.0 * sc[1] * exp(-319.54258278599701 * invT) +
          30000000.0 * sc[2] * x304 + 20000000.0 * sc[2] +
          28000000.0 * sc[3] * x483 - 5000000.0 * sc[4] * x232 +
          30000000.0 * sc[4] * x25 - 50000000.0 * sc[4] * x304 +
          20000000.0 * sc[4] * x477 - 30000000.0 * sc[4] * x483 -
          20000000.0 * sc[4] * x486 + sc[9] * x2059 + x1511 - x154 * x304 -
          x1795 * x477 + x1795 * x486 - x1877 + x1970 + x2067 + x2072 -
          x304 * x476 +
          (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
           50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
           x124 * x154 + x124 * x476 - x175 * x68) *
            (-x477 + x486) +
          (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
           5000000.0 * sc[4] * x41) *
            (-x25 + x483 + (-x68) * (-x477 + x486));
  J[102] = 30000000.0 * sc[2] * x402 + 28000000.0 * sc[3] * x494 +
           30000000.0 * sc[4] * x25 - 20000000.0 * sc[4] * x339 -
           50000000.0 * sc[4] * x402 - 5000000.0 * sc[4] * x419 +
           20000000.0 * sc[4] * x489 - 30000000.0 * sc[4] * x494 + sc[5] * x43 +
           x1511 - x154 * x402 - x175 * x25 + x1795 * x339 - x1795 * x489 +
           x1970 + x2080 + x2086 - x402 * x476 +
           (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
            50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
            x124 * x154 + x124 * x476 - x175 * x68) *
             (x339 - x489) +
           (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
            5000000.0 * sc[4] * x41) *
             (-x25 + x494 + (-x68) * (x339 - x489));
  J[116] = 12000000.0 * sc[1] * x207 + 10000000.0 * sc[2] * x207 +
           30000000.0 * sc[2] * x81 + 28000000.0 * sc[3] * x371 -
           28000000.0 * sc[3] * x373 - 5000000.0 * sc[4] * x207 -
           30000000.0 * sc[4] * x371 + 30000000.0 * sc[4] * x373 +
           20000000.0 * sc[4] * x500 - 20000000.0 * sc[4] * x506 -
           50000000.0 * sc[4] * x81 + 44000000000.0 * sc[5] * x1163 -
           20000000.0 * sc[6] * x500 + 20000000.0 * sc[6] * x506 -
           1000000.0 * x1151 + 1000000.0 * x1154 - 44000000000.0 * x1166 +
           x143 * x81 - x154 * x81 - x1795 * x500 + x1795 * x506 - x2087 +
           x2101 - x476 * x81 +
           (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
            50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
            x124 * x154 + x124 * x476 - x175 * x68) *
             (-x500 + x506) +
           (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
            5000000.0 * sc[4] * x41) *
             (x371 - x373 + (-x68) * (-x500 + x506));
  J[130] = 10000000.0 * sc[2] * x209 + 30000000.0 * sc[2] * x417 -
           28000000.0 * sc[3] * x210 + 28000000.0 * sc[3] * x378 -
           20000000.0 * sc[4] * x139 + 20000000.0 * sc[4] * x179 -
           5000000.0 * sc[4] * x209 + 30000000.0 * sc[4] * x210 -
           30000000.0 * sc[4] * x295 - 30000000.0 * sc[4] * x378 -
           50000000.0 * sc[4] * x417 + 50000000.0 * sc[5] * x287 +
           20000000.0 * sc[6] * x139 - 20000000.0 * sc[6] * x179 +
           sc[6] * x2059 - 750000.0 * x1151 + 750000.0 * x1154 - x1250 * x269 +
           x1347 + x139 * x1795 + x143 * x417 - x154 * x417 - x1565 -
           x179 * x1795 + x2111 - x417 * x476 +
           (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
            50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
            x124 * x154 + x124 * x476 - x175 * x68) *
             (x139 - x179) +
           (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
            5000000.0 * sc[4] * x41) *
             (-x210 + x378 + (-x68) * (x139 - x179));
  J[144] =
    sc[1] * x1267 + 12000000.0 * sc[1] * x207 + 30000000.0 * sc[2] * x185 +
    10000000.0 * sc[2] * x207 - 28000000.0 * sc[3] * x214 +
    28000000.0 * sc[3] * x52 + 20000000.0 * sc[4] * x161 -
    50000000.0 * sc[4] * x185 - 20000000.0 * sc[4] * x195 - sc[4] * x2017 -
    5000000.0 * sc[4] * x207 + 30000000.0 * sc[4] * x214 -
    30000000.0 * sc[4] * x52 - 20000000.0 * sc[6] * x161 +
    20000000.0 * sc[6] * x195 - 1000000.0 * x1151 + 1000000.0 * x1154 + x1347 +
    x143 * x185 - x154 * x185 - x161 * x1795 - x161 * x1962 - x175 * x214 +
    x175 * x52 + x1795 * x195 - x185 * x476 + x195 * x1962 +
    (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
     50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
     x124 * x154 + x124 * x476 - x175 * x68) *
      (-x161 + x195) +
    (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
     5000000.0 * sc[4] * x41) *
      (-x214 + x52 + (-x68) * (-x161 + x195));
  J[158] = 12000000.0 * sc[1] * x248 + 30000000.0 * sc[1] * x265 +
           20000000.0 * sc[1] * x323 + sc[2] * x1725 +
           30000000.0 * sc[2] * x359 + 28000000.0 * sc[3] * x202 +
           20000000.0 * sc[4] * x197 - 30000000.0 * sc[4] * x202 -
           20000000.0 * sc[4] * x205 - 5000000.0 * sc[4] * x248 +
           30000000.0 * sc[4] * x25 - 50000000.0 * sc[4] * x359 +
           5000000.0 * sc[5] * x228 + x1511 - x154 * x359 + x175 * x202 -
           x175 * x25 - x1795 * x197 + x1795 * x205 - x1962 * x197 +
           x1962 * x205 + x1970 + x2141 - x359 * x476 +
           (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
            50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
            x124 * x154 + x124 * x476 - x175 * x68) *
             (-x197 + x205) +
           (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
            5000000.0 * sc[4] * x41) *
             (x202 - x25 + (-x68) * (-x197 + x205));
  J[172] = 30000000.0 * sc[2] * x119 + 10000000.0 * sc[2] * x33 -
           28000000.0 * sc[3] * x282 + 28000000.0 * sc[3] * x57 -
           50000000.0 * sc[4] * x119 + 20000000.0 * sc[4] * x169 -
           20000000.0 * sc[4] * x191 + 30000000.0 * sc[4] * x282 -
           5000000.0 * sc[4] * x33 - 30000000.0 * sc[4] * x57 -
           20000000.0 * sc[6] * x169 + 20000000.0 * sc[6] * x191 + x119 * x143 -
           x119 * x154 - x119 * x476 - x169 * x1795 + x1795 * x191 + x1970 +
           x2152 +
           (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x68 +
            50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x68 - x124 * x143 +
            x124 * x154 + x124 * x476 - x175 * x68) *
             (-x169 + x191) +
           (-12000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 +
            5000000.0 * sc[4] * x41) *
             (-x282 + x57 + (-x68) * (-x169 + x191));
  J[5] = 6000000.0 * sc[1] * x321 + 12000000.0 * sc[3] * x222 -
         12000000.0 * sc[3] * x327 + 50000000.0 * sc[4] * x186 +
         5000000.0 * sc[4] * x321 - 16060000000.0 * sc[5] * x1163 +
         16060000000.0 * x1166 - x143 * x186 - x175 * x222 + x175 * x327 +
         x1962 * x346 - x1962 * x383 - x20 * x222 + x20 * x327 + x225 +
         (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 + x124 * x143 +
          x175 * x68 + x20 * x68) *
           (-x346 + x383) +
         (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
           (x222 - x327 + (-x68) * (-x346 + x383));
  J[19] = 6000000.0 * sc[1] * x15 - 12000000.0 * sc[3] * x345 +
          12000000.0 * sc[3] * x413 + sc[4] * x1301 * x861 +
          5000000.0 * sc[4] * x15 + 50000000.0 * sc[4] * x381 - sc[5] * x871 +
          sc[6] * x1275 - x128 * x1962 - x143 * x381 + x175 * x345 -
          x175 * x413 + x1962 * x400 + x20 * x345 - x20 * x413 + x2153 -
          6000000.0 * x249 - 6000000.0 * x902 +
          (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 + x124 * x143 +
           x175 * x68 + x20 * x68) *
            (x128 - x400) +
          (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
            (-x345 + x413 + (-x68) * (x128 - x400));
  J[33] = 6000000.0 * sc[1] * x28 - 12000000.0 * sc[3] * x134 +
          12000000.0 * sc[3] * x347 + 5000000.0 * sc[4] * x28 +
          50000000.0 * sc[4] * x72 - sc[5] * x1354 - 0.0357 * sc[5] * x1722 +
          x134 * x175 + x134 * x20 - x143 * x72 - x175 * x347 + x1962 * x420 -
          x1962 * x428 - x20 * x347 + x2153 +
          (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 + x124 * x143 +
           x175 * x68 + x20 * x68) *
            (-x420 + x428) +
          (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
            (-x134 + x347 + (-x68) * (-x420 + x428));
  J[47] = 6000000.0 * sc[1] * x389 + 12000000.0 * sc[3] * x351 -
          12000000.0 * sc[3] * x352 + 5000000.0 * sc[4] * x389 +
          50000000.0 * sc[4] * x392 - sc[5] * x1826 - x143 * x392 -
          x175 * x351 + x175 * x352 + x1962 * x434 - x1962 * x439 - x20 * x351 +
          x20 * x352 + x2153 + 12000000.0 * x50 +
          (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 + x124 * x143 +
           x175 * x68 + x20 * x68) *
            (-x434 + x439) +
          (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
            (x351 - x352 + (-x68) * (-x434 + x439));
  J[61] = sc[11] * x2030 + sc[1] * x1301 * x861 + 6000000.0 * sc[1] * x397 -
          12000000.0 * sc[3] * x356 + 12000000.0 * sc[3] * x84 +
          50000000.0 * sc[4] * x393 + 5000000.0 * sc[4] * x397 + sc[7] * x2026 +
          25010000.0 * sc[7] - x103 * x1962 - x143 * x393 + x175 * x356 -
          x175 * x84 + 0.0714 * x1750 + x1854 + x1962 * x445 + x20 * x356 -
          x20 * x84 + x2007 + x2153 - 5000000.0 * x249 - 50000000.0 * x299 -
          50000000.0 * x322 - 5000000.0 * x902 + x934 +
          (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 + x124 * x143 +
           x175 * x68 + x20 * x68) *
            (x103 - x445) +
          (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
            (-x356 + x84 + (-x68) * (x103 - x445));
  J[75] =
    -5000000.0 * sc[11] * x228 + 6000000.0 * sc[1] * x78 - sc[1] * x871 -
    sc[2] * x1354 - 0.0357 * sc[2] * x1722 - 12000000.0 * sc[3] * x167 -
    sc[3] * x1826 + 12000000.0 * sc[3] * x361 + 50000000.0 * sc[4] * x399 +
    5000000.0 * sc[4] * x78 - 80300000000.0 * sc[5] * x1163 - sc[7] * x43 -
    12000000.0 * sc[9] * x253 - 50000000.0 * sc[9] * x287 - x1163 * x1301 +
    80300000000.0 * x1166 - x143 * x399 - x165 * x1961 + x167 * x175 +
    x167 * x20 - x175 * x361 - x180 * x379 - x1962 * x457 + x1962 * x466 -
    x20 * x361 - 6000000.0 * x37 * x50 - 25010000.0 * x45 * x50 +
    (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 + x124 * x143 +
     x175 * x68 + x20 * x68) *
      (x457 - x466) +
    (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
      (-x167 + x361 + (-x68) * (x457 - x466));
  J[89] = sc[1] * x1275 + 6000000.0 * sc[1] * x232 + 12000000.0 * sc[3] * x483 +
          5000000.0 * sc[4] * x232 + 50000000.0 * sc[4] * x304 - x143 * x304 -
          x175 * x483 + x1877 + x1962 * x477 - x1962 * x486 - x20 * x483 +
          x2159 +
          (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 + x124 * x143 +
           x175 * x68 + x20 * x68) *
            (-x477 + x486) +
          (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
            (-x25 + x483 + (-x68) * (-x477 + x486));
  J[103] = 6000000.0 * sc[1] * x419 + 12000000.0 * sc[3] * x494 +
           sc[4] * x2026 + 50000000.0 * sc[4] * x402 +
           5000000.0 * sc[4] * x419 + 25010000.0 * sc[4] - sc[5] * x43 -
           x143 * x402 - x175 * x494 - x1962 * x339 + x1962 * x489 -
           x20 * x494 + x2159 +
           (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 +
            x124 * x143 + x175 * x68 + x20 * x68) *
             (x339 - x489) +
           (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
             (-x25 + x494 + (-x68) * (x339 - x489));
  J[117] = 6000000.0 * sc[1] * x207 + 12000000.0 * sc[3] * x371 -
           12000000.0 * sc[3] * x373 + 5000000.0 * sc[4] * x207 +
           50000000.0 * sc[4] * x81 - 44000000000.0 * sc[5] * x1163 +
           44000000000.0 * x1166 - x143 * x81 - x175 * x371 + x175 * x373 +
           x1962 * x500 - x1962 * x506 - x20 * x371 + x20 * x373 + x2087 +
           (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 +
            x124 * x143 + x175 * x68 + x20 * x68) *
             (-x500 + x506) +
           (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
             (x371 - x373 + (-x68) * (-x500 + x506));
  J[131] = 6000000.0 * sc[1] * x209 - 12000000.0 * sc[3] * x210 +
           12000000.0 * sc[3] * x378 + 5000000.0 * sc[4] * x209 +
           50000000.0 * sc[4] * x417 - 12000000.0 * sc[5] * x253 -
           50000000.0 * sc[5] * x287 - x139 * x1962 - x143 * x417 +
           x175 * x210 - x175 * x378 + x179 * x1962 + x20 * x210 - x20 * x378 +
           x2153 +
           (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 +
            x124 * x143 + x175 * x68 + x20 * x68) *
             (x139 - x179) +
           (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
             (-x210 + x378 + (-x68) * (x139 - x179));
  J[145] = 6000000.0 * sc[1] * x207 - 12000000.0 * sc[3] * x214 +
           12000000.0 * sc[3] * x52 + 50000000.0 * sc[4] * x185 +
           5000000.0 * sc[4] * x207 - x143 * x185 + x161 * x1962 + x175 * x214 -
           x175 * x52 - x195 * x1962 + x20 * x214 - x20 * x52 + x2153 +
           (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 +
            x124 * x143 + x175 * x68 + x20 * x68) *
             (-x161 + x195) +
           (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
             (-x214 + x52 + (-x68) * (-x161 + x195));
  J[159] = 6000000.0 * sc[1] * x248 + 12000000.0 * sc[3] * x202 +
           sc[4] * x2030 + 5000000.0 * sc[4] * x248 +
           50000000.0 * sc[4] * x359 - 5000000.0 * sc[5] * x228 - x143 * x359 -
           x175 * x202 + x1962 * x197 - x1962 * x205 - x20 * x202 + x2159 +
           (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 +
            x124 * x143 + x175 * x68 + x20 * x68) *
             (-x197 + x205) +
           (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
             (x202 - x25 + (-x68) * (-x197 + x205));
  J[173] = 6000000.0 * sc[1] * x33 - 12000000.0 * sc[3] * x282 +
           12000000.0 * sc[3] * x57 + 50000000.0 * sc[4] * x119 +
           5000000.0 * sc[4] * x33 - x119 * x143 + x169 * x1962 + x175 * x282 -
           x175 * x57 - x191 * x1962 + x20 * x282 - x20 * x57 + x2153 +
           (-12000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 +
            x124 * x143 + x175 * x68 + x20 * x68) *
             (-x169 + x191) +
           (-6000000.0 * sc[1] * x41 - 5000000.0 * sc[4] * x41) *
             (-x282 + x57 + (-x68) * (-x169 + x191));
  J[6] = sc[3] * x200 + 20000000.0 * sc[6] * x346 - 20000000.0 * sc[6] * x383 -
         x155 * x186 + x174 * x321 + x186 * x497 + x2164 +
         (x124 * x155 - x124 * x497) * (-x346 + x383) +
         (-x174 * x41) * (x222 - x327 + (-x68) * (-x346 + x383));
  J[20] = 375000000.0 * sc[12] * sc[3] * x1285 +
          300000000.0 * ((sc[3]) * (sc[3])) * x1285 + sc[3] * sc[5] * x1282 +
          sc[3] * x1279 - 20000000.0 * sc[6] * x128 +
          20000000.0 * sc[6] * x400 + x1277 + x15 * x174 - x155 * x381 + x2164 +
          x381 * x497 + (x124 * x155 - x124 * x497) * (x128 - x400) +
          (-x174 * x41) * (-x345 + x413 + (-x68) * (x128 - x400));
  J[34] = sc[5] * x1354 + 20000000.0 * sc[6] * x420 -
          20000000.0 * sc[6] * x428 - 20000000.0 * sc[6] - x155 * x72 +
          x174 * x28 + x2164 + x497 * x72 +
          (x124 * x155 - x124 * x497) * (-x420 + x428) +
          (-x174 * x41) * (-x134 + x347 + (-x68) * (-x420 + x428));
  J[48] = sc[0] * x200 + sc[11] * x1847 + 375000000.0 * sc[12] * sc[1] * x1285 +
          600000000.0 * sc[1] * sc[3] * x1285 + sc[1] * sc[5] * x1282 +
          sc[1] * x1279 + sc[4] * x1732 + sc[5] * x1826 +
          20000000.0 * sc[6] * x434 - 20000000.0 * sc[6] * x439 + sc[8] * x97 -
          300000000.0 * x1065 * x1389 +
          18000000.0 * x11 * exp(-452.89499922424801 * invT) - x155 * x392 +
          x174 * x389 + 7600000.0 * x379 * exp(-201.28666632188799 * invT) +
          x392 * x497 + (x124 * x155 - x124 * x497) * (-x434 + x439) +
          (-x174 * x41) * (x351 - x352 + (-x68) * (-x434 + x439));
  J[62] = sc[10] * x2017 + 20000000.0 * sc[11] * x333 + sc[3] * x1732 -
          20000000.0 * sc[6] * x103 + 20000000.0 * sc[6] * x445 +
          268000000.0 * x1409 - x155 * x393 + x174 * x397 - x1854 + x2164 +
          x393 * x497 + (x124 * x155 - x124 * x497) * (x103 - x445) +
          (-x174 * x41) * (-x356 + x84 + (-x68) * (x103 - x445));
  J[76] = sc[1] * sc[3] * x1282 + sc[2] * x1354 + sc[3] * x1826 -
          20000000.0 * sc[6] * x457 + 20000000.0 * sc[6] * x466 -
          x1065 * x1463 - x155 * x399 + x174 * x78 + x399 * x497 +
          (x124 * x155 - x124 * x497) * (x457 - x466) +
          (-x174 * x41) * (-x167 + x361 + (-x68) * (x457 - x466));
  J[90] =
    -375000000.0 * refC * sc[12] * x1389 - 300000000.0 * refC * sc[3] * x1389 -
    refC * sc[5] * x1463 - 2800000.0 * refC * x1142 * x1278 - sc[11] * x110 -
    20000000.0 * sc[2] + 20000000.0 * sc[6] * x477 - 20000000.0 * sc[6] * x486 -
    1000000.0 * sc[7] - sc[9] * x2059 - sc[9] * x286 + x1501 - x155 * x304 -
    20000000.0 * x165 + x174 * x232 - x1877 + x2164 - x30 * x379 + x304 * x497 +
    (x124 * x155 - x124 * x497) * (-x477 + x486) +
    (-x174 * x41) * (-x25 + x483 + (-x68) * (-x477 + x486));
  J[104] = -20000000.0 * sc[6] * x339 + 20000000.0 * sc[6] * x489 -
           1000000.0 * sc[6] - x155 * x402 + x174 * x419 + x2164 + x402 * x497 +
           (x124 * x155 - x124 * x497) * (x339 - x489) +
           (-x174 * x41) * (-x25 + x494 + (-x68) * (x339 - x489));
  J[118] = sc[3] * x97 + 20000000.0 * sc[6] * x500 - 20000000.0 * sc[6] * x506 -
           x155 * x81 + x174 * x207 + x2164 + x497 * x81 +
           (x124 * x155 - x124 * x497) * (-x500 + x506) +
           (-x174 * x41) * (x371 - x373 + (-x68) * (-x500 + x506));
  J[132] = -20000000.0 * sc[6] * x139 + 20000000.0 * sc[6] * x179 -
           sc[6] * x2059 - sc[6] * x286 - 2100000.0 * x1065 * x1142 +
           2100000.0 * x1147 - x155 * x417 + x174 * x209 + x417 * x497 +
           (x124 * x155 - x124 * x497) * (x139 - x179) +
           (-x174 * x41) * (-x210 + x378 + (-x68) * (x139 - x179));
  J[146] = sc[4] * x2017 + 20000000.0 * sc[6] * x161 -
           20000000.0 * sc[6] * x195 - 4200000.0 * x1065 * x1142 +
           4200000.0 * x1147 - x155 * x185 + x174 * x207 + x185 * x497 +
           (x124 * x155 - x124 * x497) * (-x161 + x195) +
           (-x174 * x41) * (-x214 + x52 + (-x68) * (-x161 + x195));
  J[160] = sc[3] * x1847 + 20000000.0 * sc[4] * x333 - sc[6] * x110 +
           20000000.0 * sc[6] * x197 - 20000000.0 * sc[6] * x205 - x155 * x359 +
           x174 * x248 + x2164 + x359 * x497 +
           (x124 * x155 - x124 * x497) * (-x197 + x205) +
           (-x174 * x41) * (x202 - x25 + (-x68) * (-x197 + x205));
  J[174] = 375000000.0 * sc[1] * sc[3] * x1285 + 20000000.0 * sc[6] * x169 -
           20000000.0 * sc[6] * x191 - 375000000.0 * x1065 * x1389 -
           x119 * x155 + x119 * x497 + x174 * x33 +
           (x124 * x155 - x124 * x497) * (-x169 + x191) +
           (-x174 * x41) * (-x282 + x57 + (-x68) * (-x169 + x191));
  J[7] = -26480000.0 * sc[7] * x186 - sc[7] * x188 + x1218 + x149 * x222 -
         x149 * x327 + x156 * x186 + x1966 - x2179 * x346 + x2179 * x383 +
         x2181 - x346 * x432 + x383 * x432 +
         (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 - x124 * x156 -
          x149 * x68 - x175 * x68) *
           (-x346 + x383) +
         (-12000000.0 * sc[1] * x41) * (x222 - x327 + (-x68) * (-x346 + x383));
  J[21] = sc[11] * x1311 - 26480000.0 * sc[7] * x381 + x1079 * x165 * x7 +
          x128 * x2179 + x128 * x432 + x1300 + x1322 - x149 * x345 +
          x149 * x413 + x156 * x381 + x1976 - x2179 * x400 + x2184 -
          x400 * x432 + x880 +
          (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 - x124 * x156 -
           x149 * x68 - x175 * x68) *
            (x128 - x400) +
          (-12000000.0 * sc[1] * x41) * (-x345 + x413 + (-x68) * (x128 - x400));
  J[35] =
    -70000000.0 * sc[0] * x134 + 70000000.0 * sc[0] * x347 -
    26480000.0 * sc[7] * x72 - 84300000.0 * sc[7] + x1322 - x134 * x149 +
    x149 * x347 + x156 * x72 + x1988 - x2179 * x420 + x2179 * x428 + x2184 -
    x220 * x420 + x220 * x428 - x420 * x432 + x428 * x432 +
    (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 - x124 * x156 -
     x149 * x68 - x175 * x68) *
      (-x420 + x428) +
    (-12000000.0 * sc[1] * x41) * (-x134 + x347 + (-x68) * (-x420 + x428));
  J[49] = 70000000.0 * sc[0] * x351 - 70000000.0 * sc[0] * x352 -
          sc[7] * x1834 - 26480000.0 * sc[7] * x392 + sc[8] * x97 + x1322 +
          x149 * x351 - x149 * x352 + x156 * x392 + x1999 - x2179 * x434 +
          x2179 * x439 + x2184 - x220 * x434 + x220 * x439 - x432 * x434 +
          x432 * x439 +
          (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 - x124 * x156 -
           x149 * x68 - x175 * x68) *
            (-x434 + x439) +
          (-12000000.0 * sc[1] * x41) * (x351 - x352 + (-x68) * (-x434 + x439));
  J[63] = -70000000.0 * sc[0] * x356 + 70000000.0 * sc[0] * x84 +
          sc[11] * x1850 - 26480000.0 * sc[7] * x393 + x103 * x2179 +
          x103 * x220 + x103 * x432 + x1322 - x149 * x356 + x149 * x84 +
          x156 * x393 + x2007 + x2028 - x2179 * x445 + x2184 - x220 * x445 -
          x432 * x445 +
          (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 - x124 * x156 -
           x149 * x68 - x175 * x68) *
            (x103 - x445) +
          (-12000000.0 * sc[1] * x41) * (-x356 + x84 + (-x68) * (x103 - x445));
  J[77] = -26480000.0 * sc[7] * x399 - sc[7] * x43 - 2.4576e+21 * x1062 * x131 +
          x1071 * (1.68 * x1067 + 12.0 * x1069) + 1.92e+16 * x1092 -
          1.92e+16 * x1102 + 2.4576e+21 * x131 * x326 -
          x1334 * (1.68 * x1067 + 12.0 * x1069) + x1480 - x149 * x167 +
          x149 * x361 + x156 * x399 + x2050 + x2179 * x457 - x2179 * x466 +
          x432 * x457 - x432 * x466 +
          (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 - x124 * x156 -
           x149 * x68 - x175 * x68) *
            (x457 - x466) +
          (-12000000.0 * sc[1] * x41) * (-x167 + x361 + (-x68) * (x457 - x466));
  J[91] = 70000000.0 * sc[0] * x483 - 26480000.0 * sc[7] * x304 -
          1000000.0 * sc[7] + x149 * x483 + x156 * x304 + x2067 - x2179 * x477 +
          x2179 * x486 + x2196 - x220 * x477 + x220 * x486 - x432 * x477 +
          x432 * x486 +
          (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 - x124 * x156 -
           x149 * x68 - x175 * x68) *
            (-x477 + x486) +
          (-12000000.0 * sc[1] * x41) * (-x25 + x483 + (-x68) * (-x477 + x486));
  J[105] =
    -sc[0] * x188 - sc[11] * x2204 - 84300000.0 * sc[2] - sc[3] * x1834 -
    sc[5] * x43 - 1000000.0 * sc[6] - 26480000.0 * sc[7] * x402 -
    9.84 * sc[7] * exp(-g_RT_qss[0] + x2176 - x260) -
    64000000.0 * sc[7] * exp(-g_RT_qss[1] + 286.83349950869001 * invT - x260) -
    x1078 + x149 * x494 + x1541 + x156 * x402 - x175 * x25 + x2086 +
    x2179 * x339 - x2179 * x489 + x2196 + 26480000.0 * x299 +
    26480000.0 * x322 + x339 * x432 - x432 * x489 +
    (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 - x124 * x156 -
     x149 * x68 - x175 * x68) *
      (x339 - x489) +
    (-12000000.0 * sc[1] * x41) * (-x25 + x494 + (-x68) * (x339 - x489));
  J[119] =
    12000000.0 * sc[1] * x207 + sc[3] * x97 - 26480000.0 * sc[7] * x81 +
    26480000.0 * sc[9] * x311 + x149 * x371 - x149 * x373 + x156 * x81 + x1563 +
    4.92 * x165 * exp(x2176) + x2087 + x2101 - x2179 * x500 + x2179 * x506 +
    x2181 + x379 * x94 - x432 * x500 + x432 * x506 +
    32000000.0 * x50 * exp(286.83349950869001 * invT) + x994 +
    (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 - x124 * x156 -
     x149 * x68 - x175 * x68) *
      (-x500 + x506) +
    (-12000000.0 * sc[1] * x41) * (x371 - x373 + (-x68) * (-x500 + x506));
  J[133] =
    -26480000.0 * sc[7] * x417 + 26480000.0 * sc[8] * x311 -
    6.144e+20 * x1062 * x131 + x1071 * (0.42 * x1067 + 3.0 * x1069) +
    4800000000000000.0 * x1092 - 4800000000000000.0 * x1102 + x1175 +
    6.144e+20 * x131 * x326 - x1334 * (0.42 * x1067 + 3.0 * x1069) +
    x139 * x2179 + x139 * x432 - x149 * x210 + x149 * x378 + x156 * x417 -
    x179 * x2179 - x179 * x432 + x2111 +
    (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 - x124 * x156 -
     x149 * x68 - x175 * x68) *
      (x139 - x179) +
    (-12000000.0 * sc[1] * x41) * (-x210 + x378 + (-x68) * (x139 - x179));
  J[147] =
    12000000.0 * sc[1] * x207 - 26480000.0 * sc[7] * x185 + x1176 -
    x149 * x214 + x149 * x52 + x156 * x185 - x161 * x1962 - x161 * x2179 -
    x161 * x432 - x175 * x214 + x175 * x52 + x195 * x1962 + x195 * x2179 +
    x195 * x432 + x2181 +
    (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 - x124 * x156 -
     x149 * x68 - x175 * x68) *
      (-x161 + x195) +
    (-12000000.0 * sc[1] * x41) * (-x214 + x52 + (-x68) * (-x161 + x195));
  J[161] = 70000000.0 * sc[0] * x202 + sc[1] * x1311 +
           12000000.0 * sc[1] * x248 + sc[4] * x1850 - sc[7] * x2204 -
           26480000.0 * sc[7] * x359 + x149 * x202 + x156 * x359 + x175 * x202 -
           x175 * x25 - x1962 * x197 + x1962 * x205 - x197 * x2179 -
           x197 * x220 - x197 * x432 + x205 * x2179 + x205 * x220 +
           x205 * x432 + x2196 +
           (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 -
            x124 * x156 - x149 * x68 - x175 * x68) *
             (-x197 + x205) +
           (-12000000.0 * sc[1] * x41) * (x202 - x25 + (-x68) * (-x197 + x205));
  J[175] =
    -26480000.0 * sc[7] * x119 + x119 * x156 - x149 * x282 + x149 * x57 -
    x169 * x2179 - x169 * x432 + x191 * x2179 + x191 * x432 + x2152 + x2184 +
    x471 +
    (-70000000.0 * sc[0] * x68 + 26480000.0 * sc[7] * x124 - x124 * x156 -
     x149 * x68 - x175 * x68) *
      (-x169 + x191) +
    (-12000000.0 * sc[1] * x41) * (-x282 + x57 + (-x68) * (-x169 + x191));
  J[8] =
    26480000.0 * sc[7] * x186 + sc[7] * x188 + x1200 - x156 * x186 -
    x222 * x87 + x2227 * x346 - x2227 * x383 + x327 * x87 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (-x346 + x383);
  J[22] =
    2.477e+21 * sc[7] * x1182 * x7 * exp(x1185) + 26480000.0 * sc[7] * x381 -
    x128 * x2227 + x1316 - x156 * x381 + x2227 * x400 + x345 * x87 -
    x413 * x87 - x880 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (x128 - x400);
  J[36] =
    26480000.0 * sc[7] * x72 - sc[8] * x1709 + x1316 + x134 * x87 - x156 * x72 +
    x2227 * x420 - x2227 * x428 - x347 * x87 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (-x420 + x428);
  J[50] =
    26480000.0 * sc[7] * x392 - sc[8] * x97 + x1316 - x156 * x392 +
    x2227 * x434 - x2227 * x439 - x351 * x87 + x352 * x87 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (-x434 + x439);
  J[64] =
    26480000.0 * sc[7] * x393 + sc[7] * x92 - x103 * x2227 + x1316 -
    x156 * x393 - x2007 + x2227 * x445 + x356 * x87 - x84 * x87 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (x103 - x445);
  J[78] =
    26480000.0 * sc[7] * x399 + sc[7] * x43 +
    2.8986751181102398e+33 * x1172 * x1179 -
    2.8986751181102398e+33 * x1172 * x1192 - x1173 * x1473 -
    1.4862e+22 * x1201 + x1204 * x1473 + 1.4862e+22 * x1210 - x156 * x399 +
    x167 * x87 - x2227 * x457 + x2227 * x466 - x361 * x87 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (x457 - x466);
  J[92] =
    26480000.0 * sc[7] * x304 + 1000000.0 * sc[7] + x1316 - x156 * x304 +
    x2227 * x477 - x2227 * x486 + x25 * x87 - x483 * x87 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (-x477 + x486);
  J[106] =
    sc[0] * x188 + sc[11] * x2204 +
    2.477e+21 * sc[1] * x1182 * x7 * exp(x1185) + sc[4] * x92 + sc[5] * x43 +
    1000000.0 * sc[6] + 26480000.0 * sc[7] * x402 +
    4.92 * sc[7] * exp(-g_RT_qss[0] + x2176 - x260) +
    32000000.0 * sc[7] * exp(-g_RT_qss[1] + 286.83349950869001 * invT - x260) +
    x1316 - x156 * x402 - x2227 * x339 + x2227 * x489 + x25 * x87 -
    26480000.0 * x299 - 26480000.0 * x322 - x494 * x87 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (x339 - x489);
  J[120] =
    -sc[2] * x1709 - sc[3] * x97 + 26480000.0 * sc[7] * x81 -
    26480000.0 * sc[9] * x311 - 2.477e+21 * x1182 * x1186 * x204 + x1200 -
    x156 * x81 - 2.46 * x165 * exp(x2176) - x2087 + x2227 * x500 -
    x2227 * x506 - x371 * x87 + x373 * x87 - x379 * x94 -
    16000000.0 * x50 * exp(286.83349950869001 * invT) - x994 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (-x500 + x506);
  J[134] =
    26480000.0 * sc[7] * x417 - 26480000.0 * sc[8] * x311 +
    7.2466877952755894e+32 * x1172 * x1179 -
    7.2466877952755894e+32 * x1172 * x1192 - x1173 * x429 - 3.7155e+21 * x1201 +
    x1204 * x429 + 3.7155e+21 * x1210 - x139 * x2227 - x156 * x417 +
    x179 * x2227 + x210 * x87 - x378 * x87 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (x139 - x179);
  J[148] =
    26480000.0 * sc[7] * x185 + x1200 - x156 * x185 + x161 * x2227 -
    x195 * x2227 + x214 * x87 - x52 * x87 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (-x161 + x195);
  J[162] =
    sc[7] * x2204 + 26480000.0 * sc[7] * x359 + x1316 - x156 * x359 +
    x197 * x2227 - x202 * x87 - x205 * x2227 + x25 * x87 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (-x197 + x205);
  J[176] =
    26480000.0 * sc[7] * x119 - x119 * x156 + x1316 + x169 * x2227 -
    x191 * x2227 + x282 * x87 - x57 * x87 +
    (-26480000.0 * sc[7] * x124 + x124 * x156 + x68 * x87) * (-x169 + x191);
  J[9] = 14000000.0 * sc[10] * x222 - 14000000.0 * sc[10] * x327 +
         30000000.0 * sc[2] * x186 + 40000000.0 * sc[3] * x222 -
         40000000.0 * sc[3] * x327 + 50000000.0 * sc[4] * x186 +
         26480000.0 * sc[7] * x186 + x1160 - 1204.0 * x1238 + 1204.0 * x17 +
         x183 * x186 + x186 * x460 + x186 * x497 + x845 + x853 +
         (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
          30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
          40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
          26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
           (-x346 + x383);
  J[23] =
    -refCinv * sc[9] * x1250 + sc[10] * x1267 - 14000000.0 * sc[10] * x345 +
    14000000.0 * sc[10] * x413 + 30000000.0 * sc[2] * x381 -
    40000000.0 * sc[3] * x345 + 40000000.0 * sc[3] * x413 +
    50000000.0 * sc[4] * x381 + 26480000.0 * sc[7] * x381 + x1296 + x2240 +
    x381 * x497 + x901 +
    (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
     40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
      (x128 - x400);
  J[37] =
    sc[10] * x10 - 14000000.0 * sc[10] * x134 + 14000000.0 * sc[10] * x347 +
    73400000.0 * sc[1] * x72 - 15000000.0 * sc[2] * x134 +
    15000000.0 * sc[2] * x347 - 40000000.0 * sc[3] * x134 +
    40000000.0 * sc[3] * x347 + 50000000.0 * sc[4] * x72 +
    26480000.0 * sc[7] * x72 - sc[9] * x1712 * x468 + x1369 + x183 * x72 +
    x2240 + x460 * x72 + x497 * x72 + 15000000.0 * x50 +
    (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
     40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
      (-x420 + x428);
  J[51] =
    14000000.0 * sc[10] * x351 - 14000000.0 * sc[10] * x352 +
    73400000.0 * sc[1] * x392 + 15000000.0 * sc[2] * x351 -
    15000000.0 * sc[2] * x352 + 30000000.0 * sc[2] * x392 +
    40000000.0 * sc[3] * x351 - 40000000.0 * sc[3] * x352 +
    50000000.0 * sc[4] * x392 + 26480000.0 * sc[7] * x392 - 3612.0 * x1238 +
    x1344 + 3612.0 * x17 - x1730 + x183 * x392 +
    7600000.0 * x379 * exp(-201.28666632188799 * invT) + x392 * x460 +
    x392 * x497 + 40000000.0 * x50 + x894 +
    (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
     40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
      (-x434 + x439);
  J[65] =
    sc[10] * x2017 - 14000000.0 * sc[10] * x356 + 14000000.0 * sc[10] * x84 +
    73400000.0 * sc[1] * x393 - 15000000.0 * sc[2] * x356 +
    30000000.0 * sc[2] * x393 + 15000000.0 * sc[2] * x84 -
    40000000.0 * sc[3] * x356 + 40000000.0 * sc[3] * x84 +
    50000000.0 * sc[4] * x393 + 26480000.0 * sc[7] * x393 -
    30000000.0 * sc[9] * x295 - x1407 + x183 * x393 + x2240 -
    28000000.0 * x270 * x271 - 50000000.0 * x299 - 50000000.0 * x322 +
    x393 * x460 + x393 * x497 +
    (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
     40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
      (x103 - x445);
  J[79] =
    -14000000.0 * sc[10] * x167 + 14000000.0 * sc[10] * x361 -
    40000000.0 * sc[3] * x167 + 40000000.0 * sc[3] * x361 +
    50000000.0 * sc[4] * x399 + 26480000.0 * sc[7] * x399 -
    12000000.0 * sc[9] * x253 - 50000000.0 * sc[9] * x287 - 3612.0 * x1238 +
    x1472 + 3612.0 * x17 + x399 * x497 + x963 +
    (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
     40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
      (x457 - x466);
  J[93] =
    14000000.0 * sc[10] * x483 + 73400000.0 * sc[1] * x304 +
    30000000.0 * sc[2] * x304 + 15000000.0 * sc[2] * x483 +
    40000000.0 * sc[3] * x483 + 50000000.0 * sc[4] * x304 +
    26480000.0 * sc[7] * x304 - sc[9] * x2059 - sc[9] * x286 + x183 * x304 +
    x2258 + x2260 + x304 * x460 + x304 * x497 +
    (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
     40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
      (-x477 + x486);
  J[107] =
    14000000.0 * sc[10] * x494 + 73400000.0 * sc[1] * x402 +
    30000000.0 * sc[2] * x402 + 15000000.0 * sc[2] * x494 +
    40000000.0 * sc[3] * x494 + 50000000.0 * sc[4] * x402 +
    2244000000000.0 * sc[5] * x402 * x450 + 26480000.0 * sc[7] * x402 + x2258 +
    x2260 - 26480000.0 * x299 - 26480000.0 * x322 +
    187000000000.0 * x402 * x450 * x6 + x402 * x497 +
    (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
     40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
      (x339 - x489);
  J[121] =
    14000000.0 * sc[10] * x371 - 14000000.0 * sc[10] * x373 +
    73400000.0 * sc[1] * x81 + 15000000.0 * sc[2] * x371 -
    15000000.0 * sc[2] * x373 + 30000000.0 * sc[2] * x81 +
    40000000.0 * sc[3] * x371 - 40000000.0 * sc[3] * x373 +
    50000000.0 * sc[4] * x81 + 26480000.0 * sc[7] * x81 -
    26480000.0 * sc[9] * x311 + x1160 - 1204.0 * x1238 + 1204.0 * x17 +
    x183 * x81 + x460 * x81 + x497 * x81 + x845 +
    (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
     40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
      (-x500 + x506);
  J[135] =
    -14000000.0 * sc[10] * x210 + 14000000.0 * sc[10] * x378 -
    14000000.0 * sc[11] * x274 - sc[2] * x1712 * x468 -
    40000000.0 * sc[3] * x210 + 40000000.0 * sc[3] * x378 -
    30000000.0 * sc[4] * x295 + 50000000.0 * sc[4] * x417 -
    12000000.0 * sc[5] * x253 - 50000000.0 * sc[5] * x287 - sc[6] * x2059 -
    sc[6] * x286 + 26480000.0 * sc[7] * x417 - 26480000.0 * sc[8] * x311 +
    x1158 - 903.0 * x1238 - x1250 * x269 - x1565 + 903.0 * x17 - x1772 + x283 +
    x417 * x497 +
    (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
     40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
      (x139 - x179);
  J[149] =
    refC * x1230 * x1712 - 14000000.0 * sc[10] * x214 +
    14000000.0 * sc[10] * x52 + sc[1] * x1267 + 73400000.0 * sc[1] * x185 +
    sc[2] * x10 + 30000000.0 * sc[2] * x185 - 15000000.0 * sc[2] * x214 +
    15000000.0 * sc[2] * x52 - 40000000.0 * sc[3] * x214 +
    40000000.0 * sc[3] * x52 + 50000000.0 * sc[4] * x185 + sc[4] * x2017 +
    26480000.0 * sc[7] * x185 + x1160 - 2107.0 * x1238 + 2107.0 * x17 +
    x183 * x185 + x185 * x460 + x185 * x497 + 14000000.0 * x50 + x845 +
    (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
     40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
      (-x161 + x195);
  J[163] =
    14000000.0 * sc[10] * x202 + 73400000.0 * sc[1] * x359 +
    15000000.0 * sc[2] * x202 + 30000000.0 * sc[2] * x359 +
    40000000.0 * sc[3] * x202 + 50000000.0 * sc[4] * x359 +
    2244000000000.0 * sc[5] * x359 * x450 + 26480000.0 * sc[7] * x359 -
    14000000.0 * sc[9] * x274 + 5070000000000000.0 * x204 * x342 * x815 +
    x2258 + x2260 + 187000000000.0 * x359 * x450 * x6 + x359 * x497 +
    (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
     40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
      (-x197 + x205);
  J[177] =
    -14000000.0 * sc[10] * x282 + 14000000.0 * sc[10] * x57 +
    73400000.0 * sc[1] * x119 + 30000000.0 * sc[2] * x119 -
    15000000.0 * sc[2] * x282 + 15000000.0 * sc[2] * x57 -
    40000000.0 * sc[3] * x282 + 40000000.0 * sc[3] * x57 +
    50000000.0 * sc[4] * x119 + 26480000.0 * sc[7] * x119 + x119 * x183 +
    x119 * x460 + x119 * x497 + x2240 +
    (-14000000.0 * sc[10] * x68 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x68 -
     40000000.0 * sc[3] * x68 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x183 - x124 * x460 - x124 * x497) *
      (-x169 + x191);
  J[10] =
    -14000000.0 * sc[10] * x222 + 14000000.0 * sc[10] * x327 +
    30000000.0 * sc[2] * x186 + 1204.0 * x1238 - 1204.0 * x17 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (-x346 + x383);
  J[24] =
    -sc[10] * x1267 - 30000000.0 * sc[10] * x308 + 14000000.0 * sc[10] * x345 -
    14000000.0 * sc[10] * x413 + 30000000.0 * sc[2] * x381 + x2281 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (x128 - x400);
  J[38] =
    -sc[10] * x10 + 14000000.0 * sc[10] * x134 - 14000000.0 * sc[10] * x347 +
    sc[9] * x1712 * x468 + x1369 + x2281 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (-x420 + x428);
  J[52] =
    -14000000.0 * sc[10] * x351 + 14000000.0 * sc[10] * x352 +
    30000000.0 * sc[2] * x392 + 3612.0 * x1238 - 3612.0 * x17 + x1730 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (-x434 + x439);
  J[66] =
    -sc[10] * x2017 + 14000000.0 * sc[10] * x356 - 14000000.0 * sc[10] * x84 +
    30000000.0 * sc[2] * x393 + x1407 + x2281 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (x103 - x445);
  J[80] =
    14000000.0 * sc[10] * x167 - 14000000.0 * sc[10] * x361 +
    30000000.0 * sc[2] * x399 + 3612.0 * x1238 - 3612.0 * x17 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (x457 - x466);
  J[94] =
    14000000.0 * sc[10] * x25 - 14000000.0 * sc[10] * x483 +
    30000000.0 * sc[2] * x304 + sc[9] * x2059 + x2281 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (-x477 + x486);
  J[108] =
    14000000.0 * sc[10] * x25 - 14000000.0 * sc[10] * x494 +
    30000000.0 * sc[2] * x402 + x2281 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (x339 - x489);
  J[122] =
    -14000000.0 * sc[10] * x371 + 14000000.0 * sc[10] * x373 +
    30000000.0 * sc[2] * x81 + 1204.0 * x1238 - 1204.0 * x17 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (-x500 + x506);
  J[136] =
    14000000.0 * sc[10] * x210 - 14000000.0 * sc[10] * x378 +
    14000000.0 * sc[11] * x274 + sc[2] * x1712 * x468 +
    30000000.0 * sc[2] * x417 + sc[6] * x2059 + 903.0 * x1238 + x1565 -
    903.0 * x17 + x1772 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (x139 - x179);
  J[150] =
    -refC * x1230 * x1712 + 14000000.0 * sc[10] * x214 -
    14000000.0 * sc[10] * x52 - sc[1] * x1267 - 30000000.0 * sc[1] * x308 -
    sc[2] * x10 + 30000000.0 * sc[2] * x185 - sc[4] * x2017 + 2107.0 * x1238 -
    2107.0 * x17 - 14000000.0 * x50 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (-x161 + x195);
  J[164] =
    -14000000.0 * sc[10] * x202 + 14000000.0 * sc[10] * x25 +
    30000000.0 * sc[2] * x359 + 14000000.0 * sc[9] * x274 + x2281 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (-x197 + x205);
  J[178] =
    14000000.0 * sc[10] * x282 - 14000000.0 * sc[10] * x57 +
    30000000.0 * sc[2] * x119 + x2281 +
    (14000000.0 * sc[10] * x68 - 30000000.0 * sc[2] * x124) * (-x169 + x191);
  J[11] =
    14000000.0 * sc[10] * x222 - 14000000.0 * sc[10] * x327 -
    20000000.0 * sc[11] * x234 + 20000000.0 * sc[1] * x321 +
    10000000.0 * sc[2] * x321 + 5000000.0 * sc[4] * x321 -
    20000000.0 * sc[6] * x346 + 20000000.0 * sc[6] * x383 +
    sc[9] * x839 * exp(x382) + x1243 + x143 * x186 + x155 * x186 + x156 * x186 +
    x174 * x321 + x186 * x443 + x186 * x476 + x827 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (-x346 + x383) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (x222 - x327 + (-x68) * (-x346 + x383));
  J[25] =
    -14000000.0 * sc[10] * x345 + 14000000.0 * sc[10] * x413 +
    20000000.0 * sc[1] * x15 + 10000000.0 * sc[2] * x15 +
    5000000.0 * sc[4] * x15 + 20000000.0 * sc[6] * x128 -
    20000000.0 * sc[6] * x400 + x1315 + x1325 + x143 * x381 + x15 * x174 +
    x155 * x381 + x156 * x381 - 20000000.0 * x249 +
    1350000000000.0 * x379 * x438 * x7 * exp(x427) + x381 * x443 + x381 * x476 +
    x890 - 20000000.0 * x902 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (x128 - x400) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (-x345 + x413 + (-x68) * (x128 - x400));
  J[39] =
    -14000000.0 * sc[10] * x134 + 14000000.0 * sc[10] * x347 - sc[11] * x1725 +
    20000000.0 * sc[1] * x28 + 5000000.0 * sc[4] * x28 + x1325 + x1378 +
    x155 * x72 + x156 * x72 + x174 * x28 + x1992 + x443 * x72 + x476 * x72 +
    x890 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (-x420 + x428) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (-x134 + x347 + (-x68) * (-x420 + x428));
  J[53] =
    14000000.0 * sc[10] * x351 - 14000000.0 * sc[10] * x352 - sc[11] * x1847 +
    20000000.0 * sc[1] * x389 + 5000000.0 * sc[4] * x389 +
    18000000.0 * x11 * exp(-452.89499922424801 * invT) + x1325 + x1405 +
    x155 * x392 + x156 * x392 + x174 * x389 + x2003 + x392 * x443 +
    x392 * x476 + x890 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (-x434 + x439) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (x351 - x352 + (-x68) * (-x434 + x439));
  J[67] =
    -14000000.0 * sc[10] * x356 + 14000000.0 * sc[10] * x84 +
    20000000.0 * sc[1] * x397 + 5000000.0 * sc[4] * x397 + x1325 + x1434 +
    x155 * x393 + x156 * x393 + x174 * x397 + x2035 - 5000000.0 * x249 +
    x379 * x475 + x393 * x443 + x393 * x476 + x890 - 5000000.0 * x902 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (x103 - x445) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (-x356 + x84 + (-x68) * (x103 - x445));
  J[81] =
    -14000000.0 * sc[10] * x167 + 14000000.0 * sc[10] * x361 -
    5000000.0 * sc[11] * x228 + 20000000.0 * sc[1] * x78 +
    5000000.0 * sc[4] * x78 + 8100000000000.0 * x1063 + x1089 * x1442 -
    1.0032110091743101e+19 * x1105 * x70 - x1442 * x75 + x1492 + x155 * x399 +
    x156 * x399 + x174 * x78 + x2055 - 3.042e+16 * x314 * x342 * x815 -
    8100000000000.0 * x314 * x826 + 1.0032110091743101e+19 * x316 * x70 +
    3.042e+16 * x342 * x837 + x399 * x443 + x399 * x476 +
    3.5867302325581401e+30 * x54 * x822 - 3.5867302325581401e+30 * x54 * x843 +
    x66 * x7 * x831 - x66 * x837 * x839 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (x457 - x466) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (-x167 + x361 + (-x68) * (x457 - x466));
  J[95] =
    14000000.0 * sc[10] * x483 - sc[11] * x110 + 20000000.0 * sc[1] * x232 +
    5000000.0 * sc[4] * x232 - 30000000.0 * sc[4] * x25 + x1325 + x1516 +
    x155 * x304 + x156 * x304 + x174 * x232 + x2072 + x2258 + x30 * x379 +
    x304 * x443 + x304 * x476 + x890 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (-x477 + x486) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (-x25 + x483 + (-x68) * (-x477 + x486));
  J[109] =
    14000000.0 * sc[10] * x494 - sc[11] * x2204 + 20000000.0 * sc[1] * x419 -
    30000000.0 * sc[4] * x25 + 5000000.0 * sc[4] * x419 + x1325 + x1539 +
    x155 * x402 + x156 * x402 + x174 * x419 + x2080 + x2258 + x402 * x443 +
    x402 * x476 + x890 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (x339 - x489) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (-x25 + x494 + (-x68) * (x339 - x489));
  J[123] =
    14000000.0 * sc[10] * x371 - 14000000.0 * sc[10] * x373 +
    20000000.0 * sc[1] * x207 + 10000000.0 * sc[2] * x207 +
    5000000.0 * sc[4] * x207 - 20000000.0 * sc[6] * x500 +
    20000000.0 * sc[6] * x506 + x143 * x81 + x155 * x81 + x156 * x81 + x1561 +
    x174 * x207 + x379 * x94 + x443 * x81 + x476 * x81 + x827 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (-x500 + x506) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (x371 - x373 + (-x68) * (-x500 + x506));
  J[137] =
    sc[0] * x839 * exp(x382) - 14000000.0 * sc[10] * x210 +
    14000000.0 * sc[10] * x378 - 14000000.0 * sc[11] * x274 +
    20000000.0 * sc[1] * x209 + 10000000.0 * sc[2] * x209 +
    5000000.0 * sc[4] * x209 + 20000000.0 * sc[6] * x139 -
    20000000.0 * sc[6] * x179 + 2025000000000.0 * x1063 + x1084 * x1089 -
    x1084 * x75 - 2.5080275229357798e+18 * x1105 * x70 + x140 + x143 * x417 +
    x155 * x417 + x156 * x417 + x174 * x209 + x223 * x833 - x223 * x841 -
    7605000000000000.0 * x314 * x342 * x815 - 2025000000000.0 * x314 * x826 +
    2.5080275229357798e+18 * x316 * x70 + 7605000000000000.0 * x342 * x837 +
    x417 * x443 + x417 * x476 + 8.9668255813953501e+29 * x54 * x822 -
    8.9668255813953501e+29 * x54 * x843 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (x139 - x179) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (-x210 + x378 + (-x68) * (x139 - x179));
  J[151] =
    -14000000.0 * sc[10] * x214 + 14000000.0 * sc[10] * x52 +
    20000000.0 * sc[1] * x207 + 10000000.0 * sc[2] * x207 +
    5000000.0 * sc[4] * x207 - 20000000.0 * sc[6] * x161 +
    20000000.0 * sc[6] * x195 + x1197 + x143 * x185 + x155 * x185 +
    x156 * x185 + x174 * x207 + x185 * x443 + x185 * x476 + 14000000.0 * x50 +
    x827 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (-x161 + x195) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (-x214 + x52 + (-x68) * (-x161 + x195));
  J[165] =
    -20000000.0 * sc[0] * x234 + 14000000.0 * sc[10] * x202 +
    20000000.0 * sc[1] * x248 - sc[2] * x1725 - sc[3] * x1847 +
    5000000.0 * sc[4] * x248 - 30000000.0 * sc[4] * x25 -
    5000000.0 * sc[5] * x228 - sc[6] * x110 - sc[7] * x2204 -
    14000000.0 * sc[9] * x274 + x1325 + x155 * x359 + x156 * x359 +
    x174 * x248 - 5070000000000000.0 * x204 * x342 * x815 -
    1350000000000.0 * x204 * x826 + x2141 + x2258 + x359 * x443 + x359 * x476 +
    x462 + x890 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (-x197 + x205) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (x202 - x25 + (-x68) * (-x197 + x205));
  J[179] =
    -14000000.0 * sc[10] * x282 + 14000000.0 * sc[10] * x57 +
    20000000.0 * sc[1] * x33 + 10000000.0 * sc[2] * x33 +
    5000000.0 * sc[4] * x33 - 20000000.0 * sc[6] * x169 +
    20000000.0 * sc[6] * x191 + x119 * x143 + x119 * x155 + x119 * x156 +
    x119 * x443 + x119 * x476 + x1217 + x174 * x33 + x890 +
    (-14000000.0 * sc[10] * x68 - 30000000.0 * sc[4] * x68 - x124 * x143 -
     x124 * x155 - x124 * x156 - x124 * x375 - x124 * x443 - x124 * x476) *
      (-x169 + x191) +
    (-20000000.0 * sc[1] * x41 - 10000000.0 * sc[2] * x41 -
     5000000.0 * sc[4] * x41 - x174 * x41 - x41 * x472) *
      (-x282 + x57 + (-x68) * (-x169 + x191));
  J[12] = 0;
  J[26] = 0;
  J[40] = 0;
  J[54] = 0;
  J[68] = 0;
  J[82] = 0;
  J[96] = 0;
  J[110] = 0;
  J[124] = 0;
  J[138] = 0;
  J[152] = 0;
  J[166] = 0;
  J[180] = 0;

  for (int k = 0; k < 13; k++) {
    J[182 + k] = (wdot_pert1[k] - wdot[k]) / (pertT);
  }
  amrex::Real c_R[13], dcRdT[13], e_RT[13];
  amrex::Real* eh_RT;
  if (consP == 1) {
    cp_R(c_R, T);
    dcvpRdT(dcRdT, T);
    eh_RT = &h_RT[0];
  } else {
    cv_R(c_R, T);
    dcvpRdT(dcRdT, T);
    speciesInternalEnergy(e_RT, T);
    eh_RT = &e_RT[0];
  }

  amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT = 0.0, dehmixdT = 0.0;
  for (int k = 0; k < 13; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[182 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 13; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 13; ++m) {
      dehmixdc += eh_RT[m] * J[k * 14 + m];
    }
    J[k * 14 + 13] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[195] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
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
  amrex::Real c[13];

  for (int k = 0; k < 13; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian(J, c, *Tp, *consP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 13; k++) {
    J[182 + k] *= 1.e-6;
    J[k * 14 + 13] *= 1.e6;
  }
}

// Transport function declarations

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENIMC(int* LENIMC)
{
  *LENIMC = 54;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENRMC(int* LENRMC)
{
  *LENRMC = 3666;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNO(int* NO)
{
  *NO = 4;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKK(int* KK)
{
  *KK = 13;
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
  WT[1] = 1.00800000E+00;
  WT[2] = 1.59990000E+01;
  WT[3] = 3.19980000E+01;
  WT[4] = 1.70070000E+01;
  WT[5] = 1.80150000E+01;
  WT[6] = 3.30060000E+01;
  WT[7] = 1.50350000E+01;
  WT[8] = 1.60430000E+01;
  WT[9] = 2.80100000E+01;
  WT[10] = 4.40090000E+01;
  WT[11] = 3.00260000E+01;
  WT[12] = 2.80140000E+01;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetEPS(amrex::Real* EPS)
{
  EPS[0] = 3.80000000E+01;
  EPS[1] = 1.45000000E+02;
  EPS[2] = 8.00000000E+01;
  EPS[3] = 1.07400000E+02;
  EPS[4] = 8.00000000E+01;
  EPS[5] = 5.72400000E+02;
  EPS[6] = 1.07400000E+02;
  EPS[7] = 1.44000000E+02;
  EPS[8] = 1.41400000E+02;
  EPS[9] = 9.81000000E+01;
  EPS[10] = 2.44000000E+02;
  EPS[11] = 4.98000000E+02;
  EPS[12] = 9.75300000E+01;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetSIG(amrex::Real* SIG)
{
  SIG[0] = 2.92000000E+00;
  SIG[1] = 2.05000000E+00;
  SIG[2] = 2.75000000E+00;
  SIG[3] = 3.45800000E+00;
  SIG[4] = 2.75000000E+00;
  SIG[5] = 2.60500000E+00;
  SIG[6] = 3.45800000E+00;
  SIG[7] = 3.80000000E+00;
  SIG[8] = 3.74600000E+00;
  SIG[9] = 3.65000000E+00;
  SIG[10] = 3.76300000E+00;
  SIG[11] = 3.59000000E+00;
  SIG[12] = 3.62100000E+00;
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
  DIP[5] = 1.84400000E+00;
  DIP[6] = 0.00000000E+00;
  DIP[7] = 0.00000000E+00;
  DIP[8] = 0.00000000E+00;
  DIP[9] = 0.00000000E+00;
  DIP[10] = 0.00000000E+00;
  DIP[11] = 0.00000000E+00;
  DIP[12] = 0.00000000E+00;
}

// the polarizability in cubic Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPOL(amrex::Real* POL)
{
  POL[0] = 7.90000000E-01;
  POL[1] = 0.00000000E+00;
  POL[2] = 0.00000000E+00;
  POL[3] = 1.60000000E+00;
  POL[4] = 0.00000000E+00;
  POL[5] = 0.00000000E+00;
  POL[6] = 0.00000000E+00;
  POL[7] = 0.00000000E+00;
  POL[8] = 2.60000000E+00;
  POL[9] = 1.95000000E+00;
  POL[10] = 2.65000000E+00;
  POL[11] = 0.00000000E+00;
  POL[12] = 1.76000000E+00;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetZROT(amrex::Real* ZROT)
{
  ZROT[0] = 2.80000000E+02;
  ZROT[1] = 0.00000000E+00;
  ZROT[2] = 0.00000000E+00;
  ZROT[3] = 3.80000000E+00;
  ZROT[4] = 0.00000000E+00;
  ZROT[5] = 4.00000000E+00;
  ZROT[6] = 1.00000000E+00;
  ZROT[7] = 0.00000000E+00;
  ZROT[8] = 1.30000000E+01;
  ZROT[9] = 1.80000000E+00;
  ZROT[10] = 2.10000000E+00;
  ZROT[11] = 2.00000000E+00;
  ZROT[12] = 4.00000000E+00;
}

// 0: monoatomic, 1: linear, 2: nonlinear
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLIN(int* NLIN)
{
  NLIN[0] = 1;
  NLIN[1] = 0;
  NLIN[2] = 0;
  NLIN[3] = 1;
  NLIN[4] = 1;
  NLIN[5] = 2;
  NLIN[6] = 2;
  NLIN[7] = 1;
  NLIN[8] = 2;
  NLIN[9] = 1;
  NLIN[10] = 1;
  NLIN[11] = 2;
  NLIN[12] = 1;
}

// Poly fits for the viscosities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFETA(amrex::Real* COFETA)
{
  COFETA[0] = -1.37549291E+01;
  COFETA[1] = 9.65530587E-01;
  COFETA[2] = -4.45720114E-02;
  COFETA[3] = 2.05871810E-03;
  COFETA[4] = -1.98744352E+01;
  COFETA[5] = 3.41660514E+00;
  COFETA[6] = -3.63206306E-01;
  COFETA[7] = 1.58671021E-02;
  COFETA[8] = -1.48001711E+01;
  COFETA[9] = 1.79491990E+00;
  COFETA[10] = -1.54008440E-01;
  COFETA[11] = 6.86719439E-03;
  COFETA[12] = -1.68118998E+01;
  COFETA[13] = 2.52362554E+00;
  COFETA[14] = -2.49309128E-01;
  COFETA[15] = 1.10211025E-02;
  COFETA[16] = -1.47696217E+01;
  COFETA[17] = 1.79491990E+00;
  COFETA[18] = -1.54008440E-01;
  COFETA[19] = 6.86719439E-03;
  COFETA[20] = -1.17770995E+01;
  COFETA[21] = -8.26744826E-01;
  COFETA[22] = 3.39009392E-01;
  COFETA[23] = -2.00674472E-02;
  COFETA[24] = -1.67963918E+01;
  COFETA[25] = 2.52362554E+00;
  COFETA[26] = -2.49309128E-01;
  COFETA[27] = 1.10211025E-02;
  COFETA[28] = -1.97154229E+01;
  COFETA[29] = 3.40193237E+00;
  COFETA[30] = -3.61470464E-01;
  COFETA[31] = 1.57990760E-02;
  COFETA[32] = -1.95453436E+01;
  COFETA[33] = 3.36385478E+00;
  COFETA[34] = -3.56948469E-01;
  COFETA[35] = 1.56210922E-02;
  COFETA[36] = -1.63031343E+01;
  COFETA[37] = 2.26143219E+00;
  COFETA[38] = -2.15114671E-01;
  COFETA[39] = 9.53461976E-03;
  COFETA[40] = -2.36749638E+01;
  COFETA[41] = 4.99775518E+00;
  COFETA[42] = -5.52687718E-01;
  COFETA[43] = 2.34353338E-02;
  COFETA[44] = -2.11136150E+01;
  COFETA[45] = 3.26961843E+00;
  COFETA[46] = -2.51355092E-01;
  COFETA[47] = 7.35605058E-03;
  COFETA[48] = -1.62526676E+01;
  COFETA[49] = 2.24839597E+00;
  COFETA[50] = -2.13428438E-01;
  COFETA[51] = 9.46192413E-03;
}

// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFLAM(amrex::Real* COFLAM)
{
  COFLAM[0] = 1.15898853E+01;
  COFLAM[1] = -1.52427526E+00;
  COFLAM[2] = 2.72840485E-01;
  COFLAM[3] = -1.03392506E-02;
  COFLAM[4] = -3.24554543E-01;
  COFLAM[5] = 3.41660514E+00;
  COFLAM[6] = -3.63206306E-01;
  COFLAM[7] = 1.58671021E-02;
  COFLAM[8] = 1.98515155E+00;
  COFLAM[9] = 1.79491990E+00;
  COFLAM[10] = -1.54008440E-01;
  COFLAM[11] = 6.86719439E-03;
  COFLAM[12] = -3.01283106E+00;
  COFLAM[13] = 3.37554996E+00;
  COFLAM[14] = -3.43353130E-01;
  COFLAM[15] = 1.51043449E-02;
  COFLAM[16] = 1.53490842E+01;
  COFLAM[17] = -3.77957911E+00;
  COFLAM[18] = 6.13516199E-01;
  COFLAM[19] = -2.72295607E-02;
  COFLAM[20] = 2.28195672E+01;
  COFLAM[21] = -8.72278727E+00;
  COFLAM[22] = 1.49300458E+00;
  COFLAM[23] = -7.41523911E-02;
  COFLAM[24] = 5.56033781E-01;
  COFLAM[25] = 1.59073629E+00;
  COFLAM[26] = -5.28054537E-02;
  COFLAM[27] = 4.07605156E-04;
  COFLAM[28] = 1.20629294E+01;
  COFLAM[29] = -3.78115762E+00;
  COFLAM[30] = 7.80383737E-01;
  COFLAM[31] = -4.14880733E-02;
  COFLAM[32] = 8.43904098E+00;
  COFLAM[33] = -2.78020359E+00;
  COFLAM[34] = 7.09313628E-01;
  COFLAM[35] = -4.04300952E-02;
  COFLAM[36] = 9.92460061E+00;
  COFLAM[37] = -2.28317864E+00;
  COFLAM[38] = 4.73113319E-01;
  COFLAM[39] = -2.40056449E-02;
  COFLAM[40] = -1.24047472E+01;
  COFLAM[41] = 6.34783072E+00;
  COFLAM[42] = -6.37857831E-01;
  COFLAM[43] = 2.37613806E-02;
  COFLAM[44] = 9.52224700E-02;
  COFLAM[45] = -1.75633398E-02;
  COFLAM[46] = 3.86704385E-01;
  COFLAM[47] = -2.84751875E-02;
  COFLAM[48] = 1.15506870E+01;
  COFLAM[49] = -2.91452034E+00;
  COFLAM[50] = 5.55043078E-01;
  COFLAM[51] = -2.75172223E-02;
}

// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFD(amrex::Real* COFD)
{
  COFD[0] = -1.02395385E+01;
  COFD[1] = 2.15403244E+00;
  COFD[2] = -6.97480266E-02;
  COFD[3] = 3.23666871E-03;
  COFD[4] = -1.11808845E+01;
  COFD[5] = 2.66936727E+00;
  COFD[6] = -1.34411514E-01;
  COFD[7] = 5.92957488E-03;
  COFD[8] = -1.06250314E+01;
  COFD[9] = 2.15849701E+00;
  COFD[10] = -6.53886401E-02;
  COFD[11] = 2.81453370E-03;
  COFD[12] = -1.15797897E+01;
  COFD[13] = 2.43235504E+00;
  COFD[14] = -1.02890179E-01;
  COFD[15] = 4.52903603E-03;
  COFD[16] = -1.06283588E+01;
  COFD[17] = 2.15849701E+00;
  COFD[18] = -6.53886401E-02;
  COFD[19] = 2.81453370E-03;
  COFD[20] = -1.68759075E+01;
  COFD[21] = 4.49460341E+00;
  COFD[22] = -3.64766179E-01;
  COFD[23] = 1.56457173E-02;
  COFD[24] = -1.15806955E+01;
  COFD[25] = 2.43235504E+00;
  COFD[26] = -1.02890179E-01;
  COFD[27] = 4.52903603E-03;
  COFD[28] = -1.22583533E+01;
  COFD[29] = 2.66488613E+00;
  COFD[30] = -1.33839930E-01;
  COFD[31] = 5.90524087E-03;
  COFD[32] = -1.22114409E+01;
  COFD[33] = 2.65215786E+00;
  COFD[34] = -1.32206766E-01;
  COFD[35] = 5.83532580E-03;
  COFD[36] = -1.13541221E+01;
  COFD[37] = 2.31999438E+00;
  COFD[38] = -8.75064804E-02;
  COFD[39] = 3.82656365E-03;
  COFD[40] = -1.35545391E+01;
  COFD[41] = 3.13878730E+00;
  COFD[42] = -1.94980335E-01;
  COFD[43] = 8.53744486E-03;
  COFD[44] = -1.57172268E+01;
  COFD[45] = 3.96062263E+00;
  COFD[46] = -2.98964970E-01;
  COFD[47] = 1.29322565E-02;
  COFD[48] = -1.13253618E+01;
  COFD[49] = 2.31195095E+00;
  COFD[50] = -8.63988037E-02;
  COFD[51] = 3.77573452E-03;
  COFD[52] = -1.11808845E+01;
  COFD[53] = 2.66936727E+00;
  COFD[54] = -1.34411514E-01;
  COFD[55] = 5.92957488E-03;
  COFD[56] = -1.43693219E+01;
  COFD[57] = 4.03992999E+00;
  COFD[58] = -3.08044800E-01;
  COFD[59] = 1.32757775E-02;
  COFD[60] = -1.31860264E+01;
  COFD[61] = 3.38003453E+00;
  COFD[62] = -2.25783856E-01;
  COFD[63] = 9.85028660E-03;
  COFD[64] = -1.43713019E+01;
  COFD[65] = 3.70920439E+00;
  COFD[66] = -2.67274113E-01;
  COFD[67] = 1.15967481E-02;
  COFD[68] = -1.31877859E+01;
  COFD[69] = 3.38003453E+00;
  COFD[70] = -2.25783856E-01;
  COFD[71] = 9.85028660E-03;
  COFD[72] = -1.93611201E+01;
  COFD[73] = 5.51579726E+00;
  COFD[74] = -4.76061961E-01;
  COFD[75] = 1.96329391E-02;
  COFD[76] = -1.43717684E+01;
  COFD[77] = 3.70920439E+00;
  COFD[78] = -2.67274113E-01;
  COFD[79] = 1.15967481E-02;
  COFD[80] = -1.53791907E+01;
  COFD[81] = 4.03496578E+00;
  COFD[82] = -3.07489899E-01;
  COFD[83] = 1.32554084E-02;
  COFD[84] = -1.53241298E+01;
  COFD[85] = 4.02249494E+00;
  COFD[86] = -3.06112211E-01;
  COFD[87] = 1.32056535E-02;
  COFD[88] = -1.40524219E+01;
  COFD[89] = 3.56261348E+00;
  COFD[90] = -2.48287981E-01;
  COFD[91] = 1.07752947E-02;
  COFD[92] = -1.72994129E+01;
  COFD[93] = 4.71931868E+00;
  COFD[94] = -3.91258152E-01;
  COFD[95] = 1.66866639E-02;
  COFD[96] = -1.95318517E+01;
  COFD[97] = 5.47046983E+00;
  COFD[98] = -4.74577605E-01;
  COFD[99] = 1.97408822E-02;
  COFD[100] = -1.40298991E+01;
  COFD[101] = 3.55837688E+00;
  COFD[102] = -2.47785790E-01;
  COFD[103] = 1.07555332E-02;
  COFD[104] = -1.06250314E+01;
  COFD[105] = 2.15849701E+00;
  COFD[106] = -6.53886401E-02;
  COFD[107] = 2.81453370E-03;
  COFD[108] = -1.31860264E+01;
  COFD[109] = 3.38003453E+00;
  COFD[110] = -2.25783856E-01;
  COFD[111] = 9.85028660E-03;
  COFD[112] = -1.29877255E+01;
  COFD[113] = 2.80841511E+00;
  COFD[114] = -1.52629888E-01;
  COFD[115] = 6.72604927E-03;
  COFD[116] = -1.40864784E+01;
  COFD[117] = 3.07458927E+00;
  COFD[118] = -1.86899591E-01;
  COFD[119] = 8.19829781E-03;
  COFD[120] = -1.30027669E+01;
  COFD[121] = 2.80841511E+00;
  COFD[122] = -1.52629888E-01;
  COFD[123] = 6.72604927E-03;
  COFD[124] = -1.91096700E+01;
  COFD[125] = 5.02608697E+00;
  COFD[126] = -4.26959993E-01;
  COFD[127] = 1.80709910E-02;
  COFD[128] = -1.40915944E+01;
  COFD[129] = 3.07458927E+00;
  COFD[130] = -1.86899591E-01;
  COFD[131] = 8.19829781E-03;
  COFD[132] = -1.48358026E+01;
  COFD[133] = 3.37160686E+00;
  COFD[134] = -2.24698307E-01;
  COFD[135] = 9.80357301E-03;
  COFD[136] = -1.47765660E+01;
  COFD[137] = 3.34934830E+00;
  COFD[138] = -2.21831026E-01;
  COFD[139] = 9.68019791E-03;
  COFD[140] = -1.39007309E+01;
  COFD[141] = 2.99164244E+00;
  COFD[142] = -1.76293106E-01;
  COFD[143] = 7.74575100E-03;
  COFD[144] = -1.67115471E+01;
  COFD[145] = 3.98859394E+00;
  COFD[146] = -3.02316219E-01;
  COFD[147] = 1.30661099E-02;
  COFD[148] = -1.91105534E+01;
  COFD[149] = 4.87977047E+00;
  COFD[150] = -4.10448693E-01;
  COFD[151] = 1.74535827E-02;
  COFD[152] = -1.38756380E+01;
  COFD[153] = 2.98558426E+00;
  COFD[154] = -1.75507216E-01;
  COFD[155] = 7.71173691E-03;
  COFD[156] = -1.15797897E+01;
  COFD[157] = 2.43235504E+00;
  COFD[158] = -1.02890179E-01;
  COFD[159] = 4.52903603E-03;
  COFD[160] = -1.43713019E+01;
  COFD[161] = 3.70920439E+00;
  COFD[162] = -2.67274113E-01;
  COFD[163] = 1.15967481E-02;
  COFD[164] = -1.40864784E+01;
  COFD[165] = 3.07458927E+00;
  COFD[166] = -1.86899591E-01;
  COFD[167] = 8.19829781E-03;
  COFD[168] = -1.53110597E+01;
  COFD[169] = 3.37317428E+00;
  COFD[170] = -2.24900439E-01;
  COFD[171] = 9.81228151E-03;
  COFD[172] = -1.41066358E+01;
  COFD[173] = 3.07458927E+00;
  COFD[174] = -1.86899591E-01;
  COFD[175] = 8.19829781E-03;
  COFD[176] = -2.10639926E+01;
  COFD[177] = 5.50980703E+00;
  COFD[178] = -4.78335492E-01;
  COFD[179] = 1.98515434E-02;
  COFD[180] = -1.53187536E+01;
  COFD[181] = 3.37317428E+00;
  COFD[182] = -2.24900439E-01;
  COFD[183] = 9.81228151E-03;
  COFD[184] = -1.60662884E+01;
  COFD[185] = 3.69744368E+00;
  COFD[186] = -2.65754489E-01;
  COFD[187] = 1.15311681E-02;
  COFD[188] = -1.59878735E+01;
  COFD[189] = 3.66478157E+00;
  COFD[190] = -2.61506432E-01;
  COFD[191] = 1.13466231E-02;
  COFD[192] = -1.50371688E+01;
  COFD[193] = 3.26249588E+00;
  COFD[194] = -2.10658287E-01;
  COFD[195] = 9.20032462E-03;
  COFD[196] = -1.81197250E+01;
  COFD[197] = 4.33684042E+00;
  COFD[198] = -3.44981265E-01;
  COFD[199] = 1.48142449E-02;
  COFD[200] = -2.05196205E+01;
  COFD[201] = 5.21473296E+00;
  COFD[202] = -4.48646311E-01;
  COFD[203] = 1.89013813E-02;
  COFD[204] = -1.50096252E+01;
  COFD[205] = 3.25515933E+00;
  COFD[206] = -2.09710110E-01;
  COFD[207] = 9.15941830E-03;
  COFD[208] = -1.06283588E+01;
  COFD[209] = 2.15849701E+00;
  COFD[210] = -6.53886401E-02;
  COFD[211] = 2.81453370E-03;
  COFD[212] = -1.31877859E+01;
  COFD[213] = 3.38003453E+00;
  COFD[214] = -2.25783856E-01;
  COFD[215] = 9.85028660E-03;
  COFD[216] = -1.30027669E+01;
  COFD[217] = 2.80841511E+00;
  COFD[218] = -1.52629888E-01;
  COFD[219] = 6.72604927E-03;
  COFD[220] = -1.41066358E+01;
  COFD[221] = 3.07458927E+00;
  COFD[222] = -1.86899591E-01;
  COFD[223] = 8.19829781E-03;
  COFD[224] = -1.30182749E+01;
  COFD[225] = 2.80841511E+00;
  COFD[226] = -1.52629888E-01;
  COFD[227] = 6.72604927E-03;
  COFD[228] = -1.91256173E+01;
  COFD[229] = 5.02608697E+00;
  COFD[230] = -4.26959993E-01;
  COFD[231] = 1.80709910E-02;
  COFD[232] = -1.41119635E+01;
  COFD[233] = 3.07458927E+00;
  COFD[234] = -1.86899591E-01;
  COFD[235] = 8.19829781E-03;
  COFD[236] = -1.48503699E+01;
  COFD[237] = 3.37160686E+00;
  COFD[238] = -2.24698307E-01;
  COFD[239] = 9.80357301E-03;
  COFD[240] = -1.47916284E+01;
  COFD[241] = 3.34934830E+00;
  COFD[242] = -2.21831026E-01;
  COFD[243] = 9.68019791E-03;
  COFD[244] = -1.39199573E+01;
  COFD[245] = 2.99164244E+00;
  COFD[246] = -1.76293106E-01;
  COFD[247] = 7.74575100E-03;
  COFD[248] = -1.67337674E+01;
  COFD[249] = 3.98859394E+00;
  COFD[250] = -3.02316219E-01;
  COFD[251] = 1.30661099E-02;
  COFD[252] = -1.91302705E+01;
  COFD[253] = 4.87977047E+00;
  COFD[254] = -4.10448693E-01;
  COFD[255] = 1.74535827E-02;
  COFD[256] = -1.38948654E+01;
  COFD[257] = 2.98558426E+00;
  COFD[258] = -1.75507216E-01;
  COFD[259] = 7.71173691E-03;
  COFD[260] = -1.68759075E+01;
  COFD[261] = 4.49460341E+00;
  COFD[262] = -3.64766179E-01;
  COFD[263] = 1.56457173E-02;
  COFD[264] = -1.93611201E+01;
  COFD[265] = 5.51579726E+00;
  COFD[266] = -4.76061961E-01;
  COFD[267] = 1.96329391E-02;
  COFD[268] = -1.91096700E+01;
  COFD[269] = 5.02608697E+00;
  COFD[270] = -4.26959993E-01;
  COFD[271] = 1.80709910E-02;
  COFD[272] = -2.10639926E+01;
  COFD[273] = 5.50980703E+00;
  COFD[274] = -4.78335492E-01;
  COFD[275] = 1.98515434E-02;
  COFD[276] = -1.91256173E+01;
  COFD[277] = 5.02608697E+00;
  COFD[278] = -4.26959993E-01;
  COFD[279] = 1.80709910E-02;
  COFD[280] = -1.31492509E+01;
  COFD[281] = 1.48004045E+00;
  COFD[282] = 1.60499956E-01;
  COFD[283] = -1.19765871E-02;
  COFD[284] = -2.04177394E+01;
  COFD[285] = 5.31457079E+00;
  COFD[286] = -4.58216496E-01;
  COFD[287] = 1.91825910E-02;
  COFD[288] = -2.10691956E+01;
  COFD[289] = 5.51604369E+00;
  COFD[290] = -4.76337283E-01;
  COFD[291] = 1.96544379E-02;
  COFD[292] = -2.12953626E+01;
  COFD[293] = 5.52385262E+00;
  COFD[294] = -4.69683817E-01;
  COFD[295] = 1.90677254E-02;
  COFD[296] = -2.08943723E+01;
  COFD[297] = 5.44718673E+00;
  COFD[298] = -4.72082974E-01;
  COFD[299] = 1.96531328E-02;
  COFD[300] = -2.12021420E+01;
  COFD[301] = 5.20775027E+00;
  COFD[302] = -4.07348285E-01;
  COFD[303] = 1.55473262E-02;
  COFD[304] = -1.87236049E+01;
  COFD[305] = 4.00967621E+00;
  COFD[306] = -2.21153539E-01;
  COFD[307] = 6.31528745E-03;
  COFD[308] = -2.08123331E+01;
  COFD[309] = 5.42470175E+00;
  COFD[310] = -4.69700438E-01;
  COFD[311] = 1.95706912E-02;
  COFD[312] = -1.15806955E+01;
  COFD[313] = 2.43235504E+00;
  COFD[314] = -1.02890179E-01;
  COFD[315] = 4.52903603E-03;
  COFD[316] = -1.43717684E+01;
  COFD[317] = 3.70920439E+00;
  COFD[318] = -2.67274113E-01;
  COFD[319] = 1.15967481E-02;
  COFD[320] = -1.40915944E+01;
  COFD[321] = 3.07458927E+00;
  COFD[322] = -1.86899591E-01;
  COFD[323] = 8.19829781E-03;
  COFD[324] = -1.53187536E+01;
  COFD[325] = 3.37317428E+00;
  COFD[326] = -2.24900439E-01;
  COFD[327] = 9.81228151E-03;
  COFD[328] = -1.41119635E+01;
  COFD[329] = 3.07458927E+00;
  COFD[330] = -1.86899591E-01;
  COFD[331] = 8.19829781E-03;
  COFD[332] = -2.04177394E+01;
  COFD[333] = 5.31457079E+00;
  COFD[334] = -4.58216496E-01;
  COFD[335] = 1.91825910E-02;
  COFD[336] = -1.53265677E+01;
  COFD[337] = 3.37317428E+00;
  COFD[338] = -2.24900439E-01;
  COFD[339] = 9.81228151E-03;
  COFD[340] = -1.60711937E+01;
  COFD[341] = 3.69744368E+00;
  COFD[342] = -2.65754489E-01;
  COFD[343] = 1.15311681E-02;
  COFD[344] = -1.59929990E+01;
  COFD[345] = 3.66478157E+00;
  COFD[346] = -2.61506432E-01;
  COFD[347] = 1.13466231E-02;
  COFD[348] = -1.50443476E+01;
  COFD[349] = 3.26249588E+00;
  COFD[350] = -2.10658287E-01;
  COFD[351] = 9.20032462E-03;
  COFD[352] = -1.81286456E+01;
  COFD[353] = 4.33684042E+00;
  COFD[354] = -3.44981265E-01;
  COFD[355] = 1.48142449E-02;
  COFD[356] = -2.05270679E+01;
  COFD[357] = 5.21473296E+00;
  COFD[358] = -4.48646311E-01;
  COFD[359] = 1.89013813E-02;
  COFD[360] = -1.50168047E+01;
  COFD[361] = 3.25515933E+00;
  COFD[362] = -2.09710110E-01;
  COFD[363] = 9.15941830E-03;
  COFD[364] = -1.22583533E+01;
  COFD[365] = 2.66488613E+00;
  COFD[366] = -1.33839930E-01;
  COFD[367] = 5.90524087E-03;
  COFD[368] = -1.53791907E+01;
  COFD[369] = 4.03496578E+00;
  COFD[370] = -3.07489899E-01;
  COFD[371] = 1.32554084E-02;
  COFD[372] = -1.48358026E+01;
  COFD[373] = 3.37160686E+00;
  COFD[374] = -2.24698307E-01;
  COFD[375] = 9.80357301E-03;
  COFD[376] = -1.60662884E+01;
  COFD[377] = 3.69744368E+00;
  COFD[378] = -2.65754489E-01;
  COFD[379] = 1.15311681E-02;
  COFD[380] = -1.48503699E+01;
  COFD[381] = 3.37160686E+00;
  COFD[382] = -2.24698307E-01;
  COFD[383] = 9.80357301E-03;
  COFD[384] = -2.10691956E+01;
  COFD[385] = 5.51604369E+00;
  COFD[386] = -4.76337283E-01;
  COFD[387] = 1.96544379E-02;
  COFD[388] = -1.60711937E+01;
  COFD[389] = 3.69744368E+00;
  COFD[390] = -2.65754489E-01;
  COFD[391] = 1.15311681E-02;
  COFD[392] = -1.69251551E+01;
  COFD[393] = 4.03028423E+00;
  COFD[394] = -3.06976249E-01;
  COFD[395] = 1.32370318E-02;
  COFD[396] = -1.68880017E+01;
  COFD[397] = 4.01759713E+00;
  COFD[398] = -3.05560618E-01;
  COFD[399] = 1.31852111E-02;
  COFD[400] = -1.57272832E+01;
  COFD[401] = 3.55756440E+00;
  COFD[402] = -2.47689003E-01;
  COFD[403] = 1.07517019E-02;
  COFD[404] = -1.90110170E+01;
  COFD[405] = 4.70842146E+00;
  COFD[406] = -3.89924637E-01;
  COFD[407] = 1.66320934E-02;
  COFD[408] = -2.12140281E+01;
  COFD[409] = 5.46309033E+00;
  COFD[410] = -4.73777184E-01;
  COFD[411] = 1.97122249E-02;
  COFD[412] = -1.57070757E+01;
  COFD[413] = 3.55327776E+00;
  COFD[414] = -2.47179250E-01;
  COFD[415] = 1.07315653E-02;
  COFD[416] = -1.22114409E+01;
  COFD[417] = 2.65215786E+00;
  COFD[418] = -1.32206766E-01;
  COFD[419] = 5.83532580E-03;
  COFD[420] = -1.53241298E+01;
  COFD[421] = 4.02249494E+00;
  COFD[422] = -3.06112211E-01;
  COFD[423] = 1.32056535E-02;
  COFD[424] = -1.47765660E+01;
  COFD[425] = 3.34934830E+00;
  COFD[426] = -2.21831026E-01;
  COFD[427] = 9.68019791E-03;
  COFD[428] = -1.59878735E+01;
  COFD[429] = 3.66478157E+00;
  COFD[430] = -2.61506432E-01;
  COFD[431] = 1.13466231E-02;
  COFD[432] = -1.47916284E+01;
  COFD[433] = 3.34934830E+00;
  COFD[434] = -2.21831026E-01;
  COFD[435] = 9.68019791E-03;
  COFD[436] = -2.12953626E+01;
  COFD[437] = 5.52385262E+00;
  COFD[438] = -4.69683817E-01;
  COFD[439] = 1.90677254E-02;
  COFD[440] = -1.59929990E+01;
  COFD[441] = 3.66478157E+00;
  COFD[442] = -2.61506432E-01;
  COFD[443] = 1.13466231E-02;
  COFD[444] = -1.68880017E+01;
  COFD[445] = 4.01759713E+00;
  COFD[446] = -3.05560618E-01;
  COFD[447] = 1.31852111E-02;
  COFD[448] = -1.68532873E+01;
  COFD[449] = 4.00572564E+00;
  COFD[450] = -3.04255586E-01;
  COFD[451] = 1.31384083E-02;
  COFD[452] = -1.56953785E+01;
  COFD[453] = 3.54443207E+00;
  COFD[454] = -2.46133003E-01;
  COFD[455] = 1.06905018E-02;
  COFD[456] = -1.89594762E+01;
  COFD[457] = 4.68735220E+00;
  COFD[458] = -3.87449514E-01;
  COFD[459] = 1.65352261E-02;
  COFD[460] = -2.11618589E+01;
  COFD[461] = 5.44451952E+00;
  COFD[462] = -4.71804032E-01;
  COFD[463] = 1.96436657E-02;
  COFD[464] = -1.56756124E+01;
  COFD[465] = 3.54035770E+00;
  COFD[466] = -2.45653442E-01;
  COFD[467] = 1.06717969E-02;
  COFD[468] = -1.13541221E+01;
  COFD[469] = 2.31999438E+00;
  COFD[470] = -8.75064804E-02;
  COFD[471] = 3.82656365E-03;
  COFD[472] = -1.40524219E+01;
  COFD[473] = 3.56261348E+00;
  COFD[474] = -2.48287981E-01;
  COFD[475] = 1.07752947E-02;
  COFD[476] = -1.39007309E+01;
  COFD[477] = 2.99164244E+00;
  COFD[478] = -1.76293106E-01;
  COFD[479] = 7.74575100E-03;
  COFD[480] = -1.50371688E+01;
  COFD[481] = 3.26249588E+00;
  COFD[482] = -2.10658287E-01;
  COFD[483] = 9.20032462E-03;
  COFD[484] = -1.39199573E+01;
  COFD[485] = 2.99164244E+00;
  COFD[486] = -1.76293106E-01;
  COFD[487] = 7.74575100E-03;
  COFD[488] = -2.08943723E+01;
  COFD[489] = 5.44718673E+00;
  COFD[490] = -4.72082974E-01;
  COFD[491] = 1.96531328E-02;
  COFD[492] = -1.50443476E+01;
  COFD[493] = 3.26249588E+00;
  COFD[494] = -2.10658287E-01;
  COFD[495] = 9.20032462E-03;
  COFD[496] = -1.57272832E+01;
  COFD[497] = 3.55756440E+00;
  COFD[498] = -2.47689003E-01;
  COFD[499] = 1.07517019E-02;
  COFD[500] = -1.56953785E+01;
  COFD[501] = 3.54443207E+00;
  COFD[502] = -2.46133003E-01;
  COFD[503] = 1.06905018E-02;
  COFD[504] = -1.48061406E+01;
  COFD[505] = 3.16912473E+00;
  COFD[506] = -1.98792456E-01;
  COFD[507] = 8.69726395E-03;
  COFD[508] = -1.77672912E+01;
  COFD[509] = 4.20234040E+00;
  COFD[510] = -3.28057658E-01;
  COFD[511] = 1.41006192E-02;
  COFD[512] = -2.02444913E+01;
  COFD[513] = 5.11785645E+00;
  COFD[514] = -4.37867828E-01;
  COFD[515] = 1.85047543E-02;
  COFD[516] = -1.47850505E+01;
  COFD[517] = 3.16433919E+00;
  COFD[518] = -1.98191564E-01;
  COFD[519] = 8.67209742E-03;
  COFD[520] = -1.35545391E+01;
  COFD[521] = 3.13878730E+00;
  COFD[522] = -1.94980335E-01;
  COFD[523] = 8.53744486E-03;
  COFD[524] = -1.72994129E+01;
  COFD[525] = 4.71931868E+00;
  COFD[526] = -3.91258152E-01;
  COFD[527] = 1.66866639E-02;
  COFD[528] = -1.67115471E+01;
  COFD[529] = 3.98859394E+00;
  COFD[530] = -3.02316219E-01;
  COFD[531] = 1.30661099E-02;
  COFD[532] = -1.81197250E+01;
  COFD[533] = 4.33684042E+00;
  COFD[534] = -3.44981265E-01;
  COFD[535] = 1.48142449E-02;
  COFD[536] = -1.67337674E+01;
  COFD[537] = 3.98859394E+00;
  COFD[538] = -3.02316219E-01;
  COFD[539] = 1.30661099E-02;
  COFD[540] = -2.12021420E+01;
  COFD[541] = 5.20775027E+00;
  COFD[542] = -4.07348285E-01;
  COFD[543] = 1.55473262E-02;
  COFD[544] = -1.81286456E+01;
  COFD[545] = 4.33684042E+00;
  COFD[546] = -3.44981265E-01;
  COFD[547] = 1.48142449E-02;
  COFD[548] = -1.90110170E+01;
  COFD[549] = 4.70842146E+00;
  COFD[550] = -3.89924637E-01;
  COFD[551] = 1.66320934E-02;
  COFD[552] = -1.89594762E+01;
  COFD[553] = 4.68735220E+00;
  COFD[554] = -3.87449514E-01;
  COFD[555] = 1.65352261E-02;
  COFD[556] = -1.77672912E+01;
  COFD[557] = 4.20234040E+00;
  COFD[558] = -3.28057658E-01;
  COFD[559] = 1.41006192E-02;
  COFD[560] = -2.10907633E+01;
  COFD[561] = 5.29211327E+00;
  COFD[562] = -4.56068366E-01;
  COFD[563] = 1.91195062E-02;
  COFD[564] = -2.20860863E+01;
  COFD[565] = 5.52171573E+00;
  COFD[566] = -4.63284984E-01;
  COFD[567] = 1.85570924E-02;
  COFD[568] = -1.77350630E+01;
  COFD[569] = 4.19328271E+00;
  COFD[570] = -3.26911461E-01;
  COFD[571] = 1.40520357E-02;
  COFD[572] = -1.57172268E+01;
  COFD[573] = 3.96062263E+00;
  COFD[574] = -2.98964970E-01;
  COFD[575] = 1.29322565E-02;
  COFD[576] = -1.95318517E+01;
  COFD[577] = 5.47046983E+00;
  COFD[578] = -4.74577605E-01;
  COFD[579] = 1.97408822E-02;
  COFD[580] = -1.91105534E+01;
  COFD[581] = 4.87977047E+00;
  COFD[582] = -4.10448693E-01;
  COFD[583] = 1.74535827E-02;
  COFD[584] = -2.05196205E+01;
  COFD[585] = 5.21473296E+00;
  COFD[586] = -4.48646311E-01;
  COFD[587] = 1.89013813E-02;
  COFD[588] = -1.91302705E+01;
  COFD[589] = 4.87977047E+00;
  COFD[590] = -4.10448693E-01;
  COFD[591] = 1.74535827E-02;
  COFD[592] = -1.87236049E+01;
  COFD[593] = 4.00967621E+00;
  COFD[594] = -2.21153539E-01;
  COFD[595] = 6.31528745E-03;
  COFD[596] = -2.05270679E+01;
  COFD[597] = 5.21473296E+00;
  COFD[598] = -4.48646311E-01;
  COFD[599] = 1.89013813E-02;
  COFD[600] = -2.12140281E+01;
  COFD[601] = 5.46309033E+00;
  COFD[602] = -4.73777184E-01;
  COFD[603] = 1.97122249E-02;
  COFD[604] = -2.11618589E+01;
  COFD[605] = 5.44451952E+00;
  COFD[606] = -4.71804032E-01;
  COFD[607] = 1.96436657E-02;
  COFD[608] = -2.02444913E+01;
  COFD[609] = 5.11785645E+00;
  COFD[610] = -4.37867828E-01;
  COFD[611] = 1.85047543E-02;
  COFD[612] = -2.20860863E+01;
  COFD[613] = 5.52171573E+00;
  COFD[614] = -4.63284984E-01;
  COFD[615] = 1.85570924E-02;
  COFD[616] = -1.99154422E+01;
  COFD[617] = 4.38041133E+00;
  COFD[618] = -2.77538214E-01;
  COFD[619] = 9.06748822E-03;
  COFD[620] = -2.02136057E+01;
  COFD[621] = 5.10993120E+00;
  COFD[622] = -4.36931630E-01;
  COFD[623] = 1.84677592E-02;
  COFD[624] = -1.13253618E+01;
  COFD[625] = 2.31195095E+00;
  COFD[626] = -8.63988037E-02;
  COFD[627] = 3.77573452E-03;
  COFD[628] = -1.40298991E+01;
  COFD[629] = 3.55837688E+00;
  COFD[630] = -2.47785790E-01;
  COFD[631] = 1.07555332E-02;
  COFD[632] = -1.38756380E+01;
  COFD[633] = 2.98558426E+00;
  COFD[634] = -1.75507216E-01;
  COFD[635] = 7.71173691E-03;
  COFD[636] = -1.50096252E+01;
  COFD[637] = 3.25515933E+00;
  COFD[638] = -2.09710110E-01;
  COFD[639] = 9.15941830E-03;
  COFD[640] = -1.38948654E+01;
  COFD[641] = 2.98558426E+00;
  COFD[642] = -1.75507216E-01;
  COFD[643] = 7.71173691E-03;
  COFD[644] = -2.08123331E+01;
  COFD[645] = 5.42470175E+00;
  COFD[646] = -4.69700438E-01;
  COFD[647] = 1.95706912E-02;
  COFD[648] = -1.50168047E+01;
  COFD[649] = 3.25515933E+00;
  COFD[650] = -2.09710110E-01;
  COFD[651] = 9.15941830E-03;
  COFD[652] = -1.57070757E+01;
  COFD[653] = 3.55327776E+00;
  COFD[654] = -2.47179250E-01;
  COFD[655] = 1.07315653E-02;
  COFD[656] = -1.56756124E+01;
  COFD[657] = 3.54035770E+00;
  COFD[658] = -2.45653442E-01;
  COFD[659] = 1.06717969E-02;
  COFD[660] = -1.47850505E+01;
  COFD[661] = 3.16433919E+00;
  COFD[662] = -1.98191564E-01;
  COFD[663] = 8.67209742E-03;
  COFD[664] = -1.77350630E+01;
  COFD[665] = 4.19328271E+00;
  COFD[666] = -3.26911461E-01;
  COFD[667] = 1.40520357E-02;
  COFD[668] = -2.02136057E+01;
  COFD[669] = 5.10993120E+00;
  COFD[670] = -4.36931630E-01;
  COFD[671] = 1.84677592E-02;
  COFD[672] = -1.47639411E+01;
  COFD[673] = 3.15955654E+00;
  COFD[674] = -1.97590757E-01;
  COFD[675] = 8.64692156E-03;
}

// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKTDIF(int* KTDIF)
{
  KTDIF[0] = 0;
  KTDIF[1] = 1;
}

// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFTD(amrex::Real* COFTD)
{
  COFTD[0] = 0.00000000E+00;
  COFTD[1] = 0.00000000E+00;
  COFTD[2] = 0.00000000E+00;
  COFTD[3] = 0.00000000E+00;
  COFTD[4] = -1.52534742E-01;
  COFTD[5] = -5.46404022E-05;
  COFTD[6] = 2.93412470E-08;
  COFTD[7] = -4.87091914E-12;
  COFTD[8] = 4.15577509E-01;
  COFTD[9] = 1.09736860E-05;
  COFTD[10] = -3.96016410E-09;
  COFTD[11] = 1.14412838E-12;
  COFTD[12] = 4.42736017E-01;
  COFTD[13] = 7.11765887E-05;
  COFTD[14] = -3.84765396E-08;
  COFTD[15] = 6.86318682E-12;
  COFTD[16] = 4.21927217E-01;
  COFTD[17] = 1.11413555E-05;
  COFTD[18] = -4.02067239E-09;
  COFTD[19] = 1.16160979E-12;
  COFTD[20] = 6.02021050E-02;
  COFTD[21] = 5.61555733E-04;
  COFTD[22] = -2.55370065E-07;
  COFTD[23] = 3.63385926E-11;
  COFTD[24] = 4.44449675E-01;
  COFTD[25] = 7.14520855E-05;
  COFTD[26] = -3.86254673E-08;
  COFTD[27] = 6.88975155E-12;
  COFTD[28] = 3.50272421E-01;
  COFTD[29] = 1.23583915E-04;
  COFTD[30] = -6.63969060E-08;
  COFTD[31] = 1.10354999E-11;
  COFTD[32] = 3.58655012E-01;
  COFTD[33] = 1.21524029E-04;
  COFTD[34] = -6.53752231E-08;
  COFTD[35] = 1.09001717E-11;
  COFTD[36] = 4.44650439E-01;
  COFTD[37] = 5.06628083E-05;
  COFTD[38] = -2.69818971E-08;
  COFTD[39] = 5.01286176E-12;
  COFTD[40] = 3.25740914E-01;
  COFTD[41] = 3.03631979E-04;
  COFTD[42] = -1.55289598E-07;
  COFTD[43] = 2.41465297E-11;
  COFTD[44] = 1.62149345E-01;
  COFTD[45] = 5.03434613E-04;
  COFTD[46] = -2.39391657E-07;
  COFTD[47] = 3.50983128E-11;
  COFTD[48] = 4.45261429E-01;
  COFTD[49] = 4.94696577E-05;
  COFTD[50] = -2.63023124E-08;
  COFTD[51] = 4.90305626E-12;
  COFTD[52] = 1.52534742E-01;
  COFTD[53] = 5.46404022E-05;
  COFTD[54] = -2.93412470E-08;
  COFTD[55] = 4.87091914E-12;
  COFTD[56] = 0.00000000E+00;
  COFTD[57] = 0.00000000E+00;
  COFTD[58] = 0.00000000E+00;
  COFTD[59] = 0.00000000E+00;
  COFTD[60] = 2.70008280E-01;
  COFTD[61] = 3.61552588E-04;
  COFTD[62] = -1.80743499E-07;
  COFTD[63] = 2.75319341E-11;
  COFTD[64] = 2.20482081E-01;
  COFTD[65] = 4.80162630E-04;
  COFTD[66] = -2.32927140E-07;
  COFTD[67] = 3.46469239E-11;
  COFTD[68] = 2.72039997E-01;
  COFTD[69] = 3.64273144E-04;
  COFTD[70] = -1.82103531E-07;
  COFTD[71] = 2.77391022E-11;
  COFTD[72] = -1.41882969E-01;
  COFTD[73] = 7.66554625E-04;
  COFTD[74] = -3.06548330E-07;
  COFTD[75] = 4.02957302E-11;
  COFTD[76] = 2.20907136E-01;
  COFTD[77] = 4.81088308E-04;
  COFTD[78] = -2.33376187E-07;
  COFTD[79] = 3.47137178E-11;
  COFTD[80] = 1.40109849E-01;
  COFTD[81] = 5.30876682E-04;
  COFTD[82] = -2.49898194E-07;
  COFTD[83] = 3.63822692E-11;
  COFTD[84] = 1.45433118E-01;
  COFTD[85] = 5.30347483E-04;
  COFTD[86] = -2.50118653E-07;
  COFTD[87] = 3.64613999E-11;
  COFTD[88] = 2.39409087E-01;
  COFTD[89] = 4.47195587E-04;
  COFTD[90] = -2.18950923E-07;
  COFTD[91] = 3.27972343E-11;
  COFTD[92] = 2.44368810E-02;
  COFTD[93] = 7.18240808E-04;
  COFTD[94] = -3.19717751E-07;
  COFTD[95] = 4.48827629E-11;
  COFTD[96] = -1.24938974E-01;
  COFTD[97] = 7.98385053E-04;
  COFTD[98] = -3.25757471E-07;
  COFTD[99] = 4.33526419E-11;
  COFTD[100] = 2.40744276E-01;
  COFTD[101] = 4.45343183E-04;
  COFTD[102] = -2.18173742E-07;
  COFTD[103] = 3.26958310E-11;
}

// compute the critical parameters for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS(
  amrex::Real* Tci, amrex::Real* ai, amrex::Real* bi, amrex::Real* acentric_i)
{

  amrex::Real EPS[13];
  amrex::Real SIG[13];
  amrex::Real wt[13];
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

  // species 1: H
  Tci[1] = 1.316 * EPS[1];
  ai[1] = (5.55 * avogadro * avogadro * EPS[1] * boltzmann * 1e-24 * SIG[1] *
           SIG[1] * SIG[1]) /
          (wt[1] * wt[1]);
  bi[1] = 0.855 * avogadro * 1e-24 * SIG[1] * SIG[1] * SIG[1] / (wt[1]);
  acentric_i[1] = 0.0;

  // species 2: O
  Tci[2] = 1.316 * EPS[2];
  ai[2] = (5.55 * avogadro * avogadro * EPS[2] * boltzmann * 1e-24 * SIG[2] *
           SIG[2] * SIG[2]) /
          (wt[2] * wt[2]);
  bi[2] = 0.855 * avogadro * 1e-24 * SIG[2] * SIG[2] * SIG[2] / (wt[2]);
  acentric_i[2] = 0.0;

  // species 3: O2
  // Imported from NIST
  Tci[3] = 154.581000;
  ai[3] = 1e6 * 0.42748 * Rcst * Rcst * Tci[3] * Tci[3] /
          (31.998800 * 31.998800 * 50.430466);
  bi[3] = 0.08664 * Rcst * Tci[3] / (31.998800 * 50.430466);
  acentric_i[3] = 0.022200;

  // species 4: OH
  Tci[4] = 1.316 * EPS[4];
  ai[4] = (5.55 * avogadro * avogadro * EPS[4] * boltzmann * 1e-24 * SIG[4] *
           SIG[4] * SIG[4]) /
          (wt[4] * wt[4]);
  bi[4] = 0.855 * avogadro * 1e-24 * SIG[4] * SIG[4] * SIG[4] / (wt[4]);
  acentric_i[4] = 0.0;

  // species 5: H2O
  // Imported from NIST
  Tci[5] = 647.096000;
  ai[5] = 1e6 * 0.42748 * Rcst * Rcst * Tci[5] * Tci[5] /
          (18.015340 * 18.015340 * 220.640000);
  bi[5] = 0.08664 * Rcst * Tci[5] / (18.015340 * 220.640000);
  acentric_i[5] = 0.344300;

  // species 6: HO2
  Tci[6] = 1.316 * EPS[6];
  ai[6] = (5.55 * avogadro * avogadro * EPS[6] * boltzmann * 1e-24 * SIG[6] *
           SIG[6] * SIG[6]) /
          (wt[6] * wt[6]);
  bi[6] = 0.855 * avogadro * 1e-24 * SIG[6] * SIG[6] * SIG[6] / (wt[6]);
  acentric_i[6] = 0.0;

  // species 7: CH3
  Tci[7] = 1.316 * EPS[7];
  ai[7] = (5.55 * avogadro * avogadro * EPS[7] * boltzmann * 1e-24 * SIG[7] *
           SIG[7] * SIG[7]) /
          (wt[7] * wt[7]);
  bi[7] = 0.855 * avogadro * 1e-24 * SIG[7] * SIG[7] * SIG[7] / (wt[7]);
  acentric_i[7] = 0.0;

  // species 8: CH4
  // Imported from NIST
  Tci[8] = 190.560000;
  ai[8] = 1e6 * 0.42748 * Rcst * Rcst * Tci[8] * Tci[8] /
          (16.043030 * 16.043030 * 45.990000);
  bi[8] = 0.08664 * Rcst * Tci[8] / (16.043030 * 45.990000);
  acentric_i[8] = 0.011000;

  // species 9: CO
  // Imported from NIST
  Tci[9] = 132.850000;
  ai[9] = 1e6 * 0.42748 * Rcst * Rcst * Tci[9] * Tci[9] /
          (28.010000 * 28.010000 * 34.940000);
  bi[9] = 0.08664 * Rcst * Tci[9] / (28.010000 * 34.940000);
  acentric_i[9] = 0.045000;

  // species 10: CO2
  // Imported from NIST
  Tci[10] = 304.120000;
  ai[10] = 1e6 * 0.42748 * Rcst * Rcst * Tci[10] * Tci[10] /
           (44.009950 * 44.009950 * 73.740000);
  bi[10] = 0.08664 * Rcst * Tci[10] / (44.009950 * 73.740000);
  acentric_i[10] = 0.225000;

  // species 11: CH2O
  Tci[11] = 1.316 * EPS[11];
  ai[11] = (5.55 * avogadro * avogadro * EPS[11] * boltzmann * 1e-24 * SIG[11] *
            SIG[11] * SIG[11]) /
           (wt[11] * wt[11]);
  bi[11] = 0.855 * avogadro * 1e-24 * SIG[11] * SIG[11] * SIG[11] / (wt[11]);
  acentric_i[11] = 0.0;

  // species 12: N2
  // Imported from NIST
  Tci[12] = 126.192000;
  ai[12] = 1e6 * 0.42748 * Rcst * Rcst * Tci[12] * Tci[12] /
           (28.013400 * 28.013400 * 33.958000);
  bi[12] = 0.08664 * Rcst * Tci[12] / (28.013400 * 33.958000);
  acentric_i[12] = 0.037200;
}

// compute the critical parameter quantities for each species for SRK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS_SRK(
  amrex::Real* sqrtOneOverTc,
  amrex::Real* sqrtAsti,
  amrex::Real* Bi,
  amrex::Real* Fomega)
{

  // species 0: H2
  // Imported from NIST
  sqrtOneOverTc[0] = 1.7369646834812e-01;
  sqrtAsti[0] = 2.4824229305911e+05;
  Bi[0] = 9.1362254169923e+00;
  Fomega[0] = 1.3798618890700e-01;

  // species 1: H
  // Computed from Lennard-Jones
  sqrtOneOverTc[1] = 7.2391579785311e-02;
  sqrtAsti[1] = 5.8451221185793e+05;
  Bi[1] = 4.4006634563410e+00;
  Fomega[1] = 4.8508000000000e-01;

  // species 2: O
  // Computed from Lennard-Jones
  sqrtOneOverTc[2] = 9.7460146948181e-02;
  sqrtAsti[2] = 4.2500156012386e+04;
  Bi[2] = 6.6930234851885e-01;
  Fomega[2] = 4.8508000000000e-01;

  // species 3: O2
  // Imported from NIST
  sqrtOneOverTc[3] = 8.0430717653241e-02;
  sqrtAsti[3] = 3.6980081924740e+04;
  Bi[3] = 6.9005216444999e-01;
  Fomega[3] = 5.1945301904908e-01;

  // species 4: OH
  // Computed from Lennard-Jones
  sqrtOneOverTc[4] = 9.7460146948181e-02;
  sqrtAsti[4] = 3.9981183985545e+04;
  Bi[4] = 6.2963299076575e-01;
  Fomega[4] = 4.8508000000000e-01;

  // species 5: H2O
  // Imported from NIST
  sqrtOneOverTc[5] = 3.9311140369715e-02;
  sqrtAsti[5] = 1.3145468419592e+05;
  Bi[5] = 1.1727204136223e+00;
  Fomega[5] = 1.0013577274636e+00;

  // species 6: HO2
  // Computed from Lennard-Jones
  sqrtOneOverTc[6] = 8.4114315417892e-02;
  sqrtAsti[6] = 3.3657825710980e+04;
  Bi[6] = 6.4505861132042e-01;
  Fomega[6] = 4.8508000000000e-01;

  // species 7: CH3
  // Computed from Lennard-Jones
  sqrtOneOverTc[7] = 7.2642504557751e-02;
  sqrtAsti[7] = 9.8558148656101e+04;
  Bi[7] = 1.8791630452967e+00;
  Fomega[7] = 4.8508000000000e-01;

  // species 8: CH4
  // Imported from NIST
  sqrtOneOverTc[8] = 7.2440948474968e-02;
  sqrtAsti[8] = 9.5215046085708e+04;
  Bi[8] = 1.8605203320660e+00;
  Fomega[8] = 5.0213035482700e-01;

  // species 9: CO
  // Imported from NIST
  sqrtOneOverTc[9] = 8.6759935530451e-02;
  sqrtAsti[9] = 4.3619324712923e+04;
  Bi[9] = 9.7786303823350e-01;
  Fomega[9] = 5.5459948367500e-01;

  // species 10: CO2
  // Imported from NIST
  sqrtOneOverTc[10] = 5.7342616962522e-02;
  sqrtAsti[10] = 4.3745610142919e+04;
  Bi[10] = 6.7506167841745e-01;
  Fomega[10] = 8.2653709187500e-01;

  // species 11: CH2O
  // Computed from Lennard-Jones
  sqrtOneOverTc[11] = 3.9062261583604e-02;
  sqrtAsti[11] = 8.4274889632487e+04;
  Bi[11] = 7.9341968505345e-01;
  Fomega[11] = 4.8508000000000e-01;

  // species 12: N2
  // Imported from NIST
  sqrtOneOverTc[12] = 8.9019282240563e-02;
  sqrtAsti[12] = 4.2022983929562e+04;
  Bi[12] = 9.5560052221563e-01;
  Fomega[12] = 5.4259343186608e-01;
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
