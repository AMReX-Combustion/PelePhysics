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

#define NUM_GAS_ELEMENTS 4      // Elements in the homogeneous phase
#define NUM_QSSA_GAS_SPECIES 13 // Species in the homogeneous phase
#define NUM_GAS_REACTIONS 73    // Reactions in the homogeneous phase

#define SITE_DENSITY 0.000000E+00 // mol/cm^2

#define NUM_SURFACE_ELEMENTS 0  // Additional elements in heterogeneous phase
#define NUM_SURFACE_SPECIES 0   // Species in the heterogeneous phase
#define NUM_SURFACE_REACTIONS 0 // Reactions in the heterogeneous phase

#define NUM_ELEMENTS (NUM_GAS_ELEMENTS + NUM_SURFACE_ELEMENTS)
#define NUM_SPECIES (NUM_QSSA_GAS_SPECIES + NUM_SURFACE_SPECIES)
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
  mm = 4;
  kk = 13;
  ii = 73;
  nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[13] = {
  4.9603174603174605e-01, // H2
  9.9206349206349209e-01, // H
  6.2503906494155881e-02, // O
  3.1251953247077940e-02, // O2
  5.8799317927912029e-02, // OH
  5.5509297807382736e-02, // H2O
  3.0297521662727988e-02, // HO2
  6.6511473229132020e-02, // CH3
  6.2332481456086768e-02, // CH4
  3.5701535166012141e-02, // CO
  2.2722624917630486e-02, // CO2
  3.3304469459801506e-02, // CH2O
  3.5696437495537948e-02, // N2
};
#endif
const amrex::Real h_global_imw[13] = {
  4.9603174603174605e-01, // H2
  9.9206349206349209e-01, // H
  6.2503906494155881e-02, // O
  3.1251953247077940e-02, // O2
  5.8799317927912029e-02, // OH
  5.5509297807382736e-02, // H2O
  3.0297521662727988e-02, // HO2
  6.6511473229132020e-02, // CH3
  6.2332481456086768e-02, // CH4
  3.5701535166012141e-02, // CO
  2.2722624917630486e-02, // CO2
  3.3304469459801506e-02, // CH2O
  3.5696437495537948e-02, // N2
};

//  molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[13] = {
  2.0160000000000000e+00, // H2
  1.0080000000000000e+00, // H
  1.5999000000000001e+01, // O
  3.1998000000000001e+01, // O2
  1.7007000000000001e+01, // OH
  1.8015000000000001e+01, // H2O
  3.3006000000000000e+01, // HO2
  1.5035000000000000e+01, // CH3
  1.6042999999999999e+01, // CH4
  2.8009999999999998e+01, // CO
  4.4009000000000000e+01, // CO2
  3.0026000000000000e+01, // CH2O
  2.8013999999999999e+01, // N2
};
#endif
const amrex::Real h_global_mw[13] = {
  2.0160000000000000e+00, // H2
  1.0080000000000000e+00, // H
  1.5999000000000001e+01, // O
  3.1998000000000001e+01, // O2
  1.7007000000000001e+01, // OH
  1.8015000000000001e+01, // H2O
  3.3006000000000000e+01, // HO2
  1.5035000000000000e+01, // CH3
  1.6042999999999999e+01, // CH4
  2.8009999999999998e+01, // CO
  4.4009000000000000e+01, // CO2
  3.0026000000000000e+01, // CH2O
  2.8013999999999999e+01, // N2
};

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_imw(amrex::Real* imw_new)
{
  imw_new[0] = 4.9603174603174605e-01;  // H2
  imw_new[1] = 9.9206349206349209e-01;  // H
  imw_new[2] = 6.2503906494155881e-02;  // O
  imw_new[3] = 3.1251953247077940e-02;  // O2
  imw_new[4] = 5.8799317927912029e-02;  // OH
  imw_new[5] = 5.5509297807382736e-02;  // H2O
  imw_new[6] = 3.0297521662727988e-02;  // HO2
  imw_new[7] = 6.6511473229132020e-02;  // CH3
  imw_new[8] = 6.2332481456086768e-02;  // CH4
  imw_new[9] = 3.5701535166012141e-02;  // CO
  imw_new[10] = 2.2722624917630486e-02; // CO2
  imw_new[11] = 3.3304469459801506e-02; // CH2O
  imw_new[12] = 3.5696437495537948e-02; // N2
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
  mw_new[0] = 2.0160000000000000e+00;  // H2
  mw_new[1] = 1.0080000000000000e+00;  // H
  mw_new[2] = 1.5999000000000001e+01;  // O
  mw_new[3] = 3.1998000000000001e+01;  // O2
  mw_new[4] = 1.7007000000000001e+01;  // OH
  mw_new[5] = 1.8015000000000001e+01;  // H2O
  mw_new[6] = 3.3006000000000000e+01;  // HO2
  mw_new[7] = 1.5035000000000000e+01;  // CH3
  mw_new[8] = 1.6042999999999999e+01;  // CH4
  mw_new[9] = 2.8009999999999998e+01;  // CO
  mw_new[10] = 4.4009000000000000e+01; // CO2
  mw_new[11] = 3.0026000000000000e+01; // CH2O
  mw_new[12] = 2.8013999999999999e+01; // N2
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;

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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
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
  const amrex::Real T3 = T * T2;

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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;

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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;

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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
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
  const amrex::Real T3 = T * T2;
  const amrex::Real T4 = T * T3;
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
  x[0] = y[0] / (2.0160000000000000e+00 * YOW);
  x[1] = y[1] / (1.0080000000000000e+00 * YOW);
  x[2] = y[2] / (1.5999000000000001e+01 * YOW);
  x[3] = y[3] / (3.1998000000000001e+01 * YOW);
  x[4] = y[4] / (1.7007000000000001e+01 * YOW);
  x[5] = y[5] / (1.8015000000000001e+01 * YOW);
  x[6] = y[6] / (3.3006000000000000e+01 * YOW);
  x[7] = y[7] / (1.5035000000000000e+01 * YOW);
  x[8] = y[8] / (1.6042999999999999e+01 * YOW);
  x[9] = y[9] / (2.8009999999999998e+01 * YOW);
  x[10] = y[10] / (4.4009000000000000e+01 * YOW);
  x[11] = y[11] / (3.0026000000000000e+01 * YOW);
  x[12] = y[12] / (2.8013999999999999e+01 * YOW);
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
  amrex::Real tmin = 300.0;  // max lower bound for thermo def
  amrex::Real tmax = 3500.0; // min upper bound for thermo def
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
  amrex::Real tmax = 3500.0; // min upper bound for thermo def
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
  XW += x[0] * 2.0160000000000000e+00;     // H2
  XW += x[1] * 1.0080000000000000e+00;     // H
  XW += x[2] * 1.5999000000000001e+01;     // O
  XW += x[3] * 3.1998000000000001e+01;     // O2
  XW += x[4] * 1.7007000000000001e+01;     // OH
  XW += x[5] * 1.8015000000000001e+01;     // H2O
  XW += x[6] * 3.3006000000000000e+01;     // HO2
  XW += x[7] * 1.5035000000000000e+01;     // CH3
  XW += x[8] * 1.6042999999999999e+01;     // CH4
  XW += x[9] * 2.8009999999999998e+01;     // CO
  XW += x[10] * 4.4009000000000000e+01;    // CO2
  XW += x[11] * 3.0026000000000000e+01;    // CH2O
  XW += x[12] * 2.8013999999999999e+01;    // N2
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
  W += c[0] * 2.0160000000000000e+00;  // H2
  W += c[1] * 1.0080000000000000e+00;  // H
  W += c[2] * 1.5999000000000001e+01;  // O
  W += c[3] * 3.1998000000000001e+01;  // O2
  W += c[4] * 1.7007000000000001e+01;  // OH
  W += c[5] * 1.8015000000000001e+01;  // H2O
  W += c[6] * 3.3006000000000000e+01;  // HO2
  W += c[7] * 1.5035000000000000e+01;  // CH3
  W += c[8] * 1.6042999999999999e+01;  // CH4
  W += c[9] * 2.8009999999999998e+01;  // CO
  W += c[10] * 4.4009000000000000e+01; // CO2
  W += c[11] * 3.0026000000000000e+01; // CH2O
  W += c[12] * 2.8013999999999999e+01; // N2

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
  XW += x[0] * 2.0160000000000000e+00;       // H2
  XW += x[1] * 1.0080000000000000e+00;       // H
  XW += x[2] * 1.5999000000000001e+01;       // O
  XW += x[3] * 3.1998000000000001e+01;       // O2
  XW += x[4] * 1.7007000000000001e+01;       // OH
  XW += x[5] * 1.8015000000000001e+01;       // H2O
  XW += x[6] * 3.3006000000000000e+01;       // HO2
  XW += x[7] * 1.5035000000000000e+01;       // CH3
  XW += x[8] * 1.6042999999999999e+01;       // CH4
  XW += x[9] * 2.8009999999999998e+01;       // CO
  XW += x[10] * 4.4009000000000000e+01;      // CO2
  XW += x[11] * 3.0026000000000000e+01;      // CH2O
  XW += x[12] * 2.8013999999999999e+01;      // N2
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
  W += c[0] * 2.0160000000000000e+00;  // H2
  W += c[1] * 1.0080000000000000e+00;  // H
  W += c[2] * 1.5999000000000001e+01;  // O
  W += c[3] * 3.1998000000000001e+01;  // O2
  W += c[4] * 1.7007000000000001e+01;  // OH
  W += c[5] * 1.8015000000000001e+01;  // H2O
  W += c[6] * 3.3006000000000000e+01;  // HO2
  W += c[7] * 1.5035000000000000e+01;  // CH3
  W += c[8] * 1.6042999999999999e+01;  // CH4
  W += c[9] * 2.8009999999999998e+01;  // CO
  W += c[10] * 4.4009000000000000e+01; // CO2
  W += c[11] * 3.0026000000000000e+01; // CH2O
  W += c[12] * 2.8013999999999999e+01; // N2

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
  amrex::Real XW = 0;                   //  see Eq 4 in CK Manual
  XW += x[0] * 2.0160000000000000e+00;  // H2
  XW += x[1] * 1.0080000000000000e+00;  // H
  XW += x[2] * 1.5999000000000001e+01;  // O
  XW += x[3] * 3.1998000000000001e+01;  // O2
  XW += x[4] * 1.7007000000000001e+01;  // OH
  XW += x[5] * 1.8015000000000001e+01;  // H2O
  XW += x[6] * 3.3006000000000000e+01;  // HO2
  XW += x[7] * 1.5035000000000000e+01;  // CH3
  XW += x[8] * 1.6042999999999999e+01;  // CH4
  XW += x[9] * 2.8009999999999998e+01;  // CO
  XW += x[10] * 4.4009000000000000e+01; // CO2
  XW += x[11] * 3.0026000000000000e+01; // CH2O
  XW += x[12] * 2.8013999999999999e+01; // N2
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
  W += c[0] * 2.0160000000000000e+00;  // H2
  W += c[1] * 1.0080000000000000e+00;  // H
  W += c[2] * 1.5999000000000001e+01;  // O
  W += c[3] * 3.1998000000000001e+01;  // O2
  W += c[4] * 1.7007000000000001e+01;  // OH
  W += c[5] * 1.8015000000000001e+01;  // H2O
  W += c[6] * 3.3006000000000000e+01;  // HO2
  W += c[7] * 1.5035000000000000e+01;  // CH3
  W += c[8] * 1.6042999999999999e+01;  // CH4
  W += c[9] * 2.8009999999999998e+01;  // CO
  W += c[10] * 4.4009000000000000e+01; // CO2
  W += c[11] * 3.0026000000000000e+01; // CH2O
  W += c[12] * 2.8013999999999999e+01; // N2

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
  XW += x[0] * 2.0160000000000000e+00;  // H2
  XW += x[1] * 1.0080000000000000e+00;  // H
  XW += x[2] * 1.5999000000000001e+01;  // O
  XW += x[3] * 3.1998000000000001e+01;  // O2
  XW += x[4] * 1.7007000000000001e+01;  // OH
  XW += x[5] * 1.8015000000000001e+01;  // H2O
  XW += x[6] * 3.3006000000000000e+01;  // HO2
  XW += x[7] * 1.5035000000000000e+01;  // CH3
  XW += x[8] * 1.6042999999999999e+01;  // CH4
  XW += x[9] * 2.8009999999999998e+01;  // CO
  XW += x[10] * 4.4009000000000000e+01; // CO2
  XW += x[11] * 3.0026000000000000e+01; // CH2O
  XW += x[12] * 2.8013999999999999e+01; // N2
  // Now compute conversion
  amrex::Real XWinv = 1.0 / XW;
  y[0] = x[0] * 2.0160000000000000e+00 * XWinv;
  y[1] = x[1] * 1.0080000000000000e+00 * XWinv;
  y[2] = x[2] * 1.5999000000000001e+01 * XWinv;
  y[3] = x[3] * 3.1998000000000001e+01 * XWinv;
  y[4] = x[4] * 1.7007000000000001e+01 * XWinv;
  y[5] = x[5] * 1.8015000000000001e+01 * XWinv;
  y[6] = x[6] * 3.3006000000000000e+01 * XWinv;
  y[7] = x[7] * 1.5035000000000000e+01 * XWinv;
  y[8] = x[8] * 1.6042999999999999e+01 * XWinv;
  y[9] = x[9] * 2.8009999999999998e+01 * XWinv;
  y[10] = x[10] * 4.4009000000000000e+01 * XWinv;
  y[11] = x[11] * 3.0026000000000000e+01 * XWinv;
  y[12] = x[12] * 2.8013999999999999e+01 * XWinv;
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
  XW += x[0] * 2.0160000000000000e+00;  // H2
  XW += x[1] * 1.0080000000000000e+00;  // H
  XW += x[2] * 1.5999000000000001e+01;  // O
  XW += x[3] * 3.1998000000000001e+01;  // O2
  XW += x[4] * 1.7007000000000001e+01;  // OH
  XW += x[5] * 1.8015000000000001e+01;  // H2O
  XW += x[6] * 3.3006000000000000e+01;  // HO2
  XW += x[7] * 1.5035000000000000e+01;  // CH3
  XW += x[8] * 1.6042999999999999e+01;  // CH4
  XW += x[9] * 2.8009999999999998e+01;  // CO
  XW += x[10] * 4.4009000000000000e+01; // CO2
  XW += x[11] * 3.0026000000000000e+01; // CH2O
  XW += x[12] * 2.8013999999999999e+01; // N2
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
  CW += c[0] * 2.0160000000000000e+00;  // H2
  CW += c[1] * 1.0080000000000000e+00;  // H
  CW += c[2] * 1.5999000000000001e+01;  // O
  CW += c[3] * 3.1998000000000001e+01;  // O2
  CW += c[4] * 1.7007000000000001e+01;  // OH
  CW += c[5] * 1.8015000000000001e+01;  // H2O
  CW += c[6] * 3.3006000000000000e+01;  // HO2
  CW += c[7] * 1.5035000000000000e+01;  // CH3
  CW += c[8] * 1.6042999999999999e+01;  // CH4
  CW += c[9] * 2.8009999999999998e+01;  // CO
  CW += c[10] * 4.4009000000000000e+01; // CO2
  CW += c[11] * 3.0026000000000000e+01; // CH2O
  CW += c[12] * 2.8013999999999999e+01; // N2
  // Now compute conversion
  amrex::Real CWinv = 1.0 / CW;
  y[0] = c[0] * 2.0160000000000000e+00 * CWinv;
  y[1] = c[1] * 1.0080000000000000e+00 * CWinv;
  y[2] = c[2] * 1.5999000000000001e+01 * CWinv;
  y[3] = c[3] * 3.1998000000000001e+01 * CWinv;
  y[4] = c[4] * 1.7007000000000001e+01 * CWinv;
  y[5] = c[5] * 1.8015000000000001e+01 * CWinv;
  y[6] = c[6] * 3.3006000000000000e+01 * CWinv;
  y[7] = c[7] * 1.5035000000000000e+01 * CWinv;
  y[8] = c[8] * 1.6042999999999999e+01 * CWinv;
  y[9] = c[9] * 2.8009999999999998e+01 * CWinv;
  y[10] = c[10] * 4.4009000000000000e+01 * CWinv;
  y[11] = c[11] * 3.0026000000000000e+01 * CWinv;
  y[12] = c[12] * 2.8013999999999999e+01 * CWinv;
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

  amrex::Real CH2OH_CH2_S = (epsilon + 1.0 * qr_co[16]) / CH2OH_denom;

  // QSSA species 1: CH2(S)

  amrex::Real CH2_S_num = epsilon + 1.0 * qr_co[6] + 1.0 * qr_co[18] +
                          1.0 * qf_co[20] + 1.0 * qr_co[29] + 1.0 * qr_co[30] +
                          1.0 * qr_co[31] + 1.0 * qr_co[33] + 1.0 * qr_co[36];
  amrex::Real CH2_S_denom = epsilon - qf_co[6] - qr_co[16] - qf_co[18] -
                            qr_co[20] - qf_co[28] - qf_co[29] - qf_co[30] -
                            qf_co[31] - qf_co[32] - qf_co[33] - qf_co[34] -
                            qf_co[35] - qf_co[36];
  amrex::Real CH2_S_rhs = -CH2_S_num / CH2_S_denom;

  amrex::Real CH2_S_CH2OH = (epsilon + 1.0 * qf_co[16]) / CH2_S_denom;
  amrex::Real CH2_S_CH2 = (epsilon + 1.0 * qr_co[28] + 1.0 * qr_co[32] +
                           1.0 * qr_co[34] + 1.0 * qr_co[35]) /
                          CH2_S_denom;

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
  amrex::Real CH2_CH2_S = (epsilon + 1.0 * qf_co[28] + 1.0 * qf_co[32] +
                           1.0 * qf_co[34] + 1.0 * qf_co[35]) /
                          CH2_denom;

  // Putting it all together
  amrex::Real H_0 = CH2_HCO;
  amrex::Real H_1 = CH2_S_CH2OH;
  amrex::Real H_2 = CH2_CH2_S / (1 - CH2OH_CH2_S * H_1);
  sc_qss[0] =
    ((CH2_rhs - HCO_rhs * H_0) - (CH2_S_rhs - CH2OH_rhs * H_1) * H_2) /
    ((1 - HCO_CH2 * H_0) - CH2_S_CH2 * H_2);

  sc_qss[1] = ((CH2_S_rhs - CH2OH_rhs * H_1) - (CH2_S_CH2 * sc_qss[0])) /
              (1 - CH2OH_CH2_S * H_1);

  sc_qss[3] = CH2OH_rhs - (+CH2OH_CH2_S * sc_qss[1]);

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
  XW += x[0] * 2.0160000000000000e+00;  // H2
  XW += x[1] * 1.0080000000000000e+00;  // H
  XW += x[2] * 1.5999000000000001e+01;  // O
  XW += x[3] * 3.1998000000000001e+01;  // O2
  XW += x[4] * 1.7007000000000001e+01;  // OH
  XW += x[5] * 1.8015000000000001e+01;  // H2O
  XW += x[6] * 3.3006000000000000e+01;  // HO2
  XW += x[7] * 1.5035000000000000e+01;  // CH3
  XW += x[8] * 1.6042999999999999e+01;  // CH4
  XW += x[9] * 2.8009999999999998e+01;  // CO
  XW += x[10] * 4.4009000000000000e+01; // CO2
  XW += x[11] * 3.0026000000000000e+01; // CH2O
  XW += x[12] * 2.8013999999999999e+01; // N2
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
  const amrex::Real x3 =
    sc[11] + sc[12] + sc[1] + sc[2] + sc[4] + sc[6] + sc[7] + 2 * sc[8];
  const amrex::Real x7 = 2 * sc[0] + 2 * sc[10] + sc[3] + 1.5 * sc[9] + x3;
  const amrex::Real x8 = 6.0 * sc[5] + x7;
  amrex::Real x11 = -g_RT[11] - g_RT[1] + g_RT_qss[3];
  amrex::Real x12 = -3286.0048277048199 * invT - 4.82 * logT;
  amrex::Real x16 = 1.27e+20 * x8 * exp(x12);
  amrex::Real x18 = log(
    exp(-4160.0 * invT) + 0.2813 * exp(-0.0097087378640776708 * T) +
    0.7187 * exp(-0.00077459333849728897 * T));
  amrex::Real x21 = -0.29097730287517898 * x18 +
                    0.43429448190325198 * log(x16 / kf_qss[2]) - 0.4;
  amrex::Real x27 =
    ((x21) * (x21)) / ((-0.55155399201712996 * x18 - 0.14 * x21 + 0.75) *
                       (-0.55155399201712996 * x18 - 0.14 * x21 + 0.75)) +
    1.0;
  amrex::Real x30 =
    exp(M_LN10 * (0.43429448190325198 * x18 / x27)) / (1.0 + x16 / kf_qss[2]);
  amrex::Real x32 = refC * x30 * exp(x11 + x12);
  amrex::Real x34 = -kf_qss[10] * sc[2] - kf_qss[14] * sc[1] -
                    kf_qss[15] * sc[1] - kf_qss[16] * sc[1] -
                    kf_qss[23] * sc[4] - kf_qss[40] * sc[3] -
                    1.27e+20 * x32 * x8 + 10.0e-13;
  const amrex::Real x38 = exp(-g_RT[1] + g_RT[5] + g_RT_qss[1] - g_RT_qss[3]);
  const amrex::Real x42 = (kf_qss[16] * sc[5] * x38 + 10.0e-13) / x34;
  amrex::Real x45 = g_RT[4] - g_RT[5];
  const amrex::Real x47 = exp(-g_RT[7] + g_RT_qss[1] - x45);
  amrex::Real x51 =
    -kf_qss[16] * sc[5] * x38 - kf_qss[18] * sc[4] - kf_qss[20] * sc[5] * x47 -
    kf_qss[28] * sc[12] - kf_qss[29] * sc[3] - kf_qss[30] * sc[3] -
    kf_qss[31] * sc[0] - kf_qss[32] * sc[5] - kf_qss[33] * sc[8] -
    kf_qss[34] * sc[9] - kf_qss[35] * sc[10] - kf_qss[36] * sc[10] -
    kf_qss[6] * sc[2] + 10.0e-13;
  amrex::Real x53 = kf_qss[16] * sc[1] + 10.0e-13;
  amrex::Real x54 = x53 / x51;
  amrex::Real x55 = -x42 * x54 + 1;
  amrex::Real x58 = exp(g_RT_qss[0] - g_RT_qss[1]);
  amrex::Real x67 = kf_qss[28] * sc[12] * x58 + kf_qss[32] * sc[5] * x58 +
                    kf_qss[34] * sc[9] * x58 + kf_qss[35] * sc[10] * x58 +
                    10.0e-13;
  const amrex::Real x69 = x67 / (x51 * x55);
  amrex::Real x70 = g_RT[7] - g_RT[8];
  amrex::Real x72 = g_RT[11] - g_RT_qss[2];
  const amrex::Real x79 = exp(g_RT[1] - g_RT[2] - g_RT_qss[0] + g_RT_qss[2]);
  amrex::Real x84 = g_RT[2] - g_RT[4];
  amrex::Real x87 = g_RT[3] - g_RT[6];
  amrex::Real x91 = -g_RT[0] + g_RT[11] + g_RT[1] - g_RT_qss[2];
  amrex::Real x93 = -717.08374877172605 * invT - 2.57 * logT;
  amrex::Real x95 = 1350000000000.0 * exp(x93);
  amrex::Real x98 = log(
    exp(-6570.0 * invT) + 0.2176 * exp(-0.0036900369003690001 * T) +
    0.7824 * exp(-0.00036297640653357502 * T));
  amrex::Real x101 = -0.29097730287517898 * x98 +
                     0.43429448190325198 * log(x8 * x95 / kf_qss[1]) - 0.4;
  amrex::Real x106 =
    ((x101) * (x101)) / ((-0.14 * x101 - 0.55155399201712996 * x98 + 0.75) *
                         (-0.14 * x101 - 0.55155399201712996 * x98 + 0.75)) +
    1.0;
  amrex::Real x109 = exp(M_LN10 * (0.43429448190325198 * x98 / x106)) /
                     (1.0 + x8 * x95 / kf_qss[1]);
  amrex::Real x117 =
    kf_qss[25] * exp(-g_RT[3] + g_RT[4] - g_RT_qss[0] + g_RT_qss[2]);
  amrex::Real x119 = -kf_qss[11] * sc[6] * exp(-x72 - x87) -
                     kf_qss[12] * sc[1] - kf_qss[13] * sc[0] * exp(-x91) -
                     kf_qss[21] * sc[4] - kf_qss[22] * sc[5] * exp(-x45 - x72) -
                     kf_qss[37] * sc[7] - kf_qss[38] * sc[8] * exp(-x70 - x72) -
                     kf_qss[39] * sc[3] - kf_qss[3] * sc[5] - kf_qss[4] * x7 -
                     kf_qss[5] * sc[1] * x79 - kf_qss[7] * sc[2] -
                     kf_qss[8] * sc[2] - kf_qss[9] * sc[4] * exp(-x72 - x84) -
                     sc[1] * x109 * x8 * x95 - sc[4] * x117 + 10.0e-13;
  const amrex::Real x124 =
    (kf_qss[25] * sc[3] + kf_qss[5] * sc[2] + 10.0e-13) / x119;
  amrex::Real x125 = 0.43429448190325198 / x8;
  amrex::Real x128 =
    ((x101) * (x101)) / ((-0.14 * x101 - 0.55155399201712996 * x98 + 0.75) *
                         (-0.14 * x101 - 0.55155399201712996 * x98 + 0.75) *
                         (-0.14 * x101 - 0.55155399201712996 * x98 + 0.75));
  x101 = x101 / ((-0.14 * x101 - 0.55155399201712996 * x98 + 0.75) *
                 (-0.14 * x101 - 0.55155399201712996 * x98 + 0.75));
  amrex::Real x133 = sc[1] * x109 * x8 * x95 * x98 / ((x106) * (x106));
  amrex::Real x136 =
    exp(M_LN10 * (0.43429448190325198 * x98 / x106)) * x8 /
    (kf_qss[1] * ((1.0 + x8 * x95 / kf_qss[1]) * (1.0 + x8 * x95 / kf_qss[1])));
  amrex::Real x139 =
    sc[1] * x136 * exp(-1434.1674975434501 * invT - 5.14 * logT);
  amrex::Real x141 =
    -2 * kf_qss[4] - 2700000000000.0 * sc[1] * x109 * exp(x93) +
    x133 * (4.0 * x101 * x125 + 0.56 * x125 * x128) + 3.645e+24 * x139;
  amrex::Real x142 = -kf_qss[13] * exp(-x91) + x141;
  amrex::Real x148 = -618.95649893980499 * invT - 3.14 * logT;
  amrex::Real x150 = 3200000000000000.0 * exp(x148);
  amrex::Real x153 = log(
    exp(-5590.0 * invT) + 0.32 * exp(-0.012820512820512799 * T) +
    0.68 * exp(-0.00050125313283207998 * T));
  amrex::Real x156 = -0.29097730287517898 * x153 +
                     0.43429448190325198 * log(x150 * x8 / kf_qss[0]) - 0.4;
  amrex::Real x161 =
    ((x156) * (x156)) / ((-0.55155399201712996 * x153 - 0.14 * x156 + 0.75) *
                         (-0.55155399201712996 * x153 - 0.14 * x156 + 0.75)) +
    1.0;
  amrex::Real x164 = exp(M_LN10 * (0.43429448190325198 * x153 / x161)) /
                     (1.0 + x150 * x8 / kf_qss[0]);
  amrex::Real x168 = -kf_qss[17] * sc[4] -
                     kf_qss[19] * sc[5] * exp(-g_RT[7] + g_RT_qss[0] - x45) -
                     kf_qss[24] * sc[6] - kf_qss[25] * sc[3] -
                     kf_qss[26] * sc[0] - kf_qss[27] * sc[8] -
                     kf_qss[28] * sc[12] * x58 - kf_qss[32] * sc[5] * x58 -
                     kf_qss[34] * sc[9] * x58 - kf_qss[35] * sc[10] * x58 -
                     kf_qss[5] * sc[2] - sc[1] * x150 * x164 * x8 + 10.0e-13;
  amrex::Real x170 = kf_qss[5] * sc[1] * x79 + sc[4] * x117 + 10.0e-13;
  amrex::Real x171 = x170 / x168;
  amrex::Real x172 = x171 *
                     (kf_qss[25] * sc[3] + kf_qss[5] * sc[2] + 10.0e-13) /
                     ((x119) * (x119));
  amrex::Real x174 = x125 * ((x156) * (x156)) /
                     ((-0.55155399201712996 * x153 - 0.14 * x156 + 0.75) *
                      (-0.55155399201712996 * x153 - 0.14 * x156 + 0.75) *
                      (-0.55155399201712996 * x153 - 0.14 * x156 + 0.75));
  x156 = x125 * x156 /
         ((-0.55155399201712996 * x153 - 0.14 * x156 + 0.75) *
          (-0.55155399201712996 * x153 - 0.14 * x156 + 0.75));
  amrex::Real x178 = sc[1] * x150 * x153 * x164 * x8 / ((x161) * (x161));
  amrex::Real x179 = exp(M_LN10 * (0.43429448190325198 * x153 / x161)) * x8 /
                     (kf_qss[0] * ((1.0 + x150 * x8 / kf_qss[0]) *
                                   (1.0 + x150 * x8 / kf_qss[0])));
  amrex::Real x186 =
    -6400000000000000.0 * sc[1] * x164 * exp(x148) +
    2.048e+31 * sc[1] * x179 * exp(-1237.91299787961 * invT - 6.28 * logT) +
    x178 * (4.0 * x156 + 0.56 * x174);
  amrex::Real x190 = x124 * x170 / ((x168) * (x168));
  amrex::Real x192 = x125 * ((x21) * (x21)) /
                     ((-0.55155399201712996 * x18 - 0.14 * x21 + 0.75) *
                      (-0.55155399201712996 * x18 - 0.14 * x21 + 0.75) *
                      (-0.55155399201712996 * x18 - 0.14 * x21 + 0.75));
  x21 = x125 * x21 /
        ((-0.55155399201712996 * x18 - 0.14 * x21 + 0.75) *
         (-0.55155399201712996 * x18 - 0.14 * x21 + 0.75));
  amrex::Real x196 = 1.27e+20 * x18 * x32 * x8 / ((x27) * (x27));
  amrex::Real x197 =
    exp(M_LN10 * (0.43429448190325198 * x18 / x27)) /
    (kf_qss[2] * ((1.0 + x16 / kf_qss[2]) * (1.0 + x16 / kf_qss[2])));
  amrex::Real x200 =
    refC * x8 * exp(-6572.0096554096399 * invT - 9.64 * logT + x11);
  amrex::Real x202 = x196 * (0.56 * x192 + 4.0 * x21) +
                     3.2258e+40 * x197 * x200 - 2.54e+20 * x32;
  amrex::Real x204 = (kf_qss[16] * sc[5] * x38 + 10.0e-13) / ((x34) * (x34));
  amrex::Real x205 = x204 * x54;
  amrex::Real x209 = x42 * x53 / ((x51) * (x51));
  amrex::Real x212 = kf_qss[28] * sc[12] + kf_qss[32] * sc[5] +
                     kf_qss[34] * sc[9] + kf_qss[35] * sc[10] + 10.0e-13;
  amrex::Real x215 = x212 * x67 / (x168 * x51 * ((x55) * (x55)));
  amrex::Real x217 = x212 * x69 / ((x168) * (x168));
  amrex::Real x218 = x212 / (x168 * x55);
  amrex::Real x220 = x218 * x67 / ((x51) * (x51));
  amrex::Real x224 = kf_qss[40] * exp(g_RT[11] - g_RT_qss[3] - x87);
  const amrex::Real x226 = exp(g_RT[11] - g_RT_qss[3] - x45);
  const amrex::Real x231 = exp(-g_RT[1] + g_RT[4] + g_RT[7] - g_RT_qss[3]);
  const amrex::Real x235 = exp(g_RT[0] + g_RT[11] - g_RT[1] - g_RT_qss[3]);
  const amrex::Real x238 = exp(g_RT[11] - g_RT_qss[3] - x84);
  amrex::Real x244 =
    kf_qss[10] * sc[11] * sc[4] * x238 + kf_qss[14] * sc[0] * sc[11] * x235 +
    kf_qss[15] * sc[4] * sc[7] * x231 + kf_qss[23] * sc[11] * sc[5] * x226 +
    sc[11] * sc[1] * x16 * x30 + sc[11] * sc[6] * x224 + 10.0e-13;
  const amrex::Real x245 = x244 / x34;
  const amrex::Real x248 = exp(-g_RT[0] + g_RT[1] + g_RT[7] - g_RT_qss[1]);
  const amrex::Real x254 = exp(-g_RT[3] + g_RT[5] + g_RT[9] - g_RT_qss[1]);
  const amrex::Real x257 =
    exp(g_RT[1] - g_RT[3] + g_RT[4] + g_RT[9] - g_RT_qss[1]);
  const amrex::Real x258 = refCinv * sc[1];
  const amrex::Real x263 = exp(-g_RT[10] + g_RT[11] + g_RT[9] - g_RT_qss[1]);
  const amrex::Real x267 = exp(g_RT[11] + g_RT[1] - g_RT[4] - g_RT_qss[1]);
  const amrex::Real x271 = -2.0 * g_RT[7] + g_RT[8];
  const amrex::Real x274 = exp(g_RT[0] - g_RT[2] + g_RT[9] - g_RT_qss[1]);
  amrex::Real x277 =
    kf_qss[18] * sc[11] * sc[1] * x267 + kf_qss[20] * sc[4] * sc[7] +
    kf_qss[29] * sc[4] * sc[9] * x257 * x258 +
    kf_qss[30] * sc[5] * sc[9] * x254 + kf_qss[31] * sc[1] * sc[7] * x248 +
    kf_qss[33] * ((sc[7]) * (sc[7])) * exp(-g_RT_qss[1] - x271) +
    kf_qss[36] * sc[11] * sc[9] * x263 + kf_qss[6] * sc[0] * sc[9] * x274 +
    10.0e-13;
  amrex::Real x278 = x245 * x54 - x277 / x51;
  amrex::Real x279 = -g_RT[9] + g_RT_qss[2];
  amrex::Real x281 = kf_qss[39] * exp(-x279 - x87);
  const amrex::Real x283 = refC * sc[11];
  amrex::Real x285 = g_RT[11] - g_RT[1] - g_RT_qss[2];
  amrex::Real x286 = x109 * exp(x285 + x93);
  amrex::Real x290 = exp(g_RT[1] - x279);
  const amrex::Real x293 = exp(-x279 - x84);
  const amrex::Real x298 = exp(g_RT[10] + g_RT[1] - g_RT[2] - g_RT_qss[2]);
  const amrex::Real x302 = exp(-x279 - x45);
  const amrex::Real x307 = exp(g_RT[0] - g_RT[1] - x279);
  const amrex::Real x311 = exp(-x279 - x70);
  amrex::Real x315 =
    kf_qss[11] * sc[11] * sc[3] + kf_qss[12] * sc[0] * sc[9] * x307 +
    kf_qss[13] * sc[11] * sc[1] + kf_qss[21] * sc[5] * sc[9] * x302 +
    kf_qss[22] * sc[11] * sc[4] + kf_qss[37] * sc[8] * sc[9] * x311 +
    kf_qss[38] * sc[11] * sc[7] + kf_qss[3] * sc[5] * sc[9] * x258 * x290 +
    kf_qss[4] * sc[9] * x258 * x290 * x7 + kf_qss[7] * sc[4] * sc[9] * x293 +
    kf_qss[8] * sc[10] * sc[1] * x298 + kf_qss[9] * sc[11] * sc[2] +
    sc[6] * sc[9] * x281 + 1350000000000.0 * x283 * x286 * x8 + 10.0e-13;
  const amrex::Real x316 = x315 / x119;
  amrex::Real x320 =
    kf_qss[26] * exp(-g_RT[0] + g_RT[1] + g_RT[7] - g_RT_qss[0]);
  const amrex::Real x324 = exp(g_RT[11] + g_RT[4] - g_RT[6] - g_RT_qss[0]);
  const amrex::Real x328 = exp(g_RT[11] + g_RT[1] - g_RT[4] - g_RT_qss[0]);
  amrex::Real x331 = -g_RT[1] + g_RT[7] - g_RT_qss[0];
  amrex::Real x335 = refC * sc[7] * x164 * exp(x148 + x331);
  amrex::Real x338 =
    kf_qss[17] * sc[11] * sc[1] * x328 + kf_qss[19] * sc[4] * sc[7] +
    kf_qss[24] * sc[11] * sc[4] * x324 +
    kf_qss[27] * ((sc[7]) * (sc[7])) * exp(-g_RT_qss[0] - x271) +
    sc[1] * sc[7] * x320 + 3200000000000000.0 * x335 * x8 + 10.0e-13;
  amrex::Real x341 = (x171 * x316 - x218 * x278 - x338 / x168) /
                     ((-x124 * x171 + 1 - x212 * x69 / x168) *
                      (-x124 * x171 + 1 - x212 * x69 / x168));
  const amrex::Real x342 =
    x341 *
    (-kf_qss[31] * x220 + x142 * x172 + x190 * (-kf_qss[26] + x186) +
     x215 * (-kf_qss[31] * x209 + x202 * x205) + x217 * (-kf_qss[26] + x186));
  amrex::Real x343 = x338 / ((x168) * (x168));
  x161 = 3200000000000000.0 * x153 * x335 * x8 / ((x161) * (x161));
  amrex::Real x345 =
    refC * sc[7] * exp(-1237.91299787961 * invT - 6.28 * logT + x331);
  amrex::Real x346 =
    -(6400000000000000.0 * refC * sc[7] * x164 * exp(x148 + x331) -
      x161 * (4.0 * x156 + 0.56 * x174) - 2.048e+31 * x179 * x345) /
    x168;
  x315 = x315 / ((x119) * (x119));
  amrex::Real x349 = x170 * x316 / ((x168) * (x168));
  x106 = 1350000000000.0 * x283 * x286 * x8 * x98 / ((x106) * (x106));
  amrex::Real x352 =
    x283 * exp(-1434.1674975434501 * invT - 5.14 * logT + x285);
  amrex::Real x354 = 2 * kf_qss[4] * sc[9] * x258 * x290 -
                     x106 * (4.0 * x101 * x125 + 0.56 * x125 * x128) -
                     3.645e+24 * x136 * x352 + 2700000000000.0 * x283 * x286;
  amrex::Real x356 = x212 * x278 / (x168 * ((x55) * (x55)));
  amrex::Real x359 = x212 * x278 / (((x168) * (x168)) * x55);
  x277 = x277 / ((x51) * (x51));
  x244 = x244 / ((x34) * (x34));
  x27 = sc[11] * sc[1] * x16 * x18 * x30 / ((x27) * (x27));
  amrex::Real x367 = sc[11] * sc[1] * x30 * exp(x12);
  amrex::Real x370 =
    sc[11] * sc[1] * x8 * exp(-6572.0096554096399 * invT - 9.64 * logT);
  amrex::Real x372 = -3.2258e+40 * x197 * x370 -
                     x27 * (0.56 * x192 + 4.0 * x21) + 2.54e+20 * x367;
  amrex::Real x374 = x245 * x53 / ((x51) * (x51));
  amrex::Real x375 = -kf_qss[31] * x277 + kf_qss[31] * x374 -
                     kf_qss[6] * sc[9] * x274 / x51 - x202 * x244 * x54 +
                     x54 * (kf_qss[14] * sc[11] * x235 + x372) / x34;
  amrex::Real x376 = (1.0 / ((-x124 * x171 + 1 - x212 * x69 / x168)));
  const amrex::Real x377 =
    x376 *
    (-x142 * x171 * x315 - x218 * x375 + x343 * (-kf_qss[26] + x186) + x346 -
     x349 * (-kf_qss[26] + x186) + x356 * (-kf_qss[31] * x209 + x202 * x205) +
     x359 * (-kf_qss[26] + x186) +
     x171 * (kf_qss[12] * sc[9] * x307 + x354) / x119);
  amrex::Real x382 = -kf_qss[4] - sc[1] * x109 * x95 +
                     x133 * (2.0 * x101 * x125 + 0.28 * x125 * x128) +
                     1.8225e+24 * x139;
  x95 = -kf_qss[12] - kf_qss[5] * x79 - x109 * x8 * x95 + x382;
  amrex::Real x385 = sc[1] * x179 * exp(-1237.91299787961 * invT - 6.28 * logT);
  amrex::Real x386 =
    -sc[1] * x150 * x164 + x178 * (2.0 * x156 + 0.28 * x174) + 1.024e+31 * x385;
  x150 = -x150 * x164 * x8 + x386;
  amrex::Real x389 = x196 * (0.28 * x192 + 2.0 * x21) +
                     1.6129e+40 * x197 * x200 - 1.27e+20 * x32;
  amrex::Real x390 = -kf_qss[14] - kf_qss[15] - kf_qss[16] + x389;
  amrex::Real x392 = -kf_qss[16] * x42 / x51 + x205 * x390;
  const amrex::Real x393 =
    x341 * (-kf_qss[5] * x124 * x79 / x168 + x150 * x190 + x150 * x217 +
            x172 * x95 + x215 * x392);
  amrex::Real x395 = -x161 * (2.0 * x156 + 0.28 * x174) -
                     1.024e+31 * x179 * x345 + 3200000000000000.0 * x335;
  amrex::Real x399 = kf_qss[4] * sc[9] * x258 * x290 -
                     x106 * (2.0 * x101 * x125 + 0.28 * x125 * x128) -
                     1.8225e+24 * x136 * x352 + 1350000000000.0 * x283 * x286;
  amrex::Real x400 =
    (kf_qss[13] * sc[11] + kf_qss[3] * refCinv * sc[5] * sc[9] * x290 +
     kf_qss[4] * refCinv * sc[9] * x290 * x7 + kf_qss[8] * sc[10] * x298 +
     x399) /
    x119;
  amrex::Real x405 = -1.6129e+40 * x197 * x370 -
                     x27 * (0.28 * x192 + 2.0 * x21) + 1.27e+20 * x367;
  const amrex::Real x407 =
    (kf_qss[16] * x245 / x51 - x244 * x390 * x54 -
     (kf_qss[18] * sc[11] * x267 + kf_qss[29] * refCinv * sc[4] * sc[9] * x257 +
      kf_qss[31] * sc[7] * x248) /
       x51 +
     x54 * (sc[11] * x16 * x30 + x405) / x34) /
    x55;
  x150 = x376 * (kf_qss[5] * x316 * x79 / x168 + x150 * x343 - x150 * x349 +
                 x150 * x359 - x171 * x315 * x95 + x171 * x400 + x356 * x392 -
                 x212 * x407 / x168 -
                 (kf_qss[17] * sc[11] * x328 + sc[7] * x320 + x395) / x168);
  amrex::Real x409 = -kf_qss[7] - kf_qss[8] + x382;
  amrex::Real x415 =
    -kf_qss[6] * x42 * x53 / ((x51) * (x51)) + x205 * (-kf_qss[10] + x389);
  const amrex::Real x416 = x341 * (-kf_qss[5] * x171 / x119 -
                                   kf_qss[6] * x218 * x67 / ((x51) * (x51)) +
                                   x172 * x409 + x190 * (-kf_qss[5] + x386) +
                                   x215 * x415 + x217 * (-kf_qss[5] + x386));
  amrex::Real x423 = kf_qss[6] * x245 * x53 / ((x51) * (x51)) -
                     kf_qss[6] * x277 - x244 * x54 * (-kf_qss[10] + x389) +
                     x405 * x54 / x34;
  const amrex::Real x424 =
    x376 *
    (-x171 * x315 * x409 - x218 * x423 + x343 * (-kf_qss[5] + x386) -
     x349 * (-kf_qss[5] + x386) + x356 * x415 + x359 * (-kf_qss[5] + x386) -
     x395 / x168 + x171 * (kf_qss[9] * sc[11] + x399) / x119);
  amrex::Real x427 = -kf_qss[29] - kf_qss[30];
  amrex::Real x428 = x205 * (-kf_qss[40] + x389) + x209 * x427;
  const amrex::Real x430 =
    x341 * (-kf_qss[25] * x171 / x119 + x172 * (-kf_qss[39] + x382) +
            x190 * (-kf_qss[25] + x386) + x215 * x428 +
            x217 * (-kf_qss[25] + x386) + x220 * x427);
  amrex::Real x434 = -x244 * x54 * (-kf_qss[40] + x389) + x277 * x427 -
                     x374 * x427 + x405 * x54 / x34;
  const amrex::Real x435 =
    x376 * (-x171 * x315 * (-kf_qss[39] + x382) - x218 * x434 +
            x343 * (-kf_qss[25] + x386) - x349 * (-kf_qss[25] + x386) +
            x356 * x428 + x359 * (-kf_qss[25] + x386) - x395 / x168 +
            x171 * (kf_qss[11] * sc[11] + x399) / x119);
  amrex::Real x437 = -kf_qss[21] - kf_qss[9] * exp(-x72 - x84) - x117 + x382;
  amrex::Real x439 = (-kf_qss[17] + x386) / ((x168) * (x168));
  amrex::Real x442 = -kf_qss[18] * x209 + x205 * (-kf_qss[23] + x389);
  const amrex::Real x443 =
    x341 * (-kf_qss[18] * x220 - x117 * x124 / x168 + x124 * x170 * x439 +
            x172 * x437 + x212 * x439 * x69 + x215 * x442);
  amrex::Real x449 =
    -kf_qss[18] * x277 + kf_qss[18] * x374 - x244 * x54 * (-kf_qss[23] + x389) -
    (kf_qss[20] * sc[7] + kf_qss[29] * sc[9] * x257 * x258) / x51 +
    x54 * (kf_qss[10] * sc[11] * x238 + kf_qss[15] * sc[7] * x231 + x405) / x34;
  x117 =
    x376 *
    (x117 * x316 / x168 - x170 * x316 * x439 - x171 * x315 * x437 +
     x212 * x278 * x439 / x55 - x218 * x449 + x343 * (-kf_qss[17] + x386) +
     x356 * x442 -
     (kf_qss[19] * sc[7] + kf_qss[24] * sc[11] * x324 + x395) / x168 +
     x171 * (kf_qss[22] * sc[11] + kf_qss[7] * sc[9] * x293 + x399) / x119);
  x170 = -kf_qss[19] * exp(-g_RT[7] + g_RT_qss[0] - x45) - kf_qss[32] * x58 -
         1.92e+16 * sc[1] * x164 * exp(x148) +
         x178 * (12.0 * x156 + 1.68 * x174) + 6.144e+31 * x385;
  x32 = x196 * (1.68 * x192 + 12.0 * x21) + 9.6774e+40 * x197 * x200 -
        7.62e+20 * x32;
  x439 = -kf_qss[16] * x38 - kf_qss[20] * x47 - kf_qss[32];
  amrex::Real x458 = -kf_qss[16] * x38 * x54 / x34 + x205 * x32 + x209 * x439;
  amrex::Real x463 =
    -kf_qss[22] * exp(-x45 - x72) - kf_qss[3] -
    8100000000000.0 * sc[1] * x109 * exp(x93) +
    1.0935e+25 * sc[1] * x136 * exp(-1434.1674975434501 * invT - 5.14 * logT) +
    x133 * (12.0 * x101 * x125 + 1.68 * x125 * x128);
  const amrex::Real x466 =
    x341 *
    (-kf_qss[32] * x218 * x58 / x51 - kf_qss[32] * x69 / x168 + x170 * x190 +
     x170 * x217 + x172 * x463 + x215 * x458 + x220 * x439);
  amrex::Real x468 =
    (kf_qss[21] * sc[9] * x302 + kf_qss[3] * sc[9] * x258 * x290 -
     x106 * (12.0 * x101 * x125 + 1.68 * x125 * x128) -
     1.0935e+25 * x136 * x352 + 8100000000000.0 * x283 * x286) /
    x119;
  amrex::Real x472 = (kf_qss[23] * sc[11] * x226 - 9.6774e+40 * x197 * x370 -
                      x27 * (1.68 * x192 + 12.0 * x21) + 7.62e+20 * x367) /
                     x34;
  amrex::Real x473 = -kf_qss[30] * sc[9] * x254 / x51 - x244 * x32 * x54 +
                     x277 * x439 - x374 * x439 + x472 * x54;
  x170 = x376 * (-kf_qss[32] * x278 / (x168 * x55) + x170 * x343 - x170 * x349 +
                 x170 * x359 - x171 * x315 * x463 + x171 * x468 - x218 * x473 +
                 x356 * x458 -
                 (-x161 * (12.0 * x156 + 1.68 * x174) -
                  6.144e+31 * x179 * x345 + 1.92e+16 * x335) /
                   x168);
  amrex::Real x476 = -kf_qss[11] * exp(-x72 - x87) + x382;
  x53 =
    x204 * x212 * x389 * x53 * x67 / (x168 * ((x51) * (x51)) * ((x55) * (x55)));
  const amrex::Real x478 = x341 * (x172 * x476 + x190 * (-kf_qss[24] + x386) +
                                   x217 * (-kf_qss[24] + x386) + x53);
  amrex::Real x484 = -x244 * x389 * x54;
  const amrex::Real x486 = (x484 + x54 * (sc[11] * x224 + x405) / x34) / x55;
  const amrex::Real x487 =
    x376 *
    (-x171 * x315 * x476 + x205 * x356 * x389 + x343 * (-kf_qss[24] + x386) -
     x349 * (-kf_qss[24] + x386) + x359 * (-kf_qss[24] + x386) -
     x212 * x486 / x168 - x395 / x168 + x171 * (sc[9] * x281 + x399) / x119);
  const amrex::Real x490 =
    x341 * (x172 * (-kf_qss[37] + x382) + x190 * x386 + x217 * x386 + x53);
  const amrex::Real x495 =
    (x484 -
     (kf_qss[20] * sc[4] + kf_qss[31] * sc[1] * x248 +
      2.0 * kf_qss[33] * ((sc[7])) * exp(-g_RT_qss[1] - x271)) /
       x51 +
     x54 * (kf_qss[15] * sc[4] * x231 + x405) / x34) /
    x55;
  x320 = x376 * (-x171 * x315 * (-kf_qss[37] + x382) + x205 * x356 * x389 +
                 x343 * x386 - x349 * x386 + x359 * x386 - x212 * x495 / x168 -
                 (kf_qss[19] * sc[4] +
                  2.0 * kf_qss[27] * ((sc[7])) * exp(-g_RT_qss[0] - x271) +
                  3200000000000000.0 * refC * x164 * x8 * exp(x148 + x331) +
                  sc[1] * x320 + x395) /
                   x168 +
                 x171 * (kf_qss[38] * sc[11] + x399) / x119);
  amrex::Real x498 = -kf_qss[38] * exp(-x70 - x72) + x141;
  const amrex::Real x501 =
    x341 *
    (-kf_qss[33] * x220 + x172 * x498 + x190 * (-kf_qss[27] + x186) +
     x215 * (-kf_qss[33] * x209 + x202 * x205) + x217 * (-kf_qss[27] + x186));
  amrex::Real x506 = -x202 * x244 * x54 + x372 * x54 / x34;
  const amrex::Real x508 =
    x376 *
    (-x171 * x315 * x498 -
     x218 * (-kf_qss[33] * x277 + kf_qss[33] * x374 + x506) +
     x343 * (-kf_qss[27] + x186) + x346 - x349 * (-kf_qss[27] + x186) +
     x356 * (-kf_qss[33] * x209 + x202 * x205) + x359 * (-kf_qss[27] + x186) +
     x171 * (kf_qss[37] * sc[9] * x311 + x354) / x119);
  x109 = -1.5 * kf_qss[4] - 2025000000000.0 * sc[1] * x109 * exp(x93) +
         x133 * (3.0 * x101 * x125 + 0.42 * x125 * x128) + 2.73375e+24 * x139;
  x133 = -kf_qss[34] * x58 - 4800000000000000.0 * sc[1] * x164 * exp(x148) +
         x178 * (3.0 * x156 + 0.42 * x174) + 1.536e+31 * x385;
  x139 = -1.905e+20 * refC * x30 * exp(x11 + x12) +
         x196 * (0.42 * x192 + 3.0 * x21) + 2.41935e+40 * x197 * x200;
  x148 = x341 * (-kf_qss[34] * x218 * x58 / x51 - kf_qss[34] * x220 -
                 kf_qss[34] * x69 / x168 + x109 * x172 + x133 * x190 +
                 x133 * x217 + x215 * (-kf_qss[34] * x209 + x139 * x205));
  x128 = (kf_qss[12] * sc[0] * x307 + kf_qss[21] * sc[5] * x302 +
          kf_qss[37] * sc[8] * x311 + kf_qss[3] * sc[5] * x258 * x290 +
          1.5 * kf_qss[4] * sc[9] * x258 * x290 + kf_qss[4] * x258 * x290 * x7 +
          kf_qss[7] * sc[4] * x293 + sc[6] * x281 -
          x106 * (3.0 * x101 * x125 + 0.42 * x125 * x128) -
          2.73375e+24 * x136 * x352 + 2025000000000.0 * x283 * x286) /
         x119;
  x101 = (-2.41935e+40 * x197 * x370 - x27 * (0.42 * x192 + 3.0 * x21) +
          1.905e+20 * x367) /
         x34;
  x136 = -kf_qss[34] * x277 + kf_qss[34] * x374 + x101 * x54 -
         x139 * x244 * x54 -
         (kf_qss[29] * sc[4] * x257 * x258 + kf_qss[30] * sc[5] * x254 +
          kf_qss[36] * sc[11] * x263 + kf_qss[6] * sc[0] * x274) /
           x51;
  x164 = x376 * (-kf_qss[34] * x278 / (x168 * x55) - x109 * x171 * x315 +
                 x128 * x171 + x133 * x343 - x133 * x349 + x133 * x359 -
                 x136 * x218 + x356 * (-kf_qss[34] * x209 + x139 * x205) -
                 (-x161 * (3.0 * x156 + 0.42 * x174) - 1.536e+31 * x179 * x345 +
                  4800000000000000.0 * x335) /
                   x168);
  x174 = -kf_qss[35] * x58 + x186;
  x156 = -kf_qss[35] - kf_qss[36];
  x178 = x341 * (-kf_qss[35] * x218 * x58 / x51 - kf_qss[35] * x69 / x168 +
                 x141 * x172 + x156 * x220 + x174 * x190 + x174 * x217 +
                 x215 * (x156 * x209 + x202 * x205));
  x179 = x376 * (-kf_qss[35] * x278 / (x168 * x55) - x141 * x171 * x315 +
                 x174 * x343 - x174 * x349 + x174 * x359 -
                 x218 * (x156 * x277 - x156 * x374 + x506) + x346 +
                 x356 * (x156 * x209 + x202 * x205) +
                 x171 * (kf_qss[8] * sc[1] * x298 + x354) / x119);
  x186 = x341 * (x172 * x382 + x190 * x386 + x217 * x386 + x53);
  x192 = (kf_qss[11] * sc[3] + kf_qss[13] * sc[1] + kf_qss[22] * sc[4] +
          kf_qss[38] * sc[7] + kf_qss[9] * sc[2] +
          1350000000000.0 * refC * x286 * x8 + x399) /
         x119;
  x21 = (kf_qss[10] * sc[4] * x238 + kf_qss[14] * sc[0] * x235 +
         kf_qss[23] * sc[5] * x226 + sc[1] * x16 * x30 + sc[6] * x224 + x405) /
        x34;
  x196 = (x21 * x54 + x484 -
          (kf_qss[18] * sc[1] * x267 + kf_qss[36] * sc[9] * x263) / x51) /
         x55;
  x197 =
    x376 *
    (x171 * x192 - x171 * x315 * x382 + x205 * x356 * x389 + x343 * x386 -
     x349 * x386 + x359 * x386 - x196 * x212 / x168 -
     (kf_qss[17] * sc[1] * x328 + kf_qss[24] * sc[4] * x324 + x395) / x168);
  x200 = -kf_qss[28] * x58 + x386;
  x172 = x341 * (-kf_qss[28] * x218 * x58 / x51 - kf_qss[28] * x220 -
                 kf_qss[28] * x69 / x168 + x172 * x382 + x190 * x200 +
                 x200 * x217 + x215 * (-kf_qss[28] * x209 + x205 * x389));
  x190 =
    x376 *
    (-kf_qss[28] * x278 / (x168 * x55) - x171 * x315 * x382 + x200 * x343 -
     x200 * x349 + x200 * x359 -
     x218 * (-kf_qss[28] * x277 + kf_qss[28] * x374 + x484 + x405 * x54 / x34) +
     x356 * (-kf_qss[28] * x209 + x205 * x389) - x395 / x168 +
     x171 * x399 / x119);
  x212 = sc_qss[0] * x67 / ((x51) * (x51));
  x215 = (-kf_qss[31] * x212 + x375) / x55;
  x217 = (-sc_qss[0] * x67 / x51 + x278) / ((x55) * (x55));
  x220 = x217 * (-kf_qss[31] * x209 + x202 * x205);
  x286 = x217 * x392;
  x290 = (-kf_qss[6] * sc_qss[0] * x67 / ((x51) * (x51)) + x423) / x55;
  x335 = x217 * x415;
  x341 = (x212 * x427 + x434) / x55;
  x343 = x217 * x428;
  x161 = (-kf_qss[18] * x212 + x449) / x55;
  x345 = x217 * x442;
  x346 = (-kf_qss[32] * sc_qss[0] * x58 / x51 + x212 * x439 + x473) / x55;
  x349 = x217 * x458;
  x106 = x205 * x217 * x389;
  x352 =
    (-kf_qss[33] * x212 - kf_qss[33] * x277 + kf_qss[33] * x374 + x506) / x55;
  x356 = x217 * (-kf_qss[33] * x209 + x202 * x205);
  x359 = (-kf_qss[34] * sc_qss[0] * x58 / x51 - kf_qss[34] * x212 + x136) / x55;
  x27 = x217 * (-kf_qss[34] * x209 + x139 * x205);
  x367 = (-kf_qss[35] * sc_qss[0] * x58 / x51 + x156 * x212 + x156 * x277 -
          x156 * x374 + x506) /
         x55;
  x370 = x217 * (x156 * x209 + x202 * x205);
  x54 = (-kf_qss[28] * sc_qss[0] * x58 / x51 - kf_qss[28] * x212 -
         kf_qss[28] * x277 + kf_qss[28] * x374 + x484 + x405 * x54 / x34) /
        x55;
  x58 = x217 * (-kf_qss[28] * x209 + x205 * x389);
  x205 = sc_qss[0] * (kf_qss[25] * sc[3] + kf_qss[5] * sc[2] + 10.0e-13) /
         ((x119) * (x119));
  x142 = x142 * x205 + x142 * x315 - (kf_qss[12] * sc[9] * x307 + x354) / x119;
  x209 = x205 * x95 + x315 * x95 - x400;
  x277 = -kf_qss[5] * sc_qss[0] / x119 + x205 * x409 + x315 * x409 -
         (kf_qss[9] * sc[11] + x399) / x119;
  x374 = -kf_qss[25] * sc_qss[0] / x119 + x205 * (-kf_qss[39] + x382) +
         x315 * (-kf_qss[39] + x382) - (kf_qss[11] * sc[11] + x399) / x119;
  x375 = x205 * x437 + x315 * x437 -
         (kf_qss[22] * sc[11] + kf_qss[7] * sc[9] * x293 + x399) / x119;
  x95 = x205 * x463 + x315 * x463 - x468;
  x281 = x205 * x476 + x315 * x476 - (sc[9] * x281 + x399) / x119;
  x385 = x205 * (-kf_qss[37] + x382) + x315 * (-kf_qss[37] + x382) -
         (kf_qss[38] * sc[11] + x399) / x119;
  x386 = x205 * x498 + x315 * x498 - (kf_qss[37] * sc[9] * x311 + x354) / x119;
  x392 = x109 * x205 + x109 * x315 - x128;
  x141 = x141 * x205 + x141 * x315 - (kf_qss[8] * sc[1] * x298 + x354) / x119;
  x354 = -x192 + x205 * x382 + x315 * x382;
  x119 = x205 * x382 + x315 * x382 - x399 / x119;
  x315 = sc_qss[1] * x202 * x204 + x202 * x244 -
         (kf_qss[14] * sc[11] * x235 + x372) / x34;
  x16 =
    sc_qss[1] * x204 * x390 + x244 * x390 - (sc[11] * x16 * x30 + x405) / x34;
  x30 = sc_qss[1] * x204 * (-kf_qss[10] + x389) + x244 * (-kf_qss[10] + x389) -
        x405 / x34;
  x382 = sc_qss[1] * x204 * (-kf_qss[40] + x389) + x244 * (-kf_qss[40] + x389) -
         x405 / x34;
  x390 = sc_qss[1] * x204 * (-kf_qss[23] + x389) + x244 * (-kf_qss[23] + x389) -
         (kf_qss[10] * sc[11] * x238 + kf_qss[15] * sc[7] * x231 + x405) / x34;
  x395 = -kf_qss[16] * sc_qss[1] * x38 / x34 + sc_qss[1] * x204 * x32 +
         x244 * x32 - x472;
  x389 = sc_qss[1] * x204 * x389 + x244 * x389;
  x224 = x389 - (sc[11] * x224 + x405) / x34;
  x399 = x389 - (kf_qss[15] * sc[4] * x231 + x405) / x34;
  x202 = sc_qss[1] * x202 * x204 + x202 * x244 - x372 / x34;
  x204 = sc_qss[1] * x139 * x204 - x101 + x139 * x244;
  x244 = -x21 + x389;
  x34 = x389 - x405 / x34;
  x91 = 23000.0 * exp(-1648.03458051046 * invT + 1.05 * logT - x91);
  x372 = sc[0] * x91;
  x400 = 80000000.0 * sc[1] * x79;
  x405 = -717.08374877172594 * invT - 2.57 * logT;
  x93 =
    1238532.11009174 * x8 * exp(-130.83633310922701 * invT - 0.48 * logT + x93);
  x409 = -0.29097730287517898 * x98 + 0.43429448190325198 * log(x93) - 0.4;
  x415 =
    ((x409) * (x409)) / ((-0.14 * x409 - 0.55155399201712996 * x98 + 0.75) *
                         (-0.14 * x409 - 0.55155399201712996 * x98 + 0.75)) +
    1.0;
  x423 = exp(M_LN10 * (0.43429448190325198 * x98 / x415)) / (x93 + 1.0);
  x427 = 1350000000000.0 * sc[1] * x423 * x8 * exp(x405);
  x428 = 6000000.0 * sc[5] * x38;
  x434 = exp(-8554.6833186802305 * invT - logT);
  x437 = 2244000000000.0 * sc[5] * x434;
  x442 = 187000000000.0 * x434 * x7;
  x12 =
    235185185185185.0 * x8 * exp(1811.57999689699 * invT - 0.454 * logT + x12);
  x449 = -0.29097730287517898 * x18 + 0.43429448190325198 * log(x12) - 0.4;
  x32 = -0.55155399201712996 * x18 - 0.14 * x449 + 0.75;
  x439 = exp(
           M_LN10 * (0.43429448190325198 * x18 /
                     (1.0 + ((x449) * (x449)) / ((x32) * (x32))))) /
         (x12 + 1.0);
  x458 = 1.27e+20 * x439 * x8;
  x463 = -3286.0048277048204 * invT - 4.82 * logT;
  x468 = refC * x458 * exp(x11 + x463);
  x472 = 39000000.0 * exp(-1781.38699694871 * invT - x72 - x84);
  x473 = sc[4] * x472;
  x476 = 100000000.0 * exp(-20128.666632188801 * invT - x72 - x87);
  x53 = sc[6] * x476;
  x484 = 7600000.0 * sc[3] * exp(-201.28666632188799 * invT);
  x498 = 13200000.0 * exp(
                        -g_RT[3] + g_RT[4] - g_RT_qss[0] + g_RT_qss[2] -
                        754.82499870708 * invT);
  x506 = sc[4] * x498;
  x109 = 18000000.0 * sc[3] * exp(-452.89499922424801 * invT);
  x133 = 25010000.0 * sc[5] * x47;
  x139 = 3430.0 * exp(224.93784961470999 * invT + 1.18 * logT - x45 - x72);
  x128 = sc[5] * x139;
  x72 = 0.00332 * exp(-2948.8496616156599 * invT + 2.81 * logT - x70 - x72);
  x101 = sc[8] * x72;
  x136 = exp(286.83349950869001 * invT);
  x174 = 32000000.0 * sc[8] * x136;
  x156 = 16000000.0 * sc[8] * x136;
  x192 = 28000000.0 *
         exp(g_RT[0] - g_RT[1] + g_RT[3] - g_RT[6] - 537.43539907944103 * invT);
  x21 = 0.05 * exp(-3165.23282791169 * invT + 2.67 * logT);
  x200 = 216.0 * sc[4] * exp(-1726.0331637101899 * invT + 1.51 * logT);
  x212 = g_RT[0] - 2.0 * g_RT[1];
  x217 = 0.5 * exp(-3638.25649376812 * invT + 2.0 * logT);
  x205 = sc[0] * x217;
  x168 = x376 * (x171 * x316 - x218 * x278 - x338 / x168);
  x171 = sc_qss[0] * x124;
  x218 = -x171 - x316;
  x338 = 660.0 * exp(
                   g_RT[0] - g_RT[1] + g_RT[7] - g_RT[8] -
                   5454.8686573231598 * invT + 1.62 * logT);
  x51 = (-sc_qss[0] * x67 / x51 + x278) / x55;
  x55 = 117906976744186.0 * x8 * exp(-2390.2791625723985 * invT - 4.92 * logT);
  x67 = log(
    exp(-10300.0 * invT) + 0.068 * exp(-0.0050761421319797002 * T) +
    0.932 * exp(-0.00064935064935064902 * T));
  x278 = -0.29097730287517898 * x67 + 0.43429448190325198 * log(x55) - 0.4;
  x376 = -0.14 * x278 - 0.55155399201712996 * x67 + 0.75;
  x389 = exp(
           M_LN10 * (0.43429448190325198 * x67 /
                     (((x278) * (x278)) / ((x376) * (x376)) + 1.0))) /
         (x55 + 1.0);
  amrex::Real x817 = -42446.325760628104 * invT - 3.42 * logT;
  const amrex::Real x819 = x283 * exp(-g_RT[0] + g_RT[11] - g_RT[9] + x817);
  const amrex::Real x821 = 5070000000000000.0 * x389 * x8 * x819;
  amrex::Real x823 = x125 * x278 / ((x376) * (x376));
  amrex::Real x824 = ((x278) * (x278)) / ((x376) * (x376) * (x376));
  const amrex::Real x826 =
    x67 / ((((x278) * (x278)) / ((x376) * (x376)) + 1.0) *
           (((x278) * (x278)) / ((x376) * (x376)) + 1.0));
  amrex::Real x827 = x826 * (0.56 * x125 * x824 + 4.0 * x823);
  const amrex::Real x831 =
    x283 * x8 *
    exp(
      -g_RT[0] + g_RT[11] - g_RT[9] - 44836.604923200604 * invT - 8.34 * logT);
  x55 = exp(
          M_LN10 * (0.43429448190325198 * x67 /
                    (((x278) * (x278)) / ((x376) * (x376)) + 1.0))) /
        ((x55 + 1.0) * (x55 + 1.0));
  x67 = sc[0] * sc[9] * exp(x817);
  x278 = 5070000000000000.0 * x389 * exp(x817);
  x376 = sc[0] * sc[9] * x278 * x8;
  const amrex::Real x844 =
    sc[0] * sc[9] * x8 * exp(-44836.604923200604 * invT - 8.34 * logT);
  const amrex::Real x848 = x376 * x827 - 1.014e+16 * x389 * x67 +
                           1.014e+16 * x389 * x819 -
                           1.19557674418605e+30 * x55 * x831 +
                           1.19557674418605e+30 * x55 * x844 - x821 * x827;
  const amrex::Real x855 =
    73400000.0 * sc[1] * x142 + 15000000.0 * sc[2] * x215 -
    15000000.0 * sc[2] * x220 - 15000000.0 * sc[9] * x274 - sc[9] * x278 * x8 -
    73400000.0 * sc[9] * x307;
  const amrex::Real x858 = 28000000.0 * exp(-537.43539907944103 * invT);
  const amrex::Real x860 =
    0.05 *
    exp(-g_RT[0] + g_RT[1] - 3165.23282791169 * invT + 2.67 * logT - x84);
  const amrex::Real x863 =
    660.0 * sc[8] * exp(-5454.8686573231598 * invT + 1.62 * logT);
  const amrex::Real x865 = exp(-logT);
  const amrex::Real x868 =
    23000.0 * exp(-1648.03458051046 * invT + 1.05 * logT);
  const amrex::Real x875 = 0.5 * exp(
                                   -g_RT[0] + g_RT[1] + g_RT[7] - g_RT_qss[0] -
                                   3638.25649376812 * invT + 2.0 * logT);
  const amrex::Real x882 =
    216.0 *
    exp(-g_RT[0] + g_RT[1] - 1726.0331637101899 * invT + 1.51 * logT - x45);
  const amrex::Real x884 = exp(-2.0 * logT);
  const amrex::Real x887 = refC * sc[0] * exp(-logT + x212);
  amrex::Real x893 = 0.28 * x125 * x824 + 2.0 * x823;
  const amrex::Real x901 =
    -sc[0] * sc[9] * x278 + x376 * x826 * x893 +
    5070000000000000.0 * x389 * x819 - 5.9778837209302301e+29 * x55 * x831 +
    5.9778837209302301e+29 * x55 * x844 - x821 * x826 * x893;
  const amrex::Real x902 =
    1000000.0 * ((sc[1]) * (sc[1])) * x865 - 1000000.0 * x887 + x901;
  const amrex::Real x908 =
    73400000.0 * sc[1] * x209 - 15000000.0 * sc[2] * x286 +
    15000000.0 * sc[2] * x407 - 73400000.0 * x171 - 73400000.0 * x316;
  const amrex::Real x909 = sc_qss[1] * x42;
  const amrex::Real x940 =
    216.0 * sc[0] * exp(-1726.0331637101899 * invT + 1.51 * logT);
  const amrex::Real x964 = x826 * (1.68 * x125 * x824 + 12.0 * x823);
  const amrex::Real x970 =
    73400000.0 * sc[1] * x95 + 15000000.0 * sc[2] * x346 -
    15000000.0 * sc[2] * x349 + x376 * x964 - 3.042e+16 * x389 * x67 +
    3.042e+16 * x389 * x819 - 3.5867302325581401e+30 * x55 * x831 +
    3.5867302325581401e+30 * x55 * x844 - x821 * x964;
  const amrex::Real x980 =
    70000000.0 * sc[0] * x106 - 15000000.0 * sc[2] * x106 + x902;
  const amrex::Real x999 =
    660.0 * sc[1] * exp(-5454.8686573231598 * invT + 1.62 * logT);
  x823 = 0.42 * x125 * x824 + 3.0 * x823;
  x824 = -15000000.0 * sc[0] * x274 - sc[0] * x278 * x8 -
         73400000.0 * sc[0] * x307 + 73400000.0 * sc[1] * x392 -
         15000000.0 * sc[2] * x27 + 15000000.0 * sc[2] * x359 +
         x376 * x823 * x826 - 7605000000000000.0 * x389 * x67 +
         7605000000000000.0 * x389 * x819 -
         8.9668255813953501e+29 * x55 * x831 +
         8.9668255813953501e+29 * x55 * x844 - x821 * x823 * x826;
  x817 = 5070000000000000.0 * refC * x389 * x8 *
         exp(-g_RT[0] + g_RT[11] - g_RT[9] + x817);
  amrex::Real x1063 =
    128000.0 * x8 * exp(-618.95649893980499 * invT - 2.34 * logT);
  amrex::Real x1064 =
    -0.29097730287517898 * x153 + 0.43429448190325198 * log(x1063) - 0.4;
  amrex::Real x1065 = -0.14 * x1064 - 0.55155399201712996 * x153 + 0.75;
  amrex::Real x1068 = ((x1064) * (x1064)) / ((x1065) * (x1065)) + 1.0;
  const amrex::Real x1071 =
    exp(M_LN10 * (0.43429448190325198 * x153 / x1068)) / (x1063 + 1.0);
  const amrex::Real x1074 =
    x168 * exp(-618.95649893980499 * invT - 3.14 * logT);
  const amrex::Real x1082 =
    ((x1064) * (x1064)) * x125 / ((x1065) * (x1065) * (x1065));
  const amrex::Real x1083 = x1064 * x125 / ((x1065) * (x1065));
  const amrex::Real x1086 =
    exp(-618.95649893980499 * invT - 3.14 * logT + x331);
  const amrex::Real x1088 = refC * sc[7] * x1071 * x1086;
  const amrex::Real x1091 = 3200000000000000.0 * x1088 * x153 * x8 *
                            (0.56 * x1082 + 4.0 * x1083) / ((x1068) * (x1068));
  amrex::Real x1094 =
    x409 / ((-0.14 * x409 - 0.55155399201712996 * x98 + 0.75) *
            (-0.14 * x409 - 0.55155399201712996 * x98 + 0.75));
  x409 =
    ((x409) * (x409)) / ((-0.14 * x409 - 0.55155399201712996 * x98 + 0.75) *
                         (-0.14 * x409 - 0.55155399201712996 * x98 + 0.75) *
                         (-0.14 * x409 - 0.55155399201712996 * x98 + 0.75));
  const amrex::Real x1096 = 4.0 * x1094 * x125 + 0.56 * x125 * x409;
  const amrex::Real x1097 = x98 / ((x415) * (x415));
  const amrex::Real x1098 = x218 * x423 * exp(x405);
  const amrex::Real x1100 = 1350000000000.0 * sc[1] * x1098 * x8;
  x1063 = exp(M_LN10 * (0.43429448190325198 * x153 / x1068)) /
          ((x1063 + 1.0) * (x1063 + 1.0));
  const amrex::Real x1105 =
    sc[1] * x168 * x8 * exp(-1237.91299787961 * invT - 5.48 * logT);
  x98 = exp(M_LN10 * (0.43429448190325198 * x98 / x415)) /
        ((x93 + 1.0) * (x93 + 1.0));
  x93 = -1565.003830652674 * invT - 5.62 * logT;
  x415 = 3.3440366972477102e+18 * sc[1] * x218 * x8 * x98 * exp(x93);
  const amrex::Real x1116 =
    3.3440366972477102e+18 * x283 * x8 * x98 * exp(x285 + x93);
  x331 = refC * sc[7] * x8 * exp(-1237.91299787961 * invT - 5.48 * logT + x331);
  x405 = x423 * exp(x285 + x405);
  x423 = 1350000000000.0 * x283 * x405 * x8;
  const amrex::Real x1133 = 3200000000000000.0 * sc[1] * x1071 * x1074 * x153 *
                            x8 * (0.56 * x1082 + 4.0 * x1083) /
                            ((x1068) * (x1068));
  const amrex::Real x1140 = 3200000000000000.0 * sc[1] * x1071 * x8 *
                            exp(-618.95649893980499 * invT - 3.14 * logT);
  const amrex::Real x1148 = exp(-g_RT[1] - 0.86 * logT - x87);
  const amrex::Real x1149 = refC * sc[6];
  const amrex::Real x1152 = exp(-0.86 * logT);
  const amrex::Real x1153 = sc[1] * sc[3];
  const amrex::Real x1156 =
    2800000.0 * x1148 * x1149 - 2800000.0 * x1152 * x1153;
  const amrex::Real x1157 = sc[1] * sc[2] * x865;
  const amrex::Real x1160 = refC * sc[4] * exp(-g_RT[1] - logT - x84);
  const amrex::Real x1162 = -1000000.0 * x1157 + 1000000.0 * x1160;
  amrex::Real x1163 = exp(g_RT[1] - 8554.6833186802305 * invT - logT - x279);
  const amrex::Real x1166 =
    -374000000000.0 * sc[9] * x1163 * x258 + 374000000000.0 * x218 * x434;
  const amrex::Real x1169 = sc[1] * sc[4] * x884;
  const amrex::Real x1171 = exp(-g_RT[1] - 2.0 * logT - x45);
  const amrex::Real x1172 = refC * sc[5] * x1171;
  const amrex::Real x1174 = sc[2] * x21 + 28000000.0 * sc[3] * x215 -
                            28000000.0 * sc[3] * x220 - 16060000000.0 * x1169 +
                            16060000000.0 * x1172;
  const amrex::Real x1179 =
    refC * sc[8] * exp(-g_RT[1] - 2262.965346123824 * invT - 8.89 * logT - x70);
  amrex::Real x1182 =
    195039370078.74 * x8 * exp(-1035.1166815603121 * invT - 4.13 * logT);
  amrex::Real x1183 = log(
    exp(-6964.0 * invT) + 0.217 * exp(-0.0135135135135135 * T) +
    0.783 * exp(-0.00034002040122407302 * T));
  amrex::Real x1185 =
    -0.29097730287517898 * x1183 + 0.43429448190325198 * log(x1182) - 0.4;
  amrex::Real x1186 = -0.55155399201712996 * x1183 - 0.14 * x1185 + 0.75;
  const amrex::Real x1192 =
    exp(
      M_LN10 * (0.43429448190325198 * x1183 /
                (((x1185) * (x1185)) / ((x1186) * (x1186)) + 1.0))) *
    x8 / ((x1182 + 1.0) * (x1182 + 1.0));
  x1182 = exp(
            M_LN10 * (0.43429448190325198 * x1183 /
                      (((x1185) * (x1185)) / ((x1186) * (x1186)) + 1.0))) /
          (x1182 + 1.0);
  const amrex::Real x1196 = -1227.8486645635121 * invT - 4.76 * logT;
  const amrex::Real x1197 = exp(-g_RT[1] + x1196 - x70);
  const amrex::Real x1201 =
    sc[1] * sc[7] * exp(-2262.965346123824 * invT - 8.89 * logT);
  amrex::Real x1204 = ((x1185) * (x1185)) / ((x1186) * (x1186) * (x1186));
  amrex::Real x1207 =
    4.0 * x1185 * x125 / ((x1186) * (x1186)) + 0.56 * x1204 * x125;
  x1183 = x1183 / ((((x1185) * (x1185)) / ((x1186) * (x1186)) + 1.0) *
                   (((x1185) * (x1185)) / ((x1186) * (x1186)) + 1.0));
  const amrex::Real x1211 = sc[1] * sc[7] * exp(x1196);
  const amrex::Real x1213 = 2.477e+21 * x1182 * x8;
  const amrex::Real x1215 = x1183 * x1211 * x1213;
  const amrex::Real x1218 = 2.477e+21 * refC * sc[8] * x1182 * x1197;
  const amrex::Real x1220 = x1183 * x1218 * x8;
  amrex::Real x1224 =
    4.954e+21 * refC * sc[8] * x1182 * x1197 -
    9.6622503937007901e+32 * x1179 * x1192 - 4.954e+21 * x1182 * x1211 +
    9.6622503937007901e+32 * x1192 * x1201 + x1207 * x1215 - x1207 * x1220;
  const amrex::Real x1225 = 70000000.0 * sc[0] * x215 -
                            70000000.0 * sc[0] * x220 + x1224 + x168 * x217 -
                            x205 * x342 + x205 * x377 + 70000000.0 * x51;
  const amrex::Real x1232 = -x245 - x909;
  amrex::Real x1233 = refC * x1232 * exp(x11 + x463);
  x12 = exp(
          M_LN10 * (0.43429448190325198 * x18 /
                    (1.0 + ((x449) * (x449)) / ((x32) * (x32))))) /
        ((x12 + 1.0) * (x12 + 1.0));
  amrex::Real x1237 =
    sc[11] * sc[1] * x8 * exp(-4760.4296585126494 * invT - 10.094 * logT);
  amrex::Real x1239 = x449 / ((x32) * (x32));
  amrex::Real x1240 = ((x449) * (x449)) / ((x32) * (x32) * (x32));
  x18 = x18 / ((1.0 + ((x449) * (x449)) / ((x32) * (x32))) *
               (1.0 + ((x449) * (x449)) / ((x32) * (x32))));
  x449 = sc[11] * exp(x463);
  x32 =
    refC * x1232 * x8 * exp(-4760.4296585126494 * invT - 10.094 * logT + x11);
  amrex::Real x1250 =
    sc[1] * x18 * x449 * x458 * (4.0 * x1239 * x125 + 0.56 * x1240 * x125) -
    2.54e+20 * sc[1] * x439 * x449 + 5.9737037037037004e+34 * x12 * x1237 -
    5.9737037037037004e+34 * x12 * x32 -
    x1233 * x18 * x458 * (4.0 * x1239 * x125 + 0.56 * x1240 * x125) +
    2.54e+20 * x1233 * x439;
  const amrex::Real x1251 =
    30000000.0 * sc[4] * x215 - 30000000.0 * sc[4] * x220 -
    20000000.0 * sc[4] * x342 + 20000000.0 * sc[4] * x377 + x1250 +
    x142 * x372 + x218 * x91 + x315 * x468;
  const amrex::Real x1259 = 28000000.0 * refCinv * sc[4] * sc[9] * x257;
  const amrex::Real x1276 = 47.6 * exp(
                                     g_RT[10] + g_RT[1] - g_RT[4] - g_RT[9] -
                                     35.225166606330397 * invT + 1.228 * logT);
  const amrex::Real x1283 = 3970000.0 * exp(-337.65838275496702 * invT);
  const amrex::Real x1285 =
    -sc[6] * x1283 - sc[6] * x858 -
    134000000.0 * sc[6] * exp(-319.54258278599701 * invT);
  const amrex::Real x1287 = 83000000.0 * exp(-7252.86180424342 * invT);
  const amrex::Real x1290 = 9380000.0 * exp(-0.76 * logT);
  const amrex::Real x1293 = exp(-1.72 * logT);
  const amrex::Real x1297 = 2800000.0 * sc[0] + 4200000.0 * sc[10] +
                            2800000.0 * sc[11] + 2800000.0 * sc[1] +
                            2800000.0 * sc[2] + 2800000.0 * sc[4] +
                            2800000.0 * sc[6] + 2800000.0 * sc[7] +
                            2800000.0 * sc[8] + 2100000.0 * sc[9];
  const amrex::Real x1301 = -375000000.0 * sc[12] * sc[3] * x1293 -
                            300000000.0 * ((sc[3]) * (sc[3])) * x1293 -
                            sc[3] * sc[5] * x1290 - sc[3] * x1152 * x1297 -
                            sc[3] * x1287;
  const amrex::Real x1303 = -2244000000000.0 * refCinv * sc[5] * sc[9] * x1163 -
                            187000000000.0 * refCinv * sc[9] * x1163 * x7 +
                            x209 * x437 + x209 * x442;
  const amrex::Real x1305 =
    -70000000.0 * sc[0] * x286 + 70000000.0 * sc[0] * x407 -
    sc[7] * x1213 * exp(x1196) - 70000000.0 * sc[7] * x248 - sc[7] * x875 +
    x150 * x205 - x205 * x393;
  const amrex::Real x1306 = 16060000000.0 * sc[0] + 22000000000.0 * sc[10] +
                            22000000000.0 * sc[3] + 80300000000.0 * sc[5] +
                            22000000000.0 * sc[9] + 22000000000.0 * x3;
  const amrex::Real x1308 = -28000000.0 * sc[3] * x286 +
                            28000000.0 * sc[3] * x407 - sc[4] * x1306 * x884 -
                            sc[4] * x860;
  const amrex::Real x1310 =
    84300000.0 * exp(g_RT[11] + g_RT[1] - g_RT[2] - g_RT[7]);
  const amrex::Real x1319 =
    -sc[11] * x1310 - 30000000.0 * sc[11] * x267 - 20000000.0 * sc[11] * x328 -
    sc[11] * x868 + 20000000.0 * sc[4] * x150 - 30000000.0 * sc[4] * x286 -
    20000000.0 * sc[4] * x393 + 30000000.0 * sc[4] * x407 + x16 * x468 +
    x209 * x372 - x449 * x458;
  amrex::Real x1320 =
    2.0 * x1185 * x125 / ((x1186) * (x1186)) + 0.28 * x1204 * x125;
  const amrex::Real x1327 = -4.83112519685039e+32 * x1179 * x1192 -
                            2.477e+21 * x1182 * x1211 +
                            4.83112519685039e+32 * x1192 * x1201 +
                            x1215 * x1320 + x1218 - x1220 * x1320;
  const amrex::Real x1335 =
    sc[1] * x18 * x449 * x458 * (2.0 * x1239 * x125 + 0.28 * x1240 * x125) -
    1.27e+20 * sc[1] * x439 * x449 + 2.9868518518518502e+34 * x12 * x1237 -
    2.9868518518518502e+34 * x12 * x32 -
    x1233 * x18 * x458 * (2.0 * x1239 * x125 + 0.28 * x1240 * x125) +
    1.27e+20 * x1233 * x439;
  const amrex::Real x1336 = x1097 * (2.0 * x1094 * x125 + 0.28 * x125 * x409);
  x285 = x283 * x8 * x98 * exp(x285 + x93);
  x98 = sc[1] * x218 * x8 * x98 * exp(x93);
  x93 = 3200000000000000.0 * sc[1] * x1071 * x1074 * x153 * x8 /
        ((x1068) * (x1068));
  x153 = 3200000000000000.0 * x1088 * x153 * x8 / ((x1068) * (x1068));
  x1068 = -187000000000.0 * sc[9] * x1163 * x258 + 187000000000.0 * x218 * x434;
  const amrex::Real x1360 = -22000000000.0 * x1169 + 22000000000.0 * x1172;
  const amrex::Real x1362 =
    -2000000.0 * ((sc[1]) * (sc[1])) * x865 -
    3200000000000000.0 * sc[1] * x1071 * x1074 -
    1350000000000.0 * sc[1] * x1098 + 4.096e+20 * x1063 * x1105 -
    4.096e+20 * x1063 * x331 + x1068 + 3200000000000000.0 * x1088 +
    x1100 * x1336 - 500000.0 * x1157 + 500000.0 * x1160 - x1336 * x423 + x1360 -
    x153 * (0.28 * x1082 + 2.0 * x1083) + 1350000000000.0 * x283 * x405 -
    1.67201834862385e+18 * x285 + 2000000.0 * x887 +
    x93 * (0.28 * x1082 + 2.0 * x1083) + 1.67201834862385e+18 * x98;
  const amrex::Real x1367 =
    3970000.0 *
    exp(-g_RT[1] + g_RT[2] + g_RT[5] - g_RT[6] - 337.65838275496702 * invT);
  const amrex::Real x1369 =
    83000000.0 *
    exp(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[4] - 7252.86180424342 * invT);
  const amrex::Real x1382 =
    30000000.0 * sc[2] * x277 - 30000000.0 * x171 - 30000000.0 * x316;
  const amrex::Real x1391 =
    30000000.0 * sc[4] * x290 - 30000000.0 * sc[4] * x335 -
    20000000.0 * sc[4] * x416 + 20000000.0 * sc[4] * x424 + 84300000.0 * sc[7] +
    x277 * x372 + x30 * x468;
  const amrex::Real x1406 = exp(-g_RT[1] - 1.72 * logT - x87);
  const amrex::Real x1411 = -375000000.0 * sc[12] * sc[1] * x1293 -
                            sc[1] * sc[5] * x1290 - sc[1] * x1152 * x1297 -
                            sc[1] * x1287 + 300000000.0 * x1149 * x1406 -
                            600000000.0 * x1153 * x1293;
  const amrex::Real x1418 =
    30000000.0 * sc[4] * x341 - 30000000.0 * sc[4] * x343 -
    20000000.0 * sc[4] * x430 + 20000000.0 * sc[4] * x435 + x372 * x374 +
    x382 * x468;
  const amrex::Real x1420 =
    47.6 * sc[9] * exp(-35.225166606330397 * invT + 1.228 * logT);
  const amrex::Real x1422 =
    sc[4] * exp(-g_RT[1] + 2.0 * g_RT[4] - g_RT[6] - 319.54258278599701 * invT);
  const amrex::Real x1425 = 500000.0 * refC * x8 * exp(-g_RT[1] - logT - x84);
  const amrex::Real x1447 =
    20000000.0 * sc[4] * x117 + 30000000.0 * sc[4] * x161 -
    30000000.0 * sc[4] * x345 - 20000000.0 * sc[4] * x443 + 20000000.0 * x168 +
    x372 * x375 + x390 * x468 + 30000000.0 * x51;
  const amrex::Real x1449 = -sc[1] * x1306 * x884 - sc[1] * x860 +
                            28000000.0 * sc[3] * x161 -
                            28000000.0 * sc[3] * x345;
  x1064 = 1.68 * ((x1064) * (x1064)) * x125 / ((x1065) * (x1065) * (x1065)) +
          12.0 * x1083;
  x1065 = 12.0 * x1094 * x125 + 1.68 * x125 * x409;
  const amrex::Real x1476 = 9380000.0 * exp(-g_RT[1] - 0.76 * logT - x87);
  const amrex::Real x1485 =
    30000000.0 * sc[2] * x95 - 2244000000000.0 * sc[9] * x1163 * x258 +
    2244000000000.0 * x218 * x434 + x437 * x95 + x442 * x95;
  const amrex::Real x1492 =
    12.0 * x1185 * x125 / ((x1186) * (x1186)) + 1.68 * x1204 * x125;
  const amrex::Real x1495 =
    1.4862e+22 * refC * sc[8] * x1182 * x1197 + 70000000.0 * sc[0] * x346 -
    70000000.0 * sc[0] * x349 - 2.8986751181102398e+33 * x1179 * x1192 -
    1.4862e+22 * x1182 * x1211 + 2.8986751181102398e+33 * x1192 * x1201 +
    x1215 * x1492 - x1220 * x1492 + x170 * x205 - x205 * x466;
  const amrex::Real x1499 = refC * x1171 * x1306 + 28000000.0 * sc[3] * x346 -
                            28000000.0 * sc[3] * x349 - 80300000000.0 * x1169 +
                            80300000000.0 * x1172;
  const amrex::Real x1506 =
    sc[1] * x18 * x449 * x458 * (12.0 * x1239 * x125 + 1.68 * x1240 * x125) -
    7.62e+20 * sc[1] * x439 * x449 + 20000000.0 * sc[4] * x170 +
    30000000.0 * sc[4] * x346 - 30000000.0 * sc[4] * x349 -
    20000000.0 * sc[4] * x466 + 1.7921111111111101e+35 * x12 * x1237 -
    1.7921111111111101e+35 * x12 * x32 -
    x1233 * x18 * x458 * (12.0 * x1239 * x125 + 1.68 * x1240 * x125) +
    7.62e+20 * x1233 * x439 + x372 * x95 + x395 * x468;
  const amrex::Real x1515 =
    -sc[1] * x1283 - sc[1] * x858 -
    134000000.0 * sc[1] * exp(-319.54258278599701 * invT);
  const amrex::Real x1521 = 375000000.0 * refC * sc[12] * x1406 +
                            300000000.0 * refC * sc[3] * x1406 +
                            refC * sc[5] * x1476 + refC * x1148 * x1297;
  const amrex::Real x1525 = -28000000.0 * sc[3] * x106;
  const amrex::Real x1530 =
    -20000000.0 * sc[4] * x478 + 30000000.0 * sc[4] * x486 +
    20000000.0 * sc[4] * x487 + x224 * x468 + x281 * x372;
  const amrex::Real x1535 = -70000000.0 * sc[0] * x106 -
                            30000000.0 * sc[4] * x106 - x106 * x428 + x1156 +
                            x1327 + x1335 + x1362;
  const amrex::Real x1553 =
    84300000.0 * sc[2] + 20000000.0 * sc[4] * x320 - 20000000.0 * sc[4] * x490 +
    30000000.0 * sc[4] * x495 + x372 * x385 + x399 * x468;
  const amrex::Real x1555 =
    70000000.0 * sc[0] * x495 - sc[1] * x1213 * exp(x1196) -
    70000000.0 * sc[1] * x248 - sc[1] * x875 + x205 * x320 - x205 * x490;
  const amrex::Real x1575 =
    30000000.0 * sc[4] * x352 - 30000000.0 * sc[4] * x356 -
    20000000.0 * sc[4] * x501 + 20000000.0 * sc[4] * x508 + x1250 +
    x202 * x468 + x372 * x386;
  const amrex::Real x1577 =
    2.477e+21 * refC * x1182 * x1197 * x8 + 70000000.0 * sc[0] * x352 -
    70000000.0 * sc[0] * x356 + x1224 - x205 * x501 + x205 * x508;
  const amrex::Real x1578 =
    -6400000000000000.0 * sc[1] * x1071 * x1074 -
    2700000000000.0 * sc[1] * x1098 - 38000000.0 * sc[1] * x202 +
    8.192e+20 * x1063 * x1105 - 8.192e+20 * x1063 * x331 +
    6400000000000000.0 * x1088 - x1091 + x1096 * x1097 * x1100 -
    x1096 * x1097 * x423 - x1116 + x1133 + x1162 + x1166 +
    2700000000000.0 * x283 * x405 + x415;
  const amrex::Real x1579 =
    47.6 * sc[4] * exp(-35.225166606330397 * invT + 1.228 * logT);
  x1094 = x1097 * (3.0 * x1094 * x125 + 0.42 * x125 * x409);
  x409 =
    x1183 * (3.0 * x1185 * x125 / ((x1186) * (x1186)) + 0.42 * x1204 * x125);
  x1185 = 3.7155e+21 * refC * sc[8] * x1182 * x1197 - 70000000.0 * sc[0] * x27 +
          70000000.0 * sc[0] * x359 - 7.2466877952755894e+32 * x1179 * x1192 -
          3.7155e+21 * x1182 * x1211 + 7.2466877952755894e+32 * x1192 * x1201 +
          x1211 * x1213 * x409 - x1218 * x409 * x8 - x148 * x205 + x164 * x205;
  x11 = 1.905e+20 * refC * x1232 * x439 * exp(x11 + x463) +
        sc[1] * x18 * x449 * x458 * (3.0 * x1239 * x125 + 0.42 * x1240 * x125) -
        1.905e+20 * sc[1] * x439 * x449 - 20000000.0 * sc[4] * x148 +
        20000000.0 * sc[4] * x164 - 30000000.0 * sc[4] * x27 +
        30000000.0 * sc[4] * x359 + 4.4802777777777799e+34 * x12 * x1237 -
        4.4802777777777799e+34 * x12 * x32 -
        x1233 * x18 * x458 * (3.0 * x1239 * x125 + 0.42 * x1240 * x125) +
        x204 * x468 + x372 * x392;
  x125 = -2244000000000.0 * refCinv * sc[1] * sc[5] * x1163 +
         30000000.0 * sc[2] * x392 - 280500000000.0 * sc[9] * x1163 * x258 -
         187000000000.0 * x1163 * x258 * x7 + 280500000000.0 * x218 * x434 +
         x392 * x437 + x392 * x442;
  x439 = 70000000.0 * sc[0] * x367 - 70000000.0 * sc[0] * x370 + x1224 -
         x178 * x205 + x179 * x205;
  x1163 = -20000000.0 * sc[4] * x178 + 20000000.0 * sc[4] * x179 +
          30000000.0 * sc[4] * x367 - 30000000.0 * sc[4] * x370 + x1250 +
          x141 * x372 + x202 * x468;
  x458 = -sc[1] * x1310 - 30000000.0 * sc[1] * x267 -
         20000000.0 * sc[1] * x328 - sc[1] * x458 * exp(x463) - sc[1] * x868 -
         20000000.0 * sc[4] * x186 + 30000000.0 * sc[4] * x196 +
         20000000.0 * sc[4] * x197 + x244 * x468 + x354 * x372;
  x463 = 70000000.0 * sc[0] * x54 - 70000000.0 * sc[0] * x58 + x1327 -
         x172 * x205 + x190 * x205;
  x1186 = -20000000.0 * sc[4] * x172 + 20000000.0 * sc[4] * x190 +
          30000000.0 * sc[4] * x54 - 30000000.0 * sc[4] * x58 + x119 * x372 +
          x1335 + x34 * x468;
  x1204 = refC * exp(-2.0 * g_RT[2] + g_RT[3] - logT);
  x1183 = sc[3] * x1204;
  x1224 = ((sc[2]) * (sc[2])) * x865;
  x1233 = exp(g_RT[10] - g_RT[2] - g_RT[9] - 1509.64999741416 * invT);
  x12 = refC * sc[10] * x1233;
  x1237 = exp(-1509.64999741416 * invT);
  x1239 = sc[2] * sc[9] * x1237;
  x1240 = -500000.0 * x1157 + 500000.0 * x1160 + 240000.0 * x1183 +
          602.0 * x12 - 240000.0 * x1224 - 602.0 * x1239;
  x18 = 2500000.0 *
        exp(g_RT[10] + g_RT[2] - g_RT[3] - g_RT[9] - 24053.756625465601 * invT);
  x449 = 1020.0 * exp(-4327.6633259205901 * invT + 1.5 * logT);
  x32 = 1204.0 * sc[0] + 2107.0 * sc[10] + 3612.0 * sc[3] + 3612.0 * sc[5] +
        903.0 * sc[9] + 602.0 * x3;
  x1250 = 2.4 * sc[0] + 3.6 * sc[10] + sc[3] + 15.4 * sc[5] + 1.75 * sc[9] + x3;
  const amrex::Real x1736 = exp(
    g_RT[2] - 2.0 * g_RT[4] + g_RT[5] + 1061.78716484796 * invT + 2.4 * logT);
  const amrex::Real x1740 = 39000000.0 * exp(-1781.38699694871 * invT);
  const amrex::Real x1744 = 2500000.0 * sc[9] * exp(-24053.756625465601 * invT);
  const amrex::Real x1746 = 20000000.0 * exp(g_RT[3] - g_RT[6] - x84);
  x84 = 1020.0 *
        exp(g_RT[7] - g_RT[8] - 4327.6633259205901 * invT + 1.5 * logT - x84);
  const amrex::Real x1762 = sc[4] * exp(1061.78716484796 * invT + 2.4 * logT);
  const amrex::Real x1785 = 2500000.0 * sc[3] * exp(-24053.756625465601 * invT);
  const amrex::Real x1811 = 13200000.0 * sc[3] * exp(-754.82499870708 * invT);
  const amrex::Real x1825 = -120000.0 * x1183 + 120000.0 * x1224;
  x70 = 1000000.0 * exp(g_RT[3] - g_RT[6] - x70);
  const amrex::Real x1839 =
    29000000.0 * exp(g_RT[3] - g_RT[6] + 251.60833290235999 * invT - x45);
  const amrex::Real x1854 = 36000.0 * exp(-4498.7569922941902 * invT);
  const amrex::Real x1860 = 100000000.0 * exp(-20128.666632188801 * invT);
  const amrex::Real x1864 =
    36000.0 *
    exp(g_RT[11] - g_RT[3] + g_RT[4] - g_RT[7] - 4498.7569922941902 * invT);
  const amrex::Real x1875 = 29000000.0 * sc[6] * exp(251.60833290235999 * invT);
  const amrex::Real x1889 = 29000000.0 * sc[4] * exp(251.60833290235999 * invT);
  const amrex::Real x1890 =
    18000000.0 * exp(g_RT[11] - g_RT_qss[3] - 452.89499922424801 * invT - x87);
  x87 = 7600000.0 * exp(-201.28666632188799 * invT - x279 - x87);
  x279 =
    56.0 *
    exp(-g_RT[7] + g_RT_qss[0] - 2727.4343286615799 * invT + 1.6 * logT - x45);
  const amrex::Real x1973 = sc[5] * x279;
  const amrex::Real x1978 = 12000000.0 * sc[1] * x315 + x133 * x215 -
                            x133 * x220 - x1973 * x342 + x1973 * x377;
  const amrex::Real x1982 = 500000.0 * x1157 - 500000.0 * x1160 + x1360;
  const amrex::Real x1988 = 12000000.0 * sc[1] * x16 - x133 * x286 +
                            x133 * x407 + x150 * x1973 - x1973 * x393 -
                            12000000.0 * x245 - 12000000.0 * x909;
  const amrex::Real x2000 = 12000000.0 * sc[1] * x30 + sc[8] * x449 +
                            x133 * x290 - x133 * x335 - x1973 * x416 +
                            x1973 * x424;
  const amrex::Real x2004 = 10000000.0 * sc[2] * x30 -
                            20000000.0 * sc[6] * x416 +
                            20000000.0 * sc[6] * x424 + x128 * x277 -
                            10000000.0 * x245 - 10000000.0 * x909;
  const amrex::Real x2011 = 12000000.0 * sc[1] * x382 + x133 * x341 -
                            x133 * x343 - x1973 * x430 + x1973 * x435;
  const amrex::Real x2015 =
    10000000.0 * sc[2] * x382 - 20000000.0 * sc[6] * x430 +
    20000000.0 * sc[6] * x435 + sc[7] * x1854 + x128 * x374;
  const amrex::Real x2019 =
    100.0 * sc[8] * exp(-1570.03599731073 * invT + 1.6 * logT);
  const amrex::Real x2029 =
    150000000.0 *
    exp(g_RT[10] + g_RT[4] - g_RT[6] - g_RT[9] - 11875.913312991401 * invT);
  const amrex::Real x2034 = 56.0 * exp(-2727.4343286615799 * invT + 1.6 * logT);
  const amrex::Real x2040 = 12000000.0 * sc[1] * x390 - sc[7] * x2034 -
                            12000000.0 * sc[7] * x231 - sc[7] * x84 -
                            25010000.0 * sc[7] + x117 * x1973 + x133 * x161 -
                            x133 * x345 - x1973 * x443;
  const amrex::Real x2044 =
    3430.0 * exp(224.93784961470999 * invT + 1.18 * logT);
  const amrex::Real x2047 =
    -sc[11] * x1864 - sc[11] * x2044 - 10000000.0 * sc[11] * x238 -
    20000000.0 * sc[11] * x324 + 10000000.0 * sc[2] * x390 +
    20000000.0 * sc[6] * x117 - 20000000.0 * sc[6] * x443 + x128 * x375;
  x45 =
    100.0 * exp(g_RT[7] - g_RT[8] - 1570.03599731073 * invT + 1.6 * logT - x45);
  const amrex::Real x2062 = 12000000.0 * sc[1] * x395 + x133 * x346 -
                            x133 * x349 + x168 * x279 + x170 * x1973 -
                            x1973 * x466 + 25010000.0 * x47 * x51;
  const amrex::Real x2067 =
    10000000.0 * sc[2] * x395 + 20000000.0 * sc[6] * x170 -
    20000000.0 * sc[6] * x466 + x128 * x95 + x139 * x218;
  const amrex::Real x2071 = 150000000.0 * exp(-11875.913312991401 * invT);
  const amrex::Real x2079 = 12000000.0 * sc[1] * x224 - x106 * x133 +
                            x133 * x486 - x1973 * x478 + x1973 * x487;
  const amrex::Real x2084 =
    10000000.0 * sc[2] * x224 - 20000000.0 * sc[6] * x478 +
    20000000.0 * sc[6] * x487 + x128 * x281 + 20000000.0 * x168;
  const amrex::Real x2092 = 10000000.0 * sc[2] * x399 + sc[3] * x1854 +
                            20000000.0 * sc[6] * x320 -
                            20000000.0 * sc[6] * x490 + x128 * x385;
  const amrex::Real x2098 = 12000000.0 * sc[1] * x399 - sc[4] * x2034 -
                            12000000.0 * sc[4] * x231 - sc[4] * x84 -
                            25010000.0 * sc[4] + x133 * x495 + x1973 * x320 -
                            x1973 * x490;
  const amrex::Real x2099 =
    100.0 * sc[4] * exp(-1570.03599731073 * invT + 1.6 * logT);
  const amrex::Real x2113 =
    sc[2] * x449 + x133 * x352 - x133 * x356 - x1973 * x501 + x1973 * x508;
  const amrex::Real x2123 = 12000000.0 * sc[1] * x204 - x133 * x27 +
                            x133 * x359 - x148 * x1973 + x164 * x1973;
  const amrex::Real x2153 =
    10000000.0 * sc[2] * x244 - sc[4] * x1864 - sc[4] * x2044 -
    10000000.0 * sc[4] * x238 - 20000000.0 * sc[4] * x324 -
    20000000.0 * sc[6] * x186 + 20000000.0 * sc[6] * x197 + x128 * x354;
  const amrex::Real x2164 = 12000000.0 * sc[1] * x34 + x133 * x54 - x133 * x58 -
                            x172 * x1973 + x190 * x1973;
  const amrex::Real x2165 = 22000000000.0 * x1169 - 22000000000.0 * x1172;
  const amrex::Real x2171 =
    -12000000.0 * sc[3] * x106 + x106 * x133 + x106 * x428 + x2165;
  const amrex::Real x2176 =
    -2800000.0 * x1148 * x1149 + 2800000.0 * x1152 * x1153;
  const amrex::Real x2187 = -4161.6018262050302 * invT + 2.0 * logT;
  const amrex::Real x2190 = 4.92 * sc[8] * exp(x2187);
  const amrex::Real x2193 =
    6400000000000000.0 * sc[1] * x1071 * x1074 - 8.192e+20 * x1063 * x1105 +
    8.192e+20 * x1063 * x331 - 6400000000000000.0 * x1088 + x1091 - x1133;
  const amrex::Real x2196 =
    3200000000000000.0 * sc[1] * x1071 * x1074 - 4.096e+20 * x1063 * x1105 +
    4.096e+20 * x1063 * x331 - 3200000000000000.0 * x1088 +
    x153 * (0.28 * x1082 + 2.0 * x1083) - x93 * (0.28 * x1082 + 2.0 * x1083);
  const amrex::Real x2208 =
    -70000000.0 * sc[0] * x106 - x106 * x174 + x1327 + x2196;
  const amrex::Real x2215 =
    0.00332 * exp(-2948.8496616156599 * invT + 2.81 * logT);
  const amrex::Real x2239 = 2.46 * sc[8] * exp(x2187);
  x1207 = -4.954e+21 * refC * sc[8] * x1182 * x1197 +
          9.6622503937007901e+32 * x1179 * x1192 + 4.954e+21 * x1182 * x1211 -
          9.6622503937007901e+32 * x1192 * x1201 - x1207 * x1215 +
          x1207 * x1220;
  x1320 = 4.83112519685039e+32 * x1179 * x1192 + 2.477e+21 * x1182 * x1211 -
          4.83112519685039e+32 * x1192 * x1201 - x1215 * x1320 - x1218 +
          x1220 * x1320;
  const amrex::Real x2252 = x1068 + 602.0 * x12 - 602.0 * x1239 + x901;
  const amrex::Real x2270 = -14000000.0 * sc[10] * x106;
  const amrex::Real x2271 =
    -15000000.0 * sc[2] * x106 - 40000000.0 * sc[3] * x106 + x2252;
  const amrex::Real x2294 = -602.0 * x12 + 602.0 * x1239;
  x827 = 2700000000000.0 * sc[1] * x1098 - x1096 * x1097 * x1100 +
         x1096 * x1097 * x423 + x1116 - 2700000000000.0 * x283 * x405 -
         x376 * x827 + 1.014e+16 * x389 * x67 - 1.014e+16 * x389 * x819 - x415 +
         1.19557674418605e+30 * x55 * x831 - 1.19557674418605e+30 * x55 * x844 +
         x821 * x827;
  x893 =
    sc[0] * sc[9] * x278 + 1350000000000.0 * sc[1] * x1098 - x1100 * x1336 +
    x1336 * x423 - 1350000000000.0 * x283 * x405 + 1.67201834862385e+18 * x285 -
    x376 * x826 * x893 - 5070000000000000.0 * x389 * x819 +
    5.9778837209302301e+29 * x55 * x831 - 5.9778837209302301e+29 * x55 * x844 +
    x821 * x826 * x893 - 1.67201834862385e+18 * x98;
  // Write the full Jacobian expression...
  J[0] = -180000.0 * refC * sc[0] * exp(-0.6 * logT + x212) -
         550000000.0 * refC * sc[10] * exp(-2.0 * logT + x212) -
         60000000.0 * refC * sc[5] * exp(-1.25 * logT + x212) -
         1000000.0 * refC * (sc[3] + sc[9] + x3) * exp(-logT + x212) -
         70000000.0 * sc[0] * x215 + 70000000.0 * sc[0] * x220 -
         20000000.0 * sc[11] * x235 +
         90000.0 * ((sc[1]) * (sc[1])) * exp(-0.6 * logT) +
         20000000.0 * sc[1] * x315 - sc[2] * x21 - sc[3] * x192 - sc[7] * x338 -
         x142 * x372 - x168 * x217 - x200 + x205 * x342 - x205 * x377 -
         x218 * x91 - 70000000.0 * x51 + x848 + x855 +
         (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
          15000000.0 * sc[2] * x69 + x124 * x372) *
           (-x342 + x377) +
         (-20000000.0 * sc[1] * x42) * (x215 - x220 + (-x69) * (-x342 + x377));
  J[14] =
    180000.0 * sc[0] * sc[1] * exp(-0.6 * logT) + 70000000.0 * sc[0] * x286 -
    70000000.0 * sc[0] * x407 + 1100000000.0 * sc[10] * sc[1] * x884 +
    sc[11] * x868 + 120000000.0 * sc[1] * sc[5] * exp(-1.25 * logT) +
    20000000.0 * sc[1] * x16 + 2000000.0 * sc[1] * x865 * (sc[3] + sc[9] + x3) +
    sc[4] * x860 + sc[5] * x882 + sc[6] * x858 + 70000000.0 * sc[7] * x248 +
    sc[7] * x875 - x150 * x205 + x205 * x393 - x209 * x372 - 20000000.0 * x245 +
    x863 + x902 + x908 - 20000000.0 * x909 +
    (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
     15000000.0 * sc[2] * x69 + x124 * x372) *
      (x150 - x393) +
    (-20000000.0 * sc[1] * x42) * (-x286 + x407 + (-x69) * (x150 - x393));
  J[28] = -sc[0] * x21 - 70000000.0 * sc[0] * x290 + 70000000.0 * sc[0] * x335 +
          73400000.0 * sc[1] * x277 + 20000000.0 * sc[1] * x30 +
          15000000.0 * sc[2] * x290 - 15000000.0 * sc[2] * x335 + x205 * x416 -
          x205 * x424 - x277 * x372 + 15000000.0 * x51 + x902 +
          (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
           15000000.0 * sc[2] * x69 + x124 * x372) *
            (-x416 + x424) +
          (-20000000.0 * sc[1] * x42) * (x290 - x335 + (-x69) * (-x416 + x424));
  J[42] = -sc[0] * x192 - 70000000.0 * sc[0] * x341 +
          70000000.0 * sc[0] * x343 + 73400000.0 * sc[1] * x374 +
          20000000.0 * sc[1] * x382 + 15000000.0 * sc[2] * x341 -
          15000000.0 * sc[2] * x343 + x205 * x430 - x205 * x435 - x372 * x374 +
          x902 +
          (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
           15000000.0 * sc[2] * x69 + x124 * x372) *
            (-x430 + x435) +
          (-20000000.0 * sc[1] * x42) * (x341 - x343 + (-x69) * (-x430 + x435));
  J[56] = -70000000.0 * sc[0] * x161 + 70000000.0 * sc[0] * x345 +
          73400000.0 * sc[1] * x375 + 20000000.0 * sc[1] * x390 + sc[1] * x860 +
          15000000.0 * sc[2] * x161 - 15000000.0 * sc[2] * x345 - x117 * x205 +
          x205 * x443 - x372 * x375 + x902 - x940 +
          (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
           15000000.0 * sc[2] * x69 + x124 * x372) *
            (x117 - x443) +
          (-20000000.0 * sc[1] * x42) * (x161 - x345 + (-x69) * (x117 - x443));
  J[70] = -60000000.0 * refC * sc[0] * exp(-1.25 * logT + x212) -
          70000000.0 * sc[0] * x346 + 70000000.0 * sc[0] * x349 +
          60000000.0 * ((sc[1]) * (sc[1])) * exp(-1.25 * logT) +
          20000000.0 * sc[1] * x395 + sc[1] * x882 - x170 * x205 + x205 * x466 -
          x372 * x95 + x970 +
          (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
           15000000.0 * sc[2] * x69 + x124 * x372) *
            (x170 - x466) +
          (-20000000.0 * sc[1] * x42) * (x346 - x349 + (-x69) * (x170 - x466));
  J[84] =
    -70000000.0 * sc[0] * x486 + 20000000.0 * sc[1] * x224 +
    73400000.0 * sc[1] * x281 + sc[1] * x858 + 15000000.0 * sc[2] * x486 +
    x205 * x478 - x205 * x487 - x281 * x372 + x980 +
    (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
     15000000.0 * sc[2] * x69 + x124 * x372) *
      (-x478 + x487) +
    (-20000000.0 * sc[1] * x42) * (-x106 + x486 + (-x69) * (-x478 + x487));
  J[98] = -sc[0] * x338 - 70000000.0 * sc[0] * x495 +
          70000000.0 * sc[1] * x248 + 73400000.0 * sc[1] * x385 +
          20000000.0 * sc[1] * x399 + sc[1] * x875 + 15000000.0 * sc[2] * x495 -
          x205 * x320 + x205 * x490 - x372 * x385 + x980 +
          (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
           15000000.0 * sc[2] * x69 + x124 * x372) *
            (x320 - x490) +
          (-20000000.0 * sc[1] * x42) * (-x106 + x495 + (-x69) * (x320 - x490));
  J[112] =
    -70000000.0 * sc[0] * x352 + 70000000.0 * sc[0] * x356 +
    2000000.0 * ((sc[1]) * (sc[1])) * x865 + 20000000.0 * sc[1] * x202 +
    73400000.0 * sc[1] * x386 + 15000000.0 * sc[2] * x352 -
    15000000.0 * sc[2] * x356 + x205 * x501 - x205 * x508 - x372 * x386 + x848 -
    2000000.0 * x887 + x999 +
    (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
     15000000.0 * sc[2] * x69 + x124 * x372) *
      (-x501 + x508) +
    (-20000000.0 * sc[1] * x42) * (x352 - x356 + (-x69) * (-x501 + x508));
  J[126] =
    70000000.0 * sc[0] * x27 - 70000000.0 * sc[0] * x359 +
    1000000.0 * ((sc[1]) * (sc[1])) * x865 + 20000000.0 * sc[1] * x204 +
    x148 * x205 - x164 * x205 - x372 * x392 + x824 - 1000000.0 * x887 +
    (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
     15000000.0 * sc[2] * x69 + x124 * x372) *
      (-x148 + x164) +
    (-20000000.0 * sc[1] * x42) * (-x27 + x359 + (-x69) * (-x148 + x164));
  J[140] =
    -550000000.0 * refC * sc[0] * exp(-2.0 * logT + x212) -
    70000000.0 * sc[0] * x367 + 70000000.0 * sc[0] * x370 +
    550000000.0 * ((sc[1]) * (sc[1])) * x884 + 73400000.0 * sc[1] * x141 +
    20000000.0 * sc[1] * x202 + 15000000.0 * sc[2] * x367 -
    15000000.0 * sc[2] * x370 - x141 * x372 + x178 * x205 - x179 * x205 + x848 +
    (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
     15000000.0 * sc[2] * x69 + x124 * x372) *
      (-x178 + x179) +
    (-20000000.0 * sc[1] * x42) * (x367 - x370 + (-x69) * (-x178 + x179));
  J[154] =
    -70000000.0 * sc[0] * x196 - 20000000.0 * sc[0] * x235 +
    20000000.0 * sc[1] * x244 + 73400000.0 * sc[1] * x354 + sc[1] * x868 +
    15000000.0 * sc[2] * x196 + x186 * x205 - x197 * x205 - x354 * x372 + x817 +
    x980 +
    (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
     15000000.0 * sc[2] * x69 + x124 * x372) *
      (-x186 + x197) +
    (-20000000.0 * sc[1] * x42) * (-x106 + x196 + (-x69) * (-x186 + x197));
  J[168] = -70000000.0 * sc[0] * x54 + 70000000.0 * sc[0] * x58 +
           73400000.0 * sc[1] * x119 + 20000000.0 * sc[1] * x34 +
           15000000.0 * sc[2] * x54 - 15000000.0 * sc[2] * x58 - x119 * x372 +
           x172 * x205 - x190 * x205 + x902 +
           (70000000.0 * sc[0] * x69 - 73400000.0 * sc[1] * x124 -
            15000000.0 * sc[2] * x69 + x124 * x372) *
             (-x172 + x190) +
           (-20000000.0 * sc[1] * x42) * (x54 - x58 + (-x69) * (-x172 + x190));
  J[1] = 360000.0 * refC * sc[0] * exp(-0.6 * logT + x212) +
         1100000000.0 * refC * sc[10] * exp(-2.0 * logT + x212) +
         120000000.0 * refC * sc[5] * exp(-1.25 * logT + x212) +
         2000000.0 * refC * (sc[3] + sc[9] + x3) * exp(-logT + x212) +
         20000000.0 * sc[11] * x235 -
         180000.0 * ((sc[1]) * (sc[1])) * exp(-0.6 * logT) -
         6400000000000000.0 * sc[1] * x1071 * x1074 -
         2700000000000.0 * sc[1] * x1098 - 73400000.0 * sc[1] * x142 -
         38000000.0 * sc[1] * x315 + 30000000.0 * sc[2] * x142 -
         80000000.0 * sc[2] * x342 + 80000000.0 * sc[2] * x377 + sc[3] * x192 +
         sc[7] * x338 + 73400000.0 * sc[9] * x307 + 8.192e+20 * x1063 * x1105 -
         8.192e+20 * x1063 * x331 + 6400000000000000.0 * x1088 - x1091 +
         x1096 * x1097 * x1100 - x1096 * x1097 * x423 - x1116 + x1133 +
         x1140 * x342 - x1140 * x377 + x1156 + x1162 + x1166 + x1174 + x1225 +
         x1251 - x142 * x400 - x142 * x427 + x142 * x437 + x142 * x442 + x200 +
         x215 * x428 - x220 * x428 + 2700000000000.0 * x283 * x405 + x415 +
         (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
          30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
          30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
          x124 * x437 - x124 * x442 - x428 * x69) *
           (-x342 + x377) +
         (38000000.0 * sc[1] * x42 - x42 * x468) *
           (x215 - x220 + (-x69) * (-x342 + x377));
  J[15] = -360000.0 * sc[0] * sc[1] * exp(-0.6 * logT) -
          2200000000.0 * sc[10] * sc[1] * x884 - sc[10] * x1276 -
          30000000.0 * sc[10] * x298 -
          240000000.0 * sc[1] * sc[5] * exp(-1.25 * logT) -
          38000000.0 * sc[1] * x16 - 73400000.0 * sc[1] * x209 -
          4000000.0 * sc[1] * x865 * (sc[3] + sc[9] + x3) +
          80000000.0 * sc[2] * x150 + 30000000.0 * sc[2] * x209 -
          80000000.0 * sc[2] * x393 - 500000.0 * sc[2] * x8 * x865 -
          sc[5] * x882 - 3200000000000000.0 * x1071 * x1074 * x8 -
          1350000000000.0 * x1098 * x8 - x1140 * x150 + x1140 * x393 + x1156 -
          x1259 + x1285 + x1301 + x1303 + x1305 + x1308 + x1319 + x1327 +
          x1335 + x1362 + 73400000.0 * x171 - x209 * x400 - x209 * x427 -
          80000000.0 * x218 * x79 + 38000000.0 * x245 - x286 * x428 +
          73400000.0 * x316 + x407 * x428 - x863 + 38000000.0 * x909 +
          (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
           30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
           30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
           x124 * x437 - x124 * x442 - x428 * x69) *
            (x150 - x393) +
          (38000000.0 * sc[1] * x42 - x42 * x468) *
            (-x286 + x407 + (-x69) * (x150 - x393));
  J[29] = sc[0] * x21 + 70000000.0 * sc[0] * x290 - 70000000.0 * sc[0] * x335 -
          73400000.0 * sc[1] * x277 - 38000000.0 * sc[1] * x30 -
          500000.0 * sc[1] * x8 * x865 - 80000000.0 * sc[2] * x416 +
          80000000.0 * sc[2] * x424 + 28000000.0 * sc[3] * x290 -
          28000000.0 * sc[3] * x335 + sc[4] * x1369 + sc[5] * x1367 +
          x1140 * x416 - x1140 * x424 + x1156 + x1327 + x1335 + x1362 + x1382 +
          x1391 + 80000000.0 * x168 - x205 * x416 + x205 * x424 - x277 * x400 -
          x277 * x427 + x277 * x437 + x277 * x442 + x290 * x428 - x335 * x428 +
          (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
           30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
           30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
           x124 * x437 - x124 * x442 - x428 * x69) *
            (-x416 + x424) +
          (38000000.0 * sc[1] * x42 - x42 * x468) *
            (x290 - x335 + (-x69) * (-x416 + x424));
  J[43] = sc[0] * x192 + 70000000.0 * sc[0] * x341 - 70000000.0 * sc[0] * x343 -
          73400000.0 * sc[1] * x374 - 38000000.0 * sc[1] * x382 +
          30000000.0 * sc[2] * x374 - 80000000.0 * sc[2] * x430 +
          80000000.0 * sc[2] * x435 + 28000000.0 * sc[3] * x341 -
          28000000.0 * sc[3] * x343 + x1140 * x430 - x1140 * x435 + x1327 +
          x1335 + x1362 + x1411 + x1418 - x205 * x430 + x205 * x435 +
          x341 * x428 - x343 * x428 - x374 * x400 - x374 * x427 + x374 * x437 +
          x374 * x442 + 28000000.0 * x51 +
          (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
           30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
           30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
           x124 * x437 - x124 * x442 - x428 * x69) *
            (-x430 + x435) +
          (38000000.0 * sc[1] * x42 - x42 * x468) *
            (x341 - x343 + (-x69) * (-x430 + x435));
  J[57] = 70000000.0 * sc[0] * x161 - 70000000.0 * sc[0] * x345 -
          73400000.0 * sc[1] * x375 - 38000000.0 * sc[1] * x390 +
          80000000.0 * sc[2] * x117 + sc[2] * x1369 +
          30000000.0 * sc[2] * x375 - 80000000.0 * sc[2] * x443 +
          12000000.0 * sc[7] * x231 - 28000000.0 * sc[9] * x257 * x258 -
          x1140 * x117 + x1140 * x443 + x1156 + x117 * x205 + x1327 + x1335 +
          x1362 + x1420 + 268000000.0 * x1422 + x1425 + x1447 + x1449 +
          x161 * x428 - x205 * x443 - x345 * x428 - x375 * x400 - x375 * x427 +
          x375 * x437 + x375 * x442 + x940 +
          (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
           30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
           30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
           x124 * x437 - x124 * x442 - x428 * x69) *
            (x117 - x443) +
          (38000000.0 * sc[1] * x42 - x42 * x468) *
            (x161 - x345 + (-x69) * (x117 - x443));
  J[71] = 120000000.0 * refC * sc[0] * exp(-1.25 * logT + x212) -
          120000000.0 * ((sc[1]) * (sc[1])) * exp(-1.25 * logT) -
          1.92e+16 * sc[1] * x1071 * x1074 - 8100000000000.0 * sc[1] * x1098 -
          38000000.0 * sc[1] * x395 - sc[1] * x882 - 73400000.0 * sc[1] * x95 +
          sc[2] * x1367 + 80000000.0 * sc[2] * x170 -
          80000000.0 * sc[2] * x466 + 2.4576e+21 * x1063 * x1105 -
          2.4576e+21 * x1063 * x331 - x1064 * x153 + x1064 * x93 +
          x1065 * x1097 * x1100 - x1065 * x1097 * x423 + 1.92e+16 * x1088 -
          x1140 * x170 + x1140 * x466 + x1149 * x1476 - x1153 * x1290 -
          3000000.0 * x1157 + 3000000.0 * x1160 + x1485 + x1495 + x1499 +
          x1506 + 8100000000000.0 * x283 * x405 -
          1.0032110091743101e+19 * x285 + x346 * x428 - x349 * x428 +
          6000000.0 * x38 * x51 - x400 * x95 - x427 * x95 +
          1.0032110091743101e+19 * x98 +
          (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
           30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
           30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
           x124 * x437 - x124 * x442 - x428 * x69) *
            (x170 - x466) +
          (38000000.0 * sc[1] * x42 - x42 * x468) *
            (x346 - x349 + (-x69) * (x170 - x466));
  J[85] = 70000000.0 * sc[0] * x486 - 38000000.0 * sc[1] * x224 -
          73400000.0 * sc[1] * x281 + 30000000.0 * sc[2] * x281 -
          80000000.0 * sc[2] * x478 + 80000000.0 * sc[2] * x487 +
          28000000.0 * sc[3] * x486 + x1140 * x478 - x1140 * x487 + x1515 +
          x1521 + x1525 + x1530 + x1535 - x205 * x478 + x205 * x487 -
          x281 * x400 - x281 * x427 + x281 * x437 + x281 * x442 + x428 * x486 +
          (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
           30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
           30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
           x124 * x437 - x124 * x442 - x428 * x69) *
            (-x478 + x487) +
          (38000000.0 * sc[1] * x42 - x42 * x468) *
            (-x106 + x486 + (-x69) * (-x478 + x487));
  J[99] = 3200000000000000.0 * refC * x1071 * x1086 * x8 + sc[0] * x338 -
          73400000.0 * sc[1] * x385 - 38000000.0 * sc[1] * x399 +
          80000000.0 * sc[2] * x320 + 30000000.0 * sc[2] * x385 -
          80000000.0 * sc[2] * x490 + 28000000.0 * sc[3] * x495 +
          12000000.0 * sc[4] * x231 + 2244000000000.0 * sc[5] * x385 * x434 -
          x1140 * x320 + x1140 * x490 + x1525 + x1535 + x1553 + x1555 -
          x385 * x400 - x385 * x427 + 187000000000.0 * x385 * x434 * x7 +
          x428 * x495 +
          (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
           30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
           30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
           x124 * x437 - x124 * x442 - x428 * x69) *
            (x320 - x490) +
          (38000000.0 * sc[1] * x42 - x42 * x468) *
            (-x106 + x495 + (-x69) * (x320 - x490));
  J[113] = -4000000.0 * ((sc[1]) * (sc[1])) * x865 - 73400000.0 * sc[1] * x386 +
           30000000.0 * sc[2] * x386 - 80000000.0 * sc[2] * x501 +
           80000000.0 * sc[2] * x508 + 28000000.0 * sc[3] * x352 -
           28000000.0 * sc[3] * x356 + x1140 * x501 - x1140 * x508 + x1156 -
           44000000000.0 * x1169 + 44000000000.0 * x1172 + x1575 + x1577 +
           x1578 + x352 * x428 - x356 * x428 - x386 * x400 - x386 * x427 +
           x386 * x437 + x386 * x442 + 4000000.0 * x887 - x999 +
           (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
            30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
            30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
            x124 * x437 - x124 * x442 - x428 * x69) *
             (-x501 + x508) +
           (38000000.0 * sc[1] * x42 - x42 * x468) *
             (x352 - x356 + (-x69) * (-x501 + x508));
  J[127] = 4800000000000000.0 * refC * sc[7] * x1071 * x1086 +
           73400000.0 * sc[0] * x307 - 2000000.0 * ((sc[1]) * (sc[1])) * x865 -
           4800000000000000.0 * sc[1] * x1071 * x1074 -
           2025000000000.0 * sc[1] * x1098 - 38000000.0 * sc[1] * x204 -
           73400000.0 * sc[1] * x392 - 80000000.0 * sc[2] * x148 +
           80000000.0 * sc[2] * x164 - 28000000.0 * sc[3] * x27 +
           28000000.0 * sc[3] * x359 - 28000000.0 * sc[4] * x257 * x258 +
           6.144e+20 * x1063 * x1105 - 6.144e+20 * x1063 * x331 +
           x1094 * x1100 - x1094 * x423 + x11 + x1140 * x148 - x1140 * x164 +
           2100000.0 * x1148 * x1149 - 2100000.0 * x1152 * x1153 -
           750000.0 * x1157 + 750000.0 * x1160 + x1185 + x125 + x1360 -
           x153 * (0.42 * x1082 + 3.0 * x1083) + x1579 - x27 * x428 +
           2025000000000.0 * x283 * x405 - 2.5080275229357798e+18 * x285 +
           x359 * x428 - x392 * x400 - x392 * x427 + 2000000.0 * x887 +
           x93 * (0.42 * x1082 + 3.0 * x1083) + 2.5080275229357798e+18 * x98 +
           (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
            30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
            30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
            x124 * x437 - x124 * x442 - x428 * x69) *
             (-x148 + x164) +
           (38000000.0 * sc[1] * x42 - x42 * x468) *
             (-x27 + x359 + (-x69) * (-x148 + x164));
  J[141] = 1100000000.0 * refC * sc[0] * exp(-2.0 * logT + x212) -
           1100000000.0 * ((sc[1]) * (sc[1])) * x884 - sc[1] * x1276 -
           73400000.0 * sc[1] * x141 - 30000000.0 * sc[1] * x298 +
           30000000.0 * sc[2] * x141 - 80000000.0 * sc[2] * x178 +
           80000000.0 * sc[2] * x179 + 28000000.0 * sc[3] * x367 -
           28000000.0 * sc[3] * x370 + x1140 * x178 - x1140 * x179 +
           4200000.0 * x1148 * x1149 - 4200000.0 * x1152 * x1153 + x1163 +
           x1360 - x141 * x400 - x141 * x427 + x141 * x437 + x141 * x442 +
           x1578 + x367 * x428 - x370 * x428 + x439 +
           (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
            30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
            30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
            x124 * x437 - x124 * x442 - x428 * x69) *
             (-x178 + x179) +
           (38000000.0 * sc[1] * x42 - x42 * x468) *
             (x367 - x370 + (-x69) * (-x178 + x179));
  J[155] = 1350000000000.0 * refC * x405 * x8 + 70000000.0 * sc[0] * x196 +
           20000000.0 * sc[0] * x235 - 38000000.0 * sc[1] * x244 -
           73400000.0 * sc[1] * x354 - 80000000.0 * sc[2] * x186 +
           80000000.0 * sc[2] * x197 + 30000000.0 * sc[2] * x354 +
           28000000.0 * sc[3] * x196 + x1140 * x186 - x1140 * x197 + x1525 +
           x1535 - x186 * x205 + x196 * x428 + x197 * x205 - x354 * x400 -
           x354 * x427 + x354 * x437 + x354 * x442 + x458 +
           (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
            30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
            30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
            x124 * x437 - x124 * x442 - x428 * x69) *
             (-x186 + x197) +
           (38000000.0 * sc[1] * x42 - x42 * x468) *
             (-x106 + x196 + (-x69) * (-x186 + x197));
  J[169] = -73400000.0 * sc[1] * x119 - 38000000.0 * sc[1] * x34 +
           30000000.0 * sc[2] * x119 - 80000000.0 * sc[2] * x172 +
           80000000.0 * sc[2] * x190 + 28000000.0 * sc[3] * x54 -
           28000000.0 * sc[3] * x58 + x1140 * x172 - x1140 * x190 +
           375000000.0 * x1149 * x1406 - 375000000.0 * x1153 * x1293 + x1186 -
           x119 * x400 - x119 * x427 + x119 * x437 + x119 * x442 + x1362 +
           x428 * x54 - x428 * x58 + x463 +
           (-70000000.0 * sc[0] * x69 + 73400000.0 * sc[1] * x124 -
            30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 -
            30000000.0 * sc[4] * x69 - x124 * x372 + x124 * x400 + x124 * x427 -
            x124 * x437 - x124 * x442 - x428 * x69) *
             (-x172 + x190) +
           (38000000.0 * sc[1] * x42 - x42 * x468) *
             (x54 - x58 + (-x69) * (-x172 + x190));
  J[2] = -60000000.0 * sc[2] * x142 - sc[2] * x21 - 15000000.0 * sc[2] * x215 +
         15000000.0 * sc[2] * x220 - 10000000.0 * sc[2] * x315 +
         80000000.0 * sc[2] * x342 - 80000000.0 * sc[2] * x377 +
         15000000.0 * sc[9] * x274 + x1162 + 576000.0 * x1183 + 1204.0 * x12 -
         576000.0 * x1224 - 1204.0 * x1239 + x142 * x400 + x142 * x473 +
         (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
          x124 * x473) *
           (-x342 + x377) +
         (10000000.0 * sc[2] * x42) * (x215 - x220 + (-x69) * (-x342 + x377));
  J[16] = 30000000.0 * sc[10] * x298 + sc[11] * x1310 -
          80000000.0 * sc[2] * x150 - 10000000.0 * sc[2] * x16 -
          60000000.0 * sc[2] * x209 + 15000000.0 * sc[2] * x286 +
          80000000.0 * sc[2] * x393 - 15000000.0 * sc[2] * x407 -
          500000.0 * sc[2] * x8 * x865 + sc[3] * x1287 + sc[4] * x860 +
          sc[6] * x1283 + x1240 + x209 * x400 + x209 * x473 +
          80000000.0 * x218 * x79 +
          (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
           x124 * x473) *
            (x150 - x393) +
          (10000000.0 * sc[2] * x42) * (-x286 + x407 + (-x69) * (x150 - x393));
  J[30] = -sc[0] * x21 - sc[10] * x18 - sc[11] * x1740 -
          500000.0 * sc[1] * x8 * x865 - 480000.0 * sc[2] * x1250 * x865 -
          60000000.0 * sc[2] * x277 - 15000000.0 * sc[2] * x290 -
          10000000.0 * sc[2] * x30 + 15000000.0 * sc[2] * x335 +
          80000000.0 * sc[2] * x416 - 80000000.0 * sc[2] * x424 -
          sc[4] * x1369 - sc[5] * x1367 - 0.0357 * sc[5] * x1736 -
          20000000.0 * sc[6] - 84300000.0 * sc[7] - sc[8] * x449 -
          sc[9] * x1237 * x32 + x1240 - 80000000.0 * x168 + 60000000.0 * x171 +
          10000000.0 * x245 + x277 * x400 + x277 * x473 + 60000000.0 * x316 -
          15000000.0 * x51 + 10000000.0 * x909 +
          (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
           x124 * x473) *
            (-x416 + x424) +
          (10000000.0 * sc[2] * x42) * (x290 - x335 + (-x69) * (-x416 + x424));
  J[44] = sc[1] * x1287 - 15000000.0 * sc[2] * x341 +
          15000000.0 * sc[2] * x343 - 60000000.0 * sc[2] * x374 -
          10000000.0 * sc[2] * x382 + 80000000.0 * sc[2] * x430 -
          80000000.0 * sc[2] * x435 + sc[4] * x1746 - 500000.0 * x1157 +
          500000.0 * x1160 + 240000.0 * x1183 + 3612.0 * x12 +
          240000.0 * x1204 * x1250 - 240000.0 * x1224 - 3612.0 * x1239 + x1744 +
          x374 * x400 + x374 * x473 +
          (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
           x124 * x473) *
            (-x430 + x435) +
          (10000000.0 * sc[2] * x42) * (x341 - x343 + (-x69) * (-x430 + x435));
  J[58] = 10000000.0 * sc[11] * x238 + sc[1] * x860 -
          80000000.0 * sc[2] * x117 - sc[2] * x1369 -
          15000000.0 * sc[2] * x161 + 15000000.0 * sc[2] * x345 -
          60000000.0 * sc[2] * x375 - 10000000.0 * sc[2] * x390 +
          80000000.0 * sc[2] * x443 + sc[3] * x1746 + sc[7] * x84 +
          30000000.0 * sc[9] * x293 + x1240 + x1425 + 0.0714 * x1762 +
          x218 * x472 + x375 * x400 + x375 * x473 +
          (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
           x124 * x473) *
            (x117 - x443) +
          (10000000.0 * sc[2] * x42) * (x161 - x345 + (-x69) * (x117 - x443));
  J[72] = -sc[2] * x1367 - 80000000.0 * sc[2] * x170 - 0.0357 * sc[2] * x1736 -
          15000000.0 * sc[2] * x346 + 15000000.0 * sc[2] * x349 -
          10000000.0 * sc[2] * x395 + 80000000.0 * sc[2] * x466 -
          60000000.0 * sc[2] * x95 - 3000000.0 * x1157 + 3000000.0 * x1160 +
          3696000.0 * x1183 + 3612.0 * x12 - 3696000.0 * x1224 -
          3612.0 * x1239 + x400 * x95 + x473 * x95 +
          (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
           x124 * x473) *
            (x170 - x466) +
          (10000000.0 * sc[2] * x42) * (x346 - x349 + (-x69) * (x170 - x466));
  J[86] = sc[1] * x1283 + 15000000.0 * sc[2] * x106 -
          10000000.0 * sc[2] * x224 - 60000000.0 * sc[2] * x281 +
          80000000.0 * sc[2] * x478 - 15000000.0 * sc[2] * x486 -
          80000000.0 * sc[2] * x487 - 20000000.0 * sc[2] + x1240 + x281 * x400 +
          x281 * x473 +
          (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
           x124 * x473) *
            (-x478 + x487) +
          (10000000.0 * sc[2] * x42) * (-x106 + x486 + (-x69) * (-x478 + x487));
  J[100] = 15000000.0 * sc[2] * x106 - 80000000.0 * sc[2] * x320 -
           60000000.0 * sc[2] * x385 - 10000000.0 * sc[2] * x399 +
           80000000.0 * sc[2] * x490 - 15000000.0 * sc[2] * x495 -
           84300000.0 * sc[2] + sc[4] * x84 + x1240 + x385 * x400 +
           x385 * x473 +
           (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
            x124 * x473) *
             (x320 - x490) +
           (10000000.0 * sc[2] * x42) * (-x106 + x495 + (-x69) * (x320 - x490));
  J[114] = -10000000.0 * sc[2] * x202 - 15000000.0 * sc[2] * x352 +
           15000000.0 * sc[2] * x356 - 60000000.0 * sc[2] * x386 -
           sc[2] * x449 + 80000000.0 * sc[2] * x501 -
           80000000.0 * sc[2] * x508 + x1162 + 480000.0 * x1183 + 1204.0 * x12 -
           480000.0 * x1224 - 1204.0 * x1239 + x386 * x400 + x386 * x473 +
           (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
            x124 * x473) *
             (-x501 + x508) +
           (10000000.0 * sc[2] * x42) * (x352 - x356 + (-x69) * (-x501 + x508));
  J[128] = 15000000.0 * sc[0] * x274 - sc[2] * x1237 * x32 +
           80000000.0 * sc[2] * x148 - 80000000.0 * sc[2] * x164 -
           10000000.0 * sc[2] * x204 + 15000000.0 * sc[2] * x27 -
           15000000.0 * sc[2] * x359 - 60000000.0 * sc[2] * x392 +
           30000000.0 * sc[4] * x293 - 750000.0 * x1157 + 750000.0 * x1160 +
           420000.0 * x1183 + 903.0 * x12 - 420000.0 * x1224 - 903.0 * x1239 +
           x1785 + x392 * x400 + x392 * x473 +
           (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
            x124 * x473) *
             (-x148 + x164) +
           (10000000.0 * sc[2] * x42) * (-x27 + x359 + (-x69) * (-x148 + x164));
  J[142] = refC * x1233 * x32 + 30000000.0 * sc[1] * x298 -
           60000000.0 * sc[2] * x141 + 80000000.0 * sc[2] * x178 -
           80000000.0 * sc[2] * x179 - sc[2] * x18 - 10000000.0 * sc[2] * x202 -
           15000000.0 * sc[2] * x367 + 15000000.0 * sc[2] * x370 + x1162 +
           864000.0 * x1183 + 2107.0 * x12 - 864000.0 * x1224 - 2107.0 * x1239 +
           x141 * x400 + x141 * x473 +
           (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
            x124 * x473) *
             (-x178 + x179) +
           (10000000.0 * sc[2] * x42) * (x367 - x370 + (-x69) * (-x178 + x179));
  J[156] =
    sc[1] * x1310 + 15000000.0 * sc[2] * x106 - sc[2] * x1740 +
    80000000.0 * sc[2] * x186 - 15000000.0 * sc[2] * x196 -
    80000000.0 * sc[2] * x197 - 10000000.0 * sc[2] * x244 -
    60000000.0 * sc[2] * x354 + 10000000.0 * sc[4] * x238 + x1240 +
    x354 * x400 + x354 * x473 +
    (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
     x124 * x473) *
      (-x186 + x197) +
    (10000000.0 * sc[2] * x42) * (-x106 + x196 + (-x69) * (-x186 + x197));
  J[170] = -60000000.0 * sc[2] * x119 + 80000000.0 * sc[2] * x172 -
           80000000.0 * sc[2] * x190 - 10000000.0 * sc[2] * x34 -
           15000000.0 * sc[2] * x54 + 15000000.0 * sc[2] * x58 + x119 * x400 +
           x119 * x473 + x1240 +
           (60000000.0 * sc[2] * x124 + 15000000.0 * sc[2] * x69 - x124 * x400 -
            x124 * x473) *
             (-x172 + x190) +
           (10000000.0 * sc[2] * x42) * (x54 - x58 + (-x69) * (-x172 + x190));
  J[3] = -sc[3] * x192 - 40000000.0 * sc[3] * x215 + 40000000.0 * sc[3] * x220 -
         x109 * x315 + x1156 - 288000.0 * x1183 + 288000.0 * x1224 -
         x142 * x484 + x142 * x506 + x142 * x53 + x1811 * x342 - x1811 * x377 +
         (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
           (-x342 + x377) +
         (x109 * x42) * (x215 - x220 + (-x69) * (-x342 + x377));
  J[17] = 40000000.0 * sc[3] * x286 - 40000000.0 * sc[3] * x407 + sc[6] * x858 -
          x109 * x16 + x1156 + x1259 + x1301 - x150 * x1811 + x1811 * x393 +
          x1825 - x209 * x484 + x209 * x506 + x209 * x53 +
          (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
            (x150 - x393) +
          (x109 * x42) * (-x286 + x407 + (-x69) * (x150 - x393));
  J[31] = sc[10] * x18 + 240000.0 * sc[2] * x1250 * x865 -
          40000000.0 * sc[3] * x290 + 40000000.0 * sc[3] * x335 +
          sc[4] * x1369 + 20000000.0 * sc[6] - x109 * x30 + x1156 +
          x1811 * x416 - x1811 * x424 + x1825 - x277 * x484 + x277 * x506 +
          x277 * x53 +
          (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
            (-x416 + x424) +
          (x109 * x42) * (x290 - x335 + (-x69) * (-x416 + x424));
  J[45] = -sc[0] * x192 - sc[11] * x1860 - 40000000.0 * sc[3] * x341 +
          40000000.0 * sc[3] * x343 - sc[4] * x1746 - sc[5] * x1839 -
          sc[7] * x1854 - sc[8] * x70 - x109 * x382 - 120000.0 * x1204 * x1250 -
          18000000.0 * x1232 * exp(-452.89499922424801 * invT) + x1411 -
          13200000.0 * x168 * exp(-754.82499870708 * invT) - x1744 +
          x1811 * x430 - x1811 * x435 + x1825 -
          7600000.0 * x218 * exp(-201.28666632188799 * invT) - x374 * x484 +
          x374 * x506 + x374 * x53 - 40000000.0 * x51 +
          (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
            (-x430 + x435) +
          (x109 * x42) * (x341 - x343 + (-x69) * (-x430 + x435));
  J[59] = sc[11] * x1864 + sc[2] * x1369 - 40000000.0 * sc[3] * x161 -
          sc[3] * x1746 + 40000000.0 * sc[3] * x345 +
          28000000.0 * sc[9] * x257 * x258 - x109 * x390 + x1156 -
          x117 * x1811 + x1811 * x443 + x1825 + x1875 + x218 * x498 -
          x375 * x484 + x375 * x506 + x375 * x53 +
          (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
            (x117 - x443) +
          (x109 * x42) * (x161 - x345 + (-x69) * (x117 - x443));
  J[73] = -sc[3] * x1839 - 40000000.0 * sc[3] * x346 +
          40000000.0 * sc[3] * x349 + 12000000.0 * sc[9] * x254 - x109 * x395 +
          x1149 * x1476 - x1153 * x1290 - 1848000.0 * x1183 +
          1848000.0 * x1224 - x170 * x1811 + x1811 * x466 - x484 * x95 +
          x506 * x95 + x53 * x95 +
          (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
            (x170 - x466) +
          (x109 * x42) * (x346 - x349 + (-x69) * (x170 - x466));
  J[87] = sc[11] * x1890 + sc[1] * x858 + 20000000.0 * sc[2] +
          40000000.0 * sc[3] * x106 - 40000000.0 * sc[3] * x486 +
          1000000.0 * sc[7] + sc[9] * x87 - x109 * x224 + x1156 + x1521 +
          x1811 * x478 - x1811 * x487 + x1825 + x1889 + x218 * x476 -
          x281 * x484 + x281 * x506 + x281 * x53 +
          (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
            (-x478 + x487) +
          (x109 * x42) * (-x106 + x486 + (-x69) * (-x478 + x487));
  J[101] = 40000000.0 * sc[3] * x106 - sc[3] * x1854 -
           40000000.0 * sc[3] * x495 + 1000000.0 * sc[6] - x109 * x399 + x1156 -
           x1811 * x320 + x1811 * x490 + x1825 - x385 * x484 + x385 * x506 +
           x385 * x53 +
           (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
             (x320 - x490) +
           (x109 * x42) * (-x106 + x495 + (-x69) * (x320 - x490));
  J[115] = -40000000.0 * sc[3] * x352 + 40000000.0 * sc[3] * x356 -
           sc[3] * x70 - x109 * x202 + x1156 - 240000.0 * x1183 +
           240000.0 * x1224 + x1811 * x501 - x1811 * x508 - x386 * x484 +
           x386 * x506 + x386 * x53 +
           (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
             (-x501 + x508) +
           (x109 * x42) * (x352 - x356 + (-x69) * (-x501 + x508));
  J[129] = 40000000.0 * sc[3] * x27 - 40000000.0 * sc[3] * x359 +
           28000000.0 * sc[4] * x257 * x258 + 12000000.0 * sc[5] * x254 +
           sc[6] * x87 - x109 * x204 + 2100000.0 * x1148 * x1149 -
           2100000.0 * x1152 * x1153 - 210000.0 * x1183 + 210000.0 * x1224 +
           x148 * x1811 - x164 * x1811 - x1785 - x392 * x484 + x392 * x506 +
           x392 * x53 +
           (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
             (-x148 + x164) +
           (x109 * x42) * (-x27 + x359 + (-x69) * (-x148 + x164));
  J[143] = sc[2] * x18 - 40000000.0 * sc[3] * x367 + 40000000.0 * sc[3] * x370 -
           x109 * x202 + 4200000.0 * x1148 * x1149 - 4200000.0 * x1152 * x1153 -
           432000.0 * x1183 + 432000.0 * x1224 - x141 * x484 + x141 * x506 +
           x141 * x53 + x178 * x1811 - x179 * x1811 +
           (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
             (-x178 + x179) +
           (x109 * x42) * (x367 - x370 + (-x69) * (-x178 + x179));
  J[157] = 40000000.0 * sc[3] * x106 - sc[3] * x1860 -
           40000000.0 * sc[3] * x196 + sc[4] * x1864 + sc[6] * x1890 -
           x109 * x244 + x1156 + x1811 * x186 - x1811 * x197 + x1825 -
           x354 * x484 + x354 * x506 + x354 * x53 +
           (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
             (-x186 + x197) +
           (x109 * x42) * (-x106 + x196 + (-x69) * (-x186 + x197));
  J[171] = -40000000.0 * sc[3] * x54 + 40000000.0 * sc[3] * x58 - x109 * x34 +
           375000000.0 * x1149 * x1406 - 375000000.0 * x1153 * x1293 -
           x119 * x484 + x119 * x506 + x119 * x53 + x172 * x1811 -
           x1811 * x190 + x1825 +
           (40000000.0 * sc[3] * x69 + x124 * x484 - x124 * x506 - x124 * x53) *
             (-x172 + x190) +
           (x109 * x42) * (x54 - x58 + (-x69) * (-x172 + x190));
  J[4] = 30000000.0 * sc[2] * x142 + 10000000.0 * sc[2] * x315 -
         50000000.0 * sc[4] * x142 - 30000000.0 * sc[4] * x215 +
         30000000.0 * sc[4] * x220 - 5000000.0 * sc[4] * x315 +
         20000000.0 * sc[4] * x342 - 20000000.0 * sc[4] * x377 -
         20000000.0 * sc[6] * x342 + 20000000.0 * sc[6] * x377 +
         1000000.0 * x1157 - 1000000.0 * x1160 + x1174 + x128 * x142 -
         x142 * x473 - x142 * x506 - x1811 * x342 + x1811 * x377 + x1978 -
         x200 +
         (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
          50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
          x124 * x473 + x124 * x506 - x133 * x69) *
           (-x342 + x377) +
         (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
          5000000.0 * sc[4] * x42) *
           (x215 - x220 + (-x69) * (-x342 + x377));
  J[18] = sc[10] * x1276 + 30000000.0 * sc[11] * x267 +
          20000000.0 * sc[11] * x328 + 10000000.0 * sc[2] * x16 +
          30000000.0 * sc[2] * x209 + 500000.0 * sc[2] * x8 * x865 +
          sc[3] * x1287 - 20000000.0 * sc[4] * x150 - 5000000.0 * sc[4] * x16 -
          50000000.0 * sc[4] * x209 + 30000000.0 * sc[4] * x286 +
          20000000.0 * sc[4] * x393 - 30000000.0 * sc[4] * x407 + sc[5] * x882 +
          20000000.0 * sc[6] * x150 - 20000000.0 * sc[6] * x393 +
          268000000.0 * sc[6] * exp(-319.54258278599701 * invT) - x1259 +
          x128 * x209 + x1308 + x150 * x1811 - x1811 * x393 + x1982 + x1988 -
          x209 * x473 - x209 * x506 +
          (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
           50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
           x124 * x473 + x124 * x506 - x133 * x69) *
            (x150 - x393) +
          (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
           5000000.0 * sc[4] * x42) *
            (-x286 + x407 + (-x69) * (x150 - x393));
  J[32] = sc[0] * x21 + sc[11] * x1740 + 500000.0 * sc[1] * x8 * x865 +
          28000000.0 * sc[3] * x290 - 28000000.0 * sc[3] * x335 -
          sc[4] * x1369 - 50000000.0 * sc[4] * x277 -
          30000000.0 * sc[4] * x290 - 5000000.0 * sc[4] * x30 +
          30000000.0 * sc[4] * x335 + 20000000.0 * sc[4] * x416 -
          20000000.0 * sc[4] * x424 + 0.0714 * sc[5] * x1736 +
          20000000.0 * sc[6] + x1382 - x1811 * x416 + x1811 * x424 + x1982 +
          x2000 + x2004 - x277 * x473 - x277 * x506 +
          (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
           50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
           x124 * x473 + x124 * x506 - x133 * x69) *
            (-x416 + x424) +
          (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
           5000000.0 * sc[4] * x42) *
            (x290 - x335 + (-x69) * (-x416 + x424));
  J[46] = sc[1] * x1287 + 30000000.0 * sc[2] * x374 +
          28000000.0 * sc[3] * x341 - 28000000.0 * sc[3] * x343 -
          sc[4] * x1746 - 30000000.0 * sc[4] * x341 +
          30000000.0 * sc[4] * x343 - 50000000.0 * sc[4] * x374 -
          5000000.0 * sc[4] * x382 + 20000000.0 * sc[4] * x430 -
          20000000.0 * sc[4] * x435 + sc[5] * x1839 +
          13200000.0 * x168 * exp(-754.82499870708 * invT) - x1811 * x430 +
          x1811 * x435 + x1982 + x2011 + x2015 - x374 * x473 - x374 * x506 +
          28000000.0 * x51 +
          (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
           50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
           x124 * x473 + x124 * x506 - x133 * x69) *
            (-x430 + x435) +
          (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
           5000000.0 * sc[4] * x42) *
            (x341 - x343 + (-x69) * (-x430 + x435));
  J[60] = -sc[10] * x2029 - sc[2] * x1369 + 30000000.0 * sc[2] * x375 -
          sc[3] * x1746 - 20000000.0 * sc[4] * x117 -
          30000000.0 * sc[4] * x161 + 30000000.0 * sc[4] * x345 -
          50000000.0 * sc[4] * x375 - 5000000.0 * sc[4] * x390 +
          20000000.0 * sc[4] * x443 - 28000000.0 * sc[9] * x257 * x258 -
          30000000.0 * sc[9] * x293 + x117 * x1811 - x1420 -
          536000000.0 * x1422 - x1425 + x1449 - 20000000.0 * x168 +
          50000000.0 * x171 - 0.1428 * x1762 - x1811 * x443 - x1875 + x1982 -
          x2019 + x2040 + x2047 - x218 * x472 - x218 * x498 + 5000000.0 * x245 +
          50000000.0 * x316 - x375 * x473 - x375 * x506 - 30000000.0 * x51 +
          5000000.0 * x909 - x940 +
          (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
           50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
           x124 * x473 + x124 * x506 - x133 * x69) *
            (x117 - x443) +
          (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
           5000000.0 * sc[4] * x42) *
            (x161 - x345 + (-x69) * (x117 - x443));
  J[74] = 5000000.0 * sc[11] * x226 + sc[1] * x882 + 0.0714 * sc[2] * x1736 +
          30000000.0 * sc[2] * x95 + sc[3] * x1839 - 20000000.0 * sc[4] * x170 -
          30000000.0 * sc[4] * x346 + 30000000.0 * sc[4] * x349 -
          5000000.0 * sc[4] * x395 + 20000000.0 * sc[4] * x466 -
          50000000.0 * sc[4] * x95 + sc[7] * x45 + 50000000.0 * sc[9] * x302 +
          3000000.0 * x1157 - 3000000.0 * x1160 + x1499 + x170 * x1811 -
          x1811 * x466 + x2062 + x2067 - x473 * x95 - x506 * x95 +
          (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
           50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
           x124 * x473 + x124 * x506 - x133 * x69) *
            (x170 - x466) +
          (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
           5000000.0 * sc[4] * x42) *
            (x346 - x349 + (-x69) * (x170 - x466));
  J[88] = 268000000.0 * sc[1] * exp(-319.54258278599701 * invT) +
          30000000.0 * sc[2] * x281 + 20000000.0 * sc[2] +
          28000000.0 * sc[3] * x486 + 30000000.0 * sc[4] * x106 -
          5000000.0 * sc[4] * x224 - 50000000.0 * sc[4] * x281 +
          20000000.0 * sc[4] * x478 - 30000000.0 * sc[4] * x486 -
          20000000.0 * sc[4] * x487 + sc[9] * x2071 + x1525 - x1811 * x478 +
          x1811 * x487 - x1889 + x1982 + x2079 + x2084 - x281 * x473 -
          x281 * x506 +
          (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
           50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
           x124 * x473 + x124 * x506 - x133 * x69) *
            (-x478 + x487) +
          (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
           5000000.0 * sc[4] * x42) *
            (-x106 + x486 + (-x69) * (-x478 + x487));
  J[102] = 30000000.0 * sc[2] * x385 + 28000000.0 * sc[3] * x495 +
           30000000.0 * sc[4] * x106 - 20000000.0 * sc[4] * x320 -
           50000000.0 * sc[4] * x385 - 5000000.0 * sc[4] * x399 +
           20000000.0 * sc[4] * x490 - 30000000.0 * sc[4] * x495 + sc[5] * x45 -
           x106 * x133 + x1525 + x1811 * x320 - x1811 * x490 + x1982 + x2092 +
           x2098 - x385 * x473 - x385 * x506 +
           (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
            50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
            x124 * x473 + x124 * x506 - x133 * x69) *
             (x320 - x490) +
           (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
            5000000.0 * sc[4] * x42) *
             (-x106 + x495 + (-x69) * (x320 - x490));
  J[116] = 12000000.0 * sc[1] * x202 + 10000000.0 * sc[2] * x202 +
           30000000.0 * sc[2] * x386 + 28000000.0 * sc[3] * x352 -
           28000000.0 * sc[3] * x356 - 5000000.0 * sc[4] * x202 -
           30000000.0 * sc[4] * x352 + 30000000.0 * sc[4] * x356 -
           50000000.0 * sc[4] * x386 + 20000000.0 * sc[4] * x501 -
           20000000.0 * sc[4] * x508 - 20000000.0 * sc[6] * x501 +
           20000000.0 * sc[6] * x508 + 1000000.0 * x1157 - 1000000.0 * x1160 -
           44000000000.0 * x1169 + 44000000000.0 * x1172 + x128 * x386 -
           x1811 * x501 + x1811 * x508 - x2099 + x2113 - x386 * x473 -
           x386 * x506 +
           (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
            50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
            x124 * x473 + x124 * x506 - x133 * x69) *
             (-x501 + x508) +
           (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
            5000000.0 * sc[4] * x42) *
             (x352 - x356 + (-x69) * (-x501 + x508));
  J[130] = 10000000.0 * sc[2] * x204 + 30000000.0 * sc[2] * x392 -
           28000000.0 * sc[3] * x27 + 28000000.0 * sc[3] * x359 +
           20000000.0 * sc[4] * x148 - 20000000.0 * sc[4] * x164 -
           5000000.0 * sc[4] * x204 - 28000000.0 * sc[4] * x257 * x258 +
           30000000.0 * sc[4] * x27 - 30000000.0 * sc[4] * x293 -
           30000000.0 * sc[4] * x359 - 50000000.0 * sc[4] * x392 +
           50000000.0 * sc[5] * x302 - 20000000.0 * sc[6] * x148 +
           20000000.0 * sc[6] * x164 + sc[6] * x2071 + 750000.0 * x1157 -
           750000.0 * x1160 + x128 * x392 + x1360 - x148 * x1811 - x1579 +
           x164 * x1811 + x2123 - x392 * x473 - x392 * x506 +
           (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
            50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
            x124 * x473 + x124 * x506 - x133 * x69) *
             (-x148 + x164) +
           (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
            5000000.0 * sc[4] * x42) *
             (-x27 + x359 + (-x69) * (-x148 + x164));
  J[144] =
    sc[1] * x1276 + 12000000.0 * sc[1] * x202 + 30000000.0 * sc[2] * x141 +
    10000000.0 * sc[2] * x202 + 28000000.0 * sc[3] * x367 -
    28000000.0 * sc[3] * x370 - 50000000.0 * sc[4] * x141 +
    20000000.0 * sc[4] * x178 - 20000000.0 * sc[4] * x179 -
    5000000.0 * sc[4] * x202 - sc[4] * x2029 - 30000000.0 * sc[4] * x367 +
    30000000.0 * sc[4] * x370 - 20000000.0 * sc[6] * x178 +
    20000000.0 * sc[6] * x179 + 1000000.0 * x1157 - 1000000.0 * x1160 +
    x128 * x141 + x133 * x367 - x133 * x370 + x1360 - x141 * x473 -
    x141 * x506 - x178 * x1811 - x178 * x1973 + x179 * x1811 + x179 * x1973 +
    (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
     50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
     x124 * x473 + x124 * x506 - x133 * x69) *
      (-x178 + x179) +
    (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
     5000000.0 * sc[4] * x42) *
      (x367 - x370 + (-x69) * (-x178 + x179));
  J[158] = 12000000.0 * sc[1] * x244 + 30000000.0 * sc[1] * x267 +
           20000000.0 * sc[1] * x328 + sc[2] * x1740 +
           30000000.0 * sc[2] * x354 + 28000000.0 * sc[3] * x196 +
           30000000.0 * sc[4] * x106 + 20000000.0 * sc[4] * x186 -
           30000000.0 * sc[4] * x196 - 20000000.0 * sc[4] * x197 -
           5000000.0 * sc[4] * x244 - 50000000.0 * sc[4] * x354 +
           5000000.0 * sc[5] * x226 - x106 * x133 + x133 * x196 + x1525 -
           x1811 * x186 + x1811 * x197 - x186 * x1973 + x197 * x1973 + x1982 +
           x2153 - x354 * x473 - x354 * x506 +
           (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
            50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
            x124 * x473 + x124 * x506 - x133 * x69) *
             (-x186 + x197) +
           (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
            5000000.0 * sc[4] * x42) *
             (-x106 + x196 + (-x69) * (-x186 + x197));
  J[172] = 30000000.0 * sc[2] * x119 + 10000000.0 * sc[2] * x34 +
           28000000.0 * sc[3] * x54 - 28000000.0 * sc[3] * x58 -
           50000000.0 * sc[4] * x119 + 20000000.0 * sc[4] * x172 -
           20000000.0 * sc[4] * x190 - 5000000.0 * sc[4] * x34 -
           30000000.0 * sc[4] * x54 + 30000000.0 * sc[4] * x58 -
           20000000.0 * sc[6] * x172 + 20000000.0 * sc[6] * x190 + x119 * x128 -
           x119 * x473 - x119 * x506 - x172 * x1811 + x1811 * x190 + x1982 +
           x2164 +
           (-30000000.0 * sc[2] * x124 - 28000000.0 * sc[3] * x69 +
            50000000.0 * sc[4] * x124 + 30000000.0 * sc[4] * x69 - x124 * x128 +
            x124 * x473 + x124 * x506 - x133 * x69) *
             (-x172 + x190) +
           (-12000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 +
            5000000.0 * sc[4] * x42) *
             (x54 - x58 + (-x69) * (-x172 + x190));
  J[5] = 6000000.0 * sc[1] * x315 + 12000000.0 * sc[3] * x215 -
         12000000.0 * sc[3] * x220 + 50000000.0 * sc[4] * x142 +
         5000000.0 * sc[4] * x315 + 16060000000.0 * x1169 -
         16060000000.0 * x1172 - x128 * x142 - x133 * x215 + x133 * x220 +
         x1973 * x342 - x1973 * x377 + x200 - x215 * x428 + x220 * x428 +
         (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 + x124 * x128 +
          x133 * x69 + x428 * x69) *
           (-x342 + x377) +
         (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
           (x215 - x220 + (-x69) * (-x342 + x377));
  J[19] = 6000000.0 * sc[1] * x16 - 12000000.0 * sc[3] * x286 +
          12000000.0 * sc[3] * x407 + sc[4] * x1306 * x884 +
          5000000.0 * sc[4] * x16 + 50000000.0 * sc[4] * x209 - sc[5] * x882 +
          sc[6] * x1283 - x128 * x209 + x133 * x286 - x133 * x407 -
          x150 * x1973 + x1973 * x393 + x2165 - 6000000.0 * x245 + x286 * x428 -
          x407 * x428 - 6000000.0 * x909 +
          (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 + x124 * x128 +
           x133 * x69 + x428 * x69) *
            (x150 - x393) +
          (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
            (-x286 + x407 + (-x69) * (x150 - x393));
  J[33] = 6000000.0 * sc[1] * x30 + 12000000.0 * sc[3] * x290 -
          12000000.0 * sc[3] * x335 + 50000000.0 * sc[4] * x277 +
          5000000.0 * sc[4] * x30 - sc[5] * x1367 - 0.0357 * sc[5] * x1736 -
          x128 * x277 - x133 * x290 + x133 * x335 + x1973 * x416 -
          x1973 * x424 + x2165 - x290 * x428 + x335 * x428 +
          (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 + x124 * x128 +
           x133 * x69 + x428 * x69) *
            (-x416 + x424) +
          (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
            (x290 - x335 + (-x69) * (-x416 + x424));
  J[47] = 6000000.0 * sc[1] * x382 + 12000000.0 * sc[3] * x341 -
          12000000.0 * sc[3] * x343 + 50000000.0 * sc[4] * x374 +
          5000000.0 * sc[4] * x382 - sc[5] * x1839 - x128 * x374 - x133 * x341 +
          x133 * x343 + x1973 * x430 - x1973 * x435 + x2165 - x341 * x428 +
          x343 * x428 + 12000000.0 * x51 +
          (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 + x124 * x128 +
           x133 * x69 + x428 * x69) *
            (-x430 + x435) +
          (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
            (x341 - x343 + (-x69) * (-x430 + x435));
  J[61] = sc[11] * x2044 + sc[1] * x1306 * x884 + 6000000.0 * sc[1] * x390 +
          12000000.0 * sc[3] * x161 - 12000000.0 * sc[3] * x345 +
          50000000.0 * sc[4] * x375 + 5000000.0 * sc[4] * x390 + sc[7] * x2034 +
          25010000.0 * sc[7] - x117 * x1973 - x128 * x375 - x133 * x161 +
          x133 * x345 - x161 * x428 - 50000000.0 * x171 + 0.0714 * x1762 +
          x1875 + x1973 * x443 + x2019 + x2165 - 5000000.0 * x245 -
          50000000.0 * x316 + x345 * x428 - 5000000.0 * x909 + x940 +
          (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 + x124 * x128 +
           x133 * x69 + x428 * x69) *
            (x117 - x443) +
          (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
            (x161 - x345 + (-x69) * (x117 - x443));
  J[75] = -refC * x1171 * x1306 - 5000000.0 * sc[11] * x226 +
          6000000.0 * sc[1] * x395 - sc[1] * x882 - sc[2] * x1367 -
          0.0357 * sc[2] * x1736 - sc[3] * x1839 + 12000000.0 * sc[3] * x346 -
          12000000.0 * sc[3] * x349 + 5000000.0 * sc[4] * x395 +
          50000000.0 * sc[4] * x95 - sc[7] * x45 - 12000000.0 * sc[9] * x254 -
          50000000.0 * sc[9] * x302 + 80300000000.0 * x1169 -
          80300000000.0 * x1172 - x128 * x95 - x133 * x346 + x133 * x349 -
          x139 * x218 - x168 * x279 - x170 * x1973 + x1973 * x466 -
          x346 * x428 + x349 * x428 - 6000000.0 * x38 * x51 -
          25010000.0 * x47 * x51 +
          (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 + x124 * x128 +
           x133 * x69 + x428 * x69) *
            (x170 - x466) +
          (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
            (x346 - x349 + (-x69) * (x170 - x466));
  J[89] = sc[1] * x1283 + 6000000.0 * sc[1] * x224 + 12000000.0 * sc[3] * x486 +
          5000000.0 * sc[4] * x224 + 50000000.0 * sc[4] * x281 - x128 * x281 -
          x133 * x486 + x1889 + x1973 * x478 - x1973 * x487 + x2171 -
          x428 * x486 +
          (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 + x124 * x128 +
           x133 * x69 + x428 * x69) *
            (-x478 + x487) +
          (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
            (-x106 + x486 + (-x69) * (-x478 + x487));
  J[103] = 6000000.0 * sc[1] * x399 + 12000000.0 * sc[3] * x495 +
           sc[4] * x2034 + 50000000.0 * sc[4] * x385 +
           5000000.0 * sc[4] * x399 + 25010000.0 * sc[4] - sc[5] * x45 -
           x128 * x385 - x133 * x495 - x1973 * x320 + x1973 * x490 + x2171 -
           x428 * x495 +
           (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 +
            x124 * x128 + x133 * x69 + x428 * x69) *
             (x320 - x490) +
           (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
             (-x106 + x495 + (-x69) * (x320 - x490));
  J[117] = 6000000.0 * sc[1] * x202 + 12000000.0 * sc[3] * x352 -
           12000000.0 * sc[3] * x356 + 5000000.0 * sc[4] * x202 +
           50000000.0 * sc[4] * x386 + 44000000000.0 * x1169 -
           44000000000.0 * x1172 - x128 * x386 - x133 * x352 + x133 * x356 +
           x1973 * x501 - x1973 * x508 + x2099 - x352 * x428 + x356 * x428 +
           (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 +
            x124 * x128 + x133 * x69 + x428 * x69) *
             (-x501 + x508) +
           (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
             (x352 - x356 + (-x69) * (-x501 + x508));
  J[131] = 6000000.0 * sc[1] * x204 - 12000000.0 * sc[3] * x27 +
           12000000.0 * sc[3] * x359 + 5000000.0 * sc[4] * x204 +
           50000000.0 * sc[4] * x392 - 12000000.0 * sc[5] * x254 -
           50000000.0 * sc[5] * x302 - x128 * x392 + x133 * x27 - x133 * x359 +
           x148 * x1973 - x164 * x1973 + x2165 + x27 * x428 - x359 * x428 +
           (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 +
            x124 * x128 + x133 * x69 + x428 * x69) *
             (-x148 + x164) +
           (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
             (-x27 + x359 + (-x69) * (-x148 + x164));
  J[145] = 6000000.0 * sc[1] * x202 + 12000000.0 * sc[3] * x367 -
           12000000.0 * sc[3] * x370 + 50000000.0 * sc[4] * x141 +
           5000000.0 * sc[4] * x202 - x128 * x141 - x133 * x367 + x133 * x370 +
           x178 * x1973 - x179 * x1973 + x2165 - x367 * x428 + x370 * x428 +
           (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 +
            x124 * x128 + x133 * x69 + x428 * x69) *
             (-x178 + x179) +
           (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
             (x367 - x370 + (-x69) * (-x178 + x179));
  J[159] = 6000000.0 * sc[1] * x244 + 12000000.0 * sc[3] * x196 +
           sc[4] * x2044 + 5000000.0 * sc[4] * x244 +
           50000000.0 * sc[4] * x354 - 5000000.0 * sc[5] * x226 - x128 * x354 -
           x133 * x196 + x186 * x1973 - x196 * x428 - x197 * x1973 + x2171 +
           (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 +
            x124 * x128 + x133 * x69 + x428 * x69) *
             (-x186 + x197) +
           (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
             (-x106 + x196 + (-x69) * (-x186 + x197));
  J[173] = 6000000.0 * sc[1] * x34 + 12000000.0 * sc[3] * x54 -
           12000000.0 * sc[3] * x58 + 50000000.0 * sc[4] * x119 +
           5000000.0 * sc[4] * x34 - x119 * x128 - x133 * x54 + x133 * x58 +
           x172 * x1973 - x190 * x1973 + x2165 - x428 * x54 + x428 * x58 +
           (-12000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 +
            x124 * x128 + x133 * x69 + x428 * x69) *
             (-x172 + x190) +
           (-6000000.0 * sc[1] * x42 - 5000000.0 * sc[4] * x42) *
             (x54 - x58 + (-x69) * (-x172 + x190));
  J[6] = sc[3] * x192 + 20000000.0 * sc[6] * x342 - 20000000.0 * sc[6] * x377 +
         x109 * x315 + x142 * x484 - x142 * x53 + x2176 +
         (-x124 * x484 + x124 * x53) * (-x342 + x377) +
         (-x109 * x42) * (x215 - x220 + (-x69) * (-x342 + x377));
  J[20] = 375000000.0 * sc[12] * sc[3] * x1293 +
          300000000.0 * ((sc[3]) * (sc[3])) * x1293 + sc[3] * sc[5] * x1290 +
          sc[3] * x1152 * x1297 - 20000000.0 * sc[6] * x150 +
          20000000.0 * sc[6] * x393 + x109 * x16 + x1285 + x209 * x484 -
          x209 * x53 + x2176 + (-x124 * x484 + x124 * x53) * (x150 - x393) +
          (-x109 * x42) * (-x286 + x407 + (-x69) * (x150 - x393));
  J[34] = sc[5] * x1367 + 20000000.0 * sc[6] * x416 -
          20000000.0 * sc[6] * x424 - 20000000.0 * sc[6] + x109 * x30 + x2176 +
          x277 * x484 - x277 * x53 +
          (-x124 * x484 + x124 * x53) * (-x416 + x424) +
          (-x109 * x42) * (x290 - x335 + (-x69) * (-x416 + x424));
  J[48] = sc[0] * x192 + sc[11] * x1860 + 375000000.0 * sc[12] * sc[1] * x1293 +
          sc[1] * sc[5] * x1290 + sc[1] * x1152 * x1297 + sc[4] * x1746 +
          sc[5] * x1839 + 20000000.0 * sc[6] * x430 -
          20000000.0 * sc[6] * x435 + sc[8] * x70 + x109 * x382 -
          300000000.0 * x1149 * x1406 + 600000000.0 * x1153 * x1293 +
          18000000.0 * x1232 * exp(-452.89499922424801 * invT) +
          7600000.0 * x218 * exp(-201.28666632188799 * invT) + x374 * x484 -
          x374 * x53 + (-x124 * x484 + x124 * x53) * (-x430 + x435) +
          (-x109 * x42) * (x341 - x343 + (-x69) * (-x430 + x435));
  J[62] = sc[10] * x2029 + 20000000.0 * sc[11] * x324 + sc[3] * x1746 -
          20000000.0 * sc[6] * x117 + 20000000.0 * sc[6] * x443 + x109 * x390 +
          268000000.0 * x1422 - x1875 + x2176 + x375 * x484 - x375 * x53 +
          (-x124 * x484 + x124 * x53) * (x117 - x443) +
          (-x109 * x42) * (x161 - x345 + (-x69) * (x117 - x443));
  J[76] = sc[2] * x1367 + sc[3] * x1839 - 20000000.0 * sc[6] * x170 +
          20000000.0 * sc[6] * x466 + x109 * x395 - x1149 * x1476 +
          x1153 * x1290 + x484 * x95 - x53 * x95 +
          (-x124 * x484 + x124 * x53) * (x170 - x466) +
          (-x109 * x42) * (x346 - x349 + (-x69) * (x170 - x466));
  J[90] =
    -375000000.0 * refC * sc[12] * x1406 - 300000000.0 * refC * sc[3] * x1406 -
    refC * sc[5] * x1476 - refC * x1148 * x1297 - sc[11] * x1890 -
    20000000.0 * sc[2] + 20000000.0 * sc[6] * x478 - 20000000.0 * sc[6] * x487 -
    1000000.0 * sc[7] - sc[9] * x2071 - sc[9] * x87 + x109 * x224 + x1515 -
    20000000.0 * x168 - x1889 + x2176 - x218 * x476 + x281 * x484 - x281 * x53 +
    (-x124 * x484 + x124 * x53) * (-x478 + x487) +
    (-x109 * x42) * (-x106 + x486 + (-x69) * (-x478 + x487));
  J[104] = -20000000.0 * sc[6] * x320 + 20000000.0 * sc[6] * x490 -
           1000000.0 * sc[6] + x109 * x399 + x2176 + x385 * x484 - x385 * x53 +
           (-x124 * x484 + x124 * x53) * (x320 - x490) +
           (-x109 * x42) * (-x106 + x495 + (-x69) * (x320 - x490));
  J[118] = sc[3] * x70 + 20000000.0 * sc[6] * x501 - 20000000.0 * sc[6] * x508 +
           x109 * x202 + x2176 + x386 * x484 - x386 * x53 +
           (-x124 * x484 + x124 * x53) * (-x501 + x508) +
           (-x109 * x42) * (x352 - x356 + (-x69) * (-x501 + x508));
  J[132] = 20000000.0 * sc[6] * x148 - 20000000.0 * sc[6] * x164 -
           sc[6] * x2071 - sc[6] * x87 + x109 * x204 -
           2100000.0 * x1148 * x1149 + 2100000.0 * x1152 * x1153 + x392 * x484 -
           x392 * x53 + (-x124 * x484 + x124 * x53) * (-x148 + x164) +
           (-x109 * x42) * (-x27 + x359 + (-x69) * (-x148 + x164));
  J[146] = sc[4] * x2029 + 20000000.0 * sc[6] * x178 -
           20000000.0 * sc[6] * x179 + x109 * x202 - 4200000.0 * x1148 * x1149 +
           4200000.0 * x1152 * x1153 + x141 * x484 - x141 * x53 +
           (-x124 * x484 + x124 * x53) * (-x178 + x179) +
           (-x109 * x42) * (x367 - x370 + (-x69) * (-x178 + x179));
  J[160] = sc[3] * x1860 + 20000000.0 * sc[4] * x324 +
           20000000.0 * sc[6] * x186 - sc[6] * x1890 -
           20000000.0 * sc[6] * x197 + x109 * x244 + x2176 + x354 * x484 -
           x354 * x53 + (-x124 * x484 + x124 * x53) * (-x186 + x197) +
           (-x109 * x42) * (-x106 + x196 + (-x69) * (-x186 + x197));
  J[174] = 20000000.0 * sc[6] * x172 - 20000000.0 * sc[6] * x190 + x109 * x34 -
           375000000.0 * x1149 * x1406 + 375000000.0 * x1153 * x1293 +
           x119 * x484 - x119 * x53 +
           (-x124 * x484 + x124 * x53) * (-x172 + x190) +
           (-x109 * x42) * (x54 - x58 + (-x69) * (-x172 + x190));
  J[7] = -26480000.0 * sc[7] * x142 - sc[7] * x338 + x101 * x142 -
         x1140 * x342 + x1140 * x377 + x1225 + x174 * x215 - x174 * x220 +
         x1978 - x2190 * x342 + x2190 * x377 + x2193 +
         (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 - x101 * x124 -
          x133 * x69 - x174 * x69) *
           (-x342 + x377) +
         (-12000000.0 * sc[1] * x42) * (x215 - x220 + (-x69) * (-x342 + x377));
  J[21] = sc[11] * x1310 - 26480000.0 * sc[7] * x209 + x101 * x209 +
          3200000000000000.0 * x1071 * x1074 * x8 + x1140 * x150 -
          x1140 * x393 + x1305 + x1327 + x150 * x2190 - x174 * x286 +
          x174 * x407 + x1988 - x2190 * x393 + x2196 + x863 +
          (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 - x101 * x124 -
           x133 * x69 - x174 * x69) *
            (x150 - x393) +
          (-12000000.0 * sc[1] * x42) * (-x286 + x407 + (-x69) * (x150 - x393));
  J[35] = 70000000.0 * sc[0] * x290 - 70000000.0 * sc[0] * x335 -
          26480000.0 * sc[7] * x277 - 84300000.0 * sc[7] + x101 * x277 -
          x1140 * x416 + x1140 * x424 + x1327 + x174 * x290 - x174 * x335 +
          x2000 - x205 * x416 + x205 * x424 - x2190 * x416 + x2190 * x424 +
          x2196 +
          (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 - x101 * x124 -
           x133 * x69 - x174 * x69) *
            (-x416 + x424) +
          (-12000000.0 * sc[1] * x42) * (x290 - x335 + (-x69) * (-x416 + x424));
  J[49] = 70000000.0 * sc[0] * x341 - 70000000.0 * sc[0] * x343 -
          sc[7] * x1854 - 26480000.0 * sc[7] * x374 + sc[8] * x70 +
          x101 * x374 - x1140 * x430 + x1140 * x435 + x1327 + x174 * x341 -
          x174 * x343 + x2011 - x205 * x430 + x205 * x435 - x2190 * x430 +
          x2190 * x435 + x2196 +
          (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 - x101 * x124 -
           x133 * x69 - x174 * x69) *
            (-x430 + x435) +
          (-12000000.0 * sc[1] * x42) * (x341 - x343 + (-x69) * (-x430 + x435));
  J[63] = 70000000.0 * sc[0] * x161 - 70000000.0 * sc[0] * x345 +
          sc[11] * x1864 - 26480000.0 * sc[7] * x375 + x101 * x375 +
          x1140 * x117 - x1140 * x443 + x117 * x205 + x117 * x2190 + x1327 +
          x161 * x174 - x174 * x345 + x2019 + x2040 - x205 * x443 -
          x2190 * x443 + x2196 +
          (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 - x101 * x124 -
           x133 * x69 - x174 * x69) *
            (x117 - x443) +
          (-12000000.0 * sc[1] * x42) * (x161 - x345 + (-x69) * (x117 - x443));
  J[77] = 1.92e+16 * sc[1] * x1071 * x1074 - sc[7] * x45 -
          26480000.0 * sc[7] * x95 + x101 * x95 - 2.4576e+21 * x1063 * x1105 +
          2.4576e+21 * x1063 * x331 + x1064 * x153 - x1064 * x93 -
          1.92e+16 * x1088 + x1140 * x170 - x1140 * x466 + x1495 +
          x170 * x2190 + x174 * x346 - x174 * x349 + x2062 - x2190 * x466 +
          (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 - x101 * x124 -
           x133 * x69 - x174 * x69) *
            (x170 - x466) +
          (-12000000.0 * sc[1] * x42) * (x346 - x349 + (-x69) * (x170 - x466));
  J[91] =
    70000000.0 * sc[0] * x486 - 26480000.0 * sc[7] * x281 - 1000000.0 * sc[7] +
    x101 * x281 - x1140 * x478 + x1140 * x487 + x174 * x486 - x205 * x478 +
    x205 * x487 + x2079 - x2190 * x478 + x2190 * x487 + x2208 +
    (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 - x101 * x124 -
     x133 * x69 - x174 * x69) *
      (-x478 + x487) +
    (-12000000.0 * sc[1] * x42) * (-x106 + x486 + (-x69) * (-x478 + x487));
  J[105] =
    -3200000000000000.0 * refC * x1071 * x1086 * x8 - sc[0] * x338 -
    sc[11] * x2215 - 84300000.0 * sc[2] - sc[3] * x1854 - sc[5] * x45 -
    1000000.0 * sc[6] - 26480000.0 * sc[7] * x385 -
    9.84 * sc[7] * exp(-g_RT_qss[0] + x2187 - x271) -
    64000000.0 * sc[7] * exp(-g_RT_qss[1] + 286.83349950869001 * invT - x271) +
    x101 * x385 - x106 * x133 + x1140 * x320 - x1140 * x490 + x1555 +
    26480000.0 * x171 + x174 * x495 + x2098 + x2190 * x320 - x2190 * x490 +
    x2208 + 26480000.0 * x316 +
    (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 - x101 * x124 -
     x133 * x69 - x174 * x69) *
      (x320 - x490) +
    (-12000000.0 * sc[1] * x42) * (-x106 + x495 + (-x69) * (x320 - x490));
  J[119] =
    12000000.0 * sc[1] * x202 + sc[3] * x70 - 26480000.0 * sc[7] * x386 +
    26480000.0 * sc[9] * x311 + x101 * x386 - x1140 * x501 + x1140 * x508 +
    32000000.0 * x136 * x51 + x1577 + 4.92 * x168 * exp(x2187) + x174 * x352 -
    x174 * x356 + x2099 + x2113 + x218 * x72 - x2190 * x501 + x2190 * x508 +
    x2193 + x999 +
    (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 - x101 * x124 -
     x133 * x69 - x174 * x69) *
      (-x501 + x508) +
    (-12000000.0 * sc[1] * x42) * (x352 - x356 + (-x69) * (-x501 + x508));
  J[133] =
    -4800000000000000.0 * refC * sc[7] * x1071 * x1086 +
    4800000000000000.0 * sc[1] * x1071 * x1074 - 26480000.0 * sc[7] * x392 +
    26480000.0 * sc[8] * x311 + x101 * x392 - 6.144e+20 * x1063 * x1105 +
    6.144e+20 * x1063 * x331 - x1140 * x148 + x1140 * x164 + x1185 -
    x148 * x2190 + x153 * (0.42 * x1082 + 3.0 * x1083) + x164 * x2190 -
    x174 * x27 + x174 * x359 + x2123 - x93 * (0.42 * x1082 + 3.0 * x1083) +
    (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 - x101 * x124 -
     x133 * x69 - x174 * x69) *
      (-x148 + x164) +
    (-12000000.0 * sc[1] * x42) * (-x27 + x359 + (-x69) * (-x148 + x164));
  J[147] =
    12000000.0 * sc[1] * x202 - 26480000.0 * sc[7] * x141 + x101 * x141 -
    x1140 * x178 + x1140 * x179 + x133 * x367 - x133 * x370 + x174 * x367 -
    x174 * x370 - x178 * x1973 - x178 * x2190 + x179 * x1973 + x179 * x2190 +
    x2193 + x439 +
    (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 - x101 * x124 -
     x133 * x69 - x174 * x69) *
      (-x178 + x179) +
    (-12000000.0 * sc[1] * x42) * (x367 - x370 + (-x69) * (-x178 + x179));
  J[161] =
    70000000.0 * sc[0] * x196 + sc[1] * x1310 + 12000000.0 * sc[1] * x244 +
    sc[4] * x1864 - sc[7] * x2215 - 26480000.0 * sc[7] * x354 + x101 * x354 -
    x106 * x133 - x1140 * x186 + x1140 * x197 + x133 * x196 + x174 * x196 -
    x186 * x1973 - x186 * x205 - x186 * x2190 + x197 * x1973 + x197 * x205 +
    x197 * x2190 + x2208 +
    (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 - x101 * x124 -
     x133 * x69 - x174 * x69) *
      (-x186 + x197) +
    (-12000000.0 * sc[1] * x42) * (-x106 + x196 + (-x69) * (-x186 + x197));
  J[175] = -26480000.0 * sc[7] * x119 + 32000000.0 * sc[8] * x136 * x54 +
           x101 * x119 - x1140 * x172 + x1140 * x190 - x172 * x2190 -
           x174 * x58 + x190 * x2190 + x2164 + x2196 + x463 +
           (-70000000.0 * sc[0] * x69 + 26480000.0 * sc[7] * x124 -
            x101 * x124 - x133 * x69 - x174 * x69) *
             (-x172 + x190) +
           (-12000000.0 * sc[1] * x42) * (x54 - x58 + (-x69) * (-x172 + x190));
  J[8] =
    26480000.0 * sc[7] * x142 + sc[7] * x338 - x101 * x142 + x1207 -
    x156 * x215 + x156 * x220 + x2239 * x342 - x2239 * x377 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (-x342 + x377);
  J[22] =
    sc[7] * x1213 * exp(x1196) + 26480000.0 * sc[7] * x209 - x101 * x209 +
    x1320 - x150 * x2239 + x156 * x286 - x156 * x407 + x2239 * x393 - x863 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (x150 - x393);
  J[36] =
    26480000.0 * sc[7] * x277 - sc[8] * x449 - x101 * x277 + x1320 -
    x156 * x290 + x156 * x335 + x2239 * x416 - x2239 * x424 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (-x416 + x424);
  J[50] =
    26480000.0 * sc[7] * x374 - sc[8] * x70 - x101 * x374 + x1320 -
    x156 * x341 + x156 * x343 + x2239 * x430 - x2239 * x435 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (-x430 + x435);
  J[64] =
    26480000.0 * sc[7] * x375 + sc[7] * x84 - x101 * x375 - x117 * x2239 +
    x1320 - x156 * x161 + x156 * x345 - x2019 + x2239 * x443 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (x117 - x443);
  J[78] =
    -1.4862e+22 * refC * sc[8] * x1182 * x1197 + sc[7] * x45 +
    26480000.0 * sc[7] * x95 - x101 * x95 +
    2.8986751181102398e+33 * x1179 * x1192 + 1.4862e+22 * x1182 * x1211 -
    2.8986751181102398e+33 * x1192 * x1201 - x1215 * x1492 + x1220 * x1492 -
    x156 * x346 + x156 * x349 - x170 * x2239 + x2239 * x466 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (x170 - x466);
  J[92] =
    26480000.0 * sc[7] * x281 + 1000000.0 * sc[7] - x101 * x281 + x106 * x156 +
    x1320 - x156 * x486 + x2239 * x478 - x2239 * x487 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (-x478 + x487);
  J[106] =
    sc[0] * x338 + sc[11] * x2215 + sc[1] * x1213 * exp(x1196) + sc[4] * x84 +
    sc[5] * x45 + 1000000.0 * sc[6] + 26480000.0 * sc[7] * x385 +
    4.92 * sc[7] * exp(-g_RT_qss[0] + x2187 - x271) +
    32000000.0 * sc[7] * exp(-g_RT_qss[1] + 286.83349950869001 * invT - x271) -
    x101 * x385 + x106 * x156 + x1320 - x156 * x495 - 26480000.0 * x171 -
    x2239 * x320 + x2239 * x490 - 26480000.0 * x316 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (x320 - x490);
  J[120] =
    -2.477e+21 * refC * x1182 * x1197 * x8 - sc[2] * x449 - sc[3] * x70 +
    26480000.0 * sc[7] * x386 - 26480000.0 * sc[9] * x311 - x101 * x386 +
    x1207 - 16000000.0 * x136 * x51 - x156 * x352 + x156 * x356 -
    2.46 * x168 * exp(x2187) - x2099 - x218 * x72 + x2239 * x501 -
    x2239 * x508 - x999 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (-x501 + x508);
  J[134] =
    -3.7155e+21 * refC * sc[8] * x1182 * x1197 + 26480000.0 * sc[7] * x392 -
    26480000.0 * sc[8] * x311 - x101 * x392 +
    7.2466877952755894e+32 * x1179 * x1192 + 3.7155e+21 * x1182 * x1211 -
    7.2466877952755894e+32 * x1192 * x1201 - x1211 * x1213 * x409 +
    x1218 * x409 * x8 + x148 * x2239 + x156 * x27 - x156 * x359 - x164 * x2239 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (-x148 + x164);
  J[148] =
    26480000.0 * sc[7] * x141 - x101 * x141 + x1207 - x156 * x367 +
    x156 * x370 + x178 * x2239 - x179 * x2239 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (-x178 + x179);
  J[162] =
    sc[7] * x2215 + 26480000.0 * sc[7] * x354 - x101 * x354 + x106 * x156 +
    x1320 - x156 * x196 + x186 * x2239 - x197 * x2239 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (-x186 + x197);
  J[176] =
    26480000.0 * sc[7] * x119 - 16000000.0 * sc[8] * x136 * x54 - x101 * x119 +
    x1320 + x156 * x58 + x172 * x2239 - x190 * x2239 +
    (-26480000.0 * sc[7] * x124 + x101 * x124 + x156 * x69) * (-x172 + x190);
  J[9] = 14000000.0 * sc[10] * x215 - 14000000.0 * sc[10] * x220 +
         30000000.0 * sc[2] * x142 + 40000000.0 * sc[3] * x215 -
         40000000.0 * sc[3] * x220 + 50000000.0 * sc[4] * x142 +
         26480000.0 * sc[7] * x142 + x1166 + 1204.0 * x12 - 1204.0 * x1239 +
         x142 * x437 + x142 * x442 + x142 * x484 + x848 + x855 +
         (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
          30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
          40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
          26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
           (-x342 + x377);
  J[23] =
    sc[10] * x1276 - 14000000.0 * sc[10] * x286 + 14000000.0 * sc[10] * x407 +
    30000000.0 * sc[2] * x209 - 40000000.0 * sc[3] * x286 +
    40000000.0 * sc[3] * x407 + 50000000.0 * sc[4] * x209 +
    26480000.0 * sc[7] * x209 - x1259 + x1303 + x209 * x484 + x2252 + x908 +
    (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
     40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
      (x150 - x393);
  J[37] =
    sc[10] * x18 + 14000000.0 * sc[10] * x290 - 14000000.0 * sc[10] * x335 +
    73400000.0 * sc[1] * x277 + 15000000.0 * sc[2] * x290 -
    15000000.0 * sc[2] * x335 + 40000000.0 * sc[3] * x290 -
    40000000.0 * sc[3] * x335 + 50000000.0 * sc[4] * x277 +
    26480000.0 * sc[7] * x277 - sc[9] * x1237 * x32 + x1382 + x2252 +
    x277 * x437 + x277 * x442 + x277 * x484 + 15000000.0 * x51 +
    (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
     40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
      (-x416 + x424);
  J[51] =
    14000000.0 * sc[10] * x341 - 14000000.0 * sc[10] * x343 +
    73400000.0 * sc[1] * x374 + 15000000.0 * sc[2] * x341 -
    15000000.0 * sc[2] * x343 + 30000000.0 * sc[2] * x374 +
    40000000.0 * sc[3] * x341 - 40000000.0 * sc[3] * x343 +
    50000000.0 * sc[4] * x374 + 26480000.0 * sc[7] * x374 + x1068 +
    3612.0 * x12 - 3612.0 * x1239 - x1744 +
    7600000.0 * x218 * exp(-201.28666632188799 * invT) + x374 * x437 +
    x374 * x442 + x374 * x484 + 40000000.0 * x51 + x901 +
    (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
     40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
      (-x430 + x435);
  J[65] =
    14000000.0 * sc[10] * x161 + sc[10] * x2029 - 14000000.0 * sc[10] * x345 +
    73400000.0 * sc[1] * x375 + 15000000.0 * sc[2] * x161 -
    15000000.0 * sc[2] * x345 + 30000000.0 * sc[2] * x375 +
    40000000.0 * sc[3] * x161 - 40000000.0 * sc[3] * x345 +
    50000000.0 * sc[4] * x375 + 26480000.0 * sc[7] * x375 -
    28000000.0 * sc[9] * x257 * x258 - 30000000.0 * sc[9] * x293 - x1420 -
    50000000.0 * x171 + x2252 - 50000000.0 * x316 + x375 * x437 + x375 * x442 +
    x375 * x484 +
    (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
     40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
      (x117 - x443);
  J[79] =
    14000000.0 * sc[10] * x346 - 14000000.0 * sc[10] * x349 +
    40000000.0 * sc[3] * x346 - 40000000.0 * sc[3] * x349 +
    50000000.0 * sc[4] * x95 + 26480000.0 * sc[7] * x95 -
    12000000.0 * sc[9] * x254 - 50000000.0 * sc[9] * x302 + 3612.0 * x12 -
    3612.0 * x1239 + x1485 + x484 * x95 + x970 +
    (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
     40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
      (x170 - x466);
  J[93] =
    14000000.0 * sc[10] * x486 + 73400000.0 * sc[1] * x281 +
    30000000.0 * sc[2] * x281 + 15000000.0 * sc[2] * x486 +
    40000000.0 * sc[3] * x486 + 50000000.0 * sc[4] * x281 +
    26480000.0 * sc[7] * x281 - sc[9] * x2071 - sc[9] * x87 + x2270 + x2271 +
    x281 * x437 + x281 * x442 + x281 * x484 +
    (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
     40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
      (-x478 + x487);
  J[107] =
    14000000.0 * sc[10] * x495 + 73400000.0 * sc[1] * x385 +
    30000000.0 * sc[2] * x385 + 15000000.0 * sc[2] * x495 +
    40000000.0 * sc[3] * x495 + 50000000.0 * sc[4] * x385 +
    2244000000000.0 * sc[5] * x385 * x434 + 26480000.0 * sc[7] * x385 -
    26480000.0 * x171 + x2270 + x2271 - 26480000.0 * x316 +
    187000000000.0 * x385 * x434 * x7 + x385 * x484 +
    (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
     40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
      (x320 - x490);
  J[121] =
    14000000.0 * sc[10] * x352 - 14000000.0 * sc[10] * x356 +
    73400000.0 * sc[1] * x386 + 15000000.0 * sc[2] * x352 -
    15000000.0 * sc[2] * x356 + 30000000.0 * sc[2] * x386 +
    40000000.0 * sc[3] * x352 - 40000000.0 * sc[3] * x356 +
    50000000.0 * sc[4] * x386 + 26480000.0 * sc[7] * x386 -
    26480000.0 * sc[9] * x311 + x1166 + 1204.0 * x12 - 1204.0 * x1239 +
    x386 * x437 + x386 * x442 + x386 * x484 + x848 +
    (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
     40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
      (-x501 + x508);
  J[135] =
    -14000000.0 * sc[10] * x27 + 14000000.0 * sc[10] * x359 -
    14000000.0 * sc[11] * x263 - sc[2] * x1237 * x32 -
    40000000.0 * sc[3] * x27 + 40000000.0 * sc[3] * x359 -
    28000000.0 * sc[4] * x257 * x258 - 30000000.0 * sc[4] * x293 +
    50000000.0 * sc[4] * x392 - 12000000.0 * sc[5] * x254 -
    50000000.0 * sc[5] * x302 - sc[6] * x2071 - sc[6] * x87 +
    26480000.0 * sc[7] * x392 - 26480000.0 * sc[8] * x311 + 903.0 * x12 -
    903.0 * x1239 + x125 - x1579 - x1785 + x392 * x484 + x824 +
    (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
     40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
      (-x148 + x164);
  J[149] =
    refC * x1233 * x32 + 14000000.0 * sc[10] * x367 -
    14000000.0 * sc[10] * x370 + sc[1] * x1276 + 73400000.0 * sc[1] * x141 +
    30000000.0 * sc[2] * x141 + sc[2] * x18 + 15000000.0 * sc[2] * x367 -
    15000000.0 * sc[2] * x370 + 40000000.0 * sc[3] * x367 -
    40000000.0 * sc[3] * x370 + 50000000.0 * sc[4] * x141 + sc[4] * x2029 +
    26480000.0 * sc[7] * x141 + x1166 + 2107.0 * x12 - 2107.0 * x1239 +
    x141 * x437 + x141 * x442 + x141 * x484 + 14000000.0 * x51 + x848 +
    (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
     40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
      (-x178 + x179);
  J[163] =
    14000000.0 * sc[10] * x196 + 73400000.0 * sc[1] * x354 +
    15000000.0 * sc[2] * x196 + 30000000.0 * sc[2] * x354 +
    40000000.0 * sc[3] * x196 + 50000000.0 * sc[4] * x354 +
    26480000.0 * sc[7] * x354 - 14000000.0 * sc[9] * x263 + x2270 + x2271 +
    x354 * x437 + x354 * x442 + x354 * x484 + x817 +
    (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
     40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
      (-x186 + x197);
  J[177] =
    14000000.0 * sc[10] * x54 - 14000000.0 * sc[10] * x58 +
    73400000.0 * sc[1] * x119 + 30000000.0 * sc[2] * x119 +
    15000000.0 * sc[2] * x54 - 15000000.0 * sc[2] * x58 +
    40000000.0 * sc[3] * x54 - 40000000.0 * sc[3] * x58 +
    50000000.0 * sc[4] * x119 + 26480000.0 * sc[7] * x119 + x119 * x437 +
    x119 * x442 + x119 * x484 + x2252 +
    (-14000000.0 * sc[10] * x69 - 73400000.0 * sc[1] * x124 -
     30000000.0 * sc[2] * x124 - 15000000.0 * sc[2] * x69 -
     40000000.0 * sc[3] * x69 - 50000000.0 * sc[4] * x124 -
     26480000.0 * sc[7] * x124 - x124 * x437 - x124 * x442 - x124 * x484) *
      (-x172 + x190);
  J[10] =
    -14000000.0 * sc[10] * x215 + 14000000.0 * sc[10] * x220 +
    30000000.0 * sc[2] * x142 - 1204.0 * x12 + 1204.0 * x1239 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (-x342 + x377);
  J[24] =
    -sc[10] * x1276 + 14000000.0 * sc[10] * x286 - 30000000.0 * sc[10] * x298 -
    14000000.0 * sc[10] * x407 + 30000000.0 * sc[2] * x209 + x2294 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (x150 - x393);
  J[38] =
    -sc[10] * x18 - 14000000.0 * sc[10] * x290 + 14000000.0 * sc[10] * x335 +
    sc[9] * x1237 * x32 + x1382 + x2294 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (-x416 + x424);
  J[52] =
    -14000000.0 * sc[10] * x341 + 14000000.0 * sc[10] * x343 +
    30000000.0 * sc[2] * x374 - 3612.0 * x12 + 3612.0 * x1239 + x1744 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (-x430 + x435);
  J[66] =
    -14000000.0 * sc[10] * x161 - sc[10] * x2029 + 14000000.0 * sc[10] * x345 +
    30000000.0 * sc[2] * x375 + x1420 + x2294 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (x117 - x443);
  J[80] =
    -14000000.0 * sc[10] * x346 + 14000000.0 * sc[10] * x349 +
    30000000.0 * sc[2] * x95 - 3612.0 * x12 + 3612.0 * x1239 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (x170 - x466);
  J[94] =
    14000000.0 * sc[10] * x106 - 14000000.0 * sc[10] * x486 +
    30000000.0 * sc[2] * x281 + sc[9] * x2071 + x2294 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (-x478 + x487);
  J[108] =
    14000000.0 * sc[10] * x106 - 14000000.0 * sc[10] * x495 +
    30000000.0 * sc[2] * x385 + x2294 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (x320 - x490);
  J[122] =
    -14000000.0 * sc[10] * x352 + 14000000.0 * sc[10] * x356 +
    30000000.0 * sc[2] * x386 - 1204.0 * x12 + 1204.0 * x1239 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (-x501 + x508);
  J[136] =
    14000000.0 * sc[10] * x27 - 14000000.0 * sc[10] * x359 +
    14000000.0 * sc[11] * x263 + sc[2] * x1237 * x32 +
    30000000.0 * sc[2] * x392 + sc[6] * x2071 - 903.0 * x12 + 903.0 * x1239 +
    x1579 + x1785 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (-x148 + x164);
  J[150] =
    -refC * x1233 * x32 - 14000000.0 * sc[10] * x367 +
    14000000.0 * sc[10] * x370 - sc[1] * x1276 - 30000000.0 * sc[1] * x298 +
    30000000.0 * sc[2] * x141 - sc[2] * x18 - sc[4] * x2029 - 2107.0 * x12 +
    2107.0 * x1239 - 14000000.0 * x51 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (-x178 + x179);
  J[164] =
    14000000.0 * sc[10] * x106 - 14000000.0 * sc[10] * x196 +
    30000000.0 * sc[2] * x354 + 14000000.0 * sc[9] * x263 + x2294 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (-x186 + x197);
  J[178] =
    -14000000.0 * sc[10] * x54 + 14000000.0 * sc[10] * x58 +
    30000000.0 * sc[2] * x119 + x2294 +
    (14000000.0 * sc[10] * x69 - 30000000.0 * sc[2] * x124) * (-x172 + x190);
  J[11] = 14000000.0 * sc[10] * x215 - 14000000.0 * sc[10] * x220 -
          20000000.0 * sc[11] * x235 + 20000000.0 * sc[1] * x315 +
          10000000.0 * sc[2] * x315 + 5000000.0 * sc[4] * x315 -
          20000000.0 * sc[6] * x342 + 20000000.0 * sc[6] * x377 +
          sc[9] * x278 * x8 + x101 * x142 + x109 * x315 + x1251 + x128 * x142 +
          x142 * x427 + x142 * x473 + x142 * x53 + x827 +
          (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
           x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
            (-x342 + x377) +
          (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
           5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
            (x215 - x220 + (-x69) * (-x342 + x377));
  J[25] = -14000000.0 * sc[10] * x286 + 14000000.0 * sc[10] * x407 +
          20000000.0 * sc[1] * x16 + 10000000.0 * sc[2] * x16 +
          5000000.0 * sc[4] * x16 + 20000000.0 * sc[6] * x150 -
          20000000.0 * sc[6] * x393 + x101 * x209 + x109 * x16 +
          1350000000000.0 * x1098 * x8 + x128 * x209 + x1319 + x1335 +
          x209 * x427 + x209 * x473 + x209 * x53 - 20000000.0 * x245 + x893 -
          20000000.0 * x909 +
          (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
           x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
            (x150 - x393) +
          (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
           5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
            (-x286 + x407 + (-x69) * (x150 - x393));
  J[39] = 14000000.0 * sc[10] * x290 - 14000000.0 * sc[10] * x335 -
          sc[11] * x1740 + 20000000.0 * sc[1] * x30 + 5000000.0 * sc[4] * x30 +
          x101 * x277 + x109 * x30 + x1335 + x1391 + x2004 + x277 * x427 +
          x277 * x473 + x277 * x53 + x893 +
          (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
           x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
            (-x416 + x424) +
          (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
           5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
            (x290 - x335 + (-x69) * (-x416 + x424));
  J[53] = 14000000.0 * sc[10] * x341 - 14000000.0 * sc[10] * x343 -
          sc[11] * x1860 + 20000000.0 * sc[1] * x382 +
          5000000.0 * sc[4] * x382 + x101 * x374 + x109 * x382 +
          18000000.0 * x1232 * exp(-452.89499922424801 * invT) + x1335 + x1418 +
          x2015 + x374 * x427 + x374 * x473 + x374 * x53 + x893 +
          (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
           x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
            (-x430 + x435) +
          (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
           5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
            (x341 - x343 + (-x69) * (-x430 + x435));
  J[67] = 14000000.0 * sc[10] * x161 - 14000000.0 * sc[10] * x345 +
          20000000.0 * sc[1] * x390 + 5000000.0 * sc[4] * x390 + x101 * x375 +
          x109 * x390 + x1335 + x1447 + x2047 + x218 * x472 - 5000000.0 * x245 +
          x375 * x427 + x375 * x473 + x375 * x53 + x893 - 5000000.0 * x909 +
          (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
           x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
            (x117 - x443) +
          (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
           5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
            (x161 - x345 + (-x69) * (x117 - x443));
  J[81] = 14000000.0 * sc[10] * x346 - 14000000.0 * sc[10] * x349 -
          5000000.0 * sc[11] * x226 + 8100000000000.0 * sc[1] * x1098 +
          20000000.0 * sc[1] * x395 + 5000000.0 * sc[4] * x395 + x101 * x95 -
          x1065 * x1097 * x1100 + x1065 * x1097 * x423 + x109 * x395 + x1506 +
          x2067 - 8100000000000.0 * x283 * x405 +
          1.0032110091743101e+19 * x285 - x376 * x964 + 3.042e+16 * x389 * x67 -
          3.042e+16 * x389 * x819 + x427 * x95 + x473 * x95 + x53 * x95 +
          3.5867302325581401e+30 * x55 * x831 -
          3.5867302325581401e+30 * x55 * x844 + x821 * x964 -
          1.0032110091743101e+19 * x98 +
          (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
           x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
            (x170 - x466) +
          (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
           5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
            (x346 - x349 + (-x69) * (x170 - x466));
  J[95] = 14000000.0 * sc[10] * x486 - sc[11] * x1890 +
          20000000.0 * sc[1] * x224 - 30000000.0 * sc[4] * x106 +
          5000000.0 * sc[4] * x224 + x101 * x281 + x109 * x224 + x1335 + x1530 +
          x2084 + x218 * x476 + x2270 + x281 * x427 + x281 * x473 + x281 * x53 +
          x893 +
          (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
           x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
            (-x478 + x487) +
          (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
           5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
            (-x106 + x486 + (-x69) * (-x478 + x487));
  J[109] =
    14000000.0 * sc[10] * x495 - sc[11] * x2215 + 20000000.0 * sc[1] * x399 -
    30000000.0 * sc[4] * x106 + 5000000.0 * sc[4] * x399 + x101 * x385 +
    x109 * x399 + x1335 + x1553 + x2092 + x2270 + x385 * x427 + x385 * x473 +
    x385 * x53 + x893 +
    (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
     x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
      (x320 - x490) +
    (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
     5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
      (-x106 + x495 + (-x69) * (x320 - x490));
  J[123] =
    14000000.0 * sc[10] * x352 - 14000000.0 * sc[10] * x356 +
    20000000.0 * sc[1] * x202 + 10000000.0 * sc[2] * x202 +
    5000000.0 * sc[4] * x202 - 20000000.0 * sc[6] * x501 +
    20000000.0 * sc[6] * x508 + x101 * x386 + x109 * x202 + x128 * x386 +
    x1575 + x218 * x72 + x386 * x427 + x386 * x473 + x386 * x53 + x827 +
    (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
     x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
      (-x501 + x508) +
    (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
     5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
      (x352 - x356 + (-x69) * (-x501 + x508));
  J[137] =
    sc[0] * x278 * x8 - 14000000.0 * sc[10] * x27 + 14000000.0 * sc[10] * x359 -
    14000000.0 * sc[11] * x263 + 2025000000000.0 * sc[1] * x1098 +
    20000000.0 * sc[1] * x204 + 10000000.0 * sc[2] * x204 +
    5000000.0 * sc[4] * x204 - 20000000.0 * sc[6] * x148 +
    20000000.0 * sc[6] * x164 + x101 * x392 + x109 * x204 - x1094 * x1100 +
    x1094 * x423 + x11 + x128 * x392 - 2025000000000.0 * x283 * x405 +
    2.5080275229357798e+18 * x285 - x376 * x823 * x826 +
    7605000000000000.0 * x389 * x67 - 7605000000000000.0 * x389 * x819 +
    x392 * x427 + x392 * x473 + x392 * x53 +
    8.9668255813953501e+29 * x55 * x831 - 8.9668255813953501e+29 * x55 * x844 +
    x821 * x823 * x826 - 2.5080275229357798e+18 * x98 +
    (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
     x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
      (-x148 + x164) +
    (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
     5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
      (-x27 + x359 + (-x69) * (-x148 + x164));
  J[151] =
    14000000.0 * sc[10] * x367 - 14000000.0 * sc[10] * x370 +
    20000000.0 * sc[1] * x202 + 10000000.0 * sc[2] * x202 +
    5000000.0 * sc[4] * x202 - 20000000.0 * sc[6] * x178 +
    20000000.0 * sc[6] * x179 + x101 * x141 + x109 * x202 + x1163 +
    x128 * x141 + x141 * x427 + x141 * x473 + x141 * x53 + 14000000.0 * x51 +
    x827 +
    (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
     x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
      (-x178 + x179) +
    (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
     5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
      (x367 - x370 + (-x69) * (-x178 + x179));
  J[165] =
    -1350000000000.0 * refC * x405 * x8 - 20000000.0 * sc[0] * x235 +
    14000000.0 * sc[10] * x196 + 20000000.0 * sc[1] * x244 - sc[2] * x1740 -
    sc[3] * x1860 - 30000000.0 * sc[4] * x106 + 5000000.0 * sc[4] * x244 -
    5000000.0 * sc[5] * x226 - sc[6] * x1890 - sc[7] * x2215 -
    14000000.0 * sc[9] * x263 + x101 * x354 + x109 * x244 + x1335 + x2153 +
    x2270 + x354 * x427 + x354 * x473 + x354 * x53 + x458 - x817 + x893 +
    (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
     x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
      (-x186 + x197) +
    (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
     5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
      (-x106 + x196 + (-x69) * (-x186 + x197));
  J[179] =
    14000000.0 * sc[10] * x54 - 14000000.0 * sc[10] * x58 +
    20000000.0 * sc[1] * x34 + 10000000.0 * sc[2] * x34 +
    5000000.0 * sc[4] * x34 - 20000000.0 * sc[6] * x172 +
    20000000.0 * sc[6] * x190 + x101 * x119 + x109 * x34 + x1186 + x119 * x128 +
    x119 * x427 + x119 * x473 + x119 * x53 + x893 +
    (-14000000.0 * sc[10] * x69 - 30000000.0 * sc[4] * x69 - x101 * x124 -
     x124 * x128 - x124 * x372 - x124 * x427 - x124 * x473 - x124 * x53) *
      (-x172 + x190) +
    (-20000000.0 * sc[1] * x42 - 10000000.0 * sc[2] * x42 -
     5000000.0 * sc[4] * x42 - x109 * x42 - x42 * x468) *
      (x54 - x58 + (-x69) * (-x172 + x190));
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
  COFD[8] = -1.43693219E+01;
  COFD[9] = 4.03992999E+00;
  COFD[10] = -3.08044800E-01;
  COFD[11] = 1.32757775E-02;
  COFD[12] = -1.06250314E+01;
  COFD[13] = 2.15849701E+00;
  COFD[14] = -6.53886401E-02;
  COFD[15] = 2.81453370E-03;
  COFD[16] = -1.31860264E+01;
  COFD[17] = 3.38003453E+00;
  COFD[18] = -2.25783856E-01;
  COFD[19] = 9.85028660E-03;
  COFD[20] = -1.29877255E+01;
  COFD[21] = 2.80841511E+00;
  COFD[22] = -1.52629888E-01;
  COFD[23] = 6.72604927E-03;
  COFD[24] = -1.15797897E+01;
  COFD[25] = 2.43235504E+00;
  COFD[26] = -1.02890179E-01;
  COFD[27] = 4.52903603E-03;
  COFD[28] = -1.43713019E+01;
  COFD[29] = 3.70920439E+00;
  COFD[30] = -2.67274113E-01;
  COFD[31] = 1.15967481E-02;
  COFD[32] = -1.40864784E+01;
  COFD[33] = 3.07458927E+00;
  COFD[34] = -1.86899591E-01;
  COFD[35] = 8.19829781E-03;
  COFD[36] = -1.53110597E+01;
  COFD[37] = 3.37317428E+00;
  COFD[38] = -2.24900439E-01;
  COFD[39] = 9.81228151E-03;
  COFD[40] = -1.06283588E+01;
  COFD[41] = 2.15849701E+00;
  COFD[42] = -6.53886401E-02;
  COFD[43] = 2.81453370E-03;
  COFD[44] = -1.31877859E+01;
  COFD[45] = 3.38003453E+00;
  COFD[46] = -2.25783856E-01;
  COFD[47] = 9.85028660E-03;
  COFD[48] = -1.30027669E+01;
  COFD[49] = 2.80841511E+00;
  COFD[50] = -1.52629888E-01;
  COFD[51] = 6.72604927E-03;
  COFD[52] = -1.41066358E+01;
  COFD[53] = 3.07458927E+00;
  COFD[54] = -1.86899591E-01;
  COFD[55] = 8.19829781E-03;
  COFD[56] = -1.30182749E+01;
  COFD[57] = 2.80841511E+00;
  COFD[58] = -1.52629888E-01;
  COFD[59] = 6.72604927E-03;
  COFD[60] = -1.68759075E+01;
  COFD[61] = 4.49460341E+00;
  COFD[62] = -3.64766179E-01;
  COFD[63] = 1.56457173E-02;
  COFD[64] = -1.93611201E+01;
  COFD[65] = 5.51579726E+00;
  COFD[66] = -4.76061961E-01;
  COFD[67] = 1.96329391E-02;
  COFD[68] = -1.91096700E+01;
  COFD[69] = 5.02608697E+00;
  COFD[70] = -4.26959993E-01;
  COFD[71] = 1.80709910E-02;
  COFD[72] = -2.10639926E+01;
  COFD[73] = 5.50980703E+00;
  COFD[74] = -4.78335492E-01;
  COFD[75] = 1.98515434E-02;
  COFD[76] = -1.91256173E+01;
  COFD[77] = 5.02608697E+00;
  COFD[78] = -4.26959993E-01;
  COFD[79] = 1.80709910E-02;
  COFD[80] = -1.31492509E+01;
  COFD[81] = 1.48004045E+00;
  COFD[82] = 1.60499956E-01;
  COFD[83] = -1.19765871E-02;
  COFD[84] = -1.15806955E+01;
  COFD[85] = 2.43235504E+00;
  COFD[86] = -1.02890179E-01;
  COFD[87] = 4.52903603E-03;
  COFD[88] = -1.43717684E+01;
  COFD[89] = 3.70920439E+00;
  COFD[90] = -2.67274113E-01;
  COFD[91] = 1.15967481E-02;
  COFD[92] = -1.40915944E+01;
  COFD[93] = 3.07458927E+00;
  COFD[94] = -1.86899591E-01;
  COFD[95] = 8.19829781E-03;
  COFD[96] = -1.53187536E+01;
  COFD[97] = 3.37317428E+00;
  COFD[98] = -2.24900439E-01;
  COFD[99] = 9.81228151E-03;
  COFD[100] = -1.41119635E+01;
  COFD[101] = 3.07458927E+00;
  COFD[102] = -1.86899591E-01;
  COFD[103] = 8.19829781E-03;
  COFD[104] = -2.04177394E+01;
  COFD[105] = 5.31457079E+00;
  COFD[106] = -4.58216496E-01;
  COFD[107] = 1.91825910E-02;
  COFD[108] = -1.53265677E+01;
  COFD[109] = 3.37317428E+00;
  COFD[110] = -2.24900439E-01;
  COFD[111] = 9.81228151E-03;
  COFD[112] = -1.22583533E+01;
  COFD[113] = 2.66488613E+00;
  COFD[114] = -1.33839930E-01;
  COFD[115] = 5.90524087E-03;
  COFD[116] = -1.53791907E+01;
  COFD[117] = 4.03496578E+00;
  COFD[118] = -3.07489899E-01;
  COFD[119] = 1.32554084E-02;
  COFD[120] = -1.48358026E+01;
  COFD[121] = 3.37160686E+00;
  COFD[122] = -2.24698307E-01;
  COFD[123] = 9.80357301E-03;
  COFD[124] = -1.60662884E+01;
  COFD[125] = 3.69744368E+00;
  COFD[126] = -2.65754489E-01;
  COFD[127] = 1.15311681E-02;
  COFD[128] = -1.48503699E+01;
  COFD[129] = 3.37160686E+00;
  COFD[130] = -2.24698307E-01;
  COFD[131] = 9.80357301E-03;
  COFD[132] = -2.10691956E+01;
  COFD[133] = 5.51604369E+00;
  COFD[134] = -4.76337283E-01;
  COFD[135] = 1.96544379E-02;
  COFD[136] = -1.60711937E+01;
  COFD[137] = 3.69744368E+00;
  COFD[138] = -2.65754489E-01;
  COFD[139] = 1.15311681E-02;
  COFD[140] = -1.69251551E+01;
  COFD[141] = 4.03028423E+00;
  COFD[142] = -3.06976249E-01;
  COFD[143] = 1.32370318E-02;
  COFD[144] = -1.22114409E+01;
  COFD[145] = 2.65215786E+00;
  COFD[146] = -1.32206766E-01;
  COFD[147] = 5.83532580E-03;
  COFD[148] = -1.53241298E+01;
  COFD[149] = 4.02249494E+00;
  COFD[150] = -3.06112211E-01;
  COFD[151] = 1.32056535E-02;
  COFD[152] = -1.47765660E+01;
  COFD[153] = 3.34934830E+00;
  COFD[154] = -2.21831026E-01;
  COFD[155] = 9.68019791E-03;
  COFD[156] = -1.59878735E+01;
  COFD[157] = 3.66478157E+00;
  COFD[158] = -2.61506432E-01;
  COFD[159] = 1.13466231E-02;
  COFD[160] = -1.47916284E+01;
  COFD[161] = 3.34934830E+00;
  COFD[162] = -2.21831026E-01;
  COFD[163] = 9.68019791E-03;
  COFD[164] = -2.12953626E+01;
  COFD[165] = 5.52385262E+00;
  COFD[166] = -4.69683817E-01;
  COFD[167] = 1.90677254E-02;
  COFD[168] = -1.59929990E+01;
  COFD[169] = 3.66478157E+00;
  COFD[170] = -2.61506432E-01;
  COFD[171] = 1.13466231E-02;
  COFD[172] = -1.68880017E+01;
  COFD[173] = 4.01759713E+00;
  COFD[174] = -3.05560618E-01;
  COFD[175] = 1.31852111E-02;
  COFD[176] = -1.68532873E+01;
  COFD[177] = 4.00572564E+00;
  COFD[178] = -3.04255586E-01;
  COFD[179] = 1.31384083E-02;
  COFD[180] = -1.13541221E+01;
  COFD[181] = 2.31999438E+00;
  COFD[182] = -8.75064804E-02;
  COFD[183] = 3.82656365E-03;
  COFD[184] = -1.40524219E+01;
  COFD[185] = 3.56261348E+00;
  COFD[186] = -2.48287981E-01;
  COFD[187] = 1.07752947E-02;
  COFD[188] = -1.39007309E+01;
  COFD[189] = 2.99164244E+00;
  COFD[190] = -1.76293106E-01;
  COFD[191] = 7.74575100E-03;
  COFD[192] = -1.50371688E+01;
  COFD[193] = 3.26249588E+00;
  COFD[194] = -2.10658287E-01;
  COFD[195] = 9.20032462E-03;
  COFD[196] = -1.39199573E+01;
  COFD[197] = 2.99164244E+00;
  COFD[198] = -1.76293106E-01;
  COFD[199] = 7.74575100E-03;
  COFD[200] = -2.08943723E+01;
  COFD[201] = 5.44718673E+00;
  COFD[202] = -4.72082974E-01;
  COFD[203] = 1.96531328E-02;
  COFD[204] = -1.50443476E+01;
  COFD[205] = 3.26249588E+00;
  COFD[206] = -2.10658287E-01;
  COFD[207] = 9.20032462E-03;
  COFD[208] = -1.57272832E+01;
  COFD[209] = 3.55756440E+00;
  COFD[210] = -2.47689003E-01;
  COFD[211] = 1.07517019E-02;
  COFD[212] = -1.56953785E+01;
  COFD[213] = 3.54443207E+00;
  COFD[214] = -2.46133003E-01;
  COFD[215] = 1.06905018E-02;
  COFD[216] = -1.48061406E+01;
  COFD[217] = 3.16912473E+00;
  COFD[218] = -1.98792456E-01;
  COFD[219] = 8.69726395E-03;
  COFD[220] = -1.35545391E+01;
  COFD[221] = 3.13878730E+00;
  COFD[222] = -1.94980335E-01;
  COFD[223] = 8.53744486E-03;
  COFD[224] = -1.72994129E+01;
  COFD[225] = 4.71931868E+00;
  COFD[226] = -3.91258152E-01;
  COFD[227] = 1.66866639E-02;
  COFD[228] = -1.67115471E+01;
  COFD[229] = 3.98859394E+00;
  COFD[230] = -3.02316219E-01;
  COFD[231] = 1.30661099E-02;
  COFD[232] = -1.81197250E+01;
  COFD[233] = 4.33684042E+00;
  COFD[234] = -3.44981265E-01;
  COFD[235] = 1.48142449E-02;
  COFD[236] = -1.67337674E+01;
  COFD[237] = 3.98859394E+00;
  COFD[238] = -3.02316219E-01;
  COFD[239] = 1.30661099E-02;
  COFD[240] = -2.12021420E+01;
  COFD[241] = 5.20775027E+00;
  COFD[242] = -4.07348285E-01;
  COFD[243] = 1.55473262E-02;
  COFD[244] = -1.81286456E+01;
  COFD[245] = 4.33684042E+00;
  COFD[246] = -3.44981265E-01;
  COFD[247] = 1.48142449E-02;
  COFD[248] = -1.90110170E+01;
  COFD[249] = 4.70842146E+00;
  COFD[250] = -3.89924637E-01;
  COFD[251] = 1.66320934E-02;
  COFD[252] = -1.89594762E+01;
  COFD[253] = 4.68735220E+00;
  COFD[254] = -3.87449514E-01;
  COFD[255] = 1.65352261E-02;
  COFD[256] = -1.77672912E+01;
  COFD[257] = 4.20234040E+00;
  COFD[258] = -3.28057658E-01;
  COFD[259] = 1.41006192E-02;
  COFD[260] = -2.10907633E+01;
  COFD[261] = 5.29211327E+00;
  COFD[262] = -4.56068366E-01;
  COFD[263] = 1.91195062E-02;
  COFD[264] = -1.57172268E+01;
  COFD[265] = 3.96062263E+00;
  COFD[266] = -2.98964970E-01;
  COFD[267] = 1.29322565E-02;
  COFD[268] = -1.95318517E+01;
  COFD[269] = 5.47046983E+00;
  COFD[270] = -4.74577605E-01;
  COFD[271] = 1.97408822E-02;
  COFD[272] = -1.91105534E+01;
  COFD[273] = 4.87977047E+00;
  COFD[274] = -4.10448693E-01;
  COFD[275] = 1.74535827E-02;
  COFD[276] = -2.05196205E+01;
  COFD[277] = 5.21473296E+00;
  COFD[278] = -4.48646311E-01;
  COFD[279] = 1.89013813E-02;
  COFD[280] = -1.91302705E+01;
  COFD[281] = 4.87977047E+00;
  COFD[282] = -4.10448693E-01;
  COFD[283] = 1.74535827E-02;
  COFD[284] = -1.87236049E+01;
  COFD[285] = 4.00967621E+00;
  COFD[286] = -2.21153539E-01;
  COFD[287] = 6.31528745E-03;
  COFD[288] = -2.05270679E+01;
  COFD[289] = 5.21473296E+00;
  COFD[290] = -4.48646311E-01;
  COFD[291] = 1.89013813E-02;
  COFD[292] = -2.12140281E+01;
  COFD[293] = 5.46309033E+00;
  COFD[294] = -4.73777184E-01;
  COFD[295] = 1.97122249E-02;
  COFD[296] = -2.11618589E+01;
  COFD[297] = 5.44451952E+00;
  COFD[298] = -4.71804032E-01;
  COFD[299] = 1.96436657E-02;
  COFD[300] = -2.02444913E+01;
  COFD[301] = 5.11785645E+00;
  COFD[302] = -4.37867828E-01;
  COFD[303] = 1.85047543E-02;
  COFD[304] = -2.20860863E+01;
  COFD[305] = 5.52171573E+00;
  COFD[306] = -4.63284984E-01;
  COFD[307] = 1.85570924E-02;
  COFD[308] = -1.99154422E+01;
  COFD[309] = 4.38041133E+00;
  COFD[310] = -2.77538214E-01;
  COFD[311] = 9.06748822E-03;
  COFD[312] = -1.13253618E+01;
  COFD[313] = 2.31195095E+00;
  COFD[314] = -8.63988037E-02;
  COFD[315] = 3.77573452E-03;
  COFD[316] = -1.40298991E+01;
  COFD[317] = 3.55837688E+00;
  COFD[318] = -2.47785790E-01;
  COFD[319] = 1.07555332E-02;
  COFD[320] = -1.38756380E+01;
  COFD[321] = 2.98558426E+00;
  COFD[322] = -1.75507216E-01;
  COFD[323] = 7.71173691E-03;
  COFD[324] = -1.50096252E+01;
  COFD[325] = 3.25515933E+00;
  COFD[326] = -2.09710110E-01;
  COFD[327] = 9.15941830E-03;
  COFD[328] = -1.38948654E+01;
  COFD[329] = 2.98558426E+00;
  COFD[330] = -1.75507216E-01;
  COFD[331] = 7.71173691E-03;
  COFD[332] = -2.08123331E+01;
  COFD[333] = 5.42470175E+00;
  COFD[334] = -4.69700438E-01;
  COFD[335] = 1.95706912E-02;
  COFD[336] = -1.50168047E+01;
  COFD[337] = 3.25515933E+00;
  COFD[338] = -2.09710110E-01;
  COFD[339] = 9.15941830E-03;
  COFD[340] = -1.57070757E+01;
  COFD[341] = 3.55327776E+00;
  COFD[342] = -2.47179250E-01;
  COFD[343] = 1.07315653E-02;
  COFD[344] = -1.56756124E+01;
  COFD[345] = 3.54035770E+00;
  COFD[346] = -2.45653442E-01;
  COFD[347] = 1.06717969E-02;
  COFD[348] = -1.47850505E+01;
  COFD[349] = 3.16433919E+00;
  COFD[350] = -1.98191564E-01;
  COFD[351] = 8.67209742E-03;
  COFD[352] = -1.77350630E+01;
  COFD[353] = 4.19328271E+00;
  COFD[354] = -3.26911461E-01;
  COFD[355] = 1.40520357E-02;
  COFD[356] = -2.02136057E+01;
  COFD[357] = 5.10993120E+00;
  COFD[358] = -4.36931630E-01;
  COFD[359] = 1.84677592E-02;
  COFD[360] = -1.47639411E+01;
  COFD[361] = 3.15955654E+00;
  COFD[362] = -1.97590757E-01;
  COFD[363] = 8.64692156E-03;
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
