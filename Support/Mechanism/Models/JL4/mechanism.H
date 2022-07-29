#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  C
1  O
2  H
3  N
*/

// Species
#define CH4_ID 0
#define O2_ID 1
#define H2O_ID 2
#define N2_ID 3
#define CO_ID 4
#define CO2_ID 5
#define H2_ID 6

#define NUM_ELEMENTS 4
#define NUM_SPECIES 7
#define NUM_IONS 0
#define NUM_REACTIONS 4

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
  *mm = 4;
  *kk = 7;
  *ii = 4;
  *nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_imw(amrex::Real* imw_new)
{
  imw_new[0] = 1.0 / 16.043000; // CH4
  imw_new[1] = 1.0 / 31.998000; // O2
  imw_new[2] = 1.0 / 18.015000; // H2O
  imw_new[3] = 1.0 / 28.014000; // N2
  imw_new[4] = 1.0 / 28.010000; // CO
  imw_new[5] = 1.0 / 44.009000; // CO2
  imw_new[6] = 1.0 / 2.016000;  // H2
}

//  molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_mw(amrex::Real* mw_new)
{
  mw_new[0] = 16.043000; // CH4
  mw_new[1] = 31.998000; // O2
  mw_new[2] = 18.015000; // H2O
  mw_new[3] = 28.014000; // N2
  mw_new[4] = 28.010000; // CO
  mw_new[5] = 44.009000; // CO2
  mw_new[6] = 2.016000;  // H2
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
    // species 0: CH4
    species[0] = +4.14987613e+00 - 1.36709788e-02 * tc[1] +
                 4.91800599e-05 * tc[2] - 4.84743026e-08 * tc[3] +
                 1.66693956e-11 * tc[4];
    // species 1: O2
    species[1] = +2.78245636e+00 - 2.99673416e-03 * tc[1] +
                 9.84730201e-06 * tc[2] - 9.68129509e-09 * tc[3] +
                 3.24372837e-12 * tc[4];
    // species 2: H2O
    species[2] = +3.19864056e+00 - 2.03643410e-03 * tc[1] +
                 6.52040211e-06 * tc[2] - 5.48797062e-09 * tc[3] +
                 1.77197817e-12 * tc[4];
    // species 3: N2
    species[3] = +2.29867700e+00 + 1.40824040e-03 * tc[1] -
                 3.96322200e-06 * tc[2] + 5.64151500e-09 * tc[3] -
                 2.44485400e-12 * tc[4];
    // species 4: CO
    species[4] = +2.57953347e+00 - 6.10353680e-04 * tc[1] +
                 1.01681433e-06 * tc[2] + 9.07005884e-10 * tc[3] -
                 9.04424499e-13 * tc[4];
    // species 5: CO2
    species[5] = +1.35677352e+00 + 8.98459677e-03 * tc[1] -
                 7.12356269e-06 * tc[2] + 2.45919022e-09 * tc[3] -
                 1.43699548e-13 * tc[4];
    // species 6: H2
    species[6] = +1.34433112e+00 + 7.98052075e-03 * tc[1] -
                 1.94781510e-05 * tc[2] + 2.01572094e-08 * tc[3] -
                 7.37611761e-12 * tc[4];
  } else {
    // species 0: CH4
    species[0] = -9.25148505e-01 + 1.33909467e-02 * tc[1] -
                 5.73285809e-06 * tc[2] + 1.22292535e-09 * tc[3] -
                 1.01815230e-13 * tc[4];
    // species 1: O2
    species[1] = +2.28253784e+00 + 1.48308754e-03 * tc[1] -
                 7.57966669e-07 * tc[2] + 2.09470555e-10 * tc[3] -
                 2.16717794e-14 * tc[4];
    // species 2: H2O
    species[2] = +2.03399249e+00 + 2.17691804e-03 * tc[1] -
                 1.64072518e-07 * tc[2] - 9.70419870e-11 * tc[3] +
                 1.68200992e-14 * tc[4];
    // species 3: N2
    species[3] = +1.92664000e+00 + 1.48797680e-03 * tc[1] -
                 5.68476000e-07 * tc[2] + 1.00970380e-10 * tc[3] -
                 6.75335100e-15 * tc[4];
    // species 4: CO
    species[4] = +1.71518561e+00 + 2.06252743e-03 * tc[1] -
                 9.98825771e-07 * tc[2] + 2.30053008e-10 * tc[3] -
                 2.03647716e-14 * tc[4];
    // species 5: CO2
    species[5] = +2.85746029e+00 + 4.41437026e-03 * tc[1] -
                 2.21481404e-06 * tc[2] + 5.23490188e-10 * tc[3] -
                 4.72084164e-14 * tc[4];
    // species 6: H2
    species[6] = +2.33727920e+00 - 4.94024731e-05 * tc[1] +
                 4.99456778e-07 * tc[2] - 1.79566394e-10 * tc[3] +
                 2.00255376e-14 * tc[4];
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
    // species 0: CH4
    species[0] = +5.14987613e+00 - 1.36709788e-02 * tc[1] +
                 4.91800599e-05 * tc[2] - 4.84743026e-08 * tc[3] +
                 1.66693956e-11 * tc[4];
    // species 1: O2
    species[1] = +3.78245636e+00 - 2.99673416e-03 * tc[1] +
                 9.84730201e-06 * tc[2] - 9.68129509e-09 * tc[3] +
                 3.24372837e-12 * tc[4];
    // species 2: H2O
    species[2] = +4.19864056e+00 - 2.03643410e-03 * tc[1] +
                 6.52040211e-06 * tc[2] - 5.48797062e-09 * tc[3] +
                 1.77197817e-12 * tc[4];
    // species 3: N2
    species[3] = +3.29867700e+00 + 1.40824040e-03 * tc[1] -
                 3.96322200e-06 * tc[2] + 5.64151500e-09 * tc[3] -
                 2.44485400e-12 * tc[4];
    // species 4: CO
    species[4] = +3.57953347e+00 - 6.10353680e-04 * tc[1] +
                 1.01681433e-06 * tc[2] + 9.07005884e-10 * tc[3] -
                 9.04424499e-13 * tc[4];
    // species 5: CO2
    species[5] = +2.35677352e+00 + 8.98459677e-03 * tc[1] -
                 7.12356269e-06 * tc[2] + 2.45919022e-09 * tc[3] -
                 1.43699548e-13 * tc[4];
    // species 6: H2
    species[6] = +2.34433112e+00 + 7.98052075e-03 * tc[1] -
                 1.94781510e-05 * tc[2] + 2.01572094e-08 * tc[3] -
                 7.37611761e-12 * tc[4];
  } else {
    // species 0: CH4
    species[0] = +7.48514950e-02 + 1.33909467e-02 * tc[1] -
                 5.73285809e-06 * tc[2] + 1.22292535e-09 * tc[3] -
                 1.01815230e-13 * tc[4];
    // species 1: O2
    species[1] = +3.28253784e+00 + 1.48308754e-03 * tc[1] -
                 7.57966669e-07 * tc[2] + 2.09470555e-10 * tc[3] -
                 2.16717794e-14 * tc[4];
    // species 2: H2O
    species[2] = +3.03399249e+00 + 2.17691804e-03 * tc[1] -
                 1.64072518e-07 * tc[2] - 9.70419870e-11 * tc[3] +
                 1.68200992e-14 * tc[4];
    // species 3: N2
    species[3] = +2.92664000e+00 + 1.48797680e-03 * tc[1] -
                 5.68476000e-07 * tc[2] + 1.00970380e-10 * tc[3] -
                 6.75335100e-15 * tc[4];
    // species 4: CO
    species[4] = +2.71518561e+00 + 2.06252743e-03 * tc[1] -
                 9.98825771e-07 * tc[2] + 2.30053008e-10 * tc[3] -
                 2.03647716e-14 * tc[4];
    // species 5: CO2
    species[5] = +3.85746029e+00 + 4.41437026e-03 * tc[1] -
                 2.21481404e-06 * tc[2] + 5.23490188e-10 * tc[3] -
                 4.72084164e-14 * tc[4];
    // species 6: H2
    species[6] = +3.33727920e+00 - 4.94024731e-05 * tc[1] +
                 4.99456778e-07 * tc[2] - 1.79566394e-10 * tc[3] +
                 2.00255376e-14 * tc[4];
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
    // species 0: CH4
    species[0] = -1.024664760000000e+04 * invT + 9.791179889999999e+00 -
                 5.149876130000000e+00 * tc[0] + 6.835489400000000e-03 * tc[1] -
                 8.196676650000000e-06 * tc[2] + 4.039525216666667e-09 * tc[3] -
                 8.334697800000000e-13 * tc[4];
    // species 1: O2
    species[1] = -1.063943560000000e+03 * invT + 1.247806300000001e-01 -
                 3.782456360000000e+00 * tc[0] + 1.498367080000000e-03 * tc[1] -
                 1.641217001666667e-06 * tc[2] + 8.067745908333334e-10 * tc[3] -
                 1.621864185000000e-13 * tc[4];
    // species 2: H2O
    species[2] = -3.029372670000000e+04 * invT + 5.047672768000000e+00 -
                 4.198640560000000e+00 * tc[0] + 1.018217050000000e-03 * tc[1] -
                 1.086733685000000e-06 * tc[2] + 4.573308850000000e-10 * tc[3] -
                 8.859890850000000e-14 * tc[4];
    // species 3: N2
    species[3] = -1.020899900000000e+03 * invT - 6.516950000000001e-01 -
                 3.298677000000000e+00 * tc[0] - 7.041202000000000e-04 * tc[1] +
                 6.605369999999999e-07 * tc[2] - 4.701262500000001e-10 * tc[3] +
                 1.222427000000000e-13 * tc[4];
    // species 4: CO
    species[4] = -1.434408600000000e+04 * invT + 7.112418999999992e-02 -
                 3.579533470000000e+00 * tc[0] + 3.051768400000000e-04 * tc[1] -
                 1.694690550000000e-07 * tc[2] - 7.558382366666667e-11 * tc[3] +
                 4.522122495000000e-14 * tc[4];
    // species 5: CO2
    species[5] = -4.837196970000000e+04 * invT - 7.544278700000000e+00 -
                 2.356773520000000e+00 * tc[0] - 4.492298385000000e-03 * tc[1] +
                 1.187260448333333e-06 * tc[2] - 2.049325183333333e-10 * tc[3] +
                 7.184977399999999e-15 * tc[4];
    // species 6: H2
    species[6] = -9.179351730000000e+02 * invT + 1.661320882000000e+00 -
                 2.344331120000000e+00 * tc[0] - 3.990260375000000e-03 * tc[1] +
                 3.246358500000000e-06 * tc[2] - 1.679767450000000e-09 * tc[3] +
                 3.688058805000000e-13 * tc[4];
  } else {
    // species 0: CH4
    species[0] = -9.468344590000001e+03 * invT - 1.836246650500000e+01 -
                 7.485149500000000e-02 * tc[0] - 6.695473350000000e-03 * tc[1] +
                 9.554763483333333e-07 * tc[2] - 1.019104458333333e-10 * tc[3] +
                 5.090761500000000e-15 * tc[4];
    // species 1: O2
    species[1] = -1.088457720000000e+03 * invT - 2.170693450000000e+00 -
                 3.282537840000000e+00 * tc[0] - 7.415437700000000e-04 * tc[1] +
                 1.263277781666667e-07 * tc[2] - 1.745587958333333e-11 * tc[3] +
                 1.083588970000000e-15 * tc[4];
    // species 2: H2O
    species[2] = -3.000429710000000e+04 * invT - 1.932777610000000e+00 -
                 3.033992490000000e+00 * tc[0] - 1.088459020000000e-03 * tc[1] +
                 2.734541966666666e-08 * tc[2] + 8.086832250000000e-12 * tc[3] -
                 8.410049600000000e-16 * tc[4];
    // species 3: N2
    species[3] = -9.227977000000000e+02 * invT - 3.053888000000000e+00 -
                 2.926640000000000e+00 * tc[0] - 7.439884000000000e-04 * tc[1] +
                 9.474600000000001e-08 * tc[2] - 8.414198333333333e-12 * tc[3] +
                 3.376675500000000e-16 * tc[4];
    // species 4: CO
    species[4] = -1.415187240000000e+04 * invT - 5.103502110000000e+00 -
                 2.715185610000000e+00 * tc[0] - 1.031263715000000e-03 * tc[1] +
                 1.664709618333334e-07 * tc[2] - 1.917108400000000e-11 * tc[3] +
                 1.018238580000000e-15 * tc[4];
    // species 5: CO2
    species[5] = -4.875916600000000e+04 * invT + 1.585822230000000e+00 -
                 3.857460290000000e+00 * tc[0] - 2.207185130000000e-03 * tc[1] +
                 3.691356733333334e-07 * tc[2] - 4.362418233333334e-11 * tc[3] +
                 2.360420820000000e-15 * tc[4];
    // species 6: H2
    species[6] = -9.501589220000000e+02 * invT + 6.542302510000000e+00 -
                 3.337279200000000e+00 * tc[0] + 2.470123655000000e-05 * tc[1] -
                 8.324279633333333e-08 * tc[2] + 1.496386616666667e-11 * tc[3] -
                 1.001276880000000e-15 * tc[4];
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
    // species 0: CH4
    species[0] = -1.02466476e+04 * invT + 8.79117989e+00 -
                 5.14987613e+00 * tc[0] + 6.83548940e-03 * tc[1] -
                 8.19667665e-06 * tc[2] + 4.03952522e-09 * tc[3] -
                 8.33469780e-13 * tc[4];
    // species 1: O2
    species[1] = -1.06394356e+03 * invT - 8.75219370e-01 -
                 3.78245636e+00 * tc[0] + 1.49836708e-03 * tc[1] -
                 1.64121700e-06 * tc[2] + 8.06774591e-10 * tc[3] -
                 1.62186418e-13 * tc[4];
    // species 2: H2O
    species[2] = -3.02937267e+04 * invT + 4.04767277e+00 -
                 4.19864056e+00 * tc[0] + 1.01821705e-03 * tc[1] -
                 1.08673369e-06 * tc[2] + 4.57330885e-10 * tc[3] -
                 8.85989085e-14 * tc[4];
    // species 3: N2
    species[3] = -1.02089990e+03 * invT - 1.65169500e+00 -
                 3.29867700e+00 * tc[0] - 7.04120200e-04 * tc[1] +
                 6.60537000e-07 * tc[2] - 4.70126250e-10 * tc[3] +
                 1.22242700e-13 * tc[4];
    // species 4: CO
    species[4] = -1.43440860e+04 * invT - 9.28875810e-01 -
                 3.57953347e+00 * tc[0] + 3.05176840e-04 * tc[1] -
                 1.69469055e-07 * tc[2] - 7.55838237e-11 * tc[3] +
                 4.52212249e-14 * tc[4];
    // species 5: CO2
    species[5] = -4.83719697e+04 * invT - 8.54427870e+00 -
                 2.35677352e+00 * tc[0] - 4.49229839e-03 * tc[1] +
                 1.18726045e-06 * tc[2] - 2.04932518e-10 * tc[3] +
                 7.18497740e-15 * tc[4];
    // species 6: H2
    species[6] = -9.17935173e+02 * invT + 6.61320882e-01 -
                 2.34433112e+00 * tc[0] - 3.99026037e-03 * tc[1] +
                 3.24635850e-06 * tc[2] - 1.67976745e-09 * tc[3] +
                 3.68805881e-13 * tc[4];
  } else {
    // species 0: CH4
    species[0] = -9.46834459e+03 * invT - 1.93624665e+01 -
                 7.48514950e-02 * tc[0] - 6.69547335e-03 * tc[1] +
                 9.55476348e-07 * tc[2] - 1.01910446e-10 * tc[3] +
                 5.09076150e-15 * tc[4];
    // species 1: O2
    species[1] = -1.08845772e+03 * invT - 3.17069345e+00 -
                 3.28253784e+00 * tc[0] - 7.41543770e-04 * tc[1] +
                 1.26327778e-07 * tc[2] - 1.74558796e-11 * tc[3] +
                 1.08358897e-15 * tc[4];
    // species 2: H2O
    species[2] = -3.00042971e+04 * invT - 2.93277761e+00 -
                 3.03399249e+00 * tc[0] - 1.08845902e-03 * tc[1] +
                 2.73454197e-08 * tc[2] + 8.08683225e-12 * tc[3] -
                 8.41004960e-16 * tc[4];
    // species 3: N2
    species[3] = -9.22797700e+02 * invT - 4.05388800e+00 -
                 2.92664000e+00 * tc[0] - 7.43988400e-04 * tc[1] +
                 9.47460000e-08 * tc[2] - 8.41419833e-12 * tc[3] +
                 3.37667550e-16 * tc[4];
    // species 4: CO
    species[4] = -1.41518724e+04 * invT - 6.10350211e+00 -
                 2.71518561e+00 * tc[0] - 1.03126372e-03 * tc[1] +
                 1.66470962e-07 * tc[2] - 1.91710840e-11 * tc[3] +
                 1.01823858e-15 * tc[4];
    // species 5: CO2
    species[5] = -4.87591660e+04 * invT + 5.85822230e-01 -
                 3.85746029e+00 * tc[0] - 2.20718513e-03 * tc[1] +
                 3.69135673e-07 * tc[2] - 4.36241823e-11 * tc[3] +
                 2.36042082e-15 * tc[4];
    // species 6: H2
    species[6] = -9.50158922e+02 * invT + 5.54230251e+00 -
                 3.33727920e+00 * tc[0] + 2.47012365e-05 * tc[1] -
                 8.32427963e-08 * tc[2] + 1.49638662e-11 * tc[3] -
                 1.00127688e-15 * tc[4];
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
    // species 0: CH4
    species[0] = +4.14987613e+00 - 6.83548940e-03 * tc[1] +
                 1.63933533e-05 * tc[2] - 1.21185757e-08 * tc[3] +
                 3.33387912e-12 * tc[4] - 1.02466476e+04 * invT;
    // species 1: O2
    species[1] = +2.78245636e+00 - 1.49836708e-03 * tc[1] +
                 3.28243400e-06 * tc[2] - 2.42032377e-09 * tc[3] +
                 6.48745674e-13 * tc[4] - 1.06394356e+03 * invT;
    // species 2: H2O
    species[2] = +3.19864056e+00 - 1.01821705e-03 * tc[1] +
                 2.17346737e-06 * tc[2] - 1.37199266e-09 * tc[3] +
                 3.54395634e-13 * tc[4] - 3.02937267e+04 * invT;
    // species 3: N2
    species[3] = +2.29867700e+00 + 7.04120200e-04 * tc[1] -
                 1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
                 4.88970800e-13 * tc[4] - 1.02089990e+03 * invT;
    // species 4: CO
    species[4] = +2.57953347e+00 - 3.05176840e-04 * tc[1] +
                 3.38938110e-07 * tc[2] + 2.26751471e-10 * tc[3] -
                 1.80884900e-13 * tc[4] - 1.43440860e+04 * invT;
    // species 5: CO2
    species[5] = +1.35677352e+00 + 4.49229839e-03 * tc[1] -
                 2.37452090e-06 * tc[2] + 6.14797555e-10 * tc[3] -
                 2.87399096e-14 * tc[4] - 4.83719697e+04 * invT;
    // species 6: H2
    species[6] = +1.34433112e+00 + 3.99026037e-03 * tc[1] -
                 6.49271700e-06 * tc[2] + 5.03930235e-09 * tc[3] -
                 1.47522352e-12 * tc[4] - 9.17935173e+02 * invT;
  } else {
    // species 0: CH4
    species[0] = -9.25148505e-01 + 6.69547335e-03 * tc[1] -
                 1.91095270e-06 * tc[2] + 3.05731338e-10 * tc[3] -
                 2.03630460e-14 * tc[4] - 9.46834459e+03 * invT;
    // species 1: O2
    species[1] = +2.28253784e+00 + 7.41543770e-04 * tc[1] -
                 2.52655556e-07 * tc[2] + 5.23676387e-11 * tc[3] -
                 4.33435588e-15 * tc[4] - 1.08845772e+03 * invT;
    // species 2: H2O
    species[2] = +2.03399249e+00 + 1.08845902e-03 * tc[1] -
                 5.46908393e-08 * tc[2] - 2.42604967e-11 * tc[3] +
                 3.36401984e-15 * tc[4] - 3.00042971e+04 * invT;
    // species 3: N2
    species[3] = +1.92664000e+00 + 7.43988400e-04 * tc[1] -
                 1.89492000e-07 * tc[2] + 2.52425950e-11 * tc[3] -
                 1.35067020e-15 * tc[4] - 9.22797700e+02 * invT;
    // species 4: CO
    species[4] = +1.71518561e+00 + 1.03126372e-03 * tc[1] -
                 3.32941924e-07 * tc[2] + 5.75132520e-11 * tc[3] -
                 4.07295432e-15 * tc[4] - 1.41518724e+04 * invT;
    // species 5: CO2
    species[5] = +2.85746029e+00 + 2.20718513e-03 * tc[1] -
                 7.38271347e-07 * tc[2] + 1.30872547e-10 * tc[3] -
                 9.44168328e-15 * tc[4] - 4.87591660e+04 * invT;
    // species 6: H2
    species[6] = +2.33727920e+00 - 2.47012365e-05 * tc[1] +
                 1.66485593e-07 * tc[2] - 4.48915985e-11 * tc[3] +
                 4.00510752e-15 * tc[4] - 9.50158922e+02 * invT;
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
    // species 0: CH4
    species[0] = +5.14987613e+00 - 6.83548940e-03 * tc[1] +
                 1.63933533e-05 * tc[2] - 1.21185757e-08 * tc[3] +
                 3.33387912e-12 * tc[4] - 1.02466476e+04 * invT;
    // species 1: O2
    species[1] = +3.78245636e+00 - 1.49836708e-03 * tc[1] +
                 3.28243400e-06 * tc[2] - 2.42032377e-09 * tc[3] +
                 6.48745674e-13 * tc[4] - 1.06394356e+03 * invT;
    // species 2: H2O
    species[2] = +4.19864056e+00 - 1.01821705e-03 * tc[1] +
                 2.17346737e-06 * tc[2] - 1.37199266e-09 * tc[3] +
                 3.54395634e-13 * tc[4] - 3.02937267e+04 * invT;
    // species 3: N2
    species[3] = +3.29867700e+00 + 7.04120200e-04 * tc[1] -
                 1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
                 4.88970800e-13 * tc[4] - 1.02089990e+03 * invT;
    // species 4: CO
    species[4] = +3.57953347e+00 - 3.05176840e-04 * tc[1] +
                 3.38938110e-07 * tc[2] + 2.26751471e-10 * tc[3] -
                 1.80884900e-13 * tc[4] - 1.43440860e+04 * invT;
    // species 5: CO2
    species[5] = +2.35677352e+00 + 4.49229839e-03 * tc[1] -
                 2.37452090e-06 * tc[2] + 6.14797555e-10 * tc[3] -
                 2.87399096e-14 * tc[4] - 4.83719697e+04 * invT;
    // species 6: H2
    species[6] = +2.34433112e+00 + 3.99026037e-03 * tc[1] -
                 6.49271700e-06 * tc[2] + 5.03930235e-09 * tc[3] -
                 1.47522352e-12 * tc[4] - 9.17935173e+02 * invT;
  } else {
    // species 0: CH4
    species[0] = +7.48514950e-02 + 6.69547335e-03 * tc[1] -
                 1.91095270e-06 * tc[2] + 3.05731338e-10 * tc[3] -
                 2.03630460e-14 * tc[4] - 9.46834459e+03 * invT;
    // species 1: O2
    species[1] = +3.28253784e+00 + 7.41543770e-04 * tc[1] -
                 2.52655556e-07 * tc[2] + 5.23676387e-11 * tc[3] -
                 4.33435588e-15 * tc[4] - 1.08845772e+03 * invT;
    // species 2: H2O
    species[2] = +3.03399249e+00 + 1.08845902e-03 * tc[1] -
                 5.46908393e-08 * tc[2] - 2.42604967e-11 * tc[3] +
                 3.36401984e-15 * tc[4] - 3.00042971e+04 * invT;
    // species 3: N2
    species[3] = +2.92664000e+00 + 7.43988400e-04 * tc[1] -
                 1.89492000e-07 * tc[2] + 2.52425950e-11 * tc[3] -
                 1.35067020e-15 * tc[4] - 9.22797700e+02 * invT;
    // species 4: CO
    species[4] = +2.71518561e+00 + 1.03126372e-03 * tc[1] -
                 3.32941924e-07 * tc[2] + 5.75132520e-11 * tc[3] -
                 4.07295432e-15 * tc[4] - 1.41518724e+04 * invT;
    // species 5: CO2
    species[5] = +3.85746029e+00 + 2.20718513e-03 * tc[1] -
                 7.38271347e-07 * tc[2] + 1.30872547e-10 * tc[3] -
                 9.44168328e-15 * tc[4] - 4.87591660e+04 * invT;
    // species 6: H2
    species[6] = +3.33727920e+00 - 2.47012365e-05 * tc[1] +
                 1.66485593e-07 * tc[2] - 4.48915985e-11 * tc[3] +
                 4.00510752e-15 * tc[4] - 9.50158922e+02 * invT;
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
    // species 0: CH4
    species[0] = +5.14987613e+00 * tc[0] - 1.36709788e-02 * tc[1] +
                 2.45900299e-05 * tc[2] - 1.61581009e-08 * tc[3] +
                 4.16734890e-12 * tc[4] - 4.64130376e+00;
    // species 1: O2
    species[1] = +3.78245636e+00 * tc[0] - 2.99673416e-03 * tc[1] +
                 4.92365101e-06 * tc[2] - 3.22709836e-09 * tc[3] +
                 8.10932092e-13 * tc[4] + 3.65767573e+00;
    // species 2: H2O
    species[2] = +4.19864056e+00 * tc[0] - 2.03643410e-03 * tc[1] +
                 3.26020105e-06 * tc[2] - 1.82932354e-09 * tc[3] +
                 4.42994543e-13 * tc[4] - 8.49032208e-01;
    // species 3: N2
    species[3] = +3.29867700e+00 * tc[0] + 1.40824040e-03 * tc[1] -
                 1.98161100e-06 * tc[2] + 1.88050500e-09 * tc[3] -
                 6.11213500e-13 * tc[4] + 3.95037200e+00;
    // species 4: CO
    species[4] = +3.57953347e+00 * tc[0] - 6.10353680e-04 * tc[1] +
                 5.08407165e-07 * tc[2] + 3.02335295e-10 * tc[3] -
                 2.26106125e-13 * tc[4] + 3.50840928e+00;
    // species 5: CO2
    species[5] = +2.35677352e+00 * tc[0] + 8.98459677e-03 * tc[1] -
                 3.56178134e-06 * tc[2] + 8.19730073e-10 * tc[3] -
                 3.59248870e-14 * tc[4] + 9.90105222e+00;
    // species 6: H2
    species[6] = +2.34433112e+00 * tc[0] + 7.98052075e-03 * tc[1] -
                 9.73907550e-06 * tc[2] + 6.71906980e-09 * tc[3] -
                 1.84402940e-12 * tc[4] + 6.83010238e-01;
  } else {
    // species 0: CH4
    species[0] = +7.48514950e-02 * tc[0] + 1.33909467e-02 * tc[1] -
                 2.86642905e-06 * tc[2] + 4.07641783e-10 * tc[3] -
                 2.54538075e-14 * tc[4] + 1.84373180e+01;
    // species 1: O2
    species[1] = +3.28253784e+00 * tc[0] + 1.48308754e-03 * tc[1] -
                 3.78983334e-07 * tc[2] + 6.98235183e-11 * tc[3] -
                 5.41794485e-15 * tc[4] + 5.45323129e+00;
    // species 2: H2O
    species[2] = +3.03399249e+00 * tc[0] + 2.17691804e-03 * tc[1] -
                 8.20362590e-08 * tc[2] - 3.23473290e-11 * tc[3] +
                 4.20502480e-15 * tc[4] + 4.96677010e+00;
    // species 3: N2
    species[3] = +2.92664000e+00 * tc[0] + 1.48797680e-03 * tc[1] -
                 2.84238000e-07 * tc[2] + 3.36567933e-11 * tc[3] -
                 1.68833775e-15 * tc[4] + 5.98052800e+00;
    // species 4: CO
    species[4] = +2.71518561e+00 * tc[0] + 2.06252743e-03 * tc[1] -
                 4.99412886e-07 * tc[2] + 7.66843360e-11 * tc[3] -
                 5.09119290e-15 * tc[4] + 7.81868772e+00;
    // species 5: CO2
    species[5] = +3.85746029e+00 * tc[0] + 4.41437026e-03 * tc[1] -
                 1.10740702e-06 * tc[2] + 1.74496729e-10 * tc[3] -
                 1.18021041e-14 * tc[4] + 2.27163806e+00;
    // species 6: H2
    species[6] = +3.33727920e+00 * tc[0] - 4.94024731e-05 * tc[1] +
                 2.49728389e-07 * tc[2] - 5.98554647e-11 * tc[3] +
                 5.00638440e-15 * tc[4] - 3.20502331e+00;
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
  amrex::Real cpor[7];                                //  temporary storage
  cp_R(cpor, tc);

  // perform dot product
  for (id = 0; id < 7; ++id) {
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
  amrex::Real cpor[7], tresult[7];                    //  temporary storage
  amrex::Real imw[7];

  get_imw(imw);

  cp_R(cpor, tc);

  for (int i = 0; i < 7; i++) {
    tresult[i] = cpor[i] * y[i] * imw[i];
  }
  for (int i = 0; i < 7; i++) {
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
  amrex::Real cvor[7];                                //  temporary storage
  cv_R(cvor, tc);

  // perform dot product
  for (id = 0; id < 7; ++id) {
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
  amrex::Real cvor[7];                                //  temporary storage
  amrex::Real imw[7];

  get_imw(imw);

  cv_R(cvor, tc);

  // multiply by y/molecularweight
  result += cvor[0] * y[0] * imw[0]; // CH4
  result += cvor[1] * y[1] * imw[1]; // O2
  result += cvor[2] * y[2] * imw[2]; // H2O
  result += cvor[3] * y[3] * imw[3]; // N2
  result += cvor[4] * y[4] * imw[4]; // CO
  result += cvor[5] * y[5] * imw[5]; // CO2
  result += cvor[6] * y[6] * imw[6]; // H2

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
  amrex::Real hml[7];                                 //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesEnthalpy(hml, tc);

  // perform dot product
  for (id = 0; id < 7; ++id) {
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
  amrex::Real hml[7], tmp[7];                         //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[7];

  get_imw(imw);

  speciesEnthalpy(hml, tc);

  int id;
  for (id = 0; id < 7; ++id) {
    tmp[id] = y[id] * hml[id] * imw[id];
  }
  for (id = 0; id < 7; ++id) {
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
  amrex::Real uml[7];                                 //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesInternalEnergy(uml, tc);

  // perform dot product
  for (id = 0; id < 7; ++id) {
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
  amrex::Real ums[7];                                 //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[7];

  get_imw(imw);

  speciesInternalEnergy(ums, tc);

  // perform dot product + scaling by wt
  result += y[0] * ums[0] * imw[0]; // CH4
  result += y[1] * ums[1] * imw[1]; // O2
  result += y[2] * ums[2] * imw[2]; // H2O
  result += y[3] * ums[3] * imw[3]; // N2
  result += y[4] * ums[4] * imw[4]; // CO
  result += y[5] * ums[5] * imw[5]; // CO2
  result += y[6] * ums[6] * imw[6]; // H2

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
  amrex::Real sor[7]; //  temporary storage
  speciesEntropy(sor, tc);

  // Compute Eq 42
  for (id = 0; id < 7; ++id) {
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
  amrex::Real sor[7];  //  temporary storage
  amrex::Real x[7];    //  need a ytx conversion
  amrex::Real YOW = 0; // See Eq 4, 6 in CK Manual
  amrex::Real imw[7];

  get_imw(imw);

  // Compute inverse of mean molecular wt first
  YOW += y[0] * imw[0]; // CH4
  YOW += y[1] * imw[1]; // O2
  YOW += y[2] * imw[2]; // H2O
  YOW += y[3] * imw[3]; // N2
  YOW += y[4] * imw[4]; // CO
  YOW += y[5] * imw[5]; // CO2
  YOW += y[6] * imw[6]; // H2
  // Now compute y to x conversion
  x[0] = y[0] / (16.043000 * YOW);
  x[1] = y[1] / (31.998000 * YOW);
  x[2] = y[2] / (18.015000 * YOW);
  x[3] = y[3] / (28.014000 * YOW);
  x[4] = y[4] / (28.010000 * YOW);
  x[5] = y[5] / (44.009000 * YOW);
  x[6] = y[6] / (2.016000 * YOW);
  speciesEntropy(sor, tc);
  // Perform computation in Eq 42 and 43
  result += x[0] * (sor[0] - log((x[0] + 1e-100)) - logPratio);
  result += x[1] * (sor[1] - log((x[1] + 1e-100)) - logPratio);
  result += x[2] * (sor[2] - log((x[2] + 1e-100)) - logPratio);
  result += x[3] * (sor[3] - log((x[3] + 1e-100)) - logPratio);
  result += x[4] * (sor[4] - log((x[4] + 1e-100)) - logPratio);
  result += x[5] * (sor[5] - log((x[5] + 1e-100)) - logPratio);
  result += x[6] * (sor[6] - log((x[6] + 1e-100)) - logPratio);
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
  XW += x[0] * 16.043000;                       // CH4
  XW += x[1] * 31.998000;                       // O2
  XW += x[2] * 18.015000;                       // H2O
  XW += x[3] * 28.014000;                       // N2
  XW += x[4] * 28.010000;                       // CO
  XW += x[5] * 44.009000;                       // CO2
  XW += x[6] * 2.016000;                        // H2
  *P = *rho * 8.31446261815324e+07 * (*T) / XW; // P = rho*R*T/W

  return;
}

// Compute P = rhoRT/W(y)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPY(amrex::Real* rho, amrex::Real* T, amrex::Real* y, amrex::Real* P)
{
  amrex::Real YOW = 0; //  for computing mean MW
  amrex::Real imw[7];

  get_imw(imw);

  YOW += y[0] * imw[0];                          // CH4
  YOW += y[1] * imw[1];                          // O2
  YOW += y[2] * imw[2];                          // H2O
  YOW += y[3] * imw[3];                          // N2
  YOW += y[4] * imw[4];                          // CO
  YOW += y[5] * imw[5];                          // CO2
  YOW += y[6] * imw[6];                          // H2
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
  W += c[0] * 16.043000; // CH4
  W += c[1] * 31.998000; // O2
  W += c[2] * 18.015000; // H2O
  W += c[3] * 28.014000; // N2
  W += c[4] * 28.010000; // CO
  W += c[5] * 44.009000; // CO2
  W += c[6] * 2.016000;  // H2

  for (id = 0; id < 7; ++id) {
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
  XW += x[0] * 16.043000;                         // CH4
  XW += x[1] * 31.998000;                         // O2
  XW += x[2] * 18.015000;                         // H2O
  XW += x[3] * 28.014000;                         // N2
  XW += x[4] * 28.010000;                         // CO
  XW += x[5] * 44.009000;                         // CO2
  XW += x[6] * 2.016000;                          // H2
  *rho = *P * XW / (8.31446261815324e+07 * (*T)); // rho = P*W/(R*T)

  return;
}

// Compute rho = P*W(y)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOY(amrex::Real* P, amrex::Real* T, amrex::Real* y, amrex::Real* rho)
{
  amrex::Real YOW = 0;
  amrex::Real tmp[7];
  amrex::Real imw[7];

  get_imw(imw);

  for (int i = 0; i < 7; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 7; i++) {
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
  W += c[0] * 16.043000; // CH4
  W += c[1] * 31.998000; // O2
  W += c[2] * 18.015000; // H2O
  W += c[3] * 28.014000; // N2
  W += c[4] * 28.010000; // CO
  W += c[5] * 44.009000; // CO2
  W += c[6] * 2.016000;  // H2

  for (id = 0; id < 7; ++id) {
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
  amrex::Real tmp[7];
  amrex::Real imw[7];

  get_imw(imw);

  for (int i = 0; i < 7; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 7; i++) {
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
  XW += x[0] * 16.043000; // CH4
  XW += x[1] * 31.998000; // O2
  XW += x[2] * 18.015000; // H2O
  XW += x[3] * 28.014000; // N2
  XW += x[4] * 28.010000; // CO
  XW += x[5] * 44.009000; // CO2
  XW += x[6] * 2.016000;  // H2
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
  W += c[0] * 16.043000; // CH4
  W += c[1] * 31.998000; // O2
  W += c[2] * 18.015000; // H2O
  W += c[3] * 28.014000; // N2
  W += c[4] * 28.010000; // CO
  W += c[5] * 44.009000; // CO2
  W += c[6] * 2.016000;  // H2

  for (id = 0; id < 7; ++id) {
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
  amrex::Real tmp[7];
  amrex::Real imw[7];

  get_imw(imw);

  for (int i = 0; i < 7; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 7; i++) {
    YOW += tmp[i];
  }

  amrex::Real YOWINV = 1.0 / YOW;

  for (int i = 0; i < 7; i++) {
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
  amrex::Real imw[7];

  get_imw(imw);

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 7; i++) {
    c[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 7; i++) {
    YOW += c[i];
  }

  // PW/RT (see Eq. 7)
  PWORT = (*P) / (YOW * 8.31446261815324e+07 * (*T));
  // Now compute conversion

  for (int i = 0; i < 7; i++) {
    c[i] = PWORT * y[i] * imw[i];
  }
  return;
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTCR(amrex::Real* rho, amrex::Real* /*T*/, amrex::Real* y, amrex::Real* c)
{
  amrex::Real imw[7];

  get_imw(imw);

  for (int i = 0; i < 7; i++) {
    c[i] = (*rho) * y[i] * imw[i];
  }
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTY(amrex::Real* x, amrex::Real* y)
{
  amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
  // Compute mean molecular wt first
  XW += x[0] * 16.043000; // CH4
  XW += x[1] * 31.998000; // O2
  XW += x[2] * 18.015000; // H2O
  XW += x[3] * 28.014000; // N2
  XW += x[4] * 28.010000; // CO
  XW += x[5] * 44.009000; // CO2
  XW += x[6] * 2.016000;  // H2
  // Now compute conversion
  amrex::Real XWinv = 1.0 / XW;
  y[0] = x[0] * 16.043000 * XWinv;
  y[1] = x[1] * 31.998000 * XWinv;
  y[2] = x[2] * 18.015000 * XWinv;
  y[3] = x[3] * 28.014000 * XWinv;
  y[4] = x[4] * 28.010000 * XWinv;
  y[5] = x[5] * 44.009000 * XWinv;
  y[6] = x[6] * 2.016000 * XWinv;

  return;
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTCP(amrex::Real* P, amrex::Real* T, amrex::Real* x, amrex::Real* c)
{
  int id;                                                  // loop counter
  amrex::Real PORT = (*P) / (8.31446261815324e+07 * (*T)); // P/RT

  // Compute conversion, see Eq 10
  for (id = 0; id < 7; ++id) {
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
  XW += x[0] * 16.043000; // CH4
  XW += x[1] * 31.998000; // O2
  XW += x[2] * 18.015000; // H2O
  XW += x[3] * 28.014000; // N2
  XW += x[4] * 28.010000; // CO
  XW += x[5] * 44.009000; // CO2
  XW += x[6] * 2.016000;  // H2
  ROW = (*rho) / XW;

  // Compute conversion, see Eq 11
  for (id = 0; id < 7; ++id) {
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
  for (id = 0; id < 7; ++id) {
    sumC += c[id];
  }

  //  See Eq 13
  amrex::Real sumCinv = 1.0 / sumC;
  for (id = 0; id < 7; ++id) {
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
  CW += c[0] * 16.043000; // CH4
  CW += c[1] * 31.998000; // O2
  CW += c[2] * 18.015000; // H2O
  CW += c[3] * 28.014000; // N2
  CW += c[4] * 28.010000; // CO
  CW += c[5] * 44.009000; // CO2
  CW += c[6] * 2.016000;  // H2
  // Now compute conversion
  amrex::Real CWinv = 1.0 / CW;
  y[0] = c[0] * 16.043000 * CWinv;
  y[1] = c[1] * 31.998000 * CWinv;
  y[2] = c[2] * 18.015000 * CWinv;
  y[3] = c[3] * 28.014000 * CWinv;
  y[4] = c[4] * 28.010000 * CWinv;
  y[5] = c[5] * 44.009000 * CWinv;
  y[6] = c[6] * 2.016000 * CWinv;

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
  for (id = 0; id < 7; ++id) {
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
  for (id = 0; id < 7; ++id) {
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
  for (id = 0; id < 7; ++id) {
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
  for (id = 0; id < 7; ++id) {
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
  for (id = 0; id < 7; ++id) {
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
  cvms[0] *= 5.182610869633635e+06; // CH4
  cvms[1] *= 2.598431970171023e+06; // O2
  cvms[2] *= 4.615299815794193e+06; // H2O
  cvms[3] *= 2.967966951578939e+06; // N2
  cvms[4] *= 2.968390795484913e+06; // CO
  cvms[5] *= 1.889264154639560e+06; // CO2
  cvms[6] *= 4.124237409798234e+07; // H2
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
  cpms[0] *= 5.182610869633635e+06; // CH4
  cpms[1] *= 2.598431970171023e+06; // O2
  cpms[2] *= 4.615299815794193e+06; // H2O
  cpms[3] *= 2.967966951578939e+06; // N2
  cpms[4] *= 2.968390795484913e+06; // CO
  cpms[5] *= 1.889264154639560e+06; // CO2
  cpms[6] *= 4.124237409798234e+07; // H2
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUMS(amrex::Real* T, amrex::Real* ums)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[7];

  get_imw(imw);

  speciesInternalEnergy(ums, tc);

  for (int i = 0; i < 7; i++) {
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
  amrex::Real imw[7];

  get_imw(imw);

  speciesEnthalpy(hms, tc);

  for (int i = 0; i < 7; i++) {
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
  sms[0] *= 5.182610869633635e+06; // CH4
  sms[1] *= 2.598431970171023e+06; // O2
  sms[2] *= 4.615299815794193e+06; // H2O
  sms[3] *= 2.967966951578939e+06; // N2
  sms[4] *= 2.968390795484913e+06; // CO
  sms[5] *= 1.889264154639560e+06; // CO2
  sms[6] *= 4.124237409798234e+07; // H2
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

  // reaction 0: 2 CH4 + O2 => 2 CO + 4 H2
  qf[0] = pow(sc[0], 0.500000) * pow(sc[1], 1.250000);
  qr[0] = 0.0;

  // reaction 1: CH4 + H2O <=> CO + 3 H2
  qf[1] = sc[0] * sc[2];
  qr[1] = sc[4] * pow(sc[6], 3.000000);

  // reaction 2: 2 H2 + O2 => 2 H2O
  qf[2] = pow(sc[1], 1.500000) * pow(sc[6], 0.250000);
  qr[2] = 0.0;

  // reaction 3: CO + H2O <=> CO2 + H2
  qf[3] = sc[2] * sc[4];
  qr[3] = sc[5] * sc[6];

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 7; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[7];
  gibbs(g_RT, tc);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 * invT;
  amrex::Real refCinv = 1 / refC;

  // Evaluate the kfs
  amrex::Real k_f, k_r, Corr;

  // reaction 0:  2 CH4 + O2 => 2 CO + 4 H2
  k_f = 219875.458149426 * exp(-(15096.4999741416) * invT);
  qf[0] *= k_f;
  qr[0] *=
    k_f *
    exp(-(
      2.000000 * g_RT[0] + g_RT[1] - 2.000000 * g_RT[4] - 4.000000 * g_RT[6])) *
    (pow(refCinv, 3.000000));
  // reaction 1:  CH4 + H2O <=> CO + 3 H2
  k_f = 300000 * exp(-(15096.4999741416) * invT);
  qf[1] *= k_f;
  qr[1] *= k_f * exp(-(g_RT[0] + g_RT[2] - g_RT[4] - 3.000000 * g_RT[6])) *
           ((refCinv * refCinv));
  // reaction 2:  2 H2 + O2 => 2 H2O
  k_f = 33993533107.7566 * exp((-1) * tc[0] - (20128.6666321888) * invT);
  qf[2] *= k_f;
  qr[2] *=
    k_f * exp(-(g_RT[1] - 2.000000 * g_RT[2] + 2.000000 * g_RT[6])) * (refC);
  // reaction 3:  CO + H2O <=> CO2 + H2
  k_f = 2750000 * exp(-(10064.3333160944) * invT);
  qf[3] *= k_f;
  qr[3] *= k_f * exp(-(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[6]));

  return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
productionRate(amrex::Real* wdot, amrex::Real* sc, amrex::Real T)
{
  const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
  const amrex::Real invT = 1.0 / tc[1];

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  const amrex::Real refC = 101325 / 8.31446 * invT;
  const amrex::Real refCinv = 1 / refC;

  for (int i = 0; i < 7; ++i) {
    wdot[i] = 0.0;
  }

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 7; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[7];
  gibbs(g_RT, tc);

  {
    // reaction 0:  2 CH4 + O2 => 2 CO + 4 H2
    const amrex::Real k_f = 219875.458149426 * exp(-(15096.4999741416) * invT);
    const amrex::Real qf = k_f * (pow(sc[0], 0.500000) * pow(sc[1], 1.250000));
    const amrex::Real qr = k_f *
                           exp(
                             -(2.000000 * g_RT[0] + g_RT[1] -
                               2.000000 * g_RT[4] - 4.000000 * g_RT[6])) *
                           (pow(refCinv, 3.000000)) * (0.0);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= 2.000000 * qdot;
    wdot[1] -= qdot;
    wdot[4] += 2.000000 * qdot;
    wdot[6] += 4.000000 * qdot;
  }

  {
    // reaction 1:  CH4 + H2O <=> CO + 3 H2
    const amrex::Real k_f = 300000 * exp(-(15096.4999741416) * invT);
    const amrex::Real qf = k_f * (sc[0] * sc[2]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] + g_RT[2] - g_RT[4] - 3.000000 * g_RT[6])) *
      ((refCinv * refCinv)) * (sc[4] * pow(sc[6], 3.000000));
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[6] += 3.000000 * qdot;
  }

  {
    // reaction 2:  2 H2 + O2 => 2 H2O
    const amrex::Real k_f =
      33993533107.7566 * exp((-1) * tc[0] - (20128.6666321888) * invT);
    const amrex::Real qf = k_f * (pow(sc[1], 1.500000) * pow(sc[6], 0.250000));
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - 2.000000 * g_RT[2] + 2.000000 * g_RT[6])) * (refC) *
      (0.0);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[2] += 2.000000 * qdot;
    wdot[6] -= 2.000000 * qdot;
  }

  {
    // reaction 3:  CO + H2O <=> CO2 + H2
    const amrex::Real k_f = 2750000 * exp(-(10064.3333160944) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[4]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[6])) * (sc[5] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
  }

  return;
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(amrex::Real* T, amrex::Real* C, amrex::Real* wdot)
{
  int id; // loop counter

  // convert to SI
  for (id = 0; id < 7; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, *T);

  // convert to chemkin units
  for (id = 0; id < 7; ++id) {
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
  amrex::Real c[7]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;
  amrex::Real imw[7];

  get_imw(imw);

  // Compute inverse of mean molecular wt first
  YOW += y[0] * imw[0]; // CH4
  YOW += y[1] * imw[1]; // O2
  YOW += y[2] * imw[2]; // H2O
  YOW += y[3] * imw[3]; // N2
  YOW += y[4] * imw[4]; // CO
  YOW += y[5] * imw[5]; // CO2
  YOW += y[6] * imw[6]; // H2
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

  // convert to chemkin units
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 7; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWXP(amrex::Real* P, amrex::Real* T, amrex::Real* x, amrex::Real* wdot)
{
  int id;           // loop counter
  amrex::Real c[7]; // temporary storage
  amrex::Real PORT =
    1e6 * (*P) /
    (8.31446261815324e+07 * (*T)); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (id = 0; id < 7; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 7; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWYR(amrex::Real* rho, amrex::Real* T, amrex::Real* y, amrex::Real* wdot)
{
  int id;           // loop counter
  amrex::Real c[7]; // temporary storage
  amrex::Real imw[7];

  get_imw(imw);

  // See Eq 8 with an extra 1e6 so c goes to SI
  c[0] = 1e6 * (*rho) * y[0] * imw[0];
  c[1] = 1e6 * (*rho) * y[1] * imw[1];
  c[2] = 1e6 * (*rho) * y[2] * imw[2];
  c[3] = 1e6 * (*rho) * y[3] * imw[3];
  c[4] = 1e6 * (*rho) * y[4] * imw[4];
  c[5] = 1e6 * (*rho) * y[5] * imw[5];
  c[6] = 1e6 * (*rho) * y[6] * imw[6];

  // call productionRate
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 7; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWXR(amrex::Real* rho, amrex::Real* T, amrex::Real* x, amrex::Real* wdot)
{
  int id;             // loop counter
  amrex::Real c[7];   // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 16.043000; // CH4
  XW += x[1] * 31.998000; // O2
  XW += x[2] * 18.015000; // H2O
  XW += x[3] * 28.014000; // N2
  XW += x[4] * 28.010000; // CO
  XW += x[5] * 44.009000; // CO2
  XW += x[6] * 2.016000;  // H2
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * (*rho) / XW;

  // Compute conversion, see Eq 11
  for (id = 0; id < 7; ++id) {
    c[id] = x[id] * ROW;
  }

  // convert to chemkin units
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 7; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int* kcharge)
{
  kcharge[0] = 0; // CH4
  kcharge[1] = 0; // O2
  kcharge[2] = 0; // H2O
  kcharge[3] = 0; // N2
  kcharge[4] = 0; // CO
  kcharge[5] = 0; // CO2
  kcharge[6] = 0; // H2
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real* zk)
{
  amrex::Real imw[7];
  get_imw(imw);

  int kchrg[7];
  CKCHRG(kchrg);

  for (int id = 0; id < 7; ++id) {
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
    // species 0: CH4
    species[0] = -1.36709788e-02 + 9.83601198e-05 * tc[1] -
                 1.45422908e-07 * tc[2] + 6.66775824e-11 * tc[3];
    // species 1: O2
    species[1] = -2.99673416e-03 + 1.96946040e-05 * tc[1] -
                 2.90438853e-08 * tc[2] + 1.29749135e-11 * tc[3];
    // species 2: H2O
    species[2] = -2.03643410e-03 + 1.30408042e-05 * tc[1] -
                 1.64639119e-08 * tc[2] + 7.08791268e-12 * tc[3];
    // species 3: N2
    species[3] = +1.40824040e-03 - 7.92644400e-06 * tc[1] +
                 1.69245450e-08 * tc[2] - 9.77941600e-12 * tc[3];
    // species 4: CO
    species[4] = -6.10353680e-04 + 2.03362866e-06 * tc[1] +
                 2.72101765e-09 * tc[2] - 3.61769800e-12 * tc[3];
    // species 5: CO2
    species[5] = +8.98459677e-03 - 1.42471254e-05 * tc[1] +
                 7.37757066e-09 * tc[2] - 5.74798192e-13 * tc[3];
    // species 6: H2
    species[6] = +7.98052075e-03 - 3.89563020e-05 * tc[1] +
                 6.04716282e-08 * tc[2] - 2.95044704e-11 * tc[3];
  } else {
    // species 0: CH4
    species[0] = +1.33909467e-02 - 1.14657162e-05 * tc[1] +
                 3.66877605e-09 * tc[2] - 4.07260920e-13 * tc[3];
    // species 1: O2
    species[1] = +1.48308754e-03 - 1.51593334e-06 * tc[1] +
                 6.28411665e-10 * tc[2] - 8.66871176e-14 * tc[3];
    // species 2: H2O
    species[2] = +2.17691804e-03 - 3.28145036e-07 * tc[1] -
                 2.91125961e-10 * tc[2] + 6.72803968e-14 * tc[3];
    // species 3: N2
    species[3] = +1.48797680e-03 - 1.13695200e-06 * tc[1] +
                 3.02911140e-10 * tc[2] - 2.70134040e-14 * tc[3];
    // species 4: CO
    species[4] = +2.06252743e-03 - 1.99765154e-06 * tc[1] +
                 6.90159024e-10 * tc[2] - 8.14590864e-14 * tc[3];
    // species 5: CO2
    species[5] = +4.41437026e-03 - 4.42962808e-06 * tc[1] +
                 1.57047056e-09 * tc[2] - 1.88833666e-13 * tc[3];
    // species 6: H2
    species[6] = -4.94024731e-05 + 9.98913556e-07 * tc[1] -
                 5.38699182e-10 * tc[2] + 8.01021504e-14 * tc[3];
  }
  return;
}

// compute an approx to the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian_precond(amrex::Real* J, amrex::Real* sc, amrex::Real T, const int HP)
{

  for (int i = 0; i < 64; i++) {
    J[i] = 0.0;
  }

  amrex::Real wdot[7];
  for (int k = 0; k < 7; k++) {
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
  for (int k = 0; k < 7; ++k) {
    mixture += sc[k];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[7];
  gibbs(g_RT, tc);

  // compute the species enthalpy
  amrex::Real h_RT[7];
  speciesEnthalpy(h_RT, tc);

  amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  amrex::Real dqdci, dcdc_fac, dqdc[7];
  amrex::Real Pr, fPr, F, k_0, logPr;
  amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
  amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
  amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const amrex::Real ln10 = log(10.0);
  const amrex::Real log10e = 1.0 / log(10.0);
  // reaction 0: 2 CH4 + O2 => 2 CO + 4 H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[0] * sc[0]) * sc[1];
  k_f = 219875.458149426 * exp(0 * tc[0] - 0.50321666580472 * (30000) * invT);
  dlnkfdT = 0 * invT + 0.50321666580472 * (30000) * invT2;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[0] -= 2 * q; // CH4
  wdot[1] -= q;     // O2
  wdot[4] += 2 * q; // CO
  wdot[6] += 4 * q; // H2
  // d()/d[CH4]
  dqdci = +k_f * 2.000000 * sc[0] * sc[1];
  J[0] += -2 * dqdci; // dwdot[CH4]/d[CH4]
  J[1] -= dqdci;      // dwdot[O2]/d[CH4]
  J[4] += 2 * dqdci;  // dwdot[CO]/d[CH4]
  J[6] += 4 * dqdci;  // dwdot[H2]/d[CH4]
  // d()/d[O2]
  dqdci = +k_f * pow(sc[0], 2.000000);
  J[8] += -2 * dqdci; // dwdot[CH4]/d[O2]
  J[9] -= dqdci;      // dwdot[O2]/d[O2]
  J[12] += 2 * dqdci; // dwdot[CO]/d[O2]
  J[14] += 4 * dqdci; // dwdot[H2]/d[O2]
  // d()/dT
  J[56] += -2 * dqdT; // dwdot[CH4]/dT
  J[57] -= dqdT;      // dwdot[O2]/dT
  J[60] += 2 * dqdT;  // dwdot[CO]/dT
  J[62] += 4 * dqdT;  // dwdot[H2]/dT

  // reaction 1: CH4 + H2O <=> CO + 3 H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[2];
  k_f = 300000 * exp(0 * tc[0] - 0.50321666580472 * (30000) * invT);
  dlnkfdT = 0 * invT + 0.50321666580472 * (30000) * invT2;
  // reverse
  phi_r = sc[4] * pow(sc[6], 3.000000);
  Kc = (refC * refC) * exp(g_RT[0] + g_RT[2] - g_RT[4] - 3.000000 * g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT =
    invT * (-(h_RT[0] + h_RT[2]) + (h_RT[4] + 3.000000 * h_RT[6]) - 2.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] -= q;     // CH4
  wdot[2] -= q;     // H2O
  wdot[4] += q;     // CO
  wdot[6] += 3 * q; // H2
  // d()/d[CH4]
  dqdci = +k_f * sc[2];
  J[0] -= dqdci;     // dwdot[CH4]/d[CH4]
  J[2] -= dqdci;     // dwdot[H2O]/d[CH4]
  J[4] += dqdci;     // dwdot[CO]/d[CH4]
  J[6] += 3 * dqdci; // dwdot[H2]/d[CH4]
  // d()/d[H2O]
  dqdci = +k_f * sc[0];
  J[16] -= dqdci;     // dwdot[CH4]/d[H2O]
  J[18] -= dqdci;     // dwdot[H2O]/d[H2O]
  J[20] += dqdci;     // dwdot[CO]/d[H2O]
  J[22] += 3 * dqdci; // dwdot[H2]/d[H2O]
  // d()/d[CO]
  dqdci = -k_r * pow(sc[6], 3.000000);
  J[32] -= dqdci;     // dwdot[CH4]/d[CO]
  J[34] -= dqdci;     // dwdot[H2O]/d[CO]
  J[36] += dqdci;     // dwdot[CO]/d[CO]
  J[38] += 3 * dqdci; // dwdot[H2]/d[CO]
  // d()/d[H2]
  dqdci = -k_r * sc[4] * 3.000000 * pow(sc[6], 2.000000);
  J[48] -= dqdci;     // dwdot[CH4]/d[H2]
  J[50] -= dqdci;     // dwdot[H2O]/d[H2]
  J[52] += dqdci;     // dwdot[CO]/d[H2]
  J[54] += 3 * dqdci; // dwdot[H2]/d[H2]
  // d()/dT
  J[56] -= dqdT;     // dwdot[CH4]/dT
  J[58] -= dqdT;     // dwdot[H2O]/dT
  J[60] += dqdT;     // dwdot[CO]/dT
  J[62] += 3 * dqdT; // dwdot[H2]/dT

  // reaction 2: 2 H2 + O2 => 2 H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * (sc[6] * sc[6]);
  k_f = 33993533107.7566 * exp(-1 * tc[0] - 0.50321666580472 * (40000) * invT);
  dlnkfdT = -1 * invT + 0.50321666580472 * (40000) * invT2;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[1] -= q;     // O2
  wdot[2] += 2 * q; // H2O
  wdot[6] -= 2 * q; // H2
  // d()/d[O2]
  dqdci = +k_f * pow(sc[6], 2.000000);
  J[9] -= dqdci;       // dwdot[O2]/d[O2]
  J[10] += 2 * dqdci;  // dwdot[H2O]/d[O2]
  J[14] += -2 * dqdci; // dwdot[H2]/d[O2]
  // d()/d[H2]
  dqdci = +k_f * sc[1] * 2.000000 * sc[6];
  J[49] -= dqdci;      // dwdot[O2]/d[H2]
  J[50] += 2 * dqdci;  // dwdot[H2O]/d[H2]
  J[54] += -2 * dqdci; // dwdot[H2]/d[H2]
  // d()/dT
  J[57] -= dqdT;      // dwdot[O2]/dT
  J[58] += 2 * dqdT;  // dwdot[H2O]/dT
  J[62] += -2 * dqdT; // dwdot[H2]/dT

  // reaction 3: CO + H2O <=> CO2 + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[4];
  k_f = 2750000 * exp(0 * tc[0] - 0.50321666580472 * (20000) * invT);
  dlnkfdT = 0 * invT + 0.50321666580472 * (20000) * invT2;
  // reverse
  phi_r = sc[5] * sc[6];
  Kc = exp(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (h_RT[5] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // H2O
  wdot[4] -= q; // CO
  wdot[5] += q; // CO2
  wdot[6] += q; // H2
  // d()/d[H2O]
  dqdci = +k_f * sc[4];
  J[18] -= dqdci; // dwdot[H2O]/d[H2O]
  J[20] -= dqdci; // dwdot[CO]/d[H2O]
  J[21] += dqdci; // dwdot[CO2]/d[H2O]
  J[22] += dqdci; // dwdot[H2]/d[H2O]
  // d()/d[CO]
  dqdci = +k_f * sc[2];
  J[34] -= dqdci; // dwdot[H2O]/d[CO]
  J[36] -= dqdci; // dwdot[CO]/d[CO]
  J[37] += dqdci; // dwdot[CO2]/d[CO]
  J[38] += dqdci; // dwdot[H2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[6];
  J[42] -= dqdci; // dwdot[H2O]/d[CO2]
  J[44] -= dqdci; // dwdot[CO]/d[CO2]
  J[45] += dqdci; // dwdot[CO2]/d[CO2]
  J[46] += dqdci; // dwdot[H2]/d[CO2]
  // d()/d[H2]
  dqdci = -k_r * sc[5];
  J[50] -= dqdci; // dwdot[H2O]/d[H2]
  J[52] -= dqdci; // dwdot[CO]/d[H2]
  J[53] += dqdci; // dwdot[CO2]/d[H2]
  J[54] += dqdci; // dwdot[H2]/d[H2]
  // d()/dT
  J[58] -= dqdT; // dwdot[H2O]/dT
  J[60] -= dqdT; // dwdot[CO]/dT
  J[61] += dqdT; // dwdot[CO2]/dT
  J[62] += dqdT; // dwdot[H2]/dT

  amrex::Real c_R[7], dcRdT[7], e_RT[7];
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
  for (int k = 0; k < 7; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[56 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 7; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 7; ++m) {
      dehmixdc += eh_RT[m] * J[k * 8 + m];
    }
    J[k * 8 + 7] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[63] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
}

// compute an approx to the reaction Jacobian (for preconditioning)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT_SIMPLIFIED(
  amrex::Real* J, amrex::Real* sc, amrex::Real* Tp, const int* HP)
{
  amrex::Real c[7];

  for (int k = 0; k < 7; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian_precond(J, c, *Tp, *HP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 7; k++) {
    J[56 + k] *= 1.e-6;
    J[k * 8 + 7] *= 1.e6;
  }

  return;
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian(amrex::Real* J, amrex::Real* sc, amrex::Real T, const int consP)
{

  for (int i = 0; i < 64; i++) {
    J[i] = 0.0;
  }

  amrex::Real wdot[7];
  for (int k = 0; k < 7; k++) {
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
  for (int k = 0; k < 7; ++k) {
    mixture += sc[k];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[7];
  gibbs(g_RT, tc);

  // compute the species enthalpy
  amrex::Real h_RT[7];
  speciesEnthalpy(h_RT, tc);

  amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  amrex::Real dqdci, dcdc_fac, dqdc[7];
  amrex::Real Pr, fPr, F, k_0, logPr;
  amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
  amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
  amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const amrex::Real ln10 = log(10.0);
  const amrex::Real log10e = 1.0 / log(10.0);
  // reaction 0: 2 CH4 + O2 => 2 CO + 4 H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[0] * sc[0]) * sc[1];
  k_f = 219875.458149426 * exp(0 * tc[0] - 0.50321666580472 * (30000) * invT);
  dlnkfdT = 0 * invT + 0.50321666580472 * (30000) * invT2;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[0] -= 2 * q; // CH4
  wdot[1] -= q;     // O2
  wdot[4] += 2 * q; // CO
  wdot[6] += 4 * q; // H2
  // d()/d[CH4]
  dqdci = +k_f * 2.000000 * sc[0] * sc[1];
  J[0] += -2 * dqdci; // dwdot[CH4]/d[CH4]
  J[1] -= dqdci;      // dwdot[O2]/d[CH4]
  J[4] += 2 * dqdci;  // dwdot[CO]/d[CH4]
  J[6] += 4 * dqdci;  // dwdot[H2]/d[CH4]
  // d()/d[O2]
  dqdci = +k_f * pow(sc[0], 2.000000);
  J[8] += -2 * dqdci; // dwdot[CH4]/d[O2]
  J[9] -= dqdci;      // dwdot[O2]/d[O2]
  J[12] += 2 * dqdci; // dwdot[CO]/d[O2]
  J[14] += 4 * dqdci; // dwdot[H2]/d[O2]
  // d()/dT
  J[56] += -2 * dqdT; // dwdot[CH4]/dT
  J[57] -= dqdT;      // dwdot[O2]/dT
  J[60] += 2 * dqdT;  // dwdot[CO]/dT
  J[62] += 4 * dqdT;  // dwdot[H2]/dT

  // reaction 1: CH4 + H2O <=> CO + 3 H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[2];
  k_f = 300000 * exp(0 * tc[0] - 0.50321666580472 * (30000) * invT);
  dlnkfdT = 0 * invT + 0.50321666580472 * (30000) * invT2;
  // reverse
  phi_r = sc[4] * pow(sc[6], 3.000000);
  Kc = (refC * refC) * exp(g_RT[0] + g_RT[2] - g_RT[4] - 3.000000 * g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT =
    invT * (-(h_RT[0] + h_RT[2]) + (h_RT[4] + 3.000000 * h_RT[6]) - 2.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[0] -= q;     // CH4
  wdot[2] -= q;     // H2O
  wdot[4] += q;     // CO
  wdot[6] += 3 * q; // H2
  // d()/d[CH4]
  dqdci = +k_f * sc[2];
  J[0] -= dqdci;     // dwdot[CH4]/d[CH4]
  J[2] -= dqdci;     // dwdot[H2O]/d[CH4]
  J[4] += dqdci;     // dwdot[CO]/d[CH4]
  J[6] += 3 * dqdci; // dwdot[H2]/d[CH4]
  // d()/d[H2O]
  dqdci = +k_f * sc[0];
  J[16] -= dqdci;     // dwdot[CH4]/d[H2O]
  J[18] -= dqdci;     // dwdot[H2O]/d[H2O]
  J[20] += dqdci;     // dwdot[CO]/d[H2O]
  J[22] += 3 * dqdci; // dwdot[H2]/d[H2O]
  // d()/d[CO]
  dqdci = -k_r * pow(sc[6], 3.000000);
  J[32] -= dqdci;     // dwdot[CH4]/d[CO]
  J[34] -= dqdci;     // dwdot[H2O]/d[CO]
  J[36] += dqdci;     // dwdot[CO]/d[CO]
  J[38] += 3 * dqdci; // dwdot[H2]/d[CO]
  // d()/d[H2]
  dqdci = -k_r * sc[4] * 3.000000 * pow(sc[6], 2.000000);
  J[48] -= dqdci;     // dwdot[CH4]/d[H2]
  J[50] -= dqdci;     // dwdot[H2O]/d[H2]
  J[52] += dqdci;     // dwdot[CO]/d[H2]
  J[54] += 3 * dqdci; // dwdot[H2]/d[H2]
  // d()/dT
  J[56] -= dqdT;     // dwdot[CH4]/dT
  J[58] -= dqdT;     // dwdot[H2O]/dT
  J[60] += dqdT;     // dwdot[CO]/dT
  J[62] += 3 * dqdT; // dwdot[H2]/dT

  // reaction 2: 2 H2 + O2 => 2 H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * (sc[6] * sc[6]);
  k_f = 33993533107.7566 * exp(-1 * tc[0] - 0.50321666580472 * (40000) * invT);
  dlnkfdT = -1 * invT + 0.50321666580472 * (40000) * invT2;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[1] -= q;     // O2
  wdot[2] += 2 * q; // H2O
  wdot[6] -= 2 * q; // H2
  // d()/d[O2]
  dqdci = +k_f * pow(sc[6], 2.000000);
  J[9] -= dqdci;       // dwdot[O2]/d[O2]
  J[10] += 2 * dqdci;  // dwdot[H2O]/d[O2]
  J[14] += -2 * dqdci; // dwdot[H2]/d[O2]
  // d()/d[H2]
  dqdci = +k_f * sc[1] * 2.000000 * sc[6];
  J[49] -= dqdci;      // dwdot[O2]/d[H2]
  J[50] += 2 * dqdci;  // dwdot[H2O]/d[H2]
  J[54] += -2 * dqdci; // dwdot[H2]/d[H2]
  // d()/dT
  J[57] -= dqdT;      // dwdot[O2]/dT
  J[58] += 2 * dqdT;  // dwdot[H2O]/dT
  J[62] += -2 * dqdT; // dwdot[H2]/dT

  // reaction 3: CO + H2O <=> CO2 + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[4];
  k_f = 2750000 * exp(0 * tc[0] - 0.50321666580472 * (20000) * invT);
  dlnkfdT = 0 * invT + 0.50321666580472 * (20000) * invT2;
  // reverse
  phi_r = sc[5] * sc[6];
  Kc = exp(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (h_RT[5] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // H2O
  wdot[4] -= q; // CO
  wdot[5] += q; // CO2
  wdot[6] += q; // H2
  // d()/d[H2O]
  dqdci = +k_f * sc[4];
  J[18] -= dqdci; // dwdot[H2O]/d[H2O]
  J[20] -= dqdci; // dwdot[CO]/d[H2O]
  J[21] += dqdci; // dwdot[CO2]/d[H2O]
  J[22] += dqdci; // dwdot[H2]/d[H2O]
  // d()/d[CO]
  dqdci = +k_f * sc[2];
  J[34] -= dqdci; // dwdot[H2O]/d[CO]
  J[36] -= dqdci; // dwdot[CO]/d[CO]
  J[37] += dqdci; // dwdot[CO2]/d[CO]
  J[38] += dqdci; // dwdot[H2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[6];
  J[42] -= dqdci; // dwdot[H2O]/d[CO2]
  J[44] -= dqdci; // dwdot[CO]/d[CO2]
  J[45] += dqdci; // dwdot[CO2]/d[CO2]
  J[46] += dqdci; // dwdot[H2]/d[CO2]
  // d()/d[H2]
  dqdci = -k_r * sc[5];
  J[50] -= dqdci; // dwdot[H2O]/d[H2]
  J[52] -= dqdci; // dwdot[CO]/d[H2]
  J[53] += dqdci; // dwdot[CO2]/d[H2]
  J[54] += dqdci; // dwdot[H2]/d[H2]
  // d()/dT
  J[58] -= dqdT; // dwdot[H2O]/dT
  J[60] -= dqdT; // dwdot[CO]/dT
  J[61] += dqdT; // dwdot[CO2]/dT
  J[62] += dqdT; // dwdot[H2]/dT

  amrex::Real c_R[7], dcRdT[7], e_RT[7];
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
  for (int k = 0; k < 7; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[56 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 7; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 7; ++m) {
      dehmixdc += eh_RT[m] * J[k * 8 + m];
    }
    J[k * 8 + 7] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[63] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT(amrex::Real* J, amrex::Real* sc, amrex::Real* Tp, const int* consP)
{
  amrex::Real c[7];

  for (int k = 0; k < 7; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian(J, c, *Tp, *consP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 7; k++) {
    J[56 + k] *= 1.e-6;
    J[k * 8 + 7] *= 1.e6;
  }

  return;
}

// Transport function declarations

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENIMC(int* LENIMC)
{
  *LENIMC = 29;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENRMC(int* LENRMC)
{
  *LENRMC = 1148;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNO(int* NO)
{
  *NO = 4;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKK(int* KK)
{
  *KK = 7;
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
  WT[0] = 1.60430000E+01;
  WT[1] = 3.19980000E+01;
  WT[2] = 1.80150000E+01;
  WT[3] = 2.80140000E+01;
  WT[4] = 2.80100000E+01;
  WT[5] = 4.40090000E+01;
  WT[6] = 2.01600000E+00;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetEPS(amrex::Real* EPS)
{
  EPS[0] = 1.41400000E+02;
  EPS[1] = 1.07400000E+02;
  EPS[2] = 5.72400000E+02;
  EPS[3] = 9.75300000E+01;
  EPS[4] = 9.81000000E+01;
  EPS[5] = 2.44000000E+02;
  EPS[6] = 3.80000000E+01;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetSIG(amrex::Real* SIG)
{
  SIG[0] = 3.74600000E+00;
  SIG[1] = 3.45800000E+00;
  SIG[2] = 2.60500000E+00;
  SIG[3] = 3.62100000E+00;
  SIG[4] = 3.65000000E+00;
  SIG[5] = 3.76300000E+00;
  SIG[6] = 2.92000000E+00;
}

// the dipole moment in Debye
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetDIP(amrex::Real* DIP)
{
  DIP[0] = 0.00000000E+00;
  DIP[1] = 0.00000000E+00;
  DIP[2] = 1.84400000E+00;
  DIP[3] = 0.00000000E+00;
  DIP[4] = 0.00000000E+00;
  DIP[5] = 0.00000000E+00;
  DIP[6] = 0.00000000E+00;
}

// the polarizability in cubic Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPOL(amrex::Real* POL)
{
  POL[0] = 2.60000000E+00;
  POL[1] = 1.60000000E+00;
  POL[2] = 0.00000000E+00;
  POL[3] = 1.76000000E+00;
  POL[4] = 1.95000000E+00;
  POL[5] = 2.65000000E+00;
  POL[6] = 7.90000000E-01;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetZROT(amrex::Real* ZROT)
{
  ZROT[0] = 1.30000000E+01;
  ZROT[1] = 3.80000000E+00;
  ZROT[2] = 4.00000000E+00;
  ZROT[3] = 4.00000000E+00;
  ZROT[4] = 1.80000000E+00;
  ZROT[5] = 2.10000000E+00;
  ZROT[6] = 2.80000000E+02;
}

// 0: monoatomic, 1: linear, 2: nonlinear
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLIN(int* NLIN)
{
  NLIN[0] = 2;
  NLIN[1] = 1;
  NLIN[2] = 2;
  NLIN[3] = 1;
  NLIN[4] = 1;
  NLIN[5] = 1;
  NLIN[6] = 1;
}

// Poly fits for the viscosities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFETA(amrex::Real* COFETA)
{
  COFETA[0] = -2.04642228E+01;
  COFETA[1] = 3.75363357E+00;
  COFETA[2] = -4.11756834E-01;
  COFETA[3] = 1.81766117E-02;
  COFETA[4] = -1.78915955E+01;
  COFETA[5] = 2.98311502E+00;
  COFETA[6] = -3.14105508E-01;
  COFETA[7] = 1.40500162E-02;
  COFETA[8] = -1.14441704E+01;
  COFETA[9] = -9.67162726E-01;
  COFETA[10] = 3.58651315E-01;
  COFETA[11] = -2.09789188E-02;
  COFETA[12] = -1.73976840E+01;
  COFETA[13] = 2.73482764E+00;
  COFETA[14] = -2.81916034E-01;
  COFETA[15] = 1.26588405E-02;
  COFETA[16] = -1.74470078E+01;
  COFETA[17] = 2.74728386E+00;
  COFETA[18] = -2.83509015E-01;
  COFETA[19] = 1.27267083E-02;
  COFETA[20] = -2.28110458E+01;
  COFETA[21] = 4.62954710E+00;
  COFETA[22] = -5.00689001E-01;
  COFETA[23] = 2.10012969E-02;
  COFETA[24] = -1.40419383E+01;
  COFETA[25] = 1.08789225E+00;
  COFETA[26] = -6.18592115E-02;
  COFETA[27] = 2.86838304E-03;
}

// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFLAM(amrex::Real* COFLAM)
{
  COFLAM[0] = 1.76957806E+01;
  COFLAM[1] = -6.71274224E+00;
  COFLAM[2] = 1.26299563E+00;
  COFLAM[3] = -6.62747841E-02;
  COFLAM[4] = 5.31587443E-01;
  COFLAM[5] = 1.87067575E+00;
  COFLAM[6] = -1.31586377E-01;
  COFLAM[7] = 5.22416988E-03;
  COFLAM[8] = 1.81612070E+01;
  COFLAM[9] = -6.74136791E+00;
  COFLAM[10] = 1.21372083E+00;
  COFLAM[11] = -6.11027797E-02;
  COFLAM[12] = 7.77701634E+00;
  COFLAM[13] = -1.30957316E+00;
  COFLAM[14] = 3.28841780E-01;
  COFLAM[15] = -1.69485282E-02;
  COFLAM[16] = 8.17515710E+00;
  COFLAM[17] = -1.53836161E+00;
  COFLAM[18] = 3.68036538E-01;
  COFLAM[19] = -1.90917319E-02;
  COFLAM[20] = -8.74830362E+00;
  COFLAM[21] = 4.79275276E+00;
  COFLAM[22] = -4.18685070E-01;
  COFLAM[23] = 1.35210258E-02;
  COFLAM[24] = 4.34727634E+00;
  COFLAM[25] = 1.55347640E+00;
  COFLAM[26] = -1.60615528E-01;
  COFLAM[27] = 9.89934243E-03;
}

// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFD(amrex::Real* COFD)
{
  COFD[0] = -1.75618461E+01;
  COFD[1] = 4.30617914E+00;
  COFD[2] = -3.46490389E-01;
  COFD[3] = 1.51071405E-02;
  COFD[4] = -1.68102679E+01;
  COFD[5] = 4.01337907E+00;
  COFD[6] = -3.10488902E-01;
  COFD[7] = 1.36288975E-02;
  COFD[8] = -1.95657595E+01;
  COFD[9] = 4.78813616E+00;
  COFD[10] = -3.65976019E-01;
  COFD[11] = 1.42215118E-02;
  COFD[12] = -1.65706931E+01;
  COFD[13] = 3.92005093E+00;
  COFD[14] = -2.99040611E-01;
  COFD[15] = 1.31607610E-02;
  COFD[16] = -1.65940113E+01;
  COFD[17] = 3.92553905E+00;
  COFD[18] = -2.99706984E-01;
  COFD[19] = 1.31876655E-02;
  COFD[20] = -1.93937242E+01;
  COFD[21] = 4.87146645E+00;
  COFD[22] = -4.13323360E-01;
  COFD[23] = 1.77408400E-02;
  COFD[24] = -1.29544979E+01;
  COFD[25] = 2.96758239E+00;
  COFD[26] = -1.76586224E-01;
  COFD[27] = 7.90559536E-03;
  COFD[28] = -1.68102679E+01;
  COFD[29] = 4.01337907E+00;
  COFD[30] = -3.10488902E-01;
  COFD[31] = 1.36288975E-02;
  COFD[32] = -1.60936459E+01;
  COFD[33] = 3.70633871E+00;
  COFD[34] = -2.71897253E-01;
  COFD[35] = 1.20097588E-02;
  COFD[36] = -1.99035490E+01;
  COFD[37] = 5.01694632E+00;
  COFD[38] = -4.08962988E-01;
  COFD[39] = 1.66143402E-02;
  COFD[40] = -1.58214949E+01;
  COFD[41] = 3.60000113E+00;
  COFD[42] = -2.58255120E-01;
  COFD[43] = 1.14251480E-02;
  COFD[44] = -1.58458185E+01;
  COFD[45] = 3.60600362E+00;
  COFD[46] = -2.59019961E-01;
  COFD[47] = 1.14576923E-02;
  COFD[48] = -1.87633989E+01;
  COFD[49] = 4.61060397E+00;
  COFD[50] = -3.83564503E-01;
  COFD[51] = 1.66168246E-02;
  COFD[52] = -1.22181330E+01;
  COFD[53] = 2.70415313E+00;
  COFD[54] = -1.41236971E-01;
  COFD[55] = 6.32236816E-03;
  COFD[56] = -1.95657595E+01;
  COFD[57] = 4.78813616E+00;
  COFD[58] = -3.65976019E-01;
  COFD[59] = 1.42215118E-02;
  COFD[60] = -1.99035490E+01;
  COFD[61] = 5.01694632E+00;
  COFD[62] = -4.08962988E-01;
  COFD[63] = 1.66143402E-02;
  COFD[64] = -1.16123739E+01;
  COFD[65] = 8.27753087E-01;
  COFD[66] = 2.52262500E-01;
  COFD[67] = -1.62567542E-02;
  COFD[68] = -1.99472346E+01;
  COFD[69] = 5.05636620E+00;
  COFD[70] = -4.17733674E-01;
  COFD[71] = 1.71403498E-02;
  COFD[72] = -1.99647324E+01;
  COFD[73] = 5.05179382E+00;
  COFD[74] = -4.16351089E-01;
  COFD[75] = 1.70488542E-02;
  COFD[76] = -1.82187534E+01;
  COFD[77] = 3.93854125E+00;
  COFD[78] = -2.28424577E-01;
  COFD[79] = 7.18603064E-03;
  COFD[80] = -1.73864192E+01;
  COFD[81] = 4.71143120E+00;
  COFD[82] = -3.95288666E-01;
  COFD[83] = 1.70702288E-02;
  COFD[84] = -1.65706931E+01;
  COFD[85] = 3.92005093E+00;
  COFD[86] = -2.99040611E-01;
  COFD[87] = 1.31607610E-02;
  COFD[88] = -1.58214949E+01;
  COFD[89] = 3.60000113E+00;
  COFD[90] = -2.58255120E-01;
  COFD[91] = 1.14251480E-02;
  COFD[92] = -1.99472346E+01;
  COFD[93] = 5.05636620E+00;
  COFD[94] = -4.17733674E-01;
  COFD[95] = 1.71403498E-02;
  COFD[96] = -1.56019684E+01;
  COFD[97] = 3.51542686E+00;
  COFD[98] = -2.47677471E-01;
  COFD[99] = 1.09841319E-02;
  COFD[100] = -1.56221646E+01;
  COFD[101] = 3.51977302E+00;
  COFD[102] = -2.48210923E-01;
  COFD[103] = 1.10059241E-02;
  COFD[104] = -1.85068911E+01;
  COFD[105] = 4.52122572E+00;
  COFD[106] = -3.73088946E-01;
  COFD[107] = 1.62076520E-02;
  COFD[108] = -1.20381551E+01;
  COFD[109] = 2.61421687E+00;
  COFD[110] = -1.28887086E-01;
  COFD[111] = 5.75609167E-03;
  COFD[112] = -1.65940113E+01;
  COFD[113] = 3.92553905E+00;
  COFD[114] = -2.99706984E-01;
  COFD[115] = 1.31876655E-02;
  COFD[116] = -1.58458185E+01;
  COFD[117] = 3.60600362E+00;
  COFD[118] = -2.59019961E-01;
  COFD[119] = 1.14576923E-02;
  COFD[120] = -1.99647324E+01;
  COFD[121] = 5.05179382E+00;
  COFD[122] = -4.16351089E-01;
  COFD[123] = 1.70488542E-02;
  COFD[124] = -1.56221646E+01;
  COFD[125] = 3.51977302E+00;
  COFD[126] = -2.48210923E-01;
  COFD[127] = 1.10059241E-02;
  COFD[128] = -1.56423496E+01;
  COFD[129] = 3.52412711E+00;
  COFD[130] = -2.48745351E-01;
  COFD[131] = 1.10277551E-02;
  COFD[132] = -1.85324272E+01;
  COFD[133] = 4.52748688E+00;
  COFD[134] = -3.73847542E-01;
  COFD[135] = 1.62384117E-02;
  COFD[136] = -1.20607836E+01;
  COFD[137] = 2.61969379E+00;
  COFD[138] = -1.29638429E-01;
  COFD[139] = 5.79050588E-03;
  COFD[140] = -1.93937242E+01;
  COFD[141] = 4.87146645E+00;
  COFD[142] = -4.13323360E-01;
  COFD[143] = 1.77408400E-02;
  COFD[144] = -1.87633989E+01;
  COFD[145] = 4.61060397E+00;
  COFD[146] = -3.83564503E-01;
  COFD[147] = 1.66168246E-02;
  COFD[148] = -1.82187534E+01;
  COFD[149] = 3.93854125E+00;
  COFD[150] = -2.28424577E-01;
  COFD[151] = 7.18603064E-03;
  COFD[152] = -1.85068911E+01;
  COFD[153] = 4.52122572E+00;
  COFD[154] = -3.73088946E-01;
  COFD[155] = 1.62076520E-02;
  COFD[156] = -1.85324272E+01;
  COFD[157] = 4.52748688E+00;
  COFD[158] = -3.73847542E-01;
  COFD[159] = 1.62384117E-02;
  COFD[160] = -2.05810575E+01;
  COFD[161] = 5.07469434E+00;
  COFD[162] = -4.25340301E-01;
  COFD[163] = 1.76800795E-02;
  COFD[164] = -1.43978813E+01;
  COFD[165] = 3.49721576E+00;
  COFD[166] = -2.45465191E-01;
  COFD[167] = 1.08948372E-02;
  COFD[168] = -1.29544979E+01;
  COFD[169] = 2.96758239E+00;
  COFD[170] = -1.76586224E-01;
  COFD[171] = 7.90559536E-03;
  COFD[172] = -1.22181330E+01;
  COFD[173] = 2.70415313E+00;
  COFD[174] = -1.41236971E-01;
  COFD[175] = 6.32236816E-03;
  COFD[176] = -1.73864192E+01;
  COFD[177] = 4.71143120E+00;
  COFD[178] = -3.95288666E-01;
  COFD[179] = 1.70702288E-02;
  COFD[180] = -1.20381551E+01;
  COFD[181] = 2.61421687E+00;
  COFD[182] = -1.28887086E-01;
  COFD[183] = 5.75609167E-03;
  COFD[184] = -1.20607836E+01;
  COFD[185] = 2.61969379E+00;
  COFD[186] = -1.29638429E-01;
  COFD[187] = 5.79050588E-03;
  COFD[188] = -1.43978813E+01;
  COFD[189] = 3.49721576E+00;
  COFD[190] = -2.45465191E-01;
  COFD[191] = 1.08948372E-02;
  COFD[192] = -1.04285243E+01;
  COFD[193] = 2.23477534E+00;
  COFD[194] = -8.11809423E-02;
  COFD[195] = 3.77342041E-03;
}

// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKTDIF(int* KTDIF)
{
  KTDIF[0] = 6;
}

// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFTD(amrex::Real* COFTD)
{
  COFTD[0] = 2.98971544E-01;
  COFTD[1] = 2.32229117E-04;
  COFTD[2] = -1.23674024E-07;
  COFTD[3] = 2.01028089E-11;
  COFTD[4] = 3.81861526E-01;
  COFTD[5] = 1.84116078E-04;
  COFTD[6] = -9.79610689E-08;
  COFTD[7] = 1.62541101E-11;
  COFTD[8] = 1.95116681E-03;
  COFTD[9] = 6.69463662E-04;
  COFTD[10] = -3.12145327E-07;
  COFTD[11] = 4.52944954E-11;
  COFTD[12] = 3.88180656E-01;
  COFTD[13] = 1.55380030E-04;
  COFTD[14] = -8.20879923E-08;
  COFTD[15] = 1.37104470E-11;
  COFTD[16] = 3.87402549E-01;
  COFTD[17] = 1.56882676E-04;
  COFTD[18] = -8.29303863E-08;
  COFTD[19] = 1.38459310E-11;
  COFTD[20] = 2.47127846E-01;
  COFTD[21] = 4.49393559E-04;
  COFTD[22] = -2.32029646E-07;
  COFTD[23] = 3.62577088E-11;
  COFTD[24] = 0.00000000E+00;
  COFTD[25] = 0.00000000E+00;
  COFTD[26] = 0.00000000E+00;
  COFTD[27] = 0.00000000E+00;
}

// compute the critical parameters for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS(
  amrex::Real* Tci, amrex::Real* ai, amrex::Real* bi, amrex::Real* acentric_i)
{

  amrex::Real EPS[7];
  amrex::Real SIG[7];
  amrex::Real wt[7];
  amrex::Real avogadro = 6.02214199e23;
  amrex::Real boltzmann = 1.3806503e-16; // we work in CGS
  amrex::Real Rcst = 83.144598;          // in bar [CGS] !

  egtransetEPS(EPS);
  egtransetSIG(SIG);
  get_mw(wt);

  // species 0: CH4
  // Imported from NIST
  Tci[0] = 190.560000;
  ai[0] = 1e6 * 0.42748 * Rcst * Rcst * Tci[0] * Tci[0] /
          (16.043030 * 16.043030 * 45.990000);
  bi[0] = 0.08664 * Rcst * Tci[0] / (16.043030 * 45.990000);
  acentric_i[0] = 0.011000;

  // species 1: O2
  // Imported from NIST
  Tci[1] = 154.581000;
  ai[1] = 1e6 * 0.42748 * Rcst * Rcst * Tci[1] * Tci[1] /
          (31.998800 * 31.998800 * 50.430466);
  bi[1] = 0.08664 * Rcst * Tci[1] / (31.998800 * 50.430466);
  acentric_i[1] = 0.022200;

  // species 2: H2O
  // Imported from NIST
  Tci[2] = 647.096000;
  ai[2] = 1e6 * 0.42748 * Rcst * Rcst * Tci[2] * Tci[2] /
          (18.015340 * 18.015340 * 220.640000);
  bi[2] = 0.08664 * Rcst * Tci[2] / (18.015340 * 220.640000);
  acentric_i[2] = 0.344300;

  // species 3: N2
  // Imported from NIST
  Tci[3] = 126.192000;
  ai[3] = 1e6 * 0.42748 * Rcst * Rcst * Tci[3] * Tci[3] /
          (28.013400 * 28.013400 * 33.958000);
  bi[3] = 0.08664 * Rcst * Tci[3] / (28.013400 * 33.958000);
  acentric_i[3] = 0.037200;

  // species 4: CO
  // Imported from NIST
  Tci[4] = 132.850000;
  ai[4] = 1e6 * 0.42748 * Rcst * Rcst * Tci[4] * Tci[4] /
          (28.010000 * 28.010000 * 34.940000);
  bi[4] = 0.08664 * Rcst * Tci[4] / (28.010000 * 34.940000);
  acentric_i[4] = 0.045000;

  // species 5: CO2
  // Imported from NIST
  Tci[5] = 304.120000;
  ai[5] = 1e6 * 0.42748 * Rcst * Rcst * Tci[5] * Tci[5] /
          (44.009950 * 44.009950 * 73.740000);
  bi[5] = 0.08664 * Rcst * Tci[5] / (44.009950 * 73.740000);
  acentric_i[5] = 0.225000;

  // species 6: H2
  // Imported from NIST
  Tci[6] = 33.145000;
  ai[6] = 1e6 * 0.42748 * Rcst * Rcst * Tci[6] * Tci[6] /
          (2.015880 * 2.015880 * 12.964000);
  bi[6] = 0.08664 * Rcst * Tci[6] / (2.015880 * 12.964000);
  acentric_i[6] = -0.219000;

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
