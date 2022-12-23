#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  N
1  H
2  O
3  C
4  E
*/

// Species
#define N2_ID 0
#define H2_ID 1
#define H_ID 2
#define O2_ID 3
#define O_ID 4
#define H2O_ID 5
#define OH_ID 6
#define H2O2_ID 7
#define HO2_ID 8
#define CO_ID 9
#define CO2_ID 10
#define CH4_ID 11
#define CH3_ID 12
#define CH2_ID 13
#define CH2S_ID 14
#define C_ID 15
#define CH_ID 16
#define CH3OH_ID 17
#define CH3O_ID 18
#define CH2OH_ID 19
#define CH2O_ID 20
#define HCO_ID 21
#define CHOp_ID 22
#define H3Op_ID 23
#define O2n_ID 24
#define E_ID 25

#define NUM_ELEMENTS 5
#define NUM_SPECIES 26
#define NUM_IONS 4
#define NUM_REACTIONS 134

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
  *mm = 5;
  *kk = 26;
  *ii = 134;
  *nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_imw(amrex::Real* imw_new)
{
  imw_new[0] = 1.0 / 28.014000;  // N2
  imw_new[1] = 1.0 / 2.016000;   // H2
  imw_new[2] = 1.0 / 1.008000;   // H
  imw_new[3] = 1.0 / 31.998000;  // O2
  imw_new[4] = 1.0 / 15.999000;  // O
  imw_new[5] = 1.0 / 18.015000;  // H2O
  imw_new[6] = 1.0 / 17.007000;  // OH
  imw_new[7] = 1.0 / 34.014000;  // H2O2
  imw_new[8] = 1.0 / 33.006000;  // HO2
  imw_new[9] = 1.0 / 28.010000;  // CO
  imw_new[10] = 1.0 / 44.009000; // CO2
  imw_new[11] = 1.0 / 16.043000; // CH4
  imw_new[12] = 1.0 / 15.035000; // CH3
  imw_new[13] = 1.0 / 14.027000; // CH2
  imw_new[14] = 1.0 / 14.027000; // CH2(S)
  imw_new[15] = 1.0 / 12.011000; // C
  imw_new[16] = 1.0 / 13.019000; // CH
  imw_new[17] = 1.0 / 32.042000; // CH3OH
  imw_new[18] = 1.0 / 31.034000; // CH3O
  imw_new[19] = 1.0 / 31.034000; // CH2OH
  imw_new[20] = 1.0 / 30.026000; // CH2O
  imw_new[21] = 1.0 / 29.018000; // HCO
  imw_new[22] = 1.0 / 29.017451; // CHO+
  imw_new[23] = 1.0 / 19.022451; // H3O+
  imw_new[24] = 1.0 / 31.998549; // O2-
  imw_new[25] = 1.0 / 0.000549;  // E
}

//  molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_mw(amrex::Real* mw_new)
{
  mw_new[0] = 28.014000;  // N2
  mw_new[1] = 2.016000;   // H2
  mw_new[2] = 1.008000;   // H
  mw_new[3] = 31.998000;  // O2
  mw_new[4] = 15.999000;  // O
  mw_new[5] = 18.015000;  // H2O
  mw_new[6] = 17.007000;  // OH
  mw_new[7] = 34.014000;  // H2O2
  mw_new[8] = 33.006000;  // HO2
  mw_new[9] = 28.010000;  // CO
  mw_new[10] = 44.009000; // CO2
  mw_new[11] = 16.043000; // CH4
  mw_new[12] = 15.035000; // CH3
  mw_new[13] = 14.027000; // CH2
  mw_new[14] = 14.027000; // CH2(S)
  mw_new[15] = 12.011000; // C
  mw_new[16] = 13.019000; // CH
  mw_new[17] = 32.042000; // CH3OH
  mw_new[18] = 31.034000; // CH3O
  mw_new[19] = 31.034000; // CH2OH
  mw_new[20] = 30.026000; // CH2O
  mw_new[21] = 29.018000; // HCO
  mw_new[22] = 29.017451; // CHO+
  mw_new[23] = 19.022451; // H3O+
  mw_new[24] = 31.998549; // O2-
  mw_new[25] = 0.000549;  // E
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

  // species with midpoint at T=1050 kelvin
  if (T < 1050) {
    // species 0: N2
    species[0] = +2.73100682e+00 - 1.83159730e-03 * tc[1] +
                 4.32324662e-06 * tc[2] - 3.04378151e-09 * tc[3] +
                 7.46071562e-13 * tc[4];
  } else {
    // species 0: N2
    species[0] = +1.81166073e+00 + 1.67067353e-03 * tc[1] -
                 6.79997428e-07 * tc[2] + 1.32881379e-10 * tc[3] -
                 1.02767442e-14 * tc[4];
  }

  // species with midpoint at T=700 kelvin
  if (T < 700) {
    // species 1: H2
    species[1] = +1.64204438e+00 + 5.49529274e-03 * tc[1] -
                 1.27163634e-05 * tc[2] + 1.28749173e-08 * tc[3] -
                 4.70027749e-12 * tc[4];
    // species 3: O2
    species[3] = +2.74403921e+00 - 2.79740147e-03 * tc[1] +
                 9.80122558e-06 * tc[2] - 1.03259643e-08 * tc[3] +
                 3.79931247e-12 * tc[4];
    // species 8: HO2
    species[8] = +2.61994299e+00 + 1.05805705e-03 * tc[1] +
                 5.06678941e-06 * tc[2] - 6.33800762e-09 * tc[3] +
                 2.41597281e-12 * tc[4];
    // species 11: CH4
    species[11] = +4.23967335e+00 - 1.46835123e-02 * tc[1] +
                  5.29732711e-05 * tc[2] - 5.41668822e-08 * tc[3] +
                  1.96318566e-11 * tc[4];
    // species 15: C
    species[15] = +1.54495192e+00 - 2.47725281e-04 * tc[1] +
                  5.48018277e-07 * tc[2] - 5.48551250e-10 * tc[3] +
                  2.04117331e-13 * tc[4];
    // species 20: CH2O
    species[20] = +3.32621296e+00 - 7.01151853e-03 * tc[1] +
                  3.15176961e-05 * tc[2] - 3.36478638e-08 * tc[3] +
                  1.23454023e-11 * tc[4];
  } else {
    // species 1: H2
    species[1] = +2.78199881e+00 - 1.01873259e-03 * tc[1] +
                 1.24226233e-06 * tc[2] - 4.19011898e-10 * tc[3] +
                 4.75543793e-14 * tc[4];
    // species 3: O2
    species[3] = +1.82012408e+00 + 2.48211357e-03 * tc[1] -
                 1.51202094e-06 * tc[2] + 4.48556201e-10 * tc[3] -
                 4.87305668e-14 * tc[4];
    // species 8: HO2
    species[8] = +2.02391889e+00 + 4.46390907e-03 * tc[1] -
                 2.23146492e-06 * tc[2] + 6.12710799e-10 * tc[3] -
                 6.64266237e-14 * tc[4];
    // species 11: CH4
    species[11] = -4.94653595e-01 + 1.23697845e-02 * tc[1] -
                  4.99807922e-06 * tc[2] + 1.04392765e-09 * tc[3] -
                  8.62897416e-14 * tc[4];
    // species 15: C
    species[15] = +1.49472531e+00 + 3.92839476e-05 * tc[1] -
                  6.70014980e-08 * tc[2] + 3.71818694e-11 * tc[3] -
                  5.07306885e-15 * tc[4];
    // species 20: CH2O
    species[20] = +3.33356520e-01 + 1.00905183e-02 * tc[1] -
                  5.12952562e-06 * tc[2] + 1.25425207e-09 * tc[3] -
                  1.19639109e-13 * tc[4];
  }

  // species with midpoint at T=860 kelvin
  if (T < 860) {
    // species 2: H
    species[2] = +1.49950544e+00 + 2.99164046e-06 * tc[1] -
                 5.92759759e-09 * tc[2] + 4.87810165e-12 * tc[3] -
                 1.45539320e-15 * tc[4];
  } else {
    // species 2: H
    species[2] = +1.50031493e+00 - 7.73406828e-07 * tc[1] +
                 6.39345346e-10 * tc[2] - 2.12551791e-13 * tc[3] +
                 2.44479191e-17 * tc[4];
  }

  // species with midpoint at T=720 kelvin
  if (T < 720) {
    // species 4: O
    species[4] = +2.14799201e+00 - 3.11174065e-03 * tc[1] +
                 6.18137897e-06 * tc[2] - 5.63808798e-09 * tc[3] +
                 1.94866016e-12 * tc[4];
  } else {
    // species 4: O
    species[4] = +1.62549143e+00 - 2.08959644e-04 * tc[1] +
                 1.33918546e-07 * tc[2] - 3.85875896e-11 * tc[3] +
                 4.38918689e-15 * tc[4];
  }

  // species with midpoint at T=1420 kelvin
  if (T < 1420) {
    // species 5: H2O
    species[5] = +3.06061172e+00 - 8.65807189e-04 * tc[1] +
                 3.24409528e-06 * tc[2] - 1.80243079e-09 * tc[3] +
                 3.32483293e-13 * tc[4];
  } else {
    // species 5: H2O
    species[5] = +1.66777075e+00 + 3.05768849e-03 * tc[1] -
                 9.00442411e-07 * tc[2] + 1.43361552e-10 * tc[3] -
                 1.00857817e-14 * tc[4];
  }

  // species with midpoint at T=1700 kelvin
  if (T < 1700) {
    // species 6: OH
    species[6] = +2.91354631e+00 - 1.66275926e-03 * tc[1] +
                 2.30920029e-06 * tc[2] - 1.02359508e-09 * tc[3] +
                 1.58829629e-13 * tc[4];
  } else {
    // species 6: OH
    species[6] = +1.49867369e+00 + 1.66635279e-03 * tc[1] -
                 6.28251516e-07 * tc[2] + 1.28346806e-10 * tc[3] -
                 1.05735894e-14 * tc[4];
  }

  // species with midpoint at T=1800 kelvin
  if (T < 1800) {
    // species 7: H2O2
    species[7] = +2.34774224e+00 + 7.05005437e-03 * tc[1] -
                 3.84522006e-06 * tc[2] + 1.16720661e-09 * tc[3] -
                 1.47618105e-13 * tc[4];
    // species 13: CH2
    species[13] = +2.76489460e+00 + 1.43839191e-03 * tc[1] +
                  4.75583077e-07 * tc[2] - 4.31788591e-10 * tc[3] +
                  7.58292874e-14 * tc[4];
    // species 17: CH3OH
    species[17] = -1.52669521e-01 + 1.63086904e-02 * tc[1] -
                  8.48344961e-06 * tc[2] + 2.29304341e-09 * tc[3] -
                  2.59952013e-13 * tc[4];
  } else {
    // species 7: H2O2
    species[7] = +3.76869639e+00 + 3.89237848e-03 * tc[1] -
                 1.21382349e-06 * tc[2] + 1.92615285e-10 * tc[3] -
                 1.22581990e-14 * tc[4];
    // species 13: CH2
    species[13] = +1.81272972e+00 + 3.55431388e-03 * tc[1] -
                  1.28768523e-06 * tc[2] + 2.21273744e-10 * tc[3] -
                  1.48738147e-14 * tc[4];
    // species 17: CH3OH
    species[17] = +1.71701530e+00 + 1.21538353e-02 * tc[1] -
                  5.02107031e-06 * tc[2] + 1.01068070e-09 * tc[3] -
                  8.18460823e-14 * tc[4];
  }

  // species with midpoint at T=960 kelvin
  if (T < 960) {
    // species 9: CO
    species[9] = +2.75723891e+00 - 2.14465241e-03 * tc[1] +
                 5.42079005e-06 * tc[2] - 4.17025963e-09 * tc[3] +
                 1.11901127e-12 * tc[4];
  } else {
    // species 9: CO
    species[9] = +1.79255381e+00 + 1.87486886e-03 * tc[1] -
                 8.59711926e-07 * tc[2] + 1.91200070e-10 * tc[3] -
                 1.67855286e-14 * tc[4];
  }

  // species with midpoint at T=1450 kelvin
  if (T < 1450) {
    // species 10: CO2
    species[10] = +1.31684347e+00 + 9.22755036e-03 * tc[1] -
                  7.75654093e-06 * tc[2] + 3.28225360e-09 * tc[3] -
                  5.48722482e-13 * tc[4];
  } else {
    // species 10: CO2
    species[10] = +3.70876468e+00 + 2.62914704e-03 * tc[1] -
                  9.30606462e-07 * tc[2] + 1.43892920e-10 * tc[3] -
                  7.62581414e-15 * tc[4];
  }

  // species with midpoint at T=1060 kelvin
  if (T < 1060) {
    // species 12: CH3
    species[12] = +2.47829310e+00 + 3.54764773e-03 * tc[1] +
                  1.47408440e-06 * tc[2] - 1.94375955e-09 * tc[3] +
                  5.21921232e-13 * tc[4];
  } else {
    // species 12: CH3
    species[12] = +1.78805104e+00 + 6.15233477e-03 * tc[1] -
                  2.21179349e-06 * tc[2] + 3.74402648e-10 * tc[3] -
                  2.48151349e-14 * tc[4];
  }

  // species with midpoint at T=970 kelvin
  if (T < 970) {
    // species 14: CH2(S)
    species[14] = +3.18185434e+00 - 2.21134314e-03 * tc[1] +
                  7.71527541e-06 * tc[2] - 5.95950381e-09 * tc[3] +
                  1.58314628e-12 * tc[4];
  } else {
    // species 14: CH2(S)
    species[14] = +1.75934299e+00 + 3.65468307e-03 * tc[1] -
                  1.35589914e-06 * tc[2] + 2.74980411e-10 * tc[3] -
                  2.36795472e-14 * tc[4];
  }

  // species with midpoint at T=1590 kelvin
  if (T < 1590) {
    // species 16: CH
    species[16] = +2.77264332e+00 - 1.58547350e-03 * tc[1] +
                  2.83512238e-06 * tc[2] - 1.36146058e-09 * tc[3] +
                  2.23995332e-13 * tc[4];
  } else {
    // species 16: CH
    species[16] = +1.27990128e+00 + 2.16985238e-03 * tc[1] -
                  7.07637885e-07 * tc[2] + 1.23973494e-10 * tc[3] -
                  9.56348515e-15 * tc[4];
  }

  // species with midpoint at T=1740 kelvin
  if (T < 1740) {
    // species 18: CH3O
    species[18] = -1.10339015e-01 + 1.70119767e-02 * tc[1] -
                  1.13807351e-05 * tc[2] + 3.88280928e-09 * tc[3] -
                  5.32841479e-13 * tc[4];
  } else {
    // species 18: CH3O
    species[18] = +4.72238062e+00 + 5.90227637e-03 * tc[1] -
                  1.80340720e-06 * tc[2] + 2.13335009e-10 * tc[3] -
                  5.61816409e-15 * tc[4];
  }

  // species with midpoint at T=1360 kelvin
  if (T < 1360) {
    // species 19: CH2OH
    species[19] = +1.34821579e+00 + 1.39600168e-02 * tc[1] -
                  1.08632206e-05 * tc[2] + 4.62498416e-09 * tc[3] -
                  8.08499162e-13 * tc[4];
  } else {
    // species 19: CH2OH
    species[19] = +4.04534950e+00 + 6.02727059e-03 * tc[1] -
                  2.11386821e-06 * tc[2] + 3.36085905e-10 * tc[3] -
                  2.00987482e-14 * tc[4];
  }

  // species with midpoint at T=770 kelvin
  if (T < 770) {
    // species 21: HCO
    species[21] = +3.03483979e+00 - 2.15836864e-03 * tc[1] +
                  1.18233875e-05 * tc[2] - 1.18459406e-08 * tc[3] +
                  4.00593954e-12 * tc[4];
  } else {
    // species 21: HCO
    species[21] = +1.60049318e+00 + 5.29278258e-03 * tc[1] -
                  2.69184211e-06 * tc[2] + 7.21357798e-10 * tc[3] -
                  7.43521409e-14 * tc[4];
  }

  // species with midpoint at T=3654.18 kelvin
  if (T < 3654.18) {
    // species 22: CHO+
    species[22] = +1.87627680e+00 + 4.94860200e-03 * tc[1] -
                  2.30853060e-06 * tc[2] + 5.11217520e-10 * tc[3] -
                  4.35206310e-14 * tc[4];
  } else {
    // species 22: CHO+
    species[22] = +5.91572650e+00 + 1.66311070e-04 * tc[1] -
                  1.84763430e-08 * tc[2] + 9.08564290e-13 * tc[3] -
                  1.64544410e-17 * tc[4];
  }

  // species with midpoint at T=5430.9 kelvin
  if (T < 5430.9) {
    // species 23: H3O+
    species[23] = +1.68852320e+00 + 5.44536540e-03 * tc[1] -
                  1.70555660e-06 * tc[2] + 2.48332930e-10 * tc[3] -
                  1.37694400e-14 * tc[4];
  } else {
    // species 23: H3O+
    species[23] = +7.95088030e+00 + 2.52816680e-04 * tc[1] -
                  2.45943970e-08 * tc[2] + 1.08783180e-12 * tc[3] -
                  1.80908370e-17 * tc[4];
  }

  // species with midpoint at T=2008.71 kelvin
  if (T < 2008.71) {
    // species 24: O2-
    species[24] = +2.10217180e+00 + 2.79808750e-03 * tc[1] -
                  2.26511260e-06 * tc[2] + 8.69165170e-10 * tc[3] -
                  1.27218840e-13 * tc[4];
  } else {
    // species 24: O2-
    species[24] = +3.25928670e+00 + 2.24680720e-04 * tc[1] -
                  5.13979550e-08 * tc[2] + 7.35459780e-12 * tc[3] -
                  3.85586520e-16 * tc[4];
  }

  // species with midpoint at T=27300 kelvin
  if (T < 27300) {
    // species 25: E
    species[25] = +1.50000000e+00 - 8.15026080e-18 * tc[1] +
                  1.30298670e-21 * tc[2] - 7.26483330e-26 * tc[3] +
                  1.30662470e-30 * tc[4];
  } else {
    // species 25: E
    species[25] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4];
  }
}

// compute Cp/R at the given temperature
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cp_R(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];

  // species with midpoint at T=1050 kelvin
  if (T < 1050) {
    // species 0: N2
    species[0] = +3.73100682e+00 - 1.83159730e-03 * tc[1] +
                 4.32324662e-06 * tc[2] - 3.04378151e-09 * tc[3] +
                 7.46071562e-13 * tc[4];
  } else {
    // species 0: N2
    species[0] = +2.81166073e+00 + 1.67067353e-03 * tc[1] -
                 6.79997428e-07 * tc[2] + 1.32881379e-10 * tc[3] -
                 1.02767442e-14 * tc[4];
  }

  // species with midpoint at T=700 kelvin
  if (T < 700) {
    // species 1: H2
    species[1] = +2.64204438e+00 + 5.49529274e-03 * tc[1] -
                 1.27163634e-05 * tc[2] + 1.28749173e-08 * tc[3] -
                 4.70027749e-12 * tc[4];
    // species 3: O2
    species[3] = +3.74403921e+00 - 2.79740147e-03 * tc[1] +
                 9.80122558e-06 * tc[2] - 1.03259643e-08 * tc[3] +
                 3.79931247e-12 * tc[4];
    // species 8: HO2
    species[8] = +3.61994299e+00 + 1.05805705e-03 * tc[1] +
                 5.06678941e-06 * tc[2] - 6.33800762e-09 * tc[3] +
                 2.41597281e-12 * tc[4];
    // species 11: CH4
    species[11] = +5.23967335e+00 - 1.46835123e-02 * tc[1] +
                  5.29732711e-05 * tc[2] - 5.41668822e-08 * tc[3] +
                  1.96318566e-11 * tc[4];
    // species 15: C
    species[15] = +2.54495192e+00 - 2.47725281e-04 * tc[1] +
                  5.48018277e-07 * tc[2] - 5.48551250e-10 * tc[3] +
                  2.04117331e-13 * tc[4];
    // species 20: CH2O
    species[20] = +4.32621296e+00 - 7.01151853e-03 * tc[1] +
                  3.15176961e-05 * tc[2] - 3.36478638e-08 * tc[3] +
                  1.23454023e-11 * tc[4];
  } else {
    // species 1: H2
    species[1] = +3.78199881e+00 - 1.01873259e-03 * tc[1] +
                 1.24226233e-06 * tc[2] - 4.19011898e-10 * tc[3] +
                 4.75543793e-14 * tc[4];
    // species 3: O2
    species[3] = +2.82012408e+00 + 2.48211357e-03 * tc[1] -
                 1.51202094e-06 * tc[2] + 4.48556201e-10 * tc[3] -
                 4.87305668e-14 * tc[4];
    // species 8: HO2
    species[8] = +3.02391889e+00 + 4.46390907e-03 * tc[1] -
                 2.23146492e-06 * tc[2] + 6.12710799e-10 * tc[3] -
                 6.64266237e-14 * tc[4];
    // species 11: CH4
    species[11] = +5.05346405e-01 + 1.23697845e-02 * tc[1] -
                  4.99807922e-06 * tc[2] + 1.04392765e-09 * tc[3] -
                  8.62897416e-14 * tc[4];
    // species 15: C
    species[15] = +2.49472531e+00 + 3.92839476e-05 * tc[1] -
                  6.70014980e-08 * tc[2] + 3.71818694e-11 * tc[3] -
                  5.07306885e-15 * tc[4];
    // species 20: CH2O
    species[20] = +1.33335652e+00 + 1.00905183e-02 * tc[1] -
                  5.12952562e-06 * tc[2] + 1.25425207e-09 * tc[3] -
                  1.19639109e-13 * tc[4];
  }

  // species with midpoint at T=860 kelvin
  if (T < 860) {
    // species 2: H
    species[2] = +2.49950544e+00 + 2.99164046e-06 * tc[1] -
                 5.92759759e-09 * tc[2] + 4.87810165e-12 * tc[3] -
                 1.45539320e-15 * tc[4];
  } else {
    // species 2: H
    species[2] = +2.50031493e+00 - 7.73406828e-07 * tc[1] +
                 6.39345346e-10 * tc[2] - 2.12551791e-13 * tc[3] +
                 2.44479191e-17 * tc[4];
  }

  // species with midpoint at T=720 kelvin
  if (T < 720) {
    // species 4: O
    species[4] = +3.14799201e+00 - 3.11174065e-03 * tc[1] +
                 6.18137897e-06 * tc[2] - 5.63808798e-09 * tc[3] +
                 1.94866016e-12 * tc[4];
  } else {
    // species 4: O
    species[4] = +2.62549143e+00 - 2.08959644e-04 * tc[1] +
                 1.33918546e-07 * tc[2] - 3.85875896e-11 * tc[3] +
                 4.38918689e-15 * tc[4];
  }

  // species with midpoint at T=1420 kelvin
  if (T < 1420) {
    // species 5: H2O
    species[5] = +4.06061172e+00 - 8.65807189e-04 * tc[1] +
                 3.24409528e-06 * tc[2] - 1.80243079e-09 * tc[3] +
                 3.32483293e-13 * tc[4];
  } else {
    // species 5: H2O
    species[5] = +2.66777075e+00 + 3.05768849e-03 * tc[1] -
                 9.00442411e-07 * tc[2] + 1.43361552e-10 * tc[3] -
                 1.00857817e-14 * tc[4];
  }

  // species with midpoint at T=1700 kelvin
  if (T < 1700) {
    // species 6: OH
    species[6] = +3.91354631e+00 - 1.66275926e-03 * tc[1] +
                 2.30920029e-06 * tc[2] - 1.02359508e-09 * tc[3] +
                 1.58829629e-13 * tc[4];
  } else {
    // species 6: OH
    species[6] = +2.49867369e+00 + 1.66635279e-03 * tc[1] -
                 6.28251516e-07 * tc[2] + 1.28346806e-10 * tc[3] -
                 1.05735894e-14 * tc[4];
  }

  // species with midpoint at T=1800 kelvin
  if (T < 1800) {
    // species 7: H2O2
    species[7] = +3.34774224e+00 + 7.05005437e-03 * tc[1] -
                 3.84522006e-06 * tc[2] + 1.16720661e-09 * tc[3] -
                 1.47618105e-13 * tc[4];
    // species 13: CH2
    species[13] = +3.76489460e+00 + 1.43839191e-03 * tc[1] +
                  4.75583077e-07 * tc[2] - 4.31788591e-10 * tc[3] +
                  7.58292874e-14 * tc[4];
    // species 17: CH3OH
    species[17] = +8.47330479e-01 + 1.63086904e-02 * tc[1] -
                  8.48344961e-06 * tc[2] + 2.29304341e-09 * tc[3] -
                  2.59952013e-13 * tc[4];
  } else {
    // species 7: H2O2
    species[7] = +4.76869639e+00 + 3.89237848e-03 * tc[1] -
                 1.21382349e-06 * tc[2] + 1.92615285e-10 * tc[3] -
                 1.22581990e-14 * tc[4];
    // species 13: CH2
    species[13] = +2.81272972e+00 + 3.55431388e-03 * tc[1] -
                  1.28768523e-06 * tc[2] + 2.21273744e-10 * tc[3] -
                  1.48738147e-14 * tc[4];
    // species 17: CH3OH
    species[17] = +2.71701530e+00 + 1.21538353e-02 * tc[1] -
                  5.02107031e-06 * tc[2] + 1.01068070e-09 * tc[3] -
                  8.18460823e-14 * tc[4];
  }

  // species with midpoint at T=960 kelvin
  if (T < 960) {
    // species 9: CO
    species[9] = +3.75723891e+00 - 2.14465241e-03 * tc[1] +
                 5.42079005e-06 * tc[2] - 4.17025963e-09 * tc[3] +
                 1.11901127e-12 * tc[4];
  } else {
    // species 9: CO
    species[9] = +2.79255381e+00 + 1.87486886e-03 * tc[1] -
                 8.59711926e-07 * tc[2] + 1.91200070e-10 * tc[3] -
                 1.67855286e-14 * tc[4];
  }

  // species with midpoint at T=1450 kelvin
  if (T < 1450) {
    // species 10: CO2
    species[10] = +2.31684347e+00 + 9.22755036e-03 * tc[1] -
                  7.75654093e-06 * tc[2] + 3.28225360e-09 * tc[3] -
                  5.48722482e-13 * tc[4];
  } else {
    // species 10: CO2
    species[10] = +4.70876468e+00 + 2.62914704e-03 * tc[1] -
                  9.30606462e-07 * tc[2] + 1.43892920e-10 * tc[3] -
                  7.62581414e-15 * tc[4];
  }

  // species with midpoint at T=1060 kelvin
  if (T < 1060) {
    // species 12: CH3
    species[12] = +3.47829310e+00 + 3.54764773e-03 * tc[1] +
                  1.47408440e-06 * tc[2] - 1.94375955e-09 * tc[3] +
                  5.21921232e-13 * tc[4];
  } else {
    // species 12: CH3
    species[12] = +2.78805104e+00 + 6.15233477e-03 * tc[1] -
                  2.21179349e-06 * tc[2] + 3.74402648e-10 * tc[3] -
                  2.48151349e-14 * tc[4];
  }

  // species with midpoint at T=970 kelvin
  if (T < 970) {
    // species 14: CH2(S)
    species[14] = +4.18185434e+00 - 2.21134314e-03 * tc[1] +
                  7.71527541e-06 * tc[2] - 5.95950381e-09 * tc[3] +
                  1.58314628e-12 * tc[4];
  } else {
    // species 14: CH2(S)
    species[14] = +2.75934299e+00 + 3.65468307e-03 * tc[1] -
                  1.35589914e-06 * tc[2] + 2.74980411e-10 * tc[3] -
                  2.36795472e-14 * tc[4];
  }

  // species with midpoint at T=1590 kelvin
  if (T < 1590) {
    // species 16: CH
    species[16] = +3.77264332e+00 - 1.58547350e-03 * tc[1] +
                  2.83512238e-06 * tc[2] - 1.36146058e-09 * tc[3] +
                  2.23995332e-13 * tc[4];
  } else {
    // species 16: CH
    species[16] = +2.27990128e+00 + 2.16985238e-03 * tc[1] -
                  7.07637885e-07 * tc[2] + 1.23973494e-10 * tc[3] -
                  9.56348515e-15 * tc[4];
  }

  // species with midpoint at T=1740 kelvin
  if (T < 1740) {
    // species 18: CH3O
    species[18] = +8.89660985e-01 + 1.70119767e-02 * tc[1] -
                  1.13807351e-05 * tc[2] + 3.88280928e-09 * tc[3] -
                  5.32841479e-13 * tc[4];
  } else {
    // species 18: CH3O
    species[18] = +5.72238062e+00 + 5.90227637e-03 * tc[1] -
                  1.80340720e-06 * tc[2] + 2.13335009e-10 * tc[3] -
                  5.61816409e-15 * tc[4];
  }

  // species with midpoint at T=1360 kelvin
  if (T < 1360) {
    // species 19: CH2OH
    species[19] = +2.34821579e+00 + 1.39600168e-02 * tc[1] -
                  1.08632206e-05 * tc[2] + 4.62498416e-09 * tc[3] -
                  8.08499162e-13 * tc[4];
  } else {
    // species 19: CH2OH
    species[19] = +5.04534950e+00 + 6.02727059e-03 * tc[1] -
                  2.11386821e-06 * tc[2] + 3.36085905e-10 * tc[3] -
                  2.00987482e-14 * tc[4];
  }

  // species with midpoint at T=770 kelvin
  if (T < 770) {
    // species 21: HCO
    species[21] = +4.03483979e+00 - 2.15836864e-03 * tc[1] +
                  1.18233875e-05 * tc[2] - 1.18459406e-08 * tc[3] +
                  4.00593954e-12 * tc[4];
  } else {
    // species 21: HCO
    species[21] = +2.60049318e+00 + 5.29278258e-03 * tc[1] -
                  2.69184211e-06 * tc[2] + 7.21357798e-10 * tc[3] -
                  7.43521409e-14 * tc[4];
  }

  // species with midpoint at T=3654.18 kelvin
  if (T < 3654.18) {
    // species 22: CHO+
    species[22] = +2.87627680e+00 + 4.94860200e-03 * tc[1] -
                  2.30853060e-06 * tc[2] + 5.11217520e-10 * tc[3] -
                  4.35206310e-14 * tc[4];
  } else {
    // species 22: CHO+
    species[22] = +6.91572650e+00 + 1.66311070e-04 * tc[1] -
                  1.84763430e-08 * tc[2] + 9.08564290e-13 * tc[3] -
                  1.64544410e-17 * tc[4];
  }

  // species with midpoint at T=5430.9 kelvin
  if (T < 5430.9) {
    // species 23: H3O+
    species[23] = +2.68852320e+00 + 5.44536540e-03 * tc[1] -
                  1.70555660e-06 * tc[2] + 2.48332930e-10 * tc[3] -
                  1.37694400e-14 * tc[4];
  } else {
    // species 23: H3O+
    species[23] = +8.95088030e+00 + 2.52816680e-04 * tc[1] -
                  2.45943970e-08 * tc[2] + 1.08783180e-12 * tc[3] -
                  1.80908370e-17 * tc[4];
  }

  // species with midpoint at T=2008.71 kelvin
  if (T < 2008.71) {
    // species 24: O2-
    species[24] = +3.10217180e+00 + 2.79808750e-03 * tc[1] -
                  2.26511260e-06 * tc[2] + 8.69165170e-10 * tc[3] -
                  1.27218840e-13 * tc[4];
  } else {
    // species 24: O2-
    species[24] = +4.25928670e+00 + 2.24680720e-04 * tc[1] -
                  5.13979550e-08 * tc[2] + 7.35459780e-12 * tc[3] -
                  3.85586520e-16 * tc[4];
  }

  // species with midpoint at T=27300 kelvin
  if (T < 27300) {
    // species 25: E
    species[25] = +2.50000000e+00 - 8.15026080e-18 * tc[1] +
                  1.30298670e-21 * tc[2] - 7.26483330e-26 * tc[3] +
                  1.30662470e-30 * tc[4];
  } else {
    // species 25: E
    species[25] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4];
  }
}

// compute the g/(RT) at the given temperature
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
gibbs(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1 / T;

  // species with midpoint at T=1050 kelvin
  if (T < 1050) {
    // species 0: N2
    species[0] = -1.062874260000000e+03 * invT + 1.562794840000000e+00 -
                 3.731006820000000e+00 * tc[0] + 9.157986500000000e-04 * tc[1] -
                 7.205411033333332e-07 * tc[2] + 2.536484591666667e-10 * tc[3] -
                 3.730357810000000e-14 * tc[4];
  } else {
    // species 0: N2
    species[0] = -8.698115790000001e+02 * invT - 3.836719770000000e+00 -
                 2.811660730000000e+00 * tc[0] - 8.353367650000000e-04 * tc[1] +
                 1.133329046666667e-07 * tc[2] - 1.107344825000000e-11 * tc[3] +
                 5.138372100000000e-16 * tc[4];
  }

  // species with midpoint at T=700 kelvin
  if (T < 700) {
    // species 1: H2
    species[1] = -9.432366140000000e+02 * invT + 3.154275482000000e+00 -
                 2.642044380000000e+00 * tc[0] - 2.747646370000000e-03 * tc[1] +
                 2.119393900000000e-06 * tc[2] - 1.072909775000000e-09 * tc[3] +
                 2.350138745000000e-13 * tc[4];
    // species 3: O2
    species[3] = -1.060698270000000e+03 * invT - 7.728724999999992e-02 -
                 3.744039210000000e+00 * tc[0] + 1.398700735000000e-03 * tc[1] -
                 1.633537596666667e-06 * tc[2] + 8.604970250000000e-10 * tc[3] -
                 1.899656235000000e-13 * tc[4];
    // species 8: HO2
    species[8] = +3.158982340000000e+02 * invT - 2.824171830000000e+00 -
                 3.619942990000000e+00 * tc[0] - 5.290285250000000e-04 * tc[1] -
                 8.444649016666667e-07 * tc[2] + 5.281673016666667e-10 * tc[3] -
                 1.207986405000000e-13 * tc[4];
    // species 11: CH4
    species[11] = -1.025263080000000e+04 * invT + 1.021617083000000e+01 -
                  5.239673350000000e+00 * tc[0] +
                  7.341756150000000e-03 * tc[1] -
                  8.828878516666666e-06 * tc[2] +
                  4.513906850000000e-09 * tc[3] - 9.815928299999999e-13 * tc[4];
    // species 15: C
    species[15] = +8.544341050000000e+04 * invT - 2.023790810000000e+00 -
                  2.544951920000000e+00 * tc[0] +
                  1.238626405000000e-04 * tc[1] -
                  9.133637950000000e-08 * tc[2] +
                  4.571260416666667e-11 * tc[3] - 1.020586655000000e-14 * tc[4];
    // species 20: CH2O
    species[20] = -1.432701690000000e+04 * invT + 1.705923940000000e+00 -
                  4.326212960000000e+00 * tc[0] +
                  3.505759265000000e-03 * tc[1] -
                  5.252949350000000e-06 * tc[2] +
                  2.803988650000000e-09 * tc[3] - 6.172701150000000e-13 * tc[4];
  } else {
    // species 1: H2
    species[1] = -1.102830230000000e+03 * invT + 9.387257909999999e+00 -
                 3.781998810000000e+00 * tc[0] + 5.093662950000000e-04 * tc[1] -
                 2.070437216666667e-07 * tc[2] + 3.491765816666666e-11 * tc[3] -
                 2.377718965000000e-15 * tc[4];
    // species 3: O2
    species[3] = -9.313501480000000e+02 * invT - 5.129021440000001e+00 -
                 2.820124080000000e+00 * tc[0] - 1.241056785000000e-03 * tc[1] +
                 2.520034900000000e-07 * tc[2] - 3.737968341666667e-11 * tc[3] +
                 2.436528340000000e-15 * tc[4];
    // species 8: HO2
    species[8] = +3.993416090000000e+02 * invT - 6.083080840000000e+00 -
                 3.023918890000000e+00 * tc[0] - 2.231954535000000e-03 * tc[1] +
                 3.719108200000000e-07 * tc[2] - 5.105923325000000e-11 * tc[3] +
                 3.321331185000000e-15 * tc[4];
    // species 11: CH4
    species[11] = -9.589825010000000e+03 * invT - 1.566993109500000e+01 -
                  5.053464050000001e-01 * tc[0] -
                  6.184892250000000e-03 * tc[1] +
                  8.330132033333333e-07 * tc[2] -
                  8.699397083333333e-11 * tc[3] + 4.314487080000000e-15 * tc[4];
    // species 15: C
    species[15] = +8.545044220000000e+04 * invT - 2.298417230000000e+00 -
                  2.494725310000000e+00 * tc[0] -
                  1.964197380000000e-05 * tc[1] +
                  1.116691633333333e-08 * tc[2] -
                  3.098489116666667e-12 * tc[3] + 2.536534425000000e-16 * tc[4];
    // species 20: CH2O
    species[20] = -1.390801700000000e+04 * invT - 1.465825788000000e+01 -
                  1.333356520000000e+00 * tc[0] -
                  5.045259150000000e-03 * tc[1] +
                  8.549209366666667e-07 * tc[2] -
                  1.045210058333333e-10 * tc[3] + 5.981955450000001e-15 * tc[4];
  }

  // species with midpoint at T=860 kelvin
  if (T < 860) {
    // species 2: H
    species[2] = +2.547378660000000e+04 * invT + 2.944079458000000e+00 -
                 2.499505440000000e+00 * tc[0] - 1.495820230000000e-06 * tc[1] +
                 9.879329316666667e-10 * tc[2] - 4.065084708333333e-13 * tc[3] +
                 7.276966000000000e-17 * tc[4];
  } else {
    // species 2: H
    species[2] = +2.547364740000000e+04 * invT + 2.948672158000000e+00 -
                 2.500314930000000e+00 * tc[0] + 3.867034140000000e-07 * tc[1] -
                 1.065575576666667e-10 * tc[2] + 1.771264925000000e-14 * tc[3] -
                 1.222395955000000e-18 * tc[4];
  }

  // species with midpoint at T=720 kelvin
  if (T < 720) {
    // species 4: O
    species[4] = +2.913091180000000e+04 * invT + 1.013526520000000e+00 -
                 3.147992010000000e+00 * tc[0] + 1.555870325000000e-03 * tc[1] -
                 1.030229828333333e-06 * tc[2] + 4.698406650000000e-10 * tc[3] -
                 9.743300800000001e-14 * tc[4];
  } else {
    // species 4: O
    species[4] = +2.920615190000000e+04 * invT - 1.858093760000000e+00 -
                 2.625491430000000e+00 * tc[0] + 1.044798220000000e-04 * tc[1] -
                 2.231975766666666e-08 * tc[2] + 3.215632466666667e-12 * tc[3] -
                 2.194593445000000e-16 * tc[4];
  }

  // species with midpoint at T=1420 kelvin
  if (T < 1420) {
    // species 5: H2O
    species[5] = -3.028313140000000e+04 * invT + 4.356762201000000e+00 -
                 4.060611720000000e+00 * tc[0] + 4.329035945000000e-04 * tc[1] -
                 5.406825466666667e-07 * tc[2] + 1.502025658333333e-10 * tc[3] -
                 1.662416465000000e-14 * tc[4];
  } else {
    // species 5: H2O
    species[5] = -2.988756450000000e+04 * invT - 4.244140560000000e+00 -
                 2.667770750000000e+00 * tc[0] - 1.528844245000000e-03 * tc[1] +
                 1.500737351666667e-07 * tc[2] - 1.194679600000000e-11 * tc[3] +
                 5.042890850000000e-16 * tc[4];
  }

  // species with midpoint at T=1700 kelvin
  if (T < 1700) {
    // species 6: OH
    species[6] = +3.400050470000000e+03 * invT + 3.708071591000000e+00 -
                 3.913546310000000e+00 * tc[0] + 8.313796300000000e-04 * tc[1] -
                 3.848667150000000e-07 * tc[2] + 8.529959000000000e-11 * tc[3] -
                 7.941481450000000e-15 * tc[4];
  } else {
    // species 6: OH
    species[6] = +3.881107160000000e+03 * invT - 5.283514930000001e+00 -
                 2.498673690000000e+00 * tc[0] - 8.331763950000000e-04 * tc[1] +
                 1.047085860000000e-07 * tc[2] - 1.069556716666667e-11 * tc[3] +
                 5.286794700000000e-16 * tc[4];
  }

  // species with midpoint at T=1800 kelvin
  if (T < 1800) {
    // species 7: H2O2
    species[7] = -1.757847850000000e+04 * invT - 3.830946270000000e+00 -
                 3.347742240000000e+00 * tc[0] - 3.525027185000000e-03 * tc[1] +
                 6.408700100000000e-07 * tc[2] - 9.726721749999999e-11 * tc[3] +
                 7.380905249999999e-15 * tc[4];
    // species 13: CH2
    species[13] = +4.586456990000000e+04 * invT + 2.275363070000000e+00 -
                  3.764894600000000e+00 * tc[0] -
                  7.191959550000000e-04 * tc[1] -
                  7.926384616666667e-08 * tc[2] +
                  3.598238258333333e-11 * tc[3] - 3.791464370000000e-15 * tc[4];
    // species 17: CH3OH
    species[17] = -2.509625440000000e+04 * invT - 1.874599922100000e+01 -
                  8.473304790000000e-01 * tc[0] -
                  8.154345200000000e-03 * tc[1] +
                  1.413908268333333e-06 * tc[2] -
                  1.910869508333333e-10 * tc[3] + 1.299760065000000e-14 * tc[4];
  } else {
    // species 7: H2O2
    species[7] = -1.809002200000000e+04 * invT + 5.280508167000000e+00 -
                 4.768696390000000e+00 * tc[0] - 1.946189240000000e-03 * tc[1] +
                 2.023039150000000e-07 * tc[2] - 1.605127375000000e-11 * tc[3] +
                 6.129099500000000e-16 * tc[4];
    // species 13: CH2
    species[13] = +4.620734920000000e+04 * invT - 3.830116800000000e+00 -
                  2.812729720000000e+00 * tc[0] -
                  1.777156940000000e-03 * tc[1] +
                  2.146142050000000e-07 * tc[2] -
                  1.843947866666667e-11 * tc[3] + 7.436907349999999e-16 * tc[4];
    // species 17: CH3OH
    species[17] = -2.576934090000000e+04 * invT - 6.757190100000001e+00 -
                  2.717015300000000e+00 * tc[0] -
                  6.076917650000000e-03 * tc[1] +
                  8.368450516666667e-07 * tc[2] -
                  8.422339166666667e-11 * tc[3] + 4.092304115000000e-15 * tc[4];
  }

  // species with midpoint at T=960 kelvin
  if (T < 960) {
    // species 9: CO
    species[9] = -1.435755300000000e+04 * invT + 9.574709200000000e-01 -
                 3.757238910000000e+00 * tc[0] + 1.072326205000000e-03 * tc[1] -
                 9.034650083333333e-07 * tc[2] + 3.475216358333333e-10 * tc[3] -
                 5.595056350000000e-14 * tc[4];
  } else {
    // species 9: CO
    species[9] = -1.417233350000000e+04 * invT - 4.621881790000000e+00 -
                 2.792553810000000e+00 * tc[0] - 9.374344300000000e-04 * tc[1] +
                 1.432853210000000e-07 * tc[2] - 1.593333916666667e-11 * tc[3] +
                 8.392764300000001e-16 * tc[4];
  }

  // species with midpoint at T=1450 kelvin
  if (T < 1450) {
    // species 10: CO2
    species[10] = -4.836260670000000e+04 * invT - 7.761779929999999e+00 -
                  2.316843470000000e+00 * tc[0] -
                  4.613775180000000e-03 * tc[1] +
                  1.292756821666666e-06 * tc[2] -
                  2.735211333333333e-10 * tc[3] + 2.743612410000000e-14 * tc[4];
  } else {
    // species 10: CO2
    species[10] = -4.905626390000000e+04 * invT + 7.058529200000000e+00 -
                  4.708764680000000e+00 * tc[0] -
                  1.314573520000000e-03 * tc[1] +
                  1.551010770000000e-07 * tc[2] -
                  1.199107666666667e-11 * tc[3] + 3.812907069999999e-16 * tc[4];
  }

  // species with midpoint at T=1060 kelvin
  if (T < 1060) {
    // species 12: CH3
    species[12] = +1.643995160000000e+04 * invT + 1.069533540000000e+00 -
                  3.478293100000000e+00 * tc[0] -
                  1.773823865000000e-03 * tc[1] -
                  2.456807333333333e-07 * tc[2] +
                  1.619799625000000e-10 * tc[3] - 2.609606160000000e-14 * tc[4];
  } else {
    // species 12: CH3
    species[12] = +1.658628290000000e+04 * invT - 2.990947140000000e+00 -
                  2.788051040000000e+00 * tc[0] -
                  3.076167385000000e-03 * tc[1] +
                  3.686322483333333e-07 * tc[2] -
                  3.120022066666666e-11 * tc[3] + 1.240756745000000e-15 * tc[4];
  }

  // species with midpoint at T=970 kelvin
  if (T < 970) {
    // species 14: CH2(S)
    species[14] = +5.036694070000000e+04 * invT + 4.884856942000000e+00 -
                  4.181854340000000e+00 * tc[0] +
                  1.105671570000000e-03 * tc[1] -
                  1.285879235000000e-06 * tc[2] +
                  4.966253174999999e-10 * tc[3] - 7.915731400000000e-14 * tc[4];
  } else {
    // species 14: CH2(S)
    species[14] = +5.064290790000000e+04 * invT - 3.357120840000000e+00 -
                  2.759342990000000e+00 * tc[0] -
                  1.827341535000000e-03 * tc[1] +
                  2.259831900000000e-07 * tc[2] -
                  2.291503425000000e-11 * tc[3] + 1.183977360000000e-15 * tc[4];
  }

  // species with midpoint at T=1590 kelvin
  if (T < 1590) {
    // species 16: CH
    species[16] = +7.063124920000001e+04 * invT + 2.893235416000000e+00 -
                  3.772643320000000e+00 * tc[0] +
                  7.927367500000000e-04 * tc[1] -
                  4.725203966666667e-07 * tc[2] +
                  1.134550483333333e-10 * tc[3] - 1.119976660000000e-14 * tc[4];
  } else {
    // species 16: CH
    species[16] = +7.110594120000000e+04 * invT - 6.493359330000000e+00 -
                  2.279901280000000e+00 * tc[0] -
                  1.084926190000000e-03 * tc[1] +
                  1.179396475000000e-07 * tc[2] -
                  1.033112450000000e-11 * tc[3] + 4.781742575000000e-16 * tc[4];
  }

  // species with midpoint at T=1740 kelvin
  if (T < 1740) {
    // species 18: CH3O
    species[18] = +1.603161210000000e+03 * invT - 1.761045241500000e+01 -
                  8.896609850000000e-01 * tc[0] -
                  8.505988350000000e-03 * tc[1] +
                  1.896789183333333e-06 * tc[2] -
                  3.235674400000000e-10 * tc[3] + 2.664207395000000e-14 * tc[4];
  } else {
    // species 18: CH3O
    species[18] = -7.862522250000001e+01 * invT + 1.321411739000000e+01 -
                  5.722380620000000e+00 * tc[0] -
                  2.951138185000000e-03 * tc[1] +
                  3.005678666666667e-07 * tc[2] -
                  1.777791741666667e-11 * tc[3] + 2.809082045000000e-16 * tc[4];
  }

  // species with midpoint at T=1360 kelvin
  if (T < 1360) {
    // species 19: CH2OH
    species[19] = -3.302221060000000e+03 * invT - 9.917981910000000e+00 -
                  2.348215790000000e+00 * tc[0] -
                  6.980008400000000e-03 * tc[1] +
                  1.810536766666667e-06 * tc[2] -
                  3.854153466666667e-10 * tc[3] + 4.042495810000000e-14 * tc[4];
  } else {
    // species 19: CH2OH
    species[19] = -4.035841430000000e+03 * invT + 6.620590230000000e+00 -
                  5.045349500000000e+00 * tc[0] -
                  3.013635295000000e-03 * tc[1] +
                  3.523113683333333e-07 * tc[2] -
                  2.800715875000000e-11 * tc[3] + 1.004937410000000e-15 * tc[4];
  }

  // species with midpoint at T=770 kelvin
  if (T < 770) {
    // species 21: HCO
    species[21] = +3.836363920000000e+03 * invT - 1.652479099999997e-01 -
                  4.034839790000000e+00 * tc[0] +
                  1.079184320000000e-03 * tc[1] -
                  1.970564583333333e-06 * tc[2] +
                  9.871617166666666e-10 * tc[3] - 2.002969770000000e-13 * tc[4];
  } else {
    // species 21: HCO
    species[21] = +4.057253300000000e+03 * invT - 8.144600120000000e+00 -
                  2.600493180000000e+00 * tc[0] -
                  2.646391290000000e-03 * tc[1] +
                  4.486403516666666e-07 * tc[2] -
                  6.011314983333334e-11 * tc[3] + 3.717607045000000e-15 * tc[4];
  }

  // species with midpoint at T=3654.18 kelvin
  if (T < 3654.18) {
    // species 22: CHO+
    species[22] = +9.913002200000000e+04 * invT - 3.625749800000000e+00 -
                  2.876276800000000e+00 * tc[0] -
                  2.474301000000000e-03 * tc[1] +
                  3.847551000000000e-07 * tc[2] -
                  4.260146000000000e-11 * tc[3] + 2.176031550000000e-15 * tc[4];
  } else {
    // species 22: CHO+
    species[22] = +9.612933300000000e+04 * invT + 2.500523650000000e+01 -
                  6.915726500000000e+00 * tc[0] -
                  8.315553500000000e-05 * tc[1] +
                  3.079390500000000e-09 * tc[2] -
                  7.571369083333333e-14 * tc[3] + 8.227220500000000e-19 * tc[4];
  }

  // species with midpoint at T=5430.9 kelvin
  if (T < 5430.9) {
    // species 23: H3O+
    species[23] = +7.089026800000000e+04 * invT - 3.672298500000000e+00 -
                  2.688523200000000e+00 * tc[0] -
                  2.722682700000000e-03 * tc[1] +
                  2.842594333333333e-07 * tc[2] -
                  2.069441083333333e-11 * tc[3] + 6.884720000000001e-16 * tc[4];
  } else {
    // species 23: H3O+
    species[23] = +6.448074900000000e+04 * invT + 4.282419930000000e+01 -
                  8.950880300000000e+00 * tc[0] -
                  1.264083400000000e-04 * tc[1] +
                  4.099066166666667e-09 * tc[2] -
                  9.065265000000000e-14 * tc[3] + 9.045418500000001e-19 * tc[4];
  }

  // species with midpoint at T=2008.71 kelvin
  if (T < 2008.71) {
    // species 24: O2-
    species[24] = -6.807479300000000e+03 * invT - 3.658730200000000e+00 -
                  3.102171800000000e+00 * tc[0] -
                  1.399043750000000e-03 * tc[1] +
                  3.775187666666667e-07 * tc[2] -
                  7.243043083333333e-11 * tc[3] + 6.360942000000000e-15 * tc[4];
  } else {
    // species 24: O2-
    species[24] = -7.242625200000000e+03 * invT + 3.783289729999999e+00 -
                  4.259286700000000e+00 * tc[0] -
                  1.123403600000000e-04 * tc[1] +
                  8.566325833333332e-09 * tc[2] -
                  6.128831500000000e-13 * tc[3] + 1.927932600000000e-17 * tc[4];
  }

  // species with midpoint at T=27300 kelvin
  if (T < 27300) {
    // species 25: E
    species[25] = -7.453750000000000e+02 * invT + 1.423397500000000e+01 -
                  2.500000000000000e+00 * tc[0] +
                  4.075130400000000e-18 * tc[1] -
                  2.171644500000000e-22 * tc[2] +
                  6.054027750000001e-27 * tc[3] - 6.533123499999999e-32 * tc[4];
  } else {
    // species 25: E
    species[25] = -7.453750000000000e+02 * invT + 1.423397500000000e+01 -
                  2.500000000000000e+00 * tc[0] -
                  0.000000000000000e+00 * tc[1] -
                  0.000000000000000e+00 * tc[2] -
                  0.000000000000000e+00 * tc[3] - 0.000000000000000e+00 * tc[4];
  }
}

// compute the a/(RT) at the given temperature
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
helmholtz(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1 / T;

  // species with midpoint at T=1050 kelvin
  if (T < 1050) {
    // species 0: N2
    species[0] = -1.06287426e+03 * invT + 5.62794840e-01 -
                 3.73100682e+00 * tc[0] + 9.15798650e-04 * tc[1] -
                 7.20541103e-07 * tc[2] + 2.53648459e-10 * tc[3] -
                 3.73035781e-14 * tc[4];
  } else {
    // species 0: N2
    species[0] = -8.69811579e+02 * invT - 4.83671977e+00 -
                 2.81166073e+00 * tc[0] - 8.35336765e-04 * tc[1] +
                 1.13332905e-07 * tc[2] - 1.10734483e-11 * tc[3] +
                 5.13837210e-16 * tc[4];
  }

  // species with midpoint at T=700 kelvin
  if (T < 700) {
    // species 1: H2
    species[1] = -9.43236614e+02 * invT + 2.15427548e+00 -
                 2.64204438e+00 * tc[0] - 2.74764637e-03 * tc[1] +
                 2.11939390e-06 * tc[2] - 1.07290977e-09 * tc[3] +
                 2.35013875e-13 * tc[4];
    // species 3: O2
    species[3] = -1.06069827e+03 * invT - 1.07728725e+00 -
                 3.74403921e+00 * tc[0] + 1.39870073e-03 * tc[1] -
                 1.63353760e-06 * tc[2] + 8.60497025e-10 * tc[3] -
                 1.89965623e-13 * tc[4];
    // species 8: HO2
    species[8] = +3.15898234e+02 * invT - 3.82417183e+00 -
                 3.61994299e+00 * tc[0] - 5.29028525e-04 * tc[1] -
                 8.44464902e-07 * tc[2] + 5.28167302e-10 * tc[3] -
                 1.20798641e-13 * tc[4];
    // species 11: CH4
    species[11] = -1.02526308e+04 * invT + 9.21617083e+00 -
                  5.23967335e+00 * tc[0] + 7.34175615e-03 * tc[1] -
                  8.82887852e-06 * tc[2] + 4.51390685e-09 * tc[3] -
                  9.81592830e-13 * tc[4];
    // species 15: C
    species[15] = +8.54434105e+04 * invT - 3.02379081e+00 -
                  2.54495192e+00 * tc[0] + 1.23862640e-04 * tc[1] -
                  9.13363795e-08 * tc[2] + 4.57126042e-11 * tc[3] -
                  1.02058665e-14 * tc[4];
    // species 20: CH2O
    species[20] = -1.43270169e+04 * invT + 7.05923940e-01 -
                  4.32621296e+00 * tc[0] + 3.50575927e-03 * tc[1] -
                  5.25294935e-06 * tc[2] + 2.80398865e-09 * tc[3] -
                  6.17270115e-13 * tc[4];
  } else {
    // species 1: H2
    species[1] = -1.10283023e+03 * invT + 8.38725791e+00 -
                 3.78199881e+00 * tc[0] + 5.09366295e-04 * tc[1] -
                 2.07043722e-07 * tc[2] + 3.49176582e-11 * tc[3] -
                 2.37771897e-15 * tc[4];
    // species 3: O2
    species[3] = -9.31350148e+02 * invT - 6.12902144e+00 -
                 2.82012408e+00 * tc[0] - 1.24105679e-03 * tc[1] +
                 2.52003490e-07 * tc[2] - 3.73796834e-11 * tc[3] +
                 2.43652834e-15 * tc[4];
    // species 8: HO2
    species[8] = +3.99341609e+02 * invT - 7.08308084e+00 -
                 3.02391889e+00 * tc[0] - 2.23195453e-03 * tc[1] +
                 3.71910820e-07 * tc[2] - 5.10592333e-11 * tc[3] +
                 3.32133119e-15 * tc[4];
    // species 11: CH4
    species[11] = -9.58982501e+03 * invT - 1.66699311e+01 -
                  5.05346405e-01 * tc[0] - 6.18489225e-03 * tc[1] +
                  8.33013203e-07 * tc[2] - 8.69939708e-11 * tc[3] +
                  4.31448708e-15 * tc[4];
    // species 15: C
    species[15] = +8.54504422e+04 * invT - 3.29841723e+00 -
                  2.49472531e+00 * tc[0] - 1.96419738e-05 * tc[1] +
                  1.11669163e-08 * tc[2] - 3.09848912e-12 * tc[3] +
                  2.53653442e-16 * tc[4];
    // species 20: CH2O
    species[20] = -1.39080170e+04 * invT - 1.56582579e+01 -
                  1.33335652e+00 * tc[0] - 5.04525915e-03 * tc[1] +
                  8.54920937e-07 * tc[2] - 1.04521006e-10 * tc[3] +
                  5.98195545e-15 * tc[4];
  }

  // species with midpoint at T=860 kelvin
  if (T < 860) {
    // species 2: H
    species[2] = +2.54737866e+04 * invT + 1.94407946e+00 -
                 2.49950544e+00 * tc[0] - 1.49582023e-06 * tc[1] +
                 9.87932932e-10 * tc[2] - 4.06508471e-13 * tc[3] +
                 7.27696600e-17 * tc[4];
  } else {
    // species 2: H
    species[2] = +2.54736474e+04 * invT + 1.94867216e+00 -
                 2.50031493e+00 * tc[0] + 3.86703414e-07 * tc[1] -
                 1.06557558e-10 * tc[2] + 1.77126492e-14 * tc[3] -
                 1.22239596e-18 * tc[4];
  }

  // species with midpoint at T=720 kelvin
  if (T < 720) {
    // species 4: O
    species[4] = +2.91309118e+04 * invT + 1.35265200e-02 -
                 3.14799201e+00 * tc[0] + 1.55587032e-03 * tc[1] -
                 1.03022983e-06 * tc[2] + 4.69840665e-10 * tc[3] -
                 9.74330080e-14 * tc[4];
  } else {
    // species 4: O
    species[4] = +2.92061519e+04 * invT - 2.85809376e+00 -
                 2.62549143e+00 * tc[0] + 1.04479822e-04 * tc[1] -
                 2.23197577e-08 * tc[2] + 3.21563247e-12 * tc[3] -
                 2.19459344e-16 * tc[4];
  }

  // species with midpoint at T=1420 kelvin
  if (T < 1420) {
    // species 5: H2O
    species[5] = -3.02831314e+04 * invT + 3.35676220e+00 -
                 4.06061172e+00 * tc[0] + 4.32903595e-04 * tc[1] -
                 5.40682547e-07 * tc[2] + 1.50202566e-10 * tc[3] -
                 1.66241647e-14 * tc[4];
  } else {
    // species 5: H2O
    species[5] = -2.98875645e+04 * invT - 5.24414056e+00 -
                 2.66777075e+00 * tc[0] - 1.52884425e-03 * tc[1] +
                 1.50073735e-07 * tc[2] - 1.19467960e-11 * tc[3] +
                 5.04289085e-16 * tc[4];
  }

  // species with midpoint at T=1700 kelvin
  if (T < 1700) {
    // species 6: OH
    species[6] = +3.40005047e+03 * invT + 2.70807159e+00 -
                 3.91354631e+00 * tc[0] + 8.31379630e-04 * tc[1] -
                 3.84866715e-07 * tc[2] + 8.52995900e-11 * tc[3] -
                 7.94148145e-15 * tc[4];
  } else {
    // species 6: OH
    species[6] = +3.88110716e+03 * invT - 6.28351493e+00 -
                 2.49867369e+00 * tc[0] - 8.33176395e-04 * tc[1] +
                 1.04708586e-07 * tc[2] - 1.06955672e-11 * tc[3] +
                 5.28679470e-16 * tc[4];
  }

  // species with midpoint at T=1800 kelvin
  if (T < 1800) {
    // species 7: H2O2
    species[7] = -1.75784785e+04 * invT - 4.83094627e+00 -
                 3.34774224e+00 * tc[0] - 3.52502719e-03 * tc[1] +
                 6.40870010e-07 * tc[2] - 9.72672175e-11 * tc[3] +
                 7.38090525e-15 * tc[4];
    // species 13: CH2
    species[13] = +4.58645699e+04 * invT + 1.27536307e+00 -
                  3.76489460e+00 * tc[0] - 7.19195955e-04 * tc[1] -
                  7.92638462e-08 * tc[2] + 3.59823826e-11 * tc[3] -
                  3.79146437e-15 * tc[4];
    // species 17: CH3OH
    species[17] = -2.50962544e+04 * invT - 1.97459992e+01 -
                  8.47330479e-01 * tc[0] - 8.15434520e-03 * tc[1] +
                  1.41390827e-06 * tc[2] - 1.91086951e-10 * tc[3] +
                  1.29976006e-14 * tc[4];
  } else {
    // species 7: H2O2
    species[7] = -1.80900220e+04 * invT + 4.28050817e+00 -
                 4.76869639e+00 * tc[0] - 1.94618924e-03 * tc[1] +
                 2.02303915e-07 * tc[2] - 1.60512737e-11 * tc[3] +
                 6.12909950e-16 * tc[4];
    // species 13: CH2
    species[13] = +4.62073492e+04 * invT - 4.83011680e+00 -
                  2.81272972e+00 * tc[0] - 1.77715694e-03 * tc[1] +
                  2.14614205e-07 * tc[2] - 1.84394787e-11 * tc[3] +
                  7.43690735e-16 * tc[4];
    // species 17: CH3OH
    species[17] = -2.57693409e+04 * invT - 7.75719010e+00 -
                  2.71701530e+00 * tc[0] - 6.07691765e-03 * tc[1] +
                  8.36845052e-07 * tc[2] - 8.42233917e-11 * tc[3] +
                  4.09230411e-15 * tc[4];
  }

  // species with midpoint at T=960 kelvin
  if (T < 960) {
    // species 9: CO
    species[9] = -1.43575530e+04 * invT - 4.25290800e-02 -
                 3.75723891e+00 * tc[0] + 1.07232621e-03 * tc[1] -
                 9.03465008e-07 * tc[2] + 3.47521636e-10 * tc[3] -
                 5.59505635e-14 * tc[4];
  } else {
    // species 9: CO
    species[9] = -1.41723335e+04 * invT - 5.62188179e+00 -
                 2.79255381e+00 * tc[0] - 9.37434430e-04 * tc[1] +
                 1.43285321e-07 * tc[2] - 1.59333392e-11 * tc[3] +
                 8.39276430e-16 * tc[4];
  }

  // species with midpoint at T=1450 kelvin
  if (T < 1450) {
    // species 10: CO2
    species[10] = -4.83626067e+04 * invT - 8.76177993e+00 -
                  2.31684347e+00 * tc[0] - 4.61377518e-03 * tc[1] +
                  1.29275682e-06 * tc[2] - 2.73521133e-10 * tc[3] +
                  2.74361241e-14 * tc[4];
  } else {
    // species 10: CO2
    species[10] = -4.90562639e+04 * invT + 6.05852920e+00 -
                  4.70876468e+00 * tc[0] - 1.31457352e-03 * tc[1] +
                  1.55101077e-07 * tc[2] - 1.19910767e-11 * tc[3] +
                  3.81290707e-16 * tc[4];
  }

  // species with midpoint at T=1060 kelvin
  if (T < 1060) {
    // species 12: CH3
    species[12] = +1.64399516e+04 * invT + 6.95335400e-02 -
                  3.47829310e+00 * tc[0] - 1.77382386e-03 * tc[1] -
                  2.45680733e-07 * tc[2] + 1.61979963e-10 * tc[3] -
                  2.60960616e-14 * tc[4];
  } else {
    // species 12: CH3
    species[12] = +1.65862829e+04 * invT - 3.99094714e+00 -
                  2.78805104e+00 * tc[0] - 3.07616739e-03 * tc[1] +
                  3.68632248e-07 * tc[2] - 3.12002207e-11 * tc[3] +
                  1.24075675e-15 * tc[4];
  }

  // species with midpoint at T=970 kelvin
  if (T < 970) {
    // species 14: CH2(S)
    species[14] = +5.03669407e+04 * invT + 3.88485694e+00 -
                  4.18185434e+00 * tc[0] + 1.10567157e-03 * tc[1] -
                  1.28587924e-06 * tc[2] + 4.96625317e-10 * tc[3] -
                  7.91573140e-14 * tc[4];
  } else {
    // species 14: CH2(S)
    species[14] = +5.06429079e+04 * invT - 4.35712084e+00 -
                  2.75934299e+00 * tc[0] - 1.82734153e-03 * tc[1] +
                  2.25983190e-07 * tc[2] - 2.29150342e-11 * tc[3] +
                  1.18397736e-15 * tc[4];
  }

  // species with midpoint at T=1590 kelvin
  if (T < 1590) {
    // species 16: CH
    species[16] = +7.06312492e+04 * invT + 1.89323542e+00 -
                  3.77264332e+00 * tc[0] + 7.92736750e-04 * tc[1] -
                  4.72520397e-07 * tc[2] + 1.13455048e-10 * tc[3] -
                  1.11997666e-14 * tc[4];
  } else {
    // species 16: CH
    species[16] = +7.11059412e+04 * invT - 7.49335933e+00 -
                  2.27990128e+00 * tc[0] - 1.08492619e-03 * tc[1] +
                  1.17939648e-07 * tc[2] - 1.03311245e-11 * tc[3] +
                  4.78174257e-16 * tc[4];
  }

  // species with midpoint at T=1740 kelvin
  if (T < 1740) {
    // species 18: CH3O
    species[18] = +1.60316121e+03 * invT - 1.86104524e+01 -
                  8.89660985e-01 * tc[0] - 8.50598835e-03 * tc[1] +
                  1.89678918e-06 * tc[2] - 3.23567440e-10 * tc[3] +
                  2.66420740e-14 * tc[4];
  } else {
    // species 18: CH3O
    species[18] = -7.86252225e+01 * invT + 1.22141174e+01 -
                  5.72238062e+00 * tc[0] - 2.95113818e-03 * tc[1] +
                  3.00567867e-07 * tc[2] - 1.77779174e-11 * tc[3] +
                  2.80908204e-16 * tc[4];
  }

  // species with midpoint at T=1360 kelvin
  if (T < 1360) {
    // species 19: CH2OH
    species[19] = -3.30222106e+03 * invT - 1.09179819e+01 -
                  2.34821579e+00 * tc[0] - 6.98000840e-03 * tc[1] +
                  1.81053677e-06 * tc[2] - 3.85415347e-10 * tc[3] +
                  4.04249581e-14 * tc[4];
  } else {
    // species 19: CH2OH
    species[19] = -4.03584143e+03 * invT + 5.62059023e+00 -
                  5.04534950e+00 * tc[0] - 3.01363529e-03 * tc[1] +
                  3.52311368e-07 * tc[2] - 2.80071587e-11 * tc[3] +
                  1.00493741e-15 * tc[4];
  }

  // species with midpoint at T=770 kelvin
  if (T < 770) {
    // species 21: HCO
    species[21] = +3.83636392e+03 * invT - 1.16524791e+00 -
                  4.03483979e+00 * tc[0] + 1.07918432e-03 * tc[1] -
                  1.97056458e-06 * tc[2] + 9.87161717e-10 * tc[3] -
                  2.00296977e-13 * tc[4];
  } else {
    // species 21: HCO
    species[21] = +4.05725330e+03 * invT - 9.14460012e+00 -
                  2.60049318e+00 * tc[0] - 2.64639129e-03 * tc[1] +
                  4.48640352e-07 * tc[2] - 6.01131498e-11 * tc[3] +
                  3.71760704e-15 * tc[4];
  }

  // species with midpoint at T=3654.18 kelvin
  if (T < 3654.18) {
    // species 22: CHO+
    species[22] = +9.91300220e+04 * invT - 4.62574980e+00 -
                  2.87627680e+00 * tc[0] - 2.47430100e-03 * tc[1] +
                  3.84755100e-07 * tc[2] - 4.26014600e-11 * tc[3] +
                  2.17603155e-15 * tc[4];
  } else {
    // species 22: CHO+
    species[22] = +9.61293330e+04 * invT + 2.40052365e+01 -
                  6.91572650e+00 * tc[0] - 8.31555350e-05 * tc[1] +
                  3.07939050e-09 * tc[2] - 7.57136908e-14 * tc[3] +
                  8.22722050e-19 * tc[4];
  }

  // species with midpoint at T=5430.9 kelvin
  if (T < 5430.9) {
    // species 23: H3O+
    species[23] = +7.08902680e+04 * invT - 4.67229850e+00 -
                  2.68852320e+00 * tc[0] - 2.72268270e-03 * tc[1] +
                  2.84259433e-07 * tc[2] - 2.06944108e-11 * tc[3] +
                  6.88472000e-16 * tc[4];
  } else {
    // species 23: H3O+
    species[23] = +6.44807490e+04 * invT + 4.18241993e+01 -
                  8.95088030e+00 * tc[0] - 1.26408340e-04 * tc[1] +
                  4.09906617e-09 * tc[2] - 9.06526500e-14 * tc[3] +
                  9.04541850e-19 * tc[4];
  }

  // species with midpoint at T=2008.71 kelvin
  if (T < 2008.71) {
    // species 24: O2-
    species[24] = -6.80747930e+03 * invT - 4.65873020e+00 -
                  3.10217180e+00 * tc[0] - 1.39904375e-03 * tc[1] +
                  3.77518767e-07 * tc[2] - 7.24304308e-11 * tc[3] +
                  6.36094200e-15 * tc[4];
  } else {
    // species 24: O2-
    species[24] = -7.24262520e+03 * invT + 2.78328973e+00 -
                  4.25928670e+00 * tc[0] - 1.12340360e-04 * tc[1] +
                  8.56632583e-09 * tc[2] - 6.12883150e-13 * tc[3] +
                  1.92793260e-17 * tc[4];
  }

  // species with midpoint at T=27300 kelvin
  if (T < 27300) {
    // species 25: E
    species[25] = -7.45375000e+02 * invT + 1.32339750e+01 -
                  2.50000000e+00 * tc[0] + 4.07513040e-18 * tc[1] -
                  2.17164450e-22 * tc[2] + 6.05402775e-27 * tc[3] -
                  6.53312350e-32 * tc[4];
  } else {
    // species 25: E
    species[25] = -7.45375000e+02 * invT + 1.32339750e+01 -
                  2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                  0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                  0.00000000e+00 * tc[4];
  }
}

// compute the e/(RT) at the given temperature
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesInternalEnergy(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1 / T;

  // species with midpoint at T=1050 kelvin
  if (T < 1050) {
    // species 0: N2
    species[0] = +2.73100682e+00 - 9.15798650e-04 * tc[1] +
                 1.44108221e-06 * tc[2] - 7.60945377e-10 * tc[3] +
                 1.49214312e-13 * tc[4] - 1.06287426e+03 * invT;
  } else {
    // species 0: N2
    species[0] = +1.81166073e+00 + 8.35336765e-04 * tc[1] -
                 2.26665809e-07 * tc[2] + 3.32203448e-11 * tc[3] -
                 2.05534884e-15 * tc[4] - 8.69811579e+02 * invT;
  }

  // species with midpoint at T=700 kelvin
  if (T < 700) {
    // species 1: H2
    species[1] = +1.64204438e+00 + 2.74764637e-03 * tc[1] -
                 4.23878780e-06 * tc[2] + 3.21872932e-09 * tc[3] -
                 9.40055498e-13 * tc[4] - 9.43236614e+02 * invT;
    // species 3: O2
    species[3] = +2.74403921e+00 - 1.39870073e-03 * tc[1] +
                 3.26707519e-06 * tc[2] - 2.58149107e-09 * tc[3] +
                 7.59862494e-13 * tc[4] - 1.06069827e+03 * invT;
    // species 8: HO2
    species[8] = +2.61994299e+00 + 5.29028525e-04 * tc[1] +
                 1.68892980e-06 * tc[2] - 1.58450190e-09 * tc[3] +
                 4.83194562e-13 * tc[4] + 3.15898234e+02 * invT;
    // species 11: CH4
    species[11] = +4.23967335e+00 - 7.34175615e-03 * tc[1] +
                  1.76577570e-05 * tc[2] - 1.35417205e-08 * tc[3] +
                  3.92637132e-12 * tc[4] - 1.02526308e+04 * invT;
    // species 15: C
    species[15] = +1.54495192e+00 - 1.23862640e-04 * tc[1] +
                  1.82672759e-07 * tc[2] - 1.37137813e-10 * tc[3] +
                  4.08234662e-14 * tc[4] + 8.54434105e+04 * invT;
    // species 20: CH2O
    species[20] = +3.32621296e+00 - 3.50575927e-03 * tc[1] +
                  1.05058987e-05 * tc[2] - 8.41196595e-09 * tc[3] +
                  2.46908046e-12 * tc[4] - 1.43270169e+04 * invT;
  } else {
    // species 1: H2
    species[1] = +2.78199881e+00 - 5.09366295e-04 * tc[1] +
                 4.14087443e-07 * tc[2] - 1.04752974e-10 * tc[3] +
                 9.51087586e-15 * tc[4] - 1.10283023e+03 * invT;
    // species 3: O2
    species[3] = +1.82012408e+00 + 1.24105679e-03 * tc[1] -
                 5.04006980e-07 * tc[2] + 1.12139050e-10 * tc[3] -
                 9.74611336e-15 * tc[4] - 9.31350148e+02 * invT;
    // species 8: HO2
    species[8] = +2.02391889e+00 + 2.23195453e-03 * tc[1] -
                 7.43821640e-07 * tc[2] + 1.53177700e-10 * tc[3] -
                 1.32853247e-14 * tc[4] + 3.99341609e+02 * invT;
    // species 11: CH4
    species[11] = -4.94653595e-01 + 6.18489225e-03 * tc[1] -
                  1.66602641e-06 * tc[2] + 2.60981913e-10 * tc[3] -
                  1.72579483e-14 * tc[4] - 9.58982501e+03 * invT;
    // species 15: C
    species[15] = +1.49472531e+00 + 1.96419738e-05 * tc[1] -
                  2.23338327e-08 * tc[2] + 9.29546735e-12 * tc[3] -
                  1.01461377e-15 * tc[4] + 8.54504422e+04 * invT;
    // species 20: CH2O
    species[20] = +3.33356520e-01 + 5.04525915e-03 * tc[1] -
                  1.70984187e-06 * tc[2] + 3.13563018e-10 * tc[3] -
                  2.39278218e-14 * tc[4] - 1.39080170e+04 * invT;
  }

  // species with midpoint at T=860 kelvin
  if (T < 860) {
    // species 2: H
    species[2] = +1.49950544e+00 + 1.49582023e-06 * tc[1] -
                 1.97586586e-09 * tc[2] + 1.21952541e-12 * tc[3] -
                 2.91078640e-16 * tc[4] + 2.54737866e+04 * invT;
  } else {
    // species 2: H
    species[2] = +1.50031493e+00 - 3.86703414e-07 * tc[1] +
                 2.13115115e-10 * tc[2] - 5.31379477e-14 * tc[3] +
                 4.88958382e-18 * tc[4] + 2.54736474e+04 * invT;
  }

  // species with midpoint at T=720 kelvin
  if (T < 720) {
    // species 4: O
    species[4] = +2.14799201e+00 - 1.55587032e-03 * tc[1] +
                 2.06045966e-06 * tc[2] - 1.40952200e-09 * tc[3] +
                 3.89732032e-13 * tc[4] + 2.91309118e+04 * invT;
  } else {
    // species 4: O
    species[4] = +1.62549143e+00 - 1.04479822e-04 * tc[1] +
                 4.46395153e-08 * tc[2] - 9.64689740e-12 * tc[3] +
                 8.77837378e-16 * tc[4] + 2.92061519e+04 * invT;
  }

  // species with midpoint at T=1420 kelvin
  if (T < 1420) {
    // species 5: H2O
    species[5] = +3.06061172e+00 - 4.32903595e-04 * tc[1] +
                 1.08136509e-06 * tc[2] - 4.50607697e-10 * tc[3] +
                 6.64966586e-14 * tc[4] - 3.02831314e+04 * invT;
  } else {
    // species 5: H2O
    species[5] = +1.66777075e+00 + 1.52884425e-03 * tc[1] -
                 3.00147470e-07 * tc[2] + 3.58403880e-11 * tc[3] -
                 2.01715634e-15 * tc[4] - 2.98875645e+04 * invT;
  }

  // species with midpoint at T=1700 kelvin
  if (T < 1700) {
    // species 6: OH
    species[6] = +2.91354631e+00 - 8.31379630e-04 * tc[1] +
                 7.69733430e-07 * tc[2] - 2.55898770e-10 * tc[3] +
                 3.17659258e-14 * tc[4] + 3.40005047e+03 * invT;
  } else {
    // species 6: OH
    species[6] = +1.49867369e+00 + 8.33176395e-04 * tc[1] -
                 2.09417172e-07 * tc[2] + 3.20867015e-11 * tc[3] -
                 2.11471788e-15 * tc[4] + 3.88110716e+03 * invT;
  }

  // species with midpoint at T=1800 kelvin
  if (T < 1800) {
    // species 7: H2O2
    species[7] = +2.34774224e+00 + 3.52502719e-03 * tc[1] -
                 1.28174002e-06 * tc[2] + 2.91801652e-10 * tc[3] -
                 2.95236210e-14 * tc[4] - 1.75784785e+04 * invT;
    // species 13: CH2
    species[13] = +2.76489460e+00 + 7.19195955e-04 * tc[1] +
                  1.58527692e-07 * tc[2] - 1.07947148e-10 * tc[3] +
                  1.51658575e-14 * tc[4] + 4.58645699e+04 * invT;
    // species 17: CH3OH
    species[17] = -1.52669521e-01 + 8.15434520e-03 * tc[1] -
                  2.82781654e-06 * tc[2] + 5.73260852e-10 * tc[3] -
                  5.19904026e-14 * tc[4] - 2.50962544e+04 * invT;
  } else {
    // species 7: H2O2
    species[7] = +3.76869639e+00 + 1.94618924e-03 * tc[1] -
                 4.04607830e-07 * tc[2] + 4.81538212e-11 * tc[3] -
                 2.45163980e-15 * tc[4] - 1.80900220e+04 * invT;
    // species 13: CH2
    species[13] = +1.81272972e+00 + 1.77715694e-03 * tc[1] -
                  4.29228410e-07 * tc[2] + 5.53184360e-11 * tc[3] -
                  2.97476294e-15 * tc[4] + 4.62073492e+04 * invT;
    // species 17: CH3OH
    species[17] = +1.71701530e+00 + 6.07691765e-03 * tc[1] -
                  1.67369010e-06 * tc[2] + 2.52670175e-10 * tc[3] -
                  1.63692165e-14 * tc[4] - 2.57693409e+04 * invT;
  }

  // species with midpoint at T=960 kelvin
  if (T < 960) {
    // species 9: CO
    species[9] = +2.75723891e+00 - 1.07232621e-03 * tc[1] +
                 1.80693002e-06 * tc[2] - 1.04256491e-09 * tc[3] +
                 2.23802254e-13 * tc[4] - 1.43575530e+04 * invT;
  } else {
    // species 9: CO
    species[9] = +1.79255381e+00 + 9.37434430e-04 * tc[1] -
                 2.86570642e-07 * tc[2] + 4.78000175e-11 * tc[3] -
                 3.35710572e-15 * tc[4] - 1.41723335e+04 * invT;
  }

  // species with midpoint at T=1450 kelvin
  if (T < 1450) {
    // species 10: CO2
    species[10] = +1.31684347e+00 + 4.61377518e-03 * tc[1] -
                  2.58551364e-06 * tc[2] + 8.20563400e-10 * tc[3] -
                  1.09744496e-13 * tc[4] - 4.83626067e+04 * invT;
  } else {
    // species 10: CO2
    species[10] = +3.70876468e+00 + 1.31457352e-03 * tc[1] -
                  3.10202154e-07 * tc[2] + 3.59732300e-11 * tc[3] -
                  1.52516283e-15 * tc[4] - 4.90562639e+04 * invT;
  }

  // species with midpoint at T=1060 kelvin
  if (T < 1060) {
    // species 12: CH3
    species[12] = +2.47829310e+00 + 1.77382386e-03 * tc[1] +
                  4.91361467e-07 * tc[2] - 4.85939888e-10 * tc[3] +
                  1.04384246e-13 * tc[4] + 1.64399516e+04 * invT;
  } else {
    // species 12: CH3
    species[12] = +1.78805104e+00 + 3.07616739e-03 * tc[1] -
                  7.37264497e-07 * tc[2] + 9.36006620e-11 * tc[3] -
                  4.96302698e-15 * tc[4] + 1.65862829e+04 * invT;
  }

  // species with midpoint at T=970 kelvin
  if (T < 970) {
    // species 14: CH2(S)
    species[14] = +3.18185434e+00 - 1.10567157e-03 * tc[1] +
                  2.57175847e-06 * tc[2] - 1.48987595e-09 * tc[3] +
                  3.16629256e-13 * tc[4] + 5.03669407e+04 * invT;
  } else {
    // species 14: CH2(S)
    species[14] = +1.75934299e+00 + 1.82734153e-03 * tc[1] -
                  4.51966380e-07 * tc[2] + 6.87451027e-11 * tc[3] -
                  4.73590944e-15 * tc[4] + 5.06429079e+04 * invT;
  }

  // species with midpoint at T=1590 kelvin
  if (T < 1590) {
    // species 16: CH
    species[16] = +2.77264332e+00 - 7.92736750e-04 * tc[1] +
                  9.45040793e-07 * tc[2] - 3.40365145e-10 * tc[3] +
                  4.47990664e-14 * tc[4] + 7.06312492e+04 * invT;
  } else {
    // species 16: CH
    species[16] = +1.27990128e+00 + 1.08492619e-03 * tc[1] -
                  2.35879295e-07 * tc[2] + 3.09933735e-11 * tc[3] -
                  1.91269703e-15 * tc[4] + 7.11059412e+04 * invT;
  }

  // species with midpoint at T=1740 kelvin
  if (T < 1740) {
    // species 18: CH3O
    species[18] = -1.10339015e-01 + 8.50598835e-03 * tc[1] -
                  3.79357837e-06 * tc[2] + 9.70702320e-10 * tc[3] -
                  1.06568296e-13 * tc[4] + 1.60316121e+03 * invT;
  } else {
    // species 18: CH3O
    species[18] = +4.72238062e+00 + 2.95113818e-03 * tc[1] -
                  6.01135733e-07 * tc[2] + 5.33337522e-11 * tc[3] -
                  1.12363282e-15 * tc[4] - 7.86252225e+01 * invT;
  }

  // species with midpoint at T=1360 kelvin
  if (T < 1360) {
    // species 19: CH2OH
    species[19] = +1.34821579e+00 + 6.98000840e-03 * tc[1] -
                  3.62107353e-06 * tc[2] + 1.15624604e-09 * tc[3] -
                  1.61699832e-13 * tc[4] - 3.30222106e+03 * invT;
  } else {
    // species 19: CH2OH
    species[19] = +4.04534950e+00 + 3.01363529e-03 * tc[1] -
                  7.04622737e-07 * tc[2] + 8.40214762e-11 * tc[3] -
                  4.01974964e-15 * tc[4] - 4.03584143e+03 * invT;
  }

  // species with midpoint at T=770 kelvin
  if (T < 770) {
    // species 21: HCO
    species[21] = +3.03483979e+00 - 1.07918432e-03 * tc[1] +
                  3.94112917e-06 * tc[2] - 2.96148515e-09 * tc[3] +
                  8.01187908e-13 * tc[4] + 3.83636392e+03 * invT;
  } else {
    // species 21: HCO
    species[21] = +1.60049318e+00 + 2.64639129e-03 * tc[1] -
                  8.97280703e-07 * tc[2] + 1.80339450e-10 * tc[3] -
                  1.48704282e-14 * tc[4] + 4.05725330e+03 * invT;
  }

  // species with midpoint at T=3654.18 kelvin
  if (T < 3654.18) {
    // species 22: CHO+
    species[22] = +1.87627680e+00 + 2.47430100e-03 * tc[1] -
                  7.69510200e-07 * tc[2] + 1.27804380e-10 * tc[3] -
                  8.70412620e-15 * tc[4] + 9.91300220e+04 * invT;
  } else {
    // species 22: CHO+
    species[22] = +5.91572650e+00 + 8.31555350e-05 * tc[1] -
                  6.15878100e-09 * tc[2] + 2.27141072e-13 * tc[3] -
                  3.29088820e-18 * tc[4] + 9.61293330e+04 * invT;
  }

  // species with midpoint at T=5430.9 kelvin
  if (T < 5430.9) {
    // species 23: H3O+
    species[23] = +1.68852320e+00 + 2.72268270e-03 * tc[1] -
                  5.68518867e-07 * tc[2] + 6.20832325e-11 * tc[3] -
                  2.75388800e-15 * tc[4] + 7.08902680e+04 * invT;
  } else {
    // species 23: H3O+
    species[23] = +7.95088030e+00 + 1.26408340e-04 * tc[1] -
                  8.19813233e-09 * tc[2] + 2.71957950e-13 * tc[3] -
                  3.61816740e-18 * tc[4] + 6.44807490e+04 * invT;
  }

  // species with midpoint at T=2008.71 kelvin
  if (T < 2008.71) {
    // species 24: O2-
    species[24] = +2.10217180e+00 + 1.39904375e-03 * tc[1] -
                  7.55037533e-07 * tc[2] + 2.17291292e-10 * tc[3] -
                  2.54437680e-14 * tc[4] - 6.80747930e+03 * invT;
  } else {
    // species 24: O2-
    species[24] = +3.25928670e+00 + 1.12340360e-04 * tc[1] -
                  1.71326517e-08 * tc[2] + 1.83864945e-12 * tc[3] -
                  7.71173040e-17 * tc[4] - 7.24262520e+03 * invT;
  }

  // species with midpoint at T=27300 kelvin
  if (T < 27300) {
    // species 25: E
    species[25] = +1.50000000e+00 - 4.07513040e-18 * tc[1] +
                  4.34328900e-22 * tc[2] - 1.81620833e-26 * tc[3] +
                  2.61324940e-31 * tc[4] - 7.45375000e+02 * invT;
  } else {
    // species 25: E
    species[25] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
  }
}

// compute the h/(RT) at the given temperature (Eq 20)
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEnthalpy(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1 / T;

  // species with midpoint at T=1050 kelvin
  if (T < 1050) {
    // species 0: N2
    species[0] = +3.73100682e+00 - 9.15798650e-04 * tc[1] +
                 1.44108221e-06 * tc[2] - 7.60945377e-10 * tc[3] +
                 1.49214312e-13 * tc[4] - 1.06287426e+03 * invT;
  } else {
    // species 0: N2
    species[0] = +2.81166073e+00 + 8.35336765e-04 * tc[1] -
                 2.26665809e-07 * tc[2] + 3.32203448e-11 * tc[3] -
                 2.05534884e-15 * tc[4] - 8.69811579e+02 * invT;
  }

  // species with midpoint at T=700 kelvin
  if (T < 700) {
    // species 1: H2
    species[1] = +2.64204438e+00 + 2.74764637e-03 * tc[1] -
                 4.23878780e-06 * tc[2] + 3.21872932e-09 * tc[3] -
                 9.40055498e-13 * tc[4] - 9.43236614e+02 * invT;
    // species 3: O2
    species[3] = +3.74403921e+00 - 1.39870073e-03 * tc[1] +
                 3.26707519e-06 * tc[2] - 2.58149107e-09 * tc[3] +
                 7.59862494e-13 * tc[4] - 1.06069827e+03 * invT;
    // species 8: HO2
    species[8] = +3.61994299e+00 + 5.29028525e-04 * tc[1] +
                 1.68892980e-06 * tc[2] - 1.58450190e-09 * tc[3] +
                 4.83194562e-13 * tc[4] + 3.15898234e+02 * invT;
    // species 11: CH4
    species[11] = +5.23967335e+00 - 7.34175615e-03 * tc[1] +
                  1.76577570e-05 * tc[2] - 1.35417205e-08 * tc[3] +
                  3.92637132e-12 * tc[4] - 1.02526308e+04 * invT;
    // species 15: C
    species[15] = +2.54495192e+00 - 1.23862640e-04 * tc[1] +
                  1.82672759e-07 * tc[2] - 1.37137813e-10 * tc[3] +
                  4.08234662e-14 * tc[4] + 8.54434105e+04 * invT;
    // species 20: CH2O
    species[20] = +4.32621296e+00 - 3.50575927e-03 * tc[1] +
                  1.05058987e-05 * tc[2] - 8.41196595e-09 * tc[3] +
                  2.46908046e-12 * tc[4] - 1.43270169e+04 * invT;
  } else {
    // species 1: H2
    species[1] = +3.78199881e+00 - 5.09366295e-04 * tc[1] +
                 4.14087443e-07 * tc[2] - 1.04752974e-10 * tc[3] +
                 9.51087586e-15 * tc[4] - 1.10283023e+03 * invT;
    // species 3: O2
    species[3] = +2.82012408e+00 + 1.24105679e-03 * tc[1] -
                 5.04006980e-07 * tc[2] + 1.12139050e-10 * tc[3] -
                 9.74611336e-15 * tc[4] - 9.31350148e+02 * invT;
    // species 8: HO2
    species[8] = +3.02391889e+00 + 2.23195453e-03 * tc[1] -
                 7.43821640e-07 * tc[2] + 1.53177700e-10 * tc[3] -
                 1.32853247e-14 * tc[4] + 3.99341609e+02 * invT;
    // species 11: CH4
    species[11] = +5.05346405e-01 + 6.18489225e-03 * tc[1] -
                  1.66602641e-06 * tc[2] + 2.60981913e-10 * tc[3] -
                  1.72579483e-14 * tc[4] - 9.58982501e+03 * invT;
    // species 15: C
    species[15] = +2.49472531e+00 + 1.96419738e-05 * tc[1] -
                  2.23338327e-08 * tc[2] + 9.29546735e-12 * tc[3] -
                  1.01461377e-15 * tc[4] + 8.54504422e+04 * invT;
    // species 20: CH2O
    species[20] = +1.33335652e+00 + 5.04525915e-03 * tc[1] -
                  1.70984187e-06 * tc[2] + 3.13563018e-10 * tc[3] -
                  2.39278218e-14 * tc[4] - 1.39080170e+04 * invT;
  }

  // species with midpoint at T=860 kelvin
  if (T < 860) {
    // species 2: H
    species[2] = +2.49950544e+00 + 1.49582023e-06 * tc[1] -
                 1.97586586e-09 * tc[2] + 1.21952541e-12 * tc[3] -
                 2.91078640e-16 * tc[4] + 2.54737866e+04 * invT;
  } else {
    // species 2: H
    species[2] = +2.50031493e+00 - 3.86703414e-07 * tc[1] +
                 2.13115115e-10 * tc[2] - 5.31379477e-14 * tc[3] +
                 4.88958382e-18 * tc[4] + 2.54736474e+04 * invT;
  }

  // species with midpoint at T=720 kelvin
  if (T < 720) {
    // species 4: O
    species[4] = +3.14799201e+00 - 1.55587032e-03 * tc[1] +
                 2.06045966e-06 * tc[2] - 1.40952200e-09 * tc[3] +
                 3.89732032e-13 * tc[4] + 2.91309118e+04 * invT;
  } else {
    // species 4: O
    species[4] = +2.62549143e+00 - 1.04479822e-04 * tc[1] +
                 4.46395153e-08 * tc[2] - 9.64689740e-12 * tc[3] +
                 8.77837378e-16 * tc[4] + 2.92061519e+04 * invT;
  }

  // species with midpoint at T=1420 kelvin
  if (T < 1420) {
    // species 5: H2O
    species[5] = +4.06061172e+00 - 4.32903595e-04 * tc[1] +
                 1.08136509e-06 * tc[2] - 4.50607697e-10 * tc[3] +
                 6.64966586e-14 * tc[4] - 3.02831314e+04 * invT;
  } else {
    // species 5: H2O
    species[5] = +2.66777075e+00 + 1.52884425e-03 * tc[1] -
                 3.00147470e-07 * tc[2] + 3.58403880e-11 * tc[3] -
                 2.01715634e-15 * tc[4] - 2.98875645e+04 * invT;
  }

  // species with midpoint at T=1700 kelvin
  if (T < 1700) {
    // species 6: OH
    species[6] = +3.91354631e+00 - 8.31379630e-04 * tc[1] +
                 7.69733430e-07 * tc[2] - 2.55898770e-10 * tc[3] +
                 3.17659258e-14 * tc[4] + 3.40005047e+03 * invT;
  } else {
    // species 6: OH
    species[6] = +2.49867369e+00 + 8.33176395e-04 * tc[1] -
                 2.09417172e-07 * tc[2] + 3.20867015e-11 * tc[3] -
                 2.11471788e-15 * tc[4] + 3.88110716e+03 * invT;
  }

  // species with midpoint at T=1800 kelvin
  if (T < 1800) {
    // species 7: H2O2
    species[7] = +3.34774224e+00 + 3.52502719e-03 * tc[1] -
                 1.28174002e-06 * tc[2] + 2.91801652e-10 * tc[3] -
                 2.95236210e-14 * tc[4] - 1.75784785e+04 * invT;
    // species 13: CH2
    species[13] = +3.76489460e+00 + 7.19195955e-04 * tc[1] +
                  1.58527692e-07 * tc[2] - 1.07947148e-10 * tc[3] +
                  1.51658575e-14 * tc[4] + 4.58645699e+04 * invT;
    // species 17: CH3OH
    species[17] = +8.47330479e-01 + 8.15434520e-03 * tc[1] -
                  2.82781654e-06 * tc[2] + 5.73260852e-10 * tc[3] -
                  5.19904026e-14 * tc[4] - 2.50962544e+04 * invT;
  } else {
    // species 7: H2O2
    species[7] = +4.76869639e+00 + 1.94618924e-03 * tc[1] -
                 4.04607830e-07 * tc[2] + 4.81538212e-11 * tc[3] -
                 2.45163980e-15 * tc[4] - 1.80900220e+04 * invT;
    // species 13: CH2
    species[13] = +2.81272972e+00 + 1.77715694e-03 * tc[1] -
                  4.29228410e-07 * tc[2] + 5.53184360e-11 * tc[3] -
                  2.97476294e-15 * tc[4] + 4.62073492e+04 * invT;
    // species 17: CH3OH
    species[17] = +2.71701530e+00 + 6.07691765e-03 * tc[1] -
                  1.67369010e-06 * tc[2] + 2.52670175e-10 * tc[3] -
                  1.63692165e-14 * tc[4] - 2.57693409e+04 * invT;
  }

  // species with midpoint at T=960 kelvin
  if (T < 960) {
    // species 9: CO
    species[9] = +3.75723891e+00 - 1.07232621e-03 * tc[1] +
                 1.80693002e-06 * tc[2] - 1.04256491e-09 * tc[3] +
                 2.23802254e-13 * tc[4] - 1.43575530e+04 * invT;
  } else {
    // species 9: CO
    species[9] = +2.79255381e+00 + 9.37434430e-04 * tc[1] -
                 2.86570642e-07 * tc[2] + 4.78000175e-11 * tc[3] -
                 3.35710572e-15 * tc[4] - 1.41723335e+04 * invT;
  }

  // species with midpoint at T=1450 kelvin
  if (T < 1450) {
    // species 10: CO2
    species[10] = +2.31684347e+00 + 4.61377518e-03 * tc[1] -
                  2.58551364e-06 * tc[2] + 8.20563400e-10 * tc[3] -
                  1.09744496e-13 * tc[4] - 4.83626067e+04 * invT;
  } else {
    // species 10: CO2
    species[10] = +4.70876468e+00 + 1.31457352e-03 * tc[1] -
                  3.10202154e-07 * tc[2] + 3.59732300e-11 * tc[3] -
                  1.52516283e-15 * tc[4] - 4.90562639e+04 * invT;
  }

  // species with midpoint at T=1060 kelvin
  if (T < 1060) {
    // species 12: CH3
    species[12] = +3.47829310e+00 + 1.77382386e-03 * tc[1] +
                  4.91361467e-07 * tc[2] - 4.85939888e-10 * tc[3] +
                  1.04384246e-13 * tc[4] + 1.64399516e+04 * invT;
  } else {
    // species 12: CH3
    species[12] = +2.78805104e+00 + 3.07616739e-03 * tc[1] -
                  7.37264497e-07 * tc[2] + 9.36006620e-11 * tc[3] -
                  4.96302698e-15 * tc[4] + 1.65862829e+04 * invT;
  }

  // species with midpoint at T=970 kelvin
  if (T < 970) {
    // species 14: CH2(S)
    species[14] = +4.18185434e+00 - 1.10567157e-03 * tc[1] +
                  2.57175847e-06 * tc[2] - 1.48987595e-09 * tc[3] +
                  3.16629256e-13 * tc[4] + 5.03669407e+04 * invT;
  } else {
    // species 14: CH2(S)
    species[14] = +2.75934299e+00 + 1.82734153e-03 * tc[1] -
                  4.51966380e-07 * tc[2] + 6.87451027e-11 * tc[3] -
                  4.73590944e-15 * tc[4] + 5.06429079e+04 * invT;
  }

  // species with midpoint at T=1590 kelvin
  if (T < 1590) {
    // species 16: CH
    species[16] = +3.77264332e+00 - 7.92736750e-04 * tc[1] +
                  9.45040793e-07 * tc[2] - 3.40365145e-10 * tc[3] +
                  4.47990664e-14 * tc[4] + 7.06312492e+04 * invT;
  } else {
    // species 16: CH
    species[16] = +2.27990128e+00 + 1.08492619e-03 * tc[1] -
                  2.35879295e-07 * tc[2] + 3.09933735e-11 * tc[3] -
                  1.91269703e-15 * tc[4] + 7.11059412e+04 * invT;
  }

  // species with midpoint at T=1740 kelvin
  if (T < 1740) {
    // species 18: CH3O
    species[18] = +8.89660985e-01 + 8.50598835e-03 * tc[1] -
                  3.79357837e-06 * tc[2] + 9.70702320e-10 * tc[3] -
                  1.06568296e-13 * tc[4] + 1.60316121e+03 * invT;
  } else {
    // species 18: CH3O
    species[18] = +5.72238062e+00 + 2.95113818e-03 * tc[1] -
                  6.01135733e-07 * tc[2] + 5.33337522e-11 * tc[3] -
                  1.12363282e-15 * tc[4] - 7.86252225e+01 * invT;
  }

  // species with midpoint at T=1360 kelvin
  if (T < 1360) {
    // species 19: CH2OH
    species[19] = +2.34821579e+00 + 6.98000840e-03 * tc[1] -
                  3.62107353e-06 * tc[2] + 1.15624604e-09 * tc[3] -
                  1.61699832e-13 * tc[4] - 3.30222106e+03 * invT;
  } else {
    // species 19: CH2OH
    species[19] = +5.04534950e+00 + 3.01363529e-03 * tc[1] -
                  7.04622737e-07 * tc[2] + 8.40214762e-11 * tc[3] -
                  4.01974964e-15 * tc[4] - 4.03584143e+03 * invT;
  }

  // species with midpoint at T=770 kelvin
  if (T < 770) {
    // species 21: HCO
    species[21] = +4.03483979e+00 - 1.07918432e-03 * tc[1] +
                  3.94112917e-06 * tc[2] - 2.96148515e-09 * tc[3] +
                  8.01187908e-13 * tc[4] + 3.83636392e+03 * invT;
  } else {
    // species 21: HCO
    species[21] = +2.60049318e+00 + 2.64639129e-03 * tc[1] -
                  8.97280703e-07 * tc[2] + 1.80339450e-10 * tc[3] -
                  1.48704282e-14 * tc[4] + 4.05725330e+03 * invT;
  }

  // species with midpoint at T=3654.18 kelvin
  if (T < 3654.18) {
    // species 22: CHO+
    species[22] = +2.87627680e+00 + 2.47430100e-03 * tc[1] -
                  7.69510200e-07 * tc[2] + 1.27804380e-10 * tc[3] -
                  8.70412620e-15 * tc[4] + 9.91300220e+04 * invT;
  } else {
    // species 22: CHO+
    species[22] = +6.91572650e+00 + 8.31555350e-05 * tc[1] -
                  6.15878100e-09 * tc[2] + 2.27141072e-13 * tc[3] -
                  3.29088820e-18 * tc[4] + 9.61293330e+04 * invT;
  }

  // species with midpoint at T=5430.9 kelvin
  if (T < 5430.9) {
    // species 23: H3O+
    species[23] = +2.68852320e+00 + 2.72268270e-03 * tc[1] -
                  5.68518867e-07 * tc[2] + 6.20832325e-11 * tc[3] -
                  2.75388800e-15 * tc[4] + 7.08902680e+04 * invT;
  } else {
    // species 23: H3O+
    species[23] = +8.95088030e+00 + 1.26408340e-04 * tc[1] -
                  8.19813233e-09 * tc[2] + 2.71957950e-13 * tc[3] -
                  3.61816740e-18 * tc[4] + 6.44807490e+04 * invT;
  }

  // species with midpoint at T=2008.71 kelvin
  if (T < 2008.71) {
    // species 24: O2-
    species[24] = +3.10217180e+00 + 1.39904375e-03 * tc[1] -
                  7.55037533e-07 * tc[2] + 2.17291292e-10 * tc[3] -
                  2.54437680e-14 * tc[4] - 6.80747930e+03 * invT;
  } else {
    // species 24: O2-
    species[24] = +4.25928670e+00 + 1.12340360e-04 * tc[1] -
                  1.71326517e-08 * tc[2] + 1.83864945e-12 * tc[3] -
                  7.71173040e-17 * tc[4] - 7.24262520e+03 * invT;
  }

  // species with midpoint at T=27300 kelvin
  if (T < 27300) {
    // species 25: E
    species[25] = +2.50000000e+00 - 4.07513040e-18 * tc[1] +
                  4.34328900e-22 * tc[2] - 1.81620833e-26 * tc[3] +
                  2.61324940e-31 * tc[4] - 7.45375000e+02 * invT;
  } else {
    // species 25: E
    species[25] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
  }
}

// compute the S/R at the given temperature (Eq 21)
// tc contains precomputed powers of T, tc[0] = log(T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEntropy(amrex::Real* species, const amrex::Real* tc)
{

  // temperature
  const amrex::Real T = tc[1];

  // species with midpoint at T=1050 kelvin
  if (T < 1050) {
    // species 0: N2
    species[0] = +3.73100682e+00 * tc[0] - 1.83159730e-03 * tc[1] +
                 2.16162331e-06 * tc[2] - 1.01459384e-09 * tc[3] +
                 1.86517890e-13 * tc[4] + 2.16821198e+00;
  } else {
    // species 0: N2
    species[0] = +2.81166073e+00 * tc[0] + 1.67067353e-03 * tc[1] -
                 3.39998714e-07 * tc[2] + 4.42937930e-11 * tc[3] -
                 2.56918605e-15 * tc[4] + 6.64838050e+00;
  }

  // species with midpoint at T=700 kelvin
  if (T < 700) {
    // species 1: H2
    species[1] = +2.64204438e+00 * tc[0] + 5.49529274e-03 * tc[1] -
                 6.35818170e-06 * tc[2] + 4.29163910e-09 * tc[3] -
                 1.17506937e-12 * tc[4] - 5.12231102e-01;
    // species 3: O2
    species[3] = +3.74403921e+00 * tc[0] - 2.79740147e-03 * tc[1] +
                 4.90061279e-06 * tc[2] - 3.44198810e-09 * tc[3] +
                 9.49828117e-13 * tc[4] + 3.82132646e+00;
    // species 8: HO2
    species[8] = +3.61994299e+00 * tc[0] + 1.05805705e-03 * tc[1] +
                 2.53339470e-06 * tc[2] - 2.11266921e-09 * tc[3] +
                 6.03993202e-13 * tc[4] + 6.44411482e+00;
    // species 11: CH4
    species[11] = +5.23967335e+00 * tc[0] - 1.46835123e-02 * tc[1] +
                  2.64866356e-05 * tc[2] - 1.80556274e-08 * tc[3] +
                  4.90796415e-12 * tc[4] - 4.97649748e+00;
    // species 15: C
    species[15] = +2.54495192e+00 * tc[0] - 2.47725281e-04 * tc[1] +
                  2.74009139e-07 * tc[2] - 1.82850417e-10 * tc[3] +
                  5.10293327e-14 * tc[4] + 4.56874273e+00;
    // species 20: CH2O
    species[20] = +4.32621296e+00 * tc[0] - 7.01151853e-03 * tc[1] +
                  1.57588480e-05 * tc[2] - 1.12159546e-08 * tc[3] +
                  3.08635057e-12 * tc[4] + 2.62028902e+00;
  } else {
    // species 1: H2
    species[1] = +3.78199881e+00 * tc[0] - 1.01873259e-03 * tc[1] +
                 6.21131165e-07 * tc[2] - 1.39670633e-10 * tc[3] +
                 1.18885948e-14 * tc[4] - 5.60525910e+00;
    // species 3: O2
    species[3] = +2.82012408e+00 * tc[0] + 2.48211357e-03 * tc[1] -
                 7.56010470e-07 * tc[2] + 1.49518734e-10 * tc[3] -
                 1.21826417e-14 * tc[4] + 7.94914552e+00;
    // species 8: HO2
    species[8] = +3.02391889e+00 * tc[0] + 4.46390907e-03 * tc[1] -
                 1.11573246e-06 * tc[2] + 2.04236933e-10 * tc[3] -
                 1.66066559e-14 * tc[4] + 9.10699973e+00;
    // species 11: CH4
    species[11] = +5.05346405e-01 * tc[0] + 1.23697845e-02 * tc[1] -
                  2.49903961e-06 * tc[2] + 3.47975883e-10 * tc[3] -
                  2.15724354e-14 * tc[4] + 1.61752775e+01;
    // species 15: C
    species[15] = +2.49472531e+00 * tc[0] + 3.92839476e-05 * tc[1] -
                  3.35007490e-08 * tc[2] + 1.23939565e-11 * tc[3] -
                  1.26826721e-15 * tc[4] + 4.79314254e+00;
    // species 20: CH2O
    species[20] = +1.33335652e+00 * tc[0] + 1.00905183e-02 * tc[1] -
                  2.56476281e-06 * tc[2] + 4.18084023e-10 * tc[3] -
                  2.99097773e-14 * tc[4] + 1.59916144e+01;
  }

  // species with midpoint at T=860 kelvin
  if (T < 860) {
    // species 2: H
    species[2] = +2.49950544e+00 * tc[0] + 2.99164046e-06 * tc[1] -
                 2.96379880e-09 * tc[2] + 1.62603388e-12 * tc[3] -
                 3.63848300e-16 * tc[4] - 4.44574018e-01;
  } else {
    // species 2: H
    species[2] = +2.50031493e+00 * tc[0] - 7.73406828e-07 * tc[1] +
                 3.19672673e-10 * tc[2] - 7.08505970e-14 * tc[3] +
                 6.11197978e-18 * tc[4] - 4.48357228e-01;
  }

  // species with midpoint at T=720 kelvin
  if (T < 720) {
    // species 4: O
    species[4] = +3.14799201e+00 * tc[0] - 3.11174065e-03 * tc[1] +
                 3.09068949e-06 * tc[2] - 1.87936266e-09 * tc[3] +
                 4.87165040e-13 * tc[4] + 2.13446549e+00;
  } else {
    // species 4: O
    species[4] = +2.62549143e+00 * tc[0] - 2.08959644e-04 * tc[1] +
                 6.69592730e-08 * tc[2] - 1.28625299e-11 * tc[3] +
                 1.09729672e-15 * tc[4] + 4.48358519e+00;
  }

  // species with midpoint at T=1420 kelvin
  if (T < 1420) {
    // species 5: H2O
    species[5] = +4.06061172e+00 * tc[0] - 8.65807189e-04 * tc[1] +
                 1.62204764e-06 * tc[2] - 6.00810263e-10 * tc[3] +
                 8.31208233e-14 * tc[4] - 2.96150481e-01;
  } else {
    // species 5: H2O
    species[5] = +2.66777075e+00 * tc[0] + 3.05768849e-03 * tc[1] -
                 4.50221205e-07 * tc[2] + 4.77871840e-11 * tc[3] -
                 2.52144543e-15 * tc[4] + 6.91191131e+00;
  }

  // species with midpoint at T=1700 kelvin
  if (T < 1700) {
    // species 6: OH
    species[6] = +3.91354631e+00 * tc[0] - 1.66275926e-03 * tc[1] +
                 1.15460014e-06 * tc[2] - 3.41198360e-10 * tc[3] +
                 3.97074073e-14 * tc[4] + 2.05474719e-01;
  } else {
    // species 6: OH
    species[6] = +2.49867369e+00 * tc[0] + 1.66635279e-03 * tc[1] -
                 3.14125758e-07 * tc[2] + 4.27822687e-11 * tc[3] -
                 2.64339735e-15 * tc[4] + 7.78218862e+00;
  }

  // species with midpoint at T=1800 kelvin
  if (T < 1800) {
    // species 7: H2O2
    species[7] = +3.34774224e+00 * tc[0] + 7.05005437e-03 * tc[1] -
                 1.92261003e-06 * tc[2] + 3.89068870e-10 * tc[3] -
                 3.69045262e-14 * tc[4] + 7.17868851e+00;
    // species 13: CH2
    species[13] = +3.76489460e+00 * tc[0] + 1.43839191e-03 * tc[1] +
                  2.37791539e-07 * tc[2] - 1.43929530e-10 * tc[3] +
                  1.89573218e-14 * tc[4] + 1.48953153e+00;
    // species 17: CH3OH
    species[17] = +8.47330479e-01 * tc[0] + 1.63086904e-02 * tc[1] -
                  4.24172480e-06 * tc[2] + 7.64347803e-10 * tc[3] -
                  6.49880033e-14 * tc[4] + 1.95933297e+01;
  } else {
    // species 7: H2O2
    species[7] = +4.76869639e+00 * tc[0] + 3.89237848e-03 * tc[1] -
                 6.06911745e-07 * tc[2] + 6.42050950e-11 * tc[3] -
                 3.06454975e-15 * tc[4] - 5.11811777e-01;
    // species 13: CH2
    species[13] = +2.81272972e+00 * tc[0] + 3.55431388e-03 * tc[1] -
                  6.43842615e-07 * tc[2] + 7.37579147e-11 * tc[3] -
                  3.71845367e-15 * tc[4] + 6.64284652e+00;
    // species 17: CH3OH
    species[17] = +2.71701530e+00 * tc[0] + 1.21538353e-02 * tc[1] -
                  2.51053516e-06 * tc[2] + 3.36893567e-10 * tc[3] -
                  2.04615206e-14 * tc[4] + 9.47420540e+00;
  }

  // species with midpoint at T=960 kelvin
  if (T < 960) {
    // species 9: CO
    species[9] = +3.75723891e+00 * tc[0] - 2.14465241e-03 * tc[1] +
                 2.71039502e-06 * tc[2] - 1.39008654e-09 * tc[3] +
                 2.79752817e-13 * tc[4] + 2.79976799e+00;
  } else {
    // species 9: CO
    species[9] = +2.79255381e+00 * tc[0] + 1.87486886e-03 * tc[1] -
                 4.29855963e-07 * tc[2] + 6.37333567e-11 * tc[3] -
                 4.19638215e-15 * tc[4] + 7.41443560e+00;
  }

  // species with midpoint at T=1450 kelvin
  if (T < 1450) {
    // species 10: CO2
    species[10] = +2.31684347e+00 * tc[0] + 9.22755036e-03 * tc[1] -
                  3.87827046e-06 * tc[2] + 1.09408453e-09 * tc[3] -
                  1.37180621e-13 * tc[4] + 1.00786234e+01;
  } else {
    // species 10: CO2
    species[10] = +4.70876468e+00 * tc[0] + 2.62914704e-03 * tc[1] -
                  4.65303231e-07 * tc[2] + 4.79643067e-11 * tc[3] -
                  1.90645353e-15 * tc[4] - 2.34976452e+00;
  }

  // species with midpoint at T=1060 kelvin
  if (T < 1060) {
    // species 12: CH3
    species[12] = +3.47829310e+00 * tc[0] + 3.54764773e-03 * tc[1] +
                  7.37042200e-07 * tc[2] - 6.47919850e-10 * tc[3] +
                  1.30480308e-13 * tc[4] + 2.40875956e+00;
  } else {
    // species 12: CH3
    species[12] = +2.78805104e+00 * tc[0] + 6.15233477e-03 * tc[1] -
                  1.10589675e-06 * tc[2] + 1.24800883e-10 * tc[3] -
                  6.20378373e-15 * tc[4] + 5.77899818e+00;
  }

  // species with midpoint at T=970 kelvin
  if (T < 970) {
    // species 14: CH2(S)
    species[14] = +4.18185434e+00 * tc[0] - 2.21134314e-03 * tc[1] +
                  3.85763771e-06 * tc[2] - 1.98650127e-09 * tc[3] +
                  3.95786570e-13 * tc[4] - 7.03002602e-01;
  } else {
    // species 14: CH2(S)
    species[14] = +2.75934299e+00 * tc[0] + 3.65468307e-03 * tc[1] -
                  6.77949570e-07 * tc[2] + 9.16601370e-11 * tc[3] -
                  5.91988680e-15 * tc[4] + 6.11646383e+00;
  }

  // species with midpoint at T=1590 kelvin
  if (T < 1590) {
    // species 16: CH
    species[16] = +3.77264332e+00 * tc[0] - 1.58547350e-03 * tc[1] +
                  1.41756119e-06 * tc[2] - 4.53820193e-10 * tc[3] +
                  5.59988330e-14 * tc[4] + 8.79407904e-01;
  } else {
    // species 16: CH
    species[16] = +2.27990128e+00 * tc[0] + 2.16985238e-03 * tc[1] -
                  3.53818943e-07 * tc[2] + 4.13244980e-11 * tc[3] -
                  2.39087129e-15 * tc[4] + 8.77326061e+00;
  }

  // species with midpoint at T=1740 kelvin
  if (T < 1740) {
    // species 18: CH3O
    species[18] = +8.89660985e-01 * tc[0] + 1.70119767e-02 * tc[1] -
                  5.69036755e-06 * tc[2] + 1.29426976e-09 * tc[3] -
                  1.33210370e-13 * tc[4] + 1.85001134e+01;
  } else {
    // species 18: CH3O
    species[18] = +5.72238062e+00 * tc[0] + 5.90227637e-03 * tc[1] -
                  9.01703600e-07 * tc[2] + 7.11116697e-11 * tc[3] -
                  1.40454102e-15 * tc[4] - 7.49173677e+00;
  }

  // species with midpoint at T=1360 kelvin
  if (T < 1360) {
    // species 19: CH2OH
    species[19] = +2.34821579e+00 * tc[0] + 1.39600168e-02 * tc[1] -
                  5.43161030e-06 * tc[2] + 1.54166139e-09 * tc[3] -
                  2.02124790e-13 * tc[4] + 1.22661977e+01;
  } else {
    // species 19: CH2OH
    species[19] = +5.04534950e+00 * tc[0] + 6.02727059e-03 * tc[1] -
                  1.05693410e-06 * tc[2] + 1.12028635e-10 * tc[3] -
                  5.02468705e-15 * tc[4] - 1.57524073e+00;
  }

  // species with midpoint at T=770 kelvin
  if (T < 770) {
    // species 21: HCO
    species[21] = +4.03483979e+00 * tc[0] - 2.15836864e-03 * tc[1] +
                  5.91169375e-06 * tc[2] - 3.94864687e-09 * tc[3] +
                  1.00148488e-12 * tc[4] + 4.20008770e+00;
  } else {
    // species 21: HCO
    species[21] = +2.60049318e+00 * tc[0] + 5.29278258e-03 * tc[1] -
                  1.34592105e-06 * tc[2] + 2.40452599e-10 * tc[3] -
                  1.85880352e-14 * tc[4] + 1.07450933e+01;
  }

  // species with midpoint at T=3654.18 kelvin
  if (T < 3654.18) {
    // species 22: CHO+
    species[22] = +2.87627680e+00 * tc[0] + 4.94860200e-03 * tc[1] -
                  1.15426530e-06 * tc[2] + 1.70405840e-10 * tc[3] -
                  1.08801578e-14 * tc[4] + 6.50202660e+00;
  } else {
    // species 22: CHO+
    species[22] = +6.91572650e+00 * tc[0] + 1.66311070e-04 * tc[1] -
                  9.23817150e-09 * tc[2] + 3.02854763e-13 * tc[3] -
                  4.11361025e-18 * tc[4] - 1.80895100e+01;
  }

  // species with midpoint at T=5430.9 kelvin
  if (T < 5430.9) {
    // species 23: H3O+
    species[23] = +2.68852320e+00 * tc[0] + 5.44536540e-03 * tc[1] -
                  8.52778300e-07 * tc[2] + 8.27776433e-11 * tc[3] -
                  3.44236000e-15 * tc[4] + 6.36082170e+00;
  } else {
    // species 23: H3O+
    species[23] = +8.95088030e+00 * tc[0] + 2.52816680e-04 * tc[1] -
                  1.22971985e-08 * tc[2] + 3.62610600e-13 * tc[3] -
                  4.52270925e-18 * tc[4] - 3.38733190e+01;
  }

  // species with midpoint at T=2008.71 kelvin
  if (T < 2008.71) {
    // species 24: O2-
    species[24] = +3.10217180e+00 * tc[0] + 2.79808750e-03 * tc[1] -
                  1.13255630e-06 * tc[2] + 2.89721723e-10 * tc[3] -
                  3.18047100e-14 * tc[4] + 6.76090200e+00;
  } else {
    // species 24: O2-
    species[24] = +4.25928670e+00 * tc[0] + 2.24680720e-04 * tc[1] -
                  2.56989775e-08 * tc[2] + 2.45153260e-12 * tc[3] -
                  9.63966300e-17 * tc[4] + 4.75996970e-01;
  }

  // species with midpoint at T=27300 kelvin
  if (T < 27300) {
    // species 25: E
    species[25] = +2.50000000e+00 * tc[0] - 8.15026080e-18 * tc[1] +
                  6.51493350e-22 * tc[2] - 2.42161110e-26 * tc[3] +
                  3.26656175e-31 * tc[4] - 1.17339750e+01;
  } else {
    // species 25: E
    species[25] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4] - 1.17339750e+01;
  }
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
  amrex::Real cpor[26];                               //  temporary storage
  cp_R(cpor, tc);

  // perform dot product
  for (id = 0; id < 26; ++id) {
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
  amrex::Real cpor[26], tresult[26];                  //  temporary storage
  amrex::Real imw[26];

  get_imw(imw);

  cp_R(cpor, tc);

  for (int i = 0; i < 26; i++) {
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
  amrex::Real cvor[26];                               //  temporary storage
  cv_R(cvor, tc);

  // perform dot product
  for (id = 0; id < 26; ++id) {
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
  amrex::Real cvor[26];                               //  temporary storage
  amrex::Real imw[26];

  get_imw(imw);

  cv_R(cvor, tc);

  // multiply by y/molecularweight
  result += cvor[0] * y[0] * imw[0];    // N2
  result += cvor[1] * y[1] * imw[1];    // H2
  result += cvor[2] * y[2] * imw[2];    // H
  result += cvor[3] * y[3] * imw[3];    // O2
  result += cvor[4] * y[4] * imw[4];    // O
  result += cvor[5] * y[5] * imw[5];    // H2O
  result += cvor[6] * y[6] * imw[6];    // OH
  result += cvor[7] * y[7] * imw[7];    // H2O2
  result += cvor[8] * y[8] * imw[8];    // HO2
  result += cvor[9] * y[9] * imw[9];    // CO
  result += cvor[10] * y[10] * imw[10]; // CO2
  result += cvor[11] * y[11] * imw[11]; // CH4
  result += cvor[12] * y[12] * imw[12]; // CH3
  result += cvor[13] * y[13] * imw[13]; // CH2
  result += cvor[14] * y[14] * imw[14]; // CH2(S)
  result += cvor[15] * y[15] * imw[15]; // C
  result += cvor[16] * y[16] * imw[16]; // CH
  result += cvor[17] * y[17] * imw[17]; // CH3OH
  result += cvor[18] * y[18] * imw[18]; // CH3O
  result += cvor[19] * y[19] * imw[19]; // CH2OH
  result += cvor[20] * y[20] * imw[20]; // CH2O
  result += cvor[21] * y[21] * imw[21]; // HCO
  result += cvor[22] * y[22] * imw[22]; // CHO+
  result += cvor[23] * y[23] * imw[23]; // H3O+
  result += cvor[24] * y[24] * imw[24]; // O2-
  result += cvor[25] * y[25] * imw[25]; // E

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
  amrex::Real hml[26];                                //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesEnthalpy(hml, tc);

  // perform dot product
  for (id = 0; id < 26; ++id) {
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
  amrex::Real hml[26], tmp[26];                       //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[26];

  get_imw(imw);

  speciesEnthalpy(hml, tc);

  int id;
  for (id = 0; id < 26; ++id) {
    tmp[id] = y[id] * hml[id] * imw[id];
  }
  for (id = 0; id < 26; ++id) {
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
  amrex::Real uml[26];                                //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  speciesInternalEnergy(uml, tc);

  // perform dot product
  for (id = 0; id < 26; ++id) {
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
  amrex::Real ums[26];                                //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[26];

  get_imw(imw);

  speciesInternalEnergy(ums, tc);

  // perform dot product + scaling by wt
  result += y[0] * ums[0] * imw[0];    // N2
  result += y[1] * ums[1] * imw[1];    // H2
  result += y[2] * ums[2] * imw[2];    // H
  result += y[3] * ums[3] * imw[3];    // O2
  result += y[4] * ums[4] * imw[4];    // O
  result += y[5] * ums[5] * imw[5];    // H2O
  result += y[6] * ums[6] * imw[6];    // OH
  result += y[7] * ums[7] * imw[7];    // H2O2
  result += y[8] * ums[8] * imw[8];    // HO2
  result += y[9] * ums[9] * imw[9];    // CO
  result += y[10] * ums[10] * imw[10]; // CO2
  result += y[11] * ums[11] * imw[11]; // CH4
  result += y[12] * ums[12] * imw[12]; // CH3
  result += y[13] * ums[13] * imw[13]; // CH2
  result += y[14] * ums[14] * imw[14]; // CH2(S)
  result += y[15] * ums[15] * imw[15]; // C
  result += y[16] * ums[16] * imw[16]; // CH
  result += y[17] * ums[17] * imw[17]; // CH3OH
  result += y[18] * ums[18] * imw[18]; // CH3O
  result += y[19] * ums[19] * imw[19]; // CH2OH
  result += y[20] * ums[20] * imw[20]; // CH2O
  result += y[21] * ums[21] * imw[21]; // HCO
  result += y[22] * ums[22] * imw[22]; // CHO+
  result += y[23] * ums[23] * imw[23]; // H3O+
  result += y[24] * ums[24] * imw[24]; // O2-
  result += y[25] * ums[25] * imw[25]; // E

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
  amrex::Real sor[26]; //  temporary storage
  speciesEntropy(sor, tc);

  // Compute Eq 42
  for (id = 0; id < 26; ++id) {
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
  amrex::Real sor[26]; //  temporary storage
  amrex::Real x[26];   //  need a ytx conversion
  amrex::Real YOW = 0; // See Eq 4, 6 in CK Manual
  amrex::Real imw[26];

  get_imw(imw);

  // Compute inverse of mean molecular wt first
  YOW += y[0] * imw[0];   // N2
  YOW += y[1] * imw[1];   // H2
  YOW += y[2] * imw[2];   // H
  YOW += y[3] * imw[3];   // O2
  YOW += y[4] * imw[4];   // O
  YOW += y[5] * imw[5];   // H2O
  YOW += y[6] * imw[6];   // OH
  YOW += y[7] * imw[7];   // H2O2
  YOW += y[8] * imw[8];   // HO2
  YOW += y[9] * imw[9];   // CO
  YOW += y[10] * imw[10]; // CO2
  YOW += y[11] * imw[11]; // CH4
  YOW += y[12] * imw[12]; // CH3
  YOW += y[13] * imw[13]; // CH2
  YOW += y[14] * imw[14]; // CH2(S)
  YOW += y[15] * imw[15]; // C
  YOW += y[16] * imw[16]; // CH
  YOW += y[17] * imw[17]; // CH3OH
  YOW += y[18] * imw[18]; // CH3O
  YOW += y[19] * imw[19]; // CH2OH
  YOW += y[20] * imw[20]; // CH2O
  YOW += y[21] * imw[21]; // HCO
  YOW += y[22] * imw[22]; // CHO+
  YOW += y[23] * imw[23]; // H3O+
  YOW += y[24] * imw[24]; // O2-
  YOW += y[25] * imw[25]; // E
  // Now compute y to x conversion
  x[0] = y[0] / (28.014000 * YOW);
  x[1] = y[1] / (2.016000 * YOW);
  x[2] = y[2] / (1.008000 * YOW);
  x[3] = y[3] / (31.998000 * YOW);
  x[4] = y[4] / (15.999000 * YOW);
  x[5] = y[5] / (18.015000 * YOW);
  x[6] = y[6] / (17.007000 * YOW);
  x[7] = y[7] / (34.014000 * YOW);
  x[8] = y[8] / (33.006000 * YOW);
  x[9] = y[9] / (28.010000 * YOW);
  x[10] = y[10] / (44.009000 * YOW);
  x[11] = y[11] / (16.043000 * YOW);
  x[12] = y[12] / (15.035000 * YOW);
  x[13] = y[13] / (14.027000 * YOW);
  x[14] = y[14] / (14.027000 * YOW);
  x[15] = y[15] / (12.011000 * YOW);
  x[16] = y[16] / (13.019000 * YOW);
  x[17] = y[17] / (32.042000 * YOW);
  x[18] = y[18] / (31.034000 * YOW);
  x[19] = y[19] / (31.034000 * YOW);
  x[20] = y[20] / (30.026000 * YOW);
  x[21] = y[21] / (29.018000 * YOW);
  x[22] = y[22] / (29.017451 * YOW);
  x[23] = y[23] / (19.022451 * YOW);
  x[24] = y[24] / (31.998549 * YOW);
  x[25] = y[25] / (0.000549 * YOW);
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
  result += x[12] * (sor[12] - log((x[12] + 1e-100)) - logPratio);
  result += x[13] * (sor[13] - log((x[13] + 1e-100)) - logPratio);
  result += x[14] * (sor[14] - log((x[14] + 1e-100)) - logPratio);
  result += x[15] * (sor[15] - log((x[15] + 1e-100)) - logPratio);
  result += x[16] * (sor[16] - log((x[16] + 1e-100)) - logPratio);
  result += x[17] * (sor[17] - log((x[17] + 1e-100)) - logPratio);
  result += x[18] * (sor[18] - log((x[18] + 1e-100)) - logPratio);
  result += x[19] * (sor[19] - log((x[19] + 1e-100)) - logPratio);
  result += x[20] * (sor[20] - log((x[20] + 1e-100)) - logPratio);
  result += x[21] * (sor[21] - log((x[21] + 1e-100)) - logPratio);
  result += x[22] * (sor[22] - log((x[22] + 1e-100)) - logPratio);
  result += x[23] * (sor[23] - log((x[23] + 1e-100)) - logPratio);
  result += x[24] * (sor[24] - log((x[24] + 1e-100)) - logPratio);
  result += x[25] * (sor[25] - log((x[25] + 1e-100)) - logPratio);
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
  XW += x[0] * 28.014000;                       // N2
  XW += x[1] * 2.016000;                        // H2
  XW += x[2] * 1.008000;                        // H
  XW += x[3] * 31.998000;                       // O2
  XW += x[4] * 15.999000;                       // O
  XW += x[5] * 18.015000;                       // H2O
  XW += x[6] * 17.007000;                       // OH
  XW += x[7] * 34.014000;                       // H2O2
  XW += x[8] * 33.006000;                       // HO2
  XW += x[9] * 28.010000;                       // CO
  XW += x[10] * 44.009000;                      // CO2
  XW += x[11] * 16.043000;                      // CH4
  XW += x[12] * 15.035000;                      // CH3
  XW += x[13] * 14.027000;                      // CH2
  XW += x[14] * 14.027000;                      // CH2(S)
  XW += x[15] * 12.011000;                      // C
  XW += x[16] * 13.019000;                      // CH
  XW += x[17] * 32.042000;                      // CH3OH
  XW += x[18] * 31.034000;                      // CH3O
  XW += x[19] * 31.034000;                      // CH2OH
  XW += x[20] * 30.026000;                      // CH2O
  XW += x[21] * 29.018000;                      // HCO
  XW += x[22] * 29.017451;                      // CHO+
  XW += x[23] * 19.022451;                      // H3O+
  XW += x[24] * 31.998549;                      // O2-
  XW += x[25] * 0.000549;                       // E
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
  amrex::Real imw[26];

  get_imw(imw);

  YOW += y[0] * imw[0];                          // N2
  YOW += y[1] * imw[1];                          // H2
  YOW += y[2] * imw[2];                          // H
  YOW += y[3] * imw[3];                          // O2
  YOW += y[4] * imw[4];                          // O
  YOW += y[5] * imw[5];                          // H2O
  YOW += y[6] * imw[6];                          // OH
  YOW += y[7] * imw[7];                          // H2O2
  YOW += y[8] * imw[8];                          // HO2
  YOW += y[9] * imw[9];                          // CO
  YOW += y[10] * imw[10];                        // CO2
  YOW += y[11] * imw[11];                        // CH4
  YOW += y[12] * imw[12];                        // CH3
  YOW += y[13] * imw[13];                        // CH2
  YOW += y[14] * imw[14];                        // CH2(S)
  YOW += y[15] * imw[15];                        // C
  YOW += y[16] * imw[16];                        // CH
  YOW += y[17] * imw[17];                        // CH3OH
  YOW += y[18] * imw[18];                        // CH3O
  YOW += y[19] * imw[19];                        // CH2OH
  YOW += y[20] * imw[20];                        // CH2O
  YOW += y[21] * imw[21];                        // HCO
  YOW += y[22] * imw[22];                        // CHO+
  YOW += y[23] * imw[23];                        // H3O+
  YOW += y[24] * imw[24];                        // O2-
  YOW += y[25] * imw[25];                        // E
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
  W += c[0] * 28.014000;  // N2
  W += c[1] * 2.016000;   // H2
  W += c[2] * 1.008000;   // H
  W += c[3] * 31.998000;  // O2
  W += c[4] * 15.999000;  // O
  W += c[5] * 18.015000;  // H2O
  W += c[6] * 17.007000;  // OH
  W += c[7] * 34.014000;  // H2O2
  W += c[8] * 33.006000;  // HO2
  W += c[9] * 28.010000;  // CO
  W += c[10] * 44.009000; // CO2
  W += c[11] * 16.043000; // CH4
  W += c[12] * 15.035000; // CH3
  W += c[13] * 14.027000; // CH2
  W += c[14] * 14.027000; // CH2(S)
  W += c[15] * 12.011000; // C
  W += c[16] * 13.019000; // CH
  W += c[17] * 32.042000; // CH3OH
  W += c[18] * 31.034000; // CH3O
  W += c[19] * 31.034000; // CH2OH
  W += c[20] * 30.026000; // CH2O
  W += c[21] * 29.018000; // HCO
  W += c[22] * 29.017451; // CHO+
  W += c[23] * 19.022451; // H3O+
  W += c[24] * 31.998549; // O2-
  W += c[25] * 0.000549;  // E

  for (id = 0; id < 26; ++id) {
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
  XW += x[0] * 28.014000;                         // N2
  XW += x[1] * 2.016000;                          // H2
  XW += x[2] * 1.008000;                          // H
  XW += x[3] * 31.998000;                         // O2
  XW += x[4] * 15.999000;                         // O
  XW += x[5] * 18.015000;                         // H2O
  XW += x[6] * 17.007000;                         // OH
  XW += x[7] * 34.014000;                         // H2O2
  XW += x[8] * 33.006000;                         // HO2
  XW += x[9] * 28.010000;                         // CO
  XW += x[10] * 44.009000;                        // CO2
  XW += x[11] * 16.043000;                        // CH4
  XW += x[12] * 15.035000;                        // CH3
  XW += x[13] * 14.027000;                        // CH2
  XW += x[14] * 14.027000;                        // CH2(S)
  XW += x[15] * 12.011000;                        // C
  XW += x[16] * 13.019000;                        // CH
  XW += x[17] * 32.042000;                        // CH3OH
  XW += x[18] * 31.034000;                        // CH3O
  XW += x[19] * 31.034000;                        // CH2OH
  XW += x[20] * 30.026000;                        // CH2O
  XW += x[21] * 29.018000;                        // HCO
  XW += x[22] * 29.017451;                        // CHO+
  XW += x[23] * 19.022451;                        // H3O+
  XW += x[24] * 31.998549;                        // O2-
  XW += x[25] * 0.000549;                         // E
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
  amrex::Real tmp[26];
  amrex::Real imw[26];

  get_imw(imw);

  for (int i = 0; i < 26; i++) {
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
  W += c[0] * 28.014000;  // N2
  W += c[1] * 2.016000;   // H2
  W += c[2] * 1.008000;   // H
  W += c[3] * 31.998000;  // O2
  W += c[4] * 15.999000;  // O
  W += c[5] * 18.015000;  // H2O
  W += c[6] * 17.007000;  // OH
  W += c[7] * 34.014000;  // H2O2
  W += c[8] * 33.006000;  // HO2
  W += c[9] * 28.010000;  // CO
  W += c[10] * 44.009000; // CO2
  W += c[11] * 16.043000; // CH4
  W += c[12] * 15.035000; // CH3
  W += c[13] * 14.027000; // CH2
  W += c[14] * 14.027000; // CH2(S)
  W += c[15] * 12.011000; // C
  W += c[16] * 13.019000; // CH
  W += c[17] * 32.042000; // CH3OH
  W += c[18] * 31.034000; // CH3O
  W += c[19] * 31.034000; // CH2OH
  W += c[20] * 30.026000; // CH2O
  W += c[21] * 29.018000; // HCO
  W += c[22] * 29.017451; // CHO+
  W += c[23] * 19.022451; // H3O+
  W += c[24] * 31.998549; // O2-
  W += c[25] * 0.000549;  // E

  for (id = 0; id < 26; ++id) {
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
  amrex::Real tmp[26];
  amrex::Real imw[26];

  get_imw(imw);

  for (int i = 0; i < 26; i++) {
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
  XW += x[0] * 28.014000;  // N2
  XW += x[1] * 2.016000;   // H2
  XW += x[2] * 1.008000;   // H
  XW += x[3] * 31.998000;  // O2
  XW += x[4] * 15.999000;  // O
  XW += x[5] * 18.015000;  // H2O
  XW += x[6] * 17.007000;  // OH
  XW += x[7] * 34.014000;  // H2O2
  XW += x[8] * 33.006000;  // HO2
  XW += x[9] * 28.010000;  // CO
  XW += x[10] * 44.009000; // CO2
  XW += x[11] * 16.043000; // CH4
  XW += x[12] * 15.035000; // CH3
  XW += x[13] * 14.027000; // CH2
  XW += x[14] * 14.027000; // CH2(S)
  XW += x[15] * 12.011000; // C
  XW += x[16] * 13.019000; // CH
  XW += x[17] * 32.042000; // CH3OH
  XW += x[18] * 31.034000; // CH3O
  XW += x[19] * 31.034000; // CH2OH
  XW += x[20] * 30.026000; // CH2O
  XW += x[21] * 29.018000; // HCO
  XW += x[22] * 29.017451; // CHO+
  XW += x[23] * 19.022451; // H3O+
  XW += x[24] * 31.998549; // O2-
  XW += x[25] * 0.000549;  // E
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
  W += c[0] * 28.014000;  // N2
  W += c[1] * 2.016000;   // H2
  W += c[2] * 1.008000;   // H
  W += c[3] * 31.998000;  // O2
  W += c[4] * 15.999000;  // O
  W += c[5] * 18.015000;  // H2O
  W += c[6] * 17.007000;  // OH
  W += c[7] * 34.014000;  // H2O2
  W += c[8] * 33.006000;  // HO2
  W += c[9] * 28.010000;  // CO
  W += c[10] * 44.009000; // CO2
  W += c[11] * 16.043000; // CH4
  W += c[12] * 15.035000; // CH3
  W += c[13] * 14.027000; // CH2
  W += c[14] * 14.027000; // CH2(S)
  W += c[15] * 12.011000; // C
  W += c[16] * 13.019000; // CH
  W += c[17] * 32.042000; // CH3OH
  W += c[18] * 31.034000; // CH3O
  W += c[19] * 31.034000; // CH2OH
  W += c[20] * 30.026000; // CH2O
  W += c[21] * 29.018000; // HCO
  W += c[22] * 29.017451; // CHO+
  W += c[23] * 19.022451; // H3O+
  W += c[24] * 31.998549; // O2-
  W += c[25] * 0.000549;  // E

  for (id = 0; id < 26; ++id) {
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
  amrex::Real tmp[26];
  amrex::Real imw[26];

  get_imw(imw);

  for (int i = 0; i < 26; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (const auto& tp : tmp) {
    YOW += tp;
  }

  amrex::Real YOWINV = 1.0 / YOW;

  for (int i = 0; i < 26; i++) {
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
  amrex::Real imw[26];

  get_imw(imw);

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 26; i++) {
    c[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 26; i++) {
    YOW += c[i];
  }

  // PW/RT (see Eq. 7)
  PWORT = (*P) / (YOW * 8.31446261815324e+07 * (*T));
  // Now compute conversion

  for (int i = 0; i < 26; i++) {
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
  amrex::Real imw[26];

  get_imw(imw);

  for (int i = 0; i < 26; i++) {
    c[i] = (*rho) * y[i] * imw[i];
  }
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTY(const amrex::Real* x, amrex::Real* y)
{
  amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
  // Compute mean molecular wt first
  XW += x[0] * 28.014000;  // N2
  XW += x[1] * 2.016000;   // H2
  XW += x[2] * 1.008000;   // H
  XW += x[3] * 31.998000;  // O2
  XW += x[4] * 15.999000;  // O
  XW += x[5] * 18.015000;  // H2O
  XW += x[6] * 17.007000;  // OH
  XW += x[7] * 34.014000;  // H2O2
  XW += x[8] * 33.006000;  // HO2
  XW += x[9] * 28.010000;  // CO
  XW += x[10] * 44.009000; // CO2
  XW += x[11] * 16.043000; // CH4
  XW += x[12] * 15.035000; // CH3
  XW += x[13] * 14.027000; // CH2
  XW += x[14] * 14.027000; // CH2(S)
  XW += x[15] * 12.011000; // C
  XW += x[16] * 13.019000; // CH
  XW += x[17] * 32.042000; // CH3OH
  XW += x[18] * 31.034000; // CH3O
  XW += x[19] * 31.034000; // CH2OH
  XW += x[20] * 30.026000; // CH2O
  XW += x[21] * 29.018000; // HCO
  XW += x[22] * 29.017451; // CHO+
  XW += x[23] * 19.022451; // H3O+
  XW += x[24] * 31.998549; // O2-
  XW += x[25] * 0.000549;  // E
  // Now compute conversion
  amrex::Real XWinv = 1.0 / XW;
  y[0] = x[0] * 28.014000 * XWinv;
  y[1] = x[1] * 2.016000 * XWinv;
  y[2] = x[2] * 1.008000 * XWinv;
  y[3] = x[3] * 31.998000 * XWinv;
  y[4] = x[4] * 15.999000 * XWinv;
  y[5] = x[5] * 18.015000 * XWinv;
  y[6] = x[6] * 17.007000 * XWinv;
  y[7] = x[7] * 34.014000 * XWinv;
  y[8] = x[8] * 33.006000 * XWinv;
  y[9] = x[9] * 28.010000 * XWinv;
  y[10] = x[10] * 44.009000 * XWinv;
  y[11] = x[11] * 16.043000 * XWinv;
  y[12] = x[12] * 15.035000 * XWinv;
  y[13] = x[13] * 14.027000 * XWinv;
  y[14] = x[14] * 14.027000 * XWinv;
  y[15] = x[15] * 12.011000 * XWinv;
  y[16] = x[16] * 13.019000 * XWinv;
  y[17] = x[17] * 32.042000 * XWinv;
  y[18] = x[18] * 31.034000 * XWinv;
  y[19] = x[19] * 31.034000 * XWinv;
  y[20] = x[20] * 30.026000 * XWinv;
  y[21] = x[21] * 29.018000 * XWinv;
  y[22] = x[22] * 29.017451 * XWinv;
  y[23] = x[23] * 19.022451 * XWinv;
  y[24] = x[24] * 31.998549 * XWinv;
  y[25] = x[25] * 0.000549 * XWinv;
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
  for (id = 0; id < 26; ++id) {
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
  XW += x[0] * 28.014000;  // N2
  XW += x[1] * 2.016000;   // H2
  XW += x[2] * 1.008000;   // H
  XW += x[3] * 31.998000;  // O2
  XW += x[4] * 15.999000;  // O
  XW += x[5] * 18.015000;  // H2O
  XW += x[6] * 17.007000;  // OH
  XW += x[7] * 34.014000;  // H2O2
  XW += x[8] * 33.006000;  // HO2
  XW += x[9] * 28.010000;  // CO
  XW += x[10] * 44.009000; // CO2
  XW += x[11] * 16.043000; // CH4
  XW += x[12] * 15.035000; // CH3
  XW += x[13] * 14.027000; // CH2
  XW += x[14] * 14.027000; // CH2(S)
  XW += x[15] * 12.011000; // C
  XW += x[16] * 13.019000; // CH
  XW += x[17] * 32.042000; // CH3OH
  XW += x[18] * 31.034000; // CH3O
  XW += x[19] * 31.034000; // CH2OH
  XW += x[20] * 30.026000; // CH2O
  XW += x[21] * 29.018000; // HCO
  XW += x[22] * 29.017451; // CHO+
  XW += x[23] * 19.022451; // H3O+
  XW += x[24] * 31.998549; // O2-
  XW += x[25] * 0.000549;  // E
  ROW = (*rho) / XW;

  // Compute conversion, see Eq 11
  for (id = 0; id < 26; ++id) {
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
  for (id = 0; id < 26; ++id) {
    sumC += c[id];
  }

  //  See Eq 13
  amrex::Real sumCinv = 1.0 / sumC;
  for (id = 0; id < 26; ++id) {
    x[id] = c[id] * sumCinv;
  }
}

// convert c[species] (molar conc) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTY(const amrex::Real* c, amrex::Real* y)
{
  amrex::Real CW = 0; // See Eq 12 in CK Manual
  // compute denominator in eq 12 first
  CW += c[0] * 28.014000;  // N2
  CW += c[1] * 2.016000;   // H2
  CW += c[2] * 1.008000;   // H
  CW += c[3] * 31.998000;  // O2
  CW += c[4] * 15.999000;  // O
  CW += c[5] * 18.015000;  // H2O
  CW += c[6] * 17.007000;  // OH
  CW += c[7] * 34.014000;  // H2O2
  CW += c[8] * 33.006000;  // HO2
  CW += c[9] * 28.010000;  // CO
  CW += c[10] * 44.009000; // CO2
  CW += c[11] * 16.043000; // CH4
  CW += c[12] * 15.035000; // CH3
  CW += c[13] * 14.027000; // CH2
  CW += c[14] * 14.027000; // CH2(S)
  CW += c[15] * 12.011000; // C
  CW += c[16] * 13.019000; // CH
  CW += c[17] * 32.042000; // CH3OH
  CW += c[18] * 31.034000; // CH3O
  CW += c[19] * 31.034000; // CH2OH
  CW += c[20] * 30.026000; // CH2O
  CW += c[21] * 29.018000; // HCO
  CW += c[22] * 29.017451; // CHO+
  CW += c[23] * 19.022451; // H3O+
  CW += c[24] * 31.998549; // O2-
  CW += c[25] * 0.000549;  // E
  // Now compute conversion
  amrex::Real CWinv = 1.0 / CW;
  y[0] = c[0] * 28.014000 * CWinv;
  y[1] = c[1] * 2.016000 * CWinv;
  y[2] = c[2] * 1.008000 * CWinv;
  y[3] = c[3] * 31.998000 * CWinv;
  y[4] = c[4] * 15.999000 * CWinv;
  y[5] = c[5] * 18.015000 * CWinv;
  y[6] = c[6] * 17.007000 * CWinv;
  y[7] = c[7] * 34.014000 * CWinv;
  y[8] = c[8] * 33.006000 * CWinv;
  y[9] = c[9] * 28.010000 * CWinv;
  y[10] = c[10] * 44.009000 * CWinv;
  y[11] = c[11] * 16.043000 * CWinv;
  y[12] = c[12] * 15.035000 * CWinv;
  y[13] = c[13] * 14.027000 * CWinv;
  y[14] = c[14] * 14.027000 * CWinv;
  y[15] = c[15] * 12.011000 * CWinv;
  y[16] = c[16] * 13.019000 * CWinv;
  y[17] = c[17] * 32.042000 * CWinv;
  y[18] = c[18] * 31.034000 * CWinv;
  y[19] = c[19] * 31.034000 * CWinv;
  y[20] = c[20] * 30.026000 * CWinv;
  y[21] = c[21] * 29.018000 * CWinv;
  y[22] = c[22] * 29.017451 * CWinv;
  y[23] = c[23] * 19.022451 * CWinv;
  y[24] = c[24] * 31.998549 * CWinv;
  y[25] = c[25] * 0.000549 * CWinv;
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
  for (id = 0; id < 26; ++id) {
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
  for (id = 0; id < 26; ++id) {
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
  for (id = 0; id < 26; ++id) {
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
  for (id = 0; id < 26; ++id) {
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
  for (id = 0; id < 26; ++id) {
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
  cvms[0] *= 2.967966951578939e+06;  // N2
  cvms[1] *= 4.124237409798234e+07;  // H2
  cvms[2] *= 8.248474819596468e+07;  // H
  cvms[3] *= 2.598431970171023e+06;  // O2
  cvms[4] *= 5.196863940342046e+06;  // O
  cvms[5] *= 4.615299815794193e+06;  // H2O
  cvms[6] *= 4.888847308845322e+06;  // OH
  cvms[7] *= 2.444423654422661e+06;  // H2O2
  cvms[8] *= 2.519076112874398e+06;  // HO2
  cvms[9] *= 2.968390795484913e+06;  // CO
  cvms[10] *= 1.889264154639560e+06; // CO2
  cvms[11] *= 5.182610869633635e+06; // CH4
  cvms[12] *= 5.530071578419182e+06; // CH3
  cvms[13] *= 5.927470320206203e+06; // CH2
  cvms[14] *= 5.927470320206203e+06; // CH2(S)
  cvms[15] *= 6.922373339566431e+06; // C
  cvms[16] *= 6.386406496776436e+06; // CH
  cvms[17] *= 2.594863809423020e+06; // CH3OH
  cvms[18] *= 2.679146297013998e+06; // CH3O
  cvms[19] *= 2.679146297013998e+06; // CH2OH
  cvms[20] *= 2.769087663409458e+06; // CH2O
  cvms[21] *= 2.865277627042952e+06; // HCO
  cvms[22] *= 2.865331795609199e+06; // CHO+
  cvms[23] *= 4.370868104502911e+06; // H3O+
  cvms[24] *= 2.598387422913829e+06; // O2-
  cvms[25] *= 1.515633818095350e+11; // E
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
  cpms[0] *= 2.967966951578939e+06;  // N2
  cpms[1] *= 4.124237409798234e+07;  // H2
  cpms[2] *= 8.248474819596468e+07;  // H
  cpms[3] *= 2.598431970171023e+06;  // O2
  cpms[4] *= 5.196863940342046e+06;  // O
  cpms[5] *= 4.615299815794193e+06;  // H2O
  cpms[6] *= 4.888847308845322e+06;  // OH
  cpms[7] *= 2.444423654422661e+06;  // H2O2
  cpms[8] *= 2.519076112874398e+06;  // HO2
  cpms[9] *= 2.968390795484913e+06;  // CO
  cpms[10] *= 1.889264154639560e+06; // CO2
  cpms[11] *= 5.182610869633635e+06; // CH4
  cpms[12] *= 5.530071578419182e+06; // CH3
  cpms[13] *= 5.927470320206203e+06; // CH2
  cpms[14] *= 5.927470320206203e+06; // CH2(S)
  cpms[15] *= 6.922373339566431e+06; // C
  cpms[16] *= 6.386406496776436e+06; // CH
  cpms[17] *= 2.594863809423020e+06; // CH3OH
  cpms[18] *= 2.679146297013998e+06; // CH3O
  cpms[19] *= 2.679146297013998e+06; // CH2OH
  cpms[20] *= 2.769087663409458e+06; // CH2O
  cpms[21] *= 2.865277627042952e+06; // HCO
  cpms[22] *= 2.865331795609199e+06; // CHO+
  cpms[23] *= 4.370868104502911e+06; // H3O+
  cpms[24] *= 2.598387422913829e+06; // O2-
  cpms[25] *= 1.515633818095350e+11; // E
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUMS(const amrex::Real* T, amrex::Real* ums)
{
  amrex::Real tT = *T; // temporary temperature
  const amrex::Real tc[5] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;         // R*T
  amrex::Real imw[26];

  get_imw(imw);

  speciesInternalEnergy(ums, tc);

  for (int i = 0; i < 26; i++) {
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
  amrex::Real imw[26];

  get_imw(imw);

  speciesEnthalpy(hms, tc);

  for (int i = 0; i < 26; i++) {
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
  sms[0] *= 2.967966951578939e+06;  // N2
  sms[1] *= 4.124237409798234e+07;  // H2
  sms[2] *= 8.248474819596468e+07;  // H
  sms[3] *= 2.598431970171023e+06;  // O2
  sms[4] *= 5.196863940342046e+06;  // O
  sms[5] *= 4.615299815794193e+06;  // H2O
  sms[6] *= 4.888847308845322e+06;  // OH
  sms[7] *= 2.444423654422661e+06;  // H2O2
  sms[8] *= 2.519076112874398e+06;  // HO2
  sms[9] *= 2.968390795484913e+06;  // CO
  sms[10] *= 1.889264154639560e+06; // CO2
  sms[11] *= 5.182610869633635e+06; // CH4
  sms[12] *= 5.530071578419182e+06; // CH3
  sms[13] *= 5.927470320206203e+06; // CH2
  sms[14] *= 5.927470320206203e+06; // CH2(S)
  sms[15] *= 6.922373339566431e+06; // C
  sms[16] *= 6.386406496776436e+06; // CH
  sms[17] *= 2.594863809423020e+06; // CH3OH
  sms[18] *= 2.679146297013998e+06; // CH3O
  sms[19] *= 2.679146297013998e+06; // CH2OH
  sms[20] *= 2.769087663409458e+06; // CH2O
  sms[21] *= 2.865277627042952e+06; // HCO
  sms[22] *= 2.865331795609199e+06; // CHO+
  sms[23] *= 4.370868104502911e+06; // H3O+
  sms[24] *= 2.598387422913829e+06; // O2-
  sms[25] *= 1.515633818095350e+11; // E
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

  // reaction 6: H2O2 (+M) <=> 2 OH (+M)
  qf[0] = sc[7];
  qr[0] = (sc[6] * sc[6]);

  // reaction 19: H + O2 (+M) <=> HO2 (+M)
  qf[1] = sc[2] * sc[3];
  qr[1] = sc[8];

  // reaction 24: CH3 + H (+M) <=> CH4 (+M)
  qf[2] = sc[2] * sc[12];
  qr[2] = sc[11];

  // reaction 43: CH2 + H (+M) <=> CH3 (+M)
  qf[3] = sc[2] * sc[13];
  qr[3] = sc[12];

  // reaction 66: CH3OH (+M) <=> CH3 + OH (+M)
  qf[4] = sc[17];
  qr[4] = sc[6] * sc[12];

  // reaction 67: CH3OH (+M) <=> CH2(S) + H2O (+M)
  qf[5] = sc[17];
  qr[5] = sc[5] * sc[14];

  // reaction 106: CH2O + H (+M) <=> CH2OH (+M)
  qf[6] = sc[2] * sc[20];
  qr[6] = sc[19];

  // reaction 107: CH3O (+M) <=> CH2O + H (+M)
  qf[7] = sc[18];
  qr[7] = sc[2] * sc[20];

  // reaction 20: CO + O (+M) <=> CO2 (+M)
  qf[8] = sc[4] * sc[9];
  qr[8] = sc[10];

  // reaction 3: H + OH + M <=> H2O + M
  qf[9] = sc[2] * sc[6];
  qr[9] = sc[5];

  // reaction 5: H + O + M <=> OH + M
  qf[10] = sc[2] * sc[4];
  qr[10] = sc[6];

  // reaction 97: HCO + M <=> CO + H + M
  qf[11] = sc[21];
  qr[11] = sc[2] * sc[9];

  // reaction 108: O + OH + M <=> HO2 + M
  qf[12] = sc[4] * sc[6];
  qr[12] = sc[8];

  // reaction 120: CH3 + O + M => CH3O + M
  qf[13] = sc[4] * sc[12];
  qr[13] = 0.0;

  // reaction 131: E + O2 + O2 <=> O2- + O2
  qf[14] = (sc[3] * sc[3]) * sc[25];
  qr[14] = sc[3] * sc[24];

  // reaction 132: E + O2 + H2O <=> O2- + H2O
  qf[15] = sc[3] * sc[5] * sc[25];
  qr[15] = sc[5] * sc[24];

  // reaction 133: E + O2 + N2 <=> O2- + N2
  qf[16] = sc[0] * sc[3] * sc[25];
  qr[16] = sc[0] * sc[24];

  // reaction 0: H2 + O <=> H + OH
  qf[17] = sc[1] * sc[4];
  qr[17] = sc[2] * sc[6];

  // reaction 1: H2 + OH <=> H + H2O
  qf[18] = sc[1] * sc[6];
  qr[18] = sc[2] * sc[5];

  // reaction 2: H + O2 <=> O + OH
  qf[19] = sc[2] * sc[3];
  qr[19] = sc[4] * sc[6];

  // reaction 4: H2O + O <=> 2 OH
  qf[20] = sc[4] * sc[5];
  qr[20] = (sc[6] * sc[6]);

  // reaction 7: H + H2O2 <=> H2O + OH
  qf[21] = sc[2] * sc[7];
  qr[21] = sc[5] * sc[6];

  // reaction 8: H + H2O2 <=> H2 + HO2
  qf[22] = sc[2] * sc[7];
  qr[22] = sc[1] * sc[8];

  // reaction 9: H2O2 + O <=> HO2 + OH
  qf[23] = sc[4] * sc[7];
  qr[23] = sc[6] * sc[8];

  // reaction 10: H2O2 + OH <=> H2O + HO2
  qf[24] = sc[6] * sc[7];
  qr[24] = sc[5] * sc[8];

  // reaction 11: H2O2 + OH <=> H2O + HO2
  qf[25] = sc[6] * sc[7];
  qr[25] = sc[5] * sc[8];

  // reaction 12: H + HO2 <=> 2 OH
  qf[26] = sc[2] * sc[8];
  qr[26] = (sc[6] * sc[6]);

  // reaction 13: H + HO2 <=> H2 + O2
  qf[27] = sc[2] * sc[8];
  qr[27] = sc[1] * sc[3];

  // reaction 14: HO2 + O <=> O2 + OH
  qf[28] = sc[4] * sc[8];
  qr[28] = sc[3] * sc[6];

  // reaction 15: HO2 + OH <=> H2O + O2
  qf[29] = sc[6] * sc[8];
  qr[29] = sc[3] * sc[5];

  // reaction 16: HO2 + OH <=> H2O + O2
  qf[30] = sc[6] * sc[8];
  qr[30] = sc[3] * sc[5];

  // reaction 17: 2 HO2 <=> H2O2 + O2
  qf[31] = (sc[8] * sc[8]);
  qr[31] = sc[3] * sc[7];

  // reaction 18: 2 HO2 <=> H2O2 + O2
  qf[32] = (sc[8] * sc[8]);
  qr[32] = sc[3] * sc[7];

  // reaction 21: CO + OH <=> CO2 + H
  qf[33] = sc[6] * sc[9];
  qr[33] = sc[2] * sc[10];

  // reaction 22: CO + OH <=> CO2 + H
  qf[34] = sc[6] * sc[9];
  qr[34] = sc[2] * sc[10];

  // reaction 23: CO + HO2 <=> CO2 + OH
  qf[35] = sc[8] * sc[9];
  qr[35] = sc[6] * sc[10];

  // reaction 25: CH4 + H <=> CH3 + H2
  qf[36] = sc[2] * sc[11];
  qr[36] = sc[1] * sc[12];

  // reaction 26: CH4 + O <=> CH3 + OH
  qf[37] = sc[4] * sc[11];
  qr[37] = sc[6] * sc[12];

  // reaction 27: CH4 + OH <=> CH3 + H2O
  qf[38] = sc[6] * sc[11];
  qr[38] = sc[5] * sc[12];

  // reaction 28: CH4 + HO2 <=> CH3 + H2O2
  qf[39] = sc[8] * sc[11];
  qr[39] = sc[7] * sc[12];

  // reaction 29: CH3 + HO2 <=> CH4 + O2
  qf[40] = sc[8] * sc[12];
  qr[40] = sc[3] * sc[11];

  // reaction 30: CH2 + CH4 <=> 2 CH3
  qf[41] = sc[11] * sc[13];
  qr[41] = (sc[12] * sc[12]);

  // reaction 31: CH2(S) + N2 <=> CH2 + N2
  qf[42] = sc[0] * sc[14];
  qr[42] = sc[0] * sc[13];

  // reaction 32: CH2(S) + H2O <=> CH2 + H2O
  qf[43] = sc[5] * sc[14];
  qr[43] = sc[5] * sc[13];

  // reaction 33: CH2(S) + CO <=> CH2 + CO
  qf[44] = sc[9] * sc[14];
  qr[44] = sc[9] * sc[13];

  // reaction 34: CH2(S) + CO2 <=> CH2 + CO2
  qf[45] = sc[10] * sc[14];
  qr[45] = sc[10] * sc[13];

  // reaction 35: CH2(S) + O2 => CO + H + OH
  qf[46] = sc[3] * sc[14];
  qr[46] = 0.0;

  // reaction 36: CH2(S) + O2 <=> CO + H2O
  qf[47] = sc[3] * sc[14];
  qr[47] = sc[5] * sc[9];

  // reaction 37: CH2(S) + O <=> CO + H2
  qf[48] = sc[4] * sc[14];
  qr[48] = sc[1] * sc[9];

  // reaction 38: CH2(S) + O <=> H + HCO
  qf[49] = sc[4] * sc[14];
  qr[49] = sc[2] * sc[21];

  // reaction 39: CH2(S) + H2 <=> CH3 + H
  qf[50] = sc[1] * sc[14];
  qr[50] = sc[2] * sc[12];

  // reaction 40: CH2(S) + H <=> CH + H2
  qf[51] = sc[2] * sc[14];
  qr[51] = sc[1] * sc[16];

  // reaction 41: CH2(S) + OH <=> CH2O + H
  qf[52] = sc[6] * sc[14];
  qr[52] = sc[2] * sc[20];

  // reaction 42: CH2(S) + CO2 <=> CH2O + CO
  qf[53] = sc[10] * sc[14];
  qr[53] = sc[9] * sc[20];

  // reaction 44: CH2 + O2 <=> HCO + OH
  qf[54] = sc[3] * sc[13];
  qr[54] = sc[6] * sc[21];

  // reaction 45: CH2 + O2 => CO2 + 2 H
  qf[55] = sc[3] * sc[13];
  qr[55] = 0.0;

  // reaction 46: CH2 + O => CO + 2 H
  qf[56] = sc[4] * sc[13];
  qr[56] = 0.0;

  // reaction 47: CH2 + H <=> CH + H2
  qf[57] = sc[2] * sc[13];
  qr[57] = sc[1] * sc[16];

  // reaction 48: CH2 + OH <=> CH + H2O
  qf[58] = sc[6] * sc[13];
  qr[58] = sc[5] * sc[16];

  // reaction 49: CH + O2 <=> HCO + O
  qf[59] = sc[3] * sc[16];
  qr[59] = sc[4] * sc[21];

  // reaction 50: CH + O <=> CO + H
  qf[60] = sc[4] * sc[16];
  qr[60] = sc[2] * sc[9];

  // reaction 51: CH + H <=> C + H2
  qf[61] = sc[2] * sc[16];
  qr[61] = sc[1] * sc[15];

  // reaction 52: CH + OH <=> H + HCO
  qf[62] = sc[6] * sc[16];
  qr[62] = sc[2] * sc[21];

  // reaction 53: CH + H2O <=> CH2O + H
  qf[63] = sc[5] * sc[16];
  qr[63] = sc[2] * sc[20];

  // reaction 54: CH + CO2 <=> CO + HCO
  qf[64] = sc[10] * sc[16];
  qr[64] = sc[9] * sc[21];

  // reaction 55: C + OH <=> CO + H
  qf[65] = sc[6] * sc[15];
  qr[65] = sc[2] * sc[9];

  // reaction 56: C + O2 <=> CO + O
  qf[66] = sc[3] * sc[15];
  qr[66] = sc[4] * sc[9];

  // reaction 57: CH3 + O2 <=> CH3O + O
  qf[67] = sc[3] * sc[12];
  qr[67] = sc[4] * sc[18];

  // reaction 58: CH3 + O2 <=> CH2O + OH
  qf[68] = sc[3] * sc[12];
  qr[68] = sc[6] * sc[20];

  // reaction 59: CH3 + O <=> CH2O + H
  qf[69] = sc[4] * sc[12];
  qr[69] = sc[2] * sc[20];

  // reaction 60: CH3 + OH <=> CH2(S) + H2O
  qf[70] = sc[6] * sc[12];
  qr[70] = sc[5] * sc[14];

  // reaction 61: CH3 + OH <=> CH2O + H2
  qf[71] = sc[6] * sc[12];
  qr[71] = sc[1] * sc[20];

  // reaction 62: CH3 + OH <=> CH2OH + H
  qf[72] = sc[6] * sc[12];
  qr[72] = sc[2] * sc[19];

  // reaction 63: CH3 + OH <=> CH3O + H
  qf[73] = sc[6] * sc[12];
  qr[73] = sc[2] * sc[18];

  // reaction 64: CH3 + OH <=> CH2 + H2O
  qf[74] = sc[6] * sc[12];
  qr[74] = sc[5] * sc[13];

  // reaction 65: CH3 + HO2 <=> CH3O + OH
  qf[75] = sc[8] * sc[12];
  qr[75] = sc[6] * sc[18];

  // reaction 68: CH3OH + H <=> CH3O + H2
  qf[76] = sc[2] * sc[17];
  qr[76] = sc[1] * sc[18];

  // reaction 69: CH3OH + H <=> CH2OH + H2
  qf[77] = sc[2] * sc[17];
  qr[77] = sc[1] * sc[19];

  // reaction 70: CH3OH + O <=> CH3O + OH
  qf[78] = sc[4] * sc[17];
  qr[78] = sc[6] * sc[18];

  // reaction 71: CH3OH + O <=> CH2OH + OH
  qf[79] = sc[4] * sc[17];
  qr[79] = sc[6] * sc[19];

  // reaction 72: CH3OH + OH <=> CH3O + H2O
  qf[80] = sc[6] * sc[17];
  qr[80] = sc[5] * sc[18];

  // reaction 73: CH3OH + OH <=> CH2OH + H2O
  qf[81] = sc[6] * sc[17];
  qr[81] = sc[5] * sc[19];

  // reaction 74: CH3OH + O2 <=> CH3O + HO2
  qf[82] = sc[3] * sc[17];
  qr[82] = sc[8] * sc[18];

  // reaction 75: CH3OH + HO2 <=> CH3O + H2O2
  qf[83] = sc[8] * sc[17];
  qr[83] = sc[7] * sc[18];

  // reaction 76: CH3OH + HO2 <=> CH2OH + H2O2
  qf[84] = sc[8] * sc[17];
  qr[84] = sc[7] * sc[19];

  // reaction 77: CH3 + CH3OH <=> CH2OH + CH4
  qf[85] = sc[12] * sc[17];
  qr[85] = sc[11] * sc[19];

  // reaction 78: CH3 + CH3OH <=> CH3O + CH4
  qf[86] = sc[12] * sc[17];
  qr[86] = sc[11] * sc[18];

  // reaction 79: CH2OH + O2 <=> CH2O + HO2
  qf[87] = sc[3] * sc[19];
  qr[87] = sc[8] * sc[20];

  // reaction 80: CH2OH + O2 <=> CH2O + HO2
  qf[88] = sc[3] * sc[19];
  qr[88] = sc[8] * sc[20];

  // reaction 81: CH2OH + H <=> CH2O + H2
  qf[89] = sc[2] * sc[19];
  qr[89] = sc[1] * sc[20];

  // reaction 82: CH2OH + HO2 <=> CH2O + H2O2
  qf[90] = sc[8] * sc[19];
  qr[90] = sc[7] * sc[20];

  // reaction 83: CH2OH + HCO <=> 2 CH2O
  qf[91] = sc[19] * sc[21];
  qr[91] = (sc[20] * sc[20]);

  // reaction 84: CH2OH + OH <=> CH2O + H2O
  qf[92] = sc[6] * sc[19];
  qr[92] = sc[5] * sc[20];

  // reaction 85: CH2OH + O <=> CH2O + OH
  qf[93] = sc[4] * sc[19];
  qr[93] = sc[6] * sc[20];

  // reaction 86: CH3O + O2 <=> CH2O + HO2
  qf[94] = sc[3] * sc[18];
  qr[94] = sc[8] * sc[20];

  // reaction 87: CH3O + H <=> CH2O + H2
  qf[95] = sc[2] * sc[18];
  qr[95] = sc[1] * sc[20];

  // reaction 88: CH3O + HO2 <=> CH2O + H2O2
  qf[96] = sc[8] * sc[18];
  qr[96] = sc[7] * sc[20];

  // reaction 89: CH3 + CH3O <=> CH2O + CH4
  qf[97] = sc[12] * sc[18];
  qr[97] = sc[11] * sc[20];

  // reaction 90: 2 CH3O <=> CH2O + CH3OH
  qf[98] = (sc[18] * sc[18]);
  qr[98] = sc[17] * sc[20];

  // reaction 91: CH2O + O <=> HCO + OH
  qf[99] = sc[4] * sc[20];
  qr[99] = sc[6] * sc[21];

  // reaction 92: CH2O + H <=> H2 + HCO
  qf[100] = sc[2] * sc[20];
  qr[100] = sc[1] * sc[21];

  // reaction 93: CH2O + OH <=> H2O + HCO
  qf[101] = sc[6] * sc[20];
  qr[101] = sc[5] * sc[21];

  // reaction 94: CH2O + HO2 <=> H2O2 + HCO
  qf[102] = sc[8] * sc[20];
  qr[102] = sc[7] * sc[21];

  // reaction 95: CH2O + CH3 <=> CH4 + HCO
  qf[103] = sc[12] * sc[20];
  qr[103] = sc[11] * sc[21];

  // reaction 96: CH2O + CH3O <=> CH3OH + HCO
  qf[104] = sc[18] * sc[20];
  qr[104] = sc[17] * sc[21];

  // reaction 98: HCO + O2 <=> CO + HO2
  qf[105] = sc[3] * sc[21];
  qr[105] = sc[8] * sc[9];

  // reaction 99: HCO + O <=> CO + OH
  qf[106] = sc[4] * sc[21];
  qr[106] = sc[6] * sc[9];

  // reaction 100: H + HCO <=> CO + H2
  qf[107] = sc[2] * sc[21];
  qr[107] = sc[1] * sc[9];

  // reaction 101: HCO + OH <=> CO + H2O
  qf[108] = sc[6] * sc[21];
  qr[108] = sc[5] * sc[9];

  // reaction 102: CH3 + HCO <=> CH4 + CO
  qf[109] = sc[12] * sc[21];
  qr[109] = sc[9] * sc[11];

  // reaction 103: 2 HCO <=> CH2O + CO
  qf[110] = (sc[21] * sc[21]);
  qr[110] = sc[9] * sc[20];

  // reaction 104: HCO + O <=> CO2 + H
  qf[111] = sc[4] * sc[21];
  qr[111] = sc[2] * sc[10];

  // reaction 105: HCO + HO2 => CO2 + H + OH
  qf[112] = sc[8] * sc[21];
  qr[112] = 0.0;

  // reaction 109: CH + OH <=> C + H2O
  qf[113] = sc[6] * sc[16];
  qr[113] = sc[5] * sc[15];

  // reaction 110: CH2 + HCO <=> CH3 + CO
  qf[114] = sc[13] * sc[21];
  qr[114] = sc[9] * sc[12];

  // reaction 111: CH2 + O <=> CO + H2
  qf[115] = sc[4] * sc[13];
  qr[115] = sc[1] * sc[9];

  // reaction 112: CH2 + OH <=> CH2O + H
  qf[116] = sc[6] * sc[13];
  qr[116] = sc[2] * sc[20];

  // reaction 113: CH2 + CO2 <=> CH2O + CO
  qf[117] = sc[10] * sc[13];
  qr[117] = sc[9] * sc[20];

  // reaction 114: CH2(S) + O <=> CO + 2 H
  qf[118] = sc[4] * sc[14];
  qr[118] = (sc[2] * sc[2]) * sc[9];

  // reaction 115: CH2(S) + CH4 <=> 2 CH3
  qf[119] = sc[11] * sc[14];
  qr[119] = (sc[12] * sc[12]);

  // reaction 116: CH3OH + H => CH3 + H2O
  qf[120] = sc[2] * sc[17];
  qr[120] = 0.0;

  // reaction 117: CH2O + O => CO2 + 2 H
  qf[121] = sc[4] * sc[20];
  qr[121] = 0.0;

  // reaction 118: CH2O + OH => CO2 + H + H2
  qf[122] = sc[6] * sc[20];
  qr[122] = 0.0;

  // reaction 119: CH3O + CO <=> CH3 + CO2
  qf[123] = sc[9] * sc[18];
  qr[123] = sc[10] * sc[12];

  // reaction 121: HCO + HO2 <=> CO + H2O2
  qf[124] = sc[8] * sc[21];
  qr[124] = sc[7] * sc[9];

  // reaction 122: CH3O + OH => CH2O + H2O
  qf[125] = sc[6] * sc[18];
  qr[125] = 0.0;

  // reaction 123: CH2OH + CH3 => CH2O + CH4
  qf[126] = sc[12] * sc[19];
  qr[126] = 0.0;

  // reaction 124: CH + O <=> CHO+ + E
  qf[127] = sc[4] * sc[16];
  qr[127] = sc[22] * sc[25];

  // reaction 125: CHO+ + H2O <=> CO + H3O+
  qf[128] = sc[5] * sc[22];
  qr[128] = sc[9] * sc[23];

  // reaction 126: E + H3O+ <=> H + H2O
  qf[129] = sc[23] * sc[25];
  qr[129] = sc[2] * sc[5];

  // reaction 127: E + H3O+ <=> 2 H + OH
  qf[130] = sc[23] * sc[25];
  qr[130] = (sc[2] * sc[2]) * sc[6];

  // reaction 128: E + H3O+ <=> H2 + OH
  qf[131] = sc[23] * sc[25];
  qr[131] = sc[1] * sc[6];

  // reaction 129: H2 + O2- <=> E + H2O2
  qf[132] = sc[1] * sc[24];
  qr[132] = sc[7] * sc[25];

  // reaction 130: H + O2- <=> E + HO2
  qf[133] = sc[2] * sc[24];
  qr[133] = sc[8] * sc[25];

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 26; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[26];
  gibbs(g_RT, tc);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 * invT;
  amrex::Real refCinv = 1 / refC;

  // Evaluate the kfs
  amrex::Real k_f, Corr;
  amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

  // reaction 0:  H2 + O <=> H + OH
  k_f = 0.0508 * exp((2.67) * tc[0] - (3166.2392612433) * invT);
  qf[17] *= k_f;
  qr[17] *= k_f * exp(-(g_RT[1] - g_RT[2] + g_RT[4] - g_RT[6]));
  // reaction 1:  H2 + OH <=> H + H2O
  k_f = 43800000 * exp(-(3517.48449397499) * invT);
  qf[18] *= k_f;
  qr[18] *= k_f * exp(-(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6]));
  // reaction 2:  H + O2 <=> O + OH
  k_f = 114000000 * exp(-(7692.16995349094) * invT);
  qf[19] *= k_f;
  qr[19] *= k_f * exp(-(g_RT[2] + g_RT[3] - g_RT[4] - g_RT[6]));
  // reaction 3:  H + OH + M <=> H2O + M
  k_f = 35000000000 * exp((-2) * tc[0]);
  Corr = mixture + (-0.27) * sc[1] + (2.65) * sc[5] + sc[11];
  qf[9] *= Corr * k_f;
  qr[9] *= Corr * k_f * exp(-(g_RT[2] - g_RT[5] + g_RT[6])) * (refC);
  // reaction 4:  H2O + O <=> 2 OH
  k_f = 67 * exp((1.704) * tc[0] - (7541.60752708217) * invT);
  qf[20] *= k_f;
  qr[20] *= k_f * exp(-(g_RT[4] + g_RT[5] - 2.000000 * g_RT[6]));
  // reaction 5:  H + O + M <=> OH + M
  k_f = 4714000 * exp((-1) * tc[0]);
  Corr =
    mixture + (1.5) * sc[1] + (11) * sc[5] + (0.5) * sc[9] + sc[10] + sc[11];
  qf[10] *= Corr * k_f;
  qr[10] *= Corr * k_f * exp(-(g_RT[2] + g_RT[4] - g_RT[6])) * (refC);
  // reaction 6:  H2O2 (+M) <=> 2 OH (+M)
  k_f = 2000000000000 * exp((0.9) * tc[0] - (24531.3092413143) * invT);
  Corr = mixture + (0.5) * sc[0] + (2.7) * sc[1] + (0.2) * sc[3] +
         (6.65) * sc[5] + (6.7) * sc[7] + (1.8) * sc[9] + (0.6) * sc[10];
  redP = Corr / k_f * 2.49e+18 * exp(-2.3 * tc[0] - (24531.3092413143) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (0.57) * exp(-tc[1] * 1e+30) + 0.43 * exp(-tc[1] * 1e-30) +
    exp(-1e-16 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[0] *= Corr * k_f;
  qr[0] *= Corr * k_f * exp(-(-2.000000 * g_RT[6] + g_RT[7])) * (refCinv);
  // reaction 7:  H + H2O2 <=> H2O + OH
  k_f = 24100000 * exp(-(1997.77016324474) * invT);
  qf[21] *= k_f;
  qr[21] *= k_f * exp(-(g_RT[2] - g_RT[5] - g_RT[6] + g_RT[7]));
  // reaction 8:  H + H2O2 <=> H2 + HO2
  k_f = 21500 * exp((1) * tc[0] - (3019.29999482832) * invT);
  qf[22] *= k_f;
  qr[22] *= k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[8]));
  // reaction 9:  H2O2 + O <=> HO2 + OH
  k_f = 9.55 * exp((2) * tc[0] - (1997.77016324474) * invT);
  qf[23] *= k_f;
  qr[23] *= k_f * exp(-(g_RT[4] - g_RT[6] + g_RT[7] - g_RT[8]));
  // reaction 10:  H2O2 + OH <=> H2O + HO2
  k_f = 1740000 * exp(-(160.022899725901) * invT);
  qf[24] *= k_f;
  qr[24] *= k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[7] - g_RT[8]));
  // reaction 11:  H2O2 + OH <=> H2O + HO2
  k_f = 75900000 * exp(-(3657.88194373451) * invT);
  qf[25] *= k_f;
  qr[25] *= k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[7] - g_RT[8]));
  // reaction 12:  H + HO2 <=> 2 OH
  k_f = 70790000 * exp(-(148.448916412392) * invT);
  qf[26] *= k_f;
  qr[26] *= k_f * exp(-(g_RT[2] - 2.000000 * g_RT[6] + g_RT[8]));
  // reaction 13:  H + HO2 <=> H2 + O2
  k_f = 11402 * exp((1.083) * tc[0] - (278.671325189338) * invT);
  qf[27] *= k_f;
  qr[27] *= k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[8]));
  // reaction 14:  HO2 + O <=> O2 + OH
  k_f = 32500000;
  qf[28] *= k_f;
  qr[28] *= k_f * exp(-(-g_RT[3] + g_RT[4] - g_RT[6] + g_RT[8]));
  // reaction 15:  HO2 + OH <=> H2O + O2
  k_f = 7000000 * exp(-(-549.995687057926) * invT);
  qf[29] *= k_f;
  qr[29] *= k_f * exp(-(-g_RT[3] - g_RT[5] + g_RT[6] + g_RT[8]));
  // reaction 16:  HO2 + OH <=> H2O + O2
  k_f = 450000000 * exp(-(5499.95687057926) * invT);
  qf[30] *= k_f;
  qr[30] *= k_f * exp(-(-g_RT[3] - g_RT[5] + g_RT[6] + g_RT[8]));
  // reaction 17:  2 HO2 <=> H2O2 + O2
  k_f = 100000000 * exp(-(5555.95482115001) * invT);
  qf[31] *= k_f;
  qr[31] *= k_f * exp(-(-g_RT[3] - g_RT[7] + 2.000000 * g_RT[8]));
  // reaction 18:  2 HO2 <=> H2O2 + O2
  k_f = 190000 * exp(-(-708.992024785586) * invT);
  qf[32] *= k_f;
  qr[32] *= k_f * exp(-(-g_RT[3] - g_RT[7] + 2.000000 * g_RT[8]));
  // reaction 19:  H + O2 (+M) <=> HO2 (+M)
  k_f = 4650000 * exp((0.44) * tc[0]);
  Corr = mixture + (0.3) * sc[1] + (9) * sc[5] + (0.9) * sc[9] +
         (2.8) * sc[10] + sc[11];
  redP = Corr / k_f * 17400000 * exp(-1.23 * tc[0]);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (0.33) * exp(-tc[1] * 1e+30) + 0.67 * exp(-tc[1] * 1e-30) +
    exp(-1e+30 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[1] *= Corr * k_f;
  qr[1] *= Corr * k_f * exp(-(g_RT[2] + g_RT[3] - g_RT[8])) * (refC);
  // reaction 20:  CO + O (+M) <=> CO2 (+M)
  k_f = 13620 * exp(-(1199.66853127845) * invT);
  Corr = mixture + sc[1] + (11) * sc[5] + (0.75) * sc[9] + (2.6) * sc[10];
  redP =
    Corr / k_f * 1170000000000 * exp(-2.79 * tc[0] - (2108.98104638758) * invT);
  Corr = redP / (1. + redP);
  qf[8] *= Corr * k_f;
  qr[8] *= Corr * k_f * exp(-(g_RT[4] + g_RT[9] - g_RT[10])) * (refC);
  // reaction 21:  CO + OH <=> CO2 + H
  k_f = 0.07015 * exp((2.053) * tc[0] - (-178.994168026739) * invT);
  qf[33] *= k_f;
  qr[33] *= k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[9] - g_RT[10]));
  // reaction 22:  CO + OH <=> CO2 + H
  k_f = 5757000 * exp((-0.664) * tc[0] - (166.967289714006) * invT);
  qf[34] *= k_f;
  qr[34] *= k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[9] - g_RT[10]));
  // reaction 23:  CO + HO2 <=> CO2 + OH
  k_f = 0.157 * exp((2.18) * tc[0] - (9027.70698453667) * invT);
  qf[35] *= k_f;
  qr[35] *= k_f * exp(-(-g_RT[6] + g_RT[8] + g_RT[9] - g_RT[10]));
  // reaction 24:  CH3 + H (+M) <=> CH4 (+M)
  k_f = 12700000000 * exp((-0.63) * tc[0] - (192.731983003208) * invT);
  Corr = mixture + sc[1] + (5) * sc[5] + (0.5) * sc[9] + sc[10] + sc[11];
  redP = Corr / k_f * 2.48e+21 * exp(-4.76 * tc[0] - (1227.84866456352) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (0.217) * exp(-tc[1] * 0.0135135135135135) +
    0.783 * exp(-tc[1] * 0.000340020401224073) + exp(-6964 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[2] *= Corr * k_f;
  qr[2] *= Corr * k_f * exp(-(g_RT[2] - g_RT[11] + g_RT[12])) * (refC);
  // reaction 25:  CH4 + H <=> CH3 + H2
  k_f = 0.614 * exp((2.5) * tc[0] - (4824.33817506985) * invT);
  qf[36] *= k_f;
  qr[36] *= k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[11] - g_RT[12]));
  // reaction 26:  CH4 + O <=> CH3 + OH
  k_f = 1020 * exp((1.5) * tc[0] - (4327.66332592059) * invT);
  qf[37] *= k_f;
  qr[37] *= k_f * exp(-(g_RT[4] - g_RT[6] + g_RT[11] - g_RT[12]));
  // reaction 27:  CH4 + OH <=> CH3 + H2O
  k_f = 0.0583 * exp((2.6) * tc[0] - (1102.04449811234) * invT);
  qf[38] *= k_f;
  qr[38] *= k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[11] - g_RT[12]));
  // reaction 28:  CH4 + HO2 <=> CH3 + H2O2
  k_f = 1.13e-05 * exp((3.74) * tc[0] - (10572.5821485572) * invT);
  qf[39] *= k_f;
  qr[39] *= k_f * exp(-(-g_RT[7] + g_RT[8] + g_RT[11] - g_RT[12]));
  // reaction 29:  CH3 + HO2 <=> CH4 + O2
  k_f = 0.116 * exp((2.23) * tc[0] - (-1520.72076406186) * invT);
  qf[40] *= k_f;
  qr[40] *= k_f * exp(-(-g_RT[3] + g_RT[8] - g_RT[11] + g_RT[12]));
  // reaction 30:  CH2 + CH4 <=> 2 CH3
  k_f = 2.46 * exp((2) * tc[0] - (4161.60182620503) * invT);
  qf[41] *= k_f;
  qr[41] *= k_f * exp(-(g_RT[11] - 2.000000 * g_RT[12] + g_RT[13]));
  // reaction 31:  CH2(S) + N2 <=> CH2 + N2
  k_f = 15000000 * exp(-(301.929999482832) * invT);
  qf[42] *= k_f;
  qr[42] *= k_f * exp(-(g_RT[0] - g_RT[0] - g_RT[13] + g_RT[14]));
  // reaction 32:  CH2(S) + H2O <=> CH2 + H2O
  k_f = 30000000;
  qf[43] *= k_f;
  qr[43] *= k_f * exp(-(g_RT[5] - g_RT[5] - g_RT[13] + g_RT[14]));
  // reaction 33:  CH2(S) + CO <=> CH2 + CO
  k_f = 9000000;
  qf[44] *= k_f;
  qr[44] *= k_f * exp(-(g_RT[9] - g_RT[9] - g_RT[13] + g_RT[14]));
  // reaction 34:  CH2(S) + CO2 <=> CH2 + CO2
  k_f = 7000000;
  qf[45] *= k_f;
  qr[45] *= k_f * exp(-(g_RT[10] - g_RT[10] - g_RT[13] + g_RT[14]));
  // reaction 35:  CH2(S) + O2 => CO + H + OH
  k_f = 28000000;
  qf[46] *= k_f;
  qr[46] *=
    k_f * exp(-(-g_RT[2] + g_RT[3] - g_RT[6] - g_RT[9] + g_RT[14])) * (refCinv);
  // reaction 36:  CH2(S) + O2 <=> CO + H2O
  k_f = 12000000;
  qf[47] *= k_f;
  qr[47] *= k_f * exp(-(g_RT[3] - g_RT[5] - g_RT[9] + g_RT[14]));
  // reaction 37:  CH2(S) + O <=> CO + H2
  k_f = 15000000;
  qf[48] *= k_f;
  qr[48] *= k_f * exp(-(-g_RT[1] + g_RT[4] - g_RT[9] + g_RT[14]));
  // reaction 38:  CH2(S) + O <=> H + HCO
  k_f = 15000000;
  qf[49] *= k_f;
  qr[49] *= k_f * exp(-(-g_RT[2] + g_RT[4] + g_RT[14] - g_RT[21]));
  // reaction 39:  CH2(S) + H2 <=> CH3 + H
  k_f = 70000000;
  qf[50] *= k_f;
  qr[50] *= k_f * exp(-(g_RT[1] - g_RT[2] - g_RT[12] + g_RT[14]));
  // reaction 40:  CH2(S) + H <=> CH + H2
  k_f = 30000000;
  qf[51] *= k_f;
  qr[51] *= k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[14] - g_RT[16]));
  // reaction 41:  CH2(S) + OH <=> CH2O + H
  k_f = 30000000;
  qf[52] *= k_f;
  qr[52] *= k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[14] - g_RT[20]));
  // reaction 42:  CH2(S) + CO2 <=> CH2O + CO
  k_f = 14000000;
  qf[53] *= k_f;
  qr[53] *= k_f * exp(-(-g_RT[9] + g_RT[10] + g_RT[14] - g_RT[20]));
  // reaction 43:  CH2 + H (+M) <=> CH3 (+M)
  k_f = 25000000000 * exp((-0.8) * tc[0]);
  Corr = mixture + sc[1] + (5) * sc[5] + (0.5) * sc[9] + sc[10] + sc[11];
  redP = Corr / k_f * 3.2e+15 * exp(-3.14 * tc[0] - (618.956498939805) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (0.32) * exp(-tc[1] * 0.0128205128205128) +
    0.68 * exp(-tc[1] * 0.00050125313283208) + exp(-5590 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[3] *= Corr * k_f;
  qr[3] *= Corr * k_f * exp(-(g_RT[2] - g_RT[12] + g_RT[13])) * (refC);
  // reaction 44:  CH2 + O2 <=> HCO + OH
  k_f = 10600000 * exp(-(754.82499870708) * invT);
  qf[54] *= k_f;
  qr[54] *= k_f * exp(-(g_RT[3] - g_RT[6] + g_RT[13] - g_RT[21]));
  // reaction 45:  CH2 + O2 => CO2 + 2 H
  k_f = 2640000 * exp(-(754.82499870708) * invT);
  qf[55] *= k_f;
  qr[55] *= k_f * exp(-(-2.000000 * g_RT[2] + g_RT[3] - g_RT[10] + g_RT[13])) *
            (refCinv);
  // reaction 46:  CH2 + O => CO + 2 H
  k_f = 50000000;
  qf[56] *= k_f;
  qr[56] *= k_f * exp(-(-2.000000 * g_RT[2] + g_RT[4] - g_RT[9] + g_RT[13])) *
            (refCinv);
  // reaction 47:  CH2 + H <=> CH + H2
  k_f = 283000000 * exp(-(369.864249366469) * invT);
  qf[57] *= k_f;
  qr[57] *= k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[13] - g_RT[16]));
  // reaction 48:  CH2 + OH <=> CH + H2O
  k_f = 11.3 * exp((2) * tc[0] - (1509.64999741416) * invT);
  qf[58] *= k_f;
  qr[58] *= k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[13] - g_RT[16]));
  // reaction 49:  CH + O2 <=> HCO + O
  k_f = 33000000;
  qf[59] *= k_f;
  qr[59] *= k_f * exp(-(g_RT[3] - g_RT[4] + g_RT[16] - g_RT[21]));
  // reaction 50:  CH + O <=> CO + H
  k_f = 57000000;
  qf[60] *= k_f;
  qr[60] *= k_f * exp(-(-g_RT[2] + g_RT[4] - g_RT[9] + g_RT[16]));
  // reaction 51:  CH + H <=> C + H2
  k_f = 110000000;
  qf[61] *= k_f;
  qr[61] *= k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[15] + g_RT[16]));
  // reaction 52:  CH + OH <=> H + HCO
  k_f = 30000000;
  qf[62] *= k_f;
  qr[62] *= k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[16] - g_RT[21]));
  // reaction 53:  CH + H2O <=> CH2O + H
  k_f = 17130000 * exp(-(-379.928582682563) * invT);
  qf[63] *= k_f;
  qr[63] *= k_f * exp(-(-g_RT[2] + g_RT[5] + g_RT[16] - g_RT[20]));
  // reaction 54:  CH + CO2 <=> CO + HCO
  k_f = 1700000 * exp(-(344.703416076233) * invT);
  qf[64] *= k_f;
  qr[64] *= k_f * exp(-(-g_RT[9] + g_RT[10] + g_RT[16] - g_RT[21]));
  // reaction 55:  C + OH <=> CO + H
  k_f = 50000000;
  qf[65] *= k_f;
  qr[65] *= k_f * exp(-(-g_RT[2] + g_RT[6] - g_RT[9] + g_RT[15]));
  // reaction 56:  C + O2 <=> CO + O
  k_f = 50000000;
  qf[66] *= k_f;
  qr[66] *= k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[9] + g_RT[15]));
  // reaction 57:  CH3 + O2 <=> CH3O + O
  k_f = 7546000 * exp(-(14251.0959755897) * invT);
  qf[67] *= k_f;
  qr[67] *= k_f * exp(-(g_RT[3] - g_RT[4] + g_RT[12] - g_RT[18]));
  // reaction 58:  CH3 + O2 <=> CH2O + OH
  k_f = 2.641e-06 * exp((3.283) * tc[0] - (4078.57107634725) * invT);
  qf[68] *= k_f;
  qr[68] *= k_f * exp(-(g_RT[3] - g_RT[6] + g_RT[12] - g_RT[20]));
  // reaction 59:  CH3 + O <=> CH2O + H
  k_f = 55400000 * exp((0.05) * tc[0] - (-68.4374665494419) * invT);
  qf[69] *= k_f;
  qr[69] *= k_f * exp(-(-g_RT[2] + g_RT[4] + g_RT[12] - g_RT[20]));
  // reaction 60:  CH3 + OH <=> CH2(S) + H2O
  k_f = 4.788e+17 * exp((-3.155) * tc[0] - (3524.02631063045) * invT);
  qf[70] *= k_f;
  qr[70] *= k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[12] - g_RT[14]));
  // reaction 61:  CH3 + OH <=> CH2O + H2
  k_f = 5374 * exp((0.287) * tc[0] - (140.900666425322) * invT);
  qf[71] *= k_f;
  qr[71] *= k_f * exp(-(-g_RT[1] + g_RT[6] + g_RT[12] - g_RT[20]));
  // reaction 62:  CH3 + OH <=> CH2OH + H
  k_f = 15250000 * exp((0.134) * tc[0] - (2838.64521180442) * invT);
  qf[72] *= k_f;
  qr[72] *= k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[12] - g_RT[19]));
  // reaction 63:  CH3 + OH <=> CH3O + H
  k_f = 1798 * exp((0.965) * tc[0] - (6068.79298960492) * invT);
  qf[73] *= k_f;
  qr[73] *= k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[12] - g_RT[18]));
  // reaction 64:  CH3 + OH <=> CH2 + H2O
  k_f = 0.04293 * exp((2.568) * tc[0] - (2011.75958655411) * invT);
  qf[74] *= k_f;
  qr[74] *= k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[12] - g_RT[13]));
  // reaction 65:  CH3 + HO2 <=> CH3O + OH
  k_f = 1000000 * exp((0.269) * tc[0] - (-345.961457740745) * invT);
  qf[75] *= k_f;
  qr[75] *= k_f * exp(-(-g_RT[6] + g_RT[8] + g_RT[12] - g_RT[18]));
  // reaction 66:  CH3OH (+M) <=> CH3 + OH (+M)
  k_f = 2.084e+18 * exp((-0.615) * tc[0] - (46567.9721835682) * invT);
  Corr = mixture;
  redP = Corr / k_f * 1.5e+37 * exp(-6.995 * tc[0] - (49311.3081588693) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (1.4748) * exp(-tc[1] * 2.81056773468241e-05) +
    -0.4748 * exp(-tc[1] * 0.000896057347670251) + exp(-9023 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[4] *= Corr * k_f;
  qr[4] *= Corr * k_f * exp(-(-g_RT[6] - g_RT[12] + g_RT[17])) * (refCinv);
  // reaction 67:  CH3OH (+M) <=> CH2(S) + H2O (+M)
  k_f = 3.121e+18 * exp((-1.017) * tc[0] - (46151.0068542825) * invT);
  Corr = mixture;
  redP =
    Corr / k_f * 1.43e+41 * exp(-8.227 * tc[0] - (50028.3415859744) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (-1.545) * exp(-tc[1] * 0.000303951367781155) +
    2.545 * exp(-tc[1] * 2.11327134404057e-05) + exp(-47110 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[5] *= Corr * k_f;
  qr[5] *= Corr * k_f * exp(-(-g_RT[5] - g_RT[14] + g_RT[17])) * (refCinv);
  // reaction 68:  CH3OH + H <=> CH3O + H2
  k_f = 0.199 * exp((2.56) * tc[0] - (5183.13165778861) * invT);
  qf[76] *= k_f;
  qr[76] *= k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[17] - g_RT[18]));
  // reaction 69:  CH3OH + H <=> CH2OH + H2
  k_f = 0.307 * exp((2.55) * tc[0] - (2737.49866197768) * invT);
  qf[77] *= k_f;
  qr[77] *= k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[17] - g_RT[19]));
  // reaction 70:  CH3OH + O <=> CH3O + OH
  k_f = 0.0388 * exp((2.5) * tc[0] - (1549.90733067854) * invT);
  qf[78] *= k_f;
  qr[78] *= k_f * exp(-(g_RT[4] - g_RT[6] + g_RT[17] - g_RT[18]));
  // reaction 71:  CH3OH + O <=> CH2OH + OH
  k_f = 0.388 * exp((2.5) * tc[0] - (1549.90733067854) * invT);
  qf[79] *= k_f;
  qr[79] *= k_f * exp(-(g_RT[4] - g_RT[6] + g_RT[17] - g_RT[19]));
  // reaction 72:  CH3OH + OH <=> CH3O + H2O
  k_f = 0.00015 * exp((3.03) * tc[0] - (-383.954316009001) * invT);
  qf[80] *= k_f;
  qr[80] *= k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[17] - g_RT[18]));
  // reaction 73:  CH3OH + OH <=> CH2OH + H2O
  k_f = 0.0308 * exp((2.65) * tc[0] - (-405.944884304667) * invT);
  qf[81] *= k_f;
  qr[81] *= k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[17] - g_RT[19]));
  // reaction 74:  CH3OH + O2 <=> CH3O + HO2
  k_f = 0.0358 * exp((2.27) * tc[0] - (21519.8091048059) * invT);
  qf[82] *= k_f;
  qr[82] *= k_f * exp(-(g_RT[3] - g_RT[8] + g_RT[17] - g_RT[18]));
  // reaction 75:  CH3OH + HO2 <=> CH3O + H2O2
  k_f = 1220000 * exp(-(10099.9107343668) * invT);
  qf[83] *= k_f;
  qr[83] *= k_f * exp(-(-g_RT[7] + g_RT[8] + g_RT[17] - g_RT[18]));
  // reaction 76:  CH3OH + HO2 <=> CH2OH + H2O2
  k_f = 32600000 * exp(-(9451.5160604774) * invT);
  qf[84] *= k_f;
  qr[84] *= k_f * exp(-(-g_RT[7] + g_RT[8] + g_RT[17] - g_RT[19]));
  // reaction 77:  CH3 + CH3OH <=> CH2OH + CH4
  k_f = 2.13e-07 * exp((3.953) * tc[0] - (3550.24389891888) * invT);
  qf[85] *= k_f;
  qr[85] *= k_f * exp(-(-g_RT[11] + g_RT[12] + g_RT[17] - g_RT[19]));
  // reaction 78:  CH3 + CH3OH <=> CH3O + CH4
  k_f = 0.00322 * exp((2.425) * tc[0] - (4317.34738427159) * invT);
  qf[86] *= k_f;
  qr[86] *= k_f * exp(-(-g_RT[11] + g_RT[12] + g_RT[17] - g_RT[18]));
  // reaction 79:  CH2OH + O2 <=> CH2O + HO2
  k_f = 1510000000 * exp((-1) * tc[0]);
  qf[87] *= k_f;
  qr[87] *= k_f * exp(-(g_RT[3] - g_RT[8] + g_RT[19] - g_RT[20]));
  // reaction 80:  CH2OH + O2 <=> CH2O + HO2
  k_f = 241000000 * exp(-(2524.63801234228) * invT);
  qf[88] *= k_f;
  qr[88] *= k_f * exp(-(g_RT[3] - g_RT[8] + g_RT[19] - g_RT[20]));
  // reaction 81:  CH2OH + H <=> CH2O + H2
  k_f = 6000000;
  qf[89] *= k_f;
  qr[89] *= k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[19] - g_RT[20]));
  // reaction 82:  CH2OH + HO2 <=> CH2O + H2O2
  k_f = 12000000;
  qf[90] *= k_f;
  qr[90] *= k_f * exp(-(-g_RT[7] + g_RT[8] + g_RT[19] - g_RT[20]));
  // reaction 83:  CH2OH + HCO <=> 2 CH2O
  k_f = 180000000;
  qf[91] *= k_f;
  qr[91] *= k_f * exp(-(g_RT[19] - 2.000000 * g_RT[20] + g_RT[21]));
  // reaction 84:  CH2OH + OH <=> CH2O + H2O
  k_f = 24000000;
  qf[92] *= k_f;
  qr[92] *= k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[19] - g_RT[20]));
  // reaction 85:  CH2OH + O <=> CH2O + OH
  k_f = 42000000;
  qf[93] *= k_f;
  qr[93] *= k_f * exp(-(g_RT[4] - g_RT[6] + g_RT[19] - g_RT[20]));
  // reaction 86:  CH3O + O2 <=> CH2O + HO2
  k_f = 4.38e-25 * exp((9.5) * tc[0] - (-2768.19487859176) * invT);
  qf[94] *= k_f;
  qr[94] *= k_f * exp(-(g_RT[3] - g_RT[8] + g_RT[18] - g_RT[20]));
  // reaction 87:  CH3O + H <=> CH2O + H2
  k_f = 20000000;
  qf[95] *= k_f;
  qr[95] *= k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[18] - g_RT[20]));
  // reaction 88:  CH3O + HO2 <=> CH2O + H2O2
  k_f = 301000;
  qf[96] *= k_f;
  qr[96] *= k_f * exp(-(-g_RT[7] + g_RT[8] + g_RT[18] - g_RT[20]));
  // reaction 89:  CH3 + CH3O <=> CH2O + CH4
  k_f = 12000000;
  qf[97] *= k_f;
  qr[97] *= k_f * exp(-(-g_RT[11] + g_RT[12] + g_RT[18] - g_RT[20]));
  // reaction 90:  2 CH3O <=> CH2O + CH3OH
  k_f = 60300000;
  qf[98] *= k_f;
  qr[98] *= k_f * exp(-(-g_RT[17] + 2.000000 * g_RT[18] - g_RT[20]));
  // reaction 91:  CH2O + O <=> HCO + OH
  k_f = 6260 * exp((1.15) * tc[0] - (1137.26966471867) * invT);
  qf[99] *= k_f;
  qr[99] *= k_f * exp(-(g_RT[4] - g_RT[6] + g_RT[20] - g_RT[21]));
  // reaction 92:  CH2O + H <=> H2 + HCO
  k_f = 57.4 * exp((1.9) * tc[0] - (1378.81366430493) * invT);
  qf[100] *= k_f;
  qr[100] *= k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[20] - g_RT[21]));
  // reaction 93:  CH2O + OH <=> H2O + HCO
  k_f = 78.2 * exp((1.63) * tc[0] - (-530.893582423979) * invT);
  qf[101] *= k_f;
  qr[101] *= k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[20] - g_RT[21]));
  // reaction 94:  CH2O + HO2 <=> H2O2 + HCO
  k_f = 0.0188 * exp((2.7) * tc[0] - (5797.05599007037) * invT);
  qf[102] *= k_f;
  qr[102] *= k_f * exp(-(-g_RT[7] + g_RT[8] + g_RT[20] - g_RT[21]));
  // reaction 95:  CH2O + CH3 <=> CH4 + HCO
  k_f = 3.83e-05 * exp((3.36) * tc[0] - (2169.87026294995) * invT);
  qf[103] *= k_f;
  qr[103] *= k_f * exp(-(-g_RT[11] + g_RT[12] + g_RT[20] - g_RT[21]));
  // reaction 96:  CH2O + CH3O <=> CH3OH + HCO
  k_f = 662000 * exp(-(1154.37903135603) * invT);
  qf[104] *= k_f;
  qr[104] *= k_f * exp(-(-g_RT[17] + g_RT[18] + g_RT[20] - g_RT[21]));
  // reaction 97:  HCO + M <=> CO + H + M
  k_f = 570000 * exp((0.66) * tc[0] - (7482.83182051618) * invT);
  Corr = mixture + sc[1] + (5) * sc[5] + (0.5) * sc[9] + sc[10] + sc[11];
  qf[11] *= Corr * k_f;
  qr[11] *= Corr * k_f * exp(-(-g_RT[2] - g_RT[9] + g_RT[21])) * (refCinv);
  // reaction 98:  HCO + O2 <=> CO + HO2
  k_f = 7580000 * exp(-(206.318832979935) * invT);
  qf[105] *= k_f;
  qr[105] *= k_f * exp(-(g_RT[3] - g_RT[8] - g_RT[9] + g_RT[21]));
  // reaction 99:  HCO + O <=> CO + OH
  k_f = 30200000;
  qf[106] *= k_f;
  qr[106] *= k_f * exp(-(g_RT[4] - g_RT[6] - g_RT[9] + g_RT[21]));
  // reaction 100:  H + HCO <=> CO + H2
  k_f = 73400000;
  qf[107] *= k_f;
  qr[107] *= k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[9] + g_RT[21]));
  // reaction 101:  HCO + OH <=> CO + H2O
  k_f = 30110000;
  qf[108] *= k_f;
  qr[108] *= k_f * exp(-(-g_RT[5] + g_RT[6] - g_RT[9] + g_RT[21]));
  // reaction 102:  CH3 + HCO <=> CH4 + CO
  k_f = 26500000;
  qf[109] *= k_f;
  qr[109] *= k_f * exp(-(-g_RT[9] - g_RT[11] + g_RT[12] + g_RT[21]));
  // reaction 103:  2 HCO <=> CH2O + CO
  k_f = 18000000;
  qf[110] *= k_f;
  qr[110] *= k_f * exp(-(-g_RT[9] - g_RT[20] + 2.000000 * g_RT[21]));
  // reaction 104:  HCO + O <=> CO2 + H
  k_f = 30000000;
  qf[111] *= k_f;
  qr[111] *= k_f * exp(-(-g_RT[2] + g_RT[4] - g_RT[10] + g_RT[21]));
  // reaction 105:  HCO + HO2 => CO2 + H + OH
  k_f = 30000000;
  qf[112] *= k_f;
  qr[112] *= k_f * exp(-(-g_RT[2] - g_RT[6] + g_RT[8] - g_RT[10] + g_RT[21])) *
             (refCinv);
  // reaction 106:  CH2O + H (+M) <=> CH2OH (+M)
  k_f = 540000 * exp((0.454) * tc[0] - (1811.57999689699) * invT);
  Corr = mixture + sc[1] + (5) * sc[5] + (0.5) * sc[9] + sc[10] + sc[11];
  redP = Corr / k_f * 1.27e+20 * exp(-4.82 * tc[0] - (3286.00482770482) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (0.2813) * exp(-tc[1] * 0.00970873786407767) +
    0.7187 * exp(-tc[1] * 0.000774593338497289) + exp(-4160 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[6] *= Corr * k_f;
  qr[6] *= Corr * k_f * exp(-(g_RT[2] - g_RT[19] + g_RT[20])) * (refC);
  // reaction 107:  CH3O (+M) <=> CH2O + H (+M)
  k_f = 68000000000000 * exp(-(13169.1801441095) * invT);
  Corr = mixture + sc[1] + (5) * sc[5] + (0.5) * sc[9] + sc[10] + sc[11];
  redP = Corr / k_f * 1.87e+19 * exp(-3 * tc[0] - (12231.6874957153) * invT);
  F = redP / (1.0 + redP);
  logPred = log10(redP);
  logFcent = log10(
    (0.1) * exp(-tc[1] * 0.0004) + 0.9 * exp(-tc[1] * 0.000769230769230769) +
    exp(-1e+99 * invT));
  troe_c = -0.4 - 0.67 * logFcent;
  troe_n = 0.75 - 1.27 * logFcent;
  troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
  F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
  Corr = F * F_troe;
  qf[7] *= Corr * k_f;
  qr[7] *= Corr * k_f * exp(-(-g_RT[2] + g_RT[18] - g_RT[20])) * (refCinv);
  // reaction 108:  O + OH + M <=> HO2 + M
  k_f = 10000;
  Corr = mixture;
  qf[12] *= Corr * k_f;
  qr[12] *= Corr * k_f * exp(-(g_RT[4] + g_RT[6] - g_RT[8])) * (refC);
  // reaction 109:  CH + OH <=> C + H2O
  k_f = 40 * exp((2) * tc[0] - (1509.64999741416) * invT);
  qf[113] *= k_f;
  qr[113] *= k_f * exp(-(-g_RT[5] + g_RT[6] - g_RT[15] + g_RT[16]));
  // reaction 110:  CH2 + HCO <=> CH3 + CO
  k_f = 20000000;
  qf[114] *= k_f;
  qr[114] *= k_f * exp(-(-g_RT[9] - g_RT[12] + g_RT[13] + g_RT[21]));
  // reaction 111:  CH2 + O <=> CO + H2
  k_f = 50000000;
  qf[115] *= k_f;
  qr[115] *= k_f * exp(-(-g_RT[1] + g_RT[4] - g_RT[9] + g_RT[13]));
  // reaction 112:  CH2 + OH <=> CH2O + H
  k_f = 30000000;
  qf[116] *= k_f;
  qr[116] *= k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[13] - g_RT[20]));
  // reaction 113:  CH2 + CO2 <=> CH2O + CO
  k_f = 110000 * exp(-(503.21666580472) * invT);
  qf[117] *= k_f;
  qr[117] *= k_f * exp(-(-g_RT[9] + g_RT[10] + g_RT[13] - g_RT[20]));
  // reaction 114:  CH2(S) + O <=> CO + 2 H
  k_f = 30000000;
  qf[118] *= k_f;
  qr[118] *= k_f * exp(-(-2.000000 * g_RT[2] + g_RT[4] - g_RT[9] + g_RT[14])) *
             (refCinv);
  // reaction 115:  CH2(S) + CH4 <=> 2 CH3
  k_f = 43000000;
  qf[119] *= k_f;
  qr[119] *= k_f * exp(-(g_RT[11] - 2.000000 * g_RT[12] + g_RT[14]));
  // reaction 116:  CH3OH + H => CH3 + H2O
  k_f = 650000 * exp(-(2667.04832876501) * invT);
  qf[120] *= k_f;
  qr[120] *= k_f * exp(-(g_RT[2] - g_RT[5] - g_RT[12] + g_RT[17]));
  // reaction 117:  CH2O + O => CO2 + 2 H
  k_f = 20000000 * exp(-(2516.0833290236) * invT);
  qf[121] *= k_f;
  qr[121] *= k_f * exp(-(-2.000000 * g_RT[2] + g_RT[4] - g_RT[10] + g_RT[20])) *
             (refCinv);
  // reaction 118:  CH2O + OH => CO2 + H + H2
  k_f = 100000;
  qf[122] *= k_f;
  qr[122] *= k_f * exp(-(-g_RT[1] - g_RT[2] + g_RT[6] - g_RT[10] + g_RT[20])) *
             (refCinv);
  // reaction 119:  CH3O + CO <=> CH3 + CO2
  k_f = 500000 * exp(-(3270.90832773068) * invT);
  qf[123] *= k_f;
  qr[123] *= k_f * exp(-(g_RT[9] - g_RT[10] - g_RT[12] + g_RT[18]));
  // reaction 120:  CH3 + O + M => CH3O + M
  k_f = 50000;
  Corr = mixture;
  qf[13] *= Corr * k_f;
  qr[13] *= Corr * k_f * exp(-(g_RT[4] + g_RT[12] - g_RT[18])) * (refC);
  // reaction 121:  HCO + HO2 <=> CO + H2O2
  k_f = 400000;
  qf[124] *= k_f;
  qr[124] *= k_f * exp(-(-g_RT[7] + g_RT[8] - g_RT[9] + g_RT[21]));
  // reaction 122:  CH3O + OH => CH2O + H2O
  k_f = 15000000;
  qf[125] *= k_f;
  qr[125] *= k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[18] - g_RT[20]));
  // reaction 123:  CH2OH + CH3 => CH2O + CH4
  k_f = 15000000;
  qf[126] *= k_f;
  qr[126] *= k_f * exp(-(-g_RT[11] + g_RT[12] + g_RT[19] - g_RT[20]));
  // reaction 124:  CH + O <=> CHO+ + E
  k_f = 1746000000000 * exp((-2.19) * tc[0] - (164.999209334044) * invT);
  qf[127] *= k_f;
  qr[127] *= k_f * exp(-(g_RT[4] + g_RT[16] - g_RT[22] - g_RT[25]));
  // reaction 125:  CHO+ + H2O <=> CO + H3O+
  k_f = 26080000000 * exp((-0.5) * tc[0]);
  qf[128] *= k_f;
  qr[128] *= k_f * exp(-(g_RT[5] - g_RT[9] + g_RT[22] - g_RT[23]));
  // reaction 126:  E + H3O+ <=> H + H2O
  k_f = 739500000000 * exp((-0.5) * tc[0]);
  qf[129] *= k_f;
  qr[129] *= k_f * exp(-(-g_RT[2] - g_RT[5] + g_RT[23] + g_RT[25]));
  // reaction 127:  E + H3O+ <=> 2 H + OH
  k_f = 3181000000000 * exp((-0.5) * tc[0]);
  qf[130] *= k_f;
  qr[130] *= k_f * exp(-(-2.000000 * g_RT[2] - g_RT[6] + g_RT[23] + g_RT[25])) *
             (refCinv);
  // reaction 128:  E + H3O+ <=> H2 + OH
  k_f = 560100000000 * exp((-0.5) * tc[0]);
  qf[131] *= k_f;
  qr[131] *= k_f * exp(-(-g_RT[1] - g_RT[6] + g_RT[23] + g_RT[25]));
  // reaction 129:  H2 + O2- <=> E + H2O2
  k_f = 602000000;
  qf[132] *= k_f;
  qr[132] *= k_f * exp(-(g_RT[1] - g_RT[7] + g_RT[24] - g_RT[25]));
  // reaction 130:  H + O2- <=> E + HO2
  k_f = 723000000;
  qf[133] *= k_f;
  qr[133] *= k_f * exp(-(g_RT[2] - g_RT[8] + g_RT[24] - g_RT[25]));
  // reaction 131:  E + O2 + O2 <=> O2- + O2
  k_f = 1520000000 * exp((-1) * tc[0] - (600.156324305341) * invT);
  qf[14] *= k_f;
  qr[14] *= k_f * exp(-(g_RT[3] - g_RT[24] + g_RT[25])) * (refC);
  // reaction 132:  E + O2 + H2O <=> O2- + H2O
  k_f = 5080000;
  qf[15] *= k_f;
  qr[15] *= k_f * exp(-(g_RT[3] - g_RT[24] + g_RT[25])) * (refC);
  // reaction 133:  E + O2 + N2 <=> O2- + N2
  k_f = 3590000000 * exp((-2) * tc[0] - (69.7558942138502) * invT);
  qf[16] *= k_f;
  qr[16] *= k_f * exp(-(g_RT[3] - g_RT[24] + g_RT[25])) * (refC);
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

  for (int i = 0; i < 26; ++i) {
    wdot[i] = 0.0;
  }

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 26; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[26];
  gibbs(g_RT, tc);

  {
    // reaction 6:  H2O2 (+M) <=> 2 OH (+M)
    const amrex::Real k_f =
      2000000000000 * exp((0.9) * tc[0] - (24531.3092413143) * invT);
    amrex::Real Corr = mixture + (0.5) * sc[0] + (2.7) * sc[1] + (0.2) * sc[3] +
                       (6.65) * sc[5] + (6.7) * sc[7] + (1.8) * sc[9] +
                       (0.6) * sc[10];
    const amrex::Real redP =
      Corr / k_f * 2.49e+18 * exp(-2.3 * tc[0] - 24531.3092413143 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.57 * exp(-tc[1] * 1e+30) + 0.43 * exp(-tc[1] * 1e-30) +
      exp(-1e-16 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[7]);
    const amrex::Real qr = Corr * k_f * exp(-(-2.000000 * g_RT[6] + g_RT[7])) *
                           (refCinv) * ((sc[6] * sc[6]));
    const amrex::Real qdot = qf - qr;
    wdot[6] += 2.000000 * qdot;
    wdot[7] -= qdot;
  }

  {
    // reaction 19:  H + O2 (+M) <=> HO2 (+M)
    const amrex::Real k_f = 4650000 * exp((0.44) * tc[0]);
    amrex::Real Corr = mixture + (0.3) * sc[1] + (9) * sc[5] + (0.9) * sc[9] +
                       (2.8) * sc[10] + sc[11];
    const amrex::Real redP = Corr / k_f * 17400000 * exp(-1.23 * tc[0]);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.33 * exp(-tc[1] * 1e+30) + 0.67 * exp(-tc[1] * 1e-30) +
      exp(-1e+30 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[2] * sc[3]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[2] + g_RT[3] - g_RT[8])) * (refC) * (sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[8] += qdot;
  }

  {
    // reaction 24:  CH3 + H (+M) <=> CH4 (+M)
    const amrex::Real k_f =
      12700000000 * exp((-0.63) * tc[0] - (192.731983003208) * invT);
    amrex::Real Corr =
      mixture + sc[1] + (5) * sc[5] + (0.5) * sc[9] + sc[10] + sc[11];
    const amrex::Real redP =
      Corr / k_f * 2.48e+21 * exp(-4.76 * tc[0] - 1227.84866456352 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.217 * exp(-tc[1] * 0.0135135135135135) +
      0.783 * exp(-tc[1] * 0.000340020401224073) + exp(-6964 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[2] * sc[12]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[2] - g_RT[11] + g_RT[12])) * (refC) * (sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;
  }

  {
    // reaction 43:  CH2 + H (+M) <=> CH3 (+M)
    const amrex::Real k_f = 25000000000 * exp((-0.8) * tc[0]);
    amrex::Real Corr =
      mixture + sc[1] + (5) * sc[5] + (0.5) * sc[9] + sc[10] + sc[11];
    const amrex::Real redP =
      Corr / k_f * 3.2e+15 * exp(-3.14 * tc[0] - 618.956498939805 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.32 * exp(-tc[1] * 0.0128205128205128) +
      0.68 * exp(-tc[1] * 0.00050125313283208) + exp(-5590 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[2] * sc[13]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[2] - g_RT[12] + g_RT[13])) * (refC) * (sc[12]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[12] += qdot;
    wdot[13] -= qdot;
  }

  {
    // reaction 66:  CH3OH (+M) <=> CH3 + OH (+M)
    const amrex::Real k_f =
      2.084e+18 * exp((-0.615) * tc[0] - (46567.9721835682) * invT);
    amrex::Real Corr = mixture;
    const amrex::Real redP =
      Corr / k_f * 1.5e+37 * exp(-6.995 * tc[0] - 49311.3081588693 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      1.4748 * exp(-tc[1] * 2.81056773468241e-05) +
      -0.4748 * exp(-tc[1] * 0.000896057347670251) + exp(-9023 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[17]);
    const amrex::Real qr = Corr * k_f * exp(-(-g_RT[6] - g_RT[12] + g_RT[17])) *
                           (refCinv) * (sc[6] * sc[12]);
    const amrex::Real qdot = qf - qr;
    wdot[6] += qdot;
    wdot[12] += qdot;
    wdot[17] -= qdot;
  }

  {
    // reaction 67:  CH3OH (+M) <=> CH2(S) + H2O (+M)
    const amrex::Real k_f =
      3.121e+18 * exp((-1.017) * tc[0] - (46151.0068542825) * invT);
    amrex::Real Corr = mixture;
    const amrex::Real redP =
      Corr / k_f * 1.43e+41 * exp(-8.227 * tc[0] - 50028.3415859744 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      -1.545 * exp(-tc[1] * 0.000303951367781155) +
      2.545 * exp(-tc[1] * 2.11327134404057e-05) + exp(-47110 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[17]);
    const amrex::Real qr = Corr * k_f * exp(-(-g_RT[5] - g_RT[14] + g_RT[17])) *
                           (refCinv) * (sc[5] * sc[14]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[14] += qdot;
    wdot[17] -= qdot;
  }

  {
    // reaction 106:  CH2O + H (+M) <=> CH2OH (+M)
    const amrex::Real k_f =
      540000 * exp((0.454) * tc[0] - (1811.57999689699) * invT);
    amrex::Real Corr =
      mixture + sc[1] + (5) * sc[5] + (0.5) * sc[9] + sc[10] + sc[11];
    const amrex::Real redP =
      Corr / k_f * 1.27e+20 * exp(-4.82 * tc[0] - 3286.00482770482 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.2813 * exp(-tc[1] * 0.00970873786407767) +
      0.7187 * exp(-tc[1] * 0.000774593338497289) + exp(-4160 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[2] * sc[20]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[2] - g_RT[19] + g_RT[20])) * (refC) * (sc[19]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[19] += qdot;
    wdot[20] -= qdot;
  }

  {
    // reaction 107:  CH3O (+M) <=> CH2O + H (+M)
    const amrex::Real k_f = 68000000000000 * exp(-(13169.1801441095) * invT);
    amrex::Real Corr =
      mixture + sc[1] + (5) * sc[5] + (0.5) * sc[9] + sc[10] + sc[11];
    const amrex::Real redP =
      Corr / k_f * 1.87e+19 * exp(-3 * tc[0] - 12231.6874957153 * invT);
    const amrex::Real F = redP / (1.0 + redP);
    const amrex::Real logPred = log10(redP);
    const amrex::Real logFcent = log10(
      0.1 * exp(-tc[1] * 0.0004) + 0.9 * exp(-tc[1] * 0.000769230769230769) +
      exp(-1e+99 * invT));
    const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
    const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
    const amrex::Real troe =
      (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
    const amrex::Real F_troe = exp(M_LN10 * logFcent / (1.0 + troe * troe));
    Corr = F * F_troe;
    const amrex::Real qf = Corr * k_f * (sc[18]);
    const amrex::Real qr = Corr * k_f * exp(-(-g_RT[2] + g_RT[18] - g_RT[20])) *
                           (refCinv) * (sc[2] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[18] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 20:  CO + O (+M) <=> CO2 (+M)
    const amrex::Real k_f = 13620 * exp(-(1199.66853127845) * invT);
    amrex::Real Corr =
      mixture + sc[1] + (11) * sc[5] + (0.75) * sc[9] + (2.6) * sc[10];
    const amrex::Real redP =
      Corr / k_f * 1170000000000 * exp(-2.79 * tc[0] - 2108.98104638758 * invT);
    Corr = redP / (1.0 + redP);
    const amrex::Real qf = Corr * k_f * (sc[4] * sc[9]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[4] + g_RT[9] - g_RT[10])) * (refC) * (sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 3:  H + OH + M <=> H2O + M
    const amrex::Real k_f = 35000000000 * exp((-2) * tc[0]);
    const amrex::Real Corr =
      mixture + (-0.27) * sc[1] + (2.65) * sc[5] + sc[11];
    const amrex::Real qf = Corr * k_f * (sc[2] * sc[6]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[2] - g_RT[5] + g_RT[6])) * (refC) * (sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;
  }

  {
    // reaction 5:  H + O + M <=> OH + M
    const amrex::Real k_f = 4714000 * exp((-1) * tc[0]);
    const amrex::Real Corr =
      mixture + (1.5) * sc[1] + (11) * sc[5] + (0.5) * sc[9] + sc[10] + sc[11];
    const amrex::Real qf = Corr * k_f * (sc[2] * sc[4]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[2] + g_RT[4] - g_RT[6])) * (refC) * (sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[4] -= qdot;
    wdot[6] += qdot;
  }

  {
    // reaction 97:  HCO + M <=> CO + H + M
    const amrex::Real k_f =
      570000 * exp((0.66) * tc[0] - (7482.83182051618) * invT);
    const amrex::Real Corr =
      mixture + sc[1] + (5) * sc[5] + (0.5) * sc[9] + sc[10] + sc[11];
    const amrex::Real qf = Corr * k_f * (sc[21]);
    const amrex::Real qr = Corr * k_f * exp(-(-g_RT[2] - g_RT[9] + g_RT[21])) *
                           (refCinv) * (sc[2] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[9] += qdot;
    wdot[21] -= qdot;
  }

  {
    // reaction 108:  O + OH + M <=> HO2 + M
    const amrex::Real k_f = 10000;
    const amrex::Real Corr = mixture;
    const amrex::Real qf = Corr * k_f * (sc[4] * sc[6]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[4] + g_RT[6] - g_RT[8])) * (refC) * (sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[6] -= qdot;
    wdot[8] += qdot;
  }

  {
    // reaction 120:  CH3 + O + M => CH3O + M
    const amrex::Real k_f = 50000;
    const amrex::Real Corr = mixture;
    const amrex::Real qf = Corr * k_f * (sc[4] * sc[12]);
    const amrex::Real qr =
      Corr * k_f * exp(-(g_RT[4] + g_RT[12] - g_RT[18])) * (refC) * (0.0);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[12] -= qdot;
    wdot[18] += qdot;
  }

  {
    // reaction 131:  E + O2 + O2 <=> O2- + O2
    const amrex::Real k_f =
      1520000000 * exp((-1) * tc[0] - (600.156324305341) * invT);
    const amrex::Real qf = k_f * ((sc[3] * sc[3]) * sc[25]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[24] + g_RT[25])) * (refC) * (sc[3] * sc[24]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[24] += qdot;
    wdot[25] -= qdot;
  }

  {
    // reaction 132:  E + O2 + H2O <=> O2- + H2O
    const amrex::Real k_f = 5080000;
    const amrex::Real qf = k_f * (sc[3] * sc[5] * sc[25]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[24] + g_RT[25])) * (refC) * (sc[5] * sc[24]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[24] += qdot;
    wdot[25] -= qdot;
  }

  {
    // reaction 133:  E + O2 + N2 <=> O2- + N2
    const amrex::Real k_f =
      3590000000 * exp((-2) * tc[0] - (69.7558942138502) * invT);
    const amrex::Real qf = k_f * (sc[0] * sc[3] * sc[25]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[24] + g_RT[25])) * (refC) * (sc[0] * sc[24]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[24] += qdot;
    wdot[25] -= qdot;
  }

  {
    // reaction 0:  H2 + O <=> H + OH
    const amrex::Real k_f =
      0.0508 * exp((2.67) * tc[0] - (3166.2392612433) * invT);
    const amrex::Real qf = k_f * (sc[1] * sc[4]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[2] + g_RT[4] - g_RT[6])) * (sc[2] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[6] += qdot;
  }

  {
    // reaction 1:  H2 + OH <=> H + H2O
    const amrex::Real k_f = 43800000 * exp(-(3517.48449397499) * invT);
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
    // reaction 2:  H + O2 <=> O + OH
    const amrex::Real k_f = 114000000 * exp(-(7692.16995349094) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[3]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] + g_RT[3] - g_RT[4] - g_RT[6])) * (sc[4] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
  }

  {
    // reaction 4:  H2O + O <=> 2 OH
    const amrex::Real k_f =
      67 * exp((1.704) * tc[0] - (7541.60752708217) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[5]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[4] + g_RT[5] - 2.000000 * g_RT[6])) * ((sc[6] * sc[6]));
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += 2.000000 * qdot;
  }

  {
    // reaction 7:  H + H2O2 <=> H2O + OH
    const amrex::Real k_f = 24100000 * exp(-(1997.77016324474) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - g_RT[5] - g_RT[6] + g_RT[7])) * (sc[5] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;
  }

  {
    // reaction 8:  H + H2O2 <=> H2 + HO2
    const amrex::Real k_f =
      21500 * exp((1) * tc[0] - (3019.29999482832) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[8])) * (sc[1] * sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;
  }

  {
    // reaction 9:  H2O2 + O <=> HO2 + OH
    const amrex::Real k_f = 9.55 * exp((2) * tc[0] - (1997.77016324474) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[4] - g_RT[6] + g_RT[7] - g_RT[8])) * (sc[6] * sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;
  }

  {
    // reaction 10:  H2O2 + OH <=> H2O + HO2
    const amrex::Real k_f = 1740000 * exp(-(160.022899725901) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[7] - g_RT[8])) * (sc[5] * sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;
  }

  {
    // reaction 11:  H2O2 + OH <=> H2O + HO2
    const amrex::Real k_f = 75900000 * exp(-(3657.88194373451) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[7]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[7] - g_RT[8])) * (sc[5] * sc[8]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;
  }

  {
    // reaction 12:  H + HO2 <=> 2 OH
    const amrex::Real k_f = 70790000 * exp(-(148.448916412392) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - 2.000000 * g_RT[6] + g_RT[8])) * ((sc[6] * sc[6]));
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[6] += 2.000000 * qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 13:  H + HO2 <=> H2 + O2
    const amrex::Real k_f =
      11402 * exp((1.083) * tc[0] - (278.671325189338) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[8])) * (sc[1] * sc[3]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 14:  HO2 + O <=> O2 + OH
    const amrex::Real k_f = 32500000;
    const amrex::Real qf = k_f * (sc[4] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[3] + g_RT[4] - g_RT[6] + g_RT[8])) * (sc[3] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 15:  HO2 + OH <=> H2O + O2
    const amrex::Real k_f = 7000000 * exp(-(-549.995687057926) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[3] - g_RT[5] + g_RT[6] + g_RT[8])) * (sc[3] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[3] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 16:  HO2 + OH <=> H2O + O2
    const amrex::Real k_f = 450000000 * exp(-(5499.95687057926) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[8]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[3] - g_RT[5] + g_RT[6] + g_RT[8])) * (sc[3] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[3] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[8] -= qdot;
  }

  {
    // reaction 17:  2 HO2 <=> H2O2 + O2
    const amrex::Real k_f = 100000000 * exp(-(5555.95482115001) * invT);
    const amrex::Real qf = k_f * ((sc[8] * sc[8]));
    const amrex::Real qr =
      k_f * exp(-(-g_RT[3] - g_RT[7] + 2.000000 * g_RT[8])) * (sc[3] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[3] += qdot;
    wdot[7] += qdot;
    wdot[8] -= 2.000000 * qdot;
  }

  {
    // reaction 18:  2 HO2 <=> H2O2 + O2
    const amrex::Real k_f = 190000 * exp(-(-708.992024785586) * invT);
    const amrex::Real qf = k_f * ((sc[8] * sc[8]));
    const amrex::Real qr =
      k_f * exp(-(-g_RT[3] - g_RT[7] + 2.000000 * g_RT[8])) * (sc[3] * sc[7]);
    const amrex::Real qdot = qf - qr;
    wdot[3] += qdot;
    wdot[7] += qdot;
    wdot[8] -= 2.000000 * qdot;
  }

  {
    // reaction 21:  CO + OH <=> CO2 + H
    const amrex::Real k_f =
      0.07015 * exp((2.053) * tc[0] - (-178.994168026739) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[9]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[9] - g_RT[10])) * (sc[2] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[6] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 22:  CO + OH <=> CO2 + H
    const amrex::Real k_f =
      5757000 * exp((-0.664) * tc[0] - (166.967289714006) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[9]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[9] - g_RT[10])) * (sc[2] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[6] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 23:  CO + HO2 <=> CO2 + OH
    const amrex::Real k_f =
      0.157 * exp((2.18) * tc[0] - (9027.70698453667) * invT);
    const amrex::Real qf = k_f * (sc[8] * sc[9]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[6] + g_RT[8] + g_RT[9] - g_RT[10])) * (sc[6] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;
  }

  {
    // reaction 25:  CH4 + H <=> CH3 + H2
    const amrex::Real k_f =
      0.614 * exp((2.5) * tc[0] - (4824.33817506985) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[11]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[11] - g_RT[12])) * (sc[1] * sc[12]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;
  }

  {
    // reaction 26:  CH4 + O <=> CH3 + OH
    const amrex::Real k_f =
      1020 * exp((1.5) * tc[0] - (4327.66332592059) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[11]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[4] - g_RT[6] + g_RT[11] - g_RT[12])) * (sc[6] * sc[12]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;
  }

  {
    // reaction 27:  CH4 + OH <=> CH3 + H2O
    const amrex::Real k_f =
      0.0583 * exp((2.6) * tc[0] - (1102.04449811234) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[11]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[11] - g_RT[12])) * (sc[5] * sc[12]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;
  }

  {
    // reaction 28:  CH4 + HO2 <=> CH3 + H2O2
    const amrex::Real k_f =
      1.13e-05 * exp((3.74) * tc[0] - (10572.5821485572) * invT);
    const amrex::Real qf = k_f * (sc[8] * sc[11]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[7] + g_RT[8] + g_RT[11] - g_RT[12])) * (sc[7] * sc[12]);
    const amrex::Real qdot = qf - qr;
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;
  }

  {
    // reaction 29:  CH3 + HO2 <=> CH4 + O2
    const amrex::Real k_f =
      0.116 * exp((2.23) * tc[0] - (-1520.72076406186) * invT);
    const amrex::Real qf = k_f * (sc[8] * sc[12]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[3] + g_RT[8] - g_RT[11] + g_RT[12])) * (sc[3] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[3] += qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;
  }

  {
    // reaction 30:  CH2 + CH4 <=> 2 CH3
    const amrex::Real k_f = 2.46 * exp((2) * tc[0] - (4161.60182620503) * invT);
    const amrex::Real qf = k_f * (sc[11] * sc[13]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[11] - 2.000000 * g_RT[12] + g_RT[13])) *
                           ((sc[12] * sc[12]));
    const amrex::Real qdot = qf - qr;
    wdot[11] -= qdot;
    wdot[12] += 2.000000 * qdot;
    wdot[13] -= qdot;
  }

  {
    // reaction 31:  CH2(S) + N2 <=> CH2 + N2
    const amrex::Real k_f = 15000000 * exp(-(301.929999482832) * invT);
    const amrex::Real qf = k_f * (sc[0] * sc[14]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] - g_RT[0] - g_RT[13] + g_RT[14])) * (sc[0] * sc[13]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[0] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;
  }

  {
    // reaction 32:  CH2(S) + H2O <=> CH2 + H2O
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[5] * sc[14]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[5] - g_RT[5] - g_RT[13] + g_RT[14])) * (sc[5] * sc[13]);
    const amrex::Real qdot = qf - qr;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;
  }

  {
    // reaction 33:  CH2(S) + CO <=> CH2 + CO
    const amrex::Real k_f = 9000000;
    const amrex::Real qf = k_f * (sc[9] * sc[14]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[9] - g_RT[9] - g_RT[13] + g_RT[14])) * (sc[9] * sc[13]);
    const amrex::Real qdot = qf - qr;
    wdot[9] -= qdot;
    wdot[9] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;
  }

  {
    // reaction 34:  CH2(S) + CO2 <=> CH2 + CO2
    const amrex::Real k_f = 7000000;
    const amrex::Real qf = k_f * (sc[10] * sc[14]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[10] - g_RT[10] - g_RT[13] + g_RT[14])) *
                           (sc[10] * sc[13]);
    const amrex::Real qdot = qf - qr;
    wdot[10] -= qdot;
    wdot[10] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;
  }

  {
    // reaction 35:  CH2(S) + O2 => CO + H + OH
    const amrex::Real k_f = 28000000;
    const amrex::Real qf = k_f * (sc[3] * sc[14]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[9] += qdot;
    wdot[14] -= qdot;
  }

  {
    // reaction 36:  CH2(S) + O2 <=> CO + H2O
    const amrex::Real k_f = 12000000;
    const amrex::Real qf = k_f * (sc[3] * sc[14]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[5] - g_RT[9] + g_RT[14])) * (sc[5] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
    wdot[14] -= qdot;
  }

  {
    // reaction 37:  CH2(S) + O <=> CO + H2
    const amrex::Real k_f = 15000000;
    const amrex::Real qf = k_f * (sc[4] * sc[14]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[4] - g_RT[9] + g_RT[14])) * (sc[1] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[9] += qdot;
    wdot[14] -= qdot;
  }

  {
    // reaction 38:  CH2(S) + O <=> H + HCO
    const amrex::Real k_f = 15000000;
    const amrex::Real qf = k_f * (sc[4] * sc[14]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[4] + g_RT[14] - g_RT[21])) * (sc[2] * sc[21]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[14] -= qdot;
    wdot[21] += qdot;
  }

  {
    // reaction 39:  CH2(S) + H2 <=> CH3 + H
    const amrex::Real k_f = 70000000;
    const amrex::Real qf = k_f * (sc[1] * sc[14]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[2] - g_RT[12] + g_RT[14])) * (sc[2] * sc[12]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[12] += qdot;
    wdot[14] -= qdot;
  }

  {
    // reaction 40:  CH2(S) + H <=> CH + H2
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[2] * sc[14]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[14] - g_RT[16])) * (sc[1] * sc[16]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[14] -= qdot;
    wdot[16] += qdot;
  }

  {
    // reaction 41:  CH2(S) + OH <=> CH2O + H
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[6] * sc[14]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[14] - g_RT[20])) * (sc[2] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[6] -= qdot;
    wdot[14] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 42:  CH2(S) + CO2 <=> CH2O + CO
    const amrex::Real k_f = 14000000;
    const amrex::Real qf = k_f * (sc[10] * sc[14]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[9] + g_RT[10] + g_RT[14] - g_RT[20])) *
                           (sc[9] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[14] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 44:  CH2 + O2 <=> HCO + OH
    const amrex::Real k_f = 10600000 * exp(-(754.82499870708) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[13]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[6] + g_RT[13] - g_RT[21])) * (sc[6] * sc[21]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[13] -= qdot;
    wdot[21] += qdot;
  }

  {
    // reaction 45:  CH2 + O2 => CO2 + 2 H
    const amrex::Real k_f = 2640000 * exp(-(754.82499870708) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[13]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[2] += 2.000000 * qdot;
    wdot[3] -= qdot;
    wdot[10] += qdot;
    wdot[13] -= qdot;
  }

  {
    // reaction 46:  CH2 + O => CO + 2 H
    const amrex::Real k_f = 50000000;
    const amrex::Real qf = k_f * (sc[4] * sc[13]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[2] += 2.000000 * qdot;
    wdot[4] -= qdot;
    wdot[9] += qdot;
    wdot[13] -= qdot;
  }

  {
    // reaction 47:  CH2 + H <=> CH + H2
    const amrex::Real k_f = 283000000 * exp(-(369.864249366469) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[13]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[13] - g_RT[16])) * (sc[1] * sc[16]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[13] -= qdot;
    wdot[16] += qdot;
  }

  {
    // reaction 48:  CH2 + OH <=> CH + H2O
    const amrex::Real k_f = 11.3 * exp((2) * tc[0] - (1509.64999741416) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[13]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[13] - g_RT[16])) * (sc[5] * sc[16]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[13] -= qdot;
    wdot[16] += qdot;
  }

  {
    // reaction 49:  CH + O2 <=> HCO + O
    const amrex::Real k_f = 33000000;
    const amrex::Real qf = k_f * (sc[3] * sc[16]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[4] + g_RT[16] - g_RT[21])) * (sc[4] * sc[21]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[16] -= qdot;
    wdot[21] += qdot;
  }

  {
    // reaction 50:  CH + O <=> CO + H
    const amrex::Real k_f = 57000000;
    const amrex::Real qf = k_f * (sc[4] * sc[16]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[4] - g_RT[9] + g_RT[16])) * (sc[2] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[9] += qdot;
    wdot[16] -= qdot;
  }

  {
    // reaction 51:  CH + H <=> C + H2
    const amrex::Real k_f = 110000000;
    const amrex::Real qf = k_f * (sc[2] * sc[16]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[15] + g_RT[16])) * (sc[1] * sc[15]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;
  }

  {
    // reaction 52:  CH + OH <=> H + HCO
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[6] * sc[16]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[16] - g_RT[21])) * (sc[2] * sc[21]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[6] -= qdot;
    wdot[16] -= qdot;
    wdot[21] += qdot;
  }

  {
    // reaction 53:  CH + H2O <=> CH2O + H
    const amrex::Real k_f = 17130000 * exp(-(-379.928582682563) * invT);
    const amrex::Real qf = k_f * (sc[5] * sc[16]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[5] + g_RT[16] - g_RT[20])) * (sc[2] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[16] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 54:  CH + CO2 <=> CO + HCO
    const amrex::Real k_f = 1700000 * exp(-(344.703416076233) * invT);
    const amrex::Real qf = k_f * (sc[10] * sc[16]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[9] + g_RT[10] + g_RT[16] - g_RT[21])) *
                           (sc[9] * sc[21]);
    const amrex::Real qdot = qf - qr;
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[16] -= qdot;
    wdot[21] += qdot;
  }

  {
    // reaction 55:  C + OH <=> CO + H
    const amrex::Real k_f = 50000000;
    const amrex::Real qf = k_f * (sc[6] * sc[15]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[6] - g_RT[9] + g_RT[15])) * (sc[2] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[6] -= qdot;
    wdot[9] += qdot;
    wdot[15] -= qdot;
  }

  {
    // reaction 56:  C + O2 <=> CO + O
    const amrex::Real k_f = 50000000;
    const amrex::Real qf = k_f * (sc[3] * sc[15]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[9] + g_RT[15])) * (sc[4] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
    wdot[15] -= qdot;
  }

  {
    // reaction 57:  CH3 + O2 <=> CH3O + O
    const amrex::Real k_f = 7546000 * exp(-(14251.0959755897) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[12]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[4] + g_RT[12] - g_RT[18])) * (sc[4] * sc[18]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[12] -= qdot;
    wdot[18] += qdot;
  }

  {
    // reaction 58:  CH3 + O2 <=> CH2O + OH
    const amrex::Real k_f =
      2.641e-06 * exp((3.283) * tc[0] - (4078.57107634725) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[12]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[6] + g_RT[12] - g_RT[20])) * (sc[6] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[12] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 59:  CH3 + O <=> CH2O + H
    const amrex::Real k_f =
      55400000 * exp((0.05) * tc[0] - (-68.4374665494419) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[12]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[4] + g_RT[12] - g_RT[20])) * (sc[2] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[12] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 60:  CH3 + OH <=> CH2(S) + H2O
    const amrex::Real k_f =
      4.788e+17 * exp((-3.155) * tc[0] - (3524.02631063045) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[12]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[12] - g_RT[14])) * (sc[5] * sc[14]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[12] -= qdot;
    wdot[14] += qdot;
  }

  {
    // reaction 61:  CH3 + OH <=> CH2O + H2
    const amrex::Real k_f =
      5374 * exp((0.287) * tc[0] - (140.900666425322) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[12]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[6] + g_RT[12] - g_RT[20])) * (sc[1] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[12] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 62:  CH3 + OH <=> CH2OH + H
    const amrex::Real k_f =
      15250000 * exp((0.134) * tc[0] - (2838.64521180442) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[12]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[12] - g_RT[19])) * (sc[2] * sc[19]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[6] -= qdot;
    wdot[12] -= qdot;
    wdot[19] += qdot;
  }

  {
    // reaction 63:  CH3 + OH <=> CH3O + H
    const amrex::Real k_f =
      1798 * exp((0.965) * tc[0] - (6068.79298960492) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[12]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[12] - g_RT[18])) * (sc[2] * sc[18]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[6] -= qdot;
    wdot[12] -= qdot;
    wdot[18] += qdot;
  }

  {
    // reaction 64:  CH3 + OH <=> CH2 + H2O
    const amrex::Real k_f =
      0.04293 * exp((2.568) * tc[0] - (2011.75958655411) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[12]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[12] - g_RT[13])) * (sc[5] * sc[13]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;
  }

  {
    // reaction 65:  CH3 + HO2 <=> CH3O + OH
    const amrex::Real k_f =
      1000000 * exp((0.269) * tc[0] - (-345.961457740745) * invT);
    const amrex::Real qf = k_f * (sc[8] * sc[12]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[6] + g_RT[8] + g_RT[12] - g_RT[18])) * (sc[6] * sc[18]);
    const amrex::Real qdot = qf - qr;
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[12] -= qdot;
    wdot[18] += qdot;
  }

  {
    // reaction 68:  CH3OH + H <=> CH3O + H2
    const amrex::Real k_f =
      0.199 * exp((2.56) * tc[0] - (5183.13165778861) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[17]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[17] - g_RT[18])) * (sc[1] * sc[18]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[17] -= qdot;
    wdot[18] += qdot;
  }

  {
    // reaction 69:  CH3OH + H <=> CH2OH + H2
    const amrex::Real k_f =
      0.307 * exp((2.55) * tc[0] - (2737.49866197768) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[17]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[17] - g_RT[19])) * (sc[1] * sc[19]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[17] -= qdot;
    wdot[19] += qdot;
  }

  {
    // reaction 70:  CH3OH + O <=> CH3O + OH
    const amrex::Real k_f =
      0.0388 * exp((2.5) * tc[0] - (1549.90733067854) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[17]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[4] - g_RT[6] + g_RT[17] - g_RT[18])) * (sc[6] * sc[18]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[17] -= qdot;
    wdot[18] += qdot;
  }

  {
    // reaction 71:  CH3OH + O <=> CH2OH + OH
    const amrex::Real k_f =
      0.388 * exp((2.5) * tc[0] - (1549.90733067854) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[17]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[4] - g_RT[6] + g_RT[17] - g_RT[19])) * (sc[6] * sc[19]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[17] -= qdot;
    wdot[19] += qdot;
  }

  {
    // reaction 72:  CH3OH + OH <=> CH3O + H2O
    const amrex::Real k_f =
      0.00015 * exp((3.03) * tc[0] - (-383.954316009001) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[17]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[17] - g_RT[18])) * (sc[5] * sc[18]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[17] -= qdot;
    wdot[18] += qdot;
  }

  {
    // reaction 73:  CH3OH + OH <=> CH2OH + H2O
    const amrex::Real k_f =
      0.0308 * exp((2.65) * tc[0] - (-405.944884304667) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[17]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[17] - g_RT[19])) * (sc[5] * sc[19]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[17] -= qdot;
    wdot[19] += qdot;
  }

  {
    // reaction 74:  CH3OH + O2 <=> CH3O + HO2
    const amrex::Real k_f =
      0.0358 * exp((2.27) * tc[0] - (21519.8091048059) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[17]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[8] + g_RT[17] - g_RT[18])) * (sc[8] * sc[18]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[8] += qdot;
    wdot[17] -= qdot;
    wdot[18] += qdot;
  }

  {
    // reaction 75:  CH3OH + HO2 <=> CH3O + H2O2
    const amrex::Real k_f = 1220000 * exp(-(10099.9107343668) * invT);
    const amrex::Real qf = k_f * (sc[8] * sc[17]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[7] + g_RT[8] + g_RT[17] - g_RT[18])) * (sc[7] * sc[18]);
    const amrex::Real qdot = qf - qr;
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[17] -= qdot;
    wdot[18] += qdot;
  }

  {
    // reaction 76:  CH3OH + HO2 <=> CH2OH + H2O2
    const amrex::Real k_f = 32600000 * exp(-(9451.5160604774) * invT);
    const amrex::Real qf = k_f * (sc[8] * sc[17]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[7] + g_RT[8] + g_RT[17] - g_RT[19])) * (sc[7] * sc[19]);
    const amrex::Real qdot = qf - qr;
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[17] -= qdot;
    wdot[19] += qdot;
  }

  {
    // reaction 77:  CH3 + CH3OH <=> CH2OH + CH4
    const amrex::Real k_f =
      2.13e-07 * exp((3.953) * tc[0] - (3550.24389891888) * invT);
    const amrex::Real qf = k_f * (sc[12] * sc[17]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[11] + g_RT[12] + g_RT[17] - g_RT[19])) *
                           (sc[11] * sc[19]);
    const amrex::Real qdot = qf - qr;
    wdot[11] += qdot;
    wdot[12] -= qdot;
    wdot[17] -= qdot;
    wdot[19] += qdot;
  }

  {
    // reaction 78:  CH3 + CH3OH <=> CH3O + CH4
    const amrex::Real k_f =
      0.00322 * exp((2.425) * tc[0] - (4317.34738427159) * invT);
    const amrex::Real qf = k_f * (sc[12] * sc[17]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[11] + g_RT[12] + g_RT[17] - g_RT[18])) *
                           (sc[11] * sc[18]);
    const amrex::Real qdot = qf - qr;
    wdot[11] += qdot;
    wdot[12] -= qdot;
    wdot[17] -= qdot;
    wdot[18] += qdot;
  }

  {
    // reaction 79:  CH2OH + O2 <=> CH2O + HO2
    const amrex::Real k_f = 1510000000 * exp((-1) * tc[0]);
    const amrex::Real qf = k_f * (sc[3] * sc[19]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[8] + g_RT[19] - g_RT[20])) * (sc[8] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[8] += qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 80:  CH2OH + O2 <=> CH2O + HO2
    const amrex::Real k_f = 241000000 * exp(-(2524.63801234228) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[19]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[8] + g_RT[19] - g_RT[20])) * (sc[8] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[8] += qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 81:  CH2OH + H <=> CH2O + H2
    const amrex::Real k_f = 6000000;
    const amrex::Real qf = k_f * (sc[2] * sc[19]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[19] - g_RT[20])) * (sc[1] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 82:  CH2OH + HO2 <=> CH2O + H2O2
    const amrex::Real k_f = 12000000;
    const amrex::Real qf = k_f * (sc[8] * sc[19]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[7] + g_RT[8] + g_RT[19] - g_RT[20])) * (sc[7] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 83:  CH2OH + HCO <=> 2 CH2O
    const amrex::Real k_f = 180000000;
    const amrex::Real qf = k_f * (sc[19] * sc[21]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[19] - 2.000000 * g_RT[20] + g_RT[21])) *
                           ((sc[20] * sc[20]));
    const amrex::Real qdot = qf - qr;
    wdot[19] -= qdot;
    wdot[20] += 2.000000 * qdot;
    wdot[21] -= qdot;
  }

  {
    // reaction 84:  CH2OH + OH <=> CH2O + H2O
    const amrex::Real k_f = 24000000;
    const amrex::Real qf = k_f * (sc[6] * sc[19]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[19] - g_RT[20])) * (sc[5] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 85:  CH2OH + O <=> CH2O + OH
    const amrex::Real k_f = 42000000;
    const amrex::Real qf = k_f * (sc[4] * sc[19]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[4] - g_RT[6] + g_RT[19] - g_RT[20])) * (sc[6] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 86:  CH3O + O2 <=> CH2O + HO2
    const amrex::Real k_f =
      4.38e-25 * exp((9.5) * tc[0] - (-2768.19487859176) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[18]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[8] + g_RT[18] - g_RT[20])) * (sc[8] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[8] += qdot;
    wdot[18] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 87:  CH3O + H <=> CH2O + H2
    const amrex::Real k_f = 20000000;
    const amrex::Real qf = k_f * (sc[2] * sc[18]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[18] - g_RT[20])) * (sc[1] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[18] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 88:  CH3O + HO2 <=> CH2O + H2O2
    const amrex::Real k_f = 301000;
    const amrex::Real qf = k_f * (sc[8] * sc[18]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[7] + g_RT[8] + g_RT[18] - g_RT[20])) * (sc[7] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[18] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 89:  CH3 + CH3O <=> CH2O + CH4
    const amrex::Real k_f = 12000000;
    const amrex::Real qf = k_f * (sc[12] * sc[18]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[11] + g_RT[12] + g_RT[18] - g_RT[20])) *
                           (sc[11] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[11] += qdot;
    wdot[12] -= qdot;
    wdot[18] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 90:  2 CH3O <=> CH2O + CH3OH
    const amrex::Real k_f = 60300000;
    const amrex::Real qf = k_f * ((sc[18] * sc[18]));
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[17] + 2.000000 * g_RT[18] - g_RT[20])) *
                           (sc[17] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[17] += qdot;
    wdot[18] -= 2.000000 * qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 91:  CH2O + O <=> HCO + OH
    const amrex::Real k_f =
      6260 * exp((1.15) * tc[0] - (1137.26966471867) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[20]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[4] - g_RT[6] + g_RT[20] - g_RT[21])) * (sc[6] * sc[21]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;
  }

  {
    // reaction 92:  CH2O + H <=> H2 + HCO
    const amrex::Real k_f =
      57.4 * exp((1.9) * tc[0] - (1378.81366430493) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[20]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] + g_RT[20] - g_RT[21])) * (sc[1] * sc[21]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;
  }

  {
    // reaction 93:  CH2O + OH <=> H2O + HCO
    const amrex::Real k_f =
      78.2 * exp((1.63) * tc[0] - (-530.893582423979) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[20]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[5] + g_RT[6] + g_RT[20] - g_RT[21])) * (sc[5] * sc[21]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;
  }

  {
    // reaction 94:  CH2O + HO2 <=> H2O2 + HCO
    const amrex::Real k_f =
      0.0188 * exp((2.7) * tc[0] - (5797.05599007037) * invT);
    const amrex::Real qf = k_f * (sc[8] * sc[20]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[7] + g_RT[8] + g_RT[20] - g_RT[21])) * (sc[7] * sc[21]);
    const amrex::Real qdot = qf - qr;
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;
  }

  {
    // reaction 95:  CH2O + CH3 <=> CH4 + HCO
    const amrex::Real k_f =
      3.83e-05 * exp((3.36) * tc[0] - (2169.87026294995) * invT);
    const amrex::Real qf = k_f * (sc[12] * sc[20]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[11] + g_RT[12] + g_RT[20] - g_RT[21])) *
                           (sc[11] * sc[21]);
    const amrex::Real qdot = qf - qr;
    wdot[11] += qdot;
    wdot[12] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;
  }

  {
    // reaction 96:  CH2O + CH3O <=> CH3OH + HCO
    const amrex::Real k_f = 662000 * exp(-(1154.37903135603) * invT);
    const amrex::Real qf = k_f * (sc[18] * sc[20]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[17] + g_RT[18] + g_RT[20] - g_RT[21])) *
                           (sc[17] * sc[21]);
    const amrex::Real qdot = qf - qr;
    wdot[17] += qdot;
    wdot[18] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;
  }

  {
    // reaction 98:  HCO + O2 <=> CO + HO2
    const amrex::Real k_f = 7580000 * exp(-(206.318832979935) * invT);
    const amrex::Real qf = k_f * (sc[3] * sc[21]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[3] - g_RT[8] - g_RT[9] + g_RT[21])) * (sc[8] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[3] -= qdot;
    wdot[8] += qdot;
    wdot[9] += qdot;
    wdot[21] -= qdot;
  }

  {
    // reaction 99:  HCO + O <=> CO + OH
    const amrex::Real k_f = 30200000;
    const amrex::Real qf = k_f * (sc[4] * sc[21]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[4] - g_RT[6] - g_RT[9] + g_RT[21])) * (sc[6] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[9] += qdot;
    wdot[21] -= qdot;
  }

  {
    // reaction 100:  H + HCO <=> CO + H2
    const amrex::Real k_f = 73400000;
    const amrex::Real qf = k_f * (sc[2] * sc[21]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[9] + g_RT[21])) * (sc[1] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[9] += qdot;
    wdot[21] -= qdot;
  }

  {
    // reaction 101:  HCO + OH <=> CO + H2O
    const amrex::Real k_f = 30110000;
    const amrex::Real qf = k_f * (sc[6] * sc[21]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[5] + g_RT[6] - g_RT[9] + g_RT[21])) * (sc[5] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[9] += qdot;
    wdot[21] -= qdot;
  }

  {
    // reaction 102:  CH3 + HCO <=> CH4 + CO
    const amrex::Real k_f = 26500000;
    const amrex::Real qf = k_f * (sc[12] * sc[21]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[9] - g_RT[11] + g_RT[12] + g_RT[21])) *
                           (sc[9] * sc[11]);
    const amrex::Real qdot = qf - qr;
    wdot[9] += qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;
    wdot[21] -= qdot;
  }

  {
    // reaction 103:  2 HCO <=> CH2O + CO
    const amrex::Real k_f = 18000000;
    const amrex::Real qf = k_f * ((sc[21] * sc[21]));
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[9] - g_RT[20] + 2.000000 * g_RT[21])) *
                           (sc[9] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[9] += qdot;
    wdot[20] += qdot;
    wdot[21] -= 2.000000 * qdot;
  }

  {
    // reaction 104:  HCO + O <=> CO2 + H
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[4] * sc[21]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[4] - g_RT[10] + g_RT[21])) * (sc[2] * sc[10]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[10] += qdot;
    wdot[21] -= qdot;
  }

  {
    // reaction 105:  HCO + HO2 => CO2 + H + OH
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[8] * sc[21]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[10] += qdot;
    wdot[21] -= qdot;
  }

  {
    // reaction 109:  CH + OH <=> C + H2O
    const amrex::Real k_f = 40 * exp((2) * tc[0] - (1509.64999741416) * invT);
    const amrex::Real qf = k_f * (sc[6] * sc[16]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[5] + g_RT[6] - g_RT[15] + g_RT[16])) * (sc[5] * sc[15]);
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;
  }

  {
    // reaction 110:  CH2 + HCO <=> CH3 + CO
    const amrex::Real k_f = 20000000;
    const amrex::Real qf = k_f * (sc[13] * sc[21]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[9] - g_RT[12] + g_RT[13] + g_RT[21])) *
                           (sc[9] * sc[12]);
    const amrex::Real qdot = qf - qr;
    wdot[9] += qdot;
    wdot[12] += qdot;
    wdot[13] -= qdot;
    wdot[21] -= qdot;
  }

  {
    // reaction 111:  CH2 + O <=> CO + H2
    const amrex::Real k_f = 50000000;
    const amrex::Real qf = k_f * (sc[4] * sc[13]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] + g_RT[4] - g_RT[9] + g_RT[13])) * (sc[1] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[9] += qdot;
    wdot[13] -= qdot;
  }

  {
    // reaction 112:  CH2 + OH <=> CH2O + H
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[6] * sc[13]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] + g_RT[6] + g_RT[13] - g_RT[20])) * (sc[2] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[6] -= qdot;
    wdot[13] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 113:  CH2 + CO2 <=> CH2O + CO
    const amrex::Real k_f = 110000 * exp(-(503.21666580472) * invT);
    const amrex::Real qf = k_f * (sc[10] * sc[13]);
    const amrex::Real qr = k_f *
                           exp(-(-g_RT[9] + g_RT[10] + g_RT[13] - g_RT[20])) *
                           (sc[9] * sc[20]);
    const amrex::Real qdot = qf - qr;
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[13] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 114:  CH2(S) + O <=> CO + 2 H
    const amrex::Real k_f = 30000000;
    const amrex::Real qf = k_f * (sc[4] * sc[14]);
    const amrex::Real qr =
      k_f * exp(-(-2.000000 * g_RT[2] + g_RT[4] - g_RT[9] + g_RT[14])) *
      (refCinv) * ((sc[2] * sc[2]) * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += 2.000000 * qdot;
    wdot[4] -= qdot;
    wdot[9] += qdot;
    wdot[14] -= qdot;
  }

  {
    // reaction 115:  CH2(S) + CH4 <=> 2 CH3
    const amrex::Real k_f = 43000000;
    const amrex::Real qf = k_f * (sc[11] * sc[14]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[11] - 2.000000 * g_RT[12] + g_RT[14])) *
                           ((sc[12] * sc[12]));
    const amrex::Real qdot = qf - qr;
    wdot[11] -= qdot;
    wdot[12] += 2.000000 * qdot;
    wdot[14] -= qdot;
  }

  {
    // reaction 116:  CH3OH + H => CH3 + H2O
    const amrex::Real k_f = 650000 * exp(-(2667.04832876501) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[17]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[12] += qdot;
    wdot[17] -= qdot;
  }

  {
    // reaction 117:  CH2O + O => CO2 + 2 H
    const amrex::Real k_f = 20000000 * exp(-(2516.0833290236) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[20]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[2] += 2.000000 * qdot;
    wdot[4] -= qdot;
    wdot[10] += qdot;
    wdot[20] -= qdot;
  }

  {
    // reaction 118:  CH2O + OH => CO2 + H + H2
    const amrex::Real k_f = 100000;
    const amrex::Real qf = k_f * (sc[6] * sc[20]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[2] += qdot;
    wdot[6] -= qdot;
    wdot[10] += qdot;
    wdot[20] -= qdot;
  }

  {
    // reaction 119:  CH3O + CO <=> CH3 + CO2
    const amrex::Real k_f = 500000 * exp(-(3270.90832773068) * invT);
    const amrex::Real qf = k_f * (sc[9] * sc[18]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[9] - g_RT[10] - g_RT[12] + g_RT[18])) *
                           (sc[10] * sc[12]);
    const amrex::Real qdot = qf - qr;
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[12] += qdot;
    wdot[18] -= qdot;
  }

  {
    // reaction 121:  HCO + HO2 <=> CO + H2O2
    const amrex::Real k_f = 400000;
    const amrex::Real qf = k_f * (sc[8] * sc[21]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[7] + g_RT[8] - g_RT[9] + g_RT[21])) * (sc[7] * sc[9]);
    const amrex::Real qdot = qf - qr;
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;
    wdot[21] -= qdot;
  }

  {
    // reaction 122:  CH3O + OH => CH2O + H2O
    const amrex::Real k_f = 15000000;
    const amrex::Real qf = k_f * (sc[6] * sc[18]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[18] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 123:  CH2OH + CH3 => CH2O + CH4
    const amrex::Real k_f = 15000000;
    const amrex::Real qf = k_f * (sc[12] * sc[19]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[11] += qdot;
    wdot[12] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;
  }

  {
    // reaction 124:  CH + O <=> CHO+ + E
    const amrex::Real k_f =
      1746000000000 * exp((-2.19) * tc[0] - (164.999209334044) * invT);
    const amrex::Real qf = k_f * (sc[4] * sc[16]);
    const amrex::Real qr = k_f *
                           exp(-(g_RT[4] + g_RT[16] - g_RT[22] - g_RT[25])) *
                           (sc[22] * sc[25]);
    const amrex::Real qdot = qf - qr;
    wdot[4] -= qdot;
    wdot[16] -= qdot;
    wdot[22] += qdot;
    wdot[25] += qdot;
  }

  {
    // reaction 125:  CHO+ + H2O <=> CO + H3O+
    const amrex::Real k_f = 26080000000 * exp((-0.5) * tc[0]);
    const amrex::Real qf = k_f * (sc[5] * sc[22]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[5] - g_RT[9] + g_RT[22] - g_RT[23])) * (sc[9] * sc[23]);
    const amrex::Real qdot = qf - qr;
    wdot[5] -= qdot;
    wdot[9] += qdot;
    wdot[22] -= qdot;
    wdot[23] += qdot;
  }

  {
    // reaction 126:  E + H3O+ <=> H + H2O
    const amrex::Real k_f = 739500000000 * exp((-0.5) * tc[0]);
    const amrex::Real qf = k_f * (sc[23] * sc[25]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[2] - g_RT[5] + g_RT[23] + g_RT[25])) * (sc[2] * sc[5]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += qdot;
    wdot[5] += qdot;
    wdot[23] -= qdot;
    wdot[25] -= qdot;
  }

  {
    // reaction 127:  E + H3O+ <=> 2 H + OH
    const amrex::Real k_f = 3181000000000 * exp((-0.5) * tc[0]);
    const amrex::Real qf = k_f * (sc[23] * sc[25]);
    const amrex::Real qr =
      k_f * exp(-(-2.000000 * g_RT[2] - g_RT[6] + g_RT[23] + g_RT[25])) *
      (refCinv) * ((sc[2] * sc[2]) * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[2] += 2.000000 * qdot;
    wdot[6] += qdot;
    wdot[23] -= qdot;
    wdot[25] -= qdot;
  }

  {
    // reaction 128:  E + H3O+ <=> H2 + OH
    const amrex::Real k_f = 560100000000 * exp((-0.5) * tc[0]);
    const amrex::Real qf = k_f * (sc[23] * sc[25]);
    const amrex::Real qr =
      k_f * exp(-(-g_RT[1] - g_RT[6] + g_RT[23] + g_RT[25])) * (sc[1] * sc[6]);
    const amrex::Real qdot = qf - qr;
    wdot[1] += qdot;
    wdot[6] += qdot;
    wdot[23] -= qdot;
    wdot[25] -= qdot;
  }

  {
    // reaction 129:  H2 + O2- <=> E + H2O2
    const amrex::Real k_f = 602000000;
    const amrex::Real qf = k_f * (sc[1] * sc[24]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[1] - g_RT[7] + g_RT[24] - g_RT[25])) * (sc[7] * sc[25]);
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[7] += qdot;
    wdot[24] -= qdot;
    wdot[25] += qdot;
  }

  {
    // reaction 130:  H + O2- <=> E + HO2
    const amrex::Real k_f = 723000000;
    const amrex::Real qf = k_f * (sc[2] * sc[24]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - g_RT[8] + g_RT[24] - g_RT[25])) * (sc[8] * sc[25]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[8] += qdot;
    wdot[24] -= qdot;
    wdot[25] += qdot;
  }
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(const amrex::Real* T, amrex::Real* C, amrex::Real* wdot)
{
  int id; // loop counter

  // convert to SI
  for (id = 0; id < 26; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, *T);

  // convert to chemkin units
  for (id = 0; id < 26; ++id) {
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
  amrex::Real c[26]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;
  amrex::Real imw[26];

  get_imw(imw);

  // Compute inverse of mean molecular wt first
  YOW += y[0] * imw[0];   // N2
  YOW += y[1] * imw[1];   // H2
  YOW += y[2] * imw[2];   // H
  YOW += y[3] * imw[3];   // O2
  YOW += y[4] * imw[4];   // O
  YOW += y[5] * imw[5];   // H2O
  YOW += y[6] * imw[6];   // OH
  YOW += y[7] * imw[7];   // H2O2
  YOW += y[8] * imw[8];   // HO2
  YOW += y[9] * imw[9];   // CO
  YOW += y[10] * imw[10]; // CO2
  YOW += y[11] * imw[11]; // CH4
  YOW += y[12] * imw[12]; // CH3
  YOW += y[13] * imw[13]; // CH2
  YOW += y[14] * imw[14]; // CH2(S)
  YOW += y[15] * imw[15]; // C
  YOW += y[16] * imw[16]; // CH
  YOW += y[17] * imw[17]; // CH3OH
  YOW += y[18] * imw[18]; // CH3O
  YOW += y[19] * imw[19]; // CH2OH
  YOW += y[20] * imw[20]; // CH2O
  YOW += y[21] * imw[21]; // HCO
  YOW += y[22] * imw[22]; // CHO+
  YOW += y[23] * imw[23]; // H3O+
  YOW += y[24] * imw[24]; // O2-
  YOW += y[25] * imw[25]; // E
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
  c[12] = PWORT * y[12] * imw[12];
  c[13] = PWORT * y[13] * imw[13];
  c[14] = PWORT * y[14] * imw[14];
  c[15] = PWORT * y[15] * imw[15];
  c[16] = PWORT * y[16] * imw[16];
  c[17] = PWORT * y[17] * imw[17];
  c[18] = PWORT * y[18] * imw[18];
  c[19] = PWORT * y[19] * imw[19];
  c[20] = PWORT * y[20] * imw[20];
  c[21] = PWORT * y[21] * imw[21];
  c[22] = PWORT * y[22] * imw[22];
  c[23] = PWORT * y[23] * imw[23];
  c[24] = PWORT * y[24] * imw[24];
  c[25] = PWORT * y[25] * imw[25];

  // convert to chemkin units
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 26; ++id) {
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
  amrex::Real c[26]; // temporary storage
  amrex::Real PORT =
    1e6 * (*P) /
    (8.31446261815324e+07 * (*T)); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (id = 0; id < 26; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 26; ++id) {
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
  amrex::Real c[26]; // temporary storage
  amrex::Real imw[26];

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
  c[12] = 1e6 * (*rho) * y[12] * imw[12];
  c[13] = 1e6 * (*rho) * y[13] * imw[13];
  c[14] = 1e6 * (*rho) * y[14] * imw[14];
  c[15] = 1e6 * (*rho) * y[15] * imw[15];
  c[16] = 1e6 * (*rho) * y[16] * imw[16];
  c[17] = 1e6 * (*rho) * y[17] * imw[17];
  c[18] = 1e6 * (*rho) * y[18] * imw[18];
  c[19] = 1e6 * (*rho) * y[19] * imw[19];
  c[20] = 1e6 * (*rho) * y[20] * imw[20];
  c[21] = 1e6 * (*rho) * y[21] * imw[21];
  c[22] = 1e6 * (*rho) * y[22] * imw[22];
  c[23] = 1e6 * (*rho) * y[23] * imw[23];
  c[24] = 1e6 * (*rho) * y[24] * imw[24];
  c[25] = 1e6 * (*rho) * y[25] * imw[25];

  // call productionRate
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 26; ++id) {
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
  amrex::Real c[26];  // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 28.014000;  // N2
  XW += x[1] * 2.016000;   // H2
  XW += x[2] * 1.008000;   // H
  XW += x[3] * 31.998000;  // O2
  XW += x[4] * 15.999000;  // O
  XW += x[5] * 18.015000;  // H2O
  XW += x[6] * 17.007000;  // OH
  XW += x[7] * 34.014000;  // H2O2
  XW += x[8] * 33.006000;  // HO2
  XW += x[9] * 28.010000;  // CO
  XW += x[10] * 44.009000; // CO2
  XW += x[11] * 16.043000; // CH4
  XW += x[12] * 15.035000; // CH3
  XW += x[13] * 14.027000; // CH2
  XW += x[14] * 14.027000; // CH2(S)
  XW += x[15] * 12.011000; // C
  XW += x[16] * 13.019000; // CH
  XW += x[17] * 32.042000; // CH3OH
  XW += x[18] * 31.034000; // CH3O
  XW += x[19] * 31.034000; // CH2OH
  XW += x[20] * 30.026000; // CH2O
  XW += x[21] * 29.018000; // HCO
  XW += x[22] * 29.017451; // CHO+
  XW += x[23] * 19.022451; // H3O+
  XW += x[24] * 31.998549; // O2-
  XW += x[25] * 0.000549;  // E
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * (*rho) / XW;

  // Compute conversion, see Eq 11
  for (id = 0; id < 26; ++id) {
    c[id] = x[id] * ROW;
  }

  // convert to chemkin units
  productionRate(wdot, c, *T);

  // convert to chemkin units
  for (id = 0; id < 26; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int* kcharge)
{
  kcharge[0] = 0;   // N2
  kcharge[1] = 0;   // H2
  kcharge[2] = 0;   // H
  kcharge[3] = 0;   // O2
  kcharge[4] = 0;   // O
  kcharge[5] = 0;   // H2O
  kcharge[6] = 0;   // OH
  kcharge[7] = 0;   // H2O2
  kcharge[8] = 0;   // HO2
  kcharge[9] = 0;   // CO
  kcharge[10] = 0;  // CO2
  kcharge[11] = 0;  // CH4
  kcharge[12] = 0;  // CH3
  kcharge[13] = 0;  // CH2
  kcharge[14] = 0;  // CH2(S)
  kcharge[15] = 0;  // C
  kcharge[16] = 0;  // CH
  kcharge[17] = 0;  // CH3OH
  kcharge[18] = 0;  // CH3O
  kcharge[19] = 0;  // CH2OH
  kcharge[20] = 0;  // CH2O
  kcharge[21] = 0;  // HCO
  kcharge[22] = 1;  // CHO+
  kcharge[23] = 1;  // H3O+
  kcharge[24] = -1; // O2-
  kcharge[25] = -1; // E
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real* zk)
{
  amrex::Real imw[26];
  get_imw(imw);

  int kchrg[26];
  CKCHRG(kchrg);

  for (int id = 0; id < 26; ++id) {
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

  // species with midpoint at T=1050 kelvin
  if (T < 1050) {
    // species 0: N2
    species[0] = -1.83159730e-03 + 8.64649324e-06 * tc[1] -
                 9.13134453e-09 * tc[2] + 2.98428625e-12 * tc[3];
  } else {
    // species 0: N2
    species[0] = +1.67067353e-03 - 1.35999486e-06 * tc[1] +
                 3.98644137e-10 * tc[2] - 4.11069768e-14 * tc[3];
  }

  // species with midpoint at T=700 kelvin
  if (T < 700) {
    // species 1: H2
    species[1] = +5.49529274e-03 - 2.54327268e-05 * tc[1] +
                 3.86247519e-08 * tc[2] - 1.88011100e-11 * tc[3];
    // species 3: O2
    species[3] = -2.79740147e-03 + 1.96024512e-05 * tc[1] -
                 3.09778929e-08 * tc[2] + 1.51972499e-11 * tc[3];
    // species 8: HO2
    species[8] = +1.05805705e-03 + 1.01335788e-05 * tc[1] -
                 1.90140229e-08 * tc[2] + 9.66389124e-12 * tc[3];
    // species 11: CH4
    species[11] = -1.46835123e-02 + 1.05946542e-04 * tc[1] -
                  1.62500647e-07 * tc[2] + 7.85274264e-11 * tc[3];
    // species 15: C
    species[15] = -2.47725281e-04 + 1.09603655e-06 * tc[1] -
                  1.64565375e-09 * tc[2] + 8.16469324e-13 * tc[3];
    // species 20: CH2O
    species[20] = -7.01151853e-03 + 6.30353922e-05 * tc[1] -
                  1.00943591e-07 * tc[2] + 4.93816092e-11 * tc[3];
  } else {
    // species 1: H2
    species[1] = -1.01873259e-03 + 2.48452466e-06 * tc[1] -
                 1.25703569e-09 * tc[2] + 1.90217517e-13 * tc[3];
    // species 3: O2
    species[3] = +2.48211357e-03 - 3.02404188e-06 * tc[1] +
                 1.34566860e-09 * tc[2] - 1.94922267e-13 * tc[3];
    // species 8: HO2
    species[8] = +4.46390907e-03 - 4.46292984e-06 * tc[1] +
                 1.83813240e-09 * tc[2] - 2.65706495e-13 * tc[3];
    // species 11: CH4
    species[11] = +1.23697845e-02 - 9.99615844e-06 * tc[1] +
                  3.13178295e-09 * tc[2] - 3.45158966e-13 * tc[3];
    // species 15: C
    species[15] = +3.92839476e-05 - 1.34002996e-07 * tc[1] +
                  1.11545608e-10 * tc[2] - 2.02922754e-14 * tc[3];
    // species 20: CH2O
    species[20] = +1.00905183e-02 - 1.02590512e-05 * tc[1] +
                  3.76275621e-09 * tc[2] - 4.78556436e-13 * tc[3];
  }

  // species with midpoint at T=860 kelvin
  if (T < 860) {
    // species 2: H
    species[2] = +2.99164046e-06 - 1.18551952e-08 * tc[1] +
                 1.46343049e-11 * tc[2] - 5.82157280e-15 * tc[3];
  } else {
    // species 2: H
    species[2] = -7.73406828e-07 + 1.27869069e-09 * tc[1] -
                 6.37655373e-13 * tc[2] + 9.77916764e-17 * tc[3];
  }

  // species with midpoint at T=720 kelvin
  if (T < 720) {
    // species 4: O
    species[4] = -3.11174065e-03 + 1.23627579e-05 * tc[1] -
                 1.69142639e-08 * tc[2] + 7.79464064e-12 * tc[3];
  } else {
    // species 4: O
    species[4] = -2.08959644e-04 + 2.67837092e-07 * tc[1] -
                 1.15762769e-10 * tc[2] + 1.75567476e-14 * tc[3];
  }

  // species with midpoint at T=1420 kelvin
  if (T < 1420) {
    // species 5: H2O
    species[5] = -8.65807189e-04 + 6.48819056e-06 * tc[1] -
                 5.40729237e-09 * tc[2] + 1.32993317e-12 * tc[3];
  } else {
    // species 5: H2O
    species[5] = +3.05768849e-03 - 1.80088482e-06 * tc[1] +
                 4.30084656e-10 * tc[2] - 4.03431268e-14 * tc[3];
  }

  // species with midpoint at T=1700 kelvin
  if (T < 1700) {
    // species 6: OH
    species[6] = -1.66275926e-03 + 4.61840058e-06 * tc[1] -
                 3.07078524e-09 * tc[2] + 6.35318516e-13 * tc[3];
  } else {
    // species 6: OH
    species[6] = +1.66635279e-03 - 1.25650303e-06 * tc[1] +
                 3.85040418e-10 * tc[2] - 4.22943576e-14 * tc[3];
  }

  // species with midpoint at T=1800 kelvin
  if (T < 1800) {
    // species 7: H2O2
    species[7] = +7.05005437e-03 - 7.69044012e-06 * tc[1] +
                 3.50161983e-09 * tc[2] - 5.90472420e-13 * tc[3];
    // species 13: CH2
    species[13] = +1.43839191e-03 + 9.51166154e-07 * tc[1] -
                  1.29536577e-09 * tc[2] + 3.03317150e-13 * tc[3];
    // species 17: CH3OH
    species[17] = +1.63086904e-02 - 1.69668992e-05 * tc[1] +
                  6.87913023e-09 * tc[2] - 1.03980805e-12 * tc[3];
  } else {
    // species 7: H2O2
    species[7] = +3.89237848e-03 - 2.42764698e-06 * tc[1] +
                 5.77845855e-10 * tc[2] - 4.90327960e-14 * tc[3];
    // species 13: CH2
    species[13] = +3.55431388e-03 - 2.57537046e-06 * tc[1] +
                  6.63821232e-10 * tc[2] - 5.94952588e-14 * tc[3];
    // species 17: CH3OH
    species[17] = +1.21538353e-02 - 1.00421406e-05 * tc[1] +
                  3.03204210e-09 * tc[2] - 3.27384329e-13 * tc[3];
  }

  // species with midpoint at T=960 kelvin
  if (T < 960) {
    // species 9: CO
    species[9] = -2.14465241e-03 + 1.08415801e-05 * tc[1] -
                 1.25107789e-08 * tc[2] + 4.47604508e-12 * tc[3];
  } else {
    // species 9: CO
    species[9] = +1.87486886e-03 - 1.71942385e-06 * tc[1] +
                 5.73600210e-10 * tc[2] - 6.71421144e-14 * tc[3];
  }

  // species with midpoint at T=1450 kelvin
  if (T < 1450) {
    // species 10: CO2
    species[10] = +9.22755036e-03 - 1.55130819e-05 * tc[1] +
                  9.84676080e-09 * tc[2] - 2.19488993e-12 * tc[3];
  } else {
    // species 10: CO2
    species[10] = +2.62914704e-03 - 1.86121292e-06 * tc[1] +
                  4.31678760e-10 * tc[2] - 3.05032566e-14 * tc[3];
  }

  // species with midpoint at T=1060 kelvin
  if (T < 1060) {
    // species 12: CH3
    species[12] = +3.54764773e-03 + 2.94816880e-06 * tc[1] -
                  5.83127865e-09 * tc[2] + 2.08768493e-12 * tc[3];
  } else {
    // species 12: CH3
    species[12] = +6.15233477e-03 - 4.42358698e-06 * tc[1] +
                  1.12320794e-09 * tc[2] - 9.92605396e-14 * tc[3];
  }

  // species with midpoint at T=970 kelvin
  if (T < 970) {
    // species 14: CH2(S)
    species[14] = -2.21134314e-03 + 1.54305508e-05 * tc[1] -
                  1.78785114e-08 * tc[2] + 6.33258512e-12 * tc[3];
  } else {
    // species 14: CH2(S)
    species[14] = +3.65468307e-03 - 2.71179828e-06 * tc[1] +
                  8.24941233e-10 * tc[2] - 9.47181888e-14 * tc[3];
  }

  // species with midpoint at T=1590 kelvin
  if (T < 1590) {
    // species 16: CH
    species[16] = -1.58547350e-03 + 5.67024476e-06 * tc[1] -
                  4.08438174e-09 * tc[2] + 8.95981328e-13 * tc[3];
  } else {
    // species 16: CH
    species[16] = +2.16985238e-03 - 1.41527577e-06 * tc[1] +
                  3.71920482e-10 * tc[2] - 3.82539406e-14 * tc[3];
  }

  // species with midpoint at T=1740 kelvin
  if (T < 1740) {
    // species 18: CH3O
    species[18] = +1.70119767e-02 - 2.27614702e-05 * tc[1] +
                  1.16484278e-08 * tc[2] - 2.13136592e-12 * tc[3];
  } else {
    // species 18: CH3O
    species[18] = +5.90227637e-03 - 3.60681440e-06 * tc[1] +
                  6.40005027e-10 * tc[2] - 2.24726564e-14 * tc[3];
  }

  // species with midpoint at T=1360 kelvin
  if (T < 1360) {
    // species 19: CH2OH
    species[19] = +1.39600168e-02 - 2.17264412e-05 * tc[1] +
                  1.38749525e-08 * tc[2] - 3.23399665e-12 * tc[3];
  } else {
    // species 19: CH2OH
    species[19] = +6.02727059e-03 - 4.22773642e-06 * tc[1] +
                  1.00825771e-09 * tc[2] - 8.03949928e-14 * tc[3];
  }

  // species with midpoint at T=770 kelvin
  if (T < 770) {
    // species 21: HCO
    species[21] = -2.15836864e-03 + 2.36467750e-05 * tc[1] -
                  3.55378218e-08 * tc[2] + 1.60237582e-11 * tc[3];
  } else {
    // species 21: HCO
    species[21] = +5.29278258e-03 - 5.38368422e-06 * tc[1] +
                  2.16407339e-09 * tc[2] - 2.97408564e-13 * tc[3];
  }

  // species with midpoint at T=3654.18 kelvin
  if (T < 3654.18) {
    // species 22: CHO+
    species[22] = +4.94860200e-03 - 4.61706120e-06 * tc[1] +
                  1.53365256e-09 * tc[2] - 1.74082524e-13 * tc[3];
  } else {
    // species 22: CHO+
    species[22] = +1.66311070e-04 - 3.69526860e-08 * tc[1] +
                  2.72569287e-12 * tc[2] - 6.58177640e-17 * tc[3];
  }

  // species with midpoint at T=5430.9 kelvin
  if (T < 5430.9) {
    // species 23: H3O+
    species[23] = +5.44536540e-03 - 3.41111320e-06 * tc[1] +
                  7.44998790e-10 * tc[2] - 5.50777600e-14 * tc[3];
  } else {
    // species 23: H3O+
    species[23] = +2.52816680e-04 - 4.91887940e-08 * tc[1] +
                  3.26349540e-12 * tc[2] - 7.23633480e-17 * tc[3];
  }

  // species with midpoint at T=2008.71 kelvin
  if (T < 2008.71) {
    // species 24: O2-
    species[24] = +2.79808750e-03 - 4.53022520e-06 * tc[1] +
                  2.60749551e-09 * tc[2] - 5.08875360e-13 * tc[3];
  } else {
    // species 24: O2-
    species[24] = +2.24680720e-04 - 1.02795910e-07 * tc[1] +
                  2.20637934e-11 * tc[2] - 1.54234608e-15 * tc[3];
  }

  // species with midpoint at T=27300 kelvin
  if (T < 27300) {
    // species 25: E
    species[25] = -8.15026080e-18 + 2.60597340e-21 * tc[1] -
                  2.17944999e-25 * tc[2] + 5.22649880e-30 * tc[3];
  } else {
    // species 25: E
    species[25] = +0.00000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3];
  }
}

// compute an approx to the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian_precond(
  amrex::Real* J, const amrex::Real* sc, const amrex::Real T, const int HP)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  for (int i = 0; i < 729; i++) {
    J[i] = 0.0;
  }

  amrex::Real wdot[26];
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
  for (int k = 0; k < 26; ++k) {
    mixture += sc[k];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[26];
  gibbs(g_RT, tc);

  // compute the species enthalpy
  amrex::Real h_RT[26];
  speciesEnthalpy(h_RT, tc);

  amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  amrex::Real dqdci, dcdc_fac, dqdc[26];
  amrex::Real Pr, fPr, F, k_0, logPr;
  amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
  amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
  amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const amrex::Real ln10 = log(10.0);
  const amrex::Real log10e = 1.0 / log(10.0);
  // reaction 6: H2O2 (+M) <=> 2 OH (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (1.5 - 1) * sc[0] + (3.7 - 1) * sc[1] + (1.2 - 1) * sc[3] +
          (7.65 - 1) * sc[5] + (7.7 - 1) * sc[7] + (2.8 - 1) * sc[9] +
          (1.6 - 1) * sc[10];
  // forward
  phi_f = sc[7];
  k_f = 2000000000000 * exp(0.9 * tc[0] - (24531.3092413143) * invT);
  dlnkfdT = 0.9 * invT + (24531.3092413143) * invT2;
  // pressure-fall-off
  k_0 = 2.49e+24 * exp(-2.3 * tc[0] - (24531.3092413143) * invT);
  Pr = 1e-6 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -2.3 * invT + (24531.3092413143) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.43) * exp(-T / 1e-30);
  Fcent2 = 0.43 * exp(-T / 1e+30);
  Fcent3 = exp(-1e-16 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT = log10e / Fcent *
                (-Fcent1 / 1e-30 - Fcent2 / 1e+30 + Fcent3 * 1e-16 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = (sc[6] * sc[6]);
  Kc = refC * exp(-2.000000 * g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[7]) + (2.000000 * h_RT[6]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[6] += 2 * q; // OH
  wdot[7] -= q;     // H2O2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = 1.5 * dcdc_fac;
  dqdc[1] = 3.7 * dcdc_fac;
  dqdc[2] = dcdc_fac;
  dqdc[3] = 1.2 * dcdc_fac;
  dqdc[4] = dcdc_fac;
  dqdc[5] = 7.65 * dcdc_fac;
  dqdc[6] = dcdc_fac - k_r * 2.000000 * sc[6];
  dqdc[7] = 7.7 * dcdc_fac + k_f;
  dqdc[8] = dcdc_fac;
  dqdc[9] = 2.8 * dcdc_fac;
  dqdc[10] = 1.6 * dcdc_fac;
  dqdc[11] = dcdc_fac;
  dqdc[12] = dcdc_fac;
  dqdc[13] = dcdc_fac;
  dqdc[14] = dcdc_fac;
  dqdc[15] = dcdc_fac;
  dqdc[16] = dcdc_fac;
  dqdc[17] = dcdc_fac;
  dqdc[18] = dcdc_fac;
  dqdc[19] = dcdc_fac;
  dqdc[20] = dcdc_fac;
  dqdc[21] = dcdc_fac;
  dqdc[22] = dcdc_fac;
  dqdc[23] = dcdc_fac;
  dqdc[24] = dcdc_fac;
  dqdc[25] = dcdc_fac;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 6] += 2 * dqdc[k];
    J[27 * k + 7] -= dqdc[k];
  }
  J[708] += 2 * dqdT; // dwdot[OH]/dT
  J[709] -= dqdT;     // dwdot[H2O2]/dT

  // reaction 19: H + O2 (+M) <=> HO2 (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (1.3 - 1) * sc[1] + (10 - 1) * sc[5] + (1.9 - 1) * sc[9] +
          (3.8 - 1) * sc[10] + sc[11];
  // forward
  phi_f = sc[2] * sc[3];
  k_f = 4650000 * exp(0.44 * tc[0]);
  dlnkfdT = 0.44 * invT;
  // pressure-fall-off
  k_0 = 1.74e+19 * exp(-1.23 * tc[0]);
  Pr = 1e-12 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -1.23 * invT;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.67) * exp(-T / 1e-30);
  Fcent2 = 0.67 * exp(-T / 1e+30);
  Fcent3 = exp(-1e+30 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT = log10e / Fcent *
                (-Fcent1 / 1e-30 - Fcent2 / 1e+30 + Fcent3 * 1e+30 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[8];
  Kc = refCinv * exp(g_RT[2] + g_RT[3] - g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[8]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[2] -= q; // H
  wdot[3] -= q; // O2
  wdot[8] += q; // HO2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = dcdc_fac;
  dqdc[1] = 1.3 * dcdc_fac;
  dqdc[2] = dcdc_fac + k_f * sc[3];
  dqdc[3] = dcdc_fac + k_f * sc[2];
  dqdc[4] = dcdc_fac;
  dqdc[5] = 10 * dcdc_fac;
  dqdc[6] = dcdc_fac;
  dqdc[7] = dcdc_fac;
  dqdc[8] = dcdc_fac - k_r;
  dqdc[9] = 1.9 * dcdc_fac;
  dqdc[10] = 3.8 * dcdc_fac;
  dqdc[11] = 2 * dcdc_fac;
  dqdc[12] = dcdc_fac;
  dqdc[13] = dcdc_fac;
  dqdc[14] = dcdc_fac;
  dqdc[15] = dcdc_fac;
  dqdc[16] = dcdc_fac;
  dqdc[17] = dcdc_fac;
  dqdc[18] = dcdc_fac;
  dqdc[19] = dcdc_fac;
  dqdc[20] = dcdc_fac;
  dqdc[21] = dcdc_fac;
  dqdc[22] = dcdc_fac;
  dqdc[23] = dcdc_fac;
  dqdc[24] = dcdc_fac;
  dqdc[25] = dcdc_fac;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 2] -= dqdc[k];
    J[27 * k + 3] -= dqdc[k];
    J[27 * k + 8] += dqdc[k];
  }
  J[704] -= dqdT; // dwdot[H]/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT

  // reaction 24: CH3 + H (+M) <=> CH4 (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha =
    mixture + sc[1] + (6 - 1) * sc[5] + (1.5 - 1) * sc[9] + sc[10] + sc[11];
  // forward
  phi_f = sc[2] * sc[12];
  k_f = 12700000000 * exp(-0.63 * tc[0] - (192.731983003208) * invT);
  dlnkfdT = -0.63 * invT + (192.731983003208) * invT2;
  // pressure-fall-off
  k_0 = 2.48e+33 * exp(-4.76 * tc[0] - (1227.84866456352) * invT);
  Pr = 1e-12 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -4.76 * invT + (1227.84866456352) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.783) * exp(-T / 74);
  Fcent2 = 0.783 * exp(-T / 2941);
  Fcent3 = exp(-6964 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT =
    log10e / Fcent * (-Fcent1 / 74 - Fcent2 / 2941 + Fcent3 * 6964 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[11];
  Kc = refCinv * exp(g_RT[2] - g_RT[11] + g_RT[12]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[12]) + (h_RT[11]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[2] -= q;  // H
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = dcdc_fac;
  dqdc[1] = 2 * dcdc_fac;
  dqdc[2] = dcdc_fac + k_f * sc[12];
  dqdc[3] = dcdc_fac;
  dqdc[4] = dcdc_fac;
  dqdc[5] = 6 * dcdc_fac;
  dqdc[6] = dcdc_fac;
  dqdc[7] = dcdc_fac;
  dqdc[8] = dcdc_fac;
  dqdc[9] = 1.5 * dcdc_fac;
  dqdc[10] = 2 * dcdc_fac;
  dqdc[11] = 2 * dcdc_fac - k_r;
  dqdc[12] = dcdc_fac + k_f * sc[2];
  dqdc[13] = dcdc_fac;
  dqdc[14] = dcdc_fac;
  dqdc[15] = dcdc_fac;
  dqdc[16] = dcdc_fac;
  dqdc[17] = dcdc_fac;
  dqdc[18] = dcdc_fac;
  dqdc[19] = dcdc_fac;
  dqdc[20] = dcdc_fac;
  dqdc[21] = dcdc_fac;
  dqdc[22] = dcdc_fac;
  dqdc[23] = dcdc_fac;
  dqdc[24] = dcdc_fac;
  dqdc[25] = dcdc_fac;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 2] -= dqdc[k];
    J[27 * k + 11] += dqdc[k];
    J[27 * k + 12] -= dqdc[k];
  }
  J[704] -= dqdT; // dwdot[H]/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT

  // reaction 43: CH2 + H (+M) <=> CH3 (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha =
    mixture + sc[1] + (6 - 1) * sc[5] + (1.5 - 1) * sc[9] + sc[10] + sc[11];
  // forward
  phi_f = sc[2] * sc[13];
  k_f = 25000000000 * exp(-0.8 * tc[0]);
  dlnkfdT = -0.8 * invT;
  // pressure-fall-off
  k_0 = 3.2e+27 * exp(-3.14 * tc[0] - (618.956498939805) * invT);
  Pr = 1e-12 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -3.14 * invT + (618.956498939805) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.68) * exp(-T / 78);
  Fcent2 = 0.68 * exp(-T / 1995);
  Fcent3 = exp(-5590 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT =
    log10e / Fcent * (-Fcent1 / 78 - Fcent2 / 1995 + Fcent3 * 5590 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[12];
  Kc = refCinv * exp(g_RT[2] - g_RT[12] + g_RT[13]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[12]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[2] -= q;  // H
  wdot[12] += q; // CH3
  wdot[13] -= q; // CH2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = dcdc_fac;
  dqdc[1] = 2 * dcdc_fac;
  dqdc[2] = dcdc_fac + k_f * sc[13];
  dqdc[3] = dcdc_fac;
  dqdc[4] = dcdc_fac;
  dqdc[5] = 6 * dcdc_fac;
  dqdc[6] = dcdc_fac;
  dqdc[7] = dcdc_fac;
  dqdc[8] = dcdc_fac;
  dqdc[9] = 1.5 * dcdc_fac;
  dqdc[10] = 2 * dcdc_fac;
  dqdc[11] = 2 * dcdc_fac;
  dqdc[12] = dcdc_fac - k_r;
  dqdc[13] = dcdc_fac + k_f * sc[2];
  dqdc[14] = dcdc_fac;
  dqdc[15] = dcdc_fac;
  dqdc[16] = dcdc_fac;
  dqdc[17] = dcdc_fac;
  dqdc[18] = dcdc_fac;
  dqdc[19] = dcdc_fac;
  dqdc[20] = dcdc_fac;
  dqdc[21] = dcdc_fac;
  dqdc[22] = dcdc_fac;
  dqdc[23] = dcdc_fac;
  dqdc[24] = dcdc_fac;
  dqdc[25] = dcdc_fac;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 2] -= dqdc[k];
    J[27 * k + 12] += dqdc[k];
    J[27 * k + 13] -= dqdc[k];
  }
  J[704] -= dqdT; // dwdot[H]/dT
  J[714] += dqdT; // dwdot[CH3]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT

  // reaction 66: CH3OH (+M) <=> CH3 + OH (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture;
  // forward
  phi_f = sc[17];
  k_f = 2.084e+18 * exp(-0.615 * tc[0] - (46567.9721835682) * invT);
  dlnkfdT = -0.615 * invT + (46567.9721835682) * invT2;
  // pressure-fall-off
  k_0 = 1.5e+43 * exp(-6.995 * tc[0] - (49311.3081588693) * invT);
  Pr = 1e-6 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -6.995 * invT + (49311.3081588693) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. + 0.4748) * exp(-T / 35580);
  Fcent2 = -0.4748 * exp(-T / 1116);
  Fcent3 = exp(-9023 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT =
    log10e / Fcent * (-Fcent1 / 35580 - Fcent2 / 1116 + Fcent3 * 9023 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[6] * sc[12];
  Kc = refC * exp(-g_RT[6] - g_RT[12] + g_RT[17]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[17]) + (h_RT[6] + h_RT[12]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[6] += q;  // OH
  wdot[12] += q; // CH3
  wdot[17] -= q; // CH3OH
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = dcdc_fac;
  dqdc[1] = dcdc_fac;
  dqdc[2] = dcdc_fac;
  dqdc[3] = dcdc_fac;
  dqdc[4] = dcdc_fac;
  dqdc[5] = dcdc_fac;
  dqdc[6] = dcdc_fac - k_r * sc[12];
  dqdc[7] = dcdc_fac;
  dqdc[8] = dcdc_fac;
  dqdc[9] = dcdc_fac;
  dqdc[10] = dcdc_fac;
  dqdc[11] = dcdc_fac;
  dqdc[12] = dcdc_fac - k_r * sc[6];
  dqdc[13] = dcdc_fac;
  dqdc[14] = dcdc_fac;
  dqdc[15] = dcdc_fac;
  dqdc[16] = dcdc_fac;
  dqdc[17] = dcdc_fac + k_f;
  dqdc[18] = dcdc_fac;
  dqdc[19] = dcdc_fac;
  dqdc[20] = dcdc_fac;
  dqdc[21] = dcdc_fac;
  dqdc[22] = dcdc_fac;
  dqdc[23] = dcdc_fac;
  dqdc[24] = dcdc_fac;
  dqdc[25] = dcdc_fac;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 6] += dqdc[k];
    J[27 * k + 12] += dqdc[k];
    J[27 * k + 17] -= dqdc[k];
  }
  J[708] += dqdT; // dwdot[OH]/dT
  J[714] += dqdT; // dwdot[CH3]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT

  // reaction 67: CH3OH (+M) <=> CH2(S) + H2O (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture;
  // forward
  phi_f = sc[17];
  k_f = 3.121e+18 * exp(-1.017 * tc[0] - (46151.0068542825) * invT);
  dlnkfdT = -1.017 * invT + (46151.0068542825) * invT2;
  // pressure-fall-off
  k_0 = 1.43e+47 * exp(-8.227 * tc[0] - (50028.3415859744) * invT);
  Pr = 1e-6 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -8.227 * invT + (50028.3415859744) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 2.545) * exp(-T / 3290);
  Fcent2 = 2.545 * exp(-T / 47320);
  Fcent3 = exp(-47110 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT =
    log10e / Fcent * (-Fcent1 / 3290 - Fcent2 / 47320 + Fcent3 * 47110 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[5] * sc[14];
  Kc = refC * exp(-g_RT[5] - g_RT[14] + g_RT[17]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[17]) + (h_RT[5] + h_RT[14]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[5] += q;  // H2O
  wdot[14] += q; // CH2(S)
  wdot[17] -= q; // CH3OH
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = dcdc_fac;
  dqdc[1] = dcdc_fac;
  dqdc[2] = dcdc_fac;
  dqdc[3] = dcdc_fac;
  dqdc[4] = dcdc_fac;
  dqdc[5] = dcdc_fac - k_r * sc[14];
  dqdc[6] = dcdc_fac;
  dqdc[7] = dcdc_fac;
  dqdc[8] = dcdc_fac;
  dqdc[9] = dcdc_fac;
  dqdc[10] = dcdc_fac;
  dqdc[11] = dcdc_fac;
  dqdc[12] = dcdc_fac;
  dqdc[13] = dcdc_fac;
  dqdc[14] = dcdc_fac - k_r * sc[5];
  dqdc[15] = dcdc_fac;
  dqdc[16] = dcdc_fac;
  dqdc[17] = dcdc_fac + k_f;
  dqdc[18] = dcdc_fac;
  dqdc[19] = dcdc_fac;
  dqdc[20] = dcdc_fac;
  dqdc[21] = dcdc_fac;
  dqdc[22] = dcdc_fac;
  dqdc[23] = dcdc_fac;
  dqdc[24] = dcdc_fac;
  dqdc[25] = dcdc_fac;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 5] += dqdc[k];
    J[27 * k + 14] += dqdc[k];
    J[27 * k + 17] -= dqdc[k];
  }
  J[707] += dqdT; // dwdot[H2O]/dT
  J[716] += dqdT; // dwdot[CH2(S)]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT

  // reaction 106: CH2O + H (+M) <=> CH2OH (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha =
    mixture + sc[1] + (6 - 1) * sc[5] + (1.5 - 1) * sc[9] + sc[10] + sc[11];
  // forward
  phi_f = sc[2] * sc[20];
  k_f = 540000 * exp(0.454 * tc[0] - (1811.57999689699) * invT);
  dlnkfdT = 0.454 * invT + (1811.57999689699) * invT2;
  // pressure-fall-off
  k_0 = 1.27e+32 * exp(-4.82 * tc[0] - (3286.00482770482) * invT);
  Pr = 1e-12 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -4.82 * invT + (3286.00482770482) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.7187) * exp(-T / 103);
  Fcent2 = 0.7187 * exp(-T / 1291);
  Fcent3 = exp(-4160 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT =
    log10e / Fcent * (-Fcent1 / 103 - Fcent2 / 1291 + Fcent3 * 4160 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[19];
  Kc = refCinv * exp(g_RT[2] - g_RT[19] + g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[20]) + (h_RT[19]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[2] -= q;  // H
  wdot[19] += q; // CH2OH
  wdot[20] -= q; // CH2O
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = dcdc_fac;
  dqdc[1] = 2 * dcdc_fac;
  dqdc[2] = dcdc_fac + k_f * sc[20];
  dqdc[3] = dcdc_fac;
  dqdc[4] = dcdc_fac;
  dqdc[5] = 6 * dcdc_fac;
  dqdc[6] = dcdc_fac;
  dqdc[7] = dcdc_fac;
  dqdc[8] = dcdc_fac;
  dqdc[9] = 1.5 * dcdc_fac;
  dqdc[10] = 2 * dcdc_fac;
  dqdc[11] = 2 * dcdc_fac;
  dqdc[12] = dcdc_fac;
  dqdc[13] = dcdc_fac;
  dqdc[14] = dcdc_fac;
  dqdc[15] = dcdc_fac;
  dqdc[16] = dcdc_fac;
  dqdc[17] = dcdc_fac;
  dqdc[18] = dcdc_fac;
  dqdc[19] = dcdc_fac - k_r;
  dqdc[20] = dcdc_fac + k_f * sc[2];
  dqdc[21] = dcdc_fac;
  dqdc[22] = dcdc_fac;
  dqdc[23] = dcdc_fac;
  dqdc[24] = dcdc_fac;
  dqdc[25] = dcdc_fac;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 2] -= dqdc[k];
    J[27 * k + 19] += dqdc[k];
    J[27 * k + 20] -= dqdc[k];
  }
  J[704] -= dqdT; // dwdot[H]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT

  // reaction 107: CH3O (+M) <=> CH2O + H (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha =
    mixture + sc[1] + (6 - 1) * sc[5] + (1.5 - 1) * sc[9] + sc[10] + sc[11];
  // forward
  phi_f = sc[18];
  k_f = 68000000000000 * exp(-(13169.1801441095) * invT);
  dlnkfdT = (13169.1801441095) * invT2;
  // pressure-fall-off
  k_0 = 1.87e+25 * exp(-3 * tc[0] - (12231.6874957153) * invT);
  Pr = 1e-6 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -3 * invT + (12231.6874957153) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.9) * exp(-T / 2500);
  Fcent2 = 0.9 * exp(-T / 1300);
  Fcent3 = exp(-1e+99 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT =
    log10e / Fcent * (-Fcent1 / 2500 - Fcent2 / 1300 + Fcent3 * 1e+99 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[2] * sc[20];
  Kc = refC * exp(-g_RT[2] + g_RT[18] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[18]) + (h_RT[2] + h_RT[20]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[2] += q;  // H
  wdot[18] -= q; // CH3O
  wdot[20] += q; // CH2O
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = dcdc_fac;
  dqdc[1] = 2 * dcdc_fac;
  dqdc[2] = dcdc_fac - k_r * sc[20];
  dqdc[3] = dcdc_fac;
  dqdc[4] = dcdc_fac;
  dqdc[5] = 6 * dcdc_fac;
  dqdc[6] = dcdc_fac;
  dqdc[7] = dcdc_fac;
  dqdc[8] = dcdc_fac;
  dqdc[9] = 1.5 * dcdc_fac;
  dqdc[10] = 2 * dcdc_fac;
  dqdc[11] = 2 * dcdc_fac;
  dqdc[12] = dcdc_fac;
  dqdc[13] = dcdc_fac;
  dqdc[14] = dcdc_fac;
  dqdc[15] = dcdc_fac;
  dqdc[16] = dcdc_fac;
  dqdc[17] = dcdc_fac;
  dqdc[18] = dcdc_fac + k_f;
  dqdc[19] = dcdc_fac;
  dqdc[20] = dcdc_fac - k_r * sc[2];
  dqdc[21] = dcdc_fac;
  dqdc[22] = dcdc_fac;
  dqdc[23] = dcdc_fac;
  dqdc[24] = dcdc_fac;
  dqdc[25] = dcdc_fac;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 2] += dqdc[k];
    J[27 * k + 18] -= dqdc[k];
    J[27 * k + 20] += dqdc[k];
  }
  J[704] += dqdT; // dwdot[H]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 20: CO + O (+M) <=> CO2 (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + sc[1] + (12 - 1) * sc[5] + (1.75 - 1) * sc[9] +
          (3.6 - 1) * sc[10];
  // forward
  phi_f = sc[4] * sc[9];
  k_f = 13620 * exp(-(1199.66853127845) * invT);
  dlnkfdT = (1199.66853127845) * invT2;
  // pressure-fall-off
  k_0 = 1.17e+24 * exp(-2.79 * tc[0] - (2108.98104638758) * invT);
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
  phi_r = sc[10];
  Kc = refCinv * exp(g_RT[4] + g_RT[9] - g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[10]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[4] -= q;  // O
  wdot[9] -= q;  // CO
  wdot[10] += q; // CO2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = dcdc_fac;
  dqdc[1] = 2 * dcdc_fac;
  dqdc[2] = dcdc_fac;
  dqdc[3] = dcdc_fac;
  dqdc[4] = dcdc_fac + k_f * sc[9];
  dqdc[5] = 12 * dcdc_fac;
  dqdc[6] = dcdc_fac;
  dqdc[7] = dcdc_fac;
  dqdc[8] = dcdc_fac;
  dqdc[9] = 1.75 * dcdc_fac + k_f * sc[4];
  dqdc[10] = 3.6 * dcdc_fac - k_r;
  dqdc[11] = dcdc_fac;
  dqdc[12] = dcdc_fac;
  dqdc[13] = dcdc_fac;
  dqdc[14] = dcdc_fac;
  dqdc[15] = dcdc_fac;
  dqdc[16] = dcdc_fac;
  dqdc[17] = dcdc_fac;
  dqdc[18] = dcdc_fac;
  dqdc[19] = dcdc_fac;
  dqdc[20] = dcdc_fac;
  dqdc[21] = dcdc_fac;
  dqdc[22] = dcdc_fac;
  dqdc[23] = dcdc_fac;
  dqdc[24] = dcdc_fac;
  dqdc[25] = dcdc_fac;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 4] -= dqdc[k];
    J[27 * k + 9] -= dqdc[k];
    J[27 * k + 10] += dqdc[k];
  }
  J[706] -= dqdT; // dwdot[O]/dT
  J[711] -= dqdT; // dwdot[CO]/dT
  J[712] += dqdT; // dwdot[CO2]/dT

  // reaction 3: H + OH + M <=> H2O + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (0.73 - 1) * sc[1] + (3.65 - 1) * sc[5] + sc[11];
  // forward
  phi_f = sc[2] * sc[6];
  k_f = 35000000000 * exp(-2 * tc[0]);
  dlnkfdT = -2 * invT;
  // reverse
  phi_r = sc[5];
  Kc = refCinv * exp(g_RT[2] - g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[5]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // H
  wdot[5] += q; // H2O
  wdot[6] -= q; // OH
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = q_nocor;
  dqdc[1] = 0.73 * q_nocor;
  dqdc[2] = q_nocor + k_f * sc[6];
  dqdc[3] = q_nocor;
  dqdc[4] = q_nocor;
  dqdc[5] = 3.65 * q_nocor - k_r;
  dqdc[6] = q_nocor + k_f * sc[2];
  dqdc[7] = q_nocor;
  dqdc[8] = q_nocor;
  dqdc[9] = q_nocor;
  dqdc[10] = q_nocor;
  dqdc[11] = 2 * q_nocor;
  dqdc[12] = q_nocor;
  dqdc[13] = q_nocor;
  dqdc[14] = q_nocor;
  dqdc[15] = q_nocor;
  dqdc[16] = q_nocor;
  dqdc[17] = q_nocor;
  dqdc[18] = q_nocor;
  dqdc[19] = q_nocor;
  dqdc[20] = q_nocor;
  dqdc[21] = q_nocor;
  dqdc[22] = q_nocor;
  dqdc[23] = q_nocor;
  dqdc[24] = q_nocor;
  dqdc[25] = q_nocor;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 2] -= dqdc[k];
    J[27 * k + 5] += dqdc[k];
    J[27 * k + 6] -= dqdc[k];
  }
  J[704] -= dqdT; // dwdot[H]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT

  // reaction 5: H + O + M <=> OH + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[1] + (12 - 1) * sc[5] + (1.5 - 1) * sc[9] +
          sc[10] + sc[11];
  // forward
  phi_f = sc[2] * sc[4];
  k_f = 4714000 * exp(-1 * tc[0]);
  dlnkfdT = -1 * invT;
  // reverse
  phi_r = sc[6];
  Kc = refCinv * exp(g_RT[2] + g_RT[4] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (h_RT[6]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // H
  wdot[4] -= q; // O
  wdot[6] += q; // OH
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = q_nocor;
  dqdc[1] = 2.5 * q_nocor;
  dqdc[2] = q_nocor + k_f * sc[4];
  dqdc[3] = q_nocor;
  dqdc[4] = q_nocor + k_f * sc[2];
  dqdc[5] = 12 * q_nocor;
  dqdc[6] = q_nocor - k_r;
  dqdc[7] = q_nocor;
  dqdc[8] = q_nocor;
  dqdc[9] = 1.5 * q_nocor;
  dqdc[10] = 2 * q_nocor;
  dqdc[11] = 2 * q_nocor;
  dqdc[12] = q_nocor;
  dqdc[13] = q_nocor;
  dqdc[14] = q_nocor;
  dqdc[15] = q_nocor;
  dqdc[16] = q_nocor;
  dqdc[17] = q_nocor;
  dqdc[18] = q_nocor;
  dqdc[19] = q_nocor;
  dqdc[20] = q_nocor;
  dqdc[21] = q_nocor;
  dqdc[22] = q_nocor;
  dqdc[23] = q_nocor;
  dqdc[24] = q_nocor;
  dqdc[25] = q_nocor;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 2] -= dqdc[k];
    J[27 * k + 4] -= dqdc[k];
    J[27 * k + 6] += dqdc[k];
  }
  J[704] -= dqdT; // dwdot[H]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT

  // reaction 97: HCO + M <=> CO + H + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha =
    mixture + sc[1] + (6 - 1) * sc[5] + (1.5 - 1) * sc[9] + sc[10] + sc[11];
  // forward
  phi_f = sc[21];
  k_f = 570000 * exp(0.66 * tc[0] - (7482.83182051618) * invT);
  dlnkfdT = 0.66 * invT + (7482.83182051618) * invT2;
  // reverse
  phi_r = sc[2] * sc[9];
  Kc = refC * exp(-g_RT[2] - g_RT[9] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[21]) + (h_RT[2] + h_RT[9]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[9] += q;  // CO
  wdot[21] -= q; // HCO
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = q_nocor;
  dqdc[1] = 2 * q_nocor;
  dqdc[2] = q_nocor - k_r * sc[9];
  dqdc[3] = q_nocor;
  dqdc[4] = q_nocor;
  dqdc[5] = 6 * q_nocor;
  dqdc[6] = q_nocor;
  dqdc[7] = q_nocor;
  dqdc[8] = q_nocor;
  dqdc[9] = 1.5 * q_nocor - k_r * sc[2];
  dqdc[10] = 2 * q_nocor;
  dqdc[11] = 2 * q_nocor;
  dqdc[12] = q_nocor;
  dqdc[13] = q_nocor;
  dqdc[14] = q_nocor;
  dqdc[15] = q_nocor;
  dqdc[16] = q_nocor;
  dqdc[17] = q_nocor;
  dqdc[18] = q_nocor;
  dqdc[19] = q_nocor;
  dqdc[20] = q_nocor;
  dqdc[21] = q_nocor + k_f;
  dqdc[22] = q_nocor;
  dqdc[23] = q_nocor;
  dqdc[24] = q_nocor;
  dqdc[25] = q_nocor;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 2] += dqdc[k];
    J[27 * k + 9] += dqdc[k];
    J[27 * k + 21] -= dqdc[k];
  }
  J[704] += dqdT; // dwdot[H]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 108: O + OH + M <=> HO2 + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture;
  // forward
  phi_f = sc[4] * sc[6];
  k_f = 10000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[8];
  Kc = refCinv * exp(g_RT[4] + g_RT[6] - g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[8]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q; // O
  wdot[6] -= q; // OH
  wdot[8] += q; // HO2
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = q_nocor;
  dqdc[1] = q_nocor;
  dqdc[2] = q_nocor;
  dqdc[3] = q_nocor;
  dqdc[4] = q_nocor + k_f * sc[6];
  dqdc[5] = q_nocor;
  dqdc[6] = q_nocor + k_f * sc[4];
  dqdc[7] = q_nocor;
  dqdc[8] = q_nocor - k_r;
  dqdc[9] = q_nocor;
  dqdc[10] = q_nocor;
  dqdc[11] = q_nocor;
  dqdc[12] = q_nocor;
  dqdc[13] = q_nocor;
  dqdc[14] = q_nocor;
  dqdc[15] = q_nocor;
  dqdc[16] = q_nocor;
  dqdc[17] = q_nocor;
  dqdc[18] = q_nocor;
  dqdc[19] = q_nocor;
  dqdc[20] = q_nocor;
  dqdc[21] = q_nocor;
  dqdc[22] = q_nocor;
  dqdc[23] = q_nocor;
  dqdc[24] = q_nocor;
  dqdc[25] = q_nocor;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 4] -= dqdc[k];
    J[27 * k + 6] -= dqdc[k];
    J[27 * k + 8] += dqdc[k];
  }
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[710] += dqdT; // dwdot[HO2]/dT

  // reaction 120: CH3 + O + M => CH3O + M
  // a third-body and non-pressure-fall-off reaction
  // Irreversible reaction in _ajac_reaction may not work
  // 3-body correction factor
  alpha = mixture;
  // forward
  phi_f = sc[4] * sc[12];
  k_f = 50000;
  dlnkfdT = 0.0;
  // rate of progress
  q_nocor = k_f * phi_f;
  q = alpha * q_nocor;
  dqdT = alpha * dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[4] -= q;  // O
  wdot[12] -= q; // CH3
  wdot[18] += q; // CH3O
  // for convenience
  k_f *= alpha;
  k_r = 0.0;
  dqdc[0] = q_nocor;
  dqdc[1] = q_nocor;
  dqdc[2] = q_nocor;
  dqdc[3] = q_nocor;
  dqdc[4] = q_nocor + k_f * sc[12];
  dqdc[5] = q_nocor;
  dqdc[6] = q_nocor;
  dqdc[7] = q_nocor;
  dqdc[8] = q_nocor;
  dqdc[9] = q_nocor;
  dqdc[10] = q_nocor;
  dqdc[11] = q_nocor;
  dqdc[12] = q_nocor + k_f * sc[4];
  dqdc[13] = q_nocor;
  dqdc[14] = q_nocor;
  dqdc[15] = q_nocor;
  dqdc[16] = q_nocor;
  dqdc[17] = q_nocor;
  dqdc[18] = q_nocor;
  dqdc[19] = q_nocor;
  dqdc[20] = q_nocor;
  dqdc[21] = q_nocor;
  dqdc[22] = q_nocor;
  dqdc[23] = q_nocor;
  dqdc[24] = q_nocor;
  dqdc[25] = q_nocor;
  for (int k = 0; k < 26; k++) {
    J[27 * k + 4] -= dqdc[k];
    J[27 * k + 12] -= dqdc[k];
    J[27 * k + 18] += dqdc[k];
  }
  J[706] -= dqdT; // dwdot[O]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 131: E + O2 + O2 <=> O2- + O2
  // a third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[3] * sc[3]) * sc[25];
  k_f = 1520000000 * exp(-1 * tc[0] - (600.156324305341) * invT);
  dlnkfdT = -1 * invT + (600.156324305341) * invT2;
  // reverse
  phi_r = sc[3] * sc[24];
  Kc = refCinv * exp(g_RT[3] - g_RT[24] + g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT =
    invT * (-(2.000000 * h_RT[3] + h_RT[25]) + (h_RT[3] + h_RT[24]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[24] += q; // O2-
  wdot[25] -= q; // E
  // d()/d[O2]
  dqdci = +k_f * 2.000000 * sc[3] * sc[25] - k_r * sc[24];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[105] += dqdci; // dwdot[O2-]/d[O2]
  J[106] -= dqdci; // dwdot[E]/d[O2]
  // d()/d[O2-]
  dqdci = -k_r * sc[3];
  J[651] -= dqdci; // dwdot[O2]/d[O2-]
  J[672] += dqdci; // dwdot[O2-]/d[O2-]
  J[673] -= dqdci; // dwdot[E]/d[O2-]
  // d()/d[E]
  dqdci = +k_f * sc[3] * sc[3];
  J[678] -= dqdci; // dwdot[O2]/d[E]
  J[699] += dqdci; // dwdot[O2-]/d[E]
  J[700] -= dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[726] += dqdT; // dwdot[O2-]/dT
  J[727] -= dqdT; // dwdot[E]/dT

  // reaction 132: E + O2 + H2O <=> O2- + H2O
  // a third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[5] * sc[25];
  k_f = 5080000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[5] * sc[24];
  Kc = refCinv * exp(g_RT[3] - g_RT[24] + g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT =
    invT * (-(h_RT[3] + h_RT[5] + h_RT[25]) + (h_RT[5] + h_RT[24]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[24] += q; // O2-
  wdot[25] -= q; // E
  // d()/d[O2]
  dqdci = +k_f * sc[5] * sc[25];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[105] += dqdci; // dwdot[O2-]/d[O2]
  J[106] -= dqdci; // dwdot[E]/d[O2]
  // d()/d[H2O]
  dqdci = +k_f * sc[3] * sc[25] - k_r * sc[24];
  J[138] -= dqdci; // dwdot[O2]/d[H2O]
  J[159] += dqdci; // dwdot[O2-]/d[H2O]
  J[160] -= dqdci; // dwdot[E]/d[H2O]
  // d()/d[O2-]
  dqdci = -k_r * sc[5];
  J[651] -= dqdci; // dwdot[O2]/d[O2-]
  J[672] += dqdci; // dwdot[O2-]/d[O2-]
  J[673] -= dqdci; // dwdot[E]/d[O2-]
  // d()/d[E]
  dqdci = +k_f * sc[3] * sc[5];
  J[678] -= dqdci; // dwdot[O2]/d[E]
  J[699] += dqdci; // dwdot[O2-]/d[E]
  J[700] -= dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[726] += dqdT; // dwdot[O2-]/dT
  J[727] -= dqdT; // dwdot[E]/dT

  // reaction 133: E + O2 + N2 <=> O2- + N2
  // a third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[3] * sc[25];
  k_f = 3590000000 * exp(-2 * tc[0] - (69.7558942138502) * invT);
  dlnkfdT = -2 * invT + (69.7558942138502) * invT2;
  // reverse
  phi_r = sc[0] * sc[24];
  Kc = refCinv * exp(g_RT[3] - g_RT[24] + g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT =
    invT * (-(h_RT[0] + h_RT[3] + h_RT[25]) + (h_RT[0] + h_RT[24]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[24] += q; // O2-
  wdot[25] -= q; // E
  // d()/d[N2]
  dqdci = +k_f * sc[3] * sc[25] - k_r * sc[24];
  J[3] -= dqdci;  // dwdot[O2]/d[N2]
  J[24] += dqdci; // dwdot[O2-]/d[N2]
  J[25] -= dqdci; // dwdot[E]/d[N2]
  // d()/d[O2]
  dqdci = +k_f * sc[0] * sc[25];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[105] += dqdci; // dwdot[O2-]/d[O2]
  J[106] -= dqdci; // dwdot[E]/d[O2]
  // d()/d[O2-]
  dqdci = -k_r * sc[0];
  J[651] -= dqdci; // dwdot[O2]/d[O2-]
  J[672] += dqdci; // dwdot[O2-]/d[O2-]
  J[673] -= dqdci; // dwdot[E]/d[O2-]
  // d()/d[E]
  dqdci = +k_f * sc[0] * sc[3];
  J[678] -= dqdci; // dwdot[O2]/d[E]
  J[699] += dqdci; // dwdot[O2-]/d[E]
  J[700] -= dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[726] += dqdT; // dwdot[O2-]/dT
  J[727] -= dqdT; // dwdot[E]/dT

  // reaction 0: H2 + O <=> H + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[4];
  k_f = 0.0508 * exp(2.67 * tc[0] - (3166.2392612433) * invT);
  dlnkfdT = 2.67 * invT + (3166.2392612433) * invT2;
  // reverse
  phi_r = sc[2] * sc[6];
  Kc = exp(g_RT[1] - g_RT[2] + g_RT[4] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[4]) + (h_RT[2] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q; // H2
  wdot[2] += q; // H
  wdot[4] -= q; // O
  wdot[6] += q; // OH
  // d()/d[H2]
  dqdci = +k_f * sc[4];
  J[28] -= dqdci; // dwdot[H2]/d[H2]
  J[29] += dqdci; // dwdot[H]/d[H2]
  J[31] -= dqdci; // dwdot[O]/d[H2]
  J[33] += dqdci; // dwdot[OH]/d[H2]
  // d()/d[H]
  dqdci = -k_r * sc[6];
  J[55] -= dqdci; // dwdot[H2]/d[H]
  J[56] += dqdci; // dwdot[H]/d[H]
  J[58] -= dqdci; // dwdot[O]/d[H]
  J[60] += dqdci; // dwdot[OH]/d[H]
  // d()/d[O]
  dqdci = +k_f * sc[1];
  J[109] -= dqdci; // dwdot[H2]/d[O]
  J[110] += dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[2];
  J[163] -= dqdci; // dwdot[H2]/d[OH]
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  // d()/dT
  J[703] -= dqdT; // dwdot[H2]/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT

  // reaction 1: H2 + OH <=> H + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[6];
  k_f = 43800000 * exp(-(3517.48449397499) * invT);
  dlnkfdT = (3517.48449397499) * invT2;
  // reverse
  phi_r = sc[2] * sc[5];
  Kc = exp(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[2] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q; // H2
  wdot[2] += q; // H
  wdot[5] += q; // H2O
  wdot[6] -= q; // OH
  // d()/d[H2]
  dqdci = +k_f * sc[6];
  J[28] -= dqdci; // dwdot[H2]/d[H2]
  J[29] += dqdci; // dwdot[H]/d[H2]
  J[32] += dqdci; // dwdot[H2O]/d[H2]
  J[33] -= dqdci; // dwdot[OH]/d[H2]
  // d()/d[H]
  dqdci = -k_r * sc[5];
  J[55] -= dqdci; // dwdot[H2]/d[H]
  J[56] += dqdci; // dwdot[H]/d[H]
  J[59] += dqdci; // dwdot[H2O]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  // d()/d[H2O]
  dqdci = -k_r * sc[2];
  J[136] -= dqdci; // dwdot[H2]/d[H2O]
  J[137] += dqdci; // dwdot[H]/d[H2O]
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[1];
  J[163] -= dqdci; // dwdot[H2]/d[OH]
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  // d()/dT
  J[703] -= dqdT; // dwdot[H2]/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT

  // reaction 2: H + O2 <=> O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[3];
  k_f = 114000000 * exp(-(7692.16995349094) * invT);
  dlnkfdT = (7692.16995349094) * invT2;
  // reverse
  phi_r = sc[4] * sc[6];
  Kc = exp(g_RT[2] + g_RT[3] - g_RT[4] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[4] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // H
  wdot[3] -= q; // O2
  wdot[4] += q; // O
  wdot[6] += q; // OH
  // d()/d[H]
  dqdci = +k_f * sc[3];
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[57] -= dqdci; // dwdot[O2]/d[H]
  J[58] += dqdci; // dwdot[O]/d[H]
  J[60] += dqdci; // dwdot[OH]/d[H]
  // d()/d[O2]
  dqdci = +k_f * sc[2];
  J[83] -= dqdci; // dwdot[H]/d[O2]
  J[84] -= dqdci; // dwdot[O2]/d[O2]
  J[85] += dqdci; // dwdot[O]/d[O2]
  J[87] += dqdci; // dwdot[OH]/d[O2]
  // d()/d[O]
  dqdci = -k_r * sc[6];
  J[110] -= dqdci; // dwdot[H]/d[O]
  J[111] -= dqdci; // dwdot[O2]/d[O]
  J[112] += dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[4];
  J[164] -= dqdci; // dwdot[H]/d[OH]
  J[165] -= dqdci; // dwdot[O2]/d[OH]
  J[166] += dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  // d()/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[706] += dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT

  // reaction 4: H2O + O <=> 2 OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[5];
  k_f = 67 * exp(1.704 * tc[0] - (7541.60752708217) * invT);
  dlnkfdT = 1.704 * invT + (7541.60752708217) * invT2;
  // reverse
  phi_r = (sc[6] * sc[6]);
  Kc = exp(g_RT[4] + g_RT[5] - 2.000000 * g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[5]) + (2.000000 * h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;     // O
  wdot[5] -= q;     // H2O
  wdot[6] += 2 * q; // OH
  // d()/d[O]
  dqdci = +k_f * sc[5];
  J[112] -= dqdci;     // dwdot[O]/d[O]
  J[113] -= dqdci;     // dwdot[H2O]/d[O]
  J[114] += 2 * dqdci; // dwdot[OH]/d[O]
  // d()/d[H2O]
  dqdci = +k_f * sc[4];
  J[139] -= dqdci;     // dwdot[O]/d[H2O]
  J[140] -= dqdci;     // dwdot[H2O]/d[H2O]
  J[141] += 2 * dqdci; // dwdot[OH]/d[H2O]
  // d()/d[OH]
  dqdci = -k_r * 2.000000 * sc[6];
  J[166] -= dqdci;     // dwdot[O]/d[OH]
  J[167] -= dqdci;     // dwdot[H2O]/d[OH]
  J[168] += 2 * dqdci; // dwdot[OH]/d[OH]
  // d()/dT
  J[706] -= dqdT;     // dwdot[O]/dT
  J[707] -= dqdT;     // dwdot[H2O]/dT
  J[708] += 2 * dqdT; // dwdot[OH]/dT

  // reaction 7: H + H2O2 <=> H2O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[7];
  k_f = 24100000 * exp(-(1997.77016324474) * invT);
  dlnkfdT = (1997.77016324474) * invT2;
  // reverse
  phi_r = sc[5] * sc[6];
  Kc = exp(g_RT[2] - g_RT[5] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[5] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // H
  wdot[5] += q; // H2O
  wdot[6] += q; // OH
  wdot[7] -= q; // H2O2
  // d()/d[H]
  dqdci = +k_f * sc[7];
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[59] += dqdci; // dwdot[H2O]/d[H]
  J[60] += dqdci; // dwdot[OH]/d[H]
  J[61] -= dqdci; // dwdot[H2O2]/d[H]
  // d()/d[H2O]
  dqdci = -k_r * sc[6];
  J[137] -= dqdci; // dwdot[H]/d[H2O]
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] += dqdci; // dwdot[OH]/d[H2O]
  J[142] -= dqdci; // dwdot[H2O2]/d[H2O]
  // d()/d[OH]
  dqdci = -k_r * sc[5];
  J[164] -= dqdci; // dwdot[H]/d[OH]
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[169] -= dqdci; // dwdot[H2O2]/d[OH]
  // d()/d[H2O2]
  dqdci = +k_f * sc[2];
  J[191] -= dqdci; // dwdot[H]/d[H2O2]
  J[194] += dqdci; // dwdot[H2O]/d[H2O2]
  J[195] += dqdci; // dwdot[OH]/d[H2O2]
  J[196] -= dqdci; // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[709] -= dqdT; // dwdot[H2O2]/dT

  // reaction 8: H + H2O2 <=> H2 + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[7];
  k_f = 21500 * exp(1 * tc[0] - (3019.29999482832) * invT);
  dlnkfdT = 1 * invT + (3019.29999482832) * invT2;
  // reverse
  phi_r = sc[1] * sc[8];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[1] + h_RT[8]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q; // H2
  wdot[2] -= q; // H
  wdot[7] -= q; // H2O2
  wdot[8] += q; // HO2
  // d()/d[H2]
  dqdci = -k_r * sc[8];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[34] -= dqdci; // dwdot[H2O2]/d[H2]
  J[35] += dqdci; // dwdot[HO2]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[7];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[61] -= dqdci; // dwdot[H2O2]/d[H]
  J[62] += dqdci; // dwdot[HO2]/d[H]
  // d()/d[H2O2]
  dqdci = +k_f * sc[2];
  J[190] += dqdci; // dwdot[H2]/d[H2O2]
  J[191] -= dqdci; // dwdot[H]/d[H2O2]
  J[196] -= dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] += dqdci; // dwdot[HO2]/d[H2O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[1];
  J[217] += dqdci; // dwdot[H2]/d[HO2]
  J[218] -= dqdci; // dwdot[H]/d[HO2]
  J[223] -= dqdci; // dwdot[H2O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[709] -= dqdT; // dwdot[H2O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT

  // reaction 9: H2O2 + O <=> HO2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[7];
  k_f = 9.55 * exp(2 * tc[0] - (1997.77016324474) * invT);
  dlnkfdT = 2 * invT + (1997.77016324474) * invT2;
  // reverse
  phi_r = sc[6] * sc[8];
  Kc = exp(g_RT[4] - g_RT[6] + g_RT[7] - g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[6] + h_RT[8]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q; // O
  wdot[6] += q; // OH
  wdot[7] -= q; // H2O2
  wdot[8] += q; // HO2
  // d()/d[O]
  dqdci = +k_f * sc[7];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[115] -= dqdci; // dwdot[H2O2]/d[O]
  J[116] += dqdci; // dwdot[HO2]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[8];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[169] -= dqdci; // dwdot[H2O2]/d[OH]
  J[170] += dqdci; // dwdot[HO2]/d[OH]
  // d()/d[H2O2]
  dqdci = +k_f * sc[4];
  J[193] -= dqdci; // dwdot[O]/d[H2O2]
  J[195] += dqdci; // dwdot[OH]/d[H2O2]
  J[196] -= dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] += dqdci; // dwdot[HO2]/d[H2O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[6];
  J[220] -= dqdci; // dwdot[O]/d[HO2]
  J[222] += dqdci; // dwdot[OH]/d[HO2]
  J[223] -= dqdci; // dwdot[H2O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[709] -= dqdT; // dwdot[H2O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT

  // reaction 10: H2O2 + OH <=> H2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[7];
  k_f = 1740000 * exp(-(160.022899725901) * invT);
  dlnkfdT = (160.022899725901) * invT2;
  // reverse
  phi_r = sc[5] * sc[8];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[7] - g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[7]) + (h_RT[5] + h_RT[8]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q; // H2O
  wdot[6] -= q; // OH
  wdot[7] -= q; // H2O2
  wdot[8] += q; // HO2
  // d()/d[H2O]
  dqdci = -k_r * sc[8];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[142] -= dqdci; // dwdot[H2O2]/d[H2O]
  J[143] += dqdci; // dwdot[HO2]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[7];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[169] -= dqdci; // dwdot[H2O2]/d[OH]
  J[170] += dqdci; // dwdot[HO2]/d[OH]
  // d()/d[H2O2]
  dqdci = +k_f * sc[6];
  J[194] += dqdci; // dwdot[H2O]/d[H2O2]
  J[195] -= dqdci; // dwdot[OH]/d[H2O2]
  J[196] -= dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] += dqdci; // dwdot[HO2]/d[H2O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[5];
  J[221] += dqdci; // dwdot[H2O]/d[HO2]
  J[222] -= dqdci; // dwdot[OH]/d[HO2]
  J[223] -= dqdci; // dwdot[H2O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[709] -= dqdT; // dwdot[H2O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT

  // reaction 11: H2O2 + OH <=> H2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[7];
  k_f = 75900000 * exp(-(3657.88194373451) * invT);
  dlnkfdT = (3657.88194373451) * invT2;
  // reverse
  phi_r = sc[5] * sc[8];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[7] - g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[7]) + (h_RT[5] + h_RT[8]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q; // H2O
  wdot[6] -= q; // OH
  wdot[7] -= q; // H2O2
  wdot[8] += q; // HO2
  // d()/d[H2O]
  dqdci = -k_r * sc[8];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[142] -= dqdci; // dwdot[H2O2]/d[H2O]
  J[143] += dqdci; // dwdot[HO2]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[7];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[169] -= dqdci; // dwdot[H2O2]/d[OH]
  J[170] += dqdci; // dwdot[HO2]/d[OH]
  // d()/d[H2O2]
  dqdci = +k_f * sc[6];
  J[194] += dqdci; // dwdot[H2O]/d[H2O2]
  J[195] -= dqdci; // dwdot[OH]/d[H2O2]
  J[196] -= dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] += dqdci; // dwdot[HO2]/d[H2O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[5];
  J[221] += dqdci; // dwdot[H2O]/d[HO2]
  J[222] -= dqdci; // dwdot[OH]/d[HO2]
  J[223] -= dqdci; // dwdot[H2O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[709] -= dqdT; // dwdot[H2O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT

  // reaction 12: H + HO2 <=> 2 OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[8];
  k_f = 70790000 * exp(-(148.448916412392) * invT);
  dlnkfdT = (148.448916412392) * invT2;
  // reverse
  phi_r = (sc[6] * sc[6]);
  Kc = exp(g_RT[2] - 2.000000 * g_RT[6] + g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (2.000000 * h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q;     // H
  wdot[6] += 2 * q; // OH
  wdot[8] -= q;     // HO2
  // d()/d[H]
  dqdci = +k_f * sc[8];
  J[56] -= dqdci;     // dwdot[H]/d[H]
  J[60] += 2 * dqdci; // dwdot[OH]/d[H]
  J[62] -= dqdci;     // dwdot[HO2]/d[H]
  // d()/d[OH]
  dqdci = -k_r * 2.000000 * sc[6];
  J[164] -= dqdci;     // dwdot[H]/d[OH]
  J[168] += 2 * dqdci; // dwdot[OH]/d[OH]
  J[170] -= dqdci;     // dwdot[HO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[2];
  J[218] -= dqdci;     // dwdot[H]/d[HO2]
  J[222] += 2 * dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci;     // dwdot[HO2]/d[HO2]
  // d()/dT
  J[704] -= dqdT;     // dwdot[H]/dT
  J[708] += 2 * dqdT; // dwdot[OH]/dT
  J[710] -= dqdT;     // dwdot[HO2]/dT

  // reaction 13: H + HO2 <=> H2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[8];
  k_f = 11402 * exp(1.083 * tc[0] - (278.671325189338) * invT);
  dlnkfdT = 1.083 * invT + (278.671325189338) * invT2;
  // reverse
  phi_r = sc[1] * sc[3];
  Kc = exp(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[1] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q; // H2
  wdot[2] -= q; // H
  wdot[3] += q; // O2
  wdot[8] -= q; // HO2
  // d()/d[H2]
  dqdci = -k_r * sc[3];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[30] += dqdci; // dwdot[O2]/d[H2]
  J[35] -= dqdci; // dwdot[HO2]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[8];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[57] += dqdci; // dwdot[O2]/d[H]
  J[62] -= dqdci; // dwdot[HO2]/d[H]
  // d()/d[O2]
  dqdci = -k_r * sc[1];
  J[82] += dqdci; // dwdot[H2]/d[O2]
  J[83] -= dqdci; // dwdot[H]/d[O2]
  J[84] += dqdci; // dwdot[O2]/d[O2]
  J[89] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[2];
  J[217] += dqdci; // dwdot[H2]/d[HO2]
  J[218] -= dqdci; // dwdot[H]/d[HO2]
  J[219] += dqdci; // dwdot[O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[705] += dqdT; // dwdot[O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT

  // reaction 14: HO2 + O <=> O2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[8];
  k_f = 32500000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[3] * sc[6];
  Kc = exp(-g_RT[3] + g_RT[4] - g_RT[6] + g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[8]) + (h_RT[3] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q; // O2
  wdot[4] -= q; // O
  wdot[6] += q; // OH
  wdot[8] -= q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[6];
  J[84] += dqdci; // dwdot[O2]/d[O2]
  J[85] -= dqdci; // dwdot[O]/d[O2]
  J[87] += dqdci; // dwdot[OH]/d[O2]
  J[89] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[O]
  dqdci = +k_f * sc[8];
  J[111] += dqdci; // dwdot[O2]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[116] -= dqdci; // dwdot[HO2]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[3];
  J[165] += dqdci; // dwdot[O2]/d[OH]
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[170] -= dqdci; // dwdot[HO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[4];
  J[219] += dqdci; // dwdot[O2]/d[HO2]
  J[220] -= dqdci; // dwdot[O]/d[HO2]
  J[222] += dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[705] += dqdT; // dwdot[O2]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT

  // reaction 15: HO2 + OH <=> H2O + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[8];
  k_f = 7000000 * exp(-(-549.995687057926) * invT);
  dlnkfdT = (-549.995687057926) * invT2;
  // reverse
  phi_r = sc[3] * sc[5];
  Kc = exp(-g_RT[3] - g_RT[5] + g_RT[6] + g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[3] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q; // O2
  wdot[5] += q; // H2O
  wdot[6] -= q; // OH
  wdot[8] -= q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[5];
  J[84] += dqdci; // dwdot[O2]/d[O2]
  J[86] += dqdci; // dwdot[H2O]/d[O2]
  J[87] -= dqdci; // dwdot[OH]/d[O2]
  J[89] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[H2O]
  dqdci = -k_r * sc[3];
  J[138] += dqdci; // dwdot[O2]/d[H2O]
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[143] -= dqdci; // dwdot[HO2]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[8];
  J[165] += dqdci; // dwdot[O2]/d[OH]
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[170] -= dqdci; // dwdot[HO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[6];
  J[219] += dqdci; // dwdot[O2]/d[HO2]
  J[221] += dqdci; // dwdot[H2O]/d[HO2]
  J[222] -= dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[705] += dqdT; // dwdot[O2]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT

  // reaction 16: HO2 + OH <=> H2O + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[8];
  k_f = 450000000 * exp(-(5499.95687057926) * invT);
  dlnkfdT = (5499.95687057926) * invT2;
  // reverse
  phi_r = sc[3] * sc[5];
  Kc = exp(-g_RT[3] - g_RT[5] + g_RT[6] + g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[3] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q; // O2
  wdot[5] += q; // H2O
  wdot[6] -= q; // OH
  wdot[8] -= q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[5];
  J[84] += dqdci; // dwdot[O2]/d[O2]
  J[86] += dqdci; // dwdot[H2O]/d[O2]
  J[87] -= dqdci; // dwdot[OH]/d[O2]
  J[89] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[H2O]
  dqdci = -k_r * sc[3];
  J[138] += dqdci; // dwdot[O2]/d[H2O]
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[143] -= dqdci; // dwdot[HO2]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[8];
  J[165] += dqdci; // dwdot[O2]/d[OH]
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[170] -= dqdci; // dwdot[HO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[6];
  J[219] += dqdci; // dwdot[O2]/d[HO2]
  J[221] += dqdci; // dwdot[H2O]/d[HO2]
  J[222] -= dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[705] += dqdT; // dwdot[O2]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT

  // reaction 17: 2 HO2 <=> H2O2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[8] * sc[8]);
  k_f = 100000000 * exp(-(5555.95482115001) * invT);
  dlnkfdT = (5555.95482115001) * invT2;
  // reverse
  phi_r = sc[3] * sc[7];
  Kc = exp(-g_RT[3] - g_RT[7] + 2.000000 * g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[8]) + (h_RT[3] + h_RT[7]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q;     // O2
  wdot[7] += q;     // H2O2
  wdot[8] -= 2 * q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[7];
  J[84] += dqdci;      // dwdot[O2]/d[O2]
  J[88] += dqdci;      // dwdot[H2O2]/d[O2]
  J[89] += -2 * dqdci; // dwdot[HO2]/d[O2]
  // d()/d[H2O2]
  dqdci = -k_r * sc[3];
  J[192] += dqdci;      // dwdot[O2]/d[H2O2]
  J[196] += dqdci;      // dwdot[H2O2]/d[H2O2]
  J[197] += -2 * dqdci; // dwdot[HO2]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * 2.000000 * sc[8];
  J[219] += dqdci;      // dwdot[O2]/d[HO2]
  J[223] += dqdci;      // dwdot[H2O2]/d[HO2]
  J[224] += -2 * dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[705] += dqdT;      // dwdot[O2]/dT
  J[709] += dqdT;      // dwdot[H2O2]/dT
  J[710] += -2 * dqdT; // dwdot[HO2]/dT

  // reaction 18: 2 HO2 <=> H2O2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[8] * sc[8]);
  k_f = 190000 * exp(-(-708.992024785586) * invT);
  dlnkfdT = (-708.992024785586) * invT2;
  // reverse
  phi_r = sc[3] * sc[7];
  Kc = exp(-g_RT[3] - g_RT[7] + 2.000000 * g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[8]) + (h_RT[3] + h_RT[7]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q;     // O2
  wdot[7] += q;     // H2O2
  wdot[8] -= 2 * q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[7];
  J[84] += dqdci;      // dwdot[O2]/d[O2]
  J[88] += dqdci;      // dwdot[H2O2]/d[O2]
  J[89] += -2 * dqdci; // dwdot[HO2]/d[O2]
  // d()/d[H2O2]
  dqdci = -k_r * sc[3];
  J[192] += dqdci;      // dwdot[O2]/d[H2O2]
  J[196] += dqdci;      // dwdot[H2O2]/d[H2O2]
  J[197] += -2 * dqdci; // dwdot[HO2]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * 2.000000 * sc[8];
  J[219] += dqdci;      // dwdot[O2]/d[HO2]
  J[223] += dqdci;      // dwdot[H2O2]/d[HO2]
  J[224] += -2 * dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[705] += dqdT;      // dwdot[O2]/dT
  J[709] += dqdT;      // dwdot[H2O2]/dT
  J[710] += -2 * dqdT; // dwdot[HO2]/dT

  // reaction 21: CO + OH <=> CO2 + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[9];
  k_f = 0.07015 * exp(2.053 * tc[0] - (-178.994168026739) * invT);
  dlnkfdT = 2.053 * invT + (-178.994168026739) * invT2;
  // reverse
  phi_r = sc[2] * sc[10];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[9] - g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[2] + h_RT[10]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[9] -= q;  // CO
  wdot[10] += q; // CO2
  // d()/d[H]
  dqdci = -k_r * sc[10];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[63] -= dqdci; // dwdot[CO]/d[H]
  J[64] += dqdci; // dwdot[CO2]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[9];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[171] -= dqdci; // dwdot[CO]/d[OH]
  J[172] += dqdci; // dwdot[CO2]/d[OH]
  // d()/d[CO]
  dqdci = +k_f * sc[6];
  J[245] += dqdci; // dwdot[H]/d[CO]
  J[249] -= dqdci; // dwdot[OH]/d[CO]
  J[252] -= dqdci; // dwdot[CO]/d[CO]
  J[253] += dqdci; // dwdot[CO2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[2];
  J[272] += dqdci; // dwdot[H]/d[CO2]
  J[276] -= dqdci; // dwdot[OH]/d[CO2]
  J[279] -= dqdci; // dwdot[CO]/d[CO2]
  J[280] += dqdci; // dwdot[CO2]/d[CO2]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[711] -= dqdT; // dwdot[CO]/dT
  J[712] += dqdT; // dwdot[CO2]/dT

  // reaction 22: CO + OH <=> CO2 + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[9];
  k_f = 5757000 * exp(-0.664 * tc[0] - (166.967289714006) * invT);
  dlnkfdT = -0.664 * invT + (166.967289714006) * invT2;
  // reverse
  phi_r = sc[2] * sc[10];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[9] - g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[2] + h_RT[10]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[9] -= q;  // CO
  wdot[10] += q; // CO2
  // d()/d[H]
  dqdci = -k_r * sc[10];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[63] -= dqdci; // dwdot[CO]/d[H]
  J[64] += dqdci; // dwdot[CO2]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[9];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[171] -= dqdci; // dwdot[CO]/d[OH]
  J[172] += dqdci; // dwdot[CO2]/d[OH]
  // d()/d[CO]
  dqdci = +k_f * sc[6];
  J[245] += dqdci; // dwdot[H]/d[CO]
  J[249] -= dqdci; // dwdot[OH]/d[CO]
  J[252] -= dqdci; // dwdot[CO]/d[CO]
  J[253] += dqdci; // dwdot[CO2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[2];
  J[272] += dqdci; // dwdot[H]/d[CO2]
  J[276] -= dqdci; // dwdot[OH]/d[CO2]
  J[279] -= dqdci; // dwdot[CO]/d[CO2]
  J[280] += dqdci; // dwdot[CO2]/d[CO2]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[711] -= dqdT; // dwdot[CO]/dT
  J[712] += dqdT; // dwdot[CO2]/dT

  // reaction 23: CO + HO2 <=> CO2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[9];
  k_f = 0.157 * exp(2.18 * tc[0] - (9027.70698453667) * invT);
  dlnkfdT = 2.18 * invT + (9027.70698453667) * invT2;
  // reverse
  phi_r = sc[6] * sc[10];
  Kc = exp(-g_RT[6] + g_RT[8] + g_RT[9] - g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[9]) + (h_RT[6] + h_RT[10]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[6] += q;  // OH
  wdot[8] -= q;  // HO2
  wdot[9] -= q;  // CO
  wdot[10] += q; // CO2
  // d()/d[OH]
  dqdci = -k_r * sc[10];
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[170] -= dqdci; // dwdot[HO2]/d[OH]
  J[171] -= dqdci; // dwdot[CO]/d[OH]
  J[172] += dqdci; // dwdot[CO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[9];
  J[222] += dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[225] -= dqdci; // dwdot[CO]/d[HO2]
  J[226] += dqdci; // dwdot[CO2]/d[HO2]
  // d()/d[CO]
  dqdci = +k_f * sc[8];
  J[249] += dqdci; // dwdot[OH]/d[CO]
  J[251] -= dqdci; // dwdot[HO2]/d[CO]
  J[252] -= dqdci; // dwdot[CO]/d[CO]
  J[253] += dqdci; // dwdot[CO2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[6];
  J[276] += dqdci; // dwdot[OH]/d[CO2]
  J[278] -= dqdci; // dwdot[HO2]/d[CO2]
  J[279] -= dqdci; // dwdot[CO]/d[CO2]
  J[280] += dqdci; // dwdot[CO2]/d[CO2]
  // d()/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[711] -= dqdT; // dwdot[CO]/dT
  J[712] += dqdT; // dwdot[CO2]/dT

  // reaction 25: CH4 + H <=> CH3 + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[11];
  k_f = 0.614 * exp(2.5 * tc[0] - (4824.33817506985) * invT);
  dlnkfdT = 2.5 * invT + (4824.33817506985) * invT2;
  // reverse
  phi_r = sc[1] * sc[12];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[11] - g_RT[12]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[11]) + (h_RT[1] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[11] -= q; // CH4
  wdot[12] += q; // CH3
  // d()/d[H2]
  dqdci = -k_r * sc[12];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[38] -= dqdci; // dwdot[CH4]/d[H2]
  J[39] += dqdci; // dwdot[CH3]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[11];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[65] -= dqdci; // dwdot[CH4]/d[H]
  J[66] += dqdci; // dwdot[CH3]/d[H]
  // d()/d[CH4]
  dqdci = +k_f * sc[2];
  J[298] += dqdci; // dwdot[H2]/d[CH4]
  J[299] -= dqdci; // dwdot[H]/d[CH4]
  J[308] -= dqdci; // dwdot[CH4]/d[CH4]
  J[309] += dqdci; // dwdot[CH3]/d[CH4]
  // d()/d[CH3]
  dqdci = -k_r * sc[1];
  J[325] += dqdci; // dwdot[H2]/d[CH3]
  J[326] -= dqdci; // dwdot[H]/d[CH3]
  J[335] -= dqdci; // dwdot[CH4]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[713] -= dqdT; // dwdot[CH4]/dT
  J[714] += dqdT; // dwdot[CH3]/dT

  // reaction 26: CH4 + O <=> CH3 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[11];
  k_f = 1020 * exp(1.5 * tc[0] - (4327.66332592059) * invT);
  dlnkfdT = 1.5 * invT + (4327.66332592059) * invT2;
  // reverse
  phi_r = sc[6] * sc[12];
  Kc = exp(g_RT[4] - g_RT[6] + g_RT[11] - g_RT[12]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[6] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[6] += q;  // OH
  wdot[11] -= q; // CH4
  wdot[12] += q; // CH3
  // d()/d[O]
  dqdci = +k_f * sc[11];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[119] -= dqdci; // dwdot[CH4]/d[O]
  J[120] += dqdci; // dwdot[CH3]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[12];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[173] -= dqdci; // dwdot[CH4]/d[OH]
  J[174] += dqdci; // dwdot[CH3]/d[OH]
  // d()/d[CH4]
  dqdci = +k_f * sc[4];
  J[301] -= dqdci; // dwdot[O]/d[CH4]
  J[303] += dqdci; // dwdot[OH]/d[CH4]
  J[308] -= dqdci; // dwdot[CH4]/d[CH4]
  J[309] += dqdci; // dwdot[CH3]/d[CH4]
  // d()/d[CH3]
  dqdci = -k_r * sc[6];
  J[328] -= dqdci; // dwdot[O]/d[CH3]
  J[330] += dqdci; // dwdot[OH]/d[CH3]
  J[335] -= dqdci; // dwdot[CH4]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[713] -= dqdT; // dwdot[CH4]/dT
  J[714] += dqdT; // dwdot[CH3]/dT

  // reaction 27: CH4 + OH <=> CH3 + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[11];
  k_f = 0.0583 * exp(2.6 * tc[0] - (1102.04449811234) * invT);
  dlnkfdT = 2.6 * invT + (1102.04449811234) * invT2;
  // reverse
  phi_r = sc[5] * sc[12];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[11] - g_RT[12]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[11]) + (h_RT[5] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[11] -= q; // CH4
  wdot[12] += q; // CH3
  // d()/d[H2O]
  dqdci = -k_r * sc[12];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[146] -= dqdci; // dwdot[CH4]/d[H2O]
  J[147] += dqdci; // dwdot[CH3]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[11];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[173] -= dqdci; // dwdot[CH4]/d[OH]
  J[174] += dqdci; // dwdot[CH3]/d[OH]
  // d()/d[CH4]
  dqdci = +k_f * sc[6];
  J[302] += dqdci; // dwdot[H2O]/d[CH4]
  J[303] -= dqdci; // dwdot[OH]/d[CH4]
  J[308] -= dqdci; // dwdot[CH4]/d[CH4]
  J[309] += dqdci; // dwdot[CH3]/d[CH4]
  // d()/d[CH3]
  dqdci = -k_r * sc[5];
  J[329] += dqdci; // dwdot[H2O]/d[CH3]
  J[330] -= dqdci; // dwdot[OH]/d[CH3]
  J[335] -= dqdci; // dwdot[CH4]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[713] -= dqdT; // dwdot[CH4]/dT
  J[714] += dqdT; // dwdot[CH3]/dT

  // reaction 28: CH4 + HO2 <=> CH3 + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[11];
  k_f = 1.13e-05 * exp(3.74 * tc[0] - (10572.5821485572) * invT);
  dlnkfdT = 3.74 * invT + (10572.5821485572) * invT2;
  // reverse
  phi_r = sc[7] * sc[12];
  Kc = exp(-g_RT[7] + g_RT[8] + g_RT[11] - g_RT[12]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[11]) + (h_RT[7] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[11] -= q; // CH4
  wdot[12] += q; // CH3
  // d()/d[H2O2]
  dqdci = -k_r * sc[12];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[200] -= dqdci; // dwdot[CH4]/d[H2O2]
  J[201] += dqdci; // dwdot[CH3]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[11];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[227] -= dqdci; // dwdot[CH4]/d[HO2]
  J[228] += dqdci; // dwdot[CH3]/d[HO2]
  // d()/d[CH4]
  dqdci = +k_f * sc[8];
  J[304] += dqdci; // dwdot[H2O2]/d[CH4]
  J[305] -= dqdci; // dwdot[HO2]/d[CH4]
  J[308] -= dqdci; // dwdot[CH4]/d[CH4]
  J[309] += dqdci; // dwdot[CH3]/d[CH4]
  // d()/d[CH3]
  dqdci = -k_r * sc[7];
  J[331] += dqdci; // dwdot[H2O2]/d[CH3]
  J[332] -= dqdci; // dwdot[HO2]/d[CH3]
  J[335] -= dqdci; // dwdot[CH4]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[713] -= dqdT; // dwdot[CH4]/dT
  J[714] += dqdT; // dwdot[CH3]/dT

  // reaction 29: CH3 + HO2 <=> CH4 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[12];
  k_f = 0.116 * exp(2.23 * tc[0] - (-1520.72076406186) * invT);
  dlnkfdT = 2.23 * invT + (-1520.72076406186) * invT2;
  // reverse
  phi_r = sc[3] * sc[11];
  Kc = exp(-g_RT[3] + g_RT[8] - g_RT[11] + g_RT[12]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[3] + h_RT[11]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q;  // O2
  wdot[8] -= q;  // HO2
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  // d()/d[O2]
  dqdci = -k_r * sc[11];
  J[84] += dqdci; // dwdot[O2]/d[O2]
  J[89] -= dqdci; // dwdot[HO2]/d[O2]
  J[92] += dqdci; // dwdot[CH4]/d[O2]
  J[93] -= dqdci; // dwdot[CH3]/d[O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[12];
  J[219] += dqdci; // dwdot[O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[227] += dqdci; // dwdot[CH4]/d[HO2]
  J[228] -= dqdci; // dwdot[CH3]/d[HO2]
  // d()/d[CH4]
  dqdci = -k_r * sc[3];
  J[300] += dqdci; // dwdot[O2]/d[CH4]
  J[305] -= dqdci; // dwdot[HO2]/d[CH4]
  J[308] += dqdci; // dwdot[CH4]/d[CH4]
  J[309] -= dqdci; // dwdot[CH3]/d[CH4]
  // d()/d[CH3]
  dqdci = +k_f * sc[8];
  J[327] += dqdci; // dwdot[O2]/d[CH3]
  J[332] -= dqdci; // dwdot[HO2]/d[CH3]
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  // d()/dT
  J[705] += dqdT; // dwdot[O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT

  // reaction 30: CH2 + CH4 <=> 2 CH3
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[11] * sc[13];
  k_f = 2.46 * exp(2 * tc[0] - (4161.60182620503) * invT);
  dlnkfdT = 2 * invT + (4161.60182620503) * invT2;
  // reverse
  phi_r = (sc[12] * sc[12]);
  Kc = exp(g_RT[11] - 2.000000 * g_RT[12] + g_RT[13]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[11] + h_RT[13]) + (2.000000 * h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[11] -= q;     // CH4
  wdot[12] += 2 * q; // CH3
  wdot[13] -= q;     // CH2
  // d()/d[CH4]
  dqdci = +k_f * sc[13];
  J[308] -= dqdci;     // dwdot[CH4]/d[CH4]
  J[309] += 2 * dqdci; // dwdot[CH3]/d[CH4]
  J[310] -= dqdci;     // dwdot[CH2]/d[CH4]
  // d()/d[CH3]
  dqdci = -k_r * 2.000000 * sc[12];
  J[335] -= dqdci;     // dwdot[CH4]/d[CH3]
  J[336] += 2 * dqdci; // dwdot[CH3]/d[CH3]
  J[337] -= dqdci;     // dwdot[CH2]/d[CH3]
  // d()/d[CH2]
  dqdci = +k_f * sc[11];
  J[362] -= dqdci;     // dwdot[CH4]/d[CH2]
  J[363] += 2 * dqdci; // dwdot[CH3]/d[CH2]
  J[364] -= dqdci;     // dwdot[CH2]/d[CH2]
  // d()/dT
  J[713] -= dqdT;     // dwdot[CH4]/dT
  J[714] += 2 * dqdT; // dwdot[CH3]/dT
  J[715] -= dqdT;     // dwdot[CH2]/dT

  // reaction 31: CH2(S) + N2 <=> CH2 + N2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[14];
  k_f = 15000000 * exp(-(301.929999482832) * invT);
  dlnkfdT = (301.929999482832) * invT2;
  // reverse
  phi_r = sc[0] * sc[13];
  Kc = exp(g_RT[0] - g_RT[0] - g_RT[13] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[14]) + (h_RT[0] + h_RT[13]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[13] += q; // CH2
  wdot[14] -= q; // CH2(S)
  // d()/d[N2]
  dqdci = +k_f * sc[14] - k_r * sc[13];
  J[13] += dqdci; // dwdot[CH2]/d[N2]
  J[14] -= dqdci; // dwdot[CH2(S)]/d[N2]
  // d()/d[CH2]
  dqdci = -k_r * sc[0];
  J[364] += dqdci; // dwdot[CH2]/d[CH2]
  J[365] -= dqdci; // dwdot[CH2(S)]/d[CH2]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[0];
  J[391] += dqdci; // dwdot[CH2]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[715] += dqdT; // dwdot[CH2]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 32: CH2(S) + H2O <=> CH2 + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[14];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[5] * sc[13];
  Kc = exp(g_RT[5] - g_RT[5] - g_RT[13] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[14]) + (h_RT[5] + h_RT[13]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[13] += q; // CH2
  wdot[14] -= q; // CH2(S)
  // d()/d[H2O]
  dqdci = +k_f * sc[14] - k_r * sc[13];
  J[148] += dqdci; // dwdot[CH2]/d[H2O]
  J[149] -= dqdci; // dwdot[CH2(S)]/d[H2O]
  // d()/d[CH2]
  dqdci = -k_r * sc[5];
  J[364] += dqdci; // dwdot[CH2]/d[CH2]
  J[365] -= dqdci; // dwdot[CH2(S)]/d[CH2]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[5];
  J[391] += dqdci; // dwdot[CH2]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[715] += dqdT; // dwdot[CH2]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 33: CH2(S) + CO <=> CH2 + CO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[9] * sc[14];
  k_f = 9000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[9] * sc[13];
  Kc = exp(g_RT[9] - g_RT[9] - g_RT[13] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[9] + h_RT[14]) + (h_RT[9] + h_RT[13]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[13] += q; // CH2
  wdot[14] -= q; // CH2(S)
  // d()/d[CO]
  dqdci = +k_f * sc[14] - k_r * sc[13];
  J[256] += dqdci; // dwdot[CH2]/d[CO]
  J[257] -= dqdci; // dwdot[CH2(S)]/d[CO]
  // d()/d[CH2]
  dqdci = -k_r * sc[9];
  J[364] += dqdci; // dwdot[CH2]/d[CH2]
  J[365] -= dqdci; // dwdot[CH2(S)]/d[CH2]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[9];
  J[391] += dqdci; // dwdot[CH2]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[715] += dqdT; // dwdot[CH2]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 34: CH2(S) + CO2 <=> CH2 + CO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[10] * sc[14];
  k_f = 7000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[10] * sc[13];
  Kc = exp(g_RT[10] - g_RT[10] - g_RT[13] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[10] + h_RT[14]) + (h_RT[10] + h_RT[13]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[13] += q; // CH2
  wdot[14] -= q; // CH2(S)
  // d()/d[CO2]
  dqdci = +k_f * sc[14] - k_r * sc[13];
  J[283] += dqdci; // dwdot[CH2]/d[CO2]
  J[284] -= dqdci; // dwdot[CH2(S)]/d[CO2]
  // d()/d[CH2]
  dqdci = -k_r * sc[10];
  J[364] += dqdci; // dwdot[CH2]/d[CH2]
  J[365] -= dqdci; // dwdot[CH2(S)]/d[CH2]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[10];
  J[391] += dqdci; // dwdot[CH2]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[715] += dqdT; // dwdot[CH2]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 35: CH2(S) + O2 => CO + H + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[14];
  k_f = 28000000;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[2] += q;  // H
  wdot[3] -= q;  // O2
  wdot[6] += q;  // OH
  wdot[9] += q;  // CO
  wdot[14] -= q; // CH2(S)
  // d()/d[O2]
  dqdci = +k_f * sc[14];
  J[83] += dqdci; // dwdot[H]/d[O2]
  J[84] -= dqdci; // dwdot[O2]/d[O2]
  J[87] += dqdci; // dwdot[OH]/d[O2]
  J[90] += dqdci; // dwdot[CO]/d[O2]
  J[95] -= dqdci; // dwdot[CH2(S)]/d[O2]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[3];
  J[380] += dqdci; // dwdot[H]/d[CH2(S)]
  J[381] -= dqdci; // dwdot[O2]/d[CH2(S)]
  J[384] += dqdci; // dwdot[OH]/d[CH2(S)]
  J[387] += dqdci; // dwdot[CO]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 36: CH2(S) + O2 <=> CO + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[14];
  k_f = 12000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[5] * sc[9];
  Kc = exp(g_RT[3] - g_RT[5] - g_RT[9] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[14]) + (h_RT[5] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[5] += q;  // H2O
  wdot[9] += q;  // CO
  wdot[14] -= q; // CH2(S)
  // d()/d[O2]
  dqdci = +k_f * sc[14];
  J[84] -= dqdci; // dwdot[O2]/d[O2]
  J[86] += dqdci; // dwdot[H2O]/d[O2]
  J[90] += dqdci; // dwdot[CO]/d[O2]
  J[95] -= dqdci; // dwdot[CH2(S)]/d[O2]
  // d()/d[H2O]
  dqdci = -k_r * sc[9];
  J[138] -= dqdci; // dwdot[O2]/d[H2O]
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[144] += dqdci; // dwdot[CO]/d[H2O]
  J[149] -= dqdci; // dwdot[CH2(S)]/d[H2O]
  // d()/d[CO]
  dqdci = -k_r * sc[5];
  J[246] -= dqdci; // dwdot[O2]/d[CO]
  J[248] += dqdci; // dwdot[H2O]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[257] -= dqdci; // dwdot[CH2(S)]/d[CO]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[3];
  J[381] -= dqdci; // dwdot[O2]/d[CH2(S)]
  J[383] += dqdci; // dwdot[H2O]/d[CH2(S)]
  J[387] += dqdci; // dwdot[CO]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 37: CH2(S) + O <=> CO + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[14];
  k_f = 15000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[9];
  Kc = exp(-g_RT[1] + g_RT[4] - g_RT[9] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[14]) + (h_RT[1] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[4] -= q;  // O
  wdot[9] += q;  // CO
  wdot[14] -= q; // CH2(S)
  // d()/d[H2]
  dqdci = -k_r * sc[9];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[31] -= dqdci; // dwdot[O]/d[H2]
  J[36] += dqdci; // dwdot[CO]/d[H2]
  J[41] -= dqdci; // dwdot[CH2(S)]/d[H2]
  // d()/d[O]
  dqdci = +k_f * sc[14];
  J[109] += dqdci; // dwdot[H2]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[117] += dqdci; // dwdot[CO]/d[O]
  J[122] -= dqdci; // dwdot[CH2(S)]/d[O]
  // d()/d[CO]
  dqdci = -k_r * sc[1];
  J[244] += dqdci; // dwdot[H2]/d[CO]
  J[247] -= dqdci; // dwdot[O]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[257] -= dqdci; // dwdot[CH2(S)]/d[CO]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[4];
  J[379] += dqdci; // dwdot[H2]/d[CH2(S)]
  J[382] -= dqdci; // dwdot[O]/d[CH2(S)]
  J[387] += dqdci; // dwdot[CO]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 38: CH2(S) + O <=> H + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[14];
  k_f = 15000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[21];
  Kc = exp(-g_RT[2] + g_RT[4] + g_RT[14] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[14]) + (h_RT[2] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[4] -= q;  // O
  wdot[14] -= q; // CH2(S)
  wdot[21] += q; // HCO
  // d()/d[H]
  dqdci = -k_r * sc[21];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[58] -= dqdci; // dwdot[O]/d[H]
  J[68] -= dqdci; // dwdot[CH2(S)]/d[H]
  J[75] += dqdci; // dwdot[HCO]/d[H]
  // d()/d[O]
  dqdci = +k_f * sc[14];
  J[110] += dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[122] -= dqdci; // dwdot[CH2(S)]/d[O]
  J[129] += dqdci; // dwdot[HCO]/d[O]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[4];
  J[380] += dqdci; // dwdot[H]/d[CH2(S)]
  J[382] -= dqdci; // dwdot[O]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  J[399] += dqdci; // dwdot[HCO]/d[CH2(S)]
  // d()/d[HCO]
  dqdci = -k_r * sc[2];
  J[569] += dqdci; // dwdot[H]/d[HCO]
  J[571] -= dqdci; // dwdot[O]/d[HCO]
  J[581] -= dqdci; // dwdot[CH2(S)]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 39: CH2(S) + H2 <=> CH3 + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[14];
  k_f = 70000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[12];
  Kc = exp(g_RT[1] - g_RT[2] - g_RT[12] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[2] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q;  // H2
  wdot[2] += q;  // H
  wdot[12] += q; // CH3
  wdot[14] -= q; // CH2(S)
  // d()/d[H2]
  dqdci = +k_f * sc[14];
  J[28] -= dqdci; // dwdot[H2]/d[H2]
  J[29] += dqdci; // dwdot[H]/d[H2]
  J[39] += dqdci; // dwdot[CH3]/d[H2]
  J[41] -= dqdci; // dwdot[CH2(S)]/d[H2]
  // d()/d[H]
  dqdci = -k_r * sc[12];
  J[55] -= dqdci; // dwdot[H2]/d[H]
  J[56] += dqdci; // dwdot[H]/d[H]
  J[66] += dqdci; // dwdot[CH3]/d[H]
  J[68] -= dqdci; // dwdot[CH2(S)]/d[H]
  // d()/d[CH3]
  dqdci = -k_r * sc[2];
  J[325] -= dqdci; // dwdot[H2]/d[CH3]
  J[326] += dqdci; // dwdot[H]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  J[338] -= dqdci; // dwdot[CH2(S)]/d[CH3]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[1];
  J[379] -= dqdci; // dwdot[H2]/d[CH2(S)]
  J[380] += dqdci; // dwdot[H]/d[CH2(S)]
  J[390] += dqdci; // dwdot[CH3]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[703] -= dqdT; // dwdot[H2]/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[714] += dqdT; // dwdot[CH3]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 40: CH2(S) + H <=> CH + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[14];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[16];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[14] - g_RT[16]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[14]) + (h_RT[1] + h_RT[16]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[14] -= q; // CH2(S)
  wdot[16] += q; // CH
  // d()/d[H2]
  dqdci = -k_r * sc[16];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[41] -= dqdci; // dwdot[CH2(S)]/d[H2]
  J[43] += dqdci; // dwdot[CH]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[14];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[68] -= dqdci; // dwdot[CH2(S)]/d[H]
  J[70] += dqdci; // dwdot[CH]/d[H]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[2];
  J[379] += dqdci; // dwdot[H2]/d[CH2(S)]
  J[380] -= dqdci; // dwdot[H]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  J[394] += dqdci; // dwdot[CH]/d[CH2(S)]
  // d()/d[CH]
  dqdci = -k_r * sc[1];
  J[433] += dqdci; // dwdot[H2]/d[CH]
  J[434] -= dqdci; // dwdot[H]/d[CH]
  J[446] -= dqdci; // dwdot[CH2(S)]/d[CH]
  J[448] += dqdci; // dwdot[CH]/d[CH]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT
  J[718] += dqdT; // dwdot[CH]/dT

  // reaction 41: CH2(S) + OH <=> CH2O + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[14];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[20];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[14] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[14]) + (h_RT[2] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[14] -= q; // CH2(S)
  wdot[20] += q; // CH2O
  // d()/d[H]
  dqdci = -k_r * sc[20];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[68] -= dqdci; // dwdot[CH2(S)]/d[H]
  J[74] += dqdci; // dwdot[CH2O]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[14];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[176] -= dqdci; // dwdot[CH2(S)]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[6];
  J[380] += dqdci; // dwdot[H]/d[CH2(S)]
  J[384] -= dqdci; // dwdot[OH]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  J[398] += dqdci; // dwdot[CH2O]/d[CH2(S)]
  // d()/d[CH2O]
  dqdci = -k_r * sc[2];
  J[542] += dqdci; // dwdot[H]/d[CH2O]
  J[546] -= dqdci; // dwdot[OH]/d[CH2O]
  J[554] -= dqdci; // dwdot[CH2(S)]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 42: CH2(S) + CO2 <=> CH2O + CO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[10] * sc[14];
  k_f = 14000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[9] * sc[20];
  Kc = exp(-g_RT[9] + g_RT[10] + g_RT[14] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[10] + h_RT[14]) + (h_RT[9] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] += q;  // CO
  wdot[10] -= q; // CO2
  wdot[14] -= q; // CH2(S)
  wdot[20] += q; // CH2O
  // d()/d[CO]
  dqdci = -k_r * sc[20];
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[253] -= dqdci; // dwdot[CO2]/d[CO]
  J[257] -= dqdci; // dwdot[CH2(S)]/d[CO]
  J[263] += dqdci; // dwdot[CH2O]/d[CO]
  // d()/d[CO2]
  dqdci = +k_f * sc[14];
  J[279] += dqdci; // dwdot[CO]/d[CO2]
  J[280] -= dqdci; // dwdot[CO2]/d[CO2]
  J[284] -= dqdci; // dwdot[CH2(S)]/d[CO2]
  J[290] += dqdci; // dwdot[CH2O]/d[CO2]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[10];
  J[387] += dqdci; // dwdot[CO]/d[CH2(S)]
  J[388] -= dqdci; // dwdot[CO2]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  J[398] += dqdci; // dwdot[CH2O]/d[CH2(S)]
  // d()/d[CH2O]
  dqdci = -k_r * sc[9];
  J[549] += dqdci; // dwdot[CO]/d[CH2O]
  J[550] -= dqdci; // dwdot[CO2]/d[CH2O]
  J[554] -= dqdci; // dwdot[CH2(S)]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[712] -= dqdT; // dwdot[CO2]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 44: CH2 + O2 <=> HCO + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[13];
  k_f = 10600000 * exp(-(754.82499870708) * invT);
  dlnkfdT = (754.82499870708) * invT2;
  // reverse
  phi_r = sc[6] * sc[21];
  Kc = exp(g_RT[3] - g_RT[6] + g_RT[13] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[13]) + (h_RT[6] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[6] += q;  // OH
  wdot[13] -= q; // CH2
  wdot[21] += q; // HCO
  // d()/d[O2]
  dqdci = +k_f * sc[13];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[87] += dqdci;  // dwdot[OH]/d[O2]
  J[94] -= dqdci;  // dwdot[CH2]/d[O2]
  J[102] += dqdci; // dwdot[HCO]/d[O2]
  // d()/d[OH]
  dqdci = -k_r * sc[21];
  J[165] -= dqdci; // dwdot[O2]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[175] -= dqdci; // dwdot[CH2]/d[OH]
  J[183] += dqdci; // dwdot[HCO]/d[OH]
  // d()/d[CH2]
  dqdci = +k_f * sc[3];
  J[354] -= dqdci; // dwdot[O2]/d[CH2]
  J[357] += dqdci; // dwdot[OH]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  J[372] += dqdci; // dwdot[HCO]/d[CH2]
  // d()/d[HCO]
  dqdci = -k_r * sc[6];
  J[570] -= dqdci; // dwdot[O2]/d[HCO]
  J[573] += dqdci; // dwdot[OH]/d[HCO]
  J[580] -= dqdci; // dwdot[CH2]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 45: CH2 + O2 => CO2 + 2 H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[13];
  k_f = 2640000 * exp(-(754.82499870708) * invT);
  dlnkfdT = (754.82499870708) * invT2;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[2] += 2 * q; // H
  wdot[3] -= q;     // O2
  wdot[10] += q;    // CO2
  wdot[13] -= q;    // CH2
  // d()/d[O2]
  dqdci = +k_f * sc[13];
  J[83] += 2 * dqdci; // dwdot[H]/d[O2]
  J[84] -= dqdci;     // dwdot[O2]/d[O2]
  J[91] += dqdci;     // dwdot[CO2]/d[O2]
  J[94] -= dqdci;     // dwdot[CH2]/d[O2]
  // d()/d[CH2]
  dqdci = +k_f * sc[3];
  J[353] += 2 * dqdci; // dwdot[H]/d[CH2]
  J[354] -= dqdci;     // dwdot[O2]/d[CH2]
  J[361] += dqdci;     // dwdot[CO2]/d[CH2]
  J[364] -= dqdci;     // dwdot[CH2]/d[CH2]
  // d()/dT
  J[704] += 2 * dqdT; // dwdot[H]/dT
  J[705] -= dqdT;     // dwdot[O2]/dT
  J[712] += dqdT;     // dwdot[CO2]/dT
  J[715] -= dqdT;     // dwdot[CH2]/dT

  // reaction 46: CH2 + O => CO + 2 H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[13];
  k_f = 50000000;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[2] += 2 * q; // H
  wdot[4] -= q;     // O
  wdot[9] += q;     // CO
  wdot[13] -= q;    // CH2
  // d()/d[O]
  dqdci = +k_f * sc[13];
  J[110] += 2 * dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci;     // dwdot[O]/d[O]
  J[117] += dqdci;     // dwdot[CO]/d[O]
  J[121] -= dqdci;     // dwdot[CH2]/d[O]
  // d()/d[CH2]
  dqdci = +k_f * sc[4];
  J[353] += 2 * dqdci; // dwdot[H]/d[CH2]
  J[355] -= dqdci;     // dwdot[O]/d[CH2]
  J[360] += dqdci;     // dwdot[CO]/d[CH2]
  J[364] -= dqdci;     // dwdot[CH2]/d[CH2]
  // d()/dT
  J[704] += 2 * dqdT; // dwdot[H]/dT
  J[706] -= dqdT;     // dwdot[O]/dT
  J[711] += dqdT;     // dwdot[CO]/dT
  J[715] -= dqdT;     // dwdot[CH2]/dT

  // reaction 47: CH2 + H <=> CH + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[13];
  k_f = 283000000 * exp(-(369.864249366469) * invT);
  dlnkfdT = (369.864249366469) * invT2;
  // reverse
  phi_r = sc[1] * sc[16];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[13] - g_RT[16]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[1] + h_RT[16]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[13] -= q; // CH2
  wdot[16] += q; // CH
  // d()/d[H2]
  dqdci = -k_r * sc[16];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[40] -= dqdci; // dwdot[CH2]/d[H2]
  J[43] += dqdci; // dwdot[CH]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[13];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[67] -= dqdci; // dwdot[CH2]/d[H]
  J[70] += dqdci; // dwdot[CH]/d[H]
  // d()/d[CH2]
  dqdci = +k_f * sc[2];
  J[352] += dqdci; // dwdot[H2]/d[CH2]
  J[353] -= dqdci; // dwdot[H]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  J[367] += dqdci; // dwdot[CH]/d[CH2]
  // d()/d[CH]
  dqdci = -k_r * sc[1];
  J[433] += dqdci; // dwdot[H2]/d[CH]
  J[434] -= dqdci; // dwdot[H]/d[CH]
  J[445] -= dqdci; // dwdot[CH2]/d[CH]
  J[448] += dqdci; // dwdot[CH]/d[CH]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT
  J[718] += dqdT; // dwdot[CH]/dT

  // reaction 48: CH2 + OH <=> CH + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[13];
  k_f = 11.3 * exp(2 * tc[0] - (1509.64999741416) * invT);
  dlnkfdT = 2 * invT + (1509.64999741416) * invT2;
  // reverse
  phi_r = sc[5] * sc[16];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[13] - g_RT[16]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[13]) + (h_RT[5] + h_RT[16]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[13] -= q; // CH2
  wdot[16] += q; // CH
  // d()/d[H2O]
  dqdci = -k_r * sc[16];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[148] -= dqdci; // dwdot[CH2]/d[H2O]
  J[151] += dqdci; // dwdot[CH]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[13];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[175] -= dqdci; // dwdot[CH2]/d[OH]
  J[178] += dqdci; // dwdot[CH]/d[OH]
  // d()/d[CH2]
  dqdci = +k_f * sc[6];
  J[356] += dqdci; // dwdot[H2O]/d[CH2]
  J[357] -= dqdci; // dwdot[OH]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  J[367] += dqdci; // dwdot[CH]/d[CH2]
  // d()/d[CH]
  dqdci = -k_r * sc[5];
  J[437] += dqdci; // dwdot[H2O]/d[CH]
  J[438] -= dqdci; // dwdot[OH]/d[CH]
  J[445] -= dqdci; // dwdot[CH2]/d[CH]
  J[448] += dqdci; // dwdot[CH]/d[CH]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT
  J[718] += dqdT; // dwdot[CH]/dT

  // reaction 49: CH + O2 <=> HCO + O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[16];
  k_f = 33000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[4] * sc[21];
  Kc = exp(g_RT[3] - g_RT[4] + g_RT[16] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[16]) + (h_RT[4] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[4] += q;  // O
  wdot[16] -= q; // CH
  wdot[21] += q; // HCO
  // d()/d[O2]
  dqdci = +k_f * sc[16];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[85] += dqdci;  // dwdot[O]/d[O2]
  J[97] -= dqdci;  // dwdot[CH]/d[O2]
  J[102] += dqdci; // dwdot[HCO]/d[O2]
  // d()/d[O]
  dqdci = -k_r * sc[21];
  J[111] -= dqdci; // dwdot[O2]/d[O]
  J[112] += dqdci; // dwdot[O]/d[O]
  J[124] -= dqdci; // dwdot[CH]/d[O]
  J[129] += dqdci; // dwdot[HCO]/d[O]
  // d()/d[CH]
  dqdci = +k_f * sc[3];
  J[435] -= dqdci; // dwdot[O2]/d[CH]
  J[436] += dqdci; // dwdot[O]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  J[453] += dqdci; // dwdot[HCO]/d[CH]
  // d()/d[HCO]
  dqdci = -k_r * sc[4];
  J[570] -= dqdci; // dwdot[O2]/d[HCO]
  J[571] += dqdci; // dwdot[O]/d[HCO]
  J[583] -= dqdci; // dwdot[CH]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[706] += dqdT; // dwdot[O]/dT
  J[718] -= dqdT; // dwdot[CH]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 50: CH + O <=> CO + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[16];
  k_f = 57000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[9];
  Kc = exp(-g_RT[2] + g_RT[4] - g_RT[9] + g_RT[16]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[16]) + (h_RT[2] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[4] -= q;  // O
  wdot[9] += q;  // CO
  wdot[16] -= q; // CH
  // d()/d[H]
  dqdci = -k_r * sc[9];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[58] -= dqdci; // dwdot[O]/d[H]
  J[63] += dqdci; // dwdot[CO]/d[H]
  J[70] -= dqdci; // dwdot[CH]/d[H]
  // d()/d[O]
  dqdci = +k_f * sc[16];
  J[110] += dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[117] += dqdci; // dwdot[CO]/d[O]
  J[124] -= dqdci; // dwdot[CH]/d[O]
  // d()/d[CO]
  dqdci = -k_r * sc[2];
  J[245] += dqdci; // dwdot[H]/d[CO]
  J[247] -= dqdci; // dwdot[O]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[259] -= dqdci; // dwdot[CH]/d[CO]
  // d()/d[CH]
  dqdci = +k_f * sc[4];
  J[434] += dqdci; // dwdot[H]/d[CH]
  J[436] -= dqdci; // dwdot[O]/d[CH]
  J[441] += dqdci; // dwdot[CO]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[718] -= dqdT; // dwdot[CH]/dT

  // reaction 51: CH + H <=> C + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[16];
  k_f = 110000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[15];
  Kc = exp(-g_RT[1] + g_RT[2] - g_RT[15] + g_RT[16]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[16]) + (h_RT[1] + h_RT[15]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[15] += q; // C
  wdot[16] -= q; // CH
  // d()/d[H2]
  dqdci = -k_r * sc[15];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[42] += dqdci; // dwdot[C]/d[H2]
  J[43] -= dqdci; // dwdot[CH]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[16];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[69] += dqdci; // dwdot[C]/d[H]
  J[70] -= dqdci; // dwdot[CH]/d[H]
  // d()/d[C]
  dqdci = -k_r * sc[1];
  J[406] += dqdci; // dwdot[H2]/d[C]
  J[407] -= dqdci; // dwdot[H]/d[C]
  J[420] += dqdci; // dwdot[C]/d[C]
  J[421] -= dqdci; // dwdot[CH]/d[C]
  // d()/d[CH]
  dqdci = +k_f * sc[2];
  J[433] += dqdci; // dwdot[H2]/d[CH]
  J[434] -= dqdci; // dwdot[H]/d[CH]
  J[447] += dqdci; // dwdot[C]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[717] += dqdT; // dwdot[C]/dT
  J[718] -= dqdT; // dwdot[CH]/dT

  // reaction 52: CH + OH <=> H + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[16];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[21];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[16] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[16]) + (h_RT[2] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[16] -= q; // CH
  wdot[21] += q; // HCO
  // d()/d[H]
  dqdci = -k_r * sc[21];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[70] -= dqdci; // dwdot[CH]/d[H]
  J[75] += dqdci; // dwdot[HCO]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[16];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[178] -= dqdci; // dwdot[CH]/d[OH]
  J[183] += dqdci; // dwdot[HCO]/d[OH]
  // d()/d[CH]
  dqdci = +k_f * sc[6];
  J[434] += dqdci; // dwdot[H]/d[CH]
  J[438] -= dqdci; // dwdot[OH]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  J[453] += dqdci; // dwdot[HCO]/d[CH]
  // d()/d[HCO]
  dqdci = -k_r * sc[2];
  J[569] += dqdci; // dwdot[H]/d[HCO]
  J[573] -= dqdci; // dwdot[OH]/d[HCO]
  J[583] -= dqdci; // dwdot[CH]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[718] -= dqdT; // dwdot[CH]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 53: CH + H2O <=> CH2O + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[16];
  k_f = 17130000 * exp(-(-379.928582682563) * invT);
  dlnkfdT = (-379.928582682563) * invT2;
  // reverse
  phi_r = sc[2] * sc[20];
  Kc = exp(-g_RT[2] + g_RT[5] + g_RT[16] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[16]) + (h_RT[2] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[5] -= q;  // H2O
  wdot[16] -= q; // CH
  wdot[20] += q; // CH2O
  // d()/d[H]
  dqdci = -k_r * sc[20];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[59] -= dqdci; // dwdot[H2O]/d[H]
  J[70] -= dqdci; // dwdot[CH]/d[H]
  J[74] += dqdci; // dwdot[CH2O]/d[H]
  // d()/d[H2O]
  dqdci = +k_f * sc[16];
  J[137] += dqdci; // dwdot[H]/d[H2O]
  J[140] -= dqdci; // dwdot[H2O]/d[H2O]
  J[151] -= dqdci; // dwdot[CH]/d[H2O]
  J[155] += dqdci; // dwdot[CH2O]/d[H2O]
  // d()/d[CH]
  dqdci = +k_f * sc[5];
  J[434] += dqdci; // dwdot[H]/d[CH]
  J[437] -= dqdci; // dwdot[H2O]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  J[452] += dqdci; // dwdot[CH2O]/d[CH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[2];
  J[542] += dqdci; // dwdot[H]/d[CH2O]
  J[545] -= dqdci; // dwdot[H2O]/d[CH2O]
  J[556] -= dqdci; // dwdot[CH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[707] -= dqdT; // dwdot[H2O]/dT
  J[718] -= dqdT; // dwdot[CH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 54: CH + CO2 <=> CO + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[10] * sc[16];
  k_f = 1700000 * exp(-(344.703416076233) * invT);
  dlnkfdT = (344.703416076233) * invT2;
  // reverse
  phi_r = sc[9] * sc[21];
  Kc = exp(-g_RT[9] + g_RT[10] + g_RT[16] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[10] + h_RT[16]) + (h_RT[9] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] += q;  // CO
  wdot[10] -= q; // CO2
  wdot[16] -= q; // CH
  wdot[21] += q; // HCO
  // d()/d[CO]
  dqdci = -k_r * sc[21];
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[253] -= dqdci; // dwdot[CO2]/d[CO]
  J[259] -= dqdci; // dwdot[CH]/d[CO]
  J[264] += dqdci; // dwdot[HCO]/d[CO]
  // d()/d[CO2]
  dqdci = +k_f * sc[16];
  J[279] += dqdci; // dwdot[CO]/d[CO2]
  J[280] -= dqdci; // dwdot[CO2]/d[CO2]
  J[286] -= dqdci; // dwdot[CH]/d[CO2]
  J[291] += dqdci; // dwdot[HCO]/d[CO2]
  // d()/d[CH]
  dqdci = +k_f * sc[10];
  J[441] += dqdci; // dwdot[CO]/d[CH]
  J[442] -= dqdci; // dwdot[CO2]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  J[453] += dqdci; // dwdot[HCO]/d[CH]
  // d()/d[HCO]
  dqdci = -k_r * sc[9];
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[577] -= dqdci; // dwdot[CO2]/d[HCO]
  J[583] -= dqdci; // dwdot[CH]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[712] -= dqdT; // dwdot[CO2]/dT
  J[718] -= dqdT; // dwdot[CH]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 55: C + OH <=> CO + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[15];
  k_f = 50000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[9];
  Kc = exp(-g_RT[2] + g_RT[6] - g_RT[9] + g_RT[15]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[15]) + (h_RT[2] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[9] += q;  // CO
  wdot[15] -= q; // C
  // d()/d[H]
  dqdci = -k_r * sc[9];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[63] += dqdci; // dwdot[CO]/d[H]
  J[69] -= dqdci; // dwdot[C]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[15];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[171] += dqdci; // dwdot[CO]/d[OH]
  J[177] -= dqdci; // dwdot[C]/d[OH]
  // d()/d[CO]
  dqdci = -k_r * sc[2];
  J[245] += dqdci; // dwdot[H]/d[CO]
  J[249] -= dqdci; // dwdot[OH]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[258] -= dqdci; // dwdot[C]/d[CO]
  // d()/d[C]
  dqdci = +k_f * sc[6];
  J[407] += dqdci; // dwdot[H]/d[C]
  J[411] -= dqdci; // dwdot[OH]/d[C]
  J[414] += dqdci; // dwdot[CO]/d[C]
  J[420] -= dqdci; // dwdot[C]/d[C]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[717] -= dqdT; // dwdot[C]/dT

  // reaction 56: C + O2 <=> CO + O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[15];
  k_f = 50000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[4] * sc[9];
  Kc = exp(g_RT[3] - g_RT[4] - g_RT[9] + g_RT[15]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[15]) + (h_RT[4] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[4] += q;  // O
  wdot[9] += q;  // CO
  wdot[15] -= q; // C
  // d()/d[O2]
  dqdci = +k_f * sc[15];
  J[84] -= dqdci; // dwdot[O2]/d[O2]
  J[85] += dqdci; // dwdot[O]/d[O2]
  J[90] += dqdci; // dwdot[CO]/d[O2]
  J[96] -= dqdci; // dwdot[C]/d[O2]
  // d()/d[O]
  dqdci = -k_r * sc[9];
  J[111] -= dqdci; // dwdot[O2]/d[O]
  J[112] += dqdci; // dwdot[O]/d[O]
  J[117] += dqdci; // dwdot[CO]/d[O]
  J[123] -= dqdci; // dwdot[C]/d[O]
  // d()/d[CO]
  dqdci = -k_r * sc[4];
  J[246] -= dqdci; // dwdot[O2]/d[CO]
  J[247] += dqdci; // dwdot[O]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[258] -= dqdci; // dwdot[C]/d[CO]
  // d()/d[C]
  dqdci = +k_f * sc[3];
  J[408] -= dqdci; // dwdot[O2]/d[C]
  J[409] += dqdci; // dwdot[O]/d[C]
  J[414] += dqdci; // dwdot[CO]/d[C]
  J[420] -= dqdci; // dwdot[C]/d[C]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[706] += dqdT; // dwdot[O]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[717] -= dqdT; // dwdot[C]/dT

  // reaction 57: CH3 + O2 <=> CH3O + O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[12];
  k_f = 7546000 * exp(-(14251.0959755897) * invT);
  dlnkfdT = (14251.0959755897) * invT2;
  // reverse
  phi_r = sc[4] * sc[18];
  Kc = exp(g_RT[3] - g_RT[4] + g_RT[12] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[12]) + (h_RT[4] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[4] += q;  // O
  wdot[12] -= q; // CH3
  wdot[18] += q; // CH3O
  // d()/d[O2]
  dqdci = +k_f * sc[12];
  J[84] -= dqdci; // dwdot[O2]/d[O2]
  J[85] += dqdci; // dwdot[O]/d[O2]
  J[93] -= dqdci; // dwdot[CH3]/d[O2]
  J[99] += dqdci; // dwdot[CH3O]/d[O2]
  // d()/d[O]
  dqdci = -k_r * sc[18];
  J[111] -= dqdci; // dwdot[O2]/d[O]
  J[112] += dqdci; // dwdot[O]/d[O]
  J[120] -= dqdci; // dwdot[CH3]/d[O]
  J[126] += dqdci; // dwdot[CH3O]/d[O]
  // d()/d[CH3]
  dqdci = +k_f * sc[3];
  J[327] -= dqdci; // dwdot[O2]/d[CH3]
  J[328] += dqdci; // dwdot[O]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[342] += dqdci; // dwdot[CH3O]/d[CH3]
  // d()/d[CH3O]
  dqdci = -k_r * sc[4];
  J[489] -= dqdci; // dwdot[O2]/d[CH3O]
  J[490] += dqdci; // dwdot[O]/d[CH3O]
  J[498] -= dqdci; // dwdot[CH3]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[706] += dqdT; // dwdot[O]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 58: CH3 + O2 <=> CH2O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[12];
  k_f = 2.641e-06 * exp(3.283 * tc[0] - (4078.57107634725) * invT);
  dlnkfdT = 3.283 * invT + (4078.57107634725) * invT2;
  // reverse
  phi_r = sc[6] * sc[20];
  Kc = exp(g_RT[3] - g_RT[6] + g_RT[12] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[12]) + (h_RT[6] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[6] += q;  // OH
  wdot[12] -= q; // CH3
  wdot[20] += q; // CH2O
  // d()/d[O2]
  dqdci = +k_f * sc[12];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[87] += dqdci;  // dwdot[OH]/d[O2]
  J[93] -= dqdci;  // dwdot[CH3]/d[O2]
  J[101] += dqdci; // dwdot[CH2O]/d[O2]
  // d()/d[OH]
  dqdci = -k_r * sc[20];
  J[165] -= dqdci; // dwdot[O2]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH3]
  dqdci = +k_f * sc[3];
  J[327] -= dqdci; // dwdot[O2]/d[CH3]
  J[330] += dqdci; // dwdot[OH]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[344] += dqdci; // dwdot[CH2O]/d[CH3]
  // d()/d[CH2O]
  dqdci = -k_r * sc[6];
  J[543] -= dqdci; // dwdot[O2]/d[CH2O]
  J[546] += dqdci; // dwdot[OH]/d[CH2O]
  J[552] -= dqdci; // dwdot[CH3]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 59: CH3 + O <=> CH2O + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[12];
  k_f = 55400000 * exp(0.05 * tc[0] - (-68.4374665494419) * invT);
  dlnkfdT = 0.05 * invT + (-68.4374665494419) * invT2;
  // reverse
  phi_r = sc[2] * sc[20];
  Kc = exp(-g_RT[2] + g_RT[4] + g_RT[12] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[12]) + (h_RT[2] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[4] -= q;  // O
  wdot[12] -= q; // CH3
  wdot[20] += q; // CH2O
  // d()/d[H]
  dqdci = -k_r * sc[20];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[58] -= dqdci; // dwdot[O]/d[H]
  J[66] -= dqdci; // dwdot[CH3]/d[H]
  J[74] += dqdci; // dwdot[CH2O]/d[H]
  // d()/d[O]
  dqdci = +k_f * sc[12];
  J[110] += dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[120] -= dqdci; // dwdot[CH3]/d[O]
  J[128] += dqdci; // dwdot[CH2O]/d[O]
  // d()/d[CH3]
  dqdci = +k_f * sc[4];
  J[326] += dqdci; // dwdot[H]/d[CH3]
  J[328] -= dqdci; // dwdot[O]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[344] += dqdci; // dwdot[CH2O]/d[CH3]
  // d()/d[CH2O]
  dqdci = -k_r * sc[2];
  J[542] += dqdci; // dwdot[H]/d[CH2O]
  J[544] -= dqdci; // dwdot[O]/d[CH2O]
  J[552] -= dqdci; // dwdot[CH3]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 60: CH3 + OH <=> CH2(S) + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[12];
  k_f = 4.788e+17 * exp(-3.155 * tc[0] - (3524.02631063045) * invT);
  dlnkfdT = -3.155 * invT + (3524.02631063045) * invT2;
  // reverse
  phi_r = sc[5] * sc[14];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[12] - g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[5] + h_RT[14]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[12] -= q; // CH3
  wdot[14] += q; // CH2(S)
  // d()/d[H2O]
  dqdci = -k_r * sc[14];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[147] -= dqdci; // dwdot[CH3]/d[H2O]
  J[149] += dqdci; // dwdot[CH2(S)]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[12];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[176] += dqdci; // dwdot[CH2(S)]/d[OH]
  // d()/d[CH3]
  dqdci = +k_f * sc[6];
  J[329] += dqdci; // dwdot[H2O]/d[CH3]
  J[330] -= dqdci; // dwdot[OH]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[338] += dqdci; // dwdot[CH2(S)]/d[CH3]
  // d()/d[CH2(S)]
  dqdci = -k_r * sc[5];
  J[383] += dqdci; // dwdot[H2O]/d[CH2(S)]
  J[384] -= dqdci; // dwdot[OH]/d[CH2(S)]
  J[390] -= dqdci; // dwdot[CH3]/d[CH2(S)]
  J[392] += dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[716] += dqdT; // dwdot[CH2(S)]/dT

  // reaction 61: CH3 + OH <=> CH2O + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[12];
  k_f = 5374 * exp(0.287 * tc[0] - (140.900666425322) * invT);
  dlnkfdT = 0.287 * invT + (140.900666425322) * invT2;
  // reverse
  phi_r = sc[1] * sc[20];
  Kc = exp(-g_RT[1] + g_RT[6] + g_RT[12] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[1] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[6] -= q;  // OH
  wdot[12] -= q; // CH3
  wdot[20] += q; // CH2O
  // d()/d[H2]
  dqdci = -k_r * sc[20];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[33] -= dqdci; // dwdot[OH]/d[H2]
  J[39] -= dqdci; // dwdot[CH3]/d[H2]
  J[47] += dqdci; // dwdot[CH2O]/d[H2]
  // d()/d[OH]
  dqdci = +k_f * sc[12];
  J[163] += dqdci; // dwdot[H2]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH3]
  dqdci = +k_f * sc[6];
  J[325] += dqdci; // dwdot[H2]/d[CH3]
  J[330] -= dqdci; // dwdot[OH]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[344] += dqdci; // dwdot[CH2O]/d[CH3]
  // d()/d[CH2O]
  dqdci = -k_r * sc[1];
  J[541] += dqdci; // dwdot[H2]/d[CH2O]
  J[546] -= dqdci; // dwdot[OH]/d[CH2O]
  J[552] -= dqdci; // dwdot[CH3]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 62: CH3 + OH <=> CH2OH + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[12];
  k_f = 15250000 * exp(0.134 * tc[0] - (2838.64521180442) * invT);
  dlnkfdT = 0.134 * invT + (2838.64521180442) * invT2;
  // reverse
  phi_r = sc[2] * sc[19];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[12] - g_RT[19]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[2] + h_RT[19]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[12] -= q; // CH3
  wdot[19] += q; // CH2OH
  // d()/d[H]
  dqdci = -k_r * sc[19];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[66] -= dqdci; // dwdot[CH3]/d[H]
  J[73] += dqdci; // dwdot[CH2OH]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[12];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[181] += dqdci; // dwdot[CH2OH]/d[OH]
  // d()/d[CH3]
  dqdci = +k_f * sc[6];
  J[326] += dqdci; // dwdot[H]/d[CH3]
  J[330] -= dqdci; // dwdot[OH]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[343] += dqdci; // dwdot[CH2OH]/d[CH3]
  // d()/d[CH2OH]
  dqdci = -k_r * sc[2];
  J[515] += dqdci; // dwdot[H]/d[CH2OH]
  J[519] -= dqdci; // dwdot[OH]/d[CH2OH]
  J[525] -= dqdci; // dwdot[CH3]/d[CH2OH]
  J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT

  // reaction 63: CH3 + OH <=> CH3O + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[12];
  k_f = 1798 * exp(0.965 * tc[0] - (6068.79298960492) * invT);
  dlnkfdT = 0.965 * invT + (6068.79298960492) * invT2;
  // reverse
  phi_r = sc[2] * sc[18];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[12] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[2] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[12] -= q; // CH3
  wdot[18] += q; // CH3O
  // d()/d[H]
  dqdci = -k_r * sc[18];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[66] -= dqdci; // dwdot[CH3]/d[H]
  J[72] += dqdci; // dwdot[CH3O]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[12];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[180] += dqdci; // dwdot[CH3O]/d[OH]
  // d()/d[CH3]
  dqdci = +k_f * sc[6];
  J[326] += dqdci; // dwdot[H]/d[CH3]
  J[330] -= dqdci; // dwdot[OH]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[342] += dqdci; // dwdot[CH3O]/d[CH3]
  // d()/d[CH3O]
  dqdci = -k_r * sc[2];
  J[488] += dqdci; // dwdot[H]/d[CH3O]
  J[492] -= dqdci; // dwdot[OH]/d[CH3O]
  J[498] -= dqdci; // dwdot[CH3]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 64: CH3 + OH <=> CH2 + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[12];
  k_f = 0.04293 * exp(2.568 * tc[0] - (2011.75958655411) * invT);
  dlnkfdT = 2.568 * invT + (2011.75958655411) * invT2;
  // reverse
  phi_r = sc[5] * sc[13];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[12] - g_RT[13]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[5] + h_RT[13]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[12] -= q; // CH3
  wdot[13] += q; // CH2
  // d()/d[H2O]
  dqdci = -k_r * sc[13];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[147] -= dqdci; // dwdot[CH3]/d[H2O]
  J[148] += dqdci; // dwdot[CH2]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[12];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[175] += dqdci; // dwdot[CH2]/d[OH]
  // d()/d[CH3]
  dqdci = +k_f * sc[6];
  J[329] += dqdci; // dwdot[H2O]/d[CH3]
  J[330] -= dqdci; // dwdot[OH]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[337] += dqdci; // dwdot[CH2]/d[CH3]
  // d()/d[CH2]
  dqdci = -k_r * sc[5];
  J[356] += dqdci; // dwdot[H2O]/d[CH2]
  J[357] -= dqdci; // dwdot[OH]/d[CH2]
  J[363] -= dqdci; // dwdot[CH3]/d[CH2]
  J[364] += dqdci; // dwdot[CH2]/d[CH2]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[715] += dqdT; // dwdot[CH2]/dT

  // reaction 65: CH3 + HO2 <=> CH3O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[12];
  k_f = 1000000 * exp(0.269 * tc[0] - (-345.961457740745) * invT);
  dlnkfdT = 0.269 * invT + (-345.961457740745) * invT2;
  // reverse
  phi_r = sc[6] * sc[18];
  Kc = exp(-g_RT[6] + g_RT[8] + g_RT[12] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[6] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[6] += q;  // OH
  wdot[8] -= q;  // HO2
  wdot[12] -= q; // CH3
  wdot[18] += q; // CH3O
  // d()/d[OH]
  dqdci = -k_r * sc[18];
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[170] -= dqdci; // dwdot[HO2]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[180] += dqdci; // dwdot[CH3O]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[12];
  J[222] += dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[228] -= dqdci; // dwdot[CH3]/d[HO2]
  J[234] += dqdci; // dwdot[CH3O]/d[HO2]
  // d()/d[CH3]
  dqdci = +k_f * sc[8];
  J[330] += dqdci; // dwdot[OH]/d[CH3]
  J[332] -= dqdci; // dwdot[HO2]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[342] += dqdci; // dwdot[CH3O]/d[CH3]
  // d()/d[CH3O]
  dqdci = -k_r * sc[6];
  J[492] += dqdci; // dwdot[OH]/d[CH3O]
  J[494] -= dqdci; // dwdot[HO2]/d[CH3O]
  J[498] -= dqdci; // dwdot[CH3]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 68: CH3OH + H <=> CH3O + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[17];
  k_f = 0.199 * exp(2.56 * tc[0] - (5183.13165778861) * invT);
  dlnkfdT = 2.56 * invT + (5183.13165778861) * invT2;
  // reverse
  phi_r = sc[1] * sc[18];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[17] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[1] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[17] -= q; // CH3OH
  wdot[18] += q; // CH3O
  // d()/d[H2]
  dqdci = -k_r * sc[18];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[44] -= dqdci; // dwdot[CH3OH]/d[H2]
  J[45] += dqdci; // dwdot[CH3O]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[17];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[71] -= dqdci; // dwdot[CH3OH]/d[H]
  J[72] += dqdci; // dwdot[CH3O]/d[H]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[2];
  J[460] += dqdci; // dwdot[H2]/d[CH3OH]
  J[461] -= dqdci; // dwdot[H]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] += dqdci; // dwdot[CH3O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = -k_r * sc[1];
  J[487] += dqdci; // dwdot[H2]/d[CH3O]
  J[488] -= dqdci; // dwdot[H]/d[CH3O]
  J[503] -= dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 69: CH3OH + H <=> CH2OH + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[17];
  k_f = 0.307 * exp(2.55 * tc[0] - (2737.49866197768) * invT);
  dlnkfdT = 2.55 * invT + (2737.49866197768) * invT2;
  // reverse
  phi_r = sc[1] * sc[19];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[17] - g_RT[19]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[1] + h_RT[19]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[17] -= q; // CH3OH
  wdot[19] += q; // CH2OH
  // d()/d[H2]
  dqdci = -k_r * sc[19];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[44] -= dqdci; // dwdot[CH3OH]/d[H2]
  J[46] += dqdci; // dwdot[CH2OH]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[17];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[71] -= dqdci; // dwdot[CH3OH]/d[H]
  J[73] += dqdci; // dwdot[CH2OH]/d[H]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[2];
  J[460] += dqdci; // dwdot[H2]/d[CH3OH]
  J[461] -= dqdci; // dwdot[H]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[478] += dqdci; // dwdot[CH2OH]/d[CH3OH]
  // d()/d[CH2OH]
  dqdci = -k_r * sc[1];
  J[514] += dqdci; // dwdot[H2]/d[CH2OH]
  J[515] -= dqdci; // dwdot[H]/d[CH2OH]
  J[530] -= dqdci; // dwdot[CH3OH]/d[CH2OH]
  J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT

  // reaction 70: CH3OH + O <=> CH3O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[17];
  k_f = 0.0388 * exp(2.5 * tc[0] - (1549.90733067854) * invT);
  dlnkfdT = 2.5 * invT + (1549.90733067854) * invT2;
  // reverse
  phi_r = sc[6] * sc[18];
  Kc = exp(g_RT[4] - g_RT[6] + g_RT[17] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[17]) + (h_RT[6] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[6] += q;  // OH
  wdot[17] -= q; // CH3OH
  wdot[18] += q; // CH3O
  // d()/d[O]
  dqdci = +k_f * sc[17];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[125] -= dqdci; // dwdot[CH3OH]/d[O]
  J[126] += dqdci; // dwdot[CH3O]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[18];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[179] -= dqdci; // dwdot[CH3OH]/d[OH]
  J[180] += dqdci; // dwdot[CH3O]/d[OH]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[4];
  J[463] -= dqdci; // dwdot[O]/d[CH3OH]
  J[465] += dqdci; // dwdot[OH]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] += dqdci; // dwdot[CH3O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = -k_r * sc[6];
  J[490] -= dqdci; // dwdot[O]/d[CH3O]
  J[492] += dqdci; // dwdot[OH]/d[CH3O]
  J[503] -= dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 71: CH3OH + O <=> CH2OH + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[17];
  k_f = 0.388 * exp(2.5 * tc[0] - (1549.90733067854) * invT);
  dlnkfdT = 2.5 * invT + (1549.90733067854) * invT2;
  // reverse
  phi_r = sc[6] * sc[19];
  Kc = exp(g_RT[4] - g_RT[6] + g_RT[17] - g_RT[19]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[17]) + (h_RT[6] + h_RT[19]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[6] += q;  // OH
  wdot[17] -= q; // CH3OH
  wdot[19] += q; // CH2OH
  // d()/d[O]
  dqdci = +k_f * sc[17];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[125] -= dqdci; // dwdot[CH3OH]/d[O]
  J[127] += dqdci; // dwdot[CH2OH]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[19];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[179] -= dqdci; // dwdot[CH3OH]/d[OH]
  J[181] += dqdci; // dwdot[CH2OH]/d[OH]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[4];
  J[463] -= dqdci; // dwdot[O]/d[CH3OH]
  J[465] += dqdci; // dwdot[OH]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[478] += dqdci; // dwdot[CH2OH]/d[CH3OH]
  // d()/d[CH2OH]
  dqdci = -k_r * sc[6];
  J[517] -= dqdci; // dwdot[O]/d[CH2OH]
  J[519] += dqdci; // dwdot[OH]/d[CH2OH]
  J[530] -= dqdci; // dwdot[CH3OH]/d[CH2OH]
  J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT

  // reaction 72: CH3OH + OH <=> CH3O + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[17];
  k_f = 0.00015 * exp(3.03 * tc[0] - (-383.954316009001) * invT);
  dlnkfdT = 3.03 * invT + (-383.954316009001) * invT2;
  // reverse
  phi_r = sc[5] * sc[18];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[17] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[17]) + (h_RT[5] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[17] -= q; // CH3OH
  wdot[18] += q; // CH3O
  // d()/d[H2O]
  dqdci = -k_r * sc[18];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[152] -= dqdci; // dwdot[CH3OH]/d[H2O]
  J[153] += dqdci; // dwdot[CH3O]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[17];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[179] -= dqdci; // dwdot[CH3OH]/d[OH]
  J[180] += dqdci; // dwdot[CH3O]/d[OH]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[6];
  J[464] += dqdci; // dwdot[H2O]/d[CH3OH]
  J[465] -= dqdci; // dwdot[OH]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] += dqdci; // dwdot[CH3O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = -k_r * sc[5];
  J[491] += dqdci; // dwdot[H2O]/d[CH3O]
  J[492] -= dqdci; // dwdot[OH]/d[CH3O]
  J[503] -= dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 73: CH3OH + OH <=> CH2OH + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[17];
  k_f = 0.0308 * exp(2.65 * tc[0] - (-405.944884304667) * invT);
  dlnkfdT = 2.65 * invT + (-405.944884304667) * invT2;
  // reverse
  phi_r = sc[5] * sc[19];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[17] - g_RT[19]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[17]) + (h_RT[5] + h_RT[19]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[17] -= q; // CH3OH
  wdot[19] += q; // CH2OH
  // d()/d[H2O]
  dqdci = -k_r * sc[19];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[152] -= dqdci; // dwdot[CH3OH]/d[H2O]
  J[154] += dqdci; // dwdot[CH2OH]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[17];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[179] -= dqdci; // dwdot[CH3OH]/d[OH]
  J[181] += dqdci; // dwdot[CH2OH]/d[OH]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[6];
  J[464] += dqdci; // dwdot[H2O]/d[CH3OH]
  J[465] -= dqdci; // dwdot[OH]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[478] += dqdci; // dwdot[CH2OH]/d[CH3OH]
  // d()/d[CH2OH]
  dqdci = -k_r * sc[5];
  J[518] += dqdci; // dwdot[H2O]/d[CH2OH]
  J[519] -= dqdci; // dwdot[OH]/d[CH2OH]
  J[530] -= dqdci; // dwdot[CH3OH]/d[CH2OH]
  J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT

  // reaction 74: CH3OH + O2 <=> CH3O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[17];
  k_f = 0.0358 * exp(2.27 * tc[0] - (21519.8091048059) * invT);
  dlnkfdT = 2.27 * invT + (21519.8091048059) * invT2;
  // reverse
  phi_r = sc[8] * sc[18];
  Kc = exp(g_RT[3] - g_RT[8] + g_RT[17] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[17]) + (h_RT[8] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[8] += q;  // HO2
  wdot[17] -= q; // CH3OH
  wdot[18] += q; // CH3O
  // d()/d[O2]
  dqdci = +k_f * sc[17];
  J[84] -= dqdci; // dwdot[O2]/d[O2]
  J[89] += dqdci; // dwdot[HO2]/d[O2]
  J[98] -= dqdci; // dwdot[CH3OH]/d[O2]
  J[99] += dqdci; // dwdot[CH3O]/d[O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[18];
  J[219] -= dqdci; // dwdot[O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  J[233] -= dqdci; // dwdot[CH3OH]/d[HO2]
  J[234] += dqdci; // dwdot[CH3O]/d[HO2]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[3];
  J[462] -= dqdci; // dwdot[O2]/d[CH3OH]
  J[467] += dqdci; // dwdot[HO2]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] += dqdci; // dwdot[CH3O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = -k_r * sc[8];
  J[489] -= dqdci; // dwdot[O2]/d[CH3O]
  J[494] += dqdci; // dwdot[HO2]/d[CH3O]
  J[503] -= dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 75: CH3OH + HO2 <=> CH3O + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[17];
  k_f = 1220000 * exp(-(10099.9107343668) * invT);
  dlnkfdT = (10099.9107343668) * invT2;
  // reverse
  phi_r = sc[7] * sc[18];
  Kc = exp(-g_RT[7] + g_RT[8] + g_RT[17] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[17]) + (h_RT[7] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[17] -= q; // CH3OH
  wdot[18] += q; // CH3O
  // d()/d[H2O2]
  dqdci = -k_r * sc[18];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[206] -= dqdci; // dwdot[CH3OH]/d[H2O2]
  J[207] += dqdci; // dwdot[CH3O]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[17];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[233] -= dqdci; // dwdot[CH3OH]/d[HO2]
  J[234] += dqdci; // dwdot[CH3O]/d[HO2]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[8];
  J[466] += dqdci; // dwdot[H2O2]/d[CH3OH]
  J[467] -= dqdci; // dwdot[HO2]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] += dqdci; // dwdot[CH3O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = -k_r * sc[7];
  J[493] += dqdci; // dwdot[H2O2]/d[CH3O]
  J[494] -= dqdci; // dwdot[HO2]/d[CH3O]
  J[503] -= dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 76: CH3OH + HO2 <=> CH2OH + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[17];
  k_f = 32600000 * exp(-(9451.5160604774) * invT);
  dlnkfdT = (9451.5160604774) * invT2;
  // reverse
  phi_r = sc[7] * sc[19];
  Kc = exp(-g_RT[7] + g_RT[8] + g_RT[17] - g_RT[19]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[17]) + (h_RT[7] + h_RT[19]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[17] -= q; // CH3OH
  wdot[19] += q; // CH2OH
  // d()/d[H2O2]
  dqdci = -k_r * sc[19];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[206] -= dqdci; // dwdot[CH3OH]/d[H2O2]
  J[208] += dqdci; // dwdot[CH2OH]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[17];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[233] -= dqdci; // dwdot[CH3OH]/d[HO2]
  J[235] += dqdci; // dwdot[CH2OH]/d[HO2]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[8];
  J[466] += dqdci; // dwdot[H2O2]/d[CH3OH]
  J[467] -= dqdci; // dwdot[HO2]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[478] += dqdci; // dwdot[CH2OH]/d[CH3OH]
  // d()/d[CH2OH]
  dqdci = -k_r * sc[7];
  J[520] += dqdci; // dwdot[H2O2]/d[CH2OH]
  J[521] -= dqdci; // dwdot[HO2]/d[CH2OH]
  J[530] -= dqdci; // dwdot[CH3OH]/d[CH2OH]
  J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT

  // reaction 77: CH3 + CH3OH <=> CH2OH + CH4
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[12] * sc[17];
  k_f = 2.13e-07 * exp(3.953 * tc[0] - (3550.24389891888) * invT);
  dlnkfdT = 3.953 * invT + (3550.24389891888) * invT2;
  // reverse
  phi_r = sc[11] * sc[19];
  Kc = exp(-g_RT[11] + g_RT[12] + g_RT[17] - g_RT[19]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[12] + h_RT[17]) + (h_RT[11] + h_RT[19]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  wdot[17] -= q; // CH3OH
  wdot[19] += q; // CH2OH
  // d()/d[CH4]
  dqdci = -k_r * sc[19];
  J[308] += dqdci; // dwdot[CH4]/d[CH4]
  J[309] -= dqdci; // dwdot[CH3]/d[CH4]
  J[314] -= dqdci; // dwdot[CH3OH]/d[CH4]
  J[316] += dqdci; // dwdot[CH2OH]/d[CH4]
  // d()/d[CH3]
  dqdci = +k_f * sc[17];
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[341] -= dqdci; // dwdot[CH3OH]/d[CH3]
  J[343] += dqdci; // dwdot[CH2OH]/d[CH3]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[12];
  J[470] += dqdci; // dwdot[CH4]/d[CH3OH]
  J[471] -= dqdci; // dwdot[CH3]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[478] += dqdci; // dwdot[CH2OH]/d[CH3OH]
  // d()/d[CH2OH]
  dqdci = -k_r * sc[11];
  J[524] += dqdci; // dwdot[CH4]/d[CH2OH]
  J[525] -= dqdci; // dwdot[CH3]/d[CH2OH]
  J[530] -= dqdci; // dwdot[CH3OH]/d[CH2OH]
  J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
  // d()/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT

  // reaction 78: CH3 + CH3OH <=> CH3O + CH4
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[12] * sc[17];
  k_f = 0.00322 * exp(2.425 * tc[0] - (4317.34738427159) * invT);
  dlnkfdT = 2.425 * invT + (4317.34738427159) * invT2;
  // reverse
  phi_r = sc[11] * sc[18];
  Kc = exp(-g_RT[11] + g_RT[12] + g_RT[17] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[12] + h_RT[17]) + (h_RT[11] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  wdot[17] -= q; // CH3OH
  wdot[18] += q; // CH3O
  // d()/d[CH4]
  dqdci = -k_r * sc[18];
  J[308] += dqdci; // dwdot[CH4]/d[CH4]
  J[309] -= dqdci; // dwdot[CH3]/d[CH4]
  J[314] -= dqdci; // dwdot[CH3OH]/d[CH4]
  J[315] += dqdci; // dwdot[CH3O]/d[CH4]
  // d()/d[CH3]
  dqdci = +k_f * sc[17];
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[341] -= dqdci; // dwdot[CH3OH]/d[CH3]
  J[342] += dqdci; // dwdot[CH3O]/d[CH3]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[12];
  J[470] += dqdci; // dwdot[CH4]/d[CH3OH]
  J[471] -= dqdci; // dwdot[CH3]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] += dqdci; // dwdot[CH3O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = -k_r * sc[11];
  J[497] += dqdci; // dwdot[CH4]/d[CH3O]
  J[498] -= dqdci; // dwdot[CH3]/d[CH3O]
  J[503] -= dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 79: CH2OH + O2 <=> CH2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[19];
  k_f = 1510000000 * exp(-1 * tc[0]);
  dlnkfdT = -1 * invT;
  // reverse
  phi_r = sc[8] * sc[20];
  Kc = exp(g_RT[3] - g_RT[8] + g_RT[19] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[19]) + (h_RT[8] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[8] += q;  // HO2
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[O2]
  dqdci = +k_f * sc[19];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[89] += dqdci;  // dwdot[HO2]/d[O2]
  J[100] -= dqdci; // dwdot[CH2OH]/d[O2]
  J[101] += dqdci; // dwdot[CH2O]/d[O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[20];
  J[219] -= dqdci; // dwdot[O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  J[235] -= dqdci; // dwdot[CH2OH]/d[HO2]
  J[236] += dqdci; // dwdot[CH2O]/d[HO2]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[3];
  J[516] -= dqdci; // dwdot[O2]/d[CH2OH]
  J[521] += dqdci; // dwdot[HO2]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[8];
  J[543] -= dqdci; // dwdot[O2]/d[CH2O]
  J[548] += dqdci; // dwdot[HO2]/d[CH2O]
  J[559] -= dqdci; // dwdot[CH2OH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 80: CH2OH + O2 <=> CH2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[19];
  k_f = 241000000 * exp(-(2524.63801234228) * invT);
  dlnkfdT = (2524.63801234228) * invT2;
  // reverse
  phi_r = sc[8] * sc[20];
  Kc = exp(g_RT[3] - g_RT[8] + g_RT[19] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[19]) + (h_RT[8] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[8] += q;  // HO2
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[O2]
  dqdci = +k_f * sc[19];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[89] += dqdci;  // dwdot[HO2]/d[O2]
  J[100] -= dqdci; // dwdot[CH2OH]/d[O2]
  J[101] += dqdci; // dwdot[CH2O]/d[O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[20];
  J[219] -= dqdci; // dwdot[O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  J[235] -= dqdci; // dwdot[CH2OH]/d[HO2]
  J[236] += dqdci; // dwdot[CH2O]/d[HO2]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[3];
  J[516] -= dqdci; // dwdot[O2]/d[CH2OH]
  J[521] += dqdci; // dwdot[HO2]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[8];
  J[543] -= dqdci; // dwdot[O2]/d[CH2O]
  J[548] += dqdci; // dwdot[HO2]/d[CH2O]
  J[559] -= dqdci; // dwdot[CH2OH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 81: CH2OH + H <=> CH2O + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[19];
  k_f = 6000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[20];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[19] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[19]) + (h_RT[1] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[H2]
  dqdci = -k_r * sc[20];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[46] -= dqdci; // dwdot[CH2OH]/d[H2]
  J[47] += dqdci; // dwdot[CH2O]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[19];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[73] -= dqdci; // dwdot[CH2OH]/d[H]
  J[74] += dqdci; // dwdot[CH2O]/d[H]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[2];
  J[514] += dqdci; // dwdot[H2]/d[CH2OH]
  J[515] -= dqdci; // dwdot[H]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[1];
  J[541] += dqdci; // dwdot[H2]/d[CH2O]
  J[542] -= dqdci; // dwdot[H]/d[CH2O]
  J[559] -= dqdci; // dwdot[CH2OH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 82: CH2OH + HO2 <=> CH2O + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[19];
  k_f = 12000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[7] * sc[20];
  Kc = exp(-g_RT[7] + g_RT[8] + g_RT[19] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[19]) + (h_RT[7] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[H2O2]
  dqdci = -k_r * sc[20];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[208] -= dqdci; // dwdot[CH2OH]/d[H2O2]
  J[209] += dqdci; // dwdot[CH2O]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[19];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[235] -= dqdci; // dwdot[CH2OH]/d[HO2]
  J[236] += dqdci; // dwdot[CH2O]/d[HO2]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[8];
  J[520] += dqdci; // dwdot[H2O2]/d[CH2OH]
  J[521] -= dqdci; // dwdot[HO2]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[7];
  J[547] += dqdci; // dwdot[H2O2]/d[CH2O]
  J[548] -= dqdci; // dwdot[HO2]/d[CH2O]
  J[559] -= dqdci; // dwdot[CH2OH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 83: CH2OH + HCO <=> 2 CH2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[19] * sc[21];
  k_f = 180000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = (sc[20] * sc[20]);
  Kc = exp(g_RT[19] - 2.000000 * g_RT[20] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[19] + h_RT[21]) + (2.000000 * h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[19] -= q;     // CH2OH
  wdot[20] += 2 * q; // CH2O
  wdot[21] -= q;     // HCO
  // d()/d[CH2OH]
  dqdci = +k_f * sc[21];
  J[532] -= dqdci;     // dwdot[CH2OH]/d[CH2OH]
  J[533] += 2 * dqdci; // dwdot[CH2O]/d[CH2OH]
  J[534] -= dqdci;     // dwdot[HCO]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * 2.000000 * sc[20];
  J[559] -= dqdci;     // dwdot[CH2OH]/d[CH2O]
  J[560] += 2 * dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] -= dqdci;     // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = +k_f * sc[19];
  J[586] -= dqdci;     // dwdot[CH2OH]/d[HCO]
  J[587] += 2 * dqdci; // dwdot[CH2O]/d[HCO]
  J[588] -= dqdci;     // dwdot[HCO]/d[HCO]
  // d()/dT
  J[721] -= dqdT;     // dwdot[CH2OH]/dT
  J[722] += 2 * dqdT; // dwdot[CH2O]/dT
  J[723] -= dqdT;     // dwdot[HCO]/dT

  // reaction 84: CH2OH + OH <=> CH2O + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[19];
  k_f = 24000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[5] * sc[20];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[19] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[19]) + (h_RT[5] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[H2O]
  dqdci = -k_r * sc[20];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[154] -= dqdci; // dwdot[CH2OH]/d[H2O]
  J[155] += dqdci; // dwdot[CH2O]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[19];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[181] -= dqdci; // dwdot[CH2OH]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[6];
  J[518] += dqdci; // dwdot[H2O]/d[CH2OH]
  J[519] -= dqdci; // dwdot[OH]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[5];
  J[545] += dqdci; // dwdot[H2O]/d[CH2O]
  J[546] -= dqdci; // dwdot[OH]/d[CH2O]
  J[559] -= dqdci; // dwdot[CH2OH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 85: CH2OH + O <=> CH2O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[19];
  k_f = 42000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[6] * sc[20];
  Kc = exp(g_RT[4] - g_RT[6] + g_RT[19] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[19]) + (h_RT[6] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[6] += q;  // OH
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[O]
  dqdci = +k_f * sc[19];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[127] -= dqdci; // dwdot[CH2OH]/d[O]
  J[128] += dqdci; // dwdot[CH2O]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[20];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[181] -= dqdci; // dwdot[CH2OH]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[4];
  J[517] -= dqdci; // dwdot[O]/d[CH2OH]
  J[519] += dqdci; // dwdot[OH]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[6];
  J[544] -= dqdci; // dwdot[O]/d[CH2O]
  J[546] += dqdci; // dwdot[OH]/d[CH2O]
  J[559] -= dqdci; // dwdot[CH2OH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 86: CH3O + O2 <=> CH2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[18];
  k_f = 4.38e-25 * exp(9.5 * tc[0] - (-2768.19487859176) * invT);
  dlnkfdT = 9.5 * invT + (-2768.19487859176) * invT2;
  // reverse
  phi_r = sc[8] * sc[20];
  Kc = exp(g_RT[3] - g_RT[8] + g_RT[18] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[18]) + (h_RT[8] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[8] += q;  // HO2
  wdot[18] -= q; // CH3O
  wdot[20] += q; // CH2O
  // d()/d[O2]
  dqdci = +k_f * sc[18];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[89] += dqdci;  // dwdot[HO2]/d[O2]
  J[99] -= dqdci;  // dwdot[CH3O]/d[O2]
  J[101] += dqdci; // dwdot[CH2O]/d[O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[20];
  J[219] -= dqdci; // dwdot[O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  J[234] -= dqdci; // dwdot[CH3O]/d[HO2]
  J[236] += dqdci; // dwdot[CH2O]/d[HO2]
  // d()/d[CH3O]
  dqdci = +k_f * sc[3];
  J[489] -= dqdci; // dwdot[O2]/d[CH3O]
  J[494] += dqdci; // dwdot[HO2]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] += dqdci; // dwdot[CH2O]/d[CH3O]
  // d()/d[CH2O]
  dqdci = -k_r * sc[8];
  J[543] -= dqdci; // dwdot[O2]/d[CH2O]
  J[548] += dqdci; // dwdot[HO2]/d[CH2O]
  J[558] -= dqdci; // dwdot[CH3O]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 87: CH3O + H <=> CH2O + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[18];
  k_f = 20000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[20];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[18] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[18]) + (h_RT[1] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[18] -= q; // CH3O
  wdot[20] += q; // CH2O
  // d()/d[H2]
  dqdci = -k_r * sc[20];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[45] -= dqdci; // dwdot[CH3O]/d[H2]
  J[47] += dqdci; // dwdot[CH2O]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[18];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[72] -= dqdci; // dwdot[CH3O]/d[H]
  J[74] += dqdci; // dwdot[CH2O]/d[H]
  // d()/d[CH3O]
  dqdci = +k_f * sc[2];
  J[487] += dqdci; // dwdot[H2]/d[CH3O]
  J[488] -= dqdci; // dwdot[H]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] += dqdci; // dwdot[CH2O]/d[CH3O]
  // d()/d[CH2O]
  dqdci = -k_r * sc[1];
  J[541] += dqdci; // dwdot[H2]/d[CH2O]
  J[542] -= dqdci; // dwdot[H]/d[CH2O]
  J[558] -= dqdci; // dwdot[CH3O]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 88: CH3O + HO2 <=> CH2O + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[18];
  k_f = 301000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[7] * sc[20];
  Kc = exp(-g_RT[7] + g_RT[8] + g_RT[18] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[18]) + (h_RT[7] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[18] -= q; // CH3O
  wdot[20] += q; // CH2O
  // d()/d[H2O2]
  dqdci = -k_r * sc[20];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[207] -= dqdci; // dwdot[CH3O]/d[H2O2]
  J[209] += dqdci; // dwdot[CH2O]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[18];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[234] -= dqdci; // dwdot[CH3O]/d[HO2]
  J[236] += dqdci; // dwdot[CH2O]/d[HO2]
  // d()/d[CH3O]
  dqdci = +k_f * sc[8];
  J[493] += dqdci; // dwdot[H2O2]/d[CH3O]
  J[494] -= dqdci; // dwdot[HO2]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] += dqdci; // dwdot[CH2O]/d[CH3O]
  // d()/d[CH2O]
  dqdci = -k_r * sc[7];
  J[547] += dqdci; // dwdot[H2O2]/d[CH2O]
  J[548] -= dqdci; // dwdot[HO2]/d[CH2O]
  J[558] -= dqdci; // dwdot[CH3O]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 89: CH3 + CH3O <=> CH2O + CH4
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[12] * sc[18];
  k_f = 12000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[11] * sc[20];
  Kc = exp(-g_RT[11] + g_RT[12] + g_RT[18] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[12] + h_RT[18]) + (h_RT[11] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  wdot[18] -= q; // CH3O
  wdot[20] += q; // CH2O
  // d()/d[CH4]
  dqdci = -k_r * sc[20];
  J[308] += dqdci; // dwdot[CH4]/d[CH4]
  J[309] -= dqdci; // dwdot[CH3]/d[CH4]
  J[315] -= dqdci; // dwdot[CH3O]/d[CH4]
  J[317] += dqdci; // dwdot[CH2O]/d[CH4]
  // d()/d[CH3]
  dqdci = +k_f * sc[18];
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[342] -= dqdci; // dwdot[CH3O]/d[CH3]
  J[344] += dqdci; // dwdot[CH2O]/d[CH3]
  // d()/d[CH3O]
  dqdci = +k_f * sc[12];
  J[497] += dqdci; // dwdot[CH4]/d[CH3O]
  J[498] -= dqdci; // dwdot[CH3]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] += dqdci; // dwdot[CH2O]/d[CH3O]
  // d()/d[CH2O]
  dqdci = -k_r * sc[11];
  J[551] += dqdci; // dwdot[CH4]/d[CH2O]
  J[552] -= dqdci; // dwdot[CH3]/d[CH2O]
  J[558] -= dqdci; // dwdot[CH3O]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 90: 2 CH3O <=> CH2O + CH3OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[18] * sc[18]);
  k_f = 60300000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[17] * sc[20];
  Kc = exp(-g_RT[17] + 2.000000 * g_RT[18] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[18]) + (h_RT[17] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[17] += q;     // CH3OH
  wdot[18] -= 2 * q; // CH3O
  wdot[20] += q;     // CH2O
  // d()/d[CH3OH]
  dqdci = -k_r * sc[20];
  J[476] += dqdci;      // dwdot[CH3OH]/d[CH3OH]
  J[477] += -2 * dqdci; // dwdot[CH3O]/d[CH3OH]
  J[479] += dqdci;      // dwdot[CH2O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = +k_f * 2.000000 * sc[18];
  J[503] += dqdci;      // dwdot[CH3OH]/d[CH3O]
  J[504] += -2 * dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] += dqdci;      // dwdot[CH2O]/d[CH3O]
  // d()/d[CH2O]
  dqdci = -k_r * sc[17];
  J[557] += dqdci;      // dwdot[CH3OH]/d[CH2O]
  J[558] += -2 * dqdci; // dwdot[CH3O]/d[CH2O]
  J[560] += dqdci;      // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[719] += dqdT;      // dwdot[CH3OH]/dT
  J[720] += -2 * dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT;      // dwdot[CH2O]/dT

  // reaction 91: CH2O + O <=> HCO + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[20];
  k_f = 6260 * exp(1.15 * tc[0] - (1137.26966471867) * invT);
  dlnkfdT = 1.15 * invT + (1137.26966471867) * invT2;
  // reverse
  phi_r = sc[6] * sc[21];
  Kc = exp(g_RT[4] - g_RT[6] + g_RT[20] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[20]) + (h_RT[6] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[6] += q;  // OH
  wdot[20] -= q; // CH2O
  wdot[21] += q; // HCO
  // d()/d[O]
  dqdci = +k_f * sc[20];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[128] -= dqdci; // dwdot[CH2O]/d[O]
  J[129] += dqdci; // dwdot[HCO]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[21];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[182] -= dqdci; // dwdot[CH2O]/d[OH]
  J[183] += dqdci; // dwdot[HCO]/d[OH]
  // d()/d[CH2O]
  dqdci = +k_f * sc[4];
  J[544] -= dqdci; // dwdot[O]/d[CH2O]
  J[546] += dqdci; // dwdot[OH]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] += dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = -k_r * sc[6];
  J[571] -= dqdci; // dwdot[O]/d[HCO]
  J[573] += dqdci; // dwdot[OH]/d[HCO]
  J[587] -= dqdci; // dwdot[CH2O]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 92: CH2O + H <=> H2 + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[20];
  k_f = 57.4 * exp(1.9 * tc[0] - (1378.81366430493) * invT);
  dlnkfdT = 1.9 * invT + (1378.81366430493) * invT2;
  // reverse
  phi_r = sc[1] * sc[21];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[20] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[20]) + (h_RT[1] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[20] -= q; // CH2O
  wdot[21] += q; // HCO
  // d()/d[H2]
  dqdci = -k_r * sc[21];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[47] -= dqdci; // dwdot[CH2O]/d[H2]
  J[48] += dqdci; // dwdot[HCO]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[20];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[74] -= dqdci; // dwdot[CH2O]/d[H]
  J[75] += dqdci; // dwdot[HCO]/d[H]
  // d()/d[CH2O]
  dqdci = +k_f * sc[2];
  J[541] += dqdci; // dwdot[H2]/d[CH2O]
  J[542] -= dqdci; // dwdot[H]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] += dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = -k_r * sc[1];
  J[568] += dqdci; // dwdot[H2]/d[HCO]
  J[569] -= dqdci; // dwdot[H]/d[HCO]
  J[587] -= dqdci; // dwdot[CH2O]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 93: CH2O + OH <=> H2O + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[20];
  k_f = 78.2 * exp(1.63 * tc[0] - (-530.893582423979) * invT);
  dlnkfdT = 1.63 * invT + (-530.893582423979) * invT2;
  // reverse
  phi_r = sc[5] * sc[21];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[20] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[20]) + (h_RT[5] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[20] -= q; // CH2O
  wdot[21] += q; // HCO
  // d()/d[H2O]
  dqdci = -k_r * sc[21];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[155] -= dqdci; // dwdot[CH2O]/d[H2O]
  J[156] += dqdci; // dwdot[HCO]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[20];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[182] -= dqdci; // dwdot[CH2O]/d[OH]
  J[183] += dqdci; // dwdot[HCO]/d[OH]
  // d()/d[CH2O]
  dqdci = +k_f * sc[6];
  J[545] += dqdci; // dwdot[H2O]/d[CH2O]
  J[546] -= dqdci; // dwdot[OH]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] += dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = -k_r * sc[5];
  J[572] += dqdci; // dwdot[H2O]/d[HCO]
  J[573] -= dqdci; // dwdot[OH]/d[HCO]
  J[587] -= dqdci; // dwdot[CH2O]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 94: CH2O + HO2 <=> H2O2 + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[20];
  k_f = 0.0188 * exp(2.7 * tc[0] - (5797.05599007037) * invT);
  dlnkfdT = 2.7 * invT + (5797.05599007037) * invT2;
  // reverse
  phi_r = sc[7] * sc[21];
  Kc = exp(-g_RT[7] + g_RT[8] + g_RT[20] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[20]) + (h_RT[7] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[20] -= q; // CH2O
  wdot[21] += q; // HCO
  // d()/d[H2O2]
  dqdci = -k_r * sc[21];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[209] -= dqdci; // dwdot[CH2O]/d[H2O2]
  J[210] += dqdci; // dwdot[HCO]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[20];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[236] -= dqdci; // dwdot[CH2O]/d[HO2]
  J[237] += dqdci; // dwdot[HCO]/d[HO2]
  // d()/d[CH2O]
  dqdci = +k_f * sc[8];
  J[547] += dqdci; // dwdot[H2O2]/d[CH2O]
  J[548] -= dqdci; // dwdot[HO2]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] += dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = -k_r * sc[7];
  J[574] += dqdci; // dwdot[H2O2]/d[HCO]
  J[575] -= dqdci; // dwdot[HO2]/d[HCO]
  J[587] -= dqdci; // dwdot[CH2O]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 95: CH2O + CH3 <=> CH4 + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[12] * sc[20];
  k_f = 3.83e-05 * exp(3.36 * tc[0] - (2169.87026294995) * invT);
  dlnkfdT = 3.36 * invT + (2169.87026294995) * invT2;
  // reverse
  phi_r = sc[11] * sc[21];
  Kc = exp(-g_RT[11] + g_RT[12] + g_RT[20] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[12] + h_RT[20]) + (h_RT[11] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  wdot[20] -= q; // CH2O
  wdot[21] += q; // HCO
  // d()/d[CH4]
  dqdci = -k_r * sc[21];
  J[308] += dqdci; // dwdot[CH4]/d[CH4]
  J[309] -= dqdci; // dwdot[CH3]/d[CH4]
  J[317] -= dqdci; // dwdot[CH2O]/d[CH4]
  J[318] += dqdci; // dwdot[HCO]/d[CH4]
  // d()/d[CH3]
  dqdci = +k_f * sc[20];
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[344] -= dqdci; // dwdot[CH2O]/d[CH3]
  J[345] += dqdci; // dwdot[HCO]/d[CH3]
  // d()/d[CH2O]
  dqdci = +k_f * sc[12];
  J[551] += dqdci; // dwdot[CH4]/d[CH2O]
  J[552] -= dqdci; // dwdot[CH3]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] += dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = -k_r * sc[11];
  J[578] += dqdci; // dwdot[CH4]/d[HCO]
  J[579] -= dqdci; // dwdot[CH3]/d[HCO]
  J[587] -= dqdci; // dwdot[CH2O]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 96: CH2O + CH3O <=> CH3OH + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[18] * sc[20];
  k_f = 662000 * exp(-(1154.37903135603) * invT);
  dlnkfdT = (1154.37903135603) * invT2;
  // reverse
  phi_r = sc[17] * sc[21];
  Kc = exp(-g_RT[17] + g_RT[18] + g_RT[20] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[18] + h_RT[20]) + (h_RT[17] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[17] += q; // CH3OH
  wdot[18] -= q; // CH3O
  wdot[20] -= q; // CH2O
  wdot[21] += q; // HCO
  // d()/d[CH3OH]
  dqdci = -k_r * sc[21];
  J[476] += dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] -= dqdci; // dwdot[CH3O]/d[CH3OH]
  J[479] -= dqdci; // dwdot[CH2O]/d[CH3OH]
  J[480] += dqdci; // dwdot[HCO]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = +k_f * sc[20];
  J[503] += dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] -= dqdci; // dwdot[CH2O]/d[CH3O]
  J[507] += dqdci; // dwdot[HCO]/d[CH3O]
  // d()/d[CH2O]
  dqdci = +k_f * sc[18];
  J[557] += dqdci; // dwdot[CH3OH]/d[CH2O]
  J[558] -= dqdci; // dwdot[CH3O]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] += dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = -k_r * sc[17];
  J[584] += dqdci; // dwdot[CH3OH]/d[HCO]
  J[585] -= dqdci; // dwdot[CH3O]/d[HCO]
  J[587] -= dqdci; // dwdot[CH2O]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[719] += dqdT; // dwdot[CH3OH]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 98: HCO + O2 <=> CO + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[21];
  k_f = 7580000 * exp(-(206.318832979935) * invT);
  dlnkfdT = (206.318832979935) * invT2;
  // reverse
  phi_r = sc[8] * sc[9];
  Kc = exp(g_RT[3] - g_RT[8] - g_RT[9] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[21]) + (h_RT[8] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[8] += q;  // HO2
  wdot[9] += q;  // CO
  wdot[21] -= q; // HCO
  // d()/d[O2]
  dqdci = +k_f * sc[21];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[89] += dqdci;  // dwdot[HO2]/d[O2]
  J[90] += dqdci;  // dwdot[CO]/d[O2]
  J[102] -= dqdci; // dwdot[HCO]/d[O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[9];
  J[219] -= dqdci; // dwdot[O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  J[225] += dqdci; // dwdot[CO]/d[HO2]
  J[237] -= dqdci; // dwdot[HCO]/d[HO2]
  // d()/d[CO]
  dqdci = -k_r * sc[8];
  J[246] -= dqdci; // dwdot[O2]/d[CO]
  J[251] += dqdci; // dwdot[HO2]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[3];
  J[570] -= dqdci; // dwdot[O2]/d[HCO]
  J[575] += dqdci; // dwdot[HO2]/d[HCO]
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 99: HCO + O <=> CO + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[21];
  k_f = 30200000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[6] * sc[9];
  Kc = exp(g_RT[4] - g_RT[6] - g_RT[9] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[21]) + (h_RT[6] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[6] += q;  // OH
  wdot[9] += q;  // CO
  wdot[21] -= q; // HCO
  // d()/d[O]
  dqdci = +k_f * sc[21];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[117] += dqdci; // dwdot[CO]/d[O]
  J[129] -= dqdci; // dwdot[HCO]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[9];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[171] += dqdci; // dwdot[CO]/d[OH]
  J[183] -= dqdci; // dwdot[HCO]/d[OH]
  // d()/d[CO]
  dqdci = -k_r * sc[6];
  J[247] -= dqdci; // dwdot[O]/d[CO]
  J[249] += dqdci; // dwdot[OH]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[4];
  J[571] -= dqdci; // dwdot[O]/d[HCO]
  J[573] += dqdci; // dwdot[OH]/d[HCO]
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 100: H + HCO <=> CO + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[21];
  k_f = 73400000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[9];
  Kc = exp(-g_RT[1] + g_RT[2] - g_RT[9] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[21]) + (h_RT[1] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[9] += q;  // CO
  wdot[21] -= q; // HCO
  // d()/d[H2]
  dqdci = -k_r * sc[9];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[36] += dqdci; // dwdot[CO]/d[H2]
  J[48] -= dqdci; // dwdot[HCO]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[21];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[63] += dqdci; // dwdot[CO]/d[H]
  J[75] -= dqdci; // dwdot[HCO]/d[H]
  // d()/d[CO]
  dqdci = -k_r * sc[1];
  J[244] += dqdci; // dwdot[H2]/d[CO]
  J[245] -= dqdci; // dwdot[H]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[2];
  J[568] += dqdci; // dwdot[H2]/d[HCO]
  J[569] -= dqdci; // dwdot[H]/d[HCO]
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 101: HCO + OH <=> CO + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[21];
  k_f = 30110000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[5] * sc[9];
  Kc = exp(-g_RT[5] + g_RT[6] - g_RT[9] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[21]) + (h_RT[5] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[9] += q;  // CO
  wdot[21] -= q; // HCO
  // d()/d[H2O]
  dqdci = -k_r * sc[9];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[144] += dqdci; // dwdot[CO]/d[H2O]
  J[156] -= dqdci; // dwdot[HCO]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[21];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[171] += dqdci; // dwdot[CO]/d[OH]
  J[183] -= dqdci; // dwdot[HCO]/d[OH]
  // d()/d[CO]
  dqdci = -k_r * sc[5];
  J[248] += dqdci; // dwdot[H2O]/d[CO]
  J[249] -= dqdci; // dwdot[OH]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[6];
  J[572] += dqdci; // dwdot[H2O]/d[HCO]
  J[573] -= dqdci; // dwdot[OH]/d[HCO]
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 102: CH3 + HCO <=> CH4 + CO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[12] * sc[21];
  k_f = 26500000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[9] * sc[11];
  Kc = exp(-g_RT[9] - g_RT[11] + g_RT[12] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[12] + h_RT[21]) + (h_RT[9] + h_RT[11]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] += q;  // CO
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  wdot[21] -= q; // HCO
  // d()/d[CO]
  dqdci = -k_r * sc[11];
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[254] += dqdci; // dwdot[CH4]/d[CO]
  J[255] -= dqdci; // dwdot[CH3]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[CH4]
  dqdci = -k_r * sc[9];
  J[306] += dqdci; // dwdot[CO]/d[CH4]
  J[308] += dqdci; // dwdot[CH4]/d[CH4]
  J[309] -= dqdci; // dwdot[CH3]/d[CH4]
  J[318] -= dqdci; // dwdot[HCO]/d[CH4]
  // d()/d[CH3]
  dqdci = +k_f * sc[21];
  J[333] += dqdci; // dwdot[CO]/d[CH3]
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[345] -= dqdci; // dwdot[HCO]/d[CH3]
  // d()/d[HCO]
  dqdci = +k_f * sc[12];
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[578] += dqdci; // dwdot[CH4]/d[HCO]
  J[579] -= dqdci; // dwdot[CH3]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 103: 2 HCO <=> CH2O + CO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[21] * sc[21]);
  k_f = 18000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[9] * sc[20];
  Kc = exp(-g_RT[9] - g_RT[20] + 2.000000 * g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[21]) + (h_RT[9] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] += q;      // CO
  wdot[20] += q;     // CH2O
  wdot[21] -= 2 * q; // HCO
  // d()/d[CO]
  dqdci = -k_r * sc[20];
  J[252] += dqdci;      // dwdot[CO]/d[CO]
  J[263] += dqdci;      // dwdot[CH2O]/d[CO]
  J[264] += -2 * dqdci; // dwdot[HCO]/d[CO]
  // d()/d[CH2O]
  dqdci = -k_r * sc[9];
  J[549] += dqdci;      // dwdot[CO]/d[CH2O]
  J[560] += dqdci;      // dwdot[CH2O]/d[CH2O]
  J[561] += -2 * dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = +k_f * 2.000000 * sc[21];
  J[576] += dqdci;      // dwdot[CO]/d[HCO]
  J[587] += dqdci;      // dwdot[CH2O]/d[HCO]
  J[588] += -2 * dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[711] += dqdT;      // dwdot[CO]/dT
  J[722] += dqdT;      // dwdot[CH2O]/dT
  J[723] += -2 * dqdT; // dwdot[HCO]/dT

  // reaction 104: HCO + O <=> CO2 + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[21];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[10];
  Kc = exp(-g_RT[2] + g_RT[4] - g_RT[10] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[21]) + (h_RT[2] + h_RT[10]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[4] -= q;  // O
  wdot[10] += q; // CO2
  wdot[21] -= q; // HCO
  // d()/d[H]
  dqdci = -k_r * sc[10];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[58] -= dqdci; // dwdot[O]/d[H]
  J[64] += dqdci; // dwdot[CO2]/d[H]
  J[75] -= dqdci; // dwdot[HCO]/d[H]
  // d()/d[O]
  dqdci = +k_f * sc[21];
  J[110] += dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[118] += dqdci; // dwdot[CO2]/d[O]
  J[129] -= dqdci; // dwdot[HCO]/d[O]
  // d()/d[CO2]
  dqdci = -k_r * sc[2];
  J[272] += dqdci; // dwdot[H]/d[CO2]
  J[274] -= dqdci; // dwdot[O]/d[CO2]
  J[280] += dqdci; // dwdot[CO2]/d[CO2]
  J[291] -= dqdci; // dwdot[HCO]/d[CO2]
  // d()/d[HCO]
  dqdci = +k_f * sc[4];
  J[569] += dqdci; // dwdot[H]/d[HCO]
  J[571] -= dqdci; // dwdot[O]/d[HCO]
  J[577] += dqdci; // dwdot[CO2]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[712] += dqdT; // dwdot[CO2]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 105: HCO + HO2 => CO2 + H + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[21];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[2] += q;  // H
  wdot[6] += q;  // OH
  wdot[8] -= q;  // HO2
  wdot[10] += q; // CO2
  wdot[21] -= q; // HCO
  // d()/d[HO2]
  dqdci = +k_f * sc[21];
  J[218] += dqdci; // dwdot[H]/d[HO2]
  J[222] += dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[226] += dqdci; // dwdot[CO2]/d[HO2]
  J[237] -= dqdci; // dwdot[HCO]/d[HO2]
  // d()/d[HCO]
  dqdci = +k_f * sc[8];
  J[569] += dqdci; // dwdot[H]/d[HCO]
  J[573] += dqdci; // dwdot[OH]/d[HCO]
  J[575] -= dqdci; // dwdot[HO2]/d[HCO]
  J[577] += dqdci; // dwdot[CO2]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[712] += dqdT; // dwdot[CO2]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 109: CH + OH <=> C + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[16];
  k_f = 40 * exp(2 * tc[0] - (1509.64999741416) * invT);
  dlnkfdT = 2 * invT + (1509.64999741416) * invT2;
  // reverse
  phi_r = sc[5] * sc[15];
  Kc = exp(-g_RT[5] + g_RT[6] - g_RT[15] + g_RT[16]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[16]) + (h_RT[5] + h_RT[15]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[15] += q; // C
  wdot[16] -= q; // CH
  // d()/d[H2O]
  dqdci = -k_r * sc[15];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[150] += dqdci; // dwdot[C]/d[H2O]
  J[151] -= dqdci; // dwdot[CH]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[16];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[177] += dqdci; // dwdot[C]/d[OH]
  J[178] -= dqdci; // dwdot[CH]/d[OH]
  // d()/d[C]
  dqdci = -k_r * sc[5];
  J[410] += dqdci; // dwdot[H2O]/d[C]
  J[411] -= dqdci; // dwdot[OH]/d[C]
  J[420] += dqdci; // dwdot[C]/d[C]
  J[421] -= dqdci; // dwdot[CH]/d[C]
  // d()/d[CH]
  dqdci = +k_f * sc[6];
  J[437] += dqdci; // dwdot[H2O]/d[CH]
  J[438] -= dqdci; // dwdot[OH]/d[CH]
  J[447] += dqdci; // dwdot[C]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[717] += dqdT; // dwdot[C]/dT
  J[718] -= dqdT; // dwdot[CH]/dT

  // reaction 110: CH2 + HCO <=> CH3 + CO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[13] * sc[21];
  k_f = 20000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[9] * sc[12];
  Kc = exp(-g_RT[9] - g_RT[12] + g_RT[13] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[13] + h_RT[21]) + (h_RT[9] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] += q;  // CO
  wdot[12] += q; // CH3
  wdot[13] -= q; // CH2
  wdot[21] -= q; // HCO
  // d()/d[CO]
  dqdci = -k_r * sc[12];
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[255] += dqdci; // dwdot[CH3]/d[CO]
  J[256] -= dqdci; // dwdot[CH2]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[CH3]
  dqdci = -k_r * sc[9];
  J[333] += dqdci; // dwdot[CO]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  J[337] -= dqdci; // dwdot[CH2]/d[CH3]
  J[345] -= dqdci; // dwdot[HCO]/d[CH3]
  // d()/d[CH2]
  dqdci = +k_f * sc[21];
  J[360] += dqdci; // dwdot[CO]/d[CH2]
  J[363] += dqdci; // dwdot[CH3]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  J[372] -= dqdci; // dwdot[HCO]/d[CH2]
  // d()/d[HCO]
  dqdci = +k_f * sc[13];
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[579] += dqdci; // dwdot[CH3]/d[HCO]
  J[580] -= dqdci; // dwdot[CH2]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[714] += dqdT; // dwdot[CH3]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 111: CH2 + O <=> CO + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[13];
  k_f = 50000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[9];
  Kc = exp(-g_RT[1] + g_RT[4] - g_RT[9] + g_RT[13]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[1] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[4] -= q;  // O
  wdot[9] += q;  // CO
  wdot[13] -= q; // CH2
  // d()/d[H2]
  dqdci = -k_r * sc[9];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[31] -= dqdci; // dwdot[O]/d[H2]
  J[36] += dqdci; // dwdot[CO]/d[H2]
  J[40] -= dqdci; // dwdot[CH2]/d[H2]
  // d()/d[O]
  dqdci = +k_f * sc[13];
  J[109] += dqdci; // dwdot[H2]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[117] += dqdci; // dwdot[CO]/d[O]
  J[121] -= dqdci; // dwdot[CH2]/d[O]
  // d()/d[CO]
  dqdci = -k_r * sc[1];
  J[244] += dqdci; // dwdot[H2]/d[CO]
  J[247] -= dqdci; // dwdot[O]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[256] -= dqdci; // dwdot[CH2]/d[CO]
  // d()/d[CH2]
  dqdci = +k_f * sc[4];
  J[352] += dqdci; // dwdot[H2]/d[CH2]
  J[355] -= dqdci; // dwdot[O]/d[CH2]
  J[360] += dqdci; // dwdot[CO]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT

  // reaction 112: CH2 + OH <=> CH2O + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[13];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[20];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[13] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[13]) + (h_RT[2] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[13] -= q; // CH2
  wdot[20] += q; // CH2O
  // d()/d[H]
  dqdci = -k_r * sc[20];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[67] -= dqdci; // dwdot[CH2]/d[H]
  J[74] += dqdci; // dwdot[CH2O]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[13];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[175] -= dqdci; // dwdot[CH2]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH2]
  dqdci = +k_f * sc[6];
  J[353] += dqdci; // dwdot[H]/d[CH2]
  J[357] -= dqdci; // dwdot[OH]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  J[371] += dqdci; // dwdot[CH2O]/d[CH2]
  // d()/d[CH2O]
  dqdci = -k_r * sc[2];
  J[542] += dqdci; // dwdot[H]/d[CH2O]
  J[546] -= dqdci; // dwdot[OH]/d[CH2O]
  J[553] -= dqdci; // dwdot[CH2]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 113: CH2 + CO2 <=> CH2O + CO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[10] * sc[13];
  k_f = 110000 * exp(-(503.21666580472) * invT);
  dlnkfdT = (503.21666580472) * invT2;
  // reverse
  phi_r = sc[9] * sc[20];
  Kc = exp(-g_RT[9] + g_RT[10] + g_RT[13] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[10] + h_RT[13]) + (h_RT[9] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] += q;  // CO
  wdot[10] -= q; // CO2
  wdot[13] -= q; // CH2
  wdot[20] += q; // CH2O
  // d()/d[CO]
  dqdci = -k_r * sc[20];
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[253] -= dqdci; // dwdot[CO2]/d[CO]
  J[256] -= dqdci; // dwdot[CH2]/d[CO]
  J[263] += dqdci; // dwdot[CH2O]/d[CO]
  // d()/d[CO2]
  dqdci = +k_f * sc[13];
  J[279] += dqdci; // dwdot[CO]/d[CO2]
  J[280] -= dqdci; // dwdot[CO2]/d[CO2]
  J[283] -= dqdci; // dwdot[CH2]/d[CO2]
  J[290] += dqdci; // dwdot[CH2O]/d[CO2]
  // d()/d[CH2]
  dqdci = +k_f * sc[10];
  J[360] += dqdci; // dwdot[CO]/d[CH2]
  J[361] -= dqdci; // dwdot[CO2]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  J[371] += dqdci; // dwdot[CH2O]/d[CH2]
  // d()/d[CH2O]
  dqdci = -k_r * sc[9];
  J[549] += dqdci; // dwdot[CO]/d[CH2O]
  J[550] -= dqdci; // dwdot[CO2]/d[CH2O]
  J[553] -= dqdci; // dwdot[CH2]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[712] -= dqdT; // dwdot[CO2]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 114: CH2(S) + O <=> CO + 2 H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[14];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = (sc[2] * sc[2]) * sc[9];
  Kc = refC * exp(-2.000000 * g_RT[2] + g_RT[4] - g_RT[9] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT =
    invT * (-(h_RT[4] + h_RT[14]) + (2.000000 * h_RT[2] + h_RT[9]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += 2 * q; // H
  wdot[4] -= q;     // O
  wdot[9] += q;     // CO
  wdot[14] -= q;    // CH2(S)
  // d()/d[H]
  dqdci = -k_r * 2.000000 * sc[2] * sc[9];
  J[56] += 2 * dqdci; // dwdot[H]/d[H]
  J[58] -= dqdci;     // dwdot[O]/d[H]
  J[63] += dqdci;     // dwdot[CO]/d[H]
  J[68] -= dqdci;     // dwdot[CH2(S)]/d[H]
  // d()/d[O]
  dqdci = +k_f * sc[14];
  J[110] += 2 * dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci;     // dwdot[O]/d[O]
  J[117] += dqdci;     // dwdot[CO]/d[O]
  J[122] -= dqdci;     // dwdot[CH2(S)]/d[O]
  // d()/d[CO]
  dqdci = -k_r * sc[2] * sc[2];
  J[245] += 2 * dqdci; // dwdot[H]/d[CO]
  J[247] -= dqdci;     // dwdot[O]/d[CO]
  J[252] += dqdci;     // dwdot[CO]/d[CO]
  J[257] -= dqdci;     // dwdot[CH2(S)]/d[CO]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[4];
  J[380] += 2 * dqdci; // dwdot[H]/d[CH2(S)]
  J[382] -= dqdci;     // dwdot[O]/d[CH2(S)]
  J[387] += dqdci;     // dwdot[CO]/d[CH2(S)]
  J[392] -= dqdci;     // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[704] += 2 * dqdT; // dwdot[H]/dT
  J[706] -= dqdT;     // dwdot[O]/dT
  J[711] += dqdT;     // dwdot[CO]/dT
  J[716] -= dqdT;     // dwdot[CH2(S)]/dT

  // reaction 115: CH2(S) + CH4 <=> 2 CH3
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[11] * sc[14];
  k_f = 43000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = (sc[12] * sc[12]);
  Kc = exp(g_RT[11] - 2.000000 * g_RT[12] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[11] + h_RT[14]) + (2.000000 * h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[11] -= q;     // CH4
  wdot[12] += 2 * q; // CH3
  wdot[14] -= q;     // CH2(S)
  // d()/d[CH4]
  dqdci = +k_f * sc[14];
  J[308] -= dqdci;     // dwdot[CH4]/d[CH4]
  J[309] += 2 * dqdci; // dwdot[CH3]/d[CH4]
  J[311] -= dqdci;     // dwdot[CH2(S)]/d[CH4]
  // d()/d[CH3]
  dqdci = -k_r * 2.000000 * sc[12];
  J[335] -= dqdci;     // dwdot[CH4]/d[CH3]
  J[336] += 2 * dqdci; // dwdot[CH3]/d[CH3]
  J[338] -= dqdci;     // dwdot[CH2(S)]/d[CH3]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[11];
  J[389] -= dqdci;     // dwdot[CH4]/d[CH2(S)]
  J[390] += 2 * dqdci; // dwdot[CH3]/d[CH2(S)]
  J[392] -= dqdci;     // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[713] -= dqdT;     // dwdot[CH4]/dT
  J[714] += 2 * dqdT; // dwdot[CH3]/dT
  J[716] -= dqdT;     // dwdot[CH2(S)]/dT

  // reaction 116: CH3OH + H => CH3 + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[17];
  k_f = 650000 * exp(-(2667.04832876501) * invT);
  dlnkfdT = (2667.04832876501) * invT2;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[2] -= q;  // H
  wdot[5] += q;  // H2O
  wdot[12] += q; // CH3
  wdot[17] -= q; // CH3OH
  // d()/d[H]
  dqdci = +k_f * sc[17];
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[59] += dqdci; // dwdot[H2O]/d[H]
  J[66] += dqdci; // dwdot[CH3]/d[H]
  J[71] -= dqdci; // dwdot[CH3OH]/d[H]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[2];
  J[461] -= dqdci; // dwdot[H]/d[CH3OH]
  J[464] += dqdci; // dwdot[H2O]/d[CH3OH]
  J[471] += dqdci; // dwdot[CH3]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  // d()/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[714] += dqdT; // dwdot[CH3]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT

  // reaction 117: CH2O + O => CO2 + 2 H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[20];
  k_f = 20000000 * exp(-(2516.0833290236) * invT);
  dlnkfdT = (2516.0833290236) * invT2;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[2] += 2 * q; // H
  wdot[4] -= q;     // O
  wdot[10] += q;    // CO2
  wdot[20] -= q;    // CH2O
  // d()/d[O]
  dqdci = +k_f * sc[20];
  J[110] += 2 * dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci;     // dwdot[O]/d[O]
  J[118] += dqdci;     // dwdot[CO2]/d[O]
  J[128] -= dqdci;     // dwdot[CH2O]/d[O]
  // d()/d[CH2O]
  dqdci = +k_f * sc[4];
  J[542] += 2 * dqdci; // dwdot[H]/d[CH2O]
  J[544] -= dqdci;     // dwdot[O]/d[CH2O]
  J[550] += dqdci;     // dwdot[CO2]/d[CH2O]
  J[560] -= dqdci;     // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[704] += 2 * dqdT; // dwdot[H]/dT
  J[706] -= dqdT;     // dwdot[O]/dT
  J[712] += dqdT;     // dwdot[CO2]/dT
  J[722] -= dqdT;     // dwdot[CH2O]/dT

  // reaction 118: CH2O + OH => CO2 + H + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[20];
  k_f = 100000;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[10] += q; // CO2
  wdot[20] -= q; // CH2O
  // d()/d[OH]
  dqdci = +k_f * sc[20];
  J[163] += dqdci; // dwdot[H2]/d[OH]
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[172] += dqdci; // dwdot[CO2]/d[OH]
  J[182] -= dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH2O]
  dqdci = +k_f * sc[6];
  J[541] += dqdci; // dwdot[H2]/d[CH2O]
  J[542] += dqdci; // dwdot[H]/d[CH2O]
  J[546] -= dqdci; // dwdot[OH]/d[CH2O]
  J[550] += dqdci; // dwdot[CO2]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[712] += dqdT; // dwdot[CO2]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT

  // reaction 119: CH3O + CO <=> CH3 + CO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[9] * sc[18];
  k_f = 500000 * exp(-(3270.90832773068) * invT);
  dlnkfdT = (3270.90832773068) * invT2;
  // reverse
  phi_r = sc[10] * sc[12];
  Kc = exp(g_RT[9] - g_RT[10] - g_RT[12] + g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[9] + h_RT[18]) + (h_RT[10] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] -= q;  // CO
  wdot[10] += q; // CO2
  wdot[12] += q; // CH3
  wdot[18] -= q; // CH3O
  // d()/d[CO]
  dqdci = +k_f * sc[18];
  J[252] -= dqdci; // dwdot[CO]/d[CO]
  J[253] += dqdci; // dwdot[CO2]/d[CO]
  J[255] += dqdci; // dwdot[CH3]/d[CO]
  J[261] -= dqdci; // dwdot[CH3O]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[12];
  J[279] -= dqdci; // dwdot[CO]/d[CO2]
  J[280] += dqdci; // dwdot[CO2]/d[CO2]
  J[282] += dqdci; // dwdot[CH3]/d[CO2]
  J[288] -= dqdci; // dwdot[CH3O]/d[CO2]
  // d()/d[CH3]
  dqdci = -k_r * sc[10];
  J[333] -= dqdci; // dwdot[CO]/d[CH3]
  J[334] += dqdci; // dwdot[CO2]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  J[342] -= dqdci; // dwdot[CH3O]/d[CH3]
  // d()/d[CH3O]
  dqdci = +k_f * sc[9];
  J[495] -= dqdci; // dwdot[CO]/d[CH3O]
  J[496] += dqdci; // dwdot[CO2]/d[CH3O]
  J[498] += dqdci; // dwdot[CH3]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[711] -= dqdT; // dwdot[CO]/dT
  J[712] += dqdT; // dwdot[CO2]/dT
  J[714] += dqdT; // dwdot[CH3]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT

  // reaction 121: HCO + HO2 <=> CO + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[21];
  k_f = 400000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[7] * sc[9];
  Kc = exp(-g_RT[7] + g_RT[8] - g_RT[9] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[21]) + (h_RT[7] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[9] += q;  // CO
  wdot[21] -= q; // HCO
  // d()/d[H2O2]
  dqdci = -k_r * sc[9];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[198] += dqdci; // dwdot[CO]/d[H2O2]
  J[210] -= dqdci; // dwdot[HCO]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[21];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[225] += dqdci; // dwdot[CO]/d[HO2]
  J[237] -= dqdci; // dwdot[HCO]/d[HO2]
  // d()/d[CO]
  dqdci = -k_r * sc[7];
  J[250] += dqdci; // dwdot[H2O2]/d[CO]
  J[251] -= dqdci; // dwdot[HO2]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[8];
  J[574] += dqdci; // dwdot[H2O2]/d[HCO]
  J[575] -= dqdci; // dwdot[HO2]/d[HCO]
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 122: CH3O + OH => CH2O + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[18];
  k_f = 15000000;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[18] -= q; // CH3O
  wdot[20] += q; // CH2O
  // d()/d[OH]
  dqdci = +k_f * sc[18];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[180] -= dqdci; // dwdot[CH3O]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH3O]
  dqdci = +k_f * sc[6];
  J[491] += dqdci; // dwdot[H2O]/d[CH3O]
  J[492] -= dqdci; // dwdot[OH]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] += dqdci; // dwdot[CH2O]/d[CH3O]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 123: CH2OH + CH3 => CH2O + CH4
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[12] * sc[19];
  k_f = 15000000;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[CH3]
  dqdci = +k_f * sc[19];
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[343] -= dqdci; // dwdot[CH2OH]/d[CH3]
  J[344] += dqdci; // dwdot[CH2O]/d[CH3]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[12];
  J[524] += dqdci; // dwdot[CH4]/d[CH2OH]
  J[525] -= dqdci; // dwdot[CH3]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 124: CH + O <=> CHO+ + E
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[16];
  k_f = 1746000000000 * exp(-2.19 * tc[0] - (164.999209334044) * invT);
  dlnkfdT = -2.19 * invT + (164.999209334044) * invT2;
  // reverse
  phi_r = sc[22] * sc[25];
  Kc = exp(g_RT[4] + g_RT[16] - g_RT[22] - g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[16]) + (h_RT[22] + h_RT[25]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[16] -= q; // CH
  wdot[22] += q; // CHO+
  wdot[25] += q; // E
  // d()/d[O]
  dqdci = +k_f * sc[16];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[124] -= dqdci; // dwdot[CH]/d[O]
  J[130] += dqdci; // dwdot[CHO+]/d[O]
  J[133] += dqdci; // dwdot[E]/d[O]
  // d()/d[CH]
  dqdci = +k_f * sc[4];
  J[436] -= dqdci; // dwdot[O]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  J[454] += dqdci; // dwdot[CHO+]/d[CH]
  J[457] += dqdci; // dwdot[E]/d[CH]
  // d()/d[CHO+]
  dqdci = -k_r * sc[25];
  J[598] -= dqdci; // dwdot[O]/d[CHO+]
  J[610] -= dqdci; // dwdot[CH]/d[CHO+]
  J[616] += dqdci; // dwdot[CHO+]/d[CHO+]
  J[619] += dqdci; // dwdot[E]/d[CHO+]
  // d()/d[E]
  dqdci = -k_r * sc[22];
  J[679] -= dqdci; // dwdot[O]/d[E]
  J[691] -= dqdci; // dwdot[CH]/d[E]
  J[697] += dqdci; // dwdot[CHO+]/d[E]
  J[700] += dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[718] -= dqdT; // dwdot[CH]/dT
  J[724] += dqdT; // dwdot[CHO+]/dT
  J[727] += dqdT; // dwdot[E]/dT

  // reaction 125: CHO+ + H2O <=> CO + H3O+
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[22];
  k_f = 26080000000 * exp(-0.5 * tc[0]);
  dlnkfdT = -0.5 * invT;
  // reverse
  phi_r = sc[9] * sc[23];
  Kc = exp(g_RT[5] - g_RT[9] + g_RT[22] - g_RT[23]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[22]) + (h_RT[9] + h_RT[23]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] -= q;  // H2O
  wdot[9] += q;  // CO
  wdot[22] -= q; // CHO+
  wdot[23] += q; // H3O+
  // d()/d[H2O]
  dqdci = +k_f * sc[22];
  J[140] -= dqdci; // dwdot[H2O]/d[H2O]
  J[144] += dqdci; // dwdot[CO]/d[H2O]
  J[157] -= dqdci; // dwdot[CHO+]/d[H2O]
  J[158] += dqdci; // dwdot[H3O+]/d[H2O]
  // d()/d[CO]
  dqdci = -k_r * sc[23];
  J[248] -= dqdci; // dwdot[H2O]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[265] -= dqdci; // dwdot[CHO+]/d[CO]
  J[266] += dqdci; // dwdot[H3O+]/d[CO]
  // d()/d[CHO+]
  dqdci = +k_f * sc[5];
  J[599] -= dqdci; // dwdot[H2O]/d[CHO+]
  J[603] += dqdci; // dwdot[CO]/d[CHO+]
  J[616] -= dqdci; // dwdot[CHO+]/d[CHO+]
  J[617] += dqdci; // dwdot[H3O+]/d[CHO+]
  // d()/d[H3O+]
  dqdci = -k_r * sc[9];
  J[626] -= dqdci; // dwdot[H2O]/d[H3O+]
  J[630] += dqdci; // dwdot[CO]/d[H3O+]
  J[643] -= dqdci; // dwdot[CHO+]/d[H3O+]
  J[644] += dqdci; // dwdot[H3O+]/d[H3O+]
  // d()/dT
  J[707] -= dqdT; // dwdot[H2O]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[724] -= dqdT; // dwdot[CHO+]/dT
  J[725] += dqdT; // dwdot[H3O+]/dT

  // reaction 126: E + H3O+ <=> H + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[23] * sc[25];
  k_f = 739500000000 * exp(-0.5 * tc[0]);
  dlnkfdT = -0.5 * invT;
  // reverse
  phi_r = sc[2] * sc[5];
  Kc = exp(-g_RT[2] - g_RT[5] + g_RT[23] + g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[23] + h_RT[25]) + (h_RT[2] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[5] += q;  // H2O
  wdot[23] -= q; // H3O+
  wdot[25] -= q; // E
  // d()/d[H]
  dqdci = -k_r * sc[5];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[59] += dqdci; // dwdot[H2O]/d[H]
  J[77] -= dqdci; // dwdot[H3O+]/d[H]
  J[79] -= dqdci; // dwdot[E]/d[H]
  // d()/d[H2O]
  dqdci = -k_r * sc[2];
  J[137] += dqdci; // dwdot[H]/d[H2O]
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[158] -= dqdci; // dwdot[H3O+]/d[H2O]
  J[160] -= dqdci; // dwdot[E]/d[H2O]
  // d()/d[H3O+]
  dqdci = +k_f * sc[25];
  J[623] += dqdci; // dwdot[H]/d[H3O+]
  J[626] += dqdci; // dwdot[H2O]/d[H3O+]
  J[644] -= dqdci; // dwdot[H3O+]/d[H3O+]
  J[646] -= dqdci; // dwdot[E]/d[H3O+]
  // d()/d[E]
  dqdci = +k_f * sc[23];
  J[677] += dqdci; // dwdot[H]/d[E]
  J[680] += dqdci; // dwdot[H2O]/d[E]
  J[698] -= dqdci; // dwdot[H3O+]/d[E]
  J[700] -= dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[725] -= dqdT; // dwdot[H3O+]/dT
  J[727] -= dqdT; // dwdot[E]/dT

  // reaction 127: E + H3O+ <=> 2 H + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[23] * sc[25];
  k_f = 3181000000000 * exp(-0.5 * tc[0]);
  dlnkfdT = -0.5 * invT;
  // reverse
  phi_r = (sc[2] * sc[2]) * sc[6];
  Kc = refC * exp(-2.000000 * g_RT[2] - g_RT[6] + g_RT[23] + g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT =
    invT * (-(h_RT[23] + h_RT[25]) + (2.000000 * h_RT[2] + h_RT[6]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += 2 * q; // H
  wdot[6] += q;     // OH
  wdot[23] -= q;    // H3O+
  wdot[25] -= q;    // E
  // d()/d[H]
  dqdci = -k_r * 2.000000 * sc[2] * sc[6];
  J[56] += 2 * dqdci; // dwdot[H]/d[H]
  J[60] += dqdci;     // dwdot[OH]/d[H]
  J[77] -= dqdci;     // dwdot[H3O+]/d[H]
  J[79] -= dqdci;     // dwdot[E]/d[H]
  // d()/d[OH]
  dqdci = -k_r * sc[2] * sc[2];
  J[164] += 2 * dqdci; // dwdot[H]/d[OH]
  J[168] += dqdci;     // dwdot[OH]/d[OH]
  J[185] -= dqdci;     // dwdot[H3O+]/d[OH]
  J[187] -= dqdci;     // dwdot[E]/d[OH]
  // d()/d[H3O+]
  dqdci = +k_f * sc[25];
  J[623] += 2 * dqdci; // dwdot[H]/d[H3O+]
  J[627] += dqdci;     // dwdot[OH]/d[H3O+]
  J[644] -= dqdci;     // dwdot[H3O+]/d[H3O+]
  J[646] -= dqdci;     // dwdot[E]/d[H3O+]
  // d()/d[E]
  dqdci = +k_f * sc[23];
  J[677] += 2 * dqdci; // dwdot[H]/d[E]
  J[681] += dqdci;     // dwdot[OH]/d[E]
  J[698] -= dqdci;     // dwdot[H3O+]/d[E]
  J[700] -= dqdci;     // dwdot[E]/d[E]
  // d()/dT
  J[704] += 2 * dqdT; // dwdot[H]/dT
  J[708] += dqdT;     // dwdot[OH]/dT
  J[725] -= dqdT;     // dwdot[H3O+]/dT
  J[727] -= dqdT;     // dwdot[E]/dT

  // reaction 128: E + H3O+ <=> H2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[23] * sc[25];
  k_f = 560100000000 * exp(-0.5 * tc[0]);
  dlnkfdT = -0.5 * invT;
  // reverse
  phi_r = sc[1] * sc[6];
  Kc = exp(-g_RT[1] - g_RT[6] + g_RT[23] + g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[23] + h_RT[25]) + (h_RT[1] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[6] += q;  // OH
  wdot[23] -= q; // H3O+
  wdot[25] -= q; // E
  // d()/d[H2]
  dqdci = -k_r * sc[6];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[33] += dqdci; // dwdot[OH]/d[H2]
  J[50] -= dqdci; // dwdot[H3O+]/d[H2]
  J[52] -= dqdci; // dwdot[E]/d[H2]
  // d()/d[OH]
  dqdci = -k_r * sc[1];
  J[163] += dqdci; // dwdot[H2]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[185] -= dqdci; // dwdot[H3O+]/d[OH]
  J[187] -= dqdci; // dwdot[E]/d[OH]
  // d()/d[H3O+]
  dqdci = +k_f * sc[25];
  J[622] += dqdci; // dwdot[H2]/d[H3O+]
  J[627] += dqdci; // dwdot[OH]/d[H3O+]
  J[644] -= dqdci; // dwdot[H3O+]/d[H3O+]
  J[646] -= dqdci; // dwdot[E]/d[H3O+]
  // d()/d[E]
  dqdci = +k_f * sc[23];
  J[676] += dqdci; // dwdot[H2]/d[E]
  J[681] += dqdci; // dwdot[OH]/d[E]
  J[698] -= dqdci; // dwdot[H3O+]/d[E]
  J[700] -= dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[725] -= dqdT; // dwdot[H3O+]/dT
  J[727] -= dqdT; // dwdot[E]/dT

  // reaction 129: H2 + O2- <=> E + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[24];
  k_f = 602000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[7] * sc[25];
  Kc = exp(g_RT[1] - g_RT[7] + g_RT[24] - g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[24]) + (h_RT[7] + h_RT[25]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q;  // H2
  wdot[7] += q;  // H2O2
  wdot[24] -= q; // O2-
  wdot[25] += q; // E
  // d()/d[H2]
  dqdci = +k_f * sc[24];
  J[28] -= dqdci; // dwdot[H2]/d[H2]
  J[34] += dqdci; // dwdot[H2O2]/d[H2]
  J[51] -= dqdci; // dwdot[O2-]/d[H2]
  J[52] += dqdci; // dwdot[E]/d[H2]
  // d()/d[H2O2]
  dqdci = -k_r * sc[25];
  J[190] -= dqdci; // dwdot[H2]/d[H2O2]
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[213] -= dqdci; // dwdot[O2-]/d[H2O2]
  J[214] += dqdci; // dwdot[E]/d[H2O2]
  // d()/d[O2-]
  dqdci = +k_f * sc[1];
  J[649] -= dqdci; // dwdot[H2]/d[O2-]
  J[655] += dqdci; // dwdot[H2O2]/d[O2-]
  J[672] -= dqdci; // dwdot[O2-]/d[O2-]
  J[673] += dqdci; // dwdot[E]/d[O2-]
  // d()/d[E]
  dqdci = -k_r * sc[7];
  J[676] -= dqdci; // dwdot[H2]/d[E]
  J[682] += dqdci; // dwdot[H2O2]/d[E]
  J[699] -= dqdci; // dwdot[O2-]/d[E]
  J[700] += dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[703] -= dqdT; // dwdot[H2]/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[726] -= dqdT; // dwdot[O2-]/dT
  J[727] += dqdT; // dwdot[E]/dT

  // reaction 130: H + O2- <=> E + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[24];
  k_f = 723000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[8] * sc[25];
  Kc = exp(g_RT[2] - g_RT[8] + g_RT[24] - g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[24]) + (h_RT[8] + h_RT[25]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q;  // H
  wdot[8] += q;  // HO2
  wdot[24] -= q; // O2-
  wdot[25] += q; // E
  // d()/d[H]
  dqdci = +k_f * sc[24];
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[62] += dqdci; // dwdot[HO2]/d[H]
  J[78] -= dqdci; // dwdot[O2-]/d[H]
  J[79] += dqdci; // dwdot[E]/d[H]
  // d()/d[HO2]
  dqdci = -k_r * sc[25];
  J[218] -= dqdci; // dwdot[H]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  J[240] -= dqdci; // dwdot[O2-]/d[HO2]
  J[241] += dqdci; // dwdot[E]/d[HO2]
  // d()/d[O2-]
  dqdci = +k_f * sc[2];
  J[650] -= dqdci; // dwdot[H]/d[O2-]
  J[656] += dqdci; // dwdot[HO2]/d[O2-]
  J[672] -= dqdci; // dwdot[O2-]/d[O2-]
  J[673] += dqdci; // dwdot[E]/d[O2-]
  // d()/d[E]
  dqdci = -k_r * sc[8];
  J[677] -= dqdci; // dwdot[H]/d[E]
  J[683] += dqdci; // dwdot[HO2]/d[E]
  J[699] -= dqdci; // dwdot[O2-]/d[E]
  J[700] += dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[710] += dqdT; // dwdot[HO2]/dT
  J[726] -= dqdT; // dwdot[O2-]/dT
  J[727] += dqdT; // dwdot[E]/dT

  amrex::Real c_R[26], dcRdT[26], e_RT[26];
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
  for (int k = 0; k < 26; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[702 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 26; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 26; ++m) {
      dehmixdc += eh_RT[m] * J[k * 27 + m];
    }
    J[k * 27 + 26] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[728] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
#else
  amrex::Abort();
#endif
}

// compute an approx to the reaction Jacobian (for preconditioning)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT_SIMPLIFIED(
  amrex::Real* J, const amrex::Real* sc, const amrex::Real* Tp, const int* HP)
{
  amrex::Real c[26];

  for (int k = 0; k < 26; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian_precond(J, c, *Tp, *HP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 26; k++) {
    J[702 + k] *= 1.e-6;
    J[k * 27 + 26] *= 1.e6;
  }
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian(
  amrex::Real* J, const amrex::Real* sc, const amrex::Real T, const int consP)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  for (int i = 0; i < 729; i++) {
    J[i] = 0.0;
  }

  amrex::Real wdot[26];
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
  for (int k = 0; k < 26; ++k) {
    mixture += sc[k];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[26];
  gibbs(g_RT, tc);

  // compute the species enthalpy
  amrex::Real h_RT[26];
  speciesEnthalpy(h_RT, tc);

  amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  amrex::Real dqdci, dcdc_fac, dqdc[26];
  amrex::Real Pr, fPr, F, k_0, logPr;
  amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
  amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
  amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const amrex::Real ln10 = log(10.0);
  const amrex::Real log10e = 1.0 / log(10.0);
  // reaction 6: H2O2 (+M) <=> 2 OH (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (1.5 - 1) * sc[0] + (3.7 - 1) * sc[1] + (1.2 - 1) * sc[3] +
          (7.65 - 1) * sc[5] + (7.7 - 1) * sc[7] + (2.8 - 1) * sc[9] +
          (1.6 - 1) * sc[10];
  // forward
  phi_f = sc[7];
  k_f = 2000000000000 * exp(0.9 * tc[0] - (24531.3092413143) * invT);
  dlnkfdT = 0.9 * invT + (24531.3092413143) * invT2;
  // pressure-fall-off
  k_0 = 2.49e+24 * exp(-2.3 * tc[0] - (24531.3092413143) * invT);
  Pr = 1e-6 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -2.3 * invT + (24531.3092413143) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.43) * exp(-T / 1e-30);
  Fcent2 = 0.43 * exp(-T / 1e+30);
  Fcent3 = exp(-1e-16 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT = log10e / Fcent *
                (-Fcent1 / 1e-30 - Fcent2 / 1e+30 + Fcent3 * 1e-16 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = (sc[6] * sc[6]);
  Kc = refC * exp(-2.000000 * g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[7]) + (2.000000 * h_RT[6]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[6] += 2 * q; // OH
  wdot[7] -= q;     // H2O2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP == 1) {
    // d()/d[N2]
    dqdci = (1.5 - 1) * dcdc_fac;
    J[6] += 2 * dqdci; // dwdot[OH]/d[N2]
    J[7] -= dqdci;     // dwdot[H2O2]/d[N2]
    // d()/d[H2]
    dqdci = (3.7 - 1) * dcdc_fac;
    J[33] += 2 * dqdci; // dwdot[OH]/d[H2]
    J[34] -= dqdci;     // dwdot[H2O2]/d[H2]
    // d()/d[O2]
    dqdci = (1.2 - 1) * dcdc_fac;
    J[87] += 2 * dqdci; // dwdot[OH]/d[O2]
    J[88] -= dqdci;     // dwdot[H2O2]/d[O2]
    // d()/d[H2O]
    dqdci = (7.65 - 1) * dcdc_fac;
    J[141] += 2 * dqdci; // dwdot[OH]/d[H2O]
    J[142] -= dqdci;     // dwdot[H2O2]/d[H2O]
    // d()/d[OH]
    dqdci = -k_r * 2.000000 * sc[6];
    J[168] += 2 * dqdci; // dwdot[OH]/d[OH]
    J[169] -= dqdci;     // dwdot[H2O2]/d[OH]
    // d()/d[H2O2]
    dqdci = (7.7 - 1) * dcdc_fac + k_f;
    J[195] += 2 * dqdci; // dwdot[OH]/d[H2O2]
    J[196] -= dqdci;     // dwdot[H2O2]/d[H2O2]
    // d()/d[CO]
    dqdci = (2.8 - 1) * dcdc_fac;
    J[249] += 2 * dqdci; // dwdot[OH]/d[CO]
    J[250] -= dqdci;     // dwdot[H2O2]/d[CO]
    // d()/d[CO2]
    dqdci = (1.6 - 1) * dcdc_fac;
    J[276] += 2 * dqdci; // dwdot[OH]/d[CO2]
    J[277] -= dqdci;     // dwdot[H2O2]/d[CO2]
  } else {
    dqdc[0] = 1.5 * dcdc_fac;
    dqdc[1] = 3.7 * dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = 1.2 * dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 7.65 * dcdc_fac;
    dqdc[6] = dcdc_fac - k_r * 2.000000 * sc[6];
    dqdc[7] = 7.7 * dcdc_fac + k_f;
    dqdc[8] = dcdc_fac;
    dqdc[9] = 2.8 * dcdc_fac;
    dqdc[10] = 1.6 * dcdc_fac;
    dqdc[11] = dcdc_fac;
    dqdc[12] = dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 6] += 2 * dqdc[k];
      J[27 * k + 7] -= dqdc[k];
    }
  }
  J[708] += 2 * dqdT; // dwdot[OH]/dT
  J[709] -= dqdT;     // dwdot[H2O2]/dT

  // reaction 19: H + O2 (+M) <=> HO2 (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (1.3 - 1) * sc[1] + (10 - 1) * sc[5] + (1.9 - 1) * sc[9] +
          (3.8 - 1) * sc[10] + sc[11];
  // forward
  phi_f = sc[2] * sc[3];
  k_f = 4650000 * exp(0.44 * tc[0]);
  dlnkfdT = 0.44 * invT;
  // pressure-fall-off
  k_0 = 1.74e+19 * exp(-1.23 * tc[0]);
  Pr = 1e-12 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -1.23 * invT;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.67) * exp(-T / 1e-30);
  Fcent2 = 0.67 * exp(-T / 1e+30);
  Fcent3 = exp(-1e+30 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT = log10e / Fcent *
                (-Fcent1 / 1e-30 - Fcent2 / 1e+30 + Fcent3 * 1e+30 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[8];
  Kc = refCinv * exp(g_RT[2] + g_RT[3] - g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[8]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[2] -= q; // H
  wdot[3] -= q; // O2
  wdot[8] += q; // HO2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (1.3 - 1) * dcdc_fac;
    J[29] -= dqdci; // dwdot[H]/d[H2]
    J[30] -= dqdci; // dwdot[O2]/d[H2]
    J[35] += dqdci; // dwdot[HO2]/d[H2]
    // d()/d[H]
    dqdci = +k_f * sc[3];
    J[56] -= dqdci; // dwdot[H]/d[H]
    J[57] -= dqdci; // dwdot[O2]/d[H]
    J[62] += dqdci; // dwdot[HO2]/d[H]
    // d()/d[O2]
    dqdci = +k_f * sc[2];
    J[83] -= dqdci; // dwdot[H]/d[O2]
    J[84] -= dqdci; // dwdot[O2]/d[O2]
    J[89] += dqdci; // dwdot[HO2]/d[O2]
    // d()/d[H2O]
    dqdci = (10 - 1) * dcdc_fac;
    J[137] -= dqdci; // dwdot[H]/d[H2O]
    J[138] -= dqdci; // dwdot[O2]/d[H2O]
    J[143] += dqdci; // dwdot[HO2]/d[H2O]
    // d()/d[HO2]
    dqdci = -k_r;
    J[218] -= dqdci; // dwdot[H]/d[HO2]
    J[219] -= dqdci; // dwdot[O2]/d[HO2]
    J[224] += dqdci; // dwdot[HO2]/d[HO2]
    // d()/d[CO]
    dqdci = (1.9 - 1) * dcdc_fac;
    J[245] -= dqdci; // dwdot[H]/d[CO]
    J[246] -= dqdci; // dwdot[O2]/d[CO]
    J[251] += dqdci; // dwdot[HO2]/d[CO]
    // d()/d[CO2]
    dqdci = (3.8 - 1) * dcdc_fac;
    J[272] -= dqdci; // dwdot[H]/d[CO2]
    J[273] -= dqdci; // dwdot[O2]/d[CO2]
    J[278] += dqdci; // dwdot[HO2]/d[CO2]
    // d()/d[CH4]
    dqdci = (2 - 1) * dcdc_fac;
    J[299] -= dqdci; // dwdot[H]/d[CH4]
    J[300] -= dqdci; // dwdot[O2]/d[CH4]
    J[305] += dqdci; // dwdot[HO2]/d[CH4]
  } else {
    dqdc[0] = dcdc_fac;
    dqdc[1] = 1.3 * dcdc_fac;
    dqdc[2] = dcdc_fac + k_f * sc[3];
    dqdc[3] = dcdc_fac + k_f * sc[2];
    dqdc[4] = dcdc_fac;
    dqdc[5] = 10 * dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac - k_r;
    dqdc[9] = 1.9 * dcdc_fac;
    dqdc[10] = 3.8 * dcdc_fac;
    dqdc[11] = 2 * dcdc_fac;
    dqdc[12] = dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 2] -= dqdc[k];
      J[27 * k + 3] -= dqdc[k];
      J[27 * k + 8] += dqdc[k];
    }
  }
  J[704] -= dqdT; // dwdot[H]/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT

  // reaction 24: CH3 + H (+M) <=> CH4 (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha =
    mixture + sc[1] + (6 - 1) * sc[5] + (1.5 - 1) * sc[9] + sc[10] + sc[11];
  // forward
  phi_f = sc[2] * sc[12];
  k_f = 12700000000 * exp(-0.63 * tc[0] - (192.731983003208) * invT);
  dlnkfdT = -0.63 * invT + (192.731983003208) * invT2;
  // pressure-fall-off
  k_0 = 2.48e+33 * exp(-4.76 * tc[0] - (1227.84866456352) * invT);
  Pr = 1e-12 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -4.76 * invT + (1227.84866456352) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.783) * exp(-T / 74);
  Fcent2 = 0.783 * exp(-T / 2941);
  Fcent3 = exp(-6964 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT =
    log10e / Fcent * (-Fcent1 / 74 - Fcent2 / 2941 + Fcent3 * 6964 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[11];
  Kc = refCinv * exp(g_RT[2] - g_RT[11] + g_RT[12]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[12]) + (h_RT[11]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[2] -= q;  // H
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2 - 1) * dcdc_fac;
    J[29] -= dqdci; // dwdot[H]/d[H2]
    J[38] += dqdci; // dwdot[CH4]/d[H2]
    J[39] -= dqdci; // dwdot[CH3]/d[H2]
    // d()/d[H]
    dqdci = +k_f * sc[12];
    J[56] -= dqdci; // dwdot[H]/d[H]
    J[65] += dqdci; // dwdot[CH4]/d[H]
    J[66] -= dqdci; // dwdot[CH3]/d[H]
    // d()/d[H2O]
    dqdci = (6 - 1) * dcdc_fac;
    J[137] -= dqdci; // dwdot[H]/d[H2O]
    J[146] += dqdci; // dwdot[CH4]/d[H2O]
    J[147] -= dqdci; // dwdot[CH3]/d[H2O]
    // d()/d[CO]
    dqdci = (1.5 - 1) * dcdc_fac;
    J[245] -= dqdci; // dwdot[H]/d[CO]
    J[254] += dqdci; // dwdot[CH4]/d[CO]
    J[255] -= dqdci; // dwdot[CH3]/d[CO]
    // d()/d[CO2]
    dqdci = (2 - 1) * dcdc_fac;
    J[272] -= dqdci; // dwdot[H]/d[CO2]
    J[281] += dqdci; // dwdot[CH4]/d[CO2]
    J[282] -= dqdci; // dwdot[CH3]/d[CO2]
    // d()/d[CH4]
    dqdci = (2 - 1) * dcdc_fac - k_r;
    J[299] -= dqdci; // dwdot[H]/d[CH4]
    J[308] += dqdci; // dwdot[CH4]/d[CH4]
    J[309] -= dqdci; // dwdot[CH3]/d[CH4]
    // d()/d[CH3]
    dqdci = +k_f * sc[2];
    J[326] -= dqdci; // dwdot[H]/d[CH3]
    J[335] += dqdci; // dwdot[CH4]/d[CH3]
    J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  } else {
    dqdc[0] = dcdc_fac;
    dqdc[1] = 2 * dcdc_fac;
    dqdc[2] = dcdc_fac + k_f * sc[12];
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 6 * dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = 1.5 * dcdc_fac;
    dqdc[10] = 2 * dcdc_fac;
    dqdc[11] = 2 * dcdc_fac - k_r;
    dqdc[12] = dcdc_fac + k_f * sc[2];
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 2] -= dqdc[k];
      J[27 * k + 11] += dqdc[k];
      J[27 * k + 12] -= dqdc[k];
    }
  }
  J[704] -= dqdT; // dwdot[H]/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT

  // reaction 43: CH2 + H (+M) <=> CH3 (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha =
    mixture + sc[1] + (6 - 1) * sc[5] + (1.5 - 1) * sc[9] + sc[10] + sc[11];
  // forward
  phi_f = sc[2] * sc[13];
  k_f = 25000000000 * exp(-0.8 * tc[0]);
  dlnkfdT = -0.8 * invT;
  // pressure-fall-off
  k_0 = 3.2e+27 * exp(-3.14 * tc[0] - (618.956498939805) * invT);
  Pr = 1e-12 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -3.14 * invT + (618.956498939805) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.68) * exp(-T / 78);
  Fcent2 = 0.68 * exp(-T / 1995);
  Fcent3 = exp(-5590 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT =
    log10e / Fcent * (-Fcent1 / 78 - Fcent2 / 1995 + Fcent3 * 5590 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[12];
  Kc = refCinv * exp(g_RT[2] - g_RT[12] + g_RT[13]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[12]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[2] -= q;  // H
  wdot[12] += q; // CH3
  wdot[13] -= q; // CH2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2 - 1) * dcdc_fac;
    J[29] -= dqdci; // dwdot[H]/d[H2]
    J[39] += dqdci; // dwdot[CH3]/d[H2]
    J[40] -= dqdci; // dwdot[CH2]/d[H2]
    // d()/d[H]
    dqdci = +k_f * sc[13];
    J[56] -= dqdci; // dwdot[H]/d[H]
    J[66] += dqdci; // dwdot[CH3]/d[H]
    J[67] -= dqdci; // dwdot[CH2]/d[H]
    // d()/d[H2O]
    dqdci = (6 - 1) * dcdc_fac;
    J[137] -= dqdci; // dwdot[H]/d[H2O]
    J[147] += dqdci; // dwdot[CH3]/d[H2O]
    J[148] -= dqdci; // dwdot[CH2]/d[H2O]
    // d()/d[CO]
    dqdci = (1.5 - 1) * dcdc_fac;
    J[245] -= dqdci; // dwdot[H]/d[CO]
    J[255] += dqdci; // dwdot[CH3]/d[CO]
    J[256] -= dqdci; // dwdot[CH2]/d[CO]
    // d()/d[CO2]
    dqdci = (2 - 1) * dcdc_fac;
    J[272] -= dqdci; // dwdot[H]/d[CO2]
    J[282] += dqdci; // dwdot[CH3]/d[CO2]
    J[283] -= dqdci; // dwdot[CH2]/d[CO2]
    // d()/d[CH4]
    dqdci = (2 - 1) * dcdc_fac;
    J[299] -= dqdci; // dwdot[H]/d[CH4]
    J[309] += dqdci; // dwdot[CH3]/d[CH4]
    J[310] -= dqdci; // dwdot[CH2]/d[CH4]
    // d()/d[CH3]
    dqdci = -k_r;
    J[326] -= dqdci; // dwdot[H]/d[CH3]
    J[336] += dqdci; // dwdot[CH3]/d[CH3]
    J[337] -= dqdci; // dwdot[CH2]/d[CH3]
    // d()/d[CH2]
    dqdci = +k_f * sc[2];
    J[353] -= dqdci; // dwdot[H]/d[CH2]
    J[363] += dqdci; // dwdot[CH3]/d[CH2]
    J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  } else {
    dqdc[0] = dcdc_fac;
    dqdc[1] = 2 * dcdc_fac;
    dqdc[2] = dcdc_fac + k_f * sc[13];
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 6 * dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = 1.5 * dcdc_fac;
    dqdc[10] = 2 * dcdc_fac;
    dqdc[11] = 2 * dcdc_fac;
    dqdc[12] = dcdc_fac - k_r;
    dqdc[13] = dcdc_fac + k_f * sc[2];
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 2] -= dqdc[k];
      J[27 * k + 12] += dqdc[k];
      J[27 * k + 13] -= dqdc[k];
    }
  }
  J[704] -= dqdT; // dwdot[H]/dT
  J[714] += dqdT; // dwdot[CH3]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT

  // reaction 66: CH3OH (+M) <=> CH3 + OH (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture;
  // forward
  phi_f = sc[17];
  k_f = 2.084e+18 * exp(-0.615 * tc[0] - (46567.9721835682) * invT);
  dlnkfdT = -0.615 * invT + (46567.9721835682) * invT2;
  // pressure-fall-off
  k_0 = 1.5e+43 * exp(-6.995 * tc[0] - (49311.3081588693) * invT);
  Pr = 1e-6 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -6.995 * invT + (49311.3081588693) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. + 0.4748) * exp(-T / 35580);
  Fcent2 = -0.4748 * exp(-T / 1116);
  Fcent3 = exp(-9023 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT =
    log10e / Fcent * (-Fcent1 / 35580 - Fcent2 / 1116 + Fcent3 * 9023 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[6] * sc[12];
  Kc = refC * exp(-g_RT[6] - g_RT[12] + g_RT[17]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[17]) + (h_RT[6] + h_RT[12]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[6] += q;  // OH
  wdot[12] += q; // CH3
  wdot[17] -= q; // CH3OH
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP == 1) {
    // d()/d[OH]
    dqdci = -k_r * sc[12];
    J[168] += dqdci; // dwdot[OH]/d[OH]
    J[174] += dqdci; // dwdot[CH3]/d[OH]
    J[179] -= dqdci; // dwdot[CH3OH]/d[OH]
    // d()/d[CH3]
    dqdci = -k_r * sc[6];
    J[330] += dqdci; // dwdot[OH]/d[CH3]
    J[336] += dqdci; // dwdot[CH3]/d[CH3]
    J[341] -= dqdci; // dwdot[CH3OH]/d[CH3]
    // d()/d[CH3OH]
    dqdci = +k_f;
    J[465] += dqdci; // dwdot[OH]/d[CH3OH]
    J[471] += dqdci; // dwdot[CH3]/d[CH3OH]
    J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  } else {
    dqdc[0] = dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = dcdc_fac - k_r * sc[12];
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    dqdc[12] = dcdc_fac - k_r * sc[6];
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac + k_f;
    dqdc[18] = dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 6] += dqdc[k];
      J[27 * k + 12] += dqdc[k];
      J[27 * k + 17] -= dqdc[k];
    }
  }
  J[708] += dqdT; // dwdot[OH]/dT
  J[714] += dqdT; // dwdot[CH3]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT

  // reaction 67: CH3OH (+M) <=> CH2(S) + H2O (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture;
  // forward
  phi_f = sc[17];
  k_f = 3.121e+18 * exp(-1.017 * tc[0] - (46151.0068542825) * invT);
  dlnkfdT = -1.017 * invT + (46151.0068542825) * invT2;
  // pressure-fall-off
  k_0 = 1.43e+47 * exp(-8.227 * tc[0] - (50028.3415859744) * invT);
  Pr = 1e-6 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -8.227 * invT + (50028.3415859744) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 2.545) * exp(-T / 3290);
  Fcent2 = 2.545 * exp(-T / 47320);
  Fcent3 = exp(-47110 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT =
    log10e / Fcent * (-Fcent1 / 3290 - Fcent2 / 47320 + Fcent3 * 47110 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[5] * sc[14];
  Kc = refC * exp(-g_RT[5] - g_RT[14] + g_RT[17]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[17]) + (h_RT[5] + h_RT[14]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[5] += q;  // H2O
  wdot[14] += q; // CH2(S)
  wdot[17] -= q; // CH3OH
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP == 1) {
    // d()/d[H2O]
    dqdci = -k_r * sc[14];
    J[140] += dqdci; // dwdot[H2O]/d[H2O]
    J[149] += dqdci; // dwdot[CH2(S)]/d[H2O]
    J[152] -= dqdci; // dwdot[CH3OH]/d[H2O]
    // d()/d[CH2(S)]
    dqdci = -k_r * sc[5];
    J[383] += dqdci; // dwdot[H2O]/d[CH2(S)]
    J[392] += dqdci; // dwdot[CH2(S)]/d[CH2(S)]
    J[395] -= dqdci; // dwdot[CH3OH]/d[CH2(S)]
    // d()/d[CH3OH]
    dqdci = +k_f;
    J[464] += dqdci; // dwdot[H2O]/d[CH3OH]
    J[473] += dqdci; // dwdot[CH2(S)]/d[CH3OH]
    J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  } else {
    dqdc[0] = dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = dcdc_fac - k_r * sc[14];
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    dqdc[12] = dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac - k_r * sc[5];
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac + k_f;
    dqdc[18] = dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 5] += dqdc[k];
      J[27 * k + 14] += dqdc[k];
      J[27 * k + 17] -= dqdc[k];
    }
  }
  J[707] += dqdT; // dwdot[H2O]/dT
  J[716] += dqdT; // dwdot[CH2(S)]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT

  // reaction 106: CH2O + H (+M) <=> CH2OH (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha =
    mixture + sc[1] + (6 - 1) * sc[5] + (1.5 - 1) * sc[9] + sc[10] + sc[11];
  // forward
  phi_f = sc[2] * sc[20];
  k_f = 540000 * exp(0.454 * tc[0] - (1811.57999689699) * invT);
  dlnkfdT = 0.454 * invT + (1811.57999689699) * invT2;
  // pressure-fall-off
  k_0 = 1.27e+32 * exp(-4.82 * tc[0] - (3286.00482770482) * invT);
  Pr = 1e-12 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -4.82 * invT + (3286.00482770482) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.7187) * exp(-T / 103);
  Fcent2 = 0.7187 * exp(-T / 1291);
  Fcent3 = exp(-4160 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT =
    log10e / Fcent * (-Fcent1 / 103 - Fcent2 / 1291 + Fcent3 * 4160 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[19];
  Kc = refCinv * exp(g_RT[2] - g_RT[19] + g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[20]) + (h_RT[19]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[2] -= q;  // H
  wdot[19] += q; // CH2OH
  wdot[20] -= q; // CH2O
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2 - 1) * dcdc_fac;
    J[29] -= dqdci; // dwdot[H]/d[H2]
    J[46] += dqdci; // dwdot[CH2OH]/d[H2]
    J[47] -= dqdci; // dwdot[CH2O]/d[H2]
    // d()/d[H]
    dqdci = +k_f * sc[20];
    J[56] -= dqdci; // dwdot[H]/d[H]
    J[73] += dqdci; // dwdot[CH2OH]/d[H]
    J[74] -= dqdci; // dwdot[CH2O]/d[H]
    // d()/d[H2O]
    dqdci = (6 - 1) * dcdc_fac;
    J[137] -= dqdci; // dwdot[H]/d[H2O]
    J[154] += dqdci; // dwdot[CH2OH]/d[H2O]
    J[155] -= dqdci; // dwdot[CH2O]/d[H2O]
    // d()/d[CO]
    dqdci = (1.5 - 1) * dcdc_fac;
    J[245] -= dqdci; // dwdot[H]/d[CO]
    J[262] += dqdci; // dwdot[CH2OH]/d[CO]
    J[263] -= dqdci; // dwdot[CH2O]/d[CO]
    // d()/d[CO2]
    dqdci = (2 - 1) * dcdc_fac;
    J[272] -= dqdci; // dwdot[H]/d[CO2]
    J[289] += dqdci; // dwdot[CH2OH]/d[CO2]
    J[290] -= dqdci; // dwdot[CH2O]/d[CO2]
    // d()/d[CH4]
    dqdci = (2 - 1) * dcdc_fac;
    J[299] -= dqdci; // dwdot[H]/d[CH4]
    J[316] += dqdci; // dwdot[CH2OH]/d[CH4]
    J[317] -= dqdci; // dwdot[CH2O]/d[CH4]
    // d()/d[CH2OH]
    dqdci = -k_r;
    J[515] -= dqdci; // dwdot[H]/d[CH2OH]
    J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
    J[533] -= dqdci; // dwdot[CH2O]/d[CH2OH]
    // d()/d[CH2O]
    dqdci = +k_f * sc[2];
    J[542] -= dqdci; // dwdot[H]/d[CH2O]
    J[559] += dqdci; // dwdot[CH2OH]/d[CH2O]
    J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  } else {
    dqdc[0] = dcdc_fac;
    dqdc[1] = 2 * dcdc_fac;
    dqdc[2] = dcdc_fac + k_f * sc[20];
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 6 * dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = 1.5 * dcdc_fac;
    dqdc[10] = 2 * dcdc_fac;
    dqdc[11] = 2 * dcdc_fac;
    dqdc[12] = dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = dcdc_fac;
    dqdc[19] = dcdc_fac - k_r;
    dqdc[20] = dcdc_fac + k_f * sc[2];
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 2] -= dqdc[k];
      J[27 * k + 19] += dqdc[k];
      J[27 * k + 20] -= dqdc[k];
    }
  }
  J[704] -= dqdT; // dwdot[H]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT

  // reaction 107: CH3O (+M) <=> CH2O + H (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha =
    mixture + sc[1] + (6 - 1) * sc[5] + (1.5 - 1) * sc[9] + sc[10] + sc[11];
  // forward
  phi_f = sc[18];
  k_f = 68000000000000 * exp(-(13169.1801441095) * invT);
  dlnkfdT = (13169.1801441095) * invT2;
  // pressure-fall-off
  k_0 = 1.87e+25 * exp(-3 * tc[0] - (12231.6874957153) * invT);
  Pr = 1e-6 * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = -3 * invT + (12231.6874957153) * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  // Troe form
  logPr = log10(Pr);
  Fcent1 = (1. - 0.9) * exp(-T / 2500);
  Fcent2 = 0.9 * exp(-T / 1300);
  Fcent3 = exp(-1e+99 * invT);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = exp(M_LN10 * logFcent * troe);
  dlogFcentdT =
    log10e / Fcent * (-Fcent1 / 2500 - Fcent2 / 1300 + Fcent3 * 1e+99 * invT2);
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  // reverse
  phi_r = sc[2] * sc[20];
  Kc = refC * exp(-g_RT[2] + g_RT[18] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[18]) + (h_RT[2] + h_RT[20]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[2] += q;  // H
  wdot[18] -= q; // CH3O
  wdot[20] += q; // CH2O
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2 - 1) * dcdc_fac;
    J[29] += dqdci; // dwdot[H]/d[H2]
    J[45] -= dqdci; // dwdot[CH3O]/d[H2]
    J[47] += dqdci; // dwdot[CH2O]/d[H2]
    // d()/d[H]
    dqdci = -k_r * sc[20];
    J[56] += dqdci; // dwdot[H]/d[H]
    J[72] -= dqdci; // dwdot[CH3O]/d[H]
    J[74] += dqdci; // dwdot[CH2O]/d[H]
    // d()/d[H2O]
    dqdci = (6 - 1) * dcdc_fac;
    J[137] += dqdci; // dwdot[H]/d[H2O]
    J[153] -= dqdci; // dwdot[CH3O]/d[H2O]
    J[155] += dqdci; // dwdot[CH2O]/d[H2O]
    // d()/d[CO]
    dqdci = (1.5 - 1) * dcdc_fac;
    J[245] += dqdci; // dwdot[H]/d[CO]
    J[261] -= dqdci; // dwdot[CH3O]/d[CO]
    J[263] += dqdci; // dwdot[CH2O]/d[CO]
    // d()/d[CO2]
    dqdci = (2 - 1) * dcdc_fac;
    J[272] += dqdci; // dwdot[H]/d[CO2]
    J[288] -= dqdci; // dwdot[CH3O]/d[CO2]
    J[290] += dqdci; // dwdot[CH2O]/d[CO2]
    // d()/d[CH4]
    dqdci = (2 - 1) * dcdc_fac;
    J[299] += dqdci; // dwdot[H]/d[CH4]
    J[315] -= dqdci; // dwdot[CH3O]/d[CH4]
    J[317] += dqdci; // dwdot[CH2O]/d[CH4]
    // d()/d[CH3O]
    dqdci = +k_f;
    J[488] += dqdci; // dwdot[H]/d[CH3O]
    J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
    J[506] += dqdci; // dwdot[CH2O]/d[CH3O]
    // d()/d[CH2O]
    dqdci = -k_r * sc[2];
    J[542] += dqdci; // dwdot[H]/d[CH2O]
    J[558] -= dqdci; // dwdot[CH3O]/d[CH2O]
    J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  } else {
    dqdc[0] = dcdc_fac;
    dqdc[1] = 2 * dcdc_fac;
    dqdc[2] = dcdc_fac - k_r * sc[20];
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 6 * dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = 1.5 * dcdc_fac;
    dqdc[10] = 2 * dcdc_fac;
    dqdc[11] = 2 * dcdc_fac;
    dqdc[12] = dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = dcdc_fac + k_f;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac - k_r * sc[2];
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 2] += dqdc[k];
      J[27 * k + 18] -= dqdc[k];
      J[27 * k + 20] += dqdc[k];
    }
  }
  J[704] += dqdT; // dwdot[H]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 20: CO + O (+M) <=> CO2 (+M)
  // a third-body and pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + sc[1] + (12 - 1) * sc[5] + (1.75 - 1) * sc[9] +
          (3.6 - 1) * sc[10];
  // forward
  phi_f = sc[4] * sc[9];
  k_f = 13620 * exp(-(1199.66853127845) * invT);
  dlnkfdT = (1199.66853127845) * invT2;
  // pressure-fall-off
  k_0 = 1.17e+24 * exp(-2.79 * tc[0] - (2108.98104638758) * invT);
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
  phi_r = sc[10];
  Kc = refCinv * exp(g_RT[4] + g_RT[9] - g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[10]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  // update wdot
  wdot[4] -= q;  // O
  wdot[9] -= q;  // CO
  wdot[10] += q; // CO2
  // for convenience
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2 - 1) * dcdc_fac;
    J[31] -= dqdci; // dwdot[O]/d[H2]
    J[36] -= dqdci; // dwdot[CO]/d[H2]
    J[37] += dqdci; // dwdot[CO2]/d[H2]
    // d()/d[O]
    dqdci = +k_f * sc[9];
    J[112] -= dqdci; // dwdot[O]/d[O]
    J[117] -= dqdci; // dwdot[CO]/d[O]
    J[118] += dqdci; // dwdot[CO2]/d[O]
    // d()/d[H2O]
    dqdci = (12 - 1) * dcdc_fac;
    J[139] -= dqdci; // dwdot[O]/d[H2O]
    J[144] -= dqdci; // dwdot[CO]/d[H2O]
    J[145] += dqdci; // dwdot[CO2]/d[H2O]
    // d()/d[CO]
    dqdci = (1.75 - 1) * dcdc_fac + k_f * sc[4];
    J[247] -= dqdci; // dwdot[O]/d[CO]
    J[252] -= dqdci; // dwdot[CO]/d[CO]
    J[253] += dqdci; // dwdot[CO2]/d[CO]
    // d()/d[CO2]
    dqdci = (3.6 - 1) * dcdc_fac - k_r;
    J[274] -= dqdci; // dwdot[O]/d[CO2]
    J[279] -= dqdci; // dwdot[CO]/d[CO2]
    J[280] += dqdci; // dwdot[CO2]/d[CO2]
  } else {
    dqdc[0] = dcdc_fac;
    dqdc[1] = 2 * dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac + k_f * sc[9];
    dqdc[5] = 12 * dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = 1.75 * dcdc_fac + k_f * sc[4];
    dqdc[10] = 3.6 * dcdc_fac - k_r;
    dqdc[11] = dcdc_fac;
    dqdc[12] = dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 4] -= dqdc[k];
      J[27 * k + 9] -= dqdc[k];
      J[27 * k + 10] += dqdc[k];
    }
  }
  J[706] -= dqdT; // dwdot[O]/dT
  J[711] -= dqdT; // dwdot[CO]/dT
  J[712] += dqdT; // dwdot[CO2]/dT

  // reaction 3: H + OH + M <=> H2O + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (0.73 - 1) * sc[1] + (3.65 - 1) * sc[5] + sc[11];
  // forward
  phi_f = sc[2] * sc[6];
  k_f = 35000000000 * exp(-2 * tc[0]);
  dlnkfdT = -2 * invT;
  // reverse
  phi_r = sc[5];
  Kc = refCinv * exp(g_RT[2] - g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[5]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // H
  wdot[5] += q; // H2O
  wdot[6] -= q; // OH
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (0.73 - 1) * q_nocor;
    J[29] -= dqdci; // dwdot[H]/d[H2]
    J[32] += dqdci; // dwdot[H2O]/d[H2]
    J[33] -= dqdci; // dwdot[OH]/d[H2]
    // d()/d[H]
    dqdci = +k_f * sc[6];
    J[56] -= dqdci; // dwdot[H]/d[H]
    J[59] += dqdci; // dwdot[H2O]/d[H]
    J[60] -= dqdci; // dwdot[OH]/d[H]
    // d()/d[H2O]
    dqdci = (3.65 - 1) * q_nocor - k_r;
    J[137] -= dqdci; // dwdot[H]/d[H2O]
    J[140] += dqdci; // dwdot[H2O]/d[H2O]
    J[141] -= dqdci; // dwdot[OH]/d[H2O]
    // d()/d[OH]
    dqdci = +k_f * sc[2];
    J[164] -= dqdci; // dwdot[H]/d[OH]
    J[167] += dqdci; // dwdot[H2O]/d[OH]
    J[168] -= dqdci; // dwdot[OH]/d[OH]
    // d()/d[CH4]
    dqdci = (2 - 1) * q_nocor;
    J[299] -= dqdci; // dwdot[H]/d[CH4]
    J[302] += dqdci; // dwdot[H2O]/d[CH4]
    J[303] -= dqdci; // dwdot[OH]/d[CH4]
  } else {
    dqdc[0] = q_nocor;
    dqdc[1] = 0.73 * q_nocor;
    dqdc[2] = q_nocor + k_f * sc[6];
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor;
    dqdc[5] = 3.65 * q_nocor - k_r;
    dqdc[6] = q_nocor + k_f * sc[2];
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = 2 * q_nocor;
    dqdc[12] = q_nocor;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = q_nocor;
    dqdc[21] = q_nocor;
    dqdc[22] = q_nocor;
    dqdc[23] = q_nocor;
    dqdc[24] = q_nocor;
    dqdc[25] = q_nocor;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 2] -= dqdc[k];
      J[27 * k + 5] += dqdc[k];
      J[27 * k + 6] -= dqdc[k];
    }
  }
  J[704] -= dqdT; // dwdot[H]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT

  // reaction 5: H + O + M <=> OH + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture + (2.5 - 1) * sc[1] + (12 - 1) * sc[5] + (1.5 - 1) * sc[9] +
          sc[10] + sc[11];
  // forward
  phi_f = sc[2] * sc[4];
  k_f = 4714000 * exp(-1 * tc[0]);
  dlnkfdT = -1 * invT;
  // reverse
  phi_r = sc[6];
  Kc = refCinv * exp(g_RT[2] + g_RT[4] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (h_RT[6]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // H
  wdot[4] -= q; // O
  wdot[6] += q; // OH
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2.5 - 1) * q_nocor;
    J[29] -= dqdci; // dwdot[H]/d[H2]
    J[31] -= dqdci; // dwdot[O]/d[H2]
    J[33] += dqdci; // dwdot[OH]/d[H2]
    // d()/d[H]
    dqdci = +k_f * sc[4];
    J[56] -= dqdci; // dwdot[H]/d[H]
    J[58] -= dqdci; // dwdot[O]/d[H]
    J[60] += dqdci; // dwdot[OH]/d[H]
    // d()/d[O]
    dqdci = +k_f * sc[2];
    J[110] -= dqdci; // dwdot[H]/d[O]
    J[112] -= dqdci; // dwdot[O]/d[O]
    J[114] += dqdci; // dwdot[OH]/d[O]
    // d()/d[H2O]
    dqdci = (12 - 1) * q_nocor;
    J[137] -= dqdci; // dwdot[H]/d[H2O]
    J[139] -= dqdci; // dwdot[O]/d[H2O]
    J[141] += dqdci; // dwdot[OH]/d[H2O]
    // d()/d[OH]
    dqdci = -k_r;
    J[164] -= dqdci; // dwdot[H]/d[OH]
    J[166] -= dqdci; // dwdot[O]/d[OH]
    J[168] += dqdci; // dwdot[OH]/d[OH]
    // d()/d[CO]
    dqdci = (1.5 - 1) * q_nocor;
    J[245] -= dqdci; // dwdot[H]/d[CO]
    J[247] -= dqdci; // dwdot[O]/d[CO]
    J[249] += dqdci; // dwdot[OH]/d[CO]
    // d()/d[CO2]
    dqdci = (2 - 1) * q_nocor;
    J[272] -= dqdci; // dwdot[H]/d[CO2]
    J[274] -= dqdci; // dwdot[O]/d[CO2]
    J[276] += dqdci; // dwdot[OH]/d[CO2]
    // d()/d[CH4]
    dqdci = (2 - 1) * q_nocor;
    J[299] -= dqdci; // dwdot[H]/d[CH4]
    J[301] -= dqdci; // dwdot[O]/d[CH4]
    J[303] += dqdci; // dwdot[OH]/d[CH4]
  } else {
    dqdc[0] = q_nocor;
    dqdc[1] = 2.5 * q_nocor;
    dqdc[2] = q_nocor + k_f * sc[4];
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor + k_f * sc[2];
    dqdc[5] = 12 * q_nocor;
    dqdc[6] = q_nocor - k_r;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = 1.5 * q_nocor;
    dqdc[10] = 2 * q_nocor;
    dqdc[11] = 2 * q_nocor;
    dqdc[12] = q_nocor;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = q_nocor;
    dqdc[21] = q_nocor;
    dqdc[22] = q_nocor;
    dqdc[23] = q_nocor;
    dqdc[24] = q_nocor;
    dqdc[25] = q_nocor;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 2] -= dqdc[k];
      J[27 * k + 4] -= dqdc[k];
      J[27 * k + 6] += dqdc[k];
    }
  }
  J[704] -= dqdT; // dwdot[H]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT

  // reaction 97: HCO + M <=> CO + H + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha =
    mixture + sc[1] + (6 - 1) * sc[5] + (1.5 - 1) * sc[9] + sc[10] + sc[11];
  // forward
  phi_f = sc[21];
  k_f = 570000 * exp(0.66 * tc[0] - (7482.83182051618) * invT);
  dlnkfdT = 0.66 * invT + (7482.83182051618) * invT2;
  // reverse
  phi_r = sc[2] * sc[9];
  Kc = refC * exp(-g_RT[2] - g_RT[9] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[21]) + (h_RT[2] + h_RT[9]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[9] += q;  // CO
  wdot[21] -= q; // HCO
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  if (consP == 1) {
    // d()/d[H2]
    dqdci = (2 - 1) * q_nocor;
    J[29] += dqdci; // dwdot[H]/d[H2]
    J[36] += dqdci; // dwdot[CO]/d[H2]
    J[48] -= dqdci; // dwdot[HCO]/d[H2]
    // d()/d[H]
    dqdci = -k_r * sc[9];
    J[56] += dqdci; // dwdot[H]/d[H]
    J[63] += dqdci; // dwdot[CO]/d[H]
    J[75] -= dqdci; // dwdot[HCO]/d[H]
    // d()/d[H2O]
    dqdci = (6 - 1) * q_nocor;
    J[137] += dqdci; // dwdot[H]/d[H2O]
    J[144] += dqdci; // dwdot[CO]/d[H2O]
    J[156] -= dqdci; // dwdot[HCO]/d[H2O]
    // d()/d[CO]
    dqdci = (1.5 - 1) * q_nocor - k_r * sc[2];
    J[245] += dqdci; // dwdot[H]/d[CO]
    J[252] += dqdci; // dwdot[CO]/d[CO]
    J[264] -= dqdci; // dwdot[HCO]/d[CO]
    // d()/d[CO2]
    dqdci = (2 - 1) * q_nocor;
    J[272] += dqdci; // dwdot[H]/d[CO2]
    J[279] += dqdci; // dwdot[CO]/d[CO2]
    J[291] -= dqdci; // dwdot[HCO]/d[CO2]
    // d()/d[CH4]
    dqdci = (2 - 1) * q_nocor;
    J[299] += dqdci; // dwdot[H]/d[CH4]
    J[306] += dqdci; // dwdot[CO]/d[CH4]
    J[318] -= dqdci; // dwdot[HCO]/d[CH4]
    // d()/d[HCO]
    dqdci = +k_f;
    J[569] += dqdci; // dwdot[H]/d[HCO]
    J[576] += dqdci; // dwdot[CO]/d[HCO]
    J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  } else {
    dqdc[0] = q_nocor;
    dqdc[1] = 2 * q_nocor;
    dqdc[2] = q_nocor - k_r * sc[9];
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor;
    dqdc[5] = 6 * q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = 1.5 * q_nocor - k_r * sc[2];
    dqdc[10] = 2 * q_nocor;
    dqdc[11] = 2 * q_nocor;
    dqdc[12] = q_nocor;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = q_nocor;
    dqdc[21] = q_nocor + k_f;
    dqdc[22] = q_nocor;
    dqdc[23] = q_nocor;
    dqdc[24] = q_nocor;
    dqdc[25] = q_nocor;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 2] += dqdc[k];
      J[27 * k + 9] += dqdc[k];
      J[27 * k + 21] -= dqdc[k];
    }
  }
  J[704] += dqdT; // dwdot[H]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 108: O + OH + M <=> HO2 + M
  // a third-body and non-pressure-fall-off reaction
  // 3-body correction factor
  alpha = mixture;
  // forward
  phi_f = sc[4] * sc[6];
  k_f = 10000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[8];
  Kc = refCinv * exp(g_RT[4] + g_RT[6] - g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[8]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q; // O
  wdot[6] -= q; // OH
  wdot[8] += q; // HO2
  // for convenience
  k_f *= alpha;
  k_r *= alpha;
  if (consP == 1) {
    // d()/d[O]
    dqdci = +k_f * sc[6];
    J[112] -= dqdci; // dwdot[O]/d[O]
    J[114] -= dqdci; // dwdot[OH]/d[O]
    J[116] += dqdci; // dwdot[HO2]/d[O]
    // d()/d[OH]
    dqdci = +k_f * sc[4];
    J[166] -= dqdci; // dwdot[O]/d[OH]
    J[168] -= dqdci; // dwdot[OH]/d[OH]
    J[170] += dqdci; // dwdot[HO2]/d[OH]
    // d()/d[HO2]
    dqdci = -k_r;
    J[220] -= dqdci; // dwdot[O]/d[HO2]
    J[222] -= dqdci; // dwdot[OH]/d[HO2]
    J[224] += dqdci; // dwdot[HO2]/d[HO2]
  } else {
    dqdc[0] = q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor + k_f * sc[6];
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor + k_f * sc[4];
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor - k_r;
    dqdc[9] = q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    dqdc[12] = q_nocor;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = q_nocor;
    dqdc[21] = q_nocor;
    dqdc[22] = q_nocor;
    dqdc[23] = q_nocor;
    dqdc[24] = q_nocor;
    dqdc[25] = q_nocor;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 4] -= dqdc[k];
      J[27 * k + 6] -= dqdc[k];
      J[27 * k + 8] += dqdc[k];
    }
  }
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[710] += dqdT; // dwdot[HO2]/dT

  // reaction 120: CH3 + O + M => CH3O + M
  // a third-body and non-pressure-fall-off reaction
  // Irreversible reaction in _ajac_reaction may not work
  // 3-body correction factor
  alpha = mixture;
  // forward
  phi_f = sc[4] * sc[12];
  k_f = 50000;
  dlnkfdT = 0.0;
  // rate of progress
  q_nocor = k_f * phi_f;
  q = alpha * q_nocor;
  dqdT = alpha * dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[4] -= q;  // O
  wdot[12] -= q; // CH3
  wdot[18] += q; // CH3O
  // for convenience
  k_f *= alpha;
  k_r = 0.0;
  if (consP == 1) {
    // d()/d[O]
    dqdci = +k_f * sc[12];
    J[112] -= dqdci; // dwdot[O]/d[O]
    J[120] -= dqdci; // dwdot[CH3]/d[O]
    J[126] += dqdci; // dwdot[CH3O]/d[O]
    // d()/d[CH3]
    dqdci = +k_f * sc[4];
    J[328] -= dqdci; // dwdot[O]/d[CH3]
    J[336] -= dqdci; // dwdot[CH3]/d[CH3]
    J[342] += dqdci; // dwdot[CH3O]/d[CH3]
  } else {
    dqdc[0] = q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor + k_f * sc[12];
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    dqdc[12] = q_nocor + k_f * sc[4];
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = q_nocor;
    dqdc[21] = q_nocor;
    dqdc[22] = q_nocor;
    dqdc[23] = q_nocor;
    dqdc[24] = q_nocor;
    dqdc[25] = q_nocor;
    for (int k = 0; k < 26; k++) {
      J[27 * k + 4] -= dqdc[k];
      J[27 * k + 12] -= dqdc[k];
      J[27 * k + 18] += dqdc[k];
    }
  }
  J[706] -= dqdT; // dwdot[O]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 131: E + O2 + O2 <=> O2- + O2
  // a third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[3] * sc[3]) * sc[25];
  k_f = 1520000000 * exp(-1 * tc[0] - (600.156324305341) * invT);
  dlnkfdT = -1 * invT + (600.156324305341) * invT2;
  // reverse
  phi_r = sc[3] * sc[24];
  Kc = refCinv * exp(g_RT[3] - g_RT[24] + g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT =
    invT * (-(2.000000 * h_RT[3] + h_RT[25]) + (h_RT[3] + h_RT[24]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[24] += q; // O2-
  wdot[25] -= q; // E
  // d()/d[O2]
  dqdci = +k_f * 2.000000 * sc[3] * sc[25] - k_r * sc[24];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[105] += dqdci; // dwdot[O2-]/d[O2]
  J[106] -= dqdci; // dwdot[E]/d[O2]
  // d()/d[O2-]
  dqdci = -k_r * sc[3];
  J[651] -= dqdci; // dwdot[O2]/d[O2-]
  J[672] += dqdci; // dwdot[O2-]/d[O2-]
  J[673] -= dqdci; // dwdot[E]/d[O2-]
  // d()/d[E]
  dqdci = +k_f * sc[3] * sc[3];
  J[678] -= dqdci; // dwdot[O2]/d[E]
  J[699] += dqdci; // dwdot[O2-]/d[E]
  J[700] -= dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[726] += dqdT; // dwdot[O2-]/dT
  J[727] -= dqdT; // dwdot[E]/dT

  // reaction 132: E + O2 + H2O <=> O2- + H2O
  // a third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[5] * sc[25];
  k_f = 5080000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[5] * sc[24];
  Kc = refCinv * exp(g_RT[3] - g_RT[24] + g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT =
    invT * (-(h_RT[3] + h_RT[5] + h_RT[25]) + (h_RT[5] + h_RT[24]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[24] += q; // O2-
  wdot[25] -= q; // E
  // d()/d[O2]
  dqdci = +k_f * sc[5] * sc[25];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[105] += dqdci; // dwdot[O2-]/d[O2]
  J[106] -= dqdci; // dwdot[E]/d[O2]
  // d()/d[H2O]
  dqdci = +k_f * sc[3] * sc[25] - k_r * sc[24];
  J[138] -= dqdci; // dwdot[O2]/d[H2O]
  J[159] += dqdci; // dwdot[O2-]/d[H2O]
  J[160] -= dqdci; // dwdot[E]/d[H2O]
  // d()/d[O2-]
  dqdci = -k_r * sc[5];
  J[651] -= dqdci; // dwdot[O2]/d[O2-]
  J[672] += dqdci; // dwdot[O2-]/d[O2-]
  J[673] -= dqdci; // dwdot[E]/d[O2-]
  // d()/d[E]
  dqdci = +k_f * sc[3] * sc[5];
  J[678] -= dqdci; // dwdot[O2]/d[E]
  J[699] += dqdci; // dwdot[O2-]/d[E]
  J[700] -= dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[726] += dqdT; // dwdot[O2-]/dT
  J[727] -= dqdT; // dwdot[E]/dT

  // reaction 133: E + O2 + N2 <=> O2- + N2
  // a third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[3] * sc[25];
  k_f = 3590000000 * exp(-2 * tc[0] - (69.7558942138502) * invT);
  dlnkfdT = -2 * invT + (69.7558942138502) * invT2;
  // reverse
  phi_r = sc[0] * sc[24];
  Kc = refCinv * exp(g_RT[3] - g_RT[24] + g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT =
    invT * (-(h_RT[0] + h_RT[3] + h_RT[25]) + (h_RT[0] + h_RT[24]) + 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[24] += q; // O2-
  wdot[25] -= q; // E
  // d()/d[N2]
  dqdci = +k_f * sc[3] * sc[25] - k_r * sc[24];
  J[3] -= dqdci;  // dwdot[O2]/d[N2]
  J[24] += dqdci; // dwdot[O2-]/d[N2]
  J[25] -= dqdci; // dwdot[E]/d[N2]
  // d()/d[O2]
  dqdci = +k_f * sc[0] * sc[25];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[105] += dqdci; // dwdot[O2-]/d[O2]
  J[106] -= dqdci; // dwdot[E]/d[O2]
  // d()/d[O2-]
  dqdci = -k_r * sc[0];
  J[651] -= dqdci; // dwdot[O2]/d[O2-]
  J[672] += dqdci; // dwdot[O2-]/d[O2-]
  J[673] -= dqdci; // dwdot[E]/d[O2-]
  // d()/d[E]
  dqdci = +k_f * sc[0] * sc[3];
  J[678] -= dqdci; // dwdot[O2]/d[E]
  J[699] += dqdci; // dwdot[O2-]/d[E]
  J[700] -= dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[726] += dqdT; // dwdot[O2-]/dT
  J[727] -= dqdT; // dwdot[E]/dT

  // reaction 0: H2 + O <=> H + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[4];
  k_f = 0.0508 * exp(2.67 * tc[0] - (3166.2392612433) * invT);
  dlnkfdT = 2.67 * invT + (3166.2392612433) * invT2;
  // reverse
  phi_r = sc[2] * sc[6];
  Kc = exp(g_RT[1] - g_RT[2] + g_RT[4] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[4]) + (h_RT[2] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q; // H2
  wdot[2] += q; // H
  wdot[4] -= q; // O
  wdot[6] += q; // OH
  // d()/d[H2]
  dqdci = +k_f * sc[4];
  J[28] -= dqdci; // dwdot[H2]/d[H2]
  J[29] += dqdci; // dwdot[H]/d[H2]
  J[31] -= dqdci; // dwdot[O]/d[H2]
  J[33] += dqdci; // dwdot[OH]/d[H2]
  // d()/d[H]
  dqdci = -k_r * sc[6];
  J[55] -= dqdci; // dwdot[H2]/d[H]
  J[56] += dqdci; // dwdot[H]/d[H]
  J[58] -= dqdci; // dwdot[O]/d[H]
  J[60] += dqdci; // dwdot[OH]/d[H]
  // d()/d[O]
  dqdci = +k_f * sc[1];
  J[109] -= dqdci; // dwdot[H2]/d[O]
  J[110] += dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[2];
  J[163] -= dqdci; // dwdot[H2]/d[OH]
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  // d()/dT
  J[703] -= dqdT; // dwdot[H2]/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT

  // reaction 1: H2 + OH <=> H + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[6];
  k_f = 43800000 * exp(-(3517.48449397499) * invT);
  dlnkfdT = (3517.48449397499) * invT2;
  // reverse
  phi_r = sc[2] * sc[5];
  Kc = exp(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[2] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q; // H2
  wdot[2] += q; // H
  wdot[5] += q; // H2O
  wdot[6] -= q; // OH
  // d()/d[H2]
  dqdci = +k_f * sc[6];
  J[28] -= dqdci; // dwdot[H2]/d[H2]
  J[29] += dqdci; // dwdot[H]/d[H2]
  J[32] += dqdci; // dwdot[H2O]/d[H2]
  J[33] -= dqdci; // dwdot[OH]/d[H2]
  // d()/d[H]
  dqdci = -k_r * sc[5];
  J[55] -= dqdci; // dwdot[H2]/d[H]
  J[56] += dqdci; // dwdot[H]/d[H]
  J[59] += dqdci; // dwdot[H2O]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  // d()/d[H2O]
  dqdci = -k_r * sc[2];
  J[136] -= dqdci; // dwdot[H2]/d[H2O]
  J[137] += dqdci; // dwdot[H]/d[H2O]
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[1];
  J[163] -= dqdci; // dwdot[H2]/d[OH]
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  // d()/dT
  J[703] -= dqdT; // dwdot[H2]/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT

  // reaction 2: H + O2 <=> O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[3];
  k_f = 114000000 * exp(-(7692.16995349094) * invT);
  dlnkfdT = (7692.16995349094) * invT2;
  // reverse
  phi_r = sc[4] * sc[6];
  Kc = exp(g_RT[2] + g_RT[3] - g_RT[4] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[4] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // H
  wdot[3] -= q; // O2
  wdot[4] += q; // O
  wdot[6] += q; // OH
  // d()/d[H]
  dqdci = +k_f * sc[3];
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[57] -= dqdci; // dwdot[O2]/d[H]
  J[58] += dqdci; // dwdot[O]/d[H]
  J[60] += dqdci; // dwdot[OH]/d[H]
  // d()/d[O2]
  dqdci = +k_f * sc[2];
  J[83] -= dqdci; // dwdot[H]/d[O2]
  J[84] -= dqdci; // dwdot[O2]/d[O2]
  J[85] += dqdci; // dwdot[O]/d[O2]
  J[87] += dqdci; // dwdot[OH]/d[O2]
  // d()/d[O]
  dqdci = -k_r * sc[6];
  J[110] -= dqdci; // dwdot[H]/d[O]
  J[111] -= dqdci; // dwdot[O2]/d[O]
  J[112] += dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[4];
  J[164] -= dqdci; // dwdot[H]/d[OH]
  J[165] -= dqdci; // dwdot[O2]/d[OH]
  J[166] += dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  // d()/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[706] += dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT

  // reaction 4: H2O + O <=> 2 OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[5];
  k_f = 67 * exp(1.704 * tc[0] - (7541.60752708217) * invT);
  dlnkfdT = 1.704 * invT + (7541.60752708217) * invT2;
  // reverse
  phi_r = (sc[6] * sc[6]);
  Kc = exp(g_RT[4] + g_RT[5] - 2.000000 * g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[5]) + (2.000000 * h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;     // O
  wdot[5] -= q;     // H2O
  wdot[6] += 2 * q; // OH
  // d()/d[O]
  dqdci = +k_f * sc[5];
  J[112] -= dqdci;     // dwdot[O]/d[O]
  J[113] -= dqdci;     // dwdot[H2O]/d[O]
  J[114] += 2 * dqdci; // dwdot[OH]/d[O]
  // d()/d[H2O]
  dqdci = +k_f * sc[4];
  J[139] -= dqdci;     // dwdot[O]/d[H2O]
  J[140] -= dqdci;     // dwdot[H2O]/d[H2O]
  J[141] += 2 * dqdci; // dwdot[OH]/d[H2O]
  // d()/d[OH]
  dqdci = -k_r * 2.000000 * sc[6];
  J[166] -= dqdci;     // dwdot[O]/d[OH]
  J[167] -= dqdci;     // dwdot[H2O]/d[OH]
  J[168] += 2 * dqdci; // dwdot[OH]/d[OH]
  // d()/dT
  J[706] -= dqdT;     // dwdot[O]/dT
  J[707] -= dqdT;     // dwdot[H2O]/dT
  J[708] += 2 * dqdT; // dwdot[OH]/dT

  // reaction 7: H + H2O2 <=> H2O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[7];
  k_f = 24100000 * exp(-(1997.77016324474) * invT);
  dlnkfdT = (1997.77016324474) * invT2;
  // reverse
  phi_r = sc[5] * sc[6];
  Kc = exp(g_RT[2] - g_RT[5] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[5] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q; // H
  wdot[5] += q; // H2O
  wdot[6] += q; // OH
  wdot[7] -= q; // H2O2
  // d()/d[H]
  dqdci = +k_f * sc[7];
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[59] += dqdci; // dwdot[H2O]/d[H]
  J[60] += dqdci; // dwdot[OH]/d[H]
  J[61] -= dqdci; // dwdot[H2O2]/d[H]
  // d()/d[H2O]
  dqdci = -k_r * sc[6];
  J[137] -= dqdci; // dwdot[H]/d[H2O]
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] += dqdci; // dwdot[OH]/d[H2O]
  J[142] -= dqdci; // dwdot[H2O2]/d[H2O]
  // d()/d[OH]
  dqdci = -k_r * sc[5];
  J[164] -= dqdci; // dwdot[H]/d[OH]
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[169] -= dqdci; // dwdot[H2O2]/d[OH]
  // d()/d[H2O2]
  dqdci = +k_f * sc[2];
  J[191] -= dqdci; // dwdot[H]/d[H2O2]
  J[194] += dqdci; // dwdot[H2O]/d[H2O2]
  J[195] += dqdci; // dwdot[OH]/d[H2O2]
  J[196] -= dqdci; // dwdot[H2O2]/d[H2O2]
  // d()/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[709] -= dqdT; // dwdot[H2O2]/dT

  // reaction 8: H + H2O2 <=> H2 + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[7];
  k_f = 21500 * exp(1 * tc[0] - (3019.29999482832) * invT);
  dlnkfdT = 1 * invT + (3019.29999482832) * invT2;
  // reverse
  phi_r = sc[1] * sc[8];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[1] + h_RT[8]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q; // H2
  wdot[2] -= q; // H
  wdot[7] -= q; // H2O2
  wdot[8] += q; // HO2
  // d()/d[H2]
  dqdci = -k_r * sc[8];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[34] -= dqdci; // dwdot[H2O2]/d[H2]
  J[35] += dqdci; // dwdot[HO2]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[7];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[61] -= dqdci; // dwdot[H2O2]/d[H]
  J[62] += dqdci; // dwdot[HO2]/d[H]
  // d()/d[H2O2]
  dqdci = +k_f * sc[2];
  J[190] += dqdci; // dwdot[H2]/d[H2O2]
  J[191] -= dqdci; // dwdot[H]/d[H2O2]
  J[196] -= dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] += dqdci; // dwdot[HO2]/d[H2O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[1];
  J[217] += dqdci; // dwdot[H2]/d[HO2]
  J[218] -= dqdci; // dwdot[H]/d[HO2]
  J[223] -= dqdci; // dwdot[H2O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[709] -= dqdT; // dwdot[H2O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT

  // reaction 9: H2O2 + O <=> HO2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[7];
  k_f = 9.55 * exp(2 * tc[0] - (1997.77016324474) * invT);
  dlnkfdT = 2 * invT + (1997.77016324474) * invT2;
  // reverse
  phi_r = sc[6] * sc[8];
  Kc = exp(g_RT[4] - g_RT[6] + g_RT[7] - g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[6] + h_RT[8]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q; // O
  wdot[6] += q; // OH
  wdot[7] -= q; // H2O2
  wdot[8] += q; // HO2
  // d()/d[O]
  dqdci = +k_f * sc[7];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[115] -= dqdci; // dwdot[H2O2]/d[O]
  J[116] += dqdci; // dwdot[HO2]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[8];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[169] -= dqdci; // dwdot[H2O2]/d[OH]
  J[170] += dqdci; // dwdot[HO2]/d[OH]
  // d()/d[H2O2]
  dqdci = +k_f * sc[4];
  J[193] -= dqdci; // dwdot[O]/d[H2O2]
  J[195] += dqdci; // dwdot[OH]/d[H2O2]
  J[196] -= dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] += dqdci; // dwdot[HO2]/d[H2O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[6];
  J[220] -= dqdci; // dwdot[O]/d[HO2]
  J[222] += dqdci; // dwdot[OH]/d[HO2]
  J[223] -= dqdci; // dwdot[H2O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[709] -= dqdT; // dwdot[H2O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT

  // reaction 10: H2O2 + OH <=> H2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[7];
  k_f = 1740000 * exp(-(160.022899725901) * invT);
  dlnkfdT = (160.022899725901) * invT2;
  // reverse
  phi_r = sc[5] * sc[8];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[7] - g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[7]) + (h_RT[5] + h_RT[8]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q; // H2O
  wdot[6] -= q; // OH
  wdot[7] -= q; // H2O2
  wdot[8] += q; // HO2
  // d()/d[H2O]
  dqdci = -k_r * sc[8];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[142] -= dqdci; // dwdot[H2O2]/d[H2O]
  J[143] += dqdci; // dwdot[HO2]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[7];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[169] -= dqdci; // dwdot[H2O2]/d[OH]
  J[170] += dqdci; // dwdot[HO2]/d[OH]
  // d()/d[H2O2]
  dqdci = +k_f * sc[6];
  J[194] += dqdci; // dwdot[H2O]/d[H2O2]
  J[195] -= dqdci; // dwdot[OH]/d[H2O2]
  J[196] -= dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] += dqdci; // dwdot[HO2]/d[H2O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[5];
  J[221] += dqdci; // dwdot[H2O]/d[HO2]
  J[222] -= dqdci; // dwdot[OH]/d[HO2]
  J[223] -= dqdci; // dwdot[H2O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[709] -= dqdT; // dwdot[H2O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT

  // reaction 11: H2O2 + OH <=> H2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[7];
  k_f = 75900000 * exp(-(3657.88194373451) * invT);
  dlnkfdT = (3657.88194373451) * invT2;
  // reverse
  phi_r = sc[5] * sc[8];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[7] - g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[7]) + (h_RT[5] + h_RT[8]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q; // H2O
  wdot[6] -= q; // OH
  wdot[7] -= q; // H2O2
  wdot[8] += q; // HO2
  // d()/d[H2O]
  dqdci = -k_r * sc[8];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[142] -= dqdci; // dwdot[H2O2]/d[H2O]
  J[143] += dqdci; // dwdot[HO2]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[7];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[169] -= dqdci; // dwdot[H2O2]/d[OH]
  J[170] += dqdci; // dwdot[HO2]/d[OH]
  // d()/d[H2O2]
  dqdci = +k_f * sc[6];
  J[194] += dqdci; // dwdot[H2O]/d[H2O2]
  J[195] -= dqdci; // dwdot[OH]/d[H2O2]
  J[196] -= dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] += dqdci; // dwdot[HO2]/d[H2O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[5];
  J[221] += dqdci; // dwdot[H2O]/d[HO2]
  J[222] -= dqdci; // dwdot[OH]/d[HO2]
  J[223] -= dqdci; // dwdot[H2O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[709] -= dqdT; // dwdot[H2O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT

  // reaction 12: H + HO2 <=> 2 OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[8];
  k_f = 70790000 * exp(-(148.448916412392) * invT);
  dlnkfdT = (148.448916412392) * invT2;
  // reverse
  phi_r = (sc[6] * sc[6]);
  Kc = exp(g_RT[2] - 2.000000 * g_RT[6] + g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (2.000000 * h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q;     // H
  wdot[6] += 2 * q; // OH
  wdot[8] -= q;     // HO2
  // d()/d[H]
  dqdci = +k_f * sc[8];
  J[56] -= dqdci;     // dwdot[H]/d[H]
  J[60] += 2 * dqdci; // dwdot[OH]/d[H]
  J[62] -= dqdci;     // dwdot[HO2]/d[H]
  // d()/d[OH]
  dqdci = -k_r * 2.000000 * sc[6];
  J[164] -= dqdci;     // dwdot[H]/d[OH]
  J[168] += 2 * dqdci; // dwdot[OH]/d[OH]
  J[170] -= dqdci;     // dwdot[HO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[2];
  J[218] -= dqdci;     // dwdot[H]/d[HO2]
  J[222] += 2 * dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci;     // dwdot[HO2]/d[HO2]
  // d()/dT
  J[704] -= dqdT;     // dwdot[H]/dT
  J[708] += 2 * dqdT; // dwdot[OH]/dT
  J[710] -= dqdT;     // dwdot[HO2]/dT

  // reaction 13: H + HO2 <=> H2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[8];
  k_f = 11402 * exp(1.083 * tc[0] - (278.671325189338) * invT);
  dlnkfdT = 1.083 * invT + (278.671325189338) * invT2;
  // reverse
  phi_r = sc[1] * sc[3];
  Kc = exp(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[1] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q; // H2
  wdot[2] -= q; // H
  wdot[3] += q; // O2
  wdot[8] -= q; // HO2
  // d()/d[H2]
  dqdci = -k_r * sc[3];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[30] += dqdci; // dwdot[O2]/d[H2]
  J[35] -= dqdci; // dwdot[HO2]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[8];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[57] += dqdci; // dwdot[O2]/d[H]
  J[62] -= dqdci; // dwdot[HO2]/d[H]
  // d()/d[O2]
  dqdci = -k_r * sc[1];
  J[82] += dqdci; // dwdot[H2]/d[O2]
  J[83] -= dqdci; // dwdot[H]/d[O2]
  J[84] += dqdci; // dwdot[O2]/d[O2]
  J[89] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[2];
  J[217] += dqdci; // dwdot[H2]/d[HO2]
  J[218] -= dqdci; // dwdot[H]/d[HO2]
  J[219] += dqdci; // dwdot[O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[705] += dqdT; // dwdot[O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT

  // reaction 14: HO2 + O <=> O2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[8];
  k_f = 32500000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[3] * sc[6];
  Kc = exp(-g_RT[3] + g_RT[4] - g_RT[6] + g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[8]) + (h_RT[3] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q; // O2
  wdot[4] -= q; // O
  wdot[6] += q; // OH
  wdot[8] -= q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[6];
  J[84] += dqdci; // dwdot[O2]/d[O2]
  J[85] -= dqdci; // dwdot[O]/d[O2]
  J[87] += dqdci; // dwdot[OH]/d[O2]
  J[89] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[O]
  dqdci = +k_f * sc[8];
  J[111] += dqdci; // dwdot[O2]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[116] -= dqdci; // dwdot[HO2]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[3];
  J[165] += dqdci; // dwdot[O2]/d[OH]
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[170] -= dqdci; // dwdot[HO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[4];
  J[219] += dqdci; // dwdot[O2]/d[HO2]
  J[220] -= dqdci; // dwdot[O]/d[HO2]
  J[222] += dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[705] += dqdT; // dwdot[O2]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT

  // reaction 15: HO2 + OH <=> H2O + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[8];
  k_f = 7000000 * exp(-(-549.995687057926) * invT);
  dlnkfdT = (-549.995687057926) * invT2;
  // reverse
  phi_r = sc[3] * sc[5];
  Kc = exp(-g_RT[3] - g_RT[5] + g_RT[6] + g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[3] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q; // O2
  wdot[5] += q; // H2O
  wdot[6] -= q; // OH
  wdot[8] -= q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[5];
  J[84] += dqdci; // dwdot[O2]/d[O2]
  J[86] += dqdci; // dwdot[H2O]/d[O2]
  J[87] -= dqdci; // dwdot[OH]/d[O2]
  J[89] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[H2O]
  dqdci = -k_r * sc[3];
  J[138] += dqdci; // dwdot[O2]/d[H2O]
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[143] -= dqdci; // dwdot[HO2]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[8];
  J[165] += dqdci; // dwdot[O2]/d[OH]
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[170] -= dqdci; // dwdot[HO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[6];
  J[219] += dqdci; // dwdot[O2]/d[HO2]
  J[221] += dqdci; // dwdot[H2O]/d[HO2]
  J[222] -= dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[705] += dqdT; // dwdot[O2]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT

  // reaction 16: HO2 + OH <=> H2O + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[8];
  k_f = 450000000 * exp(-(5499.95687057926) * invT);
  dlnkfdT = (5499.95687057926) * invT2;
  // reverse
  phi_r = sc[3] * sc[5];
  Kc = exp(-g_RT[3] - g_RT[5] + g_RT[6] + g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[3] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q; // O2
  wdot[5] += q; // H2O
  wdot[6] -= q; // OH
  wdot[8] -= q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[5];
  J[84] += dqdci; // dwdot[O2]/d[O2]
  J[86] += dqdci; // dwdot[H2O]/d[O2]
  J[87] -= dqdci; // dwdot[OH]/d[O2]
  J[89] -= dqdci; // dwdot[HO2]/d[O2]
  // d()/d[H2O]
  dqdci = -k_r * sc[3];
  J[138] += dqdci; // dwdot[O2]/d[H2O]
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[143] -= dqdci; // dwdot[HO2]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[8];
  J[165] += dqdci; // dwdot[O2]/d[OH]
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[170] -= dqdci; // dwdot[HO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[6];
  J[219] += dqdci; // dwdot[O2]/d[HO2]
  J[221] += dqdci; // dwdot[H2O]/d[HO2]
  J[222] -= dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[705] += dqdT; // dwdot[O2]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT

  // reaction 17: 2 HO2 <=> H2O2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[8] * sc[8]);
  k_f = 100000000 * exp(-(5555.95482115001) * invT);
  dlnkfdT = (5555.95482115001) * invT2;
  // reverse
  phi_r = sc[3] * sc[7];
  Kc = exp(-g_RT[3] - g_RT[7] + 2.000000 * g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[8]) + (h_RT[3] + h_RT[7]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q;     // O2
  wdot[7] += q;     // H2O2
  wdot[8] -= 2 * q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[7];
  J[84] += dqdci;      // dwdot[O2]/d[O2]
  J[88] += dqdci;      // dwdot[H2O2]/d[O2]
  J[89] += -2 * dqdci; // dwdot[HO2]/d[O2]
  // d()/d[H2O2]
  dqdci = -k_r * sc[3];
  J[192] += dqdci;      // dwdot[O2]/d[H2O2]
  J[196] += dqdci;      // dwdot[H2O2]/d[H2O2]
  J[197] += -2 * dqdci; // dwdot[HO2]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * 2.000000 * sc[8];
  J[219] += dqdci;      // dwdot[O2]/d[HO2]
  J[223] += dqdci;      // dwdot[H2O2]/d[HO2]
  J[224] += -2 * dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[705] += dqdT;      // dwdot[O2]/dT
  J[709] += dqdT;      // dwdot[H2O2]/dT
  J[710] += -2 * dqdT; // dwdot[HO2]/dT

  // reaction 18: 2 HO2 <=> H2O2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[8] * sc[8]);
  k_f = 190000 * exp(-(-708.992024785586) * invT);
  dlnkfdT = (-708.992024785586) * invT2;
  // reverse
  phi_r = sc[3] * sc[7];
  Kc = exp(-g_RT[3] - g_RT[7] + 2.000000 * g_RT[8]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[8]) + (h_RT[3] + h_RT[7]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q;     // O2
  wdot[7] += q;     // H2O2
  wdot[8] -= 2 * q; // HO2
  // d()/d[O2]
  dqdci = -k_r * sc[7];
  J[84] += dqdci;      // dwdot[O2]/d[O2]
  J[88] += dqdci;      // dwdot[H2O2]/d[O2]
  J[89] += -2 * dqdci; // dwdot[HO2]/d[O2]
  // d()/d[H2O2]
  dqdci = -k_r * sc[3];
  J[192] += dqdci;      // dwdot[O2]/d[H2O2]
  J[196] += dqdci;      // dwdot[H2O2]/d[H2O2]
  J[197] += -2 * dqdci; // dwdot[HO2]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * 2.000000 * sc[8];
  J[219] += dqdci;      // dwdot[O2]/d[HO2]
  J[223] += dqdci;      // dwdot[H2O2]/d[HO2]
  J[224] += -2 * dqdci; // dwdot[HO2]/d[HO2]
  // d()/dT
  J[705] += dqdT;      // dwdot[O2]/dT
  J[709] += dqdT;      // dwdot[H2O2]/dT
  J[710] += -2 * dqdT; // dwdot[HO2]/dT

  // reaction 21: CO + OH <=> CO2 + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[9];
  k_f = 0.07015 * exp(2.053 * tc[0] - (-178.994168026739) * invT);
  dlnkfdT = 2.053 * invT + (-178.994168026739) * invT2;
  // reverse
  phi_r = sc[2] * sc[10];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[9] - g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[2] + h_RT[10]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[9] -= q;  // CO
  wdot[10] += q; // CO2
  // d()/d[H]
  dqdci = -k_r * sc[10];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[63] -= dqdci; // dwdot[CO]/d[H]
  J[64] += dqdci; // dwdot[CO2]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[9];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[171] -= dqdci; // dwdot[CO]/d[OH]
  J[172] += dqdci; // dwdot[CO2]/d[OH]
  // d()/d[CO]
  dqdci = +k_f * sc[6];
  J[245] += dqdci; // dwdot[H]/d[CO]
  J[249] -= dqdci; // dwdot[OH]/d[CO]
  J[252] -= dqdci; // dwdot[CO]/d[CO]
  J[253] += dqdci; // dwdot[CO2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[2];
  J[272] += dqdci; // dwdot[H]/d[CO2]
  J[276] -= dqdci; // dwdot[OH]/d[CO2]
  J[279] -= dqdci; // dwdot[CO]/d[CO2]
  J[280] += dqdci; // dwdot[CO2]/d[CO2]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[711] -= dqdT; // dwdot[CO]/dT
  J[712] += dqdT; // dwdot[CO2]/dT

  // reaction 22: CO + OH <=> CO2 + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[9];
  k_f = 5757000 * exp(-0.664 * tc[0] - (166.967289714006) * invT);
  dlnkfdT = -0.664 * invT + (166.967289714006) * invT2;
  // reverse
  phi_r = sc[2] * sc[10];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[9] - g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[2] + h_RT[10]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[9] -= q;  // CO
  wdot[10] += q; // CO2
  // d()/d[H]
  dqdci = -k_r * sc[10];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[63] -= dqdci; // dwdot[CO]/d[H]
  J[64] += dqdci; // dwdot[CO2]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[9];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[171] -= dqdci; // dwdot[CO]/d[OH]
  J[172] += dqdci; // dwdot[CO2]/d[OH]
  // d()/d[CO]
  dqdci = +k_f * sc[6];
  J[245] += dqdci; // dwdot[H]/d[CO]
  J[249] -= dqdci; // dwdot[OH]/d[CO]
  J[252] -= dqdci; // dwdot[CO]/d[CO]
  J[253] += dqdci; // dwdot[CO2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[2];
  J[272] += dqdci; // dwdot[H]/d[CO2]
  J[276] -= dqdci; // dwdot[OH]/d[CO2]
  J[279] -= dqdci; // dwdot[CO]/d[CO2]
  J[280] += dqdci; // dwdot[CO2]/d[CO2]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[711] -= dqdT; // dwdot[CO]/dT
  J[712] += dqdT; // dwdot[CO2]/dT

  // reaction 23: CO + HO2 <=> CO2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[9];
  k_f = 0.157 * exp(2.18 * tc[0] - (9027.70698453667) * invT);
  dlnkfdT = 2.18 * invT + (9027.70698453667) * invT2;
  // reverse
  phi_r = sc[6] * sc[10];
  Kc = exp(-g_RT[6] + g_RT[8] + g_RT[9] - g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[9]) + (h_RT[6] + h_RT[10]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[6] += q;  // OH
  wdot[8] -= q;  // HO2
  wdot[9] -= q;  // CO
  wdot[10] += q; // CO2
  // d()/d[OH]
  dqdci = -k_r * sc[10];
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[170] -= dqdci; // dwdot[HO2]/d[OH]
  J[171] -= dqdci; // dwdot[CO]/d[OH]
  J[172] += dqdci; // dwdot[CO2]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[9];
  J[222] += dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[225] -= dqdci; // dwdot[CO]/d[HO2]
  J[226] += dqdci; // dwdot[CO2]/d[HO2]
  // d()/d[CO]
  dqdci = +k_f * sc[8];
  J[249] += dqdci; // dwdot[OH]/d[CO]
  J[251] -= dqdci; // dwdot[HO2]/d[CO]
  J[252] -= dqdci; // dwdot[CO]/d[CO]
  J[253] += dqdci; // dwdot[CO2]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[6];
  J[276] += dqdci; // dwdot[OH]/d[CO2]
  J[278] -= dqdci; // dwdot[HO2]/d[CO2]
  J[279] -= dqdci; // dwdot[CO]/d[CO2]
  J[280] += dqdci; // dwdot[CO2]/d[CO2]
  // d()/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[711] -= dqdT; // dwdot[CO]/dT
  J[712] += dqdT; // dwdot[CO2]/dT

  // reaction 25: CH4 + H <=> CH3 + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[11];
  k_f = 0.614 * exp(2.5 * tc[0] - (4824.33817506985) * invT);
  dlnkfdT = 2.5 * invT + (4824.33817506985) * invT2;
  // reverse
  phi_r = sc[1] * sc[12];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[11] - g_RT[12]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[11]) + (h_RT[1] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[11] -= q; // CH4
  wdot[12] += q; // CH3
  // d()/d[H2]
  dqdci = -k_r * sc[12];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[38] -= dqdci; // dwdot[CH4]/d[H2]
  J[39] += dqdci; // dwdot[CH3]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[11];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[65] -= dqdci; // dwdot[CH4]/d[H]
  J[66] += dqdci; // dwdot[CH3]/d[H]
  // d()/d[CH4]
  dqdci = +k_f * sc[2];
  J[298] += dqdci; // dwdot[H2]/d[CH4]
  J[299] -= dqdci; // dwdot[H]/d[CH4]
  J[308] -= dqdci; // dwdot[CH4]/d[CH4]
  J[309] += dqdci; // dwdot[CH3]/d[CH4]
  // d()/d[CH3]
  dqdci = -k_r * sc[1];
  J[325] += dqdci; // dwdot[H2]/d[CH3]
  J[326] -= dqdci; // dwdot[H]/d[CH3]
  J[335] -= dqdci; // dwdot[CH4]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[713] -= dqdT; // dwdot[CH4]/dT
  J[714] += dqdT; // dwdot[CH3]/dT

  // reaction 26: CH4 + O <=> CH3 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[11];
  k_f = 1020 * exp(1.5 * tc[0] - (4327.66332592059) * invT);
  dlnkfdT = 1.5 * invT + (4327.66332592059) * invT2;
  // reverse
  phi_r = sc[6] * sc[12];
  Kc = exp(g_RT[4] - g_RT[6] + g_RT[11] - g_RT[12]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[6] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[6] += q;  // OH
  wdot[11] -= q; // CH4
  wdot[12] += q; // CH3
  // d()/d[O]
  dqdci = +k_f * sc[11];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[119] -= dqdci; // dwdot[CH4]/d[O]
  J[120] += dqdci; // dwdot[CH3]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[12];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[173] -= dqdci; // dwdot[CH4]/d[OH]
  J[174] += dqdci; // dwdot[CH3]/d[OH]
  // d()/d[CH4]
  dqdci = +k_f * sc[4];
  J[301] -= dqdci; // dwdot[O]/d[CH4]
  J[303] += dqdci; // dwdot[OH]/d[CH4]
  J[308] -= dqdci; // dwdot[CH4]/d[CH4]
  J[309] += dqdci; // dwdot[CH3]/d[CH4]
  // d()/d[CH3]
  dqdci = -k_r * sc[6];
  J[328] -= dqdci; // dwdot[O]/d[CH3]
  J[330] += dqdci; // dwdot[OH]/d[CH3]
  J[335] -= dqdci; // dwdot[CH4]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[713] -= dqdT; // dwdot[CH4]/dT
  J[714] += dqdT; // dwdot[CH3]/dT

  // reaction 27: CH4 + OH <=> CH3 + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[11];
  k_f = 0.0583 * exp(2.6 * tc[0] - (1102.04449811234) * invT);
  dlnkfdT = 2.6 * invT + (1102.04449811234) * invT2;
  // reverse
  phi_r = sc[5] * sc[12];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[11] - g_RT[12]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[11]) + (h_RT[5] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[11] -= q; // CH4
  wdot[12] += q; // CH3
  // d()/d[H2O]
  dqdci = -k_r * sc[12];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[146] -= dqdci; // dwdot[CH4]/d[H2O]
  J[147] += dqdci; // dwdot[CH3]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[11];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[173] -= dqdci; // dwdot[CH4]/d[OH]
  J[174] += dqdci; // dwdot[CH3]/d[OH]
  // d()/d[CH4]
  dqdci = +k_f * sc[6];
  J[302] += dqdci; // dwdot[H2O]/d[CH4]
  J[303] -= dqdci; // dwdot[OH]/d[CH4]
  J[308] -= dqdci; // dwdot[CH4]/d[CH4]
  J[309] += dqdci; // dwdot[CH3]/d[CH4]
  // d()/d[CH3]
  dqdci = -k_r * sc[5];
  J[329] += dqdci; // dwdot[H2O]/d[CH3]
  J[330] -= dqdci; // dwdot[OH]/d[CH3]
  J[335] -= dqdci; // dwdot[CH4]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[713] -= dqdT; // dwdot[CH4]/dT
  J[714] += dqdT; // dwdot[CH3]/dT

  // reaction 28: CH4 + HO2 <=> CH3 + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[11];
  k_f = 1.13e-05 * exp(3.74 * tc[0] - (10572.5821485572) * invT);
  dlnkfdT = 3.74 * invT + (10572.5821485572) * invT2;
  // reverse
  phi_r = sc[7] * sc[12];
  Kc = exp(-g_RT[7] + g_RT[8] + g_RT[11] - g_RT[12]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[11]) + (h_RT[7] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[11] -= q; // CH4
  wdot[12] += q; // CH3
  // d()/d[H2O2]
  dqdci = -k_r * sc[12];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[200] -= dqdci; // dwdot[CH4]/d[H2O2]
  J[201] += dqdci; // dwdot[CH3]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[11];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[227] -= dqdci; // dwdot[CH4]/d[HO2]
  J[228] += dqdci; // dwdot[CH3]/d[HO2]
  // d()/d[CH4]
  dqdci = +k_f * sc[8];
  J[304] += dqdci; // dwdot[H2O2]/d[CH4]
  J[305] -= dqdci; // dwdot[HO2]/d[CH4]
  J[308] -= dqdci; // dwdot[CH4]/d[CH4]
  J[309] += dqdci; // dwdot[CH3]/d[CH4]
  // d()/d[CH3]
  dqdci = -k_r * sc[7];
  J[331] += dqdci; // dwdot[H2O2]/d[CH3]
  J[332] -= dqdci; // dwdot[HO2]/d[CH3]
  J[335] -= dqdci; // dwdot[CH4]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[713] -= dqdT; // dwdot[CH4]/dT
  J[714] += dqdT; // dwdot[CH3]/dT

  // reaction 29: CH3 + HO2 <=> CH4 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[12];
  k_f = 0.116 * exp(2.23 * tc[0] - (-1520.72076406186) * invT);
  dlnkfdT = 2.23 * invT + (-1520.72076406186) * invT2;
  // reverse
  phi_r = sc[3] * sc[11];
  Kc = exp(-g_RT[3] + g_RT[8] - g_RT[11] + g_RT[12]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[3] + h_RT[11]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] += q;  // O2
  wdot[8] -= q;  // HO2
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  // d()/d[O2]
  dqdci = -k_r * sc[11];
  J[84] += dqdci; // dwdot[O2]/d[O2]
  J[89] -= dqdci; // dwdot[HO2]/d[O2]
  J[92] += dqdci; // dwdot[CH4]/d[O2]
  J[93] -= dqdci; // dwdot[CH3]/d[O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[12];
  J[219] += dqdci; // dwdot[O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[227] += dqdci; // dwdot[CH4]/d[HO2]
  J[228] -= dqdci; // dwdot[CH3]/d[HO2]
  // d()/d[CH4]
  dqdci = -k_r * sc[3];
  J[300] += dqdci; // dwdot[O2]/d[CH4]
  J[305] -= dqdci; // dwdot[HO2]/d[CH4]
  J[308] += dqdci; // dwdot[CH4]/d[CH4]
  J[309] -= dqdci; // dwdot[CH3]/d[CH4]
  // d()/d[CH3]
  dqdci = +k_f * sc[8];
  J[327] += dqdci; // dwdot[O2]/d[CH3]
  J[332] -= dqdci; // dwdot[HO2]/d[CH3]
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  // d()/dT
  J[705] += dqdT; // dwdot[O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT

  // reaction 30: CH2 + CH4 <=> 2 CH3
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[11] * sc[13];
  k_f = 2.46 * exp(2 * tc[0] - (4161.60182620503) * invT);
  dlnkfdT = 2 * invT + (4161.60182620503) * invT2;
  // reverse
  phi_r = (sc[12] * sc[12]);
  Kc = exp(g_RT[11] - 2.000000 * g_RT[12] + g_RT[13]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[11] + h_RT[13]) + (2.000000 * h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[11] -= q;     // CH4
  wdot[12] += 2 * q; // CH3
  wdot[13] -= q;     // CH2
  // d()/d[CH4]
  dqdci = +k_f * sc[13];
  J[308] -= dqdci;     // dwdot[CH4]/d[CH4]
  J[309] += 2 * dqdci; // dwdot[CH3]/d[CH4]
  J[310] -= dqdci;     // dwdot[CH2]/d[CH4]
  // d()/d[CH3]
  dqdci = -k_r * 2.000000 * sc[12];
  J[335] -= dqdci;     // dwdot[CH4]/d[CH3]
  J[336] += 2 * dqdci; // dwdot[CH3]/d[CH3]
  J[337] -= dqdci;     // dwdot[CH2]/d[CH3]
  // d()/d[CH2]
  dqdci = +k_f * sc[11];
  J[362] -= dqdci;     // dwdot[CH4]/d[CH2]
  J[363] += 2 * dqdci; // dwdot[CH3]/d[CH2]
  J[364] -= dqdci;     // dwdot[CH2]/d[CH2]
  // d()/dT
  J[713] -= dqdT;     // dwdot[CH4]/dT
  J[714] += 2 * dqdT; // dwdot[CH3]/dT
  J[715] -= dqdT;     // dwdot[CH2]/dT

  // reaction 31: CH2(S) + N2 <=> CH2 + N2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[14];
  k_f = 15000000 * exp(-(301.929999482832) * invT);
  dlnkfdT = (301.929999482832) * invT2;
  // reverse
  phi_r = sc[0] * sc[13];
  Kc = exp(g_RT[0] - g_RT[0] - g_RT[13] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[14]) + (h_RT[0] + h_RT[13]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[13] += q; // CH2
  wdot[14] -= q; // CH2(S)
  // d()/d[N2]
  dqdci = +k_f * sc[14] - k_r * sc[13];
  J[13] += dqdci; // dwdot[CH2]/d[N2]
  J[14] -= dqdci; // dwdot[CH2(S)]/d[N2]
  // d()/d[CH2]
  dqdci = -k_r * sc[0];
  J[364] += dqdci; // dwdot[CH2]/d[CH2]
  J[365] -= dqdci; // dwdot[CH2(S)]/d[CH2]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[0];
  J[391] += dqdci; // dwdot[CH2]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[715] += dqdT; // dwdot[CH2]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 32: CH2(S) + H2O <=> CH2 + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[14];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[5] * sc[13];
  Kc = exp(g_RT[5] - g_RT[5] - g_RT[13] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[14]) + (h_RT[5] + h_RT[13]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[13] += q; // CH2
  wdot[14] -= q; // CH2(S)
  // d()/d[H2O]
  dqdci = +k_f * sc[14] - k_r * sc[13];
  J[148] += dqdci; // dwdot[CH2]/d[H2O]
  J[149] -= dqdci; // dwdot[CH2(S)]/d[H2O]
  // d()/d[CH2]
  dqdci = -k_r * sc[5];
  J[364] += dqdci; // dwdot[CH2]/d[CH2]
  J[365] -= dqdci; // dwdot[CH2(S)]/d[CH2]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[5];
  J[391] += dqdci; // dwdot[CH2]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[715] += dqdT; // dwdot[CH2]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 33: CH2(S) + CO <=> CH2 + CO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[9] * sc[14];
  k_f = 9000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[9] * sc[13];
  Kc = exp(g_RT[9] - g_RT[9] - g_RT[13] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[9] + h_RT[14]) + (h_RT[9] + h_RT[13]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[13] += q; // CH2
  wdot[14] -= q; // CH2(S)
  // d()/d[CO]
  dqdci = +k_f * sc[14] - k_r * sc[13];
  J[256] += dqdci; // dwdot[CH2]/d[CO]
  J[257] -= dqdci; // dwdot[CH2(S)]/d[CO]
  // d()/d[CH2]
  dqdci = -k_r * sc[9];
  J[364] += dqdci; // dwdot[CH2]/d[CH2]
  J[365] -= dqdci; // dwdot[CH2(S)]/d[CH2]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[9];
  J[391] += dqdci; // dwdot[CH2]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[715] += dqdT; // dwdot[CH2]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 34: CH2(S) + CO2 <=> CH2 + CO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[10] * sc[14];
  k_f = 7000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[10] * sc[13];
  Kc = exp(g_RT[10] - g_RT[10] - g_RT[13] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[10] + h_RT[14]) + (h_RT[10] + h_RT[13]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[13] += q; // CH2
  wdot[14] -= q; // CH2(S)
  // d()/d[CO2]
  dqdci = +k_f * sc[14] - k_r * sc[13];
  J[283] += dqdci; // dwdot[CH2]/d[CO2]
  J[284] -= dqdci; // dwdot[CH2(S)]/d[CO2]
  // d()/d[CH2]
  dqdci = -k_r * sc[10];
  J[364] += dqdci; // dwdot[CH2]/d[CH2]
  J[365] -= dqdci; // dwdot[CH2(S)]/d[CH2]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[10];
  J[391] += dqdci; // dwdot[CH2]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[715] += dqdT; // dwdot[CH2]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 35: CH2(S) + O2 => CO + H + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[14];
  k_f = 28000000;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[2] += q;  // H
  wdot[3] -= q;  // O2
  wdot[6] += q;  // OH
  wdot[9] += q;  // CO
  wdot[14] -= q; // CH2(S)
  // d()/d[O2]
  dqdci = +k_f * sc[14];
  J[83] += dqdci; // dwdot[H]/d[O2]
  J[84] -= dqdci; // dwdot[O2]/d[O2]
  J[87] += dqdci; // dwdot[OH]/d[O2]
  J[90] += dqdci; // dwdot[CO]/d[O2]
  J[95] -= dqdci; // dwdot[CH2(S)]/d[O2]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[3];
  J[380] += dqdci; // dwdot[H]/d[CH2(S)]
  J[381] -= dqdci; // dwdot[O2]/d[CH2(S)]
  J[384] += dqdci; // dwdot[OH]/d[CH2(S)]
  J[387] += dqdci; // dwdot[CO]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 36: CH2(S) + O2 <=> CO + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[14];
  k_f = 12000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[5] * sc[9];
  Kc = exp(g_RT[3] - g_RT[5] - g_RT[9] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[14]) + (h_RT[5] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[5] += q;  // H2O
  wdot[9] += q;  // CO
  wdot[14] -= q; // CH2(S)
  // d()/d[O2]
  dqdci = +k_f * sc[14];
  J[84] -= dqdci; // dwdot[O2]/d[O2]
  J[86] += dqdci; // dwdot[H2O]/d[O2]
  J[90] += dqdci; // dwdot[CO]/d[O2]
  J[95] -= dqdci; // dwdot[CH2(S)]/d[O2]
  // d()/d[H2O]
  dqdci = -k_r * sc[9];
  J[138] -= dqdci; // dwdot[O2]/d[H2O]
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[144] += dqdci; // dwdot[CO]/d[H2O]
  J[149] -= dqdci; // dwdot[CH2(S)]/d[H2O]
  // d()/d[CO]
  dqdci = -k_r * sc[5];
  J[246] -= dqdci; // dwdot[O2]/d[CO]
  J[248] += dqdci; // dwdot[H2O]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[257] -= dqdci; // dwdot[CH2(S)]/d[CO]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[3];
  J[381] -= dqdci; // dwdot[O2]/d[CH2(S)]
  J[383] += dqdci; // dwdot[H2O]/d[CH2(S)]
  J[387] += dqdci; // dwdot[CO]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 37: CH2(S) + O <=> CO + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[14];
  k_f = 15000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[9];
  Kc = exp(-g_RT[1] + g_RT[4] - g_RT[9] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[14]) + (h_RT[1] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[4] -= q;  // O
  wdot[9] += q;  // CO
  wdot[14] -= q; // CH2(S)
  // d()/d[H2]
  dqdci = -k_r * sc[9];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[31] -= dqdci; // dwdot[O]/d[H2]
  J[36] += dqdci; // dwdot[CO]/d[H2]
  J[41] -= dqdci; // dwdot[CH2(S)]/d[H2]
  // d()/d[O]
  dqdci = +k_f * sc[14];
  J[109] += dqdci; // dwdot[H2]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[117] += dqdci; // dwdot[CO]/d[O]
  J[122] -= dqdci; // dwdot[CH2(S)]/d[O]
  // d()/d[CO]
  dqdci = -k_r * sc[1];
  J[244] += dqdci; // dwdot[H2]/d[CO]
  J[247] -= dqdci; // dwdot[O]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[257] -= dqdci; // dwdot[CH2(S)]/d[CO]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[4];
  J[379] += dqdci; // dwdot[H2]/d[CH2(S)]
  J[382] -= dqdci; // dwdot[O]/d[CH2(S)]
  J[387] += dqdci; // dwdot[CO]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 38: CH2(S) + O <=> H + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[14];
  k_f = 15000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[21];
  Kc = exp(-g_RT[2] + g_RT[4] + g_RT[14] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[14]) + (h_RT[2] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[4] -= q;  // O
  wdot[14] -= q; // CH2(S)
  wdot[21] += q; // HCO
  // d()/d[H]
  dqdci = -k_r * sc[21];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[58] -= dqdci; // dwdot[O]/d[H]
  J[68] -= dqdci; // dwdot[CH2(S)]/d[H]
  J[75] += dqdci; // dwdot[HCO]/d[H]
  // d()/d[O]
  dqdci = +k_f * sc[14];
  J[110] += dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[122] -= dqdci; // dwdot[CH2(S)]/d[O]
  J[129] += dqdci; // dwdot[HCO]/d[O]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[4];
  J[380] += dqdci; // dwdot[H]/d[CH2(S)]
  J[382] -= dqdci; // dwdot[O]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  J[399] += dqdci; // dwdot[HCO]/d[CH2(S)]
  // d()/d[HCO]
  dqdci = -k_r * sc[2];
  J[569] += dqdci; // dwdot[H]/d[HCO]
  J[571] -= dqdci; // dwdot[O]/d[HCO]
  J[581] -= dqdci; // dwdot[CH2(S)]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 39: CH2(S) + H2 <=> CH3 + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[14];
  k_f = 70000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[12];
  Kc = exp(g_RT[1] - g_RT[2] - g_RT[12] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[2] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q;  // H2
  wdot[2] += q;  // H
  wdot[12] += q; // CH3
  wdot[14] -= q; // CH2(S)
  // d()/d[H2]
  dqdci = +k_f * sc[14];
  J[28] -= dqdci; // dwdot[H2]/d[H2]
  J[29] += dqdci; // dwdot[H]/d[H2]
  J[39] += dqdci; // dwdot[CH3]/d[H2]
  J[41] -= dqdci; // dwdot[CH2(S)]/d[H2]
  // d()/d[H]
  dqdci = -k_r * sc[12];
  J[55] -= dqdci; // dwdot[H2]/d[H]
  J[56] += dqdci; // dwdot[H]/d[H]
  J[66] += dqdci; // dwdot[CH3]/d[H]
  J[68] -= dqdci; // dwdot[CH2(S)]/d[H]
  // d()/d[CH3]
  dqdci = -k_r * sc[2];
  J[325] -= dqdci; // dwdot[H2]/d[CH3]
  J[326] += dqdci; // dwdot[H]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  J[338] -= dqdci; // dwdot[CH2(S)]/d[CH3]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[1];
  J[379] -= dqdci; // dwdot[H2]/d[CH2(S)]
  J[380] += dqdci; // dwdot[H]/d[CH2(S)]
  J[390] += dqdci; // dwdot[CH3]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[703] -= dqdT; // dwdot[H2]/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[714] += dqdT; // dwdot[CH3]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT

  // reaction 40: CH2(S) + H <=> CH + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[14];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[16];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[14] - g_RT[16]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[14]) + (h_RT[1] + h_RT[16]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[14] -= q; // CH2(S)
  wdot[16] += q; // CH
  // d()/d[H2]
  dqdci = -k_r * sc[16];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[41] -= dqdci; // dwdot[CH2(S)]/d[H2]
  J[43] += dqdci; // dwdot[CH]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[14];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[68] -= dqdci; // dwdot[CH2(S)]/d[H]
  J[70] += dqdci; // dwdot[CH]/d[H]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[2];
  J[379] += dqdci; // dwdot[H2]/d[CH2(S)]
  J[380] -= dqdci; // dwdot[H]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  J[394] += dqdci; // dwdot[CH]/d[CH2(S)]
  // d()/d[CH]
  dqdci = -k_r * sc[1];
  J[433] += dqdci; // dwdot[H2]/d[CH]
  J[434] -= dqdci; // dwdot[H]/d[CH]
  J[446] -= dqdci; // dwdot[CH2(S)]/d[CH]
  J[448] += dqdci; // dwdot[CH]/d[CH]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT
  J[718] += dqdT; // dwdot[CH]/dT

  // reaction 41: CH2(S) + OH <=> CH2O + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[14];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[20];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[14] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[14]) + (h_RT[2] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[14] -= q; // CH2(S)
  wdot[20] += q; // CH2O
  // d()/d[H]
  dqdci = -k_r * sc[20];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[68] -= dqdci; // dwdot[CH2(S)]/d[H]
  J[74] += dqdci; // dwdot[CH2O]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[14];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[176] -= dqdci; // dwdot[CH2(S)]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[6];
  J[380] += dqdci; // dwdot[H]/d[CH2(S)]
  J[384] -= dqdci; // dwdot[OH]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  J[398] += dqdci; // dwdot[CH2O]/d[CH2(S)]
  // d()/d[CH2O]
  dqdci = -k_r * sc[2];
  J[542] += dqdci; // dwdot[H]/d[CH2O]
  J[546] -= dqdci; // dwdot[OH]/d[CH2O]
  J[554] -= dqdci; // dwdot[CH2(S)]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 42: CH2(S) + CO2 <=> CH2O + CO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[10] * sc[14];
  k_f = 14000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[9] * sc[20];
  Kc = exp(-g_RT[9] + g_RT[10] + g_RT[14] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[10] + h_RT[14]) + (h_RT[9] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] += q;  // CO
  wdot[10] -= q; // CO2
  wdot[14] -= q; // CH2(S)
  wdot[20] += q; // CH2O
  // d()/d[CO]
  dqdci = -k_r * sc[20];
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[253] -= dqdci; // dwdot[CO2]/d[CO]
  J[257] -= dqdci; // dwdot[CH2(S)]/d[CO]
  J[263] += dqdci; // dwdot[CH2O]/d[CO]
  // d()/d[CO2]
  dqdci = +k_f * sc[14];
  J[279] += dqdci; // dwdot[CO]/d[CO2]
  J[280] -= dqdci; // dwdot[CO2]/d[CO2]
  J[284] -= dqdci; // dwdot[CH2(S)]/d[CO2]
  J[290] += dqdci; // dwdot[CH2O]/d[CO2]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[10];
  J[387] += dqdci; // dwdot[CO]/d[CH2(S)]
  J[388] -= dqdci; // dwdot[CO2]/d[CH2(S)]
  J[392] -= dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  J[398] += dqdci; // dwdot[CH2O]/d[CH2(S)]
  // d()/d[CH2O]
  dqdci = -k_r * sc[9];
  J[549] += dqdci; // dwdot[CO]/d[CH2O]
  J[550] -= dqdci; // dwdot[CO2]/d[CH2O]
  J[554] -= dqdci; // dwdot[CH2(S)]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[712] -= dqdT; // dwdot[CO2]/dT
  J[716] -= dqdT; // dwdot[CH2(S)]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 44: CH2 + O2 <=> HCO + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[13];
  k_f = 10600000 * exp(-(754.82499870708) * invT);
  dlnkfdT = (754.82499870708) * invT2;
  // reverse
  phi_r = sc[6] * sc[21];
  Kc = exp(g_RT[3] - g_RT[6] + g_RT[13] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[13]) + (h_RT[6] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[6] += q;  // OH
  wdot[13] -= q; // CH2
  wdot[21] += q; // HCO
  // d()/d[O2]
  dqdci = +k_f * sc[13];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[87] += dqdci;  // dwdot[OH]/d[O2]
  J[94] -= dqdci;  // dwdot[CH2]/d[O2]
  J[102] += dqdci; // dwdot[HCO]/d[O2]
  // d()/d[OH]
  dqdci = -k_r * sc[21];
  J[165] -= dqdci; // dwdot[O2]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[175] -= dqdci; // dwdot[CH2]/d[OH]
  J[183] += dqdci; // dwdot[HCO]/d[OH]
  // d()/d[CH2]
  dqdci = +k_f * sc[3];
  J[354] -= dqdci; // dwdot[O2]/d[CH2]
  J[357] += dqdci; // dwdot[OH]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  J[372] += dqdci; // dwdot[HCO]/d[CH2]
  // d()/d[HCO]
  dqdci = -k_r * sc[6];
  J[570] -= dqdci; // dwdot[O2]/d[HCO]
  J[573] += dqdci; // dwdot[OH]/d[HCO]
  J[580] -= dqdci; // dwdot[CH2]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 45: CH2 + O2 => CO2 + 2 H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[13];
  k_f = 2640000 * exp(-(754.82499870708) * invT);
  dlnkfdT = (754.82499870708) * invT2;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[2] += 2 * q; // H
  wdot[3] -= q;     // O2
  wdot[10] += q;    // CO2
  wdot[13] -= q;    // CH2
  // d()/d[O2]
  dqdci = +k_f * sc[13];
  J[83] += 2 * dqdci; // dwdot[H]/d[O2]
  J[84] -= dqdci;     // dwdot[O2]/d[O2]
  J[91] += dqdci;     // dwdot[CO2]/d[O2]
  J[94] -= dqdci;     // dwdot[CH2]/d[O2]
  // d()/d[CH2]
  dqdci = +k_f * sc[3];
  J[353] += 2 * dqdci; // dwdot[H]/d[CH2]
  J[354] -= dqdci;     // dwdot[O2]/d[CH2]
  J[361] += dqdci;     // dwdot[CO2]/d[CH2]
  J[364] -= dqdci;     // dwdot[CH2]/d[CH2]
  // d()/dT
  J[704] += 2 * dqdT; // dwdot[H]/dT
  J[705] -= dqdT;     // dwdot[O2]/dT
  J[712] += dqdT;     // dwdot[CO2]/dT
  J[715] -= dqdT;     // dwdot[CH2]/dT

  // reaction 46: CH2 + O => CO + 2 H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[13];
  k_f = 50000000;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[2] += 2 * q; // H
  wdot[4] -= q;     // O
  wdot[9] += q;     // CO
  wdot[13] -= q;    // CH2
  // d()/d[O]
  dqdci = +k_f * sc[13];
  J[110] += 2 * dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci;     // dwdot[O]/d[O]
  J[117] += dqdci;     // dwdot[CO]/d[O]
  J[121] -= dqdci;     // dwdot[CH2]/d[O]
  // d()/d[CH2]
  dqdci = +k_f * sc[4];
  J[353] += 2 * dqdci; // dwdot[H]/d[CH2]
  J[355] -= dqdci;     // dwdot[O]/d[CH2]
  J[360] += dqdci;     // dwdot[CO]/d[CH2]
  J[364] -= dqdci;     // dwdot[CH2]/d[CH2]
  // d()/dT
  J[704] += 2 * dqdT; // dwdot[H]/dT
  J[706] -= dqdT;     // dwdot[O]/dT
  J[711] += dqdT;     // dwdot[CO]/dT
  J[715] -= dqdT;     // dwdot[CH2]/dT

  // reaction 47: CH2 + H <=> CH + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[13];
  k_f = 283000000 * exp(-(369.864249366469) * invT);
  dlnkfdT = (369.864249366469) * invT2;
  // reverse
  phi_r = sc[1] * sc[16];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[13] - g_RT[16]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[1] + h_RT[16]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[13] -= q; // CH2
  wdot[16] += q; // CH
  // d()/d[H2]
  dqdci = -k_r * sc[16];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[40] -= dqdci; // dwdot[CH2]/d[H2]
  J[43] += dqdci; // dwdot[CH]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[13];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[67] -= dqdci; // dwdot[CH2]/d[H]
  J[70] += dqdci; // dwdot[CH]/d[H]
  // d()/d[CH2]
  dqdci = +k_f * sc[2];
  J[352] += dqdci; // dwdot[H2]/d[CH2]
  J[353] -= dqdci; // dwdot[H]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  J[367] += dqdci; // dwdot[CH]/d[CH2]
  // d()/d[CH]
  dqdci = -k_r * sc[1];
  J[433] += dqdci; // dwdot[H2]/d[CH]
  J[434] -= dqdci; // dwdot[H]/d[CH]
  J[445] -= dqdci; // dwdot[CH2]/d[CH]
  J[448] += dqdci; // dwdot[CH]/d[CH]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT
  J[718] += dqdT; // dwdot[CH]/dT

  // reaction 48: CH2 + OH <=> CH + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[13];
  k_f = 11.3 * exp(2 * tc[0] - (1509.64999741416) * invT);
  dlnkfdT = 2 * invT + (1509.64999741416) * invT2;
  // reverse
  phi_r = sc[5] * sc[16];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[13] - g_RT[16]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[13]) + (h_RT[5] + h_RT[16]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[13] -= q; // CH2
  wdot[16] += q; // CH
  // d()/d[H2O]
  dqdci = -k_r * sc[16];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[148] -= dqdci; // dwdot[CH2]/d[H2O]
  J[151] += dqdci; // dwdot[CH]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[13];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[175] -= dqdci; // dwdot[CH2]/d[OH]
  J[178] += dqdci; // dwdot[CH]/d[OH]
  // d()/d[CH2]
  dqdci = +k_f * sc[6];
  J[356] += dqdci; // dwdot[H2O]/d[CH2]
  J[357] -= dqdci; // dwdot[OH]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  J[367] += dqdci; // dwdot[CH]/d[CH2]
  // d()/d[CH]
  dqdci = -k_r * sc[5];
  J[437] += dqdci; // dwdot[H2O]/d[CH]
  J[438] -= dqdci; // dwdot[OH]/d[CH]
  J[445] -= dqdci; // dwdot[CH2]/d[CH]
  J[448] += dqdci; // dwdot[CH]/d[CH]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT
  J[718] += dqdT; // dwdot[CH]/dT

  // reaction 49: CH + O2 <=> HCO + O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[16];
  k_f = 33000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[4] * sc[21];
  Kc = exp(g_RT[3] - g_RT[4] + g_RT[16] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[16]) + (h_RT[4] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[4] += q;  // O
  wdot[16] -= q; // CH
  wdot[21] += q; // HCO
  // d()/d[O2]
  dqdci = +k_f * sc[16];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[85] += dqdci;  // dwdot[O]/d[O2]
  J[97] -= dqdci;  // dwdot[CH]/d[O2]
  J[102] += dqdci; // dwdot[HCO]/d[O2]
  // d()/d[O]
  dqdci = -k_r * sc[21];
  J[111] -= dqdci; // dwdot[O2]/d[O]
  J[112] += dqdci; // dwdot[O]/d[O]
  J[124] -= dqdci; // dwdot[CH]/d[O]
  J[129] += dqdci; // dwdot[HCO]/d[O]
  // d()/d[CH]
  dqdci = +k_f * sc[3];
  J[435] -= dqdci; // dwdot[O2]/d[CH]
  J[436] += dqdci; // dwdot[O]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  J[453] += dqdci; // dwdot[HCO]/d[CH]
  // d()/d[HCO]
  dqdci = -k_r * sc[4];
  J[570] -= dqdci; // dwdot[O2]/d[HCO]
  J[571] += dqdci; // dwdot[O]/d[HCO]
  J[583] -= dqdci; // dwdot[CH]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[706] += dqdT; // dwdot[O]/dT
  J[718] -= dqdT; // dwdot[CH]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 50: CH + O <=> CO + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[16];
  k_f = 57000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[9];
  Kc = exp(-g_RT[2] + g_RT[4] - g_RT[9] + g_RT[16]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[16]) + (h_RT[2] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[4] -= q;  // O
  wdot[9] += q;  // CO
  wdot[16] -= q; // CH
  // d()/d[H]
  dqdci = -k_r * sc[9];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[58] -= dqdci; // dwdot[O]/d[H]
  J[63] += dqdci; // dwdot[CO]/d[H]
  J[70] -= dqdci; // dwdot[CH]/d[H]
  // d()/d[O]
  dqdci = +k_f * sc[16];
  J[110] += dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[117] += dqdci; // dwdot[CO]/d[O]
  J[124] -= dqdci; // dwdot[CH]/d[O]
  // d()/d[CO]
  dqdci = -k_r * sc[2];
  J[245] += dqdci; // dwdot[H]/d[CO]
  J[247] -= dqdci; // dwdot[O]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[259] -= dqdci; // dwdot[CH]/d[CO]
  // d()/d[CH]
  dqdci = +k_f * sc[4];
  J[434] += dqdci; // dwdot[H]/d[CH]
  J[436] -= dqdci; // dwdot[O]/d[CH]
  J[441] += dqdci; // dwdot[CO]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[718] -= dqdT; // dwdot[CH]/dT

  // reaction 51: CH + H <=> C + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[16];
  k_f = 110000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[15];
  Kc = exp(-g_RT[1] + g_RT[2] - g_RT[15] + g_RT[16]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[16]) + (h_RT[1] + h_RT[15]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[15] += q; // C
  wdot[16] -= q; // CH
  // d()/d[H2]
  dqdci = -k_r * sc[15];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[42] += dqdci; // dwdot[C]/d[H2]
  J[43] -= dqdci; // dwdot[CH]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[16];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[69] += dqdci; // dwdot[C]/d[H]
  J[70] -= dqdci; // dwdot[CH]/d[H]
  // d()/d[C]
  dqdci = -k_r * sc[1];
  J[406] += dqdci; // dwdot[H2]/d[C]
  J[407] -= dqdci; // dwdot[H]/d[C]
  J[420] += dqdci; // dwdot[C]/d[C]
  J[421] -= dqdci; // dwdot[CH]/d[C]
  // d()/d[CH]
  dqdci = +k_f * sc[2];
  J[433] += dqdci; // dwdot[H2]/d[CH]
  J[434] -= dqdci; // dwdot[H]/d[CH]
  J[447] += dqdci; // dwdot[C]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[717] += dqdT; // dwdot[C]/dT
  J[718] -= dqdT; // dwdot[CH]/dT

  // reaction 52: CH + OH <=> H + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[16];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[21];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[16] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[16]) + (h_RT[2] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[16] -= q; // CH
  wdot[21] += q; // HCO
  // d()/d[H]
  dqdci = -k_r * sc[21];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[70] -= dqdci; // dwdot[CH]/d[H]
  J[75] += dqdci; // dwdot[HCO]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[16];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[178] -= dqdci; // dwdot[CH]/d[OH]
  J[183] += dqdci; // dwdot[HCO]/d[OH]
  // d()/d[CH]
  dqdci = +k_f * sc[6];
  J[434] += dqdci; // dwdot[H]/d[CH]
  J[438] -= dqdci; // dwdot[OH]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  J[453] += dqdci; // dwdot[HCO]/d[CH]
  // d()/d[HCO]
  dqdci = -k_r * sc[2];
  J[569] += dqdci; // dwdot[H]/d[HCO]
  J[573] -= dqdci; // dwdot[OH]/d[HCO]
  J[583] -= dqdci; // dwdot[CH]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[718] -= dqdT; // dwdot[CH]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 53: CH + H2O <=> CH2O + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[16];
  k_f = 17130000 * exp(-(-379.928582682563) * invT);
  dlnkfdT = (-379.928582682563) * invT2;
  // reverse
  phi_r = sc[2] * sc[20];
  Kc = exp(-g_RT[2] + g_RT[5] + g_RT[16] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[16]) + (h_RT[2] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[5] -= q;  // H2O
  wdot[16] -= q; // CH
  wdot[20] += q; // CH2O
  // d()/d[H]
  dqdci = -k_r * sc[20];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[59] -= dqdci; // dwdot[H2O]/d[H]
  J[70] -= dqdci; // dwdot[CH]/d[H]
  J[74] += dqdci; // dwdot[CH2O]/d[H]
  // d()/d[H2O]
  dqdci = +k_f * sc[16];
  J[137] += dqdci; // dwdot[H]/d[H2O]
  J[140] -= dqdci; // dwdot[H2O]/d[H2O]
  J[151] -= dqdci; // dwdot[CH]/d[H2O]
  J[155] += dqdci; // dwdot[CH2O]/d[H2O]
  // d()/d[CH]
  dqdci = +k_f * sc[5];
  J[434] += dqdci; // dwdot[H]/d[CH]
  J[437] -= dqdci; // dwdot[H2O]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  J[452] += dqdci; // dwdot[CH2O]/d[CH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[2];
  J[542] += dqdci; // dwdot[H]/d[CH2O]
  J[545] -= dqdci; // dwdot[H2O]/d[CH2O]
  J[556] -= dqdci; // dwdot[CH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[707] -= dqdT; // dwdot[H2O]/dT
  J[718] -= dqdT; // dwdot[CH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 54: CH + CO2 <=> CO + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[10] * sc[16];
  k_f = 1700000 * exp(-(344.703416076233) * invT);
  dlnkfdT = (344.703416076233) * invT2;
  // reverse
  phi_r = sc[9] * sc[21];
  Kc = exp(-g_RT[9] + g_RT[10] + g_RT[16] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[10] + h_RT[16]) + (h_RT[9] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] += q;  // CO
  wdot[10] -= q; // CO2
  wdot[16] -= q; // CH
  wdot[21] += q; // HCO
  // d()/d[CO]
  dqdci = -k_r * sc[21];
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[253] -= dqdci; // dwdot[CO2]/d[CO]
  J[259] -= dqdci; // dwdot[CH]/d[CO]
  J[264] += dqdci; // dwdot[HCO]/d[CO]
  // d()/d[CO2]
  dqdci = +k_f * sc[16];
  J[279] += dqdci; // dwdot[CO]/d[CO2]
  J[280] -= dqdci; // dwdot[CO2]/d[CO2]
  J[286] -= dqdci; // dwdot[CH]/d[CO2]
  J[291] += dqdci; // dwdot[HCO]/d[CO2]
  // d()/d[CH]
  dqdci = +k_f * sc[10];
  J[441] += dqdci; // dwdot[CO]/d[CH]
  J[442] -= dqdci; // dwdot[CO2]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  J[453] += dqdci; // dwdot[HCO]/d[CH]
  // d()/d[HCO]
  dqdci = -k_r * sc[9];
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[577] -= dqdci; // dwdot[CO2]/d[HCO]
  J[583] -= dqdci; // dwdot[CH]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[712] -= dqdT; // dwdot[CO2]/dT
  J[718] -= dqdT; // dwdot[CH]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 55: C + OH <=> CO + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[15];
  k_f = 50000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[9];
  Kc = exp(-g_RT[2] + g_RT[6] - g_RT[9] + g_RT[15]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[15]) + (h_RT[2] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[9] += q;  // CO
  wdot[15] -= q; // C
  // d()/d[H]
  dqdci = -k_r * sc[9];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[63] += dqdci; // dwdot[CO]/d[H]
  J[69] -= dqdci; // dwdot[C]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[15];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[171] += dqdci; // dwdot[CO]/d[OH]
  J[177] -= dqdci; // dwdot[C]/d[OH]
  // d()/d[CO]
  dqdci = -k_r * sc[2];
  J[245] += dqdci; // dwdot[H]/d[CO]
  J[249] -= dqdci; // dwdot[OH]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[258] -= dqdci; // dwdot[C]/d[CO]
  // d()/d[C]
  dqdci = +k_f * sc[6];
  J[407] += dqdci; // dwdot[H]/d[C]
  J[411] -= dqdci; // dwdot[OH]/d[C]
  J[414] += dqdci; // dwdot[CO]/d[C]
  J[420] -= dqdci; // dwdot[C]/d[C]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[717] -= dqdT; // dwdot[C]/dT

  // reaction 56: C + O2 <=> CO + O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[15];
  k_f = 50000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[4] * sc[9];
  Kc = exp(g_RT[3] - g_RT[4] - g_RT[9] + g_RT[15]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[15]) + (h_RT[4] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[4] += q;  // O
  wdot[9] += q;  // CO
  wdot[15] -= q; // C
  // d()/d[O2]
  dqdci = +k_f * sc[15];
  J[84] -= dqdci; // dwdot[O2]/d[O2]
  J[85] += dqdci; // dwdot[O]/d[O2]
  J[90] += dqdci; // dwdot[CO]/d[O2]
  J[96] -= dqdci; // dwdot[C]/d[O2]
  // d()/d[O]
  dqdci = -k_r * sc[9];
  J[111] -= dqdci; // dwdot[O2]/d[O]
  J[112] += dqdci; // dwdot[O]/d[O]
  J[117] += dqdci; // dwdot[CO]/d[O]
  J[123] -= dqdci; // dwdot[C]/d[O]
  // d()/d[CO]
  dqdci = -k_r * sc[4];
  J[246] -= dqdci; // dwdot[O2]/d[CO]
  J[247] += dqdci; // dwdot[O]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[258] -= dqdci; // dwdot[C]/d[CO]
  // d()/d[C]
  dqdci = +k_f * sc[3];
  J[408] -= dqdci; // dwdot[O2]/d[C]
  J[409] += dqdci; // dwdot[O]/d[C]
  J[414] += dqdci; // dwdot[CO]/d[C]
  J[420] -= dqdci; // dwdot[C]/d[C]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[706] += dqdT; // dwdot[O]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[717] -= dqdT; // dwdot[C]/dT

  // reaction 57: CH3 + O2 <=> CH3O + O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[12];
  k_f = 7546000 * exp(-(14251.0959755897) * invT);
  dlnkfdT = (14251.0959755897) * invT2;
  // reverse
  phi_r = sc[4] * sc[18];
  Kc = exp(g_RT[3] - g_RT[4] + g_RT[12] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[12]) + (h_RT[4] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[4] += q;  // O
  wdot[12] -= q; // CH3
  wdot[18] += q; // CH3O
  // d()/d[O2]
  dqdci = +k_f * sc[12];
  J[84] -= dqdci; // dwdot[O2]/d[O2]
  J[85] += dqdci; // dwdot[O]/d[O2]
  J[93] -= dqdci; // dwdot[CH3]/d[O2]
  J[99] += dqdci; // dwdot[CH3O]/d[O2]
  // d()/d[O]
  dqdci = -k_r * sc[18];
  J[111] -= dqdci; // dwdot[O2]/d[O]
  J[112] += dqdci; // dwdot[O]/d[O]
  J[120] -= dqdci; // dwdot[CH3]/d[O]
  J[126] += dqdci; // dwdot[CH3O]/d[O]
  // d()/d[CH3]
  dqdci = +k_f * sc[3];
  J[327] -= dqdci; // dwdot[O2]/d[CH3]
  J[328] += dqdci; // dwdot[O]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[342] += dqdci; // dwdot[CH3O]/d[CH3]
  // d()/d[CH3O]
  dqdci = -k_r * sc[4];
  J[489] -= dqdci; // dwdot[O2]/d[CH3O]
  J[490] += dqdci; // dwdot[O]/d[CH3O]
  J[498] -= dqdci; // dwdot[CH3]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[706] += dqdT; // dwdot[O]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 58: CH3 + O2 <=> CH2O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[12];
  k_f = 2.641e-06 * exp(3.283 * tc[0] - (4078.57107634725) * invT);
  dlnkfdT = 3.283 * invT + (4078.57107634725) * invT2;
  // reverse
  phi_r = sc[6] * sc[20];
  Kc = exp(g_RT[3] - g_RT[6] + g_RT[12] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[12]) + (h_RT[6] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[6] += q;  // OH
  wdot[12] -= q; // CH3
  wdot[20] += q; // CH2O
  // d()/d[O2]
  dqdci = +k_f * sc[12];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[87] += dqdci;  // dwdot[OH]/d[O2]
  J[93] -= dqdci;  // dwdot[CH3]/d[O2]
  J[101] += dqdci; // dwdot[CH2O]/d[O2]
  // d()/d[OH]
  dqdci = -k_r * sc[20];
  J[165] -= dqdci; // dwdot[O2]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH3]
  dqdci = +k_f * sc[3];
  J[327] -= dqdci; // dwdot[O2]/d[CH3]
  J[330] += dqdci; // dwdot[OH]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[344] += dqdci; // dwdot[CH2O]/d[CH3]
  // d()/d[CH2O]
  dqdci = -k_r * sc[6];
  J[543] -= dqdci; // dwdot[O2]/d[CH2O]
  J[546] += dqdci; // dwdot[OH]/d[CH2O]
  J[552] -= dqdci; // dwdot[CH3]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 59: CH3 + O <=> CH2O + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[12];
  k_f = 55400000 * exp(0.05 * tc[0] - (-68.4374665494419) * invT);
  dlnkfdT = 0.05 * invT + (-68.4374665494419) * invT2;
  // reverse
  phi_r = sc[2] * sc[20];
  Kc = exp(-g_RT[2] + g_RT[4] + g_RT[12] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[12]) + (h_RT[2] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[4] -= q;  // O
  wdot[12] -= q; // CH3
  wdot[20] += q; // CH2O
  // d()/d[H]
  dqdci = -k_r * sc[20];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[58] -= dqdci; // dwdot[O]/d[H]
  J[66] -= dqdci; // dwdot[CH3]/d[H]
  J[74] += dqdci; // dwdot[CH2O]/d[H]
  // d()/d[O]
  dqdci = +k_f * sc[12];
  J[110] += dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[120] -= dqdci; // dwdot[CH3]/d[O]
  J[128] += dqdci; // dwdot[CH2O]/d[O]
  // d()/d[CH3]
  dqdci = +k_f * sc[4];
  J[326] += dqdci; // dwdot[H]/d[CH3]
  J[328] -= dqdci; // dwdot[O]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[344] += dqdci; // dwdot[CH2O]/d[CH3]
  // d()/d[CH2O]
  dqdci = -k_r * sc[2];
  J[542] += dqdci; // dwdot[H]/d[CH2O]
  J[544] -= dqdci; // dwdot[O]/d[CH2O]
  J[552] -= dqdci; // dwdot[CH3]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 60: CH3 + OH <=> CH2(S) + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[12];
  k_f = 4.788e+17 * exp(-3.155 * tc[0] - (3524.02631063045) * invT);
  dlnkfdT = -3.155 * invT + (3524.02631063045) * invT2;
  // reverse
  phi_r = sc[5] * sc[14];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[12] - g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[5] + h_RT[14]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[12] -= q; // CH3
  wdot[14] += q; // CH2(S)
  // d()/d[H2O]
  dqdci = -k_r * sc[14];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[147] -= dqdci; // dwdot[CH3]/d[H2O]
  J[149] += dqdci; // dwdot[CH2(S)]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[12];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[176] += dqdci; // dwdot[CH2(S)]/d[OH]
  // d()/d[CH3]
  dqdci = +k_f * sc[6];
  J[329] += dqdci; // dwdot[H2O]/d[CH3]
  J[330] -= dqdci; // dwdot[OH]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[338] += dqdci; // dwdot[CH2(S)]/d[CH3]
  // d()/d[CH2(S)]
  dqdci = -k_r * sc[5];
  J[383] += dqdci; // dwdot[H2O]/d[CH2(S)]
  J[384] -= dqdci; // dwdot[OH]/d[CH2(S)]
  J[390] -= dqdci; // dwdot[CH3]/d[CH2(S)]
  J[392] += dqdci; // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[716] += dqdT; // dwdot[CH2(S)]/dT

  // reaction 61: CH3 + OH <=> CH2O + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[12];
  k_f = 5374 * exp(0.287 * tc[0] - (140.900666425322) * invT);
  dlnkfdT = 0.287 * invT + (140.900666425322) * invT2;
  // reverse
  phi_r = sc[1] * sc[20];
  Kc = exp(-g_RT[1] + g_RT[6] + g_RT[12] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[1] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[6] -= q;  // OH
  wdot[12] -= q; // CH3
  wdot[20] += q; // CH2O
  // d()/d[H2]
  dqdci = -k_r * sc[20];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[33] -= dqdci; // dwdot[OH]/d[H2]
  J[39] -= dqdci; // dwdot[CH3]/d[H2]
  J[47] += dqdci; // dwdot[CH2O]/d[H2]
  // d()/d[OH]
  dqdci = +k_f * sc[12];
  J[163] += dqdci; // dwdot[H2]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH3]
  dqdci = +k_f * sc[6];
  J[325] += dqdci; // dwdot[H2]/d[CH3]
  J[330] -= dqdci; // dwdot[OH]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[344] += dqdci; // dwdot[CH2O]/d[CH3]
  // d()/d[CH2O]
  dqdci = -k_r * sc[1];
  J[541] += dqdci; // dwdot[H2]/d[CH2O]
  J[546] -= dqdci; // dwdot[OH]/d[CH2O]
  J[552] -= dqdci; // dwdot[CH3]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 62: CH3 + OH <=> CH2OH + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[12];
  k_f = 15250000 * exp(0.134 * tc[0] - (2838.64521180442) * invT);
  dlnkfdT = 0.134 * invT + (2838.64521180442) * invT2;
  // reverse
  phi_r = sc[2] * sc[19];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[12] - g_RT[19]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[2] + h_RT[19]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[12] -= q; // CH3
  wdot[19] += q; // CH2OH
  // d()/d[H]
  dqdci = -k_r * sc[19];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[66] -= dqdci; // dwdot[CH3]/d[H]
  J[73] += dqdci; // dwdot[CH2OH]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[12];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[181] += dqdci; // dwdot[CH2OH]/d[OH]
  // d()/d[CH3]
  dqdci = +k_f * sc[6];
  J[326] += dqdci; // dwdot[H]/d[CH3]
  J[330] -= dqdci; // dwdot[OH]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[343] += dqdci; // dwdot[CH2OH]/d[CH3]
  // d()/d[CH2OH]
  dqdci = -k_r * sc[2];
  J[515] += dqdci; // dwdot[H]/d[CH2OH]
  J[519] -= dqdci; // dwdot[OH]/d[CH2OH]
  J[525] -= dqdci; // dwdot[CH3]/d[CH2OH]
  J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT

  // reaction 63: CH3 + OH <=> CH3O + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[12];
  k_f = 1798 * exp(0.965 * tc[0] - (6068.79298960492) * invT);
  dlnkfdT = 0.965 * invT + (6068.79298960492) * invT2;
  // reverse
  phi_r = sc[2] * sc[18];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[12] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[2] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[12] -= q; // CH3
  wdot[18] += q; // CH3O
  // d()/d[H]
  dqdci = -k_r * sc[18];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[66] -= dqdci; // dwdot[CH3]/d[H]
  J[72] += dqdci; // dwdot[CH3O]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[12];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[180] += dqdci; // dwdot[CH3O]/d[OH]
  // d()/d[CH3]
  dqdci = +k_f * sc[6];
  J[326] += dqdci; // dwdot[H]/d[CH3]
  J[330] -= dqdci; // dwdot[OH]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[342] += dqdci; // dwdot[CH3O]/d[CH3]
  // d()/d[CH3O]
  dqdci = -k_r * sc[2];
  J[488] += dqdci; // dwdot[H]/d[CH3O]
  J[492] -= dqdci; // dwdot[OH]/d[CH3O]
  J[498] -= dqdci; // dwdot[CH3]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 64: CH3 + OH <=> CH2 + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[12];
  k_f = 0.04293 * exp(2.568 * tc[0] - (2011.75958655411) * invT);
  dlnkfdT = 2.568 * invT + (2011.75958655411) * invT2;
  // reverse
  phi_r = sc[5] * sc[13];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[12] - g_RT[13]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[5] + h_RT[13]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[12] -= q; // CH3
  wdot[13] += q; // CH2
  // d()/d[H2O]
  dqdci = -k_r * sc[13];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[147] -= dqdci; // dwdot[CH3]/d[H2O]
  J[148] += dqdci; // dwdot[CH2]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[12];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[175] += dqdci; // dwdot[CH2]/d[OH]
  // d()/d[CH3]
  dqdci = +k_f * sc[6];
  J[329] += dqdci; // dwdot[H2O]/d[CH3]
  J[330] -= dqdci; // dwdot[OH]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[337] += dqdci; // dwdot[CH2]/d[CH3]
  // d()/d[CH2]
  dqdci = -k_r * sc[5];
  J[356] += dqdci; // dwdot[H2O]/d[CH2]
  J[357] -= dqdci; // dwdot[OH]/d[CH2]
  J[363] -= dqdci; // dwdot[CH3]/d[CH2]
  J[364] += dqdci; // dwdot[CH2]/d[CH2]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[715] += dqdT; // dwdot[CH2]/dT

  // reaction 65: CH3 + HO2 <=> CH3O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[12];
  k_f = 1000000 * exp(0.269 * tc[0] - (-345.961457740745) * invT);
  dlnkfdT = 0.269 * invT + (-345.961457740745) * invT2;
  // reverse
  phi_r = sc[6] * sc[18];
  Kc = exp(-g_RT[6] + g_RT[8] + g_RT[12] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[6] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[6] += q;  // OH
  wdot[8] -= q;  // HO2
  wdot[12] -= q; // CH3
  wdot[18] += q; // CH3O
  // d()/d[OH]
  dqdci = -k_r * sc[18];
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[170] -= dqdci; // dwdot[HO2]/d[OH]
  J[174] -= dqdci; // dwdot[CH3]/d[OH]
  J[180] += dqdci; // dwdot[CH3O]/d[OH]
  // d()/d[HO2]
  dqdci = +k_f * sc[12];
  J[222] += dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[228] -= dqdci; // dwdot[CH3]/d[HO2]
  J[234] += dqdci; // dwdot[CH3O]/d[HO2]
  // d()/d[CH3]
  dqdci = +k_f * sc[8];
  J[330] += dqdci; // dwdot[OH]/d[CH3]
  J[332] -= dqdci; // dwdot[HO2]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[342] += dqdci; // dwdot[CH3O]/d[CH3]
  // d()/d[CH3O]
  dqdci = -k_r * sc[6];
  J[492] += dqdci; // dwdot[OH]/d[CH3O]
  J[494] -= dqdci; // dwdot[HO2]/d[CH3O]
  J[498] -= dqdci; // dwdot[CH3]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 68: CH3OH + H <=> CH3O + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[17];
  k_f = 0.199 * exp(2.56 * tc[0] - (5183.13165778861) * invT);
  dlnkfdT = 2.56 * invT + (5183.13165778861) * invT2;
  // reverse
  phi_r = sc[1] * sc[18];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[17] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[1] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[17] -= q; // CH3OH
  wdot[18] += q; // CH3O
  // d()/d[H2]
  dqdci = -k_r * sc[18];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[44] -= dqdci; // dwdot[CH3OH]/d[H2]
  J[45] += dqdci; // dwdot[CH3O]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[17];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[71] -= dqdci; // dwdot[CH3OH]/d[H]
  J[72] += dqdci; // dwdot[CH3O]/d[H]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[2];
  J[460] += dqdci; // dwdot[H2]/d[CH3OH]
  J[461] -= dqdci; // dwdot[H]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] += dqdci; // dwdot[CH3O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = -k_r * sc[1];
  J[487] += dqdci; // dwdot[H2]/d[CH3O]
  J[488] -= dqdci; // dwdot[H]/d[CH3O]
  J[503] -= dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 69: CH3OH + H <=> CH2OH + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[17];
  k_f = 0.307 * exp(2.55 * tc[0] - (2737.49866197768) * invT);
  dlnkfdT = 2.55 * invT + (2737.49866197768) * invT2;
  // reverse
  phi_r = sc[1] * sc[19];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[17] - g_RT[19]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[1] + h_RT[19]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[17] -= q; // CH3OH
  wdot[19] += q; // CH2OH
  // d()/d[H2]
  dqdci = -k_r * sc[19];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[44] -= dqdci; // dwdot[CH3OH]/d[H2]
  J[46] += dqdci; // dwdot[CH2OH]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[17];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[71] -= dqdci; // dwdot[CH3OH]/d[H]
  J[73] += dqdci; // dwdot[CH2OH]/d[H]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[2];
  J[460] += dqdci; // dwdot[H2]/d[CH3OH]
  J[461] -= dqdci; // dwdot[H]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[478] += dqdci; // dwdot[CH2OH]/d[CH3OH]
  // d()/d[CH2OH]
  dqdci = -k_r * sc[1];
  J[514] += dqdci; // dwdot[H2]/d[CH2OH]
  J[515] -= dqdci; // dwdot[H]/d[CH2OH]
  J[530] -= dqdci; // dwdot[CH3OH]/d[CH2OH]
  J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT

  // reaction 70: CH3OH + O <=> CH3O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[17];
  k_f = 0.0388 * exp(2.5 * tc[0] - (1549.90733067854) * invT);
  dlnkfdT = 2.5 * invT + (1549.90733067854) * invT2;
  // reverse
  phi_r = sc[6] * sc[18];
  Kc = exp(g_RT[4] - g_RT[6] + g_RT[17] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[17]) + (h_RT[6] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[6] += q;  // OH
  wdot[17] -= q; // CH3OH
  wdot[18] += q; // CH3O
  // d()/d[O]
  dqdci = +k_f * sc[17];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[125] -= dqdci; // dwdot[CH3OH]/d[O]
  J[126] += dqdci; // dwdot[CH3O]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[18];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[179] -= dqdci; // dwdot[CH3OH]/d[OH]
  J[180] += dqdci; // dwdot[CH3O]/d[OH]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[4];
  J[463] -= dqdci; // dwdot[O]/d[CH3OH]
  J[465] += dqdci; // dwdot[OH]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] += dqdci; // dwdot[CH3O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = -k_r * sc[6];
  J[490] -= dqdci; // dwdot[O]/d[CH3O]
  J[492] += dqdci; // dwdot[OH]/d[CH3O]
  J[503] -= dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 71: CH3OH + O <=> CH2OH + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[17];
  k_f = 0.388 * exp(2.5 * tc[0] - (1549.90733067854) * invT);
  dlnkfdT = 2.5 * invT + (1549.90733067854) * invT2;
  // reverse
  phi_r = sc[6] * sc[19];
  Kc = exp(g_RT[4] - g_RT[6] + g_RT[17] - g_RT[19]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[17]) + (h_RT[6] + h_RT[19]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[6] += q;  // OH
  wdot[17] -= q; // CH3OH
  wdot[19] += q; // CH2OH
  // d()/d[O]
  dqdci = +k_f * sc[17];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[125] -= dqdci; // dwdot[CH3OH]/d[O]
  J[127] += dqdci; // dwdot[CH2OH]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[19];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[179] -= dqdci; // dwdot[CH3OH]/d[OH]
  J[181] += dqdci; // dwdot[CH2OH]/d[OH]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[4];
  J[463] -= dqdci; // dwdot[O]/d[CH3OH]
  J[465] += dqdci; // dwdot[OH]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[478] += dqdci; // dwdot[CH2OH]/d[CH3OH]
  // d()/d[CH2OH]
  dqdci = -k_r * sc[6];
  J[517] -= dqdci; // dwdot[O]/d[CH2OH]
  J[519] += dqdci; // dwdot[OH]/d[CH2OH]
  J[530] -= dqdci; // dwdot[CH3OH]/d[CH2OH]
  J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT

  // reaction 72: CH3OH + OH <=> CH3O + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[17];
  k_f = 0.00015 * exp(3.03 * tc[0] - (-383.954316009001) * invT);
  dlnkfdT = 3.03 * invT + (-383.954316009001) * invT2;
  // reverse
  phi_r = sc[5] * sc[18];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[17] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[17]) + (h_RT[5] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[17] -= q; // CH3OH
  wdot[18] += q; // CH3O
  // d()/d[H2O]
  dqdci = -k_r * sc[18];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[152] -= dqdci; // dwdot[CH3OH]/d[H2O]
  J[153] += dqdci; // dwdot[CH3O]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[17];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[179] -= dqdci; // dwdot[CH3OH]/d[OH]
  J[180] += dqdci; // dwdot[CH3O]/d[OH]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[6];
  J[464] += dqdci; // dwdot[H2O]/d[CH3OH]
  J[465] -= dqdci; // dwdot[OH]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] += dqdci; // dwdot[CH3O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = -k_r * sc[5];
  J[491] += dqdci; // dwdot[H2O]/d[CH3O]
  J[492] -= dqdci; // dwdot[OH]/d[CH3O]
  J[503] -= dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 73: CH3OH + OH <=> CH2OH + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[17];
  k_f = 0.0308 * exp(2.65 * tc[0] - (-405.944884304667) * invT);
  dlnkfdT = 2.65 * invT + (-405.944884304667) * invT2;
  // reverse
  phi_r = sc[5] * sc[19];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[17] - g_RT[19]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[17]) + (h_RT[5] + h_RT[19]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[17] -= q; // CH3OH
  wdot[19] += q; // CH2OH
  // d()/d[H2O]
  dqdci = -k_r * sc[19];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[152] -= dqdci; // dwdot[CH3OH]/d[H2O]
  J[154] += dqdci; // dwdot[CH2OH]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[17];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[179] -= dqdci; // dwdot[CH3OH]/d[OH]
  J[181] += dqdci; // dwdot[CH2OH]/d[OH]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[6];
  J[464] += dqdci; // dwdot[H2O]/d[CH3OH]
  J[465] -= dqdci; // dwdot[OH]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[478] += dqdci; // dwdot[CH2OH]/d[CH3OH]
  // d()/d[CH2OH]
  dqdci = -k_r * sc[5];
  J[518] += dqdci; // dwdot[H2O]/d[CH2OH]
  J[519] -= dqdci; // dwdot[OH]/d[CH2OH]
  J[530] -= dqdci; // dwdot[CH3OH]/d[CH2OH]
  J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT

  // reaction 74: CH3OH + O2 <=> CH3O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[17];
  k_f = 0.0358 * exp(2.27 * tc[0] - (21519.8091048059) * invT);
  dlnkfdT = 2.27 * invT + (21519.8091048059) * invT2;
  // reverse
  phi_r = sc[8] * sc[18];
  Kc = exp(g_RT[3] - g_RT[8] + g_RT[17] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[17]) + (h_RT[8] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[8] += q;  // HO2
  wdot[17] -= q; // CH3OH
  wdot[18] += q; // CH3O
  // d()/d[O2]
  dqdci = +k_f * sc[17];
  J[84] -= dqdci; // dwdot[O2]/d[O2]
  J[89] += dqdci; // dwdot[HO2]/d[O2]
  J[98] -= dqdci; // dwdot[CH3OH]/d[O2]
  J[99] += dqdci; // dwdot[CH3O]/d[O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[18];
  J[219] -= dqdci; // dwdot[O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  J[233] -= dqdci; // dwdot[CH3OH]/d[HO2]
  J[234] += dqdci; // dwdot[CH3O]/d[HO2]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[3];
  J[462] -= dqdci; // dwdot[O2]/d[CH3OH]
  J[467] += dqdci; // dwdot[HO2]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] += dqdci; // dwdot[CH3O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = -k_r * sc[8];
  J[489] -= dqdci; // dwdot[O2]/d[CH3O]
  J[494] += dqdci; // dwdot[HO2]/d[CH3O]
  J[503] -= dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 75: CH3OH + HO2 <=> CH3O + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[17];
  k_f = 1220000 * exp(-(10099.9107343668) * invT);
  dlnkfdT = (10099.9107343668) * invT2;
  // reverse
  phi_r = sc[7] * sc[18];
  Kc = exp(-g_RT[7] + g_RT[8] + g_RT[17] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[17]) + (h_RT[7] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[17] -= q; // CH3OH
  wdot[18] += q; // CH3O
  // d()/d[H2O2]
  dqdci = -k_r * sc[18];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[206] -= dqdci; // dwdot[CH3OH]/d[H2O2]
  J[207] += dqdci; // dwdot[CH3O]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[17];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[233] -= dqdci; // dwdot[CH3OH]/d[HO2]
  J[234] += dqdci; // dwdot[CH3O]/d[HO2]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[8];
  J[466] += dqdci; // dwdot[H2O2]/d[CH3OH]
  J[467] -= dqdci; // dwdot[HO2]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] += dqdci; // dwdot[CH3O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = -k_r * sc[7];
  J[493] += dqdci; // dwdot[H2O2]/d[CH3O]
  J[494] -= dqdci; // dwdot[HO2]/d[CH3O]
  J[503] -= dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 76: CH3OH + HO2 <=> CH2OH + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[17];
  k_f = 32600000 * exp(-(9451.5160604774) * invT);
  dlnkfdT = (9451.5160604774) * invT2;
  // reverse
  phi_r = sc[7] * sc[19];
  Kc = exp(-g_RT[7] + g_RT[8] + g_RT[17] - g_RT[19]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[17]) + (h_RT[7] + h_RT[19]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[17] -= q; // CH3OH
  wdot[19] += q; // CH2OH
  // d()/d[H2O2]
  dqdci = -k_r * sc[19];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[206] -= dqdci; // dwdot[CH3OH]/d[H2O2]
  J[208] += dqdci; // dwdot[CH2OH]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[17];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[233] -= dqdci; // dwdot[CH3OH]/d[HO2]
  J[235] += dqdci; // dwdot[CH2OH]/d[HO2]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[8];
  J[466] += dqdci; // dwdot[H2O2]/d[CH3OH]
  J[467] -= dqdci; // dwdot[HO2]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[478] += dqdci; // dwdot[CH2OH]/d[CH3OH]
  // d()/d[CH2OH]
  dqdci = -k_r * sc[7];
  J[520] += dqdci; // dwdot[H2O2]/d[CH2OH]
  J[521] -= dqdci; // dwdot[HO2]/d[CH2OH]
  J[530] -= dqdci; // dwdot[CH3OH]/d[CH2OH]
  J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT

  // reaction 77: CH3 + CH3OH <=> CH2OH + CH4
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[12] * sc[17];
  k_f = 2.13e-07 * exp(3.953 * tc[0] - (3550.24389891888) * invT);
  dlnkfdT = 3.953 * invT + (3550.24389891888) * invT2;
  // reverse
  phi_r = sc[11] * sc[19];
  Kc = exp(-g_RT[11] + g_RT[12] + g_RT[17] - g_RT[19]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[12] + h_RT[17]) + (h_RT[11] + h_RT[19]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  wdot[17] -= q; // CH3OH
  wdot[19] += q; // CH2OH
  // d()/d[CH4]
  dqdci = -k_r * sc[19];
  J[308] += dqdci; // dwdot[CH4]/d[CH4]
  J[309] -= dqdci; // dwdot[CH3]/d[CH4]
  J[314] -= dqdci; // dwdot[CH3OH]/d[CH4]
  J[316] += dqdci; // dwdot[CH2OH]/d[CH4]
  // d()/d[CH3]
  dqdci = +k_f * sc[17];
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[341] -= dqdci; // dwdot[CH3OH]/d[CH3]
  J[343] += dqdci; // dwdot[CH2OH]/d[CH3]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[12];
  J[470] += dqdci; // dwdot[CH4]/d[CH3OH]
  J[471] -= dqdci; // dwdot[CH3]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[478] += dqdci; // dwdot[CH2OH]/d[CH3OH]
  // d()/d[CH2OH]
  dqdci = -k_r * sc[11];
  J[524] += dqdci; // dwdot[CH4]/d[CH2OH]
  J[525] -= dqdci; // dwdot[CH3]/d[CH2OH]
  J[530] -= dqdci; // dwdot[CH3OH]/d[CH2OH]
  J[532] += dqdci; // dwdot[CH2OH]/d[CH2OH]
  // d()/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[721] += dqdT; // dwdot[CH2OH]/dT

  // reaction 78: CH3 + CH3OH <=> CH3O + CH4
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[12] * sc[17];
  k_f = 0.00322 * exp(2.425 * tc[0] - (4317.34738427159) * invT);
  dlnkfdT = 2.425 * invT + (4317.34738427159) * invT2;
  // reverse
  phi_r = sc[11] * sc[18];
  Kc = exp(-g_RT[11] + g_RT[12] + g_RT[17] - g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[12] + h_RT[17]) + (h_RT[11] + h_RT[18]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  wdot[17] -= q; // CH3OH
  wdot[18] += q; // CH3O
  // d()/d[CH4]
  dqdci = -k_r * sc[18];
  J[308] += dqdci; // dwdot[CH4]/d[CH4]
  J[309] -= dqdci; // dwdot[CH3]/d[CH4]
  J[314] -= dqdci; // dwdot[CH3OH]/d[CH4]
  J[315] += dqdci; // dwdot[CH3O]/d[CH4]
  // d()/d[CH3]
  dqdci = +k_f * sc[17];
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[341] -= dqdci; // dwdot[CH3OH]/d[CH3]
  J[342] += dqdci; // dwdot[CH3O]/d[CH3]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[12];
  J[470] += dqdci; // dwdot[CH4]/d[CH3OH]
  J[471] -= dqdci; // dwdot[CH3]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] += dqdci; // dwdot[CH3O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = -k_r * sc[11];
  J[497] += dqdci; // dwdot[CH4]/d[CH3O]
  J[498] -= dqdci; // dwdot[CH3]/d[CH3O]
  J[503] -= dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] += dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT
  J[720] += dqdT; // dwdot[CH3O]/dT

  // reaction 79: CH2OH + O2 <=> CH2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[19];
  k_f = 1510000000 * exp(-1 * tc[0]);
  dlnkfdT = -1 * invT;
  // reverse
  phi_r = sc[8] * sc[20];
  Kc = exp(g_RT[3] - g_RT[8] + g_RT[19] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[19]) + (h_RT[8] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[8] += q;  // HO2
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[O2]
  dqdci = +k_f * sc[19];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[89] += dqdci;  // dwdot[HO2]/d[O2]
  J[100] -= dqdci; // dwdot[CH2OH]/d[O2]
  J[101] += dqdci; // dwdot[CH2O]/d[O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[20];
  J[219] -= dqdci; // dwdot[O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  J[235] -= dqdci; // dwdot[CH2OH]/d[HO2]
  J[236] += dqdci; // dwdot[CH2O]/d[HO2]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[3];
  J[516] -= dqdci; // dwdot[O2]/d[CH2OH]
  J[521] += dqdci; // dwdot[HO2]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[8];
  J[543] -= dqdci; // dwdot[O2]/d[CH2O]
  J[548] += dqdci; // dwdot[HO2]/d[CH2O]
  J[559] -= dqdci; // dwdot[CH2OH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 80: CH2OH + O2 <=> CH2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[19];
  k_f = 241000000 * exp(-(2524.63801234228) * invT);
  dlnkfdT = (2524.63801234228) * invT2;
  // reverse
  phi_r = sc[8] * sc[20];
  Kc = exp(g_RT[3] - g_RT[8] + g_RT[19] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[19]) + (h_RT[8] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[8] += q;  // HO2
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[O2]
  dqdci = +k_f * sc[19];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[89] += dqdci;  // dwdot[HO2]/d[O2]
  J[100] -= dqdci; // dwdot[CH2OH]/d[O2]
  J[101] += dqdci; // dwdot[CH2O]/d[O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[20];
  J[219] -= dqdci; // dwdot[O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  J[235] -= dqdci; // dwdot[CH2OH]/d[HO2]
  J[236] += dqdci; // dwdot[CH2O]/d[HO2]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[3];
  J[516] -= dqdci; // dwdot[O2]/d[CH2OH]
  J[521] += dqdci; // dwdot[HO2]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[8];
  J[543] -= dqdci; // dwdot[O2]/d[CH2O]
  J[548] += dqdci; // dwdot[HO2]/d[CH2O]
  J[559] -= dqdci; // dwdot[CH2OH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 81: CH2OH + H <=> CH2O + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[19];
  k_f = 6000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[20];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[19] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[19]) + (h_RT[1] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[H2]
  dqdci = -k_r * sc[20];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[46] -= dqdci; // dwdot[CH2OH]/d[H2]
  J[47] += dqdci; // dwdot[CH2O]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[19];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[73] -= dqdci; // dwdot[CH2OH]/d[H]
  J[74] += dqdci; // dwdot[CH2O]/d[H]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[2];
  J[514] += dqdci; // dwdot[H2]/d[CH2OH]
  J[515] -= dqdci; // dwdot[H]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[1];
  J[541] += dqdci; // dwdot[H2]/d[CH2O]
  J[542] -= dqdci; // dwdot[H]/d[CH2O]
  J[559] -= dqdci; // dwdot[CH2OH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 82: CH2OH + HO2 <=> CH2O + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[19];
  k_f = 12000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[7] * sc[20];
  Kc = exp(-g_RT[7] + g_RT[8] + g_RT[19] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[19]) + (h_RT[7] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[H2O2]
  dqdci = -k_r * sc[20];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[208] -= dqdci; // dwdot[CH2OH]/d[H2O2]
  J[209] += dqdci; // dwdot[CH2O]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[19];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[235] -= dqdci; // dwdot[CH2OH]/d[HO2]
  J[236] += dqdci; // dwdot[CH2O]/d[HO2]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[8];
  J[520] += dqdci; // dwdot[H2O2]/d[CH2OH]
  J[521] -= dqdci; // dwdot[HO2]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[7];
  J[547] += dqdci; // dwdot[H2O2]/d[CH2O]
  J[548] -= dqdci; // dwdot[HO2]/d[CH2O]
  J[559] -= dqdci; // dwdot[CH2OH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 83: CH2OH + HCO <=> 2 CH2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[19] * sc[21];
  k_f = 180000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = (sc[20] * sc[20]);
  Kc = exp(g_RT[19] - 2.000000 * g_RT[20] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[19] + h_RT[21]) + (2.000000 * h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[19] -= q;     // CH2OH
  wdot[20] += 2 * q; // CH2O
  wdot[21] -= q;     // HCO
  // d()/d[CH2OH]
  dqdci = +k_f * sc[21];
  J[532] -= dqdci;     // dwdot[CH2OH]/d[CH2OH]
  J[533] += 2 * dqdci; // dwdot[CH2O]/d[CH2OH]
  J[534] -= dqdci;     // dwdot[HCO]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * 2.000000 * sc[20];
  J[559] -= dqdci;     // dwdot[CH2OH]/d[CH2O]
  J[560] += 2 * dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] -= dqdci;     // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = +k_f * sc[19];
  J[586] -= dqdci;     // dwdot[CH2OH]/d[HCO]
  J[587] += 2 * dqdci; // dwdot[CH2O]/d[HCO]
  J[588] -= dqdci;     // dwdot[HCO]/d[HCO]
  // d()/dT
  J[721] -= dqdT;     // dwdot[CH2OH]/dT
  J[722] += 2 * dqdT; // dwdot[CH2O]/dT
  J[723] -= dqdT;     // dwdot[HCO]/dT

  // reaction 84: CH2OH + OH <=> CH2O + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[19];
  k_f = 24000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[5] * sc[20];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[19] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[19]) + (h_RT[5] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[H2O]
  dqdci = -k_r * sc[20];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[154] -= dqdci; // dwdot[CH2OH]/d[H2O]
  J[155] += dqdci; // dwdot[CH2O]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[19];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[181] -= dqdci; // dwdot[CH2OH]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[6];
  J[518] += dqdci; // dwdot[H2O]/d[CH2OH]
  J[519] -= dqdci; // dwdot[OH]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[5];
  J[545] += dqdci; // dwdot[H2O]/d[CH2O]
  J[546] -= dqdci; // dwdot[OH]/d[CH2O]
  J[559] -= dqdci; // dwdot[CH2OH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 85: CH2OH + O <=> CH2O + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[19];
  k_f = 42000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[6] * sc[20];
  Kc = exp(g_RT[4] - g_RT[6] + g_RT[19] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[19]) + (h_RT[6] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[6] += q;  // OH
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[O]
  dqdci = +k_f * sc[19];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[127] -= dqdci; // dwdot[CH2OH]/d[O]
  J[128] += dqdci; // dwdot[CH2O]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[20];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[181] -= dqdci; // dwdot[CH2OH]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[4];
  J[517] -= dqdci; // dwdot[O]/d[CH2OH]
  J[519] += dqdci; // dwdot[OH]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/d[CH2O]
  dqdci = -k_r * sc[6];
  J[544] -= dqdci; // dwdot[O]/d[CH2O]
  J[546] += dqdci; // dwdot[OH]/d[CH2O]
  J[559] -= dqdci; // dwdot[CH2OH]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 86: CH3O + O2 <=> CH2O + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[18];
  k_f = 4.38e-25 * exp(9.5 * tc[0] - (-2768.19487859176) * invT);
  dlnkfdT = 9.5 * invT + (-2768.19487859176) * invT2;
  // reverse
  phi_r = sc[8] * sc[20];
  Kc = exp(g_RT[3] - g_RT[8] + g_RT[18] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[18]) + (h_RT[8] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[8] += q;  // HO2
  wdot[18] -= q; // CH3O
  wdot[20] += q; // CH2O
  // d()/d[O2]
  dqdci = +k_f * sc[18];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[89] += dqdci;  // dwdot[HO2]/d[O2]
  J[99] -= dqdci;  // dwdot[CH3O]/d[O2]
  J[101] += dqdci; // dwdot[CH2O]/d[O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[20];
  J[219] -= dqdci; // dwdot[O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  J[234] -= dqdci; // dwdot[CH3O]/d[HO2]
  J[236] += dqdci; // dwdot[CH2O]/d[HO2]
  // d()/d[CH3O]
  dqdci = +k_f * sc[3];
  J[489] -= dqdci; // dwdot[O2]/d[CH3O]
  J[494] += dqdci; // dwdot[HO2]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] += dqdci; // dwdot[CH2O]/d[CH3O]
  // d()/d[CH2O]
  dqdci = -k_r * sc[8];
  J[543] -= dqdci; // dwdot[O2]/d[CH2O]
  J[548] += dqdci; // dwdot[HO2]/d[CH2O]
  J[558] -= dqdci; // dwdot[CH3O]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 87: CH3O + H <=> CH2O + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[18];
  k_f = 20000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[20];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[18] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[18]) + (h_RT[1] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[18] -= q; // CH3O
  wdot[20] += q; // CH2O
  // d()/d[H2]
  dqdci = -k_r * sc[20];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[45] -= dqdci; // dwdot[CH3O]/d[H2]
  J[47] += dqdci; // dwdot[CH2O]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[18];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[72] -= dqdci; // dwdot[CH3O]/d[H]
  J[74] += dqdci; // dwdot[CH2O]/d[H]
  // d()/d[CH3O]
  dqdci = +k_f * sc[2];
  J[487] += dqdci; // dwdot[H2]/d[CH3O]
  J[488] -= dqdci; // dwdot[H]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] += dqdci; // dwdot[CH2O]/d[CH3O]
  // d()/d[CH2O]
  dqdci = -k_r * sc[1];
  J[541] += dqdci; // dwdot[H2]/d[CH2O]
  J[542] -= dqdci; // dwdot[H]/d[CH2O]
  J[558] -= dqdci; // dwdot[CH3O]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 88: CH3O + HO2 <=> CH2O + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[18];
  k_f = 301000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[7] * sc[20];
  Kc = exp(-g_RT[7] + g_RT[8] + g_RT[18] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[18]) + (h_RT[7] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[18] -= q; // CH3O
  wdot[20] += q; // CH2O
  // d()/d[H2O2]
  dqdci = -k_r * sc[20];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[207] -= dqdci; // dwdot[CH3O]/d[H2O2]
  J[209] += dqdci; // dwdot[CH2O]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[18];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[234] -= dqdci; // dwdot[CH3O]/d[HO2]
  J[236] += dqdci; // dwdot[CH2O]/d[HO2]
  // d()/d[CH3O]
  dqdci = +k_f * sc[8];
  J[493] += dqdci; // dwdot[H2O2]/d[CH3O]
  J[494] -= dqdci; // dwdot[HO2]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] += dqdci; // dwdot[CH2O]/d[CH3O]
  // d()/d[CH2O]
  dqdci = -k_r * sc[7];
  J[547] += dqdci; // dwdot[H2O2]/d[CH2O]
  J[548] -= dqdci; // dwdot[HO2]/d[CH2O]
  J[558] -= dqdci; // dwdot[CH3O]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 89: CH3 + CH3O <=> CH2O + CH4
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[12] * sc[18];
  k_f = 12000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[11] * sc[20];
  Kc = exp(-g_RT[11] + g_RT[12] + g_RT[18] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[12] + h_RT[18]) + (h_RT[11] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  wdot[18] -= q; // CH3O
  wdot[20] += q; // CH2O
  // d()/d[CH4]
  dqdci = -k_r * sc[20];
  J[308] += dqdci; // dwdot[CH4]/d[CH4]
  J[309] -= dqdci; // dwdot[CH3]/d[CH4]
  J[315] -= dqdci; // dwdot[CH3O]/d[CH4]
  J[317] += dqdci; // dwdot[CH2O]/d[CH4]
  // d()/d[CH3]
  dqdci = +k_f * sc[18];
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[342] -= dqdci; // dwdot[CH3O]/d[CH3]
  J[344] += dqdci; // dwdot[CH2O]/d[CH3]
  // d()/d[CH3O]
  dqdci = +k_f * sc[12];
  J[497] += dqdci; // dwdot[CH4]/d[CH3O]
  J[498] -= dqdci; // dwdot[CH3]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] += dqdci; // dwdot[CH2O]/d[CH3O]
  // d()/d[CH2O]
  dqdci = -k_r * sc[11];
  J[551] += dqdci; // dwdot[CH4]/d[CH2O]
  J[552] -= dqdci; // dwdot[CH3]/d[CH2O]
  J[558] -= dqdci; // dwdot[CH3O]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 90: 2 CH3O <=> CH2O + CH3OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[18] * sc[18]);
  k_f = 60300000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[17] * sc[20];
  Kc = exp(-g_RT[17] + 2.000000 * g_RT[18] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[18]) + (h_RT[17] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[17] += q;     // CH3OH
  wdot[18] -= 2 * q; // CH3O
  wdot[20] += q;     // CH2O
  // d()/d[CH3OH]
  dqdci = -k_r * sc[20];
  J[476] += dqdci;      // dwdot[CH3OH]/d[CH3OH]
  J[477] += -2 * dqdci; // dwdot[CH3O]/d[CH3OH]
  J[479] += dqdci;      // dwdot[CH2O]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = +k_f * 2.000000 * sc[18];
  J[503] += dqdci;      // dwdot[CH3OH]/d[CH3O]
  J[504] += -2 * dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] += dqdci;      // dwdot[CH2O]/d[CH3O]
  // d()/d[CH2O]
  dqdci = -k_r * sc[17];
  J[557] += dqdci;      // dwdot[CH3OH]/d[CH2O]
  J[558] += -2 * dqdci; // dwdot[CH3O]/d[CH2O]
  J[560] += dqdci;      // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[719] += dqdT;      // dwdot[CH3OH]/dT
  J[720] += -2 * dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT;      // dwdot[CH2O]/dT

  // reaction 91: CH2O + O <=> HCO + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[20];
  k_f = 6260 * exp(1.15 * tc[0] - (1137.26966471867) * invT);
  dlnkfdT = 1.15 * invT + (1137.26966471867) * invT2;
  // reverse
  phi_r = sc[6] * sc[21];
  Kc = exp(g_RT[4] - g_RT[6] + g_RT[20] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[20]) + (h_RT[6] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[6] += q;  // OH
  wdot[20] -= q; // CH2O
  wdot[21] += q; // HCO
  // d()/d[O]
  dqdci = +k_f * sc[20];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[128] -= dqdci; // dwdot[CH2O]/d[O]
  J[129] += dqdci; // dwdot[HCO]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[21];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[182] -= dqdci; // dwdot[CH2O]/d[OH]
  J[183] += dqdci; // dwdot[HCO]/d[OH]
  // d()/d[CH2O]
  dqdci = +k_f * sc[4];
  J[544] -= dqdci; // dwdot[O]/d[CH2O]
  J[546] += dqdci; // dwdot[OH]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] += dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = -k_r * sc[6];
  J[571] -= dqdci; // dwdot[O]/d[HCO]
  J[573] += dqdci; // dwdot[OH]/d[HCO]
  J[587] -= dqdci; // dwdot[CH2O]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 92: CH2O + H <=> H2 + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[20];
  k_f = 57.4 * exp(1.9 * tc[0] - (1378.81366430493) * invT);
  dlnkfdT = 1.9 * invT + (1378.81366430493) * invT2;
  // reverse
  phi_r = sc[1] * sc[21];
  Kc = exp(-g_RT[1] + g_RT[2] + g_RT[20] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[20]) + (h_RT[1] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[20] -= q; // CH2O
  wdot[21] += q; // HCO
  // d()/d[H2]
  dqdci = -k_r * sc[21];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[47] -= dqdci; // dwdot[CH2O]/d[H2]
  J[48] += dqdci; // dwdot[HCO]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[20];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[74] -= dqdci; // dwdot[CH2O]/d[H]
  J[75] += dqdci; // dwdot[HCO]/d[H]
  // d()/d[CH2O]
  dqdci = +k_f * sc[2];
  J[541] += dqdci; // dwdot[H2]/d[CH2O]
  J[542] -= dqdci; // dwdot[H]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] += dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = -k_r * sc[1];
  J[568] += dqdci; // dwdot[H2]/d[HCO]
  J[569] -= dqdci; // dwdot[H]/d[HCO]
  J[587] -= dqdci; // dwdot[CH2O]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 93: CH2O + OH <=> H2O + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[20];
  k_f = 78.2 * exp(1.63 * tc[0] - (-530.893582423979) * invT);
  dlnkfdT = 1.63 * invT + (-530.893582423979) * invT2;
  // reverse
  phi_r = sc[5] * sc[21];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[20] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[20]) + (h_RT[5] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[20] -= q; // CH2O
  wdot[21] += q; // HCO
  // d()/d[H2O]
  dqdci = -k_r * sc[21];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[155] -= dqdci; // dwdot[CH2O]/d[H2O]
  J[156] += dqdci; // dwdot[HCO]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[20];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[182] -= dqdci; // dwdot[CH2O]/d[OH]
  J[183] += dqdci; // dwdot[HCO]/d[OH]
  // d()/d[CH2O]
  dqdci = +k_f * sc[6];
  J[545] += dqdci; // dwdot[H2O]/d[CH2O]
  J[546] -= dqdci; // dwdot[OH]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] += dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = -k_r * sc[5];
  J[572] += dqdci; // dwdot[H2O]/d[HCO]
  J[573] -= dqdci; // dwdot[OH]/d[HCO]
  J[587] -= dqdci; // dwdot[CH2O]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 94: CH2O + HO2 <=> H2O2 + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[20];
  k_f = 0.0188 * exp(2.7 * tc[0] - (5797.05599007037) * invT);
  dlnkfdT = 2.7 * invT + (5797.05599007037) * invT2;
  // reverse
  phi_r = sc[7] * sc[21];
  Kc = exp(-g_RT[7] + g_RT[8] + g_RT[20] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[20]) + (h_RT[7] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[20] -= q; // CH2O
  wdot[21] += q; // HCO
  // d()/d[H2O2]
  dqdci = -k_r * sc[21];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[209] -= dqdci; // dwdot[CH2O]/d[H2O2]
  J[210] += dqdci; // dwdot[HCO]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[20];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[236] -= dqdci; // dwdot[CH2O]/d[HO2]
  J[237] += dqdci; // dwdot[HCO]/d[HO2]
  // d()/d[CH2O]
  dqdci = +k_f * sc[8];
  J[547] += dqdci; // dwdot[H2O2]/d[CH2O]
  J[548] -= dqdci; // dwdot[HO2]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] += dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = -k_r * sc[7];
  J[574] += dqdci; // dwdot[H2O2]/d[HCO]
  J[575] -= dqdci; // dwdot[HO2]/d[HCO]
  J[587] -= dqdci; // dwdot[CH2O]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 95: CH2O + CH3 <=> CH4 + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[12] * sc[20];
  k_f = 3.83e-05 * exp(3.36 * tc[0] - (2169.87026294995) * invT);
  dlnkfdT = 3.36 * invT + (2169.87026294995) * invT2;
  // reverse
  phi_r = sc[11] * sc[21];
  Kc = exp(-g_RT[11] + g_RT[12] + g_RT[20] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[12] + h_RT[20]) + (h_RT[11] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  wdot[20] -= q; // CH2O
  wdot[21] += q; // HCO
  // d()/d[CH4]
  dqdci = -k_r * sc[21];
  J[308] += dqdci; // dwdot[CH4]/d[CH4]
  J[309] -= dqdci; // dwdot[CH3]/d[CH4]
  J[317] -= dqdci; // dwdot[CH2O]/d[CH4]
  J[318] += dqdci; // dwdot[HCO]/d[CH4]
  // d()/d[CH3]
  dqdci = +k_f * sc[20];
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[344] -= dqdci; // dwdot[CH2O]/d[CH3]
  J[345] += dqdci; // dwdot[HCO]/d[CH3]
  // d()/d[CH2O]
  dqdci = +k_f * sc[12];
  J[551] += dqdci; // dwdot[CH4]/d[CH2O]
  J[552] -= dqdci; // dwdot[CH3]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] += dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = -k_r * sc[11];
  J[578] += dqdci; // dwdot[CH4]/d[HCO]
  J[579] -= dqdci; // dwdot[CH3]/d[HCO]
  J[587] -= dqdci; // dwdot[CH2O]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 96: CH2O + CH3O <=> CH3OH + HCO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[18] * sc[20];
  k_f = 662000 * exp(-(1154.37903135603) * invT);
  dlnkfdT = (1154.37903135603) * invT2;
  // reverse
  phi_r = sc[17] * sc[21];
  Kc = exp(-g_RT[17] + g_RT[18] + g_RT[20] - g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[18] + h_RT[20]) + (h_RT[17] + h_RT[21]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[17] += q; // CH3OH
  wdot[18] -= q; // CH3O
  wdot[20] -= q; // CH2O
  wdot[21] += q; // HCO
  // d()/d[CH3OH]
  dqdci = -k_r * sc[21];
  J[476] += dqdci; // dwdot[CH3OH]/d[CH3OH]
  J[477] -= dqdci; // dwdot[CH3O]/d[CH3OH]
  J[479] -= dqdci; // dwdot[CH2O]/d[CH3OH]
  J[480] += dqdci; // dwdot[HCO]/d[CH3OH]
  // d()/d[CH3O]
  dqdci = +k_f * sc[20];
  J[503] += dqdci; // dwdot[CH3OH]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] -= dqdci; // dwdot[CH2O]/d[CH3O]
  J[507] += dqdci; // dwdot[HCO]/d[CH3O]
  // d()/d[CH2O]
  dqdci = +k_f * sc[18];
  J[557] += dqdci; // dwdot[CH3OH]/d[CH2O]
  J[558] -= dqdci; // dwdot[CH3O]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  J[561] += dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = -k_r * sc[17];
  J[584] += dqdci; // dwdot[CH3OH]/d[HCO]
  J[585] -= dqdci; // dwdot[CH3O]/d[HCO]
  J[587] -= dqdci; // dwdot[CH2O]/d[HCO]
  J[588] += dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[719] += dqdT; // dwdot[CH3OH]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT
  J[723] += dqdT; // dwdot[HCO]/dT

  // reaction 98: HCO + O2 <=> CO + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[3] * sc[21];
  k_f = 7580000 * exp(-(206.318832979935) * invT);
  dlnkfdT = (206.318832979935) * invT2;
  // reverse
  phi_r = sc[8] * sc[9];
  Kc = exp(g_RT[3] - g_RT[8] - g_RT[9] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[21]) + (h_RT[8] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[3] -= q;  // O2
  wdot[8] += q;  // HO2
  wdot[9] += q;  // CO
  wdot[21] -= q; // HCO
  // d()/d[O2]
  dqdci = +k_f * sc[21];
  J[84] -= dqdci;  // dwdot[O2]/d[O2]
  J[89] += dqdci;  // dwdot[HO2]/d[O2]
  J[90] += dqdci;  // dwdot[CO]/d[O2]
  J[102] -= dqdci; // dwdot[HCO]/d[O2]
  // d()/d[HO2]
  dqdci = -k_r * sc[9];
  J[219] -= dqdci; // dwdot[O2]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  J[225] += dqdci; // dwdot[CO]/d[HO2]
  J[237] -= dqdci; // dwdot[HCO]/d[HO2]
  // d()/d[CO]
  dqdci = -k_r * sc[8];
  J[246] -= dqdci; // dwdot[O2]/d[CO]
  J[251] += dqdci; // dwdot[HO2]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[3];
  J[570] -= dqdci; // dwdot[O2]/d[HCO]
  J[575] += dqdci; // dwdot[HO2]/d[HCO]
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[705] -= dqdT; // dwdot[O2]/dT
  J[710] += dqdT; // dwdot[HO2]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 99: HCO + O <=> CO + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[21];
  k_f = 30200000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[6] * sc[9];
  Kc = exp(g_RT[4] - g_RT[6] - g_RT[9] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[21]) + (h_RT[6] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[6] += q;  // OH
  wdot[9] += q;  // CO
  wdot[21] -= q; // HCO
  // d()/d[O]
  dqdci = +k_f * sc[21];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[114] += dqdci; // dwdot[OH]/d[O]
  J[117] += dqdci; // dwdot[CO]/d[O]
  J[129] -= dqdci; // dwdot[HCO]/d[O]
  // d()/d[OH]
  dqdci = -k_r * sc[9];
  J[166] -= dqdci; // dwdot[O]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[171] += dqdci; // dwdot[CO]/d[OH]
  J[183] -= dqdci; // dwdot[HCO]/d[OH]
  // d()/d[CO]
  dqdci = -k_r * sc[6];
  J[247] -= dqdci; // dwdot[O]/d[CO]
  J[249] += dqdci; // dwdot[OH]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[4];
  J[571] -= dqdci; // dwdot[O]/d[HCO]
  J[573] += dqdci; // dwdot[OH]/d[HCO]
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 100: H + HCO <=> CO + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[21];
  k_f = 73400000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[9];
  Kc = exp(-g_RT[1] + g_RT[2] - g_RT[9] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[21]) + (h_RT[1] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] -= q;  // H
  wdot[9] += q;  // CO
  wdot[21] -= q; // HCO
  // d()/d[H2]
  dqdci = -k_r * sc[9];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[29] -= dqdci; // dwdot[H]/d[H2]
  J[36] += dqdci; // dwdot[CO]/d[H2]
  J[48] -= dqdci; // dwdot[HCO]/d[H2]
  // d()/d[H]
  dqdci = +k_f * sc[21];
  J[55] += dqdci; // dwdot[H2]/d[H]
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[63] += dqdci; // dwdot[CO]/d[H]
  J[75] -= dqdci; // dwdot[HCO]/d[H]
  // d()/d[CO]
  dqdci = -k_r * sc[1];
  J[244] += dqdci; // dwdot[H2]/d[CO]
  J[245] -= dqdci; // dwdot[H]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[2];
  J[568] += dqdci; // dwdot[H2]/d[HCO]
  J[569] -= dqdci; // dwdot[H]/d[HCO]
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 101: HCO + OH <=> CO + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[21];
  k_f = 30110000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[5] * sc[9];
  Kc = exp(-g_RT[5] + g_RT[6] - g_RT[9] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[21]) + (h_RT[5] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[9] += q;  // CO
  wdot[21] -= q; // HCO
  // d()/d[H2O]
  dqdci = -k_r * sc[9];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[144] += dqdci; // dwdot[CO]/d[H2O]
  J[156] -= dqdci; // dwdot[HCO]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[21];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[171] += dqdci; // dwdot[CO]/d[OH]
  J[183] -= dqdci; // dwdot[HCO]/d[OH]
  // d()/d[CO]
  dqdci = -k_r * sc[5];
  J[248] += dqdci; // dwdot[H2O]/d[CO]
  J[249] -= dqdci; // dwdot[OH]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[6];
  J[572] += dqdci; // dwdot[H2O]/d[HCO]
  J[573] -= dqdci; // dwdot[OH]/d[HCO]
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 102: CH3 + HCO <=> CH4 + CO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[12] * sc[21];
  k_f = 26500000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[9] * sc[11];
  Kc = exp(-g_RT[9] - g_RT[11] + g_RT[12] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[12] + h_RT[21]) + (h_RT[9] + h_RT[11]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] += q;  // CO
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  wdot[21] -= q; // HCO
  // d()/d[CO]
  dqdci = -k_r * sc[11];
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[254] += dqdci; // dwdot[CH4]/d[CO]
  J[255] -= dqdci; // dwdot[CH3]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[CH4]
  dqdci = -k_r * sc[9];
  J[306] += dqdci; // dwdot[CO]/d[CH4]
  J[308] += dqdci; // dwdot[CH4]/d[CH4]
  J[309] -= dqdci; // dwdot[CH3]/d[CH4]
  J[318] -= dqdci; // dwdot[HCO]/d[CH4]
  // d()/d[CH3]
  dqdci = +k_f * sc[21];
  J[333] += dqdci; // dwdot[CO]/d[CH3]
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[345] -= dqdci; // dwdot[HCO]/d[CH3]
  // d()/d[HCO]
  dqdci = +k_f * sc[12];
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[578] += dqdci; // dwdot[CH4]/d[HCO]
  J[579] -= dqdci; // dwdot[CH3]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 103: 2 HCO <=> CH2O + CO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = (sc[21] * sc[21]);
  k_f = 18000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[9] * sc[20];
  Kc = exp(-g_RT[9] - g_RT[20] + 2.000000 * g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2.000000 * h_RT[21]) + (h_RT[9] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] += q;      // CO
  wdot[20] += q;     // CH2O
  wdot[21] -= 2 * q; // HCO
  // d()/d[CO]
  dqdci = -k_r * sc[20];
  J[252] += dqdci;      // dwdot[CO]/d[CO]
  J[263] += dqdci;      // dwdot[CH2O]/d[CO]
  J[264] += -2 * dqdci; // dwdot[HCO]/d[CO]
  // d()/d[CH2O]
  dqdci = -k_r * sc[9];
  J[549] += dqdci;      // dwdot[CO]/d[CH2O]
  J[560] += dqdci;      // dwdot[CH2O]/d[CH2O]
  J[561] += -2 * dqdci; // dwdot[HCO]/d[CH2O]
  // d()/d[HCO]
  dqdci = +k_f * 2.000000 * sc[21];
  J[576] += dqdci;      // dwdot[CO]/d[HCO]
  J[587] += dqdci;      // dwdot[CH2O]/d[HCO]
  J[588] += -2 * dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[711] += dqdT;      // dwdot[CO]/dT
  J[722] += dqdT;      // dwdot[CH2O]/dT
  J[723] += -2 * dqdT; // dwdot[HCO]/dT

  // reaction 104: HCO + O <=> CO2 + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[21];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[10];
  Kc = exp(-g_RT[2] + g_RT[4] - g_RT[10] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[21]) + (h_RT[2] + h_RT[10]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[4] -= q;  // O
  wdot[10] += q; // CO2
  wdot[21] -= q; // HCO
  // d()/d[H]
  dqdci = -k_r * sc[10];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[58] -= dqdci; // dwdot[O]/d[H]
  J[64] += dqdci; // dwdot[CO2]/d[H]
  J[75] -= dqdci; // dwdot[HCO]/d[H]
  // d()/d[O]
  dqdci = +k_f * sc[21];
  J[110] += dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[118] += dqdci; // dwdot[CO2]/d[O]
  J[129] -= dqdci; // dwdot[HCO]/d[O]
  // d()/d[CO2]
  dqdci = -k_r * sc[2];
  J[272] += dqdci; // dwdot[H]/d[CO2]
  J[274] -= dqdci; // dwdot[O]/d[CO2]
  J[280] += dqdci; // dwdot[CO2]/d[CO2]
  J[291] -= dqdci; // dwdot[HCO]/d[CO2]
  // d()/d[HCO]
  dqdci = +k_f * sc[4];
  J[569] += dqdci; // dwdot[H]/d[HCO]
  J[571] -= dqdci; // dwdot[O]/d[HCO]
  J[577] += dqdci; // dwdot[CO2]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[712] += dqdT; // dwdot[CO2]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 105: HCO + HO2 => CO2 + H + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[21];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[2] += q;  // H
  wdot[6] += q;  // OH
  wdot[8] -= q;  // HO2
  wdot[10] += q; // CO2
  wdot[21] -= q; // HCO
  // d()/d[HO2]
  dqdci = +k_f * sc[21];
  J[218] += dqdci; // dwdot[H]/d[HO2]
  J[222] += dqdci; // dwdot[OH]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[226] += dqdci; // dwdot[CO2]/d[HO2]
  J[237] -= dqdci; // dwdot[HCO]/d[HO2]
  // d()/d[HCO]
  dqdci = +k_f * sc[8];
  J[569] += dqdci; // dwdot[H]/d[HCO]
  J[573] += dqdci; // dwdot[OH]/d[HCO]
  J[575] -= dqdci; // dwdot[HO2]/d[HCO]
  J[577] += dqdci; // dwdot[CO2]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[712] += dqdT; // dwdot[CO2]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 109: CH + OH <=> C + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[16];
  k_f = 40 * exp(2 * tc[0] - (1509.64999741416) * invT);
  dlnkfdT = 2 * invT + (1509.64999741416) * invT2;
  // reverse
  phi_r = sc[5] * sc[15];
  Kc = exp(-g_RT[5] + g_RT[6] - g_RT[15] + g_RT[16]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[16]) + (h_RT[5] + h_RT[15]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[15] += q; // C
  wdot[16] -= q; // CH
  // d()/d[H2O]
  dqdci = -k_r * sc[15];
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[141] -= dqdci; // dwdot[OH]/d[H2O]
  J[150] += dqdci; // dwdot[C]/d[H2O]
  J[151] -= dqdci; // dwdot[CH]/d[H2O]
  // d()/d[OH]
  dqdci = +k_f * sc[16];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[177] += dqdci; // dwdot[C]/d[OH]
  J[178] -= dqdci; // dwdot[CH]/d[OH]
  // d()/d[C]
  dqdci = -k_r * sc[5];
  J[410] += dqdci; // dwdot[H2O]/d[C]
  J[411] -= dqdci; // dwdot[OH]/d[C]
  J[420] += dqdci; // dwdot[C]/d[C]
  J[421] -= dqdci; // dwdot[CH]/d[C]
  // d()/d[CH]
  dqdci = +k_f * sc[6];
  J[437] += dqdci; // dwdot[H2O]/d[CH]
  J[438] -= dqdci; // dwdot[OH]/d[CH]
  J[447] += dqdci; // dwdot[C]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[717] += dqdT; // dwdot[C]/dT
  J[718] -= dqdT; // dwdot[CH]/dT

  // reaction 110: CH2 + HCO <=> CH3 + CO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[13] * sc[21];
  k_f = 20000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[9] * sc[12];
  Kc = exp(-g_RT[9] - g_RT[12] + g_RT[13] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[13] + h_RT[21]) + (h_RT[9] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] += q;  // CO
  wdot[12] += q; // CH3
  wdot[13] -= q; // CH2
  wdot[21] -= q; // HCO
  // d()/d[CO]
  dqdci = -k_r * sc[12];
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[255] += dqdci; // dwdot[CH3]/d[CO]
  J[256] -= dqdci; // dwdot[CH2]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[CH3]
  dqdci = -k_r * sc[9];
  J[333] += dqdci; // dwdot[CO]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  J[337] -= dqdci; // dwdot[CH2]/d[CH3]
  J[345] -= dqdci; // dwdot[HCO]/d[CH3]
  // d()/d[CH2]
  dqdci = +k_f * sc[21];
  J[360] += dqdci; // dwdot[CO]/d[CH2]
  J[363] += dqdci; // dwdot[CH3]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  J[372] -= dqdci; // dwdot[HCO]/d[CH2]
  // d()/d[HCO]
  dqdci = +k_f * sc[13];
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[579] += dqdci; // dwdot[CH3]/d[HCO]
  J[580] -= dqdci; // dwdot[CH2]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[714] += dqdT; // dwdot[CH3]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 111: CH2 + O <=> CO + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[13];
  k_f = 50000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[1] * sc[9];
  Kc = exp(-g_RT[1] + g_RT[4] - g_RT[9] + g_RT[13]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[1] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[4] -= q;  // O
  wdot[9] += q;  // CO
  wdot[13] -= q; // CH2
  // d()/d[H2]
  dqdci = -k_r * sc[9];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[31] -= dqdci; // dwdot[O]/d[H2]
  J[36] += dqdci; // dwdot[CO]/d[H2]
  J[40] -= dqdci; // dwdot[CH2]/d[H2]
  // d()/d[O]
  dqdci = +k_f * sc[13];
  J[109] += dqdci; // dwdot[H2]/d[O]
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[117] += dqdci; // dwdot[CO]/d[O]
  J[121] -= dqdci; // dwdot[CH2]/d[O]
  // d()/d[CO]
  dqdci = -k_r * sc[1];
  J[244] += dqdci; // dwdot[H2]/d[CO]
  J[247] -= dqdci; // dwdot[O]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[256] -= dqdci; // dwdot[CH2]/d[CO]
  // d()/d[CH2]
  dqdci = +k_f * sc[4];
  J[352] += dqdci; // dwdot[H2]/d[CH2]
  J[355] -= dqdci; // dwdot[O]/d[CH2]
  J[360] += dqdci; // dwdot[CO]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT

  // reaction 112: CH2 + OH <=> CH2O + H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[13];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[2] * sc[20];
  Kc = exp(-g_RT[2] + g_RT[6] + g_RT[13] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[6] + h_RT[13]) + (h_RT[2] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[13] -= q; // CH2
  wdot[20] += q; // CH2O
  // d()/d[H]
  dqdci = -k_r * sc[20];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[60] -= dqdci; // dwdot[OH]/d[H]
  J[67] -= dqdci; // dwdot[CH2]/d[H]
  J[74] += dqdci; // dwdot[CH2O]/d[H]
  // d()/d[OH]
  dqdci = +k_f * sc[13];
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[175] -= dqdci; // dwdot[CH2]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH2]
  dqdci = +k_f * sc[6];
  J[353] += dqdci; // dwdot[H]/d[CH2]
  J[357] -= dqdci; // dwdot[OH]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  J[371] += dqdci; // dwdot[CH2O]/d[CH2]
  // d()/d[CH2O]
  dqdci = -k_r * sc[2];
  J[542] += dqdci; // dwdot[H]/d[CH2O]
  J[546] -= dqdci; // dwdot[OH]/d[CH2O]
  J[553] -= dqdci; // dwdot[CH2]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 113: CH2 + CO2 <=> CH2O + CO
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[10] * sc[13];
  k_f = 110000 * exp(-(503.21666580472) * invT);
  dlnkfdT = (503.21666580472) * invT2;
  // reverse
  phi_r = sc[9] * sc[20];
  Kc = exp(-g_RT[9] + g_RT[10] + g_RT[13] - g_RT[20]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[10] + h_RT[13]) + (h_RT[9] + h_RT[20]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] += q;  // CO
  wdot[10] -= q; // CO2
  wdot[13] -= q; // CH2
  wdot[20] += q; // CH2O
  // d()/d[CO]
  dqdci = -k_r * sc[20];
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[253] -= dqdci; // dwdot[CO2]/d[CO]
  J[256] -= dqdci; // dwdot[CH2]/d[CO]
  J[263] += dqdci; // dwdot[CH2O]/d[CO]
  // d()/d[CO2]
  dqdci = +k_f * sc[13];
  J[279] += dqdci; // dwdot[CO]/d[CO2]
  J[280] -= dqdci; // dwdot[CO2]/d[CO2]
  J[283] -= dqdci; // dwdot[CH2]/d[CO2]
  J[290] += dqdci; // dwdot[CH2O]/d[CO2]
  // d()/d[CH2]
  dqdci = +k_f * sc[10];
  J[360] += dqdci; // dwdot[CO]/d[CH2]
  J[361] -= dqdci; // dwdot[CO2]/d[CH2]
  J[364] -= dqdci; // dwdot[CH2]/d[CH2]
  J[371] += dqdci; // dwdot[CH2O]/d[CH2]
  // d()/d[CH2O]
  dqdci = -k_r * sc[9];
  J[549] += dqdci; // dwdot[CO]/d[CH2O]
  J[550] -= dqdci; // dwdot[CO2]/d[CH2O]
  J[553] -= dqdci; // dwdot[CH2]/d[CH2O]
  J[560] += dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[712] -= dqdT; // dwdot[CO2]/dT
  J[715] -= dqdT; // dwdot[CH2]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 114: CH2(S) + O <=> CO + 2 H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[14];
  k_f = 30000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = (sc[2] * sc[2]) * sc[9];
  Kc = refC * exp(-2.000000 * g_RT[2] + g_RT[4] - g_RT[9] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT =
    invT * (-(h_RT[4] + h_RT[14]) + (2.000000 * h_RT[2] + h_RT[9]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += 2 * q; // H
  wdot[4] -= q;     // O
  wdot[9] += q;     // CO
  wdot[14] -= q;    // CH2(S)
  // d()/d[H]
  dqdci = -k_r * 2.000000 * sc[2] * sc[9];
  J[56] += 2 * dqdci; // dwdot[H]/d[H]
  J[58] -= dqdci;     // dwdot[O]/d[H]
  J[63] += dqdci;     // dwdot[CO]/d[H]
  J[68] -= dqdci;     // dwdot[CH2(S)]/d[H]
  // d()/d[O]
  dqdci = +k_f * sc[14];
  J[110] += 2 * dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci;     // dwdot[O]/d[O]
  J[117] += dqdci;     // dwdot[CO]/d[O]
  J[122] -= dqdci;     // dwdot[CH2(S)]/d[O]
  // d()/d[CO]
  dqdci = -k_r * sc[2] * sc[2];
  J[245] += 2 * dqdci; // dwdot[H]/d[CO]
  J[247] -= dqdci;     // dwdot[O]/d[CO]
  J[252] += dqdci;     // dwdot[CO]/d[CO]
  J[257] -= dqdci;     // dwdot[CH2(S)]/d[CO]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[4];
  J[380] += 2 * dqdci; // dwdot[H]/d[CH2(S)]
  J[382] -= dqdci;     // dwdot[O]/d[CH2(S)]
  J[387] += dqdci;     // dwdot[CO]/d[CH2(S)]
  J[392] -= dqdci;     // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[704] += 2 * dqdT; // dwdot[H]/dT
  J[706] -= dqdT;     // dwdot[O]/dT
  J[711] += dqdT;     // dwdot[CO]/dT
  J[716] -= dqdT;     // dwdot[CH2(S)]/dT

  // reaction 115: CH2(S) + CH4 <=> 2 CH3
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[11] * sc[14];
  k_f = 43000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = (sc[12] * sc[12]);
  Kc = exp(g_RT[11] - 2.000000 * g_RT[12] + g_RT[14]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[11] + h_RT[14]) + (2.000000 * h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[11] -= q;     // CH4
  wdot[12] += 2 * q; // CH3
  wdot[14] -= q;     // CH2(S)
  // d()/d[CH4]
  dqdci = +k_f * sc[14];
  J[308] -= dqdci;     // dwdot[CH4]/d[CH4]
  J[309] += 2 * dqdci; // dwdot[CH3]/d[CH4]
  J[311] -= dqdci;     // dwdot[CH2(S)]/d[CH4]
  // d()/d[CH3]
  dqdci = -k_r * 2.000000 * sc[12];
  J[335] -= dqdci;     // dwdot[CH4]/d[CH3]
  J[336] += 2 * dqdci; // dwdot[CH3]/d[CH3]
  J[338] -= dqdci;     // dwdot[CH2(S)]/d[CH3]
  // d()/d[CH2(S)]
  dqdci = +k_f * sc[11];
  J[389] -= dqdci;     // dwdot[CH4]/d[CH2(S)]
  J[390] += 2 * dqdci; // dwdot[CH3]/d[CH2(S)]
  J[392] -= dqdci;     // dwdot[CH2(S)]/d[CH2(S)]
  // d()/dT
  J[713] -= dqdT;     // dwdot[CH4]/dT
  J[714] += 2 * dqdT; // dwdot[CH3]/dT
  J[716] -= dqdT;     // dwdot[CH2(S)]/dT

  // reaction 116: CH3OH + H => CH3 + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[17];
  k_f = 650000 * exp(-(2667.04832876501) * invT);
  dlnkfdT = (2667.04832876501) * invT2;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[2] -= q;  // H
  wdot[5] += q;  // H2O
  wdot[12] += q; // CH3
  wdot[17] -= q; // CH3OH
  // d()/d[H]
  dqdci = +k_f * sc[17];
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[59] += dqdci; // dwdot[H2O]/d[H]
  J[66] += dqdci; // dwdot[CH3]/d[H]
  J[71] -= dqdci; // dwdot[CH3OH]/d[H]
  // d()/d[CH3OH]
  dqdci = +k_f * sc[2];
  J[461] -= dqdci; // dwdot[H]/d[CH3OH]
  J[464] += dqdci; // dwdot[H2O]/d[CH3OH]
  J[471] += dqdci; // dwdot[CH3]/d[CH3OH]
  J[476] -= dqdci; // dwdot[CH3OH]/d[CH3OH]
  // d()/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[714] += dqdT; // dwdot[CH3]/dT
  J[719] -= dqdT; // dwdot[CH3OH]/dT

  // reaction 117: CH2O + O => CO2 + 2 H
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[20];
  k_f = 20000000 * exp(-(2516.0833290236) * invT);
  dlnkfdT = (2516.0833290236) * invT2;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[2] += 2 * q; // H
  wdot[4] -= q;     // O
  wdot[10] += q;    // CO2
  wdot[20] -= q;    // CH2O
  // d()/d[O]
  dqdci = +k_f * sc[20];
  J[110] += 2 * dqdci; // dwdot[H]/d[O]
  J[112] -= dqdci;     // dwdot[O]/d[O]
  J[118] += dqdci;     // dwdot[CO2]/d[O]
  J[128] -= dqdci;     // dwdot[CH2O]/d[O]
  // d()/d[CH2O]
  dqdci = +k_f * sc[4];
  J[542] += 2 * dqdci; // dwdot[H]/d[CH2O]
  J[544] -= dqdci;     // dwdot[O]/d[CH2O]
  J[550] += dqdci;     // dwdot[CO2]/d[CH2O]
  J[560] -= dqdci;     // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[704] += 2 * dqdT; // dwdot[H]/dT
  J[706] -= dqdT;     // dwdot[O]/dT
  J[712] += dqdT;     // dwdot[CO2]/dT
  J[722] -= dqdT;     // dwdot[CH2O]/dT

  // reaction 118: CH2O + OH => CO2 + H + H2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[20];
  k_f = 100000;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[1] += q;  // H2
  wdot[2] += q;  // H
  wdot[6] -= q;  // OH
  wdot[10] += q; // CO2
  wdot[20] -= q; // CH2O
  // d()/d[OH]
  dqdci = +k_f * sc[20];
  J[163] += dqdci; // dwdot[H2]/d[OH]
  J[164] += dqdci; // dwdot[H]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[172] += dqdci; // dwdot[CO2]/d[OH]
  J[182] -= dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH2O]
  dqdci = +k_f * sc[6];
  J[541] += dqdci; // dwdot[H2]/d[CH2O]
  J[542] += dqdci; // dwdot[H]/d[CH2O]
  J[546] -= dqdci; // dwdot[OH]/d[CH2O]
  J[550] += dqdci; // dwdot[CO2]/d[CH2O]
  J[560] -= dqdci; // dwdot[CH2O]/d[CH2O]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[712] += dqdT; // dwdot[CO2]/dT
  J[722] -= dqdT; // dwdot[CH2O]/dT

  // reaction 119: CH3O + CO <=> CH3 + CO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[9] * sc[18];
  k_f = 500000 * exp(-(3270.90832773068) * invT);
  dlnkfdT = (3270.90832773068) * invT2;
  // reverse
  phi_r = sc[10] * sc[12];
  Kc = exp(g_RT[9] - g_RT[10] - g_RT[12] + g_RT[18]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[9] + h_RT[18]) + (h_RT[10] + h_RT[12]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[9] -= q;  // CO
  wdot[10] += q; // CO2
  wdot[12] += q; // CH3
  wdot[18] -= q; // CH3O
  // d()/d[CO]
  dqdci = +k_f * sc[18];
  J[252] -= dqdci; // dwdot[CO]/d[CO]
  J[253] += dqdci; // dwdot[CO2]/d[CO]
  J[255] += dqdci; // dwdot[CH3]/d[CO]
  J[261] -= dqdci; // dwdot[CH3O]/d[CO]
  // d()/d[CO2]
  dqdci = -k_r * sc[12];
  J[279] -= dqdci; // dwdot[CO]/d[CO2]
  J[280] += dqdci; // dwdot[CO2]/d[CO2]
  J[282] += dqdci; // dwdot[CH3]/d[CO2]
  J[288] -= dqdci; // dwdot[CH3O]/d[CO2]
  // d()/d[CH3]
  dqdci = -k_r * sc[10];
  J[333] -= dqdci; // dwdot[CO]/d[CH3]
  J[334] += dqdci; // dwdot[CO2]/d[CH3]
  J[336] += dqdci; // dwdot[CH3]/d[CH3]
  J[342] -= dqdci; // dwdot[CH3O]/d[CH3]
  // d()/d[CH3O]
  dqdci = +k_f * sc[9];
  J[495] -= dqdci; // dwdot[CO]/d[CH3O]
  J[496] += dqdci; // dwdot[CO2]/d[CH3O]
  J[498] += dqdci; // dwdot[CH3]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  // d()/dT
  J[711] -= dqdT; // dwdot[CO]/dT
  J[712] += dqdT; // dwdot[CO2]/dT
  J[714] += dqdT; // dwdot[CH3]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT

  // reaction 121: HCO + HO2 <=> CO + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[8] * sc[21];
  k_f = 400000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[7] * sc[9];
  Kc = exp(-g_RT[7] + g_RT[8] - g_RT[9] + g_RT[21]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[8] + h_RT[21]) + (h_RT[7] + h_RT[9]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[7] += q;  // H2O2
  wdot[8] -= q;  // HO2
  wdot[9] += q;  // CO
  wdot[21] -= q; // HCO
  // d()/d[H2O2]
  dqdci = -k_r * sc[9];
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[197] -= dqdci; // dwdot[HO2]/d[H2O2]
  J[198] += dqdci; // dwdot[CO]/d[H2O2]
  J[210] -= dqdci; // dwdot[HCO]/d[H2O2]
  // d()/d[HO2]
  dqdci = +k_f * sc[21];
  J[223] += dqdci; // dwdot[H2O2]/d[HO2]
  J[224] -= dqdci; // dwdot[HO2]/d[HO2]
  J[225] += dqdci; // dwdot[CO]/d[HO2]
  J[237] -= dqdci; // dwdot[HCO]/d[HO2]
  // d()/d[CO]
  dqdci = -k_r * sc[7];
  J[250] += dqdci; // dwdot[H2O2]/d[CO]
  J[251] -= dqdci; // dwdot[HO2]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[264] -= dqdci; // dwdot[HCO]/d[CO]
  // d()/d[HCO]
  dqdci = +k_f * sc[8];
  J[574] += dqdci; // dwdot[H2O2]/d[HCO]
  J[575] -= dqdci; // dwdot[HO2]/d[HCO]
  J[576] += dqdci; // dwdot[CO]/d[HCO]
  J[588] -= dqdci; // dwdot[HCO]/d[HCO]
  // d()/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[710] -= dqdT; // dwdot[HO2]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[723] -= dqdT; // dwdot[HCO]/dT

  // reaction 122: CH3O + OH => CH2O + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[6] * sc[18];
  k_f = 15000000;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[5] += q;  // H2O
  wdot[6] -= q;  // OH
  wdot[18] -= q; // CH3O
  wdot[20] += q; // CH2O
  // d()/d[OH]
  dqdci = +k_f * sc[18];
  J[167] += dqdci; // dwdot[H2O]/d[OH]
  J[168] -= dqdci; // dwdot[OH]/d[OH]
  J[180] -= dqdci; // dwdot[CH3O]/d[OH]
  J[182] += dqdci; // dwdot[CH2O]/d[OH]
  // d()/d[CH3O]
  dqdci = +k_f * sc[6];
  J[491] += dqdci; // dwdot[H2O]/d[CH3O]
  J[492] -= dqdci; // dwdot[OH]/d[CH3O]
  J[504] -= dqdci; // dwdot[CH3O]/d[CH3O]
  J[506] += dqdci; // dwdot[CH2O]/d[CH3O]
  // d()/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[708] -= dqdT; // dwdot[OH]/dT
  J[720] -= dqdT; // dwdot[CH3O]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 123: CH2OH + CH3 => CH2O + CH4
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[12] * sc[19];
  k_f = 15000000;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[11] += q; // CH4
  wdot[12] -= q; // CH3
  wdot[19] -= q; // CH2OH
  wdot[20] += q; // CH2O
  // d()/d[CH3]
  dqdci = +k_f * sc[19];
  J[335] += dqdci; // dwdot[CH4]/d[CH3]
  J[336] -= dqdci; // dwdot[CH3]/d[CH3]
  J[343] -= dqdci; // dwdot[CH2OH]/d[CH3]
  J[344] += dqdci; // dwdot[CH2O]/d[CH3]
  // d()/d[CH2OH]
  dqdci = +k_f * sc[12];
  J[524] += dqdci; // dwdot[CH4]/d[CH2OH]
  J[525] -= dqdci; // dwdot[CH3]/d[CH2OH]
  J[532] -= dqdci; // dwdot[CH2OH]/d[CH2OH]
  J[533] += dqdci; // dwdot[CH2O]/d[CH2OH]
  // d()/dT
  J[713] += dqdT; // dwdot[CH4]/dT
  J[714] -= dqdT; // dwdot[CH3]/dT
  J[721] -= dqdT; // dwdot[CH2OH]/dT
  J[722] += dqdT; // dwdot[CH2O]/dT

  // reaction 124: CH + O <=> CHO+ + E
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[4] * sc[16];
  k_f = 1746000000000 * exp(-2.19 * tc[0] - (164.999209334044) * invT);
  dlnkfdT = -2.19 * invT + (164.999209334044) * invT2;
  // reverse
  phi_r = sc[22] * sc[25];
  Kc = exp(g_RT[4] + g_RT[16] - g_RT[22] - g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4] + h_RT[16]) + (h_RT[22] + h_RT[25]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[4] -= q;  // O
  wdot[16] -= q; // CH
  wdot[22] += q; // CHO+
  wdot[25] += q; // E
  // d()/d[O]
  dqdci = +k_f * sc[16];
  J[112] -= dqdci; // dwdot[O]/d[O]
  J[124] -= dqdci; // dwdot[CH]/d[O]
  J[130] += dqdci; // dwdot[CHO+]/d[O]
  J[133] += dqdci; // dwdot[E]/d[O]
  // d()/d[CH]
  dqdci = +k_f * sc[4];
  J[436] -= dqdci; // dwdot[O]/d[CH]
  J[448] -= dqdci; // dwdot[CH]/d[CH]
  J[454] += dqdci; // dwdot[CHO+]/d[CH]
  J[457] += dqdci; // dwdot[E]/d[CH]
  // d()/d[CHO+]
  dqdci = -k_r * sc[25];
  J[598] -= dqdci; // dwdot[O]/d[CHO+]
  J[610] -= dqdci; // dwdot[CH]/d[CHO+]
  J[616] += dqdci; // dwdot[CHO+]/d[CHO+]
  J[619] += dqdci; // dwdot[E]/d[CHO+]
  // d()/d[E]
  dqdci = -k_r * sc[22];
  J[679] -= dqdci; // dwdot[O]/d[E]
  J[691] -= dqdci; // dwdot[CH]/d[E]
  J[697] += dqdci; // dwdot[CHO+]/d[E]
  J[700] += dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[706] -= dqdT; // dwdot[O]/dT
  J[718] -= dqdT; // dwdot[CH]/dT
  J[724] += dqdT; // dwdot[CHO+]/dT
  J[727] += dqdT; // dwdot[E]/dT

  // reaction 125: CHO+ + H2O <=> CO + H3O+
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[5] * sc[22];
  k_f = 26080000000 * exp(-0.5 * tc[0]);
  dlnkfdT = -0.5 * invT;
  // reverse
  phi_r = sc[9] * sc[23];
  Kc = exp(g_RT[5] - g_RT[9] + g_RT[22] - g_RT[23]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[5] + h_RT[22]) + (h_RT[9] + h_RT[23]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[5] -= q;  // H2O
  wdot[9] += q;  // CO
  wdot[22] -= q; // CHO+
  wdot[23] += q; // H3O+
  // d()/d[H2O]
  dqdci = +k_f * sc[22];
  J[140] -= dqdci; // dwdot[H2O]/d[H2O]
  J[144] += dqdci; // dwdot[CO]/d[H2O]
  J[157] -= dqdci; // dwdot[CHO+]/d[H2O]
  J[158] += dqdci; // dwdot[H3O+]/d[H2O]
  // d()/d[CO]
  dqdci = -k_r * sc[23];
  J[248] -= dqdci; // dwdot[H2O]/d[CO]
  J[252] += dqdci; // dwdot[CO]/d[CO]
  J[265] -= dqdci; // dwdot[CHO+]/d[CO]
  J[266] += dqdci; // dwdot[H3O+]/d[CO]
  // d()/d[CHO+]
  dqdci = +k_f * sc[5];
  J[599] -= dqdci; // dwdot[H2O]/d[CHO+]
  J[603] += dqdci; // dwdot[CO]/d[CHO+]
  J[616] -= dqdci; // dwdot[CHO+]/d[CHO+]
  J[617] += dqdci; // dwdot[H3O+]/d[CHO+]
  // d()/d[H3O+]
  dqdci = -k_r * sc[9];
  J[626] -= dqdci; // dwdot[H2O]/d[H3O+]
  J[630] += dqdci; // dwdot[CO]/d[H3O+]
  J[643] -= dqdci; // dwdot[CHO+]/d[H3O+]
  J[644] += dqdci; // dwdot[H3O+]/d[H3O+]
  // d()/dT
  J[707] -= dqdT; // dwdot[H2O]/dT
  J[711] += dqdT; // dwdot[CO]/dT
  J[724] -= dqdT; // dwdot[CHO+]/dT
  J[725] += dqdT; // dwdot[H3O+]/dT

  // reaction 126: E + H3O+ <=> H + H2O
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[23] * sc[25];
  k_f = 739500000000 * exp(-0.5 * tc[0]);
  dlnkfdT = -0.5 * invT;
  // reverse
  phi_r = sc[2] * sc[5];
  Kc = exp(-g_RT[2] - g_RT[5] + g_RT[23] + g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[23] + h_RT[25]) + (h_RT[2] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += q;  // H
  wdot[5] += q;  // H2O
  wdot[23] -= q; // H3O+
  wdot[25] -= q; // E
  // d()/d[H]
  dqdci = -k_r * sc[5];
  J[56] += dqdci; // dwdot[H]/d[H]
  J[59] += dqdci; // dwdot[H2O]/d[H]
  J[77] -= dqdci; // dwdot[H3O+]/d[H]
  J[79] -= dqdci; // dwdot[E]/d[H]
  // d()/d[H2O]
  dqdci = -k_r * sc[2];
  J[137] += dqdci; // dwdot[H]/d[H2O]
  J[140] += dqdci; // dwdot[H2O]/d[H2O]
  J[158] -= dqdci; // dwdot[H3O+]/d[H2O]
  J[160] -= dqdci; // dwdot[E]/d[H2O]
  // d()/d[H3O+]
  dqdci = +k_f * sc[25];
  J[623] += dqdci; // dwdot[H]/d[H3O+]
  J[626] += dqdci; // dwdot[H2O]/d[H3O+]
  J[644] -= dqdci; // dwdot[H3O+]/d[H3O+]
  J[646] -= dqdci; // dwdot[E]/d[H3O+]
  // d()/d[E]
  dqdci = +k_f * sc[23];
  J[677] += dqdci; // dwdot[H]/d[E]
  J[680] += dqdci; // dwdot[H2O]/d[E]
  J[698] -= dqdci; // dwdot[H3O+]/d[E]
  J[700] -= dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[704] += dqdT; // dwdot[H]/dT
  J[707] += dqdT; // dwdot[H2O]/dT
  J[725] -= dqdT; // dwdot[H3O+]/dT
  J[727] -= dqdT; // dwdot[E]/dT

  // reaction 127: E + H3O+ <=> 2 H + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[23] * sc[25];
  k_f = 3181000000000 * exp(-0.5 * tc[0]);
  dlnkfdT = -0.5 * invT;
  // reverse
  phi_r = (sc[2] * sc[2]) * sc[6];
  Kc = refC * exp(-2.000000 * g_RT[2] - g_RT[6] + g_RT[23] + g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT =
    invT * (-(h_RT[23] + h_RT[25]) + (2.000000 * h_RT[2] + h_RT[6]) - 1.000000);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] += 2 * q; // H
  wdot[6] += q;     // OH
  wdot[23] -= q;    // H3O+
  wdot[25] -= q;    // E
  // d()/d[H]
  dqdci = -k_r * 2.000000 * sc[2] * sc[6];
  J[56] += 2 * dqdci; // dwdot[H]/d[H]
  J[60] += dqdci;     // dwdot[OH]/d[H]
  J[77] -= dqdci;     // dwdot[H3O+]/d[H]
  J[79] -= dqdci;     // dwdot[E]/d[H]
  // d()/d[OH]
  dqdci = -k_r * sc[2] * sc[2];
  J[164] += 2 * dqdci; // dwdot[H]/d[OH]
  J[168] += dqdci;     // dwdot[OH]/d[OH]
  J[185] -= dqdci;     // dwdot[H3O+]/d[OH]
  J[187] -= dqdci;     // dwdot[E]/d[OH]
  // d()/d[H3O+]
  dqdci = +k_f * sc[25];
  J[623] += 2 * dqdci; // dwdot[H]/d[H3O+]
  J[627] += dqdci;     // dwdot[OH]/d[H3O+]
  J[644] -= dqdci;     // dwdot[H3O+]/d[H3O+]
  J[646] -= dqdci;     // dwdot[E]/d[H3O+]
  // d()/d[E]
  dqdci = +k_f * sc[23];
  J[677] += 2 * dqdci; // dwdot[H]/d[E]
  J[681] += dqdci;     // dwdot[OH]/d[E]
  J[698] -= dqdci;     // dwdot[H3O+]/d[E]
  J[700] -= dqdci;     // dwdot[E]/d[E]
  // d()/dT
  J[704] += 2 * dqdT; // dwdot[H]/dT
  J[708] += dqdT;     // dwdot[OH]/dT
  J[725] -= dqdT;     // dwdot[H3O+]/dT
  J[727] -= dqdT;     // dwdot[E]/dT

  // reaction 128: E + H3O+ <=> H2 + OH
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[23] * sc[25];
  k_f = 560100000000 * exp(-0.5 * tc[0]);
  dlnkfdT = -0.5 * invT;
  // reverse
  phi_r = sc[1] * sc[6];
  Kc = exp(-g_RT[1] - g_RT[6] + g_RT[23] + g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[23] + h_RT[25]) + (h_RT[1] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] += q;  // H2
  wdot[6] += q;  // OH
  wdot[23] -= q; // H3O+
  wdot[25] -= q; // E
  // d()/d[H2]
  dqdci = -k_r * sc[6];
  J[28] += dqdci; // dwdot[H2]/d[H2]
  J[33] += dqdci; // dwdot[OH]/d[H2]
  J[50] -= dqdci; // dwdot[H3O+]/d[H2]
  J[52] -= dqdci; // dwdot[E]/d[H2]
  // d()/d[OH]
  dqdci = -k_r * sc[1];
  J[163] += dqdci; // dwdot[H2]/d[OH]
  J[168] += dqdci; // dwdot[OH]/d[OH]
  J[185] -= dqdci; // dwdot[H3O+]/d[OH]
  J[187] -= dqdci; // dwdot[E]/d[OH]
  // d()/d[H3O+]
  dqdci = +k_f * sc[25];
  J[622] += dqdci; // dwdot[H2]/d[H3O+]
  J[627] += dqdci; // dwdot[OH]/d[H3O+]
  J[644] -= dqdci; // dwdot[H3O+]/d[H3O+]
  J[646] -= dqdci; // dwdot[E]/d[H3O+]
  // d()/d[E]
  dqdci = +k_f * sc[23];
  J[676] += dqdci; // dwdot[H2]/d[E]
  J[681] += dqdci; // dwdot[OH]/d[E]
  J[698] -= dqdci; // dwdot[H3O+]/d[E]
  J[700] -= dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[703] += dqdT; // dwdot[H2]/dT
  J[708] += dqdT; // dwdot[OH]/dT
  J[725] -= dqdT; // dwdot[H3O+]/dT
  J[727] -= dqdT; // dwdot[E]/dT

  // reaction 129: H2 + O2- <=> E + H2O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[1] * sc[24];
  k_f = 602000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[7] * sc[25];
  Kc = exp(g_RT[1] - g_RT[7] + g_RT[24] - g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[24]) + (h_RT[7] + h_RT[25]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[1] -= q;  // H2
  wdot[7] += q;  // H2O2
  wdot[24] -= q; // O2-
  wdot[25] += q; // E
  // d()/d[H2]
  dqdci = +k_f * sc[24];
  J[28] -= dqdci; // dwdot[H2]/d[H2]
  J[34] += dqdci; // dwdot[H2O2]/d[H2]
  J[51] -= dqdci; // dwdot[O2-]/d[H2]
  J[52] += dqdci; // dwdot[E]/d[H2]
  // d()/d[H2O2]
  dqdci = -k_r * sc[25];
  J[190] -= dqdci; // dwdot[H2]/d[H2O2]
  J[196] += dqdci; // dwdot[H2O2]/d[H2O2]
  J[213] -= dqdci; // dwdot[O2-]/d[H2O2]
  J[214] += dqdci; // dwdot[E]/d[H2O2]
  // d()/d[O2-]
  dqdci = +k_f * sc[1];
  J[649] -= dqdci; // dwdot[H2]/d[O2-]
  J[655] += dqdci; // dwdot[H2O2]/d[O2-]
  J[672] -= dqdci; // dwdot[O2-]/d[O2-]
  J[673] += dqdci; // dwdot[E]/d[O2-]
  // d()/d[E]
  dqdci = -k_r * sc[7];
  J[676] -= dqdci; // dwdot[H2]/d[E]
  J[682] += dqdci; // dwdot[H2O2]/d[E]
  J[699] -= dqdci; // dwdot[O2-]/d[E]
  J[700] += dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[703] -= dqdT; // dwdot[H2]/dT
  J[709] += dqdT; // dwdot[H2O2]/dT
  J[726] -= dqdT; // dwdot[O2-]/dT
  J[727] += dqdT; // dwdot[E]/dT

  // reaction 130: H + O2- <=> E + HO2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[24];
  k_f = 723000000;
  dlnkfdT = 0.0;
  // reverse
  phi_r = sc[8] * sc[25];
  Kc = exp(g_RT[2] - g_RT[8] + g_RT[24] - g_RT[25]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[24]) + (h_RT[8] + h_RT[25]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  wdot[2] -= q;  // H
  wdot[8] += q;  // HO2
  wdot[24] -= q; // O2-
  wdot[25] += q; // E
  // d()/d[H]
  dqdci = +k_f * sc[24];
  J[56] -= dqdci; // dwdot[H]/d[H]
  J[62] += dqdci; // dwdot[HO2]/d[H]
  J[78] -= dqdci; // dwdot[O2-]/d[H]
  J[79] += dqdci; // dwdot[E]/d[H]
  // d()/d[HO2]
  dqdci = -k_r * sc[25];
  J[218] -= dqdci; // dwdot[H]/d[HO2]
  J[224] += dqdci; // dwdot[HO2]/d[HO2]
  J[240] -= dqdci; // dwdot[O2-]/d[HO2]
  J[241] += dqdci; // dwdot[E]/d[HO2]
  // d()/d[O2-]
  dqdci = +k_f * sc[2];
  J[650] -= dqdci; // dwdot[H]/d[O2-]
  J[656] += dqdci; // dwdot[HO2]/d[O2-]
  J[672] -= dqdci; // dwdot[O2-]/d[O2-]
  J[673] += dqdci; // dwdot[E]/d[O2-]
  // d()/d[E]
  dqdci = -k_r * sc[8];
  J[677] -= dqdci; // dwdot[H]/d[E]
  J[683] += dqdci; // dwdot[HO2]/d[E]
  J[699] -= dqdci; // dwdot[O2-]/d[E]
  J[700] += dqdci; // dwdot[E]/d[E]
  // d()/dT
  J[704] -= dqdT; // dwdot[H]/dT
  J[710] += dqdT; // dwdot[HO2]/dT
  J[726] -= dqdT; // dwdot[O2-]/dT
  J[727] += dqdT; // dwdot[E]/dT

  amrex::Real c_R[26], dcRdT[26], e_RT[26];
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
  for (int k = 0; k < 26; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[702 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 26; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 26; ++m) {
      dehmixdc += eh_RT[m] * J[k * 27 + m];
    }
    J[k * 27 + 26] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[728] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
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
  amrex::Real c[26];

  for (int k = 0; k < 26; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian(J, c, *Tp, *consP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 26; k++) {
    J[702 + k] *= 1.e-6;
    J[k * 27 + 26] *= 1.e6;
  }
}

// Transport function declarations

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENIMC(int* LENIMC)
{
  *LENIMC = 107;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENRMC(int* LENRMC)
{
  *LENRMC = 13858;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNO(int* NO)
{
  *NO = 4;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKK(int* KK)
{
  *KK = 26;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLITE(int* NLITE)
{
  *NLITE = 3;
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
  WT[1] = 2.01600000E+00;
  WT[2] = 1.00800000E+00;
  WT[3] = 3.19980000E+01;
  WT[4] = 1.59990000E+01;
  WT[5] = 1.80150000E+01;
  WT[6] = 1.70070000E+01;
  WT[7] = 3.40140000E+01;
  WT[8] = 3.30060000E+01;
  WT[9] = 2.80100000E+01;
  WT[10] = 4.40090000E+01;
  WT[11] = 1.60430000E+01;
  WT[12] = 1.50350000E+01;
  WT[13] = 1.40270000E+01;
  WT[14] = 1.40270000E+01;
  WT[15] = 1.20110000E+01;
  WT[16] = 1.30190000E+01;
  WT[17] = 3.20420000E+01;
  WT[18] = 3.10340000E+01;
  WT[19] = 3.10340000E+01;
  WT[20] = 3.00260000E+01;
  WT[21] = 2.90180000E+01;
  WT[22] = 2.90174514E+01;
  WT[23] = 1.90224514E+01;
  WT[24] = 3.19985486E+01;
  WT[25] = 5.48579909E-04;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetEPS(amrex::Real* EPS)
{
  EPS[0] = 9.75300000E+01;
  EPS[1] = 3.80000000E+01;
  EPS[2] = 1.45000000E+02;
  EPS[3] = 1.07400000E+02;
  EPS[4] = 8.00000000E+01;
  EPS[5] = 5.72400000E+02;
  EPS[6] = 8.00000000E+01;
  EPS[7] = 1.07400000E+02;
  EPS[8] = 1.07400000E+02;
  EPS[9] = 9.81000000E+01;
  EPS[10] = 2.44000000E+02;
  EPS[11] = 1.41400000E+02;
  EPS[12] = 1.44000000E+02;
  EPS[13] = 1.44000000E+02;
  EPS[14] = 1.44000000E+02;
  EPS[15] = 7.14000000E+01;
  EPS[16] = 8.00000000E+01;
  EPS[17] = 4.81800000E+02;
  EPS[18] = 4.17000000E+02;
  EPS[19] = 4.17000000E+02;
  EPS[20] = 4.98000000E+02;
  EPS[21] = 4.98000000E+02;
  EPS[22] = 8.00000000E+01;
  EPS[23] = 1.06200000E+02;
  EPS[24] = 1.36500000E+02;
  EPS[25] = 8.50000000E+02;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetSIG(amrex::Real* SIG)
{
  SIG[0] = 3.62100000E+00;
  SIG[1] = 2.92000000E+00;
  SIG[2] = 2.05000000E+00;
  SIG[3] = 3.45800000E+00;
  SIG[4] = 2.75000000E+00;
  SIG[5] = 2.60500000E+00;
  SIG[6] = 2.75000000E+00;
  SIG[7] = 3.45800000E+00;
  SIG[8] = 3.45800000E+00;
  SIG[9] = 3.65000000E+00;
  SIG[10] = 3.76300000E+00;
  SIG[11] = 3.74600000E+00;
  SIG[12] = 3.80000000E+00;
  SIG[13] = 3.80000000E+00;
  SIG[14] = 3.80000000E+00;
  SIG[15] = 3.29800000E+00;
  SIG[16] = 2.75000000E+00;
  SIG[17] = 3.62600000E+00;
  SIG[18] = 3.69000000E+00;
  SIG[19] = 3.69000000E+00;
  SIG[20] = 3.59000000E+00;
  SIG[21] = 3.59000000E+00;
  SIG[22] = 2.75000000E+00;
  SIG[23] = 3.15000000E+00;
  SIG[24] = 3.33000000E+00;
  SIG[25] = 4.25000000E+02;
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
  DIP[5] = 1.84400053E+00;
  DIP[6] = 0.00000000E+00;
  DIP[7] = 0.00000000E+00;
  DIP[8] = 0.00000000E+00;
  DIP[9] = 0.00000000E+00;
  DIP[10] = 0.00000000E+00;
  DIP[11] = 0.00000000E+00;
  DIP[12] = 0.00000000E+00;
  DIP[13] = 0.00000000E+00;
  DIP[14] = 0.00000000E+00;
  DIP[15] = 0.00000000E+00;
  DIP[16] = 0.00000000E+00;
  DIP[17] = 0.00000000E+00;
  DIP[18] = 1.70000049E+00;
  DIP[19] = 1.70000049E+00;
  DIP[20] = 0.00000000E+00;
  DIP[21] = 0.00000000E+00;
  DIP[22] = 0.00000000E+00;
  DIP[23] = 0.00000000E+00;
  DIP[24] = 0.00000000E+00;
  DIP[25] = 0.00000000E+00;
}

// the polarizability in cubic Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPOL(amrex::Real* POL)
{
  POL[0] = 1.76000000E+00;
  POL[1] = 7.90000000E-01;
  POL[2] = 0.00000000E+00;
  POL[3] = 1.60000000E+00;
  POL[4] = 0.00000000E+00;
  POL[5] = 0.00000000E+00;
  POL[6] = 0.00000000E+00;
  POL[7] = 0.00000000E+00;
  POL[8] = 0.00000000E+00;
  POL[9] = 1.95000000E+00;
  POL[10] = 2.65000000E+00;
  POL[11] = 2.60000000E+00;
  POL[12] = 0.00000000E+00;
  POL[13] = 0.00000000E+00;
  POL[14] = 0.00000000E+00;
  POL[15] = 0.00000000E+00;
  POL[16] = 0.00000000E+00;
  POL[17] = 0.00000000E+00;
  POL[18] = 0.00000000E+00;
  POL[19] = 0.00000000E+00;
  POL[20] = 0.00000000E+00;
  POL[21] = 0.00000000E+00;
  POL[22] = 1.34100000E+00;
  POL[23] = 9.64000000E-01;
  POL[24] = 1.58100000E+00;
  POL[25] = 0.00000000E+00;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetZROT(amrex::Real* ZROT)
{
  ZROT[0] = 4.00000000E+00;
  ZROT[1] = 2.80000000E+02;
  ZROT[2] = 0.00000000E+00;
  ZROT[3] = 3.80000000E+00;
  ZROT[4] = 0.00000000E+00;
  ZROT[5] = 4.00000000E+00;
  ZROT[6] = 0.00000000E+00;
  ZROT[7] = 3.80000000E+00;
  ZROT[8] = 1.00000000E+00;
  ZROT[9] = 1.80000000E+00;
  ZROT[10] = 2.10000000E+00;
  ZROT[11] = 1.30000000E+01;
  ZROT[12] = 0.00000000E+00;
  ZROT[13] = 0.00000000E+00;
  ZROT[14] = 0.00000000E+00;
  ZROT[15] = 0.00000000E+00;
  ZROT[16] = 0.00000000E+00;
  ZROT[17] = 1.00000000E+00;
  ZROT[18] = 2.00000000E+00;
  ZROT[19] = 2.00000000E+00;
  ZROT[20] = 2.00000000E+00;
  ZROT[21] = 0.00000000E+00;
  ZROT[22] = 0.00000000E+00;
  ZROT[23] = 1.00000000E+01;
  ZROT[24] = 3.80000000E+00;
  ZROT[25] = 1.00000000E+00;
}

// 0: monoatomic, 1: linear, 2: nonlinear
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLIN(int* NLIN)
{
  NLIN[0] = 1;
  NLIN[1] = 1;
  NLIN[2] = 0;
  NLIN[3] = 1;
  NLIN[4] = 0;
  NLIN[5] = 2;
  NLIN[6] = 1;
  NLIN[7] = 2;
  NLIN[8] = 2;
  NLIN[9] = 1;
  NLIN[10] = 1;
  NLIN[11] = 2;
  NLIN[12] = 1;
  NLIN[13] = 1;
  NLIN[14] = 1;
  NLIN[15] = 0;
  NLIN[16] = 1;
  NLIN[17] = 2;
  NLIN[18] = 2;
  NLIN[19] = 2;
  NLIN[20] = 2;
  NLIN[21] = 2;
  NLIN[22] = 1;
  NLIN[23] = 2;
  NLIN[24] = 1;
  NLIN[25] = 0;
}

// Poly fits for the viscosities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFETA(amrex::Real* COFETA)
{
  COFETA[0] = -1.62526676E+01;
  COFETA[1] = 2.24839597E+00;
  COFETA[2] = -2.13428438E-01;
  COFETA[3] = 9.46192413E-03;
  COFETA[4] = -1.37549291E+01;
  COFETA[5] = 9.65530587E-01;
  COFETA[6] = -4.45720114E-02;
  COFETA[7] = 2.05871810E-03;
  COFETA[8] = -1.98744352E+01;
  COFETA[9] = 3.41660514E+00;
  COFETA[10] = -3.63206306E-01;
  COFETA[11] = 1.58671021E-02;
  COFETA[12] = -1.68118998E+01;
  COFETA[13] = 2.52362554E+00;
  COFETA[14] = -2.49309128E-01;
  COFETA[15] = 1.10211025E-02;
  COFETA[16] = -1.48001711E+01;
  COFETA[17] = 1.79491990E+00;
  COFETA[18] = -1.54008440E-01;
  COFETA[19] = 6.86719439E-03;
  COFETA[20] = -1.17770970E+01;
  COFETA[21] = -8.26746102E-01;
  COFETA[22] = 3.39009582E-01;
  COFETA[23] = -2.00674559E-02;
  COFETA[24] = -1.47696217E+01;
  COFETA[25] = 1.79491990E+00;
  COFETA[26] = -1.54008440E-01;
  COFETA[27] = 6.86719439E-03;
  COFETA[28] = -1.67813504E+01;
  COFETA[29] = 2.52362554E+00;
  COFETA[30] = -2.49309128E-01;
  COFETA[31] = 1.10211025E-02;
  COFETA[32] = -1.67963918E+01;
  COFETA[33] = 2.52362554E+00;
  COFETA[34] = -2.49309128E-01;
  COFETA[35] = 1.10211025E-02;
  COFETA[36] = -1.63031343E+01;
  COFETA[37] = 2.26143219E+00;
  COFETA[38] = -2.15114671E-01;
  COFETA[39] = 9.53461976E-03;
  COFETA[40] = -2.36749638E+01;
  COFETA[41] = 4.99775518E+00;
  COFETA[42] = -5.52687718E-01;
  COFETA[43] = 2.34353338E-02;
  COFETA[44] = -1.95453436E+01;
  COFETA[45] = 3.36385478E+00;
  COFETA[46] = -3.56948469E-01;
  COFETA[47] = 1.56210922E-02;
  COFETA[48] = -1.97154229E+01;
  COFETA[49] = 3.40193237E+00;
  COFETA[50] = -3.61470464E-01;
  COFETA[51] = 1.57990760E-02;
  COFETA[52] = -1.97501213E+01;
  COFETA[53] = 3.40193237E+00;
  COFETA[54] = -3.61470464E-01;
  COFETA[55] = 1.57990760E-02;
  COFETA[56] = -1.97501213E+01;
  COFETA[57] = 3.40193237E+00;
  COFETA[58] = -3.61470464E-01;
  COFETA[59] = 1.57990760E-02;
  COFETA[60] = -1.47250813E+01;
  COFETA[61] = 1.57011544E+00;
  COFETA[62] = -1.24217373E-01;
  COFETA[63] = 5.55161238E-03;
  COFETA[64] = -1.49032293E+01;
  COFETA[65] = 1.79491990E+00;
  COFETA[66] = -1.54008440E-01;
  COFETA[67] = 6.86719439E-03;
  COFETA[68] = -2.17058520E+01;
  COFETA[69] = 3.55238511E+00;
  COFETA[70] = -2.93340526E-01;
  COFETA[71] = 9.37852541E-03;
  COFETA[72] = -2.09501615E+01;
  COFETA[73] = 3.29835026E+00;
  COFETA[74] = -2.67151401E-01;
  COFETA[75] = 8.53547645E-03;
  COFETA[76] = -2.09501615E+01;
  COFETA[77] = 3.29835026E+00;
  COFETA[78] = -2.67151401E-01;
  COFETA[79] = 8.53547645E-03;
  COFETA[80] = -2.11136150E+01;
  COFETA[81] = 3.26961843E+00;
  COFETA[82] = -2.51355092E-01;
  COFETA[83] = 7.35605058E-03;
  COFETA[84] = -2.11306886E+01;
  COFETA[85] = 3.26961843E+00;
  COFETA[86] = -2.51355092E-01;
  COFETA[87] = 7.35605058E-03;
  COFETA[88] = -1.45024855E+01;
  COFETA[89] = 1.79491990E+00;
  COFETA[90] = -1.54008440E-01;
  COFETA[91] = 6.86719439E-03;
  COFETA[92] = -1.67955362E+01;
  COFETA[93] = 2.48907364E+00;
  COFETA[94] = -2.44793974E-01;
  COFETA[95] = 1.08243607E-02;
  COFETA[96] = -1.87177043E+01;
  COFETA[97] = 3.27402596E+00;
  COFETA[98] = -3.45827972E-01;
  COFETA[99] = 1.51622680E-02;
  COFETA[100] = -2.18499525E+01;
  COFETA[101] = -3.03985219E+00;
  COFETA[102] = 6.46768575E-01;
  COFETA[103] = -3.43986573E-02;
}

// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFLAM(amrex::Real* COFLAM)
{
  COFLAM[0] = 1.15135759E+01;
  COFLAM[1] = -2.89718927E+00;
  COFLAM[2] = 5.52227164E-01;
  COFLAM[3] = -2.73618079E-02;
  COFLAM[4] = 1.05575645E+01;
  COFLAM[5] = -1.08525230E+00;
  COFLAM[6] = 2.11035028E-01;
  COFLAM[7] = -7.45871436E-03;
  COFLAM[8] = -3.24554543E-01;
  COFLAM[9] = 3.41660514E+00;
  COFLAM[10] = -3.63206306E-01;
  COFLAM[11] = 1.58671021E-02;
  COFLAM[12] = -2.39605522E+00;
  COFLAM[13] = 3.11185492E+00;
  COFLAM[14] = -3.06053822E-01;
  COFLAM[15] = 1.33590959E-02;
  COFLAM[16] = 1.98515155E+00;
  COFLAM[17] = 1.79491990E+00;
  COFLAM[18] = -1.54008440E-01;
  COFLAM[19] = 6.86719439E-03;
  COFLAM[20] = 2.10125974E+01;
  COFLAM[21] = -7.92638744E+00;
  COFLAM[22] = 1.37729291E+00;
  COFLAM[23] = -6.86156427E-02;
  COFLAM[24] = 1.50443438E+01;
  COFLAM[25] = -3.64813621E+00;
  COFLAM[26] = 5.94792348E-01;
  COFLAM[27] = -2.63489305E-02;
  COFLAM[28] = 1.90147776E+00;
  COFLAM[29] = 9.16567599E-01;
  COFLAM[30] = 6.93878492E-02;
  COFLAM[31] = -6.39769965E-03;
  COFLAM[32] = 2.91063106E+00;
  COFLAM[33] = 5.78338476E-01;
  COFLAM[34] = 9.02531383E-02;
  COFLAM[35] = -6.19296183E-03;
  COFLAM[36] = 1.00332074E+01;
  COFLAM[37] = -2.32812074E+00;
  COFLAM[38] = 4.79111919E-01;
  COFLAM[39] = -2.42630313E-02;
  COFLAM[40] = -1.19710226E+01;
  COFLAM[41] = 6.16510101E+00;
  COFLAM[42] = -6.12547921E-01;
  COFLAM[43] = 2.26087237E-02;
  COFLAM[44] = 6.90173358E+00;
  COFLAM[45] = -2.03450657E+00;
  COFLAM[46] = 5.89314746E-01;
  COFLAM[47] = -3.40238924E-02;
  COFLAM[48] = 1.24262108E+01;
  COFLAM[49] = -3.92863760E+00;
  COFLAM[50] = 7.99986093E-01;
  COFLAM[51] = -4.23406929E-02;
  COFLAM[52] = 1.26992050E+01;
  COFLAM[53] = -3.59706429E+00;
  COFLAM[54] = 6.88908305E-01;
  COFLAM[55] = -3.47879197E-02;
  COFLAM[56] = 1.34129531E+01;
  COFLAM[57] = -3.93997381E+00;
  COFLAM[58] = 7.35902066E-01;
  COFLAM[59] = -3.65948332E-02;
  COFLAM[60] = 2.34694471E+00;
  COFLAM[61] = 1.57011544E+00;
  COFLAM[62] = -1.24217373E-01;
  COFLAM[63] = 5.55161238E-03;
  COFLAM[64] = 2.19696972E+01;
  COFLAM[65] = -6.76971206E+00;
  COFLAM[66] = 1.06194842E+00;
  COFLAM[67] = -4.90083482E-02;
  COFLAM[68] = -8.91606915E+00;
  COFLAM[69] = 3.61762165E+00;
  COFLAM[70] = -8.92916386E-02;
  COFLAM[71] = -7.46114102E-03;
  COFLAM[72] = -1.33068537E+01;
  COFLAM[73] = 5.77901128E+00;
  COFLAM[74] = -4.27810572E-01;
  COFLAM[75] = 9.36447961E-03;
  COFLAM[76] = -6.42835893E+00;
  COFLAM[77] = 3.26905843E+00;
  COFLAM[78] = -1.28184954E-01;
  COFLAM[79] = -2.45699613E-03;
  COFLAM[80] = 7.06271985E-01;
  COFLAM[81] = -2.71898947E-01;
  COFLAM[82] = 4.21512856E-01;
  COFLAM[83] = -3.00411809E-02;
  COFLAM[84] = -1.49162912E+00;
  COFLAM[85] = 1.41832455E+00;
  COFLAM[86] = 6.93252840E-02;
  COFLAM[87] = -8.61179989E-03;
  COFLAM[88] = 6.28525712E+00;
  COFLAM[89] = -8.51059745E-01;
  COFLAM[90] = 3.29554354E-01;
  COFLAM[91] = -1.93574144E-02;
  COFLAM[92] = 1.18912374E+01;
  COFLAM[93] = -3.44418514E+00;
  COFLAM[94] = 7.12558841E-01;
  COFLAM[95] = -3.73179377E-02;
  COFLAM[96] = -6.15617763E+00;
  COFLAM[97] = 4.62492578E+00;
  COFLAM[98] = -4.96911868E-01;
  COFLAM[99] = 2.09308888E-02;
  COFLAM[100] = 5.21607397E+00;
  COFLAM[101] = -3.03985219E+00;
  COFLAM[102] = 6.46768575E-01;
  COFLAM[103] = -3.43986573E-02;
}

// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFD(amrex::Real* COFD)
{
  COFD[0] = -1.47639411E+01;
  COFD[1] = 3.15955654E+00;
  COFD[2] = -1.97590757E-01;
  COFD[3] = 8.64692156E-03;
  COFD[4] = -1.13253618E+01;
  COFD[5] = 2.31195095E+00;
  COFD[6] = -8.63988037E-02;
  COFD[7] = 3.77573452E-03;
  COFD[8] = -1.40298991E+01;
  COFD[9] = 3.55837688E+00;
  COFD[10] = -2.47785790E-01;
  COFD[11] = 1.07555332E-02;
  COFD[12] = -1.50096252E+01;
  COFD[13] = 3.25515933E+00;
  COFD[14] = -2.09710110E-01;
  COFD[15] = 9.15941830E-03;
  COFD[16] = -1.38756380E+01;
  COFD[17] = 2.98558426E+00;
  COFD[18] = -1.75507216E-01;
  COFD[19] = 7.71173691E-03;
  COFD[20] = -2.08123335E+01;
  COFD[21] = 5.42470187E+00;
  COFD[22] = -4.69700452E-01;
  COFD[23] = 1.95706916E-02;
  COFD[24] = -1.38948654E+01;
  COFD[25] = 2.98558426E+00;
  COFD[26] = -1.75507216E-01;
  COFD[27] = 7.71173691E-03;
  COFD[28] = -1.50236540E+01;
  COFD[29] = 3.25515933E+00;
  COFD[30] = -2.09710110E-01;
  COFD[31] = 9.15941830E-03;
  COFD[32] = -1.50168047E+01;
  COFD[33] = 3.25515933E+00;
  COFD[34] = -2.09710110E-01;
  COFD[35] = 9.15941830E-03;
  COFD[36] = -1.47850505E+01;
  COFD[37] = 3.16433919E+00;
  COFD[38] = -1.98191564E-01;
  COFD[39] = 8.67209742E-03;
  COFD[40] = -1.77350630E+01;
  COFD[41] = 4.19328271E+00;
  COFD[42] = -3.26911461E-01;
  COFD[43] = 1.40520357E-02;
  COFD[44] = -1.56756124E+01;
  COFD[45] = 3.54035770E+00;
  COFD[46] = -2.45653442E-01;
  COFD[47] = 1.06717969E-02;
  COFD[48] = -1.57070757E+01;
  COFD[49] = 3.55327776E+00;
  COFD[50] = -2.47179250E-01;
  COFD[51] = 1.07315653E-02;
  COFD[52] = -1.56842242E+01;
  COFD[53] = 3.55327776E+00;
  COFD[54] = -2.47179250E-01;
  COFD[55] = 1.07315653E-02;
  COFD[56] = -1.56842242E+01;
  COFD[57] = 3.55327776E+00;
  COFD[58] = -2.47179250E-01;
  COFD[59] = 1.07315653E-02;
  COFD[60] = -1.36591073E+01;
  COFD[61] = 2.87846903E+00;
  COFD[62] = -1.61681103E-01;
  COFD[63] = 7.11635115E-03;
  COFD[64] = -1.38076340E+01;
  COFD[65] = 2.98558426E+00;
  COFD[66] = -1.75507216E-01;
  COFD[67] = 7.71173691E-03;
  COFD[68] = -2.01091568E+01;
  COFD[69] = 5.06508846E+00;
  COFD[70] = -4.31649695E-01;
  COFD[71] = 1.82598046E-02;
  COFD[72] = -1.98937690E+01;
  COFD[73] = 4.98725744E+00;
  COFD[74] = -4.22747564E-01;
  COFD[75] = 1.79211872E-02;
  COFD[76] = -1.98937690E+01;
  COFD[77] = 4.98725744E+00;
  COFD[78] = -4.22747564E-01;
  COFD[79] = 1.79211872E-02;
  COFD[80] = -2.02136057E+01;
  COFD[81] = 5.10993120E+00;
  COFD[82] = -4.36931630E-01;
  COFD[83] = 1.84677592E-02;
  COFD[84] = -2.02052920E+01;
  COFD[85] = 5.10993120E+00;
  COFD[86] = -4.36931630E-01;
  COFD[87] = 1.84677592E-02;
  COFD[88] = -1.40437647E+01;
  COFD[89] = 2.98558426E+00;
  COFD[90] = -1.75507216E-01;
  COFD[91] = 7.71173691E-03;
  COFD[92] = -1.47449891E+01;
  COFD[93] = 3.24099818E+00;
  COFD[94] = -2.07879868E-01;
  COFD[95] = 9.08046161E-03;
  COFD[96] = -1.56749318E+01;
  COFD[97] = 3.51447210E+00;
  COFD[98] = -2.42579007E-01;
  COFD[99] = 1.05506318E-02;
  COFD[100] = -2.46746425E+01;
  COFD[101] = 5.51584151E+00;
  COFD[102] = -4.76108968E-01;
  COFD[103] = 1.96365975E-02;
  COFD[104] = -1.13253618E+01;
  COFD[105] = 2.31195095E+00;
  COFD[106] = -8.63988037E-02;
  COFD[107] = 3.77573452E-03;
  COFD[108] = -1.02395385E+01;
  COFD[109] = 2.15403244E+00;
  COFD[110] = -6.97480266E-02;
  COFD[111] = 3.23666871E-03;
  COFD[112] = -1.11808845E+01;
  COFD[113] = 2.66936727E+00;
  COFD[114] = -1.34411514E-01;
  COFD[115] = 5.92957488E-03;
  COFD[116] = -1.15797897E+01;
  COFD[117] = 2.43235504E+00;
  COFD[118] = -1.02890179E-01;
  COFD[119] = 4.52903603E-03;
  COFD[120] = -1.06250314E+01;
  COFD[121] = 2.15849701E+00;
  COFD[122] = -6.53886401E-02;
  COFD[123] = 2.81453370E-03;
  COFD[124] = -1.68759082E+01;
  COFD[125] = 4.49460364E+00;
  COFD[126] = -3.64766207E-01;
  COFD[127] = 1.56457185E-02;
  COFD[128] = -1.06283588E+01;
  COFD[129] = 2.15849701E+00;
  COFD[130] = -6.53886401E-02;
  COFD[131] = 2.81453370E-03;
  COFD[132] = -1.15815492E+01;
  COFD[133] = 2.43235504E+00;
  COFD[134] = -1.02890179E-01;
  COFD[135] = 4.52903603E-03;
  COFD[136] = -1.15806955E+01;
  COFD[137] = 2.43235504E+00;
  COFD[138] = -1.02890179E-01;
  COFD[139] = 4.52903603E-03;
  COFD[140] = -1.13541221E+01;
  COFD[141] = 2.31999438E+00;
  COFD[142] = -8.75064804E-02;
  COFD[143] = 3.82656365E-03;
  COFD[144] = -1.35545391E+01;
  COFD[145] = 3.13878730E+00;
  COFD[146] = -1.94980335E-01;
  COFD[147] = 8.53744486E-03;
  COFD[148] = -1.22114409E+01;
  COFD[149] = 2.65215786E+00;
  COFD[150] = -1.32206766E-01;
  COFD[151] = 5.83532580E-03;
  COFD[152] = -1.22583533E+01;
  COFD[153] = 2.66488613E+00;
  COFD[154] = -1.33839930E-01;
  COFD[155] = 5.90524087E-03;
  COFD[156] = -1.22541230E+01;
  COFD[157] = 2.66488613E+00;
  COFD[158] = -1.33839930E-01;
  COFD[159] = 5.90524087E-03;
  COFD[160] = -1.22541230E+01;
  COFD[161] = 2.66488613E+00;
  COFD[162] = -1.33839930E-01;
  COFD[163] = 5.90524087E-03;
  COFD[164] = -1.06186344E+01;
  COFD[165] = 2.09311029E+00;
  COFD[166] = -5.66715295E-02;
  COFD[167] = 2.42607686E-03;
  COFD[168] = -1.06123852E+01;
  COFD[169] = 2.15849701E+00;
  COFD[170] = -6.53886401E-02;
  COFD[171] = 2.81453370E-03;
  COFD[172] = -1.56250770E+01;
  COFD[173] = 3.92232623E+00;
  COFD[174] = -2.94183188E-01;
  COFD[175] = 1.27325506E-02;
  COFD[176] = -1.55102782E+01;
  COFD[177] = 3.86994367E+00;
  COFD[178] = -2.87570272E-01;
  COFD[179] = 1.24533966E-02;
  COFD[180] = -1.55102782E+01;
  COFD[181] = 3.86994367E+00;
  COFD[182] = -2.87570272E-01;
  COFD[183] = 1.24533966E-02;
  COFD[184] = -1.57172268E+01;
  COFD[185] = 3.96062263E+00;
  COFD[186] = -2.98964970E-01;
  COFD[187] = 1.29322565E-02;
  COFD[188] = -1.57161352E+01;
  COFD[189] = 3.96062263E+00;
  COFD[190] = -2.98964970E-01;
  COFD[191] = 1.29322565E-02;
  COFD[192] = -1.06507866E+01;
  COFD[193] = 2.15849701E+00;
  COFD[194] = -6.53886401E-02;
  COFD[195] = 2.81453370E-03;
  COFD[196] = -1.14322126E+01;
  COFD[197] = 2.42102995E+00;
  COFD[198] = -1.01354401E-01;
  COFD[199] = 4.45948069E-03;
  COFD[200] = -1.20579435E+01;
  COFD[201] = 2.63303536E+00;
  COFD[202] = -1.29792632E-01;
  COFD[203] = 5.73363738E-03;
  COFD[204] = -2.18596191E+01;
  COFD[205] = 4.62023553E+00;
  COFD[206] = -3.79733248E-01;
  COFD[207] = 1.62407913E-02;
  COFD[208] = -1.40298991E+01;
  COFD[209] = 3.55837688E+00;
  COFD[210] = -2.47785790E-01;
  COFD[211] = 1.07555332E-02;
  COFD[212] = -1.11808845E+01;
  COFD[213] = 2.66936727E+00;
  COFD[214] = -1.34411514E-01;
  COFD[215] = 5.92957488E-03;
  COFD[216] = -1.43693219E+01;
  COFD[217] = 4.03992999E+00;
  COFD[218] = -3.08044800E-01;
  COFD[219] = 1.32757775E-02;
  COFD[220] = -1.43713019E+01;
  COFD[221] = 3.70920439E+00;
  COFD[222] = -2.67274113E-01;
  COFD[223] = 1.15967481E-02;
  COFD[224] = -1.31860264E+01;
  COFD[225] = 3.38003453E+00;
  COFD[226] = -2.25783856E-01;
  COFD[227] = 9.85028660E-03;
  COFD[228] = -1.93611201E+01;
  COFD[229] = 5.51579726E+00;
  COFD[230] = -4.76061961E-01;
  COFD[231] = 1.96329391E-02;
  COFD[232] = -1.31877859E+01;
  COFD[233] = 3.38003453E+00;
  COFD[234] = -2.25783856E-01;
  COFD[235] = 9.85028660E-03;
  COFD[236] = -1.43722077E+01;
  COFD[237] = 3.70920439E+00;
  COFD[238] = -2.67274113E-01;
  COFD[239] = 1.15967481E-02;
  COFD[240] = -1.43717684E+01;
  COFD[241] = 3.70920439E+00;
  COFD[242] = -2.67274113E-01;
  COFD[243] = 1.15967481E-02;
  COFD[244] = -1.40524219E+01;
  COFD[245] = 3.56261348E+00;
  COFD[246] = -2.48287981E-01;
  COFD[247] = 1.07752947E-02;
  COFD[248] = -1.72994129E+01;
  COFD[249] = 4.71931868E+00;
  COFD[250] = -3.91258152E-01;
  COFD[251] = 1.66866639E-02;
  COFD[252] = -1.53241298E+01;
  COFD[253] = 4.02249494E+00;
  COFD[254] = -3.06112211E-01;
  COFD[255] = 1.32056535E-02;
  COFD[256] = -1.53791907E+01;
  COFD[257] = 4.03496578E+00;
  COFD[258] = -3.07489899E-01;
  COFD[259] = 1.32554084E-02;
  COFD[260] = -1.53769382E+01;
  COFD[261] = 4.03496578E+00;
  COFD[262] = -3.07489899E-01;
  COFD[263] = 1.32554084E-02;
  COFD[264] = -1.53769382E+01;
  COFD[265] = 4.03496578E+00;
  COFD[266] = -3.07489899E-01;
  COFD[267] = 1.32554084E-02;
  COFD[268] = -1.30215958E+01;
  COFD[269] = 3.24041848E+00;
  COFD[270] = -2.07804876E-01;
  COFD[271] = 9.07722400E-03;
  COFD[272] = -1.31792887E+01;
  COFD[273] = 3.38003453E+00;
  COFD[274] = -2.25783856E-01;
  COFD[275] = 9.85028660E-03;
  COFD[276] = -1.94386101E+01;
  COFD[277] = 5.43670321E+00;
  COFD[278] = -4.70985776E-01;
  COFD[279] = 1.96158515E-02;
  COFD[280] = -1.90281731E+01;
  COFD[281] = 5.30278138E+00;
  COFD[282] = -4.57079018E-01;
  COFD[283] = 1.91486428E-02;
  COFD[284] = -1.90281731E+01;
  COFD[285] = 5.30278138E+00;
  COFD[286] = -4.57079018E-01;
  COFD[287] = 1.91486428E-02;
  COFD[288] = -1.95318517E+01;
  COFD[289] = 5.47046983E+00;
  COFD[290] = -4.74577605E-01;
  COFD[291] = 1.97408822E-02;
  COFD[292] = -1.95312879E+01;
  COFD[293] = 5.47046983E+00;
  COFD[294] = -4.74577605E-01;
  COFD[295] = 1.97408822E-02;
  COFD[296] = -1.31995018E+01;
  COFD[297] = 3.38003453E+00;
  COFD[298] = -2.25783856E-01;
  COFD[299] = 9.85028660E-03;
  COFD[300] = -1.41953189E+01;
  COFD[301] = 3.68995707E+00;
  COFD[302] = -2.64784442E-01;
  COFD[303] = 1.14891794E-02;
  COFD[304] = -1.51177800E+01;
  COFD[305] = 3.99904647E+00;
  COFD[306] = -3.03517220E-01;
  COFD[307] = 1.31117363E-02;
  COFD[308] = -2.50118194E+01;
  COFD[309] = 5.51533435E+00;
  COFD[310] = -4.61886312E-01;
  COFD[311] = 1.84746261E-02;
  COFD[312] = -1.50096252E+01;
  COFD[313] = 3.25515933E+00;
  COFD[314] = -2.09710110E-01;
  COFD[315] = 9.15941830E-03;
  COFD[316] = -1.15797897E+01;
  COFD[317] = 2.43235504E+00;
  COFD[318] = -1.02890179E-01;
  COFD[319] = 4.52903603E-03;
  COFD[320] = -1.43713019E+01;
  COFD[321] = 3.70920439E+00;
  COFD[322] = -2.67274113E-01;
  COFD[323] = 1.15967481E-02;
  COFD[324] = -1.53110597E+01;
  COFD[325] = 3.37317428E+00;
  COFD[326] = -2.24900439E-01;
  COFD[327] = 9.81228151E-03;
  COFD[328] = -1.40864784E+01;
  COFD[329] = 3.07458927E+00;
  COFD[330] = -1.86899591E-01;
  COFD[331] = 8.19829781E-03;
  COFD[332] = -2.10639928E+01;
  COFD[333] = 5.50980708E+00;
  COFD[334] = -4.78335495E-01;
  COFD[335] = 1.98515433E-02;
  COFD[336] = -1.41066358E+01;
  COFD[337] = 3.07458927E+00;
  COFD[338] = -1.86899591E-01;
  COFD[339] = 8.19829781E-03;
  COFD[340] = -1.53261011E+01;
  COFD[341] = 3.37317428E+00;
  COFD[342] = -2.24900439E-01;
  COFD[343] = 9.81228151E-03;
  COFD[344] = -1.53187536E+01;
  COFD[345] = 3.37317428E+00;
  COFD[346] = -2.24900439E-01;
  COFD[347] = 9.81228151E-03;
  COFD[348] = -1.50371688E+01;
  COFD[349] = 3.26249588E+00;
  COFD[350] = -2.10658287E-01;
  COFD[351] = 9.20032462E-03;
  COFD[352] = -1.81197250E+01;
  COFD[353] = 4.33684042E+00;
  COFD[354] = -3.44981265E-01;
  COFD[355] = 1.48142449E-02;
  COFD[356] = -1.59878735E+01;
  COFD[357] = 3.66478157E+00;
  COFD[358] = -2.61506432E-01;
  COFD[359] = 1.13466231E-02;
  COFD[360] = -1.60662884E+01;
  COFD[361] = 3.69744368E+00;
  COFD[362] = -2.65754489E-01;
  COFD[363] = 1.15311681E-02;
  COFD[364] = -1.60424224E+01;
  COFD[365] = 3.69744368E+00;
  COFD[366] = -2.65754489E-01;
  COFD[367] = 1.15311681E-02;
  COFD[368] = -1.60424224E+01;
  COFD[369] = 3.69744368E+00;
  COFD[370] = -2.65754489E-01;
  COFD[371] = 1.15311681E-02;
  COFD[372] = -1.38694480E+01;
  COFD[373] = 2.96800700E+00;
  COFD[374] = -1.73235360E-01;
  COFD[375] = 7.61377726E-03;
  COFD[376] = -1.40154693E+01;
  COFD[377] = 3.07458927E+00;
  COFD[378] = -1.86899591E-01;
  COFD[379] = 8.19829781E-03;
  COFD[380] = -2.04422993E+01;
  COFD[381] = 5.18126516E+00;
  COFD[382] = -4.45006101E-01;
  COFD[383] = 1.87713513E-02;
  COFD[384] = -2.02220865E+01;
  COFD[385] = 5.10234461E+00;
  COFD[386] = -4.36039243E-01;
  COFD[387] = 1.84326773E-02;
  COFD[388] = -2.02220865E+01;
  COFD[389] = 5.10234461E+00;
  COFD[390] = -4.36039243E-01;
  COFD[391] = 1.84326773E-02;
  COFD[392] = -2.05196205E+01;
  COFD[393] = 5.21473296E+00;
  COFD[394] = -4.48646311E-01;
  COFD[395] = 1.89013813E-02;
  COFD[396] = -2.05107395E+01;
  COFD[397] = 5.21473296E+00;
  COFD[398] = -4.48646311E-01;
  COFD[399] = 1.89013813E-02;
  COFD[400] = -1.42641696E+01;
  COFD[401] = 3.07458927E+00;
  COFD[402] = -1.86899591E-01;
  COFD[403] = 8.19829781E-03;
  COFD[404] = -1.50366933E+01;
  COFD[405] = 3.35941694E+00;
  COFD[406] = -2.23126552E-01;
  COFD[407] = 9.73586922E-03;
  COFD[408] = -1.59068122E+01;
  COFD[409] = 3.60186887E+00;
  COFD[410] = -2.53302622E-01;
  COFD[411] = 1.09893496E-02;
  COFD[412] = -2.47604202E+01;
  COFD[413] = 5.51868821E+00;
  COFD[414] = -4.73225401E-01;
  COFD[415] = 1.93849523E-02;
  COFD[416] = -1.38756380E+01;
  COFD[417] = 2.98558426E+00;
  COFD[418] = -1.75507216E-01;
  COFD[419] = 7.71173691E-03;
  COFD[420] = -1.06250314E+01;
  COFD[421] = 2.15849701E+00;
  COFD[422] = -6.53886401E-02;
  COFD[423] = 2.81453370E-03;
  COFD[424] = -1.31860264E+01;
  COFD[425] = 3.38003453E+00;
  COFD[426] = -2.25783856E-01;
  COFD[427] = 9.85028660E-03;
  COFD[428] = -1.40864784E+01;
  COFD[429] = 3.07458927E+00;
  COFD[430] = -1.86899591E-01;
  COFD[431] = 8.19829781E-03;
  COFD[432] = -1.29877255E+01;
  COFD[433] = 2.80841511E+00;
  COFD[434] = -1.52629888E-01;
  COFD[435] = 6.72604927E-03;
  COFD[436] = -1.91096700E+01;
  COFD[437] = 5.02608697E+00;
  COFD[438] = -4.26959993E-01;
  COFD[439] = 1.80709910E-02;
  COFD[440] = -1.30027669E+01;
  COFD[441] = 2.80841511E+00;
  COFD[442] = -1.52629888E-01;
  COFD[443] = 6.72604927E-03;
  COFD[444] = -1.40964555E+01;
  COFD[445] = 3.07458927E+00;
  COFD[446] = -1.86899591E-01;
  COFD[447] = 8.19829781E-03;
  COFD[448] = -1.40915944E+01;
  COFD[449] = 3.07458927E+00;
  COFD[450] = -1.86899591E-01;
  COFD[451] = 8.19829781E-03;
  COFD[452] = -1.39007309E+01;
  COFD[453] = 2.99164244E+00;
  COFD[454] = -1.76293106E-01;
  COFD[455] = 7.74575100E-03;
  COFD[456] = -1.67115471E+01;
  COFD[457] = 3.98859394E+00;
  COFD[458] = -3.02316219E-01;
  COFD[459] = 1.30661099E-02;
  COFD[460] = -1.47765660E+01;
  COFD[461] = 3.34934830E+00;
  COFD[462] = -2.21831026E-01;
  COFD[463] = 9.68019791E-03;
  COFD[464] = -1.48358026E+01;
  COFD[465] = 3.37160686E+00;
  COFD[466] = -2.24698307E-01;
  COFD[467] = 9.80357301E-03;
  COFD[468] = -1.48176140E+01;
  COFD[469] = 3.37160686E+00;
  COFD[470] = -2.24698307E-01;
  COFD[471] = 9.80357301E-03;
  COFD[472] = -1.48176140E+01;
  COFD[473] = 3.37160686E+00;
  COFD[474] = -2.24698307E-01;
  COFD[475] = 9.80357301E-03;
  COFD[476] = -1.28155338E+01;
  COFD[477] = 2.69967657E+00;
  COFD[478] = -1.38362256E-01;
  COFD[479] = 6.10143187E-03;
  COFD[480] = -1.29335458E+01;
  COFD[481] = 2.80841511E+00;
  COFD[482] = -1.52629888E-01;
  COFD[483] = 6.72604927E-03;
  COFD[484] = -1.90163591E+01;
  COFD[485] = 4.83982263E+00;
  COFD[486] = -4.05783981E-01;
  COFD[487] = 1.72721335E-02;
  COFD[488] = -1.84871270E+01;
  COFD[489] = 4.64907964E+00;
  COFD[490] = -3.82960584E-01;
  COFD[491] = 1.63598480E-02;
  COFD[492] = -1.84871270E+01;
  COFD[493] = 4.64907964E+00;
  COFD[494] = -3.82960584E-01;
  COFD[495] = 1.63598480E-02;
  COFD[496] = -1.91105534E+01;
  COFD[497] = 4.87977047E+00;
  COFD[498] = -4.10448693E-01;
  COFD[499] = 1.74535827E-02;
  COFD[500] = -1.91045520E+01;
  COFD[501] = 4.87977047E+00;
  COFD[502] = -4.10448693E-01;
  COFD[503] = 1.74535827E-02;
  COFD[504] = -1.31147338E+01;
  COFD[505] = 2.80841511E+00;
  COFD[506] = -1.52629888E-01;
  COFD[507] = 6.72604927E-03;
  COFD[508] = -1.38538393E+01;
  COFD[509] = 3.06403296E+00;
  COFD[510] = -1.85549879E-01;
  COFD[511] = 8.14070073E-03;
  COFD[512] = -1.46739282E+01;
  COFD[513] = 3.30683499E+00;
  COFD[514] = -2.16378602E-01;
  COFD[515] = 9.44670561E-03;
  COFD[516] = -2.42503398E+01;
  COFD[517] = 5.40793558E+00;
  COFD[518] = -4.67888803E-01;
  COFD[519] = 1.95062704E-02;
  COFD[520] = -2.08123335E+01;
  COFD[521] = 5.42470187E+00;
  COFD[522] = -4.69700452E-01;
  COFD[523] = 1.95706916E-02;
  COFD[524] = -1.68759082E+01;
  COFD[525] = 4.49460364E+00;
  COFD[526] = -3.64766207E-01;
  COFD[527] = 1.56457185E-02;
  COFD[528] = -1.93611201E+01;
  COFD[529] = 5.51579726E+00;
  COFD[530] = -4.76061961E-01;
  COFD[531] = 1.96329391E-02;
  COFD[532] = -2.10639928E+01;
  COFD[533] = 5.50980708E+00;
  COFD[534] = -4.78335495E-01;
  COFD[535] = 1.98515433E-02;
  COFD[536] = -1.91096700E+01;
  COFD[537] = 5.02608697E+00;
  COFD[538] = -4.26959993E-01;
  COFD[539] = 1.80709910E-02;
  COFD[540] = -1.31492477E+01;
  COFD[541] = 1.48003884E+00;
  COFD[542] = 1.60500200E-01;
  COFD[543] = -1.19765987E-02;
  COFD[544] = -1.91256173E+01;
  COFD[545] = 5.02608697E+00;
  COFD[546] = -4.26959993E-01;
  COFD[547] = 1.80709910E-02;
  COFD[548] = -2.04229988E+01;
  COFD[549] = 5.31457079E+00;
  COFD[550] = -4.58216496E-01;
  COFD[551] = 1.91825910E-02;
  COFD[552] = -2.04177394E+01;
  COFD[553] = 5.31457079E+00;
  COFD[554] = -4.58216496E-01;
  COFD[555] = 1.91825910E-02;
  COFD[556] = -2.08943727E+01;
  COFD[557] = 5.44718685E+00;
  COFD[558] = -4.72082987E-01;
  COFD[559] = 1.96531332E-02;
  COFD[560] = -2.12021417E+01;
  COFD[561] = 5.20775012E+00;
  COFD[562] = -4.07348260E-01;
  COFD[563] = 1.55473249E-02;
  COFD[564] = -2.12953626E+01;
  COFD[565] = 5.52385259E+00;
  COFD[566] = -4.69683808E-01;
  COFD[567] = 1.90677248E-02;
  COFD[568] = -2.10691956E+01;
  COFD[569] = 5.51604369E+00;
  COFD[570] = -4.76337283E-01;
  COFD[571] = 1.96544379E-02;
  COFD[572] = -2.10499842E+01;
  COFD[573] = 5.51604369E+00;
  COFD[574] = -4.76337283E-01;
  COFD[575] = 1.96544379E-02;
  COFD[576] = -2.10499842E+01;
  COFD[577] = 5.51604369E+00;
  COFD[578] = -4.76337283E-01;
  COFD[579] = 1.96544379E-02;
  COFD[580] = -1.88701223E+01;
  COFD[581] = 4.90944038E+00;
  COFD[582] = -4.13866957E-01;
  COFD[583] = 1.75844447E-02;
  COFD[584] = -1.90524563E+01;
  COFD[585] = 5.02608697E+00;
  COFD[586] = -4.26959993E-01;
  COFD[587] = 1.80709910E-02;
  COFD[588] = -1.89495722E+01;
  COFD[589] = 4.10779436E+00;
  COFD[590] = -2.35955221E-01;
  COFD[591] = 7.03410627E-03;
  COFD[592] = -1.75035239E+01;
  COFD[593] = 3.41230274E+00;
  COFD[594] = -1.32406923E-01;
  COFD[595] = 2.10435795E-03;
  COFD[596] = -1.75035239E+01;
  COFD[597] = 3.41230274E+00;
  COFD[598] = -1.32406923E-01;
  COFD[599] = 2.10435795E-03;
  COFD[600] = -1.87236049E+01;
  COFD[601] = 4.00967621E+00;
  COFD[602] = -2.21153539E-01;
  COFD[603] = 6.31528745E-03;
  COFD[604] = -1.87171338E+01;
  COFD[605] = 4.00967621E+00;
  COFD[606] = -2.21153539E-01;
  COFD[607] = 6.31528745E-03;
  COFD[608] = -2.04933546E+01;
  COFD[609] = 5.41191809E+00;
  COFD[610] = -4.68316503E-01;
  COFD[611] = 1.95213604E-02;
  COFD[612] = -2.07061954E+01;
  COFD[613] = 5.46880564E+00;
  COFD[614] = -4.74397826E-01;
  COFD[615] = 1.97344800E-02;
  COFD[616] = -2.12711805E+01;
  COFD[617] = 5.53137765E+00;
  COFD[618] = -4.73118467E-01;
  COFD[619] = 1.93125852E-02;
  COFD[620] = -1.78173107E+01;
  COFD[621] = 1.97328753E+00;
  COFD[622] = 7.64115525E-02;
  COFD[623] = -7.78429580E-03;
  COFD[624] = -1.38948654E+01;
  COFD[625] = 2.98558426E+00;
  COFD[626] = -1.75507216E-01;
  COFD[627] = 7.71173691E-03;
  COFD[628] = -1.06283588E+01;
  COFD[629] = 2.15849701E+00;
  COFD[630] = -6.53886401E-02;
  COFD[631] = 2.81453370E-03;
  COFD[632] = -1.31877859E+01;
  COFD[633] = 3.38003453E+00;
  COFD[634] = -2.25783856E-01;
  COFD[635] = 9.85028660E-03;
  COFD[636] = -1.41066358E+01;
  COFD[637] = 3.07458927E+00;
  COFD[638] = -1.86899591E-01;
  COFD[639] = 8.19829781E-03;
  COFD[640] = -1.30027669E+01;
  COFD[641] = 2.80841511E+00;
  COFD[642] = -1.52629888E-01;
  COFD[643] = 6.72604927E-03;
  COFD[644] = -1.91256173E+01;
  COFD[645] = 5.02608697E+00;
  COFD[646] = -4.26959993E-01;
  COFD[647] = 1.80709910E-02;
  COFD[648] = -1.30182749E+01;
  COFD[649] = 2.80841511E+00;
  COFD[650] = -1.52629888E-01;
  COFD[651] = 6.72604927E-03;
  COFD[652] = -1.41170278E+01;
  COFD[653] = 3.07458927E+00;
  COFD[654] = -1.86899591E-01;
  COFD[655] = 8.19829781E-03;
  COFD[656] = -1.41119635E+01;
  COFD[657] = 3.07458927E+00;
  COFD[658] = -1.86899591E-01;
  COFD[659] = 8.19829781E-03;
  COFD[660] = -1.39199573E+01;
  COFD[661] = 2.99164244E+00;
  COFD[662] = -1.76293106E-01;
  COFD[663] = 7.74575100E-03;
  COFD[664] = -1.67337674E+01;
  COFD[665] = 3.98859394E+00;
  COFD[666] = -3.02316219E-01;
  COFD[667] = 1.30661099E-02;
  COFD[668] = -1.47916284E+01;
  COFD[669] = 3.34934830E+00;
  COFD[670] = -2.21831026E-01;
  COFD[671] = 9.68019791E-03;
  COFD[672] = -1.48503699E+01;
  COFD[673] = 3.37160686E+00;
  COFD[674] = -2.24698307E-01;
  COFD[675] = 9.80357301E-03;
  COFD[676] = -1.48316536E+01;
  COFD[677] = 3.37160686E+00;
  COFD[678] = -2.24698307E-01;
  COFD[679] = 9.80357301E-03;
  COFD[680] = -1.48316536E+01;
  COFD[681] = 3.37160686E+00;
  COFD[682] = -2.24698307E-01;
  COFD[683] = 9.80357301E-03;
  COFD[684] = -1.28284059E+01;
  COFD[685] = 2.69967657E+00;
  COFD[686] = -1.38362256E-01;
  COFD[687] = 6.10143187E-03;
  COFD[688] = -1.29470215E+01;
  COFD[689] = 2.80841511E+00;
  COFD[690] = -1.52629888E-01;
  COFD[691] = 6.72604927E-03;
  COFD[692] = -1.90365260E+01;
  COFD[693] = 4.83982263E+00;
  COFD[694] = -4.05783981E-01;
  COFD[695] = 1.72721335E-02;
  COFD[696] = -1.85070738E+01;
  COFD[697] = 4.64907964E+00;
  COFD[698] = -3.82960584E-01;
  COFD[699] = 1.63598480E-02;
  COFD[700] = -1.85070738E+01;
  COFD[701] = 4.64907964E+00;
  COFD[702] = -3.82960584E-01;
  COFD[703] = 1.63598480E-02;
  COFD[704] = -1.91302705E+01;
  COFD[705] = 4.87977047E+00;
  COFD[706] = -4.10448693E-01;
  COFD[707] = 1.74535827E-02;
  COFD[708] = -1.91240292E+01;
  COFD[709] = 4.87977047E+00;
  COFD[710] = -4.10448693E-01;
  COFD[711] = 1.74535827E-02;
  COFD[712] = -1.31342108E+01;
  COFD[713] = 2.80841511E+00;
  COFD[714] = -1.52629888E-01;
  COFD[715] = 6.72604927E-03;
  COFD[716] = -1.38702007E+01;
  COFD[717] = 3.06403296E+00;
  COFD[718] = -1.85549879E-01;
  COFD[719] = 8.14070073E-03;
  COFD[720] = -1.46940858E+01;
  COFD[721] = 3.30683499E+00;
  COFD[722] = -2.16378602E-01;
  COFD[723] = 9.44670561E-03;
  COFD[724] = -2.42503408E+01;
  COFD[725] = 5.40793558E+00;
  COFD[726] = -4.67888803E-01;
  COFD[727] = 1.95062704E-02;
  COFD[728] = -1.50236540E+01;
  COFD[729] = 3.25515933E+00;
  COFD[730] = -2.09710110E-01;
  COFD[731] = 9.15941830E-03;
  COFD[732] = -1.15815492E+01;
  COFD[733] = 2.43235504E+00;
  COFD[734] = -1.02890179E-01;
  COFD[735] = 4.52903603E-03;
  COFD[736] = -1.43722077E+01;
  COFD[737] = 3.70920439E+00;
  COFD[738] = -2.67274113E-01;
  COFD[739] = 1.15967481E-02;
  COFD[740] = -1.53261011E+01;
  COFD[741] = 3.37317428E+00;
  COFD[742] = -2.24900439E-01;
  COFD[743] = 9.81228151E-03;
  COFD[744] = -1.40964555E+01;
  COFD[745] = 3.07458927E+00;
  COFD[746] = -1.86899591E-01;
  COFD[747] = 8.19829781E-03;
  COFD[748] = -2.04229988E+01;
  COFD[749] = 5.31457079E+00;
  COFD[750] = -4.58216496E-01;
  COFD[751] = 1.91825910E-02;
  COFD[752] = -1.41170278E+01;
  COFD[753] = 3.07458927E+00;
  COFD[754] = -1.86899591E-01;
  COFD[755] = 8.19829781E-03;
  COFD[756] = -1.53416091E+01;
  COFD[757] = 3.37317428E+00;
  COFD[758] = -2.24900439E-01;
  COFD[759] = 9.81228151E-03;
  COFD[760] = -1.53340319E+01;
  COFD[761] = 3.37317428E+00;
  COFD[762] = -2.24900439E-01;
  COFD[763] = 9.81228151E-03;
  COFD[764] = -1.50511964E+01;
  COFD[765] = 3.26249588E+00;
  COFD[766] = -2.10658287E-01;
  COFD[767] = 9.20032462E-03;
  COFD[768] = -1.81371853E+01;
  COFD[769] = 4.33684042E+00;
  COFD[770] = -3.44981265E-01;
  COFD[771] = 1.48142449E-02;
  COFD[772] = -1.59978691E+01;
  COFD[773] = 3.66478157E+00;
  COFD[774] = -2.61506432E-01;
  COFD[775] = 1.13466231E-02;
  COFD[776] = -1.60758526E+01;
  COFD[777] = 3.69744368E+00;
  COFD[778] = -2.65754489E-01;
  COFD[779] = 1.15311681E-02;
  COFD[780] = -1.60515367E+01;
  COFD[781] = 3.69744368E+00;
  COFD[782] = -2.65754489E-01;
  COFD[783] = 1.15311681E-02;
  COFD[784] = -1.60515367E+01;
  COFD[785] = 3.69744368E+00;
  COFD[786] = -2.65754489E-01;
  COFD[787] = 1.15311681E-02;
  COFD[788] = -1.38776021E+01;
  COFD[789] = 2.96800700E+00;
  COFD[790] = -1.73235360E-01;
  COFD[791] = 7.61377726E-03;
  COFD[792] = -1.40241141E+01;
  COFD[793] = 3.07458927E+00;
  COFD[794] = -1.86899591E-01;
  COFD[795] = 8.19829781E-03;
  COFD[796] = -2.04573512E+01;
  COFD[797] = 5.18126516E+00;
  COFD[798] = -4.45006101E-01;
  COFD[799] = 1.87713513E-02;
  COFD[800] = -1.99478621E+01;
  COFD[801] = 5.00532670E+00;
  COFD[802] = -4.24686880E-01;
  COFD[803] = 1.79891229E-02;
  COFD[804] = -1.99478621E+01;
  COFD[805] = 5.00532670E+00;
  COFD[806] = -4.24686880E-01;
  COFD[807] = 1.79891229E-02;
  COFD[808] = -2.05341767E+01;
  COFD[809] = 5.21473296E+00;
  COFD[810] = -4.48646311E-01;
  COFD[811] = 1.89013813E-02;
  COFD[812] = -2.05250357E+01;
  COFD[813] = 5.21473296E+00;
  COFD[814] = -4.48646311E-01;
  COFD[815] = 1.89013813E-02;
  COFD[816] = -1.42784657E+01;
  COFD[817] = 3.07458927E+00;
  COFD[818] = -1.86899591E-01;
  COFD[819] = 8.19829781E-03;
  COFD[820] = -1.50478662E+01;
  COFD[821] = 3.35941694E+00;
  COFD[822] = -2.23126552E-01;
  COFD[823] = 9.73586922E-03;
  COFD[824] = -1.59218538E+01;
  COFD[825] = 3.60186887E+00;
  COFD[826] = -2.53302622E-01;
  COFD[827] = 1.09893496E-02;
  COFD[828] = -2.47604207E+01;
  COFD[829] = 5.51868821E+00;
  COFD[830] = -4.73225401E-01;
  COFD[831] = 1.93849523E-02;
  COFD[832] = -1.50168047E+01;
  COFD[833] = 3.25515933E+00;
  COFD[834] = -2.09710110E-01;
  COFD[835] = 9.15941830E-03;
  COFD[836] = -1.15806955E+01;
  COFD[837] = 2.43235504E+00;
  COFD[838] = -1.02890179E-01;
  COFD[839] = 4.52903603E-03;
  COFD[840] = -1.43717684E+01;
  COFD[841] = 3.70920439E+00;
  COFD[842] = -2.67274113E-01;
  COFD[843] = 1.15967481E-02;
  COFD[844] = -1.53187536E+01;
  COFD[845] = 3.37317428E+00;
  COFD[846] = -2.24900439E-01;
  COFD[847] = 9.81228151E-03;
  COFD[848] = -1.40915944E+01;
  COFD[849] = 3.07458927E+00;
  COFD[850] = -1.86899591E-01;
  COFD[851] = 8.19829781E-03;
  COFD[852] = -2.04177394E+01;
  COFD[853] = 5.31457079E+00;
  COFD[854] = -4.58216496E-01;
  COFD[855] = 1.91825910E-02;
  COFD[856] = -1.41119635E+01;
  COFD[857] = 3.07458927E+00;
  COFD[858] = -1.86899591E-01;
  COFD[859] = 8.19829781E-03;
  COFD[860] = -1.53340319E+01;
  COFD[861] = 3.37317428E+00;
  COFD[862] = -2.24900439E-01;
  COFD[863] = 9.81228151E-03;
  COFD[864] = -1.53265677E+01;
  COFD[865] = 3.37317428E+00;
  COFD[866] = -2.24900439E-01;
  COFD[867] = 9.81228151E-03;
  COFD[868] = -1.50443476E+01;
  COFD[869] = 3.26249588E+00;
  COFD[870] = -2.10658287E-01;
  COFD[871] = 9.20032462E-03;
  COFD[872] = -1.81286456E+01;
  COFD[873] = 4.33684042E+00;
  COFD[874] = -3.44981265E-01;
  COFD[875] = 1.48142449E-02;
  COFD[876] = -1.59929990E+01;
  COFD[877] = 3.66478157E+00;
  COFD[878] = -2.61506432E-01;
  COFD[879] = 1.13466231E-02;
  COFD[880] = -1.60711937E+01;
  COFD[881] = 3.69744368E+00;
  COFD[882] = -2.65754489E-01;
  COFD[883] = 1.15311681E-02;
  COFD[884] = -1.60470980E+01;
  COFD[885] = 3.69744368E+00;
  COFD[886] = -2.65754489E-01;
  COFD[887] = 1.15311681E-02;
  COFD[888] = -1.60470980E+01;
  COFD[889] = 3.69744368E+00;
  COFD[890] = -2.65754489E-01;
  COFD[891] = 1.15311681E-02;
  COFD[892] = -1.38736329E+01;
  COFD[893] = 2.96800700E+00;
  COFD[894] = -1.73235360E-01;
  COFD[895] = 7.61377726E-03;
  COFD[896] = -1.40199050E+01;
  COFD[897] = 3.07458927E+00;
  COFD[898] = -1.86899591E-01;
  COFD[899] = 8.19829781E-03;
  COFD[900] = -2.04499985E+01;
  COFD[901] = 5.18126516E+00;
  COFD[902] = -4.45006101E-01;
  COFD[903] = 1.87713513E-02;
  COFD[904] = -1.99406295E+01;
  COFD[905] = 5.00532670E+00;
  COFD[906] = -4.24686880E-01;
  COFD[907] = 1.79891229E-02;
  COFD[908] = -1.99406295E+01;
  COFD[909] = 5.00532670E+00;
  COFD[910] = -4.24686880E-01;
  COFD[911] = 1.79891229E-02;
  COFD[912] = -2.05270679E+01;
  COFD[913] = 5.21473296E+00;
  COFD[914] = -4.48646311E-01;
  COFD[915] = 1.89013813E-02;
  COFD[916] = -2.05180548E+01;
  COFD[917] = 5.21473296E+00;
  COFD[918] = -4.48646311E-01;
  COFD[919] = 1.89013813E-02;
  COFD[920] = -1.42714849E+01;
  COFD[921] = 3.07458927E+00;
  COFD[922] = -1.86899591E-01;
  COFD[923] = 8.19829781E-03;
  COFD[924] = -1.50424192E+01;
  COFD[925] = 3.35941694E+00;
  COFD[926] = -2.23126552E-01;
  COFD[927] = 9.73586922E-03;
  COFD[928] = -1.59145062E+01;
  COFD[929] = 3.60186887E+00;
  COFD[930] = -2.53302622E-01;
  COFD[931] = 1.09893496E-02;
  COFD[932] = -2.47604204E+01;
  COFD[933] = 5.51868821E+00;
  COFD[934] = -4.73225401E-01;
  COFD[935] = 1.93849523E-02;
  COFD[936] = -1.47850505E+01;
  COFD[937] = 3.16433919E+00;
  COFD[938] = -1.98191564E-01;
  COFD[939] = 8.67209742E-03;
  COFD[940] = -1.13541221E+01;
  COFD[941] = 2.31999438E+00;
  COFD[942] = -8.75064804E-02;
  COFD[943] = 3.82656365E-03;
  COFD[944] = -1.40524219E+01;
  COFD[945] = 3.56261348E+00;
  COFD[946] = -2.48287981E-01;
  COFD[947] = 1.07752947E-02;
  COFD[948] = -1.50371688E+01;
  COFD[949] = 3.26249588E+00;
  COFD[950] = -2.10658287E-01;
  COFD[951] = 9.20032462E-03;
  COFD[952] = -1.39007309E+01;
  COFD[953] = 2.99164244E+00;
  COFD[954] = -1.76293106E-01;
  COFD[955] = 7.74575100E-03;
  COFD[956] = -2.08943727E+01;
  COFD[957] = 5.44718685E+00;
  COFD[958] = -4.72082987E-01;
  COFD[959] = 1.96531332E-02;
  COFD[960] = -1.39199573E+01;
  COFD[961] = 2.99164244E+00;
  COFD[962] = -1.76293106E-01;
  COFD[963] = 7.74575100E-03;
  COFD[964] = -1.50511964E+01;
  COFD[965] = 3.26249588E+00;
  COFD[966] = -2.10658287E-01;
  COFD[967] = 9.20032462E-03;
  COFD[968] = -1.50443476E+01;
  COFD[969] = 3.26249588E+00;
  COFD[970] = -2.10658287E-01;
  COFD[971] = 9.20032462E-03;
  COFD[972] = -1.48061406E+01;
  COFD[973] = 3.16912473E+00;
  COFD[974] = -1.98792456E-01;
  COFD[975] = 8.69726395E-03;
  COFD[976] = -1.77672912E+01;
  COFD[977] = 4.20234040E+00;
  COFD[978] = -3.28057658E-01;
  COFD[979] = 1.41006192E-02;
  COFD[980] = -1.56953785E+01;
  COFD[981] = 3.54443207E+00;
  COFD[982] = -2.46133003E-01;
  COFD[983] = 1.06905018E-02;
  COFD[984] = -1.57272832E+01;
  COFD[985] = 3.55756440E+00;
  COFD[986] = -2.47689003E-01;
  COFD[987] = 1.07517019E-02;
  COFD[988] = -1.57044328E+01;
  COFD[989] = 3.55756440E+00;
  COFD[990] = -2.47689003E-01;
  COFD[991] = 1.07517019E-02;
  COFD[992] = -1.57044328E+01;
  COFD[993] = 3.55756440E+00;
  COFD[994] = -2.47689003E-01;
  COFD[995] = 1.07517019E-02;
  COFD[996] = -1.36800096E+01;
  COFD[997] = 2.88311045E+00;
  COFD[998] = -1.62274683E-01;
  COFD[999] = 7.14167394E-03;
  COFD[1000] = -1.38327302E+01;
  COFD[1001] = 2.99164244E+00;
  COFD[1002] = -1.76293106E-01;
  COFD[1003] = 7.74575100E-03;
  COFD[1004] = -2.01404432E+01;
  COFD[1005] = 5.07318708E+00;
  COFD[1006] = -4.32608485E-01;
  COFD[1007] = 1.82977801E-02;
  COFD[1008] = -1.99343637E+01;
  COFD[1009] = 4.99746599E+00;
  COFD[1010] = -4.23838967E-01;
  COFD[1011] = 1.79592096E-02;
  COFD[1012] = -1.99343637E+01;
  COFD[1013] = 4.99746599E+00;
  COFD[1014] = -4.23838967E-01;
  COFD[1015] = 1.79592096E-02;
  COFD[1016] = -2.02444913E+01;
  COFD[1017] = 5.11785645E+00;
  COFD[1018] = -4.37867828E-01;
  COFD[1019] = 1.85047543E-02;
  COFD[1020] = -2.02361782E+01;
  COFD[1021] = 5.11785645E+00;
  COFD[1022] = -4.37867828E-01;
  COFD[1023] = 1.85047543E-02;
  COFD[1024] = -1.40688472E+01;
  COFD[1025] = 2.99164244E+00;
  COFD[1026] = -1.76293106E-01;
  COFD[1027] = 7.74575100E-03;
  COFD[1028] = -1.47729612E+01;
  COFD[1029] = 3.24835477E+00;
  COFD[1030] = -2.08830747E-01;
  COFD[1031] = 9.12148470E-03;
  COFD[1032] = -1.56952963E+01;
  COFD[1033] = 3.51861272E+00;
  COFD[1034] = -2.43068621E-01;
  COFD[1035] = 1.05698368E-02;
  COFD[1036] = -2.46791388E+01;
  COFD[1037] = 5.51562323E+00;
  COFD[1038] = -4.75876136E-01;
  COFD[1039] = 1.96184735E-02;
  COFD[1040] = -1.77350630E+01;
  COFD[1041] = 4.19328271E+00;
  COFD[1042] = -3.26911461E-01;
  COFD[1043] = 1.40520357E-02;
  COFD[1044] = -1.35545391E+01;
  COFD[1045] = 3.13878730E+00;
  COFD[1046] = -1.94980335E-01;
  COFD[1047] = 8.53744486E-03;
  COFD[1048] = -1.72994129E+01;
  COFD[1049] = 4.71931868E+00;
  COFD[1050] = -3.91258152E-01;
  COFD[1051] = 1.66866639E-02;
  COFD[1052] = -1.81197250E+01;
  COFD[1053] = 4.33684042E+00;
  COFD[1054] = -3.44981265E-01;
  COFD[1055] = 1.48142449E-02;
  COFD[1056] = -1.67115471E+01;
  COFD[1057] = 3.98859394E+00;
  COFD[1058] = -3.02316219E-01;
  COFD[1059] = 1.30661099E-02;
  COFD[1060] = -2.12021417E+01;
  COFD[1061] = 5.20775012E+00;
  COFD[1062] = -4.07348260E-01;
  COFD[1063] = 1.55473249E-02;
  COFD[1064] = -1.67337674E+01;
  COFD[1065] = 3.98859394E+00;
  COFD[1066] = -3.02316219E-01;
  COFD[1067] = 1.30661099E-02;
  COFD[1068] = -1.81371853E+01;
  COFD[1069] = 4.33684042E+00;
  COFD[1070] = -3.44981265E-01;
  COFD[1071] = 1.48142449E-02;
  COFD[1072] = -1.81286456E+01;
  COFD[1073] = 4.33684042E+00;
  COFD[1074] = -3.44981265E-01;
  COFD[1075] = 1.48142449E-02;
  COFD[1076] = -1.77672912E+01;
  COFD[1077] = 4.20234040E+00;
  COFD[1078] = -3.28057658E-01;
  COFD[1079] = 1.41006192E-02;
  COFD[1080] = -2.10907633E+01;
  COFD[1081] = 5.29211327E+00;
  COFD[1082] = -4.56068366E-01;
  COFD[1083] = 1.91195062E-02;
  COFD[1084] = -1.89594762E+01;
  COFD[1085] = 4.68735220E+00;
  COFD[1086] = -3.87449514E-01;
  COFD[1087] = 1.65352261E-02;
  COFD[1088] = -1.90110170E+01;
  COFD[1089] = 4.70842146E+00;
  COFD[1090] = -3.89924637E-01;
  COFD[1091] = 1.66320934E-02;
  COFD[1092] = -1.89849283E+01;
  COFD[1093] = 4.70842146E+00;
  COFD[1094] = -3.89924637E-01;
  COFD[1095] = 1.66320934E-02;
  COFD[1096] = -1.89849283E+01;
  COFD[1097] = 4.70842146E+00;
  COFD[1098] = -3.89924637E-01;
  COFD[1099] = 1.66320934E-02;
  COFD[1100] = -1.64069432E+01;
  COFD[1101] = 3.85900456E+00;
  COFD[1102] = -2.86199757E-01;
  COFD[1103] = 1.23960157E-02;
  COFD[1104] = -1.66339567E+01;
  COFD[1105] = 3.98859394E+00;
  COFD[1106] = -3.02316219E-01;
  COFD[1107] = 1.30661099E-02;
  COFD[1108] = -2.21011224E+01;
  COFD[1109] = 5.52816937E+00;
  COFD[1110] = -4.65357386E-01;
  COFD[1111] = 1.86939300E-02;
  COFD[1112] = -2.20276430E+01;
  COFD[1113] = 5.52121111E+00;
  COFD[1114] = -4.67224743E-01;
  COFD[1115] = 1.88799020E-02;
  COFD[1116] = -2.20276430E+01;
  COFD[1117] = 5.52121111E+00;
  COFD[1118] = -4.67224743E-01;
  COFD[1119] = 1.88799020E-02;
  COFD[1120] = -2.20860863E+01;
  COFD[1121] = 5.52171573E+00;
  COFD[1122] = -4.63284984E-01;
  COFD[1123] = 1.85570924E-02;
  COFD[1124] = -2.20758670E+01;
  COFD[1125] = 5.52171573E+00;
  COFD[1126] = -4.63284984E-01;
  COFD[1127] = 1.85570924E-02;
  COFD[1128] = -1.69110608E+01;
  COFD[1129] = 3.98859394E+00;
  COFD[1130] = -3.02316219E-01;
  COFD[1131] = 1.30661099E-02;
  COFD[1132] = -1.78202286E+01;
  COFD[1133] = 4.31993263E+00;
  COFD[1134] = -3.42851630E-01;
  COFD[1135] = 1.47243122E-02;
  COFD[1136] = -1.89571393E+01;
  COFD[1137] = 4.64763677E+00;
  COFD[1138] = -3.82799418E-01;
  COFD[1139] = 1.63539171E-02;
  COFD[1140] = -2.38843968E+01;
  COFD[1141] = 4.82941387E+00;
  COFD[1142] = -3.46543093E-01;
  COFD[1143] = 1.24648862E-02;
  COFD[1144] = -1.56756124E+01;
  COFD[1145] = 3.54035770E+00;
  COFD[1146] = -2.45653442E-01;
  COFD[1147] = 1.06717969E-02;
  COFD[1148] = -1.22114409E+01;
  COFD[1149] = 2.65215786E+00;
  COFD[1150] = -1.32206766E-01;
  COFD[1151] = 5.83532580E-03;
  COFD[1152] = -1.53241298E+01;
  COFD[1153] = 4.02249494E+00;
  COFD[1154] = -3.06112211E-01;
  COFD[1155] = 1.32056535E-02;
  COFD[1156] = -1.59878735E+01;
  COFD[1157] = 3.66478157E+00;
  COFD[1158] = -2.61506432E-01;
  COFD[1159] = 1.13466231E-02;
  COFD[1160] = -1.47765660E+01;
  COFD[1161] = 3.34934830E+00;
  COFD[1162] = -2.21831026E-01;
  COFD[1163] = 9.68019791E-03;
  COFD[1164] = -2.12953626E+01;
  COFD[1165] = 5.52385259E+00;
  COFD[1166] = -4.69683808E-01;
  COFD[1167] = 1.90677248E-02;
  COFD[1168] = -1.47916284E+01;
  COFD[1169] = 3.34934830E+00;
  COFD[1170] = -2.21831026E-01;
  COFD[1171] = 9.68019791E-03;
  COFD[1172] = -1.59978691E+01;
  COFD[1173] = 3.66478157E+00;
  COFD[1174] = -2.61506432E-01;
  COFD[1175] = 1.13466231E-02;
  COFD[1176] = -1.59929990E+01;
  COFD[1177] = 3.66478157E+00;
  COFD[1178] = -2.61506432E-01;
  COFD[1179] = 1.13466231E-02;
  COFD[1180] = -1.56953785E+01;
  COFD[1181] = 3.54443207E+00;
  COFD[1182] = -2.46133003E-01;
  COFD[1183] = 1.06905018E-02;
  COFD[1184] = -1.89594762E+01;
  COFD[1185] = 4.68735220E+00;
  COFD[1186] = -3.87449514E-01;
  COFD[1187] = 1.65352261E-02;
  COFD[1188] = -1.68532873E+01;
  COFD[1189] = 4.00572564E+00;
  COFD[1190] = -3.04255586E-01;
  COFD[1191] = 1.31384083E-02;
  COFD[1192] = -1.68880017E+01;
  COFD[1193] = 4.01759713E+00;
  COFD[1194] = -3.05560618E-01;
  COFD[1195] = 1.31852111E-02;
  COFD[1196] = -1.68697894E+01;
  COFD[1197] = 4.01759713E+00;
  COFD[1198] = -3.05560618E-01;
  COFD[1199] = 1.31852111E-02;
  COFD[1200] = -1.68697894E+01;
  COFD[1201] = 4.01759713E+00;
  COFD[1202] = -3.05560618E-01;
  COFD[1203] = 1.31852111E-02;
  COFD[1204] = -1.44868574E+01;
  COFD[1205] = 3.20806285E+00;
  COFD[1206] = -2.03619458E-01;
  COFD[1207] = 8.89655203E-03;
  COFD[1208] = -1.47223159E+01;
  COFD[1209] = 3.34934830E+00;
  COFD[1210] = -2.21831026E-01;
  COFD[1211] = 9.68019791E-03;
  COFD[1212] = -2.10732098E+01;
  COFD[1213] = 5.40974993E+00;
  COFD[1214] = -4.68080393E-01;
  COFD[1215] = 1.95128787E-02;
  COFD[1216] = -2.09174935E+01;
  COFD[1217] = 5.35659328E+00;
  COFD[1218] = -4.62684017E-01;
  COFD[1219] = 1.93375284E-02;
  COFD[1220] = -2.09174935E+01;
  COFD[1221] = 5.35659328E+00;
  COFD[1222] = -4.62684017E-01;
  COFD[1223] = 1.93375284E-02;
  COFD[1224] = -2.11618589E+01;
  COFD[1225] = 5.44451952E+00;
  COFD[1226] = -4.71804032E-01;
  COFD[1227] = 1.96436657E-02;
  COFD[1228] = -2.11558468E+01;
  COFD[1229] = 5.44451952E+00;
  COFD[1230] = -4.71804032E-01;
  COFD[1231] = 1.96436657E-02;
  COFD[1232] = -1.49037729E+01;
  COFD[1233] = 3.34934830E+00;
  COFD[1234] = -2.21831026E-01;
  COFD[1235] = 9.68019791E-03;
  COFD[1236] = -1.57442430E+01;
  COFD[1237] = 3.64418347E+00;
  COFD[1238] = -2.58820093E-01;
  COFD[1239] = 1.12296295E-02;
  COFD[1240] = -1.68027382E+01;
  COFD[1241] = 3.97899248E+00;
  COFD[1242] = -3.01178156E-01;
  COFD[1243] = 1.30211921E-02;
  COFD[1244] = -2.50221925E+01;
  COFD[1245] = 5.52564127E+00;
  COFD[1246] = -4.64221065E-01;
  COFD[1247] = 1.86139427E-02;
  COFD[1248] = -1.57070757E+01;
  COFD[1249] = 3.55327776E+00;
  COFD[1250] = -2.47179250E-01;
  COFD[1251] = 1.07315653E-02;
  COFD[1252] = -1.22583533E+01;
  COFD[1253] = 2.66488613E+00;
  COFD[1254] = -1.33839930E-01;
  COFD[1255] = 5.90524087E-03;
  COFD[1256] = -1.53791907E+01;
  COFD[1257] = 4.03496578E+00;
  COFD[1258] = -3.07489899E-01;
  COFD[1259] = 1.32554084E-02;
  COFD[1260] = -1.60662884E+01;
  COFD[1261] = 3.69744368E+00;
  COFD[1262] = -2.65754489E-01;
  COFD[1263] = 1.15311681E-02;
  COFD[1264] = -1.48358026E+01;
  COFD[1265] = 3.37160686E+00;
  COFD[1266] = -2.24698307E-01;
  COFD[1267] = 9.80357301E-03;
  COFD[1268] = -2.10691956E+01;
  COFD[1269] = 5.51604369E+00;
  COFD[1270] = -4.76337283E-01;
  COFD[1271] = 1.96544379E-02;
  COFD[1272] = -1.48503699E+01;
  COFD[1273] = 3.37160686E+00;
  COFD[1274] = -2.24698307E-01;
  COFD[1275] = 9.80357301E-03;
  COFD[1276] = -1.60758526E+01;
  COFD[1277] = 3.69744368E+00;
  COFD[1278] = -2.65754489E-01;
  COFD[1279] = 1.15311681E-02;
  COFD[1280] = -1.60711937E+01;
  COFD[1281] = 3.69744368E+00;
  COFD[1282] = -2.65754489E-01;
  COFD[1283] = 1.15311681E-02;
  COFD[1284] = -1.57272832E+01;
  COFD[1285] = 3.55756440E+00;
  COFD[1286] = -2.47689003E-01;
  COFD[1287] = 1.07517019E-02;
  COFD[1288] = -1.90110170E+01;
  COFD[1289] = 4.70842146E+00;
  COFD[1290] = -3.89924637E-01;
  COFD[1291] = 1.66320934E-02;
  COFD[1292] = -1.68880017E+01;
  COFD[1293] = 4.01759713E+00;
  COFD[1294] = -3.05560618E-01;
  COFD[1295] = 1.31852111E-02;
  COFD[1296] = -1.69251551E+01;
  COFD[1297] = 4.03028423E+00;
  COFD[1298] = -3.06976249E-01;
  COFD[1299] = 1.32370318E-02;
  COFD[1300] = -1.69075050E+01;
  COFD[1301] = 4.03028423E+00;
  COFD[1302] = -3.06976249E-01;
  COFD[1303] = 1.32370318E-02;
  COFD[1304] = -1.69075050E+01;
  COFD[1305] = 4.03028423E+00;
  COFD[1306] = -3.06976249E-01;
  COFD[1307] = 1.32370318E-02;
  COFD[1308] = -1.45499799E+01;
  COFD[1309] = 3.23154119E+00;
  COFD[1310] = -2.06656216E-01;
  COFD[1311] = 9.02762284E-03;
  COFD[1312] = -1.47832205E+01;
  COFD[1313] = 3.37160686E+00;
  COFD[1314] = -2.24698307E-01;
  COFD[1315] = 9.80357301E-03;
  COFD[1316] = -2.11276735E+01;
  COFD[1317] = 5.42949713E+00;
  COFD[1318] = -4.70215729E-01;
  COFD[1319] = 1.95888774E-02;
  COFD[1320] = -2.07138428E+01;
  COFD[1321] = 5.29812901E+00;
  COFD[1322] = -4.56644617E-01;
  COFD[1323] = 1.91364773E-02;
  COFD[1324] = -2.07138428E+01;
  COFD[1325] = 5.29812901E+00;
  COFD[1326] = -4.56644617E-01;
  COFD[1327] = 1.91364773E-02;
  COFD[1328] = -2.12140281E+01;
  COFD[1329] = 5.46309033E+00;
  COFD[1330] = -4.73777184E-01;
  COFD[1331] = 1.97122249E-02;
  COFD[1332] = -2.12082662E+01;
  COFD[1333] = 5.46309033E+00;
  COFD[1334] = -4.73777184E-01;
  COFD[1335] = 1.97122249E-02;
  COFD[1336] = -1.49583394E+01;
  COFD[1337] = 3.37160686E+00;
  COFD[1338] = -2.24698307E-01;
  COFD[1339] = 9.80357301E-03;
  COFD[1340] = -1.58288917E+01;
  COFD[1341] = 3.67755209E+00;
  COFD[1342] = -2.63170905E-01;
  COFD[1343] = 1.14190784E-02;
  COFD[1344] = -1.68409184E+01;
  COFD[1345] = 3.99442024E+00;
  COFD[1346] = -3.03003880E-01;
  COFD[1347] = 1.30931075E-02;
  COFD[1348] = -2.50224826E+01;
  COFD[1349] = 5.51888287E+00;
  COFD[1350] = -4.62631724E-01;
  COFD[1351] = 1.85178829E-02;
  COFD[1352] = -1.56842242E+01;
  COFD[1353] = 3.55327776E+00;
  COFD[1354] = -2.47179250E-01;
  COFD[1355] = 1.07315653E-02;
  COFD[1356] = -1.22541230E+01;
  COFD[1357] = 2.66488613E+00;
  COFD[1358] = -1.33839930E-01;
  COFD[1359] = 5.90524087E-03;
  COFD[1360] = -1.53769382E+01;
  COFD[1361] = 4.03496578E+00;
  COFD[1362] = -3.07489899E-01;
  COFD[1363] = 1.32554084E-02;
  COFD[1364] = -1.60424224E+01;
  COFD[1365] = 3.69744368E+00;
  COFD[1366] = -2.65754489E-01;
  COFD[1367] = 1.15311681E-02;
  COFD[1368] = -1.48176140E+01;
  COFD[1369] = 3.37160686E+00;
  COFD[1370] = -2.24698307E-01;
  COFD[1371] = 9.80357301E-03;
  COFD[1372] = -2.10499842E+01;
  COFD[1373] = 5.51604369E+00;
  COFD[1374] = -4.76337283E-01;
  COFD[1375] = 1.96544379E-02;
  COFD[1376] = -1.48316536E+01;
  COFD[1377] = 3.37160686E+00;
  COFD[1378] = -2.24698307E-01;
  COFD[1379] = 9.80357301E-03;
  COFD[1380] = -1.60515367E+01;
  COFD[1381] = 3.69744368E+00;
  COFD[1382] = -2.65754489E-01;
  COFD[1383] = 1.15311681E-02;
  COFD[1384] = -1.60470980E+01;
  COFD[1385] = 3.69744368E+00;
  COFD[1386] = -2.65754489E-01;
  COFD[1387] = 1.15311681E-02;
  COFD[1388] = -1.57044328E+01;
  COFD[1389] = 3.55756440E+00;
  COFD[1390] = -2.47689003E-01;
  COFD[1391] = 1.07517019E-02;
  COFD[1392] = -1.89849283E+01;
  COFD[1393] = 4.70842146E+00;
  COFD[1394] = -3.89924637E-01;
  COFD[1395] = 1.66320934E-02;
  COFD[1396] = -1.68697894E+01;
  COFD[1397] = 4.01759713E+00;
  COFD[1398] = -3.05560618E-01;
  COFD[1399] = 1.31852111E-02;
  COFD[1400] = -1.69075050E+01;
  COFD[1401] = 4.03028423E+00;
  COFD[1402] = -3.06976249E-01;
  COFD[1403] = 1.32370318E-02;
  COFD[1404] = -1.68904567E+01;
  COFD[1405] = 4.03028423E+00;
  COFD[1406] = -3.06976249E-01;
  COFD[1407] = 1.32370318E-02;
  COFD[1408] = -1.68904567E+01;
  COFD[1409] = 4.03028423E+00;
  COFD[1410] = -3.06976249E-01;
  COFD[1411] = 1.32370318E-02;
  COFD[1412] = -1.45342725E+01;
  COFD[1413] = 3.23154119E+00;
  COFD[1414] = -2.06656216E-01;
  COFD[1415] = 9.02762284E-03;
  COFD[1416] = -1.47668181E+01;
  COFD[1417] = 3.37160686E+00;
  COFD[1418] = -2.24698307E-01;
  COFD[1419] = 9.80357301E-03;
  COFD[1420] = -2.11037973E+01;
  COFD[1421] = 5.42949713E+00;
  COFD[1422] = -4.70215729E-01;
  COFD[1423] = 1.95888774E-02;
  COFD[1424] = -2.06902060E+01;
  COFD[1425] = 5.29812901E+00;
  COFD[1426] = -4.56644617E-01;
  COFD[1427] = 1.91364773E-02;
  COFD[1428] = -2.06902060E+01;
  COFD[1429] = 5.29812901E+00;
  COFD[1430] = -4.56644617E-01;
  COFD[1431] = 1.91364773E-02;
  COFD[1432] = -2.11906415E+01;
  COFD[1433] = 5.46309033E+00;
  COFD[1434] = -4.73777184E-01;
  COFD[1435] = 1.97122249E-02;
  COFD[1436] = -2.11851415E+01;
  COFD[1437] = 5.46309033E+00;
  COFD[1438] = -4.73777184E-01;
  COFD[1439] = 1.97122249E-02;
  COFD[1440] = -1.49352148E+01;
  COFD[1441] = 3.37160686E+00;
  COFD[1442] = -2.24698307E-01;
  COFD[1443] = 9.80357301E-03;
  COFD[1444] = -1.58092152E+01;
  COFD[1445] = 3.67755209E+00;
  COFD[1446] = -2.63170905E-01;
  COFD[1447] = 1.14190784E-02;
  COFD[1448] = -1.68170523E+01;
  COFD[1449] = 3.99442024E+00;
  COFD[1450] = -3.03003880E-01;
  COFD[1451] = 1.30931075E-02;
  COFD[1452] = -2.50224813E+01;
  COFD[1453] = 5.51888287E+00;
  COFD[1454] = -4.62631724E-01;
  COFD[1455] = 1.85178829E-02;
  COFD[1456] = -1.56842242E+01;
  COFD[1457] = 3.55327776E+00;
  COFD[1458] = -2.47179250E-01;
  COFD[1459] = 1.07315653E-02;
  COFD[1460] = -1.22541230E+01;
  COFD[1461] = 2.66488613E+00;
  COFD[1462] = -1.33839930E-01;
  COFD[1463] = 5.90524087E-03;
  COFD[1464] = -1.53769382E+01;
  COFD[1465] = 4.03496578E+00;
  COFD[1466] = -3.07489899E-01;
  COFD[1467] = 1.32554084E-02;
  COFD[1468] = -1.60424224E+01;
  COFD[1469] = 3.69744368E+00;
  COFD[1470] = -2.65754489E-01;
  COFD[1471] = 1.15311681E-02;
  COFD[1472] = -1.48176140E+01;
  COFD[1473] = 3.37160686E+00;
  COFD[1474] = -2.24698307E-01;
  COFD[1475] = 9.80357301E-03;
  COFD[1476] = -2.10499842E+01;
  COFD[1477] = 5.51604369E+00;
  COFD[1478] = -4.76337283E-01;
  COFD[1479] = 1.96544379E-02;
  COFD[1480] = -1.48316536E+01;
  COFD[1481] = 3.37160686E+00;
  COFD[1482] = -2.24698307E-01;
  COFD[1483] = 9.80357301E-03;
  COFD[1484] = -1.60515367E+01;
  COFD[1485] = 3.69744368E+00;
  COFD[1486] = -2.65754489E-01;
  COFD[1487] = 1.15311681E-02;
  COFD[1488] = -1.60470980E+01;
  COFD[1489] = 3.69744368E+00;
  COFD[1490] = -2.65754489E-01;
  COFD[1491] = 1.15311681E-02;
  COFD[1492] = -1.57044328E+01;
  COFD[1493] = 3.55756440E+00;
  COFD[1494] = -2.47689003E-01;
  COFD[1495] = 1.07517019E-02;
  COFD[1496] = -1.89849283E+01;
  COFD[1497] = 4.70842146E+00;
  COFD[1498] = -3.89924637E-01;
  COFD[1499] = 1.66320934E-02;
  COFD[1500] = -1.68697894E+01;
  COFD[1501] = 4.01759713E+00;
  COFD[1502] = -3.05560618E-01;
  COFD[1503] = 1.31852111E-02;
  COFD[1504] = -1.69075050E+01;
  COFD[1505] = 4.03028423E+00;
  COFD[1506] = -3.06976249E-01;
  COFD[1507] = 1.32370318E-02;
  COFD[1508] = -1.68904567E+01;
  COFD[1509] = 4.03028423E+00;
  COFD[1510] = -3.06976249E-01;
  COFD[1511] = 1.32370318E-02;
  COFD[1512] = -1.68904567E+01;
  COFD[1513] = 4.03028423E+00;
  COFD[1514] = -3.06976249E-01;
  COFD[1515] = 1.32370318E-02;
  COFD[1516] = -1.45342725E+01;
  COFD[1517] = 3.23154119E+00;
  COFD[1518] = -2.06656216E-01;
  COFD[1519] = 9.02762284E-03;
  COFD[1520] = -1.47668181E+01;
  COFD[1521] = 3.37160686E+00;
  COFD[1522] = -2.24698307E-01;
  COFD[1523] = 9.80357301E-03;
  COFD[1524] = -2.11037973E+01;
  COFD[1525] = 5.42949713E+00;
  COFD[1526] = -4.70215729E-01;
  COFD[1527] = 1.95888774E-02;
  COFD[1528] = -2.06902060E+01;
  COFD[1529] = 5.29812901E+00;
  COFD[1530] = -4.56644617E-01;
  COFD[1531] = 1.91364773E-02;
  COFD[1532] = -2.06902060E+01;
  COFD[1533] = 5.29812901E+00;
  COFD[1534] = -4.56644617E-01;
  COFD[1535] = 1.91364773E-02;
  COFD[1536] = -2.11906415E+01;
  COFD[1537] = 5.46309033E+00;
  COFD[1538] = -4.73777184E-01;
  COFD[1539] = 1.97122249E-02;
  COFD[1540] = -2.11851415E+01;
  COFD[1541] = 5.46309033E+00;
  COFD[1542] = -4.73777184E-01;
  COFD[1543] = 1.97122249E-02;
  COFD[1544] = -1.49352148E+01;
  COFD[1545] = 3.37160686E+00;
  COFD[1546] = -2.24698307E-01;
  COFD[1547] = 9.80357301E-03;
  COFD[1548] = -1.58092152E+01;
  COFD[1549] = 3.67755209E+00;
  COFD[1550] = -2.63170905E-01;
  COFD[1551] = 1.14190784E-02;
  COFD[1552] = -1.68170523E+01;
  COFD[1553] = 3.99442024E+00;
  COFD[1554] = -3.03003880E-01;
  COFD[1555] = 1.30931075E-02;
  COFD[1556] = -2.50224813E+01;
  COFD[1557] = 5.51888287E+00;
  COFD[1558] = -4.62631724E-01;
  COFD[1559] = 1.85178829E-02;
  COFD[1560] = -1.36591073E+01;
  COFD[1561] = 2.87846903E+00;
  COFD[1562] = -1.61681103E-01;
  COFD[1563] = 7.11635115E-03;
  COFD[1564] = -1.06186344E+01;
  COFD[1565] = 2.09311029E+00;
  COFD[1566] = -5.66715295E-02;
  COFD[1567] = 2.42607686E-03;
  COFD[1568] = -1.30215958E+01;
  COFD[1569] = 3.24041848E+00;
  COFD[1570] = -2.07804876E-01;
  COFD[1571] = 9.07722400E-03;
  COFD[1572] = -1.38694480E+01;
  COFD[1573] = 2.96800700E+00;
  COFD[1574] = -1.73235360E-01;
  COFD[1575] = 7.61377726E-03;
  COFD[1576] = -1.28155338E+01;
  COFD[1577] = 2.69967657E+00;
  COFD[1578] = -1.38362256E-01;
  COFD[1579] = 6.10143187E-03;
  COFD[1580] = -1.88701223E+01;
  COFD[1581] = 4.90944038E+00;
  COFD[1582] = -4.13866957E-01;
  COFD[1583] = 1.75844447E-02;
  COFD[1584] = -1.28284059E+01;
  COFD[1585] = 2.69967657E+00;
  COFD[1586] = -1.38362256E-01;
  COFD[1587] = 6.10143187E-03;
  COFD[1588] = -1.38776021E+01;
  COFD[1589] = 2.96800700E+00;
  COFD[1590] = -1.73235360E-01;
  COFD[1591] = 7.61377726E-03;
  COFD[1592] = -1.38736329E+01;
  COFD[1593] = 2.96800700E+00;
  COFD[1594] = -1.73235360E-01;
  COFD[1595] = 7.61377726E-03;
  COFD[1596] = -1.36800096E+01;
  COFD[1597] = 2.88311045E+00;
  COFD[1598] = -1.62274683E-01;
  COFD[1599] = 7.14167394E-03;
  COFD[1600] = -1.64069432E+01;
  COFD[1601] = 3.85900456E+00;
  COFD[1602] = -2.86199757E-01;
  COFD[1603] = 1.23960157E-02;
  COFD[1604] = -1.44868574E+01;
  COFD[1605] = 3.20806285E+00;
  COFD[1606] = -2.03619458E-01;
  COFD[1607] = 8.89655203E-03;
  COFD[1608] = -1.45499799E+01;
  COFD[1609] = 3.23154119E+00;
  COFD[1610] = -2.06656216E-01;
  COFD[1611] = 9.02762284E-03;
  COFD[1612] = -1.45342725E+01;
  COFD[1613] = 3.23154119E+00;
  COFD[1614] = -2.06656216E-01;
  COFD[1615] = 9.02762284E-03;
  COFD[1616] = -1.45342725E+01;
  COFD[1617] = 3.23154119E+00;
  COFD[1618] = -2.06656216E-01;
  COFD[1619] = 9.02762284E-03;
  COFD[1620] = -1.27188437E+01;
  COFD[1621] = 2.62486277E+00;
  COFD[1622] = -1.28774228E-01;
  COFD[1623] = 5.69134090E-03;
  COFD[1624] = -1.27687189E+01;
  COFD[1625] = 2.69967657E+00;
  COFD[1626] = -1.38362256E-01;
  COFD[1627] = 6.10143187E-03;
  COFD[1628] = -1.86387511E+01;
  COFD[1629] = 4.68388273E+00;
  COFD[1630] = -3.87039310E-01;
  COFD[1631] = 1.65190492E-02;
  COFD[1632] = -1.82026567E+01;
  COFD[1633] = 4.53059303E+00;
  COFD[1634] = -3.69227245E-01;
  COFD[1635] = 1.58309435E-02;
  COFD[1636] = -1.82026567E+01;
  COFD[1637] = 4.53059303E+00;
  COFD[1638] = -3.69227245E-01;
  COFD[1639] = 1.58309435E-02;
  COFD[1640] = -1.87435094E+01;
  COFD[1641] = 4.72753865E+00;
  COFD[1642] = -3.92270598E-01;
  COFD[1643] = 1.67283839E-02;
  COFD[1644] = -1.87385713E+01;
  COFD[1645] = 4.72753865E+00;
  COFD[1646] = -3.92270598E-01;
  COFD[1647] = 1.67283839E-02;
  COFD[1648] = -1.29223673E+01;
  COFD[1649] = 2.69967657E+00;
  COFD[1650] = -1.38362256E-01;
  COFD[1651] = 6.10143187E-03;
  COFD[1652] = -1.36610957E+01;
  COFD[1653] = 2.95689281E+00;
  COFD[1654] = -1.71804447E-01;
  COFD[1655] = 7.55232843E-03;
  COFD[1656] = -1.44057543E+01;
  COFD[1657] = 3.17964349E+00;
  COFD[1658] = -2.00113588E-01;
  COFD[1659] = 8.75259314E-03;
  COFD[1660] = -2.39197394E+01;
  COFD[1661] = 5.30539812E+00;
  COFD[1662] = -4.57330353E-01;
  COFD[1663] = 1.91560814E-02;
  COFD[1664] = -1.38076340E+01;
  COFD[1665] = 2.98558426E+00;
  COFD[1666] = -1.75507216E-01;
  COFD[1667] = 7.71173691E-03;
  COFD[1668] = -1.06123852E+01;
  COFD[1669] = 2.15849701E+00;
  COFD[1670] = -6.53886401E-02;
  COFD[1671] = 2.81453370E-03;
  COFD[1672] = -1.31792887E+01;
  COFD[1673] = 3.38003453E+00;
  COFD[1674] = -2.25783856E-01;
  COFD[1675] = 9.85028660E-03;
  COFD[1676] = -1.40154693E+01;
  COFD[1677] = 3.07458927E+00;
  COFD[1678] = -1.86899591E-01;
  COFD[1679] = 8.19829781E-03;
  COFD[1680] = -1.29335458E+01;
  COFD[1681] = 2.80841511E+00;
  COFD[1682] = -1.52629888E-01;
  COFD[1683] = 6.72604927E-03;
  COFD[1684] = -1.90524563E+01;
  COFD[1685] = 5.02608697E+00;
  COFD[1686] = -4.26959993E-01;
  COFD[1687] = 1.80709910E-02;
  COFD[1688] = -1.29470215E+01;
  COFD[1689] = 2.80841511E+00;
  COFD[1690] = -1.52629888E-01;
  COFD[1691] = 6.72604927E-03;
  COFD[1692] = -1.40241141E+01;
  COFD[1693] = 3.07458927E+00;
  COFD[1694] = -1.86899591E-01;
  COFD[1695] = 8.19829781E-03;
  COFD[1696] = -1.40199050E+01;
  COFD[1697] = 3.07458927E+00;
  COFD[1698] = -1.86899591E-01;
  COFD[1699] = 8.19829781E-03;
  COFD[1700] = -1.38327302E+01;
  COFD[1701] = 2.99164244E+00;
  COFD[1702] = -1.76293106E-01;
  COFD[1703] = 7.74575100E-03;
  COFD[1704] = -1.66339567E+01;
  COFD[1705] = 3.98859394E+00;
  COFD[1706] = -3.02316219E-01;
  COFD[1707] = 1.30661099E-02;
  COFD[1708] = -1.47223159E+01;
  COFD[1709] = 3.34934830E+00;
  COFD[1710] = -2.21831026E-01;
  COFD[1711] = 9.68019791E-03;
  COFD[1712] = -1.47832205E+01;
  COFD[1713] = 3.37160686E+00;
  COFD[1714] = -2.24698307E-01;
  COFD[1715] = 9.80357301E-03;
  COFD[1716] = -1.47668181E+01;
  COFD[1717] = 3.37160686E+00;
  COFD[1718] = -2.24698307E-01;
  COFD[1719] = 9.80357301E-03;
  COFD[1720] = -1.47668181E+01;
  COFD[1721] = 3.37160686E+00;
  COFD[1722] = -2.24698307E-01;
  COFD[1723] = 9.80357301E-03;
  COFD[1724] = -1.27687189E+01;
  COFD[1725] = 2.69967657E+00;
  COFD[1726] = -1.38362256E-01;
  COFD[1727] = 6.10143187E-03;
  COFD[1728] = -1.28846673E+01;
  COFD[1729] = 2.80841511E+00;
  COFD[1730] = -1.52629888E-01;
  COFD[1731] = 6.72604927E-03;
  COFD[1732] = -1.89453198E+01;
  COFD[1733] = 4.83982263E+00;
  COFD[1734] = -4.05783981E-01;
  COFD[1735] = 1.72721335E-02;
  COFD[1736] = -1.84167969E+01;
  COFD[1737] = 4.64907964E+00;
  COFD[1738] = -3.82960584E-01;
  COFD[1739] = 1.63598480E-02;
  COFD[1740] = -1.84167969E+01;
  COFD[1741] = 4.64907964E+00;
  COFD[1742] = -3.82960584E-01;
  COFD[1743] = 1.63598480E-02;
  COFD[1744] = -1.90409646E+01;
  COFD[1745] = 4.87977047E+00;
  COFD[1746] = -4.10448693E-01;
  COFD[1747] = 1.74535827E-02;
  COFD[1748] = -1.90357388E+01;
  COFD[1749] = 4.87977047E+00;
  COFD[1750] = -4.10448693E-01;
  COFD[1751] = 1.74535827E-02;
  COFD[1752] = -1.30459210E+01;
  COFD[1753] = 2.80841511E+00;
  COFD[1754] = -1.52629888E-01;
  COFD[1755] = 6.72604927E-03;
  COFD[1756] = -1.37952463E+01;
  COFD[1757] = 3.06403296E+00;
  COFD[1758] = -1.85549879E-01;
  COFD[1759] = 8.14070073E-03;
  COFD[1760] = -1.46029188E+01;
  COFD[1761] = 3.30683499E+00;
  COFD[1762] = -2.16378602E-01;
  COFD[1763] = 9.44670561E-03;
  COFD[1764] = -2.42503358E+01;
  COFD[1765] = 5.40793558E+00;
  COFD[1766] = -4.67888803E-01;
  COFD[1767] = 1.95062704E-02;
  COFD[1768] = -2.01091568E+01;
  COFD[1769] = 5.06508846E+00;
  COFD[1770] = -4.31649695E-01;
  COFD[1771] = 1.82598046E-02;
  COFD[1772] = -1.56250770E+01;
  COFD[1773] = 3.92232623E+00;
  COFD[1774] = -2.94183188E-01;
  COFD[1775] = 1.27325506E-02;
  COFD[1776] = -1.94386101E+01;
  COFD[1777] = 5.43670321E+00;
  COFD[1778] = -4.70985776E-01;
  COFD[1779] = 1.96158515E-02;
  COFD[1780] = -2.04422993E+01;
  COFD[1781] = 5.18126516E+00;
  COFD[1782] = -4.45006101E-01;
  COFD[1783] = 1.87713513E-02;
  COFD[1784] = -1.90163591E+01;
  COFD[1785] = 4.83982263E+00;
  COFD[1786] = -4.05783981E-01;
  COFD[1787] = 1.72721335E-02;
  COFD[1788] = -1.89495722E+01;
  COFD[1789] = 4.10779436E+00;
  COFD[1790] = -2.35955221E-01;
  COFD[1791] = 7.03410627E-03;
  COFD[1792] = -1.90365260E+01;
  COFD[1793] = 4.83982263E+00;
  COFD[1794] = -4.05783981E-01;
  COFD[1795] = 1.72721335E-02;
  COFD[1796] = -2.04573512E+01;
  COFD[1797] = 5.18126516E+00;
  COFD[1798] = -4.45006101E-01;
  COFD[1799] = 1.87713513E-02;
  COFD[1800] = -2.04499985E+01;
  COFD[1801] = 5.18126516E+00;
  COFD[1802] = -4.45006101E-01;
  COFD[1803] = 1.87713513E-02;
  COFD[1804] = -2.01404432E+01;
  COFD[1805] = 5.07318708E+00;
  COFD[1806] = -4.32608485E-01;
  COFD[1807] = 1.82977801E-02;
  COFD[1808] = -2.21011224E+01;
  COFD[1809] = 5.52816937E+00;
  COFD[1810] = -4.65357386E-01;
  COFD[1811] = 1.86939300E-02;
  COFD[1812] = -2.10732098E+01;
  COFD[1813] = 5.40974993E+00;
  COFD[1814] = -4.68080393E-01;
  COFD[1815] = 1.95128787E-02;
  COFD[1816] = -2.11276735E+01;
  COFD[1817] = 5.42949713E+00;
  COFD[1818] = -4.70215729E-01;
  COFD[1819] = 1.95888774E-02;
  COFD[1820] = -2.11037973E+01;
  COFD[1821] = 5.42949713E+00;
  COFD[1822] = -4.70215729E-01;
  COFD[1823] = 1.95888774E-02;
  COFD[1824] = -2.11037973E+01;
  COFD[1825] = 5.42949713E+00;
  COFD[1826] = -4.70215729E-01;
  COFD[1827] = 1.95888774E-02;
  COFD[1828] = -1.86387511E+01;
  COFD[1829] = 4.68388273E+00;
  COFD[1830] = -3.87039310E-01;
  COFD[1831] = 1.65190492E-02;
  COFD[1832] = -1.89453198E+01;
  COFD[1833] = 4.83982263E+00;
  COFD[1834] = -4.05783981E-01;
  COFD[1835] = 1.72721335E-02;
  COFD[1836] = -2.03539177E+01;
  COFD[1837] = 4.56892447E+00;
  COFD[1838] = -3.06147917E-01;
  COFD[1839] = 1.04644656E-02;
  COFD[1840] = -2.10007650E+01;
  COFD[1841] = 4.89469279E+00;
  COFD[1842] = -3.56802629E-01;
  COFD[1843] = 1.29773618E-02;
  COFD[1844] = -2.10007650E+01;
  COFD[1845] = 4.89469279E+00;
  COFD[1846] = -3.56802629E-01;
  COFD[1847] = 1.29773618E-02;
  COFD[1848] = -2.01405949E+01;
  COFD[1849] = 4.47715839E+00;
  COFD[1850] = -2.92174559E-01;
  COFD[1851] = 9.78054024E-03;
  COFD[1852] = -2.01317080E+01;
  COFD[1853] = 4.47715839E+00;
  COFD[1854] = -2.92174559E-01;
  COFD[1855] = 9.78054024E-03;
  COFD[1856] = -1.91941481E+01;
  COFD[1857] = 4.83982263E+00;
  COFD[1858] = -4.05783981E-01;
  COFD[1859] = 1.72721335E-02;
  COFD[1860] = -2.01735995E+01;
  COFD[1861] = 5.17081088E+00;
  COFD[1862] = -4.43888632E-01;
  COFD[1863] = 1.87323735E-02;
  COFD[1864] = -2.10549431E+01;
  COFD[1865] = 5.37765260E+00;
  COFD[1866] = -4.64755221E-01;
  COFD[1867] = 1.94011994E-02;
  COFD[1868] = -1.94498128E+01;
  COFD[1869] = 2.70960748E+00;
  COFD[1870] = -2.96143895E-02;
  COFD[1871] = -2.82061574E-03;
  COFD[1872] = -1.98937690E+01;
  COFD[1873] = 4.98725744E+00;
  COFD[1874] = -4.22747564E-01;
  COFD[1875] = 1.79211872E-02;
  COFD[1876] = -1.55102782E+01;
  COFD[1877] = 3.86994367E+00;
  COFD[1878] = -2.87570272E-01;
  COFD[1879] = 1.24533966E-02;
  COFD[1880] = -1.90281731E+01;
  COFD[1881] = 5.30278138E+00;
  COFD[1882] = -4.57079018E-01;
  COFD[1883] = 1.91486428E-02;
  COFD[1884] = -2.02220865E+01;
  COFD[1885] = 5.10234461E+00;
  COFD[1886] = -4.36039243E-01;
  COFD[1887] = 1.84326773E-02;
  COFD[1888] = -1.84871270E+01;
  COFD[1889] = 4.64907964E+00;
  COFD[1890] = -3.82960584E-01;
  COFD[1891] = 1.63598480E-02;
  COFD[1892] = -1.75035239E+01;
  COFD[1893] = 3.41230274E+00;
  COFD[1894] = -1.32406923E-01;
  COFD[1895] = 2.10435795E-03;
  COFD[1896] = -1.85070738E+01;
  COFD[1897] = 4.64907964E+00;
  COFD[1898] = -3.82960584E-01;
  COFD[1899] = 1.63598480E-02;
  COFD[1900] = -1.99478621E+01;
  COFD[1901] = 5.00532670E+00;
  COFD[1902] = -4.24686880E-01;
  COFD[1903] = 1.79891229E-02;
  COFD[1904] = -1.99406295E+01;
  COFD[1905] = 5.00532670E+00;
  COFD[1906] = -4.24686880E-01;
  COFD[1907] = 1.79891229E-02;
  COFD[1908] = -1.99343637E+01;
  COFD[1909] = 4.99746599E+00;
  COFD[1910] = -4.23838967E-01;
  COFD[1911] = 1.79592096E-02;
  COFD[1912] = -2.20276430E+01;
  COFD[1913] = 5.52121111E+00;
  COFD[1914] = -4.67224743E-01;
  COFD[1915] = 1.88799020E-02;
  COFD[1916] = -2.09174935E+01;
  COFD[1917] = 5.35659328E+00;
  COFD[1918] = -4.62684017E-01;
  COFD[1919] = 1.93375284E-02;
  COFD[1920] = -2.07138428E+01;
  COFD[1921] = 5.29812901E+00;
  COFD[1922] = -4.56644617E-01;
  COFD[1923] = 1.91364773E-02;
  COFD[1924] = -2.06902060E+01;
  COFD[1925] = 5.29812901E+00;
  COFD[1926] = -4.56644617E-01;
  COFD[1927] = 1.91364773E-02;
  COFD[1928] = -2.06902060E+01;
  COFD[1929] = 5.29812901E+00;
  COFD[1930] = -4.56644617E-01;
  COFD[1931] = 1.91364773E-02;
  COFD[1932] = -1.82026567E+01;
  COFD[1933] = 4.53059303E+00;
  COFD[1934] = -3.69227245E-01;
  COFD[1935] = 1.58309435E-02;
  COFD[1936] = -1.84167969E+01;
  COFD[1937] = 4.64907964E+00;
  COFD[1938] = -3.82960584E-01;
  COFD[1939] = 1.63598480E-02;
  COFD[1940] = -2.10007650E+01;
  COFD[1941] = 4.89469279E+00;
  COFD[1942] = -3.56802629E-01;
  COFD[1943] = 1.29773618E-02;
  COFD[1944] = -2.02181133E+01;
  COFD[1945] = 4.55977855E+00;
  COFD[1946] = -3.10217483E-01;
  COFD[1947] = 1.08560899E-02;
  COFD[1948] = -2.02181133E+01;
  COFD[1949] = 4.55977855E+00;
  COFD[1950] = -3.10217483E-01;
  COFD[1951] = 1.08560899E-02;
  COFD[1952] = -2.08446140E+01;
  COFD[1953] = 4.82681292E+00;
  COFD[1954] = -3.46134170E-01;
  COFD[1955] = 1.24444585E-02;
  COFD[1956] = -2.08358634E+01;
  COFD[1957] = 4.82681292E+00;
  COFD[1958] = -3.46134170E-01;
  COFD[1959] = 1.24444585E-02;
  COFD[1960] = -1.92184914E+01;
  COFD[1961] = 4.83981923E+00;
  COFD[1962] = -4.05783571E-01;
  COFD[1963] = 1.72721170E-02;
  COFD[1964] = -1.98853844E+01;
  COFD[1965] = 5.06690169E+00;
  COFD[1966] = -4.31863731E-01;
  COFD[1967] = 1.82682527E-02;
  COFD[1968] = -2.08606487E+01;
  COFD[1969] = 5.31248163E+00;
  COFD[1970] = -4.58014127E-01;
  COFD[1971] = 1.91765075E-02;
  COFD[1972] = -2.05464448E+01;
  COFD[1973] = 3.22046026E+00;
  COFD[1974] = -1.04593424E-01;
  COFD[1975] = 7.40271888E-04;
  COFD[1976] = -1.98937690E+01;
  COFD[1977] = 4.98725744E+00;
  COFD[1978] = -4.22747564E-01;
  COFD[1979] = 1.79211872E-02;
  COFD[1980] = -1.55102782E+01;
  COFD[1981] = 3.86994367E+00;
  COFD[1982] = -2.87570272E-01;
  COFD[1983] = 1.24533966E-02;
  COFD[1984] = -1.90281731E+01;
  COFD[1985] = 5.30278138E+00;
  COFD[1986] = -4.57079018E-01;
  COFD[1987] = 1.91486428E-02;
  COFD[1988] = -2.02220865E+01;
  COFD[1989] = 5.10234461E+00;
  COFD[1990] = -4.36039243E-01;
  COFD[1991] = 1.84326773E-02;
  COFD[1992] = -1.84871270E+01;
  COFD[1993] = 4.64907964E+00;
  COFD[1994] = -3.82960584E-01;
  COFD[1995] = 1.63598480E-02;
  COFD[1996] = -1.75035239E+01;
  COFD[1997] = 3.41230274E+00;
  COFD[1998] = -1.32406923E-01;
  COFD[1999] = 2.10435795E-03;
  COFD[2000] = -1.85070738E+01;
  COFD[2001] = 4.64907964E+00;
  COFD[2002] = -3.82960584E-01;
  COFD[2003] = 1.63598480E-02;
  COFD[2004] = -1.99478621E+01;
  COFD[2005] = 5.00532670E+00;
  COFD[2006] = -4.24686880E-01;
  COFD[2007] = 1.79891229E-02;
  COFD[2008] = -1.99406295E+01;
  COFD[2009] = 5.00532670E+00;
  COFD[2010] = -4.24686880E-01;
  COFD[2011] = 1.79891229E-02;
  COFD[2012] = -1.99343637E+01;
  COFD[2013] = 4.99746599E+00;
  COFD[2014] = -4.23838967E-01;
  COFD[2015] = 1.79592096E-02;
  COFD[2016] = -2.20276430E+01;
  COFD[2017] = 5.52121111E+00;
  COFD[2018] = -4.67224743E-01;
  COFD[2019] = 1.88799020E-02;
  COFD[2020] = -2.09174935E+01;
  COFD[2021] = 5.35659328E+00;
  COFD[2022] = -4.62684017E-01;
  COFD[2023] = 1.93375284E-02;
  COFD[2024] = -2.07138428E+01;
  COFD[2025] = 5.29812901E+00;
  COFD[2026] = -4.56644617E-01;
  COFD[2027] = 1.91364773E-02;
  COFD[2028] = -2.06902060E+01;
  COFD[2029] = 5.29812901E+00;
  COFD[2030] = -4.56644617E-01;
  COFD[2031] = 1.91364773E-02;
  COFD[2032] = -2.06902060E+01;
  COFD[2033] = 5.29812901E+00;
  COFD[2034] = -4.56644617E-01;
  COFD[2035] = 1.91364773E-02;
  COFD[2036] = -1.82026567E+01;
  COFD[2037] = 4.53059303E+00;
  COFD[2038] = -3.69227245E-01;
  COFD[2039] = 1.58309435E-02;
  COFD[2040] = -1.84167969E+01;
  COFD[2041] = 4.64907964E+00;
  COFD[2042] = -3.82960584E-01;
  COFD[2043] = 1.63598480E-02;
  COFD[2044] = -2.10007650E+01;
  COFD[2045] = 4.89469279E+00;
  COFD[2046] = -3.56802629E-01;
  COFD[2047] = 1.29773618E-02;
  COFD[2048] = -2.02181133E+01;
  COFD[2049] = 4.55977855E+00;
  COFD[2050] = -3.10217483E-01;
  COFD[2051] = 1.08560899E-02;
  COFD[2052] = -2.02181133E+01;
  COFD[2053] = 4.55977855E+00;
  COFD[2054] = -3.10217483E-01;
  COFD[2055] = 1.08560899E-02;
  COFD[2056] = -2.08446140E+01;
  COFD[2057] = 4.82681292E+00;
  COFD[2058] = -3.46134170E-01;
  COFD[2059] = 1.24444585E-02;
  COFD[2060] = -2.08358634E+01;
  COFD[2061] = 4.82681292E+00;
  COFD[2062] = -3.46134170E-01;
  COFD[2063] = 1.24444585E-02;
  COFD[2064] = -1.92184914E+01;
  COFD[2065] = 4.83981923E+00;
  COFD[2066] = -4.05783571E-01;
  COFD[2067] = 1.72721170E-02;
  COFD[2068] = -1.98853844E+01;
  COFD[2069] = 5.06690169E+00;
  COFD[2070] = -4.31863731E-01;
  COFD[2071] = 1.82682527E-02;
  COFD[2072] = -2.08606487E+01;
  COFD[2073] = 5.31248163E+00;
  COFD[2074] = -4.58014127E-01;
  COFD[2075] = 1.91765075E-02;
  COFD[2076] = -2.05464448E+01;
  COFD[2077] = 3.22046026E+00;
  COFD[2078] = -1.04593424E-01;
  COFD[2079] = 7.40271888E-04;
  COFD[2080] = -2.02136057E+01;
  COFD[2081] = 5.10993120E+00;
  COFD[2082] = -4.36931630E-01;
  COFD[2083] = 1.84677592E-02;
  COFD[2084] = -1.57172268E+01;
  COFD[2085] = 3.96062263E+00;
  COFD[2086] = -2.98964970E-01;
  COFD[2087] = 1.29322565E-02;
  COFD[2088] = -1.95318517E+01;
  COFD[2089] = 5.47046983E+00;
  COFD[2090] = -4.74577605E-01;
  COFD[2091] = 1.97408822E-02;
  COFD[2092] = -2.05196205E+01;
  COFD[2093] = 5.21473296E+00;
  COFD[2094] = -4.48646311E-01;
  COFD[2095] = 1.89013813E-02;
  COFD[2096] = -1.91105534E+01;
  COFD[2097] = 4.87977047E+00;
  COFD[2098] = -4.10448693E-01;
  COFD[2099] = 1.74535827E-02;
  COFD[2100] = -1.87236049E+01;
  COFD[2101] = 4.00967621E+00;
  COFD[2102] = -2.21153539E-01;
  COFD[2103] = 6.31528745E-03;
  COFD[2104] = -1.91302705E+01;
  COFD[2105] = 4.87977047E+00;
  COFD[2106] = -4.10448693E-01;
  COFD[2107] = 1.74535827E-02;
  COFD[2108] = -2.05341767E+01;
  COFD[2109] = 5.21473296E+00;
  COFD[2110] = -4.48646311E-01;
  COFD[2111] = 1.89013813E-02;
  COFD[2112] = -2.05270679E+01;
  COFD[2113] = 5.21473296E+00;
  COFD[2114] = -4.48646311E-01;
  COFD[2115] = 1.89013813E-02;
  COFD[2116] = -2.02444913E+01;
  COFD[2117] = 5.11785645E+00;
  COFD[2118] = -4.37867828E-01;
  COFD[2119] = 1.85047543E-02;
  COFD[2120] = -2.20860863E+01;
  COFD[2121] = 5.52171573E+00;
  COFD[2122] = -4.63284984E-01;
  COFD[2123] = 1.85570924E-02;
  COFD[2124] = -2.11618589E+01;
  COFD[2125] = 5.44451952E+00;
  COFD[2126] = -4.71804032E-01;
  COFD[2127] = 1.96436657E-02;
  COFD[2128] = -2.12140281E+01;
  COFD[2129] = 5.46309033E+00;
  COFD[2130] = -4.73777184E-01;
  COFD[2131] = 1.97122249E-02;
  COFD[2132] = -2.11906415E+01;
  COFD[2133] = 5.46309033E+00;
  COFD[2134] = -4.73777184E-01;
  COFD[2135] = 1.97122249E-02;
  COFD[2136] = -2.11906415E+01;
  COFD[2137] = 5.46309033E+00;
  COFD[2138] = -4.73777184E-01;
  COFD[2139] = 1.97122249E-02;
  COFD[2140] = -1.87435094E+01;
  COFD[2141] = 4.72753865E+00;
  COFD[2142] = -3.92270598E-01;
  COFD[2143] = 1.67283839E-02;
  COFD[2144] = -1.90409646E+01;
  COFD[2145] = 4.87977047E+00;
  COFD[2146] = -4.10448693E-01;
  COFD[2147] = 1.74535827E-02;
  COFD[2148] = -2.01405949E+01;
  COFD[2149] = 4.47715839E+00;
  COFD[2150] = -2.92174559E-01;
  COFD[2151] = 9.78054024E-03;
  COFD[2152] = -2.08446140E+01;
  COFD[2153] = 4.82681292E+00;
  COFD[2154] = -3.46134170E-01;
  COFD[2155] = 1.24444585E-02;
  COFD[2156] = -2.08446140E+01;
  COFD[2157] = 4.82681292E+00;
  COFD[2158] = -3.46134170E-01;
  COFD[2159] = 1.24444585E-02;
  COFD[2160] = -1.99154422E+01;
  COFD[2161] = 4.38041133E+00;
  COFD[2162] = -2.77538214E-01;
  COFD[2163] = 9.06748822E-03;
  COFD[2164] = -1.99068325E+01;
  COFD[2165] = 4.38041133E+00;
  COFD[2166] = -2.77538214E-01;
  COFD[2167] = 9.06748822E-03;
  COFD[2168] = -1.92836946E+01;
  COFD[2169] = 4.87977047E+00;
  COFD[2170] = -4.10448693E-01;
  COFD[2171] = 1.74535827E-02;
  COFD[2172] = -2.02538837E+01;
  COFD[2173] = 5.20405648E+00;
  COFD[2174] = -4.47505457E-01;
  COFD[2175] = 1.88615428E-02;
  COFD[2176] = -2.11253124E+01;
  COFD[2177] = 5.40762714E+00;
  COFD[2178] = -4.67856822E-01;
  COFD[2179] = 1.95051950E-02;
  COFD[2180] = -1.91495340E+01;
  COFD[2181] = 2.57271534E+00;
  COFD[2182] = -9.78488237E-03;
  COFD[2183] = -3.75334124E-03;
  COFD[2184] = -2.02052920E+01;
  COFD[2185] = 5.10993120E+00;
  COFD[2186] = -4.36931630E-01;
  COFD[2187] = 1.84677592E-02;
  COFD[2188] = -1.57161352E+01;
  COFD[2189] = 3.96062263E+00;
  COFD[2190] = -2.98964970E-01;
  COFD[2191] = 1.29322565E-02;
  COFD[2192] = -1.95312879E+01;
  COFD[2193] = 5.47046983E+00;
  COFD[2194] = -4.74577605E-01;
  COFD[2195] = 1.97408822E-02;
  COFD[2196] = -2.05107395E+01;
  COFD[2197] = 5.21473296E+00;
  COFD[2198] = -4.48646311E-01;
  COFD[2199] = 1.89013813E-02;
  COFD[2200] = -1.91045520E+01;
  COFD[2201] = 4.87977047E+00;
  COFD[2202] = -4.10448693E-01;
  COFD[2203] = 1.74535827E-02;
  COFD[2204] = -1.87171338E+01;
  COFD[2205] = 4.00967621E+00;
  COFD[2206] = -2.21153539E-01;
  COFD[2207] = 6.31528745E-03;
  COFD[2208] = -1.91240292E+01;
  COFD[2209] = 4.87977047E+00;
  COFD[2210] = -4.10448693E-01;
  COFD[2211] = 1.74535827E-02;
  COFD[2212] = -2.05250357E+01;
  COFD[2213] = 5.21473296E+00;
  COFD[2214] = -4.48646311E-01;
  COFD[2215] = 1.89013813E-02;
  COFD[2216] = -2.05180548E+01;
  COFD[2217] = 5.21473296E+00;
  COFD[2218] = -4.48646311E-01;
  COFD[2219] = 1.89013813E-02;
  COFD[2220] = -2.02361782E+01;
  COFD[2221] = 5.11785645E+00;
  COFD[2222] = -4.37867828E-01;
  COFD[2223] = 1.85047543E-02;
  COFD[2224] = -2.20758670E+01;
  COFD[2225] = 5.52171573E+00;
  COFD[2226] = -4.63284984E-01;
  COFD[2227] = 1.85570924E-02;
  COFD[2228] = -2.11558468E+01;
  COFD[2229] = 5.44451952E+00;
  COFD[2230] = -4.71804032E-01;
  COFD[2231] = 1.96436657E-02;
  COFD[2232] = -2.12082662E+01;
  COFD[2233] = 5.46309033E+00;
  COFD[2234] = -4.73777184E-01;
  COFD[2235] = 1.97122249E-02;
  COFD[2236] = -2.11851415E+01;
  COFD[2237] = 5.46309033E+00;
  COFD[2238] = -4.73777184E-01;
  COFD[2239] = 1.97122249E-02;
  COFD[2240] = -2.11851415E+01;
  COFD[2241] = 5.46309033E+00;
  COFD[2242] = -4.73777184E-01;
  COFD[2243] = 1.97122249E-02;
  COFD[2244] = -1.87385713E+01;
  COFD[2245] = 4.72753865E+00;
  COFD[2246] = -3.92270598E-01;
  COFD[2247] = 1.67283839E-02;
  COFD[2248] = -1.90357388E+01;
  COFD[2249] = 4.87977047E+00;
  COFD[2250] = -4.10448693E-01;
  COFD[2251] = 1.74535827E-02;
  COFD[2252] = -2.01317080E+01;
  COFD[2253] = 4.47715839E+00;
  COFD[2254] = -2.92174559E-01;
  COFD[2255] = 9.78054024E-03;
  COFD[2256] = -2.08358634E+01;
  COFD[2257] = 4.82681292E+00;
  COFD[2258] = -3.46134170E-01;
  COFD[2259] = 1.24444585E-02;
  COFD[2260] = -2.08358634E+01;
  COFD[2261] = 4.82681292E+00;
  COFD[2262] = -3.46134170E-01;
  COFD[2263] = 1.24444585E-02;
  COFD[2264] = -1.99068325E+01;
  COFD[2265] = 4.38041133E+00;
  COFD[2266] = -2.77538214E-01;
  COFD[2267] = 9.06748822E-03;
  COFD[2268] = -1.98983685E+01;
  COFD[2269] = 4.38041133E+00;
  COFD[2270] = -2.77538214E-01;
  COFD[2271] = 9.06748822E-03;
  COFD[2272] = -1.92752307E+01;
  COFD[2273] = 4.87977047E+00;
  COFD[2274] = -4.10448693E-01;
  COFD[2275] = 1.74535827E-02;
  COFD[2276] = -2.02471926E+01;
  COFD[2277] = 5.20405648E+00;
  COFD[2278] = -4.47505457E-01;
  COFD[2279] = 1.88615428E-02;
  COFD[2280] = -2.11164313E+01;
  COFD[2281] = 5.40762714E+00;
  COFD[2282] = -4.67856822E-01;
  COFD[2283] = 1.95051950E-02;
  COFD[2284] = -1.91495337E+01;
  COFD[2285] = 2.57271534E+00;
  COFD[2286] = -9.78488237E-03;
  COFD[2287] = -3.75334124E-03;
  COFD[2288] = -1.40437647E+01;
  COFD[2289] = 2.98558426E+00;
  COFD[2290] = -1.75507216E-01;
  COFD[2291] = 7.71173691E-03;
  COFD[2292] = -1.06507866E+01;
  COFD[2293] = 2.15849701E+00;
  COFD[2294] = -6.53886401E-02;
  COFD[2295] = 2.81453370E-03;
  COFD[2296] = -1.31995018E+01;
  COFD[2297] = 3.38003453E+00;
  COFD[2298] = -2.25783856E-01;
  COFD[2299] = 9.85028660E-03;
  COFD[2300] = -1.42641696E+01;
  COFD[2301] = 3.07458927E+00;
  COFD[2302] = -1.86899591E-01;
  COFD[2303] = 8.19829781E-03;
  COFD[2304] = -1.31147338E+01;
  COFD[2305] = 2.80841511E+00;
  COFD[2306] = -1.52629888E-01;
  COFD[2307] = 6.72604927E-03;
  COFD[2308] = -2.04933546E+01;
  COFD[2309] = 5.41191809E+00;
  COFD[2310] = -4.68316503E-01;
  COFD[2311] = 1.95213604E-02;
  COFD[2312] = -1.31342108E+01;
  COFD[2313] = 2.80841511E+00;
  COFD[2314] = -1.52629888E-01;
  COFD[2315] = 6.72604927E-03;
  COFD[2316] = -1.42784657E+01;
  COFD[2317] = 3.07458927E+00;
  COFD[2318] = -1.86899591E-01;
  COFD[2319] = 8.19829781E-03;
  COFD[2320] = -1.42714849E+01;
  COFD[2321] = 3.07458927E+00;
  COFD[2322] = -1.86899591E-01;
  COFD[2323] = 8.19829781E-03;
  COFD[2324] = -1.40688472E+01;
  COFD[2325] = 2.99164244E+00;
  COFD[2326] = -1.76293106E-01;
  COFD[2327] = 7.74575100E-03;
  COFD[2328] = -1.69110608E+01;
  COFD[2329] = 3.98859394E+00;
  COFD[2330] = -3.02316219E-01;
  COFD[2331] = 1.30661099E-02;
  COFD[2332] = -1.49037729E+01;
  COFD[2333] = 3.34934830E+00;
  COFD[2334] = -2.21831026E-01;
  COFD[2335] = 9.68019791E-03;
  COFD[2336] = -1.49583394E+01;
  COFD[2337] = 3.37160686E+00;
  COFD[2338] = -2.24698307E-01;
  COFD[2339] = 9.80357301E-03;
  COFD[2340] = -1.49352148E+01;
  COFD[2341] = 3.37160686E+00;
  COFD[2342] = -2.24698307E-01;
  COFD[2343] = 9.80357301E-03;
  COFD[2344] = -1.49352148E+01;
  COFD[2345] = 3.37160686E+00;
  COFD[2346] = -2.24698307E-01;
  COFD[2347] = 9.80357301E-03;
  COFD[2348] = -1.29223673E+01;
  COFD[2349] = 2.69967657E+00;
  COFD[2350] = -1.38362256E-01;
  COFD[2351] = 6.10143187E-03;
  COFD[2352] = -1.30459210E+01;
  COFD[2353] = 2.80841511E+00;
  COFD[2354] = -1.52629888E-01;
  COFD[2355] = 6.72604927E-03;
  COFD[2356] = -1.91941481E+01;
  COFD[2357] = 4.83982263E+00;
  COFD[2358] = -4.05783981E-01;
  COFD[2359] = 1.72721335E-02;
  COFD[2360] = -1.92184914E+01;
  COFD[2361] = 4.83981923E+00;
  COFD[2362] = -4.05783571E-01;
  COFD[2363] = 1.72721170E-02;
  COFD[2364] = -1.92184914E+01;
  COFD[2365] = 4.83981923E+00;
  COFD[2366] = -4.05783571E-01;
  COFD[2367] = 1.72721170E-02;
  COFD[2368] = -1.92836946E+01;
  COFD[2369] = 4.87977047E+00;
  COFD[2370] = -4.10448693E-01;
  COFD[2371] = 1.74535827E-02;
  COFD[2372] = -1.92752307E+01;
  COFD[2373] = 4.87977047E+00;
  COFD[2374] = -4.10448693E-01;
  COFD[2375] = 1.74535827E-02;
  COFD[2376] = -1.32854111E+01;
  COFD[2377] = 2.80841511E+00;
  COFD[2378] = -1.52629888E-01;
  COFD[2379] = 6.72604927E-03;
  COFD[2380] = -1.39934893E+01;
  COFD[2381] = 3.06403296E+00;
  COFD[2382] = -1.85549879E-01;
  COFD[2383] = 8.14070073E-03;
  COFD[2384] = -1.48516207E+01;
  COFD[2385] = 3.30683499E+00;
  COFD[2386] = -2.16378602E-01;
  COFD[2387] = 9.44670561E-03;
  COFD[2388] = -2.42503475E+01;
  COFD[2389] = 5.40793558E+00;
  COFD[2390] = -4.67888803E-01;
  COFD[2391] = 1.95062704E-02;
  COFD[2392] = -1.47449891E+01;
  COFD[2393] = 3.24099818E+00;
  COFD[2394] = -2.07879868E-01;
  COFD[2395] = 9.08046161E-03;
  COFD[2396] = -1.14322126E+01;
  COFD[2397] = 2.42102995E+00;
  COFD[2398] = -1.01354401E-01;
  COFD[2399] = 4.45948069E-03;
  COFD[2400] = -1.41953189E+01;
  COFD[2401] = 3.68995707E+00;
  COFD[2402] = -2.64784442E-01;
  COFD[2403] = 1.14891794E-02;
  COFD[2404] = -1.50366933E+01;
  COFD[2405] = 3.35941694E+00;
  COFD[2406] = -2.23126552E-01;
  COFD[2407] = 9.73586922E-03;
  COFD[2408] = -1.38538393E+01;
  COFD[2409] = 3.06403296E+00;
  COFD[2410] = -1.85549879E-01;
  COFD[2411] = 8.14070073E-03;
  COFD[2412] = -2.07061954E+01;
  COFD[2413] = 5.46880564E+00;
  COFD[2414] = -4.74397826E-01;
  COFD[2415] = 1.97344800E-02;
  COFD[2416] = -1.38702007E+01;
  COFD[2417] = 3.06403296E+00;
  COFD[2418] = -1.85549879E-01;
  COFD[2419] = 8.14070073E-03;
  COFD[2420] = -1.50478662E+01;
  COFD[2421] = 3.35941694E+00;
  COFD[2422] = -2.23126552E-01;
  COFD[2423] = 9.73586922E-03;
  COFD[2424] = -1.50424192E+01;
  COFD[2425] = 3.35941694E+00;
  COFD[2426] = -2.23126552E-01;
  COFD[2427] = 9.73586922E-03;
  COFD[2428] = -1.47729612E+01;
  COFD[2429] = 3.24835477E+00;
  COFD[2430] = -2.08830747E-01;
  COFD[2431] = 9.12148470E-03;
  COFD[2432] = -1.78202286E+01;
  COFD[2433] = 4.31993263E+00;
  COFD[2434] = -3.42851630E-01;
  COFD[2435] = 1.47243122E-02;
  COFD[2436] = -1.57442430E+01;
  COFD[2437] = 3.64418347E+00;
  COFD[2438] = -2.58820093E-01;
  COFD[2439] = 1.12296295E-02;
  COFD[2440] = -1.58288917E+01;
  COFD[2441] = 3.67755209E+00;
  COFD[2442] = -2.63170905E-01;
  COFD[2443] = 1.14190784E-02;
  COFD[2444] = -1.58092152E+01;
  COFD[2445] = 3.67755209E+00;
  COFD[2446] = -2.63170905E-01;
  COFD[2447] = 1.14190784E-02;
  COFD[2448] = -1.58092152E+01;
  COFD[2449] = 3.67755209E+00;
  COFD[2450] = -2.63170905E-01;
  COFD[2451] = 1.14190784E-02;
  COFD[2452] = -1.36610957E+01;
  COFD[2453] = 2.95689281E+00;
  COFD[2454] = -1.71804447E-01;
  COFD[2455] = 7.55232843E-03;
  COFD[2456] = -1.37952463E+01;
  COFD[2457] = 3.06403296E+00;
  COFD[2458] = -1.85549879E-01;
  COFD[2459] = 8.14070073E-03;
  COFD[2460] = -2.01735995E+01;
  COFD[2461] = 5.17081088E+00;
  COFD[2462] = -4.43888632E-01;
  COFD[2463] = 1.87323735E-02;
  COFD[2464] = -1.98853844E+01;
  COFD[2465] = 5.06690169E+00;
  COFD[2466] = -4.31863731E-01;
  COFD[2467] = 1.82682527E-02;
  COFD[2468] = -1.98853844E+01;
  COFD[2469] = 5.06690169E+00;
  COFD[2470] = -4.31863731E-01;
  COFD[2471] = 1.82682527E-02;
  COFD[2472] = -2.02538837E+01;
  COFD[2473] = 5.20405648E+00;
  COFD[2474] = -4.47505457E-01;
  COFD[2475] = 1.88615428E-02;
  COFD[2476] = -2.02471926E+01;
  COFD[2477] = 5.20405648E+00;
  COFD[2478] = -4.47505457E-01;
  COFD[2479] = 1.88615428E-02;
  COFD[2480] = -1.39934893E+01;
  COFD[2481] = 3.06403296E+00;
  COFD[2482] = -1.85549879E-01;
  COFD[2483] = 8.14070073E-03;
  COFD[2484] = -1.47917050E+01;
  COFD[2485] = 3.34581072E+00;
  COFD[2486] = -2.21377649E-01;
  COFD[2487] = 9.66079817E-03;
  COFD[2488] = -1.56182178E+01;
  COFD[2489] = 3.58311501E+00;
  COFD[2490] = -2.50860151E-01;
  COFD[2491] = 1.08830877E-02;
  COFD[2492] = -2.47385674E+01;
  COFD[2493] = 5.51399489E+00;
  COFD[2494] = -4.72948445E-01;
  COFD[2495] = 1.93856424E-02;
  COFD[2496] = -1.56749318E+01;
  COFD[2497] = 3.51447210E+00;
  COFD[2498] = -2.42579007E-01;
  COFD[2499] = 1.05506318E-02;
  COFD[2500] = -1.20579435E+01;
  COFD[2501] = 2.63303536E+00;
  COFD[2502] = -1.29792632E-01;
  COFD[2503] = 5.73363738E-03;
  COFD[2504] = -1.51177800E+01;
  COFD[2505] = 3.99904647E+00;
  COFD[2506] = -3.03517220E-01;
  COFD[2507] = 1.31117363E-02;
  COFD[2508] = -1.59068122E+01;
  COFD[2509] = 3.60186887E+00;
  COFD[2510] = -2.53302622E-01;
  COFD[2511] = 1.09893496E-02;
  COFD[2512] = -1.46739282E+01;
  COFD[2513] = 3.30683499E+00;
  COFD[2514] = -2.16378602E-01;
  COFD[2515] = 9.44670561E-03;
  COFD[2516] = -2.12711805E+01;
  COFD[2517] = 5.53137765E+00;
  COFD[2518] = -4.73118467E-01;
  COFD[2519] = 1.93125852E-02;
  COFD[2520] = -1.46940858E+01;
  COFD[2521] = 3.30683499E+00;
  COFD[2522] = -2.16378602E-01;
  COFD[2523] = 9.44670561E-03;
  COFD[2524] = -1.59218538E+01;
  COFD[2525] = 3.60186887E+00;
  COFD[2526] = -2.53302622E-01;
  COFD[2527] = 1.09893496E-02;
  COFD[2528] = -1.59145062E+01;
  COFD[2529] = 3.60186887E+00;
  COFD[2530] = -2.53302622E-01;
  COFD[2531] = 1.09893496E-02;
  COFD[2532] = -1.56952963E+01;
  COFD[2533] = 3.51861272E+00;
  COFD[2534] = -2.43068621E-01;
  COFD[2535] = 1.05698368E-02;
  COFD[2536] = -1.89571393E+01;
  COFD[2537] = 4.64763677E+00;
  COFD[2538] = -3.82799418E-01;
  COFD[2539] = 1.63539171E-02;
  COFD[2540] = -1.68027382E+01;
  COFD[2541] = 3.97899248E+00;
  COFD[2542] = -3.01178156E-01;
  COFD[2543] = 1.30211921E-02;
  COFD[2544] = -1.68409184E+01;
  COFD[2545] = 3.99442024E+00;
  COFD[2546] = -3.03003880E-01;
  COFD[2547] = 1.30931075E-02;
  COFD[2548] = -1.68170523E+01;
  COFD[2549] = 3.99442024E+00;
  COFD[2550] = -3.03003880E-01;
  COFD[2551] = 1.30931075E-02;
  COFD[2552] = -1.68170523E+01;
  COFD[2553] = 3.99442024E+00;
  COFD[2554] = -3.03003880E-01;
  COFD[2555] = 1.30931075E-02;
  COFD[2556] = -1.44057543E+01;
  COFD[2557] = 3.17964349E+00;
  COFD[2558] = -2.00113588E-01;
  COFD[2559] = 8.75259314E-03;
  COFD[2560] = -1.46029188E+01;
  COFD[2561] = 3.30683499E+00;
  COFD[2562] = -2.16378602E-01;
  COFD[2563] = 9.44670561E-03;
  COFD[2564] = -2.10549431E+01;
  COFD[2565] = 5.37765260E+00;
  COFD[2566] = -4.64755221E-01;
  COFD[2567] = 1.94011994E-02;
  COFD[2568] = -2.08606487E+01;
  COFD[2569] = 5.31248163E+00;
  COFD[2570] = -4.58014127E-01;
  COFD[2571] = 1.91765075E-02;
  COFD[2572] = -2.08606487E+01;
  COFD[2573] = 5.31248163E+00;
  COFD[2574] = -4.58014127E-01;
  COFD[2575] = 1.91765075E-02;
  COFD[2576] = -2.11253124E+01;
  COFD[2577] = 5.40762714E+00;
  COFD[2578] = -4.67856822E-01;
  COFD[2579] = 1.95051950E-02;
  COFD[2580] = -2.11164313E+01;
  COFD[2581] = 5.40762714E+00;
  COFD[2582] = -4.67856822E-01;
  COFD[2583] = 1.95051950E-02;
  COFD[2584] = -1.48516207E+01;
  COFD[2585] = 3.30683499E+00;
  COFD[2586] = -2.16378602E-01;
  COFD[2587] = 9.44670561E-03;
  COFD[2588] = -1.56182178E+01;
  COFD[2589] = 3.58311501E+00;
  COFD[2590] = -2.50860151E-01;
  COFD[2591] = 1.08830877E-02;
  COFD[2592] = -1.67835296E+01;
  COFD[2593] = 3.94346012E+00;
  COFD[2594] = -2.96835271E-01;
  COFD[2595] = 1.28438696E-02;
  COFD[2596] = -2.49944573E+01;
  COFD[2597] = 5.52779193E+00;
  COFD[2598] = -4.65764754E-01;
  COFD[2599] = 1.87288224E-02;
  COFD[2600] = -2.46746425E+01;
  COFD[2601] = 5.51584151E+00;
  COFD[2602] = -4.76108968E-01;
  COFD[2603] = 1.96365975E-02;
  COFD[2604] = -2.18596191E+01;
  COFD[2605] = 4.62023553E+00;
  COFD[2606] = -3.79733248E-01;
  COFD[2607] = 1.62407913E-02;
  COFD[2608] = -2.50118194E+01;
  COFD[2609] = 5.51533435E+00;
  COFD[2610] = -4.61886312E-01;
  COFD[2611] = 1.84746261E-02;
  COFD[2612] = -2.47604202E+01;
  COFD[2613] = 5.51868821E+00;
  COFD[2614] = -4.73225401E-01;
  COFD[2615] = 1.93849523E-02;
  COFD[2616] = -2.42503398E+01;
  COFD[2617] = 5.40793558E+00;
  COFD[2618] = -4.67888803E-01;
  COFD[2619] = 1.95062704E-02;
  COFD[2620] = -1.78173107E+01;
  COFD[2621] = 1.97328753E+00;
  COFD[2622] = 7.64115525E-02;
  COFD[2623] = -7.78429580E-03;
  COFD[2624] = -2.42503408E+01;
  COFD[2625] = 5.40793558E+00;
  COFD[2626] = -4.67888803E-01;
  COFD[2627] = 1.95062704E-02;
  COFD[2628] = -2.47604207E+01;
  COFD[2629] = 5.51868821E+00;
  COFD[2630] = -4.73225401E-01;
  COFD[2631] = 1.93849523E-02;
  COFD[2632] = -2.47604204E+01;
  COFD[2633] = 5.51868821E+00;
  COFD[2634] = -4.73225401E-01;
  COFD[2635] = 1.93849523E-02;
  COFD[2636] = -2.46791388E+01;
  COFD[2637] = 5.51562323E+00;
  COFD[2638] = -4.75876136E-01;
  COFD[2639] = 1.96184735E-02;
  COFD[2640] = -2.38843968E+01;
  COFD[2641] = 4.82941387E+00;
  COFD[2642] = -3.46543093E-01;
  COFD[2643] = 1.24648862E-02;
  COFD[2644] = -2.50221925E+01;
  COFD[2645] = 5.52564127E+00;
  COFD[2646] = -4.64221065E-01;
  COFD[2647] = 1.86139427E-02;
  COFD[2648] = -2.50224826E+01;
  COFD[2649] = 5.51888287E+00;
  COFD[2650] = -4.62631724E-01;
  COFD[2651] = 1.85178829E-02;
  COFD[2652] = -2.50224813E+01;
  COFD[2653] = 5.51888287E+00;
  COFD[2654] = -4.62631724E-01;
  COFD[2655] = 1.85178829E-02;
  COFD[2656] = -2.50224813E+01;
  COFD[2657] = 5.51888287E+00;
  COFD[2658] = -4.62631724E-01;
  COFD[2659] = 1.85178829E-02;
  COFD[2660] = -2.39197394E+01;
  COFD[2661] = 5.30539812E+00;
  COFD[2662] = -4.57330353E-01;
  COFD[2663] = 1.91560814E-02;
  COFD[2664] = -2.42503358E+01;
  COFD[2665] = 5.40793558E+00;
  COFD[2666] = -4.67888803E-01;
  COFD[2667] = 1.95062704E-02;
  COFD[2668] = -1.94498128E+01;
  COFD[2669] = 2.70960748E+00;
  COFD[2670] = -2.96143895E-02;
  COFD[2671] = -2.82061574E-03;
  COFD[2672] = -2.05464448E+01;
  COFD[2673] = 3.22046026E+00;
  COFD[2674] = -1.04593424E-01;
  COFD[2675] = 7.40271888E-04;
  COFD[2676] = -2.05464448E+01;
  COFD[2677] = 3.22046026E+00;
  COFD[2678] = -1.04593424E-01;
  COFD[2679] = 7.40271888E-04;
  COFD[2680] = -1.91495340E+01;
  COFD[2681] = 2.57271534E+00;
  COFD[2682] = -9.78488237E-03;
  COFD[2683] = -3.75334124E-03;
  COFD[2684] = -1.91495337E+01;
  COFD[2685] = 2.57271534E+00;
  COFD[2686] = -9.78488237E-03;
  COFD[2687] = -3.75334124E-03;
  COFD[2688] = -2.42503475E+01;
  COFD[2689] = 5.40793558E+00;
  COFD[2690] = -4.67888803E-01;
  COFD[2691] = 1.95062704E-02;
  COFD[2692] = -2.47385674E+01;
  COFD[2693] = 5.51399489E+00;
  COFD[2694] = -4.72948445E-01;
  COFD[2695] = 1.93856424E-02;
  COFD[2696] = -2.49944573E+01;
  COFD[2697] = 5.52779193E+00;
  COFD[2698] = -4.65764754E-01;
  COFD[2699] = 1.87288224E-02;
  COFD[2700] = -1.47803548E+01;
  COFD[2701] = 1.70835496E-01;
  COFD[2702] = 3.31563490E-01;
  COFD[2703] = -1.95427081E-02;
}

// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKTDIF(int* KTDIF)
{
  KTDIF[0] = 1;
  KTDIF[1] = 2;
  KTDIF[2] = 25;
}

// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFTD(amrex::Real* COFTD)
{
  COFTD[0] = 4.45261429E-01;
  COFTD[1] = 4.94696577E-05;
  COFTD[2] = -2.63023124E-08;
  COFTD[3] = 4.90305626E-12;
  COFTD[4] = 0.00000000E+00;
  COFTD[5] = 0.00000000E+00;
  COFTD[6] = 0.00000000E+00;
  COFTD[7] = 0.00000000E+00;
  COFTD[8] = -1.52534742E-01;
  COFTD[9] = -5.46404022E-05;
  COFTD[10] = 2.93412470E-08;
  COFTD[11] = -4.87091914E-12;
  COFTD[12] = 4.42736017E-01;
  COFTD[13] = 7.11765887E-05;
  COFTD[14] = -3.84765396E-08;
  COFTD[15] = 6.86318682E-12;
  COFTD[16] = 4.15577509E-01;
  COFTD[17] = 1.09736860E-05;
  COFTD[18] = -3.96016410E-09;
  COFTD[19] = 1.14412838E-12;
  COFTD[20] = 6.02020725E-02;
  COFTD[21] = 5.61555766E-04;
  COFTD[22] = -2.55370076E-07;
  COFTD[23] = 3.63385938E-11;
  COFTD[24] = 4.21927217E-01;
  COFTD[25] = 1.11413555E-05;
  COFTD[26] = -4.02067239E-09;
  COFTD[27] = 1.16160979E-12;
  COFTD[28] = 4.46067449E-01;
  COFTD[29] = 7.17121674E-05;
  COFTD[30] = -3.87660620E-08;
  COFTD[31] = 6.91482989E-12;
  COFTD[32] = 4.44449675E-01;
  COFTD[33] = 7.14520855E-05;
  COFTD[34] = -3.86254673E-08;
  COFTD[35] = 6.88975155E-12;
  COFTD[36] = 4.44650439E-01;
  COFTD[37] = 5.06628083E-05;
  COFTD[38] = -2.69818971E-08;
  COFTD[39] = 5.01286176E-12;
  COFTD[40] = 3.25740914E-01;
  COFTD[41] = 3.03631979E-04;
  COFTD[42] = -1.55289598E-07;
  COFTD[43] = 2.41465297E-11;
  COFTD[44] = 3.58655012E-01;
  COFTD[45] = 1.21524029E-04;
  COFTD[46] = -6.53752231E-08;
  COFTD[47] = 1.09001717E-11;
  COFTD[48] = 3.50272421E-01;
  COFTD[49] = 1.23583915E-04;
  COFTD[50] = -6.63969060E-08;
  COFTD[51] = 1.10354999E-11;
  COFTD[52] = 3.43456512E-01;
  COFTD[53] = 1.21179111E-04;
  COFTD[54] = -6.51048965E-08;
  COFTD[55] = 1.08207614E-11;
  COFTD[56] = 3.43456512E-01;
  COFTD[57] = 1.21179111E-04;
  COFTD[58] = -6.51048965E-08;
  COFTD[59] = 1.08207614E-11;
  COFTD[60] = 3.88529645E-01;
  COFTD[61] = -4.91887876E-06;
  COFTD[62] = 5.69539773E-09;
  COFTD[63] = -5.68222508E-13;
  COFTD[64] = 3.91826338E-01;
  COFTD[65] = 1.03465156E-05;
  COFTD[66] = -3.73383198E-09;
  COFTD[67] = 1.07873892E-12;
  COFTD[68] = 1.70978021E-01;
  COFTD[69] = 4.98264784E-04;
  COFTD[70] = -2.37738756E-07;
  COFTD[71] = 3.49395808E-11;
  COFTD[72] = 1.79391459E-01;
  COFTD[73] = 4.84490026E-04;
  COFTD[74] = -2.32134571E-07;
  COFTD[75] = 3.42178255E-11;
  COFTD[76] = 1.79391459E-01;
  COFTD[77] = 4.84490026E-04;
  COFTD[78] = -2.32134571E-07;
  COFTD[79] = 3.42178255E-11;
  COFTD[80] = 1.62149345E-01;
  COFTD[81] = 5.03434613E-04;
  COFTD[82] = -2.39391657E-07;
  COFTD[83] = 3.50983128E-11;
  COFTD[84] = 1.61391212E-01;
  COFTD[85] = 5.01080793E-04;
  COFTD[86] = -2.38272376E-07;
  COFTD[87] = 3.49342098E-11;
  COFTD[88] = 4.65846649E-01;
  COFTD[89] = 1.23010864E-05;
  COFTD[90] = -4.43919397E-09;
  COFTD[91] = 1.28252458E-12;
  COFTD[92] = 4.07197076E-01;
  COFTD[93] = 6.29805002E-05;
  COFTD[94] = -3.40117037E-08;
  COFTD[95] = 6.09207140E-12;
  COFTD[96] = 4.12032475E-01;
  COFTD[97] = 1.28802637E-04;
  COFTD[98] = -6.94527177E-08;
  COFTD[99] = 1.16546716E-11;
  COFTD[100] = -4.83575595E-02;
  COFTD[101] = -7.29902876E-04;
  COFTD[102] = 3.28156716E-07;
  COFTD[103] = -4.63569261E-11;
  COFTD[104] = 2.40744276E-01;
  COFTD[105] = 4.45343183E-04;
  COFTD[106] = -2.18173742E-07;
  COFTD[107] = 3.26958310E-11;
  COFTD[108] = 1.52534742E-01;
  COFTD[109] = 5.46404022E-05;
  COFTD[110] = -2.93412470E-08;
  COFTD[111] = 4.87091914E-12;
  COFTD[112] = 0.00000000E+00;
  COFTD[113] = 0.00000000E+00;
  COFTD[114] = 0.00000000E+00;
  COFTD[115] = 0.00000000E+00;
  COFTD[116] = 2.20482081E-01;
  COFTD[117] = 4.80162630E-04;
  COFTD[118] = -2.32927140E-07;
  COFTD[119] = 3.46469239E-11;
  COFTD[120] = 2.70008280E-01;
  COFTD[121] = 3.61552588E-04;
  COFTD[122] = -1.80743499E-07;
  COFTD[123] = 2.75319341E-11;
  COFTD[124] = -1.41882969E-01;
  COFTD[125] = 7.66554625E-04;
  COFTD[126] = -3.06548330E-07;
  COFTD[127] = 4.02957302E-11;
  COFTD[128] = 2.72039997E-01;
  COFTD[129] = 3.64273144E-04;
  COFTD[130] = -1.82103531E-07;
  COFTD[131] = 2.77391022E-11;
  COFTD[132] = 2.21307723E-01;
  COFTD[133] = 4.81960701E-04;
  COFTD[134] = -2.33799385E-07;
  COFTD[135] = 3.47766667E-11;
  COFTD[136] = 2.20907136E-01;
  COFTD[137] = 4.81088308E-04;
  COFTD[138] = -2.33376187E-07;
  COFTD[139] = 3.47137178E-11;
  COFTD[140] = 2.39409087E-01;
  COFTD[141] = 4.47195587E-04;
  COFTD[142] = -2.18950923E-07;
  COFTD[143] = 3.27972343E-11;
  COFTD[144] = 2.44368810E-02;
  COFTD[145] = 7.18240808E-04;
  COFTD[146] = -3.19717751E-07;
  COFTD[147] = 4.48827629E-11;
  COFTD[148] = 1.45433118E-01;
  COFTD[149] = 5.30347483E-04;
  COFTD[150] = -2.50118653E-07;
  COFTD[151] = 3.64613999E-11;
  COFTD[152] = 1.40109849E-01;
  COFTD[153] = 5.30876682E-04;
  COFTD[154] = -2.49898194E-07;
  COFTD[155] = 3.63822692E-11;
  COFTD[156] = 1.38759794E-01;
  COFTD[157] = 5.25761319E-04;
  COFTD[158] = -2.47490252E-07;
  COFTD[159] = 3.60317009E-11;
  COFTD[160] = 1.38759794E-01;
  COFTD[161] = 5.25761319E-04;
  COFTD[162] = -2.47490252E-07;
  COFTD[163] = 3.60317009E-11;
  COFTD[164] = 2.81033340E-01;
  COFTD[165] = 3.13353760E-04;
  COFTD[166] = -1.58475652E-07;
  COFTD[167] = 2.43852063E-11;
  COFTD[168] = 2.62294123E-01;
  COFTD[169] = 3.51223003E-04;
  COFTD[170] = -1.75579644E-07;
  COFTD[171] = 2.67453446E-11;
  COFTD[172] = -1.19542204E-01;
  COFTD[173] = 8.00169896E-04;
  COFTD[174] = -3.27980898E-07;
  COFTD[175] = 4.37706567E-11;
  COFTD[176] = -9.21021062E-02;
  COFTD[177] = 7.88085864E-04;
  COFTD[178] = -3.29283183E-07;
  COFTD[179] = 4.44614948E-11;
  COFTD[180] = -9.21021062E-02;
  COFTD[181] = 7.88085864E-04;
  COFTD[182] = -3.29283183E-07;
  COFTD[183] = 4.44614948E-11;
  COFTD[184] = -1.24938974E-01;
  COFTD[185] = 7.98385053E-04;
  COFTD[186] = -3.25757471E-07;
  COFTD[187] = 4.33526419E-11;
  COFTD[188] = -1.24647578E-01;
  COFTD[189] = 7.96522972E-04;
  COFTD[190] = -3.24997704E-07;
  COFTD[191] = 4.32515301E-11;
  COFTD[192] = 2.85751979E-01;
  COFTD[193] = 3.82634072E-04;
  COFTD[194] = -1.91282329E-07;
  COFTD[195] = 2.91372719E-11;
  COFTD[196] = 2.13713591E-01;
  COFTD[197] = 4.56510862E-04;
  COFTD[198] = -2.21708032E-07;
  COFTD[199] = 3.30064357E-11;
  COFTD[200] = 1.63366239E-01;
  COFTD[201] = 5.54240001E-04;
  COFTD[202] = -2.62337753E-07;
  COFTD[203] = 3.83388181E-11;
  COFTD[204] = 2.14687863E-01;
  COFTD[205] = -8.33408992E-04;
  COFTD[206] = 3.11551008E-07;
  COFTD[207] = -3.92419743E-11;
  COFTD[208] = -1.58494270E-01;
  COFTD[209] = 8.57382603E-04;
  COFTD[210] = -3.42931140E-07;
  COFTD[211] = 4.50831409E-11;
  COFTD[212] = 4.83575595E-02;
  COFTD[213] = 7.29902876E-04;
  COFTD[214] = -3.28156716E-07;
  COFTD[215] = 4.63569261E-11;
  COFTD[216] = -2.14687863E-01;
  COFTD[217] = 8.33408992E-04;
  COFTD[218] = -3.11551008E-07;
  COFTD[219] = 3.92419743E-11;
  COFTD[220] = -1.74367468E-01;
  COFTD[221] = 8.56399915E-04;
  COFTD[222] = -3.37494379E-07;
  COFTD[223] = 4.39622394E-11;
  COFTD[224] = -1.22051466E-01;
  COFTD[225] = 8.50463321E-04;
  COFTD[226] = -3.49880208E-07;
  COFTD[227] = 4.67982394E-11;
  COFTD[228] = -1.82074509E-01;
  COFTD[229] = 3.58494859E-04;
  COFTD[230] = -4.14797572E-08;
  COFTD[231] = -1.91123311E-12;
  COFTD[232] = -1.22051962E-01;
  COFTD[233] = 8.50466778E-04;
  COFTD[234] = -3.49881630E-07;
  COFTD[235] = 4.67984296E-11;
  COFTD[236] = -1.74367822E-01;
  COFTD[237] = 8.56401656E-04;
  COFTD[238] = -3.37495065E-07;
  COFTD[239] = 4.39623287E-11;
  COFTD[240] = -1.74367651E-01;
  COFTD[241] = 8.56400812E-04;
  COFTD[242] = -3.37494732E-07;
  COFTD[243] = 4.39622854E-11;
  COFTD[244] = -1.59489983E-01;
  COFTD[245] = 8.57405115E-04;
  COFTD[246] = -3.42640753E-07;
  COFTD[247] = 4.50206824E-11;
  COFTD[248] = -2.48073493E-01;
  COFTD[249] = 7.25180496E-04;
  COFTD[250] = -2.36929722E-07;
  COFTD[251] = 2.71093741E-11;
  COFTD[252] = -2.12064065E-01;
  COFTD[253] = 8.37240651E-04;
  COFTD[254] = -3.14517782E-07;
  COFTD[255] = 3.97420162E-11;
  COFTD[256] = -2.14134032E-01;
  COFTD[257] = 8.35097175E-04;
  COFTD[258] = -3.12605803E-07;
  COFTD[259] = 3.94096885E-11;
  COFTD[260] = -2.14132909E-01;
  COFTD[261] = 8.35092796E-04;
  COFTD[262] = -3.12604164E-07;
  COFTD[263] = 3.94094818E-11;
  COFTD[264] = -2.14132909E-01;
  COFTD[265] = 8.35092796E-04;
  COFTD[266] = -3.12604164E-07;
  COFTD[267] = 3.94094818E-11;
  COFTD[268] = -9.90513879E-02;
  COFTD[269] = 8.41283392E-04;
  COFTD[270] = -3.51342290E-07;
  COFTD[271] = 4.74261956E-11;
  COFTD[272] = -1.22049550E-01;
  COFTD[273] = 8.50449971E-04;
  COFTD[274] = -3.49874716E-07;
  COFTD[275] = 4.67975048E-11;
  COFTD[276] = -2.08192097E-01;
  COFTD[277] = 4.49774819E-04;
  COFTD[278] = -8.69551363E-08;
  COFTD[279] = 4.64889109E-12;
  COFTD[280] = -2.25054527E-01;
  COFTD[281] = 5.20024585E-04;
  COFTD[282] = -1.23031218E-07;
  COFTD[283] = 9.92192373E-12;
  COFTD[284] = -2.25054527E-01;
  COFTD[285] = 5.20024585E-04;
  COFTD[286] = -1.23031218E-07;
  COFTD[287] = 9.92192373E-12;
  COFTD[288] = -2.03684878E-01;
  COFTD[289] = 4.32867875E-04;
  COFTD[290] = -7.84216895E-08;
  COFTD[291] = 3.41085687E-12;
  COFTD[292] = -2.03684620E-01;
  COFTD[293] = 4.32867326E-04;
  COFTD[294] = -7.84215900E-08;
  COFTD[295] = 3.41085254E-12;
  COFTD[296] = -1.22055221E-01;
  COFTD[297] = 8.50489487E-04;
  COFTD[298] = -3.49890973E-07;
  COFTD[299] = 4.67996793E-11;
  COFTD[300] = -1.72580282E-01;
  COFTD[301] = 8.56644825E-04;
  COFTD[302] = -3.38190010E-07;
  COFTD[303] = 4.41015047E-11;
  COFTD[304] = -2.07905569E-01;
  COFTD[305] = 8.41099546E-04;
  COFTD[306] = -3.18089807E-07;
  COFTD[307] = 4.03679127E-11;
  COFTD[308] = 0.00000000E+00;
  COFTD[309] = 0.00000000E+00;
  COFTD[310] = 0.00000000E+00;
  COFTD[311] = 0.00000000E+00;
}

// compute the critical parameters for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS(
  amrex::Real* Tci, amrex::Real* ai, amrex::Real* bi, amrex::Real* acentric_i)
{

  amrex::Real EPS[26];
  amrex::Real SIG[26];
  amrex::Real wt[26];
  amrex::Real Rcst = 83.144598; // in bar [CGS] !
  amrex::Real avogadro = 6.02214199e23;
  amrex::Real boltzmann = 1.3806503e-16; // we work in CGS

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

  // species 1: H2
  // Imported from NIST
  Tci[1] = 33.145000;
  ai[1] = 1e6 * 0.42748 * Rcst * Rcst * Tci[1] * Tci[1] /
          (2.015880 * 2.015880 * 12.964000);
  bi[1] = 0.08664 * Rcst * Tci[1] / (2.015880 * 12.964000);
  acentric_i[1] = -0.219000;

  // species 2: H
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

  // species 4: O
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

  // species 6: OH
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

  // species 8: HO2
  Tci[8] = 1.316 * EPS[8];
  ai[8] = (5.55 * avogadro * avogadro * EPS[8] * boltzmann * 1e-24 * SIG[8] *
           SIG[8] * SIG[8]) /
          (wt[8] * wt[8]);
  bi[8] = 0.855 * avogadro * 1e-24 * SIG[8] * SIG[8] * SIG[8] / (wt[8]);
  acentric_i[8] = 0.0;

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

  // species 11: CH4
  // Imported from NIST
  Tci[11] = 190.560000;
  ai[11] = 1e6 * 0.42748 * Rcst * Rcst * Tci[11] * Tci[11] /
           (16.043030 * 16.043030 * 45.990000);
  bi[11] = 0.08664 * Rcst * Tci[11] / (16.043030 * 45.990000);
  acentric_i[11] = 0.011000;

  // species 12: CH3
  Tci[12] = 1.316 * EPS[12];
  ai[12] = (5.55 * avogadro * avogadro * EPS[12] * boltzmann * 1e-24 * SIG[12] *
            SIG[12] * SIG[12]) /
           (wt[12] * wt[12]);
  bi[12] = 0.855 * avogadro * 1e-24 * SIG[12] * SIG[12] * SIG[12] / (wt[12]);
  acentric_i[12] = 0.0;

  // species 13: CH2
  Tci[13] = 1.316 * EPS[13];
  ai[13] = (5.55 * avogadro * avogadro * EPS[13] * boltzmann * 1e-24 * SIG[13] *
            SIG[13] * SIG[13]) /
           (wt[13] * wt[13]);
  bi[13] = 0.855 * avogadro * 1e-24 * SIG[13] * SIG[13] * SIG[13] / (wt[13]);
  acentric_i[13] = 0.0;

  // species 14: CH2(S)
  Tci[14] = 1.316 * EPS[14];
  ai[14] = (5.55 * avogadro * avogadro * EPS[14] * boltzmann * 1e-24 * SIG[14] *
            SIG[14] * SIG[14]) /
           (wt[14] * wt[14]);
  bi[14] = 0.855 * avogadro * 1e-24 * SIG[14] * SIG[14] * SIG[14] / (wt[14]);
  acentric_i[14] = 0.0;

  // species 15: C
  Tci[15] = 1.316 * EPS[15];
  ai[15] = (5.55 * avogadro * avogadro * EPS[15] * boltzmann * 1e-24 * SIG[15] *
            SIG[15] * SIG[15]) /
           (wt[15] * wt[15]);
  bi[15] = 0.855 * avogadro * 1e-24 * SIG[15] * SIG[15] * SIG[15] / (wt[15]);
  acentric_i[15] = 0.0;

  // species 16: CH
  Tci[16] = 1.316 * EPS[16];
  ai[16] = (5.55 * avogadro * avogadro * EPS[16] * boltzmann * 1e-24 * SIG[16] *
            SIG[16] * SIG[16]) /
           (wt[16] * wt[16]);
  bi[16] = 0.855 * avogadro * 1e-24 * SIG[16] * SIG[16] * SIG[16] / (wt[16]);
  acentric_i[16] = 0.0;

  // species 17: CH3OH
  // Imported from NIST
  Tci[17] = 512.640000;
  ai[17] = 1e6 * 0.42748 * Rcst * Rcst * Tci[17] * Tci[17] /
           (32.042000 * 32.042000 * 80.970000);
  bi[17] = 0.08664 * Rcst * Tci[17] / (32.042000 * 80.970000);
  acentric_i[17] = 0.565000;

  // species 18: CH3O
  Tci[18] = 1.316 * EPS[18];
  ai[18] = (5.55 * avogadro * avogadro * EPS[18] * boltzmann * 1e-24 * SIG[18] *
            SIG[18] * SIG[18]) /
           (wt[18] * wt[18]);
  bi[18] = 0.855 * avogadro * 1e-24 * SIG[18] * SIG[18] * SIG[18] / (wt[18]);
  acentric_i[18] = 0.0;

  // species 19: CH2OH
  Tci[19] = 1.316 * EPS[19];
  ai[19] = (5.55 * avogadro * avogadro * EPS[19] * boltzmann * 1e-24 * SIG[19] *
            SIG[19] * SIG[19]) /
           (wt[19] * wt[19]);
  bi[19] = 0.855 * avogadro * 1e-24 * SIG[19] * SIG[19] * SIG[19] / (wt[19]);
  acentric_i[19] = 0.0;

  // species 20: CH2O
  Tci[20] = 1.316 * EPS[20];
  ai[20] = (5.55 * avogadro * avogadro * EPS[20] * boltzmann * 1e-24 * SIG[20] *
            SIG[20] * SIG[20]) /
           (wt[20] * wt[20]);
  bi[20] = 0.855 * avogadro * 1e-24 * SIG[20] * SIG[20] * SIG[20] / (wt[20]);
  acentric_i[20] = 0.0;

  // species 21: HCO
  Tci[21] = 1.316 * EPS[21];
  ai[21] = (5.55 * avogadro * avogadro * EPS[21] * boltzmann * 1e-24 * SIG[21] *
            SIG[21] * SIG[21]) /
           (wt[21] * wt[21]);
  bi[21] = 0.855 * avogadro * 1e-24 * SIG[21] * SIG[21] * SIG[21] / (wt[21]);
  acentric_i[21] = 0.0;

  // species 22: CHO+
  Tci[22] = 1.316 * EPS[22];
  ai[22] = (5.55 * avogadro * avogadro * EPS[22] * boltzmann * 1e-24 * SIG[22] *
            SIG[22] * SIG[22]) /
           (wt[22] * wt[22]);
  bi[22] = 0.855 * avogadro * 1e-24 * SIG[22] * SIG[22] * SIG[22] / (wt[22]);
  acentric_i[22] = 0.0;

  // species 23: H3O+
  Tci[23] = 1.316 * EPS[23];
  ai[23] = (5.55 * avogadro * avogadro * EPS[23] * boltzmann * 1e-24 * SIG[23] *
            SIG[23] * SIG[23]) /
           (wt[23] * wt[23]);
  bi[23] = 0.855 * avogadro * 1e-24 * SIG[23] * SIG[23] * SIG[23] / (wt[23]);
  acentric_i[23] = 0.0;

  // species 24: O2-
  Tci[24] = 1.316 * EPS[24];
  ai[24] = (5.55 * avogadro * avogadro * EPS[24] * boltzmann * 1e-24 * SIG[24] *
            SIG[24] * SIG[24]) /
           (wt[24] * wt[24]);
  bi[24] = 0.855 * avogadro * 1e-24 * SIG[24] * SIG[24] * SIG[24] / (wt[24]);
  acentric_i[24] = 0.0;

  // species 25: E
  Tci[25] = 1.316 * EPS[25];
  ai[25] = (5.55 * avogadro * avogadro * EPS[25] * boltzmann * 1e-24 * SIG[25] *
            SIG[25] * SIG[25]) /
           (wt[25] * wt[25]);
  bi[25] = 0.855 * avogadro * 1e-24 * SIG[25] * SIG[25] * SIG[25] / (wt[25]);
  acentric_i[25] = 0.0;
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
