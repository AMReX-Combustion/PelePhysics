#include "mechanism.H"
const int rmap[NUM_REACTIONS] = {
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18,
  19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
  38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
  57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < NUM_REACTIONS; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of gas species in a gas reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[NUM_GAS_REACTIONS] = {
    3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 5, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4};
  const int kiv[NUM_GAS_REACTIONS * 5] = {
    13, 1,  7,  0, 0, 7,  1,  8,  0,  0, 1,  15, 11, 0,  0, 11, 1,  16, 0,  0,
    9,  0,  11, 0, 0, 2,  3,  0,  0,  0, 1,  2,  4,  0,  0, 9,  2,  10, 0,  0,
    1,  3,  6,  0, 0, 1,  3,  6,  0,  0, 1,  3,  6,  0,  0, 1,  3,  6,  0,  0,
    1,  0,  0,  0, 0, 1,  0,  0,  0,  0, 1,  0,  0,  0,  0, 1,  0,  0,  0,  0,
    1,  4,  5,  0, 0, 15, 9,  1,  0,  0, 15, 9,  1,  0,  0, 0,  2,  1,  4,  0,
    6,  2,  3,  4, 0, 13, 2,  1,  15, 0, 14, 2,  9,  0,  0, 7,  2,  11, 1,  0,
    8,  2,  7,  4, 0, 15, 2,  9,  4,  0, 15, 2,  10, 1,  0, 11, 2,  15, 4,  0,
    16, 2,  11, 4, 0, 9,  3,  10, 2,  0, 11, 3,  15, 6,  0, 1,  3,  2,  4,  0,
    1,  6,  5,  2, 0, 1,  6,  0,  3,  0, 1,  6,  4,  0,  0, 8,  1,  7,  0,  0,
    1,  15, 9,  0, 0, 11, 1,  0,  15, 0, 16, 1,  11, 0,  0, 16, 1,  7,  4,  0,
    16, 1,  14, 5, 0, 0,  4,  1,  5,  0, 4,  5,  2,  0,  0, 6,  4,  5,  3,  0,
    13, 4,  11, 1, 0, 14, 4,  11, 1,  0, 7,  4,  13, 5,  0, 7,  4,  14, 5,  0,
    8,  4,  7,  5, 0, 9,  4,  10, 1,  0, 15, 4,  9,  5,  0, 11, 4,  5,  15, 0,
    16, 4,  11, 5, 0, 13, 6,  11, 4,  0, 7,  6,  8,  3,  0, 9,  6,  10, 4,  0,
    13, 3,  15, 4, 0, 13, 0,  7,  1,  0, 13, 8,  7,  0,  0, 14, 12, 13, 12, 0,
    14, 3,  9,  1, 4, 14, 3,  9,  5,  0, 14, 0,  7,  1,  0, 14, 5,  13, 5,  0,
    14, 8,  7,  0, 0, 14, 9,  13, 9,  0, 14, 10, 13, 10, 0, 14, 10, 11, 9,  0,
    7,  3,  11, 4, 0, 7,  15, 8,  9,  0, 11, 7,  8,  15, 0, 15, 3,  9,  6,  0,
    16, 3,  11, 6, 0};
  const int nuv[NUM_GAS_REACTIONS * 5] = {
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > NUM_GAS_REACTIONS) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 5 + j] + 1;
        nu[j] = nuv[(i - 1) * 5 + j];
      }
    }
  }
}

// Returns the progress rates of each reaction
// Given P, T, and mole fractions
void
CKKFKR(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real q_f[],
  amrex::Real q_r[])
{
  amrex::Real c[NUM_SPECIES]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // convert to SI (mol/cm^3 to mol/m^3)

  // Compute conversion, see Eq 10
  for (int id = 0; id < NUM_GAS_SPECIES; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < NUM_GAS_REACTIONS; ++id) {
    q_f[id] *= 1.0e-6;
    q_r[id] *= 1.0e-6;
  }
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void
progressRateFR(
  amrex::Real* q_f, amrex::Real* q_r, amrex::Real* sc, amrex::Real T)
{
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);
  // compute the Gibbs free energy
  amrex::Real g_RT[NUM_SPECIES];
  gibbs(g_RT, T);
  amrex::Real g_RT_qss[NUM_GAS_SPECIES];
  gibbs_qss(g_RT_qss, T);

  amrex::Real sc_qss[4];
  // Fill sc_qss here
  amrex::Real kf_qss[41], qf_qss[41], qr_qss[41];
  comp_k_f_qss(T, invT, logT, kf_qss);
  comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, T, g_RT, g_RT_qss);
  comp_sc_qss(sc_qss, qf_qss, qr_qss);
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 1.008000;  // H
  awt[1] = 15.999000; // O
  awt[2] = 12.011000; // C
  awt[3] = 14.007000; // N
}

// get atomic weight for all elements
void
CKAWT(amrex::Real* awt)
{
  atomicWeight(awt);
}

// Returns the elemental composition
// of the speciesi (mdim is num of elements)
void
CKNCF(int* ncf)
{
  int kd = 4;
  // Zero ncf
  for (int id = 0; id < kd * NUM_GAS_SPECIES; ++id) {
    ncf[id] = 0;
  }

  // H2
  ncf[0 * kd + 0] = 2; // H

  // H
  ncf[1 * kd + 0] = 1; // H

  // O
  ncf[2 * kd + 1] = 1; // O

  // O2
  ncf[3 * kd + 1] = 2; // O

  // OH
  ncf[4 * kd + 0] = 1; // H
  ncf[4 * kd + 1] = 1; // O

  // H2O
  ncf[5 * kd + 0] = 2; // H
  ncf[5 * kd + 1] = 1; // O

  // HO2
  ncf[6 * kd + 0] = 1; // H
  ncf[6 * kd + 1] = 2; // O

  // CH3
  ncf[7 * kd + 2] = 1; // C
  ncf[7 * kd + 0] = 3; // H

  // CH4
  ncf[8 * kd + 2] = 1; // C
  ncf[8 * kd + 0] = 4; // H

  // CO
  ncf[9 * kd + 2] = 1; // C
  ncf[9 * kd + 1] = 1; // O

  // CO2
  ncf[10 * kd + 2] = 1; // C
  ncf[10 * kd + 1] = 2; // O

  // CH2O
  ncf[11 * kd + 2] = 1; // C
  ncf[11 * kd + 0] = 2; // H
  ncf[11 * kd + 1] = 1; // O

  // N2
  ncf[12 * kd + 3] = 2; // N
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(NUM_ELEMENTS);
  ename[0] = "H";
  ename[1] = "O";
  ename[2] = "C";
  ename[3] = "N";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(NUM_SPECIES);
  kname[0] = "H2";
  kname[1] = "H";
  kname[2] = "O";
  kname[3] = "O2";
  kname[4] = "OH";
  kname[5] = "H2O";
  kname[6] = "HO2";
  kname[7] = "CH3";
  kname[8] = "CH4";
  kname[9] = "CO";
  kname[10] = "CO2";
  kname[11] = "CH2O";
  kname[12] = "N2";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 196> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 13> conc = {0.0};
  for (int n = 0; n < 13; n++) {
    conc[n] = 1.0 / 13.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 14; k++) {
    for (int l = 0; l < 14; l++) {
      if (Jac[14 * k + l] != 0.0) {
        nJdata_tmp = nJdata_tmp + 1;
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the system Jacobian
void
SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 196> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 13> conc = {0.0};
  for (int n = 0; n < 13; n++) {
    conc[n] = 1.0 / 13.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 14; k++) {
    for (int l = 0; l < 14; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[14 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the simplified (for preconditioning) system
// Jacobian
void
SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP)
{
  amrex::GpuArray<amrex::Real, 196> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 13> conc = {0.0};
  for (int n = 0; n < 13; n++) {
    conc[n] = 1.0 / 13.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 14; k++) {
    for (int l = 0; l < 14; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[14 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  nJdata[0] = nJdata_tmp;
}

// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base
// 0
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 196> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 13> conc = {0.0};
  for (int n = 0; n < 13; n++) {
    conc[n] = 1.0 / 13.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 14;
    int offset_col = nc * 14;
    for (int k = 0; k < 14; k++) {
      for (int l = 0; l < 14; l++) {
        if (Jac[14 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l + offset_row;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
      colPtrs[offset_col + (k + 1)] = nJdata_tmp;
    }
  }
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base
// 0
void
SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 196> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 13> conc = {0.0};
  for (int n = 0; n < 13; n++) {
    conc[n] = 1.0 / 13.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 14;
      for (int l = 0; l < 14; l++) {
        for (int k = 0; k < 14; k++) {
          if (Jac[14 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 14;
      for (int l = 0; l < 14; l++) {
        for (int k = 0; k < 14; k++) {
          if (Jac[14 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void
SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtr, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 196> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 13> conc = {0.0};
  for (int n = 0; n < 13; n++) {
    conc[n] = 1.0 / 13.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 14;
      for (int l = 0; l < 14; l++) {
        for (int k = 0; k < 14; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[14 * k + l] != 0.0) {
              colVals[nJdata_tmp - 1] = k + 1 + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 14;
      for (int l = 0; l < 14; l++) {
        for (int k = 0; k < 14; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[14 * k + l] != 0.0) {
              colVals[nJdata_tmp] = k + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// on CPU BASE 0
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP)
{
  amrex::GpuArray<amrex::Real, 196> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 13> conc = {0.0};
  for (int n = 0; n < 13; n++) {
    conc[n] = 1.0 / 13.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 14; k++) {
    for (int l = 0; l < 14; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 14 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[14 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 14 * k + l;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
    colPtrs[k + 1] = nJdata_tmp;
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base)
{
  amrex::GpuArray<amrex::Real, 196> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 13> conc = {0.0};
  for (int n = 0; n < 13; n++) {
    conc[n] = 1.0 / 13.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 14; l++) {
      for (int k = 0; k < 14; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[14 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int l = 0; l < 14; l++) {
      for (int k = 0; k < 14; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[14 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
