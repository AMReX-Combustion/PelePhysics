#include "mechanism.H"
const int rmap[73] = {
  28, 29, 31, 33, 38, 0,  1,  8,  15, 16, 17, 18, 20, 21, 22, 23, 24, 69, 70,
  2,  3,  4,  5,  6,  7,  9,  10, 11, 12, 13, 14, 19, 25, 26, 27, 30, 32, 34,
  35, 36, 37, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
  55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 71, 72};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 73; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[73] = {2, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3,
                      4, 2, 2, 2, 2, 3, 4, 4, 3, 3, 3, 4, 3, 4, 3, 4, 4, 4, 4,
                      3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3,
                      4, 5, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4};
  const int kiv[365] = {
    2,  3,  0,  0,  0, 1,  2,  4,  0,  0, 0,  2,  1,  4,  0, 6,  2,  3,  4,  0,
    7,  2,  1,  13, 0, 8,  2,  11, 0,  0, 9,  2,  14, 1,  0, 10, 2,  9,  4,  0,
    11, 2,  12, 0,  0, 13, 2,  11, 4,  0, 13, 2,  12, 1,  0, 14, 2,  13, 4,  0,
    15, 2,  14, 4,  0, 11, 3,  12, 2,  0, 14, 3,  13, 6,  0, 1,  3,  6,  0,  0,
    1,  3,  6,  0,  0, 1,  3,  6,  0,  0, 1,  3,  6,  0,  0, 1,  3,  2,  4,  0,
    1,  0,  0,  0,  0, 1,  0,  0,  0,  0, 1,  0,  0,  0,  0, 1,  0,  0,  0,  0,
    1,  4,  5,  0,  0, 1,  6,  5,  2,  0, 1,  6,  0,  3,  0, 1,  6,  4,  0,  0,
    7,  1,  9,  0,  0, 9,  1,  10, 0,  0, 10, 1,  9,  0,  0, 1,  13, 14, 0,  0,
    1,  13, 11, 0,  0, 14, 1,  15, 0,  0, 14, 1,  0,  13, 0, 15, 1,  14, 0,  0,
    15, 1,  9,  4,  0, 15, 1,  8,  5,  0, 11, 0,  14, 0,  0, 0,  4,  1,  5,  0,
    4,  5,  2,  0,  0, 6,  4,  5,  3,  0, 7,  4,  14, 1,  0, 8,  4,  14, 1,  0,
    9,  4,  7,  5,  0, 9,  4,  8,  5,  0, 10, 4,  9,  5,  0, 11, 4,  12, 1,  0,
    13, 4,  11, 5,  0, 14, 4,  5,  13, 0, 15, 4,  14, 5,  0, 7,  6,  14, 4,  0,
    9,  6,  10, 3,  0, 11, 6,  12, 4,  0, 7,  3,  13, 4,  0, 7,  0,  9,  1,  0,
    7,  10, 9,  0,  0, 8,  16, 7,  16, 0, 8,  3,  11, 1,  4, 8,  3,  11, 5,  0,
    8,  0,  9,  1,  0, 8,  5,  7,  5,  0, 8,  10, 9,  0,  0, 8,  11, 7,  11, 0,
    8,  12, 7,  12, 0, 8,  12, 14, 11, 0, 9,  3,  14, 4,  0, 9,  13, 10, 11, 0,
    14, 9,  10, 13, 0, 13, 11, 1,  0,  0, 13, 11, 1,  0,  0, 13, 3,  11, 6,  0,
    15, 3,  14, 6,  0};
  const int nuv[365] = {
    -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > 73) {
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

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void
CKKFKR(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real q_f[],
  amrex::Real q_r[])
{
  amrex::Real c[17]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 17; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 73; ++id) {
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
  const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
  amrex::Real invT = 1.0 / tc[1];
  // compute the Gibbs free energy
  amrex::Real g_RT[17];
  gibbs(g_RT, tc);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 15.999000; // O
  awt[1] = 1.008000;  // H
  awt[2] = 12.011000; // C
  awt[3] = 14.007000; // N
  awt[4] = 39.950000; // Ar
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
  int kd = 5;
  // Zero ncf
  for (int id = 0; id < kd * 17; ++id) {
    ncf[id] = 0;
  }

  // H2
  ncf[0 * kd + 1] = 2; // H

  // H
  ncf[1 * kd + 1] = 1; // H

  // O
  ncf[2 * kd + 0] = 1; // O

  // O2
  ncf[3 * kd + 0] = 2; // O

  // OH
  ncf[4 * kd + 1] = 1; // H
  ncf[4 * kd + 0] = 1; // O

  // H2O
  ncf[5 * kd + 1] = 2; // H
  ncf[5 * kd + 0] = 1; // O

  // HO2
  ncf[6 * kd + 1] = 1; // H
  ncf[6 * kd + 0] = 2; // O

  // CH2
  ncf[7 * kd + 2] = 1; // C
  ncf[7 * kd + 1] = 2; // H

  // CH2(S)
  ncf[8 * kd + 2] = 1; // C
  ncf[8 * kd + 1] = 2; // H

  // CH3
  ncf[9 * kd + 2] = 1; // C
  ncf[9 * kd + 1] = 3; // H

  // CH4
  ncf[10 * kd + 2] = 1; // C
  ncf[10 * kd + 1] = 4; // H

  // CO
  ncf[11 * kd + 2] = 1; // C
  ncf[11 * kd + 0] = 1; // O

  // CO2
  ncf[12 * kd + 2] = 1; // C
  ncf[12 * kd + 0] = 2; // O

  // HCO
  ncf[13 * kd + 2] = 1; // C
  ncf[13 * kd + 1] = 1; // H
  ncf[13 * kd + 0] = 1; // O

  // CH2O
  ncf[14 * kd + 2] = 1; // C
  ncf[14 * kd + 1] = 2; // H
  ncf[14 * kd + 0] = 1; // O

  // CH2OH
  ncf[15 * kd + 2] = 1; // C
  ncf[15 * kd + 1] = 3; // H
  ncf[15 * kd + 0] = 1; // O

  // N2
  ncf[16 * kd + 3] = 2; // N
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(5);
  ename[0] = "O";
  ename[1] = "H";
  ename[2] = "C";
  ename[3] = "N";
  ename[4] = "Ar";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(17);
  kname[0] = "H2";
  kname[1] = "H";
  kname[2] = "O";
  kname[3] = "O2";
  kname[4] = "OH";
  kname[5] = "H2O";
  kname[6] = "HO2";
  kname[7] = "CH2";
  kname[8] = "CH2(S)";
  kname[9] = "CH3";
  kname[10] = "CH4";
  kname[11] = "CO";
  kname[12] = "CO2";
  kname[13] = "HCO";
  kname[14] = "CH2O";
  kname[15] = "CH2OH";
  kname[16] = "N2";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 324> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 17> conc = {0.0};
  for (int n = 0; n < 17; n++) {
    conc[n] = 1.0 / 17.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 18; k++) {
    for (int l = 0; l < 18; l++) {
      if (Jac[18 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 324> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 17> conc = {0.0};
  for (int n = 0; n < 17; n++) {
    conc[n] = 1.0 / 17.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 18; k++) {
    for (int l = 0; l < 18; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[18 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 324> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 17> conc = {0.0};
  for (int n = 0; n < 17; n++) {
    conc[n] = 1.0 / 17.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 18; k++) {
    for (int l = 0; l < 18; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[18 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 324> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 17> conc = {0.0};
  for (int n = 0; n < 17; n++) {
    conc[n] = 1.0 / 17.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 18;
    int offset_col = nc * 18;
    for (int k = 0; k < 18; k++) {
      for (int l = 0; l < 18; l++) {
        if (Jac[18 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 324> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 17> conc = {0.0};
  for (int n = 0; n < 17; n++) {
    conc[n] = 1.0 / 17.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 18;
      for (int l = 0; l < 18; l++) {
        for (int k = 0; k < 18; k++) {
          if (Jac[18 * k + l] != 0.0) {
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
      int offset = nc * 18;
      for (int l = 0; l < 18; l++) {
        for (int k = 0; k < 18; k++) {
          if (Jac[18 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 324> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 17> conc = {0.0};
  for (int n = 0; n < 17; n++) {
    conc[n] = 1.0 / 17.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 18;
      for (int l = 0; l < 18; l++) {
        for (int k = 0; k < 18; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[18 * k + l] != 0.0) {
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
      int offset = nc * 18;
      for (int l = 0; l < 18; l++) {
        for (int k = 0; k < 18; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[18 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 324> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 17> conc = {0.0};
  for (int n = 0; n < 17; n++) {
    conc[n] = 1.0 / 17.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 18; k++) {
    for (int l = 0; l < 18; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 18 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[18 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 18 * k + l;
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
  amrex::GpuArray<amrex::Real, 324> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 17> conc = {0.0};
  for (int n = 0; n < 17; n++) {
    conc[n] = 1.0 / 17.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 18; l++) {
      for (int k = 0; k < 18; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[18 * k + l] != 0.0) {
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
    for (int l = 0; l < 18; l++) {
      for (int k = 0; k < 18; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[18 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
