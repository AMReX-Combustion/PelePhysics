#include "mechanism.H"

// Returns 0-based map of reaction order
void
GET_RMAP(int* /*_rmap*/)
{
}

// Returns a count of gas species in a gas reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int* /*ki*/, int* /*nu*/)
{
  if (i < 1) {
    // Return max num species per reaction
    nspec = 0;
  } else {
    nspec = -1;
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
  amrex::Real c[35]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 35; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void
progressRateFR(
  amrex::Real* /*q_f*/,
  amrex::Real* /*q_r*/,
  amrex::Real* /*sc*/,
  amrex::Real /*T*/)
{
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(4);
  ename[0] = "O";
  ename[1] = "H";
  ename[2] = "C";
  ename[3] = "N";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(35);
  kname[0] = "NC12H26";
  kname[1] = "H";
  kname[2] = "O";
  kname[3] = "OH";
  kname[4] = "HO2";
  kname[5] = "H2";
  kname[6] = "H2O";
  kname[7] = "H2O2";
  kname[8] = "O2";
  kname[9] = "CH3";
  kname[10] = "CH4";
  kname[11] = "CH2O";
  kname[12] = "CO";
  kname[13] = "CO2";
  kname[14] = "C2H2";
  kname[15] = "C2H4";
  kname[16] = "C2H6";
  kname[17] = "CH2CHO";
  kname[18] = "aC3H5";
  kname[19] = "C3H6";
  kname[20] = "C2H3CHO";
  kname[21] = "C4H7";
  kname[22] = "C4H81";
  kname[23] = "C5H9";
  kname[24] = "C5H10";
  kname[25] = "C6H12";
  kname[26] = "C7H14";
  kname[27] = "C8H16";
  kname[28] = "C9H18";
  kname[29] = "PXC9H19";
  kname[30] = "C10H20";
  kname[31] = "C12H24";
  kname[32] = "C12H25O2";
  kname[33] = "OC12H23OOH";
  kname[34] = "N2";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 36; k++) {
    for (int l = 0; l < 36; l++) {
      if (Jac[36 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 36; k++) {
    for (int l = 0; l < 36; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[36 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 36; k++) {
    for (int l = 0; l < 36; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[36 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 36;
    int offset_col = nc * 36;
    for (int k = 0; k < 36; k++) {
      for (int l = 0; l < 36; l++) {
        if (Jac[36 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 36;
      for (int l = 0; l < 36; l++) {
        for (int k = 0; k < 36; k++) {
          if (Jac[36 * k + l] != 0.0) {
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
      int offset = nc * 36;
      for (int l = 0; l < 36; l++) {
        for (int k = 0; k < 36; k++) {
          if (Jac[36 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 36;
      for (int l = 0; l < 36; l++) {
        for (int k = 0; k < 36; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[36 * k + l] != 0.0) {
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
      int offset = nc * 36;
      for (int l = 0; l < 36; l++) {
        for (int k = 0; k < 36; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[36 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 36; k++) {
    for (int l = 0; l < 36; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 36 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[36 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 36 * k + l;
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 36; l++) {
      for (int k = 0; k < 36; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[36 * k + l] != 0.0) {
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
    for (int l = 0; l < 36; l++) {
      for (int k = 0; k < 36; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[36 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
