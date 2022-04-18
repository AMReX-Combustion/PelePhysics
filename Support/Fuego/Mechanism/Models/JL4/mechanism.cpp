#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"

/*save atomic weights into array */
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 12.011150; /*C */
  awt[1] = 15.999400; /*O */
  awt[2] = 1.007970;  /*H */
  awt[3] = 14.006700; /*N */
}

/*get atomic weight for all elements */
void
CKAWT(amrex::Real* awt)
{
  atomicWeight(awt);
}

/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void
CKNCF(int* ncf)
{
  int id; /*loop counter */
  int kd = 4;
  /*Zero ncf */
  for (id = 0; id < kd * 7; ++id) {
    ncf[id] = 0;
  }

  /*CH4 */
  ncf[0 * kd + 0] = 1; /*C */
  ncf[0 * kd + 2] = 4; /*H */

  /*O2 */
  ncf[1 * kd + 1] = 2; /*O */

  /*H2O */
  ncf[2 * kd + 2] = 2; /*H */
  ncf[2 * kd + 1] = 1; /*O */

  /*N2 */
  ncf[3 * kd + 3] = 2; /*N */

  /*CO */
  ncf[4 * kd + 0] = 1; /*C */
  ncf[4 * kd + 1] = 1; /*O */

  /*CO2 */
  ncf[5 * kd + 0] = 1; /*C */
  ncf[5 * kd + 1] = 2; /*O */

  /*H2 */
  ncf[6 * kd + 2] = 2; /*H */
}

/* Returns the vector of strings of element names */
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(4);
  ename[0] = "C";
  ename[1] = "O";
  ename[2] = "H";
  ename[3] = "N";
}

/* Returns the vector of strings of species names */
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(7);
  kname[0] = "CH4";
  kname[1] = "O2";
  kname[2] = "H2O";
  kname[3] = "N2";
  kname[4] = "CO";
  kname[5] = "CO2";
  kname[6] = "H2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 64> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 7> conc = {0.0};
  for (int n = 0; n < 7; n++) {
    conc[n] = 1.0 / 7.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 8; k++) {
    for (int l = 0; l < 8; l++) {
      if (Jac[8 * k + l] != 0.0) {
        nJdata_tmp = nJdata_tmp + 1;
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

/*compute the sparsity pattern of the system Jacobian */
void
SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 64> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 7> conc = {0.0};
  for (int n = 0; n < 7; n++) {
    conc[n] = 1.0 / 7.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 8; k++) {
    for (int l = 0; l < 8; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[8 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

/*compute the sparsity pattern of the simplified (for preconditioning) system
 * Jacobian */
void
SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP)
{
  amrex::GpuArray<amrex::Real, 64> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 7> conc = {0.0};
  for (int n = 0; n < 7; n++) {
    conc[n] = 1.0 / 7.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 8; k++) {
    for (int l = 0; l < 8; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[8 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  nJdata[0] = nJdata_tmp;
}

/*compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0
 */
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 64> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 7> conc = {0.0};
  for (int n = 0; n < 7; n++) {
    conc[n] = 1.0 / 7.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 8;
    int offset_col = nc * 8;
    for (int k = 0; k < 8; k++) {
      for (int l = 0; l < 8; l++) {
        if (Jac[8 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l + offset_row;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
      colPtrs[offset_col + (k + 1)] = nJdata_tmp;
    }
  }
}

/*compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0
 */
void
SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 64> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 7> conc = {0.0};
  for (int n = 0; n < 7; n++) {
    conc[n] = 1.0 / 7.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 8;
      for (int l = 0; l < 8; l++) {
        for (int k = 0; k < 8; k++) {
          if (Jac[8 * k + l] != 0.0) {
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
      int offset = nc * 8;
      for (int l = 0; l < 8; l++) {
        for (int k = 0; k < 8; k++) {
          if (Jac[8 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

/*compute the sparsity pattern of the system Jacobian */
/*CSR format BASE is user choice */
void
SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtr, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 64> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 7> conc = {0.0};
  for (int n = 0; n < 7; n++) {
    conc[n] = 1.0 / 7.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 8;
      for (int l = 0; l < 8; l++) {
        for (int k = 0; k < 8; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[8 * k + l] != 0.0) {
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
      int offset = nc * 8;
      for (int l = 0; l < 8; l++) {
        for (int k = 0; k < 8; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[8 * k + l] != 0.0) {
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

/*compute the sparsity pattern of the simplified (for precond) system Jacobian
 * on CPU */
/*BASE 0 */
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP)
{
  amrex::GpuArray<amrex::Real, 64> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 7> conc = {0.0};
  for (int n = 0; n < 7; n++) {
    conc[n] = 1.0 / 7.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 8; k++) {
    for (int l = 0; l < 8; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 8 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[8 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 8 * k + l;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
    colPtrs[k + 1] = nJdata_tmp;
  }
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian
 */
/*CSR format BASE is under choice */
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base)
{
  amrex::GpuArray<amrex::Real, 64> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 7> conc = {0.0};
  for (int n = 0; n < 7; n++) {
    conc[n] = 1.0 / 7.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 8; l++) {
      for (int k = 0; k < 8; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[8 * k + l] != 0.0) {
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
    for (int l = 0; l < 8; l++) {
      for (int k = 0; k < 8; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[8 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}

#endif
