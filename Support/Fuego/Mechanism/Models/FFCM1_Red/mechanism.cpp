#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"

/*save atomic weights into array */
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 15.999400; /*O */
  awt[1] = 1.007970;  /*H */
  awt[2] = 12.011150; /*C */
  awt[3] = 14.006700; /*N */
  awt[4] = 4.002600;  /*HE */
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
  int kd = 5;
  /*Zero ncf */
  for (id = 0; id < kd * 21; ++id) {
    ncf[id] = 0;
  }

  /*N2 */
  ncf[0 * kd + 3] = 2; /*N */

  /*H2 */
  ncf[1 * kd + 1] = 2; /*H */

  /*H */
  ncf[2 * kd + 1] = 1; /*H */

  /*O */
  ncf[3 * kd + 0] = 1; /*O */

  /*O2 */
  ncf[4 * kd + 0] = 2; /*O */

  /*OH */
  ncf[5 * kd + 0] = 1; /*O */
  ncf[5 * kd + 1] = 1; /*H */

  /*H2O */
  ncf[6 * kd + 1] = 2; /*H */
  ncf[6 * kd + 0] = 1; /*O */

  /*HO2 */
  ncf[7 * kd + 1] = 1; /*H */
  ncf[7 * kd + 0] = 2; /*O */

  /*CO */
  ncf[8 * kd + 2] = 1; /*C */
  ncf[8 * kd + 0] = 1; /*O */

  /*CO2 */
  ncf[9 * kd + 2] = 1; /*C */
  ncf[9 * kd + 0] = 2; /*O */

  /*CH */
  ncf[10 * kd + 2] = 1; /*C */
  ncf[10 * kd + 1] = 1; /*H */

  /*CH2 */
  ncf[11 * kd + 2] = 1; /*C */
  ncf[11 * kd + 1] = 2; /*H */

  /*CH2(S) */
  ncf[12 * kd + 2] = 1; /*C */
  ncf[12 * kd + 1] = 2; /*H */

  /*CH3 */
  ncf[13 * kd + 2] = 1; /*C */
  ncf[13 * kd + 1] = 3; /*H */

  /*CH4 */
  ncf[14 * kd + 2] = 1; /*C */
  ncf[14 * kd + 1] = 4; /*H */

  /*HCO */
  ncf[15 * kd + 2] = 1; /*C */
  ncf[15 * kd + 1] = 1; /*H */
  ncf[15 * kd + 0] = 1; /*O */

  /*CH2O */
  ncf[16 * kd + 1] = 2; /*H */
  ncf[16 * kd + 2] = 1; /*C */
  ncf[16 * kd + 0] = 1; /*O */

  /*CH2OH */
  ncf[17 * kd + 2] = 1; /*C */
  ncf[17 * kd + 1] = 3; /*H */
  ncf[17 * kd + 0] = 1; /*O */

  /*CH3O */
  ncf[18 * kd + 2] = 1; /*C */
  ncf[18 * kd + 1] = 3; /*H */
  ncf[18 * kd + 0] = 1; /*O */

  /*CH3OH */
  ncf[19 * kd + 2] = 1; /*C */
  ncf[19 * kd + 1] = 4; /*H */
  ncf[19 * kd + 0] = 1; /*O */

  /*HE */
  ncf[20 * kd + 4] = 1; /*HE */
}

/* Returns the vector of strings of element names */
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(5);
  ename[0] = "O";
  ename[1] = "H";
  ename[2] = "C";
  ename[3] = "N";
  ename[4] = "HE";
}

/* Returns the vector of strings of species names */
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(21);
  kname[0] = "N2";
  kname[1] = "H2";
  kname[2] = "H";
  kname[3] = "O";
  kname[4] = "O2";
  kname[5] = "OH";
  kname[6] = "H2O";
  kname[7] = "HO2";
  kname[8] = "CO";
  kname[9] = "CO2";
  kname[10] = "CH";
  kname[11] = "CH2";
  kname[12] = "CH2(S)";
  kname[13] = "CH3";
  kname[14] = "CH4";
  kname[15] = "HCO";
  kname[16] = "CH2O";
  kname[17] = "CH2OH";
  kname[18] = "CH3O";
  kname[19] = "CH3OH";
  kname[20] = "HE";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 22; k++) {
    for (int l = 0; l < 22; l++) {
      if (Jac[22 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 22; k++) {
    for (int l = 0; l < 22; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[22 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 22; k++) {
    for (int l = 0; l < 22; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[22 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 22;
    int offset_col = nc * 22;
    for (int k = 0; k < 22; k++) {
      for (int l = 0; l < 22; l++) {
        if (Jac[22 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 22;
      for (int l = 0; l < 22; l++) {
        for (int k = 0; k < 22; k++) {
          if (Jac[22 * k + l] != 0.0) {
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
      int offset = nc * 22;
      for (int l = 0; l < 22; l++) {
        for (int k = 0; k < 22; k++) {
          if (Jac[22 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 22;
      for (int l = 0; l < 22; l++) {
        for (int k = 0; k < 22; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[22 * k + l] != 0.0) {
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
      int offset = nc * 22;
      for (int l = 0; l < 22; l++) {
        for (int k = 0; k < 22; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[22 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 22; k++) {
    for (int l = 0; l < 22; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 22 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[22 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 22 * k + l;
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
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 22; l++) {
      for (int k = 0; k < 22; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[22 * k + l] != 0.0) {
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
    for (int l = 0; l < 22; l++) {
      for (int k = 0; k < 22; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[22 * k + l] != 0.0) {
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
