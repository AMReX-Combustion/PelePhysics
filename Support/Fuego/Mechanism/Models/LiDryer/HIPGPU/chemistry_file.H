#ifndef CHEMISTRY_FILE_H
#define CHEMISTRY_FILE_H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <AMReX_Gpu.H>

/* PURE CPU stuff */
#if !defined(AMREX_USE_CUDA) && !defined(AMREX_USE_HIP)
namespace thermo
{

    extern double fwd_A[21], fwd_beta[21], fwd_Ea[21];
    extern double low_A[21], low_beta[21], low_Ea[21];
    extern double rev_A[21], rev_beta[21], rev_Ea[21];
    extern double troe_a[21],troe_Ts[21], troe_Tss[21], troe_Tsss[21];
    extern double sri_a[21], sri_b[21], sri_c[21], sri_d[21], sri_e[21];
    extern double activation_units[21], prefactor_units[21], phase_units[21];
    extern int is_PD[21], troe_len[21], sri_len[21], nTB[21], *TBid[21];
    extern double *TB[21];
    extern std::vector<std::vector<double>> kiv; 
    extern std::vector<std::vector<double>> nuv; 
}
/* Vectorized stuff */
void VCKYTX(int *  np, double *  y, double *  x);
void VCKHMS(int *  np, double *  T, double *  ums);
void VCKWYR(int *  np, double *  rho, double *  T,
            double *  y,
            double *  wdot);
void VCKPY(int *  np, double *  rho, double *  T, double *  y, double *  P);
void vproductionRate(int npt, double *  wdot, double *  c, double *  T);
void vcomp_k_f(int npt, double *  k_f_s, double *  tc, double *  invT);
void vcomp_gibbs(int npt, double *  g_RT, double *  tc);
void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT);
void vcomp_wdot(int npt, double *  wdot, double *  mixture, double *  sc,
                double *  k_f_s, double *  Kc_s,
                double *  tc, double *  invT, double *  T);
#endif

/* INIT and FINALIZE stuff */
void CKINIT();
void CKFINALIZE();
void atomicWeight(double *  awt);
/* MISC */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double *  rval, int * kerr, int lenline);
#if !defined(AMREX_USE_CUDA) && !defined(AMREX_USE_HIP)
void CKINU(int * i, int * nspec, int * ki, int * nu);
#endif
void CKNCF(int * ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYME(int * kname, int * lenkname);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
//void GET_CRITPARAMS(double *  Tci, double *  ai, double *  bi, double *  acentric_i);
void CKAWT(double *  awt);
/*Transport function declarations */
void egtransetLENIMC(int* LENIMC);
void egtransetLENRMC(int* LENRMC);
void egtransetNO(int* NO);
void egtransetKK(int* KK);
void egtransetNLITE(int* NLITE);
void egtransetPATM(double* PATM);
void egtransetWT(double* WT);
void egtransetEPS(double* EPS);
void egtransetSIG(double* SIG);
void egtransetDIP(double* DIP);
void egtransetPOL(double* POL);
void egtransetZROT(double* ZROT);
void egtransetNLIN(int* NLIN);
void egtransetCOFETA(double* COFETA);
void egtransetCOFLAM(double* COFLAM);
void egtransetCOFD(double* COFD);
void egtransetKTDIF(int* KTDIF);
/* SPARSE INFORMATION */
void SPARSITY_INFO(int * nJdata, int * consP, int NCELLS);
void SPARSITY_INFO_SYST(int * nJdata, int * consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int * nJdata, int * consP);
void SPARSITY_PREPROC_CSC(int * rowVals, int * colPtrs, int * consP, int NCELLS);
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base);

/* PROD RATE STUFF */
#if !defined(AMREX_USE_CUDA) && !defined(AMREX_USE_HIP)
void productionRate(double *  wdot, double *  sc, double T);
void comp_qfqr(double *  q_f, double *  q_r, double *  sc, double *  tc, double invT);
void comp_k_f(double *  tc, double invT, double *  k_f);
void comp_Kc(double *  tc, double invT, double *  Kc);
void progressRate(double *  qdot, double *  speciesConc, double T);
void progressRateFR(double *  q_f, double *  q_r, double *  speciesConc, double T);
void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r);
void CKQC(double *  T, double *  C, double *  qdot);
void CKQYP(double *  P, double *  T, double *  y, double *  qdot);
void CKQXP(double *  P, double *  T, double *  x, double *  qdot);
void CKQYR(double *  rho, double *  T, double *  y, double *  qdot);
void CKQXR(double *  rho, double *  T, double *  x, double *  qdot);
void aJacobian(double *  J, double *  sc, double T, int consP);
/* Doubly defined things */
#endif

#endif
