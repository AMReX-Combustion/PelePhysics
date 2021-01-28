#ifndef CHEMISTRY_FILE_H
#define CHEMISTRY_FILE_H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* PURE CPU stuff -- no use on GPU  */
#ifndef AMREX_USE_GPU
namespace thermo
{

    extern amrex::Real fwd_A[27], fwd_beta[27], fwd_Ea[27];
    extern amrex::Real low_A[27], low_beta[27], low_Ea[27];
    extern amrex::Real rev_A[27], rev_beta[27], rev_Ea[27];
    extern amrex::Real troe_a[27],troe_Ts[27], troe_Tss[27], troe_Tsss[27];
    extern amrex::Real sri_a[27], sri_b[27], sri_c[27], sri_d[27], sri_e[27];
    extern amrex::Real activation_units[27], prefactor_units[27], phase_units[27];
    extern int is_PD[27], troe_len[27], sri_len[27], nTB[27], *TBid[27];
    extern amrex::Real *TB[27];
}
/* Vectorized stuff  */
void VCKYTX(int *  np, amrex::Real *  y, amrex::Real *  x);
void VCKHMS(int *  np, amrex::Real *  T, amrex::Real *  ums);
void VCKWYR(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  wdot);
void VCKPY(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  P);
void vproductionRate(int npt, amrex::Real *  wdot, amrex::Real *  c, amrex::Real *  T);
void vcomp_k_f(int npt, amrex::Real *  k_f_s, amrex::Real *  tc, amrex::Real *  invT);
void vcomp_gibbs(int npt, amrex::Real *  g_RT, amrex::Real *  tc);
void vcomp_Kc(int npt, amrex::Real *  Kc_s, amrex::Real *  g_RT, amrex::Real *  invT);
void vcomp_wdot(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
                amrex::Real *  k_f_s, amrex::Real *  Kc_s,
                amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T);
/* MISC  */
/* PROD RATE STUFF  */
void productionRate(amrex::Real *  wdot, amrex::Real *  sc, amrex::Real T);
void comp_qfqr(amrex::Real *  q_f, amrex::Real *  q_r, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT);
void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f);
void comp_Kc(amrex::Real *  tc, amrex::Real invT, amrex::Real *  Kc);
void progressRate(amrex::Real *  qdot, amrex::Real *  speciesConc, amrex::Real T);
void progressRateFR(amrex::Real *  q_f, amrex::Real *  q_r, amrex::Real *  speciesConc, amrex::Real T);
void CKKFKR(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  q_f, amrex::Real *  q_r);
void CKQC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  qdot);
void CKQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot);
void CKQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot);
void CKQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot);
void CKQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot);
void aJacobian(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP);
#endif


/* ALWAYS on CPU stuff -- can have different def depending on if we are CPU or GPU based. Defined in mechanism.cpp  */
/* Transport function declarations  */
void egtransetLENIMC(int* LENIMC);
void egtransetLENRMC(int* LENRMC);
void egtransetNO(int* NO);
void egtransetKK(int* KK);
void egtransetNLITE(int* NLITE);
void egtransetPATM(amrex::Real* PATM);
void egtransetWT(amrex::Real* WT);
void egtransetEPS(amrex::Real* EPS);
void egtransetSIG(amrex::Real* SIG);
void egtransetDIP(amrex::Real* DIP);
void egtransetPOL(amrex::Real* POL);
void egtransetZROT(amrex::Real* ZROT);
void egtransetNLIN(int* NLIN);
void egtransetCOFETA(amrex::Real* COFETA);
void egtransetCOFLAM(amrex::Real* COFLAM);
void egtransetCOFD(amrex::Real* COFD);
void egtransetKTDIF(int* KTDIF);
void egtransetCOFTD(amrex::Real* COFTD);
/* INIT and FINALIZE stuff  */
void CKINIT();
void CKFINALIZE();
void atomicWeight(amrex::Real *  awt);
/* MISC  */
void CKNCF(int * ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
void CKAWT(amrex::Real *  awt);
/* SPARSE INFORMATION  */
void SPARSITY_INFO(int * nJdata, int * consP, int NCELLS);
void SPARSITY_INFO_SYST(int * nJdata, int * consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int * nJdata, int * consP);
void SPARSITY_PREPROC_CSC(int * rowVals, int * colPtrs, int * consP, int NCELLS);
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base);

/* GPU stuff that default to CPU if not compiled with GPU  */


/*A few mechanism parameters */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKINDX(int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 6;
    *kk = 13;
    *ii = 27;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/* inverse molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
    imw_new[0] = 1.0/1.007970;/*H */
    imw_new[1] = 1.0/2.015940;/*H2 */
    imw_new[2] = 1.0/15.999400;/*O */
    imw_new[3] = 1.0/17.007370;/*OH */
    imw_new[4] = 1.0/18.015340;/*H2O */
    imw_new[5] = 1.0/31.998800;/*O2 */
    imw_new[6] = 1.0/33.006770;/*HO2 */
    imw_new[7] = 1.0/34.014740;/*H2O2 */
    imw_new[8] = 1.0/28.013400;/*N2 */
    imw_new[9] = 1.0/39.948000;/*AR */
    imw_new[10] = 1.0/4.002600;/*HE */
    imw_new[11] = 1.0/28.010550;/*CO */
    imw_new[12] = 1.0/44.009950;/*CO2 */
}

/* molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_mw(amrex::Real *mw_new){
    mw_new[0] = 1.007970;/*H */
    mw_new[1] = 2.015940;/*H2 */
    mw_new[2] = 15.999400;/*O */
    mw_new[3] = 17.007370;/*OH */
    mw_new[4] = 18.015340;/*H2O */
    mw_new[5] = 31.998800;/*O2 */
    mw_new[6] = 33.006770;/*HO2 */
    mw_new[7] = 34.014740;/*H2O2 */
    mw_new[8] = 28.013400;/*N2 */
    mw_new[9] = 39.948000;/*AR */
    mw_new[10] = 4.002600;/*HE */
    mw_new[11] = 28.010550;/*CO */
    mw_new[12] = 44.009950;/*CO2 */
}



/* Returns R, Rc, Patm */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRP(amrex::Real *  ru, amrex::Real *  ruc, amrex::Real *  pa)
{
     *ru  = 8.31446261815324e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cv_R(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +2.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 2: O */
        species[2] =
            +1.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 4: H2O */
        species[4] =
            +2.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 5: O2 */
        species[5] =
            +2.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +2.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
        /*species 9: AR */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            +2.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +1.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +1.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 2: O */
        species[2] =
            +1.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 3: OH */
        species[3] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 4: H2O */
        species[4] =
            +1.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 5: O2 */
        species[5] =
            +2.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 9: AR */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            +2.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +3.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cp_R(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 2: O */
        species[2] =
            +2.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 3: OH */
        species[3] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 4: H2O */
        species[4] =
            +3.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 5: O2 */
        species[5] =
            +3.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
        /*species 9: AR */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            +3.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +2.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 2: O */
        species[2] =
            +2.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 3: OH */
        species[3] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 4: H2O */
        species[4] =
            +2.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 5: O2 */
        species[5] =
            +3.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 9: AR */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            +3.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +4.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
    }
    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gibbs(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.547163000000000e+04 * invT
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -1.012521000000000e+03 * invT
            +6.592218000000000e+00
            -3.298124000000000e+00 * tc[0]
            -4.124721000000000e-04 * tc[1]
            +1.357169166666667e-07 * tc[2]
            +7.896194999999999e-12 * tc[3]
            -2.067436000000000e-14 * tc[4];
        /*species 2: O */
        species[2] =
            +2.914764000000000e+04 * invT
            -1.756599999999997e-02
            -2.946429000000000e+00 * tc[0]
            +8.190830000000000e-04 * tc[1]
            -4.035053333333333e-07 * tc[2]
            +1.335702500000000e-10 * tc[3]
            -1.945348000000000e-14 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.346309130000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.020811000000000e+04 * invT
            +7.966090000000001e-01
            -3.386842000000000e+00 * tc[0]
            -1.737491000000000e-03 * tc[1]
            +1.059116000000000e-06 * tc[2]
            -5.807150833333333e-10 * tc[3]
            +1.253294000000000e-13 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.005249000000000e+03 * invT
            -2.821802000000000e+00
            -3.212936000000000e+00 * tc[0]
            -5.637430000000000e-04 * tc[1]
            +9.593583333333333e-08 * tc[2]
            -1.094897500000000e-10 * tc[3]
            +4.384277000000000e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.766315000000000e+04 * invT
            -3.396609000000000e+00
            -3.388754000000000e+00 * tc[0]
            -3.284613000000000e-03 * tc[1]
            +2.475021666666666e-08 * tc[2]
            +3.854838333333333e-10 * tc[3]
            -1.235757500000000e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.020900000000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
        /*species 9: AR */
        species[9] =
            -7.453750000000000e+02 * invT
            -1.866001000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            -7.453750000000000e+02 * invT
            +1.584651200000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.431054000000000e+04 * invT
            -1.586445000000000e+00
            -3.262452000000000e+00 * tc[0]
            -7.559705000000000e-04 * tc[1]
            +6.469591666666667e-07 * tc[2]
            -4.651620000000000e-10 * tc[3]
            +1.237475500000000e-13 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.837314000000000e+04 * invT
            -7.912765000000000e+00
            -2.275725000000000e+00 * tc[0]
            -4.961036000000000e-03 * tc[1]
            +1.734851666666667e-06 * tc[2]
            -5.722239166666667e-10 * tc[3]
            +1.058640000000000e-13 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.547163000000000e+04 * invT
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -8.350340000000000e+02 * invT
            +4.346533000000000e+00
            -2.991423000000000e+00 * tc[0]
            -3.500322000000000e-04 * tc[1]
            +9.389715000000000e-09 * tc[2]
            +7.692981666666667e-13 * tc[3]
            -7.913760000000000e-17 * tc[4];
        /*species 2: O */
        species[2] =
            +2.923080000000000e+04 * invT
            -2.378248000000000e+00
            -2.542060000000000e+00 * tc[0]
            +1.377531000000000e-05 * tc[1]
            +5.171338333333333e-10 * tc[2]
            -3.792555833333334e-13 * tc[3]
            +2.184026000000000e-17 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.683628750000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
        /*species 4: H2O */
        species[4] =
            -2.989921000000000e+04 * invT
            -4.190671000000000e+00
            -2.672146000000000e+00 * tc[0]
            -1.528146500000000e-03 * tc[1]
            +1.455043333333333e-07 * tc[2]
            -1.000830000000000e-11 * tc[3]
            +3.195809000000000e-16 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.233930000000000e+03 * invT
            +5.084119999999999e-01
            -3.697578000000000e+00 * tc[0]
            -3.067598500000000e-04 * tc[1]
            +2.098070000000000e-08 * tc[2]
            -1.479400833333333e-12 * tc[3]
            +5.682175000000001e-17 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.800696000000000e+04 * invT
            +4.072030000000000e+00
            -4.573167000000000e+00 * tc[0]
            -2.168068000000000e-03 * tc[1]
            +2.457815000000000e-07 * tc[2]
            -1.957420000000000e-11 * tc[3]
            +7.158270000000000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 9: AR */
        species[9] =
            -7.453750000000000e+02 * invT
            -1.866001000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            -7.453750000000000e+02 * invT
            +1.584651100000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.426835000000000e+04 * invT
            -3.083140000000000e+00
            -3.025078000000000e+00 * tc[0]
            -7.213445000000000e-04 * tc[1]
            +9.384713333333334e-08 * tc[2]
            -8.488174999999999e-12 * tc[3]
            +3.455476000000000e-16 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.896696000000000e+04 * invT
            +5.409018900000000e+00
            -4.453623000000000e+00 * tc[0]
            -1.570084500000000e-03 * tc[1]
            +2.130685000000000e-07 * tc[2]
            -1.994997500000000e-11 * tc[3]
            +8.345165000000000e-16 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void helmholtz(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.54716300e+04 * invT
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -1.01252100e+03 * invT
            +5.59221800e+00
            -3.29812400e+00 * tc[0]
            -4.12472100e-04 * tc[1]
            +1.35716917e-07 * tc[2]
            +7.89619500e-12 * tc[3]
            -2.06743600e-14 * tc[4];
        /*species 2: O */
        species[2] =
            +2.91476400e+04 * invT
            -1.01756600e+00
            -2.94642900e+00 * tc[0]
            +8.19083000e-04 * tc[1]
            -4.03505333e-07 * tc[2]
            +1.33570250e-10 * tc[3]
            -1.94534800e-14 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.34630913e+03 * invT
            +3.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.02081100e+04 * invT
            -2.03391000e-01
            -3.38684200e+00 * tc[0]
            -1.73749100e-03 * tc[1]
            +1.05911600e-06 * tc[2]
            -5.80715083e-10 * tc[3]
            +1.25329400e-13 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.00524900e+03 * invT
            -3.82180200e+00
            -3.21293600e+00 * tc[0]
            -5.63743000e-04 * tc[1]
            +9.59358333e-08 * tc[2]
            -1.09489750e-10 * tc[3]
            +4.38427700e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.94808040e+02 * invT
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.76631500e+04 * invT
            -4.39660900e+00
            -3.38875400e+00 * tc[0]
            -3.28461300e-03 * tc[1]
            +2.47502167e-08 * tc[2]
            +3.85483833e-10 * tc[3]
            -1.23575750e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.02090000e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
        /*species 9: AR */
        species[9] =
            -7.45375000e+02 * invT
            -2.86600100e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            -7.45375000e+02 * invT
            +5.84651200e-01
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.43105400e+04 * invT
            -2.58644500e+00
            -3.26245200e+00 * tc[0]
            -7.55970500e-04 * tc[1]
            +6.46959167e-07 * tc[2]
            -4.65162000e-10 * tc[3]
            +1.23747550e-13 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.83731400e+04 * invT
            -8.91276500e+00
            -2.27572500e+00 * tc[0]
            -4.96103600e-03 * tc[1]
            +1.73485167e-06 * tc[2]
            -5.72223917e-10 * tc[3]
            +1.05864000e-13 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.54716300e+04 * invT
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -8.35034000e+02 * invT
            +3.34653300e+00
            -2.99142300e+00 * tc[0]
            -3.50032200e-04 * tc[1]
            +9.38971500e-09 * tc[2]
            +7.69298167e-13 * tc[3]
            -7.91376000e-17 * tc[4];
        /*species 2: O */
        species[2] =
            +2.92308000e+04 * invT
            -3.37824800e+00
            -2.54206000e+00 * tc[0]
            +1.37753100e-05 * tc[1]
            +5.17133833e-10 * tc[2]
            -3.79255583e-13 * tc[3]
            +2.18402600e-17 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.68362875e+03 * invT
            -3.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
        /*species 4: H2O */
        species[4] =
            -2.98992100e+04 * invT
            -5.19067100e+00
            -2.67214600e+00 * tc[0]
            -1.52814650e-03 * tc[1]
            +1.45504333e-07 * tc[2]
            -1.00083000e-11 * tc[3]
            +3.19580900e-16 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.23393000e+03 * invT
            -4.91588000e-01
            -3.69757800e+00 * tc[0]
            -3.06759850e-04 * tc[1]
            +2.09807000e-08 * tc[2]
            -1.47940083e-12 * tc[3]
            +5.68217500e-17 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.11856713e+02 * invT
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.80069600e+04 * invT
            +3.07203000e+00
            -4.57316700e+00 * tc[0]
            -2.16806800e-03 * tc[1]
            +2.45781500e-07 * tc[2]
            -1.95742000e-11 * tc[3]
            +7.15827000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 9: AR */
        species[9] =
            -7.45375000e+02 * invT
            -2.86600100e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            -7.45375000e+02 * invT
            +5.84651100e-01
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.42683500e+04 * invT
            -4.08314000e+00
            -3.02507800e+00 * tc[0]
            -7.21344500e-04 * tc[1]
            +9.38471333e-08 * tc[2]
            -8.48817500e-12 * tc[3]
            +3.45547600e-16 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.89669600e+04 * invT
            +4.40901890e+00
            -4.45362300e+00 * tc[0]
            -1.57008450e-03 * tc[1]
            +2.13068500e-07 * tc[2]
            -1.99499750e-11 * tc[3]
            +8.34516500e-16 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesInternalEnergy(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +2.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 2: O */
        species[2] =
            +1.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 3: OH */
        species[3] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +2.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 5: O2 */
        species[5] =
            +2.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +2.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
        /*species 9: AR */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 10: HE */
        species[10] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 11: CO */
        species[11] =
            +2.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +1.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +1.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 2: O */
        species[2] =
            +1.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 3: OH */
        species[3] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +1.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 5: O2 */
        species[5] =
            +2.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 9: AR */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 10: HE */
        species[10] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 11: CO */
        species[11] =
            +2.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +3.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEnthalpy(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 2: O */
        species[2] =
            +2.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 3: OH */
        species[3] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +3.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 5: O2 */
        species[5] =
            +3.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
        /*species 9: AR */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 10: HE */
        species[10] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 11: CO */
        species[11] =
            +3.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +2.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 2: O */
        species[2] =
            +2.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 3: OH */
        species[3] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +2.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 5: O2 */
        species[5] =
            +3.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 9: AR */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 10: HE */
        species[10] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 11: CO */
        species[11] =
            +3.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +4.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEntropy(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00 * tc[0]
            +8.24944200e-04 * tc[1]
            -4.07150750e-07 * tc[2]
            -3.15847800e-11 * tc[3]
            +1.03371800e-13 * tc[4]
            -3.29409400e+00 ;
        /*species 2: O */
        species[2] =
            +2.94642900e+00 * tc[0]
            -1.63816600e-03 * tc[1]
            +1.21051600e-06 * tc[2]
            -5.34281000e-10 * tc[3]
            +9.72674000e-14 * tc[4]
            +2.96399500e+00 ;
        /*species 3: OH */
        species[3] =
            +4.12530561e+00 * tc[0]
            -3.22544939e-03 * tc[1]
            +3.26382346e-06 * tc[2]
            -1.93284548e-09 * tc[3]
            +5.15593447e-13 * tc[4]
            -6.90432960e-01 ;
        /*species 4: H2O */
        species[4] =
            +3.38684200e+00 * tc[0]
            +3.47498200e-03 * tc[1]
            -3.17734800e-06 * tc[2]
            +2.32286033e-09 * tc[3]
            -6.26647000e-13 * tc[4]
            +2.59023300e+00 ;
        /*species 5: O2 */
        species[5] =
            +3.21293600e+00 * tc[0]
            +1.12748600e-03 * tc[1]
            -2.87807500e-07 * tc[2]
            +4.37959000e-10 * tc[3]
            -2.19213850e-13 * tc[4]
            +6.03473800e+00 ;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00 * tc[0]
            -4.74912051e-03 * tc[1]
            +1.05791445e-05 * tc[2]
            -8.09212980e-09 * tc[3]
            +2.32306281e-12 * tc[4]
            +3.71666245e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00 * tc[0]
            +6.56922600e-03 * tc[1]
            -7.42506500e-08 * tc[2]
            -1.54193533e-09 * tc[3]
            +6.17878750e-13 * tc[4]
            +6.78536300e+00 ;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00 * tc[0]
            +1.40824000e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213750e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 9: AR */
        species[9] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600100e+00 ;
        /*species 10: HE */
        species[10] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +9.15348800e-01 ;
        /*species 11: CO */
        species[11] =
            +3.26245200e+00 * tc[0]
            +1.51194100e-03 * tc[1]
            -1.94087750e-06 * tc[2]
            +1.86064800e-09 * tc[3]
            -6.18737750e-13 * tc[4]
            +4.84889700e+00 ;
        /*species 12: CO2 */
        species[12] =
            +2.27572500e+00 * tc[0]
            +9.92207200e-03 * tc[1]
            -5.20455500e-06 * tc[2]
            +2.28889567e-09 * tc[3]
            -5.29320000e-13 * tc[4]
            +1.01884900e+01 ;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00 * tc[0]
            +7.00064400e-04 * tc[1]
            -2.81691450e-08 * tc[2]
            -3.07719267e-12 * tc[3]
            +3.95688000e-16 * tc[4]
            -1.35511000e+00 ;
        /*species 2: O */
        species[2] =
            +2.54206000e+00 * tc[0]
            -2.75506200e-05 * tc[1]
            -1.55140150e-09 * tc[2]
            +1.51702233e-12 * tc[3]
            -1.09201300e-16 * tc[4]
            +4.92030800e+00 ;
        /*species 3: OH */
        species[3] =
            +2.86472886e+00 * tc[0]
            +1.05650448e-03 * tc[1]
            -1.29541379e-07 * tc[2]
            +1.01739558e-11 * tc[3]
            -3.32989690e-16 * tc[4]
            +5.70164073e+00 ;
        /*species 4: H2O */
        species[4] =
            +2.67214600e+00 * tc[0]
            +3.05629300e-03 * tc[1]
            -4.36513000e-07 * tc[2]
            +4.00332000e-11 * tc[3]
            -1.59790450e-15 * tc[4]
            +6.86281700e+00 ;
        /*species 5: O2 */
        species[5] =
            +3.69757800e+00 * tc[0]
            +6.13519700e-04 * tc[1]
            -6.29421000e-08 * tc[2]
            +5.91760333e-12 * tc[3]
            -2.84108750e-16 * tc[4]
            +3.18916600e+00 ;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00 * tc[0]
            +2.23982013e-03 * tc[1]
            -3.16829075e-07 * tc[2]
            +3.80821233e-11 * tc[3]
            -2.69771337e-15 * tc[4]
            +3.78510215e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00 * tc[0]
            +4.33613600e-03 * tc[1]
            -7.37344500e-07 * tc[2]
            +7.82968000e-11 * tc[3]
            -3.57913500e-15 * tc[4]
            +5.01137000e-01 ;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00 * tc[0]
            +1.48797700e-03 * tc[1]
            -2.84238050e-07 * tc[2]
            +3.36568000e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 9: AR */
        species[9] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600100e+00 ;
        /*species 10: HE */
        species[10] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +9.15348900e-01 ;
        /*species 11: CO */
        species[11] =
            +3.02507800e+00 * tc[0]
            +1.44268900e-03 * tc[1]
            -2.81541400e-07 * tc[2]
            +3.39527000e-11 * tc[3]
            -1.72773800e-15 * tc[4]
            +6.10821800e+00 ;
        /*species 12: CO2 */
        species[12] =
            +4.45362300e+00 * tc[0]
            +3.14016900e-03 * tc[1]
            -6.39205500e-07 * tc[2]
            +7.97999000e-11 * tc[3]
            -4.17258250e-15 * tc[4]
            -9.55395900e-01 ;
    }
    return;
}


/*Returns the mean specific heat at CP (Eq. 33) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPBL(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  cpbl)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cpor[13]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 13; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPBS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  cpbs)
{
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cpor[13], tresult[13]; /* temporary storage */
    amrex::Real imw[13];

    get_imw(imw);

    cp_R(cpor, tc);

    for (int i = 0; i < 13; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 13; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBL(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  cvbl)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cvor[13]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 13; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  cvbs)
{
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cvor[13]; /* temporary storage */
    amrex::Real imw[13];

    get_imw(imw);

    cv_R(cvor, tc);

    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H */
    result += cvor[1]*y[1]*imw[1]; /*H2 */
    result += cvor[2]*y[2]*imw[2]; /*O */
    result += cvor[3]*y[3]*imw[3]; /*OH */
    result += cvor[4]*y[4]*imw[4]; /*H2O */
    result += cvor[5]*y[5]*imw[5]; /*O2 */
    result += cvor[6]*y[6]*imw[6]; /*HO2 */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*N2 */
    result += cvor[9]*y[9]*imw[9]; /*AR */
    result += cvor[10]*y[10]*imw[10]; /*HE */
    result += cvor[11]*y[11]*imw[11]; /*CO */
    result += cvor[12]*y[12]*imw[12]; /*CO2 */

    *cvbs = result * 8.31446261815324e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBML(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  hbml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real hml[13]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 13; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBMS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  hbms)
{
    amrex::Real result = 0;
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real hml[13], tmp[13]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[13];

    get_imw(imw);

    speciesEnthalpy(hml, tc);

    int id;
    for (id = 0; id < 13; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 13; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUBML(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  ubml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real uml[13]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 13; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUBMS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  ubms)
{
    amrex::Real result = 0;
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real ums[13]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[13];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H */
    result += y[1]*ums[1]*imw[1]; /*H2 */
    result += y[2]*ums[2]*imw[2]; /*O */
    result += y[3]*ums[3]*imw[3]; /*OH */
    result += y[4]*ums[4]*imw[4]; /*H2O */
    result += y[5]*ums[5]*imw[5]; /*O2 */
    result += y[6]*ums[6]*imw[6]; /*HO2 */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*N2 */
    result += y[9]*ums[9]*imw[9]; /*AR */
    result += y[10]*ums[10]*imw[10]; /*HE */
    result += y[11]*ums[11]*imw[11]; /*CO */
    result += y[12]*ums[12]*imw[12]; /*CO2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBML(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  sbml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    amrex::Real logPratio = log ( *P / 1013250.0 ); 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real sor[13]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 13; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31446261815324e+07;
}


/*get mixture entropy in mass units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBMS(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  sbms)
{
    amrex::Real result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    amrex::Real logPratio = log ( *P / 1013250.0 ); 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real sor[13]; /* temporary storage */
    amrex::Real x[13]; /* need a ytx conversion */
    amrex::Real YOW = 0; /*See Eq 4, 6 in CK Manual */
    amrex::Real imw[13];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(31.998800*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
    x[9] = y[9]/(39.948000*YOW); 
    x[10] = y[10]/(4.002600*YOW); 
    x[11] = y[11]/(28.010550*YOW); 
    x[12] = y[12]/(44.009950*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
    result += x[3]*(sor[3]-log((x[3]+1e-100))-logPratio);
    result += x[4]*(sor[4]-log((x[4]+1e-100))-logPratio);
    result += x[5]*(sor[5]-log((x[5]+1e-100))-logPratio);
    result += x[6]*(sor[6]-log((x[6]+1e-100))-logPratio);
    result += x[7]*(sor[7]-log((x[7]+1e-100))-logPratio);
    result += x[8]*(sor[8]-log((x[8]+1e-100))-logPratio);
    result += x[9]*(sor[9]-log((x[9]+1e-100))-logPratio);
    result += x[10]*(sor[10]-log((x[10]+1e-100))-logPratio);
    result += x[11]*(sor[11]-log((x[11]+1e-100))-logPratio);
    result += x[12]*(sor[12]-log((x[12]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.31446261815324e+07 * YOW;
}

/* get temperature given internal energy in mass units and mass fracs */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_T_GIVEN_EY(amrex::Real *  e, amrex::Real *  y, amrex::Real *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const amrex::Real tol  = 1.e-12;
#else
    const int maxiter = 200;
    const amrex::Real tol  = 1.e-6;
#endif
    amrex::Real ein  = *e;
    amrex::Real tmin = 90;/*max lower bound for thermo def */
    amrex::Real tmax = 4000;/*min upper bound for thermo def */
    amrex::Real e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */
    CKUBMS(&tmin, y, &emin);
    CKUBMS(&tmax, y, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, &cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,&e1);
        CKCVBS(&t1,y,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}

/* get temperature given enthalpy in mass units and mass fracs */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_T_GIVEN_HY(amrex::Real *  h, amrex::Real *  y, amrex::Real *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const amrex::Real tol  = 1.e-12;
#else
    const int maxiter = 200;
    const amrex::Real tol  = 1.e-6;
#endif
    amrex::Real hin  = *h;
    amrex::Real tmin = 90;/*max lower bound for thermo def */
    amrex::Real tmax = 4000;/*min upper bound for thermo def */
    amrex::Real h1,hmin,hmax,cp,t1,dt;
    int i;/* loop counter */
    CKHBMS(&tmin, y, &hmin);
    CKHBMS(&tmax, y, &hmax);
    if (hin < hmin) {
        /*Linear Extrapolation below tmin */
        CKCPBS(&tmin, y, &cp);
        *t = tmin - (hmin-hin)/cp;
        *ierr = 1;
        return;
    }
    if (hin > hmax) {
        /*Linear Extrapolation above tmax */
        CKCPBS(&tmax, y, &cp);
        *t = tmax - (hmax-hin)/cp;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKHBMS(&t1,y,&h1);
        CKCPBS(&t1,y,&cp);
        dt = (hin - h1) / cp;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}


/*Compute P = rhoRT/W(x) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPX(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  P)
{
    amrex::Real XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    *P = *rho * 8.31446261815324e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPY(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW = 0;/* for computing mean MW */
    amrex::Real imw[13];

    get_imw(imw);

    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    *P = *rho * 8.31446261815324e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(c) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPC(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  c,  amrex::Real *  P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    amrex::Real W = 0;
    amrex::Real sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*31.998800; /*O2 */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */
    W += c[9]*39.948000; /*AR */
    W += c[10]*4.002600; /*HE */
    W += c[11]*28.010550; /*CO */
    W += c[12]*44.009950; /*CO2 */

    for (id = 0; id < 13; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31446261815324e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOX(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  rho)
{
    amrex::Real XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    *rho = *P * XW / (8.31446261815324e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOY(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  rho)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[13];
    amrex::Real imw[13];

    get_imw(imw);

    for (int i = 0; i < 13; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 13; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31446261815324e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOC(amrex::Real *  P, amrex::Real *  T, amrex::Real *  c,  amrex::Real *  rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    amrex::Real W = 0;
    amrex::Real sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*31.998800; /*O2 */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */
    W += c[9]*39.948000; /*AR */
    W += c[10]*4.002600; /*HE */
    W += c[11]*28.010550; /*CO */
    W += c[12]*44.009950; /*CO2 */

    for (id = 0; id < 13; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31446261815324e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT( amrex::Real *  wt)
{
    get_mw(wt);
}

/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWY(amrex::Real *  y,  amrex::Real *  wtm)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[13];
    amrex::Real imw[13];

    get_imw(imw);

    for (int i = 0; i < 13; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 13; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWX(amrex::Real *  x,  amrex::Real *  wtm)
{
    amrex::Real XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    *wtm = XW;

    return;
}

/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWC(amrex::Real *  c,  amrex::Real *  wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    amrex::Real W = 0;
    amrex::Real sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*31.998800; /*O2 */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */
    W += c[9]*39.948000; /*AR */
    W += c[10]*4.002600; /*HE */
    W += c[11]*28.010550; /*CO */
    W += c[12]*44.009950; /*CO2 */

    for (id = 0; id < 13; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPOR(amrex::Real *  T, amrex::Real *  cpor)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHORT(amrex::Real *  T, amrex::Real *  hort)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSOR(amrex::Real *  T, amrex::Real *  sor)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}

/*convert y[species] (mass fracs) to x[species] (mole fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTX(amrex::Real *  y,  amrex::Real *  x)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[13];
    amrex::Real imw[13];

    get_imw(imw);

    for (int i = 0; i < 13; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 13; i++)
    {
        YOW += tmp[i];
    }

    amrex::Real YOWINV = 1.0/YOW;

    for (int i = 0; i < 13; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  c)
{
    amrex::Real YOW = 0;
    amrex::Real PWORT;
    amrex::Real imw[13];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 13; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 13; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446261815324e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 13; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCR(amrex::Real *  rho, amrex::Real * /*T*/, amrex::Real * y,  amrex::Real * c)
{
    amrex::Real imw[13];

    get_imw(imw);

    for (int i = 0; i < 13; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTY(amrex::Real *  x,  amrex::Real *  y)
{
    amrex::Real XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    /*Now compute conversion */
    amrex::Real XWinv = 1.0/XW;
    y[0] = x[0]*1.007970*XWinv; 
    y[1] = x[1]*2.015940*XWinv; 
    y[2] = x[2]*15.999400*XWinv; 
    y[3] = x[3]*17.007370*XWinv; 
    y[4] = x[4]*18.015340*XWinv; 
    y[5] = x[5]*31.998800*XWinv; 
    y[6] = x[6]*33.006770*XWinv; 
    y[7] = x[7]*34.014740*XWinv; 
    y[8] = x[8]*28.013400*XWinv; 
    y[9] = x[9]*39.948000*XWinv; 
    y[10] = x[10]*4.002600*XWinv; 
    y[11] = x[11]*28.010550*XWinv; 
    y[12] = x[12]*44.009950*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  c)
{
    int id; /*loop counter */
    amrex::Real PORT = (*P)/(8.31446261815324e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 13; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCR(amrex::Real *  rho, amrex::Real * /*T*/, amrex::Real *  x, amrex::Real *  c)
{
    int id; /*loop counter */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 13; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTX(amrex::Real *  c, amrex::Real *  x)
{
    int id; /*loop counter */
    amrex::Real sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 13; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    amrex::Real sumCinv = 1.0/sumC;
    for (id = 0; id < 13; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTY(amrex::Real *  c, amrex::Real *  y)
{
    amrex::Real CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*1.007970; /*H */
    CW += c[1]*2.015940; /*H2 */
    CW += c[2]*15.999400; /*O */
    CW += c[3]*17.007370; /*OH */
    CW += c[4]*18.015340; /*H2O */
    CW += c[5]*31.998800; /*O2 */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*28.013400; /*N2 */
    CW += c[9]*39.948000; /*AR */
    CW += c[10]*4.002600; /*HE */
    CW += c[11]*28.010550; /*CO */
    CW += c[12]*44.009950; /*CO2 */
    /*Now compute conversion */
    amrex::Real CWinv = 1.0/CW;
    y[0] = c[0]*1.007970*CWinv; 
    y[1] = c[1]*2.015940*CWinv; 
    y[2] = c[2]*15.999400*CWinv; 
    y[3] = c[3]*17.007370*CWinv; 
    y[4] = c[4]*18.015340*CWinv; 
    y[5] = c[5]*31.998800*CWinv; 
    y[6] = c[6]*33.006770*CWinv; 
    y[7] = c[7]*34.014740*CWinv; 
    y[8] = c[8]*28.013400*CWinv; 
    y[9] = c[9]*39.948000*CWinv; 
    y[10] = c[10]*4.002600*CWinv; 
    y[11] = c[11]*28.010550*CWinv; 
    y[12] = c[12]*44.009950*CWinv; 

    return;
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVML(amrex::Real *  T,  amrex::Real *  cvml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        cvml[id] *= 8.31446261815324e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPML(amrex::Real *  T,  amrex::Real *  cpml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        cpml[id] *= 8.31446261815324e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUML(amrex::Real *  T,  amrex::Real *  uml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHML(amrex::Real *  T,  amrex::Real *  hml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        hml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSML(amrex::Real *  T,  amrex::Real *  sml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        sml[id] *= 8.31446261815324e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVMS(amrex::Real *  T,  amrex::Real *  cvms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 8.248720317224957e+07; /*H */
    cvms[1] *= 4.124360158612479e+07; /*H2 */
    cvms[2] *= 5.196734013871295e+06; /*O */
    cvms[3] *= 4.888740950630956e+06; /*OH */
    cvms[4] *= 4.615212712140454e+06; /*H2O */
    cvms[5] *= 2.598367006935648e+06; /*O2 */
    cvms[6] *= 2.519017346487778e+06; /*HO2 */
    cvms[7] *= 2.444370475315478e+06; /*H2O2 */
    cvms[8] *= 2.968030520448514e+06; /*N2 */
    cvms[9] *= 2.081321372322329e+06; /*AR */
    cvms[10] *= 2.077265432007505e+07; /*HE */
    cvms[11] *= 2.968332509769797e+06; /*CO */
    cvms[12] *= 1.889223372931176e+06; /*CO2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPMS(amrex::Real *  T,  amrex::Real *  cpms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 8.248720317224957e+07; /*H */
    cpms[1] *= 4.124360158612479e+07; /*H2 */
    cpms[2] *= 5.196734013871295e+06; /*O */
    cpms[3] *= 4.888740950630956e+06; /*OH */
    cpms[4] *= 4.615212712140454e+06; /*H2O */
    cpms[5] *= 2.598367006935648e+06; /*O2 */
    cpms[6] *= 2.519017346487778e+06; /*HO2 */
    cpms[7] *= 2.444370475315478e+06; /*H2O2 */
    cpms[8] *= 2.968030520448514e+06; /*N2 */
    cpms[9] *= 2.081321372322329e+06; /*AR */
    cpms[10] *= 2.077265432007505e+07; /*HE */
    cpms[11] *= 2.968332509769797e+06; /*CO */
    cpms[12] *= 1.889223372931176e+06; /*CO2 */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUMS(amrex::Real *  T,  amrex::Real *  ums)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[13];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    for (int i = 0; i < 13; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHMS(amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[13];

    get_imw(imw);

    speciesEnthalpy(hms, tc);

    for (int i = 0; i < 13; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSMS(amrex::Real *  T,  amrex::Real *  sms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 8.248720317224957e+07; /*H */
    sms[1] *= 4.124360158612479e+07; /*H2 */
    sms[2] *= 5.196734013871295e+06; /*O */
    sms[3] *= 4.888740950630956e+06; /*OH */
    sms[4] *= 4.615212712140454e+06; /*H2O */
    sms[5] *= 2.598367006935648e+06; /*O2 */
    sms[6] *= 2.519017346487778e+06; /*HO2 */
    sms[7] *= 2.444370475315478e+06; /*H2O2 */
    sms[8] *= 2.968030520448514e+06; /*N2 */
    sms[9] *= 2.081321372322329e+06; /*AR */
    sms[10] *= 2.077265432007505e+07; /*HE */
    sms[11] *= 2.968332509769797e+06; /*CO */
    sms[12] *= 1.889223372931176e+06; /*CO2 */
}

#ifdef AMREX_USE_GPU

/* GPU version of productionRate: no more use of thermo namespace vectors  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qfqr(amrex::Real *  qf, amrex::Real * qr, amrex::Real * sc, amrex::Real * qss_sc, amrex::Real * tc, amrex::Real invT)
{

    /*reaction 15: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[0]*sc[5];
    qr[0] = sc[6];

    /*reaction 22: H2O2 (+M) <=> OH + OH (+M) */
    qf[1] = sc[7];
    qr[1] = sc[3]*sc[3];

    /*reaction 6: H2 + M <=> H + H + M */
    qf[2] = sc[1];
    qr[2] = sc[0]*sc[0];

    /*reaction 9: O + O + M <=> O2 + M */
    qf[3] = sc[2]*sc[2];
    qr[3] = sc[5];

    /*reaction 12: O + H + M <=> OH + M */
    qf[4] = sc[0]*sc[2];
    qr[4] = sc[3];

    /*reaction 13: H2O + M <=> H + OH + M */
    qf[5] = sc[4];
    qr[5] = sc[0]*sc[3];

    /*reaction 1: H + O2 <=> O + OH */
    qf[6] = sc[0]*sc[5];
    qr[6] = sc[2]*sc[3];

    /*reaction 2: O + H2 <=> H + OH */
    qf[7] = sc[1]*sc[2];
    qr[7] = sc[0]*sc[3];

    /*reaction 3: O + H2 <=> H + OH */
    qf[8] = sc[1]*sc[2];
    qr[8] = sc[0]*sc[3];

    /*reaction 4: H2 + OH <=> H2O + H */
    qf[9] = sc[1]*sc[3];
    qr[9] = sc[0]*sc[4];

    /*reaction 5: OH + OH <=> O + H2O */
    qf[10] = sc[3]*sc[3];
    qr[10] = sc[2]*sc[4];

    /*reaction 7: H2 + AR <=> H + H + AR */
    qf[11] = sc[1]*sc[9];
    qr[11] = sc[0]*sc[0]*sc[9];

    /*reaction 8: H2 + HE <=> H + H + HE */
    qf[12] = sc[1]*sc[10];
    qr[12] = sc[0]*sc[0]*sc[10];

    /*reaction 10: O + O + AR <=> O2 + AR */
    qf[13] = sc[2]*sc[2]*sc[9];
    qr[13] = sc[5]*sc[9];

    /*reaction 11: O + O + HE <=> O2 + HE */
    qf[14] = sc[2]*sc[2]*sc[10];
    qr[14] = sc[5]*sc[10];

    /*reaction 14: H2O + H2O <=> H + OH + H2O */
    qf[15] = sc[4]*sc[4];
    qr[15] = sc[0]*sc[3]*sc[4];

    /*reaction 16: HO2 + H <=> H2 + O2 */
    qf[16] = sc[0]*sc[6];
    qr[16] = sc[1]*sc[5];

    /*reaction 17: HO2 + H <=> OH + OH */
    qf[17] = sc[0]*sc[6];
    qr[17] = sc[3]*sc[3];

    /*reaction 18: HO2 + O <=> O2 + OH */
    qf[18] = sc[2]*sc[6];
    qr[18] = sc[3]*sc[5];

    /*reaction 19: HO2 + OH <=> H2O + O2 */
    qf[19] = sc[3]*sc[6];
    qr[19] = sc[4]*sc[5];

    /*reaction 20: HO2 + HO2 <=> H2O2 + O2 */
    qf[20] = sc[6]*sc[6];
    qr[20] = sc[5]*sc[7];

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    qf[21] = sc[6]*sc[6];
    qr[21] = sc[5]*sc[7];

    /*reaction 23: H2O2 + H <=> H2O + OH */
    qf[22] = sc[0]*sc[7];
    qr[22] = sc[3]*sc[4];

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    qf[23] = sc[0]*sc[7];
    qr[23] = sc[1]*sc[6];

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    qf[24] = sc[2]*sc[7];
    qr[24] = sc[3]*sc[6];

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    qf[25] = sc[3]*sc[7];
    qr[25] = sc[4]*sc[6];

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    qf[26] = sc[3]*sc[7];
    qr[26] = sc[4]*sc[6];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 13; ++i) {
        mixture += sc[i];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[13];
    gibbs(g_RT, tc);

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 * invT;
    amrex::Real refCinv = 1 / refC;

    /* Evaluate the kfs */
    amrex::Real k_f, k_r, Corr;
    amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

    // (0):  H + O2 <=> O + OH
    k_f = 1.0000000000000002e-06 * 104000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (15286) * invT);
    Corr  = 1.0;
    qf[6] *= Corr * k_f;
    qr[6] *= Corr * k_f / exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    // (1):  O + H2 <=> H + OH
    k_f = 1.0000000000000002e-06 * 3818000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (7948) * invT);
    Corr  = 1.0;
    qf[7] *= Corr * k_f;
    qr[7] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    // (2):  O + H2 <=> H + OH
    k_f = 1.0000000000000002e-06 * 879200000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (19170) * invT);
    Corr  = 1.0;
    qf[8] *= Corr * k_f;
    qr[8] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    // (3):  H2 + OH <=> H2O + H
    k_f = 1.0000000000000002e-06 * 216000000 
               * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    Corr  = 1.0;
    qf[9] *= Corr * k_f;
    qr[9] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4]);
    // (4):  OH + OH <=> O + H2O
    k_f = 1.0000000000000002e-06 * 33400 
               * exp(2.4199999999999999 * tc[0] - 0.50321666580471969 * (-1930) * invT);
    Corr  = 1.0;
    qf[10] *= Corr * k_f;
    qr[10] *= Corr * k_f / exp(-g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4]);
    // (5):  H2 + M <=> H + H + M
    k_f = 1.0000000000000002e-06 * 4.577e+19 
               * exp(-1.3999999999999999 * tc[0] - 0.50321666580471969 * (104380) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[10];
    qf[2] *= Corr * k_f;
    qr[2] *= Corr * k_f / (exp(-g_RT[0] - g_RT[0] + g_RT[1]) * refC);
    // (6):  H2 + AR <=> H + H + AR
    k_f = 1.0000000000000002e-06 * 5.84e+18 
               * exp(-1.1000000000000001 * tc[0] - 0.50321666580471969 * (104380) * invT);
    Corr  = 1.0;
    qf[11] *= Corr * k_f;
    qr[11] *= Corr * k_f / (exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9]) * refC);
    // (7):  H2 + HE <=> H + H + HE
    k_f = 1.0000000000000002e-06 * 5.84e+18 
               * exp(-1.1000000000000001 * tc[0] - 0.50321666580471969 * (104380) * invT);
    Corr  = 1.0;
    qf[12] *= Corr * k_f;
    qr[12] *= Corr * k_f / (exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10]) * refC);
    // (8):  O + O + M <=> O2 + M
    k_f = 1.0000000000000002e-12 * 6165000000000000 
               * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[10] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    qf[3] *= Corr * k_f;
    qr[3] *= Corr * k_f / (exp(g_RT[2] + g_RT[2] - g_RT[5]) * refCinv);
    // (9):  O + O + AR <=> O2 + AR
    k_f = 1.0000000000000002e-12 * 18860000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1788) * invT);
    Corr  = 1.0;
    qf[13] *= Corr * k_f;
    qr[13] *= Corr * k_f / (exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9]) * refCinv);
    // (10):  O + O + HE <=> O2 + HE
    k_f = 1.0000000000000002e-12 * 18860000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1788) * invT);
    Corr  = 1.0;
    qf[14] *= Corr * k_f;
    qr[14] *= Corr * k_f / (exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10]) * refCinv);
    // (11):  O + H + M <=> OH + M
    k_f = 1.0000000000000002e-12 * 4.714e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0.75 - 1)*sc[9] + ( 0.75 - 1)*sc[10] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    qf[4] *= Corr * k_f;
    qr[4] *= Corr * k_f / (exp(g_RT[0] + g_RT[2] - g_RT[3]) * refCinv);
    // (12):  H2O + M <=> H + OH + M
    k_f = 1.0000000000000002e-06 * 6.0640000000000002e+27 
               * exp(-3.3220000000000001 * tc[0] - 0.50321666580471969 * (120790) * invT);
    Corr  = mixture + ( 3 - 1)*sc[1] + ( 0 - 1)*sc[4] + ( 1.1000000000000001 - 1)*sc[10] + ( 2 - 1)*sc[8] + ( 1.5 - 1)*sc[5] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    qf[5] *= Corr * k_f;
    qr[5] *= Corr * k_f / (exp(-g_RT[0] - g_RT[3] + g_RT[4]) * refC);
    // (13):  H2O + H2O <=> H + OH + H2O
    k_f = 1.0000000000000002e-06 * 1.006e+26 
               * exp(-2.4399999999999999 * tc[0] - 0.50321666580471969 * (120180) * invT);
    Corr  = 1.0;
    qf[15] *= Corr * k_f;
    qr[15] *= Corr * k_f / (exp(-g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4]) * refC);
    // (14):  H + O2 (+M) <=> HO2 (+M)
    k_f = 1.0000000000000002e-06 * 4650840000000 
               * exp(0.44 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[1] + ( 14 - 1)*sc[4] + ( 0.78000000000000003 - 1)*sc[5] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12] + ( 0.67000000000000004 - 1)*sc[9] + ( 0.80000000000000004 - 1)*sc[10];
    redP = Corr / k_f * 1e-12 * 6.366e+20 
               * exp(-1.72  * tc[0] - 0.50321666580471969  * (524.79999999999995) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.5)*exp(-tc[1] / 1.0000000000000001e-30) 
        + 0.5 * exp(-tc[1]/1e+30)  
        + 0.);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[0] *= Corr * k_f;
    qr[0] *= Corr * k_f / (exp(g_RT[0] + g_RT[5] - g_RT[6]) * refCinv);
    // (15):  HO2 + H <=> H2 + O2
    k_f = 1.0000000000000002e-06 * 2750000 
               * exp(2.0899999999999999 * tc[0] - 0.50321666580471969 * (-1451) * invT);
    Corr  = 1.0;
    qf[16] *= Corr * k_f;
    qr[16] *= Corr * k_f / exp(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6]);
    // (16):  HO2 + H <=> OH + OH
    k_f = 1.0000000000000002e-06 * 70790000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (295) * invT);
    Corr  = 1.0;
    qf[17] *= Corr * k_f;
    qr[17] *= Corr * k_f / exp(g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6]);
    // (17):  HO2 + O <=> O2 + OH
    k_f = 1.0000000000000002e-06 * 28500000000 
               * exp(1 * tc[0] - 0.50321666580471969 * (-723.92999999999995) * invT);
    Corr  = 1.0;
    qf[18] *= Corr * k_f;
    qr[18] *= Corr * k_f / exp(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6]);
    // (18):  HO2 + OH <=> H2O + O2
    k_f = 1.0000000000000002e-06 * 28900000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-497) * invT);
    Corr  = 1.0;
    qf[19] *= Corr * k_f;
    qr[19] *= Corr * k_f / exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6]);
    // (19):  HO2 + HO2 <=> H2O2 + O2
    k_f = 1.0000000000000002e-06 * 420000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (11982) * invT);
    Corr  = 1.0;
    qf[20] *= Corr * k_f;
    qr[20] *= Corr * k_f / exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    // (20):  HO2 + HO2 <=> H2O2 + O2
    k_f = 1.0000000000000002e-06 * 130000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1629.3) * invT);
    Corr  = 1.0;
    qf[21] *= Corr * k_f;
    qr[21] *= Corr * k_f / exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    // (21):  H2O2 (+M) <=> OH + OH (+M)
    k_f = 1 * 2000000000000 
               * exp(0.90000000000000002 * tc[0] - 0.50321666580471969 * (48749) * invT);
    Corr  = mixture + ( 7.5 - 1)*sc[4] + ( 1.6000000000000001 - 1)*sc[12] + ( 1.5 - 1)*sc[8] + ( 1.2 - 1)*sc[5] + ( 0.65000000000000002 - 1)*sc[10] + ( 7.7000000000000002 - 1)*sc[7] + ( 3.7000000000000002 - 1)*sc[1] + ( 2.7999999999999998 - 1)*sc[11];
    redP = Corr / k_f * 1e-6 * 2.49e+24 
               * exp(-2.2999999999999998  * tc[0] - 0.50321666580471969  * (48749) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.42999999999999999)*exp(-tc[1] / 1.0000000000000001e-30) 
        + 0.42999999999999999 * exp(-tc[1]/1e+30)  
        + 0.);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[1] *= Corr * k_f;
    qr[1] *= Corr * k_f / (exp(-g_RT[3] - g_RT[3] + g_RT[7]) * refC);
    // (22):  H2O2 + H <=> H2O + OH
    k_f = 1.0000000000000002e-06 * 24100000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (3970) * invT);
    Corr  = 1.0;
    qf[22] *= Corr * k_f;
    qr[22] *= Corr * k_f / exp(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7]);
    // (23):  H2O2 + H <=> HO2 + H2
    k_f = 1.0000000000000002e-06 * 48200000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (7950) * invT);
    Corr  = 1.0;
    qf[23] *= Corr * k_f;
    qr[23] *= Corr * k_f / exp(g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7]);
    // (24):  H2O2 + O <=> OH + HO2
    k_f = 1.0000000000000002e-06 * 9550000 
               * exp(2 * tc[0] - 0.50321666580471969 * (3970) * invT);
    Corr  = 1.0;
    qf[24] *= Corr * k_f;
    qr[24] *= Corr * k_f / exp(g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7]);
    // (25):  H2O2 + OH <=> HO2 + H2O
    k_f = 1.0000000000000002e-06 * 1740000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (318) * invT);
    Corr  = 1.0;
    qf[25] *= Corr * k_f;
    qr[25] *= Corr * k_f / exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
    // (26):  H2O2 + OH <=> HO2 + H2O
    k_f = 1.0000000000000002e-06 * 75900000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (7270) * invT);
    Corr  = 1.0;
    qf[26] *= Corr * k_f;
    qr[26] *= Corr * k_f / exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);


    return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real * wdot, amrex::Real * sc, amrex::Real T)
{
    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];

    amrex::Real qdot, q_f[27], q_r[27];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 13; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[0] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[3] += qdot;
    wdot[3] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[2]-q_r[2];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;

    qdot = q_f[3]-q_r[3];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[0] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[0] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[7]-q_r[7];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[9] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[10] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[9] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[10] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[0] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[4] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    return;
}
#endif



/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC(amrex::Real *  T, amrex::Real *  C,  amrex::Real *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 13; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[13]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[13];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446261815324e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 
    c[7] = PWORT * y[7]*imw[7]; 
    c[8] = PWORT * y[8]*imw[8]; 
    c[9] = PWORT * y[9]*imw[9]; 
    c[10] = PWORT * y[10]*imw[10]; 
    c[11] = PWORT * y[11]*imw[11]; 
    c[12] = PWORT * y[12]*imw[12]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[13]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 13; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[13]; /*temporary storage */
    amrex::Real imw[13];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 
    c[7] = 1e6 * (*rho) * y[7]*imw[7]; 
    c[8] = 1e6 * (*rho) * y[8]*imw[8]; 
    c[9] = 1e6 * (*rho) * y[9]*imw[9]; 
    c[10] = 1e6 * (*rho) * y[10]*imw[10]; 
    c[11] = 1e6 * (*rho) * y[11]*imw[11]; 
    c[12] = 1e6 * (*rho) * y[12]*imw[12]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[13]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 13; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void dcvpRdT(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 1: H2 */
        species[1] =
            +8.24944200e-04
            -1.62860300e-06 * tc[1]
            -2.84263020e-10 * tc[2]
            +1.65394880e-12 * tc[3];
        /*species 2: O */
        species[2] =
            -1.63816600e-03
            +4.84206400e-06 * tc[1]
            -4.80852900e-09 * tc[2]
            +1.55627840e-12 * tc[3];
        /*species 3: OH */
        species[3] =
            -3.22544939e-03
            +1.30552938e-05 * tc[1]
            -1.73956093e-08 * tc[2]
            +8.24949516e-12 * tc[3];
        /*species 4: H2O */
        species[4] =
            +3.47498200e-03
            -1.27093920e-05 * tc[1]
            +2.09057430e-08 * tc[2]
            -1.00263520e-11 * tc[3];
        /*species 5: O2 */
        species[5] =
            +1.12748600e-03
            -1.15123000e-06 * tc[1]
            +3.94163100e-09 * tc[2]
            -3.50742160e-12 * tc[3];
        /*species 6: HO2 */
        species[6] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +6.56922600e-03
            -2.97002600e-07 * tc[1]
            -1.38774180e-08 * tc[2]
            +9.88606000e-12 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.40824000e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77942000e-12 * tc[3];
        /*species 9: AR */
        species[9] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 10: HE */
        species[10] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 11: CO */
        species[11] =
            +1.51194100e-03
            -7.76351000e-06 * tc[1]
            +1.67458320e-08 * tc[2]
            -9.89980400e-12 * tc[3];
        /*species 12: CO2 */
        species[12] =
            +9.92207200e-03
            -2.08182200e-05 * tc[1]
            +2.06000610e-08 * tc[2]
            -8.46912000e-12 * tc[3];
    } else {
        /*species 0: H */
        species[0] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 1: H2 */
        species[1] =
            +7.00064400e-04
            -1.12676580e-07 * tc[1]
            -2.76947340e-11 * tc[2]
            +6.33100800e-15 * tc[3];
        /*species 2: O */
        species[2] =
            -2.75506200e-05
            -6.20560600e-09 * tc[1]
            +1.36532010e-11 * tc[2]
            -1.74722080e-15 * tc[3];
        /*species 3: OH */
        species[3] =
            +1.05650448e-03
            -5.18165516e-07 * tc[1]
            +9.15656022e-11 * tc[2]
            -5.32783504e-15 * tc[3];
        /*species 4: H2O */
        species[4] =
            +3.05629300e-03
            -1.74605200e-06 * tc[1]
            +3.60298800e-10 * tc[2]
            -2.55664720e-14 * tc[3];
        /*species 5: O2 */
        species[5] =
            +6.13519700e-04
            -2.51768400e-07 * tc[1]
            +5.32584300e-11 * tc[2]
            -4.54574000e-15 * tc[3];
        /*species 6: HO2 */
        species[6] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.33613600e-03
            -2.94937800e-06 * tc[1]
            +7.04671200e-10 * tc[2]
            -5.72661600e-14 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.48797700e-03
            -1.13695220e-06 * tc[1]
            +3.02911200e-10 * tc[2]
            -2.70134040e-14 * tc[3];
        /*species 9: AR */
        species[9] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 10: HE */
        species[10] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 11: CO */
        species[11] =
            +1.44268900e-03
            -1.12616560e-06 * tc[1]
            +3.05574300e-10 * tc[2]
            -2.76438080e-14 * tc[3];
        /*species 12: CO2 */
        species[12] =
            +3.14016900e-03
            -2.55682200e-06 * tc[1]
            +7.18199100e-10 * tc[2]
            -6.67613200e-14 * tc[3];
    }
    return;
}

/*compute an approx to the reaction Jacobian */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void aJacobian_precond(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int HP)
{
    for (int i=0; i<196; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[13];
    for (int k=0; k<13; k++) {
        wdot[k] = 0.0;
    }

    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];
    amrex::Real invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 / T;
    amrex::Real refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int k = 0; k < 13; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[13];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[13];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[13];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[1] + ( 14 - 1)*sc[4] + ( 0.78000000000000003 - 1)*sc[5] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12] + ( 0.67000000000000004 - 1)*sc[9] + ( 0.80000000000000004 - 1)*sc[10];
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = 1.0000000000000002e-06 * 4650840000000
                * exp(0.44 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0.44 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* pressure-fall-off */
    k_0 = 6.366e+20 * exp(-1.72 * tc[0] - 0.50321666580471969 * (524.79999999999995) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -1.72 * invT + 0.50321666580471969 * (524.79999999999995) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.5)*exp(-T/1.0000000000000001e-30);
    Fcent2 = 0.5 * exp(-T/1e+30);
    Fcent3 = 0.;
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1.0000000000000001e-30
        -Fcent2/1e+30
        + 0.);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[5] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac + k_f*sc[5];
    dqdc[1] = 2*dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = 14*dcdc_fac;
    dqdc[5] = 0.78000000000000003*dcdc_fac + k_f*sc[0];
    dqdc[6] = dcdc_fac - k_r;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = 0.67000000000000004*dcdc_fac;
    dqdc[10] = 0.80000000000000004*dcdc_fac;
    dqdc[11] = 1.8999999999999999*dcdc_fac;
    dqdc[12] = 3.7999999999999998*dcdc_fac;
    for (int k=0; k<13; k++) {
        J[14*k+0] -= dqdc[k];
        J[14*k+5] -= dqdc[k];
        J[14*k+6] += dqdc[k];
    }
    J[182] -= dqdT; /* dwdot[H]/dT */
    J[187] -= dqdT; /* dwdot[O2]/dT */
    J[188] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 7.5 - 1)*sc[4] + ( 1.6000000000000001 - 1)*sc[12] + ( 1.5 - 1)*sc[8] + ( 1.2 - 1)*sc[5] + ( 0.65000000000000002 - 1)*sc[10] + ( 7.7000000000000002 - 1)*sc[7] + ( 3.7000000000000002 - 1)*sc[1] + ( 2.7999999999999998 - 1)*sc[11];
    /* forward */
    phi_f = sc[7];
    k_f = 1 * 2000000000000
                * exp(0.90000000000000002 * tc[0] - 0.50321666580471969 * (48749) * invT);
    dlnkfdT = 0.90000000000000002 * invT + 0.50321666580471969 *  (48749)  * invT2;
    /* pressure-fall-off */
    k_0 = 2.49e+24 * exp(-2.2999999999999998 * tc[0] - 0.50321666580471969 * (48749) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.2999999999999998 * invT + 0.50321666580471969 * (48749) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.42999999999999999)*exp(-T/1.0000000000000001e-30);
    Fcent2 = 0.42999999999999999 * exp(-T/1e+30);
    Fcent3 = 0.;
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1.0000000000000001e-30
        -Fcent2/1e+30
        + 0.);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = pow(sc[3], 2.000000);
    Kc = refC * exp(-g_RT[3] - g_RT[3] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2.000000*h_RT[3]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[3] += 2 * q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac;
    dqdc[1] = 3.7000000000000002*dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac - k_r*2.000000*sc[3];
    dqdc[4] = 7.5*dcdc_fac;
    dqdc[5] = 1.2*dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = 7.7000000000000002*dcdc_fac + k_f;
    dqdc[8] = 1.5*dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = 0.65000000000000002*dcdc_fac;
    dqdc[11] = 2.7999999999999998*dcdc_fac;
    dqdc[12] = 1.6000000000000001*dcdc_fac;
    for (int k=0; k<13; k++) {
        J[14*k+3] += 2 * dqdc[k];
        J[14*k+7] -= dqdc[k];
    }
    J[185] += 2 * dqdT; /* dwdot[OH]/dT */
    J[189] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[10];
    /* forward */
    phi_f = sc[1];
    k_f = 1.0000000000000002e-06 * 4.577e+19
                * exp(-1.3999999999999999 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.3999999999999999 * invT + 0.50321666580471969 *  (104380)  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000);
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1]) + (2.000000*h_RT[0]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = q_nocor - k_r*2.000000*sc[0];
    dqdc[1] = 2.5*q_nocor + k_f;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = 12*q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = 0.0;
    dqdc[10] = 0.0;
    dqdc[11] = 1.8999999999999999*q_nocor;
    dqdc[12] = 3.7999999999999998*q_nocor;
    for (int k=0; k<13; k++) {
        J[14*k+0] += 2 * dqdc[k];
        J[14*k+1] -= dqdc[k];
    }
    J[182] += 2 * dqdT; /* dwdot[H]/dT */
    J[183] -= dqdT; /* dwdot[H2]/dT */

    /*reaction 4: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[10] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[2], 2.000000);
    k_f = 1.0000000000000002e-12 * 6165000000000000
                * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.5 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = q_nocor;
    dqdc[1] = 2.5*q_nocor;
    dqdc[2] = q_nocor + k_f*2.000000*sc[2];
    dqdc[3] = q_nocor;
    dqdc[4] = 12*q_nocor;
    dqdc[5] = q_nocor - k_r;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = 0.0;
    dqdc[10] = 0.0;
    dqdc[11] = 1.8999999999999999*q_nocor;
    dqdc[12] = 3.7999999999999998*q_nocor;
    for (int k=0; k<13; k++) {
        J[14*k+2] += -2 * dqdc[k];
        J[14*k+5] += dqdc[k];
    }
    J[184] += -2 * dqdT; /* dwdot[O]/dT */
    J[187] += dqdT; /* dwdot[O2]/dT */

    /*reaction 5: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0.75 - 1)*sc[9] + ( 0.75 - 1)*sc[10] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = 1.0000000000000002e-12 * 4.714e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = q_nocor + k_f*sc[2];
    dqdc[1] = 2.5*q_nocor;
    dqdc[2] = q_nocor + k_f*sc[0];
    dqdc[3] = q_nocor - k_r;
    dqdc[4] = 12*q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = 0.75*q_nocor;
    dqdc[10] = 0.75*q_nocor;
    dqdc[11] = 1.8999999999999999*q_nocor;
    dqdc[12] = 3.7999999999999998*q_nocor;
    for (int k=0; k<13; k++) {
        J[14*k+0] -= dqdc[k];
        J[14*k+2] -= dqdc[k];
        J[14*k+3] += dqdc[k];
    }
    J[182] -= dqdT; /* dwdot[H]/dT */
    J[184] -= dqdT; /* dwdot[O]/dT */
    J[185] += dqdT; /* dwdot[OH]/dT */

    /*reaction 6: H2O + M <=> H + OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 3 - 1)*sc[1] + ( 0 - 1)*sc[4] + ( 1.1000000000000001 - 1)*sc[10] + ( 2 - 1)*sc[8] + ( 1.5 - 1)*sc[5] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    /* forward */
    phi_f = sc[4];
    k_f = 1.0000000000000002e-06 * 6.0640000000000002e+27
                * exp(-3.3220000000000001 * tc[0] - 0.50321666580471969 * (120790) * invT);
    dlnkfdT = -3.3220000000000001 * invT + 0.50321666580471969 *  (120790)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4]) + (h_RT[0] + h_RT[3]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = q_nocor - k_r*sc[3];
    dqdc[1] = 3*q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor - k_r*sc[0];
    dqdc[4] =  + k_f;
    dqdc[5] = 1.5*q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = 2*q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = 1.1000000000000001*q_nocor;
    dqdc[11] = 1.8999999999999999*q_nocor;
    dqdc[12] = 3.7999999999999998*q_nocor;
    for (int k=0; k<13; k++) {
        J[14*k+0] += dqdc[k];
        J[14*k+3] += dqdc[k];
        J[14*k+4] -= dqdc[k];
    }
    J[182] += dqdT; /* dwdot[H]/dT */
    J[185] += dqdT; /* dwdot[OH]/dT */
    J[186] -= dqdT; /* dwdot[H2O]/dT */

    /*reaction 7: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = 1.0000000000000002e-06 * 104000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (15286) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (15286)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[3];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[2] += dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[5] -= dqdci;                /* dwdot[O2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[3];
    J[28] -= dqdci;               /* dwdot[H]/d[O] */
    J[30] += dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[47] -= dqdci;               /* dwdot[O2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[70] -= dqdci;               /* dwdot[H]/d[O2] */
    J[72] += dqdci;               /* dwdot[O]/d[O2] */
    J[73] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] -= dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[184] += dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[187] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 8: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-06 * 3818000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7948) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (7948)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[16] -= dqdci;               /* dwdot[O]/d[H2] */
    J[17] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[28] += dqdci;               /* dwdot[H]/d[O] */
    J[29] -= dqdci;               /* dwdot[H2]/d[O] */
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 9: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-06 * 879200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (19170) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (19170)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[16] -= dqdci;               /* dwdot[O]/d[H2] */
    J[17] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[28] += dqdci;               /* dwdot[H]/d[O] */
    J[29] -= dqdci;               /* dwdot[H2]/d[O] */
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 10: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  (3430)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[4];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[0] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[3] -= dqdci;                /* dwdot[OH]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[17] -= dqdci;               /* dwdot[OH]/d[H2] */
    J[18] += dqdci;               /* dwdot[H2O]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[0];
    J[56] += dqdci;               /* dwdot[H]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 11: OH + OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[3], 2.000000);
    k_f = 1.0000000000000002e-06 * 33400
                * exp(2.4199999999999999 * tc[0] - 0.50321666580471969 * (-1930) * invT);
    dlnkfdT = 2.4199999999999999 * invT + 0.50321666580471969 *  (-1930)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= 2 * q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[30] += dqdci;               /* dwdot[O]/d[O] */
    J[31] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[32] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[3];
    J[44] += dqdci;               /* dwdot[O]/d[OH] */
    J[45] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[58] += dqdci;               /* dwdot[O]/d[H2O] */
    J[59] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[184] += dqdT;               /* dwdot[O]/dT */
    J[185] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 12: H2 + AR <=> H + H + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = 1.0000000000000002e-06 * 5.84e+18
                * exp(-1.1000000000000001 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.1000000000000001 * invT + 0.50321666580471969 *  (104380)  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000)*sc[9];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (2.000000*h_RT[0] + h_RT[9]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*2.000000*sc[0]*sc[9];
    J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[9];
    J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[1] - k_r*pow(sc[0], 2.000000);
    J[126] += 2 * dqdci;          /* dwdot[H]/d[AR] */
    J[127] -= dqdci;              /* dwdot[H2]/d[AR] */
    /* d()/dT */
    J[182] += 2 * dqdT;           /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 13: H2 + HE <=> H + H + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = 1.0000000000000002e-06 * 5.84e+18
                * exp(-1.1000000000000001 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.1000000000000001 * invT + 0.50321666580471969 *  (104380)  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000)*sc[10];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (2.000000*h_RT[0] + h_RT[10]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*2.000000*sc[0]*sc[10];
    J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[10];
    J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[HE] */
    dqdci =  + k_f*sc[1] - k_r*pow(sc[0], 2.000000);
    J[140] += 2 * dqdci;          /* dwdot[H]/d[HE] */
    J[141] -= dqdci;              /* dwdot[H2]/d[HE] */
    /* d()/dT */
    J[182] += 2 * dqdT;           /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 14: O + O + AR <=> O2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[2], 2.000000)*sc[9];
    k_f = 1.0000000000000002e-12 * 18860000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1788) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1788)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2] + h_RT[9]) + (h_RT[5] + h_RT[9]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2.000000*sc[2]*sc[9];
    J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[9];
    J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[AR] */
    dqdci =  + k_f*pow(sc[2], 2.000000) - k_r*sc[5];
    J[128] += -2 * dqdci;         /* dwdot[O]/d[AR] */
    J[131] += dqdci;              /* dwdot[O2]/d[AR] */
    /* d()/dT */
    J[184] += -2 * dqdT;          /* dwdot[O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 15: O + O + HE <=> O2 + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[2], 2.000000)*sc[10];
    k_f = 1.0000000000000002e-12 * 18860000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1788) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1788)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[10];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2] + h_RT[10]) + (h_RT[5] + h_RT[10]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2.000000*sc[2]*sc[10];
    J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[10];
    J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[HE] */
    dqdci =  + k_f*pow(sc[2], 2.000000) - k_r*sc[5];
    J[142] += -2 * dqdci;         /* dwdot[O]/d[HE] */
    J[145] += dqdci;              /* dwdot[O2]/d[HE] */
    /* d()/dT */
    J[184] += -2 * dqdT;          /* dwdot[O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-06 * 1.006e+26
                * exp(-2.4399999999999999 * tc[0] - 0.50321666580471969 * (120180) * invT);
    dlnkfdT = -2.4399999999999999 * invT + 0.50321666580471969 *  (120180)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3]*sc[4];
    Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[0] + h_RT[3] + h_RT[4]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3]*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] -= dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0]*sc[4];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  + k_f*2.000000*sc[4] - k_r*sc[0]*sc[3];
    J[56] += dqdci;               /* dwdot[H]/d[H2O] */
    J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[186] -= dqdT;               /* dwdot[H2O]/dT */

    /*reaction 17: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = 1.0000000000000002e-06 * 2750000
                * exp(2.0899999999999999 * tc[0] - 0.50321666580471969 * (-1451) * invT);
    dlnkfdT = 2.0899999999999999 * invT + 0.50321666580471969 *  (-1451)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[5] += dqdci;                /* dwdot[O2]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[14] -= dqdci;               /* dwdot[H]/d[H2] */
    J[15] += dqdci;               /* dwdot[H2]/d[H2] */
    J[19] += dqdci;               /* dwdot[O2]/d[H2] */
    J[20] -= dqdci;               /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[1];
    J[70] -= dqdci;               /* dwdot[H]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[183] += dqdT;               /* dwdot[H2]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = 1.0000000000000002e-06 * 70790000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (295) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (295)  * invT2;
    /* reverse */
    phi_r = pow(sc[3], 2.000000);
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (2.000000*h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += 2 * dqdci;            /* dwdot[OH]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[3];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[45] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[87] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[185] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 19: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = 1.0000000000000002e-06 * 28500000000
                * exp(1 * tc[0] - 0.50321666580471969 * (-723.92999999999995) * invT);
    dlnkfdT = 1 * invT + 0.50321666580471969 *  (-723.92999999999995)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    J[34] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[47] += dqdci;               /* dwdot[O2]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[3];
    J[72] -= dqdci;               /* dwdot[O]/d[O2] */
    J[73] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[86] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[87] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = 1.0000000000000002e-06 * 28900000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-497) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-497)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[47] += dqdci;               /* dwdot[O2]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[61] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[62] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[73] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[74] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (11982) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (11982)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[91] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[5];
    J[103] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[104] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[105] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[189] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1629.3) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1629.3)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[91] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[5];
    J[103] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[104] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[105] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[189] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 23: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = 1.0000000000000002e-06 * 24100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (3970)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[56] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[98] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[101] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = 1.0000000000000002e-06 * 48200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7950) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (7950)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[6];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[14] -= dqdci;               /* dwdot[H]/d[H2] */
    J[15] += dqdci;               /* dwdot[H2]/d[H2] */
    J[20] += dqdci;               /* dwdot[HO2]/d[H2] */
    J[21] -= dqdci;               /* dwdot[H2O2]/d[H2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[1];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[98] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[99] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[183] += dqdT;               /* dwdot[H2]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 9550000
                * exp(2 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (3970)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[6];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[3] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[34] += dqdci;               /* dwdot[HO2]/d[O] */
    J[35] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[86] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[87] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[100] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[101] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 1740000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (318) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (318)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[101] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 75900000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7270) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (7270)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[101] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    amrex::Real c_R[13], dcRdT[13], e_RT[13];
    amrex::Real * eh_RT;
    if (HP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 13; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[182+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 13; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 13; ++m) {
            dehmixdc += eh_RT[m]*J[k*14+m];
        }
        J[k*14+13] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[195] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

/*compute an approx to the reaction Jacobian (for preconditioning) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT_SIMPLIFIED(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * HP)
{
    amrex::Real c[13];

    for (int k=0; k<13; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<13; k++) {
        J[182+k] *= 1.e-6;
        J[k*14+13] *= 1.e6;
    }

    return;
}

#ifdef AMREX_USE_GPU
/*compute the reaction Jacobian on GPU */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void aJacobian(amrex::Real * J, amrex::Real * sc, amrex::Real T, int consP)
{

    for (int i=0; i<196; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[13];
    for (int k=0; k<13; k++) {
        wdot[k] = 0.0;
    }

    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];
    amrex::Real invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 / T;
    amrex::Real refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int k = 0; k < 13; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[13];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[13];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[13];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[1] + ( 14 - 1)*sc[4] + ( 0.78000000000000003 - 1)*sc[5] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12] + ( 0.67000000000000004 - 1)*sc[9] + ( 0.80000000000000004 - 1)*sc[10];
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = 1.0000000000000002e-06 * 4650840000000
                * exp(0.44 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0.44 * invT + 0.50321666580471969 *  0  * invT2;
    /* pressure-fall-off */
    k_0 = 6.366e+20 * exp(-1.72 * tc[0] - 0.50321666580471969 * (524.79999999999995) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -1.72 * invT + 0.50321666580471969 * (524.79999999999995) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.5)*exp(-T/1.0000000000000001e-30);
    Fcent2 = 0.5 * exp(-T/1e+30);
    Fcent3 = 0.;
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1.0000000000000001e-30
        -Fcent2/1e+30
        + 0.);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[5] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[5] -= dqdci;                /* dwdot[O2]/d[H] */
        J[6] += dqdci;                /* dwdot[HO2]/d[H] */
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[14] -= dqdci;               /* dwdot[H]/d[H2] */
        J[19] -= dqdci;               /* dwdot[O2]/d[H2] */
        J[20] += dqdci;               /* dwdot[HO2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (14 - 1)*dcdc_fac;
        J[56] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[61] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[O2] */
        dqdci = (0.78000000000000003 - 1)*dcdc_fac + k_f*sc[0];
        J[70] -= dqdci;               /* dwdot[H]/d[O2] */
        J[75] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[76] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[89] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
        /* d()/d[AR] */
        dqdci = (0.67000000000000004 - 1)*dcdc_fac;
        J[126] -= dqdci;              /* dwdot[H]/d[AR] */
        J[131] -= dqdci;              /* dwdot[O2]/d[AR] */
        J[132] += dqdci;              /* dwdot[HO2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0.80000000000000004 - 1)*dcdc_fac;
        J[140] -= dqdci;              /* dwdot[H]/d[HE] */
        J[145] -= dqdci;              /* dwdot[O2]/d[HE] */
        J[146] += dqdci;              /* dwdot[HO2]/d[HE] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*dcdc_fac;
        J[154] -= dqdci;              /* dwdot[H]/d[CO] */
        J[159] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[160] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*dcdc_fac;
        J[168] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[173] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[174] += dqdci;              /* dwdot[HO2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[5];
        dqdc[1] = 2*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = 14*dcdc_fac;
        dqdc[5] = 0.78000000000000003*dcdc_fac + k_f*sc[0];
        dqdc[6] = dcdc_fac - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = 0.67000000000000004*dcdc_fac;
        dqdc[10] = 0.80000000000000004*dcdc_fac;
        dqdc[11] = 1.8999999999999999*dcdc_fac;
        dqdc[12] = 3.7999999999999998*dcdc_fac;
        for (int k=0; k<13; k++) {
            J[14*k+0] -= dqdc[k];
            J[14*k+5] -= dqdc[k];
            J[14*k+6] += dqdc[k];
        }
    }
    J[182] -= dqdT; /* dwdot[H]/dT */
    J[187] -= dqdT; /* dwdot[O2]/dT */
    J[188] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 7.5 - 1)*sc[4] + ( 1.6000000000000001 - 1)*sc[12] + ( 1.5 - 1)*sc[8] + ( 1.2 - 1)*sc[5] + ( 0.65000000000000002 - 1)*sc[10] + ( 7.7000000000000002 - 1)*sc[7] + ( 3.7000000000000002 - 1)*sc[1] + ( 2.7999999999999998 - 1)*sc[11];
    /* forward */
    phi_f = sc[7];
    k_f = 1 * 2000000000000
                * exp(0.90000000000000002 * tc[0] - 0.50321666580471969 * (48749) * invT);
    dlnkfdT = 0.90000000000000002 * invT + 0.50321666580471969 *  48749  * invT2;
    /* pressure-fall-off */
    k_0 = 2.49e+24 * exp(-2.2999999999999998 * tc[0] - 0.50321666580471969 * (48749) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.2999999999999998 * invT + 0.50321666580471969 * (48749) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.42999999999999999)*exp(-T/1.0000000000000001e-30);
    Fcent2 = 0.42999999999999999 * exp(-T/1e+30);
    Fcent3 = 0.;
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1.0000000000000001e-30
        -Fcent2/1e+30
        + 0.);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = pow(sc[3], 2.000000);
    Kc = refC * exp(-g_RT[3] - g_RT[3] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2.000000*h_RT[3]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[3] += 2 * q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (3.7000000000000002 - 1)*dcdc_fac;
        J[17] += 2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[21] -= dqdci;               /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  - k_r*2.000000*sc[3];
        J[45] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
        J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (7.5 - 1)*dcdc_fac;
        J[59] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
        J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[O2] */
        dqdci = (1.2 - 1)*dcdc_fac;
        J[73] += 2 * dqdci;           /* dwdot[OH]/d[O2] */
        J[77] -= dqdci;               /* dwdot[H2O2]/d[O2] */
        /* d()/d[H2O2] */
        dqdci = (7.7000000000000002 - 1)*dcdc_fac + k_f;
        J[101] += 2 * dqdci;          /* dwdot[OH]/d[H2O2] */
        J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[N2] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[115] += 2 * dqdci;          /* dwdot[OH]/d[N2] */
        J[119] -= dqdci;              /* dwdot[H2O2]/d[N2] */
        /* d()/d[HE] */
        dqdci = (0.65000000000000002 - 1)*dcdc_fac;
        J[143] += 2 * dqdci;          /* dwdot[OH]/d[HE] */
        J[147] -= dqdci;              /* dwdot[H2O2]/d[HE] */
        /* d()/d[CO] */
        dqdci = (2.7999999999999998 - 1)*dcdc_fac;
        J[157] += 2 * dqdci;          /* dwdot[OH]/d[CO] */
        J[161] -= dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (1.6000000000000001 - 1)*dcdc_fac;
        J[171] += 2 * dqdci;          /* dwdot[OH]/d[CO2] */
        J[175] -= dqdci;              /* dwdot[H2O2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = 3.7000000000000002*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac - k_r*2.000000*sc[3];
        dqdc[4] = 7.5*dcdc_fac;
        dqdc[5] = 1.2*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = 7.7000000000000002*dcdc_fac + k_f;
        dqdc[8] = 1.5*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = 0.65000000000000002*dcdc_fac;
        dqdc[11] = 2.7999999999999998*dcdc_fac;
        dqdc[12] = 1.6000000000000001*dcdc_fac;
        for (int k=0; k<13; k++) {
            J[14*k+3] += 2 * dqdc[k];
            J[14*k+7] -= dqdc[k];
        }
    }
    J[185] += 2 * dqdT; /* dwdot[OH]/dT */
    J[189] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[10];
    /* forward */
    phi_f = sc[1];
    k_f = 1.0000000000000002e-06 * 4.577e+19
                * exp(-1.3999999999999999 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.3999999999999999 * invT + 0.50321666580471969 *  104380  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000);
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1]) + (2.000000*h_RT[0]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  - k_r*2.000000*sc[0];
        J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
        J[1] -= dqdci;                /* dwdot[H2]/d[H] */
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor + k_f;
        J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
        J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor;
        J[56] += 2 * dqdci;           /* dwdot[H]/d[H2O] */
        J[57] -= dqdci;               /* dwdot[H2]/d[H2O] */
        /* d()/d[AR] */
        dqdci = (0 - 1)*q_nocor;
        J[126] += 2 * dqdci;          /* dwdot[H]/d[AR] */
        J[127] -= dqdci;              /* dwdot[H2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0 - 1)*q_nocor;
        J[140] += 2 * dqdci;          /* dwdot[H]/d[HE] */
        J[141] -= dqdci;              /* dwdot[H2]/d[HE] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor;
        J[154] += 2 * dqdci;          /* dwdot[H]/d[CO] */
        J[155] -= dqdci;              /* dwdot[H2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[168] += 2 * dqdci;          /* dwdot[H]/d[CO2] */
        J[169] -= dqdci;              /* dwdot[H2]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor - k_r*2.000000*sc[0];
        dqdc[1] = 2.5*q_nocor + k_f;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = 12*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[11] = 1.8999999999999999*q_nocor;
        dqdc[12] = 3.7999999999999998*q_nocor;
        for (int k=0; k<13; k++) {
            J[14*k+0] += 2 * dqdc[k];
            J[14*k+1] -= dqdc[k];
        }
    }
    J[182] += 2 * dqdT; /* dwdot[H]/dT */
    J[183] -= dqdT; /* dwdot[H2]/dT */

    /*reaction 4: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[10] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[2], 2.000000);
    k_f = 1.0000000000000002e-12 * 6165000000000000
                * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.5 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor;
        J[16] += -2 * dqdci;          /* dwdot[O]/d[H2] */
        J[19] += dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[2];
        J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[33] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor;
        J[58] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        J[61] += dqdci;               /* dwdot[O2]/d[H2O] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        J[75] += dqdci;               /* dwdot[O2]/d[O2] */
        /* d()/d[AR] */
        dqdci = (0 - 1)*q_nocor;
        J[128] += -2 * dqdci;         /* dwdot[O]/d[AR] */
        J[131] += dqdci;              /* dwdot[O2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0 - 1)*q_nocor;
        J[142] += -2 * dqdci;         /* dwdot[O]/d[HE] */
        J[145] += dqdci;              /* dwdot[O2]/d[HE] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor;
        J[156] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[159] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[170] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[173] += dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = 2.5*q_nocor;
        dqdc[2] = q_nocor + k_f*2.000000*sc[2];
        dqdc[3] = q_nocor;
        dqdc[4] = 12*q_nocor;
        dqdc[5] = q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[11] = 1.8999999999999999*q_nocor;
        dqdc[12] = 3.7999999999999998*q_nocor;
        for (int k=0; k<13; k++) {
            J[14*k+2] += -2 * dqdc[k];
            J[14*k+5] += dqdc[k];
        }
    }
    J[184] += -2 * dqdT; /* dwdot[O]/dT */
    J[187] += dqdT; /* dwdot[O2]/dT */

    /*reaction 5: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0.75 - 1)*sc[9] + ( 0.75 - 1)*sc[10] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = 1.0000000000000002e-12 * 4.714e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[2] -= dqdci;                /* dwdot[O]/d[H] */
        J[3] += dqdci;                /* dwdot[OH]/d[H] */
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor;
        J[14] -= dqdci;               /* dwdot[H]/d[H2] */
        J[16] -= dqdci;               /* dwdot[O]/d[H2] */
        J[17] += dqdci;               /* dwdot[OH]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[0];
        J[28] -= dqdci;               /* dwdot[H]/d[O] */
        J[30] -= dqdci;               /* dwdot[O]/d[O] */
        J[31] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[42] -= dqdci;               /* dwdot[H]/d[OH] */
        J[44] -= dqdci;               /* dwdot[O]/d[OH] */
        J[45] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor;
        J[56] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[58] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
        /* d()/d[AR] */
        dqdci = (0.75 - 1)*q_nocor;
        J[126] -= dqdci;              /* dwdot[H]/d[AR] */
        J[128] -= dqdci;              /* dwdot[O]/d[AR] */
        J[129] += dqdci;              /* dwdot[OH]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0.75 - 1)*q_nocor;
        J[140] -= dqdci;              /* dwdot[H]/d[HE] */
        J[142] -= dqdci;              /* dwdot[O]/d[HE] */
        J[143] += dqdci;              /* dwdot[OH]/d[HE] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor;
        J[154] -= dqdci;              /* dwdot[H]/d[CO] */
        J[156] -= dqdci;              /* dwdot[O]/d[CO] */
        J[157] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[168] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[170] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[171] += dqdci;              /* dwdot[OH]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor + k_f*sc[2];
        dqdc[1] = 2.5*q_nocor;
        dqdc[2] = q_nocor + k_f*sc[0];
        dqdc[3] = q_nocor - k_r;
        dqdc[4] = 12*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = 0.75*q_nocor;
        dqdc[10] = 0.75*q_nocor;
        dqdc[11] = 1.8999999999999999*q_nocor;
        dqdc[12] = 3.7999999999999998*q_nocor;
        for (int k=0; k<13; k++) {
            J[14*k+0] -= dqdc[k];
            J[14*k+2] -= dqdc[k];
            J[14*k+3] += dqdc[k];
        }
    }
    J[182] -= dqdT; /* dwdot[H]/dT */
    J[184] -= dqdT; /* dwdot[O]/dT */
    J[185] += dqdT; /* dwdot[OH]/dT */

    /*reaction 6: H2O + M <=> H + OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 3 - 1)*sc[1] + ( 0 - 1)*sc[4] + ( 1.1000000000000001 - 1)*sc[10] + ( 2 - 1)*sc[8] + ( 1.5 - 1)*sc[5] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    /* forward */
    phi_f = sc[4];
    k_f = 1.0000000000000002e-06 * 6.0640000000000002e+27
                * exp(-3.3220000000000001 * tc[0] - 0.50321666580471969 * (120790) * invT);
    dlnkfdT = -3.3220000000000001 * invT + 0.50321666580471969 *  120790  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4]) + (h_RT[0] + h_RT[3]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  - k_r*sc[3];
        J[0] += dqdci;                /* dwdot[H]/d[H] */
        J[3] += dqdci;                /* dwdot[OH]/d[H] */
        J[4] -= dqdci;                /* dwdot[H2O]/d[H] */
        /* d()/d[H2] */
        dqdci = (3 - 1)*q_nocor;
        J[14] += dqdci;               /* dwdot[H]/d[H2] */
        J[17] += dqdci;               /* dwdot[OH]/d[H2] */
        J[18] -= dqdci;               /* dwdot[H2O]/d[H2] */
        /* d()/d[OH] */
        dqdci =  - k_r*sc[0];
        J[42] += dqdci;               /* dwdot[H]/d[OH] */
        J[45] += dqdci;               /* dwdot[OH]/d[OH] */
        J[46] -= dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor + k_f;
        J[56] += dqdci;               /* dwdot[H]/d[H2O] */
        J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
        J[60] -= dqdci;               /* dwdot[H2O]/d[H2O] */
        /* d()/d[O2] */
        dqdci = (1.5 - 1)*q_nocor;
        J[70] += dqdci;               /* dwdot[H]/d[O2] */
        J[73] += dqdci;               /* dwdot[OH]/d[O2] */
        J[74] -= dqdci;               /* dwdot[H2O]/d[O2] */
        /* d()/d[N2] */
        dqdci = (2 - 1)*q_nocor;
        J[112] += dqdci;              /* dwdot[H]/d[N2] */
        J[115] += dqdci;              /* dwdot[OH]/d[N2] */
        J[116] -= dqdci;              /* dwdot[H2O]/d[N2] */
        /* d()/d[HE] */
        dqdci = (1.1000000000000001 - 1)*q_nocor;
        J[140] += dqdci;              /* dwdot[H]/d[HE] */
        J[143] += dqdci;              /* dwdot[OH]/d[HE] */
        J[144] -= dqdci;              /* dwdot[H2O]/d[HE] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor;
        J[154] += dqdci;              /* dwdot[H]/d[CO] */
        J[157] += dqdci;              /* dwdot[OH]/d[CO] */
        J[158] -= dqdci;              /* dwdot[H2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[168] += dqdci;              /* dwdot[H]/d[CO2] */
        J[171] += dqdci;              /* dwdot[OH]/d[CO2] */
        J[172] -= dqdci;              /* dwdot[H2O]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor - k_r*sc[3];
        dqdc[1] = 3*q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor - k_r*sc[0];
        dqdc[4] =  + k_f;
        dqdc[5] = 1.5*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = 2*q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = 1.1000000000000001*q_nocor;
        dqdc[11] = 1.8999999999999999*q_nocor;
        dqdc[12] = 3.7999999999999998*q_nocor;
        for (int k=0; k<13; k++) {
            J[14*k+0] += dqdc[k];
            J[14*k+3] += dqdc[k];
            J[14*k+4] -= dqdc[k];
        }
    }
    J[182] += dqdT; /* dwdot[H]/dT */
    J[185] += dqdT; /* dwdot[OH]/dT */
    J[186] -= dqdT; /* dwdot[H2O]/dT */

    /*reaction 7: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = 1.0000000000000002e-06 * 104000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (15286) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  15286  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[3];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[2] += dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[5] -= dqdci;                /* dwdot[O2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[3];
    J[28] -= dqdci;               /* dwdot[H]/d[O] */
    J[30] += dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[47] -= dqdci;               /* dwdot[O2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[70] -= dqdci;               /* dwdot[H]/d[O2] */
    J[72] += dqdci;               /* dwdot[O]/d[O2] */
    J[73] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] -= dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[184] += dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[187] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 8: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-06 * 3818000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7948) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  7948  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[16] -= dqdci;               /* dwdot[O]/d[H2] */
    J[17] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[28] += dqdci;               /* dwdot[H]/d[O] */
    J[29] -= dqdci;               /* dwdot[H2]/d[O] */
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 9: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-06 * 879200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (19170) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  19170  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[16] -= dqdci;               /* dwdot[O]/d[H2] */
    J[17] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[28] += dqdci;               /* dwdot[H]/d[O] */
    J[29] -= dqdci;               /* dwdot[H2]/d[O] */
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 10: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  3430  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[4];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[0] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[3] -= dqdci;                /* dwdot[OH]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[17] -= dqdci;               /* dwdot[OH]/d[H2] */
    J[18] += dqdci;               /* dwdot[H2O]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[0];
    J[56] += dqdci;               /* dwdot[H]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 11: OH + OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[3], 2.000000);
    k_f = 1.0000000000000002e-06 * 33400
                * exp(2.4199999999999999 * tc[0] - 0.50321666580471969 * (-1930) * invT);
    dlnkfdT = 2.4199999999999999 * invT + 0.50321666580471969 *  -1930  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= 2 * q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[30] += dqdci;               /* dwdot[O]/d[O] */
    J[31] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[32] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[3];
    J[44] += dqdci;               /* dwdot[O]/d[OH] */
    J[45] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[58] += dqdci;               /* dwdot[O]/d[H2O] */
    J[59] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[184] += dqdT;               /* dwdot[O]/dT */
    J[185] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 12: H2 + AR <=> H + H + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = 1.0000000000000002e-06 * 5.84e+18
                * exp(-1.1000000000000001 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.1000000000000001 * invT + 0.50321666580471969 *  104380  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000)*sc[9];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (2.000000*h_RT[0] + h_RT[9]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*2.000000*sc[0]*sc[9];
    J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[9];
    J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[1] - k_r*pow(sc[0], 2.000000);
    J[126] += 2 * dqdci;          /* dwdot[H]/d[AR] */
    J[127] -= dqdci;              /* dwdot[H2]/d[AR] */
    /* d()/dT */
    J[182] += 2 * dqdT;           /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 13: H2 + HE <=> H + H + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = 1.0000000000000002e-06 * 5.84e+18
                * exp(-1.1000000000000001 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.1000000000000001 * invT + 0.50321666580471969 *  104380  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000)*sc[10];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (2.000000*h_RT[0] + h_RT[10]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*2.000000*sc[0]*sc[10];
    J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[10];
    J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[HE] */
    dqdci =  + k_f*sc[1] - k_r*pow(sc[0], 2.000000);
    J[140] += 2 * dqdci;          /* dwdot[H]/d[HE] */
    J[141] -= dqdci;              /* dwdot[H2]/d[HE] */
    /* d()/dT */
    J[182] += 2 * dqdT;           /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 14: O + O + AR <=> O2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[2], 2.000000)*sc[9];
    k_f = 1.0000000000000002e-12 * 18860000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1788) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1788  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2] + h_RT[9]) + (h_RT[5] + h_RT[9]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2.000000*sc[2]*sc[9];
    J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[9];
    J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[AR] */
    dqdci =  + k_f*pow(sc[2], 2.000000) - k_r*sc[5];
    J[128] += -2 * dqdci;         /* dwdot[O]/d[AR] */
    J[131] += dqdci;              /* dwdot[O2]/d[AR] */
    /* d()/dT */
    J[184] += -2 * dqdT;          /* dwdot[O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 15: O + O + HE <=> O2 + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[2], 2.000000)*sc[10];
    k_f = 1.0000000000000002e-12 * 18860000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1788) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1788  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[10];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2] + h_RT[10]) + (h_RT[5] + h_RT[10]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2.000000*sc[2]*sc[10];
    J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[10];
    J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[HE] */
    dqdci =  + k_f*pow(sc[2], 2.000000) - k_r*sc[5];
    J[142] += -2 * dqdci;         /* dwdot[O]/d[HE] */
    J[145] += dqdci;              /* dwdot[O2]/d[HE] */
    /* d()/dT */
    J[184] += -2 * dqdT;          /* dwdot[O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-06 * 1.006e+26
                * exp(-2.4399999999999999 * tc[0] - 0.50321666580471969 * (120180) * invT);
    dlnkfdT = -2.4399999999999999 * invT + 0.50321666580471969 *  120180  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3]*sc[4];
    Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[0] + h_RT[3] + h_RT[4]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3]*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] -= dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0]*sc[4];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  + k_f*2.000000*sc[4] - k_r*sc[0]*sc[3];
    J[56] += dqdci;               /* dwdot[H]/d[H2O] */
    J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[186] -= dqdT;               /* dwdot[H2O]/dT */

    /*reaction 17: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = 1.0000000000000002e-06 * 2750000
                * exp(2.0899999999999999 * tc[0] - 0.50321666580471969 * (-1451) * invT);
    dlnkfdT = 2.0899999999999999 * invT + 0.50321666580471969 *  -1451  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[5] += dqdci;                /* dwdot[O2]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[14] -= dqdci;               /* dwdot[H]/d[H2] */
    J[15] += dqdci;               /* dwdot[H2]/d[H2] */
    J[19] += dqdci;               /* dwdot[O2]/d[H2] */
    J[20] -= dqdci;               /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[1];
    J[70] -= dqdci;               /* dwdot[H]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[183] += dqdT;               /* dwdot[H2]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = 1.0000000000000002e-06 * 70790000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (295) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  295  * invT2;
    /* reverse */
    phi_r = pow(sc[3], 2.000000);
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (2.000000*h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += 2 * dqdci;            /* dwdot[OH]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[3];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[45] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[87] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[185] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 19: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = 1.0000000000000002e-06 * 28500000000
                * exp(1 * tc[0] - 0.50321666580471969 * (-723.92999999999995) * invT);
    dlnkfdT = 1 * invT + 0.50321666580471969 *  -723.92999999999995  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    J[34] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[47] += dqdci;               /* dwdot[O2]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[3];
    J[72] -= dqdci;               /* dwdot[O]/d[O2] */
    J[73] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[86] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[87] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = 1.0000000000000002e-06 * 28900000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-497) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -497  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[47] += dqdci;               /* dwdot[O2]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[61] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[62] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[73] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[74] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (11982) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  11982  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[91] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[5];
    J[103] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[104] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[105] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[189] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1629.3) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1629.3  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[91] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[5];
    J[103] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[104] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[105] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[189] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 23: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = 1.0000000000000002e-06 * 24100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3970  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[56] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[98] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[101] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = 1.0000000000000002e-06 * 48200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7950) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  7950  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[6];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[14] -= dqdci;               /* dwdot[H]/d[H2] */
    J[15] += dqdci;               /* dwdot[H2]/d[H2] */
    J[20] += dqdci;               /* dwdot[HO2]/d[H2] */
    J[21] -= dqdci;               /* dwdot[H2O2]/d[H2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[1];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[98] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[99] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[183] += dqdT;               /* dwdot[H2]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 9550000
                * exp(2 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  3970  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[6];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[3] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[34] += dqdci;               /* dwdot[HO2]/d[O] */
    J[35] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[86] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[87] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[100] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[101] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 1740000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (318) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  318  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[101] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 75900000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7270) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  7270  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[101] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    amrex::Real c_R[13], dcRdT[13], e_RT[13];
    amrex::Real * eh_RT;
    if (consP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 13; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[182+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 13; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 13; ++m) {
            dehmixdc += eh_RT[m]*J[k*14+m];
        }
        J[k*14+13] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[195] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;

return;
}
#endif


/*compute the reaction Jacobian */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * consP)
{
    amrex::Real c[13];

    for (int k=0; k<13; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<13; k++) {
        J[182+k] *= 1.e-6;
        J[k*14+13] *= 1.e6;
    }

    return;
}

/* gauss-jordan solver external routine */
/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void sgjsolve(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/) {
    amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void sgjsolve_simplified(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/) {
    amrex::Abort("sgjsolve_simplified not implemented, choose a different solver ");
}

#endif
