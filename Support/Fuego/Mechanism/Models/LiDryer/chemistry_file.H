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

    extern amrex::Real fwd_A[21], fwd_beta[21], fwd_Ea[21];
    extern amrex::Real low_A[21], low_beta[21], low_Ea[21];
    extern amrex::Real rev_A[21], rev_beta[21], rev_Ea[21];
    extern amrex::Real troe_a[21],troe_Ts[21], troe_Tss[21], troe_Tsss[21];
    extern amrex::Real sri_a[21], sri_b[21], sri_c[21], sri_d[21], sri_e[21];
    extern amrex::Real activation_units[21], prefactor_units[21], phase_units[21];
    extern int is_PD[21], troe_len[21], sri_len[21], nTB[21], *TBid[21];
    extern amrex::Real *TB[21];
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
    *mm = 3;
    *kk = 9;
    *ii = 21;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/* inverse molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
    imw_new[0] = 1.0/2.015940;/*H2 */
    imw_new[1] = 1.0/31.998800;/*O2 */
    imw_new[2] = 1.0/18.015340;/*H2O */
    imw_new[3] = 1.0/1.007970;/*H */
    imw_new[4] = 1.0/15.999400;/*O */
    imw_new[5] = 1.0/17.007370;/*OH */
    imw_new[6] = 1.0/33.006770;/*HO2 */
    imw_new[7] = 1.0/34.014740;/*H2O2 */
    imw_new[8] = 1.0/28.013400;/*N2 */
}

/* molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_mw(amrex::Real *mw_new){
    mw_new[0] = 2.015940;/*H2 */
    mw_new[1] = 31.998800;/*O2 */
    mw_new[2] = 18.015340;/*H2O */
    mw_new[3] = 1.007970;/*H */
    mw_new[4] = 15.999400;/*O */
    mw_new[5] = 17.007370;/*OH */
    mw_new[6] = 33.006770;/*HO2 */
    mw_new[7] = 34.014740;/*H2O2 */
    mw_new[8] = 28.013400;/*N2 */
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
        /*species 0: H2 */
        species[0] =
            +2.29812431e+00
            +8.24944174e-04 * tc[1]
            -8.14301529e-07 * tc[2]
            -9.47543433e-11 * tc[3]
            +4.13487224e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.21293640e+00
            +1.12748635e-03 * tc[1]
            -5.75615047e-07 * tc[2]
            +1.31387723e-09 * tc[3]
            -8.76855392e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            +2.38684249e+00
            +3.47498246e-03 * tc[1]
            -6.35469633e-06 * tc[2]
            +6.96858127e-09 * tc[3]
            -2.50658847e-12 * tc[4];
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +1.94642878e+00
            -1.63816649e-03 * tc[1]
            +2.42103170e-06 * tc[2]
            -1.60284319e-09 * tc[3]
            +3.89069636e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +2.38875365e+00
            +6.56922581e-03 * tc[1]
            -1.48501258e-07 * tc[2]
            -4.62580552e-09 * tc[3]
            +2.47151475e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142337e+00
            +7.00064411e-04 * tc[1]
            -5.63382869e-08 * tc[2]
            -9.23157818e-12 * tc[3]
            +1.58275179e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.69757819e+00
            +6.13519689e-04 * tc[1]
            -1.25884199e-07 * tc[2]
            +1.77528148e-11 * tc[3]
            -1.13643531e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            +1.67214561e+00
            +3.05629289e-03 * tc[1]
            -8.73026011e-07 * tc[2]
            +1.20099639e-10 * tc[3]
            -6.39161787e-15 * tc[4];
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +1.54205966e+00
            -2.75506191e-05 * tc[1]
            -3.10280335e-09 * tc[2]
            +4.55106742e-12 * tc[3]
            -4.36805150e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.57316685e+00
            +4.33613639e-03 * tc[1]
            -1.47468882e-06 * tc[2]
            +2.34890357e-10 * tc[3]
            -1.43165356e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
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
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00
            +8.24944174e-04 * tc[1]
            -8.14301529e-07 * tc[2]
            -9.47543433e-11 * tc[3]
            +4.13487224e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00
            +1.12748635e-03 * tc[1]
            -5.75615047e-07 * tc[2]
            +1.31387723e-09 * tc[3]
            -8.76855392e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00
            +3.47498246e-03 * tc[1]
            -6.35469633e-06 * tc[2]
            +6.96858127e-09 * tc[3]
            -2.50658847e-12 * tc[4];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.94642878e+00
            -1.63816649e-03 * tc[1]
            +2.42103170e-06 * tc[2]
            -1.60284319e-09 * tc[3]
            +3.89069636e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.38875365e+00
            +6.56922581e-03 * tc[1]
            -1.48501258e-07 * tc[2]
            -4.62580552e-09 * tc[3]
            +2.47151475e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00
            +7.00064411e-04 * tc[1]
            -5.63382869e-08 * tc[2]
            -9.23157818e-12 * tc[3]
            +1.58275179e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00
            +6.13519689e-04 * tc[1]
            -1.25884199e-07 * tc[2]
            +1.77528148e-11 * tc[3]
            -1.13643531e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00
            +3.05629289e-03 * tc[1]
            -8.73026011e-07 * tc[2]
            +1.20099639e-10 * tc[3]
            -6.39161787e-15 * tc[4];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.54205966e+00
            -2.75506191e-05 * tc[1]
            -3.10280335e-09 * tc[2]
            +4.55106742e-12 * tc[3]
            -4.36805150e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57316685e+00
            +4.33613639e-03 * tc[1]
            -1.47468882e-06 * tc[2]
            +2.34890357e-10 * tc[3]
            -1.43165356e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
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
        /*species 0: H2 */
        species[0] =
            -1.012520870000000e+03 * invT
            +6.592218400000000e+00
            -3.298124310000000e+00 * tc[0]
            -4.124720870000000e-04 * tc[1]
            +1.357169215000000e-07 * tc[2]
            +7.896195275000000e-12 * tc[3]
            -2.067436120000000e-14 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.005249020000000e+03 * invT
            -2.821801190000000e+00
            -3.212936400000000e+00 * tc[0]
            -5.637431750000000e-04 * tc[1]
            +9.593584116666666e-08 * tc[2]
            -1.094897691666667e-10 * tc[3]
            +4.384276960000000e-14 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.020811330000000e+04 * invT
            +7.966096399999998e-01
            -3.386842490000000e+00 * tc[0]
            -1.737491230000000e-03 * tc[1]
            +1.059116055000000e-06 * tc[2]
            -5.807151058333333e-10 * tc[3]
            +1.253294235000000e-13 * tc[4];
        /*species 3: H */
        species[3] =
            +2.547162700000000e+04 * invT
            +2.960117608000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.914764450000000e+04 * invT
            -1.756619999999964e-02
            -2.946428780000000e+00 * tc[0]
            +8.190832450000000e-04 * tc[1]
            -4.035052833333333e-07 * tc[2]
            +1.335702658333333e-10 * tc[3]
            -1.945348180000000e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.346309130000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
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
            -1.766314650000000e+04 * invT
            -3.396609550000000e+00
            -3.388753650000000e+00 * tc[0]
            -3.284612905000000e-03 * tc[1]
            +2.475020966666667e-08 * tc[2]
            +3.854837933333333e-10 * tc[3]
            -1.235757375000000e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.020900000000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -8.350339970000000e+02 * invT
            +4.346533540000000e+00
            -2.991423370000000e+00 * tc[0]
            -3.500322055000000e-04 * tc[1]
            +9.389714483333333e-09 * tc[2]
            +7.692981816666667e-13 * tc[3]
            -7.913758950000000e-17 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.233930180000000e+03 * invT
            +5.084126000000002e-01
            -3.697578190000000e+00 * tc[0]
            -3.067598445000000e-04 * tc[1]
            +2.098069983333333e-08 * tc[2]
            -1.479401233333333e-12 * tc[3]
            +5.682176550000000e-17 * tc[4];
        /*species 2: H2O */
        species[2] =
            -2.989920900000000e+04 * invT
            -4.190671200000001e+00
            -2.672145610000000e+00 * tc[0]
            -1.528146445000000e-03 * tc[1]
            +1.455043351666667e-07 * tc[2]
            -1.000830325000000e-11 * tc[3]
            +3.195808935000000e-16 * tc[4];
        /*species 3: H */
        species[3] =
            +2.547162700000000e+04 * invT
            +2.960117638000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.923080270000000e+04 * invT
            -2.378248450000000e+00
            -2.542059660000000e+00 * tc[0]
            +1.377530955000000e-05 * tc[1]
            +5.171338916666667e-10 * tc[2]
            -3.792556183333333e-13 * tc[3]
            +2.184025750000000e-17 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.683628750000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
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
            -1.800696090000000e+04 * invT
            +4.072029891000000e+00
            -4.573166850000000e+00 * tc[0]
            -2.168068195000000e-03 * tc[1]
            +2.457814700000000e-07 * tc[2]
            -1.957419641666667e-11 * tc[3]
            +7.158267800000000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
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
        /*species 0: H2 */
        species[0] =
            -1.01252087e+03 * invT
            +5.59221840e+00
            -3.29812431e+00 * tc[0]
            -4.12472087e-04 * tc[1]
            +1.35716922e-07 * tc[2]
            +7.89619527e-12 * tc[3]
            -2.06743612e-14 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.00524902e+03 * invT
            -3.82180119e+00
            -3.21293640e+00 * tc[0]
            -5.63743175e-04 * tc[1]
            +9.59358412e-08 * tc[2]
            -1.09489769e-10 * tc[3]
            +4.38427696e-14 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.02081133e+04 * invT
            -2.03390360e-01
            -3.38684249e+00 * tc[0]
            -1.73749123e-03 * tc[1]
            +1.05911606e-06 * tc[2]
            -5.80715106e-10 * tc[3]
            +1.25329424e-13 * tc[4];
        /*species 3: H */
        species[3] =
            +2.54716270e+04 * invT
            +1.96011761e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.91476445e+04 * invT
            -1.01756620e+00
            -2.94642878e+00 * tc[0]
            +8.19083245e-04 * tc[1]
            -4.03505283e-07 * tc[2]
            +1.33570266e-10 * tc[3]
            -1.94534818e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.34630913e+03 * invT
            +3.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
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
            -1.76631465e+04 * invT
            -4.39660955e+00
            -3.38875365e+00 * tc[0]
            -3.28461290e-03 * tc[1]
            +2.47502097e-08 * tc[2]
            +3.85483793e-10 * tc[3]
            -1.23575738e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.02090000e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -8.35033997e+02 * invT
            +3.34653354e+00
            -2.99142337e+00 * tc[0]
            -3.50032206e-04 * tc[1]
            +9.38971448e-09 * tc[2]
            +7.69298182e-13 * tc[3]
            -7.91375895e-17 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.23393018e+03 * invT
            -4.91587400e-01
            -3.69757819e+00 * tc[0]
            -3.06759845e-04 * tc[1]
            +2.09806998e-08 * tc[2]
            -1.47940123e-12 * tc[3]
            +5.68217655e-17 * tc[4];
        /*species 2: H2O */
        species[2] =
            -2.98992090e+04 * invT
            -5.19067120e+00
            -2.67214561e+00 * tc[0]
            -1.52814644e-03 * tc[1]
            +1.45504335e-07 * tc[2]
            -1.00083033e-11 * tc[3]
            +3.19580894e-16 * tc[4];
        /*species 3: H */
        species[3] =
            +2.54716270e+04 * invT
            +1.96011764e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.92308027e+04 * invT
            -3.37824845e+00
            -2.54205966e+00 * tc[0]
            +1.37753096e-05 * tc[1]
            +5.17133892e-10 * tc[2]
            -3.79255618e-13 * tc[3]
            +2.18402575e-17 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.68362875e+03 * invT
            -3.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
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
            -1.80069609e+04 * invT
            +3.07202989e+00
            -4.57316685e+00 * tc[0]
            -2.16806820e-03 * tc[1]
            +2.45781470e-07 * tc[2]
            -1.95741964e-11 * tc[3]
            +7.15826780e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
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
        /*species 0: H2 */
        species[0] =
            +2.29812431e+00
            +4.12472087e-04 * tc[1]
            -2.71433843e-07 * tc[2]
            -2.36885858e-11 * tc[3]
            +8.26974448e-14 * tc[4]
            -1.01252087e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +2.21293640e+00
            +5.63743175e-04 * tc[1]
            -1.91871682e-07 * tc[2]
            +3.28469308e-10 * tc[3]
            -1.75371078e-13 * tc[4]
            -1.00524902e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +2.38684249e+00
            +1.73749123e-03 * tc[1]
            -2.11823211e-06 * tc[2]
            +1.74214532e-09 * tc[3]
            -5.01317694e-13 * tc[4]
            -3.02081133e+04 * invT;
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +1.94642878e+00
            -8.19083245e-04 * tc[1]
            +8.07010567e-07 * tc[2]
            -4.00710797e-10 * tc[3]
            +7.78139272e-14 * tc[4]
            +2.91476445e+04 * invT;
        /*species 5: OH */
        species[5] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
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
            +2.38875365e+00
            +3.28461290e-03 * tc[1]
            -4.95004193e-08 * tc[2]
            -1.15645138e-09 * tc[3]
            +4.94302950e-13 * tc[4]
            -1.76631465e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142337e+00
            +3.50032206e-04 * tc[1]
            -1.87794290e-08 * tc[2]
            -2.30789455e-12 * tc[3]
            +3.16550358e-16 * tc[4]
            -8.35033997e+02 * invT;
        /*species 1: O2 */
        species[1] =
            +2.69757819e+00
            +3.06759845e-04 * tc[1]
            -4.19613997e-08 * tc[2]
            +4.43820370e-12 * tc[3]
            -2.27287062e-16 * tc[4]
            -1.23393018e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +1.67214561e+00
            +1.52814644e-03 * tc[1]
            -2.91008670e-07 * tc[2]
            +3.00249098e-11 * tc[3]
            -1.27832357e-15 * tc[4]
            -2.98992090e+04 * invT;
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +1.54205966e+00
            -1.37753096e-05 * tc[1]
            -1.03426778e-09 * tc[2]
            +1.13776685e-12 * tc[3]
            -8.73610300e-17 * tc[4]
            +2.92308027e+04 * invT;
        /*species 5: OH */
        species[5] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
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
            +3.57316685e+00
            +2.16806820e-03 * tc[1]
            -4.91562940e-07 * tc[2]
            +5.87225893e-11 * tc[3]
            -2.86330712e-15 * tc[4]
            -1.80069609e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
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
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00
            +4.12472087e-04 * tc[1]
            -2.71433843e-07 * tc[2]
            -2.36885858e-11 * tc[3]
            +8.26974448e-14 * tc[4]
            -1.01252087e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00
            +5.63743175e-04 * tc[1]
            -1.91871682e-07 * tc[2]
            +3.28469308e-10 * tc[3]
            -1.75371078e-13 * tc[4]
            -1.00524902e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00
            +1.73749123e-03 * tc[1]
            -2.11823211e-06 * tc[2]
            +1.74214532e-09 * tc[3]
            -5.01317694e-13 * tc[4]
            -3.02081133e+04 * invT;
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +2.94642878e+00
            -8.19083245e-04 * tc[1]
            +8.07010567e-07 * tc[2]
            -4.00710797e-10 * tc[3]
            +7.78139272e-14 * tc[4]
            +2.91476445e+04 * invT;
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
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
            +3.38875365e+00
            +3.28461290e-03 * tc[1]
            -4.95004193e-08 * tc[2]
            -1.15645138e-09 * tc[3]
            +4.94302950e-13 * tc[4]
            -1.76631465e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00
            +3.50032206e-04 * tc[1]
            -1.87794290e-08 * tc[2]
            -2.30789455e-12 * tc[3]
            +3.16550358e-16 * tc[4]
            -8.35033997e+02 * invT;
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00
            +3.06759845e-04 * tc[1]
            -4.19613997e-08 * tc[2]
            +4.43820370e-12 * tc[3]
            -2.27287062e-16 * tc[4]
            -1.23393018e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00
            +1.52814644e-03 * tc[1]
            -2.91008670e-07 * tc[2]
            +3.00249098e-11 * tc[3]
            -1.27832357e-15 * tc[4]
            -2.98992090e+04 * invT;
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +2.54205966e+00
            -1.37753096e-05 * tc[1]
            -1.03426778e-09 * tc[2]
            +1.13776685e-12 * tc[3]
            -8.73610300e-17 * tc[4]
            +2.92308027e+04 * invT;
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
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
            +4.57316685e+00
            +2.16806820e-03 * tc[1]
            -4.91562940e-07 * tc[2]
            +5.87225893e-11 * tc[3]
            -2.86330712e-15 * tc[4]
            -1.80069609e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
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
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00 * tc[0]
            +8.24944174e-04 * tc[1]
            -4.07150765e-07 * tc[2]
            -3.15847811e-11 * tc[3]
            +1.03371806e-13 * tc[4]
            -3.29409409e+00 ;
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00 * tc[0]
            +1.12748635e-03 * tc[1]
            -2.87807523e-07 * tc[2]
            +4.37959077e-10 * tc[3]
            -2.19213848e-13 * tc[4]
            +6.03473759e+00 ;
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00 * tc[0]
            +3.47498246e-03 * tc[1]
            -3.17734817e-06 * tc[2]
            +2.32286042e-09 * tc[3]
            -6.26647117e-13 * tc[4]
            +2.59023285e+00 ;
        /*species 3: H */
        species[3] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117608e-01 ;
        /*species 4: O */
        species[4] =
            +2.94642878e+00 * tc[0]
            -1.63816649e-03 * tc[1]
            +1.21051585e-06 * tc[2]
            -5.34281063e-10 * tc[3]
            +9.72674090e-14 * tc[4]
            +2.96399498e+00 ;
        /*species 5: OH */
        species[5] =
            +4.12530561e+00 * tc[0]
            -3.22544939e-03 * tc[1]
            +3.26382346e-06 * tc[2]
            -1.93284548e-09 * tc[3]
            +5.15593447e-13 * tc[4]
            -6.90432960e-01 ;
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
            +3.38875365e+00 * tc[0]
            +6.56922581e-03 * tc[1]
            -7.42506290e-08 * tc[2]
            -1.54193517e-09 * tc[3]
            +6.17878688e-13 * tc[4]
            +6.78536320e+00 ;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00 * tc[0]
            +1.40824000e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213750e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00 * tc[0]
            +7.00064411e-04 * tc[1]
            -2.81691434e-08 * tc[2]
            -3.07719273e-12 * tc[3]
            +3.95687948e-16 * tc[4]
            -1.35511017e+00 ;
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00 * tc[0]
            +6.13519689e-04 * tc[1]
            -6.29420995e-08 * tc[2]
            +5.91760493e-12 * tc[3]
            -2.84108828e-16 * tc[4]
            +3.18916559e+00 ;
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00 * tc[0]
            +3.05629289e-03 * tc[1]
            -4.36513005e-07 * tc[2]
            +4.00332130e-11 * tc[3]
            -1.59790447e-15 * tc[4]
            +6.86281681e+00 ;
        /*species 3: H */
        species[3] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117638e-01 ;
        /*species 4: O */
        species[4] =
            +2.54205966e+00 * tc[0]
            -2.75506191e-05 * tc[1]
            -1.55140167e-09 * tc[2]
            +1.51702247e-12 * tc[3]
            -1.09201287e-16 * tc[4]
            +4.92030811e+00 ;
        /*species 5: OH */
        species[5] =
            +2.86472886e+00 * tc[0]
            +1.05650448e-03 * tc[1]
            -1.29541379e-07 * tc[2]
            +1.01739558e-11 * tc[3]
            -3.32989690e-16 * tc[4]
            +5.70164073e+00 ;
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
            +4.57316685e+00 * tc[0]
            +4.33613639e-03 * tc[1]
            -7.37344410e-07 * tc[2]
            +7.82967857e-11 * tc[3]
            -3.57913390e-15 * tc[4]
            +5.01136959e-01 ;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00 * tc[0]
            +1.48797700e-03 * tc[1]
            -2.84238050e-07 * tc[2]
            +3.36568000e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
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
    amrex::Real cpor[9]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
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
    amrex::Real cpor[9], tresult[9]; /* temporary storage */
    amrex::Real imw[9];

    get_imw(imw);

    cp_R(cpor, tc);

    for (int i = 0; i < 9; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 9; i++)
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
    amrex::Real cvor[9]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
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
    amrex::Real cvor[9]; /* temporary storage */
    amrex::Real imw[9];

    get_imw(imw);

    cv_R(cvor, tc);

    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*O2 */
    result += cvor[2]*y[2]*imw[2]; /*H2O */
    result += cvor[3]*y[3]*imw[3]; /*H */
    result += cvor[4]*y[4]*imw[4]; /*O */
    result += cvor[5]*y[5]*imw[5]; /*OH */
    result += cvor[6]*y[6]*imw[6]; /*HO2 */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*N2 */

    *cvbs = result * 8.31446261815324e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBML(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  hbml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real hml[9]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
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
    amrex::Real hml[9], tmp[9]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[9];

    get_imw(imw);

    speciesEnthalpy(hml, tc);

    int id;
    for (id = 0; id < 9; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 9; ++id) {
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
    amrex::Real uml[9]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
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
    amrex::Real ums[9]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[9];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*O2 */
    result += y[2]*ums[2]*imw[2]; /*H2O */
    result += y[3]*ums[3]*imw[3]; /*H */
    result += y[4]*ums[4]*imw[4]; /*O */
    result += y[5]*ums[5]*imw[5]; /*OH */
    result += y[6]*ums[6]*imw[6]; /*HO2 */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*N2 */

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
    amrex::Real sor[9]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 9; ++id) {
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
    amrex::Real sor[9]; /* temporary storage */
    amrex::Real x[9]; /* need a ytx conversion */
    amrex::Real YOW = 0; /*See Eq 4, 6 in CK Manual */
    amrex::Real imw[9];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(18.015340*YOW); 
    x[3] = y[3]/(1.007970*YOW); 
    x[4] = y[4]/(15.999400*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
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
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *P = *rho * 8.31446261815324e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPY(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW = 0;/* for computing mean MW */
    amrex::Real imw[9];

    get_imw(imw);

    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
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
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*1.007970; /*H */
    W += c[4]*15.999400; /*O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31446261815324e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOX(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  rho)
{
    amrex::Real XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *rho = *P * XW / (8.31446261815324e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOY(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  rho)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[9];
    amrex::Real imw[9];

    get_imw(imw);

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
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
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*1.007970; /*H */
    W += c[4]*15.999400; /*O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
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
    amrex::Real tmp[9];
    amrex::Real imw[9];

    get_imw(imw);

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
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
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
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
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*1.007970; /*H */
    W += c[4]*15.999400; /*O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
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
    amrex::Real tmp[9];
    amrex::Real imw[9];

    get_imw(imw);

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    amrex::Real YOWINV = 1.0/YOW;

    for (int i = 0; i < 9; i++)
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
    amrex::Real imw[9];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 9; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446261815324e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 9; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCR(amrex::Real *  rho, amrex::Real * /*T*/, amrex::Real * y,  amrex::Real * c)
{
    amrex::Real imw[9];

    get_imw(imw);

    for (int i = 0; i < 9; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTY(amrex::Real *  x,  amrex::Real *  y)
{
    amrex::Real XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Now compute conversion */
    amrex::Real XWinv = 1.0/XW;
    y[0] = x[0]*2.015940*XWinv; 
    y[1] = x[1]*31.998800*XWinv; 
    y[2] = x[2]*18.015340*XWinv; 
    y[3] = x[3]*1.007970*XWinv; 
    y[4] = x[4]*15.999400*XWinv; 
    y[5] = x[5]*17.007370*XWinv; 
    y[6] = x[6]*33.006770*XWinv; 
    y[7] = x[7]*34.014740*XWinv; 
    y[8] = x[8]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  c)
{
    int id; /*loop counter */
    amrex::Real PORT = (*P)/(8.31446261815324e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
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
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    amrex::Real sumCinv = 1.0/sumC;
    for (id = 0; id < 9; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTY(amrex::Real *  c, amrex::Real *  y)
{
    amrex::Real CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*2.015940; /*H2 */
    CW += c[1]*31.998800; /*O2 */
    CW += c[2]*18.015340; /*H2O */
    CW += c[3]*1.007970; /*H */
    CW += c[4]*15.999400; /*O */
    CW += c[5]*17.007370; /*OH */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*28.013400; /*N2 */
    /*Now compute conversion */
    amrex::Real CWinv = 1.0/CW;
    y[0] = c[0]*2.015940*CWinv; 
    y[1] = c[1]*31.998800*CWinv; 
    y[2] = c[2]*18.015340*CWinv; 
    y[3] = c[3]*1.007970*CWinv; 
    y[4] = c[4]*15.999400*CWinv; 
    y[5] = c[5]*17.007370*CWinv; 
    y[6] = c[6]*33.006770*CWinv; 
    y[7] = c[7]*34.014740*CWinv; 
    y[8] = c[8]*28.013400*CWinv; 

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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    cvms[0] *= 4.124360158612479e+07; /*H2 */
    cvms[1] *= 2.598367006935648e+06; /*O2 */
    cvms[2] *= 4.615212712140454e+06; /*H2O */
    cvms[3] *= 8.248720317224957e+07; /*H */
    cvms[4] *= 5.196734013871295e+06; /*O */
    cvms[5] *= 4.888740950630956e+06; /*OH */
    cvms[6] *= 2.519017346487778e+06; /*HO2 */
    cvms[7] *= 2.444370475315478e+06; /*H2O2 */
    cvms[8] *= 2.968030520448514e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPMS(amrex::Real *  T,  amrex::Real *  cpms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124360158612479e+07; /*H2 */
    cpms[1] *= 2.598367006935648e+06; /*O2 */
    cpms[2] *= 4.615212712140454e+06; /*H2O */
    cpms[3] *= 8.248720317224957e+07; /*H */
    cpms[4] *= 5.196734013871295e+06; /*O */
    cpms[5] *= 4.888740950630956e+06; /*OH */
    cpms[6] *= 2.519017346487778e+06; /*HO2 */
    cpms[7] *= 2.444370475315478e+06; /*H2O2 */
    cpms[8] *= 2.968030520448514e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUMS(amrex::Real *  T,  amrex::Real *  ums)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[9];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    for (int i = 0; i < 9; i++)
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
    amrex::Real imw[9];

    get_imw(imw);

    speciesEnthalpy(hms, tc);

    for (int i = 0; i < 9; i++)
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
    sms[0] *= 4.124360158612479e+07; /*H2 */
    sms[1] *= 2.598367006935648e+06; /*O2 */
    sms[2] *= 4.615212712140454e+06; /*H2O */
    sms[3] *= 8.248720317224957e+07; /*H */
    sms[4] *= 5.196734013871295e+06; /*O */
    sms[5] *= 4.888740950630956e+06; /*OH */
    sms[6] *= 2.519017346487778e+06; /*HO2 */
    sms[7] *= 2.444370475315478e+06; /*H2O2 */
    sms[8] *= 2.968030520448514e+06; /*N2 */
}

#ifdef AMREX_USE_GPU

/* GPU version of productionRate: no more use of thermo namespace vectors  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qfqr(amrex::Real *  qf, amrex::Real * qr, amrex::Real * sc, amrex::Real * qss_sc, amrex::Real * tc, amrex::Real invT)
{

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[1]*sc[3];
    qr[0] = sc[6];

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    qf[1] = sc[7];
    qr[1] = sc[5]*sc[5];

    /*reaction 5: H2 + M <=> H + H + M */
    qf[2] = sc[0];
    qr[2] = sc[3]*sc[3];

    /*reaction 6: O + O + M <=> O2 + M */
    qf[3] = sc[4]*sc[4];
    qr[3] = sc[1];

    /*reaction 7: O + H + M <=> OH + M */
    qf[4] = sc[3]*sc[4];
    qr[4] = sc[5];

    /*reaction 8: H + OH + M <=> H2O + M */
    qf[5] = sc[3]*sc[5];
    qr[5] = sc[2];

    /*reaction 1: H + O2 <=> O + OH */
    qf[6] = sc[1]*sc[3];
    qr[6] = sc[4]*sc[5];

    /*reaction 2: O + H2 <=> H + OH */
    qf[7] = sc[0]*sc[4];
    qr[7] = sc[3]*sc[5];

    /*reaction 3: H2 + OH <=> H2O + H */
    qf[8] = sc[0]*sc[5];
    qr[8] = sc[2]*sc[3];

    /*reaction 4: O + H2O <=> OH + OH */
    qf[9] = sc[2]*sc[4];
    qr[9] = sc[5]*sc[5];

    /*reaction 10: HO2 + H <=> H2 + O2 */
    qf[10] = sc[3]*sc[6];
    qr[10] = sc[0]*sc[1];

    /*reaction 11: HO2 + H <=> OH + OH */
    qf[11] = sc[3]*sc[6];
    qr[11] = sc[5]*sc[5];

    /*reaction 12: HO2 + O <=> O2 + OH */
    qf[12] = sc[4]*sc[6];
    qr[12] = sc[1]*sc[5];

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    qf[13] = sc[5]*sc[6];
    qr[13] = sc[1]*sc[2];

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    qf[14] = sc[6]*sc[6];
    qr[14] = sc[1]*sc[7];

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    qf[15] = sc[6]*sc[6];
    qr[15] = sc[1]*sc[7];

    /*reaction 17: H2O2 + H <=> H2O + OH */
    qf[16] = sc[3]*sc[7];
    qr[16] = sc[2]*sc[5];

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    qf[17] = sc[3]*sc[7];
    qr[17] = sc[0]*sc[6];

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    qf[18] = sc[4]*sc[7];
    qr[18] = sc[5]*sc[6];

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    qf[19] = sc[5]*sc[7];
    qr[19] = sc[2]*sc[6];

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    qf[20] = sc[5]*sc[7];
    qr[20] = sc[2]*sc[6];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 9; ++i) {
        mixture += sc[i];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[9];
    gibbs(g_RT, tc);

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 * invT;
    amrex::Real refCinv = 1 / refC;

    /* Evaluate the kfs */
    amrex::Real k_f, k_r, Corr;
    amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

    // (0):  H + O2 <=> O + OH
    k_f = 1.0000000000000002e-06 * 3547000000000000 
               * exp(-0.40600000000000003 * tc[0] - 0.50321666580471969 * (16599) * invT);
    Corr  = 1.0;
    qf[6] *= Corr * k_f;
    qr[6] *= Corr * k_f / exp(g_RT[1] + g_RT[3] - g_RT[4] - g_RT[5]);
    // (1):  O + H2 <=> H + OH
    k_f = 1.0000000000000002e-06 * 50800 
               * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6290) * invT);
    Corr  = 1.0;
    qf[7] *= Corr * k_f;
    qr[7] *= Corr * k_f / exp(g_RT[0] - g_RT[3] + g_RT[4] - g_RT[5]);
    // (2):  H2 + OH <=> H2O + H
    k_f = 1.0000000000000002e-06 * 216000000 
               * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    Corr  = 1.0;
    qf[8] *= Corr * k_f;
    qr[8] *= Corr * k_f / exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    // (3):  O + H2O <=> OH + OH
    k_f = 1.0000000000000002e-06 * 2970000 
               * exp(2.02 * tc[0] - 0.50321666580471969 * (13400) * invT);
    Corr  = 1.0;
    qf[9] *= Corr * k_f;
    qr[9] *= Corr * k_f / exp(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[5]);
    // (4):  H2 + M <=> H + H + M
    k_f = 1.0000000000000002e-06 * 4.577e+19 
               * exp(-1.3999999999999999 * tc[0] - 0.50321666580471969 * (104380) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    qf[2] *= Corr * k_f;
    qr[2] *= Corr * k_f / (exp(g_RT[0] - g_RT[3] - g_RT[3]) * refC);
    // (5):  O + O + M <=> O2 + M
    k_f = 1.0000000000000002e-12 * 6165000000000000 
               * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    qf[3] *= Corr * k_f;
    qr[3] *= Corr * k_f / (exp(-g_RT[1] + g_RT[4] + g_RT[4]) * refCinv);
    // (6):  O + H + M <=> OH + M
    k_f = 1.0000000000000002e-12 * 4.714e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    qf[4] *= Corr * k_f;
    qr[4] *= Corr * k_f / (exp(g_RT[3] + g_RT[4] - g_RT[5]) * refCinv);
    // (7):  H + OH + M <=> H2O + M
    k_f = 1.0000000000000002e-12 * 3.8000000000000004e+22 
               * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    qf[5] *= Corr * k_f;
    qr[5] *= Corr * k_f / (exp(-g_RT[2] + g_RT[3] + g_RT[5]) * refCinv);
    // (8):  H + O2 (+M) <=> HO2 (+M)
    k_f = 1.0000000000000002e-06 * 1475000000000 
               * exp(0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 11 - 1)*sc[2] + ( 0.78000000000000003 - 1)*sc[1];
    redP = Corr / k_f * 1e-12 * 6.366e+20 
               * exp(-1.72  * tc[0] - 0.50321666580471969  * (524.79999999999995) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.80000000000000004)*exp(-tc[1] / 1.0000000000000001e-30) 
        + 0.80000000000000004 * exp(-tc[1]/1e+30)  
        + 0.);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[0] *= Corr * k_f;
    qr[0] *= Corr * k_f / (exp(g_RT[1] + g_RT[3] - g_RT[6]) * refCinv);
    // (9):  HO2 + H <=> H2 + O2
    k_f = 1.0000000000000002e-06 * 16600000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (823) * invT);
    Corr  = 1.0;
    qf[10] *= Corr * k_f;
    qr[10] *= Corr * k_f / exp(-g_RT[0] - g_RT[1] + g_RT[3] + g_RT[6]);
    // (10):  HO2 + H <=> OH + OH
    k_f = 1.0000000000000002e-06 * 70790000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (295) * invT);
    Corr  = 1.0;
    qf[11] *= Corr * k_f;
    qr[11] *= Corr * k_f / exp(g_RT[3] - g_RT[5] - g_RT[5] + g_RT[6]);
    // (11):  HO2 + O <=> O2 + OH
    k_f = 1.0000000000000002e-06 * 32500000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[12] *= Corr * k_f;
    qr[12] *= Corr * k_f / exp(-g_RT[1] + g_RT[4] - g_RT[5] + g_RT[6]);
    // (12):  HO2 + OH <=> H2O + O2
    k_f = 1.0000000000000002e-06 * 28900000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-497) * invT);
    Corr  = 1.0;
    qf[13] *= Corr * k_f;
    qr[13] *= Corr * k_f / exp(-g_RT[1] - g_RT[2] + g_RT[5] + g_RT[6]);
    // (13):  HO2 + HO2 <=> H2O2 + O2
    k_f = 1.0000000000000002e-06 * 420000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (11982) * invT);
    Corr  = 1.0;
    qf[14] *= Corr * k_f;
    qr[14] *= Corr * k_f / exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    // (14):  HO2 + HO2 <=> H2O2 + O2
    k_f = 1.0000000000000002e-06 * 130000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1629.3) * invT);
    Corr  = 1.0;
    qf[15] *= Corr * k_f;
    qr[15] *= Corr * k_f / exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    // (15):  H2O2 (+M) <=> OH + OH (+M)
    k_f = 1 * 295100000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (48430) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    redP = Corr / k_f * 1e-6 * 1.202e+17 
               * exp(0  * tc[0] - 0.50321666580471969  * (45500) *invT);
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
    qf[1] *= Corr * k_f;
    qr[1] *= Corr * k_f / (exp(-g_RT[5] - g_RT[5] + g_RT[7]) * refC);
    // (16):  H2O2 + H <=> H2O + OH
    k_f = 1.0000000000000002e-06 * 24100000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (3970) * invT);
    Corr  = 1.0;
    qf[16] *= Corr * k_f;
    qr[16] *= Corr * k_f / exp(-g_RT[2] + g_RT[3] - g_RT[5] + g_RT[7]);
    // (17):  H2O2 + H <=> HO2 + H2
    k_f = 1.0000000000000002e-06 * 48200000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (7950) * invT);
    Corr  = 1.0;
    qf[17] *= Corr * k_f;
    qr[17] *= Corr * k_f / exp(-g_RT[0] + g_RT[3] - g_RT[6] + g_RT[7]);
    // (18):  H2O2 + O <=> OH + HO2
    k_f = 1.0000000000000002e-06 * 9550000 
               * exp(2 * tc[0] - 0.50321666580471969 * (3970) * invT);
    Corr  = 1.0;
    qf[18] *= Corr * k_f;
    qr[18] *= Corr * k_f / exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    // (19):  H2O2 + OH <=> HO2 + H2O
    k_f = 1.0000000000000002e-06 * 1000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[19] *= Corr * k_f;
    qr[19] *= Corr * k_f / exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    // (20):  H2O2 + OH <=> HO2 + H2O
    k_f = 1.0000000000000002e-06 * 580000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (9557) * invT);
    Corr  = 1.0;
    qf[20] *= Corr * k_f;
    qr[20] *= Corr * k_f / exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);


    return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real * wdot, amrex::Real * sc, amrex::Real T)
{
    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];

    amrex::Real qdot, q_f[21], q_r[21];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 9; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[5] += qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[2]-q_r[2];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[3] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[4] -= qdot;

    qdot = q_f[4]-q_r[4];
    wdot[3] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[9]-q_r[9];
    wdot[2] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[0] += qdot;
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[6] -= qdot;

    qdot = q_f[11]-q_r[11];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[12]-q_r[12];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[13]-q_r[13];
    wdot[1] += qdot;
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[6] -= qdot;

    qdot = q_f[14]-q_r[14];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[0] += qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[2] += qdot;
    wdot[5] -= qdot;
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
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[9]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[9];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
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

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[9]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[9]; /*temporary storage */
    amrex::Real imw[9];

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

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[9]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
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
        /*species 0: H2 */
        species[0] =
            +8.24944174e-04
            -1.62860306e-06 * tc[1]
            -2.84263030e-10 * tc[2]
            +1.65394890e-12 * tc[3];
        /*species 1: O2 */
        species[1] =
            +1.12748635e-03
            -1.15123009e-06 * tc[1]
            +3.94163169e-09 * tc[2]
            -3.50742157e-12 * tc[3];
        /*species 2: H2O */
        species[2] =
            +3.47498246e-03
            -1.27093927e-05 * tc[1]
            +2.09057438e-08 * tc[2]
            -1.00263539e-11 * tc[3];
        /*species 3: H */
        species[3] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 4: O */
        species[4] =
            -1.63816649e-03
            +4.84206340e-06 * tc[1]
            -4.80852957e-09 * tc[2]
            +1.55627854e-12 * tc[3];
        /*species 5: OH */
        species[5] =
            -3.22544939e-03
            +1.30552938e-05 * tc[1]
            -1.73956093e-08 * tc[2]
            +8.24949516e-12 * tc[3];
        /*species 6: HO2 */
        species[6] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +6.56922581e-03
            -2.97002516e-07 * tc[1]
            -1.38774166e-08 * tc[2]
            +9.88605900e-12 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.40824000e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77942000e-12 * tc[3];
    } else {
        /*species 0: H2 */
        species[0] =
            +7.00064411e-04
            -1.12676574e-07 * tc[1]
            -2.76947345e-11 * tc[2]
            +6.33100716e-15 * tc[3];
        /*species 1: O2 */
        species[1] =
            +6.13519689e-04
            -2.51768398e-07 * tc[1]
            +5.32584444e-11 * tc[2]
            -4.54574124e-15 * tc[3];
        /*species 2: H2O */
        species[2] =
            +3.05629289e-03
            -1.74605202e-06 * tc[1]
            +3.60298917e-10 * tc[2]
            -2.55664715e-14 * tc[3];
        /*species 3: H */
        species[3] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 4: O */
        species[4] =
            -2.75506191e-05
            -6.20560670e-09 * tc[1]
            +1.36532023e-11 * tc[2]
            -1.74722060e-15 * tc[3];
        /*species 5: OH */
        species[5] =
            +1.05650448e-03
            -5.18165516e-07 * tc[1]
            +9.15656022e-11 * tc[2]
            -5.32783504e-15 * tc[3];
        /*species 6: HO2 */
        species[6] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.33613639e-03
            -2.94937764e-06 * tc[1]
            +7.04671071e-10 * tc[2]
            -5.72661424e-14 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.48797700e-03
            -1.13695220e-06 * tc[1]
            +3.02911200e-10 * tc[2]
            -2.70134040e-14 * tc[3];
    }
    return;
}

/*compute an approx to the reaction Jacobian */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void aJacobian_precond(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int HP)
{
    for (int i=0; i<100; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[9];
    for (int k=0; k<9; k++) {
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
    for (int k = 0; k < 9; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[9];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[9];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[9];
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
    alpha = mixture + ( 2 - 1)*sc[0] + ( 11 - 1)*sc[2] + ( 0.78000000000000003 - 1)*sc[1];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 1475000000000
                * exp(0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0.59999999999999998 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* pressure-fall-off */
    k_0 = 6.366e+20 * exp(-1.72 * tc[0] - 0.50321666580471969 * (524.79999999999995) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -1.72 * invT + 0.50321666580471969 * (524.79999999999995) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.80000000000000004)*exp(-T/1.0000000000000001e-30);
    Fcent2 = 0.80000000000000004 * exp(-T/1e+30);
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
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = 0.78000000000000003*dcdc_fac + k_f*sc[3];
    dqdc[2] = 11*dcdc_fac;
    dqdc[3] = dcdc_fac + k_f*sc[1];
    dqdc[4] = dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = dcdc_fac - k_r;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    for (int k=0; k<9; k++) {
        J[10*k+1] -= dqdc[k];
        J[10*k+3] -= dqdc[k];
        J[10*k+6] += dqdc[k];
    }
    J[91] -= dqdT; /* dwdot[O2]/dT */
    J[93] -= dqdT; /* dwdot[H]/dT */
    J[96] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    /* forward */
    phi_f = sc[7];
    k_f = 1 * 295100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (48430) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (48430)  * invT2;
    /* pressure-fall-off */
    k_0 = 1.202e+17 * exp(0 * tc[0] - 0.50321666580471969 * (45500) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = 0 * invT + 0.50321666580471969 * (45500) * invT2;
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
    phi_r = pow(sc[5], 2.000000);
    Kc = refC * exp(-g_RT[5] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2.000000*h_RT[5]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[5] += 2 * q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2.5*dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = 12*dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = dcdc_fac - k_r*2.000000*sc[5];
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac + k_f;
    dqdc[8] = dcdc_fac;
    for (int k=0; k<9; k++) {
        J[10*k+5] += 2 * dqdc[k];
        J[10*k+7] -= dqdc[k];
    }
    J[95] += 2 * dqdT; /* dwdot[OH]/dT */
    J[97] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    /* forward */
    phi_f = sc[0];
    k_f = 1.0000000000000002e-06 * 4.577e+19
                * exp(-1.3999999999999999 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.3999999999999999 * invT + 0.50321666580471969 *  (104380)  * invT2;
    /* reverse */
    phi_r = pow(sc[3], 2.000000);
    Kc = refC * exp(g_RT[0] - g_RT[3] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0]) + (2.000000*h_RT[3]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[3] += 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2.5*q_nocor + k_f;
    dqdc[1] = q_nocor;
    dqdc[2] = 12*q_nocor;
    dqdc[3] = q_nocor - k_r*2.000000*sc[3];
    dqdc[4] = q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+0] -= dqdc[k];
        J[10*k+3] += 2 * dqdc[k];
    }
    J[90] -= dqdT; /* dwdot[H2]/dT */
    J[93] += 2 * dqdT; /* dwdot[H]/dT */

    /*reaction 4: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-12 * 6165000000000000
                * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.5 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[1];
    Kc = refCinv * exp(-g_RT[1] + g_RT[4] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[1]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[4] -= 2 * q; /* O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2.5*q_nocor;
    dqdc[1] = q_nocor - k_r;
    dqdc[2] = 12*q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor + k_f*2.000000*sc[4];
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+1] += dqdc[k];
        J[10*k+4] += -2 * dqdc[k];
    }
    J[91] += dqdT; /* dwdot[O2]/dT */
    J[94] += -2 * dqdT; /* dwdot[O]/dT */

    /*reaction 5: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    /* forward */
    phi_f = sc[3]*sc[4];
    k_f = 1.0000000000000002e-12 * 4.714e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[4]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* H */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2.5*q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = 12*q_nocor;
    dqdc[3] = q_nocor + k_f*sc[4];
    dqdc[4] = q_nocor + k_f*sc[3];
    dqdc[5] = q_nocor - k_r;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+3] -= dqdc[k];
        J[10*k+4] -= dqdc[k];
        J[10*k+5] += dqdc[k];
    }
    J[93] -= dqdT; /* dwdot[H]/dT */
    J[94] -= dqdT; /* dwdot[O]/dT */
    J[95] += dqdT; /* dwdot[OH]/dT */

    /*reaction 6: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    /* forward */
    phi_f = sc[3]*sc[5];
    k_f = 1.0000000000000002e-12 * 3.8000000000000004e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[2];
    Kc = refCinv * exp(-g_RT[2] + g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[2]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[3] -= q; /* H */
    wdot[5] -= q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2.5*q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = 12*q_nocor - k_r;
    dqdc[3] = q_nocor + k_f*sc[5];
    dqdc[4] = q_nocor;
    dqdc[5] = q_nocor + k_f*sc[3];
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+2] += dqdc[k];
        J[10*k+3] -= dqdc[k];
        J[10*k+5] -= dqdc[k];
    }
    J[92] += dqdT; /* dwdot[H2O]/dT */
    J[93] -= dqdT; /* dwdot[H]/dT */
    J[95] -= dqdT; /* dwdot[OH]/dT */

    /*reaction 7: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 3547000000000000
                * exp(-0.40600000000000003 * tc[0] - 0.50321666580471969 * (16599) * invT);
    dlnkfdT = -0.40600000000000003 * invT + 0.50321666580471969 *  (16599)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[1] + g_RT[3] - g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[4] += q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[3];
    J[11] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[13] -= dqdci;               /* dwdot[H]/d[O2] */
    J[14] += dqdci;               /* dwdot[O]/d[O2] */
    J[15] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[1];
    J[31] -= dqdci;               /* dwdot[O2]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[34] += dqdci;               /* dwdot[O]/d[H] */
    J[35] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[41] -= dqdci;               /* dwdot[O2]/d[O] */
    J[43] -= dqdci;               /* dwdot[H]/d[O] */
    J[44] += dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[51] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[53] -= dqdci;               /* dwdot[H]/d[OH] */
    J[54] += dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[O2]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[94] += dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 8: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = 1.0000000000000002e-06 * 50800
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6290) * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  (6290)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[0] - g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[3] += q; /* H */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[3] += dqdci;                /* dwdot[H]/d[H2] */
    J[4] -= dqdci;                /* dwdot[O]/d[H2] */
    J[5] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[30] -= dqdci;               /* dwdot[H2]/d[H] */
    J[33] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[O]/d[H] */
    J[35] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[40] -= dqdci;               /* dwdot[H2]/d[O] */
    J[43] += dqdci;               /* dwdot[H]/d[O] */
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[50] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[53] += dqdci;               /* dwdot[H]/d[OH] */
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[93] += dqdT;                /* dwdot[H]/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 9: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  (3430)  * invT2;
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
    wdot[0] -= q; /* H2 */
    wdot[2] += q; /* H2O */
    wdot[3] += q; /* H */
    wdot[5] -= q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[2] += dqdci;                /* dwdot[H2O]/d[H2] */
    J[3] += dqdci;                /* dwdot[H]/d[H2] */
    J[5] -= dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[20] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[23] += dqdci;               /* dwdot[H]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[2];
    J[30] -= dqdci;               /* dwdot[H2]/d[H] */
    J[32] += dqdci;               /* dwdot[H2O]/d[H] */
    J[33] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[50] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[53] += dqdci;               /* dwdot[H]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[93] += dqdT;                /* dwdot[H]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */

    /*reaction 10: O + H2O <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = 1.0000000000000002e-06 * 2970000
                * exp(2.02 * tc[0] - 0.50321666580471969 * (13400) * invT);
    dlnkfdT = 2.02 * invT + 0.50321666580471969 *  (13400)  * invT2;
    /* reverse */
    phi_r = pow(sc[5], 2.000000);
    Kc = exp(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2.000000*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H2O */
    wdot[4] -= q; /* O */
    wdot[5] += 2 * q; /* OH */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[4];
    J[22] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[24] -= dqdci;               /* dwdot[O]/d[H2O] */
    J[25] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[2];
    J[42] -= dqdci;               /* dwdot[H2O]/d[O] */
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += 2 * dqdci;           /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[5];
    J[52] -= dqdci;               /* dwdot[H2O]/d[OH] */
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[92] -= dqdT;                /* dwdot[H2O]/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += 2 * dqdT;            /* dwdot[OH]/dT */

    /*reaction 11: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = 1.0000000000000002e-06 * 16600000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (823) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (823)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[1];
    Kc = exp(-g_RT[0] - g_RT[1] + g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[0] + h_RT[1]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] += q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[1];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[O2]/d[H2] */
    J[3] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[10] += dqdci;               /* dwdot[H2]/d[O2] */
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[13] -= dqdci;               /* dwdot[H]/d[O2] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[30] += dqdci;               /* dwdot[H2]/d[H] */
    J[31] += dqdci;               /* dwdot[O2]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 12: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = 1.0000000000000002e-06 * 70790000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (295) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (295)  * invT2;
    /* reverse */
    phi_r = pow(sc[5], 2.000000);
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (2.000000*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* H */
    wdot[5] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[5];
    J[53] -= dqdci;               /* dwdot[H]/d[OH] */
    J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[65] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[95] += 2 * dqdT;            /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 13: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = 1.0000000000000002e-06 * 32500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(-g_RT[1] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[5];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[14] -= dqdci;               /* dwdot[O]/d[O2] */
    J[15] += dqdci;               /* dwdot[OH]/d[O2] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[41] += dqdci;               /* dwdot[O2]/d[O] */
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[51] += dqdci;               /* dwdot[O2]/d[OH] */
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[64] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[65] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = 1.0000000000000002e-06 * 28900000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-497) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-497)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[2];
    Kc = exp(-g_RT[1] - g_RT[2] + g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (h_RT[1] + h_RT[2]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[2];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[12] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[15] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[21] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[26] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[51] += dqdci;               /* dwdot[O2]/d[OH] */
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[62] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[65] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (11982) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (11982)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[1] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[16] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[17] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[71] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1629.3) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1629.3)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[1] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[16] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[17] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[71] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 24100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (3970)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + g_RT[3] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[3] -= q; /* H */
    wdot[5] += q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[25] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[32] += dqdci;               /* dwdot[H2O]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += dqdci;               /* dwdot[OH]/d[H] */
    J[37] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[53] -= dqdci;               /* dwdot[H]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[72] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[73] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 48200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7950) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (7950)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[3] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[3] -= q; /* H */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[3] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[30] += dqdci;               /* dwdot[H2]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] += dqdci;               /* dwdot[HO2]/d[H] */
    J[37] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[70] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[73] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = 1.0000000000000002e-06 * 9550000
                * exp(2 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (3970)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    J[46] += dqdci;               /* dwdot[HO2]/d[O] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    J[56] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[65] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               /* dwdot[O]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = 1.0000000000000002e-06 * 1000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[6];
    Kc = exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[26] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[56] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[2];
    J[62] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[65] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[72] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[75] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = 1.0000000000000002e-06 * 580000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (9557) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (9557)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[6];
    Kc = exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[26] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[56] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[2];
    J[62] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[65] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[72] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[75] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    amrex::Real c_R[9], dcRdT[9], e_RT[9];
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
    for (int k = 0; k < 9; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[90+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 9; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 9; ++m) {
            dehmixdc += eh_RT[m]*J[k*10+m];
        }
        J[k*10+9] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[99] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

/*compute an approx to the reaction Jacobian (for preconditioning) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT_SIMPLIFIED(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * HP)
{
    amrex::Real c[9];

    for (int k=0; k<9; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<9; k++) {
        J[90+k] *= 1.e-6;
        J[k*10+9] *= 1.e6;
    }

    return;
}

#ifdef AMREX_USE_GPU
/*compute the reaction Jacobian on GPU */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void aJacobian(amrex::Real * J, amrex::Real * sc, amrex::Real T, int consP)
{

    for (int i=0; i<100; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[9];
    for (int k=0; k<9; k++) {
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
    for (int k = 0; k < 9; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[9];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[9];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[9];
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
    alpha = mixture + ( 2 - 1)*sc[0] + ( 11 - 1)*sc[2] + ( 0.78000000000000003 - 1)*sc[1];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 1475000000000
                * exp(0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0.59999999999999998 * invT + 0.50321666580471969 *  0  * invT2;
    /* pressure-fall-off */
    k_0 = 6.366e+20 * exp(-1.72 * tc[0] - 0.50321666580471969 * (524.79999999999995) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -1.72 * invT + 0.50321666580471969 * (524.79999999999995) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.80000000000000004)*exp(-T/1.0000000000000001e-30);
    Fcent2 = 0.80000000000000004 * exp(-T/1e+30);
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
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[O2]/d[H2] */
        J[3] -= dqdci;                /* dwdot[H]/d[H2] */
        J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
        /* d()/d[O2] */
        dqdci = (0.78000000000000003 - 1)*dcdc_fac + k_f*sc[3];
        J[11] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[13] -= dqdci;               /* dwdot[H]/d[O2] */
        J[16] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (11 - 1)*dcdc_fac;
        J[21] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[26] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[1];
        J[31] -= dqdci;               /* dwdot[O2]/d[H] */
        J[33] -= dqdci;               /* dwdot[H]/d[H] */
        J[36] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[61] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = 0.78000000000000003*dcdc_fac + k_f*sc[3];
        dqdc[2] = 11*dcdc_fac;
        dqdc[3] = dcdc_fac + k_f*sc[1];
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+3] -= dqdc[k];
            J[10*k+6] += dqdc[k];
        }
    }
    J[91] -= dqdT; /* dwdot[O2]/dT */
    J[93] -= dqdT; /* dwdot[H]/dT */
    J[96] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    /* forward */
    phi_f = sc[7];
    k_f = 1 * 295100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (48430) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  48430  * invT2;
    /* pressure-fall-off */
    k_0 = 1.202e+17 * exp(0 * tc[0] - 0.50321666580471969 * (45500) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = 0 * invT + 0.50321666580471969 * (45500) * invT2;
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
    phi_r = pow(sc[5], 2.000000);
    Kc = refC * exp(-g_RT[5] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2.000000*h_RT[5]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[5] += 2 * q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*dcdc_fac;
        J[5] += 2 * dqdci;            /* dwdot[OH]/d[H2] */
        J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*dcdc_fac;
        J[25] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
        J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[OH] */
        dqdci =  - k_r*2.000000*sc[5];
        J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
        J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O2] */
        dqdci =  + k_f;
        J[75] += 2 * dqdci;           /* dwdot[OH]/d[H2O2] */
        J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    }
    else {
        dqdc[0] = 2.5*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = 12*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac - k_r*2.000000*sc[5];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac + k_f;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+5] += 2 * dqdc[k];
            J[10*k+7] -= dqdc[k];
        }
    }
    J[95] += 2 * dqdT; /* dwdot[OH]/dT */
    J[97] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    /* forward */
    phi_f = sc[0];
    k_f = 1.0000000000000002e-06 * 4.577e+19
                * exp(-1.3999999999999999 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.3999999999999999 * invT + 0.50321666580471969 *  104380  * invT2;
    /* reverse */
    phi_r = pow(sc[3], 2.000000);
    Kc = refC * exp(g_RT[0] - g_RT[3] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0]) + (2.000000*h_RT[3]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[3] += 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor + k_f;
        J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
        J[3] += 2 * dqdci;            /* dwdot[H]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor;
        J[20] -= dqdci;               /* dwdot[H2]/d[H2O] */
        J[23] += 2 * dqdci;           /* dwdot[H]/d[H2O] */
        /* d()/d[H] */
        dqdci =  - k_r*2.000000*sc[3];
        J[30] -= dqdci;               /* dwdot[H2]/d[H] */
        J[33] += 2 * dqdci;           /* dwdot[H]/d[H] */
    }
    else {
        dqdc[0] = 2.5*q_nocor + k_f;
        dqdc[1] = q_nocor;
        dqdc[2] = 12*q_nocor;
        dqdc[3] = q_nocor - k_r*2.000000*sc[3];
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+0] -= dqdc[k];
            J[10*k+3] += 2 * dqdc[k];
        }
    }
    J[90] -= dqdT; /* dwdot[H2]/dT */
    J[93] += 2 * dqdT; /* dwdot[H]/dT */

    /*reaction 4: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-12 * 6165000000000000
                * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.5 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1];
    Kc = refCinv * exp(-g_RT[1] + g_RT[4] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[1]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[4] -= 2 * q; /* O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor;
        J[1] += dqdci;                /* dwdot[O2]/d[H2] */
        J[4] += -2 * dqdci;           /* dwdot[O]/d[H2] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[11] += dqdci;               /* dwdot[O2]/d[O2] */
        J[14] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor;
        J[21] += dqdci;               /* dwdot[O2]/d[H2O] */
        J[24] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[4];
        J[41] += dqdci;               /* dwdot[O2]/d[O] */
        J[44] += -2 * dqdci;          /* dwdot[O]/d[O] */
    }
    else {
        dqdc[0] = 2.5*q_nocor;
        dqdc[1] = q_nocor - k_r;
        dqdc[2] = 12*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f*2.000000*sc[4];
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+1] += dqdc[k];
            J[10*k+4] += -2 * dqdc[k];
        }
    }
    J[91] += dqdT; /* dwdot[O2]/dT */
    J[94] += -2 * dqdT; /* dwdot[O]/dT */

    /*reaction 5: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    /* forward */
    phi_f = sc[3]*sc[4];
    k_f = 1.0000000000000002e-12 * 4.714e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[4]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* H */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor;
        J[3] -= dqdci;                /* dwdot[H]/d[H2] */
        J[4] -= dqdci;                /* dwdot[O]/d[H2] */
        J[5] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor;
        J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[24] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[25] += dqdci;               /* dwdot[OH]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[33] -= dqdci;               /* dwdot[H]/d[H] */
        J[34] -= dqdci;               /* dwdot[O]/d[H] */
        J[35] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[3];
        J[43] -= dqdci;               /* dwdot[H]/d[O] */
        J[44] -= dqdci;               /* dwdot[O]/d[O] */
        J[45] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[53] -= dqdci;               /* dwdot[H]/d[OH] */
        J[54] -= dqdci;               /* dwdot[O]/d[OH] */
        J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    }
    else {
        dqdc[0] = 2.5*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = 12*q_nocor;
        dqdc[3] = q_nocor + k_f*sc[4];
        dqdc[4] = q_nocor + k_f*sc[3];
        dqdc[5] = q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+3] -= dqdc[k];
            J[10*k+4] -= dqdc[k];
            J[10*k+5] += dqdc[k];
        }
    }
    J[93] -= dqdT; /* dwdot[H]/dT */
    J[94] -= dqdT; /* dwdot[O]/dT */
    J[95] += dqdT; /* dwdot[OH]/dT */

    /*reaction 6: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    /* forward */
    phi_f = sc[3]*sc[5];
    k_f = 1.0000000000000002e-12 * 3.8000000000000004e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[2];
    Kc = refCinv * exp(-g_RT[2] + g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[2]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[3] -= q; /* H */
    wdot[5] -= q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor;
        J[2] += dqdci;                /* dwdot[H2O]/d[H2] */
        J[3] -= dqdci;                /* dwdot[H]/d[H2] */
        J[5] -= dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor - k_r;
        J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
        J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[32] += dqdci;               /* dwdot[H2O]/d[H] */
        J[33] -= dqdci;               /* dwdot[H]/d[H] */
        J[35] -= dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[3];
        J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
        J[53] -= dqdci;               /* dwdot[H]/d[OH] */
        J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    }
    else {
        dqdc[0] = 2.5*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = 12*q_nocor - k_r;
        dqdc[3] = q_nocor + k_f*sc[5];
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor + k_f*sc[3];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+2] += dqdc[k];
            J[10*k+3] -= dqdc[k];
            J[10*k+5] -= dqdc[k];
        }
    }
    J[92] += dqdT; /* dwdot[H2O]/dT */
    J[93] -= dqdT; /* dwdot[H]/dT */
    J[95] -= dqdT; /* dwdot[OH]/dT */

    /*reaction 7: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 3547000000000000
                * exp(-0.40600000000000003 * tc[0] - 0.50321666580471969 * (16599) * invT);
    dlnkfdT = -0.40600000000000003 * invT + 0.50321666580471969 *  16599  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[1] + g_RT[3] - g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[4] += q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[3];
    J[11] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[13] -= dqdci;               /* dwdot[H]/d[O2] */
    J[14] += dqdci;               /* dwdot[O]/d[O2] */
    J[15] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[1];
    J[31] -= dqdci;               /* dwdot[O2]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[34] += dqdci;               /* dwdot[O]/d[H] */
    J[35] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[41] -= dqdci;               /* dwdot[O2]/d[O] */
    J[43] -= dqdci;               /* dwdot[H]/d[O] */
    J[44] += dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[51] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[53] -= dqdci;               /* dwdot[H]/d[OH] */
    J[54] += dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[O2]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[94] += dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 8: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = 1.0000000000000002e-06 * 50800
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6290) * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  6290  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[0] - g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[3] += q; /* H */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[3] += dqdci;                /* dwdot[H]/d[H2] */
    J[4] -= dqdci;                /* dwdot[O]/d[H2] */
    J[5] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[30] -= dqdci;               /* dwdot[H2]/d[H] */
    J[33] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[O]/d[H] */
    J[35] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[40] -= dqdci;               /* dwdot[H2]/d[O] */
    J[43] += dqdci;               /* dwdot[H]/d[O] */
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[50] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[53] += dqdci;               /* dwdot[H]/d[OH] */
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[93] += dqdT;                /* dwdot[H]/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 9: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  3430  * invT2;
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
    wdot[0] -= q; /* H2 */
    wdot[2] += q; /* H2O */
    wdot[3] += q; /* H */
    wdot[5] -= q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[2] += dqdci;                /* dwdot[H2O]/d[H2] */
    J[3] += dqdci;                /* dwdot[H]/d[H2] */
    J[5] -= dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[20] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[23] += dqdci;               /* dwdot[H]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[2];
    J[30] -= dqdci;               /* dwdot[H2]/d[H] */
    J[32] += dqdci;               /* dwdot[H2O]/d[H] */
    J[33] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[50] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[53] += dqdci;               /* dwdot[H]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[93] += dqdT;                /* dwdot[H]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */

    /*reaction 10: O + H2O <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = 1.0000000000000002e-06 * 2970000
                * exp(2.02 * tc[0] - 0.50321666580471969 * (13400) * invT);
    dlnkfdT = 2.02 * invT + 0.50321666580471969 *  13400  * invT2;
    /* reverse */
    phi_r = pow(sc[5], 2.000000);
    Kc = exp(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2.000000*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H2O */
    wdot[4] -= q; /* O */
    wdot[5] += 2 * q; /* OH */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[4];
    J[22] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[24] -= dqdci;               /* dwdot[O]/d[H2O] */
    J[25] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[2];
    J[42] -= dqdci;               /* dwdot[H2O]/d[O] */
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += 2 * dqdci;           /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[5];
    J[52] -= dqdci;               /* dwdot[H2O]/d[OH] */
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[92] -= dqdT;                /* dwdot[H2O]/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += 2 * dqdT;            /* dwdot[OH]/dT */

    /*reaction 11: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = 1.0000000000000002e-06 * 16600000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (823) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  823  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[1];
    Kc = exp(-g_RT[0] - g_RT[1] + g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[0] + h_RT[1]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] += q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[1];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[O2]/d[H2] */
    J[3] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[10] += dqdci;               /* dwdot[H2]/d[O2] */
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[13] -= dqdci;               /* dwdot[H]/d[O2] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[30] += dqdci;               /* dwdot[H2]/d[H] */
    J[31] += dqdci;               /* dwdot[O2]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 12: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = 1.0000000000000002e-06 * 70790000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (295) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  295  * invT2;
    /* reverse */
    phi_r = pow(sc[5], 2.000000);
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (2.000000*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* H */
    wdot[5] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[5];
    J[53] -= dqdci;               /* dwdot[H]/d[OH] */
    J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[65] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[95] += 2 * dqdT;            /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 13: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = 1.0000000000000002e-06 * 32500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(-g_RT[1] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[5];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[14] -= dqdci;               /* dwdot[O]/d[O2] */
    J[15] += dqdci;               /* dwdot[OH]/d[O2] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[41] += dqdci;               /* dwdot[O2]/d[O] */
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[51] += dqdci;               /* dwdot[O2]/d[OH] */
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[64] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[65] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = 1.0000000000000002e-06 * 28900000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-497) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -497  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[2];
    Kc = exp(-g_RT[1] - g_RT[2] + g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (h_RT[1] + h_RT[2]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[2];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[12] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[15] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[21] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[26] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[51] += dqdci;               /* dwdot[O2]/d[OH] */
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[62] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[65] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (11982) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  11982  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[1] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[16] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[17] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[71] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1629.3) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1629.3  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[1] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[16] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[17] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[71] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 24100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3970  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + g_RT[3] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[3] -= q; /* H */
    wdot[5] += q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[25] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[32] += dqdci;               /* dwdot[H2O]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += dqdci;               /* dwdot[OH]/d[H] */
    J[37] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[53] -= dqdci;               /* dwdot[H]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[72] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[73] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 48200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7950) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  7950  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[3] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[3] -= q; /* H */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[3] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[30] += dqdci;               /* dwdot[H2]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] += dqdci;               /* dwdot[HO2]/d[H] */
    J[37] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[70] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[73] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = 1.0000000000000002e-06 * 9550000
                * exp(2 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  3970  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    J[46] += dqdci;               /* dwdot[HO2]/d[O] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    J[56] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[65] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               /* dwdot[O]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = 1.0000000000000002e-06 * 1000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[6];
    Kc = exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[26] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[56] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[2];
    J[62] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[65] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[72] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[75] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = 1.0000000000000002e-06 * 580000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (9557) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  9557  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[6];
    Kc = exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[26] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[56] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[2];
    J[62] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[65] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[72] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[75] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    amrex::Real c_R[9], dcRdT[9], e_RT[9];
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
    for (int k = 0; k < 9; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[90+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 9; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 9; ++m) {
            dehmixdc += eh_RT[m]*J[k*10+m];
        }
        J[k*10+9] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[99] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;

return;
}
#endif


/*compute the reaction Jacobian */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * consP)
{
    amrex::Real c[9];

    for (int k=0; k<9; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<9; k++) {
        J[90+k] *= 1.e-6;
        J[k*10+9] *= 1.e6;
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
