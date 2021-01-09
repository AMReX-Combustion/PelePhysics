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

    extern amrex::Real fwd_A[84], fwd_beta[84], fwd_Ea[84];
    extern amrex::Real low_A[84], low_beta[84], low_Ea[84];
    extern amrex::Real rev_A[84], rev_beta[84], rev_Ea[84];
    extern amrex::Real troe_a[84],troe_Ts[84], troe_Tss[84], troe_Tsss[84];
    extern amrex::Real sri_a[84], sri_b[84], sri_c[84], sri_d[84], sri_e[84];
    extern amrex::Real activation_units[84], prefactor_units[84], phase_units[84];
    extern int is_PD[84], troe_len[84], sri_len[84], nTB[84], *TBid[84];
    extern amrex::Real *TB[84];
    extern std::vector<std::vector<amrex::Real>> kiv; 
    extern std::vector<std::vector<amrex::Real>> nuv; 
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
void vcomp_wdot_1_50(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
                amrex::Real *  k_f_s, amrex::Real *  Kc_s,
                amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T);
void vcomp_wdot_51_84(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
                amrex::Real *  k_f_s, amrex::Real *  Kc_s,
                amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T);
/* MISC  */
void CKINU(int * i, int * nspec, int * ki, int * nu);
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
#endif


/* ALWAYS on CPU stuff -- can have different def depending on if we are CPU or GPU based. Defined in mechanism.cpp  */
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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKINDX(int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 5;
    *kk = 21;
    *ii = 84;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/* inverse molecular weights  */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
    imw_new[0] = 1.0/2.015940;/*H2 */
    imw_new[1] = 1.0/1.007970;/*H */
    imw_new[2] = 1.0/15.999400;/*O */
    imw_new[3] = 1.0/31.998800;/*O2 */
    imw_new[4] = 1.0/17.007370;/*OH */
    imw_new[5] = 1.0/18.015340;/*H2O */
    imw_new[6] = 1.0/33.006770;/*HO2 */
    imw_new[7] = 1.0/14.027090;/*CH2 */
    imw_new[8] = 1.0/14.027090;/*CH2(S) */
    imw_new[9] = 1.0/15.035060;/*CH3 */
    imw_new[10] = 1.0/16.043030;/*CH4 */
    imw_new[11] = 1.0/28.010550;/*CO */
    imw_new[12] = 1.0/44.009950;/*CO2 */
    imw_new[13] = 1.0/29.018520;/*HCO */
    imw_new[14] = 1.0/30.026490;/*CH2O */
    imw_new[15] = 1.0/31.034460;/*CH3O */
    imw_new[16] = 1.0/28.054180;/*C2H4 */
    imw_new[17] = 1.0/29.062150;/*C2H5 */
    imw_new[18] = 1.0/30.070120;/*C2H6 */
    imw_new[19] = 1.0/28.013400;/*N2 */
    imw_new[20] = 1.0/39.948000;/*AR */
}

/* molecular weights  */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void get_mw(amrex::Real *mw_new){
    mw_new[0] = 2.015940;/*H2 */
    mw_new[1] = 1.007970;/*H */
    mw_new[2] = 15.999400;/*O */
    mw_new[3] = 31.998800;/*O2 */
    mw_new[4] = 17.007370;/*OH */
    mw_new[5] = 18.015340;/*H2O */
    mw_new[6] = 33.006770;/*HO2 */
    mw_new[7] = 14.027090;/*CH2 */
    mw_new[8] = 14.027090;/*CH2(S) */
    mw_new[9] = 15.035060;/*CH3 */
    mw_new[10] = 16.043030;/*CH4 */
    mw_new[11] = 28.010550;/*CO */
    mw_new[12] = 44.009950;/*CO2 */
    mw_new[13] = 29.018520;/*HCO */
    mw_new[14] = 30.026490;/*CH2O */
    mw_new[15] = 31.034460;/*CH3O */
    mw_new[16] = 28.054180;/*C2H4 */
    mw_new[17] = 29.062150;/*C2H5 */
    mw_new[18] = 30.070120;/*C2H6 */
    mw_new[19] = 28.013400;/*N2 */
    mw_new[20] = 39.948000;/*AR */
}



/* Returns R, Rc, Patm */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKRP(amrex::Real *  ru, amrex::Real *  ruc, amrex::Real *  pa)
{
     *ru  = 8.31446261815324e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKNU(int * kdim,  int * nuki)
{
    int kd = (*kdim); 
    /*Zero nukis */
    for (int id = 0; id < 21 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    nuki[ 1 * kd + 0 ] += -1.000000 ;
    nuki[ 7 * kd + 0 ] += -1.000000 ;
    nuki[ 9 * kd + 0 ] += +1.000000 ;

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    nuki[ 1 * kd + 1 ] += -1.000000 ;
    nuki[ 9 * kd + 1 ] += -1.000000 ;
    nuki[ 10 * kd + 1 ] += +1.000000 ;

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    nuki[ 1 * kd + 2 ] += -1.000000 ;
    nuki[ 13 * kd + 2 ] += -1.000000 ;
    nuki[ 14 * kd + 2 ] += +1.000000 ;

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    nuki[ 1 * kd + 3 ] += -1.000000 ;
    nuki[ 14 * kd + 3 ] += -1.000000 ;
    nuki[ 15 * kd + 3 ] += +1.000000 ;

    /*reaction 5: H + C2H4 (+M) <=> C2H5 (+M) */
    nuki[ 1 * kd + 4 ] += -1.000000 ;
    nuki[ 16 * kd + 4 ] += -1.000000 ;
    nuki[ 17 * kd + 4 ] += +1.000000 ;

    /*reaction 6: H + C2H5 (+M) <=> C2H6 (+M) */
    nuki[ 1 * kd + 5 ] += -1.000000 ;
    nuki[ 17 * kd + 5 ] += -1.000000 ;
    nuki[ 18 * kd + 5 ] += +1.000000 ;

    /*reaction 7: H2 + CO (+M) <=> CH2O (+M) */
    nuki[ 0 * kd + 6 ] += -1.000000 ;
    nuki[ 11 * kd + 6 ] += -1.000000 ;
    nuki[ 14 * kd + 6 ] += +1.000000 ;

    /*reaction 8: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    nuki[ 9 * kd + 7 ] += -2.000000 ;
    nuki[ 18 * kd + 7 ] += +1.000000 ;

    /*reaction 9: O + H + M <=> OH + M */
    nuki[ 2 * kd + 8 ] += -1.000000 ;
    nuki[ 1 * kd + 8 ] += -1.000000 ;
    nuki[ 4 * kd + 8 ] += +1.000000 ;

    /*reaction 10: O + CO + M <=> CO2 + M */
    nuki[ 2 * kd + 9 ] += -1.000000 ;
    nuki[ 11 * kd + 9 ] += -1.000000 ;
    nuki[ 12 * kd + 9 ] += +1.000000 ;

    /*reaction 11: H + O2 + M <=> HO2 + M */
    nuki[ 1 * kd + 10 ] += -1.000000 ;
    nuki[ 3 * kd + 10 ] += -1.000000 ;
    nuki[ 6 * kd + 10 ] += +1.000000 ;

    /*reaction 12: 2.000000 H + M <=> H2 + M */
    nuki[ 1 * kd + 11 ] += -2.000000 ;
    nuki[ 0 * kd + 11 ] += +1.000000 ;

    /*reaction 13: H + OH + M <=> H2O + M */
    nuki[ 1 * kd + 12 ] += -1.000000 ;
    nuki[ 4 * kd + 12 ] += -1.000000 ;
    nuki[ 5 * kd + 12 ] += +1.000000 ;

    /*reaction 14: HCO + M <=> H + CO + M */
    nuki[ 13 * kd + 13 ] += -1.000000 ;
    nuki[ 1 * kd + 13 ] += +1.000000 ;
    nuki[ 11 * kd + 13 ] += +1.000000 ;

    /*reaction 15: O + H2 <=> H + OH */
    nuki[ 2 * kd + 14 ] += -1.000000 ;
    nuki[ 0 * kd + 14 ] += -1.000000 ;
    nuki[ 1 * kd + 14 ] += +1.000000 ;
    nuki[ 4 * kd + 14 ] += +1.000000 ;

    /*reaction 16: O + HO2 <=> OH + O2 */
    nuki[ 2 * kd + 15 ] += -1.000000 ;
    nuki[ 6 * kd + 15 ] += -1.000000 ;
    nuki[ 4 * kd + 15 ] += +1.000000 ;
    nuki[ 3 * kd + 15 ] += +1.000000 ;

    /*reaction 17: O + CH2 <=> H + HCO */
    nuki[ 2 * kd + 16 ] += -1.000000 ;
    nuki[ 7 * kd + 16 ] += -1.000000 ;
    nuki[ 1 * kd + 16 ] += +1.000000 ;
    nuki[ 13 * kd + 16 ] += +1.000000 ;

    /*reaction 18: O + CH2(S) <=> H + HCO */
    nuki[ 2 * kd + 17 ] += -1.000000 ;
    nuki[ 8 * kd + 17 ] += -1.000000 ;
    nuki[ 1 * kd + 17 ] += +1.000000 ;
    nuki[ 13 * kd + 17 ] += +1.000000 ;

    /*reaction 19: O + CH3 <=> H + CH2O */
    nuki[ 2 * kd + 18 ] += -1.000000 ;
    nuki[ 9 * kd + 18 ] += -1.000000 ;
    nuki[ 1 * kd + 18 ] += +1.000000 ;
    nuki[ 14 * kd + 18 ] += +1.000000 ;

    /*reaction 20: O + CH4 <=> OH + CH3 */
    nuki[ 2 * kd + 19 ] += -1.000000 ;
    nuki[ 10 * kd + 19 ] += -1.000000 ;
    nuki[ 4 * kd + 19 ] += +1.000000 ;
    nuki[ 9 * kd + 19 ] += +1.000000 ;

    /*reaction 21: O + HCO <=> OH + CO */
    nuki[ 2 * kd + 20 ] += -1.000000 ;
    nuki[ 13 * kd + 20 ] += -1.000000 ;
    nuki[ 4 * kd + 20 ] += +1.000000 ;
    nuki[ 11 * kd + 20 ] += +1.000000 ;

    /*reaction 22: O + HCO <=> H + CO2 */
    nuki[ 2 * kd + 21 ] += -1.000000 ;
    nuki[ 13 * kd + 21 ] += -1.000000 ;
    nuki[ 1 * kd + 21 ] += +1.000000 ;
    nuki[ 12 * kd + 21 ] += +1.000000 ;

    /*reaction 23: O + CH2O <=> OH + HCO */
    nuki[ 2 * kd + 22 ] += -1.000000 ;
    nuki[ 14 * kd + 22 ] += -1.000000 ;
    nuki[ 4 * kd + 22 ] += +1.000000 ;
    nuki[ 13 * kd + 22 ] += +1.000000 ;

    /*reaction 24: O + C2H4 <=> CH3 + HCO */
    nuki[ 2 * kd + 23 ] += -1.000000 ;
    nuki[ 16 * kd + 23 ] += -1.000000 ;
    nuki[ 9 * kd + 23 ] += +1.000000 ;
    nuki[ 13 * kd + 23 ] += +1.000000 ;

    /*reaction 25: O + C2H5 <=> CH3 + CH2O */
    nuki[ 2 * kd + 24 ] += -1.000000 ;
    nuki[ 17 * kd + 24 ] += -1.000000 ;
    nuki[ 9 * kd + 24 ] += +1.000000 ;
    nuki[ 14 * kd + 24 ] += +1.000000 ;

    /*reaction 26: O + C2H6 <=> OH + C2H5 */
    nuki[ 2 * kd + 25 ] += -1.000000 ;
    nuki[ 18 * kd + 25 ] += -1.000000 ;
    nuki[ 4 * kd + 25 ] += +1.000000 ;
    nuki[ 17 * kd + 25 ] += +1.000000 ;

    /*reaction 27: O2 + CO <=> O + CO2 */
    nuki[ 3 * kd + 26 ] += -1.000000 ;
    nuki[ 11 * kd + 26 ] += -1.000000 ;
    nuki[ 2 * kd + 26 ] += +1.000000 ;
    nuki[ 12 * kd + 26 ] += +1.000000 ;

    /*reaction 28: O2 + CH2O <=> HO2 + HCO */
    nuki[ 3 * kd + 27 ] += -1.000000 ;
    nuki[ 14 * kd + 27 ] += -1.000000 ;
    nuki[ 6 * kd + 27 ] += +1.000000 ;
    nuki[ 13 * kd + 27 ] += +1.000000 ;

    /*reaction 29: H + 2.000000 O2 <=> HO2 + O2 */
    nuki[ 1 * kd + 28 ] += -1.000000 ;
    nuki[ 3 * kd + 28 ] += -2.000000 ;
    nuki[ 6 * kd + 28 ] += +1.000000 ;
    nuki[ 3 * kd + 28 ] += +1.000000 ;

    /*reaction 30: H + O2 + H2O <=> HO2 + H2O */
    nuki[ 1 * kd + 29 ] += -1.000000 ;
    nuki[ 3 * kd + 29 ] += -1.000000 ;
    nuki[ 5 * kd + 29 ] += -1.000000 ;
    nuki[ 6 * kd + 29 ] += +1.000000 ;
    nuki[ 5 * kd + 29 ] += +1.000000 ;

    /*reaction 31: H + O2 + N2 <=> HO2 + N2 */
    nuki[ 1 * kd + 30 ] += -1.000000 ;
    nuki[ 3 * kd + 30 ] += -1.000000 ;
    nuki[ 19 * kd + 30 ] += -1.000000 ;
    nuki[ 6 * kd + 30 ] += +1.000000 ;
    nuki[ 19 * kd + 30 ] += +1.000000 ;

    /*reaction 32: H + O2 + AR <=> HO2 + AR */
    nuki[ 1 * kd + 31 ] += -1.000000 ;
    nuki[ 3 * kd + 31 ] += -1.000000 ;
    nuki[ 20 * kd + 31 ] += -1.000000 ;
    nuki[ 6 * kd + 31 ] += +1.000000 ;
    nuki[ 20 * kd + 31 ] += +1.000000 ;

    /*reaction 33: H + O2 <=> O + OH */
    nuki[ 1 * kd + 32 ] += -1.000000 ;
    nuki[ 3 * kd + 32 ] += -1.000000 ;
    nuki[ 2 * kd + 32 ] += +1.000000 ;
    nuki[ 4 * kd + 32 ] += +1.000000 ;

    /*reaction 34: 2.000000 H + H2 <=> 2.000000 H2 */
    nuki[ 1 * kd + 33 ] += -2.000000 ;
    nuki[ 0 * kd + 33 ] += -1.000000 ;
    nuki[ 0 * kd + 33 ] += +2.000000 ;

    /*reaction 35: 2.000000 H + H2O <=> H2 + H2O */
    nuki[ 1 * kd + 34 ] += -2.000000 ;
    nuki[ 5 * kd + 34 ] += -1.000000 ;
    nuki[ 0 * kd + 34 ] += +1.000000 ;
    nuki[ 5 * kd + 34 ] += +1.000000 ;

    /*reaction 36: 2.000000 H + CO2 <=> H2 + CO2 */
    nuki[ 1 * kd + 35 ] += -2.000000 ;
    nuki[ 12 * kd + 35 ] += -1.000000 ;
    nuki[ 0 * kd + 35 ] += +1.000000 ;
    nuki[ 12 * kd + 35 ] += +1.000000 ;

    /*reaction 37: H + HO2 <=> O2 + H2 */
    nuki[ 1 * kd + 36 ] += -1.000000 ;
    nuki[ 6 * kd + 36 ] += -1.000000 ;
    nuki[ 3 * kd + 36 ] += +1.000000 ;
    nuki[ 0 * kd + 36 ] += +1.000000 ;

    /*reaction 38: H + HO2 <=> 2.000000 OH */
    nuki[ 1 * kd + 37 ] += -1.000000 ;
    nuki[ 6 * kd + 37 ] += -1.000000 ;
    nuki[ 4 * kd + 37 ] += +2.000000 ;

    /*reaction 39: H + CH4 <=> CH3 + H2 */
    nuki[ 1 * kd + 38 ] += -1.000000 ;
    nuki[ 10 * kd + 38 ] += -1.000000 ;
    nuki[ 9 * kd + 38 ] += +1.000000 ;
    nuki[ 0 * kd + 38 ] += +1.000000 ;

    /*reaction 40: H + HCO <=> H2 + CO */
    nuki[ 1 * kd + 39 ] += -1.000000 ;
    nuki[ 13 * kd + 39 ] += -1.000000 ;
    nuki[ 0 * kd + 39 ] += +1.000000 ;
    nuki[ 11 * kd + 39 ] += +1.000000 ;

    /*reaction 41: H + CH2O <=> HCO + H2 */
    nuki[ 1 * kd + 40 ] += -1.000000 ;
    nuki[ 14 * kd + 40 ] += -1.000000 ;
    nuki[ 13 * kd + 40 ] += +1.000000 ;
    nuki[ 0 * kd + 40 ] += +1.000000 ;

    /*reaction 42: H + CH3O <=> OH + CH3 */
    nuki[ 1 * kd + 41 ] += -1.000000 ;
    nuki[ 15 * kd + 41 ] += -1.000000 ;
    nuki[ 4 * kd + 41 ] += +1.000000 ;
    nuki[ 9 * kd + 41 ] += +1.000000 ;

    /*reaction 43: H + C2H6 <=> C2H5 + H2 */
    nuki[ 1 * kd + 42 ] += -1.000000 ;
    nuki[ 18 * kd + 42 ] += -1.000000 ;
    nuki[ 17 * kd + 42 ] += +1.000000 ;
    nuki[ 0 * kd + 42 ] += +1.000000 ;

    /*reaction 44: OH + H2 <=> H + H2O */
    nuki[ 4 * kd + 43 ] += -1.000000 ;
    nuki[ 0 * kd + 43 ] += -1.000000 ;
    nuki[ 1 * kd + 43 ] += +1.000000 ;
    nuki[ 5 * kd + 43 ] += +1.000000 ;

    /*reaction 45: 2.000000 OH <=> O + H2O */
    nuki[ 4 * kd + 44 ] += -2.000000 ;
    nuki[ 2 * kd + 44 ] += +1.000000 ;
    nuki[ 5 * kd + 44 ] += +1.000000 ;

    /*reaction 46: OH + HO2 <=> O2 + H2O */
    nuki[ 4 * kd + 45 ] += -1.000000 ;
    nuki[ 6 * kd + 45 ] += -1.000000 ;
    nuki[ 3 * kd + 45 ] += +1.000000 ;
    nuki[ 5 * kd + 45 ] += +1.000000 ;

    /*reaction 47: OH + CH2 <=> H + CH2O */
    nuki[ 4 * kd + 46 ] += -1.000000 ;
    nuki[ 7 * kd + 46 ] += -1.000000 ;
    nuki[ 1 * kd + 46 ] += +1.000000 ;
    nuki[ 14 * kd + 46 ] += +1.000000 ;

    /*reaction 48: OH + CH2(S) <=> H + CH2O */
    nuki[ 4 * kd + 47 ] += -1.000000 ;
    nuki[ 8 * kd + 47 ] += -1.000000 ;
    nuki[ 1 * kd + 47 ] += +1.000000 ;
    nuki[ 14 * kd + 47 ] += +1.000000 ;

    /*reaction 49: OH + CH3 <=> CH2 + H2O */
    nuki[ 4 * kd + 48 ] += -1.000000 ;
    nuki[ 9 * kd + 48 ] += -1.000000 ;
    nuki[ 7 * kd + 48 ] += +1.000000 ;
    nuki[ 5 * kd + 48 ] += +1.000000 ;

    /*reaction 50: OH + CH3 <=> CH2(S) + H2O */
    nuki[ 4 * kd + 49 ] += -1.000000 ;
    nuki[ 9 * kd + 49 ] += -1.000000 ;
    nuki[ 8 * kd + 49 ] += +1.000000 ;
    nuki[ 5 * kd + 49 ] += +1.000000 ;

    /*reaction 51: OH + CH4 <=> CH3 + H2O */
    nuki[ 4 * kd + 50 ] += -1.000000 ;
    nuki[ 10 * kd + 50 ] += -1.000000 ;
    nuki[ 9 * kd + 50 ] += +1.000000 ;
    nuki[ 5 * kd + 50 ] += +1.000000 ;

    /*reaction 52: OH + CO <=> H + CO2 */
    nuki[ 4 * kd + 51 ] += -1.000000 ;
    nuki[ 11 * kd + 51 ] += -1.000000 ;
    nuki[ 1 * kd + 51 ] += +1.000000 ;
    nuki[ 12 * kd + 51 ] += +1.000000 ;

    /*reaction 53: OH + HCO <=> H2O + CO */
    nuki[ 4 * kd + 52 ] += -1.000000 ;
    nuki[ 13 * kd + 52 ] += -1.000000 ;
    nuki[ 5 * kd + 52 ] += +1.000000 ;
    nuki[ 11 * kd + 52 ] += +1.000000 ;

    /*reaction 54: OH + CH2O <=> HCO + H2O */
    nuki[ 4 * kd + 53 ] += -1.000000 ;
    nuki[ 14 * kd + 53 ] += -1.000000 ;
    nuki[ 13 * kd + 53 ] += +1.000000 ;
    nuki[ 5 * kd + 53 ] += +1.000000 ;

    /*reaction 55: OH + C2H6 <=> C2H5 + H2O */
    nuki[ 4 * kd + 54 ] += -1.000000 ;
    nuki[ 18 * kd + 54 ] += -1.000000 ;
    nuki[ 17 * kd + 54 ] += +1.000000 ;
    nuki[ 5 * kd + 54 ] += +1.000000 ;

    /*reaction 56: HO2 + CH2 <=> OH + CH2O */
    nuki[ 6 * kd + 55 ] += -1.000000 ;
    nuki[ 7 * kd + 55 ] += -1.000000 ;
    nuki[ 4 * kd + 55 ] += +1.000000 ;
    nuki[ 14 * kd + 55 ] += +1.000000 ;

    /*reaction 57: HO2 + CH3 <=> O2 + CH4 */
    nuki[ 6 * kd + 56 ] += -1.000000 ;
    nuki[ 9 * kd + 56 ] += -1.000000 ;
    nuki[ 3 * kd + 56 ] += +1.000000 ;
    nuki[ 10 * kd + 56 ] += +1.000000 ;

    /*reaction 58: HO2 + CH3 <=> OH + CH3O */
    nuki[ 6 * kd + 57 ] += -1.000000 ;
    nuki[ 9 * kd + 57 ] += -1.000000 ;
    nuki[ 4 * kd + 57 ] += +1.000000 ;
    nuki[ 15 * kd + 57 ] += +1.000000 ;

    /*reaction 59: HO2 + CO <=> OH + CO2 */
    nuki[ 6 * kd + 58 ] += -1.000000 ;
    nuki[ 11 * kd + 58 ] += -1.000000 ;
    nuki[ 4 * kd + 58 ] += +1.000000 ;
    nuki[ 12 * kd + 58 ] += +1.000000 ;

    /*reaction 60: CH2 + O2 <=> OH + HCO */
    nuki[ 7 * kd + 59 ] += -1.000000 ;
    nuki[ 3 * kd + 59 ] += -1.000000 ;
    nuki[ 4 * kd + 59 ] += +1.000000 ;
    nuki[ 13 * kd + 59 ] += +1.000000 ;

    /*reaction 61: CH2 + H2 <=> H + CH3 */
    nuki[ 7 * kd + 60 ] += -1.000000 ;
    nuki[ 0 * kd + 60 ] += -1.000000 ;
    nuki[ 1 * kd + 60 ] += +1.000000 ;
    nuki[ 9 * kd + 60 ] += +1.000000 ;

    /*reaction 62: CH2 + CH3 <=> H + C2H4 */
    nuki[ 7 * kd + 61 ] += -1.000000 ;
    nuki[ 9 * kd + 61 ] += -1.000000 ;
    nuki[ 1 * kd + 61 ] += +1.000000 ;
    nuki[ 16 * kd + 61 ] += +1.000000 ;

    /*reaction 63: CH2 + CH4 <=> 2.000000 CH3 */
    nuki[ 7 * kd + 62 ] += -1.000000 ;
    nuki[ 10 * kd + 62 ] += -1.000000 ;
    nuki[ 9 * kd + 62 ] += +2.000000 ;

    /*reaction 64: CH2(S) + N2 <=> CH2 + N2 */
    nuki[ 8 * kd + 63 ] += -1.000000 ;
    nuki[ 19 * kd + 63 ] += -1.000000 ;
    nuki[ 7 * kd + 63 ] += +1.000000 ;
    nuki[ 19 * kd + 63 ] += +1.000000 ;

    /*reaction 65: CH2(S) + AR <=> CH2 + AR */
    nuki[ 8 * kd + 64 ] += -1.000000 ;
    nuki[ 20 * kd + 64 ] += -1.000000 ;
    nuki[ 7 * kd + 64 ] += +1.000000 ;
    nuki[ 20 * kd + 64 ] += +1.000000 ;

    /*reaction 66: CH2(S) + O2 <=> H + OH + CO */
    nuki[ 8 * kd + 65 ] += -1.000000 ;
    nuki[ 3 * kd + 65 ] += -1.000000 ;
    nuki[ 1 * kd + 65 ] += +1.000000 ;
    nuki[ 4 * kd + 65 ] += +1.000000 ;
    nuki[ 11 * kd + 65 ] += +1.000000 ;

    /*reaction 67: CH2(S) + O2 <=> CO + H2O */
    nuki[ 8 * kd + 66 ] += -1.000000 ;
    nuki[ 3 * kd + 66 ] += -1.000000 ;
    nuki[ 11 * kd + 66 ] += +1.000000 ;
    nuki[ 5 * kd + 66 ] += +1.000000 ;

    /*reaction 68: CH2(S) + H2 <=> CH3 + H */
    nuki[ 8 * kd + 67 ] += -1.000000 ;
    nuki[ 0 * kd + 67 ] += -1.000000 ;
    nuki[ 9 * kd + 67 ] += +1.000000 ;
    nuki[ 1 * kd + 67 ] += +1.000000 ;

    /*reaction 69: CH2(S) + H2O <=> CH2 + H2O */
    nuki[ 8 * kd + 68 ] += -1.000000 ;
    nuki[ 5 * kd + 68 ] += -1.000000 ;
    nuki[ 7 * kd + 68 ] += +1.000000 ;
    nuki[ 5 * kd + 68 ] += +1.000000 ;

    /*reaction 70: CH2(S) + CH3 <=> H + C2H4 */
    nuki[ 8 * kd + 69 ] += -1.000000 ;
    nuki[ 9 * kd + 69 ] += -1.000000 ;
    nuki[ 1 * kd + 69 ] += +1.000000 ;
    nuki[ 16 * kd + 69 ] += +1.000000 ;

    /*reaction 71: CH2(S) + CH4 <=> 2.000000 CH3 */
    nuki[ 8 * kd + 70 ] += -1.000000 ;
    nuki[ 10 * kd + 70 ] += -1.000000 ;
    nuki[ 9 * kd + 70 ] += +2.000000 ;

    /*reaction 72: CH2(S) + CO <=> CH2 + CO */
    nuki[ 8 * kd + 71 ] += -1.000000 ;
    nuki[ 11 * kd + 71 ] += -1.000000 ;
    nuki[ 7 * kd + 71 ] += +1.000000 ;
    nuki[ 11 * kd + 71 ] += +1.000000 ;

    /*reaction 73: CH2(S) + CO2 <=> CH2 + CO2 */
    nuki[ 8 * kd + 72 ] += -1.000000 ;
    nuki[ 12 * kd + 72 ] += -1.000000 ;
    nuki[ 7 * kd + 72 ] += +1.000000 ;
    nuki[ 12 * kd + 72 ] += +1.000000 ;

    /*reaction 74: CH2(S) + CO2 <=> CO + CH2O */
    nuki[ 8 * kd + 73 ] += -1.000000 ;
    nuki[ 12 * kd + 73 ] += -1.000000 ;
    nuki[ 11 * kd + 73 ] += +1.000000 ;
    nuki[ 14 * kd + 73 ] += +1.000000 ;

    /*reaction 75: CH3 + O2 <=> O + CH3O */
    nuki[ 9 * kd + 74 ] += -1.000000 ;
    nuki[ 3 * kd + 74 ] += -1.000000 ;
    nuki[ 2 * kd + 74 ] += +1.000000 ;
    nuki[ 15 * kd + 74 ] += +1.000000 ;

    /*reaction 76: CH3 + O2 <=> OH + CH2O */
    nuki[ 9 * kd + 75 ] += -1.000000 ;
    nuki[ 3 * kd + 75 ] += -1.000000 ;
    nuki[ 4 * kd + 75 ] += +1.000000 ;
    nuki[ 14 * kd + 75 ] += +1.000000 ;

    /*reaction 77: 2.000000 CH3 <=> H + C2H5 */
    nuki[ 9 * kd + 76 ] += -2.000000 ;
    nuki[ 1 * kd + 76 ] += +1.000000 ;
    nuki[ 17 * kd + 76 ] += +1.000000 ;

    /*reaction 78: CH3 + HCO <=> CH4 + CO */
    nuki[ 9 * kd + 77 ] += -1.000000 ;
    nuki[ 13 * kd + 77 ] += -1.000000 ;
    nuki[ 10 * kd + 77 ] += +1.000000 ;
    nuki[ 11 * kd + 77 ] += +1.000000 ;

    /*reaction 79: CH3 + CH2O <=> HCO + CH4 */
    nuki[ 9 * kd + 78 ] += -1.000000 ;
    nuki[ 14 * kd + 78 ] += -1.000000 ;
    nuki[ 13 * kd + 78 ] += +1.000000 ;
    nuki[ 10 * kd + 78 ] += +1.000000 ;

    /*reaction 80: CH3 + C2H6 <=> C2H5 + CH4 */
    nuki[ 9 * kd + 79 ] += -1.000000 ;
    nuki[ 18 * kd + 79 ] += -1.000000 ;
    nuki[ 17 * kd + 79 ] += +1.000000 ;
    nuki[ 10 * kd + 79 ] += +1.000000 ;

    /*reaction 81: HCO + H2O <=> H + CO + H2O */
    nuki[ 13 * kd + 80 ] += -1.000000 ;
    nuki[ 5 * kd + 80 ] += -1.000000 ;
    nuki[ 1 * kd + 80 ] += +1.000000 ;
    nuki[ 11 * kd + 80 ] += +1.000000 ;
    nuki[ 5 * kd + 80 ] += +1.000000 ;

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    nuki[ 13 * kd + 81 ] += -1.000000 ;
    nuki[ 3 * kd + 81 ] += -1.000000 ;
    nuki[ 6 * kd + 81 ] += +1.000000 ;
    nuki[ 11 * kd + 81 ] += +1.000000 ;

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    nuki[ 15 * kd + 82 ] += -1.000000 ;
    nuki[ 3 * kd + 82 ] += -1.000000 ;
    nuki[ 6 * kd + 82 ] += +1.000000 ;
    nuki[ 14 * kd + 82 ] += +1.000000 ;

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    nuki[ 17 * kd + 83 ] += -1.000000 ;
    nuki[ 3 * kd + 83 ] += -1.000000 ;
    nuki[ 6 * kd + 83 ] += +1.000000 ;
    nuki[ 16 * kd + 83 ] += +1.000000 ;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cv_R(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +1.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +2.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +3.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +2.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 10: CH4 */
        species[10] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 11: CO */
        species[11] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +3.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +1.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561000e-12 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +2.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +3.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            +3.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 19: N2 */
        species[19] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 20: AR */
        species[20] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +1.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +1.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +1.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +1.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 13: HCO */
        species[13] =
            +1.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +7.60690080e-01
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +2.77079900e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +1.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +9.54656420e-01
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            +7.18815000e-02
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 19: N2 */
        species[19] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 20: AR */
        species[20] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cp_R(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +3.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +4.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +3.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 10: CH4 */
        species[10] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 11: CO */
        species[11] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +4.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +2.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561000e-12 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +3.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +4.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            +4.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 19: N2 */
        species[19] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 20: AR */
        species[20] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +2.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +2.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +2.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 10: CH4 */
        species[10] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 13: HCO */
        species[13] =
            +2.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +1.76069008e+00
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +3.77079900e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +2.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.95465642e+00
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            +1.07188150e+00
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 19: N2 */
        species[19] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 20: AR */
        species[20] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
    }
    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void gibbs(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -9.179351730000000e+02 * invT
            +1.661320882000000e+00
            -2.344331120000000e+00 * tc[0]
            -3.990260375000000e-03 * tc[1]
            +3.246358500000000e-06 * tc[2]
            -1.679767450000000e-09 * tc[3]
            +3.688058805000000e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +2.547365990000000e+04 * invT
            +2.946682853000000e+00
            -2.500000000000000e+00 * tc[0]
            -3.526664095000000e-13 * tc[1]
            +3.326532733333333e-16 * tc[2]
            -1.917346933333333e-19 * tc[3]
            +4.638661660000000e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.912225920000000e+04 * invT
            +1.116333640000000e+00
            -3.168267100000000e+00 * tc[0]
            +1.639659420000000e-03 * tc[1]
            -1.107177326666667e-06 * tc[2]
            +5.106721866666666e-10 * tc[3]
            -1.056329855000000e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.615080560000000e+03 * invT
            +4.095940888000000e+00
            -3.992015430000000e+00 * tc[0]
            +1.200658760000000e-03 * tc[1]
            -7.696564016666666e-07 * tc[2]
            +3.234277775000000e-10 * tc[3]
            -6.820573500000000e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +4.600404010000000e+04 * invT
            +2.200146820000000e+00
            -3.762678670000000e+00 * tc[0]
            -4.844360715000000e-04 * tc[1]
            -4.658164016666667e-07 * tc[2]
            +3.209092941666667e-10 * tc[3]
            -8.437085950000000e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.049681630000000e+04 * invT
            +4.967723077000000e+00
            -4.198604110000000e+00 * tc[0]
            +1.183307095000000e-03 * tc[1]
            -1.372160366666667e-06 * tc[2]
            +5.573466508333334e-10 * tc[3]
            -9.715736850000000e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.644499880000000e+04 * invT
            +2.069026070000000e+00
            -3.673590400000000e+00 * tc[0]
            -1.005475875000000e-03 * tc[1]
            -9.550364266666668e-07 * tc[2]
            +5.725978541666666e-10 * tc[3]
            -1.271928670000000e-13 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -1.024664760000000e+04 * invT
            +9.791179889999999e+00
            -5.149876130000000e+00 * tc[0]
            +6.835489400000000e-03 * tc[1]
            -8.196676650000000e-06 * tc[2]
            +4.039525216666667e-09 * tc[3]
            -8.334697800000000e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.839564960000000e+03 * invT
            +8.268134100000002e-01
            -4.221185840000000e+00 * tc[0]
            +1.621962660000000e-03 * tc[1]
            -2.296657433333333e-06 * tc[2]
            +1.109534108333333e-09 * tc[3]
            -2.168844325000000e-13 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.430895670000000e+04 * invT
            +4.190910250000000e+00
            -4.793723150000000e+00 * tc[0]
            +4.954166845000000e-03 * tc[1]
            -6.220333466666666e-06 * tc[2]
            +3.160710508333333e-09 * tc[3]
            -6.588632600000000e-13 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +9.786011000000000e+02 * invT
            -1.104597300000000e+01
            -2.106204000000000e+00 * tc[0]
            -3.608297500000000e-03 * tc[1]
            -8.897453333333333e-07 * tc[2]
            +6.148030000000000e-10 * tc[3]
            -1.037805000000000e-13 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +5.089775930000000e+03 * invT
            -1.381294799999999e-01
            -3.959201480000000e+00 * tc[0]
            +3.785261235000000e-03 * tc[1]
            -9.516504866666667e-06 * tc[2]
            +5.763239608333333e-09 * tc[3]
            -1.349421865000000e-12 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.284162650000000e+04 * invT
            -4.007435600000004e-01
            -4.306465680000000e+00 * tc[0]
            +2.093294460000000e-03 * tc[1]
            -8.285713450000000e-06 * tc[2]
            +4.992721716666666e-09 * tc[3]
            -1.152545020000000e-12 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.152220550000000e+04 * invT
            +1.624601760000000e+00
            -4.291424920000000e+00 * tc[0]
            +2.750771350000000e-03 * tc[1]
            -9.990638133333334e-06 * tc[2]
            +5.903885708333334e-09 * tc[3]
            -1.343428855000000e-12 * tc[4];
        /*species 19: N2 */
        species[19] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
        /*species 20: AR */
        species[20] =
            -7.453750000000000e+02 * invT
            -1.866000000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -9.501589220000000e+02 * invT
            +6.542302510000000e+00
            -3.337279200000000e+00 * tc[0]
            +2.470123655000000e-05 * tc[1]
            -8.324279633333333e-08 * tc[2]
            +1.496386616666667e-11 * tc[3]
            -1.001276880000000e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +2.547365990000000e+04 * invT
            +2.946682924000000e+00
            -2.500000010000000e+00 * tc[0]
            +1.154214865000000e-11 * tc[1]
            -2.692699133333334e-15 * tc[2]
            +3.945960291666667e-19 * tc[3]
            -2.490986785000000e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.921757910000000e+04 * invT
            -2.214917859999999e+00
            -2.569420780000000e+00 * tc[0]
            +4.298705685000000e-05 * tc[1]
            -6.991409816666667e-09 * tc[2]
            +8.348149916666666e-13 * tc[3]
            -6.141684549999999e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.858657000000000e+03 * invT
            -1.383808430000000e+00
            -3.092887670000000e+00 * tc[0]
            -2.742148580000000e-04 * tc[1]
            -2.108420466666667e-08 * tc[2]
            +7.328846300000000e-12 * tc[3]
            -5.870618800000000e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +4.626360400000000e+04 * invT
            -3.297092110000000e+00
            -2.874101130000000e+00 * tc[0]
            -1.828196460000000e-03 * tc[1]
            +2.348243283333333e-07 * tc[2]
            -2.168162908333333e-11 * tc[3]
            +9.386378350000000e-16 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.092599970000000e+04 * invT
            -6.334463270000000e+00
            -2.292038420000000e+00 * tc[0]
            -2.327943185000000e-03 * tc[1]
            +3.353199116666667e-07 * tc[2]
            -3.482550000000000e-11 * tc[3]
            +1.698581825000000e-15 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.677558430000000e+04 * invT
            -6.194354070000000e+00
            -2.285717720000000e+00 * tc[0]
            -3.619950185000000e-03 * tc[1]
            +4.978572466666667e-07 * tc[2]
            -4.964038700000000e-11 * tc[3]
            +2.335771970000000e-15 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -9.468344590000001e+03 * invT
            -1.836246650500000e+01
            -7.485149500000000e-02 * tc[0]
            -6.695473350000000e-03 * tc[1]
            +9.554763483333333e-07 * tc[2]
            -1.019104458333333e-10 * tc[3]
            +5.090761500000000e-15 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.011918150000000e+03 * invT
            -7.026170540000000e+00
            -2.772174380000000e+00 * tc[0]
            -2.478477630000000e-03 * tc[1]
            +4.140760216666667e-07 * tc[2]
            -4.909681483333334e-11 * tc[3]
            +2.667543555000000e-15 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.399583230000000e+04 * invT
            -1.189563292000000e+01
            -1.760690080000000e+00 * tc[0]
            -4.600000410000000e-03 * tc[1]
            +7.370980216666666e-07 * tc[2]
            -8.386767666666666e-11 * tc[3]
            +4.419278200000001e-15 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +1.278325200000000e+02 * invT
            +8.412240000000000e-01
            -3.770799000000000e+00 * tc[0]
            -3.935748500000000e-03 * tc[1]
            +4.427306666666667e-07 * tc[2]
            -3.287025833333333e-11 * tc[3]
            +1.056308000000000e-15 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +4.939886140000000e+03 * invT
            -8.269258140000002e+00
            -2.036111160000000e+00 * tc[0]
            -7.322707550000000e-03 * tc[1]
            +1.118463191666667e-06 * tc[2]
            -1.226857691666667e-10 * tc[3]
            +6.285303050000000e-15 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.285752000000000e+04 * invT
            -1.150777788000000e+01
            -1.954656420000000e+00 * tc[0]
            -8.698636100000001e-03 * tc[1]
            +1.330344446666667e-06 * tc[2]
            -1.460147408333333e-10 * tc[3]
            +7.482078800000000e-15 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.142639320000000e+04 * invT
            -1.404372920000000e+01
            -1.071881500000000e+00 * tc[0]
            -1.084263385000000e-02 * tc[1]
            +1.670934450000000e-06 * tc[2]
            -1.845100008333333e-10 * tc[3]
            +9.500144500000000e-15 * tc[4];
        /*species 19: N2 */
        species[19] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 20: AR */
        species[20] =
            -7.453750000000000e+02 * invT
            -1.866000000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void helmholtz(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -9.17935173e+02 * invT
            +6.61320882e-01
            -2.34433112e+00 * tc[0]
            -3.99026037e-03 * tc[1]
            +3.24635850e-06 * tc[2]
            -1.67976745e-09 * tc[3]
            +3.68805881e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54736599e+04 * invT
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -3.52666409e-13 * tc[1]
            +3.32653273e-16 * tc[2]
            -1.91734693e-19 * tc[3]
            +4.63866166e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.91222592e+04 * invT
            +1.16333640e-01
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.61508056e+03 * invT
            +3.09594089e+00
            -3.99201543e+00 * tc[0]
            +1.20065876e-03 * tc[1]
            -7.69656402e-07 * tc[2]
            +3.23427778e-10 * tc[3]
            -6.82057350e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.02937267e+04 * invT
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.94808040e+02 * invT
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +4.60040401e+04 * invT
            +1.20014682e+00
            -3.76267867e+00 * tc[0]
            -4.84436072e-04 * tc[1]
            -4.65816402e-07 * tc[2]
            +3.20909294e-10 * tc[3]
            -8.43708595e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.04968163e+04 * invT
            +3.96772308e+00
            -4.19860411e+00 * tc[0]
            +1.18330710e-03 * tc[1]
            -1.37216037e-06 * tc[2]
            +5.57346651e-10 * tc[3]
            -9.71573685e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.64449988e+04 * invT
            +1.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -1.02466476e+04 * invT
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.43440860e+04 * invT
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.83719697e+04 * invT
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.83956496e+03 * invT
            -1.73186590e-01
            -4.22118584e+00 * tc[0]
            +1.62196266e-03 * tc[1]
            -2.29665743e-06 * tc[2]
            +1.10953411e-09 * tc[3]
            -2.16884432e-13 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.43089567e+04 * invT
            +3.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +9.78601100e+02 * invT
            -1.20459730e+01
            -2.10620400e+00 * tc[0]
            -3.60829750e-03 * tc[1]
            -8.89745333e-07 * tc[2]
            +6.14803000e-10 * tc[3]
            -1.03780500e-13 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +5.08977593e+03 * invT
            -1.13812948e+00
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.28416265e+04 * invT
            -1.40074356e+00
            -4.30646568e+00 * tc[0]
            +2.09329446e-03 * tc[1]
            -8.28571345e-06 * tc[2]
            +4.99272172e-09 * tc[3]
            -1.15254502e-12 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.15222055e+04 * invT
            +6.24601760e-01
            -4.29142492e+00 * tc[0]
            +2.75077135e-03 * tc[1]
            -9.99063813e-06 * tc[2]
            +5.90388571e-09 * tc[3]
            -1.34342886e-12 * tc[4];
        /*species 19: N2 */
        species[19] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 20: AR */
        species[20] =
            -7.45375000e+02 * invT
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -9.50158922e+02 * invT
            +5.54230251e+00
            -3.33727920e+00 * tc[0]
            +2.47012365e-05 * tc[1]
            -8.32427963e-08 * tc[2]
            +1.49638662e-11 * tc[3]
            -1.00127688e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54736599e+04 * invT
            +1.94668292e+00
            -2.50000001e+00 * tc[0]
            +1.15421486e-11 * tc[1]
            -2.69269913e-15 * tc[2]
            +3.94596029e-19 * tc[3]
            -2.49098679e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.92175791e+04 * invT
            -3.21491786e+00
            -2.56942078e+00 * tc[0]
            +4.29870569e-05 * tc[1]
            -6.99140982e-09 * tc[2]
            +8.34814992e-13 * tc[3]
            -6.14168455e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.08845772e+03 * invT
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.85865700e+03 * invT
            -2.38380843e+00
            -3.09288767e+00 * tc[0]
            -2.74214858e-04 * tc[1]
            -2.10842047e-08 * tc[2]
            +7.32884630e-12 * tc[3]
            -5.87061880e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.00042971e+04 * invT
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.11856713e+02 * invT
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +4.62636040e+04 * invT
            -4.29709211e+00
            -2.87410113e+00 * tc[0]
            -1.82819646e-03 * tc[1]
            +2.34824328e-07 * tc[2]
            -2.16816291e-11 * tc[3]
            +9.38637835e-16 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.09259997e+04 * invT
            -7.33446327e+00
            -2.29203842e+00 * tc[0]
            -2.32794318e-03 * tc[1]
            +3.35319912e-07 * tc[2]
            -3.48255000e-11 * tc[3]
            +1.69858182e-15 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.67755843e+04 * invT
            -7.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -9.46834459e+03 * invT
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.41518724e+04 * invT
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.87591660e+04 * invT
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.01191815e+03 * invT
            -8.02617054e+00
            -2.77217438e+00 * tc[0]
            -2.47847763e-03 * tc[1]
            +4.14076022e-07 * tc[2]
            -4.90968148e-11 * tc[3]
            +2.66754356e-15 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.39958323e+04 * invT
            -1.28956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +1.27832520e+02 * invT
            -1.58776000e-01
            -3.77079900e+00 * tc[0]
            -3.93574850e-03 * tc[1]
            +4.42730667e-07 * tc[2]
            -3.28702583e-11 * tc[3]
            +1.05630800e-15 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +4.93988614e+03 * invT
            -9.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.28575200e+04 * invT
            -1.25077779e+01
            -1.95465642e+00 * tc[0]
            -8.69863610e-03 * tc[1]
            +1.33034445e-06 * tc[2]
            -1.46014741e-10 * tc[3]
            +7.48207880e-15 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.14263932e+04 * invT
            -1.50437292e+01
            -1.07188150e+00 * tc[0]
            -1.08426339e-02 * tc[1]
            +1.67093445e-06 * tc[2]
            -1.84510001e-10 * tc[3]
            +9.50014450e-15 * tc[4];
        /*species 19: N2 */
        species[19] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 20: AR */
        species[20] =
            -7.45375000e+02 * invT
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void speciesInternalEnergy(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +1.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 4: OH */
        species[4] =
            +2.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: CH2 */
        species[7] =
            +2.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
        /*species 8: CH2(S) */
        species[8] =
            +3.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 * invT;
        /*species 9: CH3 */
        species[9] =
            +2.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 10: CH4 */
        species[10] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 11: CO */
        species[11] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +3.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +3.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +1.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 16: C2H4 */
        species[16] =
            +2.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 17: C2H5 */
        species[17] =
            +3.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 * invT;
        /*species 18: C2H6 */
        species[18] =
            +3.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 19: N2 */
        species[19] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 20: AR */
        species[20] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 1: H */
        species[1] =
            +1.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +1.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 4: OH */
        species[4] =
            +2.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: CH2 */
        species[7] =
            +1.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
        /*species 8: CH2(S) */
        species[8] =
            +1.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 * invT;
        /*species 9: CH3 */
        species[9] =
            +1.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 10: CH4 */
        species[10] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 11: CO */
        species[11] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +1.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +7.60690080e-01
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +2.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 * invT;
        /*species 16: C2H4 */
        species[16] =
            +1.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 17: C2H5 */
        species[17] =
            +9.54656420e-01
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 * invT;
        /*species 18: C2H6 */
        species[18] =
            +7.18815000e-02
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 19: N2 */
        species[19] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 20: AR */
        species[20] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void speciesEnthalpy(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: CH2 */
        species[7] =
            +3.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
        /*species 8: CH2(S) */
        species[8] =
            +4.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 * invT;
        /*species 9: CH3 */
        species[9] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 10: CH4 */
        species[10] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 11: CO */
        species[11] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +4.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +2.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 16: C2H4 */
        species[16] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 17: C2H5 */
        species[17] =
            +4.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 * invT;
        /*species 18: C2H6 */
        species[18] =
            +4.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 19: N2 */
        species[19] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 20: AR */
        species[20] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: CH2 */
        species[7] =
            +2.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
        /*species 8: CH2(S) */
        species[8] =
            +2.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 * invT;
        /*species 9: CH3 */
        species[9] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 10: CH4 */
        species[10] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 11: CO */
        species[11] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +2.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +3.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 * invT;
        /*species 16: C2H4 */
        species[16] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 17: C2H5 */
        species[17] =
            +1.95465642e+00
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 * invT;
        /*species 18: C2H6 */
        species[18] =
            +1.07188150e+00
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 19: N2 */
        species[19] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 20: AR */
        species[20] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void speciesEntropy(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00 * tc[0]
            +7.98052075e-03 * tc[1]
            -9.73907550e-06 * tc[2]
            +6.71906980e-09 * tc[3]
            -1.84402940e-12 * tc[4]
            +6.83010238e-01 ;
        /*species 1: H */
        species[1] =
            +2.50000000e+00 * tc[0]
            +7.05332819e-13 * tc[1]
            -9.97959820e-16 * tc[2]
            +7.66938773e-19 * tc[3]
            -2.31933083e-22 * tc[4]
            -4.46682853e-01 ;
        /*species 2: O */
        species[2] =
            +3.16826710e+00 * tc[0]
            -3.27931884e-03 * tc[1]
            +3.32153198e-06 * tc[2]
            -2.04268875e-09 * tc[3]
            +5.28164927e-13 * tc[4]
            +2.05193346e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00 * tc[0]
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 4: OH */
        species[4] =
            +3.99201543e+00 * tc[0]
            -2.40131752e-03 * tc[1]
            +2.30896920e-06 * tc[2]
            -1.29371111e-09 * tc[3]
            +3.41028675e-13 * tc[4]
            -1.03925458e-01 ;
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00 * tc[0]
            -4.74912051e-03 * tc[1]
            +1.05791445e-05 * tc[2]
            -8.09212980e-09 * tc[3]
            +2.32306281e-12 * tc[4]
            +3.71666245e+00 ;
        /*species 7: CH2 */
        species[7] =
            +3.76267867e+00 * tc[0]
            +9.68872143e-04 * tc[1]
            +1.39744921e-06 * tc[2]
            -1.28363718e-09 * tc[3]
            +4.21854298e-13 * tc[4]
            +1.56253185e+00 ;
        /*species 8: CH2(S) */
        species[8] =
            +4.19860411e+00 * tc[0]
            -2.36661419e-03 * tc[1]
            +4.11648110e-06 * tc[2]
            -2.22938660e-09 * tc[3]
            +4.85786843e-13 * tc[4]
            -7.69118967e-01 ;
        /*species 9: CH3 */
        species[9] =
            +3.67359040e+00 * tc[0]
            +2.01095175e-03 * tc[1]
            +2.86510928e-06 * tc[2]
            -2.29039142e-09 * tc[3]
            +6.35964335e-13 * tc[4]
            +1.60456433e+00 ;
        /*species 10: CH4 */
        species[10] =
            +5.14987613e+00 * tc[0]
            -1.36709788e-02 * tc[1]
            +2.45900299e-05 * tc[2]
            -1.61581009e-08 * tc[3]
            +4.16734890e-12 * tc[4]
            -4.64130376e+00 ;
        /*species 11: CO */
        species[11] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 12: CO2 */
        species[12] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
        /*species 13: HCO */
        species[13] =
            +4.22118584e+00 * tc[0]
            -3.24392532e-03 * tc[1]
            +6.88997230e-06 * tc[2]
            -4.43813643e-09 * tc[3]
            +1.08442216e-12 * tc[4]
            +3.39437243e+00 ;
        /*species 14: CH2O */
        species[14] =
            +4.79372315e+00 * tc[0]
            -9.90833369e-03 * tc[1]
            +1.86610004e-05 * tc[2]
            -1.26428420e-08 * tc[3]
            +3.29431630e-12 * tc[4]
            +6.02812900e-01 ;
        /*species 15: CH3O */
        species[15] =
            +2.10620400e+00 * tc[0]
            +7.21659500e-03 * tc[1]
            +2.66923600e-06 * tc[2]
            -2.45921200e-09 * tc[3]
            +5.18902500e-13 * tc[4]
            +1.31521770e+01 ;
        /*species 16: C2H4 */
        species[16] =
            +3.95920148e+00 * tc[0]
            -7.57052247e-03 * tc[1]
            +2.85495146e-05 * tc[2]
            -2.30529584e-08 * tc[3]
            +6.74710933e-12 * tc[4]
            +4.09733096e+00 ;
        /*species 17: C2H5 */
        species[17] =
            +4.30646568e+00 * tc[0]
            -4.18658892e-03 * tc[1]
            +2.48571403e-05 * tc[2]
            -1.99708869e-08 * tc[3]
            +5.76272510e-12 * tc[4]
            +4.70720924e+00 ;
        /*species 18: C2H6 */
        species[18] =
            +4.29142492e+00 * tc[0]
            -5.50154270e-03 * tc[1]
            +2.99719144e-05 * tc[2]
            -2.36155428e-08 * tc[3]
            +6.71714427e-12 * tc[4]
            +2.66682316e+00 ;
        /*species 19: N2 */
        species[19] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 20: AR */
        species[20] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600000e+00 ;
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00 * tc[0]
            -4.94024731e-05 * tc[1]
            +2.49728389e-07 * tc[2]
            -5.98554647e-11 * tc[3]
            +5.00638440e-15 * tc[4]
            -3.20502331e+00 ;
        /*species 1: H */
        species[1] =
            +2.50000001e+00 * tc[0]
            -2.30842973e-11 * tc[1]
            +8.07809740e-15 * tc[2]
            -1.57838412e-18 * tc[3]
            +1.24549339e-22 * tc[4]
            -4.46682914e-01 ;
        /*species 2: O */
        species[2] =
            +2.56942078e+00 * tc[0]
            -8.59741137e-05 * tc[1]
            +2.09742295e-08 * tc[2]
            -3.33925997e-12 * tc[3]
            +3.07084227e-16 * tc[4]
            +4.78433864e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 4: OH */
        species[4] =
            +3.09288767e+00 * tc[0]
            +5.48429716e-04 * tc[1]
            +6.32526140e-08 * tc[2]
            -2.93153852e-11 * tc[3]
            +2.93530940e-15 * tc[4]
            +4.47669610e+00 ;
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00 * tc[0]
            +2.23982013e-03 * tc[1]
            -3.16829075e-07 * tc[2]
            +3.80821233e-11 * tc[3]
            -2.69771337e-15 * tc[4]
            +3.78510215e+00 ;
        /*species 7: CH2 */
        species[7] =
            +2.87410113e+00 * tc[0]
            +3.65639292e-03 * tc[1]
            -7.04472985e-07 * tc[2]
            +8.67265163e-11 * tc[3]
            -4.69318918e-15 * tc[4]
            +6.17119324e+00 ;
        /*species 8: CH2(S) */
        species[8] =
            +2.29203842e+00 * tc[0]
            +4.65588637e-03 * tc[1]
            -1.00595973e-06 * tc[2]
            +1.39302000e-10 * tc[3]
            -8.49290912e-15 * tc[4]
            +8.62650169e+00 ;
        /*species 9: CH3 */
        species[9] =
            +2.28571772e+00 * tc[0]
            +7.23990037e-03 * tc[1]
            -1.49357174e-06 * tc[2]
            +1.98561548e-10 * tc[3]
            -1.16788599e-14 * tc[4]
            +8.48007179e+00 ;
        /*species 10: CH4 */
        species[10] =
            +7.48514950e-02 * tc[0]
            +1.33909467e-02 * tc[1]
            -2.86642905e-06 * tc[2]
            +4.07641783e-10 * tc[3]
            -2.54538075e-14 * tc[4]
            +1.84373180e+01 ;
        /*species 11: CO */
        species[11] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 12: CO2 */
        species[12] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
        /*species 13: HCO */
        species[13] =
            +2.77217438e+00 * tc[0]
            +4.95695526e-03 * tc[1]
            -1.24222806e-06 * tc[2]
            +1.96387259e-10 * tc[3]
            -1.33377178e-14 * tc[4]
            +9.79834492e+00 ;
        /*species 14: CH2O */
        species[14] =
            +1.76069008e+00 * tc[0]
            +9.20000082e-03 * tc[1]
            -2.21129406e-06 * tc[2]
            +3.35470707e-10 * tc[3]
            -2.20963910e-14 * tc[4]
            +1.36563230e+01 ;
        /*species 15: CH3O */
        species[15] =
            +3.77079900e+00 * tc[0]
            +7.87149700e-03 * tc[1]
            -1.32819200e-06 * tc[2]
            +1.31481033e-10 * tc[3]
            -5.28154000e-15 * tc[4]
            +2.92957500e+00 ;
        /*species 16: C2H4 */
        species[16] =
            +2.03611116e+00 * tc[0]
            +1.46454151e-02 * tc[1]
            -3.35538958e-06 * tc[2]
            +4.90743077e-10 * tc[3]
            -3.14265152e-14 * tc[4]
            +1.03053693e+01 ;
        /*species 17: C2H5 */
        species[17] =
            +1.95465642e+00 * tc[0]
            +1.73972722e-02 * tc[1]
            -3.99103334e-06 * tc[2]
            +5.84058963e-10 * tc[3]
            -3.74103940e-14 * tc[4]
            +1.34624343e+01 ;
        /*species 18: C2H6 */
        species[18] =
            +1.07188150e+00 * tc[0]
            +2.16852677e-02 * tc[1]
            -5.01280335e-06 * tc[2]
            +7.38040003e-10 * tc[3]
            -4.75007225e-14 * tc[4]
            +1.51156107e+01 ;
        /*species 19: N2 */
        species[19] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 20: AR */
        species[20] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600000e+00 ;
    }
    return;
}


/*Returns the mean specific heat at CP (Eq. 33) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCPBL(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  cpbl)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cpor[21]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 21; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCPBS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  cpbs)
{
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cpor[21], tresult[21]; /* temporary storage */
    amrex::Real imw[21];

    get_imw(imw);

    cp_R(cpor, tc);

    for (int i = 0; i < 21; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 21; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCVBL(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  cvbl)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cvor[21]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 21; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCVBS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  cvbs)
{
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cvor[21]; /* temporary storage */
    amrex::Real imw[21];

    get_imw(imw);

    cv_R(cvor, tc);

    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*H */
    result += cvor[2]*y[2]*imw[2]; /*O */
    result += cvor[3]*y[3]*imw[3]; /*O2 */
    result += cvor[4]*y[4]*imw[4]; /*OH */
    result += cvor[5]*y[5]*imw[5]; /*H2O */
    result += cvor[6]*y[6]*imw[6]; /*HO2 */
    result += cvor[7]*y[7]*imw[7]; /*CH2 */
    result += cvor[8]*y[8]*imw[8]; /*CH2(S) */
    result += cvor[9]*y[9]*imw[9]; /*CH3 */
    result += cvor[10]*y[10]*imw[10]; /*CH4 */
    result += cvor[11]*y[11]*imw[11]; /*CO */
    result += cvor[12]*y[12]*imw[12]; /*CO2 */
    result += cvor[13]*y[13]*imw[13]; /*HCO */
    result += cvor[14]*y[14]*imw[14]; /*CH2O */
    result += cvor[15]*y[15]*imw[15]; /*CH3O */
    result += cvor[16]*y[16]*imw[16]; /*C2H4 */
    result += cvor[17]*y[17]*imw[17]; /*C2H5 */
    result += cvor[18]*y[18]*imw[18]; /*C2H6 */
    result += cvor[19]*y[19]*imw[19]; /*N2 */
    result += cvor[20]*y[20]*imw[20]; /*AR */

    *cvbs = result * 8.31446261815324e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKHBML(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  hbml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real hml[21]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 21; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKHBMS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  hbms)
{
    amrex::Real result = 0;
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real hml[21], tmp[21]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[21];

    get_imw(imw);

    speciesEnthalpy(hml, tc);

    int id;
    for (id = 0; id < 21; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 21; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKUBML(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  ubml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real uml[21]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 21; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKUBMS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  ubms)
{
    amrex::Real result = 0;
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real ums[21]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[21];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*H */
    result += y[2]*ums[2]*imw[2]; /*O */
    result += y[3]*ums[3]*imw[3]; /*O2 */
    result += y[4]*ums[4]*imw[4]; /*OH */
    result += y[5]*ums[5]*imw[5]; /*H2O */
    result += y[6]*ums[6]*imw[6]; /*HO2 */
    result += y[7]*ums[7]*imw[7]; /*CH2 */
    result += y[8]*ums[8]*imw[8]; /*CH2(S) */
    result += y[9]*ums[9]*imw[9]; /*CH3 */
    result += y[10]*ums[10]*imw[10]; /*CH4 */
    result += y[11]*ums[11]*imw[11]; /*CO */
    result += y[12]*ums[12]*imw[12]; /*CO2 */
    result += y[13]*ums[13]*imw[13]; /*HCO */
    result += y[14]*ums[14]*imw[14]; /*CH2O */
    result += y[15]*ums[15]*imw[15]; /*CH3O */
    result += y[16]*ums[16]*imw[16]; /*C2H4 */
    result += y[17]*ums[17]*imw[17]; /*C2H5 */
    result += y[18]*ums[18]*imw[18]; /*C2H6 */
    result += y[19]*ums[19]*imw[19]; /*N2 */
    result += y[20]*ums[20]*imw[20]; /*AR */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKSBML(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  sbml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    amrex::Real logPratio = log ( *P / 1013250.0 ); 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real sor[21]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 21; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31446261815324e+07;
}


/*get mixture entropy in mass units */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKSBMS(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  sbms)
{
    amrex::Real result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    amrex::Real logPratio = log ( *P / 1013250.0 ); 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real sor[21]; /* temporary storage */
    amrex::Real x[21]; /* need a ytx conversion */
    amrex::Real YOW = 0; /*See Eq 4, 6 in CK Manual */
    amrex::Real imw[21];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*CH2 */
    YOW += y[8]*imw[8]; /*CH2(S) */
    YOW += y[9]*imw[9]; /*CH3 */
    YOW += y[10]*imw[10]; /*CH4 */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    YOW += y[13]*imw[13]; /*HCO */
    YOW += y[14]*imw[14]; /*CH2O */
    YOW += y[15]*imw[15]; /*CH3O */
    YOW += y[16]*imw[16]; /*C2H4 */
    YOW += y[17]*imw[17]; /*C2H5 */
    YOW += y[18]*imw[18]; /*C2H6 */
    YOW += y[19]*imw[19]; /*N2 */
    YOW += y[20]*imw[20]; /*AR */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(14.027090*YOW); 
    x[8] = y[8]/(14.027090*YOW); 
    x[9] = y[9]/(15.035060*YOW); 
    x[10] = y[10]/(16.043030*YOW); 
    x[11] = y[11]/(28.010550*YOW); 
    x[12] = y[12]/(44.009950*YOW); 
    x[13] = y[13]/(29.018520*YOW); 
    x[14] = y[14]/(30.026490*YOW); 
    x[15] = y[15]/(31.034460*YOW); 
    x[16] = y[16]/(28.054180*YOW); 
    x[17] = y[17]/(29.062150*YOW); 
    x[18] = y[18]/(30.070120*YOW); 
    x[19] = y[19]/(28.013400*YOW); 
    x[20] = y[20]/(39.948000*YOW); 
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
    result += x[13]*(sor[13]-log((x[13]+1e-100))-logPratio);
    result += x[14]*(sor[14]-log((x[14]+1e-100))-logPratio);
    result += x[15]*(sor[15]-log((x[15]+1e-100))-logPratio);
    result += x[16]*(sor[16]-log((x[16]+1e-100))-logPratio);
    result += x[17]*(sor[17]-log((x[17]+1e-100))-logPratio);
    result += x[18]*(sor[18]-log((x[18]+1e-100))-logPratio);
    result += x[19]*(sor[19]-log((x[19]+1e-100))-logPratio);
    result += x[20]*(sor[20]-log((x[20]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.31446261815324e+07 * YOW;
}

/* get temperature given internal energy in mass units and mass fracs */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void GET_T_GIVEN_EY(amrex::Real *  e, amrex::Real *  y, amrex::Real *  t, int * ierr)
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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void GET_T_GIVEN_HY(amrex::Real *  h, amrex::Real *  y, amrex::Real *  t, int * ierr)
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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKPX(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  P)
{
    amrex::Real XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
    *P = *rho * 8.31446261815324e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKPY(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW = 0;/* for computing mean MW */
    amrex::Real imw[21];

    get_imw(imw);

    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*CH2 */
    YOW += y[8]*imw[8]; /*CH2(S) */
    YOW += y[9]*imw[9]; /*CH3 */
    YOW += y[10]*imw[10]; /*CH4 */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    YOW += y[13]*imw[13]; /*HCO */
    YOW += y[14]*imw[14]; /*CH2O */
    YOW += y[15]*imw[15]; /*CH3O */
    YOW += y[16]*imw[16]; /*C2H4 */
    YOW += y[17]*imw[17]; /*C2H5 */
    YOW += y[18]*imw[18]; /*C2H6 */
    YOW += y[19]*imw[19]; /*N2 */
    YOW += y[20]*imw[20]; /*AR */
    *P = *rho * 8.31446261815324e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(c) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKPC(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  c,  amrex::Real *  P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    amrex::Real W = 0;
    amrex::Real sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*14.027090; /*CH2 */
    W += c[8]*14.027090; /*CH2(S) */
    W += c[9]*15.035060; /*CH3 */
    W += c[10]*16.043030; /*CH4 */
    W += c[11]*28.010550; /*CO */
    W += c[12]*44.009950; /*CO2 */
    W += c[13]*29.018520; /*HCO */
    W += c[14]*30.026490; /*CH2O */
    W += c[15]*31.034460; /*CH3O */
    W += c[16]*28.054180; /*C2H4 */
    W += c[17]*29.062150; /*C2H5 */
    W += c[18]*30.070120; /*C2H6 */
    W += c[19]*28.013400; /*N2 */
    W += c[20]*39.948000; /*AR */

    for (id = 0; id < 21; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31446261815324e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKRHOX(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  rho)
{
    amrex::Real XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
    *rho = *P * XW / (8.31446261815324e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKRHOY(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  rho)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[21];
    amrex::Real imw[21];

    get_imw(imw);

    for (int i = 0; i < 21; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 21; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31446261815324e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKRHOC(amrex::Real *  P, amrex::Real *  T, amrex::Real *  c,  amrex::Real *  rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    amrex::Real W = 0;
    amrex::Real sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*14.027090; /*CH2 */
    W += c[8]*14.027090; /*CH2(S) */
    W += c[9]*15.035060; /*CH3 */
    W += c[10]*16.043030; /*CH4 */
    W += c[11]*28.010550; /*CO */
    W += c[12]*44.009950; /*CO2 */
    W += c[13]*29.018520; /*HCO */
    W += c[14]*30.026490; /*CH2O */
    W += c[15]*31.034460; /*CH3O */
    W += c[16]*28.054180; /*C2H4 */
    W += c[17]*29.062150; /*C2H5 */
    W += c[18]*30.070120; /*C2H6 */
    W += c[19]*28.013400; /*N2 */
    W += c[20]*39.948000; /*AR */

    for (id = 0; id < 21; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31446261815324e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKWT( amrex::Real *  wt)
{
    get_mw(wt);
}

/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKMMWY(amrex::Real *  y,  amrex::Real *  wtm)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[21];
    amrex::Real imw[21];

    get_imw(imw);

    for (int i = 0; i < 21; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 21; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKMMWX(amrex::Real *  x,  amrex::Real *  wtm)
{
    amrex::Real XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
    *wtm = XW;

    return;
}

/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKMMWC(amrex::Real *  c,  amrex::Real *  wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    amrex::Real W = 0;
    amrex::Real sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*14.027090; /*CH2 */
    W += c[8]*14.027090; /*CH2(S) */
    W += c[9]*15.035060; /*CH3 */
    W += c[10]*16.043030; /*CH4 */
    W += c[11]*28.010550; /*CO */
    W += c[12]*44.009950; /*CO2 */
    W += c[13]*29.018520; /*HCO */
    W += c[14]*30.026490; /*CH2O */
    W += c[15]*31.034460; /*CH3O */
    W += c[16]*28.054180; /*C2H4 */
    W += c[17]*29.062150; /*C2H5 */
    W += c[18]*30.070120; /*C2H6 */
    W += c[19]*28.013400; /*N2 */
    W += c[20]*39.948000; /*AR */

    for (id = 0; id < 21; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCPOR(amrex::Real *  T, amrex::Real *  cpor)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKHORT(amrex::Real *  T, amrex::Real *  hort)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKSOR(amrex::Real *  T, amrex::Real *  sor)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}

/*convert y[species] (mass fracs) to x[species] (mole fracs) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKYTX(amrex::Real *  y,  amrex::Real *  x)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[21];
    amrex::Real imw[21];

    get_imw(imw);

    for (int i = 0; i < 21; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 21; i++)
    {
        YOW += tmp[i];
    }

    amrex::Real YOWINV = 1.0/YOW;

    for (int i = 0; i < 21; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKYTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  c)
{
    amrex::Real YOW = 0;
    amrex::Real PWORT;
    amrex::Real imw[21];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 21; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 21; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446261815324e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 21; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKYTCR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  c)
{
    amrex::Real imw[21];

    get_imw(imw);

    for (int i = 0; i < 21; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKXTY(amrex::Real *  x,  amrex::Real *  y)
{
    amrex::Real XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
    /*Now compute conversion */
    amrex::Real XWinv = 1.0/XW;
    y[0] = x[0]*2.015940*XWinv; 
    y[1] = x[1]*1.007970*XWinv; 
    y[2] = x[2]*15.999400*XWinv; 
    y[3] = x[3]*31.998800*XWinv; 
    y[4] = x[4]*17.007370*XWinv; 
    y[5] = x[5]*18.015340*XWinv; 
    y[6] = x[6]*33.006770*XWinv; 
    y[7] = x[7]*14.027090*XWinv; 
    y[8] = x[8]*14.027090*XWinv; 
    y[9] = x[9]*15.035060*XWinv; 
    y[10] = x[10]*16.043030*XWinv; 
    y[11] = x[11]*28.010550*XWinv; 
    y[12] = x[12]*44.009950*XWinv; 
    y[13] = x[13]*29.018520*XWinv; 
    y[14] = x[14]*30.026490*XWinv; 
    y[15] = x[15]*31.034460*XWinv; 
    y[16] = x[16]*28.054180*XWinv; 
    y[17] = x[17]*29.062150*XWinv; 
    y[18] = x[18]*30.070120*XWinv; 
    y[19] = x[19]*28.013400*XWinv; 
    y[20] = x[20]*39.948000*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKXTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  c)
{
    int id; /*loop counter */
    amrex::Real PORT = (*P)/(8.31446261815324e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKXTCR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  c)
{
    int id; /*loop counter */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCTX(amrex::Real *  c, amrex::Real *  x)
{
    int id; /*loop counter */
    amrex::Real sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 21; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    amrex::Real sumCinv = 1.0/sumC;
    for (id = 0; id < 21; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCTY(amrex::Real *  c, amrex::Real *  y)
{
    amrex::Real CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*2.015940; /*H2 */
    CW += c[1]*1.007970; /*H */
    CW += c[2]*15.999400; /*O */
    CW += c[3]*31.998800; /*O2 */
    CW += c[4]*17.007370; /*OH */
    CW += c[5]*18.015340; /*H2O */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*14.027090; /*CH2 */
    CW += c[8]*14.027090; /*CH2(S) */
    CW += c[9]*15.035060; /*CH3 */
    CW += c[10]*16.043030; /*CH4 */
    CW += c[11]*28.010550; /*CO */
    CW += c[12]*44.009950; /*CO2 */
    CW += c[13]*29.018520; /*HCO */
    CW += c[14]*30.026490; /*CH2O */
    CW += c[15]*31.034460; /*CH3O */
    CW += c[16]*28.054180; /*C2H4 */
    CW += c[17]*29.062150; /*C2H5 */
    CW += c[18]*30.070120; /*C2H6 */
    CW += c[19]*28.013400; /*N2 */
    CW += c[20]*39.948000; /*AR */
    /*Now compute conversion */
    amrex::Real CWinv = 1.0/CW;
    y[0] = c[0]*2.015940*CWinv; 
    y[1] = c[1]*1.007970*CWinv; 
    y[2] = c[2]*15.999400*CWinv; 
    y[3] = c[3]*31.998800*CWinv; 
    y[4] = c[4]*17.007370*CWinv; 
    y[5] = c[5]*18.015340*CWinv; 
    y[6] = c[6]*33.006770*CWinv; 
    y[7] = c[7]*14.027090*CWinv; 
    y[8] = c[8]*14.027090*CWinv; 
    y[9] = c[9]*15.035060*CWinv; 
    y[10] = c[10]*16.043030*CWinv; 
    y[11] = c[11]*28.010550*CWinv; 
    y[12] = c[12]*44.009950*CWinv; 
    y[13] = c[13]*29.018520*CWinv; 
    y[14] = c[14]*30.026490*CWinv; 
    y[15] = c[15]*31.034460*CWinv; 
    y[16] = c[16]*28.054180*CWinv; 
    y[17] = c[17]*29.062150*CWinv; 
    y[18] = c[18]*30.070120*CWinv; 
    y[19] = c[19]*28.013400*CWinv; 
    y[20] = c[20]*39.948000*CWinv; 

    return;
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCVML(amrex::Real *  T,  amrex::Real *  cvml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        cvml[id] *= 8.31446261815324e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCPML(amrex::Real *  T,  amrex::Real *  cpml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        cpml[id] *= 8.31446261815324e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKUML(amrex::Real *  T,  amrex::Real *  uml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKHML(amrex::Real *  T,  amrex::Real *  hml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        hml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKSML(amrex::Real *  T,  amrex::Real *  sml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        sml[id] *= 8.31446261815324e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCVMS(amrex::Real *  T,  amrex::Real *  cvms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 4.124360158612479e+07; /*H2 */
    cvms[1] *= 8.248720317224957e+07; /*H */
    cvms[2] *= 5.196734013871295e+06; /*O */
    cvms[3] *= 2.598367006935648e+06; /*O2 */
    cvms[4] *= 4.888740950630956e+06; /*OH */
    cvms[5] *= 4.615212712140454e+06; /*H2O */
    cvms[6] *= 2.519017346487778e+06; /*HO2 */
    cvms[7] *= 5.927432288630956e+06; /*CH2 */
    cvms[8] *= 5.927432288630956e+06; /*CH2(S) */
    cvms[9] *= 5.530049509714786e+06; /*CH3 */
    cvms[10] *= 5.182601178301878e+06; /*CH4 */
    cvms[11] *= 2.968332509769797e+06; /*CO */
    cvms[12] *= 1.889223372931176e+06; /*CO2 */
    cvms[13] *= 2.865226282440744e+06; /*HCO */
    cvms[14] *= 2.769042474879095e+06; /*CH2O */
    cvms[15] *= 2.679106586083096e+06; /*CH3O */
    cvms[16] *= 2.963716144315478e+06; /*C2H4 */
    cvms[17] *= 2.860924817383862e+06; /*C2H5 */
    cvms[18] *= 2.765024754857393e+06; /*C2H6 */
    cvms[19] *= 2.968030520448514e+06; /*N2 */
    cvms[20] *= 2.081321372322329e+06; /*AR */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCPMS(amrex::Real *  T,  amrex::Real *  cpms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124360158612479e+07; /*H2 */
    cpms[1] *= 8.248720317224957e+07; /*H */
    cpms[2] *= 5.196734013871295e+06; /*O */
    cpms[3] *= 2.598367006935648e+06; /*O2 */
    cpms[4] *= 4.888740950630956e+06; /*OH */
    cpms[5] *= 4.615212712140454e+06; /*H2O */
    cpms[6] *= 2.519017346487778e+06; /*HO2 */
    cpms[7] *= 5.927432288630956e+06; /*CH2 */
    cpms[8] *= 5.927432288630956e+06; /*CH2(S) */
    cpms[9] *= 5.530049509714786e+06; /*CH3 */
    cpms[10] *= 5.182601178301878e+06; /*CH4 */
    cpms[11] *= 2.968332509769797e+06; /*CO */
    cpms[12] *= 1.889223372931176e+06; /*CO2 */
    cpms[13] *= 2.865226282440744e+06; /*HCO */
    cpms[14] *= 2.769042474879095e+06; /*CH2O */
    cpms[15] *= 2.679106586083096e+06; /*CH3O */
    cpms[16] *= 2.963716144315478e+06; /*C2H4 */
    cpms[17] *= 2.860924817383862e+06; /*C2H5 */
    cpms[18] *= 2.765024754857393e+06; /*C2H6 */
    cpms[19] *= 2.968030520448514e+06; /*N2 */
    cpms[20] *= 2.081321372322329e+06; /*AR */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKUMS(amrex::Real *  T,  amrex::Real *  ums)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[21];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    for (int i = 0; i < 21; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKHMS(amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[21];

    get_imw(imw);

    speciesEnthalpy(hms, tc);

    for (int i = 0; i < 21; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKSMS(amrex::Real *  T,  amrex::Real *  sms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 4.124360158612479e+07; /*H2 */
    sms[1] *= 8.248720317224957e+07; /*H */
    sms[2] *= 5.196734013871295e+06; /*O */
    sms[3] *= 2.598367006935648e+06; /*O2 */
    sms[4] *= 4.888740950630956e+06; /*OH */
    sms[5] *= 4.615212712140454e+06; /*H2O */
    sms[6] *= 2.519017346487778e+06; /*HO2 */
    sms[7] *= 5.927432288630956e+06; /*CH2 */
    sms[8] *= 5.927432288630956e+06; /*CH2(S) */
    sms[9] *= 5.530049509714786e+06; /*CH3 */
    sms[10] *= 5.182601178301878e+06; /*CH4 */
    sms[11] *= 2.968332509769797e+06; /*CO */
    sms[12] *= 1.889223372931176e+06; /*CO2 */
    sms[13] *= 2.865226282440744e+06; /*HCO */
    sms[14] *= 2.769042474879095e+06; /*CH2O */
    sms[15] *= 2.679106586083096e+06; /*CH3O */
    sms[16] *= 2.963716144315478e+06; /*C2H4 */
    sms[17] *= 2.860924817383862e+06; /*C2H5 */
    sms[18] *= 2.765024754857393e+06; /*C2H6 */
    sms[19] *= 2.968030520448514e+06; /*N2 */
    sms[20] *= 2.081321372322329e+06; /*AR */
}

#ifdef AMREX_USE_GPU

/* GPU version of productionRate: no more use of thermo namespace vectors  */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void comp_qfqr(amrex::Real *  qf, amrex::Real * qr, amrex::Real * sc, amrex::Real * qss_sc, amrex::Real * tc, amrex::Real invT)
{

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    qf[0] = sc[1]*sc[7];
    qr[0] = sc[9];

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    qf[1] = sc[1]*sc[9];
    qr[1] = sc[10];

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    qf[2] = sc[1]*sc[13];
    qr[2] = sc[14];

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    qf[3] = sc[1]*sc[14];
    qr[3] = sc[15];

    /*reaction 5: H + C2H4 (+M) <=> C2H5 (+M) */
    qf[4] = sc[1]*sc[16];
    qr[4] = sc[17];

    /*reaction 6: H + C2H5 (+M) <=> C2H6 (+M) */
    qf[5] = sc[1]*sc[17];
    qr[5] = sc[18];

    /*reaction 7: H2 + CO (+M) <=> CH2O (+M) */
    qf[6] = sc[0]*sc[11];
    qr[6] = sc[14];

    /*reaction 8: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    qf[7] = pow(sc[9], 2.000000);
    qr[7] = sc[18];

    /*reaction 9: O + H + M <=> OH + M */
    qf[8] = sc[1]*sc[2];
    qr[8] = sc[4];

    /*reaction 10: O + CO + M <=> CO2 + M */
    qf[9] = sc[2]*sc[11];
    qr[9] = sc[12];

    /*reaction 11: H + O2 + M <=> HO2 + M */
    qf[10] = sc[1]*sc[3];
    qr[10] = sc[6];

    /*reaction 12: 2.000000 H + M <=> H2 + M */
    qf[11] = pow(sc[1], 2.000000);
    qr[11] = sc[0];

    /*reaction 13: H + OH + M <=> H2O + M */
    qf[12] = sc[1]*sc[4];
    qr[12] = sc[5];

    /*reaction 14: HCO + M <=> H + CO + M */
    qf[13] = sc[13];
    qr[13] = sc[1]*sc[11];

    /*reaction 15: O + H2 <=> H + OH */
    qf[14] = sc[0]*sc[2];
    qr[14] = sc[1]*sc[4];

    /*reaction 16: O + HO2 <=> OH + O2 */
    qf[15] = sc[2]*sc[6];
    qr[15] = sc[3]*sc[4];

    /*reaction 17: O + CH2 <=> H + HCO */
    qf[16] = sc[2]*sc[7];
    qr[16] = sc[1]*sc[13];

    /*reaction 18: O + CH2(S) <=> H + HCO */
    qf[17] = sc[2]*sc[8];
    qr[17] = sc[1]*sc[13];

    /*reaction 19: O + CH3 <=> H + CH2O */
    qf[18] = sc[2]*sc[9];
    qr[18] = sc[1]*sc[14];

    /*reaction 20: O + CH4 <=> OH + CH3 */
    qf[19] = sc[2]*sc[10];
    qr[19] = sc[4]*sc[9];

    /*reaction 21: O + HCO <=> OH + CO */
    qf[20] = sc[2]*sc[13];
    qr[20] = sc[4]*sc[11];

    /*reaction 22: O + HCO <=> H + CO2 */
    qf[21] = sc[2]*sc[13];
    qr[21] = sc[1]*sc[12];

    /*reaction 23: O + CH2O <=> OH + HCO */
    qf[22] = sc[2]*sc[14];
    qr[22] = sc[4]*sc[13];

    /*reaction 24: O + C2H4 <=> CH3 + HCO */
    qf[23] = sc[2]*sc[16];
    qr[23] = sc[9]*sc[13];

    /*reaction 25: O + C2H5 <=> CH3 + CH2O */
    qf[24] = sc[2]*sc[17];
    qr[24] = sc[9]*sc[14];

    /*reaction 26: O + C2H6 <=> OH + C2H5 */
    qf[25] = sc[2]*sc[18];
    qr[25] = sc[4]*sc[17];

    /*reaction 27: O2 + CO <=> O + CO2 */
    qf[26] = sc[3]*sc[11];
    qr[26] = sc[2]*sc[12];

    /*reaction 28: O2 + CH2O <=> HO2 + HCO */
    qf[27] = sc[3]*sc[14];
    qr[27] = sc[6]*sc[13];

    /*reaction 29: H + 2.000000 O2 <=> HO2 + O2 */
    qf[28] = sc[1]*pow(sc[3], 2.000000);
    qr[28] = sc[3]*sc[6];

    /*reaction 30: H + O2 + H2O <=> HO2 + H2O */
    qf[29] = sc[1]*sc[3]*sc[5];
    qr[29] = sc[5]*sc[6];

    /*reaction 31: H + O2 + N2 <=> HO2 + N2 */
    qf[30] = sc[1]*sc[3]*sc[19];
    qr[30] = sc[6]*sc[19];

    /*reaction 32: H + O2 + AR <=> HO2 + AR */
    qf[31] = sc[1]*sc[3]*sc[20];
    qr[31] = sc[6]*sc[20];

    /*reaction 33: H + O2 <=> O + OH */
    qf[32] = sc[1]*sc[3];
    qr[32] = sc[2]*sc[4];

    /*reaction 34: 2.000000 H + H2 <=> 2.000000 H2 */
    qf[33] = sc[0]*pow(sc[1], 2.000000);
    qr[33] = pow(sc[0], 2.000000);

    /*reaction 35: 2.000000 H + H2O <=> H2 + H2O */
    qf[34] = pow(sc[1], 2.000000)*sc[5];
    qr[34] = sc[0]*sc[5];

    /*reaction 36: 2.000000 H + CO2 <=> H2 + CO2 */
    qf[35] = pow(sc[1], 2.000000)*sc[12];
    qr[35] = sc[0]*sc[12];

    /*reaction 37: H + HO2 <=> O2 + H2 */
    qf[36] = sc[1]*sc[6];
    qr[36] = sc[0]*sc[3];

    /*reaction 38: H + HO2 <=> 2.000000 OH */
    qf[37] = sc[1]*sc[6];
    qr[37] = pow(sc[4], 2.000000);

    /*reaction 39: H + CH4 <=> CH3 + H2 */
    qf[38] = sc[1]*sc[10];
    qr[38] = sc[0]*sc[9];

    /*reaction 40: H + HCO <=> H2 + CO */
    qf[39] = sc[1]*sc[13];
    qr[39] = sc[0]*sc[11];

    /*reaction 41: H + CH2O <=> HCO + H2 */
    qf[40] = sc[1]*sc[14];
    qr[40] = sc[0]*sc[13];

    /*reaction 42: H + CH3O <=> OH + CH3 */
    qf[41] = sc[1]*sc[15];
    qr[41] = sc[4]*sc[9];

    /*reaction 43: H + C2H6 <=> C2H5 + H2 */
    qf[42] = sc[1]*sc[18];
    qr[42] = sc[0]*sc[17];

    /*reaction 44: OH + H2 <=> H + H2O */
    qf[43] = sc[0]*sc[4];
    qr[43] = sc[1]*sc[5];

    /*reaction 45: 2.000000 OH <=> O + H2O */
    qf[44] = pow(sc[4], 2.000000);
    qr[44] = sc[2]*sc[5];

    /*reaction 46: OH + HO2 <=> O2 + H2O */
    qf[45] = sc[4]*sc[6];
    qr[45] = sc[3]*sc[5];

    /*reaction 47: OH + CH2 <=> H + CH2O */
    qf[46] = sc[4]*sc[7];
    qr[46] = sc[1]*sc[14];

    /*reaction 48: OH + CH2(S) <=> H + CH2O */
    qf[47] = sc[4]*sc[8];
    qr[47] = sc[1]*sc[14];

    /*reaction 49: OH + CH3 <=> CH2 + H2O */
    qf[48] = sc[4]*sc[9];
    qr[48] = sc[5]*sc[7];

    /*reaction 50: OH + CH3 <=> CH2(S) + H2O */
    qf[49] = sc[4]*sc[9];
    qr[49] = sc[5]*sc[8];

    /*reaction 51: OH + CH4 <=> CH3 + H2O */
    qf[50] = sc[4]*sc[10];
    qr[50] = sc[5]*sc[9];

    /*reaction 52: OH + CO <=> H + CO2 */
    qf[51] = sc[4]*sc[11];
    qr[51] = sc[1]*sc[12];

    /*reaction 53: OH + HCO <=> H2O + CO */
    qf[52] = sc[4]*sc[13];
    qr[52] = sc[5]*sc[11];

    /*reaction 54: OH + CH2O <=> HCO + H2O */
    qf[53] = sc[4]*sc[14];
    qr[53] = sc[5]*sc[13];

    /*reaction 55: OH + C2H6 <=> C2H5 + H2O */
    qf[54] = sc[4]*sc[18];
    qr[54] = sc[5]*sc[17];

    /*reaction 56: HO2 + CH2 <=> OH + CH2O */
    qf[55] = sc[6]*sc[7];
    qr[55] = sc[4]*sc[14];

    /*reaction 57: HO2 + CH3 <=> O2 + CH4 */
    qf[56] = sc[6]*sc[9];
    qr[56] = sc[3]*sc[10];

    /*reaction 58: HO2 + CH3 <=> OH + CH3O */
    qf[57] = sc[6]*sc[9];
    qr[57] = sc[4]*sc[15];

    /*reaction 59: HO2 + CO <=> OH + CO2 */
    qf[58] = sc[6]*sc[11];
    qr[58] = sc[4]*sc[12];

    /*reaction 60: CH2 + O2 <=> OH + HCO */
    qf[59] = sc[3]*sc[7];
    qr[59] = sc[4]*sc[13];

    /*reaction 61: CH2 + H2 <=> H + CH3 */
    qf[60] = sc[0]*sc[7];
    qr[60] = sc[1]*sc[9];

    /*reaction 62: CH2 + CH3 <=> H + C2H4 */
    qf[61] = sc[7]*sc[9];
    qr[61] = sc[1]*sc[16];

    /*reaction 63: CH2 + CH4 <=> 2.000000 CH3 */
    qf[62] = sc[7]*sc[10];
    qr[62] = pow(sc[9], 2.000000);

    /*reaction 64: CH2(S) + N2 <=> CH2 + N2 */
    qf[63] = sc[8]*sc[19];
    qr[63] = sc[7]*sc[19];

    /*reaction 65: CH2(S) + AR <=> CH2 + AR */
    qf[64] = sc[8]*sc[20];
    qr[64] = sc[7]*sc[20];

    /*reaction 66: CH2(S) + O2 <=> H + OH + CO */
    qf[65] = sc[3]*sc[8];
    qr[65] = sc[1]*sc[4]*sc[11];

    /*reaction 67: CH2(S) + O2 <=> CO + H2O */
    qf[66] = sc[3]*sc[8];
    qr[66] = sc[5]*sc[11];

    /*reaction 68: CH2(S) + H2 <=> CH3 + H */
    qf[67] = sc[0]*sc[8];
    qr[67] = sc[1]*sc[9];

    /*reaction 69: CH2(S) + H2O <=> CH2 + H2O */
    qf[68] = sc[5]*sc[8];
    qr[68] = sc[5]*sc[7];

    /*reaction 70: CH2(S) + CH3 <=> H + C2H4 */
    qf[69] = sc[8]*sc[9];
    qr[69] = sc[1]*sc[16];

    /*reaction 71: CH2(S) + CH4 <=> 2.000000 CH3 */
    qf[70] = sc[8]*sc[10];
    qr[70] = pow(sc[9], 2.000000);

    /*reaction 72: CH2(S) + CO <=> CH2 + CO */
    qf[71] = sc[8]*sc[11];
    qr[71] = sc[7]*sc[11];

    /*reaction 73: CH2(S) + CO2 <=> CH2 + CO2 */
    qf[72] = sc[8]*sc[12];
    qr[72] = sc[7]*sc[12];

    /*reaction 74: CH2(S) + CO2 <=> CO + CH2O */
    qf[73] = sc[8]*sc[12];
    qr[73] = sc[11]*sc[14];

    /*reaction 75: CH3 + O2 <=> O + CH3O */
    qf[74] = sc[3]*sc[9];
    qr[74] = sc[2]*sc[15];

    /*reaction 76: CH3 + O2 <=> OH + CH2O */
    qf[75] = sc[3]*sc[9];
    qr[75] = sc[4]*sc[14];

    /*reaction 77: 2.000000 CH3 <=> H + C2H5 */
    qf[76] = pow(sc[9], 2.000000);
    qr[76] = sc[1]*sc[17];

    /*reaction 78: CH3 + HCO <=> CH4 + CO */
    qf[77] = sc[9]*sc[13];
    qr[77] = sc[10]*sc[11];

    /*reaction 79: CH3 + CH2O <=> HCO + CH4 */
    qf[78] = sc[9]*sc[14];
    qr[78] = sc[10]*sc[13];

    /*reaction 80: CH3 + C2H6 <=> C2H5 + CH4 */
    qf[79] = sc[9]*sc[18];
    qr[79] = sc[10]*sc[17];

    /*reaction 81: HCO + H2O <=> H + CO + H2O */
    qf[80] = sc[5]*sc[13];
    qr[80] = sc[1]*sc[5]*sc[11];

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    qf[81] = sc[3]*sc[13];
    qr[81] = sc[6]*sc[11];

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    qf[82] = sc[3]*sc[15];
    qr[82] = sc[6]*sc[14];

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    qf[83] = sc[3]*sc[17];
    qr[83] = sc[6]*sc[16];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 21; ++i) {
        mixture += sc[i];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[21];
    gibbs(g_RT, tc);

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 * invT;
    amrex::Real refCinv = 1 / refC;

    /* Evaluate the kfs */
    amrex::Real k_f, k_r, Corr;
    amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

    // (0):  O + H + M <=> OH + M
    k_f = 1.0000000000000002e-12 * 5e+17 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    qf[8] *= Corr * k_f;
    qr[8] *= Corr * k_f / (exp(g_RT[1] + g_RT[2] - g_RT[4]) * refCinv);
    // (1):  O + H2 <=> H + OH
    k_f = 1.0000000000000002e-06 * 50000 
               * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6290) * invT);
    Corr  = 1.0;
    qf[14] *= Corr * k_f;
    qr[14] *= Corr * k_f / exp(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4]);
    // (2):  O + HO2 <=> OH + O2
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[15] *= Corr * k_f;
    qr[15] *= Corr * k_f / exp(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]);
    // (3):  O + CH2 <=> H + HCO
    k_f = 1.0000000000000002e-06 * 80000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[16] *= Corr * k_f;
    qr[16] *= Corr * k_f / exp(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[13]);
    // (4):  O + CH2(S) <=> H + HCO
    k_f = 1.0000000000000002e-06 * 15000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[17] *= Corr * k_f;
    qr[17] *= Corr * k_f / exp(-g_RT[1] + g_RT[2] + g_RT[8] - g_RT[13]);
    // (5):  O + CH3 <=> H + CH2O
    k_f = 1.0000000000000002e-06 * 84300000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[18] *= Corr * k_f;
    qr[18] *= Corr * k_f / exp(-g_RT[1] + g_RT[2] + g_RT[9] - g_RT[14]);
    // (6):  O + CH4 <=> OH + CH3
    k_f = 1.0000000000000002e-06 * 1020000000 
               * exp(1.5 * tc[0] - 0.50321666580471969 * (8600) * invT);
    Corr  = 1.0;
    qf[19] *= Corr * k_f;
    qr[19] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[9] + g_RT[10]);
    // (7):  O + CO + M <=> CO2 + M
    k_f = 1.0000000000000002e-12 * 602000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (3000) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[3] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 3.5 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.5 - 1)*sc[20];
    qf[9] *= Corr * k_f;
    qr[9] *= Corr * k_f / (exp(g_RT[2] + g_RT[11] - g_RT[12]) * refCinv);
    // (8):  O + HCO <=> OH + CO
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[20] *= Corr * k_f;
    qr[20] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[11] + g_RT[13]);
    // (9):  O + HCO <=> H + CO2
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[21] *= Corr * k_f;
    qr[21] *= Corr * k_f / exp(-g_RT[1] + g_RT[2] - g_RT[12] + g_RT[13]);
    // (10):  O + CH2O <=> OH + HCO
    k_f = 1.0000000000000002e-06 * 39000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (3540) * invT);
    Corr  = 1.0;
    qf[22] *= Corr * k_f;
    qr[22] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[13] + g_RT[14]);
    // (11):  O + C2H4 <=> CH3 + HCO
    k_f = 1.0000000000000002e-06 * 19200000 
               * exp(1.8300000000000001 * tc[0] - 0.50321666580471969 * (220) * invT);
    Corr  = 1.0;
    qf[23] *= Corr * k_f;
    qr[23] *= Corr * k_f / exp(g_RT[2] - g_RT[9] - g_RT[13] + g_RT[16]);
    // (12):  O + C2H5 <=> CH3 + CH2O
    k_f = 1.0000000000000002e-06 * 132000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[24] *= Corr * k_f;
    qr[24] *= Corr * k_f / exp(g_RT[2] - g_RT[9] - g_RT[14] + g_RT[17]);
    // (13):  O + C2H6 <=> OH + C2H5
    k_f = 1.0000000000000002e-06 * 89800000 
               * exp(1.9199999999999999 * tc[0] - 0.50321666580471969 * (5690) * invT);
    Corr  = 1.0;
    qf[25] *= Corr * k_f;
    qr[25] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[17] + g_RT[18]);
    // (14):  O2 + CO <=> O + CO2
    k_f = 1.0000000000000002e-06 * 2500000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (47800) * invT);
    Corr  = 1.0;
    qf[26] *= Corr * k_f;
    qr[26] *= Corr * k_f / exp(-g_RT[2] + g_RT[3] + g_RT[11] - g_RT[12]);
    // (15):  O2 + CH2O <=> HO2 + HCO
    k_f = 1.0000000000000002e-06 * 100000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (40000) * invT);
    Corr  = 1.0;
    qf[27] *= Corr * k_f;
    qr[27] *= Corr * k_f / exp(g_RT[3] - g_RT[6] - g_RT[13] + g_RT[14]);
    // (16):  H + O2 + M <=> HO2 + M
    k_f = 1.0000000000000002e-12 * 2.8e+18 
               * exp(-0.85999999999999999 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 0 - 1)*sc[3] + ( 0 - 1)*sc[5] + ( 0.75 - 1)*sc[11] + ( 1.5 - 1)*sc[12] + ( 1.5 - 1)*sc[18] + ( 0 - 1)*sc[19] + ( 0 - 1)*sc[20];
    qf[10] *= Corr * k_f;
    qr[10] *= Corr * k_f / (exp(g_RT[1] + g_RT[3] - g_RT[6]) * refCinv);
    // (17):  H + 2.000000 O2 <=> HO2 + O2
    k_f = 1.0000000000000002e-12 * 3e+20 
               * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[28] *= Corr * k_f;
    qr[28] *= Corr * k_f / (exp(g_RT[1] + 2.000000*g_RT[3] - g_RT[3] - g_RT[6]) * refCinv);
    // (18):  H + O2 + H2O <=> HO2 + H2O
    k_f = 1.0000000000000002e-12 * 9.38e+18 
               * exp(-0.76000000000000001 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[29] *= Corr * k_f;
    qr[29] *= Corr * k_f / (exp(g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6]) * refCinv);
    // (19):  H + O2 + N2 <=> HO2 + N2
    k_f = 1.0000000000000002e-12 * 3.75e+20 
               * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[30] *= Corr * k_f;
    qr[30] *= Corr * k_f / (exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[19] - g_RT[19]) * refCinv);
    // (20):  H + O2 + AR <=> HO2 + AR
    k_f = 1.0000000000000002e-12 * 7e+17 
               * exp(-0.80000000000000004 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[31] *= Corr * k_f;
    qr[31] *= Corr * k_f / (exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[20] - g_RT[20]) * refCinv);
    // (21):  H + O2 <=> O + OH
    k_f = 1.0000000000000002e-06 * 83000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (14413) * invT);
    Corr  = 1.0;
    qf[32] *= Corr * k_f;
    qr[32] *= Corr * k_f / exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4]);
    // (22):  2.000000 H + M <=> H2 + M
    k_f = 1.0000000000000002e-12 * 1e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 0 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.63 - 1)*sc[20];
    qf[11] *= Corr * k_f;
    qr[11] *= Corr * k_f / (exp(-g_RT[0] + 2.000000*g_RT[1]) * refCinv);
    // (23):  2.000000 H + H2 <=> 2.000000 H2
    k_f = 1.0000000000000002e-12 * 90000000000000000 
               * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[33] *= Corr * k_f;
    qr[33] *= Corr * k_f / (exp(g_RT[0] - 2.000000*g_RT[0] + 2.000000*g_RT[1]) * refCinv);
    // (24):  2.000000 H + H2O <=> H2 + H2O
    k_f = 1.0000000000000002e-12 * 6e+19 
               * exp(-1.25 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[34] *= Corr * k_f;
    qr[34] *= Corr * k_f / (exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[5] - g_RT[5]) * refCinv);
    // (25):  2.000000 H + CO2 <=> H2 + CO2
    k_f = 1.0000000000000002e-12 * 5.5e+20 
               * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[35] *= Corr * k_f;
    qr[35] *= Corr * k_f / (exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[12] - g_RT[12]) * refCinv);
    // (26):  H + OH + M <=> H2O + M
    k_f = 1.0000000000000002e-12 * 2.2e+22 
               * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 0.72999999999999998 - 1)*sc[0] + ( 3.6499999999999999 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 3 - 1)*sc[18] + ( 0.38 - 1)*sc[20];
    qf[12] *= Corr * k_f;
    qr[12] *= Corr * k_f / (exp(g_RT[1] + g_RT[4] - g_RT[5]) * refCinv);
    // (27):  H + HO2 <=> O2 + H2
    k_f = 1.0000000000000002e-06 * 28000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (1068) * invT);
    Corr  = 1.0;
    qf[36] *= Corr * k_f;
    qr[36] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6]);
    // (28):  H + HO2 <=> 2.000000 OH
    k_f = 1.0000000000000002e-06 * 134000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (635) * invT);
    Corr  = 1.0;
    qf[37] *= Corr * k_f;
    qr[37] *= Corr * k_f / exp(g_RT[1] - 2.000000*g_RT[4] + g_RT[6]);
    // (29):  H + CH2 (+M) <=> CH3 (+M)
    k_f = 1.0000000000000002e-06 * 25000000000000000 
               * exp(-0.80000000000000004 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 3.2000000000000002e+27 
               * exp(-3.1400000000000001  * tc[0] - 0.50321666580471969  * (1230) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.68000000000000005)*exp(-tc[1] / 78) 
        + 0.68000000000000005 * exp(-tc[1]/1995)  
        + exp(-5590 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[0] *= Corr * k_f;
    qr[0] *= Corr * k_f / (exp(g_RT[1] + g_RT[7] - g_RT[9]) * refCinv);
    // (30):  H + CH3 (+M) <=> CH4 (+M)
    k_f = 1.0000000000000002e-06 * 12700000000000000 
               * exp(-0.63 * tc[0] - 0.50321666580471969 * (383) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 2.4769999999999999e+33 
               * exp(-4.7599999999999998  * tc[0] - 0.50321666580471969  * (2440) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.78300000000000003)*exp(-tc[1] / 74) 
        + 0.78300000000000003 * exp(-tc[1]/2941)  
        + exp(-6964 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[1] *= Corr * k_f;
    qr[1] *= Corr * k_f / (exp(g_RT[1] + g_RT[9] - g_RT[10]) * refCinv);
    // (31):  H + CH4 <=> CH3 + H2
    k_f = 1.0000000000000002e-06 * 660000000 
               * exp(1.6200000000000001 * tc[0] - 0.50321666580471969 * (10840) * invT);
    Corr  = 1.0;
    qf[38] *= Corr * k_f;
    qr[38] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] - g_RT[9] + g_RT[10]);
    // (32):  H + HCO (+M) <=> CH2O (+M)
    k_f = 1.0000000000000002e-06 * 1090000000000 
               * exp(0.47999999999999998 * tc[0] - 0.50321666580471969 * (-260) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 1.35e+24 
               * exp(-2.5699999999999998  * tc[0] - 0.50321666580471969  * (1425) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.78239999999999998)*exp(-tc[1] / 271) 
        + 0.78239999999999998 * exp(-tc[1]/2755)  
        + exp(-6570 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[2] *= Corr * k_f;
    qr[2] *= Corr * k_f / (exp(g_RT[1] + g_RT[13] - g_RT[14]) * refCinv);
    // (33):  H + HCO <=> H2 + CO
    k_f = 1.0000000000000002e-06 * 73400000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[39] *= Corr * k_f;
    qr[39] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] - g_RT[11] + g_RT[13]);
    // (34):  H + CH2O (+M) <=> CH3O (+M)
    k_f = 1.0000000000000002e-06 * 540000000000 
               * exp(0.45400000000000001 * tc[0] - 0.50321666580471969 * (2600) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18];
    redP = Corr / k_f * 1e-12 * 2.2e+30 
               * exp(-4.7999999999999998  * tc[0] - 0.50321666580471969  * (5560) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.75800000000000001)*exp(-tc[1] / 94) 
        + 0.75800000000000001 * exp(-tc[1]/1555)  
        + exp(-4200 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[3] *= Corr * k_f;
    qr[3] *= Corr * k_f / (exp(g_RT[1] + g_RT[14] - g_RT[15]) * refCinv);
    // (35):  H + CH2O <=> HCO + H2
    k_f = 1.0000000000000002e-06 * 23000000000 
               * exp(1.05 * tc[0] - 0.50321666580471969 * (3275) * invT);
    Corr  = 1.0;
    qf[40] *= Corr * k_f;
    qr[40] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] - g_RT[13] + g_RT[14]);
    // (36):  H + CH3O <=> OH + CH3
    k_f = 1.0000000000000002e-06 * 32000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[41] *= Corr * k_f;
    qr[41] *= Corr * k_f / exp(g_RT[1] - g_RT[4] - g_RT[9] + g_RT[15]);
    // (37):  H + C2H4 (+M) <=> C2H5 (+M)
    k_f = 1.0000000000000002e-06 * 1080000000000 
               * exp(0.45400000000000001 * tc[0] - 0.50321666580471969 * (1820) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 1.1999999999999999e+42 
               * exp(-7.6200000000000001  * tc[0] - 0.50321666580471969  * (6970) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.97529999999999994)*exp(-tc[1] / 210) 
        + 0.97529999999999994 * exp(-tc[1]/984)  
        + exp(-4374 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[4] *= Corr * k_f;
    qr[4] *= Corr * k_f / (exp(g_RT[1] + g_RT[16] - g_RT[17]) * refCinv);
    // (38):  H + C2H5 (+M) <=> C2H6 (+M)
    k_f = 1.0000000000000002e-06 * 5.21e+17 
               * exp(-0.98999999999999999 * tc[0] - 0.50321666580471969 * (1580) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 1.9900000000000001e+41 
               * exp(-7.0800000000000001  * tc[0] - 0.50321666580471969  * (6685) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.84219999999999995)*exp(-tc[1] / 125) 
        + 0.84219999999999995 * exp(-tc[1]/2219)  
        + exp(-6882 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[5] *= Corr * k_f;
    qr[5] *= Corr * k_f / (exp(g_RT[1] + g_RT[17] - g_RT[18]) * refCinv);
    // (39):  H + C2H6 <=> C2H5 + H2
    k_f = 1.0000000000000002e-06 * 115000000 
               * exp(1.8999999999999999 * tc[0] - 0.50321666580471969 * (7530) * invT);
    Corr  = 1.0;
    qf[42] *= Corr * k_f;
    qr[42] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] - g_RT[17] + g_RT[18]);
    // (40):  H2 + CO (+M) <=> CH2O (+M)
    k_f = 1.0000000000000002e-06 * 43000000 
               * exp(1.5 * tc[0] - 0.50321666580471969 * (79600) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 5.0699999999999998e+27 
               * exp(-3.4199999999999999  * tc[0] - 0.50321666580471969  * (84350) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.93200000000000005)*exp(-tc[1] / 197) 
        + 0.93200000000000005 * exp(-tc[1]/1540)  
        + exp(-10300 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[6] *= Corr * k_f;
    qr[6] *= Corr * k_f / (exp(g_RT[0] + g_RT[11] - g_RT[14]) * refCinv);
    // (41):  OH + H2 <=> H + H2O
    k_f = 1.0000000000000002e-06 * 216000000 
               * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    Corr  = 1.0;
    qf[43] *= Corr * k_f;
    qr[43] *= Corr * k_f / exp(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5]);
    // (42):  2.000000 OH <=> O + H2O
    k_f = 1.0000000000000002e-06 * 35700 
               * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (-2110) * invT);
    Corr  = 1.0;
    qf[44] *= Corr * k_f;
    qr[44] *= Corr * k_f / exp(-g_RT[2] + 2.000000*g_RT[4] - g_RT[5]);
    // (43):  OH + HO2 <=> O2 + H2O
    k_f = 1.0000000000000002e-06 * 29000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    Corr  = 1.0;
    qf[45] *= Corr * k_f;
    qr[45] *= Corr * k_f / exp(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6]);
    // (44):  OH + CH2 <=> H + CH2O
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[46] *= Corr * k_f;
    qr[46] *= Corr * k_f / exp(-g_RT[1] + g_RT[4] + g_RT[7] - g_RT[14]);
    // (45):  OH + CH2(S) <=> H + CH2O
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[47] *= Corr * k_f;
    qr[47] *= Corr * k_f / exp(-g_RT[1] + g_RT[4] + g_RT[8] - g_RT[14]);
    // (46):  OH + CH3 <=> CH2 + H2O
    k_f = 1.0000000000000002e-06 * 56000000 
               * exp(1.6000000000000001 * tc[0] - 0.50321666580471969 * (5420) * invT);
    Corr  = 1.0;
    qf[48] *= Corr * k_f;
    qr[48] *= Corr * k_f / exp(g_RT[4] - g_RT[5] - g_RT[7] + g_RT[9]);
    // (47):  OH + CH3 <=> CH2(S) + H2O
    k_f = 1.0000000000000002e-06 * 25010000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[49] *= Corr * k_f;
    qr[49] *= Corr * k_f / exp(g_RT[4] - g_RT[5] - g_RT[8] + g_RT[9]);
    // (48):  OH + CH4 <=> CH3 + H2O
    k_f = 1.0000000000000002e-06 * 100000000 
               * exp(1.6000000000000001 * tc[0] - 0.50321666580471969 * (3120) * invT);
    Corr  = 1.0;
    qf[50] *= Corr * k_f;
    qr[50] *= Corr * k_f / exp(g_RT[4] - g_RT[5] - g_RT[9] + g_RT[10]);
    // (49):  OH + CO <=> H + CO2
    k_f = 1.0000000000000002e-06 * 47600000 
               * exp(1.228 * tc[0] - 0.50321666580471969 * (70) * invT);
    Corr  = 1.0;
    qf[51] *= Corr * k_f;
    qr[51] *= Corr * k_f / exp(-g_RT[1] + g_RT[4] + g_RT[11] - g_RT[12]);
    // (50):  OH + HCO <=> H2O + CO
    k_f = 1.0000000000000002e-06 * 50000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[52] *= Corr * k_f;
    qr[52] *= Corr * k_f / exp(g_RT[4] - g_RT[5] - g_RT[11] + g_RT[13]);
    // (51):  OH + CH2O <=> HCO + H2O
    k_f = 1.0000000000000002e-06 * 3430000000 
               * exp(1.1799999999999999 * tc[0] - 0.50321666580471969 * (-447) * invT);
    Corr  = 1.0;
    qf[53] *= Corr * k_f;
    qr[53] *= Corr * k_f / exp(g_RT[4] - g_RT[5] - g_RT[13] + g_RT[14]);
    // (52):  OH + C2H6 <=> C2H5 + H2O
    k_f = 1.0000000000000002e-06 * 3540000 
               * exp(2.1200000000000001 * tc[0] - 0.50321666580471969 * (870) * invT);
    Corr  = 1.0;
    qf[54] *= Corr * k_f;
    qr[54] *= Corr * k_f / exp(g_RT[4] - g_RT[5] - g_RT[17] + g_RT[18]);
    // (53):  HO2 + CH2 <=> OH + CH2O
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[55] *= Corr * k_f;
    qr[55] *= Corr * k_f / exp(-g_RT[4] + g_RT[6] + g_RT[7] - g_RT[14]);
    // (54):  HO2 + CH3 <=> O2 + CH4
    k_f = 1.0000000000000002e-06 * 1000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[56] *= Corr * k_f;
    qr[56] *= Corr * k_f / exp(-g_RT[3] + g_RT[6] + g_RT[9] - g_RT[10]);
    // (55):  HO2 + CH3 <=> OH + CH3O
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[57] *= Corr * k_f;
    qr[57] *= Corr * k_f / exp(-g_RT[4] + g_RT[6] + g_RT[9] - g_RT[15]);
    // (56):  HO2 + CO <=> OH + CO2
    k_f = 1.0000000000000002e-06 * 150000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (23600) * invT);
    Corr  = 1.0;
    qf[58] *= Corr * k_f;
    qr[58] *= Corr * k_f / exp(-g_RT[4] + g_RT[6] + g_RT[11] - g_RT[12]);
    // (57):  CH2 + O2 <=> OH + HCO
    k_f = 1.0000000000000002e-06 * 13200000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (1500) * invT);
    Corr  = 1.0;
    qf[59] *= Corr * k_f;
    qr[59] *= Corr * k_f / exp(g_RT[3] - g_RT[4] + g_RT[7] - g_RT[13]);
    // (58):  CH2 + H2 <=> H + CH3
    k_f = 1.0000000000000002e-06 * 500000 
               * exp(2 * tc[0] - 0.50321666580471969 * (7230) * invT);
    Corr  = 1.0;
    qf[60] *= Corr * k_f;
    qr[60] *= Corr * k_f / exp(g_RT[0] - g_RT[1] + g_RT[7] - g_RT[9]);
    // (59):  CH2 + CH3 <=> H + C2H4
    k_f = 1.0000000000000002e-06 * 40000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[61] *= Corr * k_f;
    qr[61] *= Corr * k_f / exp(-g_RT[1] + g_RT[7] + g_RT[9] - g_RT[16]);
    // (60):  CH2 + CH4 <=> 2.000000 CH3
    k_f = 1.0000000000000002e-06 * 2460000 
               * exp(2 * tc[0] - 0.50321666580471969 * (8270) * invT);
    Corr  = 1.0;
    qf[62] *= Corr * k_f;
    qr[62] *= Corr * k_f / exp(g_RT[7] - 2.000000*g_RT[9] + g_RT[10]);
    // (61):  CH2(S) + N2 <=> CH2 + N2
    k_f = 1.0000000000000002e-06 * 15000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (600) * invT);
    Corr  = 1.0;
    qf[63] *= Corr * k_f;
    qr[63] *= Corr * k_f / exp(-g_RT[7] + g_RT[8] + g_RT[19] - g_RT[19]);
    // (62):  CH2(S) + AR <=> CH2 + AR
    k_f = 1.0000000000000002e-06 * 9000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (600) * invT);
    Corr  = 1.0;
    qf[64] *= Corr * k_f;
    qr[64] *= Corr * k_f / exp(-g_RT[7] + g_RT[8] + g_RT[20] - g_RT[20]);
    // (63):  CH2(S) + O2 <=> H + OH + CO
    k_f = 1.0000000000000002e-06 * 28000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[65] *= Corr * k_f;
    qr[65] *= Corr * k_f / (exp(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[8] - g_RT[11]) * refC);
    // (64):  CH2(S) + O2 <=> CO + H2O
    k_f = 1.0000000000000002e-06 * 12000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[66] *= Corr * k_f;
    qr[66] *= Corr * k_f / exp(g_RT[3] - g_RT[5] + g_RT[8] - g_RT[11]);
    // (65):  CH2(S) + H2 <=> CH3 + H
    k_f = 1.0000000000000002e-06 * 70000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[67] *= Corr * k_f;
    qr[67] *= Corr * k_f / exp(g_RT[0] - g_RT[1] + g_RT[8] - g_RT[9]);
    // (66):  CH2(S) + H2O <=> CH2 + H2O
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[68] *= Corr * k_f;
    qr[68] *= Corr * k_f / exp(g_RT[5] - g_RT[5] - g_RT[7] + g_RT[8]);
    // (67):  CH2(S) + CH3 <=> H + C2H4
    k_f = 1.0000000000000002e-06 * 12000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-570) * invT);
    Corr  = 1.0;
    qf[69] *= Corr * k_f;
    qr[69] *= Corr * k_f / exp(-g_RT[1] + g_RT[8] + g_RT[9] - g_RT[16]);
    // (68):  CH2(S) + CH4 <=> 2.000000 CH3
    k_f = 1.0000000000000002e-06 * 16000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-570) * invT);
    Corr  = 1.0;
    qf[70] *= Corr * k_f;
    qr[70] *= Corr * k_f / exp(g_RT[8] - 2.000000*g_RT[9] + g_RT[10]);
    // (69):  CH2(S) + CO <=> CH2 + CO
    k_f = 1.0000000000000002e-06 * 9000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[71] *= Corr * k_f;
    qr[71] *= Corr * k_f / exp(-g_RT[7] + g_RT[8] + g_RT[11] - g_RT[11]);
    // (70):  CH2(S) + CO2 <=> CH2 + CO2
    k_f = 1.0000000000000002e-06 * 7000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[72] *= Corr * k_f;
    qr[72] *= Corr * k_f / exp(-g_RT[7] + g_RT[8] + g_RT[12] - g_RT[12]);
    // (71):  CH2(S) + CO2 <=> CO + CH2O
    k_f = 1.0000000000000002e-06 * 14000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[73] *= Corr * k_f;
    qr[73] *= Corr * k_f / exp(g_RT[8] - g_RT[11] + g_RT[12] - g_RT[14]);
    // (72):  CH3 + O2 <=> O + CH3O
    k_f = 1.0000000000000002e-06 * 26750000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (28800) * invT);
    Corr  = 1.0;
    qf[74] *= Corr * k_f;
    qr[74] *= Corr * k_f / exp(-g_RT[2] + g_RT[3] + g_RT[9] - g_RT[15]);
    // (73):  CH3 + O2 <=> OH + CH2O
    k_f = 1.0000000000000002e-06 * 36000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (8940) * invT);
    Corr  = 1.0;
    qf[75] *= Corr * k_f;
    qr[75] *= Corr * k_f / exp(g_RT[3] - g_RT[4] + g_RT[9] - g_RT[14]);
    // (74):  2.000000 CH3 (+M) <=> C2H6 (+M)
    k_f = 1.0000000000000002e-06 * 21200000000000000 
               * exp(-0.96999999999999997 * tc[0] - 0.50321666580471969 * (620) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 1.7700000000000001e+50 
               * exp(-9.6699999999999999  * tc[0] - 0.50321666580471969  * (6220) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.53249999999999997)*exp(-tc[1] / 151) 
        + 0.53249999999999997 * exp(-tc[1]/1038)  
        + exp(-4970 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[7] *= Corr * k_f;
    qr[7] *= Corr * k_f / (exp(2.000000*g_RT[9] - g_RT[18]) * refCinv);
    // (75):  2.000000 CH3 <=> H + C2H5
    k_f = 1.0000000000000002e-06 * 4990000000000 
               * exp(0.10000000000000001 * tc[0] - 0.50321666580471969 * (10600) * invT);
    Corr  = 1.0;
    qf[76] *= Corr * k_f;
    qr[76] *= Corr * k_f / exp(-g_RT[1] + 2.000000*g_RT[9] - g_RT[17]);
    // (76):  CH3 + HCO <=> CH4 + CO
    k_f = 1.0000000000000002e-06 * 26480000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[77] *= Corr * k_f;
    qr[77] *= Corr * k_f / exp(g_RT[9] - g_RT[10] - g_RT[11] + g_RT[13]);
    // (77):  CH3 + CH2O <=> HCO + CH4
    k_f = 1.0000000000000002e-06 * 3320 
               * exp(2.8100000000000001 * tc[0] - 0.50321666580471969 * (5860) * invT);
    Corr  = 1.0;
    qf[78] *= Corr * k_f;
    qr[78] *= Corr * k_f / exp(g_RT[9] - g_RT[10] - g_RT[13] + g_RT[14]);
    // (78):  CH3 + C2H6 <=> C2H5 + CH4
    k_f = 1.0000000000000002e-06 * 6140000 
               * exp(1.74 * tc[0] - 0.50321666580471969 * (10450) * invT);
    Corr  = 1.0;
    qf[79] *= Corr * k_f;
    qr[79] *= Corr * k_f / exp(g_RT[9] - g_RT[10] - g_RT[17] + g_RT[18]);
    // (79):  HCO + H2O <=> H + CO + H2O
    k_f = 1.0000000000000002e-06 * 2.244e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    Corr  = 1.0;
    qf[80] *= Corr * k_f;
    qr[80] *= Corr * k_f / (exp(-g_RT[1] + g_RT[5] - g_RT[5] - g_RT[11] + g_RT[13]) * refC);
    // (80):  HCO + M <=> H + CO + M
    k_f = 1.0000000000000002e-06 * 1.87e+17 
               * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 0 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18];
    qf[13] *= Corr * k_f;
    qr[13] *= Corr * k_f / (exp(-g_RT[1] - g_RT[11] + g_RT[13]) * refC);
    // (81):  HCO + O2 <=> HO2 + CO
    k_f = 1.0000000000000002e-06 * 7600000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (400) * invT);
    Corr  = 1.0;
    qf[81] *= Corr * k_f;
    qr[81] *= Corr * k_f / exp(g_RT[3] - g_RT[6] - g_RT[11] + g_RT[13]);
    // (82):  CH3O + O2 <=> HO2 + CH2O
    k_f = 1.0000000000000002e-06 * 4.2799999999999999e-13 
               * exp(7.5999999999999996 * tc[0] - 0.50321666580471969 * (-3530) * invT);
    Corr  = 1.0;
    qf[82] *= Corr * k_f;
    qr[82] *= Corr * k_f / exp(g_RT[3] - g_RT[6] - g_RT[14] + g_RT[15]);
    // (83):  C2H5 + O2 <=> HO2 + C2H4
    k_f = 1.0000000000000002e-06 * 840000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (3875) * invT);
    Corr  = 1.0;
    qf[83] *= Corr * k_f;
    qr[83] *= Corr * k_f / exp(g_RT[3] - g_RT[6] - g_RT[16] + g_RT[17]);


    return;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real * wdot, amrex::Real * sc, amrex::Real T)
{
    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];

    amrex::Real qdot, q_f[84], q_r[84];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 21; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[1] -= qdot;
    wdot[7] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[1] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[1] -= qdot;
    wdot[13] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] -= qdot;
    wdot[14] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[1] -= qdot;
    wdot[16] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[1] -= qdot;
    wdot[17] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[11] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[9] -= 2.000000 * qdot;
    wdot[18] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[2] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;

    qdot = q_f[12]-q_r[12];
    wdot[1] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[1] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[14]-q_r[14];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[16]-q_r[16];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[7] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[8] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[18]-q_r[18];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[9] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[19]-q_r[19];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[21]-q_r[21];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[12] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[22]-q_r[22];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[2] -= qdot;
    wdot[9] += qdot;
    wdot[13] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[2] -= qdot;
    wdot[9] += qdot;
    wdot[14] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[27]-q_r[27];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[1] -= qdot;
    wdot[3] -= 2.000000 * qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;

    qdot = q_f[29]-q_r[29];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;

    qdot = q_f[30]-q_r[30];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[19] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[31]-q_r[31];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[20] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[32]-q_r[32];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[0] -= qdot;
    wdot[0] += 2.000000 * qdot;
    wdot[1] -= 2.000000 * qdot;

    qdot = q_f[34]-q_r[34];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
    wdot[12] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[36]-q_r[36];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[37]-q_r[37];
    wdot[1] -= qdot;
    wdot[4] += 2.000000 * qdot;
    wdot[6] -= qdot;

    qdot = q_f[38]-q_r[38];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[39]-q_r[39];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[40]-q_r[40];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[41]-q_r[41];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[42]-q_r[42];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[43]-q_r[43];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[44]-q_r[44];
    wdot[2] += qdot;
    wdot[4] -= 2.000000 * qdot;
    wdot[5] += qdot;

    qdot = q_f[45]-q_r[45];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[46]-q_r[46];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[7] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[47]-q_r[47];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[8] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[48]-q_r[48];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[49]-q_r[49];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[8] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[50]-q_r[50];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[51]-q_r[51];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[52]-q_r[52];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[53]-q_r[53];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[54]-q_r[54];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[55]-q_r[55];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[7] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[56]-q_r[56];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[57]-q_r[57];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[9] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[58]-q_r[58];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[59]-q_r[59];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[7] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[60]-q_r[60];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[61]-q_r[61];
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[9] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[62]-q_r[62];
    wdot[7] -= qdot;
    wdot[9] += 2.000000 * qdot;
    wdot[10] -= qdot;

    qdot = q_f[63]-q_r[63];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[19] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[64]-q_r[64];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[20] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[65]-q_r[65];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[66]-q_r[66];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[67]-q_r[67];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[68]-q_r[68];
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[69]-q_r[69];
    wdot[1] += qdot;
    wdot[8] -= qdot;
    wdot[9] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[70]-q_r[70];
    wdot[8] -= qdot;
    wdot[9] += 2.000000 * qdot;
    wdot[10] -= qdot;

    qdot = q_f[71]-q_r[71];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[11] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[72]-q_r[72];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[12] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[73]-q_r[73];
    wdot[8] -= qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[74]-q_r[74];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[9] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[75]-q_r[75];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[9] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[76]-q_r[76];
    wdot[1] += qdot;
    wdot[9] -= 2.000000 * qdot;
    wdot[17] += qdot;

    qdot = q_f[77]-q_r[77];
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[78]-q_r[78];
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[79]-q_r[79];
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[80]-q_r[80];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[81]-q_r[81];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[82]-q_r[82];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[83]-q_r[83];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    return;
}
#endif



/*compute the production rate for each species */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKWC(amrex::Real *  T, amrex::Real *  C,  amrex::Real *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 21; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKWYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[21]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[21];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*CH2 */
    YOW += y[8]*imw[8]; /*CH2(S) */
    YOW += y[9]*imw[9]; /*CH3 */
    YOW += y[10]*imw[10]; /*CH4 */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    YOW += y[13]*imw[13]; /*HCO */
    YOW += y[14]*imw[14]; /*CH2O */
    YOW += y[15]*imw[15]; /*CH3O */
    YOW += y[16]*imw[16]; /*C2H4 */
    YOW += y[17]*imw[17]; /*C2H5 */
    YOW += y[18]*imw[18]; /*C2H6 */
    YOW += y[19]*imw[19]; /*N2 */
    YOW += y[20]*imw[20]; /*AR */
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
    c[13] = PWORT * y[13]*imw[13]; 
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 
    c[19] = PWORT * y[19]*imw[19]; 
    c[20] = PWORT * y[20]*imw[20]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKWXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[21]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKWYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[21]; /*temporary storage */
    amrex::Real imw[21];

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
    c[13] = 1e6 * (*rho) * y[13]*imw[13]; 
    c[14] = 1e6 * (*rho) * y[14]*imw[14]; 
    c[15] = 1e6 * (*rho) * y[15]*imw[15]; 
    c[16] = 1e6 * (*rho) * y[16]*imw[16]; 
    c[17] = 1e6 * (*rho) * y[17]*imw[17]; 
    c[18] = 1e6 * (*rho) * y[18]*imw[18]; 
    c[19] = 1e6 * (*rho) * y[19]*imw[19]; 
    c[20] = 1e6 * (*rho) * y[20]*imw[20]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKWXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[21]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void dcvpRdT(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +7.98052075e-03
            -3.89563020e-05 * tc[1]
            +6.04716282e-08 * tc[2]
            -2.95044704e-11 * tc[3];
        /*species 1: H */
        species[1] =
            +7.05332819e-13
            -3.99183928e-15 * tc[1]
            +6.90244896e-18 * tc[2]
            -3.71092933e-21 * tc[3];
        /*species 2: O */
        species[2] =
            -3.27931884e-03
            +1.32861279e-05 * tc[1]
            -1.83841987e-08 * tc[2]
            +8.45063884e-12 * tc[3];
        /*species 3: O2 */
        species[3] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 4: OH */
        species[4] =
            -2.40131752e-03
            +9.23587682e-06 * tc[1]
            -1.16434000e-08 * tc[2]
            +5.45645880e-12 * tc[3];
        /*species 5: H2O */
        species[5] =
            -2.03643410e-03
            +1.30408042e-05 * tc[1]
            -1.64639119e-08 * tc[2]
            +7.08791268e-12 * tc[3];
        /*species 6: HO2 */
        species[6] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 7: CH2 */
        species[7] =
            +9.68872143e-04
            +5.58979682e-06 * tc[1]
            -1.15527346e-08 * tc[2]
            +6.74966876e-12 * tc[3];
        /*species 8: CH2(S) */
        species[8] =
            -2.36661419e-03
            +1.64659244e-05 * tc[1]
            -2.00644794e-08 * tc[2]
            +7.77258948e-12 * tc[3];
        /*species 9: CH3 */
        species[9] =
            +2.01095175e-03
            +1.14604371e-05 * tc[1]
            -2.06135228e-08 * tc[2]
            +1.01754294e-11 * tc[3];
        /*species 10: CH4 */
        species[10] =
            -1.36709788e-02
            +9.83601198e-05 * tc[1]
            -1.45422908e-07 * tc[2]
            +6.66775824e-11 * tc[3];
        /*species 11: CO */
        species[11] =
            -6.10353680e-04
            +2.03362866e-06 * tc[1]
            +2.72101765e-09 * tc[2]
            -3.61769800e-12 * tc[3];
        /*species 12: CO2 */
        species[12] =
            +8.98459677e-03
            -1.42471254e-05 * tc[1]
            +7.37757066e-09 * tc[2]
            -5.74798192e-13 * tc[3];
        /*species 13: HCO */
        species[13] =
            -3.24392532e-03
            +2.75598892e-05 * tc[1]
            -3.99432279e-08 * tc[2]
            +1.73507546e-11 * tc[3];
        /*species 14: CH2O */
        species[14] =
            -9.90833369e-03
            +7.46440016e-05 * tc[1]
            -1.13785578e-07 * tc[2]
            +5.27090608e-11 * tc[3];
        /*species 15: CH3O */
        species[15] =
            +7.21659500e-03
            +1.06769440e-05 * tc[1]
            -2.21329080e-08 * tc[2]
            +8.30244000e-12 * tc[3];
        /*species 16: C2H4 */
        species[16] =
            -7.57052247e-03
            +1.14198058e-04 * tc[1]
            -2.07476626e-07 * tc[2]
            +1.07953749e-10 * tc[3];
        /*species 17: C2H5 */
        species[17] =
            -4.18658892e-03
            +9.94285614e-05 * tc[1]
            -1.79737982e-07 * tc[2]
            +9.22036016e-11 * tc[3];
        /*species 18: C2H6 */
        species[18] =
            -5.50154270e-03
            +1.19887658e-04 * tc[1]
            -2.12539886e-07 * tc[2]
            +1.07474308e-10 * tc[3];
        /*species 19: N2 */
        species[19] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
        /*species 20: AR */
        species[20] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
    } else {
        /*species 0: H2 */
        species[0] =
            -4.94024731e-05
            +9.98913556e-07 * tc[1]
            -5.38699182e-10 * tc[2]
            +8.01021504e-14 * tc[3];
        /*species 1: H */
        species[1] =
            -2.30842973e-11
            +3.23123896e-14 * tc[1]
            -1.42054571e-17 * tc[2]
            +1.99278943e-21 * tc[3];
        /*species 2: O */
        species[2] =
            -8.59741137e-05
            +8.38969178e-08 * tc[1]
            -3.00533397e-11 * tc[2]
            +4.91334764e-15 * tc[3];
        /*species 3: O2 */
        species[3] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 4: OH */
        species[4] =
            +5.48429716e-04
            +2.53010456e-07 * tc[1]
            -2.63838467e-10 * tc[2]
            +4.69649504e-14 * tc[3];
        /*species 5: H2O */
        species[5] =
            +2.17691804e-03
            -3.28145036e-07 * tc[1]
            -2.91125961e-10 * tc[2]
            +6.72803968e-14 * tc[3];
        /*species 6: HO2 */
        species[6] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 7: CH2 */
        species[7] =
            +3.65639292e-03
            -2.81789194e-06 * tc[1]
            +7.80538647e-10 * tc[2]
            -7.50910268e-14 * tc[3];
        /*species 8: CH2(S) */
        species[8] =
            +4.65588637e-03
            -4.02383894e-06 * tc[1]
            +1.25371800e-09 * tc[2]
            -1.35886546e-13 * tc[3];
        /*species 9: CH3 */
        species[9] =
            +7.23990037e-03
            -5.97428696e-06 * tc[1]
            +1.78705393e-09 * tc[2]
            -1.86861758e-13 * tc[3];
        /*species 10: CH4 */
        species[10] =
            +1.33909467e-02
            -1.14657162e-05 * tc[1]
            +3.66877605e-09 * tc[2]
            -4.07260920e-13 * tc[3];
        /*species 11: CO */
        species[11] =
            +2.06252743e-03
            -1.99765154e-06 * tc[1]
            +6.90159024e-10 * tc[2]
            -8.14590864e-14 * tc[3];
        /*species 12: CO2 */
        species[12] =
            +4.41437026e-03
            -4.42962808e-06 * tc[1]
            +1.57047056e-09 * tc[2]
            -1.88833666e-13 * tc[3];
        /*species 13: HCO */
        species[13] =
            +4.95695526e-03
            -4.96891226e-06 * tc[1]
            +1.76748533e-09 * tc[2]
            -2.13403484e-13 * tc[3];
        /*species 14: CH2O */
        species[14] =
            +9.20000082e-03
            -8.84517626e-06 * tc[1]
            +3.01923636e-09 * tc[2]
            -3.53542256e-13 * tc[3];
        /*species 15: CH3O */
        species[15] =
            +7.87149700e-03
            -5.31276800e-06 * tc[1]
            +1.18332930e-09 * tc[2]
            -8.45046400e-14 * tc[3];
        /*species 16: C2H4 */
        species[16] =
            +1.46454151e-02
            -1.34215583e-05 * tc[1]
            +4.41668769e-09 * tc[2]
            -5.02824244e-13 * tc[3];
        /*species 17: C2H5 */
        species[17] =
            +1.73972722e-02
            -1.59641334e-05 * tc[1]
            +5.25653067e-09 * tc[2]
            -5.98566304e-13 * tc[3];
        /*species 18: C2H6 */
        species[18] =
            +2.16852677e-02
            -2.00512134e-05 * tc[1]
            +6.64236003e-09 * tc[2]
            -7.60011560e-13 * tc[3];
        /*species 19: N2 */
        species[19] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
        /*species 20: AR */
        species[20] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
    }
    return;
}

/*compute an approx to the reaction Jacobian */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void aJacobian_precond(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int HP)
{
    for (int i=0; i<484; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[21];
    for (int k=0; k<21; k++) {
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
    for (int k = 0; k < 21; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[21];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[21];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[21];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = 1.0000000000000002e-06 * 25000000000000000
                * exp(-0.80000000000000004 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.80000000000000004 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* pressure-fall-off */
    k_0 = 3.2000000000000002e+27 * exp(-3.1400000000000001 * tc[0] - 0.50321666580471969 * (1230) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -3.1400000000000001 * invT + 0.50321666580471969 * (1230) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.68000000000000005)*exp(-T/78);
    Fcent2 = 0.68000000000000005 * exp(-T/1995);
    Fcent3 = exp(-5590 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/78
        -Fcent2/1995
        + Fcent3*5590*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[9];
    Kc = refCinv * exp(g_RT[1] + g_RT[7] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[9]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[7] -= q; /* CH2 */
    wdot[9] += q; /* CH3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = dcdc_fac + k_f*sc[7];
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 6*dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac + k_f*sc[1];
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac - k_r;
    dqdc[10] = 2*dcdc_fac;
    dqdc[11] = 1.5*dcdc_fac;
    dqdc[12] = 2*dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = 3*dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = 0.69999999999999996*dcdc_fac;
    for (int k=0; k<21; k++) {
        J[22*k+1] -= dqdc[k];
        J[22*k+7] -= dqdc[k];
        J[22*k+9] += dqdc[k];
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[469] -= dqdT; /* dwdot[CH2]/dT */
    J[471] += dqdT; /* dwdot[CH3]/dT */

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = 1.0000000000000002e-06 * 12700000000000000
                * exp(-0.63 * tc[0] - 0.50321666580471969 * (383) * invT);
    dlnkfdT = -0.63 * invT + 0.50321666580471969 *  (383)  * invT2;
    /* pressure-fall-off */
    k_0 = 2.4769999999999999e+33 * exp(-4.7599999999999998 * tc[0] - 0.50321666580471969 * (2440) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -4.7599999999999998 * invT + 0.50321666580471969 * (2440) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.78300000000000003)*exp(-T/74);
    Fcent2 = 0.78300000000000003 * exp(-T/2941);
    Fcent3 = exp(-6964 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/74
        -Fcent2/2941
        + Fcent3*6964*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[10];
    Kc = refCinv * exp(g_RT[1] + g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (h_RT[10]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = dcdc_fac + k_f*sc[9];
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 6*dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac + k_f*sc[1];
    dqdc[10] = 2*dcdc_fac - k_r;
    dqdc[11] = 1.5*dcdc_fac;
    dqdc[12] = 2*dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = 3*dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = 0.69999999999999996*dcdc_fac;
    for (int k=0; k<21; k++) {
        J[22*k+1] -= dqdc[k];
        J[22*k+9] -= dqdc[k];
        J[22*k+10] += dqdc[k];
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[471] -= dqdT; /* dwdot[CH3]/dT */
    J[472] += dqdT; /* dwdot[CH4]/dT */

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = 1.0000000000000002e-06 * 1090000000000
                * exp(0.47999999999999998 * tc[0] - 0.50321666580471969 * (-260) * invT);
    dlnkfdT = 0.47999999999999998 * invT + 0.50321666580471969 *  (-260)  * invT2;
    /* pressure-fall-off */
    k_0 = 1.35e+24 * exp(-2.5699999999999998 * tc[0] - 0.50321666580471969 * (1425) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.5699999999999998 * invT + 0.50321666580471969 * (1425) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.78239999999999998)*exp(-T/271);
    Fcent2 = 0.78239999999999998 * exp(-T/2755);
    Fcent3 = exp(-6570 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/271
        -Fcent2/2755
        + Fcent3*6570*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[14];
    Kc = refCinv * exp(g_RT[1] + g_RT[13] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[14]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[13] -= q; /* HCO */
    wdot[14] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = dcdc_fac + k_f*sc[13];
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 6*dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = 2*dcdc_fac;
    dqdc[11] = 1.5*dcdc_fac;
    dqdc[12] = 2*dcdc_fac;
    dqdc[13] = dcdc_fac + k_f*sc[1];
    dqdc[14] = dcdc_fac - k_r;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = 3*dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = 0.69999999999999996*dcdc_fac;
    for (int k=0; k<21; k++) {
        J[22*k+1] -= dqdc[k];
        J[22*k+13] -= dqdc[k];
        J[22*k+14] += dqdc[k];
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[475] -= dqdT; /* dwdot[HCO]/dT */
    J[476] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18];
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = 1.0000000000000002e-06 * 540000000000
                * exp(0.45400000000000001 * tc[0] - 0.50321666580471969 * (2600) * invT);
    dlnkfdT = 0.45400000000000001 * invT + 0.50321666580471969 *  (2600)  * invT2;
    /* pressure-fall-off */
    k_0 = 2.2e+30 * exp(-4.7999999999999998 * tc[0] - 0.50321666580471969 * (5560) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -4.7999999999999998 * invT + 0.50321666580471969 * (5560) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.75800000000000001)*exp(-T/94);
    Fcent2 = 0.75800000000000001 * exp(-T/1555);
    Fcent3 = exp(-4200 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/94
        -Fcent2/1555
        + Fcent3*4200*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[15];
    Kc = refCinv * exp(g_RT[1] + g_RT[14] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[15]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[14] -= q; /* CH2O */
    wdot[15] += q; /* CH3O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = dcdc_fac + k_f*sc[14];
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 6*dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = 2*dcdc_fac;
    dqdc[11] = 1.5*dcdc_fac;
    dqdc[12] = 2*dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac + k_f*sc[1];
    dqdc[15] = dcdc_fac - k_r;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = 3*dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    for (int k=0; k<21; k++) {
        J[22*k+1] -= dqdc[k];
        J[22*k+14] -= dqdc[k];
        J[22*k+15] += dqdc[k];
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[476] -= dqdT; /* dwdot[CH2O]/dT */
    J[477] += dqdT; /* dwdot[CH3O]/dT */

    /*reaction 5: H + C2H4 (+M) <=> C2H5 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = 1.0000000000000002e-06 * 1080000000000
                * exp(0.45400000000000001 * tc[0] - 0.50321666580471969 * (1820) * invT);
    dlnkfdT = 0.45400000000000001 * invT + 0.50321666580471969 *  (1820)  * invT2;
    /* pressure-fall-off */
    k_0 = 1.1999999999999999e+42 * exp(-7.6200000000000001 * tc[0] - 0.50321666580471969 * (6970) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -7.6200000000000001 * invT + 0.50321666580471969 * (6970) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.97529999999999994)*exp(-T/210);
    Fcent2 = 0.97529999999999994 * exp(-T/984);
    Fcent3 = exp(-4374 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/210
        -Fcent2/984
        + Fcent3*4374*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[17];
    Kc = refCinv * exp(g_RT[1] + g_RT[16] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[16]) + (h_RT[17]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[16] -= q; /* C2H4 */
    wdot[17] += q; /* C2H5 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = dcdc_fac + k_f*sc[16];
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 6*dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = 2*dcdc_fac;
    dqdc[11] = 1.5*dcdc_fac;
    dqdc[12] = 2*dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac + k_f*sc[1];
    dqdc[17] = dcdc_fac - k_r;
    dqdc[18] = 3*dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = 0.69999999999999996*dcdc_fac;
    for (int k=0; k<21; k++) {
        J[22*k+1] -= dqdc[k];
        J[22*k+16] -= dqdc[k];
        J[22*k+17] += dqdc[k];
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[478] -= dqdT; /* dwdot[C2H4]/dT */
    J[479] += dqdT; /* dwdot[C2H5]/dT */

    /*reaction 6: H + C2H5 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[17];
    k_f = 1.0000000000000002e-06 * 5.21e+17
                * exp(-0.98999999999999999 * tc[0] - 0.50321666580471969 * (1580) * invT);
    dlnkfdT = -0.98999999999999999 * invT + 0.50321666580471969 *  (1580)  * invT2;
    /* pressure-fall-off */
    k_0 = 1.9900000000000001e+41 * exp(-7.0800000000000001 * tc[0] - 0.50321666580471969 * (6685) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -7.0800000000000001 * invT + 0.50321666580471969 * (6685) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.84219999999999995)*exp(-T/125);
    Fcent2 = 0.84219999999999995 * exp(-T/2219);
    Fcent3 = exp(-6882 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/125
        -Fcent2/2219
        + Fcent3*6882*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[18];
    Kc = refCinv * exp(g_RT[1] + g_RT[17] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[17]) + (h_RT[18]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[17] -= q; /* C2H5 */
    wdot[18] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = dcdc_fac + k_f*sc[17];
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 6*dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = 2*dcdc_fac;
    dqdc[11] = 1.5*dcdc_fac;
    dqdc[12] = 2*dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac + k_f*sc[1];
    dqdc[18] = 3*dcdc_fac - k_r;
    dqdc[19] = dcdc_fac;
    dqdc[20] = 0.69999999999999996*dcdc_fac;
    for (int k=0; k<21; k++) {
        J[22*k+1] -= dqdc[k];
        J[22*k+17] -= dqdc[k];
        J[22*k+18] += dqdc[k];
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[479] -= dqdT; /* dwdot[C2H5]/dT */
    J[480] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 7: H2 + CO (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[0]*sc[11];
    k_f = 1.0000000000000002e-06 * 43000000
                * exp(1.5 * tc[0] - 0.50321666580471969 * (79600) * invT);
    dlnkfdT = 1.5 * invT + 0.50321666580471969 *  (79600)  * invT2;
    /* pressure-fall-off */
    k_0 = 5.0699999999999998e+27 * exp(-3.4199999999999999 * tc[0] - 0.50321666580471969 * (84350) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -3.4199999999999999 * invT + 0.50321666580471969 * (84350) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.93200000000000005)*exp(-T/197);
    Fcent2 = 0.93200000000000005 * exp(-T/1540);
    Fcent3 = exp(-10300 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/197
        -Fcent2/1540
        + Fcent3*10300*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[14];
    Kc = refCinv * exp(g_RT[0] + g_RT[11] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[11]) + (h_RT[14]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[11] -= q; /* CO */
    wdot[14] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac + k_f*sc[11];
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 6*dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = 2*dcdc_fac;
    dqdc[11] = 1.5*dcdc_fac + k_f*sc[0];
    dqdc[12] = 2*dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac - k_r;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = 3*dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = 0.69999999999999996*dcdc_fac;
    for (int k=0; k<21; k++) {
        J[22*k+0] -= dqdc[k];
        J[22*k+11] -= dqdc[k];
        J[22*k+14] += dqdc[k];
    }
    J[462] -= dqdT; /* dwdot[H2]/dT */
    J[473] -= dqdT; /* dwdot[CO]/dT */
    J[476] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 8: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = pow(sc[9], 2.000000);
    k_f = 1.0000000000000002e-06 * 21200000000000000
                * exp(-0.96999999999999997 * tc[0] - 0.50321666580471969 * (620) * invT);
    dlnkfdT = -0.96999999999999997 * invT + 0.50321666580471969 *  (620)  * invT2;
    /* pressure-fall-off */
    k_0 = 1.7700000000000001e+50 * exp(-9.6699999999999999 * tc[0] - 0.50321666580471969 * (6220) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -9.6699999999999999 * invT + 0.50321666580471969 * (6220) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.53249999999999997)*exp(-T/151);
    Fcent2 = 0.53249999999999997 * exp(-T/1038);
    Fcent3 = exp(-4970 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/151
        -Fcent2/1038
        + Fcent3*4970*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[18];
    Kc = refCinv * exp(2.000000*g_RT[9] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[9]) + (h_RT[18]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[9] -= 2 * q; /* CH3 */
    wdot[18] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = 6*dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac + k_f*2.000000*sc[9];
    dqdc[10] = 2*dcdc_fac;
    dqdc[11] = 1.5*dcdc_fac;
    dqdc[12] = 2*dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = 3*dcdc_fac - k_r;
    dqdc[19] = dcdc_fac;
    dqdc[20] = 0.69999999999999996*dcdc_fac;
    for (int k=0; k<21; k++) {
        J[22*k+9] += -2 * dqdc[k];
        J[22*k+18] += dqdc[k];
    }
    J[471] += -2 * dqdT; /* dwdot[CH3]/dT */
    J[480] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 9: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-12 * 5e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[4]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2*q_nocor;
    dqdc[1] = q_nocor + k_f*sc[2];
    dqdc[2] = q_nocor + k_f*sc[1];
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor - k_r;
    dqdc[5] = 6*q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = 2*q_nocor;
    dqdc[11] = 1.5*q_nocor;
    dqdc[12] = 2*q_nocor;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = 3*q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = 0.69999999999999996*q_nocor;
    for (int k=0; k<21; k++) {
        J[22*k+1] -= dqdc[k];
        J[22*k+2] -= dqdc[k];
        J[22*k+4] += dqdc[k];
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[464] -= dqdT; /* dwdot[O]/dT */
    J[466] += dqdT; /* dwdot[OH]/dT */

    /*reaction 10: O + CO + M <=> CO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[3] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 3.5 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.5 - 1)*sc[20];
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = 1.0000000000000002e-12 * 602000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (3000)  * invT2;
    /* reverse */
    phi_r = sc[12];
    Kc = refCinv * exp(g_RT[2] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[11]) + (h_RT[12]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2*q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor + k_f*sc[11];
    dqdc[3] = 6*q_nocor;
    dqdc[4] = q_nocor;
    dqdc[5] = 6*q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = 2*q_nocor;
    dqdc[11] = 1.5*q_nocor + k_f*sc[2];
    dqdc[12] = 3.5*q_nocor - k_r;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = 3*q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = 0.5*q_nocor;
    for (int k=0; k<21; k++) {
        J[22*k+2] -= dqdc[k];
        J[22*k+11] -= dqdc[k];
        J[22*k+12] += dqdc[k];
    }
    J[464] -= dqdT; /* dwdot[O]/dT */
    J[473] -= dqdT; /* dwdot[CO]/dT */
    J[474] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 11: H + O2 + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[3] + ( 0 - 1)*sc[5] + ( 0.75 - 1)*sc[11] + ( 1.5 - 1)*sc[12] + ( 1.5 - 1)*sc[18] + ( 0 - 1)*sc[19] + ( 0 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-12 * 2.8e+18
                * exp(-0.85999999999999999 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.85999999999999999 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = q_nocor;
    dqdc[1] = q_nocor + k_f*sc[3];
    dqdc[2] = q_nocor;
    dqdc[3] =  + k_f*sc[1];
    dqdc[4] = q_nocor;
    dqdc[5] = 0.0;
    dqdc[6] = q_nocor - k_r;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = 0.75*q_nocor;
    dqdc[12] = 1.5*q_nocor;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = 1.5*q_nocor;
    dqdc[19] = 0.0;
    dqdc[20] = 0.0;
    for (int k=0; k<21; k++) {
        J[22*k+1] -= dqdc[k];
        J[22*k+3] -= dqdc[k];
        J[22*k+6] += dqdc[k];
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[465] -= dqdT; /* dwdot[O2]/dT */
    J[468] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 12: 2.000000 H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 0 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.63 - 1)*sc[20];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 1e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[0];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1]) + (h_RT[0]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] =  - k_r;
    dqdc[1] = q_nocor + k_f*2.000000*sc[1];
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor;
    dqdc[5] = 0.0;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = 2*q_nocor;
    dqdc[11] = q_nocor;
    dqdc[12] = 0.0;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = 3*q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = 0.63*q_nocor;
    for (int k=0; k<21; k++) {
        J[22*k+0] += dqdc[k];
        J[22*k+1] += -2 * dqdc[k];
    }
    J[462] += dqdT; /* dwdot[H2]/dT */
    J[463] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 13: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0.72999999999999998 - 1)*sc[0] + ( 3.6499999999999999 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 3 - 1)*sc[18] + ( 0.38 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = 1.0000000000000002e-12 * 2.2e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[4]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 0.72999999999999998*q_nocor;
    dqdc[1] = q_nocor + k_f*sc[4];
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor + k_f*sc[1];
    dqdc[5] = 3.6499999999999999*q_nocor - k_r;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = 2*q_nocor;
    dqdc[11] = q_nocor;
    dqdc[12] = q_nocor;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = 3*q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = 0.38*q_nocor;
    for (int k=0; k<21; k++) {
        J[22*k+1] -= dqdc[k];
        J[22*k+4] -= dqdc[k];
        J[22*k+5] += dqdc[k];
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[466] -= dqdT; /* dwdot[OH]/dT */
    J[467] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 14: HCO + M <=> H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 0 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18];
    /* forward */
    phi_f = sc[13];
    k_f = 1.0000000000000002e-06 * 1.87e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (17000)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[11];
    Kc = refC * exp(-g_RT[1] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13]) + (h_RT[1] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2*q_nocor;
    dqdc[1] = q_nocor - k_r*sc[11];
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor;
    dqdc[5] = 0.0;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = 2*q_nocor;
    dqdc[11] = 1.5*q_nocor - k_r*sc[1];
    dqdc[12] = 2*q_nocor;
    dqdc[13] = q_nocor + k_f;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = 3*q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = q_nocor;
    for (int k=0; k<21; k++) {
        J[22*k+1] += dqdc[k];
        J[22*k+11] += dqdc[k];
        J[22*k+13] -= dqdc[k];
    }
    J[463] += dqdT; /* dwdot[H]/dT */
    J[473] += dqdT; /* dwdot[CO]/dT */
    J[475] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 15: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = 1.0000000000000002e-06 * 50000
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6290) * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  (6290)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[1] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[2] -= dqdci;                /* dwdot[O]/d[H2] */
    J[4] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[44] -= dqdci;               /* dwdot[H2]/d[O] */
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[88] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 16: O + HO2 <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[47] += dqdci;               /* dwdot[O2]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[50] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[68] -= dqdci;               /* dwdot[O]/d[O2] */
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[91] += dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[134] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 17: O + CH2 <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 80000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[7] -= q; /* CH2 */
    wdot[13] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[35] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[51] -= dqdci;               /* dwdot[CH2]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[2];
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[156] -= dqdci;              /* dwdot[O]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[167] += dqdci;              /* dwdot[HCO]/d[CH2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[293] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 18: O + CH2(S) <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = 1.0000000000000002e-06 * 15000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[8] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[8] -= q; /* CH2(S) */
    wdot[13] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[35] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[8];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[52] -= dqdci;               /* dwdot[CH2(S)]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[2];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[178] -= dqdci;              /* dwdot[O]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[189] += dqdci;              /* dwdot[HCO]/d[CH2(S)] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[294] -= dqdci;              /* dwdot[CH2(S)]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 19: O + CH3 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[9];
    k_f = 1.0000000000000002e-06 * 84300000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[9] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[9]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[9] -= q; /* CH3 */
    wdot[14] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[9];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[58] += dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[2];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[212] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[309] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[310] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[317] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 20: O + CH4 <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[10];
    k_f = 1.0000000000000002e-06 * 1020000000
                * exp(1.5 * tc[0] - 0.50321666580471969 * (8600) * invT);
    dlnkfdT = 1.5 * invT + 0.50321666580471969 *  (8600)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[9];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[10]) + (h_RT[4] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[9] += q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[53] += dqdci;               /* dwdot[CH3]/d[O] */
    J[54] -= dqdci;               /* dwdot[CH4]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[98] -= dqdci;               /* dwdot[CH4]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[222] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[224] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 21: O + HCO <=> OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[13];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[11];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[4] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[55] += dqdci;               /* dwdot[CO]/d[O] */
    J[57] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[11];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[99] += dqdci;               /* dwdot[CO]/d[OH] */
    J[101] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[4];
    J[244] -= dqdci;              /* dwdot[O]/d[CO] */
    J[246] += dqdci;              /* dwdot[OH]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[290] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 22: O + HCO <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[13];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[12] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[12] += q; /* CO2 */
    wdot[13] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[34] += dqdci;               /* dwdot[CO2]/d[H] */
    J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[56] += dqdci;               /* dwdot[CO2]/d[O] */
    J[57] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[265] += dqdci;              /* dwdot[H]/d[CO2] */
    J[266] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[277] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[298] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 23: O + CH2O <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[14];
    k_f = 1.0000000000000002e-06 * 39000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3540) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (3540)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[13];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[14]) + (h_RT[4] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[14];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    J[58] -= dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[101] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[102] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[290] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[2];
    J[310] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[312] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 24: O + C2H4 <=> CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[16];
    k_f = 1.0000000000000002e-06 * 19200000
                * exp(1.8300000000000001 * tc[0] - 0.50321666580471969 * (220) * invT);
    dlnkfdT = 1.8300000000000001 * invT + 0.50321666580471969 *  (220)  * invT2;
    /* reverse */
    phi_r = sc[9]*sc[13];
    Kc = exp(g_RT[2] - g_RT[9] - g_RT[13] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[16]) + (h_RT[9] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[9] += q; /* CH3 */
    wdot[13] += q; /* HCO */
    wdot[16] -= q; /* C2H4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] += dqdci;               /* dwdot[CH3]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    J[60] -= dqdci;               /* dwdot[C2H4]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[13];
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[211] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[214] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[9];
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[295] += dqdci;              /* dwdot[CH3]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[302] -= dqdci;              /* dwdot[C2H4]/d[HCO] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[2];
    J[354] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[361] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[365] += dqdci;              /* dwdot[HCO]/d[C2H4] */
    J[368] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[478] -= dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 25: O + C2H5 <=> CH3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = 1.0000000000000002e-06 * 132000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[9]*sc[14];
    Kc = exp(g_RT[2] - g_RT[9] - g_RT[14] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[9] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[9] += q; /* CH3 */
    wdot[14] += q; /* CH2O */
    wdot[17] -= q; /* C2H5 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] += dqdci;               /* dwdot[CH3]/d[O] */
    J[58] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[61] -= dqdci;               /* dwdot[C2H5]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[14];
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[212] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[215] -= dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[9];
    J[310] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[317] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[325] -= dqdci;              /* dwdot[C2H5]/d[CH2O] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[2];
    J[376] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[383] += dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[388] += dqdci;              /* dwdot[CH2O]/d[C2H5] */
    J[391] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */
    J[479] -= dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 26: O + C2H6 <=> OH + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = 1.0000000000000002e-06 * 89800000
                * exp(1.9199999999999999 * tc[0] - 0.50321666580471969 * (5690) * invT);
    dlnkfdT = 1.9199999999999999 * invT + 0.50321666580471969 *  (5690)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[17];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[18]) + (h_RT[4] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[18];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[61] += dqdci;               /* dwdot[C2H5]/d[O] */
    J[62] -= dqdci;               /* dwdot[C2H6]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[17];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[106] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[4];
    J[376] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[378] += dqdci;              /* dwdot[OH]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[2];
    J[398] -= dqdci;              /* dwdot[O]/d[C2H6] */
    J[400] += dqdci;              /* dwdot[OH]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 27: O2 + CO <=> O + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = 1.0000000000000002e-06 * 2500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (47800) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (47800)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[12];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[11]) + (h_RT[2] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[12];
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[47] -= dqdci;               /* dwdot[O2]/d[O] */
    J[55] -= dqdci;               /* dwdot[CO]/d[O] */
    J[56] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[68] += dqdci;               /* dwdot[O]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[77] -= dqdci;               /* dwdot[CO]/d[O2] */
    J[78] += dqdci;               /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3];
    J[244] += dqdci;              /* dwdot[O]/d[CO] */
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[2];
    J[266] += dqdci;              /* dwdot[O]/d[CO2] */
    J[267] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[473] -= dqdT;               /* dwdot[CO]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 28: O2 + CH2O <=> HO2 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[14];
    k_f = 1.0000000000000002e-06 * 100000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (40000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (40000)  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[13];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[14]) + (h_RT[6] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[79] += dqdci;               /* dwdot[HCO]/d[O2] */
    J[80] -= dqdci;               /* dwdot[CH2O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[13];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[145] += dqdci;              /* dwdot[HCO]/d[HO2] */
    J[146] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[6];
    J[289] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[292] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[3];
    J[311] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[314] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 29: H + 2.000000 O2 <=> HO2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*pow(sc[3], 2.000000);
    k_f = 1.0000000000000002e-12 * 3e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[6];
    Kc = refCinv * exp(g_RT[1] + 2.000000*g_RT[3] - g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + 2.000000*h_RT[3]) + (h_RT[3] + h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*pow(sc[3], 2.000000);
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*2.000000*sc[3] - k_r*sc[6];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 30: H + O2 + H2O <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = 1.0000000000000002e-12 * 9.38e+18
                * exp(-0.76000000000000001 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.76000000000000001 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[5]) + (h_RT[5] + h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[5];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[5];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[113] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[116] += dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 31: H + O2 + N2 <=> HO2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[19];
    k_f = 1.0000000000000002e-12 * 3.75e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[19];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[19] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[19]) + (h_RT[6] + h_RT[19]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[19];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[19];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[19];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[419] -= dqdci;              /* dwdot[H]/d[N2] */
    J[421] -= dqdci;              /* dwdot[O2]/d[N2] */
    J[424] += dqdci;              /* dwdot[HO2]/d[N2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 32: H + O2 + AR <=> HO2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[20];
    k_f = 1.0000000000000002e-12 * 7e+17
                * exp(-0.80000000000000004 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.80000000000000004 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[20];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[20] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[20]) + (h_RT[6] + h_RT[20]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[20];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[20];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[20];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[441] -= dqdci;              /* dwdot[H]/d[AR] */
    J[443] -= dqdci;              /* dwdot[O2]/d[AR] */
    J[446] += dqdci;              /* dwdot[HO2]/d[AR] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 33: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 83000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (14413) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (14413)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[24] += dqdci;               /* dwdot[O]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[45] -= dqdci;               /* dwdot[H]/d[O] */
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[47] -= dqdci;               /* dwdot[O2]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[68] += dqdci;               /* dwdot[O]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[89] -= dqdci;               /* dwdot[H]/d[OH] */
    J[90] += dqdci;               /* dwdot[O]/d[OH] */
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 34: 2.000000 H + H2 <=> 2.000000 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 90000000000000000
                * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.59999999999999998 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000);
    Kc = refCinv * exp(g_RT[0] - 2.000000*g_RT[0] + 2.000000*g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + 2.000000*h_RT[1]) + (2.000000*h_RT[0]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*2.000000*sc[0];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[0]*2.000000*sc[1];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 35: 2.000000 H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[5];
    k_f = 1.0000000000000002e-12 * 6e+19
                * exp(-1.25 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.25 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[5];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[5]) + (h_RT[0] + h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[5];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[110] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[111] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 36: 2.000000 H + CO2 <=> H2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[12];
    k_f = 1.0000000000000002e-12 * 5.5e+20
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[12];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[12] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[12]) + (h_RT[0] + h_RT[12]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[12];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[12];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[CO2] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[264] += dqdci;              /* dwdot[H2]/d[CO2] */
    J[265] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 37: H + HO2 <=> O2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 28000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (1068) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (1068)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[3] += dqdci;                /* dwdot[O2]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] += dqdci;               /* dwdot[O2]/d[H] */
    J[28] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[66] += dqdci;               /* dwdot[H2]/d[O2] */
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[132] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 38: H + HO2 <=> 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 134000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (635) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (635)  * invT2;
    /* reverse */
    phi_r = pow(sc[4], 2.000000);
    Kc = exp(g_RT[1] - 2.000000*g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (2.000000*h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[26] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[28] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[4];
    J[89] -= dqdci;               /* dwdot[H]/d[OH] */
    J[92] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[136] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[466] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 39: H + CH4 <=> CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = 1.0000000000000002e-06 * 660000000
                * exp(1.6200000000000001 * tc[0] - 0.50321666580471969 * (10840) * invT);
    dlnkfdT = 1.6200000000000001 * invT + 0.50321666580471969 *  (10840)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[9];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[0] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[9] += q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[9];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
    J[10] -= dqdci;               /* dwdot[CH4]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    J[32] -= dqdci;               /* dwdot[CH4]/d[H] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[0];
    J[198] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[1];
    J[220] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 40: H + HCO <=> H2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = 1.0000000000000002e-06 * 73400000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[11];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[0] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[11];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[11] += dqdci;               /* dwdot[CO]/d[H2] */
    J[13] -= dqdci;               /* dwdot[HCO]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[13];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[33] += dqdci;               /* dwdot[CO]/d[H] */
    J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[0];
    J[242] += dqdci;              /* dwdot[H2]/d[CO] */
    J[243] -= dqdci;              /* dwdot[H]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[286] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[287] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 41: H + CH2O <=> HCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = 1.0000000000000002e-06 * 23000000000
                * exp(1.05 * tc[0] - 0.50321666580471969 * (3275) * invT);
    dlnkfdT = 1.05 * invT + 0.50321666580471969 *  (3275)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[13];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[0] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[13];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[13] += dqdci;               /* dwdot[HCO]/d[H2] */
    J[14] -= dqdci;               /* dwdot[CH2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[14];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += dqdci;               /* dwdot[HCO]/d[H] */
    J[36] -= dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[0];
    J[286] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[287] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[1];
    J[308] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[309] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 42: H + CH3O <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[15];
    k_f = 1.0000000000000002e-06 * 32000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[9];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[9] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[15]) + (h_RT[4] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += q; /* OH */
    wdot[9] += q; /* CH3 */
    wdot[15] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[15];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    J[37] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[89] -= dqdci;               /* dwdot[H]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[103] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[213] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[1];
    J[331] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[334] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[339] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[345] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[477] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 43: H + C2H6 <=> C2H5 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[18];
    k_f = 1.0000000000000002e-06 * 115000000
                * exp(1.8999999999999999 * tc[0] - 0.50321666580471969 * (7530) * invT);
    dlnkfdT = 1.8999999999999999 * invT + 0.50321666580471969 *  (7530)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[17];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[18]) + (h_RT[0] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[17];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[17] += dqdci;               /* dwdot[C2H5]/d[H2] */
    J[18] -= dqdci;               /* dwdot[C2H6]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[39] += dqdci;               /* dwdot[C2H5]/d[H] */
    J[40] -= dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[0];
    J[374] += dqdci;              /* dwdot[H2]/d[C2H5] */
    J[375] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[1];
    J[396] += dqdci;              /* dwdot[H2]/d[C2H6] */
    J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 44: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  (3430)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
    J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[27] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[88] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[110] -= dqdci;              /* dwdot[H2]/d[H2O] */
    J[111] += dqdci;              /* dwdot[H]/d[H2O] */
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 45: 2.000000 OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-06 * 35700
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (-2110) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  (-2110)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + 2.000000*g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[4] -= 2 * q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[48] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[49] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[4];
    J[90] += dqdci;               /* dwdot[O]/d[OH] */
    J[92] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[112] += dqdci;              /* dwdot[O]/d[H2O] */
    J[114] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[466] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 46: OH + HO2 <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = 1.0000000000000002e-06 * 29000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-500)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[5];
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[70] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[91] += dqdci;               /* dwdot[O2]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[113] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[116] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[136] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[137] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 47: OH + CH2 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[7] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[7] -= q; /* CH2 */
    wdot[14] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[95] -= dqdci;               /* dwdot[CH2]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[158] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[168] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[309] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[312] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[315] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 48: OH + CH2(S) <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[8];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[8] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[8]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[8] -= q; /* CH2(S) */
    wdot[14] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[96] -= dqdci;               /* dwdot[CH2(S)]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[4];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[180] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[190] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[309] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[312] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[316] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 49: OH + CH3 <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = 1.0000000000000002e-06 * 56000000
                * exp(1.6000000000000001 * tc[0] - 0.50321666580471969 * (5420) * invT);
    dlnkfdT = 1.6000000000000001 * invT + 0.50321666580471969 *  (5420)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[7] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[7] += q; /* CH2 */
    wdot[9] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[95] += dqdci;               /* dwdot[CH2]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[7];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[117] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[119] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[5];
    J[158] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[159] += dqdci;              /* dwdot[H2O]/d[CH2] */
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[202] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[203] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[205] += dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 50: OH + CH3 <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = 1.0000000000000002e-06 * 25010000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[8];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[8] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[5] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[8] += q; /* CH2(S) */
    wdot[9] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[96] += dqdci;               /* dwdot[CH2(S)]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[118] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[119] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[5];
    J[180] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[181] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[184] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[202] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[203] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[206] += dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[470] += dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 51: OH + CH4 <=> CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = 1.0000000000000002e-06 * 100000000
                * exp(1.6000000000000001 * tc[0] - 0.50321666580471969 * (3120) * invT);
    dlnkfdT = 1.6000000000000001 * invT + 0.50321666580471969 *  (3120)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[9] += q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[97] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[98] -= dqdci;               /* dwdot[CH4]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[9];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[119] += dqdci;              /* dwdot[CH3]/d[H2O] */
    J[120] -= dqdci;              /* dwdot[CH4]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[5];
    J[202] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[203] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[4];
    J[224] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[225] += dqdci;              /* dwdot[H2O]/d[CH4] */
    J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 52: OH + CO <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = 1.0000000000000002e-06 * 47600000
                * exp(1.228 * tc[0] - 0.50321666580471969 * (70) * invT);
    dlnkfdT = 1.228 * invT + 0.50321666580471969 *  (70)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[33] -= dqdci;               /* dwdot[CO]/d[H] */
    J[34] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[99] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[100] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[4];
    J[243] += dqdci;              /* dwdot[H]/d[CO] */
    J[246] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[265] += dqdci;              /* dwdot[H]/d[CO2] */
    J[268] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[473] -= dqdT;               /* dwdot[CO]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 53: OH + HCO <=> H2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[13];
    k_f = 1.0000000000000002e-06 * 50000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[11];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[5] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[99] += dqdci;               /* dwdot[CO]/d[OH] */
    J[101] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[246] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[247] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[4];
    J[290] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[291] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 54: OH + CH2O <=> HCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[14];
    k_f = 1.0000000000000002e-06 * 3430000000
                * exp(1.1799999999999999 * tc[0] - 0.50321666580471969 * (-447) * invT);
    dlnkfdT = 1.1799999999999999 * invT + 0.50321666580471969 *  (-447)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[13];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[14]) + (h_RT[5] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[14];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[101] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[102] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[13];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[123] += dqdci;              /* dwdot[HCO]/d[H2O] */
    J[124] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[5];
    J[290] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[291] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[4];
    J[312] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[313] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 55: OH + C2H6 <=> C2H5 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[18];
    k_f = 1.0000000000000002e-06 * 3540000
                * exp(2.1200000000000001 * tc[0] - 0.50321666580471969 * (870) * invT);
    dlnkfdT = 2.1200000000000001 * invT + 0.50321666580471969 *  (870)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[17];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[18]) + (h_RT[5] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[18];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[105] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[106] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[17];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[127] += dqdci;              /* dwdot[C2H5]/d[H2O] */
    J[128] -= dqdci;              /* dwdot[C2H6]/d[H2O] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[5];
    J[378] -= dqdci;              /* dwdot[OH]/d[C2H5] */
    J[379] += dqdci;              /* dwdot[H2O]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[4];
    J[400] -= dqdci;              /* dwdot[OH]/d[C2H6] */
    J[401] += dqdci;              /* dwdot[H2O]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 56: HO2 + CH2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[7];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[7] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[7]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[7] -= q; /* CH2 */
    wdot[14] += q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[95] -= dqdci;               /* dwdot[CH2]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[7];
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[139] -= dqdci;              /* dwdot[CH2]/d[HO2] */
    J[146] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[6];
    J[158] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[160] -= dqdci;              /* dwdot[HO2]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[168] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[312] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[314] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[315] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 57: HO2 + CH3 <=> O2 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[9];
    k_f = 1.0000000000000002e-06 * 1000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[10];
    Kc = exp(-g_RT[3] + g_RT[6] + g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[3] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[10];
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    J[75] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[76] += dqdci;               /* dwdot[CH4]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[9];
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[141] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[142] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[201] += dqdci;              /* dwdot[O2]/d[CH3] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[3];
    J[223] += dqdci;              /* dwdot[O2]/d[CH4] */
    J[226] -= dqdci;              /* dwdot[HO2]/d[CH4] */
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */

    /*reaction 58: HO2 + CH3 <=> OH + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[9];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[15];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[9] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[4] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[9] -= q; /* CH3 */
    wdot[15] += q; /* CH3O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    J[103] += dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[9];
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[141] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[147] += dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[213] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[4];
    J[334] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[336] -= dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[339] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[345] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[477] += dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 59: HO2 + CO <=> OH + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[11];
    k_f = 1.0000000000000002e-06 * 150000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (23600) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (23600)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[12];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[11]) + (h_RT[4] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[12];
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[99] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[100] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[143] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[144] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[246] += dqdci;              /* dwdot[OH]/d[CO] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[4];
    J[268] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[270] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[473] -= dqdT;               /* dwdot[CO]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 60: CH2 + O2 <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 13200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (1500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (1500)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[13];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[7] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[7] -= q; /* CH2 */
    wdot[13] += q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[7];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[73] -= dqdci;               /* dwdot[CH2]/d[O2] */
    J[79] += dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[95] -= dqdci;               /* dwdot[CH2]/d[OH] */
    J[101] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[3];
    J[157] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[158] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[167] += dqdci;              /* dwdot[HCO]/d[CH2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[289] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[290] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[293] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 61: CH2 + H2 <=> H + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = 1.0000000000000002e-06 * 500000
                * exp(2 * tc[0] - 0.50321666580471969 * (7230) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (7230)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[7] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* CH2 */
    wdot[9] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[7] -= dqdci;                /* dwdot[CH2]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[0];
    J[154] -= dqdci;              /* dwdot[H2]/d[CH2] */
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] += dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[198] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 62: CH2 + CH3 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[9];
    k_f = 1.0000000000000002e-06 * 40000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[7] + g_RT[9] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[9]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* CH2 */
    wdot[9] -= q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[38] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[9];
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    J[170] += dqdci;              /* dwdot[C2H4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[214] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[353] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[359] -= dqdci;              /* dwdot[CH2]/d[C2H4] */
    J[361] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[368] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[478] += dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 63: CH2 + CH4 <=> 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[10];
    k_f = 1.0000000000000002e-06 * 2460000
                * exp(2 * tc[0] - 0.50321666580471969 * (8270) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (8270)  * invT2;
    /* reverse */
    phi_r = pow(sc[9], 2.000000);
    Kc = exp(g_RT[7] - 2.000000*g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[10]) + (2.000000*h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* CH2 */
    wdot[9] += 2 * q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[10];
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] += 2 * dqdci;          /* dwdot[CH3]/d[CH2] */
    J[164] -= dqdci;              /* dwdot[CH4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[9];
    J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[7];
    J[227] -= dqdci;              /* dwdot[CH2]/d[CH4] */
    J[229] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[471] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 64: CH2(S) + N2 <=> CH2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[19];
    k_f = 1.0000000000000002e-06 * 15000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (600) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (600)  * invT2;
    /* reverse */
    phi_r = sc[7]*sc[19];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[19] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[19]) + (h_RT[7] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[19];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[19];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[425] += dqdci;              /* dwdot[CH2]/d[N2] */
    J[426] -= dqdci;              /* dwdot[CH2(S)]/d[N2] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 65: CH2(S) + AR <=> CH2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[20];
    k_f = 1.0000000000000002e-06 * 9000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (600) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (600)  * invT2;
    /* reverse */
    phi_r = sc[7]*sc[20];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[20] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[20]) + (h_RT[7] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[20];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[20];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[447] += dqdci;              /* dwdot[CH2]/d[AR] */
    J[448] -= dqdci;              /* dwdot[CH2(S)]/d[AR] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 66: CH2(S) + O2 <=> H + OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = 1.0000000000000002e-06 * 28000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4]*sc[11];
    Kc = refC * exp(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[8] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[1] + h_RT[4] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[8] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4]*sc[11];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[33] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[67] += dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[74] -= dqdci;               /* dwdot[CH2(S)]/d[O2] */
    J[77] += dqdci;               /* dwdot[CO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1]*sc[11];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[96] -= dqdci;               /* dwdot[CH2(S)]/d[OH] */
    J[99] += dqdci;               /* dwdot[CO]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[179] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    J[180] += dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[187] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[4];
    J[243] += dqdci;              /* dwdot[H]/d[CO] */
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[246] += dqdci;              /* dwdot[OH]/d[CO] */
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 67: CH2(S) + O2 <=> CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = 1.0000000000000002e-06 * 12000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[11];
    Kc = exp(g_RT[3] - g_RT[5] + g_RT[8] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[5] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[5] += q; /* H2O */
    wdot[8] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[74] -= dqdci;               /* dwdot[CH2(S)]/d[O2] */
    J[77] += dqdci;               /* dwdot[CO]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[113] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[118] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[179] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    J[181] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[187] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[247] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 68: CH2(S) + H2 <=> CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[8];
    k_f = 1.0000000000000002e-06 * 70000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[8]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[8] -= q; /* CH2(S) */
    wdot[9] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[8];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[8] -= dqdci;                /* dwdot[CH2(S)]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[0];
    J[176] -= dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] += dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[198] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[206] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 69: CH2(S) + H2O <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[8];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(g_RT[5] - g_RT[5] - g_RT[7] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[8]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[117] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[118] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[5];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[5];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 70: CH2(S) + CH3 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[9];
    k_f = 1.0000000000000002e-06 * 12000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-570) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-570)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[8] + g_RT[9] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[9]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[8] -= q; /* CH2(S) */
    wdot[9] -= q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[38] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[9];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    J[192] += dqdci;              /* dwdot[C2H4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[206] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[214] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[353] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[360] -= dqdci;              /* dwdot[CH2(S)]/d[C2H4] */
    J[361] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[368] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[478] += dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 71: CH2(S) + CH4 <=> 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = 1.0000000000000002e-06 * 16000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-570) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-570)  * invT2;
    /* reverse */
    phi_r = pow(sc[9], 2.000000);
    Kc = exp(g_RT[8] - 2.000000*g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[10]) + (2.000000*h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= q; /* CH2(S) */
    wdot[9] += 2 * q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[10];
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] += 2 * dqdci;          /* dwdot[CH3]/d[CH2(S)] */
    J[186] -= dqdci;              /* dwdot[CH4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[9];
    J[206] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[8];
    J[228] -= dqdci;              /* dwdot[CH2(S)]/d[CH4] */
    J[229] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 72: CH2(S) + CO <=> CH2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = 1.0000000000000002e-06 * 9000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[7]*sc[11];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[11] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[11]) + (h_RT[7] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[11];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[11];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[249] += dqdci;              /* dwdot[CH2]/d[CO] */
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 73: CH2(S) + CO2 <=> CH2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = 1.0000000000000002e-06 * 7000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[7]*sc[12];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[12] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[7] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[12];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[12];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[271] += dqdci;              /* dwdot[CH2]/d[CO2] */
    J[272] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 74: CH2(S) + CO2 <=> CO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = 1.0000000000000002e-06 * 14000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[11]*sc[14];
    Kc = exp(g_RT[8] - g_RT[11] + g_RT[12] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[11] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    wdot[12] -= q; /* CO2 */
    wdot[14] += q; /* CH2O */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[12];
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[187] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[188] -= dqdci;              /* dwdot[CO2]/d[CH2(S)] */
    J[190] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[14];
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[254] -= dqdci;              /* dwdot[CO2]/d[CO] */
    J[256] += dqdci;              /* dwdot[CH2O]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[8];
    J[272] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    J[275] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    J[278] += dqdci;              /* dwdot[CH2O]/d[CO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[11];
    J[316] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[319] += dqdci;              /* dwdot[CO]/d[CH2O] */
    J[320] -= dqdci;              /* dwdot[CO2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[474] -= dqdT;               /* dwdot[CO2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 75: CH3 + O2 <=> O + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = 1.0000000000000002e-06 * 26750000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (28800) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (28800)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[15];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[9] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (h_RT[2] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[9] -= q; /* CH3 */
    wdot[15] += q; /* CH3O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[15];
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[47] -= dqdci;               /* dwdot[O2]/d[O] */
    J[53] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[59] += dqdci;               /* dwdot[CH3O]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[68] += dqdci;               /* dwdot[O]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[75] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[81] += dqdci;               /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[200] += dqdci;              /* dwdot[O]/d[CH3] */
    J[201] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[213] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[2];
    J[332] += dqdci;              /* dwdot[O]/d[CH3O] */
    J[333] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[339] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[345] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[477] += dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 76: CH3 + O2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = 1.0000000000000002e-06 * 36000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (8940) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (8940)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[9] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[9] -= q; /* CH3 */
    wdot[14] += q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[80] += dqdci;               /* dwdot[CH2O]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[201] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[212] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[311] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[312] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[317] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 77: 2.000000 CH3 <=> H + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[9], 2.000000);
    k_f = 1.0000000000000002e-06 * 4990000000000
                * exp(0.10000000000000001 * tc[0] - 0.50321666580471969 * (10600) * invT);
    dlnkfdT = 0.10000000000000001 * invT + 0.50321666580471969 *  (10600)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[17];
    Kc = exp(-g_RT[1] + 2.000000*g_RT[9] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[9]) + (h_RT[1] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[9] -= 2 * q; /* CH3 */
    wdot[17] += q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[17];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[31] += -2 * dqdci;          /* dwdot[CH3]/d[H] */
    J[39] += dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*2.000000*sc[9];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[207] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[215] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[1];
    J[375] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[383] += -2 * dqdci;         /* dwdot[CH3]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[471] += -2 * dqdT;          /* dwdot[CH3]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 78: CH3 + HCO <=> CH4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[13];
    k_f = 1.0000000000000002e-06 * 26480000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[10]*sc[11];
    Kc = exp(g_RT[9] - g_RT[10] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[13]) + (h_RT[10] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[13];
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[209] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[211] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[11];
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[231] += dqdci;              /* dwdot[CO]/d[CH4] */
    J[233] -= dqdci;              /* dwdot[HCO]/d[CH4] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[10];
    J[251] -= dqdci;              /* dwdot[CH3]/d[CO] */
    J[252] += dqdci;              /* dwdot[CH4]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[9];
    J[295] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[296] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 79: CH3 + CH2O <=> HCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[14];
    k_f = 1.0000000000000002e-06 * 3320
                * exp(2.8100000000000001 * tc[0] - 0.50321666580471969 * (5860) * invT);
    dlnkfdT = 2.8100000000000001 * invT + 0.50321666580471969 *  (5860)  * invT2;
    /* reverse */
    phi_r = sc[10]*sc[13];
    Kc = exp(g_RT[9] - g_RT[10] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[14]) + (h_RT[10] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[14];
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[211] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[212] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[13];
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[233] += dqdci;              /* dwdot[HCO]/d[CH4] */
    J[234] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[10];
    J[295] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[296] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[9];
    J[317] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[318] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 80: CH3 + C2H6 <=> C2H5 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[18];
    k_f = 1.0000000000000002e-06 * 6140000
                * exp(1.74 * tc[0] - 0.50321666580471969 * (10450) * invT);
    dlnkfdT = 1.74 * invT + 0.50321666580471969 *  (10450)  * invT2;
    /* reverse */
    phi_r = sc[10]*sc[17];
    Kc = exp(g_RT[9] - g_RT[10] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[18]) + (h_RT[10] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[18];
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[215] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[216] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[17];
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[237] += dqdci;              /* dwdot[C2H5]/d[CH4] */
    J[238] -= dqdci;              /* dwdot[C2H6]/d[CH4] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[10];
    J[383] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[384] += dqdci;              /* dwdot[CH4]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[9];
    J[405] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[406] += dqdci;              /* dwdot[CH4]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 81: HCO + H2O <=> H + CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = 1.0000000000000002e-06 * 2.244e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (17000)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5]*sc[11];
    Kc = refC * exp(-g_RT[1] + g_RT[5] - g_RT[5] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[13]) + (h_RT[1] + h_RT[5] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5]*sc[11];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[33] += dqdci;               /* dwdot[CO]/d[H] */
    J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[13] - k_r*sc[1]*sc[11];
    J[111] += dqdci;              /* dwdot[H]/d[H2O] */
    J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[5];
    J[243] += dqdci;              /* dwdot[H]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[13];
    k_f = 1.0000000000000002e-06 * 7600000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (400) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (400)  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[11];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[13]) + (h_RT[6] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[CO]/d[O2] */
    J[79] -= dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[11];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[143] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[145] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[248] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[289] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[292] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = 1.0000000000000002e-06 * 4.2799999999999999e-13
                * exp(7.5999999999999996 * tc[0] - 0.50321666580471969 * (-3530) * invT);
    dlnkfdT = 7.5999999999999996 * invT + 0.50321666580471969 *  (-3530)  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[14];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[14] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[15]) + (h_RT[6] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[14] += q; /* CH2O */
    wdot[15] -= q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[15];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[80] += dqdci;               /* dwdot[CH2O]/d[O2] */
    J[81] -= dqdci;               /* dwdot[CH3O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[14];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[146] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[147] -= dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[6];
    J[311] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[314] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[323] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[3];
    J[333] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[336] += dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[344] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[345] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */
    J[477] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[17];
    k_f = 1.0000000000000002e-06 * 840000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3875) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (3875)  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[16];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[17]) + (h_RT[6] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[16] += q; /* C2H4 */
    wdot[17] -= q; /* C2H5 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[17];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[82] += dqdci;               /* dwdot[C2H4]/d[O2] */
    J[83] -= dqdci;               /* dwdot[C2H5]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[16];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[148] += dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[149] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[6];
    J[355] -= dqdci;              /* dwdot[O2]/d[C2H4] */
    J[358] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[368] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[369] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[3];
    J[377] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[380] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[390] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[391] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[478] += dqdT;               /* dwdot[C2H4]/dT */
    J[479] -= dqdT;               /* dwdot[C2H5]/dT */

    amrex::Real c_R[21], dcRdT[21], e_RT[21];
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
    for (int k = 0; k < 21; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[462+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 21; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 21; ++m) {
            dehmixdc += eh_RT[m]*J[k*22+m];
        }
        J[k*22+21] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[483] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

/*compute an approx to the reaction Jacobian (for preconditioning) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void DWDOT_SIMPLIFIED(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * HP)
{
    amrex::Real c[21];

    for (int k=0; k<21; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<21; k++) {
        J[462+k] *= 1.e-6;
        J[k*22+21] *= 1.e6;
    }

    return;
}

#ifdef AMREX_USE_GPU
/*compute the reaction Jacobian on GPU */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void aJacobian(amrex::Real * J, amrex::Real * sc, amrex::Real T, int consP)
{

    for (int i=0; i<484; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[21];
    for (int k=0; k<21; k++) {
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
    for (int k = 0; k < 21; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[21];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[21];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[21];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = 1.0000000000000002e-06 * 25000000000000000
                * exp(-0.80000000000000004 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.80000000000000004 * invT + 0.50321666580471969 *  0  * invT2;
    /* pressure-fall-off */
    k_0 = 3.2000000000000002e+27 * exp(-3.1400000000000001 * tc[0] - 0.50321666580471969 * (1230) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -3.1400000000000001 * invT + 0.50321666580471969 * (1230) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.68000000000000005)*exp(-T/78);
    Fcent2 = 0.68000000000000005 * exp(-T/1995);
    Fcent3 = exp(-5590 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/78
        -Fcent2/1995
        + Fcent3*5590*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[9];
    Kc = refCinv * exp(g_RT[1] + g_RT[7] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[9]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[7] -= q; /* CH2 */
    wdot[9] += q; /* CH3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[7] -= dqdci;                /* dwdot[CH2]/d[H2] */
        J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[7];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
        J[31] += dqdci;               /* dwdot[CH3]/d[H] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[117] -= dqdci;              /* dwdot[CH2]/d[H2O] */
        J[119] += dqdci;              /* dwdot[CH3]/d[H2O] */
        /* d()/d[CH2] */
        dqdci =  + k_f*sc[1];
        J[155] -= dqdci;              /* dwdot[H]/d[CH2] */
        J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
        J[163] += dqdci;              /* dwdot[CH3]/d[CH2] */
        /* d()/d[CH3] */
        dqdci =  - k_r;
        J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
        J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[227] -= dqdci;              /* dwdot[CH2]/d[CH4] */
        J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[249] -= dqdci;              /* dwdot[CH2]/d[CO] */
        J[251] += dqdci;              /* dwdot[CH3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[271] -= dqdci;              /* dwdot[CH2]/d[CO2] */
        J[273] += dqdci;              /* dwdot[CH3]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[403] -= dqdci;              /* dwdot[CH2]/d[C2H6] */
        J[405] += dqdci;              /* dwdot[CH3]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[447] -= dqdci;              /* dwdot[CH2]/d[AR] */
        J[449] += dqdci;              /* dwdot[CH3]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[7];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac + k_f*sc[1];
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac - k_r;
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+7] -= dqdc[k];
            J[22*k+9] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[469] -= dqdT; /* dwdot[CH2]/dT */
    J[471] += dqdT; /* dwdot[CH3]/dT */

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = 1.0000000000000002e-06 * 12700000000000000
                * exp(-0.63 * tc[0] - 0.50321666580471969 * (383) * invT);
    dlnkfdT = -0.63 * invT + 0.50321666580471969 *  383  * invT2;
    /* pressure-fall-off */
    k_0 = 2.4769999999999999e+33 * exp(-4.7599999999999998 * tc[0] - 0.50321666580471969 * (2440) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -4.7599999999999998 * invT + 0.50321666580471969 * (2440) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.78300000000000003)*exp(-T/74);
    Fcent2 = 0.78300000000000003 * exp(-T/2941);
    Fcent3 = exp(-6964 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/74
        -Fcent2/2941
        + Fcent3*6964*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[10];
    Kc = refCinv * exp(g_RT[1] + g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (h_RT[10]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[9] -= dqdci;                /* dwdot[CH3]/d[H2] */
        J[10] += dqdci;               /* dwdot[CH4]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[9];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
        J[32] += dqdci;               /* dwdot[CH4]/d[H] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[119] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[120] += dqdci;              /* dwdot[CH4]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[1];
        J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac - k_r;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
        J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[251] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[252] += dqdci;              /* dwdot[CH4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[273] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[274] += dqdci;              /* dwdot[CH4]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[405] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
        J[406] += dqdci;              /* dwdot[CH4]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[449] -= dqdci;              /* dwdot[CH3]/d[AR] */
        J[450] += dqdci;              /* dwdot[CH4]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[9];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac + k_f*sc[1];
        dqdc[10] = 2*dcdc_fac - k_r;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+9] -= dqdc[k];
            J[22*k+10] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[471] -= dqdT; /* dwdot[CH3]/dT */
    J[472] += dqdT; /* dwdot[CH4]/dT */

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = 1.0000000000000002e-06 * 1090000000000
                * exp(0.47999999999999998 * tc[0] - 0.50321666580471969 * (-260) * invT);
    dlnkfdT = 0.47999999999999998 * invT + 0.50321666580471969 *  -260  * invT2;
    /* pressure-fall-off */
    k_0 = 1.35e+24 * exp(-2.5699999999999998 * tc[0] - 0.50321666580471969 * (1425) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.5699999999999998 * invT + 0.50321666580471969 * (1425) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.78239999999999998)*exp(-T/271);
    Fcent2 = 0.78239999999999998 * exp(-T/2755);
    Fcent3 = exp(-6570 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/271
        -Fcent2/2755
        + Fcent3*6570*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[14];
    Kc = refCinv * exp(g_RT[1] + g_RT[13] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[14]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[13] -= q; /* HCO */
    wdot[14] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[13] -= dqdci;               /* dwdot[HCO]/d[H2] */
        J[14] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[13];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
        J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        J[124] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[233] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        J[234] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
        J[256] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[277] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        J[278] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f*sc[1];
        J[287] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        J[300] += dqdci;              /* dwdot[CH2O]/d[HCO] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[309] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[321] -= dqdci;              /* dwdot[HCO]/d[CH2O] */
        J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[409] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
        J[410] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[453] -= dqdci;              /* dwdot[HCO]/d[AR] */
        J[454] += dqdci;              /* dwdot[CH2O]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[13];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac + k_f*sc[1];
        dqdc[14] = dcdc_fac - k_r;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+13] -= dqdc[k];
            J[22*k+14] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[475] -= dqdT; /* dwdot[HCO]/dT */
    J[476] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18];
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = 1.0000000000000002e-06 * 540000000000
                * exp(0.45400000000000001 * tc[0] - 0.50321666580471969 * (2600) * invT);
    dlnkfdT = 0.45400000000000001 * invT + 0.50321666580471969 *  2600  * invT2;
    /* pressure-fall-off */
    k_0 = 2.2e+30 * exp(-4.7999999999999998 * tc[0] - 0.50321666580471969 * (5560) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -4.7999999999999998 * invT + 0.50321666580471969 * (5560) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.75800000000000001)*exp(-T/94);
    Fcent2 = 0.75800000000000001 * exp(-T/1555);
    Fcent3 = exp(-4200 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/94
        -Fcent2/1555
        + Fcent3*4200*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[15];
    Kc = refCinv * exp(g_RT[1] + g_RT[14] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[15]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[14] -= q; /* CH2O */
    wdot[15] += q; /* CH3O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[14] -= dqdci;               /* dwdot[CH2O]/d[H2] */
        J[15] += dqdci;               /* dwdot[CH3O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[14];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[36] -= dqdci;               /* dwdot[CH2O]/d[H] */
        J[37] += dqdci;               /* dwdot[CH3O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[124] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
        J[125] += dqdci;              /* dwdot[CH3O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[234] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
        J[235] += dqdci;              /* dwdot[CH3O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[256] -= dqdci;              /* dwdot[CH2O]/d[CO] */
        J[257] += dqdci;              /* dwdot[CH3O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[278] -= dqdci;              /* dwdot[CH2O]/d[CO2] */
        J[279] += dqdci;              /* dwdot[CH3O]/d[CO2] */
        /* d()/d[CH2O] */
        dqdci =  + k_f*sc[1];
        J[309] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
        J[323] += dqdci;              /* dwdot[CH3O]/d[CH2O] */
        /* d()/d[CH3O] */
        dqdci =  - k_r;
        J[331] -= dqdci;              /* dwdot[H]/d[CH3O] */
        J[344] -= dqdci;              /* dwdot[CH2O]/d[CH3O] */
        J[345] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[410] -= dqdci;              /* dwdot[CH2O]/d[C2H6] */
        J[411] += dqdci;              /* dwdot[CH3O]/d[C2H6] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[14];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac + k_f*sc[1];
        dqdc[15] = dcdc_fac - k_r;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+14] -= dqdc[k];
            J[22*k+15] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[476] -= dqdT; /* dwdot[CH2O]/dT */
    J[477] += dqdT; /* dwdot[CH3O]/dT */

    /*reaction 5: H + C2H4 (+M) <=> C2H5 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = 1.0000000000000002e-06 * 1080000000000
                * exp(0.45400000000000001 * tc[0] - 0.50321666580471969 * (1820) * invT);
    dlnkfdT = 0.45400000000000001 * invT + 0.50321666580471969 *  1820  * invT2;
    /* pressure-fall-off */
    k_0 = 1.1999999999999999e+42 * exp(-7.6200000000000001 * tc[0] - 0.50321666580471969 * (6970) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -7.6200000000000001 * invT + 0.50321666580471969 * (6970) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.97529999999999994)*exp(-T/210);
    Fcent2 = 0.97529999999999994 * exp(-T/984);
    Fcent3 = exp(-4374 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/210
        -Fcent2/984
        + Fcent3*4374*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[17];
    Kc = refCinv * exp(g_RT[1] + g_RT[16] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[16]) + (h_RT[17]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[16] -= q; /* C2H4 */
    wdot[17] += q; /* C2H5 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[16] -= dqdci;               /* dwdot[C2H4]/d[H2] */
        J[17] += dqdci;               /* dwdot[C2H5]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[16];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[38] -= dqdci;               /* dwdot[C2H4]/d[H] */
        J[39] += dqdci;               /* dwdot[C2H5]/d[H] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[126] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
        J[127] += dqdci;              /* dwdot[C2H5]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[236] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
        J[237] += dqdci;              /* dwdot[C2H5]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[258] -= dqdci;              /* dwdot[C2H4]/d[CO] */
        J[259] += dqdci;              /* dwdot[C2H5]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[280] -= dqdci;              /* dwdot[C2H4]/d[CO2] */
        J[281] += dqdci;              /* dwdot[C2H5]/d[CO2] */
        /* d()/d[C2H4] */
        dqdci =  + k_f*sc[1];
        J[353] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[368] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
        J[369] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
        /* d()/d[C2H5] */
        dqdci =  - k_r;
        J[375] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[390] -= dqdci;              /* dwdot[C2H4]/d[C2H5] */
        J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[412] -= dqdci;              /* dwdot[C2H4]/d[C2H6] */
        J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[456] -= dqdci;              /* dwdot[C2H4]/d[AR] */
        J[457] += dqdci;              /* dwdot[C2H5]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[16];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac + k_f*sc[1];
        dqdc[17] = dcdc_fac - k_r;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+16] -= dqdc[k];
            J[22*k+17] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[478] -= dqdT; /* dwdot[C2H4]/dT */
    J[479] += dqdT; /* dwdot[C2H5]/dT */

    /*reaction 6: H + C2H5 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[17];
    k_f = 1.0000000000000002e-06 * 5.21e+17
                * exp(-0.98999999999999999 * tc[0] - 0.50321666580471969 * (1580) * invT);
    dlnkfdT = -0.98999999999999999 * invT + 0.50321666580471969 *  1580  * invT2;
    /* pressure-fall-off */
    k_0 = 1.9900000000000001e+41 * exp(-7.0800000000000001 * tc[0] - 0.50321666580471969 * (6685) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -7.0800000000000001 * invT + 0.50321666580471969 * (6685) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.84219999999999995)*exp(-T/125);
    Fcent2 = 0.84219999999999995 * exp(-T/2219);
    Fcent3 = exp(-6882 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/125
        -Fcent2/2219
        + Fcent3*6882*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[18];
    Kc = refCinv * exp(g_RT[1] + g_RT[17] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[17]) + (h_RT[18]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[17] -= q; /* C2H5 */
    wdot[18] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[17] -= dqdci;               /* dwdot[C2H5]/d[H2] */
        J[18] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[17];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[39] -= dqdci;               /* dwdot[C2H5]/d[H] */
        J[40] += dqdci;               /* dwdot[C2H6]/d[H] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[127] -= dqdci;              /* dwdot[C2H5]/d[H2O] */
        J[128] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[237] -= dqdci;              /* dwdot[C2H5]/d[CH4] */
        J[238] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[259] -= dqdci;              /* dwdot[C2H5]/d[CO] */
        J[260] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[281] -= dqdci;              /* dwdot[C2H5]/d[CO2] */
        J[282] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[C2H5] */
        dqdci =  + k_f*sc[1];
        J[375] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[391] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
        J[392] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac - k_r;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[413] -= dqdci;              /* dwdot[C2H5]/d[C2H6] */
        J[414] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[457] -= dqdci;              /* dwdot[C2H5]/d[AR] */
        J[458] += dqdci;              /* dwdot[C2H6]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[17];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac + k_f*sc[1];
        dqdc[18] = 3*dcdc_fac - k_r;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+17] -= dqdc[k];
            J[22*k+18] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[479] -= dqdT; /* dwdot[C2H5]/dT */
    J[480] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 7: H2 + CO (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[0]*sc[11];
    k_f = 1.0000000000000002e-06 * 43000000
                * exp(1.5 * tc[0] - 0.50321666580471969 * (79600) * invT);
    dlnkfdT = 1.5 * invT + 0.50321666580471969 *  79600  * invT2;
    /* pressure-fall-off */
    k_0 = 5.0699999999999998e+27 * exp(-3.4199999999999999 * tc[0] - 0.50321666580471969 * (84350) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -3.4199999999999999 * invT + 0.50321666580471969 * (84350) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.93200000000000005)*exp(-T/197);
    Fcent2 = 0.93200000000000005 * exp(-T/1540);
    Fcent3 = exp(-10300 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/197
        -Fcent2/1540
        + Fcent3*10300*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[14];
    Kc = refCinv * exp(g_RT[0] + g_RT[11] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[11]) + (h_RT[14]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[11] -= q; /* CO */
    wdot[14] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac + k_f*sc[11];
        J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
        J[11] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[14] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[110] -= dqdci;              /* dwdot[H2]/d[H2O] */
        J[121] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[124] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[220] -= dqdci;              /* dwdot[H2]/d[CH4] */
        J[231] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[234] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac + k_f*sc[0];
        J[242] -= dqdci;              /* dwdot[H2]/d[CO] */
        J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[256] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[264] -= dqdci;              /* dwdot[H2]/d[CO2] */
        J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[278] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[308] -= dqdci;              /* dwdot[H2]/d[CH2O] */
        J[319] -= dqdci;              /* dwdot[CO]/d[CH2O] */
        J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac;
        J[396] -= dqdci;              /* dwdot[H2]/d[C2H6] */
        J[407] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[410] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[440] -= dqdci;              /* dwdot[H2]/d[AR] */
        J[451] -= dqdci;              /* dwdot[CO]/d[AR] */
        J[454] += dqdci;              /* dwdot[CH2O]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac + k_f*sc[11];
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac + k_f*sc[0];
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac - k_r;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+0] -= dqdc[k];
            J[22*k+11] -= dqdc[k];
            J[22*k+14] += dqdc[k];
        }
    }
    J[462] -= dqdT; /* dwdot[H2]/dT */
    J[473] -= dqdT; /* dwdot[CO]/dT */
    J[476] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 8: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = pow(sc[9], 2.000000);
    k_f = 1.0000000000000002e-06 * 21200000000000000
                * exp(-0.96999999999999997 * tc[0] - 0.50321666580471969 * (620) * invT);
    dlnkfdT = -0.96999999999999997 * invT + 0.50321666580471969 *  620  * invT2;
    /* pressure-fall-off */
    k_0 = 1.7700000000000001e+50 * exp(-9.6699999999999999 * tc[0] - 0.50321666580471969 * (6220) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -9.6699999999999999 * invT + 0.50321666580471969 * (6220) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.53249999999999997)*exp(-T/151);
    Fcent2 = 0.53249999999999997 * exp(-T/1038);
    Fcent3 = exp(-4970 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/151
        -Fcent2/1038
        + Fcent3*4970*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[18];
    Kc = refCinv * exp(2.000000*g_RT[9] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[9]) + (h_RT[18]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[9] -= 2 * q; /* CH3 */
    wdot[18] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[9] += -2 * dqdci;           /* dwdot[CH3]/d[H2] */
        J[18] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[119] += -2 * dqdci;         /* dwdot[CH3]/d[H2O] */
        J[128] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*2.000000*sc[9];
        J[207] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
        J[216] += dqdci;              /* dwdot[C2H6]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[229] += -2 * dqdci;         /* dwdot[CH3]/d[CH4] */
        J[238] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[251] += -2 * dqdci;         /* dwdot[CH3]/d[CO] */
        J[260] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[273] += -2 * dqdci;         /* dwdot[CH3]/d[CO2] */
        J[282] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac - k_r;
        J[405] += -2 * dqdci;         /* dwdot[CH3]/d[C2H6] */
        J[414] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[449] += -2 * dqdci;         /* dwdot[CH3]/d[AR] */
        J[458] += dqdci;              /* dwdot[C2H6]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac + k_f*2.000000*sc[9];
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac - k_r;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+9] += -2 * dqdc[k];
            J[22*k+18] += dqdc[k];
        }
    }
    J[471] += -2 * dqdT; /* dwdot[CH3]/dT */
    J[480] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 9: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-12 * 5e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[4]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[4] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[24] -= dqdci;               /* dwdot[O]/d[H] */
        J[26] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[1];
        J[45] -= dqdci;               /* dwdot[H]/d[O] */
        J[46] -= dqdci;               /* dwdot[O]/d[O] */
        J[48] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[89] -= dqdci;               /* dwdot[H]/d[OH] */
        J[90] -= dqdci;               /* dwdot[O]/d[OH] */
        J[92] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*q_nocor;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[112] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[114] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*q_nocor;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[222] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[224] += dqdci;              /* dwdot[OH]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*q_nocor;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[244] -= dqdci;              /* dwdot[O]/d[CO] */
        J[246] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*q_nocor;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[266] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[268] += dqdci;              /* dwdot[OH]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*q_nocor;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[398] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[400] += dqdci;              /* dwdot[OH]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[442] -= dqdci;              /* dwdot[O]/d[AR] */
        J[444] += dqdci;              /* dwdot[OH]/d[AR] */
    }
    else {
        dqdc[0] = 2*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[2];
        dqdc[2] = q_nocor + k_f*sc[1];
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor - k_r;
        dqdc[5] = 6*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = 2*q_nocor;
        dqdc[11] = 1.5*q_nocor;
        dqdc[12] = 2*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = 0.69999999999999996*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+2] -= dqdc[k];
            J[22*k+4] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[464] -= dqdT; /* dwdot[O]/dT */
    J[466] += dqdT; /* dwdot[OH]/dT */

    /*reaction 10: O + CO + M <=> CO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[3] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 3.5 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.5 - 1)*sc[20];
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = 1.0000000000000002e-12 * 602000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3000  * invT2;
    /* reverse */
    phi_r = sc[12];
    Kc = refCinv * exp(g_RT[2] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[11]) + (h_RT[12]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*q_nocor;
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[11] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[12] += dqdci;               /* dwdot[CO2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[11];
        J[46] -= dqdci;               /* dwdot[O]/d[O] */
        J[55] -= dqdci;               /* dwdot[CO]/d[O] */
        J[56] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[O2] */
        dqdci = (6 - 1)*q_nocor;
        J[68] -= dqdci;               /* dwdot[O]/d[O2] */
        J[77] -= dqdci;               /* dwdot[CO]/d[O2] */
        J[78] += dqdci;               /* dwdot[CO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*q_nocor;
        J[112] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[121] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[122] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*q_nocor;
        J[222] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[231] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[232] += dqdci;              /* dwdot[CO2]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*q_nocor + k_f*sc[2];
        J[244] -= dqdci;              /* dwdot[O]/d[CO] */
        J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.5 - 1)*q_nocor - k_r;
        J[266] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*q_nocor;
        J[398] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[407] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[408] += dqdci;              /* dwdot[CO2]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.5 - 1)*q_nocor;
        J[442] -= dqdci;              /* dwdot[O]/d[AR] */
        J[451] -= dqdci;              /* dwdot[CO]/d[AR] */
        J[452] += dqdci;              /* dwdot[CO2]/d[AR] */
    }
    else {
        dqdc[0] = 2*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*sc[11];
        dqdc[3] = 6*q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = 6*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = 2*q_nocor;
        dqdc[11] = 1.5*q_nocor + k_f*sc[2];
        dqdc[12] = 3.5*q_nocor - k_r;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = 0.5*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+2] -= dqdc[k];
            J[22*k+11] -= dqdc[k];
            J[22*k+12] += dqdc[k];
        }
    }
    J[464] -= dqdT; /* dwdot[O]/dT */
    J[473] -= dqdT; /* dwdot[CO]/dT */
    J[474] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 11: H + O2 + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[3] + ( 0 - 1)*sc[5] + ( 0.75 - 1)*sc[11] + ( 1.5 - 1)*sc[12] + ( 1.5 - 1)*sc[18] + ( 0 - 1)*sc[19] + ( 0 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-12 * 2.8e+18
                * exp(-0.85999999999999999 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.85999999999999999 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[3];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[25] -= dqdci;               /* dwdot[O2]/d[H] */
        J[28] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[O2] */
        dqdci = (0 - 1)*q_nocor + k_f*sc[1];
        J[67] -= dqdci;               /* dwdot[H]/d[O2] */
        J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[113] -= dqdci;              /* dwdot[O2]/d[H2O] */
        J[116] += dqdci;              /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
        J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
        /* d()/d[CO] */
        dqdci = (0.75 - 1)*q_nocor;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[248] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (1.5 - 1)*q_nocor;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[267] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[270] += dqdci;              /* dwdot[HO2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (1.5 - 1)*q_nocor;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[399] -= dqdci;              /* dwdot[O2]/d[C2H6] */
        J[402] += dqdci;              /* dwdot[HO2]/d[C2H6] */
        /* d()/d[N2] */
        dqdci = (0 - 1)*q_nocor;
        J[419] -= dqdci;              /* dwdot[H]/d[N2] */
        J[421] -= dqdci;              /* dwdot[O2]/d[N2] */
        J[424] += dqdci;              /* dwdot[HO2]/d[N2] */
        /* d()/d[AR] */
        dqdci = (0 - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[443] -= dqdci;              /* dwdot[O2]/d[AR] */
        J[446] += dqdci;              /* dwdot[HO2]/d[AR] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*sc[3];
        dqdc[2] = q_nocor;
        dqdc[3] =  + k_f*sc[1];
        dqdc[4] = q_nocor;
        dqdc[6] = q_nocor - k_r;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = 0.75*q_nocor;
        dqdc[12] = 1.5*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 1.5*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+3] -= dqdc[k];
            J[22*k+6] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[465] -= dqdT; /* dwdot[O2]/dT */
    J[468] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 12: 2.000000 H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 0 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.63 - 1)*sc[20];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 1e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[0];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1]) + (h_RT[0]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (0 - 1)*q_nocor - k_r;
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*2.000000*sc[1];
        J[22] += dqdci;               /* dwdot[H2]/d[H] */
        J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor;
        J[110] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[111] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*q_nocor;
        J[220] += dqdci;              /* dwdot[H2]/d[CH4] */
        J[221] += -2 * dqdci;         /* dwdot[H]/d[CH4] */
        /* d()/d[CO2] */
        dqdci = (0 - 1)*q_nocor;
        J[264] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[265] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*q_nocor;
        J[396] += dqdci;              /* dwdot[H2]/d[C2H6] */
        J[397] += -2 * dqdci;         /* dwdot[H]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.63 - 1)*q_nocor;
        J[440] += dqdci;              /* dwdot[H2]/d[AR] */
        J[441] += -2 * dqdci;         /* dwdot[H]/d[AR] */
    }
    else {
        dqdc[0] =  - k_r;
        dqdc[1] = q_nocor + k_f*2.000000*sc[1];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = 2*q_nocor;
        dqdc[11] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = 0.63*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+0] += dqdc[k];
            J[22*k+1] += -2 * dqdc[k];
        }
    }
    J[462] += dqdT; /* dwdot[H2]/dT */
    J[463] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 13: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0.72999999999999998 - 1)*sc[0] + ( 3.6499999999999999 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 3 - 1)*sc[18] + ( 0.38 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = 1.0000000000000002e-12 * 2.2e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[4]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (0.72999999999999998 - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[26] -= dqdci;               /* dwdot[OH]/d[H] */
        J[27] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[1];
        J[89] -= dqdci;               /* dwdot[H]/d[OH] */
        J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (3.6499999999999999 - 1)*q_nocor - k_r;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*q_nocor;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[224] -= dqdci;              /* dwdot[OH]/d[CH4] */
        J[225] += dqdci;              /* dwdot[H2O]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*q_nocor;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[400] -= dqdci;              /* dwdot[OH]/d[C2H6] */
        J[401] += dqdci;              /* dwdot[H2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.38 - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[444] -= dqdci;              /* dwdot[OH]/d[AR] */
        J[445] += dqdci;              /* dwdot[H2O]/d[AR] */
    }
    else {
        dqdc[0] = 0.72999999999999998*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[4];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f*sc[1];
        dqdc[5] = 3.6499999999999999*q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = 2*q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = 0.38*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+4] -= dqdc[k];
            J[22*k+5] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[466] -= dqdT; /* dwdot[OH]/dT */
    J[467] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 14: HCO + M <=> H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 0 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18];
    /* forward */
    phi_f = sc[13];
    k_f = 1.0000000000000002e-06 * 1.87e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  17000  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[11];
    Kc = refC * exp(-g_RT[1] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13]) + (h_RT[1] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*q_nocor;
        J[1] += dqdci;                /* dwdot[H]/d[H2] */
        J[11] += dqdci;               /* dwdot[CO]/d[H2] */
        J[13] -= dqdci;               /* dwdot[HCO]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[11];
        J[23] += dqdci;               /* dwdot[H]/d[H] */
        J[33] += dqdci;               /* dwdot[CO]/d[H] */
        J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor;
        J[111] += dqdci;              /* dwdot[H]/d[H2O] */
        J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*q_nocor;
        J[221] += dqdci;              /* dwdot[H]/d[CH4] */
        J[231] += dqdci;              /* dwdot[CO]/d[CH4] */
        J[233] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*q_nocor - k_r*sc[1];
        J[243] += dqdci;              /* dwdot[H]/d[CO] */
        J[253] += dqdci;              /* dwdot[CO]/d[CO] */
        J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*q_nocor;
        J[265] += dqdci;              /* dwdot[H]/d[CO2] */
        J[275] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[277] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[287] += dqdci;              /* dwdot[H]/d[HCO] */
        J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*q_nocor;
        J[397] += dqdci;              /* dwdot[H]/d[C2H6] */
        J[407] += dqdci;              /* dwdot[CO]/d[C2H6] */
        J[409] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
    }
    else {
        dqdc[0] = 2*q_nocor;
        dqdc[1] = q_nocor - k_r*sc[11];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = 2*q_nocor;
        dqdc[11] = 1.5*q_nocor - k_r*sc[1];
        dqdc[12] = 2*q_nocor;
        dqdc[13] = q_nocor + k_f;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] += dqdc[k];
            J[22*k+11] += dqdc[k];
            J[22*k+13] -= dqdc[k];
        }
    }
    J[463] += dqdT; /* dwdot[H]/dT */
    J[473] += dqdT; /* dwdot[CO]/dT */
    J[475] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 15: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = 1.0000000000000002e-06 * 50000
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6290) * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  6290  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[1] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[2] -= dqdci;                /* dwdot[O]/d[H2] */
    J[4] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[44] -= dqdci;               /* dwdot[H2]/d[O] */
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[88] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 16: O + HO2 <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[47] += dqdci;               /* dwdot[O2]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[50] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[68] -= dqdci;               /* dwdot[O]/d[O2] */
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[91] += dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[134] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 17: O + CH2 <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 80000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[7] -= q; /* CH2 */
    wdot[13] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[35] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[51] -= dqdci;               /* dwdot[CH2]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[2];
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[156] -= dqdci;              /* dwdot[O]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[167] += dqdci;              /* dwdot[HCO]/d[CH2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[293] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 18: O + CH2(S) <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = 1.0000000000000002e-06 * 15000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[8] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[8] -= q; /* CH2(S) */
    wdot[13] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[35] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[8];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[52] -= dqdci;               /* dwdot[CH2(S)]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[2];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[178] -= dqdci;              /* dwdot[O]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[189] += dqdci;              /* dwdot[HCO]/d[CH2(S)] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[294] -= dqdci;              /* dwdot[CH2(S)]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 19: O + CH3 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[9];
    k_f = 1.0000000000000002e-06 * 84300000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[9] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[9]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[9] -= q; /* CH3 */
    wdot[14] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[9];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[58] += dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[2];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[212] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[309] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[310] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[317] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 20: O + CH4 <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[10];
    k_f = 1.0000000000000002e-06 * 1020000000
                * exp(1.5 * tc[0] - 0.50321666580471969 * (8600) * invT);
    dlnkfdT = 1.5 * invT + 0.50321666580471969 *  8600  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[9];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[10]) + (h_RT[4] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[9] += q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[53] += dqdci;               /* dwdot[CH3]/d[O] */
    J[54] -= dqdci;               /* dwdot[CH4]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[98] -= dqdci;               /* dwdot[CH4]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[222] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[224] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 21: O + HCO <=> OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[13];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[11];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[4] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[55] += dqdci;               /* dwdot[CO]/d[O] */
    J[57] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[11];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[99] += dqdci;               /* dwdot[CO]/d[OH] */
    J[101] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[4];
    J[244] -= dqdci;              /* dwdot[O]/d[CO] */
    J[246] += dqdci;              /* dwdot[OH]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[290] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 22: O + HCO <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[13];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[12] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[12] += q; /* CO2 */
    wdot[13] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[34] += dqdci;               /* dwdot[CO2]/d[H] */
    J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[56] += dqdci;               /* dwdot[CO2]/d[O] */
    J[57] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[265] += dqdci;              /* dwdot[H]/d[CO2] */
    J[266] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[277] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[298] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 23: O + CH2O <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[14];
    k_f = 1.0000000000000002e-06 * 39000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3540) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3540  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[13];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[14]) + (h_RT[4] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[14];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    J[58] -= dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[101] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[102] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[290] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[2];
    J[310] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[312] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 24: O + C2H4 <=> CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[16];
    k_f = 1.0000000000000002e-06 * 19200000
                * exp(1.8300000000000001 * tc[0] - 0.50321666580471969 * (220) * invT);
    dlnkfdT = 1.8300000000000001 * invT + 0.50321666580471969 *  220  * invT2;
    /* reverse */
    phi_r = sc[9]*sc[13];
    Kc = exp(g_RT[2] - g_RT[9] - g_RT[13] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[16]) + (h_RT[9] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[9] += q; /* CH3 */
    wdot[13] += q; /* HCO */
    wdot[16] -= q; /* C2H4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] += dqdci;               /* dwdot[CH3]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    J[60] -= dqdci;               /* dwdot[C2H4]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[13];
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[211] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[214] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[9];
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[295] += dqdci;              /* dwdot[CH3]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[302] -= dqdci;              /* dwdot[C2H4]/d[HCO] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[2];
    J[354] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[361] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[365] += dqdci;              /* dwdot[HCO]/d[C2H4] */
    J[368] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[478] -= dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 25: O + C2H5 <=> CH3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = 1.0000000000000002e-06 * 132000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[9]*sc[14];
    Kc = exp(g_RT[2] - g_RT[9] - g_RT[14] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[9] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[9] += q; /* CH3 */
    wdot[14] += q; /* CH2O */
    wdot[17] -= q; /* C2H5 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] += dqdci;               /* dwdot[CH3]/d[O] */
    J[58] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[61] -= dqdci;               /* dwdot[C2H5]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[14];
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[212] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[215] -= dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[9];
    J[310] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[317] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[325] -= dqdci;              /* dwdot[C2H5]/d[CH2O] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[2];
    J[376] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[383] += dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[388] += dqdci;              /* dwdot[CH2O]/d[C2H5] */
    J[391] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */
    J[479] -= dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 26: O + C2H6 <=> OH + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = 1.0000000000000002e-06 * 89800000
                * exp(1.9199999999999999 * tc[0] - 0.50321666580471969 * (5690) * invT);
    dlnkfdT = 1.9199999999999999 * invT + 0.50321666580471969 *  5690  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[17];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[18]) + (h_RT[4] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[18];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[61] += dqdci;               /* dwdot[C2H5]/d[O] */
    J[62] -= dqdci;               /* dwdot[C2H6]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[17];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[106] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[4];
    J[376] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[378] += dqdci;              /* dwdot[OH]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[2];
    J[398] -= dqdci;              /* dwdot[O]/d[C2H6] */
    J[400] += dqdci;              /* dwdot[OH]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 27: O2 + CO <=> O + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = 1.0000000000000002e-06 * 2500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (47800) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  47800  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[12];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[11]) + (h_RT[2] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[12];
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[47] -= dqdci;               /* dwdot[O2]/d[O] */
    J[55] -= dqdci;               /* dwdot[CO]/d[O] */
    J[56] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[68] += dqdci;               /* dwdot[O]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[77] -= dqdci;               /* dwdot[CO]/d[O2] */
    J[78] += dqdci;               /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3];
    J[244] += dqdci;              /* dwdot[O]/d[CO] */
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[2];
    J[266] += dqdci;              /* dwdot[O]/d[CO2] */
    J[267] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[473] -= dqdT;               /* dwdot[CO]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 28: O2 + CH2O <=> HO2 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[14];
    k_f = 1.0000000000000002e-06 * 100000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (40000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  40000  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[13];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[14]) + (h_RT[6] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[79] += dqdci;               /* dwdot[HCO]/d[O2] */
    J[80] -= dqdci;               /* dwdot[CH2O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[13];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[145] += dqdci;              /* dwdot[HCO]/d[HO2] */
    J[146] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[6];
    J[289] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[292] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[3];
    J[311] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[314] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 29: H + 2.000000 O2 <=> HO2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*pow(sc[3], 2.000000);
    k_f = 1.0000000000000002e-12 * 3e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[6];
    Kc = refCinv * exp(g_RT[1] + 2.000000*g_RT[3] - g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + 2.000000*h_RT[3]) + (h_RT[3] + h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*pow(sc[3], 2.000000);
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*2.000000*sc[3] - k_r*sc[6];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 30: H + O2 + H2O <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = 1.0000000000000002e-12 * 9.38e+18
                * exp(-0.76000000000000001 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.76000000000000001 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[5]) + (h_RT[5] + h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[5];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[5];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[113] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[116] += dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 31: H + O2 + N2 <=> HO2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[19];
    k_f = 1.0000000000000002e-12 * 3.75e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[19];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[19] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[19]) + (h_RT[6] + h_RT[19]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[19];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[19];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[19];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[419] -= dqdci;              /* dwdot[H]/d[N2] */
    J[421] -= dqdci;              /* dwdot[O2]/d[N2] */
    J[424] += dqdci;              /* dwdot[HO2]/d[N2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 32: H + O2 + AR <=> HO2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[20];
    k_f = 1.0000000000000002e-12 * 7e+17
                * exp(-0.80000000000000004 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.80000000000000004 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[20];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[20] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[20]) + (h_RT[6] + h_RT[20]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[20];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[20];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[20];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[441] -= dqdci;              /* dwdot[H]/d[AR] */
    J[443] -= dqdci;              /* dwdot[O2]/d[AR] */
    J[446] += dqdci;              /* dwdot[HO2]/d[AR] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 33: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 83000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (14413) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  14413  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[24] += dqdci;               /* dwdot[O]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[45] -= dqdci;               /* dwdot[H]/d[O] */
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[47] -= dqdci;               /* dwdot[O2]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[68] += dqdci;               /* dwdot[O]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[89] -= dqdci;               /* dwdot[H]/d[OH] */
    J[90] += dqdci;               /* dwdot[O]/d[OH] */
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 34: 2.000000 H + H2 <=> 2.000000 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 90000000000000000
                * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.59999999999999998 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000);
    Kc = refCinv * exp(g_RT[0] - 2.000000*g_RT[0] + 2.000000*g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + 2.000000*h_RT[1]) + (2.000000*h_RT[0]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*2.000000*sc[0];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[0]*2.000000*sc[1];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 35: 2.000000 H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[5];
    k_f = 1.0000000000000002e-12 * 6e+19
                * exp(-1.25 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.25 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[5];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[5]) + (h_RT[0] + h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[5];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[110] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[111] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 36: 2.000000 H + CO2 <=> H2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[12];
    k_f = 1.0000000000000002e-12 * 5.5e+20
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[12];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[12] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[12]) + (h_RT[0] + h_RT[12]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[12];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[12];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[CO2] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[264] += dqdci;              /* dwdot[H2]/d[CO2] */
    J[265] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 37: H + HO2 <=> O2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 28000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (1068) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  1068  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[3] += dqdci;                /* dwdot[O2]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] += dqdci;               /* dwdot[O2]/d[H] */
    J[28] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[66] += dqdci;               /* dwdot[H2]/d[O2] */
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[132] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 38: H + HO2 <=> 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 134000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (635) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  635  * invT2;
    /* reverse */
    phi_r = pow(sc[4], 2.000000);
    Kc = exp(g_RT[1] - 2.000000*g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (2.000000*h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[26] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[28] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[4];
    J[89] -= dqdci;               /* dwdot[H]/d[OH] */
    J[92] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[136] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[466] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 39: H + CH4 <=> CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = 1.0000000000000002e-06 * 660000000
                * exp(1.6200000000000001 * tc[0] - 0.50321666580471969 * (10840) * invT);
    dlnkfdT = 1.6200000000000001 * invT + 0.50321666580471969 *  10840  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[9];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[0] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[9] += q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[9];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
    J[10] -= dqdci;               /* dwdot[CH4]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    J[32] -= dqdci;               /* dwdot[CH4]/d[H] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[0];
    J[198] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[1];
    J[220] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 40: H + HCO <=> H2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = 1.0000000000000002e-06 * 73400000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[11];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[0] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[11];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[11] += dqdci;               /* dwdot[CO]/d[H2] */
    J[13] -= dqdci;               /* dwdot[HCO]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[13];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[33] += dqdci;               /* dwdot[CO]/d[H] */
    J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[0];
    J[242] += dqdci;              /* dwdot[H2]/d[CO] */
    J[243] -= dqdci;              /* dwdot[H]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[286] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[287] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 41: H + CH2O <=> HCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = 1.0000000000000002e-06 * 23000000000
                * exp(1.05 * tc[0] - 0.50321666580471969 * (3275) * invT);
    dlnkfdT = 1.05 * invT + 0.50321666580471969 *  3275  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[13];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[0] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[13];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[13] += dqdci;               /* dwdot[HCO]/d[H2] */
    J[14] -= dqdci;               /* dwdot[CH2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[14];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += dqdci;               /* dwdot[HCO]/d[H] */
    J[36] -= dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[0];
    J[286] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[287] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[1];
    J[308] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[309] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 42: H + CH3O <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[15];
    k_f = 1.0000000000000002e-06 * 32000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[9];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[9] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[15]) + (h_RT[4] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += q; /* OH */
    wdot[9] += q; /* CH3 */
    wdot[15] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[15];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    J[37] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[89] -= dqdci;               /* dwdot[H]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[103] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[213] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[1];
    J[331] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[334] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[339] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[345] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[477] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 43: H + C2H6 <=> C2H5 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[18];
    k_f = 1.0000000000000002e-06 * 115000000
                * exp(1.8999999999999999 * tc[0] - 0.50321666580471969 * (7530) * invT);
    dlnkfdT = 1.8999999999999999 * invT + 0.50321666580471969 *  7530  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[17];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[18]) + (h_RT[0] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[17];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[17] += dqdci;               /* dwdot[C2H5]/d[H2] */
    J[18] -= dqdci;               /* dwdot[C2H6]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[39] += dqdci;               /* dwdot[C2H5]/d[H] */
    J[40] -= dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[0];
    J[374] += dqdci;              /* dwdot[H2]/d[C2H5] */
    J[375] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[1];
    J[396] += dqdci;              /* dwdot[H2]/d[C2H6] */
    J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 44: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  3430  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
    J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[27] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[88] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[110] -= dqdci;              /* dwdot[H2]/d[H2O] */
    J[111] += dqdci;              /* dwdot[H]/d[H2O] */
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 45: 2.000000 OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-06 * 35700
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (-2110) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  -2110  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + 2.000000*g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[4] -= 2 * q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[48] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[49] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[4];
    J[90] += dqdci;               /* dwdot[O]/d[OH] */
    J[92] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[112] += dqdci;              /* dwdot[O]/d[H2O] */
    J[114] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[466] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 46: OH + HO2 <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = 1.0000000000000002e-06 * 29000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -500  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[5];
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[70] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[91] += dqdci;               /* dwdot[O2]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[113] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[116] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[136] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[137] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 47: OH + CH2 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[7] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[7] -= q; /* CH2 */
    wdot[14] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[95] -= dqdci;               /* dwdot[CH2]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[158] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[168] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[309] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[312] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[315] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 48: OH + CH2(S) <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[8];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[8] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[8]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[8] -= q; /* CH2(S) */
    wdot[14] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[96] -= dqdci;               /* dwdot[CH2(S)]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[4];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[180] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[190] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[309] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[312] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[316] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 49: OH + CH3 <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = 1.0000000000000002e-06 * 56000000
                * exp(1.6000000000000001 * tc[0] - 0.50321666580471969 * (5420) * invT);
    dlnkfdT = 1.6000000000000001 * invT + 0.50321666580471969 *  5420  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[7] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[7] += q; /* CH2 */
    wdot[9] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[95] += dqdci;               /* dwdot[CH2]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[7];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[117] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[119] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[5];
    J[158] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[159] += dqdci;              /* dwdot[H2O]/d[CH2] */
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[202] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[203] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[205] += dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 50: OH + CH3 <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = 1.0000000000000002e-06 * 25010000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[8];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[8] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[5] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[8] += q; /* CH2(S) */
    wdot[9] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[96] += dqdci;               /* dwdot[CH2(S)]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[118] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[119] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[5];
    J[180] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[181] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[184] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[202] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[203] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[206] += dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[470] += dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 51: OH + CH4 <=> CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = 1.0000000000000002e-06 * 100000000
                * exp(1.6000000000000001 * tc[0] - 0.50321666580471969 * (3120) * invT);
    dlnkfdT = 1.6000000000000001 * invT + 0.50321666580471969 *  3120  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[9] += q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[97] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[98] -= dqdci;               /* dwdot[CH4]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[9];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[119] += dqdci;              /* dwdot[CH3]/d[H2O] */
    J[120] -= dqdci;              /* dwdot[CH4]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[5];
    J[202] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[203] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[4];
    J[224] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[225] += dqdci;              /* dwdot[H2O]/d[CH4] */
    J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 52: OH + CO <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = 1.0000000000000002e-06 * 47600000
                * exp(1.228 * tc[0] - 0.50321666580471969 * (70) * invT);
    dlnkfdT = 1.228 * invT + 0.50321666580471969 *  70  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[33] -= dqdci;               /* dwdot[CO]/d[H] */
    J[34] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[99] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[100] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[4];
    J[243] += dqdci;              /* dwdot[H]/d[CO] */
    J[246] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[265] += dqdci;              /* dwdot[H]/d[CO2] */
    J[268] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[473] -= dqdT;               /* dwdot[CO]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 53: OH + HCO <=> H2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[13];
    k_f = 1.0000000000000002e-06 * 50000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[11];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[5] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[99] += dqdci;               /* dwdot[CO]/d[OH] */
    J[101] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[246] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[247] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[4];
    J[290] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[291] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 54: OH + CH2O <=> HCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[14];
    k_f = 1.0000000000000002e-06 * 3430000000
                * exp(1.1799999999999999 * tc[0] - 0.50321666580471969 * (-447) * invT);
    dlnkfdT = 1.1799999999999999 * invT + 0.50321666580471969 *  -447  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[13];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[14]) + (h_RT[5] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[14];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[101] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[102] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[13];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[123] += dqdci;              /* dwdot[HCO]/d[H2O] */
    J[124] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[5];
    J[290] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[291] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[4];
    J[312] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[313] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 55: OH + C2H6 <=> C2H5 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[18];
    k_f = 1.0000000000000002e-06 * 3540000
                * exp(2.1200000000000001 * tc[0] - 0.50321666580471969 * (870) * invT);
    dlnkfdT = 2.1200000000000001 * invT + 0.50321666580471969 *  870  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[17];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[18]) + (h_RT[5] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[18];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[105] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[106] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[17];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[127] += dqdci;              /* dwdot[C2H5]/d[H2O] */
    J[128] -= dqdci;              /* dwdot[C2H6]/d[H2O] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[5];
    J[378] -= dqdci;              /* dwdot[OH]/d[C2H5] */
    J[379] += dqdci;              /* dwdot[H2O]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[4];
    J[400] -= dqdci;              /* dwdot[OH]/d[C2H6] */
    J[401] += dqdci;              /* dwdot[H2O]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 56: HO2 + CH2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[7];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[7] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[7]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[7] -= q; /* CH2 */
    wdot[14] += q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[95] -= dqdci;               /* dwdot[CH2]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[7];
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[139] -= dqdci;              /* dwdot[CH2]/d[HO2] */
    J[146] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[6];
    J[158] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[160] -= dqdci;              /* dwdot[HO2]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[168] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[312] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[314] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[315] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 57: HO2 + CH3 <=> O2 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[9];
    k_f = 1.0000000000000002e-06 * 1000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[10];
    Kc = exp(-g_RT[3] + g_RT[6] + g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[3] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[10];
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    J[75] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[76] += dqdci;               /* dwdot[CH4]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[9];
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[141] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[142] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[201] += dqdci;              /* dwdot[O2]/d[CH3] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[3];
    J[223] += dqdci;              /* dwdot[O2]/d[CH4] */
    J[226] -= dqdci;              /* dwdot[HO2]/d[CH4] */
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */

    /*reaction 58: HO2 + CH3 <=> OH + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[9];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[15];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[9] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[4] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[9] -= q; /* CH3 */
    wdot[15] += q; /* CH3O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    J[103] += dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[9];
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[141] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[147] += dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[213] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[4];
    J[334] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[336] -= dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[339] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[345] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[477] += dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 59: HO2 + CO <=> OH + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[11];
    k_f = 1.0000000000000002e-06 * 150000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (23600) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  23600  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[12];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[11]) + (h_RT[4] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[12];
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[99] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[100] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[143] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[144] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[246] += dqdci;              /* dwdot[OH]/d[CO] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[4];
    J[268] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[270] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[473] -= dqdT;               /* dwdot[CO]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 60: CH2 + O2 <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 13200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (1500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  1500  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[13];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[7] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[7] -= q; /* CH2 */
    wdot[13] += q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[7];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[73] -= dqdci;               /* dwdot[CH2]/d[O2] */
    J[79] += dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[95] -= dqdci;               /* dwdot[CH2]/d[OH] */
    J[101] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[3];
    J[157] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[158] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[167] += dqdci;              /* dwdot[HCO]/d[CH2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[289] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[290] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[293] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 61: CH2 + H2 <=> H + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = 1.0000000000000002e-06 * 500000
                * exp(2 * tc[0] - 0.50321666580471969 * (7230) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  7230  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[7] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* CH2 */
    wdot[9] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[7] -= dqdci;                /* dwdot[CH2]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[0];
    J[154] -= dqdci;              /* dwdot[H2]/d[CH2] */
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] += dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[198] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 62: CH2 + CH3 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[9];
    k_f = 1.0000000000000002e-06 * 40000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[7] + g_RT[9] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[9]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* CH2 */
    wdot[9] -= q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[38] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[9];
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    J[170] += dqdci;              /* dwdot[C2H4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[214] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[353] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[359] -= dqdci;              /* dwdot[CH2]/d[C2H4] */
    J[361] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[368] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[478] += dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 63: CH2 + CH4 <=> 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[10];
    k_f = 1.0000000000000002e-06 * 2460000
                * exp(2 * tc[0] - 0.50321666580471969 * (8270) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  8270  * invT2;
    /* reverse */
    phi_r = pow(sc[9], 2.000000);
    Kc = exp(g_RT[7] - 2.000000*g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[10]) + (2.000000*h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* CH2 */
    wdot[9] += 2 * q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[10];
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] += 2 * dqdci;          /* dwdot[CH3]/d[CH2] */
    J[164] -= dqdci;              /* dwdot[CH4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[9];
    J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[7];
    J[227] -= dqdci;              /* dwdot[CH2]/d[CH4] */
    J[229] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[471] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 64: CH2(S) + N2 <=> CH2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[19];
    k_f = 1.0000000000000002e-06 * 15000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (600) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  600  * invT2;
    /* reverse */
    phi_r = sc[7]*sc[19];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[19] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[19]) + (h_RT[7] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[19];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[19];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[425] += dqdci;              /* dwdot[CH2]/d[N2] */
    J[426] -= dqdci;              /* dwdot[CH2(S)]/d[N2] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 65: CH2(S) + AR <=> CH2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[20];
    k_f = 1.0000000000000002e-06 * 9000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (600) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  600  * invT2;
    /* reverse */
    phi_r = sc[7]*sc[20];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[20] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[20]) + (h_RT[7] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[20];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[20];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[447] += dqdci;              /* dwdot[CH2]/d[AR] */
    J[448] -= dqdci;              /* dwdot[CH2(S)]/d[AR] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 66: CH2(S) + O2 <=> H + OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = 1.0000000000000002e-06 * 28000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4]*sc[11];
    Kc = refC * exp(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[8] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[1] + h_RT[4] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[8] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4]*sc[11];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[33] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[67] += dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[74] -= dqdci;               /* dwdot[CH2(S)]/d[O2] */
    J[77] += dqdci;               /* dwdot[CO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1]*sc[11];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[96] -= dqdci;               /* dwdot[CH2(S)]/d[OH] */
    J[99] += dqdci;               /* dwdot[CO]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[179] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    J[180] += dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[187] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[4];
    J[243] += dqdci;              /* dwdot[H]/d[CO] */
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[246] += dqdci;              /* dwdot[OH]/d[CO] */
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 67: CH2(S) + O2 <=> CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = 1.0000000000000002e-06 * 12000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[11];
    Kc = exp(g_RT[3] - g_RT[5] + g_RT[8] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[5] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[5] += q; /* H2O */
    wdot[8] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[74] -= dqdci;               /* dwdot[CH2(S)]/d[O2] */
    J[77] += dqdci;               /* dwdot[CO]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[113] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[118] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[179] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    J[181] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[187] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[247] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 68: CH2(S) + H2 <=> CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[8];
    k_f = 1.0000000000000002e-06 * 70000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[8]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[8] -= q; /* CH2(S) */
    wdot[9] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[8];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[8] -= dqdci;                /* dwdot[CH2(S)]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[0];
    J[176] -= dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] += dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[198] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[206] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 69: CH2(S) + H2O <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[8];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(g_RT[5] - g_RT[5] - g_RT[7] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[8]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[117] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[118] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[5];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[5];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 70: CH2(S) + CH3 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[9];
    k_f = 1.0000000000000002e-06 * 12000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-570) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -570  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[8] + g_RT[9] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[9]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[8] -= q; /* CH2(S) */
    wdot[9] -= q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[38] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[9];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    J[192] += dqdci;              /* dwdot[C2H4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[206] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[214] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[353] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[360] -= dqdci;              /* dwdot[CH2(S)]/d[C2H4] */
    J[361] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[368] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[478] += dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 71: CH2(S) + CH4 <=> 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = 1.0000000000000002e-06 * 16000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-570) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -570  * invT2;
    /* reverse */
    phi_r = pow(sc[9], 2.000000);
    Kc = exp(g_RT[8] - 2.000000*g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[10]) + (2.000000*h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= q; /* CH2(S) */
    wdot[9] += 2 * q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[10];
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] += 2 * dqdci;          /* dwdot[CH3]/d[CH2(S)] */
    J[186] -= dqdci;              /* dwdot[CH4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[9];
    J[206] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[8];
    J[228] -= dqdci;              /* dwdot[CH2(S)]/d[CH4] */
    J[229] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 72: CH2(S) + CO <=> CH2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = 1.0000000000000002e-06 * 9000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[7]*sc[11];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[11] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[11]) + (h_RT[7] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[11];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[11];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[249] += dqdci;              /* dwdot[CH2]/d[CO] */
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 73: CH2(S) + CO2 <=> CH2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = 1.0000000000000002e-06 * 7000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[7]*sc[12];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[12] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[7] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[12];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[12];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[271] += dqdci;              /* dwdot[CH2]/d[CO2] */
    J[272] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 74: CH2(S) + CO2 <=> CO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = 1.0000000000000002e-06 * 14000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[11]*sc[14];
    Kc = exp(g_RT[8] - g_RT[11] + g_RT[12] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[11] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    wdot[12] -= q; /* CO2 */
    wdot[14] += q; /* CH2O */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[12];
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[187] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[188] -= dqdci;              /* dwdot[CO2]/d[CH2(S)] */
    J[190] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[14];
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[254] -= dqdci;              /* dwdot[CO2]/d[CO] */
    J[256] += dqdci;              /* dwdot[CH2O]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[8];
    J[272] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    J[275] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    J[278] += dqdci;              /* dwdot[CH2O]/d[CO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[11];
    J[316] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[319] += dqdci;              /* dwdot[CO]/d[CH2O] */
    J[320] -= dqdci;              /* dwdot[CO2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[474] -= dqdT;               /* dwdot[CO2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 75: CH3 + O2 <=> O + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = 1.0000000000000002e-06 * 26750000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (28800) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  28800  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[15];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[9] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (h_RT[2] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[9] -= q; /* CH3 */
    wdot[15] += q; /* CH3O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[15];
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[47] -= dqdci;               /* dwdot[O2]/d[O] */
    J[53] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[59] += dqdci;               /* dwdot[CH3O]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[68] += dqdci;               /* dwdot[O]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[75] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[81] += dqdci;               /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[200] += dqdci;              /* dwdot[O]/d[CH3] */
    J[201] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[213] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[2];
    J[332] += dqdci;              /* dwdot[O]/d[CH3O] */
    J[333] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[339] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[345] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[477] += dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 76: CH3 + O2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = 1.0000000000000002e-06 * 36000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (8940) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  8940  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[9] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[9] -= q; /* CH3 */
    wdot[14] += q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[80] += dqdci;               /* dwdot[CH2O]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[201] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[212] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[311] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[312] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[317] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 77: 2.000000 CH3 <=> H + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[9], 2.000000);
    k_f = 1.0000000000000002e-06 * 4990000000000
                * exp(0.10000000000000001 * tc[0] - 0.50321666580471969 * (10600) * invT);
    dlnkfdT = 0.10000000000000001 * invT + 0.50321666580471969 *  10600  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[17];
    Kc = exp(-g_RT[1] + 2.000000*g_RT[9] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[9]) + (h_RT[1] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[9] -= 2 * q; /* CH3 */
    wdot[17] += q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[17];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[31] += -2 * dqdci;          /* dwdot[CH3]/d[H] */
    J[39] += dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*2.000000*sc[9];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[207] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[215] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[1];
    J[375] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[383] += -2 * dqdci;         /* dwdot[CH3]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[471] += -2 * dqdT;          /* dwdot[CH3]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 78: CH3 + HCO <=> CH4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[13];
    k_f = 1.0000000000000002e-06 * 26480000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[10]*sc[11];
    Kc = exp(g_RT[9] - g_RT[10] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[13]) + (h_RT[10] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[13];
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[209] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[211] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[11];
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[231] += dqdci;              /* dwdot[CO]/d[CH4] */
    J[233] -= dqdci;              /* dwdot[HCO]/d[CH4] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[10];
    J[251] -= dqdci;              /* dwdot[CH3]/d[CO] */
    J[252] += dqdci;              /* dwdot[CH4]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[9];
    J[295] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[296] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 79: CH3 + CH2O <=> HCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[14];
    k_f = 1.0000000000000002e-06 * 3320
                * exp(2.8100000000000001 * tc[0] - 0.50321666580471969 * (5860) * invT);
    dlnkfdT = 2.8100000000000001 * invT + 0.50321666580471969 *  5860  * invT2;
    /* reverse */
    phi_r = sc[10]*sc[13];
    Kc = exp(g_RT[9] - g_RT[10] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[14]) + (h_RT[10] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[14];
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[211] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[212] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[13];
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[233] += dqdci;              /* dwdot[HCO]/d[CH4] */
    J[234] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[10];
    J[295] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[296] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[9];
    J[317] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[318] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 80: CH3 + C2H6 <=> C2H5 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[18];
    k_f = 1.0000000000000002e-06 * 6140000
                * exp(1.74 * tc[0] - 0.50321666580471969 * (10450) * invT);
    dlnkfdT = 1.74 * invT + 0.50321666580471969 *  10450  * invT2;
    /* reverse */
    phi_r = sc[10]*sc[17];
    Kc = exp(g_RT[9] - g_RT[10] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[18]) + (h_RT[10] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[18];
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[215] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[216] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[17];
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[237] += dqdci;              /* dwdot[C2H5]/d[CH4] */
    J[238] -= dqdci;              /* dwdot[C2H6]/d[CH4] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[10];
    J[383] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[384] += dqdci;              /* dwdot[CH4]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[9];
    J[405] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[406] += dqdci;              /* dwdot[CH4]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 81: HCO + H2O <=> H + CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = 1.0000000000000002e-06 * 2.244e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  17000  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5]*sc[11];
    Kc = refC * exp(-g_RT[1] + g_RT[5] - g_RT[5] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[13]) + (h_RT[1] + h_RT[5] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5]*sc[11];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[33] += dqdci;               /* dwdot[CO]/d[H] */
    J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[13] - k_r*sc[1]*sc[11];
    J[111] += dqdci;              /* dwdot[H]/d[H2O] */
    J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[5];
    J[243] += dqdci;              /* dwdot[H]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[13];
    k_f = 1.0000000000000002e-06 * 7600000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (400) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  400  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[11];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[13]) + (h_RT[6] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[CO]/d[O2] */
    J[79] -= dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[11];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[143] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[145] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[248] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[289] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[292] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = 1.0000000000000002e-06 * 4.2799999999999999e-13
                * exp(7.5999999999999996 * tc[0] - 0.50321666580471969 * (-3530) * invT);
    dlnkfdT = 7.5999999999999996 * invT + 0.50321666580471969 *  -3530  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[14];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[14] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[15]) + (h_RT[6] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[14] += q; /* CH2O */
    wdot[15] -= q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[15];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[80] += dqdci;               /* dwdot[CH2O]/d[O2] */
    J[81] -= dqdci;               /* dwdot[CH3O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[14];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[146] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[147] -= dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[6];
    J[311] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[314] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[323] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[3];
    J[333] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[336] += dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[344] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[345] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */
    J[477] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[17];
    k_f = 1.0000000000000002e-06 * 840000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3875) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3875  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[16];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[17]) + (h_RT[6] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[16] += q; /* C2H4 */
    wdot[17] -= q; /* C2H5 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[17];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[82] += dqdci;               /* dwdot[C2H4]/d[O2] */
    J[83] -= dqdci;               /* dwdot[C2H5]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[16];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[148] += dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[149] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[6];
    J[355] -= dqdci;              /* dwdot[O2]/d[C2H4] */
    J[358] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[368] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[369] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[3];
    J[377] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[380] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[390] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[391] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[478] += dqdT;               /* dwdot[C2H4]/dT */
    J[479] -= dqdT;               /* dwdot[C2H5]/dT */

    amrex::Real c_R[21], dcRdT[21], e_RT[21];
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
    for (int k = 0; k < 21; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[462+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 21; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 21; ++m) {
            dehmixdc += eh_RT[m]*J[k*22+m];
        }
        J[k*22+21] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[483] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;

return;
}
#endif


/*compute the reaction Jacobian */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * consP)
{
    amrex::Real c[21];

    for (int k=0; k<21; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<21; k++) {
        J[462+k] *= 1.e-6;
        J[k*22+21] *= 1.e6;
    }

    return;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void sgjsolve_simplified(amrex::Real* A, amrex::Real* x, amrex::Real* b) {
    amrex::Real scratch_0 = 1.0/A[0];
    amrex::Real scratch_1 = b[0]*scratch_0;
    amrex::Real scratch_2 = A[19]*scratch_0;
    amrex::Real scratch_3 = A[2]*scratch_0;
    amrex::Real scratch_4 = A[1]*scratch_0;
    amrex::Real scratch_5 = 1.0/(-A[22]*scratch_4 + A[23]);
    amrex::Real scratch_6 = scratch_5*(-A[22]*scratch_3 + A[24]);
    amrex::Real scratch_7 = scratch_3 - scratch_4*scratch_6;
    amrex::Real scratch_8 = -A[44]*scratch_4 + A[45];
    amrex::Real scratch_9 = 1.0/(-A[44]*scratch_3 + A[46] - scratch_6*scratch_8);
    amrex::Real scratch_10 = scratch_5*(-A[22]*scratch_2 + A[41]);
    amrex::Real scratch_11 = scratch_9*(-A[44]*scratch_2 + A[63] - scratch_10*scratch_8);
    amrex::Real scratch_12 = A[3]*scratch_0;
    amrex::Real scratch_13 = scratch_5*(-A[22]*scratch_12 + A[25]);
    amrex::Real scratch_14 = scratch_9*(-A[44]*scratch_12 + A[47] - scratch_13*scratch_8);
    amrex::Real scratch_15 = scratch_12 - scratch_13*scratch_4 - scratch_14*scratch_7;
    amrex::Real scratch_16 = -A[66]*scratch_4 + A[67];
    amrex::Real scratch_17 = -A[66]*scratch_3 + A[68] - scratch_16*scratch_6;
    amrex::Real scratch_18 = 1.0/(-A[66]*scratch_12 + A[69] - scratch_13*scratch_16 - scratch_14*scratch_17);
    amrex::Real scratch_19 = scratch_18*(-A[66]*scratch_2 + A[85] - scratch_10*scratch_16 - scratch_11*scratch_17);
    amrex::Real scratch_20 = A[4]*scratch_0;
    amrex::Real scratch_21 = scratch_5*(-A[22]*scratch_20 + A[26]);
    amrex::Real scratch_22 = scratch_9*(-A[44]*scratch_20 + A[48] - scratch_21*scratch_8);
    amrex::Real scratch_23 = scratch_18*(-A[66]*scratch_20 + A[70] - scratch_16*scratch_21 - scratch_17*scratch_22);
    amrex::Real scratch_24 = -scratch_15*scratch_23 + scratch_20 - scratch_21*scratch_4 - scratch_22*scratch_7;
    amrex::Real scratch_25 = -A[88]*scratch_4 + A[89];
    amrex::Real scratch_26 = -A[88]*scratch_3 + A[90] - scratch_25*scratch_6;
    amrex::Real scratch_27 = -A[88]*scratch_12 + A[91] - scratch_13*scratch_25 - scratch_14*scratch_26;
    amrex::Real scratch_28 = 1.0/(-A[88]*scratch_20 + A[92] - scratch_21*scratch_25 - scratch_22*scratch_26 - scratch_23*scratch_27);
    amrex::Real scratch_29 = scratch_28*(A[107] - A[88]*scratch_2 - scratch_10*scratch_25 - scratch_11*scratch_26 - scratch_19*scratch_27);
    amrex::Real scratch_30 = A[5]*scratch_0;
    amrex::Real scratch_31 = scratch_5*(-A[22]*scratch_30 + A[27]);
    amrex::Real scratch_32 = scratch_9*(-A[44]*scratch_30 + A[49] - scratch_31*scratch_8);
    amrex::Real scratch_33 = scratch_18*(-A[66]*scratch_30 + A[71] - scratch_16*scratch_31 - scratch_17*scratch_32);
    amrex::Real scratch_34 = scratch_28*(-A[88]*scratch_30 + A[93] - scratch_25*scratch_31 - scratch_26*scratch_32 - scratch_27*scratch_33);
    amrex::Real scratch_35 = -scratch_15*scratch_33 - scratch_24*scratch_34 + scratch_30 - scratch_31*scratch_4 - scratch_32*scratch_7;
    amrex::Real scratch_36 = -A[110]*scratch_4 + A[111];
    amrex::Real scratch_37 = -A[110]*scratch_3 + A[112] - scratch_36*scratch_6;
    amrex::Real scratch_38 = -A[110]*scratch_12 + A[113] - scratch_13*scratch_36 - scratch_14*scratch_37;
    amrex::Real scratch_39 = -A[110]*scratch_20 + A[114] - scratch_21*scratch_36 - scratch_22*scratch_37 - scratch_23*scratch_38;
    amrex::Real scratch_40 = 1.0/(-A[110]*scratch_30 + A[115] - scratch_31*scratch_36 - scratch_32*scratch_37 - scratch_33*scratch_38 - scratch_34*scratch_39);
    amrex::Real scratch_41 = scratch_40*(-A[110]*scratch_2 + A[129] - scratch_10*scratch_36 - scratch_11*scratch_37 - scratch_19*scratch_38 - scratch_29*scratch_39);
    amrex::Real scratch_42 = A[6]*scratch_0;
    amrex::Real scratch_43 = scratch_5*(-A[22]*scratch_42 + A[28]);
    amrex::Real scratch_44 = scratch_9*(-A[44]*scratch_42 + A[50] - scratch_43*scratch_8);
    amrex::Real scratch_45 = scratch_18*(-A[66]*scratch_42 + A[72] - scratch_16*scratch_43 - scratch_17*scratch_44);
    amrex::Real scratch_46 = scratch_28*(-A[88]*scratch_42 + A[94] - scratch_25*scratch_43 - scratch_26*scratch_44 - scratch_27*scratch_45);
    amrex::Real scratch_47 = scratch_40*(-A[110]*scratch_42 + A[116] - scratch_36*scratch_43 - scratch_37*scratch_44 - scratch_38*scratch_45 - scratch_39*scratch_46);
    amrex::Real scratch_48 = -scratch_15*scratch_45 - scratch_24*scratch_46 - scratch_35*scratch_47 - scratch_4*scratch_43 + scratch_42 - scratch_44*scratch_7;
    amrex::Real scratch_49 = -A[132]*scratch_4 + A[133];
    amrex::Real scratch_50 = -A[132]*scratch_3 + A[134] - scratch_49*scratch_6;
    amrex::Real scratch_51 = -A[132]*scratch_12 + A[135] - scratch_13*scratch_49 - scratch_14*scratch_50;
    amrex::Real scratch_52 = -A[132]*scratch_20 + A[136] - scratch_21*scratch_49 - scratch_22*scratch_50 - scratch_23*scratch_51;
    amrex::Real scratch_53 = -A[132]*scratch_30 + A[137] - scratch_31*scratch_49 - scratch_32*scratch_50 - scratch_33*scratch_51 - scratch_34*scratch_52;
    amrex::Real scratch_54 = 1.0/(-A[132]*scratch_42 + A[138] - scratch_43*scratch_49 - scratch_44*scratch_50 - scratch_45*scratch_51 - scratch_46*scratch_52 - scratch_47*scratch_53);
    amrex::Real scratch_55 = scratch_54*(-A[132]*scratch_2 + A[151] - scratch_10*scratch_49 - scratch_11*scratch_50 - scratch_19*scratch_51 - scratch_29*scratch_52 - scratch_41*scratch_53);
    amrex::Real scratch_56 = A[7]*scratch_0;
    amrex::Real scratch_57 = scratch_5*(-A[22]*scratch_56 + A[29]);
    amrex::Real scratch_58 = scratch_9*(-A[44]*scratch_56 + A[51] - scratch_57*scratch_8);
    amrex::Real scratch_59 = scratch_18*(-A[66]*scratch_56 + A[73] - scratch_16*scratch_57 - scratch_17*scratch_58);
    amrex::Real scratch_60 = scratch_28*(-A[88]*scratch_56 + A[95] - scratch_25*scratch_57 - scratch_26*scratch_58 - scratch_27*scratch_59);
    amrex::Real scratch_61 = scratch_40*(-A[110]*scratch_56 + A[117] - scratch_36*scratch_57 - scratch_37*scratch_58 - scratch_38*scratch_59 - scratch_39*scratch_60);
    amrex::Real scratch_62 = scratch_54*(-A[132]*scratch_56 + A[139] - scratch_49*scratch_57 - scratch_50*scratch_58 - scratch_51*scratch_59 - scratch_52*scratch_60 - scratch_53*scratch_61);
    amrex::Real scratch_63 = -scratch_15*scratch_59 - scratch_24*scratch_60 - scratch_35*scratch_61 - scratch_4*scratch_57 - scratch_48*scratch_62 + scratch_56 - scratch_58*scratch_7;
    amrex::Real scratch_64 = -A[154]*scratch_4 + A[155];
    amrex::Real scratch_65 = -A[154]*scratch_3 + A[156] - scratch_6*scratch_64;
    amrex::Real scratch_66 = -A[154]*scratch_12 + A[157] - scratch_13*scratch_64 - scratch_14*scratch_65;
    amrex::Real scratch_67 = -A[154]*scratch_20 + A[158] - scratch_21*scratch_64 - scratch_22*scratch_65 - scratch_23*scratch_66;
    amrex::Real scratch_68 = -A[154]*scratch_30 + A[159] - scratch_31*scratch_64 - scratch_32*scratch_65 - scratch_33*scratch_66 - scratch_34*scratch_67;
    amrex::Real scratch_69 = -A[154]*scratch_42 + A[160] - scratch_43*scratch_64 - scratch_44*scratch_65 - scratch_45*scratch_66 - scratch_46*scratch_67 - scratch_47*scratch_68;
    amrex::Real scratch_70 = 1.0/(-A[154]*scratch_56 + A[161] - scratch_57*scratch_64 - scratch_58*scratch_65 - scratch_59*scratch_66 - scratch_60*scratch_67 - scratch_61*scratch_68 - scratch_62*scratch_69);
    amrex::Real scratch_71 = scratch_70*(-A[154]*scratch_2 + A[170] - scratch_10*scratch_64 - scratch_11*scratch_65 - scratch_19*scratch_66 - scratch_29*scratch_67 - scratch_41*scratch_68 - scratch_55*scratch_69);
    amrex::Real scratch_72 = A[8]*scratch_0;
    amrex::Real scratch_73 = scratch_5*(-A[22]*scratch_72 + A[30]);
    amrex::Real scratch_74 = scratch_9*(-A[44]*scratch_72 + A[52] - scratch_73*scratch_8);
    amrex::Real scratch_75 = scratch_18*(-A[66]*scratch_72 + A[74] - scratch_16*scratch_73 - scratch_17*scratch_74);
    amrex::Real scratch_76 = scratch_28*(-A[88]*scratch_72 + A[96] - scratch_25*scratch_73 - scratch_26*scratch_74 - scratch_27*scratch_75);
    amrex::Real scratch_77 = scratch_40*(-A[110]*scratch_72 + A[118] - scratch_36*scratch_73 - scratch_37*scratch_74 - scratch_38*scratch_75 - scratch_39*scratch_76);
    amrex::Real scratch_78 = scratch_54*(-A[132]*scratch_72 + A[140] - scratch_49*scratch_73 - scratch_50*scratch_74 - scratch_51*scratch_75 - scratch_52*scratch_76 - scratch_53*scratch_77);
    amrex::Real scratch_79 = scratch_70*(-A[154]*scratch_72 + A[162] - scratch_64*scratch_73 - scratch_65*scratch_74 - scratch_66*scratch_75 - scratch_67*scratch_76 - scratch_68*scratch_77 - scratch_69*scratch_78);
    amrex::Real scratch_80 = -scratch_15*scratch_75 - scratch_24*scratch_76 - scratch_35*scratch_77 - scratch_4*scratch_73 - scratch_48*scratch_78 - scratch_63*scratch_79 - scratch_7*scratch_74 + scratch_72;
    amrex::Real scratch_81 = -A[173]*scratch_4 + A[174];
    amrex::Real scratch_82 = -A[173]*scratch_3 + A[175] - scratch_6*scratch_81;
    amrex::Real scratch_83 = -A[173]*scratch_12 + A[176] - scratch_13*scratch_81 - scratch_14*scratch_82;
    amrex::Real scratch_84 = -A[173]*scratch_20 + A[177] - scratch_21*scratch_81 - scratch_22*scratch_82 - scratch_23*scratch_83;
    amrex::Real scratch_85 = -A[173]*scratch_30 + A[178] - scratch_31*scratch_81 - scratch_32*scratch_82 - scratch_33*scratch_83 - scratch_34*scratch_84;
    amrex::Real scratch_86 = -A[173]*scratch_42 - scratch_43*scratch_81 - scratch_44*scratch_82 - scratch_45*scratch_83 - scratch_46*scratch_84 - scratch_47*scratch_85;
    amrex::Real scratch_87 = -A[173]*scratch_56 + A[179] - scratch_57*scratch_81 - scratch_58*scratch_82 - scratch_59*scratch_83 - scratch_60*scratch_84 - scratch_61*scratch_85 - scratch_62*scratch_86;
    amrex::Real scratch_88 = 1.0/(-A[173]*scratch_72 + A[180] - scratch_73*scratch_81 - scratch_74*scratch_82 - scratch_75*scratch_83 - scratch_76*scratch_84 - scratch_77*scratch_85 - scratch_78*scratch_86 - scratch_79*scratch_87);
    amrex::Real scratch_89 = scratch_88*(-A[173]*scratch_2 + A[188] - scratch_10*scratch_81 - scratch_11*scratch_82 - scratch_19*scratch_83 - scratch_29*scratch_84 - scratch_41*scratch_85 - scratch_55*scratch_86 - scratch_71*scratch_87);
    amrex::Real scratch_90 = A[9]*scratch_0;
    amrex::Real scratch_91 = scratch_5*(-A[22]*scratch_90 + A[31]);
    amrex::Real scratch_92 = scratch_9*(-A[44]*scratch_90 + A[53] - scratch_8*scratch_91);
    amrex::Real scratch_93 = scratch_18*(-A[66]*scratch_90 + A[75] - scratch_16*scratch_91 - scratch_17*scratch_92);
    amrex::Real scratch_94 = scratch_28*(-A[88]*scratch_90 + A[97] - scratch_25*scratch_91 - scratch_26*scratch_92 - scratch_27*scratch_93);
    amrex::Real scratch_95 = scratch_40*(-A[110]*scratch_90 + A[119] - scratch_36*scratch_91 - scratch_37*scratch_92 - scratch_38*scratch_93 - scratch_39*scratch_94);
    amrex::Real scratch_96 = scratch_54*(-A[132]*scratch_90 + A[141] - scratch_49*scratch_91 - scratch_50*scratch_92 - scratch_51*scratch_93 - scratch_52*scratch_94 - scratch_53*scratch_95);
    amrex::Real scratch_97 = scratch_70*(-A[154]*scratch_90 + A[163] - scratch_64*scratch_91 - scratch_65*scratch_92 - scratch_66*scratch_93 - scratch_67*scratch_94 - scratch_68*scratch_95 - scratch_69*scratch_96);
    amrex::Real scratch_98 = scratch_88*(-A[173]*scratch_90 + A[181] - scratch_81*scratch_91 - scratch_82*scratch_92 - scratch_83*scratch_93 - scratch_84*scratch_94 - scratch_85*scratch_95 - scratch_86*scratch_96 - scratch_87*scratch_97);
    amrex::Real scratch_99 = -scratch_15*scratch_93 - scratch_24*scratch_94 - scratch_35*scratch_95 - scratch_4*scratch_91 - scratch_48*scratch_96 - scratch_63*scratch_97 - scratch_7*scratch_92 - scratch_80*scratch_98 + scratch_90;
    amrex::Real scratch_100 = -A[191]*scratch_4 + A[192];
    amrex::Real scratch_101 = -A[191]*scratch_3 + A[193] - scratch_100*scratch_6;
    amrex::Real scratch_102 = -A[191]*scratch_12 + A[194] - scratch_100*scratch_13 - scratch_101*scratch_14;
    amrex::Real scratch_103 = -A[191]*scratch_20 + A[195] - scratch_100*scratch_21 - scratch_101*scratch_22 - scratch_102*scratch_23;
    amrex::Real scratch_104 = -A[191]*scratch_30 + A[196] - scratch_100*scratch_31 - scratch_101*scratch_32 - scratch_102*scratch_33 - scratch_103*scratch_34;
    amrex::Real scratch_105 = -A[191]*scratch_42 + A[197] - scratch_100*scratch_43 - scratch_101*scratch_44 - scratch_102*scratch_45 - scratch_103*scratch_46 - scratch_104*scratch_47;
    amrex::Real scratch_106 = -A[191]*scratch_56 + A[198] - scratch_100*scratch_57 - scratch_101*scratch_58 - scratch_102*scratch_59 - scratch_103*scratch_60 - scratch_104*scratch_61 - scratch_105*scratch_62;
    amrex::Real scratch_107 = -A[191]*scratch_72 + A[199] - scratch_100*scratch_73 - scratch_101*scratch_74 - scratch_102*scratch_75 - scratch_103*scratch_76 - scratch_104*scratch_77 - scratch_105*scratch_78 - scratch_106*scratch_79;
    amrex::Real scratch_108 = 1.0/(-A[191]*scratch_90 + A[200] - scratch_100*scratch_91 - scratch_101*scratch_92 - scratch_102*scratch_93 - scratch_103*scratch_94 - scratch_104*scratch_95 - scratch_105*scratch_96 - scratch_106*scratch_97 - scratch_107*scratch_98);
    amrex::Real scratch_109 = scratch_108*(-A[191]*scratch_2 - scratch_10*scratch_100 - scratch_101*scratch_11 - scratch_102*scratch_19 - scratch_103*scratch_29 - scratch_104*scratch_41 - scratch_105*scratch_55 - scratch_106*scratch_71 - scratch_107*scratch_89);
    amrex::Real scratch_110 = A[10]*scratch_0;
    amrex::Real scratch_111 = scratch_5*(-A[22]*scratch_110 + A[32]);
    amrex::Real scratch_112 = scratch_9*(-A[44]*scratch_110 + A[54] - scratch_111*scratch_8);
    amrex::Real scratch_113 = scratch_18*(-A[66]*scratch_110 + A[76] - scratch_111*scratch_16 - scratch_112*scratch_17);
    amrex::Real scratch_114 = scratch_28*(-A[88]*scratch_110 + A[98] - scratch_111*scratch_25 - scratch_112*scratch_26 - scratch_113*scratch_27);
    amrex::Real scratch_115 = scratch_40*(-A[110]*scratch_110 + A[120] - scratch_111*scratch_36 - scratch_112*scratch_37 - scratch_113*scratch_38 - scratch_114*scratch_39);
    amrex::Real scratch_116 = scratch_54*(-A[132]*scratch_110 + A[142] - scratch_111*scratch_49 - scratch_112*scratch_50 - scratch_113*scratch_51 - scratch_114*scratch_52 - scratch_115*scratch_53);
    amrex::Real scratch_117 = scratch_70*(-A[154]*scratch_110 + A[164] - scratch_111*scratch_64 - scratch_112*scratch_65 - scratch_113*scratch_66 - scratch_114*scratch_67 - scratch_115*scratch_68 - scratch_116*scratch_69);
    amrex::Real scratch_118 = scratch_88*(-A[173]*scratch_110 + A[182] - scratch_111*scratch_81 - scratch_112*scratch_82 - scratch_113*scratch_83 - scratch_114*scratch_84 - scratch_115*scratch_85 - scratch_116*scratch_86 - scratch_117*scratch_87);
    amrex::Real scratch_119 = scratch_108*(-A[191]*scratch_110 + A[201] - scratch_100*scratch_111 - scratch_101*scratch_112 - scratch_102*scratch_113 - scratch_103*scratch_114 - scratch_104*scratch_115 - scratch_105*scratch_116 - scratch_106*scratch_117 - scratch_107*scratch_118);
    amrex::Real scratch_120 = scratch_110 - scratch_111*scratch_4 - scratch_112*scratch_7 - scratch_113*scratch_15 - scratch_114*scratch_24 - scratch_115*scratch_35 - scratch_116*scratch_48 - scratch_117*scratch_63 - scratch_118*scratch_80 - scratch_119*scratch_99;
    amrex::Real scratch_121 = -A[210]*scratch_4 + A[211];
    amrex::Real scratch_122 = -A[210]*scratch_3 + A[212] - scratch_121*scratch_6;
    amrex::Real scratch_123 = -A[210]*scratch_12 + A[213] - scratch_121*scratch_13 - scratch_122*scratch_14;
    amrex::Real scratch_124 = -A[210]*scratch_20 + A[214] - scratch_121*scratch_21 - scratch_122*scratch_22 - scratch_123*scratch_23;
    amrex::Real scratch_125 = -A[210]*scratch_30 + A[215] - scratch_121*scratch_31 - scratch_122*scratch_32 - scratch_123*scratch_33 - scratch_124*scratch_34;
    amrex::Real scratch_126 = -A[210]*scratch_42 + A[216] - scratch_121*scratch_43 - scratch_122*scratch_44 - scratch_123*scratch_45 - scratch_124*scratch_46 - scratch_125*scratch_47;
    amrex::Real scratch_127 = -A[210]*scratch_56 + A[217] - scratch_121*scratch_57 - scratch_122*scratch_58 - scratch_123*scratch_59 - scratch_124*scratch_60 - scratch_125*scratch_61 - scratch_126*scratch_62;
    amrex::Real scratch_128 = -A[210]*scratch_72 + A[218] - scratch_121*scratch_73 - scratch_122*scratch_74 - scratch_123*scratch_75 - scratch_124*scratch_76 - scratch_125*scratch_77 - scratch_126*scratch_78 - scratch_127*scratch_79;
    amrex::Real scratch_129 = -A[210]*scratch_90 + A[219] - scratch_121*scratch_91 - scratch_122*scratch_92 - scratch_123*scratch_93 - scratch_124*scratch_94 - scratch_125*scratch_95 - scratch_126*scratch_96 - scratch_127*scratch_97 - scratch_128*scratch_98;
    amrex::Real scratch_130 = 1.0/(-A[210]*scratch_110 + A[220] - scratch_111*scratch_121 - scratch_112*scratch_122 - scratch_113*scratch_123 - scratch_114*scratch_124 - scratch_115*scratch_125 - scratch_116*scratch_126 - scratch_117*scratch_127 - scratch_118*scratch_128 - scratch_119*scratch_129);
    amrex::Real scratch_131 = scratch_130*(-A[210]*scratch_2 - scratch_10*scratch_121 - scratch_109*scratch_129 - scratch_11*scratch_122 - scratch_123*scratch_19 - scratch_124*scratch_29 - scratch_125*scratch_41 - scratch_126*scratch_55 - scratch_127*scratch_71 - scratch_128*scratch_89);
    amrex::Real scratch_132 = A[11]*scratch_0;
    amrex::Real scratch_133 = scratch_5*(-A[22]*scratch_132 + A[33]);
    amrex::Real scratch_134 = scratch_9*(-A[44]*scratch_132 + A[55] - scratch_133*scratch_8);
    amrex::Real scratch_135 = scratch_18*(-A[66]*scratch_132 + A[77] - scratch_133*scratch_16 - scratch_134*scratch_17);
    amrex::Real scratch_136 = scratch_28*(A[99] - A[88]*scratch_132 - scratch_133*scratch_25 - scratch_134*scratch_26 - scratch_135*scratch_27);
    amrex::Real scratch_137 = scratch_40*(-A[110]*scratch_132 + A[121] - scratch_133*scratch_36 - scratch_134*scratch_37 - scratch_135*scratch_38 - scratch_136*scratch_39);
    amrex::Real scratch_138 = scratch_54*(-A[132]*scratch_132 + A[143] - scratch_133*scratch_49 - scratch_134*scratch_50 - scratch_135*scratch_51 - scratch_136*scratch_52 - scratch_137*scratch_53);
    amrex::Real scratch_139 = scratch_70*(-A[154]*scratch_132 + A[165] - scratch_133*scratch_64 - scratch_134*scratch_65 - scratch_135*scratch_66 - scratch_136*scratch_67 - scratch_137*scratch_68 - scratch_138*scratch_69);
    amrex::Real scratch_140 = scratch_88*(-A[173]*scratch_132 + A[183] - scratch_133*scratch_81 - scratch_134*scratch_82 - scratch_135*scratch_83 - scratch_136*scratch_84 - scratch_137*scratch_85 - scratch_138*scratch_86 - scratch_139*scratch_87);
    amrex::Real scratch_141 = scratch_108*(-A[191]*scratch_132 + A[202] - scratch_100*scratch_133 - scratch_101*scratch_134 - scratch_102*scratch_135 - scratch_103*scratch_136 - scratch_104*scratch_137 - scratch_105*scratch_138 - scratch_106*scratch_139 - scratch_107*scratch_140);
    amrex::Real scratch_142 = scratch_130*(-A[210]*scratch_132 + A[221] - scratch_121*scratch_133 - scratch_122*scratch_134 - scratch_123*scratch_135 - scratch_124*scratch_136 - scratch_125*scratch_137 - scratch_126*scratch_138 - scratch_127*scratch_139 - scratch_128*scratch_140 - scratch_129*scratch_141);
    amrex::Real scratch_143 = -scratch_120*scratch_142 + scratch_132 - scratch_133*scratch_4 - scratch_134*scratch_7 - scratch_135*scratch_15 - scratch_136*scratch_24 - scratch_137*scratch_35 - scratch_138*scratch_48 - scratch_139*scratch_63 - scratch_140*scratch_80 - scratch_141*scratch_99;
    amrex::Real scratch_144 = -A[227]*scratch_4 + A[228];
    amrex::Real scratch_145 = -A[227]*scratch_3 + A[229] - scratch_144*scratch_6;
    amrex::Real scratch_146 = -A[227]*scratch_12 + A[230] - scratch_13*scratch_144 - scratch_14*scratch_145;
    amrex::Real scratch_147 = -A[227]*scratch_20 + A[231] - scratch_144*scratch_21 - scratch_145*scratch_22 - scratch_146*scratch_23;
    amrex::Real scratch_148 = -A[227]*scratch_30 + A[232] - scratch_144*scratch_31 - scratch_145*scratch_32 - scratch_146*scratch_33 - scratch_147*scratch_34;
    amrex::Real scratch_149 = -A[227]*scratch_42 + A[233] - scratch_144*scratch_43 - scratch_145*scratch_44 - scratch_146*scratch_45 - scratch_147*scratch_46 - scratch_148*scratch_47;
    amrex::Real scratch_150 = -A[227]*scratch_56 + A[234] - scratch_144*scratch_57 - scratch_145*scratch_58 - scratch_146*scratch_59 - scratch_147*scratch_60 - scratch_148*scratch_61 - scratch_149*scratch_62;
    amrex::Real scratch_151 = -A[227]*scratch_72 + A[235] - scratch_144*scratch_73 - scratch_145*scratch_74 - scratch_146*scratch_75 - scratch_147*scratch_76 - scratch_148*scratch_77 - scratch_149*scratch_78 - scratch_150*scratch_79;
    amrex::Real scratch_152 = -A[227]*scratch_90 + A[236] - scratch_144*scratch_91 - scratch_145*scratch_92 - scratch_146*scratch_93 - scratch_147*scratch_94 - scratch_148*scratch_95 - scratch_149*scratch_96 - scratch_150*scratch_97 - scratch_151*scratch_98;
    amrex::Real scratch_153 = -A[227]*scratch_110 + A[237] - scratch_111*scratch_144 - scratch_112*scratch_145 - scratch_113*scratch_146 - scratch_114*scratch_147 - scratch_115*scratch_148 - scratch_116*scratch_149 - scratch_117*scratch_150 - scratch_118*scratch_151 - scratch_119*scratch_152;
    amrex::Real scratch_154 = 1.0/(-A[227]*scratch_132 + A[238] - scratch_133*scratch_144 - scratch_134*scratch_145 - scratch_135*scratch_146 - scratch_136*scratch_147 - scratch_137*scratch_148 - scratch_138*scratch_149 - scratch_139*scratch_150 - scratch_140*scratch_151 - scratch_141*scratch_152 - scratch_142*scratch_153);
    amrex::Real scratch_155 = scratch_154*(-A[227]*scratch_2 + A[246] - scratch_10*scratch_144 - scratch_109*scratch_152 - scratch_11*scratch_145 - scratch_131*scratch_153 - scratch_146*scratch_19 - scratch_147*scratch_29 - scratch_148*scratch_41 - scratch_149*scratch_55 - scratch_150*scratch_71 - scratch_151*scratch_89);
    amrex::Real scratch_156 = A[12]*scratch_0;
    amrex::Real scratch_157 = scratch_5*(-A[22]*scratch_156 + A[34]);
    amrex::Real scratch_158 = scratch_9*(-A[44]*scratch_156 + A[56] - scratch_157*scratch_8);
    amrex::Real scratch_159 = scratch_18*(-A[66]*scratch_156 + A[78] - scratch_157*scratch_16 - scratch_158*scratch_17);
    amrex::Real scratch_160 = scratch_28*(A[100] - A[88]*scratch_156 - scratch_157*scratch_25 - scratch_158*scratch_26 - scratch_159*scratch_27);
    amrex::Real scratch_161 = scratch_40*(-A[110]*scratch_156 + A[122] - scratch_157*scratch_36 - scratch_158*scratch_37 - scratch_159*scratch_38 - scratch_160*scratch_39);
    amrex::Real scratch_162 = scratch_54*(-A[132]*scratch_156 + A[144] - scratch_157*scratch_49 - scratch_158*scratch_50 - scratch_159*scratch_51 - scratch_160*scratch_52 - scratch_161*scratch_53);
    amrex::Real scratch_163 = scratch_70*(-A[154]*scratch_156 + A[166] - scratch_157*scratch_64 - scratch_158*scratch_65 - scratch_159*scratch_66 - scratch_160*scratch_67 - scratch_161*scratch_68 - scratch_162*scratch_69);
    amrex::Real scratch_164 = scratch_88*(-A[173]*scratch_156 + A[184] - scratch_157*scratch_81 - scratch_158*scratch_82 - scratch_159*scratch_83 - scratch_160*scratch_84 - scratch_161*scratch_85 - scratch_162*scratch_86 - scratch_163*scratch_87);
    amrex::Real scratch_165 = scratch_108*(-A[191]*scratch_156 - scratch_100*scratch_157 - scratch_101*scratch_158 - scratch_102*scratch_159 - scratch_103*scratch_160 - scratch_104*scratch_161 - scratch_105*scratch_162 - scratch_106*scratch_163 - scratch_107*scratch_164);
    amrex::Real scratch_166 = scratch_130*(-A[210]*scratch_156 - scratch_121*scratch_157 - scratch_122*scratch_158 - scratch_123*scratch_159 - scratch_124*scratch_160 - scratch_125*scratch_161 - scratch_126*scratch_162 - scratch_127*scratch_163 - scratch_128*scratch_164 - scratch_129*scratch_165);
    amrex::Real scratch_167 = scratch_154*(-A[227]*scratch_156 + A[239] - scratch_144*scratch_157 - scratch_145*scratch_158 - scratch_146*scratch_159 - scratch_147*scratch_160 - scratch_148*scratch_161 - scratch_149*scratch_162 - scratch_150*scratch_163 - scratch_151*scratch_164 - scratch_152*scratch_165 - scratch_153*scratch_166);
    amrex::Real scratch_168 = -scratch_120*scratch_166 - scratch_143*scratch_167 - scratch_15*scratch_159 + scratch_156 - scratch_157*scratch_4 - scratch_158*scratch_7 - scratch_160*scratch_24 - scratch_161*scratch_35 - scratch_162*scratch_48 - scratch_163*scratch_63 - scratch_164*scratch_80 - scratch_165*scratch_99;
    amrex::Real scratch_169 = -A[249]*scratch_4 + A[250];
    amrex::Real scratch_170 = -A[249]*scratch_3 + A[251] - scratch_169*scratch_6;
    amrex::Real scratch_171 = -A[249]*scratch_12 + A[252] - scratch_13*scratch_169 - scratch_14*scratch_170;
    amrex::Real scratch_172 = -A[249]*scratch_20 + A[253] - scratch_169*scratch_21 - scratch_170*scratch_22 - scratch_171*scratch_23;
    amrex::Real scratch_173 = -A[249]*scratch_30 + A[254] - scratch_169*scratch_31 - scratch_170*scratch_32 - scratch_171*scratch_33 - scratch_172*scratch_34;
    amrex::Real scratch_174 = -A[249]*scratch_42 + A[255] - scratch_169*scratch_43 - scratch_170*scratch_44 - scratch_171*scratch_45 - scratch_172*scratch_46 - scratch_173*scratch_47;
    amrex::Real scratch_175 = -A[249]*scratch_56 + A[256] - scratch_169*scratch_57 - scratch_170*scratch_58 - scratch_171*scratch_59 - scratch_172*scratch_60 - scratch_173*scratch_61 - scratch_174*scratch_62;
    amrex::Real scratch_176 = -A[249]*scratch_72 + A[257] - scratch_169*scratch_73 - scratch_170*scratch_74 - scratch_171*scratch_75 - scratch_172*scratch_76 - scratch_173*scratch_77 - scratch_174*scratch_78 - scratch_175*scratch_79;
    amrex::Real scratch_177 = -A[249]*scratch_90 + A[258] - scratch_169*scratch_91 - scratch_170*scratch_92 - scratch_171*scratch_93 - scratch_172*scratch_94 - scratch_173*scratch_95 - scratch_174*scratch_96 - scratch_175*scratch_97 - scratch_176*scratch_98;
    amrex::Real scratch_178 = -A[249]*scratch_110 + A[259] - scratch_111*scratch_169 - scratch_112*scratch_170 - scratch_113*scratch_171 - scratch_114*scratch_172 - scratch_115*scratch_173 - scratch_116*scratch_174 - scratch_117*scratch_175 - scratch_118*scratch_176 - scratch_119*scratch_177;
    amrex::Real scratch_179 = -A[249]*scratch_132 + A[260] - scratch_133*scratch_169 - scratch_134*scratch_170 - scratch_135*scratch_171 - scratch_136*scratch_172 - scratch_137*scratch_173 - scratch_138*scratch_174 - scratch_139*scratch_175 - scratch_140*scratch_176 - scratch_141*scratch_177 - scratch_142*scratch_178;
    amrex::Real scratch_180 = 1.0/(-A[249]*scratch_156 + A[261] - scratch_157*scratch_169 - scratch_158*scratch_170 - scratch_159*scratch_171 - scratch_160*scratch_172 - scratch_161*scratch_173 - scratch_162*scratch_174 - scratch_163*scratch_175 - scratch_164*scratch_176 - scratch_165*scratch_177 - scratch_166*scratch_178 - scratch_167*scratch_179);
    amrex::Real scratch_181 = scratch_180*(-A[249]*scratch_2 + A[268] - scratch_10*scratch_169 - scratch_109*scratch_177 - scratch_11*scratch_170 - scratch_131*scratch_178 - scratch_155*scratch_179 - scratch_171*scratch_19 - scratch_172*scratch_29 - scratch_173*scratch_41 - scratch_174*scratch_55 - scratch_175*scratch_71 - scratch_176*scratch_89);
    amrex::Real scratch_182 = A[13]*scratch_0;
    amrex::Real scratch_183 = scratch_5*(-A[22]*scratch_182 + A[35]);
    amrex::Real scratch_184 = scratch_9*(-A[44]*scratch_182 + A[57] - scratch_183*scratch_8);
    amrex::Real scratch_185 = scratch_18*(-A[66]*scratch_182 + A[79] - scratch_16*scratch_183 - scratch_17*scratch_184);
    amrex::Real scratch_186 = scratch_28*(A[101] - A[88]*scratch_182 - scratch_183*scratch_25 - scratch_184*scratch_26 - scratch_185*scratch_27);
    amrex::Real scratch_187 = scratch_40*(-A[110]*scratch_182 + A[123] - scratch_183*scratch_36 - scratch_184*scratch_37 - scratch_185*scratch_38 - scratch_186*scratch_39);
    amrex::Real scratch_188 = scratch_54*(-A[132]*scratch_182 + A[145] - scratch_183*scratch_49 - scratch_184*scratch_50 - scratch_185*scratch_51 - scratch_186*scratch_52 - scratch_187*scratch_53);
    amrex::Real scratch_189 = scratch_70*(-A[154]*scratch_182 + A[167] - scratch_183*scratch_64 - scratch_184*scratch_65 - scratch_185*scratch_66 - scratch_186*scratch_67 - scratch_187*scratch_68 - scratch_188*scratch_69);
    amrex::Real scratch_190 = scratch_88*(-A[173]*scratch_182 + A[185] - scratch_183*scratch_81 - scratch_184*scratch_82 - scratch_185*scratch_83 - scratch_186*scratch_84 - scratch_187*scratch_85 - scratch_188*scratch_86 - scratch_189*scratch_87);
    amrex::Real scratch_191 = scratch_108*(-A[191]*scratch_182 + A[203] - scratch_100*scratch_183 - scratch_101*scratch_184 - scratch_102*scratch_185 - scratch_103*scratch_186 - scratch_104*scratch_187 - scratch_105*scratch_188 - scratch_106*scratch_189 - scratch_107*scratch_190);
    amrex::Real scratch_192 = scratch_130*(-A[210]*scratch_182 + A[222] - scratch_121*scratch_183 - scratch_122*scratch_184 - scratch_123*scratch_185 - scratch_124*scratch_186 - scratch_125*scratch_187 - scratch_126*scratch_188 - scratch_127*scratch_189 - scratch_128*scratch_190 - scratch_129*scratch_191);
    amrex::Real scratch_193 = scratch_154*(-A[227]*scratch_182 + A[240] - scratch_144*scratch_183 - scratch_145*scratch_184 - scratch_146*scratch_185 - scratch_147*scratch_186 - scratch_148*scratch_187 - scratch_149*scratch_188 - scratch_150*scratch_189 - scratch_151*scratch_190 - scratch_152*scratch_191 - scratch_153*scratch_192);
    amrex::Real scratch_194 = scratch_180*(-A[249]*scratch_182 + A[262] - scratch_169*scratch_183 - scratch_170*scratch_184 - scratch_171*scratch_185 - scratch_172*scratch_186 - scratch_173*scratch_187 - scratch_174*scratch_188 - scratch_175*scratch_189 - scratch_176*scratch_190 - scratch_177*scratch_191 - scratch_178*scratch_192 - scratch_179*scratch_193);
    amrex::Real scratch_195 = -scratch_120*scratch_192 - scratch_143*scratch_193 - scratch_15*scratch_185 - scratch_168*scratch_194 + scratch_182 - scratch_183*scratch_4 - scratch_184*scratch_7 - scratch_186*scratch_24 - scratch_187*scratch_35 - scratch_188*scratch_48 - scratch_189*scratch_63 - scratch_190*scratch_80 - scratch_191*scratch_99;
    amrex::Real scratch_196 = -A[271]*scratch_4 + A[272];
    amrex::Real scratch_197 = -A[271]*scratch_3 + A[273] - scratch_196*scratch_6;
    amrex::Real scratch_198 = -A[271]*scratch_12 + A[274] - scratch_13*scratch_196 - scratch_14*scratch_197;
    amrex::Real scratch_199 = -A[271]*scratch_20 + A[275] - scratch_196*scratch_21 - scratch_197*scratch_22 - scratch_198*scratch_23;
    amrex::Real scratch_200 = -A[271]*scratch_30 + A[276] - scratch_196*scratch_31 - scratch_197*scratch_32 - scratch_198*scratch_33 - scratch_199*scratch_34;
    amrex::Real scratch_201 = -A[271]*scratch_42 + A[277] - scratch_196*scratch_43 - scratch_197*scratch_44 - scratch_198*scratch_45 - scratch_199*scratch_46 - scratch_200*scratch_47;
    amrex::Real scratch_202 = -A[271]*scratch_56 + A[278] - scratch_196*scratch_57 - scratch_197*scratch_58 - scratch_198*scratch_59 - scratch_199*scratch_60 - scratch_200*scratch_61 - scratch_201*scratch_62;
    amrex::Real scratch_203 = -A[271]*scratch_72 + A[279] - scratch_196*scratch_73 - scratch_197*scratch_74 - scratch_198*scratch_75 - scratch_199*scratch_76 - scratch_200*scratch_77 - scratch_201*scratch_78 - scratch_202*scratch_79;
    amrex::Real scratch_204 = -A[271]*scratch_90 + A[280] - scratch_196*scratch_91 - scratch_197*scratch_92 - scratch_198*scratch_93 - scratch_199*scratch_94 - scratch_200*scratch_95 - scratch_201*scratch_96 - scratch_202*scratch_97 - scratch_203*scratch_98;
    amrex::Real scratch_205 = -A[271]*scratch_110 + A[281] - scratch_111*scratch_196 - scratch_112*scratch_197 - scratch_113*scratch_198 - scratch_114*scratch_199 - scratch_115*scratch_200 - scratch_116*scratch_201 - scratch_117*scratch_202 - scratch_118*scratch_203 - scratch_119*scratch_204;
    amrex::Real scratch_206 = -A[271]*scratch_132 + A[282] - scratch_133*scratch_196 - scratch_134*scratch_197 - scratch_135*scratch_198 - scratch_136*scratch_199 - scratch_137*scratch_200 - scratch_138*scratch_201 - scratch_139*scratch_202 - scratch_140*scratch_203 - scratch_141*scratch_204 - scratch_142*scratch_205;
    amrex::Real scratch_207 = -A[271]*scratch_156 + A[283] - scratch_157*scratch_196 - scratch_158*scratch_197 - scratch_159*scratch_198 - scratch_160*scratch_199 - scratch_161*scratch_200 - scratch_162*scratch_201 - scratch_163*scratch_202 - scratch_164*scratch_203 - scratch_165*scratch_204 - scratch_166*scratch_205 - scratch_167*scratch_206;
    amrex::Real scratch_208 = 1.0/(-A[271]*scratch_182 + A[284] - scratch_183*scratch_196 - scratch_184*scratch_197 - scratch_185*scratch_198 - scratch_186*scratch_199 - scratch_187*scratch_200 - scratch_188*scratch_201 - scratch_189*scratch_202 - scratch_190*scratch_203 - scratch_191*scratch_204 - scratch_192*scratch_205 - scratch_193*scratch_206 - scratch_194*scratch_207);
    amrex::Real scratch_209 = scratch_208*(-A[271]*scratch_2 + A[290] - scratch_10*scratch_196 - scratch_109*scratch_204 - scratch_11*scratch_197 - scratch_131*scratch_205 - scratch_155*scratch_206 - scratch_181*scratch_207 - scratch_19*scratch_198 - scratch_199*scratch_29 - scratch_200*scratch_41 - scratch_201*scratch_55 - scratch_202*scratch_71 - scratch_203*scratch_89);
    amrex::Real scratch_210 = A[14]*scratch_0;
    amrex::Real scratch_211 = scratch_5*(-A[22]*scratch_210 + A[36]);
    amrex::Real scratch_212 = scratch_9*(-A[44]*scratch_210 + A[58] - scratch_211*scratch_8);
    amrex::Real scratch_213 = scratch_18*(-A[66]*scratch_210 + A[80] - scratch_16*scratch_211 - scratch_17*scratch_212);
    amrex::Real scratch_214 = scratch_28*(A[102] - A[88]*scratch_210 - scratch_211*scratch_25 - scratch_212*scratch_26 - scratch_213*scratch_27);
    amrex::Real scratch_215 = scratch_40*(-A[110]*scratch_210 + A[124] - scratch_211*scratch_36 - scratch_212*scratch_37 - scratch_213*scratch_38 - scratch_214*scratch_39);
    amrex::Real scratch_216 = scratch_54*(-A[132]*scratch_210 + A[146] - scratch_211*scratch_49 - scratch_212*scratch_50 - scratch_213*scratch_51 - scratch_214*scratch_52 - scratch_215*scratch_53);
    amrex::Real scratch_217 = scratch_70*(-A[154]*scratch_210 + A[168] - scratch_211*scratch_64 - scratch_212*scratch_65 - scratch_213*scratch_66 - scratch_214*scratch_67 - scratch_215*scratch_68 - scratch_216*scratch_69);
    amrex::Real scratch_218 = scratch_88*(-A[173]*scratch_210 + A[186] - scratch_211*scratch_81 - scratch_212*scratch_82 - scratch_213*scratch_83 - scratch_214*scratch_84 - scratch_215*scratch_85 - scratch_216*scratch_86 - scratch_217*scratch_87);
    amrex::Real scratch_219 = scratch_108*(-A[191]*scratch_210 + A[204] - scratch_100*scratch_211 - scratch_101*scratch_212 - scratch_102*scratch_213 - scratch_103*scratch_214 - scratch_104*scratch_215 - scratch_105*scratch_216 - scratch_106*scratch_217 - scratch_107*scratch_218);
    amrex::Real scratch_220 = scratch_130*(-A[210]*scratch_210 + A[223] - scratch_121*scratch_211 - scratch_122*scratch_212 - scratch_123*scratch_213 - scratch_124*scratch_214 - scratch_125*scratch_215 - scratch_126*scratch_216 - scratch_127*scratch_217 - scratch_128*scratch_218 - scratch_129*scratch_219);
    amrex::Real scratch_221 = scratch_154*(-A[227]*scratch_210 + A[241] - scratch_144*scratch_211 - scratch_145*scratch_212 - scratch_146*scratch_213 - scratch_147*scratch_214 - scratch_148*scratch_215 - scratch_149*scratch_216 - scratch_150*scratch_217 - scratch_151*scratch_218 - scratch_152*scratch_219 - scratch_153*scratch_220);
    amrex::Real scratch_222 = scratch_180*(-A[249]*scratch_210 + A[263] - scratch_169*scratch_211 - scratch_170*scratch_212 - scratch_171*scratch_213 - scratch_172*scratch_214 - scratch_173*scratch_215 - scratch_174*scratch_216 - scratch_175*scratch_217 - scratch_176*scratch_218 - scratch_177*scratch_219 - scratch_178*scratch_220 - scratch_179*scratch_221);
    amrex::Real scratch_223 = scratch_208*(-A[271]*scratch_210 + A[285] - scratch_196*scratch_211 - scratch_197*scratch_212 - scratch_198*scratch_213 - scratch_199*scratch_214 - scratch_200*scratch_215 - scratch_201*scratch_216 - scratch_202*scratch_217 - scratch_203*scratch_218 - scratch_204*scratch_219 - scratch_205*scratch_220 - scratch_206*scratch_221 - scratch_207*scratch_222);
    amrex::Real scratch_224 = -scratch_120*scratch_220 - scratch_143*scratch_221 - scratch_15*scratch_213 - scratch_168*scratch_222 - scratch_195*scratch_223 + scratch_210 - scratch_211*scratch_4 - scratch_212*scratch_7 - scratch_214*scratch_24 - scratch_215*scratch_35 - scratch_216*scratch_48 - scratch_217*scratch_63 - scratch_218*scratch_80 - scratch_219*scratch_99;
    amrex::Real scratch_225 = -A[293]*scratch_4 + A[294];
    amrex::Real scratch_226 = -A[293]*scratch_3 + A[295] - scratch_225*scratch_6;
    amrex::Real scratch_227 = -A[293]*scratch_12 + A[296] - scratch_13*scratch_225 - scratch_14*scratch_226;
    amrex::Real scratch_228 = -A[293]*scratch_20 + A[297] - scratch_21*scratch_225 - scratch_22*scratch_226 - scratch_227*scratch_23;
    amrex::Real scratch_229 = -A[293]*scratch_30 + A[298] - scratch_225*scratch_31 - scratch_226*scratch_32 - scratch_227*scratch_33 - scratch_228*scratch_34;
    amrex::Real scratch_230 = -A[293]*scratch_42 + A[299] - scratch_225*scratch_43 - scratch_226*scratch_44 - scratch_227*scratch_45 - scratch_228*scratch_46 - scratch_229*scratch_47;
    amrex::Real scratch_231 = -A[293]*scratch_56 + A[300] - scratch_225*scratch_57 - scratch_226*scratch_58 - scratch_227*scratch_59 - scratch_228*scratch_60 - scratch_229*scratch_61 - scratch_230*scratch_62;
    amrex::Real scratch_232 = -A[293]*scratch_72 + A[301] - scratch_225*scratch_73 - scratch_226*scratch_74 - scratch_227*scratch_75 - scratch_228*scratch_76 - scratch_229*scratch_77 - scratch_230*scratch_78 - scratch_231*scratch_79;
    amrex::Real scratch_233 = -A[293]*scratch_90 + A[302] - scratch_225*scratch_91 - scratch_226*scratch_92 - scratch_227*scratch_93 - scratch_228*scratch_94 - scratch_229*scratch_95 - scratch_230*scratch_96 - scratch_231*scratch_97 - scratch_232*scratch_98;
    amrex::Real scratch_234 = -A[293]*scratch_110 + A[303] - scratch_111*scratch_225 - scratch_112*scratch_226 - scratch_113*scratch_227 - scratch_114*scratch_228 - scratch_115*scratch_229 - scratch_116*scratch_230 - scratch_117*scratch_231 - scratch_118*scratch_232 - scratch_119*scratch_233;
    amrex::Real scratch_235 = -A[293]*scratch_132 + A[304] - scratch_133*scratch_225 - scratch_134*scratch_226 - scratch_135*scratch_227 - scratch_136*scratch_228 - scratch_137*scratch_229 - scratch_138*scratch_230 - scratch_139*scratch_231 - scratch_140*scratch_232 - scratch_141*scratch_233 - scratch_142*scratch_234;
    amrex::Real scratch_236 = -A[293]*scratch_156 + A[305] - scratch_157*scratch_225 - scratch_158*scratch_226 - scratch_159*scratch_227 - scratch_160*scratch_228 - scratch_161*scratch_229 - scratch_162*scratch_230 - scratch_163*scratch_231 - scratch_164*scratch_232 - scratch_165*scratch_233 - scratch_166*scratch_234 - scratch_167*scratch_235;
    amrex::Real scratch_237 = -A[293]*scratch_182 + A[306] - scratch_183*scratch_225 - scratch_184*scratch_226 - scratch_185*scratch_227 - scratch_186*scratch_228 - scratch_187*scratch_229 - scratch_188*scratch_230 - scratch_189*scratch_231 - scratch_190*scratch_232 - scratch_191*scratch_233 - scratch_192*scratch_234 - scratch_193*scratch_235 - scratch_194*scratch_236;
    amrex::Real scratch_238 = 1.0/(-A[293]*scratch_210 + A[307] - scratch_211*scratch_225 - scratch_212*scratch_226 - scratch_213*scratch_227 - scratch_214*scratch_228 - scratch_215*scratch_229 - scratch_216*scratch_230 - scratch_217*scratch_231 - scratch_218*scratch_232 - scratch_219*scratch_233 - scratch_220*scratch_234 - scratch_221*scratch_235 - scratch_222*scratch_236 - scratch_223*scratch_237);
    amrex::Real scratch_239 = scratch_238*(-A[293]*scratch_2 - scratch_10*scratch_225 - scratch_109*scratch_233 - scratch_11*scratch_226 - scratch_131*scratch_234 - scratch_155*scratch_235 - scratch_181*scratch_236 - scratch_19*scratch_227 - scratch_209*scratch_237 - scratch_228*scratch_29 - scratch_229*scratch_41 - scratch_230*scratch_55 - scratch_231*scratch_71 - scratch_232*scratch_89);
    amrex::Real scratch_240 = A[15]*scratch_0;
    amrex::Real scratch_241 = scratch_5*(-A[22]*scratch_240 + A[37]);
    amrex::Real scratch_242 = scratch_9*(-A[44]*scratch_240 + A[59] - scratch_241*scratch_8);
    amrex::Real scratch_243 = scratch_18*(-A[66]*scratch_240 + A[81] - scratch_16*scratch_241 - scratch_17*scratch_242);
    amrex::Real scratch_244 = scratch_28*(A[103] - A[88]*scratch_240 - scratch_241*scratch_25 - scratch_242*scratch_26 - scratch_243*scratch_27);
    amrex::Real scratch_245 = scratch_40*(-A[110]*scratch_240 + A[125] - scratch_241*scratch_36 - scratch_242*scratch_37 - scratch_243*scratch_38 - scratch_244*scratch_39);
    amrex::Real scratch_246 = scratch_54*(-A[132]*scratch_240 + A[147] - scratch_241*scratch_49 - scratch_242*scratch_50 - scratch_243*scratch_51 - scratch_244*scratch_52 - scratch_245*scratch_53);
    amrex::Real scratch_247 = scratch_70*(-A[154]*scratch_240 - scratch_241*scratch_64 - scratch_242*scratch_65 - scratch_243*scratch_66 - scratch_244*scratch_67 - scratch_245*scratch_68 - scratch_246*scratch_69);
    amrex::Real scratch_248 = scratch_88*(-A[173]*scratch_240 - scratch_241*scratch_81 - scratch_242*scratch_82 - scratch_243*scratch_83 - scratch_244*scratch_84 - scratch_245*scratch_85 - scratch_246*scratch_86 - scratch_247*scratch_87);
    amrex::Real scratch_249 = scratch_108*(-A[191]*scratch_240 + A[205] - scratch_100*scratch_241 - scratch_101*scratch_242 - scratch_102*scratch_243 - scratch_103*scratch_244 - scratch_104*scratch_245 - scratch_105*scratch_246 - scratch_106*scratch_247 - scratch_107*scratch_248);
    amrex::Real scratch_250 = scratch_130*(-A[210]*scratch_240 - scratch_121*scratch_241 - scratch_122*scratch_242 - scratch_123*scratch_243 - scratch_124*scratch_244 - scratch_125*scratch_245 - scratch_126*scratch_246 - scratch_127*scratch_247 - scratch_128*scratch_248 - scratch_129*scratch_249);
    amrex::Real scratch_251 = scratch_154*(-A[227]*scratch_240 + A[242] - scratch_144*scratch_241 - scratch_145*scratch_242 - scratch_146*scratch_243 - scratch_147*scratch_244 - scratch_148*scratch_245 - scratch_149*scratch_246 - scratch_150*scratch_247 - scratch_151*scratch_248 - scratch_152*scratch_249 - scratch_153*scratch_250);
    amrex::Real scratch_252 = scratch_180*(-A[249]*scratch_240 + A[264] - scratch_169*scratch_241 - scratch_170*scratch_242 - scratch_171*scratch_243 - scratch_172*scratch_244 - scratch_173*scratch_245 - scratch_174*scratch_246 - scratch_175*scratch_247 - scratch_176*scratch_248 - scratch_177*scratch_249 - scratch_178*scratch_250 - scratch_179*scratch_251);
    amrex::Real scratch_253 = scratch_208*(-A[271]*scratch_240 + A[286] - scratch_196*scratch_241 - scratch_197*scratch_242 - scratch_198*scratch_243 - scratch_199*scratch_244 - scratch_200*scratch_245 - scratch_201*scratch_246 - scratch_202*scratch_247 - scratch_203*scratch_248 - scratch_204*scratch_249 - scratch_205*scratch_250 - scratch_206*scratch_251 - scratch_207*scratch_252);
    amrex::Real scratch_254 = scratch_238*(-A[293]*scratch_240 + A[308] - scratch_225*scratch_241 - scratch_226*scratch_242 - scratch_227*scratch_243 - scratch_228*scratch_244 - scratch_229*scratch_245 - scratch_230*scratch_246 - scratch_231*scratch_247 - scratch_232*scratch_248 - scratch_233*scratch_249 - scratch_234*scratch_250 - scratch_235*scratch_251 - scratch_236*scratch_252 - scratch_237*scratch_253);
    amrex::Real scratch_255 = -scratch_120*scratch_250 - scratch_143*scratch_251 - scratch_15*scratch_243 - scratch_168*scratch_252 - scratch_195*scratch_253 - scratch_224*scratch_254 - scratch_24*scratch_244 + scratch_240 - scratch_241*scratch_4 - scratch_242*scratch_7 - scratch_245*scratch_35 - scratch_246*scratch_48 - scratch_247*scratch_63 - scratch_248*scratch_80 - scratch_249*scratch_99;
    amrex::Real scratch_256 = -A[311]*scratch_6 + A[312];
    amrex::Real scratch_257 = -A[311]*scratch_13 + A[313] - scratch_14*scratch_256;
    amrex::Real scratch_258 = -A[311]*scratch_21 + A[314] - scratch_22*scratch_256 - scratch_23*scratch_257;
    amrex::Real scratch_259 = -A[311]*scratch_31 - scratch_256*scratch_32 - scratch_257*scratch_33 - scratch_258*scratch_34;
    amrex::Real scratch_260 = -A[311]*scratch_43 + A[315] - scratch_256*scratch_44 - scratch_257*scratch_45 - scratch_258*scratch_46 - scratch_259*scratch_47;
    amrex::Real scratch_261 = -A[311]*scratch_57 - scratch_256*scratch_58 - scratch_257*scratch_59 - scratch_258*scratch_60 - scratch_259*scratch_61 - scratch_260*scratch_62;
    amrex::Real scratch_262 = -A[311]*scratch_73 - scratch_256*scratch_74 - scratch_257*scratch_75 - scratch_258*scratch_76 - scratch_259*scratch_77 - scratch_260*scratch_78 - scratch_261*scratch_79;
    amrex::Real scratch_263 = -A[311]*scratch_91 + A[316] - scratch_256*scratch_92 - scratch_257*scratch_93 - scratch_258*scratch_94 - scratch_259*scratch_95 - scratch_260*scratch_96 - scratch_261*scratch_97 - scratch_262*scratch_98;
    amrex::Real scratch_264 = -A[311]*scratch_111 - scratch_112*scratch_256 - scratch_113*scratch_257 - scratch_114*scratch_258 - scratch_115*scratch_259 - scratch_116*scratch_260 - scratch_117*scratch_261 - scratch_118*scratch_262 - scratch_119*scratch_263;
    amrex::Real scratch_265 = -A[311]*scratch_133 - scratch_134*scratch_256 - scratch_135*scratch_257 - scratch_136*scratch_258 - scratch_137*scratch_259 - scratch_138*scratch_260 - scratch_139*scratch_261 - scratch_140*scratch_262 - scratch_141*scratch_263 - scratch_142*scratch_264;
    amrex::Real scratch_266 = -A[311]*scratch_157 - scratch_158*scratch_256 - scratch_159*scratch_257 - scratch_160*scratch_258 - scratch_161*scratch_259 - scratch_162*scratch_260 - scratch_163*scratch_261 - scratch_164*scratch_262 - scratch_165*scratch_263 - scratch_166*scratch_264 - scratch_167*scratch_265;
    amrex::Real scratch_267 = -A[311]*scratch_183 - scratch_184*scratch_256 - scratch_185*scratch_257 - scratch_186*scratch_258 - scratch_187*scratch_259 - scratch_188*scratch_260 - scratch_189*scratch_261 - scratch_190*scratch_262 - scratch_191*scratch_263 - scratch_192*scratch_264 - scratch_193*scratch_265 - scratch_194*scratch_266;
    amrex::Real scratch_268 = -A[311]*scratch_211 + A[317] - scratch_212*scratch_256 - scratch_213*scratch_257 - scratch_214*scratch_258 - scratch_215*scratch_259 - scratch_216*scratch_260 - scratch_217*scratch_261 - scratch_218*scratch_262 - scratch_219*scratch_263 - scratch_220*scratch_264 - scratch_221*scratch_265 - scratch_222*scratch_266 - scratch_223*scratch_267;
    amrex::Real scratch_269 = 1.0/(-A[311]*scratch_241 + A[318] - scratch_242*scratch_256 - scratch_243*scratch_257 - scratch_244*scratch_258 - scratch_245*scratch_259 - scratch_246*scratch_260 - scratch_247*scratch_261 - scratch_248*scratch_262 - scratch_249*scratch_263 - scratch_250*scratch_264 - scratch_251*scratch_265 - scratch_252*scratch_266 - scratch_253*scratch_267 - scratch_254*scratch_268);
    amrex::Real scratch_270 = scratch_269*(-A[311]*scratch_10 - scratch_109*scratch_263 - scratch_11*scratch_256 - scratch_131*scratch_264 - scratch_155*scratch_265 - scratch_181*scratch_266 - scratch_19*scratch_257 - scratch_209*scratch_267 - scratch_239*scratch_268 - scratch_258*scratch_29 - scratch_259*scratch_41 - scratch_260*scratch_55 - scratch_261*scratch_71 - scratch_262*scratch_89);
    amrex::Real scratch_271 = A[16]*scratch_0;
    amrex::Real scratch_272 = scratch_5*(-A[22]*scratch_271 + A[38]);
    amrex::Real scratch_273 = scratch_9*(-A[44]*scratch_271 + A[60] - scratch_272*scratch_8);
    amrex::Real scratch_274 = scratch_18*(-A[66]*scratch_271 + A[82] - scratch_16*scratch_272 - scratch_17*scratch_273);
    amrex::Real scratch_275 = scratch_28*(A[104] - A[88]*scratch_271 - scratch_25*scratch_272 - scratch_26*scratch_273 - scratch_27*scratch_274);
    amrex::Real scratch_276 = scratch_40*(-A[110]*scratch_271 + A[126] - scratch_272*scratch_36 - scratch_273*scratch_37 - scratch_274*scratch_38 - scratch_275*scratch_39);
    amrex::Real scratch_277 = scratch_54*(-A[132]*scratch_271 + A[148] - scratch_272*scratch_49 - scratch_273*scratch_50 - scratch_274*scratch_51 - scratch_275*scratch_52 - scratch_276*scratch_53);
    amrex::Real scratch_278 = scratch_70*(-A[154]*scratch_271 + A[169] - scratch_272*scratch_64 - scratch_273*scratch_65 - scratch_274*scratch_66 - scratch_275*scratch_67 - scratch_276*scratch_68 - scratch_277*scratch_69);
    amrex::Real scratch_279 = scratch_88*(-A[173]*scratch_271 + A[187] - scratch_272*scratch_81 - scratch_273*scratch_82 - scratch_274*scratch_83 - scratch_275*scratch_84 - scratch_276*scratch_85 - scratch_277*scratch_86 - scratch_278*scratch_87);
    amrex::Real scratch_280 = scratch_108*(-A[191]*scratch_271 + A[206] - scratch_100*scratch_272 - scratch_101*scratch_273 - scratch_102*scratch_274 - scratch_103*scratch_275 - scratch_104*scratch_276 - scratch_105*scratch_277 - scratch_106*scratch_278 - scratch_107*scratch_279);
    amrex::Real scratch_281 = scratch_130*(-A[210]*scratch_271 - scratch_121*scratch_272 - scratch_122*scratch_273 - scratch_123*scratch_274 - scratch_124*scratch_275 - scratch_125*scratch_276 - scratch_126*scratch_277 - scratch_127*scratch_278 - scratch_128*scratch_279 - scratch_129*scratch_280);
    amrex::Real scratch_282 = scratch_154*(-A[227]*scratch_271 + A[243] - scratch_144*scratch_272 - scratch_145*scratch_273 - scratch_146*scratch_274 - scratch_147*scratch_275 - scratch_148*scratch_276 - scratch_149*scratch_277 - scratch_150*scratch_278 - scratch_151*scratch_279 - scratch_152*scratch_280 - scratch_153*scratch_281);
    amrex::Real scratch_283 = scratch_180*(-A[249]*scratch_271 + A[265] - scratch_169*scratch_272 - scratch_170*scratch_273 - scratch_171*scratch_274 - scratch_172*scratch_275 - scratch_173*scratch_276 - scratch_174*scratch_277 - scratch_175*scratch_278 - scratch_176*scratch_279 - scratch_177*scratch_280 - scratch_178*scratch_281 - scratch_179*scratch_282);
    amrex::Real scratch_284 = scratch_208*(-A[271]*scratch_271 + A[287] - scratch_196*scratch_272 - scratch_197*scratch_273 - scratch_198*scratch_274 - scratch_199*scratch_275 - scratch_200*scratch_276 - scratch_201*scratch_277 - scratch_202*scratch_278 - scratch_203*scratch_279 - scratch_204*scratch_280 - scratch_205*scratch_281 - scratch_206*scratch_282 - scratch_207*scratch_283);
    amrex::Real scratch_285 = scratch_238*(-A[293]*scratch_271 - scratch_225*scratch_272 - scratch_226*scratch_273 - scratch_227*scratch_274 - scratch_228*scratch_275 - scratch_229*scratch_276 - scratch_230*scratch_277 - scratch_231*scratch_278 - scratch_232*scratch_279 - scratch_233*scratch_280 - scratch_234*scratch_281 - scratch_235*scratch_282 - scratch_236*scratch_283 - scratch_237*scratch_284);
    amrex::Real scratch_286 = scratch_269*(-A[311]*scratch_272 - scratch_256*scratch_273 - scratch_257*scratch_274 - scratch_258*scratch_275 - scratch_259*scratch_276 - scratch_260*scratch_277 - scratch_261*scratch_278 - scratch_262*scratch_279 - scratch_263*scratch_280 - scratch_264*scratch_281 - scratch_265*scratch_282 - scratch_266*scratch_283 - scratch_267*scratch_284 - scratch_268*scratch_285);
    amrex::Real scratch_287 = -scratch_120*scratch_281 - scratch_143*scratch_282 - scratch_15*scratch_274 - scratch_168*scratch_283 - scratch_195*scratch_284 - scratch_224*scratch_285 - scratch_24*scratch_275 - scratch_255*scratch_286 + scratch_271 - scratch_272*scratch_4 - scratch_273*scratch_7 - scratch_276*scratch_35 - scratch_277*scratch_48 - scratch_278*scratch_63 - scratch_279*scratch_80 - scratch_280*scratch_99;
    amrex::Real scratch_288 = -A[320]*scratch_6 + A[321];
    amrex::Real scratch_289 = -A[320]*scratch_13 + A[322] - scratch_14*scratch_288;
    amrex::Real scratch_290 = -A[320]*scratch_21 - scratch_22*scratch_288 - scratch_23*scratch_289;
    amrex::Real scratch_291 = -A[320]*scratch_31 - scratch_288*scratch_32 - scratch_289*scratch_33 - scratch_290*scratch_34;
    amrex::Real scratch_292 = -A[320]*scratch_43 + A[323] - scratch_288*scratch_44 - scratch_289*scratch_45 - scratch_290*scratch_46 - scratch_291*scratch_47;
    amrex::Real scratch_293 = -A[320]*scratch_57 + A[324] - scratch_288*scratch_58 - scratch_289*scratch_59 - scratch_290*scratch_60 - scratch_291*scratch_61 - scratch_292*scratch_62;
    amrex::Real scratch_294 = -A[320]*scratch_73 + A[325] - scratch_288*scratch_74 - scratch_289*scratch_75 - scratch_290*scratch_76 - scratch_291*scratch_77 - scratch_292*scratch_78 - scratch_293*scratch_79;
    amrex::Real scratch_295 = -A[320]*scratch_91 + A[326] - scratch_288*scratch_92 - scratch_289*scratch_93 - scratch_290*scratch_94 - scratch_291*scratch_95 - scratch_292*scratch_96 - scratch_293*scratch_97 - scratch_294*scratch_98;
    amrex::Real scratch_296 = -A[320]*scratch_111 - scratch_112*scratch_288 - scratch_113*scratch_289 - scratch_114*scratch_290 - scratch_115*scratch_291 - scratch_116*scratch_292 - scratch_117*scratch_293 - scratch_118*scratch_294 - scratch_119*scratch_295;
    amrex::Real scratch_297 = -A[320]*scratch_133 - scratch_134*scratch_288 - scratch_135*scratch_289 - scratch_136*scratch_290 - scratch_137*scratch_291 - scratch_138*scratch_292 - scratch_139*scratch_293 - scratch_140*scratch_294 - scratch_141*scratch_295 - scratch_142*scratch_296;
    amrex::Real scratch_298 = -A[320]*scratch_157 - scratch_158*scratch_288 - scratch_159*scratch_289 - scratch_160*scratch_290 - scratch_161*scratch_291 - scratch_162*scratch_292 - scratch_163*scratch_293 - scratch_164*scratch_294 - scratch_165*scratch_295 - scratch_166*scratch_296 - scratch_167*scratch_297;
    amrex::Real scratch_299 = -A[320]*scratch_183 + A[327] - scratch_184*scratch_288 - scratch_185*scratch_289 - scratch_186*scratch_290 - scratch_187*scratch_291 - scratch_188*scratch_292 - scratch_189*scratch_293 - scratch_190*scratch_294 - scratch_191*scratch_295 - scratch_192*scratch_296 - scratch_193*scratch_297 - scratch_194*scratch_298;
    amrex::Real scratch_300 = -A[320]*scratch_211 - scratch_212*scratch_288 - scratch_213*scratch_289 - scratch_214*scratch_290 - scratch_215*scratch_291 - scratch_216*scratch_292 - scratch_217*scratch_293 - scratch_218*scratch_294 - scratch_219*scratch_295 - scratch_220*scratch_296 - scratch_221*scratch_297 - scratch_222*scratch_298 - scratch_223*scratch_299;
    amrex::Real scratch_301 = -A[320]*scratch_241 - scratch_242*scratch_288 - scratch_243*scratch_289 - scratch_244*scratch_290 - scratch_245*scratch_291 - scratch_246*scratch_292 - scratch_247*scratch_293 - scratch_248*scratch_294 - scratch_249*scratch_295 - scratch_250*scratch_296 - scratch_251*scratch_297 - scratch_252*scratch_298 - scratch_253*scratch_299 - scratch_254*scratch_300;
    amrex::Real scratch_302 = 1.0/(-A[320]*scratch_272 + A[328] - scratch_273*scratch_288 - scratch_274*scratch_289 - scratch_275*scratch_290 - scratch_276*scratch_291 - scratch_277*scratch_292 - scratch_278*scratch_293 - scratch_279*scratch_294 - scratch_280*scratch_295 - scratch_281*scratch_296 - scratch_282*scratch_297 - scratch_283*scratch_298 - scratch_284*scratch_299 - scratch_285*scratch_300 - scratch_286*scratch_301);
    amrex::Real scratch_303 = scratch_302*(-A[320]*scratch_10 - scratch_109*scratch_295 - scratch_11*scratch_288 - scratch_131*scratch_296 - scratch_155*scratch_297 - scratch_181*scratch_298 - scratch_19*scratch_289 - scratch_209*scratch_299 - scratch_239*scratch_300 - scratch_270*scratch_301 - scratch_29*scratch_290 - scratch_291*scratch_41 - scratch_292*scratch_55 - scratch_293*scratch_71 - scratch_294*scratch_89);
    amrex::Real scratch_304 = A[17]*scratch_0;
    amrex::Real scratch_305 = scratch_5*(-A[22]*scratch_304 + A[39]);
    amrex::Real scratch_306 = scratch_9*(-A[44]*scratch_304 + A[61] - scratch_305*scratch_8);
    amrex::Real scratch_307 = scratch_18*(-A[66]*scratch_304 + A[83] - scratch_16*scratch_305 - scratch_17*scratch_306);
    amrex::Real scratch_308 = scratch_28*(A[105] - A[88]*scratch_304 - scratch_25*scratch_305 - scratch_26*scratch_306 - scratch_27*scratch_307);
    amrex::Real scratch_309 = scratch_40*(-A[110]*scratch_304 + A[127] - scratch_305*scratch_36 - scratch_306*scratch_37 - scratch_307*scratch_38 - scratch_308*scratch_39);
    amrex::Real scratch_310 = scratch_54*(-A[132]*scratch_304 + A[149] - scratch_305*scratch_49 - scratch_306*scratch_50 - scratch_307*scratch_51 - scratch_308*scratch_52 - scratch_309*scratch_53);
    amrex::Real scratch_311 = scratch_70*(-A[154]*scratch_304 - scratch_305*scratch_64 - scratch_306*scratch_65 - scratch_307*scratch_66 - scratch_308*scratch_67 - scratch_309*scratch_68 - scratch_310*scratch_69);
    amrex::Real scratch_312 = scratch_88*(-A[173]*scratch_304 - scratch_305*scratch_81 - scratch_306*scratch_82 - scratch_307*scratch_83 - scratch_308*scratch_84 - scratch_309*scratch_85 - scratch_310*scratch_86 - scratch_311*scratch_87);
    amrex::Real scratch_313 = scratch_108*(-A[191]*scratch_304 + A[207] - scratch_100*scratch_305 - scratch_101*scratch_306 - scratch_102*scratch_307 - scratch_103*scratch_308 - scratch_104*scratch_309 - scratch_105*scratch_310 - scratch_106*scratch_311 - scratch_107*scratch_312);
    amrex::Real scratch_314 = scratch_130*(-A[210]*scratch_304 + A[224] - scratch_121*scratch_305 - scratch_122*scratch_306 - scratch_123*scratch_307 - scratch_124*scratch_308 - scratch_125*scratch_309 - scratch_126*scratch_310 - scratch_127*scratch_311 - scratch_128*scratch_312 - scratch_129*scratch_313);
    amrex::Real scratch_315 = scratch_154*(-A[227]*scratch_304 + A[244] - scratch_144*scratch_305 - scratch_145*scratch_306 - scratch_146*scratch_307 - scratch_147*scratch_308 - scratch_148*scratch_309 - scratch_149*scratch_310 - scratch_150*scratch_311 - scratch_151*scratch_312 - scratch_152*scratch_313 - scratch_153*scratch_314);
    amrex::Real scratch_316 = scratch_180*(-A[249]*scratch_304 + A[266] - scratch_169*scratch_305 - scratch_170*scratch_306 - scratch_171*scratch_307 - scratch_172*scratch_308 - scratch_173*scratch_309 - scratch_174*scratch_310 - scratch_175*scratch_311 - scratch_176*scratch_312 - scratch_177*scratch_313 - scratch_178*scratch_314 - scratch_179*scratch_315);
    amrex::Real scratch_317 = scratch_208*(-A[271]*scratch_304 + A[288] - scratch_196*scratch_305 - scratch_197*scratch_306 - scratch_198*scratch_307 - scratch_199*scratch_308 - scratch_200*scratch_309 - scratch_201*scratch_310 - scratch_202*scratch_311 - scratch_203*scratch_312 - scratch_204*scratch_313 - scratch_205*scratch_314 - scratch_206*scratch_315 - scratch_207*scratch_316);
    amrex::Real scratch_318 = scratch_238*(-A[293]*scratch_304 + A[309] - scratch_225*scratch_305 - scratch_226*scratch_306 - scratch_227*scratch_307 - scratch_228*scratch_308 - scratch_229*scratch_309 - scratch_230*scratch_310 - scratch_231*scratch_311 - scratch_232*scratch_312 - scratch_233*scratch_313 - scratch_234*scratch_314 - scratch_235*scratch_315 - scratch_236*scratch_316 - scratch_237*scratch_317);
    amrex::Real scratch_319 = scratch_269*(-A[311]*scratch_305 - scratch_256*scratch_306 - scratch_257*scratch_307 - scratch_258*scratch_308 - scratch_259*scratch_309 - scratch_260*scratch_310 - scratch_261*scratch_311 - scratch_262*scratch_312 - scratch_263*scratch_313 - scratch_264*scratch_314 - scratch_265*scratch_315 - scratch_266*scratch_316 - scratch_267*scratch_317 - scratch_268*scratch_318);
    amrex::Real scratch_320 = scratch_302*(-A[320]*scratch_305 + A[329] - scratch_288*scratch_306 - scratch_289*scratch_307 - scratch_290*scratch_308 - scratch_291*scratch_309 - scratch_292*scratch_310 - scratch_293*scratch_311 - scratch_294*scratch_312 - scratch_295*scratch_313 - scratch_296*scratch_314 - scratch_297*scratch_315 - scratch_298*scratch_316 - scratch_299*scratch_317 - scratch_300*scratch_318 - scratch_301*scratch_319);
    amrex::Real scratch_321 = -scratch_120*scratch_314 - scratch_143*scratch_315 - scratch_15*scratch_307 - scratch_168*scratch_316 - scratch_195*scratch_317 - scratch_224*scratch_318 - scratch_24*scratch_308 - scratch_255*scratch_319 - scratch_287*scratch_320 + scratch_304 - scratch_305*scratch_4 - scratch_306*scratch_7 - scratch_309*scratch_35 - scratch_310*scratch_48 - scratch_311*scratch_63 - scratch_312*scratch_80 - scratch_313*scratch_99;
    amrex::Real scratch_322 = -A[331]*scratch_4 + A[332];
    amrex::Real scratch_323 = -A[331]*scratch_3 + A[333] - scratch_322*scratch_6;
    amrex::Real scratch_324 = -A[331]*scratch_12 + A[334] - scratch_13*scratch_322 - scratch_14*scratch_323;
    amrex::Real scratch_325 = -A[331]*scratch_20 + A[335] - scratch_21*scratch_322 - scratch_22*scratch_323 - scratch_23*scratch_324;
    amrex::Real scratch_326 = -A[331]*scratch_30 + A[336] - scratch_31*scratch_322 - scratch_32*scratch_323 - scratch_324*scratch_33 - scratch_325*scratch_34;
    amrex::Real scratch_327 = -A[331]*scratch_42 + A[337] - scratch_322*scratch_43 - scratch_323*scratch_44 - scratch_324*scratch_45 - scratch_325*scratch_46 - scratch_326*scratch_47;
    amrex::Real scratch_328 = -A[331]*scratch_56 - scratch_322*scratch_57 - scratch_323*scratch_58 - scratch_324*scratch_59 - scratch_325*scratch_60 - scratch_326*scratch_61 - scratch_327*scratch_62;
    amrex::Real scratch_329 = -A[331]*scratch_72 - scratch_322*scratch_73 - scratch_323*scratch_74 - scratch_324*scratch_75 - scratch_325*scratch_76 - scratch_326*scratch_77 - scratch_327*scratch_78 - scratch_328*scratch_79;
    amrex::Real scratch_330 = -A[331]*scratch_90 + A[338] - scratch_322*scratch_91 - scratch_323*scratch_92 - scratch_324*scratch_93 - scratch_325*scratch_94 - scratch_326*scratch_95 - scratch_327*scratch_96 - scratch_328*scratch_97 - scratch_329*scratch_98;
    amrex::Real scratch_331 = -A[331]*scratch_110 + A[339] - scratch_111*scratch_322 - scratch_112*scratch_323 - scratch_113*scratch_324 - scratch_114*scratch_325 - scratch_115*scratch_326 - scratch_116*scratch_327 - scratch_117*scratch_328 - scratch_118*scratch_329 - scratch_119*scratch_330;
    amrex::Real scratch_332 = -A[331]*scratch_132 - scratch_133*scratch_322 - scratch_134*scratch_323 - scratch_135*scratch_324 - scratch_136*scratch_325 - scratch_137*scratch_326 - scratch_138*scratch_327 - scratch_139*scratch_328 - scratch_140*scratch_329 - scratch_141*scratch_330 - scratch_142*scratch_331;
    amrex::Real scratch_333 = -A[331]*scratch_156 - scratch_157*scratch_322 - scratch_158*scratch_323 - scratch_159*scratch_324 - scratch_160*scratch_325 - scratch_161*scratch_326 - scratch_162*scratch_327 - scratch_163*scratch_328 - scratch_164*scratch_329 - scratch_165*scratch_330 - scratch_166*scratch_331 - scratch_167*scratch_332;
    amrex::Real scratch_334 = -A[331]*scratch_182 - scratch_183*scratch_322 - scratch_184*scratch_323 - scratch_185*scratch_324 - scratch_186*scratch_325 - scratch_187*scratch_326 - scratch_188*scratch_327 - scratch_189*scratch_328 - scratch_190*scratch_329 - scratch_191*scratch_330 - scratch_192*scratch_331 - scratch_193*scratch_332 - scratch_194*scratch_333;
    amrex::Real scratch_335 = -A[331]*scratch_210 + A[340] - scratch_211*scratch_322 - scratch_212*scratch_323 - scratch_213*scratch_324 - scratch_214*scratch_325 - scratch_215*scratch_326 - scratch_216*scratch_327 - scratch_217*scratch_328 - scratch_218*scratch_329 - scratch_219*scratch_330 - scratch_220*scratch_331 - scratch_221*scratch_332 - scratch_222*scratch_333 - scratch_223*scratch_334;
    amrex::Real scratch_336 = -A[331]*scratch_240 - scratch_241*scratch_322 - scratch_242*scratch_323 - scratch_243*scratch_324 - scratch_244*scratch_325 - scratch_245*scratch_326 - scratch_246*scratch_327 - scratch_247*scratch_328 - scratch_248*scratch_329 - scratch_249*scratch_330 - scratch_250*scratch_331 - scratch_251*scratch_332 - scratch_252*scratch_333 - scratch_253*scratch_334 - scratch_254*scratch_335;
    amrex::Real scratch_337 = -A[331]*scratch_271 + A[341] - scratch_272*scratch_322 - scratch_273*scratch_323 - scratch_274*scratch_324 - scratch_275*scratch_325 - scratch_276*scratch_326 - scratch_277*scratch_327 - scratch_278*scratch_328 - scratch_279*scratch_329 - scratch_280*scratch_330 - scratch_281*scratch_331 - scratch_282*scratch_332 - scratch_283*scratch_333 - scratch_284*scratch_334 - scratch_285*scratch_335 - scratch_286*scratch_336;
    amrex::Real scratch_338 = 1.0/(-A[331]*scratch_304 + A[342] - scratch_305*scratch_322 - scratch_306*scratch_323 - scratch_307*scratch_324 - scratch_308*scratch_325 - scratch_309*scratch_326 - scratch_310*scratch_327 - scratch_311*scratch_328 - scratch_312*scratch_329 - scratch_313*scratch_330 - scratch_314*scratch_331 - scratch_315*scratch_332 - scratch_316*scratch_333 - scratch_317*scratch_334 - scratch_318*scratch_335 - scratch_319*scratch_336 - scratch_320*scratch_337);
    amrex::Real scratch_339 = scratch_338*(-A[331]*scratch_2 - scratch_10*scratch_322 - scratch_109*scratch_330 - scratch_11*scratch_323 - scratch_131*scratch_331 - scratch_155*scratch_332 - scratch_181*scratch_333 - scratch_19*scratch_324 - scratch_209*scratch_334 - scratch_239*scratch_335 - scratch_270*scratch_336 - scratch_29*scratch_325 - scratch_303*scratch_337 - scratch_326*scratch_41 - scratch_327*scratch_55 - scratch_328*scratch_71 - scratch_329*scratch_89);
    amrex::Real scratch_340 = -A[345]*scratch_4 + A[346];
    amrex::Real scratch_341 = -A[345]*scratch_3 + A[347] - scratch_340*scratch_6;
    amrex::Real scratch_342 = -A[345]*scratch_12 - scratch_13*scratch_340 - scratch_14*scratch_341;
    amrex::Real scratch_343 = -A[345]*scratch_20 + A[348] - scratch_21*scratch_340 - scratch_22*scratch_341 - scratch_23*scratch_342;
    amrex::Real scratch_344 = -A[345]*scratch_30 + A[349] - scratch_31*scratch_340 - scratch_32*scratch_341 - scratch_33*scratch_342 - scratch_34*scratch_343;
    amrex::Real scratch_345 = -A[345]*scratch_42 - scratch_340*scratch_43 - scratch_341*scratch_44 - scratch_342*scratch_45 - scratch_343*scratch_46 - scratch_344*scratch_47;
    amrex::Real scratch_346 = -A[345]*scratch_56 - scratch_340*scratch_57 - scratch_341*scratch_58 - scratch_342*scratch_59 - scratch_343*scratch_60 - scratch_344*scratch_61 - scratch_345*scratch_62;
    amrex::Real scratch_347 = -A[345]*scratch_72 - scratch_340*scratch_73 - scratch_341*scratch_74 - scratch_342*scratch_75 - scratch_343*scratch_76 - scratch_344*scratch_77 - scratch_345*scratch_78 - scratch_346*scratch_79;
    amrex::Real scratch_348 = -A[345]*scratch_90 + A[350] - scratch_340*scratch_91 - scratch_341*scratch_92 - scratch_342*scratch_93 - scratch_343*scratch_94 - scratch_344*scratch_95 - scratch_345*scratch_96 - scratch_346*scratch_97 - scratch_347*scratch_98;
    amrex::Real scratch_349 = -A[345]*scratch_110 + A[351] - scratch_111*scratch_340 - scratch_112*scratch_341 - scratch_113*scratch_342 - scratch_114*scratch_343 - scratch_115*scratch_344 - scratch_116*scratch_345 - scratch_117*scratch_346 - scratch_118*scratch_347 - scratch_119*scratch_348;
    amrex::Real scratch_350 = -A[345]*scratch_132 - scratch_133*scratch_340 - scratch_134*scratch_341 - scratch_135*scratch_342 - scratch_136*scratch_343 - scratch_137*scratch_344 - scratch_138*scratch_345 - scratch_139*scratch_346 - scratch_140*scratch_347 - scratch_141*scratch_348 - scratch_142*scratch_349;
    amrex::Real scratch_351 = -A[345]*scratch_156 - scratch_157*scratch_340 - scratch_158*scratch_341 - scratch_159*scratch_342 - scratch_160*scratch_343 - scratch_161*scratch_344 - scratch_162*scratch_345 - scratch_163*scratch_346 - scratch_164*scratch_347 - scratch_165*scratch_348 - scratch_166*scratch_349 - scratch_167*scratch_350;
    amrex::Real scratch_352 = -A[345]*scratch_182 - scratch_183*scratch_340 - scratch_184*scratch_341 - scratch_185*scratch_342 - scratch_186*scratch_343 - scratch_187*scratch_344 - scratch_188*scratch_345 - scratch_189*scratch_346 - scratch_190*scratch_347 - scratch_191*scratch_348 - scratch_192*scratch_349 - scratch_193*scratch_350 - scratch_194*scratch_351;
    amrex::Real scratch_353 = -A[345]*scratch_210 - scratch_211*scratch_340 - scratch_212*scratch_341 - scratch_213*scratch_342 - scratch_214*scratch_343 - scratch_215*scratch_344 - scratch_216*scratch_345 - scratch_217*scratch_346 - scratch_218*scratch_347 - scratch_219*scratch_348 - scratch_220*scratch_349 - scratch_221*scratch_350 - scratch_222*scratch_351 - scratch_223*scratch_352;
    amrex::Real scratch_354 = -A[345]*scratch_240 - scratch_241*scratch_340 - scratch_242*scratch_341 - scratch_243*scratch_342 - scratch_244*scratch_343 - scratch_245*scratch_344 - scratch_246*scratch_345 - scratch_247*scratch_346 - scratch_248*scratch_347 - scratch_249*scratch_348 - scratch_250*scratch_349 - scratch_251*scratch_350 - scratch_252*scratch_351 - scratch_253*scratch_352 - scratch_254*scratch_353;
    amrex::Real scratch_355 = -A[345]*scratch_271 - scratch_272*scratch_340 - scratch_273*scratch_341 - scratch_274*scratch_342 - scratch_275*scratch_343 - scratch_276*scratch_344 - scratch_277*scratch_345 - scratch_278*scratch_346 - scratch_279*scratch_347 - scratch_280*scratch_348 - scratch_281*scratch_349 - scratch_282*scratch_350 - scratch_283*scratch_351 - scratch_284*scratch_352 - scratch_285*scratch_353 - scratch_286*scratch_354;
    amrex::Real scratch_356 = -A[345]*scratch_304 + A[352] - scratch_305*scratch_340 - scratch_306*scratch_341 - scratch_307*scratch_342 - scratch_308*scratch_343 - scratch_309*scratch_344 - scratch_310*scratch_345 - scratch_311*scratch_346 - scratch_312*scratch_347 - scratch_313*scratch_348 - scratch_314*scratch_349 - scratch_315*scratch_350 - scratch_316*scratch_351 - scratch_317*scratch_352 - scratch_318*scratch_353 - scratch_319*scratch_354 - scratch_320*scratch_355;
    amrex::Real scratch_357 = -A[345]*scratch_2 - scratch_10*scratch_340 - scratch_109*scratch_348 - scratch_11*scratch_341 - scratch_131*scratch_349 - scratch_155*scratch_350 - scratch_181*scratch_351 - scratch_19*scratch_342 - scratch_209*scratch_352 - scratch_239*scratch_353 - scratch_270*scratch_354 - scratch_29*scratch_343 - scratch_303*scratch_355 - scratch_339*scratch_356 - scratch_344*scratch_41 - scratch_345*scratch_55 - scratch_346*scratch_71 - scratch_347*scratch_89;
    amrex::Real scratch_358 = A[18]*scratch_0;
    amrex::Real scratch_359 = scratch_5*(-A[22]*scratch_358 + A[40]);
    amrex::Real scratch_360 = scratch_9*(-A[44]*scratch_358 + A[62] - scratch_359*scratch_8);
    amrex::Real scratch_361 = scratch_18*(-A[66]*scratch_358 + A[84] - scratch_16*scratch_359 - scratch_17*scratch_360);
    amrex::Real scratch_362 = scratch_28*(A[106] - A[88]*scratch_358 - scratch_25*scratch_359 - scratch_26*scratch_360 - scratch_27*scratch_361);
    amrex::Real scratch_363 = scratch_40*(-A[110]*scratch_358 + A[128] - scratch_359*scratch_36 - scratch_360*scratch_37 - scratch_361*scratch_38 - scratch_362*scratch_39);
    amrex::Real scratch_364 = scratch_54*(-A[132]*scratch_358 + A[150] - scratch_359*scratch_49 - scratch_360*scratch_50 - scratch_361*scratch_51 - scratch_362*scratch_52 - scratch_363*scratch_53);
    amrex::Real scratch_365 = scratch_70*(-A[154]*scratch_358 - scratch_359*scratch_64 - scratch_360*scratch_65 - scratch_361*scratch_66 - scratch_362*scratch_67 - scratch_363*scratch_68 - scratch_364*scratch_69);
    amrex::Real scratch_366 = scratch_88*(-A[173]*scratch_358 - scratch_359*scratch_81 - scratch_360*scratch_82 - scratch_361*scratch_83 - scratch_362*scratch_84 - scratch_363*scratch_85 - scratch_364*scratch_86 - scratch_365*scratch_87);
    amrex::Real scratch_367 = scratch_108*(-A[191]*scratch_358 + A[208] - scratch_100*scratch_359 - scratch_101*scratch_360 - scratch_102*scratch_361 - scratch_103*scratch_362 - scratch_104*scratch_363 - scratch_105*scratch_364 - scratch_106*scratch_365 - scratch_107*scratch_366);
    amrex::Real scratch_368 = scratch_130*(-A[210]*scratch_358 + A[225] - scratch_121*scratch_359 - scratch_122*scratch_360 - scratch_123*scratch_361 - scratch_124*scratch_362 - scratch_125*scratch_363 - scratch_126*scratch_364 - scratch_127*scratch_365 - scratch_128*scratch_366 - scratch_129*scratch_367);
    amrex::Real scratch_369 = scratch_154*(-A[227]*scratch_358 + A[245] - scratch_144*scratch_359 - scratch_145*scratch_360 - scratch_146*scratch_361 - scratch_147*scratch_362 - scratch_148*scratch_363 - scratch_149*scratch_364 - scratch_150*scratch_365 - scratch_151*scratch_366 - scratch_152*scratch_367 - scratch_153*scratch_368);
    amrex::Real scratch_370 = scratch_180*(-A[249]*scratch_358 + A[267] - scratch_169*scratch_359 - scratch_170*scratch_360 - scratch_171*scratch_361 - scratch_172*scratch_362 - scratch_173*scratch_363 - scratch_174*scratch_364 - scratch_175*scratch_365 - scratch_176*scratch_366 - scratch_177*scratch_367 - scratch_178*scratch_368 - scratch_179*scratch_369);
    amrex::Real scratch_371 = scratch_208*(-A[271]*scratch_358 + A[289] - scratch_196*scratch_359 - scratch_197*scratch_360 - scratch_198*scratch_361 - scratch_199*scratch_362 - scratch_200*scratch_363 - scratch_201*scratch_364 - scratch_202*scratch_365 - scratch_203*scratch_366 - scratch_204*scratch_367 - scratch_205*scratch_368 - scratch_206*scratch_369 - scratch_207*scratch_370);
    amrex::Real scratch_372 = scratch_238*(-A[293]*scratch_358 - scratch_225*scratch_359 - scratch_226*scratch_360 - scratch_227*scratch_361 - scratch_228*scratch_362 - scratch_229*scratch_363 - scratch_230*scratch_364 - scratch_231*scratch_365 - scratch_232*scratch_366 - scratch_233*scratch_367 - scratch_234*scratch_368 - scratch_235*scratch_369 - scratch_236*scratch_370 - scratch_237*scratch_371);
    amrex::Real scratch_373 = scratch_269*(-A[311]*scratch_359 - scratch_256*scratch_360 - scratch_257*scratch_361 - scratch_258*scratch_362 - scratch_259*scratch_363 - scratch_260*scratch_364 - scratch_261*scratch_365 - scratch_262*scratch_366 - scratch_263*scratch_367 - scratch_264*scratch_368 - scratch_265*scratch_369 - scratch_266*scratch_370 - scratch_267*scratch_371 - scratch_268*scratch_372);
    amrex::Real scratch_374 = scratch_302*(-A[320]*scratch_359 - scratch_288*scratch_360 - scratch_289*scratch_361 - scratch_290*scratch_362 - scratch_291*scratch_363 - scratch_292*scratch_364 - scratch_293*scratch_365 - scratch_294*scratch_366 - scratch_295*scratch_367 - scratch_296*scratch_368 - scratch_297*scratch_369 - scratch_298*scratch_370 - scratch_299*scratch_371 - scratch_300*scratch_372 - scratch_301*scratch_373);
    amrex::Real scratch_375 = scratch_338*(-A[331]*scratch_358 + A[343] - scratch_322*scratch_359 - scratch_323*scratch_360 - scratch_324*scratch_361 - scratch_325*scratch_362 - scratch_326*scratch_363 - scratch_327*scratch_364 - scratch_328*scratch_365 - scratch_329*scratch_366 - scratch_330*scratch_367 - scratch_331*scratch_368 - scratch_332*scratch_369 - scratch_333*scratch_370 - scratch_334*scratch_371 - scratch_335*scratch_372 - scratch_336*scratch_373 - scratch_337*scratch_374);
    amrex::Real scratch_376 = 1.0/(-A[345]*scratch_358 + A[353] - scratch_340*scratch_359 - scratch_341*scratch_360 - scratch_342*scratch_361 - scratch_343*scratch_362 - scratch_344*scratch_363 - scratch_345*scratch_364 - scratch_346*scratch_365 - scratch_347*scratch_366 - scratch_348*scratch_367 - scratch_349*scratch_368 - scratch_350*scratch_369 - scratch_351*scratch_370 - scratch_352*scratch_371 - scratch_353*scratch_372 - scratch_354*scratch_373 - scratch_355*scratch_374 - scratch_356*scratch_375);
    amrex::Real scratch_377 = -scratch_120*scratch_368 - scratch_143*scratch_369 - scratch_15*scratch_361 - scratch_168*scratch_370 - scratch_195*scratch_371 - scratch_224*scratch_372 - scratch_24*scratch_362 - scratch_255*scratch_373 - scratch_287*scratch_374 - scratch_321*scratch_375 - scratch_35*scratch_363 + scratch_358 - scratch_359*scratch_4 - scratch_360*scratch_7 - scratch_364*scratch_48 - scratch_365*scratch_63 - scratch_366*scratch_80 - scratch_367*scratch_99;
    amrex::Real scratch_378 = scratch_376*scratch_377;
    amrex::Real scratch_379 = b[19]/A[355];
    amrex::Real scratch_380 = A[20]*scratch_0;
    amrex::Real scratch_381 = scratch_5*(-A[22]*scratch_380 + A[42]);
    amrex::Real scratch_382 = scratch_9*(-A[44]*scratch_380 + A[64] - scratch_381*scratch_8);
    amrex::Real scratch_383 = scratch_18*(-A[66]*scratch_380 + A[86] - scratch_16*scratch_381 - scratch_17*scratch_382);
    amrex::Real scratch_384 = scratch_28*(A[108] - A[88]*scratch_380 - scratch_25*scratch_381 - scratch_26*scratch_382 - scratch_27*scratch_383);
    amrex::Real scratch_385 = scratch_40*(-A[110]*scratch_380 + A[130] - scratch_36*scratch_381 - scratch_37*scratch_382 - scratch_38*scratch_383 - scratch_384*scratch_39);
    amrex::Real scratch_386 = scratch_54*(-A[132]*scratch_380 + A[152] - scratch_381*scratch_49 - scratch_382*scratch_50 - scratch_383*scratch_51 - scratch_384*scratch_52 - scratch_385*scratch_53);
    amrex::Real scratch_387 = scratch_70*(-A[154]*scratch_380 + A[171] - scratch_381*scratch_64 - scratch_382*scratch_65 - scratch_383*scratch_66 - scratch_384*scratch_67 - scratch_385*scratch_68 - scratch_386*scratch_69);
    amrex::Real scratch_388 = scratch_88*(-A[173]*scratch_380 + A[189] - scratch_381*scratch_81 - scratch_382*scratch_82 - scratch_383*scratch_83 - scratch_384*scratch_84 - scratch_385*scratch_85 - scratch_386*scratch_86 - scratch_387*scratch_87);
    amrex::Real scratch_389 = scratch_108*(-A[191]*scratch_380 - scratch_100*scratch_381 - scratch_101*scratch_382 - scratch_102*scratch_383 - scratch_103*scratch_384 - scratch_104*scratch_385 - scratch_105*scratch_386 - scratch_106*scratch_387 - scratch_107*scratch_388);
    amrex::Real scratch_390 = scratch_130*(-A[210]*scratch_380 - scratch_121*scratch_381 - scratch_122*scratch_382 - scratch_123*scratch_383 - scratch_124*scratch_384 - scratch_125*scratch_385 - scratch_126*scratch_386 - scratch_127*scratch_387 - scratch_128*scratch_388 - scratch_129*scratch_389);
    amrex::Real scratch_391 = scratch_154*(-A[227]*scratch_380 + A[247] - scratch_144*scratch_381 - scratch_145*scratch_382 - scratch_146*scratch_383 - scratch_147*scratch_384 - scratch_148*scratch_385 - scratch_149*scratch_386 - scratch_150*scratch_387 - scratch_151*scratch_388 - scratch_152*scratch_389 - scratch_153*scratch_390);
    amrex::Real scratch_392 = scratch_180*(-A[249]*scratch_380 + A[269] - scratch_169*scratch_381 - scratch_170*scratch_382 - scratch_171*scratch_383 - scratch_172*scratch_384 - scratch_173*scratch_385 - scratch_174*scratch_386 - scratch_175*scratch_387 - scratch_176*scratch_388 - scratch_177*scratch_389 - scratch_178*scratch_390 - scratch_179*scratch_391);
    amrex::Real scratch_393 = scratch_208*(-A[271]*scratch_380 + A[291] - scratch_196*scratch_381 - scratch_197*scratch_382 - scratch_198*scratch_383 - scratch_199*scratch_384 - scratch_200*scratch_385 - scratch_201*scratch_386 - scratch_202*scratch_387 - scratch_203*scratch_388 - scratch_204*scratch_389 - scratch_205*scratch_390 - scratch_206*scratch_391 - scratch_207*scratch_392);
    amrex::Real scratch_394 = scratch_238*(-A[293]*scratch_380 - scratch_225*scratch_381 - scratch_226*scratch_382 - scratch_227*scratch_383 - scratch_228*scratch_384 - scratch_229*scratch_385 - scratch_230*scratch_386 - scratch_231*scratch_387 - scratch_232*scratch_388 - scratch_233*scratch_389 - scratch_234*scratch_390 - scratch_235*scratch_391 - scratch_236*scratch_392 - scratch_237*scratch_393);
    amrex::Real scratch_395 = scratch_269*(-A[311]*scratch_381 - scratch_256*scratch_382 - scratch_257*scratch_383 - scratch_258*scratch_384 - scratch_259*scratch_385 - scratch_260*scratch_386 - scratch_261*scratch_387 - scratch_262*scratch_388 - scratch_263*scratch_389 - scratch_264*scratch_390 - scratch_265*scratch_391 - scratch_266*scratch_392 - scratch_267*scratch_393 - scratch_268*scratch_394);
    amrex::Real scratch_396 = scratch_302*(-A[320]*scratch_381 - scratch_288*scratch_382 - scratch_289*scratch_383 - scratch_290*scratch_384 - scratch_291*scratch_385 - scratch_292*scratch_386 - scratch_293*scratch_387 - scratch_294*scratch_388 - scratch_295*scratch_389 - scratch_296*scratch_390 - scratch_297*scratch_391 - scratch_298*scratch_392 - scratch_299*scratch_393 - scratch_300*scratch_394 - scratch_301*scratch_395);
    amrex::Real scratch_397 = scratch_338*(-A[331]*scratch_380 - scratch_322*scratch_381 - scratch_323*scratch_382 - scratch_324*scratch_383 - scratch_325*scratch_384 - scratch_326*scratch_385 - scratch_327*scratch_386 - scratch_328*scratch_387 - scratch_329*scratch_388 - scratch_330*scratch_389 - scratch_331*scratch_390 - scratch_332*scratch_391 - scratch_333*scratch_392 - scratch_334*scratch_393 - scratch_335*scratch_394 - scratch_336*scratch_395 - scratch_337*scratch_396);
    amrex::Real scratch_398 = -A[345]*scratch_380 - scratch_340*scratch_381 - scratch_341*scratch_382 - scratch_342*scratch_383 - scratch_343*scratch_384 - scratch_344*scratch_385 - scratch_345*scratch_386 - scratch_346*scratch_387 - scratch_347*scratch_388 - scratch_348*scratch_389 - scratch_349*scratch_390 - scratch_350*scratch_391 - scratch_351*scratch_392 - scratch_352*scratch_393 - scratch_353*scratch_394 - scratch_354*scratch_395 - scratch_355*scratch_396 - scratch_356*scratch_397;
    amrex::Real scratch_399 = b[20]/A[356];
    amrex::Real scratch_400 = scratch_5*(-A[22]*scratch_1 + b[1]);
    amrex::Real scratch_401 = scratch_9*(-A[44]*scratch_1 + b[2] - scratch_400*scratch_8);
    amrex::Real scratch_402 = scratch_18*(-A[66]*scratch_1 + b[3] - scratch_16*scratch_400 - scratch_17*scratch_401);
    amrex::Real scratch_403 = scratch_28*(-A[88]*scratch_1 + b[4] - scratch_25*scratch_400 - scratch_26*scratch_401 - scratch_27*scratch_402);
    amrex::Real scratch_404 = scratch_40*(-A[110]*scratch_1 + b[5] - scratch_36*scratch_400 - scratch_37*scratch_401 - scratch_38*scratch_402 - scratch_39*scratch_403);
    amrex::Real scratch_405 = scratch_54*(-A[132]*scratch_1 + b[6] - scratch_400*scratch_49 - scratch_401*scratch_50 - scratch_402*scratch_51 - scratch_403*scratch_52 - scratch_404*scratch_53);
    amrex::Real scratch_406 = scratch_70*(-A[154]*scratch_1 + b[7] - scratch_400*scratch_64 - scratch_401*scratch_65 - scratch_402*scratch_66 - scratch_403*scratch_67 - scratch_404*scratch_68 - scratch_405*scratch_69);
    amrex::Real scratch_407 = scratch_88*(-A[173]*scratch_1 + b[8] - scratch_400*scratch_81 - scratch_401*scratch_82 - scratch_402*scratch_83 - scratch_403*scratch_84 - scratch_404*scratch_85 - scratch_405*scratch_86 - scratch_406*scratch_87);
    amrex::Real scratch_408 = scratch_108*(-A[191]*scratch_1 + b[9] - scratch_100*scratch_400 - scratch_101*scratch_401 - scratch_102*scratch_402 - scratch_103*scratch_403 - scratch_104*scratch_404 - scratch_105*scratch_405 - scratch_106*scratch_406 - scratch_107*scratch_407);
    amrex::Real scratch_409 = scratch_130*(-A[210]*scratch_1 + b[10] - scratch_121*scratch_400 - scratch_122*scratch_401 - scratch_123*scratch_402 - scratch_124*scratch_403 - scratch_125*scratch_404 - scratch_126*scratch_405 - scratch_127*scratch_406 - scratch_128*scratch_407 - scratch_129*scratch_408);
    amrex::Real scratch_410 = scratch_154*(-A[227]*scratch_1 + b[11] - scratch_144*scratch_400 - scratch_145*scratch_401 - scratch_146*scratch_402 - scratch_147*scratch_403 - scratch_148*scratch_404 - scratch_149*scratch_405 - scratch_150*scratch_406 - scratch_151*scratch_407 - scratch_152*scratch_408 - scratch_153*scratch_409);
    amrex::Real scratch_411 = scratch_180*(-A[249]*scratch_1 + b[12] - scratch_169*scratch_400 - scratch_170*scratch_401 - scratch_171*scratch_402 - scratch_172*scratch_403 - scratch_173*scratch_404 - scratch_174*scratch_405 - scratch_175*scratch_406 - scratch_176*scratch_407 - scratch_177*scratch_408 - scratch_178*scratch_409 - scratch_179*scratch_410);
    amrex::Real scratch_412 = scratch_208*(-A[271]*scratch_1 + b[13] - scratch_196*scratch_400 - scratch_197*scratch_401 - scratch_198*scratch_402 - scratch_199*scratch_403 - scratch_200*scratch_404 - scratch_201*scratch_405 - scratch_202*scratch_406 - scratch_203*scratch_407 - scratch_204*scratch_408 - scratch_205*scratch_409 - scratch_206*scratch_410 - scratch_207*scratch_411);
    amrex::Real scratch_413 = scratch_238*(-A[293]*scratch_1 + b[14] - scratch_225*scratch_400 - scratch_226*scratch_401 - scratch_227*scratch_402 - scratch_228*scratch_403 - scratch_229*scratch_404 - scratch_230*scratch_405 - scratch_231*scratch_406 - scratch_232*scratch_407 - scratch_233*scratch_408 - scratch_234*scratch_409 - scratch_235*scratch_410 - scratch_236*scratch_411 - scratch_237*scratch_412);
    amrex::Real scratch_414 = scratch_269*(-A[311]*scratch_400 + b[15] - scratch_256*scratch_401 - scratch_257*scratch_402 - scratch_258*scratch_403 - scratch_259*scratch_404 - scratch_260*scratch_405 - scratch_261*scratch_406 - scratch_262*scratch_407 - scratch_263*scratch_408 - scratch_264*scratch_409 - scratch_265*scratch_410 - scratch_266*scratch_411 - scratch_267*scratch_412 - scratch_268*scratch_413);
    amrex::Real scratch_415 = scratch_302*(-A[320]*scratch_400 + b[16] - scratch_288*scratch_401 - scratch_289*scratch_402 - scratch_290*scratch_403 - scratch_291*scratch_404 - scratch_292*scratch_405 - scratch_293*scratch_406 - scratch_294*scratch_407 - scratch_295*scratch_408 - scratch_296*scratch_409 - scratch_297*scratch_410 - scratch_298*scratch_411 - scratch_299*scratch_412 - scratch_300*scratch_413 - scratch_301*scratch_414);
    amrex::Real scratch_416 = scratch_338*(-A[331]*scratch_1 + b[17] - scratch_322*scratch_400 - scratch_323*scratch_401 - scratch_324*scratch_402 - scratch_325*scratch_403 - scratch_326*scratch_404 - scratch_327*scratch_405 - scratch_328*scratch_406 - scratch_329*scratch_407 - scratch_330*scratch_408 - scratch_331*scratch_409 - scratch_332*scratch_410 - scratch_333*scratch_411 - scratch_334*scratch_412 - scratch_335*scratch_413 - scratch_336*scratch_414 - scratch_337*scratch_415);
    amrex::Real scratch_417 = scratch_376*(-A[345]*scratch_1 + b[18] - scratch_340*scratch_400 - scratch_341*scratch_401 - scratch_342*scratch_402 - scratch_343*scratch_403 - scratch_344*scratch_404 - scratch_345*scratch_405 - scratch_346*scratch_406 - scratch_347*scratch_407 - scratch_348*scratch_408 - scratch_349*scratch_409 - scratch_350*scratch_410 - scratch_351*scratch_411 - scratch_352*scratch_412 - scratch_353*scratch_413 - scratch_354*scratch_414 - scratch_355*scratch_415 - scratch_356*scratch_416);
    amrex::Real scratch_418 = A[21]*scratch_0;
    amrex::Real scratch_419 = scratch_5*(-A[22]*scratch_418 + A[43]);
    amrex::Real scratch_420 = scratch_9*(-A[44]*scratch_418 + A[65] - scratch_419*scratch_8);
    amrex::Real scratch_421 = scratch_18*(-A[66]*scratch_418 + A[87] - scratch_16*scratch_419 - scratch_17*scratch_420);
    amrex::Real scratch_422 = scratch_28*(A[109] - A[88]*scratch_418 - scratch_25*scratch_419 - scratch_26*scratch_420 - scratch_27*scratch_421);
    amrex::Real scratch_423 = scratch_40*(-A[110]*scratch_418 + A[131] - scratch_36*scratch_419 - scratch_37*scratch_420 - scratch_38*scratch_421 - scratch_39*scratch_422);
    amrex::Real scratch_424 = scratch_54*(-A[132]*scratch_418 + A[153] - scratch_419*scratch_49 - scratch_420*scratch_50 - scratch_421*scratch_51 - scratch_422*scratch_52 - scratch_423*scratch_53);
    amrex::Real scratch_425 = scratch_70*(-A[154]*scratch_418 + A[172] - scratch_419*scratch_64 - scratch_420*scratch_65 - scratch_421*scratch_66 - scratch_422*scratch_67 - scratch_423*scratch_68 - scratch_424*scratch_69);
    amrex::Real scratch_426 = scratch_88*(-A[173]*scratch_418 + A[190] - scratch_419*scratch_81 - scratch_420*scratch_82 - scratch_421*scratch_83 - scratch_422*scratch_84 - scratch_423*scratch_85 - scratch_424*scratch_86 - scratch_425*scratch_87);
    amrex::Real scratch_427 = scratch_108*(-A[191]*scratch_418 + A[209] - scratch_100*scratch_419 - scratch_101*scratch_420 - scratch_102*scratch_421 - scratch_103*scratch_422 - scratch_104*scratch_423 - scratch_105*scratch_424 - scratch_106*scratch_425 - scratch_107*scratch_426);
    amrex::Real scratch_428 = scratch_130*(-A[210]*scratch_418 + A[226] - scratch_121*scratch_419 - scratch_122*scratch_420 - scratch_123*scratch_421 - scratch_124*scratch_422 - scratch_125*scratch_423 - scratch_126*scratch_424 - scratch_127*scratch_425 - scratch_128*scratch_426 - scratch_129*scratch_427);
    amrex::Real scratch_429 = scratch_154*(-A[227]*scratch_418 + A[248] - scratch_144*scratch_419 - scratch_145*scratch_420 - scratch_146*scratch_421 - scratch_147*scratch_422 - scratch_148*scratch_423 - scratch_149*scratch_424 - scratch_150*scratch_425 - scratch_151*scratch_426 - scratch_152*scratch_427 - scratch_153*scratch_428);
    amrex::Real scratch_430 = scratch_180*(-A[249]*scratch_418 + A[270] - scratch_169*scratch_419 - scratch_170*scratch_420 - scratch_171*scratch_421 - scratch_172*scratch_422 - scratch_173*scratch_423 - scratch_174*scratch_424 - scratch_175*scratch_425 - scratch_176*scratch_426 - scratch_177*scratch_427 - scratch_178*scratch_428 - scratch_179*scratch_429);
    amrex::Real scratch_431 = scratch_208*(-A[271]*scratch_418 + A[292] - scratch_196*scratch_419 - scratch_197*scratch_420 - scratch_198*scratch_421 - scratch_199*scratch_422 - scratch_200*scratch_423 - scratch_201*scratch_424 - scratch_202*scratch_425 - scratch_203*scratch_426 - scratch_204*scratch_427 - scratch_205*scratch_428 - scratch_206*scratch_429 - scratch_207*scratch_430);
    amrex::Real scratch_432 = scratch_238*(-A[293]*scratch_418 + A[310] - scratch_225*scratch_419 - scratch_226*scratch_420 - scratch_227*scratch_421 - scratch_228*scratch_422 - scratch_229*scratch_423 - scratch_230*scratch_424 - scratch_231*scratch_425 - scratch_232*scratch_426 - scratch_233*scratch_427 - scratch_234*scratch_428 - scratch_235*scratch_429 - scratch_236*scratch_430 - scratch_237*scratch_431);
    amrex::Real scratch_433 = scratch_269*(-A[311]*scratch_419 + A[319] - scratch_256*scratch_420 - scratch_257*scratch_421 - scratch_258*scratch_422 - scratch_259*scratch_423 - scratch_260*scratch_424 - scratch_261*scratch_425 - scratch_262*scratch_426 - scratch_263*scratch_427 - scratch_264*scratch_428 - scratch_265*scratch_429 - scratch_266*scratch_430 - scratch_267*scratch_431 - scratch_268*scratch_432);
    amrex::Real scratch_434 = scratch_302*(-A[320]*scratch_419 + A[330] - scratch_288*scratch_420 - scratch_289*scratch_421 - scratch_290*scratch_422 - scratch_291*scratch_423 - scratch_292*scratch_424 - scratch_293*scratch_425 - scratch_294*scratch_426 - scratch_295*scratch_427 - scratch_296*scratch_428 - scratch_297*scratch_429 - scratch_298*scratch_430 - scratch_299*scratch_431 - scratch_300*scratch_432 - scratch_301*scratch_433);
    amrex::Real scratch_435 = scratch_338*(-A[331]*scratch_418 + A[344] - scratch_322*scratch_419 - scratch_323*scratch_420 - scratch_324*scratch_421 - scratch_325*scratch_422 - scratch_326*scratch_423 - scratch_327*scratch_424 - scratch_328*scratch_425 - scratch_329*scratch_426 - scratch_330*scratch_427 - scratch_331*scratch_428 - scratch_332*scratch_429 - scratch_333*scratch_430 - scratch_334*scratch_431 - scratch_335*scratch_432 - scratch_336*scratch_433 - scratch_337*scratch_434);
    amrex::Real scratch_436 = -A[345]*scratch_418 + A[354] - scratch_340*scratch_419 - scratch_341*scratch_420 - scratch_342*scratch_421 - scratch_343*scratch_422 - scratch_344*scratch_423 - scratch_345*scratch_424 - scratch_346*scratch_425 - scratch_347*scratch_426 - scratch_348*scratch_427 - scratch_349*scratch_428 - scratch_350*scratch_429 - scratch_351*scratch_430 - scratch_352*scratch_431 - scratch_353*scratch_432 - scratch_354*scratch_433 - scratch_355*scratch_434 - scratch_356*scratch_435;
    amrex::Real scratch_437 = -A[357]*scratch_4 + A[358];
    amrex::Real scratch_438 = -A[357]*scratch_3 + A[359] - scratch_437*scratch_6;
    amrex::Real scratch_439 = -A[357]*scratch_12 + A[360] - scratch_13*scratch_437 - scratch_14*scratch_438;
    amrex::Real scratch_440 = -A[357]*scratch_20 + A[361] - scratch_21*scratch_437 - scratch_22*scratch_438 - scratch_23*scratch_439;
    amrex::Real scratch_441 = -A[357]*scratch_30 + A[362] - scratch_31*scratch_437 - scratch_32*scratch_438 - scratch_33*scratch_439 - scratch_34*scratch_440;
    amrex::Real scratch_442 = -A[357]*scratch_42 + A[363] - scratch_43*scratch_437 - scratch_438*scratch_44 - scratch_439*scratch_45 - scratch_440*scratch_46 - scratch_441*scratch_47;
    amrex::Real scratch_443 = -A[357]*scratch_56 + A[364] - scratch_437*scratch_57 - scratch_438*scratch_58 - scratch_439*scratch_59 - scratch_440*scratch_60 - scratch_441*scratch_61 - scratch_442*scratch_62;
    amrex::Real scratch_444 = -A[357]*scratch_72 + A[365] - scratch_437*scratch_73 - scratch_438*scratch_74 - scratch_439*scratch_75 - scratch_440*scratch_76 - scratch_441*scratch_77 - scratch_442*scratch_78 - scratch_443*scratch_79;
    amrex::Real scratch_445 = -A[357]*scratch_90 + A[366] - scratch_437*scratch_91 - scratch_438*scratch_92 - scratch_439*scratch_93 - scratch_440*scratch_94 - scratch_441*scratch_95 - scratch_442*scratch_96 - scratch_443*scratch_97 - scratch_444*scratch_98;
    amrex::Real scratch_446 = -A[357]*scratch_110 + A[367] - scratch_111*scratch_437 - scratch_112*scratch_438 - scratch_113*scratch_439 - scratch_114*scratch_440 - scratch_115*scratch_441 - scratch_116*scratch_442 - scratch_117*scratch_443 - scratch_118*scratch_444 - scratch_119*scratch_445;
    amrex::Real scratch_447 = -A[357]*scratch_132 + A[368] - scratch_133*scratch_437 - scratch_134*scratch_438 - scratch_135*scratch_439 - scratch_136*scratch_440 - scratch_137*scratch_441 - scratch_138*scratch_442 - scratch_139*scratch_443 - scratch_140*scratch_444 - scratch_141*scratch_445 - scratch_142*scratch_446;
    amrex::Real scratch_448 = -A[357]*scratch_156 + A[369] - scratch_157*scratch_437 - scratch_158*scratch_438 - scratch_159*scratch_439 - scratch_160*scratch_440 - scratch_161*scratch_441 - scratch_162*scratch_442 - scratch_163*scratch_443 - scratch_164*scratch_444 - scratch_165*scratch_445 - scratch_166*scratch_446 - scratch_167*scratch_447;
    amrex::Real scratch_449 = -A[357]*scratch_182 + A[370] - scratch_183*scratch_437 - scratch_184*scratch_438 - scratch_185*scratch_439 - scratch_186*scratch_440 - scratch_187*scratch_441 - scratch_188*scratch_442 - scratch_189*scratch_443 - scratch_190*scratch_444 - scratch_191*scratch_445 - scratch_192*scratch_446 - scratch_193*scratch_447 - scratch_194*scratch_448;
    amrex::Real scratch_450 = -A[357]*scratch_210 + A[371] - scratch_211*scratch_437 - scratch_212*scratch_438 - scratch_213*scratch_439 - scratch_214*scratch_440 - scratch_215*scratch_441 - scratch_216*scratch_442 - scratch_217*scratch_443 - scratch_218*scratch_444 - scratch_219*scratch_445 - scratch_220*scratch_446 - scratch_221*scratch_447 - scratch_222*scratch_448 - scratch_223*scratch_449;
    amrex::Real scratch_451 = -A[357]*scratch_240 + A[372] - scratch_241*scratch_437 - scratch_242*scratch_438 - scratch_243*scratch_439 - scratch_244*scratch_440 - scratch_245*scratch_441 - scratch_246*scratch_442 - scratch_247*scratch_443 - scratch_248*scratch_444 - scratch_249*scratch_445 - scratch_250*scratch_446 - scratch_251*scratch_447 - scratch_252*scratch_448 - scratch_253*scratch_449 - scratch_254*scratch_450;
    amrex::Real scratch_452 = -A[357]*scratch_271 + A[373] - scratch_272*scratch_437 - scratch_273*scratch_438 - scratch_274*scratch_439 - scratch_275*scratch_440 - scratch_276*scratch_441 - scratch_277*scratch_442 - scratch_278*scratch_443 - scratch_279*scratch_444 - scratch_280*scratch_445 - scratch_281*scratch_446 - scratch_282*scratch_447 - scratch_283*scratch_448 - scratch_284*scratch_449 - scratch_285*scratch_450 - scratch_286*scratch_451;
    amrex::Real scratch_453 = -A[357]*scratch_304 + A[374] - scratch_305*scratch_437 - scratch_306*scratch_438 - scratch_307*scratch_439 - scratch_308*scratch_440 - scratch_309*scratch_441 - scratch_310*scratch_442 - scratch_311*scratch_443 - scratch_312*scratch_444 - scratch_313*scratch_445 - scratch_314*scratch_446 - scratch_315*scratch_447 - scratch_316*scratch_448 - scratch_317*scratch_449 - scratch_318*scratch_450 - scratch_319*scratch_451 - scratch_320*scratch_452;
    amrex::Real scratch_454 = -A[357]*scratch_358 + A[375] - scratch_359*scratch_437 - scratch_360*scratch_438 - scratch_361*scratch_439 - scratch_362*scratch_440 - scratch_363*scratch_441 - scratch_364*scratch_442 - scratch_365*scratch_443 - scratch_366*scratch_444 - scratch_367*scratch_445 - scratch_368*scratch_446 - scratch_369*scratch_447 - scratch_370*scratch_448 - scratch_371*scratch_449 - scratch_372*scratch_450 - scratch_373*scratch_451 - scratch_374*scratch_452 - scratch_375*scratch_453;
    amrex::Real scratch_455 = scratch_376*scratch_454;
    amrex::Real scratch_456 = (-A[357]*scratch_1 + b[21] - scratch_379*(-A[357]*scratch_2 + A[376] - scratch_10*scratch_437 - scratch_109*scratch_445 - scratch_11*scratch_438 - scratch_131*scratch_446 - scratch_155*scratch_447 - scratch_181*scratch_448 - scratch_19*scratch_439 - scratch_209*scratch_449 - scratch_239*scratch_450 - scratch_270*scratch_451 - scratch_29*scratch_440 - scratch_303*scratch_452 - scratch_339*scratch_453 - scratch_357*scratch_455 - scratch_41*scratch_441 - scratch_442*scratch_55 - scratch_443*scratch_71 - scratch_444*scratch_89) - scratch_399*(-A[357]*scratch_380 + A[377] - scratch_381*scratch_437 - scratch_382*scratch_438 - scratch_383*scratch_439 - scratch_384*scratch_440 - scratch_385*scratch_441 - scratch_386*scratch_442 - scratch_387*scratch_443 - scratch_388*scratch_444 - scratch_389*scratch_445 - scratch_390*scratch_446 - scratch_391*scratch_447 - scratch_392*scratch_448 - scratch_393*scratch_449 - scratch_394*scratch_450 - scratch_395*scratch_451 - scratch_396*scratch_452 - scratch_397*scratch_453 - scratch_398*scratch_455) - scratch_400*scratch_437 - scratch_401*scratch_438 - scratch_402*scratch_439 - scratch_403*scratch_440 - scratch_404*scratch_441 - scratch_405*scratch_442 - scratch_406*scratch_443 - scratch_407*scratch_444 - scratch_408*scratch_445 - scratch_409*scratch_446 - scratch_410*scratch_447 - scratch_411*scratch_448 - scratch_412*scratch_449 - scratch_413*scratch_450 - scratch_414*scratch_451 - scratch_415*scratch_452 - scratch_416*scratch_453 - scratch_417*scratch_454)/(-A[357]*scratch_418 + A[378] - scratch_419*scratch_437 - scratch_420*scratch_438 - scratch_421*scratch_439 - scratch_422*scratch_440 - scratch_423*scratch_441 - scratch_424*scratch_442 - scratch_425*scratch_443 - scratch_426*scratch_444 - scratch_427*scratch_445 - scratch_428*scratch_446 - scratch_429*scratch_447 - scratch_430*scratch_448 - scratch_431*scratch_449 - scratch_432*scratch_450 - scratch_433*scratch_451 - scratch_434*scratch_452 - scratch_435*scratch_453 - scratch_436*scratch_455);
    amrex::Real scratch_457 = scratch_13 - scratch_14*scratch_6;
    amrex::Real scratch_458 = scratch_21 - scratch_22*scratch_6 - scratch_23*scratch_457;
    amrex::Real scratch_459 = scratch_31 - scratch_32*scratch_6 - scratch_33*scratch_457 - scratch_34*scratch_458;
    amrex::Real scratch_460 = scratch_43 - scratch_44*scratch_6 - scratch_45*scratch_457 - scratch_458*scratch_46 - scratch_459*scratch_47;
    amrex::Real scratch_461 = -scratch_457*scratch_59 - scratch_458*scratch_60 - scratch_459*scratch_61 - scratch_460*scratch_62 + scratch_57 - scratch_58*scratch_6;
    amrex::Real scratch_462 = -scratch_457*scratch_75 - scratch_458*scratch_76 - scratch_459*scratch_77 - scratch_460*scratch_78 - scratch_461*scratch_79 - scratch_6*scratch_74 + scratch_73;
    amrex::Real scratch_463 = -scratch_457*scratch_93 - scratch_458*scratch_94 - scratch_459*scratch_95 - scratch_460*scratch_96 - scratch_461*scratch_97 - scratch_462*scratch_98 - scratch_6*scratch_92 + scratch_91;
    amrex::Real scratch_464 = scratch_111 - scratch_112*scratch_6 - scratch_113*scratch_457 - scratch_114*scratch_458 - scratch_115*scratch_459 - scratch_116*scratch_460 - scratch_117*scratch_461 - scratch_118*scratch_462 - scratch_119*scratch_463;
    amrex::Real scratch_465 = scratch_133 - scratch_134*scratch_6 - scratch_135*scratch_457 - scratch_136*scratch_458 - scratch_137*scratch_459 - scratch_138*scratch_460 - scratch_139*scratch_461 - scratch_140*scratch_462 - scratch_141*scratch_463 - scratch_142*scratch_464;
    amrex::Real scratch_466 = scratch_157 - scratch_158*scratch_6 - scratch_159*scratch_457 - scratch_160*scratch_458 - scratch_161*scratch_459 - scratch_162*scratch_460 - scratch_163*scratch_461 - scratch_164*scratch_462 - scratch_165*scratch_463 - scratch_166*scratch_464 - scratch_167*scratch_465;
    amrex::Real scratch_467 = scratch_183 - scratch_184*scratch_6 - scratch_185*scratch_457 - scratch_186*scratch_458 - scratch_187*scratch_459 - scratch_188*scratch_460 - scratch_189*scratch_461 - scratch_190*scratch_462 - scratch_191*scratch_463 - scratch_192*scratch_464 - scratch_193*scratch_465 - scratch_194*scratch_466;
    amrex::Real scratch_468 = scratch_211 - scratch_212*scratch_6 - scratch_213*scratch_457 - scratch_214*scratch_458 - scratch_215*scratch_459 - scratch_216*scratch_460 - scratch_217*scratch_461 - scratch_218*scratch_462 - scratch_219*scratch_463 - scratch_220*scratch_464 - scratch_221*scratch_465 - scratch_222*scratch_466 - scratch_223*scratch_467;
    amrex::Real scratch_469 = scratch_241 - scratch_242*scratch_6 - scratch_243*scratch_457 - scratch_244*scratch_458 - scratch_245*scratch_459 - scratch_246*scratch_460 - scratch_247*scratch_461 - scratch_248*scratch_462 - scratch_249*scratch_463 - scratch_250*scratch_464 - scratch_251*scratch_465 - scratch_252*scratch_466 - scratch_253*scratch_467 - scratch_254*scratch_468;
    amrex::Real scratch_470 = scratch_272 - scratch_273*scratch_6 - scratch_274*scratch_457 - scratch_275*scratch_458 - scratch_276*scratch_459 - scratch_277*scratch_460 - scratch_278*scratch_461 - scratch_279*scratch_462 - scratch_280*scratch_463 - scratch_281*scratch_464 - scratch_282*scratch_465 - scratch_283*scratch_466 - scratch_284*scratch_467 - scratch_285*scratch_468 - scratch_286*scratch_469;
    amrex::Real scratch_471 = scratch_305 - scratch_306*scratch_6 - scratch_307*scratch_457 - scratch_308*scratch_458 - scratch_309*scratch_459 - scratch_310*scratch_460 - scratch_311*scratch_461 - scratch_312*scratch_462 - scratch_313*scratch_463 - scratch_314*scratch_464 - scratch_315*scratch_465 - scratch_316*scratch_466 - scratch_317*scratch_467 - scratch_318*scratch_468 - scratch_319*scratch_469 - scratch_320*scratch_470;
    amrex::Real scratch_472 = scratch_359 - scratch_360*scratch_6 - scratch_361*scratch_457 - scratch_362*scratch_458 - scratch_363*scratch_459 - scratch_364*scratch_460 - scratch_365*scratch_461 - scratch_366*scratch_462 - scratch_367*scratch_463 - scratch_368*scratch_464 - scratch_369*scratch_465 - scratch_370*scratch_466 - scratch_371*scratch_467 - scratch_372*scratch_468 - scratch_373*scratch_469 - scratch_374*scratch_470 - scratch_375*scratch_471;
    amrex::Real scratch_473 = scratch_376*scratch_472;
    amrex::Real scratch_474 = -scratch_14*scratch_23 + scratch_22;
    amrex::Real scratch_475 = -scratch_14*scratch_33 + scratch_32 - scratch_34*scratch_474;
    amrex::Real scratch_476 = -scratch_14*scratch_45 + scratch_44 - scratch_46*scratch_474 - scratch_47*scratch_475;
    amrex::Real scratch_477 = -scratch_14*scratch_59 - scratch_474*scratch_60 - scratch_475*scratch_61 - scratch_476*scratch_62 + scratch_58;
    amrex::Real scratch_478 = -scratch_14*scratch_75 - scratch_474*scratch_76 - scratch_475*scratch_77 - scratch_476*scratch_78 - scratch_477*scratch_79 + scratch_74;
    amrex::Real scratch_479 = -scratch_14*scratch_93 - scratch_474*scratch_94 - scratch_475*scratch_95 - scratch_476*scratch_96 - scratch_477*scratch_97 - scratch_478*scratch_98 + scratch_92;
    amrex::Real scratch_480 = scratch_112 - scratch_113*scratch_14 - scratch_114*scratch_474 - scratch_115*scratch_475 - scratch_116*scratch_476 - scratch_117*scratch_477 - scratch_118*scratch_478 - scratch_119*scratch_479;
    amrex::Real scratch_481 = scratch_134 - scratch_135*scratch_14 - scratch_136*scratch_474 - scratch_137*scratch_475 - scratch_138*scratch_476 - scratch_139*scratch_477 - scratch_140*scratch_478 - scratch_141*scratch_479 - scratch_142*scratch_480;
    amrex::Real scratch_482 = -scratch_14*scratch_159 + scratch_158 - scratch_160*scratch_474 - scratch_161*scratch_475 - scratch_162*scratch_476 - scratch_163*scratch_477 - scratch_164*scratch_478 - scratch_165*scratch_479 - scratch_166*scratch_480 - scratch_167*scratch_481;
    amrex::Real scratch_483 = -scratch_14*scratch_185 + scratch_184 - scratch_186*scratch_474 - scratch_187*scratch_475 - scratch_188*scratch_476 - scratch_189*scratch_477 - scratch_190*scratch_478 - scratch_191*scratch_479 - scratch_192*scratch_480 - scratch_193*scratch_481 - scratch_194*scratch_482;
    amrex::Real scratch_484 = -scratch_14*scratch_213 + scratch_212 - scratch_214*scratch_474 - scratch_215*scratch_475 - scratch_216*scratch_476 - scratch_217*scratch_477 - scratch_218*scratch_478 - scratch_219*scratch_479 - scratch_220*scratch_480 - scratch_221*scratch_481 - scratch_222*scratch_482 - scratch_223*scratch_483;
    amrex::Real scratch_485 = -scratch_14*scratch_243 + scratch_242 - scratch_244*scratch_474 - scratch_245*scratch_475 - scratch_246*scratch_476 - scratch_247*scratch_477 - scratch_248*scratch_478 - scratch_249*scratch_479 - scratch_250*scratch_480 - scratch_251*scratch_481 - scratch_252*scratch_482 - scratch_253*scratch_483 - scratch_254*scratch_484;
    amrex::Real scratch_486 = -scratch_14*scratch_274 + scratch_273 - scratch_275*scratch_474 - scratch_276*scratch_475 - scratch_277*scratch_476 - scratch_278*scratch_477 - scratch_279*scratch_478 - scratch_280*scratch_479 - scratch_281*scratch_480 - scratch_282*scratch_481 - scratch_283*scratch_482 - scratch_284*scratch_483 - scratch_285*scratch_484 - scratch_286*scratch_485;
    amrex::Real scratch_487 = -scratch_14*scratch_307 + scratch_306 - scratch_308*scratch_474 - scratch_309*scratch_475 - scratch_310*scratch_476 - scratch_311*scratch_477 - scratch_312*scratch_478 - scratch_313*scratch_479 - scratch_314*scratch_480 - scratch_315*scratch_481 - scratch_316*scratch_482 - scratch_317*scratch_483 - scratch_318*scratch_484 - scratch_319*scratch_485 - scratch_320*scratch_486;
    amrex::Real scratch_488 = -scratch_14*scratch_361 + scratch_360 - scratch_362*scratch_474 - scratch_363*scratch_475 - scratch_364*scratch_476 - scratch_365*scratch_477 - scratch_366*scratch_478 - scratch_367*scratch_479 - scratch_368*scratch_480 - scratch_369*scratch_481 - scratch_370*scratch_482 - scratch_371*scratch_483 - scratch_372*scratch_484 - scratch_373*scratch_485 - scratch_374*scratch_486 - scratch_375*scratch_487;
    amrex::Real scratch_489 = scratch_376*scratch_488;
    amrex::Real scratch_490 = -scratch_23*scratch_34 + scratch_33;
    amrex::Real scratch_491 = -scratch_23*scratch_46 + scratch_45 - scratch_47*scratch_490;
    amrex::Real scratch_492 = -scratch_23*scratch_60 - scratch_490*scratch_61 - scratch_491*scratch_62 + scratch_59;
    amrex::Real scratch_493 = -scratch_23*scratch_76 - scratch_490*scratch_77 - scratch_491*scratch_78 - scratch_492*scratch_79 + scratch_75;
    amrex::Real scratch_494 = -scratch_23*scratch_94 - scratch_490*scratch_95 - scratch_491*scratch_96 - scratch_492*scratch_97 - scratch_493*scratch_98 + scratch_93;
    amrex::Real scratch_495 = scratch_113 - scratch_114*scratch_23 - scratch_115*scratch_490 - scratch_116*scratch_491 - scratch_117*scratch_492 - scratch_118*scratch_493 - scratch_119*scratch_494;
    amrex::Real scratch_496 = scratch_135 - scratch_136*scratch_23 - scratch_137*scratch_490 - scratch_138*scratch_491 - scratch_139*scratch_492 - scratch_140*scratch_493 - scratch_141*scratch_494 - scratch_142*scratch_495;
    amrex::Real scratch_497 = scratch_159 - scratch_160*scratch_23 - scratch_161*scratch_490 - scratch_162*scratch_491 - scratch_163*scratch_492 - scratch_164*scratch_493 - scratch_165*scratch_494 - scratch_166*scratch_495 - scratch_167*scratch_496;
    amrex::Real scratch_498 = scratch_185 - scratch_186*scratch_23 - scratch_187*scratch_490 - scratch_188*scratch_491 - scratch_189*scratch_492 - scratch_190*scratch_493 - scratch_191*scratch_494 - scratch_192*scratch_495 - scratch_193*scratch_496 - scratch_194*scratch_497;
    amrex::Real scratch_499 = scratch_213 - scratch_214*scratch_23 - scratch_215*scratch_490 - scratch_216*scratch_491 - scratch_217*scratch_492 - scratch_218*scratch_493 - scratch_219*scratch_494 - scratch_220*scratch_495 - scratch_221*scratch_496 - scratch_222*scratch_497 - scratch_223*scratch_498;
    amrex::Real scratch_500 = -scratch_23*scratch_244 + scratch_243 - scratch_245*scratch_490 - scratch_246*scratch_491 - scratch_247*scratch_492 - scratch_248*scratch_493 - scratch_249*scratch_494 - scratch_250*scratch_495 - scratch_251*scratch_496 - scratch_252*scratch_497 - scratch_253*scratch_498 - scratch_254*scratch_499;
    amrex::Real scratch_501 = -scratch_23*scratch_275 + scratch_274 - scratch_276*scratch_490 - scratch_277*scratch_491 - scratch_278*scratch_492 - scratch_279*scratch_493 - scratch_280*scratch_494 - scratch_281*scratch_495 - scratch_282*scratch_496 - scratch_283*scratch_497 - scratch_284*scratch_498 - scratch_285*scratch_499 - scratch_286*scratch_500;
    amrex::Real scratch_502 = -scratch_23*scratch_308 + scratch_307 - scratch_309*scratch_490 - scratch_310*scratch_491 - scratch_311*scratch_492 - scratch_312*scratch_493 - scratch_313*scratch_494 - scratch_314*scratch_495 - scratch_315*scratch_496 - scratch_316*scratch_497 - scratch_317*scratch_498 - scratch_318*scratch_499 - scratch_319*scratch_500 - scratch_320*scratch_501;
    amrex::Real scratch_503 = -scratch_23*scratch_362 + scratch_361 - scratch_363*scratch_490 - scratch_364*scratch_491 - scratch_365*scratch_492 - scratch_366*scratch_493 - scratch_367*scratch_494 - scratch_368*scratch_495 - scratch_369*scratch_496 - scratch_370*scratch_497 - scratch_371*scratch_498 - scratch_372*scratch_499 - scratch_373*scratch_500 - scratch_374*scratch_501 - scratch_375*scratch_502;
    amrex::Real scratch_504 = scratch_376*scratch_503;
    amrex::Real scratch_505 = -scratch_34*scratch_47 + scratch_46;
    amrex::Real scratch_506 = -scratch_34*scratch_61 - scratch_505*scratch_62 + scratch_60;
    amrex::Real scratch_507 = -scratch_34*scratch_77 - scratch_505*scratch_78 - scratch_506*scratch_79 + scratch_76;
    amrex::Real scratch_508 = -scratch_34*scratch_95 - scratch_505*scratch_96 - scratch_506*scratch_97 - scratch_507*scratch_98 + scratch_94;
    amrex::Real scratch_509 = scratch_114 - scratch_115*scratch_34 - scratch_116*scratch_505 - scratch_117*scratch_506 - scratch_118*scratch_507 - scratch_119*scratch_508;
    amrex::Real scratch_510 = scratch_136 - scratch_137*scratch_34 - scratch_138*scratch_505 - scratch_139*scratch_506 - scratch_140*scratch_507 - scratch_141*scratch_508 - scratch_142*scratch_509;
    amrex::Real scratch_511 = scratch_160 - scratch_161*scratch_34 - scratch_162*scratch_505 - scratch_163*scratch_506 - scratch_164*scratch_507 - scratch_165*scratch_508 - scratch_166*scratch_509 - scratch_167*scratch_510;
    amrex::Real scratch_512 = scratch_186 - scratch_187*scratch_34 - scratch_188*scratch_505 - scratch_189*scratch_506 - scratch_190*scratch_507 - scratch_191*scratch_508 - scratch_192*scratch_509 - scratch_193*scratch_510 - scratch_194*scratch_511;
    amrex::Real scratch_513 = scratch_214 - scratch_215*scratch_34 - scratch_216*scratch_505 - scratch_217*scratch_506 - scratch_218*scratch_507 - scratch_219*scratch_508 - scratch_220*scratch_509 - scratch_221*scratch_510 - scratch_222*scratch_511 - scratch_223*scratch_512;
    amrex::Real scratch_514 = scratch_244 - scratch_245*scratch_34 - scratch_246*scratch_505 - scratch_247*scratch_506 - scratch_248*scratch_507 - scratch_249*scratch_508 - scratch_250*scratch_509 - scratch_251*scratch_510 - scratch_252*scratch_511 - scratch_253*scratch_512 - scratch_254*scratch_513;
    amrex::Real scratch_515 = scratch_275 - scratch_276*scratch_34 - scratch_277*scratch_505 - scratch_278*scratch_506 - scratch_279*scratch_507 - scratch_280*scratch_508 - scratch_281*scratch_509 - scratch_282*scratch_510 - scratch_283*scratch_511 - scratch_284*scratch_512 - scratch_285*scratch_513 - scratch_286*scratch_514;
    amrex::Real scratch_516 = scratch_308 - scratch_309*scratch_34 - scratch_310*scratch_505 - scratch_311*scratch_506 - scratch_312*scratch_507 - scratch_313*scratch_508 - scratch_314*scratch_509 - scratch_315*scratch_510 - scratch_316*scratch_511 - scratch_317*scratch_512 - scratch_318*scratch_513 - scratch_319*scratch_514 - scratch_320*scratch_515;
    amrex::Real scratch_517 = -scratch_34*scratch_363 + scratch_362 - scratch_364*scratch_505 - scratch_365*scratch_506 - scratch_366*scratch_507 - scratch_367*scratch_508 - scratch_368*scratch_509 - scratch_369*scratch_510 - scratch_370*scratch_511 - scratch_371*scratch_512 - scratch_372*scratch_513 - scratch_373*scratch_514 - scratch_374*scratch_515 - scratch_375*scratch_516;
    amrex::Real scratch_518 = scratch_376*scratch_517;
    amrex::Real scratch_519 = -scratch_47*scratch_62 + scratch_61;
    amrex::Real scratch_520 = -scratch_47*scratch_78 - scratch_519*scratch_79 + scratch_77;
    amrex::Real scratch_521 = -scratch_47*scratch_96 - scratch_519*scratch_97 - scratch_520*scratch_98 + scratch_95;
    amrex::Real scratch_522 = scratch_115 - scratch_116*scratch_47 - scratch_117*scratch_519 - scratch_118*scratch_520 - scratch_119*scratch_521;
    amrex::Real scratch_523 = scratch_137 - scratch_138*scratch_47 - scratch_139*scratch_519 - scratch_140*scratch_520 - scratch_141*scratch_521 - scratch_142*scratch_522;
    amrex::Real scratch_524 = scratch_161 - scratch_162*scratch_47 - scratch_163*scratch_519 - scratch_164*scratch_520 - scratch_165*scratch_521 - scratch_166*scratch_522 - scratch_167*scratch_523;
    amrex::Real scratch_525 = scratch_187 - scratch_188*scratch_47 - scratch_189*scratch_519 - scratch_190*scratch_520 - scratch_191*scratch_521 - scratch_192*scratch_522 - scratch_193*scratch_523 - scratch_194*scratch_524;
    amrex::Real scratch_526 = scratch_215 - scratch_216*scratch_47 - scratch_217*scratch_519 - scratch_218*scratch_520 - scratch_219*scratch_521 - scratch_220*scratch_522 - scratch_221*scratch_523 - scratch_222*scratch_524 - scratch_223*scratch_525;
    amrex::Real scratch_527 = scratch_245 - scratch_246*scratch_47 - scratch_247*scratch_519 - scratch_248*scratch_520 - scratch_249*scratch_521 - scratch_250*scratch_522 - scratch_251*scratch_523 - scratch_252*scratch_524 - scratch_253*scratch_525 - scratch_254*scratch_526;
    amrex::Real scratch_528 = scratch_276 - scratch_277*scratch_47 - scratch_278*scratch_519 - scratch_279*scratch_520 - scratch_280*scratch_521 - scratch_281*scratch_522 - scratch_282*scratch_523 - scratch_283*scratch_524 - scratch_284*scratch_525 - scratch_285*scratch_526 - scratch_286*scratch_527;
    amrex::Real scratch_529 = scratch_309 - scratch_310*scratch_47 - scratch_311*scratch_519 - scratch_312*scratch_520 - scratch_313*scratch_521 - scratch_314*scratch_522 - scratch_315*scratch_523 - scratch_316*scratch_524 - scratch_317*scratch_525 - scratch_318*scratch_526 - scratch_319*scratch_527 - scratch_320*scratch_528;
    amrex::Real scratch_530 = scratch_363 - scratch_364*scratch_47 - scratch_365*scratch_519 - scratch_366*scratch_520 - scratch_367*scratch_521 - scratch_368*scratch_522 - scratch_369*scratch_523 - scratch_370*scratch_524 - scratch_371*scratch_525 - scratch_372*scratch_526 - scratch_373*scratch_527 - scratch_374*scratch_528 - scratch_375*scratch_529;
    amrex::Real scratch_531 = scratch_376*scratch_530;
    amrex::Real scratch_532 = -scratch_62*scratch_79 + scratch_78;
    amrex::Real scratch_533 = -scratch_532*scratch_98 - scratch_62*scratch_97 + scratch_96;
    amrex::Real scratch_534 = scratch_116 - scratch_117*scratch_62 - scratch_118*scratch_532 - scratch_119*scratch_533;
    amrex::Real scratch_535 = scratch_138 - scratch_139*scratch_62 - scratch_140*scratch_532 - scratch_141*scratch_533 - scratch_142*scratch_534;
    amrex::Real scratch_536 = scratch_162 - scratch_163*scratch_62 - scratch_164*scratch_532 - scratch_165*scratch_533 - scratch_166*scratch_534 - scratch_167*scratch_535;
    amrex::Real scratch_537 = scratch_188 - scratch_189*scratch_62 - scratch_190*scratch_532 - scratch_191*scratch_533 - scratch_192*scratch_534 - scratch_193*scratch_535 - scratch_194*scratch_536;
    amrex::Real scratch_538 = scratch_216 - scratch_217*scratch_62 - scratch_218*scratch_532 - scratch_219*scratch_533 - scratch_220*scratch_534 - scratch_221*scratch_535 - scratch_222*scratch_536 - scratch_223*scratch_537;
    amrex::Real scratch_539 = scratch_246 - scratch_247*scratch_62 - scratch_248*scratch_532 - scratch_249*scratch_533 - scratch_250*scratch_534 - scratch_251*scratch_535 - scratch_252*scratch_536 - scratch_253*scratch_537 - scratch_254*scratch_538;
    amrex::Real scratch_540 = scratch_277 - scratch_278*scratch_62 - scratch_279*scratch_532 - scratch_280*scratch_533 - scratch_281*scratch_534 - scratch_282*scratch_535 - scratch_283*scratch_536 - scratch_284*scratch_537 - scratch_285*scratch_538 - scratch_286*scratch_539;
    amrex::Real scratch_541 = scratch_310 - scratch_311*scratch_62 - scratch_312*scratch_532 - scratch_313*scratch_533 - scratch_314*scratch_534 - scratch_315*scratch_535 - scratch_316*scratch_536 - scratch_317*scratch_537 - scratch_318*scratch_538 - scratch_319*scratch_539 - scratch_320*scratch_540;
    amrex::Real scratch_542 = scratch_364 - scratch_365*scratch_62 - scratch_366*scratch_532 - scratch_367*scratch_533 - scratch_368*scratch_534 - scratch_369*scratch_535 - scratch_370*scratch_536 - scratch_371*scratch_537 - scratch_372*scratch_538 - scratch_373*scratch_539 - scratch_374*scratch_540 - scratch_375*scratch_541;
    amrex::Real scratch_543 = scratch_376*scratch_542;
    amrex::Real scratch_544 = -scratch_79*scratch_98 + scratch_97;
    amrex::Real scratch_545 = scratch_117 - scratch_118*scratch_79 - scratch_119*scratch_544;
    amrex::Real scratch_546 = scratch_139 - scratch_140*scratch_79 - scratch_141*scratch_544 - scratch_142*scratch_545;
    amrex::Real scratch_547 = scratch_163 - scratch_164*scratch_79 - scratch_165*scratch_544 - scratch_166*scratch_545 - scratch_167*scratch_546;
    amrex::Real scratch_548 = scratch_189 - scratch_190*scratch_79 - scratch_191*scratch_544 - scratch_192*scratch_545 - scratch_193*scratch_546 - scratch_194*scratch_547;
    amrex::Real scratch_549 = scratch_217 - scratch_218*scratch_79 - scratch_219*scratch_544 - scratch_220*scratch_545 - scratch_221*scratch_546 - scratch_222*scratch_547 - scratch_223*scratch_548;
    amrex::Real scratch_550 = scratch_247 - scratch_248*scratch_79 - scratch_249*scratch_544 - scratch_250*scratch_545 - scratch_251*scratch_546 - scratch_252*scratch_547 - scratch_253*scratch_548 - scratch_254*scratch_549;
    amrex::Real scratch_551 = scratch_278 - scratch_279*scratch_79 - scratch_280*scratch_544 - scratch_281*scratch_545 - scratch_282*scratch_546 - scratch_283*scratch_547 - scratch_284*scratch_548 - scratch_285*scratch_549 - scratch_286*scratch_550;
    amrex::Real scratch_552 = scratch_311 - scratch_312*scratch_79 - scratch_313*scratch_544 - scratch_314*scratch_545 - scratch_315*scratch_546 - scratch_316*scratch_547 - scratch_317*scratch_548 - scratch_318*scratch_549 - scratch_319*scratch_550 - scratch_320*scratch_551;
    amrex::Real scratch_553 = scratch_365 - scratch_366*scratch_79 - scratch_367*scratch_544 - scratch_368*scratch_545 - scratch_369*scratch_546 - scratch_370*scratch_547 - scratch_371*scratch_548 - scratch_372*scratch_549 - scratch_373*scratch_550 - scratch_374*scratch_551 - scratch_375*scratch_552;
    amrex::Real scratch_554 = scratch_376*scratch_553;
    amrex::Real scratch_555 = scratch_118 - scratch_119*scratch_98;
    amrex::Real scratch_556 = scratch_140 - scratch_141*scratch_98 - scratch_142*scratch_555;
    amrex::Real scratch_557 = scratch_164 - scratch_165*scratch_98 - scratch_166*scratch_555 - scratch_167*scratch_556;
    amrex::Real scratch_558 = scratch_190 - scratch_191*scratch_98 - scratch_192*scratch_555 - scratch_193*scratch_556 - scratch_194*scratch_557;
    amrex::Real scratch_559 = scratch_218 - scratch_219*scratch_98 - scratch_220*scratch_555 - scratch_221*scratch_556 - scratch_222*scratch_557 - scratch_223*scratch_558;
    amrex::Real scratch_560 = scratch_248 - scratch_249*scratch_98 - scratch_250*scratch_555 - scratch_251*scratch_556 - scratch_252*scratch_557 - scratch_253*scratch_558 - scratch_254*scratch_559;
    amrex::Real scratch_561 = scratch_279 - scratch_280*scratch_98 - scratch_281*scratch_555 - scratch_282*scratch_556 - scratch_283*scratch_557 - scratch_284*scratch_558 - scratch_285*scratch_559 - scratch_286*scratch_560;
    amrex::Real scratch_562 = scratch_312 - scratch_313*scratch_98 - scratch_314*scratch_555 - scratch_315*scratch_556 - scratch_316*scratch_557 - scratch_317*scratch_558 - scratch_318*scratch_559 - scratch_319*scratch_560 - scratch_320*scratch_561;
    amrex::Real scratch_563 = scratch_366 - scratch_367*scratch_98 - scratch_368*scratch_555 - scratch_369*scratch_556 - scratch_370*scratch_557 - scratch_371*scratch_558 - scratch_372*scratch_559 - scratch_373*scratch_560 - scratch_374*scratch_561 - scratch_375*scratch_562;
    amrex::Real scratch_564 = scratch_376*scratch_563;
    amrex::Real scratch_565 = -scratch_119*scratch_142 + scratch_141;
    amrex::Real scratch_566 = -scratch_119*scratch_166 + scratch_165 - scratch_167*scratch_565;
    amrex::Real scratch_567 = -scratch_119*scratch_192 + scratch_191 - scratch_193*scratch_565 - scratch_194*scratch_566;
    amrex::Real scratch_568 = -scratch_119*scratch_220 + scratch_219 - scratch_221*scratch_565 - scratch_222*scratch_566 - scratch_223*scratch_567;
    amrex::Real scratch_569 = -scratch_119*scratch_250 + scratch_249 - scratch_251*scratch_565 - scratch_252*scratch_566 - scratch_253*scratch_567 - scratch_254*scratch_568;
    amrex::Real scratch_570 = -scratch_119*scratch_281 + scratch_280 - scratch_282*scratch_565 - scratch_283*scratch_566 - scratch_284*scratch_567 - scratch_285*scratch_568 - scratch_286*scratch_569;
    amrex::Real scratch_571 = -scratch_119*scratch_314 + scratch_313 - scratch_315*scratch_565 - scratch_316*scratch_566 - scratch_317*scratch_567 - scratch_318*scratch_568 - scratch_319*scratch_569 - scratch_320*scratch_570;
    amrex::Real scratch_572 = -scratch_119*scratch_368 + scratch_367 - scratch_369*scratch_565 - scratch_370*scratch_566 - scratch_371*scratch_567 - scratch_372*scratch_568 - scratch_373*scratch_569 - scratch_374*scratch_570 - scratch_375*scratch_571;
    amrex::Real scratch_573 = scratch_376*scratch_572;
    amrex::Real scratch_574 = -scratch_142*scratch_167 + scratch_166;
    amrex::Real scratch_575 = -scratch_142*scratch_193 + scratch_192 - scratch_194*scratch_574;
    amrex::Real scratch_576 = -scratch_142*scratch_221 + scratch_220 - scratch_222*scratch_574 - scratch_223*scratch_575;
    amrex::Real scratch_577 = -scratch_142*scratch_251 + scratch_250 - scratch_252*scratch_574 - scratch_253*scratch_575 - scratch_254*scratch_576;
    amrex::Real scratch_578 = -scratch_142*scratch_282 + scratch_281 - scratch_283*scratch_574 - scratch_284*scratch_575 - scratch_285*scratch_576 - scratch_286*scratch_577;
    amrex::Real scratch_579 = -scratch_142*scratch_315 + scratch_314 - scratch_316*scratch_574 - scratch_317*scratch_575 - scratch_318*scratch_576 - scratch_319*scratch_577 - scratch_320*scratch_578;
    amrex::Real scratch_580 = -scratch_142*scratch_369 + scratch_368 - scratch_370*scratch_574 - scratch_371*scratch_575 - scratch_372*scratch_576 - scratch_373*scratch_577 - scratch_374*scratch_578 - scratch_375*scratch_579;
    amrex::Real scratch_581 = scratch_376*scratch_580;
    amrex::Real scratch_582 = -scratch_167*scratch_194 + scratch_193;
    amrex::Real scratch_583 = -scratch_167*scratch_222 + scratch_221 - scratch_223*scratch_582;
    amrex::Real scratch_584 = -scratch_167*scratch_252 + scratch_251 - scratch_253*scratch_582 - scratch_254*scratch_583;
    amrex::Real scratch_585 = -scratch_167*scratch_283 + scratch_282 - scratch_284*scratch_582 - scratch_285*scratch_583 - scratch_286*scratch_584;
    amrex::Real scratch_586 = -scratch_167*scratch_316 + scratch_315 - scratch_317*scratch_582 - scratch_318*scratch_583 - scratch_319*scratch_584 - scratch_320*scratch_585;
    amrex::Real scratch_587 = -scratch_167*scratch_370 + scratch_369 - scratch_371*scratch_582 - scratch_372*scratch_583 - scratch_373*scratch_584 - scratch_374*scratch_585 - scratch_375*scratch_586;
    amrex::Real scratch_588 = scratch_376*scratch_587;
    amrex::Real scratch_589 = -scratch_194*scratch_223 + scratch_222;
    amrex::Real scratch_590 = -scratch_194*scratch_253 + scratch_252 - scratch_254*scratch_589;
    amrex::Real scratch_591 = -scratch_194*scratch_284 + scratch_283 - scratch_285*scratch_589 - scratch_286*scratch_590;
    amrex::Real scratch_592 = -scratch_194*scratch_317 + scratch_316 - scratch_318*scratch_589 - scratch_319*scratch_590 - scratch_320*scratch_591;
    amrex::Real scratch_593 = -scratch_194*scratch_371 + scratch_370 - scratch_372*scratch_589 - scratch_373*scratch_590 - scratch_374*scratch_591 - scratch_375*scratch_592;
    amrex::Real scratch_594 = scratch_376*scratch_593;
    amrex::Real scratch_595 = -scratch_223*scratch_254 + scratch_253;
    amrex::Real scratch_596 = -scratch_223*scratch_285 + scratch_284 - scratch_286*scratch_595;
    amrex::Real scratch_597 = -scratch_223*scratch_318 + scratch_317 - scratch_319*scratch_595 - scratch_320*scratch_596;
    amrex::Real scratch_598 = -scratch_223*scratch_372 + scratch_371 - scratch_373*scratch_595 - scratch_374*scratch_596 - scratch_375*scratch_597;
    amrex::Real scratch_599 = scratch_376*scratch_598;
    amrex::Real scratch_600 = -scratch_254*scratch_286 + scratch_285;
    amrex::Real scratch_601 = -scratch_254*scratch_319 + scratch_318 - scratch_320*scratch_600;
    amrex::Real scratch_602 = -scratch_254*scratch_373 + scratch_372 - scratch_374*scratch_600 - scratch_375*scratch_601;
    amrex::Real scratch_603 = scratch_376*scratch_602;
    amrex::Real scratch_604 = -scratch_286*scratch_320 + scratch_319;
    amrex::Real scratch_605 = -scratch_286*scratch_374 + scratch_373 - scratch_375*scratch_604;
    amrex::Real scratch_606 = scratch_376*scratch_605;
    amrex::Real scratch_607 = -scratch_320*scratch_375 + scratch_374;
    amrex::Real scratch_608 = scratch_376*scratch_607;
    amrex::Real scratch_609 = scratch_375*scratch_376;

    x[0] = scratch_1 - scratch_120*scratch_409 - scratch_143*scratch_410 - scratch_15*scratch_402 - scratch_168*scratch_411 - scratch_195*scratch_412 - scratch_224*scratch_413 - scratch_24*scratch_403 - scratch_255*scratch_414 - scratch_287*scratch_415 - scratch_321*scratch_416 - scratch_35*scratch_404 - scratch_377*scratch_417 - scratch_379*(-scratch_10*scratch_4 - scratch_109*scratch_99 - scratch_11*scratch_7 - scratch_120*scratch_131 - scratch_143*scratch_155 - scratch_15*scratch_19 - scratch_168*scratch_181 - scratch_195*scratch_209 + scratch_2 - scratch_224*scratch_239 - scratch_24*scratch_29 - scratch_255*scratch_270 - scratch_287*scratch_303 - scratch_321*scratch_339 - scratch_35*scratch_41 - scratch_357*scratch_378 - scratch_48*scratch_55 - scratch_63*scratch_71 - scratch_80*scratch_89) - scratch_399*(-scratch_120*scratch_390 - scratch_143*scratch_391 - scratch_15*scratch_383 - scratch_168*scratch_392 - scratch_195*scratch_393 - scratch_224*scratch_394 - scratch_24*scratch_384 - scratch_255*scratch_395 - scratch_287*scratch_396 - scratch_321*scratch_397 - scratch_35*scratch_385 - scratch_378*scratch_398 + scratch_380 - scratch_381*scratch_4 - scratch_382*scratch_7 - scratch_386*scratch_48 - scratch_387*scratch_63 - scratch_388*scratch_80 - scratch_389*scratch_99) - scratch_4*scratch_400 - scratch_401*scratch_7 - scratch_405*scratch_48 - scratch_406*scratch_63 - scratch_407*scratch_80 - scratch_408*scratch_99 - scratch_456*(-scratch_120*scratch_428 - scratch_143*scratch_429 - scratch_15*scratch_421 - scratch_168*scratch_430 - scratch_195*scratch_431 - scratch_224*scratch_432 - scratch_24*scratch_422 - scratch_255*scratch_433 - scratch_287*scratch_434 - scratch_321*scratch_435 - scratch_35*scratch_423 - scratch_378*scratch_436 - scratch_4*scratch_419 + scratch_418 - scratch_420*scratch_7 - scratch_424*scratch_48 - scratch_425*scratch_63 - scratch_426*scratch_80 - scratch_427*scratch_99);
    x[1] = -scratch_379*(scratch_10 - scratch_109*scratch_463 - scratch_11*scratch_6 - scratch_131*scratch_464 - scratch_155*scratch_465 - scratch_181*scratch_466 - scratch_19*scratch_457 - scratch_209*scratch_467 - scratch_239*scratch_468 - scratch_270*scratch_469 - scratch_29*scratch_458 - scratch_303*scratch_470 - scratch_339*scratch_471 - scratch_357*scratch_473 - scratch_41*scratch_459 - scratch_460*scratch_55 - scratch_461*scratch_71 - scratch_462*scratch_89) - scratch_399*(scratch_381 - scratch_382*scratch_6 - scratch_383*scratch_457 - scratch_384*scratch_458 - scratch_385*scratch_459 - scratch_386*scratch_460 - scratch_387*scratch_461 - scratch_388*scratch_462 - scratch_389*scratch_463 - scratch_390*scratch_464 - scratch_391*scratch_465 - scratch_392*scratch_466 - scratch_393*scratch_467 - scratch_394*scratch_468 - scratch_395*scratch_469 - scratch_396*scratch_470 - scratch_397*scratch_471 - scratch_398*scratch_473) + scratch_400 - scratch_401*scratch_6 - scratch_402*scratch_457 - scratch_403*scratch_458 - scratch_404*scratch_459 - scratch_405*scratch_460 - scratch_406*scratch_461 - scratch_407*scratch_462 - scratch_408*scratch_463 - scratch_409*scratch_464 - scratch_410*scratch_465 - scratch_411*scratch_466 - scratch_412*scratch_467 - scratch_413*scratch_468 - scratch_414*scratch_469 - scratch_415*scratch_470 - scratch_416*scratch_471 - scratch_417*scratch_472 - scratch_456*(scratch_419 - scratch_420*scratch_6 - scratch_421*scratch_457 - scratch_422*scratch_458 - scratch_423*scratch_459 - scratch_424*scratch_460 - scratch_425*scratch_461 - scratch_426*scratch_462 - scratch_427*scratch_463 - scratch_428*scratch_464 - scratch_429*scratch_465 - scratch_430*scratch_466 - scratch_431*scratch_467 - scratch_432*scratch_468 - scratch_433*scratch_469 - scratch_434*scratch_470 - scratch_435*scratch_471 - scratch_436*scratch_473);
    x[2] = -scratch_14*scratch_402 - scratch_379*(-scratch_109*scratch_479 + scratch_11 - scratch_131*scratch_480 - scratch_14*scratch_19 - scratch_155*scratch_481 - scratch_181*scratch_482 - scratch_209*scratch_483 - scratch_239*scratch_484 - scratch_270*scratch_485 - scratch_29*scratch_474 - scratch_303*scratch_486 - scratch_339*scratch_487 - scratch_357*scratch_489 - scratch_41*scratch_475 - scratch_476*scratch_55 - scratch_477*scratch_71 - scratch_478*scratch_89) - scratch_399*(-scratch_14*scratch_383 + scratch_382 - scratch_384*scratch_474 - scratch_385*scratch_475 - scratch_386*scratch_476 - scratch_387*scratch_477 - scratch_388*scratch_478 - scratch_389*scratch_479 - scratch_390*scratch_480 - scratch_391*scratch_481 - scratch_392*scratch_482 - scratch_393*scratch_483 - scratch_394*scratch_484 - scratch_395*scratch_485 - scratch_396*scratch_486 - scratch_397*scratch_487 - scratch_398*scratch_489) + scratch_401 - scratch_403*scratch_474 - scratch_404*scratch_475 - scratch_405*scratch_476 - scratch_406*scratch_477 - scratch_407*scratch_478 - scratch_408*scratch_479 - scratch_409*scratch_480 - scratch_410*scratch_481 - scratch_411*scratch_482 - scratch_412*scratch_483 - scratch_413*scratch_484 - scratch_414*scratch_485 - scratch_415*scratch_486 - scratch_416*scratch_487 - scratch_417*scratch_488 - scratch_456*(-scratch_14*scratch_421 + scratch_420 - scratch_422*scratch_474 - scratch_423*scratch_475 - scratch_424*scratch_476 - scratch_425*scratch_477 - scratch_426*scratch_478 - scratch_427*scratch_479 - scratch_428*scratch_480 - scratch_429*scratch_481 - scratch_430*scratch_482 - scratch_431*scratch_483 - scratch_432*scratch_484 - scratch_433*scratch_485 - scratch_434*scratch_486 - scratch_435*scratch_487 - scratch_436*scratch_489);
    x[3] = -scratch_23*scratch_403 - scratch_379*(-scratch_109*scratch_494 - scratch_131*scratch_495 - scratch_155*scratch_496 - scratch_181*scratch_497 + scratch_19 - scratch_209*scratch_498 - scratch_23*scratch_29 - scratch_239*scratch_499 - scratch_270*scratch_500 - scratch_303*scratch_501 - scratch_339*scratch_502 - scratch_357*scratch_504 - scratch_41*scratch_490 - scratch_491*scratch_55 - scratch_492*scratch_71 - scratch_493*scratch_89) - scratch_399*(-scratch_23*scratch_384 + scratch_383 - scratch_385*scratch_490 - scratch_386*scratch_491 - scratch_387*scratch_492 - scratch_388*scratch_493 - scratch_389*scratch_494 - scratch_390*scratch_495 - scratch_391*scratch_496 - scratch_392*scratch_497 - scratch_393*scratch_498 - scratch_394*scratch_499 - scratch_395*scratch_500 - scratch_396*scratch_501 - scratch_397*scratch_502 - scratch_398*scratch_504) + scratch_402 - scratch_404*scratch_490 - scratch_405*scratch_491 - scratch_406*scratch_492 - scratch_407*scratch_493 - scratch_408*scratch_494 - scratch_409*scratch_495 - scratch_410*scratch_496 - scratch_411*scratch_497 - scratch_412*scratch_498 - scratch_413*scratch_499 - scratch_414*scratch_500 - scratch_415*scratch_501 - scratch_416*scratch_502 - scratch_417*scratch_503 - scratch_456*(-scratch_23*scratch_422 + scratch_421 - scratch_423*scratch_490 - scratch_424*scratch_491 - scratch_425*scratch_492 - scratch_426*scratch_493 - scratch_427*scratch_494 - scratch_428*scratch_495 - scratch_429*scratch_496 - scratch_430*scratch_497 - scratch_431*scratch_498 - scratch_432*scratch_499 - scratch_433*scratch_500 - scratch_434*scratch_501 - scratch_435*scratch_502 - scratch_436*scratch_504);
    x[4] = -scratch_34*scratch_404 - scratch_379*(-scratch_109*scratch_508 - scratch_131*scratch_509 - scratch_155*scratch_510 - scratch_181*scratch_511 - scratch_209*scratch_512 - scratch_239*scratch_513 - scratch_270*scratch_514 + scratch_29 - scratch_303*scratch_515 - scratch_339*scratch_516 - scratch_34*scratch_41 - scratch_357*scratch_518 - scratch_505*scratch_55 - scratch_506*scratch_71 - scratch_507*scratch_89) - scratch_399*(-scratch_34*scratch_385 + scratch_384 - scratch_386*scratch_505 - scratch_387*scratch_506 - scratch_388*scratch_507 - scratch_389*scratch_508 - scratch_390*scratch_509 - scratch_391*scratch_510 - scratch_392*scratch_511 - scratch_393*scratch_512 - scratch_394*scratch_513 - scratch_395*scratch_514 - scratch_396*scratch_515 - scratch_397*scratch_516 - scratch_398*scratch_518) + scratch_403 - scratch_405*scratch_505 - scratch_406*scratch_506 - scratch_407*scratch_507 - scratch_408*scratch_508 - scratch_409*scratch_509 - scratch_410*scratch_510 - scratch_411*scratch_511 - scratch_412*scratch_512 - scratch_413*scratch_513 - scratch_414*scratch_514 - scratch_415*scratch_515 - scratch_416*scratch_516 - scratch_417*scratch_517 - scratch_456*(-scratch_34*scratch_423 + scratch_422 - scratch_424*scratch_505 - scratch_425*scratch_506 - scratch_426*scratch_507 - scratch_427*scratch_508 - scratch_428*scratch_509 - scratch_429*scratch_510 - scratch_430*scratch_511 - scratch_431*scratch_512 - scratch_432*scratch_513 - scratch_433*scratch_514 - scratch_434*scratch_515 - scratch_435*scratch_516 - scratch_436*scratch_518);
    x[5] = -scratch_379*(-scratch_109*scratch_521 - scratch_131*scratch_522 - scratch_155*scratch_523 - scratch_181*scratch_524 - scratch_209*scratch_525 - scratch_239*scratch_526 - scratch_270*scratch_527 - scratch_303*scratch_528 - scratch_339*scratch_529 - scratch_357*scratch_531 + scratch_41 - scratch_47*scratch_55 - scratch_519*scratch_71 - scratch_520*scratch_89) - scratch_399*(scratch_385 - scratch_386*scratch_47 - scratch_387*scratch_519 - scratch_388*scratch_520 - scratch_389*scratch_521 - scratch_390*scratch_522 - scratch_391*scratch_523 - scratch_392*scratch_524 - scratch_393*scratch_525 - scratch_394*scratch_526 - scratch_395*scratch_527 - scratch_396*scratch_528 - scratch_397*scratch_529 - scratch_398*scratch_531) + scratch_404 - scratch_405*scratch_47 - scratch_406*scratch_519 - scratch_407*scratch_520 - scratch_408*scratch_521 - scratch_409*scratch_522 - scratch_410*scratch_523 - scratch_411*scratch_524 - scratch_412*scratch_525 - scratch_413*scratch_526 - scratch_414*scratch_527 - scratch_415*scratch_528 - scratch_416*scratch_529 - scratch_417*scratch_530 - scratch_456*(scratch_423 - scratch_424*scratch_47 - scratch_425*scratch_519 - scratch_426*scratch_520 - scratch_427*scratch_521 - scratch_428*scratch_522 - scratch_429*scratch_523 - scratch_430*scratch_524 - scratch_431*scratch_525 - scratch_432*scratch_526 - scratch_433*scratch_527 - scratch_434*scratch_528 - scratch_435*scratch_529 - scratch_436*scratch_531);
    x[6] = -scratch_379*(-scratch_109*scratch_533 - scratch_131*scratch_534 - scratch_155*scratch_535 - scratch_181*scratch_536 - scratch_209*scratch_537 - scratch_239*scratch_538 - scratch_270*scratch_539 - scratch_303*scratch_540 - scratch_339*scratch_541 - scratch_357*scratch_543 - scratch_532*scratch_89 + scratch_55 - scratch_62*scratch_71) - scratch_399*(scratch_386 - scratch_387*scratch_62 - scratch_388*scratch_532 - scratch_389*scratch_533 - scratch_390*scratch_534 - scratch_391*scratch_535 - scratch_392*scratch_536 - scratch_393*scratch_537 - scratch_394*scratch_538 - scratch_395*scratch_539 - scratch_396*scratch_540 - scratch_397*scratch_541 - scratch_398*scratch_543) + scratch_405 - scratch_406*scratch_62 - scratch_407*scratch_532 - scratch_408*scratch_533 - scratch_409*scratch_534 - scratch_410*scratch_535 - scratch_411*scratch_536 - scratch_412*scratch_537 - scratch_413*scratch_538 - scratch_414*scratch_539 - scratch_415*scratch_540 - scratch_416*scratch_541 - scratch_417*scratch_542 - scratch_456*(scratch_424 - scratch_425*scratch_62 - scratch_426*scratch_532 - scratch_427*scratch_533 - scratch_428*scratch_534 - scratch_429*scratch_535 - scratch_430*scratch_536 - scratch_431*scratch_537 - scratch_432*scratch_538 - scratch_433*scratch_539 - scratch_434*scratch_540 - scratch_435*scratch_541 - scratch_436*scratch_543);
    x[7] = -scratch_379*(-scratch_109*scratch_544 - scratch_131*scratch_545 - scratch_155*scratch_546 - scratch_181*scratch_547 - scratch_209*scratch_548 - scratch_239*scratch_549 - scratch_270*scratch_550 - scratch_303*scratch_551 - scratch_339*scratch_552 - scratch_357*scratch_554 + scratch_71 - scratch_79*scratch_89) - scratch_399*(scratch_387 - scratch_388*scratch_79 - scratch_389*scratch_544 - scratch_390*scratch_545 - scratch_391*scratch_546 - scratch_392*scratch_547 - scratch_393*scratch_548 - scratch_394*scratch_549 - scratch_395*scratch_550 - scratch_396*scratch_551 - scratch_397*scratch_552 - scratch_398*scratch_554) + scratch_406 - scratch_407*scratch_79 - scratch_408*scratch_544 - scratch_409*scratch_545 - scratch_410*scratch_546 - scratch_411*scratch_547 - scratch_412*scratch_548 - scratch_413*scratch_549 - scratch_414*scratch_550 - scratch_415*scratch_551 - scratch_416*scratch_552 - scratch_417*scratch_553 - scratch_456*(scratch_425 - scratch_426*scratch_79 - scratch_427*scratch_544 - scratch_428*scratch_545 - scratch_429*scratch_546 - scratch_430*scratch_547 - scratch_431*scratch_548 - scratch_432*scratch_549 - scratch_433*scratch_550 - scratch_434*scratch_551 - scratch_435*scratch_552 - scratch_436*scratch_554);
    x[8] = -scratch_379*(-scratch_109*scratch_98 - scratch_131*scratch_555 - scratch_155*scratch_556 - scratch_181*scratch_557 - scratch_209*scratch_558 - scratch_239*scratch_559 - scratch_270*scratch_560 - scratch_303*scratch_561 - scratch_339*scratch_562 - scratch_357*scratch_564 + scratch_89) - scratch_399*(scratch_388 - scratch_389*scratch_98 - scratch_390*scratch_555 - scratch_391*scratch_556 - scratch_392*scratch_557 - scratch_393*scratch_558 - scratch_394*scratch_559 - scratch_395*scratch_560 - scratch_396*scratch_561 - scratch_397*scratch_562 - scratch_398*scratch_564) + scratch_407 - scratch_408*scratch_98 - scratch_409*scratch_555 - scratch_410*scratch_556 - scratch_411*scratch_557 - scratch_412*scratch_558 - scratch_413*scratch_559 - scratch_414*scratch_560 - scratch_415*scratch_561 - scratch_416*scratch_562 - scratch_417*scratch_563 - scratch_456*(scratch_426 - scratch_427*scratch_98 - scratch_428*scratch_555 - scratch_429*scratch_556 - scratch_430*scratch_557 - scratch_431*scratch_558 - scratch_432*scratch_559 - scratch_433*scratch_560 - scratch_434*scratch_561 - scratch_435*scratch_562 - scratch_436*scratch_564);
    x[9] = -scratch_119*scratch_409 - scratch_379*(scratch_109 - scratch_119*scratch_131 - scratch_155*scratch_565 - scratch_181*scratch_566 - scratch_209*scratch_567 - scratch_239*scratch_568 - scratch_270*scratch_569 - scratch_303*scratch_570 - scratch_339*scratch_571 - scratch_357*scratch_573) - scratch_399*(-scratch_119*scratch_390 + scratch_389 - scratch_391*scratch_565 - scratch_392*scratch_566 - scratch_393*scratch_567 - scratch_394*scratch_568 - scratch_395*scratch_569 - scratch_396*scratch_570 - scratch_397*scratch_571 - scratch_398*scratch_573) + scratch_408 - scratch_410*scratch_565 - scratch_411*scratch_566 - scratch_412*scratch_567 - scratch_413*scratch_568 - scratch_414*scratch_569 - scratch_415*scratch_570 - scratch_416*scratch_571 - scratch_417*scratch_572 - scratch_456*(-scratch_119*scratch_428 + scratch_427 - scratch_429*scratch_565 - scratch_430*scratch_566 - scratch_431*scratch_567 - scratch_432*scratch_568 - scratch_433*scratch_569 - scratch_434*scratch_570 - scratch_435*scratch_571 - scratch_436*scratch_573);
    x[10] = -scratch_142*scratch_410 - scratch_379*(scratch_131 - scratch_142*scratch_155 - scratch_181*scratch_574 - scratch_209*scratch_575 - scratch_239*scratch_576 - scratch_270*scratch_577 - scratch_303*scratch_578 - scratch_339*scratch_579 - scratch_357*scratch_581) - scratch_399*(-scratch_142*scratch_391 + scratch_390 - scratch_392*scratch_574 - scratch_393*scratch_575 - scratch_394*scratch_576 - scratch_395*scratch_577 - scratch_396*scratch_578 - scratch_397*scratch_579 - scratch_398*scratch_581) + scratch_409 - scratch_411*scratch_574 - scratch_412*scratch_575 - scratch_413*scratch_576 - scratch_414*scratch_577 - scratch_415*scratch_578 - scratch_416*scratch_579 - scratch_417*scratch_580 - scratch_456*(-scratch_142*scratch_429 + scratch_428 - scratch_430*scratch_574 - scratch_431*scratch_575 - scratch_432*scratch_576 - scratch_433*scratch_577 - scratch_434*scratch_578 - scratch_435*scratch_579 - scratch_436*scratch_581);
    x[11] = -scratch_167*scratch_411 - scratch_379*(scratch_155 - scratch_167*scratch_181 - scratch_209*scratch_582 - scratch_239*scratch_583 - scratch_270*scratch_584 - scratch_303*scratch_585 - scratch_339*scratch_586 - scratch_357*scratch_588) - scratch_399*(-scratch_167*scratch_392 + scratch_391 - scratch_393*scratch_582 - scratch_394*scratch_583 - scratch_395*scratch_584 - scratch_396*scratch_585 - scratch_397*scratch_586 - scratch_398*scratch_588) + scratch_410 - scratch_412*scratch_582 - scratch_413*scratch_583 - scratch_414*scratch_584 - scratch_415*scratch_585 - scratch_416*scratch_586 - scratch_417*scratch_587 - scratch_456*(-scratch_167*scratch_430 + scratch_429 - scratch_431*scratch_582 - scratch_432*scratch_583 - scratch_433*scratch_584 - scratch_434*scratch_585 - scratch_435*scratch_586 - scratch_436*scratch_588);
    x[12] = -scratch_194*scratch_412 - scratch_379*(scratch_181 - scratch_194*scratch_209 - scratch_239*scratch_589 - scratch_270*scratch_590 - scratch_303*scratch_591 - scratch_339*scratch_592 - scratch_357*scratch_594) - scratch_399*(-scratch_194*scratch_393 + scratch_392 - scratch_394*scratch_589 - scratch_395*scratch_590 - scratch_396*scratch_591 - scratch_397*scratch_592 - scratch_398*scratch_594) + scratch_411 - scratch_413*scratch_589 - scratch_414*scratch_590 - scratch_415*scratch_591 - scratch_416*scratch_592 - scratch_417*scratch_593 - scratch_456*(-scratch_194*scratch_431 + scratch_430 - scratch_432*scratch_589 - scratch_433*scratch_590 - scratch_434*scratch_591 - scratch_435*scratch_592 - scratch_436*scratch_594);
    x[13] = -scratch_223*scratch_413 - scratch_379*(scratch_209 - scratch_223*scratch_239 - scratch_270*scratch_595 - scratch_303*scratch_596 - scratch_339*scratch_597 - scratch_357*scratch_599) - scratch_399*(-scratch_223*scratch_394 + scratch_393 - scratch_395*scratch_595 - scratch_396*scratch_596 - scratch_397*scratch_597 - scratch_398*scratch_599) + scratch_412 - scratch_414*scratch_595 - scratch_415*scratch_596 - scratch_416*scratch_597 - scratch_417*scratch_598 - scratch_456*(-scratch_223*scratch_432 + scratch_431 - scratch_433*scratch_595 - scratch_434*scratch_596 - scratch_435*scratch_597 - scratch_436*scratch_599);
    x[14] = -scratch_254*scratch_414 - scratch_379*(scratch_239 - scratch_254*scratch_270 - scratch_303*scratch_600 - scratch_339*scratch_601 - scratch_357*scratch_603) - scratch_399*(-scratch_254*scratch_395 + scratch_394 - scratch_396*scratch_600 - scratch_397*scratch_601 - scratch_398*scratch_603) + scratch_413 - scratch_415*scratch_600 - scratch_416*scratch_601 - scratch_417*scratch_602 - scratch_456*(-scratch_254*scratch_433 + scratch_432 - scratch_434*scratch_600 - scratch_435*scratch_601 - scratch_436*scratch_603);
    x[15] = -scratch_286*scratch_415 - scratch_379*(scratch_270 - scratch_286*scratch_303 - scratch_339*scratch_604 - scratch_357*scratch_606) - scratch_399*(-scratch_286*scratch_396 + scratch_395 - scratch_397*scratch_604 - scratch_398*scratch_606) + scratch_414 - scratch_416*scratch_604 - scratch_417*scratch_605 - scratch_456*(-scratch_286*scratch_434 + scratch_433 - scratch_435*scratch_604 - scratch_436*scratch_606);
    x[16] = -scratch_320*scratch_416 - scratch_379*(scratch_303 - scratch_320*scratch_339 - scratch_357*scratch_608) - scratch_399*(-scratch_320*scratch_397 + scratch_396 - scratch_398*scratch_608) + scratch_415 - scratch_417*scratch_607 - scratch_456*(-scratch_320*scratch_435 + scratch_434 - scratch_436*scratch_608);
    x[17] = -scratch_375*scratch_417 - scratch_379*(scratch_339 - scratch_357*scratch_609) - scratch_399*(scratch_397 - scratch_398*scratch_609) + scratch_416 - scratch_456*(scratch_435 - scratch_436*scratch_609);
    x[18] = -scratch_357*scratch_376*scratch_379 - scratch_376*scratch_398*scratch_399 - scratch_376*scratch_436*scratch_456 + scratch_417;
    x[19] = scratch_379;
    x[20] = scratch_399;
    x[21] = scratch_456;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void sgjsolve(amrex::Real* A, amrex::Real* x, amrex::Real* b) {
    amrex::Real scratch_0 = 1.0/A[0];
    amrex::Real scratch_1 = b[0]*scratch_0;
    amrex::Real scratch_2 = A[19]*scratch_0;
    amrex::Real scratch_3 = A[2]*scratch_0;
    amrex::Real scratch_4 = A[1]*scratch_0;
    amrex::Real scratch_5 = 1.0/(-A[22]*scratch_4 + A[23]);
    amrex::Real scratch_6 = scratch_5*(-A[22]*scratch_3 + A[24]);
    amrex::Real scratch_7 = scratch_3 - scratch_4*scratch_6;
    amrex::Real scratch_8 = -A[44]*scratch_4 + A[45];
    amrex::Real scratch_9 = 1.0/(-A[44]*scratch_3 + A[46] - scratch_6*scratch_8);
    amrex::Real scratch_10 = scratch_5*(-A[22]*scratch_2 + A[41]);
    amrex::Real scratch_11 = scratch_9*(-A[44]*scratch_2 + A[63] - scratch_10*scratch_8);
    amrex::Real scratch_12 = A[3]*scratch_0;
    amrex::Real scratch_13 = scratch_5*(-A[22]*scratch_12 + A[25]);
    amrex::Real scratch_14 = scratch_9*(-A[44]*scratch_12 + A[47] - scratch_13*scratch_8);
    amrex::Real scratch_15 = scratch_12 - scratch_13*scratch_4 - scratch_14*scratch_7;
    amrex::Real scratch_16 = -A[66]*scratch_4 + A[67];
    amrex::Real scratch_17 = -A[66]*scratch_3 + A[68] - scratch_16*scratch_6;
    amrex::Real scratch_18 = 1.0/(-A[66]*scratch_12 + A[69] - scratch_13*scratch_16 - scratch_14*scratch_17);
    amrex::Real scratch_19 = scratch_18*(-A[66]*scratch_2 + A[85] - scratch_10*scratch_16 - scratch_11*scratch_17);
    amrex::Real scratch_20 = A[4]*scratch_0;
    amrex::Real scratch_21 = scratch_5*(-A[22]*scratch_20 + A[26]);
    amrex::Real scratch_22 = scratch_9*(-A[44]*scratch_20 + A[48] - scratch_21*scratch_8);
    amrex::Real scratch_23 = scratch_18*(-A[66]*scratch_20 + A[70] - scratch_16*scratch_21 - scratch_17*scratch_22);
    amrex::Real scratch_24 = -scratch_15*scratch_23 + scratch_20 - scratch_21*scratch_4 - scratch_22*scratch_7;
    amrex::Real scratch_25 = -A[88]*scratch_4 + A[89];
    amrex::Real scratch_26 = -A[88]*scratch_3 + A[90] - scratch_25*scratch_6;
    amrex::Real scratch_27 = -A[88]*scratch_12 + A[91] - scratch_13*scratch_25 - scratch_14*scratch_26;
    amrex::Real scratch_28 = 1.0/(-A[88]*scratch_20 + A[92] - scratch_21*scratch_25 - scratch_22*scratch_26 - scratch_23*scratch_27);
    amrex::Real scratch_29 = scratch_28*(A[107] - A[88]*scratch_2 - scratch_10*scratch_25 - scratch_11*scratch_26 - scratch_19*scratch_27);
    amrex::Real scratch_30 = A[5]*scratch_0;
    amrex::Real scratch_31 = scratch_5*(-A[22]*scratch_30 + A[27]);
    amrex::Real scratch_32 = scratch_9*(-A[44]*scratch_30 + A[49] - scratch_31*scratch_8);
    amrex::Real scratch_33 = scratch_18*(-A[66]*scratch_30 + A[71] - scratch_16*scratch_31 - scratch_17*scratch_32);
    amrex::Real scratch_34 = scratch_28*(-A[88]*scratch_30 + A[93] - scratch_25*scratch_31 - scratch_26*scratch_32 - scratch_27*scratch_33);
    amrex::Real scratch_35 = -scratch_15*scratch_33 - scratch_24*scratch_34 + scratch_30 - scratch_31*scratch_4 - scratch_32*scratch_7;
    amrex::Real scratch_36 = -A[110]*scratch_4 + A[111];
    amrex::Real scratch_37 = -A[110]*scratch_3 + A[112] - scratch_36*scratch_6;
    amrex::Real scratch_38 = -A[110]*scratch_12 + A[113] - scratch_13*scratch_36 - scratch_14*scratch_37;
    amrex::Real scratch_39 = -A[110]*scratch_20 + A[114] - scratch_21*scratch_36 - scratch_22*scratch_37 - scratch_23*scratch_38;
    amrex::Real scratch_40 = 1.0/(-A[110]*scratch_30 + A[115] - scratch_31*scratch_36 - scratch_32*scratch_37 - scratch_33*scratch_38 - scratch_34*scratch_39);
    amrex::Real scratch_41 = scratch_40*(-A[110]*scratch_2 + A[129] - scratch_10*scratch_36 - scratch_11*scratch_37 - scratch_19*scratch_38 - scratch_29*scratch_39);
    amrex::Real scratch_42 = A[6]*scratch_0;
    amrex::Real scratch_43 = scratch_5*(-A[22]*scratch_42 + A[28]);
    amrex::Real scratch_44 = scratch_9*(-A[44]*scratch_42 + A[50] - scratch_43*scratch_8);
    amrex::Real scratch_45 = scratch_18*(-A[66]*scratch_42 + A[72] - scratch_16*scratch_43 - scratch_17*scratch_44);
    amrex::Real scratch_46 = scratch_28*(-A[88]*scratch_42 + A[94] - scratch_25*scratch_43 - scratch_26*scratch_44 - scratch_27*scratch_45);
    amrex::Real scratch_47 = scratch_40*(-A[110]*scratch_42 + A[116] - scratch_36*scratch_43 - scratch_37*scratch_44 - scratch_38*scratch_45 - scratch_39*scratch_46);
    amrex::Real scratch_48 = -scratch_15*scratch_45 - scratch_24*scratch_46 - scratch_35*scratch_47 - scratch_4*scratch_43 + scratch_42 - scratch_44*scratch_7;
    amrex::Real scratch_49 = -A[132]*scratch_4 + A[133];
    amrex::Real scratch_50 = -A[132]*scratch_3 + A[134] - scratch_49*scratch_6;
    amrex::Real scratch_51 = -A[132]*scratch_12 + A[135] - scratch_13*scratch_49 - scratch_14*scratch_50;
    amrex::Real scratch_52 = -A[132]*scratch_20 + A[136] - scratch_21*scratch_49 - scratch_22*scratch_50 - scratch_23*scratch_51;
    amrex::Real scratch_53 = -A[132]*scratch_30 + A[137] - scratch_31*scratch_49 - scratch_32*scratch_50 - scratch_33*scratch_51 - scratch_34*scratch_52;
    amrex::Real scratch_54 = 1.0/(-A[132]*scratch_42 + A[138] - scratch_43*scratch_49 - scratch_44*scratch_50 - scratch_45*scratch_51 - scratch_46*scratch_52 - scratch_47*scratch_53);
    amrex::Real scratch_55 = scratch_54*(-A[132]*scratch_2 + A[151] - scratch_10*scratch_49 - scratch_11*scratch_50 - scratch_19*scratch_51 - scratch_29*scratch_52 - scratch_41*scratch_53);
    amrex::Real scratch_56 = A[7]*scratch_0;
    amrex::Real scratch_57 = scratch_5*(-A[22]*scratch_56 + A[29]);
    amrex::Real scratch_58 = scratch_9*(-A[44]*scratch_56 + A[51] - scratch_57*scratch_8);
    amrex::Real scratch_59 = scratch_18*(-A[66]*scratch_56 + A[73] - scratch_16*scratch_57 - scratch_17*scratch_58);
    amrex::Real scratch_60 = scratch_28*(-A[88]*scratch_56 + A[95] - scratch_25*scratch_57 - scratch_26*scratch_58 - scratch_27*scratch_59);
    amrex::Real scratch_61 = scratch_40*(-A[110]*scratch_56 + A[117] - scratch_36*scratch_57 - scratch_37*scratch_58 - scratch_38*scratch_59 - scratch_39*scratch_60);
    amrex::Real scratch_62 = scratch_54*(-A[132]*scratch_56 + A[139] - scratch_49*scratch_57 - scratch_50*scratch_58 - scratch_51*scratch_59 - scratch_52*scratch_60 - scratch_53*scratch_61);
    amrex::Real scratch_63 = -scratch_15*scratch_59 - scratch_24*scratch_60 - scratch_35*scratch_61 - scratch_4*scratch_57 - scratch_48*scratch_62 + scratch_56 - scratch_58*scratch_7;
    amrex::Real scratch_64 = -A[154]*scratch_4 + A[155];
    amrex::Real scratch_65 = -A[154]*scratch_3 + A[156] - scratch_6*scratch_64;
    amrex::Real scratch_66 = -A[154]*scratch_12 + A[157] - scratch_13*scratch_64 - scratch_14*scratch_65;
    amrex::Real scratch_67 = -A[154]*scratch_20 + A[158] - scratch_21*scratch_64 - scratch_22*scratch_65 - scratch_23*scratch_66;
    amrex::Real scratch_68 = -A[154]*scratch_30 + A[159] - scratch_31*scratch_64 - scratch_32*scratch_65 - scratch_33*scratch_66 - scratch_34*scratch_67;
    amrex::Real scratch_69 = -A[154]*scratch_42 + A[160] - scratch_43*scratch_64 - scratch_44*scratch_65 - scratch_45*scratch_66 - scratch_46*scratch_67 - scratch_47*scratch_68;
    amrex::Real scratch_70 = 1.0/(-A[154]*scratch_56 + A[161] - scratch_57*scratch_64 - scratch_58*scratch_65 - scratch_59*scratch_66 - scratch_60*scratch_67 - scratch_61*scratch_68 - scratch_62*scratch_69);
    amrex::Real scratch_71 = scratch_70*(-A[154]*scratch_2 + A[173] - scratch_10*scratch_64 - scratch_11*scratch_65 - scratch_19*scratch_66 - scratch_29*scratch_67 - scratch_41*scratch_68 - scratch_55*scratch_69);
    amrex::Real scratch_72 = A[8]*scratch_0;
    amrex::Real scratch_73 = scratch_5*(-A[22]*scratch_72 + A[30]);
    amrex::Real scratch_74 = scratch_9*(-A[44]*scratch_72 + A[52] - scratch_73*scratch_8);
    amrex::Real scratch_75 = scratch_18*(-A[66]*scratch_72 + A[74] - scratch_16*scratch_73 - scratch_17*scratch_74);
    amrex::Real scratch_76 = scratch_28*(-A[88]*scratch_72 + A[96] - scratch_25*scratch_73 - scratch_26*scratch_74 - scratch_27*scratch_75);
    amrex::Real scratch_77 = scratch_40*(-A[110]*scratch_72 + A[118] - scratch_36*scratch_73 - scratch_37*scratch_74 - scratch_38*scratch_75 - scratch_39*scratch_76);
    amrex::Real scratch_78 = scratch_54*(-A[132]*scratch_72 + A[140] - scratch_49*scratch_73 - scratch_50*scratch_74 - scratch_51*scratch_75 - scratch_52*scratch_76 - scratch_53*scratch_77);
    amrex::Real scratch_79 = scratch_70*(-A[154]*scratch_72 + A[162] - scratch_64*scratch_73 - scratch_65*scratch_74 - scratch_66*scratch_75 - scratch_67*scratch_76 - scratch_68*scratch_77 - scratch_69*scratch_78);
    amrex::Real scratch_80 = -scratch_15*scratch_75 - scratch_24*scratch_76 - scratch_35*scratch_77 - scratch_4*scratch_73 - scratch_48*scratch_78 - scratch_63*scratch_79 - scratch_7*scratch_74 + scratch_72;
    amrex::Real scratch_81 = -A[176]*scratch_4 + A[177];
    amrex::Real scratch_82 = -A[176]*scratch_3 + A[178] - scratch_6*scratch_81;
    amrex::Real scratch_83 = -A[176]*scratch_12 + A[179] - scratch_13*scratch_81 - scratch_14*scratch_82;
    amrex::Real scratch_84 = -A[176]*scratch_20 + A[180] - scratch_21*scratch_81 - scratch_22*scratch_82 - scratch_23*scratch_83;
    amrex::Real scratch_85 = -A[176]*scratch_30 + A[181] - scratch_31*scratch_81 - scratch_32*scratch_82 - scratch_33*scratch_83 - scratch_34*scratch_84;
    amrex::Real scratch_86 = -A[176]*scratch_42 - scratch_43*scratch_81 - scratch_44*scratch_82 - scratch_45*scratch_83 - scratch_46*scratch_84 - scratch_47*scratch_85;
    amrex::Real scratch_87 = -A[176]*scratch_56 + A[182] - scratch_57*scratch_81 - scratch_58*scratch_82 - scratch_59*scratch_83 - scratch_60*scratch_84 - scratch_61*scratch_85 - scratch_62*scratch_86;
    amrex::Real scratch_88 = 1.0/(-A[176]*scratch_72 + A[183] - scratch_73*scratch_81 - scratch_74*scratch_82 - scratch_75*scratch_83 - scratch_76*scratch_84 - scratch_77*scratch_85 - scratch_78*scratch_86 - scratch_79*scratch_87);
    amrex::Real scratch_89 = scratch_88*(-A[176]*scratch_2 + A[191] - scratch_10*scratch_81 - scratch_11*scratch_82 - scratch_19*scratch_83 - scratch_29*scratch_84 - scratch_41*scratch_85 - scratch_55*scratch_86 - scratch_71*scratch_87);
    amrex::Real scratch_90 = A[9]*scratch_0;
    amrex::Real scratch_91 = scratch_5*(-A[22]*scratch_90 + A[31]);
    amrex::Real scratch_92 = scratch_9*(-A[44]*scratch_90 + A[53] - scratch_8*scratch_91);
    amrex::Real scratch_93 = scratch_18*(-A[66]*scratch_90 + A[75] - scratch_16*scratch_91 - scratch_17*scratch_92);
    amrex::Real scratch_94 = scratch_28*(-A[88]*scratch_90 + A[97] - scratch_25*scratch_91 - scratch_26*scratch_92 - scratch_27*scratch_93);
    amrex::Real scratch_95 = scratch_40*(-A[110]*scratch_90 + A[119] - scratch_36*scratch_91 - scratch_37*scratch_92 - scratch_38*scratch_93 - scratch_39*scratch_94);
    amrex::Real scratch_96 = scratch_54*(-A[132]*scratch_90 + A[141] - scratch_49*scratch_91 - scratch_50*scratch_92 - scratch_51*scratch_93 - scratch_52*scratch_94 - scratch_53*scratch_95);
    amrex::Real scratch_97 = scratch_70*(-A[154]*scratch_90 + A[163] - scratch_64*scratch_91 - scratch_65*scratch_92 - scratch_66*scratch_93 - scratch_67*scratch_94 - scratch_68*scratch_95 - scratch_69*scratch_96);
    amrex::Real scratch_98 = scratch_88*(-A[176]*scratch_90 + A[184] - scratch_81*scratch_91 - scratch_82*scratch_92 - scratch_83*scratch_93 - scratch_84*scratch_94 - scratch_85*scratch_95 - scratch_86*scratch_96 - scratch_87*scratch_97);
    amrex::Real scratch_99 = -scratch_15*scratch_93 - scratch_24*scratch_94 - scratch_35*scratch_95 - scratch_4*scratch_91 - scratch_48*scratch_96 - scratch_63*scratch_97 - scratch_7*scratch_92 - scratch_80*scratch_98 + scratch_90;
    amrex::Real scratch_100 = -A[194]*scratch_4 + A[195];
    amrex::Real scratch_101 = -A[194]*scratch_3 + A[196] - scratch_100*scratch_6;
    amrex::Real scratch_102 = -A[194]*scratch_12 + A[197] - scratch_100*scratch_13 - scratch_101*scratch_14;
    amrex::Real scratch_103 = -A[194]*scratch_20 + A[198] - scratch_100*scratch_21 - scratch_101*scratch_22 - scratch_102*scratch_23;
    amrex::Real scratch_104 = -A[194]*scratch_30 + A[199] - scratch_100*scratch_31 - scratch_101*scratch_32 - scratch_102*scratch_33 - scratch_103*scratch_34;
    amrex::Real scratch_105 = -A[194]*scratch_42 + A[200] - scratch_100*scratch_43 - scratch_101*scratch_44 - scratch_102*scratch_45 - scratch_103*scratch_46 - scratch_104*scratch_47;
    amrex::Real scratch_106 = -A[194]*scratch_56 + A[201] - scratch_100*scratch_57 - scratch_101*scratch_58 - scratch_102*scratch_59 - scratch_103*scratch_60 - scratch_104*scratch_61 - scratch_105*scratch_62;
    amrex::Real scratch_107 = -A[194]*scratch_72 + A[202] - scratch_100*scratch_73 - scratch_101*scratch_74 - scratch_102*scratch_75 - scratch_103*scratch_76 - scratch_104*scratch_77 - scratch_105*scratch_78 - scratch_106*scratch_79;
    amrex::Real scratch_108 = 1.0/(-A[194]*scratch_90 + A[203] - scratch_100*scratch_91 - scratch_101*scratch_92 - scratch_102*scratch_93 - scratch_103*scratch_94 - scratch_104*scratch_95 - scratch_105*scratch_96 - scratch_106*scratch_97 - scratch_107*scratch_98);
    amrex::Real scratch_109 = scratch_108*(-A[194]*scratch_2 + A[213] - scratch_10*scratch_100 - scratch_101*scratch_11 - scratch_102*scratch_19 - scratch_103*scratch_29 - scratch_104*scratch_41 - scratch_105*scratch_55 - scratch_106*scratch_71 - scratch_107*scratch_89);
    amrex::Real scratch_110 = A[10]*scratch_0;
    amrex::Real scratch_111 = scratch_5*(-A[22]*scratch_110 + A[32]);
    amrex::Real scratch_112 = scratch_9*(-A[44]*scratch_110 + A[54] - scratch_111*scratch_8);
    amrex::Real scratch_113 = scratch_18*(-A[66]*scratch_110 + A[76] - scratch_111*scratch_16 - scratch_112*scratch_17);
    amrex::Real scratch_114 = scratch_28*(-A[88]*scratch_110 + A[98] - scratch_111*scratch_25 - scratch_112*scratch_26 - scratch_113*scratch_27);
    amrex::Real scratch_115 = scratch_40*(-A[110]*scratch_110 + A[120] - scratch_111*scratch_36 - scratch_112*scratch_37 - scratch_113*scratch_38 - scratch_114*scratch_39);
    amrex::Real scratch_116 = scratch_54*(-A[132]*scratch_110 + A[142] - scratch_111*scratch_49 - scratch_112*scratch_50 - scratch_113*scratch_51 - scratch_114*scratch_52 - scratch_115*scratch_53);
    amrex::Real scratch_117 = scratch_70*(-A[154]*scratch_110 + A[164] - scratch_111*scratch_64 - scratch_112*scratch_65 - scratch_113*scratch_66 - scratch_114*scratch_67 - scratch_115*scratch_68 - scratch_116*scratch_69);
    amrex::Real scratch_118 = scratch_88*(-A[176]*scratch_110 + A[185] - scratch_111*scratch_81 - scratch_112*scratch_82 - scratch_113*scratch_83 - scratch_114*scratch_84 - scratch_115*scratch_85 - scratch_116*scratch_86 - scratch_117*scratch_87);
    amrex::Real scratch_119 = scratch_108*(-A[194]*scratch_110 + A[204] - scratch_100*scratch_111 - scratch_101*scratch_112 - scratch_102*scratch_113 - scratch_103*scratch_114 - scratch_104*scratch_115 - scratch_105*scratch_116 - scratch_106*scratch_117 - scratch_107*scratch_118);
    amrex::Real scratch_120 = scratch_110 - scratch_111*scratch_4 - scratch_112*scratch_7 - scratch_113*scratch_15 - scratch_114*scratch_24 - scratch_115*scratch_35 - scratch_116*scratch_48 - scratch_117*scratch_63 - scratch_118*scratch_80 - scratch_119*scratch_99;
    amrex::Real scratch_121 = -A[216]*scratch_4 + A[217];
    amrex::Real scratch_122 = -A[216]*scratch_3 + A[218] - scratch_121*scratch_6;
    amrex::Real scratch_123 = -A[216]*scratch_12 + A[219] - scratch_121*scratch_13 - scratch_122*scratch_14;
    amrex::Real scratch_124 = -A[216]*scratch_20 + A[220] - scratch_121*scratch_21 - scratch_122*scratch_22 - scratch_123*scratch_23;
    amrex::Real scratch_125 = -A[216]*scratch_30 + A[221] - scratch_121*scratch_31 - scratch_122*scratch_32 - scratch_123*scratch_33 - scratch_124*scratch_34;
    amrex::Real scratch_126 = -A[216]*scratch_42 + A[222] - scratch_121*scratch_43 - scratch_122*scratch_44 - scratch_123*scratch_45 - scratch_124*scratch_46 - scratch_125*scratch_47;
    amrex::Real scratch_127 = -A[216]*scratch_56 + A[223] - scratch_121*scratch_57 - scratch_122*scratch_58 - scratch_123*scratch_59 - scratch_124*scratch_60 - scratch_125*scratch_61 - scratch_126*scratch_62;
    amrex::Real scratch_128 = -A[216]*scratch_72 + A[224] - scratch_121*scratch_73 - scratch_122*scratch_74 - scratch_123*scratch_75 - scratch_124*scratch_76 - scratch_125*scratch_77 - scratch_126*scratch_78 - scratch_127*scratch_79;
    amrex::Real scratch_129 = -A[216]*scratch_90 + A[225] - scratch_121*scratch_91 - scratch_122*scratch_92 - scratch_123*scratch_93 - scratch_124*scratch_94 - scratch_125*scratch_95 - scratch_126*scratch_96 - scratch_127*scratch_97 - scratch_128*scratch_98;
    amrex::Real scratch_130 = 1.0/(-A[216]*scratch_110 + A[226] - scratch_111*scratch_121 - scratch_112*scratch_122 - scratch_113*scratch_123 - scratch_114*scratch_124 - scratch_115*scratch_125 - scratch_116*scratch_126 - scratch_117*scratch_127 - scratch_118*scratch_128 - scratch_119*scratch_129);
    amrex::Real scratch_131 = scratch_130*(-A[216]*scratch_2 + A[235] - scratch_10*scratch_121 - scratch_109*scratch_129 - scratch_11*scratch_122 - scratch_123*scratch_19 - scratch_124*scratch_29 - scratch_125*scratch_41 - scratch_126*scratch_55 - scratch_127*scratch_71 - scratch_128*scratch_89);
    amrex::Real scratch_132 = A[11]*scratch_0;
    amrex::Real scratch_133 = scratch_5*(-A[22]*scratch_132 + A[33]);
    amrex::Real scratch_134 = scratch_9*(-A[44]*scratch_132 + A[55] - scratch_133*scratch_8);
    amrex::Real scratch_135 = scratch_18*(-A[66]*scratch_132 + A[77] - scratch_133*scratch_16 - scratch_134*scratch_17);
    amrex::Real scratch_136 = scratch_28*(A[99] - A[88]*scratch_132 - scratch_133*scratch_25 - scratch_134*scratch_26 - scratch_135*scratch_27);
    amrex::Real scratch_137 = scratch_40*(-A[110]*scratch_132 + A[121] - scratch_133*scratch_36 - scratch_134*scratch_37 - scratch_135*scratch_38 - scratch_136*scratch_39);
    amrex::Real scratch_138 = scratch_54*(-A[132]*scratch_132 + A[143] - scratch_133*scratch_49 - scratch_134*scratch_50 - scratch_135*scratch_51 - scratch_136*scratch_52 - scratch_137*scratch_53);
    amrex::Real scratch_139 = scratch_70*(-A[154]*scratch_132 + A[165] - scratch_133*scratch_64 - scratch_134*scratch_65 - scratch_135*scratch_66 - scratch_136*scratch_67 - scratch_137*scratch_68 - scratch_138*scratch_69);
    amrex::Real scratch_140 = scratch_88*(-A[176]*scratch_132 + A[186] - scratch_133*scratch_81 - scratch_134*scratch_82 - scratch_135*scratch_83 - scratch_136*scratch_84 - scratch_137*scratch_85 - scratch_138*scratch_86 - scratch_139*scratch_87);
    amrex::Real scratch_141 = scratch_108*(-A[194]*scratch_132 + A[205] - scratch_100*scratch_133 - scratch_101*scratch_134 - scratch_102*scratch_135 - scratch_103*scratch_136 - scratch_104*scratch_137 - scratch_105*scratch_138 - scratch_106*scratch_139 - scratch_107*scratch_140);
    amrex::Real scratch_142 = scratch_130*(-A[216]*scratch_132 + A[227] - scratch_121*scratch_133 - scratch_122*scratch_134 - scratch_123*scratch_135 - scratch_124*scratch_136 - scratch_125*scratch_137 - scratch_126*scratch_138 - scratch_127*scratch_139 - scratch_128*scratch_140 - scratch_129*scratch_141);
    amrex::Real scratch_143 = -scratch_120*scratch_142 + scratch_132 - scratch_133*scratch_4 - scratch_134*scratch_7 - scratch_135*scratch_15 - scratch_136*scratch_24 - scratch_137*scratch_35 - scratch_138*scratch_48 - scratch_139*scratch_63 - scratch_140*scratch_80 - scratch_141*scratch_99;
    amrex::Real scratch_144 = -A[238]*scratch_4 + A[239];
    amrex::Real scratch_145 = -A[238]*scratch_3 + A[240] - scratch_144*scratch_6;
    amrex::Real scratch_146 = -A[238]*scratch_12 + A[241] - scratch_13*scratch_144 - scratch_14*scratch_145;
    amrex::Real scratch_147 = -A[238]*scratch_20 + A[242] - scratch_144*scratch_21 - scratch_145*scratch_22 - scratch_146*scratch_23;
    amrex::Real scratch_148 = -A[238]*scratch_30 + A[243] - scratch_144*scratch_31 - scratch_145*scratch_32 - scratch_146*scratch_33 - scratch_147*scratch_34;
    amrex::Real scratch_149 = -A[238]*scratch_42 + A[244] - scratch_144*scratch_43 - scratch_145*scratch_44 - scratch_146*scratch_45 - scratch_147*scratch_46 - scratch_148*scratch_47;
    amrex::Real scratch_150 = -A[238]*scratch_56 + A[245] - scratch_144*scratch_57 - scratch_145*scratch_58 - scratch_146*scratch_59 - scratch_147*scratch_60 - scratch_148*scratch_61 - scratch_149*scratch_62;
    amrex::Real scratch_151 = -A[238]*scratch_72 + A[246] - scratch_144*scratch_73 - scratch_145*scratch_74 - scratch_146*scratch_75 - scratch_147*scratch_76 - scratch_148*scratch_77 - scratch_149*scratch_78 - scratch_150*scratch_79;
    amrex::Real scratch_152 = -A[238]*scratch_90 + A[247] - scratch_144*scratch_91 - scratch_145*scratch_92 - scratch_146*scratch_93 - scratch_147*scratch_94 - scratch_148*scratch_95 - scratch_149*scratch_96 - scratch_150*scratch_97 - scratch_151*scratch_98;
    amrex::Real scratch_153 = -A[238]*scratch_110 + A[248] - scratch_111*scratch_144 - scratch_112*scratch_145 - scratch_113*scratch_146 - scratch_114*scratch_147 - scratch_115*scratch_148 - scratch_116*scratch_149 - scratch_117*scratch_150 - scratch_118*scratch_151 - scratch_119*scratch_152;
    amrex::Real scratch_154 = 1.0/(-A[238]*scratch_132 + A[249] - scratch_133*scratch_144 - scratch_134*scratch_145 - scratch_135*scratch_146 - scratch_136*scratch_147 - scratch_137*scratch_148 - scratch_138*scratch_149 - scratch_139*scratch_150 - scratch_140*scratch_151 - scratch_141*scratch_152 - scratch_142*scratch_153);
    amrex::Real scratch_155 = scratch_154*(-A[238]*scratch_2 + A[257] - scratch_10*scratch_144 - scratch_109*scratch_152 - scratch_11*scratch_145 - scratch_131*scratch_153 - scratch_146*scratch_19 - scratch_147*scratch_29 - scratch_148*scratch_41 - scratch_149*scratch_55 - scratch_150*scratch_71 - scratch_151*scratch_89);
    amrex::Real scratch_156 = A[12]*scratch_0;
    amrex::Real scratch_157 = scratch_5*(-A[22]*scratch_156 + A[34]);
    amrex::Real scratch_158 = scratch_9*(-A[44]*scratch_156 + A[56] - scratch_157*scratch_8);
    amrex::Real scratch_159 = scratch_18*(-A[66]*scratch_156 + A[78] - scratch_157*scratch_16 - scratch_158*scratch_17);
    amrex::Real scratch_160 = scratch_28*(A[100] - A[88]*scratch_156 - scratch_157*scratch_25 - scratch_158*scratch_26 - scratch_159*scratch_27);
    amrex::Real scratch_161 = scratch_40*(-A[110]*scratch_156 + A[122] - scratch_157*scratch_36 - scratch_158*scratch_37 - scratch_159*scratch_38 - scratch_160*scratch_39);
    amrex::Real scratch_162 = scratch_54*(-A[132]*scratch_156 + A[144] - scratch_157*scratch_49 - scratch_158*scratch_50 - scratch_159*scratch_51 - scratch_160*scratch_52 - scratch_161*scratch_53);
    amrex::Real scratch_163 = scratch_70*(-A[154]*scratch_156 + A[166] - scratch_157*scratch_64 - scratch_158*scratch_65 - scratch_159*scratch_66 - scratch_160*scratch_67 - scratch_161*scratch_68 - scratch_162*scratch_69);
    amrex::Real scratch_164 = scratch_88*(-A[176]*scratch_156 + A[187] - scratch_157*scratch_81 - scratch_158*scratch_82 - scratch_159*scratch_83 - scratch_160*scratch_84 - scratch_161*scratch_85 - scratch_162*scratch_86 - scratch_163*scratch_87);
    amrex::Real scratch_165 = scratch_108*(-A[194]*scratch_156 + A[206] - scratch_100*scratch_157 - scratch_101*scratch_158 - scratch_102*scratch_159 - scratch_103*scratch_160 - scratch_104*scratch_161 - scratch_105*scratch_162 - scratch_106*scratch_163 - scratch_107*scratch_164);
    amrex::Real scratch_166 = scratch_130*(-A[216]*scratch_156 + A[228] - scratch_121*scratch_157 - scratch_122*scratch_158 - scratch_123*scratch_159 - scratch_124*scratch_160 - scratch_125*scratch_161 - scratch_126*scratch_162 - scratch_127*scratch_163 - scratch_128*scratch_164 - scratch_129*scratch_165);
    amrex::Real scratch_167 = scratch_154*(-A[238]*scratch_156 + A[250] - scratch_144*scratch_157 - scratch_145*scratch_158 - scratch_146*scratch_159 - scratch_147*scratch_160 - scratch_148*scratch_161 - scratch_149*scratch_162 - scratch_150*scratch_163 - scratch_151*scratch_164 - scratch_152*scratch_165 - scratch_153*scratch_166);
    amrex::Real scratch_168 = -scratch_120*scratch_166 - scratch_143*scratch_167 - scratch_15*scratch_159 + scratch_156 - scratch_157*scratch_4 - scratch_158*scratch_7 - scratch_160*scratch_24 - scratch_161*scratch_35 - scratch_162*scratch_48 - scratch_163*scratch_63 - scratch_164*scratch_80 - scratch_165*scratch_99;
    amrex::Real scratch_169 = -A[260]*scratch_4 + A[261];
    amrex::Real scratch_170 = -A[260]*scratch_3 + A[262] - scratch_169*scratch_6;
    amrex::Real scratch_171 = -A[260]*scratch_12 + A[263] - scratch_13*scratch_169 - scratch_14*scratch_170;
    amrex::Real scratch_172 = -A[260]*scratch_20 + A[264] - scratch_169*scratch_21 - scratch_170*scratch_22 - scratch_171*scratch_23;
    amrex::Real scratch_173 = -A[260]*scratch_30 + A[265] - scratch_169*scratch_31 - scratch_170*scratch_32 - scratch_171*scratch_33 - scratch_172*scratch_34;
    amrex::Real scratch_174 = -A[260]*scratch_42 + A[266] - scratch_169*scratch_43 - scratch_170*scratch_44 - scratch_171*scratch_45 - scratch_172*scratch_46 - scratch_173*scratch_47;
    amrex::Real scratch_175 = -A[260]*scratch_56 + A[267] - scratch_169*scratch_57 - scratch_170*scratch_58 - scratch_171*scratch_59 - scratch_172*scratch_60 - scratch_173*scratch_61 - scratch_174*scratch_62;
    amrex::Real scratch_176 = -A[260]*scratch_72 + A[268] - scratch_169*scratch_73 - scratch_170*scratch_74 - scratch_171*scratch_75 - scratch_172*scratch_76 - scratch_173*scratch_77 - scratch_174*scratch_78 - scratch_175*scratch_79;
    amrex::Real scratch_177 = -A[260]*scratch_90 + A[269] - scratch_169*scratch_91 - scratch_170*scratch_92 - scratch_171*scratch_93 - scratch_172*scratch_94 - scratch_173*scratch_95 - scratch_174*scratch_96 - scratch_175*scratch_97 - scratch_176*scratch_98;
    amrex::Real scratch_178 = -A[260]*scratch_110 + A[270] - scratch_111*scratch_169 - scratch_112*scratch_170 - scratch_113*scratch_171 - scratch_114*scratch_172 - scratch_115*scratch_173 - scratch_116*scratch_174 - scratch_117*scratch_175 - scratch_118*scratch_176 - scratch_119*scratch_177;
    amrex::Real scratch_179 = -A[260]*scratch_132 + A[271] - scratch_133*scratch_169 - scratch_134*scratch_170 - scratch_135*scratch_171 - scratch_136*scratch_172 - scratch_137*scratch_173 - scratch_138*scratch_174 - scratch_139*scratch_175 - scratch_140*scratch_176 - scratch_141*scratch_177 - scratch_142*scratch_178;
    amrex::Real scratch_180 = 1.0/(-A[260]*scratch_156 + A[272] - scratch_157*scratch_169 - scratch_158*scratch_170 - scratch_159*scratch_171 - scratch_160*scratch_172 - scratch_161*scratch_173 - scratch_162*scratch_174 - scratch_163*scratch_175 - scratch_164*scratch_176 - scratch_165*scratch_177 - scratch_166*scratch_178 - scratch_167*scratch_179);
    amrex::Real scratch_181 = scratch_180*(-A[260]*scratch_2 + A[279] - scratch_10*scratch_169 - scratch_109*scratch_177 - scratch_11*scratch_170 - scratch_131*scratch_178 - scratch_155*scratch_179 - scratch_171*scratch_19 - scratch_172*scratch_29 - scratch_173*scratch_41 - scratch_174*scratch_55 - scratch_175*scratch_71 - scratch_176*scratch_89);
    amrex::Real scratch_182 = A[13]*scratch_0;
    amrex::Real scratch_183 = scratch_5*(-A[22]*scratch_182 + A[35]);
    amrex::Real scratch_184 = scratch_9*(-A[44]*scratch_182 + A[57] - scratch_183*scratch_8);
    amrex::Real scratch_185 = scratch_18*(-A[66]*scratch_182 + A[79] - scratch_16*scratch_183 - scratch_17*scratch_184);
    amrex::Real scratch_186 = scratch_28*(A[101] - A[88]*scratch_182 - scratch_183*scratch_25 - scratch_184*scratch_26 - scratch_185*scratch_27);
    amrex::Real scratch_187 = scratch_40*(-A[110]*scratch_182 + A[123] - scratch_183*scratch_36 - scratch_184*scratch_37 - scratch_185*scratch_38 - scratch_186*scratch_39);
    amrex::Real scratch_188 = scratch_54*(-A[132]*scratch_182 + A[145] - scratch_183*scratch_49 - scratch_184*scratch_50 - scratch_185*scratch_51 - scratch_186*scratch_52 - scratch_187*scratch_53);
    amrex::Real scratch_189 = scratch_70*(-A[154]*scratch_182 + A[167] - scratch_183*scratch_64 - scratch_184*scratch_65 - scratch_185*scratch_66 - scratch_186*scratch_67 - scratch_187*scratch_68 - scratch_188*scratch_69);
    amrex::Real scratch_190 = scratch_88*(-A[176]*scratch_182 + A[188] - scratch_183*scratch_81 - scratch_184*scratch_82 - scratch_185*scratch_83 - scratch_186*scratch_84 - scratch_187*scratch_85 - scratch_188*scratch_86 - scratch_189*scratch_87);
    amrex::Real scratch_191 = scratch_108*(-A[194]*scratch_182 + A[207] - scratch_100*scratch_183 - scratch_101*scratch_184 - scratch_102*scratch_185 - scratch_103*scratch_186 - scratch_104*scratch_187 - scratch_105*scratch_188 - scratch_106*scratch_189 - scratch_107*scratch_190);
    amrex::Real scratch_192 = scratch_130*(-A[216]*scratch_182 + A[229] - scratch_121*scratch_183 - scratch_122*scratch_184 - scratch_123*scratch_185 - scratch_124*scratch_186 - scratch_125*scratch_187 - scratch_126*scratch_188 - scratch_127*scratch_189 - scratch_128*scratch_190 - scratch_129*scratch_191);
    amrex::Real scratch_193 = scratch_154*(-A[238]*scratch_182 + A[251] - scratch_144*scratch_183 - scratch_145*scratch_184 - scratch_146*scratch_185 - scratch_147*scratch_186 - scratch_148*scratch_187 - scratch_149*scratch_188 - scratch_150*scratch_189 - scratch_151*scratch_190 - scratch_152*scratch_191 - scratch_153*scratch_192);
    amrex::Real scratch_194 = scratch_180*(-A[260]*scratch_182 + A[273] - scratch_169*scratch_183 - scratch_170*scratch_184 - scratch_171*scratch_185 - scratch_172*scratch_186 - scratch_173*scratch_187 - scratch_174*scratch_188 - scratch_175*scratch_189 - scratch_176*scratch_190 - scratch_177*scratch_191 - scratch_178*scratch_192 - scratch_179*scratch_193);
    amrex::Real scratch_195 = -scratch_120*scratch_192 - scratch_143*scratch_193 - scratch_15*scratch_185 - scratch_168*scratch_194 + scratch_182 - scratch_183*scratch_4 - scratch_184*scratch_7 - scratch_186*scratch_24 - scratch_187*scratch_35 - scratch_188*scratch_48 - scratch_189*scratch_63 - scratch_190*scratch_80 - scratch_191*scratch_99;
    amrex::Real scratch_196 = -A[282]*scratch_4 + A[283];
    amrex::Real scratch_197 = -A[282]*scratch_3 + A[284] - scratch_196*scratch_6;
    amrex::Real scratch_198 = -A[282]*scratch_12 + A[285] - scratch_13*scratch_196 - scratch_14*scratch_197;
    amrex::Real scratch_199 = -A[282]*scratch_20 + A[286] - scratch_196*scratch_21 - scratch_197*scratch_22 - scratch_198*scratch_23;
    amrex::Real scratch_200 = -A[282]*scratch_30 + A[287] - scratch_196*scratch_31 - scratch_197*scratch_32 - scratch_198*scratch_33 - scratch_199*scratch_34;
    amrex::Real scratch_201 = -A[282]*scratch_42 + A[288] - scratch_196*scratch_43 - scratch_197*scratch_44 - scratch_198*scratch_45 - scratch_199*scratch_46 - scratch_200*scratch_47;
    amrex::Real scratch_202 = -A[282]*scratch_56 + A[289] - scratch_196*scratch_57 - scratch_197*scratch_58 - scratch_198*scratch_59 - scratch_199*scratch_60 - scratch_200*scratch_61 - scratch_201*scratch_62;
    amrex::Real scratch_203 = -A[282]*scratch_72 + A[290] - scratch_196*scratch_73 - scratch_197*scratch_74 - scratch_198*scratch_75 - scratch_199*scratch_76 - scratch_200*scratch_77 - scratch_201*scratch_78 - scratch_202*scratch_79;
    amrex::Real scratch_204 = -A[282]*scratch_90 + A[291] - scratch_196*scratch_91 - scratch_197*scratch_92 - scratch_198*scratch_93 - scratch_199*scratch_94 - scratch_200*scratch_95 - scratch_201*scratch_96 - scratch_202*scratch_97 - scratch_203*scratch_98;
    amrex::Real scratch_205 = -A[282]*scratch_110 + A[292] - scratch_111*scratch_196 - scratch_112*scratch_197 - scratch_113*scratch_198 - scratch_114*scratch_199 - scratch_115*scratch_200 - scratch_116*scratch_201 - scratch_117*scratch_202 - scratch_118*scratch_203 - scratch_119*scratch_204;
    amrex::Real scratch_206 = -A[282]*scratch_132 + A[293] - scratch_133*scratch_196 - scratch_134*scratch_197 - scratch_135*scratch_198 - scratch_136*scratch_199 - scratch_137*scratch_200 - scratch_138*scratch_201 - scratch_139*scratch_202 - scratch_140*scratch_203 - scratch_141*scratch_204 - scratch_142*scratch_205;
    amrex::Real scratch_207 = -A[282]*scratch_156 + A[294] - scratch_157*scratch_196 - scratch_158*scratch_197 - scratch_159*scratch_198 - scratch_160*scratch_199 - scratch_161*scratch_200 - scratch_162*scratch_201 - scratch_163*scratch_202 - scratch_164*scratch_203 - scratch_165*scratch_204 - scratch_166*scratch_205 - scratch_167*scratch_206;
    amrex::Real scratch_208 = 1.0/(-A[282]*scratch_182 + A[295] - scratch_183*scratch_196 - scratch_184*scratch_197 - scratch_185*scratch_198 - scratch_186*scratch_199 - scratch_187*scratch_200 - scratch_188*scratch_201 - scratch_189*scratch_202 - scratch_190*scratch_203 - scratch_191*scratch_204 - scratch_192*scratch_205 - scratch_193*scratch_206 - scratch_194*scratch_207);
    amrex::Real scratch_209 = scratch_208*(-A[282]*scratch_2 + A[301] - scratch_10*scratch_196 - scratch_109*scratch_204 - scratch_11*scratch_197 - scratch_131*scratch_205 - scratch_155*scratch_206 - scratch_181*scratch_207 - scratch_19*scratch_198 - scratch_199*scratch_29 - scratch_200*scratch_41 - scratch_201*scratch_55 - scratch_202*scratch_71 - scratch_203*scratch_89);
    amrex::Real scratch_210 = A[14]*scratch_0;
    amrex::Real scratch_211 = scratch_5*(-A[22]*scratch_210 + A[36]);
    amrex::Real scratch_212 = scratch_9*(-A[44]*scratch_210 + A[58] - scratch_211*scratch_8);
    amrex::Real scratch_213 = scratch_18*(-A[66]*scratch_210 + A[80] - scratch_16*scratch_211 - scratch_17*scratch_212);
    amrex::Real scratch_214 = scratch_28*(A[102] - A[88]*scratch_210 - scratch_211*scratch_25 - scratch_212*scratch_26 - scratch_213*scratch_27);
    amrex::Real scratch_215 = scratch_40*(-A[110]*scratch_210 + A[124] - scratch_211*scratch_36 - scratch_212*scratch_37 - scratch_213*scratch_38 - scratch_214*scratch_39);
    amrex::Real scratch_216 = scratch_54*(-A[132]*scratch_210 + A[146] - scratch_211*scratch_49 - scratch_212*scratch_50 - scratch_213*scratch_51 - scratch_214*scratch_52 - scratch_215*scratch_53);
    amrex::Real scratch_217 = scratch_70*(-A[154]*scratch_210 + A[168] - scratch_211*scratch_64 - scratch_212*scratch_65 - scratch_213*scratch_66 - scratch_214*scratch_67 - scratch_215*scratch_68 - scratch_216*scratch_69);
    amrex::Real scratch_218 = scratch_88*(-A[176]*scratch_210 + A[189] - scratch_211*scratch_81 - scratch_212*scratch_82 - scratch_213*scratch_83 - scratch_214*scratch_84 - scratch_215*scratch_85 - scratch_216*scratch_86 - scratch_217*scratch_87);
    amrex::Real scratch_219 = scratch_108*(-A[194]*scratch_210 + A[208] - scratch_100*scratch_211 - scratch_101*scratch_212 - scratch_102*scratch_213 - scratch_103*scratch_214 - scratch_104*scratch_215 - scratch_105*scratch_216 - scratch_106*scratch_217 - scratch_107*scratch_218);
    amrex::Real scratch_220 = scratch_130*(-A[216]*scratch_210 + A[230] - scratch_121*scratch_211 - scratch_122*scratch_212 - scratch_123*scratch_213 - scratch_124*scratch_214 - scratch_125*scratch_215 - scratch_126*scratch_216 - scratch_127*scratch_217 - scratch_128*scratch_218 - scratch_129*scratch_219);
    amrex::Real scratch_221 = scratch_154*(-A[238]*scratch_210 + A[252] - scratch_144*scratch_211 - scratch_145*scratch_212 - scratch_146*scratch_213 - scratch_147*scratch_214 - scratch_148*scratch_215 - scratch_149*scratch_216 - scratch_150*scratch_217 - scratch_151*scratch_218 - scratch_152*scratch_219 - scratch_153*scratch_220);
    amrex::Real scratch_222 = scratch_180*(-A[260]*scratch_210 + A[274] - scratch_169*scratch_211 - scratch_170*scratch_212 - scratch_171*scratch_213 - scratch_172*scratch_214 - scratch_173*scratch_215 - scratch_174*scratch_216 - scratch_175*scratch_217 - scratch_176*scratch_218 - scratch_177*scratch_219 - scratch_178*scratch_220 - scratch_179*scratch_221);
    amrex::Real scratch_223 = scratch_208*(-A[282]*scratch_210 + A[296] - scratch_196*scratch_211 - scratch_197*scratch_212 - scratch_198*scratch_213 - scratch_199*scratch_214 - scratch_200*scratch_215 - scratch_201*scratch_216 - scratch_202*scratch_217 - scratch_203*scratch_218 - scratch_204*scratch_219 - scratch_205*scratch_220 - scratch_206*scratch_221 - scratch_207*scratch_222);
    amrex::Real scratch_224 = -scratch_120*scratch_220 - scratch_143*scratch_221 - scratch_15*scratch_213 - scratch_168*scratch_222 - scratch_195*scratch_223 + scratch_210 - scratch_211*scratch_4 - scratch_212*scratch_7 - scratch_214*scratch_24 - scratch_215*scratch_35 - scratch_216*scratch_48 - scratch_217*scratch_63 - scratch_218*scratch_80 - scratch_219*scratch_99;
    amrex::Real scratch_225 = -A[304]*scratch_4 + A[305];
    amrex::Real scratch_226 = -A[304]*scratch_3 + A[306] - scratch_225*scratch_6;
    amrex::Real scratch_227 = -A[304]*scratch_12 + A[307] - scratch_13*scratch_225 - scratch_14*scratch_226;
    amrex::Real scratch_228 = -A[304]*scratch_20 + A[308] - scratch_21*scratch_225 - scratch_22*scratch_226 - scratch_227*scratch_23;
    amrex::Real scratch_229 = -A[304]*scratch_30 + A[309] - scratch_225*scratch_31 - scratch_226*scratch_32 - scratch_227*scratch_33 - scratch_228*scratch_34;
    amrex::Real scratch_230 = -A[304]*scratch_42 + A[310] - scratch_225*scratch_43 - scratch_226*scratch_44 - scratch_227*scratch_45 - scratch_228*scratch_46 - scratch_229*scratch_47;
    amrex::Real scratch_231 = -A[304]*scratch_56 + A[311] - scratch_225*scratch_57 - scratch_226*scratch_58 - scratch_227*scratch_59 - scratch_228*scratch_60 - scratch_229*scratch_61 - scratch_230*scratch_62;
    amrex::Real scratch_232 = -A[304]*scratch_72 + A[312] - scratch_225*scratch_73 - scratch_226*scratch_74 - scratch_227*scratch_75 - scratch_228*scratch_76 - scratch_229*scratch_77 - scratch_230*scratch_78 - scratch_231*scratch_79;
    amrex::Real scratch_233 = -A[304]*scratch_90 + A[313] - scratch_225*scratch_91 - scratch_226*scratch_92 - scratch_227*scratch_93 - scratch_228*scratch_94 - scratch_229*scratch_95 - scratch_230*scratch_96 - scratch_231*scratch_97 - scratch_232*scratch_98;
    amrex::Real scratch_234 = -A[304]*scratch_110 + A[314] - scratch_111*scratch_225 - scratch_112*scratch_226 - scratch_113*scratch_227 - scratch_114*scratch_228 - scratch_115*scratch_229 - scratch_116*scratch_230 - scratch_117*scratch_231 - scratch_118*scratch_232 - scratch_119*scratch_233;
    amrex::Real scratch_235 = -A[304]*scratch_132 + A[315] - scratch_133*scratch_225 - scratch_134*scratch_226 - scratch_135*scratch_227 - scratch_136*scratch_228 - scratch_137*scratch_229 - scratch_138*scratch_230 - scratch_139*scratch_231 - scratch_140*scratch_232 - scratch_141*scratch_233 - scratch_142*scratch_234;
    amrex::Real scratch_236 = -A[304]*scratch_156 + A[316] - scratch_157*scratch_225 - scratch_158*scratch_226 - scratch_159*scratch_227 - scratch_160*scratch_228 - scratch_161*scratch_229 - scratch_162*scratch_230 - scratch_163*scratch_231 - scratch_164*scratch_232 - scratch_165*scratch_233 - scratch_166*scratch_234 - scratch_167*scratch_235;
    amrex::Real scratch_237 = -A[304]*scratch_182 + A[317] - scratch_183*scratch_225 - scratch_184*scratch_226 - scratch_185*scratch_227 - scratch_186*scratch_228 - scratch_187*scratch_229 - scratch_188*scratch_230 - scratch_189*scratch_231 - scratch_190*scratch_232 - scratch_191*scratch_233 - scratch_192*scratch_234 - scratch_193*scratch_235 - scratch_194*scratch_236;
    amrex::Real scratch_238 = 1.0/(-A[304]*scratch_210 + A[318] - scratch_211*scratch_225 - scratch_212*scratch_226 - scratch_213*scratch_227 - scratch_214*scratch_228 - scratch_215*scratch_229 - scratch_216*scratch_230 - scratch_217*scratch_231 - scratch_218*scratch_232 - scratch_219*scratch_233 - scratch_220*scratch_234 - scratch_221*scratch_235 - scratch_222*scratch_236 - scratch_223*scratch_237);
    amrex::Real scratch_239 = scratch_238*(-A[304]*scratch_2 + A[323] - scratch_10*scratch_225 - scratch_109*scratch_233 - scratch_11*scratch_226 - scratch_131*scratch_234 - scratch_155*scratch_235 - scratch_181*scratch_236 - scratch_19*scratch_227 - scratch_209*scratch_237 - scratch_228*scratch_29 - scratch_229*scratch_41 - scratch_230*scratch_55 - scratch_231*scratch_71 - scratch_232*scratch_89);
    amrex::Real scratch_240 = A[15]*scratch_0;
    amrex::Real scratch_241 = scratch_5*(-A[22]*scratch_240 + A[37]);
    amrex::Real scratch_242 = scratch_9*(-A[44]*scratch_240 + A[59] - scratch_241*scratch_8);
    amrex::Real scratch_243 = scratch_18*(-A[66]*scratch_240 + A[81] - scratch_16*scratch_241 - scratch_17*scratch_242);
    amrex::Real scratch_244 = scratch_28*(A[103] - A[88]*scratch_240 - scratch_241*scratch_25 - scratch_242*scratch_26 - scratch_243*scratch_27);
    amrex::Real scratch_245 = scratch_40*(-A[110]*scratch_240 + A[125] - scratch_241*scratch_36 - scratch_242*scratch_37 - scratch_243*scratch_38 - scratch_244*scratch_39);
    amrex::Real scratch_246 = scratch_54*(-A[132]*scratch_240 + A[147] - scratch_241*scratch_49 - scratch_242*scratch_50 - scratch_243*scratch_51 - scratch_244*scratch_52 - scratch_245*scratch_53);
    amrex::Real scratch_247 = scratch_70*(-A[154]*scratch_240 + A[169] - scratch_241*scratch_64 - scratch_242*scratch_65 - scratch_243*scratch_66 - scratch_244*scratch_67 - scratch_245*scratch_68 - scratch_246*scratch_69);
    amrex::Real scratch_248 = scratch_88*(-A[176]*scratch_240 - scratch_241*scratch_81 - scratch_242*scratch_82 - scratch_243*scratch_83 - scratch_244*scratch_84 - scratch_245*scratch_85 - scratch_246*scratch_86 - scratch_247*scratch_87);
    amrex::Real scratch_249 = scratch_108*(-A[194]*scratch_240 + A[209] - scratch_100*scratch_241 - scratch_101*scratch_242 - scratch_102*scratch_243 - scratch_103*scratch_244 - scratch_104*scratch_245 - scratch_105*scratch_246 - scratch_106*scratch_247 - scratch_107*scratch_248);
    amrex::Real scratch_250 = scratch_130*(-A[216]*scratch_240 + A[231] - scratch_121*scratch_241 - scratch_122*scratch_242 - scratch_123*scratch_243 - scratch_124*scratch_244 - scratch_125*scratch_245 - scratch_126*scratch_246 - scratch_127*scratch_247 - scratch_128*scratch_248 - scratch_129*scratch_249);
    amrex::Real scratch_251 = scratch_154*(-A[238]*scratch_240 + A[253] - scratch_144*scratch_241 - scratch_145*scratch_242 - scratch_146*scratch_243 - scratch_147*scratch_244 - scratch_148*scratch_245 - scratch_149*scratch_246 - scratch_150*scratch_247 - scratch_151*scratch_248 - scratch_152*scratch_249 - scratch_153*scratch_250);
    amrex::Real scratch_252 = scratch_180*(-A[260]*scratch_240 + A[275] - scratch_169*scratch_241 - scratch_170*scratch_242 - scratch_171*scratch_243 - scratch_172*scratch_244 - scratch_173*scratch_245 - scratch_174*scratch_246 - scratch_175*scratch_247 - scratch_176*scratch_248 - scratch_177*scratch_249 - scratch_178*scratch_250 - scratch_179*scratch_251);
    amrex::Real scratch_253 = scratch_208*(-A[282]*scratch_240 + A[297] - scratch_196*scratch_241 - scratch_197*scratch_242 - scratch_198*scratch_243 - scratch_199*scratch_244 - scratch_200*scratch_245 - scratch_201*scratch_246 - scratch_202*scratch_247 - scratch_203*scratch_248 - scratch_204*scratch_249 - scratch_205*scratch_250 - scratch_206*scratch_251 - scratch_207*scratch_252);
    amrex::Real scratch_254 = scratch_238*(-A[304]*scratch_240 + A[319] - scratch_225*scratch_241 - scratch_226*scratch_242 - scratch_227*scratch_243 - scratch_228*scratch_244 - scratch_229*scratch_245 - scratch_230*scratch_246 - scratch_231*scratch_247 - scratch_232*scratch_248 - scratch_233*scratch_249 - scratch_234*scratch_250 - scratch_235*scratch_251 - scratch_236*scratch_252 - scratch_237*scratch_253);
    amrex::Real scratch_255 = -scratch_120*scratch_250 - scratch_143*scratch_251 - scratch_15*scratch_243 - scratch_168*scratch_252 - scratch_195*scratch_253 - scratch_224*scratch_254 - scratch_24*scratch_244 + scratch_240 - scratch_241*scratch_4 - scratch_242*scratch_7 - scratch_245*scratch_35 - scratch_246*scratch_48 - scratch_247*scratch_63 - scratch_248*scratch_80 - scratch_249*scratch_99;
    amrex::Real scratch_256 = -A[326]*scratch_4 + A[327];
    amrex::Real scratch_257 = -A[326]*scratch_3 + A[328] - scratch_256*scratch_6;
    amrex::Real scratch_258 = -A[326]*scratch_12 + A[329] - scratch_13*scratch_256 - scratch_14*scratch_257;
    amrex::Real scratch_259 = -A[326]*scratch_20 + A[330] - scratch_21*scratch_256 - scratch_22*scratch_257 - scratch_23*scratch_258;
    amrex::Real scratch_260 = -A[326]*scratch_30 + A[331] - scratch_256*scratch_31 - scratch_257*scratch_32 - scratch_258*scratch_33 - scratch_259*scratch_34;
    amrex::Real scratch_261 = -A[326]*scratch_42 + A[332] - scratch_256*scratch_43 - scratch_257*scratch_44 - scratch_258*scratch_45 - scratch_259*scratch_46 - scratch_260*scratch_47;
    amrex::Real scratch_262 = -A[326]*scratch_56 + A[333] - scratch_256*scratch_57 - scratch_257*scratch_58 - scratch_258*scratch_59 - scratch_259*scratch_60 - scratch_260*scratch_61 - scratch_261*scratch_62;
    amrex::Real scratch_263 = -A[326]*scratch_72 + A[334] - scratch_256*scratch_73 - scratch_257*scratch_74 - scratch_258*scratch_75 - scratch_259*scratch_76 - scratch_260*scratch_77 - scratch_261*scratch_78 - scratch_262*scratch_79;
    amrex::Real scratch_264 = -A[326]*scratch_90 + A[335] - scratch_256*scratch_91 - scratch_257*scratch_92 - scratch_258*scratch_93 - scratch_259*scratch_94 - scratch_260*scratch_95 - scratch_261*scratch_96 - scratch_262*scratch_97 - scratch_263*scratch_98;
    amrex::Real scratch_265 = -A[326]*scratch_110 + A[336] - scratch_111*scratch_256 - scratch_112*scratch_257 - scratch_113*scratch_258 - scratch_114*scratch_259 - scratch_115*scratch_260 - scratch_116*scratch_261 - scratch_117*scratch_262 - scratch_118*scratch_263 - scratch_119*scratch_264;
    amrex::Real scratch_266 = -A[326]*scratch_132 + A[337] - scratch_133*scratch_256 - scratch_134*scratch_257 - scratch_135*scratch_258 - scratch_136*scratch_259 - scratch_137*scratch_260 - scratch_138*scratch_261 - scratch_139*scratch_262 - scratch_140*scratch_263 - scratch_141*scratch_264 - scratch_142*scratch_265;
    amrex::Real scratch_267 = -A[326]*scratch_156 + A[338] - scratch_157*scratch_256 - scratch_158*scratch_257 - scratch_159*scratch_258 - scratch_160*scratch_259 - scratch_161*scratch_260 - scratch_162*scratch_261 - scratch_163*scratch_262 - scratch_164*scratch_263 - scratch_165*scratch_264 - scratch_166*scratch_265 - scratch_167*scratch_266;
    amrex::Real scratch_268 = -A[326]*scratch_182 + A[339] - scratch_183*scratch_256 - scratch_184*scratch_257 - scratch_185*scratch_258 - scratch_186*scratch_259 - scratch_187*scratch_260 - scratch_188*scratch_261 - scratch_189*scratch_262 - scratch_190*scratch_263 - scratch_191*scratch_264 - scratch_192*scratch_265 - scratch_193*scratch_266 - scratch_194*scratch_267;
    amrex::Real scratch_269 = -A[326]*scratch_210 + A[340] - scratch_211*scratch_256 - scratch_212*scratch_257 - scratch_213*scratch_258 - scratch_214*scratch_259 - scratch_215*scratch_260 - scratch_216*scratch_261 - scratch_217*scratch_262 - scratch_218*scratch_263 - scratch_219*scratch_264 - scratch_220*scratch_265 - scratch_221*scratch_266 - scratch_222*scratch_267 - scratch_223*scratch_268;
    amrex::Real scratch_270 = 1.0/(-A[326]*scratch_240 + A[341] - scratch_241*scratch_256 - scratch_242*scratch_257 - scratch_243*scratch_258 - scratch_244*scratch_259 - scratch_245*scratch_260 - scratch_246*scratch_261 - scratch_247*scratch_262 - scratch_248*scratch_263 - scratch_249*scratch_264 - scratch_250*scratch_265 - scratch_251*scratch_266 - scratch_252*scratch_267 - scratch_253*scratch_268 - scratch_254*scratch_269);
    amrex::Real scratch_271 = scratch_270*(-A[326]*scratch_2 + A[345] - scratch_10*scratch_256 - scratch_109*scratch_264 - scratch_11*scratch_257 - scratch_131*scratch_265 - scratch_155*scratch_266 - scratch_181*scratch_267 - scratch_19*scratch_258 - scratch_209*scratch_268 - scratch_239*scratch_269 - scratch_259*scratch_29 - scratch_260*scratch_41 - scratch_261*scratch_55 - scratch_262*scratch_71 - scratch_263*scratch_89);
    amrex::Real scratch_272 = A[16]*scratch_0;
    amrex::Real scratch_273 = scratch_5*(-A[22]*scratch_272 + A[38]);
    amrex::Real scratch_274 = scratch_9*(-A[44]*scratch_272 + A[60] - scratch_273*scratch_8);
    amrex::Real scratch_275 = scratch_18*(-A[66]*scratch_272 + A[82] - scratch_16*scratch_273 - scratch_17*scratch_274);
    amrex::Real scratch_276 = scratch_28*(A[104] - A[88]*scratch_272 - scratch_25*scratch_273 - scratch_26*scratch_274 - scratch_27*scratch_275);
    amrex::Real scratch_277 = scratch_40*(-A[110]*scratch_272 + A[126] - scratch_273*scratch_36 - scratch_274*scratch_37 - scratch_275*scratch_38 - scratch_276*scratch_39);
    amrex::Real scratch_278 = scratch_54*(-A[132]*scratch_272 + A[148] - scratch_273*scratch_49 - scratch_274*scratch_50 - scratch_275*scratch_51 - scratch_276*scratch_52 - scratch_277*scratch_53);
    amrex::Real scratch_279 = scratch_70*(-A[154]*scratch_272 + A[170] - scratch_273*scratch_64 - scratch_274*scratch_65 - scratch_275*scratch_66 - scratch_276*scratch_67 - scratch_277*scratch_68 - scratch_278*scratch_69);
    amrex::Real scratch_280 = scratch_88*(-A[176]*scratch_272 + A[190] - scratch_273*scratch_81 - scratch_274*scratch_82 - scratch_275*scratch_83 - scratch_276*scratch_84 - scratch_277*scratch_85 - scratch_278*scratch_86 - scratch_279*scratch_87);
    amrex::Real scratch_281 = scratch_108*(-A[194]*scratch_272 + A[210] - scratch_100*scratch_273 - scratch_101*scratch_274 - scratch_102*scratch_275 - scratch_103*scratch_276 - scratch_104*scratch_277 - scratch_105*scratch_278 - scratch_106*scratch_279 - scratch_107*scratch_280);
    amrex::Real scratch_282 = scratch_130*(-A[216]*scratch_272 + A[232] - scratch_121*scratch_273 - scratch_122*scratch_274 - scratch_123*scratch_275 - scratch_124*scratch_276 - scratch_125*scratch_277 - scratch_126*scratch_278 - scratch_127*scratch_279 - scratch_128*scratch_280 - scratch_129*scratch_281);
    amrex::Real scratch_283 = scratch_154*(-A[238]*scratch_272 + A[254] - scratch_144*scratch_273 - scratch_145*scratch_274 - scratch_146*scratch_275 - scratch_147*scratch_276 - scratch_148*scratch_277 - scratch_149*scratch_278 - scratch_150*scratch_279 - scratch_151*scratch_280 - scratch_152*scratch_281 - scratch_153*scratch_282);
    amrex::Real scratch_284 = scratch_180*(-A[260]*scratch_272 + A[276] - scratch_169*scratch_273 - scratch_170*scratch_274 - scratch_171*scratch_275 - scratch_172*scratch_276 - scratch_173*scratch_277 - scratch_174*scratch_278 - scratch_175*scratch_279 - scratch_176*scratch_280 - scratch_177*scratch_281 - scratch_178*scratch_282 - scratch_179*scratch_283);
    amrex::Real scratch_285 = scratch_208*(-A[282]*scratch_272 + A[298] - scratch_196*scratch_273 - scratch_197*scratch_274 - scratch_198*scratch_275 - scratch_199*scratch_276 - scratch_200*scratch_277 - scratch_201*scratch_278 - scratch_202*scratch_279 - scratch_203*scratch_280 - scratch_204*scratch_281 - scratch_205*scratch_282 - scratch_206*scratch_283 - scratch_207*scratch_284);
    amrex::Real scratch_286 = scratch_238*(-A[304]*scratch_272 + A[320] - scratch_225*scratch_273 - scratch_226*scratch_274 - scratch_227*scratch_275 - scratch_228*scratch_276 - scratch_229*scratch_277 - scratch_230*scratch_278 - scratch_231*scratch_279 - scratch_232*scratch_280 - scratch_233*scratch_281 - scratch_234*scratch_282 - scratch_235*scratch_283 - scratch_236*scratch_284 - scratch_237*scratch_285);
    amrex::Real scratch_287 = scratch_270*(-A[326]*scratch_272 + A[342] - scratch_256*scratch_273 - scratch_257*scratch_274 - scratch_258*scratch_275 - scratch_259*scratch_276 - scratch_260*scratch_277 - scratch_261*scratch_278 - scratch_262*scratch_279 - scratch_263*scratch_280 - scratch_264*scratch_281 - scratch_265*scratch_282 - scratch_266*scratch_283 - scratch_267*scratch_284 - scratch_268*scratch_285 - scratch_269*scratch_286);
    amrex::Real scratch_288 = -scratch_120*scratch_282 - scratch_143*scratch_283 - scratch_15*scratch_275 - scratch_168*scratch_284 - scratch_195*scratch_285 - scratch_224*scratch_286 - scratch_24*scratch_276 - scratch_255*scratch_287 + scratch_272 - scratch_273*scratch_4 - scratch_274*scratch_7 - scratch_277*scratch_35 - scratch_278*scratch_48 - scratch_279*scratch_63 - scratch_280*scratch_80 - scratch_281*scratch_99;
    amrex::Real scratch_289 = -A[348]*scratch_4 + A[349];
    amrex::Real scratch_290 = -A[348]*scratch_3 + A[350] - scratch_289*scratch_6;
    amrex::Real scratch_291 = -A[348]*scratch_12 + A[351] - scratch_13*scratch_289 - scratch_14*scratch_290;
    amrex::Real scratch_292 = -A[348]*scratch_20 + A[352] - scratch_21*scratch_289 - scratch_22*scratch_290 - scratch_23*scratch_291;
    amrex::Real scratch_293 = -A[348]*scratch_30 + A[353] - scratch_289*scratch_31 - scratch_290*scratch_32 - scratch_291*scratch_33 - scratch_292*scratch_34;
    amrex::Real scratch_294 = -A[348]*scratch_42 + A[354] - scratch_289*scratch_43 - scratch_290*scratch_44 - scratch_291*scratch_45 - scratch_292*scratch_46 - scratch_293*scratch_47;
    amrex::Real scratch_295 = -A[348]*scratch_56 + A[355] - scratch_289*scratch_57 - scratch_290*scratch_58 - scratch_291*scratch_59 - scratch_292*scratch_60 - scratch_293*scratch_61 - scratch_294*scratch_62;
    amrex::Real scratch_296 = -A[348]*scratch_72 + A[356] - scratch_289*scratch_73 - scratch_290*scratch_74 - scratch_291*scratch_75 - scratch_292*scratch_76 - scratch_293*scratch_77 - scratch_294*scratch_78 - scratch_295*scratch_79;
    amrex::Real scratch_297 = -A[348]*scratch_90 + A[357] - scratch_289*scratch_91 - scratch_290*scratch_92 - scratch_291*scratch_93 - scratch_292*scratch_94 - scratch_293*scratch_95 - scratch_294*scratch_96 - scratch_295*scratch_97 - scratch_296*scratch_98;
    amrex::Real scratch_298 = -A[348]*scratch_110 + A[358] - scratch_111*scratch_289 - scratch_112*scratch_290 - scratch_113*scratch_291 - scratch_114*scratch_292 - scratch_115*scratch_293 - scratch_116*scratch_294 - scratch_117*scratch_295 - scratch_118*scratch_296 - scratch_119*scratch_297;
    amrex::Real scratch_299 = -A[348]*scratch_132 + A[359] - scratch_133*scratch_289 - scratch_134*scratch_290 - scratch_135*scratch_291 - scratch_136*scratch_292 - scratch_137*scratch_293 - scratch_138*scratch_294 - scratch_139*scratch_295 - scratch_140*scratch_296 - scratch_141*scratch_297 - scratch_142*scratch_298;
    amrex::Real scratch_300 = -A[348]*scratch_156 + A[360] - scratch_157*scratch_289 - scratch_158*scratch_290 - scratch_159*scratch_291 - scratch_160*scratch_292 - scratch_161*scratch_293 - scratch_162*scratch_294 - scratch_163*scratch_295 - scratch_164*scratch_296 - scratch_165*scratch_297 - scratch_166*scratch_298 - scratch_167*scratch_299;
    amrex::Real scratch_301 = -A[348]*scratch_182 + A[361] - scratch_183*scratch_289 - scratch_184*scratch_290 - scratch_185*scratch_291 - scratch_186*scratch_292 - scratch_187*scratch_293 - scratch_188*scratch_294 - scratch_189*scratch_295 - scratch_190*scratch_296 - scratch_191*scratch_297 - scratch_192*scratch_298 - scratch_193*scratch_299 - scratch_194*scratch_300;
    amrex::Real scratch_302 = -A[348]*scratch_210 + A[362] - scratch_211*scratch_289 - scratch_212*scratch_290 - scratch_213*scratch_291 - scratch_214*scratch_292 - scratch_215*scratch_293 - scratch_216*scratch_294 - scratch_217*scratch_295 - scratch_218*scratch_296 - scratch_219*scratch_297 - scratch_220*scratch_298 - scratch_221*scratch_299 - scratch_222*scratch_300 - scratch_223*scratch_301;
    amrex::Real scratch_303 = -A[348]*scratch_240 + A[363] - scratch_241*scratch_289 - scratch_242*scratch_290 - scratch_243*scratch_291 - scratch_244*scratch_292 - scratch_245*scratch_293 - scratch_246*scratch_294 - scratch_247*scratch_295 - scratch_248*scratch_296 - scratch_249*scratch_297 - scratch_250*scratch_298 - scratch_251*scratch_299 - scratch_252*scratch_300 - scratch_253*scratch_301 - scratch_254*scratch_302;
    amrex::Real scratch_304 = 1.0/(-A[348]*scratch_272 + A[364] - scratch_273*scratch_289 - scratch_274*scratch_290 - scratch_275*scratch_291 - scratch_276*scratch_292 - scratch_277*scratch_293 - scratch_278*scratch_294 - scratch_279*scratch_295 - scratch_280*scratch_296 - scratch_281*scratch_297 - scratch_282*scratch_298 - scratch_283*scratch_299 - scratch_284*scratch_300 - scratch_285*scratch_301 - scratch_286*scratch_302 - scratch_287*scratch_303);
    amrex::Real scratch_305 = scratch_304*(-A[348]*scratch_2 + A[367] - scratch_10*scratch_289 - scratch_109*scratch_297 - scratch_11*scratch_290 - scratch_131*scratch_298 - scratch_155*scratch_299 - scratch_181*scratch_300 - scratch_19*scratch_291 - scratch_209*scratch_301 - scratch_239*scratch_302 - scratch_271*scratch_303 - scratch_29*scratch_292 - scratch_293*scratch_41 - scratch_294*scratch_55 - scratch_295*scratch_71 - scratch_296*scratch_89);
    amrex::Real scratch_306 = A[17]*scratch_0;
    amrex::Real scratch_307 = scratch_5*(-A[22]*scratch_306 + A[39]);
    amrex::Real scratch_308 = scratch_9*(-A[44]*scratch_306 + A[61] - scratch_307*scratch_8);
    amrex::Real scratch_309 = scratch_18*(-A[66]*scratch_306 + A[83] - scratch_16*scratch_307 - scratch_17*scratch_308);
    amrex::Real scratch_310 = scratch_28*(A[105] - A[88]*scratch_306 - scratch_25*scratch_307 - scratch_26*scratch_308 - scratch_27*scratch_309);
    amrex::Real scratch_311 = scratch_40*(-A[110]*scratch_306 + A[127] - scratch_307*scratch_36 - scratch_308*scratch_37 - scratch_309*scratch_38 - scratch_310*scratch_39);
    amrex::Real scratch_312 = scratch_54*(-A[132]*scratch_306 + A[149] - scratch_307*scratch_49 - scratch_308*scratch_50 - scratch_309*scratch_51 - scratch_310*scratch_52 - scratch_311*scratch_53);
    amrex::Real scratch_313 = scratch_70*(-A[154]*scratch_306 + A[171] - scratch_307*scratch_64 - scratch_308*scratch_65 - scratch_309*scratch_66 - scratch_310*scratch_67 - scratch_311*scratch_68 - scratch_312*scratch_69);
    amrex::Real scratch_314 = scratch_88*(-A[176]*scratch_306 - scratch_307*scratch_81 - scratch_308*scratch_82 - scratch_309*scratch_83 - scratch_310*scratch_84 - scratch_311*scratch_85 - scratch_312*scratch_86 - scratch_313*scratch_87);
    amrex::Real scratch_315 = scratch_108*(-A[194]*scratch_306 + A[211] - scratch_100*scratch_307 - scratch_101*scratch_308 - scratch_102*scratch_309 - scratch_103*scratch_310 - scratch_104*scratch_311 - scratch_105*scratch_312 - scratch_106*scratch_313 - scratch_107*scratch_314);
    amrex::Real scratch_316 = scratch_130*(-A[216]*scratch_306 + A[233] - scratch_121*scratch_307 - scratch_122*scratch_308 - scratch_123*scratch_309 - scratch_124*scratch_310 - scratch_125*scratch_311 - scratch_126*scratch_312 - scratch_127*scratch_313 - scratch_128*scratch_314 - scratch_129*scratch_315);
    amrex::Real scratch_317 = scratch_154*(-A[238]*scratch_306 + A[255] - scratch_144*scratch_307 - scratch_145*scratch_308 - scratch_146*scratch_309 - scratch_147*scratch_310 - scratch_148*scratch_311 - scratch_149*scratch_312 - scratch_150*scratch_313 - scratch_151*scratch_314 - scratch_152*scratch_315 - scratch_153*scratch_316);
    amrex::Real scratch_318 = scratch_180*(-A[260]*scratch_306 + A[277] - scratch_169*scratch_307 - scratch_170*scratch_308 - scratch_171*scratch_309 - scratch_172*scratch_310 - scratch_173*scratch_311 - scratch_174*scratch_312 - scratch_175*scratch_313 - scratch_176*scratch_314 - scratch_177*scratch_315 - scratch_178*scratch_316 - scratch_179*scratch_317);
    amrex::Real scratch_319 = scratch_208*(-A[282]*scratch_306 + A[299] - scratch_196*scratch_307 - scratch_197*scratch_308 - scratch_198*scratch_309 - scratch_199*scratch_310 - scratch_200*scratch_311 - scratch_201*scratch_312 - scratch_202*scratch_313 - scratch_203*scratch_314 - scratch_204*scratch_315 - scratch_205*scratch_316 - scratch_206*scratch_317 - scratch_207*scratch_318);
    amrex::Real scratch_320 = scratch_238*(-A[304]*scratch_306 + A[321] - scratch_225*scratch_307 - scratch_226*scratch_308 - scratch_227*scratch_309 - scratch_228*scratch_310 - scratch_229*scratch_311 - scratch_230*scratch_312 - scratch_231*scratch_313 - scratch_232*scratch_314 - scratch_233*scratch_315 - scratch_234*scratch_316 - scratch_235*scratch_317 - scratch_236*scratch_318 - scratch_237*scratch_319);
    amrex::Real scratch_321 = scratch_270*(-A[326]*scratch_306 + A[343] - scratch_256*scratch_307 - scratch_257*scratch_308 - scratch_258*scratch_309 - scratch_259*scratch_310 - scratch_260*scratch_311 - scratch_261*scratch_312 - scratch_262*scratch_313 - scratch_263*scratch_314 - scratch_264*scratch_315 - scratch_265*scratch_316 - scratch_266*scratch_317 - scratch_267*scratch_318 - scratch_268*scratch_319 - scratch_269*scratch_320);
    amrex::Real scratch_322 = scratch_304*(-A[348]*scratch_306 + A[365] - scratch_289*scratch_307 - scratch_290*scratch_308 - scratch_291*scratch_309 - scratch_292*scratch_310 - scratch_293*scratch_311 - scratch_294*scratch_312 - scratch_295*scratch_313 - scratch_296*scratch_314 - scratch_297*scratch_315 - scratch_298*scratch_316 - scratch_299*scratch_317 - scratch_300*scratch_318 - scratch_301*scratch_319 - scratch_302*scratch_320 - scratch_303*scratch_321);
    amrex::Real scratch_323 = -scratch_120*scratch_316 - scratch_143*scratch_317 - scratch_15*scratch_309 - scratch_168*scratch_318 - scratch_195*scratch_319 - scratch_224*scratch_320 - scratch_24*scratch_310 - scratch_255*scratch_321 - scratch_288*scratch_322 + scratch_306 - scratch_307*scratch_4 - scratch_308*scratch_7 - scratch_311*scratch_35 - scratch_312*scratch_48 - scratch_313*scratch_63 - scratch_314*scratch_80 - scratch_315*scratch_99;
    amrex::Real scratch_324 = -A[370]*scratch_4 + A[371];
    amrex::Real scratch_325 = -A[370]*scratch_3 + A[372] - scratch_324*scratch_6;
    amrex::Real scratch_326 = -A[370]*scratch_12 + A[373] - scratch_13*scratch_324 - scratch_14*scratch_325;
    amrex::Real scratch_327 = -A[370]*scratch_20 + A[374] - scratch_21*scratch_324 - scratch_22*scratch_325 - scratch_23*scratch_326;
    amrex::Real scratch_328 = -A[370]*scratch_30 + A[375] - scratch_31*scratch_324 - scratch_32*scratch_325 - scratch_326*scratch_33 - scratch_327*scratch_34;
    amrex::Real scratch_329 = -A[370]*scratch_42 + A[376] - scratch_324*scratch_43 - scratch_325*scratch_44 - scratch_326*scratch_45 - scratch_327*scratch_46 - scratch_328*scratch_47;
    amrex::Real scratch_330 = -A[370]*scratch_56 + A[377] - scratch_324*scratch_57 - scratch_325*scratch_58 - scratch_326*scratch_59 - scratch_327*scratch_60 - scratch_328*scratch_61 - scratch_329*scratch_62;
    amrex::Real scratch_331 = -A[370]*scratch_72 + A[378] - scratch_324*scratch_73 - scratch_325*scratch_74 - scratch_326*scratch_75 - scratch_327*scratch_76 - scratch_328*scratch_77 - scratch_329*scratch_78 - scratch_330*scratch_79;
    amrex::Real scratch_332 = -A[370]*scratch_90 + A[379] - scratch_324*scratch_91 - scratch_325*scratch_92 - scratch_326*scratch_93 - scratch_327*scratch_94 - scratch_328*scratch_95 - scratch_329*scratch_96 - scratch_330*scratch_97 - scratch_331*scratch_98;
    amrex::Real scratch_333 = -A[370]*scratch_110 + A[380] - scratch_111*scratch_324 - scratch_112*scratch_325 - scratch_113*scratch_326 - scratch_114*scratch_327 - scratch_115*scratch_328 - scratch_116*scratch_329 - scratch_117*scratch_330 - scratch_118*scratch_331 - scratch_119*scratch_332;
    amrex::Real scratch_334 = -A[370]*scratch_132 + A[381] - scratch_133*scratch_324 - scratch_134*scratch_325 - scratch_135*scratch_326 - scratch_136*scratch_327 - scratch_137*scratch_328 - scratch_138*scratch_329 - scratch_139*scratch_330 - scratch_140*scratch_331 - scratch_141*scratch_332 - scratch_142*scratch_333;
    amrex::Real scratch_335 = -A[370]*scratch_156 + A[382] - scratch_157*scratch_324 - scratch_158*scratch_325 - scratch_159*scratch_326 - scratch_160*scratch_327 - scratch_161*scratch_328 - scratch_162*scratch_329 - scratch_163*scratch_330 - scratch_164*scratch_331 - scratch_165*scratch_332 - scratch_166*scratch_333 - scratch_167*scratch_334;
    amrex::Real scratch_336 = -A[370]*scratch_182 + A[383] - scratch_183*scratch_324 - scratch_184*scratch_325 - scratch_185*scratch_326 - scratch_186*scratch_327 - scratch_187*scratch_328 - scratch_188*scratch_329 - scratch_189*scratch_330 - scratch_190*scratch_331 - scratch_191*scratch_332 - scratch_192*scratch_333 - scratch_193*scratch_334 - scratch_194*scratch_335;
    amrex::Real scratch_337 = -A[370]*scratch_210 + A[384] - scratch_211*scratch_324 - scratch_212*scratch_325 - scratch_213*scratch_326 - scratch_214*scratch_327 - scratch_215*scratch_328 - scratch_216*scratch_329 - scratch_217*scratch_330 - scratch_218*scratch_331 - scratch_219*scratch_332 - scratch_220*scratch_333 - scratch_221*scratch_334 - scratch_222*scratch_335 - scratch_223*scratch_336;
    amrex::Real scratch_338 = -A[370]*scratch_240 + A[385] - scratch_241*scratch_324 - scratch_242*scratch_325 - scratch_243*scratch_326 - scratch_244*scratch_327 - scratch_245*scratch_328 - scratch_246*scratch_329 - scratch_247*scratch_330 - scratch_248*scratch_331 - scratch_249*scratch_332 - scratch_250*scratch_333 - scratch_251*scratch_334 - scratch_252*scratch_335 - scratch_253*scratch_336 - scratch_254*scratch_337;
    amrex::Real scratch_339 = -A[370]*scratch_272 + A[386] - scratch_273*scratch_324 - scratch_274*scratch_325 - scratch_275*scratch_326 - scratch_276*scratch_327 - scratch_277*scratch_328 - scratch_278*scratch_329 - scratch_279*scratch_330 - scratch_280*scratch_331 - scratch_281*scratch_332 - scratch_282*scratch_333 - scratch_283*scratch_334 - scratch_284*scratch_335 - scratch_285*scratch_336 - scratch_286*scratch_337 - scratch_287*scratch_338;
    amrex::Real scratch_340 = 1.0/(-A[370]*scratch_306 + A[387] - scratch_307*scratch_324 - scratch_308*scratch_325 - scratch_309*scratch_326 - scratch_310*scratch_327 - scratch_311*scratch_328 - scratch_312*scratch_329 - scratch_313*scratch_330 - scratch_314*scratch_331 - scratch_315*scratch_332 - scratch_316*scratch_333 - scratch_317*scratch_334 - scratch_318*scratch_335 - scratch_319*scratch_336 - scratch_320*scratch_337 - scratch_321*scratch_338 - scratch_322*scratch_339);
    amrex::Real scratch_341 = scratch_340*(-A[370]*scratch_2 + A[389] - scratch_10*scratch_324 - scratch_109*scratch_332 - scratch_11*scratch_325 - scratch_131*scratch_333 - scratch_155*scratch_334 - scratch_181*scratch_335 - scratch_19*scratch_326 - scratch_209*scratch_336 - scratch_239*scratch_337 - scratch_271*scratch_338 - scratch_29*scratch_327 - scratch_305*scratch_339 - scratch_328*scratch_41 - scratch_329*scratch_55 - scratch_330*scratch_71 - scratch_331*scratch_89);
    amrex::Real scratch_342 = -A[392]*scratch_4 + A[393];
    amrex::Real scratch_343 = -A[392]*scratch_3 + A[394] - scratch_342*scratch_6;
    amrex::Real scratch_344 = -A[392]*scratch_12 + A[395] - scratch_13*scratch_342 - scratch_14*scratch_343;
    amrex::Real scratch_345 = -A[392]*scratch_20 + A[396] - scratch_21*scratch_342 - scratch_22*scratch_343 - scratch_23*scratch_344;
    amrex::Real scratch_346 = -A[392]*scratch_30 + A[397] - scratch_31*scratch_342 - scratch_32*scratch_343 - scratch_33*scratch_344 - scratch_34*scratch_345;
    amrex::Real scratch_347 = -A[392]*scratch_42 + A[398] - scratch_342*scratch_43 - scratch_343*scratch_44 - scratch_344*scratch_45 - scratch_345*scratch_46 - scratch_346*scratch_47;
    amrex::Real scratch_348 = -A[392]*scratch_56 + A[399] - scratch_342*scratch_57 - scratch_343*scratch_58 - scratch_344*scratch_59 - scratch_345*scratch_60 - scratch_346*scratch_61 - scratch_347*scratch_62;
    amrex::Real scratch_349 = -A[392]*scratch_72 + A[400] - scratch_342*scratch_73 - scratch_343*scratch_74 - scratch_344*scratch_75 - scratch_345*scratch_76 - scratch_346*scratch_77 - scratch_347*scratch_78 - scratch_348*scratch_79;
    amrex::Real scratch_350 = -A[392]*scratch_90 + A[401] - scratch_342*scratch_91 - scratch_343*scratch_92 - scratch_344*scratch_93 - scratch_345*scratch_94 - scratch_346*scratch_95 - scratch_347*scratch_96 - scratch_348*scratch_97 - scratch_349*scratch_98;
    amrex::Real scratch_351 = -A[392]*scratch_110 + A[402] - scratch_111*scratch_342 - scratch_112*scratch_343 - scratch_113*scratch_344 - scratch_114*scratch_345 - scratch_115*scratch_346 - scratch_116*scratch_347 - scratch_117*scratch_348 - scratch_118*scratch_349 - scratch_119*scratch_350;
    amrex::Real scratch_352 = -A[392]*scratch_132 + A[403] - scratch_133*scratch_342 - scratch_134*scratch_343 - scratch_135*scratch_344 - scratch_136*scratch_345 - scratch_137*scratch_346 - scratch_138*scratch_347 - scratch_139*scratch_348 - scratch_140*scratch_349 - scratch_141*scratch_350 - scratch_142*scratch_351;
    amrex::Real scratch_353 = -A[392]*scratch_156 + A[404] - scratch_157*scratch_342 - scratch_158*scratch_343 - scratch_159*scratch_344 - scratch_160*scratch_345 - scratch_161*scratch_346 - scratch_162*scratch_347 - scratch_163*scratch_348 - scratch_164*scratch_349 - scratch_165*scratch_350 - scratch_166*scratch_351 - scratch_167*scratch_352;
    amrex::Real scratch_354 = -A[392]*scratch_182 + A[405] - scratch_183*scratch_342 - scratch_184*scratch_343 - scratch_185*scratch_344 - scratch_186*scratch_345 - scratch_187*scratch_346 - scratch_188*scratch_347 - scratch_189*scratch_348 - scratch_190*scratch_349 - scratch_191*scratch_350 - scratch_192*scratch_351 - scratch_193*scratch_352 - scratch_194*scratch_353;
    amrex::Real scratch_355 = -A[392]*scratch_210 + A[406] - scratch_211*scratch_342 - scratch_212*scratch_343 - scratch_213*scratch_344 - scratch_214*scratch_345 - scratch_215*scratch_346 - scratch_216*scratch_347 - scratch_217*scratch_348 - scratch_218*scratch_349 - scratch_219*scratch_350 - scratch_220*scratch_351 - scratch_221*scratch_352 - scratch_222*scratch_353 - scratch_223*scratch_354;
    amrex::Real scratch_356 = -A[392]*scratch_240 + A[407] - scratch_241*scratch_342 - scratch_242*scratch_343 - scratch_243*scratch_344 - scratch_244*scratch_345 - scratch_245*scratch_346 - scratch_246*scratch_347 - scratch_247*scratch_348 - scratch_248*scratch_349 - scratch_249*scratch_350 - scratch_250*scratch_351 - scratch_251*scratch_352 - scratch_252*scratch_353 - scratch_253*scratch_354 - scratch_254*scratch_355;
    amrex::Real scratch_357 = -A[392]*scratch_272 + A[408] - scratch_273*scratch_342 - scratch_274*scratch_343 - scratch_275*scratch_344 - scratch_276*scratch_345 - scratch_277*scratch_346 - scratch_278*scratch_347 - scratch_279*scratch_348 - scratch_280*scratch_349 - scratch_281*scratch_350 - scratch_282*scratch_351 - scratch_283*scratch_352 - scratch_284*scratch_353 - scratch_285*scratch_354 - scratch_286*scratch_355 - scratch_287*scratch_356;
    amrex::Real scratch_358 = -A[392]*scratch_306 + A[409] - scratch_307*scratch_342 - scratch_308*scratch_343 - scratch_309*scratch_344 - scratch_310*scratch_345 - scratch_311*scratch_346 - scratch_312*scratch_347 - scratch_313*scratch_348 - scratch_314*scratch_349 - scratch_315*scratch_350 - scratch_316*scratch_351 - scratch_317*scratch_352 - scratch_318*scratch_353 - scratch_319*scratch_354 - scratch_320*scratch_355 - scratch_321*scratch_356 - scratch_322*scratch_357;
    amrex::Real scratch_359 = -A[392]*scratch_2 + A[411] - scratch_10*scratch_342 - scratch_109*scratch_350 - scratch_11*scratch_343 - scratch_131*scratch_351 - scratch_155*scratch_352 - scratch_181*scratch_353 - scratch_19*scratch_344 - scratch_209*scratch_354 - scratch_239*scratch_355 - scratch_271*scratch_356 - scratch_29*scratch_345 - scratch_305*scratch_357 - scratch_341*scratch_358 - scratch_346*scratch_41 - scratch_347*scratch_55 - scratch_348*scratch_71 - scratch_349*scratch_89;
    amrex::Real scratch_360 = A[18]*scratch_0;
    amrex::Real scratch_361 = scratch_5*(-A[22]*scratch_360 + A[40]);
    amrex::Real scratch_362 = scratch_9*(-A[44]*scratch_360 + A[62] - scratch_361*scratch_8);
    amrex::Real scratch_363 = scratch_18*(-A[66]*scratch_360 + A[84] - scratch_16*scratch_361 - scratch_17*scratch_362);
    amrex::Real scratch_364 = scratch_28*(A[106] - A[88]*scratch_360 - scratch_25*scratch_361 - scratch_26*scratch_362 - scratch_27*scratch_363);
    amrex::Real scratch_365 = scratch_40*(-A[110]*scratch_360 + A[128] - scratch_36*scratch_361 - scratch_362*scratch_37 - scratch_363*scratch_38 - scratch_364*scratch_39);
    amrex::Real scratch_366 = scratch_54*(-A[132]*scratch_360 + A[150] - scratch_361*scratch_49 - scratch_362*scratch_50 - scratch_363*scratch_51 - scratch_364*scratch_52 - scratch_365*scratch_53);
    amrex::Real scratch_367 = scratch_70*(-A[154]*scratch_360 + A[172] - scratch_361*scratch_64 - scratch_362*scratch_65 - scratch_363*scratch_66 - scratch_364*scratch_67 - scratch_365*scratch_68 - scratch_366*scratch_69);
    amrex::Real scratch_368 = scratch_88*(-A[176]*scratch_360 - scratch_361*scratch_81 - scratch_362*scratch_82 - scratch_363*scratch_83 - scratch_364*scratch_84 - scratch_365*scratch_85 - scratch_366*scratch_86 - scratch_367*scratch_87);
    amrex::Real scratch_369 = scratch_108*(-A[194]*scratch_360 + A[212] - scratch_100*scratch_361 - scratch_101*scratch_362 - scratch_102*scratch_363 - scratch_103*scratch_364 - scratch_104*scratch_365 - scratch_105*scratch_366 - scratch_106*scratch_367 - scratch_107*scratch_368);
    amrex::Real scratch_370 = scratch_130*(-A[216]*scratch_360 + A[234] - scratch_121*scratch_361 - scratch_122*scratch_362 - scratch_123*scratch_363 - scratch_124*scratch_364 - scratch_125*scratch_365 - scratch_126*scratch_366 - scratch_127*scratch_367 - scratch_128*scratch_368 - scratch_129*scratch_369);
    amrex::Real scratch_371 = scratch_154*(-A[238]*scratch_360 + A[256] - scratch_144*scratch_361 - scratch_145*scratch_362 - scratch_146*scratch_363 - scratch_147*scratch_364 - scratch_148*scratch_365 - scratch_149*scratch_366 - scratch_150*scratch_367 - scratch_151*scratch_368 - scratch_152*scratch_369 - scratch_153*scratch_370);
    amrex::Real scratch_372 = scratch_180*(-A[260]*scratch_360 + A[278] - scratch_169*scratch_361 - scratch_170*scratch_362 - scratch_171*scratch_363 - scratch_172*scratch_364 - scratch_173*scratch_365 - scratch_174*scratch_366 - scratch_175*scratch_367 - scratch_176*scratch_368 - scratch_177*scratch_369 - scratch_178*scratch_370 - scratch_179*scratch_371);
    amrex::Real scratch_373 = scratch_208*(-A[282]*scratch_360 + A[300] - scratch_196*scratch_361 - scratch_197*scratch_362 - scratch_198*scratch_363 - scratch_199*scratch_364 - scratch_200*scratch_365 - scratch_201*scratch_366 - scratch_202*scratch_367 - scratch_203*scratch_368 - scratch_204*scratch_369 - scratch_205*scratch_370 - scratch_206*scratch_371 - scratch_207*scratch_372);
    amrex::Real scratch_374 = scratch_238*(-A[304]*scratch_360 + A[322] - scratch_225*scratch_361 - scratch_226*scratch_362 - scratch_227*scratch_363 - scratch_228*scratch_364 - scratch_229*scratch_365 - scratch_230*scratch_366 - scratch_231*scratch_367 - scratch_232*scratch_368 - scratch_233*scratch_369 - scratch_234*scratch_370 - scratch_235*scratch_371 - scratch_236*scratch_372 - scratch_237*scratch_373);
    amrex::Real scratch_375 = scratch_270*(-A[326]*scratch_360 + A[344] - scratch_256*scratch_361 - scratch_257*scratch_362 - scratch_258*scratch_363 - scratch_259*scratch_364 - scratch_260*scratch_365 - scratch_261*scratch_366 - scratch_262*scratch_367 - scratch_263*scratch_368 - scratch_264*scratch_369 - scratch_265*scratch_370 - scratch_266*scratch_371 - scratch_267*scratch_372 - scratch_268*scratch_373 - scratch_269*scratch_374);
    amrex::Real scratch_376 = scratch_304*(-A[348]*scratch_360 + A[366] - scratch_289*scratch_361 - scratch_290*scratch_362 - scratch_291*scratch_363 - scratch_292*scratch_364 - scratch_293*scratch_365 - scratch_294*scratch_366 - scratch_295*scratch_367 - scratch_296*scratch_368 - scratch_297*scratch_369 - scratch_298*scratch_370 - scratch_299*scratch_371 - scratch_300*scratch_372 - scratch_301*scratch_373 - scratch_302*scratch_374 - scratch_303*scratch_375);
    amrex::Real scratch_377 = scratch_340*(-A[370]*scratch_360 + A[388] - scratch_324*scratch_361 - scratch_325*scratch_362 - scratch_326*scratch_363 - scratch_327*scratch_364 - scratch_328*scratch_365 - scratch_329*scratch_366 - scratch_330*scratch_367 - scratch_331*scratch_368 - scratch_332*scratch_369 - scratch_333*scratch_370 - scratch_334*scratch_371 - scratch_335*scratch_372 - scratch_336*scratch_373 - scratch_337*scratch_374 - scratch_338*scratch_375 - scratch_339*scratch_376);
    amrex::Real scratch_378 = 1.0/(-A[392]*scratch_360 + A[410] - scratch_342*scratch_361 - scratch_343*scratch_362 - scratch_344*scratch_363 - scratch_345*scratch_364 - scratch_346*scratch_365 - scratch_347*scratch_366 - scratch_348*scratch_367 - scratch_349*scratch_368 - scratch_350*scratch_369 - scratch_351*scratch_370 - scratch_352*scratch_371 - scratch_353*scratch_372 - scratch_354*scratch_373 - scratch_355*scratch_374 - scratch_356*scratch_375 - scratch_357*scratch_376 - scratch_358*scratch_377);
    amrex::Real scratch_379 = -scratch_120*scratch_370 - scratch_143*scratch_371 - scratch_15*scratch_363 - scratch_168*scratch_372 - scratch_195*scratch_373 - scratch_224*scratch_374 - scratch_24*scratch_364 - scratch_255*scratch_375 - scratch_288*scratch_376 - scratch_323*scratch_377 - scratch_35*scratch_365 + scratch_360 - scratch_361*scratch_4 - scratch_362*scratch_7 - scratch_366*scratch_48 - scratch_367*scratch_63 - scratch_368*scratch_80 - scratch_369*scratch_99;
    amrex::Real scratch_380 = scratch_378*scratch_379;
    amrex::Real scratch_381 = b[19]/A[414];
    amrex::Real scratch_382 = A[20]*scratch_0;
    amrex::Real scratch_383 = scratch_5*(-A[22]*scratch_382 + A[42]);
    amrex::Real scratch_384 = scratch_9*(-A[44]*scratch_382 + A[64] - scratch_383*scratch_8);
    amrex::Real scratch_385 = scratch_18*(-A[66]*scratch_382 + A[86] - scratch_16*scratch_383 - scratch_17*scratch_384);
    amrex::Real scratch_386 = scratch_28*(A[108] - A[88]*scratch_382 - scratch_25*scratch_383 - scratch_26*scratch_384 - scratch_27*scratch_385);
    amrex::Real scratch_387 = scratch_40*(-A[110]*scratch_382 + A[130] - scratch_36*scratch_383 - scratch_37*scratch_384 - scratch_38*scratch_385 - scratch_386*scratch_39);
    amrex::Real scratch_388 = scratch_54*(-A[132]*scratch_382 + A[152] - scratch_383*scratch_49 - scratch_384*scratch_50 - scratch_385*scratch_51 - scratch_386*scratch_52 - scratch_387*scratch_53);
    amrex::Real scratch_389 = scratch_70*(-A[154]*scratch_382 + A[174] - scratch_383*scratch_64 - scratch_384*scratch_65 - scratch_385*scratch_66 - scratch_386*scratch_67 - scratch_387*scratch_68 - scratch_388*scratch_69);
    amrex::Real scratch_390 = scratch_88*(-A[176]*scratch_382 + A[192] - scratch_383*scratch_81 - scratch_384*scratch_82 - scratch_385*scratch_83 - scratch_386*scratch_84 - scratch_387*scratch_85 - scratch_388*scratch_86 - scratch_389*scratch_87);
    amrex::Real scratch_391 = scratch_108*(-A[194]*scratch_382 + A[214] - scratch_100*scratch_383 - scratch_101*scratch_384 - scratch_102*scratch_385 - scratch_103*scratch_386 - scratch_104*scratch_387 - scratch_105*scratch_388 - scratch_106*scratch_389 - scratch_107*scratch_390);
    amrex::Real scratch_392 = scratch_130*(-A[216]*scratch_382 + A[236] - scratch_121*scratch_383 - scratch_122*scratch_384 - scratch_123*scratch_385 - scratch_124*scratch_386 - scratch_125*scratch_387 - scratch_126*scratch_388 - scratch_127*scratch_389 - scratch_128*scratch_390 - scratch_129*scratch_391);
    amrex::Real scratch_393 = scratch_154*(-A[238]*scratch_382 + A[258] - scratch_144*scratch_383 - scratch_145*scratch_384 - scratch_146*scratch_385 - scratch_147*scratch_386 - scratch_148*scratch_387 - scratch_149*scratch_388 - scratch_150*scratch_389 - scratch_151*scratch_390 - scratch_152*scratch_391 - scratch_153*scratch_392);
    amrex::Real scratch_394 = scratch_180*(-A[260]*scratch_382 + A[280] - scratch_169*scratch_383 - scratch_170*scratch_384 - scratch_171*scratch_385 - scratch_172*scratch_386 - scratch_173*scratch_387 - scratch_174*scratch_388 - scratch_175*scratch_389 - scratch_176*scratch_390 - scratch_177*scratch_391 - scratch_178*scratch_392 - scratch_179*scratch_393);
    amrex::Real scratch_395 = scratch_208*(-A[282]*scratch_382 + A[302] - scratch_196*scratch_383 - scratch_197*scratch_384 - scratch_198*scratch_385 - scratch_199*scratch_386 - scratch_200*scratch_387 - scratch_201*scratch_388 - scratch_202*scratch_389 - scratch_203*scratch_390 - scratch_204*scratch_391 - scratch_205*scratch_392 - scratch_206*scratch_393 - scratch_207*scratch_394);
    amrex::Real scratch_396 = scratch_238*(-A[304]*scratch_382 + A[324] - scratch_225*scratch_383 - scratch_226*scratch_384 - scratch_227*scratch_385 - scratch_228*scratch_386 - scratch_229*scratch_387 - scratch_230*scratch_388 - scratch_231*scratch_389 - scratch_232*scratch_390 - scratch_233*scratch_391 - scratch_234*scratch_392 - scratch_235*scratch_393 - scratch_236*scratch_394 - scratch_237*scratch_395);
    amrex::Real scratch_397 = scratch_270*(-A[326]*scratch_382 + A[346] - scratch_256*scratch_383 - scratch_257*scratch_384 - scratch_258*scratch_385 - scratch_259*scratch_386 - scratch_260*scratch_387 - scratch_261*scratch_388 - scratch_262*scratch_389 - scratch_263*scratch_390 - scratch_264*scratch_391 - scratch_265*scratch_392 - scratch_266*scratch_393 - scratch_267*scratch_394 - scratch_268*scratch_395 - scratch_269*scratch_396);
    amrex::Real scratch_398 = scratch_304*(-A[348]*scratch_382 + A[368] - scratch_289*scratch_383 - scratch_290*scratch_384 - scratch_291*scratch_385 - scratch_292*scratch_386 - scratch_293*scratch_387 - scratch_294*scratch_388 - scratch_295*scratch_389 - scratch_296*scratch_390 - scratch_297*scratch_391 - scratch_298*scratch_392 - scratch_299*scratch_393 - scratch_300*scratch_394 - scratch_301*scratch_395 - scratch_302*scratch_396 - scratch_303*scratch_397);
    amrex::Real scratch_399 = scratch_340*(-A[370]*scratch_382 + A[390] - scratch_324*scratch_383 - scratch_325*scratch_384 - scratch_326*scratch_385 - scratch_327*scratch_386 - scratch_328*scratch_387 - scratch_329*scratch_388 - scratch_330*scratch_389 - scratch_331*scratch_390 - scratch_332*scratch_391 - scratch_333*scratch_392 - scratch_334*scratch_393 - scratch_335*scratch_394 - scratch_336*scratch_395 - scratch_337*scratch_396 - scratch_338*scratch_397 - scratch_339*scratch_398);
    amrex::Real scratch_400 = -A[392]*scratch_382 + A[412] - scratch_342*scratch_383 - scratch_343*scratch_384 - scratch_344*scratch_385 - scratch_345*scratch_386 - scratch_346*scratch_387 - scratch_347*scratch_388 - scratch_348*scratch_389 - scratch_349*scratch_390 - scratch_350*scratch_391 - scratch_351*scratch_392 - scratch_352*scratch_393 - scratch_353*scratch_394 - scratch_354*scratch_395 - scratch_355*scratch_396 - scratch_356*scratch_397 - scratch_357*scratch_398 - scratch_358*scratch_399;
    amrex::Real scratch_401 = b[20]/A[415];
    amrex::Real scratch_402 = scratch_5*(-A[22]*scratch_1 + b[1]);
    amrex::Real scratch_403 = scratch_9*(-A[44]*scratch_1 + b[2] - scratch_402*scratch_8);
    amrex::Real scratch_404 = scratch_18*(-A[66]*scratch_1 + b[3] - scratch_16*scratch_402 - scratch_17*scratch_403);
    amrex::Real scratch_405 = scratch_28*(-A[88]*scratch_1 + b[4] - scratch_25*scratch_402 - scratch_26*scratch_403 - scratch_27*scratch_404);
    amrex::Real scratch_406 = scratch_40*(-A[110]*scratch_1 + b[5] - scratch_36*scratch_402 - scratch_37*scratch_403 - scratch_38*scratch_404 - scratch_39*scratch_405);
    amrex::Real scratch_407 = scratch_54*(-A[132]*scratch_1 + b[6] - scratch_402*scratch_49 - scratch_403*scratch_50 - scratch_404*scratch_51 - scratch_405*scratch_52 - scratch_406*scratch_53);
    amrex::Real scratch_408 = scratch_70*(-A[154]*scratch_1 + b[7] - scratch_402*scratch_64 - scratch_403*scratch_65 - scratch_404*scratch_66 - scratch_405*scratch_67 - scratch_406*scratch_68 - scratch_407*scratch_69);
    amrex::Real scratch_409 = scratch_88*(-A[176]*scratch_1 + b[8] - scratch_402*scratch_81 - scratch_403*scratch_82 - scratch_404*scratch_83 - scratch_405*scratch_84 - scratch_406*scratch_85 - scratch_407*scratch_86 - scratch_408*scratch_87);
    amrex::Real scratch_410 = scratch_108*(-A[194]*scratch_1 + b[9] - scratch_100*scratch_402 - scratch_101*scratch_403 - scratch_102*scratch_404 - scratch_103*scratch_405 - scratch_104*scratch_406 - scratch_105*scratch_407 - scratch_106*scratch_408 - scratch_107*scratch_409);
    amrex::Real scratch_411 = scratch_130*(-A[216]*scratch_1 + b[10] - scratch_121*scratch_402 - scratch_122*scratch_403 - scratch_123*scratch_404 - scratch_124*scratch_405 - scratch_125*scratch_406 - scratch_126*scratch_407 - scratch_127*scratch_408 - scratch_128*scratch_409 - scratch_129*scratch_410);
    amrex::Real scratch_412 = scratch_154*(-A[238]*scratch_1 + b[11] - scratch_144*scratch_402 - scratch_145*scratch_403 - scratch_146*scratch_404 - scratch_147*scratch_405 - scratch_148*scratch_406 - scratch_149*scratch_407 - scratch_150*scratch_408 - scratch_151*scratch_409 - scratch_152*scratch_410 - scratch_153*scratch_411);
    amrex::Real scratch_413 = scratch_180*(-A[260]*scratch_1 + b[12] - scratch_169*scratch_402 - scratch_170*scratch_403 - scratch_171*scratch_404 - scratch_172*scratch_405 - scratch_173*scratch_406 - scratch_174*scratch_407 - scratch_175*scratch_408 - scratch_176*scratch_409 - scratch_177*scratch_410 - scratch_178*scratch_411 - scratch_179*scratch_412);
    amrex::Real scratch_414 = scratch_208*(-A[282]*scratch_1 + b[13] - scratch_196*scratch_402 - scratch_197*scratch_403 - scratch_198*scratch_404 - scratch_199*scratch_405 - scratch_200*scratch_406 - scratch_201*scratch_407 - scratch_202*scratch_408 - scratch_203*scratch_409 - scratch_204*scratch_410 - scratch_205*scratch_411 - scratch_206*scratch_412 - scratch_207*scratch_413);
    amrex::Real scratch_415 = scratch_238*(-A[304]*scratch_1 + b[14] - scratch_225*scratch_402 - scratch_226*scratch_403 - scratch_227*scratch_404 - scratch_228*scratch_405 - scratch_229*scratch_406 - scratch_230*scratch_407 - scratch_231*scratch_408 - scratch_232*scratch_409 - scratch_233*scratch_410 - scratch_234*scratch_411 - scratch_235*scratch_412 - scratch_236*scratch_413 - scratch_237*scratch_414);
    amrex::Real scratch_416 = scratch_270*(-A[326]*scratch_1 + b[15] - scratch_256*scratch_402 - scratch_257*scratch_403 - scratch_258*scratch_404 - scratch_259*scratch_405 - scratch_260*scratch_406 - scratch_261*scratch_407 - scratch_262*scratch_408 - scratch_263*scratch_409 - scratch_264*scratch_410 - scratch_265*scratch_411 - scratch_266*scratch_412 - scratch_267*scratch_413 - scratch_268*scratch_414 - scratch_269*scratch_415);
    amrex::Real scratch_417 = scratch_304*(-A[348]*scratch_1 + b[16] - scratch_289*scratch_402 - scratch_290*scratch_403 - scratch_291*scratch_404 - scratch_292*scratch_405 - scratch_293*scratch_406 - scratch_294*scratch_407 - scratch_295*scratch_408 - scratch_296*scratch_409 - scratch_297*scratch_410 - scratch_298*scratch_411 - scratch_299*scratch_412 - scratch_300*scratch_413 - scratch_301*scratch_414 - scratch_302*scratch_415 - scratch_303*scratch_416);
    amrex::Real scratch_418 = scratch_340*(-A[370]*scratch_1 + b[17] - scratch_324*scratch_402 - scratch_325*scratch_403 - scratch_326*scratch_404 - scratch_327*scratch_405 - scratch_328*scratch_406 - scratch_329*scratch_407 - scratch_330*scratch_408 - scratch_331*scratch_409 - scratch_332*scratch_410 - scratch_333*scratch_411 - scratch_334*scratch_412 - scratch_335*scratch_413 - scratch_336*scratch_414 - scratch_337*scratch_415 - scratch_338*scratch_416 - scratch_339*scratch_417);
    amrex::Real scratch_419 = scratch_378*(-A[392]*scratch_1 + b[18] - scratch_342*scratch_402 - scratch_343*scratch_403 - scratch_344*scratch_404 - scratch_345*scratch_405 - scratch_346*scratch_406 - scratch_347*scratch_407 - scratch_348*scratch_408 - scratch_349*scratch_409 - scratch_350*scratch_410 - scratch_351*scratch_411 - scratch_352*scratch_412 - scratch_353*scratch_413 - scratch_354*scratch_414 - scratch_355*scratch_415 - scratch_356*scratch_416 - scratch_357*scratch_417 - scratch_358*scratch_418);
    amrex::Real scratch_420 = A[21]*scratch_0;
    amrex::Real scratch_421 = scratch_5*(-A[22]*scratch_420 + A[43]);
    amrex::Real scratch_422 = scratch_9*(-A[44]*scratch_420 + A[65] - scratch_421*scratch_8);
    amrex::Real scratch_423 = scratch_18*(-A[66]*scratch_420 + A[87] - scratch_16*scratch_421 - scratch_17*scratch_422);
    amrex::Real scratch_424 = scratch_28*(A[109] - A[88]*scratch_420 - scratch_25*scratch_421 - scratch_26*scratch_422 - scratch_27*scratch_423);
    amrex::Real scratch_425 = scratch_40*(-A[110]*scratch_420 + A[131] - scratch_36*scratch_421 - scratch_37*scratch_422 - scratch_38*scratch_423 - scratch_39*scratch_424);
    amrex::Real scratch_426 = scratch_54*(-A[132]*scratch_420 + A[153] - scratch_421*scratch_49 - scratch_422*scratch_50 - scratch_423*scratch_51 - scratch_424*scratch_52 - scratch_425*scratch_53);
    amrex::Real scratch_427 = scratch_70*(-A[154]*scratch_420 + A[175] - scratch_421*scratch_64 - scratch_422*scratch_65 - scratch_423*scratch_66 - scratch_424*scratch_67 - scratch_425*scratch_68 - scratch_426*scratch_69);
    amrex::Real scratch_428 = scratch_88*(-A[176]*scratch_420 + A[193] - scratch_421*scratch_81 - scratch_422*scratch_82 - scratch_423*scratch_83 - scratch_424*scratch_84 - scratch_425*scratch_85 - scratch_426*scratch_86 - scratch_427*scratch_87);
    amrex::Real scratch_429 = scratch_108*(-A[194]*scratch_420 + A[215] - scratch_100*scratch_421 - scratch_101*scratch_422 - scratch_102*scratch_423 - scratch_103*scratch_424 - scratch_104*scratch_425 - scratch_105*scratch_426 - scratch_106*scratch_427 - scratch_107*scratch_428);
    amrex::Real scratch_430 = scratch_130*(-A[216]*scratch_420 + A[237] - scratch_121*scratch_421 - scratch_122*scratch_422 - scratch_123*scratch_423 - scratch_124*scratch_424 - scratch_125*scratch_425 - scratch_126*scratch_426 - scratch_127*scratch_427 - scratch_128*scratch_428 - scratch_129*scratch_429);
    amrex::Real scratch_431 = scratch_154*(-A[238]*scratch_420 + A[259] - scratch_144*scratch_421 - scratch_145*scratch_422 - scratch_146*scratch_423 - scratch_147*scratch_424 - scratch_148*scratch_425 - scratch_149*scratch_426 - scratch_150*scratch_427 - scratch_151*scratch_428 - scratch_152*scratch_429 - scratch_153*scratch_430);
    amrex::Real scratch_432 = scratch_180*(-A[260]*scratch_420 + A[281] - scratch_169*scratch_421 - scratch_170*scratch_422 - scratch_171*scratch_423 - scratch_172*scratch_424 - scratch_173*scratch_425 - scratch_174*scratch_426 - scratch_175*scratch_427 - scratch_176*scratch_428 - scratch_177*scratch_429 - scratch_178*scratch_430 - scratch_179*scratch_431);
    amrex::Real scratch_433 = scratch_208*(-A[282]*scratch_420 + A[303] - scratch_196*scratch_421 - scratch_197*scratch_422 - scratch_198*scratch_423 - scratch_199*scratch_424 - scratch_200*scratch_425 - scratch_201*scratch_426 - scratch_202*scratch_427 - scratch_203*scratch_428 - scratch_204*scratch_429 - scratch_205*scratch_430 - scratch_206*scratch_431 - scratch_207*scratch_432);
    amrex::Real scratch_434 = scratch_238*(-A[304]*scratch_420 + A[325] - scratch_225*scratch_421 - scratch_226*scratch_422 - scratch_227*scratch_423 - scratch_228*scratch_424 - scratch_229*scratch_425 - scratch_230*scratch_426 - scratch_231*scratch_427 - scratch_232*scratch_428 - scratch_233*scratch_429 - scratch_234*scratch_430 - scratch_235*scratch_431 - scratch_236*scratch_432 - scratch_237*scratch_433);
    amrex::Real scratch_435 = scratch_270*(-A[326]*scratch_420 + A[347] - scratch_256*scratch_421 - scratch_257*scratch_422 - scratch_258*scratch_423 - scratch_259*scratch_424 - scratch_260*scratch_425 - scratch_261*scratch_426 - scratch_262*scratch_427 - scratch_263*scratch_428 - scratch_264*scratch_429 - scratch_265*scratch_430 - scratch_266*scratch_431 - scratch_267*scratch_432 - scratch_268*scratch_433 - scratch_269*scratch_434);
    amrex::Real scratch_436 = scratch_304*(-A[348]*scratch_420 + A[369] - scratch_289*scratch_421 - scratch_290*scratch_422 - scratch_291*scratch_423 - scratch_292*scratch_424 - scratch_293*scratch_425 - scratch_294*scratch_426 - scratch_295*scratch_427 - scratch_296*scratch_428 - scratch_297*scratch_429 - scratch_298*scratch_430 - scratch_299*scratch_431 - scratch_300*scratch_432 - scratch_301*scratch_433 - scratch_302*scratch_434 - scratch_303*scratch_435);
    amrex::Real scratch_437 = scratch_340*(-A[370]*scratch_420 + A[391] - scratch_324*scratch_421 - scratch_325*scratch_422 - scratch_326*scratch_423 - scratch_327*scratch_424 - scratch_328*scratch_425 - scratch_329*scratch_426 - scratch_330*scratch_427 - scratch_331*scratch_428 - scratch_332*scratch_429 - scratch_333*scratch_430 - scratch_334*scratch_431 - scratch_335*scratch_432 - scratch_336*scratch_433 - scratch_337*scratch_434 - scratch_338*scratch_435 - scratch_339*scratch_436);
    amrex::Real scratch_438 = -A[392]*scratch_420 + A[413] - scratch_342*scratch_421 - scratch_343*scratch_422 - scratch_344*scratch_423 - scratch_345*scratch_424 - scratch_346*scratch_425 - scratch_347*scratch_426 - scratch_348*scratch_427 - scratch_349*scratch_428 - scratch_350*scratch_429 - scratch_351*scratch_430 - scratch_352*scratch_431 - scratch_353*scratch_432 - scratch_354*scratch_433 - scratch_355*scratch_434 - scratch_356*scratch_435 - scratch_357*scratch_436 - scratch_358*scratch_437;
    amrex::Real scratch_439 = -A[416]*scratch_4 + A[417];
    amrex::Real scratch_440 = -A[416]*scratch_3 + A[418] - scratch_439*scratch_6;
    amrex::Real scratch_441 = -A[416]*scratch_12 + A[419] - scratch_13*scratch_439 - scratch_14*scratch_440;
    amrex::Real scratch_442 = -A[416]*scratch_20 + A[420] - scratch_21*scratch_439 - scratch_22*scratch_440 - scratch_23*scratch_441;
    amrex::Real scratch_443 = -A[416]*scratch_30 + A[421] - scratch_31*scratch_439 - scratch_32*scratch_440 - scratch_33*scratch_441 - scratch_34*scratch_442;
    amrex::Real scratch_444 = -A[416]*scratch_42 + A[422] - scratch_43*scratch_439 - scratch_44*scratch_440 - scratch_441*scratch_45 - scratch_442*scratch_46 - scratch_443*scratch_47;
    amrex::Real scratch_445 = -A[416]*scratch_56 + A[423] - scratch_439*scratch_57 - scratch_440*scratch_58 - scratch_441*scratch_59 - scratch_442*scratch_60 - scratch_443*scratch_61 - scratch_444*scratch_62;
    amrex::Real scratch_446 = -A[416]*scratch_72 + A[424] - scratch_439*scratch_73 - scratch_440*scratch_74 - scratch_441*scratch_75 - scratch_442*scratch_76 - scratch_443*scratch_77 - scratch_444*scratch_78 - scratch_445*scratch_79;
    amrex::Real scratch_447 = -A[416]*scratch_90 + A[425] - scratch_439*scratch_91 - scratch_440*scratch_92 - scratch_441*scratch_93 - scratch_442*scratch_94 - scratch_443*scratch_95 - scratch_444*scratch_96 - scratch_445*scratch_97 - scratch_446*scratch_98;
    amrex::Real scratch_448 = -A[416]*scratch_110 + A[426] - scratch_111*scratch_439 - scratch_112*scratch_440 - scratch_113*scratch_441 - scratch_114*scratch_442 - scratch_115*scratch_443 - scratch_116*scratch_444 - scratch_117*scratch_445 - scratch_118*scratch_446 - scratch_119*scratch_447;
    amrex::Real scratch_449 = -A[416]*scratch_132 + A[427] - scratch_133*scratch_439 - scratch_134*scratch_440 - scratch_135*scratch_441 - scratch_136*scratch_442 - scratch_137*scratch_443 - scratch_138*scratch_444 - scratch_139*scratch_445 - scratch_140*scratch_446 - scratch_141*scratch_447 - scratch_142*scratch_448;
    amrex::Real scratch_450 = -A[416]*scratch_156 + A[428] - scratch_157*scratch_439 - scratch_158*scratch_440 - scratch_159*scratch_441 - scratch_160*scratch_442 - scratch_161*scratch_443 - scratch_162*scratch_444 - scratch_163*scratch_445 - scratch_164*scratch_446 - scratch_165*scratch_447 - scratch_166*scratch_448 - scratch_167*scratch_449;
    amrex::Real scratch_451 = -A[416]*scratch_182 + A[429] - scratch_183*scratch_439 - scratch_184*scratch_440 - scratch_185*scratch_441 - scratch_186*scratch_442 - scratch_187*scratch_443 - scratch_188*scratch_444 - scratch_189*scratch_445 - scratch_190*scratch_446 - scratch_191*scratch_447 - scratch_192*scratch_448 - scratch_193*scratch_449 - scratch_194*scratch_450;
    amrex::Real scratch_452 = -A[416]*scratch_210 + A[430] - scratch_211*scratch_439 - scratch_212*scratch_440 - scratch_213*scratch_441 - scratch_214*scratch_442 - scratch_215*scratch_443 - scratch_216*scratch_444 - scratch_217*scratch_445 - scratch_218*scratch_446 - scratch_219*scratch_447 - scratch_220*scratch_448 - scratch_221*scratch_449 - scratch_222*scratch_450 - scratch_223*scratch_451;
    amrex::Real scratch_453 = -A[416]*scratch_240 + A[431] - scratch_241*scratch_439 - scratch_242*scratch_440 - scratch_243*scratch_441 - scratch_244*scratch_442 - scratch_245*scratch_443 - scratch_246*scratch_444 - scratch_247*scratch_445 - scratch_248*scratch_446 - scratch_249*scratch_447 - scratch_250*scratch_448 - scratch_251*scratch_449 - scratch_252*scratch_450 - scratch_253*scratch_451 - scratch_254*scratch_452;
    amrex::Real scratch_454 = -A[416]*scratch_272 + A[432] - scratch_273*scratch_439 - scratch_274*scratch_440 - scratch_275*scratch_441 - scratch_276*scratch_442 - scratch_277*scratch_443 - scratch_278*scratch_444 - scratch_279*scratch_445 - scratch_280*scratch_446 - scratch_281*scratch_447 - scratch_282*scratch_448 - scratch_283*scratch_449 - scratch_284*scratch_450 - scratch_285*scratch_451 - scratch_286*scratch_452 - scratch_287*scratch_453;
    amrex::Real scratch_455 = -A[416]*scratch_306 + A[433] - scratch_307*scratch_439 - scratch_308*scratch_440 - scratch_309*scratch_441 - scratch_310*scratch_442 - scratch_311*scratch_443 - scratch_312*scratch_444 - scratch_313*scratch_445 - scratch_314*scratch_446 - scratch_315*scratch_447 - scratch_316*scratch_448 - scratch_317*scratch_449 - scratch_318*scratch_450 - scratch_319*scratch_451 - scratch_320*scratch_452 - scratch_321*scratch_453 - scratch_322*scratch_454;
    amrex::Real scratch_456 = -A[416]*scratch_360 + A[434] - scratch_361*scratch_439 - scratch_362*scratch_440 - scratch_363*scratch_441 - scratch_364*scratch_442 - scratch_365*scratch_443 - scratch_366*scratch_444 - scratch_367*scratch_445 - scratch_368*scratch_446 - scratch_369*scratch_447 - scratch_370*scratch_448 - scratch_371*scratch_449 - scratch_372*scratch_450 - scratch_373*scratch_451 - scratch_374*scratch_452 - scratch_375*scratch_453 - scratch_376*scratch_454 - scratch_377*scratch_455;
    amrex::Real scratch_457 = scratch_378*scratch_456;
    amrex::Real scratch_458 = (-A[416]*scratch_1 + b[21] - scratch_381*(-A[416]*scratch_2 + A[435] - scratch_10*scratch_439 - scratch_109*scratch_447 - scratch_11*scratch_440 - scratch_131*scratch_448 - scratch_155*scratch_449 - scratch_181*scratch_450 - scratch_19*scratch_441 - scratch_209*scratch_451 - scratch_239*scratch_452 - scratch_271*scratch_453 - scratch_29*scratch_442 - scratch_305*scratch_454 - scratch_341*scratch_455 - scratch_359*scratch_457 - scratch_41*scratch_443 - scratch_444*scratch_55 - scratch_445*scratch_71 - scratch_446*scratch_89) - scratch_401*(-A[416]*scratch_382 + A[436] - scratch_383*scratch_439 - scratch_384*scratch_440 - scratch_385*scratch_441 - scratch_386*scratch_442 - scratch_387*scratch_443 - scratch_388*scratch_444 - scratch_389*scratch_445 - scratch_390*scratch_446 - scratch_391*scratch_447 - scratch_392*scratch_448 - scratch_393*scratch_449 - scratch_394*scratch_450 - scratch_395*scratch_451 - scratch_396*scratch_452 - scratch_397*scratch_453 - scratch_398*scratch_454 - scratch_399*scratch_455 - scratch_400*scratch_457) - scratch_402*scratch_439 - scratch_403*scratch_440 - scratch_404*scratch_441 - scratch_405*scratch_442 - scratch_406*scratch_443 - scratch_407*scratch_444 - scratch_408*scratch_445 - scratch_409*scratch_446 - scratch_410*scratch_447 - scratch_411*scratch_448 - scratch_412*scratch_449 - scratch_413*scratch_450 - scratch_414*scratch_451 - scratch_415*scratch_452 - scratch_416*scratch_453 - scratch_417*scratch_454 - scratch_418*scratch_455 - scratch_419*scratch_456)/(-A[416]*scratch_420 + A[437] - scratch_421*scratch_439 - scratch_422*scratch_440 - scratch_423*scratch_441 - scratch_424*scratch_442 - scratch_425*scratch_443 - scratch_426*scratch_444 - scratch_427*scratch_445 - scratch_428*scratch_446 - scratch_429*scratch_447 - scratch_430*scratch_448 - scratch_431*scratch_449 - scratch_432*scratch_450 - scratch_433*scratch_451 - scratch_434*scratch_452 - scratch_435*scratch_453 - scratch_436*scratch_454 - scratch_437*scratch_455 - scratch_438*scratch_457);
    amrex::Real scratch_459 = scratch_13 - scratch_14*scratch_6;
    amrex::Real scratch_460 = scratch_21 - scratch_22*scratch_6 - scratch_23*scratch_459;
    amrex::Real scratch_461 = scratch_31 - scratch_32*scratch_6 - scratch_33*scratch_459 - scratch_34*scratch_460;
    amrex::Real scratch_462 = scratch_43 - scratch_44*scratch_6 - scratch_45*scratch_459 - scratch_46*scratch_460 - scratch_461*scratch_47;
    amrex::Real scratch_463 = -scratch_459*scratch_59 - scratch_460*scratch_60 - scratch_461*scratch_61 - scratch_462*scratch_62 + scratch_57 - scratch_58*scratch_6;
    amrex::Real scratch_464 = -scratch_459*scratch_75 - scratch_460*scratch_76 - scratch_461*scratch_77 - scratch_462*scratch_78 - scratch_463*scratch_79 - scratch_6*scratch_74 + scratch_73;
    amrex::Real scratch_465 = -scratch_459*scratch_93 - scratch_460*scratch_94 - scratch_461*scratch_95 - scratch_462*scratch_96 - scratch_463*scratch_97 - scratch_464*scratch_98 - scratch_6*scratch_92 + scratch_91;
    amrex::Real scratch_466 = scratch_111 - scratch_112*scratch_6 - scratch_113*scratch_459 - scratch_114*scratch_460 - scratch_115*scratch_461 - scratch_116*scratch_462 - scratch_117*scratch_463 - scratch_118*scratch_464 - scratch_119*scratch_465;
    amrex::Real scratch_467 = scratch_133 - scratch_134*scratch_6 - scratch_135*scratch_459 - scratch_136*scratch_460 - scratch_137*scratch_461 - scratch_138*scratch_462 - scratch_139*scratch_463 - scratch_140*scratch_464 - scratch_141*scratch_465 - scratch_142*scratch_466;
    amrex::Real scratch_468 = scratch_157 - scratch_158*scratch_6 - scratch_159*scratch_459 - scratch_160*scratch_460 - scratch_161*scratch_461 - scratch_162*scratch_462 - scratch_163*scratch_463 - scratch_164*scratch_464 - scratch_165*scratch_465 - scratch_166*scratch_466 - scratch_167*scratch_467;
    amrex::Real scratch_469 = scratch_183 - scratch_184*scratch_6 - scratch_185*scratch_459 - scratch_186*scratch_460 - scratch_187*scratch_461 - scratch_188*scratch_462 - scratch_189*scratch_463 - scratch_190*scratch_464 - scratch_191*scratch_465 - scratch_192*scratch_466 - scratch_193*scratch_467 - scratch_194*scratch_468;
    amrex::Real scratch_470 = scratch_211 - scratch_212*scratch_6 - scratch_213*scratch_459 - scratch_214*scratch_460 - scratch_215*scratch_461 - scratch_216*scratch_462 - scratch_217*scratch_463 - scratch_218*scratch_464 - scratch_219*scratch_465 - scratch_220*scratch_466 - scratch_221*scratch_467 - scratch_222*scratch_468 - scratch_223*scratch_469;
    amrex::Real scratch_471 = scratch_241 - scratch_242*scratch_6 - scratch_243*scratch_459 - scratch_244*scratch_460 - scratch_245*scratch_461 - scratch_246*scratch_462 - scratch_247*scratch_463 - scratch_248*scratch_464 - scratch_249*scratch_465 - scratch_250*scratch_466 - scratch_251*scratch_467 - scratch_252*scratch_468 - scratch_253*scratch_469 - scratch_254*scratch_470;
    amrex::Real scratch_472 = scratch_273 - scratch_274*scratch_6 - scratch_275*scratch_459 - scratch_276*scratch_460 - scratch_277*scratch_461 - scratch_278*scratch_462 - scratch_279*scratch_463 - scratch_280*scratch_464 - scratch_281*scratch_465 - scratch_282*scratch_466 - scratch_283*scratch_467 - scratch_284*scratch_468 - scratch_285*scratch_469 - scratch_286*scratch_470 - scratch_287*scratch_471;
    amrex::Real scratch_473 = scratch_307 - scratch_308*scratch_6 - scratch_309*scratch_459 - scratch_310*scratch_460 - scratch_311*scratch_461 - scratch_312*scratch_462 - scratch_313*scratch_463 - scratch_314*scratch_464 - scratch_315*scratch_465 - scratch_316*scratch_466 - scratch_317*scratch_467 - scratch_318*scratch_468 - scratch_319*scratch_469 - scratch_320*scratch_470 - scratch_321*scratch_471 - scratch_322*scratch_472;
    amrex::Real scratch_474 = scratch_361 - scratch_362*scratch_6 - scratch_363*scratch_459 - scratch_364*scratch_460 - scratch_365*scratch_461 - scratch_366*scratch_462 - scratch_367*scratch_463 - scratch_368*scratch_464 - scratch_369*scratch_465 - scratch_370*scratch_466 - scratch_371*scratch_467 - scratch_372*scratch_468 - scratch_373*scratch_469 - scratch_374*scratch_470 - scratch_375*scratch_471 - scratch_376*scratch_472 - scratch_377*scratch_473;
    amrex::Real scratch_475 = scratch_378*scratch_474;
    amrex::Real scratch_476 = -scratch_14*scratch_23 + scratch_22;
    amrex::Real scratch_477 = -scratch_14*scratch_33 + scratch_32 - scratch_34*scratch_476;
    amrex::Real scratch_478 = -scratch_14*scratch_45 + scratch_44 - scratch_46*scratch_476 - scratch_47*scratch_477;
    amrex::Real scratch_479 = -scratch_14*scratch_59 - scratch_476*scratch_60 - scratch_477*scratch_61 - scratch_478*scratch_62 + scratch_58;
    amrex::Real scratch_480 = -scratch_14*scratch_75 - scratch_476*scratch_76 - scratch_477*scratch_77 - scratch_478*scratch_78 - scratch_479*scratch_79 + scratch_74;
    amrex::Real scratch_481 = -scratch_14*scratch_93 - scratch_476*scratch_94 - scratch_477*scratch_95 - scratch_478*scratch_96 - scratch_479*scratch_97 - scratch_480*scratch_98 + scratch_92;
    amrex::Real scratch_482 = scratch_112 - scratch_113*scratch_14 - scratch_114*scratch_476 - scratch_115*scratch_477 - scratch_116*scratch_478 - scratch_117*scratch_479 - scratch_118*scratch_480 - scratch_119*scratch_481;
    amrex::Real scratch_483 = scratch_134 - scratch_135*scratch_14 - scratch_136*scratch_476 - scratch_137*scratch_477 - scratch_138*scratch_478 - scratch_139*scratch_479 - scratch_140*scratch_480 - scratch_141*scratch_481 - scratch_142*scratch_482;
    amrex::Real scratch_484 = -scratch_14*scratch_159 + scratch_158 - scratch_160*scratch_476 - scratch_161*scratch_477 - scratch_162*scratch_478 - scratch_163*scratch_479 - scratch_164*scratch_480 - scratch_165*scratch_481 - scratch_166*scratch_482 - scratch_167*scratch_483;
    amrex::Real scratch_485 = -scratch_14*scratch_185 + scratch_184 - scratch_186*scratch_476 - scratch_187*scratch_477 - scratch_188*scratch_478 - scratch_189*scratch_479 - scratch_190*scratch_480 - scratch_191*scratch_481 - scratch_192*scratch_482 - scratch_193*scratch_483 - scratch_194*scratch_484;
    amrex::Real scratch_486 = -scratch_14*scratch_213 + scratch_212 - scratch_214*scratch_476 - scratch_215*scratch_477 - scratch_216*scratch_478 - scratch_217*scratch_479 - scratch_218*scratch_480 - scratch_219*scratch_481 - scratch_220*scratch_482 - scratch_221*scratch_483 - scratch_222*scratch_484 - scratch_223*scratch_485;
    amrex::Real scratch_487 = -scratch_14*scratch_243 + scratch_242 - scratch_244*scratch_476 - scratch_245*scratch_477 - scratch_246*scratch_478 - scratch_247*scratch_479 - scratch_248*scratch_480 - scratch_249*scratch_481 - scratch_250*scratch_482 - scratch_251*scratch_483 - scratch_252*scratch_484 - scratch_253*scratch_485 - scratch_254*scratch_486;
    amrex::Real scratch_488 = -scratch_14*scratch_275 + scratch_274 - scratch_276*scratch_476 - scratch_277*scratch_477 - scratch_278*scratch_478 - scratch_279*scratch_479 - scratch_280*scratch_480 - scratch_281*scratch_481 - scratch_282*scratch_482 - scratch_283*scratch_483 - scratch_284*scratch_484 - scratch_285*scratch_485 - scratch_286*scratch_486 - scratch_287*scratch_487;
    amrex::Real scratch_489 = -scratch_14*scratch_309 + scratch_308 - scratch_310*scratch_476 - scratch_311*scratch_477 - scratch_312*scratch_478 - scratch_313*scratch_479 - scratch_314*scratch_480 - scratch_315*scratch_481 - scratch_316*scratch_482 - scratch_317*scratch_483 - scratch_318*scratch_484 - scratch_319*scratch_485 - scratch_320*scratch_486 - scratch_321*scratch_487 - scratch_322*scratch_488;
    amrex::Real scratch_490 = -scratch_14*scratch_363 + scratch_362 - scratch_364*scratch_476 - scratch_365*scratch_477 - scratch_366*scratch_478 - scratch_367*scratch_479 - scratch_368*scratch_480 - scratch_369*scratch_481 - scratch_370*scratch_482 - scratch_371*scratch_483 - scratch_372*scratch_484 - scratch_373*scratch_485 - scratch_374*scratch_486 - scratch_375*scratch_487 - scratch_376*scratch_488 - scratch_377*scratch_489;
    amrex::Real scratch_491 = scratch_378*scratch_490;
    amrex::Real scratch_492 = -scratch_23*scratch_34 + scratch_33;
    amrex::Real scratch_493 = -scratch_23*scratch_46 + scratch_45 - scratch_47*scratch_492;
    amrex::Real scratch_494 = -scratch_23*scratch_60 - scratch_492*scratch_61 - scratch_493*scratch_62 + scratch_59;
    amrex::Real scratch_495 = -scratch_23*scratch_76 - scratch_492*scratch_77 - scratch_493*scratch_78 - scratch_494*scratch_79 + scratch_75;
    amrex::Real scratch_496 = -scratch_23*scratch_94 - scratch_492*scratch_95 - scratch_493*scratch_96 - scratch_494*scratch_97 - scratch_495*scratch_98 + scratch_93;
    amrex::Real scratch_497 = scratch_113 - scratch_114*scratch_23 - scratch_115*scratch_492 - scratch_116*scratch_493 - scratch_117*scratch_494 - scratch_118*scratch_495 - scratch_119*scratch_496;
    amrex::Real scratch_498 = scratch_135 - scratch_136*scratch_23 - scratch_137*scratch_492 - scratch_138*scratch_493 - scratch_139*scratch_494 - scratch_140*scratch_495 - scratch_141*scratch_496 - scratch_142*scratch_497;
    amrex::Real scratch_499 = scratch_159 - scratch_160*scratch_23 - scratch_161*scratch_492 - scratch_162*scratch_493 - scratch_163*scratch_494 - scratch_164*scratch_495 - scratch_165*scratch_496 - scratch_166*scratch_497 - scratch_167*scratch_498;
    amrex::Real scratch_500 = scratch_185 - scratch_186*scratch_23 - scratch_187*scratch_492 - scratch_188*scratch_493 - scratch_189*scratch_494 - scratch_190*scratch_495 - scratch_191*scratch_496 - scratch_192*scratch_497 - scratch_193*scratch_498 - scratch_194*scratch_499;
    amrex::Real scratch_501 = scratch_213 - scratch_214*scratch_23 - scratch_215*scratch_492 - scratch_216*scratch_493 - scratch_217*scratch_494 - scratch_218*scratch_495 - scratch_219*scratch_496 - scratch_220*scratch_497 - scratch_221*scratch_498 - scratch_222*scratch_499 - scratch_223*scratch_500;
    amrex::Real scratch_502 = -scratch_23*scratch_244 + scratch_243 - scratch_245*scratch_492 - scratch_246*scratch_493 - scratch_247*scratch_494 - scratch_248*scratch_495 - scratch_249*scratch_496 - scratch_250*scratch_497 - scratch_251*scratch_498 - scratch_252*scratch_499 - scratch_253*scratch_500 - scratch_254*scratch_501;
    amrex::Real scratch_503 = -scratch_23*scratch_276 + scratch_275 - scratch_277*scratch_492 - scratch_278*scratch_493 - scratch_279*scratch_494 - scratch_280*scratch_495 - scratch_281*scratch_496 - scratch_282*scratch_497 - scratch_283*scratch_498 - scratch_284*scratch_499 - scratch_285*scratch_500 - scratch_286*scratch_501 - scratch_287*scratch_502;
    amrex::Real scratch_504 = -scratch_23*scratch_310 + scratch_309 - scratch_311*scratch_492 - scratch_312*scratch_493 - scratch_313*scratch_494 - scratch_314*scratch_495 - scratch_315*scratch_496 - scratch_316*scratch_497 - scratch_317*scratch_498 - scratch_318*scratch_499 - scratch_319*scratch_500 - scratch_320*scratch_501 - scratch_321*scratch_502 - scratch_322*scratch_503;
    amrex::Real scratch_505 = -scratch_23*scratch_364 + scratch_363 - scratch_365*scratch_492 - scratch_366*scratch_493 - scratch_367*scratch_494 - scratch_368*scratch_495 - scratch_369*scratch_496 - scratch_370*scratch_497 - scratch_371*scratch_498 - scratch_372*scratch_499 - scratch_373*scratch_500 - scratch_374*scratch_501 - scratch_375*scratch_502 - scratch_376*scratch_503 - scratch_377*scratch_504;
    amrex::Real scratch_506 = scratch_378*scratch_505;
    amrex::Real scratch_507 = -scratch_34*scratch_47 + scratch_46;
    amrex::Real scratch_508 = -scratch_34*scratch_61 - scratch_507*scratch_62 + scratch_60;
    amrex::Real scratch_509 = -scratch_34*scratch_77 - scratch_507*scratch_78 - scratch_508*scratch_79 + scratch_76;
    amrex::Real scratch_510 = -scratch_34*scratch_95 - scratch_507*scratch_96 - scratch_508*scratch_97 - scratch_509*scratch_98 + scratch_94;
    amrex::Real scratch_511 = scratch_114 - scratch_115*scratch_34 - scratch_116*scratch_507 - scratch_117*scratch_508 - scratch_118*scratch_509 - scratch_119*scratch_510;
    amrex::Real scratch_512 = scratch_136 - scratch_137*scratch_34 - scratch_138*scratch_507 - scratch_139*scratch_508 - scratch_140*scratch_509 - scratch_141*scratch_510 - scratch_142*scratch_511;
    amrex::Real scratch_513 = scratch_160 - scratch_161*scratch_34 - scratch_162*scratch_507 - scratch_163*scratch_508 - scratch_164*scratch_509 - scratch_165*scratch_510 - scratch_166*scratch_511 - scratch_167*scratch_512;
    amrex::Real scratch_514 = scratch_186 - scratch_187*scratch_34 - scratch_188*scratch_507 - scratch_189*scratch_508 - scratch_190*scratch_509 - scratch_191*scratch_510 - scratch_192*scratch_511 - scratch_193*scratch_512 - scratch_194*scratch_513;
    amrex::Real scratch_515 = scratch_214 - scratch_215*scratch_34 - scratch_216*scratch_507 - scratch_217*scratch_508 - scratch_218*scratch_509 - scratch_219*scratch_510 - scratch_220*scratch_511 - scratch_221*scratch_512 - scratch_222*scratch_513 - scratch_223*scratch_514;
    amrex::Real scratch_516 = scratch_244 - scratch_245*scratch_34 - scratch_246*scratch_507 - scratch_247*scratch_508 - scratch_248*scratch_509 - scratch_249*scratch_510 - scratch_250*scratch_511 - scratch_251*scratch_512 - scratch_252*scratch_513 - scratch_253*scratch_514 - scratch_254*scratch_515;
    amrex::Real scratch_517 = scratch_276 - scratch_277*scratch_34 - scratch_278*scratch_507 - scratch_279*scratch_508 - scratch_280*scratch_509 - scratch_281*scratch_510 - scratch_282*scratch_511 - scratch_283*scratch_512 - scratch_284*scratch_513 - scratch_285*scratch_514 - scratch_286*scratch_515 - scratch_287*scratch_516;
    amrex::Real scratch_518 = scratch_310 - scratch_311*scratch_34 - scratch_312*scratch_507 - scratch_313*scratch_508 - scratch_314*scratch_509 - scratch_315*scratch_510 - scratch_316*scratch_511 - scratch_317*scratch_512 - scratch_318*scratch_513 - scratch_319*scratch_514 - scratch_320*scratch_515 - scratch_321*scratch_516 - scratch_322*scratch_517;
    amrex::Real scratch_519 = -scratch_34*scratch_365 + scratch_364 - scratch_366*scratch_507 - scratch_367*scratch_508 - scratch_368*scratch_509 - scratch_369*scratch_510 - scratch_370*scratch_511 - scratch_371*scratch_512 - scratch_372*scratch_513 - scratch_373*scratch_514 - scratch_374*scratch_515 - scratch_375*scratch_516 - scratch_376*scratch_517 - scratch_377*scratch_518;
    amrex::Real scratch_520 = scratch_378*scratch_519;
    amrex::Real scratch_521 = -scratch_47*scratch_62 + scratch_61;
    amrex::Real scratch_522 = -scratch_47*scratch_78 - scratch_521*scratch_79 + scratch_77;
    amrex::Real scratch_523 = -scratch_47*scratch_96 - scratch_521*scratch_97 - scratch_522*scratch_98 + scratch_95;
    amrex::Real scratch_524 = scratch_115 - scratch_116*scratch_47 - scratch_117*scratch_521 - scratch_118*scratch_522 - scratch_119*scratch_523;
    amrex::Real scratch_525 = scratch_137 - scratch_138*scratch_47 - scratch_139*scratch_521 - scratch_140*scratch_522 - scratch_141*scratch_523 - scratch_142*scratch_524;
    amrex::Real scratch_526 = scratch_161 - scratch_162*scratch_47 - scratch_163*scratch_521 - scratch_164*scratch_522 - scratch_165*scratch_523 - scratch_166*scratch_524 - scratch_167*scratch_525;
    amrex::Real scratch_527 = scratch_187 - scratch_188*scratch_47 - scratch_189*scratch_521 - scratch_190*scratch_522 - scratch_191*scratch_523 - scratch_192*scratch_524 - scratch_193*scratch_525 - scratch_194*scratch_526;
    amrex::Real scratch_528 = scratch_215 - scratch_216*scratch_47 - scratch_217*scratch_521 - scratch_218*scratch_522 - scratch_219*scratch_523 - scratch_220*scratch_524 - scratch_221*scratch_525 - scratch_222*scratch_526 - scratch_223*scratch_527;
    amrex::Real scratch_529 = scratch_245 - scratch_246*scratch_47 - scratch_247*scratch_521 - scratch_248*scratch_522 - scratch_249*scratch_523 - scratch_250*scratch_524 - scratch_251*scratch_525 - scratch_252*scratch_526 - scratch_253*scratch_527 - scratch_254*scratch_528;
    amrex::Real scratch_530 = scratch_277 - scratch_278*scratch_47 - scratch_279*scratch_521 - scratch_280*scratch_522 - scratch_281*scratch_523 - scratch_282*scratch_524 - scratch_283*scratch_525 - scratch_284*scratch_526 - scratch_285*scratch_527 - scratch_286*scratch_528 - scratch_287*scratch_529;
    amrex::Real scratch_531 = scratch_311 - scratch_312*scratch_47 - scratch_313*scratch_521 - scratch_314*scratch_522 - scratch_315*scratch_523 - scratch_316*scratch_524 - scratch_317*scratch_525 - scratch_318*scratch_526 - scratch_319*scratch_527 - scratch_320*scratch_528 - scratch_321*scratch_529 - scratch_322*scratch_530;
    amrex::Real scratch_532 = scratch_365 - scratch_366*scratch_47 - scratch_367*scratch_521 - scratch_368*scratch_522 - scratch_369*scratch_523 - scratch_370*scratch_524 - scratch_371*scratch_525 - scratch_372*scratch_526 - scratch_373*scratch_527 - scratch_374*scratch_528 - scratch_375*scratch_529 - scratch_376*scratch_530 - scratch_377*scratch_531;
    amrex::Real scratch_533 = scratch_378*scratch_532;
    amrex::Real scratch_534 = -scratch_62*scratch_79 + scratch_78;
    amrex::Real scratch_535 = -scratch_534*scratch_98 - scratch_62*scratch_97 + scratch_96;
    amrex::Real scratch_536 = scratch_116 - scratch_117*scratch_62 - scratch_118*scratch_534 - scratch_119*scratch_535;
    amrex::Real scratch_537 = scratch_138 - scratch_139*scratch_62 - scratch_140*scratch_534 - scratch_141*scratch_535 - scratch_142*scratch_536;
    amrex::Real scratch_538 = scratch_162 - scratch_163*scratch_62 - scratch_164*scratch_534 - scratch_165*scratch_535 - scratch_166*scratch_536 - scratch_167*scratch_537;
    amrex::Real scratch_539 = scratch_188 - scratch_189*scratch_62 - scratch_190*scratch_534 - scratch_191*scratch_535 - scratch_192*scratch_536 - scratch_193*scratch_537 - scratch_194*scratch_538;
    amrex::Real scratch_540 = scratch_216 - scratch_217*scratch_62 - scratch_218*scratch_534 - scratch_219*scratch_535 - scratch_220*scratch_536 - scratch_221*scratch_537 - scratch_222*scratch_538 - scratch_223*scratch_539;
    amrex::Real scratch_541 = scratch_246 - scratch_247*scratch_62 - scratch_248*scratch_534 - scratch_249*scratch_535 - scratch_250*scratch_536 - scratch_251*scratch_537 - scratch_252*scratch_538 - scratch_253*scratch_539 - scratch_254*scratch_540;
    amrex::Real scratch_542 = scratch_278 - scratch_279*scratch_62 - scratch_280*scratch_534 - scratch_281*scratch_535 - scratch_282*scratch_536 - scratch_283*scratch_537 - scratch_284*scratch_538 - scratch_285*scratch_539 - scratch_286*scratch_540 - scratch_287*scratch_541;
    amrex::Real scratch_543 = scratch_312 - scratch_313*scratch_62 - scratch_314*scratch_534 - scratch_315*scratch_535 - scratch_316*scratch_536 - scratch_317*scratch_537 - scratch_318*scratch_538 - scratch_319*scratch_539 - scratch_320*scratch_540 - scratch_321*scratch_541 - scratch_322*scratch_542;
    amrex::Real scratch_544 = scratch_366 - scratch_367*scratch_62 - scratch_368*scratch_534 - scratch_369*scratch_535 - scratch_370*scratch_536 - scratch_371*scratch_537 - scratch_372*scratch_538 - scratch_373*scratch_539 - scratch_374*scratch_540 - scratch_375*scratch_541 - scratch_376*scratch_542 - scratch_377*scratch_543;
    amrex::Real scratch_545 = scratch_378*scratch_544;
    amrex::Real scratch_546 = -scratch_79*scratch_98 + scratch_97;
    amrex::Real scratch_547 = scratch_117 - scratch_118*scratch_79 - scratch_119*scratch_546;
    amrex::Real scratch_548 = scratch_139 - scratch_140*scratch_79 - scratch_141*scratch_546 - scratch_142*scratch_547;
    amrex::Real scratch_549 = scratch_163 - scratch_164*scratch_79 - scratch_165*scratch_546 - scratch_166*scratch_547 - scratch_167*scratch_548;
    amrex::Real scratch_550 = scratch_189 - scratch_190*scratch_79 - scratch_191*scratch_546 - scratch_192*scratch_547 - scratch_193*scratch_548 - scratch_194*scratch_549;
    amrex::Real scratch_551 = scratch_217 - scratch_218*scratch_79 - scratch_219*scratch_546 - scratch_220*scratch_547 - scratch_221*scratch_548 - scratch_222*scratch_549 - scratch_223*scratch_550;
    amrex::Real scratch_552 = scratch_247 - scratch_248*scratch_79 - scratch_249*scratch_546 - scratch_250*scratch_547 - scratch_251*scratch_548 - scratch_252*scratch_549 - scratch_253*scratch_550 - scratch_254*scratch_551;
    amrex::Real scratch_553 = scratch_279 - scratch_280*scratch_79 - scratch_281*scratch_546 - scratch_282*scratch_547 - scratch_283*scratch_548 - scratch_284*scratch_549 - scratch_285*scratch_550 - scratch_286*scratch_551 - scratch_287*scratch_552;
    amrex::Real scratch_554 = scratch_313 - scratch_314*scratch_79 - scratch_315*scratch_546 - scratch_316*scratch_547 - scratch_317*scratch_548 - scratch_318*scratch_549 - scratch_319*scratch_550 - scratch_320*scratch_551 - scratch_321*scratch_552 - scratch_322*scratch_553;
    amrex::Real scratch_555 = scratch_367 - scratch_368*scratch_79 - scratch_369*scratch_546 - scratch_370*scratch_547 - scratch_371*scratch_548 - scratch_372*scratch_549 - scratch_373*scratch_550 - scratch_374*scratch_551 - scratch_375*scratch_552 - scratch_376*scratch_553 - scratch_377*scratch_554;
    amrex::Real scratch_556 = scratch_378*scratch_555;
    amrex::Real scratch_557 = scratch_118 - scratch_119*scratch_98;
    amrex::Real scratch_558 = scratch_140 - scratch_141*scratch_98 - scratch_142*scratch_557;
    amrex::Real scratch_559 = scratch_164 - scratch_165*scratch_98 - scratch_166*scratch_557 - scratch_167*scratch_558;
    amrex::Real scratch_560 = scratch_190 - scratch_191*scratch_98 - scratch_192*scratch_557 - scratch_193*scratch_558 - scratch_194*scratch_559;
    amrex::Real scratch_561 = scratch_218 - scratch_219*scratch_98 - scratch_220*scratch_557 - scratch_221*scratch_558 - scratch_222*scratch_559 - scratch_223*scratch_560;
    amrex::Real scratch_562 = scratch_248 - scratch_249*scratch_98 - scratch_250*scratch_557 - scratch_251*scratch_558 - scratch_252*scratch_559 - scratch_253*scratch_560 - scratch_254*scratch_561;
    amrex::Real scratch_563 = scratch_280 - scratch_281*scratch_98 - scratch_282*scratch_557 - scratch_283*scratch_558 - scratch_284*scratch_559 - scratch_285*scratch_560 - scratch_286*scratch_561 - scratch_287*scratch_562;
    amrex::Real scratch_564 = scratch_314 - scratch_315*scratch_98 - scratch_316*scratch_557 - scratch_317*scratch_558 - scratch_318*scratch_559 - scratch_319*scratch_560 - scratch_320*scratch_561 - scratch_321*scratch_562 - scratch_322*scratch_563;
    amrex::Real scratch_565 = scratch_368 - scratch_369*scratch_98 - scratch_370*scratch_557 - scratch_371*scratch_558 - scratch_372*scratch_559 - scratch_373*scratch_560 - scratch_374*scratch_561 - scratch_375*scratch_562 - scratch_376*scratch_563 - scratch_377*scratch_564;
    amrex::Real scratch_566 = scratch_378*scratch_565;
    amrex::Real scratch_567 = -scratch_119*scratch_142 + scratch_141;
    amrex::Real scratch_568 = -scratch_119*scratch_166 + scratch_165 - scratch_167*scratch_567;
    amrex::Real scratch_569 = -scratch_119*scratch_192 + scratch_191 - scratch_193*scratch_567 - scratch_194*scratch_568;
    amrex::Real scratch_570 = -scratch_119*scratch_220 + scratch_219 - scratch_221*scratch_567 - scratch_222*scratch_568 - scratch_223*scratch_569;
    amrex::Real scratch_571 = -scratch_119*scratch_250 + scratch_249 - scratch_251*scratch_567 - scratch_252*scratch_568 - scratch_253*scratch_569 - scratch_254*scratch_570;
    amrex::Real scratch_572 = -scratch_119*scratch_282 + scratch_281 - scratch_283*scratch_567 - scratch_284*scratch_568 - scratch_285*scratch_569 - scratch_286*scratch_570 - scratch_287*scratch_571;
    amrex::Real scratch_573 = -scratch_119*scratch_316 + scratch_315 - scratch_317*scratch_567 - scratch_318*scratch_568 - scratch_319*scratch_569 - scratch_320*scratch_570 - scratch_321*scratch_571 - scratch_322*scratch_572;
    amrex::Real scratch_574 = -scratch_119*scratch_370 + scratch_369 - scratch_371*scratch_567 - scratch_372*scratch_568 - scratch_373*scratch_569 - scratch_374*scratch_570 - scratch_375*scratch_571 - scratch_376*scratch_572 - scratch_377*scratch_573;
    amrex::Real scratch_575 = scratch_378*scratch_574;
    amrex::Real scratch_576 = -scratch_142*scratch_167 + scratch_166;
    amrex::Real scratch_577 = -scratch_142*scratch_193 + scratch_192 - scratch_194*scratch_576;
    amrex::Real scratch_578 = -scratch_142*scratch_221 + scratch_220 - scratch_222*scratch_576 - scratch_223*scratch_577;
    amrex::Real scratch_579 = -scratch_142*scratch_251 + scratch_250 - scratch_252*scratch_576 - scratch_253*scratch_577 - scratch_254*scratch_578;
    amrex::Real scratch_580 = -scratch_142*scratch_283 + scratch_282 - scratch_284*scratch_576 - scratch_285*scratch_577 - scratch_286*scratch_578 - scratch_287*scratch_579;
    amrex::Real scratch_581 = -scratch_142*scratch_317 + scratch_316 - scratch_318*scratch_576 - scratch_319*scratch_577 - scratch_320*scratch_578 - scratch_321*scratch_579 - scratch_322*scratch_580;
    amrex::Real scratch_582 = -scratch_142*scratch_371 + scratch_370 - scratch_372*scratch_576 - scratch_373*scratch_577 - scratch_374*scratch_578 - scratch_375*scratch_579 - scratch_376*scratch_580 - scratch_377*scratch_581;
    amrex::Real scratch_583 = scratch_378*scratch_582;
    amrex::Real scratch_584 = -scratch_167*scratch_194 + scratch_193;
    amrex::Real scratch_585 = -scratch_167*scratch_222 + scratch_221 - scratch_223*scratch_584;
    amrex::Real scratch_586 = -scratch_167*scratch_252 + scratch_251 - scratch_253*scratch_584 - scratch_254*scratch_585;
    amrex::Real scratch_587 = -scratch_167*scratch_284 + scratch_283 - scratch_285*scratch_584 - scratch_286*scratch_585 - scratch_287*scratch_586;
    amrex::Real scratch_588 = -scratch_167*scratch_318 + scratch_317 - scratch_319*scratch_584 - scratch_320*scratch_585 - scratch_321*scratch_586 - scratch_322*scratch_587;
    amrex::Real scratch_589 = -scratch_167*scratch_372 + scratch_371 - scratch_373*scratch_584 - scratch_374*scratch_585 - scratch_375*scratch_586 - scratch_376*scratch_587 - scratch_377*scratch_588;
    amrex::Real scratch_590 = scratch_378*scratch_589;
    amrex::Real scratch_591 = -scratch_194*scratch_223 + scratch_222;
    amrex::Real scratch_592 = -scratch_194*scratch_253 + scratch_252 - scratch_254*scratch_591;
    amrex::Real scratch_593 = -scratch_194*scratch_285 + scratch_284 - scratch_286*scratch_591 - scratch_287*scratch_592;
    amrex::Real scratch_594 = -scratch_194*scratch_319 + scratch_318 - scratch_320*scratch_591 - scratch_321*scratch_592 - scratch_322*scratch_593;
    amrex::Real scratch_595 = -scratch_194*scratch_373 + scratch_372 - scratch_374*scratch_591 - scratch_375*scratch_592 - scratch_376*scratch_593 - scratch_377*scratch_594;
    amrex::Real scratch_596 = scratch_378*scratch_595;
    amrex::Real scratch_597 = -scratch_223*scratch_254 + scratch_253;
    amrex::Real scratch_598 = -scratch_223*scratch_286 + scratch_285 - scratch_287*scratch_597;
    amrex::Real scratch_599 = -scratch_223*scratch_320 + scratch_319 - scratch_321*scratch_597 - scratch_322*scratch_598;
    amrex::Real scratch_600 = -scratch_223*scratch_374 + scratch_373 - scratch_375*scratch_597 - scratch_376*scratch_598 - scratch_377*scratch_599;
    amrex::Real scratch_601 = scratch_378*scratch_600;
    amrex::Real scratch_602 = -scratch_254*scratch_287 + scratch_286;
    amrex::Real scratch_603 = -scratch_254*scratch_321 + scratch_320 - scratch_322*scratch_602;
    amrex::Real scratch_604 = -scratch_254*scratch_375 + scratch_374 - scratch_376*scratch_602 - scratch_377*scratch_603;
    amrex::Real scratch_605 = scratch_378*scratch_604;
    amrex::Real scratch_606 = -scratch_287*scratch_322 + scratch_321;
    amrex::Real scratch_607 = -scratch_287*scratch_376 + scratch_375 - scratch_377*scratch_606;
    amrex::Real scratch_608 = scratch_378*scratch_607;
    amrex::Real scratch_609 = -scratch_322*scratch_377 + scratch_376;
    amrex::Real scratch_610 = scratch_378*scratch_609;
    amrex::Real scratch_611 = scratch_377*scratch_378;

    x[0] = scratch_1 - scratch_120*scratch_411 - scratch_143*scratch_412 - scratch_15*scratch_404 - scratch_168*scratch_413 - scratch_195*scratch_414 - scratch_224*scratch_415 - scratch_24*scratch_405 - scratch_255*scratch_416 - scratch_288*scratch_417 - scratch_323*scratch_418 - scratch_35*scratch_406 - scratch_379*scratch_419 - scratch_381*(-scratch_10*scratch_4 - scratch_109*scratch_99 - scratch_11*scratch_7 - scratch_120*scratch_131 - scratch_143*scratch_155 - scratch_15*scratch_19 - scratch_168*scratch_181 - scratch_195*scratch_209 + scratch_2 - scratch_224*scratch_239 - scratch_24*scratch_29 - scratch_255*scratch_271 - scratch_288*scratch_305 - scratch_323*scratch_341 - scratch_35*scratch_41 - scratch_359*scratch_380 - scratch_48*scratch_55 - scratch_63*scratch_71 - scratch_80*scratch_89) - scratch_4*scratch_402 - scratch_401*(-scratch_120*scratch_392 - scratch_143*scratch_393 - scratch_15*scratch_385 - scratch_168*scratch_394 - scratch_195*scratch_395 - scratch_224*scratch_396 - scratch_24*scratch_386 - scratch_255*scratch_397 - scratch_288*scratch_398 - scratch_323*scratch_399 - scratch_35*scratch_387 - scratch_380*scratch_400 + scratch_382 - scratch_383*scratch_4 - scratch_384*scratch_7 - scratch_388*scratch_48 - scratch_389*scratch_63 - scratch_390*scratch_80 - scratch_391*scratch_99) - scratch_403*scratch_7 - scratch_407*scratch_48 - scratch_408*scratch_63 - scratch_409*scratch_80 - scratch_410*scratch_99 - scratch_458*(-scratch_120*scratch_430 - scratch_143*scratch_431 - scratch_15*scratch_423 - scratch_168*scratch_432 - scratch_195*scratch_433 - scratch_224*scratch_434 - scratch_24*scratch_424 - scratch_255*scratch_435 - scratch_288*scratch_436 - scratch_323*scratch_437 - scratch_35*scratch_425 - scratch_380*scratch_438 - scratch_4*scratch_421 + scratch_420 - scratch_422*scratch_7 - scratch_426*scratch_48 - scratch_427*scratch_63 - scratch_428*scratch_80 - scratch_429*scratch_99);
    x[1] = -scratch_381*(scratch_10 - scratch_109*scratch_465 - scratch_11*scratch_6 - scratch_131*scratch_466 - scratch_155*scratch_467 - scratch_181*scratch_468 - scratch_19*scratch_459 - scratch_209*scratch_469 - scratch_239*scratch_470 - scratch_271*scratch_471 - scratch_29*scratch_460 - scratch_305*scratch_472 - scratch_341*scratch_473 - scratch_359*scratch_475 - scratch_41*scratch_461 - scratch_462*scratch_55 - scratch_463*scratch_71 - scratch_464*scratch_89) - scratch_401*(scratch_383 - scratch_384*scratch_6 - scratch_385*scratch_459 - scratch_386*scratch_460 - scratch_387*scratch_461 - scratch_388*scratch_462 - scratch_389*scratch_463 - scratch_390*scratch_464 - scratch_391*scratch_465 - scratch_392*scratch_466 - scratch_393*scratch_467 - scratch_394*scratch_468 - scratch_395*scratch_469 - scratch_396*scratch_470 - scratch_397*scratch_471 - scratch_398*scratch_472 - scratch_399*scratch_473 - scratch_400*scratch_475) + scratch_402 - scratch_403*scratch_6 - scratch_404*scratch_459 - scratch_405*scratch_460 - scratch_406*scratch_461 - scratch_407*scratch_462 - scratch_408*scratch_463 - scratch_409*scratch_464 - scratch_410*scratch_465 - scratch_411*scratch_466 - scratch_412*scratch_467 - scratch_413*scratch_468 - scratch_414*scratch_469 - scratch_415*scratch_470 - scratch_416*scratch_471 - scratch_417*scratch_472 - scratch_418*scratch_473 - scratch_419*scratch_474 - scratch_458*(scratch_421 - scratch_422*scratch_6 - scratch_423*scratch_459 - scratch_424*scratch_460 - scratch_425*scratch_461 - scratch_426*scratch_462 - scratch_427*scratch_463 - scratch_428*scratch_464 - scratch_429*scratch_465 - scratch_430*scratch_466 - scratch_431*scratch_467 - scratch_432*scratch_468 - scratch_433*scratch_469 - scratch_434*scratch_470 - scratch_435*scratch_471 - scratch_436*scratch_472 - scratch_437*scratch_473 - scratch_438*scratch_475);
    x[2] = -scratch_14*scratch_404 - scratch_381*(-scratch_109*scratch_481 + scratch_11 - scratch_131*scratch_482 - scratch_14*scratch_19 - scratch_155*scratch_483 - scratch_181*scratch_484 - scratch_209*scratch_485 - scratch_239*scratch_486 - scratch_271*scratch_487 - scratch_29*scratch_476 - scratch_305*scratch_488 - scratch_341*scratch_489 - scratch_359*scratch_491 - scratch_41*scratch_477 - scratch_478*scratch_55 - scratch_479*scratch_71 - scratch_480*scratch_89) - scratch_401*(-scratch_14*scratch_385 + scratch_384 - scratch_386*scratch_476 - scratch_387*scratch_477 - scratch_388*scratch_478 - scratch_389*scratch_479 - scratch_390*scratch_480 - scratch_391*scratch_481 - scratch_392*scratch_482 - scratch_393*scratch_483 - scratch_394*scratch_484 - scratch_395*scratch_485 - scratch_396*scratch_486 - scratch_397*scratch_487 - scratch_398*scratch_488 - scratch_399*scratch_489 - scratch_400*scratch_491) + scratch_403 - scratch_405*scratch_476 - scratch_406*scratch_477 - scratch_407*scratch_478 - scratch_408*scratch_479 - scratch_409*scratch_480 - scratch_410*scratch_481 - scratch_411*scratch_482 - scratch_412*scratch_483 - scratch_413*scratch_484 - scratch_414*scratch_485 - scratch_415*scratch_486 - scratch_416*scratch_487 - scratch_417*scratch_488 - scratch_418*scratch_489 - scratch_419*scratch_490 - scratch_458*(-scratch_14*scratch_423 + scratch_422 - scratch_424*scratch_476 - scratch_425*scratch_477 - scratch_426*scratch_478 - scratch_427*scratch_479 - scratch_428*scratch_480 - scratch_429*scratch_481 - scratch_430*scratch_482 - scratch_431*scratch_483 - scratch_432*scratch_484 - scratch_433*scratch_485 - scratch_434*scratch_486 - scratch_435*scratch_487 - scratch_436*scratch_488 - scratch_437*scratch_489 - scratch_438*scratch_491);
    x[3] = -scratch_23*scratch_405 - scratch_381*(-scratch_109*scratch_496 - scratch_131*scratch_497 - scratch_155*scratch_498 - scratch_181*scratch_499 + scratch_19 - scratch_209*scratch_500 - scratch_23*scratch_29 - scratch_239*scratch_501 - scratch_271*scratch_502 - scratch_305*scratch_503 - scratch_341*scratch_504 - scratch_359*scratch_506 - scratch_41*scratch_492 - scratch_493*scratch_55 - scratch_494*scratch_71 - scratch_495*scratch_89) - scratch_401*(-scratch_23*scratch_386 + scratch_385 - scratch_387*scratch_492 - scratch_388*scratch_493 - scratch_389*scratch_494 - scratch_390*scratch_495 - scratch_391*scratch_496 - scratch_392*scratch_497 - scratch_393*scratch_498 - scratch_394*scratch_499 - scratch_395*scratch_500 - scratch_396*scratch_501 - scratch_397*scratch_502 - scratch_398*scratch_503 - scratch_399*scratch_504 - scratch_400*scratch_506) + scratch_404 - scratch_406*scratch_492 - scratch_407*scratch_493 - scratch_408*scratch_494 - scratch_409*scratch_495 - scratch_410*scratch_496 - scratch_411*scratch_497 - scratch_412*scratch_498 - scratch_413*scratch_499 - scratch_414*scratch_500 - scratch_415*scratch_501 - scratch_416*scratch_502 - scratch_417*scratch_503 - scratch_418*scratch_504 - scratch_419*scratch_505 - scratch_458*(-scratch_23*scratch_424 + scratch_423 - scratch_425*scratch_492 - scratch_426*scratch_493 - scratch_427*scratch_494 - scratch_428*scratch_495 - scratch_429*scratch_496 - scratch_430*scratch_497 - scratch_431*scratch_498 - scratch_432*scratch_499 - scratch_433*scratch_500 - scratch_434*scratch_501 - scratch_435*scratch_502 - scratch_436*scratch_503 - scratch_437*scratch_504 - scratch_438*scratch_506);
    x[4] = -scratch_34*scratch_406 - scratch_381*(-scratch_109*scratch_510 - scratch_131*scratch_511 - scratch_155*scratch_512 - scratch_181*scratch_513 - scratch_209*scratch_514 - scratch_239*scratch_515 - scratch_271*scratch_516 + scratch_29 - scratch_305*scratch_517 - scratch_34*scratch_41 - scratch_341*scratch_518 - scratch_359*scratch_520 - scratch_507*scratch_55 - scratch_508*scratch_71 - scratch_509*scratch_89) - scratch_401*(-scratch_34*scratch_387 + scratch_386 - scratch_388*scratch_507 - scratch_389*scratch_508 - scratch_390*scratch_509 - scratch_391*scratch_510 - scratch_392*scratch_511 - scratch_393*scratch_512 - scratch_394*scratch_513 - scratch_395*scratch_514 - scratch_396*scratch_515 - scratch_397*scratch_516 - scratch_398*scratch_517 - scratch_399*scratch_518 - scratch_400*scratch_520) + scratch_405 - scratch_407*scratch_507 - scratch_408*scratch_508 - scratch_409*scratch_509 - scratch_410*scratch_510 - scratch_411*scratch_511 - scratch_412*scratch_512 - scratch_413*scratch_513 - scratch_414*scratch_514 - scratch_415*scratch_515 - scratch_416*scratch_516 - scratch_417*scratch_517 - scratch_418*scratch_518 - scratch_419*scratch_519 - scratch_458*(-scratch_34*scratch_425 + scratch_424 - scratch_426*scratch_507 - scratch_427*scratch_508 - scratch_428*scratch_509 - scratch_429*scratch_510 - scratch_430*scratch_511 - scratch_431*scratch_512 - scratch_432*scratch_513 - scratch_433*scratch_514 - scratch_434*scratch_515 - scratch_435*scratch_516 - scratch_436*scratch_517 - scratch_437*scratch_518 - scratch_438*scratch_520);
    x[5] = -scratch_381*(-scratch_109*scratch_523 - scratch_131*scratch_524 - scratch_155*scratch_525 - scratch_181*scratch_526 - scratch_209*scratch_527 - scratch_239*scratch_528 - scratch_271*scratch_529 - scratch_305*scratch_530 - scratch_341*scratch_531 - scratch_359*scratch_533 + scratch_41 - scratch_47*scratch_55 - scratch_521*scratch_71 - scratch_522*scratch_89) - scratch_401*(scratch_387 - scratch_388*scratch_47 - scratch_389*scratch_521 - scratch_390*scratch_522 - scratch_391*scratch_523 - scratch_392*scratch_524 - scratch_393*scratch_525 - scratch_394*scratch_526 - scratch_395*scratch_527 - scratch_396*scratch_528 - scratch_397*scratch_529 - scratch_398*scratch_530 - scratch_399*scratch_531 - scratch_400*scratch_533) + scratch_406 - scratch_407*scratch_47 - scratch_408*scratch_521 - scratch_409*scratch_522 - scratch_410*scratch_523 - scratch_411*scratch_524 - scratch_412*scratch_525 - scratch_413*scratch_526 - scratch_414*scratch_527 - scratch_415*scratch_528 - scratch_416*scratch_529 - scratch_417*scratch_530 - scratch_418*scratch_531 - scratch_419*scratch_532 - scratch_458*(scratch_425 - scratch_426*scratch_47 - scratch_427*scratch_521 - scratch_428*scratch_522 - scratch_429*scratch_523 - scratch_430*scratch_524 - scratch_431*scratch_525 - scratch_432*scratch_526 - scratch_433*scratch_527 - scratch_434*scratch_528 - scratch_435*scratch_529 - scratch_436*scratch_530 - scratch_437*scratch_531 - scratch_438*scratch_533);
    x[6] = -scratch_381*(-scratch_109*scratch_535 - scratch_131*scratch_536 - scratch_155*scratch_537 - scratch_181*scratch_538 - scratch_209*scratch_539 - scratch_239*scratch_540 - scratch_271*scratch_541 - scratch_305*scratch_542 - scratch_341*scratch_543 - scratch_359*scratch_545 - scratch_534*scratch_89 + scratch_55 - scratch_62*scratch_71) - scratch_401*(scratch_388 - scratch_389*scratch_62 - scratch_390*scratch_534 - scratch_391*scratch_535 - scratch_392*scratch_536 - scratch_393*scratch_537 - scratch_394*scratch_538 - scratch_395*scratch_539 - scratch_396*scratch_540 - scratch_397*scratch_541 - scratch_398*scratch_542 - scratch_399*scratch_543 - scratch_400*scratch_545) + scratch_407 - scratch_408*scratch_62 - scratch_409*scratch_534 - scratch_410*scratch_535 - scratch_411*scratch_536 - scratch_412*scratch_537 - scratch_413*scratch_538 - scratch_414*scratch_539 - scratch_415*scratch_540 - scratch_416*scratch_541 - scratch_417*scratch_542 - scratch_418*scratch_543 - scratch_419*scratch_544 - scratch_458*(scratch_426 - scratch_427*scratch_62 - scratch_428*scratch_534 - scratch_429*scratch_535 - scratch_430*scratch_536 - scratch_431*scratch_537 - scratch_432*scratch_538 - scratch_433*scratch_539 - scratch_434*scratch_540 - scratch_435*scratch_541 - scratch_436*scratch_542 - scratch_437*scratch_543 - scratch_438*scratch_545);
    x[7] = -scratch_381*(-scratch_109*scratch_546 - scratch_131*scratch_547 - scratch_155*scratch_548 - scratch_181*scratch_549 - scratch_209*scratch_550 - scratch_239*scratch_551 - scratch_271*scratch_552 - scratch_305*scratch_553 - scratch_341*scratch_554 - scratch_359*scratch_556 + scratch_71 - scratch_79*scratch_89) - scratch_401*(scratch_389 - scratch_390*scratch_79 - scratch_391*scratch_546 - scratch_392*scratch_547 - scratch_393*scratch_548 - scratch_394*scratch_549 - scratch_395*scratch_550 - scratch_396*scratch_551 - scratch_397*scratch_552 - scratch_398*scratch_553 - scratch_399*scratch_554 - scratch_400*scratch_556) + scratch_408 - scratch_409*scratch_79 - scratch_410*scratch_546 - scratch_411*scratch_547 - scratch_412*scratch_548 - scratch_413*scratch_549 - scratch_414*scratch_550 - scratch_415*scratch_551 - scratch_416*scratch_552 - scratch_417*scratch_553 - scratch_418*scratch_554 - scratch_419*scratch_555 - scratch_458*(scratch_427 - scratch_428*scratch_79 - scratch_429*scratch_546 - scratch_430*scratch_547 - scratch_431*scratch_548 - scratch_432*scratch_549 - scratch_433*scratch_550 - scratch_434*scratch_551 - scratch_435*scratch_552 - scratch_436*scratch_553 - scratch_437*scratch_554 - scratch_438*scratch_556);
    x[8] = -scratch_381*(-scratch_109*scratch_98 - scratch_131*scratch_557 - scratch_155*scratch_558 - scratch_181*scratch_559 - scratch_209*scratch_560 - scratch_239*scratch_561 - scratch_271*scratch_562 - scratch_305*scratch_563 - scratch_341*scratch_564 - scratch_359*scratch_566 + scratch_89) - scratch_401*(scratch_390 - scratch_391*scratch_98 - scratch_392*scratch_557 - scratch_393*scratch_558 - scratch_394*scratch_559 - scratch_395*scratch_560 - scratch_396*scratch_561 - scratch_397*scratch_562 - scratch_398*scratch_563 - scratch_399*scratch_564 - scratch_400*scratch_566) + scratch_409 - scratch_410*scratch_98 - scratch_411*scratch_557 - scratch_412*scratch_558 - scratch_413*scratch_559 - scratch_414*scratch_560 - scratch_415*scratch_561 - scratch_416*scratch_562 - scratch_417*scratch_563 - scratch_418*scratch_564 - scratch_419*scratch_565 - scratch_458*(scratch_428 - scratch_429*scratch_98 - scratch_430*scratch_557 - scratch_431*scratch_558 - scratch_432*scratch_559 - scratch_433*scratch_560 - scratch_434*scratch_561 - scratch_435*scratch_562 - scratch_436*scratch_563 - scratch_437*scratch_564 - scratch_438*scratch_566);
    x[9] = -scratch_119*scratch_411 - scratch_381*(scratch_109 - scratch_119*scratch_131 - scratch_155*scratch_567 - scratch_181*scratch_568 - scratch_209*scratch_569 - scratch_239*scratch_570 - scratch_271*scratch_571 - scratch_305*scratch_572 - scratch_341*scratch_573 - scratch_359*scratch_575) - scratch_401*(-scratch_119*scratch_392 + scratch_391 - scratch_393*scratch_567 - scratch_394*scratch_568 - scratch_395*scratch_569 - scratch_396*scratch_570 - scratch_397*scratch_571 - scratch_398*scratch_572 - scratch_399*scratch_573 - scratch_400*scratch_575) + scratch_410 - scratch_412*scratch_567 - scratch_413*scratch_568 - scratch_414*scratch_569 - scratch_415*scratch_570 - scratch_416*scratch_571 - scratch_417*scratch_572 - scratch_418*scratch_573 - scratch_419*scratch_574 - scratch_458*(-scratch_119*scratch_430 + scratch_429 - scratch_431*scratch_567 - scratch_432*scratch_568 - scratch_433*scratch_569 - scratch_434*scratch_570 - scratch_435*scratch_571 - scratch_436*scratch_572 - scratch_437*scratch_573 - scratch_438*scratch_575);
    x[10] = -scratch_142*scratch_412 - scratch_381*(scratch_131 - scratch_142*scratch_155 - scratch_181*scratch_576 - scratch_209*scratch_577 - scratch_239*scratch_578 - scratch_271*scratch_579 - scratch_305*scratch_580 - scratch_341*scratch_581 - scratch_359*scratch_583) - scratch_401*(-scratch_142*scratch_393 + scratch_392 - scratch_394*scratch_576 - scratch_395*scratch_577 - scratch_396*scratch_578 - scratch_397*scratch_579 - scratch_398*scratch_580 - scratch_399*scratch_581 - scratch_400*scratch_583) + scratch_411 - scratch_413*scratch_576 - scratch_414*scratch_577 - scratch_415*scratch_578 - scratch_416*scratch_579 - scratch_417*scratch_580 - scratch_418*scratch_581 - scratch_419*scratch_582 - scratch_458*(-scratch_142*scratch_431 + scratch_430 - scratch_432*scratch_576 - scratch_433*scratch_577 - scratch_434*scratch_578 - scratch_435*scratch_579 - scratch_436*scratch_580 - scratch_437*scratch_581 - scratch_438*scratch_583);
    x[11] = -scratch_167*scratch_413 - scratch_381*(scratch_155 - scratch_167*scratch_181 - scratch_209*scratch_584 - scratch_239*scratch_585 - scratch_271*scratch_586 - scratch_305*scratch_587 - scratch_341*scratch_588 - scratch_359*scratch_590) - scratch_401*(-scratch_167*scratch_394 + scratch_393 - scratch_395*scratch_584 - scratch_396*scratch_585 - scratch_397*scratch_586 - scratch_398*scratch_587 - scratch_399*scratch_588 - scratch_400*scratch_590) + scratch_412 - scratch_414*scratch_584 - scratch_415*scratch_585 - scratch_416*scratch_586 - scratch_417*scratch_587 - scratch_418*scratch_588 - scratch_419*scratch_589 - scratch_458*(-scratch_167*scratch_432 + scratch_431 - scratch_433*scratch_584 - scratch_434*scratch_585 - scratch_435*scratch_586 - scratch_436*scratch_587 - scratch_437*scratch_588 - scratch_438*scratch_590);
    x[12] = -scratch_194*scratch_414 - scratch_381*(scratch_181 - scratch_194*scratch_209 - scratch_239*scratch_591 - scratch_271*scratch_592 - scratch_305*scratch_593 - scratch_341*scratch_594 - scratch_359*scratch_596) - scratch_401*(-scratch_194*scratch_395 + scratch_394 - scratch_396*scratch_591 - scratch_397*scratch_592 - scratch_398*scratch_593 - scratch_399*scratch_594 - scratch_400*scratch_596) + scratch_413 - scratch_415*scratch_591 - scratch_416*scratch_592 - scratch_417*scratch_593 - scratch_418*scratch_594 - scratch_419*scratch_595 - scratch_458*(-scratch_194*scratch_433 + scratch_432 - scratch_434*scratch_591 - scratch_435*scratch_592 - scratch_436*scratch_593 - scratch_437*scratch_594 - scratch_438*scratch_596);
    x[13] = -scratch_223*scratch_415 - scratch_381*(scratch_209 - scratch_223*scratch_239 - scratch_271*scratch_597 - scratch_305*scratch_598 - scratch_341*scratch_599 - scratch_359*scratch_601) - scratch_401*(-scratch_223*scratch_396 + scratch_395 - scratch_397*scratch_597 - scratch_398*scratch_598 - scratch_399*scratch_599 - scratch_400*scratch_601) + scratch_414 - scratch_416*scratch_597 - scratch_417*scratch_598 - scratch_418*scratch_599 - scratch_419*scratch_600 - scratch_458*(-scratch_223*scratch_434 + scratch_433 - scratch_435*scratch_597 - scratch_436*scratch_598 - scratch_437*scratch_599 - scratch_438*scratch_601);
    x[14] = -scratch_254*scratch_416 - scratch_381*(scratch_239 - scratch_254*scratch_271 - scratch_305*scratch_602 - scratch_341*scratch_603 - scratch_359*scratch_605) - scratch_401*(-scratch_254*scratch_397 + scratch_396 - scratch_398*scratch_602 - scratch_399*scratch_603 - scratch_400*scratch_605) + scratch_415 - scratch_417*scratch_602 - scratch_418*scratch_603 - scratch_419*scratch_604 - scratch_458*(-scratch_254*scratch_435 + scratch_434 - scratch_436*scratch_602 - scratch_437*scratch_603 - scratch_438*scratch_605);
    x[15] = -scratch_287*scratch_417 - scratch_381*(scratch_271 - scratch_287*scratch_305 - scratch_341*scratch_606 - scratch_359*scratch_608) - scratch_401*(-scratch_287*scratch_398 + scratch_397 - scratch_399*scratch_606 - scratch_400*scratch_608) + scratch_416 - scratch_418*scratch_606 - scratch_419*scratch_607 - scratch_458*(-scratch_287*scratch_436 + scratch_435 - scratch_437*scratch_606 - scratch_438*scratch_608);
    x[16] = -scratch_322*scratch_418 - scratch_381*(scratch_305 - scratch_322*scratch_341 - scratch_359*scratch_610) - scratch_401*(-scratch_322*scratch_399 + scratch_398 - scratch_400*scratch_610) + scratch_417 - scratch_419*scratch_609 - scratch_458*(-scratch_322*scratch_437 + scratch_436 - scratch_438*scratch_610);
    x[17] = -scratch_377*scratch_419 - scratch_381*(scratch_341 - scratch_359*scratch_611) - scratch_401*(scratch_399 - scratch_400*scratch_611) + scratch_418 - scratch_458*(scratch_437 - scratch_438*scratch_611);
    x[18] = -scratch_359*scratch_378*scratch_381 - scratch_378*scratch_400*scratch_401 - scratch_378*scratch_438*scratch_458 + scratch_419;
    x[19] = scratch_381;
    x[20] = scratch_401;
    x[21] = scratch_458;
}

#endif
