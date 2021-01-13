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

    extern amrex::Real fwd_A[167], fwd_beta[167], fwd_Ea[167];
    extern amrex::Real low_A[167], low_beta[167], low_Ea[167];
    extern amrex::Real rev_A[167], rev_beta[167], rev_Ea[167];
    extern amrex::Real troe_a[167],troe_Ts[167], troe_Tss[167], troe_Tsss[167];
    extern amrex::Real sri_a[167], sri_b[167], sri_c[167], sri_d[167], sri_e[167];
    extern amrex::Real activation_units[167], prefactor_units[167], phase_units[167];
    extern int is_PD[167], troe_len[167], sri_len[167], nTB[167], *TBid[167];
    extern amrex::Real *TB[167];
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
void vcomp_wdot_51_100(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
                amrex::Real *  k_f_s, amrex::Real *  Kc_s,
                amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T);
void vcomp_wdot_101_150(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
                amrex::Real *  k_f_s, amrex::Real *  Kc_s,
                amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T);
void vcomp_wdot_151_167(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKINDX(int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 4;
    *kk = 34;
    *ii = 167;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/* inverse molecular weights  */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
    imw_new[0] = 1.0/28.013400;/*N2 */
    imw_new[1] = 1.0/15.999400;/*O */
    imw_new[2] = 1.0/1.007970;/*H */
    imw_new[3] = 1.0/17.007370;/*OH */
    imw_new[4] = 1.0/2.015940;/*H2 */
    imw_new[5] = 1.0/31.998800;/*O2 */
    imw_new[6] = 1.0/18.015340;/*H2O */
    imw_new[7] = 1.0/34.014740;/*H2O2 */
    imw_new[8] = 1.0/33.006770;/*HO2 */
    imw_new[9] = 1.0/14.027090;/*CH2GSG */
    imw_new[10] = 1.0/30.026490;/*CH2O */
    imw_new[11] = 1.0/15.035060;/*CH3 */
    imw_new[12] = 1.0/28.010550;/*CO */
    imw_new[13] = 1.0/31.034460;/*CH3O */
    imw_new[14] = 1.0/29.062150;/*C2H5 */
    imw_new[15] = 1.0/16.043030;/*CH4 */
    imw_new[16] = 1.0/28.054180;/*C2H4 */
    imw_new[17] = 1.0/30.070120;/*C2H6 */
    imw_new[18] = 1.0/44.009950;/*CO2 */
    imw_new[19] = 1.0/29.018520;/*HCO */
    imw_new[20] = 1.0/47.033860;/*CH3O2 */
    imw_new[21] = 1.0/48.041830;/*CH3O2H */
    imw_new[22] = 1.0/26.038240;/*C2H2 */
    imw_new[23] = 1.0/41.029670;/*HCCO */
    imw_new[24] = 1.0/27.046210;/*C2H3 */
    imw_new[25] = 1.0/43.045610;/*CH2CHO */
    imw_new[26] = 1.0/41.073300;/*C3H5XA */
    imw_new[27] = 1.0/42.081270;/*C3H6 */
    imw_new[28] = 1.0/57.072700;/*C3H5O */
    imw_new[29] = 1.0/43.089240;/*IXC3H7 */
    imw_new[30] = 1.0/43.089240;/*NXC3H7 */
    imw_new[31] = 1.0/44.097210;/*C3H8 */
    imw_new[32] = 1.0/75.088040;/*IXC3H7O2 */
    imw_new[33] = 1.0/75.088040;/*NXC3H7O2 */
}

/* molecular weights  */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void get_mw(amrex::Real *mw_new){
    mw_new[0] = 28.013400;/*N2 */
    mw_new[1] = 15.999400;/*O */
    mw_new[2] = 1.007970;/*H */
    mw_new[3] = 17.007370;/*OH */
    mw_new[4] = 2.015940;/*H2 */
    mw_new[5] = 31.998800;/*O2 */
    mw_new[6] = 18.015340;/*H2O */
    mw_new[7] = 34.014740;/*H2O2 */
    mw_new[8] = 33.006770;/*HO2 */
    mw_new[9] = 14.027090;/*CH2GSG */
    mw_new[10] = 30.026490;/*CH2O */
    mw_new[11] = 15.035060;/*CH3 */
    mw_new[12] = 28.010550;/*CO */
    mw_new[13] = 31.034460;/*CH3O */
    mw_new[14] = 29.062150;/*C2H5 */
    mw_new[15] = 16.043030;/*CH4 */
    mw_new[16] = 28.054180;/*C2H4 */
    mw_new[17] = 30.070120;/*C2H6 */
    mw_new[18] = 44.009950;/*CO2 */
    mw_new[19] = 29.018520;/*HCO */
    mw_new[20] = 47.033860;/*CH3O2 */
    mw_new[21] = 48.041830;/*CH3O2H */
    mw_new[22] = 26.038240;/*C2H2 */
    mw_new[23] = 41.029670;/*HCCO */
    mw_new[24] = 27.046210;/*C2H3 */
    mw_new[25] = 43.045610;/*CH2CHO */
    mw_new[26] = 41.073300;/*C3H5XA */
    mw_new[27] = 42.081270;/*C3H6 */
    mw_new[28] = 57.072700;/*C3H5O */
    mw_new[29] = 43.089240;/*IXC3H7 */
    mw_new[30] = 43.089240;/*NXC3H7 */
    mw_new[31] = 44.097210;/*C3H8 */
    mw_new[32] = 75.088040;/*IXC3H7O2 */
    mw_new[33] = 75.088040;/*NXC3H7O2 */
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
    for (int id = 0; id < 34 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    nuki[ 3 * kd + 0 ] += -2 ;
    nuki[ 7 * kd + 0 ] += +1 ;

    /*reaction 2: H + O2 (+M) <=> HO2 (+M) */
    nuki[ 2 * kd + 1 ] += -1 ;
    nuki[ 5 * kd + 1 ] += -1 ;
    nuki[ 8 * kd + 1 ] += +1 ;

    /*reaction 3: CH3 + H (+M) => CH4 (+M) */
    nuki[ 11 * kd + 2 ] += -1 ;
    nuki[ 2 * kd + 2 ] += -1 ;
    nuki[ 15 * kd + 2 ] += +1 ;

    /*reaction 4: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    nuki[ 11 * kd + 3 ] += -2 ;
    nuki[ 17 * kd + 3 ] += +1 ;

    /*reaction 5: CO + O (+M) => CO2 (+M) */
    nuki[ 12 * kd + 4 ] += -1 ;
    nuki[ 1 * kd + 4 ] += -1 ;
    nuki[ 18 * kd + 4 ] += +1 ;

    /*reaction 6: CH3O (+M) => CH2O + H (+M) */
    nuki[ 13 * kd + 5 ] += -1 ;
    nuki[ 10 * kd + 5 ] += +1 ;
    nuki[ 2 * kd + 5 ] += +1 ;

    /*reaction 7: C2H3 (+M) => H + C2H2 (+M) */
    nuki[ 24 * kd + 6 ] += -1 ;
    nuki[ 2 * kd + 6 ] += +1 ;
    nuki[ 22 * kd + 6 ] += +1 ;

    /*reaction 8: H + C2H4 (+M) <=> C2H5 (+M) */
    nuki[ 2 * kd + 7 ] += -1 ;
    nuki[ 16 * kd + 7 ] += -1 ;
    nuki[ 14 * kd + 7 ] += +1 ;

    /*reaction 9: C2H4 (+M) => C2H2 + H2 (+M) */
    nuki[ 16 * kd + 8 ] += -1 ;
    nuki[ 22 * kd + 8 ] += +1 ;
    nuki[ 4 * kd + 8 ] += +1 ;

    /*reaction 10: O + H + M => OH + M */
    nuki[ 1 * kd + 9 ] += -1 ;
    nuki[ 2 * kd + 9 ] += -1 ;
    nuki[ 3 * kd + 9 ] += +1 ;

    /*reaction 11: 2.000000 O + M => O2 + M */
    nuki[ 1 * kd + 10 ] += -2 ;
    nuki[ 5 * kd + 10 ] += +1 ;

    /*reaction 12: H + OH + M => H2O + M */
    nuki[ 2 * kd + 11 ] += -1 ;
    nuki[ 3 * kd + 11 ] += -1 ;
    nuki[ 6 * kd + 11 ] += +1 ;

    /*reaction 13: HCO + M => H + CO + M */
    nuki[ 19 * kd + 12 ] += -1 ;
    nuki[ 2 * kd + 12 ] += +1 ;
    nuki[ 12 * kd + 12 ] += +1 ;

    /*reaction 14: H + CO + M => HCO + M */
    nuki[ 2 * kd + 13 ] += -1 ;
    nuki[ 12 * kd + 13 ] += -1 ;
    nuki[ 19 * kd + 13 ] += +1 ;

    /*reaction 15: CH3O2 + M => CH3 + O2 + M */
    nuki[ 20 * kd + 14 ] += -1 ;
    nuki[ 11 * kd + 14 ] += +1 ;
    nuki[ 5 * kd + 14 ] += +1 ;

    /*reaction 16: CH3 + O2 + M => CH3O2 + M */
    nuki[ 11 * kd + 15 ] += -1 ;
    nuki[ 5 * kd + 15 ] += -1 ;
    nuki[ 20 * kd + 15 ] += +1 ;

    /*reaction 17: O + H2 => H + OH */
    nuki[ 1 * kd + 16 ] += -1 ;
    nuki[ 4 * kd + 16 ] += -1 ;
    nuki[ 2 * kd + 16 ] += +1 ;
    nuki[ 3 * kd + 16 ] += +1 ;

    /*reaction 18: H + OH => O + H2 */
    nuki[ 2 * kd + 17 ] += -1 ;
    nuki[ 3 * kd + 17 ] += -1 ;
    nuki[ 1 * kd + 17 ] += +1 ;
    nuki[ 4 * kd + 17 ] += +1 ;

    /*reaction 19: OH + H2 => H + H2O */
    nuki[ 3 * kd + 18 ] += -1 ;
    nuki[ 4 * kd + 18 ] += -1 ;
    nuki[ 2 * kd + 18 ] += +1 ;
    nuki[ 6 * kd + 18 ] += +1 ;

    /*reaction 20: H + H2O => OH + H2 */
    nuki[ 2 * kd + 19 ] += -1 ;
    nuki[ 6 * kd + 19 ] += -1 ;
    nuki[ 3 * kd + 19 ] += +1 ;
    nuki[ 4 * kd + 19 ] += +1 ;

    /*reaction 21: O + H2O => 2.000000 OH */
    nuki[ 1 * kd + 20 ] += -1 ;
    nuki[ 6 * kd + 20 ] += -1 ;
    nuki[ 3 * kd + 20 ] += +2 ;

    /*reaction 22: 2.000000 OH => O + H2O */
    nuki[ 3 * kd + 21 ] += -2 ;
    nuki[ 1 * kd + 21 ] += +1 ;
    nuki[ 6 * kd + 21 ] += +1 ;

    /*reaction 23: H + O2 => O + OH */
    nuki[ 2 * kd + 22 ] += -1 ;
    nuki[ 5 * kd + 22 ] += -1 ;
    nuki[ 1 * kd + 22 ] += +1 ;
    nuki[ 3 * kd + 22 ] += +1 ;

    /*reaction 24: O + OH => H + O2 */
    nuki[ 1 * kd + 23 ] += -1 ;
    nuki[ 3 * kd + 23 ] += -1 ;
    nuki[ 2 * kd + 23 ] += +1 ;
    nuki[ 5 * kd + 23 ] += +1 ;

    /*reaction 25: HO2 + H => 2.000000 OH */
    nuki[ 8 * kd + 24 ] += -1 ;
    nuki[ 2 * kd + 24 ] += -1 ;
    nuki[ 3 * kd + 24 ] += +2 ;

    /*reaction 26: 2.000000 HO2 => H2O2 + O2 */
    nuki[ 8 * kd + 25 ] += -2 ;
    nuki[ 7 * kd + 25 ] += +1 ;
    nuki[ 5 * kd + 25 ] += +1 ;

    /*reaction 27: HO2 + H => H2 + O2 */
    nuki[ 8 * kd + 26 ] += -1 ;
    nuki[ 2 * kd + 26 ] += -1 ;
    nuki[ 4 * kd + 26 ] += +1 ;
    nuki[ 5 * kd + 26 ] += +1 ;

    /*reaction 28: HO2 + OH => H2O + O2 */
    nuki[ 8 * kd + 27 ] += -1 ;
    nuki[ 3 * kd + 27 ] += -1 ;
    nuki[ 6 * kd + 27 ] += +1 ;
    nuki[ 5 * kd + 27 ] += +1 ;

    /*reaction 29: H2O + O2 => HO2 + OH */
    nuki[ 6 * kd + 28 ] += -1 ;
    nuki[ 5 * kd + 28 ] += -1 ;
    nuki[ 8 * kd + 28 ] += +1 ;
    nuki[ 3 * kd + 28 ] += +1 ;

    /*reaction 30: 2.000000 HO2 => H2O2 + O2 */
    nuki[ 8 * kd + 29 ] += -2 ;
    nuki[ 7 * kd + 29 ] += +1 ;
    nuki[ 5 * kd + 29 ] += +1 ;

    /*reaction 31: HO2 + O => OH + O2 */
    nuki[ 8 * kd + 30 ] += -1 ;
    nuki[ 1 * kd + 30 ] += -1 ;
    nuki[ 3 * kd + 30 ] += +1 ;
    nuki[ 5 * kd + 30 ] += +1 ;

    /*reaction 32: OH + O2 => HO2 + O */
    nuki[ 3 * kd + 31 ] += -1 ;
    nuki[ 5 * kd + 31 ] += -1 ;
    nuki[ 8 * kd + 31 ] += +1 ;
    nuki[ 1 * kd + 31 ] += +1 ;

    /*reaction 33: H2O2 + OH => H2O + HO2 */
    nuki[ 7 * kd + 32 ] += -1 ;
    nuki[ 3 * kd + 32 ] += -1 ;
    nuki[ 6 * kd + 32 ] += +1 ;
    nuki[ 8 * kd + 32 ] += +1 ;

    /*reaction 34: H2O2 + H => H2O + OH */
    nuki[ 7 * kd + 33 ] += -1 ;
    nuki[ 2 * kd + 33 ] += -1 ;
    nuki[ 6 * kd + 33 ] += +1 ;
    nuki[ 3 * kd + 33 ] += +1 ;

    /*reaction 35: H2O2 + OH => H2O + HO2 */
    nuki[ 7 * kd + 34 ] += -1 ;
    nuki[ 3 * kd + 34 ] += -1 ;
    nuki[ 6 * kd + 34 ] += +1 ;
    nuki[ 8 * kd + 34 ] += +1 ;

    /*reaction 36: H2O + HO2 => H2O2 + OH */
    nuki[ 6 * kd + 35 ] += -1 ;
    nuki[ 8 * kd + 35 ] += -1 ;
    nuki[ 7 * kd + 35 ] += +1 ;
    nuki[ 3 * kd + 35 ] += +1 ;

    /*reaction 37: H2O2 + O => OH + HO2 */
    nuki[ 7 * kd + 36 ] += -1 ;
    nuki[ 1 * kd + 36 ] += -1 ;
    nuki[ 3 * kd + 36 ] += +1 ;
    nuki[ 8 * kd + 36 ] += +1 ;

    /*reaction 38: H2O2 + H => H2 + HO2 */
    nuki[ 7 * kd + 37 ] += -1 ;
    nuki[ 2 * kd + 37 ] += -1 ;
    nuki[ 4 * kd + 37 ] += +1 ;
    nuki[ 8 * kd + 37 ] += +1 ;

    /*reaction 39: H2 + HO2 => H2O2 + H */
    nuki[ 4 * kd + 38 ] += -1 ;
    nuki[ 8 * kd + 38 ] += -1 ;
    nuki[ 7 * kd + 38 ] += +1 ;
    nuki[ 2 * kd + 38 ] += +1 ;

    /*reaction 40: CH2GSG + OH => CH2O + H */
    nuki[ 9 * kd + 39 ] += -1 ;
    nuki[ 3 * kd + 39 ] += -1 ;
    nuki[ 10 * kd + 39 ] += +1 ;
    nuki[ 2 * kd + 39 ] += +1 ;

    /*reaction 41: CH2GSG + H2 => CH3 + H */
    nuki[ 9 * kd + 40 ] += -1 ;
    nuki[ 4 * kd + 40 ] += -1 ;
    nuki[ 11 * kd + 40 ] += +1 ;
    nuki[ 2 * kd + 40 ] += +1 ;

    /*reaction 42: CH3 + H => CH2GSG + H2 */
    nuki[ 11 * kd + 41 ] += -1 ;
    nuki[ 2 * kd + 41 ] += -1 ;
    nuki[ 9 * kd + 41 ] += +1 ;
    nuki[ 4 * kd + 41 ] += +1 ;

    /*reaction 43: CH2GSG + O2 => CO + OH + H */
    nuki[ 9 * kd + 42 ] += -1 ;
    nuki[ 5 * kd + 42 ] += -1 ;
    nuki[ 12 * kd + 42 ] += +1 ;
    nuki[ 3 * kd + 42 ] += +1 ;
    nuki[ 2 * kd + 42 ] += +1 ;

    /*reaction 44: CH2GSG + O => CO + 2.000000 H */
    nuki[ 9 * kd + 43 ] += -1 ;
    nuki[ 1 * kd + 43 ] += -1 ;
    nuki[ 12 * kd + 43 ] += +1 ;
    nuki[ 2 * kd + 43 ] += +2 ;

    /*reaction 45: CH3 + HO2 => CH3O + OH */
    nuki[ 11 * kd + 44 ] += -1 ;
    nuki[ 8 * kd + 44 ] += -1 ;
    nuki[ 13 * kd + 44 ] += +1 ;
    nuki[ 3 * kd + 44 ] += +1 ;

    /*reaction 46: CH3 + O2 => CH2O + OH */
    nuki[ 11 * kd + 45 ] += -1 ;
    nuki[ 5 * kd + 45 ] += -1 ;
    nuki[ 10 * kd + 45 ] += +1 ;
    nuki[ 3 * kd + 45 ] += +1 ;

    /*reaction 47: CH3 + O2 => CH3O + O */
    nuki[ 11 * kd + 46 ] += -1 ;
    nuki[ 5 * kd + 46 ] += -1 ;
    nuki[ 13 * kd + 46 ] += +1 ;
    nuki[ 1 * kd + 46 ] += +1 ;

    /*reaction 48: CH3O + O => CH3 + O2 */
    nuki[ 13 * kd + 47 ] += -1 ;
    nuki[ 1 * kd + 47 ] += -1 ;
    nuki[ 11 * kd + 47 ] += +1 ;
    nuki[ 5 * kd + 47 ] += +1 ;

    /*reaction 49: 2.000000 CH3 <=> H + C2H5 */
    nuki[ 11 * kd + 48 ] += -2 ;
    nuki[ 2 * kd + 48 ] += +1 ;
    nuki[ 14 * kd + 48 ] += +1 ;

    /*reaction 50: CH3 + HO2 => CH4 + O2 */
    nuki[ 11 * kd + 49 ] += -1 ;
    nuki[ 8 * kd + 49 ] += -1 ;
    nuki[ 15 * kd + 49 ] += +1 ;
    nuki[ 5 * kd + 49 ] += +1 ;

    /*reaction 51: CH3 + O => CH2O + H */
    nuki[ 11 * kd + 50 ] += -1 ;
    nuki[ 1 * kd + 50 ] += -1 ;
    nuki[ 10 * kd + 50 ] += +1 ;
    nuki[ 2 * kd + 50 ] += +1 ;

    /*reaction 52: CH3 + OH => CH2O + H2 */
    nuki[ 11 * kd + 51 ] += -1 ;
    nuki[ 3 * kd + 51 ] += -1 ;
    nuki[ 10 * kd + 51 ] += +1 ;
    nuki[ 4 * kd + 51 ] += +1 ;

    /*reaction 53: CH3 + OH => CH2GSG + H2O */
    nuki[ 11 * kd + 52 ] += -1 ;
    nuki[ 3 * kd + 52 ] += -1 ;
    nuki[ 9 * kd + 52 ] += +1 ;
    nuki[ 6 * kd + 52 ] += +1 ;

    /*reaction 54: CH2GSG + H2O => CH3 + OH */
    nuki[ 9 * kd + 53 ] += -1 ;
    nuki[ 6 * kd + 53 ] += -1 ;
    nuki[ 11 * kd + 53 ] += +1 ;
    nuki[ 3 * kd + 53 ] += +1 ;

    /*reaction 55: CH3 + H2O2 => CH4 + HO2 */
    nuki[ 11 * kd + 54 ] += -1 ;
    nuki[ 7 * kd + 54 ] += -1 ;
    nuki[ 15 * kd + 54 ] += +1 ;
    nuki[ 8 * kd + 54 ] += +1 ;

    /*reaction 56: CH2GSG + CH3 => C2H4 + H */
    nuki[ 9 * kd + 55 ] += -1 ;
    nuki[ 11 * kd + 55 ] += -1 ;
    nuki[ 16 * kd + 55 ] += +1 ;
    nuki[ 2 * kd + 55 ] += +1 ;

    /*reaction 57: CH2GSG + CH4 => 2.000000 CH3 */
    nuki[ 9 * kd + 56 ] += -1 ;
    nuki[ 15 * kd + 56 ] += -1 ;
    nuki[ 11 * kd + 56 ] += +2 ;

    /*reaction 58: 2.000000 CH3 => CH2GSG + CH4 */
    nuki[ 11 * kd + 57 ] += -2 ;
    nuki[ 9 * kd + 57 ] += +1 ;
    nuki[ 15 * kd + 57 ] += +1 ;

    /*reaction 59: CH4 + O => CH3 + OH */
    nuki[ 15 * kd + 58 ] += -1 ;
    nuki[ 1 * kd + 58 ] += -1 ;
    nuki[ 11 * kd + 58 ] += +1 ;
    nuki[ 3 * kd + 58 ] += +1 ;

    /*reaction 60: CH3 + OH => CH4 + O */
    nuki[ 11 * kd + 59 ] += -1 ;
    nuki[ 3 * kd + 59 ] += -1 ;
    nuki[ 15 * kd + 59 ] += +1 ;
    nuki[ 1 * kd + 59 ] += +1 ;

    /*reaction 61: CH4 + H => CH3 + H2 */
    nuki[ 15 * kd + 60 ] += -1 ;
    nuki[ 2 * kd + 60 ] += -1 ;
    nuki[ 11 * kd + 60 ] += +1 ;
    nuki[ 4 * kd + 60 ] += +1 ;

    /*reaction 62: CH3 + H2 => CH4 + H */
    nuki[ 11 * kd + 61 ] += -1 ;
    nuki[ 4 * kd + 61 ] += -1 ;
    nuki[ 15 * kd + 61 ] += +1 ;
    nuki[ 2 * kd + 61 ] += +1 ;

    /*reaction 63: CH4 + OH => CH3 + H2O */
    nuki[ 15 * kd + 62 ] += -1 ;
    nuki[ 3 * kd + 62 ] += -1 ;
    nuki[ 11 * kd + 62 ] += +1 ;
    nuki[ 6 * kd + 62 ] += +1 ;

    /*reaction 64: CH3 + H2O => CH4 + OH */
    nuki[ 11 * kd + 63 ] += -1 ;
    nuki[ 6 * kd + 63 ] += -1 ;
    nuki[ 15 * kd + 63 ] += +1 ;
    nuki[ 3 * kd + 63 ] += +1 ;

    /*reaction 65: CO + HO2 => CO2 + OH */
    nuki[ 12 * kd + 64 ] += -1 ;
    nuki[ 8 * kd + 64 ] += -1 ;
    nuki[ 18 * kd + 64 ] += +1 ;
    nuki[ 3 * kd + 64 ] += +1 ;

    /*reaction 66: CO + O2 => CO2 + O */
    nuki[ 12 * kd + 65 ] += -1 ;
    nuki[ 5 * kd + 65 ] += -1 ;
    nuki[ 18 * kd + 65 ] += +1 ;
    nuki[ 1 * kd + 65 ] += +1 ;

    /*reaction 67: CO2 + O => CO + O2 */
    nuki[ 18 * kd + 66 ] += -1 ;
    nuki[ 1 * kd + 66 ] += -1 ;
    nuki[ 12 * kd + 66 ] += +1 ;
    nuki[ 5 * kd + 66 ] += +1 ;

    /*reaction 68: CO + OH => CO2 + H */
    nuki[ 12 * kd + 67 ] += -1 ;
    nuki[ 3 * kd + 67 ] += -1 ;
    nuki[ 18 * kd + 67 ] += +1 ;
    nuki[ 2 * kd + 67 ] += +1 ;

    /*reaction 69: CO2 + H => CO + OH */
    nuki[ 18 * kd + 68 ] += -1 ;
    nuki[ 2 * kd + 68 ] += -1 ;
    nuki[ 12 * kd + 68 ] += +1 ;
    nuki[ 3 * kd + 68 ] += +1 ;

    /*reaction 70: HCO + CH3 => CH4 + CO */
    nuki[ 19 * kd + 69 ] += -1 ;
    nuki[ 11 * kd + 69 ] += -1 ;
    nuki[ 15 * kd + 69 ] += +1 ;
    nuki[ 12 * kd + 69 ] += +1 ;

    /*reaction 71: HCO + H => CO + H2 */
    nuki[ 19 * kd + 70 ] += -1 ;
    nuki[ 2 * kd + 70 ] += -1 ;
    nuki[ 12 * kd + 70 ] += +1 ;
    nuki[ 4 * kd + 70 ] += +1 ;

    /*reaction 72: HCO + O2 => CO + HO2 */
    nuki[ 19 * kd + 71 ] += -1 ;
    nuki[ 5 * kd + 71 ] += -1 ;
    nuki[ 12 * kd + 71 ] += +1 ;
    nuki[ 8 * kd + 71 ] += +1 ;

    /*reaction 73: HCO + O => CO + OH */
    nuki[ 19 * kd + 72 ] += -1 ;
    nuki[ 1 * kd + 72 ] += -1 ;
    nuki[ 12 * kd + 72 ] += +1 ;
    nuki[ 3 * kd + 72 ] += +1 ;

    /*reaction 74: HCO + O => CO2 + H */
    nuki[ 19 * kd + 73 ] += -1 ;
    nuki[ 1 * kd + 73 ] += -1 ;
    nuki[ 18 * kd + 73 ] += +1 ;
    nuki[ 2 * kd + 73 ] += +1 ;

    /*reaction 75: HCO + OH => CO + H2O */
    nuki[ 19 * kd + 74 ] += -1 ;
    nuki[ 3 * kd + 74 ] += -1 ;
    nuki[ 12 * kd + 74 ] += +1 ;
    nuki[ 6 * kd + 74 ] += +1 ;

    /*reaction 76: CH2O + H => HCO + H2 */
    nuki[ 10 * kd + 75 ] += -1 ;
    nuki[ 2 * kd + 75 ] += -1 ;
    nuki[ 19 * kd + 75 ] += +1 ;
    nuki[ 4 * kd + 75 ] += +1 ;

    /*reaction 77: CH2O + O2 => HCO + HO2 */
    nuki[ 10 * kd + 76 ] += -1 ;
    nuki[ 5 * kd + 76 ] += -1 ;
    nuki[ 19 * kd + 76 ] += +1 ;
    nuki[ 8 * kd + 76 ] += +1 ;

    /*reaction 78: CH2O + OH => HCO + H2O */
    nuki[ 10 * kd + 77 ] += -1 ;
    nuki[ 3 * kd + 77 ] += -1 ;
    nuki[ 19 * kd + 77 ] += +1 ;
    nuki[ 6 * kd + 77 ] += +1 ;

    /*reaction 79: CH2O + HO2 => HCO + H2O2 */
    nuki[ 10 * kd + 78 ] += -1 ;
    nuki[ 8 * kd + 78 ] += -1 ;
    nuki[ 19 * kd + 78 ] += +1 ;
    nuki[ 7 * kd + 78 ] += +1 ;

    /*reaction 80: CH2O + O => HCO + OH */
    nuki[ 10 * kd + 79 ] += -1 ;
    nuki[ 1 * kd + 79 ] += -1 ;
    nuki[ 19 * kd + 79 ] += +1 ;
    nuki[ 3 * kd + 79 ] += +1 ;

    /*reaction 81: CH2O + CH3 => HCO + CH4 */
    nuki[ 10 * kd + 80 ] += -1 ;
    nuki[ 11 * kd + 80 ] += -1 ;
    nuki[ 19 * kd + 80 ] += +1 ;
    nuki[ 15 * kd + 80 ] += +1 ;

    /*reaction 82: CH3O + O2 => CH2O + HO2 */
    nuki[ 13 * kd + 81 ] += -1 ;
    nuki[ 5 * kd + 81 ] += -1 ;
    nuki[ 10 * kd + 81 ] += +1 ;
    nuki[ 8 * kd + 81 ] += +1 ;

    /*reaction 83: CH2GSG + CO2 => CH2O + CO */
    nuki[ 9 * kd + 82 ] += -1 ;
    nuki[ 18 * kd + 82 ] += -1 ;
    nuki[ 10 * kd + 82 ] += +1 ;
    nuki[ 12 * kd + 82 ] += +1 ;

    /*reaction 84: CH3O2 + CH3 => 2.000000 CH3O */
    nuki[ 20 * kd + 83 ] += -1 ;
    nuki[ 11 * kd + 83 ] += -1 ;
    nuki[ 13 * kd + 83 ] += +2 ;

    /*reaction 85: CH3O2 + HO2 => CH3O2H + O2 */
    nuki[ 20 * kd + 84 ] += -1 ;
    nuki[ 8 * kd + 84 ] += -1 ;
    nuki[ 21 * kd + 84 ] += +1 ;
    nuki[ 5 * kd + 84 ] += +1 ;

    /*reaction 86: 2.000000 CH3O2 => O2 + 2.000000 CH3O */
    nuki[ 20 * kd + 85 ] += -2 ;
    nuki[ 5 * kd + 85 ] += +1 ;
    nuki[ 13 * kd + 85 ] += +2 ;

    /*reaction 87: CH3O2 + CH2O => CH3O2H + HCO */
    nuki[ 20 * kd + 86 ] += -1 ;
    nuki[ 10 * kd + 86 ] += -1 ;
    nuki[ 21 * kd + 86 ] += +1 ;
    nuki[ 19 * kd + 86 ] += +1 ;

    /*reaction 88: CH3O2H => CH3O + OH */
    nuki[ 21 * kd + 87 ] += -1 ;
    nuki[ 13 * kd + 87 ] += +1 ;
    nuki[ 3 * kd + 87 ] += +1 ;

    /*reaction 89: CH3O + OH => CH3O2H */
    nuki[ 13 * kd + 88 ] += -1 ;
    nuki[ 3 * kd + 88 ] += -1 ;
    nuki[ 21 * kd + 88 ] += +1 ;

    /*reaction 90: C2H2 + O2 => HCCO + OH */
    nuki[ 22 * kd + 89 ] += -1 ;
    nuki[ 5 * kd + 89 ] += -1 ;
    nuki[ 23 * kd + 89 ] += +1 ;
    nuki[ 3 * kd + 89 ] += +1 ;

    /*reaction 91: C2H2 + O => HCCO + H */
    nuki[ 22 * kd + 90 ] += -1 ;
    nuki[ 1 * kd + 90 ] += -1 ;
    nuki[ 23 * kd + 90 ] += +1 ;
    nuki[ 2 * kd + 90 ] += +1 ;

    /*reaction 92: C2H3 + O2 => CH2CHO + O */
    nuki[ 24 * kd + 91 ] += -1 ;
    nuki[ 5 * kd + 91 ] += -1 ;
    nuki[ 25 * kd + 91 ] += +1 ;
    nuki[ 1 * kd + 91 ] += +1 ;

    /*reaction 93: C2H3 + H => C2H2 + H2 */
    nuki[ 24 * kd + 92 ] += -1 ;
    nuki[ 2 * kd + 92 ] += -1 ;
    nuki[ 22 * kd + 92 ] += +1 ;
    nuki[ 4 * kd + 92 ] += +1 ;

    /*reaction 94: C2H3 + O2 => CH2O + HCO */
    nuki[ 24 * kd + 93 ] += -1 ;
    nuki[ 5 * kd + 93 ] += -1 ;
    nuki[ 10 * kd + 93 ] += +1 ;
    nuki[ 19 * kd + 93 ] += +1 ;

    /*reaction 95: C2H3 + O2 => C2H2 + HO2 */
    nuki[ 24 * kd + 94 ] += -1 ;
    nuki[ 5 * kd + 94 ] += -1 ;
    nuki[ 22 * kd + 94 ] += +1 ;
    nuki[ 8 * kd + 94 ] += +1 ;

    /*reaction 96: C2H4 + O => CH3 + HCO */
    nuki[ 16 * kd + 95 ] += -1 ;
    nuki[ 1 * kd + 95 ] += -1 ;
    nuki[ 11 * kd + 95 ] += +1 ;
    nuki[ 19 * kd + 95 ] += +1 ;

    /*reaction 97: C2H4 + OH => C2H3 + H2O */
    nuki[ 16 * kd + 96 ] += -1 ;
    nuki[ 3 * kd + 96 ] += -1 ;
    nuki[ 24 * kd + 96 ] += +1 ;
    nuki[ 6 * kd + 96 ] += +1 ;

    /*reaction 98: C2H3 + H2O => C2H4 + OH */
    nuki[ 24 * kd + 97 ] += -1 ;
    nuki[ 6 * kd + 97 ] += -1 ;
    nuki[ 16 * kd + 97 ] += +1 ;
    nuki[ 3 * kd + 97 ] += +1 ;

    /*reaction 99: C2H4 + H => C2H3 + H2 */
    nuki[ 16 * kd + 98 ] += -1 ;
    nuki[ 2 * kd + 98 ] += -1 ;
    nuki[ 24 * kd + 98 ] += +1 ;
    nuki[ 4 * kd + 98 ] += +1 ;

    /*reaction 100: C2H3 + H2 => C2H4 + H */
    nuki[ 24 * kd + 99 ] += -1 ;
    nuki[ 4 * kd + 99 ] += -1 ;
    nuki[ 16 * kd + 99 ] += +1 ;
    nuki[ 2 * kd + 99 ] += +1 ;

    /*reaction 101: C2H4 + CH3O2 => C2H3 + CH3O2H */
    nuki[ 16 * kd + 100 ] += -1 ;
    nuki[ 20 * kd + 100 ] += -1 ;
    nuki[ 24 * kd + 100 ] += +1 ;
    nuki[ 21 * kd + 100 ] += +1 ;

    /*reaction 102: C2H4 + CH3 => C2H3 + CH4 */
    nuki[ 16 * kd + 101 ] += -1 ;
    nuki[ 11 * kd + 101 ] += -1 ;
    nuki[ 24 * kd + 101 ] += +1 ;
    nuki[ 15 * kd + 101 ] += +1 ;

    /*reaction 103: C2H3 + CH4 => C2H4 + CH3 */
    nuki[ 24 * kd + 102 ] += -1 ;
    nuki[ 15 * kd + 102 ] += -1 ;
    nuki[ 16 * kd + 102 ] += +1 ;
    nuki[ 11 * kd + 102 ] += +1 ;

    /*reaction 104: C2H4 + O => CH2CHO + H */
    nuki[ 16 * kd + 103 ] += -1 ;
    nuki[ 1 * kd + 103 ] += -1 ;
    nuki[ 25 * kd + 103 ] += +1 ;
    nuki[ 2 * kd + 103 ] += +1 ;

    /*reaction 105: CH3 + C2H5 => CH4 + C2H4 */
    nuki[ 11 * kd + 104 ] += -1 ;
    nuki[ 14 * kd + 104 ] += -1 ;
    nuki[ 15 * kd + 104 ] += +1 ;
    nuki[ 16 * kd + 104 ] += +1 ;

    /*reaction 106: C2H5 + O2 => C2H4 + HO2 */
    nuki[ 14 * kd + 105 ] += -1 ;
    nuki[ 5 * kd + 105 ] += -1 ;
    nuki[ 16 * kd + 105 ] += +1 ;
    nuki[ 8 * kd + 105 ] += +1 ;

    /*reaction 107: C2H4 + HO2 => C2H5 + O2 */
    nuki[ 16 * kd + 106 ] += -1 ;
    nuki[ 8 * kd + 106 ] += -1 ;
    nuki[ 14 * kd + 106 ] += +1 ;
    nuki[ 5 * kd + 106 ] += +1 ;

    /*reaction 108: C2H5 + HO2 => C2H6 + O2 */
    nuki[ 14 * kd + 107 ] += -1 ;
    nuki[ 8 * kd + 107 ] += -1 ;
    nuki[ 17 * kd + 107 ] += +1 ;
    nuki[ 5 * kd + 107 ] += +1 ;

    /*reaction 109: H + C2H5 => C2H6 */
    nuki[ 2 * kd + 108 ] += -1 ;
    nuki[ 14 * kd + 108 ] += -1 ;
    nuki[ 17 * kd + 108 ] += +1 ;

    /*reaction 110: C2H6 + H => C2H5 + H2 */
    nuki[ 17 * kd + 109 ] += -1 ;
    nuki[ 2 * kd + 109 ] += -1 ;
    nuki[ 14 * kd + 109 ] += +1 ;
    nuki[ 4 * kd + 109 ] += +1 ;

    /*reaction 111: C2H5 + H2 => C2H6 + H */
    nuki[ 14 * kd + 110 ] += -1 ;
    nuki[ 4 * kd + 110 ] += -1 ;
    nuki[ 17 * kd + 110 ] += +1 ;
    nuki[ 2 * kd + 110 ] += +1 ;

    /*reaction 112: C2H6 + OH => C2H5 + H2O */
    nuki[ 17 * kd + 111 ] += -1 ;
    nuki[ 3 * kd + 111 ] += -1 ;
    nuki[ 14 * kd + 111 ] += +1 ;
    nuki[ 6 * kd + 111 ] += +1 ;

    /*reaction 113: CH2GSG + C2H6 => CH3 + C2H5 */
    nuki[ 9 * kd + 112 ] += -1 ;
    nuki[ 17 * kd + 112 ] += -1 ;
    nuki[ 11 * kd + 112 ] += +1 ;
    nuki[ 14 * kd + 112 ] += +1 ;

    /*reaction 114: C2H6 + O => C2H5 + OH */
    nuki[ 17 * kd + 113 ] += -1 ;
    nuki[ 1 * kd + 113 ] += -1 ;
    nuki[ 14 * kd + 113 ] += +1 ;
    nuki[ 3 * kd + 113 ] += +1 ;

    /*reaction 115: C2H6 + CH3 => C2H5 + CH4 */
    nuki[ 17 * kd + 114 ] += -1 ;
    nuki[ 11 * kd + 114 ] += -1 ;
    nuki[ 14 * kd + 114 ] += +1 ;
    nuki[ 15 * kd + 114 ] += +1 ;

    /*reaction 116: HCCO + O => H + 2.000000 CO */
    nuki[ 23 * kd + 115 ] += -1 ;
    nuki[ 1 * kd + 115 ] += -1 ;
    nuki[ 2 * kd + 115 ] += +1 ;
    nuki[ 12 * kd + 115 ] += +2 ;

    /*reaction 117: HCCO + O2 => CO2 + HCO */
    nuki[ 23 * kd + 116 ] += -1 ;
    nuki[ 5 * kd + 116 ] += -1 ;
    nuki[ 18 * kd + 116 ] += +1 ;
    nuki[ 19 * kd + 116 ] += +1 ;

    /*reaction 118: HCCO + OH => 2.000000 HCO */
    nuki[ 23 * kd + 117 ] += -1 ;
    nuki[ 3 * kd + 117 ] += -1 ;
    nuki[ 19 * kd + 117 ] += +2 ;

    /*reaction 119: HCCO + H => CH2GSG + CO */
    nuki[ 23 * kd + 118 ] += -1 ;
    nuki[ 2 * kd + 118 ] += -1 ;
    nuki[ 9 * kd + 118 ] += +1 ;
    nuki[ 12 * kd + 118 ] += +1 ;

    /*reaction 120: CH2GSG + CO => HCCO + H */
    nuki[ 9 * kd + 119 ] += -1 ;
    nuki[ 12 * kd + 119 ] += -1 ;
    nuki[ 23 * kd + 119 ] += +1 ;
    nuki[ 2 * kd + 119 ] += +1 ;

    /*reaction 121: CH2CHO + O2 => CH2O + CO + OH */
    nuki[ 25 * kd + 120 ] += -1 ;
    nuki[ 5 * kd + 120 ] += -1 ;
    nuki[ 10 * kd + 120 ] += +1 ;
    nuki[ 12 * kd + 120 ] += +1 ;
    nuki[ 3 * kd + 120 ] += +1 ;

    /*reaction 122: C3H5XA + CH2O => C3H6 + HCO */
    nuki[ 26 * kd + 121 ] += -1 ;
    nuki[ 10 * kd + 121 ] += -1 ;
    nuki[ 27 * kd + 121 ] += +1 ;
    nuki[ 19 * kd + 121 ] += +1 ;

    /*reaction 123: C3H6 + HCO => C3H5XA + CH2O */
    nuki[ 27 * kd + 122 ] += -1 ;
    nuki[ 19 * kd + 122 ] += -1 ;
    nuki[ 26 * kd + 122 ] += +1 ;
    nuki[ 10 * kd + 122 ] += +1 ;

    /*reaction 124: C3H5XA + HO2 => C3H5O + OH */
    nuki[ 26 * kd + 123 ] += -1 ;
    nuki[ 8 * kd + 123 ] += -1 ;
    nuki[ 28 * kd + 123 ] += +1 ;
    nuki[ 3 * kd + 123 ] += +1 ;

    /*reaction 125: C3H5XA + CH3O2 => C3H5O + CH3O */
    nuki[ 26 * kd + 124 ] += -1 ;
    nuki[ 20 * kd + 124 ] += -1 ;
    nuki[ 28 * kd + 124 ] += +1 ;
    nuki[ 13 * kd + 124 ] += +1 ;

    /*reaction 126: C3H5XA + HO2 => C3H6 + O2 */
    nuki[ 26 * kd + 125 ] += -1 ;
    nuki[ 8 * kd + 125 ] += -1 ;
    nuki[ 27 * kd + 125 ] += +1 ;
    nuki[ 5 * kd + 125 ] += +1 ;

    /*reaction 127: C3H5XA + O2 => CH2CHO + CH2O */
    nuki[ 26 * kd + 126 ] += -1 ;
    nuki[ 5 * kd + 126 ] += -1 ;
    nuki[ 25 * kd + 126 ] += +1 ;
    nuki[ 10 * kd + 126 ] += +1 ;

    /*reaction 128: C3H5XA + O2 => C2H2 + CH2O + OH */
    nuki[ 26 * kd + 127 ] += -1 ;
    nuki[ 5 * kd + 127 ] += -1 ;
    nuki[ 22 * kd + 127 ] += +1 ;
    nuki[ 10 * kd + 127 ] += +1 ;
    nuki[ 3 * kd + 127 ] += +1 ;

    /*reaction 129: C3H5XA + H => C3H6 */
    nuki[ 26 * kd + 128 ] += -1 ;
    nuki[ 2 * kd + 128 ] += -1 ;
    nuki[ 27 * kd + 128 ] += +1 ;

    /*reaction 130: C3H5XA => C2H2 + CH3 */
    nuki[ 26 * kd + 129 ] += -1 ;
    nuki[ 22 * kd + 129 ] += +1 ;
    nuki[ 11 * kd + 129 ] += +1 ;

    /*reaction 131: C3H6 => C2H3 + CH3 */
    nuki[ 27 * kd + 130 ] += -1 ;
    nuki[ 24 * kd + 130 ] += +1 ;
    nuki[ 11 * kd + 130 ] += +1 ;

    /*reaction 132: C2H3 + CH3 => C3H6 */
    nuki[ 24 * kd + 131 ] += -1 ;
    nuki[ 11 * kd + 131 ] += -1 ;
    nuki[ 27 * kd + 131 ] += +1 ;

    /*reaction 133: C3H6 + O => C2H5 + HCO */
    nuki[ 27 * kd + 132 ] += -1 ;
    nuki[ 1 * kd + 132 ] += -1 ;
    nuki[ 14 * kd + 132 ] += +1 ;
    nuki[ 19 * kd + 132 ] += +1 ;

    /*reaction 134: C3H6 + H => C2H4 + CH3 */
    nuki[ 27 * kd + 133 ] += -1 ;
    nuki[ 2 * kd + 133 ] += -1 ;
    nuki[ 16 * kd + 133 ] += +1 ;
    nuki[ 11 * kd + 133 ] += +1 ;

    /*reaction 135: C2H4 + CH3 => C3H6 + H */
    nuki[ 16 * kd + 134 ] += -1 ;
    nuki[ 11 * kd + 134 ] += -1 ;
    nuki[ 27 * kd + 134 ] += +1 ;
    nuki[ 2 * kd + 134 ] += +1 ;

    /*reaction 136: C3H6 + H => C3H5XA + H2 */
    nuki[ 27 * kd + 135 ] += -1 ;
    nuki[ 2 * kd + 135 ] += -1 ;
    nuki[ 26 * kd + 135 ] += +1 ;
    nuki[ 4 * kd + 135 ] += +1 ;

    /*reaction 137: C3H5XA + H2 => C3H6 + H */
    nuki[ 26 * kd + 136 ] += -1 ;
    nuki[ 4 * kd + 136 ] += -1 ;
    nuki[ 27 * kd + 136 ] += +1 ;
    nuki[ 2 * kd + 136 ] += +1 ;

    /*reaction 138: C3H6 + OH => C3H5XA + H2O */
    nuki[ 27 * kd + 137 ] += -1 ;
    nuki[ 3 * kd + 137 ] += -1 ;
    nuki[ 26 * kd + 137 ] += +1 ;
    nuki[ 6 * kd + 137 ] += +1 ;

    /*reaction 139: C3H6 + O => C3H5XA + OH */
    nuki[ 27 * kd + 138 ] += -1 ;
    nuki[ 1 * kd + 138 ] += -1 ;
    nuki[ 26 * kd + 138 ] += +1 ;
    nuki[ 3 * kd + 138 ] += +1 ;

    /*reaction 140: C3H6 + CH3 => C3H5XA + CH4 */
    nuki[ 27 * kd + 139 ] += -1 ;
    nuki[ 11 * kd + 139 ] += -1 ;
    nuki[ 26 * kd + 139 ] += +1 ;
    nuki[ 15 * kd + 139 ] += +1 ;

    /*reaction 141: IXC3H7 + O2 => C3H6 + HO2 */
    nuki[ 29 * kd + 140 ] += -1 ;
    nuki[ 5 * kd + 140 ] += -1 ;
    nuki[ 27 * kd + 140 ] += +1 ;
    nuki[ 8 * kd + 140 ] += +1 ;

    /*reaction 142: IXC3H7 => H + C3H6 */
    nuki[ 29 * kd + 141 ] += -1 ;
    nuki[ 2 * kd + 141 ] += +1 ;
    nuki[ 27 * kd + 141 ] += +1 ;

    /*reaction 143: H + C3H6 => IXC3H7 */
    nuki[ 2 * kd + 142 ] += -1 ;
    nuki[ 27 * kd + 142 ] += -1 ;
    nuki[ 29 * kd + 142 ] += +1 ;

    /*reaction 144: NXC3H7 => CH3 + C2H4 */
    nuki[ 30 * kd + 143 ] += -1 ;
    nuki[ 11 * kd + 143 ] += +1 ;
    nuki[ 16 * kd + 143 ] += +1 ;

    /*reaction 145: CH3 + C2H4 => NXC3H7 */
    nuki[ 11 * kd + 144 ] += -1 ;
    nuki[ 16 * kd + 144 ] += -1 ;
    nuki[ 30 * kd + 144 ] += +1 ;

    /*reaction 146: NXC3H7 + HO2 => C3H8 + O2 */
    nuki[ 30 * kd + 145 ] += -1 ;
    nuki[ 8 * kd + 145 ] += -1 ;
    nuki[ 31 * kd + 145 ] += +1 ;
    nuki[ 5 * kd + 145 ] += +1 ;

    /*reaction 147: NXC3H7 + O2 => C3H6 + HO2 */
    nuki[ 30 * kd + 146 ] += -1 ;
    nuki[ 5 * kd + 146 ] += -1 ;
    nuki[ 27 * kd + 146 ] += +1 ;
    nuki[ 8 * kd + 146 ] += +1 ;

    /*reaction 148: NXC3H7 => H + C3H6 */
    nuki[ 30 * kd + 147 ] += -1 ;
    nuki[ 2 * kd + 147 ] += +1 ;
    nuki[ 27 * kd + 147 ] += +1 ;

    /*reaction 149: H + C3H6 => NXC3H7 */
    nuki[ 2 * kd + 148 ] += -1 ;
    nuki[ 27 * kd + 148 ] += -1 ;
    nuki[ 30 * kd + 148 ] += +1 ;

    /*reaction 150: C3H8 + OH => NXC3H7 + H2O */
    nuki[ 31 * kd + 149 ] += -1 ;
    nuki[ 3 * kd + 149 ] += -1 ;
    nuki[ 30 * kd + 149 ] += +1 ;
    nuki[ 6 * kd + 149 ] += +1 ;

    /*reaction 151: C3H8 + HO2 => NXC3H7 + H2O2 */
    nuki[ 31 * kd + 150 ] += -1 ;
    nuki[ 8 * kd + 150 ] += -1 ;
    nuki[ 30 * kd + 150 ] += +1 ;
    nuki[ 7 * kd + 150 ] += +1 ;

    /*reaction 152: H + C3H8 <=> H2 + NXC3H7 */
    nuki[ 2 * kd + 151 ] += -1 ;
    nuki[ 31 * kd + 151 ] += -1 ;
    nuki[ 4 * kd + 151 ] += +1 ;
    nuki[ 30 * kd + 151 ] += +1 ;

    /*reaction 153: C3H8 + OH => IXC3H7 + H2O */
    nuki[ 31 * kd + 152 ] += -1 ;
    nuki[ 3 * kd + 152 ] += -1 ;
    nuki[ 29 * kd + 152 ] += +1 ;
    nuki[ 6 * kd + 152 ] += +1 ;

    /*reaction 154: CH3 + C3H8 => CH4 + IXC3H7 */
    nuki[ 11 * kd + 153 ] += -1 ;
    nuki[ 31 * kd + 153 ] += -1 ;
    nuki[ 15 * kd + 153 ] += +1 ;
    nuki[ 29 * kd + 153 ] += +1 ;

    /*reaction 155: CH3 + C3H8 => CH4 + NXC3H7 */
    nuki[ 11 * kd + 154 ] += -1 ;
    nuki[ 31 * kd + 154 ] += -1 ;
    nuki[ 15 * kd + 154 ] += +1 ;
    nuki[ 30 * kd + 154 ] += +1 ;

    /*reaction 156: C3H8 + O => IXC3H7 + OH */
    nuki[ 31 * kd + 155 ] += -1 ;
    nuki[ 1 * kd + 155 ] += -1 ;
    nuki[ 29 * kd + 155 ] += +1 ;
    nuki[ 3 * kd + 155 ] += +1 ;

    /*reaction 157: C3H8 + HO2 => IXC3H7 + H2O2 */
    nuki[ 31 * kd + 156 ] += -1 ;
    nuki[ 8 * kd + 156 ] += -1 ;
    nuki[ 29 * kd + 156 ] += +1 ;
    nuki[ 7 * kd + 156 ] += +1 ;

    /*reaction 158: C3H8 + O => NXC3H7 + OH */
    nuki[ 31 * kd + 157 ] += -1 ;
    nuki[ 1 * kd + 157 ] += -1 ;
    nuki[ 30 * kd + 157 ] += +1 ;
    nuki[ 3 * kd + 157 ] += +1 ;

    /*reaction 159: C3H8 + O2 => IXC3H7 + HO2 */
    nuki[ 31 * kd + 158 ] += -1 ;
    nuki[ 5 * kd + 158 ] += -1 ;
    nuki[ 29 * kd + 158 ] += +1 ;
    nuki[ 8 * kd + 158 ] += +1 ;

    /*reaction 160: IXC3H7 + HO2 => C3H8 + O2 */
    nuki[ 29 * kd + 159 ] += -1 ;
    nuki[ 8 * kd + 159 ] += -1 ;
    nuki[ 31 * kd + 159 ] += +1 ;
    nuki[ 5 * kd + 159 ] += +1 ;

    /*reaction 161: H + C3H8 => H2 + IXC3H7 */
    nuki[ 2 * kd + 160 ] += -1 ;
    nuki[ 31 * kd + 160 ] += -1 ;
    nuki[ 4 * kd + 160 ] += +1 ;
    nuki[ 29 * kd + 160 ] += +1 ;

    /*reaction 162: H2 + IXC3H7 => H + C3H8 */
    nuki[ 4 * kd + 161 ] += -1 ;
    nuki[ 29 * kd + 161 ] += -1 ;
    nuki[ 2 * kd + 161 ] += +1 ;
    nuki[ 31 * kd + 161 ] += +1 ;

    /*reaction 163: C3H5O => C2H3 + CH2O */
    nuki[ 28 * kd + 162 ] += -1 ;
    nuki[ 24 * kd + 162 ] += +1 ;
    nuki[ 10 * kd + 162 ] += +1 ;

    /*reaction 164: IXC3H7O2 => IXC3H7 + O2 */
    nuki[ 32 * kd + 163 ] += -1 ;
    nuki[ 29 * kd + 163 ] += +1 ;
    nuki[ 5 * kd + 163 ] += +1 ;

    /*reaction 165: IXC3H7 + O2 => IXC3H7O2 */
    nuki[ 29 * kd + 164 ] += -1 ;
    nuki[ 5 * kd + 164 ] += -1 ;
    nuki[ 32 * kd + 164 ] += +1 ;

    /*reaction 166: NXC3H7O2 => NXC3H7 + O2 */
    nuki[ 33 * kd + 165 ] += -1 ;
    nuki[ 30 * kd + 165 ] += +1 ;
    nuki[ 5 * kd + 165 ] += +1 ;

    /*reaction 167: NXC3H7 + O2 => NXC3H7O2 */
    nuki[ 30 * kd + 166 ] += -1 ;
    nuki[ 5 * kd + 166 ] += -1 ;
    nuki[ 33 * kd + 166 ] += +1 ;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cv_R(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: N2 */
        species[0] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
        /*species 1: O */
        species[1] =
            +1.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 2: H */
        species[2] =
            +1.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: H2 */
        species[4] =
            +2.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 5: O2 */
        species[5] =
            +2.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 6: H2O */
        species[6] =
            +2.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +2.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 9: CH2GSG */
        species[9] =
            +2.97126500e+00
            -1.69908900e-04 * tc[1]
            +1.02536900e-06 * tc[2]
            +2.49255100e-09 * tc[3]
            -1.98126600e-12 * tc[4];
        /*species 10: CH2O */
        species[10] =
            +6.52731000e-01
            +1.26314400e-02 * tc[1]
            -1.88816800e-05 * tc[2]
            +2.05003100e-08 * tc[3]
            -8.41323700e-12 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +1.43044300e+00
            +1.11241000e-02 * tc[1]
            -1.68022000e-05 * tc[2]
            +1.62182900e-08 * tc[3]
            -5.86495300e-12 * tc[4];
        /*species 12: CO */
        species[12] =
            +2.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +1.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561100e-12 * tc[4];
        /*species 14: C2H5 */
        species[14] =
            +1.69070200e+00
            +8.71913300e-03 * tc[1]
            +4.41983900e-06 * tc[2]
            +9.33870300e-10 * tc[3]
            -3.92777300e-12 * tc[4];
        /*species 15: CH4 */
        species[15] =
            -2.21258500e-01
            +1.74766800e-02 * tc[1]
            -2.78340900e-05 * tc[2]
            +3.04970800e-08 * tc[3]
            -1.22393100e-11 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            -1.86148800e+00
            +2.79616300e-02 * tc[1]
            -3.38867700e-05 * tc[2]
            +2.78515200e-08 * tc[3]
            -9.73787900e-12 * tc[4];
        /*species 18: CO2 */
        species[18] =
            +1.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
        /*species 19: HCO */
        species[19] =
            +1.89833000e+00
            +6.19914700e-03 * tc[1]
            -9.62308400e-06 * tc[2]
            +1.08982500e-08 * tc[3]
            -4.57488500e-12 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +1.01356200e+00
            +1.51904500e-02 * tc[1]
            -1.61631900e-05 * tc[2]
            +9.07899200e-09 * tc[3]
            -1.91274600e-12 * tc[4];
        /*species 23: HCCO */
        species[23] =
            +4.04796500e+00
            +4.45347800e-03 * tc[1]
            +2.26828300e-07 * tc[2]
            -1.48209500e-09 * tc[3]
            +2.25074200e-13 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +1.45927600e+00
            +7.37147600e-03 * tc[1]
            +2.10987300e-06 * tc[2]
            -1.32164200e-09 * tc[3]
            -1.18478400e-12 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +2.40906200e+00
            +1.07385700e-02 * tc[1]
            +1.89149200e-06 * tc[2]
            -7.15858300e-09 * tc[3]
            +2.86738500e-12 * tc[4];
        /*species 29: IXC3H7 */
        species[29] =
            +7.13300000e-01
            +2.54261600e-02 * tc[1]
            +1.58080800e-06 * tc[2]
            -1.82128600e-08 * tc[3]
            +8.82771000e-12 * tc[4];
        /*species 30: NXC3H7 */
        species[30] =
            +9.22537000e-01
            +2.47892700e-02 * tc[1]
            +1.81024900e-06 * tc[2]
            -1.78326600e-08 * tc[3]
            +8.58299600e-12 * tc[4];
        /*species 31: C3H8 */
        species[31] =
            -1.03079200e-01
            +2.66898600e-02 * tc[1]
            +5.43142500e-06 * tc[2]
            -2.12600100e-08 * tc[3]
            +9.24333000e-12 * tc[4];
    } else {
        /*species 0: N2 */
        species[0] =
            +1.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 1: O */
        species[1] =
            +1.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 2: H */
        species[2] =
            +1.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: H2 */
        species[4] =
            +1.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 5: O2 */
        species[5] =
            +2.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 6: H2O */
        species[6] =
            +1.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 9: CH2GSG */
        species[9] =
            +2.55288900e+00
            +2.06678800e-03 * tc[1]
            -1.91411600e-07 * tc[2]
            -1.10467300e-10 * tc[3]
            +2.02135000e-14 * tc[4];
        /*species 10: CH2O */
        species[10] =
            +1.99560600e+00
            +6.68132100e-03 * tc[1]
            -2.62895500e-06 * tc[2]
            +4.73715300e-10 * tc[3]
            -3.21251700e-14 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +1.84405200e+00
            +6.13797400e-03 * tc[1]
            -2.23034500e-06 * tc[2]
            +3.78516100e-10 * tc[3]
            -2.45215900e-14 * tc[4];
        /*species 12: CO */
        species[12] =
            +2.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +2.77080000e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 14: C2H5 */
        species[14] =
            +6.19048000e+00
            +6.48407700e-03 * tc[1]
            -6.42806500e-07 * tc[2]
            -2.34787900e-10 * tc[3]
            +3.88087700e-14 * tc[4];
        /*species 15: CH4 */
        species[15] =
            +6.83479000e-01
            +1.02372400e-02 * tc[1]
            -3.87512900e-06 * tc[2]
            +6.78558500e-10 * tc[3]
            -4.50342300e-14 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +2.52841900e+00
            +1.14851800e-02 * tc[1]
            -4.41838500e-06 * tc[2]
            +7.84460100e-10 * tc[3]
            -5.26684800e-14 * tc[4];
        /*species 18: CO2 */
        species[18] =
            +3.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
        /*species 19: HCO */
        species[19] =
            +2.55727100e+00
            +3.34557300e-03 * tc[1]
            -1.33500600e-06 * tc[2]
            +2.47057300e-10 * tc[3]
            -1.71385100e-14 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +3.43677000e+00
            +5.37603900e-03 * tc[1]
            -1.91281700e-06 * tc[2]
            +3.28637900e-10 * tc[3]
            -2.15671000e-14 * tc[4];
        /*species 23: HCCO */
        species[23] =
            +5.75807300e+00
            +2.00040000e-03 * tc[1]
            -2.02760700e-07 * tc[2]
            -1.04113200e-10 * tc[3]
            +1.96516500e-14 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +4.93346800e+00
            +4.01774600e-03 * tc[1]
            -3.96674000e-07 * tc[2]
            -1.44126700e-10 * tc[3]
            +2.37864400e-14 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +4.97567000e+00
            +8.13059100e-03 * tc[1]
            -2.74362400e-06 * tc[2]
            +4.07030400e-10 * tc[3]
            -2.17601700e-14 * tc[4];
        /*species 29: IXC3H7 */
        species[29] =
            +7.06336900e+00
            +1.57448800e-02 * tc[1]
            -5.18239200e-06 * tc[2]
            +7.47724500e-10 * tc[3]
            -3.85442200e-14 * tc[4];
        /*species 30: NXC3H7 */
        species[30] =
            +6.97829100e+00
            +1.57611300e-02 * tc[1]
            -5.17324300e-06 * tc[2]
            +7.44389200e-10 * tc[3]
            -3.82497800e-14 * tc[4];
        /*species 31: C3H8 */
        species[31] =
            +6.52521700e+00
            +1.88903400e-02 * tc[1]
            -6.28392400e-06 * tc[2]
            +9.17937300e-10 * tc[3]
            -4.81241000e-14 * tc[4];
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 3: OH */
        species[3] =
            +2.43586219e+00
            +2.02235804e-04 * tc[1]
            -1.13546412e-07 * tc[2]
            +2.42445149e-10 * tc[3]
            -7.43651031e-14 * tc[4];
    } else {
        /*species 3: OH */
        species[3] =
            +1.62599754e+00
            +1.31992406e-03 * tc[1]
            -3.59724670e-07 * tc[2]
            +4.25630800e-11 * tc[3]
            -1.82048016e-15 * tc[4];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 8: HO2 */
        species[8] =
            +2.18310656e+00
            +3.66767950e-03 * tc[1]
            -9.32385122e-07 * tc[2]
            -3.25852919e-10 * tc[3]
            +1.51139912e-13 * tc[4];
        /*species 21: CH3O2H */
        species[21] =
            +2.23442817e+00
            +1.90129767e-02 * tc[1]
            -1.13386287e-05 * tc[2]
            +3.40306653e-09 * tc[3]
            -4.11830222e-13 * tc[4];
    } else {
        /*species 8: HO2 */
        species[8] =
            +3.10547423e+00
            +2.38452835e-03 * tc[1]
            -8.06347989e-07 * tc[2]
            +1.24191723e-10 * tc[3]
            -7.16400108e-15 * tc[4];
        /*species 21: CH3O2H */
        species[21] =
            +7.43117091e+00
            +8.06817909e-03 * tc[1]
            -2.77094921e-06 * tc[2]
            +4.31332243e-10 * tc[3]
            -2.50692146e-14 * tc[4];
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 17: C2H6 */
        species[17] =
            -1.02528543e+00
            +2.40764754e-02 * tc[1]
            -1.11893472e-05 * tc[2]
            +2.08340901e-09 * tc[3]
            -5.29868616e-14 * tc[4];
        /*species 33: NXC3H7O2 */
        species[33] =
            +1.10731492e+00
            +3.96164986e-02 * tc[1]
            -2.49491599e-05 * tc[2]
            +8.59450300e-09 * tc[3]
            -1.31240330e-12 * tc[4];
    } else {
        /*species 17: C2H6 */
        species[17] =
            +5.10683385e+00
            +1.29236361e-02 * tc[1]
            -4.42527196e-06 * tc[2]
            +6.87391726e-10 * tc[3]
            -3.98901732e-14 * tc[4];
        /*species 33: NXC3H7O2 */
        species[33] =
            +1.16327059e+01
            +1.69910726e-02 * tc[1]
            -5.88866873e-06 * tc[2]
            +9.22195396e-10 * tc[3]
            -5.38230675e-14 * tc[4];
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 20: CH3O2 */
        species[20] =
            +3.26146906e+00
            +1.00873599e-02 * tc[1]
            -3.21506184e-06 * tc[2]
            +2.09409267e-10 * tc[3]
            +4.18339103e-14 * tc[4];
    } else {
        /*species 20: CH3O2 */
        species[20] =
            +4.95787891e+00
            +7.90728626e-03 * tc[1]
            -2.68246234e-06 * tc[2]
            +4.13891337e-10 * tc[3]
            -2.39007330e-14 * tc[4];
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 26: C3H5XA */
        species[26] =
            -1.52913196e+00
            +3.34559100e-02 * tc[1]
            -2.53401027e-05 * tc[2]
            +1.02865754e-08 * tc[3]
            -1.73258340e-12 * tc[4];
    } else {
        /*species 26: C3H5XA */
        species[26] =
            +7.45883958e+00
            +1.12695483e-02 * tc[1]
            -3.83792864e-06 * tc[2]
            +5.94059119e-10 * tc[3]
            -3.43918030e-14 * tc[4];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 27: C3H6 */
        species[27] =
            -6.05384556e-01
            +2.89107662e-02 * tc[1]
            -1.54886808e-05 * tc[2]
            +3.88814209e-09 * tc[3]
            -3.37890352e-13 * tc[4];
        /*species 32: IXC3H7O2 */
        species[32] =
            +4.99416390e-01
            +4.43081205e-02 * tc[1]
            -3.22414456e-05 * tc[2]
            +1.29687136e-08 * tc[3]
            -2.23370569e-12 * tc[4];
    } else {
        /*species 27: C3H6 */
        species[27] =
            +7.01595958e+00
            +1.37023634e-02 * tc[1]
            -4.66249733e-06 * tc[2]
            +7.21254402e-10 * tc[3]
            -4.17370126e-14 * tc[4];
        /*species 32: IXC3H7O2 */
        species[32] =
            +1.22493493e+01
            +1.64082190e-02 * tc[1]
            -5.67432062e-06 * tc[2]
            +8.87336340e-10 * tc[3]
            -5.17361535e-14 * tc[4];
    }

    /*species with midpoint at T=1380 kelvin */
    if (T < 1380) {
        /*species 28: C3H5O */
        species[28] =
            +1.98225820e-01
            +3.05579837e-02 * tc[1]
            -1.80630276e-05 * tc[2]
            +4.86150033e-09 * tc[3]
            -4.19854562e-13 * tc[4];
    } else {
        /*species 28: C3H5O */
        species[28] =
            +9.25517520e+00
            +1.14983720e-02 * tc[1]
            -3.84645659e-06 * tc[2]
            +5.88910346e-10 * tc[3]
            -3.38557923e-14 * tc[4];
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
        /*species 0: N2 */
        species[0] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
        /*species 1: O */
        species[1] =
            +2.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 2: H */
        species[2] =
            +2.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: H2 */
        species[4] =
            +3.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 5: O2 */
        species[5] =
            +3.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 6: H2O */
        species[6] =
            +3.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 9: CH2GSG */
        species[9] =
            +3.97126500e+00
            -1.69908900e-04 * tc[1]
            +1.02536900e-06 * tc[2]
            +2.49255100e-09 * tc[3]
            -1.98126600e-12 * tc[4];
        /*species 10: CH2O */
        species[10] =
            +1.65273100e+00
            +1.26314400e-02 * tc[1]
            -1.88816800e-05 * tc[2]
            +2.05003100e-08 * tc[3]
            -8.41323700e-12 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +2.43044300e+00
            +1.11241000e-02 * tc[1]
            -1.68022000e-05 * tc[2]
            +1.62182900e-08 * tc[3]
            -5.86495300e-12 * tc[4];
        /*species 12: CO */
        species[12] =
            +3.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +2.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561100e-12 * tc[4];
        /*species 14: C2H5 */
        species[14] =
            +2.69070200e+00
            +8.71913300e-03 * tc[1]
            +4.41983900e-06 * tc[2]
            +9.33870300e-10 * tc[3]
            -3.92777300e-12 * tc[4];
        /*species 15: CH4 */
        species[15] =
            +7.78741500e-01
            +1.74766800e-02 * tc[1]
            -2.78340900e-05 * tc[2]
            +3.04970800e-08 * tc[3]
            -1.22393100e-11 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            -8.61488000e-01
            +2.79616300e-02 * tc[1]
            -3.38867700e-05 * tc[2]
            +2.78515200e-08 * tc[3]
            -9.73787900e-12 * tc[4];
        /*species 18: CO2 */
        species[18] =
            +2.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
        /*species 19: HCO */
        species[19] =
            +2.89833000e+00
            +6.19914700e-03 * tc[1]
            -9.62308400e-06 * tc[2]
            +1.08982500e-08 * tc[3]
            -4.57488500e-12 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.01356200e+00
            +1.51904500e-02 * tc[1]
            -1.61631900e-05 * tc[2]
            +9.07899200e-09 * tc[3]
            -1.91274600e-12 * tc[4];
        /*species 23: HCCO */
        species[23] =
            +5.04796500e+00
            +4.45347800e-03 * tc[1]
            +2.26828300e-07 * tc[2]
            -1.48209500e-09 * tc[3]
            +2.25074200e-13 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +2.45927600e+00
            +7.37147600e-03 * tc[1]
            +2.10987300e-06 * tc[2]
            -1.32164200e-09 * tc[3]
            -1.18478400e-12 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +3.40906200e+00
            +1.07385700e-02 * tc[1]
            +1.89149200e-06 * tc[2]
            -7.15858300e-09 * tc[3]
            +2.86738500e-12 * tc[4];
        /*species 29: IXC3H7 */
        species[29] =
            +1.71330000e+00
            +2.54261600e-02 * tc[1]
            +1.58080800e-06 * tc[2]
            -1.82128600e-08 * tc[3]
            +8.82771000e-12 * tc[4];
        /*species 30: NXC3H7 */
        species[30] =
            +1.92253700e+00
            +2.47892700e-02 * tc[1]
            +1.81024900e-06 * tc[2]
            -1.78326600e-08 * tc[3]
            +8.58299600e-12 * tc[4];
        /*species 31: C3H8 */
        species[31] =
            +8.96920800e-01
            +2.66898600e-02 * tc[1]
            +5.43142500e-06 * tc[2]
            -2.12600100e-08 * tc[3]
            +9.24333000e-12 * tc[4];
    } else {
        /*species 0: N2 */
        species[0] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 1: O */
        species[1] =
            +2.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 2: H */
        species[2] =
            +2.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: H2 */
        species[4] =
            +2.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 5: O2 */
        species[5] =
            +3.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 6: H2O */
        species[6] =
            +2.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 9: CH2GSG */
        species[9] =
            +3.55288900e+00
            +2.06678800e-03 * tc[1]
            -1.91411600e-07 * tc[2]
            -1.10467300e-10 * tc[3]
            +2.02135000e-14 * tc[4];
        /*species 10: CH2O */
        species[10] =
            +2.99560600e+00
            +6.68132100e-03 * tc[1]
            -2.62895500e-06 * tc[2]
            +4.73715300e-10 * tc[3]
            -3.21251700e-14 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +2.84405200e+00
            +6.13797400e-03 * tc[1]
            -2.23034500e-06 * tc[2]
            +3.78516100e-10 * tc[3]
            -2.45215900e-14 * tc[4];
        /*species 12: CO */
        species[12] =
            +3.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +3.77080000e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 14: C2H5 */
        species[14] =
            +7.19048000e+00
            +6.48407700e-03 * tc[1]
            -6.42806500e-07 * tc[2]
            -2.34787900e-10 * tc[3]
            +3.88087700e-14 * tc[4];
        /*species 15: CH4 */
        species[15] =
            +1.68347900e+00
            +1.02372400e-02 * tc[1]
            -3.87512900e-06 * tc[2]
            +6.78558500e-10 * tc[3]
            -4.50342300e-14 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +3.52841900e+00
            +1.14851800e-02 * tc[1]
            -4.41838500e-06 * tc[2]
            +7.84460100e-10 * tc[3]
            -5.26684800e-14 * tc[4];
        /*species 18: CO2 */
        species[18] =
            +4.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
        /*species 19: HCO */
        species[19] =
            +3.55727100e+00
            +3.34557300e-03 * tc[1]
            -1.33500600e-06 * tc[2]
            +2.47057300e-10 * tc[3]
            -1.71385100e-14 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +4.43677000e+00
            +5.37603900e-03 * tc[1]
            -1.91281700e-06 * tc[2]
            +3.28637900e-10 * tc[3]
            -2.15671000e-14 * tc[4];
        /*species 23: HCCO */
        species[23] =
            +6.75807300e+00
            +2.00040000e-03 * tc[1]
            -2.02760700e-07 * tc[2]
            -1.04113200e-10 * tc[3]
            +1.96516500e-14 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +5.93346800e+00
            +4.01774600e-03 * tc[1]
            -3.96674000e-07 * tc[2]
            -1.44126700e-10 * tc[3]
            +2.37864400e-14 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +5.97567000e+00
            +8.13059100e-03 * tc[1]
            -2.74362400e-06 * tc[2]
            +4.07030400e-10 * tc[3]
            -2.17601700e-14 * tc[4];
        /*species 29: IXC3H7 */
        species[29] =
            +8.06336900e+00
            +1.57448800e-02 * tc[1]
            -5.18239200e-06 * tc[2]
            +7.47724500e-10 * tc[3]
            -3.85442200e-14 * tc[4];
        /*species 30: NXC3H7 */
        species[30] =
            +7.97829100e+00
            +1.57611300e-02 * tc[1]
            -5.17324300e-06 * tc[2]
            +7.44389200e-10 * tc[3]
            -3.82497800e-14 * tc[4];
        /*species 31: C3H8 */
        species[31] =
            +7.52521700e+00
            +1.88903400e-02 * tc[1]
            -6.28392400e-06 * tc[2]
            +9.17937300e-10 * tc[3]
            -4.81241000e-14 * tc[4];
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 3: OH */
        species[3] =
            +3.43586219e+00
            +2.02235804e-04 * tc[1]
            -1.13546412e-07 * tc[2]
            +2.42445149e-10 * tc[3]
            -7.43651031e-14 * tc[4];
    } else {
        /*species 3: OH */
        species[3] =
            +2.62599754e+00
            +1.31992406e-03 * tc[1]
            -3.59724670e-07 * tc[2]
            +4.25630800e-11 * tc[3]
            -1.82048016e-15 * tc[4];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 8: HO2 */
        species[8] =
            +3.18310656e+00
            +3.66767950e-03 * tc[1]
            -9.32385122e-07 * tc[2]
            -3.25852919e-10 * tc[3]
            +1.51139912e-13 * tc[4];
        /*species 21: CH3O2H */
        species[21] =
            +3.23442817e+00
            +1.90129767e-02 * tc[1]
            -1.13386287e-05 * tc[2]
            +3.40306653e-09 * tc[3]
            -4.11830222e-13 * tc[4];
    } else {
        /*species 8: HO2 */
        species[8] =
            +4.10547423e+00
            +2.38452835e-03 * tc[1]
            -8.06347989e-07 * tc[2]
            +1.24191723e-10 * tc[3]
            -7.16400108e-15 * tc[4];
        /*species 21: CH3O2H */
        species[21] =
            +8.43117091e+00
            +8.06817909e-03 * tc[1]
            -2.77094921e-06 * tc[2]
            +4.31332243e-10 * tc[3]
            -2.50692146e-14 * tc[4];
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 17: C2H6 */
        species[17] =
            -2.52854344e-02
            +2.40764754e-02 * tc[1]
            -1.11893472e-05 * tc[2]
            +2.08340901e-09 * tc[3]
            -5.29868616e-14 * tc[4];
        /*species 33: NXC3H7O2 */
        species[33] =
            +2.10731492e+00
            +3.96164986e-02 * tc[1]
            -2.49491599e-05 * tc[2]
            +8.59450300e-09 * tc[3]
            -1.31240330e-12 * tc[4];
    } else {
        /*species 17: C2H6 */
        species[17] =
            +6.10683385e+00
            +1.29236361e-02 * tc[1]
            -4.42527196e-06 * tc[2]
            +6.87391726e-10 * tc[3]
            -3.98901732e-14 * tc[4];
        /*species 33: NXC3H7O2 */
        species[33] =
            +1.26327059e+01
            +1.69910726e-02 * tc[1]
            -5.88866873e-06 * tc[2]
            +9.22195396e-10 * tc[3]
            -5.38230675e-14 * tc[4];
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 20: CH3O2 */
        species[20] =
            +4.26146906e+00
            +1.00873599e-02 * tc[1]
            -3.21506184e-06 * tc[2]
            +2.09409267e-10 * tc[3]
            +4.18339103e-14 * tc[4];
    } else {
        /*species 20: CH3O2 */
        species[20] =
            +5.95787891e+00
            +7.90728626e-03 * tc[1]
            -2.68246234e-06 * tc[2]
            +4.13891337e-10 * tc[3]
            -2.39007330e-14 * tc[4];
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 26: C3H5XA */
        species[26] =
            -5.29131958e-01
            +3.34559100e-02 * tc[1]
            -2.53401027e-05 * tc[2]
            +1.02865754e-08 * tc[3]
            -1.73258340e-12 * tc[4];
    } else {
        /*species 26: C3H5XA */
        species[26] =
            +8.45883958e+00
            +1.12695483e-02 * tc[1]
            -3.83792864e-06 * tc[2]
            +5.94059119e-10 * tc[3]
            -3.43918030e-14 * tc[4];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 27: C3H6 */
        species[27] =
            +3.94615444e-01
            +2.89107662e-02 * tc[1]
            -1.54886808e-05 * tc[2]
            +3.88814209e-09 * tc[3]
            -3.37890352e-13 * tc[4];
        /*species 32: IXC3H7O2 */
        species[32] =
            +1.49941639e+00
            +4.43081205e-02 * tc[1]
            -3.22414456e-05 * tc[2]
            +1.29687136e-08 * tc[3]
            -2.23370569e-12 * tc[4];
    } else {
        /*species 27: C3H6 */
        species[27] =
            +8.01595958e+00
            +1.37023634e-02 * tc[1]
            -4.66249733e-06 * tc[2]
            +7.21254402e-10 * tc[3]
            -4.17370126e-14 * tc[4];
        /*species 32: IXC3H7O2 */
        species[32] =
            +1.32493493e+01
            +1.64082190e-02 * tc[1]
            -5.67432062e-06 * tc[2]
            +8.87336340e-10 * tc[3]
            -5.17361535e-14 * tc[4];
    }

    /*species with midpoint at T=1380 kelvin */
    if (T < 1380) {
        /*species 28: C3H5O */
        species[28] =
            +1.19822582e+00
            +3.05579837e-02 * tc[1]
            -1.80630276e-05 * tc[2]
            +4.86150033e-09 * tc[3]
            -4.19854562e-13 * tc[4];
    } else {
        /*species 28: C3H5O */
        species[28] =
            +1.02551752e+01
            +1.14983720e-02 * tc[1]
            -3.84645659e-06 * tc[2]
            +5.88910346e-10 * tc[3]
            -3.38557923e-14 * tc[4];
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
        /*species 0: N2 */
        species[0] =
            -1.020900000000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
        /*species 1: O */
        species[1] =
            +2.914764000000000e+04 * invT
            -1.756599999999997e-02
            -2.946429000000000e+00 * tc[0]
            +8.190830000000000e-04 * tc[1]
            -4.035053333333333e-07 * tc[2]
            +1.335702500000000e-10 * tc[3]
            -1.945348000000000e-14 * tc[4];
        /*species 2: H */
        species[2] =
            +2.547474660000000e+04 * invT
            +2.966385537000000e+00
            -2.501044220000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 4: H2 */
        species[4] =
            -1.012521000000000e+03 * invT
            +6.592218000000000e+00
            -3.298124000000000e+00 * tc[0]
            -4.124721000000000e-04 * tc[1]
            +1.357169166666667e-07 * tc[2]
            +7.896194999999999e-12 * tc[3]
            -2.067436000000000e-14 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.005249000000000e+03 * invT
            -2.821802000000000e+00
            -3.212936000000000e+00 * tc[0]
            -5.637430000000000e-04 * tc[1]
            +9.593583333333333e-08 * tc[2]
            -1.094897500000000e-10 * tc[3]
            +4.384277000000000e-14 * tc[4];
        /*species 6: H2O */
        species[6] =
            -3.020811000000000e+04 * invT
            +7.966090000000001e-01
            -3.386842000000000e+00 * tc[0]
            -1.737491000000000e-03 * tc[1]
            +1.059116000000000e-06 * tc[2]
            -5.807150833333333e-10 * tc[3]
            +1.253294000000000e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.766315000000000e+04 * invT
            -3.396609000000000e+00
            -3.388754000000000e+00 * tc[0]
            -3.284613000000000e-03 * tc[1]
            +2.475021666666666e-08 * tc[2]
            +3.854838333333333e-10 * tc[3]
            -1.235757500000000e-13 * tc[4];
        /*species 9: CH2GSG */
        species[9] =
            +4.989368000000000e+04 * invT
            +3.913732930000000e+00
            -3.971265000000000e+00 * tc[0]
            +8.495445000000000e-05 * tc[1]
            -1.708948333333333e-07 * tc[2]
            -2.077125833333333e-10 * tc[3]
            +9.906330000000000e-14 * tc[4];
        /*species 10: CH2O */
        species[10] =
            -1.486540000000000e+04 * invT
            -1.213208900000000e+01
            -1.652731000000000e+00 * tc[0]
            -6.315720000000000e-03 * tc[1]
            +3.146946666666667e-06 * tc[2]
            -1.708359166666667e-09 * tc[3]
            +4.206618500000000e-13 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +1.642378000000000e+04 * invT
            -4.359351000000000e+00
            -2.430443000000000e+00 * tc[0]
            -5.562050000000000e-03 * tc[1]
            +2.800366666666666e-06 * tc[2]
            -1.351524166666667e-09 * tc[3]
            +2.932476500000000e-13 * tc[4];
        /*species 12: CO */
        species[12] =
            -1.431054000000000e+04 * invT
            -1.586445000000000e+00
            -3.262452000000000e+00 * tc[0]
            -7.559705000000000e-04 * tc[1]
            +6.469591666666667e-07 * tc[2]
            -4.651620000000000e-10 * tc[3]
            +1.237475500000000e-13 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +9.786011000000000e+02 * invT
            -1.104597600000000e+01
            -2.106204000000000e+00 * tc[0]
            -3.608297500000000e-03 * tc[1]
            -8.897453333333333e-07 * tc[2]
            +6.148030000000000e-10 * tc[3]
            -1.037805500000000e-13 * tc[4];
        /*species 14: C2H5 */
        species[14] =
            +1.287040000000000e+04 * invT
            -9.447498000000000e+00
            -2.690702000000000e+00 * tc[0]
            -4.359566500000000e-03 * tc[1]
            -7.366398333333332e-07 * tc[2]
            -7.782252500000001e-11 * tc[3]
            +1.963886500000000e-13 * tc[4];
        /*species 15: CH4 */
        species[15] =
            -9.825228999999999e+03 * invT
            -1.294344850000000e+01
            -7.787415000000000e-01 * tc[0]
            -8.738340000000001e-03 * tc[1]
            +4.639015000000000e-06 * tc[2]
            -2.541423333333333e-09 * tc[3]
            +6.119655000000000e-13 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +5.573046000000000e+03 * invT
            -2.507297800000000e+01
            +8.614880000000000e-01 * tc[0]
            -1.398081500000000e-02 * tc[1]
            +5.647795000000000e-06 * tc[2]
            -2.320960000000000e-09 * tc[3]
            +4.868939500000000e-13 * tc[4];
        /*species 18: CO2 */
        species[18] =
            -4.837314000000000e+04 * invT
            -7.912765000000000e+00
            -2.275725000000000e+00 * tc[0]
            -4.961036000000000e-03 * tc[1]
            +1.734851666666667e-06 * tc[2]
            -5.722239166666667e-10 * tc[3]
            +1.058640000000000e-13 * tc[4];
        /*species 19: HCO */
        species[19] =
            +4.159922000000000e+03 * invT
            -6.085284000000000e+00
            -2.898330000000000e+00 * tc[0]
            -3.099573500000000e-03 * tc[1]
            +1.603847333333333e-06 * tc[2]
            -9.081875000000000e-10 * tc[3]
            +2.287442500000000e-13 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.612444000000000e+04 * invT
            -6.791815999999999e+00
            -2.013562000000000e+00 * tc[0]
            -7.595225000000000e-03 * tc[1]
            +2.693865000000000e-06 * tc[2]
            -7.565826666666667e-10 * tc[3]
            +9.563730000000000e-14 * tc[4];
        /*species 23: HCCO */
        species[23] =
            +1.965892000000000e+04 * invT
            +4.566121099999999e+00
            -5.047965000000000e+00 * tc[0]
            -2.226739000000000e-03 * tc[1]
            -3.780471666666667e-08 * tc[2]
            +1.235079166666667e-10 * tc[3]
            -1.125371000000000e-14 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +3.335225000000000e+04 * invT
            -9.096924000000001e+00
            -2.459276000000000e+00 * tc[0]
            -3.685738000000000e-03 * tc[1]
            -3.516455000000000e-07 * tc[2]
            +1.101368333333333e-10 * tc[3]
            +5.923920000000000e-14 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +1.521477000000000e+03 * invT
            -6.149227999999999e+00
            -3.409062000000000e+00 * tc[0]
            -5.369285000000000e-03 * tc[1]
            -3.152486666666667e-07 * tc[2]
            +5.965485833333333e-10 * tc[3]
            -1.433692500000000e-13 * tc[4];
        /*species 29: IXC3H7 */
        species[29] =
            +7.535809000000000e+03 * invT
            -1.126571000000000e+01
            -1.713300000000000e+00 * tc[0]
            -1.271308000000000e-02 * tc[1]
            -2.634680000000000e-07 * tc[2]
            +1.517738333333333e-09 * tc[3]
            -4.413855000000000e-13 * tc[4];
        /*species 30: NXC3H7 */
        species[30] =
            +9.713281000000001e+03 * invT
            -1.207017300000000e+01
            -1.922537000000000e+00 * tc[0]
            -1.239463500000000e-02 * tc[1]
            -3.017081666666666e-07 * tc[2]
            +1.486055000000000e-09 * tc[3]
            -4.291498000000000e-13 * tc[4];
        /*species 31: C3H8 */
        species[31] =
            -1.395492000000000e+04 * invT
            -1.845840920000000e+01
            -8.969208000000000e-01 * tc[0]
            -1.334493000000000e-02 * tc[1]
            -9.052375000000000e-07 * tc[2]
            +1.771667500000000e-09 * tc[3]
            -4.621665000000000e-13 * tc[4];
    } else {
        /*species 0: N2 */
        species[0] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 1: O */
        species[1] =
            +2.923080000000000e+04 * invT
            -2.378248000000000e+00
            -2.542060000000000e+00 * tc[0]
            +1.377531000000000e-05 * tc[1]
            +5.171338333333333e-10 * tc[2]
            -3.792555833333334e-13 * tc[3]
            +2.184026000000000e-17 * tc[4];
        /*species 2: H */
        species[2] =
            +2.547474660000000e+04 * invT
            +2.966385537000000e+00
            -2.501044220000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 4: H2 */
        species[4] =
            -8.350340000000000e+02 * invT
            +4.346533000000000e+00
            -2.991423000000000e+00 * tc[0]
            -3.500322000000000e-04 * tc[1]
            +9.389715000000000e-09 * tc[2]
            +7.692981666666667e-13 * tc[3]
            -7.913760000000000e-17 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.233930000000000e+03 * invT
            +5.084119999999999e-01
            -3.697578000000000e+00 * tc[0]
            -3.067598500000000e-04 * tc[1]
            +2.098070000000000e-08 * tc[2]
            -1.479400833333333e-12 * tc[3]
            +5.682175000000001e-17 * tc[4];
        /*species 6: H2O */
        species[6] =
            -2.989921000000000e+04 * invT
            -4.190671000000000e+00
            -2.672146000000000e+00 * tc[0]
            -1.528146500000000e-03 * tc[1]
            +1.455043333333333e-07 * tc[2]
            -1.000830000000000e-11 * tc[3]
            +3.195809000000000e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.800696000000000e+04 * invT
            +4.072030000000000e+00
            -4.573167000000000e+00 * tc[0]
            -2.168068000000000e-03 * tc[1]
            +2.457815000000000e-07 * tc[2]
            -1.957420000000000e-11 * tc[3]
            +7.158270000000000e-16 * tc[4];
        /*species 9: CH2GSG */
        species[9] =
            +4.984975000000000e+04 * invT
            +1.866319000000000e+00
            -3.552889000000000e+00 * tc[0]
            -1.033394000000000e-03 * tc[1]
            +3.190193333333333e-08 * tc[2]
            +9.205608333333333e-12 * tc[3]
            -1.010675000000000e-15 * tc[4];
        /*species 10: CH2O */
        species[10] =
            -1.532037000000000e+04 * invT
            -3.916966000000000e+00
            -2.995606000000000e+00 * tc[0]
            -3.340660500000000e-03 * tc[1]
            +4.381591666666666e-07 * tc[2]
            -3.947627500000000e-11 * tc[3]
            +1.606258500000000e-15 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +1.643781000000000e+04 * invT
            -2.608645000000000e+00
            -2.844052000000000e+00 * tc[0]
            -3.068987000000000e-03 * tc[1]
            +3.717241666666666e-07 * tc[2]
            -3.154300833333333e-11 * tc[3]
            +1.226079500000000e-15 * tc[4];
        /*species 12: CO */
        species[12] =
            -1.426835000000000e+04 * invT
            -3.083140000000000e+00
            -3.025078000000000e+00 * tc[0]
            -7.213445000000000e-04 * tc[1]
            +9.384713333333334e-08 * tc[2]
            -8.488174999999999e-12 * tc[3]
            +3.455476000000000e-16 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +1.278325000000000e+02 * invT
            +8.412250000000001e-01
            -3.770800000000000e+00 * tc[0]
            -3.935748500000000e-03 * tc[1]
            +4.427306666666667e-07 * tc[2]
            -3.287025833333333e-11 * tc[3]
            +1.056308000000000e-15 * tc[4];
        /*species 14: C2H5 */
        species[14] =
            +1.067455000000000e+04 * invT
            +2.197137000000000e+01
            -7.190480000000000e+00 * tc[0]
            -3.242038500000000e-03 * tc[1]
            +1.071344166666667e-07 * tc[2]
            +1.956565833333333e-11 * tc[3]
            -1.940438500000000e-15 * tc[4];
        /*species 15: CH4 */
        species[15] =
            -1.008079000000000e+04 * invT
            -7.939916000000000e+00
            -1.683479000000000e+00 * tc[0]
            -5.118620000000000e-03 * tc[1]
            +6.458548333333333e-07 * tc[2]
            -5.654654166666667e-11 * tc[3]
            +2.251711500000000e-15 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +4.428289000000000e+03 * invT
            +1.298030000000000e+00
            -3.528419000000000e+00 * tc[0]
            -5.742590000000000e-03 * tc[1]
            +7.363975000000000e-07 * tc[2]
            -6.537167500000001e-11 * tc[3]
            +2.633424000000000e-15 * tc[4];
        /*species 18: CO2 */
        species[18] =
            -4.896696000000000e+04 * invT
            +5.409018900000000e+00
            -4.453623000000000e+00 * tc[0]
            -1.570084500000000e-03 * tc[1]
            +2.130685000000000e-07 * tc[2]
            -1.994997500000000e-11 * tc[3]
            +8.345165000000000e-16 * tc[4];
        /*species 19: HCO */
        species[19] =
            +3.916324000000000e+03 * invT
            -1.995028000000000e+00
            -3.557271000000000e+00 * tc[0]
            -1.672786500000000e-03 * tc[1]
            +2.225010000000000e-07 * tc[2]
            -2.058810833333333e-11 * tc[3]
            +8.569255000000000e-16 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.566766000000000e+04 * invT
            +7.237108000000000e+00
            -4.436770000000000e+00 * tc[0]
            -2.688019500000000e-03 * tc[1]
            +3.188028333333333e-07 * tc[2]
            -2.738649166666667e-11 * tc[3]
            +1.078355000000000e-15 * tc[4];
        /*species 23: HCCO */
        species[23] =
            +1.901513000000000e+04 * invT
            +1.582933500000000e+01
            -6.758073000000000e+00 * tc[0]
            -1.000200000000000e-03 * tc[1]
            +3.379345000000000e-08 * tc[2]
            +8.676100000000000e-12 * tc[3]
            -9.825825000000000e-16 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +3.185435000000000e+04 * invT
            +1.446378100000000e+01
            -5.933468000000000e+00 * tc[0]
            -2.008873000000000e-03 * tc[1]
            +6.611233333333333e-08 * tc[2]
            +1.201055833333333e-11 * tc[3]
            -1.189322000000000e-15 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +4.903218000000000e+02 * invT
            +1.102092100000000e+01
            -5.975670000000000e+00 * tc[0]
            -4.065295500000000e-03 * tc[1]
            +4.572706666666667e-07 * tc[2]
            -3.391920000000000e-11 * tc[3]
            +1.088008500000000e-15 * tc[4];
        /*species 29: IXC3H7 */
        species[29] =
            +5.313871000000000e+03 * invT
            +2.998983900000000e+01
            -8.063369000000000e+00 * tc[0]
            -7.872440000000000e-03 * tc[1]
            +8.637320000000000e-07 * tc[2]
            -6.231037499999999e-11 * tc[3]
            +1.927211000000000e-15 * tc[4];
        /*species 30: NXC3H7 */
        species[30] =
            +7.579402000000000e+03 * invT
            +2.733440100000000e+01
            -7.978291000000000e+00 * tc[0]
            -7.880565000000001e-03 * tc[1]
            +8.622071666666666e-07 * tc[2]
            -6.203243333333333e-11 * tc[3]
            +1.912489000000000e-15 * tc[4];
        /*species 31: C3H8 */
        species[31] =
            -1.646455000000000e+04 * invT
            +2.536911700000000e+01
            -7.525217000000000e+00 * tc[0]
            -9.445169999999999e-03 * tc[1]
            +1.047320666666667e-06 * tc[2]
            -7.649477500000000e-11 * tc[3]
            +2.406205000000000e-15 * tc[4];
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 3: OH */
        species[3] =
            +3.743212520000000e+03 * invT
            +9.857209199999999e-01
            -3.435862190000000e+00 * tc[0]
            -1.011179020000000e-04 * tc[1]
            +1.892440200000000e-08 * tc[2]
            -2.020376241666666e-11 * tc[3]
            +3.718255155000000e-15 * tc[4];
    } else {
        /*species 3: OH */
        species[3] =
            +4.120853740000000e+03 * invT
            -4.480675530000000e+00
            -2.625997540000000e+00 * tc[0]
            -6.599620300000000e-04 * tc[1]
            +5.995411166666667e-08 * tc[2]
            -3.546923333333333e-12 * tc[3]
            +9.102400800000000e-17 * tc[4];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 8: HO2 */
        species[8] =
            +8.091810130000000e+02 * invT
            -5.210604430000000e+00
            -3.183106560000000e+00 * tc[0]
            -1.833839750000000e-03 * tc[1]
            +1.553975203333333e-07 * tc[2]
            +2.715440991666667e-11 * tc[3]
            -7.556995600000000e-15 * tc[4];
        /*species 21: CH3O2H */
        species[21] =
            -1.771979260000000e+04 * invT
            -6.021811320000000e+00
            -3.234428170000000e+00 * tc[0]
            -9.506488350000000e-03 * tc[1]
            +1.889771450000000e-06 * tc[2]
            -2.835888775000000e-10 * tc[3]
            +2.059151110000000e-14 * tc[4];
    } else {
        /*species 8: HO2 */
        species[8] =
            +3.981276890000000e+02 * invT
            +9.803158699999996e-01
            -4.105474230000000e+00 * tc[0]
            -1.192264175000000e-03 * tc[1]
            +1.343913315000000e-07 * tc[2]
            -1.034931025000000e-11 * tc[3]
            +3.582000540000000e-16 * tc[4];
        /*species 21: CH3O2H */
        species[21] =
            -1.966787710000000e+04 * invT
            +2.754823381000000e+01
            -8.431170910000001e+00 * tc[0]
            -4.034089545000000e-03 * tc[1]
            +4.618248683333333e-07 * tc[2]
            -3.594435358333333e-11 * tc[3]
            +1.253460730000000e-15 * tc[4];
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 17: C2H6 */
        species[17] =
            -1.123455340000000e+04 * invT
            -2.119016043440000e+01
            +2.528543440000000e-02 * tc[0]
            -1.203823770000000e-02 * tc[1]
            +1.864891200000000e-06 * tc[2]
            -1.736174175000000e-10 * tc[3]
            +2.649343080000000e-15 * tc[4];
        /*species 33: NXC3H7O2 */
        species[33] =
            -7.937455670000000e+03 * invT
            -1.680095988000000e+01
            -2.107314920000000e+00 * tc[0]
            -1.980824930000000e-02 * tc[1]
            +4.158193316666667e-06 * tc[2]
            -7.162085833333334e-10 * tc[3]
            +6.562016500000000e-14 * tc[4];
    } else {
        /*species 17: C2H6 */
        species[17] =
            -1.375000140000000e+04 * invT
            +1.911495885000000e+01
            -6.106833850000000e+00 * tc[0]
            -6.461818050000000e-03 * tc[1]
            +7.375453266666666e-07 * tc[2]
            -5.728264383333333e-11 * tc[3]
            +1.994508660000000e-15 * tc[4];
        /*species 33: NXC3H7O2 */
        species[33] =
            -1.191946520000000e+04 * invT
            +5.116763560000000e+01
            -1.263270590000000e+01 * tc[0]
            -8.495536300000000e-03 * tc[1]
            +9.814447883333334e-07 * tc[2]
            -7.684961633333334e-11 * tc[3]
            +2.691153375000000e-15 * tc[4];
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 20: CH3O2 */
        species[20] =
            -6.843942590000000e+02 * invT
            -9.018341400000001e-01
            -4.261469060000000e+00 * tc[0]
            -5.043679950000000e-03 * tc[1]
            +5.358436400000000e-07 * tc[2]
            -1.745077225000000e-11 * tc[3]
            -2.091695515000000e-15 * tc[4];
    } else {
        /*species 20: CH3O2 */
        species[20] =
            -1.535748380000000e+03 * invT
            +1.067751777000000e+01
            -5.957878910000000e+00 * tc[0]
            -3.953643130000000e-03 * tc[1]
            +4.470770566666667e-07 * tc[2]
            -3.449094475000000e-11 * tc[3]
            +1.195036650000000e-15 * tc[4];
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 26: C3H5XA */
        species[26] =
            +1.938342260000000e+04 * invT
            -2.583584505800000e+01
            +5.291319580000000e-01 * tc[0]
            -1.672795500000000e-02 * tc[1]
            +4.223350450000001e-06 * tc[2]
            -8.572146166666666e-10 * tc[3]
            +8.662917000000000e-14 * tc[4];
    } else {
        /*species 26: C3H5XA */
        species[26] =
            +1.635760920000000e+04 * invT
            +3.103978458000000e+01
            -8.458839579999999e+00 * tc[0]
            -5.634774150000000e-03 * tc[1]
            +6.396547733333333e-07 * tc[2]
            -4.950492658333333e-11 * tc[3]
            +1.719590150000000e-15 * tc[4];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 27: C3H6 */
        species[27] =
            +1.066881640000000e+03 * invT
            -2.150575815600000e+01
            -3.946154440000000e-01 * tc[0]
            -1.445538310000000e-02 * tc[1]
            +2.581446800000000e-06 * tc[2]
            -3.240118408333333e-10 * tc[3]
            +1.689451760000000e-14 * tc[4];
        /*species 32: IXC3H7O2 */
        species[32] =
            -1.025879800000000e+04 * invT
            -1.873423261000000e+01
            -1.499416390000000e+00 * tc[0]
            -2.215406025000000e-02 * tc[1]
            +5.373574266666667e-06 * tc[2]
            -1.080726133333333e-09 * tc[3]
            +1.116852845000000e-13 * tc[4];
    } else {
        /*species 27: C3H6 */
        species[27] =
            -1.878212710000000e+03 * invT
            +2.803202638000000e+01
            -8.015959580000001e+00 * tc[0]
            -6.851181700000000e-03 * tc[1]
            +7.770828883333333e-07 * tc[2]
            -6.010453350000000e-11 * tc[3]
            +2.086850630000000e-15 * tc[4];
        /*species 32: IXC3H7O2 */
        species[32] =
            -1.441098550000000e+04 * invT
            +5.615597060000000e+01
            -1.324934930000000e+01 * tc[0]
            -8.204109500000001e-03 * tc[1]
            +9.457201033333333e-07 * tc[2]
            -7.394469499999999e-11 * tc[3]
            +2.586807675000000e-15 * tc[4];
    }

    /*species with midpoint at T=1380 kelvin */
    if (T < 1380) {
        /*species 28: C3H5O */
        species[28] =
            +9.582177840000000e+03 * invT
            -2.035839628000000e+01
            -1.198225820000000e+00 * tc[0]
            -1.527899185000000e-02 * tc[1]
            +3.010504600000000e-06 * tc[2]
            -4.051250275000000e-10 * tc[3]
            +2.099272810000000e-14 * tc[4];
    } else {
        /*species 28: C3H5O */
        species[28] =
            +6.265608100000000e+03 * invT
            +3.802067940000000e+01
            -1.025517520000000e+01 * tc[0]
            -5.749186000000000e-03 * tc[1]
            +6.410760983333333e-07 * tc[2]
            -4.907586216666667e-11 * tc[3]
            +1.692789615000000e-15 * tc[4];
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
        /*species 0: N2 */
        species[0] =
            -1.02090000e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
        /*species 1: O */
        species[1] =
            +2.91476400e+04 * invT
            -1.01756600e+00
            -2.94642900e+00 * tc[0]
            +8.19083000e-04 * tc[1]
            -4.03505333e-07 * tc[2]
            +1.33570250e-10 * tc[3]
            -1.94534800e-14 * tc[4];
        /*species 2: H */
        species[2] =
            +2.54747466e+04 * invT
            +1.96638554e+00
            -2.50104422e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 4: H2 */
        species[4] =
            -1.01252100e+03 * invT
            +5.59221800e+00
            -3.29812400e+00 * tc[0]
            -4.12472100e-04 * tc[1]
            +1.35716917e-07 * tc[2]
            +7.89619500e-12 * tc[3]
            -2.06743600e-14 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.00524900e+03 * invT
            -3.82180200e+00
            -3.21293600e+00 * tc[0]
            -5.63743000e-04 * tc[1]
            +9.59358333e-08 * tc[2]
            -1.09489750e-10 * tc[3]
            +4.38427700e-14 * tc[4];
        /*species 6: H2O */
        species[6] =
            -3.02081100e+04 * invT
            -2.03391000e-01
            -3.38684200e+00 * tc[0]
            -1.73749100e-03 * tc[1]
            +1.05911600e-06 * tc[2]
            -5.80715083e-10 * tc[3]
            +1.25329400e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.76631500e+04 * invT
            -4.39660900e+00
            -3.38875400e+00 * tc[0]
            -3.28461300e-03 * tc[1]
            +2.47502167e-08 * tc[2]
            +3.85483833e-10 * tc[3]
            -1.23575750e-13 * tc[4];
        /*species 9: CH2GSG */
        species[9] =
            +4.98936800e+04 * invT
            +2.91373293e+00
            -3.97126500e+00 * tc[0]
            +8.49544500e-05 * tc[1]
            -1.70894833e-07 * tc[2]
            -2.07712583e-10 * tc[3]
            +9.90633000e-14 * tc[4];
        /*species 10: CH2O */
        species[10] =
            -1.48654000e+04 * invT
            -1.31320890e+01
            -1.65273100e+00 * tc[0]
            -6.31572000e-03 * tc[1]
            +3.14694667e-06 * tc[2]
            -1.70835917e-09 * tc[3]
            +4.20661850e-13 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +1.64237800e+04 * invT
            -5.35935100e+00
            -2.43044300e+00 * tc[0]
            -5.56205000e-03 * tc[1]
            +2.80036667e-06 * tc[2]
            -1.35152417e-09 * tc[3]
            +2.93247650e-13 * tc[4];
        /*species 12: CO */
        species[12] =
            -1.43105400e+04 * invT
            -2.58644500e+00
            -3.26245200e+00 * tc[0]
            -7.55970500e-04 * tc[1]
            +6.46959167e-07 * tc[2]
            -4.65162000e-10 * tc[3]
            +1.23747550e-13 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +9.78601100e+02 * invT
            -1.20459760e+01
            -2.10620400e+00 * tc[0]
            -3.60829750e-03 * tc[1]
            -8.89745333e-07 * tc[2]
            +6.14803000e-10 * tc[3]
            -1.03780550e-13 * tc[4];
        /*species 14: C2H5 */
        species[14] =
            +1.28704000e+04 * invT
            -1.04474980e+01
            -2.69070200e+00 * tc[0]
            -4.35956650e-03 * tc[1]
            -7.36639833e-07 * tc[2]
            -7.78225250e-11 * tc[3]
            +1.96388650e-13 * tc[4];
        /*species 15: CH4 */
        species[15] =
            -9.82522900e+03 * invT
            -1.39434485e+01
            -7.78741500e-01 * tc[0]
            -8.73834000e-03 * tc[1]
            +4.63901500e-06 * tc[2]
            -2.54142333e-09 * tc[3]
            +6.11965500e-13 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +5.57304600e+03 * invT
            -2.60729780e+01
            +8.61488000e-01 * tc[0]
            -1.39808150e-02 * tc[1]
            +5.64779500e-06 * tc[2]
            -2.32096000e-09 * tc[3]
            +4.86893950e-13 * tc[4];
        /*species 18: CO2 */
        species[18] =
            -4.83731400e+04 * invT
            -8.91276500e+00
            -2.27572500e+00 * tc[0]
            -4.96103600e-03 * tc[1]
            +1.73485167e-06 * tc[2]
            -5.72223917e-10 * tc[3]
            +1.05864000e-13 * tc[4];
        /*species 19: HCO */
        species[19] =
            +4.15992200e+03 * invT
            -7.08528400e+00
            -2.89833000e+00 * tc[0]
            -3.09957350e-03 * tc[1]
            +1.60384733e-06 * tc[2]
            -9.08187500e-10 * tc[3]
            +2.28744250e-13 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.61244400e+04 * invT
            -7.79181600e+00
            -2.01356200e+00 * tc[0]
            -7.59522500e-03 * tc[1]
            +2.69386500e-06 * tc[2]
            -7.56582667e-10 * tc[3]
            +9.56373000e-14 * tc[4];
        /*species 23: HCCO */
        species[23] =
            +1.96589200e+04 * invT
            +3.56612110e+00
            -5.04796500e+00 * tc[0]
            -2.22673900e-03 * tc[1]
            -3.78047167e-08 * tc[2]
            +1.23507917e-10 * tc[3]
            -1.12537100e-14 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +3.33522500e+04 * invT
            -1.00969240e+01
            -2.45927600e+00 * tc[0]
            -3.68573800e-03 * tc[1]
            -3.51645500e-07 * tc[2]
            +1.10136833e-10 * tc[3]
            +5.92392000e-14 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +1.52147700e+03 * invT
            -7.14922800e+00
            -3.40906200e+00 * tc[0]
            -5.36928500e-03 * tc[1]
            -3.15248667e-07 * tc[2]
            +5.96548583e-10 * tc[3]
            -1.43369250e-13 * tc[4];
        /*species 29: IXC3H7 */
        species[29] =
            +7.53580900e+03 * invT
            -1.22657100e+01
            -1.71330000e+00 * tc[0]
            -1.27130800e-02 * tc[1]
            -2.63468000e-07 * tc[2]
            +1.51773833e-09 * tc[3]
            -4.41385500e-13 * tc[4];
        /*species 30: NXC3H7 */
        species[30] =
            +9.71328100e+03 * invT
            -1.30701730e+01
            -1.92253700e+00 * tc[0]
            -1.23946350e-02 * tc[1]
            -3.01708167e-07 * tc[2]
            +1.48605500e-09 * tc[3]
            -4.29149800e-13 * tc[4];
        /*species 31: C3H8 */
        species[31] =
            -1.39549200e+04 * invT
            -1.94584092e+01
            -8.96920800e-01 * tc[0]
            -1.33449300e-02 * tc[1]
            -9.05237500e-07 * tc[2]
            +1.77166750e-09 * tc[3]
            -4.62166500e-13 * tc[4];
    } else {
        /*species 0: N2 */
        species[0] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 1: O */
        species[1] =
            +2.92308000e+04 * invT
            -3.37824800e+00
            -2.54206000e+00 * tc[0]
            +1.37753100e-05 * tc[1]
            +5.17133833e-10 * tc[2]
            -3.79255583e-13 * tc[3]
            +2.18402600e-17 * tc[4];
        /*species 2: H */
        species[2] =
            +2.54747466e+04 * invT
            +1.96638554e+00
            -2.50104422e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 4: H2 */
        species[4] =
            -8.35034000e+02 * invT
            +3.34653300e+00
            -2.99142300e+00 * tc[0]
            -3.50032200e-04 * tc[1]
            +9.38971500e-09 * tc[2]
            +7.69298167e-13 * tc[3]
            -7.91376000e-17 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.23393000e+03 * invT
            -4.91588000e-01
            -3.69757800e+00 * tc[0]
            -3.06759850e-04 * tc[1]
            +2.09807000e-08 * tc[2]
            -1.47940083e-12 * tc[3]
            +5.68217500e-17 * tc[4];
        /*species 6: H2O */
        species[6] =
            -2.98992100e+04 * invT
            -5.19067100e+00
            -2.67214600e+00 * tc[0]
            -1.52814650e-03 * tc[1]
            +1.45504333e-07 * tc[2]
            -1.00083000e-11 * tc[3]
            +3.19580900e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.80069600e+04 * invT
            +3.07203000e+00
            -4.57316700e+00 * tc[0]
            -2.16806800e-03 * tc[1]
            +2.45781500e-07 * tc[2]
            -1.95742000e-11 * tc[3]
            +7.15827000e-16 * tc[4];
        /*species 9: CH2GSG */
        species[9] =
            +4.98497500e+04 * invT
            +8.66319000e-01
            -3.55288900e+00 * tc[0]
            -1.03339400e-03 * tc[1]
            +3.19019333e-08 * tc[2]
            +9.20560833e-12 * tc[3]
            -1.01067500e-15 * tc[4];
        /*species 10: CH2O */
        species[10] =
            -1.53203700e+04 * invT
            -4.91696600e+00
            -2.99560600e+00 * tc[0]
            -3.34066050e-03 * tc[1]
            +4.38159167e-07 * tc[2]
            -3.94762750e-11 * tc[3]
            +1.60625850e-15 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +1.64378100e+04 * invT
            -3.60864500e+00
            -2.84405200e+00 * tc[0]
            -3.06898700e-03 * tc[1]
            +3.71724167e-07 * tc[2]
            -3.15430083e-11 * tc[3]
            +1.22607950e-15 * tc[4];
        /*species 12: CO */
        species[12] =
            -1.42683500e+04 * invT
            -4.08314000e+00
            -3.02507800e+00 * tc[0]
            -7.21344500e-04 * tc[1]
            +9.38471333e-08 * tc[2]
            -8.48817500e-12 * tc[3]
            +3.45547600e-16 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +1.27832500e+02 * invT
            -1.58775000e-01
            -3.77080000e+00 * tc[0]
            -3.93574850e-03 * tc[1]
            +4.42730667e-07 * tc[2]
            -3.28702583e-11 * tc[3]
            +1.05630800e-15 * tc[4];
        /*species 14: C2H5 */
        species[14] =
            +1.06745500e+04 * invT
            +2.09713700e+01
            -7.19048000e+00 * tc[0]
            -3.24203850e-03 * tc[1]
            +1.07134417e-07 * tc[2]
            +1.95656583e-11 * tc[3]
            -1.94043850e-15 * tc[4];
        /*species 15: CH4 */
        species[15] =
            -1.00807900e+04 * invT
            -8.93991600e+00
            -1.68347900e+00 * tc[0]
            -5.11862000e-03 * tc[1]
            +6.45854833e-07 * tc[2]
            -5.65465417e-11 * tc[3]
            +2.25171150e-15 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +4.42828900e+03 * invT
            +2.98030000e-01
            -3.52841900e+00 * tc[0]
            -5.74259000e-03 * tc[1]
            +7.36397500e-07 * tc[2]
            -6.53716750e-11 * tc[3]
            +2.63342400e-15 * tc[4];
        /*species 18: CO2 */
        species[18] =
            -4.89669600e+04 * invT
            +4.40901890e+00
            -4.45362300e+00 * tc[0]
            -1.57008450e-03 * tc[1]
            +2.13068500e-07 * tc[2]
            -1.99499750e-11 * tc[3]
            +8.34516500e-16 * tc[4];
        /*species 19: HCO */
        species[19] =
            +3.91632400e+03 * invT
            -2.99502800e+00
            -3.55727100e+00 * tc[0]
            -1.67278650e-03 * tc[1]
            +2.22501000e-07 * tc[2]
            -2.05881083e-11 * tc[3]
            +8.56925500e-16 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.56676600e+04 * invT
            +6.23710800e+00
            -4.43677000e+00 * tc[0]
            -2.68801950e-03 * tc[1]
            +3.18802833e-07 * tc[2]
            -2.73864917e-11 * tc[3]
            +1.07835500e-15 * tc[4];
        /*species 23: HCCO */
        species[23] =
            +1.90151300e+04 * invT
            +1.48293350e+01
            -6.75807300e+00 * tc[0]
            -1.00020000e-03 * tc[1]
            +3.37934500e-08 * tc[2]
            +8.67610000e-12 * tc[3]
            -9.82582500e-16 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +3.18543500e+04 * invT
            +1.34637810e+01
            -5.93346800e+00 * tc[0]
            -2.00887300e-03 * tc[1]
            +6.61123333e-08 * tc[2]
            +1.20105583e-11 * tc[3]
            -1.18932200e-15 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +4.90321800e+02 * invT
            +1.00209210e+01
            -5.97567000e+00 * tc[0]
            -4.06529550e-03 * tc[1]
            +4.57270667e-07 * tc[2]
            -3.39192000e-11 * tc[3]
            +1.08800850e-15 * tc[4];
        /*species 29: IXC3H7 */
        species[29] =
            +5.31387100e+03 * invT
            +2.89898390e+01
            -8.06336900e+00 * tc[0]
            -7.87244000e-03 * tc[1]
            +8.63732000e-07 * tc[2]
            -6.23103750e-11 * tc[3]
            +1.92721100e-15 * tc[4];
        /*species 30: NXC3H7 */
        species[30] =
            +7.57940200e+03 * invT
            +2.63344010e+01
            -7.97829100e+00 * tc[0]
            -7.88056500e-03 * tc[1]
            +8.62207167e-07 * tc[2]
            -6.20324333e-11 * tc[3]
            +1.91248900e-15 * tc[4];
        /*species 31: C3H8 */
        species[31] =
            -1.64645500e+04 * invT
            +2.43691170e+01
            -7.52521700e+00 * tc[0]
            -9.44517000e-03 * tc[1]
            +1.04732067e-06 * tc[2]
            -7.64947750e-11 * tc[3]
            +2.40620500e-15 * tc[4];
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 3: OH */
        species[3] =
            +3.74321252e+03 * invT
            -1.42790800e-02
            -3.43586219e+00 * tc[0]
            -1.01117902e-04 * tc[1]
            +1.89244020e-08 * tc[2]
            -2.02037624e-11 * tc[3]
            +3.71825515e-15 * tc[4];
    } else {
        /*species 3: OH */
        species[3] =
            +4.12085374e+03 * invT
            -5.48067553e+00
            -2.62599754e+00 * tc[0]
            -6.59962030e-04 * tc[1]
            +5.99541117e-08 * tc[2]
            -3.54692333e-12 * tc[3]
            +9.10240080e-17 * tc[4];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 8: HO2 */
        species[8] =
            +8.09181013e+02 * invT
            -6.21060443e+00
            -3.18310656e+00 * tc[0]
            -1.83383975e-03 * tc[1]
            +1.55397520e-07 * tc[2]
            +2.71544099e-11 * tc[3]
            -7.55699560e-15 * tc[4];
        /*species 21: CH3O2H */
        species[21] =
            -1.77197926e+04 * invT
            -7.02181132e+00
            -3.23442817e+00 * tc[0]
            -9.50648835e-03 * tc[1]
            +1.88977145e-06 * tc[2]
            -2.83588878e-10 * tc[3]
            +2.05915111e-14 * tc[4];
    } else {
        /*species 8: HO2 */
        species[8] =
            +3.98127689e+02 * invT
            -1.96841300e-02
            -4.10547423e+00 * tc[0]
            -1.19226417e-03 * tc[1]
            +1.34391332e-07 * tc[2]
            -1.03493103e-11 * tc[3]
            +3.58200054e-16 * tc[4];
        /*species 21: CH3O2H */
        species[21] =
            -1.96678771e+04 * invT
            +2.65482338e+01
            -8.43117091e+00 * tc[0]
            -4.03408955e-03 * tc[1]
            +4.61824868e-07 * tc[2]
            -3.59443536e-11 * tc[3]
            +1.25346073e-15 * tc[4];
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 17: C2H6 */
        species[17] =
            -1.12345534e+04 * invT
            -2.21901604e+01
            +2.52854344e-02 * tc[0]
            -1.20382377e-02 * tc[1]
            +1.86489120e-06 * tc[2]
            -1.73617417e-10 * tc[3]
            +2.64934308e-15 * tc[4];
        /*species 33: NXC3H7O2 */
        species[33] =
            -7.93745567e+03 * invT
            -1.78009599e+01
            -2.10731492e+00 * tc[0]
            -1.98082493e-02 * tc[1]
            +4.15819332e-06 * tc[2]
            -7.16208583e-10 * tc[3]
            +6.56201650e-14 * tc[4];
    } else {
        /*species 17: C2H6 */
        species[17] =
            -1.37500014e+04 * invT
            +1.81149589e+01
            -6.10683385e+00 * tc[0]
            -6.46181805e-03 * tc[1]
            +7.37545327e-07 * tc[2]
            -5.72826438e-11 * tc[3]
            +1.99450866e-15 * tc[4];
        /*species 33: NXC3H7O2 */
        species[33] =
            -1.19194652e+04 * invT
            +5.01676356e+01
            -1.26327059e+01 * tc[0]
            -8.49553630e-03 * tc[1]
            +9.81444788e-07 * tc[2]
            -7.68496163e-11 * tc[3]
            +2.69115337e-15 * tc[4];
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 20: CH3O2 */
        species[20] =
            -6.84394259e+02 * invT
            -1.90183414e+00
            -4.26146906e+00 * tc[0]
            -5.04367995e-03 * tc[1]
            +5.35843640e-07 * tc[2]
            -1.74507723e-11 * tc[3]
            -2.09169552e-15 * tc[4];
    } else {
        /*species 20: CH3O2 */
        species[20] =
            -1.53574838e+03 * invT
            +9.67751777e+00
            -5.95787891e+00 * tc[0]
            -3.95364313e-03 * tc[1]
            +4.47077057e-07 * tc[2]
            -3.44909447e-11 * tc[3]
            +1.19503665e-15 * tc[4];
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 26: C3H5XA */
        species[26] =
            +1.93834226e+04 * invT
            -2.68358451e+01
            +5.29131958e-01 * tc[0]
            -1.67279550e-02 * tc[1]
            +4.22335045e-06 * tc[2]
            -8.57214617e-10 * tc[3]
            +8.66291700e-14 * tc[4];
    } else {
        /*species 26: C3H5XA */
        species[26] =
            +1.63576092e+04 * invT
            +3.00397846e+01
            -8.45883958e+00 * tc[0]
            -5.63477415e-03 * tc[1]
            +6.39654773e-07 * tc[2]
            -4.95049266e-11 * tc[3]
            +1.71959015e-15 * tc[4];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 27: C3H6 */
        species[27] =
            +1.06688164e+03 * invT
            -2.25057582e+01
            -3.94615444e-01 * tc[0]
            -1.44553831e-02 * tc[1]
            +2.58144680e-06 * tc[2]
            -3.24011841e-10 * tc[3]
            +1.68945176e-14 * tc[4];
        /*species 32: IXC3H7O2 */
        species[32] =
            -1.02587980e+04 * invT
            -1.97342326e+01
            -1.49941639e+00 * tc[0]
            -2.21540602e-02 * tc[1]
            +5.37357427e-06 * tc[2]
            -1.08072613e-09 * tc[3]
            +1.11685285e-13 * tc[4];
    } else {
        /*species 27: C3H6 */
        species[27] =
            -1.87821271e+03 * invT
            +2.70320264e+01
            -8.01595958e+00 * tc[0]
            -6.85118170e-03 * tc[1]
            +7.77082888e-07 * tc[2]
            -6.01045335e-11 * tc[3]
            +2.08685063e-15 * tc[4];
        /*species 32: IXC3H7O2 */
        species[32] =
            -1.44109855e+04 * invT
            +5.51559706e+01
            -1.32493493e+01 * tc[0]
            -8.20410950e-03 * tc[1]
            +9.45720103e-07 * tc[2]
            -7.39446950e-11 * tc[3]
            +2.58680767e-15 * tc[4];
    }

    /*species with midpoint at T=1380 kelvin */
    if (T < 1380) {
        /*species 28: C3H5O */
        species[28] =
            +9.58217784e+03 * invT
            -2.13583963e+01
            -1.19822582e+00 * tc[0]
            -1.52789918e-02 * tc[1]
            +3.01050460e-06 * tc[2]
            -4.05125028e-10 * tc[3]
            +2.09927281e-14 * tc[4];
    } else {
        /*species 28: C3H5O */
        species[28] =
            +6.26560810e+03 * invT
            +3.70206794e+01
            -1.02551752e+01 * tc[0]
            -5.74918600e-03 * tc[1]
            +6.41076098e-07 * tc[2]
            -4.90758622e-11 * tc[3]
            +1.69278961e-15 * tc[4];
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
        /*species 0: N2 */
        species[0] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
        /*species 1: O */
        species[1] =
            +1.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 2: H */
        species[2] =
            +1.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54747466e+04 * invT;
        /*species 4: H2 */
        species[4] =
            +2.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 5: O2 */
        species[5] =
            +2.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
        /*species 6: H2O */
        species[6] =
            +2.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 7: H2O2 */
        species[7] =
            +2.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 9: CH2GSG */
        species[9] =
            +2.97126500e+00
            -8.49544500e-05 * tc[1]
            +3.41789667e-07 * tc[2]
            +6.23137750e-10 * tc[3]
            -3.96253200e-13 * tc[4]
            +4.98936800e+04 * invT;
        /*species 10: CH2O */
        species[10] =
            +6.52731000e-01
            +6.31572000e-03 * tc[1]
            -6.29389333e-06 * tc[2]
            +5.12507750e-09 * tc[3]
            -1.68264740e-12 * tc[4]
            -1.48654000e+04 * invT;
        /*species 11: CH3 */
        species[11] =
            +1.43044300e+00
            +5.56205000e-03 * tc[1]
            -5.60073333e-06 * tc[2]
            +4.05457250e-09 * tc[3]
            -1.17299060e-12 * tc[4]
            +1.64237800e+04 * invT;
        /*species 12: CO */
        species[12] =
            +2.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 13: CH3O */
        species[13] =
            +1.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122200e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 14: C2H5 */
        species[14] =
            +1.69070200e+00
            +4.35956650e-03 * tc[1]
            +1.47327967e-06 * tc[2]
            +2.33467575e-10 * tc[3]
            -7.85554600e-13 * tc[4]
            +1.28704000e+04 * invT;
        /*species 15: CH4 */
        species[15] =
            -2.21258500e-01
            +8.73834000e-03 * tc[1]
            -9.27803000e-06 * tc[2]
            +7.62427000e-09 * tc[3]
            -2.44786200e-12 * tc[4]
            -9.82522900e+03 * invT;
        /*species 16: C2H4 */
        species[16] =
            -1.86148800e+00
            +1.39808150e-02 * tc[1]
            -1.12955900e-05 * tc[2]
            +6.96288000e-09 * tc[3]
            -1.94757580e-12 * tc[4]
            +5.57304600e+03 * invT;
        /*species 18: CO2 */
        species[18] =
            +1.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
        /*species 19: HCO */
        species[19] =
            +1.89833000e+00
            +3.09957350e-03 * tc[1]
            -3.20769467e-06 * tc[2]
            +2.72456250e-09 * tc[3]
            -9.14977000e-13 * tc[4]
            +4.15992200e+03 * invT;
        /*species 22: C2H2 */
        species[22] =
            +1.01356200e+00
            +7.59522500e-03 * tc[1]
            -5.38773000e-06 * tc[2]
            +2.26974800e-09 * tc[3]
            -3.82549200e-13 * tc[4]
            +2.61244400e+04 * invT;
        /*species 23: HCCO */
        species[23] =
            +4.04796500e+00
            +2.22673900e-03 * tc[1]
            +7.56094333e-08 * tc[2]
            -3.70523750e-10 * tc[3]
            +4.50148400e-14 * tc[4]
            +1.96589200e+04 * invT;
        /*species 24: C2H3 */
        species[24] =
            +1.45927600e+00
            +3.68573800e-03 * tc[1]
            +7.03291000e-07 * tc[2]
            -3.30410500e-10 * tc[3]
            -2.36956800e-13 * tc[4]
            +3.33522500e+04 * invT;
        /*species 25: CH2CHO */
        species[25] =
            +2.40906200e+00
            +5.36928500e-03 * tc[1]
            +6.30497333e-07 * tc[2]
            -1.78964575e-09 * tc[3]
            +5.73477000e-13 * tc[4]
            +1.52147700e+03 * invT;
        /*species 29: IXC3H7 */
        species[29] =
            +7.13300000e-01
            +1.27130800e-02 * tc[1]
            +5.26936000e-07 * tc[2]
            -4.55321500e-09 * tc[3]
            +1.76554200e-12 * tc[4]
            +7.53580900e+03 * invT;
        /*species 30: NXC3H7 */
        species[30] =
            +9.22537000e-01
            +1.23946350e-02 * tc[1]
            +6.03416333e-07 * tc[2]
            -4.45816500e-09 * tc[3]
            +1.71659920e-12 * tc[4]
            +9.71328100e+03 * invT;
        /*species 31: C3H8 */
        species[31] =
            -1.03079200e-01
            +1.33449300e-02 * tc[1]
            +1.81047500e-06 * tc[2]
            -5.31500250e-09 * tc[3]
            +1.84866600e-12 * tc[4]
            -1.39549200e+04 * invT;
    } else {
        /*species 0: N2 */
        species[0] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 1: O */
        species[1] =
            +1.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 2: H */
        species[2] =
            +1.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54747466e+04 * invT;
        /*species 4: H2 */
        species[4] =
            +1.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 5: O2 */
        species[5] =
            +2.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
        /*species 6: H2O */
        species[6] =
            +1.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 9: CH2GSG */
        species[9] =
            +2.55288900e+00
            +1.03339400e-03 * tc[1]
            -6.38038667e-08 * tc[2]
            -2.76168250e-11 * tc[3]
            +4.04270000e-15 * tc[4]
            +4.98497500e+04 * invT;
        /*species 10: CH2O */
        species[10] =
            +1.99560600e+00
            +3.34066050e-03 * tc[1]
            -8.76318333e-07 * tc[2]
            +1.18428825e-10 * tc[3]
            -6.42503400e-15 * tc[4]
            -1.53203700e+04 * invT;
        /*species 11: CH3 */
        species[11] =
            +1.84405200e+00
            +3.06898700e-03 * tc[1]
            -7.43448333e-07 * tc[2]
            +9.46290250e-11 * tc[3]
            -4.90431800e-15 * tc[4]
            +1.64378100e+04 * invT;
        /*species 12: CO */
        species[12] =
            +2.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 13: CH3O */
        species[13] =
            +2.77080000e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832500e+02 * invT;
        /*species 14: C2H5 */
        species[14] =
            +6.19048000e+00
            +3.24203850e-03 * tc[1]
            -2.14268833e-07 * tc[2]
            -5.86969750e-11 * tc[3]
            +7.76175400e-15 * tc[4]
            +1.06745500e+04 * invT;
        /*species 15: CH4 */
        species[15] =
            +6.83479000e-01
            +5.11862000e-03 * tc[1]
            -1.29170967e-06 * tc[2]
            +1.69639625e-10 * tc[3]
            -9.00684600e-15 * tc[4]
            -1.00807900e+04 * invT;
        /*species 16: C2H4 */
        species[16] =
            +2.52841900e+00
            +5.74259000e-03 * tc[1]
            -1.47279500e-06 * tc[2]
            +1.96115025e-10 * tc[3]
            -1.05336960e-14 * tc[4]
            +4.42828900e+03 * invT;
        /*species 18: CO2 */
        species[18] =
            +3.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
        /*species 19: HCO */
        species[19] =
            +2.55727100e+00
            +1.67278650e-03 * tc[1]
            -4.45002000e-07 * tc[2]
            +6.17643250e-11 * tc[3]
            -3.42770200e-15 * tc[4]
            +3.91632400e+03 * invT;
        /*species 22: C2H2 */
        species[22] =
            +3.43677000e+00
            +2.68801950e-03 * tc[1]
            -6.37605667e-07 * tc[2]
            +8.21594750e-11 * tc[3]
            -4.31342000e-15 * tc[4]
            +2.56676600e+04 * invT;
        /*species 23: HCCO */
        species[23] =
            +5.75807300e+00
            +1.00020000e-03 * tc[1]
            -6.75869000e-08 * tc[2]
            -2.60283000e-11 * tc[3]
            +3.93033000e-15 * tc[4]
            +1.90151300e+04 * invT;
        /*species 24: C2H3 */
        species[24] =
            +4.93346800e+00
            +2.00887300e-03 * tc[1]
            -1.32224667e-07 * tc[2]
            -3.60316750e-11 * tc[3]
            +4.75728800e-15 * tc[4]
            +3.18543500e+04 * invT;
        /*species 25: CH2CHO */
        species[25] =
            +4.97567000e+00
            +4.06529550e-03 * tc[1]
            -9.14541333e-07 * tc[2]
            +1.01757600e-10 * tc[3]
            -4.35203400e-15 * tc[4]
            +4.90321800e+02 * invT;
        /*species 29: IXC3H7 */
        species[29] =
            +7.06336900e+00
            +7.87244000e-03 * tc[1]
            -1.72746400e-06 * tc[2]
            +1.86931125e-10 * tc[3]
            -7.70884400e-15 * tc[4]
            +5.31387100e+03 * invT;
        /*species 30: NXC3H7 */
        species[30] =
            +6.97829100e+00
            +7.88056500e-03 * tc[1]
            -1.72441433e-06 * tc[2]
            +1.86097300e-10 * tc[3]
            -7.64995600e-15 * tc[4]
            +7.57940200e+03 * invT;
        /*species 31: C3H8 */
        species[31] =
            +6.52521700e+00
            +9.44517000e-03 * tc[1]
            -2.09464133e-06 * tc[2]
            +2.29484325e-10 * tc[3]
            -9.62482000e-15 * tc[4]
            -1.64645500e+04 * invT;
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 3: OH */
        species[3] =
            +2.43586219e+00
            +1.01117902e-04 * tc[1]
            -3.78488040e-08 * tc[2]
            +6.06112872e-11 * tc[3]
            -1.48730206e-14 * tc[4]
            +3.74321252e+03 * invT;
    } else {
        /*species 3: OH */
        species[3] =
            +1.62599754e+00
            +6.59962030e-04 * tc[1]
            -1.19908223e-07 * tc[2]
            +1.06407700e-11 * tc[3]
            -3.64096032e-16 * tc[4]
            +4.12085374e+03 * invT;
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 8: HO2 */
        species[8] =
            +2.18310656e+00
            +1.83383975e-03 * tc[1]
            -3.10795041e-07 * tc[2]
            -8.14632298e-11 * tc[3]
            +3.02279824e-14 * tc[4]
            +8.09181013e+02 * invT;
        /*species 21: CH3O2H */
        species[21] =
            +2.23442817e+00
            +9.50648835e-03 * tc[1]
            -3.77954290e-06 * tc[2]
            +8.50766632e-10 * tc[3]
            -8.23660444e-14 * tc[4]
            -1.77197926e+04 * invT;
    } else {
        /*species 8: HO2 */
        species[8] =
            +3.10547423e+00
            +1.19226417e-03 * tc[1]
            -2.68782663e-07 * tc[2]
            +3.10479308e-11 * tc[3]
            -1.43280022e-15 * tc[4]
            +3.98127689e+02 * invT;
        /*species 21: CH3O2H */
        species[21] =
            +7.43117091e+00
            +4.03408955e-03 * tc[1]
            -9.23649737e-07 * tc[2]
            +1.07833061e-10 * tc[3]
            -5.01384292e-15 * tc[4]
            -1.96678771e+04 * invT;
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 17: C2H6 */
        species[17] =
            -1.02528543e+00
            +1.20382377e-02 * tc[1]
            -3.72978240e-06 * tc[2]
            +5.20852252e-10 * tc[3]
            -1.05973723e-14 * tc[4]
            -1.12345534e+04 * invT;
        /*species 33: NXC3H7O2 */
        species[33] =
            +1.10731492e+00
            +1.98082493e-02 * tc[1]
            -8.31638663e-06 * tc[2]
            +2.14862575e-09 * tc[3]
            -2.62480660e-13 * tc[4]
            -7.93745567e+03 * invT;
    } else {
        /*species 17: C2H6 */
        species[17] =
            +5.10683385e+00
            +6.46181805e-03 * tc[1]
            -1.47509065e-06 * tc[2]
            +1.71847932e-10 * tc[3]
            -7.97803464e-15 * tc[4]
            -1.37500014e+04 * invT;
        /*species 33: NXC3H7O2 */
        species[33] =
            +1.16327059e+01
            +8.49553630e-03 * tc[1]
            -1.96288958e-06 * tc[2]
            +2.30548849e-10 * tc[3]
            -1.07646135e-14 * tc[4]
            -1.19194652e+04 * invT;
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 20: CH3O2 */
        species[20] =
            +3.26146906e+00
            +5.04367995e-03 * tc[1]
            -1.07168728e-06 * tc[2]
            +5.23523168e-11 * tc[3]
            +8.36678206e-15 * tc[4]
            -6.84394259e+02 * invT;
    } else {
        /*species 20: CH3O2 */
        species[20] =
            +4.95787891e+00
            +3.95364313e-03 * tc[1]
            -8.94154113e-07 * tc[2]
            +1.03472834e-10 * tc[3]
            -4.78014660e-15 * tc[4]
            -1.53574838e+03 * invT;
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 26: C3H5XA */
        species[26] =
            -1.52913196e+00
            +1.67279550e-02 * tc[1]
            -8.44670090e-06 * tc[2]
            +2.57164385e-09 * tc[3]
            -3.46516680e-13 * tc[4]
            +1.93834226e+04 * invT;
    } else {
        /*species 26: C3H5XA */
        species[26] =
            +7.45883958e+00
            +5.63477415e-03 * tc[1]
            -1.27930955e-06 * tc[2]
            +1.48514780e-10 * tc[3]
            -6.87836060e-15 * tc[4]
            +1.63576092e+04 * invT;
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 27: C3H6 */
        species[27] =
            -6.05384556e-01
            +1.44553831e-02 * tc[1]
            -5.16289360e-06 * tc[2]
            +9.72035522e-10 * tc[3]
            -6.75780704e-14 * tc[4]
            +1.06688164e+03 * invT;
        /*species 32: IXC3H7O2 */
        species[32] =
            +4.99416390e-01
            +2.21540602e-02 * tc[1]
            -1.07471485e-05 * tc[2]
            +3.24217840e-09 * tc[3]
            -4.46741138e-13 * tc[4]
            -1.02587980e+04 * invT;
    } else {
        /*species 27: C3H6 */
        species[27] =
            +7.01595958e+00
            +6.85118170e-03 * tc[1]
            -1.55416578e-06 * tc[2]
            +1.80313601e-10 * tc[3]
            -8.34740252e-15 * tc[4]
            -1.87821271e+03 * invT;
        /*species 32: IXC3H7O2 */
        species[32] =
            +1.22493493e+01
            +8.20410950e-03 * tc[1]
            -1.89144021e-06 * tc[2]
            +2.21834085e-10 * tc[3]
            -1.03472307e-14 * tc[4]
            -1.44109855e+04 * invT;
    }

    /*species with midpoint at T=1380 kelvin */
    if (T < 1380) {
        /*species 28: C3H5O */
        species[28] =
            +1.98225820e-01
            +1.52789918e-02 * tc[1]
            -6.02100920e-06 * tc[2]
            +1.21537508e-09 * tc[3]
            -8.39709124e-14 * tc[4]
            +9.58217784e+03 * invT;
    } else {
        /*species 28: C3H5O */
        species[28] =
            +9.25517520e+00
            +5.74918600e-03 * tc[1]
            -1.28215220e-06 * tc[2]
            +1.47227587e-10 * tc[3]
            -6.77115846e-15 * tc[4]
            +6.26560810e+03 * invT;
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
        /*species 0: N2 */
        species[0] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
        /*species 1: O */
        species[1] =
            +2.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 2: H */
        species[2] =
            +2.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54747466e+04 * invT;
        /*species 4: H2 */
        species[4] =
            +3.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 5: O2 */
        species[5] =
            +3.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
        /*species 6: H2O */
        species[6] =
            +3.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 9: CH2GSG */
        species[9] =
            +3.97126500e+00
            -8.49544500e-05 * tc[1]
            +3.41789667e-07 * tc[2]
            +6.23137750e-10 * tc[3]
            -3.96253200e-13 * tc[4]
            +4.98936800e+04 * invT;
        /*species 10: CH2O */
        species[10] =
            +1.65273100e+00
            +6.31572000e-03 * tc[1]
            -6.29389333e-06 * tc[2]
            +5.12507750e-09 * tc[3]
            -1.68264740e-12 * tc[4]
            -1.48654000e+04 * invT;
        /*species 11: CH3 */
        species[11] =
            +2.43044300e+00
            +5.56205000e-03 * tc[1]
            -5.60073333e-06 * tc[2]
            +4.05457250e-09 * tc[3]
            -1.17299060e-12 * tc[4]
            +1.64237800e+04 * invT;
        /*species 12: CO */
        species[12] =
            +3.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 13: CH3O */
        species[13] =
            +2.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122200e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 14: C2H5 */
        species[14] =
            +2.69070200e+00
            +4.35956650e-03 * tc[1]
            +1.47327967e-06 * tc[2]
            +2.33467575e-10 * tc[3]
            -7.85554600e-13 * tc[4]
            +1.28704000e+04 * invT;
        /*species 15: CH4 */
        species[15] =
            +7.78741500e-01
            +8.73834000e-03 * tc[1]
            -9.27803000e-06 * tc[2]
            +7.62427000e-09 * tc[3]
            -2.44786200e-12 * tc[4]
            -9.82522900e+03 * invT;
        /*species 16: C2H4 */
        species[16] =
            -8.61488000e-01
            +1.39808150e-02 * tc[1]
            -1.12955900e-05 * tc[2]
            +6.96288000e-09 * tc[3]
            -1.94757580e-12 * tc[4]
            +5.57304600e+03 * invT;
        /*species 18: CO2 */
        species[18] =
            +2.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
        /*species 19: HCO */
        species[19] =
            +2.89833000e+00
            +3.09957350e-03 * tc[1]
            -3.20769467e-06 * tc[2]
            +2.72456250e-09 * tc[3]
            -9.14977000e-13 * tc[4]
            +4.15992200e+03 * invT;
        /*species 22: C2H2 */
        species[22] =
            +2.01356200e+00
            +7.59522500e-03 * tc[1]
            -5.38773000e-06 * tc[2]
            +2.26974800e-09 * tc[3]
            -3.82549200e-13 * tc[4]
            +2.61244400e+04 * invT;
        /*species 23: HCCO */
        species[23] =
            +5.04796500e+00
            +2.22673900e-03 * tc[1]
            +7.56094333e-08 * tc[2]
            -3.70523750e-10 * tc[3]
            +4.50148400e-14 * tc[4]
            +1.96589200e+04 * invT;
        /*species 24: C2H3 */
        species[24] =
            +2.45927600e+00
            +3.68573800e-03 * tc[1]
            +7.03291000e-07 * tc[2]
            -3.30410500e-10 * tc[3]
            -2.36956800e-13 * tc[4]
            +3.33522500e+04 * invT;
        /*species 25: CH2CHO */
        species[25] =
            +3.40906200e+00
            +5.36928500e-03 * tc[1]
            +6.30497333e-07 * tc[2]
            -1.78964575e-09 * tc[3]
            +5.73477000e-13 * tc[4]
            +1.52147700e+03 * invT;
        /*species 29: IXC3H7 */
        species[29] =
            +1.71330000e+00
            +1.27130800e-02 * tc[1]
            +5.26936000e-07 * tc[2]
            -4.55321500e-09 * tc[3]
            +1.76554200e-12 * tc[4]
            +7.53580900e+03 * invT;
        /*species 30: NXC3H7 */
        species[30] =
            +1.92253700e+00
            +1.23946350e-02 * tc[1]
            +6.03416333e-07 * tc[2]
            -4.45816500e-09 * tc[3]
            +1.71659920e-12 * tc[4]
            +9.71328100e+03 * invT;
        /*species 31: C3H8 */
        species[31] =
            +8.96920800e-01
            +1.33449300e-02 * tc[1]
            +1.81047500e-06 * tc[2]
            -5.31500250e-09 * tc[3]
            +1.84866600e-12 * tc[4]
            -1.39549200e+04 * invT;
    } else {
        /*species 0: N2 */
        species[0] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 1: O */
        species[1] =
            +2.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 2: H */
        species[2] =
            +2.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54747466e+04 * invT;
        /*species 4: H2 */
        species[4] =
            +2.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 5: O2 */
        species[5] =
            +3.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
        /*species 6: H2O */
        species[6] =
            +2.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 9: CH2GSG */
        species[9] =
            +3.55288900e+00
            +1.03339400e-03 * tc[1]
            -6.38038667e-08 * tc[2]
            -2.76168250e-11 * tc[3]
            +4.04270000e-15 * tc[4]
            +4.98497500e+04 * invT;
        /*species 10: CH2O */
        species[10] =
            +2.99560600e+00
            +3.34066050e-03 * tc[1]
            -8.76318333e-07 * tc[2]
            +1.18428825e-10 * tc[3]
            -6.42503400e-15 * tc[4]
            -1.53203700e+04 * invT;
        /*species 11: CH3 */
        species[11] =
            +2.84405200e+00
            +3.06898700e-03 * tc[1]
            -7.43448333e-07 * tc[2]
            +9.46290250e-11 * tc[3]
            -4.90431800e-15 * tc[4]
            +1.64378100e+04 * invT;
        /*species 12: CO */
        species[12] =
            +3.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 13: CH3O */
        species[13] =
            +3.77080000e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832500e+02 * invT;
        /*species 14: C2H5 */
        species[14] =
            +7.19048000e+00
            +3.24203850e-03 * tc[1]
            -2.14268833e-07 * tc[2]
            -5.86969750e-11 * tc[3]
            +7.76175400e-15 * tc[4]
            +1.06745500e+04 * invT;
        /*species 15: CH4 */
        species[15] =
            +1.68347900e+00
            +5.11862000e-03 * tc[1]
            -1.29170967e-06 * tc[2]
            +1.69639625e-10 * tc[3]
            -9.00684600e-15 * tc[4]
            -1.00807900e+04 * invT;
        /*species 16: C2H4 */
        species[16] =
            +3.52841900e+00
            +5.74259000e-03 * tc[1]
            -1.47279500e-06 * tc[2]
            +1.96115025e-10 * tc[3]
            -1.05336960e-14 * tc[4]
            +4.42828900e+03 * invT;
        /*species 18: CO2 */
        species[18] =
            +4.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
        /*species 19: HCO */
        species[19] =
            +3.55727100e+00
            +1.67278650e-03 * tc[1]
            -4.45002000e-07 * tc[2]
            +6.17643250e-11 * tc[3]
            -3.42770200e-15 * tc[4]
            +3.91632400e+03 * invT;
        /*species 22: C2H2 */
        species[22] =
            +4.43677000e+00
            +2.68801950e-03 * tc[1]
            -6.37605667e-07 * tc[2]
            +8.21594750e-11 * tc[3]
            -4.31342000e-15 * tc[4]
            +2.56676600e+04 * invT;
        /*species 23: HCCO */
        species[23] =
            +6.75807300e+00
            +1.00020000e-03 * tc[1]
            -6.75869000e-08 * tc[2]
            -2.60283000e-11 * tc[3]
            +3.93033000e-15 * tc[4]
            +1.90151300e+04 * invT;
        /*species 24: C2H3 */
        species[24] =
            +5.93346800e+00
            +2.00887300e-03 * tc[1]
            -1.32224667e-07 * tc[2]
            -3.60316750e-11 * tc[3]
            +4.75728800e-15 * tc[4]
            +3.18543500e+04 * invT;
        /*species 25: CH2CHO */
        species[25] =
            +5.97567000e+00
            +4.06529550e-03 * tc[1]
            -9.14541333e-07 * tc[2]
            +1.01757600e-10 * tc[3]
            -4.35203400e-15 * tc[4]
            +4.90321800e+02 * invT;
        /*species 29: IXC3H7 */
        species[29] =
            +8.06336900e+00
            +7.87244000e-03 * tc[1]
            -1.72746400e-06 * tc[2]
            +1.86931125e-10 * tc[3]
            -7.70884400e-15 * tc[4]
            +5.31387100e+03 * invT;
        /*species 30: NXC3H7 */
        species[30] =
            +7.97829100e+00
            +7.88056500e-03 * tc[1]
            -1.72441433e-06 * tc[2]
            +1.86097300e-10 * tc[3]
            -7.64995600e-15 * tc[4]
            +7.57940200e+03 * invT;
        /*species 31: C3H8 */
        species[31] =
            +7.52521700e+00
            +9.44517000e-03 * tc[1]
            -2.09464133e-06 * tc[2]
            +2.29484325e-10 * tc[3]
            -9.62482000e-15 * tc[4]
            -1.64645500e+04 * invT;
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 3: OH */
        species[3] =
            +3.43586219e+00
            +1.01117902e-04 * tc[1]
            -3.78488040e-08 * tc[2]
            +6.06112872e-11 * tc[3]
            -1.48730206e-14 * tc[4]
            +3.74321252e+03 * invT;
    } else {
        /*species 3: OH */
        species[3] =
            +2.62599754e+00
            +6.59962030e-04 * tc[1]
            -1.19908223e-07 * tc[2]
            +1.06407700e-11 * tc[3]
            -3.64096032e-16 * tc[4]
            +4.12085374e+03 * invT;
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 8: HO2 */
        species[8] =
            +3.18310656e+00
            +1.83383975e-03 * tc[1]
            -3.10795041e-07 * tc[2]
            -8.14632298e-11 * tc[3]
            +3.02279824e-14 * tc[4]
            +8.09181013e+02 * invT;
        /*species 21: CH3O2H */
        species[21] =
            +3.23442817e+00
            +9.50648835e-03 * tc[1]
            -3.77954290e-06 * tc[2]
            +8.50766632e-10 * tc[3]
            -8.23660444e-14 * tc[4]
            -1.77197926e+04 * invT;
    } else {
        /*species 8: HO2 */
        species[8] =
            +4.10547423e+00
            +1.19226417e-03 * tc[1]
            -2.68782663e-07 * tc[2]
            +3.10479308e-11 * tc[3]
            -1.43280022e-15 * tc[4]
            +3.98127689e+02 * invT;
        /*species 21: CH3O2H */
        species[21] =
            +8.43117091e+00
            +4.03408955e-03 * tc[1]
            -9.23649737e-07 * tc[2]
            +1.07833061e-10 * tc[3]
            -5.01384292e-15 * tc[4]
            -1.96678771e+04 * invT;
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 17: C2H6 */
        species[17] =
            -2.52854344e-02
            +1.20382377e-02 * tc[1]
            -3.72978240e-06 * tc[2]
            +5.20852252e-10 * tc[3]
            -1.05973723e-14 * tc[4]
            -1.12345534e+04 * invT;
        /*species 33: NXC3H7O2 */
        species[33] =
            +2.10731492e+00
            +1.98082493e-02 * tc[1]
            -8.31638663e-06 * tc[2]
            +2.14862575e-09 * tc[3]
            -2.62480660e-13 * tc[4]
            -7.93745567e+03 * invT;
    } else {
        /*species 17: C2H6 */
        species[17] =
            +6.10683385e+00
            +6.46181805e-03 * tc[1]
            -1.47509065e-06 * tc[2]
            +1.71847932e-10 * tc[3]
            -7.97803464e-15 * tc[4]
            -1.37500014e+04 * invT;
        /*species 33: NXC3H7O2 */
        species[33] =
            +1.26327059e+01
            +8.49553630e-03 * tc[1]
            -1.96288958e-06 * tc[2]
            +2.30548849e-10 * tc[3]
            -1.07646135e-14 * tc[4]
            -1.19194652e+04 * invT;
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 20: CH3O2 */
        species[20] =
            +4.26146906e+00
            +5.04367995e-03 * tc[1]
            -1.07168728e-06 * tc[2]
            +5.23523168e-11 * tc[3]
            +8.36678206e-15 * tc[4]
            -6.84394259e+02 * invT;
    } else {
        /*species 20: CH3O2 */
        species[20] =
            +5.95787891e+00
            +3.95364313e-03 * tc[1]
            -8.94154113e-07 * tc[2]
            +1.03472834e-10 * tc[3]
            -4.78014660e-15 * tc[4]
            -1.53574838e+03 * invT;
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 26: C3H5XA */
        species[26] =
            -5.29131958e-01
            +1.67279550e-02 * tc[1]
            -8.44670090e-06 * tc[2]
            +2.57164385e-09 * tc[3]
            -3.46516680e-13 * tc[4]
            +1.93834226e+04 * invT;
    } else {
        /*species 26: C3H5XA */
        species[26] =
            +8.45883958e+00
            +5.63477415e-03 * tc[1]
            -1.27930955e-06 * tc[2]
            +1.48514780e-10 * tc[3]
            -6.87836060e-15 * tc[4]
            +1.63576092e+04 * invT;
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 27: C3H6 */
        species[27] =
            +3.94615444e-01
            +1.44553831e-02 * tc[1]
            -5.16289360e-06 * tc[2]
            +9.72035522e-10 * tc[3]
            -6.75780704e-14 * tc[4]
            +1.06688164e+03 * invT;
        /*species 32: IXC3H7O2 */
        species[32] =
            +1.49941639e+00
            +2.21540602e-02 * tc[1]
            -1.07471485e-05 * tc[2]
            +3.24217840e-09 * tc[3]
            -4.46741138e-13 * tc[4]
            -1.02587980e+04 * invT;
    } else {
        /*species 27: C3H6 */
        species[27] =
            +8.01595958e+00
            +6.85118170e-03 * tc[1]
            -1.55416578e-06 * tc[2]
            +1.80313601e-10 * tc[3]
            -8.34740252e-15 * tc[4]
            -1.87821271e+03 * invT;
        /*species 32: IXC3H7O2 */
        species[32] =
            +1.32493493e+01
            +8.20410950e-03 * tc[1]
            -1.89144021e-06 * tc[2]
            +2.21834085e-10 * tc[3]
            -1.03472307e-14 * tc[4]
            -1.44109855e+04 * invT;
    }

    /*species with midpoint at T=1380 kelvin */
    if (T < 1380) {
        /*species 28: C3H5O */
        species[28] =
            +1.19822582e+00
            +1.52789918e-02 * tc[1]
            -6.02100920e-06 * tc[2]
            +1.21537508e-09 * tc[3]
            -8.39709124e-14 * tc[4]
            +9.58217784e+03 * invT;
    } else {
        /*species 28: C3H5O */
        species[28] =
            +1.02551752e+01
            +5.74918600e-03 * tc[1]
            -1.28215220e-06 * tc[2]
            +1.47227587e-10 * tc[3]
            -6.77115846e-15 * tc[4]
            +6.26560810e+03 * invT;
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
        /*species 0: N2 */
        species[0] =
            +3.29867700e+00 * tc[0]
            +1.40824000e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213750e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 1: O */
        species[1] =
            +2.94642900e+00 * tc[0]
            -1.63816600e-03 * tc[1]
            +1.21051600e-06 * tc[2]
            -5.34281000e-10 * tc[3]
            +9.72674000e-14 * tc[4]
            +2.96399500e+00 ;
        /*species 2: H */
        species[2] =
            +2.50104422e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.65341317e-01 ;
        /*species 4: H2 */
        species[4] =
            +3.29812400e+00 * tc[0]
            +8.24944200e-04 * tc[1]
            -4.07150750e-07 * tc[2]
            -3.15847800e-11 * tc[3]
            +1.03371800e-13 * tc[4]
            -3.29409400e+00 ;
        /*species 5: O2 */
        species[5] =
            +3.21293600e+00 * tc[0]
            +1.12748600e-03 * tc[1]
            -2.87807500e-07 * tc[2]
            +4.37959000e-10 * tc[3]
            -2.19213850e-13 * tc[4]
            +6.03473800e+00 ;
        /*species 6: H2O */
        species[6] =
            +3.38684200e+00 * tc[0]
            +3.47498200e-03 * tc[1]
            -3.17734800e-06 * tc[2]
            +2.32286033e-09 * tc[3]
            -6.26647000e-13 * tc[4]
            +2.59023300e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00 * tc[0]
            +6.56922600e-03 * tc[1]
            -7.42506500e-08 * tc[2]
            -1.54193533e-09 * tc[3]
            +6.17878750e-13 * tc[4]
            +6.78536300e+00 ;
        /*species 9: CH2GSG */
        species[9] =
            +3.97126500e+00 * tc[0]
            -1.69908900e-04 * tc[1]
            +5.12684500e-07 * tc[2]
            +8.30850333e-10 * tc[3]
            -4.95316500e-13 * tc[4]
            +5.75320700e-02 ;
        /*species 10: CH2O */
        species[10] =
            +1.65273100e+00 * tc[0]
            +1.26314400e-02 * tc[1]
            -9.44084000e-06 * tc[2]
            +6.83343667e-09 * tc[3]
            -2.10330925e-12 * tc[4]
            +1.37848200e+01 ;
        /*species 11: CH3 */
        species[11] =
            +2.43044300e+00 * tc[0]
            +1.11241000e-02 * tc[1]
            -8.40110000e-06 * tc[2]
            +5.40609667e-09 * tc[3]
            -1.46623825e-12 * tc[4]
            +6.78979400e+00 ;
        /*species 12: CO */
        species[12] =
            +3.26245200e+00 * tc[0]
            +1.51194100e-03 * tc[1]
            -1.94087750e-06 * tc[2]
            +1.86064800e-09 * tc[3]
            -6.18737750e-13 * tc[4]
            +4.84889700e+00 ;
        /*species 13: CH3O */
        species[13] =
            +2.10620400e+00 * tc[0]
            +7.21659500e-03 * tc[1]
            +2.66923600e-06 * tc[2]
            -2.45921200e-09 * tc[3]
            +5.18902750e-13 * tc[4]
            +1.31521800e+01 ;
        /*species 14: C2H5 */
        species[14] =
            +2.69070200e+00 * tc[0]
            +8.71913300e-03 * tc[1]
            +2.20991950e-06 * tc[2]
            +3.11290100e-10 * tc[3]
            -9.81943250e-13 * tc[4]
            +1.21382000e+01 ;
        /*species 15: CH4 */
        species[15] =
            +7.78741500e-01 * tc[0]
            +1.74766800e-02 * tc[1]
            -1.39170450e-05 * tc[2]
            +1.01656933e-08 * tc[3]
            -3.05982750e-12 * tc[4]
            +1.37221900e+01 ;
        /*species 16: C2H4 */
        species[16] =
            -8.61488000e-01 * tc[0]
            +2.79616300e-02 * tc[1]
            -1.69433850e-05 * tc[2]
            +9.28384000e-09 * tc[3]
            -2.43446975e-12 * tc[4]
            +2.42114900e+01 ;
        /*species 18: CO2 */
        species[18] =
            +2.27572500e+00 * tc[0]
            +9.92207200e-03 * tc[1]
            -5.20455500e-06 * tc[2]
            +2.28889567e-09 * tc[3]
            -5.29320000e-13 * tc[4]
            +1.01884900e+01 ;
        /*species 19: HCO */
        species[19] =
            +2.89833000e+00 * tc[0]
            +6.19914700e-03 * tc[1]
            -4.81154200e-06 * tc[2]
            +3.63275000e-09 * tc[3]
            -1.14372125e-12 * tc[4]
            +8.98361400e+00 ;
        /*species 22: C2H2 */
        species[22] =
            +2.01356200e+00 * tc[0]
            +1.51904500e-02 * tc[1]
            -8.08159500e-06 * tc[2]
            +3.02633067e-09 * tc[3]
            -4.78186500e-13 * tc[4]
            +8.80537800e+00 ;
        /*species 23: HCCO */
        species[23] =
            +5.04796500e+00 * tc[0]
            +4.45347800e-03 * tc[1]
            +1.13414150e-07 * tc[2]
            -4.94031667e-10 * tc[3]
            +5.62685500e-14 * tc[4]
            +4.81843900e-01 ;
        /*species 24: C2H3 */
        species[24] =
            +2.45927600e+00 * tc[0]
            +7.37147600e-03 * tc[1]
            +1.05493650e-06 * tc[2]
            -4.40547333e-10 * tc[3]
            -2.96196000e-13 * tc[4]
            +1.15562000e+01 ;
        /*species 25: CH2CHO */
        species[25] =
            +3.40906200e+00 * tc[0]
            +1.07385700e-02 * tc[1]
            +9.45746000e-07 * tc[2]
            -2.38619433e-09 * tc[3]
            +7.16846250e-13 * tc[4]
            +9.55829000e+00 ;
        /*species 29: IXC3H7 */
        species[29] =
            +1.71330000e+00 * tc[0]
            +2.54261600e-02 * tc[1]
            +7.90404000e-07 * tc[2]
            -6.07095333e-09 * tc[3]
            +2.20692750e-12 * tc[4]
            +1.29790100e+01 ;
        /*species 30: NXC3H7 */
        species[30] =
            +1.92253700e+00 * tc[0]
            +2.47892700e-02 * tc[1]
            +9.05124500e-07 * tc[2]
            -5.94422000e-09 * tc[3]
            +2.14574900e-12 * tc[4]
            +1.39927100e+01 ;
        /*species 31: C3H8 */
        species[31] =
            +8.96920800e-01 * tc[0]
            +2.66898600e-02 * tc[1]
            +2.71571250e-06 * tc[2]
            -7.08667000e-09 * tc[3]
            +2.31083250e-12 * tc[4]
            +1.93553300e+01 ;
    } else {
        /*species 0: N2 */
        species[0] =
            +2.92664000e+00 * tc[0]
            +1.48797700e-03 * tc[1]
            -2.84238050e-07 * tc[2]
            +3.36568000e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 1: O */
        species[1] =
            +2.54206000e+00 * tc[0]
            -2.75506200e-05 * tc[1]
            -1.55140150e-09 * tc[2]
            +1.51702233e-12 * tc[3]
            -1.09201300e-16 * tc[4]
            +4.92030800e+00 ;
        /*species 2: H */
        species[2] =
            +2.50104422e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.65341317e-01 ;
        /*species 4: H2 */
        species[4] =
            +2.99142300e+00 * tc[0]
            +7.00064400e-04 * tc[1]
            -2.81691450e-08 * tc[2]
            -3.07719267e-12 * tc[3]
            +3.95688000e-16 * tc[4]
            -1.35511000e+00 ;
        /*species 5: O2 */
        species[5] =
            +3.69757800e+00 * tc[0]
            +6.13519700e-04 * tc[1]
            -6.29421000e-08 * tc[2]
            +5.91760333e-12 * tc[3]
            -2.84108750e-16 * tc[4]
            +3.18916600e+00 ;
        /*species 6: H2O */
        species[6] =
            +2.67214600e+00 * tc[0]
            +3.05629300e-03 * tc[1]
            -4.36513000e-07 * tc[2]
            +4.00332000e-11 * tc[3]
            -1.59790450e-15 * tc[4]
            +6.86281700e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00 * tc[0]
            +4.33613600e-03 * tc[1]
            -7.37344500e-07 * tc[2]
            +7.82968000e-11 * tc[3]
            -3.57913500e-15 * tc[4]
            +5.01137000e-01 ;
        /*species 9: CH2GSG */
        species[9] =
            +3.55288900e+00 * tc[0]
            +2.06678800e-03 * tc[1]
            -9.57058000e-08 * tc[2]
            -3.68224333e-11 * tc[3]
            +5.05337500e-15 * tc[4]
            +1.68657000e+00 ;
        /*species 10: CH2O */
        species[10] =
            +2.99560600e+00 * tc[0]
            +6.68132100e-03 * tc[1]
            -1.31447750e-06 * tc[2]
            +1.57905100e-10 * tc[3]
            -8.03129250e-15 * tc[4]
            +6.91257200e+00 ;
        /*species 11: CH3 */
        species[11] =
            +2.84405200e+00 * tc[0]
            +6.13797400e-03 * tc[1]
            -1.11517250e-06 * tc[2]
            +1.26172033e-10 * tc[3]
            -6.13039750e-15 * tc[4]
            +5.45269700e+00 ;
        /*species 12: CO */
        species[12] =
            +3.02507800e+00 * tc[0]
            +1.44268900e-03 * tc[1]
            -2.81541400e-07 * tc[2]
            +3.39527000e-11 * tc[3]
            -1.72773800e-15 * tc[4]
            +6.10821800e+00 ;
        /*species 13: CH3O */
        species[13] =
            +3.77080000e+00 * tc[0]
            +7.87149700e-03 * tc[1]
            -1.32819200e-06 * tc[2]
            +1.31481033e-10 * tc[3]
            -5.28154000e-15 * tc[4]
            +2.92957500e+00 ;
        /*species 14: C2H5 */
        species[14] =
            +7.19048000e+00 * tc[0]
            +6.48407700e-03 * tc[1]
            -3.21403250e-07 * tc[2]
            -7.82626333e-11 * tc[3]
            +9.70219250e-15 * tc[4]
            -1.47808900e+01 ;
        /*species 15: CH4 */
        species[15] =
            +1.68347900e+00 * tc[0]
            +1.02372400e-02 * tc[1]
            -1.93756450e-06 * tc[2]
            +2.26186167e-10 * tc[3]
            -1.12585575e-14 * tc[4]
            +9.62339500e+00 ;
        /*species 16: C2H4 */
        species[16] =
            +3.52841900e+00 * tc[0]
            +1.14851800e-02 * tc[1]
            -2.20919250e-06 * tc[2]
            +2.61486700e-10 * tc[3]
            -1.31671200e-14 * tc[4]
            +2.23038900e+00 ;
        /*species 18: CO2 */
        species[18] =
            +4.45362300e+00 * tc[0]
            +3.14016900e-03 * tc[1]
            -6.39205500e-07 * tc[2]
            +7.97999000e-11 * tc[3]
            -4.17258250e-15 * tc[4]
            -9.55395900e-01 ;
        /*species 19: HCO */
        species[19] =
            +3.55727100e+00 * tc[0]
            +3.34557300e-03 * tc[1]
            -6.67503000e-07 * tc[2]
            +8.23524333e-11 * tc[3]
            -4.28462750e-15 * tc[4]
            +5.55229900e+00 ;
        /*species 22: C2H2 */
        species[22] =
            +4.43677000e+00 * tc[0]
            +5.37603900e-03 * tc[1]
            -9.56408500e-07 * tc[2]
            +1.09545967e-10 * tc[3]
            -5.39177500e-15 * tc[4]
            -2.80033800e+00 ;
        /*species 23: HCCO */
        species[23] =
            +6.75807300e+00 * tc[0]
            +2.00040000e-03 * tc[1]
            -1.01380350e-07 * tc[2]
            -3.47044000e-11 * tc[3]
            +4.91291250e-15 * tc[4]
            -9.07126200e+00 ;
        /*species 24: C2H3 */
        species[24] =
            +5.93346800e+00 * tc[0]
            +4.01774600e-03 * tc[1]
            -1.98337000e-07 * tc[2]
            -4.80422333e-11 * tc[3]
            +5.94661000e-15 * tc[4]
            -8.53031300e+00 ;
        /*species 25: CH2CHO */
        species[25] =
            +5.97567000e+00 * tc[0]
            +8.13059100e-03 * tc[1]
            -1.37181200e-06 * tc[2]
            +1.35676800e-10 * tc[3]
            -5.44004250e-15 * tc[4]
            -5.04525100e+00 ;
        /*species 29: IXC3H7 */
        species[29] =
            +8.06336900e+00 * tc[0]
            +1.57448800e-02 * tc[1]
            -2.59119600e-06 * tc[2]
            +2.49241500e-10 * tc[3]
            -9.63605500e-15 * tc[4]
            -2.19264700e+01 ;
        /*species 30: NXC3H7 */
        species[30] =
            +7.97829100e+00 * tc[0]
            +1.57611300e-02 * tc[1]
            -2.58662150e-06 * tc[2]
            +2.48129733e-10 * tc[3]
            -9.56244500e-15 * tc[4]
            -1.93561100e+01 ;
        /*species 31: C3H8 */
        species[31] =
            +7.52521700e+00 * tc[0]
            +1.88903400e-02 * tc[1]
            -3.14196200e-06 * tc[2]
            +3.05979100e-10 * tc[3]
            -1.20310250e-14 * tc[4]
            -1.78439000e+01 ;
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 3: OH */
        species[3] =
            +3.43586219e+00 * tc[0]
            +2.02235804e-04 * tc[1]
            -5.67732060e-08 * tc[2]
            +8.08150497e-11 * tc[3]
            -1.85912758e-14 * tc[4]
            +2.45014127e+00 ;
    } else {
        /*species 3: OH */
        species[3] =
            +2.62599754e+00 * tc[0]
            +1.31992406e-03 * tc[1]
            -1.79862335e-07 * tc[2]
            +1.41876933e-11 * tc[3]
            -4.55120040e-16 * tc[4]
            +7.10667307e+00 ;
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 8: HO2 */
        species[8] =
            +3.18310656e+00 * tc[0]
            +3.66767950e-03 * tc[1]
            -4.66192561e-07 * tc[2]
            -1.08617640e-10 * tc[3]
            +3.77849780e-14 * tc[4]
            +8.39371099e+00 ;
        /*species 21: CH3O2H */
        species[21] =
            +3.23442817e+00 * tc[0]
            +1.90129767e-02 * tc[1]
            -5.66931435e-06 * tc[2]
            +1.13435551e-09 * tc[3]
            -1.02957555e-13 * tc[4]
            +9.25623949e+00 ;
    } else {
        /*species 8: HO2 */
        species[8] =
            +4.10547423e+00 * tc[0]
            +2.38452835e-03 * tc[1]
            -4.03173995e-07 * tc[2]
            +4.13972410e-11 * tc[3]
            -1.79100027e-15 * tc[4]
            +3.12515836e+00 ;
        /*species 21: CH3O2H */
        species[21] =
            +8.43117091e+00 * tc[0]
            +8.06817909e-03 * tc[1]
            -1.38547461e-06 * tc[2]
            +1.43777414e-10 * tc[3]
            -6.26730365e-15 * tc[4]
            -1.91170629e+01 ;
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 17: C2H6 */
        species[17] =
            -2.52854344e-02 * tc[0]
            +2.40764754e-02 * tc[1]
            -5.59467360e-06 * tc[2]
            +6.94469670e-10 * tc[3]
            -1.32467154e-14 * tc[4]
            +2.11648750e+01 ;
        /*species 33: NXC3H7O2 */
        species[33] =
            +2.10731492e+00 * tc[0]
            +3.96164986e-02 * tc[1]
            -1.24745800e-05 * tc[2]
            +2.86483433e-09 * tc[3]
            -3.28100825e-13 * tc[4]
            +1.89082748e+01 ;
    } else {
        /*species 17: C2H6 */
        species[17] =
            +6.10683385e+00 * tc[0]
            +1.29236361e-02 * tc[1]
            -2.21263598e-06 * tc[2]
            +2.29130575e-10 * tc[3]
            -9.97254330e-15 * tc[4]
            -1.30081250e+01 ;
        /*species 33: NXC3H7O2 */
        species[33] =
            +1.26327059e+01 * tc[0]
            +1.69910726e-02 * tc[1]
            -2.94433436e-06 * tc[2]
            +3.07398465e-10 * tc[3]
            -1.34557669e-14 * tc[4]
            -3.85349297e+01 ;
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 20: CH3O2 */
        species[20] =
            +4.26146906e+00 * tc[0]
            +1.00873599e-02 * tc[1]
            -1.60753092e-06 * tc[2]
            +6.98030890e-11 * tc[3]
            +1.04584776e-14 * tc[4]
            +5.16330320e+00 ;
    } else {
        /*species 20: CH3O2 */
        species[20] =
            +5.95787891e+00 * tc[0]
            +7.90728626e-03 * tc[1]
            -1.34123117e-06 * tc[2]
            +1.37963779e-10 * tc[3]
            -5.97518325e-15 * tc[4]
            -4.71963886e+00 ;
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 26: C3H5XA */
        species[26] =
            -5.29131958e-01 * tc[0]
            +3.34559100e-02 * tc[1]
            -1.26700514e-05 * tc[2]
            +3.42885847e-09 * tc[3]
            -4.33145850e-13 * tc[4]
            +2.53067131e+01 ;
    } else {
        /*species 26: C3H5XA */
        species[26] =
            +8.45883958e+00 * tc[0]
            +1.12695483e-02 * tc[1]
            -1.91896432e-06 * tc[2]
            +1.98019706e-10 * tc[3]
            -8.59795075e-15 * tc[4]
            -2.25809450e+01 ;
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 27: C3H6 */
        species[27] =
            +3.94615444e-01 * tc[0]
            +2.89107662e-02 * tc[1]
            -7.74434040e-06 * tc[2]
            +1.29604736e-09 * tc[3]
            -8.44725880e-14 * tc[4]
            +2.19003736e+01 ;
        /*species 32: IXC3H7O2 */
        species[32] =
            +1.49941639e+00 * tc[0]
            +4.43081205e-02 * tc[1]
            -1.61207228e-05 * tc[2]
            +4.32290453e-09 * tc[3]
            -5.58426423e-13 * tc[4]
            +2.02336490e+01 ;
    } else {
        /*species 27: C3H6 */
        species[27] =
            +8.01595958e+00 * tc[0]
            +1.37023634e-02 * tc[1]
            -2.33124867e-06 * tc[2]
            +2.40418134e-10 * tc[3]
            -1.04342532e-14 * tc[4]
            -2.00160668e+01 ;
        /*species 32: IXC3H7O2 */
        species[32] =
            +1.32493493e+01 * tc[0]
            +1.64082190e-02 * tc[1]
            -2.83716031e-06 * tc[2]
            +2.95778780e-10 * tc[3]
            -1.29340384e-14 * tc[4]
            -4.29066213e+01 ;
    }

    /*species with midpoint at T=1380 kelvin */
    if (T < 1380) {
        /*species 28: C3H5O */
        species[28] =
            +1.19822582e+00 * tc[0]
            +3.05579837e-02 * tc[1]
            -9.03151380e-06 * tc[2]
            +1.62050011e-09 * tc[3]
            -1.04963640e-13 * tc[4]
            +2.15566221e+01 ;
    } else {
        /*species 28: C3H5O */
        species[28] =
            +1.02551752e+01 * tc[0]
            +1.14983720e-02 * tc[1]
            -1.92322829e-06 * tc[2]
            +1.96303449e-10 * tc[3]
            -8.46394807e-15 * tc[4]
            -2.77655042e+01 ;
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
    amrex::Real cpor[34]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 34; ++id) {
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
    amrex::Real cpor[34], tresult[34]; /* temporary storage */
    amrex::Real imw[34];

    get_imw(imw);

    cp_R(cpor, tc);

    for (int i = 0; i < 34; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 34; i++)
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
    amrex::Real cvor[34]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 34; ++id) {
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
    amrex::Real cvor[34]; /* temporary storage */
    amrex::Real imw[34];

    get_imw(imw);

    cv_R(cvor, tc);

    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*N2 */
    result += cvor[1]*y[1]*imw[1]; /*O */
    result += cvor[2]*y[2]*imw[2]; /*H */
    result += cvor[3]*y[3]*imw[3]; /*OH */
    result += cvor[4]*y[4]*imw[4]; /*H2 */
    result += cvor[5]*y[5]*imw[5]; /*O2 */
    result += cvor[6]*y[6]*imw[6]; /*H2O */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*HO2 */
    result += cvor[9]*y[9]*imw[9]; /*CH2GSG */
    result += cvor[10]*y[10]*imw[10]; /*CH2O */
    result += cvor[11]*y[11]*imw[11]; /*CH3 */
    result += cvor[12]*y[12]*imw[12]; /*CO */
    result += cvor[13]*y[13]*imw[13]; /*CH3O */
    result += cvor[14]*y[14]*imw[14]; /*C2H5 */
    result += cvor[15]*y[15]*imw[15]; /*CH4 */
    result += cvor[16]*y[16]*imw[16]; /*C2H4 */
    result += cvor[17]*y[17]*imw[17]; /*C2H6 */
    result += cvor[18]*y[18]*imw[18]; /*CO2 */
    result += cvor[19]*y[19]*imw[19]; /*HCO */
    result += cvor[20]*y[20]*imw[20]; /*CH3O2 */
    result += cvor[21]*y[21]*imw[21]; /*CH3O2H */
    result += cvor[22]*y[22]*imw[22]; /*C2H2 */
    result += cvor[23]*y[23]*imw[23]; /*HCCO */
    result += cvor[24]*y[24]*imw[24]; /*C2H3 */
    result += cvor[25]*y[25]*imw[25]; /*CH2CHO */
    result += cvor[26]*y[26]*imw[26]; /*C3H5XA */
    result += cvor[27]*y[27]*imw[27]; /*C3H6 */
    result += cvor[28]*y[28]*imw[28]; /*C3H5O */
    result += cvor[29]*y[29]*imw[29]; /*IXC3H7 */
    result += cvor[30]*y[30]*imw[30]; /*NXC3H7 */
    result += cvor[31]*y[31]*imw[31]; /*C3H8 */
    result += cvor[32]*y[32]*imw[32]; /*IXC3H7O2 */
    result += cvor[33]*y[33]*imw[33]; /*NXC3H7O2 */

    *cvbs = result * 8.31446261815324e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKHBML(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  hbml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real hml[34]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 34; ++id) {
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
    amrex::Real hml[34], tmp[34]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[34];

    get_imw(imw);

    speciesEnthalpy(hml, tc);

    int id;
    for (id = 0; id < 34; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 34; ++id) {
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
    amrex::Real uml[34]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 34; ++id) {
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
    amrex::Real ums[34]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[34];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*N2 */
    result += y[1]*ums[1]*imw[1]; /*O */
    result += y[2]*ums[2]*imw[2]; /*H */
    result += y[3]*ums[3]*imw[3]; /*OH */
    result += y[4]*ums[4]*imw[4]; /*H2 */
    result += y[5]*ums[5]*imw[5]; /*O2 */
    result += y[6]*ums[6]*imw[6]; /*H2O */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*HO2 */
    result += y[9]*ums[9]*imw[9]; /*CH2GSG */
    result += y[10]*ums[10]*imw[10]; /*CH2O */
    result += y[11]*ums[11]*imw[11]; /*CH3 */
    result += y[12]*ums[12]*imw[12]; /*CO */
    result += y[13]*ums[13]*imw[13]; /*CH3O */
    result += y[14]*ums[14]*imw[14]; /*C2H5 */
    result += y[15]*ums[15]*imw[15]; /*CH4 */
    result += y[16]*ums[16]*imw[16]; /*C2H4 */
    result += y[17]*ums[17]*imw[17]; /*C2H6 */
    result += y[18]*ums[18]*imw[18]; /*CO2 */
    result += y[19]*ums[19]*imw[19]; /*HCO */
    result += y[20]*ums[20]*imw[20]; /*CH3O2 */
    result += y[21]*ums[21]*imw[21]; /*CH3O2H */
    result += y[22]*ums[22]*imw[22]; /*C2H2 */
    result += y[23]*ums[23]*imw[23]; /*HCCO */
    result += y[24]*ums[24]*imw[24]; /*C2H3 */
    result += y[25]*ums[25]*imw[25]; /*CH2CHO */
    result += y[26]*ums[26]*imw[26]; /*C3H5XA */
    result += y[27]*ums[27]*imw[27]; /*C3H6 */
    result += y[28]*ums[28]*imw[28]; /*C3H5O */
    result += y[29]*ums[29]*imw[29]; /*IXC3H7 */
    result += y[30]*ums[30]*imw[30]; /*NXC3H7 */
    result += y[31]*ums[31]*imw[31]; /*C3H8 */
    result += y[32]*ums[32]*imw[32]; /*IXC3H7O2 */
    result += y[33]*ums[33]*imw[33]; /*NXC3H7O2 */

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
    amrex::Real sor[34]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 34; ++id) {
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
    amrex::Real sor[34]; /* temporary storage */
    amrex::Real x[34]; /* need a ytx conversion */
    amrex::Real YOW = 0; /*See Eq 4, 6 in CK Manual */
    amrex::Real imw[34];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*N2 */
    YOW += y[1]*imw[1]; /*O */
    YOW += y[2]*imw[2]; /*H */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2 */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*CH2GSG */
    YOW += y[10]*imw[10]; /*CH2O */
    YOW += y[11]*imw[11]; /*CH3 */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CH3O */
    YOW += y[14]*imw[14]; /*C2H5 */
    YOW += y[15]*imw[15]; /*CH4 */
    YOW += y[16]*imw[16]; /*C2H4 */
    YOW += y[17]*imw[17]; /*C2H6 */
    YOW += y[18]*imw[18]; /*CO2 */
    YOW += y[19]*imw[19]; /*HCO */
    YOW += y[20]*imw[20]; /*CH3O2 */
    YOW += y[21]*imw[21]; /*CH3O2H */
    YOW += y[22]*imw[22]; /*C2H2 */
    YOW += y[23]*imw[23]; /*HCCO */
    YOW += y[24]*imw[24]; /*C2H3 */
    YOW += y[25]*imw[25]; /*CH2CHO */
    YOW += y[26]*imw[26]; /*C3H5XA */
    YOW += y[27]*imw[27]; /*C3H6 */
    YOW += y[28]*imw[28]; /*C3H5O */
    YOW += y[29]*imw[29]; /*IXC3H7 */
    YOW += y[30]*imw[30]; /*NXC3H7 */
    YOW += y[31]*imw[31]; /*C3H8 */
    YOW += y[32]*imw[32]; /*IXC3H7O2 */
    YOW += y[33]*imw[33]; /*NXC3H7O2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(28.013400*YOW); 
    x[1] = y[1]/(15.999400*YOW); 
    x[2] = y[2]/(1.007970*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(2.015940*YOW); 
    x[5] = y[5]/(31.998800*YOW); 
    x[6] = y[6]/(18.015340*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(14.027090*YOW); 
    x[10] = y[10]/(30.026490*YOW); 
    x[11] = y[11]/(15.035060*YOW); 
    x[12] = y[12]/(28.010550*YOW); 
    x[13] = y[13]/(31.034460*YOW); 
    x[14] = y[14]/(29.062150*YOW); 
    x[15] = y[15]/(16.043030*YOW); 
    x[16] = y[16]/(28.054180*YOW); 
    x[17] = y[17]/(30.070120*YOW); 
    x[18] = y[18]/(44.009950*YOW); 
    x[19] = y[19]/(29.018520*YOW); 
    x[20] = y[20]/(47.033860*YOW); 
    x[21] = y[21]/(48.041830*YOW); 
    x[22] = y[22]/(26.038240*YOW); 
    x[23] = y[23]/(41.029670*YOW); 
    x[24] = y[24]/(27.046210*YOW); 
    x[25] = y[25]/(43.045610*YOW); 
    x[26] = y[26]/(41.073300*YOW); 
    x[27] = y[27]/(42.081270*YOW); 
    x[28] = y[28]/(57.072700*YOW); 
    x[29] = y[29]/(43.089240*YOW); 
    x[30] = y[30]/(43.089240*YOW); 
    x[31] = y[31]/(44.097210*YOW); 
    x[32] = y[32]/(75.088040*YOW); 
    x[33] = y[33]/(75.088040*YOW); 
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
    result += x[21]*(sor[21]-log((x[21]+1e-100))-logPratio);
    result += x[22]*(sor[22]-log((x[22]+1e-100))-logPratio);
    result += x[23]*(sor[23]-log((x[23]+1e-100))-logPratio);
    result += x[24]*(sor[24]-log((x[24]+1e-100))-logPratio);
    result += x[25]*(sor[25]-log((x[25]+1e-100))-logPratio);
    result += x[26]*(sor[26]-log((x[26]+1e-100))-logPratio);
    result += x[27]*(sor[27]-log((x[27]+1e-100))-logPratio);
    result += x[28]*(sor[28]-log((x[28]+1e-100))-logPratio);
    result += x[29]*(sor[29]-log((x[29]+1e-100))-logPratio);
    result += x[30]*(sor[30]-log((x[30]+1e-100))-logPratio);
    result += x[31]*(sor[31]-log((x[31]+1e-100))-logPratio);
    result += x[32]*(sor[32]-log((x[32]+1e-100))-logPratio);
    result += x[33]*(sor[33]-log((x[33]+1e-100))-logPratio);
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
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*2.015940; /*H2 */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*14.027090; /*CH2GSG */
    XW += x[10]*30.026490; /*CH2O */
    XW += x[11]*15.035060; /*CH3 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*29.062150; /*C2H5 */
    XW += x[15]*16.043030; /*CH4 */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*30.070120; /*C2H6 */
    XW += x[18]*44.009950; /*CO2 */
    XW += x[19]*29.018520; /*HCO */
    XW += x[20]*47.033860; /*CH3O2 */
    XW += x[21]*48.041830; /*CH3O2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*41.029670; /*HCCO */
    XW += x[24]*27.046210; /*C2H3 */
    XW += x[25]*43.045610; /*CH2CHO */
    XW += x[26]*41.073300; /*C3H5XA */
    XW += x[27]*42.081270; /*C3H6 */
    XW += x[28]*57.072700; /*C3H5O */
    XW += x[29]*43.089240; /*IXC3H7 */
    XW += x[30]*43.089240; /*NXC3H7 */
    XW += x[31]*44.097210; /*C3H8 */
    XW += x[32]*75.088040; /*IXC3H7O2 */
    XW += x[33]*75.088040; /*NXC3H7O2 */
    *P = *rho * 8.31446261815324e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKPY(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW = 0;/* for computing mean MW */
    amrex::Real imw[34];

    get_imw(imw);

    YOW += y[0]*imw[0]; /*N2 */
    YOW += y[1]*imw[1]; /*O */
    YOW += y[2]*imw[2]; /*H */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2 */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*CH2GSG */
    YOW += y[10]*imw[10]; /*CH2O */
    YOW += y[11]*imw[11]; /*CH3 */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CH3O */
    YOW += y[14]*imw[14]; /*C2H5 */
    YOW += y[15]*imw[15]; /*CH4 */
    YOW += y[16]*imw[16]; /*C2H4 */
    YOW += y[17]*imw[17]; /*C2H6 */
    YOW += y[18]*imw[18]; /*CO2 */
    YOW += y[19]*imw[19]; /*HCO */
    YOW += y[20]*imw[20]; /*CH3O2 */
    YOW += y[21]*imw[21]; /*CH3O2H */
    YOW += y[22]*imw[22]; /*C2H2 */
    YOW += y[23]*imw[23]; /*HCCO */
    YOW += y[24]*imw[24]; /*C2H3 */
    YOW += y[25]*imw[25]; /*CH2CHO */
    YOW += y[26]*imw[26]; /*C3H5XA */
    YOW += y[27]*imw[27]; /*C3H6 */
    YOW += y[28]*imw[28]; /*C3H5O */
    YOW += y[29]*imw[29]; /*IXC3H7 */
    YOW += y[30]*imw[30]; /*NXC3H7 */
    YOW += y[31]*imw[31]; /*C3H8 */
    YOW += y[32]*imw[32]; /*IXC3H7O2 */
    YOW += y[33]*imw[33]; /*NXC3H7O2 */
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
    W += c[0]*28.013400; /*N2 */
    W += c[1]*15.999400; /*O */
    W += c[2]*1.007970; /*H */
    W += c[3]*17.007370; /*OH */
    W += c[4]*2.015940; /*H2 */
    W += c[5]*31.998800; /*O2 */
    W += c[6]*18.015340; /*H2O */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*14.027090; /*CH2GSG */
    W += c[10]*30.026490; /*CH2O */
    W += c[11]*15.035060; /*CH3 */
    W += c[12]*28.010550; /*CO */
    W += c[13]*31.034460; /*CH3O */
    W += c[14]*29.062150; /*C2H5 */
    W += c[15]*16.043030; /*CH4 */
    W += c[16]*28.054180; /*C2H4 */
    W += c[17]*30.070120; /*C2H6 */
    W += c[18]*44.009950; /*CO2 */
    W += c[19]*29.018520; /*HCO */
    W += c[20]*47.033860; /*CH3O2 */
    W += c[21]*48.041830; /*CH3O2H */
    W += c[22]*26.038240; /*C2H2 */
    W += c[23]*41.029670; /*HCCO */
    W += c[24]*27.046210; /*C2H3 */
    W += c[25]*43.045610; /*CH2CHO */
    W += c[26]*41.073300; /*C3H5XA */
    W += c[27]*42.081270; /*C3H6 */
    W += c[28]*57.072700; /*C3H5O */
    W += c[29]*43.089240; /*IXC3H7 */
    W += c[30]*43.089240; /*NXC3H7 */
    W += c[31]*44.097210; /*C3H8 */
    W += c[32]*75.088040; /*IXC3H7O2 */
    W += c[33]*75.088040; /*NXC3H7O2 */

    for (id = 0; id < 34; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31446261815324e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKRHOX(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  rho)
{
    amrex::Real XW = 0;/* To hold mean molecular wt */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*2.015940; /*H2 */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*14.027090; /*CH2GSG */
    XW += x[10]*30.026490; /*CH2O */
    XW += x[11]*15.035060; /*CH3 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*29.062150; /*C2H5 */
    XW += x[15]*16.043030; /*CH4 */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*30.070120; /*C2H6 */
    XW += x[18]*44.009950; /*CO2 */
    XW += x[19]*29.018520; /*HCO */
    XW += x[20]*47.033860; /*CH3O2 */
    XW += x[21]*48.041830; /*CH3O2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*41.029670; /*HCCO */
    XW += x[24]*27.046210; /*C2H3 */
    XW += x[25]*43.045610; /*CH2CHO */
    XW += x[26]*41.073300; /*C3H5XA */
    XW += x[27]*42.081270; /*C3H6 */
    XW += x[28]*57.072700; /*C3H5O */
    XW += x[29]*43.089240; /*IXC3H7 */
    XW += x[30]*43.089240; /*NXC3H7 */
    XW += x[31]*44.097210; /*C3H8 */
    XW += x[32]*75.088040; /*IXC3H7O2 */
    XW += x[33]*75.088040; /*NXC3H7O2 */
    *rho = *P * XW / (8.31446261815324e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKRHOY(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  rho)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[34];
    amrex::Real imw[34];

    get_imw(imw);

    for (int i = 0; i < 34; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 34; i++)
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
    W += c[0]*28.013400; /*N2 */
    W += c[1]*15.999400; /*O */
    W += c[2]*1.007970; /*H */
    W += c[3]*17.007370; /*OH */
    W += c[4]*2.015940; /*H2 */
    W += c[5]*31.998800; /*O2 */
    W += c[6]*18.015340; /*H2O */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*14.027090; /*CH2GSG */
    W += c[10]*30.026490; /*CH2O */
    W += c[11]*15.035060; /*CH3 */
    W += c[12]*28.010550; /*CO */
    W += c[13]*31.034460; /*CH3O */
    W += c[14]*29.062150; /*C2H5 */
    W += c[15]*16.043030; /*CH4 */
    W += c[16]*28.054180; /*C2H4 */
    W += c[17]*30.070120; /*C2H6 */
    W += c[18]*44.009950; /*CO2 */
    W += c[19]*29.018520; /*HCO */
    W += c[20]*47.033860; /*CH3O2 */
    W += c[21]*48.041830; /*CH3O2H */
    W += c[22]*26.038240; /*C2H2 */
    W += c[23]*41.029670; /*HCCO */
    W += c[24]*27.046210; /*C2H3 */
    W += c[25]*43.045610; /*CH2CHO */
    W += c[26]*41.073300; /*C3H5XA */
    W += c[27]*42.081270; /*C3H6 */
    W += c[28]*57.072700; /*C3H5O */
    W += c[29]*43.089240; /*IXC3H7 */
    W += c[30]*43.089240; /*NXC3H7 */
    W += c[31]*44.097210; /*C3H8 */
    W += c[32]*75.088040; /*IXC3H7O2 */
    W += c[33]*75.088040; /*NXC3H7O2 */

    for (id = 0; id < 34; ++id) {
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
    amrex::Real tmp[34];
    amrex::Real imw[34];

    get_imw(imw);

    for (int i = 0; i < 34; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 34; i++)
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
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*2.015940; /*H2 */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*14.027090; /*CH2GSG */
    XW += x[10]*30.026490; /*CH2O */
    XW += x[11]*15.035060; /*CH3 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*29.062150; /*C2H5 */
    XW += x[15]*16.043030; /*CH4 */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*30.070120; /*C2H6 */
    XW += x[18]*44.009950; /*CO2 */
    XW += x[19]*29.018520; /*HCO */
    XW += x[20]*47.033860; /*CH3O2 */
    XW += x[21]*48.041830; /*CH3O2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*41.029670; /*HCCO */
    XW += x[24]*27.046210; /*C2H3 */
    XW += x[25]*43.045610; /*CH2CHO */
    XW += x[26]*41.073300; /*C3H5XA */
    XW += x[27]*42.081270; /*C3H6 */
    XW += x[28]*57.072700; /*C3H5O */
    XW += x[29]*43.089240; /*IXC3H7 */
    XW += x[30]*43.089240; /*NXC3H7 */
    XW += x[31]*44.097210; /*C3H8 */
    XW += x[32]*75.088040; /*IXC3H7O2 */
    XW += x[33]*75.088040; /*NXC3H7O2 */
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
    W += c[0]*28.013400; /*N2 */
    W += c[1]*15.999400; /*O */
    W += c[2]*1.007970; /*H */
    W += c[3]*17.007370; /*OH */
    W += c[4]*2.015940; /*H2 */
    W += c[5]*31.998800; /*O2 */
    W += c[6]*18.015340; /*H2O */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*14.027090; /*CH2GSG */
    W += c[10]*30.026490; /*CH2O */
    W += c[11]*15.035060; /*CH3 */
    W += c[12]*28.010550; /*CO */
    W += c[13]*31.034460; /*CH3O */
    W += c[14]*29.062150; /*C2H5 */
    W += c[15]*16.043030; /*CH4 */
    W += c[16]*28.054180; /*C2H4 */
    W += c[17]*30.070120; /*C2H6 */
    W += c[18]*44.009950; /*CO2 */
    W += c[19]*29.018520; /*HCO */
    W += c[20]*47.033860; /*CH3O2 */
    W += c[21]*48.041830; /*CH3O2H */
    W += c[22]*26.038240; /*C2H2 */
    W += c[23]*41.029670; /*HCCO */
    W += c[24]*27.046210; /*C2H3 */
    W += c[25]*43.045610; /*CH2CHO */
    W += c[26]*41.073300; /*C3H5XA */
    W += c[27]*42.081270; /*C3H6 */
    W += c[28]*57.072700; /*C3H5O */
    W += c[29]*43.089240; /*IXC3H7 */
    W += c[30]*43.089240; /*NXC3H7 */
    W += c[31]*44.097210; /*C3H8 */
    W += c[32]*75.088040; /*IXC3H7O2 */
    W += c[33]*75.088040; /*NXC3H7O2 */

    for (id = 0; id < 34; ++id) {
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
    amrex::Real tmp[34];
    amrex::Real imw[34];

    get_imw(imw);

    for (int i = 0; i < 34; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 34; i++)
    {
        YOW += tmp[i];
    }

    amrex::Real YOWINV = 1.0/YOW;

    for (int i = 0; i < 34; i++)
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
    amrex::Real imw[34];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 34; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 34; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446261815324e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 34; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKYTCR(amrex::Real *  rho, amrex::Real * /* T */, amrex::Real *  y,  amrex::Real *  c)
{
    amrex::Real imw[34];

    get_imw(imw);

    for (int i = 0; i < 34; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKXTY(amrex::Real *  x,  amrex::Real *  y)
{
    amrex::Real XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*2.015940; /*H2 */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*14.027090; /*CH2GSG */
    XW += x[10]*30.026490; /*CH2O */
    XW += x[11]*15.035060; /*CH3 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*29.062150; /*C2H5 */
    XW += x[15]*16.043030; /*CH4 */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*30.070120; /*C2H6 */
    XW += x[18]*44.009950; /*CO2 */
    XW += x[19]*29.018520; /*HCO */
    XW += x[20]*47.033860; /*CH3O2 */
    XW += x[21]*48.041830; /*CH3O2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*41.029670; /*HCCO */
    XW += x[24]*27.046210; /*C2H3 */
    XW += x[25]*43.045610; /*CH2CHO */
    XW += x[26]*41.073300; /*C3H5XA */
    XW += x[27]*42.081270; /*C3H6 */
    XW += x[28]*57.072700; /*C3H5O */
    XW += x[29]*43.089240; /*IXC3H7 */
    XW += x[30]*43.089240; /*NXC3H7 */
    XW += x[31]*44.097210; /*C3H8 */
    XW += x[32]*75.088040; /*IXC3H7O2 */
    XW += x[33]*75.088040; /*NXC3H7O2 */
    /*Now compute conversion */
    amrex::Real XWinv = 1.0/XW;
    y[0] = x[0]*28.013400*XWinv; 
    y[1] = x[1]*15.999400*XWinv; 
    y[2] = x[2]*1.007970*XWinv; 
    y[3] = x[3]*17.007370*XWinv; 
    y[4] = x[4]*2.015940*XWinv; 
    y[5] = x[5]*31.998800*XWinv; 
    y[6] = x[6]*18.015340*XWinv; 
    y[7] = x[7]*34.014740*XWinv; 
    y[8] = x[8]*33.006770*XWinv; 
    y[9] = x[9]*14.027090*XWinv; 
    y[10] = x[10]*30.026490*XWinv; 
    y[11] = x[11]*15.035060*XWinv; 
    y[12] = x[12]*28.010550*XWinv; 
    y[13] = x[13]*31.034460*XWinv; 
    y[14] = x[14]*29.062150*XWinv; 
    y[15] = x[15]*16.043030*XWinv; 
    y[16] = x[16]*28.054180*XWinv; 
    y[17] = x[17]*30.070120*XWinv; 
    y[18] = x[18]*44.009950*XWinv; 
    y[19] = x[19]*29.018520*XWinv; 
    y[20] = x[20]*47.033860*XWinv; 
    y[21] = x[21]*48.041830*XWinv; 
    y[22] = x[22]*26.038240*XWinv; 
    y[23] = x[23]*41.029670*XWinv; 
    y[24] = x[24]*27.046210*XWinv; 
    y[25] = x[25]*43.045610*XWinv; 
    y[26] = x[26]*41.073300*XWinv; 
    y[27] = x[27]*42.081270*XWinv; 
    y[28] = x[28]*57.072700*XWinv; 
    y[29] = x[29]*43.089240*XWinv; 
    y[30] = x[30]*43.089240*XWinv; 
    y[31] = x[31]*44.097210*XWinv; 
    y[32] = x[32]*75.088040*XWinv; 
    y[33] = x[33]*75.088040*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKXTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  c)
{
    int id; /*loop counter */
    amrex::Real PORT = (*P)/(8.31446261815324e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 34; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKXTCR(amrex::Real *  rho, amrex::Real * /* T */, amrex::Real *  x, amrex::Real *  c)
{
    int id; /*loop counter */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*2.015940; /*H2 */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*14.027090; /*CH2GSG */
    XW += x[10]*30.026490; /*CH2O */
    XW += x[11]*15.035060; /*CH3 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*29.062150; /*C2H5 */
    XW += x[15]*16.043030; /*CH4 */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*30.070120; /*C2H6 */
    XW += x[18]*44.009950; /*CO2 */
    XW += x[19]*29.018520; /*HCO */
    XW += x[20]*47.033860; /*CH3O2 */
    XW += x[21]*48.041830; /*CH3O2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*41.029670; /*HCCO */
    XW += x[24]*27.046210; /*C2H3 */
    XW += x[25]*43.045610; /*CH2CHO */
    XW += x[26]*41.073300; /*C3H5XA */
    XW += x[27]*42.081270; /*C3H6 */
    XW += x[28]*57.072700; /*C3H5O */
    XW += x[29]*43.089240; /*IXC3H7 */
    XW += x[30]*43.089240; /*NXC3H7 */
    XW += x[31]*44.097210; /*C3H8 */
    XW += x[32]*75.088040; /*IXC3H7O2 */
    XW += x[33]*75.088040; /*NXC3H7O2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 34; ++id) {
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
    for (id = 0; id < 34; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    amrex::Real sumCinv = 1.0/sumC;
    for (id = 0; id < 34; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCTY(amrex::Real *  c, amrex::Real *  y)
{
    amrex::Real CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*28.013400; /*N2 */
    CW += c[1]*15.999400; /*O */
    CW += c[2]*1.007970; /*H */
    CW += c[3]*17.007370; /*OH */
    CW += c[4]*2.015940; /*H2 */
    CW += c[5]*31.998800; /*O2 */
    CW += c[6]*18.015340; /*H2O */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*33.006770; /*HO2 */
    CW += c[9]*14.027090; /*CH2GSG */
    CW += c[10]*30.026490; /*CH2O */
    CW += c[11]*15.035060; /*CH3 */
    CW += c[12]*28.010550; /*CO */
    CW += c[13]*31.034460; /*CH3O */
    CW += c[14]*29.062150; /*C2H5 */
    CW += c[15]*16.043030; /*CH4 */
    CW += c[16]*28.054180; /*C2H4 */
    CW += c[17]*30.070120; /*C2H6 */
    CW += c[18]*44.009950; /*CO2 */
    CW += c[19]*29.018520; /*HCO */
    CW += c[20]*47.033860; /*CH3O2 */
    CW += c[21]*48.041830; /*CH3O2H */
    CW += c[22]*26.038240; /*C2H2 */
    CW += c[23]*41.029670; /*HCCO */
    CW += c[24]*27.046210; /*C2H3 */
    CW += c[25]*43.045610; /*CH2CHO */
    CW += c[26]*41.073300; /*C3H5XA */
    CW += c[27]*42.081270; /*C3H6 */
    CW += c[28]*57.072700; /*C3H5O */
    CW += c[29]*43.089240; /*IXC3H7 */
    CW += c[30]*43.089240; /*NXC3H7 */
    CW += c[31]*44.097210; /*C3H8 */
    CW += c[32]*75.088040; /*IXC3H7O2 */
    CW += c[33]*75.088040; /*NXC3H7O2 */
    /*Now compute conversion */
    amrex::Real CWinv = 1.0/CW;
    y[0] = c[0]*28.013400*CWinv; 
    y[1] = c[1]*15.999400*CWinv; 
    y[2] = c[2]*1.007970*CWinv; 
    y[3] = c[3]*17.007370*CWinv; 
    y[4] = c[4]*2.015940*CWinv; 
    y[5] = c[5]*31.998800*CWinv; 
    y[6] = c[6]*18.015340*CWinv; 
    y[7] = c[7]*34.014740*CWinv; 
    y[8] = c[8]*33.006770*CWinv; 
    y[9] = c[9]*14.027090*CWinv; 
    y[10] = c[10]*30.026490*CWinv; 
    y[11] = c[11]*15.035060*CWinv; 
    y[12] = c[12]*28.010550*CWinv; 
    y[13] = c[13]*31.034460*CWinv; 
    y[14] = c[14]*29.062150*CWinv; 
    y[15] = c[15]*16.043030*CWinv; 
    y[16] = c[16]*28.054180*CWinv; 
    y[17] = c[17]*30.070120*CWinv; 
    y[18] = c[18]*44.009950*CWinv; 
    y[19] = c[19]*29.018520*CWinv; 
    y[20] = c[20]*47.033860*CWinv; 
    y[21] = c[21]*48.041830*CWinv; 
    y[22] = c[22]*26.038240*CWinv; 
    y[23] = c[23]*41.029670*CWinv; 
    y[24] = c[24]*27.046210*CWinv; 
    y[25] = c[25]*43.045610*CWinv; 
    y[26] = c[26]*41.073300*CWinv; 
    y[27] = c[27]*42.081270*CWinv; 
    y[28] = c[28]*57.072700*CWinv; 
    y[29] = c[29]*43.089240*CWinv; 
    y[30] = c[30]*43.089240*CWinv; 
    y[31] = c[31]*44.097210*CWinv; 
    y[32] = c[32]*75.088040*CWinv; 
    y[33] = c[33]*75.088040*CWinv; 

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
    for (id = 0; id < 34; ++id) {
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
    for (id = 0; id < 34; ++id) {
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
    for (id = 0; id < 34; ++id) {
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
    for (id = 0; id < 34; ++id) {
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
    for (id = 0; id < 34; ++id) {
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
    cvms[0] *= 2.968030520448514e+06; /*N2 */
    cvms[1] *= 5.196734013871295e+06; /*O */
    cvms[2] *= 8.248720317224957e+07; /*H */
    cvms[3] *= 4.888740950630956e+06; /*OH */
    cvms[4] *= 4.124360158612479e+07; /*H2 */
    cvms[5] *= 2.598367006935648e+06; /*O2 */
    cvms[6] *= 4.615212712140454e+06; /*H2O */
    cvms[7] *= 2.444370475315478e+06; /*H2O2 */
    cvms[8] *= 2.519017346487778e+06; /*HO2 */
    cvms[9] *= 5.927432288630956e+06; /*CH2GSG */
    cvms[10] *= 2.769042474879095e+06; /*CH2O */
    cvms[11] *= 5.530049509714786e+06; /*CH3 */
    cvms[12] *= 2.968332509769797e+06; /*CO */
    cvms[13] *= 2.679106586083096e+06; /*CH3O */
    cvms[14] *= 2.860924817383862e+06; /*C2H5 */
    cvms[15] *= 5.182601178301878e+06; /*CH4 */
    cvms[16] *= 2.963716144315478e+06; /*C2H4 */
    cvms[17] *= 2.765024754857393e+06; /*C2H6 */
    cvms[18] *= 1.889223372931176e+06; /*CO2 */
    cvms[19] *= 2.865226282440744e+06; /*HCO */
    cvms[20] *= 1.767761059405552e+06; /*CH3O2 */
    cvms[21] *= 1.730671504010825e+06; /*CH3O2H */
    cvms[22] *= 3.193173815954242e+06; /*C2H2 */
    cvms[23] *= 2.026451252996488e+06; /*HCCO */
    cvms[24] *= 3.074169215632519e+06; /*C2H3 */
    cvms[25] *= 1.931547170118681e+06; /*CH2CHO */
    cvms[26] *= 2.024298660724422e+06; /*C3H5XA */
    cvms[27] *= 1.975810762876985e+06; /*C3H6 */
    cvms[28] *= 1.456819568401922e+06; /*C3H5O */
    cvms[29] *= 1.929591382478141e+06; /*IXC3H7 */
    cvms[30] *= 1.929591382478141e+06; /*NXC3H7 */
    cvms[31] *= 1.885484958833731e+06; /*C3H8 */
    cvms[32] *= 1.107295198829699e+06; /*IXC3H7O2 */
    cvms[33] *= 1.107295198829699e+06; /*NXC3H7O2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKCPMS(amrex::Real *  T,  amrex::Real *  cpms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 2.968030520448514e+06; /*N2 */
    cpms[1] *= 5.196734013871295e+06; /*O */
    cpms[2] *= 8.248720317224957e+07; /*H */
    cpms[3] *= 4.888740950630956e+06; /*OH */
    cpms[4] *= 4.124360158612479e+07; /*H2 */
    cpms[5] *= 2.598367006935648e+06; /*O2 */
    cpms[6] *= 4.615212712140454e+06; /*H2O */
    cpms[7] *= 2.444370475315478e+06; /*H2O2 */
    cpms[8] *= 2.519017346487778e+06; /*HO2 */
    cpms[9] *= 5.927432288630956e+06; /*CH2GSG */
    cpms[10] *= 2.769042474879095e+06; /*CH2O */
    cpms[11] *= 5.530049509714786e+06; /*CH3 */
    cpms[12] *= 2.968332509769797e+06; /*CO */
    cpms[13] *= 2.679106586083096e+06; /*CH3O */
    cpms[14] *= 2.860924817383862e+06; /*C2H5 */
    cpms[15] *= 5.182601178301878e+06; /*CH4 */
    cpms[16] *= 2.963716144315478e+06; /*C2H4 */
    cpms[17] *= 2.765024754857393e+06; /*C2H6 */
    cpms[18] *= 1.889223372931176e+06; /*CO2 */
    cpms[19] *= 2.865226282440744e+06; /*HCO */
    cpms[20] *= 1.767761059405552e+06; /*CH3O2 */
    cpms[21] *= 1.730671504010825e+06; /*CH3O2H */
    cpms[22] *= 3.193173815954242e+06; /*C2H2 */
    cpms[23] *= 2.026451252996488e+06; /*HCCO */
    cpms[24] *= 3.074169215632519e+06; /*C2H3 */
    cpms[25] *= 1.931547170118681e+06; /*CH2CHO */
    cpms[26] *= 2.024298660724422e+06; /*C3H5XA */
    cpms[27] *= 1.975810762876985e+06; /*C3H6 */
    cpms[28] *= 1.456819568401922e+06; /*C3H5O */
    cpms[29] *= 1.929591382478141e+06; /*IXC3H7 */
    cpms[30] *= 1.929591382478141e+06; /*NXC3H7 */
    cpms[31] *= 1.885484958833731e+06; /*C3H8 */
    cpms[32] *= 1.107295198829699e+06; /*IXC3H7O2 */
    cpms[33] *= 1.107295198829699e+06; /*NXC3H7O2 */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKUMS(amrex::Real *  T,  amrex::Real *  ums)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[34];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    for (int i = 0; i < 34; i++)
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
    amrex::Real imw[34];

    get_imw(imw);

    speciesEnthalpy(hms, tc);

    for (int i = 0; i < 34; i++)
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
    sms[0] *= 2.968030520448514e+06; /*N2 */
    sms[1] *= 5.196734013871295e+06; /*O */
    sms[2] *= 8.248720317224957e+07; /*H */
    sms[3] *= 4.888740950630956e+06; /*OH */
    sms[4] *= 4.124360158612479e+07; /*H2 */
    sms[5] *= 2.598367006935648e+06; /*O2 */
    sms[6] *= 4.615212712140454e+06; /*H2O */
    sms[7] *= 2.444370475315478e+06; /*H2O2 */
    sms[8] *= 2.519017346487778e+06; /*HO2 */
    sms[9] *= 5.927432288630956e+06; /*CH2GSG */
    sms[10] *= 2.769042474879095e+06; /*CH2O */
    sms[11] *= 5.530049509714786e+06; /*CH3 */
    sms[12] *= 2.968332509769797e+06; /*CO */
    sms[13] *= 2.679106586083096e+06; /*CH3O */
    sms[14] *= 2.860924817383862e+06; /*C2H5 */
    sms[15] *= 5.182601178301878e+06; /*CH4 */
    sms[16] *= 2.963716144315478e+06; /*C2H4 */
    sms[17] *= 2.765024754857393e+06; /*C2H6 */
    sms[18] *= 1.889223372931176e+06; /*CO2 */
    sms[19] *= 2.865226282440744e+06; /*HCO */
    sms[20] *= 1.767761059405552e+06; /*CH3O2 */
    sms[21] *= 1.730671504010825e+06; /*CH3O2H */
    sms[22] *= 3.193173815954242e+06; /*C2H2 */
    sms[23] *= 2.026451252996488e+06; /*HCCO */
    sms[24] *= 3.074169215632519e+06; /*C2H3 */
    sms[25] *= 1.931547170118681e+06; /*CH2CHO */
    sms[26] *= 2.024298660724422e+06; /*C3H5XA */
    sms[27] *= 1.975810762876985e+06; /*C3H6 */
    sms[28] *= 1.456819568401922e+06; /*C3H5O */
    sms[29] *= 1.929591382478141e+06; /*IXC3H7 */
    sms[30] *= 1.929591382478141e+06; /*NXC3H7 */
    sms[31] *= 1.885484958833731e+06; /*C3H8 */
    sms[32] *= 1.107295198829699e+06; /*IXC3H7O2 */
    sms[33] *= 1.107295198829699e+06; /*NXC3H7O2 */
}

#ifdef AMREX_USE_GPU

/* GPU version of productionRate: no more use of thermo namespace vectors  */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void comp_qfqr(amrex::Real *  qf, amrex::Real * qr, amrex::Real * sc, amrex::Real * qss_sc, amrex::Real * tc, amrex::Real invT)
{

    /*reaction 7: 2.000000 OH (+M) <=> H2O2 (+M) */
    qf[0] = pow(sc[3], 2.000000);
    qr[0] = sc[7];

    /*reaction 13: H + O2 (+M) <=> HO2 (+M) */
    qf[1] = sc[2]*sc[5];
    qr[1] = sc[8];

    /*reaction 46: CH3 + H (+M) => CH4 (+M) */
    qf[2] = sc[2]*sc[11];
    qr[2] = 0.0;

    /*reaction 47: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    qf[3] = pow(sc[11], 2.000000);
    qr[3] = sc[17];

    /*reaction 61: CO + O (+M) => CO2 (+M) */
    qf[4] = sc[1]*sc[12];
    qr[4] = 0.0;

    /*reaction 76: CH3O (+M) => CH2O + H (+M) */
    qf[5] = sc[13];
    qr[5] = 0.0;

    /*reaction 93: C2H3 (+M) => H + C2H2 (+M) */
    qf[6] = sc[24];
    qr[6] = 0.0;

    /*reaction 95: H + C2H4 (+M) <=> C2H5 (+M) */
    qf[7] = sc[2]*sc[16];
    qr[7] = sc[14];

    /*reaction 104: C2H4 (+M) => C2H2 + H2 (+M) */
    qf[8] = sc[16];
    qr[8] = 0.0;

    /*reaction 1: O + H + M => OH + M */
    qf[9] = sc[1]*sc[2];
    qr[9] = 0.0;

    /*reaction 4: 2.000000 O + M => O2 + M */
    qf[10] = pow(sc[1], 2.000000);
    qr[10] = 0.0;

    /*reaction 8: H + OH + M => H2O + M */
    qf[11] = sc[2]*sc[3];
    qr[11] = 0.0;

    /*reaction 67: HCO + M => H + CO + M */
    qf[12] = sc[19];
    qr[12] = 0.0;

    /*reaction 68: H + CO + M => HCO + M */
    qf[13] = sc[2]*sc[12];
    qr[13] = 0.0;

    /*reaction 83: CH3O2 + M => CH3 + O2 + M */
    qf[14] = sc[20];
    qr[14] = 0.0;

    /*reaction 84: CH3 + O2 + M => CH3O2 + M */
    qf[15] = sc[5]*sc[11];
    qr[15] = 0.0;

    /*reaction 2: O + H2 => H + OH */
    qf[16] = sc[1]*sc[4];
    qr[16] = 0.0;

    /*reaction 3: H + OH => O + H2 */
    qf[17] = sc[2]*sc[3];
    qr[17] = 0.0;

    /*reaction 5: OH + H2 => H + H2O */
    qf[18] = sc[3]*sc[4];
    qr[18] = 0.0;

    /*reaction 6: H + H2O => OH + H2 */
    qf[19] = sc[2]*sc[6];
    qr[19] = 0.0;

    /*reaction 9: O + H2O => 2.000000 OH */
    qf[20] = sc[1]*sc[6];
    qr[20] = 0.0;

    /*reaction 10: 2.000000 OH => O + H2O */
    qf[21] = pow(sc[3], 2.000000);
    qr[21] = 0.0;

    /*reaction 11: H + O2 => O + OH */
    qf[22] = sc[2]*sc[5];
    qr[22] = 0.0;

    /*reaction 12: O + OH => H + O2 */
    qf[23] = sc[1]*sc[3];
    qr[23] = 0.0;

    /*reaction 14: HO2 + H => 2.000000 OH */
    qf[24] = sc[2]*sc[8];
    qr[24] = 0.0;

    /*reaction 15: 2.000000 HO2 => H2O2 + O2 */
    qf[25] = pow(sc[8], 2.000000);
    qr[25] = 0.0;

    /*reaction 16: HO2 + H => H2 + O2 */
    qf[26] = sc[2]*sc[8];
    qr[26] = 0.0;

    /*reaction 17: HO2 + OH => H2O + O2 */
    qf[27] = sc[3]*sc[8];
    qr[27] = 0.0;

    /*reaction 18: H2O + O2 => HO2 + OH */
    qf[28] = sc[5]*sc[6];
    qr[28] = 0.0;

    /*reaction 19: 2.000000 HO2 => H2O2 + O2 */
    qf[29] = pow(sc[8], 2.000000);
    qr[29] = 0.0;

    /*reaction 20: HO2 + O => OH + O2 */
    qf[30] = sc[1]*sc[8];
    qr[30] = 0.0;

    /*reaction 21: OH + O2 => HO2 + O */
    qf[31] = sc[3]*sc[5];
    qr[31] = 0.0;

    /*reaction 22: H2O2 + OH => H2O + HO2 */
    qf[32] = sc[3]*sc[7];
    qr[32] = 0.0;

    /*reaction 23: H2O2 + H => H2O + OH */
    qf[33] = sc[2]*sc[7];
    qr[33] = 0.0;

    /*reaction 24: H2O2 + OH => H2O + HO2 */
    qf[34] = sc[3]*sc[7];
    qr[34] = 0.0;

    /*reaction 25: H2O + HO2 => H2O2 + OH */
    qf[35] = sc[6]*sc[8];
    qr[35] = 0.0;

    /*reaction 26: H2O2 + O => OH + HO2 */
    qf[36] = sc[1]*sc[7];
    qr[36] = 0.0;

    /*reaction 27: H2O2 + H => H2 + HO2 */
    qf[37] = sc[2]*sc[7];
    qr[37] = 0.0;

    /*reaction 28: H2 + HO2 => H2O2 + H */
    qf[38] = sc[4]*sc[8];
    qr[38] = 0.0;

    /*reaction 29: CH2GSG + OH => CH2O + H */
    qf[39] = sc[3]*sc[9];
    qr[39] = 0.0;

    /*reaction 30: CH2GSG + H2 => CH3 + H */
    qf[40] = sc[4]*sc[9];
    qr[40] = 0.0;

    /*reaction 31: CH3 + H => CH2GSG + H2 */
    qf[41] = sc[2]*sc[11];
    qr[41] = 0.0;

    /*reaction 32: CH2GSG + O2 => CO + OH + H */
    qf[42] = sc[5]*sc[9];
    qr[42] = 0.0;

    /*reaction 33: CH2GSG + O => CO + 2.000000 H */
    qf[43] = sc[1]*sc[9];
    qr[43] = 0.0;

    /*reaction 34: CH3 + HO2 => CH3O + OH */
    qf[44] = sc[8]*sc[11];
    qr[44] = 0.0;

    /*reaction 35: CH3 + O2 => CH2O + OH */
    qf[45] = sc[5]*sc[11];
    qr[45] = 0.0;

    /*reaction 36: CH3 + O2 => CH3O + O */
    qf[46] = sc[5]*sc[11];
    qr[46] = 0.0;

    /*reaction 37: CH3O + O => CH3 + O2 */
    qf[47] = sc[1]*sc[13];
    qr[47] = 0.0;

    /*reaction 38: 2.000000 CH3 <=> H + C2H5 */
    qf[48] = pow(sc[11], 2.000000);
    qr[48] = sc[2]*sc[14];

    /*reaction 39: CH3 + HO2 => CH4 + O2 */
    qf[49] = sc[8]*sc[11];
    qr[49] = 0.0;

    /*reaction 40: CH3 + O => CH2O + H */
    qf[50] = sc[1]*sc[11];
    qr[50] = 0.0;

    /*reaction 41: CH3 + OH => CH2O + H2 */
    qf[51] = sc[3]*sc[11];
    qr[51] = 0.0;

    /*reaction 42: CH3 + OH => CH2GSG + H2O */
    qf[52] = sc[3]*sc[11];
    qr[52] = 0.0;

    /*reaction 43: CH2GSG + H2O => CH3 + OH */
    qf[53] = sc[6]*sc[9];
    qr[53] = 0.0;

    /*reaction 44: CH3 + H2O2 => CH4 + HO2 */
    qf[54] = sc[7]*sc[11];
    qr[54] = 0.0;

    /*reaction 45: CH2GSG + CH3 => C2H4 + H */
    qf[55] = sc[9]*sc[11];
    qr[55] = 0.0;

    /*reaction 48: CH2GSG + CH4 => 2.000000 CH3 */
    qf[56] = sc[9]*sc[15];
    qr[56] = 0.0;

    /*reaction 49: 2.000000 CH3 => CH2GSG + CH4 */
    qf[57] = pow(sc[11], 2.000000);
    qr[57] = 0.0;

    /*reaction 50: CH4 + O => CH3 + OH */
    qf[58] = sc[1]*sc[15];
    qr[58] = 0.0;

    /*reaction 51: CH3 + OH => CH4 + O */
    qf[59] = sc[3]*sc[11];
    qr[59] = 0.0;

    /*reaction 52: CH4 + H => CH3 + H2 */
    qf[60] = sc[2]*sc[15];
    qr[60] = 0.0;

    /*reaction 53: CH3 + H2 => CH4 + H */
    qf[61] = sc[4]*sc[11];
    qr[61] = 0.0;

    /*reaction 54: CH4 + OH => CH3 + H2O */
    qf[62] = sc[3]*sc[15];
    qr[62] = 0.0;

    /*reaction 55: CH3 + H2O => CH4 + OH */
    qf[63] = sc[6]*sc[11];
    qr[63] = 0.0;

    /*reaction 56: CO + HO2 => CO2 + OH */
    qf[64] = sc[8]*sc[12];
    qr[64] = 0.0;

    /*reaction 57: CO + O2 => CO2 + O */
    qf[65] = sc[5]*sc[12];
    qr[65] = 0.0;

    /*reaction 58: CO2 + O => CO + O2 */
    qf[66] = sc[1]*sc[18];
    qr[66] = 0.0;

    /*reaction 59: CO + OH => CO2 + H */
    qf[67] = sc[3]*sc[12];
    qr[67] = 0.0;

    /*reaction 60: CO2 + H => CO + OH */
    qf[68] = sc[2]*sc[18];
    qr[68] = 0.0;

    /*reaction 62: HCO + CH3 => CH4 + CO */
    qf[69] = sc[11]*sc[19];
    qr[69] = 0.0;

    /*reaction 63: HCO + H => CO + H2 */
    qf[70] = sc[2]*sc[19];
    qr[70] = 0.0;

    /*reaction 64: HCO + O2 => CO + HO2 */
    qf[71] = sc[5]*sc[19];
    qr[71] = 0.0;

    /*reaction 65: HCO + O => CO + OH */
    qf[72] = sc[1]*sc[19];
    qr[72] = 0.0;

    /*reaction 66: HCO + O => CO2 + H */
    qf[73] = sc[1]*sc[19];
    qr[73] = 0.0;

    /*reaction 69: HCO + OH => CO + H2O */
    qf[74] = sc[3]*sc[19];
    qr[74] = 0.0;

    /*reaction 70: CH2O + H => HCO + H2 */
    qf[75] = sc[2]*sc[10];
    qr[75] = 0.0;

    /*reaction 71: CH2O + O2 => HCO + HO2 */
    qf[76] = sc[5]*sc[10];
    qr[76] = 0.0;

    /*reaction 72: CH2O + OH => HCO + H2O */
    qf[77] = sc[3]*sc[10];
    qr[77] = 0.0;

    /*reaction 73: CH2O + HO2 => HCO + H2O2 */
    qf[78] = sc[8]*sc[10];
    qr[78] = 0.0;

    /*reaction 74: CH2O + O => HCO + OH */
    qf[79] = sc[1]*sc[10];
    qr[79] = 0.0;

    /*reaction 75: CH2O + CH3 => HCO + CH4 */
    qf[80] = sc[10]*sc[11];
    qr[80] = 0.0;

    /*reaction 77: CH3O + O2 => CH2O + HO2 */
    qf[81] = sc[5]*sc[13];
    qr[81] = 0.0;

    /*reaction 78: CH2GSG + CO2 => CH2O + CO */
    qf[82] = sc[9]*sc[18];
    qr[82] = 0.0;

    /*reaction 79: CH3O2 + CH3 => 2.000000 CH3O */
    qf[83] = sc[11]*sc[20];
    qr[83] = 0.0;

    /*reaction 80: CH3O2 + HO2 => CH3O2H + O2 */
    qf[84] = sc[8]*sc[20];
    qr[84] = 0.0;

    /*reaction 81: 2.000000 CH3O2 => O2 + 2.000000 CH3O */
    qf[85] = pow(sc[20], 2.000000);
    qr[85] = 0.0;

    /*reaction 82: CH3O2 + CH2O => CH3O2H + HCO */
    qf[86] = sc[10]*sc[20];
    qr[86] = 0.0;

    /*reaction 85: CH3O2H => CH3O + OH */
    qf[87] = sc[21];
    qr[87] = 0.0;

    /*reaction 86: CH3O + OH => CH3O2H */
    qf[88] = sc[3]*sc[13];
    qr[88] = 0.0;

    /*reaction 87: C2H2 + O2 => HCCO + OH */
    qf[89] = sc[5]*sc[22];
    qr[89] = 0.0;

    /*reaction 88: C2H2 + O => HCCO + H */
    qf[90] = sc[1]*sc[22];
    qr[90] = 0.0;

    /*reaction 89: C2H3 + O2 => CH2CHO + O */
    qf[91] = sc[5]*sc[24];
    qr[91] = 0.0;

    /*reaction 90: C2H3 + H => C2H2 + H2 */
    qf[92] = sc[2]*sc[24];
    qr[92] = 0.0;

    /*reaction 91: C2H3 + O2 => CH2O + HCO */
    qf[93] = sc[5]*sc[24];
    qr[93] = 0.0;

    /*reaction 92: C2H3 + O2 => C2H2 + HO2 */
    qf[94] = sc[5]*sc[24];
    qr[94] = 0.0;

    /*reaction 94: C2H4 + O => CH3 + HCO */
    qf[95] = sc[1]*sc[16];
    qr[95] = 0.0;

    /*reaction 96: C2H4 + OH => C2H3 + H2O */
    qf[96] = sc[3]*sc[16];
    qr[96] = 0.0;

    /*reaction 97: C2H3 + H2O => C2H4 + OH */
    qf[97] = sc[6]*sc[24];
    qr[97] = 0.0;

    /*reaction 98: C2H4 + H => C2H3 + H2 */
    qf[98] = sc[2]*sc[16];
    qr[98] = 0.0;

    /*reaction 99: C2H3 + H2 => C2H4 + H */
    qf[99] = sc[4]*sc[24];
    qr[99] = 0.0;

    /*reaction 100: C2H4 + CH3O2 => C2H3 + CH3O2H */
    qf[100] = sc[16]*sc[20];
    qr[100] = 0.0;

    /*reaction 101: C2H4 + CH3 => C2H3 + CH4 */
    qf[101] = sc[11]*sc[16];
    qr[101] = 0.0;

    /*reaction 102: C2H3 + CH4 => C2H4 + CH3 */
    qf[102] = sc[15]*sc[24];
    qr[102] = 0.0;

    /*reaction 103: C2H4 + O => CH2CHO + H */
    qf[103] = sc[1]*sc[16];
    qr[103] = 0.0;

    /*reaction 105: CH3 + C2H5 => CH4 + C2H4 */
    qf[104] = sc[11]*sc[14];
    qr[104] = 0.0;

    /*reaction 106: C2H5 + O2 => C2H4 + HO2 */
    qf[105] = sc[5]*sc[14];
    qr[105] = 0.0;

    /*reaction 107: C2H4 + HO2 => C2H5 + O2 */
    qf[106] = sc[8]*sc[16];
    qr[106] = 0.0;

    /*reaction 108: C2H5 + HO2 => C2H6 + O2 */
    qf[107] = sc[8]*sc[14];
    qr[107] = 0.0;

    /*reaction 109: H + C2H5 => C2H6 */
    qf[108] = sc[2]*sc[14];
    qr[108] = 0.0;

    /*reaction 110: C2H6 + H => C2H5 + H2 */
    qf[109] = sc[2]*sc[17];
    qr[109] = 0.0;

    /*reaction 111: C2H5 + H2 => C2H6 + H */
    qf[110] = sc[4]*sc[14];
    qr[110] = 0.0;

    /*reaction 112: C2H6 + OH => C2H5 + H2O */
    qf[111] = sc[3]*sc[17];
    qr[111] = 0.0;

    /*reaction 113: CH2GSG + C2H6 => CH3 + C2H5 */
    qf[112] = sc[9]*sc[17];
    qr[112] = 0.0;

    /*reaction 114: C2H6 + O => C2H5 + OH */
    qf[113] = sc[1]*sc[17];
    qr[113] = 0.0;

    /*reaction 115: C2H6 + CH3 => C2H5 + CH4 */
    qf[114] = sc[11]*sc[17];
    qr[114] = 0.0;

    /*reaction 116: HCCO + O => H + 2.000000 CO */
    qf[115] = sc[1]*sc[23];
    qr[115] = 0.0;

    /*reaction 117: HCCO + O2 => CO2 + HCO */
    qf[116] = sc[5]*sc[23];
    qr[116] = 0.0;

    /*reaction 118: HCCO + OH => 2.000000 HCO */
    qf[117] = sc[3]*sc[23];
    qr[117] = 0.0;

    /*reaction 119: HCCO + H => CH2GSG + CO */
    qf[118] = sc[2]*sc[23];
    qr[118] = 0.0;

    /*reaction 120: CH2GSG + CO => HCCO + H */
    qf[119] = sc[9]*sc[12];
    qr[119] = 0.0;

    /*reaction 121: CH2CHO + O2 => CH2O + CO + OH */
    qf[120] = sc[5]*sc[25];
    qr[120] = 0.0;

    /*reaction 122: C3H5XA + CH2O => C3H6 + HCO */
    qf[121] = sc[10]*sc[26];
    qr[121] = 0.0;

    /*reaction 123: C3H6 + HCO => C3H5XA + CH2O */
    qf[122] = sc[19]*sc[27];
    qr[122] = 0.0;

    /*reaction 124: C3H5XA + HO2 => C3H5O + OH */
    qf[123] = sc[8]*sc[26];
    qr[123] = 0.0;

    /*reaction 125: C3H5XA + CH3O2 => C3H5O + CH3O */
    qf[124] = sc[20]*sc[26];
    qr[124] = 0.0;

    /*reaction 126: C3H5XA + HO2 => C3H6 + O2 */
    qf[125] = sc[8]*sc[26];
    qr[125] = 0.0;

    /*reaction 127: C3H5XA + O2 => CH2CHO + CH2O */
    qf[126] = sc[5]*sc[26];
    qr[126] = 0.0;

    /*reaction 128: C3H5XA + O2 => C2H2 + CH2O + OH */
    qf[127] = sc[5]*sc[26];
    qr[127] = 0.0;

    /*reaction 129: C3H5XA + H => C3H6 */
    qf[128] = sc[2]*sc[26];
    qr[128] = 0.0;

    /*reaction 130: C3H5XA => C2H2 + CH3 */
    qf[129] = sc[26];
    qr[129] = 0.0;

    /*reaction 131: C3H6 => C2H3 + CH3 */
    qf[130] = sc[27];
    qr[130] = 0.0;

    /*reaction 132: C2H3 + CH3 => C3H6 */
    qf[131] = sc[11]*sc[24];
    qr[131] = 0.0;

    /*reaction 133: C3H6 + O => C2H5 + HCO */
    qf[132] = sc[1]*sc[27];
    qr[132] = 0.0;

    /*reaction 134: C3H6 + H => C2H4 + CH3 */
    qf[133] = sc[2]*sc[27];
    qr[133] = 0.0;

    /*reaction 135: C2H4 + CH3 => C3H6 + H */
    qf[134] = sc[11]*sc[16];
    qr[134] = 0.0;

    /*reaction 136: C3H6 + H => C3H5XA + H2 */
    qf[135] = sc[2]*sc[27];
    qr[135] = 0.0;

    /*reaction 137: C3H5XA + H2 => C3H6 + H */
    qf[136] = sc[4]*sc[26];
    qr[136] = 0.0;

    /*reaction 138: C3H6 + OH => C3H5XA + H2O */
    qf[137] = sc[3]*sc[27];
    qr[137] = 0.0;

    /*reaction 139: C3H6 + O => C3H5XA + OH */
    qf[138] = sc[1]*sc[27];
    qr[138] = 0.0;

    /*reaction 140: C3H6 + CH3 => C3H5XA + CH4 */
    qf[139] = sc[11]*sc[27];
    qr[139] = 0.0;

    /*reaction 141: IXC3H7 + O2 => C3H6 + HO2 */
    qf[140] = sc[5]*sc[29];
    qr[140] = 0.0;

    /*reaction 142: IXC3H7 => H + C3H6 */
    qf[141] = sc[29];
    qr[141] = 0.0;

    /*reaction 143: H + C3H6 => IXC3H7 */
    qf[142] = sc[2]*sc[27];
    qr[142] = 0.0;

    /*reaction 144: NXC3H7 => CH3 + C2H4 */
    qf[143] = sc[30];
    qr[143] = 0.0;

    /*reaction 145: CH3 + C2H4 => NXC3H7 */
    qf[144] = sc[11]*sc[16];
    qr[144] = 0.0;

    /*reaction 146: NXC3H7 + HO2 => C3H8 + O2 */
    qf[145] = sc[8]*sc[30];
    qr[145] = 0.0;

    /*reaction 147: NXC3H7 + O2 => C3H6 + HO2 */
    qf[146] = sc[5]*sc[30];
    qr[146] = 0.0;

    /*reaction 148: NXC3H7 => H + C3H6 */
    qf[147] = sc[30];
    qr[147] = 0.0;

    /*reaction 149: H + C3H6 => NXC3H7 */
    qf[148] = sc[2]*sc[27];
    qr[148] = 0.0;

    /*reaction 150: C3H8 + OH => NXC3H7 + H2O */
    qf[149] = sc[3]*sc[31];
    qr[149] = 0.0;

    /*reaction 151: C3H8 + HO2 => NXC3H7 + H2O2 */
    qf[150] = sc[8]*sc[31];
    qr[150] = 0.0;

    /*reaction 152: H + C3H8 <=> H2 + NXC3H7 */
    qf[151] = sc[2]*sc[31];
    qr[151] = sc[4]*sc[30];

    /*reaction 153: C3H8 + OH => IXC3H7 + H2O */
    qf[152] = sc[3]*sc[31];
    qr[152] = 0.0;

    /*reaction 154: CH3 + C3H8 => CH4 + IXC3H7 */
    qf[153] = sc[11]*sc[31];
    qr[153] = 0.0;

    /*reaction 155: CH3 + C3H8 => CH4 + NXC3H7 */
    qf[154] = sc[11]*sc[31];
    qr[154] = 0.0;

    /*reaction 156: C3H8 + O => IXC3H7 + OH */
    qf[155] = sc[1]*sc[31];
    qr[155] = 0.0;

    /*reaction 157: C3H8 + HO2 => IXC3H7 + H2O2 */
    qf[156] = sc[8]*sc[31];
    qr[156] = 0.0;

    /*reaction 158: C3H8 + O => NXC3H7 + OH */
    qf[157] = sc[1]*sc[31];
    qr[157] = 0.0;

    /*reaction 159: C3H8 + O2 => IXC3H7 + HO2 */
    qf[158] = sc[5]*sc[31];
    qr[158] = 0.0;

    /*reaction 160: IXC3H7 + HO2 => C3H8 + O2 */
    qf[159] = sc[8]*sc[29];
    qr[159] = 0.0;

    /*reaction 161: H + C3H8 => H2 + IXC3H7 */
    qf[160] = sc[2]*sc[31];
    qr[160] = 0.0;

    /*reaction 162: H2 + IXC3H7 => H + C3H8 */
    qf[161] = sc[4]*sc[29];
    qr[161] = 0.0;

    /*reaction 163: C3H5O => C2H3 + CH2O */
    qf[162] = sc[28];
    qr[162] = 0.0;

    /*reaction 164: IXC3H7O2 => IXC3H7 + O2 */
    qf[163] = sc[32];
    qr[163] = 0.0;

    /*reaction 165: IXC3H7 + O2 => IXC3H7O2 */
    qf[164] = sc[5]*sc[29];
    qr[164] = 0.0;

    /*reaction 166: NXC3H7O2 => NXC3H7 + O2 */
    qf[165] = sc[33];
    qr[165] = 0.0;

    /*reaction 167: NXC3H7 + O2 => NXC3H7O2 */
    qf[166] = sc[5]*sc[30];
    qr[166] = 0.0;

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 34; ++i) {
        mixture += sc[i];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[34];
    gibbs(g_RT, tc);

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 * invT;
    amrex::Real refCinv = 1 / refC;

    /* Evaluate the kfs */
    amrex::Real k_f, k_r, Corr;
    amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

    // (0):  O + H + M => OH + M
    k_f = 1.0000000000000002e-12 * 4.72e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    qf[9] *= Corr * k_f;
    qr[9] *= Corr * k_f / (exp(g_RT[1] + g_RT[2] - g_RT[3]) * refCinv);
    // (1):  O + H2 => H + OH
    k_f = 1.0000000000000002e-06 * 50800 
               * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6292.0699999999997) * invT);
    Corr  = 1.0;
    qf[16] *= Corr * k_f;
    qr[16] *= Corr * k_f / exp(g_RT[1] - g_RT[2] - g_RT[3] + g_RT[4]);
    // (2):  H + OH => O + H2
    k_f = 1.0000000000000002e-06 * 22310 
               * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (4196.9399999999996) * invT);
    Corr  = 1.0;
    qf[17] *= Corr * k_f;
    qr[17] *= Corr * k_f / exp(-g_RT[1] + g_RT[2] + g_RT[3] - g_RT[4]);
    // (3):  2.000000 O + M => O2 + M
    k_f = 1.0000000000000002e-12 * 6170000000000000 
               * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    qf[10] *= Corr * k_f;
    qr[10] *= Corr * k_f / (exp(2.000000*g_RT[1] - g_RT[5]) * refCinv);
    // (4):  OH + H2 => H + H2O
    k_f = 1.0000000000000002e-06 * 216000000 
               * exp(1.51 * tc[0] - 0.50321666580471969 * (3429.9699999999998) * invT);
    Corr  = 1.0;
    qf[18] *= Corr * k_f;
    qr[18] *= Corr * k_f / exp(-g_RT[2] + g_RT[3] + g_RT[4] - g_RT[6]);
    // (5):  H + H2O => OH + H2
    k_f = 1.0000000000000002e-06 * 935200000 
               * exp(1.51 * tc[0] - 0.50321666580471969 * (18580.07) * invT);
    Corr  = 1.0;
    qf[19] *= Corr * k_f;
    qr[19] *= Corr * k_f / exp(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]);
    // (6):  2.000000 OH (+M) <=> H2O2 (+M)
    k_f = 1.0000000000000002e-06 * 123600000000000 
               * exp(-0.37 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    redP = Corr / k_f * 1e-12 * 3.0410000000000002e+30 
               * exp(-4.6299999999999999  * tc[0] - 0.50321666580471969  * (2049) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.46999999999999997)*exp(-tc[1] / 100) 
        + 0.46999999999999997 * exp(-tc[1]/2000)  
        + exp(-1000000000000000 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[0] *= Corr * k_f;
    qr[0] *= Corr * k_f / (exp(2.000000*g_RT[3] - g_RT[7]) * refCinv);
    // (7):  H + OH + M => H2O + M
    k_f = 1.0000000000000002e-12 * 2.2499999999999999e+22 
               * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    qf[11] *= Corr * k_f;
    qr[11] *= Corr * k_f / (exp(g_RT[2] + g_RT[3] - g_RT[6]) * refCinv);
    // (8):  O + H2O => 2.000000 OH
    k_f = 1.0000000000000002e-06 * 2970000 
               * exp(2.02 * tc[0] - 0.50321666580471969 * (13400.1) * invT);
    Corr  = 1.0;
    qf[20] *= Corr * k_f;
    qr[20] *= Corr * k_f / exp(g_RT[1] - 2.000000*g_RT[3] + g_RT[6]);
    // (9):  2.000000 OH => O + H2O
    k_f = 1.0000000000000002e-06 * 301300 
               * exp(2.02 * tc[0] - 0.50321666580471969 * (-3849.9000000000001) * invT);
    Corr  = 1.0;
    qf[21] *= Corr * k_f;
    qr[21] *= Corr * k_f / exp(-g_RT[1] + 2.000000*g_RT[3] - g_RT[6]);
    // (10):  H + O2 => O + OH
    k_f = 1.0000000000000002e-06 * 197000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (16539.91) * invT);
    Corr  = 1.0;
    qf[22] *= Corr * k_f;
    qr[22] *= Corr * k_f / exp(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[5]);
    // (11):  O + OH => H + O2
    k_f = 1.0000000000000002e-06 * 15550000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (424.94999999999999) * invT);
    Corr  = 1.0;
    qf[23] *= Corr * k_f;
    qr[23] *= Corr * k_f / exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[5]);
    // (12):  H + O2 (+M) <=> HO2 (+M)
    k_f = 1.0000000000000002e-06 * 1475000000000 
               * exp(0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    redP = Corr / k_f * 1e-12 * 35000000000000000 
               * exp(-0.40999999999999998  * tc[0] - 0.50321666580471969  * (-1115.9200000000001) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
         0. 
        + 0.5 * exp(-tc[1]/1e+30)  
        + exp(-1e+100 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[1] *= Corr * k_f;
    qr[1] *= Corr * k_f / (exp(g_RT[2] + g_RT[5] - g_RT[8]) * refCinv);
    // (13):  HO2 + H => 2.000000 OH
    k_f = 1.0000000000000002e-06 * 70800000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (299.94999999999999) * invT);
    Corr  = 1.0;
    qf[24] *= Corr * k_f;
    qr[24] *= Corr * k_f / exp(g_RT[2] - 2.000000*g_RT[3] + g_RT[8]);
    // (14):  2.000000 HO2 => H2O2 + O2
    k_f = 1.0000000000000002e-06 * 420000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (11979.92) * invT);
    Corr  = 1.0;
    qf[25] *= Corr * k_f;
    qr[25] *= Corr * k_f / exp(-g_RT[5] - g_RT[7] + 2.000000*g_RT[8]);
    // (15):  HO2 + H => H2 + O2
    k_f = 1.0000000000000002e-06 * 16600000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (820.02999999999997) * invT);
    Corr  = 1.0;
    qf[26] *= Corr * k_f;
    qr[26] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[5] + g_RT[8]);
    // (16):  HO2 + OH => H2O + O2
    k_f = 1.0000000000000002e-06 * 28900000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    Corr  = 1.0;
    qf[27] *= Corr * k_f;
    qr[27] *= Corr * k_f / exp(g_RT[3] - g_RT[5] - g_RT[6] + g_RT[8]);
    // (17):  H2O + O2 => HO2 + OH
    k_f = 1.0000000000000002e-06 * 6888000000000000 
               * exp(-0.33000000000000002 * tc[0] - 0.50321666580471969 * (72140.059999999998) * invT);
    Corr  = 1.0;
    qf[28] *= Corr * k_f;
    qr[28] *= Corr * k_f / exp(-g_RT[3] + g_RT[5] + g_RT[6] - g_RT[8]);
    // (18):  2.000000 HO2 => H2O2 + O2
    k_f = 1.0000000000000002e-06 * 130000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1629.0599999999999) * invT);
    Corr  = 1.0;
    qf[29] *= Corr * k_f;
    qr[29] *= Corr * k_f / exp(-g_RT[5] - g_RT[7] + 2.000000*g_RT[8]);
    // (19):  HO2 + O => OH + O2
    k_f = 1.0000000000000002e-06 * 32500000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[30] *= Corr * k_f;
    qr[30] *= Corr * k_f / exp(g_RT[1] - g_RT[3] - g_RT[5] + g_RT[8]);
    // (20):  OH + O2 => HO2 + O
    k_f = 1.0000000000000002e-06 * 785700000000000 
               * exp(-0.33000000000000002 * tc[0] - 0.50321666580471969 * (55390.059999999998) * invT);
    Corr  = 1.0;
    qf[31] *= Corr * k_f;
    qr[31] *= Corr * k_f / exp(-g_RT[1] + g_RT[3] + g_RT[5] - g_RT[8]);
    // (21):  H2O2 + OH => H2O + HO2
    k_f = 1.0000000000000002e-06 * 1000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[32] *= Corr * k_f;
    qr[32] *= Corr * k_f / exp(g_RT[3] - g_RT[6] + g_RT[7] - g_RT[8]);
    // (22):  H2O2 + H => H2O + OH
    k_f = 1.0000000000000002e-06 * 24100000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (3969.8899999999999) * invT);
    Corr  = 1.0;
    qf[33] *= Corr * k_f;
    qr[33] *= Corr * k_f / exp(g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7]);
    // (23):  H2O2 + OH => H2O + HO2
    k_f = 1.0000000000000002e-06 * 580000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (9559.9899999999998) * invT);
    Corr  = 1.0;
    qf[34] *= Corr * k_f;
    qr[34] *= Corr * k_f / exp(g_RT[3] - g_RT[6] + g_RT[7] - g_RT[8]);
    // (24):  H2O + HO2 => H2O2 + OH
    k_f = 1.0000000000000002e-06 * 97710000000000 
               * exp(0.33000000000000002 * tc[0] - 0.50321666580471969 * (41020.080000000002) * invT);
    Corr  = 1.0;
    qf[35] *= Corr * k_f;
    qr[35] *= Corr * k_f / exp(-g_RT[3] + g_RT[6] - g_RT[7] + g_RT[8]);
    // (25):  H2O2 + O => OH + HO2
    k_f = 1.0000000000000002e-06 * 9550000 
               * exp(2 * tc[0] - 0.50321666580471969 * (3969.8899999999999) * invT);
    Corr  = 1.0;
    qf[36] *= Corr * k_f;
    qr[36] *= Corr * k_f / exp(g_RT[1] - g_RT[3] + g_RT[7] - g_RT[8]);
    // (26):  H2O2 + H => H2 + HO2
    k_f = 1.0000000000000002e-06 * 48200000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (7950.0500000000002) * invT);
    Corr  = 1.0;
    qf[37] *= Corr * k_f;
    qr[37] *= Corr * k_f / exp(g_RT[2] - g_RT[4] + g_RT[7] - g_RT[8]);
    // (27):  H2 + HO2 => H2O2 + H
    k_f = 1.0000000000000002e-06 * 1875000000000 
               * exp(0.33000000000000002 * tc[0] - 0.50321666580471969 * (24260.040000000001) * invT);
    Corr  = 1.0;
    qf[38] *= Corr * k_f;
    qr[38] *= Corr * k_f / exp(-g_RT[2] + g_RT[4] - g_RT[7] + g_RT[8]);
    // (28):  CH2GSG + OH => CH2O + H
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[39] *= Corr * k_f;
    qr[39] *= Corr * k_f / exp(-g_RT[2] + g_RT[3] + g_RT[9] - g_RT[10]);
    // (29):  CH2GSG + H2 => CH3 + H
    k_f = 1.0000000000000002e-06 * 70000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[40] *= Corr * k_f;
    qr[40] *= Corr * k_f / exp(-g_RT[2] + g_RT[4] + g_RT[9] - g_RT[11]);
    // (30):  CH3 + H => CH2GSG + H2
    k_f = 1.0000000000000002e-06 * 2.482e+17 
               * exp(-0.89000000000000001 * tc[0] - 0.50321666580471969 * (16130.02) * invT);
    Corr  = 1.0;
    qf[41] *= Corr * k_f;
    qr[41] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[9] + g_RT[11]);
    // (31):  CH2GSG + O2 => CO + OH + H
    k_f = 1.0000000000000002e-06 * 70000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[42] *= Corr * k_f;
    qr[42] *= Corr * k_f / (exp(-g_RT[2] - g_RT[3] + g_RT[5] + g_RT[9] - g_RT[12]) * refC);
    // (32):  CH2GSG + O => CO + 2.000000 H
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[43] *= Corr * k_f;
    qr[43] *= Corr * k_f / (exp(g_RT[1] - 2.000000*g_RT[2] + g_RT[9] - g_RT[12]) * refC);
    // (33):  CH3 + HO2 => CH3O + OH
    k_f = 1.0000000000000002e-06 * 11000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[44] *= Corr * k_f;
    qr[44] *= Corr * k_f / exp(-g_RT[3] + g_RT[8] + g_RT[11] - g_RT[13]);
    // (34):  CH3 + O2 => CH2O + OH
    k_f = 1.0000000000000002e-06 * 747000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (14250) * invT);
    Corr  = 1.0;
    qf[45] *= Corr * k_f;
    qr[45] *= Corr * k_f / exp(-g_RT[3] + g_RT[5] - g_RT[10] + g_RT[11]);
    // (35):  CH3 + O2 => CH3O + O
    k_f = 1.0000000000000002e-06 * 1.995e+18 
               * exp(-1.5700000000000001 * tc[0] - 0.50321666580471969 * (29210.09) * invT);
    Corr  = 1.0;
    qf[46] *= Corr * k_f;
    qr[46] *= Corr * k_f / exp(-g_RT[1] + g_RT[5] + g_RT[11] - g_RT[13]);
    // (36):  CH3O + O => CH3 + O2
    k_f = 1.0000000000000002e-06 * 3.585e+18 
               * exp(-1.5900000000000001 * tc[0] - 0.50321666580471969 * (-1630.98) * invT);
    Corr  = 1.0;
    qf[47] *= Corr * k_f;
    qr[47] *= Corr * k_f / exp(g_RT[1] - g_RT[5] - g_RT[11] + g_RT[13]);
    // (37):  2.000000 CH3 <=> H + C2H5
    k_f = 1.0000000000000002e-06 * 6840000000000 
               * exp(0.10000000000000001 * tc[0] - 0.50321666580471969 * (10599.9) * invT);
    Corr  = 1.0;
    qf[48] *= Corr * k_f;
    qr[48] *= Corr * k_f / exp(-g_RT[2] + 2.000000*g_RT[11] - g_RT[14]);
    // (38):  CH3 + HO2 => CH4 + O2
    k_f = 1.0000000000000002e-06 * 3600000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[49] *= Corr * k_f;
    qr[49] *= Corr * k_f / exp(-g_RT[5] + g_RT[8] + g_RT[11] - g_RT[15]);
    // (39):  CH3 + O => CH2O + H
    k_f = 1.0000000000000002e-06 * 80000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[50] *= Corr * k_f;
    qr[50] *= Corr * k_f / exp(g_RT[1] - g_RT[2] - g_RT[10] + g_RT[11]);
    // (40):  CH3 + OH => CH2O + H2
    k_f = 1.0000000000000002e-06 * 22500000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (4299.9499999999998) * invT);
    Corr  = 1.0;
    qf[51] *= Corr * k_f;
    qr[51] *= Corr * k_f / exp(g_RT[3] - g_RT[4] - g_RT[10] + g_RT[11]);
    // (41):  CH3 + OH => CH2GSG + H2O
    k_f = 1.0000000000000002e-06 * 26500000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (2185.9499999999998) * invT);
    Corr  = 1.0;
    qf[52] *= Corr * k_f;
    qr[52] *= Corr * k_f / exp(g_RT[3] - g_RT[6] - g_RT[9] + g_RT[11]);
    // (42):  CH2GSG + H2O => CH3 + OH
    k_f = 1.0000000000000002e-06 * 32360000000 
               * exp(0.89000000000000001 * tc[0] - 0.50321666580471969 * (1211.04) * invT);
    Corr  = 1.0;
    qf[53] *= Corr * k_f;
    qr[53] *= Corr * k_f / exp(-g_RT[3] + g_RT[6] + g_RT[9] - g_RT[11]);
    // (43):  CH3 + H2O2 => CH4 + HO2
    k_f = 1.0000000000000002e-06 * 336500000000 
               * exp(-0.33000000000000002 * tc[0] - 0.50321666580471969 * (2501.9099999999999) * invT);
    Corr  = 1.0;
    qf[54] *= Corr * k_f;
    qr[54] *= Corr * k_f / exp(g_RT[7] - g_RT[8] + g_RT[11] - g_RT[15]);
    // (44):  CH2GSG + CH3 => C2H4 + H
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[55] *= Corr * k_f;
    qr[55] *= Corr * k_f / exp(-g_RT[2] + g_RT[9] + g_RT[11] - g_RT[16]);
    // (45):  CH3 + H (+M) => CH4 (+M)
    k_f = 1.0000000000000002e-06 * 2138000000000000 
               * exp(-0.40000000000000002 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[4] + ( 5 - 1)*sc[6] + ( 3 - 1)*sc[18] + ( 2 - 1)*sc[12];
    redP = Corr / k_f * 1e-12 * 3.31e+30 
               * exp(-4  * tc[0] - 0.50321666580471969  * (2108.0300000000002) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
         0. 
        + 0. 
        + exp(-40 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[2] *= Corr * k_f;
    qr[2] *= Corr * k_f / (exp(g_RT[2] + g_RT[11] - g_RT[15]) * refCinv);
    // (46):  2.000000 CH3 (+M) <=> C2H6 (+M)
    k_f = 1.0000000000000002e-06 * 92140000000000000 
               * exp(-1.1699999999999999 * tc[0] - 0.50321666580471969 * (635.75999999999999) * invT);
    Corr  = mixture + ( 2 - 1)*sc[4] + ( 5 - 1)*sc[6] + ( 3 - 1)*sc[18] + ( 2 - 1)*sc[12];
    redP = Corr / k_f * 1e-12 * 1.135e+36 
               * exp(-5.2460000000000004  * tc[0] - 0.50321666580471969  * (1705.0699999999999) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.40500000000000003)*exp(-tc[1] / 1120) 
        + 0.40500000000000003 * exp(-tc[1]/69.599999999999994)  
        + exp(-1000000000000000 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[3] *= Corr * k_f;
    qr[3] *= Corr * k_f / (exp(2.000000*g_RT[11] - g_RT[17]) * refCinv);
    // (47):  CH2GSG + CH4 => 2.000000 CH3
    k_f = 1.0000000000000002e-06 * 40000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[56] *= Corr * k_f;
    qr[56] *= Corr * k_f / exp(g_RT[9] - 2.000000*g_RT[11] + g_RT[15]);
    // (48):  2.000000 CH3 => CH2GSG + CH4
    k_f = 1.0000000000000002e-06 * 5429000000000000 
               * exp(-0.89000000000000001 * tc[0] - 0.50321666580471969 * (15650.1) * invT);
    Corr  = 1.0;
    qf[57] *= Corr * k_f;
    qr[57] *= Corr * k_f / exp(-g_RT[9] + 2.000000*g_RT[11] - g_RT[15]);
    // (49):  CH4 + O => CH3 + OH
    k_f = 1.0000000000000002e-06 * 3150000000000 
               * exp(0.5 * tc[0] - 0.50321666580471969 * (10289.91) * invT);
    Corr  = 1.0;
    qf[58] *= Corr * k_f;
    qr[58] *= Corr * k_f / exp(g_RT[1] - g_RT[3] - g_RT[11] + g_RT[15]);
    // (50):  CH3 + OH => CH4 + O
    k_f = 1.0000000000000002e-06 * 52960000000 
               * exp(0.5 * tc[0] - 0.50321666580471969 * (7715.1099999999997) * invT);
    Corr  = 1.0;
    qf[59] *= Corr * k_f;
    qr[59] *= Corr * k_f / exp(-g_RT[1] + g_RT[3] + g_RT[11] - g_RT[15]);
    // (51):  CH4 + H => CH3 + H2
    k_f = 1.0000000000000002e-06 * 17270 
               * exp(3 * tc[0] - 0.50321666580471969 * (8223.9500000000007) * invT);
    Corr  = 1.0;
    qf[60] *= Corr * k_f;
    qr[60] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[11] + g_RT[15]);
    // (52):  CH3 + H2 => CH4 + H
    k_f = 1.0000000000000002e-06 * 661 
               * exp(3 * tc[0] - 0.50321666580471969 * (7744.0200000000004) * invT);
    Corr  = 1.0;
    qf[61] *= Corr * k_f;
    qr[61] *= Corr * k_f / exp(-g_RT[2] + g_RT[4] + g_RT[11] - g_RT[15]);
    // (53):  CH4 + OH => CH3 + H2O
    k_f = 1.0000000000000002e-06 * 193000 
               * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (2106.1199999999999) * invT);
    Corr  = 1.0;
    qf[62] *= Corr * k_f;
    qr[62] *= Corr * k_f / exp(g_RT[3] - g_RT[6] - g_RT[11] + g_RT[15]);
    // (54):  CH3 + H2O => CH4 + OH
    k_f = 1.0000000000000002e-06 * 31990 
               * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (16780.110000000001) * invT);
    Corr  = 1.0;
    qf[63] *= Corr * k_f;
    qr[63] *= Corr * k_f / exp(-g_RT[3] + g_RT[6] + g_RT[11] - g_RT[15]);
    // (55):  CO + HO2 => CO2 + OH
    k_f = 1.0000000000000002e-06 * 30100000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (23000) * invT);
    Corr  = 1.0;
    qf[64] *= Corr * k_f;
    qr[64] *= Corr * k_f / exp(-g_RT[3] + g_RT[8] + g_RT[12] - g_RT[18]);
    // (56):  CO + O2 => CO2 + O
    k_f = 1.0000000000000002e-06 * 16200000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (47700.050000000003) * invT);
    Corr  = 1.0;
    qf[65] *= Corr * k_f;
    qr[65] *= Corr * k_f / exp(-g_RT[1] + g_RT[5] + g_RT[12] - g_RT[18]);
    // (57):  CO2 + O => CO + O2
    k_f = 1.0000000000000002e-06 * 143300000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (53919.93) * invT);
    Corr  = 1.0;
    qf[66] *= Corr * k_f;
    qr[66] *= Corr * k_f / exp(g_RT[1] - g_RT[5] - g_RT[12] + g_RT[18]);
    // (58):  CO + OH => CO2 + H
    k_f = 1.0000000000000002e-06 * 140000 
               * exp(1.95 * tc[0] - 0.50321666580471969 * (-1347.04) * invT);
    Corr  = 1.0;
    qf[67] *= Corr * k_f;
    qr[67] *= Corr * k_f / exp(-g_RT[2] + g_RT[3] + g_RT[12] - g_RT[18]);
    // (59):  CO2 + H => CO + OH
    k_f = 1.0000000000000002e-06 * 15680000 
               * exp(1.95 * tc[0] - 0.50321666580471969 * (20989.959999999999) * invT);
    Corr  = 1.0;
    qf[68] *= Corr * k_f;
    qr[68] *= Corr * k_f / exp(g_RT[2] - g_RT[3] - g_RT[12] + g_RT[18]);
    // (60):  CO + O (+M) => CO2 (+M)
    k_f = 1.0000000000000002e-06 * 18000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (2384.0799999999999) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    redP = Corr / k_f * 1e-12 * 1.35e+24 
               * exp(-2.7879999999999998  * tc[0] - 0.50321666580471969  * (4190.9700000000003) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-1)*exp(-tc[1] / 1) 
        + 1 * exp(-tc[1]/10000000)  
        + exp(-10000000 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[4] *= Corr * k_f;
    qr[4] *= Corr * k_f / (exp(g_RT[1] + g_RT[12] - g_RT[18]) * refCinv);
    // (61):  HCO + CH3 => CH4 + CO
    k_f = 1.0000000000000002e-06 * 121000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[69] *= Corr * k_f;
    qr[69] *= Corr * k_f / exp(g_RT[11] - g_RT[12] - g_RT[15] + g_RT[19]);
    // (62):  HCO + H => CO + H2
    k_f = 1.0000000000000002e-06 * 73400000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[70] *= Corr * k_f;
    qr[70] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[12] + g_RT[19]);
    // (63):  HCO + O2 => CO + HO2
    k_f = 1.0000000000000002e-06 * 7580000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (409.88999999999999) * invT);
    Corr  = 1.0;
    qf[71] *= Corr * k_f;
    qr[71] *= Corr * k_f / exp(g_RT[5] - g_RT[8] - g_RT[12] + g_RT[19]);
    // (64):  HCO + O => CO + OH
    k_f = 1.0000000000000002e-06 * 30200000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[72] *= Corr * k_f;
    qr[72] *= Corr * k_f / exp(g_RT[1] - g_RT[3] - g_RT[12] + g_RT[19]);
    // (65):  HCO + O => CO2 + H
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[73] *= Corr * k_f;
    qr[73] *= Corr * k_f / exp(g_RT[1] - g_RT[2] - g_RT[18] + g_RT[19]);
    // (66):  HCO + M => H + CO + M
    k_f = 1.0000000000000002e-06 * 1.86e+17 
               * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[4] + ( 6 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    qf[12] *= Corr * k_f;
    qr[12] *= Corr * k_f / (exp(-g_RT[2] - g_RT[12] + g_RT[19]) * refC);
    // (67):  H + CO + M => HCO + M
    k_f = 1.0000000000000002e-12 * 64670000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-441.92000000000002) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[4] + ( 6 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    qf[13] *= Corr * k_f;
    qr[13] *= Corr * k_f / (exp(g_RT[2] + g_RT[12] - g_RT[19]) * refCinv);
    // (68):  HCO + OH => CO + H2O
    k_f = 1.0000000000000002e-06 * 102000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[74] *= Corr * k_f;
    qr[74] *= Corr * k_f / exp(g_RT[3] - g_RT[6] - g_RT[12] + g_RT[19]);
    // (69):  CH2O + H => HCO + H2
    k_f = 1.0000000000000002e-06 * 933400000 
               * exp(1.5 * tc[0] - 0.50321666580471969 * (2976.0999999999999) * invT);
    Corr  = 1.0;
    qf[75] *= Corr * k_f;
    qr[75] *= Corr * k_f / exp(g_RT[2] - g_RT[4] + g_RT[10] - g_RT[19]);
    // (70):  CH2O + O2 => HCO + HO2
    k_f = 1.0000000000000002e-06 * 20500000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (38950.050000000003) * invT);
    Corr  = 1.0;
    qf[76] *= Corr * k_f;
    qr[76] *= Corr * k_f / exp(g_RT[5] - g_RT[8] + g_RT[10] - g_RT[19]);
    // (71):  CH2O + OH => HCO + H2O
    k_f = 1.0000000000000002e-06 * 3430000000 
               * exp(1.1799999999999999 * tc[0] - 0.50321666580471969 * (-446.94) * invT);
    Corr  = 1.0;
    qf[77] *= Corr * k_f;
    qr[77] *= Corr * k_f / exp(g_RT[3] - g_RT[6] + g_RT[10] - g_RT[19]);
    // (72):  CH2O + HO2 => HCO + H2O2
    k_f = 1.0000000000000002e-06 * 0.0058199999999999997 
               * exp(4.5300000000000002 * tc[0] - 0.50321666580471969 * (6556.8800000000001) * invT);
    Corr  = 1.0;
    qf[78] *= Corr * k_f;
    qr[78] *= Corr * k_f / exp(-g_RT[7] + g_RT[8] + g_RT[10] - g_RT[19]);
    // (73):  CH2O + O => HCO + OH
    k_f = 1.0000000000000002e-06 * 416000000000 
               * exp(0.56999999999999995 * tc[0] - 0.50321666580471969 * (2761.9499999999998) * invT);
    Corr  = 1.0;
    qf[79] *= Corr * k_f;
    qr[79] *= Corr * k_f / exp(g_RT[1] - g_RT[3] + g_RT[10] - g_RT[19]);
    // (74):  CH2O + CH3 => HCO + CH4
    k_f = 1.0000000000000002e-06 * 3.636e-06 
               * exp(5.4199999999999999 * tc[0] - 0.50321666580471969 * (998.09000000000003) * invT);
    Corr  = 1.0;
    qf[80] *= Corr * k_f;
    qr[80] *= Corr * k_f / exp(g_RT[10] + g_RT[11] - g_RT[15] - g_RT[19]);
    // (75):  CH3O (+M) => CH2O + H (+M)
    k_f = 1 * 54500000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (13500) * invT);
    Corr  = mixture;
    redP = Corr / k_f * 1e-6 * 2.344e+25 
               * exp(-2.7000000000000002  * tc[0] - 0.50321666580471969  * (30599.900000000001) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-1)*exp(-tc[1] / 1) 
        + 1 * exp(-tc[1]/10000000)  
        + exp(-10000000 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[5] *= Corr * k_f;
    qr[5] *= Corr * k_f / (exp(-g_RT[2] - g_RT[10] + g_RT[13]) * refC);
    // (76):  CH3O + O2 => CH2O + HO2
    k_f = 1.0000000000000002e-06 * 55000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (2424) * invT);
    Corr  = 1.0;
    qf[81] *= Corr * k_f;
    qr[81] *= Corr * k_f / exp(g_RT[5] - g_RT[8] - g_RT[10] + g_RT[13]);
    // (77):  CH2GSG + CO2 => CH2O + CO
    k_f = 1.0000000000000002e-06 * 3000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[82] *= Corr * k_f;
    qr[82] *= Corr * k_f / exp(g_RT[9] - g_RT[10] - g_RT[12] + g_RT[18]);
    // (78):  CH3O2 + CH3 => 2.000000 CH3O
    k_f = 1.0000000000000002e-06 * 7000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1000) * invT);
    Corr  = 1.0;
    qf[83] *= Corr * k_f;
    qr[83] *= Corr * k_f / exp(g_RT[11] - 2.000000*g_RT[13] + g_RT[20]);
    // (79):  CH3O2 + HO2 => CH3O2H + O2
    k_f = 1.0000000000000002e-06 * 17500000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-3275.0999999999999) * invT);
    Corr  = 1.0;
    qf[84] *= Corr * k_f;
    qr[84] *= Corr * k_f / exp(-g_RT[5] + g_RT[8] + g_RT[20] - g_RT[21]);
    // (80):  2.000000 CH3O2 => O2 + 2.000000 CH3O
    k_f = 1.0000000000000002e-06 * 14000000000000000 
               * exp(-1.6100000000000001 * tc[0] - 0.50321666580471969 * (1859.9400000000001) * invT);
    Corr  = 1.0;
    qf[85] *= Corr * k_f;
    qr[85] *= Corr * k_f / (exp(-g_RT[5] - 2.000000*g_RT[13] + 2.000000*g_RT[20]) * refC);
    // (81):  CH3O2 + CH2O => CH3O2H + HCO
    k_f = 1.0000000000000002e-06 * 1990000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (11659.889999999999) * invT);
    Corr  = 1.0;
    qf[86] *= Corr * k_f;
    qr[86] *= Corr * k_f / exp(g_RT[10] - g_RT[19] + g_RT[20] - g_RT[21]);
    // (82):  CH3O2 + M => CH3 + O2 + M
    k_f = 1.0000000000000002e-06 * 4.3430000000000002e+27 
               * exp(-3.4199999999999999 * tc[0] - 0.50321666580471969 * (30469.889999999999) * invT);
    Corr  = mixture;
    qf[14] *= Corr * k_f;
    qr[14] *= Corr * k_f / (exp(-g_RT[5] - g_RT[11] + g_RT[20]) * refC);
    // (83):  CH3 + O2 + M => CH3O2 + M
    k_f = 1.0000000000000002e-12 * 5.4400000000000001e+25 
               * exp(-3.2999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture;
    qf[15] *= Corr * k_f;
    qr[15] *= Corr * k_f / (exp(g_RT[5] + g_RT[11] - g_RT[20]) * refCinv);
    // (84):  CH3O2H => CH3O + OH
    k_f = 1 * 631000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (42299.949999999997) * invT);
    Corr  = 1.0;
    qf[87] *= Corr * k_f;
    qr[87] *= Corr * k_f / (exp(-g_RT[3] - g_RT[13] + g_RT[21]) * refC);
    // (85):  CH3O + OH => CH3O2H
    k_f = 1.0000000000000002e-06 * 116600000000 
               * exp(0.59999999999999998 * tc[0] - 0.50321666580471969 * (-1771.03) * invT);
    Corr  = 1.0;
    qf[88] *= Corr * k_f;
    qr[88] *= Corr * k_f / (exp(g_RT[3] + g_RT[13] - g_RT[21]) * refCinv);
    // (86):  C2H2 + O2 => HCCO + OH
    k_f = 1.0000000000000002e-06 * 200000000 
               * exp(1.5 * tc[0] - 0.50321666580471969 * (30099.900000000001) * invT);
    Corr  = 1.0;
    qf[89] *= Corr * k_f;
    qr[89] *= Corr * k_f / exp(-g_RT[3] + g_RT[5] + g_RT[22] - g_RT[23]);
    // (87):  C2H2 + O => HCCO + H
    k_f = 1.0000000000000002e-06 * 14300000 
               * exp(2 * tc[0] - 0.50321666580471969 * (1900.0999999999999) * invT);
    Corr  = 1.0;
    qf[90] *= Corr * k_f;
    qr[90] *= Corr * k_f / exp(g_RT[1] - g_RT[2] + g_RT[22] - g_RT[23]);
    // (88):  C2H3 + O2 => CH2CHO + O
    k_f = 1.0000000000000002e-06 * 350000000000000 
               * exp(-0.60999999999999999 * tc[0] - 0.50321666580471969 * (5260.04) * invT);
    Corr  = 1.0;
    qf[91] *= Corr * k_f;
    qr[91] *= Corr * k_f / exp(-g_RT[1] + g_RT[5] + g_RT[24] - g_RT[25]);
    // (89):  C2H3 + H => C2H2 + H2
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (2500) * invT);
    Corr  = 1.0;
    qf[92] *= Corr * k_f;
    qr[92] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[22] + g_RT[24]);
    // (90):  C2H3 + O2 => CH2O + HCO
    k_f = 1.0000000000000002e-06 * 1.6999999999999999e+29 
               * exp(-5.3099999999999996 * tc[0] - 0.50321666580471969 * (6500) * invT);
    Corr  = 1.0;
    qf[93] *= Corr * k_f;
    qr[93] *= Corr * k_f / exp(g_RT[5] - g_RT[10] - g_RT[19] + g_RT[24]);
    // (91):  C2H3 + O2 => C2H2 + HO2
    k_f = 1.0000000000000002e-06 * 2.12e-06 
               * exp(6 * tc[0] - 0.50321666580471969 * (9483.9899999999998) * invT);
    Corr  = 1.0;
    qf[94] *= Corr * k_f;
    qr[94] *= Corr * k_f / exp(g_RT[5] - g_RT[8] - g_RT[22] + g_RT[24]);
    // (92):  C2H3 (+M) => H + C2H2 (+M)
    k_f = 1 * 16060000000 
               * exp(1.028 * tc[0] - 0.50321666580471969 * (40503.589999999997) * invT);
    Corr  = mixture + ( 2 - 1)*sc[4] + ( 5 - 1)*sc[6] + ( 3 - 1)*sc[18] + ( 2 - 1)*sc[12];
    redP = Corr / k_f * 1e-6 * 1.164e+39 
               * exp(-6.8209999999999997  * tc[0] - 0.50321666580471969  * (44491.629999999997) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
         0. 
        + 1 * exp(-tc[1]/675)  
        + exp(-1000000000000000 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[6] *= Corr * k_f;
    qr[6] *= Corr * k_f / (exp(-g_RT[2] - g_RT[22] + g_RT[24]) * refC);
    // (93):  C2H4 + O => CH3 + HCO
    k_f = 1.0000000000000002e-06 * 10200000 
               * exp(1.8799999999999999 * tc[0] - 0.50321666580471969 * (179.02000000000001) * invT);
    Corr  = 1.0;
    qf[95] *= Corr * k_f;
    qr[95] *= Corr * k_f / exp(g_RT[1] - g_RT[11] + g_RT[16] - g_RT[19]);
    // (94):  H + C2H4 (+M) <=> C2H5 (+M)
    k_f = 1.0000000000000002e-06 * 1081000000000 
               * exp(0.45000000000000001 * tc[0] - 0.50321666580471969 * (1821.9400000000001) * invT);
    Corr  = mixture;
    redP = Corr / k_f * 1e-12 * 1.1120000000000001e+34 
               * exp(-5  * tc[0] - 0.50321666580471969  * (4447.8999999999996) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
         0. 
        + 1 * exp(-tc[1]/95)  
        + exp(-200 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[7] *= Corr * k_f;
    qr[7] *= Corr * k_f / (exp(g_RT[2] - g_RT[14] + g_RT[16]) * refCinv);
    // (95):  C2H4 + OH => C2H3 + H2O
    k_f = 1.0000000000000002e-06 * 20500000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (5950.0500000000002) * invT);
    Corr  = 1.0;
    qf[96] *= Corr * k_f;
    qr[96] *= Corr * k_f / exp(g_RT[3] - g_RT[6] + g_RT[16] - g_RT[24]);
    // (96):  C2H3 + H2O => C2H4 + OH
    k_f = 1.0000000000000002e-06 * 6033000000000000 
               * exp(-0.82999999999999996 * tc[0] - 0.50321666580471969 * (21760.040000000001) * invT);
    Corr  = 1.0;
    qf[97] *= Corr * k_f;
    qr[97] *= Corr * k_f / exp(-g_RT[3] + g_RT[6] - g_RT[16] + g_RT[24]);
    // (97):  C2H4 + H => C2H3 + H2
    k_f = 1.0000000000000002e-06 * 0.0084200000000000004 
               * exp(4.6200000000000001 * tc[0] - 0.50321666580471969 * (2582.9299999999998) * invT);
    Corr  = 1.0;
    qf[98] *= Corr * k_f;
    qr[98] *= Corr * k_f / exp(g_RT[2] - g_RT[4] + g_RT[16] - g_RT[24]);
    // (98):  C2H3 + H2 => C2H4 + H
    k_f = 1.0000000000000002e-06 * 0.57230000000000003 
               * exp(3.79 * tc[0] - 0.50321666580471969 * (3233.0300000000002) * invT);
    Corr  = 1.0;
    qf[99] *= Corr * k_f;
    qr[99] *= Corr * k_f / exp(-g_RT[2] + g_RT[4] - g_RT[16] + g_RT[24]);
    // (99):  C2H4 + CH3O2 => C2H3 + CH3O2H
    k_f = 1.0000000000000002e-06 * 2230000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (17190.009999999998) * invT);
    Corr  = 1.0;
    qf[100] *= Corr * k_f;
    qr[100] *= Corr * k_f / exp(g_RT[16] + g_RT[20] - g_RT[21] - g_RT[24]);
    // (100):  C2H4 + CH3 => C2H3 + CH4
    k_f = 1.0000000000000002e-06 * 6.6200000000000001 
               * exp(3.7000000000000002 * tc[0] - 0.50321666580471969 * (9500) * invT);
    Corr  = 1.0;
    qf[101] *= Corr * k_f;
    qr[101] *= Corr * k_f / exp(g_RT[11] - g_RT[15] + g_RT[16] - g_RT[24]);
    // (101):  C2H3 + CH4 => C2H4 + CH3
    k_f = 1.0000000000000002e-06 * 1.4399999999999999 
               * exp(4.0199999999999996 * tc[0] - 0.50321666580471969 * (5472.04) * invT);
    Corr  = 1.0;
    qf[102] *= Corr * k_f;
    qr[102] *= Corr * k_f / exp(-g_RT[11] + g_RT[15] - g_RT[16] + g_RT[24]);
    // (102):  C2H4 + O => CH2CHO + H
    k_f = 1.0000000000000002e-06 * 3390000 
               * exp(1.8799999999999999 * tc[0] - 0.50321666580471969 * (179.02000000000001) * invT);
    Corr  = 1.0;
    qf[103] *= Corr * k_f;
    qr[103] *= Corr * k_f / exp(g_RT[1] - g_RT[2] + g_RT[16] - g_RT[25]);
    // (103):  C2H4 (+M) => C2H2 + H2 (+M)
    k_f = 1 * 18000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (76000) * invT);
    Corr  = mixture;
    redP = Corr / k_f * 1e-6 * 1500000000000000 
               * exp(0  * tc[0] - 0.50321666580471969  * (55440.010000000002) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-1)*exp(-tc[1] / 1) 
        + 1 * exp(-tc[1]/10000000)  
        + exp(-10000000 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[8] *= Corr * k_f;
    qr[8] *= Corr * k_f / (exp(-g_RT[4] + g_RT[16] - g_RT[22]) * refC);
    // (104):  CH3 + C2H5 => CH4 + C2H4
    k_f = 1.0000000000000002e-06 * 19500000000000 
               * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[104] *= Corr * k_f;
    qr[104] *= Corr * k_f / exp(g_RT[11] + g_RT[14] - g_RT[15] - g_RT[16]);
    // (105):  C2H5 + O2 => C2H4 + HO2
    k_f = 1.0000000000000002e-06 * 1.22e+30 
               * exp(-5.7599999999999998 * tc[0] - 0.50321666580471969 * (10099.9) * invT);
    Corr  = 1.0;
    qf[105] *= Corr * k_f;
    qr[105] *= Corr * k_f / exp(g_RT[5] - g_RT[8] + g_RT[14] - g_RT[16]);
    // (106):  C2H4 + HO2 => C2H5 + O2
    k_f = 1.0000000000000002e-06 * 1.259e+30 
               * exp(-5.6299999999999999 * tc[0] - 0.50321666580471969 * (22299.950000000001) * invT);
    Corr  = 1.0;
    qf[106] *= Corr * k_f;
    qr[106] *= Corr * k_f / exp(-g_RT[5] + g_RT[8] - g_RT[14] + g_RT[16]);
    // (107):  C2H5 + HO2 => C2H6 + O2
    k_f = 1.0000000000000002e-06 * 267900000 
               * exp(0.89000000000000001 * tc[0] - 0.50321666580471969 * (-1922.0799999999999) * invT);
    Corr  = 1.0;
    qf[107] *= Corr * k_f;
    qr[107] *= Corr * k_f / exp(-g_RT[5] + g_RT[8] + g_RT[14] - g_RT[17]);
    // (108):  H + C2H5 => C2H6
    k_f = 1.0000000000000002e-06 * 583100000000 
               * exp(0.59899999999999998 * tc[0] - 0.50321666580471969 * (-2913) * invT);
    Corr  = 1.0;
    qf[108] *= Corr * k_f;
    qr[108] *= Corr * k_f / (exp(g_RT[2] + g_RT[14] - g_RT[17]) * refCinv);
    // (109):  C2H6 + H => C2H5 + H2
    k_f = 1.0000000000000002e-06 * 554 
               * exp(3.5 * tc[0] - 0.50321666580471969 * (5167.0699999999997) * invT);
    Corr  = 1.0;
    qf[109] *= Corr * k_f;
    qr[109] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[14] + g_RT[17]);
    // (110):  C2H5 + H2 => C2H6 + H
    k_f = 1.0000000000000002e-06 * 0.13550000000000001 
               * exp(4.0599999999999996 * tc[0] - 0.50321666580471969 * (8857.0699999999997) * invT);
    Corr  = 1.0;
    qf[110] *= Corr * k_f;
    qr[110] *= Corr * k_f / exp(-g_RT[2] + g_RT[4] + g_RT[14] - g_RT[17]);
    // (111):  C2H6 + OH => C2H5 + H2O
    k_f = 1.0000000000000002e-06 * 58000000 
               * exp(1.73 * tc[0] - 0.50321666580471969 * (1159.8900000000001) * invT);
    Corr  = 1.0;
    qf[111] *= Corr * k_f;
    qr[111] *= Corr * k_f / exp(g_RT[3] - g_RT[6] - g_RT[14] + g_RT[17]);
    // (112):  CH2GSG + C2H6 => CH3 + C2H5
    k_f = 1.0000000000000002e-06 * 120000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[112] *= Corr * k_f;
    qr[112] *= Corr * k_f / exp(g_RT[9] - g_RT[11] - g_RT[14] + g_RT[17]);
    // (113):  C2H6 + O => C2H5 + OH
    k_f = 1.0000000000000002e-06 * 13000000 
               * exp(2.1299999999999999 * tc[0] - 0.50321666580471969 * (5190.0100000000002) * invT);
    Corr  = 1.0;
    qf[113] *= Corr * k_f;
    qr[113] *= Corr * k_f / exp(g_RT[1] - g_RT[3] - g_RT[14] + g_RT[17]);
    // (114):  C2H6 + CH3 => C2H5 + CH4
    k_f = 1.0000000000000002e-06 * 1.5099999999999999e-07 
               * exp(6 * tc[0] - 0.50321666580471969 * (6047.0799999999999) * invT);
    Corr  = 1.0;
    qf[114] *= Corr * k_f;
    qr[114] *= Corr * k_f / exp(g_RT[11] - g_RT[14] - g_RT[15] + g_RT[17]);
    // (115):  HCCO + O => H + 2.000000 CO
    k_f = 1.0000000000000002e-06 * 80000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[115] *= Corr * k_f;
    qr[115] *= Corr * k_f / (exp(g_RT[1] - g_RT[2] - 2.000000*g_RT[12] + g_RT[23]) * refC);
    // (116):  HCCO + O2 => CO2 + HCO
    k_f = 1.0000000000000002e-06 * 240000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-853.97000000000003) * invT);
    Corr  = 1.0;
    qf[116] *= Corr * k_f;
    qr[116] *= Corr * k_f / exp(g_RT[5] - g_RT[18] - g_RT[19] + g_RT[23]);
    // (117):  HCCO + OH => 2.000000 HCO
    k_f = 1.0000000000000002e-06 * 10000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[117] *= Corr * k_f;
    qr[117] *= Corr * k_f / exp(g_RT[3] - 2.000000*g_RT[19] + g_RT[23]);
    // (118):  HCCO + H => CH2GSG + CO
    k_f = 1.0000000000000002e-06 * 110000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[118] *= Corr * k_f;
    qr[118] *= Corr * k_f / exp(g_RT[2] - g_RT[9] - g_RT[12] + g_RT[23]);
    // (119):  CH2GSG + CO => HCCO + H
    k_f = 1.0000000000000002e-06 * 2046000000000 
               * exp(0.89000000000000001 * tc[0] - 0.50321666580471969 * (27830.07) * invT);
    Corr  = 1.0;
    qf[119] *= Corr * k_f;
    qr[119] *= Corr * k_f / exp(-g_RT[2] + g_RT[9] + g_RT[12] - g_RT[23]);
    // (120):  CH2CHO + O2 => CH2O + CO + OH
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (4200.0500000000002) * invT);
    Corr  = 1.0;
    qf[120] *= Corr * k_f;
    qr[120] *= Corr * k_f / (exp(-g_RT[3] + g_RT[5] - g_RT[10] - g_RT[12] + g_RT[25]) * refC);
    // (121):  C3H5XA + CH2O => C3H6 + HCO
    k_f = 1.0000000000000002e-06 * 630000000 
               * exp(1.8999999999999999 * tc[0] - 0.50321666580471969 * (18190.009999999998) * invT);
    Corr  = 1.0;
    qf[121] *= Corr * k_f;
    qr[121] *= Corr * k_f / exp(g_RT[10] - g_RT[19] + g_RT[26] - g_RT[27]);
    // (122):  C3H6 + HCO => C3H5XA + CH2O
    k_f = 1.0000000000000002e-06 * 109700000 
               * exp(1.8899999999999999 * tc[0] - 0.50321666580471969 * (15840.110000000001) * invT);
    Corr  = 1.0;
    qf[122] *= Corr * k_f;
    qr[122] *= Corr * k_f / exp(-g_RT[10] + g_RT[19] - g_RT[26] + g_RT[27]);
    // (123):  C3H5XA + HO2 => C3H5O + OH
    k_f = 1.0000000000000002e-06 * 7000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1000) * invT);
    Corr  = 1.0;
    qf[123] *= Corr * k_f;
    qr[123] *= Corr * k_f / exp(-g_RT[3] + g_RT[8] + g_RT[26] - g_RT[28]);
    // (124):  C3H5XA + CH3O2 => C3H5O + CH3O
    k_f = 1.0000000000000002e-06 * 7000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1000) * invT);
    Corr  = 1.0;
    qf[124] *= Corr * k_f;
    qr[124] *= Corr * k_f / exp(-g_RT[13] + g_RT[20] + g_RT[26] - g_RT[28]);
    // (125):  C3H5XA + HO2 => C3H6 + O2
    k_f = 1.0000000000000002e-06 * 33320000000 
               * exp(0.34000000000000002 * tc[0] - 0.50321666580471969 * (-555.92999999999995) * invT);
    Corr  = 1.0;
    qf[125] *= Corr * k_f;
    qr[125] *= Corr * k_f / exp(-g_RT[5] + g_RT[8] + g_RT[26] - g_RT[27]);
    // (126):  C3H5XA + O2 => CH2CHO + CH2O
    k_f = 1.0000000000000002e-06 * 7140000000000000 
               * exp(-1.21 * tc[0] - 0.50321666580471969 * (21049.950000000001) * invT);
    Corr  = 1.0;
    qf[126] *= Corr * k_f;
    qr[126] *= Corr * k_f / exp(g_RT[5] - g_RT[10] - g_RT[25] + g_RT[26]);
    // (127):  C3H5XA + O2 => C2H2 + CH2O + OH
    k_f = 1.0000000000000002e-06 * 9.7200000000000003e+29 
               * exp(-5.71 * tc[0] - 0.50321666580471969 * (21450.049999999999) * invT);
    Corr  = 1.0;
    qf[127] *= Corr * k_f;
    qr[127] *= Corr * k_f / (exp(-g_RT[3] + g_RT[5] - g_RT[10] - g_RT[22] + g_RT[26]) * refC);
    // (128):  C3H5XA + H => C3H6
    k_f = 1.0000000000000002e-06 * 4.8869999999999999e+56 
               * exp(-12.25 * tc[0] - 0.50321666580471969 * (28080.07) * invT);
    Corr  = 1.0;
    qf[128] *= Corr * k_f;
    qr[128] *= Corr * k_f / (exp(g_RT[2] + g_RT[26] - g_RT[27]) * refCinv);
    // (129):  C3H5XA => C2H2 + CH3
    k_f = 1 * 2.397e+48 
               * exp(-9.9000000000000004 * tc[0] - 0.50321666580471969 * (82080.070000000007) * invT);
    Corr  = 1.0;
    qf[129] *= Corr * k_f;
    qr[129] *= Corr * k_f / (exp(-g_RT[11] - g_RT[22] + g_RT[26]) * refC);
    // (130):  C3H6 => C2H3 + CH3
    k_f = 1 * 2.7300000000000001e+62 
               * exp(-13.279999999999999 * tc[0] - 0.50321666580471969 * (123200.05) * invT);
    Corr  = 1.0;
    qf[130] *= Corr * k_f;
    qr[130] *= Corr * k_f / (exp(-g_RT[11] - g_RT[24] + g_RT[27]) * refC);
    // (131):  C2H3 + CH3 => C3H6
    k_f = 1.0000000000000002e-06 * 4.7119999999999996e+59 
               * exp(-13.19 * tc[0] - 0.50321666580471969 * (29539.91) * invT);
    Corr  = 1.0;
    qf[131] *= Corr * k_f;
    qr[131] *= Corr * k_f / (exp(g_RT[11] + g_RT[24] - g_RT[27]) * refCinv);
    // (132):  C3H6 + O => C2H5 + HCO
    k_f = 1.0000000000000002e-06 * 15800000 
               * exp(1.76 * tc[0] - 0.50321666580471969 * (-1216.0599999999999) * invT);
    Corr  = 1.0;
    qf[132] *= Corr * k_f;
    qr[132] *= Corr * k_f / exp(g_RT[1] - g_RT[14] - g_RT[19] + g_RT[27]);
    // (133):  C3H6 + H => C2H4 + CH3
    k_f = 1.0000000000000002e-06 * 4.8299999999999998e+33 
               * exp(-5.8099999999999996 * tc[0] - 0.50321666580471969 * (18500) * invT);
    Corr  = 1.0;
    qf[133] *= Corr * k_f;
    qr[133] *= Corr * k_f / exp(g_RT[2] - g_RT[11] - g_RT[16] + g_RT[27]);
    // (134):  C2H4 + CH3 => C3H6 + H
    k_f = 1.0000000000000002e-06 * 2.3130000000000001e+33 
               * exp(-5.9000000000000004 * tc[0] - 0.50321666580471969 * (31619.98) * invT);
    Corr  = 1.0;
    qf[134] *= Corr * k_f;
    qr[134] *= Corr * k_f / exp(-g_RT[2] + g_RT[11] + g_RT[16] - g_RT[27]);
    // (135):  C3H6 + H => C3H5XA + H2
    k_f = 1.0000000000000002e-06 * 173000 
               * exp(2.5 * tc[0] - 0.50321666580471969 * (2492.1100000000001) * invT);
    Corr  = 1.0;
    qf[135] *= Corr * k_f;
    qr[135] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[26] + g_RT[27]);
    // (136):  C3H5XA + H2 => C3H6 + H
    k_f = 1.0000000000000002e-06 * 79330 
               * exp(2.5099999999999998 * tc[0] - 0.50321666580471969 * (19520.080000000002) * invT);
    Corr  = 1.0;
    qf[136] *= Corr * k_f;
    qr[136] *= Corr * k_f / exp(-g_RT[2] + g_RT[4] + g_RT[26] - g_RT[27]);
    // (137):  C3H6 + OH => C3H5XA + H2O
    k_f = 1.0000000000000002e-06 * 3120000 
               * exp(2 * tc[0] - 0.50321666580471969 * (-298.04000000000002) * invT);
    Corr  = 1.0;
    qf[137] *= Corr * k_f;
    qr[137] *= Corr * k_f / exp(g_RT[3] - g_RT[6] - g_RT[26] + g_RT[27]);
    // (138):  C3H6 + O => C3H5XA + OH
    k_f = 1.0000000000000002e-06 * 524000000000 
               * exp(0.69999999999999996 * tc[0] - 0.50321666580471969 * (5884.0799999999999) * invT);
    Corr  = 1.0;
    qf[138] *= Corr * k_f;
    qr[138] *= Corr * k_f / exp(g_RT[1] - g_RT[3] - g_RT[26] + g_RT[27]);
    // (139):  C3H6 + CH3 => C3H5XA + CH4
    k_f = 1.0000000000000002e-06 * 2.21 
               * exp(3.5 * tc[0] - 0.50321666580471969 * (5674.9499999999998) * invT);
    Corr  = 1.0;
    qf[139] *= Corr * k_f;
    qr[139] *= Corr * k_f / exp(g_RT[11] - g_RT[15] - g_RT[26] + g_RT[27]);
    // (140):  IXC3H7 + O2 => C3H6 + HO2
    k_f = 1.0000000000000002e-06 * 450000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (5020.0799999999999) * invT);
    Corr  = 1.0;
    qf[140] *= Corr * k_f;
    qr[140] *= Corr * k_f / exp(g_RT[5] - g_RT[8] - g_RT[27] + g_RT[29]);
    // (141):  IXC3H7 => H + C3H6
    k_f = 1 * 8.569e+18 
               * exp(-1.5700000000000001 * tc[0] - 0.50321666580471969 * (40340.110000000001) * invT);
    Corr  = 1.0;
    qf[141] *= Corr * k_f;
    qr[141] *= Corr * k_f / (exp(-g_RT[2] - g_RT[27] + g_RT[29]) * refC);
    // (142):  H + C3H6 => IXC3H7
    k_f = 1.0000000000000002e-06 * 13000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (1559.99) * invT);
    Corr  = 1.0;
    qf[142] *= Corr * k_f;
    qr[142] *= Corr * k_f / (exp(g_RT[2] + g_RT[27] - g_RT[29]) * refCinv);
    // (143):  NXC3H7 => CH3 + C2H4
    k_f = 1 * 228400000000000 
               * exp(-0.55000000000000004 * tc[0] - 0.50321666580471969 * (28400.099999999999) * invT);
    Corr  = 1.0;
    qf[143] *= Corr * k_f;
    qr[143] *= Corr * k_f / (exp(-g_RT[11] - g_RT[16] + g_RT[30]) * refC);
    // (144):  CH3 + C2H4 => NXC3H7
    k_f = 1.0000000000000002e-06 * 410000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (7204.1099999999997) * invT);
    Corr  = 1.0;
    qf[144] *= Corr * k_f;
    qr[144] *= Corr * k_f / (exp(g_RT[11] + g_RT[16] - g_RT[30]) * refCinv);
    // (145):  NXC3H7 + HO2 => C3H8 + O2
    k_f = 1.0000000000000002e-06 * 2080000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[145] *= Corr * k_f;
    qr[145] *= Corr * k_f / exp(-g_RT[5] + g_RT[8] + g_RT[30] - g_RT[31]);
    // (146):  NXC3H7 + O2 => C3H6 + HO2
    k_f = 1.0000000000000002e-06 * 300000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (3000) * invT);
    Corr  = 1.0;
    qf[146] *= Corr * k_f;
    qr[146] *= Corr * k_f / exp(g_RT[5] - g_RT[8] - g_RT[27] + g_RT[30]);
    // (147):  NXC3H7 => H + C3H6
    k_f = 1 * 2667000000000000 
               * exp(-0.64000000000000001 * tc[0] - 0.50321666580471969 * (36820.029999999999) * invT);
    Corr  = 1.0;
    qf[147] *= Corr * k_f;
    qr[147] *= Corr * k_f / (exp(-g_RT[2] - g_RT[27] + g_RT[30]) * refC);
    // (148):  H + C3H6 => NXC3H7
    k_f = 1.0000000000000002e-06 * 10000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (2500) * invT);
    Corr  = 1.0;
    qf[148] *= Corr * k_f;
    qr[148] *= Corr * k_f / (exp(g_RT[2] + g_RT[27] - g_RT[30]) * refCinv);
    // (149):  C3H8 + OH => NXC3H7 + H2O
    k_f = 1.0000000000000002e-06 * 10540000000 
               * exp(0.96999999999999997 * tc[0] - 0.50321666580471969 * (1586.04) * invT);
    Corr  = 1.0;
    qf[149] *= Corr * k_f;
    qr[149] *= Corr * k_f / exp(g_RT[3] - g_RT[6] - g_RT[30] + g_RT[31]);
    // (150):  C3H8 + HO2 => NXC3H7 + H2O2
    k_f = 1.0000000000000002e-06 * 16800000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (20429.970000000001) * invT);
    Corr  = 1.0;
    qf[150] *= Corr * k_f;
    qr[150] *= Corr * k_f / exp(-g_RT[7] + g_RT[8] - g_RT[30] + g_RT[31]);
    // (151):  H + C3H8 <=> H2 + NXC3H7
    k_f = 1.0000000000000002e-06 * 3972000 
               * exp(2.75 * tc[0] - 0.50321666580471969 * (6756.6899999999996) * invT);
    Corr  = 1.0;
    qf[151] *= Corr * k_f;
    qr[151] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[30] + g_RT[31]);
    // (152):  C3H8 + OH => IXC3H7 + H2O
    k_f = 1.0000000000000002e-06 * 46700000 
               * exp(1.6100000000000001 * tc[0] - 0.50321666580471969 * (-34.890000000000001) * invT);
    Corr  = 1.0;
    qf[152] *= Corr * k_f;
    qr[152] *= Corr * k_f / exp(g_RT[3] - g_RT[6] - g_RT[29] + g_RT[31]);
    // (153):  CH3 + C3H8 => CH4 + IXC3H7
    k_f = 1.0000000000000002e-06 * 398000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (9500) * invT);
    Corr  = 1.0;
    qf[153] *= Corr * k_f;
    qr[153] *= Corr * k_f / exp(g_RT[11] - g_RT[15] - g_RT[29] + g_RT[31]);
    // (154):  CH3 + C3H8 => CH4 + NXC3H7
    k_f = 1.0000000000000002e-06 * 1290000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (11599.9) * invT);
    Corr  = 1.0;
    qf[154] *= Corr * k_f;
    qr[154] *= Corr * k_f / exp(g_RT[11] - g_RT[15] - g_RT[30] + g_RT[31]);
    // (155):  C3H8 + O => IXC3H7 + OH
    k_f = 1.0000000000000002e-06 * 28100000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (5200.0500000000002) * invT);
    Corr  = 1.0;
    qf[155] *= Corr * k_f;
    qr[155] *= Corr * k_f / exp(g_RT[1] - g_RT[3] - g_RT[29] + g_RT[31]);
    // (156):  C3H8 + HO2 => IXC3H7 + H2O2
    k_f = 1.0000000000000002e-06 * 5600000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (17700.049999999999) * invT);
    Corr  = 1.0;
    qf[156] *= Corr * k_f;
    qr[156] *= Corr * k_f / exp(-g_RT[7] + g_RT[8] - g_RT[29] + g_RT[31]);
    // (157):  C3H8 + O => NXC3H7 + OH
    k_f = 1.0000000000000002e-06 * 113000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (7849.8999999999996) * invT);
    Corr  = 1.0;
    qf[157] *= Corr * k_f;
    qr[157] *= Corr * k_f / exp(g_RT[1] - g_RT[3] - g_RT[30] + g_RT[31]);
    // (158):  C3H8 + O2 => IXC3H7 + HO2
    k_f = 1.0000000000000002e-06 * 40000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (47500) * invT);
    Corr  = 1.0;
    qf[158] *= Corr * k_f;
    qr[158] *= Corr * k_f / exp(g_RT[5] - g_RT[8] - g_RT[29] + g_RT[31]);
    // (159):  IXC3H7 + HO2 => C3H8 + O2
    k_f = 1.0000000000000002e-06 * 2080000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[159] *= Corr * k_f;
    qr[159] *= Corr * k_f / exp(-g_RT[5] + g_RT[8] + g_RT[29] - g_RT[31]);
    // (160):  H + C3H8 => H2 + IXC3H7
    k_f = 1.0000000000000002e-06 * 1300000 
               * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (4471.0799999999999) * invT);
    Corr  = 1.0;
    qf[160] *= Corr * k_f;
    qr[160] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[29] + g_RT[31]);
    // (161):  H2 + IXC3H7 => H + C3H8
    k_f = 1.0000000000000002e-06 * 470900 
               * exp(2.1499999999999999 * tc[0] - 0.50321666580471969 * (12179.969999999999) * invT);
    Corr  = 1.0;
    qf[161] *= Corr * k_f;
    qr[161] *= Corr * k_f / exp(-g_RT[2] + g_RT[4] + g_RT[29] - g_RT[31]);
    // (162):  C3H5O => C2H3 + CH2O
    k_f = 1 * 2028000000000 
               * exp(0.089999999999999997 * tc[0] - 0.50321666580471969 * (23559.990000000002) * invT);
    Corr  = 1.0;
    qf[162] *= Corr * k_f;
    qr[162] *= Corr * k_f / (exp(-g_RT[10] - g_RT[24] + g_RT[28]) * refC);
    // (163):  IXC3H7O2 => IXC3H7 + O2
    k_f = 1 * 2.803e+17 
               * exp(-0.62 * tc[0] - 0.50321666580471969 * (36039.910000000003) * invT);
    Corr  = 1.0;
    qf[163] *= Corr * k_f;
    qr[163] *= Corr * k_f / (exp(-g_RT[5] - g_RT[29] + g_RT[32]) * refC);
    // (164):  IXC3H7 + O2 => IXC3H7O2
    k_f = 1.0000000000000002e-06 * 7540000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[164] *= Corr * k_f;
    qr[164] *= Corr * k_f / (exp(g_RT[5] + g_RT[29] - g_RT[32]) * refCinv);
    // (165):  NXC3H7O2 => NXC3H7 + O2
    k_f = 1 * 3.364e+19 
               * exp(-1.3200000000000001 * tc[0] - 0.50321666580471969 * (35760.040000000001) * invT);
    Corr  = 1.0;
    qf[165] *= Corr * k_f;
    qr[165] *= Corr * k_f / (exp(-g_RT[5] - g_RT[30] + g_RT[33]) * refC);
    // (166):  NXC3H7 + O2 => NXC3H7O2
    k_f = 1.0000000000000002e-06 * 4520000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[166] *= Corr * k_f;
    qr[166] *= Corr * k_f / (exp(g_RT[5] + g_RT[30] - g_RT[33]) * refCinv);


    return;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real * wdot, amrex::Real * sc, amrex::Real T)
{
    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];

    amrex::Real qdot, q_f[167], q_r[167];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 34; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[3] -= 2.000000 * qdot;
    wdot[7] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[2] -= qdot;
    wdot[5] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[2] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[11] -= 2.000000 * qdot;
    wdot[17] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[1] -= qdot;
    wdot[12] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[2] += qdot;
    wdot[10] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[2] += qdot;
    wdot[22] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[7]-q_r[7];
    wdot[2] -= qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[8]-q_r[8];
    wdot[4] += qdot;
    wdot[16] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[1] -= 2.000000 * qdot;
    wdot[5] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[2] += qdot;
    wdot[12] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[13]-q_r[13];
    wdot[2] -= qdot;
    wdot[12] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[15]-q_r[15];
    wdot[5] -= qdot;
    wdot[11] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[18]-q_r[18];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[19]-q_r[19];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[1] -= qdot;
    wdot[3] += 2.000000 * qdot;
    wdot[6] -= qdot;

    qdot = q_f[21]-q_r[21];
    wdot[1] += qdot;
    wdot[3] -= 2.000000 * qdot;
    wdot[6] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[24]-q_r[24];
    wdot[2] -= qdot;
    wdot[3] += 2.000000 * qdot;
    wdot[8] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[8] -= 2.000000 * qdot;

    qdot = q_f[26]-q_r[26];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[27]-q_r[27];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[6] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[29]-q_r[29];
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[8] -= 2.000000 * qdot;

    qdot = q_f[30]-q_r[30];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[31]-q_r[31];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[32]-q_r[32];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[34]-q_r[34];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[36]-q_r[36];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[37]-q_r[37];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[38]-q_r[38];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[39]-q_r[39];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[40]-q_r[40];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[9] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[41]-q_r[41];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[42]-q_r[42];
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[9] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[43]-q_r[43];
    wdot[1] -= qdot;
    wdot[2] += 2.000000 * qdot;
    wdot[9] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[44]-q_r[44];
    wdot[3] += qdot;
    wdot[8] -= qdot;
    wdot[11] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[45]-q_r[45];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[46]-q_r[46];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[11] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[47]-q_r[47];
    wdot[1] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[48]-q_r[48];
    wdot[2] += qdot;
    wdot[11] -= 2.000000 * qdot;
    wdot[14] += qdot;

    qdot = q_f[49]-q_r[49];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[50]-q_r[50];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[51]-q_r[51];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[52]-q_r[52];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[9] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[53]-q_r[53];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[9] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[54]-q_r[54];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[55]-q_r[55];
    wdot[2] += qdot;
    wdot[9] -= qdot;
    wdot[11] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[56]-q_r[56];
    wdot[9] -= qdot;
    wdot[11] += 2.000000 * qdot;
    wdot[15] -= qdot;

    qdot = q_f[57]-q_r[57];
    wdot[9] += qdot;
    wdot[11] -= 2.000000 * qdot;
    wdot[15] += qdot;

    qdot = q_f[58]-q_r[58];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[11] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[59]-q_r[59];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[60]-q_r[60];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[11] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[61]-q_r[61];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[62]-q_r[62];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[11] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[63]-q_r[63];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[64]-q_r[64];
    wdot[3] += qdot;
    wdot[8] -= qdot;
    wdot[12] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[65]-q_r[65];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[12] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[66]-q_r[66];
    wdot[1] -= qdot;
    wdot[5] += qdot;
    wdot[12] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[67]-q_r[67];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[12] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[68]-q_r[68];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[12] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[69]-q_r[69];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[15] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[70]-q_r[70];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[12] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[71]-q_r[71];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[12] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[72]-q_r[72];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[12] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[73]-q_r[73];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[74]-q_r[74];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[12] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[75]-q_r[75];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[10] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[76]-q_r[76];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[77]-q_r[77];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[10] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[78]-q_r[78];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[10] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[79]-q_r[79];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[10] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[80]-q_r[80];
    wdot[10] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;
    wdot[19] += qdot;

    qdot = q_f[81]-q_r[81];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[10] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[82]-q_r[82];
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[12] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[83]-q_r[83];
    wdot[11] -= qdot;
    wdot[13] += 2.000000 * qdot;
    wdot[20] -= qdot;

    qdot = q_f[84]-q_r[84];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[85]-q_r[85];
    wdot[5] += qdot;
    wdot[13] += 2.000000 * qdot;
    wdot[20] -= 2.000000 * qdot;

    qdot = q_f[86]-q_r[86];
    wdot[10] -= qdot;
    wdot[19] += qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[87]-q_r[87];
    wdot[3] += qdot;
    wdot[13] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[88]-q_r[88];
    wdot[3] -= qdot;
    wdot[13] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[89]-q_r[89];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[22] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[90]-q_r[90];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[22] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[91]-q_r[91];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[24] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[92]-q_r[92];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[22] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[93]-q_r[93];
    wdot[5] -= qdot;
    wdot[10] += qdot;
    wdot[19] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[94]-q_r[94];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[22] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[95]-q_r[95];
    wdot[1] -= qdot;
    wdot[11] += qdot;
    wdot[16] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[96]-q_r[96];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[16] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[97]-q_r[97];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[16] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[98]-q_r[98];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[16] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[99]-q_r[99];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[16] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[100]-q_r[100];
    wdot[16] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;
    wdot[24] += qdot;

    qdot = q_f[101]-q_r[101];
    wdot[11] -= qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[102]-q_r[102];
    wdot[11] += qdot;
    wdot[15] -= qdot;
    wdot[16] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[103]-q_r[103];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[16] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[104]-q_r[104];
    wdot[11] -= qdot;
    wdot[14] -= qdot;
    wdot[15] += qdot;
    wdot[16] += qdot;

    qdot = q_f[105]-q_r[105];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[14] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[106]-q_r[106];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[107]-q_r[107];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[14] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[108]-q_r[108];
    wdot[2] -= qdot;
    wdot[14] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[109]-q_r[109];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[14] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[110]-q_r[110];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[14] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[111]-q_r[111];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[14] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[112]-q_r[112];
    wdot[9] -= qdot;
    wdot[11] += qdot;
    wdot[14] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[113]-q_r[113];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[14] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[114]-q_r[114];
    wdot[11] -= qdot;
    wdot[14] += qdot;
    wdot[15] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[115]-q_r[115];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[12] += 2.000000 * qdot;
    wdot[23] -= qdot;

    qdot = q_f[116]-q_r[116];
    wdot[5] -= qdot;
    wdot[18] += qdot;
    wdot[19] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[117]-q_r[117];
    wdot[3] -= qdot;
    wdot[19] += 2.000000 * qdot;
    wdot[23] -= qdot;

    qdot = q_f[118]-q_r[118];
    wdot[2] -= qdot;
    wdot[9] += qdot;
    wdot[12] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[119]-q_r[119];
    wdot[2] += qdot;
    wdot[9] -= qdot;
    wdot[12] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[120]-q_r[120];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[10] += qdot;
    wdot[12] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[121]-q_r[121];
    wdot[10] -= qdot;
    wdot[19] += qdot;
    wdot[26] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[122]-q_r[122];
    wdot[10] += qdot;
    wdot[19] -= qdot;
    wdot[26] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[123]-q_r[123];
    wdot[3] += qdot;
    wdot[8] -= qdot;
    wdot[26] -= qdot;
    wdot[28] += qdot;

    qdot = q_f[124]-q_r[124];
    wdot[13] += qdot;
    wdot[20] -= qdot;
    wdot[26] -= qdot;
    wdot[28] += qdot;

    qdot = q_f[125]-q_r[125];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[26] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[126]-q_r[126];
    wdot[5] -= qdot;
    wdot[10] += qdot;
    wdot[25] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[127]-q_r[127];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[10] += qdot;
    wdot[22] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[128]-q_r[128];
    wdot[2] -= qdot;
    wdot[26] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[129]-q_r[129];
    wdot[11] += qdot;
    wdot[22] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[130]-q_r[130];
    wdot[11] += qdot;
    wdot[24] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[131]-q_r[131];
    wdot[11] -= qdot;
    wdot[24] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[132]-q_r[132];
    wdot[1] -= qdot;
    wdot[14] += qdot;
    wdot[19] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[133]-q_r[133];
    wdot[2] -= qdot;
    wdot[11] += qdot;
    wdot[16] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[134]-q_r[134];
    wdot[2] += qdot;
    wdot[11] -= qdot;
    wdot[16] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[135]-q_r[135];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[26] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[136]-q_r[136];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[26] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[137]-q_r[137];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[26] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[138]-q_r[138];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[26] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[139]-q_r[139];
    wdot[11] -= qdot;
    wdot[15] += qdot;
    wdot[26] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[140]-q_r[140];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[27] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[141]-q_r[141];
    wdot[2] += qdot;
    wdot[27] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[142]-q_r[142];
    wdot[2] -= qdot;
    wdot[27] -= qdot;
    wdot[29] += qdot;

    qdot = q_f[143]-q_r[143];
    wdot[11] += qdot;
    wdot[16] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[144]-q_r[144];
    wdot[11] -= qdot;
    wdot[16] -= qdot;
    wdot[30] += qdot;

    qdot = q_f[145]-q_r[145];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[30] -= qdot;
    wdot[31] += qdot;

    qdot = q_f[146]-q_r[146];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[27] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[147]-q_r[147];
    wdot[2] += qdot;
    wdot[27] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[148]-q_r[148];
    wdot[2] -= qdot;
    wdot[27] -= qdot;
    wdot[30] += qdot;

    qdot = q_f[149]-q_r[149];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[30] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[150]-q_r[150];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[30] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[151]-q_r[151];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[30] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[152]-q_r[152];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[29] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[153]-q_r[153];
    wdot[11] -= qdot;
    wdot[15] += qdot;
    wdot[29] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[154]-q_r[154];
    wdot[11] -= qdot;
    wdot[15] += qdot;
    wdot[30] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[155]-q_r[155];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[29] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[156]-q_r[156];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[29] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[157]-q_r[157];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[30] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[158]-q_r[158];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[29] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[159]-q_r[159];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[29] -= qdot;
    wdot[31] += qdot;

    qdot = q_f[160]-q_r[160];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[29] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[161]-q_r[161];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[29] -= qdot;
    wdot[31] += qdot;

    qdot = q_f[162]-q_r[162];
    wdot[10] += qdot;
    wdot[24] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[163]-q_r[163];
    wdot[5] += qdot;
    wdot[29] += qdot;
    wdot[32] -= qdot;

    qdot = q_f[164]-q_r[164];
    wdot[5] -= qdot;
    wdot[29] -= qdot;
    wdot[32] += qdot;

    qdot = q_f[165]-q_r[165];
    wdot[5] += qdot;
    wdot[30] += qdot;
    wdot[33] -= qdot;

    qdot = q_f[166]-q_r[166];
    wdot[5] -= qdot;
    wdot[30] -= qdot;
    wdot[33] += qdot;

    return;
}
#endif



/*compute the production rate for each species */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKWC(amrex::Real *  T, amrex::Real *  C,  amrex::Real *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 34; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 34; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKWYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[34]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[34];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*N2 */
    YOW += y[1]*imw[1]; /*O */
    YOW += y[2]*imw[2]; /*H */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2 */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*CH2GSG */
    YOW += y[10]*imw[10]; /*CH2O */
    YOW += y[11]*imw[11]; /*CH3 */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CH3O */
    YOW += y[14]*imw[14]; /*C2H5 */
    YOW += y[15]*imw[15]; /*CH4 */
    YOW += y[16]*imw[16]; /*C2H4 */
    YOW += y[17]*imw[17]; /*C2H6 */
    YOW += y[18]*imw[18]; /*CO2 */
    YOW += y[19]*imw[19]; /*HCO */
    YOW += y[20]*imw[20]; /*CH3O2 */
    YOW += y[21]*imw[21]; /*CH3O2H */
    YOW += y[22]*imw[22]; /*C2H2 */
    YOW += y[23]*imw[23]; /*HCCO */
    YOW += y[24]*imw[24]; /*C2H3 */
    YOW += y[25]*imw[25]; /*CH2CHO */
    YOW += y[26]*imw[26]; /*C3H5XA */
    YOW += y[27]*imw[27]; /*C3H6 */
    YOW += y[28]*imw[28]; /*C3H5O */
    YOW += y[29]*imw[29]; /*IXC3H7 */
    YOW += y[30]*imw[30]; /*NXC3H7 */
    YOW += y[31]*imw[31]; /*C3H8 */
    YOW += y[32]*imw[32]; /*IXC3H7O2 */
    YOW += y[33]*imw[33]; /*NXC3H7O2 */
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
    c[21] = PWORT * y[21]*imw[21]; 
    c[22] = PWORT * y[22]*imw[22]; 
    c[23] = PWORT * y[23]*imw[23]; 
    c[24] = PWORT * y[24]*imw[24]; 
    c[25] = PWORT * y[25]*imw[25]; 
    c[26] = PWORT * y[26]*imw[26]; 
    c[27] = PWORT * y[27]*imw[27]; 
    c[28] = PWORT * y[28]*imw[28]; 
    c[29] = PWORT * y[29]*imw[29]; 
    c[30] = PWORT * y[30]*imw[30]; 
    c[31] = PWORT * y[31]*imw[31]; 
    c[32] = PWORT * y[32]*imw[32]; 
    c[33] = PWORT * y[33]*imw[33]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 34; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKWXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[34]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 34; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 34; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKWYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[34]; /*temporary storage */
    amrex::Real imw[34];

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
    c[21] = 1e6 * (*rho) * y[21]*imw[21]; 
    c[22] = 1e6 * (*rho) * y[22]*imw[22]; 
    c[23] = 1e6 * (*rho) * y[23]*imw[23]; 
    c[24] = 1e6 * (*rho) * y[24]*imw[24]; 
    c[25] = 1e6 * (*rho) * y[25]*imw[25]; 
    c[26] = 1e6 * (*rho) * y[26]*imw[26]; 
    c[27] = 1e6 * (*rho) * y[27]*imw[27]; 
    c[28] = 1e6 * (*rho) * y[28]*imw[28]; 
    c[29] = 1e6 * (*rho) * y[29]*imw[29]; 
    c[30] = 1e6 * (*rho) * y[30]*imw[30]; 
    c[31] = 1e6 * (*rho) * y[31]*imw[31]; 
    c[32] = 1e6 * (*rho) * y[32]*imw[32]; 
    c[33] = 1e6 * (*rho) * y[33]*imw[33]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 34; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CKWXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[34]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*2.015940; /*H2 */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*14.027090; /*CH2GSG */
    XW += x[10]*30.026490; /*CH2O */
    XW += x[11]*15.035060; /*CH3 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*29.062150; /*C2H5 */
    XW += x[15]*16.043030; /*CH4 */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*30.070120; /*C2H6 */
    XW += x[18]*44.009950; /*CO2 */
    XW += x[19]*29.018520; /*HCO */
    XW += x[20]*47.033860; /*CH3O2 */
    XW += x[21]*48.041830; /*CH3O2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*41.029670; /*HCCO */
    XW += x[24]*27.046210; /*C2H3 */
    XW += x[25]*43.045610; /*CH2CHO */
    XW += x[26]*41.073300; /*C3H5XA */
    XW += x[27]*42.081270; /*C3H6 */
    XW += x[28]*57.072700; /*C3H5O */
    XW += x[29]*43.089240; /*IXC3H7 */
    XW += x[30]*43.089240; /*NXC3H7 */
    XW += x[31]*44.097210; /*C3H8 */
    XW += x[32]*75.088040; /*IXC3H7O2 */
    XW += x[33]*75.088040; /*NXC3H7O2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 34; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 34; ++id) {
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
        /*species 0: N2 */
        species[0] =
            +1.40824000e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77942000e-12 * tc[3];
        /*species 1: O */
        species[1] =
            -1.63816600e-03
            +4.84206400e-06 * tc[1]
            -4.80852900e-09 * tc[2]
            +1.55627840e-12 * tc[3];
        /*species 2: H */
        species[2] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 4: H2 */
        species[4] =
            +8.24944200e-04
            -1.62860300e-06 * tc[1]
            -2.84263020e-10 * tc[2]
            +1.65394880e-12 * tc[3];
        /*species 5: O2 */
        species[5] =
            +1.12748600e-03
            -1.15123000e-06 * tc[1]
            +3.94163100e-09 * tc[2]
            -3.50742160e-12 * tc[3];
        /*species 6: H2O */
        species[6] =
            +3.47498200e-03
            -1.27093920e-05 * tc[1]
            +2.09057430e-08 * tc[2]
            -1.00263520e-11 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +6.56922600e-03
            -2.97002600e-07 * tc[1]
            -1.38774180e-08 * tc[2]
            +9.88606000e-12 * tc[3];
        /*species 9: CH2GSG */
        species[9] =
            -1.69908900e-04
            +2.05073800e-06 * tc[1]
            +7.47765300e-09 * tc[2]
            -7.92506400e-12 * tc[3];
        /*species 10: CH2O */
        species[10] =
            +1.26314400e-02
            -3.77633600e-05 * tc[1]
            +6.15009300e-08 * tc[2]
            -3.36529480e-11 * tc[3];
        /*species 11: CH3 */
        species[11] =
            +1.11241000e-02
            -3.36044000e-05 * tc[1]
            +4.86548700e-08 * tc[2]
            -2.34598120e-11 * tc[3];
        /*species 12: CO */
        species[12] =
            +1.51194100e-03
            -7.76351000e-06 * tc[1]
            +1.67458320e-08 * tc[2]
            -9.89980400e-12 * tc[3];
        /*species 13: CH3O */
        species[13] =
            +7.21659500e-03
            +1.06769440e-05 * tc[1]
            -2.21329080e-08 * tc[2]
            +8.30244400e-12 * tc[3];
        /*species 14: C2H5 */
        species[14] =
            +8.71913300e-03
            +8.83967800e-06 * tc[1]
            +2.80161090e-09 * tc[2]
            -1.57110920e-11 * tc[3];
        /*species 15: CH4 */
        species[15] =
            +1.74766800e-02
            -5.56681800e-05 * tc[1]
            +9.14912400e-08 * tc[2]
            -4.89572400e-11 * tc[3];
        /*species 16: C2H4 */
        species[16] =
            +2.79616300e-02
            -6.77735400e-05 * tc[1]
            +8.35545600e-08 * tc[2]
            -3.89515160e-11 * tc[3];
        /*species 18: CO2 */
        species[18] =
            +9.92207200e-03
            -2.08182200e-05 * tc[1]
            +2.06000610e-08 * tc[2]
            -8.46912000e-12 * tc[3];
        /*species 19: HCO */
        species[19] =
            +6.19914700e-03
            -1.92461680e-05 * tc[1]
            +3.26947500e-08 * tc[2]
            -1.82995400e-11 * tc[3];
        /*species 22: C2H2 */
        species[22] =
            +1.51904500e-02
            -3.23263800e-05 * tc[1]
            +2.72369760e-08 * tc[2]
            -7.65098400e-12 * tc[3];
        /*species 23: HCCO */
        species[23] =
            +4.45347800e-03
            +4.53656600e-07 * tc[1]
            -4.44628500e-09 * tc[2]
            +9.00296800e-13 * tc[3];
        /*species 24: C2H3 */
        species[24] =
            +7.37147600e-03
            +4.21974600e-06 * tc[1]
            -3.96492600e-09 * tc[2]
            -4.73913600e-12 * tc[3];
        /*species 25: CH2CHO */
        species[25] =
            +1.07385700e-02
            +3.78298400e-06 * tc[1]
            -2.14757490e-08 * tc[2]
            +1.14695400e-11 * tc[3];
        /*species 29: IXC3H7 */
        species[29] =
            +2.54261600e-02
            +3.16161600e-06 * tc[1]
            -5.46385800e-08 * tc[2]
            +3.53108400e-11 * tc[3];
        /*species 30: NXC3H7 */
        species[30] =
            +2.47892700e-02
            +3.62049800e-06 * tc[1]
            -5.34979800e-08 * tc[2]
            +3.43319840e-11 * tc[3];
        /*species 31: C3H8 */
        species[31] =
            +2.66898600e-02
            +1.08628500e-05 * tc[1]
            -6.37800300e-08 * tc[2]
            +3.69733200e-11 * tc[3];
    } else {
        /*species 0: N2 */
        species[0] =
            +1.48797700e-03
            -1.13695220e-06 * tc[1]
            +3.02911200e-10 * tc[2]
            -2.70134040e-14 * tc[3];
        /*species 1: O */
        species[1] =
            -2.75506200e-05
            -6.20560600e-09 * tc[1]
            +1.36532010e-11 * tc[2]
            -1.74722080e-15 * tc[3];
        /*species 2: H */
        species[2] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 4: H2 */
        species[4] =
            +7.00064400e-04
            -1.12676580e-07 * tc[1]
            -2.76947340e-11 * tc[2]
            +6.33100800e-15 * tc[3];
        /*species 5: O2 */
        species[5] =
            +6.13519700e-04
            -2.51768400e-07 * tc[1]
            +5.32584300e-11 * tc[2]
            -4.54574000e-15 * tc[3];
        /*species 6: H2O */
        species[6] =
            +3.05629300e-03
            -1.74605200e-06 * tc[1]
            +3.60298800e-10 * tc[2]
            -2.55664720e-14 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.33613600e-03
            -2.94937800e-06 * tc[1]
            +7.04671200e-10 * tc[2]
            -5.72661600e-14 * tc[3];
        /*species 9: CH2GSG */
        species[9] =
            +2.06678800e-03
            -3.82823200e-07 * tc[1]
            -3.31401900e-10 * tc[2]
            +8.08540000e-14 * tc[3];
        /*species 10: CH2O */
        species[10] =
            +6.68132100e-03
            -5.25791000e-06 * tc[1]
            +1.42114590e-09 * tc[2]
            -1.28500680e-13 * tc[3];
        /*species 11: CH3 */
        species[11] =
            +6.13797400e-03
            -4.46069000e-06 * tc[1]
            +1.13554830e-09 * tc[2]
            -9.80863600e-14 * tc[3];
        /*species 12: CO */
        species[12] =
            +1.44268900e-03
            -1.12616560e-06 * tc[1]
            +3.05574300e-10 * tc[2]
            -2.76438080e-14 * tc[3];
        /*species 13: CH3O */
        species[13] =
            +7.87149700e-03
            -5.31276800e-06 * tc[1]
            +1.18332930e-09 * tc[2]
            -8.45046400e-14 * tc[3];
        /*species 14: C2H5 */
        species[14] =
            +6.48407700e-03
            -1.28561300e-06 * tc[1]
            -7.04363700e-10 * tc[2]
            +1.55235080e-13 * tc[3];
        /*species 15: CH4 */
        species[15] =
            +1.02372400e-02
            -7.75025800e-06 * tc[1]
            +2.03567550e-09 * tc[2]
            -1.80136920e-13 * tc[3];
        /*species 16: C2H4 */
        species[16] =
            +1.14851800e-02
            -8.83677000e-06 * tc[1]
            +2.35338030e-09 * tc[2]
            -2.10673920e-13 * tc[3];
        /*species 18: CO2 */
        species[18] =
            +3.14016900e-03
            -2.55682200e-06 * tc[1]
            +7.18199100e-10 * tc[2]
            -6.67613200e-14 * tc[3];
        /*species 19: HCO */
        species[19] =
            +3.34557300e-03
            -2.67001200e-06 * tc[1]
            +7.41171900e-10 * tc[2]
            -6.85540400e-14 * tc[3];
        /*species 22: C2H2 */
        species[22] =
            +5.37603900e-03
            -3.82563400e-06 * tc[1]
            +9.85913700e-10 * tc[2]
            -8.62684000e-14 * tc[3];
        /*species 23: HCCO */
        species[23] =
            +2.00040000e-03
            -4.05521400e-07 * tc[1]
            -3.12339600e-10 * tc[2]
            +7.86066000e-14 * tc[3];
        /*species 24: C2H3 */
        species[24] =
            +4.01774600e-03
            -7.93348000e-07 * tc[1]
            -4.32380100e-10 * tc[2]
            +9.51457600e-14 * tc[3];
        /*species 25: CH2CHO */
        species[25] =
            +8.13059100e-03
            -5.48724800e-06 * tc[1]
            +1.22109120e-09 * tc[2]
            -8.70406800e-14 * tc[3];
        /*species 29: IXC3H7 */
        species[29] =
            +1.57448800e-02
            -1.03647840e-05 * tc[1]
            +2.24317350e-09 * tc[2]
            -1.54176880e-13 * tc[3];
        /*species 30: NXC3H7 */
        species[30] =
            +1.57611300e-02
            -1.03464860e-05 * tc[1]
            +2.23316760e-09 * tc[2]
            -1.52999120e-13 * tc[3];
        /*species 31: C3H8 */
        species[31] =
            +1.88903400e-02
            -1.25678480e-05 * tc[1]
            +2.75381190e-09 * tc[2]
            -1.92496400e-13 * tc[3];
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 3: OH */
        species[3] =
            +2.02235804e-04
            -2.27092824e-07 * tc[1]
            +7.27335447e-10 * tc[2]
            -2.97460412e-13 * tc[3];
    } else {
        /*species 3: OH */
        species[3] =
            +1.31992406e-03
            -7.19449340e-07 * tc[1]
            +1.27689240e-10 * tc[2]
            -7.28192064e-15 * tc[3];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 8: HO2 */
        species[8] =
            +3.66767950e-03
            -1.86477024e-06 * tc[1]
            -9.77558757e-10 * tc[2]
            +6.04559648e-13 * tc[3];
        /*species 21: CH3O2H */
        species[21] =
            +1.90129767e-02
            -2.26772574e-05 * tc[1]
            +1.02091996e-08 * tc[2]
            -1.64732089e-12 * tc[3];
    } else {
        /*species 8: HO2 */
        species[8] =
            +2.38452835e-03
            -1.61269598e-06 * tc[1]
            +3.72575169e-10 * tc[2]
            -2.86560043e-14 * tc[3];
        /*species 21: CH3O2H */
        species[21] =
            +8.06817909e-03
            -5.54189842e-06 * tc[1]
            +1.29399673e-09 * tc[2]
            -1.00276858e-13 * tc[3];
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 17: C2H6 */
        species[17] =
            +2.40764754e-02
            -2.23786944e-05 * tc[1]
            +6.25022703e-09 * tc[2]
            -2.11947446e-13 * tc[3];
        /*species 33: NXC3H7O2 */
        species[33] =
            +3.96164986e-02
            -4.98983198e-05 * tc[1]
            +2.57835090e-08 * tc[2]
            -5.24961320e-12 * tc[3];
    } else {
        /*species 17: C2H6 */
        species[17] =
            +1.29236361e-02
            -8.85054392e-06 * tc[1]
            +2.06217518e-09 * tc[2]
            -1.59560693e-13 * tc[3];
        /*species 33: NXC3H7O2 */
        species[33] =
            +1.69910726e-02
            -1.17773375e-05 * tc[1]
            +2.76658619e-09 * tc[2]
            -2.15292270e-13 * tc[3];
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 20: CH3O2 */
        species[20] =
            +1.00873599e-02
            -6.43012368e-06 * tc[1]
            +6.28227801e-10 * tc[2]
            +1.67335641e-13 * tc[3];
    } else {
        /*species 20: CH3O2 */
        species[20] =
            +7.90728626e-03
            -5.36492468e-06 * tc[1]
            +1.24167401e-09 * tc[2]
            -9.56029320e-14 * tc[3];
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 26: C3H5XA */
        species[26] =
            +3.34559100e-02
            -5.06802054e-05 * tc[1]
            +3.08597262e-08 * tc[2]
            -6.93033360e-12 * tc[3];
    } else {
        /*species 26: C3H5XA */
        species[26] =
            +1.12695483e-02
            -7.67585728e-06 * tc[1]
            +1.78217736e-09 * tc[2]
            -1.37567212e-13 * tc[3];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 27: C3H6 */
        species[27] =
            +2.89107662e-02
            -3.09773616e-05 * tc[1]
            +1.16644263e-08 * tc[2]
            -1.35156141e-12 * tc[3];
        /*species 32: IXC3H7O2 */
        species[32] =
            +4.43081205e-02
            -6.44828912e-05 * tc[1]
            +3.89061408e-08 * tc[2]
            -8.93482276e-12 * tc[3];
    } else {
        /*species 27: C3H6 */
        species[27] =
            +1.37023634e-02
            -9.32499466e-06 * tc[1]
            +2.16376321e-09 * tc[2]
            -1.66948050e-13 * tc[3];
        /*species 32: IXC3H7O2 */
        species[32] =
            +1.64082190e-02
            -1.13486412e-05 * tc[1]
            +2.66200902e-09 * tc[2]
            -2.06944614e-13 * tc[3];
    }

    /*species with midpoint at T=1380 kelvin */
    if (T < 1380) {
        /*species 28: C3H5O */
        species[28] =
            +3.05579837e-02
            -3.61260552e-05 * tc[1]
            +1.45845010e-08 * tc[2]
            -1.67941825e-12 * tc[3];
    } else {
        /*species 28: C3H5O */
        species[28] =
            +1.14983720e-02
            -7.69291318e-06 * tc[1]
            +1.76673104e-09 * tc[2]
            -1.35423169e-13 * tc[3];
    }
    return;
}

/*compute an approx to the reaction Jacobian */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void aJacobian_precond(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int HP)
{
    for (int i=0; i<1225; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[34];
    for (int k=0; k<34; k++) {
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
    for (int k = 0; k < 34; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[34];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[34];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[34];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[3], 2.000000);
    k_f = 1.0000000000000002e-06 * 123600000000000
                * exp(-0.37 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.37 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* pressure-fall-off */
    k_0 = 3.0410000000000002e+30 * exp(-4.6299999999999999 * tc[0] - 0.50321666580471969 * (2049) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -4.6299999999999999 * invT + 0.50321666580471969 * (2049) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.46999999999999997)*exp(-T/100);
    Fcent2 = 0.46999999999999997 * exp(-T/2000);
    Fcent3 = exp(-1000000000000000 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/100
        -Fcent2/2000
        + Fcent3*1000000000000000*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[7];
    Kc = refCinv * exp(2.000000*g_RT[3] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[3]) + (h_RT[7]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[3] -= 2 * q; /* OH */
    wdot[7] += q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac + k_f*2.000000*sc[3];
    dqdc[4] = 2.5*dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = 12*dcdc_fac;
    dqdc[7] = dcdc_fac - k_r;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    dqdc[12] = 1.8999999999999999*dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = 3.7999999999999998*dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    dqdc[26] = dcdc_fac;
    dqdc[27] = dcdc_fac;
    dqdc[28] = dcdc_fac;
    dqdc[29] = dcdc_fac;
    dqdc[30] = dcdc_fac;
    dqdc[31] = dcdc_fac;
    dqdc[32] = dcdc_fac;
    dqdc[33] = dcdc_fac;
    for (int k=0; k<34; k++) {
        J[35*k+3] += -2 * dqdc[k];
        J[35*k+7] += dqdc[k];
    }
    J[1193] += -2 * dqdT; /* dwdot[OH]/dT */
    J[1197] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 2: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = sc[2]*sc[5];
    k_f = 1.0000000000000002e-06 * 1475000000000
                * exp(0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0.59999999999999998 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* pressure-fall-off */
    k_0 = 35000000000000000 * exp(-0.40999999999999998 * tc[0] - 0.50321666580471969 * (-1115.9200000000001) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -0.40999999999999998 * invT + 0.50321666580471969 * (-1115.9200000000001) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = 0.;
    Fcent2 = 0.5 * exp(-T/1e+30);
    Fcent3 = exp(-1e+100 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        +0.
        -Fcent2/1e+30
        + Fcent3*1e+100*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[8];
    Kc = refCinv * exp(g_RT[2] + g_RT[5] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[8]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac + k_f*sc[5];
    dqdc[3] = dcdc_fac;
    dqdc[4] = 2.5*dcdc_fac;
    dqdc[5] = dcdc_fac + k_f*sc[2];
    dqdc[6] = 12*dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac - k_r;
    dqdc[9] = dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    dqdc[12] = 1.8999999999999999*dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = 3.7999999999999998*dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    dqdc[26] = dcdc_fac;
    dqdc[27] = dcdc_fac;
    dqdc[28] = dcdc_fac;
    dqdc[29] = dcdc_fac;
    dqdc[30] = dcdc_fac;
    dqdc[31] = dcdc_fac;
    dqdc[32] = dcdc_fac;
    dqdc[33] = dcdc_fac;
    for (int k=0; k<34; k++) {
        J[35*k+2] -= dqdc[k];
        J[35*k+5] -= dqdc[k];
        J[35*k+8] += dqdc[k];
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1195] -= dqdT; /* dwdot[O2]/dT */
    J[1198] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 3: CH3 + H (+M) => CH4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[4] + ( 5 - 1)*sc[6] + ( 3 - 1)*sc[18] + ( 2 - 1)*sc[12];
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = 1.0000000000000002e-06 * 2138000000000000
                * exp(-0.40000000000000002 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.40000000000000002 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* pressure-fall-off */
    k_0 = 3.31e+30 * exp(-4 * tc[0] - 0.50321666580471969 * (2108.0300000000002) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -4 * invT + 0.50321666580471969 * (2108.0300000000002) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = 0.;
    Fcent2 = 0.;
    Fcent3 = exp(-40 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        +0.
        +0.
        + Fcent3*40*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac + k_f*sc[11];
    dqdc[3] = dcdc_fac;
    dqdc[4] = 2*dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = 5*dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac + k_f*sc[2];
    dqdc[12] = 2*dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = 3*dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    dqdc[26] = dcdc_fac;
    dqdc[27] = dcdc_fac;
    dqdc[28] = dcdc_fac;
    dqdc[29] = dcdc_fac;
    dqdc[30] = dcdc_fac;
    dqdc[31] = dcdc_fac;
    dqdc[32] = dcdc_fac;
    dqdc[33] = dcdc_fac;
    for (int k=0; k<34; k++) {
        J[35*k+2] -= dqdc[k];
        J[35*k+11] -= dqdc[k];
        J[35*k+15] += dqdc[k];
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1201] -= dqdT; /* dwdot[CH3]/dT */
    J[1205] += dqdT; /* dwdot[CH4]/dT */

    /*reaction 4: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[4] + ( 5 - 1)*sc[6] + ( 3 - 1)*sc[18] + ( 2 - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[11], 2.000000);
    k_f = 1.0000000000000002e-06 * 92140000000000000
                * exp(-1.1699999999999999 * tc[0] - 0.50321666580471969 * (635.75999999999999) * invT);
    dlnkfdT = -1.1699999999999999 * invT + 0.50321666580471969 *  (635.75999999999999)  * invT2;
    /* pressure-fall-off */
    k_0 = 1.135e+36 * exp(-5.2460000000000004 * tc[0] - 0.50321666580471969 * (1705.0699999999999) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -5.2460000000000004 * invT + 0.50321666580471969 * (1705.0699999999999) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.40500000000000003)*exp(-T/1120);
    Fcent2 = 0.40500000000000003 * exp(-T/69.599999999999994);
    Fcent3 = exp(-1000000000000000 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1120
        -Fcent2/69.599999999999994
        + Fcent3*1000000000000000*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[17];
    Kc = refCinv * exp(2.000000*g_RT[11] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[11]) + (h_RT[17]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[11] -= 2 * q; /* CH3 */
    wdot[17] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = 2*dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = 5*dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac + k_f*2.000000*sc[11];
    dqdc[12] = 2*dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac - k_r;
    dqdc[18] = 3*dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    dqdc[26] = dcdc_fac;
    dqdc[27] = dcdc_fac;
    dqdc[28] = dcdc_fac;
    dqdc[29] = dcdc_fac;
    dqdc[30] = dcdc_fac;
    dqdc[31] = dcdc_fac;
    dqdc[32] = dcdc_fac;
    dqdc[33] = dcdc_fac;
    for (int k=0; k<34; k++) {
        J[35*k+11] += -2 * dqdc[k];
        J[35*k+17] += dqdc[k];
    }
    J[1201] += -2 * dqdT; /* dwdot[CH3]/dT */
    J[1207] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 5: CO + O (+M) => CO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = sc[1]*sc[12];
    k_f = 1.0000000000000002e-06 * 18000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (2384.0799999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (2384.0799999999999)  * invT2;
    /* pressure-fall-off */
    k_0 = 1.35e+24 * exp(-2.7879999999999998 * tc[0] - 0.50321666580471969 * (4190.9700000000003) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.7879999999999998 * invT + 0.50321666580471969 * (4190.9700000000003) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-1)*exp(-T/1);
    Fcent2 = 1 * exp(-T/10000000);
    Fcent3 = exp(-10000000 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1
        -Fcent2/10000000
        + Fcent3*10000000*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[12] -= q; /* CO */
    wdot[18] += q; /* CO2 */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac;
    dqdc[1] = dcdc_fac + k_f*sc[12];
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = 2.5*dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = 12*dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    dqdc[12] = 1.8999999999999999*dcdc_fac + k_f*sc[1];
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = 3.7999999999999998*dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    dqdc[26] = dcdc_fac;
    dqdc[27] = dcdc_fac;
    dqdc[28] = dcdc_fac;
    dqdc[29] = dcdc_fac;
    dqdc[30] = dcdc_fac;
    dqdc[31] = dcdc_fac;
    dqdc[32] = dcdc_fac;
    dqdc[33] = dcdc_fac;
    for (int k=0; k<34; k++) {
        J[35*k+1] -= dqdc[k];
        J[35*k+12] -= dqdc[k];
        J[35*k+18] += dqdc[k];
    }
    J[1191] -= dqdT; /* dwdot[O]/dT */
    J[1202] -= dqdT; /* dwdot[CO]/dT */
    J[1208] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 6: CH3O (+M) => CH2O + H (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[13];
    k_f = 1 * 54500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (13500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (13500)  * invT2;
    /* pressure-fall-off */
    k_0 = 2.344e+25 * exp(-2.7000000000000002 * tc[0] - 0.50321666580471969 * (30599.900000000001) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.7000000000000002 * invT + 0.50321666580471969 * (30599.900000000001) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-1)*exp(-T/1);
    Fcent2 = 1 * exp(-T/10000000);
    Fcent3 = exp(-10000000 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1
        -Fcent2/10000000
        + Fcent3*10000000*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[10] += q; /* CH2O */
    wdot[13] -= q; /* CH3O */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    dqdc[12] = dcdc_fac;
    dqdc[13] = dcdc_fac + k_f;
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
    dqdc[26] = dcdc_fac;
    dqdc[27] = dcdc_fac;
    dqdc[28] = dcdc_fac;
    dqdc[29] = dcdc_fac;
    dqdc[30] = dcdc_fac;
    dqdc[31] = dcdc_fac;
    dqdc[32] = dcdc_fac;
    dqdc[33] = dcdc_fac;
    for (int k=0; k<34; k++) {
        J[35*k+2] += dqdc[k];
        J[35*k+10] += dqdc[k];
        J[35*k+13] -= dqdc[k];
    }
    J[1192] += dqdT; /* dwdot[H]/dT */
    J[1200] += dqdT; /* dwdot[CH2O]/dT */
    J[1203] -= dqdT; /* dwdot[CH3O]/dT */

    /*reaction 7: C2H3 (+M) => H + C2H2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[4] + ( 5 - 1)*sc[6] + ( 3 - 1)*sc[18] + ( 2 - 1)*sc[12];
    /* forward */
    phi_f = sc[24];
    k_f = 1 * 16060000000
                * exp(1.028 * tc[0] - 0.50321666580471969 * (40503.589999999997) * invT);
    dlnkfdT = 1.028 * invT + 0.50321666580471969 *  (40503.589999999997)  * invT2;
    /* pressure-fall-off */
    k_0 = 1.164e+39 * exp(-6.8209999999999997 * tc[0] - 0.50321666580471969 * (44491.629999999997) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -6.8209999999999997 * invT + 0.50321666580471969 * (44491.629999999997) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = 0.;
    Fcent2 = 1 * exp(-T/675);
    Fcent3 = exp(-1000000000000000 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        +0.
        -Fcent2/675
        + Fcent3*1000000000000000*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[22] += q; /* C2H2 */
    wdot[24] -= q; /* C2H3 */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = 2*dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = 5*dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    dqdc[12] = 2*dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac;
    dqdc[17] = dcdc_fac;
    dqdc[18] = 3*dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac + k_f;
    dqdc[25] = dcdc_fac;
    dqdc[26] = dcdc_fac;
    dqdc[27] = dcdc_fac;
    dqdc[28] = dcdc_fac;
    dqdc[29] = dcdc_fac;
    dqdc[30] = dcdc_fac;
    dqdc[31] = dcdc_fac;
    dqdc[32] = dcdc_fac;
    dqdc[33] = dcdc_fac;
    for (int k=0; k<34; k++) {
        J[35*k+2] += dqdc[k];
        J[35*k+22] += dqdc[k];
        J[35*k+24] -= dqdc[k];
    }
    J[1192] += dqdT; /* dwdot[H]/dT */
    J[1212] += dqdT; /* dwdot[C2H2]/dT */
    J[1214] -= dqdT; /* dwdot[C2H3]/dT */

    /*reaction 8: H + C2H4 (+M) <=> C2H5 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[2]*sc[16];
    k_f = 1.0000000000000002e-06 * 1081000000000
                * exp(0.45000000000000001 * tc[0] - 0.50321666580471969 * (1821.9400000000001) * invT);
    dlnkfdT = 0.45000000000000001 * invT + 0.50321666580471969 *  (1821.9400000000001)  * invT2;
    /* pressure-fall-off */
    k_0 = 1.1120000000000001e+34 * exp(-5 * tc[0] - 0.50321666580471969 * (4447.8999999999996) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -5 * invT + 0.50321666580471969 * (4447.8999999999996) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = 0.;
    Fcent2 = 1 * exp(-T/95);
    Fcent3 = exp(-200 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        +0.
        -Fcent2/95
        + Fcent3*200*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[14];
    Kc = refCinv * exp(g_RT[2] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[16]) + (h_RT[14]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[14] += q; /* C2H5 */
    wdot[16] -= q; /* C2H4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac + k_f*sc[16];
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    dqdc[12] = dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac - k_r;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac + k_f*sc[2];
    dqdc[17] = dcdc_fac;
    dqdc[18] = dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    dqdc[26] = dcdc_fac;
    dqdc[27] = dcdc_fac;
    dqdc[28] = dcdc_fac;
    dqdc[29] = dcdc_fac;
    dqdc[30] = dcdc_fac;
    dqdc[31] = dcdc_fac;
    dqdc[32] = dcdc_fac;
    dqdc[33] = dcdc_fac;
    for (int k=0; k<34; k++) {
        J[35*k+2] -= dqdc[k];
        J[35*k+14] += dqdc[k];
        J[35*k+16] -= dqdc[k];
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1204] += dqdT; /* dwdot[C2H5]/dT */
    J[1206] -= dqdT; /* dwdot[C2H4]/dT */

    /*reaction 9: C2H4 (+M) => C2H2 + H2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[16];
    k_f = 1 * 18000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (76000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (76000)  * invT2;
    /* pressure-fall-off */
    k_0 = 1500000000000000 * exp(0 * tc[0] - 0.50321666580471969 * (55440.010000000002) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = 0 * invT + 0.50321666580471969 * (55440.010000000002) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-1)*exp(-T/1);
    Fcent2 = 1 * exp(-T/10000000);
    Fcent3 = exp(-10000000 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1
        -Fcent2/10000000
        + Fcent3*10000000*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[4] += q; /* H2 */
    wdot[16] -= q; /* C2H4 */
    wdot[22] += q; /* C2H2 */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    dqdc[12] = dcdc_fac;
    dqdc[13] = dcdc_fac;
    dqdc[14] = dcdc_fac;
    dqdc[15] = dcdc_fac;
    dqdc[16] = dcdc_fac + k_f;
    dqdc[17] = dcdc_fac;
    dqdc[18] = dcdc_fac;
    dqdc[19] = dcdc_fac;
    dqdc[20] = dcdc_fac;
    dqdc[21] = dcdc_fac;
    dqdc[22] = dcdc_fac;
    dqdc[23] = dcdc_fac;
    dqdc[24] = dcdc_fac;
    dqdc[25] = dcdc_fac;
    dqdc[26] = dcdc_fac;
    dqdc[27] = dcdc_fac;
    dqdc[28] = dcdc_fac;
    dqdc[29] = dcdc_fac;
    dqdc[30] = dcdc_fac;
    dqdc[31] = dcdc_fac;
    dqdc[32] = dcdc_fac;
    dqdc[33] = dcdc_fac;
    for (int k=0; k<34; k++) {
        J[35*k+4] += dqdc[k];
        J[35*k+16] -= dqdc[k];
        J[35*k+22] += dqdc[k];
    }
    J[1194] += dqdT; /* dwdot[H2]/dT */
    J[1206] -= dqdT; /* dwdot[C2H4]/dT */
    J[1212] += dqdT; /* dwdot[C2H2]/dT */

    /*reaction 10: O + H + M => OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-12 * 4.72e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    dqdc[0] = q_nocor;
    dqdc[1] = q_nocor + k_f*sc[2];
    dqdc[2] = q_nocor + k_f*sc[1];
    dqdc[3] = q_nocor;
    dqdc[4] = 2.5*q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = 12*q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    dqdc[12] = 1.8999999999999999*q_nocor;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = 3.7999999999999998*q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = q_nocor;
    dqdc[21] = q_nocor;
    dqdc[22] = q_nocor;
    dqdc[23] = q_nocor;
    dqdc[24] = q_nocor;
    dqdc[25] = q_nocor;
    dqdc[26] = q_nocor;
    dqdc[27] = q_nocor;
    dqdc[28] = q_nocor;
    dqdc[29] = q_nocor;
    dqdc[30] = q_nocor;
    dqdc[31] = q_nocor;
    dqdc[32] = q_nocor;
    dqdc[33] = q_nocor;
    for (int k=0; k<34; k++) {
        J[35*k+1] -= dqdc[k];
        J[35*k+2] -= dqdc[k];
        J[35*k+3] += dqdc[k];
    }
    J[1191] -= dqdT; /* dwdot[O]/dT */
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1193] += dqdT; /* dwdot[OH]/dT */

    /*reaction 11: 2.000000 O + M => O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 6170000000000000
                * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.5 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    dqdc[0] = q_nocor;
    dqdc[1] = q_nocor + k_f*2.000000*sc[1];
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = 2.5*q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = 12*q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    dqdc[12] = 1.8999999999999999*q_nocor;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = 3.7999999999999998*q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = q_nocor;
    dqdc[21] = q_nocor;
    dqdc[22] = q_nocor;
    dqdc[23] = q_nocor;
    dqdc[24] = q_nocor;
    dqdc[25] = q_nocor;
    dqdc[26] = q_nocor;
    dqdc[27] = q_nocor;
    dqdc[28] = q_nocor;
    dqdc[29] = q_nocor;
    dqdc[30] = q_nocor;
    dqdc[31] = q_nocor;
    dqdc[32] = q_nocor;
    dqdc[33] = q_nocor;
    for (int k=0; k<34; k++) {
        J[35*k+1] += -2 * dqdc[k];
        J[35*k+5] += dqdc[k];
    }
    J[1191] += -2 * dqdT; /* dwdot[O]/dT */
    J[1195] += dqdT; /* dwdot[O2]/dT */

    /*reaction 12: H + OH + M => H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = sc[2]*sc[3];
    k_f = 1.0000000000000002e-12 * 2.2499999999999999e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    dqdc[0] = q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor + k_f*sc[3];
    dqdc[3] = q_nocor + k_f*sc[2];
    dqdc[4] = 2.5*q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = 12*q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    dqdc[12] = 1.8999999999999999*q_nocor;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = 3.7999999999999998*q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = q_nocor;
    dqdc[21] = q_nocor;
    dqdc[22] = q_nocor;
    dqdc[23] = q_nocor;
    dqdc[24] = q_nocor;
    dqdc[25] = q_nocor;
    dqdc[26] = q_nocor;
    dqdc[27] = q_nocor;
    dqdc[28] = q_nocor;
    dqdc[29] = q_nocor;
    dqdc[30] = q_nocor;
    dqdc[31] = q_nocor;
    dqdc[32] = q_nocor;
    dqdc[33] = q_nocor;
    for (int k=0; k<34; k++) {
        J[35*k+2] -= dqdc[k];
        J[35*k+3] -= dqdc[k];
        J[35*k+6] += dqdc[k];
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1193] -= dqdT; /* dwdot[OH]/dT */
    J[1196] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 13: HCO + M => H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 6 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = sc[19];
    k_f = 1.0000000000000002e-06 * 1.86e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (17000)  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    dqdc[0] = q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = 2.5*q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = 6*q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    dqdc[12] = 1.8999999999999999*q_nocor;
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = 3.7999999999999998*q_nocor;
    dqdc[19] = q_nocor + k_f;
    dqdc[20] = q_nocor;
    dqdc[21] = q_nocor;
    dqdc[22] = q_nocor;
    dqdc[23] = q_nocor;
    dqdc[24] = q_nocor;
    dqdc[25] = q_nocor;
    dqdc[26] = q_nocor;
    dqdc[27] = q_nocor;
    dqdc[28] = q_nocor;
    dqdc[29] = q_nocor;
    dqdc[30] = q_nocor;
    dqdc[31] = q_nocor;
    dqdc[32] = q_nocor;
    dqdc[33] = q_nocor;
    for (int k=0; k<34; k++) {
        J[35*k+2] += dqdc[k];
        J[35*k+12] += dqdc[k];
        J[35*k+19] -= dqdc[k];
    }
    J[1192] += dqdT; /* dwdot[H]/dT */
    J[1202] += dqdT; /* dwdot[CO]/dT */
    J[1209] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 14: H + CO + M => HCO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 6 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = sc[2]*sc[12];
    k_f = 1.0000000000000002e-12 * 64670000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-441.92000000000002) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-441.92000000000002)  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[12] -= q; /* CO */
    wdot[19] += q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    dqdc[0] = q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor + k_f*sc[12];
    dqdc[3] = q_nocor;
    dqdc[4] = 2.5*q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = 6*q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    dqdc[12] = 1.8999999999999999*q_nocor + k_f*sc[2];
    dqdc[13] = q_nocor;
    dqdc[14] = q_nocor;
    dqdc[15] = q_nocor;
    dqdc[16] = q_nocor;
    dqdc[17] = q_nocor;
    dqdc[18] = 3.7999999999999998*q_nocor;
    dqdc[19] = q_nocor;
    dqdc[20] = q_nocor;
    dqdc[21] = q_nocor;
    dqdc[22] = q_nocor;
    dqdc[23] = q_nocor;
    dqdc[24] = q_nocor;
    dqdc[25] = q_nocor;
    dqdc[26] = q_nocor;
    dqdc[27] = q_nocor;
    dqdc[28] = q_nocor;
    dqdc[29] = q_nocor;
    dqdc[30] = q_nocor;
    dqdc[31] = q_nocor;
    dqdc[32] = q_nocor;
    dqdc[33] = q_nocor;
    for (int k=0; k<34; k++) {
        J[35*k+2] -= dqdc[k];
        J[35*k+12] -= dqdc[k];
        J[35*k+19] += dqdc[k];
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1202] -= dqdT; /* dwdot[CO]/dT */
    J[1209] += dqdT; /* dwdot[HCO]/dT */

    /*reaction 15: CH3O2 + M => CH3 + O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[20];
    k_f = 1.0000000000000002e-06 * 4.3430000000000002e+27
                * exp(-3.4199999999999999 * tc[0] - 0.50321666580471969 * (30469.889999999999) * invT);
    dlnkfdT = -3.4199999999999999 * invT + 0.50321666580471969 *  (30469.889999999999)  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[11] += q; /* CH3 */
    wdot[20] -= q; /* CH3O2 */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    dqdc[0] = q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
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
    dqdc[20] = q_nocor + k_f;
    dqdc[21] = q_nocor;
    dqdc[22] = q_nocor;
    dqdc[23] = q_nocor;
    dqdc[24] = q_nocor;
    dqdc[25] = q_nocor;
    dqdc[26] = q_nocor;
    dqdc[27] = q_nocor;
    dqdc[28] = q_nocor;
    dqdc[29] = q_nocor;
    dqdc[30] = q_nocor;
    dqdc[31] = q_nocor;
    dqdc[32] = q_nocor;
    dqdc[33] = q_nocor;
    for (int k=0; k<34; k++) {
        J[35*k+5] += dqdc[k];
        J[35*k+11] += dqdc[k];
        J[35*k+20] -= dqdc[k];
    }
    J[1195] += dqdT; /* dwdot[O2]/dT */
    J[1201] += dqdT; /* dwdot[CH3]/dT */
    J[1210] -= dqdT; /* dwdot[CH3O2]/dT */

    /*reaction 16: CH3 + O2 + M => CH3O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = 1.0000000000000002e-12 * 5.4400000000000001e+25
                * exp(-3.2999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -3.2999999999999998 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[11] -= q; /* CH3 */
    wdot[20] += q; /* CH3O2 */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    dqdc[0] = q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor;
    dqdc[5] = q_nocor + k_f*sc[11];
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor + k_f*sc[5];
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
    dqdc[26] = q_nocor;
    dqdc[27] = q_nocor;
    dqdc[28] = q_nocor;
    dqdc[29] = q_nocor;
    dqdc[30] = q_nocor;
    dqdc[31] = q_nocor;
    dqdc[32] = q_nocor;
    dqdc[33] = q_nocor;
    for (int k=0; k<34; k++) {
        J[35*k+5] -= dqdc[k];
        J[35*k+11] -= dqdc[k];
        J[35*k+20] += dqdc[k];
    }
    J[1195] -= dqdT; /* dwdot[O2]/dT */
    J[1201] -= dqdT; /* dwdot[CH3]/dT */
    J[1210] += dqdT; /* dwdot[CH3O2]/dT */

    /*reaction 17: O + H2 => H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = 1.0000000000000002e-06 * 50800
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6292.0699999999997) * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  (6292.0699999999997)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[4];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[39] -= dqdci;               /* dwdot[H2]/d[O] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[1];
    J[141] -= dqdci;              /* dwdot[O]/d[H2] */
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[143] += dqdci;              /* dwdot[OH]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */

    /*reaction 18: H + OH => O + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[3];
    k_f = 1.0000000000000002e-06 * 22310
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (4196.9399999999996) * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  (4196.9399999999996)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[2] -= q; /* H */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3];
    J[71] += dqdci;               /* dwdot[O]/d[H] */
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] -= dqdci;               /* dwdot[OH]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[2];
    J[106] += dqdci;              /* dwdot[O]/d[OH] */
    J[107] -= dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[109] += dqdci;              /* dwdot[H2]/d[OH] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */

    /*reaction 19: OH + H2 => H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[4];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3429.9699999999998) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  (3429.9699999999998)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[4] -= q; /* H2 */
    wdot[6] += q; /* H2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[4];
    J[107] += dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[109] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[143] -= dqdci;              /* dwdot[OH]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[146] += dqdci;              /* dwdot[H2O]/d[H2] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */

    /*reaction 20: H + H2O => OH + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = 1.0000000000000002e-06 * 935200000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (18580.07) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  (18580.07)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] += q; /* H2 */
    wdot[6] -= q; /* H2O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[OH]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[76] -= dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[2];
    J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[214] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */

    /*reaction 21: O + H2O => 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 2970000
                * exp(2.02 * tc[0] - 0.50321666580471969 * (13400.1) * invT);
    dlnkfdT = 2.02 * invT + 0.50321666580471969 *  (13400.1)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += 2 * q; /* OH */
    wdot[6] -= q; /* H2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += 2 * dqdci;           /* dwdot[OH]/d[O] */
    J[41] -= dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1];
    J[211] -= dqdci;              /* dwdot[O]/d[H2O] */
    J[213] += 2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += 2 * dqdT;          /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */

    /*reaction 22: 2.000000 OH => O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[3], 2.000000);
    k_f = 1.0000000000000002e-06 * 301300
                * exp(2.02 * tc[0] - 0.50321666580471969 * (-3849.9000000000001) * invT);
    dlnkfdT = 2.02 * invT + 0.50321666580471969 *  (-3849.9000000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[3] -= 2 * q; /* OH */
    wdot[6] += q; /* H2O */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[3];
    J[106] += dqdci;              /* dwdot[O]/d[OH] */
    J[108] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1193] += -2 * dqdT;         /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */

    /*reaction 23: H + O2 => O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[5];
    k_f = 1.0000000000000002e-06 * 197000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (16539.91) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (16539.91)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[5];
    J[71] += dqdci;               /* dwdot[O]/d[H] */
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[OH]/d[H] */
    J[75] -= dqdci;               /* dwdot[O2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[2];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[177] -= dqdci;              /* dwdot[H]/d[O2] */
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */

    /*reaction 24: O + OH => H + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 15550000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (424.94999999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (424.94999999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[3];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[38] -= dqdci;               /* dwdot[OH]/d[O] */
    J[40] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[106] -= dqdci;              /* dwdot[O]/d[OH] */
    J[107] += dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[110] += dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */

    /*reaction 25: HO2 + H => 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = 1.0000000000000002e-06 * 70800000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (299.94999999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (299.94999999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += 2 * q; /* OH */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[78] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[282] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[283] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += 2 * dqdT;          /* dwdot[OH]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 26: 2.000000 HO2 => H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[8], 2.000000);
    k_f = 1.0000000000000002e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (11979.92) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (11979.92)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= 2 * q; /* HO2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[8];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += -2 * dqdT;         /* dwdot[HO2]/dT */

    /*reaction 27: HO2 + H => H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = 1.0000000000000002e-06 * 16600000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (820.02999999999997) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (820.02999999999997)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[75] += dqdci;               /* dwdot[O2]/d[H] */
    J[78] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[282] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[284] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 28: HO2 + OH => H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = 1.0000000000000002e-06 * 28900000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-500)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[6] += q; /* H2O */
    wdot[8] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[110] += dqdci;              /* dwdot[O2]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[113] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[283] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[286] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 29: H2O + O2 => HO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = 1.0000000000000002e-06 * 6888000000000000
                * exp(-0.33000000000000002 * tc[0] - 0.50321666580471969 * (72140.059999999998) * invT);
    dlnkfdT = -0.33000000000000002 * invT + 0.50321666580471969 *  (72140.059999999998)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[6] -= q; /* H2O */
    wdot[8] += q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[6];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[181] -= dqdci;              /* dwdot[H2O]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[5];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[215] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[218] += dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 30: 2.000000 HO2 => H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[8], 2.000000);
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1629.0599999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1629.0599999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= 2 * q; /* HO2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[8];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += -2 * dqdT;         /* dwdot[HO2]/dT */

    /*reaction 31: HO2 + O => OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[8];
    k_f = 1.0000000000000002e-06 * 32500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[8];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[40] += dqdci;               /* dwdot[O2]/d[O] */
    J[43] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[281] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 32: OH + O2 => HO2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[5];
    k_f = 1.0000000000000002e-06 * 785700000000000
                * exp(-0.33000000000000002 * tc[0] - 0.50321666580471969 * (55390.059999999998) * invT);
    dlnkfdT = -0.33000000000000002 * invT + 0.50321666580471969 *  (55390.059999999998)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[3] -= q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[5];
    J[106] += dqdci;              /* dwdot[O]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[110] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[113] += dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[3];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[178] -= dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 33: H2O2 + OH => H2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 1000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[112] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    J[113] += dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[248] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[251] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 34: H2O2 + H => H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 24100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3969.8899999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (3969.8899999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[OH]/d[H] */
    J[76] += dqdci;               /* dwdot[H2O]/d[H] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[247] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[248] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[251] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 35: H2O2 + OH => H2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 580000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (9559.9899999999998) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (9559.9899999999998)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[112] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    J[113] += dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[248] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[251] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 36: H2O + HO2 => H2O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[8];
    k_f = 1.0000000000000002e-06 * 97710000000000
                * exp(0.33000000000000002 * tc[0] - 0.50321666580471969 * (41020.080000000002) * invT);
    dlnkfdT = 0.33000000000000002 * invT + 0.50321666580471969 *  (41020.080000000002)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* H2O */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[8];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O2]/d[H2O] */
    J[218] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[6];
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[286] -= dqdci;              /* dwdot[H2O]/d[HO2] */
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 37: H2O2 + O => OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = 1.0000000000000002e-06 * 9550000
                * exp(2 * tc[0] - 0.50321666580471969 * (3969.8899999999999) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (3969.8899999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[42] -= dqdci;               /* dwdot[H2O2]/d[O] */
    J[43] += dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[246] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[248] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 38: H2O2 + H => H2 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 48200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7950.0500000000002) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (7950.0500000000002)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H] */
    J[78] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[247] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[249] += dqdci;              /* dwdot[H2]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 39: H2 + HO2 => H2O2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[8];
    k_f = 1.0000000000000002e-06 * 1875000000000
                * exp(0.33000000000000002 * tc[0] - 0.50321666580471969 * (24260.040000000001) * invT);
    dlnkfdT = 0.33000000000000002 * invT + 0.50321666580471969 *  (24260.040000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[8];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[147] += dqdci;              /* dwdot[H2O2]/d[H2] */
    J[148] -= dqdci;              /* dwdot[HO2]/d[H2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[282] += dqdci;              /* dwdot[H]/d[HO2] */
    J[284] -= dqdci;              /* dwdot[H2]/d[HO2] */
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 40: CH2GSG + OH => CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[9] -= q; /* CH2GSG */
    wdot[10] += q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[107] += dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[114] -= dqdci;              /* dwdot[CH2GSG]/d[OH] */
    J[115] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[3];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[318] -= dqdci;              /* dwdot[OH]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[325] += dqdci;              /* dwdot[CH2O]/d[CH2GSG] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 41: CH2GSG + H2 => CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = 1.0000000000000002e-06 * 70000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[9];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[149] -= dqdci;              /* dwdot[CH2GSG]/d[H2] */
    J[151] += dqdci;              /* dwdot[CH3]/d[H2] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[4];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[319] -= dqdci;              /* dwdot[H2]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] += dqdci;              /* dwdot[CH3]/d[CH2GSG] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */

    /*reaction 42: CH3 + H => CH2GSG + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = 1.0000000000000002e-06 * 2.482e+17
                * exp(-0.89000000000000001 * tc[0] - 0.50321666580471969 * (16130.02) * invT);
    dlnkfdT = -0.89000000000000001 * invT + 0.50321666580471969 *  (16130.02)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[9] += q; /* CH2GSG */
    wdot[11] -= q; /* CH3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[11];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[79] += dqdci;               /* dwdot[CH2GSG]/d[H] */
    J[81] -= dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[2];
    J[387] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[389] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[394] += dqdci;              /* dwdot[CH2GSG]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1199] += dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 43: CH2GSG + O2 => CO + OH + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[9];
    k_f = 1.0000000000000002e-06 * 70000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[9] -= q; /* CH2GSG */
    wdot[12] += q; /* CO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[177] += dqdci;              /* dwdot[H]/d[O2] */
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[184] -= dqdci;              /* dwdot[CH2GSG]/d[O2] */
    J[187] += dqdci;              /* dwdot[CO]/d[O2] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[5];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[318] += dqdci;              /* dwdot[OH]/d[CH2GSG] */
    J[320] -= dqdci;              /* dwdot[O2]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[327] += dqdci;              /* dwdot[CO]/d[CH2GSG] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 44: CH2GSG + O => CO + 2.000000 H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += 2 * q; /* H */
    wdot[9] -= q; /* CH2GSG */
    wdot[12] += q; /* CO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[9];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += 2 * dqdci;           /* dwdot[H]/d[O] */
    J[44] -= dqdci;               /* dwdot[CH2GSG]/d[O] */
    J[47] += dqdci;               /* dwdot[CO]/d[O] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[1];
    J[316] -= dqdci;              /* dwdot[O]/d[CH2GSG] */
    J[317] += 2 * dqdci;          /* dwdot[H]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[327] += dqdci;              /* dwdot[CO]/d[CH2GSG] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += 2 * dqdT;          /* dwdot[H]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 45: CH3 + HO2 => CH3O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = 1.0000000000000002e-06 * 11000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[11] -= q; /* CH3 */
    wdot[13] += q; /* CH3O */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[291] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[293] += dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[388] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[393] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[398] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1203] += dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 46: CH3 + O2 => CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = 1.0000000000000002e-06 * 747000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (14250) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (14250)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[11] -= q; /* CH3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[186] -= dqdci;              /* dwdot[CH3]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[5];
    J[388] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[390] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[395] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 47: CH3 + O2 => CH3O + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = 1.0000000000000002e-06 * 1.995e+18
                * exp(-1.5700000000000001 * tc[0] - 0.50321666580471969 * (29210.09) * invT);
    dlnkfdT = -1.5700000000000001 * invT + 0.50321666580471969 *  (29210.09)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[11] -= q; /* CH3 */
    wdot[13] += q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[186] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[188] += dqdci;              /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[5];
    J[386] += dqdci;              /* dwdot[O]/d[CH3] */
    J[390] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[398] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1203] += dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 48: CH3O + O => CH3 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = 1.0000000000000002e-06 * 3.585e+18
                * exp(-1.5900000000000001 * tc[0] - 0.50321666580471969 * (-1630.98) * invT);
    dlnkfdT = -1.5900000000000001 * invT + 0.50321666580471969 *  (-1630.98)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[5] += q; /* O2 */
    wdot[11] += q; /* CH3 */
    wdot[13] -= q; /* CH3O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[40] += dqdci;               /* dwdot[O2]/d[O] */
    J[46] += dqdci;               /* dwdot[CH3]/d[O] */
    J[48] -= dqdci;               /* dwdot[CH3O]/d[O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[1];
    J[456] -= dqdci;              /* dwdot[O]/d[CH3O] */
    J[460] += dqdci;              /* dwdot[O2]/d[CH3O] */
    J[466] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[468] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1203] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 49: 2.000000 CH3 <=> H + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[11], 2.000000);
    k_f = 1.0000000000000002e-06 * 6840000000000
                * exp(0.10000000000000001 * tc[0] - 0.50321666580471969 * (10599.9) * invT);
    dlnkfdT = 0.10000000000000001 * invT + 0.50321666580471969 *  (10599.9)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[14];
    Kc = exp(-g_RT[2] + 2.000000*g_RT[11] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[11]) + (h_RT[2] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[11] -= 2 * q; /* CH3 */
    wdot[14] += q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[72] += dqdci;               /* dwdot[H]/d[H] */
    J[81] += -2 * dqdci;          /* dwdot[CH3]/d[H] */
    J[84] += dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*2.000000*sc[11];
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[396] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[399] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[2];
    J[492] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[501] += -2 * dqdci;         /* dwdot[CH3]/d[C2H5] */
    J[504] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1201] += -2 * dqdT;         /* dwdot[CH3]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */

    /*reaction 50: CH3 + HO2 => CH4 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = 1.0000000000000002e-06 * 3600000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[291] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[295] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[390] += dqdci;              /* dwdot[O2]/d[CH3] */
    J[393] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 51: CH3 + O => CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = 1.0000000000000002e-06 * 80000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[10] += q; /* CH2O */
    wdot[11] -= q; /* CH3 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[11];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[45] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[46] -= dqdci;               /* dwdot[CH3]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[1];
    J[386] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[395] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 52: CH3 + OH => CH2O + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = 1.0000000000000002e-06 * 22500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (4299.9499999999998) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (4299.9499999999998)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2 */
    wdot[10] += q; /* CH2O */
    wdot[11] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[109] += dqdci;              /* dwdot[H2]/d[OH] */
    J[115] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[116] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[388] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[389] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[395] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 53: CH3 + OH => CH2GSG + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = 1.0000000000000002e-06 * 26500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (2185.9499999999998) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (2185.9499999999998)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[9] += q; /* CH2GSG */
    wdot[11] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[114] += dqdci;              /* dwdot[CH2GSG]/d[OH] */
    J[116] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[388] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[391] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[394] += dqdci;              /* dwdot[CH2GSG]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1199] += dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 54: CH2GSG + H2O => CH3 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[9];
    k_f = 1.0000000000000002e-06 * 32360000000
                * exp(0.89000000000000001 * tc[0] - 0.50321666580471969 * (1211.04) * invT);
    dlnkfdT = 0.89000000000000001 * invT + 0.50321666580471969 *  (1211.04)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* H2O */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] += q; /* CH3 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[9];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[219] -= dqdci;              /* dwdot[CH2GSG]/d[H2O] */
    J[221] += dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[6];
    J[318] += dqdci;              /* dwdot[OH]/d[CH2GSG] */
    J[321] -= dqdci;              /* dwdot[H2O]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] += dqdci;              /* dwdot[CH3]/d[CH2GSG] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */

    /*reaction 55: CH3 + H2O2 => CH4 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[11];
    k_f = 1.0000000000000002e-06 * 336500000000
                * exp(-0.33000000000000002 * tc[0] - 0.50321666580471969 * (2501.9099999999999) * invT);
    dlnkfdT = -0.33000000000000002 * invT + 0.50321666580471969 *  (2501.9099999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[11];
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[256] -= dqdci;              /* dwdot[CH3]/d[H2O2] */
    J[260] += dqdci;              /* dwdot[CH4]/d[H2O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[392] -= dqdci;              /* dwdot[H2O2]/d[CH3] */
    J[393] += dqdci;              /* dwdot[HO2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 56: CH2GSG + CH3 => C2H4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[11];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] -= q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[11];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] -= dqdci;              /* dwdot[CH3]/d[CH2GSG] */
    J[331] += dqdci;              /* dwdot[C2H4]/d[CH2GSG] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[9];
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[394] -= dqdci;              /* dwdot[CH2GSG]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[401] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 57: CH2GSG + CH4 => 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[15];
    k_f = 1.0000000000000002e-06 * 40000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] += 2 * q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[15];
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] += 2 * dqdci;          /* dwdot[CH3]/d[CH2GSG] */
    J[330] -= dqdci;              /* dwdot[CH4]/d[CH2GSG] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[9];
    J[534] -= dqdci;              /* dwdot[CH2GSG]/d[CH4] */
    J[536] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += 2 * dqdT;          /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 58: 2.000000 CH3 => CH2GSG + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[11], 2.000000);
    k_f = 1.0000000000000002e-06 * 5429000000000000
                * exp(-0.89000000000000001 * tc[0] - 0.50321666580471969 * (15650.1) * invT);
    dlnkfdT = -0.89000000000000001 * invT + 0.50321666580471969 *  (15650.1)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] += q; /* CH2GSG */
    wdot[11] -= 2 * q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[CH3] */
    dqdci =  + k_f*2.000000*sc[11];
    J[394] += dqdci;              /* dwdot[CH2GSG]/d[CH3] */
    J[396] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1199] += dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += -2 * dqdT;         /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 59: CH4 + O => CH3 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[15];
    k_f = 1.0000000000000002e-06 * 3150000000000
                * exp(0.5 * tc[0] - 0.50321666580471969 * (10289.91) * invT);
    dlnkfdT = 0.5 * invT + 0.50321666580471969 *  (10289.91)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[11] += q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[15];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[46] += dqdci;               /* dwdot[CH3]/d[O] */
    J[50] -= dqdci;               /* dwdot[CH4]/d[O] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[1];
    J[526] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[528] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[536] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 60: CH3 + OH => CH4 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = 1.0000000000000002e-06 * 52960000000
                * exp(0.5 * tc[0] - 0.50321666580471969 * (7715.1099999999997) * invT);
    dlnkfdT = 0.5 * invT + 0.50321666580471969 *  (7715.1099999999997)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[3] -= q; /* OH */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[106] += dqdci;              /* dwdot[O]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[120] += dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[386] += dqdci;              /* dwdot[O]/d[CH3] */
    J[388] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 61: CH4 + H => CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[15];
    k_f = 1.0000000000000002e-06 * 17270
                * exp(3 * tc[0] - 0.50321666580471969 * (8223.9500000000007) * invT);
    dlnkfdT = 3 * invT + 0.50321666580471969 *  (8223.9500000000007)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[11] += q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[15];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[81] += dqdci;               /* dwdot[CH3]/d[H] */
    J[85] -= dqdci;               /* dwdot[CH4]/d[H] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[527] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[529] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[536] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 62: CH3 + H2 => CH4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = 1.0000000000000002e-06 * 661
                * exp(3 * tc[0] - 0.50321666580471969 * (7744.0200000000004) * invT);
    dlnkfdT = 3 * invT + 0.50321666580471969 *  (7744.0200000000004)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[11];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[151] -= dqdci;              /* dwdot[CH3]/d[H2] */
    J[155] += dqdci;              /* dwdot[CH4]/d[H2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[389] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 63: CH4 + OH => CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = 1.0000000000000002e-06 * 193000
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (2106.1199999999999) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  (2106.1199999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[11] += q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[15];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[116] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[120] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[3];
    J[528] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[531] += dqdci;              /* dwdot[H2O]/d[CH4] */
    J[536] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 64: CH3 + H2O => CH4 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[11];
    k_f = 1.0000000000000002e-06 * 31990
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (16780.110000000001) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  (16780.110000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* H2O */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[11];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[221] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    J[225] += dqdci;              /* dwdot[CH4]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[388] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[391] -= dqdci;              /* dwdot[H2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 65: CO + HO2 => CO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = 1.0000000000000002e-06 * 30100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (23000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (23000)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[12] -= q; /* CO */
    wdot[18] += q; /* CO2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[12];
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[292] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[298] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[8];
    J[423] += dqdci;              /* dwdot[OH]/d[CO] */
    J[428] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[438] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1202] -= dqdT;              /* dwdot[CO]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 66: CO + O2 => CO2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[12];
    k_f = 1.0000000000000002e-06 * 16200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (47700.050000000003) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (47700.050000000003)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[12] -= q; /* CO */
    wdot[18] += q; /* CO2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[12];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[187] -= dqdci;              /* dwdot[CO]/d[O2] */
    J[193] += dqdci;              /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[5];
    J[421] += dqdci;              /* dwdot[O]/d[CO] */
    J[425] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[438] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1202] -= dqdT;              /* dwdot[CO]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 67: CO2 + O => CO + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[18];
    k_f = 1.0000000000000002e-06 * 143300000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (53919.93) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (53919.93)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[5] += q; /* O2 */
    wdot[12] += q; /* CO */
    wdot[18] -= q; /* CO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[18];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[40] += dqdci;               /* dwdot[O2]/d[O] */
    J[47] += dqdci;               /* dwdot[CO]/d[O] */
    J[53] -= dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[1];
    J[631] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[635] += dqdci;              /* dwdot[O2]/d[CO2] */
    J[642] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[648] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1208] -= dqdT;              /* dwdot[CO2]/dT */

    /*reaction 68: CO + OH => CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[12];
    k_f = 1.0000000000000002e-06 * 140000
                * exp(1.95 * tc[0] - 0.50321666580471969 * (-1347.04) * invT);
    dlnkfdT = 1.95 * invT + 0.50321666580471969 *  (-1347.04)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[12] -= q; /* CO */
    wdot[18] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[107] += dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[117] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[123] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3];
    J[422] += dqdci;              /* dwdot[H]/d[CO] */
    J[423] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[438] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1202] -= dqdT;              /* dwdot[CO]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 69: CO2 + H => CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = 1.0000000000000002e-06 * 15680000
                * exp(1.95 * tc[0] - 0.50321666580471969 * (20989.959999999999) * invT);
    dlnkfdT = 1.95 * invT + 0.50321666580471969 *  (20989.959999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[12] += q; /* CO */
    wdot[18] -= q; /* CO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[OH]/d[H] */
    J[82] += dqdci;               /* dwdot[CO]/d[H] */
    J[88] -= dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[2];
    J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
    J[633] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[642] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[648] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1208] -= dqdT;              /* dwdot[CO2]/dT */

    /*reaction 70: HCO + CH3 => CH4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[19];
    k_f = 1.0000000000000002e-06 * 121000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[12] += q; /* CO */
    wdot[15] += q; /* CH4 */
    wdot[19] -= q; /* HCO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[19];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[397] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[404] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[11];
    J[676] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[680] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 71: HCO + H => CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[19];
    k_f = 1.0000000000000002e-06 * 73400000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[19];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[82] += dqdci;               /* dwdot[CO]/d[H] */
    J[89] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[667] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[669] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 72: HCO + O2 => CO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[19];
    k_f = 1.0000000000000002e-06 * 7580000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (409.88999999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (409.88999999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[19];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[187] += dqdci;              /* dwdot[CO]/d[O2] */
    J[194] -= dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[670] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[673] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 73: HCO + O => CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = 1.0000000000000002e-06 * 30200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[19];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[47] += dqdci;               /* dwdot[CO]/d[O] */
    J[54] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[666] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[668] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 74: HCO + O => CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[18] += q; /* CO2 */
    wdot[19] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[19];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[53] += dqdci;               /* dwdot[CO2]/d[O] */
    J[54] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[666] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[667] += dqdci;              /* dwdot[H]/d[HCO] */
    J[683] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 75: HCO + OH => CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[19];
    k_f = 1.0000000000000002e-06 * 102000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[19];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[117] += dqdci;              /* dwdot[CO]/d[OH] */
    J[124] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[668] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[671] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 76: CH2O + H => HCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[10];
    k_f = 1.0000000000000002e-06 * 933400000
                * exp(1.5 * tc[0] - 0.50321666580471969 * (2976.0999999999999) * invT);
    dlnkfdT = 1.5 * invT + 0.50321666580471969 *  (2976.0999999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[80] -= dqdci;               /* dwdot[CH2O]/d[H] */
    J[89] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[2];
    J[352] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[354] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 77: CH2O + O2 => HCO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[10];
    k_f = 1.0000000000000002e-06 * 20500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (38950.050000000003) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (38950.050000000003)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[185] -= dqdci;              /* dwdot[CH2O]/d[O2] */
    J[194] += dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[5];
    J[355] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[358] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 78: CH2O + OH => HCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[10];
    k_f = 1.0000000000000002e-06 * 3430000000
                * exp(1.1799999999999999 * tc[0] - 0.50321666580471969 * (-446.94) * invT);
    dlnkfdT = 1.1799999999999999 * invT + 0.50321666580471969 *  (-446.94)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[115] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    J[124] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[3];
    J[353] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[356] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 79: CH2O + HO2 => HCO + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = 1.0000000000000002e-06 * 0.0058199999999999997
                * exp(4.5300000000000002 * tc[0] - 0.50321666580471969 * (6556.8800000000001) * invT);
    dlnkfdT = 4.5300000000000002 * invT + 0.50321666580471969 *  (6556.8800000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[10];
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[290] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[8];
    J[357] += dqdci;              /* dwdot[H2O2]/d[CH2O] */
    J[358] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 80: CH2O + O => HCO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = 1.0000000000000002e-06 * 416000000000
                * exp(0.56999999999999995 * tc[0] - 0.50321666580471969 * (2761.9499999999998) * invT);
    dlnkfdT = 0.56999999999999995 * invT + 0.50321666580471969 *  (2761.9499999999998)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[45] -= dqdci;               /* dwdot[CH2O]/d[O] */
    J[54] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[1];
    J[351] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[353] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 81: CH2O + CH3 => HCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[11];
    k_f = 1.0000000000000002e-06 * 3.636e-06
                * exp(5.4199999999999999 * tc[0] - 0.50321666580471969 * (998.09000000000003) * invT);
    dlnkfdT = 5.4199999999999999 * invT + 0.50321666580471969 *  (998.09000000000003)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] -= q; /* CH2O */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[19] += q; /* HCO */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[11];
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[361] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[365] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[10];
    J[395] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[404] += dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 82: CH3O + O2 => CH2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = 1.0000000000000002e-06 * 55000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (2424) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (2424)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[10] += q; /* CH2O */
    wdot[13] -= q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[188] -= dqdci;              /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[5];
    J[460] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[463] += dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[465] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[468] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1203] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 83: CH2GSG + CO2 => CH2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[18];
    k_f = 1.0000000000000002e-06 * 3000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] -= q; /* CH2GSG */
    wdot[10] += q; /* CH2O */
    wdot[12] += q; /* CO */
    wdot[18] -= q; /* CO2 */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[18];
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[325] += dqdci;              /* dwdot[CH2O]/d[CH2GSG] */
    J[327] += dqdci;              /* dwdot[CO]/d[CH2GSG] */
    J[333] -= dqdci;              /* dwdot[CO2]/d[CH2GSG] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[9];
    J[639] -= dqdci;              /* dwdot[CH2GSG]/d[CO2] */
    J[640] += dqdci;              /* dwdot[CH2O]/d[CO2] */
    J[642] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[648] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1208] -= dqdT;              /* dwdot[CO2]/dT */

    /*reaction 84: CH3O2 + CH3 => 2.000000 CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[20];
    k_f = 1.0000000000000002e-06 * 7000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1000)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[13] += 2 * q; /* CH3O */
    wdot[20] -= q; /* CH3O2 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[20];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[398] += 2 * dqdci;          /* dwdot[CH3O]/d[CH3] */
    J[405] -= dqdci;              /* dwdot[CH3O2]/d[CH3] */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[11];
    J[711] -= dqdci;              /* dwdot[CH3]/d[CH3O2] */
    J[713] += 2 * dqdci;          /* dwdot[CH3O]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1203] += 2 * dqdT;          /* dwdot[CH3O]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */

    /*reaction 85: CH3O2 + HO2 => CH3O2H + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[20];
    k_f = 1.0000000000000002e-06 * 17500000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-3275.0999999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-3275.0999999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[20] -= q; /* CH3O2 */
    wdot[21] += q; /* CH3O2H */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[20];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[300] -= dqdci;              /* dwdot[CH3O2]/d[HO2] */
    J[301] += dqdci;              /* dwdot[CH3O2H]/d[HO2] */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[8];
    J[705] += dqdci;              /* dwdot[O2]/d[CH3O2] */
    J[708] -= dqdci;              /* dwdot[HO2]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    J[721] += dqdci;              /* dwdot[CH3O2H]/d[CH3O2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */
    J[1211] += dqdT;              /* dwdot[CH3O2H]/dT */

    /*reaction 86: 2.000000 CH3O2 => O2 + 2.000000 CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[20], 2.000000);
    k_f = 1.0000000000000002e-06 * 14000000000000000
                * exp(-1.6100000000000001 * tc[0] - 0.50321666580471969 * (1859.9400000000001) * invT);
    dlnkfdT = -1.6100000000000001 * invT + 0.50321666580471969 *  (1859.9400000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[13] += 2 * q; /* CH3O */
    wdot[20] -= 2 * q; /* CH3O2 */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*2.000000*sc[20];
    J[705] += dqdci;              /* dwdot[O2]/d[CH3O2] */
    J[713] += 2 * dqdci;          /* dwdot[CH3O]/d[CH3O2] */
    J[720] += -2 * dqdci;         /* dwdot[CH3O2]/d[CH3O2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1203] += 2 * dqdT;          /* dwdot[CH3O]/dT */
    J[1210] += -2 * dqdT;         /* dwdot[CH3O2]/dT */

    /*reaction 87: CH3O2 + CH2O => CH3O2H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[20];
    k_f = 1.0000000000000002e-06 * 1990000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (11659.889999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (11659.889999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    wdot[20] -= q; /* CH3O2 */
    wdot[21] += q; /* CH3O2H */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[20];
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[370] -= dqdci;              /* dwdot[CH3O2]/d[CH2O] */
    J[371] += dqdci;              /* dwdot[CH3O2H]/d[CH2O] */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[10];
    J[710] -= dqdci;              /* dwdot[CH2O]/d[CH3O2] */
    J[719] += dqdci;              /* dwdot[HCO]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    J[721] += dqdci;              /* dwdot[CH3O2H]/d[CH3O2] */
    /* d()/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */
    J[1211] += dqdT;              /* dwdot[CH3O2H]/dT */

    /*reaction 88: CH3O2H => CH3O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[21];
    k_f = 1 * 631000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (42299.949999999997) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (42299.949999999997)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[13] += q; /* CH3O */
    wdot[21] -= q; /* CH3O2H */
    /* d()/d[CH3O2H] */
    dqdci =  + k_f;
    J[738] += dqdci;              /* dwdot[OH]/d[CH3O2H] */
    J[748] += dqdci;              /* dwdot[CH3O]/d[CH3O2H] */
    J[756] -= dqdci;              /* dwdot[CH3O2H]/d[CH3O2H] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1203] += dqdT;              /* dwdot[CH3O]/dT */
    J[1211] -= dqdT;              /* dwdot[CH3O2H]/dT */

    /*reaction 89: CH3O + OH => CH3O2H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[13];
    k_f = 1.0000000000000002e-06 * 116600000000
                * exp(0.59999999999999998 * tc[0] - 0.50321666580471969 * (-1771.03) * invT);
    dlnkfdT = 0.59999999999999998 * invT + 0.50321666580471969 *  (-1771.03)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[13] -= q; /* CH3O */
    wdot[21] += q; /* CH3O2H */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[118] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    J[126] += dqdci;              /* dwdot[CH3O2H]/d[OH] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[3];
    J[458] -= dqdci;              /* dwdot[OH]/d[CH3O] */
    J[468] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[476] += dqdci;              /* dwdot[CH3O2H]/d[CH3O] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1203] -= dqdT;              /* dwdot[CH3O]/dT */
    J[1211] += dqdT;              /* dwdot[CH3O2H]/dT */

    /*reaction 90: C2H2 + O2 => HCCO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[22];
    k_f = 1.0000000000000002e-06 * 200000000
                * exp(1.5 * tc[0] - 0.50321666580471969 * (30099.900000000001) * invT);
    dlnkfdT = 1.5 * invT + 0.50321666580471969 *  (30099.900000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[22] -= q; /* C2H2 */
    wdot[23] += q; /* HCCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[22];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[197] -= dqdci;              /* dwdot[C2H2]/d[O2] */
    J[198] += dqdci;              /* dwdot[HCCO]/d[O2] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[5];
    J[773] += dqdci;              /* dwdot[OH]/d[C2H2] */
    J[775] -= dqdci;              /* dwdot[O2]/d[C2H2] */
    J[792] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[793] += dqdci;              /* dwdot[HCCO]/d[C2H2] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1212] -= dqdT;              /* dwdot[C2H2]/dT */
    J[1213] += dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 91: C2H2 + O => HCCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[22];
    k_f = 1.0000000000000002e-06 * 14300000
                * exp(2 * tc[0] - 0.50321666580471969 * (1900.0999999999999) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (1900.0999999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[22] -= q; /* C2H2 */
    wdot[23] += q; /* HCCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[22];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[57] -= dqdci;               /* dwdot[C2H2]/d[O] */
    J[58] += dqdci;               /* dwdot[HCCO]/d[O] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[1];
    J[771] -= dqdci;              /* dwdot[O]/d[C2H2] */
    J[772] += dqdci;              /* dwdot[H]/d[C2H2] */
    J[792] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[793] += dqdci;              /* dwdot[HCCO]/d[C2H2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1212] -= dqdT;              /* dwdot[C2H2]/dT */
    J[1213] += dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 92: C2H3 + O2 => CH2CHO + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[24];
    k_f = 1.0000000000000002e-06 * 350000000000000
                * exp(-0.60999999999999999 * tc[0] - 0.50321666580471969 * (5260.04) * invT);
    dlnkfdT = -0.60999999999999999 * invT + 0.50321666580471969 *  (5260.04)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[24] -= q; /* C2H3 */
    wdot[25] += q; /* CH2CHO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[24];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[199] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    J[200] += dqdci;              /* dwdot[CH2CHO]/d[O2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[5];
    J[841] += dqdci;              /* dwdot[O]/d[C2H3] */
    J[845] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[865] += dqdci;              /* dwdot[CH2CHO]/d[C2H3] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */
    J[1215] += dqdT;              /* dwdot[CH2CHO]/dT */

    /*reaction 93: C2H3 + H => C2H2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[24];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (2500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (2500)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[22] += q; /* C2H2 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[24];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[92] += dqdci;               /* dwdot[C2H2]/d[H] */
    J[94] -= dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[2];
    J[842] -= dqdci;              /* dwdot[H]/d[C2H3] */
    J[844] += dqdci;              /* dwdot[H2]/d[C2H3] */
    J[862] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1212] += dqdT;              /* dwdot[C2H2]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 94: C2H3 + O2 => CH2O + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[24];
    k_f = 1.0000000000000002e-06 * 1.6999999999999999e+29
                * exp(-5.3099999999999996 * tc[0] - 0.50321666580471969 * (6500) * invT);
    dlnkfdT = -5.3099999999999996 * invT + 0.50321666580471969 *  (6500)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[19] += q; /* HCO */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[24];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[194] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[199] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[5];
    J[845] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[850] += dqdci;              /* dwdot[CH2O]/d[C2H3] */
    J[859] += dqdci;              /* dwdot[HCO]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 95: C2H3 + O2 => C2H2 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[24];
    k_f = 1.0000000000000002e-06 * 2.12e-06
                * exp(6 * tc[0] - 0.50321666580471969 * (9483.9899999999998) * invT);
    dlnkfdT = 6 * invT + 0.50321666580471969 *  (9483.9899999999998)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[22] += q; /* C2H2 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[24];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[197] += dqdci;              /* dwdot[C2H2]/d[O2] */
    J[199] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[5];
    J[845] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[848] += dqdci;              /* dwdot[HO2]/d[C2H3] */
    J[862] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1212] += dqdT;              /* dwdot[C2H2]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 96: C2H4 + O => CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = 1.0000000000000002e-06 * 10200000
                * exp(1.8799999999999999 * tc[0] - 0.50321666580471969 * (179.02000000000001) * invT);
    dlnkfdT = 1.8799999999999999 * invT + 0.50321666580471969 *  (179.02000000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[11] += q; /* CH3 */
    wdot[16] -= q; /* C2H4 */
    wdot[19] += q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[46] += dqdci;               /* dwdot[CH3]/d[O] */
    J[51] -= dqdci;               /* dwdot[C2H4]/d[O] */
    J[54] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[1];
    J[561] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[571] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[579] += dqdci;              /* dwdot[HCO]/d[C2H4] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 97: C2H4 + OH => C2H3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[16];
    k_f = 1.0000000000000002e-06 * 20500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (5950.0500000000002) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (5950.0500000000002)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[16] -= q; /* C2H4 */
    wdot[24] += q; /* C2H3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[16];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[121] -= dqdci;              /* dwdot[C2H4]/d[OH] */
    J[129] += dqdci;              /* dwdot[C2H3]/d[OH] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[3];
    J[563] -= dqdci;              /* dwdot[OH]/d[C2H4] */
    J[566] += dqdci;              /* dwdot[H2O]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[584] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 98: C2H3 + H2O => C2H4 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[24];
    k_f = 1.0000000000000002e-06 * 6033000000000000
                * exp(-0.82999999999999996 * tc[0] - 0.50321666580471969 * (21760.040000000001) * invT);
    dlnkfdT = -0.82999999999999996 * invT + 0.50321666580471969 *  (21760.040000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* H2O */
    wdot[16] += q; /* C2H4 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[24];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[226] += dqdci;              /* dwdot[C2H4]/d[H2O] */
    J[234] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[6];
    J[843] += dqdci;              /* dwdot[OH]/d[C2H3] */
    J[846] -= dqdci;              /* dwdot[H2O]/d[C2H3] */
    J[856] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 99: C2H4 + H => C2H3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[16];
    k_f = 1.0000000000000002e-06 * 0.0084200000000000004
                * exp(4.6200000000000001 * tc[0] - 0.50321666580471969 * (2582.9299999999998) * invT);
    dlnkfdT = 4.6200000000000001 * invT + 0.50321666580471969 *  (2582.9299999999998)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[16] -= q; /* C2H4 */
    wdot[24] += q; /* C2H3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[16];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[86] -= dqdci;               /* dwdot[C2H4]/d[H] */
    J[94] += dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[2];
    J[562] -= dqdci;              /* dwdot[H]/d[C2H4] */
    J[564] += dqdci;              /* dwdot[H2]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[584] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 100: C2H3 + H2 => C2H4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[24];
    k_f = 1.0000000000000002e-06 * 0.57230000000000003
                * exp(3.79 * tc[0] - 0.50321666580471969 * (3233.0300000000002) * invT);
    dlnkfdT = 3.79 * invT + 0.50321666580471969 *  (3233.0300000000002)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[16] += q; /* C2H4 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[24];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[156] += dqdci;              /* dwdot[C2H4]/d[H2] */
    J[164] -= dqdci;              /* dwdot[C2H3]/d[H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[4];
    J[842] += dqdci;              /* dwdot[H]/d[C2H3] */
    J[844] -= dqdci;              /* dwdot[H2]/d[C2H3] */
    J[856] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 101: C2H4 + CH3O2 => C2H3 + CH3O2H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[16]*sc[20];
    k_f = 1.0000000000000002e-06 * 2230000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (17190.009999999998) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (17190.009999999998)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[16] -= q; /* C2H4 */
    wdot[20] -= q; /* CH3O2 */
    wdot[21] += q; /* CH3O2H */
    wdot[24] += q; /* C2H3 */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[20];
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[580] -= dqdci;              /* dwdot[CH3O2]/d[C2H4] */
    J[581] += dqdci;              /* dwdot[CH3O2H]/d[C2H4] */
    J[584] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[16];
    J[716] -= dqdci;              /* dwdot[C2H4]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    J[721] += dqdci;              /* dwdot[CH3O2H]/d[CH3O2] */
    J[724] += dqdci;              /* dwdot[C2H3]/d[CH3O2] */
    /* d()/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */
    J[1211] += dqdT;              /* dwdot[CH3O2H]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 102: C2H4 + CH3 => C2H3 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[16];
    k_f = 1.0000000000000002e-06 * 6.6200000000000001
                * exp(3.7000000000000002 * tc[0] - 0.50321666580471969 * (9500) * invT);
    dlnkfdT = 3.7000000000000002 * invT + 0.50321666580471969 *  (9500)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[16] -= q; /* C2H4 */
    wdot[24] += q; /* C2H3 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[401] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[409] += dqdci;              /* dwdot[C2H3]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[11];
    J[571] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[575] += dqdci;              /* dwdot[CH4]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[584] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 103: C2H3 + CH4 => C2H4 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[15]*sc[24];
    k_f = 1.0000000000000002e-06 * 1.4399999999999999
                * exp(4.0199999999999996 * tc[0] - 0.50321666580471969 * (5472.04) * invT);
    dlnkfdT = 4.0199999999999996 * invT + 0.50321666580471969 *  (5472.04)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] += q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    wdot[16] += q; /* C2H4 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[24];
    J[536] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    J[541] += dqdci;              /* dwdot[C2H4]/d[CH4] */
    J[549] -= dqdci;              /* dwdot[C2H3]/d[CH4] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[15];
    J[851] += dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[855] -= dqdci;              /* dwdot[CH4]/d[C2H3] */
    J[856] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 104: C2H4 + O => CH2CHO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = 1.0000000000000002e-06 * 3390000
                * exp(1.8799999999999999 * tc[0] - 0.50321666580471969 * (179.02000000000001) * invT);
    dlnkfdT = 1.8799999999999999 * invT + 0.50321666580471969 *  (179.02000000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[16] -= q; /* C2H4 */
    wdot[25] += q; /* CH2CHO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[51] -= dqdci;               /* dwdot[C2H4]/d[O] */
    J[60] += dqdci;               /* dwdot[CH2CHO]/d[O] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[1];
    J[561] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[562] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[585] += dqdci;              /* dwdot[CH2CHO]/d[C2H4] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1215] += dqdT;              /* dwdot[CH2CHO]/dT */

    /*reaction 105: CH3 + C2H5 => CH4 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[14];
    k_f = 1.0000000000000002e-06 * 19500000000000
                * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.5 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[14] -= q; /* C2H5 */
    wdot[15] += q; /* CH4 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[14];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[399] -= dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[401] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[11];
    J[501] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[505] += dqdci;              /* dwdot[CH4]/d[C2H5] */
    J[506] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 106: C2H5 + O2 => C2H4 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[14];
    k_f = 1.0000000000000002e-06 * 1.22e+30
                * exp(-5.7599999999999998 * tc[0] - 0.50321666580471969 * (10099.9) * invT);
    dlnkfdT = -5.7599999999999998 * invT + 0.50321666580471969 *  (10099.9)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[14] -= q; /* C2H5 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[189] -= dqdci;              /* dwdot[C2H5]/d[O2] */
    J[191] += dqdci;              /* dwdot[C2H4]/d[O2] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[5];
    J[495] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[498] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[506] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 107: C2H4 + HO2 => C2H5 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[16];
    k_f = 1.0000000000000002e-06 * 1.259e+30
                * exp(-5.6299999999999999 * tc[0] - 0.50321666580471969 * (22299.950000000001) * invT);
    dlnkfdT = -5.6299999999999999 * invT + 0.50321666580471969 *  (22299.950000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[14] += q; /* C2H5 */
    wdot[16] -= q; /* C2H4 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[16];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[294] += dqdci;              /* dwdot[C2H5]/d[HO2] */
    J[296] -= dqdci;              /* dwdot[C2H4]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[8];
    J[565] += dqdci;              /* dwdot[O2]/d[C2H4] */
    J[568] -= dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[574] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 108: C2H5 + HO2 => C2H6 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[14];
    k_f = 1.0000000000000002e-06 * 267900000
                * exp(0.89000000000000001 * tc[0] - 0.50321666580471969 * (-1922.0799999999999) * invT);
    dlnkfdT = 0.89000000000000001 * invT + 0.50321666580471969 *  (-1922.0799999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[14] -= q; /* C2H5 */
    wdot[17] += q; /* C2H6 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[14];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[294] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    J[297] += dqdci;              /* dwdot[C2H6]/d[HO2] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[8];
    J[495] += dqdci;              /* dwdot[O2]/d[C2H5] */
    J[498] -= dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[507] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1207] += dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 109: H + C2H5 => C2H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[14];
    k_f = 1.0000000000000002e-06 * 583100000000
                * exp(0.59899999999999998 * tc[0] - 0.50321666580471969 * (-2913) * invT);
    dlnkfdT = 0.59899999999999998 * invT + 0.50321666580471969 *  (-2913)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[14] -= q; /* C2H5 */
    wdot[17] += q; /* C2H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[14];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[84] -= dqdci;               /* dwdot[C2H5]/d[H] */
    J[87] += dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[2];
    J[492] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[507] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1207] += dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 110: C2H6 + H => C2H5 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = 1.0000000000000002e-06 * 554
                * exp(3.5 * tc[0] - 0.50321666580471969 * (5167.0699999999997) * invT);
    dlnkfdT = 3.5 * invT + 0.50321666580471969 *  (5167.0699999999997)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[14] += q; /* C2H5 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[17];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[84] += dqdci;               /* dwdot[C2H5]/d[H] */
    J[87] -= dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[2];
    J[597] -= dqdci;              /* dwdot[H]/d[C2H6] */
    J[599] += dqdci;              /* dwdot[H2]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 111: C2H5 + H2 => C2H6 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[14];
    k_f = 1.0000000000000002e-06 * 0.13550000000000001
                * exp(4.0599999999999996 * tc[0] - 0.50321666580471969 * (8857.0699999999997) * invT);
    dlnkfdT = 4.0599999999999996 * invT + 0.50321666580471969 *  (8857.0699999999997)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[14] -= q; /* C2H5 */
    wdot[17] += q; /* C2H6 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[14];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[154] -= dqdci;              /* dwdot[C2H5]/d[H2] */
    J[157] += dqdci;              /* dwdot[C2H6]/d[H2] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[4];
    J[492] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[494] -= dqdci;              /* dwdot[H2]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[507] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1207] += dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 112: C2H6 + OH => C2H5 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[17];
    k_f = 1.0000000000000002e-06 * 58000000
                * exp(1.73 * tc[0] - 0.50321666580471969 * (1159.8900000000001) * invT);
    dlnkfdT = 1.73 * invT + 0.50321666580471969 *  (1159.8900000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[14] += q; /* C2H5 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[17];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[119] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[122] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[3];
    J[598] -= dqdci;              /* dwdot[OH]/d[C2H6] */
    J[601] += dqdci;              /* dwdot[H2O]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 113: CH2GSG + C2H6 => CH3 + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[17];
    k_f = 1.0000000000000002e-06 * 120000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] += q; /* CH3 */
    wdot[14] += q; /* C2H5 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[17];
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] += dqdci;              /* dwdot[CH3]/d[CH2GSG] */
    J[329] += dqdci;              /* dwdot[C2H5]/d[CH2GSG] */
    J[332] -= dqdci;              /* dwdot[C2H6]/d[CH2GSG] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[9];
    J[604] -= dqdci;              /* dwdot[CH2GSG]/d[C2H6] */
    J[606] += dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 114: C2H6 + O => C2H5 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[17];
    k_f = 1.0000000000000002e-06 * 13000000
                * exp(2.1299999999999999 * tc[0] - 0.50321666580471969 * (5190.0100000000002) * invT);
    dlnkfdT = 2.1299999999999999 * invT + 0.50321666580471969 *  (5190.0100000000002)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[14] += q; /* C2H5 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[49] += dqdci;               /* dwdot[C2H5]/d[O] */
    J[52] -= dqdci;               /* dwdot[C2H6]/d[O] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[1];
    J[596] -= dqdci;              /* dwdot[O]/d[C2H6] */
    J[598] += dqdci;              /* dwdot[OH]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 115: C2H6 + CH3 => C2H5 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[17];
    k_f = 1.0000000000000002e-06 * 1.5099999999999999e-07
                * exp(6 * tc[0] - 0.50321666580471969 * (6047.0799999999999) * invT);
    dlnkfdT = 6 * invT + 0.50321666580471969 *  (6047.0799999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[14] += q; /* C2H5 */
    wdot[15] += q; /* CH4 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[17];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[399] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[402] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[11];
    J[606] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[610] += dqdci;              /* dwdot[CH4]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 116: HCCO + O => H + 2.000000 CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[23];
    k_f = 1.0000000000000002e-06 * 80000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[12] += 2 * q; /* CO */
    wdot[23] -= q; /* HCCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[23];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[47] += 2 * dqdci;           /* dwdot[CO]/d[O] */
    J[58] -= dqdci;               /* dwdot[HCCO]/d[O] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[1];
    J[806] -= dqdci;              /* dwdot[O]/d[HCCO] */
    J[807] += dqdci;              /* dwdot[H]/d[HCCO] */
    J[817] += 2 * dqdci;          /* dwdot[CO]/d[HCCO] */
    J[828] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1202] += 2 * dqdT;          /* dwdot[CO]/dT */
    J[1213] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 117: HCCO + O2 => CO2 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[23];
    k_f = 1.0000000000000002e-06 * 240000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-853.97000000000003) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-853.97000000000003)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[18] += q; /* CO2 */
    wdot[19] += q; /* HCO */
    wdot[23] -= q; /* HCCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[23];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[193] += dqdci;              /* dwdot[CO2]/d[O2] */
    J[194] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[198] -= dqdci;              /* dwdot[HCCO]/d[O2] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[5];
    J[810] -= dqdci;              /* dwdot[O2]/d[HCCO] */
    J[823] += dqdci;              /* dwdot[CO2]/d[HCCO] */
    J[824] += dqdci;              /* dwdot[HCO]/d[HCCO] */
    J[828] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1213] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 118: HCCO + OH => 2.000000 HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[23];
    k_f = 1.0000000000000002e-06 * 10000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[19] += 2 * q; /* HCO */
    wdot[23] -= q; /* HCCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[23];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[124] += 2 * dqdci;          /* dwdot[HCO]/d[OH] */
    J[128] -= dqdci;              /* dwdot[HCCO]/d[OH] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[3];
    J[808] -= dqdci;              /* dwdot[OH]/d[HCCO] */
    J[824] += 2 * dqdci;          /* dwdot[HCO]/d[HCCO] */
    J[828] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1209] += 2 * dqdT;          /* dwdot[HCO]/dT */
    J[1213] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 119: HCCO + H => CH2GSG + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[23];
    k_f = 1.0000000000000002e-06 * 110000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[9] += q; /* CH2GSG */
    wdot[12] += q; /* CO */
    wdot[23] -= q; /* HCCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[23];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[79] += dqdci;               /* dwdot[CH2GSG]/d[H] */
    J[82] += dqdci;               /* dwdot[CO]/d[H] */
    J[93] -= dqdci;               /* dwdot[HCCO]/d[H] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[2];
    J[807] -= dqdci;              /* dwdot[H]/d[HCCO] */
    J[814] += dqdci;              /* dwdot[CH2GSG]/d[HCCO] */
    J[817] += dqdci;              /* dwdot[CO]/d[HCCO] */
    J[828] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1199] += dqdT;              /* dwdot[CH2GSG]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1213] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 120: CH2GSG + CO => HCCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[12];
    k_f = 1.0000000000000002e-06 * 2046000000000
                * exp(0.89000000000000001 * tc[0] - 0.50321666580471969 * (27830.07) * invT);
    dlnkfdT = 0.89000000000000001 * invT + 0.50321666580471969 *  (27830.07)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[9] -= q; /* CH2GSG */
    wdot[12] -= q; /* CO */
    wdot[23] += q; /* HCCO */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[12];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[327] -= dqdci;              /* dwdot[CO]/d[CH2GSG] */
    J[338] += dqdci;              /* dwdot[HCCO]/d[CH2GSG] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[9];
    J[422] += dqdci;              /* dwdot[H]/d[CO] */
    J[429] -= dqdci;              /* dwdot[CH2GSG]/d[CO] */
    J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[443] += dqdci;              /* dwdot[HCCO]/d[CO] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1202] -= dqdT;              /* dwdot[CO]/dT */
    J[1213] += dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 121: CH2CHO + O2 => CH2O + CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[25];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (4200.0500000000002) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (4200.0500000000002)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[12] += q; /* CO */
    wdot[25] -= q; /* CH2CHO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[25];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[187] += dqdci;              /* dwdot[CO]/d[O2] */
    J[200] -= dqdci;              /* dwdot[CH2CHO]/d[O2] */
    /* d()/d[CH2CHO] */
    dqdci =  + k_f*sc[5];
    J[878] += dqdci;              /* dwdot[OH]/d[CH2CHO] */
    J[880] -= dqdci;              /* dwdot[O2]/d[CH2CHO] */
    J[885] += dqdci;              /* dwdot[CH2O]/d[CH2CHO] */
    J[887] += dqdci;              /* dwdot[CO]/d[CH2CHO] */
    J[900] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1215] -= dqdT;              /* dwdot[CH2CHO]/dT */

    /*reaction 122: C3H5XA + CH2O => C3H6 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[26];
    k_f = 1.0000000000000002e-06 * 630000000
                * exp(1.8999999999999999 * tc[0] - 0.50321666580471969 * (18190.009999999998) * invT);
    dlnkfdT = 1.8999999999999999 * invT + 0.50321666580471969 *  (18190.009999999998)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    wdot[26] -= q; /* C3H5XA */
    wdot[27] += q; /* C3H6 */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[26];
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[376] -= dqdci;              /* dwdot[C3H5XA]/d[CH2O] */
    J[377] += dqdci;              /* dwdot[C3H6]/d[CH2O] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[10];
    J[920] -= dqdci;              /* dwdot[CH2O]/d[C3H5XA] */
    J[929] += dqdci;              /* dwdot[HCO]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[937] += dqdci;              /* dwdot[C3H6]/d[C3H5XA] */
    /* d()/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 123: C3H6 + HCO => C3H5XA + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[19]*sc[27];
    k_f = 1.0000000000000002e-06 * 109700000
                * exp(1.8899999999999999 * tc[0] - 0.50321666580471969 * (15840.110000000001) * invT);
    dlnkfdT = 1.8899999999999999 * invT + 0.50321666580471969 *  (15840.110000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] += q; /* CH2O */
    wdot[19] -= q; /* HCO */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[27];
    J[675] += dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[691] += dqdci;              /* dwdot[C3H5XA]/d[HCO] */
    J[692] -= dqdci;              /* dwdot[C3H6]/d[HCO] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[19];
    J[955] += dqdci;              /* dwdot[CH2O]/d[C3H6] */
    J[964] -= dqdci;              /* dwdot[HCO]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 124: C3H5XA + HO2 => C3H5O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[26];
    k_f = 1.0000000000000002e-06 * 7000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1000)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[26] -= q; /* C3H5XA */
    wdot[28] += q; /* C3H5O */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[26];
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[306] -= dqdci;              /* dwdot[C3H5XA]/d[HO2] */
    J[308] += dqdci;              /* dwdot[C3H5O]/d[HO2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[8];
    J[913] += dqdci;              /* dwdot[OH]/d[C3H5XA] */
    J[918] -= dqdci;              /* dwdot[HO2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[938] += dqdci;              /* dwdot[C3H5O]/d[C3H5XA] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1218] += dqdT;              /* dwdot[C3H5O]/dT */

    /*reaction 125: C3H5XA + CH3O2 => C3H5O + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[20]*sc[26];
    k_f = 1.0000000000000002e-06 * 7000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1000)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[13] += q; /* CH3O */
    wdot[20] -= q; /* CH3O2 */
    wdot[26] -= q; /* C3H5XA */
    wdot[28] += q; /* C3H5O */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[26];
    J[713] += dqdci;              /* dwdot[CH3O]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    J[726] -= dqdci;              /* dwdot[C3H5XA]/d[CH3O2] */
    J[728] += dqdci;              /* dwdot[C3H5O]/d[CH3O2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[20];
    J[923] += dqdci;              /* dwdot[CH3O]/d[C3H5XA] */
    J[930] -= dqdci;              /* dwdot[CH3O2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[938] += dqdci;              /* dwdot[C3H5O]/d[C3H5XA] */
    /* d()/dT */
    J[1203] += dqdT;              /* dwdot[CH3O]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1218] += dqdT;              /* dwdot[C3H5O]/dT */

    /*reaction 126: C3H5XA + HO2 => C3H6 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[26];
    k_f = 1.0000000000000002e-06 * 33320000000
                * exp(0.34000000000000002 * tc[0] - 0.50321666580471969 * (-555.92999999999995) * invT);
    dlnkfdT = 0.34000000000000002 * invT + 0.50321666580471969 *  (-555.92999999999995)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[26] -= q; /* C3H5XA */
    wdot[27] += q; /* C3H6 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[26];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[306] -= dqdci;              /* dwdot[C3H5XA]/d[HO2] */
    J[307] += dqdci;              /* dwdot[C3H6]/d[HO2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[8];
    J[915] += dqdci;              /* dwdot[O2]/d[C3H5XA] */
    J[918] -= dqdci;              /* dwdot[HO2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[937] += dqdci;              /* dwdot[C3H6]/d[C3H5XA] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 127: C3H5XA + O2 => CH2CHO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[26];
    k_f = 1.0000000000000002e-06 * 7140000000000000
                * exp(-1.21 * tc[0] - 0.50321666580471969 * (21049.950000000001) * invT);
    dlnkfdT = -1.21 * invT + 0.50321666580471969 *  (21049.950000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[25] += q; /* CH2CHO */
    wdot[26] -= q; /* C3H5XA */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[26];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[200] += dqdci;              /* dwdot[CH2CHO]/d[O2] */
    J[201] -= dqdci;              /* dwdot[C3H5XA]/d[O2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[5];
    J[915] -= dqdci;              /* dwdot[O2]/d[C3H5XA] */
    J[920] += dqdci;              /* dwdot[CH2O]/d[C3H5XA] */
    J[935] += dqdci;              /* dwdot[CH2CHO]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1215] += dqdT;              /* dwdot[CH2CHO]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */

    /*reaction 128: C3H5XA + O2 => C2H2 + CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[26];
    k_f = 1.0000000000000002e-06 * 9.7200000000000003e+29
                * exp(-5.71 * tc[0] - 0.50321666580471969 * (21450.049999999999) * invT);
    dlnkfdT = -5.71 * invT + 0.50321666580471969 *  (21450.049999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[22] += q; /* C2H2 */
    wdot[26] -= q; /* C3H5XA */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[26];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[197] += dqdci;              /* dwdot[C2H2]/d[O2] */
    J[201] -= dqdci;              /* dwdot[C3H5XA]/d[O2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[5];
    J[913] += dqdci;              /* dwdot[OH]/d[C3H5XA] */
    J[915] -= dqdci;              /* dwdot[O2]/d[C3H5XA] */
    J[920] += dqdci;              /* dwdot[CH2O]/d[C3H5XA] */
    J[932] += dqdci;              /* dwdot[C2H2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1212] += dqdT;              /* dwdot[C2H2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */

    /*reaction 129: C3H5XA + H => C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[26];
    k_f = 1.0000000000000002e-06 * 4.8869999999999999e+56
                * exp(-12.25 * tc[0] - 0.50321666580471969 * (28080.07) * invT);
    dlnkfdT = -12.25 * invT + 0.50321666580471969 *  (28080.07)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[26] -= q; /* C3H5XA */
    wdot[27] += q; /* C3H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[26];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[96] -= dqdci;               /* dwdot[C3H5XA]/d[H] */
    J[97] += dqdci;               /* dwdot[C3H6]/d[H] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[2];
    J[912] -= dqdci;              /* dwdot[H]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[937] += dqdci;              /* dwdot[C3H6]/d[C3H5XA] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 130: C3H5XA => C2H2 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[26];
    k_f = 1 * 2.397e+48
                * exp(-9.9000000000000004 * tc[0] - 0.50321666580471969 * (82080.070000000007) * invT);
    dlnkfdT = -9.9000000000000004 * invT + 0.50321666580471969 *  (82080.070000000007)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] += q; /* CH3 */
    wdot[22] += q; /* C2H2 */
    wdot[26] -= q; /* C3H5XA */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f;
    J[921] += dqdci;              /* dwdot[CH3]/d[C3H5XA] */
    J[932] += dqdci;              /* dwdot[C2H2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    /* d()/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1212] += dqdT;              /* dwdot[C2H2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */

    /*reaction 131: C3H6 => C2H3 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[27];
    k_f = 1 * 2.7300000000000001e+62
                * exp(-13.279999999999999 * tc[0] - 0.50321666580471969 * (123200.05) * invT);
    dlnkfdT = -13.279999999999999 * invT + 0.50321666580471969 *  (123200.05)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] += q; /* CH3 */
    wdot[24] += q; /* C2H3 */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[C3H6] */
    dqdci =  + k_f;
    J[956] += dqdci;              /* dwdot[CH3]/d[C3H6] */
    J[969] += dqdci;              /* dwdot[C2H3]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 132: C2H3 + CH3 => C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[24];
    k_f = 1.0000000000000002e-06 * 4.7119999999999996e+59
                * exp(-13.19 * tc[0] - 0.50321666580471969 * (29539.91) * invT);
    dlnkfdT = -13.19 * invT + 0.50321666580471969 *  (29539.91)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[24] -= q; /* C2H3 */
    wdot[27] += q; /* C3H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[24];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] -= dqdci;              /* dwdot[C2H3]/d[CH3] */
    J[412] += dqdci;              /* dwdot[C3H6]/d[CH3] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[11];
    J[851] -= dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[867] += dqdci;              /* dwdot[C3H6]/d[C2H3] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 133: C3H6 + O => C2H5 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[27];
    k_f = 1.0000000000000002e-06 * 15800000
                * exp(1.76 * tc[0] - 0.50321666580471969 * (-1216.0599999999999) * invT);
    dlnkfdT = 1.76 * invT + 0.50321666580471969 *  (-1216.0599999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[14] += q; /* C2H5 */
    wdot[19] += q; /* HCO */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[27];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[49] += dqdci;               /* dwdot[C2H5]/d[O] */
    J[54] += dqdci;               /* dwdot[HCO]/d[O] */
    J[62] -= dqdci;               /* dwdot[C3H6]/d[O] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[1];
    J[946] -= dqdci;              /* dwdot[O]/d[C3H6] */
    J[959] += dqdci;              /* dwdot[C2H5]/d[C3H6] */
    J[964] += dqdci;              /* dwdot[HCO]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 134: C3H6 + H => C2H4 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = 1.0000000000000002e-06 * 4.8299999999999998e+33
                * exp(-5.8099999999999996 * tc[0] - 0.50321666580471969 * (18500) * invT);
    dlnkfdT = -5.8099999999999996 * invT + 0.50321666580471969 *  (18500)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[11] += q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[81] += dqdci;               /* dwdot[CH3]/d[H] */
    J[86] += dqdci;               /* dwdot[C2H4]/d[H] */
    J[97] -= dqdci;               /* dwdot[C3H6]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[2];
    J[947] -= dqdci;              /* dwdot[H]/d[C3H6] */
    J[956] += dqdci;              /* dwdot[CH3]/d[C3H6] */
    J[961] += dqdci;              /* dwdot[C2H4]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 135: C2H4 + CH3 => C3H6 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[16];
    k_f = 1.0000000000000002e-06 * 2.3130000000000001e+33
                * exp(-5.9000000000000004 * tc[0] - 0.50321666580471969 * (31619.98) * invT);
    dlnkfdT = -5.9000000000000004 * invT + 0.50321666580471969 *  (31619.98)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[11] -= q; /* CH3 */
    wdot[16] -= q; /* C2H4 */
    wdot[27] += q; /* C3H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[401] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[412] += dqdci;              /* dwdot[C3H6]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[11];
    J[562] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[571] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[587] += dqdci;              /* dwdot[C3H6]/d[C2H4] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 136: C3H6 + H => C3H5XA + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = 1.0000000000000002e-06 * 173000
                * exp(2.5 * tc[0] - 0.50321666580471969 * (2492.1100000000001) * invT);
    dlnkfdT = 2.5 * invT + 0.50321666580471969 *  (2492.1100000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[96] += dqdci;               /* dwdot[C3H5XA]/d[H] */
    J[97] -= dqdci;               /* dwdot[C3H6]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[2];
    J[947] -= dqdci;              /* dwdot[H]/d[C3H6] */
    J[949] += dqdci;              /* dwdot[H2]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 137: C3H5XA + H2 => C3H6 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[26];
    k_f = 1.0000000000000002e-06 * 79330
                * exp(2.5099999999999998 * tc[0] - 0.50321666580471969 * (19520.080000000002) * invT);
    dlnkfdT = 2.5099999999999998 * invT + 0.50321666580471969 *  (19520.080000000002)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[26] -= q; /* C3H5XA */
    wdot[27] += q; /* C3H6 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[26];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[166] -= dqdci;              /* dwdot[C3H5XA]/d[H2] */
    J[167] += dqdci;              /* dwdot[C3H6]/d[H2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[4];
    J[912] += dqdci;              /* dwdot[H]/d[C3H5XA] */
    J[914] -= dqdci;              /* dwdot[H2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[937] += dqdci;              /* dwdot[C3H6]/d[C3H5XA] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 138: C3H6 + OH => C3H5XA + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[27];
    k_f = 1.0000000000000002e-06 * 3120000
                * exp(2 * tc[0] - 0.50321666580471969 * (-298.04000000000002) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (-298.04000000000002)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[27];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[131] += dqdci;              /* dwdot[C3H5XA]/d[OH] */
    J[132] -= dqdci;              /* dwdot[C3H6]/d[OH] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[3];
    J[948] -= dqdci;              /* dwdot[OH]/d[C3H6] */
    J[951] += dqdci;              /* dwdot[H2O]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 139: C3H6 + O => C3H5XA + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[27];
    k_f = 1.0000000000000002e-06 * 524000000000
                * exp(0.69999999999999996 * tc[0] - 0.50321666580471969 * (5884.0799999999999) * invT);
    dlnkfdT = 0.69999999999999996 * invT + 0.50321666580471969 *  (5884.0799999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[27];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[61] += dqdci;               /* dwdot[C3H5XA]/d[O] */
    J[62] -= dqdci;               /* dwdot[C3H6]/d[O] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[1];
    J[946] -= dqdci;              /* dwdot[O]/d[C3H6] */
    J[948] += dqdci;              /* dwdot[OH]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 140: C3H6 + CH3 => C3H5XA + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[27];
    k_f = 1.0000000000000002e-06 * 2.21
                * exp(3.5 * tc[0] - 0.50321666580471969 * (5674.9499999999998) * invT);
    dlnkfdT = 3.5 * invT + 0.50321666580471969 *  (5674.9499999999998)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[27];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[411] += dqdci;              /* dwdot[C3H5XA]/d[CH3] */
    J[412] -= dqdci;              /* dwdot[C3H6]/d[CH3] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[11];
    J[956] -= dqdci;              /* dwdot[CH3]/d[C3H6] */
    J[960] += dqdci;              /* dwdot[CH4]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 141: IXC3H7 + O2 => C3H6 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[29];
    k_f = 1.0000000000000002e-06 * 450000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (5020.0799999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (5020.0799999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[27] += q; /* C3H6 */
    wdot[29] -= q; /* IXC3H7 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[29];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[202] += dqdci;              /* dwdot[C3H6]/d[O2] */
    J[204] -= dqdci;              /* dwdot[IXC3H7]/d[O2] */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f*sc[5];
    J[1020] -= dqdci;             /* dwdot[O2]/d[IXC3H7] */
    J[1023] += dqdci;             /* dwdot[HO2]/d[IXC3H7] */
    J[1042] += dqdci;             /* dwdot[C3H6]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */

    /*reaction 142: IXC3H7 => H + C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[29];
    k_f = 1 * 8.569e+18
                * exp(-1.5700000000000001 * tc[0] - 0.50321666580471969 * (40340.110000000001) * invT);
    dlnkfdT = -1.5700000000000001 * invT + 0.50321666580471969 *  (40340.110000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[27] += q; /* C3H6 */
    wdot[29] -= q; /* IXC3H7 */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f;
    J[1017] += dqdci;             /* dwdot[H]/d[IXC3H7] */
    J[1042] += dqdci;             /* dwdot[C3H6]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */

    /*reaction 143: H + C3H6 => IXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = 1.0000000000000002e-06 * 13000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (1559.99) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (1559.99)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[27] -= q; /* C3H6 */
    wdot[29] += q; /* IXC3H7 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[97] -= dqdci;               /* dwdot[C3H6]/d[H] */
    J[99] += dqdci;               /* dwdot[IXC3H7]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[2];
    J[947] -= dqdci;              /* dwdot[H]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    J[974] += dqdci;              /* dwdot[IXC3H7]/d[C3H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */

    /*reaction 144: NXC3H7 => CH3 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[30];
    k_f = 1 * 228400000000000
                * exp(-0.55000000000000004 * tc[0] - 0.50321666580471969 * (28400.099999999999) * invT);
    dlnkfdT = -0.55000000000000004 * invT + 0.50321666580471969 *  (28400.099999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] += q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    wdot[30] -= q; /* NXC3H7 */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f;
    J[1061] += dqdci;             /* dwdot[CH3]/d[NXC3H7] */
    J[1066] += dqdci;             /* dwdot[C2H4]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    /* d()/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 145: CH3 + C2H4 => NXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[16];
    k_f = 1.0000000000000002e-06 * 410000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7204.1099999999997) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (7204.1099999999997)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[16] -= q; /* C2H4 */
    wdot[30] += q; /* NXC3H7 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[401] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[415] += dqdci;              /* dwdot[NXC3H7]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[11];
    J[571] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[590] += dqdci;              /* dwdot[NXC3H7]/d[C2H4] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 146: NXC3H7 + HO2 => C3H8 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[30];
    k_f = 1.0000000000000002e-06 * 2080000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[30] -= q; /* NXC3H7 */
    wdot[31] += q; /* C3H8 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[30];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[310] -= dqdci;              /* dwdot[NXC3H7]/d[HO2] */
    J[311] += dqdci;              /* dwdot[C3H8]/d[HO2] */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f*sc[8];
    J[1055] += dqdci;             /* dwdot[O2]/d[NXC3H7] */
    J[1058] -= dqdci;             /* dwdot[HO2]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    J[1081] += dqdci;             /* dwdot[C3H8]/d[NXC3H7] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] += dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 147: NXC3H7 + O2 => C3H6 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[30];
    k_f = 1.0000000000000002e-06 * 300000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (3000)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[27] += q; /* C3H6 */
    wdot[30] -= q; /* NXC3H7 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[30];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[202] += dqdci;              /* dwdot[C3H6]/d[O2] */
    J[205] -= dqdci;              /* dwdot[NXC3H7]/d[O2] */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f*sc[5];
    J[1055] -= dqdci;             /* dwdot[O2]/d[NXC3H7] */
    J[1058] += dqdci;             /* dwdot[HO2]/d[NXC3H7] */
    J[1077] += dqdci;             /* dwdot[C3H6]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 148: NXC3H7 => H + C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[30];
    k_f = 1 * 2667000000000000
                * exp(-0.64000000000000001 * tc[0] - 0.50321666580471969 * (36820.029999999999) * invT);
    dlnkfdT = -0.64000000000000001 * invT + 0.50321666580471969 *  (36820.029999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[27] += q; /* C3H6 */
    wdot[30] -= q; /* NXC3H7 */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f;
    J[1052] += dqdci;             /* dwdot[H]/d[NXC3H7] */
    J[1077] += dqdci;             /* dwdot[C3H6]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 149: H + C3H6 => NXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = 1.0000000000000002e-06 * 10000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (2500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (2500)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[27] -= q; /* C3H6 */
    wdot[30] += q; /* NXC3H7 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[97] -= dqdci;               /* dwdot[C3H6]/d[H] */
    J[100] += dqdci;              /* dwdot[NXC3H7]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[2];
    J[947] -= dqdci;              /* dwdot[H]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    J[975] += dqdci;              /* dwdot[NXC3H7]/d[C3H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 150: C3H8 + OH => NXC3H7 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[31];
    k_f = 1.0000000000000002e-06 * 10540000000
                * exp(0.96999999999999997 * tc[0] - 0.50321666580471969 * (1586.04) * invT);
    dlnkfdT = 0.96999999999999997 * invT + 0.50321666580471969 *  (1586.04)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[31];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[135] += dqdci;              /* dwdot[NXC3H7]/d[OH] */
    J[136] -= dqdci;              /* dwdot[C3H8]/d[OH] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[3];
    J[1088] -= dqdci;             /* dwdot[OH]/d[C3H8] */
    J[1091] += dqdci;             /* dwdot[H2O]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 151: C3H8 + HO2 => NXC3H7 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[31];
    k_f = 1.0000000000000002e-06 * 16800000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (20429.970000000001) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (20429.970000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[31];
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[310] += dqdci;              /* dwdot[NXC3H7]/d[HO2] */
    J[311] -= dqdci;              /* dwdot[C3H8]/d[HO2] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[8];
    J[1092] += dqdci;             /* dwdot[H2O2]/d[C3H8] */
    J[1093] -= dqdci;             /* dwdot[HO2]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 152: H + C3H8 <=> H2 + NXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[31];
    k_f = 1.0000000000000002e-06 * 3972000
                * exp(2.75 * tc[0] - 0.50321666580471969 * (6756.6899999999996) * invT);
    dlnkfdT = 2.75 * invT + 0.50321666580471969 *  (6756.6899999999996)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[30];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[30] + g_RT[31]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[31]) + (h_RT[4] + h_RT[30]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[31];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[100] += dqdci;              /* dwdot[NXC3H7]/d[H] */
    J[101] -= dqdci;              /* dwdot[C3H8]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[30];
    J[142] -= dqdci;              /* dwdot[H]/d[H2] */
    J[144] += dqdci;              /* dwdot[H2]/d[H2] */
    J[170] += dqdci;              /* dwdot[NXC3H7]/d[H2] */
    J[171] -= dqdci;              /* dwdot[C3H8]/d[H2] */
    /* d()/d[NXC3H7] */
    dqdci =  - k_r*sc[4];
    J[1052] -= dqdci;             /* dwdot[H]/d[NXC3H7] */
    J[1054] += dqdci;             /* dwdot[H2]/d[NXC3H7] */
    J[1080] += dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    J[1081] -= dqdci;             /* dwdot[C3H8]/d[NXC3H7] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[2];
    J[1087] -= dqdci;             /* dwdot[H]/d[C3H8] */
    J[1089] += dqdci;             /* dwdot[H2]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 153: C3H8 + OH => IXC3H7 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[31];
    k_f = 1.0000000000000002e-06 * 46700000
                * exp(1.6100000000000001 * tc[0] - 0.50321666580471969 * (-34.890000000000001) * invT);
    dlnkfdT = 1.6100000000000001 * invT + 0.50321666580471969 *  (-34.890000000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[31];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[134] += dqdci;              /* dwdot[IXC3H7]/d[OH] */
    J[136] -= dqdci;              /* dwdot[C3H8]/d[OH] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[3];
    J[1088] -= dqdci;             /* dwdot[OH]/d[C3H8] */
    J[1091] += dqdci;             /* dwdot[H2O]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 154: CH3 + C3H8 => CH4 + IXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[31];
    k_f = 1.0000000000000002e-06 * 398000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (9500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (9500)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[31];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[414] += dqdci;              /* dwdot[IXC3H7]/d[CH3] */
    J[416] -= dqdci;              /* dwdot[C3H8]/d[CH3] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[11];
    J[1096] -= dqdci;             /* dwdot[CH3]/d[C3H8] */
    J[1100] += dqdci;             /* dwdot[CH4]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 155: CH3 + C3H8 => CH4 + NXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[31];
    k_f = 1.0000000000000002e-06 * 1290000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (11599.9) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (11599.9)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[31];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[415] += dqdci;              /* dwdot[NXC3H7]/d[CH3] */
    J[416] -= dqdci;              /* dwdot[C3H8]/d[CH3] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[11];
    J[1096] -= dqdci;             /* dwdot[CH3]/d[C3H8] */
    J[1100] += dqdci;             /* dwdot[CH4]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 156: C3H8 + O => IXC3H7 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[31];
    k_f = 1.0000000000000002e-06 * 28100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (5200.0500000000002) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (5200.0500000000002)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[31];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[64] += dqdci;               /* dwdot[IXC3H7]/d[O] */
    J[66] -= dqdci;               /* dwdot[C3H8]/d[O] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[1];
    J[1086] -= dqdci;             /* dwdot[O]/d[C3H8] */
    J[1088] += dqdci;             /* dwdot[OH]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 157: C3H8 + HO2 => IXC3H7 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[31];
    k_f = 1.0000000000000002e-06 * 5600000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (17700.049999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (17700.049999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[31];
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[309] += dqdci;              /* dwdot[IXC3H7]/d[HO2] */
    J[311] -= dqdci;              /* dwdot[C3H8]/d[HO2] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[8];
    J[1092] += dqdci;             /* dwdot[H2O2]/d[C3H8] */
    J[1093] -= dqdci;             /* dwdot[HO2]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 158: C3H8 + O => NXC3H7 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[31];
    k_f = 1.0000000000000002e-06 * 113000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7849.8999999999996) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (7849.8999999999996)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[31];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[65] += dqdci;               /* dwdot[NXC3H7]/d[O] */
    J[66] -= dqdci;               /* dwdot[C3H8]/d[O] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[1];
    J[1086] -= dqdci;             /* dwdot[O]/d[C3H8] */
    J[1088] += dqdci;             /* dwdot[OH]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 159: C3H8 + O2 => IXC3H7 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[31];
    k_f = 1.0000000000000002e-06 * 40000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (47500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (47500)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[31];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[204] += dqdci;              /* dwdot[IXC3H7]/d[O2] */
    J[206] -= dqdci;              /* dwdot[C3H8]/d[O2] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[5];
    J[1090] -= dqdci;             /* dwdot[O2]/d[C3H8] */
    J[1093] += dqdci;             /* dwdot[HO2]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 160: IXC3H7 + HO2 => C3H8 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[29];
    k_f = 1.0000000000000002e-06 * 2080000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[29] -= q; /* IXC3H7 */
    wdot[31] += q; /* C3H8 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[29];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[309] -= dqdci;              /* dwdot[IXC3H7]/d[HO2] */
    J[311] += dqdci;              /* dwdot[C3H8]/d[HO2] */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f*sc[8];
    J[1020] += dqdci;             /* dwdot[O2]/d[IXC3H7] */
    J[1023] -= dqdci;             /* dwdot[HO2]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    J[1046] += dqdci;             /* dwdot[C3H8]/d[IXC3H7] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] += dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 161: H + C3H8 => H2 + IXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[31];
    k_f = 1.0000000000000002e-06 * 1300000
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (4471.0799999999999) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  (4471.0799999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[31];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[99] += dqdci;               /* dwdot[IXC3H7]/d[H] */
    J[101] -= dqdci;              /* dwdot[C3H8]/d[H] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[2];
    J[1087] -= dqdci;             /* dwdot[H]/d[C3H8] */
    J[1089] += dqdci;             /* dwdot[H2]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 162: H2 + IXC3H7 => H + C3H8 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[29];
    k_f = 1.0000000000000002e-06 * 470900
                * exp(2.1499999999999999 * tc[0] - 0.50321666580471969 * (12179.969999999999) * invT);
    dlnkfdT = 2.1499999999999999 * invT + 0.50321666580471969 *  (12179.969999999999)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[29] -= q; /* IXC3H7 */
    wdot[31] += q; /* C3H8 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[29];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[169] -= dqdci;              /* dwdot[IXC3H7]/d[H2] */
    J[171] += dqdci;              /* dwdot[C3H8]/d[H2] */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f*sc[4];
    J[1017] += dqdci;             /* dwdot[H]/d[IXC3H7] */
    J[1019] -= dqdci;             /* dwdot[H2]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    J[1046] += dqdci;             /* dwdot[C3H8]/d[IXC3H7] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] += dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 163: C3H5O => C2H3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[28];
    k_f = 1 * 2028000000000
                * exp(0.089999999999999997 * tc[0] - 0.50321666580471969 * (23559.990000000002) * invT);
    dlnkfdT = 0.089999999999999997 * invT + 0.50321666580471969 *  (23559.990000000002)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] += q; /* CH2O */
    wdot[24] += q; /* C2H3 */
    wdot[28] -= q; /* C3H5O */
    /* d()/d[C3H5O] */
    dqdci =  + k_f;
    J[990] += dqdci;              /* dwdot[CH2O]/d[C3H5O] */
    J[1004] += dqdci;             /* dwdot[C2H3]/d[C3H5O] */
    J[1008] -= dqdci;             /* dwdot[C3H5O]/d[C3H5O] */
    /* d()/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */
    J[1218] -= dqdT;              /* dwdot[C3H5O]/dT */

    /*reaction 164: IXC3H7O2 => IXC3H7 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[32];
    k_f = 1 * 2.803e+17
                * exp(-0.62 * tc[0] - 0.50321666580471969 * (36039.910000000003) * invT);
    dlnkfdT = -0.62 * invT + 0.50321666580471969 *  (36039.910000000003)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[29] += q; /* IXC3H7 */
    wdot[32] -= q; /* IXC3H7O2 */
    /* d()/d[IXC3H7O2] */
    dqdci =  + k_f;
    J[1125] += dqdci;             /* dwdot[O2]/d[IXC3H7O2] */
    J[1149] += dqdci;             /* dwdot[IXC3H7]/d[IXC3H7O2] */
    J[1152] -= dqdci;             /* dwdot[IXC3H7O2]/d[IXC3H7O2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1222] -= dqdT;              /* dwdot[IXC3H7O2]/dT */

    /*reaction 165: IXC3H7 + O2 => IXC3H7O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[29];
    k_f = 1.0000000000000002e-06 * 7540000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[29] -= q; /* IXC3H7 */
    wdot[32] += q; /* IXC3H7O2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[29];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[204] -= dqdci;              /* dwdot[IXC3H7]/d[O2] */
    J[207] += dqdci;              /* dwdot[IXC3H7O2]/d[O2] */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f*sc[5];
    J[1020] -= dqdci;             /* dwdot[O2]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    J[1047] += dqdci;             /* dwdot[IXC3H7O2]/d[IXC3H7] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */
    J[1222] += dqdT;              /* dwdot[IXC3H7O2]/dT */

    /*reaction 166: NXC3H7O2 => NXC3H7 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[33];
    k_f = 1 * 3.364e+19
                * exp(-1.3200000000000001 * tc[0] - 0.50321666580471969 * (35760.040000000001) * invT);
    dlnkfdT = -1.3200000000000001 * invT + 0.50321666580471969 *  (35760.040000000001)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[30] += q; /* NXC3H7 */
    wdot[33] -= q; /* NXC3H7O2 */
    /* d()/d[NXC3H7O2] */
    dqdci =  + k_f;
    J[1160] += dqdci;             /* dwdot[O2]/d[NXC3H7O2] */
    J[1185] += dqdci;             /* dwdot[NXC3H7]/d[NXC3H7O2] */
    J[1188] -= dqdci;             /* dwdot[NXC3H7O2]/d[NXC3H7O2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1223] -= dqdT;              /* dwdot[NXC3H7O2]/dT */

    /*reaction 167: NXC3H7 + O2 => NXC3H7O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[30];
    k_f = 1.0000000000000002e-06 * 4520000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[30] -= q; /* NXC3H7 */
    wdot[33] += q; /* NXC3H7O2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[30];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[205] -= dqdci;              /* dwdot[NXC3H7]/d[O2] */
    J[208] += dqdci;              /* dwdot[NXC3H7O2]/d[O2] */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f*sc[5];
    J[1055] -= dqdci;             /* dwdot[O2]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    J[1083] += dqdci;             /* dwdot[NXC3H7O2]/d[NXC3H7] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */
    J[1223] += dqdT;              /* dwdot[NXC3H7O2]/dT */

    amrex::Real c_R[34], dcRdT[34], e_RT[34];
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
    for (int k = 0; k < 34; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[1190+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 34; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 34; ++m) {
            dehmixdc += eh_RT[m]*J[k*35+m];
        }
        J[k*35+34] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[1224] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

/*compute an approx to the reaction Jacobian (for preconditioning) */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void DWDOT_SIMPLIFIED(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * HP)
{
    amrex::Real c[34];

    for (int k=0; k<34; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<34; k++) {
        J[1190+k] *= 1.e-6;
        J[k*35+34] *= 1.e6;
    }

    return;
}

#ifdef AMREX_USE_GPU
/*compute the reaction Jacobian on GPU */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void aJacobian(amrex::Real * J, amrex::Real * sc, amrex::Real T, int consP)
{

    for (int i=0; i<1225; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[34];
    for (int k=0; k<34; k++) {
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
    for (int k = 0; k < 34; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[34];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[34];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[34];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[3], 2.000000);
    k_f = 1.0000000000000002e-06 * 123600000000000
                * exp(-0.37 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.37 * invT + 0.50321666580471969 *  0  * invT2;
    /* pressure-fall-off */
    k_0 = 3.0410000000000002e+30 * exp(-4.6299999999999999 * tc[0] - 0.50321666580471969 * (2049) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -4.6299999999999999 * invT + 0.50321666580471969 * (2049) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.46999999999999997)*exp(-T/100);
    Fcent2 = 0.46999999999999997 * exp(-T/2000);
    Fcent3 = exp(-1000000000000000 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/100
        -Fcent2/2000
        + Fcent3*1000000000000000*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[7];
    Kc = refCinv * exp(2.000000*g_RT[3] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[3]) + (h_RT[7]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[3] -= 2 * q; /* OH */
    wdot[7] += q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[OH] */
        dqdci =  + k_f*2.000000*sc[3];
        J[108] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
        J[112] += dqdci;              /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*dcdc_fac;
        J[143] += -2 * dqdci;         /* dwdot[OH]/d[H2] */
        J[147] += dqdci;              /* dwdot[H2O2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*dcdc_fac;
        J[213] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
        J[217] += dqdci;              /* dwdot[H2O2]/d[H2O] */
        /* d()/d[H2O2] */
        dqdci =  - k_r;
        J[248] += -2 * dqdci;         /* dwdot[OH]/d[H2O2] */
        J[252] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*dcdc_fac;
        J[423] += -2 * dqdci;         /* dwdot[OH]/d[CO] */
        J[427] += dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*dcdc_fac;
        J[633] += -2 * dqdci;         /* dwdot[OH]/d[CO2] */
        J[637] += dqdci;              /* dwdot[H2O2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac + k_f*2.000000*sc[3];
        dqdc[4] = 2.5*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = 12*dcdc_fac;
        dqdc[7] = dcdc_fac - k_r;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = 1.8999999999999999*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3.7999999999999998*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        for (int k=0; k<34; k++) {
            J[35*k+3] += -2 * dqdc[k];
            J[35*k+7] += dqdc[k];
        }
    }
    J[1193] += -2 * dqdT; /* dwdot[OH]/dT */
    J[1197] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 2: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = sc[2]*sc[5];
    k_f = 1.0000000000000002e-06 * 1475000000000
                * exp(0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0.59999999999999998 * invT + 0.50321666580471969 *  0  * invT2;
    /* pressure-fall-off */
    k_0 = 35000000000000000 * exp(-0.40999999999999998 * tc[0] - 0.50321666580471969 * (-1115.9200000000001) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -0.40999999999999998 * invT + 0.50321666580471969 * (-1115.9200000000001) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = 0.;
    Fcent2 = 0.5 * exp(-T/1e+30);
    Fcent3 = exp(-1e+100 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        0.
        -Fcent2/1e+30
        + Fcent3*1e+100*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[8];
    Kc = refCinv * exp(g_RT[2] + g_RT[5] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[8]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[72] -= dqdci;               /* dwdot[H]/d[H] */
        J[75] -= dqdci;               /* dwdot[O2]/d[H] */
        J[78] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*dcdc_fac;
        J[142] -= dqdci;              /* dwdot[H]/d[H2] */
        J[145] -= dqdci;              /* dwdot[O2]/d[H2] */
        J[148] += dqdci;              /* dwdot[HO2]/d[H2] */
        /* d()/d[O2] */
        dqdci =  + k_f*sc[2];
        J[177] -= dqdci;              /* dwdot[H]/d[O2] */
        J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
        J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*dcdc_fac;
        J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[215] -= dqdci;              /* dwdot[O2]/d[H2O] */
        J[218] += dqdci;              /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[282] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[285] -= dqdci;              /* dwdot[O2]/d[HO2] */
        J[288] += dqdci;              /* dwdot[HO2]/d[HO2] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*dcdc_fac;
        J[422] -= dqdci;              /* dwdot[H]/d[CO] */
        J[425] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[428] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*dcdc_fac;
        J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[635] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[638] += dqdci;              /* dwdot[HO2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[5];
        dqdc[3] = dcdc_fac;
        dqdc[4] = 2.5*dcdc_fac;
        dqdc[5] = dcdc_fac + k_f*sc[2];
        dqdc[6] = 12*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac - k_r;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = 1.8999999999999999*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3.7999999999999998*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        for (int k=0; k<34; k++) {
            J[35*k+2] -= dqdc[k];
            J[35*k+5] -= dqdc[k];
            J[35*k+8] += dqdc[k];
        }
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1195] -= dqdT; /* dwdot[O2]/dT */
    J[1198] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 3: CH3 + H (+M) => CH4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[4] + ( 5 - 1)*sc[6] + ( 3 - 1)*sc[18] + ( 2 - 1)*sc[12];
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = 1.0000000000000002e-06 * 2138000000000000
                * exp(-0.40000000000000002 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.40000000000000002 * invT + 0.50321666580471969 *  0  * invT2;
    /* pressure-fall-off */
    k_0 = 3.31e+30 * exp(-4 * tc[0] - 0.50321666580471969 * (2108.0300000000002) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -4 * invT + 0.50321666580471969 * (2108.0300000000002) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = 0.;
    Fcent2 = 0.;
    Fcent3 = exp(-40 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        0.
        0.
        + Fcent3*40*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[11];
        J[72] -= dqdci;               /* dwdot[H]/d[H] */
        J[81] -= dqdci;               /* dwdot[CH3]/d[H] */
        J[85] += dqdci;               /* dwdot[CH4]/d[H] */
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[142] -= dqdci;              /* dwdot[H]/d[H2] */
        J[151] -= dqdci;              /* dwdot[CH3]/d[H2] */
        J[155] += dqdci;              /* dwdot[CH4]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (5 - 1)*dcdc_fac;
        J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[221] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[225] += dqdci;              /* dwdot[CH4]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[2];
        J[387] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
        /* d()/d[CO] */
        dqdci = (2 - 1)*dcdc_fac;
        J[422] -= dqdci;              /* dwdot[H]/d[CO] */
        J[431] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[435] += dqdci;              /* dwdot[CH4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3 - 1)*dcdc_fac;
        J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[641] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[645] += dqdci;              /* dwdot[CH4]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[11];
        dqdc[3] = dcdc_fac;
        dqdc[4] = 2*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = 5*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac + k_f*sc[2];
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        for (int k=0; k<34; k++) {
            J[35*k+2] -= dqdc[k];
            J[35*k+11] -= dqdc[k];
            J[35*k+15] += dqdc[k];
        }
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1201] -= dqdT; /* dwdot[CH3]/dT */
    J[1205] += dqdT; /* dwdot[CH4]/dT */

    /*reaction 4: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[4] + ( 5 - 1)*sc[6] + ( 3 - 1)*sc[18] + ( 2 - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[11], 2.000000);
    k_f = 1.0000000000000002e-06 * 92140000000000000
                * exp(-1.1699999999999999 * tc[0] - 0.50321666580471969 * (635.75999999999999) * invT);
    dlnkfdT = -1.1699999999999999 * invT + 0.50321666580471969 *  635.75999999999999  * invT2;
    /* pressure-fall-off */
    k_0 = 1.135e+36 * exp(-5.2460000000000004 * tc[0] - 0.50321666580471969 * (1705.0699999999999) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -5.2460000000000004 * invT + 0.50321666580471969 * (1705.0699999999999) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.40500000000000003)*exp(-T/1120);
    Fcent2 = 0.40500000000000003 * exp(-T/69.599999999999994);
    Fcent3 = exp(-1000000000000000 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1120
        -Fcent2/69.599999999999994
        + Fcent3*1000000000000000*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[17];
    Kc = refCinv * exp(2.000000*g_RT[11] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[11]) + (h_RT[17]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[11] -= 2 * q; /* CH3 */
    wdot[17] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[151] += -2 * dqdci;         /* dwdot[CH3]/d[H2] */
        J[157] += dqdci;              /* dwdot[C2H6]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (5 - 1)*dcdc_fac;
        J[221] += -2 * dqdci;         /* dwdot[CH3]/d[H2O] */
        J[227] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*2.000000*sc[11];
        J[396] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
        J[402] += dqdci;              /* dwdot[C2H6]/d[CH3] */
        /* d()/d[CO] */
        dqdci = (2 - 1)*dcdc_fac;
        J[431] += -2 * dqdci;         /* dwdot[CH3]/d[CO] */
        J[437] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[C2H6] */
        dqdci =  - k_r;
        J[606] += -2 * dqdci;         /* dwdot[CH3]/d[C2H6] */
        J[612] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[CO2] */
        dqdci = (3 - 1)*dcdc_fac;
        J[641] += -2 * dqdci;         /* dwdot[CH3]/d[CO2] */
        J[647] += dqdci;              /* dwdot[C2H6]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = 2*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = 5*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac + k_f*2.000000*sc[11];
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac - k_r;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        for (int k=0; k<34; k++) {
            J[35*k+11] += -2 * dqdc[k];
            J[35*k+17] += dqdc[k];
        }
    }
    J[1201] += -2 * dqdT; /* dwdot[CH3]/dT */
    J[1207] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 5: CO + O (+M) => CO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = sc[1]*sc[12];
    k_f = 1.0000000000000002e-06 * 18000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (2384.0799999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  2384.0799999999999  * invT2;
    /* pressure-fall-off */
    k_0 = 1.35e+24 * exp(-2.7879999999999998 * tc[0] - 0.50321666580471969 * (4190.9700000000003) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.7879999999999998 * invT + 0.50321666580471969 * (4190.9700000000003) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-1)*exp(-T/1);
    Fcent2 = 1 * exp(-T/10000000);
    Fcent3 = exp(-10000000 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1
        -Fcent2/10000000
        + Fcent3*10000000*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[12] -= q; /* CO */
    wdot[18] += q; /* CO2 */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[O] */
        dqdci =  + k_f*sc[12];
        J[36] -= dqdci;               /* dwdot[O]/d[O] */
        J[47] -= dqdci;               /* dwdot[CO]/d[O] */
        J[53] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*dcdc_fac;
        J[141] -= dqdci;              /* dwdot[O]/d[H2] */
        J[152] -= dqdci;              /* dwdot[CO]/d[H2] */
        J[158] += dqdci;              /* dwdot[CO2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*dcdc_fac;
        J[211] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[222] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[228] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*dcdc_fac + k_f*sc[1];
        J[421] -= dqdci;              /* dwdot[O]/d[CO] */
        J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[438] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*dcdc_fac;
        J[631] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[642] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[648] += dqdci;              /* dwdot[CO2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[12];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = 2.5*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = 12*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = 1.8999999999999999*dcdc_fac + k_f*sc[1];
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3.7999999999999998*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        for (int k=0; k<34; k++) {
            J[35*k+1] -= dqdc[k];
            J[35*k+12] -= dqdc[k];
            J[35*k+18] += dqdc[k];
        }
    }
    J[1191] -= dqdT; /* dwdot[O]/dT */
    J[1202] -= dqdT; /* dwdot[CO]/dT */
    J[1208] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 6: CH3O (+M) => CH2O + H (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[13];
    k_f = 1 * 54500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (13500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  13500  * invT2;
    /* pressure-fall-off */
    k_0 = 2.344e+25 * exp(-2.7000000000000002 * tc[0] - 0.50321666580471969 * (30599.900000000001) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.7000000000000002 * invT + 0.50321666580471969 * (30599.900000000001) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-1)*exp(-T/1);
    Fcent2 = 1 * exp(-T/10000000);
    Fcent3 = exp(-10000000 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1
        -Fcent2/10000000
        + Fcent3*10000000*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[10] += q; /* CH2O */
    wdot[13] -= q; /* CH3O */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[CH3O] */
        dqdci =  + k_f;
        J[457] += dqdci;              /* dwdot[H]/d[CH3O] */
        J[465] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
        J[468] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac + k_f;
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
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        for (int k=0; k<34; k++) {
            J[35*k+2] += dqdc[k];
            J[35*k+10] += dqdc[k];
            J[35*k+13] -= dqdc[k];
        }
    }
    J[1192] += dqdT; /* dwdot[H]/dT */
    J[1200] += dqdT; /* dwdot[CH2O]/dT */
    J[1203] -= dqdT; /* dwdot[CH3O]/dT */

    /*reaction 7: C2H3 (+M) => H + C2H2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[4] + ( 5 - 1)*sc[6] + ( 3 - 1)*sc[18] + ( 2 - 1)*sc[12];
    /* forward */
    phi_f = sc[24];
    k_f = 1 * 16060000000
                * exp(1.028 * tc[0] - 0.50321666580471969 * (40503.589999999997) * invT);
    dlnkfdT = 1.028 * invT + 0.50321666580471969 *  40503.589999999997  * invT2;
    /* pressure-fall-off */
    k_0 = 1.164e+39 * exp(-6.8209999999999997 * tc[0] - 0.50321666580471969 * (44491.629999999997) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -6.8209999999999997 * invT + 0.50321666580471969 * (44491.629999999997) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = 0.;
    Fcent2 = 1 * exp(-T/675);
    Fcent3 = exp(-1000000000000000 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        0.
        -Fcent2/675
        + Fcent3*1000000000000000*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[22] += q; /* C2H2 */
    wdot[24] -= q; /* C2H3 */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[142] += dqdci;              /* dwdot[H]/d[H2] */
        J[162] += dqdci;              /* dwdot[C2H2]/d[H2] */
        J[164] -= dqdci;              /* dwdot[C2H3]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (5 - 1)*dcdc_fac;
        J[212] += dqdci;              /* dwdot[H]/d[H2O] */
        J[232] += dqdci;              /* dwdot[C2H2]/d[H2O] */
        J[234] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (2 - 1)*dcdc_fac;
        J[422] += dqdci;              /* dwdot[H]/d[CO] */
        J[442] += dqdci;              /* dwdot[C2H2]/d[CO] */
        J[444] -= dqdci;              /* dwdot[C2H3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3 - 1)*dcdc_fac;
        J[632] += dqdci;              /* dwdot[H]/d[CO2] */
        J[652] += dqdci;              /* dwdot[C2H2]/d[CO2] */
        J[654] -= dqdci;              /* dwdot[C2H3]/d[CO2] */
        /* d()/d[C2H3] */
        dqdci =  + k_f;
        J[842] += dqdci;              /* dwdot[H]/d[C2H3] */
        J[862] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
        J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = 2*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = 5*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac + k_f;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        for (int k=0; k<34; k++) {
            J[35*k+2] += dqdc[k];
            J[35*k+22] += dqdc[k];
            J[35*k+24] -= dqdc[k];
        }
    }
    J[1192] += dqdT; /* dwdot[H]/dT */
    J[1212] += dqdT; /* dwdot[C2H2]/dT */
    J[1214] -= dqdT; /* dwdot[C2H3]/dT */

    /*reaction 8: H + C2H4 (+M) <=> C2H5 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[2]*sc[16];
    k_f = 1.0000000000000002e-06 * 1081000000000
                * exp(0.45000000000000001 * tc[0] - 0.50321666580471969 * (1821.9400000000001) * invT);
    dlnkfdT = 0.45000000000000001 * invT + 0.50321666580471969 *  1821.9400000000001  * invT2;
    /* pressure-fall-off */
    k_0 = 1.1120000000000001e+34 * exp(-5 * tc[0] - 0.50321666580471969 * (4447.8999999999996) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -5 * invT + 0.50321666580471969 * (4447.8999999999996) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = 0.;
    Fcent2 = 1 * exp(-T/95);
    Fcent3 = exp(-200 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        0.
        -Fcent2/95
        + Fcent3*200*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[14];
    Kc = refCinv * exp(g_RT[2] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[16]) + (h_RT[14]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[14] += q; /* C2H5 */
    wdot[16] -= q; /* C2H4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[16];
        J[72] -= dqdci;               /* dwdot[H]/d[H] */
        J[84] += dqdci;               /* dwdot[C2H5]/d[H] */
        J[86] -= dqdci;               /* dwdot[C2H4]/d[H] */
        /* d()/d[C2H5] */
        dqdci =  - k_r;
        J[492] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[504] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
        J[506] -= dqdci;              /* dwdot[C2H4]/d[C2H5] */
        /* d()/d[C2H4] */
        dqdci =  + k_f*sc[2];
        J[562] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[574] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
        J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[16];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac - k_r;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac + k_f*sc[2];
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        for (int k=0; k<34; k++) {
            J[35*k+2] -= dqdc[k];
            J[35*k+14] += dqdc[k];
            J[35*k+16] -= dqdc[k];
        }
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1204] += dqdT; /* dwdot[C2H5]/dT */
    J[1206] -= dqdT; /* dwdot[C2H4]/dT */

    /*reaction 9: C2H4 (+M) => C2H2 + H2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[16];
    k_f = 1 * 18000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (76000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  76000  * invT2;
    /* pressure-fall-off */
    k_0 = 1500000000000000 * exp(0 * tc[0] - 0.50321666580471969 * (55440.010000000002) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = 0 * invT + 0.50321666580471969 * (55440.010000000002) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-1)*exp(-T/1);
    Fcent2 = 1 * exp(-T/10000000);
    Fcent3 = exp(-10000000 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1
        -Fcent2/10000000
        + Fcent3*10000000*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[4] += q; /* H2 */
    wdot[16] -= q; /* C2H4 */
    wdot[22] += q; /* C2H2 */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[C2H4] */
        dqdci =  + k_f;
        J[564] += dqdci;              /* dwdot[H2]/d[C2H4] */
        J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
        J[582] += dqdci;              /* dwdot[C2H2]/d[C2H4] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac + k_f;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        for (int k=0; k<34; k++) {
            J[35*k+4] += dqdc[k];
            J[35*k+16] -= dqdc[k];
            J[35*k+22] += dqdc[k];
        }
    }
    J[1194] += dqdT; /* dwdot[H2]/dT */
    J[1206] -= dqdT; /* dwdot[C2H4]/dT */
    J[1212] += dqdT; /* dwdot[C2H2]/dT */

    /*reaction 10: O + H + M => OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-12 * 4.72e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[O] */
        dqdci =  + k_f*sc[2];
        J[36] -= dqdci;               /* dwdot[O]/d[O] */
        J[37] -= dqdci;               /* dwdot[H]/d[O] */
        J[38] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[1];
        J[71] -= dqdci;               /* dwdot[O]/d[H] */
        J[72] -= dqdci;               /* dwdot[H]/d[H] */
        J[73] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor;
        J[141] -= dqdci;              /* dwdot[O]/d[H2] */
        J[142] -= dqdci;              /* dwdot[H]/d[H2] */
        J[143] += dqdci;              /* dwdot[OH]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor;
        J[211] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor;
        J[421] -= dqdci;              /* dwdot[O]/d[CO] */
        J[422] -= dqdci;              /* dwdot[H]/d[CO] */
        J[423] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[631] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[633] += dqdci;              /* dwdot[OH]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*sc[2];
        dqdc[2] = q_nocor + k_f*sc[1];
        dqdc[3] = q_nocor;
        dqdc[4] = 2.5*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = 12*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = 1.8999999999999999*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3.7999999999999998*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        for (int k=0; k<34; k++) {
            J[35*k+1] -= dqdc[k];
            J[35*k+2] -= dqdc[k];
            J[35*k+3] += dqdc[k];
        }
    }
    J[1191] -= dqdT; /* dwdot[O]/dT */
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1193] += dqdT; /* dwdot[OH]/dT */

    /*reaction 11: 2.000000 O + M => O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 6170000000000000
                * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.5 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[1];
        J[36] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[40] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor;
        J[141] += -2 * dqdci;         /* dwdot[O]/d[H2] */
        J[145] += dqdci;              /* dwdot[O2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor;
        J[211] += -2 * dqdci;         /* dwdot[O]/d[H2O] */
        J[215] += dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor;
        J[421] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[425] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[631] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[635] += dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*2.000000*sc[1];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = 2.5*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = 12*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = 1.8999999999999999*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3.7999999999999998*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        for (int k=0; k<34; k++) {
            J[35*k+1] += -2 * dqdc[k];
            J[35*k+5] += dqdc[k];
        }
    }
    J[1191] += -2 * dqdT; /* dwdot[O]/dT */
    J[1195] += dqdT; /* dwdot[O2]/dT */

    /*reaction 12: H + OH + M => H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 12 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = sc[2]*sc[3];
    k_f = 1.0000000000000002e-12 * 2.2499999999999999e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[3];
        J[72] -= dqdci;               /* dwdot[H]/d[H] */
        J[73] -= dqdci;               /* dwdot[OH]/d[H] */
        J[76] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[2];
        J[107] -= dqdci;              /* dwdot[H]/d[OH] */
        J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
        J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor;
        J[142] -= dqdci;              /* dwdot[H]/d[H2] */
        J[143] -= dqdci;              /* dwdot[OH]/d[H2] */
        J[146] += dqdci;              /* dwdot[H2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor;
        J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[213] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[216] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor;
        J[422] -= dqdci;              /* dwdot[H]/d[CO] */
        J[423] -= dqdci;              /* dwdot[OH]/d[CO] */
        J[426] += dqdci;              /* dwdot[H2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[633] -= dqdci;              /* dwdot[OH]/d[CO2] */
        J[636] += dqdci;              /* dwdot[H2O]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*sc[3];
        dqdc[3] = q_nocor + k_f*sc[2];
        dqdc[4] = 2.5*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = 12*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = 1.8999999999999999*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3.7999999999999998*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        for (int k=0; k<34; k++) {
            J[35*k+2] -= dqdc[k];
            J[35*k+3] -= dqdc[k];
            J[35*k+6] += dqdc[k];
        }
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1193] -= dqdT; /* dwdot[OH]/dT */
    J[1196] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 13: HCO + M => H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 6 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = sc[19];
    k_f = 1.0000000000000002e-06 * 1.86e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  17000  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor;
        J[142] += dqdci;              /* dwdot[H]/d[H2] */
        J[152] += dqdci;              /* dwdot[CO]/d[H2] */
        J[159] -= dqdci;              /* dwdot[HCO]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*q_nocor;
        J[212] += dqdci;              /* dwdot[H]/d[H2O] */
        J[222] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[229] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor;
        J[422] += dqdci;              /* dwdot[H]/d[CO] */
        J[432] += dqdci;              /* dwdot[CO]/d[CO] */
        J[439] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[632] += dqdci;              /* dwdot[H]/d[CO2] */
        J[642] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[649] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[667] += dqdci;              /* dwdot[H]/d[HCO] */
        J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = 2.5*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = 6*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = 1.8999999999999999*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3.7999999999999998*q_nocor;
        dqdc[19] = q_nocor + k_f;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        for (int k=0; k<34; k++) {
            J[35*k+2] += dqdc[k];
            J[35*k+12] += dqdc[k];
            J[35*k+19] -= dqdc[k];
        }
    }
    J[1192] += dqdT; /* dwdot[H]/dT */
    J[1202] += dqdT; /* dwdot[CO]/dT */
    J[1209] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 14: H + CO + M => HCO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[4] + ( 6 - 1)*sc[6] + ( 3.7999999999999998 - 1)*sc[18] + ( 1.8999999999999999 - 1)*sc[12];
    /* forward */
    phi_f = sc[2]*sc[12];
    k_f = 1.0000000000000002e-12 * 64670000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-441.92000000000002) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -441.92000000000002  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[12] -= q; /* CO */
    wdot[19] += q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[12];
        J[72] -= dqdci;               /* dwdot[H]/d[H] */
        J[82] -= dqdci;               /* dwdot[CO]/d[H] */
        J[89] += dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor;
        J[142] -= dqdci;              /* dwdot[H]/d[H2] */
        J[152] -= dqdci;              /* dwdot[CO]/d[H2] */
        J[159] += dqdci;              /* dwdot[HCO]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*q_nocor;
        J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[222] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[229] += dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor + k_f*sc[2];
        J[422] -= dqdci;              /* dwdot[H]/d[CO] */
        J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[439] += dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[642] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[649] += dqdci;              /* dwdot[HCO]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*sc[12];
        dqdc[3] = q_nocor;
        dqdc[4] = 2.5*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = 6*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = 1.8999999999999999*q_nocor + k_f*sc[2];
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3.7999999999999998*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        for (int k=0; k<34; k++) {
            J[35*k+2] -= dqdc[k];
            J[35*k+12] -= dqdc[k];
            J[35*k+19] += dqdc[k];
        }
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1202] -= dqdT; /* dwdot[CO]/dT */
    J[1209] += dqdT; /* dwdot[HCO]/dT */

    /*reaction 15: CH3O2 + M => CH3 + O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[20];
    k_f = 1.0000000000000002e-06 * 4.3430000000000002e+27
                * exp(-3.4199999999999999 * tc[0] - 0.50321666580471969 * (30469.889999999999) * invT);
    dlnkfdT = -3.4199999999999999 * invT + 0.50321666580471969 *  30469.889999999999  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[11] += q; /* CH3 */
    wdot[20] -= q; /* CH3O2 */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[CH3O2] */
        dqdci =  + k_f;
        J[705] += dqdci;              /* dwdot[O2]/d[CH3O2] */
        J[711] += dqdci;              /* dwdot[CH3]/d[CH3O2] */
        J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
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
        dqdc[20] = q_nocor + k_f;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        for (int k=0; k<34; k++) {
            J[35*k+5] += dqdc[k];
            J[35*k+11] += dqdc[k];
            J[35*k+20] -= dqdc[k];
        }
    }
    J[1195] += dqdT; /* dwdot[O2]/dT */
    J[1201] += dqdT; /* dwdot[CH3]/dT */
    J[1210] -= dqdT; /* dwdot[CH3O2]/dT */

    /*reaction 16: CH3 + O2 + M => CH3O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: irreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = 1.0000000000000002e-12 * 5.4400000000000001e+25
                * exp(-3.2999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -3.2999999999999998 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[11] -= q; /* CH3 */
    wdot[20] += q; /* CH3O2 */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[O2] */
        dqdci =  + k_f*sc[11];
        J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
        J[186] -= dqdci;              /* dwdot[CH3]/d[O2] */
        J[195] += dqdci;              /* dwdot[CH3O2]/d[O2] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[5];
        J[390] -= dqdci;              /* dwdot[O2]/d[CH3] */
        J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[405] += dqdci;              /* dwdot[CH3O2]/d[CH3] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor + k_f*sc[11];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor + k_f*sc[5];
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
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        for (int k=0; k<34; k++) {
            J[35*k+5] -= dqdc[k];
            J[35*k+11] -= dqdc[k];
            J[35*k+20] += dqdc[k];
        }
    }
    J[1195] -= dqdT; /* dwdot[O2]/dT */
    J[1201] -= dqdT; /* dwdot[CH3]/dT */
    J[1210] += dqdT; /* dwdot[CH3O2]/dT */

    /*reaction 17: O + H2 => H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = 1.0000000000000002e-06 * 50800
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6292.0699999999997) * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  6292.0699999999997  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[4];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[39] -= dqdci;               /* dwdot[H2]/d[O] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[1];
    J[141] -= dqdci;              /* dwdot[O]/d[H2] */
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[143] += dqdci;              /* dwdot[OH]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */

    /*reaction 18: H + OH => O + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[3];
    k_f = 1.0000000000000002e-06 * 22310
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (4196.9399999999996) * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  4196.9399999999996  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[2] -= q; /* H */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3];
    J[71] += dqdci;               /* dwdot[O]/d[H] */
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] -= dqdci;               /* dwdot[OH]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[2];
    J[106] += dqdci;              /* dwdot[O]/d[OH] */
    J[107] -= dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[109] += dqdci;              /* dwdot[H2]/d[OH] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */

    /*reaction 19: OH + H2 => H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[4];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3429.9699999999998) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  3429.9699999999998  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[4] -= q; /* H2 */
    wdot[6] += q; /* H2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[4];
    J[107] += dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[109] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[143] -= dqdci;              /* dwdot[OH]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[146] += dqdci;              /* dwdot[H2O]/d[H2] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */

    /*reaction 20: H + H2O => OH + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = 1.0000000000000002e-06 * 935200000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (18580.07) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  18580.07  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] += q; /* H2 */
    wdot[6] -= q; /* H2O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[OH]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[76] -= dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[2];
    J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[214] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */

    /*reaction 21: O + H2O => 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 2970000
                * exp(2.02 * tc[0] - 0.50321666580471969 * (13400.1) * invT);
    dlnkfdT = 2.02 * invT + 0.50321666580471969 *  13400.1  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += 2 * q; /* OH */
    wdot[6] -= q; /* H2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += 2 * dqdci;           /* dwdot[OH]/d[O] */
    J[41] -= dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1];
    J[211] -= dqdci;              /* dwdot[O]/d[H2O] */
    J[213] += 2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += 2 * dqdT;          /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */

    /*reaction 22: 2.000000 OH => O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[3], 2.000000);
    k_f = 1.0000000000000002e-06 * 301300
                * exp(2.02 * tc[0] - 0.50321666580471969 * (-3849.9000000000001) * invT);
    dlnkfdT = 2.02 * invT + 0.50321666580471969 *  -3849.9000000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[3] -= 2 * q; /* OH */
    wdot[6] += q; /* H2O */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[3];
    J[106] += dqdci;              /* dwdot[O]/d[OH] */
    J[108] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1193] += -2 * dqdT;         /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */

    /*reaction 23: H + O2 => O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[5];
    k_f = 1.0000000000000002e-06 * 197000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (16539.91) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  16539.91  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[5];
    J[71] += dqdci;               /* dwdot[O]/d[H] */
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[OH]/d[H] */
    J[75] -= dqdci;               /* dwdot[O2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[2];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[177] -= dqdci;              /* dwdot[H]/d[O2] */
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */

    /*reaction 24: O + OH => H + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 15550000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (424.94999999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  424.94999999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[3];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[38] -= dqdci;               /* dwdot[OH]/d[O] */
    J[40] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[106] -= dqdci;              /* dwdot[O]/d[OH] */
    J[107] += dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[110] += dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */

    /*reaction 25: HO2 + H => 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = 1.0000000000000002e-06 * 70800000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (299.94999999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  299.94999999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += 2 * q; /* OH */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[78] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[282] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[283] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += 2 * dqdT;          /* dwdot[OH]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 26: 2.000000 HO2 => H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[8], 2.000000);
    k_f = 1.0000000000000002e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (11979.92) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  11979.92  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= 2 * q; /* HO2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[8];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += -2 * dqdT;         /* dwdot[HO2]/dT */

    /*reaction 27: HO2 + H => H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = 1.0000000000000002e-06 * 16600000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (820.02999999999997) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  820.02999999999997  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[75] += dqdci;               /* dwdot[O2]/d[H] */
    J[78] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[282] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[284] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 28: HO2 + OH => H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = 1.0000000000000002e-06 * 28900000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -500  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[6] += q; /* H2O */
    wdot[8] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[110] += dqdci;              /* dwdot[O2]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[113] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[283] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[286] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 29: H2O + O2 => HO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = 1.0000000000000002e-06 * 6888000000000000
                * exp(-0.33000000000000002 * tc[0] - 0.50321666580471969 * (72140.059999999998) * invT);
    dlnkfdT = -0.33000000000000002 * invT + 0.50321666580471969 *  72140.059999999998  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[6] -= q; /* H2O */
    wdot[8] += q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[6];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[181] -= dqdci;              /* dwdot[H2O]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[5];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[215] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[218] += dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 30: 2.000000 HO2 => H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[8], 2.000000);
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1629.0599999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1629.0599999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= 2 * q; /* HO2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[8];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += -2 * dqdT;         /* dwdot[HO2]/dT */

    /*reaction 31: HO2 + O => OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[8];
    k_f = 1.0000000000000002e-06 * 32500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[8];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[40] += dqdci;               /* dwdot[O2]/d[O] */
    J[43] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[281] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 32: OH + O2 => HO2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[5];
    k_f = 1.0000000000000002e-06 * 785700000000000
                * exp(-0.33000000000000002 * tc[0] - 0.50321666580471969 * (55390.059999999998) * invT);
    dlnkfdT = -0.33000000000000002 * invT + 0.50321666580471969 *  55390.059999999998  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[3] -= q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[5];
    J[106] += dqdci;              /* dwdot[O]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[110] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[113] += dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[3];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[178] -= dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 33: H2O2 + OH => H2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 1000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[112] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    J[113] += dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[248] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[251] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 34: H2O2 + H => H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 24100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3969.8899999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3969.8899999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[OH]/d[H] */
    J[76] += dqdci;               /* dwdot[H2O]/d[H] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[247] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[248] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[251] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 35: H2O2 + OH => H2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 580000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (9559.9899999999998) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  9559.9899999999998  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[112] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    J[113] += dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[248] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[251] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 36: H2O + HO2 => H2O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[8];
    k_f = 1.0000000000000002e-06 * 97710000000000
                * exp(0.33000000000000002 * tc[0] - 0.50321666580471969 * (41020.080000000002) * invT);
    dlnkfdT = 0.33000000000000002 * invT + 0.50321666580471969 *  41020.080000000002  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* H2O */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[8];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O2]/d[H2O] */
    J[218] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[6];
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[286] -= dqdci;              /* dwdot[H2O]/d[HO2] */
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 37: H2O2 + O => OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = 1.0000000000000002e-06 * 9550000
                * exp(2 * tc[0] - 0.50321666580471969 * (3969.8899999999999) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  3969.8899999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[42] -= dqdci;               /* dwdot[H2O2]/d[O] */
    J[43] += dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[246] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[248] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 38: H2O2 + H => H2 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 48200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7950.0500000000002) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  7950.0500000000002  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H] */
    J[78] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[247] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[249] += dqdci;              /* dwdot[H2]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 39: H2 + HO2 => H2O2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[8];
    k_f = 1.0000000000000002e-06 * 1875000000000
                * exp(0.33000000000000002 * tc[0] - 0.50321666580471969 * (24260.040000000001) * invT);
    dlnkfdT = 0.33000000000000002 * invT + 0.50321666580471969 *  24260.040000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[8];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[147] += dqdci;              /* dwdot[H2O2]/d[H2] */
    J[148] -= dqdci;              /* dwdot[HO2]/d[H2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[282] += dqdci;              /* dwdot[H]/d[HO2] */
    J[284] -= dqdci;              /* dwdot[H2]/d[HO2] */
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 40: CH2GSG + OH => CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[9] -= q; /* CH2GSG */
    wdot[10] += q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[107] += dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[114] -= dqdci;              /* dwdot[CH2GSG]/d[OH] */
    J[115] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[3];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[318] -= dqdci;              /* dwdot[OH]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[325] += dqdci;              /* dwdot[CH2O]/d[CH2GSG] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 41: CH2GSG + H2 => CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = 1.0000000000000002e-06 * 70000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[9];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[149] -= dqdci;              /* dwdot[CH2GSG]/d[H2] */
    J[151] += dqdci;              /* dwdot[CH3]/d[H2] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[4];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[319] -= dqdci;              /* dwdot[H2]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] += dqdci;              /* dwdot[CH3]/d[CH2GSG] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */

    /*reaction 42: CH3 + H => CH2GSG + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = 1.0000000000000002e-06 * 2.482e+17
                * exp(-0.89000000000000001 * tc[0] - 0.50321666580471969 * (16130.02) * invT);
    dlnkfdT = -0.89000000000000001 * invT + 0.50321666580471969 *  16130.02  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[9] += q; /* CH2GSG */
    wdot[11] -= q; /* CH3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[11];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[79] += dqdci;               /* dwdot[CH2GSG]/d[H] */
    J[81] -= dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[2];
    J[387] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[389] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[394] += dqdci;              /* dwdot[CH2GSG]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1199] += dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 43: CH2GSG + O2 => CO + OH + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[9];
    k_f = 1.0000000000000002e-06 * 70000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[9] -= q; /* CH2GSG */
    wdot[12] += q; /* CO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[177] += dqdci;              /* dwdot[H]/d[O2] */
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[184] -= dqdci;              /* dwdot[CH2GSG]/d[O2] */
    J[187] += dqdci;              /* dwdot[CO]/d[O2] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[5];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[318] += dqdci;              /* dwdot[OH]/d[CH2GSG] */
    J[320] -= dqdci;              /* dwdot[O2]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[327] += dqdci;              /* dwdot[CO]/d[CH2GSG] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 44: CH2GSG + O => CO + 2.000000 H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += 2 * q; /* H */
    wdot[9] -= q; /* CH2GSG */
    wdot[12] += q; /* CO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[9];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += 2 * dqdci;           /* dwdot[H]/d[O] */
    J[44] -= dqdci;               /* dwdot[CH2GSG]/d[O] */
    J[47] += dqdci;               /* dwdot[CO]/d[O] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[1];
    J[316] -= dqdci;              /* dwdot[O]/d[CH2GSG] */
    J[317] += 2 * dqdci;          /* dwdot[H]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[327] += dqdci;              /* dwdot[CO]/d[CH2GSG] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += 2 * dqdT;          /* dwdot[H]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 45: CH3 + HO2 => CH3O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = 1.0000000000000002e-06 * 11000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[11] -= q; /* CH3 */
    wdot[13] += q; /* CH3O */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[291] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[293] += dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[388] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[393] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[398] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1203] += dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 46: CH3 + O2 => CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = 1.0000000000000002e-06 * 747000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (14250) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  14250  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[11] -= q; /* CH3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[186] -= dqdci;              /* dwdot[CH3]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[5];
    J[388] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[390] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[395] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 47: CH3 + O2 => CH3O + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = 1.0000000000000002e-06 * 1.995e+18
                * exp(-1.5700000000000001 * tc[0] - 0.50321666580471969 * (29210.09) * invT);
    dlnkfdT = -1.5700000000000001 * invT + 0.50321666580471969 *  29210.09  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[11] -= q; /* CH3 */
    wdot[13] += q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[186] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[188] += dqdci;              /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[5];
    J[386] += dqdci;              /* dwdot[O]/d[CH3] */
    J[390] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[398] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1203] += dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 48: CH3O + O => CH3 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = 1.0000000000000002e-06 * 3.585e+18
                * exp(-1.5900000000000001 * tc[0] - 0.50321666580471969 * (-1630.98) * invT);
    dlnkfdT = -1.5900000000000001 * invT + 0.50321666580471969 *  -1630.98  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[5] += q; /* O2 */
    wdot[11] += q; /* CH3 */
    wdot[13] -= q; /* CH3O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[40] += dqdci;               /* dwdot[O2]/d[O] */
    J[46] += dqdci;               /* dwdot[CH3]/d[O] */
    J[48] -= dqdci;               /* dwdot[CH3O]/d[O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[1];
    J[456] -= dqdci;              /* dwdot[O]/d[CH3O] */
    J[460] += dqdci;              /* dwdot[O2]/d[CH3O] */
    J[466] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[468] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1203] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 49: 2.000000 CH3 <=> H + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[11], 2.000000);
    k_f = 1.0000000000000002e-06 * 6840000000000
                * exp(0.10000000000000001 * tc[0] - 0.50321666580471969 * (10599.9) * invT);
    dlnkfdT = 0.10000000000000001 * invT + 0.50321666580471969 *  10599.9  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[14];
    Kc = exp(-g_RT[2] + 2.000000*g_RT[11] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[11]) + (h_RT[2] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[11] -= 2 * q; /* CH3 */
    wdot[14] += q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[72] += dqdci;               /* dwdot[H]/d[H] */
    J[81] += -2 * dqdci;          /* dwdot[CH3]/d[H] */
    J[84] += dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*2.000000*sc[11];
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[396] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[399] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[2];
    J[492] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[501] += -2 * dqdci;         /* dwdot[CH3]/d[C2H5] */
    J[504] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1201] += -2 * dqdT;         /* dwdot[CH3]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */

    /*reaction 50: CH3 + HO2 => CH4 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = 1.0000000000000002e-06 * 3600000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[291] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[295] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[390] += dqdci;              /* dwdot[O2]/d[CH3] */
    J[393] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 51: CH3 + O => CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = 1.0000000000000002e-06 * 80000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[10] += q; /* CH2O */
    wdot[11] -= q; /* CH3 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[11];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[45] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[46] -= dqdci;               /* dwdot[CH3]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[1];
    J[386] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[395] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 52: CH3 + OH => CH2O + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = 1.0000000000000002e-06 * 22500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (4299.9499999999998) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  4299.9499999999998  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2 */
    wdot[10] += q; /* CH2O */
    wdot[11] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[109] += dqdci;              /* dwdot[H2]/d[OH] */
    J[115] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[116] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[388] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[389] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[395] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 53: CH3 + OH => CH2GSG + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = 1.0000000000000002e-06 * 26500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (2185.9499999999998) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  2185.9499999999998  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[9] += q; /* CH2GSG */
    wdot[11] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[114] += dqdci;              /* dwdot[CH2GSG]/d[OH] */
    J[116] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[388] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[391] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[394] += dqdci;              /* dwdot[CH2GSG]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1199] += dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 54: CH2GSG + H2O => CH3 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[9];
    k_f = 1.0000000000000002e-06 * 32360000000
                * exp(0.89000000000000001 * tc[0] - 0.50321666580471969 * (1211.04) * invT);
    dlnkfdT = 0.89000000000000001 * invT + 0.50321666580471969 *  1211.04  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* H2O */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] += q; /* CH3 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[9];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[219] -= dqdci;              /* dwdot[CH2GSG]/d[H2O] */
    J[221] += dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[6];
    J[318] += dqdci;              /* dwdot[OH]/d[CH2GSG] */
    J[321] -= dqdci;              /* dwdot[H2O]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] += dqdci;              /* dwdot[CH3]/d[CH2GSG] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */

    /*reaction 55: CH3 + H2O2 => CH4 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[11];
    k_f = 1.0000000000000002e-06 * 336500000000
                * exp(-0.33000000000000002 * tc[0] - 0.50321666580471969 * (2501.9099999999999) * invT);
    dlnkfdT = -0.33000000000000002 * invT + 0.50321666580471969 *  2501.9099999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[11];
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[256] -= dqdci;              /* dwdot[CH3]/d[H2O2] */
    J[260] += dqdci;              /* dwdot[CH4]/d[H2O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[392] -= dqdci;              /* dwdot[H2O2]/d[CH3] */
    J[393] += dqdci;              /* dwdot[HO2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 56: CH2GSG + CH3 => C2H4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[11];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] -= q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[11];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] -= dqdci;              /* dwdot[CH3]/d[CH2GSG] */
    J[331] += dqdci;              /* dwdot[C2H4]/d[CH2GSG] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[9];
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[394] -= dqdci;              /* dwdot[CH2GSG]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[401] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 57: CH2GSG + CH4 => 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[15];
    k_f = 1.0000000000000002e-06 * 40000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] += 2 * q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[15];
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] += 2 * dqdci;          /* dwdot[CH3]/d[CH2GSG] */
    J[330] -= dqdci;              /* dwdot[CH4]/d[CH2GSG] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[9];
    J[534] -= dqdci;              /* dwdot[CH2GSG]/d[CH4] */
    J[536] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += 2 * dqdT;          /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 58: 2.000000 CH3 => CH2GSG + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[11], 2.000000);
    k_f = 1.0000000000000002e-06 * 5429000000000000
                * exp(-0.89000000000000001 * tc[0] - 0.50321666580471969 * (15650.1) * invT);
    dlnkfdT = -0.89000000000000001 * invT + 0.50321666580471969 *  15650.1  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] += q; /* CH2GSG */
    wdot[11] -= 2 * q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[CH3] */
    dqdci =  + k_f*2.000000*sc[11];
    J[394] += dqdci;              /* dwdot[CH2GSG]/d[CH3] */
    J[396] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1199] += dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += -2 * dqdT;         /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 59: CH4 + O => CH3 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[15];
    k_f = 1.0000000000000002e-06 * 3150000000000
                * exp(0.5 * tc[0] - 0.50321666580471969 * (10289.91) * invT);
    dlnkfdT = 0.5 * invT + 0.50321666580471969 *  10289.91  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[11] += q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[15];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[46] += dqdci;               /* dwdot[CH3]/d[O] */
    J[50] -= dqdci;               /* dwdot[CH4]/d[O] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[1];
    J[526] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[528] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[536] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 60: CH3 + OH => CH4 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = 1.0000000000000002e-06 * 52960000000
                * exp(0.5 * tc[0] - 0.50321666580471969 * (7715.1099999999997) * invT);
    dlnkfdT = 0.5 * invT + 0.50321666580471969 *  7715.1099999999997  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[3] -= q; /* OH */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[106] += dqdci;              /* dwdot[O]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[120] += dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[386] += dqdci;              /* dwdot[O]/d[CH3] */
    J[388] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 61: CH4 + H => CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[15];
    k_f = 1.0000000000000002e-06 * 17270
                * exp(3 * tc[0] - 0.50321666580471969 * (8223.9500000000007) * invT);
    dlnkfdT = 3 * invT + 0.50321666580471969 *  8223.9500000000007  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[11] += q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[15];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[81] += dqdci;               /* dwdot[CH3]/d[H] */
    J[85] -= dqdci;               /* dwdot[CH4]/d[H] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[527] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[529] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[536] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 62: CH3 + H2 => CH4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = 1.0000000000000002e-06 * 661
                * exp(3 * tc[0] - 0.50321666580471969 * (7744.0200000000004) * invT);
    dlnkfdT = 3 * invT + 0.50321666580471969 *  7744.0200000000004  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[11];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[151] -= dqdci;              /* dwdot[CH3]/d[H2] */
    J[155] += dqdci;              /* dwdot[CH4]/d[H2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[389] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 63: CH4 + OH => CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = 1.0000000000000002e-06 * 193000
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (2106.1199999999999) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  2106.1199999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[11] += q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[15];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[116] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[120] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[3];
    J[528] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[531] += dqdci;              /* dwdot[H2O]/d[CH4] */
    J[536] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 64: CH3 + H2O => CH4 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[11];
    k_f = 1.0000000000000002e-06 * 31990
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (16780.110000000001) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  16780.110000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* H2O */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[11];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[221] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    J[225] += dqdci;              /* dwdot[CH4]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[388] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[391] -= dqdci;              /* dwdot[H2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 65: CO + HO2 => CO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = 1.0000000000000002e-06 * 30100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (23000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  23000  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[12] -= q; /* CO */
    wdot[18] += q; /* CO2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[12];
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[292] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[298] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[8];
    J[423] += dqdci;              /* dwdot[OH]/d[CO] */
    J[428] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[438] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1202] -= dqdT;              /* dwdot[CO]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 66: CO + O2 => CO2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[12];
    k_f = 1.0000000000000002e-06 * 16200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (47700.050000000003) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  47700.050000000003  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[12] -= q; /* CO */
    wdot[18] += q; /* CO2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[12];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[187] -= dqdci;              /* dwdot[CO]/d[O2] */
    J[193] += dqdci;              /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[5];
    J[421] += dqdci;              /* dwdot[O]/d[CO] */
    J[425] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[438] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1202] -= dqdT;              /* dwdot[CO]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 67: CO2 + O => CO + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[18];
    k_f = 1.0000000000000002e-06 * 143300000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (53919.93) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  53919.93  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[5] += q; /* O2 */
    wdot[12] += q; /* CO */
    wdot[18] -= q; /* CO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[18];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[40] += dqdci;               /* dwdot[O2]/d[O] */
    J[47] += dqdci;               /* dwdot[CO]/d[O] */
    J[53] -= dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[1];
    J[631] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[635] += dqdci;              /* dwdot[O2]/d[CO2] */
    J[642] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[648] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1208] -= dqdT;              /* dwdot[CO2]/dT */

    /*reaction 68: CO + OH => CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[12];
    k_f = 1.0000000000000002e-06 * 140000
                * exp(1.95 * tc[0] - 0.50321666580471969 * (-1347.04) * invT);
    dlnkfdT = 1.95 * invT + 0.50321666580471969 *  -1347.04  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[12] -= q; /* CO */
    wdot[18] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[107] += dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[117] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[123] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3];
    J[422] += dqdci;              /* dwdot[H]/d[CO] */
    J[423] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[438] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1202] -= dqdT;              /* dwdot[CO]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 69: CO2 + H => CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = 1.0000000000000002e-06 * 15680000
                * exp(1.95 * tc[0] - 0.50321666580471969 * (20989.959999999999) * invT);
    dlnkfdT = 1.95 * invT + 0.50321666580471969 *  20989.959999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[12] += q; /* CO */
    wdot[18] -= q; /* CO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[OH]/d[H] */
    J[82] += dqdci;               /* dwdot[CO]/d[H] */
    J[88] -= dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[2];
    J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
    J[633] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[642] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[648] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1208] -= dqdT;              /* dwdot[CO2]/dT */

    /*reaction 70: HCO + CH3 => CH4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[19];
    k_f = 1.0000000000000002e-06 * 121000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[12] += q; /* CO */
    wdot[15] += q; /* CH4 */
    wdot[19] -= q; /* HCO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[19];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[397] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[404] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[11];
    J[676] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[680] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 71: HCO + H => CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[19];
    k_f = 1.0000000000000002e-06 * 73400000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[19];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[82] += dqdci;               /* dwdot[CO]/d[H] */
    J[89] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[667] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[669] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 72: HCO + O2 => CO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[19];
    k_f = 1.0000000000000002e-06 * 7580000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (409.88999999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  409.88999999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[19];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[187] += dqdci;              /* dwdot[CO]/d[O2] */
    J[194] -= dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[670] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[673] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 73: HCO + O => CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = 1.0000000000000002e-06 * 30200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[19];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[47] += dqdci;               /* dwdot[CO]/d[O] */
    J[54] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[666] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[668] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 74: HCO + O => CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[18] += q; /* CO2 */
    wdot[19] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[19];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[53] += dqdci;               /* dwdot[CO2]/d[O] */
    J[54] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[666] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[667] += dqdci;              /* dwdot[H]/d[HCO] */
    J[683] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 75: HCO + OH => CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[19];
    k_f = 1.0000000000000002e-06 * 102000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[19];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[117] += dqdci;              /* dwdot[CO]/d[OH] */
    J[124] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[668] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[671] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 76: CH2O + H => HCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[10];
    k_f = 1.0000000000000002e-06 * 933400000
                * exp(1.5 * tc[0] - 0.50321666580471969 * (2976.0999999999999) * invT);
    dlnkfdT = 1.5 * invT + 0.50321666580471969 *  2976.0999999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[80] -= dqdci;               /* dwdot[CH2O]/d[H] */
    J[89] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[2];
    J[352] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[354] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 77: CH2O + O2 => HCO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[10];
    k_f = 1.0000000000000002e-06 * 20500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (38950.050000000003) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  38950.050000000003  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[185] -= dqdci;              /* dwdot[CH2O]/d[O2] */
    J[194] += dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[5];
    J[355] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[358] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 78: CH2O + OH => HCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[10];
    k_f = 1.0000000000000002e-06 * 3430000000
                * exp(1.1799999999999999 * tc[0] - 0.50321666580471969 * (-446.94) * invT);
    dlnkfdT = 1.1799999999999999 * invT + 0.50321666580471969 *  -446.94  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[115] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    J[124] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[3];
    J[353] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[356] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 79: CH2O + HO2 => HCO + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = 1.0000000000000002e-06 * 0.0058199999999999997
                * exp(4.5300000000000002 * tc[0] - 0.50321666580471969 * (6556.8800000000001) * invT);
    dlnkfdT = 4.5300000000000002 * invT + 0.50321666580471969 *  6556.8800000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[10];
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[290] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[8];
    J[357] += dqdci;              /* dwdot[H2O2]/d[CH2O] */
    J[358] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 80: CH2O + O => HCO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = 1.0000000000000002e-06 * 416000000000
                * exp(0.56999999999999995 * tc[0] - 0.50321666580471969 * (2761.9499999999998) * invT);
    dlnkfdT = 0.56999999999999995 * invT + 0.50321666580471969 *  2761.9499999999998  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[45] -= dqdci;               /* dwdot[CH2O]/d[O] */
    J[54] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[1];
    J[351] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[353] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 81: CH2O + CH3 => HCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[11];
    k_f = 1.0000000000000002e-06 * 3.636e-06
                * exp(5.4199999999999999 * tc[0] - 0.50321666580471969 * (998.09000000000003) * invT);
    dlnkfdT = 5.4199999999999999 * invT + 0.50321666580471969 *  998.09000000000003  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] -= q; /* CH2O */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[19] += q; /* HCO */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[11];
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[361] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[365] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[10];
    J[395] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[404] += dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 82: CH3O + O2 => CH2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = 1.0000000000000002e-06 * 55000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (2424) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  2424  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[10] += q; /* CH2O */
    wdot[13] -= q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[188] -= dqdci;              /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[5];
    J[460] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[463] += dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[465] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[468] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1203] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 83: CH2GSG + CO2 => CH2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[18];
    k_f = 1.0000000000000002e-06 * 3000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] -= q; /* CH2GSG */
    wdot[10] += q; /* CH2O */
    wdot[12] += q; /* CO */
    wdot[18] -= q; /* CO2 */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[18];
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[325] += dqdci;              /* dwdot[CH2O]/d[CH2GSG] */
    J[327] += dqdci;              /* dwdot[CO]/d[CH2GSG] */
    J[333] -= dqdci;              /* dwdot[CO2]/d[CH2GSG] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[9];
    J[639] -= dqdci;              /* dwdot[CH2GSG]/d[CO2] */
    J[640] += dqdci;              /* dwdot[CH2O]/d[CO2] */
    J[642] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[648] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1208] -= dqdT;              /* dwdot[CO2]/dT */

    /*reaction 84: CH3O2 + CH3 => 2.000000 CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[20];
    k_f = 1.0000000000000002e-06 * 7000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1000  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[13] += 2 * q; /* CH3O */
    wdot[20] -= q; /* CH3O2 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[20];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[398] += 2 * dqdci;          /* dwdot[CH3O]/d[CH3] */
    J[405] -= dqdci;              /* dwdot[CH3O2]/d[CH3] */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[11];
    J[711] -= dqdci;              /* dwdot[CH3]/d[CH3O2] */
    J[713] += 2 * dqdci;          /* dwdot[CH3O]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1203] += 2 * dqdT;          /* dwdot[CH3O]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */

    /*reaction 85: CH3O2 + HO2 => CH3O2H + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[20];
    k_f = 1.0000000000000002e-06 * 17500000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-3275.0999999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -3275.0999999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[20] -= q; /* CH3O2 */
    wdot[21] += q; /* CH3O2H */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[20];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[300] -= dqdci;              /* dwdot[CH3O2]/d[HO2] */
    J[301] += dqdci;              /* dwdot[CH3O2H]/d[HO2] */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[8];
    J[705] += dqdci;              /* dwdot[O2]/d[CH3O2] */
    J[708] -= dqdci;              /* dwdot[HO2]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    J[721] += dqdci;              /* dwdot[CH3O2H]/d[CH3O2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */
    J[1211] += dqdT;              /* dwdot[CH3O2H]/dT */

    /*reaction 86: 2.000000 CH3O2 => O2 + 2.000000 CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[20], 2.000000);
    k_f = 1.0000000000000002e-06 * 14000000000000000
                * exp(-1.6100000000000001 * tc[0] - 0.50321666580471969 * (1859.9400000000001) * invT);
    dlnkfdT = -1.6100000000000001 * invT + 0.50321666580471969 *  1859.9400000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[13] += 2 * q; /* CH3O */
    wdot[20] -= 2 * q; /* CH3O2 */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*2.000000*sc[20];
    J[705] += dqdci;              /* dwdot[O2]/d[CH3O2] */
    J[713] += 2 * dqdci;          /* dwdot[CH3O]/d[CH3O2] */
    J[720] += -2 * dqdci;         /* dwdot[CH3O2]/d[CH3O2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1203] += 2 * dqdT;          /* dwdot[CH3O]/dT */
    J[1210] += -2 * dqdT;         /* dwdot[CH3O2]/dT */

    /*reaction 87: CH3O2 + CH2O => CH3O2H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[20];
    k_f = 1.0000000000000002e-06 * 1990000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (11659.889999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  11659.889999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    wdot[20] -= q; /* CH3O2 */
    wdot[21] += q; /* CH3O2H */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[20];
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[370] -= dqdci;              /* dwdot[CH3O2]/d[CH2O] */
    J[371] += dqdci;              /* dwdot[CH3O2H]/d[CH2O] */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[10];
    J[710] -= dqdci;              /* dwdot[CH2O]/d[CH3O2] */
    J[719] += dqdci;              /* dwdot[HCO]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    J[721] += dqdci;              /* dwdot[CH3O2H]/d[CH3O2] */
    /* d()/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */
    J[1211] += dqdT;              /* dwdot[CH3O2H]/dT */

    /*reaction 88: CH3O2H => CH3O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[21];
    k_f = 1 * 631000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (42299.949999999997) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  42299.949999999997  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[13] += q; /* CH3O */
    wdot[21] -= q; /* CH3O2H */
    /* d()/d[CH3O2H] */
    dqdci =  + k_f;
    J[738] += dqdci;              /* dwdot[OH]/d[CH3O2H] */
    J[748] += dqdci;              /* dwdot[CH3O]/d[CH3O2H] */
    J[756] -= dqdci;              /* dwdot[CH3O2H]/d[CH3O2H] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1203] += dqdT;              /* dwdot[CH3O]/dT */
    J[1211] -= dqdT;              /* dwdot[CH3O2H]/dT */

    /*reaction 89: CH3O + OH => CH3O2H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[13];
    k_f = 1.0000000000000002e-06 * 116600000000
                * exp(0.59999999999999998 * tc[0] - 0.50321666580471969 * (-1771.03) * invT);
    dlnkfdT = 0.59999999999999998 * invT + 0.50321666580471969 *  -1771.03  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[13] -= q; /* CH3O */
    wdot[21] += q; /* CH3O2H */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[118] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    J[126] += dqdci;              /* dwdot[CH3O2H]/d[OH] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[3];
    J[458] -= dqdci;              /* dwdot[OH]/d[CH3O] */
    J[468] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[476] += dqdci;              /* dwdot[CH3O2H]/d[CH3O] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1203] -= dqdT;              /* dwdot[CH3O]/dT */
    J[1211] += dqdT;              /* dwdot[CH3O2H]/dT */

    /*reaction 90: C2H2 + O2 => HCCO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[22];
    k_f = 1.0000000000000002e-06 * 200000000
                * exp(1.5 * tc[0] - 0.50321666580471969 * (30099.900000000001) * invT);
    dlnkfdT = 1.5 * invT + 0.50321666580471969 *  30099.900000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[22] -= q; /* C2H2 */
    wdot[23] += q; /* HCCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[22];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[197] -= dqdci;              /* dwdot[C2H2]/d[O2] */
    J[198] += dqdci;              /* dwdot[HCCO]/d[O2] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[5];
    J[773] += dqdci;              /* dwdot[OH]/d[C2H2] */
    J[775] -= dqdci;              /* dwdot[O2]/d[C2H2] */
    J[792] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[793] += dqdci;              /* dwdot[HCCO]/d[C2H2] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1212] -= dqdT;              /* dwdot[C2H2]/dT */
    J[1213] += dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 91: C2H2 + O => HCCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[22];
    k_f = 1.0000000000000002e-06 * 14300000
                * exp(2 * tc[0] - 0.50321666580471969 * (1900.0999999999999) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  1900.0999999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[22] -= q; /* C2H2 */
    wdot[23] += q; /* HCCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[22];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[57] -= dqdci;               /* dwdot[C2H2]/d[O] */
    J[58] += dqdci;               /* dwdot[HCCO]/d[O] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[1];
    J[771] -= dqdci;              /* dwdot[O]/d[C2H2] */
    J[772] += dqdci;              /* dwdot[H]/d[C2H2] */
    J[792] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[793] += dqdci;              /* dwdot[HCCO]/d[C2H2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1212] -= dqdT;              /* dwdot[C2H2]/dT */
    J[1213] += dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 92: C2H3 + O2 => CH2CHO + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[24];
    k_f = 1.0000000000000002e-06 * 350000000000000
                * exp(-0.60999999999999999 * tc[0] - 0.50321666580471969 * (5260.04) * invT);
    dlnkfdT = -0.60999999999999999 * invT + 0.50321666580471969 *  5260.04  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[24] -= q; /* C2H3 */
    wdot[25] += q; /* CH2CHO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[24];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[199] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    J[200] += dqdci;              /* dwdot[CH2CHO]/d[O2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[5];
    J[841] += dqdci;              /* dwdot[O]/d[C2H3] */
    J[845] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[865] += dqdci;              /* dwdot[CH2CHO]/d[C2H3] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */
    J[1215] += dqdT;              /* dwdot[CH2CHO]/dT */

    /*reaction 93: C2H3 + H => C2H2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[24];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (2500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  2500  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[22] += q; /* C2H2 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[24];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[92] += dqdci;               /* dwdot[C2H2]/d[H] */
    J[94] -= dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[2];
    J[842] -= dqdci;              /* dwdot[H]/d[C2H3] */
    J[844] += dqdci;              /* dwdot[H2]/d[C2H3] */
    J[862] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1212] += dqdT;              /* dwdot[C2H2]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 94: C2H3 + O2 => CH2O + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[24];
    k_f = 1.0000000000000002e-06 * 1.6999999999999999e+29
                * exp(-5.3099999999999996 * tc[0] - 0.50321666580471969 * (6500) * invT);
    dlnkfdT = -5.3099999999999996 * invT + 0.50321666580471969 *  6500  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[19] += q; /* HCO */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[24];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[194] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[199] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[5];
    J[845] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[850] += dqdci;              /* dwdot[CH2O]/d[C2H3] */
    J[859] += dqdci;              /* dwdot[HCO]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 95: C2H3 + O2 => C2H2 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[24];
    k_f = 1.0000000000000002e-06 * 2.12e-06
                * exp(6 * tc[0] - 0.50321666580471969 * (9483.9899999999998) * invT);
    dlnkfdT = 6 * invT + 0.50321666580471969 *  9483.9899999999998  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[22] += q; /* C2H2 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[24];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[197] += dqdci;              /* dwdot[C2H2]/d[O2] */
    J[199] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[5];
    J[845] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[848] += dqdci;              /* dwdot[HO2]/d[C2H3] */
    J[862] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1212] += dqdT;              /* dwdot[C2H2]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 96: C2H4 + O => CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = 1.0000000000000002e-06 * 10200000
                * exp(1.8799999999999999 * tc[0] - 0.50321666580471969 * (179.02000000000001) * invT);
    dlnkfdT = 1.8799999999999999 * invT + 0.50321666580471969 *  179.02000000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[11] += q; /* CH3 */
    wdot[16] -= q; /* C2H4 */
    wdot[19] += q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[46] += dqdci;               /* dwdot[CH3]/d[O] */
    J[51] -= dqdci;               /* dwdot[C2H4]/d[O] */
    J[54] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[1];
    J[561] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[571] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[579] += dqdci;              /* dwdot[HCO]/d[C2H4] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 97: C2H4 + OH => C2H3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[16];
    k_f = 1.0000000000000002e-06 * 20500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (5950.0500000000002) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  5950.0500000000002  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[16] -= q; /* C2H4 */
    wdot[24] += q; /* C2H3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[16];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[121] -= dqdci;              /* dwdot[C2H4]/d[OH] */
    J[129] += dqdci;              /* dwdot[C2H3]/d[OH] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[3];
    J[563] -= dqdci;              /* dwdot[OH]/d[C2H4] */
    J[566] += dqdci;              /* dwdot[H2O]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[584] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 98: C2H3 + H2O => C2H4 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[24];
    k_f = 1.0000000000000002e-06 * 6033000000000000
                * exp(-0.82999999999999996 * tc[0] - 0.50321666580471969 * (21760.040000000001) * invT);
    dlnkfdT = -0.82999999999999996 * invT + 0.50321666580471969 *  21760.040000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* H2O */
    wdot[16] += q; /* C2H4 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[24];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[226] += dqdci;              /* dwdot[C2H4]/d[H2O] */
    J[234] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[6];
    J[843] += dqdci;              /* dwdot[OH]/d[C2H3] */
    J[846] -= dqdci;              /* dwdot[H2O]/d[C2H3] */
    J[856] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 99: C2H4 + H => C2H3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[16];
    k_f = 1.0000000000000002e-06 * 0.0084200000000000004
                * exp(4.6200000000000001 * tc[0] - 0.50321666580471969 * (2582.9299999999998) * invT);
    dlnkfdT = 4.6200000000000001 * invT + 0.50321666580471969 *  2582.9299999999998  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[16] -= q; /* C2H4 */
    wdot[24] += q; /* C2H3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[16];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[86] -= dqdci;               /* dwdot[C2H4]/d[H] */
    J[94] += dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[2];
    J[562] -= dqdci;              /* dwdot[H]/d[C2H4] */
    J[564] += dqdci;              /* dwdot[H2]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[584] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 100: C2H3 + H2 => C2H4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[24];
    k_f = 1.0000000000000002e-06 * 0.57230000000000003
                * exp(3.79 * tc[0] - 0.50321666580471969 * (3233.0300000000002) * invT);
    dlnkfdT = 3.79 * invT + 0.50321666580471969 *  3233.0300000000002  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[16] += q; /* C2H4 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[24];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[156] += dqdci;              /* dwdot[C2H4]/d[H2] */
    J[164] -= dqdci;              /* dwdot[C2H3]/d[H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[4];
    J[842] += dqdci;              /* dwdot[H]/d[C2H3] */
    J[844] -= dqdci;              /* dwdot[H2]/d[C2H3] */
    J[856] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 101: C2H4 + CH3O2 => C2H3 + CH3O2H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[16]*sc[20];
    k_f = 1.0000000000000002e-06 * 2230000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (17190.009999999998) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  17190.009999999998  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[16] -= q; /* C2H4 */
    wdot[20] -= q; /* CH3O2 */
    wdot[21] += q; /* CH3O2H */
    wdot[24] += q; /* C2H3 */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[20];
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[580] -= dqdci;              /* dwdot[CH3O2]/d[C2H4] */
    J[581] += dqdci;              /* dwdot[CH3O2H]/d[C2H4] */
    J[584] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[16];
    J[716] -= dqdci;              /* dwdot[C2H4]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    J[721] += dqdci;              /* dwdot[CH3O2H]/d[CH3O2] */
    J[724] += dqdci;              /* dwdot[C2H3]/d[CH3O2] */
    /* d()/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */
    J[1211] += dqdT;              /* dwdot[CH3O2H]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 102: C2H4 + CH3 => C2H3 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[16];
    k_f = 1.0000000000000002e-06 * 6.6200000000000001
                * exp(3.7000000000000002 * tc[0] - 0.50321666580471969 * (9500) * invT);
    dlnkfdT = 3.7000000000000002 * invT + 0.50321666580471969 *  9500  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[16] -= q; /* C2H4 */
    wdot[24] += q; /* C2H3 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[401] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[409] += dqdci;              /* dwdot[C2H3]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[11];
    J[571] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[575] += dqdci;              /* dwdot[CH4]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[584] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 103: C2H3 + CH4 => C2H4 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[15]*sc[24];
    k_f = 1.0000000000000002e-06 * 1.4399999999999999
                * exp(4.0199999999999996 * tc[0] - 0.50321666580471969 * (5472.04) * invT);
    dlnkfdT = 4.0199999999999996 * invT + 0.50321666580471969 *  5472.04  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] += q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    wdot[16] += q; /* C2H4 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[24];
    J[536] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    J[541] += dqdci;              /* dwdot[C2H4]/d[CH4] */
    J[549] -= dqdci;              /* dwdot[C2H3]/d[CH4] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[15];
    J[851] += dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[855] -= dqdci;              /* dwdot[CH4]/d[C2H3] */
    J[856] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 104: C2H4 + O => CH2CHO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = 1.0000000000000002e-06 * 3390000
                * exp(1.8799999999999999 * tc[0] - 0.50321666580471969 * (179.02000000000001) * invT);
    dlnkfdT = 1.8799999999999999 * invT + 0.50321666580471969 *  179.02000000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[16] -= q; /* C2H4 */
    wdot[25] += q; /* CH2CHO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[51] -= dqdci;               /* dwdot[C2H4]/d[O] */
    J[60] += dqdci;               /* dwdot[CH2CHO]/d[O] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[1];
    J[561] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[562] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[585] += dqdci;              /* dwdot[CH2CHO]/d[C2H4] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1215] += dqdT;              /* dwdot[CH2CHO]/dT */

    /*reaction 105: CH3 + C2H5 => CH4 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[14];
    k_f = 1.0000000000000002e-06 * 19500000000000
                * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.5 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[14] -= q; /* C2H5 */
    wdot[15] += q; /* CH4 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[14];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[399] -= dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[401] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[11];
    J[501] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[505] += dqdci;              /* dwdot[CH4]/d[C2H5] */
    J[506] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 106: C2H5 + O2 => C2H4 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[14];
    k_f = 1.0000000000000002e-06 * 1.22e+30
                * exp(-5.7599999999999998 * tc[0] - 0.50321666580471969 * (10099.9) * invT);
    dlnkfdT = -5.7599999999999998 * invT + 0.50321666580471969 *  10099.9  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[14] -= q; /* C2H5 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[189] -= dqdci;              /* dwdot[C2H5]/d[O2] */
    J[191] += dqdci;              /* dwdot[C2H4]/d[O2] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[5];
    J[495] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[498] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[506] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 107: C2H4 + HO2 => C2H5 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[16];
    k_f = 1.0000000000000002e-06 * 1.259e+30
                * exp(-5.6299999999999999 * tc[0] - 0.50321666580471969 * (22299.950000000001) * invT);
    dlnkfdT = -5.6299999999999999 * invT + 0.50321666580471969 *  22299.950000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[14] += q; /* C2H5 */
    wdot[16] -= q; /* C2H4 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[16];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[294] += dqdci;              /* dwdot[C2H5]/d[HO2] */
    J[296] -= dqdci;              /* dwdot[C2H4]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[8];
    J[565] += dqdci;              /* dwdot[O2]/d[C2H4] */
    J[568] -= dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[574] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 108: C2H5 + HO2 => C2H6 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[14];
    k_f = 1.0000000000000002e-06 * 267900000
                * exp(0.89000000000000001 * tc[0] - 0.50321666580471969 * (-1922.0799999999999) * invT);
    dlnkfdT = 0.89000000000000001 * invT + 0.50321666580471969 *  -1922.0799999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[14] -= q; /* C2H5 */
    wdot[17] += q; /* C2H6 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[14];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[294] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    J[297] += dqdci;              /* dwdot[C2H6]/d[HO2] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[8];
    J[495] += dqdci;              /* dwdot[O2]/d[C2H5] */
    J[498] -= dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[507] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1207] += dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 109: H + C2H5 => C2H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[14];
    k_f = 1.0000000000000002e-06 * 583100000000
                * exp(0.59899999999999998 * tc[0] - 0.50321666580471969 * (-2913) * invT);
    dlnkfdT = 0.59899999999999998 * invT + 0.50321666580471969 *  -2913  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[14] -= q; /* C2H5 */
    wdot[17] += q; /* C2H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[14];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[84] -= dqdci;               /* dwdot[C2H5]/d[H] */
    J[87] += dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[2];
    J[492] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[507] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1207] += dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 110: C2H6 + H => C2H5 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = 1.0000000000000002e-06 * 554
                * exp(3.5 * tc[0] - 0.50321666580471969 * (5167.0699999999997) * invT);
    dlnkfdT = 3.5 * invT + 0.50321666580471969 *  5167.0699999999997  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[14] += q; /* C2H5 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[17];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[84] += dqdci;               /* dwdot[C2H5]/d[H] */
    J[87] -= dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[2];
    J[597] -= dqdci;              /* dwdot[H]/d[C2H6] */
    J[599] += dqdci;              /* dwdot[H2]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 111: C2H5 + H2 => C2H6 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[14];
    k_f = 1.0000000000000002e-06 * 0.13550000000000001
                * exp(4.0599999999999996 * tc[0] - 0.50321666580471969 * (8857.0699999999997) * invT);
    dlnkfdT = 4.0599999999999996 * invT + 0.50321666580471969 *  8857.0699999999997  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[14] -= q; /* C2H5 */
    wdot[17] += q; /* C2H6 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[14];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[154] -= dqdci;              /* dwdot[C2H5]/d[H2] */
    J[157] += dqdci;              /* dwdot[C2H6]/d[H2] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[4];
    J[492] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[494] -= dqdci;              /* dwdot[H2]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[507] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1207] += dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 112: C2H6 + OH => C2H5 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[17];
    k_f = 1.0000000000000002e-06 * 58000000
                * exp(1.73 * tc[0] - 0.50321666580471969 * (1159.8900000000001) * invT);
    dlnkfdT = 1.73 * invT + 0.50321666580471969 *  1159.8900000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[14] += q; /* C2H5 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[17];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[119] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[122] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[3];
    J[598] -= dqdci;              /* dwdot[OH]/d[C2H6] */
    J[601] += dqdci;              /* dwdot[H2O]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 113: CH2GSG + C2H6 => CH3 + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[17];
    k_f = 1.0000000000000002e-06 * 120000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] += q; /* CH3 */
    wdot[14] += q; /* C2H5 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[17];
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] += dqdci;              /* dwdot[CH3]/d[CH2GSG] */
    J[329] += dqdci;              /* dwdot[C2H5]/d[CH2GSG] */
    J[332] -= dqdci;              /* dwdot[C2H6]/d[CH2GSG] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[9];
    J[604] -= dqdci;              /* dwdot[CH2GSG]/d[C2H6] */
    J[606] += dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 114: C2H6 + O => C2H5 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[17];
    k_f = 1.0000000000000002e-06 * 13000000
                * exp(2.1299999999999999 * tc[0] - 0.50321666580471969 * (5190.0100000000002) * invT);
    dlnkfdT = 2.1299999999999999 * invT + 0.50321666580471969 *  5190.0100000000002  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[14] += q; /* C2H5 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[49] += dqdci;               /* dwdot[C2H5]/d[O] */
    J[52] -= dqdci;               /* dwdot[C2H6]/d[O] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[1];
    J[596] -= dqdci;              /* dwdot[O]/d[C2H6] */
    J[598] += dqdci;              /* dwdot[OH]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 115: C2H6 + CH3 => C2H5 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[17];
    k_f = 1.0000000000000002e-06 * 1.5099999999999999e-07
                * exp(6 * tc[0] - 0.50321666580471969 * (6047.0799999999999) * invT);
    dlnkfdT = 6 * invT + 0.50321666580471969 *  6047.0799999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[14] += q; /* C2H5 */
    wdot[15] += q; /* CH4 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[17];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[399] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[402] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[11];
    J[606] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[610] += dqdci;              /* dwdot[CH4]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 116: HCCO + O => H + 2.000000 CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[23];
    k_f = 1.0000000000000002e-06 * 80000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[12] += 2 * q; /* CO */
    wdot[23] -= q; /* HCCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[23];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[47] += 2 * dqdci;           /* dwdot[CO]/d[O] */
    J[58] -= dqdci;               /* dwdot[HCCO]/d[O] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[1];
    J[806] -= dqdci;              /* dwdot[O]/d[HCCO] */
    J[807] += dqdci;              /* dwdot[H]/d[HCCO] */
    J[817] += 2 * dqdci;          /* dwdot[CO]/d[HCCO] */
    J[828] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1202] += 2 * dqdT;          /* dwdot[CO]/dT */
    J[1213] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 117: HCCO + O2 => CO2 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[23];
    k_f = 1.0000000000000002e-06 * 240000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-853.97000000000003) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -853.97000000000003  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[18] += q; /* CO2 */
    wdot[19] += q; /* HCO */
    wdot[23] -= q; /* HCCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[23];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[193] += dqdci;              /* dwdot[CO2]/d[O2] */
    J[194] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[198] -= dqdci;              /* dwdot[HCCO]/d[O2] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[5];
    J[810] -= dqdci;              /* dwdot[O2]/d[HCCO] */
    J[823] += dqdci;              /* dwdot[CO2]/d[HCCO] */
    J[824] += dqdci;              /* dwdot[HCO]/d[HCCO] */
    J[828] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1213] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 118: HCCO + OH => 2.000000 HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[23];
    k_f = 1.0000000000000002e-06 * 10000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[19] += 2 * q; /* HCO */
    wdot[23] -= q; /* HCCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[23];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[124] += 2 * dqdci;          /* dwdot[HCO]/d[OH] */
    J[128] -= dqdci;              /* dwdot[HCCO]/d[OH] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[3];
    J[808] -= dqdci;              /* dwdot[OH]/d[HCCO] */
    J[824] += 2 * dqdci;          /* dwdot[HCO]/d[HCCO] */
    J[828] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1209] += 2 * dqdT;          /* dwdot[HCO]/dT */
    J[1213] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 119: HCCO + H => CH2GSG + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[23];
    k_f = 1.0000000000000002e-06 * 110000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[9] += q; /* CH2GSG */
    wdot[12] += q; /* CO */
    wdot[23] -= q; /* HCCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[23];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[79] += dqdci;               /* dwdot[CH2GSG]/d[H] */
    J[82] += dqdci;               /* dwdot[CO]/d[H] */
    J[93] -= dqdci;               /* dwdot[HCCO]/d[H] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[2];
    J[807] -= dqdci;              /* dwdot[H]/d[HCCO] */
    J[814] += dqdci;              /* dwdot[CH2GSG]/d[HCCO] */
    J[817] += dqdci;              /* dwdot[CO]/d[HCCO] */
    J[828] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1199] += dqdT;              /* dwdot[CH2GSG]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1213] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 120: CH2GSG + CO => HCCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[12];
    k_f = 1.0000000000000002e-06 * 2046000000000
                * exp(0.89000000000000001 * tc[0] - 0.50321666580471969 * (27830.07) * invT);
    dlnkfdT = 0.89000000000000001 * invT + 0.50321666580471969 *  27830.07  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[9] -= q; /* CH2GSG */
    wdot[12] -= q; /* CO */
    wdot[23] += q; /* HCCO */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[12];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[327] -= dqdci;              /* dwdot[CO]/d[CH2GSG] */
    J[338] += dqdci;              /* dwdot[HCCO]/d[CH2GSG] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[9];
    J[422] += dqdci;              /* dwdot[H]/d[CO] */
    J[429] -= dqdci;              /* dwdot[CH2GSG]/d[CO] */
    J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[443] += dqdci;              /* dwdot[HCCO]/d[CO] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1202] -= dqdT;              /* dwdot[CO]/dT */
    J[1213] += dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 121: CH2CHO + O2 => CH2O + CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[25];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (4200.0500000000002) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  4200.0500000000002  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[12] += q; /* CO */
    wdot[25] -= q; /* CH2CHO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[25];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[187] += dqdci;              /* dwdot[CO]/d[O2] */
    J[200] -= dqdci;              /* dwdot[CH2CHO]/d[O2] */
    /* d()/d[CH2CHO] */
    dqdci =  + k_f*sc[5];
    J[878] += dqdci;              /* dwdot[OH]/d[CH2CHO] */
    J[880] -= dqdci;              /* dwdot[O2]/d[CH2CHO] */
    J[885] += dqdci;              /* dwdot[CH2O]/d[CH2CHO] */
    J[887] += dqdci;              /* dwdot[CO]/d[CH2CHO] */
    J[900] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1215] -= dqdT;              /* dwdot[CH2CHO]/dT */

    /*reaction 122: C3H5XA + CH2O => C3H6 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[26];
    k_f = 1.0000000000000002e-06 * 630000000
                * exp(1.8999999999999999 * tc[0] - 0.50321666580471969 * (18190.009999999998) * invT);
    dlnkfdT = 1.8999999999999999 * invT + 0.50321666580471969 *  18190.009999999998  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    wdot[26] -= q; /* C3H5XA */
    wdot[27] += q; /* C3H6 */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[26];
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[376] -= dqdci;              /* dwdot[C3H5XA]/d[CH2O] */
    J[377] += dqdci;              /* dwdot[C3H6]/d[CH2O] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[10];
    J[920] -= dqdci;              /* dwdot[CH2O]/d[C3H5XA] */
    J[929] += dqdci;              /* dwdot[HCO]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[937] += dqdci;              /* dwdot[C3H6]/d[C3H5XA] */
    /* d()/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 123: C3H6 + HCO => C3H5XA + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[19]*sc[27];
    k_f = 1.0000000000000002e-06 * 109700000
                * exp(1.8899999999999999 * tc[0] - 0.50321666580471969 * (15840.110000000001) * invT);
    dlnkfdT = 1.8899999999999999 * invT + 0.50321666580471969 *  15840.110000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] += q; /* CH2O */
    wdot[19] -= q; /* HCO */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[27];
    J[675] += dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[691] += dqdci;              /* dwdot[C3H5XA]/d[HCO] */
    J[692] -= dqdci;              /* dwdot[C3H6]/d[HCO] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[19];
    J[955] += dqdci;              /* dwdot[CH2O]/d[C3H6] */
    J[964] -= dqdci;              /* dwdot[HCO]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 124: C3H5XA + HO2 => C3H5O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[26];
    k_f = 1.0000000000000002e-06 * 7000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1000  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[26] -= q; /* C3H5XA */
    wdot[28] += q; /* C3H5O */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[26];
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[306] -= dqdci;              /* dwdot[C3H5XA]/d[HO2] */
    J[308] += dqdci;              /* dwdot[C3H5O]/d[HO2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[8];
    J[913] += dqdci;              /* dwdot[OH]/d[C3H5XA] */
    J[918] -= dqdci;              /* dwdot[HO2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[938] += dqdci;              /* dwdot[C3H5O]/d[C3H5XA] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1218] += dqdT;              /* dwdot[C3H5O]/dT */

    /*reaction 125: C3H5XA + CH3O2 => C3H5O + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[20]*sc[26];
    k_f = 1.0000000000000002e-06 * 7000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1000  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[13] += q; /* CH3O */
    wdot[20] -= q; /* CH3O2 */
    wdot[26] -= q; /* C3H5XA */
    wdot[28] += q; /* C3H5O */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[26];
    J[713] += dqdci;              /* dwdot[CH3O]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    J[726] -= dqdci;              /* dwdot[C3H5XA]/d[CH3O2] */
    J[728] += dqdci;              /* dwdot[C3H5O]/d[CH3O2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[20];
    J[923] += dqdci;              /* dwdot[CH3O]/d[C3H5XA] */
    J[930] -= dqdci;              /* dwdot[CH3O2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[938] += dqdci;              /* dwdot[C3H5O]/d[C3H5XA] */
    /* d()/dT */
    J[1203] += dqdT;              /* dwdot[CH3O]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1218] += dqdT;              /* dwdot[C3H5O]/dT */

    /*reaction 126: C3H5XA + HO2 => C3H6 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[26];
    k_f = 1.0000000000000002e-06 * 33320000000
                * exp(0.34000000000000002 * tc[0] - 0.50321666580471969 * (-555.92999999999995) * invT);
    dlnkfdT = 0.34000000000000002 * invT + 0.50321666580471969 *  -555.92999999999995  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[26] -= q; /* C3H5XA */
    wdot[27] += q; /* C3H6 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[26];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[306] -= dqdci;              /* dwdot[C3H5XA]/d[HO2] */
    J[307] += dqdci;              /* dwdot[C3H6]/d[HO2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[8];
    J[915] += dqdci;              /* dwdot[O2]/d[C3H5XA] */
    J[918] -= dqdci;              /* dwdot[HO2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[937] += dqdci;              /* dwdot[C3H6]/d[C3H5XA] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 127: C3H5XA + O2 => CH2CHO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[26];
    k_f = 1.0000000000000002e-06 * 7140000000000000
                * exp(-1.21 * tc[0] - 0.50321666580471969 * (21049.950000000001) * invT);
    dlnkfdT = -1.21 * invT + 0.50321666580471969 *  21049.950000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[25] += q; /* CH2CHO */
    wdot[26] -= q; /* C3H5XA */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[26];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[200] += dqdci;              /* dwdot[CH2CHO]/d[O2] */
    J[201] -= dqdci;              /* dwdot[C3H5XA]/d[O2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[5];
    J[915] -= dqdci;              /* dwdot[O2]/d[C3H5XA] */
    J[920] += dqdci;              /* dwdot[CH2O]/d[C3H5XA] */
    J[935] += dqdci;              /* dwdot[CH2CHO]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1215] += dqdT;              /* dwdot[CH2CHO]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */

    /*reaction 128: C3H5XA + O2 => C2H2 + CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[26];
    k_f = 1.0000000000000002e-06 * 9.7200000000000003e+29
                * exp(-5.71 * tc[0] - 0.50321666580471969 * (21450.049999999999) * invT);
    dlnkfdT = -5.71 * invT + 0.50321666580471969 *  21450.049999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[22] += q; /* C2H2 */
    wdot[26] -= q; /* C3H5XA */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[26];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[197] += dqdci;              /* dwdot[C2H2]/d[O2] */
    J[201] -= dqdci;              /* dwdot[C3H5XA]/d[O2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[5];
    J[913] += dqdci;              /* dwdot[OH]/d[C3H5XA] */
    J[915] -= dqdci;              /* dwdot[O2]/d[C3H5XA] */
    J[920] += dqdci;              /* dwdot[CH2O]/d[C3H5XA] */
    J[932] += dqdci;              /* dwdot[C2H2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1212] += dqdT;              /* dwdot[C2H2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */

    /*reaction 129: C3H5XA + H => C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[26];
    k_f = 1.0000000000000002e-06 * 4.8869999999999999e+56
                * exp(-12.25 * tc[0] - 0.50321666580471969 * (28080.07) * invT);
    dlnkfdT = -12.25 * invT + 0.50321666580471969 *  28080.07  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[26] -= q; /* C3H5XA */
    wdot[27] += q; /* C3H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[26];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[96] -= dqdci;               /* dwdot[C3H5XA]/d[H] */
    J[97] += dqdci;               /* dwdot[C3H6]/d[H] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[2];
    J[912] -= dqdci;              /* dwdot[H]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[937] += dqdci;              /* dwdot[C3H6]/d[C3H5XA] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 130: C3H5XA => C2H2 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[26];
    k_f = 1 * 2.397e+48
                * exp(-9.9000000000000004 * tc[0] - 0.50321666580471969 * (82080.070000000007) * invT);
    dlnkfdT = -9.9000000000000004 * invT + 0.50321666580471969 *  82080.070000000007  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] += q; /* CH3 */
    wdot[22] += q; /* C2H2 */
    wdot[26] -= q; /* C3H5XA */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f;
    J[921] += dqdci;              /* dwdot[CH3]/d[C3H5XA] */
    J[932] += dqdci;              /* dwdot[C2H2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    /* d()/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1212] += dqdT;              /* dwdot[C2H2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */

    /*reaction 131: C3H6 => C2H3 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[27];
    k_f = 1 * 2.7300000000000001e+62
                * exp(-13.279999999999999 * tc[0] - 0.50321666580471969 * (123200.05) * invT);
    dlnkfdT = -13.279999999999999 * invT + 0.50321666580471969 *  123200.05  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] += q; /* CH3 */
    wdot[24] += q; /* C2H3 */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[C3H6] */
    dqdci =  + k_f;
    J[956] += dqdci;              /* dwdot[CH3]/d[C3H6] */
    J[969] += dqdci;              /* dwdot[C2H3]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 132: C2H3 + CH3 => C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[24];
    k_f = 1.0000000000000002e-06 * 4.7119999999999996e+59
                * exp(-13.19 * tc[0] - 0.50321666580471969 * (29539.91) * invT);
    dlnkfdT = -13.19 * invT + 0.50321666580471969 *  29539.91  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[24] -= q; /* C2H3 */
    wdot[27] += q; /* C3H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[24];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] -= dqdci;              /* dwdot[C2H3]/d[CH3] */
    J[412] += dqdci;              /* dwdot[C3H6]/d[CH3] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[11];
    J[851] -= dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[867] += dqdci;              /* dwdot[C3H6]/d[C2H3] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 133: C3H6 + O => C2H5 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[27];
    k_f = 1.0000000000000002e-06 * 15800000
                * exp(1.76 * tc[0] - 0.50321666580471969 * (-1216.0599999999999) * invT);
    dlnkfdT = 1.76 * invT + 0.50321666580471969 *  -1216.0599999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[14] += q; /* C2H5 */
    wdot[19] += q; /* HCO */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[27];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[49] += dqdci;               /* dwdot[C2H5]/d[O] */
    J[54] += dqdci;               /* dwdot[HCO]/d[O] */
    J[62] -= dqdci;               /* dwdot[C3H6]/d[O] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[1];
    J[946] -= dqdci;              /* dwdot[O]/d[C3H6] */
    J[959] += dqdci;              /* dwdot[C2H5]/d[C3H6] */
    J[964] += dqdci;              /* dwdot[HCO]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 134: C3H6 + H => C2H4 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = 1.0000000000000002e-06 * 4.8299999999999998e+33
                * exp(-5.8099999999999996 * tc[0] - 0.50321666580471969 * (18500) * invT);
    dlnkfdT = -5.8099999999999996 * invT + 0.50321666580471969 *  18500  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[11] += q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[81] += dqdci;               /* dwdot[CH3]/d[H] */
    J[86] += dqdci;               /* dwdot[C2H4]/d[H] */
    J[97] -= dqdci;               /* dwdot[C3H6]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[2];
    J[947] -= dqdci;              /* dwdot[H]/d[C3H6] */
    J[956] += dqdci;              /* dwdot[CH3]/d[C3H6] */
    J[961] += dqdci;              /* dwdot[C2H4]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 135: C2H4 + CH3 => C3H6 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[16];
    k_f = 1.0000000000000002e-06 * 2.3130000000000001e+33
                * exp(-5.9000000000000004 * tc[0] - 0.50321666580471969 * (31619.98) * invT);
    dlnkfdT = -5.9000000000000004 * invT + 0.50321666580471969 *  31619.98  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[11] -= q; /* CH3 */
    wdot[16] -= q; /* C2H4 */
    wdot[27] += q; /* C3H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[401] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[412] += dqdci;              /* dwdot[C3H6]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[11];
    J[562] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[571] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[587] += dqdci;              /* dwdot[C3H6]/d[C2H4] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 136: C3H6 + H => C3H5XA + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = 1.0000000000000002e-06 * 173000
                * exp(2.5 * tc[0] - 0.50321666580471969 * (2492.1100000000001) * invT);
    dlnkfdT = 2.5 * invT + 0.50321666580471969 *  2492.1100000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[96] += dqdci;               /* dwdot[C3H5XA]/d[H] */
    J[97] -= dqdci;               /* dwdot[C3H6]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[2];
    J[947] -= dqdci;              /* dwdot[H]/d[C3H6] */
    J[949] += dqdci;              /* dwdot[H2]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 137: C3H5XA + H2 => C3H6 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[26];
    k_f = 1.0000000000000002e-06 * 79330
                * exp(2.5099999999999998 * tc[0] - 0.50321666580471969 * (19520.080000000002) * invT);
    dlnkfdT = 2.5099999999999998 * invT + 0.50321666580471969 *  19520.080000000002  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[26] -= q; /* C3H5XA */
    wdot[27] += q; /* C3H6 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[26];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[166] -= dqdci;              /* dwdot[C3H5XA]/d[H2] */
    J[167] += dqdci;              /* dwdot[C3H6]/d[H2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[4];
    J[912] += dqdci;              /* dwdot[H]/d[C3H5XA] */
    J[914] -= dqdci;              /* dwdot[H2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[937] += dqdci;              /* dwdot[C3H6]/d[C3H5XA] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 138: C3H6 + OH => C3H5XA + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[27];
    k_f = 1.0000000000000002e-06 * 3120000
                * exp(2 * tc[0] - 0.50321666580471969 * (-298.04000000000002) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  -298.04000000000002  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[27];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[131] += dqdci;              /* dwdot[C3H5XA]/d[OH] */
    J[132] -= dqdci;              /* dwdot[C3H6]/d[OH] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[3];
    J[948] -= dqdci;              /* dwdot[OH]/d[C3H6] */
    J[951] += dqdci;              /* dwdot[H2O]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 139: C3H6 + O => C3H5XA + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[27];
    k_f = 1.0000000000000002e-06 * 524000000000
                * exp(0.69999999999999996 * tc[0] - 0.50321666580471969 * (5884.0799999999999) * invT);
    dlnkfdT = 0.69999999999999996 * invT + 0.50321666580471969 *  5884.0799999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[27];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[61] += dqdci;               /* dwdot[C3H5XA]/d[O] */
    J[62] -= dqdci;               /* dwdot[C3H6]/d[O] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[1];
    J[946] -= dqdci;              /* dwdot[O]/d[C3H6] */
    J[948] += dqdci;              /* dwdot[OH]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 140: C3H6 + CH3 => C3H5XA + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[27];
    k_f = 1.0000000000000002e-06 * 2.21
                * exp(3.5 * tc[0] - 0.50321666580471969 * (5674.9499999999998) * invT);
    dlnkfdT = 3.5 * invT + 0.50321666580471969 *  5674.9499999999998  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[27];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[411] += dqdci;              /* dwdot[C3H5XA]/d[CH3] */
    J[412] -= dqdci;              /* dwdot[C3H6]/d[CH3] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[11];
    J[956] -= dqdci;              /* dwdot[CH3]/d[C3H6] */
    J[960] += dqdci;              /* dwdot[CH4]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 141: IXC3H7 + O2 => C3H6 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[29];
    k_f = 1.0000000000000002e-06 * 450000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (5020.0799999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  5020.0799999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[27] += q; /* C3H6 */
    wdot[29] -= q; /* IXC3H7 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[29];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[202] += dqdci;              /* dwdot[C3H6]/d[O2] */
    J[204] -= dqdci;              /* dwdot[IXC3H7]/d[O2] */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f*sc[5];
    J[1020] -= dqdci;             /* dwdot[O2]/d[IXC3H7] */
    J[1023] += dqdci;             /* dwdot[HO2]/d[IXC3H7] */
    J[1042] += dqdci;             /* dwdot[C3H6]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */

    /*reaction 142: IXC3H7 => H + C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[29];
    k_f = 1 * 8.569e+18
                * exp(-1.5700000000000001 * tc[0] - 0.50321666580471969 * (40340.110000000001) * invT);
    dlnkfdT = -1.5700000000000001 * invT + 0.50321666580471969 *  40340.110000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[27] += q; /* C3H6 */
    wdot[29] -= q; /* IXC3H7 */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f;
    J[1017] += dqdci;             /* dwdot[H]/d[IXC3H7] */
    J[1042] += dqdci;             /* dwdot[C3H6]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */

    /*reaction 143: H + C3H6 => IXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = 1.0000000000000002e-06 * 13000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (1559.99) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  1559.99  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[27] -= q; /* C3H6 */
    wdot[29] += q; /* IXC3H7 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[97] -= dqdci;               /* dwdot[C3H6]/d[H] */
    J[99] += dqdci;               /* dwdot[IXC3H7]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[2];
    J[947] -= dqdci;              /* dwdot[H]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    J[974] += dqdci;              /* dwdot[IXC3H7]/d[C3H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */

    /*reaction 144: NXC3H7 => CH3 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[30];
    k_f = 1 * 228400000000000
                * exp(-0.55000000000000004 * tc[0] - 0.50321666580471969 * (28400.099999999999) * invT);
    dlnkfdT = -0.55000000000000004 * invT + 0.50321666580471969 *  28400.099999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] += q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    wdot[30] -= q; /* NXC3H7 */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f;
    J[1061] += dqdci;             /* dwdot[CH3]/d[NXC3H7] */
    J[1066] += dqdci;             /* dwdot[C2H4]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    /* d()/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 145: CH3 + C2H4 => NXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[16];
    k_f = 1.0000000000000002e-06 * 410000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7204.1099999999997) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  7204.1099999999997  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[16] -= q; /* C2H4 */
    wdot[30] += q; /* NXC3H7 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[401] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[415] += dqdci;              /* dwdot[NXC3H7]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[11];
    J[571] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[590] += dqdci;              /* dwdot[NXC3H7]/d[C2H4] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 146: NXC3H7 + HO2 => C3H8 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[30];
    k_f = 1.0000000000000002e-06 * 2080000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[30] -= q; /* NXC3H7 */
    wdot[31] += q; /* C3H8 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[30];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[310] -= dqdci;              /* dwdot[NXC3H7]/d[HO2] */
    J[311] += dqdci;              /* dwdot[C3H8]/d[HO2] */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f*sc[8];
    J[1055] += dqdci;             /* dwdot[O2]/d[NXC3H7] */
    J[1058] -= dqdci;             /* dwdot[HO2]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    J[1081] += dqdci;             /* dwdot[C3H8]/d[NXC3H7] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] += dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 147: NXC3H7 + O2 => C3H6 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[30];
    k_f = 1.0000000000000002e-06 * 300000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3000  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[27] += q; /* C3H6 */
    wdot[30] -= q; /* NXC3H7 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[30];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[202] += dqdci;              /* dwdot[C3H6]/d[O2] */
    J[205] -= dqdci;              /* dwdot[NXC3H7]/d[O2] */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f*sc[5];
    J[1055] -= dqdci;             /* dwdot[O2]/d[NXC3H7] */
    J[1058] += dqdci;             /* dwdot[HO2]/d[NXC3H7] */
    J[1077] += dqdci;             /* dwdot[C3H6]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 148: NXC3H7 => H + C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[30];
    k_f = 1 * 2667000000000000
                * exp(-0.64000000000000001 * tc[0] - 0.50321666580471969 * (36820.029999999999) * invT);
    dlnkfdT = -0.64000000000000001 * invT + 0.50321666580471969 *  36820.029999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[27] += q; /* C3H6 */
    wdot[30] -= q; /* NXC3H7 */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f;
    J[1052] += dqdci;             /* dwdot[H]/d[NXC3H7] */
    J[1077] += dqdci;             /* dwdot[C3H6]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 149: H + C3H6 => NXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = 1.0000000000000002e-06 * 10000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (2500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  2500  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[27] -= q; /* C3H6 */
    wdot[30] += q; /* NXC3H7 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[97] -= dqdci;               /* dwdot[C3H6]/d[H] */
    J[100] += dqdci;              /* dwdot[NXC3H7]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[2];
    J[947] -= dqdci;              /* dwdot[H]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    J[975] += dqdci;              /* dwdot[NXC3H7]/d[C3H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 150: C3H8 + OH => NXC3H7 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[31];
    k_f = 1.0000000000000002e-06 * 10540000000
                * exp(0.96999999999999997 * tc[0] - 0.50321666580471969 * (1586.04) * invT);
    dlnkfdT = 0.96999999999999997 * invT + 0.50321666580471969 *  1586.04  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[31];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[135] += dqdci;              /* dwdot[NXC3H7]/d[OH] */
    J[136] -= dqdci;              /* dwdot[C3H8]/d[OH] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[3];
    J[1088] -= dqdci;             /* dwdot[OH]/d[C3H8] */
    J[1091] += dqdci;             /* dwdot[H2O]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 151: C3H8 + HO2 => NXC3H7 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[31];
    k_f = 1.0000000000000002e-06 * 16800000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (20429.970000000001) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  20429.970000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[31];
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[310] += dqdci;              /* dwdot[NXC3H7]/d[HO2] */
    J[311] -= dqdci;              /* dwdot[C3H8]/d[HO2] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[8];
    J[1092] += dqdci;             /* dwdot[H2O2]/d[C3H8] */
    J[1093] -= dqdci;             /* dwdot[HO2]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 152: H + C3H8 <=> H2 + NXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[31];
    k_f = 1.0000000000000002e-06 * 3972000
                * exp(2.75 * tc[0] - 0.50321666580471969 * (6756.6899999999996) * invT);
    dlnkfdT = 2.75 * invT + 0.50321666580471969 *  6756.6899999999996  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[30];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[30] + g_RT[31]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[31]) + (h_RT[4] + h_RT[30]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[31];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[100] += dqdci;              /* dwdot[NXC3H7]/d[H] */
    J[101] -= dqdci;              /* dwdot[C3H8]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[30];
    J[142] -= dqdci;              /* dwdot[H]/d[H2] */
    J[144] += dqdci;              /* dwdot[H2]/d[H2] */
    J[170] += dqdci;              /* dwdot[NXC3H7]/d[H2] */
    J[171] -= dqdci;              /* dwdot[C3H8]/d[H2] */
    /* d()/d[NXC3H7] */
    dqdci =  - k_r*sc[4];
    J[1052] -= dqdci;             /* dwdot[H]/d[NXC3H7] */
    J[1054] += dqdci;             /* dwdot[H2]/d[NXC3H7] */
    J[1080] += dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    J[1081] -= dqdci;             /* dwdot[C3H8]/d[NXC3H7] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[2];
    J[1087] -= dqdci;             /* dwdot[H]/d[C3H8] */
    J[1089] += dqdci;             /* dwdot[H2]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 153: C3H8 + OH => IXC3H7 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[31];
    k_f = 1.0000000000000002e-06 * 46700000
                * exp(1.6100000000000001 * tc[0] - 0.50321666580471969 * (-34.890000000000001) * invT);
    dlnkfdT = 1.6100000000000001 * invT + 0.50321666580471969 *  -34.890000000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[31];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[134] += dqdci;              /* dwdot[IXC3H7]/d[OH] */
    J[136] -= dqdci;              /* dwdot[C3H8]/d[OH] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[3];
    J[1088] -= dqdci;             /* dwdot[OH]/d[C3H8] */
    J[1091] += dqdci;             /* dwdot[H2O]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 154: CH3 + C3H8 => CH4 + IXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[31];
    k_f = 1.0000000000000002e-06 * 398000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (9500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  9500  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[31];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[414] += dqdci;              /* dwdot[IXC3H7]/d[CH3] */
    J[416] -= dqdci;              /* dwdot[C3H8]/d[CH3] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[11];
    J[1096] -= dqdci;             /* dwdot[CH3]/d[C3H8] */
    J[1100] += dqdci;             /* dwdot[CH4]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 155: CH3 + C3H8 => CH4 + NXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[31];
    k_f = 1.0000000000000002e-06 * 1290000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (11599.9) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  11599.9  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[31];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[415] += dqdci;              /* dwdot[NXC3H7]/d[CH3] */
    J[416] -= dqdci;              /* dwdot[C3H8]/d[CH3] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[11];
    J[1096] -= dqdci;             /* dwdot[CH3]/d[C3H8] */
    J[1100] += dqdci;             /* dwdot[CH4]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 156: C3H8 + O => IXC3H7 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[31];
    k_f = 1.0000000000000002e-06 * 28100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (5200.0500000000002) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  5200.0500000000002  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[31];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[64] += dqdci;               /* dwdot[IXC3H7]/d[O] */
    J[66] -= dqdci;               /* dwdot[C3H8]/d[O] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[1];
    J[1086] -= dqdci;             /* dwdot[O]/d[C3H8] */
    J[1088] += dqdci;             /* dwdot[OH]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 157: C3H8 + HO2 => IXC3H7 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[31];
    k_f = 1.0000000000000002e-06 * 5600000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (17700.049999999999) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  17700.049999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[31];
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[309] += dqdci;              /* dwdot[IXC3H7]/d[HO2] */
    J[311] -= dqdci;              /* dwdot[C3H8]/d[HO2] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[8];
    J[1092] += dqdci;             /* dwdot[H2O2]/d[C3H8] */
    J[1093] -= dqdci;             /* dwdot[HO2]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 158: C3H8 + O => NXC3H7 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[31];
    k_f = 1.0000000000000002e-06 * 113000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7849.8999999999996) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  7849.8999999999996  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[31];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[65] += dqdci;               /* dwdot[NXC3H7]/d[O] */
    J[66] -= dqdci;               /* dwdot[C3H8]/d[O] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[1];
    J[1086] -= dqdci;             /* dwdot[O]/d[C3H8] */
    J[1088] += dqdci;             /* dwdot[OH]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 159: C3H8 + O2 => IXC3H7 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[31];
    k_f = 1.0000000000000002e-06 * 40000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (47500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  47500  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[31];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[204] += dqdci;              /* dwdot[IXC3H7]/d[O2] */
    J[206] -= dqdci;              /* dwdot[C3H8]/d[O2] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[5];
    J[1090] -= dqdci;             /* dwdot[O2]/d[C3H8] */
    J[1093] += dqdci;             /* dwdot[HO2]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 160: IXC3H7 + HO2 => C3H8 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[29];
    k_f = 1.0000000000000002e-06 * 2080000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[29] -= q; /* IXC3H7 */
    wdot[31] += q; /* C3H8 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[29];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[309] -= dqdci;              /* dwdot[IXC3H7]/d[HO2] */
    J[311] += dqdci;              /* dwdot[C3H8]/d[HO2] */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f*sc[8];
    J[1020] += dqdci;             /* dwdot[O2]/d[IXC3H7] */
    J[1023] -= dqdci;             /* dwdot[HO2]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    J[1046] += dqdci;             /* dwdot[C3H8]/d[IXC3H7] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] += dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 161: H + C3H8 => H2 + IXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[31];
    k_f = 1.0000000000000002e-06 * 1300000
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (4471.0799999999999) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  4471.0799999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[31];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[99] += dqdci;               /* dwdot[IXC3H7]/d[H] */
    J[101] -= dqdci;              /* dwdot[C3H8]/d[H] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[2];
    J[1087] -= dqdci;             /* dwdot[H]/d[C3H8] */
    J[1089] += dqdci;             /* dwdot[H2]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 162: H2 + IXC3H7 => H + C3H8 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[29];
    k_f = 1.0000000000000002e-06 * 470900
                * exp(2.1499999999999999 * tc[0] - 0.50321666580471969 * (12179.969999999999) * invT);
    dlnkfdT = 2.1499999999999999 * invT + 0.50321666580471969 *  12179.969999999999  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[29] -= q; /* IXC3H7 */
    wdot[31] += q; /* C3H8 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[29];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[169] -= dqdci;              /* dwdot[IXC3H7]/d[H2] */
    J[171] += dqdci;              /* dwdot[C3H8]/d[H2] */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f*sc[4];
    J[1017] += dqdci;             /* dwdot[H]/d[IXC3H7] */
    J[1019] -= dqdci;             /* dwdot[H2]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    J[1046] += dqdci;             /* dwdot[C3H8]/d[IXC3H7] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] += dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 163: C3H5O => C2H3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[28];
    k_f = 1 * 2028000000000
                * exp(0.089999999999999997 * tc[0] - 0.50321666580471969 * (23559.990000000002) * invT);
    dlnkfdT = 0.089999999999999997 * invT + 0.50321666580471969 *  23559.990000000002  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] += q; /* CH2O */
    wdot[24] += q; /* C2H3 */
    wdot[28] -= q; /* C3H5O */
    /* d()/d[C3H5O] */
    dqdci =  + k_f;
    J[990] += dqdci;              /* dwdot[CH2O]/d[C3H5O] */
    J[1004] += dqdci;             /* dwdot[C2H3]/d[C3H5O] */
    J[1008] -= dqdci;             /* dwdot[C3H5O]/d[C3H5O] */
    /* d()/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */
    J[1218] -= dqdT;              /* dwdot[C3H5O]/dT */

    /*reaction 164: IXC3H7O2 => IXC3H7 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[32];
    k_f = 1 * 2.803e+17
                * exp(-0.62 * tc[0] - 0.50321666580471969 * (36039.910000000003) * invT);
    dlnkfdT = -0.62 * invT + 0.50321666580471969 *  36039.910000000003  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[29] += q; /* IXC3H7 */
    wdot[32] -= q; /* IXC3H7O2 */
    /* d()/d[IXC3H7O2] */
    dqdci =  + k_f;
    J[1125] += dqdci;             /* dwdot[O2]/d[IXC3H7O2] */
    J[1149] += dqdci;             /* dwdot[IXC3H7]/d[IXC3H7O2] */
    J[1152] -= dqdci;             /* dwdot[IXC3H7O2]/d[IXC3H7O2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1222] -= dqdT;              /* dwdot[IXC3H7O2]/dT */

    /*reaction 165: IXC3H7 + O2 => IXC3H7O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[29];
    k_f = 1.0000000000000002e-06 * 7540000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[29] -= q; /* IXC3H7 */
    wdot[32] += q; /* IXC3H7O2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[29];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[204] -= dqdci;              /* dwdot[IXC3H7]/d[O2] */
    J[207] += dqdci;              /* dwdot[IXC3H7O2]/d[O2] */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f*sc[5];
    J[1020] -= dqdci;             /* dwdot[O2]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    J[1047] += dqdci;             /* dwdot[IXC3H7O2]/d[IXC3H7] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */
    J[1222] += dqdT;              /* dwdot[IXC3H7O2]/dT */

    /*reaction 166: NXC3H7O2 => NXC3H7 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[33];
    k_f = 1 * 3.364e+19
                * exp(-1.3200000000000001 * tc[0] - 0.50321666580471969 * (35760.040000000001) * invT);
    dlnkfdT = -1.3200000000000001 * invT + 0.50321666580471969 *  35760.040000000001  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[30] += q; /* NXC3H7 */
    wdot[33] -= q; /* NXC3H7O2 */
    /* d()/d[NXC3H7O2] */
    dqdci =  + k_f;
    J[1160] += dqdci;             /* dwdot[O2]/d[NXC3H7O2] */
    J[1185] += dqdci;             /* dwdot[NXC3H7]/d[NXC3H7O2] */
    J[1188] -= dqdci;             /* dwdot[NXC3H7O2]/d[NXC3H7O2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1223] -= dqdT;              /* dwdot[NXC3H7O2]/dT */

    /*reaction 167: NXC3H7 + O2 => NXC3H7O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[30];
    k_f = 1.0000000000000002e-06 * 4520000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[30] -= q; /* NXC3H7 */
    wdot[33] += q; /* NXC3H7O2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[30];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[205] -= dqdci;              /* dwdot[NXC3H7]/d[O2] */
    J[208] += dqdci;              /* dwdot[NXC3H7O2]/d[O2] */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f*sc[5];
    J[1055] -= dqdci;             /* dwdot[O2]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    J[1083] += dqdci;             /* dwdot[NXC3H7O2]/d[NXC3H7] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */
    J[1223] += dqdT;              /* dwdot[NXC3H7O2]/dT */

    amrex::Real c_R[34], dcRdT[34], e_RT[34];
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
    for (int k = 0; k < 34; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[1190+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 34; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 34; ++m) {
            dehmixdc += eh_RT[m]*J[k*35+m];
        }
        J[k*35+34] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[1224] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;

return;
}
#endif


/*compute the reaction Jacobian */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * consP)
{
    amrex::Real c[34];

    for (int k=0; k<34; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<34; k++) {
        J[1190+k] *= 1.e-6;
        J[k*35+34] *= 1.e6;
    }

    return;
}

/* gauss-jordan solver external routine */
/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void sgjsolve(amrex::Real* A, amrex::Real* x, amrex::Real* b) {
    amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void sgjsolve_simplified(amrex::Real* A, amrex::Real* x, amrex::Real* b) {
    amrex::Abort("sgjsolve_simplified not implemented, choose a different solver ");
}

#endif
