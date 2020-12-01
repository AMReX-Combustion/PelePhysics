#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <AMReX_Gpu.H>

extern "C"
{
AMREX_GPU_HOST_DEVICE void get_imw(amrex::Real imw_new[]);
AMREX_GPU_HOST_DEVICE void get_mw(amrex::Real mw_new[]);
void egtransetEPS(amrex::Real *  EPS);
void egtransetSIG(amrex::Real* SIG);
void atomicWeight(amrex::Real *  awt);
void molecularWeight(amrex::Real *  wt);
AMREX_GPU_HOST_DEVICE void gibbs(amrex::Real *  species, amrex::Real *  tc);
AMREX_GPU_HOST_DEVICE void helmholtz(amrex::Real *  species, amrex::Real *  tc);
AMREX_GPU_HOST_DEVICE void speciesInternalEnergy(amrex::Real *  species, amrex::Real *  tc);
AMREX_GPU_HOST_DEVICE void speciesEnthalpy(amrex::Real *  species, amrex::Real *  tc);
AMREX_GPU_HOST_DEVICE void speciesEntropy(amrex::Real *  species, amrex::Real *  tc);
AMREX_GPU_HOST_DEVICE void cp_R(amrex::Real *  species, amrex::Real *  tc);
AMREX_GPU_HOST_DEVICE void cv_R(amrex::Real *  species, amrex::Real *  tc);
void equilibriumConstants(amrex::Real *  kc, amrex::Real *  g_RT, amrex::Real T);
AMREX_GPU_HOST_DEVICE void productionRate(amrex::Real *  wdot, amrex::Real *  sc, amrex::Real T);
AMREX_GPU_HOST_DEVICE void comp_qfqr(amrex::Real *  q_f, amrex::Real *  q_r, amrex::Real *  sc, amrex::Real *  tc, amrex::Real invT);
#ifndef AMREX_USE_CUDA
void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f);
void comp_Kc(amrex::Real *  tc, amrex::Real invT, amrex::Real *  Kc);
#endif
AMREX_GPU_HOST_DEVICE void progressRate(amrex::Real *  qdot, amrex::Real *  speciesConc, amrex::Real T);
AMREX_GPU_HOST_DEVICE void progressRateFR(amrex::Real *  q_f, amrex::Real *  q_r, amrex::Real *  speciesConc, amrex::Real T);
AMREX_GPU_HOST_DEVICE void CKINIT();
AMREX_GPU_HOST_DEVICE void CKFINALIZE();
#ifndef AMREX_USE_CUDA
void GET_REACTION_MAP(int *  rmap);
void SetAllDefaults();
#endif
void CKINDX(int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, amrex::Real *  rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, amrex::Real *  rval, int * kerr, int lenline, int lenkray);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYME(int * kname, int * lenkname);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(amrex::Real *  ru, amrex::Real *  ruc, amrex::Real *  pa);
void CKPX(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  P);
AMREX_GPU_HOST_DEVICE void CKPY(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  P);
void CKPC(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  c, amrex::Real *  P);
void CKRHOX(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  rho);
AMREX_GPU_HOST_DEVICE void CKRHOY(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  rho);
void CKRHOC(amrex::Real *  P, amrex::Real *  T, amrex::Real *  c, amrex::Real *  rho);
void CKWT(amrex::Real *  wt);
void CKAWT(amrex::Real *  awt);
AMREX_GPU_HOST_DEVICE void CKMMWY(amrex::Real *  y, amrex::Real *  wtm);
void CKMMWX(amrex::Real *  x, amrex::Real *  wtm);
void CKMMWC(amrex::Real *  c, amrex::Real *  wtm);
AMREX_GPU_HOST_DEVICE void CKYTX(amrex::Real *  y, amrex::Real *  x);
void CKYTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  c);
AMREX_GPU_HOST_DEVICE void CKYTCR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  c);
AMREX_GPU_HOST_DEVICE void CKXTY(amrex::Real *  x, amrex::Real *  y);
void CKXTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  c);
void CKXTCR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  c);
void CKCTX(amrex::Real *  c, amrex::Real *  x);
void CKCTY(amrex::Real *  c, amrex::Real *  y);
void CKCPOR(amrex::Real *  T, amrex::Real *  cpor);
void CKHORT(amrex::Real *  T, amrex::Real *  hort);
void CKSOR(amrex::Real *  T, amrex::Real *  sor);
void CKCVML(amrex::Real *  T, amrex::Real *  cvml);
void CKCPML(amrex::Real *  T, amrex::Real *  cvml);
void CKUML(amrex::Real *  T, amrex::Real *  uml);
void CKHML(amrex::Real *  T, amrex::Real *  uml);
void CKGML(amrex::Real *  T, amrex::Real *  gml);
void CKAML(amrex::Real *  T, amrex::Real *  aml);
void CKSML(amrex::Real *  T, amrex::Real *  sml);
AMREX_GPU_HOST_DEVICE void CKCVMS(amrex::Real *  T, amrex::Real *  cvms);
AMREX_GPU_HOST_DEVICE void CKCPMS(amrex::Real *  T, amrex::Real *  cvms);
AMREX_GPU_HOST_DEVICE void CKUMS(amrex::Real *  T, amrex::Real *  ums);
AMREX_GPU_HOST_DEVICE void CKHMS(amrex::Real *  T, amrex::Real *  ums);
void CKGMS(amrex::Real *  T, amrex::Real *  gms);
void CKAMS(amrex::Real *  T, amrex::Real *  ams);
void CKSMS(amrex::Real *  T, amrex::Real *  sms);
void CKCPBL(amrex::Real *  T, amrex::Real *  x, amrex::Real *  cpbl);
AMREX_GPU_HOST_DEVICE void CKCPBS(amrex::Real *  T, amrex::Real *  y, amrex::Real *  cpbs);
void CKCVBL(amrex::Real *  T, amrex::Real *  x, amrex::Real *  cpbl);
AMREX_GPU_HOST_DEVICE void CKCVBS(amrex::Real *  T, amrex::Real *  y, amrex::Real *  cpbs);
void CKHBML(amrex::Real *  T, amrex::Real *  x, amrex::Real *  hbml);
AMREX_GPU_HOST_DEVICE void CKHBMS(amrex::Real *  T, amrex::Real *  y, amrex::Real *  hbms);
void CKUBML(amrex::Real *  T, amrex::Real *  x, amrex::Real *  ubml);
AMREX_GPU_HOST_DEVICE void CKUBMS(amrex::Real *  T, amrex::Real *  y, amrex::Real *  ubms);
void CKSBML(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  sbml);
void CKSBMS(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  sbms);
void CKGBML(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  gbml);
void CKGBMS(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  gbms);
void CKABML(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  abml);
void CKABMS(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  abms);
AMREX_GPU_HOST_DEVICE void CKWC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  wdot);
void CKWYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  wdot);
void CKWXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  wdot);
AMREX_GPU_HOST_DEVICE void CKWYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  wdot);
void CKWXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  wdot);
void CKQC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  qdot);
void CKKFKR(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  q_f, amrex::Real *  q_r);
void CKQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot);
void CKQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot);
void CKQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot);
void CKQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot);
void CKNU(int * kdim, int * nuki);
#ifndef AMREX_USE_CUDA
void CKINU(int * i, int * nspec, int * ki, int * nu);
#endif
void CKNCF(int * ncf);
void CKABE(amrex::Real *  a, amrex::Real *  b, amrex::Real *  e );
void CKEQC(amrex::Real *  T, amrex::Real *  C , amrex::Real *  eqcon );
void CKEQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  eqcon);
void CKEQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  eqcon);
void CKEQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  eqcon);
void CKEQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  eqcon);
AMREX_GPU_HOST_DEVICE void DWDOT(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  T, int * consP);
AMREX_GPU_HOST_DEVICE void DWDOT_SIMPLIFIED(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * HP);
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO(int * nJdata, int * consP, int NCELLS);
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO_SYST(int * nJdata, int * consP, int NCELLS);
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO_SYST_SIMPLIFIED(int * nJdata, int * consP);
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_CSC(int * rowVals, int * colPtrs, int * consP, int NCELLS);
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP);
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base);
AMREX_GPU_HOST_DEVICE void aJacobian(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP);
AMREX_GPU_HOST_DEVICE void aJacobian_precond(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int HP);
AMREX_GPU_HOST_DEVICE void dcvpRdT(amrex::Real *  species, amrex::Real *  tc);
AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_EY(amrex::Real *  e, amrex::Real *  y, amrex::Real *  t, int *ierr);
AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_HY(amrex::Real *  h, amrex::Real *  y, amrex::Real *  t, int *ierr);
void GET_CRITPARAMS(amrex::Real *  Tci, amrex::Real *  ai, amrex::Real *  bi, amrex::Real *  acentric_i);
/*vector version */
void VCKYTX(int *  np, amrex::Real *  y, amrex::Real *  x);
void VCKHMS(int *  np, amrex::Real *  T, amrex::Real *  ums);
void VCKWYR(int *  np, amrex::Real *  rho, amrex::Real *  T,
            amrex::Real *  y,
            amrex::Real *  wdot);
#ifndef AMREX_USE_CUDA
void vproductionRate(int npt, amrex::Real *  wdot, amrex::Real *  c, amrex::Real *  T);
void VCKPY(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  P);
void vcomp_k_f(int npt, amrex::Real *  k_f_s, amrex::Real *  tc, amrex::Real *  invT);
void vcomp_gibbs(int npt, amrex::Real *  g_RT, amrex::Real *  tc);
void vcomp_Kc(int npt, amrex::Real *  Kc_s, amrex::Real *  g_RT, amrex::Real *  invT);
void vcomp_wdot_1_50(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
                amrex::Real *  k_f_s, amrex::Real *  Kc_s,
                amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T);
void vcomp_wdot_51_84(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
                amrex::Real *  k_f_s, amrex::Real *  Kc_s,
                amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T);
#endif
/*Transport function declarations */
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
/*gauss-jordan solver external routine */
AMREX_GPU_HOST_DEVICE void sgjsolve(amrex::Real* A, amrex::Real* x, amrex::Real* b);
AMREX_GPU_HOST_DEVICE void sgjsolve_simplified(amrex::Real* A, amrex::Real* x, amrex::Real* b);
}

#ifndef AMREX_USE_CUDA
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

    extern amrex::Real fwd_A_DEF[84], fwd_beta_DEF[84], fwd_Ea_DEF[84];
    extern amrex::Real low_A_DEF[84], low_beta_DEF[84], low_Ea_DEF[84];
    extern amrex::Real rev_A_DEF[84], rev_beta_DEF[84], rev_Ea_DEF[84];
    extern amrex::Real troe_a_DEF[84],troe_Ts_DEF[84], troe_Tss_DEF[84], troe_Tsss_DEF[84];
    extern amrex::Real sri_a_DEF[84], sri_b_DEF[84], sri_c_DEF[84], sri_d_DEF[84], sri_e_DEF[84];
    extern amrex::Real activation_units_DEF[84], prefactor_units_DEF[84], phase_units_DEF[84];
    extern int is_PD_DEF[84], troe_len_DEF[84], sri_len_DEF[84], nTB_DEF[84], *TBid_DEF[84];
    extern amrex::Real *TB_DEF[84];
    extern std::vector<int> rxn_map;
}
#endif
