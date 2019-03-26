
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if defined(BL_FORT_USE_UPPERCASE)
#define CKINDX CKINDX
#define CKINIT CKINIT
#define CKFINALIZE CKFINALIZE
#define CKXNUM CKXNUM
#define CKSYME CKSYME
#define CKSYMS CKSYMS
#define CKRP CKRP
#define CKPX CKPX
#define CKPY CKPY
#define CKPC CKPC
#define CKRHOX CKRHOX
#define CKRHOY CKRHOY
#define CKRHOC CKRHOC
#define CKWT CKWT
#define CKAWT CKAWT
#define CKMMWY CKMMWY
#define CKMMWX CKMMWX
#define CKMMWC CKMMWC
#define CKYTX CKYTX
#define CKYTCP CKYTCP
#define CKYTCR CKYTCR
#define CKXTY CKXTY
#define CKXTCP CKXTCP
#define CKXTCR CKXTCR
#define CKCTX CKCTX
#define CKCTY CKCTY
#define CKCPOR CKCPOR
#define CKHORT CKHORT
#define CKSOR CKSOR
#define CKCVML CKCVML
#define CKCPML CKCPML
#define CKUML CKUML
#define CKHML CKHML
#define CKGML CKGML
#define CKAML CKAML
#define CKSML CKSML
#define CKCVMS CKCVMS
#define CKCPMS CKCPMS
#define CKUMS CKUMS
#define CKHMS CKHMS
#define CKGMS CKGMS
#define CKAMS CKAMS
#define CKSMS CKSMS
#define CKCPBL CKCPBL
#define CKCPBS CKCPBS
#define CKCVBL CKCVBL
#define CKCVBS CKCVBS
#define CKHBML CKHBML
#define CKHBMS CKHBMS
#define CKUBML CKUBML
#define CKUBMS CKUBMS
#define CKSBML CKSBML
#define CKSBMS CKSBMS
#define CKGBML CKGBML
#define CKGBMS CKGBMS
#define CKABML CKABML
#define CKABMS CKABMS
#define CKWC CKWC
#define CKWYP CKWYP
#define CKWXP CKWXP
#define CKWYR CKWYR
#define CKWXR CKWXR
#define CKQC CKQC
#define CKKFKR CKKFKR
#define CKQYP CKQYP
#define CKQXP CKQXP
#define CKQYR CKQYR
#define CKQXR CKQXR
#define CKNU CKNU
#define CKNCF CKNCF
#define CKABE CKABE
#define CKEQC CKEQC
#define CKEQYP CKEQYP
#define CKEQXP CKEQXP
#define CKEQYR CKEQYR
#define CKEQXR CKEQXR
#define DWDOT DWDOT
#define DWDOT_PRECOND DWDOT_PRECOND
#define SPARSITY_INFO SPARSITY_INFO
#define SPARSITY_INFO_PRECOND SPARSITY_INFO_PRECOND
#define SPARSITY_PREPROC SPARSITY_PREPROC
#define SPARSITY_PREPROC_PRECOND SPARSITY_PREPROC_PRECOND
#define VCKHMS VCKHMS
#define VCKPY VCKPY
#define VCKWYR VCKWYR
#define VCKYTX VCKYTX
#define GET_T_GIVEN_EY GET_T_GIVEN_EY
#define GET_T_GIVEN_HY GET_T_GIVEN_HY
#define GET_REACTION_MAP GET_REACTION_MAP
#define GET_CRITPARAMS GET_CRITPARAMS
#elif defined(BL_FORT_USE_LOWERCASE)
#define CKINDX ckindx
#define CKINIT ckinit
#define CKFINALIZE ckfinalize
#define CKXNUM ckxnum
#define CKSYME cksyme
#define CKSYMS cksyms
#define CKRP ckrp
#define CKPX ckpx
#define CKPY ckpy
#define CKPC ckpc
#define CKRHOX ckrhox
#define CKRHOY ckrhoy
#define CKRHOC ckrhoc
#define CKWT ckwt
#define CKAWT ckawt
#define CKMMWY ckmmwy
#define CKMMWX ckmmwx
#define CKMMWC ckmmwc
#define CKYTX ckytx
#define CKYTCP ckytcp
#define CKYTCR ckytcr
#define CKXTY ckxty
#define CKXTCP ckxtcp
#define CKXTCR ckxtcr
#define CKCTX ckctx
#define CKCTY ckcty
#define CKCPOR ckcpor
#define CKHORT ckhort
#define CKSOR cksor
#define CKCVML ckcvml
#define CKCPML ckcpml
#define CKUML ckuml
#define CKHML ckhml
#define CKGML ckgml
#define CKAML ckaml
#define CKSML cksml
#define CKCVMS ckcvms
#define CKCPMS ckcpms
#define CKUMS ckums
#define CKHMS ckhms
#define CKGMS ckgms
#define CKAMS ckams
#define CKSMS cksms
#define CKCPBL ckcpbl
#define CKCPBS ckcpbs
#define CKCVBL ckcvbl
#define CKCVBS ckcvbs
#define CKHBML ckhbml
#define CKHBMS ckhbms
#define CKUBML ckubml
#define CKUBMS ckubms
#define CKSBML cksbml
#define CKSBMS cksbms
#define CKGBML ckgbml
#define CKGBMS ckgbms
#define CKABML ckabml
#define CKABMS ckabms
#define CKWC ckwc
#define CKWYP ckwyp
#define CKWXP ckwxp
#define CKWYR ckwyr
#define CKWXR ckwxr
#define CKQC ckqc
#define CKKFKR ckkfkr
#define CKQYP ckqyp
#define CKQXP ckqxp
#define CKQYR ckqyr
#define CKQXR ckqxr
#define CKNU cknu
#define CKNCF ckncf
#define CKABE ckabe
#define CKEQC ckeqc
#define CKEQYP ckeqyp
#define CKEQXP ckeqxp
#define CKEQYR ckeqyr
#define CKEQXR ckeqxr
#define DWDOT dwdot
#define DWDOT_PRECOND dwdot_precond
#define SPARSITY_INFO sparsity_info
#define SPARSITY_INFO_PRECOND sparsity_info_precond
#define SPARSITY_PREPROC sparsity_preproc
#define SPARSITY_PREPROC_PRECOND sparsity_preproc_precond
#define VCKHMS vckhms
#define VCKPY vckpy
#define VCKWYR vckwyr
#define VCKYTX vckytx
#define GET_T_GIVEN_EY get_t_given_ey
#define GET_T_GIVEN_HY get_t_given_hy
#define GET_REACTION_MAP get_reaction_map
#define GET_CRITPARAMS get_critparams
#elif defined(BL_FORT_USE_UNDERSCORE)
#define CKINDX ckindx_
#define CKINIT ckinit_
#define CKFINALIZE ckfinalize_
#define CKXNUM ckxnum_
#define CKSYME cksyme_
#define CKSYMS cksyms_
#define CKRP ckrp_
#define CKPX ckpx_
#define CKPY ckpy_
#define CKPC ckpc_
#define CKRHOX ckrhox_
#define CKRHOY ckrhoy_
#define CKRHOC ckrhoc_
#define CKWT ckwt_
#define CKAWT ckawt_
#define CKMMWY ckmmwy_
#define CKMMWX ckmmwx_
#define CKMMWC ckmmwc_
#define CKYTX ckytx_
#define CKYTCP ckytcp_
#define CKYTCR ckytcr_
#define CKXTY ckxty_
#define CKXTCP ckxtcp_
#define CKXTCR ckxtcr_
#define CKCTX ckctx_
#define CKCTY ckcty_
#define CKCPOR ckcpor_
#define CKHORT ckhort_
#define CKSOR cksor_
#define CKCVML ckcvml_
#define CKCPML ckcpml_
#define CKUML ckuml_
#define CKHML ckhml_
#define CKGML ckgml_
#define CKAML ckaml_
#define CKSML cksml_
#define CKCVMS ckcvms_
#define CKCPMS ckcpms_
#define CKUMS ckums_
#define CKHMS ckhms_
#define CKGMS ckgms_
#define CKAMS ckams_
#define CKSMS cksms_
#define CKCPBL ckcpbl_
#define CKCPBS ckcpbs_
#define CKCVBL ckcvbl_
#define CKCVBS ckcvbs_
#define CKHBML ckhbml_
#define CKHBMS ckhbms_
#define CKUBML ckubml_
#define CKUBMS ckubms_
#define CKSBML cksbml_
#define CKSBMS cksbms_
#define CKGBML ckgbml_
#define CKGBMS ckgbms_
#define CKABML ckabml_
#define CKABMS ckabms_
#define CKWC ckwc_
#define CKWYP ckwyp_
#define CKWXP ckwxp_
#define CKWYR ckwyr_
#define CKWXR ckwxr_
#define CKQC ckqc_
#define CKKFKR ckkfkr_
#define CKQYP ckqyp_
#define CKQXP ckqxp_
#define CKQYR ckqyr_
#define CKQXR ckqxr_
#define CKNU cknu_
#define CKNCF ckncf_
#define CKABE ckabe_
#define CKEQC ckeqc_
#define CKEQYP ckeqyp_
#define CKEQXP ckeqxp_
#define CKEQYR ckeqyr_
#define CKEQXR ckeqxr_
#define DWDOT dwdot_
#define DWDOT_PRECOND dwdot_precond_
#define SPARSITY_INFO sparsity_info_
#define SPARSITY_INFO_PRECOND sparsity_info_precond_ 
#define SPARSITY_PREPROC sparsity_preproc_
#define SPARSITY_PREPROC_PRECOND sparsity_preproc_precond_
#define VCKHMS vckhms_
#define VCKPY vckpy_
#define VCKWYR vckwyr_
#define VCKYTX vckytx_
#define GET_T_GIVEN_EY get_t_given_ey_
#define GET_T_GIVEN_HY get_t_given_hy_
#define GET_REACTION_MAP get_reaction_map_
#define GET_CRITPARAMS get_critparams_
#endif

/*function declarations */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
extern "C"
{
void egtransetEPS(double *  EPS);
void egtransetSIG(double* SIG);
void atomicWeight(double *  awt);
void molecularWeight(double *  wt);
void gibbs(double *  species, double *  tc);
void helmholtz(double *  species, double *  tc);
void speciesInternalEnergy(double *  species, double *  tc);
void speciesEnthalpy(double *  species, double *  tc);
void speciesEntropy(double *  species, double *  tc);
void cp_R(double *  species, double *  tc);
void cv_R(double *  species, double *  tc);
void equilibriumConstants(double *  kc, double *  g_RT, double T);
void productionRate(double *  wdot, double *  sc, double T);
void comp_k_f(double *  tc, double invT, double *  k_f);
void comp_Kc(double *  tc, double invT, double *  Kc);
void comp_qfqr(double *  q_f, double *  q_r, double *  sc, double *  tc, double invT);
void progressRate(double *  qdot, double *  speciesConc, double T);
void progressRateFR(double *  q_f, double *  q_r, double *  speciesConc, double T);
void CKINIT();
void CKFINALIZE();
void CKINDX(int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double *  rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double *  rval, int * kerr, int lenline, int lenkray);
void CKSYME(int * kname, int * lenkname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(double *  ru, double *  ruc, double *  pa);
void CKPX(double *  rho, double *  T, double *  x, double *  P);
void CKPY(double *  rho, double *  T, double *  y, double *  P);
void CKPC(double *  rho, double *  T, double *  c, double *  P);
void CKRHOX(double *  P, double *  T, double *  x, double *  rho);
void CKRHOY(double *  P, double *  T, double *  y, double *  rho);
void CKRHOC(double *  P, double *  T, double *  c, double *  rho);
void CKWT(double *  wt);
void CKAWT(double *  awt);
void CKMMWY(double *  y, double *  wtm);
void CKMMWX(double *  x, double *  wtm);
void CKMMWC(double *  c, double *  wtm);
void CKYTX(double *  y, double *  x);
void CKYTCP(double *  P, double *  T, double *  y, double *  c);
void CKYTCR(double *  rho, double *  T, double *  y, double *  c);
void CKXTY(double *  x, double *  y);
void CKXTCP(double *  P, double *  T, double *  x, double *  c);
void CKXTCR(double *  rho, double *  T, double *  x, double *  c);
void CKCTX(double *  c, double *  x);
void CKCTY(double *  c, double *  y);
void CKCPOR(double *  T, double *  cpor);
void CKHORT(double *  T, double *  hort);
void CKSOR(double *  T, double *  sor);
void CKCVML(double *  T, double *  cvml);
void CKCPML(double *  T, double *  cvml);
void CKUML(double *  T, double *  uml);
void CKHML(double *  T, double *  uml);
void CKGML(double *  T, double *  gml);
void CKAML(double *  T, double *  aml);
void CKSML(double *  T, double *  sml);
void CKCVMS(double *  T, double *  cvms);
void CKCPMS(double *  T, double *  cvms);
void CKUMS(double *  T, double *  ums);
void CKHMS(double *  T, double *  ums);
void CKGMS(double *  T, double *  gms);
void CKAMS(double *  T, double *  ams);
void CKSMS(double *  T, double *  sms);
void CKCPBL(double *  T, double *  x, double *  cpbl);
void CKCPBS(double *  T, double *  y, double *  cpbs);
void CKCVBL(double *  T, double *  x, double *  cpbl);
void CKCVBS(double *  T, double *  y, double *  cpbs);
void CKHBML(double *  T, double *  x, double *  hbml);
void CKHBMS(double *  T, double *  y, double *  hbms);
void CKUBML(double *  T, double *  x, double *  ubml);
void CKUBMS(double *  T, double *  y, double *  ubms);
void CKSBML(double *  P, double *  T, double *  x, double *  sbml);
void CKSBMS(double *  P, double *  T, double *  y, double *  sbms);
void CKGBML(double *  P, double *  T, double *  x, double *  gbml);
void CKGBMS(double *  P, double *  T, double *  y, double *  gbms);
void CKABML(double *  P, double *  T, double *  x, double *  abml);
void CKABMS(double *  P, double *  T, double *  y, double *  abms);
void CKWC(double *  T, double *  C, double *  wdot);
void CKWYP(double *  P, double *  T, double *  y, double *  wdot);
void CKWXP(double *  P, double *  T, double *  x, double *  wdot);
void CKWYR(double *  rho, double *  T, double *  y, double *  wdot);
void CKWXR(double *  rho, double *  T, double *  x, double *  wdot);
void CKQC(double *  T, double *  C, double *  qdot);
void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r);
void CKQYP(double *  P, double *  T, double *  y, double *  qdot);
void CKQXP(double *  P, double *  T, double *  x, double *  qdot);
void CKQYR(double *  rho, double *  T, double *  y, double *  qdot);
void CKQXR(double *  rho, double *  T, double *  x, double *  qdot);
void CKNU(int * kdim, int * nuki);
void CKNCF(int * mdim, int * ncf);
void CKABE(double *  a, double *  b, double *  e );
void CKEQC(double *  T, double *  C , double *  eqcon );
void CKEQYP(double *  P, double *  T, double *  y, double *  eqcon);
void CKEQXP(double *  P, double *  T, double *  x, double *  eqcon);
void CKEQYR(double *  rho, double *  T, double *  y, double *  eqcon);
void CKEQXR(double *  rho, double *  T, double *  x, double *  eqcon);
void DWDOT(double *  J, double *  sc, double *  T, int * consP);
void DWDOT_PRECOND(double *  J, double *  sc, double *  Tp, int * HP);
void SPARSITY_INFO(int * nJdata, int * consP, int NCELLS);
void SPARSITY_INFO_PRECOND(int * nJdata, int * consP);
void SPARSITY_PREPROC(int * rowVals, int * colPtrs, int * consP, int NCELLS);
void SPARSITY_PREPROC_PRECOND(int * rowVals, int * colPtrs, int * consP);
void aJacobian(double *  J, double *  sc, double T, int consP);
void aJacobian_precond(double *  J, double *  sc, double T, int HP);
void dcvpRdT(double *  species, double *  tc);
void GET_T_GIVEN_EY(double *  e, double *  y, double *  t, int *ierr);
void GET_T_GIVEN_HY(double *  h, double *  y, double *  t, int *ierr);
void GET_REACTION_MAP(int *  rmap);
/*vector version */
void vproductionRate(int npt, double *  wdot, double *  c, double *  T);
void VCKHMS(int *  np, double *  T, double *  ums);
void VCKPY(int *  np, double *  rho, double *  T, double *  y, double *  P);
void VCKWYR(int *  np, double *  rho, double *  T,
            double *  y,
            double *  wdot);
void VCKYTX(int *  np, double *  y, double *  x);
void vcomp_k_f(int npt, double *  k_f_s, double *  tc, double *  invT);
void vcomp_gibbs(int npt, double *  g_RT, double *  tc);
void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT);
void GET_CRITPARAMS(double *  Tci, double *  ai, double *  bi, double *  acentric_i);
void vcomp_wdot(int npt, double *  wdot, double *  mixture, double *  sc,
                double *  k_f_s, double *  Kc_s,
                double *  tc, double *  invT, double *  T);

/* Inverse molecular weights */
static const double imw[12] = {
    1.0 / 2.015940,  /*H2 */
    1.0 / 31.998800,  /*O2 */
    1.0 / 15.999400,  /*O */
    1.0 / 17.007370,  /*OH */
    1.0 / 18.015340,  /*H2O */
    1.0 / 1.007970,  /*H */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 28.010550,  /*CO */
    1.0 / 44.009950,  /*CO2 */
    1.0 / 29.018520,  /*HCO */
    1.0 / 28.013400};  /*N2 */



static double fwd_A[29], fwd_beta[29], fwd_Ea[29];
static double low_A[29], low_beta[29], low_Ea[29];
static double rev_A[29], rev_beta[29], rev_Ea[29];
static double troe_a[29],troe_Ts[29], troe_Tss[29], troe_Tsss[29];
static double sri_a[29], sri_b[29], sri_c[29], sri_d[29], sri_e[29];
static double activation_units[29], prefactor_units[29], phase_units[29];
static int is_PD[29], troe_len[29], sri_len[29], nTB[29], *TBid[29];
static double *TB[29];

static double fwd_A_DEF[29], fwd_beta_DEF[29], fwd_Ea_DEF[29];
static double low_A_DEF[29], low_beta_DEF[29], low_Ea_DEF[29];
static double rev_A_DEF[29], rev_beta_DEF[29], rev_Ea_DEF[29];
static double troe_a_DEF[29],troe_Ts_DEF[29], troe_Tss_DEF[29], troe_Tsss_DEF[29];
static double sri_a_DEF[29], sri_b_DEF[29], sri_c_DEF[29], sri_d_DEF[29], sri_e_DEF[29];
static double activation_units_DEF[29], prefactor_units_DEF[29], phase_units_DEF[29];
static int is_PD_DEF[29], troe_len_DEF[29], sri_len_DEF[29], nTB_DEF[29], *TBid_DEF[29];
static double *TB_DEF[29];
static int rxn_map[29] = {8,9,10,11,3,4,5,6,0,12,13,14,15,16,17,1,18,19,20,21,22,2,23,24,25,7,26,27,28};

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<29; ++i) {
        rmap[i] = rxn_map[i];
    }
}


#include <ReactionData.H>
double* GetParamPtr(int                reaction_id,
                    REACTION_PARAMETER param_id,
                    int                species_id,
                    int                get_default)
{
  double* ret = 0;
  if (reaction_id<0 || reaction_id>=29) {
    printf("Bad reaction id = %d",reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id<0 || species_id>=12) {
      printf("GetParamPtr: Bad species id = %d",species_id);
      abort();
    }
    if (get_default) {
      for (int i=0; i<nTB_DEF[mrid]; ++i) {
        if (species_id == TBid_DEF[mrid][i]) {
          ret = &(TB_DEF[mrid][i]);
        }
      }
    }
    else {
      for (int i=0; i<nTB[mrid]; ++i) {
        if (species_id == TBid[mrid][i]) {
          ret = &(TB[mrid][i]);
        }
      }
    }
    if (ret == 0) {
      printf("GetParamPtr: No TB for reaction id = %d",reaction_id);
      abort();
    }
  }
  else {
    if (     param_id == FWD_A)     {ret = (get_default ? &(fwd_A_DEF[mrid]) : &(fwd_A[mrid]));}
      else if (param_id == FWD_BETA)  {ret = (get_default ? &(fwd_beta_DEF[mrid]) : &(fwd_beta[mrid]));}
      else if (param_id == FWD_EA)    {ret = (get_default ? &(fwd_Ea_DEF[mrid]) : &(fwd_Ea[mrid]));}
      else if (param_id == LOW_A)     {ret = (get_default ? &(low_A_DEF[mrid]) : &(low_A[mrid]));}
      else if (param_id == LOW_BETA)  {ret = (get_default ? &(low_beta_DEF[mrid]) : &(low_beta[mrid]));}
      else if (param_id == LOW_EA)    {ret = (get_default ? &(low_Ea_DEF[mrid]) : &(low_Ea[mrid]));}
      else if (param_id == REV_A)     {ret = (get_default ? &(rev_A_DEF[mrid]) : &(rev_A[mrid]));}
      else if (param_id == REV_BETA)  {ret = (get_default ? &(rev_beta_DEF[mrid]) : &(rev_beta[mrid]));}
      else if (param_id == REV_EA)    {ret = (get_default ? &(rev_Ea_DEF[mrid]) : &(rev_Ea[mrid]));}
      else if (param_id == TROE_A)    {ret = (get_default ? &(troe_a_DEF[mrid]) : &(troe_a[mrid]));}
      else if (param_id == TROE_TS)   {ret = (get_default ? &(troe_Ts_DEF[mrid]) : &(troe_Ts[mrid]));}
      else if (param_id == TROE_TSS)  {ret = (get_default ? &(troe_Tss_DEF[mrid]) : &(troe_Tss[mrid]));}
      else if (param_id == TROE_TSSS) {ret = (get_default ? &(troe_Tsss_DEF[mrid]) : &(troe_Tsss[mrid]));}
      else if (param_id == SRI_A)     {ret = (get_default ? &(sri_a_DEF[mrid]) : &(sri_a[mrid]));}
      else if (param_id == SRI_B)     {ret = (get_default ? &(sri_b_DEF[mrid]) : &(sri_b[mrid]));}
      else if (param_id == SRI_C)     {ret = (get_default ? &(sri_c_DEF[mrid]) : &(sri_c[mrid]));}
      else if (param_id == SRI_D)     {ret = (get_default ? &(sri_d_DEF[mrid]) : &(sri_d[mrid]));}
      else if (param_id == SRI_E)     {ret = (get_default ? &(sri_e_DEF[mrid]) : &(sri_e[mrid]));}
    else {
      printf("GetParamPtr: Unknown parameter id");
      abort();
    }
  }
  return ret;
}

void ResetAllParametersToDefault()
{
    for (int i=0; i<29; i++) {
        if (nTB[i] != 0) {
            nTB[i] = 0;
            free(TB[i]);
            free(TBid[i]);
        }

        fwd_A[i]    = fwd_A_DEF[i];
        fwd_beta[i] = fwd_beta_DEF[i];
        fwd_Ea[i]   = fwd_Ea_DEF[i];

        low_A[i]    = low_A_DEF[i];
        low_beta[i] = low_beta_DEF[i];
        low_Ea[i]   = low_Ea_DEF[i];

        rev_A[i]    = rev_A_DEF[i];
        rev_beta[i] = rev_beta_DEF[i];
        rev_Ea[i]   = rev_Ea_DEF[i];

        troe_a[i]    = troe_a_DEF[i];
        troe_Ts[i]   = troe_Ts_DEF[i];
        troe_Tss[i]  = troe_Tss_DEF[i];
        troe_Tsss[i] = troe_Tsss_DEF[i];

        sri_a[i] = sri_a_DEF[i];
        sri_b[i] = sri_b_DEF[i];
        sri_c[i] = sri_c_DEF[i];
        sri_d[i] = sri_d_DEF[i];
        sri_e[i] = sri_e_DEF[i];

        is_PD[i]    = is_PD_DEF[i];
        troe_len[i] = troe_len_DEF[i];
        sri_len[i]  = sri_len_DEF[i];

        activation_units[i] = activation_units_DEF[i];
        prefactor_units[i]  = prefactor_units_DEF[i];
        phase_units[i]      = phase_units_DEF[i];

        nTB[i]  = nTB_DEF[i];
        if (nTB[i] != 0) {
           TB[i] = (double *) malloc(sizeof(double) * nTB[i]);
           TBid[i] = (int *) malloc(sizeof(int) * nTB[i]);
           for (int j=0; j<nTB[i]; j++) {
             TB[i][j] = TB_DEF[i][j];
             TBid[i][j] = TBid_DEF[i][j];
           }
        }
    }
}

void SetAllDefaults()
{
    for (int i=0; i<29; i++) {
        if (nTB_DEF[i] != 0) {
            nTB_DEF[i] = 0;
            free(TB_DEF[i]);
            free(TBid_DEF[i]);
        }

        fwd_A_DEF[i]    = fwd_A[i];
        fwd_beta_DEF[i] = fwd_beta[i];
        fwd_Ea_DEF[i]   = fwd_Ea[i];

        low_A_DEF[i]    = low_A[i];
        low_beta_DEF[i] = low_beta[i];
        low_Ea_DEF[i]   = low_Ea[i];

        rev_A_DEF[i]    = rev_A[i];
        rev_beta_DEF[i] = rev_beta[i];
        rev_Ea_DEF[i]   = rev_Ea[i];

        troe_a_DEF[i]    = troe_a[i];
        troe_Ts_DEF[i]   = troe_Ts[i];
        troe_Tss_DEF[i]  = troe_Tss[i];
        troe_Tsss_DEF[i] = troe_Tsss[i];

        sri_a_DEF[i] = sri_a[i];
        sri_b_DEF[i] = sri_b[i];
        sri_c_DEF[i] = sri_c[i];
        sri_d_DEF[i] = sri_d[i];
        sri_e_DEF[i] = sri_e[i];

        is_PD_DEF[i]    = is_PD[i];
        troe_len_DEF[i] = troe_len[i];
        sri_len_DEF[i]  = sri_len[i];

        activation_units_DEF[i] = activation_units[i];
        prefactor_units_DEF[i]  = prefactor_units[i];
        phase_units_DEF[i]      = phase_units[i];

        nTB_DEF[i]  = nTB[i];
        if (nTB_DEF[i] != 0) {
           TB_DEF[i] = (double *) malloc(sizeof(double) * nTB_DEF[i]);
           TBid_DEF[i] = (int *) malloc(sizeof(int) * nTB_DEF[i]);
           for (int j=0; j<nTB_DEF[i]; j++) {
             TB_DEF[i][j] = TB[i][j];
             TBid_DEF[i][j] = TBid[i][j];
           }
        }
    }
}

/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<29; ++i) {
    free(TB[i]); TB[i] = 0; 
    free(TBid[i]); TBid[i] = 0;
    nTB[i] = 0;

    free(TB_DEF[i]); TB_DEF[i] = 0; 
    free(TBid_DEF[i]); TBid_DEF[i] = 0;
    nTB_DEF[i] = 0;
  }
}

/* Initializes parameter database */
void CKINIT()
{
    // (0):  H + O2 <=> O + OH
    fwd_A[8]     = 3547000000000000;
    fwd_beta[8]  = -0.40600000000000003;
    fwd_Ea[8]    = 16599;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = 1e-12;
    is_PD[8] = 0;
    nTB[8] = 0;

    // (1):  O + H2 <=> H + OH
    fwd_A[9]     = 50800;
    fwd_beta[9]  = 2.6699999999999999;
    fwd_Ea[9]    = 6290;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = 1e-12;
    is_PD[9] = 0;
    nTB[9] = 0;

    // (2):  H2 + OH <=> H2O + H
    fwd_A[10]     = 216000000;
    fwd_beta[10]  = 1.51;
    fwd_Ea[10]    = 3430;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = 1e-12;
    is_PD[10] = 0;
    nTB[10] = 0;

    // (3):  O + H2O <=> OH + OH
    fwd_A[11]     = 2970000;
    fwd_beta[11]  = 2.02;
    fwd_Ea[11]    = 13400;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = 1e-12;
    is_PD[11] = 0;
    nTB[11] = 0;

    // (4):  H2 + M <=> H + H + M
    fwd_A[3]     = 4.577e+19;
    fwd_beta[3]  = -1.3999999999999999;
    fwd_Ea[3]    = 104380;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = 1e-6;
    is_PD[3] = 0;
    nTB[3] = 4;
    TB[3] = (double *) malloc(4 * sizeof(double));
    TBid[3] = (int *) malloc(4 * sizeof(int));
    TBid[3][0] = 0; TB[3][0] = 2.5; // H2
    TBid[3][1] = 4; TB[3][1] = 12; // H2O
    TBid[3][2] = 8; TB[3][2] = 1.8999999999999999; // CO
    TBid[3][3] = 9; TB[3][3] = 3.7999999999999998; // CO2

    // (5):  O + O + M <=> O2 + M
    fwd_A[4]     = 6165000000000000;
    fwd_beta[4]  = -0.5;
    fwd_Ea[4]    = 0;
    prefactor_units[4]  = 1.0000000000000002e-12;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = 1e-12;
    is_PD[4] = 0;
    nTB[4] = 4;
    TB[4] = (double *) malloc(4 * sizeof(double));
    TBid[4] = (int *) malloc(4 * sizeof(int));
    TBid[4][0] = 0; TB[4][0] = 2.5; // H2
    TBid[4][1] = 4; TB[4][1] = 12; // H2O
    TBid[4][2] = 8; TB[4][2] = 1.8999999999999999; // CO
    TBid[4][3] = 9; TB[4][3] = 3.7999999999999998; // CO2

    // (6):  O + H + M <=> OH + M
    fwd_A[5]     = 4.714e+18;
    fwd_beta[5]  = -1;
    fwd_Ea[5]    = 0;
    prefactor_units[5]  = 1.0000000000000002e-12;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = 1e-12;
    is_PD[5] = 0;
    nTB[5] = 4;
    TB[5] = (double *) malloc(4 * sizeof(double));
    TBid[5] = (int *) malloc(4 * sizeof(int));
    TBid[5][0] = 0; TB[5][0] = 2.5; // H2
    TBid[5][1] = 4; TB[5][1] = 12; // H2O
    TBid[5][2] = 8; TB[5][2] = 1.8999999999999999; // CO
    TBid[5][3] = 9; TB[5][3] = 3.7999999999999998; // CO2

    // (7):  H + OH + M <=> H2O + M
    fwd_A[6]     = 3.8000000000000004e+22;
    fwd_beta[6]  = -2;
    fwd_Ea[6]    = 0;
    prefactor_units[6]  = 1.0000000000000002e-12;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = 1e-12;
    is_PD[6] = 0;
    nTB[6] = 4;
    TB[6] = (double *) malloc(4 * sizeof(double));
    TBid[6] = (int *) malloc(4 * sizeof(int));
    TBid[6][0] = 0; TB[6][0] = 2.5; // H2
    TBid[6][1] = 4; TB[6][1] = 12; // H2O
    TBid[6][2] = 8; TB[6][2] = 1.8999999999999999; // CO
    TBid[6][3] = 9; TB[6][3] = 3.7999999999999998; // CO2

    // (8):  H + O2 (+M) <=> HO2 (+M)
    fwd_A[0]     = 1475000000000;
    fwd_beta[0]  = 0.59999999999999998;
    fwd_Ea[0]    = 0;
    low_A[0]     = 6.366e+20;
    low_beta[0]  = -1.72;
    low_Ea[0]    = 524.79999999999995;
    troe_a[0]    = 0.80000000000000004;
    troe_Tsss[0] = 1.0000000000000001e-30;
    troe_Ts[0]   = 1e+30;
    troe_len[0]  = 3;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = 1e-12;
    is_PD[0] = 1;
    nTB[0] = 5;
    TB[0] = (double *) malloc(5 * sizeof(double));
    TBid[0] = (int *) malloc(5 * sizeof(int));
    TBid[0][0] = 0; TB[0][0] = 2; // H2
    TBid[0][1] = 4; TB[0][1] = 11; // H2O
    TBid[0][2] = 1; TB[0][2] = 0.78000000000000003; // O2
    TBid[0][3] = 8; TB[0][3] = 1.8999999999999999; // CO
    TBid[0][4] = 9; TB[0][4] = 3.7999999999999998; // CO2

    // (9):  HO2 + H <=> H2 + O2
    fwd_A[12]     = 16600000000000;
    fwd_beta[12]  = 0;
    fwd_Ea[12]    = 823;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = 1e-12;
    is_PD[12] = 0;
    nTB[12] = 0;

    // (10):  HO2 + H <=> OH + OH
    fwd_A[13]     = 70790000000000;
    fwd_beta[13]  = 0;
    fwd_Ea[13]    = 295;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = 1e-12;
    is_PD[13] = 0;
    nTB[13] = 0;

    // (11):  HO2 + O <=> O2 + OH
    fwd_A[14]     = 32500000000000;
    fwd_beta[14]  = 0;
    fwd_Ea[14]    = 0;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = 1e-12;
    is_PD[14] = 0;
    nTB[14] = 0;

    // (12):  HO2 + OH <=> H2O + O2
    fwd_A[15]     = 28900000000000;
    fwd_beta[15]  = 0;
    fwd_Ea[15]    = -497;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = 1e-12;
    is_PD[15] = 0;
    nTB[15] = 0;

    // (13):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[16]     = 420000000000000;
    fwd_beta[16]  = 0;
    fwd_Ea[16]    = 11982;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = 1e-12;
    is_PD[16] = 0;
    nTB[16] = 0;

    // (14):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[17]     = 130000000000;
    fwd_beta[17]  = 0;
    fwd_Ea[17]    = -1629.3;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = 1e-12;
    is_PD[17] = 0;
    nTB[17] = 0;

    // (15):  H2O2 (+M) <=> OH + OH (+M)
    fwd_A[1]     = 295100000000000;
    fwd_beta[1]  = 0;
    fwd_Ea[1]    = 48430;
    low_A[1]     = 1.202e+17;
    low_beta[1]  = 0;
    low_Ea[1]    = 45500;
    troe_a[1]    = 0.5;
    troe_Tsss[1] = 1.0000000000000001e-30;
    troe_Ts[1]   = 1e+30;
    troe_len[1]  = 3;
    prefactor_units[1]  = 1;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = 1e-6;
    is_PD[1] = 1;
    nTB[1] = 4;
    TB[1] = (double *) malloc(4 * sizeof(double));
    TBid[1] = (int *) malloc(4 * sizeof(int));
    TBid[1][0] = 0; TB[1][0] = 2.5; // H2
    TBid[1][1] = 4; TB[1][1] = 12; // H2O
    TBid[1][2] = 8; TB[1][2] = 1.8999999999999999; // CO
    TBid[1][3] = 9; TB[1][3] = 3.7999999999999998; // CO2

    // (16):  H2O2 + H <=> H2O + OH
    fwd_A[18]     = 24100000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 3970;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = 1e-12;
    is_PD[18] = 0;
    nTB[18] = 0;

    // (17):  H2O2 + H <=> HO2 + H2
    fwd_A[19]     = 48200000000000;
    fwd_beta[19]  = 0;
    fwd_Ea[19]    = 7950;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = 1e-12;
    is_PD[19] = 0;
    nTB[19] = 0;

    // (18):  H2O2 + O <=> OH + HO2
    fwd_A[20]     = 9550000;
    fwd_beta[20]  = 2;
    fwd_Ea[20]    = 3970;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = 1e-12;
    is_PD[20] = 0;
    nTB[20] = 0;

    // (19):  H2O2 + OH <=> HO2 + H2O
    fwd_A[21]     = 1000000000000;
    fwd_beta[21]  = 0;
    fwd_Ea[21]    = 0;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = 1e-12;
    is_PD[21] = 0;
    nTB[21] = 0;

    // (20):  H2O2 + OH <=> HO2 + H2O
    fwd_A[22]     = 580000000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 9557;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = 1e-12;
    is_PD[22] = 0;
    nTB[22] = 0;

    // (21):  CO + O (+M) <=> CO2 (+M)
    fwd_A[2]     = 18000000000;
    fwd_beta[2]  = 0;
    fwd_Ea[2]    = 2384;
    low_A[2]     = 1.5500000000000001e+24;
    low_beta[2]  = -2.79;
    low_Ea[2]    = 4191;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = 1e-12;
    is_PD[2] = 1;
    nTB[2] = 4;
    TB[2] = (double *) malloc(4 * sizeof(double));
    TBid[2] = (int *) malloc(4 * sizeof(int));
    TBid[2][0] = 0; TB[2][0] = 2.5; // H2
    TBid[2][1] = 4; TB[2][1] = 12; // H2O
    TBid[2][2] = 8; TB[2][2] = 1.8999999999999999; // CO
    TBid[2][3] = 9; TB[2][3] = 3.7999999999999998; // CO2

    // (22):  CO + O2 <=> CO2 + O
    fwd_A[23]     = 2530000000000;
    fwd_beta[23]  = 0;
    fwd_Ea[23]    = 47700;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = 1e-12;
    is_PD[23] = 0;
    nTB[23] = 0;

    // (23):  CO + HO2 <=> CO2 + OH
    fwd_A[24]     = 30100000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 23000;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = 1e-12;
    is_PD[24] = 0;
    nTB[24] = 0;

    // (24):  CO + OH <=> CO2 + H
    fwd_A[25]     = 222900;
    fwd_beta[25]  = 1.8899999999999999;
    fwd_Ea[25]    = -1158.7;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = 1e-12;
    is_PD[25] = 0;
    nTB[25] = 0;

    // (25):  HCO + M <=> H + CO + M
    fwd_A[7]     = 474850000000;
    fwd_beta[7]  = 0.65900000000000003;
    fwd_Ea[7]    = 14874;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = 1e-6;
    is_PD[7] = 0;
    nTB[7] = 4;
    TB[7] = (double *) malloc(4 * sizeof(double));
    TBid[7] = (int *) malloc(4 * sizeof(int));
    TBid[7][0] = 0; TB[7][0] = 2.5; // H2
    TBid[7][1] = 4; TB[7][1] = 6; // H2O
    TBid[7][2] = 8; TB[7][2] = 1.8999999999999999; // CO
    TBid[7][3] = 9; TB[7][3] = 3.7999999999999998; // CO2

    // (26):  HCO + O2 <=> CO + HO2
    fwd_A[26]     = 7580000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 410;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = 1e-12;
    is_PD[26] = 0;
    nTB[26] = 0;

    // (27):  HCO + H <=> CO + H2
    fwd_A[27]     = 72300000000000;
    fwd_beta[27]  = 0;
    fwd_Ea[27]    = 0;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = 1e-12;
    is_PD[27] = 0;
    nTB[27] = 0;

    // (28):  HCO + O <=> CO2 + H
    fwd_A[28]     = 30000000000000;
    fwd_beta[28]  = 0;
    fwd_Ea[28]    = 0;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = 1e-12;
    is_PD[28] = 0;
    nTB[28] = 0;

    SetAllDefaults();
}



/*A few mechanism parameters */
void CKINDX(int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 4;
    *kk = 12;
    *ii = 29;
    *nfit = -1; /*Why do you need this anyway ?  */
}



/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double *  rval, int * kerr, int lenline )
{
    int n,i; /*Loop Counters */
    char cstr[1000];
    char *saveptr;
    char *p; /*String Tokens */
    /* Strip Comments  */
    for (i=0; i<lenline; ++i) {
        if (line[i]=='!') {
            break;
        }
        cstr[i] = line[i];
    }
    cstr[i] = '\0';

    p = strtok_r(cstr," ", &saveptr);
    if (!p) {
        *nval = 0;
        *kerr = 1;
        return;
    }
    for (n=0; n<*nexp; ++n) {
        rval[n] = atof(p);
        p = strtok_r(NULL, " ", &saveptr);
        if (!p) break;
    }
    *nval = n+1;
    if (*nval < *nexp) *kerr = 1;
    return;
}


/* cksnum... for parsing strings  */
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double *  rval, int * kerr, int lenline, int lenkray)
{
    /*Not done yet ... */
}


/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*4; i++) {
        kname[i] = ' ';
    }

    /* C  */
    kname[ 0*lenkname + 0 ] = 'C';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 2*lenkname + 0 ] = 'O';
    kname[ 2*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 3*lenkname + 0 ] = 'N';
    kname[ 3*lenkname + 1 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*12; i++) {
        kname[i] = ' ';
    }

    /* H2  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* O2  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = '2';
    kname[ 1*lenkname + 2 ] = ' ';

    /* O  */
    kname[ 2*lenkname + 0 ] = 'O';
    kname[ 2*lenkname + 1 ] = ' ';

    /* OH  */
    kname[ 3*lenkname + 0 ] = 'O';
    kname[ 3*lenkname + 1 ] = 'H';
    kname[ 3*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 4*lenkname + 0 ] = 'H';
    kname[ 4*lenkname + 1 ] = '2';
    kname[ 4*lenkname + 2 ] = 'O';
    kname[ 4*lenkname + 3 ] = ' ';

    /* H  */
    kname[ 5*lenkname + 0 ] = 'H';
    kname[ 5*lenkname + 1 ] = ' ';

    /* HO2  */
    kname[ 6*lenkname + 0 ] = 'H';
    kname[ 6*lenkname + 1 ] = 'O';
    kname[ 6*lenkname + 2 ] = '2';
    kname[ 6*lenkname + 3 ] = ' ';

    /* H2O2  */
    kname[ 7*lenkname + 0 ] = 'H';
    kname[ 7*lenkname + 1 ] = '2';
    kname[ 7*lenkname + 2 ] = 'O';
    kname[ 7*lenkname + 3 ] = '2';
    kname[ 7*lenkname + 4 ] = ' ';

    /* CO  */
    kname[ 8*lenkname + 0 ] = 'C';
    kname[ 8*lenkname + 1 ] = 'O';
    kname[ 8*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 9*lenkname + 0 ] = 'C';
    kname[ 9*lenkname + 1 ] = 'O';
    kname[ 9*lenkname + 2 ] = '2';
    kname[ 9*lenkname + 3 ] = ' ';

    /* HCO  */
    kname[ 10*lenkname + 0 ] = 'H';
    kname[ 10*lenkname + 1 ] = 'C';
    kname[ 10*lenkname + 2 ] = 'O';
    kname[ 10*lenkname + 3 ] = ' ';

    /* N2  */
    kname[ 11*lenkname + 0 ] = 'N';
    kname[ 11*lenkname + 1 ] = '2';
    kname[ 11*lenkname + 2 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(double *  ru, double *  ruc, double *  pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double *  rho, double *  T, double *  x, double *  P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*1.007970; /*H */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*44.009950; /*CO2 */
    XW += x[10]*29.018520; /*HCO */
    XW += x[11]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*H */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*CO */
    YOW += y[9]*imw[9]; /*CO2 */
    YOW += y[10]*imw[10]; /*HCO */
    YOW += y[11]*imw[11]; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<12; n++) {
        for (int i=0; i<(*np); i++) {
            YOW[i] += y[n*(*np)+i] * imw[n];
        }
    }

    for (int i=0; i<(*np); i++) {
        P[i] = rho[i] * 8.31451e+07 * T[i] * YOW[i]; /*P = rho*R*T/W */
    }

    return;
}


/*Compute P = rhoRT/W(c) */
void CKPC(double *  rho, double *  T, double *  c,  double *  P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*1.007970; /*H */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.010550; /*CO */
    W += c[9]*44.009950; /*CO2 */
    W += c[10]*29.018520; /*HCO */
    W += c[11]*28.013400; /*N2 */

    for (id = 0; id < 12; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double *  P, double *  T, double *  x,  double *  rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*1.007970; /*H */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*44.009950; /*CO2 */
    XW += x[10]*29.018520; /*HCO */
    XW += x[11]*28.013400; /*N2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double *  P, double *  T, double *  y,  double *  rho)
{
    double YOW = 0;
    double tmp[12];

    for (int i = 0; i < 12; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 12; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31451e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double *  P, double *  T, double *  c,  double *  rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*1.007970; /*H */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.010550; /*CO */
    W += c[9]*44.009950; /*CO2 */
    W += c[10]*29.018520; /*HCO */
    W += c[11]*28.013400; /*N2 */

    for (id = 0; id < 12; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31451e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT( double *  wt)
{
    molecularWeight(wt);
}


/*get atomic weight for all elements */
void CKAWT( double *  awt)
{
    atomicWeight(awt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWY(double *  y,  double *  wtm)
{
    double YOW = 0;
    double tmp[12];

    for (int i = 0; i < 12; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 12; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double *  x,  double *  wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*1.007970; /*H */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*44.009950; /*CO2 */
    XW += x[10]*29.018520; /*HCO */
    XW += x[11]*28.013400; /*N2 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void CKMMWC(double *  c,  double *  wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*1.007970; /*H */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.010550; /*CO */
    W += c[9]*44.009950; /*CO2 */
    W += c[10]*29.018520; /*HCO */
    W += c[11]*28.013400; /*N2 */

    for (id = 0; id < 12; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
void CKYTX(double *  y,  double *  x)
{
    double YOW = 0;
    double tmp[12];

    for (int i = 0; i < 12; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 12; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 12; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int *  np, double *  y,  double *  x)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<12; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<12; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double *  P, double *  T, double *  y,  double *  c)
{
    double YOW = 0;
    double PWORT;

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 12; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 12; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 12; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double *  rho, double *  T, double *  y,  double *  c)
{
    for (int i = 0; i < 12; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double *  x,  double *  y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*1.007970; /*H */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*44.009950; /*CO2 */
    XW += x[10]*29.018520; /*HCO */
    XW += x[11]*28.013400; /*N2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*2.015940*XWinv; 
    y[1] = x[1]*31.998800*XWinv; 
    y[2] = x[2]*15.999400*XWinv; 
    y[3] = x[3]*17.007370*XWinv; 
    y[4] = x[4]*18.015340*XWinv; 
    y[5] = x[5]*1.007970*XWinv; 
    y[6] = x[6]*33.006770*XWinv; 
    y[7] = x[7]*34.014740*XWinv; 
    y[8] = x[8]*28.010550*XWinv; 
    y[9] = x[9]*44.009950*XWinv; 
    y[10] = x[10]*29.018520*XWinv; 
    y[11] = x[11]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double *  P, double *  T, double *  x,  double *  c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 12; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCR(double *  rho, double *  T, double *  x, double *  c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*1.007970; /*H */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*44.009950; /*CO2 */
    XW += x[10]*29.018520; /*HCO */
    XW += x[11]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 12; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
void CKCTX(double *  c, double *  x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 12; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 12; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double *  c, double *  y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*2.015940; /*H2 */
    CW += c[1]*31.998800; /*O2 */
    CW += c[2]*15.999400; /*O */
    CW += c[3]*17.007370; /*OH */
    CW += c[4]*18.015340; /*H2O */
    CW += c[5]*1.007970; /*H */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*28.010550; /*CO */
    CW += c[9]*44.009950; /*CO2 */
    CW += c[10]*29.018520; /*HCO */
    CW += c[11]*28.013400; /*N2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*2.015940*CWinv; 
    y[1] = c[1]*31.998800*CWinv; 
    y[2] = c[2]*15.999400*CWinv; 
    y[3] = c[3]*17.007370*CWinv; 
    y[4] = c[4]*18.015340*CWinv; 
    y[5] = c[5]*1.007970*CWinv; 
    y[6] = c[6]*33.006770*CWinv; 
    y[7] = c[7]*34.014740*CWinv; 
    y[8] = c[8]*28.010550*CWinv; 
    y[9] = c[9]*44.009950*CWinv; 
    y[10] = c[10]*29.018520*CWinv; 
    y[11] = c[11]*28.013400*CWinv; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void CKCPOR(double *  T, double *  cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void CKHORT(double *  T, double *  hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void CKSOR(double *  T, double *  sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void CKCVML(double *  T,  double *  cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        cvml[id] *= 8.31451e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void CKCPML(double *  T,  double *  cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        cpml[id] *= 8.31451e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double *  T,  double *  uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void CKHML(double *  T,  double *  hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void CKGML(double *  T,  double *  gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void CKAML(double *  T,  double *  aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void CKSML(double *  T,  double *  sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        sml[id] *= 8.31451e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void CKCVMS(double *  T,  double *  cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 4.124383662212169e+07; /*H2 */
    cvms[1] *= 2.598381814318037e+06; /*O2 */
    cvms[2] *= 5.196763628636074e+06; /*O */
    cvms[3] *= 4.888768810227566e+06; /*OH */
    cvms[4] *= 4.615239012974499e+06; /*H2O */
    cvms[5] *= 8.248767324424338e+07; /*H */
    cvms[6] *= 2.519031701678171e+06; /*HO2 */
    cvms[7] *= 2.444384405113783e+06; /*H2O2 */
    cvms[8] *= 2.968349425484326e+06; /*CO */
    cvms[9] *= 1.889234139098090e+06; /*CO2 */
    cvms[10] *= 2.865242610581105e+06; /*HCO */
    cvms[11] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double *  T,  double *  cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124383662212169e+07; /*H2 */
    cpms[1] *= 2.598381814318037e+06; /*O2 */
    cpms[2] *= 5.196763628636074e+06; /*O */
    cpms[3] *= 4.888768810227566e+06; /*OH */
    cpms[4] *= 4.615239012974499e+06; /*H2O */
    cpms[5] *= 8.248767324424338e+07; /*H */
    cpms[6] *= 2.519031701678171e+06; /*HO2 */
    cpms[7] *= 2.444384405113783e+06; /*H2O2 */
    cpms[8] *= 2.968349425484326e+06; /*CO */
    cpms[9] *= 1.889234139098090e+06; /*CO2 */
    cpms[10] *= 2.865242610581105e+06; /*HCO */
    cpms[11] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double *  T,  double *  ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 12; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double *  T,  double *  hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 12; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, double *  T,  double *  hms)
{
    double tc[5], h[12];

    for (int i=0; i<(*np); i++) {
        tc[0] = 0.0;
        tc[1] = T[i];
        tc[2] = T[i]*T[i];
        tc[3] = T[i]*T[i]*T[i];
        tc[4] = T[i]*T[i]*T[i]*T[i];

        speciesEnthalpy(h, tc);

        hms[0*(*np)+i] = h[0];
        hms[1*(*np)+i] = h[1];
        hms[2*(*np)+i] = h[2];
        hms[3*(*np)+i] = h[3];
        hms[4*(*np)+i] = h[4];
        hms[5*(*np)+i] = h[5];
        hms[6*(*np)+i] = h[6];
        hms[7*(*np)+i] = h[7];
        hms[8*(*np)+i] = h[8];
        hms[9*(*np)+i] = h[9];
        hms[10*(*np)+i] = h[10];
        hms[11*(*np)+i] = h[11];
    }

    for (int n=0; n<12; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31451e+07 * T[i] * imw[n];
        }
    }
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double *  T,  double *  gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gms, tc);
    for (int i = 0; i < 12; i++)
    {
        gms[i] *= RT*imw[i];
    }
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double *  T,  double *  ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(ams, tc);
    for (int i = 0; i < 12; i++)
    {
        ams[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double *  T,  double *  sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 4.124383662212169e+07; /*H2 */
    sms[1] *= 2.598381814318037e+06; /*O2 */
    sms[2] *= 5.196763628636074e+06; /*O */
    sms[3] *= 4.888768810227566e+06; /*OH */
    sms[4] *= 4.615239012974499e+06; /*H2O */
    sms[5] *= 8.248767324424338e+07; /*H */
    sms[6] *= 2.519031701678171e+06; /*HO2 */
    sms[7] *= 2.444384405113783e+06; /*H2O2 */
    sms[8] *= 2.968349425484326e+06; /*CO */
    sms[9] *= 1.889234139098090e+06; /*CO2 */
    sms[10] *= 2.865242610581105e+06; /*HCO */
    sms[11] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *  T, double *  x,  double *  cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[12]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 12; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void CKCPBS(double *  T, double *  y,  double *  cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[12], tresult[12]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 12; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 12; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *  T, double *  x,  double *  cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[12]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 12; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
void CKCVBS(double *  T, double *  y,  double *  cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[12]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*O2 */
    result += cvor[2]*y[2]*imw[2]; /*O */
    result += cvor[3]*y[3]*imw[3]; /*OH */
    result += cvor[4]*y[4]*imw[4]; /*H2O */
    result += cvor[5]*y[5]*imw[5]; /*H */
    result += cvor[6]*y[6]*imw[6]; /*HO2 */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*CO */
    result += cvor[9]*y[9]*imw[9]; /*CO2 */
    result += cvor[10]*y[10]*imw[10]; /*HCO */
    result += cvor[11]*y[11]*imw[11]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *  T, double *  x,  double *  hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[12]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 12; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
void CKHBMS(double *  T, double *  y,  double *  hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[12], tmp[12]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 12; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 12; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *  T, double *  x,  double *  ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[12]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 12; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
void CKUBMS(double *  T, double *  y,  double *  ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[12]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*O2 */
    result += y[2]*ums[2]*imw[2]; /*O */
    result += y[3]*ums[3]*imw[3]; /*OH */
    result += y[4]*ums[4]*imw[4]; /*H2O */
    result += y[5]*ums[5]*imw[5]; /*H */
    result += y[6]*ums[6]*imw[6]; /*HO2 */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*CO */
    result += y[9]*ums[9]*imw[9]; /*CO2 */
    result += y[10]*ums[10]*imw[10]; /*HCO */
    result += y[11]*ums[11]*imw[11]; /*N2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
void CKSBML(double *  P, double *  T, double *  x,  double *  sbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[12]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 12; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31451e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double *  P, double *  T, double *  y,  double *  sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[12]; /* temporary storage */
    double x[12]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*H */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*CO */
    YOW += y[9]*imw[9]; /*CO2 */
    YOW += y[10]*imw[10]; /*HCO */
    YOW += y[11]*imw[11]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(1.007970*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.010550*YOW); 
    x[9] = y[9]/(44.009950*YOW); 
    x[10] = y[10]/(29.018520*YOW); 
    x[11] = y[11]/(28.013400*YOW); 
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
    /*Scale by R/W */
    *sbms = result * 8.31451e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void CKGBML(double *  P, double *  T, double *  x,  double *  gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[12]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 12; ++id) {
        result += x[id]*(gort[id]+log((x[id]+1e-100))+logPratio);
    }

    *gbml = result * RT;
}


/*Returns mixture gibbs free energy in mass units */
void CKGBMS(double *  P, double *  T, double *  y,  double *  gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[12]; /* temporary storage */
    double x[12]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*H */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*CO */
    YOW += y[9]*imw[9]; /*CO2 */
    YOW += y[10]*imw[10]; /*HCO */
    YOW += y[11]*imw[11]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(1.007970*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.010550*YOW); 
    x[9] = y[9]/(44.009950*YOW); 
    x[10] = y[10]/(29.018520*YOW); 
    x[11] = y[11]/(28.013400*YOW); 
    gibbs(gort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(gort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(gort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(gort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(gort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(gort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(gort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(gort[6]+log((x[6]+1e-100))+logPratio);
    result += x[7]*(gort[7]+log((x[7]+1e-100))+logPratio);
    result += x[8]*(gort[8]+log((x[8]+1e-100))+logPratio);
    result += x[9]*(gort[9]+log((x[9]+1e-100))+logPratio);
    result += x[10]*(gort[10]+log((x[10]+1e-100))+logPratio);
    result += x[11]*(gort[11]+log((x[11]+1e-100))+logPratio);
    /*Scale by RT/W */
    *gbms = result * RT * YOW;
}


/*Returns mean helmholtz free energy in molar units */
void CKABML(double *  P, double *  T, double *  x,  double *  abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[12]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 12; ++id) {
        result += x[id]*(aort[id]+log((x[id]+1e-100))+logPratio);
    }

    *abml = result * RT;
}


/*Returns mixture helmholtz free energy in mass units */
void CKABMS(double *  P, double *  T, double *  y,  double *  abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[12]; /* temporary storage */
    double x[12]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*H */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*CO */
    YOW += y[9]*imw[9]; /*CO2 */
    YOW += y[10]*imw[10]; /*HCO */
    YOW += y[11]*imw[11]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(1.007970*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.010550*YOW); 
    x[9] = y[9]/(44.009950*YOW); 
    x[10] = y[10]/(29.018520*YOW); 
    x[11] = y[11]/(28.013400*YOW); 
    helmholtz(aort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(aort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(aort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(aort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(aort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(aort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(aort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(aort[6]+log((x[6]+1e-100))+logPratio);
    result += x[7]*(aort[7]+log((x[7]+1e-100))+logPratio);
    result += x[8]*(aort[8]+log((x[8]+1e-100))+logPratio);
    result += x[9]*(aort[9]+log((x[9]+1e-100))+logPratio);
    result += x[10]*(aort[10]+log((x[10]+1e-100))+logPratio);
    result += x[11]*(aort[11]+log((x[11]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double *  T, double *  C,  double *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 12; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double *  P, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[12]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*H */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*CO */
    YOW += y[9]*imw[9]; /*CO2 */
    YOW += y[10]*imw[10]; /*HCO */
    YOW += y[11]*imw[11]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
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

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double *  P, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[12]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 12; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double *  rho, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[12]; /*temporary storage */
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

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int *  np, double *  rho, double *  T,
	    double *  y,
	    double *  wdot)
{
    double c[12*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<12; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<12*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double *  rho, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[12]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*1.007970; /*H */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*44.009950; /*CO2 */
    XW += x[10]*29.018520; /*HCO */
    XW += x[11]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 12; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double *  T, double *  C, double *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 12; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 12; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 29; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r)
{
    int id; /*loop counter */
    double c[12]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 12; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double *  P, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[12]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*H */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*CO */
    YOW += y[9]*imw[9]; /*CO2 */
    YOW += y[10]*imw[10]; /*HCO */
    YOW += y[11]*imw[11]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
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

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double *  P, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[12]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 12; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double *  rho, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[12]; /*temporary storage */
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

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double *  rho, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[12]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*1.007970; /*H */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*44.009950; /*CO2 */
    XW += x[10]*29.018520; /*HCO */
    XW += x[11]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 12; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim,  int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 12 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    nuki[ 5 * kd + 0 ] += -1 ;
    nuki[ 1 * kd + 0 ] += -1 ;
    nuki[ 6 * kd + 0 ] += +1 ;

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    nuki[ 7 * kd + 1 ] += -1 ;
    nuki[ 3 * kd + 1 ] += +1 ;
    nuki[ 3 * kd + 1 ] += +1 ;

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    nuki[ 8 * kd + 2 ] += -1 ;
    nuki[ 2 * kd + 2 ] += -1 ;
    nuki[ 9 * kd + 2 ] += +1 ;

    /*reaction 4: H2 + M <=> H + H + M */
    nuki[ 0 * kd + 3 ] += -1 ;
    nuki[ 5 * kd + 3 ] += +1 ;
    nuki[ 5 * kd + 3 ] += +1 ;

    /*reaction 5: O + O + M <=> O2 + M */
    nuki[ 2 * kd + 4 ] += -1 ;
    nuki[ 2 * kd + 4 ] += -1 ;
    nuki[ 1 * kd + 4 ] += +1 ;

    /*reaction 6: O + H + M <=> OH + M */
    nuki[ 2 * kd + 5 ] += -1 ;
    nuki[ 5 * kd + 5 ] += -1 ;
    nuki[ 3 * kd + 5 ] += +1 ;

    /*reaction 7: H + OH + M <=> H2O + M */
    nuki[ 5 * kd + 6 ] += -1 ;
    nuki[ 3 * kd + 6 ] += -1 ;
    nuki[ 4 * kd + 6 ] += +1 ;

    /*reaction 8: HCO + M <=> H + CO + M */
    nuki[ 10 * kd + 7 ] += -1 ;
    nuki[ 5 * kd + 7 ] += +1 ;
    nuki[ 8 * kd + 7 ] += +1 ;

    /*reaction 9: H + O2 <=> O + OH */
    nuki[ 5 * kd + 8 ] += -1 ;
    nuki[ 1 * kd + 8 ] += -1 ;
    nuki[ 2 * kd + 8 ] += +1 ;
    nuki[ 3 * kd + 8 ] += +1 ;

    /*reaction 10: O + H2 <=> H + OH */
    nuki[ 2 * kd + 9 ] += -1 ;
    nuki[ 0 * kd + 9 ] += -1 ;
    nuki[ 5 * kd + 9 ] += +1 ;
    nuki[ 3 * kd + 9 ] += +1 ;

    /*reaction 11: H2 + OH <=> H2O + H */
    nuki[ 0 * kd + 10 ] += -1 ;
    nuki[ 3 * kd + 10 ] += -1 ;
    nuki[ 4 * kd + 10 ] += +1 ;
    nuki[ 5 * kd + 10 ] += +1 ;

    /*reaction 12: O + H2O <=> OH + OH */
    nuki[ 2 * kd + 11 ] += -1 ;
    nuki[ 4 * kd + 11 ] += -1 ;
    nuki[ 3 * kd + 11 ] += +1 ;
    nuki[ 3 * kd + 11 ] += +1 ;

    /*reaction 13: HO2 + H <=> H2 + O2 */
    nuki[ 6 * kd + 12 ] += -1 ;
    nuki[ 5 * kd + 12 ] += -1 ;
    nuki[ 0 * kd + 12 ] += +1 ;
    nuki[ 1 * kd + 12 ] += +1 ;

    /*reaction 14: HO2 + H <=> OH + OH */
    nuki[ 6 * kd + 13 ] += -1 ;
    nuki[ 5 * kd + 13 ] += -1 ;
    nuki[ 3 * kd + 13 ] += +1 ;
    nuki[ 3 * kd + 13 ] += +1 ;

    /*reaction 15: HO2 + O <=> O2 + OH */
    nuki[ 6 * kd + 14 ] += -1 ;
    nuki[ 2 * kd + 14 ] += -1 ;
    nuki[ 1 * kd + 14 ] += +1 ;
    nuki[ 3 * kd + 14 ] += +1 ;

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    nuki[ 6 * kd + 15 ] += -1 ;
    nuki[ 3 * kd + 15 ] += -1 ;
    nuki[ 4 * kd + 15 ] += +1 ;
    nuki[ 1 * kd + 15 ] += +1 ;

    /*reaction 17: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 6 * kd + 16 ] += -1 ;
    nuki[ 6 * kd + 16 ] += -1 ;
    nuki[ 7 * kd + 16 ] += +1 ;
    nuki[ 1 * kd + 16 ] += +1 ;

    /*reaction 18: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 6 * kd + 17 ] += -1 ;
    nuki[ 6 * kd + 17 ] += -1 ;
    nuki[ 7 * kd + 17 ] += +1 ;
    nuki[ 1 * kd + 17 ] += +1 ;

    /*reaction 19: H2O2 + H <=> H2O + OH */
    nuki[ 7 * kd + 18 ] += -1 ;
    nuki[ 5 * kd + 18 ] += -1 ;
    nuki[ 4 * kd + 18 ] += +1 ;
    nuki[ 3 * kd + 18 ] += +1 ;

    /*reaction 20: H2O2 + H <=> HO2 + H2 */
    nuki[ 7 * kd + 19 ] += -1 ;
    nuki[ 5 * kd + 19 ] += -1 ;
    nuki[ 6 * kd + 19 ] += +1 ;
    nuki[ 0 * kd + 19 ] += +1 ;

    /*reaction 21: H2O2 + O <=> OH + HO2 */
    nuki[ 7 * kd + 20 ] += -1 ;
    nuki[ 2 * kd + 20 ] += -1 ;
    nuki[ 3 * kd + 20 ] += +1 ;
    nuki[ 6 * kd + 20 ] += +1 ;

    /*reaction 22: H2O2 + OH <=> HO2 + H2O */
    nuki[ 7 * kd + 21 ] += -1 ;
    nuki[ 3 * kd + 21 ] += -1 ;
    nuki[ 6 * kd + 21 ] += +1 ;
    nuki[ 4 * kd + 21 ] += +1 ;

    /*reaction 23: H2O2 + OH <=> HO2 + H2O */
    nuki[ 7 * kd + 22 ] += -1 ;
    nuki[ 3 * kd + 22 ] += -1 ;
    nuki[ 6 * kd + 22 ] += +1 ;
    nuki[ 4 * kd + 22 ] += +1 ;

    /*reaction 24: CO + O2 <=> CO2 + O */
    nuki[ 8 * kd + 23 ] += -1 ;
    nuki[ 1 * kd + 23 ] += -1 ;
    nuki[ 9 * kd + 23 ] += +1 ;
    nuki[ 2 * kd + 23 ] += +1 ;

    /*reaction 25: CO + HO2 <=> CO2 + OH */
    nuki[ 8 * kd + 24 ] += -1 ;
    nuki[ 6 * kd + 24 ] += -1 ;
    nuki[ 9 * kd + 24 ] += +1 ;
    nuki[ 3 * kd + 24 ] += +1 ;

    /*reaction 26: CO + OH <=> CO2 + H */
    nuki[ 8 * kd + 25 ] += -1 ;
    nuki[ 3 * kd + 25 ] += -1 ;
    nuki[ 9 * kd + 25 ] += +1 ;
    nuki[ 5 * kd + 25 ] += +1 ;

    /*reaction 27: HCO + O2 <=> CO + HO2 */
    nuki[ 10 * kd + 26 ] += -1 ;
    nuki[ 1 * kd + 26 ] += -1 ;
    nuki[ 8 * kd + 26 ] += +1 ;
    nuki[ 6 * kd + 26 ] += +1 ;

    /*reaction 28: HCO + H <=> CO + H2 */
    nuki[ 10 * kd + 27 ] += -1 ;
    nuki[ 5 * kd + 27 ] += -1 ;
    nuki[ 8 * kd + 27 ] += +1 ;
    nuki[ 0 * kd + 27 ] += +1 ;

    /*reaction 29: HCO + O <=> CO2 + H */
    nuki[ 10 * kd + 28 ] += -1 ;
    nuki[ 2 * kd + 28 ] += -1 ;
    nuki[ 9 * kd + 28 ] += +1 ;
    nuki[ 5 * kd + 28 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim,  int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 12; ++ id) {
         ncf[id] = 0; 
    }

    /*H2 */
    ncf[ 0 * kd + 1 ] = 2; /*H */

    /*O2 */
    ncf[ 1 * kd + 2 ] = 2; /*O */

    /*O */
    ncf[ 2 * kd + 2 ] = 1; /*O */

    /*OH */
    ncf[ 3 * kd + 2 ] = 1; /*O */
    ncf[ 3 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 4 * kd + 1 ] = 2; /*H */
    ncf[ 4 * kd + 2 ] = 1; /*O */

    /*H */
    ncf[ 5 * kd + 1 ] = 1; /*H */

    /*HO2 */
    ncf[ 6 * kd + 1 ] = 1; /*H */
    ncf[ 6 * kd + 2 ] = 2; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 1 ] = 2; /*H */
    ncf[ 7 * kd + 2 ] = 2; /*O */

    /*CO */
    ncf[ 8 * kd + 0 ] = 1; /*C */
    ncf[ 8 * kd + 2 ] = 1; /*O */

    /*CO2 */
    ncf[ 9 * kd + 0 ] = 1; /*C */
    ncf[ 9 * kd + 2 ] = 2; /*O */

    /*HCO */
    ncf[ 10 * kd + 1 ] = 1; /*H */
    ncf[ 10 * kd + 0 ] = 1; /*C */
    ncf[ 10 * kd + 2 ] = 1; /*O */

    /*N2 */
    ncf[ 11 * kd + 3 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE( double *  a, double *  b, double *  e)
{
    for (int i=0; i<29; ++i) {
        a[i] = fwd_A[i];
        b[i] = fwd_beta[i];
        e[i] = fwd_Ea[i];
    }

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double *  T, double *  C, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[12]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H2 + M <=> H + H + M */
    eqcon[3] *= 1e-06; 

    /*reaction 5: O + O + M <=> O2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + OH + M <=> H2O + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> H + CO + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + H2O <=> OH + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + H <=> H2 + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + H <=> OH + OH */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + O <=> O2 + OH */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + H <=> H2O + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + H <=> HO2 + H2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + O <=> OH + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: CO + O2 <=> CO2 + O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: CO + HO2 <=> CO2 + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: CO + OH <=> CO2 + H */
    /*eqcon[25] *= 1;  */

    /*reaction 27: HCO + O2 <=> CO + HO2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: HCO + H <=> CO + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HCO + O <=> CO2 + H */
    /*eqcon[28] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double *  P, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[12]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H2 + M <=> H + H + M */
    eqcon[3] *= 1e-06; 

    /*reaction 5: O + O + M <=> O2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + OH + M <=> H2O + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> H + CO + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + H2O <=> OH + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + H <=> H2 + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + H <=> OH + OH */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + O <=> O2 + OH */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + H <=> H2O + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + H <=> HO2 + H2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + O <=> OH + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: CO + O2 <=> CO2 + O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: CO + HO2 <=> CO2 + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: CO + OH <=> CO2 + H */
    /*eqcon[25] *= 1;  */

    /*reaction 27: HCO + O2 <=> CO + HO2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: HCO + H <=> CO + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HCO + O <=> CO2 + H */
    /*eqcon[28] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double *  P, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[12]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H2 + M <=> H + H + M */
    eqcon[3] *= 1e-06; 

    /*reaction 5: O + O + M <=> O2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + OH + M <=> H2O + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> H + CO + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + H2O <=> OH + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + H <=> H2 + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + H <=> OH + OH */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + O <=> O2 + OH */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + H <=> H2O + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + H <=> HO2 + H2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + O <=> OH + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: CO + O2 <=> CO2 + O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: CO + HO2 <=> CO2 + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: CO + OH <=> CO2 + H */
    /*eqcon[25] *= 1;  */

    /*reaction 27: HCO + O2 <=> CO + HO2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: HCO + H <=> CO + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HCO + O <=> CO2 + H */
    /*eqcon[28] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double *  rho, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[12]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H2 + M <=> H + H + M */
    eqcon[3] *= 1e-06; 

    /*reaction 5: O + O + M <=> O2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + OH + M <=> H2O + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> H + CO + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + H2O <=> OH + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + H <=> H2 + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + H <=> OH + OH */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + O <=> O2 + OH */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + H <=> H2O + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + H <=> HO2 + H2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + O <=> OH + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: CO + O2 <=> CO2 + O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: CO + HO2 <=> CO2 + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: CO + OH <=> CO2 + H */
    /*eqcon[25] *= 1;  */

    /*reaction 27: HCO + O2 <=> CO + HO2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: HCO + H <=> CO + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HCO + O <=> CO2 + H */
    /*eqcon[28] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double *  rho, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[12]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H2 + M <=> H + H + M */
    eqcon[3] *= 1e-06; 

    /*reaction 5: O + O + M <=> O2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + OH + M <=> H2O + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> H + CO + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + H2O <=> OH + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + H <=> H2 + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + H <=> OH + OH */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + O <=> O2 + OH */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + H <=> H2O + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + H <=> HO2 + H2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + O <=> OH + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: CO + O2 <=> CO2 + O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: CO + HO2 <=> CO2 + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: CO + OH <=> CO2 + H */
    /*eqcon[25] *= 1;  */

    /*reaction 27: HCO + O2 <=> CO + HO2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: HCO + H <=> CO + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HCO + O <=> CO2 + H */
    /*eqcon[28] *= 1;  */
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[29];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[29];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif


/*compute the production rate for each species */
void productionRate(double *  wdot, double *  sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double qdot, q_f[29], q_r[29];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 12; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[1] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[3] += qdot;
    wdot[3] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[2]-q_r[2];
    wdot[2] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[0] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[2] -= qdot;

    qdot = q_f[5]-q_r[5];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[7]-q_r[7];
    wdot[5] += qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[8]-q_r[8];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[9]-q_r[9];
    wdot[0] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[0] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;

    qdot = q_f[12]-q_r[12];
    wdot[0] += qdot;
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[6] -= qdot;

    qdot = q_f[13]-q_r[13];
    wdot[3] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[6] -= qdot;

    qdot = q_f[14]-q_r[14];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[15]-q_r[15];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[16]-q_r[16];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[18]-q_r[18];
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[7] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[0] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[21]-q_r[21];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[22]-q_r[22];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[24]-q_r[24];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[25]-q_r[25];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[26]-q_r[26];
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[27]-q_r[27];
    wdot[0] += qdot;
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    return;
}

void comp_k_f(double *  tc, double invT, double *  k_f)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<29; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double *  tc, double invT, double *  Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[12];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[1] + g_RT[5] - g_RT[6];
    Kc[1] = -g_RT[3] - g_RT[3] + g_RT[7];
    Kc[2] = g_RT[2] + g_RT[8] - g_RT[9];
    Kc[3] = g_RT[0] - g_RT[5] - g_RT[5];
    Kc[4] = -g_RT[1] + g_RT[2] + g_RT[2];
    Kc[5] = g_RT[2] - g_RT[3] + g_RT[5];
    Kc[6] = g_RT[3] - g_RT[4] + g_RT[5];
    Kc[7] = -g_RT[5] - g_RT[8] + g_RT[10];
    Kc[8] = g_RT[1] - g_RT[2] - g_RT[3] + g_RT[5];
    Kc[9] = g_RT[0] + g_RT[2] - g_RT[3] - g_RT[5];
    Kc[10] = g_RT[0] + g_RT[3] - g_RT[4] - g_RT[5];
    Kc[11] = g_RT[2] - g_RT[3] - g_RT[3] + g_RT[4];
    Kc[12] = -g_RT[0] - g_RT[1] + g_RT[5] + g_RT[6];
    Kc[13] = -g_RT[3] - g_RT[3] + g_RT[5] + g_RT[6];
    Kc[14] = -g_RT[1] + g_RT[2] - g_RT[3] + g_RT[6];
    Kc[15] = -g_RT[1] + g_RT[3] - g_RT[4] + g_RT[6];
    Kc[16] = -g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7];
    Kc[17] = -g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7];
    Kc[18] = -g_RT[3] - g_RT[4] + g_RT[5] + g_RT[7];
    Kc[19] = -g_RT[0] + g_RT[5] - g_RT[6] + g_RT[7];
    Kc[20] = g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7];
    Kc[21] = g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7];
    Kc[22] = g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7];
    Kc[23] = g_RT[1] - g_RT[2] + g_RT[8] - g_RT[9];
    Kc[24] = -g_RT[3] + g_RT[6] + g_RT[8] - g_RT[9];
    Kc[25] = g_RT[3] - g_RT[5] + g_RT[8] - g_RT[9];
    Kc[26] = g_RT[1] - g_RT[6] - g_RT[8] + g_RT[10];
    Kc[27] = -g_RT[0] + g_RT[5] - g_RT[8] + g_RT[10];
    Kc[28] = g_RT[2] - g_RT[5] - g_RT[9] + g_RT[10];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<29; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refC;
    Kc[2] *= refCinv;
    Kc[3] *= refC;
    Kc[4] *= refCinv;
    Kc[5] *= refCinv;
    Kc[6] *= refCinv;
    Kc[7] *= refC;

    return;
}

void comp_qfqr(double *  qf, double *  qr, double *  sc, double *  tc, double invT)
{

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[1]*sc[5];
    qr[0] = sc[6];

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    qf[1] = sc[7];
    qr[1] = sc[3]*sc[3];

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    qf[2] = sc[2]*sc[8];
    qr[2] = sc[9];

    /*reaction 4: H2 + M <=> H + H + M */
    qf[3] = sc[0];
    qr[3] = sc[5]*sc[5];

    /*reaction 5: O + O + M <=> O2 + M */
    qf[4] = sc[2]*sc[2];
    qr[4] = sc[1];

    /*reaction 6: O + H + M <=> OH + M */
    qf[5] = sc[2]*sc[5];
    qr[5] = sc[3];

    /*reaction 7: H + OH + M <=> H2O + M */
    qf[6] = sc[3]*sc[5];
    qr[6] = sc[4];

    /*reaction 8: HCO + M <=> H + CO + M */
    qf[7] = sc[10];
    qr[7] = sc[5]*sc[8];

    /*reaction 9: H + O2 <=> O + OH */
    qf[8] = sc[1]*sc[5];
    qr[8] = sc[2]*sc[3];

    /*reaction 10: O + H2 <=> H + OH */
    qf[9] = sc[0]*sc[2];
    qr[9] = sc[3]*sc[5];

    /*reaction 11: H2 + OH <=> H2O + H */
    qf[10] = sc[0]*sc[3];
    qr[10] = sc[4]*sc[5];

    /*reaction 12: O + H2O <=> OH + OH */
    qf[11] = sc[2]*sc[4];
    qr[11] = sc[3]*sc[3];

    /*reaction 13: HO2 + H <=> H2 + O2 */
    qf[12] = sc[5]*sc[6];
    qr[12] = sc[0]*sc[1];

    /*reaction 14: HO2 + H <=> OH + OH */
    qf[13] = sc[5]*sc[6];
    qr[13] = sc[3]*sc[3];

    /*reaction 15: HO2 + O <=> O2 + OH */
    qf[14] = sc[2]*sc[6];
    qr[14] = sc[1]*sc[3];

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    qf[15] = sc[3]*sc[6];
    qr[15] = sc[1]*sc[4];

    /*reaction 17: HO2 + HO2 <=> H2O2 + O2 */
    qf[16] = sc[6]*sc[6];
    qr[16] = sc[1]*sc[7];

    /*reaction 18: HO2 + HO2 <=> H2O2 + O2 */
    qf[17] = sc[6]*sc[6];
    qr[17] = sc[1]*sc[7];

    /*reaction 19: H2O2 + H <=> H2O + OH */
    qf[18] = sc[5]*sc[7];
    qr[18] = sc[3]*sc[4];

    /*reaction 20: H2O2 + H <=> HO2 + H2 */
    qf[19] = sc[5]*sc[7];
    qr[19] = sc[0]*sc[6];

    /*reaction 21: H2O2 + O <=> OH + HO2 */
    qf[20] = sc[2]*sc[7];
    qr[20] = sc[3]*sc[6];

    /*reaction 22: H2O2 + OH <=> HO2 + H2O */
    qf[21] = sc[3]*sc[7];
    qr[21] = sc[4]*sc[6];

    /*reaction 23: H2O2 + OH <=> HO2 + H2O */
    qf[22] = sc[3]*sc[7];
    qr[22] = sc[4]*sc[6];

    /*reaction 24: CO + O2 <=> CO2 + O */
    qf[23] = sc[1]*sc[8];
    qr[23] = sc[2]*sc[9];

    /*reaction 25: CO + HO2 <=> CO2 + OH */
    qf[24] = sc[6]*sc[8];
    qr[24] = sc[3]*sc[9];

    /*reaction 26: CO + OH <=> CO2 + H */
    qf[25] = sc[3]*sc[8];
    qr[25] = sc[5]*sc[9];

    /*reaction 27: HCO + O2 <=> CO + HO2 */
    qf[26] = sc[1]*sc[10];
    qr[26] = sc[6]*sc[8];

    /*reaction 28: HCO + H <=> CO + H2 */
    qf[27] = sc[5]*sc[10];
    qr[27] = sc[0]*sc[8];

    /*reaction 29: HCO + O <=> CO2 + H */
    qf[28] = sc[2]*sc[10];
    qr[28] = sc[5]*sc[9];

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 12; ++i) {
        mixture += sc[i];
    }

    double Corr[29];
    for (int i = 0; i < 29; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        double alpha[2];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[4] + (TB[0][2] - 1)*sc[1] + (TB[0][3] - 1)*sc[8] + (TB[0][4] - 1)*sc[9];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[4] + (TB[1][2] - 1)*sc[8] + (TB[1][3] - 1)*sc[9];
        for (int i=0; i<2; i++)
        {
            double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;
            redP = alpha[i-0] / k_f_save[i] * phase_units[i] * low_A[i] * exp(low_beta[i] * tc[0] - activation_units[i] * low_Ea[i] *invT);
            F = redP / (1.0 + redP);
            logPred = log10(redP);
            logFcent = log10(
                (fabs(troe_Tsss[i]) > 1.e-100 ? (1.-troe_a[i])*exp(-T/troe_Tsss[i]) : 0.) 
                + (fabs(troe_Ts[i]) > 1.e-100 ? troe_a[i] * exp(-T/troe_Ts[i]) : 0.) 
                + (troe_len[i] == 4 ? exp(-troe_Tss[i] * invT) : 0.) );
            troe_c = -.4 - .67 * logFcent;
            troe_n = .75 - 1.27 * logFcent;
            troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
            F_troe = pow(10., logFcent / (1.0 + troe*troe));
            Corr[i] = F * F_troe;
        }
    }

    /* Lindemann */
    {
        double alpha;
        alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[4] + (TB[2][2] - 1)*sc[8] + (TB[2][3] - 1)*sc[9];
        double redP = alpha / k_f_save[2] * phase_units[2] * low_A[2] * exp(low_beta[2] * tc[0] - activation_units[2] * low_Ea[2] * invT);
        Corr[2] = redP / (1. + redP);
    }

    /* simple three-body correction */
    {
        double alpha;
        alpha = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[4] + (TB[3][2] - 1)*sc[8] + (TB[3][3] - 1)*sc[9];
        Corr[3] = alpha;
        alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[4] + (TB[4][2] - 1)*sc[8] + (TB[4][3] - 1)*sc[9];
        Corr[4] = alpha;
        alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[4] + (TB[5][2] - 1)*sc[8] + (TB[5][3] - 1)*sc[9];
        Corr[5] = alpha;
        alpha = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[4] + (TB[6][2] - 1)*sc[8] + (TB[6][3] - 1)*sc[9];
        Corr[6] = alpha;
        alpha = mixture + (TB[7][0] - 1)*sc[0] + (TB[7][1] - 1)*sc[4] + (TB[7][2] - 1)*sc[8] + (TB[7][3] - 1)*sc[9];
        Corr[7] = alpha;
    }

    for (int i=0; i<29; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double *  wdot, double *  sc, double *  T)
{
    double k_f_s[29*npt], Kc_s[29*npt], mixture[npt], g_RT[12*npt];
    double tc[5*npt], invT[npt];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        tc[0*npt+i] = log(T[i]);
        tc[1*npt+i] = T[i];
        tc[2*npt+i] = T[i]*T[i];
        tc[3*npt+i] = T[i]*T[i]*T[i];
        tc[4*npt+i] = T[i]*T[i]*T[i]*T[i];
        invT[i] = 1.0 / T[i];
    }

    for (int i=0; i<npt; i++) {
        mixture[i] = 0.0;
    }

    for (int n=0; n<12; n++) {
        for (int i=0; i<npt; i++) {
            mixture[i] += sc[n*npt+i];
            wdot[n*npt+i] = 0.0;
        }
    }

    vcomp_k_f(npt, k_f_s, tc, invT);

    vcomp_gibbs(npt, g_RT, tc);

    vcomp_Kc(npt, Kc_s, g_RT, invT);

    vcomp_wdot(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
}

void vcomp_k_f(int npt, double *  k_f_s, double *  tc, double *  invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        k_f_s[0*npt+i] = prefactor_units[0] * fwd_A[0] * exp(fwd_beta[0] * tc[i] - activation_units[0] * fwd_Ea[0] * invT[i]);
        k_f_s[1*npt+i] = prefactor_units[1] * fwd_A[1] * exp(fwd_beta[1] * tc[i] - activation_units[1] * fwd_Ea[1] * invT[i]);
        k_f_s[2*npt+i] = prefactor_units[2] * fwd_A[2] * exp(fwd_beta[2] * tc[i] - activation_units[2] * fwd_Ea[2] * invT[i]);
        k_f_s[3*npt+i] = prefactor_units[3] * fwd_A[3] * exp(fwd_beta[3] * tc[i] - activation_units[3] * fwd_Ea[3] * invT[i]);
        k_f_s[4*npt+i] = prefactor_units[4] * fwd_A[4] * exp(fwd_beta[4] * tc[i] - activation_units[4] * fwd_Ea[4] * invT[i]);
        k_f_s[5*npt+i] = prefactor_units[5] * fwd_A[5] * exp(fwd_beta[5] * tc[i] - activation_units[5] * fwd_Ea[5] * invT[i]);
        k_f_s[6*npt+i] = prefactor_units[6] * fwd_A[6] * exp(fwd_beta[6] * tc[i] - activation_units[6] * fwd_Ea[6] * invT[i]);
        k_f_s[7*npt+i] = prefactor_units[7] * fwd_A[7] * exp(fwd_beta[7] * tc[i] - activation_units[7] * fwd_Ea[7] * invT[i]);
        k_f_s[8*npt+i] = prefactor_units[8] * fwd_A[8] * exp(fwd_beta[8] * tc[i] - activation_units[8] * fwd_Ea[8] * invT[i]);
        k_f_s[9*npt+i] = prefactor_units[9] * fwd_A[9] * exp(fwd_beta[9] * tc[i] - activation_units[9] * fwd_Ea[9] * invT[i]);
        k_f_s[10*npt+i] = prefactor_units[10] * fwd_A[10] * exp(fwd_beta[10] * tc[i] - activation_units[10] * fwd_Ea[10] * invT[i]);
        k_f_s[11*npt+i] = prefactor_units[11] * fwd_A[11] * exp(fwd_beta[11] * tc[i] - activation_units[11] * fwd_Ea[11] * invT[i]);
        k_f_s[12*npt+i] = prefactor_units[12] * fwd_A[12] * exp(fwd_beta[12] * tc[i] - activation_units[12] * fwd_Ea[12] * invT[i]);
        k_f_s[13*npt+i] = prefactor_units[13] * fwd_A[13] * exp(fwd_beta[13] * tc[i] - activation_units[13] * fwd_Ea[13] * invT[i]);
        k_f_s[14*npt+i] = prefactor_units[14] * fwd_A[14] * exp(fwd_beta[14] * tc[i] - activation_units[14] * fwd_Ea[14] * invT[i]);
        k_f_s[15*npt+i] = prefactor_units[15] * fwd_A[15] * exp(fwd_beta[15] * tc[i] - activation_units[15] * fwd_Ea[15] * invT[i]);
        k_f_s[16*npt+i] = prefactor_units[16] * fwd_A[16] * exp(fwd_beta[16] * tc[i] - activation_units[16] * fwd_Ea[16] * invT[i]);
        k_f_s[17*npt+i] = prefactor_units[17] * fwd_A[17] * exp(fwd_beta[17] * tc[i] - activation_units[17] * fwd_Ea[17] * invT[i]);
        k_f_s[18*npt+i] = prefactor_units[18] * fwd_A[18] * exp(fwd_beta[18] * tc[i] - activation_units[18] * fwd_Ea[18] * invT[i]);
        k_f_s[19*npt+i] = prefactor_units[19] * fwd_A[19] * exp(fwd_beta[19] * tc[i] - activation_units[19] * fwd_Ea[19] * invT[i]);
        k_f_s[20*npt+i] = prefactor_units[20] * fwd_A[20] * exp(fwd_beta[20] * tc[i] - activation_units[20] * fwd_Ea[20] * invT[i]);
        k_f_s[21*npt+i] = prefactor_units[21] * fwd_A[21] * exp(fwd_beta[21] * tc[i] - activation_units[21] * fwd_Ea[21] * invT[i]);
        k_f_s[22*npt+i] = prefactor_units[22] * fwd_A[22] * exp(fwd_beta[22] * tc[i] - activation_units[22] * fwd_Ea[22] * invT[i]);
        k_f_s[23*npt+i] = prefactor_units[23] * fwd_A[23] * exp(fwd_beta[23] * tc[i] - activation_units[23] * fwd_Ea[23] * invT[i]);
        k_f_s[24*npt+i] = prefactor_units[24] * fwd_A[24] * exp(fwd_beta[24] * tc[i] - activation_units[24] * fwd_Ea[24] * invT[i]);
        k_f_s[25*npt+i] = prefactor_units[25] * fwd_A[25] * exp(fwd_beta[25] * tc[i] - activation_units[25] * fwd_Ea[25] * invT[i]);
        k_f_s[26*npt+i] = prefactor_units[26] * fwd_A[26] * exp(fwd_beta[26] * tc[i] - activation_units[26] * fwd_Ea[26] * invT[i]);
        k_f_s[27*npt+i] = prefactor_units[27] * fwd_A[27] * exp(fwd_beta[27] * tc[i] - activation_units[27] * fwd_Ea[27] * invT[i]);
        k_f_s[28*npt+i] = prefactor_units[28] * fwd_A[28] * exp(fwd_beta[28] * tc[i] - activation_units[28] * fwd_Ea[28] * invT[i]);
    }
}

void vcomp_gibbs(int npt, double *  g_RT, double *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[12];
        tg[0] = tc[0*npt+i];
        tg[1] = tc[1*npt+i];
        tg[2] = tc[2*npt+i];
        tg[3] = tc[3*npt+i];
        tg[4] = tc[4*npt+i];

        gibbs(g, tg);

        g_RT[0*npt+i] = g[0];
        g_RT[1*npt+i] = g[1];
        g_RT[2*npt+i] = g[2];
        g_RT[3*npt+i] = g[3];
        g_RT[4*npt+i] = g[4];
        g_RT[5*npt+i] = g[5];
        g_RT[6*npt+i] = g[6];
        g_RT[7*npt+i] = g[7];
        g_RT[8*npt+i] = g[8];
        g_RT[9*npt+i] = g[9];
        g_RT[10*npt+i] = g[10];
        g_RT[11*npt+i] = g[11];
    }
}

void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = (101325. / 8.31451) * invT[i];
        double refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[1*npt+i] = refC * exp((g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[3*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[8*npt+i]) - (g_RT[9*npt+i]));
        Kc_s[3*npt+i] = refC * exp((g_RT[0*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[2*npt+i]) - (g_RT[1*npt+i]));
        Kc_s[5*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[5*npt+i]) - (g_RT[3*npt+i]));
        Kc_s[6*npt+i] = refCinv * exp((g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[4*npt+i]));
        Kc_s[7*npt+i] = refC * exp((g_RT[10*npt+i]) - (g_RT[5*npt+i] + g_RT[8*npt+i]));
        Kc_s[8*npt+i] = exp((g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[3*npt+i]));
        Kc_s[9*npt+i] = exp((g_RT[0*npt+i] + g_RT[2*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[10*npt+i] = exp((g_RT[0*npt+i] + g_RT[3*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[11*npt+i] = exp((g_RT[2*npt+i] + g_RT[4*npt+i]) - (g_RT[3*npt+i] + g_RT[3*npt+i]));
        Kc_s[12*npt+i] = exp((g_RT[5*npt+i] + g_RT[6*npt+i]) - (g_RT[0*npt+i] + g_RT[1*npt+i]));
        Kc_s[13*npt+i] = exp((g_RT[5*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[3*npt+i]));
        Kc_s[14*npt+i] = exp((g_RT[2*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[3*npt+i]));
        Kc_s[15*npt+i] = exp((g_RT[3*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[7*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[7*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[0*npt+i] + g_RT[6*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[6*npt+i]));
        Kc_s[21*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[1*npt+i] + g_RT[8*npt+i]) - (g_RT[2*npt+i] + g_RT[9*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[6*npt+i] + g_RT[8*npt+i]) - (g_RT[3*npt+i] + g_RT[9*npt+i]));
        Kc_s[25*npt+i] = exp((g_RT[3*npt+i] + g_RT[8*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[1*npt+i] + g_RT[10*npt+i]) - (g_RT[6*npt+i] + g_RT[8*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[5*npt+i] + g_RT[10*npt+i]) - (g_RT[0*npt+i] + g_RT[8*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[2*npt+i] + g_RT[10*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
    }
}

void vcomp_wdot(int npt, double *  wdot, double *  mixture, double *  sc,
		double *  k_f_s, double *  Kc_s,
		double *  tc, double *  invT, double *  T)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        double alpha;
        double redP, F;
        double logPred;
        double logFcent, troe_c, troe_n, troe, F_troe;

        /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
        phi_f = sc[1*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[0*npt+i] + (TB[0][1] - 1)*sc[4*npt+i] + (TB[0][2] - 1)*sc[1*npt+i] + (TB[0][3] - 1)*sc[8*npt+i] + (TB[0][4] - 1)*sc[9*npt+i];
        k_f = k_f_s[0*npt+i];
        redP = alpha / k_f * phase_units[0] * low_A[0] * exp(low_beta[0] * tc[i] - activation_units[0] * low_Ea[0] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[0]) > 1.e-100 ? (1.-troe_a[0])*exp(-T[i]/troe_Tsss[0]) : 0.) 
            + (fabs(troe_Ts[0]) > 1.e-100 ? troe_a[0] * exp(-T[i]/troe_Ts[0]) : 0.) 
            + (troe_len[0] == 4 ? exp(-troe_Tss[0] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
        phi_f = sc[7*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[0*npt+i] + (TB[1][1] - 1)*sc[4*npt+i] + (TB[1][2] - 1)*sc[8*npt+i] + (TB[1][3] - 1)*sc[9*npt+i];
        k_f = k_f_s[1*npt+i];
        redP = alpha / k_f * phase_units[1] * low_A[1] * exp(low_beta[1] * tc[i] - activation_units[1] * low_Ea[1] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[1]) > 1.e-100 ? (1.-troe_a[1])*exp(-T[i]/troe_Tsss[1]) : 0.) 
            + (fabs(troe_Ts[1]) > 1.e-100 ? troe_a[1] * exp(-T[i]/troe_Ts[1]) : 0.) 
            + (troe_len[1] == 4 ? exp(-troe_Tss[1] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[3*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 3: CO + O (+M) <=> CO2 (+M) */
        phi_f = sc[2*npt+i]*sc[8*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[0*npt+i] + (TB[2][1] - 1)*sc[4*npt+i] + (TB[2][2] - 1)*sc[8*npt+i] + (TB[2][3] - 1)*sc[9*npt+i];
        k_f = k_f_s[2*npt+i];
        redP = alpha / k_f * phase_units[2] * low_A[2] * exp(low_beta[2] * tc[i] - activation_units[2] * low_Ea[2] * invT[i]);
        F = redP / (1 + redP);
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 4: H2 + M <=> H + H + M */
        phi_f = sc[0*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[0*npt+i] + (TB[3][1] - 1)*sc[4*npt+i] + (TB[3][2] - 1)*sc[8*npt+i] + (TB[3][3] - 1)*sc[9*npt+i];
        k_f = alpha * k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 5: O + O + M <=> O2 + M */
        phi_f = sc[2*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[0*npt+i] + (TB[4][1] - 1)*sc[4*npt+i] + (TB[4][2] - 1)*sc[8*npt+i] + (TB[4][3] - 1)*sc[9*npt+i];
        k_f = alpha * k_f_s[4*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;

        /*reaction 6: O + H + M <=> OH + M */
        phi_f = sc[2*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[0*npt+i] + (TB[5][1] - 1)*sc[4*npt+i] + (TB[5][2] - 1)*sc[8*npt+i] + (TB[5][3] - 1)*sc[9*npt+i];
        k_f = alpha * k_f_s[5*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 7: H + OH + M <=> H2O + M */
        phi_f = sc[3*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[6][0] - 1)*sc[0*npt+i] + (TB[6][1] - 1)*sc[4*npt+i] + (TB[6][2] - 1)*sc[8*npt+i] + (TB[6][3] - 1)*sc[9*npt+i];
        k_f = alpha * k_f_s[6*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 8: HCO + M <=> H + CO + M */
        phi_f = sc[10*npt+i];
        alpha = mixture[i] + (TB[7][0] - 1)*sc[0*npt+i] + (TB[7][1] - 1)*sc[4*npt+i] + (TB[7][2] - 1)*sc[8*npt+i] + (TB[7][3] - 1)*sc[9*npt+i];
        k_f = alpha * k_f_s[7*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[8*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 9: H + O2 <=> O + OH */
        phi_f = sc[1*npt+i]*sc[5*npt+i];
        k_f = k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[3*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 10: O + H2 <=> H + OH */
        phi_f = sc[0*npt+i]*sc[2*npt+i];
        k_f = k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 11: H2 + OH <=> H2O + H */
        phi_f = sc[0*npt+i]*sc[3*npt+i];
        k_f = k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[5*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 12: O + H2O <=> OH + OH */
        phi_f = sc[2*npt+i]*sc[4*npt+i];
        k_f = k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[3*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;

        /*reaction 13: HO2 + H <=> H2 + O2 */
        phi_f = sc[5*npt+i]*sc[6*npt+i];
        k_f = k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[1*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 14: HO2 + H <=> OH + OH */
        phi_f = sc[5*npt+i]*sc[6*npt+i];
        k_f = k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[3*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 15: HO2 + O <=> O2 + OH */
        phi_f = sc[2*npt+i]*sc[6*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[3*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 16: HO2 + OH <=> H2O + O2 */
        phi_f = sc[3*npt+i]*sc[6*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 17: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[7*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 18: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[7*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 19: H2O2 + H <=> H2O + OH */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[4*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 20: H2O2 + H <=> HO2 + H2 */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[6*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 21: H2O2 + O <=> OH + HO2 */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[6*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 22: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 23: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 24: CO + O2 <=> CO2 + O */
        phi_f = sc[1*npt+i]*sc[8*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[9*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 25: CO + HO2 <=> CO2 + OH */
        phi_f = sc[6*npt+i]*sc[8*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[9*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 26: CO + OH <=> CO2 + H */
        phi_f = sc[3*npt+i]*sc[8*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[9*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 27: HCO + O2 <=> CO + HO2 */
        phi_f = sc[1*npt+i]*sc[10*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[8*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 28: HCO + H <=> CO + H2 */
        phi_f = sc[5*npt+i]*sc[10*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[8*npt+i];
        Kc = Kc_s[27*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 29: HCO + O <=> CO2 + H */
        phi_f = sc[2*npt+i]*sc[10*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[9*npt+i];
        Kc = Kc_s[28*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
    }
}

/*compute an approx to the reaction Jacobian */
void DWDOT_PRECOND(double *  J, double *  sc, double *  Tp, int * HP)
{
    double c[12];

    for (int k=0; k<12; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<12; k++) {
        J[156+k] *= 1.e-6;
        J[k*13+12] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
void DWDOT(double *  J, double *  sc, double *  Tp, int * consP)
{
    double c[12];

    for (int k=0; k<12; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<12; k++) {
        J[156+k] *= 1.e-6;
        J[k*13+12] *= 1.e6;
    }

    return;
}

/*compute the sparsity pattern Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    double c[12];
    double J[169];

    for (int k=0; k<12; k++) {
        c[k] = 1.0/ 12.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<13; k++) {
        for (int l=0; l<13; l++) {
            if(J[ 13 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}



/*compute the sparsity pattern of simplified Jacobian */
void SPARSITY_INFO_PRECOND( int * nJdata, int * consP)
{
    double c[12];
    double J[169];

    for (int k=0; k<12; k++) {
        c[k] = 1.0/ 12.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<13; k++) {
        for (int l=0; l<13; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J[ 13 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    nJdata[0] = nJdata_tmp;

    return;
}


/*compute the sparsity pattern of the simplified precond Jacobian */
void SPARSITY_PREPROC_PRECOND(int * rowVals, int * colPtrs, int * consP)
{
    double c[12];
    double J[169];

    for (int k=0; k<12; k++) {
        c[k] = 1.0/ 12.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<13; k++) {
        for (int l=0; l<13; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J[13*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        colPtrs[k+1] = nJdata_tmp;
    }

    return;
}
/*compute the sparsity pattern of the Jacobian */
void SPARSITY_PREPROC(int *  rowVals, int *  colPtrs, int * consP, int NCELLS)
{
    double c[12];
    double J[169];
    int offset_row;
    int offset_col;

    for (int k=0; k<12; k++) {
        c[k] = 1.0/ 12.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        offset_row = nc * 13;
        offset_col = nc * 13;
        for (int k=0; k<13; k++) {
            for (int l=0; l<13; l++) {
                if(J[13*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l + offset_row; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
            colPtrs[offset_col + (k + 1)] = nJdata_tmp;
        }
    }

    return;
}

/*compute the reaction Jacobian */
void aJacobian(double *  J, double *  sc, double T, int consP)
{
    for (int i=0; i<169; i++) {
        J[i] = 0.0;
    }

    double wdot[12];
    for (int k=0; k<12; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 12; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[12];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[12];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[12];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[4] + (TB[0][2] - 1)*sc[1] + (TB[0][3] - 1)*sc[8] + (TB[0][4] - 1)*sc[9];
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = prefactor_units[0] * fwd_A[0]
                * exp(fwd_beta[0] * tc[0] - activation_units[0] * fwd_Ea[0] * invT);
    dlnkfdT = fwd_beta[0] * invT + activation_units[0] * fwd_Ea[0] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[0] * exp(low_beta[0] * tc[0] - activation_units[0] * low_Ea[0] * invT);
    Pr = phase_units[0] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[0] * invT + activation_units[0] * low_Ea[0] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[0]) > 1.e-100 ? (1.-troe_a[0])*exp(-T/troe_Tsss[0]) : 0.);
    Fcent2 = (fabs(troe_Ts[0]) > 1.e-100 ? troe_a[0] * exp(-T/troe_Ts[0]) : 0.);
    Fcent3 = (troe_len[0] == 4 ? exp(-troe_Tss[0] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[0]) > 1.e-100 ? -Fcent1/troe_Tsss[0] : 0.)
      + (fabs(troe_Ts[0]) > 1.e-100 ? -Fcent2/troe_Ts[0] : 0.)
      + (troe_len[0] == 4 ? Fcent3*troe_Tss[0]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[5] -= q; /* H */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[O2]/d[H2] */
        J[5] -= dqdci;                /* dwdot[H]/d[H2] */
        J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
        /* d()/d[O2] */
        dqdci = (TB[0][2] - 1)*dcdc_fac + k_f*sc[5];
        J[14] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[18] -= dqdci;               /* dwdot[H]/d[O2] */
        J[19] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[53] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[57] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[58] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[1];
        J[66] -= dqdci;               /* dwdot[O2]/d[H] */
        J[70] -= dqdci;               /* dwdot[H]/d[H] */
        J[71] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[79] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[83] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[84] += dqdci;               /* dwdot[HO2]/d[HO2] */
        /* d()/d[CO] */
        dqdci = (TB[0][3] - 1)*dcdc_fac;
        J[105] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[109] -= dqdci;              /* dwdot[H]/d[CO] */
        J[110] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[0][4] - 1)*dcdc_fac;
        J[118] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[122] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[123] += dqdci;              /* dwdot[HO2]/d[CO2] */
    }
    else {
        dqdc[0] = TB[0][0]*dcdc_fac;
        dqdc[1] = TB[0][2]*dcdc_fac + k_f*sc[5];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[0][1]*dcdc_fac;
        dqdc[5] = dcdc_fac + k_f*sc[1];
        dqdc[6] = dcdc_fac - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[0][3]*dcdc_fac;
        dqdc[9] = TB[0][4]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        for (int k=0; k<12; k++) {
            J[13*k+1] -= dqdc[k];
            J[13*k+5] -= dqdc[k];
            J[13*k+6] += dqdc[k];
        }
    }
    J[157] -= dqdT; /* dwdot[O2]/dT */
    J[161] -= dqdT; /* dwdot[H]/dT */
    J[162] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[4] + (TB[1][2] - 1)*sc[8] + (TB[1][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[7];
    k_f = prefactor_units[1] * fwd_A[1]
                * exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
    dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[1] * exp(low_beta[1] * tc[0] - activation_units[1] * low_Ea[1] * invT);
    Pr = phase_units[1] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[1] * invT + activation_units[1] * low_Ea[1] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[1]) > 1.e-100 ? (1.-troe_a[1])*exp(-T/troe_Tsss[1]) : 0.);
    Fcent2 = (fabs(troe_Ts[1]) > 1.e-100 ? troe_a[1] * exp(-T/troe_Ts[1]) : 0.);
    Fcent3 = (troe_len[1] == 4 ? exp(-troe_Tss[1] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[1]) > 1.e-100 ? -Fcent1/troe_Tsss[1] : 0.)
      + (fabs(troe_Ts[1]) > 1.e-100 ? -Fcent2/troe_Ts[1] : 0.)
      + (troe_len[1] == 4 ? Fcent3*troe_Tss[1]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[3]*sc[3];
    Kc = refC * exp(-g_RT[3] - g_RT[3] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2*h_RT[3]) - 1);
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
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[3] += 2 * dqdci;            /* dwdot[OH]/d[H2] */
        J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  - k_r*2*sc[3];
        J[42] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
        J[46] -= dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[1][1] - 1)*dcdc_fac;
        J[55] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
        J[59] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[H2O2] */
        dqdci =  + k_f;
        J[94] += 2 * dqdci;           /* dwdot[OH]/d[H2O2] */
        J[98] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[CO] */
        dqdci = (TB[1][2] - 1)*dcdc_fac;
        J[107] += 2 * dqdci;          /* dwdot[OH]/d[CO] */
        J[111] -= dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[1][3] - 1)*dcdc_fac;
        J[120] += 2 * dqdci;          /* dwdot[OH]/d[CO2] */
        J[124] -= dqdci;              /* dwdot[H2O2]/d[CO2] */
    }
    else {
        dqdc[0] = TB[1][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac - k_r*2*sc[3];
        dqdc[4] = TB[1][1]*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac + k_f;
        dqdc[8] = TB[1][2]*dcdc_fac;
        dqdc[9] = TB[1][3]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        for (int k=0; k<12; k++) {
            J[13*k+3] += 2 * dqdc[k];
            J[13*k+7] -= dqdc[k];
        }
    }
    J[159] += 2 * dqdT; /* dwdot[OH]/dT */
    J[163] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[4] + (TB[2][2] - 1)*sc[8] + (TB[2][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[2] * exp(low_beta[2] * tc[0] - activation_units[2] * low_Ea[2] * invT);
    Pr = phase_units[2] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[2] * invT + activation_units[2] * low_Ea[2] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Lindemann form */
    F = 1.0;
    dlogFdlogPr = 0.0;
    dlogFdT = 0.0;
    /* reverse */
    phi_r = sc[9];
    Kc = refCinv * exp(g_RT[2] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[9]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[2][0] - 1)*dcdc_fac;
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[8] -= dqdci;                /* dwdot[CO]/d[H2] */
        J[9] += dqdci;                /* dwdot[CO2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[8];
        J[28] -= dqdci;               /* dwdot[O]/d[O] */
        J[34] -= dqdci;               /* dwdot[CO]/d[O] */
        J[35] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*dcdc_fac;
        J[54] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[60] -= dqdci;               /* dwdot[CO]/d[H2O] */
        J[61] += dqdci;               /* dwdot[CO2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[2][2] - 1)*dcdc_fac + k_f*sc[2];
        J[106] -= dqdci;              /* dwdot[O]/d[CO] */
        J[112] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[113] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[2][3] - 1)*dcdc_fac - k_r;
        J[119] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[125] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[126] += dqdci;              /* dwdot[CO2]/d[CO2] */
    }
    else {
        dqdc[0] = TB[2][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[8];
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[2][1]*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[2][2]*dcdc_fac + k_f*sc[2];
        dqdc[9] = TB[2][3]*dcdc_fac - k_r;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        for (int k=0; k<12; k++) {
            J[13*k+2] -= dqdc[k];
            J[13*k+8] -= dqdc[k];
            J[13*k+9] += dqdc[k];
        }
    }
    J[158] -= dqdT; /* dwdot[O]/dT */
    J[164] -= dqdT; /* dwdot[CO]/dT */
    J[165] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 4: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[4] + (TB[3][2] - 1)*sc[8] + (TB[3][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[0];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[5];
    Kc = refC * exp(g_RT[0] - g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0]) + (2*h_RT[5]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[5] += 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*q_nocor + k_f;
        J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
        J[5] += 2 * dqdci;            /* dwdot[H]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*q_nocor;
        J[52] -= dqdci;               /* dwdot[H2]/d[H2O] */
        J[57] += 2 * dqdci;           /* dwdot[H]/d[H2O] */
        /* d()/d[H] */
        dqdci =  - k_r*2*sc[5];
        J[65] -= dqdci;               /* dwdot[H2]/d[H] */
        J[70] += 2 * dqdci;           /* dwdot[H]/d[H] */
        /* d()/d[CO] */
        dqdci = (TB[3][2] - 1)*q_nocor;
        J[104] -= dqdci;              /* dwdot[H2]/d[CO] */
        J[109] += 2 * dqdci;          /* dwdot[H]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[3][3] - 1)*q_nocor;
        J[117] -= dqdci;              /* dwdot[H2]/d[CO2] */
        J[122] += 2 * dqdci;          /* dwdot[H]/d[CO2] */
    }
    else {
        dqdc[0] = TB[3][0]*q_nocor + k_f;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = TB[3][1]*q_nocor;
        dqdc[5] = q_nocor - k_r*2*sc[5];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[3][2]*q_nocor;
        dqdc[9] = TB[3][3]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        for (int k=0; k<12; k++) {
            J[13*k+0] -= dqdc[k];
            J[13*k+5] += 2 * dqdc[k];
        }
    }
    J[156] -= dqdT; /* dwdot[H2]/dT */
    J[161] += 2 * dqdT; /* dwdot[H]/dT */

    /*reaction 5: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[4] + (TB[4][2] - 1)*sc[8] + (TB[4][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[2]*sc[2];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[1];
    Kc = refCinv * exp(-g_RT[1] + g_RT[2] + g_RT[2]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[2]) + (h_RT[1]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[2] -= 2 * q; /* O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[4][0] - 1)*q_nocor;
        J[1] += dqdci;                /* dwdot[O2]/d[H2] */
        J[2] += -2 * dqdci;           /* dwdot[O]/d[H2] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[14] += dqdci;               /* dwdot[O2]/d[O2] */
        J[15] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        /* d()/d[O] */
        dqdci =  + k_f*2*sc[2];
        J[27] += dqdci;               /* dwdot[O2]/d[O] */
        J[28] += -2 * dqdci;          /* dwdot[O]/d[O] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*q_nocor;
        J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
        J[54] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[4][2] - 1)*q_nocor;
        J[105] += dqdci;              /* dwdot[O2]/d[CO] */
        J[106] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[4][3] - 1)*q_nocor;
        J[118] += dqdci;              /* dwdot[O2]/d[CO2] */
        J[119] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
    }
    else {
        dqdc[0] = TB[4][0]*q_nocor;
        dqdc[1] = q_nocor - k_r;
        dqdc[2] = q_nocor + k_f*2*sc[2];
        dqdc[3] = q_nocor;
        dqdc[4] = TB[4][1]*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[4][2]*q_nocor;
        dqdc[9] = TB[4][3]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        for (int k=0; k<12; k++) {
            J[13*k+1] += dqdc[k];
            J[13*k+2] += -2 * dqdc[k];
        }
    }
    J[157] += dqdT; /* dwdot[O2]/dT */
    J[158] += -2 * dqdT; /* dwdot[O]/dT */

    /*reaction 6: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[4] + (TB[5][2] - 1)*sc[8] + (TB[5][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[2]*sc[5];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[3]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[5][0] - 1)*q_nocor;
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[3] += dqdci;                /* dwdot[OH]/d[H2] */
        J[5] -= dqdci;                /* dwdot[H]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[5];
        J[28] -= dqdci;               /* dwdot[O]/d[O] */
        J[29] += dqdci;               /* dwdot[OH]/d[O] */
        J[31] -= dqdci;               /* dwdot[H]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[41] -= dqdci;               /* dwdot[O]/d[OH] */
        J[42] += dqdci;               /* dwdot[OH]/d[OH] */
        J[44] -= dqdci;               /* dwdot[H]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*q_nocor;
        J[54] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[55] += dqdci;               /* dwdot[OH]/d[H2O] */
        J[57] -= dqdci;               /* dwdot[H]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[67] -= dqdci;               /* dwdot[O]/d[H] */
        J[68] += dqdci;               /* dwdot[OH]/d[H] */
        J[70] -= dqdci;               /* dwdot[H]/d[H] */
        /* d()/d[CO] */
        dqdci = (TB[5][2] - 1)*q_nocor;
        J[106] -= dqdci;              /* dwdot[O]/d[CO] */
        J[107] += dqdci;              /* dwdot[OH]/d[CO] */
        J[109] -= dqdci;              /* dwdot[H]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[5][3] - 1)*q_nocor;
        J[119] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[120] += dqdci;              /* dwdot[OH]/d[CO2] */
        J[122] -= dqdci;              /* dwdot[H]/d[CO2] */
    }
    else {
        dqdc[0] = TB[5][0]*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*sc[5];
        dqdc[3] = q_nocor - k_r;
        dqdc[4] = TB[5][1]*q_nocor;
        dqdc[5] = q_nocor + k_f*sc[2];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[5][2]*q_nocor;
        dqdc[9] = TB[5][3]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        for (int k=0; k<12; k++) {
            J[13*k+2] -= dqdc[k];
            J[13*k+3] += dqdc[k];
            J[13*k+5] -= dqdc[k];
        }
    }
    J[158] -= dqdT; /* dwdot[O]/dT */
    J[159] += dqdT; /* dwdot[OH]/dT */
    J[161] -= dqdT; /* dwdot[H]/dT */

    /*reaction 7: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[4] + (TB[6][2] - 1)*sc[8] + (TB[6][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[3]*sc[5];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[3] - g_RT[4] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[4]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] -= q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[6][0] - 1)*q_nocor;
        J[3] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[4] += dqdci;                /* dwdot[H2O]/d[H2] */
        J[5] -= dqdci;                /* dwdot[H]/d[H2] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[5];
        J[42] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[43] += dqdci;               /* dwdot[H2O]/d[OH] */
        J[44] -= dqdci;               /* dwdot[H]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[6][1] - 1)*q_nocor - k_r;
        J[55] -= dqdci;               /* dwdot[OH]/d[H2O] */
        J[56] += dqdci;               /* dwdot[H2O]/d[H2O] */
        J[57] -= dqdci;               /* dwdot[H]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[3];
        J[68] -= dqdci;               /* dwdot[OH]/d[H] */
        J[69] += dqdci;               /* dwdot[H2O]/d[H] */
        J[70] -= dqdci;               /* dwdot[H]/d[H] */
        /* d()/d[CO] */
        dqdci = (TB[6][2] - 1)*q_nocor;
        J[107] -= dqdci;              /* dwdot[OH]/d[CO] */
        J[108] += dqdci;              /* dwdot[H2O]/d[CO] */
        J[109] -= dqdci;              /* dwdot[H]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[6][3] - 1)*q_nocor;
        J[120] -= dqdci;              /* dwdot[OH]/d[CO2] */
        J[121] += dqdci;              /* dwdot[H2O]/d[CO2] */
        J[122] -= dqdci;              /* dwdot[H]/d[CO2] */
    }
    else {
        dqdc[0] = TB[6][0]*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor + k_f*sc[5];
        dqdc[4] = TB[6][1]*q_nocor - k_r;
        dqdc[5] = q_nocor + k_f*sc[3];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[6][2]*q_nocor;
        dqdc[9] = TB[6][3]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        for (int k=0; k<12; k++) {
            J[13*k+3] -= dqdc[k];
            J[13*k+4] += dqdc[k];
            J[13*k+5] -= dqdc[k];
        }
    }
    J[159] -= dqdT; /* dwdot[OH]/dT */
    J[160] += dqdT; /* dwdot[H2O]/dT */
    J[161] -= dqdT; /* dwdot[H]/dT */

    /*reaction 8: HCO + M <=> H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[7][0] - 1)*sc[0] + (TB[7][1] - 1)*sc[4] + (TB[7][2] - 1)*sc[8] + (TB[7][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[10];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[8];
    Kc = refC * exp(-g_RT[5] - g_RT[8] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10]) + (h_RT[5] + h_RT[8]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* H */
    wdot[8] += q; /* CO */
    wdot[10] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[7][0] - 1)*q_nocor;
        J[5] += dqdci;                /* dwdot[H]/d[H2] */
        J[8] += dqdci;                /* dwdot[CO]/d[H2] */
        J[10] -= dqdci;               /* dwdot[HCO]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[7][1] - 1)*q_nocor;
        J[57] += dqdci;               /* dwdot[H]/d[H2O] */
        J[60] += dqdci;               /* dwdot[CO]/d[H2O] */
        J[62] -= dqdci;               /* dwdot[HCO]/d[H2O] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[8];
        J[70] += dqdci;               /* dwdot[H]/d[H] */
        J[73] += dqdci;               /* dwdot[CO]/d[H] */
        J[75] -= dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[CO] */
        dqdci = (TB[7][2] - 1)*q_nocor - k_r*sc[5];
        J[109] += dqdci;              /* dwdot[H]/d[CO] */
        J[112] += dqdci;              /* dwdot[CO]/d[CO] */
        J[114] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[7][3] - 1)*q_nocor;
        J[122] += dqdci;              /* dwdot[H]/d[CO2] */
        J[125] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[127] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[135] += dqdci;              /* dwdot[H]/d[HCO] */
        J[138] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[140] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    }
    else {
        dqdc[0] = TB[7][0]*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = TB[7][1]*q_nocor;
        dqdc[5] = q_nocor - k_r*sc[8];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[7][2]*q_nocor - k_r*sc[5];
        dqdc[9] = TB[7][3]*q_nocor;
        dqdc[10] = q_nocor + k_f;
        dqdc[11] = q_nocor;
        for (int k=0; k<12; k++) {
            J[13*k+5] += dqdc[k];
            J[13*k+8] += dqdc[k];
            J[13*k+10] -= dqdc[k];
        }
    }
    J[161] += dqdT; /* dwdot[H]/dT */
    J[164] += dqdT; /* dwdot[CO]/dT */
    J[166] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 9: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[3];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[2] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[2] += q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* H */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[5];
    J[14] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[15] += dqdci;               /* dwdot[O]/d[O2] */
    J[16] += dqdci;               /* dwdot[OH]/d[O2] */
    J[18] -= dqdci;               /* dwdot[H]/d[O2] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[3];
    J[27] -= dqdci;               /* dwdot[O2]/d[O] */
    J[28] += dqdci;               /* dwdot[O]/d[O] */
    J[29] += dqdci;               /* dwdot[OH]/d[O] */
    J[31] -= dqdci;               /* dwdot[H]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[40] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[41] += dqdci;               /* dwdot[O]/d[OH] */
    J[42] += dqdci;               /* dwdot[OH]/d[OH] */
    J[44] -= dqdci;               /* dwdot[H]/d[OH] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[1];
    J[66] -= dqdci;               /* dwdot[O2]/d[H] */
    J[67] += dqdci;               /* dwdot[O]/d[H] */
    J[68] += dqdci;               /* dwdot[OH]/d[H] */
    J[70] -= dqdci;               /* dwdot[H]/d[H] */
    /* d()/dT */
    J[157] -= dqdT;               /* dwdot[O2]/dT */
    J[158] += dqdT;               /* dwdot[O]/dT */
    J[159] += dqdT;               /* dwdot[OH]/dT */
    J[161] -= dqdT;               /* dwdot[H]/dT */

    /*reaction 10: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[0] + g_RT[2] - g_RT[3] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] += q; /* H */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[2] -= dqdci;                /* dwdot[O]/d[H2] */
    J[3] += dqdci;                /* dwdot[OH]/d[H2] */
    J[5] += dqdci;                /* dwdot[H]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[26] -= dqdci;               /* dwdot[H2]/d[O] */
    J[28] -= dqdci;               /* dwdot[O]/d[O] */
    J[29] += dqdci;               /* dwdot[OH]/d[O] */
    J[31] += dqdci;               /* dwdot[H]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[39] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[41] -= dqdci;               /* dwdot[O]/d[OH] */
    J[42] += dqdci;               /* dwdot[OH]/d[OH] */
    J[44] += dqdci;               /* dwdot[H]/d[OH] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[65] -= dqdci;               /* dwdot[H2]/d[H] */
    J[67] -= dqdci;               /* dwdot[O]/d[H] */
    J[68] += dqdci;               /* dwdot[OH]/d[H] */
    J[70] += dqdci;               /* dwdot[H]/d[H] */
    /* d()/dT */
    J[156] -= dqdT;               /* dwdot[H2]/dT */
    J[158] -= dqdT;               /* dwdot[O]/dT */
    J[159] += dqdT;               /* dwdot[OH]/dT */
    J[161] += dqdT;               /* dwdot[H]/dT */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[3];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[0] + g_RT[3] - g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[3]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] += q; /* H */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[3] -= dqdci;                /* dwdot[OH]/d[H2] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H2] */
    J[5] += dqdci;                /* dwdot[H]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[39] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[42] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[43] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[44] += dqdci;               /* dwdot[H]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[52] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[55] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[56] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[57] += dqdci;               /* dwdot[H]/d[H2O] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[65] -= dqdci;               /* dwdot[H2]/d[H] */
    J[68] -= dqdci;               /* dwdot[OH]/d[H] */
    J[69] += dqdci;               /* dwdot[H2O]/d[H] */
    J[70] += dqdci;               /* dwdot[H]/d[H] */
    /* d()/dT */
    J[156] -= dqdT;               /* dwdot[H2]/dT */
    J[159] -= dqdT;               /* dwdot[OH]/dT */
    J[160] += dqdT;               /* dwdot[H2O]/dT */
    J[161] += dqdT;               /* dwdot[H]/dT */

    /*reaction 12: O + H2O <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[3];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[3] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2*h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += 2 * q; /* OH */
    wdot[4] -= q; /* H2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[4];
    J[28] -= dqdci;               /* dwdot[O]/d[O] */
    J[29] += 2 * dqdci;           /* dwdot[OH]/d[O] */
    J[30] -= dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[3];
    J[41] -= dqdci;               /* dwdot[O]/d[OH] */
    J[42] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[2];
    J[54] -= dqdci;               /* dwdot[O]/d[H2O] */
    J[55] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
    J[56] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[158] -= dqdT;               /* dwdot[O]/dT */
    J[159] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[160] -= dqdT;               /* dwdot[H2O]/dT */

    /*reaction 13: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[1];
    Kc = exp(-g_RT[0] - g_RT[1] + g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (h_RT[0] + h_RT[1]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] += q; /* O2 */
    wdot[5] -= q; /* H */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[1];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[O2]/d[H2] */
    J[5] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[13] += dqdci;               /* dwdot[H2]/d[O2] */
    J[14] += dqdci;               /* dwdot[O2]/d[O2] */
    J[18] -= dqdci;               /* dwdot[H]/d[O2] */
    J[19] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[65] += dqdci;               /* dwdot[H2]/d[H] */
    J[66] += dqdci;               /* dwdot[O2]/d[H] */
    J[70] -= dqdci;               /* dwdot[H]/d[H] */
    J[71] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[78] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[79] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[83] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[84] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[156] += dqdT;               /* dwdot[H2]/dT */
    J[157] += dqdT;               /* dwdot[O2]/dT */
    J[161] -= dqdT;               /* dwdot[H]/dT */
    J[162] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 14: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[3];
    Kc = exp(-g_RT[3] - g_RT[3] + g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (2*h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += 2 * q; /* OH */
    wdot[5] -= q; /* H */
    wdot[6] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[3];
    J[42] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[44] -= dqdci;               /* dwdot[H]/d[OH] */
    J[45] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[68] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[70] -= dqdci;               /* dwdot[H]/d[H] */
    J[71] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[81] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[83] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[84] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[159] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[161] -= dqdT;               /* dwdot[H]/dT */
    J[162] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 15: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[3];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[1] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[3];
    J[14] += dqdci;               /* dwdot[O2]/d[O2] */
    J[15] -= dqdci;               /* dwdot[O]/d[O2] */
    J[16] += dqdci;               /* dwdot[OH]/d[O2] */
    J[19] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[27] += dqdci;               /* dwdot[O2]/d[O] */
    J[28] -= dqdci;               /* dwdot[O]/d[O] */
    J[29] += dqdci;               /* dwdot[OH]/d[O] */
    J[32] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[40] += dqdci;               /* dwdot[O2]/d[OH] */
    J[41] -= dqdci;               /* dwdot[O]/d[OH] */
    J[42] += dqdci;               /* dwdot[OH]/d[OH] */
    J[45] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[79] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[80] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[81] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[84] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[157] += dqdT;               /* dwdot[O2]/dT */
    J[158] -= dqdT;               /* dwdot[O]/dT */
    J[159] += dqdT;               /* dwdot[OH]/dT */
    J[162] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4];
    Kc = exp(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[1] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[14] += dqdci;               /* dwdot[O2]/d[O2] */
    J[16] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[17] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[19] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[40] += dqdci;               /* dwdot[O2]/d[OH] */
    J[42] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[43] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[45] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[55] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[56] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[58] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[79] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[81] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[82] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[84] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[157] += dqdT;               /* dwdot[O2]/dT */
    J[159] -= dqdT;               /* dwdot[OH]/dT */
    J[160] += dqdT;               /* dwdot[H2O]/dT */
    J[162] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 17: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[1] + h_RT[7]));
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
    J[14] += dqdci;               /* dwdot[O2]/d[O2] */
    J[19] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[20] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[6];
    J[79] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[84] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[92] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[97] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[98] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[157] += dqdT;               /* dwdot[O2]/dT */
    J[162] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[163] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 18: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[1] + h_RT[7]));
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
    J[14] += dqdci;               /* dwdot[O2]/d[O2] */
    J[19] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[20] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[6];
    J[79] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[84] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[92] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[97] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[98] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[157] += dqdT;               /* dwdot[O2]/dT */
    J[162] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[163] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 19: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(-g_RT[3] - g_RT[4] + g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] -= q; /* H */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[42] += dqdci;               /* dwdot[OH]/d[OH] */
    J[43] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[44] -= dqdci;               /* dwdot[H]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[55] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[56] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[59] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[68] += dqdci;               /* dwdot[OH]/d[H] */
    J[69] += dqdci;               /* dwdot[H2O]/d[H] */
    J[70] -= dqdci;               /* dwdot[H]/d[H] */
    J[72] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[94] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[95] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[96] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[98] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[159] += dqdT;               /* dwdot[OH]/dT */
    J[160] += dqdT;               /* dwdot[H2O]/dT */
    J[161] -= dqdT;               /* dwdot[H]/dT */
    J[163] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 20: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[5] -= q; /* H */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[5] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[65] += dqdci;               /* dwdot[H2]/d[H] */
    J[70] -= dqdci;               /* dwdot[H]/d[H] */
    J[71] += dqdci;               /* dwdot[HO2]/d[H] */
    J[72] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[78] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[83] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[84] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[91] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[96] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[97] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[98] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[156] += dqdT;               /* dwdot[H2]/dT */
    J[161] -= dqdT;               /* dwdot[H]/dT */
    J[162] += dqdT;               /* dwdot[HO2]/dT */
    J[163] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 21: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
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
    J[28] -= dqdci;               /* dwdot[O]/d[O] */
    J[29] += dqdci;               /* dwdot[OH]/d[O] */
    J[32] += dqdci;               /* dwdot[HO2]/d[O] */
    J[33] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[41] -= dqdci;               /* dwdot[O]/d[OH] */
    J[42] += dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[80] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[81] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[84] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[93] -= dqdci;               /* dwdot[O]/d[H2O2] */
    J[94] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[97] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[98] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[158] -= dqdT;               /* dwdot[O]/dT */
    J[159] += dqdT;               /* dwdot[OH]/dT */
    J[162] += dqdT;               /* dwdot[HO2]/dT */
    J[163] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 22: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
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
    J[42] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[43] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[45] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[55] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[56] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[58] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[59] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[81] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[82] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[84] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[94] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[95] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[97] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[98] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[159] -= dqdT;               /* dwdot[OH]/dT */
    J[160] += dqdT;               /* dwdot[H2O]/dT */
    J[162] += dqdT;               /* dwdot[HO2]/dT */
    J[163] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 23: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
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
    J[42] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[43] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[45] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[55] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[56] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[58] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[59] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[81] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[82] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[84] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[94] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[95] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[97] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[98] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[159] -= dqdT;               /* dwdot[OH]/dT */
    J[160] += dqdT;               /* dwdot[H2O]/dT */
    J[162] += dqdT;               /* dwdot[HO2]/dT */
    J[163] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 24: CO + O2 <=> CO2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[8];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[9];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[8]) + (h_RT[2] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[2] += q; /* O */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[14] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[15] += dqdci;               /* dwdot[O]/d[O2] */
    J[21] -= dqdci;               /* dwdot[CO]/d[O2] */
    J[22] += dqdci;               /* dwdot[CO2]/d[O2] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[9];
    J[27] -= dqdci;               /* dwdot[O2]/d[O] */
    J[28] += dqdci;               /* dwdot[O]/d[O] */
    J[34] -= dqdci;               /* dwdot[CO]/d[O] */
    J[35] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[1];
    J[105] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[106] += dqdci;              /* dwdot[O]/d[CO] */
    J[112] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[113] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[2];
    J[118] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[119] += dqdci;              /* dwdot[O]/d[CO2] */
    J[125] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[126] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[157] -= dqdT;               /* dwdot[O2]/dT */
    J[158] += dqdT;               /* dwdot[O]/dT */
    J[164] -= dqdT;               /* dwdot[CO]/dT */
    J[165] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 25: CO + HO2 <=> CO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[8];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[9];
    Kc = exp(-g_RT[3] + g_RT[6] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[3] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[42] += dqdci;               /* dwdot[OH]/d[OH] */
    J[45] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[48] += dqdci;               /* dwdot[CO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[8];
    J[81] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[84] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    J[86] -= dqdci;               /* dwdot[CO]/d[HO2] */
    J[87] += dqdci;               /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[107] += dqdci;              /* dwdot[OH]/d[CO] */
    J[110] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[112] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[113] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[3];
    J[120] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[123] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[125] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[126] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[159] += dqdT;               /* dwdot[OH]/dT */
    J[162] -= dqdT;               /* dwdot[HO2]/dT */
    J[164] -= dqdT;               /* dwdot[CO]/dT */
    J[165] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 26: CO + OH <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(g_RT[3] - g_RT[5] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[5] += q; /* H */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[42] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[44] += dqdci;               /* dwdot[H]/d[OH] */
    J[47] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[48] += dqdci;               /* dwdot[CO2]/d[OH] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[68] -= dqdci;               /* dwdot[OH]/d[H] */
    J[70] += dqdci;               /* dwdot[H]/d[H] */
    J[73] -= dqdci;               /* dwdot[CO]/d[H] */
    J[74] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3];
    J[107] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[109] += dqdci;              /* dwdot[H]/d[CO] */
    J[112] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[113] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[5];
    J[120] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[122] += dqdci;              /* dwdot[H]/d[CO2] */
    J[125] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[126] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[159] -= dqdT;               /* dwdot[OH]/dT */
    J[161] += dqdT;               /* dwdot[H]/dT */
    J[164] -= dqdT;               /* dwdot[CO]/dT */
    J[165] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 27: HCO + O2 <=> CO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[8];
    Kc = exp(g_RT[1] - g_RT[6] - g_RT[8] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[6] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[8] += q; /* CO */
    wdot[10] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[14] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[19] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[21] += dqdci;               /* dwdot[CO]/d[O2] */
    J[23] -= dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[8];
    J[79] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[84] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[86] += dqdci;               /* dwdot[CO]/d[HO2] */
    J[88] -= dqdci;               /* dwdot[HCO]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[105] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[110] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[112] += dqdci;              /* dwdot[CO]/d[CO] */
    J[114] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[131] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[136] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[138] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[140] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[157] -= dqdT;               /* dwdot[O2]/dT */
    J[162] += dqdT;               /* dwdot[HO2]/dT */
    J[164] += dqdT;               /* dwdot[CO]/dT */
    J[166] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 28: HCO + H <=> CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[10];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[8];
    Kc = exp(-g_RT[0] + g_RT[5] - g_RT[8] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[10]) + (h_RT[0] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[5] -= q; /* H */
    wdot[8] += q; /* CO */
    wdot[10] -= q; /* HCO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[8];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[5] -= dqdci;                /* dwdot[H]/d[H2] */
    J[8] += dqdci;                /* dwdot[CO]/d[H2] */
    J[10] -= dqdci;               /* dwdot[HCO]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[65] += dqdci;               /* dwdot[H2]/d[H] */
    J[70] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[CO]/d[H] */
    J[75] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[0];
    J[104] += dqdci;              /* dwdot[H2]/d[CO] */
    J[109] -= dqdci;              /* dwdot[H]/d[CO] */
    J[112] += dqdci;              /* dwdot[CO]/d[CO] */
    J[114] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[130] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[135] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[138] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[140] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[156] += dqdT;               /* dwdot[H2]/dT */
    J[161] -= dqdT;               /* dwdot[H]/dT */
    J[164] += dqdT;               /* dwdot[CO]/dT */
    J[166] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 29: HCO + O <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[10];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(g_RT[2] - g_RT[5] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[10]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[5] += q; /* H */
    wdot[9] += q; /* CO2 */
    wdot[10] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[28] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[H]/d[O] */
    J[35] += dqdci;               /* dwdot[CO2]/d[O] */
    J[36] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[67] -= dqdci;               /* dwdot[O]/d[H] */
    J[70] += dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[CO2]/d[H] */
    J[75] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[5];
    J[119] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[122] += dqdci;              /* dwdot[H]/d[CO2] */
    J[126] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[127] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[132] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[135] += dqdci;              /* dwdot[H]/d[HCO] */
    J[139] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[140] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[158] -= dqdT;               /* dwdot[O]/dT */
    J[161] += dqdT;               /* dwdot[H]/dT */
    J[165] += dqdT;               /* dwdot[CO2]/dT */
    J[166] -= dqdT;               /* dwdot[HCO]/dT */

    double c_R[12], dcRdT[12], e_RT[12];
    double * eh_RT;
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

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 12; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[156+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 12; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 12; ++m) {
            dehmixdc += eh_RT[m]*J[k*13+m];
        }
        J[k*13+12] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[168] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

/*compute an approx to the reaction Jacobian */
void aJacobian_precond(double *  J, double *  sc, double T, int HP)
{
    for (int i=0; i<169; i++) {
        J[i] = 0.0;
    }

    double wdot[12];
    for (int k=0; k<12; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 12; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[12];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[12];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[12];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[4] + (TB[0][2] - 1)*sc[1] + (TB[0][3] - 1)*sc[8] + (TB[0][4] - 1)*sc[9];
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = prefactor_units[0] * fwd_A[0]
                * exp(fwd_beta[0] * tc[0] - activation_units[0] * fwd_Ea[0] * invT);
    dlnkfdT = fwd_beta[0] * invT + activation_units[0] * fwd_Ea[0] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[0] * exp(low_beta[0] * tc[0] - activation_units[0] * low_Ea[0] * invT);
    Pr = phase_units[0] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[0] * invT + activation_units[0] * low_Ea[0] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[0]) > 1.e-100 ? (1.-troe_a[0])*exp(-T/troe_Tsss[0]) : 0.);
    Fcent2 = (fabs(troe_Ts[0]) > 1.e-100 ? troe_a[0] * exp(-T/troe_Ts[0]) : 0.);
    Fcent3 = (troe_len[0] == 4 ? exp(-troe_Tss[0] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[0]) > 1.e-100 ? -Fcent1/troe_Tsss[0] : 0.)
      + (fabs(troe_Ts[0]) > 1.e-100 ? -Fcent2/troe_Ts[0] : 0.)
      + (troe_len[0] == 4 ? Fcent3*troe_Tss[0]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[5] -= q; /* H */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = TB[0][0]*dcdc_fac;
    dqdc[1] = TB[0][2]*dcdc_fac + k_f*sc[5];
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = TB[0][1]*dcdc_fac;
    dqdc[5] = dcdc_fac + k_f*sc[1];
    dqdc[6] = dcdc_fac - k_r;
    dqdc[7] = dcdc_fac;
    dqdc[8] = TB[0][3]*dcdc_fac;
    dqdc[9] = TB[0][4]*dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    for (int k=0; k<12; k++) {
        J[13*k+1] -= dqdc[k];
        J[13*k+5] -= dqdc[k];
        J[13*k+6] += dqdc[k];
    }
    J[157] -= dqdT; /* dwdot[O2]/dT */
    J[161] -= dqdT; /* dwdot[H]/dT */
    J[162] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[4] + (TB[1][2] - 1)*sc[8] + (TB[1][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[7];
    k_f = prefactor_units[1] * fwd_A[1]
                * exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
    dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[1] * exp(low_beta[1] * tc[0] - activation_units[1] * low_Ea[1] * invT);
    Pr = phase_units[1] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[1] * invT + activation_units[1] * low_Ea[1] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[1]) > 1.e-100 ? (1.-troe_a[1])*exp(-T/troe_Tsss[1]) : 0.);
    Fcent2 = (fabs(troe_Ts[1]) > 1.e-100 ? troe_a[1] * exp(-T/troe_Ts[1]) : 0.);
    Fcent3 = (troe_len[1] == 4 ? exp(-troe_Tss[1] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[1]) > 1.e-100 ? -Fcent1/troe_Tsss[1] : 0.)
      + (fabs(troe_Ts[1]) > 1.e-100 ? -Fcent2/troe_Ts[1] : 0.)
      + (troe_len[1] == 4 ? Fcent3*troe_Tss[1]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[3]*sc[3];
    Kc = refC * exp(-g_RT[3] - g_RT[3] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2*h_RT[3]) - 1);
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
    dqdc[0] = TB[1][0]*dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac - k_r*2*sc[3];
    dqdc[4] = TB[1][1]*dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac + k_f;
    dqdc[8] = TB[1][2]*dcdc_fac;
    dqdc[9] = TB[1][3]*dcdc_fac;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    for (int k=0; k<12; k++) {
        J[13*k+3] += 2 * dqdc[k];
        J[13*k+7] -= dqdc[k];
    }
    J[159] += 2 * dqdT; /* dwdot[OH]/dT */
    J[163] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[4] + (TB[2][2] - 1)*sc[8] + (TB[2][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[2] * exp(low_beta[2] * tc[0] - activation_units[2] * low_Ea[2] * invT);
    Pr = phase_units[2] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[2] * invT + activation_units[2] * low_Ea[2] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Lindemann form */
    F = 1.0;
    dlogFdlogPr = 0.0;
    dlogFdT = 0.0;
    /* reverse */
    phi_r = sc[9];
    Kc = refCinv * exp(g_RT[2] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[9]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = TB[2][0]*dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac + k_f*sc[8];
    dqdc[3] = dcdc_fac;
    dqdc[4] = TB[2][1]*dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = TB[2][2]*dcdc_fac + k_f*sc[2];
    dqdc[9] = TB[2][3]*dcdc_fac - k_r;
    dqdc[10] = dcdc_fac;
    dqdc[11] = dcdc_fac;
    for (int k=0; k<12; k++) {
        J[13*k+2] -= dqdc[k];
        J[13*k+8] -= dqdc[k];
        J[13*k+9] += dqdc[k];
    }
    J[158] -= dqdT; /* dwdot[O]/dT */
    J[164] -= dqdT; /* dwdot[CO]/dT */
    J[165] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 4: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[4] + (TB[3][2] - 1)*sc[8] + (TB[3][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[0];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[5];
    Kc = refC * exp(g_RT[0] - g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0]) + (2*h_RT[5]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[5] += 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = TB[3][0]*q_nocor + k_f;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = TB[3][1]*q_nocor;
    dqdc[5] = q_nocor - k_r*2*sc[5];
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = TB[3][2]*q_nocor;
    dqdc[9] = TB[3][3]*q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    for (int k=0; k<12; k++) {
        J[13*k+0] -= dqdc[k];
        J[13*k+5] += 2 * dqdc[k];
    }
    J[156] -= dqdT; /* dwdot[H2]/dT */
    J[161] += 2 * dqdT; /* dwdot[H]/dT */

    /*reaction 5: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[4] + (TB[4][2] - 1)*sc[8] + (TB[4][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[2]*sc[2];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[1];
    Kc = refCinv * exp(-g_RT[1] + g_RT[2] + g_RT[2]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[2]) + (h_RT[1]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[2] -= 2 * q; /* O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = TB[4][0]*q_nocor;
    dqdc[1] = q_nocor - k_r;
    dqdc[2] = q_nocor + k_f*2*sc[2];
    dqdc[3] = q_nocor;
    dqdc[4] = TB[4][1]*q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = TB[4][2]*q_nocor;
    dqdc[9] = TB[4][3]*q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    for (int k=0; k<12; k++) {
        J[13*k+1] += dqdc[k];
        J[13*k+2] += -2 * dqdc[k];
    }
    J[157] += dqdT; /* dwdot[O2]/dT */
    J[158] += -2 * dqdT; /* dwdot[O]/dT */

    /*reaction 6: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[4] + (TB[5][2] - 1)*sc[8] + (TB[5][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[2]*sc[5];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[3]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = TB[5][0]*q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor + k_f*sc[5];
    dqdc[3] = q_nocor - k_r;
    dqdc[4] = TB[5][1]*q_nocor;
    dqdc[5] = q_nocor + k_f*sc[2];
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = TB[5][2]*q_nocor;
    dqdc[9] = TB[5][3]*q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    for (int k=0; k<12; k++) {
        J[13*k+2] -= dqdc[k];
        J[13*k+3] += dqdc[k];
        J[13*k+5] -= dqdc[k];
    }
    J[158] -= dqdT; /* dwdot[O]/dT */
    J[159] += dqdT; /* dwdot[OH]/dT */
    J[161] -= dqdT; /* dwdot[H]/dT */

    /*reaction 7: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[4] + (TB[6][2] - 1)*sc[8] + (TB[6][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[3]*sc[5];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[3] - g_RT[4] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[4]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] -= q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = TB[6][0]*q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor + k_f*sc[5];
    dqdc[4] = TB[6][1]*q_nocor - k_r;
    dqdc[5] = q_nocor + k_f*sc[3];
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = TB[6][2]*q_nocor;
    dqdc[9] = TB[6][3]*q_nocor;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    for (int k=0; k<12; k++) {
        J[13*k+3] -= dqdc[k];
        J[13*k+4] += dqdc[k];
        J[13*k+5] -= dqdc[k];
    }
    J[159] -= dqdT; /* dwdot[OH]/dT */
    J[160] += dqdT; /* dwdot[H2O]/dT */
    J[161] -= dqdT; /* dwdot[H]/dT */

    /*reaction 8: HCO + M <=> H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[7][0] - 1)*sc[0] + (TB[7][1] - 1)*sc[4] + (TB[7][2] - 1)*sc[8] + (TB[7][3] - 1)*sc[9];
    /* forward */
    phi_f = sc[10];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[8];
    Kc = refC * exp(-g_RT[5] - g_RT[8] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10]) + (h_RT[5] + h_RT[8]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* H */
    wdot[8] += q; /* CO */
    wdot[10] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = TB[7][0]*q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = TB[7][1]*q_nocor;
    dqdc[5] = q_nocor - k_r*sc[8];
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = TB[7][2]*q_nocor - k_r*sc[5];
    dqdc[9] = TB[7][3]*q_nocor;
    dqdc[10] = q_nocor + k_f;
    dqdc[11] = q_nocor;
    for (int k=0; k<12; k++) {
        J[13*k+5] += dqdc[k];
        J[13*k+8] += dqdc[k];
        J[13*k+10] -= dqdc[k];
    }
    J[161] += dqdT; /* dwdot[H]/dT */
    J[164] += dqdT; /* dwdot[CO]/dT */
    J[166] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 9: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[3];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[2] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[2] += q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* H */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[5];
    J[14] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[15] += dqdci;               /* dwdot[O]/d[O2] */
    J[16] += dqdci;               /* dwdot[OH]/d[O2] */
    J[18] -= dqdci;               /* dwdot[H]/d[O2] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[3];
    J[27] -= dqdci;               /* dwdot[O2]/d[O] */
    J[28] += dqdci;               /* dwdot[O]/d[O] */
    J[29] += dqdci;               /* dwdot[OH]/d[O] */
    J[31] -= dqdci;               /* dwdot[H]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[40] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[41] += dqdci;               /* dwdot[O]/d[OH] */
    J[42] += dqdci;               /* dwdot[OH]/d[OH] */
    J[44] -= dqdci;               /* dwdot[H]/d[OH] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[1];
    J[66] -= dqdci;               /* dwdot[O2]/d[H] */
    J[67] += dqdci;               /* dwdot[O]/d[H] */
    J[68] += dqdci;               /* dwdot[OH]/d[H] */
    J[70] -= dqdci;               /* dwdot[H]/d[H] */
    /* d()/dT */
    J[157] -= dqdT;               /* dwdot[O2]/dT */
    J[158] += dqdT;               /* dwdot[O]/dT */
    J[159] += dqdT;               /* dwdot[OH]/dT */
    J[161] -= dqdT;               /* dwdot[H]/dT */

    /*reaction 10: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[0] + g_RT[2] - g_RT[3] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] += q; /* H */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[2] -= dqdci;                /* dwdot[O]/d[H2] */
    J[3] += dqdci;                /* dwdot[OH]/d[H2] */
    J[5] += dqdci;                /* dwdot[H]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[26] -= dqdci;               /* dwdot[H2]/d[O] */
    J[28] -= dqdci;               /* dwdot[O]/d[O] */
    J[29] += dqdci;               /* dwdot[OH]/d[O] */
    J[31] += dqdci;               /* dwdot[H]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[39] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[41] -= dqdci;               /* dwdot[O]/d[OH] */
    J[42] += dqdci;               /* dwdot[OH]/d[OH] */
    J[44] += dqdci;               /* dwdot[H]/d[OH] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[65] -= dqdci;               /* dwdot[H2]/d[H] */
    J[67] -= dqdci;               /* dwdot[O]/d[H] */
    J[68] += dqdci;               /* dwdot[OH]/d[H] */
    J[70] += dqdci;               /* dwdot[H]/d[H] */
    /* d()/dT */
    J[156] -= dqdT;               /* dwdot[H2]/dT */
    J[158] -= dqdT;               /* dwdot[O]/dT */
    J[159] += dqdT;               /* dwdot[OH]/dT */
    J[161] += dqdT;               /* dwdot[H]/dT */

    /*reaction 11: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[3];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[0] + g_RT[3] - g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[3]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] += q; /* H */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[3] -= dqdci;                /* dwdot[OH]/d[H2] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H2] */
    J[5] += dqdci;                /* dwdot[H]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[39] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[42] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[43] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[44] += dqdci;               /* dwdot[H]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[52] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[55] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[56] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[57] += dqdci;               /* dwdot[H]/d[H2O] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[65] -= dqdci;               /* dwdot[H2]/d[H] */
    J[68] -= dqdci;               /* dwdot[OH]/d[H] */
    J[69] += dqdci;               /* dwdot[H2O]/d[H] */
    J[70] += dqdci;               /* dwdot[H]/d[H] */
    /* d()/dT */
    J[156] -= dqdT;               /* dwdot[H2]/dT */
    J[159] -= dqdT;               /* dwdot[OH]/dT */
    J[160] += dqdT;               /* dwdot[H2O]/dT */
    J[161] += dqdT;               /* dwdot[H]/dT */

    /*reaction 12: O + H2O <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[3];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[3] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2*h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += 2 * q; /* OH */
    wdot[4] -= q; /* H2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[4];
    J[28] -= dqdci;               /* dwdot[O]/d[O] */
    J[29] += 2 * dqdci;           /* dwdot[OH]/d[O] */
    J[30] -= dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[3];
    J[41] -= dqdci;               /* dwdot[O]/d[OH] */
    J[42] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[2];
    J[54] -= dqdci;               /* dwdot[O]/d[H2O] */
    J[55] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
    J[56] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[158] -= dqdT;               /* dwdot[O]/dT */
    J[159] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[160] -= dqdT;               /* dwdot[H2O]/dT */

    /*reaction 13: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[1];
    Kc = exp(-g_RT[0] - g_RT[1] + g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (h_RT[0] + h_RT[1]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] += q; /* O2 */
    wdot[5] -= q; /* H */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[1];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[O2]/d[H2] */
    J[5] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[13] += dqdci;               /* dwdot[H2]/d[O2] */
    J[14] += dqdci;               /* dwdot[O2]/d[O2] */
    J[18] -= dqdci;               /* dwdot[H]/d[O2] */
    J[19] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[65] += dqdci;               /* dwdot[H2]/d[H] */
    J[66] += dqdci;               /* dwdot[O2]/d[H] */
    J[70] -= dqdci;               /* dwdot[H]/d[H] */
    J[71] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[78] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[79] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[83] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[84] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[156] += dqdT;               /* dwdot[H2]/dT */
    J[157] += dqdT;               /* dwdot[O2]/dT */
    J[161] -= dqdT;               /* dwdot[H]/dT */
    J[162] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 14: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[3];
    Kc = exp(-g_RT[3] - g_RT[3] + g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (2*h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += 2 * q; /* OH */
    wdot[5] -= q; /* H */
    wdot[6] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[3];
    J[42] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[44] -= dqdci;               /* dwdot[H]/d[OH] */
    J[45] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[68] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[70] -= dqdci;               /* dwdot[H]/d[H] */
    J[71] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[81] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[83] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[84] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[159] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[161] -= dqdT;               /* dwdot[H]/dT */
    J[162] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 15: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[3];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[1] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[3];
    J[14] += dqdci;               /* dwdot[O2]/d[O2] */
    J[15] -= dqdci;               /* dwdot[O]/d[O2] */
    J[16] += dqdci;               /* dwdot[OH]/d[O2] */
    J[19] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[27] += dqdci;               /* dwdot[O2]/d[O] */
    J[28] -= dqdci;               /* dwdot[O]/d[O] */
    J[29] += dqdci;               /* dwdot[OH]/d[O] */
    J[32] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[40] += dqdci;               /* dwdot[O2]/d[OH] */
    J[41] -= dqdci;               /* dwdot[O]/d[OH] */
    J[42] += dqdci;               /* dwdot[OH]/d[OH] */
    J[45] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[79] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[80] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[81] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[84] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[157] += dqdT;               /* dwdot[O2]/dT */
    J[158] -= dqdT;               /* dwdot[O]/dT */
    J[159] += dqdT;               /* dwdot[OH]/dT */
    J[162] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4];
    Kc = exp(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[1] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[14] += dqdci;               /* dwdot[O2]/d[O2] */
    J[16] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[17] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[19] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[40] += dqdci;               /* dwdot[O2]/d[OH] */
    J[42] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[43] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[45] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[55] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[56] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[58] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[79] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[81] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[82] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[84] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[157] += dqdT;               /* dwdot[O2]/dT */
    J[159] -= dqdT;               /* dwdot[OH]/dT */
    J[160] += dqdT;               /* dwdot[H2O]/dT */
    J[162] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 17: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[1] + h_RT[7]));
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
    J[14] += dqdci;               /* dwdot[O2]/d[O2] */
    J[19] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[20] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[6];
    J[79] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[84] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[92] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[97] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[98] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[157] += dqdT;               /* dwdot[O2]/dT */
    J[162] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[163] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 18: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[1] + h_RT[7]));
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
    J[14] += dqdci;               /* dwdot[O2]/d[O2] */
    J[19] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[20] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[6];
    J[79] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[84] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[92] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[97] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[98] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[157] += dqdT;               /* dwdot[O2]/dT */
    J[162] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[163] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 19: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(-g_RT[3] - g_RT[4] + g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] -= q; /* H */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[42] += dqdci;               /* dwdot[OH]/d[OH] */
    J[43] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[44] -= dqdci;               /* dwdot[H]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[55] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[56] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[59] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[68] += dqdci;               /* dwdot[OH]/d[H] */
    J[69] += dqdci;               /* dwdot[H2O]/d[H] */
    J[70] -= dqdci;               /* dwdot[H]/d[H] */
    J[72] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[94] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[95] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[96] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[98] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[159] += dqdT;               /* dwdot[OH]/dT */
    J[160] += dqdT;               /* dwdot[H2O]/dT */
    J[161] -= dqdT;               /* dwdot[H]/dT */
    J[163] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 20: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[5] -= q; /* H */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[5] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[65] += dqdci;               /* dwdot[H2]/d[H] */
    J[70] -= dqdci;               /* dwdot[H]/d[H] */
    J[71] += dqdci;               /* dwdot[HO2]/d[H] */
    J[72] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[78] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[83] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[84] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[91] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[96] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[97] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[98] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[156] += dqdT;               /* dwdot[H2]/dT */
    J[161] -= dqdT;               /* dwdot[H]/dT */
    J[162] += dqdT;               /* dwdot[HO2]/dT */
    J[163] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 21: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
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
    J[28] -= dqdci;               /* dwdot[O]/d[O] */
    J[29] += dqdci;               /* dwdot[OH]/d[O] */
    J[32] += dqdci;               /* dwdot[HO2]/d[O] */
    J[33] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[41] -= dqdci;               /* dwdot[O]/d[OH] */
    J[42] += dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[80] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[81] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[84] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[93] -= dqdci;               /* dwdot[O]/d[H2O2] */
    J[94] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[97] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[98] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[158] -= dqdT;               /* dwdot[O]/dT */
    J[159] += dqdT;               /* dwdot[OH]/dT */
    J[162] += dqdT;               /* dwdot[HO2]/dT */
    J[163] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 22: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
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
    J[42] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[43] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[45] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[55] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[56] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[58] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[59] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[81] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[82] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[84] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[94] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[95] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[97] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[98] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[159] -= dqdT;               /* dwdot[OH]/dT */
    J[160] += dqdT;               /* dwdot[H2O]/dT */
    J[162] += dqdT;               /* dwdot[HO2]/dT */
    J[163] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 23: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
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
    J[42] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[43] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[45] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[55] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[56] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[58] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[59] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[81] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[82] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[84] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[85] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[94] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[95] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[97] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[98] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[159] -= dqdT;               /* dwdot[OH]/dT */
    J[160] += dqdT;               /* dwdot[H2O]/dT */
    J[162] += dqdT;               /* dwdot[HO2]/dT */
    J[163] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 24: CO + O2 <=> CO2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[8];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[9];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[8]) + (h_RT[2] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[2] += q; /* O */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[14] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[15] += dqdci;               /* dwdot[O]/d[O2] */
    J[21] -= dqdci;               /* dwdot[CO]/d[O2] */
    J[22] += dqdci;               /* dwdot[CO2]/d[O2] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[9];
    J[27] -= dqdci;               /* dwdot[O2]/d[O] */
    J[28] += dqdci;               /* dwdot[O]/d[O] */
    J[34] -= dqdci;               /* dwdot[CO]/d[O] */
    J[35] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[1];
    J[105] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[106] += dqdci;              /* dwdot[O]/d[CO] */
    J[112] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[113] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[2];
    J[118] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[119] += dqdci;              /* dwdot[O]/d[CO2] */
    J[125] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[126] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[157] -= dqdT;               /* dwdot[O2]/dT */
    J[158] += dqdT;               /* dwdot[O]/dT */
    J[164] -= dqdT;               /* dwdot[CO]/dT */
    J[165] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 25: CO + HO2 <=> CO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[8];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[9];
    Kc = exp(-g_RT[3] + g_RT[6] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[3] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[42] += dqdci;               /* dwdot[OH]/d[OH] */
    J[45] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[48] += dqdci;               /* dwdot[CO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[8];
    J[81] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[84] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    J[86] -= dqdci;               /* dwdot[CO]/d[HO2] */
    J[87] += dqdci;               /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[107] += dqdci;              /* dwdot[OH]/d[CO] */
    J[110] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[112] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[113] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[3];
    J[120] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[123] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[125] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[126] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[159] += dqdT;               /* dwdot[OH]/dT */
    J[162] -= dqdT;               /* dwdot[HO2]/dT */
    J[164] -= dqdT;               /* dwdot[CO]/dT */
    J[165] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 26: CO + OH <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(g_RT[3] - g_RT[5] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[5] += q; /* H */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[42] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[44] += dqdci;               /* dwdot[H]/d[OH] */
    J[47] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[48] += dqdci;               /* dwdot[CO2]/d[OH] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[68] -= dqdci;               /* dwdot[OH]/d[H] */
    J[70] += dqdci;               /* dwdot[H]/d[H] */
    J[73] -= dqdci;               /* dwdot[CO]/d[H] */
    J[74] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3];
    J[107] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[109] += dqdci;              /* dwdot[H]/d[CO] */
    J[112] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[113] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[5];
    J[120] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[122] += dqdci;              /* dwdot[H]/d[CO2] */
    J[125] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[126] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[159] -= dqdT;               /* dwdot[OH]/dT */
    J[161] += dqdT;               /* dwdot[H]/dT */
    J[164] -= dqdT;               /* dwdot[CO]/dT */
    J[165] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 27: HCO + O2 <=> CO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[8];
    Kc = exp(g_RT[1] - g_RT[6] - g_RT[8] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[6] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[8] += q; /* CO */
    wdot[10] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[14] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[19] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[21] += dqdci;               /* dwdot[CO]/d[O2] */
    J[23] -= dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[8];
    J[79] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[84] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[86] += dqdci;               /* dwdot[CO]/d[HO2] */
    J[88] -= dqdci;               /* dwdot[HCO]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[105] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[110] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[112] += dqdci;              /* dwdot[CO]/d[CO] */
    J[114] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[131] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[136] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[138] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[140] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[157] -= dqdT;               /* dwdot[O2]/dT */
    J[162] += dqdT;               /* dwdot[HO2]/dT */
    J[164] += dqdT;               /* dwdot[CO]/dT */
    J[166] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 28: HCO + H <=> CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[10];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[8];
    Kc = exp(-g_RT[0] + g_RT[5] - g_RT[8] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[10]) + (h_RT[0] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[5] -= q; /* H */
    wdot[8] += q; /* CO */
    wdot[10] -= q; /* HCO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[8];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[5] -= dqdci;                /* dwdot[H]/d[H2] */
    J[8] += dqdci;                /* dwdot[CO]/d[H2] */
    J[10] -= dqdci;               /* dwdot[HCO]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[65] += dqdci;               /* dwdot[H2]/d[H] */
    J[70] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[CO]/d[H] */
    J[75] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[0];
    J[104] += dqdci;              /* dwdot[H2]/d[CO] */
    J[109] -= dqdci;              /* dwdot[H]/d[CO] */
    J[112] += dqdci;              /* dwdot[CO]/d[CO] */
    J[114] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[130] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[135] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[138] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[140] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[156] += dqdT;               /* dwdot[H2]/dT */
    J[161] -= dqdT;               /* dwdot[H]/dT */
    J[164] += dqdT;               /* dwdot[CO]/dT */
    J[166] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 29: HCO + O <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[10];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(g_RT[2] - g_RT[5] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[10]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[5] += q; /* H */
    wdot[9] += q; /* CO2 */
    wdot[10] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[28] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[H]/d[O] */
    J[35] += dqdci;               /* dwdot[CO2]/d[O] */
    J[36] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[67] -= dqdci;               /* dwdot[O]/d[H] */
    J[70] += dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[CO2]/d[H] */
    J[75] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[5];
    J[119] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[122] += dqdci;              /* dwdot[H]/d[CO2] */
    J[126] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[127] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[132] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[135] += dqdci;              /* dwdot[H]/d[HCO] */
    J[139] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[140] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[158] -= dqdT;               /* dwdot[O]/dT */
    J[161] += dqdT;               /* dwdot[H]/dT */
    J[165] += dqdT;               /* dwdot[CO2]/dT */
    J[166] -= dqdT;               /* dwdot[HCO]/dT */

    double c_R[12], dcRdT[12], e_RT[12];
    double * eh_RT;
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

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 12; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[156+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 12; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 12; ++m) {
            dehmixdc += eh_RT[m]*J[k*13+m];
        }
        J[k*13+12] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[168] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void dcvpRdT(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +8.24944200e-04
            -1.62860300e-06 * tc[1]
            -2.84263020e-10 * tc[2]
            +1.65394880e-12 * tc[3];
        /*species 1: O2 */
        species[1] =
            +1.12748600e-03
            -1.15123000e-06 * tc[1]
            +3.94163100e-09 * tc[2]
            -3.50742160e-12 * tc[3];
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
        /*species 5: H */
        species[5] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
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
        /*species 8: CO */
        species[8] =
            +1.51194100e-03
            -7.76351000e-06 * tc[1]
            +1.67458320e-08 * tc[2]
            -9.89980400e-12 * tc[3];
        /*species 9: CO2 */
        species[9] =
            +9.92207200e-03
            -2.08182200e-05 * tc[1]
            +2.06000610e-08 * tc[2]
            -8.46912000e-12 * tc[3];
        /*species 10: HCO */
        species[10] =
            +6.19914700e-03
            -1.92461680e-05 * tc[1]
            +3.26947500e-08 * tc[2]
            -1.82995400e-11 * tc[3];
        /*species 11: N2 */
        species[11] =
            +1.40824000e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77942000e-12 * tc[3];
    } else {
        /*species 0: H2 */
        species[0] =
            +7.00064400e-04
            -1.12676580e-07 * tc[1]
            -2.76947340e-11 * tc[2]
            +6.33100800e-15 * tc[3];
        /*species 1: O2 */
        species[1] =
            +6.13519700e-04
            -2.51768400e-07 * tc[1]
            +5.32584300e-11 * tc[2]
            -4.54574000e-15 * tc[3];
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
        /*species 5: H */
        species[5] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
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
        /*species 8: CO */
        species[8] =
            +1.44268900e-03
            -1.12616560e-06 * tc[1]
            +3.05574300e-10 * tc[2]
            -2.76438080e-14 * tc[3];
        /*species 9: CO2 */
        species[9] =
            +3.14016900e-03
            -2.55682200e-06 * tc[1]
            +7.18199100e-10 * tc[2]
            -6.67613200e-14 * tc[3];
        /*species 10: HCO */
        species[10] =
            +3.34557300e-03
            -2.67001200e-06 * tc[1]
            +7.41171900e-10 * tc[2]
            -6.85540400e-14 * tc[3];
        /*species 11: N2 */
        species[11] =
            +1.48797700e-03
            -1.13695220e-06 * tc[1]
            +3.02911200e-10 * tc[2]
            -2.70134040e-14 * tc[3];
    }
    return;
}


/*compute the progress rate for each reaction */
void progressRate(double *  qdot, double *  sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double q_f[29], q_r[29];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 29; ++i) {
        qdot[i] = q_f[i] - q_r[i];
    }

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double *  q_f, double *  q_r, double *  sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    comp_qfqr(q_f, q_r, sc, tc, invT);

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *  kc, double *  g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    kc[0] = 1.0 / (refC) * exp((g_RT[5] + g_RT[1]) - (g_RT[6]));

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    kc[1] = refC * exp((g_RT[7]) - (g_RT[3] + g_RT[3]));

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    kc[2] = 1.0 / (refC) * exp((g_RT[8] + g_RT[2]) - (g_RT[9]));

    /*reaction 4: H2 + M <=> H + H + M */
    kc[3] = refC * exp((g_RT[0]) - (g_RT[5] + g_RT[5]));

    /*reaction 5: O + O + M <=> O2 + M */
    kc[4] = 1.0 / (refC) * exp((g_RT[2] + g_RT[2]) - (g_RT[1]));

    /*reaction 6: O + H + M <=> OH + M */
    kc[5] = 1.0 / (refC) * exp((g_RT[2] + g_RT[5]) - (g_RT[3]));

    /*reaction 7: H + OH + M <=> H2O + M */
    kc[6] = 1.0 / (refC) * exp((g_RT[5] + g_RT[3]) - (g_RT[4]));

    /*reaction 8: HCO + M <=> H + CO + M */
    kc[7] = refC * exp((g_RT[10]) - (g_RT[5] + g_RT[8]));

    /*reaction 9: H + O2 <=> O + OH */
    kc[8] = exp((g_RT[5] + g_RT[1]) - (g_RT[2] + g_RT[3]));

    /*reaction 10: O + H2 <=> H + OH */
    kc[9] = exp((g_RT[2] + g_RT[0]) - (g_RT[5] + g_RT[3]));

    /*reaction 11: H2 + OH <=> H2O + H */
    kc[10] = exp((g_RT[0] + g_RT[3]) - (g_RT[4] + g_RT[5]));

    /*reaction 12: O + H2O <=> OH + OH */
    kc[11] = exp((g_RT[2] + g_RT[4]) - (g_RT[3] + g_RT[3]));

    /*reaction 13: HO2 + H <=> H2 + O2 */
    kc[12] = exp((g_RT[6] + g_RT[5]) - (g_RT[0] + g_RT[1]));

    /*reaction 14: HO2 + H <=> OH + OH */
    kc[13] = exp((g_RT[6] + g_RT[5]) - (g_RT[3] + g_RT[3]));

    /*reaction 15: HO2 + O <=> O2 + OH */
    kc[14] = exp((g_RT[6] + g_RT[2]) - (g_RT[1] + g_RT[3]));

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    kc[15] = exp((g_RT[6] + g_RT[3]) - (g_RT[4] + g_RT[1]));

    /*reaction 17: HO2 + HO2 <=> H2O2 + O2 */
    kc[16] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));

    /*reaction 18: HO2 + HO2 <=> H2O2 + O2 */
    kc[17] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));

    /*reaction 19: H2O2 + H <=> H2O + OH */
    kc[18] = exp((g_RT[7] + g_RT[5]) - (g_RT[4] + g_RT[3]));

    /*reaction 20: H2O2 + H <=> HO2 + H2 */
    kc[19] = exp((g_RT[7] + g_RT[5]) - (g_RT[6] + g_RT[0]));

    /*reaction 21: H2O2 + O <=> OH + HO2 */
    kc[20] = exp((g_RT[7] + g_RT[2]) - (g_RT[3] + g_RT[6]));

    /*reaction 22: H2O2 + OH <=> HO2 + H2O */
    kc[21] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[4]));

    /*reaction 23: H2O2 + OH <=> HO2 + H2O */
    kc[22] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[4]));

    /*reaction 24: CO + O2 <=> CO2 + O */
    kc[23] = exp((g_RT[8] + g_RT[1]) - (g_RT[9] + g_RT[2]));

    /*reaction 25: CO + HO2 <=> CO2 + OH */
    kc[24] = exp((g_RT[8] + g_RT[6]) - (g_RT[9] + g_RT[3]));

    /*reaction 26: CO + OH <=> CO2 + H */
    kc[25] = exp((g_RT[8] + g_RT[3]) - (g_RT[9] + g_RT[5]));

    /*reaction 27: HCO + O2 <=> CO + HO2 */
    kc[26] = exp((g_RT[10] + g_RT[1]) - (g_RT[8] + g_RT[6]));

    /*reaction 28: HCO + H <=> CO + H2 */
    kc[27] = exp((g_RT[10] + g_RT[5]) - (g_RT[8] + g_RT[0]));

    /*reaction 29: HCO + O <=> CO2 + H */
    kc[28] = exp((g_RT[10] + g_RT[2]) - (g_RT[9] + g_RT[5]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -1.012521000000000e+03 * invT
            +6.592218000000000e+00
            -3.298124000000000e+00 * tc[0]
            -4.124721000000000e-04 * tc[1]
            +1.357169166666667e-07 * tc[2]
            +7.896194999999999e-12 * tc[3]
            -2.067436000000000e-14 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.005249000000000e+03 * invT
            -2.821802000000000e+00
            -3.212936000000000e+00 * tc[0]
            -5.637430000000000e-04 * tc[1]
            +9.593583333333333e-08 * tc[2]
            -1.094897500000000e-10 * tc[3]
            +4.384277000000000e-14 * tc[4];
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
        /*species 5: H */
        species[5] =
            +2.547163000000000e+04 * invT
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
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
        /*species 8: CO */
        species[8] =
            -1.431054000000000e+04 * invT
            -1.586445000000000e+00
            -3.262452000000000e+00 * tc[0]
            -7.559705000000000e-04 * tc[1]
            +6.469591666666667e-07 * tc[2]
            -4.651620000000000e-10 * tc[3]
            +1.237475500000000e-13 * tc[4];
        /*species 9: CO2 */
        species[9] =
            -4.837314000000000e+04 * invT
            -7.912765000000000e+00
            -2.275725000000000e+00 * tc[0]
            -4.961036000000000e-03 * tc[1]
            +1.734851666666667e-06 * tc[2]
            -5.722239166666667e-10 * tc[3]
            +1.058640000000000e-13 * tc[4];
        /*species 10: HCO */
        species[10] =
            +4.159922000000000e+03 * invT
            -6.085284000000000e+00
            -2.898330000000000e+00 * tc[0]
            -3.099573500000000e-03 * tc[1]
            +1.603847333333333e-06 * tc[2]
            -9.081875000000000e-10 * tc[3]
            +2.287442500000000e-13 * tc[4];
        /*species 11: N2 */
        species[11] =
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
            -8.350340000000000e+02 * invT
            +4.346533000000000e+00
            -2.991423000000000e+00 * tc[0]
            -3.500322000000000e-04 * tc[1]
            +9.389715000000000e-09 * tc[2]
            +7.692981666666667e-13 * tc[3]
            -7.913760000000000e-17 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.233930000000000e+03 * invT
            +5.084119999999999e-01
            -3.697578000000000e+00 * tc[0]
            -3.067598500000000e-04 * tc[1]
            +2.098070000000000e-08 * tc[2]
            -1.479400833333333e-12 * tc[3]
            +5.682175000000001e-17 * tc[4];
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
        /*species 5: H */
        species[5] =
            +2.547163000000000e+04 * invT
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
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
        /*species 8: CO */
        species[8] =
            -1.426835000000000e+04 * invT
            -3.083140000000000e+00
            -3.025078000000000e+00 * tc[0]
            -7.213445000000000e-04 * tc[1]
            +9.384713333333334e-08 * tc[2]
            -8.488174999999999e-12 * tc[3]
            +3.455476000000000e-16 * tc[4];
        /*species 9: CO2 */
        species[9] =
            -4.896696000000000e+04 * invT
            +5.409018900000000e+00
            -4.453623000000000e+00 * tc[0]
            -1.570084500000000e-03 * tc[1]
            +2.130685000000000e-07 * tc[2]
            -1.994997500000000e-11 * tc[3]
            +8.345165000000000e-16 * tc[4];
        /*species 10: HCO */
        species[10] =
            +3.916324000000000e+03 * invT
            -1.995028000000000e+00
            -3.557271000000000e+00 * tc[0]
            -1.672786500000000e-03 * tc[1]
            +2.225010000000000e-07 * tc[2]
            -2.058810833333333e-11 * tc[3]
            +8.569255000000000e-16 * tc[4];
        /*species 11: N2 */
        species[11] =
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
void helmholtz(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -1.01252100e+03 * invT
            +5.59221800e+00
            -3.29812400e+00 * tc[0]
            -4.12472100e-04 * tc[1]
            +1.35716917e-07 * tc[2]
            +7.89619500e-12 * tc[3]
            -2.06743600e-14 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.00524900e+03 * invT
            -3.82180200e+00
            -3.21293600e+00 * tc[0]
            -5.63743000e-04 * tc[1]
            +9.59358333e-08 * tc[2]
            -1.09489750e-10 * tc[3]
            +4.38427700e-14 * tc[4];
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
        /*species 5: H */
        species[5] =
            +2.54716300e+04 * invT
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
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
        /*species 8: CO */
        species[8] =
            -1.43105400e+04 * invT
            -2.58644500e+00
            -3.26245200e+00 * tc[0]
            -7.55970500e-04 * tc[1]
            +6.46959167e-07 * tc[2]
            -4.65162000e-10 * tc[3]
            +1.23747550e-13 * tc[4];
        /*species 9: CO2 */
        species[9] =
            -4.83731400e+04 * invT
            -8.91276500e+00
            -2.27572500e+00 * tc[0]
            -4.96103600e-03 * tc[1]
            +1.73485167e-06 * tc[2]
            -5.72223917e-10 * tc[3]
            +1.05864000e-13 * tc[4];
        /*species 10: HCO */
        species[10] =
            +4.15992200e+03 * invT
            -7.08528400e+00
            -2.89833000e+00 * tc[0]
            -3.09957350e-03 * tc[1]
            +1.60384733e-06 * tc[2]
            -9.08187500e-10 * tc[3]
            +2.28744250e-13 * tc[4];
        /*species 11: N2 */
        species[11] =
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
            -8.35034000e+02 * invT
            +3.34653300e+00
            -2.99142300e+00 * tc[0]
            -3.50032200e-04 * tc[1]
            +9.38971500e-09 * tc[2]
            +7.69298167e-13 * tc[3]
            -7.91376000e-17 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.23393000e+03 * invT
            -4.91588000e-01
            -3.69757800e+00 * tc[0]
            -3.06759850e-04 * tc[1]
            +2.09807000e-08 * tc[2]
            -1.47940083e-12 * tc[3]
            +5.68217500e-17 * tc[4];
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
        /*species 5: H */
        species[5] =
            +2.54716300e+04 * invT
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
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
        /*species 8: CO */
        species[8] =
            -1.42683500e+04 * invT
            -4.08314000e+00
            -3.02507800e+00 * tc[0]
            -7.21344500e-04 * tc[1]
            +9.38471333e-08 * tc[2]
            -8.48817500e-12 * tc[3]
            +3.45547600e-16 * tc[4];
        /*species 9: CO2 */
        species[9] =
            -4.89669600e+04 * invT
            +4.40901890e+00
            -4.45362300e+00 * tc[0]
            -1.57008450e-03 * tc[1]
            +2.13068500e-07 * tc[2]
            -1.99499750e-11 * tc[3]
            +8.34516500e-16 * tc[4];
        /*species 10: HCO */
        species[10] =
            +3.91632400e+03 * invT
            -2.99502800e+00
            -3.55727100e+00 * tc[0]
            -1.67278650e-03 * tc[1]
            +2.22501000e-07 * tc[2]
            -2.05881083e-11 * tc[3]
            +8.56925500e-16 * tc[4];
        /*species 11: N2 */
        species[11] =
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


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cv_R(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
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
        /*species 5: H */
        species[5] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
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
        /*species 8: CO */
        species[8] =
            +2.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 9: CO2 */
        species[9] =
            +1.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
        /*species 10: HCO */
        species[10] =
            +1.89833000e+00
            +6.19914700e-03 * tc[1]
            -9.62308400e-06 * tc[2]
            +1.08982500e-08 * tc[3]
            -4.57488500e-12 * tc[4];
        /*species 11: N2 */
        species[11] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
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
        /*species 5: H */
        species[5] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
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
        /*species 8: CO */
        species[8] =
            +2.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 9: CO2 */
        species[9] =
            +3.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
        /*species 10: HCO */
        species[10] =
            +2.55727100e+00
            +3.34557300e-03 * tc[1]
            -1.33500600e-06 * tc[2]
            +2.47057300e-10 * tc[3]
            -1.71385100e-14 * tc[4];
        /*species 11: N2 */
        species[11] =
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
void cp_R(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +3.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
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
        /*species 5: H */
        species[5] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
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
        /*species 8: CO */
        species[8] =
            +3.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 9: CO2 */
        species[9] =
            +2.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
        /*species 10: HCO */
        species[10] =
            +2.89833000e+00
            +6.19914700e-03 * tc[1]
            -9.62308400e-06 * tc[2]
            +1.08982500e-08 * tc[3]
            -4.57488500e-12 * tc[4];
        /*species 11: N2 */
        species[11] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
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
        /*species 5: H */
        species[5] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
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
        /*species 8: CO */
        species[8] =
            +3.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 9: CO2 */
        species[9] =
            +4.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
        /*species 10: HCO */
        species[10] =
            +3.55727100e+00
            +3.34557300e-03 * tc[1]
            -1.33500600e-06 * tc[2]
            +2.47057300e-10 * tc[3]
            -1.71385100e-14 * tc[4];
        /*species 11: N2 */
        species[11] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesInternalEnergy(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +2.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
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
        /*species 5: H */
        species[5] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
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
        /*species 8: CO */
        species[8] =
            +2.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 9: CO2 */
        species[9] =
            +1.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
        /*species 10: HCO */
        species[10] =
            +1.89833000e+00
            +3.09957350e-03 * tc[1]
            -3.20769467e-06 * tc[2]
            +2.72456250e-09 * tc[3]
            -9.14977000e-13 * tc[4]
            +4.15992200e+03 * invT;
        /*species 11: N2 */
        species[11] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 1: O2 */
        species[1] =
            +2.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
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
        /*species 5: H */
        species[5] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
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
        /*species 8: CO */
        species[8] =
            +2.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 9: CO2 */
        species[9] =
            +3.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
        /*species 10: HCO */
        species[10] =
            +2.55727100e+00
            +1.67278650e-03 * tc[1]
            -4.45002000e-07 * tc[2]
            +6.17643250e-11 * tc[3]
            -3.42770200e-15 * tc[4]
            +3.91632400e+03 * invT;
        /*species 11: N2 */
        species[11] =
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
void speciesEnthalpy(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +3.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +3.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
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
        /*species 5: H */
        species[5] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
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
        /*species 8: CO */
        species[8] =
            +3.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 9: CO2 */
        species[9] =
            +2.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
        /*species 10: HCO */
        species[10] =
            +2.89833000e+00
            +3.09957350e-03 * tc[1]
            -3.20769467e-06 * tc[2]
            +2.72456250e-09 * tc[3]
            -9.14977000e-13 * tc[4]
            +4.15992200e+03 * invT;
        /*species 11: N2 */
        species[11] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 1: O2 */
        species[1] =
            +3.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
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
        /*species 5: H */
        species[5] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
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
        /*species 8: CO */
        species[8] =
            +3.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 9: CO2 */
        species[9] =
            +4.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
        /*species 10: HCO */
        species[10] =
            +3.55727100e+00
            +1.67278650e-03 * tc[1]
            -4.45002000e-07 * tc[2]
            +6.17643250e-11 * tc[3]
            -3.42770200e-15 * tc[4]
            +3.91632400e+03 * invT;
        /*species 11: N2 */
        species[11] =
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
void speciesEntropy(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +3.29812400e+00 * tc[0]
            +8.24944200e-04 * tc[1]
            -4.07150750e-07 * tc[2]
            -3.15847800e-11 * tc[3]
            +1.03371800e-13 * tc[4]
            -3.29409400e+00 ;
        /*species 1: O2 */
        species[1] =
            +3.21293600e+00 * tc[0]
            +1.12748600e-03 * tc[1]
            -2.87807500e-07 * tc[2]
            +4.37959000e-10 * tc[3]
            -2.19213850e-13 * tc[4]
            +6.03473800e+00 ;
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
        /*species 5: H */
        species[5] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
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
        /*species 8: CO */
        species[8] =
            +3.26245200e+00 * tc[0]
            +1.51194100e-03 * tc[1]
            -1.94087750e-06 * tc[2]
            +1.86064800e-09 * tc[3]
            -6.18737750e-13 * tc[4]
            +4.84889700e+00 ;
        /*species 9: CO2 */
        species[9] =
            +2.27572500e+00 * tc[0]
            +9.92207200e-03 * tc[1]
            -5.20455500e-06 * tc[2]
            +2.28889567e-09 * tc[3]
            -5.29320000e-13 * tc[4]
            +1.01884900e+01 ;
        /*species 10: HCO */
        species[10] =
            +2.89833000e+00 * tc[0]
            +6.19914700e-03 * tc[1]
            -4.81154200e-06 * tc[2]
            +3.63275000e-09 * tc[3]
            -1.14372125e-12 * tc[4]
            +8.98361400e+00 ;
        /*species 11: N2 */
        species[11] =
            +3.29867700e+00 * tc[0]
            +1.40824000e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213750e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142300e+00 * tc[0]
            +7.00064400e-04 * tc[1]
            -2.81691450e-08 * tc[2]
            -3.07719267e-12 * tc[3]
            +3.95688000e-16 * tc[4]
            -1.35511000e+00 ;
        /*species 1: O2 */
        species[1] =
            +3.69757800e+00 * tc[0]
            +6.13519700e-04 * tc[1]
            -6.29421000e-08 * tc[2]
            +5.91760333e-12 * tc[3]
            -2.84108750e-16 * tc[4]
            +3.18916600e+00 ;
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
        /*species 5: H */
        species[5] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
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
        /*species 8: CO */
        species[8] =
            +3.02507800e+00 * tc[0]
            +1.44268900e-03 * tc[1]
            -2.81541400e-07 * tc[2]
            +3.39527000e-11 * tc[3]
            -1.72773800e-15 * tc[4]
            +6.10821800e+00 ;
        /*species 9: CO2 */
        species[9] =
            +4.45362300e+00 * tc[0]
            +3.14016900e-03 * tc[1]
            -6.39205500e-07 * tc[2]
            +7.97999000e-11 * tc[3]
            -4.17258250e-15 * tc[4]
            -9.55395900e-01 ;
        /*species 10: HCO */
        species[10] =
            +3.55727100e+00 * tc[0]
            +3.34557300e-03 * tc[1]
            -6.67503000e-07 * tc[2]
            +8.23524333e-11 * tc[3]
            -4.28462750e-15 * tc[4]
            +5.55229900e+00 ;
        /*species 11: N2 */
        species[11] =
            +2.92664000e+00 * tc[0]
            +1.48797700e-03 * tc[1]
            -2.84238050e-07 * tc[2]
            +3.36568000e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double *  wt)
{
    wt[0] = 2.015940; /*H2 */
    wt[1] = 31.998800; /*O2 */
    wt[2] = 15.999400; /*O */
    wt[3] = 17.007370; /*OH */
    wt[4] = 18.015340; /*H2O */
    wt[5] = 1.007970; /*H */
    wt[6] = 33.006770; /*HO2 */
    wt[7] = 34.014740; /*H2O2 */
    wt[8] = 28.010550; /*CO */
    wt[9] = 44.009950; /*CO2 */
    wt[10] = 29.018520; /*HCO */
    wt[11] = 28.013400; /*N2 */

    return;
}


/*save atomic weights into array */
void atomicWeight(double *  awt)
{
    awt[0] = 12.011150; /*C */
    awt[1] = 1.007970; /*H */
    awt[2] = 15.999400; /*O */
    awt[3] = 14.006700; /*N */

    return;
}
/* get temperature given internal energy in mass units and mass fracs */
void GET_T_GIVEN_EY(double *  e, double *  y, double *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double ein  = *e;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double e1,emin,emax,cv,t1,dt;
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
void GET_T_GIVEN_HY(double *  h, double *  y, double *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double hin  = *h;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double h1,hmin,hmax,cp,t1,dt;
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


/*compute the critical parameters for each species */
void GET_CRITPARAMS(double *  Tci, double *  ai, double *  bi, double *  acentric_i)
{

    double   EPS[12];
    double   SIG[12];
    double    wt[12];
    double avogadro = 6.02214199e23;
    double boltzmann = 1.3806503e-16; //we work in CGS
    double Rcst = 83.144598; //in bar [CGS] !

    egtransetEPS(EPS);
    egtransetSIG(SIG);
    molecularWeight(wt);

    /*species 0: H2 */
    /*Imported from NIST */
    Tci[0] = 33.145000 ; 
    ai[0] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[0],2.0) / (pow(2.015880,2.0) * 12.964000); 
    bi[0] = 0.08664 * Rcst * Tci[0] / (2.015880 * 12.964000); 
    acentric_i[0] = -0.219000 ;

    /*species 1: O2 */
    /*Imported from NIST */
    Tci[1] = 154.581000 ; 
    ai[1] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[1],2.0) / (pow(31.998800,2.0) * 50.430466); 
    bi[1] = 0.08664 * Rcst * Tci[1] / (31.998800 * 50.430466); 
    acentric_i[1] = 0.022200 ;

    /*species 2: O */
    Tci[2] = 1.316 * EPS[2] ; 
    ai[2] = (5.55 * pow(avogadro,2.0) * EPS[2]*boltzmann * pow(1e-8*SIG[2],3.0) ) / (pow(wt[2],2.0)); 
    bi[2] = 0.855 * avogadro * pow(1e-8*SIG[2],3.0) / (wt[2]); 
    acentric_i[2] = 0.0 ;

    /*species 3: OH */
    Tci[3] = 1.316 * EPS[3] ; 
    ai[3] = (5.55 * pow(avogadro,2.0) * EPS[3]*boltzmann * pow(1e-8*SIG[3],3.0) ) / (pow(wt[3],2.0)); 
    bi[3] = 0.855 * avogadro * pow(1e-8*SIG[3],3.0) / (wt[3]); 
    acentric_i[3] = 0.0 ;

    /*species 4: H2O */
    /*Imported from NIST */
    Tci[4] = 647.096000 ; 
    ai[4] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[4],2.0) / (pow(18.015340,2.0) * 220.640000); 
    bi[4] = 0.08664 * Rcst * Tci[4] / (18.015340 * 220.640000); 
    acentric_i[4] = 0.344300 ;

    /*species 5: H */
    Tci[5] = 1.316 * EPS[5] ; 
    ai[5] = (5.55 * pow(avogadro,2.0) * EPS[5]*boltzmann * pow(1e-8*SIG[5],3.0) ) / (pow(wt[5],2.0)); 
    bi[5] = 0.855 * avogadro * pow(1e-8*SIG[5],3.0) / (wt[5]); 
    acentric_i[5] = 0.0 ;

    /*species 6: HO2 */
    Tci[6] = 1.316 * EPS[6] ; 
    ai[6] = (5.55 * pow(avogadro,2.0) * EPS[6]*boltzmann * pow(1e-8*SIG[6],3.0) ) / (pow(wt[6],2.0)); 
    bi[6] = 0.855 * avogadro * pow(1e-8*SIG[6],3.0) / (wt[6]); 
    acentric_i[6] = 0.0 ;

    /*species 7: H2O2 */
    Tci[7] = 1.316 * EPS[7] ; 
    ai[7] = (5.55 * pow(avogadro,2.0) * EPS[7]*boltzmann * pow(1e-8*SIG[7],3.0) ) / (pow(wt[7],2.0)); 
    bi[7] = 0.855 * avogadro * pow(1e-8*SIG[7],3.0) / (wt[7]); 
    acentric_i[7] = 0.0 ;

    /*species 8: CO */
    /*Imported from NIST */
    Tci[8] = 132.850000 ; 
    ai[8] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[8],2.0) / (pow(28.010000,2.0) * 34.940000); 
    bi[8] = 0.08664 * Rcst * Tci[8] / (28.010000 * 34.940000); 
    acentric_i[8] = 0.045000 ;

    /*species 9: CO2 */
    /*Imported from NIST */
    Tci[9] = 304.120000 ; 
    ai[9] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[9],2.0) / (pow(44.009950,2.0) * 73.740000); 
    bi[9] = 0.08664 * Rcst * Tci[9] / (44.009950 * 73.740000); 
    acentric_i[9] = 0.225000 ;

    /*species 10: HCO */
    Tci[10] = 1.316 * EPS[10] ; 
    ai[10] = (5.55 * pow(avogadro,2.0) * EPS[10]*boltzmann * pow(1e-8*SIG[10],3.0) ) / (pow(wt[10],2.0)); 
    bi[10] = 0.855 * avogadro * pow(1e-8*SIG[10],3.0) / (wt[10]); 
    acentric_i[10] = 0.0 ;

    /*species 11: N2 */
    /*Imported from NIST */
    Tci[11] = 126.192000 ; 
    ai[11] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[11],2.0) / (pow(28.013400,2.0) * 33.958000); 
    bi[11] = 0.08664 * Rcst * Tci[11] / (28.013400 * 33.958000); 
    acentric_i[11] = 0.037200 ;

    return;
}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENIMC EGTRANSETLENIMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENIMC egtransetlenimc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENIMC egtransetlenimc_
#endif
void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 50;}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 3156;}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNO EGTRANSETNO
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNO egtransetno
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNO egtransetno_
#endif
void egtransetNO(int* NO ) {
    *NO = 4;}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKK EGTRANSETKK
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKK egtransetkk
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKK egtransetkk_
#endif
void egtransetKK(int* KK ) {
    *KK = 12;}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLITE EGTRANSETNLITE
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLITE egtransetnlite
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLITE egtransetnlite_
#endif
void egtransetNLITE(int* NLITE ) {
    *NLITE = 2;}


/*Patm in ergs/cm3 */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPATM EGTRANSETPATM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPATM egtransetpatm
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPATM egtransetpatm_
#endif
void egtransetPATM(double* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetWT EGTRANSETWT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetWT egtransetwt
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetWT egtransetwt_
#endif
void egtransetWT(double* WT ) {
    WT[0] = 2.01594000E+00;
    WT[1] = 3.19988000E+01;
    WT[2] = 1.59994000E+01;
    WT[3] = 1.70073700E+01;
    WT[4] = 1.80153400E+01;
    WT[5] = 1.00797000E+00;
    WT[6] = 3.30067700E+01;
    WT[7] = 3.40147400E+01;
    WT[8] = 2.80105500E+01;
    WT[9] = 4.40099500E+01;
    WT[10] = 2.90185200E+01;
    WT[11] = 2.80134000E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
void egtransetEPS(double* EPS ) {
    EPS[1] = 1.07400000E+02;
    EPS[9] = 2.44000000E+02;
    EPS[8] = 9.81000000E+01;
    EPS[11] = 9.75300000E+01;
    EPS[3] = 8.00000000E+01;
    EPS[0] = 3.80000000E+01;
    EPS[4] = 5.72400000E+02;
    EPS[6] = 1.07400000E+02;
    EPS[5] = 1.45000000E+02;
    EPS[10] = 4.98000000E+02;
    EPS[2] = 8.00000000E+01;
    EPS[7] = 1.07400000E+02;
}


/*the lennard-jones collision diameter in Angstroms */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
void egtransetSIG(double* SIG ) {
    SIG[1] = 3.45800000E+00;
    SIG[9] = 3.76300000E+00;
    SIG[8] = 3.65000000E+00;
    SIG[11] = 3.62100000E+00;
    SIG[3] = 2.75000000E+00;
    SIG[0] = 2.92000000E+00;
    SIG[4] = 2.60500000E+00;
    SIG[6] = 3.45800000E+00;
    SIG[5] = 2.05000000E+00;
    SIG[10] = 3.59000000E+00;
    SIG[2] = 2.75000000E+00;
    SIG[7] = 3.45800000E+00;
}


/*the dipole moment in Debye */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetDIP EGTRANSETDIP
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetDIP egtransetdip
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetDIP egtransetdip_
#endif
void egtransetDIP(double* DIP ) {
    DIP[1] = 0.00000000E+00;
    DIP[9] = 0.00000000E+00;
    DIP[8] = 0.00000000E+00;
    DIP[11] = 0.00000000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[0] = 0.00000000E+00;
    DIP[4] = 1.84400000E+00;
    DIP[6] = 0.00000000E+00;
    DIP[5] = 0.00000000E+00;
    DIP[10] = 0.00000000E+00;
    DIP[2] = 0.00000000E+00;
    DIP[7] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPOL EGTRANSETPOL
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPOL egtransetpol
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPOL egtransetpol_
#endif
void egtransetPOL(double* POL ) {
    POL[1] = 1.60000000E+00;
    POL[9] = 2.65000000E+00;
    POL[8] = 1.95000000E+00;
    POL[11] = 1.76000000E+00;
    POL[3] = 0.00000000E+00;
    POL[0] = 7.90000000E-01;
    POL[4] = 0.00000000E+00;
    POL[6] = 0.00000000E+00;
    POL[5] = 0.00000000E+00;
    POL[10] = 0.00000000E+00;
    POL[2] = 0.00000000E+00;
    POL[7] = 0.00000000E+00;
}


/*the rotational relaxation collision number at 298 K */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetZROT EGTRANSETZROT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetZROT egtransetzrot
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetZROT egtransetzrot_
#endif
void egtransetZROT(double* ZROT ) {
    ZROT[1] = 3.80000000E+00;
    ZROT[9] = 2.10000000E+00;
    ZROT[8] = 1.80000000E+00;
    ZROT[11] = 4.00000000E+00;
    ZROT[3] = 0.00000000E+00;
    ZROT[0] = 2.80000000E+02;
    ZROT[4] = 4.00000000E+00;
    ZROT[6] = 1.00000000E+00;
    ZROT[5] = 0.00000000E+00;
    ZROT[10] = 0.00000000E+00;
    ZROT[2] = 0.00000000E+00;
    ZROT[7] = 3.80000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLIN EGTRANSETNLIN
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLIN egtransetnlin
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLIN egtransetnlin_
#endif
void egtransetNLIN(int* NLIN) {
    NLIN[1] = 1;
    NLIN[9] = 1;
    NLIN[8] = 1;
    NLIN[11] = 1;
    NLIN[3] = 1;
    NLIN[0] = 1;
    NLIN[4] = 2;
    NLIN[6] = 2;
    NLIN[5] = 0;
    NLIN[10] = 2;
    NLIN[2] = 0;
    NLIN[7] = 2;
}


/*Poly fits for the viscosities, dim NO*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFETA EGTRANSETCOFETA
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFETA egtransetcofeta
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFETA egtransetcofeta_
#endif
void egtransetCOFETA(double* COFETA) {
    COFETA[0] = -1.37549435E+01;
    COFETA[1] = 9.65530587E-01;
    COFETA[2] = -4.45720114E-02;
    COFETA[3] = 2.05871810E-03;
    COFETA[4] = -1.68118868E+01;
    COFETA[5] = 2.52362554E+00;
    COFETA[6] = -2.49309128E-01;
    COFETA[7] = 1.10211025E-02;
    COFETA[8] = -1.48001581E+01;
    COFETA[9] = 1.79491990E+00;
    COFETA[10] = -1.54008440E-01;
    COFETA[11] = 6.86719439E-03;
    COFETA[12] = -1.47696103E+01;
    COFETA[13] = 1.79491990E+00;
    COFETA[14] = -1.54008440E-01;
    COFETA[15] = 6.86719439E-03;
    COFETA[16] = -1.17770937E+01;
    COFETA[17] = -8.26742721E-01;
    COFETA[18] = 3.39009079E-01;
    COFETA[19] = -2.00674327E-02;
    COFETA[20] = -1.98744496E+01;
    COFETA[21] = 3.41660514E+00;
    COFETA[22] = -3.63206306E-01;
    COFETA[23] = 1.58671021E-02;
    COFETA[24] = -1.67963797E+01;
    COFETA[25] = 2.52362554E+00;
    COFETA[26] = -2.49309128E-01;
    COFETA[27] = 1.10211025E-02;
    COFETA[28] = -1.67813391E+01;
    COFETA[29] = 2.52362554E+00;
    COFETA[30] = -2.49309128E-01;
    COFETA[31] = 1.10211025E-02;
    COFETA[32] = -1.63031240E+01;
    COFETA[33] = 2.26143219E+00;
    COFETA[34] = -2.15114671E-01;
    COFETA[35] = 9.53461976E-03;
    COFETA[36] = -2.36749526E+01;
    COFETA[37] = 4.99775518E+00;
    COFETA[38] = -5.52687718E-01;
    COFETA[39] = 2.34353338E-02;
    COFETA[40] = -2.11306792E+01;
    COFETA[41] = 3.26961843E+00;
    COFETA[42] = -2.51355092E-01;
    COFETA[43] = 7.35605058E-03;
    COFETA[44] = -1.62526779E+01;
    COFETA[45] = 2.24839597E+00;
    COFETA[46] = -2.13428438E-01;
    COFETA[47] = 9.46192413E-03;
}


/*Poly fits for the conductivities, dim NO*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
void egtransetCOFLAM(double* COFLAM) {
    COFLAM[0] = 1.11035511E+01;
    COFLAM[1] = -1.31883912E+00;
    COFLAM[2] = 2.44042473E-01;
    COFLAM[3] = -8.99836359E-03;
    COFLAM[4] = -2.51296685E+00;
    COFLAM[5] = 3.15165687E+00;
    COFLAM[6] = -3.10007661E-01;
    COFLAM[7] = 1.34522321E-02;
    COFLAM[8] = 1.98513952E+00;
    COFLAM[9] = 1.79491990E+00;
    COFLAM[10] = -1.54008440E-01;
    COFLAM[11] = 6.86719439E-03;
    COFLAM[12] = 1.60618734E+01;
    COFLAM[13] = -4.10626869E+00;
    COFLAM[14] = 6.63571339E-01;
    COFLAM[15] = -2.97906324E-02;
    COFLAM[16] = 2.21730522E+01;
    COFLAM[17] = -8.46935675E+00;
    COFLAM[18] = 1.46153820E+00;
    COFLAM[19] = -7.29502441E-02;
    COFLAM[20] = -3.24539191E-01;
    COFLAM[21] = 3.41660514E+00;
    COFLAM[22] = -3.63206306E-01;
    COFLAM[23] = 1.58671021E-02;
    COFLAM[24] = 5.56023763E-01;
    COFLAM[25] = 1.59073590E+00;
    COFLAM[26] = -5.28053839E-02;
    COFLAM[27] = 4.07601571E-04;
    COFLAM[28] = 1.48801095E+00;
    COFLAM[29] = 1.06176238E+00;
    COFLAM[30] = 5.72195640E-02;
    COFLAM[31] = -6.38391491E-03;
    COFLAM[32] = 9.94279402E+00;
    COFLAM[33] = -2.29161875E+00;
    COFLAM[34] = 4.74393585E-01;
    COFLAM[35] = -2.40686534E-02;
    COFLAM[36] = -1.21375633E+01;
    COFLAM[37] = 6.23624427E+00;
    COFLAM[38] = -6.22471402E-01;
    COFLAM[39] = 2.30613281E-02;
    COFLAM[40] = 1.70861874E+00;
    COFLAM[41] = -1.65305270E-01;
    COFLAM[42] = 3.30690261E-01;
    COFLAM[43] = -2.30097345E-02;
    COFLAM[44] = 1.15507419E+01;
    COFLAM[45] = -2.91453917E+00;
    COFLAM[46] = 5.55045765E-01;
    COFLAM[47] = -2.75173485E-02;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFD EGTRANSETCOFD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFD egtransetcofd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFD egtransetcofd_
#endif
void egtransetCOFD(double* COFD) {
    COFD[0] = -1.02395222E+01;
    COFD[1] = 2.15403244E+00;
    COFD[2] = -6.97480266E-02;
    COFD[3] = 3.23666871E-03;
    COFD[4] = -1.15797750E+01;
    COFD[5] = 2.43235504E+00;
    COFD[6] = -1.02890179E-01;
    COFD[7] = 4.52903603E-03;
    COFD[8] = -1.06250182E+01;
    COFD[9] = 2.15849701E+00;
    COFD[10] = -6.53886401E-02;
    COFD[11] = 2.81453370E-03;
    COFD[12] = -1.06283453E+01;
    COFD[13] = 2.15849701E+00;
    COFD[14] = -6.53886401E-02;
    COFD[15] = 2.81453370E-03;
    COFD[16] = -1.68758926E+01;
    COFD[17] = 4.49460303E+00;
    COFD[18] = -3.64766132E-01;
    COFD[19] = 1.56457153E-02;
    COFD[20] = -1.11808682E+01;
    COFD[21] = 2.66936727E+00;
    COFD[22] = -1.34411514E-01;
    COFD[23] = 5.92957488E-03;
    COFD[24] = -1.15806808E+01;
    COFD[25] = 2.43235504E+00;
    COFD[26] = -1.02890179E-01;
    COFD[27] = 4.52903603E-03;
    COFD[28] = -1.15815344E+01;
    COFD[29] = 2.43235504E+00;
    COFD[30] = -1.02890179E-01;
    COFD[31] = 4.52903603E-03;
    COFD[32] = -1.13541075E+01;
    COFD[33] = 2.31999438E+00;
    COFD[34] = -8.75064804E-02;
    COFD[35] = 3.82656365E-03;
    COFD[36] = -1.35545239E+01;
    COFD[37] = 3.13878730E+00;
    COFD[38] = -1.94980335E-01;
    COFD[39] = 8.53744486E-03;
    COFD[40] = -1.57161204E+01;
    COFD[41] = 3.96062263E+00;
    COFD[42] = -2.98964970E-01;
    COFD[43] = 1.29322565E-02;
    COFD[44] = -1.13253458E+01;
    COFD[45] = 2.31195095E+00;
    COFD[46] = -8.63988037E-02;
    COFD[47] = 3.77573452E-03;
    COFD[48] = -1.15797750E+01;
    COFD[49] = 2.43235504E+00;
    COFD[50] = -1.02890179E-01;
    COFD[51] = 4.52903603E-03;
    COFD[52] = -1.53110708E+01;
    COFD[53] = 3.37317428E+00;
    COFD[54] = -2.24900439E-01;
    COFD[55] = 9.81228151E-03;
    COFD[56] = -1.40864894E+01;
    COFD[57] = 3.07458927E+00;
    COFD[58] = -1.86899591E-01;
    COFD[59] = 8.19829781E-03;
    COFD[60] = -1.41066459E+01;
    COFD[61] = 3.07458927E+00;
    COFD[62] = -1.86899591E-01;
    COFD[63] = 8.19829781E-03;
    COFD[64] = -2.10640014E+01;
    COFD[65] = 5.50980695E+00;
    COFD[66] = -4.78335488E-01;
    COFD[67] = 1.98515434E-02;
    COFD[68] = -1.43712864E+01;
    COFD[69] = 3.70920439E+00;
    COFD[70] = -2.67274113E-01;
    COFD[71] = 1.15967481E-02;
    COFD[72] = -1.53187643E+01;
    COFD[73] = 3.37317428E+00;
    COFD[74] = -2.24900439E-01;
    COFD[75] = 9.81228151E-03;
    COFD[76] = -1.53261114E+01;
    COFD[77] = 3.37317428E+00;
    COFD[78] = -2.24900439E-01;
    COFD[79] = 9.81228151E-03;
    COFD[80] = -1.50371784E+01;
    COFD[81] = 3.26249588E+00;
    COFD[82] = -2.10658287E-01;
    COFD[83] = 9.20032462E-03;
    COFD[84] = -1.81197354E+01;
    COFD[85] = 4.33684042E+00;
    COFD[86] = -3.44981265E-01;
    COFD[87] = 1.48142449E-02;
    COFD[88] = -2.05107487E+01;
    COFD[89] = 5.21473296E+00;
    COFD[90] = -4.48646311E-01;
    COFD[91] = 1.89013813E-02;
    COFD[92] = -1.50096240E+01;
    COFD[93] = 3.25515933E+00;
    COFD[94] = -2.09710110E-01;
    COFD[95] = 9.15941830E-03;
    COFD[96] = -1.06250182E+01;
    COFD[97] = 2.15849701E+00;
    COFD[98] = -6.53886401E-02;
    COFD[99] = 2.81453370E-03;
    COFD[100] = -1.40864894E+01;
    COFD[101] = 3.07458927E+00;
    COFD[102] = -1.86899591E-01;
    COFD[103] = 8.19829781E-03;
    COFD[104] = -1.29877365E+01;
    COFD[105] = 2.80841511E+00;
    COFD[106] = -1.52629888E-01;
    COFD[107] = 6.72604927E-03;
    COFD[108] = -1.30027772E+01;
    COFD[109] = 2.80841511E+00;
    COFD[110] = -1.52629888E-01;
    COFD[111] = 6.72604927E-03;
    COFD[112] = -1.91096797E+01;
    COFD[113] = 5.02608697E+00;
    COFD[114] = -4.26959993E-01;
    COFD[115] = 1.80709910E-02;
    COFD[116] = -1.31860117E+01;
    COFD[117] = 3.38003453E+00;
    COFD[118] = -2.25783856E-01;
    COFD[119] = 9.85028660E-03;
    COFD[120] = -1.40916052E+01;
    COFD[121] = 3.07458927E+00;
    COFD[122] = -1.86899591E-01;
    COFD[123] = 8.19829781E-03;
    COFD[124] = -1.40964661E+01;
    COFD[125] = 3.07458927E+00;
    COFD[126] = -1.86899591E-01;
    COFD[127] = 8.19829781E-03;
    COFD[128] = -1.39007410E+01;
    COFD[129] = 2.99164244E+00;
    COFD[130] = -1.76293106E-01;
    COFD[131] = 7.74575100E-03;
    COFD[132] = -1.67115577E+01;
    COFD[133] = 3.98859394E+00;
    COFD[134] = -3.02316219E-01;
    COFD[135] = 1.30661099E-02;
    COFD[136] = -1.91045619E+01;
    COFD[137] = 4.87977047E+00;
    COFD[138] = -4.10448693E-01;
    COFD[139] = 1.74535827E-02;
    COFD[140] = -1.38756407E+01;
    COFD[141] = 2.98558426E+00;
    COFD[142] = -1.75507216E-01;
    COFD[143] = 7.71173691E-03;
    COFD[144] = -1.06283453E+01;
    COFD[145] = 2.15849701E+00;
    COFD[146] = -6.53886401E-02;
    COFD[147] = 2.81453370E-03;
    COFD[148] = -1.41066459E+01;
    COFD[149] = 3.07458927E+00;
    COFD[150] = -1.86899591E-01;
    COFD[151] = 8.19829781E-03;
    COFD[152] = -1.30027772E+01;
    COFD[153] = 2.80841511E+00;
    COFD[154] = -1.52629888E-01;
    COFD[155] = 6.72604927E-03;
    COFD[156] = -1.30182843E+01;
    COFD[157] = 2.80841511E+00;
    COFD[158] = -1.52629888E-01;
    COFD[159] = 6.72604927E-03;
    COFD[160] = -1.91256261E+01;
    COFD[161] = 5.02608697E+00;
    COFD[162] = -4.26959993E-01;
    COFD[163] = 1.80709910E-02;
    COFD[164] = -1.31877711E+01;
    COFD[165] = 3.38003453E+00;
    COFD[166] = -2.25783856E-01;
    COFD[167] = 9.85028660E-03;
    COFD[168] = -1.41119732E+01;
    COFD[169] = 3.07458927E+00;
    COFD[170] = -1.86899591E-01;
    COFD[171] = 8.19829781E-03;
    COFD[172] = -1.41170372E+01;
    COFD[173] = 3.07458927E+00;
    COFD[174] = -1.86899591E-01;
    COFD[175] = 8.19829781E-03;
    COFD[176] = -1.39199663E+01;
    COFD[177] = 2.99164244E+00;
    COFD[178] = -1.76293106E-01;
    COFD[179] = 7.74575100E-03;
    COFD[180] = -1.67337768E+01;
    COFD[181] = 3.98859394E+00;
    COFD[182] = -3.02316219E-01;
    COFD[183] = 1.30661099E-02;
    COFD[184] = -1.91240379E+01;
    COFD[185] = 4.87977047E+00;
    COFD[186] = -4.10448693E-01;
    COFD[187] = 1.74535827E-02;
    COFD[188] = -1.38948667E+01;
    COFD[189] = 2.98558426E+00;
    COFD[190] = -1.75507216E-01;
    COFD[191] = 7.71173691E-03;
    COFD[192] = -1.68758926E+01;
    COFD[193] = 4.49460303E+00;
    COFD[194] = -3.64766132E-01;
    COFD[195] = 1.56457153E-02;
    COFD[196] = -2.10640014E+01;
    COFD[197] = 5.50980695E+00;
    COFD[198] = -4.78335488E-01;
    COFD[199] = 1.98515434E-02;
    COFD[200] = -1.91096797E+01;
    COFD[201] = 5.02608697E+00;
    COFD[202] = -4.26959993E-01;
    COFD[203] = 1.80709910E-02;
    COFD[204] = -1.91256261E+01;
    COFD[205] = 5.02608697E+00;
    COFD[206] = -4.26959993E-01;
    COFD[207] = 1.80709910E-02;
    COFD[208] = -1.31492641E+01;
    COFD[209] = 1.48004311E+00;
    COFD[210] = 1.60499553E-01;
    COFD[211] = -1.19765679E-02;
    COFD[212] = -1.93611051E+01;
    COFD[213] = 5.51579726E+00;
    COFD[214] = -4.76061961E-01;
    COFD[215] = 1.96329391E-02;
    COFD[216] = -2.04177482E+01;
    COFD[217] = 5.31457079E+00;
    COFD[218] = -4.58216496E-01;
    COFD[219] = 1.91825910E-02;
    COFD[220] = -2.04230073E+01;
    COFD[221] = 5.31457079E+00;
    COFD[222] = -4.58216496E-01;
    COFD[223] = 1.91825910E-02;
    COFD[224] = -2.08943798E+01;
    COFD[225] = 5.44718652E+00;
    COFD[226] = -4.72082953E-01;
    COFD[227] = 1.96531321E-02;
    COFD[228] = -2.12021508E+01;
    COFD[229] = 5.20775052E+00;
    COFD[230] = -4.07348327E-01;
    COFD[231] = 1.55473283E-02;
    COFD[232] = -1.87171417E+01;
    COFD[233] = 4.00967621E+00;
    COFD[234] = -2.21153539E-01;
    COFD[235] = 6.31528745E-03;
    COFD[236] = -2.08123325E+01;
    COFD[237] = 5.42470154E+00;
    COFD[238] = -4.69700416E-01;
    COFD[239] = 1.95706904E-02;
    COFD[240] = -1.11808682E+01;
    COFD[241] = 2.66936727E+00;
    COFD[242] = -1.34411514E-01;
    COFD[243] = 5.92957488E-03;
    COFD[244] = -1.43712864E+01;
    COFD[245] = 3.70920439E+00;
    COFD[246] = -2.67274113E-01;
    COFD[247] = 1.15967481E-02;
    COFD[248] = -1.31860117E+01;
    COFD[249] = 3.38003453E+00;
    COFD[250] = -2.25783856E-01;
    COFD[251] = 9.85028660E-03;
    COFD[252] = -1.31877711E+01;
    COFD[253] = 3.38003453E+00;
    COFD[254] = -2.25783856E-01;
    COFD[255] = 9.85028660E-03;
    COFD[256] = -1.93611051E+01;
    COFD[257] = 5.51579726E+00;
    COFD[258] = -4.76061961E-01;
    COFD[259] = 1.96329391E-02;
    COFD[260] = -1.43693056E+01;
    COFD[261] = 4.03992999E+00;
    COFD[262] = -3.08044800E-01;
    COFD[263] = 1.32757775E-02;
    COFD[264] = -1.43717529E+01;
    COFD[265] = 3.70920439E+00;
    COFD[266] = -2.67274113E-01;
    COFD[267] = 1.15967481E-02;
    COFD[268] = -1.43721922E+01;
    COFD[269] = 3.70920439E+00;
    COFD[270] = -2.67274113E-01;
    COFD[271] = 1.15967481E-02;
    COFD[272] = -1.40524065E+01;
    COFD[273] = 3.56261348E+00;
    COFD[274] = -2.48287981E-01;
    COFD[275] = 1.07752947E-02;
    COFD[276] = -1.72993972E+01;
    COFD[277] = 4.71931868E+00;
    COFD[278] = -3.91258152E-01;
    COFD[279] = 1.66866639E-02;
    COFD[280] = -1.95312724E+01;
    COFD[281] = 5.47046983E+00;
    COFD[282] = -4.74577605E-01;
    COFD[283] = 1.97408822E-02;
    COFD[284] = -1.40298830E+01;
    COFD[285] = 3.55837688E+00;
    COFD[286] = -2.47785790E-01;
    COFD[287] = 1.07555332E-02;
    COFD[288] = -1.15806808E+01;
    COFD[289] = 2.43235504E+00;
    COFD[290] = -1.02890179E-01;
    COFD[291] = 4.52903603E-03;
    COFD[292] = -1.53187643E+01;
    COFD[293] = 3.37317428E+00;
    COFD[294] = -2.24900439E-01;
    COFD[295] = 9.81228151E-03;
    COFD[296] = -1.40916052E+01;
    COFD[297] = 3.07458927E+00;
    COFD[298] = -1.86899591E-01;
    COFD[299] = 8.19829781E-03;
    COFD[300] = -1.41119732E+01;
    COFD[301] = 3.07458927E+00;
    COFD[302] = -1.86899591E-01;
    COFD[303] = 8.19829781E-03;
    COFD[304] = -2.04177482E+01;
    COFD[305] = 5.31457079E+00;
    COFD[306] = -4.58216496E-01;
    COFD[307] = 1.91825910E-02;
    COFD[308] = -1.43717529E+01;
    COFD[309] = 3.70920439E+00;
    COFD[310] = -2.67274113E-01;
    COFD[311] = 1.15967481E-02;
    COFD[312] = -1.53265780E+01;
    COFD[313] = 3.37317428E+00;
    COFD[314] = -2.24900439E-01;
    COFD[315] = 9.81228151E-03;
    COFD[316] = -1.53340417E+01;
    COFD[317] = 3.37317428E+00;
    COFD[318] = -2.24900439E-01;
    COFD[319] = 9.81228151E-03;
    COFD[320] = -1.50443569E+01;
    COFD[321] = 3.26249588E+00;
    COFD[322] = -2.10658287E-01;
    COFD[323] = 9.20032462E-03;
    COFD[324] = -1.81286555E+01;
    COFD[325] = 4.33684042E+00;
    COFD[326] = -3.44981265E-01;
    COFD[327] = 1.48142449E-02;
    COFD[328] = -2.05180636E+01;
    COFD[329] = 5.21473296E+00;
    COFD[330] = -4.48646311E-01;
    COFD[331] = 1.89013813E-02;
    COFD[332] = -1.50168028E+01;
    COFD[333] = 3.25515933E+00;
    COFD[334] = -2.09710110E-01;
    COFD[335] = 9.15941830E-03;
    COFD[336] = -1.15815344E+01;
    COFD[337] = 2.43235504E+00;
    COFD[338] = -1.02890179E-01;
    COFD[339] = 4.52903603E-03;
    COFD[340] = -1.53261114E+01;
    COFD[341] = 3.37317428E+00;
    COFD[342] = -2.24900439E-01;
    COFD[343] = 9.81228151E-03;
    COFD[344] = -1.40964661E+01;
    COFD[345] = 3.07458927E+00;
    COFD[346] = -1.86899591E-01;
    COFD[347] = 8.19829781E-03;
    COFD[348] = -1.41170372E+01;
    COFD[349] = 3.07458927E+00;
    COFD[350] = -1.86899591E-01;
    COFD[351] = 8.19829781E-03;
    COFD[352] = -2.04230073E+01;
    COFD[353] = 5.31457079E+00;
    COFD[354] = -4.58216496E-01;
    COFD[355] = 1.91825910E-02;
    COFD[356] = -1.43721922E+01;
    COFD[357] = 3.70920439E+00;
    COFD[358] = -2.67274113E-01;
    COFD[359] = 1.15967481E-02;
    COFD[360] = -1.53340417E+01;
    COFD[361] = 3.37317428E+00;
    COFD[362] = -2.24900439E-01;
    COFD[363] = 9.81228151E-03;
    COFD[364] = -1.53416186E+01;
    COFD[365] = 3.37317428E+00;
    COFD[366] = -2.24900439E-01;
    COFD[367] = 9.81228151E-03;
    COFD[368] = -1.50512053E+01;
    COFD[369] = 3.26249588E+00;
    COFD[370] = -2.10658287E-01;
    COFD[371] = 9.20032462E-03;
    COFD[372] = -1.81371948E+01;
    COFD[373] = 4.33684042E+00;
    COFD[374] = -3.44981265E-01;
    COFD[375] = 1.48142449E-02;
    COFD[376] = -2.05250441E+01;
    COFD[377] = 5.21473296E+00;
    COFD[378] = -4.48646311E-01;
    COFD[379] = 1.89013813E-02;
    COFD[380] = -1.50236516E+01;
    COFD[381] = 3.25515933E+00;
    COFD[382] = -2.09710110E-01;
    COFD[383] = 9.15941830E-03;
    COFD[384] = -1.13541075E+01;
    COFD[385] = 2.31999438E+00;
    COFD[386] = -8.75064804E-02;
    COFD[387] = 3.82656365E-03;
    COFD[388] = -1.50371784E+01;
    COFD[389] = 3.26249588E+00;
    COFD[390] = -2.10658287E-01;
    COFD[391] = 9.20032462E-03;
    COFD[392] = -1.39007410E+01;
    COFD[393] = 2.99164244E+00;
    COFD[394] = -1.76293106E-01;
    COFD[395] = 7.74575100E-03;
    COFD[396] = -1.39199663E+01;
    COFD[397] = 2.99164244E+00;
    COFD[398] = -1.76293106E-01;
    COFD[399] = 7.74575100E-03;
    COFD[400] = -2.08943798E+01;
    COFD[401] = 5.44718652E+00;
    COFD[402] = -4.72082953E-01;
    COFD[403] = 1.96531321E-02;
    COFD[404] = -1.40524065E+01;
    COFD[405] = 3.56261348E+00;
    COFD[406] = -2.48287981E-01;
    COFD[407] = 1.07752947E-02;
    COFD[408] = -1.50443569E+01;
    COFD[409] = 3.26249588E+00;
    COFD[410] = -2.10658287E-01;
    COFD[411] = 9.20032462E-03;
    COFD[412] = -1.50512053E+01;
    COFD[413] = 3.26249588E+00;
    COFD[414] = -2.10658287E-01;
    COFD[415] = 9.20032462E-03;
    COFD[416] = -1.48061490E+01;
    COFD[417] = 3.16912473E+00;
    COFD[418] = -1.98792456E-01;
    COFD[419] = 8.69726395E-03;
    COFD[420] = -1.77673000E+01;
    COFD[421] = 4.20234040E+00;
    COFD[422] = -3.28057658E-01;
    COFD[423] = 1.41006192E-02;
    COFD[424] = -2.02361862E+01;
    COFD[425] = 5.11785645E+00;
    COFD[426] = -4.37867828E-01;
    COFD[427] = 1.85047543E-02;
    COFD[428] = -1.47850486E+01;
    COFD[429] = 3.16433919E+00;
    COFD[430] = -1.98191564E-01;
    COFD[431] = 8.67209742E-03;
    COFD[432] = -1.35545239E+01;
    COFD[433] = 3.13878730E+00;
    COFD[434] = -1.94980335E-01;
    COFD[435] = 8.53744486E-03;
    COFD[436] = -1.81197354E+01;
    COFD[437] = 4.33684042E+00;
    COFD[438] = -3.44981265E-01;
    COFD[439] = 1.48142449E-02;
    COFD[440] = -1.67115577E+01;
    COFD[441] = 3.98859394E+00;
    COFD[442] = -3.02316219E-01;
    COFD[443] = 1.30661099E-02;
    COFD[444] = -1.67337768E+01;
    COFD[445] = 3.98859394E+00;
    COFD[446] = -3.02316219E-01;
    COFD[447] = 1.30661099E-02;
    COFD[448] = -2.12021508E+01;
    COFD[449] = 5.20775052E+00;
    COFD[450] = -4.07348327E-01;
    COFD[451] = 1.55473283E-02;
    COFD[452] = -1.72993972E+01;
    COFD[453] = 4.71931868E+00;
    COFD[454] = -3.91258152E-01;
    COFD[455] = 1.66866639E-02;
    COFD[456] = -1.81286555E+01;
    COFD[457] = 4.33684042E+00;
    COFD[458] = -3.44981265E-01;
    COFD[459] = 1.48142449E-02;
    COFD[460] = -1.81371948E+01;
    COFD[461] = 4.33684042E+00;
    COFD[462] = -3.44981265E-01;
    COFD[463] = 1.48142449E-02;
    COFD[464] = -1.77673000E+01;
    COFD[465] = 4.20234040E+00;
    COFD[466] = -3.28057658E-01;
    COFD[467] = 1.41006192E-02;
    COFD[468] = -2.10907727E+01;
    COFD[469] = 5.29211327E+00;
    COFD[470] = -4.56068366E-01;
    COFD[471] = 1.91195062E-02;
    COFD[472] = -2.20758752E+01;
    COFD[473] = 5.52171573E+00;
    COFD[474] = -4.63284984E-01;
    COFD[475] = 1.85570924E-02;
    COFD[476] = -1.77350592E+01;
    COFD[477] = 4.19328271E+00;
    COFD[478] = -3.26911461E-01;
    COFD[479] = 1.40520357E-02;
    COFD[480] = -1.57161204E+01;
    COFD[481] = 3.96062263E+00;
    COFD[482] = -2.98964970E-01;
    COFD[483] = 1.29322565E-02;
    COFD[484] = -2.05107487E+01;
    COFD[485] = 5.21473296E+00;
    COFD[486] = -4.48646311E-01;
    COFD[487] = 1.89013813E-02;
    COFD[488] = -1.91045619E+01;
    COFD[489] = 4.87977047E+00;
    COFD[490] = -4.10448693E-01;
    COFD[491] = 1.74535827E-02;
    COFD[492] = -1.91240379E+01;
    COFD[493] = 4.87977047E+00;
    COFD[494] = -4.10448693E-01;
    COFD[495] = 1.74535827E-02;
    COFD[496] = -1.87171417E+01;
    COFD[497] = 4.00967621E+00;
    COFD[498] = -2.21153539E-01;
    COFD[499] = 6.31528745E-03;
    COFD[500] = -1.95312724E+01;
    COFD[501] = 5.47046983E+00;
    COFD[502] = -4.74577605E-01;
    COFD[503] = 1.97408822E-02;
    COFD[504] = -2.05180636E+01;
    COFD[505] = 5.21473296E+00;
    COFD[506] = -4.48646311E-01;
    COFD[507] = 1.89013813E-02;
    COFD[508] = -2.05250441E+01;
    COFD[509] = 5.21473296E+00;
    COFD[510] = -4.48646311E-01;
    COFD[511] = 1.89013813E-02;
    COFD[512] = -2.02361862E+01;
    COFD[513] = 5.11785645E+00;
    COFD[514] = -4.37867828E-01;
    COFD[515] = 1.85047543E-02;
    COFD[516] = -2.20758752E+01;
    COFD[517] = 5.52171573E+00;
    COFD[518] = -4.63284984E-01;
    COFD[519] = 1.85570924E-02;
    COFD[520] = -1.98983761E+01;
    COFD[521] = 4.38041133E+00;
    COFD[522] = -2.77538214E-01;
    COFD[523] = 9.06748822E-03;
    COFD[524] = -2.02052895E+01;
    COFD[525] = 5.10993120E+00;
    COFD[526] = -4.36931630E-01;
    COFD[527] = 1.84677592E-02;
    COFD[528] = -1.13253458E+01;
    COFD[529] = 2.31195095E+00;
    COFD[530] = -8.63988037E-02;
    COFD[531] = 3.77573452E-03;
    COFD[532] = -1.50096240E+01;
    COFD[533] = 3.25515933E+00;
    COFD[534] = -2.09710110E-01;
    COFD[535] = 9.15941830E-03;
    COFD[536] = -1.38756407E+01;
    COFD[537] = 2.98558426E+00;
    COFD[538] = -1.75507216E-01;
    COFD[539] = 7.71173691E-03;
    COFD[540] = -1.38948667E+01;
    COFD[541] = 2.98558426E+00;
    COFD[542] = -1.75507216E-01;
    COFD[543] = 7.71173691E-03;
    COFD[544] = -2.08123325E+01;
    COFD[545] = 5.42470154E+00;
    COFD[546] = -4.69700416E-01;
    COFD[547] = 1.95706904E-02;
    COFD[548] = -1.40298830E+01;
    COFD[549] = 3.55837688E+00;
    COFD[550] = -2.47785790E-01;
    COFD[551] = 1.07555332E-02;
    COFD[552] = -1.50168028E+01;
    COFD[553] = 3.25515933E+00;
    COFD[554] = -2.09710110E-01;
    COFD[555] = 9.15941830E-03;
    COFD[556] = -1.50236516E+01;
    COFD[557] = 3.25515933E+00;
    COFD[558] = -2.09710110E-01;
    COFD[559] = 9.15941830E-03;
    COFD[560] = -1.47850486E+01;
    COFD[561] = 3.16433919E+00;
    COFD[562] = -1.98191564E-01;
    COFD[563] = 8.67209742E-03;
    COFD[564] = -1.77350592E+01;
    COFD[565] = 4.19328271E+00;
    COFD[566] = -3.26911461E-01;
    COFD[567] = 1.40520357E-02;
    COFD[568] = -2.02052895E+01;
    COFD[569] = 5.10993120E+00;
    COFD[570] = -4.36931630E-01;
    COFD[571] = 1.84677592E-02;
    COFD[572] = -1.47639290E+01;
    COFD[573] = 3.15955654E+00;
    COFD[574] = -1.97590757E-01;
    COFD[575] = 8.64692156E-03;
}


/*List of specs with small weight, dim NLITE */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKTDIF EGTRANSETKTDIF
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKTDIF egtransetktdif
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKTDIF egtransetktdif_
#endif
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 1;
    KTDIF[1] = 6;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFTD EGTRANSETCOFTD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFTD egtransetcoftd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFTD egtransetcoftd_
#endif
void egtransetCOFTD(double* COFTD) {
    COFTD[0] = 0.00000000E+00;
    COFTD[1] = 0.00000000E+00;
    COFTD[2] = 0.00000000E+00;
    COFTD[3] = 0.00000000E+00;
    COFTD[4] = 4.42739084E-01;
    COFTD[5] = 7.11770818E-05;
    COFTD[6] = -3.84768062E-08;
    COFTD[7] = 6.86323437E-12;
    COFTD[8] = 4.15583337E-01;
    COFTD[9] = 1.09738399E-05;
    COFTD[10] = -3.96021963E-09;
    COFTD[11] = 1.14414443E-12;
    COFTD[12] = 4.21932443E-01;
    COFTD[13] = 1.11414935E-05;
    COFTD[14] = -4.02072219E-09;
    COFTD[15] = 1.16162418E-12;
    COFTD[16] = 6.02028221E-02;
    COFTD[17] = 5.61561867E-04;
    COFTD[18] = -2.55372862E-07;
    COFTD[19] = 3.63389913E-11;
    COFTD[20] = -1.52534742E-01;
    COFTD[21] = -5.46404022E-05;
    COFTD[22] = 2.93412470E-08;
    COFTD[23] = -4.87091914E-12;
    COFTD[24] = 4.44452569E-01;
    COFTD[25] = 7.14525507E-05;
    COFTD[26] = -3.86257187E-08;
    COFTD[27] = 6.88979640E-12;
    COFTD[28] = 4.46070183E-01;
    COFTD[29] = 7.17126069E-05;
    COFTD[30] = -3.87662996E-08;
    COFTD[31] = 6.91487226E-12;
    COFTD[32] = 4.44653617E-01;
    COFTD[33] = 5.06631704E-05;
    COFTD[34] = -2.69820900E-08;
    COFTD[35] = 5.01289759E-12;
    COFTD[36] = 3.25742450E-01;
    COFTD[37] = 3.03633411E-04;
    COFTD[38] = -1.55290330E-07;
    COFTD[39] = 2.41466436E-11;
    COFTD[40] = 1.61392287E-01;
    COFTD[41] = 5.01084129E-04;
    COFTD[42] = -2.38273963E-07;
    COFTD[43] = 3.49344424E-11;
    COFTD[44] = 4.45261966E-01;
    COFTD[45] = 4.94697174E-05;
    COFTD[46] = -2.63023442E-08;
    COFTD[47] = 4.90306217E-12;
    COFTD[48] = 1.52534742E-01;
    COFTD[49] = 5.46404022E-05;
    COFTD[50] = -2.93412470E-08;
    COFTD[51] = 4.87091914E-12;
    COFTD[52] = 2.20482843E-01;
    COFTD[53] = 4.80164288E-04;
    COFTD[54] = -2.32927944E-07;
    COFTD[55] = 3.46470436E-11;
    COFTD[56] = 2.70010150E-01;
    COFTD[57] = 3.61555093E-04;
    COFTD[58] = -1.80744752E-07;
    COFTD[59] = 2.75321248E-11;
    COFTD[60] = 2.72041664E-01;
    COFTD[61] = 3.64275376E-04;
    COFTD[62] = -1.82104647E-07;
    COFTD[63] = 2.77392722E-11;
    COFTD[64] = -1.41883744E-01;
    COFTD[65] = 7.66558810E-04;
    COFTD[66] = -3.06550003E-07;
    COFTD[67] = 4.02959502E-11;
    COFTD[68] = 0.00000000E+00;
    COFTD[69] = 0.00000000E+00;
    COFTD[70] = 0.00000000E+00;
    COFTD[71] = 0.00000000E+00;
    COFTD[72] = 2.20907853E-01;
    COFTD[73] = 4.81089870E-04;
    COFTD[74] = -2.33376944E-07;
    COFTD[75] = 3.47138305E-11;
    COFTD[76] = 2.21308399E-01;
    COFTD[77] = 4.81962174E-04;
    COFTD[78] = -2.33800100E-07;
    COFTD[79] = 3.47767730E-11;
    COFTD[80] = 2.39409939E-01;
    COFTD[81] = 4.47197179E-04;
    COFTD[82] = -2.18951702E-07;
    COFTD[83] = 3.27973510E-11;
    COFTD[84] = 2.44369385E-02;
    COFTD[85] = 7.18242498E-04;
    COFTD[86] = -3.19718504E-07;
    COFTD[87] = 4.48828685E-11;
    COFTD[88] = -1.24647991E-01;
    COFTD[89] = 7.96525614E-04;
    COFTD[90] = -3.24998782E-07;
    COFTD[91] = 4.32516736E-11;
    COFTD[92] = 2.40744421E-01;
    COFTD[93] = 4.45343451E-04;
    COFTD[94] = -2.18173874E-07;
    COFTD[95] = 3.26958506E-11;
}

}

/* End of file  */




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
! Mechanism for CO/H2 oxidation by Hemanth Kolla
! based on the skeletal mechanism of Evatt R. Hawkes and Ramanan Sankaran
!
! Only difference from skeletal mechanism is the inclusion of H2O2 and
! the reactions corresponding to its formation and consumption
! 
! The skeletal mechanism is derived from a detailed C1 mechanism by J. Li
!
! Original notes on the skeletal mechanism and the detailed mechanisms follow
!
! ------------------------------------------------------------------------------
! Skeletal Mechanism for CO/H2 oxidation
! by Evatt R. Hawkes and Ramanan Sankaran
!
! Reduced from complete C1 Mechanism that is published in:
!
! J. Li, PhD Thesis, 
! Mechanical and Aerospace Engineering Department, 
! Princeton University, Princeton NJ.  November 2004. Thesis No. 3122-T.
!
! http://www.princeton.edu/~combust/database/files/symposium/C1_Mechanism.zip
!
! At the time of writing, a publication to IJCK is in preparation 
! by the authors of the complete mechanism.
!
! This mechanism was reduced specifically for the purpose of the
! Direct Numerical Simulations performed in
! Hawkes, E.R., Sankaran, R., Sutherland, J.C., and Chen, J.H. (2006)
! Proc. Combust. Inst. 31, to appear.
!
! It was validated by comparison with the full mechanism in several cases:
! 1. freely propagating premixed flames in a range of equivalence ratios,
! 2. opposed-flow non-premixed flames in a range of strains up to extinction,
! 3. homogeneous ignition calculations for mixtures of fuel and oxidizer streams 
!    and equilibrium products,
! 4. two-dimensional DNS featuring extinction and local reignition.
! In all cases the agreement was excellent.
!
! However, the mechanism is validated ONLY for the specific conditions of the 
! DNS and is not expected to be valid in general.
!
! The following changes (only) were made to the complete mechanism:
! 1) Only the species H2 O2 O OH H2O H HO2 CO CO2 HCO N2 were retained. 
!    All other species and reactions involving these species were removed.  Note
!    this includes all C containing species other than those essential for CO 
!    oxidation: CO, CO2 and HCO.  For the atmospheric pressure of the simulation,
!    H2O2 was also found to be unimportant and was removed.
! 2) It was found HCO has only a minor importance, and its reaction rates were 
!    dominated by a few key reactions.  These reactions (below) were retained
!    and all others neglected.
!
! Steady state assumptions were investigated and found to be reasonably accurate
! for fully burning conditions but it was found they increased stiffness, hence 
! they were not employed.  For reference, a steady state approximation for HCO
! and HO2 may be expected to perform reasonably well if the added stiffness can
! be tackled.  However, note the HO2 steady state assumption will degrade the 
! prediction of ignition at longer ignition delay times.
!
! ---------------------------------------------------------------------------------
!
! Notes on the original mechanism by its authors:
!
! Authors: J. Li, Z. Zhao, A. Kazakov, F.L. Dryer, 
! Address: Dept. of Mechanical and Aerospace Engineering, 
! Princeton University, Princeton, NJ 08544
!
! This C1 mechanism is based on the CH3OH mechanism of Held and Dryer (IJCK,1998, 30, 805)
! with following important revision:
! 1.  H2/O2 subset is updated to that of Li et al. (IJCK, in press, 2004)
! 2.  CO + OH = CO2 + H is optimized to fit the literature experimental result 
! 3.  HCO + M = H + CO + M is optimized to fit the literature experimental result
! 4.  CH3 + HO2 = CH3O + OH is modified to match Scire's value at 1000 K 
! 5.  CH3 + HO2 = CH4 + H is taken from Scire (IJCK, 2001, 33, 75)
! 6.  CH3 + O2 = CH2O + OH is taken from Scire (2002, Ph.D. thesis)
! 7.  CH2O + HO2 = HCO + H2O2 is from Eiteneer et al. (JPC A. 1998, 102, 5196)
! 8.  CH2O + H = HCO + H2 is from Irdam et al. (IJCK 1993, 25, 285)
! 9.  CH2O + M reactions are from Friedrichs et al.(IJCK 2004, 36, 157)
! 10. CH3OH decomposition reactions are taken from GRI-3.0 (1999)
! 11. CH2OH + HCO = CH2O + CH2O is taken from Friedrichs et al. (IJCK, 2004, 36, 157)
! 12. CH2OH + HCO = CH3OH + CO is changed to keep the branching ratio with the above reaction
! 13. HCOOH reactions are not included since it is not important and has large uncertainties
! 14. dHf of OH is adjusted to 8.91 kcal/mol (Ruscic et al. JPC A. 2002, 106, 2727)
! 15. thermochemical data of CH2OH is fit to Johnson & Hudgens' table (JPC 1996, 100, 19874)


ELEMENTS
C H O N
END

SPECIES
H2 O2 O OH H2O H HO2 H2O2 CO CO2 HCO N2 
END

TRANS ALL
H2                 1    38.000     2.920     0.000     0.790   280.000          
O2                 1   107.400     3.458     0.000     1.600     3.800          
O                  0    80.000     2.750     0.000     0.000     0.000          
OH                 1    80.000     2.750     0.000     0.000     0.000          
H2O                2   572.400     2.605     1.844     0.000     4.000          
H                  0   145.000     2.050     0.000     0.000     0.000          
HO2                2   107.400     3.458     0.000     0.000     1.000          
H2O2               2   107.400     3.458     0.000     0.000     3.800          
CO                 1    98.100     3.650     0.000     1.950     1.800          
CO2                1   244.000     3.763     0.000     2.650     2.100          
HCO                2   498.000     3.590     0.000     0.000     0.000          
N2                 1    97.530     3.621     0.000     1.760     4.000          
END



REACTIONS

! ************ H2-O2 Chain Reactions **********************

! Hessler, J. Phys. Chem. A, 102:4517 (1998)
H+O2=O+OH                 3.547e+15 -0.406  1.6599E+4

! Sutherland et al., 21st Symposium, p. 929 (1986)
O+H2=H+OH                 0.508E+05  2.67  0.629E+04

! Michael and Sutherland, J. Phys. Chem. 92:3853 (1988)
H2+OH=H2O+H               0.216E+09  1.51  0.343E+04

! Sutherland et al., 23rd Symposium, p. 51 (1990)
O+H2O=OH+OH               2.97e+06   2.02  1.34e+4

! *************** H2-O2 Dissociation Reactions ******************

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2+M=H+H+M                4.577E+19 -1.40  1.0438E+05
   H2/2.5/ H2O/12/
   CO/1.9/ CO2/3.8/

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
O+O+M=O2+M                6.165E+15 -0.50  0.000E+00
   H2/2.5/ H2O/12/
   CO/1.9/ CO2/3.8/


! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
O+H+M=OH+M                4.714E+18 -1.00  0.000E+00
   H2/2.5/ H2O/12/
   CO/1.9/ CO2/3.8/

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
!H+OH+M=H2O+M              2.212E+22 -2.00  0.000E+00
H+OH+M=H2O+M               3.800E+22 -2.00  0.000E+00
   H2/2.5/ H2O/12/
   CO/1.9/ CO2/3.8/

!************** Formation and Consumption of HO2******************

! Cobos et al., J. Phys. Chem. 89:342 (1985) for kinf
! Michael, et al., J. Phys. Chem. A, 106:5297 (2002) for k0

!******************************************************************************
! MAIN BATH GAS IS N2 (comment this reaction otherwise)
!
 H+O2(+M)=HO2(+M)      1.475E+12  0.60  0.00E+00
     LOW/6.366E+20  -1.72  5.248E+02/
     TROE/0.8  1E-30  1E+30/
     H2/2.0/ H2O/11./ O2/0.78/ CO/1.9/ CO2/3.8/

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) [modified]
HO2+H=H2+O2               1.66E+13   0.00   0.823E+03

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) [modified]
HO2+H=OH+OH               7.079E+13   0.00   2.95E+02

! Baulch et al., J. Phys. Chem. Ref Data, 21:411 (1992)
HO2+O=O2+OH               0.325E+14  0.00   0.00E+00

! Keyser, J. Phys. Chem. 92:1193 (1988)
HO2+OH=H2O+O2             2.890E+13  0.00 -4.970E+02

! ***************Formation and Consumption of H2O2******************

! Hippler et al., J. Chem. Phys. 93:1755 (1990)
HO2+HO2=H2O2+O2            4.200e+14  0.00  1.1982e+04
  DUPLICATE
HO2+HO2=H2O2+O2            1.300e+11  0.00 -1.6293e+3
  DUPLICATE

! Brouwer et al., J. Chem. Phys. 86:6171 (1987) for kinf
! Warnatz, J. in Combustion chemistry (1984) for k0
H2O2(+M)=OH+OH(+M)         2.951e+14   0.00  4.843E+04
  LOW/1.202E+17  0.00  4.55E+04/
  TROE/0.5 1E-30 1E+30/
  H2/2.5/ H2O/12/
  CO/1.9/ CO2/3.8/


! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2O2+H=H2O+OH             0.241E+14  0.00  0.397E+04

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2O2+H=HO2+H2             0.482E+14  0.00  0.795E+04

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2O2+O=OH+HO2             9.550E+06  2.00  3.970E+03

! Hippler and Troe, J. Chem. Phys. Lett. 192:333 (1992)
H2O2+OH=HO2+H2O           1.000E+12  0.00  0.000
    DUPLICATE
H2O2+OH=HO2+H2O           5.800E+14  0.00  9.557E+03
    DUPLICATE

!************** CO/HCO REACTIONS *****************

! Troe, 15th Symposium
CO+O(+M)=CO2(+M)           1.80E+10  0.00  2384.
! Fit of Westmoreland, AiChe J., 1986, rel. to N2 - Tim adjusted from MTA's
! rate constant, which was rel to Ar.
   LOW/1.55E+24 -2.79  4191./
   H2/2.5/ H2O/12/ CO/1.9/ CO2/3.8/

! Tsang and Hampson, JPC Ref. Data, 15:1087 (1986)
CO+O2=CO2+O               0.253E+13  0.00  0.477E+05

! This rate constant is modified per an updated value for HO2+HO2=H2O2+OH
CO+HO2=CO2+OH             3.01E+13   0.00   2.30E+04

! This study (2004) by matching literature experiment results
CO+OH=CO2+H                    2.229E+05  1.89  -1158.7

! This study (2004) by matching literature experiment results
HCO+M=H+CO+M           4.7485E+11  0.659  1.4874E+04
H2/2.5/ H2O/6/ CO/1.9/ CO2/3.8/

! Timonen et al., JPC, 92:651 (1988)
HCO+O2=CO+HO2             0.758E+13  0.00  0.410E+03 

! Timonen et al., JPC, 91:692 (1987)
HCO+H=CO+H2               0.723E+14  0.00  0.000E+00  

! All reactions from Tsang and Hampson, JPC Ref. Data, 15:1087 (1986)
HCO + O = CO2 + H                 3.000E+13  0.00  0.000E+00

END

\\
\\
\\  This is the therm file
\\
\\
THERMO ALL
0300.00   1000.00 5000.00
HO2               L 5/89H   1O   2   00   00G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00 1.00021620E+04    4
CO                121286C   1O   1          G  0300.00   5000.00  1000.00      1
 0.03025078E+02 0.01442689E-01-0.05630828E-05 0.01018581E-08-0.06910952E-13    2
-0.01426835E+06 0.06108218E+02 0.03262452E+02 0.01511941E-01-0.03881755E-04    3
 0.05581944E-07-0.02474951E-10-0.01431054E+06 0.04848897E+02                   4
CO2               121286C   1O   2          G  0300.00   5000.00  1000.00      1
 0.04453623E+02 0.03140169E-01-0.01278411E-04 0.02393997E-08-0.01669033E-12    2
-0.04896696E+06-0.09553959E+01 0.02275725E+02 0.09922072E-01-0.01040911E-03    3
 0.06866687E-07-0.02117280E-10-0.04837314E+06 0.01018849E+03                   4
H                 120186H   1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 0.02547163E+06-0.04601176E+01 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.02547163E+06-0.04601176E+01                   4
H2                121286H   2               G  0300.00   5000.00  1000.00      1
 0.02991423E+02 0.07000644E-02-0.05633829E-06-0.09231578E-10 0.01582752E-13    2
-0.08350340E+04-0.01355110E+02 0.03298124E+02 0.08249442E-02-0.08143015E-05    3
-0.09475434E-09 0.04134872E-11-0.01012521E+05-0.03294094E+02                   4
H2O                20387H   2O   1          G  0300.00   5000.00  1000.00      1
 0.02672146E+02 0.03056293E-01-0.08730260E-05 0.01200996E-08-0.06391618E-13    2
-0.02989921E+06 0.06862817E+02 0.03386842E+02 0.03474982E-01-0.06354696E-04    3
 0.06968581E-07-0.02506588E-10-0.03020811E+06 0.02590233E+02                   4
H2O2              120186H   2O   2          G  0300.00   5000.00  1000.00      1
 0.04573167E+02 0.04336136E-01-0.01474689E-04 0.02348904E-08-0.01431654E-12    2
-0.01800696E+06 0.05011370E+01 0.03388754E+02 0.06569226E-01-0.01485013E-05    3
-0.04625806E-07 0.02471515E-10-0.01766315E+06 0.06785363E+02                   4
HCO               121286H   1C   1O   1     G  0300.00   5000.00  1000.00      1
 0.03557271E+02 0.03345573E-01-0.01335006E-04 0.02470573E-08-0.01713851E-12    2
 0.03916324E+05 0.05552299E+02 0.02898330E+02 0.06199147E-01-0.09623084E-04    3
 0.01089825E-06-0.04574885E-10 0.04159922E+05 0.08983614E+02                   4
O                 120186O   1               G  0300.00   5000.00  1000.00      1
 0.02542060E+02-0.02755062E-03-0.03102803E-07 0.04551067E-10-0.04368052E-14    2
 0.02923080E+06 0.04920308E+02 0.02946429E+02-0.01638166E-01 0.02421032E-04    3
-0.01602843E-07 0.03890696E-11 0.02914764E+06 0.02963995E+02                   4
O2                121386O   2               G  0300.00   5000.00  1000.00      1
 0.03697578E+02 0.06135197E-02-0.01258842E-05 0.01775281E-09-0.01136435E-13    2
-0.01233930E+05 0.03189166E+02 0.03212936E+02 0.01127486E-01-0.05756150E-05    3
 0.01313877E-07-0.08768554E-11-0.01005249E+05 0.06034738E+02                   4
N2                121286N   2               G  0300.00   5000.00  1000.00      1
 0.02926640E+02 0.01487977E-01-0.05684761E-05 0.01009704E-08-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.01408240E-01-0.03963222E-04    3
 0.05641515E-07-0.02444855E-10-0.01020900E+05 0.03950372E+02                   4
OH                S 9/01O   1H   1    0    0G   200.000  6000.000 1000.        1
 2.86472886E+00 1.05650448E-03-2.59082758E-07 3.05218674E-11-1.33195876E-15    2
 3.68362875E+03 5.70164073E+00 4.12530561E+00-3.22544939E-03 6.52764691E-06    3
-5.79853643E-09 2.06237379E-12 3.34630913E+03-6.90432960E-01 4.51532273E+03    4
END

#endif
