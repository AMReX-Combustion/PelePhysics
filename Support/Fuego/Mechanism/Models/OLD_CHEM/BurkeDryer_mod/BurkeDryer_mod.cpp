
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
extern "C" {
void egtransetEPS(double* EPS);
void egtransetSIG(double* SIG);
void atomicWeight(double* awt);
void molecularWeight(double* wt);
void gibbs(double* species, double* tc);
void helmholtz(double* species, double* tc);
void speciesInternalEnergy(double* species, double* tc);
void speciesEnthalpy(double* species, double* tc);
void speciesEntropy(double* species, double* tc);
void cp_R(double* species, double* tc);
void cv_R(double* species, double* tc);
void equilibriumConstants(double* kc, double* g_RT, double T);
void productionRate(double* wdot, double* sc, double T);
void comp_k_f(double* tc, double invT, double* k_f);
void comp_Kc(double* tc, double invT, double* Kc);
void comp_qfqr(double* q_f, double* q_r, double* sc, double* tc, double invT);
void progressRate(double* qdot, double* speciesConc, double T);
void progressRateFR(double* q_f, double* q_r, double* speciesConc, double T);
void CKINIT();
void CKFINALIZE();
void CKINDX(int* mm, int* kk, int* ii, int* nfit);
void CKXNUM(
  char* line,
  int* nexp,
  int* lout,
  int* nval,
  double* rval,
  int* kerr,
  int lenline);
void CKSNUM(
  char* line,
  int* nexp,
  int* lout,
  char* kray,
  int* nn,
  int* knum,
  int* nval,
  double* rval,
  int* kerr,
  int lenline,
  int lenkray);
void CKSYME(int* kname, int* lenkname);
void CKSYMS(int* kname, int* lenkname);
void CKRP(double* ru, double* ruc, double* pa);
void CKPX(double* rho, double* T, double* x, double* P);
void CKPY(double* rho, double* T, double* y, double* P);
void CKPC(double* rho, double* T, double* c, double* P);
void CKRHOX(double* P, double* T, double* x, double* rho);
void CKRHOY(double* P, double* T, double* y, double* rho);
void CKRHOC(double* P, double* T, double* c, double* rho);
void CKWT(double* wt);
void CKAWT(double* awt);
void CKMMWY(double* y, double* wtm);
void CKMMWX(double* x, double* wtm);
void CKMMWC(double* c, double* wtm);
void CKYTX(double* y, double* x);
void CKYTCP(double* P, double* T, double* y, double* c);
void CKYTCR(double* rho, double* T, double* y, double* c);
void CKXTY(double* x, double* y);
void CKXTCP(double* P, double* T, double* x, double* c);
void CKXTCR(double* rho, double* T, double* x, double* c);
void CKCTX(double* c, double* x);
void CKCTY(double* c, double* y);
void CKCPOR(double* T, double* cpor);
void CKHORT(double* T, double* hort);
void CKSOR(double* T, double* sor);
void CKCVML(double* T, double* cvml);
void CKCPML(double* T, double* cvml);
void CKUML(double* T, double* uml);
void CKHML(double* T, double* uml);
void CKGML(double* T, double* gml);
void CKAML(double* T, double* aml);
void CKSML(double* T, double* sml);
void CKCVMS(double* T, double* cvms);
void CKCPMS(double* T, double* cvms);
void CKUMS(double* T, double* ums);
void CKHMS(double* T, double* ums);
void CKGMS(double* T, double* gms);
void CKAMS(double* T, double* ams);
void CKSMS(double* T, double* sms);
void CKCPBL(double* T, double* x, double* cpbl);
void CKCPBS(double* T, double* y, double* cpbs);
void CKCVBL(double* T, double* x, double* cpbl);
void CKCVBS(double* T, double* y, double* cpbs);
void CKHBML(double* T, double* x, double* hbml);
void CKHBMS(double* T, double* y, double* hbms);
void CKUBML(double* T, double* x, double* ubml);
void CKUBMS(double* T, double* y, double* ubms);
void CKSBML(double* P, double* T, double* x, double* sbml);
void CKSBMS(double* P, double* T, double* y, double* sbms);
void CKGBML(double* P, double* T, double* x, double* gbml);
void CKGBMS(double* P, double* T, double* y, double* gbms);
void CKABML(double* P, double* T, double* x, double* abml);
void CKABMS(double* P, double* T, double* y, double* abms);
void CKWC(double* T, double* C, double* wdot);
void CKWYP(double* P, double* T, double* y, double* wdot);
void CKWXP(double* P, double* T, double* x, double* wdot);
void CKWYR(double* rho, double* T, double* y, double* wdot);
void CKWXR(double* rho, double* T, double* x, double* wdot);
void CKQC(double* T, double* C, double* qdot);
void CKKFKR(double* P, double* T, double* x, double* q_f, double* q_r);
void CKQYP(double* P, double* T, double* y, double* qdot);
void CKQXP(double* P, double* T, double* x, double* qdot);
void CKQYR(double* rho, double* T, double* y, double* qdot);
void CKQXR(double* rho, double* T, double* x, double* qdot);
void CKNU(int* kdim, int* nuki);
void CKNCF(int* mdim, int* ncf);
void CKABE(double* a, double* b, double* e);
void CKEQC(double* T, double* C, double* eqcon);
void CKEQYP(double* P, double* T, double* y, double* eqcon);
void CKEQXP(double* P, double* T, double* x, double* eqcon);
void CKEQYR(double* rho, double* T, double* y, double* eqcon);
void CKEQXR(double* rho, double* T, double* x, double* eqcon);
void DWDOT(double* J, double* sc, double* T, int* consP);
void DWDOT_PRECOND(double* J, double* sc, double* Tp, int* HP);
void SPARSITY_INFO(int* nJdata, int* consP, int NCELLS);
void SPARSITY_INFO_PRECOND(int* nJdata, int* consP);
void SPARSITY_PREPROC(int* rowVals, int* colPtrs, int* consP, int NCELLS);
void SPARSITY_PREPROC_PRECOND(int* rowVals, int* colPtrs, int* consP);
void aJacobian(double* J, double* sc, double T, int consP);
void aJacobian_precond(double* J, double* sc, double T, int HP);
void dcvpRdT(double* species, double* tc);
void GET_T_GIVEN_EY(double* e, double* y, double* t, int* ierr);
void GET_T_GIVEN_HY(double* h, double* y, double* t, int* ierr);
void GET_REACTION_MAP(int* rmap);
/*vector version */
void vproductionRate(int npt, double* wdot, double* c, double* T);
void VCKHMS(int* np, double* T, double* ums);
void VCKPY(int* np, double* rho, double* T, double* y, double* P);
void VCKWYR(int* np, double* rho, double* T, double* y, double* wdot);
void VCKYTX(int* np, double* y, double* x);
void vcomp_k_f(int npt, double* k_f_s, double* tc, double* invT);
void vcomp_gibbs(int npt, double* g_RT, double* tc);
void vcomp_Kc(int npt, double* Kc_s, double* g_RT, double* invT);
void GET_CRITPARAMS(double* Tci, double* ai, double* bi, double* acentric_i);
void vcomp_wdot(
  int npt,
  double* wdot,
  double* mixture,
  double* sc,
  double* k_f_s,
  double* Kc_s,
  double* tc,
  double* invT,
  double* T);

/* Inverse molecular weights */
static const double imw[11] = {
  1.0 / 1.007970,  /*H */
  1.0 / 2.015940,  /*H2 */
  1.0 / 15.999400, /*O */
  1.0 / 17.007370, /*OH */
  1.0 / 18.015340, /*H2O */
  1.0 / 31.998800, /*O2 */
  1.0 / 33.006770, /*HO2 */
  1.0 / 34.014740, /*H2O2 */
  1.0 / 28.013400, /*N2 */
  1.0 / 39.948000, /*AR */
  1.0 / 4.002600}; /*HE */

static double fwd_A[29], fwd_beta[29], fwd_Ea[29];
static double low_A[29], low_beta[29], low_Ea[29];
static double rev_A[29], rev_beta[29], rev_Ea[29];
static double troe_a[29], troe_Ts[29], troe_Tss[29], troe_Tsss[29];
static double sri_a[29], sri_b[29], sri_c[29], sri_d[29], sri_e[29];
static double activation_units[29], prefactor_units[29], phase_units[29];
static int is_PD[29], troe_len[29], sri_len[29], nTB[29], *TBid[29];
static double* TB[29];

static double fwd_A_DEF[29], fwd_beta_DEF[29], fwd_Ea_DEF[29];
static double low_A_DEF[29], low_beta_DEF[29], low_Ea_DEF[29];
static double rev_A_DEF[29], rev_beta_DEF[29], rev_Ea_DEF[29];
static double troe_a_DEF[29], troe_Ts_DEF[29], troe_Tss_DEF[29],
  troe_Tsss_DEF[29];
static double sri_a_DEF[29], sri_b_DEF[29], sri_c_DEF[29], sri_d_DEF[29],
  sri_e_DEF[29];
static double activation_units_DEF[29], prefactor_units_DEF[29],
  phase_units_DEF[29];
static int is_PD_DEF[29], troe_len_DEF[29], sri_len_DEF[29], nTB_DEF[29],
  *TBid_DEF[29];
static double* TB_DEF[29];
static int rxn_map[29] = {7,  8, 9,  10, 11, 2,  12, 13, 3,  14,
                          15, 4, 5,  16, 0,  17, 18, 19, 20, 21,
                          22, 1, 23, 24, 25, 26, 27, 28, 6};

void
GET_REACTION_MAP(int* rmap)
{
  for (int i = 0; i < 29; ++i) {
    rmap[i] = rxn_map[i];
  }
}

#include <ReactionData.H>
double*
GetParamPtr(
  int reaction_id, REACTION_PARAMETER param_id, int species_id, int get_default)
{
  double* ret = 0;
  if (reaction_id < 0 || reaction_id >= 29) {
    printf("Bad reaction id = %d", reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id < 0 || species_id >= 11) {
      printf("GetParamPtr: Bad species id = %d", species_id);
      abort();
    }
    if (get_default) {
      for (int i = 0; i < nTB_DEF[mrid]; ++i) {
        if (species_id == TBid_DEF[mrid][i]) {
          ret = &(TB_DEF[mrid][i]);
        }
      }
    } else {
      for (int i = 0; i < nTB[mrid]; ++i) {
        if (species_id == TBid[mrid][i]) {
          ret = &(TB[mrid][i]);
        }
      }
    }
    if (ret == 0) {
      printf("GetParamPtr: No TB for reaction id = %d", reaction_id);
      abort();
    }
  } else {
    if (param_id == FWD_A) {
      ret = (get_default ? &(fwd_A_DEF[mrid]) : &(fwd_A[mrid]));
    } else if (param_id == FWD_BETA) {
      ret = (get_default ? &(fwd_beta_DEF[mrid]) : &(fwd_beta[mrid]));
    } else if (param_id == FWD_EA) {
      ret = (get_default ? &(fwd_Ea_DEF[mrid]) : &(fwd_Ea[mrid]));
    } else if (param_id == LOW_A) {
      ret = (get_default ? &(low_A_DEF[mrid]) : &(low_A[mrid]));
    } else if (param_id == LOW_BETA) {
      ret = (get_default ? &(low_beta_DEF[mrid]) : &(low_beta[mrid]));
    } else if (param_id == LOW_EA) {
      ret = (get_default ? &(low_Ea_DEF[mrid]) : &(low_Ea[mrid]));
    } else if (param_id == REV_A) {
      ret = (get_default ? &(rev_A_DEF[mrid]) : &(rev_A[mrid]));
    } else if (param_id == REV_BETA) {
      ret = (get_default ? &(rev_beta_DEF[mrid]) : &(rev_beta[mrid]));
    } else if (param_id == REV_EA) {
      ret = (get_default ? &(rev_Ea_DEF[mrid]) : &(rev_Ea[mrid]));
    } else if (param_id == TROE_A) {
      ret = (get_default ? &(troe_a_DEF[mrid]) : &(troe_a[mrid]));
    } else if (param_id == TROE_TS) {
      ret = (get_default ? &(troe_Ts_DEF[mrid]) : &(troe_Ts[mrid]));
    } else if (param_id == TROE_TSS) {
      ret = (get_default ? &(troe_Tss_DEF[mrid]) : &(troe_Tss[mrid]));
    } else if (param_id == TROE_TSSS) {
      ret = (get_default ? &(troe_Tsss_DEF[mrid]) : &(troe_Tsss[mrid]));
    } else if (param_id == SRI_A) {
      ret = (get_default ? &(sri_a_DEF[mrid]) : &(sri_a[mrid]));
    } else if (param_id == SRI_B) {
      ret = (get_default ? &(sri_b_DEF[mrid]) : &(sri_b[mrid]));
    } else if (param_id == SRI_C) {
      ret = (get_default ? &(sri_c_DEF[mrid]) : &(sri_c[mrid]));
    } else if (param_id == SRI_D) {
      ret = (get_default ? &(sri_d_DEF[mrid]) : &(sri_d[mrid]));
    } else if (param_id == SRI_E) {
      ret = (get_default ? &(sri_e_DEF[mrid]) : &(sri_e[mrid]));
    } else {
      printf("GetParamPtr: Unknown parameter id");
      abort();
    }
  }
  return ret;
}

void
ResetAllParametersToDefault()
{
  for (int i = 0; i < 29; i++) {
    if (nTB[i] != 0) {
      nTB[i] = 0;
      free(TB[i]);
      free(TBid[i]);
    }

    fwd_A[i] = fwd_A_DEF[i];
    fwd_beta[i] = fwd_beta_DEF[i];
    fwd_Ea[i] = fwd_Ea_DEF[i];

    low_A[i] = low_A_DEF[i];
    low_beta[i] = low_beta_DEF[i];
    low_Ea[i] = low_Ea_DEF[i];

    rev_A[i] = rev_A_DEF[i];
    rev_beta[i] = rev_beta_DEF[i];
    rev_Ea[i] = rev_Ea_DEF[i];

    troe_a[i] = troe_a_DEF[i];
    troe_Ts[i] = troe_Ts_DEF[i];
    troe_Tss[i] = troe_Tss_DEF[i];
    troe_Tsss[i] = troe_Tsss_DEF[i];

    sri_a[i] = sri_a_DEF[i];
    sri_b[i] = sri_b_DEF[i];
    sri_c[i] = sri_c_DEF[i];
    sri_d[i] = sri_d_DEF[i];
    sri_e[i] = sri_e_DEF[i];

    is_PD[i] = is_PD_DEF[i];
    troe_len[i] = troe_len_DEF[i];
    sri_len[i] = sri_len_DEF[i];

    activation_units[i] = activation_units_DEF[i];
    prefactor_units[i] = prefactor_units_DEF[i];
    phase_units[i] = phase_units_DEF[i];

    nTB[i] = nTB_DEF[i];
    if (nTB[i] != 0) {
      TB[i] = (double*)malloc(sizeof(double) * nTB[i]);
      TBid[i] = (int*)malloc(sizeof(int) * nTB[i]);
      for (int j = 0; j < nTB[i]; j++) {
        TB[i][j] = TB_DEF[i][j];
        TBid[i][j] = TBid_DEF[i][j];
      }
    }
  }
}

void
SetAllDefaults()
{
  for (int i = 0; i < 29; i++) {
    if (nTB_DEF[i] != 0) {
      nTB_DEF[i] = 0;
      free(TB_DEF[i]);
      free(TBid_DEF[i]);
    }

    fwd_A_DEF[i] = fwd_A[i];
    fwd_beta_DEF[i] = fwd_beta[i];
    fwd_Ea_DEF[i] = fwd_Ea[i];

    low_A_DEF[i] = low_A[i];
    low_beta_DEF[i] = low_beta[i];
    low_Ea_DEF[i] = low_Ea[i];

    rev_A_DEF[i] = rev_A[i];
    rev_beta_DEF[i] = rev_beta[i];
    rev_Ea_DEF[i] = rev_Ea[i];

    troe_a_DEF[i] = troe_a[i];
    troe_Ts_DEF[i] = troe_Ts[i];
    troe_Tss_DEF[i] = troe_Tss[i];
    troe_Tsss_DEF[i] = troe_Tsss[i];

    sri_a_DEF[i] = sri_a[i];
    sri_b_DEF[i] = sri_b[i];
    sri_c_DEF[i] = sri_c[i];
    sri_d_DEF[i] = sri_d[i];
    sri_e_DEF[i] = sri_e[i];

    is_PD_DEF[i] = is_PD[i];
    troe_len_DEF[i] = troe_len[i];
    sri_len_DEF[i] = sri_len[i];

    activation_units_DEF[i] = activation_units[i];
    prefactor_units_DEF[i] = prefactor_units[i];
    phase_units_DEF[i] = phase_units[i];

    nTB_DEF[i] = nTB[i];
    if (nTB_DEF[i] != 0) {
      TB_DEF[i] = (double*)malloc(sizeof(double) * nTB_DEF[i]);
      TBid_DEF[i] = (int*)malloc(sizeof(int) * nTB_DEF[i]);
      for (int j = 0; j < nTB_DEF[i]; j++) {
        TB_DEF[i][j] = TB[i][j];
        TBid_DEF[i][j] = TBid[i][j];
      }
    }
  }
}

/* Finalizes parameter database */
void
CKFINALIZE()
{
  for (int i = 0; i < 29; ++i) {
    free(TB[i]);
    TB[i] = 0;
    free(TBid[i]);
    TBid[i] = 0;
    nTB[i] = 0;

    free(TB_DEF[i]);
    TB_DEF[i] = 0;
    free(TBid_DEF[i]);
    TBid_DEF[i] = 0;
    nTB_DEF[i] = 0;
  }
}

/* Initializes parameter database */
void
CKINIT()
{
  // (0):  H + O2 <=> O + OH
  fwd_A[7] = 104000000000000;
  fwd_beta[7] = 0;
  fwd_Ea[7] = 15286;
  prefactor_units[7] = 1.0000000000000002e-06;
  activation_units[7] = 0.50321666580471969;
  phase_units[7] = 1e-12;
  is_PD[7] = 0;
  nTB[7] = 0;

  // (1):  O + H2 <=> H + OH
  fwd_A[8] = 3818000000000;
  fwd_beta[8] = 0;
  fwd_Ea[8] = 7948;
  prefactor_units[8] = 1.0000000000000002e-06;
  activation_units[8] = 0.50321666580471969;
  phase_units[8] = 1e-12;
  is_PD[8] = 0;
  nTB[8] = 0;

  // (2):  O + H2 <=> H + OH
  fwd_A[9] = 879200000000000;
  fwd_beta[9] = 0;
  fwd_Ea[9] = 19170;
  prefactor_units[9] = 1.0000000000000002e-06;
  activation_units[9] = 0.50321666580471969;
  phase_units[9] = 1e-12;
  is_PD[9] = 0;
  nTB[9] = 0;

  // (3):  H2 + OH <=> H2O + H
  fwd_A[10] = 216000000;
  fwd_beta[10] = 1.51;
  fwd_Ea[10] = 3430;
  prefactor_units[10] = 1.0000000000000002e-06;
  activation_units[10] = 0.50321666580471969;
  phase_units[10] = 1e-12;
  is_PD[10] = 0;
  nTB[10] = 0;

  // (4):  OH + OH <=> O + H2O
  fwd_A[11] = 33400;
  fwd_beta[11] = 2.4199999999999999;
  fwd_Ea[11] = -1930;
  prefactor_units[11] = 1.0000000000000002e-06;
  activation_units[11] = 0.50321666580471969;
  phase_units[11] = 1e-12;
  is_PD[11] = 0;
  nTB[11] = 0;

  // (5):  H2 + M <=> H + H + M
  fwd_A[2] = 4.577e+19;
  fwd_beta[2] = -1.3999999999999999;
  fwd_Ea[2] = 104380;
  prefactor_units[2] = 1.0000000000000002e-06;
  activation_units[2] = 0.50321666580471969;
  phase_units[2] = 1e-6;
  is_PD[2] = 0;
  nTB[2] = 4;
  TB[2] = (double*)malloc(4 * sizeof(double));
  TBid[2] = (int*)malloc(4 * sizeof(int));
  TBid[2][0] = 1;
  TB[2][0] = 2.5; // H2
  TBid[2][1] = 4;
  TB[2][1] = 12; // H2O
  TBid[2][2] = 9;
  TB[2][2] = 0; // AR
  TBid[2][3] = 10;
  TB[2][3] = 0; // HE

  // (6):  H2 + AR <=> H + H + AR
  fwd_A[12] = 5.84e+18;
  fwd_beta[12] = -1.1000000000000001;
  fwd_Ea[12] = 104380;
  prefactor_units[12] = 1.0000000000000002e-06;
  activation_units[12] = 0.50321666580471969;
  phase_units[12] = 1e-12;
  is_PD[12] = 0;
  nTB[12] = 0;

  // (7):  H2 + HE <=> H + H + HE
  fwd_A[13] = 5.84e+18;
  fwd_beta[13] = -1.1000000000000001;
  fwd_Ea[13] = 104380;
  prefactor_units[13] = 1.0000000000000002e-06;
  activation_units[13] = 0.50321666580471969;
  phase_units[13] = 1e-12;
  is_PD[13] = 0;
  nTB[13] = 0;

  // (8):  O + O + M <=> O2 + M
  fwd_A[3] = 6165000000000000;
  fwd_beta[3] = -0.5;
  fwd_Ea[3] = 0;
  prefactor_units[3] = 1.0000000000000002e-12;
  activation_units[3] = 0.50321666580471969;
  phase_units[3] = 1e-12;
  is_PD[3] = 0;
  nTB[3] = 4;
  TB[3] = (double*)malloc(4 * sizeof(double));
  TBid[3] = (int*)malloc(4 * sizeof(int));
  TBid[3][0] = 1;
  TB[3][0] = 2.5; // H2
  TBid[3][1] = 4;
  TB[3][1] = 12; // H2O
  TBid[3][2] = 9;
  TB[3][2] = 0; // AR
  TBid[3][3] = 10;
  TB[3][3] = 0; // HE

  // (9):  O + O + AR <=> O2 + AR
  fwd_A[14] = 18860000000000;
  fwd_beta[14] = 0;
  fwd_Ea[14] = -1788;
  prefactor_units[14] = 1.0000000000000002e-12;
  activation_units[14] = 0.50321666580471969;
  phase_units[14] = 1e-18;
  is_PD[14] = 0;
  nTB[14] = 0;

  // (10):  O + O + HE <=> O2 + HE
  fwd_A[15] = 18860000000000;
  fwd_beta[15] = 0;
  fwd_Ea[15] = -1788;
  prefactor_units[15] = 1.0000000000000002e-12;
  activation_units[15] = 0.50321666580471969;
  phase_units[15] = 1e-18;
  is_PD[15] = 0;
  nTB[15] = 0;

  // (11):  O + H + M <=> OH + M
  fwd_A[4] = 4.714e+18;
  fwd_beta[4] = -1;
  fwd_Ea[4] = 0;
  prefactor_units[4] = 1.0000000000000002e-12;
  activation_units[4] = 0.50321666580471969;
  phase_units[4] = 1e-12;
  is_PD[4] = 0;
  nTB[4] = 4;
  TB[4] = (double*)malloc(4 * sizeof(double));
  TBid[4] = (int*)malloc(4 * sizeof(int));
  TBid[4][0] = 1;
  TB[4][0] = 2.5; // H2
  TBid[4][1] = 4;
  TB[4][1] = 12; // H2O
  TBid[4][2] = 9;
  TB[4][2] = 0.75; // AR
  TBid[4][3] = 10;
  TB[4][3] = 0.75; // HE

  // (12):  H2O + M <=> H + OH + M
  fwd_A[5] = 6.0640000000000002e+27;
  fwd_beta[5] = -3.3220000000000001;
  fwd_Ea[5] = 120790;
  prefactor_units[5] = 1.0000000000000002e-06;
  activation_units[5] = 0.50321666580471969;
  phase_units[5] = 1e-6;
  is_PD[5] = 0;
  nTB[5] = 5;
  TB[5] = (double*)malloc(5 * sizeof(double));
  TBid[5] = (int*)malloc(5 * sizeof(int));
  TBid[5][0] = 1;
  TB[5][0] = 3; // H2
  TBid[5][1] = 4;
  TB[5][1] = 0; // H2O
  TBid[5][2] = 10;
  TB[5][2] = 1.1000000000000001; // HE
  TBid[5][3] = 8;
  TB[5][3] = 2; // N2
  TBid[5][4] = 5;
  TB[5][4] = 1.5; // O2

  // (13):  H2O + H2O <=> H + OH + H2O
  fwd_A[16] = 1.006e+26;
  fwd_beta[16] = -2.4399999999999999;
  fwd_Ea[16] = 120180;
  prefactor_units[16] = 1.0000000000000002e-06;
  activation_units[16] = 0.50321666580471969;
  phase_units[16] = 1e-12;
  is_PD[16] = 0;
  nTB[16] = 0;

  // (14):  H + O2 (+M) <=> HO2 (+M)
  fwd_A[0] = 4650840000000;
  fwd_beta[0] = 0.44;
  fwd_Ea[0] = 0;
  low_A[0] = 6.366e+20;
  low_beta[0] = -1.72;
  low_Ea[0] = 524.79999999999995;
  troe_a[0] = 0.5;
  troe_Tsss[0] = 1.0000000000000001e-30;
  troe_Ts[0] = 1e+30;
  troe_len[0] = 3;
  prefactor_units[0] = 1.0000000000000002e-06;
  activation_units[0] = 0.50321666580471969;
  phase_units[0] = 1e-12;
  is_PD[0] = 1;
  nTB[0] = 5;
  TB[0] = (double*)malloc(5 * sizeof(double));
  TBid[0] = (int*)malloc(5 * sizeof(int));
  TBid[0][0] = 1;
  TB[0][0] = 2; // H2
  TBid[0][1] = 4;
  TB[0][1] = 14; // H2O
  TBid[0][2] = 5;
  TB[0][2] = 0.78000000000000003; // O2
  TBid[0][3] = 9;
  TB[0][3] = 0.67000000000000004; // AR
  TBid[0][4] = 10;
  TB[0][4] = 0.80000000000000004; // HE

  // (15):  HO2 + H <=> H2 + O2
  fwd_A[17] = 2750000;
  fwd_beta[17] = 2.0899999999999999;
  fwd_Ea[17] = -1451;
  prefactor_units[17] = 1.0000000000000002e-06;
  activation_units[17] = 0.50321666580471969;
  phase_units[17] = 1e-12;
  is_PD[17] = 0;
  nTB[17] = 0;

  // (16):  HO2 + H <=> OH + OH
  fwd_A[18] = 70790000000000;
  fwd_beta[18] = 0;
  fwd_Ea[18] = 295;
  prefactor_units[18] = 1.0000000000000002e-06;
  activation_units[18] = 0.50321666580471969;
  phase_units[18] = 1e-12;
  is_PD[18] = 0;
  nTB[18] = 0;

  // (17):  HO2 + O <=> O2 + OH
  fwd_A[19] = 28500000000;
  fwd_beta[19] = 1;
  fwd_Ea[19] = -723.92999999999995;
  prefactor_units[19] = 1.0000000000000002e-06;
  activation_units[19] = 0.50321666580471969;
  phase_units[19] = 1e-12;
  is_PD[19] = 0;
  nTB[19] = 0;

  // (18):  HO2 + OH <=> H2O + O2
  fwd_A[20] = 28900000000000;
  fwd_beta[20] = 0;
  fwd_Ea[20] = -497;
  prefactor_units[20] = 1.0000000000000002e-06;
  activation_units[20] = 0.50321666580471969;
  phase_units[20] = 1e-12;
  is_PD[20] = 0;
  nTB[20] = 0;

  // (19):  HO2 + HO2 <=> H2O2 + O2
  fwd_A[21] = 420000000000000;
  fwd_beta[21] = 0;
  fwd_Ea[21] = 11982;
  prefactor_units[21] = 1.0000000000000002e-06;
  activation_units[21] = 0.50321666580471969;
  phase_units[21] = 1e-12;
  is_PD[21] = 0;
  nTB[21] = 0;

  // (20):  HO2 + HO2 <=> H2O2 + O2
  fwd_A[22] = 130000000000;
  fwd_beta[22] = 0;
  fwd_Ea[22] = -1629.3;
  prefactor_units[22] = 1.0000000000000002e-06;
  activation_units[22] = 0.50321666580471969;
  phase_units[22] = 1e-12;
  is_PD[22] = 0;
  nTB[22] = 0;

  // (21):  H2O2 (+M) <=> OH + OH (+M)
  fwd_A[1] = 2000000000000;
  fwd_beta[1] = 0.90000000000000002;
  fwd_Ea[1] = 48749;
  low_A[1] = 2.49e+24;
  low_beta[1] = -2.2999999999999998;
  low_Ea[1] = 48749;
  troe_a[1] = 0.42999999999999999;
  troe_Tsss[1] = 1.0000000000000001e-30;
  troe_Ts[1] = 1e+30;
  troe_len[1] = 3;
  prefactor_units[1] = 1;
  activation_units[1] = 0.50321666580471969;
  phase_units[1] = 1e-6;
  is_PD[1] = 1;
  nTB[1] = 6;
  TB[1] = (double*)malloc(6 * sizeof(double));
  TBid[1] = (int*)malloc(6 * sizeof(int));
  TBid[1][0] = 4;
  TB[1][0] = 7.5; // H2O
  TBid[1][1] = 8;
  TB[1][1] = 1.5; // N2
  TBid[1][2] = 5;
  TB[1][2] = 1.2; // O2
  TBid[1][3] = 10;
  TB[1][3] = 0.65000000000000002; // HE
  TBid[1][4] = 7;
  TB[1][4] = 7.7000000000000002; // H2O2
  TBid[1][5] = 1;
  TB[1][5] = 3.7000000000000002; // H2

  // (22):  H2O2 + H <=> H2O + OH
  fwd_A[23] = 24100000000000;
  fwd_beta[23] = 0;
  fwd_Ea[23] = 3970;
  prefactor_units[23] = 1.0000000000000002e-06;
  activation_units[23] = 0.50321666580471969;
  phase_units[23] = 1e-12;
  is_PD[23] = 0;
  nTB[23] = 0;

  // (23):  H2O2 + H <=> HO2 + H2
  fwd_A[24] = 48200000000000;
  fwd_beta[24] = 0;
  fwd_Ea[24] = 7950;
  prefactor_units[24] = 1.0000000000000002e-06;
  activation_units[24] = 0.50321666580471969;
  phase_units[24] = 1e-12;
  is_PD[24] = 0;
  nTB[24] = 0;

  // (24):  H2O2 + O <=> OH + HO2
  fwd_A[25] = 9550000;
  fwd_beta[25] = 2;
  fwd_Ea[25] = 3970;
  prefactor_units[25] = 1.0000000000000002e-06;
  activation_units[25] = 0.50321666580471969;
  phase_units[25] = 1e-12;
  is_PD[25] = 0;
  nTB[25] = 0;

  // (25):  H2O2 + OH <=> HO2 + H2O
  fwd_A[26] = 1740000000000;
  fwd_beta[26] = 0;
  fwd_Ea[26] = 318;
  prefactor_units[26] = 1.0000000000000002e-06;
  activation_units[26] = 0.50321666580471969;
  phase_units[26] = 1e-12;
  is_PD[26] = 0;
  nTB[26] = 0;

  // (26):  H2O2 + OH <=> HO2 + H2O
  fwd_A[27] = 75900000000000;
  fwd_beta[27] = 0;
  fwd_Ea[27] = 7270;
  prefactor_units[27] = 1.0000000000000002e-06;
  activation_units[27] = 0.50321666580471969;
  phase_units[27] = 1e-12;
  is_PD[27] = 0;
  nTB[27] = 0;

  // (27):  HO2 + H <=> O + H2O
  fwd_A[28] = 3970000000000;
  fwd_beta[28] = 0;
  fwd_Ea[28] = 671;
  prefactor_units[28] = 1.0000000000000002e-06;
  activation_units[28] = 0.50321666580471969;
  phase_units[28] = 1e-12;
  is_PD[28] = 0;
  nTB[28] = 0;

  // (28):  O + OH + M <=> HO2 + M
  fwd_A[6] = 8000000000000000;
  fwd_beta[6] = 0;
  fwd_Ea[6] = 0;
  prefactor_units[6] = 1.0000000000000002e-12;
  activation_units[6] = 0.50321666580471969;
  phase_units[6] = 1e-12;
  is_PD[6] = 0;
  nTB[6] = 4;
  TB[6] = (double*)malloc(4 * sizeof(double));
  TBid[6] = (int*)malloc(4 * sizeof(int));
  TBid[6][0] = 1;
  TB[6][0] = 2; // H2
  TBid[6][1] = 4;
  TB[6][1] = 12; // H2O
  TBid[6][2] = 9;
  TB[6][2] = 0.69999999999999996; // AR
  TBid[6][3] = 10;
  TB[6][3] = 0.69999999999999996; // HE

  SetAllDefaults();
}

/*A few mechanism parameters */
void
CKINDX(int* mm, int* kk, int* ii, int* nfit)
{
  *mm = 6;
  *kk = 11;
  *ii = 29;
  *nfit = -1; /*Why do you need this anyway ?  */
}

/* ckxnum... for parsing strings  */
void
CKXNUM(
  char* line,
  int* nexp,
  int* lout,
  int* nval,
  double* rval,
  int* kerr,
  int lenline)
{
  int n, i; /*Loop Counters */
  char cstr[1000];
  char* saveptr;
  char* p; /*String Tokens */
  /* Strip Comments  */
  for (i = 0; i < lenline; ++i) {
    if (line[i] == '!') {
      break;
    }
    cstr[i] = line[i];
  }
  cstr[i] = '\0';

  p = strtok_r(cstr, " ", &saveptr);
  if (!p) {
    *nval = 0;
    *kerr = 1;
    return;
  }
  for (n = 0; n < *nexp; ++n) {
    rval[n] = atof(p);
    p = strtok_r(NULL, " ", &saveptr);
    if (!p)
      break;
  }
  *nval = n + 1;
  if (*nval < *nexp)
    *kerr = 1;
  return;
}

/* cksnum... for parsing strings  */
void
CKSNUM(
  char* line,
  int* nexp,
  int* lout,
  char* kray,
  int* nn,
  int* knum,
  int* nval,
  double* rval,
  int* kerr,
  int lenline,
  int lenkray)
{
  /*Not done yet ... */
}

/* Returns the char strings of element names */
void
CKSYME(int* kname, int* plenkname)
{
  int i; /*Loop Counter */
  int lenkname = *plenkname;
  /*clear kname */
  for (i = 0; i < lenkname * 6; i++) {
    kname[i] = ' ';
  }

  /* H  */
  kname[0 * lenkname + 0] = 'H';
  kname[0 * lenkname + 1] = ' ';

  /* O  */
  kname[1 * lenkname + 0] = 'O';
  kname[1 * lenkname + 1] = ' ';

  /* N  */
  kname[2 * lenkname + 0] = 'N';
  kname[2 * lenkname + 1] = ' ';

  /* AR  */
  kname[3 * lenkname + 0] = 'A';
  kname[3 * lenkname + 1] = 'R';
  kname[3 * lenkname + 2] = ' ';

  /* HE  */
  kname[4 * lenkname + 0] = 'H';
  kname[4 * lenkname + 1] = 'E';
  kname[4 * lenkname + 2] = ' ';

  /* C  */
  kname[5 * lenkname + 0] = 'C';
  kname[5 * lenkname + 1] = ' ';
}

/* Returns the char strings of species names */
void
CKSYMS(int* kname, int* plenkname)
{
  int i; /*Loop Counter */
  int lenkname = *plenkname;
  /*clear kname */
  for (i = 0; i < lenkname * 11; i++) {
    kname[i] = ' ';
  }

  /* H  */
  kname[0 * lenkname + 0] = 'H';
  kname[0 * lenkname + 1] = ' ';

  /* H2  */
  kname[1 * lenkname + 0] = 'H';
  kname[1 * lenkname + 1] = '2';
  kname[1 * lenkname + 2] = ' ';

  /* O  */
  kname[2 * lenkname + 0] = 'O';
  kname[2 * lenkname + 1] = ' ';

  /* OH  */
  kname[3 * lenkname + 0] = 'O';
  kname[3 * lenkname + 1] = 'H';
  kname[3 * lenkname + 2] = ' ';

  /* H2O  */
  kname[4 * lenkname + 0] = 'H';
  kname[4 * lenkname + 1] = '2';
  kname[4 * lenkname + 2] = 'O';
  kname[4 * lenkname + 3] = ' ';

  /* O2  */
  kname[5 * lenkname + 0] = 'O';
  kname[5 * lenkname + 1] = '2';
  kname[5 * lenkname + 2] = ' ';

  /* HO2  */
  kname[6 * lenkname + 0] = 'H';
  kname[6 * lenkname + 1] = 'O';
  kname[6 * lenkname + 2] = '2';
  kname[6 * lenkname + 3] = ' ';

  /* H2O2  */
  kname[7 * lenkname + 0] = 'H';
  kname[7 * lenkname + 1] = '2';
  kname[7 * lenkname + 2] = 'O';
  kname[7 * lenkname + 3] = '2';
  kname[7 * lenkname + 4] = ' ';

  /* N2  */
  kname[8 * lenkname + 0] = 'N';
  kname[8 * lenkname + 1] = '2';
  kname[8 * lenkname + 2] = ' ';

  /* AR  */
  kname[9 * lenkname + 0] = 'A';
  kname[9 * lenkname + 1] = 'R';
  kname[9 * lenkname + 2] = ' ';

  /* HE  */
  kname[10 * lenkname + 0] = 'H';
  kname[10 * lenkname + 1] = 'E';
  kname[10 * lenkname + 2] = ' ';
}

/* Returns R, Rc, Patm */
void
CKRP(double* ru, double* ruc, double* pa)
{
  *ru = 8.31451e+07;
  *ruc = 1.98721558317399615845;
  *pa = 1.01325e+06;
}

/*Compute P = rhoRT/W(x) */
void
CKPX(double* rho, double* T, double* x, double* P)
{
  double XW = 0;                       /* To hold mean molecular wt */
  XW += x[0] * 1.007970;               /*H */
  XW += x[1] * 2.015940;               /*H2 */
  XW += x[2] * 15.999400;              /*O */
  XW += x[3] * 17.007370;              /*OH */
  XW += x[4] * 18.015340;              /*H2O */
  XW += x[5] * 31.998800;              /*O2 */
  XW += x[6] * 33.006770;              /*HO2 */
  XW += x[7] * 34.014740;              /*H2O2 */
  XW += x[8] * 28.013400;              /*N2 */
  XW += x[9] * 39.948000;              /*AR */
  XW += x[10] * 4.002600;              /*HE */
  *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

  return;
}

/*Compute P = rhoRT/W(y) */
void
CKPY(double* rho, double* T, double* y, double* P)
{
  double YOW = 0;                       /* for computing mean MW */
  YOW += y[0] * imw[0];                 /*H */
  YOW += y[1] * imw[1];                 /*H2 */
  YOW += y[2] * imw[2];                 /*O */
  YOW += y[3] * imw[3];                 /*OH */
  YOW += y[4] * imw[4];                 /*H2O */
  YOW += y[5] * imw[5];                 /*O2 */
  YOW += y[6] * imw[6];                 /*HO2 */
  YOW += y[7] * imw[7];                 /*H2O2 */
  YOW += y[8] * imw[8];                 /*N2 */
  YOW += y[9] * imw[9];                 /*AR */
  YOW += y[10] * imw[10];               /*HE */
  *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

  return;
}

/*Compute P = rhoRT/W(y) */
void
VCKPY(int* np, double* rho, double* T, double* y, double* P)
{
  double YOW[*np];
  for (int i = 0; i < (*np); i++) {
    YOW[i] = 0.0;
  }

  for (int n = 0; n < 11; n++) {
    for (int i = 0; i < (*np); i++) {
      YOW[i] += y[n * (*np) + i] * imw[n];
    }
  }

  for (int i = 0; i < (*np); i++) {
    P[i] = rho[i] * 8.31451e+07 * T[i] * YOW[i]; /*P = rho*R*T/W */
  }

  return;
}

/*Compute P = rhoRT/W(c) */
void
CKPC(double* rho, double* T, double* c, double* P)
{
  int id; /*loop counter */
  /*See Eq 5 in CK Manual */
  double W = 0;
  double sumC = 0;
  W += c[0] * 1.007970;  /*H */
  W += c[1] * 2.015940;  /*H2 */
  W += c[2] * 15.999400; /*O */
  W += c[3] * 17.007370; /*OH */
  W += c[4] * 18.015340; /*H2O */
  W += c[5] * 31.998800; /*O2 */
  W += c[6] * 33.006770; /*HO2 */
  W += c[7] * 34.014740; /*H2O2 */
  W += c[8] * 28.013400; /*N2 */
  W += c[9] * 39.948000; /*AR */
  W += c[10] * 4.002600; /*HE */

  for (id = 0; id < 11; ++id) {
    sumC += c[id];
  }
  *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

  return;
}

/*Compute rho = PW(x)/RT */
void
CKRHOX(double* P, double* T, double* x, double* rho)
{
  double XW = 0;                         /* To hold mean molecular wt */
  XW += x[0] * 1.007970;                 /*H */
  XW += x[1] * 2.015940;                 /*H2 */
  XW += x[2] * 15.999400;                /*O */
  XW += x[3] * 17.007370;                /*OH */
  XW += x[4] * 18.015340;                /*H2O */
  XW += x[5] * 31.998800;                /*O2 */
  XW += x[6] * 33.006770;                /*HO2 */
  XW += x[7] * 34.014740;                /*H2O2 */
  XW += x[8] * 28.013400;                /*N2 */
  XW += x[9] * 39.948000;                /*AR */
  XW += x[10] * 4.002600;                /*HE */
  *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

  return;
}

/*Compute rho = P*W(y)/RT */
void
CKRHOY(double* P, double* T, double* y, double* rho)
{
  double YOW = 0;
  double tmp[11];

  for (int i = 0; i < 11; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 11; i++) {
    YOW += tmp[i];
  }

  *rho = *P / (8.31451e+07 * (*T) * YOW); /*rho = P*W/(R*T) */
  return;
}

/*Compute rho = P*W(c)/(R*T) */
void
CKRHOC(double* P, double* T, double* c, double* rho)
{
  int id; /*loop counter */
  /*See Eq 5 in CK Manual */
  double W = 0;
  double sumC = 0;
  W += c[0] * 1.007970;  /*H */
  W += c[1] * 2.015940;  /*H2 */
  W += c[2] * 15.999400; /*O */
  W += c[3] * 17.007370; /*OH */
  W += c[4] * 18.015340; /*H2O */
  W += c[5] * 31.998800; /*O2 */
  W += c[6] * 33.006770; /*HO2 */
  W += c[7] * 34.014740; /*H2O2 */
  W += c[8] * 28.013400; /*N2 */
  W += c[9] * 39.948000; /*AR */
  W += c[10] * 4.002600; /*HE */

  for (id = 0; id < 11; ++id) {
    sumC += c[id];
  }
  *rho = *P * W / (sumC * (*T) * 8.31451e+07); /*rho = PW/(R*T) */

  return;
}

/*get molecular weight for all species */
void
CKWT(double* wt)
{
  molecularWeight(wt);
}

/*get atomic weight for all elements */
void
CKAWT(double* awt)
{
  atomicWeight(awt);
}

/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void
CKMMWY(double* y, double* wtm)
{
  double YOW = 0;
  double tmp[11];

  for (int i = 0; i < 11; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 11; i++) {
    YOW += tmp[i];
  }

  *wtm = 1.0 / YOW;
  return;
}

/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void
CKMMWX(double* x, double* wtm)
{
  double XW = 0;          /* see Eq 4 in CK Manual */
  XW += x[0] * 1.007970;  /*H */
  XW += x[1] * 2.015940;  /*H2 */
  XW += x[2] * 15.999400; /*O */
  XW += x[3] * 17.007370; /*OH */
  XW += x[4] * 18.015340; /*H2O */
  XW += x[5] * 31.998800; /*O2 */
  XW += x[6] * 33.006770; /*HO2 */
  XW += x[7] * 34.014740; /*H2O2 */
  XW += x[8] * 28.013400; /*N2 */
  XW += x[9] * 39.948000; /*AR */
  XW += x[10] * 4.002600; /*HE */
  *wtm = XW;

  return;
}

/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void
CKMMWC(double* c, double* wtm)
{
  int id; /*loop counter */
  /*See Eq 5 in CK Manual */
  double W = 0;
  double sumC = 0;
  W += c[0] * 1.007970;  /*H */
  W += c[1] * 2.015940;  /*H2 */
  W += c[2] * 15.999400; /*O */
  W += c[3] * 17.007370; /*OH */
  W += c[4] * 18.015340; /*H2O */
  W += c[5] * 31.998800; /*O2 */
  W += c[6] * 33.006770; /*HO2 */
  W += c[7] * 34.014740; /*H2O2 */
  W += c[8] * 28.013400; /*N2 */
  W += c[9] * 39.948000; /*AR */
  W += c[10] * 4.002600; /*HE */

  for (id = 0; id < 11; ++id) {
    sumC += c[id];
  }
  /* CK provides no guard against divison by zero */
  *wtm = W / sumC;

  return;
}

/*convert y[species] (mass fracs) to x[species] (mole fracs) */
void
CKYTX(double* y, double* x)
{
  double YOW = 0;
  double tmp[11];

  for (int i = 0; i < 11; i++) {
    tmp[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 11; i++) {
    YOW += tmp[i];
  }

  double YOWINV = 1.0 / YOW;

  for (int i = 0; i < 11; i++) {
    x[i] = y[i] * imw[i] * YOWINV;
  }
  return;
}

/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void
VCKYTX(int* np, double* y, double* x)
{
  double YOW[*np];
  for (int i = 0; i < (*np); i++) {
    YOW[i] = 0.0;
  }

  for (int n = 0; n < 11; n++) {
    for (int i = 0; i < (*np); i++) {
      x[n * (*np) + i] = y[n * (*np) + i] * imw[n];
      YOW[i] += x[n * (*np) + i];
    }
  }

  for (int i = 0; i < (*np); i++) {
    YOW[i] = 1.0 / YOW[i];
  }

  for (int n = 0; n < 11; n++) {
    for (int i = 0; i < (*np); i++) {
      x[n * (*np) + i] *= YOW[i];
    }
  }
}

/*convert y[species] (mass fracs) to c[species] (molar conc) */
void
CKYTCP(double* P, double* T, double* y, double* c)
{
  double YOW = 0;
  double PWORT;

  /*Compute inverse of mean molecular wt first */
  for (int i = 0; i < 11; i++) {
    c[i] = y[i] * imw[i];
  }
  for (int i = 0; i < 11; i++) {
    YOW += c[i];
  }

  /*PW/RT (see Eq. 7) */
  PWORT = (*P) / (YOW * 8.31451e+07 * (*T));
  /*Now compute conversion */

  for (int i = 0; i < 11; i++) {
    c[i] = PWORT * y[i] * imw[i];
  }
  return;
}

/*convert y[species] (mass fracs) to c[species] (molar conc) */
void
CKYTCR(double* rho, double* T, double* y, double* c)
{
  for (int i = 0; i < 11; i++) {
    c[i] = (*rho) * y[i] * imw[i];
  }
}

/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void
CKXTY(double* x, double* y)
{
  double XW = 0; /*See Eq 4, 9 in CK Manual */
  /*Compute mean molecular wt first */
  XW += x[0] * 1.007970;  /*H */
  XW += x[1] * 2.015940;  /*H2 */
  XW += x[2] * 15.999400; /*O */
  XW += x[3] * 17.007370; /*OH */
  XW += x[4] * 18.015340; /*H2O */
  XW += x[5] * 31.998800; /*O2 */
  XW += x[6] * 33.006770; /*HO2 */
  XW += x[7] * 34.014740; /*H2O2 */
  XW += x[8] * 28.013400; /*N2 */
  XW += x[9] * 39.948000; /*AR */
  XW += x[10] * 4.002600; /*HE */
  /*Now compute conversion */
  double XWinv = 1.0 / XW;
  y[0] = x[0] * 1.007970 * XWinv;
  y[1] = x[1] * 2.015940 * XWinv;
  y[2] = x[2] * 15.999400 * XWinv;
  y[3] = x[3] * 17.007370 * XWinv;
  y[4] = x[4] * 18.015340 * XWinv;
  y[5] = x[5] * 31.998800 * XWinv;
  y[6] = x[6] * 33.006770 * XWinv;
  y[7] = x[7] * 34.014740 * XWinv;
  y[8] = x[8] * 28.013400 * XWinv;
  y[9] = x[9] * 39.948000 * XWinv;
  y[10] = x[10] * 4.002600 * XWinv;

  return;
}

/*convert x[species] (mole fracs) to c[species] (molar conc) */
void
CKXTCP(double* P, double* T, double* x, double* c)
{
  int id;                                    /*loop counter */
  double PORT = (*P) / (8.31451e+07 * (*T)); /*P/RT */

  /*Compute conversion, see Eq 10 */
  for (id = 0; id < 11; ++id) {
    c[id] = x[id] * PORT;
  }

  return;
}

/*convert x[species] (mole fracs) to c[species] (molar conc) */
void
CKXTCR(double* rho, double* T, double* x, double* c)
{
  int id;        /*loop counter */
  double XW = 0; /*See Eq 4, 11 in CK Manual */
  double ROW;
  /*Compute mean molecular wt first */
  XW += x[0] * 1.007970;  /*H */
  XW += x[1] * 2.015940;  /*H2 */
  XW += x[2] * 15.999400; /*O */
  XW += x[3] * 17.007370; /*OH */
  XW += x[4] * 18.015340; /*H2O */
  XW += x[5] * 31.998800; /*O2 */
  XW += x[6] * 33.006770; /*HO2 */
  XW += x[7] * 34.014740; /*H2O2 */
  XW += x[8] * 28.013400; /*N2 */
  XW += x[9] * 39.948000; /*AR */
  XW += x[10] * 4.002600; /*HE */
  ROW = (*rho) / XW;

  /*Compute conversion, see Eq 11 */
  for (id = 0; id < 11; ++id) {
    c[id] = x[id] * ROW;
  }

  return;
}

/*convert c[species] (molar conc) to x[species] (mole fracs) */
void
CKCTX(double* c, double* x)
{
  int id; /*loop counter */
  double sumC = 0;

  /*compute sum of c  */
  for (id = 0; id < 11; ++id) {
    sumC += c[id];
  }

  /* See Eq 13  */
  double sumCinv = 1.0 / sumC;
  for (id = 0; id < 11; ++id) {
    x[id] = c[id] * sumCinv;
  }

  return;
}

/*convert c[species] (molar conc) to y[species] (mass fracs) */
void
CKCTY(double* c, double* y)
{
  double CW = 0; /*See Eq 12 in CK Manual */
  /*compute denominator in eq 12 first */
  CW += c[0] * 1.007970;  /*H */
  CW += c[1] * 2.015940;  /*H2 */
  CW += c[2] * 15.999400; /*O */
  CW += c[3] * 17.007370; /*OH */
  CW += c[4] * 18.015340; /*H2O */
  CW += c[5] * 31.998800; /*O2 */
  CW += c[6] * 33.006770; /*HO2 */
  CW += c[7] * 34.014740; /*H2O2 */
  CW += c[8] * 28.013400; /*N2 */
  CW += c[9] * 39.948000; /*AR */
  CW += c[10] * 4.002600; /*HE */
  /*Now compute conversion */
  double CWinv = 1.0 / CW;
  y[0] = c[0] * 1.007970 * CWinv;
  y[1] = c[1] * 2.015940 * CWinv;
  y[2] = c[2] * 15.999400 * CWinv;
  y[3] = c[3] * 17.007370 * CWinv;
  y[4] = c[4] * 18.015340 * CWinv;
  y[5] = c[5] * 31.998800 * CWinv;
  y[6] = c[6] * 33.006770 * CWinv;
  y[7] = c[7] * 34.014740 * CWinv;
  y[8] = c[8] * 28.013400 * CWinv;
  y[9] = c[9] * 39.948000 * CWinv;
  y[10] = c[10] * 4.002600 * CWinv;

  return;
}

/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void
CKCPOR(double* T, double* cpor)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  cp_R(cpor, tc);
}

/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void
CKHORT(double* T, double* hort)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  speciesEnthalpy(hort, tc);
}

/*get S/R as a function of T  */
/*for all species (Eq 21) */
void
CKSOR(double* T, double* sor)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT}; /*temperature cache */
  speciesEntropy(sor, tc);
}

/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void
CKCVML(double* T, double* cvml)
{
  int id;         /*loop counter */
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  cv_R(cvml, tc);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    cvml[id] *= 8.31451e+07;
  }
}

/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void
CKCPML(double* T, double* cpml)
{
  int id;         /*loop counter */
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  cp_R(cpml, tc);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    cpml[id] *= 8.31451e+07;
  }
}

/*get internal energy as a function  */
/*of T for all species (molar units) */
void
CKUML(double* T, double* uml)
{
  int id;         /*loop counter */
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  double RT = 8.31451e+07 * tT;                       /*R*T */
  speciesInternalEnergy(uml, tc);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    uml[id] *= RT;
  }
}

/*get enthalpy as a function  */
/*of T for all species (molar units) */
void
CKHML(double* T, double* hml)
{
  int id;         /*loop counter */
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  double RT = 8.31451e+07 * tT;                       /*R*T */
  speciesEnthalpy(hml, tc);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    hml[id] *= RT;
  }
}

/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void
CKGML(double* T, double* gml)
{
  int id;         /*loop counter */
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT};         /*temperature cache */
  double RT = 8.31451e+07 * tT; /*R*T */
  gibbs(gml, tc);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    gml[id] *= RT;
  }
}

/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void
CKAML(double* T, double* aml)
{
  int id;         /*loop counter */
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT};         /*temperature cache */
  double RT = 8.31451e+07 * tT; /*R*T */
  helmholtz(aml, tc);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    aml[id] *= RT;
  }
}

/*Returns the standard-state entropies in molar units */
void
CKSML(double* T, double* sml)
{
  int id;         /*loop counter */
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT}; /*temperature cache */
  speciesEntropy(sml, tc);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    sml[id] *= 8.31451e+07;
  }
}

/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void
CKCVMS(double* T, double* cvms)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  cv_R(cvms, tc);
  /*multiply by R/molecularweight */
  cvms[0] *= 8.248767324424338e+07;  /*H */
  cvms[1] *= 4.124383662212169e+07;  /*H2 */
  cvms[2] *= 5.196763628636074e+06;  /*O */
  cvms[3] *= 4.888768810227566e+06;  /*OH */
  cvms[4] *= 4.615239012974499e+06;  /*H2O */
  cvms[5] *= 2.598381814318037e+06;  /*O2 */
  cvms[6] *= 2.519031701678171e+06;  /*HO2 */
  cvms[7] *= 2.444384405113783e+06;  /*H2O2 */
  cvms[8] *= 2.968047434442088e+06;  /*N2 */
  cvms[9] *= 2.081333233203164e+06;  /*AR */
  cvms[10] *= 2.077277269774646e+07; /*HE */
}

/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void
CKCPMS(double* T, double* cpms)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  cp_R(cpms, tc);
  /*multiply by R/molecularweight */
  cpms[0] *= 8.248767324424338e+07;  /*H */
  cpms[1] *= 4.124383662212169e+07;  /*H2 */
  cpms[2] *= 5.196763628636074e+06;  /*O */
  cpms[3] *= 4.888768810227566e+06;  /*OH */
  cpms[4] *= 4.615239012974499e+06;  /*H2O */
  cpms[5] *= 2.598381814318037e+06;  /*O2 */
  cpms[6] *= 2.519031701678171e+06;  /*HO2 */
  cpms[7] *= 2.444384405113783e+06;  /*H2O2 */
  cpms[8] *= 2.968047434442088e+06;  /*N2 */
  cpms[9] *= 2.081333233203164e+06;  /*AR */
  cpms[10] *= 2.077277269774646e+07; /*HE */
}

/*Returns internal energy in mass units (Eq 30.) */
void
CKUMS(double* T, double* ums)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  double RT = 8.31451e+07 * tT;                       /*R*T */
  speciesInternalEnergy(ums, tc);
  for (int i = 0; i < 11; i++) {
    ums[i] *= RT * imw[i];
  }
}

/*Returns enthalpy in mass units (Eq 27.) */
void
CKHMS(double* T, double* hms)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  double RT = 8.31451e+07 * tT;                       /*R*T */
  speciesEnthalpy(hms, tc);
  for (int i = 0; i < 11; i++) {
    hms[i] *= RT * imw[i];
  }
}

/*Returns enthalpy in mass units (Eq 27.) */
void
VCKHMS(int* np, double* T, double* hms)
{
  double tc[5], h[11];

  for (int i = 0; i < (*np); i++) {
    tc[0] = 0.0;
    tc[1] = T[i];
    tc[2] = T[i] * T[i];
    tc[3] = T[i] * T[i] * T[i];
    tc[4] = T[i] * T[i] * T[i] * T[i];

    speciesEnthalpy(h, tc);

    hms[0 * (*np) + i] = h[0];
    hms[1 * (*np) + i] = h[1];
    hms[2 * (*np) + i] = h[2];
    hms[3 * (*np) + i] = h[3];
    hms[4 * (*np) + i] = h[4];
    hms[5 * (*np) + i] = h[5];
    hms[6 * (*np) + i] = h[6];
    hms[7 * (*np) + i] = h[7];
    hms[8 * (*np) + i] = h[8];
    hms[9 * (*np) + i] = h[9];
    hms[10 * (*np) + i] = h[10];
  }

  for (int n = 0; n < 11; n++) {
    for (int i = 0; i < (*np); i++) {
      hms[n * (*np) + i] *= 8.31451e+07 * T[i] * imw[n];
    }
  }
}

/*Returns gibbs in mass units (Eq 31.) */
void
CKGMS(double* T, double* gms)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT};         /*temperature cache */
  double RT = 8.31451e+07 * tT; /*R*T */
  gibbs(gms, tc);
  for (int i = 0; i < 11; i++) {
    gms[i] *= RT * imw[i];
  }
}

/*Returns helmholtz in mass units (Eq 32.) */
void
CKAMS(double* T, double* ams)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT};         /*temperature cache */
  double RT = 8.31451e+07 * tT; /*R*T */
  helmholtz(ams, tc);
  for (int i = 0; i < 11; i++) {
    ams[i] *= RT * imw[i];
  }
}

/*Returns the entropies in mass units (Eq 28.) */
void
CKSMS(double* T, double* sms)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT}; /*temperature cache */
  speciesEntropy(sms, tc);
  /*multiply by R/molecularweight */
  sms[0] *= 8.248767324424338e+07;  /*H */
  sms[1] *= 4.124383662212169e+07;  /*H2 */
  sms[2] *= 5.196763628636074e+06;  /*O */
  sms[3] *= 4.888768810227566e+06;  /*OH */
  sms[4] *= 4.615239012974499e+06;  /*H2O */
  sms[5] *= 2.598381814318037e+06;  /*O2 */
  sms[6] *= 2.519031701678171e+06;  /*HO2 */
  sms[7] *= 2.444384405113783e+06;  /*H2O2 */
  sms[8] *= 2.968047434442088e+06;  /*N2 */
  sms[9] *= 2.081333233203164e+06;  /*AR */
  sms[10] *= 2.077277269774646e+07; /*HE */
}

/*Returns the mean specific heat at CP (Eq. 33) */
void
CKCPBL(double* T, double* x, double* cpbl)
{
  int id; /*loop counter */
  double result = 0;
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  double cpor[11];                                    /* temporary storage */
  cp_R(cpor, tc);

  /*perform dot product */
  for (id = 0; id < 11; ++id) {
    result += x[id] * cpor[id];
  }

  *cpbl = result * 8.31451e+07;
}

/*Returns the mean specific heat at CP (Eq. 34) */
void
CKCPBS(double* T, double* y, double* cpbs)
{
  double result = 0;
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  double cpor[11], tresult[11];                       /* temporary storage */
  cp_R(cpor, tc);
  for (int i = 0; i < 11; i++) {
    tresult[i] = cpor[i] * y[i] * imw[i];
  }
  for (int i = 0; i < 11; i++) {
    result += tresult[i];
  }

  *cpbs = result * 8.31451e+07;
}

/*Returns the mean specific heat at CV (Eq. 35) */
void
CKCVBL(double* T, double* x, double* cvbl)
{
  int id; /*loop counter */
  double result = 0;
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  double cvor[11];                                    /* temporary storage */
  cv_R(cvor, tc);

  /*perform dot product */
  for (id = 0; id < 11; ++id) {
    result += x[id] * cvor[id];
  }

  *cvbl = result * 8.31451e+07;
}

/*Returns the mean specific heat at CV (Eq. 36) */
void
CKCVBS(double* T, double* y, double* cvbs)
{
  double result = 0;
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  double cvor[11];                                    /* temporary storage */
  cv_R(cvor, tc);
  /*multiply by y/molecularweight */
  result += cvor[0] * y[0] * imw[0];    /*H */
  result += cvor[1] * y[1] * imw[1];    /*H2 */
  result += cvor[2] * y[2] * imw[2];    /*O */
  result += cvor[3] * y[3] * imw[3];    /*OH */
  result += cvor[4] * y[4] * imw[4];    /*H2O */
  result += cvor[5] * y[5] * imw[5];    /*O2 */
  result += cvor[6] * y[6] * imw[6];    /*HO2 */
  result += cvor[7] * y[7] * imw[7];    /*H2O2 */
  result += cvor[8] * y[8] * imw[8];    /*N2 */
  result += cvor[9] * y[9] * imw[9];    /*AR */
  result += cvor[10] * y[10] * imw[10]; /*HE */

  *cvbs = result * 8.31451e+07;
}

/*Returns the mean enthalpy of the mixture in molar units */
void
CKHBML(double* T, double* x, double* hbml)
{
  int id; /*loop counter */
  double result = 0;
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  double hml[11];                                     /* temporary storage */
  double RT = 8.31451e+07 * tT;                       /*R*T */
  speciesEnthalpy(hml, tc);

  /*perform dot product */
  for (id = 0; id < 11; ++id) {
    result += x[id] * hml[id];
  }

  *hbml = result * RT;
}

/*Returns mean enthalpy of mixture in mass units */
void
CKHBMS(double* T, double* y, double* hbms)
{
  double result = 0;
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  double hml[11], tmp[11];                            /* temporary storage */
  double RT = 8.31451e+07 * tT;                       /*R*T */
  speciesEnthalpy(hml, tc);
  int id;
  for (id = 0; id < 11; ++id) {
    tmp[id] = y[id] * hml[id] * imw[id];
  }
  for (id = 0; id < 11; ++id) {
    result += tmp[id];
  }

  *hbms = result * RT;
}

/*get mean internal energy in molar units */
void
CKUBML(double* T, double* x, double* ubml)
{
  int id; /*loop counter */
  double result = 0;
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  double uml[11];               /* temporary energy array */
  double RT = 8.31451e+07 * tT; /*R*T */
  speciesInternalEnergy(uml, tc);

  /*perform dot product */
  for (id = 0; id < 11; ++id) {
    result += x[id] * uml[id];
  }

  *ubml = result * RT;
}

/*get mean internal energy in mass units */
void
CKUBMS(double* T, double* y, double* ubms)
{
  double result = 0;
  double tT = *T; /*temporary temperature */
  double tc[] = {
    0, tT, tT * tT, tT * tT * tT, tT * tT * tT * tT}; /*temperature cache */
  double ums[11];               /* temporary energy array */
  double RT = 8.31451e+07 * tT; /*R*T */
  speciesInternalEnergy(ums, tc);
  /*perform dot product + scaling by wt */
  result += y[0] * ums[0] * imw[0];    /*H */
  result += y[1] * ums[1] * imw[1];    /*H2 */
  result += y[2] * ums[2] * imw[2];    /*O */
  result += y[3] * ums[3] * imw[3];    /*OH */
  result += y[4] * ums[4] * imw[4];    /*H2O */
  result += y[5] * ums[5] * imw[5];    /*O2 */
  result += y[6] * ums[6] * imw[6];    /*HO2 */
  result += y[7] * ums[7] * imw[7];    /*H2O2 */
  result += y[8] * ums[8] * imw[8];    /*N2 */
  result += y[9] * ums[9] * imw[9];    /*AR */
  result += y[10] * ums[10] * imw[10]; /*HE */

  *ubms = result * RT;
}

/*get mixture entropy in molar units */
void
CKSBML(double* P, double* T, double* x, double* sbml)
{
  int id; /*loop counter */
  double result = 0;
  /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
  double logPratio = log(*P / 1013250.0);
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT}; /*temperature cache */
  double sor[11];       /* temporary storage */
  speciesEntropy(sor, tc);

  /*Compute Eq 42 */
  for (id = 0; id < 11; ++id) {
    result += x[id] * (sor[id] - log((x[id] + 1e-100)) - logPratio);
  }

  *sbml = result * 8.31451e+07;
}

/*get mixture entropy in mass units */
void
CKSBMS(double* P, double* T, double* y, double* sbms)
{
  double result = 0;
  /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
  double logPratio = log(*P / 1013250.0);
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT}; /*temperature cache */
  double sor[11];       /* temporary storage */
  double x[11];         /* need a ytx conversion */
  double YOW = 0;       /*See Eq 4, 6 in CK Manual */
  /*Compute inverse of mean molecular wt first */
  YOW += y[0] * imw[0];   /*H */
  YOW += y[1] * imw[1];   /*H2 */
  YOW += y[2] * imw[2];   /*O */
  YOW += y[3] * imw[3];   /*OH */
  YOW += y[4] * imw[4];   /*H2O */
  YOW += y[5] * imw[5];   /*O2 */
  YOW += y[6] * imw[6];   /*HO2 */
  YOW += y[7] * imw[7];   /*H2O2 */
  YOW += y[8] * imw[8];   /*N2 */
  YOW += y[9] * imw[9];   /*AR */
  YOW += y[10] * imw[10]; /*HE */
  /*Now compute y to x conversion */
  x[0] = y[0] / (1.007970 * YOW);
  x[1] = y[1] / (2.015940 * YOW);
  x[2] = y[2] / (15.999400 * YOW);
  x[3] = y[3] / (17.007370 * YOW);
  x[4] = y[4] / (18.015340 * YOW);
  x[5] = y[5] / (31.998800 * YOW);
  x[6] = y[6] / (33.006770 * YOW);
  x[7] = y[7] / (34.014740 * YOW);
  x[8] = y[8] / (28.013400 * YOW);
  x[9] = y[9] / (39.948000 * YOW);
  x[10] = y[10] / (4.002600 * YOW);
  speciesEntropy(sor, tc);
  /*Perform computation in Eq 42 and 43 */
  result += x[0] * (sor[0] - log((x[0] + 1e-100)) - logPratio);
  result += x[1] * (sor[1] - log((x[1] + 1e-100)) - logPratio);
  result += x[2] * (sor[2] - log((x[2] + 1e-100)) - logPratio);
  result += x[3] * (sor[3] - log((x[3] + 1e-100)) - logPratio);
  result += x[4] * (sor[4] - log((x[4] + 1e-100)) - logPratio);
  result += x[5] * (sor[5] - log((x[5] + 1e-100)) - logPratio);
  result += x[6] * (sor[6] - log((x[6] + 1e-100)) - logPratio);
  result += x[7] * (sor[7] - log((x[7] + 1e-100)) - logPratio);
  result += x[8] * (sor[8] - log((x[8] + 1e-100)) - logPratio);
  result += x[9] * (sor[9] - log((x[9] + 1e-100)) - logPratio);
  result += x[10] * (sor[10] - log((x[10] + 1e-100)) - logPratio);
  /*Scale by R/W */
  *sbms = result * 8.31451e+07 * YOW;
}

/*Returns mean gibbs free energy in molar units */
void
CKGBML(double* P, double* T, double* x, double* gbml)
{
  int id; /*loop counter */
  double result = 0;
  /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
  double logPratio = log(*P / 1013250.0);
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT};         /*temperature cache */
  double RT = 8.31451e+07 * tT; /*R*T */
  double gort[11];              /* temporary storage */
  /*Compute g/RT */
  gibbs(gort, tc);

  /*Compute Eq 44 */
  for (id = 0; id < 11; ++id) {
    result += x[id] * (gort[id] + log((x[id] + 1e-100)) + logPratio);
  }

  *gbml = result * RT;
}

/*Returns mixture gibbs free energy in mass units */
void
CKGBMS(double* P, double* T, double* y, double* gbms)
{
  double result = 0;
  /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
  double logPratio = log(*P / 1013250.0);
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT};         /*temperature cache */
  double RT = 8.31451e+07 * tT; /*R*T */
  double gort[11];              /* temporary storage */
  double x[11];                 /* need a ytx conversion */
  double YOW = 0;               /*To hold 1/molecularweight */
  /*Compute inverse of mean molecular wt first */
  YOW += y[0] * imw[0];   /*H */
  YOW += y[1] * imw[1];   /*H2 */
  YOW += y[2] * imw[2];   /*O */
  YOW += y[3] * imw[3];   /*OH */
  YOW += y[4] * imw[4];   /*H2O */
  YOW += y[5] * imw[5];   /*O2 */
  YOW += y[6] * imw[6];   /*HO2 */
  YOW += y[7] * imw[7];   /*H2O2 */
  YOW += y[8] * imw[8];   /*N2 */
  YOW += y[9] * imw[9];   /*AR */
  YOW += y[10] * imw[10]; /*HE */
  /*Now compute y to x conversion */
  x[0] = y[0] / (1.007970 * YOW);
  x[1] = y[1] / (2.015940 * YOW);
  x[2] = y[2] / (15.999400 * YOW);
  x[3] = y[3] / (17.007370 * YOW);
  x[4] = y[4] / (18.015340 * YOW);
  x[5] = y[5] / (31.998800 * YOW);
  x[6] = y[6] / (33.006770 * YOW);
  x[7] = y[7] / (34.014740 * YOW);
  x[8] = y[8] / (28.013400 * YOW);
  x[9] = y[9] / (39.948000 * YOW);
  x[10] = y[10] / (4.002600 * YOW);
  gibbs(gort, tc);
  /*Perform computation in Eq 44 */
  result += x[0] * (gort[0] + log((x[0] + 1e-100)) + logPratio);
  result += x[1] * (gort[1] + log((x[1] + 1e-100)) + logPratio);
  result += x[2] * (gort[2] + log((x[2] + 1e-100)) + logPratio);
  result += x[3] * (gort[3] + log((x[3] + 1e-100)) + logPratio);
  result += x[4] * (gort[4] + log((x[4] + 1e-100)) + logPratio);
  result += x[5] * (gort[5] + log((x[5] + 1e-100)) + logPratio);
  result += x[6] * (gort[6] + log((x[6] + 1e-100)) + logPratio);
  result += x[7] * (gort[7] + log((x[7] + 1e-100)) + logPratio);
  result += x[8] * (gort[8] + log((x[8] + 1e-100)) + logPratio);
  result += x[9] * (gort[9] + log((x[9] + 1e-100)) + logPratio);
  result += x[10] * (gort[10] + log((x[10] + 1e-100)) + logPratio);
  /*Scale by RT/W */
  *gbms = result * RT * YOW;
}

/*Returns mean helmholtz free energy in molar units */
void
CKABML(double* P, double* T, double* x, double* abml)
{
  int id; /*loop counter */
  double result = 0;
  /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
  double logPratio = log(*P / 1013250.0);
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT};         /*temperature cache */
  double RT = 8.31451e+07 * tT; /*R*T */
  double aort[11];              /* temporary storage */
  /*Compute g/RT */
  helmholtz(aort, tc);

  /*Compute Eq 44 */
  for (id = 0; id < 11; ++id) {
    result += x[id] * (aort[id] + log((x[id] + 1e-100)) + logPratio);
  }

  *abml = result * RT;
}

/*Returns mixture helmholtz free energy in mass units */
void
CKABMS(double* P, double* T, double* y, double* abms)
{
  double result = 0;
  /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
  double logPratio = log(*P / 1013250.0);
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT};         /*temperature cache */
  double RT = 8.31451e+07 * tT; /*R*T */
  double aort[11];              /* temporary storage */
  double x[11];                 /* need a ytx conversion */
  double YOW = 0;               /*To hold 1/molecularweight */
  /*Compute inverse of mean molecular wt first */
  YOW += y[0] * imw[0];   /*H */
  YOW += y[1] * imw[1];   /*H2 */
  YOW += y[2] * imw[2];   /*O */
  YOW += y[3] * imw[3];   /*OH */
  YOW += y[4] * imw[4];   /*H2O */
  YOW += y[5] * imw[5];   /*O2 */
  YOW += y[6] * imw[6];   /*HO2 */
  YOW += y[7] * imw[7];   /*H2O2 */
  YOW += y[8] * imw[8];   /*N2 */
  YOW += y[9] * imw[9];   /*AR */
  YOW += y[10] * imw[10]; /*HE */
  /*Now compute y to x conversion */
  x[0] = y[0] / (1.007970 * YOW);
  x[1] = y[1] / (2.015940 * YOW);
  x[2] = y[2] / (15.999400 * YOW);
  x[3] = y[3] / (17.007370 * YOW);
  x[4] = y[4] / (18.015340 * YOW);
  x[5] = y[5] / (31.998800 * YOW);
  x[6] = y[6] / (33.006770 * YOW);
  x[7] = y[7] / (34.014740 * YOW);
  x[8] = y[8] / (28.013400 * YOW);
  x[9] = y[9] / (39.948000 * YOW);
  x[10] = y[10] / (4.002600 * YOW);
  helmholtz(aort, tc);
  /*Perform computation in Eq 44 */
  result += x[0] * (aort[0] + log((x[0] + 1e-100)) + logPratio);
  result += x[1] * (aort[1] + log((x[1] + 1e-100)) + logPratio);
  result += x[2] * (aort[2] + log((x[2] + 1e-100)) + logPratio);
  result += x[3] * (aort[3] + log((x[3] + 1e-100)) + logPratio);
  result += x[4] * (aort[4] + log((x[4] + 1e-100)) + logPratio);
  result += x[5] * (aort[5] + log((x[5] + 1e-100)) + logPratio);
  result += x[6] * (aort[6] + log((x[6] + 1e-100)) + logPratio);
  result += x[7] * (aort[7] + log((x[7] + 1e-100)) + logPratio);
  result += x[8] * (aort[8] + log((x[8] + 1e-100)) + logPratio);
  result += x[9] * (aort[9] + log((x[9] + 1e-100)) + logPratio);
  result += x[10] * (aort[10] + log((x[10] + 1e-100)) + logPratio);
  /*Scale by RT/W */
  *abms = result * RT * YOW;
}

/*compute the production rate for each species */
void
CKWC(double* T, double* C, double* wdot)
{
  int id; /*loop counter */

  /*convert to SI */
  for (id = 0; id < 11; ++id) {
    C[id] *= 1.0e6;
  }

  /*convert to chemkin units */
  productionRate(wdot, C, *T);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    C[id] *= 1.0e-6;
    wdot[id] *= 1.0e-6;
  }
}

/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void
CKWYP(double* P, double* T, double* y, double* wdot)
{
  int id;       /*loop counter */
  double c[11]; /*temporary storage */
  double YOW = 0;
  double PWORT;
  /*Compute inverse of mean molecular wt first */
  YOW += y[0] * imw[0];   /*H */
  YOW += y[1] * imw[1];   /*H2 */
  YOW += y[2] * imw[2];   /*O */
  YOW += y[3] * imw[3];   /*OH */
  YOW += y[4] * imw[4];   /*H2O */
  YOW += y[5] * imw[5];   /*O2 */
  YOW += y[6] * imw[6];   /*HO2 */
  YOW += y[7] * imw[7];   /*H2O2 */
  YOW += y[8] * imw[8];   /*N2 */
  YOW += y[9] * imw[9];   /*AR */
  YOW += y[10] * imw[10]; /*HE */
  /*PW/RT (see Eq. 7) */
  PWORT = (*P) / (YOW * 8.31451e+07 * (*T));
  /*multiply by 1e6 so c goes to SI */
  PWORT *= 1e6;
  /*Now compute conversion (and go to SI) */
  c[0] = PWORT * y[0] * imw[0];
  c[1] = PWORT * y[1] * imw[1];
  c[2] = PWORT * y[2] * imw[2];
  c[3] = PWORT * y[3] * imw[3];
  c[4] = PWORT * y[4] * imw[4];
  c[5] = PWORT * y[5] * imw[5];
  c[6] = PWORT * y[6] * imw[6];
  c[7] = PWORT * y[7] * imw[7];
  c[8] = PWORT * y[8] * imw[8];
  c[9] = PWORT * y[9] * imw[9];
  c[10] = PWORT * y[10] * imw[10];

  /*convert to chemkin units */
  productionRate(wdot, c, *T);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void
CKWXP(double* P, double* T, double* x, double* wdot)
{
  int id;       /*loop counter */
  double c[11]; /*temporary storage */
  double PORT =
    1e6 * (*P) / (8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

  /*Compute conversion, see Eq 10 */
  for (id = 0; id < 11; ++id) {
    c[id] = x[id] * PORT;
  }

  /*convert to chemkin units */
  productionRate(wdot, c, *T);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void
CKWYR(double* rho, double* T, double* y, double* wdot)
{
  int id;       /*loop counter */
  double c[11]; /*temporary storage */
  /*See Eq 8 with an extra 1e6 so c goes to SI */
  c[0] = 1e6 * (*rho) * y[0] * imw[0];
  c[1] = 1e6 * (*rho) * y[1] * imw[1];
  c[2] = 1e6 * (*rho) * y[2] * imw[2];
  c[3] = 1e6 * (*rho) * y[3] * imw[3];
  c[4] = 1e6 * (*rho) * y[4] * imw[4];
  c[5] = 1e6 * (*rho) * y[5] * imw[5];
  c[6] = 1e6 * (*rho) * y[6] * imw[6];
  c[7] = 1e6 * (*rho) * y[7] * imw[7];
  c[8] = 1e6 * (*rho) * y[8] * imw[8];
  c[9] = 1e6 * (*rho) * y[9] * imw[9];
  c[10] = 1e6 * (*rho) * y[10] * imw[10];

  /*call productionRate */
  productionRate(wdot, c, *T);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void
VCKWYR(int* np, double* rho, double* T, double* y, double* wdot)
{
  double c[11 * (*np)]; /*temporary storage */
  /*See Eq 8 with an extra 1e6 so c goes to SI */
  for (int n = 0; n < 11; n++) {
    for (int i = 0; i < (*np); i++) {
      c[n * (*np) + i] = 1.0e6 * rho[i] * y[n * (*np) + i] * imw[n];
    }
  }

  /*call productionRate */
  vproductionRate(*np, wdot, c, T);

  /*convert to chemkin units */
  for (int i = 0; i < 11 * (*np); i++) {
    wdot[i] *= 1.0e-6;
  }
}

/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void
CKWXR(double* rho, double* T, double* x, double* wdot)
{
  int id;        /*loop counter */
  double c[11];  /*temporary storage */
  double XW = 0; /*See Eq 4, 11 in CK Manual */
  double ROW;
  /*Compute mean molecular wt first */
  XW += x[0] * 1.007970;  /*H */
  XW += x[1] * 2.015940;  /*H2 */
  XW += x[2] * 15.999400; /*O */
  XW += x[3] * 17.007370; /*OH */
  XW += x[4] * 18.015340; /*H2O */
  XW += x[5] * 31.998800; /*O2 */
  XW += x[6] * 33.006770; /*HO2 */
  XW += x[7] * 34.014740; /*H2O2 */
  XW += x[8] * 28.013400; /*N2 */
  XW += x[9] * 39.948000; /*AR */
  XW += x[10] * 4.002600; /*HE */
  /*Extra 1e6 factor to take c to SI */
  ROW = 1e6 * (*rho) / XW;

  /*Compute conversion, see Eq 11 */
  for (id = 0; id < 11; ++id) {
    c[id] = x[id] * ROW;
  }

  /*convert to chemkin units */
  productionRate(wdot, c, *T);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

/*Returns the rate of progress for each reaction */
void
CKQC(double* T, double* C, double* qdot)
{
  int id; /*loop counter */

  /*convert to SI */
  for (id = 0; id < 11; ++id) {
    C[id] *= 1.0e6;
  }

  /*convert to chemkin units */
  progressRate(qdot, C, *T);

  /*convert to chemkin units */
  for (id = 0; id < 11; ++id) {
    C[id] *= 1.0e-6;
  }

  for (id = 0; id < 29; ++id) {
    qdot[id] *= 1.0e-6;
  }
}

/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void
CKKFKR(double* P, double* T, double* x, double* q_f, double* q_r)
{
  int id;       /*loop counter */
  double c[11]; /*temporary storage */
  double PORT =
    1e6 * (*P) / (8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

  /*Compute conversion, see Eq 10 */
  for (id = 0; id < 11; ++id) {
    c[id] = x[id] * PORT;
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
void
CKQYP(double* P, double* T, double* y, double* qdot)
{
  int id;       /*loop counter */
  double c[11]; /*temporary storage */
  double YOW = 0;
  double PWORT;
  /*Compute inverse of mean molecular wt first */
  YOW += y[0] * imw[0];   /*H */
  YOW += y[1] * imw[1];   /*H2 */
  YOW += y[2] * imw[2];   /*O */
  YOW += y[3] * imw[3];   /*OH */
  YOW += y[4] * imw[4];   /*H2O */
  YOW += y[5] * imw[5];   /*O2 */
  YOW += y[6] * imw[6];   /*HO2 */
  YOW += y[7] * imw[7];   /*H2O2 */
  YOW += y[8] * imw[8];   /*N2 */
  YOW += y[9] * imw[9];   /*AR */
  YOW += y[10] * imw[10]; /*HE */
  /*PW/RT (see Eq. 7) */
  PWORT = (*P) / (YOW * 8.31451e+07 * (*T));
  /*multiply by 1e6 so c goes to SI */
  PWORT *= 1e6;
  /*Now compute conversion (and go to SI) */
  c[0] = PWORT * y[0] * imw[0];
  c[1] = PWORT * y[1] * imw[1];
  c[2] = PWORT * y[2] * imw[2];
  c[3] = PWORT * y[3] * imw[3];
  c[4] = PWORT * y[4] * imw[4];
  c[5] = PWORT * y[5] * imw[5];
  c[6] = PWORT * y[6] * imw[6];
  c[7] = PWORT * y[7] * imw[7];
  c[8] = PWORT * y[8] * imw[8];
  c[9] = PWORT * y[9] * imw[9];
  c[10] = PWORT * y[10] * imw[10];

  /*convert to chemkin units */
  progressRate(qdot, c, *T);

  /*convert to chemkin units */
  for (id = 0; id < 29; ++id) {
    qdot[id] *= 1.0e-6;
  }
}

/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void
CKQXP(double* P, double* T, double* x, double* qdot)
{
  int id;       /*loop counter */
  double c[11]; /*temporary storage */
  double PORT =
    1e6 * (*P) / (8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

  /*Compute conversion, see Eq 10 */
  for (id = 0; id < 11; ++id) {
    c[id] = x[id] * PORT;
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
void
CKQYR(double* rho, double* T, double* y, double* qdot)
{
  int id;       /*loop counter */
  double c[11]; /*temporary storage */
  /*See Eq 8 with an extra 1e6 so c goes to SI */
  c[0] = 1e6 * (*rho) * y[0] * imw[0];
  c[1] = 1e6 * (*rho) * y[1] * imw[1];
  c[2] = 1e6 * (*rho) * y[2] * imw[2];
  c[3] = 1e6 * (*rho) * y[3] * imw[3];
  c[4] = 1e6 * (*rho) * y[4] * imw[4];
  c[5] = 1e6 * (*rho) * y[5] * imw[5];
  c[6] = 1e6 * (*rho) * y[6] * imw[6];
  c[7] = 1e6 * (*rho) * y[7] * imw[7];
  c[8] = 1e6 * (*rho) * y[8] * imw[8];
  c[9] = 1e6 * (*rho) * y[9] * imw[9];
  c[10] = 1e6 * (*rho) * y[10] * imw[10];

  /*call progressRate */
  progressRate(qdot, c, *T);

  /*convert to chemkin units */
  for (id = 0; id < 29; ++id) {
    qdot[id] *= 1.0e-6;
  }
}

/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void
CKQXR(double* rho, double* T, double* x, double* qdot)
{
  int id;        /*loop counter */
  double c[11];  /*temporary storage */
  double XW = 0; /*See Eq 4, 11 in CK Manual */
  double ROW;
  /*Compute mean molecular wt first */
  XW += x[0] * 1.007970;  /*H */
  XW += x[1] * 2.015940;  /*H2 */
  XW += x[2] * 15.999400; /*O */
  XW += x[3] * 17.007370; /*OH */
  XW += x[4] * 18.015340; /*H2O */
  XW += x[5] * 31.998800; /*O2 */
  XW += x[6] * 33.006770; /*HO2 */
  XW += x[7] * 34.014740; /*H2O2 */
  XW += x[8] * 28.013400; /*N2 */
  XW += x[9] * 39.948000; /*AR */
  XW += x[10] * 4.002600; /*HE */
  /*Extra 1e6 factor to take c to SI */
  ROW = 1e6 * (*rho) / XW;

  /*Compute conversion, see Eq 11 */
  for (id = 0; id < 11; ++id) {
    c[id] = x[id] * ROW;
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
void
CKNU(int* kdim, int* nuki)
{
  int id; /*loop counter */
  int kd = (*kdim);
  /*Zero nuki */
  for (id = 0; id < 11 * kd; ++id) {
    nuki[id] = 0;
  }

  /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
  nuki[0 * kd + 0] += -1;
  nuki[5 * kd + 0] += -1;
  nuki[6 * kd + 0] += +1;

  /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
  nuki[7 * kd + 1] += -1;
  nuki[3 * kd + 1] += +1;
  nuki[3 * kd + 1] += +1;

  /*reaction 3: H2 + M <=> H + H + M */
  nuki[1 * kd + 2] += -1;
  nuki[0 * kd + 2] += +1;
  nuki[0 * kd + 2] += +1;

  /*reaction 4: O + O + M <=> O2 + M */
  nuki[2 * kd + 3] += -1;
  nuki[2 * kd + 3] += -1;
  nuki[5 * kd + 3] += +1;

  /*reaction 5: O + H + M <=> OH + M */
  nuki[2 * kd + 4] += -1;
  nuki[0 * kd + 4] += -1;
  nuki[3 * kd + 4] += +1;

  /*reaction 6: H2O + M <=> H + OH + M */
  nuki[4 * kd + 5] += -1;
  nuki[0 * kd + 5] += +1;
  nuki[3 * kd + 5] += +1;

  /*reaction 7: O + OH + M <=> HO2 + M */
  nuki[2 * kd + 6] += -1;
  nuki[3 * kd + 6] += -1;
  nuki[6 * kd + 6] += +1;

  /*reaction 8: H + O2 <=> O + OH */
  nuki[0 * kd + 7] += -1;
  nuki[5 * kd + 7] += -1;
  nuki[2 * kd + 7] += +1;
  nuki[3 * kd + 7] += +1;

  /*reaction 9: O + H2 <=> H + OH */
  nuki[2 * kd + 8] += -1;
  nuki[1 * kd + 8] += -1;
  nuki[0 * kd + 8] += +1;
  nuki[3 * kd + 8] += +1;

  /*reaction 10: O + H2 <=> H + OH */
  nuki[2 * kd + 9] += -1;
  nuki[1 * kd + 9] += -1;
  nuki[0 * kd + 9] += +1;
  nuki[3 * kd + 9] += +1;

  /*reaction 11: H2 + OH <=> H2O + H */
  nuki[1 * kd + 10] += -1;
  nuki[3 * kd + 10] += -1;
  nuki[4 * kd + 10] += +1;
  nuki[0 * kd + 10] += +1;

  /*reaction 12: OH + OH <=> O + H2O */
  nuki[3 * kd + 11] += -1;
  nuki[3 * kd + 11] += -1;
  nuki[2 * kd + 11] += +1;
  nuki[4 * kd + 11] += +1;

  /*reaction 13: H2 + AR <=> H + H + AR */
  nuki[1 * kd + 12] += -1;
  nuki[9 * kd + 12] += -1;
  nuki[0 * kd + 12] += +1;
  nuki[0 * kd + 12] += +1;
  nuki[9 * kd + 12] += +1;

  /*reaction 14: H2 + HE <=> H + H + HE */
  nuki[1 * kd + 13] += -1;
  nuki[10 * kd + 13] += -1;
  nuki[0 * kd + 13] += +1;
  nuki[0 * kd + 13] += +1;
  nuki[10 * kd + 13] += +1;

  /*reaction 15: O + O + AR <=> O2 + AR */
  nuki[2 * kd + 14] += -1;
  nuki[2 * kd + 14] += -1;
  nuki[9 * kd + 14] += -1;
  nuki[5 * kd + 14] += +1;
  nuki[9 * kd + 14] += +1;

  /*reaction 16: O + O + HE <=> O2 + HE */
  nuki[2 * kd + 15] += -1;
  nuki[2 * kd + 15] += -1;
  nuki[10 * kd + 15] += -1;
  nuki[5 * kd + 15] += +1;
  nuki[10 * kd + 15] += +1;

  /*reaction 17: H2O + H2O <=> H + OH + H2O */
  nuki[4 * kd + 16] += -1;
  nuki[4 * kd + 16] += -1;
  nuki[0 * kd + 16] += +1;
  nuki[3 * kd + 16] += +1;
  nuki[4 * kd + 16] += +1;

  /*reaction 18: HO2 + H <=> H2 + O2 */
  nuki[6 * kd + 17] += -1;
  nuki[0 * kd + 17] += -1;
  nuki[1 * kd + 17] += +1;
  nuki[5 * kd + 17] += +1;

  /*reaction 19: HO2 + H <=> OH + OH */
  nuki[6 * kd + 18] += -1;
  nuki[0 * kd + 18] += -1;
  nuki[3 * kd + 18] += +1;
  nuki[3 * kd + 18] += +1;

  /*reaction 20: HO2 + O <=> O2 + OH */
  nuki[6 * kd + 19] += -1;
  nuki[2 * kd + 19] += -1;
  nuki[5 * kd + 19] += +1;
  nuki[3 * kd + 19] += +1;

  /*reaction 21: HO2 + OH <=> H2O + O2 */
  nuki[6 * kd + 20] += -1;
  nuki[3 * kd + 20] += -1;
  nuki[4 * kd + 20] += +1;
  nuki[5 * kd + 20] += +1;

  /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
  nuki[6 * kd + 21] += -1;
  nuki[6 * kd + 21] += -1;
  nuki[7 * kd + 21] += +1;
  nuki[5 * kd + 21] += +1;

  /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
  nuki[6 * kd + 22] += -1;
  nuki[6 * kd + 22] += -1;
  nuki[7 * kd + 22] += +1;
  nuki[5 * kd + 22] += +1;

  /*reaction 24: H2O2 + H <=> H2O + OH */
  nuki[7 * kd + 23] += -1;
  nuki[0 * kd + 23] += -1;
  nuki[4 * kd + 23] += +1;
  nuki[3 * kd + 23] += +1;

  /*reaction 25: H2O2 + H <=> HO2 + H2 */
  nuki[7 * kd + 24] += -1;
  nuki[0 * kd + 24] += -1;
  nuki[6 * kd + 24] += +1;
  nuki[1 * kd + 24] += +1;

  /*reaction 26: H2O2 + O <=> OH + HO2 */
  nuki[7 * kd + 25] += -1;
  nuki[2 * kd + 25] += -1;
  nuki[3 * kd + 25] += +1;
  nuki[6 * kd + 25] += +1;

  /*reaction 27: H2O2 + OH <=> HO2 + H2O */
  nuki[7 * kd + 26] += -1;
  nuki[3 * kd + 26] += -1;
  nuki[6 * kd + 26] += +1;
  nuki[4 * kd + 26] += +1;

  /*reaction 28: H2O2 + OH <=> HO2 + H2O */
  nuki[7 * kd + 27] += -1;
  nuki[3 * kd + 27] += -1;
  nuki[6 * kd + 27] += +1;
  nuki[4 * kd + 27] += +1;

  /*reaction 29: HO2 + H <=> O + H2O */
  nuki[6 * kd + 28] += -1;
  nuki[0 * kd + 28] += -1;
  nuki[2 * kd + 28] += +1;
  nuki[4 * kd + 28] += +1;
}

/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void
CKNCF(int* mdim, int* ncf)
{
  int id; /*loop counter */
  int kd = (*mdim);
  /*Zero ncf */
  for (id = 0; id < kd * 11; ++id) {
    ncf[id] = 0;
  }

  /*H */
  ncf[0 * kd + 0] = 1; /*H */

  /*H2 */
  ncf[1 * kd + 0] = 2; /*H */

  /*O */
  ncf[2 * kd + 1] = 1; /*O */

  /*OH */
  ncf[3 * kd + 1] = 1; /*O */
  ncf[3 * kd + 0] = 1; /*H */

  /*H2O */
  ncf[4 * kd + 0] = 2; /*H */
  ncf[4 * kd + 1] = 1; /*O */

  /*O2 */
  ncf[5 * kd + 1] = 2; /*O */

  /*HO2 */
  ncf[6 * kd + 0] = 1; /*H */
  ncf[6 * kd + 1] = 2; /*O */

  /*H2O2 */
  ncf[7 * kd + 0] = 2; /*H */
  ncf[7 * kd + 1] = 2; /*O */

  /*N2 */
  ncf[8 * kd + 2] = 2; /*N */

  /*AR */
  ncf[9 * kd + 3] = 1; /*AR */

  /*HE */
  ncf[10 * kd + 4] = 1; /*HE */
}

/*Returns the arrehenius coefficients  */
/*for all reactions */
void
CKABE(double* a, double* b, double* e)
{
  for (int i = 0; i < 29; ++i) {
    a[i] = fwd_A[i];
    b[i] = fwd_beta[i];
    e[i] = fwd_Ea[i];
  }

  return;
}

/*Returns the equil constants for each reaction */
void
CKEQC(double* T, double* C, double* eqcon)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT}; /*temperature cache */
  double gort[11];      /* temporary storage */

  /*compute the Gibbs free energy */
  gibbs(gort, tc);

  /*compute the equilibrium constants */
  equilibriumConstants(eqcon, gort, tT);

  /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
  eqcon[0] *= 1e+06;

  /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
  eqcon[1] *= 1e-06;

  /*reaction 3: H2 + M <=> H + H + M */
  eqcon[2] *= 1e-06;

  /*reaction 4: O + O + M <=> O2 + M */
  eqcon[3] *= 1e+06;

  /*reaction 5: O + H + M <=> OH + M */
  eqcon[4] *= 1e+06;

  /*reaction 6: H2O + M <=> H + OH + M */
  eqcon[5] *= 1e-06;

  /*reaction 7: O + OH + M <=> HO2 + M */
  eqcon[6] *= 1e+06;

  /*reaction 8: H + O2 <=> O + OH */
  /*eqcon[7] *= 1;  */

  /*reaction 9: O + H2 <=> H + OH */
  /*eqcon[8] *= 1;  */

  /*reaction 10: O + H2 <=> H + OH */
  /*eqcon[9] *= 1;  */

  /*reaction 11: H2 + OH <=> H2O + H */
  /*eqcon[10] *= 1;  */

  /*reaction 12: OH + OH <=> O + H2O */
  /*eqcon[11] *= 1;  */

  /*reaction 13: H2 + AR <=> H + H + AR */
  eqcon[12] *= 1e-06;

  /*reaction 14: H2 + HE <=> H + H + HE */
  eqcon[13] *= 1e-06;

  /*reaction 15: O + O + AR <=> O2 + AR */
  eqcon[14] *= 1e+06;

  /*reaction 16: O + O + HE <=> O2 + HE */
  eqcon[15] *= 1e+06;

  /*reaction 17: H2O + H2O <=> H + OH + H2O */
  eqcon[16] *= 1e-06;

  /*reaction 18: HO2 + H <=> H2 + O2 */
  /*eqcon[17] *= 1;  */

  /*reaction 19: HO2 + H <=> OH + OH */
  /*eqcon[18] *= 1;  */

  /*reaction 20: HO2 + O <=> O2 + OH */
  /*eqcon[19] *= 1;  */

  /*reaction 21: HO2 + OH <=> H2O + O2 */
  /*eqcon[20] *= 1;  */

  /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
  /*eqcon[21] *= 1;  */

  /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
  /*eqcon[22] *= 1;  */

  /*reaction 24: H2O2 + H <=> H2O + OH */
  /*eqcon[23] *= 1;  */

  /*reaction 25: H2O2 + H <=> HO2 + H2 */
  /*eqcon[24] *= 1;  */

  /*reaction 26: H2O2 + O <=> OH + HO2 */
  /*eqcon[25] *= 1;  */

  /*reaction 27: H2O2 + OH <=> HO2 + H2O */
  /*eqcon[26] *= 1;  */

  /*reaction 28: H2O2 + OH <=> HO2 + H2O */
  /*eqcon[27] *= 1;  */

  /*reaction 29: HO2 + H <=> O + H2O */
  /*eqcon[28] *= 1;  */
}

/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void
CKEQYP(double* P, double* T, double* y, double* eqcon)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT}; /*temperature cache */
  double gort[11];      /* temporary storage */

  /*compute the Gibbs free energy */
  gibbs(gort, tc);

  /*compute the equilibrium constants */
  equilibriumConstants(eqcon, gort, tT);

  /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
  eqcon[0] *= 1e+06;

  /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
  eqcon[1] *= 1e-06;

  /*reaction 3: H2 + M <=> H + H + M */
  eqcon[2] *= 1e-06;

  /*reaction 4: O + O + M <=> O2 + M */
  eqcon[3] *= 1e+06;

  /*reaction 5: O + H + M <=> OH + M */
  eqcon[4] *= 1e+06;

  /*reaction 6: H2O + M <=> H + OH + M */
  eqcon[5] *= 1e-06;

  /*reaction 7: O + OH + M <=> HO2 + M */
  eqcon[6] *= 1e+06;

  /*reaction 8: H + O2 <=> O + OH */
  /*eqcon[7] *= 1;  */

  /*reaction 9: O + H2 <=> H + OH */
  /*eqcon[8] *= 1;  */

  /*reaction 10: O + H2 <=> H + OH */
  /*eqcon[9] *= 1;  */

  /*reaction 11: H2 + OH <=> H2O + H */
  /*eqcon[10] *= 1;  */

  /*reaction 12: OH + OH <=> O + H2O */
  /*eqcon[11] *= 1;  */

  /*reaction 13: H2 + AR <=> H + H + AR */
  eqcon[12] *= 1e-06;

  /*reaction 14: H2 + HE <=> H + H + HE */
  eqcon[13] *= 1e-06;

  /*reaction 15: O + O + AR <=> O2 + AR */
  eqcon[14] *= 1e+06;

  /*reaction 16: O + O + HE <=> O2 + HE */
  eqcon[15] *= 1e+06;

  /*reaction 17: H2O + H2O <=> H + OH + H2O */
  eqcon[16] *= 1e-06;

  /*reaction 18: HO2 + H <=> H2 + O2 */
  /*eqcon[17] *= 1;  */

  /*reaction 19: HO2 + H <=> OH + OH */
  /*eqcon[18] *= 1;  */

  /*reaction 20: HO2 + O <=> O2 + OH */
  /*eqcon[19] *= 1;  */

  /*reaction 21: HO2 + OH <=> H2O + O2 */
  /*eqcon[20] *= 1;  */

  /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
  /*eqcon[21] *= 1;  */

  /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
  /*eqcon[22] *= 1;  */

  /*reaction 24: H2O2 + H <=> H2O + OH */
  /*eqcon[23] *= 1;  */

  /*reaction 25: H2O2 + H <=> HO2 + H2 */
  /*eqcon[24] *= 1;  */

  /*reaction 26: H2O2 + O <=> OH + HO2 */
  /*eqcon[25] *= 1;  */

  /*reaction 27: H2O2 + OH <=> HO2 + H2O */
  /*eqcon[26] *= 1;  */

  /*reaction 28: H2O2 + OH <=> HO2 + H2O */
  /*eqcon[27] *= 1;  */

  /*reaction 29: HO2 + H <=> O + H2O */
  /*eqcon[28] *= 1;  */
}

/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void
CKEQXP(double* P, double* T, double* x, double* eqcon)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT}; /*temperature cache */
  double gort[11];      /* temporary storage */

  /*compute the Gibbs free energy */
  gibbs(gort, tc);

  /*compute the equilibrium constants */
  equilibriumConstants(eqcon, gort, tT);

  /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
  eqcon[0] *= 1e+06;

  /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
  eqcon[1] *= 1e-06;

  /*reaction 3: H2 + M <=> H + H + M */
  eqcon[2] *= 1e-06;

  /*reaction 4: O + O + M <=> O2 + M */
  eqcon[3] *= 1e+06;

  /*reaction 5: O + H + M <=> OH + M */
  eqcon[4] *= 1e+06;

  /*reaction 6: H2O + M <=> H + OH + M */
  eqcon[5] *= 1e-06;

  /*reaction 7: O + OH + M <=> HO2 + M */
  eqcon[6] *= 1e+06;

  /*reaction 8: H + O2 <=> O + OH */
  /*eqcon[7] *= 1;  */

  /*reaction 9: O + H2 <=> H + OH */
  /*eqcon[8] *= 1;  */

  /*reaction 10: O + H2 <=> H + OH */
  /*eqcon[9] *= 1;  */

  /*reaction 11: H2 + OH <=> H2O + H */
  /*eqcon[10] *= 1;  */

  /*reaction 12: OH + OH <=> O + H2O */
  /*eqcon[11] *= 1;  */

  /*reaction 13: H2 + AR <=> H + H + AR */
  eqcon[12] *= 1e-06;

  /*reaction 14: H2 + HE <=> H + H + HE */
  eqcon[13] *= 1e-06;

  /*reaction 15: O + O + AR <=> O2 + AR */
  eqcon[14] *= 1e+06;

  /*reaction 16: O + O + HE <=> O2 + HE */
  eqcon[15] *= 1e+06;

  /*reaction 17: H2O + H2O <=> H + OH + H2O */
  eqcon[16] *= 1e-06;

  /*reaction 18: HO2 + H <=> H2 + O2 */
  /*eqcon[17] *= 1;  */

  /*reaction 19: HO2 + H <=> OH + OH */
  /*eqcon[18] *= 1;  */

  /*reaction 20: HO2 + O <=> O2 + OH */
  /*eqcon[19] *= 1;  */

  /*reaction 21: HO2 + OH <=> H2O + O2 */
  /*eqcon[20] *= 1;  */

  /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
  /*eqcon[21] *= 1;  */

  /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
  /*eqcon[22] *= 1;  */

  /*reaction 24: H2O2 + H <=> H2O + OH */
  /*eqcon[23] *= 1;  */

  /*reaction 25: H2O2 + H <=> HO2 + H2 */
  /*eqcon[24] *= 1;  */

  /*reaction 26: H2O2 + O <=> OH + HO2 */
  /*eqcon[25] *= 1;  */

  /*reaction 27: H2O2 + OH <=> HO2 + H2O */
  /*eqcon[26] *= 1;  */

  /*reaction 28: H2O2 + OH <=> HO2 + H2O */
  /*eqcon[27] *= 1;  */

  /*reaction 29: HO2 + H <=> O + H2O */
  /*eqcon[28] *= 1;  */
}

/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void
CKEQYR(double* rho, double* T, double* y, double* eqcon)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT}; /*temperature cache */
  double gort[11];      /* temporary storage */

  /*compute the Gibbs free energy */
  gibbs(gort, tc);

  /*compute the equilibrium constants */
  equilibriumConstants(eqcon, gort, tT);

  /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
  eqcon[0] *= 1e+06;

  /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
  eqcon[1] *= 1e-06;

  /*reaction 3: H2 + M <=> H + H + M */
  eqcon[2] *= 1e-06;

  /*reaction 4: O + O + M <=> O2 + M */
  eqcon[3] *= 1e+06;

  /*reaction 5: O + H + M <=> OH + M */
  eqcon[4] *= 1e+06;

  /*reaction 6: H2O + M <=> H + OH + M */
  eqcon[5] *= 1e-06;

  /*reaction 7: O + OH + M <=> HO2 + M */
  eqcon[6] *= 1e+06;

  /*reaction 8: H + O2 <=> O + OH */
  /*eqcon[7] *= 1;  */

  /*reaction 9: O + H2 <=> H + OH */
  /*eqcon[8] *= 1;  */

  /*reaction 10: O + H2 <=> H + OH */
  /*eqcon[9] *= 1;  */

  /*reaction 11: H2 + OH <=> H2O + H */
  /*eqcon[10] *= 1;  */

  /*reaction 12: OH + OH <=> O + H2O */
  /*eqcon[11] *= 1;  */

  /*reaction 13: H2 + AR <=> H + H + AR */
  eqcon[12] *= 1e-06;

  /*reaction 14: H2 + HE <=> H + H + HE */
  eqcon[13] *= 1e-06;

  /*reaction 15: O + O + AR <=> O2 + AR */
  eqcon[14] *= 1e+06;

  /*reaction 16: O + O + HE <=> O2 + HE */
  eqcon[15] *= 1e+06;

  /*reaction 17: H2O + H2O <=> H + OH + H2O */
  eqcon[16] *= 1e-06;

  /*reaction 18: HO2 + H <=> H2 + O2 */
  /*eqcon[17] *= 1;  */

  /*reaction 19: HO2 + H <=> OH + OH */
  /*eqcon[18] *= 1;  */

  /*reaction 20: HO2 + O <=> O2 + OH */
  /*eqcon[19] *= 1;  */

  /*reaction 21: HO2 + OH <=> H2O + O2 */
  /*eqcon[20] *= 1;  */

  /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
  /*eqcon[21] *= 1;  */

  /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
  /*eqcon[22] *= 1;  */

  /*reaction 24: H2O2 + H <=> H2O + OH */
  /*eqcon[23] *= 1;  */

  /*reaction 25: H2O2 + H <=> HO2 + H2 */
  /*eqcon[24] *= 1;  */

  /*reaction 26: H2O2 + O <=> OH + HO2 */
  /*eqcon[25] *= 1;  */

  /*reaction 27: H2O2 + OH <=> HO2 + H2O */
  /*eqcon[26] *= 1;  */

  /*reaction 28: H2O2 + OH <=> HO2 + H2O */
  /*eqcon[27] *= 1;  */

  /*reaction 29: HO2 + H <=> O + H2O */
  /*eqcon[28] *= 1;  */
}

/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void
CKEQXR(double* rho, double* T, double* x, double* eqcon)
{
  double tT = *T; /*temporary temperature */
  double tc[] = {
    log(tT), tT, tT * tT, tT * tT * tT,
    tT * tT * tT * tT}; /*temperature cache */
  double gort[11];      /* temporary storage */

  /*compute the Gibbs free energy */
  gibbs(gort, tc);

  /*compute the equilibrium constants */
  equilibriumConstants(eqcon, gort, tT);

  /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
  eqcon[0] *= 1e+06;

  /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
  eqcon[1] *= 1e-06;

  /*reaction 3: H2 + M <=> H + H + M */
  eqcon[2] *= 1e-06;

  /*reaction 4: O + O + M <=> O2 + M */
  eqcon[3] *= 1e+06;

  /*reaction 5: O + H + M <=> OH + M */
  eqcon[4] *= 1e+06;

  /*reaction 6: H2O + M <=> H + OH + M */
  eqcon[5] *= 1e-06;

  /*reaction 7: O + OH + M <=> HO2 + M */
  eqcon[6] *= 1e+06;

  /*reaction 8: H + O2 <=> O + OH */
  /*eqcon[7] *= 1;  */

  /*reaction 9: O + H2 <=> H + OH */
  /*eqcon[8] *= 1;  */

  /*reaction 10: O + H2 <=> H + OH */
  /*eqcon[9] *= 1;  */

  /*reaction 11: H2 + OH <=> H2O + H */
  /*eqcon[10] *= 1;  */

  /*reaction 12: OH + OH <=> O + H2O */
  /*eqcon[11] *= 1;  */

  /*reaction 13: H2 + AR <=> H + H + AR */
  eqcon[12] *= 1e-06;

  /*reaction 14: H2 + HE <=> H + H + HE */
  eqcon[13] *= 1e-06;

  /*reaction 15: O + O + AR <=> O2 + AR */
  eqcon[14] *= 1e+06;

  /*reaction 16: O + O + HE <=> O2 + HE */
  eqcon[15] *= 1e+06;

  /*reaction 17: H2O + H2O <=> H + OH + H2O */
  eqcon[16] *= 1e-06;

  /*reaction 18: HO2 + H <=> H2 + O2 */
  /*eqcon[17] *= 1;  */

  /*reaction 19: HO2 + H <=> OH + OH */
  /*eqcon[18] *= 1;  */

  /*reaction 20: HO2 + O <=> O2 + OH */
  /*eqcon[19] *= 1;  */

  /*reaction 21: HO2 + OH <=> H2O + O2 */
  /*eqcon[20] *= 1;  */

  /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
  /*eqcon[21] *= 1;  */

  /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
  /*eqcon[22] *= 1;  */

  /*reaction 24: H2O2 + H <=> H2O + OH */
  /*eqcon[23] *= 1;  */

  /*reaction 25: H2O2 + H <=> HO2 + H2 */
  /*eqcon[24] *= 1;  */

  /*reaction 26: H2O2 + O <=> OH + HO2 */
  /*eqcon[25] *= 1;  */

  /*reaction 27: H2O2 + OH <=> HO2 + H2O */
  /*eqcon[26] *= 1;  */

  /*reaction 28: H2O2 + OH <=> HO2 + H2O */
  /*eqcon[27] *= 1;  */

  /*reaction 29: HO2 + H <=> O + H2O */
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
void
productionRate(double* wdot, double* sc, double T)
{
  double tc[] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; /*temperature cache */
  double invT = 1.0 / tc[1];

  if (T != T_save) {
    T_save = T;
    comp_k_f(tc, invT, k_f_save);
    comp_Kc(tc, invT, Kc_save);
  }

  double qdot, q_f[29], q_r[29];
  comp_qfqr(q_f, q_r, sc, tc, invT);

  for (int i = 0; i < 11; ++i) {
    wdot[i] = 0.0;
  }

  qdot = q_f[0] - q_r[0];
  wdot[0] -= qdot;
  wdot[5] -= qdot;
  wdot[6] += qdot;

  qdot = q_f[1] - q_r[1];
  wdot[3] += qdot;
  wdot[3] += qdot;
  wdot[7] -= qdot;

  qdot = q_f[2] - q_r[2];
  wdot[0] += qdot;
  wdot[0] += qdot;
  wdot[1] -= qdot;

  qdot = q_f[3] - q_r[3];
  wdot[2] -= qdot;
  wdot[2] -= qdot;
  wdot[5] += qdot;

  qdot = q_f[4] - q_r[4];
  wdot[0] -= qdot;
  wdot[2] -= qdot;
  wdot[3] += qdot;

  qdot = q_f[5] - q_r[5];
  wdot[0] += qdot;
  wdot[3] += qdot;
  wdot[4] -= qdot;

  qdot = q_f[6] - q_r[6];
  wdot[2] -= qdot;
  wdot[3] -= qdot;
  wdot[6] += qdot;

  qdot = q_f[7] - q_r[7];
  wdot[0] -= qdot;
  wdot[2] += qdot;
  wdot[3] += qdot;
  wdot[5] -= qdot;

  qdot = q_f[8] - q_r[8];
  wdot[0] += qdot;
  wdot[1] -= qdot;
  wdot[2] -= qdot;
  wdot[3] += qdot;

  qdot = q_f[9] - q_r[9];
  wdot[0] += qdot;
  wdot[1] -= qdot;
  wdot[2] -= qdot;
  wdot[3] += qdot;

  qdot = q_f[10] - q_r[10];
  wdot[0] += qdot;
  wdot[1] -= qdot;
  wdot[3] -= qdot;
  wdot[4] += qdot;

  qdot = q_f[11] - q_r[11];
  wdot[2] += qdot;
  wdot[3] -= qdot;
  wdot[3] -= qdot;
  wdot[4] += qdot;

  qdot = q_f[12] - q_r[12];
  wdot[0] += qdot;
  wdot[0] += qdot;
  wdot[1] -= qdot;
  wdot[9] -= qdot;
  wdot[9] += qdot;

  qdot = q_f[13] - q_r[13];
  wdot[0] += qdot;
  wdot[0] += qdot;
  wdot[1] -= qdot;
  wdot[10] -= qdot;
  wdot[10] += qdot;

  qdot = q_f[14] - q_r[14];
  wdot[2] -= qdot;
  wdot[2] -= qdot;
  wdot[5] += qdot;
  wdot[9] -= qdot;
  wdot[9] += qdot;

  qdot = q_f[15] - q_r[15];
  wdot[2] -= qdot;
  wdot[2] -= qdot;
  wdot[5] += qdot;
  wdot[10] -= qdot;
  wdot[10] += qdot;

  qdot = q_f[16] - q_r[16];
  wdot[0] += qdot;
  wdot[3] += qdot;
  wdot[4] -= qdot;
  wdot[4] -= qdot;
  wdot[4] += qdot;

  qdot = q_f[17] - q_r[17];
  wdot[0] -= qdot;
  wdot[1] += qdot;
  wdot[5] += qdot;
  wdot[6] -= qdot;

  qdot = q_f[18] - q_r[18];
  wdot[0] -= qdot;
  wdot[3] += qdot;
  wdot[3] += qdot;
  wdot[6] -= qdot;

  qdot = q_f[19] - q_r[19];
  wdot[2] -= qdot;
  wdot[3] += qdot;
  wdot[5] += qdot;
  wdot[6] -= qdot;

  qdot = q_f[20] - q_r[20];
  wdot[3] -= qdot;
  wdot[4] += qdot;
  wdot[5] += qdot;
  wdot[6] -= qdot;

  qdot = q_f[21] - q_r[21];
  wdot[5] += qdot;
  wdot[6] -= qdot;
  wdot[6] -= qdot;
  wdot[7] += qdot;

  qdot = q_f[22] - q_r[22];
  wdot[5] += qdot;
  wdot[6] -= qdot;
  wdot[6] -= qdot;
  wdot[7] += qdot;

  qdot = q_f[23] - q_r[23];
  wdot[0] -= qdot;
  wdot[3] += qdot;
  wdot[4] += qdot;
  wdot[7] -= qdot;

  qdot = q_f[24] - q_r[24];
  wdot[0] -= qdot;
  wdot[1] += qdot;
  wdot[6] += qdot;
  wdot[7] -= qdot;

  qdot = q_f[25] - q_r[25];
  wdot[2] -= qdot;
  wdot[3] += qdot;
  wdot[6] += qdot;
  wdot[7] -= qdot;

  qdot = q_f[26] - q_r[26];
  wdot[3] -= qdot;
  wdot[4] += qdot;
  wdot[6] += qdot;
  wdot[7] -= qdot;

  qdot = q_f[27] - q_r[27];
  wdot[3] -= qdot;
  wdot[4] += qdot;
  wdot[6] += qdot;
  wdot[7] -= qdot;

  qdot = q_f[28] - q_r[28];
  wdot[0] -= qdot;
  wdot[2] += qdot;
  wdot[4] += qdot;
  wdot[6] -= qdot;

  return;
}

void
comp_k_f(double* tc, double invT, double* k_f)
{
#ifdef __INTEL_COMPILER
#pragma simd
#endif
  for (int i = 0; i < 29; ++i) {
    k_f[i] = prefactor_units[i] * fwd_A[i] *
             exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
  };
  return;
}

void
comp_Kc(double* tc, double invT, double* Kc)
{
  /*compute the Gibbs free energy */
  double g_RT[11];
  gibbs(g_RT, tc);

  Kc[0] = g_RT[0] + g_RT[5] - g_RT[6];
  Kc[1] = -g_RT[3] - g_RT[3] + g_RT[7];
  Kc[2] = -g_RT[0] - g_RT[0] + g_RT[1];
  Kc[3] = g_RT[2] + g_RT[2] - g_RT[5];
  Kc[4] = g_RT[0] + g_RT[2] - g_RT[3];
  Kc[5] = -g_RT[0] - g_RT[3] + g_RT[4];
  Kc[6] = g_RT[2] + g_RT[3] - g_RT[6];
  Kc[7] = g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5];
  Kc[8] = -g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3];
  Kc[9] = -g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3];
  Kc[10] = -g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4];
  Kc[11] = -g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4];
  Kc[12] = -g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9];
  Kc[13] = -g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10];
  Kc[14] = g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9];
  Kc[15] = g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10];
  Kc[16] = -g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4];
  Kc[17] = g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6];
  Kc[18] = g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6];
  Kc[19] = g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6];
  Kc[20] = g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6];
  Kc[21] = -g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7];
  Kc[22] = -g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7];
  Kc[23] = g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7];
  Kc[24] = g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7];
  Kc[25] = g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7];
  Kc[26] = g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7];
  Kc[27] = g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7];
  Kc[28] = g_RT[0] - g_RT[2] - g_RT[4] + g_RT[6];

#ifdef __INTEL_COMPILER
#pragma simd
#endif
  for (int i = 0; i < 29; ++i) {
    Kc[i] = exp(Kc[i]);
  };

  /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
  double refC = 101325 / 8.31451 * invT;
  double refCinv = 1 / refC;

  Kc[0] *= refCinv;
  Kc[1] *= refC;
  Kc[2] *= refC;
  Kc[3] *= refCinv;
  Kc[4] *= refCinv;
  Kc[5] *= refC;
  Kc[6] *= refCinv;
  Kc[12] *= refC;
  Kc[13] *= refC;
  Kc[14] *= refCinv;
  Kc[15] *= refCinv;
  Kc[16] *= refC;

  return;
}

void
comp_qfqr(double* qf, double* qr, double* sc, double* tc, double invT)
{

  /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
  qf[0] = sc[0] * sc[5];
  qr[0] = sc[6];

  /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
  qf[1] = sc[7];
  qr[1] = sc[3] * sc[3];

  /*reaction 3: H2 + M <=> H + H + M */
  qf[2] = sc[1];
  qr[2] = sc[0] * sc[0];

  /*reaction 4: O + O + M <=> O2 + M */
  qf[3] = sc[2] * sc[2];
  qr[3] = sc[5];

  /*reaction 5: O + H + M <=> OH + M */
  qf[4] = sc[0] * sc[2];
  qr[4] = sc[3];

  /*reaction 6: H2O + M <=> H + OH + M */
  qf[5] = sc[4];
  qr[5] = sc[0] * sc[3];

  /*reaction 7: O + OH + M <=> HO2 + M */
  qf[6] = sc[2] * sc[3];
  qr[6] = sc[6];

  /*reaction 8: H + O2 <=> O + OH */
  qf[7] = sc[0] * sc[5];
  qr[7] = sc[2] * sc[3];

  /*reaction 9: O + H2 <=> H + OH */
  qf[8] = sc[1] * sc[2];
  qr[8] = sc[0] * sc[3];

  /*reaction 10: O + H2 <=> H + OH */
  qf[9] = sc[1] * sc[2];
  qr[9] = sc[0] * sc[3];

  /*reaction 11: H2 + OH <=> H2O + H */
  qf[10] = sc[1] * sc[3];
  qr[10] = sc[0] * sc[4];

  /*reaction 12: OH + OH <=> O + H2O */
  qf[11] = sc[3] * sc[3];
  qr[11] = sc[2] * sc[4];

  /*reaction 13: H2 + AR <=> H + H + AR */
  qf[12] = sc[1] * sc[9];
  qr[12] = sc[0] * sc[0] * sc[9];

  /*reaction 14: H2 + HE <=> H + H + HE */
  qf[13] = sc[1] * sc[10];
  qr[13] = sc[0] * sc[0] * sc[10];

  /*reaction 15: O + O + AR <=> O2 + AR */
  qf[14] = sc[2] * sc[2] * sc[9];
  qr[14] = sc[5] * sc[9];

  /*reaction 16: O + O + HE <=> O2 + HE */
  qf[15] = sc[2] * sc[2] * sc[10];
  qr[15] = sc[5] * sc[10];

  /*reaction 17: H2O + H2O <=> H + OH + H2O */
  qf[16] = sc[4] * sc[4];
  qr[16] = sc[0] * sc[3] * sc[4];

  /*reaction 18: HO2 + H <=> H2 + O2 */
  qf[17] = sc[0] * sc[6];
  qr[17] = sc[1] * sc[5];

  /*reaction 19: HO2 + H <=> OH + OH */
  qf[18] = sc[0] * sc[6];
  qr[18] = sc[3] * sc[3];

  /*reaction 20: HO2 + O <=> O2 + OH */
  qf[19] = sc[2] * sc[6];
  qr[19] = sc[3] * sc[5];

  /*reaction 21: HO2 + OH <=> H2O + O2 */
  qf[20] = sc[3] * sc[6];
  qr[20] = sc[4] * sc[5];

  /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
  qf[21] = sc[6] * sc[6];
  qr[21] = sc[5] * sc[7];

  /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
  qf[22] = sc[6] * sc[6];
  qr[22] = sc[5] * sc[7];

  /*reaction 24: H2O2 + H <=> H2O + OH */
  qf[23] = sc[0] * sc[7];
  qr[23] = sc[3] * sc[4];

  /*reaction 25: H2O2 + H <=> HO2 + H2 */
  qf[24] = sc[0] * sc[7];
  qr[24] = sc[1] * sc[6];

  /*reaction 26: H2O2 + O <=> OH + HO2 */
  qf[25] = sc[2] * sc[7];
  qr[25] = sc[3] * sc[6];

  /*reaction 27: H2O2 + OH <=> HO2 + H2O */
  qf[26] = sc[3] * sc[7];
  qr[26] = sc[4] * sc[6];

  /*reaction 28: H2O2 + OH <=> HO2 + H2O */
  qf[27] = sc[3] * sc[7];
  qr[27] = sc[4] * sc[6];

  /*reaction 29: HO2 + H <=> O + H2O */
  qf[28] = sc[0] * sc[6];
  qr[28] = sc[2] * sc[4];

  double T = tc[1];

  /*compute the mixture concentration */
  double mixture = 0.0;
  for (int i = 0; i < 11; ++i) {
    mixture += sc[i];
  }

  double Corr[29];
  for (int i = 0; i < 29; ++i) {
    Corr[i] = 1.0;
  }

  /* troe */
  {
    double alpha[2];
    alpha[0] = mixture + (TB[0][0] - 1) * sc[1] + (TB[0][1] - 1) * sc[4] +
               (TB[0][2] - 1) * sc[5] + (TB[0][3] - 1) * sc[9] +
               (TB[0][4] - 1) * sc[10];
    alpha[1] = mixture + (TB[1][0] - 1) * sc[4] + (TB[1][1] - 1) * sc[8] +
               (TB[1][2] - 1) * sc[5] + (TB[1][3] - 1) * sc[10] +
               (TB[1][4] - 1) * sc[7] + (TB[1][5] - 1) * sc[1];
    for (int i = 0; i < 2; i++) {
      double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;
      redP = alpha[i - 0] / k_f_save[i] * phase_units[i] * low_A[i] *
             exp(low_beta[i] * tc[0] - activation_units[i] * low_Ea[i] * invT);
      F = redP / (1.0 + redP);
      logPred = log10(redP);
      logFcent = log10(
        (fabs(troe_Tsss[i]) > 1.e-100
           ? (1. - troe_a[i]) * exp(-T / troe_Tsss[i])
           : 0.) +
        (fabs(troe_Ts[i]) > 1.e-100 ? troe_a[i] * exp(-T / troe_Ts[i]) : 0.) +
        (troe_len[i] == 4 ? exp(-troe_Tss[i] * invT) : 0.));
      troe_c = -.4 - .67 * logFcent;
      troe_n = .75 - 1.27 * logFcent;
      troe = (troe_c + logPred) / (troe_n - .14 * (troe_c + logPred));
      F_troe = pow(10., logFcent / (1.0 + troe * troe));
      Corr[i] = F * F_troe;
    }
  }

  /* simple three-body correction */
  {
    double alpha;
    alpha = mixture + (TB[2][0] - 1) * sc[1] + (TB[2][1] - 1) * sc[4] +
            (TB[2][2] - 1) * sc[9] + (TB[2][3] - 1) * sc[10];
    Corr[2] = alpha;
    alpha = mixture + (TB[3][0] - 1) * sc[1] + (TB[3][1] - 1) * sc[4] +
            (TB[3][2] - 1) * sc[9] + (TB[3][3] - 1) * sc[10];
    Corr[3] = alpha;
    alpha = mixture + (TB[4][0] - 1) * sc[1] + (TB[4][1] - 1) * sc[4] +
            (TB[4][2] - 1) * sc[9] + (TB[4][3] - 1) * sc[10];
    Corr[4] = alpha;
    alpha = mixture + (TB[5][0] - 1) * sc[1] + (TB[5][1] - 1) * sc[4] +
            (TB[5][2] - 1) * sc[10] + (TB[5][3] - 1) * sc[8] +
            (TB[5][4] - 1) * sc[5];
    Corr[5] = alpha;
    alpha = mixture + (TB[6][0] - 1) * sc[1] + (TB[6][1] - 1) * sc[4] +
            (TB[6][2] - 1) * sc[9] + (TB[6][3] - 1) * sc[10];
    Corr[6] = alpha;
  }

  for (int i = 0; i < 29; i++) {
    qf[i] *= Corr[i] * k_f_save[i];
    qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
  }

  return;
}

/*compute the production rate for each species */
void
vproductionRate(int npt, double* wdot, double* sc, double* T)
{
  double k_f_s[29 * npt], Kc_s[29 * npt], mixture[npt], g_RT[11 * npt];
  double tc[5 * npt], invT[npt];

#ifdef __INTEL_COMPILER
#pragma simd
#endif
  for (int i = 0; i < npt; i++) {
    tc[0 * npt + i] = log(T[i]);
    tc[1 * npt + i] = T[i];
    tc[2 * npt + i] = T[i] * T[i];
    tc[3 * npt + i] = T[i] * T[i] * T[i];
    tc[4 * npt + i] = T[i] * T[i] * T[i] * T[i];
    invT[i] = 1.0 / T[i];
  }

  for (int i = 0; i < npt; i++) {
    mixture[i] = 0.0;
  }

  for (int n = 0; n < 11; n++) {
    for (int i = 0; i < npt; i++) {
      mixture[i] += sc[n * npt + i];
      wdot[n * npt + i] = 0.0;
    }
  }

  vcomp_k_f(npt, k_f_s, tc, invT);

  vcomp_gibbs(npt, g_RT, tc);

  vcomp_Kc(npt, Kc_s, g_RT, invT);

  vcomp_wdot(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
}

void
vcomp_k_f(int npt, double* k_f_s, double* tc, double* invT)
{
#ifdef __INTEL_COMPILER
#pragma simd
#endif
  for (int i = 0; i < npt; i++) {
    k_f_s[0 * npt + i] =
      prefactor_units[0] * fwd_A[0] *
      exp(fwd_beta[0] * tc[i] - activation_units[0] * fwd_Ea[0] * invT[i]);
    k_f_s[1 * npt + i] =
      prefactor_units[1] * fwd_A[1] *
      exp(fwd_beta[1] * tc[i] - activation_units[1] * fwd_Ea[1] * invT[i]);
    k_f_s[2 * npt + i] =
      prefactor_units[2] * fwd_A[2] *
      exp(fwd_beta[2] * tc[i] - activation_units[2] * fwd_Ea[2] * invT[i]);
    k_f_s[3 * npt + i] =
      prefactor_units[3] * fwd_A[3] *
      exp(fwd_beta[3] * tc[i] - activation_units[3] * fwd_Ea[3] * invT[i]);
    k_f_s[4 * npt + i] =
      prefactor_units[4] * fwd_A[4] *
      exp(fwd_beta[4] * tc[i] - activation_units[4] * fwd_Ea[4] * invT[i]);
    k_f_s[5 * npt + i] =
      prefactor_units[5] * fwd_A[5] *
      exp(fwd_beta[5] * tc[i] - activation_units[5] * fwd_Ea[5] * invT[i]);
    k_f_s[6 * npt + i] =
      prefactor_units[6] * fwd_A[6] *
      exp(fwd_beta[6] * tc[i] - activation_units[6] * fwd_Ea[6] * invT[i]);
    k_f_s[7 * npt + i] =
      prefactor_units[7] * fwd_A[7] *
      exp(fwd_beta[7] * tc[i] - activation_units[7] * fwd_Ea[7] * invT[i]);
    k_f_s[8 * npt + i] =
      prefactor_units[8] * fwd_A[8] *
      exp(fwd_beta[8] * tc[i] - activation_units[8] * fwd_Ea[8] * invT[i]);
    k_f_s[9 * npt + i] =
      prefactor_units[9] * fwd_A[9] *
      exp(fwd_beta[9] * tc[i] - activation_units[9] * fwd_Ea[9] * invT[i]);
    k_f_s[10 * npt + i] =
      prefactor_units[10] * fwd_A[10] *
      exp(fwd_beta[10] * tc[i] - activation_units[10] * fwd_Ea[10] * invT[i]);
    k_f_s[11 * npt + i] =
      prefactor_units[11] * fwd_A[11] *
      exp(fwd_beta[11] * tc[i] - activation_units[11] * fwd_Ea[11] * invT[i]);
    k_f_s[12 * npt + i] =
      prefactor_units[12] * fwd_A[12] *
      exp(fwd_beta[12] * tc[i] - activation_units[12] * fwd_Ea[12] * invT[i]);
    k_f_s[13 * npt + i] =
      prefactor_units[13] * fwd_A[13] *
      exp(fwd_beta[13] * tc[i] - activation_units[13] * fwd_Ea[13] * invT[i]);
    k_f_s[14 * npt + i] =
      prefactor_units[14] * fwd_A[14] *
      exp(fwd_beta[14] * tc[i] - activation_units[14] * fwd_Ea[14] * invT[i]);
    k_f_s[15 * npt + i] =
      prefactor_units[15] * fwd_A[15] *
      exp(fwd_beta[15] * tc[i] - activation_units[15] * fwd_Ea[15] * invT[i]);
    k_f_s[16 * npt + i] =
      prefactor_units[16] * fwd_A[16] *
      exp(fwd_beta[16] * tc[i] - activation_units[16] * fwd_Ea[16] * invT[i]);
    k_f_s[17 * npt + i] =
      prefactor_units[17] * fwd_A[17] *
      exp(fwd_beta[17] * tc[i] - activation_units[17] * fwd_Ea[17] * invT[i]);
    k_f_s[18 * npt + i] =
      prefactor_units[18] * fwd_A[18] *
      exp(fwd_beta[18] * tc[i] - activation_units[18] * fwd_Ea[18] * invT[i]);
    k_f_s[19 * npt + i] =
      prefactor_units[19] * fwd_A[19] *
      exp(fwd_beta[19] * tc[i] - activation_units[19] * fwd_Ea[19] * invT[i]);
    k_f_s[20 * npt + i] =
      prefactor_units[20] * fwd_A[20] *
      exp(fwd_beta[20] * tc[i] - activation_units[20] * fwd_Ea[20] * invT[i]);
    k_f_s[21 * npt + i] =
      prefactor_units[21] * fwd_A[21] *
      exp(fwd_beta[21] * tc[i] - activation_units[21] * fwd_Ea[21] * invT[i]);
    k_f_s[22 * npt + i] =
      prefactor_units[22] * fwd_A[22] *
      exp(fwd_beta[22] * tc[i] - activation_units[22] * fwd_Ea[22] * invT[i]);
    k_f_s[23 * npt + i] =
      prefactor_units[23] * fwd_A[23] *
      exp(fwd_beta[23] * tc[i] - activation_units[23] * fwd_Ea[23] * invT[i]);
    k_f_s[24 * npt + i] =
      prefactor_units[24] * fwd_A[24] *
      exp(fwd_beta[24] * tc[i] - activation_units[24] * fwd_Ea[24] * invT[i]);
    k_f_s[25 * npt + i] =
      prefactor_units[25] * fwd_A[25] *
      exp(fwd_beta[25] * tc[i] - activation_units[25] * fwd_Ea[25] * invT[i]);
    k_f_s[26 * npt + i] =
      prefactor_units[26] * fwd_A[26] *
      exp(fwd_beta[26] * tc[i] - activation_units[26] * fwd_Ea[26] * invT[i]);
    k_f_s[27 * npt + i] =
      prefactor_units[27] * fwd_A[27] *
      exp(fwd_beta[27] * tc[i] - activation_units[27] * fwd_Ea[27] * invT[i]);
    k_f_s[28 * npt + i] =
      prefactor_units[28] * fwd_A[28] *
      exp(fwd_beta[28] * tc[i] - activation_units[28] * fwd_Ea[28] * invT[i]);
  }
}

void
vcomp_gibbs(int npt, double* g_RT, double* tc)
{
  /*compute the Gibbs free energy */
  for (int i = 0; i < npt; i++) {
    double tg[5], g[11];
    tg[0] = tc[0 * npt + i];
    tg[1] = tc[1 * npt + i];
    tg[2] = tc[2 * npt + i];
    tg[3] = tc[3 * npt + i];
    tg[4] = tc[4 * npt + i];

    gibbs(g, tg);

    g_RT[0 * npt + i] = g[0];
    g_RT[1 * npt + i] = g[1];
    g_RT[2 * npt + i] = g[2];
    g_RT[3 * npt + i] = g[3];
    g_RT[4 * npt + i] = g[4];
    g_RT[5 * npt + i] = g[5];
    g_RT[6 * npt + i] = g[6];
    g_RT[7 * npt + i] = g[7];
    g_RT[8 * npt + i] = g[8];
    g_RT[9 * npt + i] = g[9];
    g_RT[10 * npt + i] = g[10];
  }
}

void
vcomp_Kc(int npt, double* Kc_s, double* g_RT, double* invT)
{
#ifdef __INTEL_COMPILER
#pragma simd
#endif
  for (int i = 0; i < npt; i++) {
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = (101325. / 8.31451) * invT[i];
    double refCinv = 1.0 / refC;

    Kc_s[0 * npt + i] =
      refCinv *
      exp((g_RT[0 * npt + i] + g_RT[5 * npt + i]) - (g_RT[6 * npt + i]));
    Kc_s[1 * npt + i] =
      refC * exp((g_RT[7 * npt + i]) - (g_RT[3 * npt + i] + g_RT[3 * npt + i]));
    Kc_s[2 * npt + i] =
      refC * exp((g_RT[1 * npt + i]) - (g_RT[0 * npt + i] + g_RT[0 * npt + i]));
    Kc_s[3 * npt + i] =
      refCinv *
      exp((g_RT[2 * npt + i] + g_RT[2 * npt + i]) - (g_RT[5 * npt + i]));
    Kc_s[4 * npt + i] =
      refCinv *
      exp((g_RT[0 * npt + i] + g_RT[2 * npt + i]) - (g_RT[3 * npt + i]));
    Kc_s[5 * npt + i] =
      refC * exp((g_RT[4 * npt + i]) - (g_RT[0 * npt + i] + g_RT[3 * npt + i]));
    Kc_s[6 * npt + i] =
      refCinv *
      exp((g_RT[2 * npt + i] + g_RT[3 * npt + i]) - (g_RT[6 * npt + i]));
    Kc_s[7 * npt + i] = exp(
      (g_RT[0 * npt + i] + g_RT[5 * npt + i]) -
      (g_RT[2 * npt + i] + g_RT[3 * npt + i]));
    Kc_s[8 * npt + i] = exp(
      (g_RT[1 * npt + i] + g_RT[2 * npt + i]) -
      (g_RT[0 * npt + i] + g_RT[3 * npt + i]));
    Kc_s[9 * npt + i] = exp(
      (g_RT[1 * npt + i] + g_RT[2 * npt + i]) -
      (g_RT[0 * npt + i] + g_RT[3 * npt + i]));
    Kc_s[10 * npt + i] = exp(
      (g_RT[1 * npt + i] + g_RT[3 * npt + i]) -
      (g_RT[0 * npt + i] + g_RT[4 * npt + i]));
    Kc_s[11 * npt + i] = exp(
      (g_RT[3 * npt + i] + g_RT[3 * npt + i]) -
      (g_RT[2 * npt + i] + g_RT[4 * npt + i]));
    Kc_s[12 * npt + i] =
      refC * exp(
               (g_RT[1 * npt + i] + g_RT[9 * npt + i]) -
               (g_RT[0 * npt + i] + g_RT[0 * npt + i] + g_RT[9 * npt + i]));
    Kc_s[13 * npt + i] =
      refC * exp(
               (g_RT[1 * npt + i] + g_RT[10 * npt + i]) -
               (g_RT[0 * npt + i] + g_RT[0 * npt + i] + g_RT[10 * npt + i]));
    Kc_s[14 * npt + i] =
      refCinv * exp(
                  (g_RT[2 * npt + i] + g_RT[2 * npt + i] + g_RT[9 * npt + i]) -
                  (g_RT[5 * npt + i] + g_RT[9 * npt + i]));
    Kc_s[15 * npt + i] =
      refCinv * exp(
                  (g_RT[2 * npt + i] + g_RT[2 * npt + i] + g_RT[10 * npt + i]) -
                  (g_RT[5 * npt + i] + g_RT[10 * npt + i]));
    Kc_s[16 * npt + i] =
      refC * exp(
               (g_RT[4 * npt + i] + g_RT[4 * npt + i]) -
               (g_RT[0 * npt + i] + g_RT[3 * npt + i] + g_RT[4 * npt + i]));
    Kc_s[17 * npt + i] = exp(
      (g_RT[0 * npt + i] + g_RT[6 * npt + i]) -
      (g_RT[1 * npt + i] + g_RT[5 * npt + i]));
    Kc_s[18 * npt + i] = exp(
      (g_RT[0 * npt + i] + g_RT[6 * npt + i]) -
      (g_RT[3 * npt + i] + g_RT[3 * npt + i]));
    Kc_s[19 * npt + i] = exp(
      (g_RT[2 * npt + i] + g_RT[6 * npt + i]) -
      (g_RT[3 * npt + i] + g_RT[5 * npt + i]));
    Kc_s[20 * npt + i] = exp(
      (g_RT[3 * npt + i] + g_RT[6 * npt + i]) -
      (g_RT[4 * npt + i] + g_RT[5 * npt + i]));
    Kc_s[21 * npt + i] = exp(
      (g_RT[6 * npt + i] + g_RT[6 * npt + i]) -
      (g_RT[5 * npt + i] + g_RT[7 * npt + i]));
    Kc_s[22 * npt + i] = exp(
      (g_RT[6 * npt + i] + g_RT[6 * npt + i]) -
      (g_RT[5 * npt + i] + g_RT[7 * npt + i]));
    Kc_s[23 * npt + i] = exp(
      (g_RT[0 * npt + i] + g_RT[7 * npt + i]) -
      (g_RT[3 * npt + i] + g_RT[4 * npt + i]));
    Kc_s[24 * npt + i] = exp(
      (g_RT[0 * npt + i] + g_RT[7 * npt + i]) -
      (g_RT[1 * npt + i] + g_RT[6 * npt + i]));
    Kc_s[25 * npt + i] = exp(
      (g_RT[2 * npt + i] + g_RT[7 * npt + i]) -
      (g_RT[3 * npt + i] + g_RT[6 * npt + i]));
    Kc_s[26 * npt + i] = exp(
      (g_RT[3 * npt + i] + g_RT[7 * npt + i]) -
      (g_RT[4 * npt + i] + g_RT[6 * npt + i]));
    Kc_s[27 * npt + i] = exp(
      (g_RT[3 * npt + i] + g_RT[7 * npt + i]) -
      (g_RT[4 * npt + i] + g_RT[6 * npt + i]));
    Kc_s[28 * npt + i] = exp(
      (g_RT[0 * npt + i] + g_RT[6 * npt + i]) -
      (g_RT[2 * npt + i] + g_RT[4 * npt + i]));
  }
}

void
vcomp_wdot(
  int npt,
  double* wdot,
  double* mixture,
  double* sc,
  double* k_f_s,
  double* Kc_s,
  double* tc,
  double* invT,
  double* T)
{
#ifdef __INTEL_COMPILER
#pragma simd
#endif
  for (int i = 0; i < npt; i++) {
    double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
    double alpha;
    double redP, F;
    double logPred;
    double logFcent, troe_c, troe_n, troe, F_troe;

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[0 * npt + i] * sc[5 * npt + i];
    alpha =
      mixture[i] + (TB[0][0] - 1) * sc[1 * npt + i] +
      (TB[0][1] - 1) * sc[4 * npt + i] + (TB[0][2] - 1) * sc[5 * npt + i] +
      (TB[0][3] - 1) * sc[9 * npt + i] + (TB[0][4] - 1) * sc[10 * npt + i];
    k_f = k_f_s[0 * npt + i];
    redP = alpha / k_f * phase_units[0] * low_A[0] *
           exp(low_beta[0] * tc[i] - activation_units[0] * low_Ea[0] * invT[i]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10(
      (fabs(troe_Tsss[0]) > 1.e-100
         ? (1. - troe_a[0]) * exp(-T[i] / troe_Tsss[0])
         : 0.) +
      (fabs(troe_Ts[0]) > 1.e-100 ? troe_a[0] * exp(-T[i] / troe_Ts[0]) : 0.) +
      (troe_len[0] == 4 ? exp(-troe_Tss[0] * invT[i]) : 0.));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14 * (troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe * troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[6 * npt + i];
    Kc = Kc_s[0 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] -= qdot;
    wdot[5 * npt + i] -= qdot;
    wdot[6 * npt + i] += qdot;

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    phi_f = sc[7 * npt + i];
    alpha = mixture[i] + (TB[1][0] - 1) * sc[4 * npt + i] +
            (TB[1][1] - 1) * sc[8 * npt + i] +
            (TB[1][2] - 1) * sc[5 * npt + i] +
            (TB[1][3] - 1) * sc[10 * npt + i] +
            (TB[1][4] - 1) * sc[7 * npt + i] + (TB[1][5] - 1) * sc[1 * npt + i];
    k_f = k_f_s[1 * npt + i];
    redP = alpha / k_f * phase_units[1] * low_A[1] *
           exp(low_beta[1] * tc[i] - activation_units[1] * low_Ea[1] * invT[i]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10(
      (fabs(troe_Tsss[1]) > 1.e-100
         ? (1. - troe_a[1]) * exp(-T[i] / troe_Tsss[1])
         : 0.) +
      (fabs(troe_Ts[1]) > 1.e-100 ? troe_a[1] * exp(-T[i] / troe_Ts[1]) : 0.) +
      (troe_len[1] == 4 ? exp(-troe_Tss[1] * invT[i]) : 0.));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14 * (troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe * troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[3 * npt + i] * sc[3 * npt + i];
    Kc = Kc_s[1 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3 * npt + i] += qdot;
    wdot[3 * npt + i] += qdot;
    wdot[7 * npt + i] -= qdot;

    /*reaction 3: H2 + M <=> H + H + M */
    phi_f = sc[1 * npt + i];
    alpha = mixture[i] + (TB[2][0] - 1) * sc[1 * npt + i] +
            (TB[2][1] - 1) * sc[4 * npt + i] +
            (TB[2][2] - 1) * sc[9 * npt + i] +
            (TB[2][3] - 1) * sc[10 * npt + i];
    k_f = alpha * k_f_s[2 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[0 * npt + i] * sc[0 * npt + i];
    Kc = Kc_s[2 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] += qdot;
    wdot[0 * npt + i] += qdot;
    wdot[1 * npt + i] -= qdot;

    /*reaction 4: O + O + M <=> O2 + M */
    phi_f = sc[2 * npt + i] * sc[2 * npt + i];
    alpha = mixture[i] + (TB[3][0] - 1) * sc[1 * npt + i] +
            (TB[3][1] - 1) * sc[4 * npt + i] +
            (TB[3][2] - 1) * sc[9 * npt + i] +
            (TB[3][3] - 1) * sc[10 * npt + i];
    k_f = alpha * k_f_s[3 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[5 * npt + i];
    Kc = Kc_s[3 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2 * npt + i] -= qdot;
    wdot[2 * npt + i] -= qdot;
    wdot[5 * npt + i] += qdot;

    /*reaction 5: O + H + M <=> OH + M */
    phi_f = sc[0 * npt + i] * sc[2 * npt + i];
    alpha = mixture[i] + (TB[4][0] - 1) * sc[1 * npt + i] +
            (TB[4][1] - 1) * sc[4 * npt + i] +
            (TB[4][2] - 1) * sc[9 * npt + i] +
            (TB[4][3] - 1) * sc[10 * npt + i];
    k_f = alpha * k_f_s[4 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[3 * npt + i];
    Kc = Kc_s[4 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] -= qdot;
    wdot[2 * npt + i] -= qdot;
    wdot[3 * npt + i] += qdot;

    /*reaction 6: H2O + M <=> H + OH + M */
    phi_f = sc[4 * npt + i];
    alpha = mixture[i] + (TB[5][0] - 1) * sc[1 * npt + i] +
            (TB[5][1] - 1) * sc[4 * npt + i] +
            (TB[5][2] - 1) * sc[10 * npt + i] +
            (TB[5][3] - 1) * sc[8 * npt + i] + (TB[5][4] - 1) * sc[5 * npt + i];
    k_f = alpha * k_f_s[5 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[0 * npt + i] * sc[3 * npt + i];
    Kc = Kc_s[5 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] += qdot;
    wdot[3 * npt + i] += qdot;
    wdot[4 * npt + i] -= qdot;

    /*reaction 7: O + OH + M <=> HO2 + M */
    phi_f = sc[2 * npt + i] * sc[3 * npt + i];
    alpha = mixture[i] + (TB[6][0] - 1) * sc[1 * npt + i] +
            (TB[6][1] - 1) * sc[4 * npt + i] +
            (TB[6][2] - 1) * sc[9 * npt + i] +
            (TB[6][3] - 1) * sc[10 * npt + i];
    k_f = alpha * k_f_s[6 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[6 * npt + i];
    Kc = Kc_s[6 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2 * npt + i] -= qdot;
    wdot[3 * npt + i] -= qdot;
    wdot[6 * npt + i] += qdot;

    /*reaction 8: H + O2 <=> O + OH */
    phi_f = sc[0 * npt + i] * sc[5 * npt + i];
    k_f = k_f_s[7 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[2 * npt + i] * sc[3 * npt + i];
    Kc = Kc_s[7 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] -= qdot;
    wdot[2 * npt + i] += qdot;
    wdot[3 * npt + i] += qdot;
    wdot[5 * npt + i] -= qdot;

    /*reaction 9: O + H2 <=> H + OH */
    phi_f = sc[1 * npt + i] * sc[2 * npt + i];
    k_f = k_f_s[8 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[0 * npt + i] * sc[3 * npt + i];
    Kc = Kc_s[8 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] += qdot;
    wdot[1 * npt + i] -= qdot;
    wdot[2 * npt + i] -= qdot;
    wdot[3 * npt + i] += qdot;

    /*reaction 10: O + H2 <=> H + OH */
    phi_f = sc[1 * npt + i] * sc[2 * npt + i];
    k_f = k_f_s[9 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[0 * npt + i] * sc[3 * npt + i];
    Kc = Kc_s[9 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] += qdot;
    wdot[1 * npt + i] -= qdot;
    wdot[2 * npt + i] -= qdot;
    wdot[3 * npt + i] += qdot;

    /*reaction 11: H2 + OH <=> H2O + H */
    phi_f = sc[1 * npt + i] * sc[3 * npt + i];
    k_f = k_f_s[10 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[0 * npt + i] * sc[4 * npt + i];
    Kc = Kc_s[10 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] += qdot;
    wdot[1 * npt + i] -= qdot;
    wdot[3 * npt + i] -= qdot;
    wdot[4 * npt + i] += qdot;

    /*reaction 12: OH + OH <=> O + H2O */
    phi_f = sc[3 * npt + i] * sc[3 * npt + i];
    k_f = k_f_s[11 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[2 * npt + i] * sc[4 * npt + i];
    Kc = Kc_s[11 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2 * npt + i] += qdot;
    wdot[3 * npt + i] -= qdot;
    wdot[3 * npt + i] -= qdot;
    wdot[4 * npt + i] += qdot;

    /*reaction 13: H2 + AR <=> H + H + AR */
    phi_f = sc[1 * npt + i] * sc[9 * npt + i];
    k_f = k_f_s[12 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[0 * npt + i] * sc[0 * npt + i] * sc[9 * npt + i];
    Kc = Kc_s[12 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] += qdot;
    wdot[0 * npt + i] += qdot;
    wdot[1 * npt + i] -= qdot;
    wdot[9 * npt + i] -= qdot;
    wdot[9 * npt + i] += qdot;

    /*reaction 14: H2 + HE <=> H + H + HE */
    phi_f = sc[1 * npt + i] * sc[10 * npt + i];
    k_f = k_f_s[13 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[0 * npt + i] * sc[0 * npt + i] * sc[10 * npt + i];
    Kc = Kc_s[13 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] += qdot;
    wdot[0 * npt + i] += qdot;
    wdot[1 * npt + i] -= qdot;
    wdot[10 * npt + i] -= qdot;
    wdot[10 * npt + i] += qdot;

    /*reaction 15: O + O + AR <=> O2 + AR */
    phi_f = sc[2 * npt + i] * sc[2 * npt + i] * sc[9 * npt + i];
    k_f = k_f_s[14 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[5 * npt + i] * sc[9 * npt + i];
    Kc = Kc_s[14 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2 * npt + i] -= qdot;
    wdot[2 * npt + i] -= qdot;
    wdot[5 * npt + i] += qdot;
    wdot[9 * npt + i] -= qdot;
    wdot[9 * npt + i] += qdot;

    /*reaction 16: O + O + HE <=> O2 + HE */
    phi_f = sc[2 * npt + i] * sc[2 * npt + i] * sc[10 * npt + i];
    k_f = k_f_s[15 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[5 * npt + i] * sc[10 * npt + i];
    Kc = Kc_s[15 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2 * npt + i] -= qdot;
    wdot[2 * npt + i] -= qdot;
    wdot[5 * npt + i] += qdot;
    wdot[10 * npt + i] -= qdot;
    wdot[10 * npt + i] += qdot;

    /*reaction 17: H2O + H2O <=> H + OH + H2O */
    phi_f = sc[4 * npt + i] * sc[4 * npt + i];
    k_f = k_f_s[16 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[0 * npt + i] * sc[3 * npt + i] * sc[4 * npt + i];
    Kc = Kc_s[16 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] += qdot;
    wdot[3 * npt + i] += qdot;
    wdot[4 * npt + i] -= qdot;
    wdot[4 * npt + i] -= qdot;
    wdot[4 * npt + i] += qdot;

    /*reaction 18: HO2 + H <=> H2 + O2 */
    phi_f = sc[0 * npt + i] * sc[6 * npt + i];
    k_f = k_f_s[17 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[1 * npt + i] * sc[5 * npt + i];
    Kc = Kc_s[17 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] -= qdot;
    wdot[1 * npt + i] += qdot;
    wdot[5 * npt + i] += qdot;
    wdot[6 * npt + i] -= qdot;

    /*reaction 19: HO2 + H <=> OH + OH */
    phi_f = sc[0 * npt + i] * sc[6 * npt + i];
    k_f = k_f_s[18 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[3 * npt + i] * sc[3 * npt + i];
    Kc = Kc_s[18 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] -= qdot;
    wdot[3 * npt + i] += qdot;
    wdot[3 * npt + i] += qdot;
    wdot[6 * npt + i] -= qdot;

    /*reaction 20: HO2 + O <=> O2 + OH */
    phi_f = sc[2 * npt + i] * sc[6 * npt + i];
    k_f = k_f_s[19 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[3 * npt + i] * sc[5 * npt + i];
    Kc = Kc_s[19 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2 * npt + i] -= qdot;
    wdot[3 * npt + i] += qdot;
    wdot[5 * npt + i] += qdot;
    wdot[6 * npt + i] -= qdot;

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    phi_f = sc[3 * npt + i] * sc[6 * npt + i];
    k_f = k_f_s[20 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[4 * npt + i] * sc[5 * npt + i];
    Kc = Kc_s[20 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3 * npt + i] -= qdot;
    wdot[4 * npt + i] += qdot;
    wdot[5 * npt + i] += qdot;
    wdot[6 * npt + i] -= qdot;

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[6 * npt + i] * sc[6 * npt + i];
    k_f = k_f_s[21 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[5 * npt + i] * sc[7 * npt + i];
    Kc = Kc_s[21 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5 * npt + i] += qdot;
    wdot[6 * npt + i] -= qdot;
    wdot[6 * npt + i] -= qdot;
    wdot[7 * npt + i] += qdot;

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[6 * npt + i] * sc[6 * npt + i];
    k_f = k_f_s[22 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[5 * npt + i] * sc[7 * npt + i];
    Kc = Kc_s[22 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5 * npt + i] += qdot;
    wdot[6 * npt + i] -= qdot;
    wdot[6 * npt + i] -= qdot;
    wdot[7 * npt + i] += qdot;

    /*reaction 24: H2O2 + H <=> H2O + OH */
    phi_f = sc[0 * npt + i] * sc[7 * npt + i];
    k_f = k_f_s[23 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[3 * npt + i] * sc[4 * npt + i];
    Kc = Kc_s[23 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] -= qdot;
    wdot[3 * npt + i] += qdot;
    wdot[4 * npt + i] += qdot;
    wdot[7 * npt + i] -= qdot;

    /*reaction 25: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[0 * npt + i] * sc[7 * npt + i];
    k_f = k_f_s[24 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[1 * npt + i] * sc[6 * npt + i];
    Kc = Kc_s[24 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] -= qdot;
    wdot[1 * npt + i] += qdot;
    wdot[6 * npt + i] += qdot;
    wdot[7 * npt + i] -= qdot;

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    phi_f = sc[2 * npt + i] * sc[7 * npt + i];
    k_f = k_f_s[25 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[3 * npt + i] * sc[6 * npt + i];
    Kc = Kc_s[25 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2 * npt + i] -= qdot;
    wdot[3 * npt + i] += qdot;
    wdot[6 * npt + i] += qdot;
    wdot[7 * npt + i] -= qdot;

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[3 * npt + i] * sc[7 * npt + i];
    k_f = k_f_s[26 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[4 * npt + i] * sc[6 * npt + i];
    Kc = Kc_s[26 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3 * npt + i] -= qdot;
    wdot[4 * npt + i] += qdot;
    wdot[6 * npt + i] += qdot;
    wdot[7 * npt + i] -= qdot;

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[3 * npt + i] * sc[7 * npt + i];
    k_f = k_f_s[27 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[4 * npt + i] * sc[6 * npt + i];
    Kc = Kc_s[27 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3 * npt + i] -= qdot;
    wdot[4 * npt + i] += qdot;
    wdot[6 * npt + i] += qdot;
    wdot[7 * npt + i] -= qdot;

    /*reaction 29: HO2 + H <=> O + H2O */
    phi_f = sc[0 * npt + i] * sc[6 * npt + i];
    k_f = k_f_s[28 * npt + i];
    q_f = phi_f * k_f;
    phi_r = sc[2 * npt + i] * sc[4 * npt + i];
    Kc = Kc_s[28 * npt + i];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0 * npt + i] -= qdot;
    wdot[2 * npt + i] += qdot;
    wdot[4 * npt + i] += qdot;
    wdot[6 * npt + i] -= qdot;
  }
}

/*compute an approx to the reaction Jacobian */
void
DWDOT_PRECOND(double* J, double* sc, double* Tp, int* HP)
{
  double c[11];

  for (int k = 0; k < 11; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian_precond(J, c, *Tp, *HP);

  /* dwdot[k]/dT */
  /* dTdot/d[X] */
  for (int k = 0; k < 11; k++) {
    J[132 + k] *= 1.e-6;
    J[k * 12 + 11] *= 1.e6;
  }

  return;
}

/*compute the reaction Jacobian */
void
DWDOT(double* J, double* sc, double* Tp, int* consP)
{
  double c[11];

  for (int k = 0; k < 11; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian(J, c, *Tp, *consP);

  /* dwdot[k]/dT */
  /* dTdot/d[X] */
  for (int k = 0; k < 11; k++) {
    J[132 + k] *= 1.e-6;
    J[k * 12 + 11] *= 1.e6;
  }

  return;
}

/*compute the sparsity pattern Jacobian */
void
SPARSITY_INFO(int* nJdata, int* consP, int NCELLS)
{
  double c[11];
  double J[144];

  for (int k = 0; k < 11; k++) {
    c[k] = 1.0 / 11.000000;
  }

  aJacobian(J, c, 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 12; k++) {
    for (int l = 0; l < 12; l++) {
      if (J[12 * k + l] != 0.0) {
        nJdata_tmp = nJdata_tmp + 1;
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;

  return;
}

/*compute the sparsity pattern of simplified Jacobian */
void
SPARSITY_INFO_PRECOND(int* nJdata, int* consP)
{
  double c[11];
  double J[144];

  for (int k = 0; k < 11; k++) {
    c[k] = 1.0 / 11.000000;
  }

  aJacobian_precond(J, c, 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 12; k++) {
    for (int l = 0; l < 12; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (J[12 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  nJdata[0] = nJdata_tmp;

  return;
}

/*compute the sparsity pattern of the simplified precond Jacobian */
void
SPARSITY_PREPROC_PRECOND(int* rowVals, int* colPtrs, int* consP)
{
  double c[11];
  double J[144];

  for (int k = 0; k < 11; k++) {
    c[k] = 1.0 / 11.000000;
  }

  aJacobian_precond(J, c, 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 12; k++) {
    for (int l = 0; l < 12; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (J[12 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
    colPtrs[k + 1] = nJdata_tmp;
  }

  return;
}
/*compute the sparsity pattern of the Jacobian */
void
SPARSITY_PREPROC(int* rowVals, int* colPtrs, int* consP, int NCELLS)
{
  double c[11];
  double J[144];
  int offset_row;
  int offset_col;

  for (int k = 0; k < 11; k++) {
    c[k] = 1.0 / 11.000000;
  }

  aJacobian(J, c, 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    offset_row = nc * 12;
    offset_col = nc * 12;
    for (int k = 0; k < 12; k++) {
      for (int l = 0; l < 12; l++) {
        if (J[12 * k + l] != 0.0) {
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
void
aJacobian(double* J, double* sc, double T, int consP)
{
  for (int i = 0; i < 144; i++) {
    J[i] = 0.0;
  }

  double wdot[11];
  for (int k = 0; k < 11; k++) {
    wdot[k] = 0.0;
  }

  double tc[] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; /*temperature cache */
  double invT = 1.0 / tc[1];
  double invT2 = invT * invT;

  /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
  double refC = 101325 / 8.31451 / T;
  double refCinv = 1.0 / refC;

  /*compute the mixture concentration */
  double mixture = 0.0;
  for (int k = 0; k < 11; ++k) {
    mixture += sc[k];
  }

  /*compute the Gibbs free energy */
  double g_RT[11];
  gibbs(g_RT, tc);

  /*compute the species enthalpy */
  double h_RT[11];
  speciesEnthalpy(h_RT, tc);

  double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  double dqdci, dcdc_fac, dqdc[11];
  double Pr, fPr, F, k_0, logPr;
  double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  double Fcent1, Fcent2, Fcent3, Fcent;
  double dlogFdc, dlogFdn, dlogFdcn_fac;
  double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const double ln10 = log(10.0);
  const double log10e = 1.0 / log(10.0);
  /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
  /*a pressure-fall-off reaction */
  /* also 3-body */
  /* 3-body correction factor */
  alpha = mixture + (TB[0][0] - 1) * sc[1] + (TB[0][1] - 1) * sc[4] +
          (TB[0][2] - 1) * sc[5] + (TB[0][3] - 1) * sc[9] +
          (TB[0][4] - 1) * sc[10];
  /* forward */
  phi_f = sc[0] * sc[5];
  k_f = prefactor_units[0] * fwd_A[0] *
        exp(fwd_beta[0] * tc[0] - activation_units[0] * fwd_Ea[0] * invT);
  dlnkfdT = fwd_beta[0] * invT + activation_units[0] * fwd_Ea[0] * invT2;
  /* pressure-fall-off */
  k_0 = low_A[0] *
        exp(low_beta[0] * tc[0] - activation_units[0] * low_Ea[0] * invT);
  Pr = phase_units[0] * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = low_beta[0] * invT + activation_units[0] * low_Ea[0] * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  /* Troe form */
  logPr = log10(Pr);
  Fcent1 =
    (fabs(troe_Tsss[0]) > 1.e-100 ? (1. - troe_a[0]) * exp(-T / troe_Tsss[0])
                                  : 0.);
  Fcent2 = (fabs(troe_Ts[0]) > 1.e-100 ? troe_a[0] * exp(-T / troe_Ts[0]) : 0.);
  Fcent3 = (troe_len[0] == 4 ? exp(-troe_Tss[0] * invT) : 0.);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = pow(10.0, logFcent * troe);
  dlogFcentdT = log10e / Fcent *
                ((fabs(troe_Tsss[0]) > 1.e-100 ? -Fcent1 / troe_Tsss[0] : 0.) +
                 (fabs(troe_Ts[0]) > 1.e-100 ? -Fcent2 / troe_Ts[0] : 0.) +
                 (troe_len[0] == 4 ? Fcent3 * troe_Tss[0] * invT2 : 0.));
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  /* reverse */
  phi_r = sc[6];
  Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[6]) + 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[5] -= q; /* O2 */
  wdot[6] += q; /* HO2 */
  /* for convenience */
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP) {
    /* d()/d[H] */
    dqdci = +k_f * sc[5];
    J[0] -= dqdci; /* dwdot[H]/d[H] */
    J[5] -= dqdci; /* dwdot[O2]/d[H] */
    J[6] += dqdci; /* dwdot[HO2]/d[H] */
    /* d()/d[H2] */
    dqdci = (TB[0][0] - 1) * dcdc_fac;
    J[12] -= dqdci; /* dwdot[H]/d[H2] */
    J[17] -= dqdci; /* dwdot[O2]/d[H2] */
    J[18] += dqdci; /* dwdot[HO2]/d[H2] */
    /* d()/d[H2O] */
    dqdci = (TB[0][1] - 1) * dcdc_fac;
    J[48] -= dqdci; /* dwdot[H]/d[H2O] */
    J[53] -= dqdci; /* dwdot[O2]/d[H2O] */
    J[54] += dqdci; /* dwdot[HO2]/d[H2O] */
    /* d()/d[O2] */
    dqdci = (TB[0][2] - 1) * dcdc_fac + k_f * sc[0];
    J[60] -= dqdci; /* dwdot[H]/d[O2] */
    J[65] -= dqdci; /* dwdot[O2]/d[O2] */
    J[66] += dqdci; /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci = -k_r;
    J[72] -= dqdci; /* dwdot[H]/d[HO2] */
    J[77] -= dqdci; /* dwdot[O2]/d[HO2] */
    J[78] += dqdci; /* dwdot[HO2]/d[HO2] */
    /* d()/d[AR] */
    dqdci = (TB[0][3] - 1) * dcdc_fac;
    J[108] -= dqdci; /* dwdot[H]/d[AR] */
    J[113] -= dqdci; /* dwdot[O2]/d[AR] */
    J[114] += dqdci; /* dwdot[HO2]/d[AR] */
    /* d()/d[HE] */
    dqdci = (TB[0][4] - 1) * dcdc_fac;
    J[120] -= dqdci; /* dwdot[H]/d[HE] */
    J[125] -= dqdci; /* dwdot[O2]/d[HE] */
    J[126] += dqdci; /* dwdot[HO2]/d[HE] */
  } else {
    dqdc[0] = dcdc_fac + k_f * sc[5];
    dqdc[1] = TB[0][0] * dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = TB[0][1] * dcdc_fac;
    dqdc[5] = TB[0][2] * dcdc_fac + k_f * sc[0];
    dqdc[6] = dcdc_fac - k_r;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = TB[0][3] * dcdc_fac;
    dqdc[10] = TB[0][4] * dcdc_fac;
    for (int k = 0; k < 11; k++) {
      J[12 * k + 0] -= dqdc[k];
      J[12 * k + 5] -= dqdc[k];
      J[12 * k + 6] += dqdc[k];
    }
  }
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[137] -= dqdT; /* dwdot[O2]/dT */
  J[138] += dqdT; /* dwdot[HO2]/dT */

  /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
  /*a pressure-fall-off reaction */
  /* also 3-body */
  /* 3-body correction factor */
  alpha = mixture + (TB[1][0] - 1) * sc[4] + (TB[1][1] - 1) * sc[8] +
          (TB[1][2] - 1) * sc[5] + (TB[1][3] - 1) * sc[10] +
          (TB[1][4] - 1) * sc[7] + (TB[1][5] - 1) * sc[1];
  /* forward */
  phi_f = sc[7];
  k_f = prefactor_units[1] * fwd_A[1] *
        exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
  dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
  /* pressure-fall-off */
  k_0 = low_A[1] *
        exp(low_beta[1] * tc[0] - activation_units[1] * low_Ea[1] * invT);
  Pr = phase_units[1] * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = low_beta[1] * invT + activation_units[1] * low_Ea[1] * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  /* Troe form */
  logPr = log10(Pr);
  Fcent1 =
    (fabs(troe_Tsss[1]) > 1.e-100 ? (1. - troe_a[1]) * exp(-T / troe_Tsss[1])
                                  : 0.);
  Fcent2 = (fabs(troe_Ts[1]) > 1.e-100 ? troe_a[1] * exp(-T / troe_Ts[1]) : 0.);
  Fcent3 = (troe_len[1] == 4 ? exp(-troe_Tss[1] * invT) : 0.);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = pow(10.0, logFcent * troe);
  dlogFcentdT = log10e / Fcent *
                ((fabs(troe_Tsss[1]) > 1.e-100 ? -Fcent1 / troe_Tsss[1] : 0.) +
                 (fabs(troe_Ts[1]) > 1.e-100 ? -Fcent2 / troe_Ts[1] : 0.) +
                 (troe_len[1] == 4 ? Fcent3 * troe_Tss[1] * invT2 : 0.));
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  /* reverse */
  phi_r = sc[3] * sc[3];
  Kc = refC * exp(-g_RT[3] - g_RT[3] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[7]) + (2 * h_RT[3]) - 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  /* update wdot */
  wdot[3] += 2 * q; /* OH */
  wdot[7] -= q;     /* H2O2 */
  /* for convenience */
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = q / alpha * (1.0 / (Pr + 1.0) + dlogFdlogPr);
  if (consP) {
    /* d()/d[H2] */
    dqdci = (TB[1][5] - 1) * dcdc_fac;
    J[15] += 2 * dqdci; /* dwdot[OH]/d[H2] */
    J[19] -= dqdci;     /* dwdot[H2O2]/d[H2] */
    /* d()/d[OH] */
    dqdci = -k_r * 2 * sc[3];
    J[39] += 2 * dqdci; /* dwdot[OH]/d[OH] */
    J[43] -= dqdci;     /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci = (TB[1][0] - 1) * dcdc_fac;
    J[51] += 2 * dqdci; /* dwdot[OH]/d[H2O] */
    J[55] -= dqdci;     /* dwdot[H2O2]/d[H2O] */
    /* d()/d[O2] */
    dqdci = (TB[1][2] - 1) * dcdc_fac;
    J[63] += 2 * dqdci; /* dwdot[OH]/d[O2] */
    J[67] -= dqdci;     /* dwdot[H2O2]/d[O2] */
    /* d()/d[H2O2] */
    dqdci = (TB[1][4] - 1) * dcdc_fac + k_f;
    J[87] += 2 * dqdci; /* dwdot[OH]/d[H2O2] */
    J[91] -= dqdci;     /* dwdot[H2O2]/d[H2O2] */
    /* d()/d[N2] */
    dqdci = (TB[1][1] - 1) * dcdc_fac;
    J[99] += 2 * dqdci; /* dwdot[OH]/d[N2] */
    J[103] -= dqdci;    /* dwdot[H2O2]/d[N2] */
    /* d()/d[HE] */
    dqdci = (TB[1][3] - 1) * dcdc_fac;
    J[123] += 2 * dqdci; /* dwdot[OH]/d[HE] */
    J[127] -= dqdci;     /* dwdot[H2O2]/d[HE] */
  } else {
    dqdc[0] = dcdc_fac;
    dqdc[1] = TB[1][5] * dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac - k_r * 2 * sc[3];
    dqdc[4] = TB[1][0] * dcdc_fac;
    dqdc[5] = TB[1][2] * dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = TB[1][4] * dcdc_fac + k_f;
    dqdc[8] = TB[1][1] * dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = TB[1][3] * dcdc_fac;
    for (int k = 0; k < 11; k++) {
      J[12 * k + 3] += 2 * dqdc[k];
      J[12 * k + 7] -= dqdc[k];
    }
  }
  J[135] += 2 * dqdT; /* dwdot[OH]/dT */
  J[139] -= dqdT;     /* dwdot[H2O2]/dT */

  /*reaction 3: H2 + M <=> H + H + M */
  /*a third-body and non-pressure-fall-off reaction */
  /* 3-body correction factor */
  alpha = mixture + (TB[2][0] - 1) * sc[1] + (TB[2][1] - 1) * sc[4] +
          (TB[2][2] - 1) * sc[9] + (TB[2][3] - 1) * sc[10];
  /* forward */
  phi_f = sc[1];
  k_f = prefactor_units[2] * fwd_A[2] *
        exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
  dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[0];
  Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1]) + (2 * h_RT[0]) - 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += 2 * q; /* H */
  wdot[1] -= q;     /* H2 */
  /* for convenience */
  k_f *= alpha;
  k_r *= alpha;
  if (consP) {
    /* d()/d[H] */
    dqdci = -k_r * 2 * sc[0];
    J[0] += 2 * dqdci; /* dwdot[H]/d[H] */
    J[1] -= dqdci;     /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci = (TB[2][0] - 1) * q_nocor + k_f;
    J[12] += 2 * dqdci; /* dwdot[H]/d[H2] */
    J[13] -= dqdci;     /* dwdot[H2]/d[H2] */
    /* d()/d[H2O] */
    dqdci = (TB[2][1] - 1) * q_nocor;
    J[48] += 2 * dqdci; /* dwdot[H]/d[H2O] */
    J[49] -= dqdci;     /* dwdot[H2]/d[H2O] */
    /* d()/d[AR] */
    dqdci = (TB[2][2] - 1) * q_nocor;
    J[108] += 2 * dqdci; /* dwdot[H]/d[AR] */
    J[109] -= dqdci;     /* dwdot[H2]/d[AR] */
    /* d()/d[HE] */
    dqdci = (TB[2][3] - 1) * q_nocor;
    J[120] += 2 * dqdci; /* dwdot[H]/d[HE] */
    J[121] -= dqdci;     /* dwdot[H2]/d[HE] */
  } else {
    dqdc[0] = q_nocor - k_r * 2 * sc[0];
    dqdc[1] = TB[2][0] * q_nocor + k_f;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = TB[2][1] * q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = TB[2][2] * q_nocor;
    dqdc[10] = TB[2][3] * q_nocor;
    for (int k = 0; k < 11; k++) {
      J[12 * k + 0] += 2 * dqdc[k];
      J[12 * k + 1] -= dqdc[k];
    }
  }
  J[132] += 2 * dqdT; /* dwdot[H]/dT */
  J[133] -= dqdT;     /* dwdot[H2]/dT */

  /*reaction 4: O + O + M <=> O2 + M */
  /*a third-body and non-pressure-fall-off reaction */
  /* 3-body correction factor */
  alpha = mixture + (TB[3][0] - 1) * sc[1] + (TB[3][1] - 1) * sc[4] +
          (TB[3][2] - 1) * sc[9] + (TB[3][3] - 1) * sc[10];
  /* forward */
  phi_f = sc[2] * sc[2];
  k_f = prefactor_units[3] * fwd_A[3] *
        exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
  dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
  /* reverse */
  phi_r = sc[5];
  Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[2]) + (h_RT[5]) + 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] -= 2 * q; /* O */
  wdot[5] += q;     /* O2 */
  /* for convenience */
  k_f *= alpha;
  k_r *= alpha;
  if (consP) {
    /* d()/d[H2] */
    dqdci = (TB[3][0] - 1) * q_nocor;
    J[14] += -2 * dqdci; /* dwdot[O]/d[H2] */
    J[17] += dqdci;      /* dwdot[O2]/d[H2] */
    /* d()/d[O] */
    dqdci = +k_f * 2 * sc[2];
    J[26] += -2 * dqdci; /* dwdot[O]/d[O] */
    J[29] += dqdci;      /* dwdot[O2]/d[O] */
    /* d()/d[H2O] */
    dqdci = (TB[3][1] - 1) * q_nocor;
    J[50] += -2 * dqdci; /* dwdot[O]/d[H2O] */
    J[53] += dqdci;      /* dwdot[O2]/d[H2O] */
    /* d()/d[O2] */
    dqdci = -k_r;
    J[62] += -2 * dqdci; /* dwdot[O]/d[O2] */
    J[65] += dqdci;      /* dwdot[O2]/d[O2] */
    /* d()/d[AR] */
    dqdci = (TB[3][2] - 1) * q_nocor;
    J[110] += -2 * dqdci; /* dwdot[O]/d[AR] */
    J[113] += dqdci;      /* dwdot[O2]/d[AR] */
    /* d()/d[HE] */
    dqdci = (TB[3][3] - 1) * q_nocor;
    J[122] += -2 * dqdci; /* dwdot[O]/d[HE] */
    J[125] += dqdci;      /* dwdot[O2]/d[HE] */
  } else {
    dqdc[0] = q_nocor;
    dqdc[1] = TB[3][0] * q_nocor;
    dqdc[2] = q_nocor + k_f * 2 * sc[2];
    dqdc[3] = q_nocor;
    dqdc[4] = TB[3][1] * q_nocor;
    dqdc[5] = q_nocor - k_r;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = TB[3][2] * q_nocor;
    dqdc[10] = TB[3][3] * q_nocor;
    for (int k = 0; k < 11; k++) {
      J[12 * k + 2] += -2 * dqdc[k];
      J[12 * k + 5] += dqdc[k];
    }
  }
  J[134] += -2 * dqdT; /* dwdot[O]/dT */
  J[137] += dqdT;      /* dwdot[O2]/dT */

  /*reaction 5: O + H + M <=> OH + M */
  /*a third-body and non-pressure-fall-off reaction */
  /* 3-body correction factor */
  alpha = mixture + (TB[4][0] - 1) * sc[1] + (TB[4][1] - 1) * sc[4] +
          (TB[4][2] - 1) * sc[9] + (TB[4][3] - 1) * sc[10];
  /* forward */
  phi_f = sc[0] * sc[2];
  k_f = prefactor_units[4] * fwd_A[4] *
        exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
  dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
  /* reverse */
  phi_r = sc[3];
  Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[3]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3]) + 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[2] -= q; /* O */
  wdot[3] += q; /* OH */
  /* for convenience */
  k_f *= alpha;
  k_r *= alpha;
  if (consP) {
    /* d()/d[H] */
    dqdci = +k_f * sc[2];
    J[0] -= dqdci; /* dwdot[H]/d[H] */
    J[2] -= dqdci; /* dwdot[O]/d[H] */
    J[3] += dqdci; /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci = (TB[4][0] - 1) * q_nocor;
    J[12] -= dqdci; /* dwdot[H]/d[H2] */
    J[14] -= dqdci; /* dwdot[O]/d[H2] */
    J[15] += dqdci; /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci = +k_f * sc[0];
    J[24] -= dqdci; /* dwdot[H]/d[O] */
    J[26] -= dqdci; /* dwdot[O]/d[O] */
    J[27] += dqdci; /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci = -k_r;
    J[36] -= dqdci; /* dwdot[H]/d[OH] */
    J[38] -= dqdci; /* dwdot[O]/d[OH] */
    J[39] += dqdci; /* dwdot[OH]/d[OH] */
    /* d()/d[H2O] */
    dqdci = (TB[4][1] - 1) * q_nocor;
    J[48] -= dqdci; /* dwdot[H]/d[H2O] */
    J[50] -= dqdci; /* dwdot[O]/d[H2O] */
    J[51] += dqdci; /* dwdot[OH]/d[H2O] */
    /* d()/d[AR] */
    dqdci = (TB[4][2] - 1) * q_nocor;
    J[108] -= dqdci; /* dwdot[H]/d[AR] */
    J[110] -= dqdci; /* dwdot[O]/d[AR] */
    J[111] += dqdci; /* dwdot[OH]/d[AR] */
    /* d()/d[HE] */
    dqdci = (TB[4][3] - 1) * q_nocor;
    J[120] -= dqdci; /* dwdot[H]/d[HE] */
    J[122] -= dqdci; /* dwdot[O]/d[HE] */
    J[123] += dqdci; /* dwdot[OH]/d[HE] */
  } else {
    dqdc[0] = q_nocor + k_f * sc[2];
    dqdc[1] = TB[4][0] * q_nocor;
    dqdc[2] = q_nocor + k_f * sc[0];
    dqdc[3] = q_nocor - k_r;
    dqdc[4] = TB[4][1] * q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = TB[4][2] * q_nocor;
    dqdc[10] = TB[4][3] * q_nocor;
    for (int k = 0; k < 11; k++) {
      J[12 * k + 0] -= dqdc[k];
      J[12 * k + 2] -= dqdc[k];
      J[12 * k + 3] += dqdc[k];
    }
  }
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[134] -= dqdT; /* dwdot[O]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */

  /*reaction 6: H2O + M <=> H + OH + M */
  /*a third-body and non-pressure-fall-off reaction */
  /* 3-body correction factor */
  alpha = mixture + (TB[5][0] - 1) * sc[1] + (TB[5][1] - 1) * sc[4] +
          (TB[5][2] - 1) * sc[10] + (TB[5][3] - 1) * sc[8] +
          (TB[5][4] - 1) * sc[5];
  /* forward */
  phi_f = sc[4];
  k_f = prefactor_units[5] * fwd_A[5] *
        exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
  dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[3];
  Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4]) + (h_RT[0] + h_RT[3]) - 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += q; /* H */
  wdot[3] += q; /* OH */
  wdot[4] -= q; /* H2O */
  /* for convenience */
  k_f *= alpha;
  k_r *= alpha;
  if (consP) {
    /* d()/d[H] */
    dqdci = -k_r * sc[3];
    J[0] += dqdci; /* dwdot[H]/d[H] */
    J[3] += dqdci; /* dwdot[OH]/d[H] */
    J[4] -= dqdci; /* dwdot[H2O]/d[H] */
    /* d()/d[H2] */
    dqdci = (TB[5][0] - 1) * q_nocor;
    J[12] += dqdci; /* dwdot[H]/d[H2] */
    J[15] += dqdci; /* dwdot[OH]/d[H2] */
    J[16] -= dqdci; /* dwdot[H2O]/d[H2] */
    /* d()/d[OH] */
    dqdci = -k_r * sc[0];
    J[36] += dqdci; /* dwdot[H]/d[OH] */
    J[39] += dqdci; /* dwdot[OH]/d[OH] */
    J[40] -= dqdci; /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci = (TB[5][1] - 1) * q_nocor + k_f;
    J[48] += dqdci; /* dwdot[H]/d[H2O] */
    J[51] += dqdci; /* dwdot[OH]/d[H2O] */
    J[52] -= dqdci; /* dwdot[H2O]/d[H2O] */
    /* d()/d[O2] */
    dqdci = (TB[5][4] - 1) * q_nocor;
    J[60] += dqdci; /* dwdot[H]/d[O2] */
    J[63] += dqdci; /* dwdot[OH]/d[O2] */
    J[64] -= dqdci; /* dwdot[H2O]/d[O2] */
    /* d()/d[N2] */
    dqdci = (TB[5][3] - 1) * q_nocor;
    J[96] += dqdci;  /* dwdot[H]/d[N2] */
    J[99] += dqdci;  /* dwdot[OH]/d[N2] */
    J[100] -= dqdci; /* dwdot[H2O]/d[N2] */
    /* d()/d[HE] */
    dqdci = (TB[5][2] - 1) * q_nocor;
    J[120] += dqdci; /* dwdot[H]/d[HE] */
    J[123] += dqdci; /* dwdot[OH]/d[HE] */
    J[124] -= dqdci; /* dwdot[H2O]/d[HE] */
  } else {
    dqdc[0] = q_nocor - k_r * sc[3];
    dqdc[1] = TB[5][0] * q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor - k_r * sc[0];
    dqdc[4] = TB[5][1] * q_nocor + k_f;
    dqdc[5] = TB[5][4] * q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = TB[5][3] * q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = TB[5][2] * q_nocor;
    for (int k = 0; k < 11; k++) {
      J[12 * k + 0] += dqdc[k];
      J[12 * k + 3] += dqdc[k];
      J[12 * k + 4] -= dqdc[k];
    }
  }
  J[132] += dqdT; /* dwdot[H]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */
  J[136] -= dqdT; /* dwdot[H2O]/dT */

  /*reaction 7: O + OH + M <=> HO2 + M */
  /*a third-body and non-pressure-fall-off reaction */
  /* 3-body correction factor */
  alpha = mixture + (TB[6][0] - 1) * sc[1] + (TB[6][1] - 1) * sc[4] +
          (TB[6][2] - 1) * sc[9] + (TB[6][3] - 1) * sc[10];
  /* forward */
  phi_f = sc[2] * sc[3];
  k_f = prefactor_units[6] * fwd_A[6] *
        exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
  dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
  /* reverse */
  phi_r = sc[6];
  Kc = refCinv * exp(g_RT[2] + g_RT[3] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[6]) + 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] -= q; /* O */
  wdot[3] -= q; /* OH */
  wdot[6] += q; /* HO2 */
  /* for convenience */
  k_f *= alpha;
  k_r *= alpha;
  if (consP) {
    /* d()/d[H2] */
    dqdci = (TB[6][0] - 1) * q_nocor;
    J[14] -= dqdci; /* dwdot[O]/d[H2] */
    J[15] -= dqdci; /* dwdot[OH]/d[H2] */
    J[18] += dqdci; /* dwdot[HO2]/d[H2] */
    /* d()/d[O] */
    dqdci = +k_f * sc[3];
    J[26] -= dqdci; /* dwdot[O]/d[O] */
    J[27] -= dqdci; /* dwdot[OH]/d[O] */
    J[30] += dqdci; /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci = +k_f * sc[2];
    J[38] -= dqdci; /* dwdot[O]/d[OH] */
    J[39] -= dqdci; /* dwdot[OH]/d[OH] */
    J[42] += dqdci; /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci = (TB[6][1] - 1) * q_nocor;
    J[50] -= dqdci; /* dwdot[O]/d[H2O] */
    J[51] -= dqdci; /* dwdot[OH]/d[H2O] */
    J[54] += dqdci; /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci = -k_r;
    J[74] -= dqdci; /* dwdot[O]/d[HO2] */
    J[75] -= dqdci; /* dwdot[OH]/d[HO2] */
    J[78] += dqdci; /* dwdot[HO2]/d[HO2] */
    /* d()/d[AR] */
    dqdci = (TB[6][2] - 1) * q_nocor;
    J[110] -= dqdci; /* dwdot[O]/d[AR] */
    J[111] -= dqdci; /* dwdot[OH]/d[AR] */
    J[114] += dqdci; /* dwdot[HO2]/d[AR] */
    /* d()/d[HE] */
    dqdci = (TB[6][3] - 1) * q_nocor;
    J[122] -= dqdci; /* dwdot[O]/d[HE] */
    J[123] -= dqdci; /* dwdot[OH]/d[HE] */
    J[126] += dqdci; /* dwdot[HO2]/d[HE] */
  } else {
    dqdc[0] = q_nocor;
    dqdc[1] = TB[6][0] * q_nocor;
    dqdc[2] = q_nocor + k_f * sc[3];
    dqdc[3] = q_nocor + k_f * sc[2];
    dqdc[4] = TB[6][1] * q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor - k_r;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = TB[6][2] * q_nocor;
    dqdc[10] = TB[6][3] * q_nocor;
    for (int k = 0; k < 11; k++) {
      J[12 * k + 2] -= dqdc[k];
      J[12 * k + 3] -= dqdc[k];
      J[12 * k + 6] += dqdc[k];
    }
  }
  J[134] -= dqdT; /* dwdot[O]/dT */
  J[135] -= dqdT; /* dwdot[OH]/dT */
  J[138] += dqdT; /* dwdot[HO2]/dT */

  /*reaction 8: H + O2 <=> O + OH */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[0] * sc[5];
  k_f = prefactor_units[7] * fwd_A[7] *
        exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
  dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
  /* reverse */
  phi_r = sc[2] * sc[3];
  Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[2] += q; /* O */
  wdot[3] += q; /* OH */
  wdot[5] -= q; /* O2 */
  /* d()/d[H] */
  dqdci = +k_f * sc[5];
  J[0] -= dqdci; /* dwdot[H]/d[H] */
  J[2] += dqdci; /* dwdot[O]/d[H] */
  J[3] += dqdci; /* dwdot[OH]/d[H] */
  J[5] -= dqdci; /* dwdot[O2]/d[H] */
  /* d()/d[O] */
  dqdci = -k_r * sc[3];
  J[24] -= dqdci; /* dwdot[H]/d[O] */
  J[26] += dqdci; /* dwdot[O]/d[O] */
  J[27] += dqdci; /* dwdot[OH]/d[O] */
  J[29] -= dqdci; /* dwdot[O2]/d[O] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[2];
  J[36] -= dqdci; /* dwdot[H]/d[OH] */
  J[38] += dqdci; /* dwdot[O]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  J[41] -= dqdci; /* dwdot[O2]/d[OH] */
  /* d()/d[O2] */
  dqdci = +k_f * sc[0];
  J[60] -= dqdci; /* dwdot[H]/d[O2] */
  J[62] += dqdci; /* dwdot[O]/d[O2] */
  J[63] += dqdci; /* dwdot[OH]/d[O2] */
  J[65] -= dqdci; /* dwdot[O2]/d[O2] */
  /* d()/dT */
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[134] += dqdT; /* dwdot[O]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */
  J[137] -= dqdT; /* dwdot[O2]/dT */

  /*reaction 9: O + H2 <=> H + OH */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[1] * sc[2];
  k_f = prefactor_units[8] * fwd_A[8] *
        exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
  dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[3];
  Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += q; /* H */
  wdot[1] -= q; /* H2 */
  wdot[2] -= q; /* O */
  wdot[3] += q; /* OH */
  /* d()/d[H] */
  dqdci = -k_r * sc[3];
  J[0] += dqdci; /* dwdot[H]/d[H] */
  J[1] -= dqdci; /* dwdot[H2]/d[H] */
  J[2] -= dqdci; /* dwdot[O]/d[H] */
  J[3] += dqdci; /* dwdot[OH]/d[H] */
  /* d()/d[H2] */
  dqdci = +k_f * sc[2];
  J[12] += dqdci; /* dwdot[H]/d[H2] */
  J[13] -= dqdci; /* dwdot[H2]/d[H2] */
  J[14] -= dqdci; /* dwdot[O]/d[H2] */
  J[15] += dqdci; /* dwdot[OH]/d[H2] */
  /* d()/d[O] */
  dqdci = +k_f * sc[1];
  J[24] += dqdci; /* dwdot[H]/d[O] */
  J[25] -= dqdci; /* dwdot[H2]/d[O] */
  J[26] -= dqdci; /* dwdot[O]/d[O] */
  J[27] += dqdci; /* dwdot[OH]/d[O] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[0];
  J[36] += dqdci; /* dwdot[H]/d[OH] */
  J[37] -= dqdci; /* dwdot[H2]/d[OH] */
  J[38] -= dqdci; /* dwdot[O]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  /* d()/dT */
  J[132] += dqdT; /* dwdot[H]/dT */
  J[133] -= dqdT; /* dwdot[H2]/dT */
  J[134] -= dqdT; /* dwdot[O]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */

  /*reaction 10: O + H2 <=> H + OH */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[1] * sc[2];
  k_f = prefactor_units[9] * fwd_A[9] *
        exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
  dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[3];
  Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += q; /* H */
  wdot[1] -= q; /* H2 */
  wdot[2] -= q; /* O */
  wdot[3] += q; /* OH */
  /* d()/d[H] */
  dqdci = -k_r * sc[3];
  J[0] += dqdci; /* dwdot[H]/d[H] */
  J[1] -= dqdci; /* dwdot[H2]/d[H] */
  J[2] -= dqdci; /* dwdot[O]/d[H] */
  J[3] += dqdci; /* dwdot[OH]/d[H] */
  /* d()/d[H2] */
  dqdci = +k_f * sc[2];
  J[12] += dqdci; /* dwdot[H]/d[H2] */
  J[13] -= dqdci; /* dwdot[H2]/d[H2] */
  J[14] -= dqdci; /* dwdot[O]/d[H2] */
  J[15] += dqdci; /* dwdot[OH]/d[H2] */
  /* d()/d[O] */
  dqdci = +k_f * sc[1];
  J[24] += dqdci; /* dwdot[H]/d[O] */
  J[25] -= dqdci; /* dwdot[H2]/d[O] */
  J[26] -= dqdci; /* dwdot[O]/d[O] */
  J[27] += dqdci; /* dwdot[OH]/d[O] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[0];
  J[36] += dqdci; /* dwdot[H]/d[OH] */
  J[37] -= dqdci; /* dwdot[H2]/d[OH] */
  J[38] -= dqdci; /* dwdot[O]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  /* d()/dT */
  J[132] += dqdT; /* dwdot[H]/dT */
  J[133] -= dqdT; /* dwdot[H2]/dT */
  J[134] -= dqdT; /* dwdot[O]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */

  /*reaction 11: H2 + OH <=> H2O + H */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[1] * sc[3];
  k_f = prefactor_units[10] * fwd_A[10] *
        exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
  dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[4];
  Kc = exp(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[0] + h_RT[4]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += q; /* H */
  wdot[1] -= q; /* H2 */
  wdot[3] -= q; /* OH */
  wdot[4] += q; /* H2O */
  /* d()/d[H] */
  dqdci = -k_r * sc[4];
  J[0] += dqdci; /* dwdot[H]/d[H] */
  J[1] -= dqdci; /* dwdot[H2]/d[H] */
  J[3] -= dqdci; /* dwdot[OH]/d[H] */
  J[4] += dqdci; /* dwdot[H2O]/d[H] */
  /* d()/d[H2] */
  dqdci = +k_f * sc[3];
  J[12] += dqdci; /* dwdot[H]/d[H2] */
  J[13] -= dqdci; /* dwdot[H2]/d[H2] */
  J[15] -= dqdci; /* dwdot[OH]/d[H2] */
  J[16] += dqdci; /* dwdot[H2O]/d[H2] */
  /* d()/d[OH] */
  dqdci = +k_f * sc[1];
  J[36] += dqdci; /* dwdot[H]/d[OH] */
  J[37] -= dqdci; /* dwdot[H2]/d[OH] */
  J[39] -= dqdci; /* dwdot[OH]/d[OH] */
  J[40] += dqdci; /* dwdot[H2O]/d[OH] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[0];
  J[48] += dqdci; /* dwdot[H]/d[H2O] */
  J[49] -= dqdci; /* dwdot[H2]/d[H2O] */
  J[51] -= dqdci; /* dwdot[OH]/d[H2O] */
  J[52] += dqdci; /* dwdot[H2O]/d[H2O] */
  /* d()/dT */
  J[132] += dqdT; /* dwdot[H]/dT */
  J[133] -= dqdT; /* dwdot[H2]/dT */
  J[135] -= dqdT; /* dwdot[OH]/dT */
  J[136] += dqdT; /* dwdot[H2O]/dT */

  /*reaction 12: OH + OH <=> O + H2O */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[3] * sc[3];
  k_f = prefactor_units[11] * fwd_A[11] *
        exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
  dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
  /* reverse */
  phi_r = sc[2] * sc[4];
  Kc = exp(-g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[3]) + (h_RT[2] + h_RT[4]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] += q;     /* O */
  wdot[3] -= 2 * q; /* OH */
  wdot[4] += q;     /* H2O */
  /* d()/d[O] */
  dqdci = -k_r * sc[4];
  J[26] += dqdci;      /* dwdot[O]/d[O] */
  J[27] += -2 * dqdci; /* dwdot[OH]/d[O] */
  J[28] += dqdci;      /* dwdot[H2O]/d[O] */
  /* d()/d[OH] */
  dqdci = +k_f * 2 * sc[3];
  J[38] += dqdci;      /* dwdot[O]/d[OH] */
  J[39] += -2 * dqdci; /* dwdot[OH]/d[OH] */
  J[40] += dqdci;      /* dwdot[H2O]/d[OH] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[2];
  J[50] += dqdci;      /* dwdot[O]/d[H2O] */
  J[51] += -2 * dqdci; /* dwdot[OH]/d[H2O] */
  J[52] += dqdci;      /* dwdot[H2O]/d[H2O] */
  /* d()/dT */
  J[134] += dqdT;      /* dwdot[O]/dT */
  J[135] += -2 * dqdT; /* dwdot[OH]/dT */
  J[136] += dqdT;      /* dwdot[H2O]/dT */

  /*reaction 13: H2 + AR <=> H + H + AR */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[1] * sc[9];
  k_f = prefactor_units[12] * fwd_A[12] *
        exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
  dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[0] * sc[9];
  Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (2 * h_RT[0] + h_RT[9]) - 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += 2 * q; /* H */
  wdot[1] -= q;     /* H2 */
  /* d()/d[H] */
  dqdci = -k_r * 2 * sc[0] * sc[9];
  J[0] += 2 * dqdci; /* dwdot[H]/d[H] */
  J[1] -= dqdci;     /* dwdot[H2]/d[H] */
  /* d()/d[H2] */
  dqdci = +k_f * sc[9];
  J[12] += 2 * dqdci; /* dwdot[H]/d[H2] */
  J[13] -= dqdci;     /* dwdot[H2]/d[H2] */
  /* d()/d[AR] */
  dqdci = +k_f * sc[1] - k_r * sc[0] * sc[0];
  J[108] += 2 * dqdci; /* dwdot[H]/d[AR] */
  J[109] -= dqdci;     /* dwdot[H2]/d[AR] */
  /* d()/dT */
  J[132] += 2 * dqdT; /* dwdot[H]/dT */
  J[133] -= dqdT;     /* dwdot[H2]/dT */

  /*reaction 14: H2 + HE <=> H + H + HE */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[1] * sc[10];
  k_f = prefactor_units[13] * fwd_A[13] *
        exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
  dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[0] * sc[10];
  Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (2 * h_RT[0] + h_RT[10]) - 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += 2 * q; /* H */
  wdot[1] -= q;     /* H2 */
  /* d()/d[H] */
  dqdci = -k_r * 2 * sc[0] * sc[10];
  J[0] += 2 * dqdci; /* dwdot[H]/d[H] */
  J[1] -= dqdci;     /* dwdot[H2]/d[H] */
  /* d()/d[H2] */
  dqdci = +k_f * sc[10];
  J[12] += 2 * dqdci; /* dwdot[H]/d[H2] */
  J[13] -= dqdci;     /* dwdot[H2]/d[H2] */
  /* d()/d[HE] */
  dqdci = +k_f * sc[1] - k_r * sc[0] * sc[0];
  J[120] += 2 * dqdci; /* dwdot[H]/d[HE] */
  J[121] -= dqdci;     /* dwdot[H2]/d[HE] */
  /* d()/dT */
  J[132] += 2 * dqdT; /* dwdot[H]/dT */
  J[133] -= dqdT;     /* dwdot[H2]/dT */

  /*reaction 15: O + O + AR <=> O2 + AR */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[2] * sc[2] * sc[9];
  k_f = prefactor_units[14] * fwd_A[14] *
        exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
  dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
  /* reverse */
  phi_r = sc[5] * sc[9];
  Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[2] + h_RT[9]) + (h_RT[5] + h_RT[9]) + 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] -= 2 * q; /* O */
  wdot[5] += q;     /* O2 */
  /* d()/d[O] */
  dqdci = +k_f * 2 * sc[2] * sc[9];
  J[26] += -2 * dqdci; /* dwdot[O]/d[O] */
  J[29] += dqdci;      /* dwdot[O2]/d[O] */
  /* d()/d[O2] */
  dqdci = -k_r * sc[9];
  J[62] += -2 * dqdci; /* dwdot[O]/d[O2] */
  J[65] += dqdci;      /* dwdot[O2]/d[O2] */
  /* d()/d[AR] */
  dqdci = +k_f * sc[2] * sc[2] - k_r * sc[5];
  J[110] += -2 * dqdci; /* dwdot[O]/d[AR] */
  J[113] += dqdci;      /* dwdot[O2]/d[AR] */
  /* d()/dT */
  J[134] += -2 * dqdT; /* dwdot[O]/dT */
  J[137] += dqdT;      /* dwdot[O2]/dT */

  /*reaction 16: O + O + HE <=> O2 + HE */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[2] * sc[2] * sc[10];
  k_f = prefactor_units[15] * fwd_A[15] *
        exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
  dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
  /* reverse */
  phi_r = sc[5] * sc[10];
  Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[2] + h_RT[10]) + (h_RT[5] + h_RT[10]) + 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] -= 2 * q; /* O */
  wdot[5] += q;     /* O2 */
  /* d()/d[O] */
  dqdci = +k_f * 2 * sc[2] * sc[10];
  J[26] += -2 * dqdci; /* dwdot[O]/d[O] */
  J[29] += dqdci;      /* dwdot[O2]/d[O] */
  /* d()/d[O2] */
  dqdci = -k_r * sc[10];
  J[62] += -2 * dqdci; /* dwdot[O]/d[O2] */
  J[65] += dqdci;      /* dwdot[O2]/d[O2] */
  /* d()/d[HE] */
  dqdci = +k_f * sc[2] * sc[2] - k_r * sc[5];
  J[122] += -2 * dqdci; /* dwdot[O]/d[HE] */
  J[125] += dqdci;      /* dwdot[O2]/d[HE] */
  /* d()/dT */
  J[134] += -2 * dqdT; /* dwdot[O]/dT */
  J[137] += dqdT;      /* dwdot[O2]/dT */

  /*reaction 17: H2O + H2O <=> H + OH + H2O */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[4] * sc[4];
  k_f = prefactor_units[16] * fwd_A[16] *
        exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
  dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[3] * sc[4];
  Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[4]) + (h_RT[0] + h_RT[3] + h_RT[4]) - 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += q; /* H */
  wdot[3] += q; /* OH */
  wdot[4] -= q; /* H2O */
  /* d()/d[H] */
  dqdci = -k_r * sc[3] * sc[4];
  J[0] += dqdci; /* dwdot[H]/d[H] */
  J[3] += dqdci; /* dwdot[OH]/d[H] */
  J[4] -= dqdci; /* dwdot[H2O]/d[H] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[0] * sc[4];
  J[36] += dqdci; /* dwdot[H]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  J[40] -= dqdci; /* dwdot[H2O]/d[OH] */
  /* d()/d[H2O] */
  dqdci = +k_f * 2 * sc[4] - k_r * sc[0] * sc[3];
  J[48] += dqdci; /* dwdot[H]/d[H2O] */
  J[51] += dqdci; /* dwdot[OH]/d[H2O] */
  J[52] -= dqdci; /* dwdot[H2O]/d[H2O] */
  /* d()/dT */
  J[132] += dqdT; /* dwdot[H]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */
  J[136] -= dqdT; /* dwdot[H2O]/dT */

  /*reaction 18: HO2 + H <=> H2 + O2 */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[0] * sc[6];
  k_f = prefactor_units[17] * fwd_A[17] *
        exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
  dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
  /* reverse */
  phi_r = sc[1] * sc[5];
  Kc = exp(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[1] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[1] += q; /* H2 */
  wdot[5] += q; /* O2 */
  wdot[6] -= q; /* HO2 */
  /* d()/d[H] */
  dqdci = +k_f * sc[6];
  J[0] -= dqdci; /* dwdot[H]/d[H] */
  J[1] += dqdci; /* dwdot[H2]/d[H] */
  J[5] += dqdci; /* dwdot[O2]/d[H] */
  J[6] -= dqdci; /* dwdot[HO2]/d[H] */
  /* d()/d[H2] */
  dqdci = -k_r * sc[5];
  J[12] -= dqdci; /* dwdot[H]/d[H2] */
  J[13] += dqdci; /* dwdot[H2]/d[H2] */
  J[17] += dqdci; /* dwdot[O2]/d[H2] */
  J[18] -= dqdci; /* dwdot[HO2]/d[H2] */
  /* d()/d[O2] */
  dqdci = -k_r * sc[1];
  J[60] -= dqdci; /* dwdot[H]/d[O2] */
  J[61] += dqdci; /* dwdot[H2]/d[O2] */
  J[65] += dqdci; /* dwdot[O2]/d[O2] */
  J[66] -= dqdci; /* dwdot[HO2]/d[O2] */
  /* d()/d[HO2] */
  dqdci = +k_f * sc[0];
  J[72] -= dqdci; /* dwdot[H]/d[HO2] */
  J[73] += dqdci; /* dwdot[H2]/d[HO2] */
  J[77] += dqdci; /* dwdot[O2]/d[HO2] */
  J[78] -= dqdci; /* dwdot[HO2]/d[HO2] */
  /* d()/dT */
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[133] += dqdT; /* dwdot[H2]/dT */
  J[137] += dqdT; /* dwdot[O2]/dT */
  J[138] -= dqdT; /* dwdot[HO2]/dT */

  /*reaction 19: HO2 + H <=> OH + OH */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[0] * sc[6];
  k_f = prefactor_units[18] * fwd_A[18] *
        exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
  dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
  /* reverse */
  phi_r = sc[3] * sc[3];
  Kc = exp(g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (2 * h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q;     /* H */
  wdot[3] += 2 * q; /* OH */
  wdot[6] -= q;     /* HO2 */
  /* d()/d[H] */
  dqdci = +k_f * sc[6];
  J[0] -= dqdci;     /* dwdot[H]/d[H] */
  J[3] += 2 * dqdci; /* dwdot[OH]/d[H] */
  J[6] -= dqdci;     /* dwdot[HO2]/d[H] */
  /* d()/d[OH] */
  dqdci = -k_r * 2 * sc[3];
  J[36] -= dqdci;     /* dwdot[H]/d[OH] */
  J[39] += 2 * dqdci; /* dwdot[OH]/d[OH] */
  J[42] -= dqdci;     /* dwdot[HO2]/d[OH] */
  /* d()/d[HO2] */
  dqdci = +k_f * sc[0];
  J[72] -= dqdci;     /* dwdot[H]/d[HO2] */
  J[75] += 2 * dqdci; /* dwdot[OH]/d[HO2] */
  J[78] -= dqdci;     /* dwdot[HO2]/d[HO2] */
  /* d()/dT */
  J[132] -= dqdT;     /* dwdot[H]/dT */
  J[135] += 2 * dqdT; /* dwdot[OH]/dT */
  J[138] -= dqdT;     /* dwdot[HO2]/dT */

  /*reaction 20: HO2 + O <=> O2 + OH */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[2] * sc[6];
  k_f = prefactor_units[19] * fwd_A[19] *
        exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
  dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
  /* reverse */
  phi_r = sc[3] * sc[5];
  Kc = exp(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] -= q; /* O */
  wdot[3] += q; /* OH */
  wdot[5] += q; /* O2 */
  wdot[6] -= q; /* HO2 */
  /* d()/d[O] */
  dqdci = +k_f * sc[6];
  J[26] -= dqdci; /* dwdot[O]/d[O] */
  J[27] += dqdci; /* dwdot[OH]/d[O] */
  J[29] += dqdci; /* dwdot[O2]/d[O] */
  J[30] -= dqdci; /* dwdot[HO2]/d[O] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[5];
  J[38] -= dqdci; /* dwdot[O]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  J[41] += dqdci; /* dwdot[O2]/d[OH] */
  J[42] -= dqdci; /* dwdot[HO2]/d[OH] */
  /* d()/d[O2] */
  dqdci = -k_r * sc[3];
  J[62] -= dqdci; /* dwdot[O]/d[O2] */
  J[63] += dqdci; /* dwdot[OH]/d[O2] */
  J[65] += dqdci; /* dwdot[O2]/d[O2] */
  J[66] -= dqdci; /* dwdot[HO2]/d[O2] */
  /* d()/d[HO2] */
  dqdci = +k_f * sc[2];
  J[74] -= dqdci; /* dwdot[O]/d[HO2] */
  J[75] += dqdci; /* dwdot[OH]/d[HO2] */
  J[77] += dqdci; /* dwdot[O2]/d[HO2] */
  J[78] -= dqdci; /* dwdot[HO2]/d[HO2] */
  /* d()/dT */
  J[134] -= dqdT; /* dwdot[O]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */
  J[137] += dqdT; /* dwdot[O2]/dT */
  J[138] -= dqdT; /* dwdot[HO2]/dT */

  /*reaction 21: HO2 + OH <=> H2O + O2 */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[3] * sc[6];
  k_f = prefactor_units[20] * fwd_A[20] *
        exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
  dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
  /* reverse */
  phi_r = sc[4] * sc[5];
  Kc = exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[4] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[3] -= q; /* OH */
  wdot[4] += q; /* H2O */
  wdot[5] += q; /* O2 */
  wdot[6] -= q; /* HO2 */
  /* d()/d[OH] */
  dqdci = +k_f * sc[6];
  J[39] -= dqdci; /* dwdot[OH]/d[OH] */
  J[40] += dqdci; /* dwdot[H2O]/d[OH] */
  J[41] += dqdci; /* dwdot[O2]/d[OH] */
  J[42] -= dqdci; /* dwdot[HO2]/d[OH] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[5];
  J[51] -= dqdci; /* dwdot[OH]/d[H2O] */
  J[52] += dqdci; /* dwdot[H2O]/d[H2O] */
  J[53] += dqdci; /* dwdot[O2]/d[H2O] */
  J[54] -= dqdci; /* dwdot[HO2]/d[H2O] */
  /* d()/d[O2] */
  dqdci = -k_r * sc[4];
  J[63] -= dqdci; /* dwdot[OH]/d[O2] */
  J[64] += dqdci; /* dwdot[H2O]/d[O2] */
  J[65] += dqdci; /* dwdot[O2]/d[O2] */
  J[66] -= dqdci; /* dwdot[HO2]/d[O2] */
  /* d()/d[HO2] */
  dqdci = +k_f * sc[3];
  J[75] -= dqdci; /* dwdot[OH]/d[HO2] */
  J[76] += dqdci; /* dwdot[H2O]/d[HO2] */
  J[77] += dqdci; /* dwdot[O2]/d[HO2] */
  J[78] -= dqdci; /* dwdot[HO2]/d[HO2] */
  /* d()/dT */
  J[135] -= dqdT; /* dwdot[OH]/dT */
  J[136] += dqdT; /* dwdot[H2O]/dT */
  J[137] += dqdT; /* dwdot[O2]/dT */
  J[138] -= dqdT; /* dwdot[HO2]/dT */

  /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[6] * sc[6];
  k_f = prefactor_units[21] * fwd_A[21] *
        exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
  dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
  /* reverse */
  phi_r = sc[5] * sc[7];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[6]) + (h_RT[5] + h_RT[7]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[5] += q;     /* O2 */
  wdot[6] -= 2 * q; /* HO2 */
  wdot[7] += q;     /* H2O2 */
  /* d()/d[O2] */
  dqdci = -k_r * sc[7];
  J[65] += dqdci;      /* dwdot[O2]/d[O2] */
  J[66] += -2 * dqdci; /* dwdot[HO2]/d[O2] */
  J[67] += dqdci;      /* dwdot[H2O2]/d[O2] */
  /* d()/d[HO2] */
  dqdci = +k_f * 2 * sc[6];
  J[77] += dqdci;      /* dwdot[O2]/d[HO2] */
  J[78] += -2 * dqdci; /* dwdot[HO2]/d[HO2] */
  J[79] += dqdci;      /* dwdot[H2O2]/d[HO2] */
  /* d()/d[H2O2] */
  dqdci = -k_r * sc[5];
  J[89] += dqdci;      /* dwdot[O2]/d[H2O2] */
  J[90] += -2 * dqdci; /* dwdot[HO2]/d[H2O2] */
  J[91] += dqdci;      /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[137] += dqdT;      /* dwdot[O2]/dT */
  J[138] += -2 * dqdT; /* dwdot[HO2]/dT */
  J[139] += dqdT;      /* dwdot[H2O2]/dT */

  /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[6] * sc[6];
  k_f = prefactor_units[22] * fwd_A[22] *
        exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
  dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
  /* reverse */
  phi_r = sc[5] * sc[7];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[6]) + (h_RT[5] + h_RT[7]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[5] += q;     /* O2 */
  wdot[6] -= 2 * q; /* HO2 */
  wdot[7] += q;     /* H2O2 */
  /* d()/d[O2] */
  dqdci = -k_r * sc[7];
  J[65] += dqdci;      /* dwdot[O2]/d[O2] */
  J[66] += -2 * dqdci; /* dwdot[HO2]/d[O2] */
  J[67] += dqdci;      /* dwdot[H2O2]/d[O2] */
  /* d()/d[HO2] */
  dqdci = +k_f * 2 * sc[6];
  J[77] += dqdci;      /* dwdot[O2]/d[HO2] */
  J[78] += -2 * dqdci; /* dwdot[HO2]/d[HO2] */
  J[79] += dqdci;      /* dwdot[H2O2]/d[HO2] */
  /* d()/d[H2O2] */
  dqdci = -k_r * sc[5];
  J[89] += dqdci;      /* dwdot[O2]/d[H2O2] */
  J[90] += -2 * dqdci; /* dwdot[HO2]/d[H2O2] */
  J[91] += dqdci;      /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[137] += dqdT;      /* dwdot[O2]/dT */
  J[138] += -2 * dqdT; /* dwdot[HO2]/dT */
  J[139] += dqdT;      /* dwdot[H2O2]/dT */

  /*reaction 24: H2O2 + H <=> H2O + OH */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[0] * sc[7];
  k_f = prefactor_units[23] * fwd_A[23] *
        exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
  dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
  /* reverse */
  phi_r = sc[3] * sc[4];
  Kc = exp(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[3] + h_RT[4]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[3] += q; /* OH */
  wdot[4] += q; /* H2O */
  wdot[7] -= q; /* H2O2 */
  /* d()/d[H] */
  dqdci = +k_f * sc[7];
  J[0] -= dqdci; /* dwdot[H]/d[H] */
  J[3] += dqdci; /* dwdot[OH]/d[H] */
  J[4] += dqdci; /* dwdot[H2O]/d[H] */
  J[7] -= dqdci; /* dwdot[H2O2]/d[H] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[4];
  J[36] -= dqdci; /* dwdot[H]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  J[40] += dqdci; /* dwdot[H2O]/d[OH] */
  J[43] -= dqdci; /* dwdot[H2O2]/d[OH] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[3];
  J[48] -= dqdci; /* dwdot[H]/d[H2O] */
  J[51] += dqdci; /* dwdot[OH]/d[H2O] */
  J[52] += dqdci; /* dwdot[H2O]/d[H2O] */
  J[55] -= dqdci; /* dwdot[H2O2]/d[H2O] */
  /* d()/d[H2O2] */
  dqdci = +k_f * sc[0];
  J[84] -= dqdci; /* dwdot[H]/d[H2O2] */
  J[87] += dqdci; /* dwdot[OH]/d[H2O2] */
  J[88] += dqdci; /* dwdot[H2O]/d[H2O2] */
  J[91] -= dqdci; /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */
  J[136] += dqdT; /* dwdot[H2O]/dT */
  J[139] -= dqdT; /* dwdot[H2O2]/dT */

  /*reaction 25: H2O2 + H <=> HO2 + H2 */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[0] * sc[7];
  k_f = prefactor_units[24] * fwd_A[24] *
        exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
  dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
  /* reverse */
  phi_r = sc[1] * sc[6];
  Kc = exp(g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[1] += q; /* H2 */
  wdot[6] += q; /* HO2 */
  wdot[7] -= q; /* H2O2 */
  /* d()/d[H] */
  dqdci = +k_f * sc[7];
  J[0] -= dqdci; /* dwdot[H]/d[H] */
  J[1] += dqdci; /* dwdot[H2]/d[H] */
  J[6] += dqdci; /* dwdot[HO2]/d[H] */
  J[7] -= dqdci; /* dwdot[H2O2]/d[H] */
  /* d()/d[H2] */
  dqdci = -k_r * sc[6];
  J[12] -= dqdci; /* dwdot[H]/d[H2] */
  J[13] += dqdci; /* dwdot[H2]/d[H2] */
  J[18] += dqdci; /* dwdot[HO2]/d[H2] */
  J[19] -= dqdci; /* dwdot[H2O2]/d[H2] */
  /* d()/d[HO2] */
  dqdci = -k_r * sc[1];
  J[72] -= dqdci; /* dwdot[H]/d[HO2] */
  J[73] += dqdci; /* dwdot[H2]/d[HO2] */
  J[78] += dqdci; /* dwdot[HO2]/d[HO2] */
  J[79] -= dqdci; /* dwdot[H2O2]/d[HO2] */
  /* d()/d[H2O2] */
  dqdci = +k_f * sc[0];
  J[84] -= dqdci; /* dwdot[H]/d[H2O2] */
  J[85] += dqdci; /* dwdot[H2]/d[H2O2] */
  J[90] += dqdci; /* dwdot[HO2]/d[H2O2] */
  J[91] -= dqdci; /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[133] += dqdT; /* dwdot[H2]/dT */
  J[138] += dqdT; /* dwdot[HO2]/dT */
  J[139] -= dqdT; /* dwdot[H2O2]/dT */

  /*reaction 26: H2O2 + O <=> OH + HO2 */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[2] * sc[7];
  k_f = prefactor_units[25] * fwd_A[25] *
        exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
  dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
  /* reverse */
  phi_r = sc[3] * sc[6];
  Kc = exp(g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[3] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] -= q; /* O */
  wdot[3] += q; /* OH */
  wdot[6] += q; /* HO2 */
  wdot[7] -= q; /* H2O2 */
  /* d()/d[O] */
  dqdci = +k_f * sc[7];
  J[26] -= dqdci; /* dwdot[O]/d[O] */
  J[27] += dqdci; /* dwdot[OH]/d[O] */
  J[30] += dqdci; /* dwdot[HO2]/d[O] */
  J[31] -= dqdci; /* dwdot[H2O2]/d[O] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[6];
  J[38] -= dqdci; /* dwdot[O]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  J[42] += dqdci; /* dwdot[HO2]/d[OH] */
  J[43] -= dqdci; /* dwdot[H2O2]/d[OH] */
  /* d()/d[HO2] */
  dqdci = -k_r * sc[3];
  J[74] -= dqdci; /* dwdot[O]/d[HO2] */
  J[75] += dqdci; /* dwdot[OH]/d[HO2] */
  J[78] += dqdci; /* dwdot[HO2]/d[HO2] */
  J[79] -= dqdci; /* dwdot[H2O2]/d[HO2] */
  /* d()/d[H2O2] */
  dqdci = +k_f * sc[2];
  J[86] -= dqdci; /* dwdot[O]/d[H2O2] */
  J[87] += dqdci; /* dwdot[OH]/d[H2O2] */
  J[90] += dqdci; /* dwdot[HO2]/d[H2O2] */
  J[91] -= dqdci; /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[134] -= dqdT; /* dwdot[O]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */
  J[138] += dqdT; /* dwdot[HO2]/dT */
  J[139] -= dqdT; /* dwdot[H2O2]/dT */

  /*reaction 27: H2O2 + OH <=> HO2 + H2O */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[3] * sc[7];
  k_f = prefactor_units[26] * fwd_A[26] *
        exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
  dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
  /* reverse */
  phi_r = sc[4] * sc[6];
  Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[3] -= q; /* OH */
  wdot[4] += q; /* H2O */
  wdot[6] += q; /* HO2 */
  wdot[7] -= q; /* H2O2 */
  /* d()/d[OH] */
  dqdci = +k_f * sc[7];
  J[39] -= dqdci; /* dwdot[OH]/d[OH] */
  J[40] += dqdci; /* dwdot[H2O]/d[OH] */
  J[42] += dqdci; /* dwdot[HO2]/d[OH] */
  J[43] -= dqdci; /* dwdot[H2O2]/d[OH] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[6];
  J[51] -= dqdci; /* dwdot[OH]/d[H2O] */
  J[52] += dqdci; /* dwdot[H2O]/d[H2O] */
  J[54] += dqdci; /* dwdot[HO2]/d[H2O] */
  J[55] -= dqdci; /* dwdot[H2O2]/d[H2O] */
  /* d()/d[HO2] */
  dqdci = -k_r * sc[4];
  J[75] -= dqdci; /* dwdot[OH]/d[HO2] */
  J[76] += dqdci; /* dwdot[H2O]/d[HO2] */
  J[78] += dqdci; /* dwdot[HO2]/d[HO2] */
  J[79] -= dqdci; /* dwdot[H2O2]/d[HO2] */
  /* d()/d[H2O2] */
  dqdci = +k_f * sc[3];
  J[87] -= dqdci; /* dwdot[OH]/d[H2O2] */
  J[88] += dqdci; /* dwdot[H2O]/d[H2O2] */
  J[90] += dqdci; /* dwdot[HO2]/d[H2O2] */
  J[91] -= dqdci; /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[135] -= dqdT; /* dwdot[OH]/dT */
  J[136] += dqdT; /* dwdot[H2O]/dT */
  J[138] += dqdT; /* dwdot[HO2]/dT */
  J[139] -= dqdT; /* dwdot[H2O2]/dT */

  /*reaction 28: H2O2 + OH <=> HO2 + H2O */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[3] * sc[7];
  k_f = prefactor_units[27] * fwd_A[27] *
        exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
  dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
  /* reverse */
  phi_r = sc[4] * sc[6];
  Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[3] -= q; /* OH */
  wdot[4] += q; /* H2O */
  wdot[6] += q; /* HO2 */
  wdot[7] -= q; /* H2O2 */
  /* d()/d[OH] */
  dqdci = +k_f * sc[7];
  J[39] -= dqdci; /* dwdot[OH]/d[OH] */
  J[40] += dqdci; /* dwdot[H2O]/d[OH] */
  J[42] += dqdci; /* dwdot[HO2]/d[OH] */
  J[43] -= dqdci; /* dwdot[H2O2]/d[OH] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[6];
  J[51] -= dqdci; /* dwdot[OH]/d[H2O] */
  J[52] += dqdci; /* dwdot[H2O]/d[H2O] */
  J[54] += dqdci; /* dwdot[HO2]/d[H2O] */
  J[55] -= dqdci; /* dwdot[H2O2]/d[H2O] */
  /* d()/d[HO2] */
  dqdci = -k_r * sc[4];
  J[75] -= dqdci; /* dwdot[OH]/d[HO2] */
  J[76] += dqdci; /* dwdot[H2O]/d[HO2] */
  J[78] += dqdci; /* dwdot[HO2]/d[HO2] */
  J[79] -= dqdci; /* dwdot[H2O2]/d[HO2] */
  /* d()/d[H2O2] */
  dqdci = +k_f * sc[3];
  J[87] -= dqdci; /* dwdot[OH]/d[H2O2] */
  J[88] += dqdci; /* dwdot[H2O]/d[H2O2] */
  J[90] += dqdci; /* dwdot[HO2]/d[H2O2] */
  J[91] -= dqdci; /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[135] -= dqdT; /* dwdot[OH]/dT */
  J[136] += dqdT; /* dwdot[H2O]/dT */
  J[138] += dqdT; /* dwdot[HO2]/dT */
  J[139] -= dqdT; /* dwdot[H2O2]/dT */

  /*reaction 29: HO2 + H <=> O + H2O */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[0] * sc[6];
  k_f = prefactor_units[28] * fwd_A[28] *
        exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
  dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
  /* reverse */
  phi_r = sc[2] * sc[4];
  Kc = exp(g_RT[0] - g_RT[2] - g_RT[4] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[2] + h_RT[4]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[2] += q; /* O */
  wdot[4] += q; /* H2O */
  wdot[6] -= q; /* HO2 */
  /* d()/d[H] */
  dqdci = +k_f * sc[6];
  J[0] -= dqdci; /* dwdot[H]/d[H] */
  J[2] += dqdci; /* dwdot[O]/d[H] */
  J[4] += dqdci; /* dwdot[H2O]/d[H] */
  J[6] -= dqdci; /* dwdot[HO2]/d[H] */
  /* d()/d[O] */
  dqdci = -k_r * sc[4];
  J[24] -= dqdci; /* dwdot[H]/d[O] */
  J[26] += dqdci; /* dwdot[O]/d[O] */
  J[28] += dqdci; /* dwdot[H2O]/d[O] */
  J[30] -= dqdci; /* dwdot[HO2]/d[O] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[2];
  J[48] -= dqdci; /* dwdot[H]/d[H2O] */
  J[50] += dqdci; /* dwdot[O]/d[H2O] */
  J[52] += dqdci; /* dwdot[H2O]/d[H2O] */
  J[54] -= dqdci; /* dwdot[HO2]/d[H2O] */
  /* d()/d[HO2] */
  dqdci = +k_f * sc[0];
  J[72] -= dqdci; /* dwdot[H]/d[HO2] */
  J[74] += dqdci; /* dwdot[O]/d[HO2] */
  J[76] += dqdci; /* dwdot[H2O]/d[HO2] */
  J[78] -= dqdci; /* dwdot[HO2]/d[HO2] */
  /* d()/dT */
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[134] += dqdT; /* dwdot[O]/dT */
  J[136] += dqdT; /* dwdot[H2O]/dT */
  J[138] -= dqdT; /* dwdot[HO2]/dT */

  double c_R[11], dcRdT[11], e_RT[11];
  double* eh_RT;
  if (consP) {
    cp_R(c_R, tc);
    dcvpRdT(dcRdT, tc);
    eh_RT = &h_RT[0];
  } else {
    cv_R(c_R, tc);
    dcvpRdT(dcRdT, tc);
    speciesInternalEnergy(e_RT, tc);
    eh_RT = &e_RT[0];
  }

  double cmix = 0.0, ehmix = 0.0, dcmixdT = 0.0, dehmixdT = 0.0;
  for (int k = 0; k < 11; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[132 + k];
  }

  double cmixinv = 1.0 / cmix;
  double tmp1 = ehmix * cmixinv;
  double tmp3 = cmixinv * T;
  double tmp2 = tmp1 * tmp3;
  double dehmixdc;
  /* dTdot/d[X] */
  for (int k = 0; k < 11; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 11; ++m) {
      dehmixdc += eh_RT[m] * J[k * 12 + m];
    }
    J[k * 12 + 11] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  /* dTdot/dT */
  J[143] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
}

/*compute an approx to the reaction Jacobian */
void
aJacobian_precond(double* J, double* sc, double T, int HP)
{
  for (int i = 0; i < 144; i++) {
    J[i] = 0.0;
  }

  double wdot[11];
  for (int k = 0; k < 11; k++) {
    wdot[k] = 0.0;
  }

  double tc[] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; /*temperature cache */
  double invT = 1.0 / tc[1];
  double invT2 = invT * invT;

  /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
  double refC = 101325 / 8.31451 / T;
  double refCinv = 1.0 / refC;

  /*compute the mixture concentration */
  double mixture = 0.0;
  for (int k = 0; k < 11; ++k) {
    mixture += sc[k];
  }

  /*compute the Gibbs free energy */
  double g_RT[11];
  gibbs(g_RT, tc);

  /*compute the species enthalpy */
  double h_RT[11];
  speciesEnthalpy(h_RT, tc);

  double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  double dqdci, dcdc_fac, dqdc[11];
  double Pr, fPr, F, k_0, logPr;
  double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  double Fcent1, Fcent2, Fcent3, Fcent;
  double dlogFdc, dlogFdn, dlogFdcn_fac;
  double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const double ln10 = log(10.0);
  const double log10e = 1.0 / log(10.0);
  /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
  /*a pressure-fall-off reaction */
  /* also 3-body */
  /* 3-body correction factor */
  alpha = mixture + (TB[0][0] - 1) * sc[1] + (TB[0][1] - 1) * sc[4] +
          (TB[0][2] - 1) * sc[5] + (TB[0][3] - 1) * sc[9] +
          (TB[0][4] - 1) * sc[10];
  /* forward */
  phi_f = sc[0] * sc[5];
  k_f = prefactor_units[0] * fwd_A[0] *
        exp(fwd_beta[0] * tc[0] - activation_units[0] * fwd_Ea[0] * invT);
  dlnkfdT = fwd_beta[0] * invT + activation_units[0] * fwd_Ea[0] * invT2;
  /* pressure-fall-off */
  k_0 = low_A[0] *
        exp(low_beta[0] * tc[0] - activation_units[0] * low_Ea[0] * invT);
  Pr = phase_units[0] * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = low_beta[0] * invT + activation_units[0] * low_Ea[0] * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  /* Troe form */
  logPr = log10(Pr);
  Fcent1 =
    (fabs(troe_Tsss[0]) > 1.e-100 ? (1. - troe_a[0]) * exp(-T / troe_Tsss[0])
                                  : 0.);
  Fcent2 = (fabs(troe_Ts[0]) > 1.e-100 ? troe_a[0] * exp(-T / troe_Ts[0]) : 0.);
  Fcent3 = (troe_len[0] == 4 ? exp(-troe_Tss[0] * invT) : 0.);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = pow(10.0, logFcent * troe);
  dlogFcentdT = log10e / Fcent *
                ((fabs(troe_Tsss[0]) > 1.e-100 ? -Fcent1 / troe_Tsss[0] : 0.) +
                 (fabs(troe_Ts[0]) > 1.e-100 ? -Fcent2 / troe_Ts[0] : 0.) +
                 (troe_len[0] == 4 ? Fcent3 * troe_Tss[0] * invT2 : 0.));
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  /* reverse */
  phi_r = sc[6];
  Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[6]) + 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[5] -= q; /* O2 */
  wdot[6] += q; /* HO2 */
  /* for convenience */
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = dcdc_fac + k_f * sc[5];
  dqdc[1] = TB[0][0] * dcdc_fac;
  dqdc[2] = dcdc_fac;
  dqdc[3] = dcdc_fac;
  dqdc[4] = TB[0][1] * dcdc_fac;
  dqdc[5] = TB[0][2] * dcdc_fac + k_f * sc[0];
  dqdc[6] = dcdc_fac - k_r;
  dqdc[7] = dcdc_fac;
  dqdc[8] = dcdc_fac;
  dqdc[9] = TB[0][3] * dcdc_fac;
  dqdc[10] = TB[0][4] * dcdc_fac;
  for (int k = 0; k < 11; k++) {
    J[12 * k + 0] -= dqdc[k];
    J[12 * k + 5] -= dqdc[k];
    J[12 * k + 6] += dqdc[k];
  }
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[137] -= dqdT; /* dwdot[O2]/dT */
  J[138] += dqdT; /* dwdot[HO2]/dT */

  /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
  /*a pressure-fall-off reaction */
  /* also 3-body */
  /* 3-body correction factor */
  alpha = mixture + (TB[1][0] - 1) * sc[4] + (TB[1][1] - 1) * sc[8] +
          (TB[1][2] - 1) * sc[5] + (TB[1][3] - 1) * sc[10] +
          (TB[1][4] - 1) * sc[7] + (TB[1][5] - 1) * sc[1];
  /* forward */
  phi_f = sc[7];
  k_f = prefactor_units[1] * fwd_A[1] *
        exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
  dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
  /* pressure-fall-off */
  k_0 = low_A[1] *
        exp(low_beta[1] * tc[0] - activation_units[1] * low_Ea[1] * invT);
  Pr = phase_units[1] * alpha / k_f * k_0;
  fPr = Pr / (1.0 + Pr);
  dlnk0dT = low_beta[1] * invT + activation_units[1] * low_Ea[1] * invT2;
  dlogPrdT = log10e * (dlnk0dT - dlnkfdT);
  dlogfPrdT = dlogPrdT / (1.0 + Pr);
  /* Troe form */
  logPr = log10(Pr);
  Fcent1 =
    (fabs(troe_Tsss[1]) > 1.e-100 ? (1. - troe_a[1]) * exp(-T / troe_Tsss[1])
                                  : 0.);
  Fcent2 = (fabs(troe_Ts[1]) > 1.e-100 ? troe_a[1] * exp(-T / troe_Ts[1]) : 0.);
  Fcent3 = (troe_len[1] == 4 ? exp(-troe_Tss[1] * invT) : 0.);
  Fcent = Fcent1 + Fcent2 + Fcent3;
  logFcent = log10(Fcent);
  troe_c = -.4 - .67 * logFcent;
  troe_n = .75 - 1.27 * logFcent;
  troePr_den = 1.0 / (troe_n - .14 * (troe_c + logPr));
  troePr = (troe_c + logPr) * troePr_den;
  troe = 1.0 / (1.0 + troePr * troePr);
  F = pow(10.0, logFcent * troe);
  dlogFcentdT = log10e / Fcent *
                ((fabs(troe_Tsss[1]) > 1.e-100 ? -Fcent1 / troe_Tsss[1] : 0.) +
                 (fabs(troe_Ts[1]) > 1.e-100 ? -Fcent2 / troe_Ts[1] : 0.) +
                 (troe_len[1] == 4 ? Fcent3 * troe_Tss[1] * invT2 : 0.));
  dlogFdcn_fac = 2.0 * logFcent * troe * troe * troePr * troePr_den;
  dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
  dlogFdn = dlogFdcn_fac * troePr;
  dlogFdlogPr = dlogFdc;
  dlogFdT = dlogFcentdT * (troe - 0.67 * dlogFdc - 1.27 * dlogFdn) +
            dlogFdlogPr * dlogPrdT;
  /* reverse */
  phi_r = sc[3] * sc[3];
  Kc = refC * exp(-g_RT[3] - g_RT[3] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[7]) + (2 * h_RT[3]) - 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  Corr = fPr * F;
  q = Corr * q_nocor;
  dlnCorrdT = ln10 * (dlogfPrdT + dlogFdT);
  dqdT = Corr * (dlnkfdT * k_f * phi_f - dkrdT * phi_r) + dlnCorrdT * q;
  /* update wdot */
  wdot[3] += 2 * q; /* OH */
  wdot[7] -= q;     /* H2O2 */
  /* for convenience */
  k_f *= Corr;
  k_r *= Corr;
  dcdc_fac = 0.0;
  dqdc[0] = dcdc_fac;
  dqdc[1] = TB[1][5] * dcdc_fac;
  dqdc[2] = dcdc_fac;
  dqdc[3] = dcdc_fac - k_r * 2 * sc[3];
  dqdc[4] = TB[1][0] * dcdc_fac;
  dqdc[5] = TB[1][2] * dcdc_fac;
  dqdc[6] = dcdc_fac;
  dqdc[7] = TB[1][4] * dcdc_fac + k_f;
  dqdc[8] = TB[1][1] * dcdc_fac;
  dqdc[9] = dcdc_fac;
  dqdc[10] = TB[1][3] * dcdc_fac;
  for (int k = 0; k < 11; k++) {
    J[12 * k + 3] += 2 * dqdc[k];
    J[12 * k + 7] -= dqdc[k];
  }
  J[135] += 2 * dqdT; /* dwdot[OH]/dT */
  J[139] -= dqdT;     /* dwdot[H2O2]/dT */

  /*reaction 3: H2 + M <=> H + H + M */
  /*a third-body and non-pressure-fall-off reaction */
  /* 3-body correction factor */
  alpha = mixture + (TB[2][0] - 1) * sc[1] + (TB[2][1] - 1) * sc[4] +
          (TB[2][2] - 1) * sc[9] + (TB[2][3] - 1) * sc[10];
  /* forward */
  phi_f = sc[1];
  k_f = prefactor_units[2] * fwd_A[2] *
        exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
  dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[0];
  Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1]) + (2 * h_RT[0]) - 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += 2 * q; /* H */
  wdot[1] -= q;     /* H2 */
  /* for convenience */
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = q_nocor - k_r * 2 * sc[0];
  dqdc[1] = TB[2][0] * q_nocor + k_f;
  dqdc[2] = q_nocor;
  dqdc[3] = q_nocor;
  dqdc[4] = TB[2][1] * q_nocor;
  dqdc[5] = q_nocor;
  dqdc[6] = q_nocor;
  dqdc[7] = q_nocor;
  dqdc[8] = q_nocor;
  dqdc[9] = TB[2][2] * q_nocor;
  dqdc[10] = TB[2][3] * q_nocor;
  for (int k = 0; k < 11; k++) {
    J[12 * k + 0] += 2 * dqdc[k];
    J[12 * k + 1] -= dqdc[k];
  }
  J[132] += 2 * dqdT; /* dwdot[H]/dT */
  J[133] -= dqdT;     /* dwdot[H2]/dT */

  /*reaction 4: O + O + M <=> O2 + M */
  /*a third-body and non-pressure-fall-off reaction */
  /* 3-body correction factor */
  alpha = mixture + (TB[3][0] - 1) * sc[1] + (TB[3][1] - 1) * sc[4] +
          (TB[3][2] - 1) * sc[9] + (TB[3][3] - 1) * sc[10];
  /* forward */
  phi_f = sc[2] * sc[2];
  k_f = prefactor_units[3] * fwd_A[3] *
        exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
  dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
  /* reverse */
  phi_r = sc[5];
  Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[2]) + (h_RT[5]) + 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] -= 2 * q; /* O */
  wdot[5] += q;     /* O2 */
  /* for convenience */
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = q_nocor;
  dqdc[1] = TB[3][0] * q_nocor;
  dqdc[2] = q_nocor + k_f * 2 * sc[2];
  dqdc[3] = q_nocor;
  dqdc[4] = TB[3][1] * q_nocor;
  dqdc[5] = q_nocor - k_r;
  dqdc[6] = q_nocor;
  dqdc[7] = q_nocor;
  dqdc[8] = q_nocor;
  dqdc[9] = TB[3][2] * q_nocor;
  dqdc[10] = TB[3][3] * q_nocor;
  for (int k = 0; k < 11; k++) {
    J[12 * k + 2] += -2 * dqdc[k];
    J[12 * k + 5] += dqdc[k];
  }
  J[134] += -2 * dqdT; /* dwdot[O]/dT */
  J[137] += dqdT;      /* dwdot[O2]/dT */

  /*reaction 5: O + H + M <=> OH + M */
  /*a third-body and non-pressure-fall-off reaction */
  /* 3-body correction factor */
  alpha = mixture + (TB[4][0] - 1) * sc[1] + (TB[4][1] - 1) * sc[4] +
          (TB[4][2] - 1) * sc[9] + (TB[4][3] - 1) * sc[10];
  /* forward */
  phi_f = sc[0] * sc[2];
  k_f = prefactor_units[4] * fwd_A[4] *
        exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
  dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
  /* reverse */
  phi_r = sc[3];
  Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[3]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3]) + 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[2] -= q; /* O */
  wdot[3] += q; /* OH */
  /* for convenience */
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = q_nocor + k_f * sc[2];
  dqdc[1] = TB[4][0] * q_nocor;
  dqdc[2] = q_nocor + k_f * sc[0];
  dqdc[3] = q_nocor - k_r;
  dqdc[4] = TB[4][1] * q_nocor;
  dqdc[5] = q_nocor;
  dqdc[6] = q_nocor;
  dqdc[7] = q_nocor;
  dqdc[8] = q_nocor;
  dqdc[9] = TB[4][2] * q_nocor;
  dqdc[10] = TB[4][3] * q_nocor;
  for (int k = 0; k < 11; k++) {
    J[12 * k + 0] -= dqdc[k];
    J[12 * k + 2] -= dqdc[k];
    J[12 * k + 3] += dqdc[k];
  }
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[134] -= dqdT; /* dwdot[O]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */

  /*reaction 6: H2O + M <=> H + OH + M */
  /*a third-body and non-pressure-fall-off reaction */
  /* 3-body correction factor */
  alpha = mixture + (TB[5][0] - 1) * sc[1] + (TB[5][1] - 1) * sc[4] +
          (TB[5][2] - 1) * sc[10] + (TB[5][3] - 1) * sc[8] +
          (TB[5][4] - 1) * sc[5];
  /* forward */
  phi_f = sc[4];
  k_f = prefactor_units[5] * fwd_A[5] *
        exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
  dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[3];
  Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[4]) + (h_RT[0] + h_RT[3]) - 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += q; /* H */
  wdot[3] += q; /* OH */
  wdot[4] -= q; /* H2O */
  /* for convenience */
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = q_nocor - k_r * sc[3];
  dqdc[1] = TB[5][0] * q_nocor;
  dqdc[2] = q_nocor;
  dqdc[3] = q_nocor - k_r * sc[0];
  dqdc[4] = TB[5][1] * q_nocor + k_f;
  dqdc[5] = TB[5][4] * q_nocor;
  dqdc[6] = q_nocor;
  dqdc[7] = q_nocor;
  dqdc[8] = TB[5][3] * q_nocor;
  dqdc[9] = q_nocor;
  dqdc[10] = TB[5][2] * q_nocor;
  for (int k = 0; k < 11; k++) {
    J[12 * k + 0] += dqdc[k];
    J[12 * k + 3] += dqdc[k];
    J[12 * k + 4] -= dqdc[k];
  }
  J[132] += dqdT; /* dwdot[H]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */
  J[136] -= dqdT; /* dwdot[H2O]/dT */

  /*reaction 7: O + OH + M <=> HO2 + M */
  /*a third-body and non-pressure-fall-off reaction */
  /* 3-body correction factor */
  alpha = mixture + (TB[6][0] - 1) * sc[1] + (TB[6][1] - 1) * sc[4] +
          (TB[6][2] - 1) * sc[9] + (TB[6][3] - 1) * sc[10];
  /* forward */
  phi_f = sc[2] * sc[3];
  k_f = prefactor_units[6] * fwd_A[6] *
        exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
  dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
  /* reverse */
  phi_r = sc[6];
  Kc = refCinv * exp(g_RT[2] + g_RT[3] - g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[6]) + 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q_nocor = k_f * phi_f - k_r * phi_r;
  q = alpha * q_nocor;
  dqdT = alpha * (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] -= q; /* O */
  wdot[3] -= q; /* OH */
  wdot[6] += q; /* HO2 */
  /* for convenience */
  k_f *= alpha;
  k_r *= alpha;
  dqdc[0] = q_nocor;
  dqdc[1] = TB[6][0] * q_nocor;
  dqdc[2] = q_nocor + k_f * sc[3];
  dqdc[3] = q_nocor + k_f * sc[2];
  dqdc[4] = TB[6][1] * q_nocor;
  dqdc[5] = q_nocor;
  dqdc[6] = q_nocor - k_r;
  dqdc[7] = q_nocor;
  dqdc[8] = q_nocor;
  dqdc[9] = TB[6][2] * q_nocor;
  dqdc[10] = TB[6][3] * q_nocor;
  for (int k = 0; k < 11; k++) {
    J[12 * k + 2] -= dqdc[k];
    J[12 * k + 3] -= dqdc[k];
    J[12 * k + 6] += dqdc[k];
  }
  J[134] -= dqdT; /* dwdot[O]/dT */
  J[135] -= dqdT; /* dwdot[OH]/dT */
  J[138] += dqdT; /* dwdot[HO2]/dT */

  /*reaction 8: H + O2 <=> O + OH */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[0] * sc[5];
  k_f = prefactor_units[7] * fwd_A[7] *
        exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
  dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
  /* reverse */
  phi_r = sc[2] * sc[3];
  Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[2] += q; /* O */
  wdot[3] += q; /* OH */
  wdot[5] -= q; /* O2 */
  /* d()/d[H] */
  dqdci = +k_f * sc[5];
  J[0] -= dqdci; /* dwdot[H]/d[H] */
  J[2] += dqdci; /* dwdot[O]/d[H] */
  J[3] += dqdci; /* dwdot[OH]/d[H] */
  J[5] -= dqdci; /* dwdot[O2]/d[H] */
  /* d()/d[O] */
  dqdci = -k_r * sc[3];
  J[24] -= dqdci; /* dwdot[H]/d[O] */
  J[26] += dqdci; /* dwdot[O]/d[O] */
  J[27] += dqdci; /* dwdot[OH]/d[O] */
  J[29] -= dqdci; /* dwdot[O2]/d[O] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[2];
  J[36] -= dqdci; /* dwdot[H]/d[OH] */
  J[38] += dqdci; /* dwdot[O]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  J[41] -= dqdci; /* dwdot[O2]/d[OH] */
  /* d()/d[O2] */
  dqdci = +k_f * sc[0];
  J[60] -= dqdci; /* dwdot[H]/d[O2] */
  J[62] += dqdci; /* dwdot[O]/d[O2] */
  J[63] += dqdci; /* dwdot[OH]/d[O2] */
  J[65] -= dqdci; /* dwdot[O2]/d[O2] */
  /* d()/dT */
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[134] += dqdT; /* dwdot[O]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */
  J[137] -= dqdT; /* dwdot[O2]/dT */

  /*reaction 9: O + H2 <=> H + OH */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[1] * sc[2];
  k_f = prefactor_units[8] * fwd_A[8] *
        exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
  dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[3];
  Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += q; /* H */
  wdot[1] -= q; /* H2 */
  wdot[2] -= q; /* O */
  wdot[3] += q; /* OH */
  /* d()/d[H] */
  dqdci = -k_r * sc[3];
  J[0] += dqdci; /* dwdot[H]/d[H] */
  J[1] -= dqdci; /* dwdot[H2]/d[H] */
  J[2] -= dqdci; /* dwdot[O]/d[H] */
  J[3] += dqdci; /* dwdot[OH]/d[H] */
  /* d()/d[H2] */
  dqdci = +k_f * sc[2];
  J[12] += dqdci; /* dwdot[H]/d[H2] */
  J[13] -= dqdci; /* dwdot[H2]/d[H2] */
  J[14] -= dqdci; /* dwdot[O]/d[H2] */
  J[15] += dqdci; /* dwdot[OH]/d[H2] */
  /* d()/d[O] */
  dqdci = +k_f * sc[1];
  J[24] += dqdci; /* dwdot[H]/d[O] */
  J[25] -= dqdci; /* dwdot[H2]/d[O] */
  J[26] -= dqdci; /* dwdot[O]/d[O] */
  J[27] += dqdci; /* dwdot[OH]/d[O] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[0];
  J[36] += dqdci; /* dwdot[H]/d[OH] */
  J[37] -= dqdci; /* dwdot[H2]/d[OH] */
  J[38] -= dqdci; /* dwdot[O]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  /* d()/dT */
  J[132] += dqdT; /* dwdot[H]/dT */
  J[133] -= dqdT; /* dwdot[H2]/dT */
  J[134] -= dqdT; /* dwdot[O]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */

  /*reaction 10: O + H2 <=> H + OH */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[1] * sc[2];
  k_f = prefactor_units[9] * fwd_A[9] *
        exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
  dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[3];
  Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += q; /* H */
  wdot[1] -= q; /* H2 */
  wdot[2] -= q; /* O */
  wdot[3] += q; /* OH */
  /* d()/d[H] */
  dqdci = -k_r * sc[3];
  J[0] += dqdci; /* dwdot[H]/d[H] */
  J[1] -= dqdci; /* dwdot[H2]/d[H] */
  J[2] -= dqdci; /* dwdot[O]/d[H] */
  J[3] += dqdci; /* dwdot[OH]/d[H] */
  /* d()/d[H2] */
  dqdci = +k_f * sc[2];
  J[12] += dqdci; /* dwdot[H]/d[H2] */
  J[13] -= dqdci; /* dwdot[H2]/d[H2] */
  J[14] -= dqdci; /* dwdot[O]/d[H2] */
  J[15] += dqdci; /* dwdot[OH]/d[H2] */
  /* d()/d[O] */
  dqdci = +k_f * sc[1];
  J[24] += dqdci; /* dwdot[H]/d[O] */
  J[25] -= dqdci; /* dwdot[H2]/d[O] */
  J[26] -= dqdci; /* dwdot[O]/d[O] */
  J[27] += dqdci; /* dwdot[OH]/d[O] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[0];
  J[36] += dqdci; /* dwdot[H]/d[OH] */
  J[37] -= dqdci; /* dwdot[H2]/d[OH] */
  J[38] -= dqdci; /* dwdot[O]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  /* d()/dT */
  J[132] += dqdT; /* dwdot[H]/dT */
  J[133] -= dqdT; /* dwdot[H2]/dT */
  J[134] -= dqdT; /* dwdot[O]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */

  /*reaction 11: H2 + OH <=> H2O + H */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[1] * sc[3];
  k_f = prefactor_units[10] * fwd_A[10] *
        exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
  dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[4];
  Kc = exp(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[0] + h_RT[4]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += q; /* H */
  wdot[1] -= q; /* H2 */
  wdot[3] -= q; /* OH */
  wdot[4] += q; /* H2O */
  /* d()/d[H] */
  dqdci = -k_r * sc[4];
  J[0] += dqdci; /* dwdot[H]/d[H] */
  J[1] -= dqdci; /* dwdot[H2]/d[H] */
  J[3] -= dqdci; /* dwdot[OH]/d[H] */
  J[4] += dqdci; /* dwdot[H2O]/d[H] */
  /* d()/d[H2] */
  dqdci = +k_f * sc[3];
  J[12] += dqdci; /* dwdot[H]/d[H2] */
  J[13] -= dqdci; /* dwdot[H2]/d[H2] */
  J[15] -= dqdci; /* dwdot[OH]/d[H2] */
  J[16] += dqdci; /* dwdot[H2O]/d[H2] */
  /* d()/d[OH] */
  dqdci = +k_f * sc[1];
  J[36] += dqdci; /* dwdot[H]/d[OH] */
  J[37] -= dqdci; /* dwdot[H2]/d[OH] */
  J[39] -= dqdci; /* dwdot[OH]/d[OH] */
  J[40] += dqdci; /* dwdot[H2O]/d[OH] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[0];
  J[48] += dqdci; /* dwdot[H]/d[H2O] */
  J[49] -= dqdci; /* dwdot[H2]/d[H2O] */
  J[51] -= dqdci; /* dwdot[OH]/d[H2O] */
  J[52] += dqdci; /* dwdot[H2O]/d[H2O] */
  /* d()/dT */
  J[132] += dqdT; /* dwdot[H]/dT */
  J[133] -= dqdT; /* dwdot[H2]/dT */
  J[135] -= dqdT; /* dwdot[OH]/dT */
  J[136] += dqdT; /* dwdot[H2O]/dT */

  /*reaction 12: OH + OH <=> O + H2O */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[3] * sc[3];
  k_f = prefactor_units[11] * fwd_A[11] *
        exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
  dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
  /* reverse */
  phi_r = sc[2] * sc[4];
  Kc = exp(-g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[3]) + (h_RT[2] + h_RT[4]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] += q;     /* O */
  wdot[3] -= 2 * q; /* OH */
  wdot[4] += q;     /* H2O */
  /* d()/d[O] */
  dqdci = -k_r * sc[4];
  J[26] += dqdci;      /* dwdot[O]/d[O] */
  J[27] += -2 * dqdci; /* dwdot[OH]/d[O] */
  J[28] += dqdci;      /* dwdot[H2O]/d[O] */
  /* d()/d[OH] */
  dqdci = +k_f * 2 * sc[3];
  J[38] += dqdci;      /* dwdot[O]/d[OH] */
  J[39] += -2 * dqdci; /* dwdot[OH]/d[OH] */
  J[40] += dqdci;      /* dwdot[H2O]/d[OH] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[2];
  J[50] += dqdci;      /* dwdot[O]/d[H2O] */
  J[51] += -2 * dqdci; /* dwdot[OH]/d[H2O] */
  J[52] += dqdci;      /* dwdot[H2O]/d[H2O] */
  /* d()/dT */
  J[134] += dqdT;      /* dwdot[O]/dT */
  J[135] += -2 * dqdT; /* dwdot[OH]/dT */
  J[136] += dqdT;      /* dwdot[H2O]/dT */

  /*reaction 13: H2 + AR <=> H + H + AR */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[1] * sc[9];
  k_f = prefactor_units[12] * fwd_A[12] *
        exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
  dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[0] * sc[9];
  Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (2 * h_RT[0] + h_RT[9]) - 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += 2 * q; /* H */
  wdot[1] -= q;     /* H2 */
  /* d()/d[H] */
  dqdci = -k_r * 2 * sc[0] * sc[9];
  J[0] += 2 * dqdci; /* dwdot[H]/d[H] */
  J[1] -= dqdci;     /* dwdot[H2]/d[H] */
  /* d()/d[H2] */
  dqdci = +k_f * sc[9];
  J[12] += 2 * dqdci; /* dwdot[H]/d[H2] */
  J[13] -= dqdci;     /* dwdot[H2]/d[H2] */
  /* d()/d[AR] */
  dqdci = +k_f * sc[1] - k_r * sc[0] * sc[0];
  J[108] += 2 * dqdci; /* dwdot[H]/d[AR] */
  J[109] -= dqdci;     /* dwdot[H2]/d[AR] */
  /* d()/dT */
  J[132] += 2 * dqdT; /* dwdot[H]/dT */
  J[133] -= dqdT;     /* dwdot[H2]/dT */

  /*reaction 14: H2 + HE <=> H + H + HE */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[1] * sc[10];
  k_f = prefactor_units[13] * fwd_A[13] *
        exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
  dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[0] * sc[10];
  Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (2 * h_RT[0] + h_RT[10]) - 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += 2 * q; /* H */
  wdot[1] -= q;     /* H2 */
  /* d()/d[H] */
  dqdci = -k_r * 2 * sc[0] * sc[10];
  J[0] += 2 * dqdci; /* dwdot[H]/d[H] */
  J[1] -= dqdci;     /* dwdot[H2]/d[H] */
  /* d()/d[H2] */
  dqdci = +k_f * sc[10];
  J[12] += 2 * dqdci; /* dwdot[H]/d[H2] */
  J[13] -= dqdci;     /* dwdot[H2]/d[H2] */
  /* d()/d[HE] */
  dqdci = +k_f * sc[1] - k_r * sc[0] * sc[0];
  J[120] += 2 * dqdci; /* dwdot[H]/d[HE] */
  J[121] -= dqdci;     /* dwdot[H2]/d[HE] */
  /* d()/dT */
  J[132] += 2 * dqdT; /* dwdot[H]/dT */
  J[133] -= dqdT;     /* dwdot[H2]/dT */

  /*reaction 15: O + O + AR <=> O2 + AR */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[2] * sc[2] * sc[9];
  k_f = prefactor_units[14] * fwd_A[14] *
        exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
  dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
  /* reverse */
  phi_r = sc[5] * sc[9];
  Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[2] + h_RT[9]) + (h_RT[5] + h_RT[9]) + 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] -= 2 * q; /* O */
  wdot[5] += q;     /* O2 */
  /* d()/d[O] */
  dqdci = +k_f * 2 * sc[2] * sc[9];
  J[26] += -2 * dqdci; /* dwdot[O]/d[O] */
  J[29] += dqdci;      /* dwdot[O2]/d[O] */
  /* d()/d[O2] */
  dqdci = -k_r * sc[9];
  J[62] += -2 * dqdci; /* dwdot[O]/d[O2] */
  J[65] += dqdci;      /* dwdot[O2]/d[O2] */
  /* d()/d[AR] */
  dqdci = +k_f * sc[2] * sc[2] - k_r * sc[5];
  J[110] += -2 * dqdci; /* dwdot[O]/d[AR] */
  J[113] += dqdci;      /* dwdot[O2]/d[AR] */
  /* d()/dT */
  J[134] += -2 * dqdT; /* dwdot[O]/dT */
  J[137] += dqdT;      /* dwdot[O2]/dT */

  /*reaction 16: O + O + HE <=> O2 + HE */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[2] * sc[2] * sc[10];
  k_f = prefactor_units[15] * fwd_A[15] *
        exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
  dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
  /* reverse */
  phi_r = sc[5] * sc[10];
  Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[2] + h_RT[10]) + (h_RT[5] + h_RT[10]) + 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] -= 2 * q; /* O */
  wdot[5] += q;     /* O2 */
  /* d()/d[O] */
  dqdci = +k_f * 2 * sc[2] * sc[10];
  J[26] += -2 * dqdci; /* dwdot[O]/d[O] */
  J[29] += dqdci;      /* dwdot[O2]/d[O] */
  /* d()/d[O2] */
  dqdci = -k_r * sc[10];
  J[62] += -2 * dqdci; /* dwdot[O]/d[O2] */
  J[65] += dqdci;      /* dwdot[O2]/d[O2] */
  /* d()/d[HE] */
  dqdci = +k_f * sc[2] * sc[2] - k_r * sc[5];
  J[122] += -2 * dqdci; /* dwdot[O]/d[HE] */
  J[125] += dqdci;      /* dwdot[O2]/d[HE] */
  /* d()/dT */
  J[134] += -2 * dqdT; /* dwdot[O]/dT */
  J[137] += dqdT;      /* dwdot[O2]/dT */

  /*reaction 17: H2O + H2O <=> H + OH + H2O */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[4] * sc[4];
  k_f = prefactor_units[16] * fwd_A[16] *
        exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
  dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
  /* reverse */
  phi_r = sc[0] * sc[3] * sc[4];
  Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[4]) + (h_RT[0] + h_RT[3] + h_RT[4]) - 1);
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] += q; /* H */
  wdot[3] += q; /* OH */
  wdot[4] -= q; /* H2O */
  /* d()/d[H] */
  dqdci = -k_r * sc[3] * sc[4];
  J[0] += dqdci; /* dwdot[H]/d[H] */
  J[3] += dqdci; /* dwdot[OH]/d[H] */
  J[4] -= dqdci; /* dwdot[H2O]/d[H] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[0] * sc[4];
  J[36] += dqdci; /* dwdot[H]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  J[40] -= dqdci; /* dwdot[H2O]/d[OH] */
  /* d()/d[H2O] */
  dqdci = +k_f * 2 * sc[4] - k_r * sc[0] * sc[3];
  J[48] += dqdci; /* dwdot[H]/d[H2O] */
  J[51] += dqdci; /* dwdot[OH]/d[H2O] */
  J[52] -= dqdci; /* dwdot[H2O]/d[H2O] */
  /* d()/dT */
  J[132] += dqdT; /* dwdot[H]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */
  J[136] -= dqdT; /* dwdot[H2O]/dT */

  /*reaction 18: HO2 + H <=> H2 + O2 */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[0] * sc[6];
  k_f = prefactor_units[17] * fwd_A[17] *
        exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
  dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
  /* reverse */
  phi_r = sc[1] * sc[5];
  Kc = exp(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[1] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[1] += q; /* H2 */
  wdot[5] += q; /* O2 */
  wdot[6] -= q; /* HO2 */
  /* d()/d[H] */
  dqdci = +k_f * sc[6];
  J[0] -= dqdci; /* dwdot[H]/d[H] */
  J[1] += dqdci; /* dwdot[H2]/d[H] */
  J[5] += dqdci; /* dwdot[O2]/d[H] */
  J[6] -= dqdci; /* dwdot[HO2]/d[H] */
  /* d()/d[H2] */
  dqdci = -k_r * sc[5];
  J[12] -= dqdci; /* dwdot[H]/d[H2] */
  J[13] += dqdci; /* dwdot[H2]/d[H2] */
  J[17] += dqdci; /* dwdot[O2]/d[H2] */
  J[18] -= dqdci; /* dwdot[HO2]/d[H2] */
  /* d()/d[O2] */
  dqdci = -k_r * sc[1];
  J[60] -= dqdci; /* dwdot[H]/d[O2] */
  J[61] += dqdci; /* dwdot[H2]/d[O2] */
  J[65] += dqdci; /* dwdot[O2]/d[O2] */
  J[66] -= dqdci; /* dwdot[HO2]/d[O2] */
  /* d()/d[HO2] */
  dqdci = +k_f * sc[0];
  J[72] -= dqdci; /* dwdot[H]/d[HO2] */
  J[73] += dqdci; /* dwdot[H2]/d[HO2] */
  J[77] += dqdci; /* dwdot[O2]/d[HO2] */
  J[78] -= dqdci; /* dwdot[HO2]/d[HO2] */
  /* d()/dT */
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[133] += dqdT; /* dwdot[H2]/dT */
  J[137] += dqdT; /* dwdot[O2]/dT */
  J[138] -= dqdT; /* dwdot[HO2]/dT */

  /*reaction 19: HO2 + H <=> OH + OH */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[0] * sc[6];
  k_f = prefactor_units[18] * fwd_A[18] *
        exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
  dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
  /* reverse */
  phi_r = sc[3] * sc[3];
  Kc = exp(g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (2 * h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q;     /* H */
  wdot[3] += 2 * q; /* OH */
  wdot[6] -= q;     /* HO2 */
  /* d()/d[H] */
  dqdci = +k_f * sc[6];
  J[0] -= dqdci;     /* dwdot[H]/d[H] */
  J[3] += 2 * dqdci; /* dwdot[OH]/d[H] */
  J[6] -= dqdci;     /* dwdot[HO2]/d[H] */
  /* d()/d[OH] */
  dqdci = -k_r * 2 * sc[3];
  J[36] -= dqdci;     /* dwdot[H]/d[OH] */
  J[39] += 2 * dqdci; /* dwdot[OH]/d[OH] */
  J[42] -= dqdci;     /* dwdot[HO2]/d[OH] */
  /* d()/d[HO2] */
  dqdci = +k_f * sc[0];
  J[72] -= dqdci;     /* dwdot[H]/d[HO2] */
  J[75] += 2 * dqdci; /* dwdot[OH]/d[HO2] */
  J[78] -= dqdci;     /* dwdot[HO2]/d[HO2] */
  /* d()/dT */
  J[132] -= dqdT;     /* dwdot[H]/dT */
  J[135] += 2 * dqdT; /* dwdot[OH]/dT */
  J[138] -= dqdT;     /* dwdot[HO2]/dT */

  /*reaction 20: HO2 + O <=> O2 + OH */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[2] * sc[6];
  k_f = prefactor_units[19] * fwd_A[19] *
        exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
  dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
  /* reverse */
  phi_r = sc[3] * sc[5];
  Kc = exp(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] -= q; /* O */
  wdot[3] += q; /* OH */
  wdot[5] += q; /* O2 */
  wdot[6] -= q; /* HO2 */
  /* d()/d[O] */
  dqdci = +k_f * sc[6];
  J[26] -= dqdci; /* dwdot[O]/d[O] */
  J[27] += dqdci; /* dwdot[OH]/d[O] */
  J[29] += dqdci; /* dwdot[O2]/d[O] */
  J[30] -= dqdci; /* dwdot[HO2]/d[O] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[5];
  J[38] -= dqdci; /* dwdot[O]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  J[41] += dqdci; /* dwdot[O2]/d[OH] */
  J[42] -= dqdci; /* dwdot[HO2]/d[OH] */
  /* d()/d[O2] */
  dqdci = -k_r * sc[3];
  J[62] -= dqdci; /* dwdot[O]/d[O2] */
  J[63] += dqdci; /* dwdot[OH]/d[O2] */
  J[65] += dqdci; /* dwdot[O2]/d[O2] */
  J[66] -= dqdci; /* dwdot[HO2]/d[O2] */
  /* d()/d[HO2] */
  dqdci = +k_f * sc[2];
  J[74] -= dqdci; /* dwdot[O]/d[HO2] */
  J[75] += dqdci; /* dwdot[OH]/d[HO2] */
  J[77] += dqdci; /* dwdot[O2]/d[HO2] */
  J[78] -= dqdci; /* dwdot[HO2]/d[HO2] */
  /* d()/dT */
  J[134] -= dqdT; /* dwdot[O]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */
  J[137] += dqdT; /* dwdot[O2]/dT */
  J[138] -= dqdT; /* dwdot[HO2]/dT */

  /*reaction 21: HO2 + OH <=> H2O + O2 */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[3] * sc[6];
  k_f = prefactor_units[20] * fwd_A[20] *
        exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
  dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
  /* reverse */
  phi_r = sc[4] * sc[5];
  Kc = exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[4] + h_RT[5]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[3] -= q; /* OH */
  wdot[4] += q; /* H2O */
  wdot[5] += q; /* O2 */
  wdot[6] -= q; /* HO2 */
  /* d()/d[OH] */
  dqdci = +k_f * sc[6];
  J[39] -= dqdci; /* dwdot[OH]/d[OH] */
  J[40] += dqdci; /* dwdot[H2O]/d[OH] */
  J[41] += dqdci; /* dwdot[O2]/d[OH] */
  J[42] -= dqdci; /* dwdot[HO2]/d[OH] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[5];
  J[51] -= dqdci; /* dwdot[OH]/d[H2O] */
  J[52] += dqdci; /* dwdot[H2O]/d[H2O] */
  J[53] += dqdci; /* dwdot[O2]/d[H2O] */
  J[54] -= dqdci; /* dwdot[HO2]/d[H2O] */
  /* d()/d[O2] */
  dqdci = -k_r * sc[4];
  J[63] -= dqdci; /* dwdot[OH]/d[O2] */
  J[64] += dqdci; /* dwdot[H2O]/d[O2] */
  J[65] += dqdci; /* dwdot[O2]/d[O2] */
  J[66] -= dqdci; /* dwdot[HO2]/d[O2] */
  /* d()/d[HO2] */
  dqdci = +k_f * sc[3];
  J[75] -= dqdci; /* dwdot[OH]/d[HO2] */
  J[76] += dqdci; /* dwdot[H2O]/d[HO2] */
  J[77] += dqdci; /* dwdot[O2]/d[HO2] */
  J[78] -= dqdci; /* dwdot[HO2]/d[HO2] */
  /* d()/dT */
  J[135] -= dqdT; /* dwdot[OH]/dT */
  J[136] += dqdT; /* dwdot[H2O]/dT */
  J[137] += dqdT; /* dwdot[O2]/dT */
  J[138] -= dqdT; /* dwdot[HO2]/dT */

  /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[6] * sc[6];
  k_f = prefactor_units[21] * fwd_A[21] *
        exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
  dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
  /* reverse */
  phi_r = sc[5] * sc[7];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[6]) + (h_RT[5] + h_RT[7]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[5] += q;     /* O2 */
  wdot[6] -= 2 * q; /* HO2 */
  wdot[7] += q;     /* H2O2 */
  /* d()/d[O2] */
  dqdci = -k_r * sc[7];
  J[65] += dqdci;      /* dwdot[O2]/d[O2] */
  J[66] += -2 * dqdci; /* dwdot[HO2]/d[O2] */
  J[67] += dqdci;      /* dwdot[H2O2]/d[O2] */
  /* d()/d[HO2] */
  dqdci = +k_f * 2 * sc[6];
  J[77] += dqdci;      /* dwdot[O2]/d[HO2] */
  J[78] += -2 * dqdci; /* dwdot[HO2]/d[HO2] */
  J[79] += dqdci;      /* dwdot[H2O2]/d[HO2] */
  /* d()/d[H2O2] */
  dqdci = -k_r * sc[5];
  J[89] += dqdci;      /* dwdot[O2]/d[H2O2] */
  J[90] += -2 * dqdci; /* dwdot[HO2]/d[H2O2] */
  J[91] += dqdci;      /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[137] += dqdT;      /* dwdot[O2]/dT */
  J[138] += -2 * dqdT; /* dwdot[HO2]/dT */
  J[139] += dqdT;      /* dwdot[H2O2]/dT */

  /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[6] * sc[6];
  k_f = prefactor_units[22] * fwd_A[22] *
        exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
  dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
  /* reverse */
  phi_r = sc[5] * sc[7];
  Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(2 * h_RT[6]) + (h_RT[5] + h_RT[7]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[5] += q;     /* O2 */
  wdot[6] -= 2 * q; /* HO2 */
  wdot[7] += q;     /* H2O2 */
  /* d()/d[O2] */
  dqdci = -k_r * sc[7];
  J[65] += dqdci;      /* dwdot[O2]/d[O2] */
  J[66] += -2 * dqdci; /* dwdot[HO2]/d[O2] */
  J[67] += dqdci;      /* dwdot[H2O2]/d[O2] */
  /* d()/d[HO2] */
  dqdci = +k_f * 2 * sc[6];
  J[77] += dqdci;      /* dwdot[O2]/d[HO2] */
  J[78] += -2 * dqdci; /* dwdot[HO2]/d[HO2] */
  J[79] += dqdci;      /* dwdot[H2O2]/d[HO2] */
  /* d()/d[H2O2] */
  dqdci = -k_r * sc[5];
  J[89] += dqdci;      /* dwdot[O2]/d[H2O2] */
  J[90] += -2 * dqdci; /* dwdot[HO2]/d[H2O2] */
  J[91] += dqdci;      /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[137] += dqdT;      /* dwdot[O2]/dT */
  J[138] += -2 * dqdT; /* dwdot[HO2]/dT */
  J[139] += dqdT;      /* dwdot[H2O2]/dT */

  /*reaction 24: H2O2 + H <=> H2O + OH */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[0] * sc[7];
  k_f = prefactor_units[23] * fwd_A[23] *
        exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
  dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
  /* reverse */
  phi_r = sc[3] * sc[4];
  Kc = exp(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[3] + h_RT[4]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[3] += q; /* OH */
  wdot[4] += q; /* H2O */
  wdot[7] -= q; /* H2O2 */
  /* d()/d[H] */
  dqdci = +k_f * sc[7];
  J[0] -= dqdci; /* dwdot[H]/d[H] */
  J[3] += dqdci; /* dwdot[OH]/d[H] */
  J[4] += dqdci; /* dwdot[H2O]/d[H] */
  J[7] -= dqdci; /* dwdot[H2O2]/d[H] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[4];
  J[36] -= dqdci; /* dwdot[H]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  J[40] += dqdci; /* dwdot[H2O]/d[OH] */
  J[43] -= dqdci; /* dwdot[H2O2]/d[OH] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[3];
  J[48] -= dqdci; /* dwdot[H]/d[H2O] */
  J[51] += dqdci; /* dwdot[OH]/d[H2O] */
  J[52] += dqdci; /* dwdot[H2O]/d[H2O] */
  J[55] -= dqdci; /* dwdot[H2O2]/d[H2O] */
  /* d()/d[H2O2] */
  dqdci = +k_f * sc[0];
  J[84] -= dqdci; /* dwdot[H]/d[H2O2] */
  J[87] += dqdci; /* dwdot[OH]/d[H2O2] */
  J[88] += dqdci; /* dwdot[H2O]/d[H2O2] */
  J[91] -= dqdci; /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */
  J[136] += dqdT; /* dwdot[H2O]/dT */
  J[139] -= dqdT; /* dwdot[H2O2]/dT */

  /*reaction 25: H2O2 + H <=> HO2 + H2 */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[0] * sc[7];
  k_f = prefactor_units[24] * fwd_A[24] *
        exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
  dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
  /* reverse */
  phi_r = sc[1] * sc[6];
  Kc = exp(g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[1] += q; /* H2 */
  wdot[6] += q; /* HO2 */
  wdot[7] -= q; /* H2O2 */
  /* d()/d[H] */
  dqdci = +k_f * sc[7];
  J[0] -= dqdci; /* dwdot[H]/d[H] */
  J[1] += dqdci; /* dwdot[H2]/d[H] */
  J[6] += dqdci; /* dwdot[HO2]/d[H] */
  J[7] -= dqdci; /* dwdot[H2O2]/d[H] */
  /* d()/d[H2] */
  dqdci = -k_r * sc[6];
  J[12] -= dqdci; /* dwdot[H]/d[H2] */
  J[13] += dqdci; /* dwdot[H2]/d[H2] */
  J[18] += dqdci; /* dwdot[HO2]/d[H2] */
  J[19] -= dqdci; /* dwdot[H2O2]/d[H2] */
  /* d()/d[HO2] */
  dqdci = -k_r * sc[1];
  J[72] -= dqdci; /* dwdot[H]/d[HO2] */
  J[73] += dqdci; /* dwdot[H2]/d[HO2] */
  J[78] += dqdci; /* dwdot[HO2]/d[HO2] */
  J[79] -= dqdci; /* dwdot[H2O2]/d[HO2] */
  /* d()/d[H2O2] */
  dqdci = +k_f * sc[0];
  J[84] -= dqdci; /* dwdot[H]/d[H2O2] */
  J[85] += dqdci; /* dwdot[H2]/d[H2O2] */
  J[90] += dqdci; /* dwdot[HO2]/d[H2O2] */
  J[91] -= dqdci; /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[133] += dqdT; /* dwdot[H2]/dT */
  J[138] += dqdT; /* dwdot[HO2]/dT */
  J[139] -= dqdT; /* dwdot[H2O2]/dT */

  /*reaction 26: H2O2 + O <=> OH + HO2 */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[2] * sc[7];
  k_f = prefactor_units[25] * fwd_A[25] *
        exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
  dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
  /* reverse */
  phi_r = sc[3] * sc[6];
  Kc = exp(g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[3] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[2] -= q; /* O */
  wdot[3] += q; /* OH */
  wdot[6] += q; /* HO2 */
  wdot[7] -= q; /* H2O2 */
  /* d()/d[O] */
  dqdci = +k_f * sc[7];
  J[26] -= dqdci; /* dwdot[O]/d[O] */
  J[27] += dqdci; /* dwdot[OH]/d[O] */
  J[30] += dqdci; /* dwdot[HO2]/d[O] */
  J[31] -= dqdci; /* dwdot[H2O2]/d[O] */
  /* d()/d[OH] */
  dqdci = -k_r * sc[6];
  J[38] -= dqdci; /* dwdot[O]/d[OH] */
  J[39] += dqdci; /* dwdot[OH]/d[OH] */
  J[42] += dqdci; /* dwdot[HO2]/d[OH] */
  J[43] -= dqdci; /* dwdot[H2O2]/d[OH] */
  /* d()/d[HO2] */
  dqdci = -k_r * sc[3];
  J[74] -= dqdci; /* dwdot[O]/d[HO2] */
  J[75] += dqdci; /* dwdot[OH]/d[HO2] */
  J[78] += dqdci; /* dwdot[HO2]/d[HO2] */
  J[79] -= dqdci; /* dwdot[H2O2]/d[HO2] */
  /* d()/d[H2O2] */
  dqdci = +k_f * sc[2];
  J[86] -= dqdci; /* dwdot[O]/d[H2O2] */
  J[87] += dqdci; /* dwdot[OH]/d[H2O2] */
  J[90] += dqdci; /* dwdot[HO2]/d[H2O2] */
  J[91] -= dqdci; /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[134] -= dqdT; /* dwdot[O]/dT */
  J[135] += dqdT; /* dwdot[OH]/dT */
  J[138] += dqdT; /* dwdot[HO2]/dT */
  J[139] -= dqdT; /* dwdot[H2O2]/dT */

  /*reaction 27: H2O2 + OH <=> HO2 + H2O */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[3] * sc[7];
  k_f = prefactor_units[26] * fwd_A[26] *
        exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
  dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
  /* reverse */
  phi_r = sc[4] * sc[6];
  Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[3] -= q; /* OH */
  wdot[4] += q; /* H2O */
  wdot[6] += q; /* HO2 */
  wdot[7] -= q; /* H2O2 */
  /* d()/d[OH] */
  dqdci = +k_f * sc[7];
  J[39] -= dqdci; /* dwdot[OH]/d[OH] */
  J[40] += dqdci; /* dwdot[H2O]/d[OH] */
  J[42] += dqdci; /* dwdot[HO2]/d[OH] */
  J[43] -= dqdci; /* dwdot[H2O2]/d[OH] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[6];
  J[51] -= dqdci; /* dwdot[OH]/d[H2O] */
  J[52] += dqdci; /* dwdot[H2O]/d[H2O] */
  J[54] += dqdci; /* dwdot[HO2]/d[H2O] */
  J[55] -= dqdci; /* dwdot[H2O2]/d[H2O] */
  /* d()/d[HO2] */
  dqdci = -k_r * sc[4];
  J[75] -= dqdci; /* dwdot[OH]/d[HO2] */
  J[76] += dqdci; /* dwdot[H2O]/d[HO2] */
  J[78] += dqdci; /* dwdot[HO2]/d[HO2] */
  J[79] -= dqdci; /* dwdot[H2O2]/d[HO2] */
  /* d()/d[H2O2] */
  dqdci = +k_f * sc[3];
  J[87] -= dqdci; /* dwdot[OH]/d[H2O2] */
  J[88] += dqdci; /* dwdot[H2O]/d[H2O2] */
  J[90] += dqdci; /* dwdot[HO2]/d[H2O2] */
  J[91] -= dqdci; /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[135] -= dqdT; /* dwdot[OH]/dT */
  J[136] += dqdT; /* dwdot[H2O]/dT */
  J[138] += dqdT; /* dwdot[HO2]/dT */
  J[139] -= dqdT; /* dwdot[H2O2]/dT */

  /*reaction 28: H2O2 + OH <=> HO2 + H2O */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[3] * sc[7];
  k_f = prefactor_units[27] * fwd_A[27] *
        exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
  dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
  /* reverse */
  phi_r = sc[4] * sc[6];
  Kc = exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[6]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[3] -= q; /* OH */
  wdot[4] += q; /* H2O */
  wdot[6] += q; /* HO2 */
  wdot[7] -= q; /* H2O2 */
  /* d()/d[OH] */
  dqdci = +k_f * sc[7];
  J[39] -= dqdci; /* dwdot[OH]/d[OH] */
  J[40] += dqdci; /* dwdot[H2O]/d[OH] */
  J[42] += dqdci; /* dwdot[HO2]/d[OH] */
  J[43] -= dqdci; /* dwdot[H2O2]/d[OH] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[6];
  J[51] -= dqdci; /* dwdot[OH]/d[H2O] */
  J[52] += dqdci; /* dwdot[H2O]/d[H2O] */
  J[54] += dqdci; /* dwdot[HO2]/d[H2O] */
  J[55] -= dqdci; /* dwdot[H2O2]/d[H2O] */
  /* d()/d[HO2] */
  dqdci = -k_r * sc[4];
  J[75] -= dqdci; /* dwdot[OH]/d[HO2] */
  J[76] += dqdci; /* dwdot[H2O]/d[HO2] */
  J[78] += dqdci; /* dwdot[HO2]/d[HO2] */
  J[79] -= dqdci; /* dwdot[H2O2]/d[HO2] */
  /* d()/d[H2O2] */
  dqdci = +k_f * sc[3];
  J[87] -= dqdci; /* dwdot[OH]/d[H2O2] */
  J[88] += dqdci; /* dwdot[H2O]/d[H2O2] */
  J[90] += dqdci; /* dwdot[HO2]/d[H2O2] */
  J[91] -= dqdci; /* dwdot[H2O2]/d[H2O2] */
  /* d()/dT */
  J[135] -= dqdT; /* dwdot[OH]/dT */
  J[136] += dqdT; /* dwdot[H2O]/dT */
  J[138] += dqdT; /* dwdot[HO2]/dT */
  J[139] -= dqdT; /* dwdot[H2O2]/dT */

  /*reaction 29: HO2 + H <=> O + H2O */
  /*a non-third-body and non-pressure-fall-off reaction */
  /* forward */
  phi_f = sc[0] * sc[6];
  k_f = prefactor_units[28] * fwd_A[28] *
        exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
  dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
  /* reverse */
  phi_r = sc[2] * sc[4];
  Kc = exp(g_RT[0] - g_RT[2] - g_RT[4] + g_RT[6]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[2] + h_RT[4]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  /* rate of progress */
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  /* update wdot */
  wdot[0] -= q; /* H */
  wdot[2] += q; /* O */
  wdot[4] += q; /* H2O */
  wdot[6] -= q; /* HO2 */
  /* d()/d[H] */
  dqdci = +k_f * sc[6];
  J[0] -= dqdci; /* dwdot[H]/d[H] */
  J[2] += dqdci; /* dwdot[O]/d[H] */
  J[4] += dqdci; /* dwdot[H2O]/d[H] */
  J[6] -= dqdci; /* dwdot[HO2]/d[H] */
  /* d()/d[O] */
  dqdci = -k_r * sc[4];
  J[24] -= dqdci; /* dwdot[H]/d[O] */
  J[26] += dqdci; /* dwdot[O]/d[O] */
  J[28] += dqdci; /* dwdot[H2O]/d[O] */
  J[30] -= dqdci; /* dwdot[HO2]/d[O] */
  /* d()/d[H2O] */
  dqdci = -k_r * sc[2];
  J[48] -= dqdci; /* dwdot[H]/d[H2O] */
  J[50] += dqdci; /* dwdot[O]/d[H2O] */
  J[52] += dqdci; /* dwdot[H2O]/d[H2O] */
  J[54] -= dqdci; /* dwdot[HO2]/d[H2O] */
  /* d()/d[HO2] */
  dqdci = +k_f * sc[0];
  J[72] -= dqdci; /* dwdot[H]/d[HO2] */
  J[74] += dqdci; /* dwdot[O]/d[HO2] */
  J[76] += dqdci; /* dwdot[H2O]/d[HO2] */
  J[78] -= dqdci; /* dwdot[HO2]/d[HO2] */
  /* d()/dT */
  J[132] -= dqdT; /* dwdot[H]/dT */
  J[134] += dqdT; /* dwdot[O]/dT */
  J[136] += dqdT; /* dwdot[H2O]/dT */
  J[138] -= dqdT; /* dwdot[HO2]/dT */

  double c_R[11], dcRdT[11], e_RT[11];
  double* eh_RT;
  if (HP) {
    cp_R(c_R, tc);
    dcvpRdT(dcRdT, tc);
    eh_RT = &h_RT[0];
  } else {
    cv_R(c_R, tc);
    dcvpRdT(dcRdT, tc);
    speciesInternalEnergy(e_RT, tc);
    eh_RT = &e_RT[0];
  }

  double cmix = 0.0, ehmix = 0.0, dcmixdT = 0.0, dehmixdT = 0.0;
  for (int k = 0; k < 11; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[132 + k];
  }

  double cmixinv = 1.0 / cmix;
  double tmp1 = ehmix * cmixinv;
  double tmp3 = cmixinv * T;
  double tmp2 = tmp1 * tmp3;
  double dehmixdc;
  /* dTdot/d[X] */
  for (int k = 0; k < 11; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 11; ++m) {
      dehmixdc += eh_RT[m] * J[k * 12 + m];
    }
    J[k * 12 + 11] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  /* dTdot/dT */
  J[143] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
}

/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void
dcvpRdT(double* species, double* tc)
{

  /*temperature */
  double T = tc[1];

  /*species with midpoint at T=1000 kelvin */
  if (T < 1000) {
    /*species 0: H */
    species[0] = +0.00000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3];
    /*species 1: H2 */
    species[1] = +8.24944200e-04 - 1.62860300e-06 * tc[1] -
                 2.84263020e-10 * tc[2] + 1.65394880e-12 * tc[3];
    /*species 2: O */
    species[2] = -1.63816600e-03 + 4.84206400e-06 * tc[1] -
                 4.80852900e-09 * tc[2] + 1.55627840e-12 * tc[3];
    /*species 3: OH */
    species[3] = -3.22544939e-03 + 1.30552938e-05 * tc[1] -
                 1.73956093e-08 * tc[2] + 8.24949516e-12 * tc[3];
    /*species 4: H2O */
    species[4] = +3.47498200e-03 - 1.27093920e-05 * tc[1] +
                 2.09057430e-08 * tc[2] - 1.00263520e-11 * tc[3];
    /*species 5: O2 */
    species[5] = +1.12748600e-03 - 1.15123000e-06 * tc[1] +
                 3.94163100e-09 * tc[2] - 3.50742160e-12 * tc[3];
    /*species 6: HO2 */
    species[6] = -4.74912051e-03 + 4.23165782e-05 * tc[1] -
                 7.28291682e-08 * tc[2] + 3.71690050e-11 * tc[3];
    /*species 7: H2O2 */
    species[7] = +6.56922600e-03 - 2.97002600e-07 * tc[1] -
                 1.38774180e-08 * tc[2] + 9.88606000e-12 * tc[3];
    /*species 8: N2 */
    species[8] = +1.40824000e-03 - 7.92644400e-06 * tc[1] +
                 1.69245450e-08 * tc[2] - 9.77942000e-12 * tc[3];
    /*species 9: AR */
    species[9] = +0.00000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3];
    /*species 10: HE */
    species[10] = +0.00000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3];
  } else {
    /*species 0: H */
    species[0] = +0.00000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3];
    /*species 1: H2 */
    species[1] = +7.00064400e-04 - 1.12676580e-07 * tc[1] -
                 2.76947340e-11 * tc[2] + 6.33100800e-15 * tc[3];
    /*species 2: O */
    species[2] = -2.75506200e-05 - 6.20560600e-09 * tc[1] +
                 1.36532010e-11 * tc[2] - 1.74722080e-15 * tc[3];
    /*species 3: OH */
    species[3] = +1.05650448e-03 - 5.18165516e-07 * tc[1] +
                 9.15656022e-11 * tc[2] - 5.32783504e-15 * tc[3];
    /*species 4: H2O */
    species[4] = +3.05629300e-03 - 1.74605200e-06 * tc[1] +
                 3.60298800e-10 * tc[2] - 2.55664720e-14 * tc[3];
    /*species 5: O2 */
    species[5] = +6.13519700e-04 - 2.51768400e-07 * tc[1] +
                 5.32584300e-11 * tc[2] - 4.54574000e-15 * tc[3];
    /*species 6: HO2 */
    species[6] = +2.23982013e-03 - 1.26731630e-06 * tc[1] +
                 3.42739110e-10 * tc[2] - 4.31634140e-14 * tc[3];
    /*species 7: H2O2 */
    species[7] = +4.33613600e-03 - 2.94937800e-06 * tc[1] +
                 7.04671200e-10 * tc[2] - 5.72661600e-14 * tc[3];
    /*species 8: N2 */
    species[8] = +1.48797700e-03 - 1.13695220e-06 * tc[1] +
                 3.02911200e-10 * tc[2] - 2.70134040e-14 * tc[3];
    /*species 9: AR */
    species[9] = +0.00000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3];
    /*species 10: HE */
    species[10] = +0.00000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3];
  }
  return;
}

/*compute the progress rate for each reaction */
void
progressRate(double* qdot, double* sc, double T)
{
  double tc[] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; /*temperature cache */
  double invT = 1.0 / tc[1];

  if (T != T_save) {
    T_save = T;
    comp_k_f(tc, invT, k_f_save);
    comp_Kc(tc, invT, Kc_save);
  }

  double q_f[29], q_r[29];
  comp_qfqr(q_f, q_r, sc, tc, invT);

  for (int i = 0; i < 29; ++i) {
    qdot[i] = q_f[i] - q_r[i];
  }

  return;
}

/*compute the progress rate for each reaction */
void
progressRateFR(double* q_f, double* q_r, double* sc, double T)
{
  double tc[] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; /*temperature cache */
  double invT = 1.0 / tc[1];

  if (T != T_save) {
    T_save = T;
    comp_k_f(tc, invT, k_f_save);
    comp_Kc(tc, invT, Kc_save);
  }

  comp_qfqr(q_f, q_r, sc, tc, invT);

  return;
}

/*compute the equilibrium constants for each reaction */
void
equilibriumConstants(double* kc, double* g_RT, double T)
{
  /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
  double refC = 101325 / 8.31451 / T;

  /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
  kc[0] = 1.0 / (refC)*exp((g_RT[0] + g_RT[5]) - (g_RT[6]));

  /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
  kc[1] = refC * exp((g_RT[7]) - (g_RT[3] + g_RT[3]));

  /*reaction 3: H2 + M <=> H + H + M */
  kc[2] = refC * exp((g_RT[1]) - (g_RT[0] + g_RT[0]));

  /*reaction 4: O + O + M <=> O2 + M */
  kc[3] = 1.0 / (refC)*exp((g_RT[2] + g_RT[2]) - (g_RT[5]));

  /*reaction 5: O + H + M <=> OH + M */
  kc[4] = 1.0 / (refC)*exp((g_RT[2] + g_RT[0]) - (g_RT[3]));

  /*reaction 6: H2O + M <=> H + OH + M */
  kc[5] = refC * exp((g_RT[4]) - (g_RT[0] + g_RT[3]));

  /*reaction 7: O + OH + M <=> HO2 + M */
  kc[6] = 1.0 / (refC)*exp((g_RT[2] + g_RT[3]) - (g_RT[6]));

  /*reaction 8: H + O2 <=> O + OH */
  kc[7] = exp((g_RT[0] + g_RT[5]) - (g_RT[2] + g_RT[3]));

  /*reaction 9: O + H2 <=> H + OH */
  kc[8] = exp((g_RT[2] + g_RT[1]) - (g_RT[0] + g_RT[3]));

  /*reaction 10: O + H2 <=> H + OH */
  kc[9] = exp((g_RT[2] + g_RT[1]) - (g_RT[0] + g_RT[3]));

  /*reaction 11: H2 + OH <=> H2O + H */
  kc[10] = exp((g_RT[1] + g_RT[3]) - (g_RT[4] + g_RT[0]));

  /*reaction 12: OH + OH <=> O + H2O */
  kc[11] = exp((g_RT[3] + g_RT[3]) - (g_RT[2] + g_RT[4]));

  /*reaction 13: H2 + AR <=> H + H + AR */
  kc[12] = refC * exp((g_RT[1] + g_RT[9]) - (g_RT[0] + g_RT[0] + g_RT[9]));

  /*reaction 14: H2 + HE <=> H + H + HE */
  kc[13] = refC * exp((g_RT[1] + g_RT[10]) - (g_RT[0] + g_RT[0] + g_RT[10]));

  /*reaction 15: O + O + AR <=> O2 + AR */
  kc[14] =
    1.0 / (refC)*exp((g_RT[2] + g_RT[2] + g_RT[9]) - (g_RT[5] + g_RT[9]));

  /*reaction 16: O + O + HE <=> O2 + HE */
  kc[15] =
    1.0 / (refC)*exp((g_RT[2] + g_RT[2] + g_RT[10]) - (g_RT[5] + g_RT[10]));

  /*reaction 17: H2O + H2O <=> H + OH + H2O */
  kc[16] = refC * exp((g_RT[4] + g_RT[4]) - (g_RT[0] + g_RT[3] + g_RT[4]));

  /*reaction 18: HO2 + H <=> H2 + O2 */
  kc[17] = exp((g_RT[6] + g_RT[0]) - (g_RT[1] + g_RT[5]));

  /*reaction 19: HO2 + H <=> OH + OH */
  kc[18] = exp((g_RT[6] + g_RT[0]) - (g_RT[3] + g_RT[3]));

  /*reaction 20: HO2 + O <=> O2 + OH */
  kc[19] = exp((g_RT[6] + g_RT[2]) - (g_RT[5] + g_RT[3]));

  /*reaction 21: HO2 + OH <=> H2O + O2 */
  kc[20] = exp((g_RT[6] + g_RT[3]) - (g_RT[4] + g_RT[5]));

  /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
  kc[21] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[5]));

  /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
  kc[22] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[5]));

  /*reaction 24: H2O2 + H <=> H2O + OH */
  kc[23] = exp((g_RT[7] + g_RT[0]) - (g_RT[4] + g_RT[3]));

  /*reaction 25: H2O2 + H <=> HO2 + H2 */
  kc[24] = exp((g_RT[7] + g_RT[0]) - (g_RT[6] + g_RT[1]));

  /*reaction 26: H2O2 + O <=> OH + HO2 */
  kc[25] = exp((g_RT[7] + g_RT[2]) - (g_RT[3] + g_RT[6]));

  /*reaction 27: H2O2 + OH <=> HO2 + H2O */
  kc[26] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[4]));

  /*reaction 28: H2O2 + OH <=> HO2 + H2O */
  kc[27] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[4]));

  /*reaction 29: HO2 + H <=> O + H2O */
  kc[28] = exp((g_RT[6] + g_RT[0]) - (g_RT[2] + g_RT[4]));

  return;
}

/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void
gibbs(double* species, double* tc)
{

  /*temperature */
  double T = tc[1];
  double invT = 1 / T;

  /*species with midpoint at T=1000 kelvin */
  if (T < 1000) {
    /*species 0: H */
    species[0] = +2.547163000000000e+04 * invT + 2.960117600000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    /*species 1: H2 */
    species[1] = -1.012521000000000e+03 * invT + 6.592218000000000e+00 -
                 3.298124000000000e+00 * tc[0] - 4.124721000000000e-04 * tc[1] +
                 1.357169166666667e-07 * tc[2] + 7.896194999999999e-12 * tc[3] -
                 2.067436000000000e-14 * tc[4];
    /*species 2: O */
    species[2] = +2.914764000000000e+04 * invT - 1.756599999999997e-02 -
                 2.946429000000000e+00 * tc[0] + 8.190830000000000e-04 * tc[1] -
                 4.035053333333333e-07 * tc[2] + 1.335702500000000e-10 * tc[3] -
                 1.945348000000000e-14 * tc[4];
    /*species 3: OH */
    species[3] = +3.346309130000000e+03 * invT + 4.815738570000000e+00 -
                 4.125305610000000e+00 * tc[0] + 1.612724695000000e-03 * tc[1] -
                 1.087941151666667e-06 * tc[2] + 4.832113691666666e-10 * tc[3] -
                 1.031186895000000e-13 * tc[4];
    /*species 4: H2O */
    species[4] = -3.020811000000000e+04 * invT + 7.966090000000001e-01 -
                 3.386842000000000e+00 * tc[0] - 1.737491000000000e-03 * tc[1] +
                 1.059116000000000e-06 * tc[2] - 5.807150833333333e-10 * tc[3] +
                 1.253294000000000e-13 * tc[4];
    /*species 5: O2 */
    species[5] = -1.005249000000000e+03 * invT - 2.821802000000000e+00 -
                 3.212936000000000e+00 * tc[0] - 5.637430000000000e-04 * tc[1] +
                 9.593583333333333e-08 * tc[2] - 1.094897500000000e-10 * tc[3] +
                 4.384277000000000e-14 * tc[4];
    /*species 6: HO2 */
    species[6] = +2.948080400000000e+02 * invT + 5.851355599999999e-01 -
                 4.301798010000000e+00 * tc[0] + 2.374560255000000e-03 * tc[1] -
                 3.526381516666666e-06 * tc[2] + 2.023032450000000e-09 * tc[3] -
                 4.646125620000001e-13 * tc[4];
    /*species 7: H2O2 */
    species[7] = -1.766315000000000e+04 * invT - 3.396609000000000e+00 -
                 3.388754000000000e+00 * tc[0] - 3.284613000000000e-03 * tc[1] +
                 2.475021666666666e-08 * tc[2] + 3.854838333333333e-10 * tc[3] -
                 1.235757500000000e-13 * tc[4];
    /*species 8: N2 */
    species[8] = -1.020900000000000e+03 * invT - 6.516950000000001e-01 -
                 3.298677000000000e+00 * tc[0] - 7.041200000000000e-04 * tc[1] +
                 6.605369999999999e-07 * tc[2] - 4.701262500000001e-10 * tc[3] +
                 1.222427500000000e-13 * tc[4];
    /*species 9: AR */
    species[9] = -7.453750000000000e+02 * invT - 1.866001000000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    /*species 10: HE */
    species[10] = -7.453750000000000e+02 * invT + 1.584651200000000e+00 -
                  2.500000000000000e+00 * tc[0] -
                  0.000000000000000e+00 * tc[1] -
                  0.000000000000000e+00 * tc[2] -
                  0.000000000000000e+00 * tc[3] - 0.000000000000000e+00 * tc[4];
  } else {
    /*species 0: H */
    species[0] = +2.547163000000000e+04 * invT + 2.960117600000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    /*species 1: H2 */
    species[1] = -8.350340000000000e+02 * invT + 4.346533000000000e+00 -
                 2.991423000000000e+00 * tc[0] - 3.500322000000000e-04 * tc[1] +
                 9.389715000000000e-09 * tc[2] + 7.692981666666667e-13 * tc[3] -
                 7.913760000000000e-17 * tc[4];
    /*species 2: O */
    species[2] = +2.923080000000000e+04 * invT - 2.378248000000000e+00 -
                 2.542060000000000e+00 * tc[0] + 1.377531000000000e-05 * tc[1] +
                 5.171338333333333e-10 * tc[2] - 3.792555833333334e-13 * tc[3] +
                 2.184026000000000e-17 * tc[4];
    /*species 3: OH */
    species[3] = +3.683628750000000e+03 * invT - 2.836911870000000e+00 -
                 2.864728860000000e+00 * tc[0] - 5.282522400000000e-04 * tc[1] +
                 4.318045966666667e-08 * tc[2] - 2.543488950000000e-12 * tc[3] +
                 6.659793800000000e-17 * tc[4];
    /*species 4: H2O */
    species[4] = -2.989921000000000e+04 * invT - 4.190671000000000e+00 -
                 2.672146000000000e+00 * tc[0] - 1.528146500000000e-03 * tc[1] +
                 1.455043333333333e-07 * tc[2] - 1.000830000000000e-11 * tc[3] +
                 3.195809000000000e-16 * tc[4];
    /*species 5: O2 */
    species[5] = -1.233930000000000e+03 * invT + 5.084119999999999e-01 -
                 3.697578000000000e+00 * tc[0] - 3.067598500000000e-04 * tc[1] +
                 2.098070000000000e-08 * tc[2] - 1.479400833333333e-12 * tc[3] +
                 5.682175000000001e-17 * tc[4];
    /*species 6: HO2 */
    species[6] = +1.118567130000000e+02 * invT + 2.321087500000001e-01 -
                 4.017210900000000e+00 * tc[0] - 1.119910065000000e-03 * tc[1] +
                 1.056096916666667e-07 * tc[2] - 9.520530833333334e-12 * tc[3] +
                 5.395426750000000e-16 * tc[4];
    /*species 7: H2O2 */
    species[7] = -1.800696000000000e+04 * invT + 4.072030000000000e+00 -
                 4.573167000000000e+00 * tc[0] - 2.168068000000000e-03 * tc[1] +
                 2.457815000000000e-07 * tc[2] - 1.957420000000000e-11 * tc[3] +
                 7.158270000000000e-16 * tc[4];
    /*species 8: N2 */
    species[8] = -9.227977000000000e+02 * invT - 3.053888000000000e+00 -
                 2.926640000000000e+00 * tc[0] - 7.439885000000000e-04 * tc[1] +
                 9.474601666666666e-08 * tc[2] - 8.414199999999999e-12 * tc[3] +
                 3.376675500000000e-16 * tc[4];
    /*species 9: AR */
    species[9] = -7.453750000000000e+02 * invT - 1.866001000000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    /*species 10: HE */
    species[10] = -7.453750000000000e+02 * invT + 1.584651100000000e+00 -
                  2.500000000000000e+00 * tc[0] -
                  0.000000000000000e+00 * tc[1] -
                  0.000000000000000e+00 * tc[2] -
                  0.000000000000000e+00 * tc[3] - 0.000000000000000e+00 * tc[4];
  }
  return;
}

/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void
helmholtz(double* species, double* tc)
{

  /*temperature */
  double T = tc[1];
  double invT = 1 / T;

  /*species with midpoint at T=1000 kelvin */
  if (T < 1000) {
    /*species 0: H */
    species[0] = +2.54716300e+04 * invT + 1.96011760e+00 -
                 2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                 0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                 0.00000000e+00 * tc[4];
    /*species 1: H2 */
    species[1] = -1.01252100e+03 * invT + 5.59221800e+00 -
                 3.29812400e+00 * tc[0] - 4.12472100e-04 * tc[1] +
                 1.35716917e-07 * tc[2] + 7.89619500e-12 * tc[3] -
                 2.06743600e-14 * tc[4];
    /*species 2: O */
    species[2] = +2.91476400e+04 * invT - 1.01756600e+00 -
                 2.94642900e+00 * tc[0] + 8.19083000e-04 * tc[1] -
                 4.03505333e-07 * tc[2] + 1.33570250e-10 * tc[3] -
                 1.94534800e-14 * tc[4];
    /*species 3: OH */
    species[3] = +3.34630913e+03 * invT + 3.81573857e+00 -
                 4.12530561e+00 * tc[0] + 1.61272470e-03 * tc[1] -
                 1.08794115e-06 * tc[2] + 4.83211369e-10 * tc[3] -
                 1.03118689e-13 * tc[4];
    /*species 4: H2O */
    species[4] = -3.02081100e+04 * invT - 2.03391000e-01 -
                 3.38684200e+00 * tc[0] - 1.73749100e-03 * tc[1] +
                 1.05911600e-06 * tc[2] - 5.80715083e-10 * tc[3] +
                 1.25329400e-13 * tc[4];
    /*species 5: O2 */
    species[5] = -1.00524900e+03 * invT - 3.82180200e+00 -
                 3.21293600e+00 * tc[0] - 5.63743000e-04 * tc[1] +
                 9.59358333e-08 * tc[2] - 1.09489750e-10 * tc[3] +
                 4.38427700e-14 * tc[4];
    /*species 6: HO2 */
    species[6] = +2.94808040e+02 * invT - 4.14864440e-01 -
                 4.30179801e+00 * tc[0] + 2.37456025e-03 * tc[1] -
                 3.52638152e-06 * tc[2] + 2.02303245e-09 * tc[3] -
                 4.64612562e-13 * tc[4];
    /*species 7: H2O2 */
    species[7] = -1.76631500e+04 * invT - 4.39660900e+00 -
                 3.38875400e+00 * tc[0] - 3.28461300e-03 * tc[1] +
                 2.47502167e-08 * tc[2] + 3.85483833e-10 * tc[3] -
                 1.23575750e-13 * tc[4];
    /*species 8: N2 */
    species[8] = -1.02090000e+03 * invT - 1.65169500e+00 -
                 3.29867700e+00 * tc[0] - 7.04120000e-04 * tc[1] +
                 6.60537000e-07 * tc[2] - 4.70126250e-10 * tc[3] +
                 1.22242750e-13 * tc[4];
    /*species 9: AR */
    species[9] = -7.45375000e+02 * invT - 2.86600100e+00 -
                 2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                 0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                 0.00000000e+00 * tc[4];
    /*species 10: HE */
    species[10] = -7.45375000e+02 * invT + 5.84651200e-01 -
                  2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                  0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                  0.00000000e+00 * tc[4];
  } else {
    /*species 0: H */
    species[0] = +2.54716300e+04 * invT + 1.96011760e+00 -
                 2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                 0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                 0.00000000e+00 * tc[4];
    /*species 1: H2 */
    species[1] = -8.35034000e+02 * invT + 3.34653300e+00 -
                 2.99142300e+00 * tc[0] - 3.50032200e-04 * tc[1] +
                 9.38971500e-09 * tc[2] + 7.69298167e-13 * tc[3] -
                 7.91376000e-17 * tc[4];
    /*species 2: O */
    species[2] = +2.92308000e+04 * invT - 3.37824800e+00 -
                 2.54206000e+00 * tc[0] + 1.37753100e-05 * tc[1] +
                 5.17133833e-10 * tc[2] - 3.79255583e-13 * tc[3] +
                 2.18402600e-17 * tc[4];
    /*species 3: OH */
    species[3] = +3.68362875e+03 * invT - 3.83691187e+00 -
                 2.86472886e+00 * tc[0] - 5.28252240e-04 * tc[1] +
                 4.31804597e-08 * tc[2] - 2.54348895e-12 * tc[3] +
                 6.65979380e-17 * tc[4];
    /*species 4: H2O */
    species[4] = -2.98992100e+04 * invT - 5.19067100e+00 -
                 2.67214600e+00 * tc[0] - 1.52814650e-03 * tc[1] +
                 1.45504333e-07 * tc[2] - 1.00083000e-11 * tc[3] +
                 3.19580900e-16 * tc[4];
    /*species 5: O2 */
    species[5] = -1.23393000e+03 * invT - 4.91588000e-01 -
                 3.69757800e+00 * tc[0] - 3.06759850e-04 * tc[1] +
                 2.09807000e-08 * tc[2] - 1.47940083e-12 * tc[3] +
                 5.68217500e-17 * tc[4];
    /*species 6: HO2 */
    species[6] = +1.11856713e+02 * invT - 7.67891250e-01 -
                 4.01721090e+00 * tc[0] - 1.11991006e-03 * tc[1] +
                 1.05609692e-07 * tc[2] - 9.52053083e-12 * tc[3] +
                 5.39542675e-16 * tc[4];
    /*species 7: H2O2 */
    species[7] = -1.80069600e+04 * invT + 3.07203000e+00 -
                 4.57316700e+00 * tc[0] - 2.16806800e-03 * tc[1] +
                 2.45781500e-07 * tc[2] - 1.95742000e-11 * tc[3] +
                 7.15827000e-16 * tc[4];
    /*species 8: N2 */
    species[8] = -9.22797700e+02 * invT - 4.05388800e+00 -
                 2.92664000e+00 * tc[0] - 7.43988500e-04 * tc[1] +
                 9.47460167e-08 * tc[2] - 8.41420000e-12 * tc[3] +
                 3.37667550e-16 * tc[4];
    /*species 9: AR */
    species[9] = -7.45375000e+02 * invT - 2.86600100e+00 -
                 2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                 0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                 0.00000000e+00 * tc[4];
    /*species 10: HE */
    species[10] = -7.45375000e+02 * invT + 5.84651100e-01 -
                  2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                  0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                  0.00000000e+00 * tc[4];
  }
  return;
}

/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void
cv_R(double* species, double* tc)
{

  /*temperature */
  double T = tc[1];

  /*species with midpoint at T=1000 kelvin */
  if (T < 1000) {
    /*species 0: H */
    species[0] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4];
    /*species 1: H2 */
    species[1] = +2.29812400e+00 + 8.24944200e-04 * tc[1] -
                 8.14301500e-07 * tc[2] - 9.47543400e-11 * tc[3] +
                 4.13487200e-13 * tc[4];
    /*species 2: O */
    species[2] = +1.94642900e+00 - 1.63816600e-03 * tc[1] +
                 2.42103200e-06 * tc[2] - 1.60284300e-09 * tc[3] +
                 3.89069600e-13 * tc[4];
    /*species 3: OH */
    species[3] = +3.12530561e+00 - 3.22544939e-03 * tc[1] +
                 6.52764691e-06 * tc[2] - 5.79853643e-09 * tc[3] +
                 2.06237379e-12 * tc[4];
    /*species 4: H2O */
    species[4] = +2.38684200e+00 + 3.47498200e-03 * tc[1] -
                 6.35469600e-06 * tc[2] + 6.96858100e-09 * tc[3] -
                 2.50658800e-12 * tc[4];
    /*species 5: O2 */
    species[5] = +2.21293600e+00 + 1.12748600e-03 * tc[1] -
                 5.75615000e-07 * tc[2] + 1.31387700e-09 * tc[3] -
                 8.76855400e-13 * tc[4];
    /*species 6: HO2 */
    species[6] = +3.30179801e+00 - 4.74912051e-03 * tc[1] +
                 2.11582891e-05 * tc[2] - 2.42763894e-08 * tc[3] +
                 9.29225124e-12 * tc[4];
    /*species 7: H2O2 */
    species[7] = +2.38875400e+00 + 6.56922600e-03 * tc[1] -
                 1.48501300e-07 * tc[2] - 4.62580600e-09 * tc[3] +
                 2.47151500e-12 * tc[4];
    /*species 8: N2 */
    species[8] = +2.29867700e+00 + 1.40824000e-03 * tc[1] -
                 3.96322200e-06 * tc[2] + 5.64151500e-09 * tc[3] -
                 2.44485500e-12 * tc[4];
    /*species 9: AR */
    species[9] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4];
    /*species 10: HE */
    species[10] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4];
  } else {
    /*species 0: H */
    species[0] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4];
    /*species 1: H2 */
    species[1] = +1.99142300e+00 + 7.00064400e-04 * tc[1] -
                 5.63382900e-08 * tc[2] - 9.23157800e-12 * tc[3] +
                 1.58275200e-15 * tc[4];
    /*species 2: O */
    species[2] = +1.54206000e+00 - 2.75506200e-05 * tc[1] -
                 3.10280300e-09 * tc[2] + 4.55106700e-12 * tc[3] -
                 4.36805200e-16 * tc[4];
    /*species 3: OH */
    species[3] = +1.86472886e+00 + 1.05650448e-03 * tc[1] -
                 2.59082758e-07 * tc[2] + 3.05218674e-11 * tc[3] -
                 1.33195876e-15 * tc[4];
    /*species 4: H2O */
    species[4] = +1.67214600e+00 + 3.05629300e-03 * tc[1] -
                 8.73026000e-07 * tc[2] + 1.20099600e-10 * tc[3] -
                 6.39161800e-15 * tc[4];
    /*species 5: O2 */
    species[5] = +2.69757800e+00 + 6.13519700e-04 * tc[1] -
                 1.25884200e-07 * tc[2] + 1.77528100e-11 * tc[3] -
                 1.13643500e-15 * tc[4];
    /*species 6: HO2 */
    species[6] = +3.01721090e+00 + 2.23982013e-03 * tc[1] -
                 6.33658150e-07 * tc[2] + 1.14246370e-10 * tc[3] -
                 1.07908535e-14 * tc[4];
    /*species 7: H2O2 */
    species[7] = +3.57316700e+00 + 4.33613600e-03 * tc[1] -
                 1.47468900e-06 * tc[2] + 2.34890400e-10 * tc[3] -
                 1.43165400e-14 * tc[4];
    /*species 8: N2 */
    species[8] = +1.92664000e+00 + 1.48797700e-03 * tc[1] -
                 5.68476100e-07 * tc[2] + 1.00970400e-10 * tc[3] -
                 6.75335100e-15 * tc[4];
    /*species 9: AR */
    species[9] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4];
    /*species 10: HE */
    species[10] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4];
  }
  return;
}

/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void
cp_R(double* species, double* tc)
{

  /*temperature */
  double T = tc[1];

  /*species with midpoint at T=1000 kelvin */
  if (T < 1000) {
    /*species 0: H */
    species[0] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4];
    /*species 1: H2 */
    species[1] = +3.29812400e+00 + 8.24944200e-04 * tc[1] -
                 8.14301500e-07 * tc[2] - 9.47543400e-11 * tc[3] +
                 4.13487200e-13 * tc[4];
    /*species 2: O */
    species[2] = +2.94642900e+00 - 1.63816600e-03 * tc[1] +
                 2.42103200e-06 * tc[2] - 1.60284300e-09 * tc[3] +
                 3.89069600e-13 * tc[4];
    /*species 3: OH */
    species[3] = +4.12530561e+00 - 3.22544939e-03 * tc[1] +
                 6.52764691e-06 * tc[2] - 5.79853643e-09 * tc[3] +
                 2.06237379e-12 * tc[4];
    /*species 4: H2O */
    species[4] = +3.38684200e+00 + 3.47498200e-03 * tc[1] -
                 6.35469600e-06 * tc[2] + 6.96858100e-09 * tc[3] -
                 2.50658800e-12 * tc[4];
    /*species 5: O2 */
    species[5] = +3.21293600e+00 + 1.12748600e-03 * tc[1] -
                 5.75615000e-07 * tc[2] + 1.31387700e-09 * tc[3] -
                 8.76855400e-13 * tc[4];
    /*species 6: HO2 */
    species[6] = +4.30179801e+00 - 4.74912051e-03 * tc[1] +
                 2.11582891e-05 * tc[2] - 2.42763894e-08 * tc[3] +
                 9.29225124e-12 * tc[4];
    /*species 7: H2O2 */
    species[7] = +3.38875400e+00 + 6.56922600e-03 * tc[1] -
                 1.48501300e-07 * tc[2] - 4.62580600e-09 * tc[3] +
                 2.47151500e-12 * tc[4];
    /*species 8: N2 */
    species[8] = +3.29867700e+00 + 1.40824000e-03 * tc[1] -
                 3.96322200e-06 * tc[2] + 5.64151500e-09 * tc[3] -
                 2.44485500e-12 * tc[4];
    /*species 9: AR */
    species[9] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4];
    /*species 10: HE */
    species[10] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4];
  } else {
    /*species 0: H */
    species[0] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4];
    /*species 1: H2 */
    species[1] = +2.99142300e+00 + 7.00064400e-04 * tc[1] -
                 5.63382900e-08 * tc[2] - 9.23157800e-12 * tc[3] +
                 1.58275200e-15 * tc[4];
    /*species 2: O */
    species[2] = +2.54206000e+00 - 2.75506200e-05 * tc[1] -
                 3.10280300e-09 * tc[2] + 4.55106700e-12 * tc[3] -
                 4.36805200e-16 * tc[4];
    /*species 3: OH */
    species[3] = +2.86472886e+00 + 1.05650448e-03 * tc[1] -
                 2.59082758e-07 * tc[2] + 3.05218674e-11 * tc[3] -
                 1.33195876e-15 * tc[4];
    /*species 4: H2O */
    species[4] = +2.67214600e+00 + 3.05629300e-03 * tc[1] -
                 8.73026000e-07 * tc[2] + 1.20099600e-10 * tc[3] -
                 6.39161800e-15 * tc[4];
    /*species 5: O2 */
    species[5] = +3.69757800e+00 + 6.13519700e-04 * tc[1] -
                 1.25884200e-07 * tc[2] + 1.77528100e-11 * tc[3] -
                 1.13643500e-15 * tc[4];
    /*species 6: HO2 */
    species[6] = +4.01721090e+00 + 2.23982013e-03 * tc[1] -
                 6.33658150e-07 * tc[2] + 1.14246370e-10 * tc[3] -
                 1.07908535e-14 * tc[4];
    /*species 7: H2O2 */
    species[7] = +4.57316700e+00 + 4.33613600e-03 * tc[1] -
                 1.47468900e-06 * tc[2] + 2.34890400e-10 * tc[3] -
                 1.43165400e-14 * tc[4];
    /*species 8: N2 */
    species[8] = +2.92664000e+00 + 1.48797700e-03 * tc[1] -
                 5.68476100e-07 * tc[2] + 1.00970400e-10 * tc[3] -
                 6.75335100e-15 * tc[4];
    /*species 9: AR */
    species[9] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4];
    /*species 10: HE */
    species[10] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4];
  }
  return;
}

/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void
speciesInternalEnergy(double* species, double* tc)
{

  /*temperature */
  double T = tc[1];
  double invT = 1 / T;

  /*species with midpoint at T=1000 kelvin */
  if (T < 1000) {
    /*species 0: H */
    species[0] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] + 2.54716300e+04 * invT;
    /*species 1: H2 */
    species[1] = +2.29812400e+00 + 4.12472100e-04 * tc[1] -
                 2.71433833e-07 * tc[2] - 2.36885850e-11 * tc[3] +
                 8.26974400e-14 * tc[4] - 1.01252100e+03 * invT;
    /*species 2: O */
    species[2] = +1.94642900e+00 - 8.19083000e-04 * tc[1] +
                 8.07010667e-07 * tc[2] - 4.00710750e-10 * tc[3] +
                 7.78139200e-14 * tc[4] + 2.91476400e+04 * invT;
    /*species 3: OH */
    species[3] = +3.12530561e+00 - 1.61272470e-03 * tc[1] +
                 2.17588230e-06 * tc[2] - 1.44963411e-09 * tc[3] +
                 4.12474758e-13 * tc[4] + 3.34630913e+03 * invT;
    /*species 4: H2O */
    species[4] = +2.38684200e+00 + 1.73749100e-03 * tc[1] -
                 2.11823200e-06 * tc[2] + 1.74214525e-09 * tc[3] -
                 5.01317600e-13 * tc[4] - 3.02081100e+04 * invT;
    /*species 5: O2 */
    species[5] = +2.21293600e+00 + 5.63743000e-04 * tc[1] -
                 1.91871667e-07 * tc[2] + 3.28469250e-10 * tc[3] -
                 1.75371080e-13 * tc[4] - 1.00524900e+03 * invT;
    /*species 6: HO2 */
    species[6] = +3.30179801e+00 - 2.37456025e-03 * tc[1] +
                 7.05276303e-06 * tc[2] - 6.06909735e-09 * tc[3] +
                 1.85845025e-12 * tc[4] + 2.94808040e+02 * invT;
    /*species 7: H2O2 */
    species[7] = +2.38875400e+00 + 3.28461300e-03 * tc[1] -
                 4.95004333e-08 * tc[2] - 1.15645150e-09 * tc[3] +
                 4.94303000e-13 * tc[4] - 1.76631500e+04 * invT;
    /*species 8: N2 */
    species[8] = +2.29867700e+00 + 7.04120000e-04 * tc[1] -
                 1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
                 4.88971000e-13 * tc[4] - 1.02090000e+03 * invT;
    /*species 9: AR */
    species[9] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
    /*species 10: HE */
    species[10] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
  } else {
    /*species 0: H */
    species[0] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] + 2.54716300e+04 * invT;
    /*species 1: H2 */
    species[1] = +1.99142300e+00 + 3.50032200e-04 * tc[1] -
                 1.87794300e-08 * tc[2] - 2.30789450e-12 * tc[3] +
                 3.16550400e-16 * tc[4] - 8.35034000e+02 * invT;
    /*species 2: O */
    species[2] = +1.54206000e+00 - 1.37753100e-05 * tc[1] -
                 1.03426767e-09 * tc[2] + 1.13776675e-12 * tc[3] -
                 8.73610400e-17 * tc[4] + 2.92308000e+04 * invT;
    /*species 3: OH */
    species[3] = +1.86472886e+00 + 5.28252240e-04 * tc[1] -
                 8.63609193e-08 * tc[2] + 7.63046685e-12 * tc[3] -
                 2.66391752e-16 * tc[4] + 3.68362875e+03 * invT;
    /*species 4: H2O */
    species[4] = +1.67214600e+00 + 1.52814650e-03 * tc[1] -
                 2.91008667e-07 * tc[2] + 3.00249000e-11 * tc[3] -
                 1.27832360e-15 * tc[4] - 2.98992100e+04 * invT;
    /*species 5: O2 */
    species[5] = +2.69757800e+00 + 3.06759850e-04 * tc[1] -
                 4.19614000e-08 * tc[2] + 4.43820250e-12 * tc[3] -
                 2.27287000e-16 * tc[4] - 1.23393000e+03 * invT;
    /*species 6: HO2 */
    species[6] = +3.01721090e+00 + 1.11991006e-03 * tc[1] -
                 2.11219383e-07 * tc[2] + 2.85615925e-11 * tc[3] -
                 2.15817070e-15 * tc[4] + 1.11856713e+02 * invT;
    /*species 7: H2O2 */
    species[7] = +3.57316700e+00 + 2.16806800e-03 * tc[1] -
                 4.91563000e-07 * tc[2] + 5.87226000e-11 * tc[3] -
                 2.86330800e-15 * tc[4] - 1.80069600e+04 * invT;
    /*species 8: N2 */
    species[8] = +1.92664000e+00 + 7.43988500e-04 * tc[1] -
                 1.89492033e-07 * tc[2] + 2.52426000e-11 * tc[3] -
                 1.35067020e-15 * tc[4] - 9.22797700e+02 * invT;
    /*species 9: AR */
    species[9] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
    /*species 10: HE */
    species[10] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
  }
  return;
}

/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void
speciesEnthalpy(double* species, double* tc)
{

  /*temperature */
  double T = tc[1];
  double invT = 1 / T;

  /*species with midpoint at T=1000 kelvin */
  if (T < 1000) {
    /*species 0: H */
    species[0] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] + 2.54716300e+04 * invT;
    /*species 1: H2 */
    species[1] = +3.29812400e+00 + 4.12472100e-04 * tc[1] -
                 2.71433833e-07 * tc[2] - 2.36885850e-11 * tc[3] +
                 8.26974400e-14 * tc[4] - 1.01252100e+03 * invT;
    /*species 2: O */
    species[2] = +2.94642900e+00 - 8.19083000e-04 * tc[1] +
                 8.07010667e-07 * tc[2] - 4.00710750e-10 * tc[3] +
                 7.78139200e-14 * tc[4] + 2.91476400e+04 * invT;
    /*species 3: OH */
    species[3] = +4.12530561e+00 - 1.61272470e-03 * tc[1] +
                 2.17588230e-06 * tc[2] - 1.44963411e-09 * tc[3] +
                 4.12474758e-13 * tc[4] + 3.34630913e+03 * invT;
    /*species 4: H2O */
    species[4] = +3.38684200e+00 + 1.73749100e-03 * tc[1] -
                 2.11823200e-06 * tc[2] + 1.74214525e-09 * tc[3] -
                 5.01317600e-13 * tc[4] - 3.02081100e+04 * invT;
    /*species 5: O2 */
    species[5] = +3.21293600e+00 + 5.63743000e-04 * tc[1] -
                 1.91871667e-07 * tc[2] + 3.28469250e-10 * tc[3] -
                 1.75371080e-13 * tc[4] - 1.00524900e+03 * invT;
    /*species 6: HO2 */
    species[6] = +4.30179801e+00 - 2.37456025e-03 * tc[1] +
                 7.05276303e-06 * tc[2] - 6.06909735e-09 * tc[3] +
                 1.85845025e-12 * tc[4] + 2.94808040e+02 * invT;
    /*species 7: H2O2 */
    species[7] = +3.38875400e+00 + 3.28461300e-03 * tc[1] -
                 4.95004333e-08 * tc[2] - 1.15645150e-09 * tc[3] +
                 4.94303000e-13 * tc[4] - 1.76631500e+04 * invT;
    /*species 8: N2 */
    species[8] = +3.29867700e+00 + 7.04120000e-04 * tc[1] -
                 1.32107400e-06 * tc[2] + 1.41037875e-09 * tc[3] -
                 4.88971000e-13 * tc[4] - 1.02090000e+03 * invT;
    /*species 9: AR */
    species[9] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
    /*species 10: HE */
    species[10] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
  } else {
    /*species 0: H */
    species[0] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] + 2.54716300e+04 * invT;
    /*species 1: H2 */
    species[1] = +2.99142300e+00 + 3.50032200e-04 * tc[1] -
                 1.87794300e-08 * tc[2] - 2.30789450e-12 * tc[3] +
                 3.16550400e-16 * tc[4] - 8.35034000e+02 * invT;
    /*species 2: O */
    species[2] = +2.54206000e+00 - 1.37753100e-05 * tc[1] -
                 1.03426767e-09 * tc[2] + 1.13776675e-12 * tc[3] -
                 8.73610400e-17 * tc[4] + 2.92308000e+04 * invT;
    /*species 3: OH */
    species[3] = +2.86472886e+00 + 5.28252240e-04 * tc[1] -
                 8.63609193e-08 * tc[2] + 7.63046685e-12 * tc[3] -
                 2.66391752e-16 * tc[4] + 3.68362875e+03 * invT;
    /*species 4: H2O */
    species[4] = +2.67214600e+00 + 1.52814650e-03 * tc[1] -
                 2.91008667e-07 * tc[2] + 3.00249000e-11 * tc[3] -
                 1.27832360e-15 * tc[4] - 2.98992100e+04 * invT;
    /*species 5: O2 */
    species[5] = +3.69757800e+00 + 3.06759850e-04 * tc[1] -
                 4.19614000e-08 * tc[2] + 4.43820250e-12 * tc[3] -
                 2.27287000e-16 * tc[4] - 1.23393000e+03 * invT;
    /*species 6: HO2 */
    species[6] = +4.01721090e+00 + 1.11991006e-03 * tc[1] -
                 2.11219383e-07 * tc[2] + 2.85615925e-11 * tc[3] -
                 2.15817070e-15 * tc[4] + 1.11856713e+02 * invT;
    /*species 7: H2O2 */
    species[7] = +4.57316700e+00 + 2.16806800e-03 * tc[1] -
                 4.91563000e-07 * tc[2] + 5.87226000e-11 * tc[3] -
                 2.86330800e-15 * tc[4] - 1.80069600e+04 * invT;
    /*species 8: N2 */
    species[8] = +2.92664000e+00 + 7.43988500e-04 * tc[1] -
                 1.89492033e-07 * tc[2] + 2.52426000e-11 * tc[3] -
                 1.35067020e-15 * tc[4] - 9.22797700e+02 * invT;
    /*species 9: AR */
    species[9] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
    /*species 10: HE */
    species[10] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
  }
  return;
}

/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void
speciesEntropy(double* species, double* tc)
{

  /*temperature */
  double T = tc[1];

  /*species with midpoint at T=1000 kelvin */
  if (T < 1000) {
    /*species 0: H */
    species[0] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 4.60117600e-01;
    /*species 1: H2 */
    species[1] = +3.29812400e+00 * tc[0] + 8.24944200e-04 * tc[1] -
                 4.07150750e-07 * tc[2] - 3.15847800e-11 * tc[3] +
                 1.03371800e-13 * tc[4] - 3.29409400e+00;
    /*species 2: O */
    species[2] = +2.94642900e+00 * tc[0] - 1.63816600e-03 * tc[1] +
                 1.21051600e-06 * tc[2] - 5.34281000e-10 * tc[3] +
                 9.72674000e-14 * tc[4] + 2.96399500e+00;
    /*species 3: OH */
    species[3] = +4.12530561e+00 * tc[0] - 3.22544939e-03 * tc[1] +
                 3.26382346e-06 * tc[2] - 1.93284548e-09 * tc[3] +
                 5.15593447e-13 * tc[4] - 6.90432960e-01;
    /*species 4: H2O */
    species[4] = +3.38684200e+00 * tc[0] + 3.47498200e-03 * tc[1] -
                 3.17734800e-06 * tc[2] + 2.32286033e-09 * tc[3] -
                 6.26647000e-13 * tc[4] + 2.59023300e+00;
    /*species 5: O2 */
    species[5] = +3.21293600e+00 * tc[0] + 1.12748600e-03 * tc[1] -
                 2.87807500e-07 * tc[2] + 4.37959000e-10 * tc[3] -
                 2.19213850e-13 * tc[4] + 6.03473800e+00;
    /*species 6: HO2 */
    species[6] = +4.30179801e+00 * tc[0] - 4.74912051e-03 * tc[1] +
                 1.05791445e-05 * tc[2] - 8.09212980e-09 * tc[3] +
                 2.32306281e-12 * tc[4] + 3.71666245e+00;
    /*species 7: H2O2 */
    species[7] = +3.38875400e+00 * tc[0] + 6.56922600e-03 * tc[1] -
                 7.42506500e-08 * tc[2] - 1.54193533e-09 * tc[3] +
                 6.17878750e-13 * tc[4] + 6.78536300e+00;
    /*species 8: N2 */
    species[8] = +3.29867700e+00 * tc[0] + 1.40824000e-03 * tc[1] -
                 1.98161100e-06 * tc[2] + 1.88050500e-09 * tc[3] -
                 6.11213750e-13 * tc[4] + 3.95037200e+00;
    /*species 9: AR */
    species[9] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] + 4.36600100e+00;
    /*species 10: HE */
    species[10] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4] + 9.15348800e-01;
  } else {
    /*species 0: H */
    species[0] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 4.60117600e-01;
    /*species 1: H2 */
    species[1] = +2.99142300e+00 * tc[0] + 7.00064400e-04 * tc[1] -
                 2.81691450e-08 * tc[2] - 3.07719267e-12 * tc[3] +
                 3.95688000e-16 * tc[4] - 1.35511000e+00;
    /*species 2: O */
    species[2] = +2.54206000e+00 * tc[0] - 2.75506200e-05 * tc[1] -
                 1.55140150e-09 * tc[2] + 1.51702233e-12 * tc[3] -
                 1.09201300e-16 * tc[4] + 4.92030800e+00;
    /*species 3: OH */
    species[3] = +2.86472886e+00 * tc[0] + 1.05650448e-03 * tc[1] -
                 1.29541379e-07 * tc[2] + 1.01739558e-11 * tc[3] -
                 3.32989690e-16 * tc[4] + 5.70164073e+00;
    /*species 4: H2O */
    species[4] = +2.67214600e+00 * tc[0] + 3.05629300e-03 * tc[1] -
                 4.36513000e-07 * tc[2] + 4.00332000e-11 * tc[3] -
                 1.59790450e-15 * tc[4] + 6.86281700e+00;
    /*species 5: O2 */
    species[5] = +3.69757800e+00 * tc[0] + 6.13519700e-04 * tc[1] -
                 6.29421000e-08 * tc[2] + 5.91760333e-12 * tc[3] -
                 2.84108750e-16 * tc[4] + 3.18916600e+00;
    /*species 6: HO2 */
    species[6] = +4.01721090e+00 * tc[0] + 2.23982013e-03 * tc[1] -
                 3.16829075e-07 * tc[2] + 3.80821233e-11 * tc[3] -
                 2.69771337e-15 * tc[4] + 3.78510215e+00;
    /*species 7: H2O2 */
    species[7] = +4.57316700e+00 * tc[0] + 4.33613600e-03 * tc[1] -
                 7.37344500e-07 * tc[2] + 7.82968000e-11 * tc[3] -
                 3.57913500e-15 * tc[4] + 5.01137000e-01;
    /*species 8: N2 */
    species[8] = +2.92664000e+00 * tc[0] + 1.48797700e-03 * tc[1] -
                 2.84238050e-07 * tc[2] + 3.36568000e-11 * tc[3] -
                 1.68833775e-15 * tc[4] + 5.98052800e+00;
    /*species 9: AR */
    species[9] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] + 4.36600100e+00;
    /*species 10: HE */
    species[10] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                  0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                  0.00000000e+00 * tc[4] + 9.15348900e-01;
  }
  return;
}

/*save molecular weights into array */
void
molecularWeight(double* wt)
{
  wt[0] = 1.007970;  /*H */
  wt[1] = 2.015940;  /*H2 */
  wt[2] = 15.999400; /*O */
  wt[3] = 17.007370; /*OH */
  wt[4] = 18.015340; /*H2O */
  wt[5] = 31.998800; /*O2 */
  wt[6] = 33.006770; /*HO2 */
  wt[7] = 34.014740; /*H2O2 */
  wt[8] = 28.013400; /*N2 */
  wt[9] = 39.948000; /*AR */
  wt[10] = 4.002600; /*HE */

  return;
}

/*save atomic weights into array */
void
atomicWeight(double* awt)
{
  awt[0] = 1.007970;  /*H */
  awt[1] = 15.999400; /*O */
  awt[2] = 14.006700; /*N */
  awt[3] = 39.948000; /*AR */
  awt[4] = 4.002600;  /*HE */
  awt[5] = 12.011150; /*C */

  return;
}
/* get temperature given internal energy in mass units and mass fracs */
void
GET_T_GIVEN_EY(double* e, double* y, double* t, int* ierr)
{
#ifdef CONVERGENCE
  const int maxiter = 5000;
  const double tol = 1.e-12;
#else
  const int maxiter = 200;
  const double tol = 1.e-6;
#endif
  double ein = *e;
  double tmin = 90;   /*max lower bound for thermo def */
  double tmax = 4000; /*min upper bound for thermo def */
  double e1, emin, emax, cv, t1, dt;
  int i; /* loop counter */
  CKUBMS(&tmin, y, &emin);
  CKUBMS(&tmax, y, &emax);
  if (ein < emin) {
    /*Linear Extrapolation below tmin */
    CKCVBS(&tmin, y, &cv);
    *t = tmin - (emin - ein) / cv;
    *ierr = 1;
    return;
  }
  if (ein > emax) {
    /*Linear Extrapolation above tmax */
    CKCVBS(&tmax, y, &cv);
    *t = tmax - (emax - ein) / cv;
    *ierr = 1;
    return;
  }
  t1 = *t;
  if (t1 < tmin || t1 > tmax) {
    t1 = tmin + (tmax - tmin) / (emax - emin) * (ein - emin);
  }
  for (i = 0; i < maxiter; ++i) {
    CKUBMS(&t1, y, &e1);
    CKCVBS(&t1, y, &cv);
    dt = (ein - e1) / cv;
    if (dt > 100.) {
      dt = 100.;
    } else if (dt < -100.) {
      dt = -100.;
    } else if (fabs(dt) < tol)
      break;
    else if (t1 + dt == t1)
      break;
    t1 += dt;
  }
  *t = t1;
  *ierr = 0;
  return;
}

/* get temperature given enthalpy in mass units and mass fracs */
void
GET_T_GIVEN_HY(double* h, double* y, double* t, int* ierr)
{
#ifdef CONVERGENCE
  const int maxiter = 5000;
  const double tol = 1.e-12;
#else
  const int maxiter = 200;
  const double tol = 1.e-6;
#endif
  double hin = *h;
  double tmin = 90;   /*max lower bound for thermo def */
  double tmax = 4000; /*min upper bound for thermo def */
  double h1, hmin, hmax, cp, t1, dt;
  int i; /* loop counter */
  CKHBMS(&tmin, y, &hmin);
  CKHBMS(&tmax, y, &hmax);
  if (hin < hmin) {
    /*Linear Extrapolation below tmin */
    CKCPBS(&tmin, y, &cp);
    *t = tmin - (hmin - hin) / cp;
    *ierr = 1;
    return;
  }
  if (hin > hmax) {
    /*Linear Extrapolation above tmax */
    CKCPBS(&tmax, y, &cp);
    *t = tmax - (hmax - hin) / cp;
    *ierr = 1;
    return;
  }
  t1 = *t;
  if (t1 < tmin || t1 > tmax) {
    t1 = tmin + (tmax - tmin) / (hmax - hmin) * (hin - hmin);
  }
  for (i = 0; i < maxiter; ++i) {
    CKHBMS(&t1, y, &h1);
    CKCPBS(&t1, y, &cp);
    dt = (hin - h1) / cp;
    if (dt > 100.) {
      dt = 100.;
    } else if (dt < -100.) {
      dt = -100.;
    } else if (fabs(dt) < tol)
      break;
    else if (t1 + dt == t1)
      break;
    t1 += dt;
  }
  *t = t1;
  *ierr = 0;
  return;
}

/*compute the critical parameters for each species */
void
GET_CRITPARAMS(double* Tci, double* ai, double* bi, double* acentric_i)
{

  double EPS[11];
  double SIG[11];
  double wt[11];
  double avogadro = 6.02214199e23;
  double boltzmann = 1.3806503e-16; // we work in CGS
  double Rcst = 83.144598;          // in bar [CGS] !

  egtransetEPS(EPS);
  egtransetSIG(SIG);
  molecularWeight(wt);

  /*species 0: H */
  Tci[0] = 1.316 * EPS[0];
  ai[0] =
    (5.55 * pow(avogadro, 2.0) * EPS[0] * boltzmann * pow(1e-8 * SIG[0], 3.0)) /
    (pow(wt[0], 2.0));
  bi[0] = 0.855 * avogadro * pow(1e-8 * SIG[0], 3.0) / (wt[0]);
  acentric_i[0] = 0.0;

  /*species 1: H2 */
  /*Imported from NIST */
  Tci[1] = 33.145000;
  ai[1] = 1e6 * 0.42748 * pow(Rcst, 2.0) * pow(Tci[1], 2.0) /
          (pow(2.015880, 2.0) * 12.964000);
  bi[1] = 0.08664 * Rcst * Tci[1] / (2.015880 * 12.964000);
  acentric_i[1] = -0.219000;

  /*species 2: O */
  Tci[2] = 1.316 * EPS[2];
  ai[2] =
    (5.55 * pow(avogadro, 2.0) * EPS[2] * boltzmann * pow(1e-8 * SIG[2], 3.0)) /
    (pow(wt[2], 2.0));
  bi[2] = 0.855 * avogadro * pow(1e-8 * SIG[2], 3.0) / (wt[2]);
  acentric_i[2] = 0.0;

  /*species 3: OH */
  Tci[3] = 1.316 * EPS[3];
  ai[3] =
    (5.55 * pow(avogadro, 2.0) * EPS[3] * boltzmann * pow(1e-8 * SIG[3], 3.0)) /
    (pow(wt[3], 2.0));
  bi[3] = 0.855 * avogadro * pow(1e-8 * SIG[3], 3.0) / (wt[3]);
  acentric_i[3] = 0.0;

  /*species 4: H2O */
  /*Imported from NIST */
  Tci[4] = 647.096000;
  ai[4] = 1e6 * 0.42748 * pow(Rcst, 2.0) * pow(Tci[4], 2.0) /
          (pow(18.015340, 2.0) * 220.640000);
  bi[4] = 0.08664 * Rcst * Tci[4] / (18.015340 * 220.640000);
  acentric_i[4] = 0.344300;

  /*species 5: O2 */
  /*Imported from NIST */
  Tci[5] = 154.581000;
  ai[5] = 1e6 * 0.42748 * pow(Rcst, 2.0) * pow(Tci[5], 2.0) /
          (pow(31.998800, 2.0) * 50.430466);
  bi[5] = 0.08664 * Rcst * Tci[5] / (31.998800 * 50.430466);
  acentric_i[5] = 0.022200;

  /*species 6: HO2 */
  Tci[6] = 1.316 * EPS[6];
  ai[6] =
    (5.55 * pow(avogadro, 2.0) * EPS[6] * boltzmann * pow(1e-8 * SIG[6], 3.0)) /
    (pow(wt[6], 2.0));
  bi[6] = 0.855 * avogadro * pow(1e-8 * SIG[6], 3.0) / (wt[6]);
  acentric_i[6] = 0.0;

  /*species 7: H2O2 */
  Tci[7] = 1.316 * EPS[7];
  ai[7] =
    (5.55 * pow(avogadro, 2.0) * EPS[7] * boltzmann * pow(1e-8 * SIG[7], 3.0)) /
    (pow(wt[7], 2.0));
  bi[7] = 0.855 * avogadro * pow(1e-8 * SIG[7], 3.0) / (wt[7]);
  acentric_i[7] = 0.0;

  /*species 8: N2 */
  /*Imported from NIST */
  Tci[8] = 126.192000;
  ai[8] = 1e6 * 0.42748 * pow(Rcst, 2.0) * pow(Tci[8], 2.0) /
          (pow(28.013400, 2.0) * 33.958000);
  bi[8] = 0.08664 * Rcst * Tci[8] / (28.013400 * 33.958000);
  acentric_i[8] = 0.037200;

  /*species 9: AR */
  /*Imported from NIST */
  Tci[9] = 150.860000;
  ai[9] = 1e6 * 0.42748 * pow(Rcst, 2.0) * pow(Tci[9], 2.0) /
          (pow(39.948000, 2.0) * 48.980000);
  bi[9] = 0.08664 * Rcst * Tci[9] / (39.948000 * 48.980000);
  acentric_i[9] = -0.002000;

  /*species 10: HE */
  Tci[10] = 1.316 * EPS[10];
  ai[10] = (5.55 * pow(avogadro, 2.0) * EPS[10] * boltzmann *
            pow(1e-8 * SIG[10], 3.0)) /
           (pow(wt[10], 2.0));
  bi[10] = 0.855 * avogadro * pow(1e-8 * SIG[10], 3.0) / (wt[10]);
  acentric_i[10] = 0.0;

  return;
}

#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENIMC EGTRANSETLENIMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENIMC egtransetlenimc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENIMC egtransetlenimc_
#endif
void
egtransetLENIMC(int* LENIMC)
{
  *LENIMC = 47;
}

#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void
egtransetLENRMC(int* LENRMC)
{
  *LENRMC = 2728;
}

#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNO EGTRANSETNO
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNO egtransetno
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNO egtransetno_
#endif
void
egtransetNO(int* NO)
{
  *NO = 4;
}

#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKK EGTRANSETKK
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKK egtransetkk
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKK egtransetkk_
#endif
void
egtransetKK(int* KK)
{
  *KK = 11;
}

#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLITE EGTRANSETNLITE
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLITE egtransetnlite
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLITE egtransetnlite_
#endif
void
egtransetNLITE(int* NLITE)
{
  *NLITE = 3;
}

/*Patm in ergs/cm3 */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPATM EGTRANSETPATM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPATM egtransetpatm
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPATM egtransetpatm_
#endif
void
egtransetPATM(double* PATM)
{
  *PATM = 0.1013250000000000E+07;
}

/*the molecular weights in g/mol */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetWT EGTRANSETWT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetWT egtransetwt
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetWT egtransetwt_
#endif
void
egtransetWT(double* WT)
{
  WT[0] = 1.00797000E+00;
  WT[1] = 2.01594000E+00;
  WT[2] = 1.59994000E+01;
  WT[3] = 1.70073700E+01;
  WT[4] = 1.80153400E+01;
  WT[5] = 3.19988000E+01;
  WT[6] = 3.30067700E+01;
  WT[7] = 3.40147400E+01;
  WT[8] = 2.80134000E+01;
  WT[9] = 3.99480000E+01;
  WT[10] = 4.00260000E+00;
}

/*the lennard-jones potential well depth eps/kb in K */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
void
egtransetEPS(double* EPS)
{
  EPS[6] = 1.07400000E+02;
  EPS[7] = 1.07400000E+02;
  EPS[9] = 1.36500000E+02;
  EPS[0] = 1.45000000E+02;
  EPS[1] = 3.80000000E+01;
  EPS[8] = 9.75300000E+01;
  EPS[2] = 8.00000000E+01;
  EPS[3] = 8.00000000E+01;
  EPS[4] = 5.72400000E+02;
  EPS[5] = 1.07400000E+02;
  EPS[10] = 1.02000000E+01;
}

/*the lennard-jones collision diameter in Angstroms */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
void
egtransetSIG(double* SIG)
{
  SIG[6] = 3.45800000E+00;
  SIG[7] = 3.45800000E+00;
  SIG[9] = 3.33000000E+00;
  SIG[0] = 2.05000000E+00;
  SIG[1] = 2.92000000E+00;
  SIG[8] = 3.62100000E+00;
  SIG[2] = 2.75000000E+00;
  SIG[3] = 2.75000000E+00;
  SIG[4] = 2.60500000E+00;
  SIG[5] = 3.45800000E+00;
  SIG[10] = 2.57600000E+00;
}

/*the dipole moment in Debye */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetDIP EGTRANSETDIP
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetDIP egtransetdip
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetDIP egtransetdip_
#endif
void
egtransetDIP(double* DIP)
{
  DIP[6] = 0.00000000E+00;
  DIP[7] = 0.00000000E+00;
  DIP[9] = 0.00000000E+00;
  DIP[0] = 0.00000000E+00;
  DIP[1] = 0.00000000E+00;
  DIP[8] = 0.00000000E+00;
  DIP[2] = 0.00000000E+00;
  DIP[3] = 0.00000000E+00;
  DIP[4] = 1.84400000E+00;
  DIP[5] = 0.00000000E+00;
  DIP[10] = 0.00000000E+00;
}

/*the polarizability in cubic Angstroms */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPOL EGTRANSETPOL
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPOL egtransetpol
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPOL egtransetpol_
#endif
void
egtransetPOL(double* POL)
{
  POL[6] = 0.00000000E+00;
  POL[7] = 0.00000000E+00;
  POL[9] = 0.00000000E+00;
  POL[0] = 0.00000000E+00;
  POL[1] = 7.90000000E-01;
  POL[8] = 1.76000000E+00;
  POL[2] = 0.00000000E+00;
  POL[3] = 0.00000000E+00;
  POL[4] = 0.00000000E+00;
  POL[5] = 1.60000000E+00;
  POL[10] = 0.00000000E+00;
}

/*the rotational relaxation collision number at 298 K */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetZROT EGTRANSETZROT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetZROT egtransetzrot
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetZROT egtransetzrot_
#endif
void
egtransetZROT(double* ZROT)
{
  ZROT[6] = 1.00000000E+00;
  ZROT[7] = 3.80000000E+00;
  ZROT[9] = 0.00000000E+00;
  ZROT[0] = 0.00000000E+00;
  ZROT[1] = 2.80000000E+02;
  ZROT[8] = 4.00000000E+00;
  ZROT[2] = 0.00000000E+00;
  ZROT[3] = 0.00000000E+00;
  ZROT[4] = 4.00000000E+00;
  ZROT[5] = 3.80000000E+00;
  ZROT[10] = 0.00000000E+00;
}

/*0: monoatomic, 1: linear, 2: nonlinear */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLIN EGTRANSETNLIN
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLIN egtransetnlin
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLIN egtransetnlin_
#endif
void
egtransetNLIN(int* NLIN)
{
  NLIN[6] = 2;
  NLIN[7] = 2;
  NLIN[9] = 0;
  NLIN[0] = 0;
  NLIN[1] = 1;
  NLIN[8] = 1;
  NLIN[2] = 0;
  NLIN[3] = 1;
  NLIN[4] = 2;
  NLIN[5] = 1;
  NLIN[10] = 0;
}

/*Poly fits for the viscosities, dim NO*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFETA EGTRANSETCOFETA
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFETA egtransetcofeta
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFETA egtransetcofeta_
#endif
void
egtransetCOFETA(double* COFETA)
{
  COFETA[0] = -1.98744496E+01;
  COFETA[1] = 3.41660514E+00;
  COFETA[2] = -3.63206306E-01;
  COFETA[3] = 1.58671021E-02;
  COFETA[4] = -1.37549435E+01;
  COFETA[5] = 9.65530587E-01;
  COFETA[6] = -4.45720114E-02;
  COFETA[7] = 2.05871810E-03;
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
  COFETA[20] = -1.68118868E+01;
  COFETA[21] = 2.52362554E+00;
  COFETA[22] = -2.49309128E-01;
  COFETA[23] = 1.10211025E-02;
  COFETA[24] = -1.67963797E+01;
  COFETA[25] = 2.52362554E+00;
  COFETA[26] = -2.49309128E-01;
  COFETA[27] = 1.10211025E-02;
  COFETA[28] = -1.67813391E+01;
  COFETA[29] = 2.52362554E+00;
  COFETA[30] = -2.49309128E-01;
  COFETA[31] = 1.10211025E-02;
  COFETA[32] = -1.62526779E+01;
  COFETA[33] = 2.24839597E+00;
  COFETA[34] = -2.13428438E-01;
  COFETA[35] = 9.46192413E-03;
  COFETA[36] = -1.86067598E+01;
  COFETA[37] = 3.27402596E+00;
  COFETA[38] = -3.45827972E-01;
  COFETA[39] = 1.51622680E-02;
  COFETA[40] = -1.11555213E+01;
  COFETA[41] = 2.18772782E-01;
  COFETA[42] = 5.60263799E-02;
  COFETA[43] = -2.36018246E-03;
}

/*Poly fits for the conductivities, dim NO*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
void
egtransetCOFLAM(double* COFLAM)
{
  COFLAM[0] = -3.24539191E-01;
  COFLAM[1] = 3.41660514E+00;
  COFLAM[2] = -3.63206306E-01;
  COFLAM[3] = 1.58671021E-02;
  COFLAM[4] = 1.11035511E+01;
  COFLAM[5] = -1.31883912E+00;
  COFLAM[6] = 2.44042473E-01;
  COFLAM[7] = -8.99836359E-03;
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
  COFLAM[20] = -2.51296685E+00;
  COFLAM[21] = 3.15165687E+00;
  COFLAM[22] = -3.10007661E-01;
  COFLAM[23] = 1.34522321E-02;
  COFLAM[24] = 5.56023763E-01;
  COFLAM[25] = 1.59073590E+00;
  COFLAM[26] = -5.28053839E-02;
  COFLAM[27] = 4.07601571E-04;
  COFLAM[28] = 1.48801095E+00;
  COFLAM[29] = 1.06176238E+00;
  COFLAM[30] = 5.72195640E-02;
  COFLAM[31] = -6.38391491E-03;
  COFLAM[32] = 1.15507419E+01;
  COFLAM[33] = -2.91453917E+00;
  COFLAM[34] = 5.55045765E-01;
  COFLAM[35] = -2.75173485E-02;
  COFLAM[36] = -2.73648952E+00;
  COFLAM[37] = 3.27402596E+00;
  COFLAM[38] = -3.45827972E-01;
  COFLAM[39] = 1.51622680E-02;
  COFLAM[40] = 7.01538340E+00;
  COFLAM[41] = 2.18772782E-01;
  COFLAM[42] = 5.60263799E-02;
  COFLAM[43] = -2.36018246E-03;
}

/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFD EGTRANSETCOFD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFD egtransetcofd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFD egtransetcofd_
#endif
void
egtransetCOFD(double* COFD)
{
  COFD[0] = -1.43693056E+01;
  COFD[1] = 4.03992999E+00;
  COFD[2] = -3.08044800E-01;
  COFD[3] = 1.32757775E-02;
  COFD[4] = -1.11808682E+01;
  COFD[5] = 2.66936727E+00;
  COFD[6] = -1.34411514E-01;
  COFD[7] = 5.92957488E-03;
  COFD[8] = -1.31860117E+01;
  COFD[9] = 3.38003453E+00;
  COFD[10] = -2.25783856E-01;
  COFD[11] = 9.85028660E-03;
  COFD[12] = -1.31877711E+01;
  COFD[13] = 3.38003453E+00;
  COFD[14] = -2.25783856E-01;
  COFD[15] = 9.85028660E-03;
  COFD[16] = -1.93611051E+01;
  COFD[17] = 5.51579726E+00;
  COFD[18] = -4.76061961E-01;
  COFD[19] = 1.96329391E-02;
  COFD[20] = -1.43712864E+01;
  COFD[21] = 3.70920439E+00;
  COFD[22] = -2.67274113E-01;
  COFD[23] = 1.15967481E-02;
  COFD[24] = -1.43717529E+01;
  COFD[25] = 3.70920439E+00;
  COFD[26] = -2.67274113E-01;
  COFD[27] = 1.15967481E-02;
  COFD[28] = -1.43721922E+01;
  COFD[29] = 3.70920439E+00;
  COFD[30] = -2.67274113E-01;
  COFD[31] = 1.15967481E-02;
  COFD[32] = -1.40298830E+01;
  COFD[33] = 3.55837688E+00;
  COFD[34] = -2.47785790E-01;
  COFD[35] = 1.07555332E-02;
  COFD[36] = -1.51208119E+01;
  COFD[37] = 3.99904647E+00;
  COFD[38] = -3.03517220E-01;
  COFD[39] = 1.31117363E-02;
  COFD[40] = -9.71338331E+00;
  COFD[41] = 2.17561180E+00;
  COFD[42] = -7.28270090E-02;
  COFD[43] = 3.38302182E-03;
  COFD[44] = -1.11808682E+01;
  COFD[45] = 2.66936727E+00;
  COFD[46] = -1.34411514E-01;
  COFD[47] = 5.92957488E-03;
  COFD[48] = -1.02395222E+01;
  COFD[49] = 2.15403244E+00;
  COFD[50] = -6.97480266E-02;
  COFD[51] = 3.23666871E-03;
  COFD[52] = -1.06250182E+01;
  COFD[53] = 2.15849701E+00;
  COFD[54] = -6.53886401E-02;
  COFD[55] = 2.81453370E-03;
  COFD[56] = -1.06283453E+01;
  COFD[57] = 2.15849701E+00;
  COFD[58] = -6.53886401E-02;
  COFD[59] = 2.81453370E-03;
  COFD[60] = -1.68758926E+01;
  COFD[61] = 4.49460303E+00;
  COFD[62] = -3.64766132E-01;
  COFD[63] = 1.56457153E-02;
  COFD[64] = -1.15797750E+01;
  COFD[65] = 2.43235504E+00;
  COFD[66] = -1.02890179E-01;
  COFD[67] = 4.52903603E-03;
  COFD[68] = -1.15806808E+01;
  COFD[69] = 2.43235504E+00;
  COFD[70] = -1.02890179E-01;
  COFD[71] = 4.52903603E-03;
  COFD[72] = -1.15815344E+01;
  COFD[73] = 2.43235504E+00;
  COFD[74] = -1.02890179E-01;
  COFD[75] = 4.52903603E-03;
  COFD[76] = -1.13253458E+01;
  COFD[77] = 2.31195095E+00;
  COFD[78] = -8.63988037E-02;
  COFD[79] = 3.77573452E-03;
  COFD[80] = -1.20638601E+01;
  COFD[81] = 2.63303536E+00;
  COFD[82] = -1.29792632E-01;
  COFD[83] = 5.73363738E-03;
  COFD[84] = -9.86429034E+00;
  COFD[85] = 2.05348746E+00;
  COFD[86] = -5.90289007E-02;
  COFD[87] = 2.89596157E-03;
  COFD[88] = -1.31860117E+01;
  COFD[89] = 3.38003453E+00;
  COFD[90] = -2.25783856E-01;
  COFD[91] = 9.85028660E-03;
  COFD[92] = -1.06250182E+01;
  COFD[93] = 2.15849701E+00;
  COFD[94] = -6.53886401E-02;
  COFD[95] = 2.81453370E-03;
  COFD[96] = -1.29877365E+01;
  COFD[97] = 2.80841511E+00;
  COFD[98] = -1.52629888E-01;
  COFD[99] = 6.72604927E-03;
  COFD[100] = -1.30027772E+01;
  COFD[101] = 2.80841511E+00;
  COFD[102] = -1.52629888E-01;
  COFD[103] = 6.72604927E-03;
  COFD[104] = -1.91096797E+01;
  COFD[105] = 5.02608697E+00;
  COFD[106] = -4.26959993E-01;
  COFD[107] = 1.80709910E-02;
  COFD[108] = -1.40864894E+01;
  COFD[109] = 3.07458927E+00;
  COFD[110] = -1.86899591E-01;
  COFD[111] = 8.19829781E-03;
  COFD[112] = -1.40916052E+01;
  COFD[113] = 3.07458927E+00;
  COFD[114] = -1.86899591E-01;
  COFD[115] = 8.19829781E-03;
  COFD[116] = -1.40964661E+01;
  COFD[117] = 3.07458927E+00;
  COFD[118] = -1.86899591E-01;
  COFD[119] = 8.19829781E-03;
  COFD[120] = -1.38756407E+01;
  COFD[121] = 2.98558426E+00;
  COFD[122] = -1.75507216E-01;
  COFD[123] = 7.71173691E-03;
  COFD[124] = -1.47082523E+01;
  COFD[125] = 3.30683499E+00;
  COFD[126] = -2.16378602E-01;
  COFD[127] = 9.44670561E-03;
  COFD[128] = -9.70779324E+00;
  COFD[129] = 1.77912272E+00;
  COFD[130] = -1.67349571E-02;
  COFD[131] = 7.45446845E-04;
  COFD[132] = -1.31877711E+01;
  COFD[133] = 3.38003453E+00;
  COFD[134] = -2.25783856E-01;
  COFD[135] = 9.85028660E-03;
  COFD[136] = -1.06283453E+01;
  COFD[137] = 2.15849701E+00;
  COFD[138] = -6.53886401E-02;
  COFD[139] = 2.81453370E-03;
  COFD[140] = -1.30027772E+01;
  COFD[141] = 2.80841511E+00;
  COFD[142] = -1.52629888E-01;
  COFD[143] = 6.72604927E-03;
  COFD[144] = -1.30182843E+01;
  COFD[145] = 2.80841511E+00;
  COFD[146] = -1.52629888E-01;
  COFD[147] = 6.72604927E-03;
  COFD[148] = -1.91256261E+01;
  COFD[149] = 5.02608697E+00;
  COFD[150] = -4.26959993E-01;
  COFD[151] = 1.80709910E-02;
  COFD[152] = -1.41066459E+01;
  COFD[153] = 3.07458927E+00;
  COFD[154] = -1.86899591E-01;
  COFD[155] = 8.19829781E-03;
  COFD[156] = -1.41119732E+01;
  COFD[157] = 3.07458927E+00;
  COFD[158] = -1.86899591E-01;
  COFD[159] = 8.19829781E-03;
  COFD[160] = -1.41170372E+01;
  COFD[161] = 3.07458927E+00;
  COFD[162] = -1.86899591E-01;
  COFD[163] = 8.19829781E-03;
  COFD[164] = -1.38948667E+01;
  COFD[165] = 2.98558426E+00;
  COFD[166] = -1.75507216E-01;
  COFD[167] = 7.71173691E-03;
  COFD[168] = -1.47298720E+01;
  COFD[169] = 3.30683499E+00;
  COFD[170] = -2.16378602E-01;
  COFD[171] = 9.44670561E-03;
  COFD[172] = -9.71375861E+00;
  COFD[173] = 1.77912272E+00;
  COFD[174] = -1.67349571E-02;
  COFD[175] = 7.45446845E-04;
  COFD[176] = -1.93611051E+01;
  COFD[177] = 5.51579726E+00;
  COFD[178] = -4.76061961E-01;
  COFD[179] = 1.96329391E-02;
  COFD[180] = -1.68758926E+01;
  COFD[181] = 4.49460303E+00;
  COFD[182] = -3.64766132E-01;
  COFD[183] = 1.56457153E-02;
  COFD[184] = -1.91096797E+01;
  COFD[185] = 5.02608697E+00;
  COFD[186] = -4.26959993E-01;
  COFD[187] = 1.80709910E-02;
  COFD[188] = -1.91256261E+01;
  COFD[189] = 5.02608697E+00;
  COFD[190] = -4.26959993E-01;
  COFD[191] = 1.80709910E-02;
  COFD[192] = -1.31492641E+01;
  COFD[193] = 1.48004311E+00;
  COFD[194] = 1.60499553E-01;
  COFD[195] = -1.19765679E-02;
  COFD[196] = -2.10640014E+01;
  COFD[197] = 5.50980695E+00;
  COFD[198] = -4.78335488E-01;
  COFD[199] = 1.98515434E-02;
  COFD[200] = -2.04177482E+01;
  COFD[201] = 5.31457079E+00;
  COFD[202] = -4.58216496E-01;
  COFD[203] = 1.91825910E-02;
  COFD[204] = -2.04230073E+01;
  COFD[205] = 5.31457079E+00;
  COFD[206] = -4.58216496E-01;
  COFD[207] = 1.91825910E-02;
  COFD[208] = -2.08123325E+01;
  COFD[209] = 5.42470154E+00;
  COFD[210] = -4.69700416E-01;
  COFD[211] = 1.95706904E-02;
  COFD[212] = -2.10785324E+01;
  COFD[213] = 5.51573149E+00;
  COFD[214] = -4.78177665E-01;
  COFD[215] = 1.98082796E-02;
  COFD[216] = -1.21950642E+01;
  COFD[217] = 2.72222246E+00;
  COFD[218] = -1.41335602E-01;
  COFD[219] = 6.23222872E-03;
  COFD[220] = -1.43712864E+01;
  COFD[221] = 3.70920439E+00;
  COFD[222] = -2.67274113E-01;
  COFD[223] = 1.15967481E-02;
  COFD[224] = -1.15797750E+01;
  COFD[225] = 2.43235504E+00;
  COFD[226] = -1.02890179E-01;
  COFD[227] = 4.52903603E-03;
  COFD[228] = -1.40864894E+01;
  COFD[229] = 3.07458927E+00;
  COFD[230] = -1.86899591E-01;
  COFD[231] = 8.19829781E-03;
  COFD[232] = -1.41066459E+01;
  COFD[233] = 3.07458927E+00;
  COFD[234] = -1.86899591E-01;
  COFD[235] = 8.19829781E-03;
  COFD[236] = -2.10640014E+01;
  COFD[237] = 5.50980695E+00;
  COFD[238] = -4.78335488E-01;
  COFD[239] = 1.98515434E-02;
  COFD[240] = -1.53110708E+01;
  COFD[241] = 3.37317428E+00;
  COFD[242] = -2.24900439E-01;
  COFD[243] = 9.81228151E-03;
  COFD[244] = -1.53187643E+01;
  COFD[245] = 3.37317428E+00;
  COFD[246] = -2.24900439E-01;
  COFD[247] = 9.81228151E-03;
  COFD[248] = -1.53261114E+01;
  COFD[249] = 3.37317428E+00;
  COFD[250] = -2.24900439E-01;
  COFD[251] = 9.81228151E-03;
  COFD[252] = -1.50096240E+01;
  COFD[253] = 3.25515933E+00;
  COFD[254] = -2.09710110E-01;
  COFD[255] = 9.15941830E-03;
  COFD[256] = -1.59592184E+01;
  COFD[257] = 3.60186887E+00;
  COFD[258] = -2.53302622E-01;
  COFD[259] = 1.09893496E-02;
  COFD[260] = -1.03310318E+01;
  COFD[261] = 1.90522472E+00;
  COFD[262] = -3.44812795E-02;
  COFD[263] = 1.57640018E-03;
  COFD[264] = -1.43717529E+01;
  COFD[265] = 3.70920439E+00;
  COFD[266] = -2.67274113E-01;
  COFD[267] = 1.15967481E-02;
  COFD[268] = -1.15806808E+01;
  COFD[269] = 2.43235504E+00;
  COFD[270] = -1.02890179E-01;
  COFD[271] = 4.52903603E-03;
  COFD[272] = -1.40916052E+01;
  COFD[273] = 3.07458927E+00;
  COFD[274] = -1.86899591E-01;
  COFD[275] = 8.19829781E-03;
  COFD[276] = -1.41119732E+01;
  COFD[277] = 3.07458927E+00;
  COFD[278] = -1.86899591E-01;
  COFD[279] = 8.19829781E-03;
  COFD[280] = -2.04177482E+01;
  COFD[281] = 5.31457079E+00;
  COFD[282] = -4.58216496E-01;
  COFD[283] = 1.91825910E-02;
  COFD[284] = -1.53187643E+01;
  COFD[285] = 3.37317428E+00;
  COFD[286] = -2.24900439E-01;
  COFD[287] = 9.81228151E-03;
  COFD[288] = -1.53265780E+01;
  COFD[289] = 3.37317428E+00;
  COFD[290] = -2.24900439E-01;
  COFD[291] = 9.81228151E-03;
  COFD[292] = -1.53340417E+01;
  COFD[293] = 3.37317428E+00;
  COFD[294] = -2.24900439E-01;
  COFD[295] = 9.81228151E-03;
  COFD[296] = -1.50168028E+01;
  COFD[297] = 3.25515933E+00;
  COFD[298] = -2.09710110E-01;
  COFD[299] = 9.15941830E-03;
  COFD[300] = -1.59677692E+01;
  COFD[301] = 3.60186887E+00;
  COFD[302] = -2.53302622E-01;
  COFD[303] = 1.09893496E-02;
  COFD[304] = -1.03327323E+01;
  COFD[305] = 1.90522472E+00;
  COFD[306] = -3.44812795E-02;
  COFD[307] = 1.57640018E-03;
  COFD[308] = -1.43721922E+01;
  COFD[309] = 3.70920439E+00;
  COFD[310] = -2.67274113E-01;
  COFD[311] = 1.15967481E-02;
  COFD[312] = -1.15815344E+01;
  COFD[313] = 2.43235504E+00;
  COFD[314] = -1.02890179E-01;
  COFD[315] = 4.52903603E-03;
  COFD[316] = -1.40964661E+01;
  COFD[317] = 3.07458927E+00;
  COFD[318] = -1.86899591E-01;
  COFD[319] = 8.19829781E-03;
  COFD[320] = -1.41170372E+01;
  COFD[321] = 3.07458927E+00;
  COFD[322] = -1.86899591E-01;
  COFD[323] = 8.19829781E-03;
  COFD[324] = -2.04230073E+01;
  COFD[325] = 5.31457079E+00;
  COFD[326] = -4.58216496E-01;
  COFD[327] = 1.91825910E-02;
  COFD[328] = -1.53261114E+01;
  COFD[329] = 3.37317428E+00;
  COFD[330] = -2.24900439E-01;
  COFD[331] = 9.81228151E-03;
  COFD[332] = -1.53340417E+01;
  COFD[333] = 3.37317428E+00;
  COFD[334] = -2.24900439E-01;
  COFD[335] = 9.81228151E-03;
  COFD[336] = -1.53416186E+01;
  COFD[337] = 3.37317428E+00;
  COFD[338] = -2.24900439E-01;
  COFD[339] = 9.81228151E-03;
  COFD[340] = -1.50236516E+01;
  COFD[341] = 3.25515933E+00;
  COFD[342] = -2.09710110E-01;
  COFD[343] = 9.15941830E-03;
  COFD[344] = -1.59759490E+01;
  COFD[345] = 3.60186887E+00;
  COFD[346] = -2.53302622E-01;
  COFD[347] = 1.09893496E-02;
  COFD[348] = -1.03343373E+01;
  COFD[349] = 1.90522472E+00;
  COFD[350] = -3.44812795E-02;
  COFD[351] = 1.57640018E-03;
  COFD[352] = -1.40298830E+01;
  COFD[353] = 3.55837688E+00;
  COFD[354] = -2.47785790E-01;
  COFD[355] = 1.07555332E-02;
  COFD[356] = -1.13253458E+01;
  COFD[357] = 2.31195095E+00;
  COFD[358] = -8.63988037E-02;
  COFD[359] = 3.77573452E-03;
  COFD[360] = -1.38756407E+01;
  COFD[361] = 2.98558426E+00;
  COFD[362] = -1.75507216E-01;
  COFD[363] = 7.71173691E-03;
  COFD[364] = -1.38948667E+01;
  COFD[365] = 2.98558426E+00;
  COFD[366] = -1.75507216E-01;
  COFD[367] = 7.71173691E-03;
  COFD[368] = -2.08123325E+01;
  COFD[369] = 5.42470154E+00;
  COFD[370] = -4.69700416E-01;
  COFD[371] = 1.95706904E-02;
  COFD[372] = -1.50096240E+01;
  COFD[373] = 3.25515933E+00;
  COFD[374] = -2.09710110E-01;
  COFD[375] = 9.15941830E-03;
  COFD[376] = -1.50168028E+01;
  COFD[377] = 3.25515933E+00;
  COFD[378] = -2.09710110E-01;
  COFD[379] = 9.15941830E-03;
  COFD[380] = -1.50236516E+01;
  COFD[381] = 3.25515933E+00;
  COFD[382] = -2.09710110E-01;
  COFD[383] = 9.15941830E-03;
  COFD[384] = -1.47639290E+01;
  COFD[385] = 3.15955654E+00;
  COFD[386] = -1.97590757E-01;
  COFD[387] = 8.64692156E-03;
  COFD[388] = -1.57236706E+01;
  COFD[389] = 3.51447210E+00;
  COFD[390] = -2.42579007E-01;
  COFD[391] = 1.05506318E-02;
  COFD[392] = -1.01976409E+01;
  COFD[393] = 1.83188320E+00;
  COFD[394] = -2.40547456E-02;
  COFD[395] = 1.08399898E-03;
  COFD[396] = -1.51208119E+01;
  COFD[397] = 3.99904647E+00;
  COFD[398] = -3.03517220E-01;
  COFD[399] = 1.31117363E-02;
  COFD[400] = -1.20638601E+01;
  COFD[401] = 2.63303536E+00;
  COFD[402] = -1.29792632E-01;
  COFD[403] = 5.73363738E-03;
  COFD[404] = -1.47082523E+01;
  COFD[405] = 3.30683499E+00;
  COFD[406] = -2.16378602E-01;
  COFD[407] = 9.44670561E-03;
  COFD[408] = -1.47298720E+01;
  COFD[409] = 3.30683499E+00;
  COFD[410] = -2.16378602E-01;
  COFD[411] = 9.44670561E-03;
  COFD[412] = -2.10785324E+01;
  COFD[413] = 5.51573149E+00;
  COFD[414] = -4.78177665E-01;
  COFD[415] = 1.98082796E-02;
  COFD[416] = -1.59592184E+01;
  COFD[417] = 3.60186887E+00;
  COFD[418] = -2.53302622E-01;
  COFD[419] = 1.09893496E-02;
  COFD[420] = -1.59677692E+01;
  COFD[421] = 3.60186887E+00;
  COFD[422] = -2.53302622E-01;
  COFD[423] = 1.09893496E-02;
  COFD[424] = -1.59759490E+01;
  COFD[425] = 3.60186887E+00;
  COFD[426] = -2.53302622E-01;
  COFD[427] = 1.09893496E-02;
  COFD[428] = -1.57236706E+01;
  COFD[429] = 3.51447210E+00;
  COFD[430] = -2.42579007E-01;
  COFD[431] = 1.05506318E-02;
  COFD[432] = -1.68944722E+01;
  COFD[433] = 3.94346012E+00;
  COFD[434] = -2.96835271E-01;
  COFD[435] = 1.28438696E-02;
  COFD[436] = -1.08140177E+01;
  COFD[437] = 2.11737538E+00;
  COFD[438] = -6.46167749E-02;
  COFD[439] = 2.99827695E-03;
  COFD[440] = -9.71338331E+00;
  COFD[441] = 2.17561180E+00;
  COFD[442] = -7.28270090E-02;
  COFD[443] = 3.38302182E-03;
  COFD[444] = -9.86429034E+00;
  COFD[445] = 2.05348746E+00;
  COFD[446] = -5.90289007E-02;
  COFD[447] = 2.89596157E-03;
  COFD[448] = -9.70779324E+00;
  COFD[449] = 1.77912272E+00;
  COFD[450] = -1.67349571E-02;
  COFD[451] = 7.45446845E-04;
  COFD[452] = -9.71375861E+00;
  COFD[453] = 1.77912272E+00;
  COFD[454] = -1.67349571E-02;
  COFD[455] = 7.45446845E-04;
  COFD[456] = -1.21950642E+01;
  COFD[457] = 2.72222246E+00;
  COFD[458] = -1.41335602E-01;
  COFD[459] = 6.23222872E-03;
  COFD[460] = -1.03310318E+01;
  COFD[461] = 1.90522472E+00;
  COFD[462] = -3.44812795E-02;
  COFD[463] = 1.57640018E-03;
  COFD[464] = -1.03327323E+01;
  COFD[465] = 1.90522472E+00;
  COFD[466] = -3.44812795E-02;
  COFD[467] = 1.57640018E-03;
  COFD[468] = -1.03343373E+01;
  COFD[469] = 1.90522472E+00;
  COFD[470] = -3.44812795E-02;
  COFD[471] = 1.57640018E-03;
  COFD[472] = -1.01976409E+01;
  COFD[473] = 1.83188320E+00;
  COFD[474] = -2.40547456E-02;
  COFD[475] = 1.08399898E-03;
  COFD[476] = -1.08140177E+01;
  COFD[477] = 2.11737538E+00;
  COFD[478] = -6.46167749E-02;
  COFD[479] = 2.99827695E-03;
  COFD[480] = -7.72963289E+00;
  COFD[481] = 1.13864728E+00;
  COFD[482] = 7.22991035E-02;
  COFD[483] = -3.32491895E-03;
}

/*List of specs with small weight, dim NLITE */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKTDIF EGTRANSETKTDIF
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKTDIF egtransetktdif
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKTDIF egtransetktdif_
#endif
void
egtransetKTDIF(int* KTDIF)
{
  KTDIF[0] = 1;
  KTDIF[1] = 2;
  KTDIF[2] = 11;
}

/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFTD EGTRANSETCOFTD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFTD egtransetcoftd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFTD egtransetcoftd_
#endif
void
egtransetCOFTD(double* COFTD)
{
  COFTD[0] = 0.00000000E+00;
  COFTD[1] = 0.00000000E+00;
  COFTD[2] = 0.00000000E+00;
  COFTD[3] = 0.00000000E+00;
  COFTD[4] = 1.52534742E-01;
  COFTD[5] = 5.46404022E-05;
  COFTD[6] = -2.93412470E-08;
  COFTD[7] = 4.87091914E-12;
  COFTD[8] = 2.70010150E-01;
  COFTD[9] = 3.61555093E-04;
  COFTD[10] = -1.80744752E-07;
  COFTD[11] = 2.75321248E-11;
  COFTD[12] = 2.72041664E-01;
  COFTD[13] = 3.64275376E-04;
  COFTD[14] = -1.82104647E-07;
  COFTD[15] = 2.77392722E-11;
  COFTD[16] = -1.41883744E-01;
  COFTD[17] = 7.66558810E-04;
  COFTD[18] = -3.06550003E-07;
  COFTD[19] = 4.02959502E-11;
  COFTD[20] = 2.20482843E-01;
  COFTD[21] = 4.80164288E-04;
  COFTD[22] = -2.32927944E-07;
  COFTD[23] = 3.46470436E-11;
  COFTD[24] = 2.20907853E-01;
  COFTD[25] = 4.81089870E-04;
  COFTD[26] = -2.33376944E-07;
  COFTD[27] = 3.47138305E-11;
  COFTD[28] = 2.21308399E-01;
  COFTD[29] = 4.81962174E-04;
  COFTD[30] = -2.33800100E-07;
  COFTD[31] = 3.47767730E-11;
  COFTD[32] = 2.40744421E-01;
  COFTD[33] = 4.45343451E-04;
  COFTD[34] = -2.18173874E-07;
  COFTD[35] = 3.26958506E-11;
  COFTD[36] = 1.65429221E-01;
  COFTD[37] = 5.61238922E-04;
  COFTD[38] = -2.65650544E-07;
  COFTD[39] = 3.88229592E-11;
  COFTD[40] = 3.40762433E-01;
  COFTD[41] = -4.04057756E-05;
  COFTD[42] = 3.27879533E-08;
  COFTD[43] = -6.27093812E-12;
  COFTD[44] = -1.52534742E-01;
  COFTD[45] = -5.46404022E-05;
  COFTD[46] = 2.93412470E-08;
  COFTD[47] = -4.87091914E-12;
  COFTD[48] = 0.00000000E+00;
  COFTD[49] = 0.00000000E+00;
  COFTD[50] = 0.00000000E+00;
  COFTD[51] = 0.00000000E+00;
  COFTD[52] = 4.15583337E-01;
  COFTD[53] = 1.09738399E-05;
  COFTD[54] = -3.96021963E-09;
  COFTD[55] = 1.14414443E-12;
  COFTD[56] = 4.21932443E-01;
  COFTD[57] = 1.11414935E-05;
  COFTD[58] = -4.02072219E-09;
  COFTD[59] = 1.16162418E-12;
  COFTD[60] = 6.02028221E-02;
  COFTD[61] = 5.61561867E-04;
  COFTD[62] = -2.55372862E-07;
  COFTD[63] = 3.63389913E-11;
  COFTD[64] = 4.42739084E-01;
  COFTD[65] = 7.11770818E-05;
  COFTD[66] = -3.84768062E-08;
  COFTD[67] = 6.86323437E-12;
  COFTD[68] = 4.44452569E-01;
  COFTD[69] = 7.14525507E-05;
  COFTD[70] = -3.86257187E-08;
  COFTD[71] = 6.88979640E-12;
  COFTD[72] = 4.46070183E-01;
  COFTD[73] = 7.17126069E-05;
  COFTD[74] = -3.87662996E-08;
  COFTD[75] = 6.91487226E-12;
  COFTD[76] = 4.45261966E-01;
  COFTD[77] = 4.94697174E-05;
  COFTD[78] = -2.63023442E-08;
  COFTD[79] = 4.90306217E-12;
  COFTD[80] = 4.22530228E-01;
  COFTD[81] = 1.32084268E-04;
  COFTD[82] = -7.12222323E-08;
  COFTD[83] = 1.19516090E-11;
  COFTD[84] = 1.61613664E-01;
  COFTD[85] = 4.74155340E-05;
  COFTD[86] = -1.67115247E-08;
  COFTD[87] = -1.88982125E-12;
  COFTD[88] = -3.40762433E-01;
  COFTD[89] = 4.04057756E-05;
  COFTD[90] = -3.27879533E-08;
  COFTD[91] = 6.27093812E-12;
  COFTD[92] = -1.61613664E-01;
  COFTD[93] = -4.74155340E-05;
  COFTD[94] = 1.67115247E-08;
  COFTD[95] = 1.88982125E-12;
  COFTD[96] = 3.31587939E-01;
  COFTD[97] = -1.96388078E-05;
  COFTD[98] = 3.02388828E-08;
  COFTD[99] = -8.44998018E-12;
  COFTD[100] = 3.42203127E-01;
  COFTD[101] = -2.02675087E-05;
  COFTD[102] = 3.12069259E-08;
  COFTD[103] = -8.72049099E-12;
  COFTD[104] = 2.84983505E-01;
  COFTD[105] = 1.15460005E-04;
  COFTD[106] = -6.17197869E-08;
  COFTD[107] = 1.01504212E-11;
  COFTD[108] = 4.40220831E-01;
  COFTD[109] = -4.83717413E-05;
  COFTD[110] = 4.66088897E-08;
  COFTD[111] = -1.02768430E-11;
  COFTD[112] = 4.43649137E-01;
  COFTD[113] = -4.87484458E-05;
  COFTD[114] = 4.69718656E-08;
  COFTD[115] = -1.03568760E-11;
  COFTD[116] = 4.46895651E-01;
  COFTD[117] = -4.91051748E-05;
  COFTD[118] = 4.73155940E-08;
  COFTD[119] = -1.04326650E-11;
  COFTD[120] = 4.22009934E-01;
  COFTD[121] = -4.14042334E-05;
  COFTD[122] = 4.38751613E-08;
  COFTD[123] = -1.02860246E-11;
  COFTD[124] = 4.66315159E-01;
  COFTD[125] = -5.60150425E-05;
  COFTD[126] = 4.65987669E-08;
  COFTD[127] = -9.13646318E-12;
  COFTD[128] = 0.00000000E+00;
  COFTD[129] = 0.00000000E+00;
  COFTD[130] = 0.00000000E+00;
  COFTD[131] = 0.00000000E+00;
}
}

/* End of file  */

#if 0




\\
\\
\\  This is the mechanism file
\\
\\
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
!
!                     ----- H2 Kinetic Mechanism -----
!                     -----   Version 6-10-2011  -----
!
! (c) Burke, Chaos, Ju, Dryer, and Klippenstein; Princeton University, 2011.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  HOW TO USE THIS MECHANISM:
!
! (*) Due to limitations of CHEMKIN-II format (specifically, an inability to
!     implement temperature-dependent collision efficiencies in falloff
!     reactions) and the lack of fundamental understanding of the mixing rules
!     for the falloff reactions with the bath gases that have different
!     broadening factors, the present implementation represents a compromise
!     (approximate) formulation.  As a consequence,
!
!     PRIOR TO ITS USE IN THE CALCULATIONS, THIS FILE HAS TO BE MODIFIED.
!     DEPENDING ON WHAT BATH GAS (DILUTANT) IS MOST ABUNDANT IN YOUR SYSTEM
!     (THE PRESENT CHOICES ARE N2, AR, OR HE),  YOU  SHOULD UNCOMMENT THE
!     CORRESPONDING BLOCK FOR THE REACTION H+O2(+M)=HO2(+M), AND COMMENT THE
!     BLOCK FOR OTHER DILUTANT(S).  AS GIVEN, THE MAIN DILUTANT IS SET TO BE N2.
!
!
!  HOW TO REFERENCE THIS MECHANISM:
!
!     M.P. Burke, M. Chaos, Y. Ju, F.L. Dryer, S.J. Klippenstein
!        "Comprehensive H2/O2 Kinetic Model for High-Pressure Combustion,"
!        Int. J. Chem. Kinet. (2011).
!
!  FUTURE REVISIONS/UPDATES MAY BE FOUND ON THE FUELS AND COMBUSTION RESEARCH LABORATORY
!  WEBSITE: < http://www.princeton.edu/mae/people/faculty/dryer/homepage/combustion_lab/ >
!
!
!  HOW TO CONTACT THE AUTHORS:
!
!     Dr. Michael P. Burke
!     R122 Building 200
!     Chemical Sciences and Engineering Division
!     Argonne National Laboratory
!     Argonne, IL 60439
!     Email: mpburke@anl.gov
!
!     Prof. Frederick L. Dryer
!     D-329D Engineering Quadrangle
!     Mechanical and Aerospace Engineering
!     Princeton University
!     Princeton, NJ 08544
!     Phone: 609-258-5206
!     Lab:   609-258-0316
!     FAX:   609-258-1939
!     Email: fldryer@princeton.edu
!
!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
!
ELEMENTS
H O N AR HE C
END

SPECIES
H        H2       O        OH
H2O      O2       HO2      H2O2     
N2       AR       HE       
END

!*********************************************************************************

TRANS ALL
H                  0   145.000     2.050     0.000     0.000     0.000          
H2                 1    38.000     2.920     0.000     0.790   280.000          
O                  0    80.000     2.750     0.000     0.000     0.000          
OH                 1    80.000     2.750     0.000     0.000     0.000          
H2O                2   572.400     2.605     1.844     0.000     4.000          
O2                 1   107.400     3.458     0.000     1.600     3.800          
HO2                2   107.400     3.458     0.000     0.000     1.000          
H2O2               2   107.400     3.458     0.000     0.000     3.800          
N2                 1    97.530     3.621     0.000     1.760     4.000          
AR                 0   136.500     3.330     0.000     0.000     0.000          
HE                 0    10.200     2.576     0.000     0.000     0.000          
END

THERMO ALL
0300.00  1000.00  5000.00
H                 120186H   1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 0.02547163E+06-0.04601176E+01 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.02547163E+06-0.04601176E+01                   4
H2                121286H   2               G  0300.00   5000.00  1000.00      1
 0.02991423E+02 0.07000644E-02-0.05633829E-06-0.09231578E-10 0.01582752E-13    2
-0.08350340E+04-0.01355110E+02 0.03298124E+02 0.08249442E-02-0.08143015E-05    3
-0.09475434E-09 0.04134872E-11-0.01012521E+05-0.03294094E+02                   4
O                 120186O   1               G  0300.00   5000.00  1000.00      1
 0.02542060E+02-0.02755062E-03-0.03102803E-07 0.04551067E-10-0.04368052E-14    2
 0.02923080E+06 0.04920308E+02 0.02946429E+02-0.01638166E-01 0.02421032E-04    3
-0.01602843E-07 0.03890696E-11 0.02914764E+06 0.02963995E+02                   4
OH                S 9/01O   1H   1    0    0G   200.000  6000.000 1000.        1
 2.86472886E+00 1.05650448E-03-2.59082758E-07 3.05218674E-11-1.33195876E-15    2
 3.68362875E+03 5.70164073E+00 4.12530561E+00-3.22544939E-03 6.52764691E-06    3
-5.79853643E-09 2.06237379E-12 3.34630913E+03-6.90432960E-01 4.51532273E+03    4
H2O                20387H   2O   1          G  0300.00   5000.00  1000.00      1
 0.02672146E+02 0.03056293E-01-0.08730260E-05 0.01200996E-08-0.06391618E-13    2
-0.02989921E+06 0.06862817E+02 0.03386842E+02 0.03474982E-01-0.06354696E-04    3
 0.06968581E-07-0.02506588E-10-0.03020811E+06 0.02590233E+02                   4
O2                121386O   2               G  0300.00   5000.00  1000.00      1
 0.03697578E+02 0.06135197E-02-0.01258842E-05 0.01775281E-09-0.01136435E-13    2
-0.01233930E+05 0.03189166E+02 0.03212936E+02 0.01127486E-01-0.05756150E-05    3
 0.01313877E-07-0.08768554E-11-0.01005249E+05 0.06034738E+02                   4
HO2               L 5/89H   1O   2   00   00G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00 1.00021620E+04    4
H2O2              120186H   2O   2          G  0300.00   5000.00  1000.00      1
 0.04573167E+02 0.04336136E-01-0.01474689E-04 0.02348904E-08-0.01431654E-12    2
-0.01800696E+06 0.05011370E+01 0.03388754E+02 0.06569226E-01-0.01485013E-05    3
-0.04625806E-07 0.02471515E-10-0.01766315E+06 0.06785363E+02                   4
N2                121286N   2               G  0300.00   5000.00  1000.00      1
 0.02926640E+02 0.01487977E-01-0.05684761E-05 0.01009704E-08-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.01408240E-01-0.03963222E-04    3
 0.05641515E-07-0.02444855E-10-0.01020900E+05 0.03950372E+02                   4
AR                120186AR  1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.04366001E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366001E+02                   4
HE                120186HE  1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.09153489E+01 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.09153488E+01                   4
END

!*********************************************************************************

REACTIONS

!======================
!H2-O2 Chain Reactions
!======================

! Hong et al., Proc. Comb. Inst. 33:309-316 (2011)
H+O2 = O+OH                                 	1.04E+14   0.00  1.5286E+04

! Baulch et al., J. Phys. Chem. Ref. Data, 21:411 (1992)
O+H2 = H+OH						3.818E+12  0.00  7.948E+03
   DUPLICATE
O+H2 = H+OH						8.792E+14  0.00  1.917E+04
   DUPLICATE

! Michael and Sutherland, J. Phys. Chem. 92:3853 (1988)
H2+OH = H2O+H						0.216E+09  1.51  0.343E+04

! Baulch et al., J. Phys. Chem. Ref. Data, 21:411 (1992)
OH+OH = O+H2O						3.34E+04   2.42  -1.93E+03

!============================
!H2-O2 Dissociation Reactions
!============================

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) 
H2+M = H+H+M						4.577E+19 -1.40  1.0438E+05
   H2/2.5/ H2O/12/
   AR/0.0/ HE/0.0/

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) 
H2+AR = H+H+AR                              	5.840E+18 -1.10  1.0438E+05
H2+HE = H+H+HE                              	5.840E+18 -1.10  1.0438E+05

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) 
O+O+M = O2+M                                	6.165E+15 -0.50  0.000E+00
   H2/2.5/ H2O/12/
   AR/0.0/ HE/0.0/

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) 
O+O+AR = O2+AR                              	1.886E+13  0.00 -1.788E+03
O+O+HE = O2+HE                              	1.886E+13  0.00 -1.788E+03

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) 
O+H+M = OH+M                                	4.714E+18 -1.00  0.000E+00
   H2/2.5/  H2O/12/
   AR/0.75/ HE/0.75/

! Srinivasan and Michael, Int. J. Chem. Kinetic. 38 (2006)
! Rate constant is for Ar with efficiencies from Michael et al., J. Phys. Chem. A, 106 (2002)
H2O+M = H+OH+M                              	6.064E+27 -3.322 1.2079E+05
   H2/3.0/  H2O/0.0/
   HE/1.10/ N2/2.00/
   O2/1.5/
! Efficiencies for CO and CO2 taken from Li et al., Int. J. Chem. Kinet. 36:566-575 (2004)

! Srinivasan and Michael, Int. J. Chem. Kinetic. 38 (2006)
H2O+H2O = H+OH+H2O                          	1.006E+26 -2.44  1.2018E+05

!=================================
! Formation and consumption of HO2
!=================================

! High-pressure limit from Troe, Proc. Comb. Inst. 28:1463-1469 (2000)
! Low-pressure  limit from Michael et al., J. Phys. Chem. A 106:5297-5313
! Centering factors from Fernandes et al., Phys. Chem. Chem. Phys. 10:4313-4321 (2008)
!=================================================================================
! MAIN BATH GAS IS N2 (comment this reaction otherwise)
!
H+O2(+M) = HO2(+M)                          	4.65084E+12  0.44  0.000E+00
   LOW/6.366E+20 -1.72  5.248E+02/
   TROE/0.5  1E-30  1E+30/
   H2/2.0/ H2O/14/ O2/0.78/ AR/0.67/ HE/0.8/
!=================================================================================
! MAIN BATH GAS IS AR OR HE (comment this reaction otherwise)
!
!H+O2(+M) = HO2(+M)                         	4.65084E+12  0.44  0.000E+00
!   LOW/9.042E+19 -1.50  4.922E+02/
!   TROE/0.5 1E-30  1E+30/
!   H2/3.0/ H2O/21/ O2/1.1/ CO/2.7/ CO2/5.4/ HE/1.2/ N2/1.5/
!=================================================================================

! Michael et al., Proc. Comb. Inst. 28:1471 (2000)
!HO2+H = H2+O2                                 	3.659E+06  2.09 -1.451E+03
!Scaled by 0.75
HO2+H = H2+O2                                 	2.750E+06  2.09 -1.451E+03

! Mueller et al., Int. J. Chem. Kinetic. 31:113 (1999) 
HO2+H = OH+OH                               	7.079E+13  0.00  2.950E+02

! Fernandez-Ramos and Varandas, J. Phys. Chem. A 106:4077-4083 (2002)
!HO2+O = O2+OH                               	4.750E+10  1.00 -7.2393E+02
!Scaled by 0.60
HO2+O = O2+OH                               	2.850E+10  1.00 -7.2393E+02

! Keyser, J. Phys. Chem. 92:1193 (1988)
HO2+OH = H2O+O2                             	2.890E+13  0.00 -4.970E+02

!=====================================
!Formation and Consumption of H2O2
!=====================================

! Hippler et al., J. Chem. Phys. 93:1755 (1990)
HO2+HO2 = H2O2+O2                           	4.200E+14  0.00  1.1982E+04
   DUPLICATE
HO2+HO2 = H2O2+O2                           	1.300E+11  0.00 -1.6293E+03
   DUPLICATE

! Troe, Combust. Flame,  158:594-601 (2011)
! Rate constant is for Ar
H2O2(+M) = OH+OH(+M)            			2.00E+12   0.90  4.8749E+04
   LOW/2.49E+24 -2.30 4.8749E+04/
   TROE/0.43 1E-30 1E+30/
   H2O/7.5/ 
   N2/1.5/  O2/1.2/
   HE/0.65/ H2O2/7.7/
! Efficiencies for H2 and CO taken from Li et al., Int. J. Chem. Kinet. 36:566-575 (2004)
   H2/3.7/ 

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2O2+H = H2O+OH                             	2.410E+13  0.00  3.970E+03
H2O2+H = HO2+H2                             	4.820E+13  0.00  7.950E+03
H2O2+O = OH+HO2                             	9.550E+06  2.00  3.970E+03

! Hong et al., J. Phys. Chem. A  114 (2010) 5718-5727
H2O2+OH = HO2+H2O                           	1.740E+12  0.00  3.180E+02
   DUPLICATE
H2O2+OH = HO2+H2O                           	7.590E+13  0.00  7.270E+03
   DUPLICATE

!  JBB added reactions  X1 and X6 from paper

HO2+H = O+H2O                                   3.970E+12  0.00  6.710E+02

O+OH+M = HO2+M                                  8.000E+15  0.00  0.000E+00
   H2/2.0/  H2O/12./ AR/0.7/  HE/0.7/  



END

\\
\\
\\  This is the therm file
\\
\\
THERMO ALL
300.0 1000.0 5000.0
H                 120186H   1               G  0300.00   5000.00  1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 2.54716270E+04-4.60117638E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 2.54716270E+04-4.60117608E-01                   4
O                 120186O   1               G  0300.00   5000.00  1000.00      1
 2.54205966E+00-2.75506191E-05-3.10280335E-09 4.55106742E-12-4.36805150E-16    2
 2.92308027E+04 4.92030811E+00 2.94642878E+00-1.63816649E-03 2.42103170E-06    3
-1.60284319E-09 3.89069636E-13 2.91476445E+04 2.96399498E+00                   4
OH                S 9/01O   1H   1    0    0G   200.000  6000.000 1000.        1
 2.86472886E+00 1.05650448E-03-2.59082758E-07 3.05218674E-11-1.33195876E-15    2
 3.68362875E+03 5.70164073E+00 4.12530561E+00-3.22544939E-03 6.52764691E-06    3
-5.79853643E-09 2.06237379E-12 3.34630913E+03-6.90432960E-01 4.51532273E+03    4
H2                121286H   2               G  0300.00   5000.00  1000.00      1
 2.99142337E+00 7.00064411E-04-5.63382869E-08-9.23157818E-12 1.58275179E-15    2
-8.35033997E+02-1.35511017E+00 3.29812431E+00 8.24944174E-04-8.14301529E-07    3
-9.47543433E-11 4.13487224E-13-1.01252087E+03-3.29409409E+00                   4
O2                121386O   2               G  0300.00   5000.00  1000.00      1
 3.69757819E+00 6.13519689E-04-1.25884199E-07 1.77528148E-11-1.13643531E-15    2
-1.23393018E+03 3.18916559E+00 3.21293640E+00 1.12748635E-03-5.75615047E-07    3
 1.31387723E-09-8.76855392E-13-1.00524902E+03 6.03473759E+00                   4
H2O                20387H   2O   1          G  0300.00   5000.00  1000.00      1
 2.67214561E+00 3.05629289E-03-8.73026011E-07 1.20099639E-10-6.39161787E-15    2
-2.98992090E+04 6.86281681E+00 3.38684249E+00 3.47498246E-03-6.35469633E-06    3
 6.96858127E-09-2.50658847E-12-3.02081133E+04 2.59023285E+00                   4
HO2               L 5/89H   1O   2   00   00G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00 1.00021620E+04    4
H2O2              120186H   2O   2          G  0300.00   5000.00  1000.00      1
 4.57316685E+00 4.33613639E-03-1.47468882E-06 2.34890357E-10-1.43165356E-14    2
-1.80069609E+04 5.01136959E-01 3.38875365E+00 6.56922581E-03-1.48501258E-07    3
-4.62580552E-09 2.47151475E-12-1.76631465E+04 6.78536320E+00                   4
N2                121286N   2               G  0300.00   5000.00  1000.00      1
 0.02926640E+02 0.01487977E-01-0.05684761E-05 0.01009704E-08-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.01408240E-01-0.03963222E-04    3
 0.05641515E-07-0.02444855E-10-0.01020900E+05 0.03950372E+02                   4
END

#endif
