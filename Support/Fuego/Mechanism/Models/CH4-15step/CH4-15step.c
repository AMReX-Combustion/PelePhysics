
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
/* #define CKWC CKWC */
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
/* #define CKWC ckwc */
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
/* #define CKWC ckwc_ */
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
void egtransetEPS(double *  EPS);
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
void egtransetSIG(double* SIG);
void atomicWeight(double * restrict awt);
void molecularWeight(double * restrict wt);
void gibbs(double * restrict species, double * restrict tc);
void helmholtz(double * restrict species, double * restrict tc);
void speciesInternalEnergy(double * restrict species, double * restrict tc);
void speciesEnthalpy(double * restrict species, double * restrict tc);
void speciesEntropy(double * restrict species, double * restrict tc);
void cp_R(double * restrict species, double * restrict tc);
void cv_R(double * restrict species, double * restrict tc);
void equilibriumConstants(double * restrict kc, double * restrict g_RT, double T);
void productionRate(double * restrict wdot, double * restrict sc, double T);
void comp_k_f(double * restrict tc, double invT, double * restrict k_f);
void comp_Kc(double * restrict tc, double invT, double * restrict Kc);
void comp_qfqr(double * restrict q_f, double * restrict q_r, double * restrict sc, double * restrict tc, double invT);
void progressRate(double * restrict qdot, double * restrict speciesConc, double T);
void progressRateFR(double * restrict q_f, double * restrict q_r, double * restrict speciesConc, double T);
void CKINIT();
void CKFINALIZE();
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict rval, int * kerr, int lenline, int lenkray);
void CKSYME(int * kname, int * lenkname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(int * ickwrk, double * restrict rckwrk, double * restrict ru, double * restrict ruc, double * restrict pa);
void CKPX(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict P);
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P);
void CKPC(double * restrict rho, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict P);
void CKRHOX(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKRHOC(double * restrict P, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKWT(int * iwrk, double * restrict rwrk, double * restrict wt);
void CKAWT(int * iwrk, double * restrict rwrk, double * restrict awt);
void CKMMWY(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKMMWX(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKMMWC(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKYTX(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x);
void CKYTCP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c);
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c);
void CKXTY(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y);
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c);
void CKXTCR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c);
void CKCTX(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict x);
void CKCTY(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict y);
void CKCPOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpor);
void CKHORT(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hort);
void CKSOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sor);
void CKCVML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml);
void CKCPML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml);
void CKUML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml);
void CKHML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml);
void CKGML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gml);
void CKAML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict aml);
void CKSML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sml);
void CKCVMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms);
void CKCPMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms);
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void CKHMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void CKGMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gms);
void CKAMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ams);
void CKSMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sms);
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl);
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs);
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl);
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs);
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml);
void CKHBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict hbms);
void CKUBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict ubml);
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms);
void CKSBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict sbml);
void CKSBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict sbms);
void CKGBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict gbml);
void CKGBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict gbms);
void CKABML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict abml);
void CKABMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict abms);
/* void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot); */
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r);
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKNU(int * kdim, int * iwrk, double * restrict rwrk, int * nuki);
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf);
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e );
void CKEQC(double * restrict T, double * restrict C , int * iwrk, double * restrict rwrk, double * restrict eqcon );
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void DWDOT(double * restrict J, double * restrict sc, double * restrict T, int * consP);
void aJacobian(double * restrict J, double * restrict sc, double T, int consP);
void dcvpRdT(double * restrict species, double * restrict tc);
void GET_T_GIVEN_EY(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int *ierr);
void GET_T_GIVEN_HY(double * restrict h, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int *ierr);
void GET_REACTION_MAP(int * restrict rmap);
/*vector version */
void vproductionRate(int npt, double * restrict wdot, double * restrict c, double * restrict T);
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void VCKPY(int * restrict np, double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P);
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
            double * restrict y, int * restrict iwrk, double * restrict rwrk,
            double * restrict wdot);
void VCKYTX(int * restrict np, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x);
void vcomp_k_f(int npt, double * restrict k_f_s, double * restrict tc, double * restrict invT);
void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc);
void vcomp_Kc(int npt, double * restrict Kc_s, double * restrict g_RT, double * restrict invT);
void GET_CRITPARAMS(double * restrict Tci, double * restrict ai, double * restrict bi, double * restrict acentric_i);
void vcomp_wdot(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
                double * restrict k_f_s, double * restrict Kc_s,
                double * restrict tc, double * restrict invT, double * restrict T);

/* Inverse molecular weights */
static const double imw[19] = {
    1.0 / 2.015940,  /*H2 */
    1.0 / 1.007970,  /*H */
    1.0 / 31.998800,  /*O2 */
    1.0 / 17.007370,  /*OH */
    1.0 / 18.015340,  /*H2O */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 15.035060,  /*CH3 */
    1.0 / 16.043030,  /*CH4 */
    1.0 / 28.010550,  /*CO */
    1.0 / 44.009950,  /*CO2 */
    1.0 / 30.026490,  /*CH2O */
    1.0 / 26.038240,  /*C2H2 */
    1.0 / 28.054180,  /*C2H4 */
    1.0 / 30.070120,  /*C2H6 */
    1.0 / 17.030610,  /*NH3 */
    1.0 / 30.006100,  /*NO */
    1.0 / 27.025820,  /*HCN */
    1.0 / 28.013400};  /*N2 */



static double fwd_A[0], fwd_beta[0], fwd_Ea[0];
static double low_A[0], low_beta[0], low_Ea[0];
static double rev_A[0], rev_beta[0], rev_Ea[0];
static double troe_a[0],troe_Ts[0], troe_Tss[0], troe_Tsss[0];
static double sri_a[0], sri_b[0], sri_c[0], sri_d[0], sri_e[0];
static double activation_units[0], prefactor_units[0], phase_units[0];
static int is_PD[0], troe_len[0], sri_len[0], nTB[0], *TBid[0];
static double *TB[0];

static double fwd_A_DEF[0], fwd_beta_DEF[0], fwd_Ea_DEF[0];
static double low_A_DEF[0], low_beta_DEF[0], low_Ea_DEF[0];
static double rev_A_DEF[0], rev_beta_DEF[0], rev_Ea_DEF[0];
static double troe_a_DEF[0],troe_Ts_DEF[0], troe_Tss_DEF[0], troe_Tsss_DEF[0];
static double sri_a_DEF[0], sri_b_DEF[0], sri_c_DEF[0], sri_d_DEF[0], sri_e_DEF[0];
static double activation_units_DEF[0], prefactor_units_DEF[0], phase_units_DEF[0];
static int is_PD_DEF[0], troe_len_DEF[0], sri_len_DEF[0], nTB_DEF[0], *TBid_DEF[0];
static double *TB_DEF[0];
static int rxn_map[0] = {};

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<0; ++i) {
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
  if (reaction_id<0 || reaction_id>=0) {
    printf("Bad reaction id = %d",reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id<0 || species_id>=19) {
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
    for (int i=0; i<0; i++) {
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
    for (int i=0; i<0; i++) {
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
  for (int i=0; i<0; ++i) {
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
    SetAllDefaults();
}



/*A few mechanism parameters */
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 4;
    *kk = 19;
    *ii = 0;
    *nfit = -1; /*Why do you need this anyway ?  */
}



/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict rval, int * kerr, int lenline )
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
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict rval, int * kerr, int lenline, int lenkray)
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

    /* O  */
    kname[ 0*lenkname + 0 ] = 'O';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* C  */
    kname[ 2*lenkname + 0 ] = 'C';
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
    for (i=0; i<lenkname*19; i++) {
        kname[i] = ' ';
    }

    /* H2  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* O2  */
    kname[ 2*lenkname + 0 ] = 'O';
    kname[ 2*lenkname + 1 ] = '2';
    kname[ 2*lenkname + 2 ] = ' ';

    /* OH  */
    kname[ 3*lenkname + 0 ] = 'O';
    kname[ 3*lenkname + 1 ] = 'H';
    kname[ 3*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 4*lenkname + 0 ] = 'H';
    kname[ 4*lenkname + 1 ] = '2';
    kname[ 4*lenkname + 2 ] = 'O';
    kname[ 4*lenkname + 3 ] = ' ';

    /* HO2  */
    kname[ 5*lenkname + 0 ] = 'H';
    kname[ 5*lenkname + 1 ] = 'O';
    kname[ 5*lenkname + 2 ] = '2';
    kname[ 5*lenkname + 3 ] = ' ';

    /* H2O2  */
    kname[ 6*lenkname + 0 ] = 'H';
    kname[ 6*lenkname + 1 ] = '2';
    kname[ 6*lenkname + 2 ] = 'O';
    kname[ 6*lenkname + 3 ] = '2';
    kname[ 6*lenkname + 4 ] = ' ';

    /* CH3  */
    kname[ 7*lenkname + 0 ] = 'C';
    kname[ 7*lenkname + 1 ] = 'H';
    kname[ 7*lenkname + 2 ] = '3';
    kname[ 7*lenkname + 3 ] = ' ';

    /* CH4  */
    kname[ 8*lenkname + 0 ] = 'C';
    kname[ 8*lenkname + 1 ] = 'H';
    kname[ 8*lenkname + 2 ] = '4';
    kname[ 8*lenkname + 3 ] = ' ';

    /* CO  */
    kname[ 9*lenkname + 0 ] = 'C';
    kname[ 9*lenkname + 1 ] = 'O';
    kname[ 9*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 10*lenkname + 0 ] = 'C';
    kname[ 10*lenkname + 1 ] = 'O';
    kname[ 10*lenkname + 2 ] = '2';
    kname[ 10*lenkname + 3 ] = ' ';

    /* CH2O  */
    kname[ 11*lenkname + 0 ] = 'C';
    kname[ 11*lenkname + 1 ] = 'H';
    kname[ 11*lenkname + 2 ] = '2';
    kname[ 11*lenkname + 3 ] = 'O';
    kname[ 11*lenkname + 4 ] = ' ';

    /* C2H2  */
    kname[ 12*lenkname + 0 ] = 'C';
    kname[ 12*lenkname + 1 ] = '2';
    kname[ 12*lenkname + 2 ] = 'H';
    kname[ 12*lenkname + 3 ] = '2';
    kname[ 12*lenkname + 4 ] = ' ';

    /* C2H4  */
    kname[ 13*lenkname + 0 ] = 'C';
    kname[ 13*lenkname + 1 ] = '2';
    kname[ 13*lenkname + 2 ] = 'H';
    kname[ 13*lenkname + 3 ] = '4';
    kname[ 13*lenkname + 4 ] = ' ';

    /* C2H6  */
    kname[ 14*lenkname + 0 ] = 'C';
    kname[ 14*lenkname + 1 ] = '2';
    kname[ 14*lenkname + 2 ] = 'H';
    kname[ 14*lenkname + 3 ] = '6';
    kname[ 14*lenkname + 4 ] = ' ';

    /* NH3  */
    kname[ 15*lenkname + 0 ] = 'N';
    kname[ 15*lenkname + 1 ] = 'H';
    kname[ 15*lenkname + 2 ] = '3';
    kname[ 15*lenkname + 3 ] = ' ';

    /* NO  */
    kname[ 16*lenkname + 0 ] = 'N';
    kname[ 16*lenkname + 1 ] = 'O';
    kname[ 16*lenkname + 2 ] = ' ';

    /* HCN  */
    kname[ 17*lenkname + 0 ] = 'H';
    kname[ 17*lenkname + 1 ] = 'C';
    kname[ 17*lenkname + 2 ] = 'N';
    kname[ 17*lenkname + 3 ] = ' ';

    /* N2  */
    kname[ 18*lenkname + 0 ] = 'N';
    kname[ 18*lenkname + 1 ] = '2';
    kname[ 18*lenkname + 2 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(int * ickwrk, double * restrict rckwrk, double * restrict ru, double * restrict ruc, double * restrict pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*31.998800; /*O2 */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*15.035060; /*CH3 */
    XW += x[8]*16.043030; /*CH4 */
    XW += x[9]*28.010550; /*CO */
    XW += x[10]*44.009950; /*CO2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*26.038240; /*C2H2 */
    XW += x[13]*28.054180; /*C2H4 */
    XW += x[14]*30.070120; /*C2H6 */
    XW += x[15]*17.030610; /*NH3 */
    XW += x[16]*30.006100; /*NO */
    XW += x[17]*27.025820; /*HCN */
    XW += x[18]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O2 */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*HO2 */
    YOW += y[6]*imw[6]; /*H2O2 */
    YOW += y[7]*imw[7]; /*CH3 */
    YOW += y[8]*imw[8]; /*CH4 */
    YOW += y[9]*imw[9]; /*CO */
    YOW += y[10]*imw[10]; /*CO2 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*C2H2 */
    YOW += y[13]*imw[13]; /*C2H4 */
    YOW += y[14]*imw[14]; /*C2H6 */
    YOW += y[15]*imw[15]; /*NH3 */
    YOW += y[16]*imw[16]; /*NO */
    YOW += y[17]*imw[17]; /*HCN */
    YOW += y[18]*imw[18]; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int * restrict np, double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<19; n++) {
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
void CKPC(double * restrict rho, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*31.998800; /*O2 */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*33.006770; /*HO2 */
    W += c[6]*34.014740; /*H2O2 */
    W += c[7]*15.035060; /*CH3 */
    W += c[8]*16.043030; /*CH4 */
    W += c[9]*28.010550; /*CO */
    W += c[10]*44.009950; /*CO2 */
    W += c[11]*30.026490; /*CH2O */
    W += c[12]*26.038240; /*C2H2 */
    W += c[13]*28.054180; /*C2H4 */
    W += c[14]*30.070120; /*C2H6 */
    W += c[15]*17.030610; /*NH3 */
    W += c[16]*30.006100; /*NO */
    W += c[17]*27.025820; /*HCN */
    W += c[18]*28.013400; /*N2 */

    for (id = 0; id < 19; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*31.998800; /*O2 */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*15.035060; /*CH3 */
    XW += x[8]*16.043030; /*CH4 */
    XW += x[9]*28.010550; /*CO */
    XW += x[10]*44.009950; /*CO2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*26.038240; /*C2H2 */
    XW += x[13]*28.054180; /*C2H4 */
    XW += x[14]*30.070120; /*C2H6 */
    XW += x[15]*17.030610; /*NH3 */
    XW += x[16]*30.006100; /*NO */
    XW += x[17]*27.025820; /*HCN */
    XW += x[18]*28.013400; /*N2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double YOW = 0;
    double tmp[19];

    for (int i = 0; i < 19; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 19; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31451e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double * restrict P, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*31.998800; /*O2 */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*33.006770; /*HO2 */
    W += c[6]*34.014740; /*H2O2 */
    W += c[7]*15.035060; /*CH3 */
    W += c[8]*16.043030; /*CH4 */
    W += c[9]*28.010550; /*CO */
    W += c[10]*44.009950; /*CO2 */
    W += c[11]*30.026490; /*CH2O */
    W += c[12]*26.038240; /*C2H2 */
    W += c[13]*28.054180; /*C2H4 */
    W += c[14]*30.070120; /*C2H6 */
    W += c[15]*17.030610; /*NH3 */
    W += c[16]*30.006100; /*NO */
    W += c[17]*27.025820; /*HCN */
    W += c[18]*28.013400; /*N2 */

    for (id = 0; id < 19; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31451e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT(int * iwrk, double * restrict rwrk, double * restrict wt)
{
    molecularWeight(wt);
}


/*get atomic weight for all elements */
void CKAWT(int * iwrk, double * restrict rwrk, double * restrict awt)
{
    atomicWeight(awt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWY(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    double YOW = 0;
    double tmp[19];

    for (int i = 0; i < 19; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 19; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*31.998800; /*O2 */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*15.035060; /*CH3 */
    XW += x[8]*16.043030; /*CH4 */
    XW += x[9]*28.010550; /*CO */
    XW += x[10]*44.009950; /*CO2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*26.038240; /*C2H2 */
    XW += x[13]*28.054180; /*C2H4 */
    XW += x[14]*30.070120; /*C2H6 */
    XW += x[15]*17.030610; /*NH3 */
    XW += x[16]*30.006100; /*NO */
    XW += x[17]*27.025820; /*HCN */
    XW += x[18]*28.013400; /*N2 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void CKMMWC(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*31.998800; /*O2 */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*33.006770; /*HO2 */
    W += c[6]*34.014740; /*H2O2 */
    W += c[7]*15.035060; /*CH3 */
    W += c[8]*16.043030; /*CH4 */
    W += c[9]*28.010550; /*CO */
    W += c[10]*44.009950; /*CO2 */
    W += c[11]*30.026490; /*CH2O */
    W += c[12]*26.038240; /*C2H2 */
    W += c[13]*28.054180; /*C2H4 */
    W += c[14]*30.070120; /*C2H6 */
    W += c[15]*17.030610; /*NH3 */
    W += c[16]*30.006100; /*NO */
    W += c[17]*27.025820; /*HCN */
    W += c[18]*28.013400; /*N2 */

    for (id = 0; id < 19; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
void CKYTX(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)
{
    double YOW = 0;
    double tmp[19];

    for (int i = 0; i < 19; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 19; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 19; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int * restrict np, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<19; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<19; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    double YOW = 0;
    double PWORT;

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 19; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 19; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 19; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    for (int i = 0; i < 19; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*31.998800; /*O2 */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*15.035060; /*CH3 */
    XW += x[8]*16.043030; /*CH4 */
    XW += x[9]*28.010550; /*CO */
    XW += x[10]*44.009950; /*CO2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*26.038240; /*C2H2 */
    XW += x[13]*28.054180; /*C2H4 */
    XW += x[14]*30.070120; /*C2H6 */
    XW += x[15]*17.030610; /*NH3 */
    XW += x[16]*30.006100; /*NO */
    XW += x[17]*27.025820; /*HCN */
    XW += x[18]*28.013400; /*N2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*2.015940*XWinv; 
    y[1] = x[1]*1.007970*XWinv; 
    y[2] = x[2]*31.998800*XWinv; 
    y[3] = x[3]*17.007370*XWinv; 
    y[4] = x[4]*18.015340*XWinv; 
    y[5] = x[5]*33.006770*XWinv; 
    y[6] = x[6]*34.014740*XWinv; 
    y[7] = x[7]*15.035060*XWinv; 
    y[8] = x[8]*16.043030*XWinv; 
    y[9] = x[9]*28.010550*XWinv; 
    y[10] = x[10]*44.009950*XWinv; 
    y[11] = x[11]*30.026490*XWinv; 
    y[12] = x[12]*26.038240*XWinv; 
    y[13] = x[13]*28.054180*XWinv; 
    y[14] = x[14]*30.070120*XWinv; 
    y[15] = x[15]*17.030610*XWinv; 
    y[16] = x[16]*30.006100*XWinv; 
    y[17] = x[17]*27.025820*XWinv; 
    y[18] = x[18]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 19; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*31.998800; /*O2 */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*15.035060; /*CH3 */
    XW += x[8]*16.043030; /*CH4 */
    XW += x[9]*28.010550; /*CO */
    XW += x[10]*44.009950; /*CO2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*26.038240; /*C2H2 */
    XW += x[13]*28.054180; /*C2H4 */
    XW += x[14]*30.070120; /*C2H6 */
    XW += x[15]*17.030610; /*NH3 */
    XW += x[16]*30.006100; /*NO */
    XW += x[17]*27.025820; /*HCN */
    XW += x[18]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 19; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
void CKCTX(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 19; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 19; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*2.015940; /*H2 */
    CW += c[1]*1.007970; /*H */
    CW += c[2]*31.998800; /*O2 */
    CW += c[3]*17.007370; /*OH */
    CW += c[4]*18.015340; /*H2O */
    CW += c[5]*33.006770; /*HO2 */
    CW += c[6]*34.014740; /*H2O2 */
    CW += c[7]*15.035060; /*CH3 */
    CW += c[8]*16.043030; /*CH4 */
    CW += c[9]*28.010550; /*CO */
    CW += c[10]*44.009950; /*CO2 */
    CW += c[11]*30.026490; /*CH2O */
    CW += c[12]*26.038240; /*C2H2 */
    CW += c[13]*28.054180; /*C2H4 */
    CW += c[14]*30.070120; /*C2H6 */
    CW += c[15]*17.030610; /*NH3 */
    CW += c[16]*30.006100; /*NO */
    CW += c[17]*27.025820; /*HCN */
    CW += c[18]*28.013400; /*N2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*2.015940*CWinv; 
    y[1] = c[1]*1.007970*CWinv; 
    y[2] = c[2]*31.998800*CWinv; 
    y[3] = c[3]*17.007370*CWinv; 
    y[4] = c[4]*18.015340*CWinv; 
    y[5] = c[5]*33.006770*CWinv; 
    y[6] = c[6]*34.014740*CWinv; 
    y[7] = c[7]*15.035060*CWinv; 
    y[8] = c[8]*16.043030*CWinv; 
    y[9] = c[9]*28.010550*CWinv; 
    y[10] = c[10]*44.009950*CWinv; 
    y[11] = c[11]*30.026490*CWinv; 
    y[12] = c[12]*26.038240*CWinv; 
    y[13] = c[13]*28.054180*CWinv; 
    y[14] = c[14]*30.070120*CWinv; 
    y[15] = c[15]*17.030610*CWinv; 
    y[16] = c[16]*30.006100*CWinv; 
    y[17] = c[17]*27.025820*CWinv; 
    y[18] = c[18]*28.013400*CWinv; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void CKCPOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void CKHORT(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void CKSOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void CKCVML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 19; ++id) {
        cvml[id] *= 8.31451e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void CKCPML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 19; ++id) {
        cpml[id] *= 8.31451e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 19; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void CKHML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 19; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void CKGML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 19; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void CKAML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 19; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void CKSML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 19; ++id) {
        sml[id] *= 8.31451e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void CKCVMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 4.124383662212169e+07; /*H2 */
    cvms[1] *= 8.248767324424338e+07; /*H */
    cvms[2] *= 2.598381814318037e+06; /*O2 */
    cvms[3] *= 4.888768810227566e+06; /*OH */
    cvms[4] *= 4.615239012974499e+06; /*H2O */
    cvms[5] *= 2.519031701678171e+06; /*HO2 */
    cvms[6] *= 2.444384405113783e+06; /*H2O2 */
    cvms[7] *= 5.530081023953346e+06; /*CH3 */
    cvms[8] *= 5.182630712527496e+06; /*CH4 */
    cvms[9] *= 2.968349425484326e+06; /*CO */
    cvms[10] *= 1.889234139098090e+06; /*CO2 */
    cvms[11] *= 2.769058254894261e+06; /*CH2O */
    cvms[12] *= 3.193192012977835e+06; /*C2H2 */
    cvms[13] *= 2.963733033722604e+06; /*C2H4 */
    cvms[14] *= 2.765040511976673e+06; /*C2H6 */
    cvms[15] *= 4.882097587813943e+06; /*NH3 */
    cvms[16] *= 2.770939908885194e+06; /*NO */
    cvms[17] *= 3.076506096762281e+06; /*HCN */
    cvms[18] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124383662212169e+07; /*H2 */
    cpms[1] *= 8.248767324424338e+07; /*H */
    cpms[2] *= 2.598381814318037e+06; /*O2 */
    cpms[3] *= 4.888768810227566e+06; /*OH */
    cpms[4] *= 4.615239012974499e+06; /*H2O */
    cpms[5] *= 2.519031701678171e+06; /*HO2 */
    cpms[6] *= 2.444384405113783e+06; /*H2O2 */
    cpms[7] *= 5.530081023953346e+06; /*CH3 */
    cpms[8] *= 5.182630712527496e+06; /*CH4 */
    cpms[9] *= 2.968349425484326e+06; /*CO */
    cpms[10] *= 1.889234139098090e+06; /*CO2 */
    cpms[11] *= 2.769058254894261e+06; /*CH2O */
    cpms[12] *= 3.193192012977835e+06; /*C2H2 */
    cpms[13] *= 2.963733033722604e+06; /*C2H4 */
    cpms[14] *= 2.765040511976673e+06; /*C2H6 */
    cpms[15] *= 4.882097587813943e+06; /*NH3 */
    cpms[16] *= 2.770939908885194e+06; /*NO */
    cpms[17] *= 3.076506096762281e+06; /*HCN */
    cpms[18] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 19; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 19; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tc[5], h[19];

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
        hms[12*(*np)+i] = h[12];
        hms[13*(*np)+i] = h[13];
        hms[14*(*np)+i] = h[14];
        hms[15*(*np)+i] = h[15];
        hms[16*(*np)+i] = h[16];
        hms[17*(*np)+i] = h[17];
        hms[18*(*np)+i] = h[18];
    }

    for (int n=0; n<19; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31451e+07 * T[i] * imw[n];
        }
    }
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gms, tc);
    for (int i = 0; i < 19; i++)
    {
        gms[i] *= RT*imw[i];
    }
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(ams, tc);
    for (int i = 0; i < 19; i++)
    {
        ams[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 4.124383662212169e+07; /*H2 */
    sms[1] *= 8.248767324424338e+07; /*H */
    sms[2] *= 2.598381814318037e+06; /*O2 */
    sms[3] *= 4.888768810227566e+06; /*OH */
    sms[4] *= 4.615239012974499e+06; /*H2O */
    sms[5] *= 2.519031701678171e+06; /*HO2 */
    sms[6] *= 2.444384405113783e+06; /*H2O2 */
    sms[7] *= 5.530081023953346e+06; /*CH3 */
    sms[8] *= 5.182630712527496e+06; /*CH4 */
    sms[9] *= 2.968349425484326e+06; /*CO */
    sms[10] *= 1.889234139098090e+06; /*CO2 */
    sms[11] *= 2.769058254894261e+06; /*CH2O */
    sms[12] *= 3.193192012977835e+06; /*C2H2 */
    sms[13] *= 2.963733033722604e+06; /*C2H4 */
    sms[14] *= 2.765040511976673e+06; /*C2H6 */
    sms[15] *= 4.882097587813943e+06; /*NH3 */
    sms[16] *= 2.770939908885194e+06; /*NO */
    sms[17] *= 3.076506096762281e+06; /*HCN */
    sms[18] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[19]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 19; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[19], tresult[19]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 19; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 19; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[19]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 19; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[19]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*H */
    result += cvor[2]*y[2]*imw[2]; /*O2 */
    result += cvor[3]*y[3]*imw[3]; /*OH */
    result += cvor[4]*y[4]*imw[4]; /*H2O */
    result += cvor[5]*y[5]*imw[5]; /*HO2 */
    result += cvor[6]*y[6]*imw[6]; /*H2O2 */
    result += cvor[7]*y[7]*imw[7]; /*CH3 */
    result += cvor[8]*y[8]*imw[8]; /*CH4 */
    result += cvor[9]*y[9]*imw[9]; /*CO */
    result += cvor[10]*y[10]*imw[10]; /*CO2 */
    result += cvor[11]*y[11]*imw[11]; /*CH2O */
    result += cvor[12]*y[12]*imw[12]; /*C2H2 */
    result += cvor[13]*y[13]*imw[13]; /*C2H4 */
    result += cvor[14]*y[14]*imw[14]; /*C2H6 */
    result += cvor[15]*y[15]*imw[15]; /*NH3 */
    result += cvor[16]*y[16]*imw[16]; /*NO */
    result += cvor[17]*y[17]*imw[17]; /*HCN */
    result += cvor[18]*y[18]*imw[18]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[19]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 19; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
void CKHBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[19], tmp[19]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 19; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 19; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[19]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 19; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[19]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*H */
    result += y[2]*ums[2]*imw[2]; /*O2 */
    result += y[3]*ums[3]*imw[3]; /*OH */
    result += y[4]*ums[4]*imw[4]; /*H2O */
    result += y[5]*ums[5]*imw[5]; /*HO2 */
    result += y[6]*ums[6]*imw[6]; /*H2O2 */
    result += y[7]*ums[7]*imw[7]; /*CH3 */
    result += y[8]*ums[8]*imw[8]; /*CH4 */
    result += y[9]*ums[9]*imw[9]; /*CO */
    result += y[10]*ums[10]*imw[10]; /*CO2 */
    result += y[11]*ums[11]*imw[11]; /*CH2O */
    result += y[12]*ums[12]*imw[12]; /*C2H2 */
    result += y[13]*ums[13]*imw[13]; /*C2H4 */
    result += y[14]*ums[14]*imw[14]; /*C2H6 */
    result += y[15]*ums[15]*imw[15]; /*NH3 */
    result += y[16]*ums[16]*imw[16]; /*NO */
    result += y[17]*ums[17]*imw[17]; /*HCN */
    result += y[18]*ums[18]*imw[18]; /*N2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
void CKSBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict sbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[19]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 19; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31451e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[19]; /* temporary storage */
    double x[19]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O2 */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*HO2 */
    YOW += y[6]*imw[6]; /*H2O2 */
    YOW += y[7]*imw[7]; /*CH3 */
    YOW += y[8]*imw[8]; /*CH4 */
    YOW += y[9]*imw[9]; /*CO */
    YOW += y[10]*imw[10]; /*CO2 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*C2H2 */
    YOW += y[13]*imw[13]; /*C2H4 */
    YOW += y[14]*imw[14]; /*C2H6 */
    YOW += y[15]*imw[15]; /*NH3 */
    YOW += y[16]*imw[16]; /*NO */
    YOW += y[17]*imw[17]; /*HCN */
    YOW += y[18]*imw[18]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(31.998800*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(33.006770*YOW); 
    x[6] = y[6]/(34.014740*YOW); 
    x[7] = y[7]/(15.035060*YOW); 
    x[8] = y[8]/(16.043030*YOW); 
    x[9] = y[9]/(28.010550*YOW); 
    x[10] = y[10]/(44.009950*YOW); 
    x[11] = y[11]/(30.026490*YOW); 
    x[12] = y[12]/(26.038240*YOW); 
    x[13] = y[13]/(28.054180*YOW); 
    x[14] = y[14]/(30.070120*YOW); 
    x[15] = y[15]/(17.030610*YOW); 
    x[16] = y[16]/(30.006100*YOW); 
    x[17] = y[17]/(27.025820*YOW); 
    x[18] = y[18]/(28.013400*YOW); 
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
    /*Scale by R/W */
    *sbms = result * 8.31451e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void CKGBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[19]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 19; ++id) {
        result += x[id]*(gort[id]+log((x[id]+1e-100))+logPratio);
    }

    *gbml = result * RT;
}


/*Returns mixture gibbs free energy in mass units */
void CKGBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[19]; /* temporary storage */
    double x[19]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O2 */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*HO2 */
    YOW += y[6]*imw[6]; /*H2O2 */
    YOW += y[7]*imw[7]; /*CH3 */
    YOW += y[8]*imw[8]; /*CH4 */
    YOW += y[9]*imw[9]; /*CO */
    YOW += y[10]*imw[10]; /*CO2 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*C2H2 */
    YOW += y[13]*imw[13]; /*C2H4 */
    YOW += y[14]*imw[14]; /*C2H6 */
    YOW += y[15]*imw[15]; /*NH3 */
    YOW += y[16]*imw[16]; /*NO */
    YOW += y[17]*imw[17]; /*HCN */
    YOW += y[18]*imw[18]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(31.998800*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(33.006770*YOW); 
    x[6] = y[6]/(34.014740*YOW); 
    x[7] = y[7]/(15.035060*YOW); 
    x[8] = y[8]/(16.043030*YOW); 
    x[9] = y[9]/(28.010550*YOW); 
    x[10] = y[10]/(44.009950*YOW); 
    x[11] = y[11]/(30.026490*YOW); 
    x[12] = y[12]/(26.038240*YOW); 
    x[13] = y[13]/(28.054180*YOW); 
    x[14] = y[14]/(30.070120*YOW); 
    x[15] = y[15]/(17.030610*YOW); 
    x[16] = y[16]/(30.006100*YOW); 
    x[17] = y[17]/(27.025820*YOW); 
    x[18] = y[18]/(28.013400*YOW); 
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
    result += x[12]*(gort[12]+log((x[12]+1e-100))+logPratio);
    result += x[13]*(gort[13]+log((x[13]+1e-100))+logPratio);
    result += x[14]*(gort[14]+log((x[14]+1e-100))+logPratio);
    result += x[15]*(gort[15]+log((x[15]+1e-100))+logPratio);
    result += x[16]*(gort[16]+log((x[16]+1e-100))+logPratio);
    result += x[17]*(gort[17]+log((x[17]+1e-100))+logPratio);
    result += x[18]*(gort[18]+log((x[18]+1e-100))+logPratio);
    /*Scale by RT/W */
    *gbms = result * RT * YOW;
}


/*Returns mean helmholtz free energy in molar units */
void CKABML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[19]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 19; ++id) {
        result += x[id]*(aort[id]+log((x[id]+1e-100))+logPratio);
    }

    *abml = result * RT;
}


/*Returns mixture helmholtz free energy in mass units */
void CKABMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[19]; /* temporary storage */
    double x[19]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O2 */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*HO2 */
    YOW += y[6]*imw[6]; /*H2O2 */
    YOW += y[7]*imw[7]; /*CH3 */
    YOW += y[8]*imw[8]; /*CH4 */
    YOW += y[9]*imw[9]; /*CO */
    YOW += y[10]*imw[10]; /*CO2 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*C2H2 */
    YOW += y[13]*imw[13]; /*C2H4 */
    YOW += y[14]*imw[14]; /*C2H6 */
    YOW += y[15]*imw[15]; /*NH3 */
    YOW += y[16]*imw[16]; /*NO */
    YOW += y[17]*imw[17]; /*HCN */
    YOW += y[18]*imw[18]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(31.998800*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(33.006770*YOW); 
    x[6] = y[6]/(34.014740*YOW); 
    x[7] = y[7]/(15.035060*YOW); 
    x[8] = y[8]/(16.043030*YOW); 
    x[9] = y[9]/(28.010550*YOW); 
    x[10] = y[10]/(44.009950*YOW); 
    x[11] = y[11]/(30.026490*YOW); 
    x[12] = y[12]/(26.038240*YOW); 
    x[13] = y[13]/(28.054180*YOW); 
    x[14] = y[14]/(30.070120*YOW); 
    x[15] = y[15]/(17.030610*YOW); 
    x[16] = y[16]/(30.006100*YOW); 
    x[17] = y[17]/(27.025820*YOW); 
    x[18] = y[18]/(28.013400*YOW); 
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
    result += x[12]*(aort[12]+log((x[12]+1e-100))+logPratio);
    result += x[13]*(aort[13]+log((x[13]+1e-100))+logPratio);
    result += x[14]*(aort[14]+log((x[14]+1e-100))+logPratio);
    result += x[15]*(aort[15]+log((x[15]+1e-100))+logPratio);
    result += x[16]*(aort[16]+log((x[16]+1e-100))+logPratio);
    result += x[17]*(aort[17]+log((x[17]+1e-100))+logPratio);
    result += x[18]*(aort[18]+log((x[18]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
/* void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot) */
/* { */
/*     int id; /\*loop counter *\/ */

/*     /\*convert to SI *\/ */
/*     for (id = 0; id < 19; ++id) { */
/*         C[id] *= 1.0e6; */
/*     } */

/*     /\*convert to chemkin units *\/ */
/*     productionRate(wdot, C, *T); */

/*     /\*convert to chemkin units *\/ */
/*     for (id = 0; id < 19; ++id) { */
/*         C[id] *= 1.0e-6; */
/*         wdot[id] *= 1.0e-6; */
/*     } */
/* } */


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[19]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O2 */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*HO2 */
    YOW += y[6]*imw[6]; /*H2O2 */
    YOW += y[7]*imw[7]; /*CH3 */
    YOW += y[8]*imw[8]; /*CH4 */
    YOW += y[9]*imw[9]; /*CO */
    YOW += y[10]*imw[10]; /*CO2 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*C2H2 */
    YOW += y[13]*imw[13]; /*C2H4 */
    YOW += y[14]*imw[14]; /*C2H6 */
    YOW += y[15]*imw[15]; /*NH3 */
    YOW += y[16]*imw[16]; /*NO */
    YOW += y[17]*imw[17]; /*HCN */
    YOW += y[18]*imw[18]; /*N2 */
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
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 19; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[19]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 19; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 19; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[19]; /*temporary storage */
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

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 19; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
	    double * restrict y, int * restrict iwrk, double * restrict rwrk,
	    double * restrict wdot)
{
    double c[19*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<19; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<19*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[19]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*31.998800; /*O2 */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*15.035060; /*CH3 */
    XW += x[8]*16.043030; /*CH4 */
    XW += x[9]*28.010550; /*CO */
    XW += x[10]*44.009950; /*CO2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*26.038240; /*C2H2 */
    XW += x[13]*28.054180; /*C2H4 */
    XW += x[14]*30.070120; /*C2H6 */
    XW += x[15]*17.030610; /*NH3 */
    XW += x[16]*30.006100; /*NO */
    XW += x[17]*27.025820; /*HCN */
    XW += x[18]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 19; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 19; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 19; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 19; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r)
{
    int id; /*loop counter */
    double c[19]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 19; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[19]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O2 */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*HO2 */
    YOW += y[6]*imw[6]; /*H2O2 */
    YOW += y[7]*imw[7]; /*CH3 */
    YOW += y[8]*imw[8]; /*CH4 */
    YOW += y[9]*imw[9]; /*CO */
    YOW += y[10]*imw[10]; /*CO2 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*C2H2 */
    YOW += y[13]*imw[13]; /*C2H4 */
    YOW += y[14]*imw[14]; /*C2H6 */
    YOW += y[15]*imw[15]; /*NH3 */
    YOW += y[16]*imw[16]; /*NO */
    YOW += y[17]*imw[17]; /*HCN */
    YOW += y[18]*imw[18]; /*N2 */
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
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[19]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 19; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[19]; /*temporary storage */
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

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[19]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*31.998800; /*O2 */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*15.035060; /*CH3 */
    XW += x[8]*16.043030; /*CH4 */
    XW += x[9]*28.010550; /*CO */
    XW += x[10]*44.009950; /*CO2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*26.038240; /*C2H2 */
    XW += x[13]*28.054180; /*C2H4 */
    XW += x[14]*30.070120; /*C2H6 */
    XW += x[15]*17.030610; /*NH3 */
    XW += x[16]*30.006100; /*NO */
    XW += x[17]*27.025820; /*HCN */
    XW += x[18]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 19; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim, int * iwrk, double * restrict rwrk, int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 19 * kd; ++ id) {
         nuki[id] = 0; 
    }
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 19; ++ id) {
         ncf[id] = 0; 
    }

    /*H2 */
    ncf[ 0 * kd + 1 ] = 2; /*H */

    /*H */
    ncf[ 1 * kd + 1 ] = 1; /*H */

    /*O2 */
    ncf[ 2 * kd + 0 ] = 2; /*O */

    /*OH */
    ncf[ 3 * kd + 0 ] = 1; /*O */
    ncf[ 3 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 4 * kd + 1 ] = 2; /*H */
    ncf[ 4 * kd + 0 ] = 1; /*O */

    /*HO2 */
    ncf[ 5 * kd + 1 ] = 1; /*H */
    ncf[ 5 * kd + 0 ] = 2; /*O */

    /*H2O2 */
    ncf[ 6 * kd + 1 ] = 2; /*H */
    ncf[ 6 * kd + 0 ] = 2; /*O */

    /*CH3 */
    ncf[ 7 * kd + 2 ] = 1; /*C */
    ncf[ 7 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 8 * kd + 2 ] = 1; /*C */
    ncf[ 8 * kd + 1 ] = 4; /*H */

    /*CO */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 0 ] = 2; /*O */

    /*CH2O */
    ncf[ 11 * kd + 1 ] = 2; /*H */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 0 ] = 1; /*O */

    /*C2H2 */
    ncf[ 12 * kd + 2 ] = 2; /*C */
    ncf[ 12 * kd + 1 ] = 2; /*H */

    /*C2H4 */
    ncf[ 13 * kd + 2 ] = 2; /*C */
    ncf[ 13 * kd + 1 ] = 4; /*H */

    /*C2H6 */
    ncf[ 14 * kd + 2 ] = 2; /*C */
    ncf[ 14 * kd + 1 ] = 6; /*H */

    /*NH3 */
    ncf[ 15 * kd + 3 ] = 1; /*N */
    ncf[ 15 * kd + 1 ] = 3; /*H */

    /*NO */
    ncf[ 16 * kd + 3 ] = 1; /*N */
    ncf[ 16 * kd + 0 ] = 1; /*O */

    /*HCN */
    ncf[ 17 * kd + 1 ] = 1; /*H */
    ncf[ 17 * kd + 2 ] = 1; /*C */
    ncf[ 17 * kd + 3 ] = 1; /*N */

    /*N2 */
    ncf[ 18 * kd + 3 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e)
{
    for (int i=0; i<0; ++i) {
        a[i] = fwd_A[i];
        b[i] = fwd_beta[i];
        e[i] = fwd_Ea[i];
    }

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[19]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[19]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[19]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[19]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[19]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[0];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[0];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif


/*compute the production rate for each species */
void productionRate(double * restrict wdot, double * restrict sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double qdot, q_f[0], q_r[0];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 19; ++i) {
        wdot[i] = 0.0;
    }

    return;
}

void comp_k_f(double * restrict tc, double invT, double * restrict k_f)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<0; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double * restrict tc, double invT, double * restrict Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[19];
    gibbs(g_RT, tc);


#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<0; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;


    return;
}

void comp_qfqr(double * restrict qf, double * restrict qr, double * restrict sc, double * restrict tc, double invT)
{

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 19; ++i) {
        mixture += sc[i];
    }

    double Corr[0];
    for (int i = 0; i < 0; ++i) {
        Corr[i] = 1.0;
    }

    for (int i=0; i<0; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[0*npt], Kc_s[0*npt], mixture[npt], g_RT[19*npt];
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

    for (int n=0; n<19; n++) {
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

void vcomp_k_f(int npt, double * restrict k_f_s, double * restrict tc, double * restrict invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
    }
}

void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[19];
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
        g_RT[12*npt+i] = g[12];
        g_RT[13*npt+i] = g[13];
        g_RT[14*npt+i] = g[14];
        g_RT[15*npt+i] = g[15];
        g_RT[16*npt+i] = g[16];
        g_RT[17*npt+i] = g[17];
        g_RT[18*npt+i] = g[18];
    }
}

void vcomp_Kc(int npt, double * restrict Kc_s, double * restrict g_RT, double * restrict invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = (101325. / 8.31451) * invT[i];
        double refCinv = 1.0 / refC;

    }
}

void vcomp_wdot(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
		double * restrict k_f_s, double * restrict Kc_s,
		double * restrict tc, double * restrict invT, double * restrict T)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
    }
}

/*compute the reaction Jacobian */
void DWDOT(double * restrict J, double * restrict sc, double * restrict Tp, int * consP)
{
    double c[19];

    for (int k=0; k<19; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    for (int k=0; k<19; k++) {
        J[380+k] *= 1.e-6;
    }

    /* dTdot/d[X] */
    for (int k=0; k<19; k++) {
        J[k*20+19] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
void aJacobian(double * restrict J, double * restrict sc, double T, int consP)
{
    for (int i=0; i<400; i++) {
        J[i] = 0.0;
    }

    double wdot[19];
    for (int k=0; k<19; k++) {
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
    for (int k = 0; k < 19; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[19];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[19];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[19];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    double c_R[19], dcRdT[19], e_RT[19];
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
    for (int k = 0; k < 19; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[380+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 19; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 19; ++m) {
            dehmixdc += eh_RT[m]*J[k*20+m];
        }
        J[k*20+19] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[399] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void dcvpRdT(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

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
        /*species 2: O2 */
        species[2] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 3: OH */
        species[3] =
            -2.40131752e-03
            +9.23587682e-06 * tc[1]
            -1.16434000e-08 * tc[2]
            +5.45645880e-12 * tc[3];
        /*species 4: H2O */
        species[4] =
            -2.03643410e-03
            +1.30408042e-05 * tc[1]
            -1.64639119e-08 * tc[2]
            +7.08791268e-12 * tc[3];
        /*species 5: HO2 */
        species[5] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 6: H2O2 */
        species[6] =
            -5.42822417e-04
            +3.34671402e-05 * tc[1]
            -6.47312439e-08 * tc[2]
            +3.44981745e-11 * tc[3];
        /*species 7: CH3 */
        species[7] =
            +2.01095175e-03
            +1.14604371e-05 * tc[1]
            -2.06135228e-08 * tc[2]
            +1.01754294e-11 * tc[3];
        /*species 8: CH4 */
        species[8] =
            -1.36709788e-02
            +9.83601198e-05 * tc[1]
            -1.45422908e-07 * tc[2]
            +6.66775824e-11 * tc[3];
        /*species 9: CO */
        species[9] =
            -6.10353680e-04
            +2.03362866e-06 * tc[1]
            +2.72101765e-09 * tc[2]
            -3.61769800e-12 * tc[3];
        /*species 10: CO2 */
        species[10] =
            +8.98459677e-03
            -1.42471254e-05 * tc[1]
            +7.37757066e-09 * tc[2]
            -5.74798192e-13 * tc[3];
        /*species 11: CH2O */
        species[11] =
            -9.90833369e-03
            +7.46440016e-05 * tc[1]
            -1.13785578e-07 * tc[2]
            +5.27090608e-11 * tc[3];
        /*species 12: C2H2 */
        species[12] =
            +2.33615629e-02
            -7.10343630e-05 * tc[1]
            +8.40457311e-08 * tc[2]
            -3.40029190e-11 * tc[3];
        /*species 13: C2H4 */
        species[13] =
            -7.57052247e-03
            +1.14198058e-04 * tc[1]
            -2.07476626e-07 * tc[2]
            +1.07953749e-10 * tc[3];
        /*species 14: C2H6 */
        species[14] =
            -5.50154270e-03
            +1.19887658e-04 * tc[1]
            -2.12539886e-07 * tc[2]
            +1.07474308e-10 * tc[3];
        /*species 15: NH3 */
        species[15] =
            -4.66052300e-03
            +4.34370260e-05 * tc[1]
            -6.84266610e-08 * tc[2]
            +3.30552184e-11 * tc[3];
        /*species 16: NO */
        species[16] =
            -4.63897600e-03
            +2.20820440e-05 * tc[1]
            -2.80084062e-08 * tc[2]
            +1.12143080e-11 * tc[3];
        /*species 17: HCN */
        species[17] =
            +1.00511700e-02
            -2.67035260e-05 * tc[1]
            +3.02770470e-08 * tc[2]
            -1.20356112e-11 * tc[3];
        /*species 18: N2 */
        species[18] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
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
        /*species 2: O2 */
        species[2] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 3: OH */
        species[3] =
            +5.48429716e-04
            +2.53010456e-07 * tc[1]
            -2.63838467e-10 * tc[2]
            +4.69649504e-14 * tc[3];
        /*species 4: H2O */
        species[4] =
            +2.17691804e-03
            -3.28145036e-07 * tc[1]
            -2.91125961e-10 * tc[2]
            +6.72803968e-14 * tc[3];
        /*species 5: HO2 */
        species[5] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 6: H2O2 */
        species[6] =
            +4.90831694e-03
            -3.80278450e-06 * tc[1]
            +1.11355796e-09 * tc[2]
            -1.15163322e-13 * tc[3];
        /*species 7: CH3 */
        species[7] =
            +7.23990037e-03
            -5.97428696e-06 * tc[1]
            +1.78705393e-09 * tc[2]
            -1.86861758e-13 * tc[3];
        /*species 8: CH4 */
        species[8] =
            +1.33909467e-02
            -1.14657162e-05 * tc[1]
            +3.66877605e-09 * tc[2]
            -4.07260920e-13 * tc[3];
        /*species 9: CO */
        species[9] =
            +2.06252743e-03
            -1.99765154e-06 * tc[1]
            +6.90159024e-10 * tc[2]
            -8.14590864e-14 * tc[3];
        /*species 10: CO2 */
        species[10] =
            +4.41437026e-03
            -4.42962808e-06 * tc[1]
            +1.57047056e-09 * tc[2]
            -1.88833666e-13 * tc[3];
        /*species 11: CH2O */
        species[11] =
            +9.20000082e-03
            -8.84517626e-06 * tc[1]
            +3.01923636e-09 * tc[2]
            -3.53542256e-13 * tc[3];
        /*species 12: C2H2 */
        species[12] =
            +5.96166664e-03
            -4.74589704e-06 * tc[1]
            +1.40223651e-09 * tc[2]
            -1.44494085e-13 * tc[3];
        /*species 13: C2H4 */
        species[13] =
            +1.46454151e-02
            -1.34215583e-05 * tc[1]
            +4.41668769e-09 * tc[2]
            -5.02824244e-13 * tc[3];
        /*species 14: C2H6 */
        species[14] =
            +2.16852677e-02
            -2.00512134e-05 * tc[1]
            +6.64236003e-09 * tc[2]
            -7.60011560e-13 * tc[3];
        /*species 15: NH3 */
        species[15] =
            +5.66625600e-03
            -3.45573520e-06 * tc[1]
            +7.16014830e-10 * tc[2]
            -5.03151440e-14 * tc[3];
        /*species 16: NO */
        species[16] =
            +1.19110430e-03
            -8.58340960e-07 * tc[1]
            +2.08373007e-10 * tc[2]
            -1.61344396e-14 * tc[3];
        /*species 17: HCN */
        species[17] =
            +3.14642280e-03
            -2.12643700e-06 * tc[1]
            +4.98592710e-10 * tc[2]
            -3.91990280e-14 * tc[3];
        /*species 18: N2 */
        species[18] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
    }
    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * restrict qdot, double * restrict sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double q_f[0], q_r[0];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 0; ++i) {
        qdot[i] = q_f[i] - q_r[i];
    }

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double * restrict q_f, double * restrict q_r, double * restrict sc, double T)
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
void equilibriumConstants(double * restrict kc, double * restrict g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

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
        /*species 2: O2 */
        species[2] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.615080560000000e+03 * invT
            +4.095940888000000e+00
            -3.992015430000000e+00 * tc[0]
            +1.200658760000000e-03 * tc[1]
            -7.696564016666666e-07 * tc[2]
            +3.234277775000000e-10 * tc[3]
            -6.820573500000000e-14 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            -1.770258210000000e+04 * invT
            +8.410619499999998e-01
            -4.276112690000000e+00 * tc[0]
            +2.714112085000000e-04 * tc[1]
            -2.788928350000000e-06 * tc[2]
            +1.798090108333333e-09 * tc[3]
            -4.312271815000000e-13 * tc[4];
        /*species 7: CH3 */
        species[7] =
            +1.644499880000000e+04 * invT
            +2.069026070000000e+00
            -3.673590400000000e+00 * tc[0]
            -1.005475875000000e-03 * tc[1]
            -9.550364266666668e-07 * tc[2]
            +5.725978541666666e-10 * tc[3]
            -1.271928670000000e-13 * tc[4];
        /*species 8: CH4 */
        species[8] =
            -1.024664760000000e+04 * invT
            +9.791179889999999e+00
            -5.149876130000000e+00 * tc[0]
            +6.835489400000000e-03 * tc[1]
            -8.196676650000000e-06 * tc[2]
            +4.039525216666667e-09 * tc[3]
            -8.334697800000000e-13 * tc[4];
        /*species 9: CO */
        species[9] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 10: CO2 */
        species[10] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
        /*species 11: CH2O */
        species[11] =
            -1.430895670000000e+04 * invT
            +4.190910250000000e+00
            -4.793723150000000e+00 * tc[0]
            +4.954166845000000e-03 * tc[1]
            -6.220333466666666e-06 * tc[2]
            +3.160710508333333e-09 * tc[3]
            -6.588632600000000e-13 * tc[4];
        /*species 12: C2H2 */
        species[12] =
            +2.642898070000000e+04 * invT
            -1.313102400600000e+01
            -8.086810940000000e-01 * tc[0]
            -1.168078145000000e-02 * tc[1]
            +5.919530250000000e-06 * tc[2]
            -2.334603641666667e-09 * tc[3]
            +4.250364870000000e-13 * tc[4];
        /*species 13: C2H4 */
        species[13] =
            +5.089775930000000e+03 * invT
            -1.381294799999999e-01
            -3.959201480000000e+00 * tc[0]
            +3.785261235000000e-03 * tc[1]
            -9.516504866666667e-06 * tc[2]
            +5.763239608333333e-09 * tc[3]
            -1.349421865000000e-12 * tc[4];
        /*species 14: C2H6 */
        species[14] =
            -1.152220550000000e+04 * invT
            +1.624601760000000e+00
            -4.291424920000000e+00 * tc[0]
            +2.750771350000000e-03 * tc[1]
            -9.990638133333334e-06 * tc[2]
            +5.903885708333334e-09 * tc[3]
            -1.343428855000000e-12 * tc[4];
        /*species 15: NH3 */
        species[15] =
            -6.741728500000000e+03 * invT
            +4.911400170000000e+00
            -4.286027400000000e+00 * tc[0]
            +2.330261500000000e-03 * tc[1]
            -3.619752166666667e-06 * tc[2]
            +1.900740583333333e-09 * tc[3]
            -4.131902300000000e-13 * tc[4];
        /*species 16: NO */
        species[16] =
            +9.844623000000000e+03 * invT
            +1.937629900000000e+00
            -4.218476300000000e+00 * tc[0]
            +2.319488000000000e-03 * tc[1]
            -1.840170333333333e-06 * tc[2]
            +7.780112833333333e-10 * tc[3]
            -1.401788500000000e-13 * tc[4];
        /*species 17: HCN */
        species[17] =
            +1.471263300000000e+04 * invT
            -6.657453300000000e+00
            -2.258988600000000e+00 * tc[0]
            -5.025585000000000e-03 * tc[1]
            +2.225293833333333e-06 * tc[2]
            -8.410290833333334e-10 * tc[3]
            +1.504451400000000e-13 * tc[4];
        /*species 18: N2 */
        species[18] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
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
        /*species 2: O2 */
        species[2] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.858657000000000e+03 * invT
            -1.383808430000000e+00
            -3.092887670000000e+00 * tc[0]
            -2.742148580000000e-04 * tc[1]
            -2.108420466666667e-08 * tc[2]
            +7.328846300000000e-12 * tc[3]
            -5.870618800000000e-16 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            -1.786178770000000e+04 * invT
            +1.248846229999999e+00
            -4.165002850000000e+00 * tc[0]
            -2.454158470000000e-03 * tc[1]
            +3.168987083333333e-07 * tc[2]
            -3.093216550000000e-11 * tc[3]
            +1.439541525000000e-15 * tc[4];
        /*species 7: CH3 */
        species[7] =
            +1.677558430000000e+04 * invT
            -6.194354070000000e+00
            -2.285717720000000e+00 * tc[0]
            -3.619950185000000e-03 * tc[1]
            +4.978572466666667e-07 * tc[2]
            -4.964038700000000e-11 * tc[3]
            +2.335771970000000e-15 * tc[4];
        /*species 8: CH4 */
        species[8] =
            -9.468344590000001e+03 * invT
            -1.836246650500000e+01
            -7.485149500000000e-02 * tc[0]
            -6.695473350000000e-03 * tc[1]
            +9.554763483333333e-07 * tc[2]
            -1.019104458333333e-10 * tc[3]
            +5.090761500000000e-15 * tc[4];
        /*species 9: CO */
        species[9] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 10: CO2 */
        species[10] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
        /*species 11: CH2O */
        species[11] =
            -1.399583230000000e+04 * invT
            -1.189563292000000e+01
            -1.760690080000000e+00 * tc[0]
            -4.600000410000000e-03 * tc[1]
            +7.370980216666666e-07 * tc[2]
            -8.386767666666666e-11 * tc[3]
            +4.419278200000001e-15 * tc[4];
        /*species 12: C2H2 */
        species[12] =
            +2.593599920000000e+04 * invT
            +5.377850850000001e+00
            -4.147569640000000e+00 * tc[0]
            -2.980833320000000e-03 * tc[1]
            +3.954914200000000e-07 * tc[2]
            -3.895101425000000e-11 * tc[3]
            +1.806176065000000e-15 * tc[4];
        /*species 13: C2H4 */
        species[13] =
            +4.939886140000000e+03 * invT
            -8.269258140000002e+00
            -2.036111160000000e+00 * tc[0]
            -7.322707550000000e-03 * tc[1]
            +1.118463191666667e-06 * tc[2]
            -1.226857691666667e-10 * tc[3]
            +6.285303050000000e-15 * tc[4];
        /*species 14: C2H6 */
        species[14] =
            -1.142639320000000e+04 * invT
            -1.404372920000000e+01
            -1.071881500000000e+00 * tc[0]
            -1.084263385000000e-02 * tc[1]
            +1.670934450000000e-06 * tc[2]
            -1.845100008333333e-10 * tc[3]
            +9.500144500000000e-15 * tc[4];
        /*species 15: NH3 */
        species[15] =
            -6.544695800000000e+03 * invT
            -3.931840700000000e+00
            -2.634452100000000e+00 * tc[0]
            -2.833128000000000e-03 * tc[1]
            +2.879779333333333e-07 * tc[2]
            -1.988930083333333e-11 * tc[3]
            +6.289392999999999e-16 * tc[4];
        /*species 16: NO */
        species[16] =
            +9.920974600000000e+03 * invT
            -3.108697100000001e+00
            -3.260605600000000e+00 * tc[0]
            -5.955521500000000e-04 * tc[1]
            +7.152841333333333e-08 * tc[2]
            -5.788139083333334e-12 * tc[3]
            +2.016804950000000e-16 * tc[4];
        /*species 17: HCN */
        species[17] =
            +1.440729200000000e+04 * invT
            +2.226779100000000e+00
            -3.802239200000000e+00 * tc[0]
            -1.573211400000000e-03 * tc[1]
            +1.772030833333333e-07 * tc[2]
            -1.384979750000000e-11 * tc[3]
            +4.899878500000000e-16 * tc[4];
        /*species 18: N2 */
        species[18] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void helmholtz(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

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
        /*species 2: O2 */
        species[2] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.61508056e+03 * invT
            +3.09594089e+00
            -3.99201543e+00 * tc[0]
            +1.20065876e-03 * tc[1]
            -7.69656402e-07 * tc[2]
            +3.23427778e-10 * tc[3]
            -6.82057350e-14 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.02937267e+04 * invT
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +2.94808040e+02 * invT
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            -1.77025821e+04 * invT
            -1.58938050e-01
            -4.27611269e+00 * tc[0]
            +2.71411208e-04 * tc[1]
            -2.78892835e-06 * tc[2]
            +1.79809011e-09 * tc[3]
            -4.31227182e-13 * tc[4];
        /*species 7: CH3 */
        species[7] =
            +1.64449988e+04 * invT
            +1.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 8: CH4 */
        species[8] =
            -1.02466476e+04 * invT
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 9: CO */
        species[9] =
            -1.43440860e+04 * invT
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 10: CO2 */
        species[10] =
            -4.83719697e+04 * invT
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 11: CH2O */
        species[11] =
            -1.43089567e+04 * invT
            +3.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 12: C2H2 */
        species[12] =
            +2.64289807e+04 * invT
            -1.41310240e+01
            -8.08681094e-01 * tc[0]
            -1.16807815e-02 * tc[1]
            +5.91953025e-06 * tc[2]
            -2.33460364e-09 * tc[3]
            +4.25036487e-13 * tc[4];
        /*species 13: C2H4 */
        species[13] =
            +5.08977593e+03 * invT
            -1.13812948e+00
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 14: C2H6 */
        species[14] =
            -1.15222055e+04 * invT
            +6.24601760e-01
            -4.29142492e+00 * tc[0]
            +2.75077135e-03 * tc[1]
            -9.99063813e-06 * tc[2]
            +5.90388571e-09 * tc[3]
            -1.34342886e-12 * tc[4];
        /*species 15: NH3 */
        species[15] =
            -6.74172850e+03 * invT
            +3.91140017e+00
            -4.28602740e+00 * tc[0]
            +2.33026150e-03 * tc[1]
            -3.61975217e-06 * tc[2]
            +1.90074058e-09 * tc[3]
            -4.13190230e-13 * tc[4];
        /*species 16: NO */
        species[16] =
            +9.84462300e+03 * invT
            +9.37629900e-01
            -4.21847630e+00 * tc[0]
            +2.31948800e-03 * tc[1]
            -1.84017033e-06 * tc[2]
            +7.78011283e-10 * tc[3]
            -1.40178850e-13 * tc[4];
        /*species 17: HCN */
        species[17] =
            +1.47126330e+04 * invT
            -7.65745330e+00
            -2.25898860e+00 * tc[0]
            -5.02558500e-03 * tc[1]
            +2.22529383e-06 * tc[2]
            -8.41029083e-10 * tc[3]
            +1.50445140e-13 * tc[4];
        /*species 18: N2 */
        species[18] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
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
        /*species 2: O2 */
        species[2] =
            -1.08845772e+03 * invT
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.85865700e+03 * invT
            -2.38380843e+00
            -3.09288767e+00 * tc[0]
            -2.74214858e-04 * tc[1]
            -2.10842047e-08 * tc[2]
            +7.32884630e-12 * tc[3]
            -5.87061880e-16 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.00042971e+04 * invT
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +1.11856713e+02 * invT
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            -1.78617877e+04 * invT
            +2.48846230e-01
            -4.16500285e+00 * tc[0]
            -2.45415847e-03 * tc[1]
            +3.16898708e-07 * tc[2]
            -3.09321655e-11 * tc[3]
            +1.43954153e-15 * tc[4];
        /*species 7: CH3 */
        species[7] =
            +1.67755843e+04 * invT
            -7.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 8: CH4 */
        species[8] =
            -9.46834459e+03 * invT
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 9: CO */
        species[9] =
            -1.41518724e+04 * invT
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 10: CO2 */
        species[10] =
            -4.87591660e+04 * invT
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 11: CH2O */
        species[11] =
            -1.39958323e+04 * invT
            -1.28956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 12: C2H2 */
        species[12] =
            +2.59359992e+04 * invT
            +4.37785085e+00
            -4.14756964e+00 * tc[0]
            -2.98083332e-03 * tc[1]
            +3.95491420e-07 * tc[2]
            -3.89510143e-11 * tc[3]
            +1.80617607e-15 * tc[4];
        /*species 13: C2H4 */
        species[13] =
            +4.93988614e+03 * invT
            -9.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 14: C2H6 */
        species[14] =
            -1.14263932e+04 * invT
            -1.50437292e+01
            -1.07188150e+00 * tc[0]
            -1.08426339e-02 * tc[1]
            +1.67093445e-06 * tc[2]
            -1.84510001e-10 * tc[3]
            +9.50014450e-15 * tc[4];
        /*species 15: NH3 */
        species[15] =
            -6.54469580e+03 * invT
            -4.93184070e+00
            -2.63445210e+00 * tc[0]
            -2.83312800e-03 * tc[1]
            +2.87977933e-07 * tc[2]
            -1.98893008e-11 * tc[3]
            +6.28939300e-16 * tc[4];
        /*species 16: NO */
        species[16] =
            +9.92097460e+03 * invT
            -4.10869710e+00
            -3.26060560e+00 * tc[0]
            -5.95552150e-04 * tc[1]
            +7.15284133e-08 * tc[2]
            -5.78813908e-12 * tc[3]
            +2.01680495e-16 * tc[4];
        /*species 17: HCN */
        species[17] =
            +1.44072920e+04 * invT
            +1.22677910e+00
            -3.80223920e+00 * tc[0]
            -1.57321140e-03 * tc[1]
            +1.77203083e-07 * tc[2]
            -1.38497975e-11 * tc[3]
            +4.89987850e-16 * tc[4];
        /*species 18: N2 */
        species[18] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cv_R(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

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
        /*species 2: O2 */
        species[2] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 3: OH */
        species[3] =
            +2.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 4: H2O */
        species[4] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            +3.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 7: CH3 */
        species[7] =
            +2.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 8: CH4 */
        species[8] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 9: CO */
        species[9] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 10: CO2 */
        species[10] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 11: CH2O */
        species[11] =
            +3.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 12: C2H2 */
        species[12] =
            -1.91318906e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 13: C2H4 */
        species[13] =
            +2.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 14: C2H6 */
        species[14] =
            +3.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 15: NH3 */
        species[15] =
            +3.28602740e+00
            -4.66052300e-03 * tc[1]
            +2.17185130e-05 * tc[2]
            -2.28088870e-08 * tc[3]
            +8.26380460e-12 * tc[4];
        /*species 16: NO */
        species[16] =
            +3.21847630e+00
            -4.63897600e-03 * tc[1]
            +1.10410220e-05 * tc[2]
            -9.33613540e-09 * tc[3]
            +2.80357700e-12 * tc[4];
        /*species 17: HCN */
        species[17] =
            +1.25898860e+00
            +1.00511700e-02 * tc[1]
            -1.33517630e-05 * tc[2]
            +1.00923490e-08 * tc[3]
            -3.00890280e-12 * tc[4];
        /*species 18: N2 */
        species[18] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
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
        /*species 2: O2 */
        species[2] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 3: OH */
        species[3] =
            +2.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 4: H2O */
        species[4] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            +3.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 7: CH3 */
        species[7] =
            +1.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 8: CH4 */
        species[8] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 9: CO */
        species[9] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 10: CO2 */
        species[10] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 11: CH2O */
        species[11] =
            +7.60690080e-01
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 12: C2H2 */
        species[12] =
            +3.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 13: C2H4 */
        species[13] =
            +1.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 14: C2H6 */
        species[14] =
            +7.18815000e-02
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 15: NH3 */
        species[15] =
            +1.63445210e+00
            +5.66625600e-03 * tc[1]
            -1.72786760e-06 * tc[2]
            +2.38671610e-10 * tc[3]
            -1.25787860e-14 * tc[4];
        /*species 16: NO */
        species[16] =
            +2.26060560e+00
            +1.19110430e-03 * tc[1]
            -4.29170480e-07 * tc[2]
            +6.94576690e-11 * tc[3]
            -4.03360990e-15 * tc[4];
        /*species 17: HCN */
        species[17] =
            +2.80223920e+00
            +3.14642280e-03 * tc[1]
            -1.06321850e-06 * tc[2]
            +1.66197570e-10 * tc[3]
            -9.79975700e-15 * tc[4];
        /*species 18: N2 */
        species[18] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cp_R(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

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
        /*species 2: O2 */
        species[2] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 4: H2O */
        species[4] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            +4.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 7: CH3 */
        species[7] =
            +3.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 8: CH4 */
        species[8] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 9: CO */
        species[9] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 10: CO2 */
        species[10] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 11: CH2O */
        species[11] =
            +4.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 12: C2H2 */
        species[12] =
            +8.08681094e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 13: C2H4 */
        species[13] =
            +3.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 14: C2H6 */
        species[14] =
            +4.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 15: NH3 */
        species[15] =
            +4.28602740e+00
            -4.66052300e-03 * tc[1]
            +2.17185130e-05 * tc[2]
            -2.28088870e-08 * tc[3]
            +8.26380460e-12 * tc[4];
        /*species 16: NO */
        species[16] =
            +4.21847630e+00
            -4.63897600e-03 * tc[1]
            +1.10410220e-05 * tc[2]
            -9.33613540e-09 * tc[3]
            +2.80357700e-12 * tc[4];
        /*species 17: HCN */
        species[17] =
            +2.25898860e+00
            +1.00511700e-02 * tc[1]
            -1.33517630e-05 * tc[2]
            +1.00923490e-08 * tc[3]
            -3.00890280e-12 * tc[4];
        /*species 18: N2 */
        species[18] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
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
        /*species 2: O2 */
        species[2] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 4: H2O */
        species[4] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            +4.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 7: CH3 */
        species[7] =
            +2.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 8: CH4 */
        species[8] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 9: CO */
        species[9] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 10: CO2 */
        species[10] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 11: CH2O */
        species[11] =
            +1.76069008e+00
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 12: C2H2 */
        species[12] =
            +4.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 13: C2H4 */
        species[13] =
            +2.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 14: C2H6 */
        species[14] =
            +1.07188150e+00
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 15: NH3 */
        species[15] =
            +2.63445210e+00
            +5.66625600e-03 * tc[1]
            -1.72786760e-06 * tc[2]
            +2.38671610e-10 * tc[3]
            -1.25787860e-14 * tc[4];
        /*species 16: NO */
        species[16] =
            +3.26060560e+00
            +1.19110430e-03 * tc[1]
            -4.29170480e-07 * tc[2]
            +6.94576690e-11 * tc[3]
            -4.03360990e-15 * tc[4];
        /*species 17: HCN */
        species[17] =
            +3.80223920e+00
            +3.14642280e-03 * tc[1]
            -1.06321850e-06 * tc[2]
            +1.66197570e-10 * tc[3]
            -9.79975700e-15 * tc[4];
        /*species 18: N2 */
        species[18] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesInternalEnergy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

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
        /*species 2: O2 */
        species[2] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 3: OH */
        species[3] =
            +2.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 5: HO2 */
        species[5] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 6: H2O2 */
        species[6] =
            +3.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 7: CH3 */
        species[7] =
            +2.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 8: CH4 */
        species[8] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 9: CO */
        species[9] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 10: CO2 */
        species[10] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 11: CH2O */
        species[11] =
            +3.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 12: C2H2 */
        species[12] =
            -1.91318906e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 * invT;
        /*species 13: C2H4 */
        species[13] =
            +2.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 14: C2H6 */
        species[14] =
            +3.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 15: NH3 */
        species[15] =
            +3.28602740e+00
            -2.33026150e-03 * tc[1]
            +7.23950433e-06 * tc[2]
            -5.70222175e-09 * tc[3]
            +1.65276092e-12 * tc[4]
            -6.74172850e+03 * invT;
        /*species 16: NO */
        species[16] =
            +3.21847630e+00
            -2.31948800e-03 * tc[1]
            +3.68034067e-06 * tc[2]
            -2.33403385e-09 * tc[3]
            +5.60715400e-13 * tc[4]
            +9.84462300e+03 * invT;
        /*species 17: HCN */
        species[17] =
            +1.25898860e+00
            +5.02558500e-03 * tc[1]
            -4.45058767e-06 * tc[2]
            +2.52308725e-09 * tc[3]
            -6.01780560e-13 * tc[4]
            +1.47126330e+04 * invT;
        /*species 18: N2 */
        species[18] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
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
        /*species 2: O2 */
        species[2] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 3: OH */
        species[3] =
            +2.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 5: HO2 */
        species[5] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 6: H2O2 */
        species[6] =
            +3.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 7: CH3 */
        species[7] =
            +1.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 8: CH4 */
        species[8] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 9: CO */
        species[9] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 10: CO2 */
        species[10] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 11: CH2O */
        species[11] =
            +7.60690080e-01
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 12: C2H2 */
        species[12] =
            +3.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 * invT;
        /*species 13: C2H4 */
        species[13] =
            +1.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 14: C2H6 */
        species[14] =
            +7.18815000e-02
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 15: NH3 */
        species[15] =
            +1.63445210e+00
            +2.83312800e-03 * tc[1]
            -5.75955867e-07 * tc[2]
            +5.96679025e-11 * tc[3]
            -2.51575720e-15 * tc[4]
            -6.54469580e+03 * invT;
        /*species 16: NO */
        species[16] =
            +2.26060560e+00
            +5.95552150e-04 * tc[1]
            -1.43056827e-07 * tc[2]
            +1.73644173e-11 * tc[3]
            -8.06721980e-16 * tc[4]
            +9.92097460e+03 * invT;
        /*species 17: HCN */
        species[17] =
            +2.80223920e+00
            +1.57321140e-03 * tc[1]
            -3.54406167e-07 * tc[2]
            +4.15493925e-11 * tc[3]
            -1.95995140e-15 * tc[4]
            +1.44072920e+04 * invT;
        /*species 18: N2 */
        species[18] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

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
        /*species 2: O2 */
        species[2] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 3: OH */
        species[3] =
            +3.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 5: HO2 */
        species[5] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 6: H2O2 */
        species[6] =
            +4.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 7: CH3 */
        species[7] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 8: CH4 */
        species[8] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 9: CO */
        species[9] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 10: CO2 */
        species[10] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 11: CH2O */
        species[11] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 12: C2H2 */
        species[12] =
            +8.08681094e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 * invT;
        /*species 13: C2H4 */
        species[13] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 14: C2H6 */
        species[14] =
            +4.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 15: NH3 */
        species[15] =
            +4.28602740e+00
            -2.33026150e-03 * tc[1]
            +7.23950433e-06 * tc[2]
            -5.70222175e-09 * tc[3]
            +1.65276092e-12 * tc[4]
            -6.74172850e+03 * invT;
        /*species 16: NO */
        species[16] =
            +4.21847630e+00
            -2.31948800e-03 * tc[1]
            +3.68034067e-06 * tc[2]
            -2.33403385e-09 * tc[3]
            +5.60715400e-13 * tc[4]
            +9.84462300e+03 * invT;
        /*species 17: HCN */
        species[17] =
            +2.25898860e+00
            +5.02558500e-03 * tc[1]
            -4.45058767e-06 * tc[2]
            +2.52308725e-09 * tc[3]
            -6.01780560e-13 * tc[4]
            +1.47126330e+04 * invT;
        /*species 18: N2 */
        species[18] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
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
        /*species 2: O2 */
        species[2] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 3: OH */
        species[3] =
            +3.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 5: HO2 */
        species[5] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 6: H2O2 */
        species[6] =
            +4.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 7: CH3 */
        species[7] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 8: CH4 */
        species[8] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 9: CO */
        species[9] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 10: CO2 */
        species[10] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 11: CH2O */
        species[11] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 12: C2H2 */
        species[12] =
            +4.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 * invT;
        /*species 13: C2H4 */
        species[13] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 14: C2H6 */
        species[14] =
            +1.07188150e+00
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 15: NH3 */
        species[15] =
            +2.63445210e+00
            +2.83312800e-03 * tc[1]
            -5.75955867e-07 * tc[2]
            +5.96679025e-11 * tc[3]
            -2.51575720e-15 * tc[4]
            -6.54469580e+03 * invT;
        /*species 16: NO */
        species[16] =
            +3.26060560e+00
            +5.95552150e-04 * tc[1]
            -1.43056827e-07 * tc[2]
            +1.73644173e-11 * tc[3]
            -8.06721980e-16 * tc[4]
            +9.92097460e+03 * invT;
        /*species 17: HCN */
        species[17] =
            +3.80223920e+00
            +1.57321140e-03 * tc[1]
            -3.54406167e-07 * tc[2]
            +4.15493925e-11 * tc[3]
            -1.95995140e-15 * tc[4]
            +1.44072920e+04 * invT;
        /*species 18: N2 */
        species[18] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEntropy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

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
        /*species 2: O2 */
        species[2] =
            +3.78245636e+00 * tc[0]
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 3: OH */
        species[3] =
            +3.99201543e+00 * tc[0]
            -2.40131752e-03 * tc[1]
            +2.30896920e-06 * tc[2]
            -1.29371111e-09 * tc[3]
            +3.41028675e-13 * tc[4]
            -1.03925458e-01 ;
        /*species 4: H2O */
        species[4] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
        /*species 5: HO2 */
        species[5] =
            +4.30179801e+00 * tc[0]
            -4.74912051e-03 * tc[1]
            +1.05791445e-05 * tc[2]
            -8.09212980e-09 * tc[3]
            +2.32306281e-12 * tc[4]
            +3.71666245e+00 ;
        /*species 6: H2O2 */
        species[6] =
            +4.27611269e+00 * tc[0]
            -5.42822417e-04 * tc[1]
            +8.36678505e-06 * tc[2]
            -7.19236043e-09 * tc[3]
            +2.15613591e-12 * tc[4]
            +3.43505074e+00 ;
        /*species 7: CH3 */
        species[7] =
            +3.67359040e+00 * tc[0]
            +2.01095175e-03 * tc[1]
            +2.86510928e-06 * tc[2]
            -2.29039142e-09 * tc[3]
            +6.35964335e-13 * tc[4]
            +1.60456433e+00 ;
        /*species 8: CH4 */
        species[8] =
            +5.14987613e+00 * tc[0]
            -1.36709788e-02 * tc[1]
            +2.45900299e-05 * tc[2]
            -1.61581009e-08 * tc[3]
            +4.16734890e-12 * tc[4]
            -4.64130376e+00 ;
        /*species 9: CO */
        species[9] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 10: CO2 */
        species[10] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
        /*species 11: CH2O */
        species[11] =
            +4.79372315e+00 * tc[0]
            -9.90833369e-03 * tc[1]
            +1.86610004e-05 * tc[2]
            -1.26428420e-08 * tc[3]
            +3.29431630e-12 * tc[4]
            +6.02812900e-01 ;
        /*species 12: C2H2 */
        species[12] =
            +8.08681094e-01 * tc[0]
            +2.33615629e-02 * tc[1]
            -1.77585907e-05 * tc[2]
            +9.33841457e-09 * tc[3]
            -2.12518243e-12 * tc[4]
            +1.39397051e+01 ;
        /*species 13: C2H4 */
        species[13] =
            +3.95920148e+00 * tc[0]
            -7.57052247e-03 * tc[1]
            +2.85495146e-05 * tc[2]
            -2.30529584e-08 * tc[3]
            +6.74710933e-12 * tc[4]
            +4.09733096e+00 ;
        /*species 14: C2H6 */
        species[14] =
            +4.29142492e+00 * tc[0]
            -5.50154270e-03 * tc[1]
            +2.99719144e-05 * tc[2]
            -2.36155428e-08 * tc[3]
            +6.71714427e-12 * tc[4]
            +2.66682316e+00 ;
        /*species 15: NH3 */
        species[15] =
            +4.28602740e+00 * tc[0]
            -4.66052300e-03 * tc[1]
            +1.08592565e-05 * tc[2]
            -7.60296233e-09 * tc[3]
            +2.06595115e-12 * tc[4]
            -6.25372770e-01 ;
        /*species 16: NO */
        species[16] =
            +4.21847630e+00 * tc[0]
            -4.63897600e-03 * tc[1]
            +5.52051100e-06 * tc[2]
            -3.11204513e-09 * tc[3]
            +7.00894250e-13 * tc[4]
            +2.28084640e+00 ;
        /*species 17: HCN */
        species[17] =
            +2.25898860e+00 * tc[0]
            +1.00511700e-02 * tc[1]
            -6.67588150e-06 * tc[2]
            +3.36411633e-09 * tc[3]
            -7.52225700e-13 * tc[4]
            +8.91644190e+00 ;
        /*species 18: N2 */
        species[18] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
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
        /*species 2: O2 */
        species[2] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 3: OH */
        species[3] =
            +3.09288767e+00 * tc[0]
            +5.48429716e-04 * tc[1]
            +6.32526140e-08 * tc[2]
            -2.93153852e-11 * tc[3]
            +2.93530940e-15 * tc[4]
            +4.47669610e+00 ;
        /*species 4: H2O */
        species[4] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
        /*species 5: HO2 */
        species[5] =
            +4.01721090e+00 * tc[0]
            +2.23982013e-03 * tc[1]
            -3.16829075e-07 * tc[2]
            +3.80821233e-11 * tc[3]
            -2.69771337e-15 * tc[4]
            +3.78510215e+00 ;
        /*species 6: H2O2 */
        species[6] =
            +4.16500285e+00 * tc[0]
            +4.90831694e-03 * tc[1]
            -9.50696125e-07 * tc[2]
            +1.23728662e-10 * tc[3]
            -7.19770763e-15 * tc[4]
            +2.91615662e+00 ;
        /*species 7: CH3 */
        species[7] =
            +2.28571772e+00 * tc[0]
            +7.23990037e-03 * tc[1]
            -1.49357174e-06 * tc[2]
            +1.98561548e-10 * tc[3]
            -1.16788599e-14 * tc[4]
            +8.48007179e+00 ;
        /*species 8: CH4 */
        species[8] =
            +7.48514950e-02 * tc[0]
            +1.33909467e-02 * tc[1]
            -2.86642905e-06 * tc[2]
            +4.07641783e-10 * tc[3]
            -2.54538075e-14 * tc[4]
            +1.84373180e+01 ;
        /*species 9: CO */
        species[9] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 10: CO2 */
        species[10] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
        /*species 11: CH2O */
        species[11] =
            +1.76069008e+00 * tc[0]
            +9.20000082e-03 * tc[1]
            -2.21129406e-06 * tc[2]
            +3.35470707e-10 * tc[3]
            -2.20963910e-14 * tc[4]
            +1.36563230e+01 ;
        /*species 12: C2H2 */
        species[12] =
            +4.14756964e+00 * tc[0]
            +5.96166664e-03 * tc[1]
            -1.18647426e-06 * tc[2]
            +1.55804057e-10 * tc[3]
            -9.03088033e-15 * tc[4]
            -1.23028121e+00 ;
        /*species 13: C2H4 */
        species[13] =
            +2.03611116e+00 * tc[0]
            +1.46454151e-02 * tc[1]
            -3.35538958e-06 * tc[2]
            +4.90743077e-10 * tc[3]
            -3.14265152e-14 * tc[4]
            +1.03053693e+01 ;
        /*species 14: C2H6 */
        species[14] =
            +1.07188150e+00 * tc[0]
            +2.16852677e-02 * tc[1]
            -5.01280335e-06 * tc[2]
            +7.38040003e-10 * tc[3]
            -4.75007225e-14 * tc[4]
            +1.51156107e+01 ;
        /*species 15: NH3 */
        species[15] =
            +2.63445210e+00 * tc[0]
            +5.66625600e-03 * tc[1]
            -8.63933800e-07 * tc[2]
            +7.95572033e-11 * tc[3]
            -3.14469650e-15 * tc[4]
            +6.56629280e+00 ;
        /*species 16: NO */
        species[16] =
            +3.26060560e+00 * tc[0]
            +1.19110430e-03 * tc[1]
            -2.14585240e-07 * tc[2]
            +2.31525563e-11 * tc[3]
            -1.00840247e-15 * tc[4]
            +6.36930270e+00 ;
        /*species 17: HCN */
        species[17] =
            +3.80223920e+00 * tc[0]
            +3.14642280e-03 * tc[1]
            -5.31609250e-07 * tc[2]
            +5.53991900e-11 * tc[3]
            -2.44993925e-15 * tc[4]
            +1.57546010e+00 ;
        /*species 18: N2 */
        species[18] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double * restrict wt)
{
    wt[0] = 2.015940; /*H2 */
    wt[1] = 1.007970; /*H */
    wt[2] = 31.998800; /*O2 */
    wt[3] = 17.007370; /*OH */
    wt[4] = 18.015340; /*H2O */
    wt[5] = 33.006770; /*HO2 */
    wt[6] = 34.014740; /*H2O2 */
    wt[7] = 15.035060; /*CH3 */
    wt[8] = 16.043030; /*CH4 */
    wt[9] = 28.010550; /*CO */
    wt[10] = 44.009950; /*CO2 */
    wt[11] = 30.026490; /*CH2O */
    wt[12] = 26.038240; /*C2H2 */
    wt[13] = 28.054180; /*C2H4 */
    wt[14] = 30.070120; /*C2H6 */
    wt[15] = 17.030610; /*NH3 */
    wt[16] = 30.006100; /*NO */
    wt[17] = 27.025820; /*HCN */
    wt[18] = 28.013400; /*N2 */

    return;
}


/*save atomic weights into array */
void atomicWeight(double * restrict awt)
{
    awt[0] = 15.999400; /*O */
    awt[1] = 1.007970; /*H */
    awt[2] = 12.011150; /*C */
    awt[3] = 14.006700; /*N */

    return;
}
/* get temperature given internal energy in mass units and mass fracs */
void GET_T_GIVEN_EY(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int * ierr)
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
    CKUBMS(&tmin, y, iwrk, rwrk, &emin);
    CKUBMS(&tmax, y, iwrk, rwrk, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, iwrk, rwrk, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, iwrk, rwrk, &cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,iwrk,rwrk,&e1);
        CKCVBS(&t1,y,iwrk,rwrk,&cv);
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
void GET_T_GIVEN_HY(double * restrict h, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int * ierr)
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
    CKHBMS(&tmin, y, iwrk, rwrk, &hmin);
    CKHBMS(&tmax, y, iwrk, rwrk, &hmax);
    if (hin < hmin) {
        /*Linear Extrapolation below tmin */
        CKCPBS(&tmin, y, iwrk, rwrk, &cp);
        *t = tmin - (hmin-hin)/cp;
        *ierr = 1;
        return;
    }
    if (hin > hmax) {
        /*Linear Extrapolation above tmax */
        CKCPBS(&tmax, y, iwrk, rwrk, &cp);
        *t = tmax - (hmax-hin)/cp;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKHBMS(&t1,y,iwrk,rwrk,&h1);
        CKCPBS(&t1,y,iwrk,rwrk,&cp);
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
void GET_CRITPARAMS(double * restrict Tci, double * restrict ai, double * restrict bi, double * restrict acentric_i)
{

    double   EPS[19];
    double   SIG[19];
    double    wt[19];
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

    /*species 1: H */
    Tci[1] = 1.316 * EPS[1] ; 
    ai[1] = (5.55 * pow(avogadro,2.0) * EPS[1]*boltzmann * pow(1e-8*SIG[1],3.0) ) / (pow(wt[1],2.0)); 
    bi[1] = 0.855 * avogadro * pow(1e-8*SIG[1],3.0) / (wt[1]); 
    acentric_i[1] = 0.0 ;

    /*species 2: O2 */
    /*Imported from NIST */
    Tci[2] = 154.581000 ; 
    ai[2] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[2],2.0) / (pow(31.998800,2.0) * 50.430466); 
    bi[2] = 0.08664 * Rcst * Tci[2] / (31.998800 * 50.430466); 
    acentric_i[2] = 0.022200 ;

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

    /*species 5: HO2 */
    Tci[5] = 1.316 * EPS[5] ; 
    ai[5] = (5.55 * pow(avogadro,2.0) * EPS[5]*boltzmann * pow(1e-8*SIG[5],3.0) ) / (pow(wt[5],2.0)); 
    bi[5] = 0.855 * avogadro * pow(1e-8*SIG[5],3.0) / (wt[5]); 
    acentric_i[5] = 0.0 ;

    /*species 6: H2O2 */
    Tci[6] = 1.316 * EPS[6] ; 
    ai[6] = (5.55 * pow(avogadro,2.0) * EPS[6]*boltzmann * pow(1e-8*SIG[6],3.0) ) / (pow(wt[6],2.0)); 
    bi[6] = 0.855 * avogadro * pow(1e-8*SIG[6],3.0) / (wt[6]); 
    acentric_i[6] = 0.0 ;

    /*species 7: CH3 */
    Tci[7] = 1.316 * EPS[7] ; 
    ai[7] = (5.55 * pow(avogadro,2.0) * EPS[7]*boltzmann * pow(1e-8*SIG[7],3.0) ) / (pow(wt[7],2.0)); 
    bi[7] = 0.855 * avogadro * pow(1e-8*SIG[7],3.0) / (wt[7]); 
    acentric_i[7] = 0.0 ;

    /*species 8: CH4 */
    /*Imported from NIST */
    Tci[8] = 190.560000 ; 
    ai[8] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[8],2.0) / (pow(16.043030,2.0) * 45.990000); 
    bi[8] = 0.08664 * Rcst * Tci[8] / (16.043030 * 45.990000); 
    acentric_i[8] = 0.011000 ;

    /*species 9: CO */
    /*Imported from NIST */
    Tci[9] = 132.850000 ; 
    ai[9] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[9],2.0) / (pow(28.010000,2.0) * 34.940000); 
    bi[9] = 0.08664 * Rcst * Tci[9] / (28.010000 * 34.940000); 
    acentric_i[9] = 0.045000 ;

    /*species 10: CO2 */
    /*Imported from NIST */
    Tci[10] = 304.120000 ; 
    ai[10] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[10],2.0) / (pow(44.009950,2.0) * 73.740000); 
    bi[10] = 0.08664 * Rcst * Tci[10] / (44.009950 * 73.740000); 
    acentric_i[10] = 0.225000 ;

    /*species 11: CH2O */
    Tci[11] = 1.316 * EPS[11] ; 
    ai[11] = (5.55 * pow(avogadro,2.0) * EPS[11]*boltzmann * pow(1e-8*SIG[11],3.0) ) / (pow(wt[11],2.0)); 
    bi[11] = 0.855 * avogadro * pow(1e-8*SIG[11],3.0) / (wt[11]); 
    acentric_i[11] = 0.0 ;

    /*species 12: C2H2 */
    /*Imported from NIST */
    Tci[12] = 308.300000 ; 
    ai[12] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[12],2.0) / (pow(26.038000,2.0) * 61.140000); 
    bi[12] = 0.08664 * Rcst * Tci[12] / (26.038000 * 61.140000); 
    acentric_i[12] = 0.189000 ;

    /*species 13: C2H4 */
    /*Imported from NIST */
    Tci[13] = 282.340000 ; 
    ai[13] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[13],2.0) / (pow(28.054000,2.0) * 50.410000); 
    bi[13] = 0.08664 * Rcst * Tci[13] / (28.054000 * 50.410000); 
    acentric_i[13] = 0.087000 ;

    /*species 14: C2H6 */
    /*Imported from NIST */
    Tci[14] = 305.320000 ; 
    ai[14] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[14],2.0) / (pow(30.070120,2.0) * 48.720000); 
    bi[14] = 0.08664 * Rcst * Tci[14] / (30.070120 * 48.720000); 
    acentric_i[14] = 0.099000 ;

    /*species 15: NH3 */
    Tci[15] = 1.316 * EPS[15] ; 
    ai[15] = (5.55 * pow(avogadro,2.0) * EPS[15]*boltzmann * pow(1e-8*SIG[15],3.0) ) / (pow(wt[15],2.0)); 
    bi[15] = 0.855 * avogadro * pow(1e-8*SIG[15],3.0) / (wt[15]); 
    acentric_i[15] = 0.0 ;

    /*species 16: NO */
    /*Imported from NIST */
    Tci[16] = 180.000000 ; 
    ai[16] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[16],2.0) / (pow(30.006000,2.0) * 64.800000); 
    bi[16] = 0.08664 * Rcst * Tci[16] / (30.006000 * 64.800000); 
    acentric_i[16] = 0.582000 ;

    /*species 17: HCN */
    Tci[17] = 1.316 * EPS[17] ; 
    ai[17] = (5.55 * pow(avogadro,2.0) * EPS[17]*boltzmann * pow(1e-8*SIG[17],3.0) ) / (pow(wt[17],2.0)); 
    bi[17] = 0.855 * avogadro * pow(1e-8*SIG[17],3.0) / (wt[17]); 
    acentric_i[17] = 0.0 ;

    /*species 18: N2 */
    /*Imported from NIST */
    Tci[18] = 126.192000 ; 
    ai[18] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[18],2.0) / (pow(28.013400,2.0) * 33.958000); 
    bi[18] = 0.08664 * Rcst * Tci[18] / (28.013400 * 33.958000); 
    acentric_i[18] = 0.037200 ;

    return;
}

/* End of file  */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENIMC EGTRANSETLENIMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENIMC egtransetlenimc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENIMC egtransetlenimc_
#endif
void egtransetLENIMC(int* LENIMC) {
  *LENIMC =           70;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void egtransetLENRMC(int* LENRMC) {
  *LENRMC =         6086;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNO EGTRANSETNO
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNO egtransetno
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNO egtransetno_
#endif
void egtransetNO(int* NO) {
  *NO =            4;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKK EGTRANSETKK
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKK egtransetkk
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKK egtransetkk_
#endif
void egtransetKK(int* KK) {
  *KK =           17;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLITE EGTRANSETNLITE
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLITE egtransetnlite
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLITE egtransetnlite_
#endif
void egtransetNLITE(int* NLITE) {
  *NLITE =            2;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPATM EGTRANSETPATM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPATM egtransetpatm
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPATM egtransetpatm_
#endif
void egtransetPATM(double* PATM) {
  *PATM =   0.1013250000000000E+07;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetWT EGTRANSETWT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetWT egtransetwt
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetWT egtransetwt_
#endif
void egtransetWT(double* WT) {
  WT[           0] =   0.2015939950942993E+01;
  WT[           1] =   0.1007969975471497E+01;
  WT[           2] =   0.3199880027770996E+02;
  WT[           3] =   0.1700737011432648E+02;
  WT[           4] =   0.1801534008979797E+02;
  WT[           5] =   0.3300677025318146E+02;
  WT[           6] =   0.3401474022865295E+02;
  WT[           7] =   0.1503506028652191E+02;
  WT[           8] =   0.1604303026199341E+02;
  WT[           9] =   0.2801055049896240E+02;
  WT[          10] =   0.4400995063781738E+02;
  WT[          11] =   0.3002649044990540E+02;
  WT[          12] =   0.2603824067115784E+02;
  WT[          13] =   0.2805418062210083E+02;
  WT[          14] =   0.3007012057304382E+02;
  WT[          15] =   0.1703060948848724E+02;
  WT[          16] =   0.1400669956207275E+02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
void egtransetEPS(double* EPS) {
  EPS[           0] =   0.3800000000000000E+02;
  EPS[           1] =   0.1450000000000000E+03;
  EPS[           2] =   0.1074000000000000E+03;
  EPS[           3] =   0.8000000000000000E+02;
  EPS[           4] =   0.5724000000000000E+03;
  EPS[           5] =   0.1074000000000000E+03;
  EPS[           6] =   0.1074000000000000E+03;
  EPS[           7] =   0.1440000000000000E+03;
  EPS[           8] =   0.1414000000000000E+03;
  EPS[           9] =   0.9809999999999999E+02;
  EPS[          10] =   0.2440000000000000E+03;
  EPS[          11] =   0.4980000000000000E+03;
  EPS[          12] =   0.2090000000000000E+03;
  EPS[          13] =   0.2808000000000000E+03;
  EPS[          14] =   0.2523000000000000E+03;
  EPS[          15] =   0.4810000000000000E+03;
  EPS[          16] =   0.7140000000000001E+02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
void egtransetSIG(double* SIG) {
  SIG[           0] =   0.2920000000000000E+01;
  SIG[           1] =   0.2050000000000000E+01;
  SIG[           2] =   0.3458000000000000E+01;
  SIG[           3] =   0.2750000000000000E+01;
  SIG[           4] =   0.2605000000000000E+01;
  SIG[           5] =   0.3458000000000000E+01;
  SIG[           6] =   0.3458000000000000E+01;
  SIG[           7] =   0.3800000000000000E+01;
  SIG[           8] =   0.3746000000000000E+01;
  SIG[           9] =   0.3650000000000000E+01;
  SIG[          10] =   0.3763000000000000E+01;
  SIG[          11] =   0.3590000000000000E+01;
  SIG[          12] =   0.4100000000000000E+01;
  SIG[          13] =   0.3971000000000000E+01;
  SIG[          14] =   0.4302000000000000E+01;
  SIG[          15] =   0.2920000000000000E+01;
  SIG[          16] =   0.3298000000000000E+01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetDIP EGTRANSETDIP
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetDIP egtransetdip
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetDIP egtransetdip_
#endif
void egtransetDIP(double* DIP) {
  DIP[           0] =   0.0000000000000000E+00;
  DIP[           1] =   0.0000000000000000E+00;
  DIP[           2] =   0.0000000000000000E+00;
  DIP[           3] =   0.0000000000000000E+00;
  DIP[           4] =   0.1844000000000000E+01;
  DIP[           5] =   0.0000000000000000E+00;
  DIP[           6] =   0.0000000000000000E+00;
  DIP[           7] =   0.0000000000000000E+00;
  DIP[           8] =   0.0000000000000000E+00;
  DIP[           9] =   0.0000000000000000E+00;
  DIP[          10] =   0.0000000000000000E+00;
  DIP[          11] =   0.0000000000000000E+00;
  DIP[          12] =   0.0000000000000000E+00;
  DIP[          13] =   0.0000000000000000E+00;
  DIP[          14] =   0.0000000000000000E+00;
  DIP[          15] =   0.1470000000000000E+01;
  DIP[          16] =   0.0000000000000000E+00;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPOL EGTRANSETPOL
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPOL egtransetpol
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPOL egtransetpol_
#endif
void egtransetPOL(double* POL) {
  POL[           0] =   0.7900000000000000E+00;
  POL[           1] =   0.0000000000000000E+00;
  POL[           2] =   0.1600000000000000E+01;
  POL[           3] =   0.0000000000000000E+00;
  POL[           4] =   0.0000000000000000E+00;
  POL[           5] =   0.0000000000000000E+00;
  POL[           6] =   0.0000000000000000E+00;
  POL[           7] =   0.0000000000000000E+00;
  POL[           8] =   0.2600000000000000E+01;
  POL[           9] =   0.1950000000000000E+01;
  POL[          10] =   0.2650000000000000E+01;
  POL[          11] =   0.0000000000000000E+00;
  POL[          12] =   0.0000000000000000E+00;
  POL[          13] =   0.0000000000000000E+00;
  POL[          14] =   0.0000000000000000E+00;
  POL[          15] =   0.0000000000000000E+00;
  POL[          16] =   0.0000000000000000E+00;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetZROT EGTRANSETZROT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetZROT egtransetzrot
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetZROT egtransetzrot_
#endif
void egtransetZROT(double* ZROT) {
  ZROT[           0] =   0.2800000000000000E+03;
  ZROT[           1] =   0.0000000000000000E+00;
  ZROT[           2] =   0.3800000000000000E+01;
  ZROT[           3] =   0.0000000000000000E+00;
  ZROT[           4] =   0.4000000000000000E+01;
  ZROT[           5] =   0.1000000000000000E+01;
  ZROT[           6] =   0.3800000000000000E+01;
  ZROT[           7] =   0.0000000000000000E+00;
  ZROT[           8] =   0.1300000000000000E+02;
  ZROT[           9] =   0.1800000000000000E+01;
  ZROT[          10] =   0.2100000000000000E+01;
  ZROT[          11] =   0.2000000000000000E+01;
  ZROT[          12] =   0.2500000000000000E+01;
  ZROT[          13] =   0.1500000000000000E+01;
  ZROT[          14] =   0.1500000000000000E+01;
  ZROT[          15] =   0.1000000000000000E+02;
  ZROT[          16] =   0.0000000000000000E+00;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLIN EGTRANSETNLIN
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLIN egtransetnlin
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLIN egtransetnlin_
#endif
void egtransetNLIN(int* NLIN) {
  NLIN[           0] =            1;
  NLIN[           1] =            0;
  NLIN[           2] =            1;
  NLIN[           3] =            1;
  NLIN[           4] =            2;
  NLIN[           5] =            2;
  NLIN[           6] =            2;
  NLIN[           7] =            1;
  NLIN[           8] =            2;
  NLIN[           9] =            1;
  NLIN[          10] =            1;
  NLIN[          11] =            2;
  NLIN[          12] =            1;
  NLIN[          13] =            2;
  NLIN[          14] =            2;
  NLIN[          15] =            2;
  NLIN[          16] =            0;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
void egtransetCOFLAM(double* COFLAM) {
  COFLAM[           0] =   0.4338938909590841E+01;
  COFLAM[           1] =   0.1557129190255232E+01;
  COFLAM[           2] =  -0.1611482772000911E+00;
  COFLAM[           3] =   0.9924831871819461E-02;
  COFLAM[           4] =  -0.1277787249296346E+01;
  COFLAM[           5] =   0.3820054734917150E+01;
  COFLAM[           6] =  -0.4198233798779158E+00;
  COFLAM[           7] =   0.1850208398735055E-01;
  COFLAM[           8] =   0.5325069412912258E+00;
  COFLAM[           9] =   0.1869914528199901E+01;
  COFLAM[          10] =  -0.1314362832820327E+00;
  COFLAM[          11] =   0.5215120926063267E-02;
  COFLAM[          12] =   0.9093696928027398E+01;
  COFLAM[          13] =  -0.1120775592691738E+01;
  COFLAM[          14] =   0.2389873914448641E+00;
  COFLAM[          15] =  -0.9739555323602451E-02;
  COFLAM[          16] =   0.1823299080601272E+02;
  COFLAM[          17] =  -0.6774533358601863E+01;
  COFLAM[          18] =   0.1218719282523822E+01;
  COFLAM[          19] =  -0.6135010154070134E-01;
  COFLAM[          20] =   0.3326317596296849E+01;
  COFLAM[          21] =   0.4132335554940667E+00;
  COFLAM[          22] =   0.1130570823333984E+00;
  COFLAM[          23] =  -0.7338162598088820E-02;
  COFLAM[          24] =   0.2915877425653172E+01;
  COFLAM[          25] =   0.4587330741194222E+00;
  COFLAM[          26] =   0.1387858388240978E+00;
  COFLAM[          27] =  -0.9951989822150117E-02;
  COFLAM[          28] =   0.1076771824317687E+02;
  COFLAM[          29] =  -0.3230785742845311E+01;
  COFLAM[          30] =   0.7028613656761165E+00;
  COFLAM[          31] =  -0.3786794222679780E-01;
  COFLAM[          32] =   0.1770109820287010E+02;
  COFLAM[          33] =  -0.6715184244214581E+01;
  COFLAM[          34] =   0.1263352812473795E+01;
  COFLAM[          35] =  -0.6629210623948648E-01;
  COFLAM[          36] =   0.8169761711928375E+01;
  COFLAM[          37] =  -0.1536354590024332E+01;
  COFLAM[          38] =   0.3677852316201770E+00;
  COFLAM[          39] =  -0.1908151636745373E-01;
  COFLAM[          40] =  -0.8741249063434607E+01;
  COFLAM[          41] =   0.4789682737843632E+01;
  COFLAM[          42] =  -0.4182583551792970E+00;
  COFLAM[          43] =   0.1350148479415579E-01;
  COFLAM[          44] =   0.1432241064124226E+02;
  COFLAM[          45] =  -0.6063127607975618E+01;
  COFLAM[          46] =   0.1238070730356745E+01;
  COFLAM[          47] =  -0.6822324216312974E-01;
  COFLAM[          48] =  -0.1080654352351750E+02;
  COFLAM[          49] =   0.5871580561001701E+01;
  COFLAM[          50] =  -0.5858009165203976E+00;
  COFLAM[          51] =   0.2242112198055643E-01;
  COFLAM[          52] =  -0.4034386834241054E+01;
  COFLAM[          53] =   0.1906845175798028E+01;
  COFLAM[          54] =   0.1176313823060724E+00;
  COFLAM[          55] =  -0.1610854731203796E-01;
  COFLAM[          56] =  -0.2116070664089148E+01;
  COFLAM[          57] =   0.9961322379018972E+00;
  COFLAM[          58] =   0.2612385762947667E+00;
  COFLAM[          59] =  -0.2335852593974065E-01;
  COFLAM[          60] =   0.1756249168236376E+02;
  COFLAM[          61] =  -0.6844841344645484E+01;
  COFLAM[          62] =   0.1296574793469918E+01;
  COFLAM[          63] =  -0.6830286748937907E-01;
  COFLAM[          64] =   0.1334087340220970E+01;
  COFLAM[          65] =   0.1967564865885963E+01;
  COFLAM[          66] =  -0.1801563268324271E+00;
  COFLAM[          67] =   0.8161858212342093E-02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFETA EGTRANSETCOFETA
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFETA egtransetcofeta
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFETA egtransetcofeta_
#endif
void egtransetCOFETA(double* COFETA) {
  COFETA[           0] =  -0.1405365243967108E+02;
  COFETA[           1] =   0.1093067315260161E+01;
  COFETA[           2] =  -0.6261814846952836E-01;
  COFETA[           3] =   0.2905012497485487E-02;
  COFETA[           4] =  -0.2082770253853309E+02;
  COFETA[           5] =   0.3820054734917064E+01;
  COFETA[           6] =  -0.4198233798779034E+00;
  COFETA[           7] =   0.1850208398734995E-01;
  COFETA[           8] =  -0.1788787183918792E+02;
  COFETA[           9] =   0.2980942390384383E+01;
  COFETA[          10] =  -0.3137212871769931E+00;
  COFETA[          11] =   0.1402845871347659E-01;
  COFETA[          12] =  -0.1575591740228496E+02;
  COFETA[          13] =   0.2214372539700379E+01;
  COFETA[          14] =  -0.2131246920584266E+00;
  COFETA[          15] =   0.9629058498016490E-02;
  COFETA[          16] =  -0.1132071455875007E+02;
  COFETA[          17] =  -0.1024187591396996E+01;
  COFETA[          18] =   0.3672551958155785E+00;
  COFETA[          19] =  -0.2140495811549247E-01;
  COFETA[          20] =  -0.1787236469476987E+02;
  COFETA[          21] =   0.2980942390384401E+01;
  COFETA[          22] =  -0.3137212871769956E+00;
  COFETA[          23] =   0.1402845871347671E-01;
  COFETA[          24] =  -0.1785732406117402E+02;
  COFETA[          25] =   0.2980942390384431E+01;
  COFETA[          26] =  -0.3137212871770002E+00;
  COFETA[          27] =   0.1402845871347693E-01;
  COFETA[          28] =  -0.2064817961278261E+02;
  COFETA[          29] =   0.3796664365547116E+01;
  COFETA[          30] =  -0.4168591838616045E+00;
  COFETA[          31] =   0.1837669926760990E-01;
  COFETA[          32] =  -0.2044795093120063E+02;
  COFETA[          33] =   0.3746140041781780E+01;
  COFETA[          34] =  -0.4106325045512437E+00;
  COFETA[          35] =   0.1812122592865459E-01;
  COFETA[          36] =  -0.1745394886003354E+02;
  COFETA[          37] =   0.2749779734238378E+01;
  COFETA[          38] =  -0.2838005922559860E+00;
  COFETA[          39] =   0.1273750763474062E-01;
  COFETA[          40] =  -0.2279534579079101E+02;
  COFETA[          41] =   0.4622798720286126E+01;
  COFETA[          42] =  -0.4997444888565643E+00;
  COFETA[          43] =   0.2095793339101399E-01;
  COFETA[          44] =  -0.1574565551467780E+02;
  COFETA[          45] =   0.9879045622370226E+00;
  COFETA[          46] =   0.7005046117242672E-01;
  COFETA[          47] =  -0.7652880399373136E-02;
  COFETA[          48] =  -0.2285054221476539E+02;
  COFETA[          49] =   0.4573684124304381E+01;
  COFETA[          50] =  -0.5044627764440257E+00;
  COFETA[          51] =   0.2162047229656920E-01;
  COFETA[          52] =  -0.2294894266555372E+02;
  COFETA[          53] =   0.4436068808351478E+01;
  COFETA[          54] =  -0.4620857825400307E+00;
  COFETA[          55] =   0.1877687773259413E-01;
  COFETA[          56] =  -0.2324874285445497E+02;
  COFETA[          57] =   0.4594844249509825E+01;
  COFETA[          58] =  -0.4931221039570846E+00;
  COFETA[          59] =   0.2054776726722393E-01;
  COFETA[          60] =  -0.1419898917157052E+02;
  COFETA[          61] =   0.3629587083500767E+00;
  COFETA[          62] =   0.1568120182105773E+00;
  COFETA[          63] =  -0.1151267609047790E-01;
  COFETA[          64] =  -0.1558423057668358E+02;
  COFETA[          65] =   0.1967564865886029E+01;
  COFETA[          66] =  -0.1801563268324368E+00;
  COFETA[          67] =   0.8161858212342565E-02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFD EGTRANSETCOFD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFD egtransetcofd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFD egtransetcofd_
#endif
void egtransetCOFD(double* COFD) {
  COFD[           0] =  -0.1041995820481732E+02;
  COFD[           1] =   0.2234441574336671E+01;
  COFD[           2] =  -0.8113717628974293E-01;
  COFD[           3] =   0.3771369728593280E-02;
  COFD[           4] =  -0.1246996149641938E+02;
  COFD[           5] =   0.3218029211436469E+01;
  COFD[           6] =  -0.2110489360993846E+00;
  COFD[           7] =   0.9480877224339204E-02;
  COFD[           8] =  -0.1309173178181482E+02;
  COFD[           9] =   0.3085278377844914E+01;
  COFD[          10] =  -0.1946824097426824E+00;
  COFD[          11] =   0.8807381240537916E-02;
  COFD[          12] =  -0.1175032758512814E+02;
  COFD[          13] =   0.2642045353939402E+01;
  COFD[          14] =  -0.1332637450759486E+00;
  COFD[          15] =   0.5973546569611733E-02;
  COFD[          16] =  -0.1757366663670086E+02;
  COFD[          17] =   0.4778618797909657E+01;
  COFD[          18] =  -0.3996041319366184E+00;
  COFD[          19] =   0.1705904378841811E-01;
  COFD[          20] =  -0.1309613046501959E+02;
  COFD[          21] =   0.3086801219615145E+01;
  COFD[          22] =  -0.1948957472387874E+00;
  COFD[          23] =   0.8817294833761086E-02;
  COFD[          24] =  -0.1310029609771552E+02;
  COFD[          25] =   0.3088245187916686E+01;
  COFD[          26] =  -0.1950980328845662E+00;
  COFD[          27] =   0.8826694775135658E-02;
  COFD[          28] =  -0.1401582227027964E+02;
  COFD[          29] =   0.3415426664334749E+01;
  COFD[          30] =  -0.2385853115052575E+00;
  COFD[          31] =   0.1075469512610122E-01;
  COFD[          32] =  -0.1394568249588876E+02;
  COFD[          33] =   0.3392701225291731E+01;
  COFD[          34] =  -0.2355263017693934E+00;
  COFD[          35] =   0.1061747669734484E-01;
  COFD[          36] =  -0.1288505703505633E+02;
  COFD[          37] =   0.2980129383938292E+01;
  COFD[          38] =  -0.1801545917674105E+00;
  COFD[          39] =   0.8138117217492757E-02;
  COFD[          40] =  -0.1558361674716745E+02;
  COFD[          41] =   0.3987098568284434E+01;
  COFD[          42] =  -0.3109866036899985E+00;
  COFD[          43] =   0.1380795341338578E-01;
  COFD[          44] =  -0.1688564690991911E+02;
  COFD[          45] =   0.4415865473479337E+01;
  COFD[          46] =  -0.3564104485793251E+00;
  COFD[          47] =   0.1533855847922275E-01;
  COFD[          48] =  -0.1548207029093904E+02;
  COFD[          49] =   0.3939084188623417E+01;
  COFD[          50] =  -0.3072761108436687E+00;
  COFD[          51] =   0.1375770270895890E-01;
  COFD[          52] =  -0.1606152931877502E+02;
  COFD[          53] =   0.4146430359164761E+01;
  COFD[          54] =  -0.3312210923439469E+00;
  COFD[          55] =   0.1466430884967795E-01;
  COFD[          56] =  -0.1596236229610858E+02;
  COFD[          57] =   0.4079842750882008E+01;
  COFD[          58] =  -0.3231960602071978E+00;
  COFD[          59] =   0.1434305777680907E-01;
  COFD[          60] =  -0.1682806034335371E+02;
  COFD[          61] =   0.4485450500929320E+01;
  COFD[          62] =  -0.3655695289846634E+00;
  COFD[          63] =   0.1574614737371121E-01;
  COFD[          64] =  -0.1163276178298392E+02;
  COFD[          65] =   0.2530205071288130E+01;
  COFD[          66] =  -0.1183776350934922E+00;
  COFD[          67] =   0.5312988352557955E-02;
  COFD[          68] =  -0.1246996149641938E+02;
  COFD[          69] =   0.3218029211436469E+01;
  COFD[          70] =  -0.2110489360993846E+00;
  COFD[          71] =   0.9480877224339204E-02;
  COFD[          72] =  -0.1523325703608638E+02;
  COFD[          73] =   0.4393426193693418E+01;
  COFD[          74] =  -0.3557313428561681E+00;
  COFD[          75] =   0.1541554803417400E-01;
  COFD[          76] =  -0.1615679289714481E+02;
  COFD[          77] =   0.4416835307133869E+01;
  COFD[          78] =  -0.3582737234170646E+00;
  COFD[          79] =   0.1548663572743747E-01;
  COFD[          80] =  -0.1489056810719442E+02;
  COFD[          81] =   0.4079672697171075E+01;
  COFD[          82] =  -0.3197193422647769E+00;
  COFD[          83] =   0.1404463618351796E-01;
  COFD[          84] =  -0.1615075966577466E+02;
  COFD[          85] =   0.4076548482999812E+01;
  COFD[          86] =  -0.2637510941604614E+00;
  COFD[          87] =   0.9371841080052064E-02;
  COFD[          88] =  -0.1616058613782880E+02;
  COFD[          89] =   0.4418087193226073E+01;
  COFD[          90] =  -0.3584223903851448E+00;
  COFD[          91] =   0.1549243249303281E-01;
  COFD[          92] =  -0.1616421015670871E+02;
  COFD[          93] =   0.4419285820914456E+01;
  COFD[          94] =  -0.3585647608630650E+00;
  COFD[          95] =   0.1549798503551403E-01;
  COFD[          96] =  -0.1668220221378553E+02;
  COFD[          97] =   0.4517665522196412E+01;
  COFD[          98] =  -0.3643177404654742E+00;
  COFD[          99] =   0.1544861119982853E-01;
  COFD[         100] =  -0.1661124059439949E+02;
  COFD[         101] =   0.4499983815327449E+01;
  COFD[         102] =  -0.3624203284397300E+00;
  COFD[         103] =   0.1538350856590451E-01;
  COFD[         104] =  -0.1640877269708571E+02;
  COFD[         105] =   0.4519498926662536E+01;
  COFD[         106] =  -0.3750329779403845E+00;
  COFD[         107] =   0.1635855421757920E-01;
  COFD[         108] =  -0.1664237108262573E+02;
  COFD[         109] =   0.4321207774376617E+01;
  COFD[         110] =  -0.3181315549506706E+00;
  COFD[         111] =   0.1258308183870717E-01;
  COFD[         112] =  -0.1567567878102999E+02;
  COFD[         113] =   0.3706302510623138E+01;
  COFD[         114] =  -0.2096236244277931E+00;
  COFD[         115] =   0.6748353057110522E-02;
  COFD[         116] =  -0.1685361380845761E+02;
  COFD[         117] =   0.4411251776453823E+01;
  COFD[         118] =  -0.3350756131548513E+00;
  COFD[         119] =   0.1353157937321782E-01;
  COFD[         120] =  -0.1635823810097287E+02;
  COFD[         121] =   0.4123977007004703E+01;
  COFD[         122] =  -0.2854728908229589E+00;
  COFD[         123] =   0.1088433914065527E-01;
  COFD[         124] =  -0.1646195534760534E+02;
  COFD[         125] =   0.4141822899037704E+01;
  COFD[         126] =  -0.2892308013569180E+00;
  COFD[         127] =   0.1109614462742981E-01;
  COFD[         128] =  -0.1626359848098758E+02;
  COFD[         129] =   0.4116788683843708E+01;
  COFD[         130] =  -0.2733289330957227E+00;
  COFD[         131] =   0.9943035956817482E-02;
  COFD[         132] =  -0.1522605770447405E+02;
  COFD[         133] =   0.4154863253653862E+01;
  COFD[         134] =  -0.3321404775716639E+00;
  COFD[         135] =   0.1469592836326692E-01;
  COFD[         136] =  -0.1309173178181482E+02;
  COFD[         137] =   0.3085278377844914E+01;
  COFD[         138] =  -0.1946824097426824E+00;
  COFD[         139] =   0.8807381240537916E-02;
  COFD[         140] =  -0.1615679289714481E+02;
  COFD[         141] =   0.4416835307133869E+01;
  COFD[         142] =  -0.3582737234170646E+00;
  COFD[         143] =   0.1548663572743747E-01;
  COFD[         144] =  -0.1633859102926806E+02;
  COFD[         145] =   0.3804061718963007E+01;
  COFD[         146] =  -0.2844545503177112E+00;
  COFD[         147] =   0.1254577088876504E-01;
  COFD[         148] =  -0.1549441379638191E+02;
  COFD[         149] =   0.3661937775617350E+01;
  COFD[         150] =  -0.2686571425887965E+00;
  COFD[         151] =   0.1197472750180553E-01;
  COFD[         152] =  -0.1922081039841695E+02;
  COFD[         153] =   0.4711535347815869E+01;
  COFD[         154] =  -0.3619579159000997E+00;
  COFD[         155] =   0.1426886476331241E-01;
  COFD[         156] =  -0.1634649054694938E+02;
  COFD[         157] =   0.3804144246517877E+01;
  COFD[         158] =  -0.2844651891354975E+00;
  COFD[         159] =   0.1254622626209081E-01;
  COFD[         160] =  -0.1635443062264980E+02;
  COFD[         161] =   0.3804381702146611E+01;
  COFD[         162] =  -0.2844957999638338E+00;
  COFD[         163] =   0.1254753649276844E-01;
  COFD[         164] =  -0.1710012222365139E+02;
  COFD[         165] =   0.4123808736773267E+01;
  COFD[         166] =  -0.3236376370297240E+00;
  COFD[         167] =   0.1414381001550686E-01;
  COFD[         168] =  -0.1705130515010214E+02;
  COFD[         169] =   0.4103994885241071E+01;
  COFD[         170] =  -0.3212499365446575E+00;
  COFD[         171] =   0.1404820978963865E-01;
  COFD[         172] =  -0.1615768146554276E+02;
  COFD[         173] =   0.3734150883385773E+01;
  COFD[         174] =  -0.2760989842105895E+00;
  COFD[         175] =   0.1221438663500826E-01;
  COFD[         176] =  -0.1871735010457793E+02;
  COFD[         177] =   0.4573421885700214E+01;
  COFD[         178] =  -0.3757214086933051E+00;
  COFD[         179] =   0.1614451502251515E-01;
  COFD[         180] =  -0.1977364410044468E+02;
  COFD[         181] =   0.4876531592391508E+01;
  COFD[         182] =  -0.3981741984755758E+00;
  COFD[         183] =   0.1643476290195535E-01;
  COFD[         184] =  -0.1830650296367939E+02;
  COFD[         185] =   0.4452394507404697E+01;
  COFD[         186] =  -0.3624991200767191E+00;
  COFD[         187] =   0.1567286266089603E-01;
  COFD[         188] =  -0.1904759289384530E+02;
  COFD[         189] =   0.4700929573683383E+01;
  COFD[         190] =  -0.3903608191209737E+00;
  COFD[         191] =   0.1670229359855150E-01;
  COFD[         192] =  -0.1883799365221673E+02;
  COFD[         193] =   0.4594611640098751E+01;
  COFD[         194] =  -0.3781168670160998E+00;
  COFD[         195] =   0.1623467419613562E-01;
  COFD[         196] =  -0.1934710272078301E+02;
  COFD[         197] =   0.4822779400529180E+01;
  COFD[         198] =  -0.3871228449439424E+00;
  COFD[         199] =   0.1578686341252481E-01;
  COFD[         200] =  -0.1527049784849730E+02;
  COFD[         201] =   0.3538365317994180E+01;
  COFD[         202] =  -0.2526814234121028E+00;
  COFD[         203] =   0.1128573168209360E-01;
  COFD[         204] =  -0.1175032758512814E+02;
  COFD[         205] =   0.2642045353939402E+01;
  COFD[         206] =  -0.1332637450759486E+00;
  COFD[         207] =   0.5973546569611733E-02;
  COFD[         208] =  -0.1489056810719442E+02;
  COFD[         209] =   0.4079672697171075E+01;
  COFD[         210] =  -0.3197193422647769E+00;
  COFD[         211] =   0.1404463618351796E-01;
  COFD[         212] =  -0.1549441379638191E+02;
  COFD[         213] =   0.3661937775617350E+01;
  COFD[         214] =  -0.2686571425887965E+00;
  COFD[         215] =   0.1197472750180553E-01;
  COFD[         216] =  -0.1408727594805337E+02;
  COFD[         217] =   0.3263692027745076E+01;
  COFD[         218] =  -0.2164586669743292E+00;
  COFD[         219] =   0.9693817632062094E-02;
  COFD[         220] =  -0.1876923157317771E+02;
  COFD[         221] =   0.4853854184816717E+01;
  COFD[         222] =  -0.3998497982364926E+00;
  COFD[         223] =   0.1669227468405473E-01;
  COFD[         224] =  -0.1551249966394157E+02;
  COFD[         225] =   0.3667301435894143E+01;
  COFD[         226] =  -0.2693892845913942E+00;
  COFD[         227] =   0.1200791993248851E-01;
  COFD[         228] =  -0.1553013778710249E+02;
  COFD[         229] =   0.3672587405633829E+01;
  COFD[         230] =  -0.2701107938194129E+00;
  COFD[         231] =   0.1204062933869048E-01;
  COFD[         232] =  -0.1596628118290846E+02;
  COFD[         233] =   0.3841273000042274E+01;
  COFD[         234] =  -0.2898297027626218E+00;
  COFD[         235] =   0.1280193859902638E-01;
  COFD[         236] =  -0.1591206444980951E+02;
  COFD[         237] =   0.3820917505869285E+01;
  COFD[         238] =  -0.2872227297192593E+00;
  COFD[         239] =   0.1269037691099738E-01;
  COFD[         240] =  -0.1524511379647520E+02;
  COFD[         241] =   0.3553053597674808E+01;
  COFD[         242] =  -0.2545212965294182E+00;
  COFD[         243] =   0.1136271070888101E-01;
  COFD[         244] =  -0.1768659705517047E+02;
  COFD[         245] =   0.4367905899188474E+01;
  COFD[         246] =  -0.3516645879639785E+00;
  COFD[         247] =   0.1520186726409712E-01;
  COFD[         248] =  -0.1891518602118661E+02;
  COFD[         249] =   0.4756739681751712E+01;
  COFD[         250] =  -0.3886685070061102E+00;
  COFD[         251] =   0.1625748441332711E-01;
  COFD[         252] =  -0.1730907989620245E+02;
  COFD[         253] =   0.4243749180870270E+01;
  COFD[         254] =  -0.3387466685616417E+00;
  COFD[         255] =   0.1477905643158847E-01;
  COFD[         256] =  -0.1795985745949001E+02;
  COFD[         257] =   0.4461091096847464E+01;
  COFD[         258] =  -0.3624212969270659E+00;
  COFD[         259] =   0.1561553586311572E-01;
  COFD[         260] =  -0.1782798069130013E+02;
  COFD[         261] =   0.4381516462622495E+01;
  COFD[         262] =  -0.3532970918753254E+00;
  COFD[         263] =   0.1526804966457564E-01;
  COFD[         264] =  -0.1863707722655403E+02;
  COFD[         265] =   0.4806992114833483E+01;
  COFD[         266] =  -0.3984196894018535E+00;
  COFD[         267] =   0.1681856268295304E-01;
  COFD[         268] =  -0.1395220381010925E+02;
  COFD[         269] =   0.3159688299904867E+01;
  COFD[         270] =  -0.2028455188834762E+00;
  COFD[         271] =   0.9099190221367019E-02;
  COFD[         272] =  -0.1757366663670086E+02;
  COFD[         273] =   0.4778618797909657E+01;
  COFD[         274] =  -0.3996041319366184E+00;
  COFD[         275] =   0.1705904378841811E-01;
  COFD[         276] =  -0.1615075966577466E+02;
  COFD[         277] =   0.4076548482999812E+01;
  COFD[         278] =  -0.2637510941604614E+00;
  COFD[         279] =   0.9371841080052064E-02;
  COFD[         280] =  -0.1922081039841695E+02;
  COFD[         281] =   0.4711535347815869E+01;
  COFD[         282] =  -0.3619579159000997E+00;
  COFD[         283] =   0.1426886476331241E-01;
  COFD[         284] =  -0.1876923157317771E+02;
  COFD[         285] =   0.4853854184816717E+01;
  COFD[         286] =  -0.3998497982364926E+00;
  COFD[         287] =   0.1669227468405473E-01;
  COFD[         288] =  -0.1159394277352764E+02;
  COFD[         289] =   0.8274415074083117E+00;
  COFD[         290] =   0.2507810870995229E+00;
  COFD[         291] =  -0.1608920912964349E-01;
  COFD[         292] =  -0.1918278409127273E+02;
  COFD[         293] =   0.4757442212360504E+01;
  COFD[         294] =  -0.3754949936858373E+00;
  COFD[         295] =   0.1515864657024127E-01;
  COFD[         296] =  -0.1917291211841810E+02;
  COFD[         297] =   0.4750317181782021E+01;
  COFD[         298] =  -0.3743997734248102E+00;
  COFD[         299] =   0.1510405742415517E-01;
  COFD[         300] =  -0.1924307146615256E+02;
  COFD[         301] =   0.4719749979653547E+01;
  COFD[         302] =  -0.3616169791891675E+00;
  COFD[         303] =   0.1420601032870719E-01;
  COFD[         304] =  -0.1899850780979201E+02;
  COFD[         305] =   0.4543165455349725E+01;
  COFD[         306] =  -0.3288657155661528E+00;
  COFD[         307] =   0.1239537851241695E-01;
  COFD[         308] =  -0.1935427072301072E+02;
  COFD[         309] =   0.4779618521455131E+01;
  COFD[         310] =  -0.3741503724582112E+00;
  COFD[         311] =   0.1493019768764606E-01;
  COFD[         312] =  -0.1719414340129819E+02;
  COFD[         313] =   0.3490298819499314E+01;
  COFD[         314] =  -0.1625329575445097E+00;
  COFD[         315] =   0.4043757780510787E-02;
  COFD[         316] =  -0.1468058225757294E+02;
  COFD[         317] =   0.2298576535112391E+01;
  COFD[         318] =   0.1810975986268222E-01;
  COFD[         319] =  -0.4731609614535192E-02;
  COFD[         320] =  -0.1868831150678477E+02;
  COFD[         321] =   0.4247585976166683E+01;
  COFD[         322] =  -0.2815735628073077E+00;
  COFD[         323] =   0.1000602366593966E-01;
  COFD[         324] =  -0.1767910948230725E+02;
  COFD[         325] =   0.3722086580781685E+01;
  COFD[         326] =  -0.1980254146986777E+00;
  COFD[         327] =   0.5799038224028695E-02;
  COFD[         328] =  -0.1813805288510309E+02;
  COFD[         329] =   0.3907124675584541E+01;
  COFD[         330] =  -0.2273991439446390E+00;
  COFD[         331] =   0.7280601459305962E-02;
  COFD[         332] =  -0.1316791411191276E+02;
  COFD[         333] =   0.1635596276910710E+01;
  COFD[         334] =   0.1248562922809231E+00;
  COFD[         335] =  -0.1002143991895203E-01;
  COFD[         336] =  -0.1875637020961290E+02;
  COFD[         337] =   0.4821872691746870E+01;
  COFD[         338] =  -0.3987313751514457E+00;
  COFD[         339] =   0.1676685214877852E-01;
  COFD[         340] =  -0.1309613046501959E+02;
  COFD[         341] =   0.3086801219615145E+01;
  COFD[         342] =  -0.1948957472387874E+00;
  COFD[         343] =   0.8817294833761086E-02;
  COFD[         344] =  -0.1616058613782880E+02;
  COFD[         345] =   0.4418087193226073E+01;
  COFD[         346] =  -0.3584223903851448E+00;
  COFD[         347] =   0.1549243249303281E-01;
  COFD[         348] =  -0.1634649054694938E+02;
  COFD[         349] =   0.3804144246517877E+01;
  COFD[         350] =  -0.2844651891354975E+00;
  COFD[         351] =   0.1254622626209081E-01;
  COFD[         352] =  -0.1551249966394157E+02;
  COFD[         353] =   0.3667301435894143E+01;
  COFD[         354] =  -0.2693892845913942E+00;
  COFD[         355] =   0.1200791993248851E-01;
  COFD[         356] =  -0.1918278409127273E+02;
  COFD[         357] =   0.4757442212360504E+01;
  COFD[         358] =  -0.3754949936858373E+00;
  COFD[         359] =   0.1515864657024127E-01;
  COFD[         360] =  -0.1635409817368614E+02;
  COFD[         361] =   0.3804061718963002E+01;
  COFD[         362] =  -0.2844545503177105E+00;
  COFD[         363] =   0.1254577088876500E-01;
  COFD[         364] =  -0.1636175579013954E+02;
  COFD[         365] =   0.3804139357106159E+01;
  COFD[         366] =  -0.2844645588307625E+00;
  COFD[         367] =   0.1254619928317991E-01;
  COFD[         368] =  -0.1710991586467964E+02;
  COFD[         369] =   0.4125646326528441E+01;
  COFD[         370] =  -0.3238554740258663E+00;
  COFD[         371] =   0.1415227419625766E-01;
  COFD[         372] =  -0.1706120101779605E+02;
  COFD[         373] =   0.4105791991868936E+01;
  COFD[         374] =  -0.3214635787029052E+00;
  COFD[         375] =   0.1405654084975138E-01;
  COFD[         376] =  -0.1616672146196606E+02;
  COFD[         377] =   0.3734919571093601E+01;
  COFD[         378] =  -0.2762018861637152E+00;
  COFD[         379] =   0.1221896274296108E-01;
  COFD[         380] =  -0.1872670165824808E+02;
  COFD[         381] =   0.4573903479896379E+01;
  COFD[         382] =  -0.3758340786206700E+00;
  COFD[         383] =   0.1615161324401359E-01;
  COFD[         384] =  -0.1978049689264223E+02;
  COFD[         385] =   0.4876281709421428E+01;
  COFD[         386] =  -0.3981401304053240E+00;
  COFD[         387] =   0.1643325343977014E-01;
  COFD[         388] =  -0.1831322697803230E+02;
  COFD[         389] =   0.4452261161272578E+01;
  COFD[         390] =  -0.3624715246876605E+00;
  COFD[         391] =   0.1567121551911240E-01;
  COFD[         392] =  -0.1905419842832635E+02;
  COFD[         393] =   0.4700662819460902E+01;
  COFD[         394] =  -0.3903206776226419E+00;
  COFD[         395] =   0.1670032573134212E-01;
  COFD[         396] =  -0.1884529530697174E+02;
  COFD[         397] =   0.4594634055023398E+01;
  COFD[         398] =  -0.3781333057516071E+00;
  COFD[         399] =   0.1623599168682821E-01;
  COFD[         400] =  -0.1919185322703664E+02;
  COFD[         401] =   0.4784033288807873E+01;
  COFD[         402] =  -0.3849459064100618E+00;
  COFD[         403] =   0.1580648990064550E-01;
  COFD[         404] =  -0.1528586419146581E+02;
  COFD[         405] =   0.3542872339374690E+01;
  COFD[         406] =  -0.2532972672017028E+00;
  COFD[         407] =   0.1131367837559019E-01;
  COFD[         408] =  -0.1310029609771552E+02;
  COFD[         409] =   0.3088245187916686E+01;
  COFD[         410] =  -0.1950980328845662E+00;
  COFD[         411] =   0.8826694775135658E-02;
  COFD[         412] =  -0.1616421015670871E+02;
  COFD[         413] =   0.4419285820914456E+01;
  COFD[         414] =  -0.3585647608630650E+00;
  COFD[         415] =   0.1549798503551403E-01;
  COFD[         416] =  -0.1635443062264980E+02;
  COFD[         417] =   0.3804381702146611E+01;
  COFD[         418] =  -0.2844957999638338E+00;
  COFD[         419] =   0.1254753649276844E-01;
  COFD[         420] =  -0.1553013778710249E+02;
  COFD[         421] =   0.3672587405633829E+01;
  COFD[         422] =  -0.2701107938194129E+00;
  COFD[         423] =   0.1204062933869048E-01;
  COFD[         424] =  -0.1917291211841810E+02;
  COFD[         425] =   0.4750317181782021E+01;
  COFD[         426] =  -0.3743997734248102E+00;
  COFD[         427] =   0.1510405742415517E-01;
  COFD[         428] =  -0.1636175579013954E+02;
  COFD[         429] =   0.3804139357106159E+01;
  COFD[         430] =  -0.2844645588307625E+00;
  COFD[         431] =   0.1254619928317991E-01;
  COFD[         432] =  -0.1636913880728208E+02;
  COFD[         433] =   0.3804061718963007E+01;
  COFD[         434] =  -0.2844545503177112E+00;
  COFD[         435] =   0.1254577088876504E-01;
  COFD[         436] =  -0.1711940058927942E+02;
  COFD[         437] =   0.4127460426680499E+01;
  COFD[         438] =  -0.3240705310922555E+00;
  COFD[         439] =   0.1416063080699452E-01;
  COFD[         440] =  -0.1707080539569448E+02;
  COFD[         441] =   0.4107575519112949E+01;
  COFD[         442] =  -0.3216756148320791E+00;
  COFD[         443] =   0.1406480989264605E-01;
  COFD[         444] =  -0.1617584857259722E+02;
  COFD[         445] =   0.3735860038531880E+01;
  COFD[         446] =  -0.2763277124196010E+00;
  COFD[         447] =   0.1222455543181157E-01;
  COFD[         448] =  -0.1873560076291526E+02;
  COFD[         449] =   0.4574325021189105E+01;
  COFD[         450] =  -0.3759338652876284E+00;
  COFD[         451] =   0.1615792922437573E-01;
  COFD[         452] =  -0.1978644159772173E+02;
  COFD[         453] =   0.4875758120324897E+01;
  COFD[         454] =  -0.3980631994813040E+00;
  COFD[         455] =   0.1642957450530744E-01;
  COFD[         456] =  -0.1831968967523493E+02;
  COFD[         457] =   0.4452134774757628E+01;
  COFD[         458] =  -0.3624420830935042E+00;
  COFD[         459] =   0.1566936578675393E-01;
  COFD[         460] =  -0.1906036212420135E+02;
  COFD[         461] =   0.4700323443044125E+01;
  COFD[         462] =  -0.3902665913786998E+00;
  COFD[         463] =   0.1669755465043006E-01;
  COFD[         464] =  -0.1885218260686505E+02;
  COFD[         465] =   0.4594596529989258E+01;
  COFD[         466] =  -0.3781371531513473E+00;
  COFD[         467] =   0.1623655159949430E-01;
  COFD[         468] =  -0.1918517200292939E+02;
  COFD[         469] =   0.4778400021399422E+01;
  COFD[         470] =  -0.3840658172050402E+00;
  COFD[         471] =   0.1576203256152660E-01;
  COFD[         472] =  -0.1530075071186979E+02;
  COFD[         473] =   0.3547276678457322E+01;
  COFD[         474] =  -0.2538990508478524E+00;
  COFD[         475] =   0.1134098602286849E-01;
  COFD[         476] =  -0.1401582227027964E+02;
  COFD[         477] =   0.3415426664334749E+01;
  COFD[         478] =  -0.2385853115052575E+00;
  COFD[         479] =   0.1075469512610122E-01;
  COFD[         480] =  -0.1668220221378553E+02;
  COFD[         481] =   0.4517665522196412E+01;
  COFD[         482] =  -0.3643177404654742E+00;
  COFD[         483] =   0.1544861119982853E-01;
  COFD[         484] =  -0.1710012222365139E+02;
  COFD[         485] =   0.4123808736773267E+01;
  COFD[         486] =  -0.3236376370297240E+00;
  COFD[         487] =   0.1414381001550686E-01;
  COFD[         488] =  -0.1596628118290846E+02;
  COFD[         489] =   0.3841273000042274E+01;
  COFD[         490] =  -0.2898297027626218E+00;
  COFD[         491] =   0.1280193859902638E-01;
  COFD[         492] =  -0.1924307146615256E+02;
  COFD[         493] =   0.4719749979653547E+01;
  COFD[         494] =  -0.3616169791891675E+00;
  COFD[         495] =   0.1420601032870719E-01;
  COFD[         496] =  -0.1710991586467964E+02;
  COFD[         497] =   0.4125646326528441E+01;
  COFD[         498] =  -0.3238554740258663E+00;
  COFD[         499] =   0.1415227419625766E-01;
  COFD[         500] =  -0.1711940058927942E+02;
  COFD[         501] =   0.4127460426680499E+01;
  COFD[         502] =  -0.3240705310922555E+00;
  COFD[         503] =   0.1416063080699452E-01;
  COFD[         504] =  -0.1776945810612460E+02;
  COFD[         505] =   0.4375557614241441E+01;
  COFD[         506] =  -0.3535212257934470E+00;
  COFD[         507] =   0.1532421336526862E-01;
  COFD[         508] =  -0.1770845512212784E+02;
  COFD[         509] =   0.4353002298484032E+01;
  COFD[         510] =  -0.3507494100853503E+00;
  COFD[         511] =   0.1521052013931193E-01;
  COFD[         512] =  -0.1695356550477676E+02;
  COFD[         513] =   0.4068213336589554E+01;
  COFD[         514] =  -0.3177710408064736E+00;
  COFD[         515] =   0.1394701157712223E-01;
  COFD[         516] =  -0.1913993542632300E+02;
  COFD[         517] =   0.4730620081114692E+01;
  COFD[         518] =  -0.3883890489649657E+00;
  COFD[         519] =   0.1637150502726411E-01;
  COFD[         520] =  -0.1958734834972716E+02;
  COFD[         521] =   0.4743227463757541E+01;
  COFD[         522] =  -0.3686474607497347E+00;
  COFD[         523] =   0.1466232449273593E-01;
  COFD[         524] =  -0.1892132288088230E+02;
  COFD[         525] =   0.4686257416133524E+01;
  COFD[         526] =  -0.3876178620425971E+00;
  COFD[         527] =   0.1654627973246352E-01;
  COFD[         528] =  -0.1934401530563503E+02;
  COFD[         529] =   0.4788242539311199E+01;
  COFD[         530] =  -0.3933254249328604E+00;
  COFD[         531] =   0.1648505821153991E-01;
  COFD[         532] =  -0.1926514266945471E+02;
  COFD[         533] =   0.4743728688652060E+01;
  COFD[         534] =  -0.3899270700585402E+00;
  COFD[         535] =   0.1643466452124135E-01;
  COFD[         536] =  -0.1944043258564383E+02;
  COFD[         537] =   0.4828743160693618E+01;
  COFD[         538] =  -0.3828021863328961E+00;
  COFD[         539] =   0.1540186987470838E-01;
  COFD[         540] =  -0.1573476318984406E+02;
  COFD[         541] =   0.3712839423735899E+01;
  COFD[         542] =  -0.2734220687823376E+00;
  COFD[         543] =   0.1210212757075583E-01;
  COFD[         544] =  -0.1394568249588876E+02;
  COFD[         545] =   0.3392701225291731E+01;
  COFD[         546] =  -0.2355263017693934E+00;
  COFD[         547] =   0.1061747669734484E-01;
  COFD[         548] =  -0.1661124059439949E+02;
  COFD[         549] =   0.4499983815327449E+01;
  COFD[         550] =  -0.3624203284397300E+00;
  COFD[         551] =   0.1538350856590451E-01;
  COFD[         552] =  -0.1705130515010214E+02;
  COFD[         553] =   0.4103994885241071E+01;
  COFD[         554] =  -0.3212499365446575E+00;
  COFD[         555] =   0.1404820978963865E-01;
  COFD[         556] =  -0.1591206444980951E+02;
  COFD[         557] =   0.3820917505869285E+01;
  COFD[         558] =  -0.2872227297192593E+00;
  COFD[         559] =   0.1269037691099738E-01;
  COFD[         560] =  -0.1899850780979201E+02;
  COFD[         561] =   0.4543165455349725E+01;
  COFD[         562] =  -0.3288657155661528E+00;
  COFD[         563] =   0.1239537851241695E-01;
  COFD[         564] =  -0.1706120101779605E+02;
  COFD[         565] =   0.4105791991868936E+01;
  COFD[         566] =  -0.3214635787029052E+00;
  COFD[         567] =   0.1405654084975138E-01;
  COFD[         568] =  -0.1707080539569448E+02;
  COFD[         569] =   0.4107575519112949E+01;
  COFD[         570] =  -0.3216756148320791E+00;
  COFD[         571] =   0.1406480989264605E-01;
  COFD[         572] =  -0.1770845512212784E+02;
  COFD[         573] =   0.4353002298484032E+01;
  COFD[         574] =  -0.3507494100853503E+00;
  COFD[         575] =   0.1521052013931193E-01;
  COFD[         576] =  -0.1766523778708329E+02;
  COFD[         577] =   0.4337952155636497E+01;
  COFD[         578] =  -0.3490538655747895E+00;
  COFD[         579] =   0.1514793254733080E-01;
  COFD[         580] =  -0.1690087953694522E+02;
  COFD[         581] =   0.4046945935955795E+01;
  COFD[         582] =  -0.3151742473102278E+00;
  COFD[         583] =   0.1384162361172659E-01;
  COFD[         584] =  -0.1913011857829392E+02;
  COFD[         585] =   0.4728539523227080E+01;
  COFD[         586] =  -0.3887605628634140E+00;
  COFD[         587] =   0.1641451966125198E-01;
  COFD[         588] =  -0.1962180466510501E+02;
  COFD[         589] =   0.4762356031601866E+01;
  COFD[         590] =  -0.3720986023672735E+00;
  COFD[         591] =   0.1485121157223651E-01;
  COFD[         592] =  -0.1887963312848053E+02;
  COFD[         593] =   0.4671422222663214E+01;
  COFD[         594] =  -0.3861362810478654E+00;
  COFD[         595] =   0.1650068107616750E-01;
  COFD[         596] =  -0.1933004429336219E+02;
  COFD[         597] =   0.4785693509382203E+01;
  COFD[         598] =  -0.3936420075604487E+00;
  COFD[         599] =   0.1652573779841039E-01;
  COFD[         600] =  -0.1926485085357470E+02;
  COFD[         601] =   0.4746666060969545E+01;
  COFD[         602] =  -0.3910387566978661E+00;
  COFD[         603] =   0.1651352760050653E-01;
  COFD[         604] =  -0.1941730561245821E+02;
  COFD[         605] =   0.4788578363438019E+01;
  COFD[         606] =  -0.3737994268951498E+00;
  COFD[         607] =   0.1486050997611076E-01;
  COFD[         608] =  -0.1569001110440693E+02;
  COFD[         609] =   0.3696584277188425E+01;
  COFD[         610] =  -0.2713867054547939E+00;
  COFD[         611] =   0.1201712146112774E-01;
  COFD[         612] =  -0.1288505703505633E+02;
  COFD[         613] =   0.2980129383938292E+01;
  COFD[         614] =  -0.1801545917674105E+00;
  COFD[         615] =   0.8138117217492757E-02;
  COFD[         616] =  -0.1640877269708571E+02;
  COFD[         617] =   0.4519498926662536E+01;
  COFD[         618] =  -0.3750329779403845E+00;
  COFD[         619] =   0.1635855421757920E-01;
  COFD[         620] =  -0.1615768146554276E+02;
  COFD[         621] =   0.3734150883385773E+01;
  COFD[         622] =  -0.2760989842105895E+00;
  COFD[         623] =   0.1221438663500826E-01;
  COFD[         624] =  -0.1524511379647520E+02;
  COFD[         625] =   0.3553053597674808E+01;
  COFD[         626] =  -0.2545212965294182E+00;
  COFD[         627] =   0.1136271070888101E-01;
  COFD[         628] =  -0.1935427072301072E+02;
  COFD[         629] =   0.4779618521455131E+01;
  COFD[         630] =  -0.3741503724582112E+00;
  COFD[         631] =   0.1493019768764606E-01;
  COFD[         632] =  -0.1616672146196606E+02;
  COFD[         633] =   0.3734919571093601E+01;
  COFD[         634] =  -0.2762018861637152E+00;
  COFD[         635] =   0.1221896274296108E-01;
  COFD[         636] =  -0.1617584857259722E+02;
  COFD[         637] =   0.3735860038531880E+01;
  COFD[         638] =  -0.2763277124196010E+00;
  COFD[         639] =   0.1222455543181157E-01;
  COFD[         640] =  -0.1695356550477676E+02;
  COFD[         641] =   0.4068213336589554E+01;
  COFD[         642] =  -0.3177710408064736E+00;
  COFD[         643] =   0.1394701157712223E-01;
  COFD[         644] =  -0.1690087953694522E+02;
  COFD[         645] =   0.4046945935955795E+01;
  COFD[         646] =  -0.3151742473102278E+00;
  COFD[         647] =   0.1384162361172659E-01;
  COFD[         648] =  -0.1594065528219861E+02;
  COFD[         649] =   0.3646923650815759E+01;
  COFD[         650] =  -0.2651273478800260E+00;
  COFD[         651] =   0.1175431470196917E-01;
  COFD[         652] =  -0.1856069485825941E+02;
  COFD[         653] =   0.4522343094041945E+01;
  COFD[         654] =  -0.3705240910138461E+00;
  COFD[         655] =   0.1597653817363581E-01;
  COFD[         656] =  -0.1966898676800112E+02;
  COFD[         657] =   0.4852357825322780E+01;
  COFD[         658] =  -0.3976302929205434E+00;
  COFD[         659] =   0.1651454081234691E-01;
  COFD[         660] =  -0.1805980237913599E+02;
  COFD[         661] =   0.4360474000056099E+01;
  COFD[         662] =  -0.3516356776422001E+00;
  COFD[         663] =   0.1524545596600636E-01;
  COFD[         664] =  -0.1878468158135861E+02;
  COFD[         665] =   0.4605034135660564E+01;
  COFD[         666] =  -0.3792399254380860E+00;
  COFD[         667] =   0.1627418668274144E-01;
  COFD[         668] =  -0.1869955819260176E+02;
  COFD[         669] =   0.4551222122230642E+01;
  COFD[         670] =  -0.3741499971521470E+00;
  COFD[         671] =   0.1613010818331805E-01;
  COFD[         672] =  -0.1934447945708759E+02;
  COFD[         673] =   0.4835867901258055E+01;
  COFD[         674] =  -0.3919040198993948E+00;
  COFD[         675] =   0.1611935633959806E-01;
  COFD[         676] =  -0.1499765769172950E+02;
  COFD[         677] =   0.3419142531858579E+01;
  COFD[         678] =  -0.2369480887688776E+00;
  COFD[         679] =   0.1059345659536257E-01;
  COFD[         680] =  -0.1558361674716745E+02;
  COFD[         681] =   0.3987098568284434E+01;
  COFD[         682] =  -0.3109866036899985E+00;
  COFD[         683] =   0.1380795341338578E-01;
  COFD[         684] =  -0.1664237108262573E+02;
  COFD[         685] =   0.4321207774376617E+01;
  COFD[         686] =  -0.3181315549506706E+00;
  COFD[         687] =   0.1258308183870717E-01;
  COFD[         688] =  -0.1871735010457793E+02;
  COFD[         689] =   0.4573421885700214E+01;
  COFD[         690] =  -0.3757214086933051E+00;
  COFD[         691] =   0.1614451502251515E-01;
  COFD[         692] =  -0.1768659705517047E+02;
  COFD[         693] =   0.4367905899188474E+01;
  COFD[         694] =  -0.3516645879639785E+00;
  COFD[         695] =   0.1520186726409712E-01;
  COFD[         696] =  -0.1719414340129819E+02;
  COFD[         697] =   0.3490298819499314E+01;
  COFD[         698] =  -0.1625329575445097E+00;
  COFD[         699] =   0.4043757780510787E-02;
  COFD[         700] =  -0.1872670165824808E+02;
  COFD[         701] =   0.4573903479896379E+01;
  COFD[         702] =  -0.3758340786206700E+00;
  COFD[         703] =   0.1615161324401359E-01;
  COFD[         704] =  -0.1873560076291526E+02;
  COFD[         705] =   0.4574325021189105E+01;
  COFD[         706] =  -0.3759338652876284E+00;
  COFD[         707] =   0.1615792922437573E-01;
  COFD[         708] =  -0.1913993542632300E+02;
  COFD[         709] =   0.4730620081114692E+01;
  COFD[         710] =  -0.3883890489649657E+00;
  COFD[         711] =   0.1637150502726411E-01;
  COFD[         712] =  -0.1913011857829392E+02;
  COFD[         713] =   0.4728539523227080E+01;
  COFD[         714] =  -0.3887605628634140E+00;
  COFD[         715] =   0.1641451966125198E-01;
  COFD[         716] =  -0.1856069485825941E+02;
  COFD[         717] =   0.4522343094041945E+01;
  COFD[         718] =  -0.3705240910138461E+00;
  COFD[         719] =   0.1597653817363581E-01;
  COFD[         720] =  -0.2016579860319706E+02;
  COFD[         721] =   0.4877618917687775E+01;
  COFD[         722] =  -0.3948460761426021E+00;
  COFD[         723] =   0.1615166256107893E-01;
  COFD[         724] =  -0.1928397124960116E+02;
  COFD[         725] =   0.4312110001899904E+01;
  COFD[         726] =  -0.2907322088872635E+00;
  COFD[         727] =   0.1043236734543387E-01;
  COFD[         728] =  -0.1989540745008587E+02;
  COFD[         729] =   0.4834853238941779E+01;
  COFD[         730] =  -0.3932322467650001E+00;
  COFD[         731] =   0.1623552876595314E-01;
  COFD[         732] =  -0.2005638433296784E+02;
  COFD[         733] =   0.4810333449107644E+01;
  COFD[         734] =  -0.3802752927546432E+00;
  COFD[         735] =   0.1528555688791561E-01;
  COFD[         736] =  -0.2016849050066446E+02;
  COFD[         737] =   0.4851147993946014E+01;
  COFD[         738] =  -0.3897545968990269E+00;
  COFD[         739] =   0.1586327845280228E-01;
  COFD[         740] =  -0.1833717188007337E+02;
  COFD[         741] =   0.4039707350126803E+01;
  COFD[         742] =  -0.2488565303881858E+00;
  COFD[         743] =   0.8361653792794556E-02;
  COFD[         744] =  -0.1761289252565630E+02;
  COFD[         745] =   0.4327706010171974E+01;
  COFD[         746] =  -0.3490300604788930E+00;
  COFD[         747] =   0.1519713199892770E-01;
  COFD[         748] =  -0.1688564690991911E+02;
  COFD[         749] =   0.4415865473479337E+01;
  COFD[         750] =  -0.3564104485793251E+00;
  COFD[         751] =   0.1533855847922275E-01;
  COFD[         752] =  -0.1567567878102999E+02;
  COFD[         753] =   0.3706302510623138E+01;
  COFD[         754] =  -0.2096236244277931E+00;
  COFD[         755] =   0.6748353057110522E-02;
  COFD[         756] =  -0.1977364410044468E+02;
  COFD[         757] =   0.4876531592391508E+01;
  COFD[         758] =  -0.3981741984755758E+00;
  COFD[         759] =   0.1643476290195535E-01;
  COFD[         760] =  -0.1891518602118661E+02;
  COFD[         761] =   0.4756739681751712E+01;
  COFD[         762] =  -0.3886685070061102E+00;
  COFD[         763] =   0.1625748441332711E-01;
  COFD[         764] =  -0.1468058225757294E+02;
  COFD[         765] =   0.2298576535112391E+01;
  COFD[         766] =   0.1810975986268222E-01;
  COFD[         767] =  -0.4731609614535192E-02;
  COFD[         768] =  -0.1978049689264223E+02;
  COFD[         769] =   0.4876281709421428E+01;
  COFD[         770] =  -0.3981401304053240E+00;
  COFD[         771] =   0.1643325343977014E-01;
  COFD[         772] =  -0.1978644159772173E+02;
  COFD[         773] =   0.4875758120324897E+01;
  COFD[         774] =  -0.3980631994813040E+00;
  COFD[         775] =   0.1642957450530744E-01;
  COFD[         776] =  -0.1958734834972716E+02;
  COFD[         777] =   0.4743227463757541E+01;
  COFD[         778] =  -0.3686474607497347E+00;
  COFD[         779] =   0.1466232449273593E-01;
  COFD[         780] =  -0.1962180466510501E+02;
  COFD[         781] =   0.4762356031601866E+01;
  COFD[         782] =  -0.3720986023672735E+00;
  COFD[         783] =   0.1485121157223651E-01;
  COFD[         784] =  -0.1966898676800112E+02;
  COFD[         785] =   0.4852357825322780E+01;
  COFD[         786] =  -0.3976302929205434E+00;
  COFD[         787] =   0.1651454081234691E-01;
  COFD[         788] =  -0.1928397124960116E+02;
  COFD[         789] =   0.4312110001899904E+01;
  COFD[         790] =  -0.2907322088872635E+00;
  COFD[         791] =   0.1043236734543387E-01;
  COFD[         792] =  -0.1603997666374001E+02;
  COFD[         793] =   0.2734208547161103E+01;
  COFD[         794] =  -0.4628114204507498E-01;
  COFD[         795] =  -0.1672840630334733E-02;
  COFD[         796] =  -0.1961392673202087E+02;
  COFD[         797] =   0.4524886156774725E+01;
  COFD[         798] =  -0.3261423023221592E+00;
  COFD[         799] =   0.1226835649915157E-01;
  COFD[         800] =  -0.1890313414494545E+02;
  COFD[         801] =   0.4121436960338011E+01;
  COFD[         802] =  -0.2594362965642011E+00;
  COFD[         803] =   0.8823603190561257E-02;
  COFD[         804] =  -0.1932385614900669E+02;
  COFD[         805] =   0.4295315571390907E+01;
  COFD[         806] =  -0.2874898966822066E+00;
  COFD[         807] =   0.1025198769300411E-01;
  COFD[         808] =  -0.1580770444028599E+02;
  COFD[         809] =   0.2799524263213967E+01;
  COFD[         810] =  -0.5666575630034967E-01;
  COFD[         811] =  -0.1137221314445150E-02;
  COFD[         812] =  -0.1888662285046225E+02;
  COFD[         813] =   0.4736242637853255E+01;
  COFD[         814] =  -0.3893357605389257E+00;
  COFD[         815] =   0.1642422988902124E-01;
  COFD[         816] =  -0.1548207029093904E+02;
  COFD[         817] =   0.3939084188623417E+01;
  COFD[         818] =  -0.3072761108436687E+00;
  COFD[         819] =   0.1375770270895890E-01;
  COFD[         820] =  -0.1685361380845761E+02;
  COFD[         821] =   0.4411251776453823E+01;
  COFD[         822] =  -0.3350756131548513E+00;
  COFD[         823] =   0.1353157937321782E-01;
  COFD[         824] =  -0.1830650296367939E+02;
  COFD[         825] =   0.4452394507404697E+01;
  COFD[         826] =  -0.3624991200767191E+00;
  COFD[         827] =   0.1567286266089603E-01;
  COFD[         828] =  -0.1730907989620245E+02;
  COFD[         829] =   0.4243749180870270E+01;
  COFD[         830] =  -0.3387466685616417E+00;
  COFD[         831] =   0.1477905643158847E-01;
  COFD[         832] =  -0.1868831150678477E+02;
  COFD[         833] =   0.4247585976166683E+01;
  COFD[         834] =  -0.2815735628073077E+00;
  COFD[         835] =   0.1000602366593966E-01;
  COFD[         836] =  -0.1831322697803230E+02;
  COFD[         837] =   0.4452261161272578E+01;
  COFD[         838] =  -0.3624715246876605E+00;
  COFD[         839] =   0.1567121551911240E-01;
  COFD[         840] =  -0.1831968967523493E+02;
  COFD[         841] =   0.4452134774757628E+01;
  COFD[         842] =  -0.3624420830935042E+00;
  COFD[         843] =   0.1566936578675393E-01;
  COFD[         844] =  -0.1892132288088230E+02;
  COFD[         845] =   0.4686257416133524E+01;
  COFD[         846] =  -0.3876178620425971E+00;
  COFD[         847] =   0.1654627973246352E-01;
  COFD[         848] =  -0.1887963312848053E+02;
  COFD[         849] =   0.4671422222663214E+01;
  COFD[         850] =  -0.3861362810478654E+00;
  COFD[         851] =   0.1650068107616750E-01;
  COFD[         852] =  -0.1805980237913599E+02;
  COFD[         853] =   0.4360474000056099E+01;
  COFD[         854] =  -0.3516356776422001E+00;
  COFD[         855] =   0.1524545596600636E-01;
  COFD[         856] =  -0.1989540745008587E+02;
  COFD[         857] =   0.4834853238941779E+01;
  COFD[         858] =  -0.3932322467650001E+00;
  COFD[         859] =   0.1623552876595314E-01;
  COFD[         860] =  -0.1961392673202087E+02;
  COFD[         861] =   0.4524886156774725E+01;
  COFD[         862] =  -0.3261423023221592E+00;
  COFD[         863] =   0.1226835649915157E-01;
  COFD[         864] =  -0.1977017209130380E+02;
  COFD[         865] =   0.4847119688172202E+01;
  COFD[         866] =  -0.4003461404722360E+00;
  COFD[         867] =   0.1676917774316241E-01;
  COFD[         868] =  -0.2005922850303636E+02;
  COFD[         869] =   0.4881576631488546E+01;
  COFD[         870] =  -0.3958945479488925E+00;
  COFD[         871] =   0.1621878544379618E-01;
  COFD[         872] =  -0.2004561228963652E+02;
  COFD[         873] =   0.4868786155445560E+01;
  COFD[         874] =  -0.3974311326078929E+00;
  COFD[         875] =   0.1641238256860971E-01;
  COFD[         876] =  -0.1912297200146966E+02;
  COFD[         877] =   0.4463559909783323E+01;
  COFD[         878] =  -0.3178210889061098E+00;
  COFD[         879] =   0.1189219041826154E-01;
  COFD[         880] =  -0.1703791763653080E+02;
  COFD[         881] =   0.4110925857281930E+01;
  COFD[         882] =  -0.3221728744644342E+00;
  COFD[         883] =   0.1408821528952375E-01;
  COFD[         884] =  -0.1606152931877502E+02;
  COFD[         885] =   0.4146430359164761E+01;
  COFD[         886] =  -0.3312210923439469E+00;
  COFD[         887] =   0.1466430884967795E-01;
  COFD[         888] =  -0.1635823810097287E+02;
  COFD[         889] =   0.4123977007004703E+01;
  COFD[         890] =  -0.2854728908229589E+00;
  COFD[         891] =   0.1088433914065527E-01;
  COFD[         892] =  -0.1904759289384530E+02;
  COFD[         893] =   0.4700929573683383E+01;
  COFD[         894] =  -0.3903608191209737E+00;
  COFD[         895] =   0.1670229359855150E-01;
  COFD[         896] =  -0.1795985745949001E+02;
  COFD[         897] =   0.4461091096847464E+01;
  COFD[         898] =  -0.3624212969270659E+00;
  COFD[         899] =   0.1561553586311572E-01;
  COFD[         900] =  -0.1767910948230725E+02;
  COFD[         901] =   0.3722086580781685E+01;
  COFD[         902] =  -0.1980254146986777E+00;
  COFD[         903] =   0.5799038224028695E-02;
  COFD[         904] =  -0.1905419842832635E+02;
  COFD[         905] =   0.4700662819460902E+01;
  COFD[         906] =  -0.3903206776226419E+00;
  COFD[         907] =   0.1670032573134212E-01;
  COFD[         908] =  -0.1906036212420135E+02;
  COFD[         909] =   0.4700323443044125E+01;
  COFD[         910] =  -0.3902665913786998E+00;
  COFD[         911] =   0.1669755465043006E-01;
  COFD[         912] =  -0.1934401530563503E+02;
  COFD[         913] =   0.4788242539311199E+01;
  COFD[         914] =  -0.3933254249328604E+00;
  COFD[         915] =   0.1648505821153991E-01;
  COFD[         916] =  -0.1933004429336219E+02;
  COFD[         917] =   0.4785693509382203E+01;
  COFD[         918] =  -0.3936420075604487E+00;
  COFD[         919] =   0.1652573779841039E-01;
  COFD[         920] =  -0.1878468158135861E+02;
  COFD[         921] =   0.4605034135660564E+01;
  COFD[         922] =  -0.3792399254380860E+00;
  COFD[         923] =   0.1627418668274144E-01;
  COFD[         924] =  -0.2005638433296784E+02;
  COFD[         925] =   0.4810333449107644E+01;
  COFD[         926] =  -0.3802752927546432E+00;
  COFD[         927] =   0.1528555688791561E-01;
  COFD[         928] =  -0.1890313414494545E+02;
  COFD[         929] =   0.4121436960338011E+01;
  COFD[         930] =  -0.2594362965642011E+00;
  COFD[         931] =   0.8823603190561257E-02;
  COFD[         932] =  -0.2005922850303636E+02;
  COFD[         933] =   0.4881576631488546E+01;
  COFD[         934] =  -0.3958945479488925E+00;
  COFD[         935] =   0.1621878544379618E-01;
  COFD[         936] =  -0.1999234499569811E+02;
  COFD[         937] =   0.4757309526477115E+01;
  COFD[         938] =  -0.3684627400545891E+00;
  COFD[         939] =   0.1457942336558470E-01;
  COFD[         940] =  -0.2015885283044970E+02;
  COFD[         941] =   0.4822728985557825E+01;
  COFD[         942] =  -0.3812000138893480E+00;
  COFD[         943] =   0.1529926403243838E-01;
  COFD[         944] =  -0.1846533503760062E+02;
  COFD[         945] =   0.4088215083590887E+01;
  COFD[         946] =  -0.2553371937407484E+00;
  COFD[         947] =   0.8657586372083866E-02;
  COFD[         948] =  -0.1776457707271082E+02;
  COFD[         949] =   0.4365715498242292E+01;
  COFD[         950] =  -0.3517547303197583E+00;
  COFD[         951] =   0.1522394939239619E-01;
  COFD[         952] =  -0.1596236229610858E+02;
  COFD[         953] =   0.4079842750882008E+01;
  COFD[         954] =  -0.3231960602071978E+00;
  COFD[         955] =   0.1434305777680907E-01;
  COFD[         956] =  -0.1646195534760534E+02;
  COFD[         957] =   0.4141822899037704E+01;
  COFD[         958] =  -0.2892308013569180E+00;
  COFD[         959] =   0.1109614462742981E-01;
  COFD[         960] =  -0.1883799365221673E+02;
  COFD[         961] =   0.4594611640098751E+01;
  COFD[         962] =  -0.3781168670160998E+00;
  COFD[         963] =   0.1623467419613562E-01;
  COFD[         964] =  -0.1782798069130013E+02;
  COFD[         965] =   0.4381516462622495E+01;
  COFD[         966] =  -0.3532970918753254E+00;
  COFD[         967] =   0.1526804966457564E-01;
  COFD[         968] =  -0.1813805288510309E+02;
  COFD[         969] =   0.3907124675584541E+01;
  COFD[         970] =  -0.2273991439446390E+00;
  COFD[         971] =   0.7280601459305962E-02;
  COFD[         972] =  -0.1884529530697174E+02;
  COFD[         973] =   0.4594634055023398E+01;
  COFD[         974] =  -0.3781333057516071E+00;
  COFD[         975] =   0.1623599168682821E-01;
  COFD[         976] =  -0.1885218260686505E+02;
  COFD[         977] =   0.4594596529989258E+01;
  COFD[         978] =  -0.3781371531513473E+00;
  COFD[         979] =   0.1623655159949430E-01;
  COFD[         980] =  -0.1926514266945471E+02;
  COFD[         981] =   0.4743728688652060E+01;
  COFD[         982] =  -0.3899270700585402E+00;
  COFD[         983] =   0.1643466452124135E-01;
  COFD[         984] =  -0.1926485085357470E+02;
  COFD[         985] =   0.4746666060969545E+01;
  COFD[         986] =  -0.3910387566978661E+00;
  COFD[         987] =   0.1651352760050653E-01;
  COFD[         988] =  -0.1869955819260176E+02;
  COFD[         989] =   0.4551222122230642E+01;
  COFD[         990] =  -0.3741499971521470E+00;
  COFD[         991] =   0.1613010818331805E-01;
  COFD[         992] =  -0.2016849050066446E+02;
  COFD[         993] =   0.4851147993946014E+01;
  COFD[         994] =  -0.3897545968990269E+00;
  COFD[         995] =   0.1586327845280228E-01;
  COFD[         996] =  -0.1932385614900669E+02;
  COFD[         997] =   0.4295315571390907E+01;
  COFD[         998] =  -0.2874898966822066E+00;
  COFD[         999] =   0.1025198769300411E-01;
  COFD[        1000] =  -0.2004561228963652E+02;
  COFD[        1001] =   0.4868786155445560E+01;
  COFD[        1002] =  -0.3974311326078929E+00;
  COFD[        1003] =   0.1641238256860971E-01;
  COFD[        1004] =  -0.2015885283044970E+02;
  COFD[        1005] =   0.4822728985557825E+01;
  COFD[        1006] =  -0.3812000138893480E+00;
  COFD[        1007] =   0.1529926403243838E-01;
  COFD[        1008] =  -0.2024931478123324E+02;
  COFD[        1009] =   0.4856999247129295E+01;
  COFD[        1010] =  -0.3896182557523967E+00;
  COFD[        1011] =   0.1582188160996442E-01;
  COFD[        1012] =  -0.1871694128239326E+02;
  COFD[        1013] =   0.4183824876778368E+01;
  COFD[        1014] =  -0.2717303313431965E+00;
  COFD[        1015] =   0.9517448825356330E-02;
  COFD[        1016] =  -0.1776451185181849E+02;
  COFD[        1017] =   0.4347765208885912E+01;
  COFD[        1018] =  -0.3515930667599710E+00;
  COFD[        1019] =   0.1530764063925441E-01;
  COFD[        1020] =  -0.1682806034335371E+02;
  COFD[        1021] =   0.4485450500929320E+01;
  COFD[        1022] =  -0.3655695289846634E+00;
  COFD[        1023] =   0.1574614737371121E-01;
  COFD[        1024] =  -0.1626359848098758E+02;
  COFD[        1025] =   0.4116788683843708E+01;
  COFD[        1026] =  -0.2733289330957227E+00;
  COFD[        1027] =   0.9943035956817482E-02;
  COFD[        1028] =  -0.1934710272078301E+02;
  COFD[        1029] =   0.4822779400529180E+01;
  COFD[        1030] =  -0.3871228449439424E+00;
  COFD[        1031] =   0.1578686341252481E-01;
  COFD[        1032] =  -0.1863707722655403E+02;
  COFD[        1033] =   0.4806992114833483E+01;
  COFD[        1034] =  -0.3984196894018535E+00;
  COFD[        1035] =   0.1681856268295304E-01;
  COFD[        1036] =  -0.1316791411191276E+02;
  COFD[        1037] =   0.1635596276910710E+01;
  COFD[        1038] =   0.1248562922809231E+00;
  COFD[        1039] =  -0.1002143991895203E-01;
  COFD[        1040] =  -0.1919185322703664E+02;
  COFD[        1041] =   0.4784033288807873E+01;
  COFD[        1042] =  -0.3849459064100618E+00;
  COFD[        1043] =   0.1580648990064550E-01;
  COFD[        1044] =  -0.1918517200292939E+02;
  COFD[        1045] =   0.4778400021399422E+01;
  COFD[        1046] =  -0.3840658172050402E+00;
  COFD[        1047] =   0.1576203256152660E-01;
  COFD[        1048] =  -0.1944043258564383E+02;
  COFD[        1049] =   0.4828743160693618E+01;
  COFD[        1050] =  -0.3828021863328961E+00;
  COFD[        1051] =   0.1540186987470838E-01;
  COFD[        1052] =  -0.1941730561245821E+02;
  COFD[        1053] =   0.4788578363438019E+01;
  COFD[        1054] =  -0.3737994268951498E+00;
  COFD[        1055] =   0.1486050997611076E-01;
  COFD[        1056] =  -0.1934447945708759E+02;
  COFD[        1057] =   0.4835867901258055E+01;
  COFD[        1058] =  -0.3919040198993948E+00;
  COFD[        1059] =   0.1611935633959806E-01;
  COFD[        1060] =  -0.1833717188007337E+02;
  COFD[        1061] =   0.4039707350126803E+01;
  COFD[        1062] =  -0.2488565303881858E+00;
  COFD[        1063] =   0.8361653792794556E-02;
  COFD[        1064] =  -0.1580770444028599E+02;
  COFD[        1065] =   0.2799524263213967E+01;
  COFD[        1066] =  -0.5666575630034967E-01;
  COFD[        1067] =  -0.1137221314445150E-02;
  COFD[        1068] =  -0.1912297200146966E+02;
  COFD[        1069] =   0.4463559909783323E+01;
  COFD[        1070] =  -0.3178210889061098E+00;
  COFD[        1071] =   0.1189219041826154E-01;
  COFD[        1072] =  -0.1846533503760062E+02;
  COFD[        1073] =   0.4088215083590887E+01;
  COFD[        1074] =  -0.2553371937407484E+00;
  COFD[        1075] =   0.8657586372083866E-02;
  COFD[        1076] =  -0.1871694128239326E+02;
  COFD[        1077] =   0.4183824876778368E+01;
  COFD[        1078] =  -0.2717303313431965E+00;
  COFD[        1079] =   0.9517448825356330E-02;
  COFD[        1080] =  -0.1451532645593394E+02;
  COFD[        1081] =   0.2310965406901932E+01;
  COFD[        1082] =   0.1872493495986990E-01;
  COFD[        1083] =  -0.4810190587113683E-02;
  COFD[        1084] =  -0.1862683104563700E+02;
  COFD[        1085] =   0.4778876685436832E+01;
  COFD[        1086] =  -0.3977561689282376E+00;
  COFD[        1087] =   0.1691190324024170E-01;
  COFD[        1088] =  -0.1163276178298392E+02;
  COFD[        1089] =   0.2530205071288130E+01;
  COFD[        1090] =  -0.1183776350934922E+00;
  COFD[        1091] =   0.5312988352557955E-02;
  COFD[        1092] =  -0.1522605770447405E+02;
  COFD[        1093] =   0.4154863253653862E+01;
  COFD[        1094] =  -0.3321404775716639E+00;
  COFD[        1095] =   0.1469592836326692E-01;
  COFD[        1096] =  -0.1527049784849730E+02;
  COFD[        1097] =   0.3538365317994180E+01;
  COFD[        1098] =  -0.2526814234121028E+00;
  COFD[        1099] =   0.1128573168209360E-01;
  COFD[        1100] =  -0.1395220381010925E+02;
  COFD[        1101] =   0.3159688299904867E+01;
  COFD[        1102] =  -0.2028455188834762E+00;
  COFD[        1103] =   0.9099190221367019E-02;
  COFD[        1104] =  -0.1875637020961290E+02;
  COFD[        1105] =   0.4821872691746870E+01;
  COFD[        1106] =  -0.3987313751514457E+00;
  COFD[        1107] =   0.1676685214877852E-01;
  COFD[        1108] =  -0.1528586419146581E+02;
  COFD[        1109] =   0.3542872339374690E+01;
  COFD[        1110] =  -0.2532972672017028E+00;
  COFD[        1111] =   0.1131367837559019E-01;
  COFD[        1112] =  -0.1530075071186979E+02;
  COFD[        1113] =   0.3547276678457322E+01;
  COFD[        1114] =  -0.2538990508478524E+00;
  COFD[        1115] =   0.1134098602286849E-01;
  COFD[        1116] =  -0.1573476318984406E+02;
  COFD[        1117] =   0.3712839423735899E+01;
  COFD[        1118] =  -0.2734220687823376E+00;
  COFD[        1119] =   0.1210212757075583E-01;
  COFD[        1120] =  -0.1569001110440693E+02;
  COFD[        1121] =   0.3696584277188425E+01;
  COFD[        1122] =  -0.2713867054547939E+00;
  COFD[        1123] =   0.1201712146112774E-01;
  COFD[        1124] =  -0.1499765769172950E+02;
  COFD[        1125] =   0.3419142531858579E+01;
  COFD[        1126] =  -0.2369480887688776E+00;
  COFD[        1127] =   0.1059345659536257E-01;
  COFD[        1128] =  -0.1761289252565630E+02;
  COFD[        1129] =   0.4327706010171974E+01;
  COFD[        1130] =  -0.3490300604788930E+00;
  COFD[        1131] =   0.1519713199892770E-01;
  COFD[        1132] =  -0.1888662285046225E+02;
  COFD[        1133] =   0.4736242637853255E+01;
  COFD[        1134] =  -0.3893357605389257E+00;
  COFD[        1135] =   0.1642422988902124E-01;
  COFD[        1136] =  -0.1703791763653080E+02;
  COFD[        1137] =   0.4110925857281930E+01;
  COFD[        1138] =  -0.3221728744644342E+00;
  COFD[        1139] =   0.1408821528952375E-01;
  COFD[        1140] =  -0.1776457707271082E+02;
  COFD[        1141] =   0.4365715498242292E+01;
  COFD[        1142] =  -0.3517547303197583E+00;
  COFD[        1143] =   0.1522394939239619E-01;
  COFD[        1144] =  -0.1776451185181849E+02;
  COFD[        1145] =   0.4347765208885912E+01;
  COFD[        1146] =  -0.3515930667599710E+00;
  COFD[        1147] =   0.1530764063925441E-01;
  COFD[        1148] =  -0.1862683104563700E+02;
  COFD[        1149] =   0.4778876685436832E+01;
  COFD[        1150] =  -0.3977561689282376E+00;
  COFD[        1151] =   0.1691190324024170E-01;
  COFD[        1152] =  -0.1377266831742670E+02;
  COFD[        1153] =   0.3040578783767133E+01;
  COFD[        1154] =  -0.1869754178307637E+00;
  COFD[        1155] =   0.8394057441551481E-02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKTDIF EGTRANSETKTDIF
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKTDIF egtransetktdif
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKTDIF egtransetktdif_
#endif
void egtransetKTDIF(int* KTDIF) {
  KTDIF[           0] =            1;
  KTDIF[           1] =            2;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFTD EGTRANSETCOFTD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFTD egtransetcoftd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFTD egtransetcoftd_
#endif
void egtransetCOFTD(double* COFTD) {
  COFTD[           0] =   0.0000000000000000E+00;
  COFTD[           1] =   0.0000000000000000E+00;
  COFTD[           2] =   0.0000000000000000E+00;
  COFTD[           3] =   0.0000000000000000E+00;
  COFTD[           4] =  -0.1267157358260665E+00;
  COFTD[           5] =  -0.1025304929454753E-03;
  COFTD[           6] =   0.5456049480958656E-07;
  COFTD[           7] =  -0.8851811492383965E-11;
  COFTD[           8] =   0.3818642928024142E+00;
  COFTD[           9] =   0.1841173676342166E-03;
  COFTD[          10] =  -0.9796176034833310E-07;
  COFTD[          11] =   0.1625422476308086E-10;
  COFTD[          12] =   0.3754754585061315E+00;
  COFTD[          13] =   0.9735779997691615E-04;
  COFTD[          14] =  -0.4943935849608963E-07;
  COFTD[          15] =   0.8333052020411856E-11;
  COFTD[          16] =   0.1951583360602999E-02;
  COFTD[          17] =   0.6694708047799631E-03;
  COFTD[          18] =  -0.3121487146752215E-06;
  COFTD[          19] =   0.4529499206739157E-10;
  COFTD[          20] =   0.3833421804573635E+00;
  COFTD[          21] =   0.1848299369679226E-03;
  COFTD[          22] =  -0.9834089104739220E-07;
  COFTD[          23] =   0.1631713171345801E-10;
  COFTD[          24] =   0.3847373793384187E+00;
  COFTD[          25] =   0.1855026375847330E-03;
  COFTD[          26] =  -0.9869880913766751E-07;
  COFTD[          27] =   0.1637651897911576E-10;
  COFTD[          28] =   0.2912627324548178E+00;
  COFTD[          29] =   0.2330450691123891E-03;
  COFTD[          30] =  -0.1240401108265922E-06;
  COFTD[          31] =   0.2013458793291140E-10;
  COFTD[          32] =   0.2989740640904121E+00;
  COFTD[          33] =   0.2322310123534718E-03;
  COFTD[          34] =  -0.1236750374915925E-06;
  COFTD[          35] =   0.2010297352068701E-10;
  COFTD[          36] =   0.3874054384854359E+00;
  COFTD[          37] =   0.1568838091395251E-03;
  COFTD[          38] =  -0.8293099098207665E-07;
  COFTD[          39] =   0.1384603182952203E-10;
  COFTD[          40] =   0.2471291249759435E+00;
  COFTD[          41] =   0.4493957116165469E-03;
  COFTD[          42] =  -0.2320307603798581E-06;
  COFTD[          43] =   0.3625788270902565E-10;
  COFTD[          44] =   0.8921833323260459E-01;
  COFTD[          45] =   0.6385973457264696E-03;
  COFTD[          46] =  -0.3105267185106736E-06;
  COFTD[          47] =   0.4632182714912095E-10;
  COFTD[          48] =   0.2613831106363483E+00;
  COFTD[          49] =   0.3741001421151409E-03;
  COFTD[          50] =  -0.1952756438202893E-06;
  COFTD[          51] =   0.3084446425883240E-10;
  COFTD[          52] =   0.2066211976818736E+00;
  COFTD[          53] =   0.4698391340071242E-03;
  COFTD[          54] =  -0.2399965709388047E-06;
  COFTD[          55] =   0.3714866667604425E-10;
  COFTD[          56] =   0.2301782357423144E+00;
  COFTD[          57] =   0.4411239723719602E-03;
  COFTD[          58] =  -0.2271927867378409E-06;
  COFTD[          59] =   0.3542100728590062E-10;
  COFTD[          60] =   0.5978875014515522E-01;
  COFTD[          61] =   0.6004068868355234E-03;
  COFTD[          62] =  -0.2889528616835724E-06;
  COFTD[          63] =   0.4280270387409600E-10;
  COFTD[          64] =   0.3671545750621259E+00;
  COFTD[          65] =   0.7074871807561138E-04;
  COFTD[          66] =  -0.3401390884461344E-07;
  COFTD[          67] =   0.5718817519581169E-11;
  COFTD[          68] =   0.1267157358260665E+00;
  COFTD[          69] =   0.1025304929454753E-03;
  COFTD[          70] =  -0.5456049480958656E-07;
  COFTD[          71] =   0.8851811492383965E-11;
  COFTD[          72] =   0.0000000000000000E+00;
  COFTD[          73] =   0.0000000000000000E+00;
  COFTD[          74] =   0.0000000000000000E+00;
  COFTD[          75] =   0.0000000000000000E+00;
  COFTD[          76] =   0.1397457852745939E+00;
  COFTD[          77] =   0.6298108591719154E-03;
  COFTD[          78] =  -0.3116940342118423E-06;
  COFTD[          79] =   0.4707558637232602E-10;
  COFTD[          80] =   0.1945605631423317E+00;
  COFTD[          81] =   0.5079167441644920E-03;
  COFTD[          82] =  -0.2577206587518775E-06;
  COFTD[          83] =   0.3967209139227169E-10;
  COFTD[          84] =  -0.1609284370381552E+00;
  COFTD[          85] =   0.8016856147508327E-03;
  COFTD[          86] =  -0.3249766388129991E-06;
  COFTD[          87] =   0.4319581912733095E-10;
  COFTD[          88] =   0.1400151641499396E+00;
  COFTD[          89] =   0.6310249046659604E-03;
  COFTD[          90] =  -0.3122948665605463E-06;
  COFTD[          91] =   0.4716633092314242E-10;
  COFTD[          92] =   0.1402690373253105E+00;
  COFTD[          93] =   0.6321690685660495E-03;
  COFTD[          94] =  -0.3128611144373887E-06;
  COFTD[          95] =   0.4725185213275378E-10;
  COFTD[          96] =   0.6875565922217064E-01;
  COFTD[          97] =   0.6631049392258758E-03;
  COFTD[          98] =  -0.3194849555055448E-06;
  COFTD[          99] =   0.4736101303399774E-10;
  COFTD[         100] =   0.7315173371679451E-01;
  COFTD[         101] =   0.6642955892285024E-03;
  COFTD[         102] =  -0.3206112835377731E-06;
  COFTD[         103] =   0.4758315870108564E-10;
  COFTD[         104] =   0.1587269043383360E+00;
  COFTD[         105] =   0.5967534238686169E-03;
  COFTD[         106] =  -0.2976737939540804E-06;
  COFTD[         107] =   0.4521938974583591E-10;
  COFTD[         108] =  -0.3849651678153934E-01;
  COFTD[         109] =   0.8347949814005065E-03;
  COFTD[         110] =  -0.3810312564452182E-06;
  COFTD[         111] =   0.5455319374182929E-10;
  COFTD[         112] =  -0.1529552889766811E+00;
  COFTD[         113] =   0.8501473210897147E-03;
  COFTD[         114] =  -0.3529415371202759E-06;
  COFTD[         115] =   0.4763453033863902E-10;
  COFTD[         116] =  -0.6479936638152356E-02;
  COFTD[         117] =   0.7836066032604117E-03;
  COFTD[         118] =  -0.3637347512791469E-06;
  COFTD[         119] =   0.5263072586686018E-10;
  COFTD[         120] =  -0.6410089132874891E-01;
  COFTD[         121] =   0.8311312698059269E-03;
  COFTD[         122] =  -0.3732877529348841E-06;
  COFTD[         123] =   0.5291055177769992E-10;
  COFTD[         124] =  -0.4419302611803710E-01;
  COFTD[         125] =   0.8219675155377782E-03;
  COFTD[         126] =  -0.3737680415530062E-06;
  COFTD[         127] =   0.5338781013168227E-10;
  COFTD[         128] =  -0.1412627145163882E+00;
  COFTD[         129] =   0.8094401725653276E-03;
  COFTD[         130] =  -0.3379215246058519E-06;
  COFTD[         131] =   0.4576856099561897E-10;
  COFTD[         132] =   0.2126392909101736E+00;
  COFTD[         133] =   0.4604824955661828E-03;
  COFTD[         134] =  -0.2357734960762928E-06;
  COFTD[         135] =   0.3656872132522026E-10;
};




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
ELEMENTS
O   H   C   N
END
SPECIES
H2   H    O2   OH   H2O  HO2  H2O2 CH3  CH4  CO   CO2  CH2O C2H2 C2H4 C2H6 NH3 NO   HCN  N2
END
REACTIONS
END

\\
\\
\\  This is the therm file
\\
\\
THERMO
   300.000  1000.000  5000.000
! GRI-Mech Version 3.0 Thermodynamics released 7/30/99
! NASA Polynomial format for CHEMKIN-II
! see README file for disclaimer
O                 L 1/90O   1               G   200.000  3500.000  1000.000    1
 2.56942078E+00-8.59741137E-05 4.19484589E-08-1.00177799E-11 1.22833691E-15    2
 2.92175791E+04 4.78433864E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00                   4
O2                TPIS89O   2               G   200.000  3500.000  1000.000    1
 3.28253784E+00 1.48308754E-03-7.57966669E-07 2.09470555E-10-2.16717794E-14    2
-1.08845772E+03 5.45323129E+00 3.78245636E+00-2.99673416E-03 9.84730201E-06    3
-9.68129509E-09 3.24372837E-12-1.06394356E+03 3.65767573E+00                   4
H                 L 7/88H   1               G   200.000  3500.000  1000.000    1
 2.50000001E+00-2.30842973E-11 1.61561948E-14-4.73515235E-18 4.98197357E-22    2
 2.54736599E+04-4.46682914E-01 2.50000000E+00 7.05332819E-13-1.99591964E-15    3
 2.30081632E-18-9.27732332E-22 2.54736599E+04-4.46682853E-01                   4
H2                TPIS78H   2               G   200.000  3500.000  1000.000    1
 3.33727920E+00-4.94024731E-05 4.99456778E-07-1.79566394E-10 2.00255376E-14    2
-9.50158922E+02-3.20502331E+00 2.34433112E+00 7.98052075E-03-1.94781510E-05    3
 2.01572094E-08-7.37611761E-12-9.17935173E+02 6.83010238E-01                   4
OH                RUS 78O   1H   1          G   200.000  3500.000  1000.000    1
 3.09288767E+00 5.48429716E-04 1.26505228E-07-8.79461556E-11 1.17412376E-14    2
 3.85865700E+03 4.47669610E+00 3.99201543E+00-2.40131752E-03 4.61793841E-06    3
-3.88113333E-09 1.36411470E-12 3.61508056E+03-1.03925458E-01                   4
H2O               L 8/89H   2O   1          G   200.000  3500.000  1000.000    1
 3.03399249E+00 2.17691804E-03-1.64072518E-07-9.70419870E-11 1.68200992E-14    2
-3.00042971E+04 4.96677010E+00 4.19864056E+00-2.03643410E-03 6.52040211E-06    3
-5.48797062E-09 1.77197817E-12-3.02937267E+04-8.49032208E-01                   4
HO2               L 5/89H   1O   2          G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00                   4
H2O2              L 7/88H   2O   2          G   200.000  3500.000  1000.000    1
 4.16500285E+00 4.90831694E-03-1.90139225E-06 3.71185986E-10-2.87908305E-14    2
-1.78617877E+04 2.91615662E+00 4.27611269E+00-5.42822417E-04 1.67335701E-05    3
-2.15770813E-08 8.62454363E-12-1.77025821E+04 3.43505074E+00                   4
C                 L11/88C   1               G   200.000  3500.000  1000.000    1
 2.49266888E+00 4.79889284E-05-7.24335020E-08 3.74291029E-11-4.87277893E-15    2
 8.54512953E+04 4.80150373E+00 2.55423955E+00-3.21537724E-04 7.33792245E-07    3
-7.32234889E-10 2.66521446E-13 8.54438832E+04 4.53130848E+00                   4
CH                TPIS79C   1H   1          G   200.000  3500.000  1000.000    1
 2.87846473E+00 9.70913681E-04 1.44445655E-07-1.30687849E-10 1.76079383E-14    2
 7.10124364E+04 5.48497999E+00 3.48981665E+00 3.23835541E-04-1.68899065E-06    3
 3.16217327E-09-1.40609067E-12 7.07972934E+04 2.08401108E+00                   4
CH2               L S/93C   1H   2          G   200.000  3500.000  1000.000    1
 2.87410113E+00 3.65639292E-03-1.40894597E-06 2.60179549E-10-1.87727567E-14    2
 4.62636040E+04 6.17119324E+00 3.76267867E+00 9.68872143E-04 2.79489841E-06    3
-3.85091153E-09 1.68741719E-12 4.60040401E+04 1.56253185E+00                   4
CH2(S)            L S/93C   1H   2          G   200.000  3500.000  1000.000    1
 2.29203842E+00 4.65588637E-03-2.01191947E-06 4.17906000E-10-3.39716365E-14    2
 5.09259997E+04 8.62650169E+00 4.19860411E+00-2.36661419E-03 8.23296220E-06    3
-6.68815981E-09 1.94314737E-12 5.04968163E+04-7.69118967E-01                   4
CH3               L11/89C   1H   3          G   200.000  3500.000  1000.000    1
 2.28571772E+00 7.23990037E-03-2.98714348E-06 5.95684644E-10-4.67154394E-14    2
 1.67755843E+04 8.48007179E+00 3.67359040E+00 2.01095175E-03 5.73021856E-06    3
-6.87117425E-09 2.54385734E-12 1.64449988E+04 1.60456433E+00                   4
CH4               L 8/88C   1H   4          G   200.000  3500.000  1000.000    1
 7.48514950E-02 1.33909467E-02-5.73285809E-06 1.22292535E-09-1.01815230E-13    2
-9.46834459E+03 1.84373180E+01 5.14987613E+00-1.36709788E-02 4.91800599E-05    3
-4.84743026E-08 1.66693956E-11-1.02466476E+04-4.64130376E+00                   4
CO                TPIS79C   1O   1          G   200.000  3500.000  1000.000    1
 2.71518561E+00 2.06252743E-03-9.98825771E-07 2.30053008E-10-2.03647716E-14    2
-1.41518724E+04 7.81868772E+00 3.57953347E+00-6.10353680E-04 1.01681433E-06    3
 9.07005884E-10-9.04424499E-13-1.43440860E+04 3.50840928E+00                   4
CO2               L 7/88C   1O   2          G   200.000  3500.000  1000.000    1
 3.85746029E+00 4.41437026E-03-2.21481404E-06 5.23490188E-10-4.72084164E-14    2
-4.87591660E+04 2.27163806E+00 2.35677352E+00 8.98459677E-03-7.12356269E-06    3
 2.45919022E-09-1.43699548E-13-4.83719697E+04 9.90105222E+00                   4
HCO               L12/89H   1C   1O   1     G   200.000  3500.000  1000.000    1
 2.77217438E+00 4.95695526E-03-2.48445613E-06 5.89161778E-10-5.33508711E-14    2
 4.01191815E+03 9.79834492E+00 4.22118584E+00-3.24392532E-03 1.37799446E-05    3
-1.33144093E-08 4.33768865E-12 3.83956496E+03 3.39437243E+00                   4
CH2O              L 8/88H   2C   1O   1     G   200.000  3500.000  1000.000    1
 1.76069008E+00 9.20000082E-03-4.42258813E-06 1.00641212E-09-8.83855640E-14    2
-1.39958323E+04 1.36563230E+01 4.79372315E+00-9.90833369E-03 3.73220008E-05    3
-3.79285261E-08 1.31772652E-11-1.43089567E+04 6.02812900E-01                   4
CH2OH             GUNL93C   1H   3O   1     G   200.000  3500.000  1000.000    1
 3.69266569E+00 8.64576797E-03-3.75101120E-06 7.87234636E-10-6.48554201E-14    2
-3.24250627E+03 5.81043215E+00 3.86388918E+00 5.59672304E-03 5.93271791E-06    3
-1.04532012E-08 4.36967278E-12-3.19391367E+03 5.47302243E+00                   4
CH3O              121686C   1H   3O   1     G   300.00   3000.00   1000.000    1
 0.03770799E+02 0.07871497E-01-0.02656384E-04 0.03944431E-08-0.02112616E-12    2
 0.12783252E+03 0.02929575E+02 0.02106204E+02 0.07216595E-01 0.05338472E-04    3
-0.07377636E-07 0.02075610E-10 0.09786011E+04 0.13152177E+02                   4
CH3OH             L 8/88C   1H   4O   1     G   200.000  3500.000  1000.000    1
 1.78970791E+00 1.40938292E-02-6.36500835E-06 1.38171085E-09-1.17060220E-13    2
-2.53748747E+04 1.45023623E+01 5.71539582E+00-1.52309129E-02 6.52441155E-05    3
-7.10806889E-08 2.61352698E-11-2.56427656E+04-1.50409823E+00                   4
C2H               L 1/91C   2H   1          G   200.000  3500.000  1000.000    1
 3.16780652E+00 4.75221902E-03-1.83787077E-06 3.04190252E-10-1.77232770E-14    2
 6.71210650E+04 6.63589475E+00 2.88965733E+00 1.34099611E-02-2.84769501E-05    3
 2.94791045E-08-1.09331511E-11 6.68393932E+04 6.22296438E+00                   4
C2H2              L 1/91C   2H   2          G   200.000  3500.000  1000.000    1
 4.14756964E+00 5.96166664E-03-2.37294852E-06 4.67412171E-10-3.61235213E-14    2
 2.59359992E+04-1.23028121E+00 8.08681094E-01 2.33615629E-02-3.55171815E-05    3
 2.80152437E-08-8.50072974E-12 2.64289807E+04 1.39397051E+01                   4
C2H3              L 2/92C   2H   3          G   200.000  3500.000  1000.000    1
 3.01672400E+00 1.03302292E-02-4.68082349E-06 1.01763288E-09-8.62607041E-14    2
 3.46128739E+04 7.78732378E+00 3.21246645E+00 1.51479162E-03 2.59209412E-05    3
-3.57657847E-08 1.47150873E-11 3.48598468E+04 8.51054025E+00                   4
C2H4              L 1/91C   2H   4          G   200.000  3500.000  1000.000    1
 2.03611116E+00 1.46454151E-02-6.71077915E-06 1.47222923E-09-1.25706061E-13    2
 4.93988614E+03 1.03053693E+01 3.95920148E+00-7.57052247E-03 5.70990292E-05    3
-6.91588753E-08 2.69884373E-11 5.08977593E+03 4.09733096E+00                   4
C2H5              L12/92C   2H   5          G   200.000  3500.000  1000.000    1
 1.95465642E+00 1.73972722E-02-7.98206668E-06 1.75217689E-09-1.49641576E-13    2
 1.28575200E+04 1.34624343E+01 4.30646568E+00-4.18658892E-03 4.97142807E-05    3
-5.99126606E-08 2.30509004E-11 1.28416265E+04 4.70720924E+00                   4
C2H6              L 8/88C   2H   6          G   200.000  3500.000  1000.000    1
 1.07188150E+00 2.16852677E-02-1.00256067E-05 2.21412001E-09-1.90002890E-13    2
-1.14263932E+04 1.51156107E+01 4.29142492E+00-5.50154270E-03 5.99438288E-05    3
-7.08466285E-08 2.68685771E-11-1.15222055E+04 2.66682316E+00                   4
CH2CO             L 5/90C   2H   2O   1     G   200.000  3500.000  1000.000    1
 4.51129732E+00 9.00359745E-03-4.16939635E-06 9.23345882E-10-7.94838201E-14    2
-7.55105311E+03 6.32247205E-01 2.13583630E+00 1.81188721E-02-1.73947474E-05    3
 9.34397568E-09-2.01457615E-12-7.04291804E+03 1.22156480E+01                   4
HCCO              SRIC91H   1C   2O   1     G   300.00   4000.00   1000.000    1
 0.56282058E+01 0.40853401E-02-0.15934547E-05 0.28626052E-09-0.19407832E-13    2
 0.19327215E+05-0.39302595E+01 0.22517214E+01 0.17655021E-01-0.23729101E-04    3
 0.17275759E-07-0.50664811E-11 0.20059449E+05 0.12490417E+02                   4
HCCOH              SRI91C   2O   1H   2     G   300.000  5000.000  1000.000    1
 0.59238291E+01 0.67923600E-02-0.25658564E-05 0.44987841E-09-0.29940101E-13    2
 0.72646260E+04-0.76017742E+01 0.12423733E+01 0.31072201E-01-0.50866864E-04    3
 0.43137131E-07-0.14014594E-10 0.80316143E+04 0.13874319E+02                   4
H2CN               41687H   2C   1N   1     G   300.00   4000.000  1000.000    1
 0.52097030E+01 0.29692911E-02-0.28555891E-06-0.16355500E-09 0.30432589E-13    2
 0.27677109E+05-0.44444780E+01 0.28516610E+01 0.56952331E-02 0.10711400E-05    3
-0.16226120E-08-0.23511081E-12 0.28637820E+05 0.89927511E+01                   4
HCN               GRI/98H   1C   1N   1     G   200.000  6000.000  1000.000    1
 0.38022392E+01 0.31464228E-02-0.10632185E-05 0.16619757E-09-0.97997570E-14    2
 0.14407292E+05 0.15754601E+01 0.22589886E+01 0.10051170E-01-0.13351763E-04    3
 0.10092349E-07-0.30089028E-11 0.14712633E+05 0.89164419E+01                   4
HNO               And93 H   1N   1O   1     G   200.000  6000.000  1000.000    1
 0.29792509E+01 0.34944059E-02-0.78549778E-06 0.57479594E-10-0.19335916E-15    2
 0.11750582E+05 0.86063728E+01 0.45334916E+01-0.56696171E-02 0.18473207E-04    3
-0.17137094E-07 0.55454573E-11 0.11548297E+05 0.17498417E+01                   4
N                 L 6/88N   1               G   200.000  6000.000  1000.000    1
 0.24159429E+01 0.17489065E-03-0.11902369E-06 0.30226245E-10-0.20360982E-14    2
 0.56133773E+05 0.46496096E+01 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.56104637E+05 0.41939087E+01                   4
NNH               T07/93N   2H   1          G   200.000  6000.000  1000.000    1
 0.37667544E+01 0.28915082E-02-0.10416620E-05 0.16842594E-09-0.10091896E-13    2
 0.28650697E+05 0.44705067E+01 0.43446927E+01-0.48497072E-02 0.20059459E-04    3
-0.21726464E-07 0.79469539E-11 0.28791973E+05 0.29779410E+01                   4
N2O               L 7/88N   2O   1          G   200.000  6000.000  1000.000    1
 0.48230729E+01 0.26270251E-02-0.95850874E-06 0.16000712E-09-0.97752303E-14    2
 0.80734048E+04-0.22017207E+01 0.22571502E+01 0.11304728E-01-0.13671319E-04    3
 0.96819806E-08-0.29307182E-11 0.87417744E+04 0.10757992E+02                   4
NH                And94 N   1H   1          G   200.000  6000.000  1000.000    1
 0.27836928E+01 0.13298430E-02-0.42478047E-06 0.78348501E-10-0.55044470E-14    2
 0.42120848E+05 0.57407799E+01 0.34929085E+01 0.31179198E-03-0.14890484E-05    3
 0.24816442E-08-0.10356967E-11 0.41880629E+05 0.18483278E+01                   4
NH2               And89 N   1H   2          G   200.000  6000.000  1000.000    1
 0.28347421E+01 0.32073082E-02-0.93390804E-06 0.13702953E-09-0.79206144E-14    2
 0.22171957E+05 0.65204163E+01 0.42040029E+01-0.21061385E-02 0.71068348E-05    3
-0.56115197E-08 0.16440717E-11 0.21885910E+05-0.14184248E+00                   4
NH3               J 6/77N   1H   3          G   200.000  6000.000  1000.000    1
 0.26344521E+01 0.56662560E-02-0.17278676E-05 0.23867161E-09-0.12578786E-13    2
-0.65446958E+04 0.65662928E+01 0.42860274E+01-0.46605230E-02 0.21718513E-04    3
-0.22808887E-07 0.82638046E-11-0.67417285E+04-0.62537277E+00                   4
NO                RUS 78N   1O   1          G   200.000  6000.000  1000.000    1
 0.32606056E+01 0.11911043E-02-0.42917048E-06 0.69457669E-10-0.40336099E-14    2
 0.99209746E+04 0.63693027E+01 0.42184763E+01-0.46389760E-02 0.11041022E-04    3
-0.93361354E-08 0.28035770E-11 0.98446230E+04 0.22808464E+01                   4
NO2               L 7/88N   1O   2          G   200.000  6000.000  1000.000    1
 0.48847542E+01 0.21723956E-02-0.82806906E-06 0.15747510E-09-0.10510895E-13    2
 0.23164983E+04-0.11741695E+00 0.39440312E+01-0.15854290E-02 0.16657812E-04    3
-0.20475426E-07 0.78350564E-11 0.28966179E+04 0.63119917E+01                   4
HCNO              BDEA94H   1N   1C   1O   1G   300.000  5000.000  1382.000    1
 6.59860456E+00 3.02778626E-03-1.07704346E-06 1.71666528E-10-1.01439391E-14    2
 1.79661339E+04-1.03306599E+01 2.64727989E+00 1.27505342E-02-1.04794236E-05    3
 4.41432836E-09-7.57521466E-13 1.92990252E+04 1.07332972E+01                   4
HOCN              BDEA94H   1N   1C   1O   1G   300.000  5000.000  1368.000    1
 5.89784885E+00 3.16789393E-03-1.11801064E-06 1.77243144E-10-1.04339177E-14    2
-3.70653331E+03-6.18167825E+00 3.78604952E+00 6.88667922E-03-3.21487864E-06    3
 5.17195767E-10 1.19360788E-14-2.82698400E+03 5.63292162E+00                   4
HNCO              BDEA94H   1N   1C   1O   1G   300.000  5000.000  1478.000    1
 6.22395134E+00 3.17864004E-03-1.09378755E-06 1.70735163E-10-9.95021955E-15    2
-1.66599344E+04-8.38224741E+00 3.63096317E+00 7.30282357E-03-2.28050003E-06    3
-6.61271298E-10 3.62235752E-13-1.55873636E+04 6.19457727E+00                   4
NCO               EA 93 N   1C   1O   1     G   200.000  6000.000  1000.000    1
 0.51521845E+01 0.23051761E-02-0.88033153E-06 0.14789098E-09-0.90977996E-14    2
 0.14004123E+05-0.25442660E+01 0.28269308E+01 0.88051688E-02-0.83866134E-05    3
 0.48016964E-08-0.13313595E-11 0.14682477E+05 0.95504646E+01                   4
CN                HBH92 C   1N   1          G   200.000  6000.000  1000.000    1
 0.37459805E+01 0.43450775E-04 0.29705984E-06-0.68651806E-10 0.44134173E-14    2
 0.51536188E+05 0.27867601E+01 0.36129351E+01-0.95551327E-03 0.21442977E-05    3
-0.31516323E-09-0.46430356E-12 0.51708340E+05 0.39804995E+01                   4
HCNN              SRI/94C   1N   2H   1     G   300.000  5000.000  1000.000    1
 0.58946362E+01 0.39895959E-02-0.15982380E-05 0.29249395E-09-0.20094686E-13    2
 0.53452941E+05-0.51030502E+01 0.25243194E+01 0.15960619E-01-0.18816354E-04    3
 0.12125540E-07-0.32357378E-11 0.54261984E+05 0.11675870E+02                   4
N2                121286N   2               G   300.000  5000.000  1000.000    1
 0.02926640E+02 0.14879768E-02-0.05684760E-05 0.10097038E-09-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.14082404E-02-0.03963222E-04    3
 0.05641515E-07-0.02444854E-10-0.10208999E+04 0.03950372E+02                   4
AR                120186AR  1               G   300.000  5000.000  1000.000    1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.04366000E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366000E+02                   4
C3H8              L 4/85C   3H   8          G   300.000  5000.000  1000.000    1
 0.75341368E+01 0.18872239E-01-0.62718491E-05 0.91475649E-09-0.47838069E-13    2
-0.16467516E+05-0.17892349E+02 0.93355381E+00 0.26424579E-01 0.61059727E-05    3
-0.21977499E-07 0.95149253E-11-0.13958520E+05 0.19201691E+02                   4
C3H7              L 9/84C   3H   7          G   300.000  5000.000  1000.000    1
 0.77026987E+01 0.16044203E-01-0.52833220E-05 0.76298590E-09-0.39392284E-13    2
 0.82984336E+04-0.15480180E+02 0.10515518E+01 0.25991980E-01 0.23800540E-05    3
-0.19609569E-07 0.93732470E-11 0.10631863E+05 0.21122559E+02                   4
CH3CHO            L 8/88C   2H   4O   1     G   200.000  6000.000  1000.000    1
 0.54041108E+01 0.11723059E-01-0.42263137E-05 0.68372451E-09-0.40984863E-13    2
-0.22593122E+05-0.34807917E+01 0.47294595E+01-0.31932858E-02 0.47534921E-04    3
-0.57458611E-07 0.21931112E-10-0.21572878E+05 0.41030159E+01                   4
CH2CHO            SAND86O   1H   3C   2     G   300.000  5000.000  1000.000    1
 0.05975670E+02 0.08130591E-01-0.02743624E-04 0.04070304E-08-0.02176017E-12    2
 0.04903218E+04-0.05045251E+02 0.03409062E+02 0.10738574E-01 0.01891492E-04    3
-0.07158583E-07 0.02867385E-10 0.15214766E+04 0.09558290E+02                   4
END





\\
\\
\\  This is the tran file
\\
\\
AR                 0   136.500     3.330     0.000     0.000     0.000
C                  0    71.400     3.298     0.000     0.000     0.000 ! *
C2                 1    97.530     3.621     0.000     1.760     4.000
C2O                1   232.400     3.828     0.000     0.000     1.000 ! *
CN2                1   232.400     3.828     0.000     0.000     1.000 ! OIS
C2H                1   209.000     4.100     0.000     0.000     2.500
C2H2               1   209.000     4.100     0.000     0.000     2.500
C2H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *
C2H3               2   209.000     4.100     0.000     0.000     1.000 ! *
C2H4               2   280.800     3.971     0.000     0.000     1.500
C2H5               2   252.300     4.302     0.000     0.000     1.500
C2H6               2   252.300     4.302     0.000     0.000     1.500
C2N                1   232.400     3.828     0.000     0.000     1.000 ! OIS
C2N2               1   349.000     4.361     0.000     0.000     1.000 ! OIS
C3H2               2   209.000     4.100     0.000     0.000     1.000 ! *
C3H4               1   252.000     4.760     0.000     0.000     1.000
C3H6               2   266.800     4.982     0.000     0.000     1.000
C3H7               2   266.800     4.982     0.000     0.000     1.000
C4H6               2   357.000     5.180     0.000     0.000     1.000
I*C3H7             2   266.800     4.982     0.000     0.000     1.000
N*C3H7             2   266.800     4.982     0.000     0.000     1.000
C3H8               2   266.800     4.982     0.000     0.000     1.000
C4H                1   357.000     5.180     0.000     0.000     1.000
C4H2               1   357.000     5.180     0.000     0.000     1.000
C4H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *
C4H8               2   357.000     5.176     0.000     0.000     1.000
C4H9               2   357.000     5.176     0.000     0.000     1.000
I*C4H9             2   357.000     5.176     0.000     0.000     1.000
C5H2               1   357.000     5.180     0.000     0.000     1.000
C5H3               1   357.000     5.180     0.000     0.000     1.000
C6H2               1   357.000     5.180     0.000     0.000     1.000
C6H5               2   412.300     5.349     0.000     0.000     1.000 ! JAM
C6H5O              2   450.000     5.500     0.000     0.000     1.000 ! JAM
C5H5OH             2   450.000     5.500     0.000     0.000     1.000 ! JAM
C6H6               2   412.300     5.349     0.000     0.000     1.000 ! SVE
C6H7               2   412.300     5.349     0.000     0.000     1.000 ! JAM
CH                 1    80.000     2.750     0.000     0.000     0.000
CH2                1   144.000     3.800     0.000     0.000     0.000
CH2(S)             1   144.000     3.800     0.000     0.000     0.000
CH2*               1   144.000     3.800     0.000     0.000     0.000
CH2CHCCH           2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CHCCH2          2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CHCH2           2   260.000     4.850     0.000     0.000     1.000 ! JAM
CH2CHCHCH          2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CHCHCH2         2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CO              2   436.000     3.970     0.000     0.000     2.000
CH2O               2   498.000     3.590     0.000     0.000     2.000
CH2OH              2   417.000     3.690     1.700     0.000     2.000
CH3                1   144.000     3.800     0.000     0.000     0.000
CH3CC              2   252.000     4.760     0.000     0.000     1.000 ! JAM
CH3CCCH2           2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH3CCCH3           2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH3CCH2            2   260.000     4.850     0.000     0.000     1.000 ! JAM
CH3CHCH            2   260.000     4.850     0.000     0.000     1.000 ! JAM
CH3CH2CCH          2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH3CHO             2   436.000     3.970     0.000     0.000     2.000
CH2CHO             2   436.000     3.970     0.000     0.000     2.000
CH3CO              2   436.000     3.970     0.000     0.000     2.000
CH3O               2   417.000     3.690     1.700     0.000     2.000
CH3OH              2   481.800     3.626     0.000     0.000     1.000 ! SVE
CH4                2   141.400     3.746     0.000     2.600    13.000
CH4O               2   417.000     3.690     1.700     0.000     2.000
CN                 1    75.000     3.856     0.000     0.000     1.000 ! OIS
CNC                1   232.400     3.828     0.000     0.000     1.000 ! OIS
CNN                1   232.400     3.828     0.000     0.000     1.000 ! OIS
CO                 1    98.100     3.650     0.000     1.950     1.800
CO2                1   244.000     3.763     0.000     2.650     2.100
H                  0   145.000     2.050     0.000     0.000     0.000
H2C4O              2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2                 1    38.000     2.920     0.000     0.790   280.000
H2CCCCH            2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2CCCCH2           2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2CCCH             2   252.000     4.760     0.000     0.000     1.000 ! JAM
H2CN               1   569.000     3.630     0.000     0.000     1.000 ! os/jm
H2NO               2   116.700     3.492     0.000     0.000     1.000 ! JAM
H2O                2   572.400     2.605     1.844     0.000     4.000
H2O2               2   107.400     3.458     0.000     0.000     3.800
HC2N2              1   349.000     4.361     0.000     0.000     1.000 ! OIS
HCCHCCH            2   357.000     5.180     0.000     0.000     1.000 ! JAM
HCCO               2   150.000     2.500     0.000     0.000     1.000 ! *
HCNN               2   150.000     2.500     0.000     0.000     1.000 ! *
HCCOH              2   436.000     3.970     0.000     0.000     2.000
HCN                1   569.000     3.630     0.000     0.000     1.000 ! OIS
HCO                2   498.000     3.590     0.000     0.000     0.000
HE                 0    10.200     2.576     0.000     0.000     0.000 ! *
HCNO               2   232.400     3.828     0.000     0.000     1.000 ! JAM
HOCN               2   232.400     3.828     0.000     0.000     1.000 ! JAM
HNCO               2   232.400     3.828     0.000     0.000     1.000 ! OIS
HNNO               2   232.400     3.828     0.000     0.000     1.000 ! *
HNO                2   116.700     3.492     0.000     0.000     1.000 ! *
HNOH               2   116.700     3.492     0.000     0.000     1.000 ! JAM
HO2                2   107.400     3.458     0.000     0.000     1.000 ! *
N                  0    71.400     3.298     0.000     0.000     0.000 ! *
N2                 1    97.530     3.621     0.000     1.760     4.000
N2H2               2    71.400     3.798     0.000     0.000     1.000 ! *
N2H3               2   200.000     3.900     0.000     0.000     1.000 ! *
N2H4               2   205.000     4.230     0.000     4.260     1.500
N2O                1   232.400     3.828     0.000     0.000     1.000 ! *
NCN                1   232.400     3.828     0.000     0.000     1.000 ! OIS
NCO                1   232.400     3.828     0.000     0.000     1.000 ! OIS
NH                 1    80.000     2.650     0.000     0.000     4.000
NH2                2    80.000     2.650     0.000     2.260     4.000
NH3                2   481.000     2.920     1.470     0.000    10.000
NNH                2    71.400     3.798     0.000     0.000     1.000 ! *
NO                 1    97.530     3.621     0.000     1.760     4.000
NCNO               2   232.400     3.828     0.000     0.000     1.000 ! OIS
NO2                2   200.000     3.500     0.000     0.000     1.000 ! *
O                  0    80.000     2.750     0.000     0.000     0.000
O2                 1   107.400     3.458     0.000     1.600     3.800
OH                 1    80.000     2.750     0.000     0.000     0.000

#endif
