/*  -*- C -*-  */
/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *                                   Marc Day
 *                    Lawrence Berkeley National Laboratory
 *                      (C) 1998-2003  All Rights Reserved
 *
 * <LicenseText>
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#if defined(BL_FORT_USE_UPPERCASE)
#define CKINDX CKINDX
#define CKINIT CKINIT
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
#elif defined(BL_FORT_USE_LOWERCASE)
#define CKINDX ckindx
#define CKINIT ckinit
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
#elif defined(BL_FORT_USE_UNDERSCORE)
#define CKINDX ckindx_
#define CKINIT ckinit_
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
#endif

/*function declarations */
extern "C" {
void molecularWeight(double * wt);
void gibbs(double * species, double * tc);
void helmholtz(double * species, double * tc);
void speciesInternalEnergy(double * species, double * tc);
void speciesEnthalpy(double * species, double * tc);
void speciesEntropy(double * species, double * tc);
void cp_R(double * species, double * tc);
void cv_R(double * species, double * tc);
void equilibriumConstants(double * kc, double * g_RT, double T);
void productionRate(double * wdot, double * sc, double T);
void progressRate(double * qdot, double * speciesConc, double T);
void CKINDX(int * iwrk, double *rwrk, int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * rval, int * kerr, int lenline, int lenkray);
void CKSYME(int * kname, int * lenkname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(int * ickwrk, double * rckwrk, double * ru, double * ruc, double * pa);
void CKPX(double * rho, double * T, double * x, int * iwrk, double *rwrk, double * P);
void CKPY(double * rho, double * T, double * y, int * iwrk, double *rwrk, double * P);
void CKPC(double * rho, double * T, double * c, int * iwrk, double *rwrk, double * P);
void CKRHOX(double * P, double * T, double * x, int * iwrk, double *rwrk, double * rho);
void CKRHOY(double * P, double * T, double * y, int * iwrk, double *rwrk, double * rho);
void CKRHOC(double * P, double * T, double * c, int * iwrk, double *rwrk, double * rho);
void CKWT(int * iwrk, double *rwrk, double * wt);
void CKMMWY(double * y, int * iwrk, double * rwrk, double * wtm);
void CKMMWX(double * x, int * iwrk, double * rwrk, double * wtm);
void CKMMWC(double * c, int * iwrk, double * rwrk, double * wtm);
void CKYTX(double * y, int * iwrk, double * rwrk, double * x);
void CKYTCP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * c);
void CKYTCR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * c);
void CKXTY(double * x, int * iwrk, double * rwrk, double * y);
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c);
void CKXTCR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * c);
void CKCTX(double * c, int * iwrk, double * rwrk, double * x);
void CKCTY(double * c, int * iwrk, double * rwrk, double * y);
void CKCPOR(double * T, int * iwrk, double * rwrk, double * cpor);
void CKHORT(double * T, int * iwrk, double * rwrk, double * hort);
void CKSOR(double * T, int * iwrk, double * rwrk, double * sor);
void CKCVML(double * T, int * iwrk, double * rwrk, double * cvml);
void CKCPML(double * T, int * iwrk, double * rwrk, double * cvml);
void CKUML(double * T, int * iwrk, double * rwrk, double * uml);
void CKHML(double * T, int * iwrk, double * rwrk, double * uml);
void CKGML(double * T, int * iwrk, double * rwrk, double * gml);
void CKAML(double * T, int * iwrk, double * rwrk, double * aml);
void CKSML(double * T, int * iwrk, double * rwrk, double * sml);
void CKCVMS(double * T, int * iwrk, double * rwrk, double * cvms);
void CKCPMS(double * T, int * iwrk, double * rwrk, double * cvms);
void CKUMS(double * T, int * iwrk, double * rwrk, double * ums);
void CKHMS(double * T, int * iwrk, double * rwrk, double * ums);
void CKGMS(double * T, int * iwrk, double * rwrk, double * gms);
void CKAMS(double * T, int * iwrk, double * rwrk, double * ams);
void CKSMS(double * T, int * iwrk, double * rwrk, double * sms);
void CKCPBL(double * T, double * x, int * iwrk, double * rwrk, double * cpbl);
void CKCPBS(double * T, double * y, int * iwrk, double * rwrk, double * cpbs);
void CKCVBL(double * T, double * x, int * iwrk, double * rwrk, double * cpbl);
void CKCVBS(double * T, double * y, int * iwrk, double * rwrk, double * cpbs);
void CKHBML(double * T, double * x, int * iwrk, double * rwrk, double * hbml);
void CKHBMS(double * T, double * y, int * iwrk, double * rwrk, double * hbms);
void CKUBML(double * T, double * x, int * iwrk, double * rwrk, double * ubml);
void CKUBMS(double * T, double * y, int * iwrk, double * rwrk, double * ubms);
void CKSBML(double * P, double * T, double * x, int * iwrk, double * rwrk, double * sbml);
void CKSBMS(double * P, double * T, double * y, int * iwrk, double * rwrk, double * sbms);
void CKGBML(double * P, double * T, double * x, int * iwrk, double * rwrk, double * gbml);
void CKGBMS(double * P, double * T, double * y, int * iwrk, double * rwrk, double * gbms);
void CKABML(double * P, double * T, double * x, int * iwrk, double * rwrk, double * abml);
void CKABMS(double * P, double * T, double * y, int * iwrk, double * rwrk, double * abms);
void CKWC(double * T, double * C, int * iwrk, double *rwrk, double * wdot);
void CKWYP(double * P, double * T, double * y, int * iwrk, double *rwrk, double * wdot);
void CKWXP(double * P, double * T, double * x, int * iwrk, double *rwrk, double * wdot);
void CKWYR(double * rho, double * T, double * y, int * iwrk, double *rwrk, double * wdot);
void CKWXR(double * rho, double * T, double * x, int * iwrk, double *rwrk, double * wdot);
void CKQC(double * T, double * C, int * iwrk, double *rwrk, double * qdot);
void CKQYP(double * P, double * T, double * y, int * iwrk, double *rwrk, double * qdot);
void CKQXP(double * P, double * T, double * x, int * iwrk, double *rwrk, double * qdot);
void CKQYR(double * rho, double * T, double * y, int * iwrk, double *rwrk, double * qdot);
void CKQXR(double * rho, double * T, double * x, int * iwrk, double *rwrk, double * qdot);
void CKNU(int * kdim, int * iwrk, double *rwrk, int * nuki);
void CKNCF(int * mdim, int * iwrk, double *rwrk, int * ncf);
void CKABE(int * iwrk, double *rwrk, double * a, double * b, double * e );
void CKEQC(double * T, double * C , int * iwrk, double *rwrk, double * eqcon );
void CKEQYP(double * P, double * T, double * y, int * iwrk, double *rwrk, double * eqcon);
void CKEQXP(double * P, double * T, double * x, int * iwrk, double *rwrk, double * eqcon);
void CKEQYR(double * rho, double * T, double * y, int * iwrk, double *rwrk, double * eqcon);
void CKEQXR(double * rho, double * T, double * x, int * iwrk, double *rwrk, double * eqcon);
int  feeytt_(double * e, double * y, int * iwrk, double *rwrk, double * t);
void fephity_(double * phi, int * iwrk, double *rwrk, double * y);
void feytphi_(double * y, int * iwrk, double *rwrk, double * phi);
void fectyr_(double * c, double * rho, int * iwrk, double *rwrk, double * y);
void fecvrhs_(double * time, double * phi, double * phidot, double * rckwrk, int * ickwrk);
int fecvdim_();
void fezndrhs_(double * time, double * z, double * zdot, double * rckwrk, int * ickwrk);
int feznddim_();
char* femechfile_();
char* fesymname_(int sn);
int fesymnum_(const char* s1);
}


/*A few mechanism parameters */
void CKINDX(int * iwrk, double * rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 5;
    *kk = 39;
    *ii = 175;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/*Dummy ckinit */
void fginit_(int * leniwk, int * lenrwk, int * lencwk, int * linc, int * lout, int * ickwrk, double * rckwrk, char * cckwrk )
{
    if ((*lout) != 0) {
        printf(" ***       Congratulations       *** \n");
        printf(" * You are using the Fuego Library * \n");
        printf(" *****    Say NO to cklib.f    ***** \n");
    }
}


/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * rval, int * kerr, int lenline )
{
    int n,i; /*Loop Counters */
    char *p; /*String Tokens */
    char cstr[1000];
    /* Strip Comments  */
    for (i=0; i<lenline; ++i) {
        if (line[i]=='!') {
            cstr[i] = '\0';
            break;
        }
        cstr[i] = line[i];
    }

    p = strtok(cstr," ");
    if (!p) {
        *nval = 0;
        *kerr = 1;
        return;
    }
    for (n=0; n<*nexp; ++n) {
        rval[n] = atof(p);
        p = strtok(NULL, " ");
        if (!p) break;
    }
    *nval = n+1;
    if (*nval < *nexp) *kerr = 1;
    return;
}


/* cksnum... for parsing strings  */
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * rval, int * kerr, int lenline, int lenkray)
{
    /*Not done yet ... */
}


/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*5; i++) {
        kname[i] = ' ';
    }

    /* N  */
    kname[ 0*lenkname + 0 ] = 'N';
    kname[ 0*lenkname + 1 ] = ' ';

    /* AR  */
    kname[ 1*lenkname + 0 ] = 'A';
    kname[ 1*lenkname + 1 ] = 'R';

    /* H  */
    kname[ 2*lenkname + 0 ] = 'H';
    kname[ 2*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 3*lenkname + 0 ] = 'O';
    kname[ 3*lenkname + 1 ] = ' ';

    /* C  */
    kname[ 4*lenkname + 0 ] = 'C';
    kname[ 4*lenkname + 1 ] = ' ';
}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*39; i++) {
        kname[i] = ' ';
    }

    /* N2  */
    kname[ 0*lenkname + 0 ] = 'N';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* AR  */
    kname[ 1*lenkname + 0 ] = 'A';
    kname[ 1*lenkname + 1 ] = 'R';
    kname[ 1*lenkname + 2 ] = ' ';

    /* H  */
    kname[ 2*lenkname + 0 ] = 'H';
    kname[ 2*lenkname + 1 ] = ' ';

    /* O2  */
    kname[ 3*lenkname + 0 ] = 'O';
    kname[ 3*lenkname + 1 ] = '2';
    kname[ 3*lenkname + 2 ] = ' ';

    /* OH  */
    kname[ 4*lenkname + 0 ] = 'O';
    kname[ 4*lenkname + 1 ] = 'H';
    kname[ 4*lenkname + 2 ] = ' ';

    /* O  */
    kname[ 5*lenkname + 0 ] = 'O';
    kname[ 5*lenkname + 1 ] = ' ';

    /* H2  */
    kname[ 6*lenkname + 0 ] = 'H';
    kname[ 6*lenkname + 1 ] = '2';
    kname[ 6*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 7*lenkname + 0 ] = 'H';
    kname[ 7*lenkname + 1 ] = '2';
    kname[ 7*lenkname + 2 ] = 'O';
    kname[ 7*lenkname + 3 ] = ' ';

    /* HO2  */
    kname[ 8*lenkname + 0 ] = 'H';
    kname[ 8*lenkname + 1 ] = 'O';
    kname[ 8*lenkname + 2 ] = '2';
    kname[ 8*lenkname + 3 ] = ' ';

    /* H2O2  */
    kname[ 9*lenkname + 0 ] = 'H';
    kname[ 9*lenkname + 1 ] = '2';
    kname[ 9*lenkname + 2 ] = 'O';
    kname[ 9*lenkname + 3 ] = '2';
    kname[ 9*lenkname + 4 ] = ' ';

    /* CO  */
    kname[ 10*lenkname + 0 ] = 'C';
    kname[ 10*lenkname + 1 ] = 'O';
    kname[ 10*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 11*lenkname + 0 ] = 'C';
    kname[ 11*lenkname + 1 ] = 'O';
    kname[ 11*lenkname + 2 ] = '2';
    kname[ 11*lenkname + 3 ] = ' ';

    /* HCO  */
    kname[ 12*lenkname + 0 ] = 'H';
    kname[ 12*lenkname + 1 ] = 'C';
    kname[ 12*lenkname + 2 ] = 'O';
    kname[ 12*lenkname + 3 ] = ' ';

    /* CH2O  */
    kname[ 13*lenkname + 0 ] = 'C';
    kname[ 13*lenkname + 1 ] = 'H';
    kname[ 13*lenkname + 2 ] = '2';
    kname[ 13*lenkname + 3 ] = 'O';
    kname[ 13*lenkname + 4 ] = ' ';

    /* CH4  */
    kname[ 14*lenkname + 0 ] = 'C';
    kname[ 14*lenkname + 1 ] = 'H';
    kname[ 14*lenkname + 2 ] = '4';
    kname[ 14*lenkname + 3 ] = ' ';

    /* CH3  */
    kname[ 15*lenkname + 0 ] = 'C';
    kname[ 15*lenkname + 1 ] = 'H';
    kname[ 15*lenkname + 2 ] = '3';
    kname[ 15*lenkname + 3 ] = ' ';

    /* T-CH2  */
    kname[ 16*lenkname + 0 ] = 'T';
    kname[ 16*lenkname + 1 ] = '-';
    kname[ 16*lenkname + 2 ] = 'C';
    kname[ 16*lenkname + 3 ] = 'H';
    kname[ 16*lenkname + 4 ] = '2';
    kname[ 16*lenkname + 5 ] = ' ';

    /* S-CH2  */
    kname[ 17*lenkname + 0 ] = 'S';
    kname[ 17*lenkname + 1 ] = '-';
    kname[ 17*lenkname + 2 ] = 'C';
    kname[ 17*lenkname + 3 ] = 'H';
    kname[ 17*lenkname + 4 ] = '2';
    kname[ 17*lenkname + 5 ] = ' ';

    /* C2H4  */
    kname[ 18*lenkname + 0 ] = 'C';
    kname[ 18*lenkname + 1 ] = '2';
    kname[ 18*lenkname + 2 ] = 'H';
    kname[ 18*lenkname + 3 ] = '4';
    kname[ 18*lenkname + 4 ] = ' ';

    /* CH3O  */
    kname[ 19*lenkname + 0 ] = 'C';
    kname[ 19*lenkname + 1 ] = 'H';
    kname[ 19*lenkname + 2 ] = '3';
    kname[ 19*lenkname + 3 ] = 'O';
    kname[ 19*lenkname + 4 ] = ' ';

    /* C2H5  */
    kname[ 20*lenkname + 0 ] = 'C';
    kname[ 20*lenkname + 1 ] = '2';
    kname[ 20*lenkname + 2 ] = 'H';
    kname[ 20*lenkname + 3 ] = '5';
    kname[ 20*lenkname + 4 ] = ' ';

    /* C2H6  */
    kname[ 21*lenkname + 0 ] = 'C';
    kname[ 21*lenkname + 1 ] = '2';
    kname[ 21*lenkname + 2 ] = 'H';
    kname[ 21*lenkname + 3 ] = '6';
    kname[ 21*lenkname + 4 ] = ' ';

    /* CH  */
    kname[ 22*lenkname + 0 ] = 'C';
    kname[ 22*lenkname + 1 ] = 'H';
    kname[ 22*lenkname + 2 ] = ' ';

    /* C2H2  */
    kname[ 23*lenkname + 0 ] = 'C';
    kname[ 23*lenkname + 1 ] = '2';
    kname[ 23*lenkname + 2 ] = 'H';
    kname[ 23*lenkname + 3 ] = '2';
    kname[ 23*lenkname + 4 ] = ' ';

    /* C2H3  */
    kname[ 24*lenkname + 0 ] = 'C';
    kname[ 24*lenkname + 1 ] = '2';
    kname[ 24*lenkname + 2 ] = 'H';
    kname[ 24*lenkname + 3 ] = '3';
    kname[ 24*lenkname + 4 ] = ' ';

    /* CH2CHO  */
    kname[ 25*lenkname + 0 ] = 'C';
    kname[ 25*lenkname + 1 ] = 'H';
    kname[ 25*lenkname + 2 ] = '2';
    kname[ 25*lenkname + 3 ] = 'C';
    kname[ 25*lenkname + 4 ] = 'H';
    kname[ 25*lenkname + 5 ] = 'O';
    kname[ 25*lenkname + 6 ] = ' ';

    /* C2H4O  */
    kname[ 26*lenkname + 0 ] = 'C';
    kname[ 26*lenkname + 1 ] = '2';
    kname[ 26*lenkname + 2 ] = 'H';
    kname[ 26*lenkname + 3 ] = '4';
    kname[ 26*lenkname + 4 ] = 'O';
    kname[ 26*lenkname + 5 ] = ' ';

    /* CH2CO  */
    kname[ 27*lenkname + 0 ] = 'C';
    kname[ 27*lenkname + 1 ] = 'H';
    kname[ 27*lenkname + 2 ] = '2';
    kname[ 27*lenkname + 3 ] = 'C';
    kname[ 27*lenkname + 4 ] = 'O';
    kname[ 27*lenkname + 5 ] = ' ';

    /* HCCO  */
    kname[ 28*lenkname + 0 ] = 'H';
    kname[ 28*lenkname + 1 ] = 'C';
    kname[ 28*lenkname + 2 ] = 'C';
    kname[ 28*lenkname + 3 ] = 'O';
    kname[ 28*lenkname + 4 ] = ' ';

    /* C2H  */
    kname[ 29*lenkname + 0 ] = 'C';
    kname[ 29*lenkname + 1 ] = '2';
    kname[ 29*lenkname + 2 ] = 'H';
    kname[ 29*lenkname + 3 ] = ' ';

    /* CH2OH  */
    kname[ 30*lenkname + 0 ] = 'C';
    kname[ 30*lenkname + 1 ] = 'H';
    kname[ 30*lenkname + 2 ] = '2';
    kname[ 30*lenkname + 3 ] = 'O';
    kname[ 30*lenkname + 4 ] = 'H';
    kname[ 30*lenkname + 5 ] = ' ';

    /* CH3OH  */
    kname[ 31*lenkname + 0 ] = 'C';
    kname[ 31*lenkname + 1 ] = 'H';
    kname[ 31*lenkname + 2 ] = '3';
    kname[ 31*lenkname + 3 ] = 'O';
    kname[ 31*lenkname + 4 ] = 'H';
    kname[ 31*lenkname + 5 ] = ' ';

    /* C3H4  */
    kname[ 32*lenkname + 0 ] = 'C';
    kname[ 32*lenkname + 1 ] = '3';
    kname[ 32*lenkname + 2 ] = 'H';
    kname[ 32*lenkname + 3 ] = '4';
    kname[ 32*lenkname + 4 ] = ' ';

    /* C3H3  */
    kname[ 33*lenkname + 0 ] = 'C';
    kname[ 33*lenkname + 1 ] = '3';
    kname[ 33*lenkname + 2 ] = 'H';
    kname[ 33*lenkname + 3 ] = '3';
    kname[ 33*lenkname + 4 ] = ' ';

    /* C3H5  */
    kname[ 34*lenkname + 0 ] = 'C';
    kname[ 34*lenkname + 1 ] = '3';
    kname[ 34*lenkname + 2 ] = 'H';
    kname[ 34*lenkname + 3 ] = '5';
    kname[ 34*lenkname + 4 ] = ' ';

    /* C3H6  */
    kname[ 35*lenkname + 0 ] = 'C';
    kname[ 35*lenkname + 1 ] = '3';
    kname[ 35*lenkname + 2 ] = 'H';
    kname[ 35*lenkname + 3 ] = '6';
    kname[ 35*lenkname + 4 ] = ' ';

    /* C3H8  */
    kname[ 36*lenkname + 0 ] = 'C';
    kname[ 36*lenkname + 1 ] = '3';
    kname[ 36*lenkname + 2 ] = 'H';
    kname[ 36*lenkname + 3 ] = '8';
    kname[ 36*lenkname + 4 ] = ' ';

    /* I-C3H7  */
    kname[ 37*lenkname + 0 ] = 'I';
    kname[ 37*lenkname + 1 ] = '-';
    kname[ 37*lenkname + 2 ] = 'C';
    kname[ 37*lenkname + 3 ] = '3';
    kname[ 37*lenkname + 4 ] = 'H';
    kname[ 37*lenkname + 5 ] = '7';
    kname[ 37*lenkname + 6 ] = ' ';

    /* N-C3H7  */
    kname[ 38*lenkname + 0 ] = 'N';
    kname[ 38*lenkname + 1 ] = '-';
    kname[ 38*lenkname + 2 ] = 'C';
    kname[ 38*lenkname + 3 ] = '3';
    kname[ 38*lenkname + 4 ] = 'H';
    kname[ 38*lenkname + 5 ] = '7';
    kname[ 38*lenkname + 6 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(int * ickwrk, double * rckwrk, double * ru, double * ruc, double * pa)
{
     *ru  = 8.314e+07; 
     *ruc = 1.987; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*39.948000; /*AR */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*2.015940; /*H2 */
    XW += x[7]*18.015340; /*H2O */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*34.014740; /*H2O2 */
    XW += x[10]*28.010400; /*CO */
    XW += x[11]*44.009800; /*CO2 */
    XW += x[12]*29.018370; /*HCO */
    XW += x[13]*30.026340; /*CH2O */
    XW += x[14]*16.042880; /*CH4 */
    XW += x[15]*15.034910; /*CH3 */
    XW += x[16]*14.026940; /*T-CH2 */
    XW += x[17]*14.026940; /*S-CH2 */
    XW += x[18]*28.053880; /*C2H4 */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*29.061850; /*C2H5 */
    XW += x[21]*30.069820; /*C2H6 */
    XW += x[22]*13.018970; /*CH */
    XW += x[23]*26.037940; /*C2H2 */
    XW += x[24]*27.045910; /*C2H3 */
    XW += x[25]*43.045310; /*CH2CHO */
    XW += x[26]*44.053280; /*C2H4O */
    XW += x[27]*42.037340; /*CH2CO */
    XW += x[28]*41.029370; /*HCCO */
    XW += x[29]*25.029970; /*C2H */
    XW += x[30]*31.034310; /*CH2OH */
    XW += x[31]*32.042280; /*CH3OH */
    XW += x[32]*40.064880; /*C3H4 */
    XW += x[33]*39.056910; /*C3H3 */
    XW += x[34]*41.072850; /*C3H5 */
    XW += x[35]*42.080820; /*C3H6 */
    XW += x[36]*44.096760; /*C3H8 */
    XW += x[37]*43.088790; /*I-C3H7 */
    XW += x[38]*43.088790; /*N-C3H7 */
    *P = *rho * 8.314e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]/28.013400; /*N2 */
    YOW += y[1]/39.948000; /*AR */
    YOW += y[2]/1.007970; /*H */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/15.999400; /*O */
    YOW += y[6]/2.015940; /*H2 */
    YOW += y[7]/18.015340; /*H2O */
    YOW += y[8]/33.006770; /*HO2 */
    YOW += y[9]/34.014740; /*H2O2 */
    YOW += y[10]/28.010400; /*CO */
    YOW += y[11]/44.009800; /*CO2 */
    YOW += y[12]/29.018370; /*HCO */
    YOW += y[13]/30.026340; /*CH2O */
    YOW += y[14]/16.042880; /*CH4 */
    YOW += y[15]/15.034910; /*CH3 */
    YOW += y[16]/14.026940; /*T-CH2 */
    YOW += y[17]/14.026940; /*S-CH2 */
    YOW += y[18]/28.053880; /*C2H4 */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/29.061850; /*C2H5 */
    YOW += y[21]/30.069820; /*C2H6 */
    YOW += y[22]/13.018970; /*CH */
    YOW += y[23]/26.037940; /*C2H2 */
    YOW += y[24]/27.045910; /*C2H3 */
    YOW += y[25]/43.045310; /*CH2CHO */
    YOW += y[26]/44.053280; /*C2H4O */
    YOW += y[27]/42.037340; /*CH2CO */
    YOW += y[28]/41.029370; /*HCCO */
    YOW += y[29]/25.029970; /*C2H */
    YOW += y[30]/31.034310; /*CH2OH */
    YOW += y[31]/32.042280; /*CH3OH */
    YOW += y[32]/40.064880; /*C3H4 */
    YOW += y[33]/39.056910; /*C3H3 */
    YOW += y[34]/41.072850; /*C3H5 */
    YOW += y[35]/42.080820; /*C3H6 */
    YOW += y[36]/44.096760; /*C3H8 */
    YOW += y[37]/43.088790; /*I-C3H7 */
    YOW += y[38]/43.088790; /*N-C3H7 */
    *P = *rho * 8.314e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(c) */
void CKPC(double * rho, double * T, double * c, int * iwrk, double * rwrk, double * P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*28.013400; /*N2 */
    W += c[1]*39.948000; /*AR */
    W += c[2]*1.007970; /*H */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*15.999400; /*O */
    W += c[6]*2.015940; /*H2 */
    W += c[7]*18.015340; /*H2O */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*34.014740; /*H2O2 */
    W += c[10]*28.010400; /*CO */
    W += c[11]*44.009800; /*CO2 */
    W += c[12]*29.018370; /*HCO */
    W += c[13]*30.026340; /*CH2O */
    W += c[14]*16.042880; /*CH4 */
    W += c[15]*15.034910; /*CH3 */
    W += c[16]*14.026940; /*T-CH2 */
    W += c[17]*14.026940; /*S-CH2 */
    W += c[18]*28.053880; /*C2H4 */
    W += c[19]*31.034310; /*CH3O */
    W += c[20]*29.061850; /*C2H5 */
    W += c[21]*30.069820; /*C2H6 */
    W += c[22]*13.018970; /*CH */
    W += c[23]*26.037940; /*C2H2 */
    W += c[24]*27.045910; /*C2H3 */
    W += c[25]*43.045310; /*CH2CHO */
    W += c[26]*44.053280; /*C2H4O */
    W += c[27]*42.037340; /*CH2CO */
    W += c[28]*41.029370; /*HCCO */
    W += c[29]*25.029970; /*C2H */
    W += c[30]*31.034310; /*CH2OH */
    W += c[31]*32.042280; /*CH3OH */
    W += c[32]*40.064880; /*C3H4 */
    W += c[33]*39.056910; /*C3H3 */
    W += c[34]*41.072850; /*C3H5 */
    W += c[35]*42.080820; /*C3H6 */
    W += c[36]*44.096760; /*C3H8 */
    W += c[37]*43.088790; /*I-C3H7 */
    W += c[38]*43.088790; /*N-C3H7 */

    for (id = 0; id < 39; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.314e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * P, double * T, double * x, int * iwrk, double * rwrk, double * rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*39.948000; /*AR */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*2.015940; /*H2 */
    XW += x[7]*18.015340; /*H2O */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*34.014740; /*H2O2 */
    XW += x[10]*28.010400; /*CO */
    XW += x[11]*44.009800; /*CO2 */
    XW += x[12]*29.018370; /*HCO */
    XW += x[13]*30.026340; /*CH2O */
    XW += x[14]*16.042880; /*CH4 */
    XW += x[15]*15.034910; /*CH3 */
    XW += x[16]*14.026940; /*T-CH2 */
    XW += x[17]*14.026940; /*S-CH2 */
    XW += x[18]*28.053880; /*C2H4 */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*29.061850; /*C2H5 */
    XW += x[21]*30.069820; /*C2H6 */
    XW += x[22]*13.018970; /*CH */
    XW += x[23]*26.037940; /*C2H2 */
    XW += x[24]*27.045910; /*C2H3 */
    XW += x[25]*43.045310; /*CH2CHO */
    XW += x[26]*44.053280; /*C2H4O */
    XW += x[27]*42.037340; /*CH2CO */
    XW += x[28]*41.029370; /*HCCO */
    XW += x[29]*25.029970; /*C2H */
    XW += x[30]*31.034310; /*CH2OH */
    XW += x[31]*32.042280; /*CH3OH */
    XW += x[32]*40.064880; /*C3H4 */
    XW += x[33]*39.056910; /*C3H3 */
    XW += x[34]*41.072850; /*C3H5 */
    XW += x[35]*42.080820; /*C3H6 */
    XW += x[36]*44.096760; /*C3H8 */
    XW += x[37]*43.088790; /*I-C3H7 */
    XW += x[38]*43.088790; /*N-C3H7 */
    *rho = *P * XW / (8.314e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * P, double * T, double * y, int * iwrk, double * rwrk, double * rho)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]/28.013400; /*N2 */
    YOW += y[1]/39.948000; /*AR */
    YOW += y[2]/1.007970; /*H */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/15.999400; /*O */
    YOW += y[6]/2.015940; /*H2 */
    YOW += y[7]/18.015340; /*H2O */
    YOW += y[8]/33.006770; /*HO2 */
    YOW += y[9]/34.014740; /*H2O2 */
    YOW += y[10]/28.010400; /*CO */
    YOW += y[11]/44.009800; /*CO2 */
    YOW += y[12]/29.018370; /*HCO */
    YOW += y[13]/30.026340; /*CH2O */
    YOW += y[14]/16.042880; /*CH4 */
    YOW += y[15]/15.034910; /*CH3 */
    YOW += y[16]/14.026940; /*T-CH2 */
    YOW += y[17]/14.026940; /*S-CH2 */
    YOW += y[18]/28.053880; /*C2H4 */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/29.061850; /*C2H5 */
    YOW += y[21]/30.069820; /*C2H6 */
    YOW += y[22]/13.018970; /*CH */
    YOW += y[23]/26.037940; /*C2H2 */
    YOW += y[24]/27.045910; /*C2H3 */
    YOW += y[25]/43.045310; /*CH2CHO */
    YOW += y[26]/44.053280; /*C2H4O */
    YOW += y[27]/42.037340; /*CH2CO */
    YOW += y[28]/41.029370; /*HCCO */
    YOW += y[29]/25.029970; /*C2H */
    YOW += y[30]/31.034310; /*CH2OH */
    YOW += y[31]/32.042280; /*CH3OH */
    YOW += y[32]/40.064880; /*C3H4 */
    YOW += y[33]/39.056910; /*C3H3 */
    YOW += y[34]/41.072850; /*C3H5 */
    YOW += y[35]/42.080820; /*C3H6 */
    YOW += y[36]/44.096760; /*C3H8 */
    YOW += y[37]/43.088790; /*I-C3H7 */
    YOW += y[38]/43.088790; /*N-C3H7 */
    *rho = *P / (8.314e+07 * (*T) * YOW); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double * P, double * T, double * c, int * iwrk, double * rwrk, double * rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*28.013400; /*N2 */
    W += c[1]*39.948000; /*AR */
    W += c[2]*1.007970; /*H */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*15.999400; /*O */
    W += c[6]*2.015940; /*H2 */
    W += c[7]*18.015340; /*H2O */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*34.014740; /*H2O2 */
    W += c[10]*28.010400; /*CO */
    W += c[11]*44.009800; /*CO2 */
    W += c[12]*29.018370; /*HCO */
    W += c[13]*30.026340; /*CH2O */
    W += c[14]*16.042880; /*CH4 */
    W += c[15]*15.034910; /*CH3 */
    W += c[16]*14.026940; /*T-CH2 */
    W += c[17]*14.026940; /*S-CH2 */
    W += c[18]*28.053880; /*C2H4 */
    W += c[19]*31.034310; /*CH3O */
    W += c[20]*29.061850; /*C2H5 */
    W += c[21]*30.069820; /*C2H6 */
    W += c[22]*13.018970; /*CH */
    W += c[23]*26.037940; /*C2H2 */
    W += c[24]*27.045910; /*C2H3 */
    W += c[25]*43.045310; /*CH2CHO */
    W += c[26]*44.053280; /*C2H4O */
    W += c[27]*42.037340; /*CH2CO */
    W += c[28]*41.029370; /*HCCO */
    W += c[29]*25.029970; /*C2H */
    W += c[30]*31.034310; /*CH2OH */
    W += c[31]*32.042280; /*CH3OH */
    W += c[32]*40.064880; /*C3H4 */
    W += c[33]*39.056910; /*C3H3 */
    W += c[34]*41.072850; /*C3H5 */
    W += c[35]*42.080820; /*C3H6 */
    W += c[36]*44.096760; /*C3H8 */
    W += c[37]*43.088790; /*I-C3H7 */
    W += c[38]*43.088790; /*N-C3H7 */

    for (id = 0; id < 39; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.314e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT(int * iwrk, double * rwrk, double * wt)
{
    molecularWeight(wt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWY(double *y, int * iwrk, double * rwrk, double * wtm)
{
    double YOW = 0;/* see Eq 3 in CK Manual */
    YOW += y[0]/28.013400; /*N2 */
    YOW += y[1]/39.948000; /*AR */
    YOW += y[2]/1.007970; /*H */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/15.999400; /*O */
    YOW += y[6]/2.015940; /*H2 */
    YOW += y[7]/18.015340; /*H2O */
    YOW += y[8]/33.006770; /*HO2 */
    YOW += y[9]/34.014740; /*H2O2 */
    YOW += y[10]/28.010400; /*CO */
    YOW += y[11]/44.009800; /*CO2 */
    YOW += y[12]/29.018370; /*HCO */
    YOW += y[13]/30.026340; /*CH2O */
    YOW += y[14]/16.042880; /*CH4 */
    YOW += y[15]/15.034910; /*CH3 */
    YOW += y[16]/14.026940; /*T-CH2 */
    YOW += y[17]/14.026940; /*S-CH2 */
    YOW += y[18]/28.053880; /*C2H4 */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/29.061850; /*C2H5 */
    YOW += y[21]/30.069820; /*C2H6 */
    YOW += y[22]/13.018970; /*CH */
    YOW += y[23]/26.037940; /*C2H2 */
    YOW += y[24]/27.045910; /*C2H3 */
    YOW += y[25]/43.045310; /*CH2CHO */
    YOW += y[26]/44.053280; /*C2H4O */
    YOW += y[27]/42.037340; /*CH2CO */
    YOW += y[28]/41.029370; /*HCCO */
    YOW += y[29]/25.029970; /*C2H */
    YOW += y[30]/31.034310; /*CH2OH */
    YOW += y[31]/32.042280; /*CH3OH */
    YOW += y[32]/40.064880; /*C3H4 */
    YOW += y[33]/39.056910; /*C3H3 */
    YOW += y[34]/41.072850; /*C3H5 */
    YOW += y[35]/42.080820; /*C3H6 */
    YOW += y[36]/44.096760; /*C3H8 */
    YOW += y[37]/43.088790; /*I-C3H7 */
    YOW += y[38]/43.088790; /*N-C3H7 */
    *wtm = 1.0 / YOW;

    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double *x, int * iwrk, double * rwrk, double * wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*39.948000; /*AR */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*2.015940; /*H2 */
    XW += x[7]*18.015340; /*H2O */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*34.014740; /*H2O2 */
    XW += x[10]*28.010400; /*CO */
    XW += x[11]*44.009800; /*CO2 */
    XW += x[12]*29.018370; /*HCO */
    XW += x[13]*30.026340; /*CH2O */
    XW += x[14]*16.042880; /*CH4 */
    XW += x[15]*15.034910; /*CH3 */
    XW += x[16]*14.026940; /*T-CH2 */
    XW += x[17]*14.026940; /*S-CH2 */
    XW += x[18]*28.053880; /*C2H4 */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*29.061850; /*C2H5 */
    XW += x[21]*30.069820; /*C2H6 */
    XW += x[22]*13.018970; /*CH */
    XW += x[23]*26.037940; /*C2H2 */
    XW += x[24]*27.045910; /*C2H3 */
    XW += x[25]*43.045310; /*CH2CHO */
    XW += x[26]*44.053280; /*C2H4O */
    XW += x[27]*42.037340; /*CH2CO */
    XW += x[28]*41.029370; /*HCCO */
    XW += x[29]*25.029970; /*C2H */
    XW += x[30]*31.034310; /*CH2OH */
    XW += x[31]*32.042280; /*CH3OH */
    XW += x[32]*40.064880; /*C3H4 */
    XW += x[33]*39.056910; /*C3H3 */
    XW += x[34]*41.072850; /*C3H5 */
    XW += x[35]*42.080820; /*C3H6 */
    XW += x[36]*44.096760; /*C3H8 */
    XW += x[37]*43.088790; /*I-C3H7 */
    XW += x[38]*43.088790; /*N-C3H7 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void CKMMWC(double *c, int * iwrk, double * rwrk, double * wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*28.013400; /*N2 */
    W += c[1]*39.948000; /*AR */
    W += c[2]*1.007970; /*H */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*15.999400; /*O */
    W += c[6]*2.015940; /*H2 */
    W += c[7]*18.015340; /*H2O */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*34.014740; /*H2O2 */
    W += c[10]*28.010400; /*CO */
    W += c[11]*44.009800; /*CO2 */
    W += c[12]*29.018370; /*HCO */
    W += c[13]*30.026340; /*CH2O */
    W += c[14]*16.042880; /*CH4 */
    W += c[15]*15.034910; /*CH3 */
    W += c[16]*14.026940; /*T-CH2 */
    W += c[17]*14.026940; /*S-CH2 */
    W += c[18]*28.053880; /*C2H4 */
    W += c[19]*31.034310; /*CH3O */
    W += c[20]*29.061850; /*C2H5 */
    W += c[21]*30.069820; /*C2H6 */
    W += c[22]*13.018970; /*CH */
    W += c[23]*26.037940; /*C2H2 */
    W += c[24]*27.045910; /*C2H3 */
    W += c[25]*43.045310; /*CH2CHO */
    W += c[26]*44.053280; /*C2H4O */
    W += c[27]*42.037340; /*CH2CO */
    W += c[28]*41.029370; /*HCCO */
    W += c[29]*25.029970; /*C2H */
    W += c[30]*31.034310; /*CH2OH */
    W += c[31]*32.042280; /*CH3OH */
    W += c[32]*40.064880; /*C3H4 */
    W += c[33]*39.056910; /*C3H3 */
    W += c[34]*41.072850; /*C3H5 */
    W += c[35]*42.080820; /*C3H6 */
    W += c[36]*44.096760; /*C3H8 */
    W += c[37]*43.088790; /*I-C3H7 */
    W += c[38]*43.088790; /*N-C3H7 */

    for (id = 0; id < 39; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
void CKYTX(double * y, int * iwrk, double * rwrk, double * x)
{
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/28.013400; /*N2 */
    YOW += y[1]/39.948000; /*AR */
    YOW += y[2]/1.007970; /*H */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/15.999400; /*O */
    YOW += y[6]/2.015940; /*H2 */
    YOW += y[7]/18.015340; /*H2O */
    YOW += y[8]/33.006770; /*HO2 */
    YOW += y[9]/34.014740; /*H2O2 */
    YOW += y[10]/28.010400; /*CO */
    YOW += y[11]/44.009800; /*CO2 */
    YOW += y[12]/29.018370; /*HCO */
    YOW += y[13]/30.026340; /*CH2O */
    YOW += y[14]/16.042880; /*CH4 */
    YOW += y[15]/15.034910; /*CH3 */
    YOW += y[16]/14.026940; /*T-CH2 */
    YOW += y[17]/14.026940; /*S-CH2 */
    YOW += y[18]/28.053880; /*C2H4 */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/29.061850; /*C2H5 */
    YOW += y[21]/30.069820; /*C2H6 */
    YOW += y[22]/13.018970; /*CH */
    YOW += y[23]/26.037940; /*C2H2 */
    YOW += y[24]/27.045910; /*C2H3 */
    YOW += y[25]/43.045310; /*CH2CHO */
    YOW += y[26]/44.053280; /*C2H4O */
    YOW += y[27]/42.037340; /*CH2CO */
    YOW += y[28]/41.029370; /*HCCO */
    YOW += y[29]/25.029970; /*C2H */
    YOW += y[30]/31.034310; /*CH2OH */
    YOW += y[31]/32.042280; /*CH3OH */
    YOW += y[32]/40.064880; /*C3H4 */
    YOW += y[33]/39.056910; /*C3H3 */
    YOW += y[34]/41.072850; /*C3H5 */
    YOW += y[35]/42.080820; /*C3H6 */
    YOW += y[36]/44.096760; /*C3H8 */
    YOW += y[37]/43.088790; /*I-C3H7 */
    YOW += y[38]/43.088790; /*N-C3H7 */
    /*Now compute conversion */
    x[0] = y[0]/(28.013400*YOW); 
    x[1] = y[1]/(39.948000*YOW); 
    x[2] = y[2]/(1.007970*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(2.015940*YOW); 
    x[7] = y[7]/(18.015340*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(34.014740*YOW); 
    x[10] = y[10]/(28.010400*YOW); 
    x[11] = y[11]/(44.009800*YOW); 
    x[12] = y[12]/(29.018370*YOW); 
    x[13] = y[13]/(30.026340*YOW); 
    x[14] = y[14]/(16.042880*YOW); 
    x[15] = y[15]/(15.034910*YOW); 
    x[16] = y[16]/(14.026940*YOW); 
    x[17] = y[17]/(14.026940*YOW); 
    x[18] = y[18]/(28.053880*YOW); 
    x[19] = y[19]/(31.034310*YOW); 
    x[20] = y[20]/(29.061850*YOW); 
    x[21] = y[21]/(30.069820*YOW); 
    x[22] = y[22]/(13.018970*YOW); 
    x[23] = y[23]/(26.037940*YOW); 
    x[24] = y[24]/(27.045910*YOW); 
    x[25] = y[25]/(43.045310*YOW); 
    x[26] = y[26]/(44.053280*YOW); 
    x[27] = y[27]/(42.037340*YOW); 
    x[28] = y[28]/(41.029370*YOW); 
    x[29] = y[29]/(25.029970*YOW); 
    x[30] = y[30]/(31.034310*YOW); 
    x[31] = y[31]/(32.042280*YOW); 
    x[32] = y[32]/(40.064880*YOW); 
    x[33] = y[33]/(39.056910*YOW); 
    x[34] = y[34]/(41.072850*YOW); 
    x[35] = y[35]/(42.080820*YOW); 
    x[36] = y[36]/(44.096760*YOW); 
    x[37] = y[37]/(43.088790*YOW); 
    x[38] = y[38]/(43.088790*YOW); 

    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * c)
{
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/28.013400; /*N2 */
    YOW += y[1]/39.948000; /*AR */
    YOW += y[2]/1.007970; /*H */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/15.999400; /*O */
    YOW += y[6]/2.015940; /*H2 */
    YOW += y[7]/18.015340; /*H2O */
    YOW += y[8]/33.006770; /*HO2 */
    YOW += y[9]/34.014740; /*H2O2 */
    YOW += y[10]/28.010400; /*CO */
    YOW += y[11]/44.009800; /*CO2 */
    YOW += y[12]/29.018370; /*HCO */
    YOW += y[13]/30.026340; /*CH2O */
    YOW += y[14]/16.042880; /*CH4 */
    YOW += y[15]/15.034910; /*CH3 */
    YOW += y[16]/14.026940; /*T-CH2 */
    YOW += y[17]/14.026940; /*S-CH2 */
    YOW += y[18]/28.053880; /*C2H4 */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/29.061850; /*C2H5 */
    YOW += y[21]/30.069820; /*C2H6 */
    YOW += y[22]/13.018970; /*CH */
    YOW += y[23]/26.037940; /*C2H2 */
    YOW += y[24]/27.045910; /*C2H3 */
    YOW += y[25]/43.045310; /*CH2CHO */
    YOW += y[26]/44.053280; /*C2H4O */
    YOW += y[27]/42.037340; /*CH2CO */
    YOW += y[28]/41.029370; /*HCCO */
    YOW += y[29]/25.029970; /*C2H */
    YOW += y[30]/31.034310; /*CH2OH */
    YOW += y[31]/32.042280; /*CH3OH */
    YOW += y[32]/40.064880; /*C3H4 */
    YOW += y[33]/39.056910; /*C3H3 */
    YOW += y[34]/41.072850; /*C3H5 */
    YOW += y[35]/42.080820; /*C3H6 */
    YOW += y[36]/44.096760; /*C3H8 */
    YOW += y[37]/43.088790; /*I-C3H7 */
    YOW += y[38]/43.088790; /*N-C3H7 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.314e+07 * (*T)); 
    /*Now compute conversion */
    c[0] = PWORT * y[0]/28.013400; 
    c[1] = PWORT * y[1]/39.948000; 
    c[2] = PWORT * y[2]/1.007970; 
    c[3] = PWORT * y[3]/31.998800; 
    c[4] = PWORT * y[4]/17.007370; 
    c[5] = PWORT * y[5]/15.999400; 
    c[6] = PWORT * y[6]/2.015940; 
    c[7] = PWORT * y[7]/18.015340; 
    c[8] = PWORT * y[8]/33.006770; 
    c[9] = PWORT * y[9]/34.014740; 
    c[10] = PWORT * y[10]/28.010400; 
    c[11] = PWORT * y[11]/44.009800; 
    c[12] = PWORT * y[12]/29.018370; 
    c[13] = PWORT * y[13]/30.026340; 
    c[14] = PWORT * y[14]/16.042880; 
    c[15] = PWORT * y[15]/15.034910; 
    c[16] = PWORT * y[16]/14.026940; 
    c[17] = PWORT * y[17]/14.026940; 
    c[18] = PWORT * y[18]/28.053880; 
    c[19] = PWORT * y[19]/31.034310; 
    c[20] = PWORT * y[20]/29.061850; 
    c[21] = PWORT * y[21]/30.069820; 
    c[22] = PWORT * y[22]/13.018970; 
    c[23] = PWORT * y[23]/26.037940; 
    c[24] = PWORT * y[24]/27.045910; 
    c[25] = PWORT * y[25]/43.045310; 
    c[26] = PWORT * y[26]/44.053280; 
    c[27] = PWORT * y[27]/42.037340; 
    c[28] = PWORT * y[28]/41.029370; 
    c[29] = PWORT * y[29]/25.029970; 
    c[30] = PWORT * y[30]/31.034310; 
    c[31] = PWORT * y[31]/32.042280; 
    c[32] = PWORT * y[32]/40.064880; 
    c[33] = PWORT * y[33]/39.056910; 
    c[34] = PWORT * y[34]/41.072850; 
    c[35] = PWORT * y[35]/42.080820; 
    c[36] = PWORT * y[36]/44.096760; 
    c[37] = PWORT * y[37]/43.088790; 
    c[38] = PWORT * y[38]/43.088790; 

    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * c)
{
    /*See Eq 8 (Temperature not used) */
    c[0] = (*rho) * y[0]/28.013400; 
    c[1] = (*rho) * y[1]/39.948000; 
    c[2] = (*rho) * y[2]/1.007970; 
    c[3] = (*rho) * y[3]/31.998800; 
    c[4] = (*rho) * y[4]/17.007370; 
    c[5] = (*rho) * y[5]/15.999400; 
    c[6] = (*rho) * y[6]/2.015940; 
    c[7] = (*rho) * y[7]/18.015340; 
    c[8] = (*rho) * y[8]/33.006770; 
    c[9] = (*rho) * y[9]/34.014740; 
    c[10] = (*rho) * y[10]/28.010400; 
    c[11] = (*rho) * y[11]/44.009800; 
    c[12] = (*rho) * y[12]/29.018370; 
    c[13] = (*rho) * y[13]/30.026340; 
    c[14] = (*rho) * y[14]/16.042880; 
    c[15] = (*rho) * y[15]/15.034910; 
    c[16] = (*rho) * y[16]/14.026940; 
    c[17] = (*rho) * y[17]/14.026940; 
    c[18] = (*rho) * y[18]/28.053880; 
    c[19] = (*rho) * y[19]/31.034310; 
    c[20] = (*rho) * y[20]/29.061850; 
    c[21] = (*rho) * y[21]/30.069820; 
    c[22] = (*rho) * y[22]/13.018970; 
    c[23] = (*rho) * y[23]/26.037940; 
    c[24] = (*rho) * y[24]/27.045910; 
    c[25] = (*rho) * y[25]/43.045310; 
    c[26] = (*rho) * y[26]/44.053280; 
    c[27] = (*rho) * y[27]/42.037340; 
    c[28] = (*rho) * y[28]/41.029370; 
    c[29] = (*rho) * y[29]/25.029970; 
    c[30] = (*rho) * y[30]/31.034310; 
    c[31] = (*rho) * y[31]/32.042280; 
    c[32] = (*rho) * y[32]/40.064880; 
    c[33] = (*rho) * y[33]/39.056910; 
    c[34] = (*rho) * y[34]/41.072850; 
    c[35] = (*rho) * y[35]/42.080820; 
    c[36] = (*rho) * y[36]/44.096760; 
    c[37] = (*rho) * y[37]/43.088790; 
    c[38] = (*rho) * y[38]/43.088790; 

    return;
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * x, int * iwrk, double * rwrk, double * y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*39.948000; /*AR */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*2.015940; /*H2 */
    XW += x[7]*18.015340; /*H2O */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*34.014740; /*H2O2 */
    XW += x[10]*28.010400; /*CO */
    XW += x[11]*44.009800; /*CO2 */
    XW += x[12]*29.018370; /*HCO */
    XW += x[13]*30.026340; /*CH2O */
    XW += x[14]*16.042880; /*CH4 */
    XW += x[15]*15.034910; /*CH3 */
    XW += x[16]*14.026940; /*T-CH2 */
    XW += x[17]*14.026940; /*S-CH2 */
    XW += x[18]*28.053880; /*C2H4 */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*29.061850; /*C2H5 */
    XW += x[21]*30.069820; /*C2H6 */
    XW += x[22]*13.018970; /*CH */
    XW += x[23]*26.037940; /*C2H2 */
    XW += x[24]*27.045910; /*C2H3 */
    XW += x[25]*43.045310; /*CH2CHO */
    XW += x[26]*44.053280; /*C2H4O */
    XW += x[27]*42.037340; /*CH2CO */
    XW += x[28]*41.029370; /*HCCO */
    XW += x[29]*25.029970; /*C2H */
    XW += x[30]*31.034310; /*CH2OH */
    XW += x[31]*32.042280; /*CH3OH */
    XW += x[32]*40.064880; /*C3H4 */
    XW += x[33]*39.056910; /*C3H3 */
    XW += x[34]*41.072850; /*C3H5 */
    XW += x[35]*42.080820; /*C3H6 */
    XW += x[36]*44.096760; /*C3H8 */
    XW += x[37]*43.088790; /*I-C3H7 */
    XW += x[38]*43.088790; /*N-C3H7 */
    /*Now compute conversion */
    y[0] = x[0]*28.013400/XW; 
    y[1] = x[1]*39.948000/XW; 
    y[2] = x[2]*1.007970/XW; 
    y[3] = x[3]*31.998800/XW; 
    y[4] = x[4]*17.007370/XW; 
    y[5] = x[5]*15.999400/XW; 
    y[6] = x[6]*2.015940/XW; 
    y[7] = x[7]*18.015340/XW; 
    y[8] = x[8]*33.006770/XW; 
    y[9] = x[9]*34.014740/XW; 
    y[10] = x[10]*28.010400/XW; 
    y[11] = x[11]*44.009800/XW; 
    y[12] = x[12]*29.018370/XW; 
    y[13] = x[13]*30.026340/XW; 
    y[14] = x[14]*16.042880/XW; 
    y[15] = x[15]*15.034910/XW; 
    y[16] = x[16]*14.026940/XW; 
    y[17] = x[17]*14.026940/XW; 
    y[18] = x[18]*28.053880/XW; 
    y[19] = x[19]*31.034310/XW; 
    y[20] = x[20]*29.061850/XW; 
    y[21] = x[21]*30.069820/XW; 
    y[22] = x[22]*13.018970/XW; 
    y[23] = x[23]*26.037940/XW; 
    y[24] = x[24]*27.045910/XW; 
    y[25] = x[25]*43.045310/XW; 
    y[26] = x[26]*44.053280/XW; 
    y[27] = x[27]*42.037340/XW; 
    y[28] = x[28]*41.029370/XW; 
    y[29] = x[29]*25.029970/XW; 
    y[30] = x[30]*31.034310/XW; 
    y[31] = x[31]*32.042280/XW; 
    y[32] = x[32]*40.064880/XW; 
    y[33] = x[33]*39.056910/XW; 
    y[34] = x[34]*41.072850/XW; 
    y[35] = x[35]*42.080820/XW; 
    y[36] = x[36]*44.096760/XW; 
    y[37] = x[37]*43.088790/XW; 
    y[38] = x[38]*43.088790/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.314e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 39; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*39.948000; /*AR */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*2.015940; /*H2 */
    XW += x[7]*18.015340; /*H2O */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*34.014740; /*H2O2 */
    XW += x[10]*28.010400; /*CO */
    XW += x[11]*44.009800; /*CO2 */
    XW += x[12]*29.018370; /*HCO */
    XW += x[13]*30.026340; /*CH2O */
    XW += x[14]*16.042880; /*CH4 */
    XW += x[15]*15.034910; /*CH3 */
    XW += x[16]*14.026940; /*T-CH2 */
    XW += x[17]*14.026940; /*S-CH2 */
    XW += x[18]*28.053880; /*C2H4 */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*29.061850; /*C2H5 */
    XW += x[21]*30.069820; /*C2H6 */
    XW += x[22]*13.018970; /*CH */
    XW += x[23]*26.037940; /*C2H2 */
    XW += x[24]*27.045910; /*C2H3 */
    XW += x[25]*43.045310; /*CH2CHO */
    XW += x[26]*44.053280; /*C2H4O */
    XW += x[27]*42.037340; /*CH2CO */
    XW += x[28]*41.029370; /*HCCO */
    XW += x[29]*25.029970; /*C2H */
    XW += x[30]*31.034310; /*CH2OH */
    XW += x[31]*32.042280; /*CH3OH */
    XW += x[32]*40.064880; /*C3H4 */
    XW += x[33]*39.056910; /*C3H3 */
    XW += x[34]*41.072850; /*C3H5 */
    XW += x[35]*42.080820; /*C3H6 */
    XW += x[36]*44.096760; /*C3H8 */
    XW += x[37]*43.088790; /*I-C3H7 */
    XW += x[38]*43.088790; /*N-C3H7 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 39; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
void CKCTX(double * c, int * iwrk, double * rwrk, double * x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 39; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    for (id = 0; id < 39; ++id) {
        x[id] = c[id]/sumC;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double * c, int * iwrk, double * rwrk, double * y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*28.013400; /*N2 */
    CW += c[1]*39.948000; /*AR */
    CW += c[2]*1.007970; /*H */
    CW += c[3]*31.998800; /*O2 */
    CW += c[4]*17.007370; /*OH */
    CW += c[5]*15.999400; /*O */
    CW += c[6]*2.015940; /*H2 */
    CW += c[7]*18.015340; /*H2O */
    CW += c[8]*33.006770; /*HO2 */
    CW += c[9]*34.014740; /*H2O2 */
    CW += c[10]*28.010400; /*CO */
    CW += c[11]*44.009800; /*CO2 */
    CW += c[12]*29.018370; /*HCO */
    CW += c[13]*30.026340; /*CH2O */
    CW += c[14]*16.042880; /*CH4 */
    CW += c[15]*15.034910; /*CH3 */
    CW += c[16]*14.026940; /*T-CH2 */
    CW += c[17]*14.026940; /*S-CH2 */
    CW += c[18]*28.053880; /*C2H4 */
    CW += c[19]*31.034310; /*CH3O */
    CW += c[20]*29.061850; /*C2H5 */
    CW += c[21]*30.069820; /*C2H6 */
    CW += c[22]*13.018970; /*CH */
    CW += c[23]*26.037940; /*C2H2 */
    CW += c[24]*27.045910; /*C2H3 */
    CW += c[25]*43.045310; /*CH2CHO */
    CW += c[26]*44.053280; /*C2H4O */
    CW += c[27]*42.037340; /*CH2CO */
    CW += c[28]*41.029370; /*HCCO */
    CW += c[29]*25.029970; /*C2H */
    CW += c[30]*31.034310; /*CH2OH */
    CW += c[31]*32.042280; /*CH3OH */
    CW += c[32]*40.064880; /*C3H4 */
    CW += c[33]*39.056910; /*C3H3 */
    CW += c[34]*41.072850; /*C3H5 */
    CW += c[35]*42.080820; /*C3H6 */
    CW += c[36]*44.096760; /*C3H8 */
    CW += c[37]*43.088790; /*I-C3H7 */
    CW += c[38]*43.088790; /*N-C3H7 */
    /*Now compute conversion */
    y[0] = c[0]*28.013400/CW; 
    y[1] = c[1]*39.948000/CW; 
    y[2] = c[2]*1.007970/CW; 
    y[3] = c[3]*31.998800/CW; 
    y[4] = c[4]*17.007370/CW; 
    y[5] = c[5]*15.999400/CW; 
    y[6] = c[6]*2.015940/CW; 
    y[7] = c[7]*18.015340/CW; 
    y[8] = c[8]*33.006770/CW; 
    y[9] = c[9]*34.014740/CW; 
    y[10] = c[10]*28.010400/CW; 
    y[11] = c[11]*44.009800/CW; 
    y[12] = c[12]*29.018370/CW; 
    y[13] = c[13]*30.026340/CW; 
    y[14] = c[14]*16.042880/CW; 
    y[15] = c[15]*15.034910/CW; 
    y[16] = c[16]*14.026940/CW; 
    y[17] = c[17]*14.026940/CW; 
    y[18] = c[18]*28.053880/CW; 
    y[19] = c[19]*31.034310/CW; 
    y[20] = c[20]*29.061850/CW; 
    y[21] = c[21]*30.069820/CW; 
    y[22] = c[22]*13.018970/CW; 
    y[23] = c[23]*26.037940/CW; 
    y[24] = c[24]*27.045910/CW; 
    y[25] = c[25]*43.045310/CW; 
    y[26] = c[26]*44.053280/CW; 
    y[27] = c[27]*42.037340/CW; 
    y[28] = c[28]*41.029370/CW; 
    y[29] = c[29]*25.029970/CW; 
    y[30] = c[30]*31.034310/CW; 
    y[31] = c[31]*32.042280/CW; 
    y[32] = c[32]*40.064880/CW; 
    y[33] = c[33]*39.056910/CW; 
    y[34] = c[34]*41.072850/CW; 
    y[35] = c[35]*42.080820/CW; 
    y[36] = c[36]*44.096760/CW; 
    y[37] = c[37]*43.088790/CW; 
    y[38] = c[38]*43.088790/CW; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void CKCPOR(double *T, int * iwrk, double * rwrk, double * cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void CKHORT(double *T, int * iwrk, double * rwrk, double * hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void CKSOR(double *T, int * iwrk, double * rwrk, double * sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void CKCVML(double *T, int * iwrk, double * rwrk, double * cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        cvml[id] *= 8.314e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void CKCPML(double *T, int * iwrk, double * rwrk, double * cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        cpml[id] *= 8.314e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double *T, int * iwrk, double * rwrk, double * uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void CKHML(double *T, int * iwrk, double * rwrk, double * hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void CKGML(double *T, int * iwrk, double * rwrk, double * gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void CKAML(double *T, int * iwrk, double * rwrk, double * aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void CKSML(double *T, int * iwrk, double * rwrk, double * sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        sml[id] *= 8.314e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void CKCVMS(double *T, int * iwrk, double * rwrk, double * cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 2967865.378712; /*N2 */
    cvms[1] *= 2081205.567237; /*AR */
    cvms[2] *= 82482613.569848; /*H */
    cvms[3] *= 2598222.433341; /*O2 */
    cvms[4] *= 4888468.940230; /*OH */
    cvms[5] *= 5196444.866683; /*O */
    cvms[6] *= 41241306.784924; /*H2 */
    cvms[7] *= 4614955.920899; /*H2O */
    cvms[8] *= 2518877.187922; /*HO2 */
    cvms[9] *= 2444234.470115; /*H2O2 */
    cvms[10] *= 2968183.246223; /*CO */
    cvms[11] *= 1889124.694954; /*CO2 */
    cvms[12] *= 2865081.670680; /*HCO */
    cvms[13] *= 2768902.237169; /*CH2O */
    cvms[14] *= 5182361.271792; /*CH4 */
    cvms[15] *= 5529796.985815; /*CH3 */
    cvms[16] *= 5927165.867966; /*T-CH2 */
    cvms[17] *= 5927165.867966; /*S-CH2 */
    cvms[18] *= 2963582.933983; /*C2H4 */
    cvms[19] *= 2678970.468491; /*CH3O */
    cvms[20] *= 2860795.166171; /*C2H5 */
    cvms[21] *= 2764898.492908; /*C2H6 */
    cvms[22] *= 6386065.871570; /*CH */
    cvms[23] *= 3193032.935785; /*C2H2 */
    cvms[24] *= 3074032.265877; /*C2H3 */
    cvms[25] *= 1931453.159473; /*CH2CHO */
    cvms[26] *= 1887260.154068; /*C2H4O */
    cvms[27] *= 1977765.481831; /*CH2CO */
    cvms[28] *= 2026353.317148; /*HCCO */
    cvms[29] *= 3321618.044289; /*C2H */
    cvms[30] *= 2678970.468491; /*CH2OH */
    cvms[31] *= 2594696.756910; /*CH3OH */
    cvms[32] *= 2075134.132437; /*C3H4 */
    cvms[33] *= 2128688.623857; /*C3H3 */
    cvms[34] *= 2024208.205664; /*C3H5 */
    cvms[35] *= 1975721.955989; /*C3H6 */
    cvms[36] *= 1885399.290107; /*C3H8 */
    cvms[37] *= 1929504.170342; /*I-C3H7 */
    cvms[38] *= 1929504.170342; /*N-C3H7 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double *T, int * iwrk, double * rwrk, double * cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 2967865.378712; /*N2 */
    cpms[1] *= 2081205.567237; /*AR */
    cpms[2] *= 82482613.569848; /*H */
    cpms[3] *= 2598222.433341; /*O2 */
    cpms[4] *= 4888468.940230; /*OH */
    cpms[5] *= 5196444.866683; /*O */
    cpms[6] *= 41241306.784924; /*H2 */
    cpms[7] *= 4614955.920899; /*H2O */
    cpms[8] *= 2518877.187922; /*HO2 */
    cpms[9] *= 2444234.470115; /*H2O2 */
    cpms[10] *= 2968183.246223; /*CO */
    cpms[11] *= 1889124.694954; /*CO2 */
    cpms[12] *= 2865081.670680; /*HCO */
    cpms[13] *= 2768902.237169; /*CH2O */
    cpms[14] *= 5182361.271792; /*CH4 */
    cpms[15] *= 5529796.985815; /*CH3 */
    cpms[16] *= 5927165.867966; /*T-CH2 */
    cpms[17] *= 5927165.867966; /*S-CH2 */
    cpms[18] *= 2963582.933983; /*C2H4 */
    cpms[19] *= 2678970.468491; /*CH3O */
    cpms[20] *= 2860795.166171; /*C2H5 */
    cpms[21] *= 2764898.492908; /*C2H6 */
    cpms[22] *= 6386065.871570; /*CH */
    cpms[23] *= 3193032.935785; /*C2H2 */
    cpms[24] *= 3074032.265877; /*C2H3 */
    cpms[25] *= 1931453.159473; /*CH2CHO */
    cpms[26] *= 1887260.154068; /*C2H4O */
    cpms[27] *= 1977765.481831; /*CH2CO */
    cpms[28] *= 2026353.317148; /*HCCO */
    cpms[29] *= 3321618.044289; /*C2H */
    cpms[30] *= 2678970.468491; /*CH2OH */
    cpms[31] *= 2594696.756910; /*CH3OH */
    cpms[32] *= 2075134.132437; /*C3H4 */
    cpms[33] *= 2128688.623857; /*C3H3 */
    cpms[34] *= 2024208.205664; /*C3H5 */
    cpms[35] *= 1975721.955989; /*C3H6 */
    cpms[36] *= 1885399.290107; /*C3H8 */
    cpms[37] *= 1929504.170342; /*I-C3H7 */
    cpms[38] *= 1929504.170342; /*N-C3H7 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double *T, int * iwrk, double * rwrk, double * ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    ums[0] *= RT/28.013400; /*N2 */
    ums[1] *= RT/39.948000; /*AR */
    ums[2] *= RT/1.007970; /*H */
    ums[3] *= RT/31.998800; /*O2 */
    ums[4] *= RT/17.007370; /*OH */
    ums[5] *= RT/15.999400; /*O */
    ums[6] *= RT/2.015940; /*H2 */
    ums[7] *= RT/18.015340; /*H2O */
    ums[8] *= RT/33.006770; /*HO2 */
    ums[9] *= RT/34.014740; /*H2O2 */
    ums[10] *= RT/28.010400; /*CO */
    ums[11] *= RT/44.009800; /*CO2 */
    ums[12] *= RT/29.018370; /*HCO */
    ums[13] *= RT/30.026340; /*CH2O */
    ums[14] *= RT/16.042880; /*CH4 */
    ums[15] *= RT/15.034910; /*CH3 */
    ums[16] *= RT/14.026940; /*T-CH2 */
    ums[17] *= RT/14.026940; /*S-CH2 */
    ums[18] *= RT/28.053880; /*C2H4 */
    ums[19] *= RT/31.034310; /*CH3O */
    ums[20] *= RT/29.061850; /*C2H5 */
    ums[21] *= RT/30.069820; /*C2H6 */
    ums[22] *= RT/13.018970; /*CH */
    ums[23] *= RT/26.037940; /*C2H2 */
    ums[24] *= RT/27.045910; /*C2H3 */
    ums[25] *= RT/43.045310; /*CH2CHO */
    ums[26] *= RT/44.053280; /*C2H4O */
    ums[27] *= RT/42.037340; /*CH2CO */
    ums[28] *= RT/41.029370; /*HCCO */
    ums[29] *= RT/25.029970; /*C2H */
    ums[30] *= RT/31.034310; /*CH2OH */
    ums[31] *= RT/32.042280; /*CH3OH */
    ums[32] *= RT/40.064880; /*C3H4 */
    ums[33] *= RT/39.056910; /*C3H3 */
    ums[34] *= RT/41.072850; /*C3H5 */
    ums[35] *= RT/42.080820; /*C3H6 */
    ums[36] *= RT/44.096760; /*C3H8 */
    ums[37] *= RT/43.088790; /*I-C3H7 */
    ums[38] *= RT/43.088790; /*N-C3H7 */
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double *T, int * iwrk, double * rwrk, double * hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    hms[0] *= RT/28.013400; /*N2 */
    hms[1] *= RT/39.948000; /*AR */
    hms[2] *= RT/1.007970; /*H */
    hms[3] *= RT/31.998800; /*O2 */
    hms[4] *= RT/17.007370; /*OH */
    hms[5] *= RT/15.999400; /*O */
    hms[6] *= RT/2.015940; /*H2 */
    hms[7] *= RT/18.015340; /*H2O */
    hms[8] *= RT/33.006770; /*HO2 */
    hms[9] *= RT/34.014740; /*H2O2 */
    hms[10] *= RT/28.010400; /*CO */
    hms[11] *= RT/44.009800; /*CO2 */
    hms[12] *= RT/29.018370; /*HCO */
    hms[13] *= RT/30.026340; /*CH2O */
    hms[14] *= RT/16.042880; /*CH4 */
    hms[15] *= RT/15.034910; /*CH3 */
    hms[16] *= RT/14.026940; /*T-CH2 */
    hms[17] *= RT/14.026940; /*S-CH2 */
    hms[18] *= RT/28.053880; /*C2H4 */
    hms[19] *= RT/31.034310; /*CH3O */
    hms[20] *= RT/29.061850; /*C2H5 */
    hms[21] *= RT/30.069820; /*C2H6 */
    hms[22] *= RT/13.018970; /*CH */
    hms[23] *= RT/26.037940; /*C2H2 */
    hms[24] *= RT/27.045910; /*C2H3 */
    hms[25] *= RT/43.045310; /*CH2CHO */
    hms[26] *= RT/44.053280; /*C2H4O */
    hms[27] *= RT/42.037340; /*CH2CO */
    hms[28] *= RT/41.029370; /*HCCO */
    hms[29] *= RT/25.029970; /*C2H */
    hms[30] *= RT/31.034310; /*CH2OH */
    hms[31] *= RT/32.042280; /*CH3OH */
    hms[32] *= RT/40.064880; /*C3H4 */
    hms[33] *= RT/39.056910; /*C3H3 */
    hms[34] *= RT/41.072850; /*C3H5 */
    hms[35] *= RT/42.080820; /*C3H6 */
    hms[36] *= RT/44.096760; /*C3H8 */
    hms[37] *= RT/43.088790; /*I-C3H7 */
    hms[38] *= RT/43.088790; /*N-C3H7 */
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double *T, int * iwrk, double * rwrk, double * gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    gibbs(gms, tc);
    gms[0] *= RT/28.013400; /*N2 */
    gms[1] *= RT/39.948000; /*AR */
    gms[2] *= RT/1.007970; /*H */
    gms[3] *= RT/31.998800; /*O2 */
    gms[4] *= RT/17.007370; /*OH */
    gms[5] *= RT/15.999400; /*O */
    gms[6] *= RT/2.015940; /*H2 */
    gms[7] *= RT/18.015340; /*H2O */
    gms[8] *= RT/33.006770; /*HO2 */
    gms[9] *= RT/34.014740; /*H2O2 */
    gms[10] *= RT/28.010400; /*CO */
    gms[11] *= RT/44.009800; /*CO2 */
    gms[12] *= RT/29.018370; /*HCO */
    gms[13] *= RT/30.026340; /*CH2O */
    gms[14] *= RT/16.042880; /*CH4 */
    gms[15] *= RT/15.034910; /*CH3 */
    gms[16] *= RT/14.026940; /*T-CH2 */
    gms[17] *= RT/14.026940; /*S-CH2 */
    gms[18] *= RT/28.053880; /*C2H4 */
    gms[19] *= RT/31.034310; /*CH3O */
    gms[20] *= RT/29.061850; /*C2H5 */
    gms[21] *= RT/30.069820; /*C2H6 */
    gms[22] *= RT/13.018970; /*CH */
    gms[23] *= RT/26.037940; /*C2H2 */
    gms[24] *= RT/27.045910; /*C2H3 */
    gms[25] *= RT/43.045310; /*CH2CHO */
    gms[26] *= RT/44.053280; /*C2H4O */
    gms[27] *= RT/42.037340; /*CH2CO */
    gms[28] *= RT/41.029370; /*HCCO */
    gms[29] *= RT/25.029970; /*C2H */
    gms[30] *= RT/31.034310; /*CH2OH */
    gms[31] *= RT/32.042280; /*CH3OH */
    gms[32] *= RT/40.064880; /*C3H4 */
    gms[33] *= RT/39.056910; /*C3H3 */
    gms[34] *= RT/41.072850; /*C3H5 */
    gms[35] *= RT/42.080820; /*C3H6 */
    gms[36] *= RT/44.096760; /*C3H8 */
    gms[37] *= RT/43.088790; /*I-C3H7 */
    gms[38] *= RT/43.088790; /*N-C3H7 */
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double *T, int * iwrk, double * rwrk, double * ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    helmholtz(ams, tc);
    ams[0] *= RT/28.013400; /*N2 */
    ams[1] *= RT/39.948000; /*AR */
    ams[2] *= RT/1.007970; /*H */
    ams[3] *= RT/31.998800; /*O2 */
    ams[4] *= RT/17.007370; /*OH */
    ams[5] *= RT/15.999400; /*O */
    ams[6] *= RT/2.015940; /*H2 */
    ams[7] *= RT/18.015340; /*H2O */
    ams[8] *= RT/33.006770; /*HO2 */
    ams[9] *= RT/34.014740; /*H2O2 */
    ams[10] *= RT/28.010400; /*CO */
    ams[11] *= RT/44.009800; /*CO2 */
    ams[12] *= RT/29.018370; /*HCO */
    ams[13] *= RT/30.026340; /*CH2O */
    ams[14] *= RT/16.042880; /*CH4 */
    ams[15] *= RT/15.034910; /*CH3 */
    ams[16] *= RT/14.026940; /*T-CH2 */
    ams[17] *= RT/14.026940; /*S-CH2 */
    ams[18] *= RT/28.053880; /*C2H4 */
    ams[19] *= RT/31.034310; /*CH3O */
    ams[20] *= RT/29.061850; /*C2H5 */
    ams[21] *= RT/30.069820; /*C2H6 */
    ams[22] *= RT/13.018970; /*CH */
    ams[23] *= RT/26.037940; /*C2H2 */
    ams[24] *= RT/27.045910; /*C2H3 */
    ams[25] *= RT/43.045310; /*CH2CHO */
    ams[26] *= RT/44.053280; /*C2H4O */
    ams[27] *= RT/42.037340; /*CH2CO */
    ams[28] *= RT/41.029370; /*HCCO */
    ams[29] *= RT/25.029970; /*C2H */
    ams[30] *= RT/31.034310; /*CH2OH */
    ams[31] *= RT/32.042280; /*CH3OH */
    ams[32] *= RT/40.064880; /*C3H4 */
    ams[33] *= RT/39.056910; /*C3H3 */
    ams[34] *= RT/41.072850; /*C3H5 */
    ams[35] *= RT/42.080820; /*C3H6 */
    ams[36] *= RT/44.096760; /*C3H8 */
    ams[37] *= RT/43.088790; /*I-C3H7 */
    ams[38] *= RT/43.088790; /*N-C3H7 */
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double *T, int * iwrk, double * rwrk, double * sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 2967865.378712; /*N2 */
    sms[1] *= 2081205.567237; /*AR */
    sms[2] *= 82482613.569848; /*H */
    sms[3] *= 2598222.433341; /*O2 */
    sms[4] *= 4888468.940230; /*OH */
    sms[5] *= 5196444.866683; /*O */
    sms[6] *= 41241306.784924; /*H2 */
    sms[7] *= 4614955.920899; /*H2O */
    sms[8] *= 2518877.187922; /*HO2 */
    sms[9] *= 2444234.470115; /*H2O2 */
    sms[10] *= 2968183.246223; /*CO */
    sms[11] *= 1889124.694954; /*CO2 */
    sms[12] *= 2865081.670680; /*HCO */
    sms[13] *= 2768902.237169; /*CH2O */
    sms[14] *= 5182361.271792; /*CH4 */
    sms[15] *= 5529796.985815; /*CH3 */
    sms[16] *= 5927165.867966; /*T-CH2 */
    sms[17] *= 5927165.867966; /*S-CH2 */
    sms[18] *= 2963582.933983; /*C2H4 */
    sms[19] *= 2678970.468491; /*CH3O */
    sms[20] *= 2860795.166171; /*C2H5 */
    sms[21] *= 2764898.492908; /*C2H6 */
    sms[22] *= 6386065.871570; /*CH */
    sms[23] *= 3193032.935785; /*C2H2 */
    sms[24] *= 3074032.265877; /*C2H3 */
    sms[25] *= 1931453.159473; /*CH2CHO */
    sms[26] *= 1887260.154068; /*C2H4O */
    sms[27] *= 1977765.481831; /*CH2CO */
    sms[28] *= 2026353.317148; /*HCCO */
    sms[29] *= 3321618.044289; /*C2H */
    sms[30] *= 2678970.468491; /*CH2OH */
    sms[31] *= 2594696.756910; /*CH3OH */
    sms[32] *= 2075134.132437; /*C3H4 */
    sms[33] *= 2128688.623857; /*C3H3 */
    sms[34] *= 2024208.205664; /*C3H5 */
    sms[35] *= 1975721.955989; /*C3H6 */
    sms[36] *= 1885399.290107; /*C3H8 */
    sms[37] *= 1929504.170342; /*I-C3H7 */
    sms[38] *= 1929504.170342; /*N-C3H7 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *T, double *x, int * iwrk, double * rwrk, double * cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[39]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 39; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.314e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void CKCPBS(double *T, double *y, int * iwrk, double * rwrk, double * cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[39]; /* temporary storage */
    cp_R(cpor, tc);
    /*multiply by y/molecularweight */
    result += cpor[0]*y[0]/28.0134; /*N2 */
    result += cpor[1]*y[1]/39.948; /*AR */
    result += cpor[2]*y[2]/1.00797; /*H */
    result += cpor[3]*y[3]/31.9988; /*O2 */
    result += cpor[4]*y[4]/17.0074; /*OH */
    result += cpor[5]*y[5]/15.9994; /*O */
    result += cpor[6]*y[6]/2.01594; /*H2 */
    result += cpor[7]*y[7]/18.0153; /*H2O */
    result += cpor[8]*y[8]/33.0068; /*HO2 */
    result += cpor[9]*y[9]/34.0147; /*H2O2 */
    result += cpor[10]*y[10]/28.0104; /*CO */
    result += cpor[11]*y[11]/44.0098; /*CO2 */
    result += cpor[12]*y[12]/29.0184; /*HCO */
    result += cpor[13]*y[13]/30.0263; /*CH2O */
    result += cpor[14]*y[14]/16.0429; /*CH4 */
    result += cpor[15]*y[15]/15.0349; /*CH3 */
    result += cpor[16]*y[16]/14.0269; /*T-CH2 */
    result += cpor[17]*y[17]/14.0269; /*S-CH2 */
    result += cpor[18]*y[18]/28.0539; /*C2H4 */
    result += cpor[19]*y[19]/31.0343; /*CH3O */
    result += cpor[20]*y[20]/29.0618; /*C2H5 */
    result += cpor[21]*y[21]/30.0698; /*C2H6 */
    result += cpor[22]*y[22]/13.019; /*CH */
    result += cpor[23]*y[23]/26.0379; /*C2H2 */
    result += cpor[24]*y[24]/27.0459; /*C2H3 */
    result += cpor[25]*y[25]/43.0453; /*CH2CHO */
    result += cpor[26]*y[26]/44.0533; /*C2H4O */
    result += cpor[27]*y[27]/42.0373; /*CH2CO */
    result += cpor[28]*y[28]/41.0294; /*HCCO */
    result += cpor[29]*y[29]/25.03; /*C2H */
    result += cpor[30]*y[30]/31.0343; /*CH2OH */
    result += cpor[31]*y[31]/32.0423; /*CH3OH */
    result += cpor[32]*y[32]/40.0649; /*C3H4 */
    result += cpor[33]*y[33]/39.0569; /*C3H3 */
    result += cpor[34]*y[34]/41.0729; /*C3H5 */
    result += cpor[35]*y[35]/42.0808; /*C3H6 */
    result += cpor[36]*y[36]/44.0968; /*C3H8 */
    result += cpor[37]*y[37]/43.0888; /*I-C3H7 */
    result += cpor[38]*y[38]/43.0888; /*N-C3H7 */

    *cpbs = result * 8.314e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *T, double *x, int * iwrk, double * rwrk, double * cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[39]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 39; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.314e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
void CKCVBS(double *T, double *y, int * iwrk, double * rwrk, double * cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[39]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]/28.0134; /*N2 */
    result += cvor[1]*y[1]/39.948; /*AR */
    result += cvor[2]*y[2]/1.00797; /*H */
    result += cvor[3]*y[3]/31.9988; /*O2 */
    result += cvor[4]*y[4]/17.0074; /*OH */
    result += cvor[5]*y[5]/15.9994; /*O */
    result += cvor[6]*y[6]/2.01594; /*H2 */
    result += cvor[7]*y[7]/18.0153; /*H2O */
    result += cvor[8]*y[8]/33.0068; /*HO2 */
    result += cvor[9]*y[9]/34.0147; /*H2O2 */
    result += cvor[10]*y[10]/28.0104; /*CO */
    result += cvor[11]*y[11]/44.0098; /*CO2 */
    result += cvor[12]*y[12]/29.0184; /*HCO */
    result += cvor[13]*y[13]/30.0263; /*CH2O */
    result += cvor[14]*y[14]/16.0429; /*CH4 */
    result += cvor[15]*y[15]/15.0349; /*CH3 */
    result += cvor[16]*y[16]/14.0269; /*T-CH2 */
    result += cvor[17]*y[17]/14.0269; /*S-CH2 */
    result += cvor[18]*y[18]/28.0539; /*C2H4 */
    result += cvor[19]*y[19]/31.0343; /*CH3O */
    result += cvor[20]*y[20]/29.0618; /*C2H5 */
    result += cvor[21]*y[21]/30.0698; /*C2H6 */
    result += cvor[22]*y[22]/13.019; /*CH */
    result += cvor[23]*y[23]/26.0379; /*C2H2 */
    result += cvor[24]*y[24]/27.0459; /*C2H3 */
    result += cvor[25]*y[25]/43.0453; /*CH2CHO */
    result += cvor[26]*y[26]/44.0533; /*C2H4O */
    result += cvor[27]*y[27]/42.0373; /*CH2CO */
    result += cvor[28]*y[28]/41.0294; /*HCCO */
    result += cvor[29]*y[29]/25.03; /*C2H */
    result += cvor[30]*y[30]/31.0343; /*CH2OH */
    result += cvor[31]*y[31]/32.0423; /*CH3OH */
    result += cvor[32]*y[32]/40.0649; /*C3H4 */
    result += cvor[33]*y[33]/39.0569; /*C3H3 */
    result += cvor[34]*y[34]/41.0729; /*C3H5 */
    result += cvor[35]*y[35]/42.0808; /*C3H6 */
    result += cvor[36]*y[36]/44.0968; /*C3H8 */
    result += cvor[37]*y[37]/43.0888; /*I-C3H7 */
    result += cvor[38]*y[38]/43.0888; /*N-C3H7 */

    *cvbs = result * 8.314e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *T, double *x, int * iwrk, double * rwrk, double * hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[39]; /* temporary storage */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 39; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
void CKHBMS(double *T, double *y, int * iwrk, double * rwrk, double * hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[39]; /* temporary storage */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*hml[0]/28.013400; /*N2 */
    result += y[1]*hml[1]/39.948000; /*AR */
    result += y[2]*hml[2]/1.007970; /*H */
    result += y[3]*hml[3]/31.998800; /*O2 */
    result += y[4]*hml[4]/17.007370; /*OH */
    result += y[5]*hml[5]/15.999400; /*O */
    result += y[6]*hml[6]/2.015940; /*H2 */
    result += y[7]*hml[7]/18.015340; /*H2O */
    result += y[8]*hml[8]/33.006770; /*HO2 */
    result += y[9]*hml[9]/34.014740; /*H2O2 */
    result += y[10]*hml[10]/28.010400; /*CO */
    result += y[11]*hml[11]/44.009800; /*CO2 */
    result += y[12]*hml[12]/29.018370; /*HCO */
    result += y[13]*hml[13]/30.026340; /*CH2O */
    result += y[14]*hml[14]/16.042880; /*CH4 */
    result += y[15]*hml[15]/15.034910; /*CH3 */
    result += y[16]*hml[16]/14.026940; /*T-CH2 */
    result += y[17]*hml[17]/14.026940; /*S-CH2 */
    result += y[18]*hml[18]/28.053880; /*C2H4 */
    result += y[19]*hml[19]/31.034310; /*CH3O */
    result += y[20]*hml[20]/29.061850; /*C2H5 */
    result += y[21]*hml[21]/30.069820; /*C2H6 */
    result += y[22]*hml[22]/13.018970; /*CH */
    result += y[23]*hml[23]/26.037940; /*C2H2 */
    result += y[24]*hml[24]/27.045910; /*C2H3 */
    result += y[25]*hml[25]/43.045310; /*CH2CHO */
    result += y[26]*hml[26]/44.053280; /*C2H4O */
    result += y[27]*hml[27]/42.037340; /*CH2CO */
    result += y[28]*hml[28]/41.029370; /*HCCO */
    result += y[29]*hml[29]/25.029970; /*C2H */
    result += y[30]*hml[30]/31.034310; /*CH2OH */
    result += y[31]*hml[31]/32.042280; /*CH3OH */
    result += y[32]*hml[32]/40.064880; /*C3H4 */
    result += y[33]*hml[33]/39.056910; /*C3H3 */
    result += y[34]*hml[34]/41.072850; /*C3H5 */
    result += y[35]*hml[35]/42.080820; /*C3H6 */
    result += y[36]*hml[36]/44.096760; /*C3H8 */
    result += y[37]*hml[37]/43.088790; /*I-C3H7 */
    result += y[38]*hml[38]/43.088790; /*N-C3H7 */

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *T, double *x, int * iwrk, double * rwrk, double * ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[39]; /* temporary energy array */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 39; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
void CKUBMS(double *T, double *y, int * iwrk, double * rwrk, double * ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[39]; /* temporary energy array */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]/28.013400; /*N2 */
    result += y[1]*ums[1]/39.948000; /*AR */
    result += y[2]*ums[2]/1.007970; /*H */
    result += y[3]*ums[3]/31.998800; /*O2 */
    result += y[4]*ums[4]/17.007370; /*OH */
    result += y[5]*ums[5]/15.999400; /*O */
    result += y[6]*ums[6]/2.015940; /*H2 */
    result += y[7]*ums[7]/18.015340; /*H2O */
    result += y[8]*ums[8]/33.006770; /*HO2 */
    result += y[9]*ums[9]/34.014740; /*H2O2 */
    result += y[10]*ums[10]/28.010400; /*CO */
    result += y[11]*ums[11]/44.009800; /*CO2 */
    result += y[12]*ums[12]/29.018370; /*HCO */
    result += y[13]*ums[13]/30.026340; /*CH2O */
    result += y[14]*ums[14]/16.042880; /*CH4 */
    result += y[15]*ums[15]/15.034910; /*CH3 */
    result += y[16]*ums[16]/14.026940; /*T-CH2 */
    result += y[17]*ums[17]/14.026940; /*S-CH2 */
    result += y[18]*ums[18]/28.053880; /*C2H4 */
    result += y[19]*ums[19]/31.034310; /*CH3O */
    result += y[20]*ums[20]/29.061850; /*C2H5 */
    result += y[21]*ums[21]/30.069820; /*C2H6 */
    result += y[22]*ums[22]/13.018970; /*CH */
    result += y[23]*ums[23]/26.037940; /*C2H2 */
    result += y[24]*ums[24]/27.045910; /*C2H3 */
    result += y[25]*ums[25]/43.045310; /*CH2CHO */
    result += y[26]*ums[26]/44.053280; /*C2H4O */
    result += y[27]*ums[27]/42.037340; /*CH2CO */
    result += y[28]*ums[28]/41.029370; /*HCCO */
    result += y[29]*ums[29]/25.029970; /*C2H */
    result += y[30]*ums[30]/31.034310; /*CH2OH */
    result += y[31]*ums[31]/32.042280; /*CH3OH */
    result += y[32]*ums[32]/40.064880; /*C3H4 */
    result += y[33]*ums[33]/39.056910; /*C3H3 */
    result += y[34]*ums[34]/41.072850; /*C3H5 */
    result += y[35]*ums[35]/42.080820; /*C3H6 */
    result += y[36]*ums[36]/44.096760; /*C3H8 */
    result += y[37]*ums[37]/43.088790; /*I-C3H7 */
    result += y[38]*ums[38]/43.088790; /*N-C3H7 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
void CKSBML(double *P, double *T, double *x, int * iwrk, double * rwrk, double * sbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[39]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 39; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.314e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double *P, double *T, double *y, int * iwrk, double * rwrk, double * sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[39]; /* temporary storage */
    double x[39]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/28.013400; /*N2 */
    YOW += y[1]/39.948000; /*AR */
    YOW += y[2]/1.007970; /*H */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/15.999400; /*O */
    YOW += y[6]/2.015940; /*H2 */
    YOW += y[7]/18.015340; /*H2O */
    YOW += y[8]/33.006770; /*HO2 */
    YOW += y[9]/34.014740; /*H2O2 */
    YOW += y[10]/28.010400; /*CO */
    YOW += y[11]/44.009800; /*CO2 */
    YOW += y[12]/29.018370; /*HCO */
    YOW += y[13]/30.026340; /*CH2O */
    YOW += y[14]/16.042880; /*CH4 */
    YOW += y[15]/15.034910; /*CH3 */
    YOW += y[16]/14.026940; /*T-CH2 */
    YOW += y[17]/14.026940; /*S-CH2 */
    YOW += y[18]/28.053880; /*C2H4 */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/29.061850; /*C2H5 */
    YOW += y[21]/30.069820; /*C2H6 */
    YOW += y[22]/13.018970; /*CH */
    YOW += y[23]/26.037940; /*C2H2 */
    YOW += y[24]/27.045910; /*C2H3 */
    YOW += y[25]/43.045310; /*CH2CHO */
    YOW += y[26]/44.053280; /*C2H4O */
    YOW += y[27]/42.037340; /*CH2CO */
    YOW += y[28]/41.029370; /*HCCO */
    YOW += y[29]/25.029970; /*C2H */
    YOW += y[30]/31.034310; /*CH2OH */
    YOW += y[31]/32.042280; /*CH3OH */
    YOW += y[32]/40.064880; /*C3H4 */
    YOW += y[33]/39.056910; /*C3H3 */
    YOW += y[34]/41.072850; /*C3H5 */
    YOW += y[35]/42.080820; /*C3H6 */
    YOW += y[36]/44.096760; /*C3H8 */
    YOW += y[37]/43.088790; /*I-C3H7 */
    YOW += y[38]/43.088790; /*N-C3H7 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(28.013400*YOW); 
    x[1] = y[1]/(39.948000*YOW); 
    x[2] = y[2]/(1.007970*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(2.015940*YOW); 
    x[7] = y[7]/(18.015340*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(34.014740*YOW); 
    x[10] = y[10]/(28.010400*YOW); 
    x[11] = y[11]/(44.009800*YOW); 
    x[12] = y[12]/(29.018370*YOW); 
    x[13] = y[13]/(30.026340*YOW); 
    x[14] = y[14]/(16.042880*YOW); 
    x[15] = y[15]/(15.034910*YOW); 
    x[16] = y[16]/(14.026940*YOW); 
    x[17] = y[17]/(14.026940*YOW); 
    x[18] = y[18]/(28.053880*YOW); 
    x[19] = y[19]/(31.034310*YOW); 
    x[20] = y[20]/(29.061850*YOW); 
    x[21] = y[21]/(30.069820*YOW); 
    x[22] = y[22]/(13.018970*YOW); 
    x[23] = y[23]/(26.037940*YOW); 
    x[24] = y[24]/(27.045910*YOW); 
    x[25] = y[25]/(43.045310*YOW); 
    x[26] = y[26]/(44.053280*YOW); 
    x[27] = y[27]/(42.037340*YOW); 
    x[28] = y[28]/(41.029370*YOW); 
    x[29] = y[29]/(25.029970*YOW); 
    x[30] = y[30]/(31.034310*YOW); 
    x[31] = y[31]/(32.042280*YOW); 
    x[32] = y[32]/(40.064880*YOW); 
    x[33] = y[33]/(39.056910*YOW); 
    x[34] = y[34]/(41.072850*YOW); 
    x[35] = y[35]/(42.080820*YOW); 
    x[36] = y[36]/(44.096760*YOW); 
    x[37] = y[37]/(43.088790*YOW); 
    x[38] = y[38]/(43.088790*YOW); 
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
    result += x[34]*(sor[34]-log((x[34]+1e-100))-logPratio);
    result += x[35]*(sor[35]-log((x[35]+1e-100))-logPratio);
    result += x[36]*(sor[36]-log((x[36]+1e-100))-logPratio);
    result += x[37]*(sor[37]-log((x[37]+1e-100))-logPratio);
    result += x[38]*(sor[38]-log((x[38]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.314e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void CKGBML(double *P, double *T, double *x, int * iwrk, double * rwrk, double * gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    double gort[39]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 39; ++id) {
        result += x[id]*(gort[id]+log((x[id]+1e-100))+logPratio);
    }

    *gbml = result * RT;
}


/*Returns mixture gibbs free energy in mass units */
void CKGBMS(double *P, double *T, double *y, int * iwrk, double * rwrk, double * gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    double gort[39]; /* temporary storage */
    double x[39]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/28.013400; /*N2 */
    YOW += y[1]/39.948000; /*AR */
    YOW += y[2]/1.007970; /*H */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/15.999400; /*O */
    YOW += y[6]/2.015940; /*H2 */
    YOW += y[7]/18.015340; /*H2O */
    YOW += y[8]/33.006770; /*HO2 */
    YOW += y[9]/34.014740; /*H2O2 */
    YOW += y[10]/28.010400; /*CO */
    YOW += y[11]/44.009800; /*CO2 */
    YOW += y[12]/29.018370; /*HCO */
    YOW += y[13]/30.026340; /*CH2O */
    YOW += y[14]/16.042880; /*CH4 */
    YOW += y[15]/15.034910; /*CH3 */
    YOW += y[16]/14.026940; /*T-CH2 */
    YOW += y[17]/14.026940; /*S-CH2 */
    YOW += y[18]/28.053880; /*C2H4 */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/29.061850; /*C2H5 */
    YOW += y[21]/30.069820; /*C2H6 */
    YOW += y[22]/13.018970; /*CH */
    YOW += y[23]/26.037940; /*C2H2 */
    YOW += y[24]/27.045910; /*C2H3 */
    YOW += y[25]/43.045310; /*CH2CHO */
    YOW += y[26]/44.053280; /*C2H4O */
    YOW += y[27]/42.037340; /*CH2CO */
    YOW += y[28]/41.029370; /*HCCO */
    YOW += y[29]/25.029970; /*C2H */
    YOW += y[30]/31.034310; /*CH2OH */
    YOW += y[31]/32.042280; /*CH3OH */
    YOW += y[32]/40.064880; /*C3H4 */
    YOW += y[33]/39.056910; /*C3H3 */
    YOW += y[34]/41.072850; /*C3H5 */
    YOW += y[35]/42.080820; /*C3H6 */
    YOW += y[36]/44.096760; /*C3H8 */
    YOW += y[37]/43.088790; /*I-C3H7 */
    YOW += y[38]/43.088790; /*N-C3H7 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(28.013400*YOW); 
    x[1] = y[1]/(39.948000*YOW); 
    x[2] = y[2]/(1.007970*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(2.015940*YOW); 
    x[7] = y[7]/(18.015340*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(34.014740*YOW); 
    x[10] = y[10]/(28.010400*YOW); 
    x[11] = y[11]/(44.009800*YOW); 
    x[12] = y[12]/(29.018370*YOW); 
    x[13] = y[13]/(30.026340*YOW); 
    x[14] = y[14]/(16.042880*YOW); 
    x[15] = y[15]/(15.034910*YOW); 
    x[16] = y[16]/(14.026940*YOW); 
    x[17] = y[17]/(14.026940*YOW); 
    x[18] = y[18]/(28.053880*YOW); 
    x[19] = y[19]/(31.034310*YOW); 
    x[20] = y[20]/(29.061850*YOW); 
    x[21] = y[21]/(30.069820*YOW); 
    x[22] = y[22]/(13.018970*YOW); 
    x[23] = y[23]/(26.037940*YOW); 
    x[24] = y[24]/(27.045910*YOW); 
    x[25] = y[25]/(43.045310*YOW); 
    x[26] = y[26]/(44.053280*YOW); 
    x[27] = y[27]/(42.037340*YOW); 
    x[28] = y[28]/(41.029370*YOW); 
    x[29] = y[29]/(25.029970*YOW); 
    x[30] = y[30]/(31.034310*YOW); 
    x[31] = y[31]/(32.042280*YOW); 
    x[32] = y[32]/(40.064880*YOW); 
    x[33] = y[33]/(39.056910*YOW); 
    x[34] = y[34]/(41.072850*YOW); 
    x[35] = y[35]/(42.080820*YOW); 
    x[36] = y[36]/(44.096760*YOW); 
    x[37] = y[37]/(43.088790*YOW); 
    x[38] = y[38]/(43.088790*YOW); 
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
    result += x[19]*(gort[19]+log((x[19]+1e-100))+logPratio);
    result += x[20]*(gort[20]+log((x[20]+1e-100))+logPratio);
    result += x[21]*(gort[21]+log((x[21]+1e-100))+logPratio);
    result += x[22]*(gort[22]+log((x[22]+1e-100))+logPratio);
    result += x[23]*(gort[23]+log((x[23]+1e-100))+logPratio);
    result += x[24]*(gort[24]+log((x[24]+1e-100))+logPratio);
    result += x[25]*(gort[25]+log((x[25]+1e-100))+logPratio);
    result += x[26]*(gort[26]+log((x[26]+1e-100))+logPratio);
    result += x[27]*(gort[27]+log((x[27]+1e-100))+logPratio);
    result += x[28]*(gort[28]+log((x[28]+1e-100))+logPratio);
    result += x[29]*(gort[29]+log((x[29]+1e-100))+logPratio);
    result += x[30]*(gort[30]+log((x[30]+1e-100))+logPratio);
    result += x[31]*(gort[31]+log((x[31]+1e-100))+logPratio);
    result += x[32]*(gort[32]+log((x[32]+1e-100))+logPratio);
    result += x[33]*(gort[33]+log((x[33]+1e-100))+logPratio);
    result += x[34]*(gort[34]+log((x[34]+1e-100))+logPratio);
    result += x[35]*(gort[35]+log((x[35]+1e-100))+logPratio);
    result += x[36]*(gort[36]+log((x[36]+1e-100))+logPratio);
    result += x[37]*(gort[37]+log((x[37]+1e-100))+logPratio);
    result += x[38]*(gort[38]+log((x[38]+1e-100))+logPratio);
    /*Scale by RT/W */
    *gbms = result * RT * YOW;
}


/*Returns mean helmholtz free energy in molar units */
void CKABML(double *P, double *T, double *x, int * iwrk, double * rwrk, double * abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    double aort[39]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 39; ++id) {
        result += x[id]*(aort[id]+log((x[id]+1e-100))+logPratio);
    }

    *abml = result * RT;
}


/*Returns mixture helmholtz free energy in mass units */
void CKABMS(double *P, double *T, double *y, int * iwrk, double * rwrk, double * abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    double aort[39]; /* temporary storage */
    double x[39]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/28.013400; /*N2 */
    YOW += y[1]/39.948000; /*AR */
    YOW += y[2]/1.007970; /*H */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/15.999400; /*O */
    YOW += y[6]/2.015940; /*H2 */
    YOW += y[7]/18.015340; /*H2O */
    YOW += y[8]/33.006770; /*HO2 */
    YOW += y[9]/34.014740; /*H2O2 */
    YOW += y[10]/28.010400; /*CO */
    YOW += y[11]/44.009800; /*CO2 */
    YOW += y[12]/29.018370; /*HCO */
    YOW += y[13]/30.026340; /*CH2O */
    YOW += y[14]/16.042880; /*CH4 */
    YOW += y[15]/15.034910; /*CH3 */
    YOW += y[16]/14.026940; /*T-CH2 */
    YOW += y[17]/14.026940; /*S-CH2 */
    YOW += y[18]/28.053880; /*C2H4 */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/29.061850; /*C2H5 */
    YOW += y[21]/30.069820; /*C2H6 */
    YOW += y[22]/13.018970; /*CH */
    YOW += y[23]/26.037940; /*C2H2 */
    YOW += y[24]/27.045910; /*C2H3 */
    YOW += y[25]/43.045310; /*CH2CHO */
    YOW += y[26]/44.053280; /*C2H4O */
    YOW += y[27]/42.037340; /*CH2CO */
    YOW += y[28]/41.029370; /*HCCO */
    YOW += y[29]/25.029970; /*C2H */
    YOW += y[30]/31.034310; /*CH2OH */
    YOW += y[31]/32.042280; /*CH3OH */
    YOW += y[32]/40.064880; /*C3H4 */
    YOW += y[33]/39.056910; /*C3H3 */
    YOW += y[34]/41.072850; /*C3H5 */
    YOW += y[35]/42.080820; /*C3H6 */
    YOW += y[36]/44.096760; /*C3H8 */
    YOW += y[37]/43.088790; /*I-C3H7 */
    YOW += y[38]/43.088790; /*N-C3H7 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(28.013400*YOW); 
    x[1] = y[1]/(39.948000*YOW); 
    x[2] = y[2]/(1.007970*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(2.015940*YOW); 
    x[7] = y[7]/(18.015340*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(34.014740*YOW); 
    x[10] = y[10]/(28.010400*YOW); 
    x[11] = y[11]/(44.009800*YOW); 
    x[12] = y[12]/(29.018370*YOW); 
    x[13] = y[13]/(30.026340*YOW); 
    x[14] = y[14]/(16.042880*YOW); 
    x[15] = y[15]/(15.034910*YOW); 
    x[16] = y[16]/(14.026940*YOW); 
    x[17] = y[17]/(14.026940*YOW); 
    x[18] = y[18]/(28.053880*YOW); 
    x[19] = y[19]/(31.034310*YOW); 
    x[20] = y[20]/(29.061850*YOW); 
    x[21] = y[21]/(30.069820*YOW); 
    x[22] = y[22]/(13.018970*YOW); 
    x[23] = y[23]/(26.037940*YOW); 
    x[24] = y[24]/(27.045910*YOW); 
    x[25] = y[25]/(43.045310*YOW); 
    x[26] = y[26]/(44.053280*YOW); 
    x[27] = y[27]/(42.037340*YOW); 
    x[28] = y[28]/(41.029370*YOW); 
    x[29] = y[29]/(25.029970*YOW); 
    x[30] = y[30]/(31.034310*YOW); 
    x[31] = y[31]/(32.042280*YOW); 
    x[32] = y[32]/(40.064880*YOW); 
    x[33] = y[33]/(39.056910*YOW); 
    x[34] = y[34]/(41.072850*YOW); 
    x[35] = y[35]/(42.080820*YOW); 
    x[36] = y[36]/(44.096760*YOW); 
    x[37] = y[37]/(43.088790*YOW); 
    x[38] = y[38]/(43.088790*YOW); 
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
    result += x[19]*(aort[19]+log((x[19]+1e-100))+logPratio);
    result += x[20]*(aort[20]+log((x[20]+1e-100))+logPratio);
    result += x[21]*(aort[21]+log((x[21]+1e-100))+logPratio);
    result += x[22]*(aort[22]+log((x[22]+1e-100))+logPratio);
    result += x[23]*(aort[23]+log((x[23]+1e-100))+logPratio);
    result += x[24]*(aort[24]+log((x[24]+1e-100))+logPratio);
    result += x[25]*(aort[25]+log((x[25]+1e-100))+logPratio);
    result += x[26]*(aort[26]+log((x[26]+1e-100))+logPratio);
    result += x[27]*(aort[27]+log((x[27]+1e-100))+logPratio);
    result += x[28]*(aort[28]+log((x[28]+1e-100))+logPratio);
    result += x[29]*(aort[29]+log((x[29]+1e-100))+logPratio);
    result += x[30]*(aort[30]+log((x[30]+1e-100))+logPratio);
    result += x[31]*(aort[31]+log((x[31]+1e-100))+logPratio);
    result += x[32]*(aort[32]+log((x[32]+1e-100))+logPratio);
    result += x[33]*(aort[33]+log((x[33]+1e-100))+logPratio);
    result += x[34]*(aort[34]+log((x[34]+1e-100))+logPratio);
    result += x[35]*(aort[35]+log((x[35]+1e-100))+logPratio);
    result += x[36]*(aort[36]+log((x[36]+1e-100))+logPratio);
    result += x[37]*(aort[37]+log((x[37]+1e-100))+logPratio);
    result += x[38]*(aort[38]+log((x[38]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * T, double * C, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 39; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[39]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/28.013400; /*N2 */
    YOW += y[1]/39.948000; /*AR */
    YOW += y[2]/1.007970; /*H */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/15.999400; /*O */
    YOW += y[6]/2.015940; /*H2 */
    YOW += y[7]/18.015340; /*H2O */
    YOW += y[8]/33.006770; /*HO2 */
    YOW += y[9]/34.014740; /*H2O2 */
    YOW += y[10]/28.010400; /*CO */
    YOW += y[11]/44.009800; /*CO2 */
    YOW += y[12]/29.018370; /*HCO */
    YOW += y[13]/30.026340; /*CH2O */
    YOW += y[14]/16.042880; /*CH4 */
    YOW += y[15]/15.034910; /*CH3 */
    YOW += y[16]/14.026940; /*T-CH2 */
    YOW += y[17]/14.026940; /*S-CH2 */
    YOW += y[18]/28.053880; /*C2H4 */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/29.061850; /*C2H5 */
    YOW += y[21]/30.069820; /*C2H6 */
    YOW += y[22]/13.018970; /*CH */
    YOW += y[23]/26.037940; /*C2H2 */
    YOW += y[24]/27.045910; /*C2H3 */
    YOW += y[25]/43.045310; /*CH2CHO */
    YOW += y[26]/44.053280; /*C2H4O */
    YOW += y[27]/42.037340; /*CH2CO */
    YOW += y[28]/41.029370; /*HCCO */
    YOW += y[29]/25.029970; /*C2H */
    YOW += y[30]/31.034310; /*CH2OH */
    YOW += y[31]/32.042280; /*CH3OH */
    YOW += y[32]/40.064880; /*C3H4 */
    YOW += y[33]/39.056910; /*C3H3 */
    YOW += y[34]/41.072850; /*C3H5 */
    YOW += y[35]/42.080820; /*C3H6 */
    YOW += y[36]/44.096760; /*C3H8 */
    YOW += y[37]/43.088790; /*I-C3H7 */
    YOW += y[38]/43.088790; /*N-C3H7 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.314e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]/28.013400; 
    c[1] = PWORT * y[1]/39.948000; 
    c[2] = PWORT * y[2]/1.007970; 
    c[3] = PWORT * y[3]/31.998800; 
    c[4] = PWORT * y[4]/17.007370; 
    c[5] = PWORT * y[5]/15.999400; 
    c[6] = PWORT * y[6]/2.015940; 
    c[7] = PWORT * y[7]/18.015340; 
    c[8] = PWORT * y[8]/33.006770; 
    c[9] = PWORT * y[9]/34.014740; 
    c[10] = PWORT * y[10]/28.010400; 
    c[11] = PWORT * y[11]/44.009800; 
    c[12] = PWORT * y[12]/29.018370; 
    c[13] = PWORT * y[13]/30.026340; 
    c[14] = PWORT * y[14]/16.042880; 
    c[15] = PWORT * y[15]/15.034910; 
    c[16] = PWORT * y[16]/14.026940; 
    c[17] = PWORT * y[17]/14.026940; 
    c[18] = PWORT * y[18]/28.053880; 
    c[19] = PWORT * y[19]/31.034310; 
    c[20] = PWORT * y[20]/29.061850; 
    c[21] = PWORT * y[21]/30.069820; 
    c[22] = PWORT * y[22]/13.018970; 
    c[23] = PWORT * y[23]/26.037940; 
    c[24] = PWORT * y[24]/27.045910; 
    c[25] = PWORT * y[25]/43.045310; 
    c[26] = PWORT * y[26]/44.053280; 
    c[27] = PWORT * y[27]/42.037340; 
    c[28] = PWORT * y[28]/41.029370; 
    c[29] = PWORT * y[29]/25.029970; 
    c[30] = PWORT * y[30]/31.034310; 
    c[31] = PWORT * y[31]/32.042280; 
    c[32] = PWORT * y[32]/40.064880; 
    c[33] = PWORT * y[33]/39.056910; 
    c[34] = PWORT * y[34]/41.072850; 
    c[35] = PWORT * y[35]/42.080820; 
    c[36] = PWORT * y[36]/44.096760; 
    c[37] = PWORT * y[37]/43.088790; 
    c[38] = PWORT * y[38]/43.088790; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[39]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 39; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[39]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/28.013400; 
    c[1] = 1e6 * (*rho) * y[1]/39.948000; 
    c[2] = 1e6 * (*rho) * y[2]/1.007970; 
    c[3] = 1e6 * (*rho) * y[3]/31.998800; 
    c[4] = 1e6 * (*rho) * y[4]/17.007370; 
    c[5] = 1e6 * (*rho) * y[5]/15.999400; 
    c[6] = 1e6 * (*rho) * y[6]/2.015940; 
    c[7] = 1e6 * (*rho) * y[7]/18.015340; 
    c[8] = 1e6 * (*rho) * y[8]/33.006770; 
    c[9] = 1e6 * (*rho) * y[9]/34.014740; 
    c[10] = 1e6 * (*rho) * y[10]/28.010400; 
    c[11] = 1e6 * (*rho) * y[11]/44.009800; 
    c[12] = 1e6 * (*rho) * y[12]/29.018370; 
    c[13] = 1e6 * (*rho) * y[13]/30.026340; 
    c[14] = 1e6 * (*rho) * y[14]/16.042880; 
    c[15] = 1e6 * (*rho) * y[15]/15.034910; 
    c[16] = 1e6 * (*rho) * y[16]/14.026940; 
    c[17] = 1e6 * (*rho) * y[17]/14.026940; 
    c[18] = 1e6 * (*rho) * y[18]/28.053880; 
    c[19] = 1e6 * (*rho) * y[19]/31.034310; 
    c[20] = 1e6 * (*rho) * y[20]/29.061850; 
    c[21] = 1e6 * (*rho) * y[21]/30.069820; 
    c[22] = 1e6 * (*rho) * y[22]/13.018970; 
    c[23] = 1e6 * (*rho) * y[23]/26.037940; 
    c[24] = 1e6 * (*rho) * y[24]/27.045910; 
    c[25] = 1e6 * (*rho) * y[25]/43.045310; 
    c[26] = 1e6 * (*rho) * y[26]/44.053280; 
    c[27] = 1e6 * (*rho) * y[27]/42.037340; 
    c[28] = 1e6 * (*rho) * y[28]/41.029370; 
    c[29] = 1e6 * (*rho) * y[29]/25.029970; 
    c[30] = 1e6 * (*rho) * y[30]/31.034310; 
    c[31] = 1e6 * (*rho) * y[31]/32.042280; 
    c[32] = 1e6 * (*rho) * y[32]/40.064880; 
    c[33] = 1e6 * (*rho) * y[33]/39.056910; 
    c[34] = 1e6 * (*rho) * y[34]/41.072850; 
    c[35] = 1e6 * (*rho) * y[35]/42.080820; 
    c[36] = 1e6 * (*rho) * y[36]/44.096760; 
    c[37] = 1e6 * (*rho) * y[37]/43.088790; 
    c[38] = 1e6 * (*rho) * y[38]/43.088790; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[39]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*39.948000; /*AR */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*2.015940; /*H2 */
    XW += x[7]*18.015340; /*H2O */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*34.014740; /*H2O2 */
    XW += x[10]*28.010400; /*CO */
    XW += x[11]*44.009800; /*CO2 */
    XW += x[12]*29.018370; /*HCO */
    XW += x[13]*30.026340; /*CH2O */
    XW += x[14]*16.042880; /*CH4 */
    XW += x[15]*15.034910; /*CH3 */
    XW += x[16]*14.026940; /*T-CH2 */
    XW += x[17]*14.026940; /*S-CH2 */
    XW += x[18]*28.053880; /*C2H4 */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*29.061850; /*C2H5 */
    XW += x[21]*30.069820; /*C2H6 */
    XW += x[22]*13.018970; /*CH */
    XW += x[23]*26.037940; /*C2H2 */
    XW += x[24]*27.045910; /*C2H3 */
    XW += x[25]*43.045310; /*CH2CHO */
    XW += x[26]*44.053280; /*C2H4O */
    XW += x[27]*42.037340; /*CH2CO */
    XW += x[28]*41.029370; /*HCCO */
    XW += x[29]*25.029970; /*C2H */
    XW += x[30]*31.034310; /*CH2OH */
    XW += x[31]*32.042280; /*CH3OH */
    XW += x[32]*40.064880; /*C3H4 */
    XW += x[33]*39.056910; /*C3H3 */
    XW += x[34]*41.072850; /*C3H5 */
    XW += x[35]*42.080820; /*C3H6 */
    XW += x[36]*44.096760; /*C3H8 */
    XW += x[37]*43.088790; /*I-C3H7 */
    XW += x[38]*43.088790; /*N-C3H7 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 39; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * T, double * C, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 39; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 175; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[39]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/28.013400; /*N2 */
    YOW += y[1]/39.948000; /*AR */
    YOW += y[2]/1.007970; /*H */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/15.999400; /*O */
    YOW += y[6]/2.015940; /*H2 */
    YOW += y[7]/18.015340; /*H2O */
    YOW += y[8]/33.006770; /*HO2 */
    YOW += y[9]/34.014740; /*H2O2 */
    YOW += y[10]/28.010400; /*CO */
    YOW += y[11]/44.009800; /*CO2 */
    YOW += y[12]/29.018370; /*HCO */
    YOW += y[13]/30.026340; /*CH2O */
    YOW += y[14]/16.042880; /*CH4 */
    YOW += y[15]/15.034910; /*CH3 */
    YOW += y[16]/14.026940; /*T-CH2 */
    YOW += y[17]/14.026940; /*S-CH2 */
    YOW += y[18]/28.053880; /*C2H4 */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/29.061850; /*C2H5 */
    YOW += y[21]/30.069820; /*C2H6 */
    YOW += y[22]/13.018970; /*CH */
    YOW += y[23]/26.037940; /*C2H2 */
    YOW += y[24]/27.045910; /*C2H3 */
    YOW += y[25]/43.045310; /*CH2CHO */
    YOW += y[26]/44.053280; /*C2H4O */
    YOW += y[27]/42.037340; /*CH2CO */
    YOW += y[28]/41.029370; /*HCCO */
    YOW += y[29]/25.029970; /*C2H */
    YOW += y[30]/31.034310; /*CH2OH */
    YOW += y[31]/32.042280; /*CH3OH */
    YOW += y[32]/40.064880; /*C3H4 */
    YOW += y[33]/39.056910; /*C3H3 */
    YOW += y[34]/41.072850; /*C3H5 */
    YOW += y[35]/42.080820; /*C3H6 */
    YOW += y[36]/44.096760; /*C3H8 */
    YOW += y[37]/43.088790; /*I-C3H7 */
    YOW += y[38]/43.088790; /*N-C3H7 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.314e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]/28.013400; 
    c[1] = PWORT * y[1]/39.948000; 
    c[2] = PWORT * y[2]/1.007970; 
    c[3] = PWORT * y[3]/31.998800; 
    c[4] = PWORT * y[4]/17.007370; 
    c[5] = PWORT * y[5]/15.999400; 
    c[6] = PWORT * y[6]/2.015940; 
    c[7] = PWORT * y[7]/18.015340; 
    c[8] = PWORT * y[8]/33.006770; 
    c[9] = PWORT * y[9]/34.014740; 
    c[10] = PWORT * y[10]/28.010400; 
    c[11] = PWORT * y[11]/44.009800; 
    c[12] = PWORT * y[12]/29.018370; 
    c[13] = PWORT * y[13]/30.026340; 
    c[14] = PWORT * y[14]/16.042880; 
    c[15] = PWORT * y[15]/15.034910; 
    c[16] = PWORT * y[16]/14.026940; 
    c[17] = PWORT * y[17]/14.026940; 
    c[18] = PWORT * y[18]/28.053880; 
    c[19] = PWORT * y[19]/31.034310; 
    c[20] = PWORT * y[20]/29.061850; 
    c[21] = PWORT * y[21]/30.069820; 
    c[22] = PWORT * y[22]/13.018970; 
    c[23] = PWORT * y[23]/26.037940; 
    c[24] = PWORT * y[24]/27.045910; 
    c[25] = PWORT * y[25]/43.045310; 
    c[26] = PWORT * y[26]/44.053280; 
    c[27] = PWORT * y[27]/42.037340; 
    c[28] = PWORT * y[28]/41.029370; 
    c[29] = PWORT * y[29]/25.029970; 
    c[30] = PWORT * y[30]/31.034310; 
    c[31] = PWORT * y[31]/32.042280; 
    c[32] = PWORT * y[32]/40.064880; 
    c[33] = PWORT * y[33]/39.056910; 
    c[34] = PWORT * y[34]/41.072850; 
    c[35] = PWORT * y[35]/42.080820; 
    c[36] = PWORT * y[36]/44.096760; 
    c[37] = PWORT * y[37]/43.088790; 
    c[38] = PWORT * y[38]/43.088790; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 175; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[39]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 39; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 175; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[39]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/28.013400; 
    c[1] = 1e6 * (*rho) * y[1]/39.948000; 
    c[2] = 1e6 * (*rho) * y[2]/1.007970; 
    c[3] = 1e6 * (*rho) * y[3]/31.998800; 
    c[4] = 1e6 * (*rho) * y[4]/17.007370; 
    c[5] = 1e6 * (*rho) * y[5]/15.999400; 
    c[6] = 1e6 * (*rho) * y[6]/2.015940; 
    c[7] = 1e6 * (*rho) * y[7]/18.015340; 
    c[8] = 1e6 * (*rho) * y[8]/33.006770; 
    c[9] = 1e6 * (*rho) * y[9]/34.014740; 
    c[10] = 1e6 * (*rho) * y[10]/28.010400; 
    c[11] = 1e6 * (*rho) * y[11]/44.009800; 
    c[12] = 1e6 * (*rho) * y[12]/29.018370; 
    c[13] = 1e6 * (*rho) * y[13]/30.026340; 
    c[14] = 1e6 * (*rho) * y[14]/16.042880; 
    c[15] = 1e6 * (*rho) * y[15]/15.034910; 
    c[16] = 1e6 * (*rho) * y[16]/14.026940; 
    c[17] = 1e6 * (*rho) * y[17]/14.026940; 
    c[18] = 1e6 * (*rho) * y[18]/28.053880; 
    c[19] = 1e6 * (*rho) * y[19]/31.034310; 
    c[20] = 1e6 * (*rho) * y[20]/29.061850; 
    c[21] = 1e6 * (*rho) * y[21]/30.069820; 
    c[22] = 1e6 * (*rho) * y[22]/13.018970; 
    c[23] = 1e6 * (*rho) * y[23]/26.037940; 
    c[24] = 1e6 * (*rho) * y[24]/27.045910; 
    c[25] = 1e6 * (*rho) * y[25]/43.045310; 
    c[26] = 1e6 * (*rho) * y[26]/44.053280; 
    c[27] = 1e6 * (*rho) * y[27]/42.037340; 
    c[28] = 1e6 * (*rho) * y[28]/41.029370; 
    c[29] = 1e6 * (*rho) * y[29]/25.029970; 
    c[30] = 1e6 * (*rho) * y[30]/31.034310; 
    c[31] = 1e6 * (*rho) * y[31]/32.042280; 
    c[32] = 1e6 * (*rho) * y[32]/40.064880; 
    c[33] = 1e6 * (*rho) * y[33]/39.056910; 
    c[34] = 1e6 * (*rho) * y[34]/41.072850; 
    c[35] = 1e6 * (*rho) * y[35]/42.080820; 
    c[36] = 1e6 * (*rho) * y[36]/44.096760; 
    c[37] = 1e6 * (*rho) * y[37]/43.088790; 
    c[38] = 1e6 * (*rho) * y[38]/43.088790; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 175; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[39]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*39.948000; /*AR */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*2.015940; /*H2 */
    XW += x[7]*18.015340; /*H2O */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*34.014740; /*H2O2 */
    XW += x[10]*28.010400; /*CO */
    XW += x[11]*44.009800; /*CO2 */
    XW += x[12]*29.018370; /*HCO */
    XW += x[13]*30.026340; /*CH2O */
    XW += x[14]*16.042880; /*CH4 */
    XW += x[15]*15.034910; /*CH3 */
    XW += x[16]*14.026940; /*T-CH2 */
    XW += x[17]*14.026940; /*S-CH2 */
    XW += x[18]*28.053880; /*C2H4 */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*29.061850; /*C2H5 */
    XW += x[21]*30.069820; /*C2H6 */
    XW += x[22]*13.018970; /*CH */
    XW += x[23]*26.037940; /*C2H2 */
    XW += x[24]*27.045910; /*C2H3 */
    XW += x[25]*43.045310; /*CH2CHO */
    XW += x[26]*44.053280; /*C2H4O */
    XW += x[27]*42.037340; /*CH2CO */
    XW += x[28]*41.029370; /*HCCO */
    XW += x[29]*25.029970; /*C2H */
    XW += x[30]*31.034310; /*CH2OH */
    XW += x[31]*32.042280; /*CH3OH */
    XW += x[32]*40.064880; /*C3H4 */
    XW += x[33]*39.056910; /*C3H3 */
    XW += x[34]*41.072850; /*C3H5 */
    XW += x[35]*42.080820; /*C3H6 */
    XW += x[36]*44.096760; /*C3H8 */
    XW += x[37]*43.088790; /*I-C3H7 */
    XW += x[38]*43.088790; /*N-C3H7 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 39; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 175; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim, int * iwrk, double * rwrk, int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 39 * 175; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + O2 <=> OH + O */
    nuki[ 2 * kd + 0 ] = -1 ;
    nuki[ 3 * kd + 0 ] = -1 ;
    nuki[ 4 * kd + 0 ] = +1 ;
    nuki[ 5 * kd + 0 ] = +1 ;

    /*reaction 2: H2 + O <=> OH + H */
    nuki[ 6 * kd + 1 ] = -1 ;
    nuki[ 5 * kd + 1 ] = -1 ;
    nuki[ 4 * kd + 1 ] = +1 ;
    nuki[ 2 * kd + 1 ] = +1 ;

    /*reaction 3: H2 + OH <=> H2O + H */
    nuki[ 6 * kd + 2 ] = -1 ;
    nuki[ 4 * kd + 2 ] = -1 ;
    nuki[ 7 * kd + 2 ] = +1 ;
    nuki[ 2 * kd + 2 ] = +1 ;

    /*reaction 4: H2O + O <=> 2 OH */
    nuki[ 7 * kd + 3 ] = -1 ;
    nuki[ 5 * kd + 3 ] = -1 ;
    nuki[ 4 * kd + 3 ] = +2 ;

    /*reaction 5: H + O + M <=> OH + M */
    nuki[ 2 * kd + 4 ] = -1 ;
    nuki[ 5 * kd + 4 ] = -1 ;
    nuki[ 4 * kd + 4 ] = +1 ;

    /*reaction 6: H2 + O2 <=> 2 OH */
    nuki[ 6 * kd + 5 ] = -1 ;
    nuki[ 3 * kd + 5 ] = -1 ;
    nuki[ 4 * kd + 5 ] = +2 ;

    /*reaction 7: 2 H + M <=> H2 + M */
    nuki[ 2 * kd + 6 ] = -2 ;
    nuki[ 6 * kd + 6 ] = +1 ;

    /*reaction 8: H + OH + M <=> H2O + M */
    nuki[ 2 * kd + 7 ] = -1 ;
    nuki[ 4 * kd + 7 ] = -1 ;
    nuki[ 7 * kd + 7 ] = +1 ;

    /*reaction 9: 2 O + M <=> O2 + M */
    nuki[ 5 * kd + 8 ] = -2 ;
    nuki[ 3 * kd + 8 ] = +1 ;

    /*reaction 10: H + O2 (+M) <=> HO2 (+M) */
    nuki[ 2 * kd + 9 ] = -1 ;
    nuki[ 3 * kd + 9 ] = -1 ;
    nuki[ 8 * kd + 9 ] = +1 ;

    /*reaction 11: O + OH + M <=> HO2 + M */
    nuki[ 5 * kd + 10 ] = -1 ;
    nuki[ 4 * kd + 10 ] = -1 ;
    nuki[ 8 * kd + 10 ] = +1 ;

    /*reaction 12: HO2 + H <=> 2 OH */
    nuki[ 8 * kd + 11 ] = -1 ;
    nuki[ 2 * kd + 11 ] = -1 ;
    nuki[ 4 * kd + 11 ] = +2 ;

    /*reaction 13: HO2 + H <=> H2 + O2 */
    nuki[ 8 * kd + 12 ] = -1 ;
    nuki[ 2 * kd + 12 ] = -1 ;
    nuki[ 6 * kd + 12 ] = +1 ;
    nuki[ 3 * kd + 12 ] = +1 ;

    /*reaction 14: HO2 + H <=> H2O + O */
    nuki[ 8 * kd + 13 ] = -1 ;
    nuki[ 2 * kd + 13 ] = -1 ;
    nuki[ 7 * kd + 13 ] = +1 ;
    nuki[ 5 * kd + 13 ] = +1 ;

    /*reaction 15: HO2 + O <=> OH + O2 */
    nuki[ 8 * kd + 14 ] = -1 ;
    nuki[ 5 * kd + 14 ] = -1 ;
    nuki[ 4 * kd + 14 ] = +1 ;
    nuki[ 3 * kd + 14 ] = +1 ;

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    nuki[ 8 * kd + 15 ] = -1 ;
    nuki[ 4 * kd + 15 ] = -1 ;
    nuki[ 7 * kd + 15 ] = +1 ;
    nuki[ 3 * kd + 15 ] = +1 ;

    /*reaction 17: 2 OH (+M) <=> H2O2 (+M) */
    nuki[ 4 * kd + 16 ] = -2 ;
    nuki[ 9 * kd + 16 ] = +1 ;

    /*reaction 18: 2 HO2 <=> H2O2 + O2 */
    nuki[ 8 * kd + 17 ] = -2 ;
    nuki[ 9 * kd + 17 ] = +1 ;
    nuki[ 3 * kd + 17 ] = +1 ;

    /*reaction 19: H2O2 + H <=> HO2 + H2 */
    nuki[ 9 * kd + 18 ] = -1 ;
    nuki[ 2 * kd + 18 ] = -1 ;
    nuki[ 8 * kd + 18 ] = +1 ;
    nuki[ 6 * kd + 18 ] = +1 ;

    /*reaction 20: H2O2 + H <=> H2O + OH */
    nuki[ 9 * kd + 19 ] = -1 ;
    nuki[ 2 * kd + 19 ] = -1 ;
    nuki[ 7 * kd + 19 ] = +1 ;
    nuki[ 4 * kd + 19 ] = +1 ;

    /*reaction 21: H2O2 + OH <=> H2O + HO2 */
    nuki[ 9 * kd + 20 ] = -1 ;
    nuki[ 4 * kd + 20 ] = -1 ;
    nuki[ 7 * kd + 20 ] = +1 ;
    nuki[ 8 * kd + 20 ] = +1 ;

    /*reaction 22: H2O2 + O <=> HO2 + OH */
    nuki[ 9 * kd + 21 ] = -1 ;
    nuki[ 5 * kd + 21 ] = -1 ;
    nuki[ 8 * kd + 21 ] = +1 ;
    nuki[ 4 * kd + 21 ] = +1 ;

    /*reaction 23: CO + OH <=> CO2 + H */
    nuki[ 10 * kd + 22 ] = -1 ;
    nuki[ 4 * kd + 22 ] = -1 ;
    nuki[ 11 * kd + 22 ] = +1 ;
    nuki[ 2 * kd + 22 ] = +1 ;

    /*reaction 24: CO + HO2 <=> CO2 + OH */
    nuki[ 10 * kd + 23 ] = -1 ;
    nuki[ 8 * kd + 23 ] = -1 ;
    nuki[ 11 * kd + 23 ] = +1 ;
    nuki[ 4 * kd + 23 ] = +1 ;

    /*reaction 25: HCO + M <=> CO + H + M */
    nuki[ 12 * kd + 24 ] = -1 ;
    nuki[ 10 * kd + 24 ] = +1 ;
    nuki[ 2 * kd + 24 ] = +1 ;

    /*reaction 26: HCO + H <=> CO + H2 */
    nuki[ 12 * kd + 25 ] = -1 ;
    nuki[ 2 * kd + 25 ] = -1 ;
    nuki[ 10 * kd + 25 ] = +1 ;
    nuki[ 6 * kd + 25 ] = +1 ;

    /*reaction 27: HCO + O <=> CO + OH */
    nuki[ 12 * kd + 26 ] = -1 ;
    nuki[ 5 * kd + 26 ] = -1 ;
    nuki[ 10 * kd + 26 ] = +1 ;
    nuki[ 4 * kd + 26 ] = +1 ;

    /*reaction 28: HCO + O <=> CO2 + H */
    nuki[ 12 * kd + 27 ] = -1 ;
    nuki[ 5 * kd + 27 ] = -1 ;
    nuki[ 11 * kd + 27 ] = +1 ;
    nuki[ 2 * kd + 27 ] = +1 ;

    /*reaction 29: HCO + OH <=> CO + H2O */
    nuki[ 12 * kd + 28 ] = -1 ;
    nuki[ 4 * kd + 28 ] = -1 ;
    nuki[ 10 * kd + 28 ] = +1 ;
    nuki[ 7 * kd + 28 ] = +1 ;

    /*reaction 30: HCO + O2 <=> CO + HO2 */
    nuki[ 12 * kd + 29 ] = -1 ;
    nuki[ 3 * kd + 29 ] = -1 ;
    nuki[ 10 * kd + 29 ] = +1 ;
    nuki[ 8 * kd + 29 ] = +1 ;

    /*reaction 31: CH2O + M <=> HCO + H + M */
    nuki[ 13 * kd + 30 ] = -1 ;
    nuki[ 12 * kd + 30 ] = +1 ;
    nuki[ 2 * kd + 30 ] = +1 ;

    /*reaction 32: CH2O + H <=> HCO + H2 */
    nuki[ 13 * kd + 31 ] = -1 ;
    nuki[ 2 * kd + 31 ] = -1 ;
    nuki[ 12 * kd + 31 ] = +1 ;
    nuki[ 6 * kd + 31 ] = +1 ;

    /*reaction 33: CH2O + O <=> HCO + OH */
    nuki[ 13 * kd + 32 ] = -1 ;
    nuki[ 5 * kd + 32 ] = -1 ;
    nuki[ 12 * kd + 32 ] = +1 ;
    nuki[ 4 * kd + 32 ] = +1 ;

    /*reaction 34: CH2O + OH <=> HCO + H2O */
    nuki[ 13 * kd + 33 ] = -1 ;
    nuki[ 4 * kd + 33 ] = -1 ;
    nuki[ 12 * kd + 33 ] = +1 ;
    nuki[ 7 * kd + 33 ] = +1 ;

    /*reaction 35: CH4 + H <=> H2 + CH3 */
    nuki[ 14 * kd + 34 ] = -1 ;
    nuki[ 2 * kd + 34 ] = -1 ;
    nuki[ 6 * kd + 34 ] = +1 ;
    nuki[ 15 * kd + 34 ] = +1 ;

    /*reaction 36: CH4 + OH <=> H2O + CH3 */
    nuki[ 14 * kd + 35 ] = -1 ;
    nuki[ 4 * kd + 35 ] = -1 ;
    nuki[ 7 * kd + 35 ] = +1 ;
    nuki[ 15 * kd + 35 ] = +1 ;

    /*reaction 37: CH4 + O <=> CH3 + OH */
    nuki[ 14 * kd + 36 ] = -1 ;
    nuki[ 5 * kd + 36 ] = -1 ;
    nuki[ 15 * kd + 36 ] = +1 ;
    nuki[ 4 * kd + 36 ] = +1 ;

    /*reaction 38: CH4 + O2 <=> CH3 + HO2 */
    nuki[ 14 * kd + 37 ] = -1 ;
    nuki[ 3 * kd + 37 ] = -1 ;
    nuki[ 15 * kd + 37 ] = +1 ;
    nuki[ 8 * kd + 37 ] = +1 ;

    /*reaction 39: CH4 + HO2 <=> CH3 + H2O2 */
    nuki[ 14 * kd + 38 ] = -1 ;
    nuki[ 8 * kd + 38 ] = -1 ;
    nuki[ 15 * kd + 38 ] = +1 ;
    nuki[ 9 * kd + 38 ] = +1 ;

    /*reaction 40: CH3 + H <=> T-CH2 + H2 */
    nuki[ 15 * kd + 39 ] = -1 ;
    nuki[ 2 * kd + 39 ] = -1 ;
    nuki[ 16 * kd + 39 ] = +1 ;
    nuki[ 6 * kd + 39 ] = +1 ;

    /*reaction 41: CH3 + H <=> S-CH2 + H2 */
    nuki[ 15 * kd + 40 ] = -1 ;
    nuki[ 2 * kd + 40 ] = -1 ;
    nuki[ 17 * kd + 40 ] = +1 ;
    nuki[ 6 * kd + 40 ] = +1 ;

    /*reaction 42: CH3 + OH <=> S-CH2 + H2O */
    nuki[ 15 * kd + 41 ] = -1 ;
    nuki[ 4 * kd + 41 ] = -1 ;
    nuki[ 17 * kd + 41 ] = +1 ;
    nuki[ 7 * kd + 41 ] = +1 ;

    /*reaction 43: CH3 + O <=> CH2O + H */
    nuki[ 15 * kd + 42 ] = -1 ;
    nuki[ 5 * kd + 42 ] = -1 ;
    nuki[ 13 * kd + 42 ] = +1 ;
    nuki[ 2 * kd + 42 ] = +1 ;

    /*reaction 44: CH3 + T-CH2 <=> C2H4 + H */
    nuki[ 15 * kd + 43 ] = -1 ;
    nuki[ 16 * kd + 43 ] = -1 ;
    nuki[ 18 * kd + 43 ] = +1 ;
    nuki[ 2 * kd + 43 ] = +1 ;

    /*reaction 45: CH3 + HO2 <=> CH3O + OH */
    nuki[ 15 * kd + 44 ] = -1 ;
    nuki[ 8 * kd + 44 ] = -1 ;
    nuki[ 19 * kd + 44 ] = +1 ;
    nuki[ 4 * kd + 44 ] = +1 ;

    /*reaction 46: CH3 + O2 <=> CH2O + OH */
    nuki[ 15 * kd + 45 ] = -1 ;
    nuki[ 3 * kd + 45 ] = -1 ;
    nuki[ 13 * kd + 45 ] = +1 ;
    nuki[ 4 * kd + 45 ] = +1 ;

    /*reaction 47: CH3 + O2 <=> CH3O + O */
    nuki[ 15 * kd + 46 ] = -1 ;
    nuki[ 3 * kd + 46 ] = -1 ;
    nuki[ 19 * kd + 46 ] = +1 ;
    nuki[ 5 * kd + 46 ] = +1 ;

    /*reaction 48: 2 CH3 <=> C2H4 + H2 */
    nuki[ 15 * kd + 47 ] = -2 ;
    nuki[ 18 * kd + 47 ] = +1 ;
    nuki[ 6 * kd + 47 ] = +1 ;

    /*reaction 49: 2 CH3 <=> C2H5 + H */
    nuki[ 15 * kd + 48 ] = -2 ;
    nuki[ 20 * kd + 48 ] = +1 ;
    nuki[ 2 * kd + 48 ] = +1 ;

    /*reaction 50: CH3 + H (+M) <=> CH4 (+M) */
    nuki[ 15 * kd + 49 ] = -1 ;
    nuki[ 2 * kd + 49 ] = -1 ;
    nuki[ 14 * kd + 49 ] = +1 ;

    /*reaction 51: 2 CH3 (+M) <=> C2H6 (+M) */
    nuki[ 15 * kd + 50 ] = -2 ;
    nuki[ 21 * kd + 50 ] = +1 ;

    /*reaction 52: S-CH2 + OH <=> CH2O + H */
    nuki[ 17 * kd + 51 ] = -1 ;
    nuki[ 4 * kd + 51 ] = -1 ;
    nuki[ 13 * kd + 51 ] = +1 ;
    nuki[ 2 * kd + 51 ] = +1 ;

    /*reaction 53: S-CH2 + O2 <=> CO + OH + H */
    nuki[ 17 * kd + 52 ] = -1 ;
    nuki[ 3 * kd + 52 ] = -1 ;
    nuki[ 10 * kd + 52 ] = +1 ;
    nuki[ 4 * kd + 52 ] = +1 ;
    nuki[ 2 * kd + 52 ] = +1 ;

    /*reaction 54: S-CH2 + CO2 <=> CO + CH2O */
    nuki[ 17 * kd + 53 ] = -1 ;
    nuki[ 11 * kd + 53 ] = -1 ;
    nuki[ 10 * kd + 53 ] = +1 ;
    nuki[ 13 * kd + 53 ] = +1 ;

    /*reaction 55: S-CH2 + M <=> T-CH2 + M */
    nuki[ 17 * kd + 54 ] = -1 ;
    nuki[ 16 * kd + 54 ] = +1 ;

    /*reaction 56: T-CH2 + H <=> CH + H2 */
    nuki[ 16 * kd + 55 ] = -1 ;
    nuki[ 2 * kd + 55 ] = -1 ;
    nuki[ 22 * kd + 55 ] = +1 ;
    nuki[ 6 * kd + 55 ] = +1 ;

    /*reaction 57: T-CH2 + OH <=> CH2O + H */
    nuki[ 16 * kd + 56 ] = -1 ;
    nuki[ 4 * kd + 56 ] = -1 ;
    nuki[ 13 * kd + 56 ] = +1 ;
    nuki[ 2 * kd + 56 ] = +1 ;

    /*reaction 58: T-CH2 + OH <=> CH + H2O */
    nuki[ 16 * kd + 57 ] = -1 ;
    nuki[ 4 * kd + 57 ] = -1 ;
    nuki[ 22 * kd + 57 ] = +1 ;
    nuki[ 7 * kd + 57 ] = +1 ;

    /*reaction 59: T-CH2 + O <=> CO + 2 H */
    nuki[ 16 * kd + 58 ] = -1 ;
    nuki[ 5 * kd + 58 ] = -1 ;
    nuki[ 10 * kd + 58 ] = +1 ;
    nuki[ 2 * kd + 58 ] = +2 ;

    /*reaction 60: T-CH2 + O <=> CO + H2 */
    nuki[ 16 * kd + 59 ] = -1 ;
    nuki[ 5 * kd + 59 ] = -1 ;
    nuki[ 10 * kd + 59 ] = +1 ;
    nuki[ 6 * kd + 59 ] = +1 ;

    /*reaction 61: T-CH2 + O2 <=> CO2 + H2 */
    nuki[ 16 * kd + 60 ] = -1 ;
    nuki[ 3 * kd + 60 ] = -1 ;
    nuki[ 11 * kd + 60 ] = +1 ;
    nuki[ 6 * kd + 60 ] = +1 ;

    /*reaction 62: T-CH2 + O2 <=> CO + OH + H */
    nuki[ 16 * kd + 61 ] = -1 ;
    nuki[ 3 * kd + 61 ] = -1 ;
    nuki[ 10 * kd + 61 ] = +1 ;
    nuki[ 4 * kd + 61 ] = +1 ;
    nuki[ 2 * kd + 61 ] = +1 ;

    /*reaction 63: 2 T-CH2 <=> C2H2 + 2 H */
    nuki[ 16 * kd + 62 ] = -2 ;
    nuki[ 23 * kd + 62 ] = +1 ;
    nuki[ 2 * kd + 62 ] = +2 ;

    /*reaction 64: CH + O <=> CO + H */
    nuki[ 22 * kd + 63 ] = -1 ;
    nuki[ 5 * kd + 63 ] = -1 ;
    nuki[ 10 * kd + 63 ] = +1 ;
    nuki[ 2 * kd + 63 ] = +1 ;

    /*reaction 65: CH + O2 <=> HCO + O */
    nuki[ 22 * kd + 64 ] = -1 ;
    nuki[ 3 * kd + 64 ] = -1 ;
    nuki[ 12 * kd + 64 ] = +1 ;
    nuki[ 5 * kd + 64 ] = +1 ;

    /*reaction 66: CH + H2O <=> CH2O + H */
    nuki[ 22 * kd + 65 ] = -1 ;
    nuki[ 7 * kd + 65 ] = -1 ;
    nuki[ 13 * kd + 65 ] = +1 ;
    nuki[ 2 * kd + 65 ] = +1 ;

    /*reaction 67: CH + CO2 <=> HCO + CO */
    nuki[ 22 * kd + 66 ] = -1 ;
    nuki[ 11 * kd + 66 ] = -1 ;
    nuki[ 12 * kd + 66 ] = +1 ;
    nuki[ 10 * kd + 66 ] = +1 ;

    /*reaction 68: CH3O + H <=> CH2O + H2 */
    nuki[ 19 * kd + 67 ] = -1 ;
    nuki[ 2 * kd + 67 ] = -1 ;
    nuki[ 13 * kd + 67 ] = +1 ;
    nuki[ 6 * kd + 67 ] = +1 ;

    /*reaction 69: CH3O + H <=> S-CH2 + H2O */
    nuki[ 19 * kd + 68 ] = -1 ;
    nuki[ 2 * kd + 68 ] = -1 ;
    nuki[ 17 * kd + 68 ] = +1 ;
    nuki[ 7 * kd + 68 ] = +1 ;

    /*reaction 70: CH3O + OH <=> CH2O + H2O */
    nuki[ 19 * kd + 69 ] = -1 ;
    nuki[ 4 * kd + 69 ] = -1 ;
    nuki[ 13 * kd + 69 ] = +1 ;
    nuki[ 7 * kd + 69 ] = +1 ;

    /*reaction 71: CH3O + O <=> OH + CH2O */
    nuki[ 19 * kd + 70 ] = -1 ;
    nuki[ 5 * kd + 70 ] = -1 ;
    nuki[ 4 * kd + 70 ] = +1 ;
    nuki[ 13 * kd + 70 ] = +1 ;

    /*reaction 72: CH3O + O2 <=> CH2O + HO2 */
    nuki[ 19 * kd + 71 ] = -1 ;
    nuki[ 3 * kd + 71 ] = -1 ;
    nuki[ 13 * kd + 71 ] = +1 ;
    nuki[ 8 * kd + 71 ] = +1 ;

    /*reaction 73: CH3O + M <=> CH2O + H + M */
    nuki[ 19 * kd + 72 ] = -1 ;
    nuki[ 13 * kd + 72 ] = +1 ;
    nuki[ 2 * kd + 72 ] = +1 ;

    /*reaction 74: C2H6 + H <=> C2H5 + H2 */
    nuki[ 21 * kd + 73 ] = -1 ;
    nuki[ 2 * kd + 73 ] = -1 ;
    nuki[ 20 * kd + 73 ] = +1 ;
    nuki[ 6 * kd + 73 ] = +1 ;

    /*reaction 75: C2H6 + O <=> C2H5 + OH */
    nuki[ 21 * kd + 74 ] = -1 ;
    nuki[ 5 * kd + 74 ] = -1 ;
    nuki[ 20 * kd + 74 ] = +1 ;
    nuki[ 4 * kd + 74 ] = +1 ;

    /*reaction 76: C2H6 + OH <=> C2H5 + H2O */
    nuki[ 21 * kd + 75 ] = -1 ;
    nuki[ 4 * kd + 75 ] = -1 ;
    nuki[ 20 * kd + 75 ] = +1 ;
    nuki[ 7 * kd + 75 ] = +1 ;

    /*reaction 77: C2H6 + CH3 <=> C2H5 + CH4 */
    nuki[ 21 * kd + 76 ] = -1 ;
    nuki[ 15 * kd + 76 ] = -1 ;
    nuki[ 20 * kd + 76 ] = +1 ;
    nuki[ 14 * kd + 76 ] = +1 ;

    /*reaction 78: C2H6 (+M) <=> C2H5 + H (+M) */
    nuki[ 21 * kd + 77 ] = -1 ;
    nuki[ 20 * kd + 77 ] = +1 ;
    nuki[ 2 * kd + 77 ] = +1 ;

    /*reaction 79: C2H5 + H <=> C2H4 + H2 */
    nuki[ 20 * kd + 78 ] = -1 ;
    nuki[ 2 * kd + 78 ] = -1 ;
    nuki[ 18 * kd + 78 ] = +1 ;
    nuki[ 6 * kd + 78 ] = +1 ;

    /*reaction 80: C2H5 + O <=> C2H4 + OH */
    nuki[ 20 * kd + 79 ] = -1 ;
    nuki[ 5 * kd + 79 ] = -1 ;
    nuki[ 18 * kd + 79 ] = +1 ;
    nuki[ 4 * kd + 79 ] = +1 ;

    /*reaction 81: C2H5 + O <=> CH3 + CH2O */
    nuki[ 20 * kd + 80 ] = -1 ;
    nuki[ 5 * kd + 80 ] = -1 ;
    nuki[ 15 * kd + 80 ] = +1 ;
    nuki[ 13 * kd + 80 ] = +1 ;

    /*reaction 82: C2H5 + O2 <=> C2H4 + HO2 */
    nuki[ 20 * kd + 81 ] = -1 ;
    nuki[ 3 * kd + 81 ] = -1 ;
    nuki[ 18 * kd + 81 ] = +1 ;
    nuki[ 8 * kd + 81 ] = +1 ;

    /*reaction 83: C2H5 (+M) <=> C2H4 + H (+M) */
    nuki[ 20 * kd + 82 ] = -1 ;
    nuki[ 18 * kd + 82 ] = +1 ;
    nuki[ 2 * kd + 82 ] = +1 ;

    /*reaction 84: C2H4 + H <=> C2H3 + H2 */
    nuki[ 18 * kd + 83 ] = -1 ;
    nuki[ 2 * kd + 83 ] = -1 ;
    nuki[ 24 * kd + 83 ] = +1 ;
    nuki[ 6 * kd + 83 ] = +1 ;

    /*reaction 85: C2H4 + OH <=> C2H3 + H2O */
    nuki[ 18 * kd + 84 ] = -1 ;
    nuki[ 4 * kd + 84 ] = -1 ;
    nuki[ 24 * kd + 84 ] = +1 ;
    nuki[ 7 * kd + 84 ] = +1 ;

    /*reaction 86: C2H4 + O <=> CH3 + HCO */
    nuki[ 18 * kd + 85 ] = -1 ;
    nuki[ 5 * kd + 85 ] = -1 ;
    nuki[ 15 * kd + 85 ] = +1 ;
    nuki[ 12 * kd + 85 ] = +1 ;

    /*reaction 87: C2H4 + O <=> CH2CHO + H */
    nuki[ 18 * kd + 86 ] = -1 ;
    nuki[ 5 * kd + 86 ] = -1 ;
    nuki[ 25 * kd + 86 ] = +1 ;
    nuki[ 2 * kd + 86 ] = +1 ;

    /*reaction 88: 2 C2H4 <=> C2H3 + C2H5 */
    nuki[ 18 * kd + 87 ] = -2 ;
    nuki[ 24 * kd + 87 ] = +1 ;
    nuki[ 20 * kd + 87 ] = +1 ;

    /*reaction 89: C2H4 + O2 <=> C2H3 + HO2 */
    nuki[ 18 * kd + 88 ] = -1 ;
    nuki[ 3 * kd + 88 ] = -1 ;
    nuki[ 24 * kd + 88 ] = +1 ;
    nuki[ 8 * kd + 88 ] = +1 ;

    /*reaction 90: C2H4 + HO2 <=> C2H4O + OH */
    nuki[ 18 * kd + 89 ] = -1 ;
    nuki[ 8 * kd + 89 ] = -1 ;
    nuki[ 26 * kd + 89 ] = +1 ;
    nuki[ 4 * kd + 89 ] = +1 ;

    /*reaction 91: C2H4O + HO2 <=> CH3 + CO + H2O2 */
    nuki[ 26 * kd + 90 ] = -1 ;
    nuki[ 8 * kd + 90 ] = -1 ;
    nuki[ 15 * kd + 90 ] = +1 ;
    nuki[ 10 * kd + 90 ] = +1 ;
    nuki[ 9 * kd + 90 ] = +1 ;

    /*reaction 92: C2H4 + M <=> C2H3 + H + M */
    nuki[ 18 * kd + 91 ] = -1 ;
    nuki[ 24 * kd + 91 ] = +1 ;
    nuki[ 2 * kd + 91 ] = +1 ;

    /*reaction 93: C2H4 + M <=> C2H2 + H2 + M */
    nuki[ 18 * kd + 92 ] = -1 ;
    nuki[ 23 * kd + 92 ] = +1 ;
    nuki[ 6 * kd + 92 ] = +1 ;

    /*reaction 94: C2H3 + H <=> C2H2 + H2 */
    nuki[ 24 * kd + 93 ] = -1 ;
    nuki[ 2 * kd + 93 ] = -1 ;
    nuki[ 23 * kd + 93 ] = +1 ;
    nuki[ 6 * kd + 93 ] = +1 ;

    /*reaction 95: C2H3 (+M) <=> C2H2 + H (+M) */
    nuki[ 24 * kd + 94 ] = -1 ;
    nuki[ 23 * kd + 94 ] = +1 ;
    nuki[ 2 * kd + 94 ] = +1 ;

    /*reaction 96: C2H3 + O2 <=> CH2O + HCO */
    nuki[ 24 * kd + 95 ] = -1 ;
    nuki[ 3 * kd + 95 ] = -1 ;
    nuki[ 13 * kd + 95 ] = +1 ;
    nuki[ 12 * kd + 95 ] = +1 ;

    /*reaction 97: C2H3 + O2 <=> CH2CHO + O */
    nuki[ 24 * kd + 96 ] = -1 ;
    nuki[ 3 * kd + 96 ] = -1 ;
    nuki[ 25 * kd + 96 ] = +1 ;
    nuki[ 5 * kd + 96 ] = +1 ;

    /*reaction 98: C2H3 + O2 <=> C2H2 + HO2 */
    nuki[ 24 * kd + 97 ] = -1 ;
    nuki[ 3 * kd + 97 ] = -1 ;
    nuki[ 23 * kd + 97 ] = +1 ;
    nuki[ 8 * kd + 97 ] = +1 ;

    /*reaction 99: CH2CHO <=> CH2CO + H */
    nuki[ 25 * kd + 98 ] = -1 ;
    nuki[ 27 * kd + 98 ] = +1 ;
    nuki[ 2 * kd + 98 ] = +1 ;

    /*reaction 100: C2H2 + O <=> HCCO + H */
    nuki[ 23 * kd + 99 ] = -1 ;
    nuki[ 5 * kd + 99 ] = -1 ;
    nuki[ 28 * kd + 99 ] = +1 ;
    nuki[ 2 * kd + 99 ] = +1 ;

    /*reaction 101: C2H2 + O <=> T-CH2 + CO */
    nuki[ 23 * kd + 100 ] = -1 ;
    nuki[ 5 * kd + 100 ] = -1 ;
    nuki[ 16 * kd + 100 ] = +1 ;
    nuki[ 10 * kd + 100 ] = +1 ;

    /*reaction 102: C2H2 + O2 <=> CH2O + CO */
    nuki[ 23 * kd + 101 ] = -1 ;
    nuki[ 3 * kd + 101 ] = -1 ;
    nuki[ 13 * kd + 101 ] = +1 ;
    nuki[ 10 * kd + 101 ] = +1 ;

    /*reaction 103: C2H2 + OH <=> CH2CO + H */
    nuki[ 23 * kd + 102 ] = -1 ;
    nuki[ 4 * kd + 102 ] = -1 ;
    nuki[ 27 * kd + 102 ] = +1 ;
    nuki[ 2 * kd + 102 ] = +1 ;

    /*reaction 104: C2H2 + OH <=> C2H + H2O */
    nuki[ 23 * kd + 103 ] = -1 ;
    nuki[ 4 * kd + 103 ] = -1 ;
    nuki[ 29 * kd + 103 ] = +1 ;
    nuki[ 7 * kd + 103 ] = +1 ;

    /*reaction 105: CH2CO + H <=> CH3 + CO */
    nuki[ 27 * kd + 104 ] = -1 ;
    nuki[ 2 * kd + 104 ] = -1 ;
    nuki[ 15 * kd + 104 ] = +1 ;
    nuki[ 10 * kd + 104 ] = +1 ;

    /*reaction 106: CH2CO + O <=> T-CH2 + CO2 */
    nuki[ 27 * kd + 105 ] = -1 ;
    nuki[ 5 * kd + 105 ] = -1 ;
    nuki[ 16 * kd + 105 ] = +1 ;
    nuki[ 11 * kd + 105 ] = +1 ;

    /*reaction 107: CH2CO + O <=> HCCO + OH */
    nuki[ 27 * kd + 106 ] = -1 ;
    nuki[ 5 * kd + 106 ] = -1 ;
    nuki[ 28 * kd + 106 ] = +1 ;
    nuki[ 4 * kd + 106 ] = +1 ;

    /*reaction 108: CH2CO + CH3 <=> C2H5 + CO */
    nuki[ 27 * kd + 107 ] = -1 ;
    nuki[ 15 * kd + 107 ] = -1 ;
    nuki[ 20 * kd + 107 ] = +1 ;
    nuki[ 10 * kd + 107 ] = +1 ;

    /*reaction 109: HCCO + H <=> S-CH2 + CO */
    nuki[ 28 * kd + 108 ] = -1 ;
    nuki[ 2 * kd + 108 ] = -1 ;
    nuki[ 17 * kd + 108 ] = +1 ;
    nuki[ 10 * kd + 108 ] = +1 ;

    /*reaction 110: HCCO + OH <=> HCO + CO + H */
    nuki[ 28 * kd + 109 ] = -1 ;
    nuki[ 4 * kd + 109 ] = -1 ;
    nuki[ 12 * kd + 109 ] = +1 ;
    nuki[ 10 * kd + 109 ] = +1 ;
    nuki[ 2 * kd + 109 ] = +1 ;

    /*reaction 111: HCCO + O <=> 2 CO + H */
    nuki[ 28 * kd + 110 ] = -1 ;
    nuki[ 5 * kd + 110 ] = -1 ;
    nuki[ 10 * kd + 110 ] = +2 ;
    nuki[ 2 * kd + 110 ] = +1 ;

    /*reaction 112: HCCO + O2 <=> 2 CO + OH */
    nuki[ 28 * kd + 111 ] = -1 ;
    nuki[ 3 * kd + 111 ] = -1 ;
    nuki[ 10 * kd + 111 ] = +2 ;
    nuki[ 4 * kd + 111 ] = +1 ;

    /*reaction 113: HCCO + O2 <=> CO2 + CO + H */
    nuki[ 28 * kd + 112 ] = -1 ;
    nuki[ 3 * kd + 112 ] = -1 ;
    nuki[ 11 * kd + 112 ] = +1 ;
    nuki[ 10 * kd + 112 ] = +1 ;
    nuki[ 2 * kd + 112 ] = +1 ;

    /*reaction 114: C2H + OH <=> HCCO + H */
    nuki[ 29 * kd + 113 ] = -1 ;
    nuki[ 4 * kd + 113 ] = -1 ;
    nuki[ 28 * kd + 113 ] = +1 ;
    nuki[ 2 * kd + 113 ] = +1 ;

    /*reaction 115: C2H + O <=> CO + CH */
    nuki[ 29 * kd + 114 ] = -1 ;
    nuki[ 5 * kd + 114 ] = -1 ;
    nuki[ 10 * kd + 114 ] = +1 ;
    nuki[ 22 * kd + 114 ] = +1 ;

    /*reaction 116: C2H + O2 <=> HCCO + O */
    nuki[ 29 * kd + 115 ] = -1 ;
    nuki[ 3 * kd + 115 ] = -1 ;
    nuki[ 28 * kd + 115 ] = +1 ;
    nuki[ 5 * kd + 115 ] = +1 ;

    /*reaction 117: C2H + O2 <=> CH + CO2 */
    nuki[ 29 * kd + 116 ] = -1 ;
    nuki[ 3 * kd + 116 ] = -1 ;
    nuki[ 22 * kd + 116 ] = +1 ;
    nuki[ 11 * kd + 116 ] = +1 ;

    /*reaction 118: C2H + O2 <=> HCO + CO */
    nuki[ 29 * kd + 117 ] = -1 ;
    nuki[ 3 * kd + 117 ] = -1 ;
    nuki[ 12 * kd + 117 ] = +1 ;
    nuki[ 10 * kd + 117 ] = +1 ;

    /*reaction 119: CH2OH + H <=> CH2O + H2 */
    nuki[ 30 * kd + 118 ] = -1 ;
    nuki[ 2 * kd + 118 ] = -1 ;
    nuki[ 13 * kd + 118 ] = +1 ;
    nuki[ 6 * kd + 118 ] = +1 ;

    /*reaction 120: CH2OH + H <=> CH3 + OH */
    nuki[ 30 * kd + 119 ] = -1 ;
    nuki[ 2 * kd + 119 ] = -1 ;
    nuki[ 15 * kd + 119 ] = +1 ;
    nuki[ 4 * kd + 119 ] = +1 ;

    /*reaction 121: CH2OH + OH <=> CH2O + H2O */
    nuki[ 30 * kd + 120 ] = -1 ;
    nuki[ 4 * kd + 120 ] = -1 ;
    nuki[ 13 * kd + 120 ] = +1 ;
    nuki[ 7 * kd + 120 ] = +1 ;

    /*reaction 122: CH2OH + O2 <=> CH2O + HO2 */
    nuki[ 30 * kd + 121 ] = -1 ;
    nuki[ 3 * kd + 121 ] = -1 ;
    nuki[ 13 * kd + 121 ] = +1 ;
    nuki[ 8 * kd + 121 ] = +1 ;

    /*reaction 123: CH2OH + M <=> CH2O + H + M */
    nuki[ 30 * kd + 122 ] = -1 ;
    nuki[ 13 * kd + 122 ] = +1 ;
    nuki[ 2 * kd + 122 ] = +1 ;

    /*reaction 124: CH3O + M <=> CH2OH + M */
    nuki[ 19 * kd + 123 ] = -1 ;
    nuki[ 30 * kd + 123 ] = +1 ;

    /*reaction 125: CH2CO + OH <=> CH2OH + CO */
    nuki[ 27 * kd + 124 ] = -1 ;
    nuki[ 4 * kd + 124 ] = -1 ;
    nuki[ 30 * kd + 124 ] = +1 ;
    nuki[ 10 * kd + 124 ] = +1 ;

    /*reaction 126: CH3OH + OH <=> CH2OH + H2O */
    nuki[ 31 * kd + 125 ] = -1 ;
    nuki[ 4 * kd + 125 ] = -1 ;
    nuki[ 30 * kd + 125 ] = +1 ;
    nuki[ 7 * kd + 125 ] = +1 ;

    /*reaction 127: CH3OH + OH <=> CH3O + H2O */
    nuki[ 31 * kd + 126 ] = -1 ;
    nuki[ 4 * kd + 126 ] = -1 ;
    nuki[ 19 * kd + 126 ] = +1 ;
    nuki[ 7 * kd + 126 ] = +1 ;

    /*reaction 128: CH3OH + H <=> CH2OH + H2 */
    nuki[ 31 * kd + 127 ] = -1 ;
    nuki[ 2 * kd + 127 ] = -1 ;
    nuki[ 30 * kd + 127 ] = +1 ;
    nuki[ 6 * kd + 127 ] = +1 ;

    /*reaction 129: CH3OH + H <=> CH3O + H2 */
    nuki[ 31 * kd + 128 ] = -1 ;
    nuki[ 2 * kd + 128 ] = -1 ;
    nuki[ 19 * kd + 128 ] = +1 ;
    nuki[ 6 * kd + 128 ] = +1 ;

    /*reaction 130: CH3OH + O <=> CH2OH + OH */
    nuki[ 31 * kd + 129 ] = -1 ;
    nuki[ 5 * kd + 129 ] = -1 ;
    nuki[ 30 * kd + 129 ] = +1 ;
    nuki[ 4 * kd + 129 ] = +1 ;

    /*reaction 131: CH3OH + HO2 <=> CH2OH + H2O2 */
    nuki[ 31 * kd + 130 ] = -1 ;
    nuki[ 8 * kd + 130 ] = -1 ;
    nuki[ 30 * kd + 130 ] = +1 ;
    nuki[ 9 * kd + 130 ] = +1 ;

    /*reaction 132: CH3OH + O2 <=> CH2OH + HO2 */
    nuki[ 31 * kd + 131 ] = -1 ;
    nuki[ 3 * kd + 131 ] = -1 ;
    nuki[ 30 * kd + 131 ] = +1 ;
    nuki[ 8 * kd + 131 ] = +1 ;

    /*reaction 133: C3H4 + O <=> C2H4 + CO */
    nuki[ 32 * kd + 132 ] = -1 ;
    nuki[ 5 * kd + 132 ] = -1 ;
    nuki[ 18 * kd + 132 ] = +1 ;
    nuki[ 10 * kd + 132 ] = +1 ;

    /*reaction 134: CH3 + C2H2 <=> C3H4 + H */
    nuki[ 15 * kd + 133 ] = -1 ;
    nuki[ 23 * kd + 133 ] = -1 ;
    nuki[ 32 * kd + 133 ] = +1 ;
    nuki[ 2 * kd + 133 ] = +1 ;

    /*reaction 135: C3H4 + O <=> HCCO + CH3 */
    nuki[ 32 * kd + 134 ] = -1 ;
    nuki[ 5 * kd + 134 ] = -1 ;
    nuki[ 28 * kd + 134 ] = +1 ;
    nuki[ 15 * kd + 134 ] = +1 ;

    /*reaction 136: C3H3 + H (+M) <=> C3H4 (+M) */
    nuki[ 33 * kd + 135 ] = -1 ;
    nuki[ 2 * kd + 135 ] = -1 ;
    nuki[ 32 * kd + 135 ] = +1 ;

    /*reaction 137: C3H3 + HO2 <=> C3H4 + O2 */
    nuki[ 33 * kd + 136 ] = -1 ;
    nuki[ 8 * kd + 136 ] = -1 ;
    nuki[ 32 * kd + 136 ] = +1 ;
    nuki[ 3 * kd + 136 ] = +1 ;

    /*reaction 138: C3H4 + OH <=> C3H3 + H2O */
    nuki[ 32 * kd + 137 ] = -1 ;
    nuki[ 4 * kd + 137 ] = -1 ;
    nuki[ 33 * kd + 137 ] = +1 ;
    nuki[ 7 * kd + 137 ] = +1 ;

    /*reaction 139: C3H3 + O2 <=> CH2CO + HCO */
    nuki[ 33 * kd + 138 ] = -1 ;
    nuki[ 3 * kd + 138 ] = -1 ;
    nuki[ 27 * kd + 138 ] = +1 ;
    nuki[ 12 * kd + 138 ] = +1 ;

    /*reaction 140: C3H4 + H (+M) <=> C3H5 (+M) */
    nuki[ 32 * kd + 139 ] = -1 ;
    nuki[ 2 * kd + 139 ] = -1 ;
    nuki[ 34 * kd + 139 ] = +1 ;

    /*reaction 141: C3H5 + H <=> C3H4 + H2 */
    nuki[ 34 * kd + 140 ] = -1 ;
    nuki[ 2 * kd + 140 ] = -1 ;
    nuki[ 32 * kd + 140 ] = +1 ;
    nuki[ 6 * kd + 140 ] = +1 ;

    /*reaction 142: C3H5 + O2 <=> C3H4 + HO2 */
    nuki[ 34 * kd + 141 ] = -1 ;
    nuki[ 3 * kd + 141 ] = -1 ;
    nuki[ 32 * kd + 141 ] = +1 ;
    nuki[ 8 * kd + 141 ] = +1 ;

    /*reaction 143: C3H5 + CH3 <=> C3H4 + CH4 */
    nuki[ 34 * kd + 142 ] = -1 ;
    nuki[ 15 * kd + 142 ] = -1 ;
    nuki[ 32 * kd + 142 ] = +1 ;
    nuki[ 14 * kd + 142 ] = +1 ;

    /*reaction 144: C2H2 + CH3 (+M) <=> C3H5 (+M) */
    nuki[ 23 * kd + 143 ] = -1 ;
    nuki[ 15 * kd + 143 ] = -1 ;
    nuki[ 34 * kd + 143 ] = +1 ;

    /*reaction 145: C3H5 + OH <=> C3H4 + H2O */
    nuki[ 34 * kd + 144 ] = -1 ;
    nuki[ 4 * kd + 144 ] = -1 ;
    nuki[ 32 * kd + 144 ] = +1 ;
    nuki[ 7 * kd + 144 ] = +1 ;

    /*reaction 146: C3H3 + HCO <=> C3H4 + CO */
    nuki[ 33 * kd + 145 ] = -1 ;
    nuki[ 12 * kd + 145 ] = -1 ;
    nuki[ 32 * kd + 145 ] = +1 ;
    nuki[ 10 * kd + 145 ] = +1 ;

    /*reaction 147: C3H3 + HO2 <=> OH + CO + C2H3 */
    nuki[ 33 * kd + 146 ] = -1 ;
    nuki[ 8 * kd + 146 ] = -1 ;
    nuki[ 4 * kd + 146 ] = +1 ;
    nuki[ 10 * kd + 146 ] = +1 ;
    nuki[ 24 * kd + 146 ] = +1 ;

    /*reaction 148: C3H4 + O2 <=> CH3 + HCO + CO */
    nuki[ 32 * kd + 147 ] = -1 ;
    nuki[ 3 * kd + 147 ] = -1 ;
    nuki[ 15 * kd + 147 ] = +1 ;
    nuki[ 12 * kd + 147 ] = +1 ;
    nuki[ 10 * kd + 147 ] = +1 ;

    /*reaction 149: C3H6 + O <=> C2H5 + HCO */
    nuki[ 35 * kd + 148 ] = -1 ;
    nuki[ 5 * kd + 148 ] = -1 ;
    nuki[ 20 * kd + 148 ] = +1 ;
    nuki[ 12 * kd + 148 ] = +1 ;

    /*reaction 150: C3H6 + OH <=> C3H5 + H2O */
    nuki[ 35 * kd + 149 ] = -1 ;
    nuki[ 4 * kd + 149 ] = -1 ;
    nuki[ 34 * kd + 149 ] = +1 ;
    nuki[ 7 * kd + 149 ] = +1 ;

    /*reaction 151: C3H6 + O <=> CH2CO + CH3 + H */
    nuki[ 35 * kd + 150 ] = -1 ;
    nuki[ 5 * kd + 150 ] = -1 ;
    nuki[ 27 * kd + 150 ] = +1 ;
    nuki[ 15 * kd + 150 ] = +1 ;
    nuki[ 2 * kd + 150 ] = +1 ;

    /*reaction 152: C3H6 + H <=> C3H5 + H2 */
    nuki[ 35 * kd + 151 ] = -1 ;
    nuki[ 2 * kd + 151 ] = -1 ;
    nuki[ 34 * kd + 151 ] = +1 ;
    nuki[ 6 * kd + 151 ] = +1 ;

    /*reaction 153: C3H5 + H (+M) <=> C3H6 (+M) */
    nuki[ 34 * kd + 152 ] = -1 ;
    nuki[ 2 * kd + 152 ] = -1 ;
    nuki[ 35 * kd + 152 ] = +1 ;

    /*reaction 154: C3H5 + HO2 <=> C3H6 + O2 */
    nuki[ 34 * kd + 153 ] = -1 ;
    nuki[ 8 * kd + 153 ] = -1 ;
    nuki[ 35 * kd + 153 ] = +1 ;
    nuki[ 3 * kd + 153 ] = +1 ;

    /*reaction 155: C3H5 + HO2 <=> OH + C2H3 + CH2O */
    nuki[ 34 * kd + 154 ] = -1 ;
    nuki[ 8 * kd + 154 ] = -1 ;
    nuki[ 4 * kd + 154 ] = +1 ;
    nuki[ 24 * kd + 154 ] = +1 ;
    nuki[ 13 * kd + 154 ] = +1 ;

    /*reaction 156: C2H3 + CH3 (+M) <=> C3H6 (+M) */
    nuki[ 24 * kd + 155 ] = -1 ;
    nuki[ 15 * kd + 155 ] = -1 ;
    nuki[ 35 * kd + 155 ] = +1 ;

    /*reaction 157: C3H6 + H <=> C2H4 + CH3 */
    nuki[ 35 * kd + 156 ] = -1 ;
    nuki[ 2 * kd + 156 ] = -1 ;
    nuki[ 18 * kd + 156 ] = +1 ;
    nuki[ 15 * kd + 156 ] = +1 ;

    /*reaction 158: CH3 + C2H3 <=> C3H5 + H */
    nuki[ 15 * kd + 157 ] = -1 ;
    nuki[ 24 * kd + 157 ] = -1 ;
    nuki[ 34 * kd + 157 ] = +1 ;
    nuki[ 2 * kd + 157 ] = +1 ;

    /*reaction 159: C3H8 (+M) <=> CH3 + C2H5 (+M) */
    nuki[ 36 * kd + 158 ] = -1 ;
    nuki[ 15 * kd + 158 ] = +1 ;
    nuki[ 20 * kd + 158 ] = +1 ;

    /*reaction 160: C3H8 + O2 <=> I-C3H7 + HO2 */
    nuki[ 36 * kd + 159 ] = -1 ;
    nuki[ 3 * kd + 159 ] = -1 ;
    nuki[ 37 * kd + 159 ] = +1 ;
    nuki[ 8 * kd + 159 ] = +1 ;

    /*reaction 161: C3H8 + O2 <=> N-C3H7 + HO2 */
    nuki[ 36 * kd + 160 ] = -1 ;
    nuki[ 3 * kd + 160 ] = -1 ;
    nuki[ 38 * kd + 160 ] = +1 ;
    nuki[ 8 * kd + 160 ] = +1 ;

    /*reaction 162: C3H8 + H <=> I-C3H7 + H2 */
    nuki[ 36 * kd + 161 ] = -1 ;
    nuki[ 2 * kd + 161 ] = -1 ;
    nuki[ 37 * kd + 161 ] = +1 ;
    nuki[ 6 * kd + 161 ] = +1 ;

    /*reaction 163: C3H8 + H <=> N-C3H7 + H2 */
    nuki[ 36 * kd + 162 ] = -1 ;
    nuki[ 2 * kd + 162 ] = -1 ;
    nuki[ 38 * kd + 162 ] = +1 ;
    nuki[ 6 * kd + 162 ] = +1 ;

    /*reaction 164: C3H8 + O <=> I-C3H7 + OH */
    nuki[ 36 * kd + 163 ] = -1 ;
    nuki[ 5 * kd + 163 ] = -1 ;
    nuki[ 37 * kd + 163 ] = +1 ;
    nuki[ 4 * kd + 163 ] = +1 ;

    /*reaction 165: C3H8 + O <=> N-C3H7 + OH */
    nuki[ 36 * kd + 164 ] = -1 ;
    nuki[ 5 * kd + 164 ] = -1 ;
    nuki[ 38 * kd + 164 ] = +1 ;
    nuki[ 4 * kd + 164 ] = +1 ;

    /*reaction 166: C3H8 + OH <=> N-C3H7 + H2O */
    nuki[ 36 * kd + 165 ] = -1 ;
    nuki[ 4 * kd + 165 ] = -1 ;
    nuki[ 38 * kd + 165 ] = +1 ;
    nuki[ 7 * kd + 165 ] = +1 ;

    /*reaction 167: C3H8 + OH <=> I-C3H7 + H2O */
    nuki[ 36 * kd + 166 ] = -1 ;
    nuki[ 4 * kd + 166 ] = -1 ;
    nuki[ 37 * kd + 166 ] = +1 ;
    nuki[ 7 * kd + 166 ] = +1 ;

    /*reaction 168: C3H8 + HO2 <=> I-C3H7 + H2O2 */
    nuki[ 36 * kd + 167 ] = -1 ;
    nuki[ 8 * kd + 167 ] = -1 ;
    nuki[ 37 * kd + 167 ] = +1 ;
    nuki[ 9 * kd + 167 ] = +1 ;

    /*reaction 169: C3H8 + HO2 <=> N-C3H7 + H2O2 */
    nuki[ 36 * kd + 168 ] = -1 ;
    nuki[ 8 * kd + 168 ] = -1 ;
    nuki[ 38 * kd + 168 ] = +1 ;
    nuki[ 9 * kd + 168 ] = +1 ;

    /*reaction 170: I-C3H7 + C3H8 <=> N-C3H7 + C3H8 */
    nuki[ 37 * kd + 169 ] = -1 ;
    nuki[ 36 * kd + 169 ] = -1 ;
    nuki[ 38 * kd + 169 ] = +1 ;
    nuki[ 36 * kd + 169 ] = +1 ;

    /*reaction 171: C3H6 + H (+M) <=> I-C3H7 (+M) */
    nuki[ 35 * kd + 170 ] = -1 ;
    nuki[ 2 * kd + 170 ] = -1 ;
    nuki[ 37 * kd + 170 ] = +1 ;

    /*reaction 172: I-C3H7 + O2 <=> C3H6 + HO2 */
    nuki[ 37 * kd + 171 ] = -1 ;
    nuki[ 3 * kd + 171 ] = -1 ;
    nuki[ 35 * kd + 171 ] = +1 ;
    nuki[ 8 * kd + 171 ] = +1 ;

    /*reaction 173: N-C3H7 (+M) <=> CH3 + C2H4 (+M) */
    nuki[ 38 * kd + 172 ] = -1 ;
    nuki[ 15 * kd + 172 ] = +1 ;
    nuki[ 18 * kd + 172 ] = +1 ;

    /*reaction 174: H + C3H6 (+M) <=> N-C3H7 (+M) */
    nuki[ 2 * kd + 173 ] = -1 ;
    nuki[ 35 * kd + 173 ] = -1 ;
    nuki[ 38 * kd + 173 ] = +1 ;

    /*reaction 175: N-C3H7 + O2 <=> C3H6 + HO2 */
    nuki[ 38 * kd + 174 ] = -1 ;
    nuki[ 3 * kd + 174 ] = -1 ;
    nuki[ 35 * kd + 174 ] = +1 ;
    nuki[ 8 * kd + 174 ] = +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < 5 * 39; ++ id) {
         ncf[id] = 0; 
    }

    /*N2 */
    ncf[ 0 * kd + 0 ] = 2; /*N */

    /*AR */
    ncf[ 1 * kd + 1 ] = 1; /*AR */

    /*H */
    ncf[ 2 * kd + 2 ] = 1; /*H */

    /*O2 */
    ncf[ 3 * kd + 3 ] = 2; /*O */

    /*OH */
    ncf[ 4 * kd + 3 ] = 1; /*O */
    ncf[ 4 * kd + 2 ] = 1; /*H */

    /*O */
    ncf[ 5 * kd + 3 ] = 1; /*O */

    /*H2 */
    ncf[ 6 * kd + 2 ] = 2; /*H */

    /*H2O */
    ncf[ 7 * kd + 2 ] = 2; /*H */
    ncf[ 7 * kd + 3 ] = 1; /*O */

    /*HO2 */
    ncf[ 8 * kd + 2 ] = 1; /*H */
    ncf[ 8 * kd + 3 ] = 2; /*O */

    /*H2O2 */
    ncf[ 9 * kd + 2 ] = 2; /*H */
    ncf[ 9 * kd + 3 ] = 2; /*O */

    /*CO */
    ncf[ 10 * kd + 4 ] = 1; /*C */
    ncf[ 10 * kd + 3 ] = 1; /*O */

    /*CO2 */
    ncf[ 11 * kd + 4 ] = 1; /*C */
    ncf[ 11 * kd + 3 ] = 2; /*O */

    /*HCO */
    ncf[ 12 * kd + 2 ] = 1; /*H */
    ncf[ 12 * kd + 4 ] = 1; /*C */
    ncf[ 12 * kd + 3 ] = 1; /*O */

    /*CH2O */
    ncf[ 13 * kd + 2 ] = 2; /*H */
    ncf[ 13 * kd + 4 ] = 1; /*C */
    ncf[ 13 * kd + 3 ] = 1; /*O */

    /*CH4 */
    ncf[ 14 * kd + 4 ] = 1; /*C */
    ncf[ 14 * kd + 2 ] = 4; /*H */

    /*CH3 */
    ncf[ 15 * kd + 4 ] = 1; /*C */
    ncf[ 15 * kd + 2 ] = 3; /*H */

    /*T-CH2 */
    ncf[ 16 * kd + 4 ] = 1; /*C */
    ncf[ 16 * kd + 2 ] = 2; /*H */

    /*S-CH2 */
    ncf[ 17 * kd + 4 ] = 1; /*C */
    ncf[ 17 * kd + 2 ] = 2; /*H */

    /*C2H4 */
    ncf[ 18 * kd + 4 ] = 2; /*C */
    ncf[ 18 * kd + 2 ] = 4; /*H */

    /*CH3O */
    ncf[ 19 * kd + 4 ] = 1; /*C */
    ncf[ 19 * kd + 2 ] = 3; /*H */
    ncf[ 19 * kd + 3 ] = 1; /*O */

    /*C2H5 */
    ncf[ 20 * kd + 4 ] = 2; /*C */
    ncf[ 20 * kd + 2 ] = 5; /*H */

    /*C2H6 */
    ncf[ 21 * kd + 4 ] = 2; /*C */
    ncf[ 21 * kd + 2 ] = 6; /*H */

    /*CH */
    ncf[ 22 * kd + 4 ] = 1; /*C */
    ncf[ 22 * kd + 2 ] = 1; /*H */

    /*C2H2 */
    ncf[ 23 * kd + 4 ] = 2; /*C */
    ncf[ 23 * kd + 2 ] = 2; /*H */

    /*C2H3 */
    ncf[ 24 * kd + 4 ] = 2; /*C */
    ncf[ 24 * kd + 2 ] = 3; /*H */

    /*CH2CHO */
    ncf[ 25 * kd + 3 ] = 1; /*O */
    ncf[ 25 * kd + 2 ] = 3; /*H */
    ncf[ 25 * kd + 4 ] = 2; /*C */

    /*C2H4O */
    ncf[ 26 * kd + 4 ] = 2; /*C */
    ncf[ 26 * kd + 2 ] = 4; /*H */
    ncf[ 26 * kd + 3 ] = 1; /*O */

    /*CH2CO */
    ncf[ 27 * kd + 4 ] = 2; /*C */
    ncf[ 27 * kd + 2 ] = 2; /*H */
    ncf[ 27 * kd + 3 ] = 1; /*O */

    /*HCCO */
    ncf[ 28 * kd + 2 ] = 1; /*H */
    ncf[ 28 * kd + 4 ] = 2; /*C */
    ncf[ 28 * kd + 3 ] = 1; /*O */

    /*C2H */
    ncf[ 29 * kd + 4 ] = 2; /*C */
    ncf[ 29 * kd + 2 ] = 1; /*H */

    /*CH2OH */
    ncf[ 30 * kd + 4 ] = 1; /*C */
    ncf[ 30 * kd + 2 ] = 3; /*H */
    ncf[ 30 * kd + 3 ] = 1; /*O */

    /*CH3OH */
    ncf[ 31 * kd + 4 ] = 1; /*C */
    ncf[ 31 * kd + 2 ] = 4; /*H */
    ncf[ 31 * kd + 3 ] = 1; /*O */

    /*C3H4 */
    ncf[ 32 * kd + 4 ] = 3; /*C */
    ncf[ 32 * kd + 2 ] = 4; /*H */

    /*C3H3 */
    ncf[ 33 * kd + 4 ] = 3; /*C */
    ncf[ 33 * kd + 2 ] = 3; /*H */

    /*C3H5 */
    ncf[ 34 * kd + 4 ] = 3; /*C */
    ncf[ 34 * kd + 2 ] = 5; /*H */

    /*C3H6 */
    ncf[ 35 * kd + 4 ] = 3; /*C */
    ncf[ 35 * kd + 2 ] = 6; /*H */

    /*C3H8 */
    ncf[ 36 * kd + 4 ] = 3; /*C */
    ncf[ 36 * kd + 2 ] = 8; /*H */

    /*I-C3H7 */
    ncf[ 37 * kd + 4 ] = 3; /*C */
    ncf[ 37 * kd + 2 ] = 7; /*H */

    /*N-C3H7 */
    ncf[ 38 * kd + 4 ] = 3; /*C */
    ncf[ 38 * kd + 2 ] = 7; /*H */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * rwrk, double * a, double * b, double * e)
{

    /*reaction 1: H + O2 <=> OH + O */
    a[0] = 3.52e+16;
    b[0] = -0.7;
    e[0] = 17069.8;

    /*reaction 2: H2 + O <=> OH + H */
    a[1] = 50600;
    b[1] = 2.67;
    e[1] = 6290.63;

    /*reaction 3: H2 + OH <=> H2O + H */
    a[2] = 1.17e+09;
    b[2] = 1.3;
    e[2] = 3635.28;

    /*reaction 4: H2O + O <=> 2 OH */
    a[3] = 7.6;
    b[3] = 3.84;
    e[3] = 12779.6;

    /*reaction 5: H + O + M <=> OH + M */
    a[4] = 6.2e+16;
    b[4] = -0.6;
    e[4] = 0;

    /*reaction 6: H2 + O2 <=> 2 OH */
    a[5] = 1.7e+13;
    b[5] = 0;
    e[5] = 47813.1;

    /*reaction 7: 2 H + M <=> H2 + M */
    a[6] = 7.2e+17;
    b[6] = -1;
    e[6] = 0;

    /*reaction 8: H + OH + M <=> H2O + M */
    a[7] = 3.8e+22;
    b[7] = -2;
    e[7] = 0;

    /*reaction 9: 2 O + M <=> O2 + M */
    a[8] = 6.17e+15;
    b[8] = -0.5;
    e[8] = 0;

    /*reaction 10: H + O2 (+M) <=> HO2 (+M) */
    a[9] = 4.65e+12;
    b[9] = 0.44;
    e[9] = 0;

    /*reaction 11: O + OH + M <=> HO2 + M */
    a[10] = 1e+16;
    b[10] = 0;
    e[10] = 0;

    /*reaction 12: HO2 + H <=> 2 OH */
    a[11] = 7.08e+13;
    b[11] = 0;
    e[11] = 299.95;

    /*reaction 13: HO2 + H <=> H2 + O2 */
    a[12] = 4.28e+13;
    b[12] = 0;
    e[12] = 1410.13;

    /*reaction 14: HO2 + H <=> H2O + O */
    a[13] = 3.1e+13;
    b[13] = 0;
    e[13] = 1720.84;

    /*reaction 15: HO2 + O <=> OH + O2 */
    a[14] = 2e+13;
    b[14] = 0;
    e[14] = 0;

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    a[15] = 2.89e+13;
    b[15] = 0;
    e[15] = -497.13;

    /*reaction 17: 2 OH (+M) <=> H2O2 (+M) */
    a[16] = 7.4e+13;
    b[16] = -0.37;
    e[16] = 0;

    /*reaction 18: 2 HO2 <=> H2O2 + O2 */
    a[17] = 3.02e+12;
    b[17] = 0;
    e[17] = 1386.23;

    /*reaction 19: H2O2 + H <=> HO2 + H2 */
    a[18] = 4.79e+13;
    b[18] = 0;
    e[18] = 7958.89;

    /*reaction 20: H2O2 + H <=> H2O + OH */
    a[19] = 1e+13;
    b[19] = 0;
    e[19] = 3585.09;

    /*reaction 21: H2O2 + OH <=> H2O + HO2 */
    a[20] = 7.08e+12;
    b[20] = 0;
    e[20] = 1434.03;

    /*reaction 22: H2O2 + O <=> HO2 + OH */
    a[21] = 9.63e+06;
    b[21] = 2;
    e[21] = 3991.4;

    /*reaction 23: CO + OH <=> CO2 + H */
    a[22] = 4.4e+06;
    b[22] = 1.5;
    e[22] = -740.92;

    /*reaction 24: CO + HO2 <=> CO2 + OH */
    a[23] = 6.03e+13;
    b[23] = 0;
    e[23] = 22944.5;

    /*reaction 25: HCO + M <=> CO + H + M */
    a[24] = 1.86e+17;
    b[24] = -1;
    e[24] = 17000.5;

    /*reaction 26: HCO + H <=> CO + H2 */
    a[25] = 1e+14;
    b[25] = 0;
    e[25] = 0;

    /*reaction 27: HCO + O <=> CO + OH */
    a[26] = 3e+13;
    b[26] = 0;
    e[26] = 0;

    /*reaction 28: HCO + O <=> CO2 + H */
    a[27] = 3e+13;
    b[27] = 0;
    e[27] = 0;

    /*reaction 29: HCO + OH <=> CO + H2O */
    a[28] = 5.02e+13;
    b[28] = 0;
    e[28] = 0;

    /*reaction 30: HCO + O2 <=> CO + HO2 */
    a[29] = 3e+12;
    b[29] = 0;
    e[29] = 0;

    /*reaction 31: CH2O + M <=> HCO + H + M */
    a[30] = 6.26e+16;
    b[30] = 0;
    e[30] = 77915.9;

    /*reaction 32: CH2O + H <=> HCO + H2 */
    a[31] = 1.26e+08;
    b[31] = 1.62;
    e[31] = 2165.39;

    /*reaction 33: CH2O + O <=> HCO + OH */
    a[32] = 3.5e+13;
    b[32] = 0;
    e[32] = 3513.38;

    /*reaction 34: CH2O + OH <=> HCO + H2O */
    a[33] = 3.9e+10;
    b[33] = 0.89;
    e[33] = 406.31;

    /*reaction 35: CH4 + H <=> H2 + CH3 */
    a[34] = 13000;
    b[34] = 3;
    e[34] = 8037.76;

    /*reaction 36: CH4 + OH <=> H2O + CH3 */
    a[35] = 1.6e+07;
    b[35] = 1.83;
    e[35] = 2782.03;

    /*reaction 37: CH4 + O <=> CH3 + OH */
    a[36] = 1.9e+09;
    b[36] = 1.44;
    e[36] = 8675.91;

    /*reaction 38: CH4 + O2 <=> CH3 + HO2 */
    a[37] = 3.98e+13;
    b[37] = 0;
    e[37] = 56890.5;

    /*reaction 39: CH4 + HO2 <=> CH3 + H2O2 */
    a[38] = 9.03e+12;
    b[38] = 0;
    e[38] = 24641.5;

    /*reaction 40: CH3 + H <=> T-CH2 + H2 */
    a[39] = 1.8e+14;
    b[39] = 0;
    e[39] = 15105.2;

    /*reaction 41: CH3 + H <=> S-CH2 + H2 */
    a[40] = 1.55e+14;
    b[40] = 0;
    e[40] = 13479.9;

    /*reaction 42: CH3 + OH <=> S-CH2 + H2O */
    a[41] = 1e+13;
    b[41] = 0;
    e[41] = 2502.39;

    /*reaction 43: CH3 + O <=> CH2O + H */
    a[42] = 8.43e+13;
    b[42] = 0;
    e[42] = 0;

    /*reaction 44: CH3 + T-CH2 <=> C2H4 + H */
    a[43] = 4.22e+13;
    b[43] = 0;
    e[43] = 0;

    /*reaction 45: CH3 + HO2 <=> CH3O + OH */
    a[44] = 5e+12;
    b[44] = 0;
    e[44] = 0;

    /*reaction 46: CH3 + O2 <=> CH2O + OH */
    a[45] = 3.3e+11;
    b[45] = 0;
    e[45] = 8941.2;

    /*reaction 47: CH3 + O2 <=> CH3O + O */
    a[46] = 1.33e+14;
    b[46] = 0;
    e[46] = 31405.3;

    /*reaction 48: 2 CH3 <=> C2H4 + H2 */
    a[47] = 1e+14;
    b[47] = 0;
    e[47] = 32002.9;

    /*reaction 49: 2 CH3 <=> C2H5 + H */
    a[48] = 3.16e+13;
    b[48] = 0;
    e[48] = 14698.9;

    /*reaction 50: CH3 + H (+M) <=> CH4 (+M) */
    a[49] = 2.11e+14;
    b[49] = 0;
    e[49] = 0;

    /*reaction 51: 2 CH3 (+M) <=> C2H6 (+M) */
    a[50] = 1.81e+13;
    b[50] = 0;
    e[50] = 0;

    /*reaction 52: S-CH2 + OH <=> CH2O + H */
    a[51] = 3e+13;
    b[51] = 0;
    e[51] = 0;

    /*reaction 53: S-CH2 + O2 <=> CO + OH + H */
    a[52] = 3.13e+13;
    b[52] = 0;
    e[52] = 0;

    /*reaction 54: S-CH2 + CO2 <=> CO + CH2O */
    a[53] = 3e+12;
    b[53] = 0;
    e[53] = 0;

    /*reaction 55: S-CH2 + M <=> T-CH2 + M */
    a[54] = 6e+12;
    b[54] = 0;
    e[54] = 0;

    /*reaction 56: T-CH2 + H <=> CH + H2 */
    a[55] = 6.02e+12;
    b[55] = 0;
    e[55] = -1787.76;

    /*reaction 57: T-CH2 + OH <=> CH2O + H */
    a[56] = 2.5e+13;
    b[56] = 0;
    e[56] = 0;

    /*reaction 58: T-CH2 + OH <=> CH + H2O */
    a[57] = 1.13e+07;
    b[57] = 2;
    e[57] = 2999.52;

    /*reaction 59: T-CH2 + O <=> CO + 2 H */
    a[58] = 8e+13;
    b[58] = 0;
    e[58] = 0;

    /*reaction 60: T-CH2 + O <=> CO + H2 */
    a[59] = 4e+13;
    b[59] = 0;
    e[59] = 0;

    /*reaction 61: T-CH2 + O2 <=> CO2 + H2 */
    a[60] = 2.63e+12;
    b[60] = 0;
    e[60] = 1491.4;

    /*reaction 62: T-CH2 + O2 <=> CO + OH + H */
    a[61] = 6.58e+12;
    b[61] = 0;
    e[61] = 1491.4;

    /*reaction 63: 2 T-CH2 <=> C2H2 + 2 H */
    a[62] = 1e+14;
    b[62] = 0;
    e[62] = 0;

    /*reaction 64: CH + O <=> CO + H */
    a[63] = 4e+13;
    b[63] = 0;
    e[63] = 0;

    /*reaction 65: CH + O2 <=> HCO + O */
    a[64] = 1.77e+11;
    b[64] = 0.76;
    e[64] = -478.01;

    /*reaction 66: CH + H2O <=> CH2O + H */
    a[65] = 1.17e+15;
    b[65] = -0.75;
    e[65] = 0;

    /*reaction 67: CH + CO2 <=> HCO + CO */
    a[66] = 48;
    b[66] = 3.22;
    e[66] = -3226.58;

    /*reaction 68: CH3O + H <=> CH2O + H2 */
    a[67] = 2e+13;
    b[67] = 0;
    e[67] = 0;

    /*reaction 69: CH3O + H <=> S-CH2 + H2O */
    a[68] = 1.6e+13;
    b[68] = 0;
    e[68] = 0;

    /*reaction 70: CH3O + OH <=> CH2O + H2O */
    a[69] = 5e+12;
    b[69] = 0;
    e[69] = 0;

    /*reaction 71: CH3O + O <=> OH + CH2O */
    a[70] = 1e+13;
    b[70] = 0;
    e[70] = 0;

    /*reaction 72: CH3O + O2 <=> CH2O + HO2 */
    a[71] = 4.28e-13;
    b[71] = 7.6;
    e[71] = -3537.28;

    /*reaction 73: CH3O + M <=> CH2O + H + M */
    a[72] = 1e+13;
    b[72] = 0;
    e[72] = 13503.8;

    /*reaction 74: C2H6 + H <=> C2H5 + H2 */
    a[73] = 540;
    b[73] = 3.5;
    e[73] = 5210.33;

    /*reaction 75: C2H6 + O <=> C2H5 + OH */
    a[74] = 1.4;
    b[74] = 4.3;
    e[74] = 2772.47;

    /*reaction 76: C2H6 + OH <=> C2H5 + H2O */
    a[75] = 2.2e+07;
    b[75] = 1.9;
    e[75] = 1123.33;

    /*reaction 77: C2H6 + CH3 <=> C2H5 + CH4 */
    a[76] = 0.55;
    b[76] = 4;
    e[76] = 8293.5;

    /*reaction 78: C2H6 (+M) <=> C2H5 + H (+M) */
    a[77] = 8.85e+20;
    b[77] = -1.23;
    e[77] = 102223;

    /*reaction 79: C2H5 + H <=> C2H4 + H2 */
    a[78] = 3e+13;
    b[78] = 0;
    e[78] = 0;

    /*reaction 80: C2H5 + O <=> C2H4 + OH */
    a[79] = 3.06e+13;
    b[79] = 0;
    e[79] = 0;

    /*reaction 81: C2H5 + O <=> CH3 + CH2O */
    a[80] = 4.24e+13;
    b[80] = 0;
    e[80] = 0;

    /*reaction 82: C2H5 + O2 <=> C2H4 + HO2 */
    a[81] = 2e+12;
    b[81] = 0;
    e[81] = 4995.22;

    /*reaction 83: C2H5 (+M) <=> C2H4 + H (+M) */
    a[82] = 1.11e+10;
    b[82] = 1.037;
    e[82] = 36768.6;

    /*reaction 84: C2H4 + H <=> C2H3 + H2 */
    a[83] = 4.49e+07;
    b[83] = 2.12;
    e[83] = 13360.4;

    /*reaction 85: C2H4 + OH <=> C2H3 + H2O */
    a[84] = 553000;
    b[84] = 2.31;
    e[84] = 2963.67;

    /*reaction 86: C2H4 + O <=> CH3 + HCO */
    a[85] = 2.25e+06;
    b[85] = 2.08;
    e[85] = 0;

    /*reaction 87: C2H4 + O <=> CH2CHO + H */
    a[86] = 1.21e+06;
    b[86] = 2.08;
    e[86] = 0;

    /*reaction 88: 2 C2H4 <=> C2H3 + C2H5 */
    a[87] = 5.01e+14;
    b[87] = 0;
    e[87] = 64700.1;

    /*reaction 89: C2H4 + O2 <=> C2H3 + HO2 */
    a[88] = 4.22e+13;
    b[88] = 0;
    e[88] = 57623.1;

    /*reaction 90: C2H4 + HO2 <=> C2H4O + OH */
    a[89] = 2.23e+12;
    b[89] = 0;
    e[89] = 17189.3;

    /*reaction 91: C2H4O + HO2 <=> CH3 + CO + H2O2 */
    a[90] = 4e+12;
    b[90] = 0;
    e[90] = 17007.7;

    /*reaction 92: C2H4 + M <=> C2H3 + H + M */
    a[91] = 2.6e+17;
    b[91] = 0;
    e[91] = 96568.1;

    /*reaction 93: C2H4 + M <=> C2H2 + H2 + M */
    a[92] = 3.5e+16;
    b[92] = 0;
    e[92] = 71532;

    /*reaction 94: C2H3 + H <=> C2H2 + H2 */
    a[93] = 1.21e+13;
    b[93] = 0;
    e[93] = 0;

    /*reaction 95: C2H3 (+M) <=> C2H2 + H (+M) */
    a[94] = 6.38e+09;
    b[94] = 1;
    e[94] = 37626.7;

    /*reaction 96: C2H3 + O2 <=> CH2O + HCO */
    a[95] = 1.7e+29;
    b[95] = -5.312;
    e[95] = 6503.11;

    /*reaction 97: C2H3 + O2 <=> CH2CHO + O */
    a[96] = 7e+14;
    b[96] = -0.611;
    e[96] = 5262.43;

    /*reaction 98: C2H3 + O2 <=> C2H2 + HO2 */
    a[97] = 5.19e+15;
    b[97] = -1.26;
    e[97] = 3312.62;

    /*reaction 99: CH2CHO <=> CH2CO + H */
    a[98] = 1.047e+37;
    b[98] = -7.189;
    e[98] = 44340.3;

    /*reaction 100: C2H2 + O <=> HCCO + H */
    a[99] = 4e+14;
    b[99] = 0;
    e[99] = 10659.7;

    /*reaction 101: C2H2 + O <=> T-CH2 + CO */
    a[100] = 1.6e+14;
    b[100] = 0;
    e[100] = 9894.84;

    /*reaction 102: C2H2 + O2 <=> CH2O + CO */
    a[101] = 4.6e+15;
    b[101] = -0.54;
    e[101] = 44933.1;

    /*reaction 103: C2H2 + OH <=> CH2CO + H */
    a[102] = 1.9e+07;
    b[102] = 1.7;
    e[102] = 999.04;

    /*reaction 104: C2H2 + OH <=> C2H + H2O */
    a[103] = 3.37e+07;
    b[103] = 2;
    e[103] = 14001;

    /*reaction 105: CH2CO + H <=> CH3 + CO */
    a[104] = 1.5e+09;
    b[104] = 1.43;
    e[104] = 2688.81;

    /*reaction 106: CH2CO + O <=> T-CH2 + CO2 */
    a[105] = 2e+13;
    b[105] = 0;
    e[105] = 2294.46;

    /*reaction 107: CH2CO + O <=> HCCO + OH */
    a[106] = 1e+13;
    b[106] = 0;
    e[106] = 2000.48;

    /*reaction 108: CH2CO + CH3 <=> C2H5 + CO */
    a[107] = 9e+10;
    b[107] = 0;
    e[107] = 0;

    /*reaction 109: HCCO + H <=> S-CH2 + CO */
    a[108] = 1.5e+14;
    b[108] = 0;
    e[108] = 0;

    /*reaction 110: HCCO + OH <=> HCO + CO + H */
    a[109] = 2e+12;
    b[109] = 0;
    e[109] = 0;

    /*reaction 111: HCCO + O <=> 2 CO + H */
    a[110] = 9.64e+13;
    b[110] = 0;
    e[110] = 0;

    /*reaction 112: HCCO + O2 <=> 2 CO + OH */
    a[111] = 2.88e+07;
    b[111] = 1.7;
    e[111] = 1001.43;

    /*reaction 113: HCCO + O2 <=> CO2 + CO + H */
    a[112] = 1.4e+07;
    b[112] = 1.7;
    e[112] = 1001.43;

    /*reaction 114: C2H + OH <=> HCCO + H */
    a[113] = 2e+13;
    b[113] = 0;
    e[113] = 0;

    /*reaction 115: C2H + O <=> CO + CH */
    a[114] = 1.02e+13;
    b[114] = 0;
    e[114] = 0;

    /*reaction 116: C2H + O2 <=> HCCO + O */
    a[115] = 6.02e+11;
    b[115] = 0;
    e[115] = 0;

    /*reaction 117: C2H + O2 <=> CH + CO2 */
    a[116] = 4.5e+15;
    b[116] = 0;
    e[116] = 25095.6;

    /*reaction 118: C2H + O2 <=> HCO + CO */
    a[117] = 2.41e+12;
    b[117] = 0;
    e[117] = 0;

    /*reaction 119: CH2OH + H <=> CH2O + H2 */
    a[118] = 3e+13;
    b[118] = 0;
    e[118] = 0;

    /*reaction 120: CH2OH + H <=> CH3 + OH */
    a[119] = 1.75e+14;
    b[119] = 0;
    e[119] = 2796.37;

    /*reaction 121: CH2OH + OH <=> CH2O + H2O */
    a[120] = 2.4e+13;
    b[120] = 0;
    e[120] = 0;

    /*reaction 122: CH2OH + O2 <=> CH2O + HO2 */
    a[121] = 5e+12;
    b[121] = 0;
    e[121] = 0;

    /*reaction 123: CH2OH + M <=> CH2O + H + M */
    a[122] = 5e+13;
    b[122] = 0;
    e[122] = 25119.5;

    /*reaction 124: CH3O + M <=> CH2OH + M */
    a[123] = 1e+14;
    b[123] = 0;
    e[123] = 19120.5;

    /*reaction 125: CH2CO + OH <=> CH2OH + CO */
    a[124] = 1.02e+13;
    b[124] = 0;
    e[124] = 0;

    /*reaction 126: CH3OH + OH <=> CH2OH + H2O */
    a[125] = 1.44e+06;
    b[125] = 2;
    e[125] = -838.91;

    /*reaction 127: CH3OH + OH <=> CH3O + H2O */
    a[126] = 6.3e+06;
    b[126] = 2;
    e[126] = 1505.74;

    /*reaction 128: CH3OH + H <=> CH2OH + H2 */
    a[127] = 1.64e+07;
    b[127] = 2;
    e[127] = 4517.21;

    /*reaction 129: CH3OH + H <=> CH3O + H2 */
    a[128] = 3.83e+07;
    b[128] = 2;
    e[128] = 5855.64;

    /*reaction 130: CH3OH + O <=> CH2OH + OH */
    a[129] = 1e+13;
    b[129] = 0;
    e[129] = 4684.51;

    /*reaction 131: CH3OH + HO2 <=> CH2OH + H2O2 */
    a[130] = 6.2e+12;
    b[130] = 0;
    e[130] = 19383.4;

    /*reaction 132: CH3OH + O2 <=> CH2OH + HO2 */
    a[131] = 2e+13;
    b[131] = 0;
    e[131] = 44933.1;

    /*reaction 133: C3H4 + O <=> C2H4 + CO */
    a[132] = 2e+07;
    b[132] = 1.8;
    e[132] = 1000;

    /*reaction 134: CH3 + C2H2 <=> C3H4 + H */
    a[133] = 2.56e+09;
    b[133] = 1.1;
    e[133] = 13643.9;

    /*reaction 135: C3H4 + O <=> HCCO + CH3 */
    a[134] = 7.3e+12;
    b[134] = 0;
    e[134] = 2250;

    /*reaction 136: C3H3 + H (+M) <=> C3H4 (+M) */
    a[135] = 3e+13;
    b[135] = 0;
    e[135] = 0;

    /*reaction 137: C3H3 + HO2 <=> C3H4 + O2 */
    a[136] = 2.5e+12;
    b[136] = 0;
    e[136] = 0;

    /*reaction 138: C3H4 + OH <=> C3H3 + H2O */
    a[137] = 5.3e+06;
    b[137] = 2;
    e[137] = 2000;

    /*reaction 139: C3H3 + O2 <=> CH2CO + HCO */
    a[138] = 3e+10;
    b[138] = 0;
    e[138] = 2868.07;

    /*reaction 140: C3H4 + H (+M) <=> C3H5 (+M) */
    a[139] = 4e+13;
    b[139] = 0;
    e[139] = 0;

    /*reaction 141: C3H5 + H <=> C3H4 + H2 */
    a[140] = 1.8e+13;
    b[140] = 0;
    e[140] = 0;

    /*reaction 142: C3H5 + O2 <=> C3H4 + HO2 */
    a[141] = 4.99e+15;
    b[141] = -1.4;
    e[141] = 22428.1;

    /*reaction 143: C3H5 + CH3 <=> C3H4 + CH4 */
    a[142] = 3e+12;
    b[142] = -0.32;
    e[142] = -130.98;

    /*reaction 144: C2H2 + CH3 (+M) <=> C3H5 (+M) */
    a[143] = 6e+08;
    b[143] = 0;
    e[143] = 0;

    /*reaction 145: C3H5 + OH <=> C3H4 + H2O */
    a[144] = 6e+12;
    b[144] = 0;
    e[144] = 0;

    /*reaction 146: C3H3 + HCO <=> C3H4 + CO */
    a[145] = 2.5e+13;
    b[145] = 0;
    e[145] = 0;

    /*reaction 147: C3H3 + HO2 <=> OH + CO + C2H3 */
    a[146] = 8e+11;
    b[146] = 0;
    e[146] = 0;

    /*reaction 148: C3H4 + O2 <=> CH3 + HCO + CO */
    a[147] = 4e+14;
    b[147] = 0;
    e[147] = 41826;

    /*reaction 149: C3H6 + O <=> C2H5 + HCO */
    a[148] = 3.5e+07;
    b[148] = 1.65;
    e[148] = -972.75;

    /*reaction 150: C3H6 + OH <=> C3H5 + H2O */
    a[149] = 3.1e+06;
    b[149] = 2;
    e[149] = -298.28;

    /*reaction 151: C3H6 + O <=> CH2CO + CH3 + H */
    a[150] = 1.2e+08;
    b[150] = 1.65;
    e[150] = 327.44;

    /*reaction 152: C3H6 + H <=> C3H5 + H2 */
    a[151] = 170000;
    b[151] = 2.5;
    e[151] = 2492.83;

    /*reaction 153: C3H5 + H (+M) <=> C3H6 (+M) */
    a[152] = 2e+14;
    b[152] = 0;
    e[152] = 0;

    /*reaction 154: C3H5 + HO2 <=> C3H6 + O2 */
    a[153] = 2.66e+12;
    b[153] = 0;
    e[153] = 0;

    /*reaction 155: C3H5 + HO2 <=> OH + C2H3 + CH2O */
    a[154] = 3e+12;
    b[154] = 0;
    e[154] = 0;

    /*reaction 156: C2H3 + CH3 (+M) <=> C3H6 (+M) */
    a[155] = 2.5e+13;
    b[155] = 0;
    e[155] = 0;

    /*reaction 157: C3H6 + H <=> C2H4 + CH3 */
    a[156] = 1.6e+22;
    b[156] = -2.39;
    e[156] = 11185.5;

    /*reaction 158: CH3 + C2H3 <=> C3H5 + H */
    a[157] = 1.5e+24;
    b[157] = -2.83;
    e[157] = 18618.5;

    /*reaction 159: C3H8 (+M) <=> CH3 + C2H5 (+M) */
    a[158] = 1.1e+17;
    b[158] = 0;
    e[158] = 84392.9;

    /*reaction 160: C3H8 + O2 <=> I-C3H7 + HO2 */
    a[159] = 4e+13;
    b[159] = 0;
    e[159] = 47500;

    /*reaction 161: C3H8 + O2 <=> N-C3H7 + HO2 */
    a[160] = 4e+13;
    b[160] = 0;
    e[160] = 50932.1;

    /*reaction 162: C3H8 + H <=> I-C3H7 + H2 */
    a[161] = 1.3e+06;
    b[161] = 2.4;
    e[161] = 4471.08;

    /*reaction 163: C3H8 + H <=> N-C3H7 + H2 */
    a[162] = 1.33e+06;
    b[162] = 2.54;
    e[162] = 6761.47;

    /*reaction 164: C3H8 + O <=> I-C3H7 + OH */
    a[163] = 47600;
    b[163] = 2.71;
    e[163] = 2107.31;

    /*reaction 165: C3H8 + O <=> N-C3H7 + OH */
    a[164] = 190000;
    b[164] = 2.68;
    e[164] = 3718.45;

    /*reaction 166: C3H8 + OH <=> N-C3H7 + H2O */
    a[165] = 1400;
    b[165] = 2.66;
    e[165] = 527.25;

    /*reaction 167: C3H8 + OH <=> I-C3H7 + H2O */
    a[166] = 27000;
    b[166] = 2.39;
    e[166] = 393.16;

    /*reaction 168: C3H8 + HO2 <=> I-C3H7 + H2O2 */
    a[167] = 9640;
    b[167] = 2.6;
    e[167] = 13917.3;

    /*reaction 169: C3H8 + HO2 <=> N-C3H7 + H2O2 */
    a[168] = 47600;
    b[168] = 2.55;
    e[168] = 16491.4;

    /*reaction 170: I-C3H7 + C3H8 <=> N-C3H7 + C3H8 */
    a[169] = 0.0084;
    b[169] = 4.2;
    e[169] = 8675.91;

    /*reaction 171: C3H6 + H (+M) <=> I-C3H7 (+M) */
    a[170] = 1.33e+13;
    b[170] = 0;
    e[170] = 1560.71;

    /*reaction 172: I-C3H7 + O2 <=> C3H6 + HO2 */
    a[171] = 1.3e+11;
    b[171] = 0;
    e[171] = 0;

    /*reaction 173: N-C3H7 (+M) <=> CH3 + C2H4 (+M) */
    a[172] = 1.23e+13;
    b[172] = -0.1;
    e[172] = 30210.3;

    /*reaction 174: H + C3H6 (+M) <=> N-C3H7 (+M) */
    a[173] = 1.33e+13;
    b[173] = 0;
    e[173] = 3260.04;

    /*reaction 175: N-C3H7 + O2 <=> C3H6 + HO2 */
    a[174] = 9e+10;
    b[174] = 0;
    e[174] = 0;

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * T, double * C, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[39]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> OH + O */
    /*eqcon[0] *= 1;  */

    /*reaction 2: H2 + O <=> OH + H */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */

    /*reaction 4: H2O + O <=> 2 OH */
    /*eqcon[3] *= 1;  */

    /*reaction 5: H + O + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2 + O2 <=> 2 OH */
    /*eqcon[5] *= 1;  */

    /*reaction 7: 2 H + M <=> H2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + OH + M <=> H2O + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: 2 O + M <=> O2 + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 (+M) <=> HO2 (+M) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: O + OH + M <=> HO2 + M */
    eqcon[10] *= 1e+06; 

    /*reaction 12: HO2 + H <=> 2 OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + H <=> H2 + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + H <=> H2O + O */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + O <=> OH + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[16] *= 1e+06; 

    /*reaction 18: 2 HO2 <=> H2O2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + H <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + H <=> H2O + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2O2 + O <=> HO2 + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: CO + OH <=> CO2 + H */
    /*eqcon[22] *= 1;  */

    /*reaction 24: CO + HO2 <=> CO2 + OH */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HCO + M <=> CO + H + M */
    eqcon[24] *= 1e-06; 

    /*reaction 26: HCO + H <=> CO + H2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: HCO + O <=> CO + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: HCO + O <=> CO2 + H */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HCO + OH <=> CO + H2O */
    /*eqcon[28] *= 1;  */

    /*reaction 30: HCO + O2 <=> CO + HO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CH2O + M <=> HCO + H + M */
    eqcon[30] *= 1e-06; 

    /*reaction 32: CH2O + H <=> HCO + H2 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: CH2O + O <=> HCO + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: CH2O + OH <=> HCO + H2O */
    /*eqcon[33] *= 1;  */

    /*reaction 35: CH4 + H <=> H2 + CH3 */
    /*eqcon[34] *= 1;  */

    /*reaction 36: CH4 + OH <=> H2O + CH3 */
    /*eqcon[35] *= 1;  */

    /*reaction 37: CH4 + O <=> CH3 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: CH4 + O2 <=> CH3 + HO2 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: CH4 + HO2 <=> CH3 + H2O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: CH3 + H <=> T-CH2 + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: CH3 + H <=> S-CH2 + H2 */
    /*eqcon[40] *= 1;  */

    /*reaction 42: CH3 + OH <=> S-CH2 + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: CH3 + O <=> CH2O + H */
    /*eqcon[42] *= 1;  */

    /*reaction 44: CH3 + T-CH2 <=> C2H4 + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: CH3 + HO2 <=> CH3O + OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: CH3 + O2 <=> CH2O + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: CH3 + O2 <=> CH3O + O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: 2 CH3 <=> C2H4 + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: 2 CH3 <=> C2H5 + H */
    /*eqcon[48] *= 1;  */

    /*reaction 50: CH3 + H (+M) <=> CH4 (+M) */
    eqcon[49] *= 1e+06; 

    /*reaction 51: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[50] *= 1e+06; 

    /*reaction 52: S-CH2 + OH <=> CH2O + H */
    /*eqcon[51] *= 1;  */

    /*reaction 53: S-CH2 + O2 <=> CO + OH + H */
    eqcon[52] *= 1e-06; 

    /*reaction 54: S-CH2 + CO2 <=> CO + CH2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: S-CH2 + M <=> T-CH2 + M */
    /*eqcon[54] *= 1;  */

    /*reaction 56: T-CH2 + H <=> CH + H2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: T-CH2 + OH <=> CH2O + H */
    /*eqcon[56] *= 1;  */

    /*reaction 58: T-CH2 + OH <=> CH + H2O */
    /*eqcon[57] *= 1;  */

    /*reaction 59: T-CH2 + O <=> CO + 2 H */
    eqcon[58] *= 1e-06; 

    /*reaction 60: T-CH2 + O <=> CO + H2 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: T-CH2 + O2 <=> CO2 + H2 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: T-CH2 + O2 <=> CO + OH + H */
    eqcon[61] *= 1e-06; 

    /*reaction 63: 2 T-CH2 <=> C2H2 + 2 H */
    eqcon[62] *= 1e-06; 

    /*reaction 64: CH + O <=> CO + H */
    /*eqcon[63] *= 1;  */

    /*reaction 65: CH + O2 <=> HCO + O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: CH + H2O <=> CH2O + H */
    /*eqcon[65] *= 1;  */

    /*reaction 67: CH + CO2 <=> HCO + CO */
    /*eqcon[66] *= 1;  */

    /*reaction 68: CH3O + H <=> CH2O + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: CH3O + H <=> S-CH2 + H2O */
    /*eqcon[68] *= 1;  */

    /*reaction 70: CH3O + OH <=> CH2O + H2O */
    /*eqcon[69] *= 1;  */

    /*reaction 71: CH3O + O <=> OH + CH2O */
    /*eqcon[70] *= 1;  */

    /*reaction 72: CH3O + O2 <=> CH2O + HO2 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: CH3O + M <=> CH2O + H + M */
    eqcon[72] *= 1e-06; 

    /*reaction 74: C2H6 + H <=> C2H5 + H2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: C2H6 + O <=> C2H5 + OH */
    /*eqcon[74] *= 1;  */

    /*reaction 76: C2H6 + OH <=> C2H5 + H2O */
    /*eqcon[75] *= 1;  */

    /*reaction 77: C2H6 + CH3 <=> C2H5 + CH4 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: C2H6 (+M) <=> C2H5 + H (+M) */
    eqcon[77] *= 1e-06; 

    /*reaction 79: C2H5 + H <=> C2H4 + H2 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: C2H5 + O <=> C2H4 + OH */
    /*eqcon[79] *= 1;  */

    /*reaction 81: C2H5 + O <=> CH3 + CH2O */
    /*eqcon[80] *= 1;  */

    /*reaction 82: C2H5 + O2 <=> C2H4 + HO2 */
    /*eqcon[81] *= 1;  */

    /*reaction 83: C2H5 (+M) <=> C2H4 + H (+M) */
    eqcon[82] *= 1e-06; 

    /*reaction 84: C2H4 + H <=> C2H3 + H2 */
    /*eqcon[83] *= 1;  */

    /*reaction 85: C2H4 + OH <=> C2H3 + H2O */
    /*eqcon[84] *= 1;  */

    /*reaction 86: C2H4 + O <=> CH3 + HCO */
    /*eqcon[85] *= 1;  */

    /*reaction 87: C2H4 + O <=> CH2CHO + H */
    /*eqcon[86] *= 1;  */

    /*reaction 88: 2 C2H4 <=> C2H3 + C2H5 */
    /*eqcon[87] *= 1;  */

    /*reaction 89: C2H4 + O2 <=> C2H3 + HO2 */
    /*eqcon[88] *= 1;  */

    /*reaction 90: C2H4 + HO2 <=> C2H4O + OH */
    /*eqcon[89] *= 1;  */

    /*reaction 91: C2H4O + HO2 <=> CH3 + CO + H2O2 */
    eqcon[90] *= 1e-06; 

    /*reaction 92: C2H4 + M <=> C2H3 + H + M */
    eqcon[91] *= 1e-06; 

    /*reaction 93: C2H4 + M <=> C2H2 + H2 + M */
    eqcon[92] *= 1e-06; 

    /*reaction 94: C2H3 + H <=> C2H2 + H2 */
    /*eqcon[93] *= 1;  */

    /*reaction 95: C2H3 (+M) <=> C2H2 + H (+M) */
    eqcon[94] *= 1e-06; 

    /*reaction 96: C2H3 + O2 <=> CH2O + HCO */
    /*eqcon[95] *= 1;  */

    /*reaction 97: C2H3 + O2 <=> CH2CHO + O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: C2H3 + O2 <=> C2H2 + HO2 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: CH2CHO <=> CH2CO + H */
    eqcon[98] *= 1e-06; 

    /*reaction 100: C2H2 + O <=> HCCO + H */
    /*eqcon[99] *= 1;  */

    /*reaction 101: C2H2 + O <=> T-CH2 + CO */
    /*eqcon[100] *= 1;  */

    /*reaction 102: C2H2 + O2 <=> CH2O + CO */
    /*eqcon[101] *= 1;  */

    /*reaction 103: C2H2 + OH <=> CH2CO + H */
    /*eqcon[102] *= 1;  */

    /*reaction 104: C2H2 + OH <=> C2H + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: CH2CO + H <=> CH3 + CO */
    /*eqcon[104] *= 1;  */

    /*reaction 106: CH2CO + O <=> T-CH2 + CO2 */
    /*eqcon[105] *= 1;  */

    /*reaction 107: CH2CO + O <=> HCCO + OH */
    /*eqcon[106] *= 1;  */

    /*reaction 108: CH2CO + CH3 <=> C2H5 + CO */
    /*eqcon[107] *= 1;  */

    /*reaction 109: HCCO + H <=> S-CH2 + CO */
    /*eqcon[108] *= 1;  */

    /*reaction 110: HCCO + OH <=> HCO + CO + H */
    eqcon[109] *= 1e-06; 

    /*reaction 111: HCCO + O <=> 2 CO + H */
    eqcon[110] *= 1e-06; 

    /*reaction 112: HCCO + O2 <=> 2 CO + OH */
    eqcon[111] *= 1e-06; 

    /*reaction 113: HCCO + O2 <=> CO2 + CO + H */
    eqcon[112] *= 1e-06; 

    /*reaction 114: C2H + OH <=> HCCO + H */
    /*eqcon[113] *= 1;  */

    /*reaction 115: C2H + O <=> CO + CH */
    /*eqcon[114] *= 1;  */

    /*reaction 116: C2H + O2 <=> HCCO + O */
    /*eqcon[115] *= 1;  */

    /*reaction 117: C2H + O2 <=> CH + CO2 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: C2H + O2 <=> HCO + CO */
    /*eqcon[117] *= 1;  */

    /*reaction 119: CH2OH + H <=> CH2O + H2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: CH2OH + H <=> CH3 + OH */
    /*eqcon[119] *= 1;  */

    /*reaction 121: CH2OH + OH <=> CH2O + H2O */
    /*eqcon[120] *= 1;  */

    /*reaction 122: CH2OH + O2 <=> CH2O + HO2 */
    /*eqcon[121] *= 1;  */

    /*reaction 123: CH2OH + M <=> CH2O + H + M */
    eqcon[122] *= 1e-06; 

    /*reaction 124: CH3O + M <=> CH2OH + M */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH2CO + OH <=> CH2OH + CO */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH3OH + OH <=> CH2OH + H2O */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH3OH + OH <=> CH3O + H2O */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH3OH + H <=> CH2OH + H2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH3OH + H <=> CH3O + H2 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH3OH + O <=> CH2OH + OH */
    /*eqcon[129] *= 1;  */

    /*reaction 131: CH3OH + HO2 <=> CH2OH + H2O2 */
    /*eqcon[130] *= 1;  */

    /*reaction 132: CH3OH + O2 <=> CH2OH + HO2 */
    /*eqcon[131] *= 1;  */

    /*reaction 133: C3H4 + O <=> C2H4 + CO */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH3 + C2H2 <=> C3H4 + H */
    /*eqcon[133] *= 1;  */

    /*reaction 135: C3H4 + O <=> HCCO + CH3 */
    /*eqcon[134] *= 1;  */

    /*reaction 136: C3H3 + H (+M) <=> C3H4 (+M) */
    eqcon[135] *= 1e+06; 

    /*reaction 137: C3H3 + HO2 <=> C3H4 + O2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: C3H4 + OH <=> C3H3 + H2O */
    /*eqcon[137] *= 1;  */

    /*reaction 139: C3H3 + O2 <=> CH2CO + HCO */
    /*eqcon[138] *= 1;  */

    /*reaction 140: C3H4 + H (+M) <=> C3H5 (+M) */
    eqcon[139] *= 1e+06; 

    /*reaction 141: C3H5 + H <=> C3H4 + H2 */
    /*eqcon[140] *= 1;  */

    /*reaction 142: C3H5 + O2 <=> C3H4 + HO2 */
    /*eqcon[141] *= 1;  */

    /*reaction 143: C3H5 + CH3 <=> C3H4 + CH4 */
    /*eqcon[142] *= 1;  */

    /*reaction 144: C2H2 + CH3 (+M) <=> C3H5 (+M) */
    eqcon[143] *= 1e+06; 

    /*reaction 145: C3H5 + OH <=> C3H4 + H2O */
    /*eqcon[144] *= 1;  */

    /*reaction 146: C3H3 + HCO <=> C3H4 + CO */
    /*eqcon[145] *= 1;  */

    /*reaction 147: C3H3 + HO2 <=> OH + CO + C2H3 */
    eqcon[146] *= 1e-06; 

    /*reaction 148: C3H4 + O2 <=> CH3 + HCO + CO */
    eqcon[147] *= 1e-06; 

    /*reaction 149: C3H6 + O <=> C2H5 + HCO */
    /*eqcon[148] *= 1;  */

    /*reaction 150: C3H6 + OH <=> C3H5 + H2O */
    /*eqcon[149] *= 1;  */

    /*reaction 151: C3H6 + O <=> CH2CO + CH3 + H */
    eqcon[150] *= 1e-06; 

    /*reaction 152: C3H6 + H <=> C3H5 + H2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: C3H5 + H (+M) <=> C3H6 (+M) */
    eqcon[152] *= 1e+06; 

    /*reaction 154: C3H5 + HO2 <=> C3H6 + O2 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: C3H5 + HO2 <=> OH + C2H3 + CH2O */
    eqcon[154] *= 1e-06; 

    /*reaction 156: C2H3 + CH3 (+M) <=> C3H6 (+M) */
    eqcon[155] *= 1e+06; 

    /*reaction 157: C3H6 + H <=> C2H4 + CH3 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: CH3 + C2H3 <=> C3H5 + H */
    /*eqcon[157] *= 1;  */

    /*reaction 159: C3H8 (+M) <=> CH3 + C2H5 (+M) */
    eqcon[158] *= 1e-06; 

    /*reaction 160: C3H8 + O2 <=> I-C3H7 + HO2 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: C3H8 + O2 <=> N-C3H7 + HO2 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: C3H8 + H <=> I-C3H7 + H2 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: C3H8 + H <=> N-C3H7 + H2 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: C3H8 + O <=> I-C3H7 + OH */
    /*eqcon[163] *= 1;  */

    /*reaction 165: C3H8 + O <=> N-C3H7 + OH */
    /*eqcon[164] *= 1;  */

    /*reaction 166: C3H8 + OH <=> N-C3H7 + H2O */
    /*eqcon[165] *= 1;  */

    /*reaction 167: C3H8 + OH <=> I-C3H7 + H2O */
    /*eqcon[166] *= 1;  */

    /*reaction 168: C3H8 + HO2 <=> I-C3H7 + H2O2 */
    /*eqcon[167] *= 1;  */

    /*reaction 169: C3H8 + HO2 <=> N-C3H7 + H2O2 */
    /*eqcon[168] *= 1;  */

    /*reaction 170: I-C3H7 + C3H8 <=> N-C3H7 + C3H8 */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C3H6 + H (+M) <=> I-C3H7 (+M) */
    eqcon[170] *= 1e+06; 

    /*reaction 172: I-C3H7 + O2 <=> C3H6 + HO2 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: N-C3H7 (+M) <=> CH3 + C2H4 (+M) */
    eqcon[172] *= 1e-06; 

    /*reaction 174: H + C3H6 (+M) <=> N-C3H7 (+M) */
    eqcon[173] *= 1e+06; 

    /*reaction 175: N-C3H7 + O2 <=> C3H6 + HO2 */
    /*eqcon[174] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[39]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> OH + O */
    /*eqcon[0] *= 1;  */

    /*reaction 2: H2 + O <=> OH + H */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */

    /*reaction 4: H2O + O <=> 2 OH */
    /*eqcon[3] *= 1;  */

    /*reaction 5: H + O + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2 + O2 <=> 2 OH */
    /*eqcon[5] *= 1;  */

    /*reaction 7: 2 H + M <=> H2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + OH + M <=> H2O + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: 2 O + M <=> O2 + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 (+M) <=> HO2 (+M) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: O + OH + M <=> HO2 + M */
    eqcon[10] *= 1e+06; 

    /*reaction 12: HO2 + H <=> 2 OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + H <=> H2 + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + H <=> H2O + O */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + O <=> OH + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[16] *= 1e+06; 

    /*reaction 18: 2 HO2 <=> H2O2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + H <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + H <=> H2O + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2O2 + O <=> HO2 + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: CO + OH <=> CO2 + H */
    /*eqcon[22] *= 1;  */

    /*reaction 24: CO + HO2 <=> CO2 + OH */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HCO + M <=> CO + H + M */
    eqcon[24] *= 1e-06; 

    /*reaction 26: HCO + H <=> CO + H2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: HCO + O <=> CO + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: HCO + O <=> CO2 + H */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HCO + OH <=> CO + H2O */
    /*eqcon[28] *= 1;  */

    /*reaction 30: HCO + O2 <=> CO + HO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CH2O + M <=> HCO + H + M */
    eqcon[30] *= 1e-06; 

    /*reaction 32: CH2O + H <=> HCO + H2 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: CH2O + O <=> HCO + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: CH2O + OH <=> HCO + H2O */
    /*eqcon[33] *= 1;  */

    /*reaction 35: CH4 + H <=> H2 + CH3 */
    /*eqcon[34] *= 1;  */

    /*reaction 36: CH4 + OH <=> H2O + CH3 */
    /*eqcon[35] *= 1;  */

    /*reaction 37: CH4 + O <=> CH3 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: CH4 + O2 <=> CH3 + HO2 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: CH4 + HO2 <=> CH3 + H2O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: CH3 + H <=> T-CH2 + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: CH3 + H <=> S-CH2 + H2 */
    /*eqcon[40] *= 1;  */

    /*reaction 42: CH3 + OH <=> S-CH2 + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: CH3 + O <=> CH2O + H */
    /*eqcon[42] *= 1;  */

    /*reaction 44: CH3 + T-CH2 <=> C2H4 + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: CH3 + HO2 <=> CH3O + OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: CH3 + O2 <=> CH2O + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: CH3 + O2 <=> CH3O + O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: 2 CH3 <=> C2H4 + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: 2 CH3 <=> C2H5 + H */
    /*eqcon[48] *= 1;  */

    /*reaction 50: CH3 + H (+M) <=> CH4 (+M) */
    eqcon[49] *= 1e+06; 

    /*reaction 51: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[50] *= 1e+06; 

    /*reaction 52: S-CH2 + OH <=> CH2O + H */
    /*eqcon[51] *= 1;  */

    /*reaction 53: S-CH2 + O2 <=> CO + OH + H */
    eqcon[52] *= 1e-06; 

    /*reaction 54: S-CH2 + CO2 <=> CO + CH2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: S-CH2 + M <=> T-CH2 + M */
    /*eqcon[54] *= 1;  */

    /*reaction 56: T-CH2 + H <=> CH + H2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: T-CH2 + OH <=> CH2O + H */
    /*eqcon[56] *= 1;  */

    /*reaction 58: T-CH2 + OH <=> CH + H2O */
    /*eqcon[57] *= 1;  */

    /*reaction 59: T-CH2 + O <=> CO + 2 H */
    eqcon[58] *= 1e-06; 

    /*reaction 60: T-CH2 + O <=> CO + H2 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: T-CH2 + O2 <=> CO2 + H2 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: T-CH2 + O2 <=> CO + OH + H */
    eqcon[61] *= 1e-06; 

    /*reaction 63: 2 T-CH2 <=> C2H2 + 2 H */
    eqcon[62] *= 1e-06; 

    /*reaction 64: CH + O <=> CO + H */
    /*eqcon[63] *= 1;  */

    /*reaction 65: CH + O2 <=> HCO + O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: CH + H2O <=> CH2O + H */
    /*eqcon[65] *= 1;  */

    /*reaction 67: CH + CO2 <=> HCO + CO */
    /*eqcon[66] *= 1;  */

    /*reaction 68: CH3O + H <=> CH2O + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: CH3O + H <=> S-CH2 + H2O */
    /*eqcon[68] *= 1;  */

    /*reaction 70: CH3O + OH <=> CH2O + H2O */
    /*eqcon[69] *= 1;  */

    /*reaction 71: CH3O + O <=> OH + CH2O */
    /*eqcon[70] *= 1;  */

    /*reaction 72: CH3O + O2 <=> CH2O + HO2 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: CH3O + M <=> CH2O + H + M */
    eqcon[72] *= 1e-06; 

    /*reaction 74: C2H6 + H <=> C2H5 + H2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: C2H6 + O <=> C2H5 + OH */
    /*eqcon[74] *= 1;  */

    /*reaction 76: C2H6 + OH <=> C2H5 + H2O */
    /*eqcon[75] *= 1;  */

    /*reaction 77: C2H6 + CH3 <=> C2H5 + CH4 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: C2H6 (+M) <=> C2H5 + H (+M) */
    eqcon[77] *= 1e-06; 

    /*reaction 79: C2H5 + H <=> C2H4 + H2 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: C2H5 + O <=> C2H4 + OH */
    /*eqcon[79] *= 1;  */

    /*reaction 81: C2H5 + O <=> CH3 + CH2O */
    /*eqcon[80] *= 1;  */

    /*reaction 82: C2H5 + O2 <=> C2H4 + HO2 */
    /*eqcon[81] *= 1;  */

    /*reaction 83: C2H5 (+M) <=> C2H4 + H (+M) */
    eqcon[82] *= 1e-06; 

    /*reaction 84: C2H4 + H <=> C2H3 + H2 */
    /*eqcon[83] *= 1;  */

    /*reaction 85: C2H4 + OH <=> C2H3 + H2O */
    /*eqcon[84] *= 1;  */

    /*reaction 86: C2H4 + O <=> CH3 + HCO */
    /*eqcon[85] *= 1;  */

    /*reaction 87: C2H4 + O <=> CH2CHO + H */
    /*eqcon[86] *= 1;  */

    /*reaction 88: 2 C2H4 <=> C2H3 + C2H5 */
    /*eqcon[87] *= 1;  */

    /*reaction 89: C2H4 + O2 <=> C2H3 + HO2 */
    /*eqcon[88] *= 1;  */

    /*reaction 90: C2H4 + HO2 <=> C2H4O + OH */
    /*eqcon[89] *= 1;  */

    /*reaction 91: C2H4O + HO2 <=> CH3 + CO + H2O2 */
    eqcon[90] *= 1e-06; 

    /*reaction 92: C2H4 + M <=> C2H3 + H + M */
    eqcon[91] *= 1e-06; 

    /*reaction 93: C2H4 + M <=> C2H2 + H2 + M */
    eqcon[92] *= 1e-06; 

    /*reaction 94: C2H3 + H <=> C2H2 + H2 */
    /*eqcon[93] *= 1;  */

    /*reaction 95: C2H3 (+M) <=> C2H2 + H (+M) */
    eqcon[94] *= 1e-06; 

    /*reaction 96: C2H3 + O2 <=> CH2O + HCO */
    /*eqcon[95] *= 1;  */

    /*reaction 97: C2H3 + O2 <=> CH2CHO + O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: C2H3 + O2 <=> C2H2 + HO2 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: CH2CHO <=> CH2CO + H */
    eqcon[98] *= 1e-06; 

    /*reaction 100: C2H2 + O <=> HCCO + H */
    /*eqcon[99] *= 1;  */

    /*reaction 101: C2H2 + O <=> T-CH2 + CO */
    /*eqcon[100] *= 1;  */

    /*reaction 102: C2H2 + O2 <=> CH2O + CO */
    /*eqcon[101] *= 1;  */

    /*reaction 103: C2H2 + OH <=> CH2CO + H */
    /*eqcon[102] *= 1;  */

    /*reaction 104: C2H2 + OH <=> C2H + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: CH2CO + H <=> CH3 + CO */
    /*eqcon[104] *= 1;  */

    /*reaction 106: CH2CO + O <=> T-CH2 + CO2 */
    /*eqcon[105] *= 1;  */

    /*reaction 107: CH2CO + O <=> HCCO + OH */
    /*eqcon[106] *= 1;  */

    /*reaction 108: CH2CO + CH3 <=> C2H5 + CO */
    /*eqcon[107] *= 1;  */

    /*reaction 109: HCCO + H <=> S-CH2 + CO */
    /*eqcon[108] *= 1;  */

    /*reaction 110: HCCO + OH <=> HCO + CO + H */
    eqcon[109] *= 1e-06; 

    /*reaction 111: HCCO + O <=> 2 CO + H */
    eqcon[110] *= 1e-06; 

    /*reaction 112: HCCO + O2 <=> 2 CO + OH */
    eqcon[111] *= 1e-06; 

    /*reaction 113: HCCO + O2 <=> CO2 + CO + H */
    eqcon[112] *= 1e-06; 

    /*reaction 114: C2H + OH <=> HCCO + H */
    /*eqcon[113] *= 1;  */

    /*reaction 115: C2H + O <=> CO + CH */
    /*eqcon[114] *= 1;  */

    /*reaction 116: C2H + O2 <=> HCCO + O */
    /*eqcon[115] *= 1;  */

    /*reaction 117: C2H + O2 <=> CH + CO2 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: C2H + O2 <=> HCO + CO */
    /*eqcon[117] *= 1;  */

    /*reaction 119: CH2OH + H <=> CH2O + H2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: CH2OH + H <=> CH3 + OH */
    /*eqcon[119] *= 1;  */

    /*reaction 121: CH2OH + OH <=> CH2O + H2O */
    /*eqcon[120] *= 1;  */

    /*reaction 122: CH2OH + O2 <=> CH2O + HO2 */
    /*eqcon[121] *= 1;  */

    /*reaction 123: CH2OH + M <=> CH2O + H + M */
    eqcon[122] *= 1e-06; 

    /*reaction 124: CH3O + M <=> CH2OH + M */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH2CO + OH <=> CH2OH + CO */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH3OH + OH <=> CH2OH + H2O */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH3OH + OH <=> CH3O + H2O */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH3OH + H <=> CH2OH + H2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH3OH + H <=> CH3O + H2 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH3OH + O <=> CH2OH + OH */
    /*eqcon[129] *= 1;  */

    /*reaction 131: CH3OH + HO2 <=> CH2OH + H2O2 */
    /*eqcon[130] *= 1;  */

    /*reaction 132: CH3OH + O2 <=> CH2OH + HO2 */
    /*eqcon[131] *= 1;  */

    /*reaction 133: C3H4 + O <=> C2H4 + CO */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH3 + C2H2 <=> C3H4 + H */
    /*eqcon[133] *= 1;  */

    /*reaction 135: C3H4 + O <=> HCCO + CH3 */
    /*eqcon[134] *= 1;  */

    /*reaction 136: C3H3 + H (+M) <=> C3H4 (+M) */
    eqcon[135] *= 1e+06; 

    /*reaction 137: C3H3 + HO2 <=> C3H4 + O2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: C3H4 + OH <=> C3H3 + H2O */
    /*eqcon[137] *= 1;  */

    /*reaction 139: C3H3 + O2 <=> CH2CO + HCO */
    /*eqcon[138] *= 1;  */

    /*reaction 140: C3H4 + H (+M) <=> C3H5 (+M) */
    eqcon[139] *= 1e+06; 

    /*reaction 141: C3H5 + H <=> C3H4 + H2 */
    /*eqcon[140] *= 1;  */

    /*reaction 142: C3H5 + O2 <=> C3H4 + HO2 */
    /*eqcon[141] *= 1;  */

    /*reaction 143: C3H5 + CH3 <=> C3H4 + CH4 */
    /*eqcon[142] *= 1;  */

    /*reaction 144: C2H2 + CH3 (+M) <=> C3H5 (+M) */
    eqcon[143] *= 1e+06; 

    /*reaction 145: C3H5 + OH <=> C3H4 + H2O */
    /*eqcon[144] *= 1;  */

    /*reaction 146: C3H3 + HCO <=> C3H4 + CO */
    /*eqcon[145] *= 1;  */

    /*reaction 147: C3H3 + HO2 <=> OH + CO + C2H3 */
    eqcon[146] *= 1e-06; 

    /*reaction 148: C3H4 + O2 <=> CH3 + HCO + CO */
    eqcon[147] *= 1e-06; 

    /*reaction 149: C3H6 + O <=> C2H5 + HCO */
    /*eqcon[148] *= 1;  */

    /*reaction 150: C3H6 + OH <=> C3H5 + H2O */
    /*eqcon[149] *= 1;  */

    /*reaction 151: C3H6 + O <=> CH2CO + CH3 + H */
    eqcon[150] *= 1e-06; 

    /*reaction 152: C3H6 + H <=> C3H5 + H2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: C3H5 + H (+M) <=> C3H6 (+M) */
    eqcon[152] *= 1e+06; 

    /*reaction 154: C3H5 + HO2 <=> C3H6 + O2 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: C3H5 + HO2 <=> OH + C2H3 + CH2O */
    eqcon[154] *= 1e-06; 

    /*reaction 156: C2H3 + CH3 (+M) <=> C3H6 (+M) */
    eqcon[155] *= 1e+06; 

    /*reaction 157: C3H6 + H <=> C2H4 + CH3 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: CH3 + C2H3 <=> C3H5 + H */
    /*eqcon[157] *= 1;  */

    /*reaction 159: C3H8 (+M) <=> CH3 + C2H5 (+M) */
    eqcon[158] *= 1e-06; 

    /*reaction 160: C3H8 + O2 <=> I-C3H7 + HO2 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: C3H8 + O2 <=> N-C3H7 + HO2 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: C3H8 + H <=> I-C3H7 + H2 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: C3H8 + H <=> N-C3H7 + H2 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: C3H8 + O <=> I-C3H7 + OH */
    /*eqcon[163] *= 1;  */

    /*reaction 165: C3H8 + O <=> N-C3H7 + OH */
    /*eqcon[164] *= 1;  */

    /*reaction 166: C3H8 + OH <=> N-C3H7 + H2O */
    /*eqcon[165] *= 1;  */

    /*reaction 167: C3H8 + OH <=> I-C3H7 + H2O */
    /*eqcon[166] *= 1;  */

    /*reaction 168: C3H8 + HO2 <=> I-C3H7 + H2O2 */
    /*eqcon[167] *= 1;  */

    /*reaction 169: C3H8 + HO2 <=> N-C3H7 + H2O2 */
    /*eqcon[168] *= 1;  */

    /*reaction 170: I-C3H7 + C3H8 <=> N-C3H7 + C3H8 */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C3H6 + H (+M) <=> I-C3H7 (+M) */
    eqcon[170] *= 1e+06; 

    /*reaction 172: I-C3H7 + O2 <=> C3H6 + HO2 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: N-C3H7 (+M) <=> CH3 + C2H4 (+M) */
    eqcon[172] *= 1e-06; 

    /*reaction 174: H + C3H6 (+M) <=> N-C3H7 (+M) */
    eqcon[173] *= 1e+06; 

    /*reaction 175: N-C3H7 + O2 <=> C3H6 + HO2 */
    /*eqcon[174] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[39]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> OH + O */
    /*eqcon[0] *= 1;  */

    /*reaction 2: H2 + O <=> OH + H */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */

    /*reaction 4: H2O + O <=> 2 OH */
    /*eqcon[3] *= 1;  */

    /*reaction 5: H + O + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2 + O2 <=> 2 OH */
    /*eqcon[5] *= 1;  */

    /*reaction 7: 2 H + M <=> H2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + OH + M <=> H2O + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: 2 O + M <=> O2 + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 (+M) <=> HO2 (+M) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: O + OH + M <=> HO2 + M */
    eqcon[10] *= 1e+06; 

    /*reaction 12: HO2 + H <=> 2 OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + H <=> H2 + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + H <=> H2O + O */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + O <=> OH + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[16] *= 1e+06; 

    /*reaction 18: 2 HO2 <=> H2O2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + H <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + H <=> H2O + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2O2 + O <=> HO2 + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: CO + OH <=> CO2 + H */
    /*eqcon[22] *= 1;  */

    /*reaction 24: CO + HO2 <=> CO2 + OH */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HCO + M <=> CO + H + M */
    eqcon[24] *= 1e-06; 

    /*reaction 26: HCO + H <=> CO + H2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: HCO + O <=> CO + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: HCO + O <=> CO2 + H */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HCO + OH <=> CO + H2O */
    /*eqcon[28] *= 1;  */

    /*reaction 30: HCO + O2 <=> CO + HO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CH2O + M <=> HCO + H + M */
    eqcon[30] *= 1e-06; 

    /*reaction 32: CH2O + H <=> HCO + H2 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: CH2O + O <=> HCO + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: CH2O + OH <=> HCO + H2O */
    /*eqcon[33] *= 1;  */

    /*reaction 35: CH4 + H <=> H2 + CH3 */
    /*eqcon[34] *= 1;  */

    /*reaction 36: CH4 + OH <=> H2O + CH3 */
    /*eqcon[35] *= 1;  */

    /*reaction 37: CH4 + O <=> CH3 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: CH4 + O2 <=> CH3 + HO2 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: CH4 + HO2 <=> CH3 + H2O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: CH3 + H <=> T-CH2 + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: CH3 + H <=> S-CH2 + H2 */
    /*eqcon[40] *= 1;  */

    /*reaction 42: CH3 + OH <=> S-CH2 + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: CH3 + O <=> CH2O + H */
    /*eqcon[42] *= 1;  */

    /*reaction 44: CH3 + T-CH2 <=> C2H4 + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: CH3 + HO2 <=> CH3O + OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: CH3 + O2 <=> CH2O + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: CH3 + O2 <=> CH3O + O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: 2 CH3 <=> C2H4 + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: 2 CH3 <=> C2H5 + H */
    /*eqcon[48] *= 1;  */

    /*reaction 50: CH3 + H (+M) <=> CH4 (+M) */
    eqcon[49] *= 1e+06; 

    /*reaction 51: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[50] *= 1e+06; 

    /*reaction 52: S-CH2 + OH <=> CH2O + H */
    /*eqcon[51] *= 1;  */

    /*reaction 53: S-CH2 + O2 <=> CO + OH + H */
    eqcon[52] *= 1e-06; 

    /*reaction 54: S-CH2 + CO2 <=> CO + CH2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: S-CH2 + M <=> T-CH2 + M */
    /*eqcon[54] *= 1;  */

    /*reaction 56: T-CH2 + H <=> CH + H2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: T-CH2 + OH <=> CH2O + H */
    /*eqcon[56] *= 1;  */

    /*reaction 58: T-CH2 + OH <=> CH + H2O */
    /*eqcon[57] *= 1;  */

    /*reaction 59: T-CH2 + O <=> CO + 2 H */
    eqcon[58] *= 1e-06; 

    /*reaction 60: T-CH2 + O <=> CO + H2 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: T-CH2 + O2 <=> CO2 + H2 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: T-CH2 + O2 <=> CO + OH + H */
    eqcon[61] *= 1e-06; 

    /*reaction 63: 2 T-CH2 <=> C2H2 + 2 H */
    eqcon[62] *= 1e-06; 

    /*reaction 64: CH + O <=> CO + H */
    /*eqcon[63] *= 1;  */

    /*reaction 65: CH + O2 <=> HCO + O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: CH + H2O <=> CH2O + H */
    /*eqcon[65] *= 1;  */

    /*reaction 67: CH + CO2 <=> HCO + CO */
    /*eqcon[66] *= 1;  */

    /*reaction 68: CH3O + H <=> CH2O + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: CH3O + H <=> S-CH2 + H2O */
    /*eqcon[68] *= 1;  */

    /*reaction 70: CH3O + OH <=> CH2O + H2O */
    /*eqcon[69] *= 1;  */

    /*reaction 71: CH3O + O <=> OH + CH2O */
    /*eqcon[70] *= 1;  */

    /*reaction 72: CH3O + O2 <=> CH2O + HO2 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: CH3O + M <=> CH2O + H + M */
    eqcon[72] *= 1e-06; 

    /*reaction 74: C2H6 + H <=> C2H5 + H2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: C2H6 + O <=> C2H5 + OH */
    /*eqcon[74] *= 1;  */

    /*reaction 76: C2H6 + OH <=> C2H5 + H2O */
    /*eqcon[75] *= 1;  */

    /*reaction 77: C2H6 + CH3 <=> C2H5 + CH4 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: C2H6 (+M) <=> C2H5 + H (+M) */
    eqcon[77] *= 1e-06; 

    /*reaction 79: C2H5 + H <=> C2H4 + H2 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: C2H5 + O <=> C2H4 + OH */
    /*eqcon[79] *= 1;  */

    /*reaction 81: C2H5 + O <=> CH3 + CH2O */
    /*eqcon[80] *= 1;  */

    /*reaction 82: C2H5 + O2 <=> C2H4 + HO2 */
    /*eqcon[81] *= 1;  */

    /*reaction 83: C2H5 (+M) <=> C2H4 + H (+M) */
    eqcon[82] *= 1e-06; 

    /*reaction 84: C2H4 + H <=> C2H3 + H2 */
    /*eqcon[83] *= 1;  */

    /*reaction 85: C2H4 + OH <=> C2H3 + H2O */
    /*eqcon[84] *= 1;  */

    /*reaction 86: C2H4 + O <=> CH3 + HCO */
    /*eqcon[85] *= 1;  */

    /*reaction 87: C2H4 + O <=> CH2CHO + H */
    /*eqcon[86] *= 1;  */

    /*reaction 88: 2 C2H4 <=> C2H3 + C2H5 */
    /*eqcon[87] *= 1;  */

    /*reaction 89: C2H4 + O2 <=> C2H3 + HO2 */
    /*eqcon[88] *= 1;  */

    /*reaction 90: C2H4 + HO2 <=> C2H4O + OH */
    /*eqcon[89] *= 1;  */

    /*reaction 91: C2H4O + HO2 <=> CH3 + CO + H2O2 */
    eqcon[90] *= 1e-06; 

    /*reaction 92: C2H4 + M <=> C2H3 + H + M */
    eqcon[91] *= 1e-06; 

    /*reaction 93: C2H4 + M <=> C2H2 + H2 + M */
    eqcon[92] *= 1e-06; 

    /*reaction 94: C2H3 + H <=> C2H2 + H2 */
    /*eqcon[93] *= 1;  */

    /*reaction 95: C2H3 (+M) <=> C2H2 + H (+M) */
    eqcon[94] *= 1e-06; 

    /*reaction 96: C2H3 + O2 <=> CH2O + HCO */
    /*eqcon[95] *= 1;  */

    /*reaction 97: C2H3 + O2 <=> CH2CHO + O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: C2H3 + O2 <=> C2H2 + HO2 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: CH2CHO <=> CH2CO + H */
    eqcon[98] *= 1e-06; 

    /*reaction 100: C2H2 + O <=> HCCO + H */
    /*eqcon[99] *= 1;  */

    /*reaction 101: C2H2 + O <=> T-CH2 + CO */
    /*eqcon[100] *= 1;  */

    /*reaction 102: C2H2 + O2 <=> CH2O + CO */
    /*eqcon[101] *= 1;  */

    /*reaction 103: C2H2 + OH <=> CH2CO + H */
    /*eqcon[102] *= 1;  */

    /*reaction 104: C2H2 + OH <=> C2H + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: CH2CO + H <=> CH3 + CO */
    /*eqcon[104] *= 1;  */

    /*reaction 106: CH2CO + O <=> T-CH2 + CO2 */
    /*eqcon[105] *= 1;  */

    /*reaction 107: CH2CO + O <=> HCCO + OH */
    /*eqcon[106] *= 1;  */

    /*reaction 108: CH2CO + CH3 <=> C2H5 + CO */
    /*eqcon[107] *= 1;  */

    /*reaction 109: HCCO + H <=> S-CH2 + CO */
    /*eqcon[108] *= 1;  */

    /*reaction 110: HCCO + OH <=> HCO + CO + H */
    eqcon[109] *= 1e-06; 

    /*reaction 111: HCCO + O <=> 2 CO + H */
    eqcon[110] *= 1e-06; 

    /*reaction 112: HCCO + O2 <=> 2 CO + OH */
    eqcon[111] *= 1e-06; 

    /*reaction 113: HCCO + O2 <=> CO2 + CO + H */
    eqcon[112] *= 1e-06; 

    /*reaction 114: C2H + OH <=> HCCO + H */
    /*eqcon[113] *= 1;  */

    /*reaction 115: C2H + O <=> CO + CH */
    /*eqcon[114] *= 1;  */

    /*reaction 116: C2H + O2 <=> HCCO + O */
    /*eqcon[115] *= 1;  */

    /*reaction 117: C2H + O2 <=> CH + CO2 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: C2H + O2 <=> HCO + CO */
    /*eqcon[117] *= 1;  */

    /*reaction 119: CH2OH + H <=> CH2O + H2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: CH2OH + H <=> CH3 + OH */
    /*eqcon[119] *= 1;  */

    /*reaction 121: CH2OH + OH <=> CH2O + H2O */
    /*eqcon[120] *= 1;  */

    /*reaction 122: CH2OH + O2 <=> CH2O + HO2 */
    /*eqcon[121] *= 1;  */

    /*reaction 123: CH2OH + M <=> CH2O + H + M */
    eqcon[122] *= 1e-06; 

    /*reaction 124: CH3O + M <=> CH2OH + M */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH2CO + OH <=> CH2OH + CO */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH3OH + OH <=> CH2OH + H2O */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH3OH + OH <=> CH3O + H2O */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH3OH + H <=> CH2OH + H2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH3OH + H <=> CH3O + H2 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH3OH + O <=> CH2OH + OH */
    /*eqcon[129] *= 1;  */

    /*reaction 131: CH3OH + HO2 <=> CH2OH + H2O2 */
    /*eqcon[130] *= 1;  */

    /*reaction 132: CH3OH + O2 <=> CH2OH + HO2 */
    /*eqcon[131] *= 1;  */

    /*reaction 133: C3H4 + O <=> C2H4 + CO */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH3 + C2H2 <=> C3H4 + H */
    /*eqcon[133] *= 1;  */

    /*reaction 135: C3H4 + O <=> HCCO + CH3 */
    /*eqcon[134] *= 1;  */

    /*reaction 136: C3H3 + H (+M) <=> C3H4 (+M) */
    eqcon[135] *= 1e+06; 

    /*reaction 137: C3H3 + HO2 <=> C3H4 + O2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: C3H4 + OH <=> C3H3 + H2O */
    /*eqcon[137] *= 1;  */

    /*reaction 139: C3H3 + O2 <=> CH2CO + HCO */
    /*eqcon[138] *= 1;  */

    /*reaction 140: C3H4 + H (+M) <=> C3H5 (+M) */
    eqcon[139] *= 1e+06; 

    /*reaction 141: C3H5 + H <=> C3H4 + H2 */
    /*eqcon[140] *= 1;  */

    /*reaction 142: C3H5 + O2 <=> C3H4 + HO2 */
    /*eqcon[141] *= 1;  */

    /*reaction 143: C3H5 + CH3 <=> C3H4 + CH4 */
    /*eqcon[142] *= 1;  */

    /*reaction 144: C2H2 + CH3 (+M) <=> C3H5 (+M) */
    eqcon[143] *= 1e+06; 

    /*reaction 145: C3H5 + OH <=> C3H4 + H2O */
    /*eqcon[144] *= 1;  */

    /*reaction 146: C3H3 + HCO <=> C3H4 + CO */
    /*eqcon[145] *= 1;  */

    /*reaction 147: C3H3 + HO2 <=> OH + CO + C2H3 */
    eqcon[146] *= 1e-06; 

    /*reaction 148: C3H4 + O2 <=> CH3 + HCO + CO */
    eqcon[147] *= 1e-06; 

    /*reaction 149: C3H6 + O <=> C2H5 + HCO */
    /*eqcon[148] *= 1;  */

    /*reaction 150: C3H6 + OH <=> C3H5 + H2O */
    /*eqcon[149] *= 1;  */

    /*reaction 151: C3H6 + O <=> CH2CO + CH3 + H */
    eqcon[150] *= 1e-06; 

    /*reaction 152: C3H6 + H <=> C3H5 + H2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: C3H5 + H (+M) <=> C3H6 (+M) */
    eqcon[152] *= 1e+06; 

    /*reaction 154: C3H5 + HO2 <=> C3H6 + O2 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: C3H5 + HO2 <=> OH + C2H3 + CH2O */
    eqcon[154] *= 1e-06; 

    /*reaction 156: C2H3 + CH3 (+M) <=> C3H6 (+M) */
    eqcon[155] *= 1e+06; 

    /*reaction 157: C3H6 + H <=> C2H4 + CH3 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: CH3 + C2H3 <=> C3H5 + H */
    /*eqcon[157] *= 1;  */

    /*reaction 159: C3H8 (+M) <=> CH3 + C2H5 (+M) */
    eqcon[158] *= 1e-06; 

    /*reaction 160: C3H8 + O2 <=> I-C3H7 + HO2 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: C3H8 + O2 <=> N-C3H7 + HO2 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: C3H8 + H <=> I-C3H7 + H2 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: C3H8 + H <=> N-C3H7 + H2 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: C3H8 + O <=> I-C3H7 + OH */
    /*eqcon[163] *= 1;  */

    /*reaction 165: C3H8 + O <=> N-C3H7 + OH */
    /*eqcon[164] *= 1;  */

    /*reaction 166: C3H8 + OH <=> N-C3H7 + H2O */
    /*eqcon[165] *= 1;  */

    /*reaction 167: C3H8 + OH <=> I-C3H7 + H2O */
    /*eqcon[166] *= 1;  */

    /*reaction 168: C3H8 + HO2 <=> I-C3H7 + H2O2 */
    /*eqcon[167] *= 1;  */

    /*reaction 169: C3H8 + HO2 <=> N-C3H7 + H2O2 */
    /*eqcon[168] *= 1;  */

    /*reaction 170: I-C3H7 + C3H8 <=> N-C3H7 + C3H8 */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C3H6 + H (+M) <=> I-C3H7 (+M) */
    eqcon[170] *= 1e+06; 

    /*reaction 172: I-C3H7 + O2 <=> C3H6 + HO2 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: N-C3H7 (+M) <=> CH3 + C2H4 (+M) */
    eqcon[172] *= 1e-06; 

    /*reaction 174: H + C3H6 (+M) <=> N-C3H7 (+M) */
    eqcon[173] *= 1e+06; 

    /*reaction 175: N-C3H7 + O2 <=> C3H6 + HO2 */
    /*eqcon[174] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[39]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> OH + O */
    /*eqcon[0] *= 1;  */

    /*reaction 2: H2 + O <=> OH + H */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */

    /*reaction 4: H2O + O <=> 2 OH */
    /*eqcon[3] *= 1;  */

    /*reaction 5: H + O + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2 + O2 <=> 2 OH */
    /*eqcon[5] *= 1;  */

    /*reaction 7: 2 H + M <=> H2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + OH + M <=> H2O + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: 2 O + M <=> O2 + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 (+M) <=> HO2 (+M) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: O + OH + M <=> HO2 + M */
    eqcon[10] *= 1e+06; 

    /*reaction 12: HO2 + H <=> 2 OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + H <=> H2 + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + H <=> H2O + O */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + O <=> OH + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[16] *= 1e+06; 

    /*reaction 18: 2 HO2 <=> H2O2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + H <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + H <=> H2O + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2O2 + O <=> HO2 + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: CO + OH <=> CO2 + H */
    /*eqcon[22] *= 1;  */

    /*reaction 24: CO + HO2 <=> CO2 + OH */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HCO + M <=> CO + H + M */
    eqcon[24] *= 1e-06; 

    /*reaction 26: HCO + H <=> CO + H2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: HCO + O <=> CO + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: HCO + O <=> CO2 + H */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HCO + OH <=> CO + H2O */
    /*eqcon[28] *= 1;  */

    /*reaction 30: HCO + O2 <=> CO + HO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CH2O + M <=> HCO + H + M */
    eqcon[30] *= 1e-06; 

    /*reaction 32: CH2O + H <=> HCO + H2 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: CH2O + O <=> HCO + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: CH2O + OH <=> HCO + H2O */
    /*eqcon[33] *= 1;  */

    /*reaction 35: CH4 + H <=> H2 + CH3 */
    /*eqcon[34] *= 1;  */

    /*reaction 36: CH4 + OH <=> H2O + CH3 */
    /*eqcon[35] *= 1;  */

    /*reaction 37: CH4 + O <=> CH3 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: CH4 + O2 <=> CH3 + HO2 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: CH4 + HO2 <=> CH3 + H2O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: CH3 + H <=> T-CH2 + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: CH3 + H <=> S-CH2 + H2 */
    /*eqcon[40] *= 1;  */

    /*reaction 42: CH3 + OH <=> S-CH2 + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: CH3 + O <=> CH2O + H */
    /*eqcon[42] *= 1;  */

    /*reaction 44: CH3 + T-CH2 <=> C2H4 + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: CH3 + HO2 <=> CH3O + OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: CH3 + O2 <=> CH2O + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: CH3 + O2 <=> CH3O + O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: 2 CH3 <=> C2H4 + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: 2 CH3 <=> C2H5 + H */
    /*eqcon[48] *= 1;  */

    /*reaction 50: CH3 + H (+M) <=> CH4 (+M) */
    eqcon[49] *= 1e+06; 

    /*reaction 51: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[50] *= 1e+06; 

    /*reaction 52: S-CH2 + OH <=> CH2O + H */
    /*eqcon[51] *= 1;  */

    /*reaction 53: S-CH2 + O2 <=> CO + OH + H */
    eqcon[52] *= 1e-06; 

    /*reaction 54: S-CH2 + CO2 <=> CO + CH2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: S-CH2 + M <=> T-CH2 + M */
    /*eqcon[54] *= 1;  */

    /*reaction 56: T-CH2 + H <=> CH + H2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: T-CH2 + OH <=> CH2O + H */
    /*eqcon[56] *= 1;  */

    /*reaction 58: T-CH2 + OH <=> CH + H2O */
    /*eqcon[57] *= 1;  */

    /*reaction 59: T-CH2 + O <=> CO + 2 H */
    eqcon[58] *= 1e-06; 

    /*reaction 60: T-CH2 + O <=> CO + H2 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: T-CH2 + O2 <=> CO2 + H2 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: T-CH2 + O2 <=> CO + OH + H */
    eqcon[61] *= 1e-06; 

    /*reaction 63: 2 T-CH2 <=> C2H2 + 2 H */
    eqcon[62] *= 1e-06; 

    /*reaction 64: CH + O <=> CO + H */
    /*eqcon[63] *= 1;  */

    /*reaction 65: CH + O2 <=> HCO + O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: CH + H2O <=> CH2O + H */
    /*eqcon[65] *= 1;  */

    /*reaction 67: CH + CO2 <=> HCO + CO */
    /*eqcon[66] *= 1;  */

    /*reaction 68: CH3O + H <=> CH2O + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: CH3O + H <=> S-CH2 + H2O */
    /*eqcon[68] *= 1;  */

    /*reaction 70: CH3O + OH <=> CH2O + H2O */
    /*eqcon[69] *= 1;  */

    /*reaction 71: CH3O + O <=> OH + CH2O */
    /*eqcon[70] *= 1;  */

    /*reaction 72: CH3O + O2 <=> CH2O + HO2 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: CH3O + M <=> CH2O + H + M */
    eqcon[72] *= 1e-06; 

    /*reaction 74: C2H6 + H <=> C2H5 + H2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: C2H6 + O <=> C2H5 + OH */
    /*eqcon[74] *= 1;  */

    /*reaction 76: C2H6 + OH <=> C2H5 + H2O */
    /*eqcon[75] *= 1;  */

    /*reaction 77: C2H6 + CH3 <=> C2H5 + CH4 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: C2H6 (+M) <=> C2H5 + H (+M) */
    eqcon[77] *= 1e-06; 

    /*reaction 79: C2H5 + H <=> C2H4 + H2 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: C2H5 + O <=> C2H4 + OH */
    /*eqcon[79] *= 1;  */

    /*reaction 81: C2H5 + O <=> CH3 + CH2O */
    /*eqcon[80] *= 1;  */

    /*reaction 82: C2H5 + O2 <=> C2H4 + HO2 */
    /*eqcon[81] *= 1;  */

    /*reaction 83: C2H5 (+M) <=> C2H4 + H (+M) */
    eqcon[82] *= 1e-06; 

    /*reaction 84: C2H4 + H <=> C2H3 + H2 */
    /*eqcon[83] *= 1;  */

    /*reaction 85: C2H4 + OH <=> C2H3 + H2O */
    /*eqcon[84] *= 1;  */

    /*reaction 86: C2H4 + O <=> CH3 + HCO */
    /*eqcon[85] *= 1;  */

    /*reaction 87: C2H4 + O <=> CH2CHO + H */
    /*eqcon[86] *= 1;  */

    /*reaction 88: 2 C2H4 <=> C2H3 + C2H5 */
    /*eqcon[87] *= 1;  */

    /*reaction 89: C2H4 + O2 <=> C2H3 + HO2 */
    /*eqcon[88] *= 1;  */

    /*reaction 90: C2H4 + HO2 <=> C2H4O + OH */
    /*eqcon[89] *= 1;  */

    /*reaction 91: C2H4O + HO2 <=> CH3 + CO + H2O2 */
    eqcon[90] *= 1e-06; 

    /*reaction 92: C2H4 + M <=> C2H3 + H + M */
    eqcon[91] *= 1e-06; 

    /*reaction 93: C2H4 + M <=> C2H2 + H2 + M */
    eqcon[92] *= 1e-06; 

    /*reaction 94: C2H3 + H <=> C2H2 + H2 */
    /*eqcon[93] *= 1;  */

    /*reaction 95: C2H3 (+M) <=> C2H2 + H (+M) */
    eqcon[94] *= 1e-06; 

    /*reaction 96: C2H3 + O2 <=> CH2O + HCO */
    /*eqcon[95] *= 1;  */

    /*reaction 97: C2H3 + O2 <=> CH2CHO + O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: C2H3 + O2 <=> C2H2 + HO2 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: CH2CHO <=> CH2CO + H */
    eqcon[98] *= 1e-06; 

    /*reaction 100: C2H2 + O <=> HCCO + H */
    /*eqcon[99] *= 1;  */

    /*reaction 101: C2H2 + O <=> T-CH2 + CO */
    /*eqcon[100] *= 1;  */

    /*reaction 102: C2H2 + O2 <=> CH2O + CO */
    /*eqcon[101] *= 1;  */

    /*reaction 103: C2H2 + OH <=> CH2CO + H */
    /*eqcon[102] *= 1;  */

    /*reaction 104: C2H2 + OH <=> C2H + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: CH2CO + H <=> CH3 + CO */
    /*eqcon[104] *= 1;  */

    /*reaction 106: CH2CO + O <=> T-CH2 + CO2 */
    /*eqcon[105] *= 1;  */

    /*reaction 107: CH2CO + O <=> HCCO + OH */
    /*eqcon[106] *= 1;  */

    /*reaction 108: CH2CO + CH3 <=> C2H5 + CO */
    /*eqcon[107] *= 1;  */

    /*reaction 109: HCCO + H <=> S-CH2 + CO */
    /*eqcon[108] *= 1;  */

    /*reaction 110: HCCO + OH <=> HCO + CO + H */
    eqcon[109] *= 1e-06; 

    /*reaction 111: HCCO + O <=> 2 CO + H */
    eqcon[110] *= 1e-06; 

    /*reaction 112: HCCO + O2 <=> 2 CO + OH */
    eqcon[111] *= 1e-06; 

    /*reaction 113: HCCO + O2 <=> CO2 + CO + H */
    eqcon[112] *= 1e-06; 

    /*reaction 114: C2H + OH <=> HCCO + H */
    /*eqcon[113] *= 1;  */

    /*reaction 115: C2H + O <=> CO + CH */
    /*eqcon[114] *= 1;  */

    /*reaction 116: C2H + O2 <=> HCCO + O */
    /*eqcon[115] *= 1;  */

    /*reaction 117: C2H + O2 <=> CH + CO2 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: C2H + O2 <=> HCO + CO */
    /*eqcon[117] *= 1;  */

    /*reaction 119: CH2OH + H <=> CH2O + H2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: CH2OH + H <=> CH3 + OH */
    /*eqcon[119] *= 1;  */

    /*reaction 121: CH2OH + OH <=> CH2O + H2O */
    /*eqcon[120] *= 1;  */

    /*reaction 122: CH2OH + O2 <=> CH2O + HO2 */
    /*eqcon[121] *= 1;  */

    /*reaction 123: CH2OH + M <=> CH2O + H + M */
    eqcon[122] *= 1e-06; 

    /*reaction 124: CH3O + M <=> CH2OH + M */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH2CO + OH <=> CH2OH + CO */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH3OH + OH <=> CH2OH + H2O */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH3OH + OH <=> CH3O + H2O */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH3OH + H <=> CH2OH + H2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH3OH + H <=> CH3O + H2 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH3OH + O <=> CH2OH + OH */
    /*eqcon[129] *= 1;  */

    /*reaction 131: CH3OH + HO2 <=> CH2OH + H2O2 */
    /*eqcon[130] *= 1;  */

    /*reaction 132: CH3OH + O2 <=> CH2OH + HO2 */
    /*eqcon[131] *= 1;  */

    /*reaction 133: C3H4 + O <=> C2H4 + CO */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH3 + C2H2 <=> C3H4 + H */
    /*eqcon[133] *= 1;  */

    /*reaction 135: C3H4 + O <=> HCCO + CH3 */
    /*eqcon[134] *= 1;  */

    /*reaction 136: C3H3 + H (+M) <=> C3H4 (+M) */
    eqcon[135] *= 1e+06; 

    /*reaction 137: C3H3 + HO2 <=> C3H4 + O2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: C3H4 + OH <=> C3H3 + H2O */
    /*eqcon[137] *= 1;  */

    /*reaction 139: C3H3 + O2 <=> CH2CO + HCO */
    /*eqcon[138] *= 1;  */

    /*reaction 140: C3H4 + H (+M) <=> C3H5 (+M) */
    eqcon[139] *= 1e+06; 

    /*reaction 141: C3H5 + H <=> C3H4 + H2 */
    /*eqcon[140] *= 1;  */

    /*reaction 142: C3H5 + O2 <=> C3H4 + HO2 */
    /*eqcon[141] *= 1;  */

    /*reaction 143: C3H5 + CH3 <=> C3H4 + CH4 */
    /*eqcon[142] *= 1;  */

    /*reaction 144: C2H2 + CH3 (+M) <=> C3H5 (+M) */
    eqcon[143] *= 1e+06; 

    /*reaction 145: C3H5 + OH <=> C3H4 + H2O */
    /*eqcon[144] *= 1;  */

    /*reaction 146: C3H3 + HCO <=> C3H4 + CO */
    /*eqcon[145] *= 1;  */

    /*reaction 147: C3H3 + HO2 <=> OH + CO + C2H3 */
    eqcon[146] *= 1e-06; 

    /*reaction 148: C3H4 + O2 <=> CH3 + HCO + CO */
    eqcon[147] *= 1e-06; 

    /*reaction 149: C3H6 + O <=> C2H5 + HCO */
    /*eqcon[148] *= 1;  */

    /*reaction 150: C3H6 + OH <=> C3H5 + H2O */
    /*eqcon[149] *= 1;  */

    /*reaction 151: C3H6 + O <=> CH2CO + CH3 + H */
    eqcon[150] *= 1e-06; 

    /*reaction 152: C3H6 + H <=> C3H5 + H2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: C3H5 + H (+M) <=> C3H6 (+M) */
    eqcon[152] *= 1e+06; 

    /*reaction 154: C3H5 + HO2 <=> C3H6 + O2 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: C3H5 + HO2 <=> OH + C2H3 + CH2O */
    eqcon[154] *= 1e-06; 

    /*reaction 156: C2H3 + CH3 (+M) <=> C3H6 (+M) */
    eqcon[155] *= 1e+06; 

    /*reaction 157: C3H6 + H <=> C2H4 + CH3 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: CH3 + C2H3 <=> C3H5 + H */
    /*eqcon[157] *= 1;  */

    /*reaction 159: C3H8 (+M) <=> CH3 + C2H5 (+M) */
    eqcon[158] *= 1e-06; 

    /*reaction 160: C3H8 + O2 <=> I-C3H7 + HO2 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: C3H8 + O2 <=> N-C3H7 + HO2 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: C3H8 + H <=> I-C3H7 + H2 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: C3H8 + H <=> N-C3H7 + H2 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: C3H8 + O <=> I-C3H7 + OH */
    /*eqcon[163] *= 1;  */

    /*reaction 165: C3H8 + O <=> N-C3H7 + OH */
    /*eqcon[164] *= 1;  */

    /*reaction 166: C3H8 + OH <=> N-C3H7 + H2O */
    /*eqcon[165] *= 1;  */

    /*reaction 167: C3H8 + OH <=> I-C3H7 + H2O */
    /*eqcon[166] *= 1;  */

    /*reaction 168: C3H8 + HO2 <=> I-C3H7 + H2O2 */
    /*eqcon[167] *= 1;  */

    /*reaction 169: C3H8 + HO2 <=> N-C3H7 + H2O2 */
    /*eqcon[168] *= 1;  */

    /*reaction 170: I-C3H7 + C3H8 <=> N-C3H7 + C3H8 */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C3H6 + H (+M) <=> I-C3H7 (+M) */
    eqcon[170] *= 1e+06; 

    /*reaction 172: I-C3H7 + O2 <=> C3H6 + HO2 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: N-C3H7 (+M) <=> CH3 + C2H4 (+M) */
    eqcon[172] *= 1e-06; 

    /*reaction 174: H + C3H6 (+M) <=> N-C3H7 (+M) */
    eqcon[173] *= 1e+06; 

    /*reaction 175: N-C3H7 + O2 <=> C3H6 + HO2 */
    /*eqcon[174] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[39]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> OH + O */
    /*eqcon[0] *= 1;  */

    /*reaction 2: H2 + O <=> OH + H */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */

    /*reaction 4: H2O + O <=> 2 OH */
    /*eqcon[3] *= 1;  */

    /*reaction 5: H + O + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2 + O2 <=> 2 OH */
    /*eqcon[5] *= 1;  */

    /*reaction 7: 2 H + M <=> H2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + OH + M <=> H2O + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: 2 O + M <=> O2 + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 (+M) <=> HO2 (+M) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: O + OH + M <=> HO2 + M */
    eqcon[10] *= 1e+06; 

    /*reaction 12: HO2 + H <=> 2 OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + H <=> H2 + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + H <=> H2O + O */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + O <=> OH + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[16] *= 1e+06; 

    /*reaction 18: 2 HO2 <=> H2O2 + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + H <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + H <=> H2O + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2O2 + O <=> HO2 + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: CO + OH <=> CO2 + H */
    /*eqcon[22] *= 1;  */

    /*reaction 24: CO + HO2 <=> CO2 + OH */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HCO + M <=> CO + H + M */
    eqcon[24] *= 1e-06; 

    /*reaction 26: HCO + H <=> CO + H2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: HCO + O <=> CO + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: HCO + O <=> CO2 + H */
    /*eqcon[27] *= 1;  */

    /*reaction 29: HCO + OH <=> CO + H2O */
    /*eqcon[28] *= 1;  */

    /*reaction 30: HCO + O2 <=> CO + HO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CH2O + M <=> HCO + H + M */
    eqcon[30] *= 1e-06; 

    /*reaction 32: CH2O + H <=> HCO + H2 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: CH2O + O <=> HCO + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: CH2O + OH <=> HCO + H2O */
    /*eqcon[33] *= 1;  */

    /*reaction 35: CH4 + H <=> H2 + CH3 */
    /*eqcon[34] *= 1;  */

    /*reaction 36: CH4 + OH <=> H2O + CH3 */
    /*eqcon[35] *= 1;  */

    /*reaction 37: CH4 + O <=> CH3 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: CH4 + O2 <=> CH3 + HO2 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: CH4 + HO2 <=> CH3 + H2O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: CH3 + H <=> T-CH2 + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: CH3 + H <=> S-CH2 + H2 */
    /*eqcon[40] *= 1;  */

    /*reaction 42: CH3 + OH <=> S-CH2 + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: CH3 + O <=> CH2O + H */
    /*eqcon[42] *= 1;  */

    /*reaction 44: CH3 + T-CH2 <=> C2H4 + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: CH3 + HO2 <=> CH3O + OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: CH3 + O2 <=> CH2O + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: CH3 + O2 <=> CH3O + O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: 2 CH3 <=> C2H4 + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: 2 CH3 <=> C2H5 + H */
    /*eqcon[48] *= 1;  */

    /*reaction 50: CH3 + H (+M) <=> CH4 (+M) */
    eqcon[49] *= 1e+06; 

    /*reaction 51: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[50] *= 1e+06; 

    /*reaction 52: S-CH2 + OH <=> CH2O + H */
    /*eqcon[51] *= 1;  */

    /*reaction 53: S-CH2 + O2 <=> CO + OH + H */
    eqcon[52] *= 1e-06; 

    /*reaction 54: S-CH2 + CO2 <=> CO + CH2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: S-CH2 + M <=> T-CH2 + M */
    /*eqcon[54] *= 1;  */

    /*reaction 56: T-CH2 + H <=> CH + H2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: T-CH2 + OH <=> CH2O + H */
    /*eqcon[56] *= 1;  */

    /*reaction 58: T-CH2 + OH <=> CH + H2O */
    /*eqcon[57] *= 1;  */

    /*reaction 59: T-CH2 + O <=> CO + 2 H */
    eqcon[58] *= 1e-06; 

    /*reaction 60: T-CH2 + O <=> CO + H2 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: T-CH2 + O2 <=> CO2 + H2 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: T-CH2 + O2 <=> CO + OH + H */
    eqcon[61] *= 1e-06; 

    /*reaction 63: 2 T-CH2 <=> C2H2 + 2 H */
    eqcon[62] *= 1e-06; 

    /*reaction 64: CH + O <=> CO + H */
    /*eqcon[63] *= 1;  */

    /*reaction 65: CH + O2 <=> HCO + O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: CH + H2O <=> CH2O + H */
    /*eqcon[65] *= 1;  */

    /*reaction 67: CH + CO2 <=> HCO + CO */
    /*eqcon[66] *= 1;  */

    /*reaction 68: CH3O + H <=> CH2O + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: CH3O + H <=> S-CH2 + H2O */
    /*eqcon[68] *= 1;  */

    /*reaction 70: CH3O + OH <=> CH2O + H2O */
    /*eqcon[69] *= 1;  */

    /*reaction 71: CH3O + O <=> OH + CH2O */
    /*eqcon[70] *= 1;  */

    /*reaction 72: CH3O + O2 <=> CH2O + HO2 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: CH3O + M <=> CH2O + H + M */
    eqcon[72] *= 1e-06; 

    /*reaction 74: C2H6 + H <=> C2H5 + H2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: C2H6 + O <=> C2H5 + OH */
    /*eqcon[74] *= 1;  */

    /*reaction 76: C2H6 + OH <=> C2H5 + H2O */
    /*eqcon[75] *= 1;  */

    /*reaction 77: C2H6 + CH3 <=> C2H5 + CH4 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: C2H6 (+M) <=> C2H5 + H (+M) */
    eqcon[77] *= 1e-06; 

    /*reaction 79: C2H5 + H <=> C2H4 + H2 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: C2H5 + O <=> C2H4 + OH */
    /*eqcon[79] *= 1;  */

    /*reaction 81: C2H5 + O <=> CH3 + CH2O */
    /*eqcon[80] *= 1;  */

    /*reaction 82: C2H5 + O2 <=> C2H4 + HO2 */
    /*eqcon[81] *= 1;  */

    /*reaction 83: C2H5 (+M) <=> C2H4 + H (+M) */
    eqcon[82] *= 1e-06; 

    /*reaction 84: C2H4 + H <=> C2H3 + H2 */
    /*eqcon[83] *= 1;  */

    /*reaction 85: C2H4 + OH <=> C2H3 + H2O */
    /*eqcon[84] *= 1;  */

    /*reaction 86: C2H4 + O <=> CH3 + HCO */
    /*eqcon[85] *= 1;  */

    /*reaction 87: C2H4 + O <=> CH2CHO + H */
    /*eqcon[86] *= 1;  */

    /*reaction 88: 2 C2H4 <=> C2H3 + C2H5 */
    /*eqcon[87] *= 1;  */

    /*reaction 89: C2H4 + O2 <=> C2H3 + HO2 */
    /*eqcon[88] *= 1;  */

    /*reaction 90: C2H4 + HO2 <=> C2H4O + OH */
    /*eqcon[89] *= 1;  */

    /*reaction 91: C2H4O + HO2 <=> CH3 + CO + H2O2 */
    eqcon[90] *= 1e-06; 

    /*reaction 92: C2H4 + M <=> C2H3 + H + M */
    eqcon[91] *= 1e-06; 

    /*reaction 93: C2H4 + M <=> C2H2 + H2 + M */
    eqcon[92] *= 1e-06; 

    /*reaction 94: C2H3 + H <=> C2H2 + H2 */
    /*eqcon[93] *= 1;  */

    /*reaction 95: C2H3 (+M) <=> C2H2 + H (+M) */
    eqcon[94] *= 1e-06; 

    /*reaction 96: C2H3 + O2 <=> CH2O + HCO */
    /*eqcon[95] *= 1;  */

    /*reaction 97: C2H3 + O2 <=> CH2CHO + O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: C2H3 + O2 <=> C2H2 + HO2 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: CH2CHO <=> CH2CO + H */
    eqcon[98] *= 1e-06; 

    /*reaction 100: C2H2 + O <=> HCCO + H */
    /*eqcon[99] *= 1;  */

    /*reaction 101: C2H2 + O <=> T-CH2 + CO */
    /*eqcon[100] *= 1;  */

    /*reaction 102: C2H2 + O2 <=> CH2O + CO */
    /*eqcon[101] *= 1;  */

    /*reaction 103: C2H2 + OH <=> CH2CO + H */
    /*eqcon[102] *= 1;  */

    /*reaction 104: C2H2 + OH <=> C2H + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: CH2CO + H <=> CH3 + CO */
    /*eqcon[104] *= 1;  */

    /*reaction 106: CH2CO + O <=> T-CH2 + CO2 */
    /*eqcon[105] *= 1;  */

    /*reaction 107: CH2CO + O <=> HCCO + OH */
    /*eqcon[106] *= 1;  */

    /*reaction 108: CH2CO + CH3 <=> C2H5 + CO */
    /*eqcon[107] *= 1;  */

    /*reaction 109: HCCO + H <=> S-CH2 + CO */
    /*eqcon[108] *= 1;  */

    /*reaction 110: HCCO + OH <=> HCO + CO + H */
    eqcon[109] *= 1e-06; 

    /*reaction 111: HCCO + O <=> 2 CO + H */
    eqcon[110] *= 1e-06; 

    /*reaction 112: HCCO + O2 <=> 2 CO + OH */
    eqcon[111] *= 1e-06; 

    /*reaction 113: HCCO + O2 <=> CO2 + CO + H */
    eqcon[112] *= 1e-06; 

    /*reaction 114: C2H + OH <=> HCCO + H */
    /*eqcon[113] *= 1;  */

    /*reaction 115: C2H + O <=> CO + CH */
    /*eqcon[114] *= 1;  */

    /*reaction 116: C2H + O2 <=> HCCO + O */
    /*eqcon[115] *= 1;  */

    /*reaction 117: C2H + O2 <=> CH + CO2 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: C2H + O2 <=> HCO + CO */
    /*eqcon[117] *= 1;  */

    /*reaction 119: CH2OH + H <=> CH2O + H2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: CH2OH + H <=> CH3 + OH */
    /*eqcon[119] *= 1;  */

    /*reaction 121: CH2OH + OH <=> CH2O + H2O */
    /*eqcon[120] *= 1;  */

    /*reaction 122: CH2OH + O2 <=> CH2O + HO2 */
    /*eqcon[121] *= 1;  */

    /*reaction 123: CH2OH + M <=> CH2O + H + M */
    eqcon[122] *= 1e-06; 

    /*reaction 124: CH3O + M <=> CH2OH + M */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH2CO + OH <=> CH2OH + CO */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH3OH + OH <=> CH2OH + H2O */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH3OH + OH <=> CH3O + H2O */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH3OH + H <=> CH2OH + H2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH3OH + H <=> CH3O + H2 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH3OH + O <=> CH2OH + OH */
    /*eqcon[129] *= 1;  */

    /*reaction 131: CH3OH + HO2 <=> CH2OH + H2O2 */
    /*eqcon[130] *= 1;  */

    /*reaction 132: CH3OH + O2 <=> CH2OH + HO2 */
    /*eqcon[131] *= 1;  */

    /*reaction 133: C3H4 + O <=> C2H4 + CO */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH3 + C2H2 <=> C3H4 + H */
    /*eqcon[133] *= 1;  */

    /*reaction 135: C3H4 + O <=> HCCO + CH3 */
    /*eqcon[134] *= 1;  */

    /*reaction 136: C3H3 + H (+M) <=> C3H4 (+M) */
    eqcon[135] *= 1e+06; 

    /*reaction 137: C3H3 + HO2 <=> C3H4 + O2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: C3H4 + OH <=> C3H3 + H2O */
    /*eqcon[137] *= 1;  */

    /*reaction 139: C3H3 + O2 <=> CH2CO + HCO */
    /*eqcon[138] *= 1;  */

    /*reaction 140: C3H4 + H (+M) <=> C3H5 (+M) */
    eqcon[139] *= 1e+06; 

    /*reaction 141: C3H5 + H <=> C3H4 + H2 */
    /*eqcon[140] *= 1;  */

    /*reaction 142: C3H5 + O2 <=> C3H4 + HO2 */
    /*eqcon[141] *= 1;  */

    /*reaction 143: C3H5 + CH3 <=> C3H4 + CH4 */
    /*eqcon[142] *= 1;  */

    /*reaction 144: C2H2 + CH3 (+M) <=> C3H5 (+M) */
    eqcon[143] *= 1e+06; 

    /*reaction 145: C3H5 + OH <=> C3H4 + H2O */
    /*eqcon[144] *= 1;  */

    /*reaction 146: C3H3 + HCO <=> C3H4 + CO */
    /*eqcon[145] *= 1;  */

    /*reaction 147: C3H3 + HO2 <=> OH + CO + C2H3 */
    eqcon[146] *= 1e-06; 

    /*reaction 148: C3H4 + O2 <=> CH3 + HCO + CO */
    eqcon[147] *= 1e-06; 

    /*reaction 149: C3H6 + O <=> C2H5 + HCO */
    /*eqcon[148] *= 1;  */

    /*reaction 150: C3H6 + OH <=> C3H5 + H2O */
    /*eqcon[149] *= 1;  */

    /*reaction 151: C3H6 + O <=> CH2CO + CH3 + H */
    eqcon[150] *= 1e-06; 

    /*reaction 152: C3H6 + H <=> C3H5 + H2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: C3H5 + H (+M) <=> C3H6 (+M) */
    eqcon[152] *= 1e+06; 

    /*reaction 154: C3H5 + HO2 <=> C3H6 + O2 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: C3H5 + HO2 <=> OH + C2H3 + CH2O */
    eqcon[154] *= 1e-06; 

    /*reaction 156: C2H3 + CH3 (+M) <=> C3H6 (+M) */
    eqcon[155] *= 1e+06; 

    /*reaction 157: C3H6 + H <=> C2H4 + CH3 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: CH3 + C2H3 <=> C3H5 + H */
    /*eqcon[157] *= 1;  */

    /*reaction 159: C3H8 (+M) <=> CH3 + C2H5 (+M) */
    eqcon[158] *= 1e-06; 

    /*reaction 160: C3H8 + O2 <=> I-C3H7 + HO2 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: C3H8 + O2 <=> N-C3H7 + HO2 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: C3H8 + H <=> I-C3H7 + H2 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: C3H8 + H <=> N-C3H7 + H2 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: C3H8 + O <=> I-C3H7 + OH */
    /*eqcon[163] *= 1;  */

    /*reaction 165: C3H8 + O <=> N-C3H7 + OH */
    /*eqcon[164] *= 1;  */

    /*reaction 166: C3H8 + OH <=> N-C3H7 + H2O */
    /*eqcon[165] *= 1;  */

    /*reaction 167: C3H8 + OH <=> I-C3H7 + H2O */
    /*eqcon[166] *= 1;  */

    /*reaction 168: C3H8 + HO2 <=> I-C3H7 + H2O2 */
    /*eqcon[167] *= 1;  */

    /*reaction 169: C3H8 + HO2 <=> N-C3H7 + H2O2 */
    /*eqcon[168] *= 1;  */

    /*reaction 170: I-C3H7 + C3H8 <=> N-C3H7 + C3H8 */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C3H6 + H (+M) <=> I-C3H7 (+M) */
    eqcon[170] *= 1e+06; 

    /*reaction 172: I-C3H7 + O2 <=> C3H6 + HO2 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: N-C3H7 (+M) <=> CH3 + C2H4 (+M) */
    eqcon[172] *= 1e-06; 

    /*reaction 174: H + C3H6 (+M) <=> N-C3H7 (+M) */
    eqcon[173] *= 1e+06; 

    /*reaction 175: N-C3H7 + O2 <=> C3H6 + HO2 */
    /*eqcon[174] *= 1;  */
}


/*compute the production rate for each species */
void productionRate(double * wdot, double * sc, double T)
{
    double qdot;

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[39];                /*Gibbs free energy */
    double Kc;                      /*equilibrium constant */
    double k_f;                     /*forward reaction rate */
    double k_r;                     /*reverse reaction rate */
    double q_f;                     /*forward progress rate */
    double q_r;                     /*reverse progress rate */
    double phi_f;                   /*forward phase space factor */
    double phi_r;                   /*reverse phase space factor */
    double alpha;                   /*enhancement */
    double redP;                    /*reduced pressure */
    double logPred;                 /*log of above */
    double F;                       /*fallof rate enhancement */

    double F_troe;                  /*TROE intermediate */
    double logFcent;                /*TROE intermediate */
    double troe;                    /*TROE intermediate */
    double troe_c;                  /*TROE intermediate */
    double troe_n;                  /*TROE intermediate */

    double X;                       /*SRI intermediate */
    double F_sri;                   /*SRI intermediate */
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.314 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 39; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*zero out wdot */
    for (id = 0; id < 39; ++id) {
        wdot[id] = 0.0;
    }

    /*reaction 1: H + O2 <=> OH + O */
    phi_f = sc[2]*sc[3];
    k_f = 1e-06 * 3.52e+16*exp(-0.7*tc[0]-8590.73/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = exp((g_RT[2] + g_RT[3]) - (g_RT[4] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 2: H2 + O <=> OH + H */
    phi_f = sc[6]*sc[5];
    k_f = 1e-06 * 50600*exp(2.67*tc[0]-3165.89/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[2];
    Kc = exp((g_RT[6] + g_RT[5]) - (g_RT[4] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 3: H2 + OH <=> H2O + H */
    phi_f = sc[6]*sc[4];
    k_f = 1e-06 * 1.17e+09*exp(1.3*tc[0]-1829.53/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[2];
    Kc = exp((g_RT[6] + g_RT[4]) - (g_RT[7] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 4: H2O + O <=> 2 OH */
    phi_f = sc[7]*sc[5];
    k_f = 1e-06 * 7.6*exp(3.84*tc[0]-6431.63/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = exp((g_RT[7] + g_RT[5]) - (2 * g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[4] += 2 * qdot;

    /*reaction 5: H + O + M <=> OH + M */
    phi_f = sc[2]*sc[5];
    alpha = mixture + 1.5*sc[6] + 11*sc[7] + 0.9*sc[10] + 2.8*sc[11];
    k_f = 1e-12 * alpha * 6.2e+16*exp(-0.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[5]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 6: H2 + O2 <=> 2 OH */
    phi_f = sc[6]*sc[3];
    k_f = 1e-06 * 1.7e+13*exp(-24063/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = exp((g_RT[6] + g_RT[3]) - (2 * g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 2 * qdot;

    /*reaction 7: 2 H + M <=> H2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.5*sc[6] + 15.3*sc[7] + 0.9*sc[10] + 2.8*sc[11];
    k_f = 1e-12 * alpha * 7.2e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = 1.0 / (refC) * exp((2 * g_RT[2]) - (g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 2 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 8: H + OH + M <=> H2O + M */
    phi_f = sc[2]*sc[4];
    alpha = mixture + 1.5*sc[6] + 11*sc[7] + 0.9*sc[10] + 2.8*sc[11];
    k_f = 1e-12 * alpha * 3.8e+22*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[7];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[4]) - (g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 9: 2 O + M <=> O2 + M */
    phi_f = sc[5]*sc[5];
    alpha = mixture + 1.5*sc[6] + 11*sc[7] + 0.9*sc[10] + 2.8*sc[11];
    k_f = 1e-12 * alpha * 6.17e+15*exp(-0.5*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = 1.0 / (refC) * exp((2 * g_RT[5]) - (g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 2 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 10: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[2]*sc[3];
    alpha = mixture + -0.5*sc[1] + -0.7*sc[3] + 6*sc[7] + -0.25*sc[10] + 0.5*sc[11] + 0.5*sc[21];
    k_f = 1e-06 * 4.65e+12*exp(0.44*tc[0]);
    redP = 1e-12 * alpha / k_f * 2.6e+19*exp(-1.2*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((1*exp(T/-345))+ (0*exp(T/-1))+ (exp(-345/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[8];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[3]) - (g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 11: O + OH + M <=> HO2 + M */
    phi_f = sc[5]*sc[4];
    alpha = mixture;
    k_f = 1e-12 * alpha * 1e+16;
    q_f = phi_f * k_f;
    phi_r = sc[8];
    Kc = 1.0 / (refC) * exp((g_RT[5] + g_RT[4]) - (g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 12: HO2 + H <=> 2 OH */
    phi_f = sc[8]*sc[2];
    k_f = 1e-06 * 7.08e+13*exp(-150.956/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = exp((g_RT[8] + g_RT[2]) - (2 * g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[4] += 2 * qdot;

    /*reaction 13: HO2 + H <=> H2 + O2 */
    phi_f = sc[8]*sc[2];
    k_f = 1e-06 * 4.28e+13*exp(-709.678/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = exp((g_RT[8] + g_RT[2]) - (g_RT[6] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 14: HO2 + H <=> H2O + O */
    phi_f = sc[8]*sc[2];
    k_f = 1e-06 * 3.1e+13*exp(-866.049/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = exp((g_RT[8] + g_RT[2]) - (g_RT[7] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 15: HO2 + O <=> OH + O2 */
    phi_f = sc[8]*sc[5];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = exp((g_RT[8] + g_RT[5]) - (g_RT[4] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    phi_f = sc[8]*sc[4];
    k_f = 1e-06 * 2.89e+13*exp(+250.191/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[3];
    Kc = exp((g_RT[8] + g_RT[4]) - (g_RT[7] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 17: 2 OH (+M) <=> H2O2 (+M) */
    phi_f = sc[4]*sc[4];
    alpha = mixture + -0.3*sc[1] + sc[6] + 5*sc[7] + 0.5*sc[10] + sc[11] + sc[14] + 2*sc[21];
    k_f = 1e-06 * 7.4e+13*exp(-0.37*tc[0]);
    redP = 1e-12 * alpha / k_f * 1.34e+17*exp(-0.584*tc[0]+1154.74/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.265*exp(T/-94))+ (0.735*exp(T/-1756))+ (exp(-5182/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[9];
    Kc = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 2 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 18: 2 HO2 <=> H2O2 + O2 */
    phi_f = sc[8]*sc[8];
    k_f = 1e-06 * 3.02e+12*exp(-697.65/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[3];
    Kc = exp((2 * g_RT[8]) - (g_RT[9] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 2 * qdot;
    wdot[9] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 19: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[9]*sc[2];
    k_f = 1e-06 * 4.79e+13*exp(-4005.48/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[6];
    Kc = exp((g_RT[9] + g_RT[2]) - (g_RT[8] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 20: H2O2 + H <=> H2O + OH */
    phi_f = sc[9]*sc[2];
    k_f = 1e-06 * 1e+13*exp(-1804.27/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = exp((g_RT[9] + g_RT[2]) - (g_RT[7] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 21: H2O2 + OH <=> H2O + HO2 */
    phi_f = sc[9]*sc[4];
    k_f = 1e-06 * 7.08e+12*exp(-721.706/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[8];
    Kc = exp((g_RT[9] + g_RT[4]) - (g_RT[7] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 22: H2O2 + O <=> HO2 + OH */
    phi_f = sc[9]*sc[5];
    k_f = 1e-06 * 9.63e+06*exp(2*tc[0]-2008.76/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[4];
    Kc = exp((g_RT[9] + g_RT[5]) - (g_RT[8] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 23: CO + OH <=> CO2 + H */
    phi_f = sc[10]*sc[4];
    k_f = 1e-06 * 4.4e+06*exp(1.5*tc[0]+372.884/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[2];
    Kc = exp((g_RT[10] + g_RT[4]) - (g_RT[11] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 24: CO + HO2 <=> CO2 + OH */
    phi_f = sc[10]*sc[8];
    k_f = 1e-06 * 6.03e+13*exp(-11547.3/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[4];
    Kc = exp((g_RT[10] + g_RT[8]) - (g_RT[11] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 25: HCO + M <=> CO + H + M */
    phi_f = sc[12];
    alpha = mixture + 0.9*sc[6] + 11*sc[7] + 1.5*sc[10] + 1.5*sc[11];
    k_f = 1e-06 * alpha * 1.86e+17*exp(-1*tc[0]-8555.85/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[2];
    Kc = refC * exp((g_RT[12]) - (g_RT[10] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 26: HCO + H <=> CO + H2 */
    phi_f = sc[12]*sc[2];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[6];
    Kc = exp((g_RT[12] + g_RT[2]) - (g_RT[10] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 27: HCO + O <=> CO + OH */
    phi_f = sc[12]*sc[5];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[4];
    Kc = exp((g_RT[12] + g_RT[5]) - (g_RT[10] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 28: HCO + O <=> CO2 + H */
    phi_f = sc[12]*sc[5];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[2];
    Kc = exp((g_RT[12] + g_RT[5]) - (g_RT[11] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 29: HCO + OH <=> CO + H2O */
    phi_f = sc[12]*sc[4];
    k_f = 1e-06 * 5.02e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[7];
    Kc = exp((g_RT[12] + g_RT[4]) - (g_RT[10] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 30: HCO + O2 <=> CO + HO2 */
    phi_f = sc[12]*sc[3];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[8];
    Kc = exp((g_RT[12] + g_RT[3]) - (g_RT[10] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 31: CH2O + M <=> HCO + H + M */
    phi_f = sc[13];
    alpha = mixture + 1.5*sc[6] + 15.3*sc[7] + 0.9*sc[10] + 2.8*sc[11];
    k_f = 1e-06 * alpha * 6.26e+16*exp(-39212.8/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[2];
    Kc = refC * exp((g_RT[13]) - (g_RT[12] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 32: CH2O + H <=> HCO + H2 */
    phi_f = sc[13]*sc[2];
    k_f = 1e-06 * 1.26e+08*exp(1.62*tc[0]-1089.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[6];
    Kc = exp((g_RT[13] + g_RT[2]) - (g_RT[12] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 33: CH2O + O <=> HCO + OH */
    phi_f = sc[13]*sc[5];
    k_f = 1e-06 * 3.5e+13*exp(-1768.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[4];
    Kc = exp((g_RT[13] + g_RT[5]) - (g_RT[12] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 34: CH2O + OH <=> HCO + H2O */
    phi_f = sc[13]*sc[4];
    k_f = 1e-06 * 3.9e+10*exp(0.89*tc[0]-204.484/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[7];
    Kc = exp((g_RT[13] + g_RT[4]) - (g_RT[12] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 35: CH4 + H <=> H2 + CH3 */
    phi_f = sc[14]*sc[2];
    k_f = 1e-06 * 13000*exp(3*tc[0]-4045.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[15];
    Kc = exp((g_RT[14] + g_RT[2]) - (g_RT[6] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[14] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 36: CH4 + OH <=> H2O + CH3 */
    phi_f = sc[14]*sc[4];
    k_f = 1e-06 * 1.6e+07*exp(1.83*tc[0]-1400.12/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[15];
    Kc = exp((g_RT[14] + g_RT[4]) - (g_RT[7] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[14] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 37: CH4 + O <=> CH3 + OH */
    phi_f = sc[14]*sc[5];
    k_f = 1e-06 * 1.9e+09*exp(1.44*tc[0]-4366.34/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[4];
    Kc = exp((g_RT[14] + g_RT[5]) - (g_RT[15] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[14] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 38: CH4 + O2 <=> CH3 + HO2 */
    phi_f = sc[14]*sc[3];
    k_f = 1e-06 * 3.98e+13*exp(-28631.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[8];
    Kc = exp((g_RT[14] + g_RT[3]) - (g_RT[15] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[14] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 39: CH4 + HO2 <=> CH3 + H2O2 */
    phi_f = sc[14]*sc[8];
    k_f = 1e-06 * 9.03e+12*exp(-12401.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[9];
    Kc = exp((g_RT[14] + g_RT[8]) - (g_RT[15] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[14] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 40: CH3 + H <=> T-CH2 + H2 */
    phi_f = sc[15]*sc[2];
    k_f = 1e-06 * 1.8e+14*exp(-7601.99/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[6];
    Kc = exp((g_RT[15] + g_RT[2]) - (g_RT[16] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 41: CH3 + H <=> S-CH2 + H2 */
    phi_f = sc[15]*sc[2];
    k_f = 1e-06 * 1.55e+14*exp(-6784.06/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[6];
    Kc = exp((g_RT[15] + g_RT[2]) - (g_RT[17] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 42: CH3 + OH <=> S-CH2 + H2O */
    phi_f = sc[15]*sc[4];
    k_f = 1e-06 * 1e+13*exp(-1259.38/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[7];
    Kc = exp((g_RT[15] + g_RT[4]) - (g_RT[17] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 43: CH3 + O <=> CH2O + H */
    phi_f = sc[15]*sc[5];
    k_f = 1e-06 * 8.43e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = exp((g_RT[15] + g_RT[5]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 44: CH3 + T-CH2 <=> C2H4 + H */
    phi_f = sc[15]*sc[16];
    k_f = 1e-06 * 4.22e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[2];
    Kc = exp((g_RT[15] + g_RT[16]) - (g_RT[18] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 45: CH3 + HO2 <=> CH3O + OH */
    phi_f = sc[15]*sc[8];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[4];
    Kc = exp((g_RT[15] + g_RT[8]) - (g_RT[19] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 46: CH3 + O2 <=> CH2O + OH */
    phi_f = sc[15]*sc[3];
    k_f = 1e-06 * 3.3e+11*exp(-4499.85/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[4];
    Kc = exp((g_RT[15] + g_RT[3]) - (g_RT[13] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 47: CH3 + O2 <=> CH3O + O */
    phi_f = sc[15]*sc[3];
    k_f = 1e-06 * 1.33e+14*exp(-15805.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[5];
    Kc = exp((g_RT[15] + g_RT[3]) - (g_RT[19] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 48: 2 CH3 <=> C2H4 + H2 */
    phi_f = sc[15]*sc[15];
    k_f = 1e-06 * 1e+14*exp(-16106.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[6];
    Kc = exp((2 * g_RT[15]) - (g_RT[18] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 2 * qdot;
    wdot[18] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 49: 2 CH3 <=> C2H5 + H */
    phi_f = sc[15]*sc[15];
    k_f = 1e-06 * 3.16e+13*exp(-7397.51/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[2];
    Kc = exp((2 * g_RT[15]) - (g_RT[20] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 2 * qdot;
    wdot[20] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 50: CH3 + H (+M) <=> CH4 (+M) */
    phi_f = sc[15]*sc[2];
    alpha = mixture;
    k_f = 1e-06 * 2.11e+14;
    redP = 1e-12 * alpha / k_f * 6.26e+23*exp(-1.8*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.63*exp(T/-3315))+ (0.37*exp(T/-61)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[15] + g_RT[2]) - (g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 51: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[15]*sc[15];
    alpha = mixture;
    k_f = 1e-06 * 1.81e+13;
    redP = 1e-12 * alpha / k_f * 1.27e+41*exp(-7*tc[0]-1390.49/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.38*exp(T/-73))+ (0.62*exp(T/-1180)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[21];
    Kc = 1.0 / (refC) * exp((2 * g_RT[15]) - (g_RT[21]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 2 * qdot;
    wdot[21] += 1 * qdot;

    /*reaction 52: S-CH2 + OH <=> CH2O + H */
    phi_f = sc[17]*sc[4];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = exp((g_RT[17] + g_RT[4]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 53: S-CH2 + O2 <=> CO + OH + H */
    phi_f = sc[17]*sc[3];
    k_f = 1e-06 * 3.13e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[4]*sc[2];
    Kc = refC * exp((g_RT[17] + g_RT[3]) - (g_RT[10] + g_RT[4] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 54: S-CH2 + CO2 <=> CO + CH2O */
    phi_f = sc[17]*sc[11];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[13];
    Kc = exp((g_RT[17] + g_RT[11]) - (g_RT[10] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 55: S-CH2 + M <=> T-CH2 + M */
    phi_f = sc[17];
    alpha = mixture + 1.4*sc[6] + 14.4*sc[7] + 0.8*sc[10] + 2.6*sc[11];
    k_f = 1e-06 * alpha * 6e+12;
    q_f = phi_f * k_f;
    phi_r = sc[16];
    Kc = exp((g_RT[17]) - (g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 56: T-CH2 + H <=> CH + H2 */
    phi_f = sc[16]*sc[2];
    k_f = 1e-06 * 6.02e+12*exp(+899.728/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[22]*sc[6];
    Kc = exp((g_RT[16] + g_RT[2]) - (g_RT[22] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[22] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 57: T-CH2 + OH <=> CH2O + H */
    phi_f = sc[16]*sc[4];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = exp((g_RT[16] + g_RT[4]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 58: T-CH2 + OH <=> CH + H2O */
    phi_f = sc[16]*sc[4];
    k_f = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.57/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[22]*sc[7];
    Kc = exp((g_RT[16] + g_RT[4]) - (g_RT[22] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[22] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 59: T-CH2 + O <=> CO + 2 H */
    phi_f = sc[16]*sc[5];
    k_f = 1e-06 * 8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[2]*sc[2];
    Kc = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[10] + 2 * g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[2] += 2 * qdot;

    /*reaction 60: T-CH2 + O <=> CO + H2 */
    phi_f = sc[16]*sc[5];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[6];
    Kc = exp((g_RT[16] + g_RT[5]) - (g_RT[10] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 61: T-CH2 + O2 <=> CO2 + H2 */
    phi_f = sc[16]*sc[3];
    k_f = 1e-06 * 2.63e+12*exp(-750.579/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[6];
    Kc = exp((g_RT[16] + g_RT[3]) - (g_RT[11] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 62: T-CH2 + O2 <=> CO + OH + H */
    phi_f = sc[16]*sc[3];
    k_f = 1e-06 * 6.58e+12*exp(-750.579/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[4]*sc[2];
    Kc = refC * exp((g_RT[16] + g_RT[3]) - (g_RT[10] + g_RT[4] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 63: 2 T-CH2 <=> C2H2 + 2 H */
    phi_f = sc[16]*sc[16];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[2]*sc[2];
    Kc = refC * exp((2 * g_RT[16]) - (g_RT[23] + 2 * g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 2 * qdot;
    wdot[23] += 1 * qdot;
    wdot[2] += 2 * qdot;

    /*reaction 64: CH + O <=> CO + H */
    phi_f = sc[22]*sc[5];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[2];
    Kc = exp((g_RT[22] + g_RT[5]) - (g_RT[10] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[22] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 65: CH + O2 <=> HCO + O */
    phi_f = sc[22]*sc[3];
    k_f = 1e-06 * 1.77e+11*exp(0.76*tc[0]+240.569/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[5];
    Kc = exp((g_RT[22] + g_RT[3]) - (g_RT[12] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[22] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 66: CH + H2O <=> CH2O + H */
    phi_f = sc[22]*sc[7];
    k_f = 1e-06 * 1.17e+15*exp(-0.75*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = exp((g_RT[22] + g_RT[7]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[22] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 67: CH + CO2 <=> HCO + CO */
    phi_f = sc[22]*sc[11];
    k_f = 1e-06 * 48*exp(3.22*tc[0]+1623.84/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[10];
    Kc = exp((g_RT[22] + g_RT[11]) - (g_RT[12] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[22] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 68: CH3O + H <=> CH2O + H2 */
    phi_f = sc[19]*sc[2];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[6];
    Kc = exp((g_RT[19] + g_RT[2]) - (g_RT[13] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 69: CH3O + H <=> S-CH2 + H2O */
    phi_f = sc[19]*sc[2];
    k_f = 1e-06 * 1.6e+13;
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[7];
    Kc = exp((g_RT[19] + g_RT[2]) - (g_RT[17] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 70: CH3O + OH <=> CH2O + H2O */
    phi_f = sc[19]*sc[4];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[7];
    Kc = exp((g_RT[19] + g_RT[4]) - (g_RT[13] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 71: CH3O + O <=> OH + CH2O */
    phi_f = sc[19]*sc[5];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = exp((g_RT[19] + g_RT[5]) - (g_RT[4] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 72: CH3O + O2 <=> CH2O + HO2 */
    phi_f = sc[19]*sc[3];
    k_f = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1780.21/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[8];
    Kc = exp((g_RT[19] + g_RT[3]) - (g_RT[13] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 73: CH3O + M <=> CH2O + H + M */
    phi_f = sc[19];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1e+13*exp(-6796.08/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = refC * exp((g_RT[19]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 74: C2H6 + H <=> C2H5 + H2 */
    phi_f = sc[21]*sc[2];
    k_f = 1e-06 * 540*exp(3.5*tc[0]-2622.21/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[6];
    Kc = exp((g_RT[21] + g_RT[2]) - (g_RT[20] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 75: C2H6 + O <=> C2H5 + OH */
    phi_f = sc[21]*sc[5];
    k_f = 1e-06 * 1.4*exp(4.3*tc[0]-1395.3/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[4];
    Kc = exp((g_RT[21] + g_RT[5]) - (g_RT[20] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 76: C2H6 + OH <=> C2H5 + H2O */
    phi_f = sc[21]*sc[4];
    k_f = 1e-06 * 2.2e+07*exp(1.9*tc[0]-565.34/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[7];
    Kc = exp((g_RT[21] + g_RT[4]) - (g_RT[20] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 77: C2H6 + CH3 <=> C2H5 + CH4 */
    phi_f = sc[21]*sc[15];
    k_f = 1e-06 * 0.55*exp(4*tc[0]-4173.88/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[14];
    Kc = exp((g_RT[21] + g_RT[15]) - (g_RT[20] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 78: C2H6 (+M) <=> C2H5 + H (+M) */
    phi_f = sc[21];
    alpha = mixture;
    k_f = 1 * 8.85e+20*exp(-1.23*tc[0]-51445.8/tc[1]);
    redP = 1e-12 * alpha / k_f * 4.9e+42*exp(-6.43*tc[0]-53935.7/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.16*exp(T/-125))+ (0.84*exp(T/-2219))+ (exp(-6882/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[2];
    Kc = refC * exp((g_RT[21]) - (g_RT[20] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 79: C2H5 + H <=> C2H4 + H2 */
    phi_f = sc[20]*sc[2];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[6];
    Kc = exp((g_RT[20] + g_RT[2]) - (g_RT[18] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 80: C2H5 + O <=> C2H4 + OH */
    phi_f = sc[20]*sc[5];
    k_f = 1e-06 * 3.06e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[4];
    Kc = exp((g_RT[20] + g_RT[5]) - (g_RT[18] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 81: C2H5 + O <=> CH3 + CH2O */
    phi_f = sc[20]*sc[5];
    k_f = 1e-06 * 4.24e+13;
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[13];
    Kc = exp((g_RT[20] + g_RT[5]) - (g_RT[15] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 82: C2H5 + O2 <=> C2H4 + HO2 */
    phi_f = sc[20]*sc[3];
    k_f = 1e-06 * 2e+12*exp(-2513.95/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[8];
    Kc = exp((g_RT[20] + g_RT[3]) - (g_RT[18] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 83: C2H5 (+M) <=> C2H4 + H (+M) */
    phi_f = sc[20];
    alpha = mixture;
    k_f = 1 * 1.11e+10*exp(1.037*tc[0]-18504.6/tc[1]);
    redP = 1e-12 * alpha / k_f * 3.99e+33*exp(-4.99*tc[0]-20130.9/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.832*exp(T/-1203)) /*+ (0.168*exp(T/-0))*/);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[2];
    Kc = refC * exp((g_RT[20]) - (g_RT[18] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 84: C2H4 + H <=> C2H3 + H2 */
    phi_f = sc[18]*sc[2];
    k_f = 1e-06 * 4.49e+07*exp(2.12*tc[0]-6723.92/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[6];
    Kc = exp((g_RT[18] + g_RT[2]) - (g_RT[24] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[24] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 85: C2H4 + OH <=> C2H3 + H2O */
    phi_f = sc[18]*sc[4];
    k_f = 1e-06 * 553000*exp(2.31*tc[0]-1491.53/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[7];
    Kc = exp((g_RT[18] + g_RT[4]) - (g_RT[24] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[24] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 86: C2H4 + O <=> CH3 + HCO */
    phi_f = sc[18]*sc[5];
    k_f = 1e-06 * 2.25e+06*exp(2.08*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[12];
    Kc = exp((g_RT[18] + g_RT[5]) - (g_RT[15] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 87: C2H4 + O <=> CH2CHO + H */
    phi_f = sc[18]*sc[5];
    k_f = 1e-06 * 1.21e+06*exp(2.08*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[2];
    Kc = exp((g_RT[18] + g_RT[5]) - (g_RT[25] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 88: 2 C2H4 <=> C2H3 + C2H5 */
    phi_f = sc[18]*sc[18];
    k_f = 1e-06 * 5.01e+14*exp(-32561.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[20];
    Kc = exp((2 * g_RT[18]) - (g_RT[24] + g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 2 * qdot;
    wdot[24] += 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 89: C2H4 + O2 <=> C2H3 + HO2 */
    phi_f = sc[18]*sc[3];
    k_f = 1e-06 * 4.22e+13*exp(-29000/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[8];
    Kc = exp((g_RT[18] + g_RT[3]) - (g_RT[24] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[24] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 90: C2H4 + HO2 <=> C2H4O + OH */
    phi_f = sc[18]*sc[8];
    k_f = 1e-06 * 2.23e+12*exp(-8650.88/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[26]*sc[4];
    Kc = exp((g_RT[18] + g_RT[8]) - (g_RT[26] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[26] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 91: C2H4O + HO2 <=> CH3 + CO + H2O2 */
    phi_f = sc[26]*sc[8];
    k_f = 1e-06 * 4e+12*exp(-8559.46/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[10]*sc[9];
    Kc = refC * exp((g_RT[26] + g_RT[8]) - (g_RT[15] + g_RT[10] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[26] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 92: C2H4 + M <=> C2H3 + H + M */
    phi_f = sc[18];
    alpha = mixture;
    k_f = 1e-06 * alpha * 2.6e+17*exp(-48600/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[2];
    Kc = refC * exp((g_RT[18]) - (g_RT[24] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[24] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 93: C2H4 + M <=> C2H2 + H2 + M */
    phi_f = sc[18];
    alpha = mixture;
    k_f = 1e-06 * alpha * 3.5e+16*exp(-36000/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[6];
    Kc = refC * exp((g_RT[18]) - (g_RT[23] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 94: C2H3 + H <=> C2H2 + H2 */
    phi_f = sc[24]*sc[2];
    k_f = 1e-06 * 1.21e+13;
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[6];
    Kc = exp((g_RT[24] + g_RT[2]) - (g_RT[23] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 95: C2H3 (+M) <=> C2H2 + H (+M) */
    phi_f = sc[24];
    alpha = mixture;
    k_f = 1 * 6.38e+09*exp(1*tc[0]-18936.4/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.51e+14*exp(0.1*tc[0]-16450.1/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.7*exp(T/-1e+30))+ (0.3*exp(T/-1e-30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[2];
    Kc = refC * exp((g_RT[24]) - (g_RT[23] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 96: C2H3 + O2 <=> CH2O + HCO */
    phi_f = sc[24]*sc[3];
    k_f = 1e-06 * 1.7e+29*exp(-5.312*tc[0]-3272.83/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[12];
    Kc = exp((g_RT[24] + g_RT[3]) - (g_RT[13] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 97: C2H3 + O2 <=> CH2CHO + O */
    phi_f = sc[24]*sc[3];
    k_f = 1e-06 * 7e+14*exp(-0.611*tc[0]-2648.43/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[5];
    Kc = exp((g_RT[24] + g_RT[3]) - (g_RT[25] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 98: C2H3 + O2 <=> C2H2 + HO2 */
    phi_f = sc[24]*sc[3];
    k_f = 1e-06 * 5.19e+15*exp(-1.26*tc[0]-1667.15/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[8];
    Kc = exp((g_RT[24] + g_RT[3]) - (g_RT[23] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 99: CH2CHO <=> CH2CO + H */
    phi_f = sc[25];
    k_f = 1 * 1.047e+37*exp(-7.189*tc[0]-22315.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[2];
    Kc = refC * exp((g_RT[25]) - (g_RT[27] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[25] -= 1 * qdot;
    wdot[27] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 100: C2H2 + O <=> HCCO + H */
    phi_f = sc[23]*sc[5];
    k_f = 1e-06 * 4e+14*exp(-5364.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[2];
    Kc = exp((g_RT[23] + g_RT[5]) - (g_RT[28] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[28] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 101: C2H2 + O <=> T-CH2 + CO */
    phi_f = sc[23]*sc[5];
    k_f = 1e-06 * 1.6e+14*exp(-4979.79/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[10];
    Kc = exp((g_RT[23] + g_RT[5]) - (g_RT[16] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 102: C2H2 + O2 <=> CH2O + CO */
    phi_f = sc[23]*sc[3];
    k_f = 1e-06 * 4.6e+15*exp(-0.54*tc[0]-22613.5/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[10];
    Kc = exp((g_RT[23] + g_RT[3]) - (g_RT[13] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 103: C2H2 + OH <=> CH2CO + H */
    phi_f = sc[23]*sc[4];
    k_f = 1e-06 * 1.9e+07*exp(1.7*tc[0]-502.788/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[2];
    Kc = exp((g_RT[23] + g_RT[4]) - (g_RT[27] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[27] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 104: C2H2 + OH <=> C2H + H2O */
    phi_f = sc[23]*sc[4];
    k_f = 1e-06 * 3.37e+07*exp(2*tc[0]-7046.28/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[7];
    Kc = exp((g_RT[23] + g_RT[4]) - (g_RT[29] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[29] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 105: CH2CO + H <=> CH3 + CO */
    phi_f = sc[27]*sc[2];
    k_f = 1e-06 * 1.5e+09*exp(1.43*tc[0]-1353.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[10];
    Kc = exp((g_RT[27] + g_RT[2]) - (g_RT[15] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[27] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 106: CH2CO + O <=> T-CH2 + CO2 */
    phi_f = sc[27]*sc[5];
    k_f = 1e-06 * 2e+13*exp(-1154.74/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[11];
    Kc = exp((g_RT[27] + g_RT[5]) - (g_RT[16] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[27] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 107: CH2CO + O <=> HCCO + OH */
    phi_f = sc[27]*sc[5];
    k_f = 1e-06 * 1e+13*exp(-1006.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[4];
    Kc = exp((g_RT[27] + g_RT[5]) - (g_RT[28] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[27] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[28] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 108: CH2CO + CH3 <=> C2H5 + CO */
    phi_f = sc[27]*sc[15];
    k_f = 1e-06 * 9e+10;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[10];
    Kc = exp((g_RT[27] + g_RT[15]) - (g_RT[20] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[27] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 109: HCCO + H <=> S-CH2 + CO */
    phi_f = sc[28]*sc[2];
    k_f = 1e-06 * 1.5e+14;
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[10];
    Kc = exp((g_RT[28] + g_RT[2]) - (g_RT[17] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 110: HCCO + OH <=> HCO + CO + H */
    phi_f = sc[28]*sc[4];
    k_f = 1e-06 * 2e+12;
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[10]*sc[2];
    Kc = refC * exp((g_RT[28] + g_RT[4]) - (g_RT[12] + g_RT[10] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 111: HCCO + O <=> 2 CO + H */
    phi_f = sc[28]*sc[5];
    k_f = 1e-06 * 9.64e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[10]*sc[2];
    Kc = refC * exp((g_RT[28] + g_RT[5]) - (2 * g_RT[10] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[10] += 2 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 112: HCCO + O2 <=> 2 CO + OH */
    phi_f = sc[28]*sc[3];
    k_f = 1e-06 * 2.88e+07*exp(1.7*tc[0]-503.991/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[10]*sc[4];
    Kc = refC * exp((g_RT[28] + g_RT[3]) - (2 * g_RT[10] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[10] += 2 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 113: HCCO + O2 <=> CO2 + CO + H */
    phi_f = sc[28]*sc[3];
    k_f = 1e-06 * 1.4e+07*exp(1.7*tc[0]-503.991/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[10]*sc[2];
    Kc = refC * exp((g_RT[28] + g_RT[3]) - (g_RT[11] + g_RT[10] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 114: C2H + OH <=> HCCO + H */
    phi_f = sc[29]*sc[4];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[2];
    Kc = exp((g_RT[29] + g_RT[4]) - (g_RT[28] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[29] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[28] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 115: C2H + O <=> CO + CH */
    phi_f = sc[29]*sc[5];
    k_f = 1e-06 * 1.02e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[22];
    Kc = exp((g_RT[29] + g_RT[5]) - (g_RT[10] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[29] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 116: C2H + O2 <=> HCCO + O */
    phi_f = sc[29]*sc[3];
    k_f = 1e-06 * 6.02e+11;
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[5];
    Kc = exp((g_RT[29] + g_RT[3]) - (g_RT[28] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[29] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[28] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 117: C2H + O2 <=> CH + CO2 */
    phi_f = sc[29]*sc[3];
    k_f = 1e-06 * 4.5e+15*exp(-12629.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[22]*sc[11];
    Kc = exp((g_RT[29] + g_RT[3]) - (g_RT[22] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[29] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[22] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 118: C2H + O2 <=> HCO + CO */
    phi_f = sc[29]*sc[3];
    k_f = 1e-06 * 2.41e+12;
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[10];
    Kc = exp((g_RT[29] + g_RT[3]) - (g_RT[12] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[29] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 119: CH2OH + H <=> CH2O + H2 */
    phi_f = sc[30]*sc[2];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[6];
    Kc = exp((g_RT[30] + g_RT[2]) - (g_RT[13] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 120: CH2OH + H <=> CH3 + OH */
    phi_f = sc[30]*sc[2];
    k_f = 1e-06 * 1.75e+14*exp(-1407.33/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[4];
    Kc = exp((g_RT[30] + g_RT[2]) - (g_RT[15] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 121: CH2OH + OH <=> CH2O + H2O */
    phi_f = sc[30]*sc[4];
    k_f = 1e-06 * 2.4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[7];
    Kc = exp((g_RT[30] + g_RT[4]) - (g_RT[13] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 122: CH2OH + O2 <=> CH2O + HO2 */
    phi_f = sc[30]*sc[3];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[8];
    Kc = exp((g_RT[30] + g_RT[3]) - (g_RT[13] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 123: CH2OH + M <=> CH2O + H + M */
    phi_f = sc[30];
    alpha = mixture + 1.4*sc[6] + 14.4*sc[7] + 0.8*sc[10] + 2.6*sc[11];
    k_f = 1e-06 * alpha * 5e+13*exp(-12641.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = refC * exp((g_RT[30]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 124: CH3O + M <=> CH2OH + M */
    phi_f = sc[19];
    alpha = mixture + 1.5*sc[6] + 11*sc[7] + 0.9*sc[10] + 2.8*sc[11];
    k_f = 1e-06 * alpha * 1e+14*exp(-9622.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30];
    Kc = exp((g_RT[19]) - (g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[30] += 1 * qdot;

    /*reaction 125: CH2CO + OH <=> CH2OH + CO */
    phi_f = sc[27]*sc[4];
    k_f = 1e-06 * 1.02e+13;
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[10];
    Kc = exp((g_RT[27] + g_RT[4]) - (g_RT[30] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[27] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 126: CH3OH + OH <=> CH2OH + H2O */
    phi_f = sc[31]*sc[4];
    k_f = 1e-06 * 1.44e+06*exp(2*tc[0]+422.199/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[7];
    Kc = exp((g_RT[31] + g_RT[4]) - (g_RT[30] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 127: CH3OH + OH <=> CH3O + H2O */
    phi_f = sc[31]*sc[4];
    k_f = 1e-06 * 6.3e+06*exp(2*tc[0]-757.796/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[7];
    Kc = exp((g_RT[31] + g_RT[4]) - (g_RT[19] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 128: CH3OH + H <=> CH2OH + H2 */
    phi_f = sc[31]*sc[2];
    k_f = 1e-06 * 1.64e+07*exp(2*tc[0]-2273.38/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[6];
    Kc = exp((g_RT[31] + g_RT[2]) - (g_RT[30] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 129: CH3OH + H <=> CH3O + H2 */
    phi_f = sc[31]*sc[2];
    k_f = 1e-06 * 3.83e+07*exp(2*tc[0]-2946.98/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[6];
    Kc = exp((g_RT[31] + g_RT[2]) - (g_RT[19] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 130: CH3OH + O <=> CH2OH + OH */
    phi_f = sc[31]*sc[5];
    k_f = 1e-06 * 1e+13*exp(-2357.58/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[4];
    Kc = exp((g_RT[31] + g_RT[5]) - (g_RT[30] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 131: CH3OH + HO2 <=> CH2OH + H2O2 */
    phi_f = sc[31]*sc[8];
    k_f = 1e-06 * 6.2e+12*exp(-9755.09/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[9];
    Kc = exp((g_RT[31] + g_RT[8]) - (g_RT[30] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 132: CH3OH + O2 <=> CH2OH + HO2 */
    phi_f = sc[31]*sc[3];
    k_f = 1e-06 * 2e+13*exp(-22613.5/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[8];
    Kc = exp((g_RT[31] + g_RT[3]) - (g_RT[30] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 133: C3H4 + O <=> C2H4 + CO */
    phi_f = sc[32]*sc[5];
    k_f = 1e-06 * 2e+07*exp(1.8*tc[0]-503.271/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[10];
    Kc = exp((g_RT[32] + g_RT[5]) - (g_RT[18] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 134: CH3 + C2H2 <=> C3H4 + H */
    phi_f = sc[15]*sc[23];
    k_f = 1e-06 * 2.56e+09*exp(1.1*tc[0]-6866.57/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[2];
    Kc = exp((g_RT[15] + g_RT[23]) - (g_RT[32] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[23] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 135: C3H4 + O <=> HCCO + CH3 */
    phi_f = sc[32]*sc[5];
    k_f = 1e-06 * 7.3e+12*exp(-1132.36/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[15];
    Kc = exp((g_RT[32] + g_RT[5]) - (g_RT[28] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[28] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 136: C3H3 + H (+M) <=> C3H4 (+M) */
    phi_f = sc[33]*sc[2];
    alpha = mixture;
    k_f = 1e-06 * 3e+13;
    redP = 1e-12 * alpha / k_f * 9e+15*exp(1*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-1e+30))+ (0.5*exp(T/-1e-30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[32];
    Kc = 1.0 / (refC) * exp((g_RT[33] + g_RT[2]) - (g_RT[32]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[32] += 1 * qdot;

    /*reaction 137: C3H3 + HO2 <=> C3H4 + O2 */
    phi_f = sc[33]*sc[8];
    k_f = 1e-06 * 2.5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[3];
    Kc = exp((g_RT[33] + g_RT[8]) - (g_RT[32] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 138: C3H4 + OH <=> C3H3 + H2O */
    phi_f = sc[32]*sc[4];
    k_f = 1e-06 * 5.3e+06*exp(2*tc[0]-1006.54/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[33]*sc[7];
    Kc = exp((g_RT[32] + g_RT[4]) - (g_RT[33] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[33] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 139: C3H3 + O2 <=> CH2CO + HCO */
    phi_f = sc[33]*sc[3];
    k_f = 1e-06 * 3e+10*exp(-1443.42/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[12];
    Kc = exp((g_RT[33] + g_RT[3]) - (g_RT[27] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[27] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 140: C3H4 + H (+M) <=> C3H5 (+M) */
    phi_f = sc[32]*sc[2];
    alpha = mixture;
    k_f = 1e-06 * 4e+13;
    redP = 1e-12 * alpha / k_f * 3e+24*exp(-2*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2*exp(T/-1e+30))+ (0.8*exp(T/-1e-30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[34];
    Kc = 1.0 / (refC) * exp((g_RT[32] + g_RT[2]) - (g_RT[34]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[34] += 1 * qdot;

    /*reaction 141: C3H5 + H <=> C3H4 + H2 */
    phi_f = sc[34]*sc[2];
    k_f = 1e-06 * 1.8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[6];
    Kc = exp((g_RT[34] + g_RT[2]) - (g_RT[32] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 142: C3H5 + O2 <=> C3H4 + HO2 */
    phi_f = sc[34]*sc[3];
    k_f = 1e-06 * 4.99e+15*exp(-1.4*tc[0]-11287.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[8];
    Kc = exp((g_RT[34] + g_RT[3]) - (g_RT[32] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 143: C3H5 + CH3 <=> C3H4 + CH4 */
    phi_f = sc[34]*sc[15];
    k_f = 1e-06 * 3e+12*exp(-0.32*tc[0]+65.9185/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[14];
    Kc = exp((g_RT[34] + g_RT[15]) - (g_RT[32] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 144: C2H2 + CH3 (+M) <=> C3H5 (+M) */
    phi_f = sc[23]*sc[15];
    alpha = mixture;
    k_f = 1e-06 * 6e+08;
    redP = 1e-12 * alpha / k_f * 2e+09*exp(1*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-1e+30))+ (0.5*exp(T/-1e-30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[34];
    Kc = 1.0 / (refC) * exp((g_RT[23] + g_RT[15]) - (g_RT[34]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[34] += 1 * qdot;

    /*reaction 145: C3H5 + OH <=> C3H4 + H2O */
    phi_f = sc[34]*sc[4];
    k_f = 1e-06 * 6e+12;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[7];
    Kc = exp((g_RT[34] + g_RT[4]) - (g_RT[32] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 146: C3H3 + HCO <=> C3H4 + CO */
    phi_f = sc[33]*sc[12];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[10];
    Kc = exp((g_RT[33] + g_RT[12]) - (g_RT[32] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 147: C3H3 + HO2 <=> OH + CO + C2H3 */
    phi_f = sc[33]*sc[8];
    k_f = 1e-06 * 8e+11;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[10]*sc[24];
    Kc = refC * exp((g_RT[33] + g_RT[8]) - (g_RT[4] + g_RT[10] + g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 148: C3H4 + O2 <=> CH3 + HCO + CO */
    phi_f = sc[32]*sc[3];
    k_f = 1e-06 * 4e+14*exp(-21049.8/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[12]*sc[10];
    Kc = refC * exp((g_RT[32] + g_RT[3]) - (g_RT[15] + g_RT[12] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 149: C3H6 + O <=> C2H5 + HCO */
    phi_f = sc[35]*sc[5];
    k_f = 1e-06 * 3.5e+07*exp(1.65*tc[0]+489.557/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[12];
    Kc = exp((g_RT[35] + g_RT[5]) - (g_RT[20] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 150: C3H6 + OH <=> C3H5 + H2O */
    phi_f = sc[35]*sc[4];
    k_f = 1e-06 * 3.1e+06*exp(2*tc[0]+150.116/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[7];
    Kc = exp((g_RT[35] + g_RT[4]) - (g_RT[34] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[34] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 151: C3H6 + O <=> CH2CO + CH3 + H */
    phi_f = sc[35]*sc[5];
    k_f = 1e-06 * 1.2e+08*exp(1.65*tc[0]-164.791/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[15]*sc[2];
    Kc = refC * exp((g_RT[35] + g_RT[5]) - (g_RT[27] + g_RT[15] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[27] += 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 152: C3H6 + H <=> C3H5 + H2 */
    phi_f = sc[35]*sc[2];
    k_f = 1e-06 * 170000*exp(2.5*tc[0]-1254.57/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[6];
    Kc = exp((g_RT[35] + g_RT[2]) - (g_RT[34] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[34] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 153: C3H5 + H (+M) <=> C3H6 (+M) */
    phi_f = sc[34]*sc[2];
    alpha = mixture + -0.3*sc[1] + sc[6] + 5*sc[7] + 0.5*sc[10] + sc[11] + sc[14] + 2*sc[21];
    k_f = 1e-06 * 2e+14;
    redP = 1e-12 * alpha / k_f * 1.33e+60*exp(-12*tc[0]-3003.51/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.98*exp(T/-1097))+ (0.02*exp(T/-1097))+ (exp(-6860/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[35];
    Kc = 1.0 / (refC) * exp((g_RT[34] + g_RT[2]) - (g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[35] += 1 * qdot;

    /*reaction 154: C3H5 + HO2 <=> C3H6 + O2 */
    phi_f = sc[34]*sc[8];
    k_f = 1e-06 * 2.66e+12;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[3];
    Kc = exp((g_RT[34] + g_RT[8]) - (g_RT[35] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 155: C3H5 + HO2 <=> OH + C2H3 + CH2O */
    phi_f = sc[34]*sc[8];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[24]*sc[13];
    Kc = refC * exp((g_RT[34] + g_RT[8]) - (g_RT[4] + g_RT[24] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[24] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 156: C2H3 + CH3 (+M) <=> C3H6 (+M) */
    phi_f = sc[24]*sc[15];
    alpha = mixture + -0.3*sc[1] + sc[6] + 5*sc[7] + 0.5*sc[10] + sc[11] + sc[14] + 2*sc[21];
    k_f = 1e-06 * 2.5e+13;
    redP = 1e-12 * alpha / k_f * 4.27e+58*exp(-11.94*tc[0]-4917.24/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.825*exp(T/-1341))+ (0.175*exp(T/-60000))+ (exp(-10140/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[35];
    Kc = 1.0 / (refC) * exp((g_RT[24] + g_RT[15]) - (g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[35] += 1 * qdot;

    /*reaction 157: C3H6 + H <=> C2H4 + CH3 */
    phi_f = sc[35]*sc[2];
    k_f = 1e-06 * 1.6e+22*exp(-2.39*tc[0]-5629.33/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[15];
    Kc = exp((g_RT[35] + g_RT[2]) - (g_RT[18] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 158: CH3 + C2H3 <=> C3H5 + H */
    phi_f = sc[15]*sc[24];
    k_f = 1e-06 * 1.5e+24*exp(-2.83*tc[0]-9370.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[2];
    Kc = exp((g_RT[15] + g_RT[24]) - (g_RT[34] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[34] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 159: C3H8 (+M) <=> CH3 + C2H5 (+M) */
    phi_f = sc[36];
    alpha = mixture;
    k_f = 1 * 1.1e+17*exp(-42472.5/tc[1]);
    redP = 1e-12 * alpha / k_f * 7.83e+18*exp(-32701.8/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.24*exp(T/-1946))+ (0.76*exp(T/-38)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[20];
    Kc = refC * exp((g_RT[36]) - (g_RT[15] + g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 160: C3H8 + O2 <=> I-C3H7 + HO2 */
    phi_f = sc[36]*sc[3];
    k_f = 1e-06 * 4e+13*exp(-23905.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[8];
    Kc = exp((g_RT[36] + g_RT[3]) - (g_RT[37] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[37] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 161: C3H8 + O2 <=> N-C3H7 + HO2 */
    phi_f = sc[36]*sc[3];
    k_f = 1e-06 * 4e+13*exp(-25632.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[8];
    Kc = exp((g_RT[36] + g_RT[3]) - (g_RT[38] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[38] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 162: C3H8 + H <=> I-C3H7 + H2 */
    phi_f = sc[36]*sc[2];
    k_f = 1e-06 * 1.3e+06*exp(2.4*tc[0]-2250.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[6];
    Kc = exp((g_RT[36] + g_RT[2]) - (g_RT[37] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[37] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 163: C3H8 + H <=> N-C3H7 + H2 */
    phi_f = sc[36]*sc[2];
    k_f = 1e-06 * 1.33e+06*exp(2.54*tc[0]-3402.85/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[6];
    Kc = exp((g_RT[36] + g_RT[2]) - (g_RT[38] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[38] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 164: C3H8 + O <=> I-C3H7 + OH */
    phi_f = sc[36]*sc[5];
    k_f = 1e-06 * 47600*exp(2.71*tc[0]-1060.55/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[4];
    Kc = exp((g_RT[36] + g_RT[5]) - (g_RT[37] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[37] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 165: C3H8 + O <=> N-C3H7 + OH */
    phi_f = sc[36]*sc[5];
    k_f = 1e-06 * 190000*exp(2.68*tc[0]-1871.39/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[4];
    Kc = exp((g_RT[36] + g_RT[5]) - (g_RT[38] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[38] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 166: C3H8 + OH <=> N-C3H7 + H2O */
    phi_f = sc[36]*sc[4];
    k_f = 1e-06 * 1400*exp(2.66*tc[0]-265.35/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[7];
    Kc = exp((g_RT[36] + g_RT[4]) - (g_RT[38] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[38] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 167: C3H8 + OH <=> I-C3H7 + H2O */
    phi_f = sc[36]*sc[4];
    k_f = 1e-06 * 27000*exp(2.39*tc[0]-197.866/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[7];
    Kc = exp((g_RT[36] + g_RT[4]) - (g_RT[37] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[37] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 168: C3H8 + HO2 <=> I-C3H7 + H2O2 */
    phi_f = sc[36]*sc[8];
    k_f = 1e-06 * 9640*exp(2.6*tc[0]-7004.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[9];
    Kc = exp((g_RT[36] + g_RT[8]) - (g_RT[37] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[37] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 169: C3H8 + HO2 <=> N-C3H7 + H2O2 */
    phi_f = sc[36]*sc[8];
    k_f = 1e-06 * 47600*exp(2.55*tc[0]-8299.65/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[9];
    Kc = exp((g_RT[36] + g_RT[8]) - (g_RT[38] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[38] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 170: I-C3H7 + C3H8 <=> N-C3H7 + C3H8 */
    phi_f = sc[37]*sc[36];
    k_f = 1e-06 * 0.0084*exp(4.2*tc[0]-4366.34/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[36];
    Kc = exp((g_RT[37] + g_RT[36]) - (g_RT[38] + g_RT[36]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[37] -= 1 * qdot;
    wdot[36] -= 1 * qdot;
    wdot[38] += 1 * qdot;
    wdot[36] += 1 * qdot;

    /*reaction 171: C3H6 + H (+M) <=> I-C3H7 (+M) */
    phi_f = sc[35]*sc[2];
    alpha = mixture + -0.3*sc[1] + sc[6] + 5*sc[7] + 0.5*sc[10] + sc[11] + sc[14] + 2*sc[21];
    k_f = 1e-06 * 1.33e+13*exp(-785.46/tc[1]);
    redP = 1e-12 * alpha / k_f * 8.7e+42*exp(-7.5*tc[0]-2381.64/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0*exp(T/-1000))+ (1*exp(T/-645.4))+ (exp(-6844.3/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[37];
    Kc = 1.0 / (refC) * exp((g_RT[35] + g_RT[2]) - (g_RT[37]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[37] += 1 * qdot;

    /*reaction 172: I-C3H7 + O2 <=> C3H6 + HO2 */
    phi_f = sc[37]*sc[3];
    k_f = 1e-06 * 1.3e+11;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[8];
    Kc = exp((g_RT[37] + g_RT[3]) - (g_RT[35] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[37] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 173: N-C3H7 (+M) <=> CH3 + C2H4 (+M) */
    phi_f = sc[38];
    alpha = mixture;
    k_f = 1 * 1.23e+13*exp(-0.1*tc[0]-15204/tc[1]);
    redP = 1e-12 * alpha / k_f * 5.49e+49*exp(-10*tc[0]-18006.6/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((2.17*exp(T/-251))+ /*(-1.17*exp(T/-0))+*/ (exp(-1185/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[18];
    Kc = refC * exp((g_RT[38]) - (g_RT[15] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[38] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 174: H + C3H6 (+M) <=> N-C3H7 (+M) */
    phi_f = sc[2]*sc[35];
    alpha = mixture + -0.3*sc[1] + sc[6] + 5*sc[7] + 0.5*sc[10] + sc[11] + sc[14] + 2*sc[21];
    k_f = 1e-06 * 1.33e+13*exp(-1640.68/tc[1]);
    redP = 1e-12 * alpha / k_f * 6.26e+38*exp(-6.66*tc[0]-3523.14/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0*exp(T/-1000))+ (1*exp(T/-1310))+ (exp(-48097/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[38];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[35]) - (g_RT[38]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[38] += 1 * qdot;

    /*reaction 175: N-C3H7 + O2 <=> C3H6 + HO2 */
    phi_f = sc[38]*sc[3];
    k_f = 1e-06 * 9e+10;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[8];
    Kc = exp((g_RT[38] + g_RT[3]) - (g_RT[35] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[38] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[8] += 1 * qdot;

    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * qdot, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[39];                /*Gibbs free energy */
    double Kc;                      /*equilibrium constant */
    double k_f;                     /*forward reaction rate */
    double k_r;                     /*reverse reaction rate */
    double q_f;                     /*forward progress rate */
    double q_r;                     /*reverse progress rate */
    double phi_f;                   /*forward phase space factor */
    double phi_r;                   /*reverse phase space factor */
    double alpha;                   /*enhancement */
    double redP;                    /*reduced pressure */
    double logPred;                 /*log of above */
    double F;                       /*fallof rate enhancement */

    double F_troe;                  /*TROE intermediate */
    double logFcent;                /*TROE intermediate */
    double troe;                    /*TROE intermediate */
    double troe_c;                  /*TROE intermediate */
    double troe_n;                  /*TROE intermediate */

    double X;                       /*SRI intermediate */
    double F_sri;                   /*SRI intermediate */
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.314 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 39; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*reaction 1: H + O2 <=> OH + O */
    phi_f = sc[2]*sc[3];
    k_f = 1e-06 * 3.52e+16*exp(-0.7*tc[0]-8590.73/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = exp((g_RT[2] + g_RT[3]) - (g_RT[4] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: H2 + O <=> OH + H */
    phi_f = sc[6]*sc[5];
    k_f = 1e-06 * 50600*exp(2.67*tc[0]-3165.89/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[2];
    Kc = exp((g_RT[6] + g_RT[5]) - (g_RT[4] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[1] = q_f - q_r;

    /*reaction 3: H2 + OH <=> H2O + H */
    phi_f = sc[6]*sc[4];
    k_f = 1e-06 * 1.17e+09*exp(1.3*tc[0]-1829.53/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[2];
    Kc = exp((g_RT[6] + g_RT[4]) - (g_RT[7] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[2] = q_f - q_r;

    /*reaction 4: H2O + O <=> 2 OH */
    phi_f = sc[7]*sc[5];
    k_f = 1e-06 * 7.6*exp(3.84*tc[0]-6431.63/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = exp((g_RT[7] + g_RT[5]) - (2 * g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[3] = q_f - q_r;

    /*reaction 5: H + O + M <=> OH + M */
    phi_f = sc[2]*sc[5];
    alpha = mixture + 1.5*sc[6] + 11*sc[7] + 0.9*sc[10] + 2.8*sc[11];
    k_f = 1e-12 * alpha * 6.2e+16*exp(-0.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[5]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[4] = q_f - q_r;

    /*reaction 6: H2 + O2 <=> 2 OH */
    phi_f = sc[6]*sc[3];
    k_f = 1e-06 * 1.7e+13*exp(-24063/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = exp((g_RT[6] + g_RT[3]) - (2 * g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[5] = q_f - q_r;

    /*reaction 7: 2 H + M <=> H2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.5*sc[6] + 15.3*sc[7] + 0.9*sc[10] + 2.8*sc[11];
    k_f = 1e-12 * alpha * 7.2e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = 1.0 / (refC) * exp((2 * g_RT[2]) - (g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[6] = q_f - q_r;

    /*reaction 8: H + OH + M <=> H2O + M */
    phi_f = sc[2]*sc[4];
    alpha = mixture + 1.5*sc[6] + 11*sc[7] + 0.9*sc[10] + 2.8*sc[11];
    k_f = 1e-12 * alpha * 3.8e+22*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[7];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[4]) - (g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[7] = q_f - q_r;

    /*reaction 9: 2 O + M <=> O2 + M */
    phi_f = sc[5]*sc[5];
    alpha = mixture + 1.5*sc[6] + 11*sc[7] + 0.9*sc[10] + 2.8*sc[11];
    k_f = 1e-12 * alpha * 6.17e+15*exp(-0.5*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = 1.0 / (refC) * exp((2 * g_RT[5]) - (g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[8] = q_f - q_r;

    /*reaction 10: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[2]*sc[3];
    alpha = mixture + -0.5*sc[1] + -0.7*sc[3] + 6*sc[7] + -0.25*sc[10] + 0.5*sc[11] + 0.5*sc[21];
    k_f = 1e-06 * 4.65e+12*exp(0.44*tc[0]);
    redP = 1e-12 * alpha / k_f * 2.6e+19*exp(-1.2*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((1*exp(T/-345))+ (0*exp(T/-1))+ (exp(-345/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[8];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[3]) - (g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[9] = q_f - q_r;

    /*reaction 11: O + OH + M <=> HO2 + M */
    phi_f = sc[5]*sc[4];
    alpha = mixture;
    k_f = 1e-12 * alpha * 1e+16;
    q_f = phi_f * k_f;
    phi_r = sc[8];
    Kc = 1.0 / (refC) * exp((g_RT[5] + g_RT[4]) - (g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[10] = q_f - q_r;

    /*reaction 12: HO2 + H <=> 2 OH */
    phi_f = sc[8]*sc[2];
    k_f = 1e-06 * 7.08e+13*exp(-150.956/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = exp((g_RT[8] + g_RT[2]) - (2 * g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[11] = q_f - q_r;

    /*reaction 13: HO2 + H <=> H2 + O2 */
    phi_f = sc[8]*sc[2];
    k_f = 1e-06 * 4.28e+13*exp(-709.678/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = exp((g_RT[8] + g_RT[2]) - (g_RT[6] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[12] = q_f - q_r;

    /*reaction 14: HO2 + H <=> H2O + O */
    phi_f = sc[8]*sc[2];
    k_f = 1e-06 * 3.1e+13*exp(-866.049/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = exp((g_RT[8] + g_RT[2]) - (g_RT[7] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[13] = q_f - q_r;

    /*reaction 15: HO2 + O <=> OH + O2 */
    phi_f = sc[8]*sc[5];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = exp((g_RT[8] + g_RT[5]) - (g_RT[4] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[14] = q_f - q_r;

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    phi_f = sc[8]*sc[4];
    k_f = 1e-06 * 2.89e+13*exp(+250.191/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[3];
    Kc = exp((g_RT[8] + g_RT[4]) - (g_RT[7] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[15] = q_f - q_r;

    /*reaction 17: 2 OH (+M) <=> H2O2 (+M) */
    phi_f = sc[4]*sc[4];
    alpha = mixture + -0.3*sc[1] + sc[6] + 5*sc[7] + 0.5*sc[10] + sc[11] + sc[14] + 2*sc[21];
    k_f = 1e-06 * 7.4e+13*exp(-0.37*tc[0]);
    redP = 1e-12 * alpha / k_f * 1.34e+17*exp(-0.584*tc[0]+1154.74/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.265*exp(T/-94))+ (0.735*exp(T/-1756))+ (exp(-5182/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[9];
    Kc = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[16] = q_f - q_r;

    /*reaction 18: 2 HO2 <=> H2O2 + O2 */
    phi_f = sc[8]*sc[8];
    k_f = 1e-06 * 3.02e+12*exp(-697.65/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[3];
    Kc = exp((2 * g_RT[8]) - (g_RT[9] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[17] = q_f - q_r;

    /*reaction 19: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[9]*sc[2];
    k_f = 1e-06 * 4.79e+13*exp(-4005.48/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[6];
    Kc = exp((g_RT[9] + g_RT[2]) - (g_RT[8] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[18] = q_f - q_r;

    /*reaction 20: H2O2 + H <=> H2O + OH */
    phi_f = sc[9]*sc[2];
    k_f = 1e-06 * 1e+13*exp(-1804.27/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = exp((g_RT[9] + g_RT[2]) - (g_RT[7] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[19] = q_f - q_r;

    /*reaction 21: H2O2 + OH <=> H2O + HO2 */
    phi_f = sc[9]*sc[4];
    k_f = 1e-06 * 7.08e+12*exp(-721.706/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[8];
    Kc = exp((g_RT[9] + g_RT[4]) - (g_RT[7] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[20] = q_f - q_r;

    /*reaction 22: H2O2 + O <=> HO2 + OH */
    phi_f = sc[9]*sc[5];
    k_f = 1e-06 * 9.63e+06*exp(2*tc[0]-2008.76/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[4];
    Kc = exp((g_RT[9] + g_RT[5]) - (g_RT[8] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[21] = q_f - q_r;

    /*reaction 23: CO + OH <=> CO2 + H */
    phi_f = sc[10]*sc[4];
    k_f = 1e-06 * 4.4e+06*exp(1.5*tc[0]+372.884/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[2];
    Kc = exp((g_RT[10] + g_RT[4]) - (g_RT[11] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[22] = q_f - q_r;

    /*reaction 24: CO + HO2 <=> CO2 + OH */
    phi_f = sc[10]*sc[8];
    k_f = 1e-06 * 6.03e+13*exp(-11547.3/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[4];
    Kc = exp((g_RT[10] + g_RT[8]) - (g_RT[11] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[23] = q_f - q_r;

    /*reaction 25: HCO + M <=> CO + H + M */
    phi_f = sc[12];
    alpha = mixture + 0.9*sc[6] + 11*sc[7] + 1.5*sc[10] + 1.5*sc[11];
    k_f = 1e-06 * alpha * 1.86e+17*exp(-1*tc[0]-8555.85/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[2];
    Kc = refC * exp((g_RT[12]) - (g_RT[10] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[24] = q_f - q_r;

    /*reaction 26: HCO + H <=> CO + H2 */
    phi_f = sc[12]*sc[2];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[6];
    Kc = exp((g_RT[12] + g_RT[2]) - (g_RT[10] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[25] = q_f - q_r;

    /*reaction 27: HCO + O <=> CO + OH */
    phi_f = sc[12]*sc[5];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[4];
    Kc = exp((g_RT[12] + g_RT[5]) - (g_RT[10] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[26] = q_f - q_r;

    /*reaction 28: HCO + O <=> CO2 + H */
    phi_f = sc[12]*sc[5];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[2];
    Kc = exp((g_RT[12] + g_RT[5]) - (g_RT[11] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[27] = q_f - q_r;

    /*reaction 29: HCO + OH <=> CO + H2O */
    phi_f = sc[12]*sc[4];
    k_f = 1e-06 * 5.02e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[7];
    Kc = exp((g_RT[12] + g_RT[4]) - (g_RT[10] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[28] = q_f - q_r;

    /*reaction 30: HCO + O2 <=> CO + HO2 */
    phi_f = sc[12]*sc[3];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[8];
    Kc = exp((g_RT[12] + g_RT[3]) - (g_RT[10] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[29] = q_f - q_r;

    /*reaction 31: CH2O + M <=> HCO + H + M */
    phi_f = sc[13];
    alpha = mixture + 1.5*sc[6] + 15.3*sc[7] + 0.9*sc[10] + 2.8*sc[11];
    k_f = 1e-06 * alpha * 6.26e+16*exp(-39212.8/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[2];
    Kc = refC * exp((g_RT[13]) - (g_RT[12] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[30] = q_f - q_r;

    /*reaction 32: CH2O + H <=> HCO + H2 */
    phi_f = sc[13]*sc[2];
    k_f = 1e-06 * 1.26e+08*exp(1.62*tc[0]-1089.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[6];
    Kc = exp((g_RT[13] + g_RT[2]) - (g_RT[12] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[31] = q_f - q_r;

    /*reaction 33: CH2O + O <=> HCO + OH */
    phi_f = sc[13]*sc[5];
    k_f = 1e-06 * 3.5e+13*exp(-1768.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[4];
    Kc = exp((g_RT[13] + g_RT[5]) - (g_RT[12] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[32] = q_f - q_r;

    /*reaction 34: CH2O + OH <=> HCO + H2O */
    phi_f = sc[13]*sc[4];
    k_f = 1e-06 * 3.9e+10*exp(0.89*tc[0]-204.484/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[7];
    Kc = exp((g_RT[13] + g_RT[4]) - (g_RT[12] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[33] = q_f - q_r;

    /*reaction 35: CH4 + H <=> H2 + CH3 */
    phi_f = sc[14]*sc[2];
    k_f = 1e-06 * 13000*exp(3*tc[0]-4045.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[15];
    Kc = exp((g_RT[14] + g_RT[2]) - (g_RT[6] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[34] = q_f - q_r;

    /*reaction 36: CH4 + OH <=> H2O + CH3 */
    phi_f = sc[14]*sc[4];
    k_f = 1e-06 * 1.6e+07*exp(1.83*tc[0]-1400.12/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[15];
    Kc = exp((g_RT[14] + g_RT[4]) - (g_RT[7] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[35] = q_f - q_r;

    /*reaction 37: CH4 + O <=> CH3 + OH */
    phi_f = sc[14]*sc[5];
    k_f = 1e-06 * 1.9e+09*exp(1.44*tc[0]-4366.34/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[4];
    Kc = exp((g_RT[14] + g_RT[5]) - (g_RT[15] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[36] = q_f - q_r;

    /*reaction 38: CH4 + O2 <=> CH3 + HO2 */
    phi_f = sc[14]*sc[3];
    k_f = 1e-06 * 3.98e+13*exp(-28631.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[8];
    Kc = exp((g_RT[14] + g_RT[3]) - (g_RT[15] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[37] = q_f - q_r;

    /*reaction 39: CH4 + HO2 <=> CH3 + H2O2 */
    phi_f = sc[14]*sc[8];
    k_f = 1e-06 * 9.03e+12*exp(-12401.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[9];
    Kc = exp((g_RT[14] + g_RT[8]) - (g_RT[15] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[38] = q_f - q_r;

    /*reaction 40: CH3 + H <=> T-CH2 + H2 */
    phi_f = sc[15]*sc[2];
    k_f = 1e-06 * 1.8e+14*exp(-7601.99/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[6];
    Kc = exp((g_RT[15] + g_RT[2]) - (g_RT[16] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[39] = q_f - q_r;

    /*reaction 41: CH3 + H <=> S-CH2 + H2 */
    phi_f = sc[15]*sc[2];
    k_f = 1e-06 * 1.55e+14*exp(-6784.06/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[6];
    Kc = exp((g_RT[15] + g_RT[2]) - (g_RT[17] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[40] = q_f - q_r;

    /*reaction 42: CH3 + OH <=> S-CH2 + H2O */
    phi_f = sc[15]*sc[4];
    k_f = 1e-06 * 1e+13*exp(-1259.38/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[7];
    Kc = exp((g_RT[15] + g_RT[4]) - (g_RT[17] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[41] = q_f - q_r;

    /*reaction 43: CH3 + O <=> CH2O + H */
    phi_f = sc[15]*sc[5];
    k_f = 1e-06 * 8.43e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = exp((g_RT[15] + g_RT[5]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[42] = q_f - q_r;

    /*reaction 44: CH3 + T-CH2 <=> C2H4 + H */
    phi_f = sc[15]*sc[16];
    k_f = 1e-06 * 4.22e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[2];
    Kc = exp((g_RT[15] + g_RT[16]) - (g_RT[18] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[43] = q_f - q_r;

    /*reaction 45: CH3 + HO2 <=> CH3O + OH */
    phi_f = sc[15]*sc[8];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[4];
    Kc = exp((g_RT[15] + g_RT[8]) - (g_RT[19] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[44] = q_f - q_r;

    /*reaction 46: CH3 + O2 <=> CH2O + OH */
    phi_f = sc[15]*sc[3];
    k_f = 1e-06 * 3.3e+11*exp(-4499.85/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[4];
    Kc = exp((g_RT[15] + g_RT[3]) - (g_RT[13] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[45] = q_f - q_r;

    /*reaction 47: CH3 + O2 <=> CH3O + O */
    phi_f = sc[15]*sc[3];
    k_f = 1e-06 * 1.33e+14*exp(-15805.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[5];
    Kc = exp((g_RT[15] + g_RT[3]) - (g_RT[19] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[46] = q_f - q_r;

    /*reaction 48: 2 CH3 <=> C2H4 + H2 */
    phi_f = sc[15]*sc[15];
    k_f = 1e-06 * 1e+14*exp(-16106.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[6];
    Kc = exp((2 * g_RT[15]) - (g_RT[18] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[47] = q_f - q_r;

    /*reaction 49: 2 CH3 <=> C2H5 + H */
    phi_f = sc[15]*sc[15];
    k_f = 1e-06 * 3.16e+13*exp(-7397.51/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[2];
    Kc = exp((2 * g_RT[15]) - (g_RT[20] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[48] = q_f - q_r;

    /*reaction 50: CH3 + H (+M) <=> CH4 (+M) */
    phi_f = sc[15]*sc[2];
    alpha = mixture;
    k_f = 1e-06 * 2.11e+14;
    redP = 1e-12 * alpha / k_f * 6.26e+23*exp(-1.8*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.63*exp(T/-3315))+ (0.37*exp(T/-61)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[15] + g_RT[2]) - (g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[49] = q_f - q_r;

    /*reaction 51: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[15]*sc[15];
    alpha = mixture;
    k_f = 1e-06 * 1.81e+13;
    redP = 1e-12 * alpha / k_f * 1.27e+41*exp(-7*tc[0]-1390.49/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.38*exp(T/-73))+ (0.62*exp(T/-1180)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[21];
    Kc = 1.0 / (refC) * exp((2 * g_RT[15]) - (g_RT[21]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[50] = q_f - q_r;

    /*reaction 52: S-CH2 + OH <=> CH2O + H */
    phi_f = sc[17]*sc[4];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = exp((g_RT[17] + g_RT[4]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[51] = q_f - q_r;

    /*reaction 53: S-CH2 + O2 <=> CO + OH + H */
    phi_f = sc[17]*sc[3];
    k_f = 1e-06 * 3.13e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[4]*sc[2];
    Kc = refC * exp((g_RT[17] + g_RT[3]) - (g_RT[10] + g_RT[4] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[52] = q_f - q_r;

    /*reaction 54: S-CH2 + CO2 <=> CO + CH2O */
    phi_f = sc[17]*sc[11];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[13];
    Kc = exp((g_RT[17] + g_RT[11]) - (g_RT[10] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[53] = q_f - q_r;

    /*reaction 55: S-CH2 + M <=> T-CH2 + M */
    phi_f = sc[17];
    alpha = mixture + 1.4*sc[6] + 14.4*sc[7] + 0.8*sc[10] + 2.6*sc[11];
    k_f = 1e-06 * alpha * 6e+12;
    q_f = phi_f * k_f;
    phi_r = sc[16];
    Kc = exp((g_RT[17]) - (g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[54] = q_f - q_r;

    /*reaction 56: T-CH2 + H <=> CH + H2 */
    phi_f = sc[16]*sc[2];
    k_f = 1e-06 * 6.02e+12*exp(+899.728/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[22]*sc[6];
    Kc = exp((g_RT[16] + g_RT[2]) - (g_RT[22] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[55] = q_f - q_r;

    /*reaction 57: T-CH2 + OH <=> CH2O + H */
    phi_f = sc[16]*sc[4];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = exp((g_RT[16] + g_RT[4]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[56] = q_f - q_r;

    /*reaction 58: T-CH2 + OH <=> CH + H2O */
    phi_f = sc[16]*sc[4];
    k_f = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.57/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[22]*sc[7];
    Kc = exp((g_RT[16] + g_RT[4]) - (g_RT[22] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[57] = q_f - q_r;

    /*reaction 59: T-CH2 + O <=> CO + 2 H */
    phi_f = sc[16]*sc[5];
    k_f = 1e-06 * 8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[2]*sc[2];
    Kc = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[10] + 2 * g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[58] = q_f - q_r;

    /*reaction 60: T-CH2 + O <=> CO + H2 */
    phi_f = sc[16]*sc[5];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[6];
    Kc = exp((g_RT[16] + g_RT[5]) - (g_RT[10] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[59] = q_f - q_r;

    /*reaction 61: T-CH2 + O2 <=> CO2 + H2 */
    phi_f = sc[16]*sc[3];
    k_f = 1e-06 * 2.63e+12*exp(-750.579/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[6];
    Kc = exp((g_RT[16] + g_RT[3]) - (g_RT[11] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[60] = q_f - q_r;

    /*reaction 62: T-CH2 + O2 <=> CO + OH + H */
    phi_f = sc[16]*sc[3];
    k_f = 1e-06 * 6.58e+12*exp(-750.579/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[4]*sc[2];
    Kc = refC * exp((g_RT[16] + g_RT[3]) - (g_RT[10] + g_RT[4] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[61] = q_f - q_r;

    /*reaction 63: 2 T-CH2 <=> C2H2 + 2 H */
    phi_f = sc[16]*sc[16];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[2]*sc[2];
    Kc = refC * exp((2 * g_RT[16]) - (g_RT[23] + 2 * g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[62] = q_f - q_r;

    /*reaction 64: CH + O <=> CO + H */
    phi_f = sc[22]*sc[5];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[2];
    Kc = exp((g_RT[22] + g_RT[5]) - (g_RT[10] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[63] = q_f - q_r;

    /*reaction 65: CH + O2 <=> HCO + O */
    phi_f = sc[22]*sc[3];
    k_f = 1e-06 * 1.77e+11*exp(0.76*tc[0]+240.569/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[5];
    Kc = exp((g_RT[22] + g_RT[3]) - (g_RT[12] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[64] = q_f - q_r;

    /*reaction 66: CH + H2O <=> CH2O + H */
    phi_f = sc[22]*sc[7];
    k_f = 1e-06 * 1.17e+15*exp(-0.75*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = exp((g_RT[22] + g_RT[7]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[65] = q_f - q_r;

    /*reaction 67: CH + CO2 <=> HCO + CO */
    phi_f = sc[22]*sc[11];
    k_f = 1e-06 * 48*exp(3.22*tc[0]+1623.84/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[10];
    Kc = exp((g_RT[22] + g_RT[11]) - (g_RT[12] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[66] = q_f - q_r;

    /*reaction 68: CH3O + H <=> CH2O + H2 */
    phi_f = sc[19]*sc[2];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[6];
    Kc = exp((g_RT[19] + g_RT[2]) - (g_RT[13] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[67] = q_f - q_r;

    /*reaction 69: CH3O + H <=> S-CH2 + H2O */
    phi_f = sc[19]*sc[2];
    k_f = 1e-06 * 1.6e+13;
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[7];
    Kc = exp((g_RT[19] + g_RT[2]) - (g_RT[17] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[68] = q_f - q_r;

    /*reaction 70: CH3O + OH <=> CH2O + H2O */
    phi_f = sc[19]*sc[4];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[7];
    Kc = exp((g_RT[19] + g_RT[4]) - (g_RT[13] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[69] = q_f - q_r;

    /*reaction 71: CH3O + O <=> OH + CH2O */
    phi_f = sc[19]*sc[5];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = exp((g_RT[19] + g_RT[5]) - (g_RT[4] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[70] = q_f - q_r;

    /*reaction 72: CH3O + O2 <=> CH2O + HO2 */
    phi_f = sc[19]*sc[3];
    k_f = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1780.21/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[8];
    Kc = exp((g_RT[19] + g_RT[3]) - (g_RT[13] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[71] = q_f - q_r;

    /*reaction 73: CH3O + M <=> CH2O + H + M */
    phi_f = sc[19];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1e+13*exp(-6796.08/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = refC * exp((g_RT[19]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[72] = q_f - q_r;

    /*reaction 74: C2H6 + H <=> C2H5 + H2 */
    phi_f = sc[21]*sc[2];
    k_f = 1e-06 * 540*exp(3.5*tc[0]-2622.21/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[6];
    Kc = exp((g_RT[21] + g_RT[2]) - (g_RT[20] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[73] = q_f - q_r;

    /*reaction 75: C2H6 + O <=> C2H5 + OH */
    phi_f = sc[21]*sc[5];
    k_f = 1e-06 * 1.4*exp(4.3*tc[0]-1395.3/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[4];
    Kc = exp((g_RT[21] + g_RT[5]) - (g_RT[20] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[74] = q_f - q_r;

    /*reaction 76: C2H6 + OH <=> C2H5 + H2O */
    phi_f = sc[21]*sc[4];
    k_f = 1e-06 * 2.2e+07*exp(1.9*tc[0]-565.34/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[7];
    Kc = exp((g_RT[21] + g_RT[4]) - (g_RT[20] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[75] = q_f - q_r;

    /*reaction 77: C2H6 + CH3 <=> C2H5 + CH4 */
    phi_f = sc[21]*sc[15];
    k_f = 1e-06 * 0.55*exp(4*tc[0]-4173.88/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[14];
    Kc = exp((g_RT[21] + g_RT[15]) - (g_RT[20] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[76] = q_f - q_r;

    /*reaction 78: C2H6 (+M) <=> C2H5 + H (+M) */
    phi_f = sc[21];
    alpha = mixture;
    k_f = 1 * 8.85e+20*exp(-1.23*tc[0]-51445.8/tc[1]);
    redP = 1e-12 * alpha / k_f * 4.9e+42*exp(-6.43*tc[0]-53935.7/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.16*exp(T/-125))+ (0.84*exp(T/-2219))+ (exp(-6882/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[2];
    Kc = refC * exp((g_RT[21]) - (g_RT[20] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[77] = q_f - q_r;

    /*reaction 79: C2H5 + H <=> C2H4 + H2 */
    phi_f = sc[20]*sc[2];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[6];
    Kc = exp((g_RT[20] + g_RT[2]) - (g_RT[18] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[78] = q_f - q_r;

    /*reaction 80: C2H5 + O <=> C2H4 + OH */
    phi_f = sc[20]*sc[5];
    k_f = 1e-06 * 3.06e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[4];
    Kc = exp((g_RT[20] + g_RT[5]) - (g_RT[18] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[79] = q_f - q_r;

    /*reaction 81: C2H5 + O <=> CH3 + CH2O */
    phi_f = sc[20]*sc[5];
    k_f = 1e-06 * 4.24e+13;
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[13];
    Kc = exp((g_RT[20] + g_RT[5]) - (g_RT[15] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[80] = q_f - q_r;

    /*reaction 82: C2H5 + O2 <=> C2H4 + HO2 */
    phi_f = sc[20]*sc[3];
    k_f = 1e-06 * 2e+12*exp(-2513.95/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[8];
    Kc = exp((g_RT[20] + g_RT[3]) - (g_RT[18] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[81] = q_f - q_r;

    /*reaction 83: C2H5 (+M) <=> C2H4 + H (+M) */
    phi_f = sc[20];
    alpha = mixture;
    k_f = 1 * 1.11e+10*exp(1.037*tc[0]-18504.6/tc[1]);
    redP = 1e-12 * alpha / k_f * 3.99e+33*exp(-4.99*tc[0]-20130.9/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.832*exp(T/-1203))/*+ (0.168*exp(T/-0))*/);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[2];
    Kc = refC * exp((g_RT[20]) - (g_RT[18] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[82] = q_f - q_r;

    /*reaction 84: C2H4 + H <=> C2H3 + H2 */
    phi_f = sc[18]*sc[2];
    k_f = 1e-06 * 4.49e+07*exp(2.12*tc[0]-6723.92/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[6];
    Kc = exp((g_RT[18] + g_RT[2]) - (g_RT[24] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[83] = q_f - q_r;

    /*reaction 85: C2H4 + OH <=> C2H3 + H2O */
    phi_f = sc[18]*sc[4];
    k_f = 1e-06 * 553000*exp(2.31*tc[0]-1491.53/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[7];
    Kc = exp((g_RT[18] + g_RT[4]) - (g_RT[24] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[84] = q_f - q_r;

    /*reaction 86: C2H4 + O <=> CH3 + HCO */
    phi_f = sc[18]*sc[5];
    k_f = 1e-06 * 2.25e+06*exp(2.08*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[12];
    Kc = exp((g_RT[18] + g_RT[5]) - (g_RT[15] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[85] = q_f - q_r;

    /*reaction 87: C2H4 + O <=> CH2CHO + H */
    phi_f = sc[18]*sc[5];
    k_f = 1e-06 * 1.21e+06*exp(2.08*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[2];
    Kc = exp((g_RT[18] + g_RT[5]) - (g_RT[25] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[86] = q_f - q_r;

    /*reaction 88: 2 C2H4 <=> C2H3 + C2H5 */
    phi_f = sc[18]*sc[18];
    k_f = 1e-06 * 5.01e+14*exp(-32561.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[20];
    Kc = exp((2 * g_RT[18]) - (g_RT[24] + g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[87] = q_f - q_r;

    /*reaction 89: C2H4 + O2 <=> C2H3 + HO2 */
    phi_f = sc[18]*sc[3];
    k_f = 1e-06 * 4.22e+13*exp(-29000/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[8];
    Kc = exp((g_RT[18] + g_RT[3]) - (g_RT[24] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[88] = q_f - q_r;

    /*reaction 90: C2H4 + HO2 <=> C2H4O + OH */
    phi_f = sc[18]*sc[8];
    k_f = 1e-06 * 2.23e+12*exp(-8650.88/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[26]*sc[4];
    Kc = exp((g_RT[18] + g_RT[8]) - (g_RT[26] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[89] = q_f - q_r;

    /*reaction 91: C2H4O + HO2 <=> CH3 + CO + H2O2 */
    phi_f = sc[26]*sc[8];
    k_f = 1e-06 * 4e+12*exp(-8559.46/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[10]*sc[9];
    Kc = refC * exp((g_RT[26] + g_RT[8]) - (g_RT[15] + g_RT[10] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[90] = q_f - q_r;

    /*reaction 92: C2H4 + M <=> C2H3 + H + M */
    phi_f = sc[18];
    alpha = mixture;
    k_f = 1e-06 * alpha * 2.6e+17*exp(-48600/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[2];
    Kc = refC * exp((g_RT[18]) - (g_RT[24] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[91] = q_f - q_r;

    /*reaction 93: C2H4 + M <=> C2H2 + H2 + M */
    phi_f = sc[18];
    alpha = mixture;
    k_f = 1e-06 * alpha * 3.5e+16*exp(-36000/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[6];
    Kc = refC * exp((g_RT[18]) - (g_RT[23] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[92] = q_f - q_r;

    /*reaction 94: C2H3 + H <=> C2H2 + H2 */
    phi_f = sc[24]*sc[2];
    k_f = 1e-06 * 1.21e+13;
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[6];
    Kc = exp((g_RT[24] + g_RT[2]) - (g_RT[23] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[93] = q_f - q_r;

    /*reaction 95: C2H3 (+M) <=> C2H2 + H (+M) */
    phi_f = sc[24];
    alpha = mixture;
    k_f = 1 * 6.38e+09*exp(1*tc[0]-18936.4/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.51e+14*exp(0.1*tc[0]-16450.1/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.7*exp(T/-1e+30))+ (0.3*exp(T/-1e-30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[2];
    Kc = refC * exp((g_RT[24]) - (g_RT[23] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[94] = q_f - q_r;

    /*reaction 96: C2H3 + O2 <=> CH2O + HCO */
    phi_f = sc[24]*sc[3];
    k_f = 1e-06 * 1.7e+29*exp(-5.312*tc[0]-3272.83/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[12];
    Kc = exp((g_RT[24] + g_RT[3]) - (g_RT[13] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[95] = q_f - q_r;

    /*reaction 97: C2H3 + O2 <=> CH2CHO + O */
    phi_f = sc[24]*sc[3];
    k_f = 1e-06 * 7e+14*exp(-0.611*tc[0]-2648.43/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[5];
    Kc = exp((g_RT[24] + g_RT[3]) - (g_RT[25] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[96] = q_f - q_r;

    /*reaction 98: C2H3 + O2 <=> C2H2 + HO2 */
    phi_f = sc[24]*sc[3];
    k_f = 1e-06 * 5.19e+15*exp(-1.26*tc[0]-1667.15/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[8];
    Kc = exp((g_RT[24] + g_RT[3]) - (g_RT[23] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[97] = q_f - q_r;

    /*reaction 99: CH2CHO <=> CH2CO + H */
    phi_f = sc[25];
    k_f = 1 * 1.047e+37*exp(-7.189*tc[0]-22315.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[2];
    Kc = refC * exp((g_RT[25]) - (g_RT[27] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[98] = q_f - q_r;

    /*reaction 100: C2H2 + O <=> HCCO + H */
    phi_f = sc[23]*sc[5];
    k_f = 1e-06 * 4e+14*exp(-5364.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[2];
    Kc = exp((g_RT[23] + g_RT[5]) - (g_RT[28] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[99] = q_f - q_r;

    /*reaction 101: C2H2 + O <=> T-CH2 + CO */
    phi_f = sc[23]*sc[5];
    k_f = 1e-06 * 1.6e+14*exp(-4979.79/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[10];
    Kc = exp((g_RT[23] + g_RT[5]) - (g_RT[16] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[100] = q_f - q_r;

    /*reaction 102: C2H2 + O2 <=> CH2O + CO */
    phi_f = sc[23]*sc[3];
    k_f = 1e-06 * 4.6e+15*exp(-0.54*tc[0]-22613.5/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[10];
    Kc = exp((g_RT[23] + g_RT[3]) - (g_RT[13] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[101] = q_f - q_r;

    /*reaction 103: C2H2 + OH <=> CH2CO + H */
    phi_f = sc[23]*sc[4];
    k_f = 1e-06 * 1.9e+07*exp(1.7*tc[0]-502.788/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[2];
    Kc = exp((g_RT[23] + g_RT[4]) - (g_RT[27] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[102] = q_f - q_r;

    /*reaction 104: C2H2 + OH <=> C2H + H2O */
    phi_f = sc[23]*sc[4];
    k_f = 1e-06 * 3.37e+07*exp(2*tc[0]-7046.28/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[7];
    Kc = exp((g_RT[23] + g_RT[4]) - (g_RT[29] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[103] = q_f - q_r;

    /*reaction 105: CH2CO + H <=> CH3 + CO */
    phi_f = sc[27]*sc[2];
    k_f = 1e-06 * 1.5e+09*exp(1.43*tc[0]-1353.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[10];
    Kc = exp((g_RT[27] + g_RT[2]) - (g_RT[15] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[104] = q_f - q_r;

    /*reaction 106: CH2CO + O <=> T-CH2 + CO2 */
    phi_f = sc[27]*sc[5];
    k_f = 1e-06 * 2e+13*exp(-1154.74/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[11];
    Kc = exp((g_RT[27] + g_RT[5]) - (g_RT[16] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[105] = q_f - q_r;

    /*reaction 107: CH2CO + O <=> HCCO + OH */
    phi_f = sc[27]*sc[5];
    k_f = 1e-06 * 1e+13*exp(-1006.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[4];
    Kc = exp((g_RT[27] + g_RT[5]) - (g_RT[28] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[106] = q_f - q_r;

    /*reaction 108: CH2CO + CH3 <=> C2H5 + CO */
    phi_f = sc[27]*sc[15];
    k_f = 1e-06 * 9e+10;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[10];
    Kc = exp((g_RT[27] + g_RT[15]) - (g_RT[20] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[107] = q_f - q_r;

    /*reaction 109: HCCO + H <=> S-CH2 + CO */
    phi_f = sc[28]*sc[2];
    k_f = 1e-06 * 1.5e+14;
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[10];
    Kc = exp((g_RT[28] + g_RT[2]) - (g_RT[17] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[108] = q_f - q_r;

    /*reaction 110: HCCO + OH <=> HCO + CO + H */
    phi_f = sc[28]*sc[4];
    k_f = 1e-06 * 2e+12;
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[10]*sc[2];
    Kc = refC * exp((g_RT[28] + g_RT[4]) - (g_RT[12] + g_RT[10] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[109] = q_f - q_r;

    /*reaction 111: HCCO + O <=> 2 CO + H */
    phi_f = sc[28]*sc[5];
    k_f = 1e-06 * 9.64e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[10]*sc[2];
    Kc = refC * exp((g_RT[28] + g_RT[5]) - (2 * g_RT[10] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[110] = q_f - q_r;

    /*reaction 112: HCCO + O2 <=> 2 CO + OH */
    phi_f = sc[28]*sc[3];
    k_f = 1e-06 * 2.88e+07*exp(1.7*tc[0]-503.991/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[10]*sc[4];
    Kc = refC * exp((g_RT[28] + g_RT[3]) - (2 * g_RT[10] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[111] = q_f - q_r;

    /*reaction 113: HCCO + O2 <=> CO2 + CO + H */
    phi_f = sc[28]*sc[3];
    k_f = 1e-06 * 1.4e+07*exp(1.7*tc[0]-503.991/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[10]*sc[2];
    Kc = refC * exp((g_RT[28] + g_RT[3]) - (g_RT[11] + g_RT[10] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[112] = q_f - q_r;

    /*reaction 114: C2H + OH <=> HCCO + H */
    phi_f = sc[29]*sc[4];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[2];
    Kc = exp((g_RT[29] + g_RT[4]) - (g_RT[28] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[113] = q_f - q_r;

    /*reaction 115: C2H + O <=> CO + CH */
    phi_f = sc[29]*sc[5];
    k_f = 1e-06 * 1.02e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[22];
    Kc = exp((g_RT[29] + g_RT[5]) - (g_RT[10] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[114] = q_f - q_r;

    /*reaction 116: C2H + O2 <=> HCCO + O */
    phi_f = sc[29]*sc[3];
    k_f = 1e-06 * 6.02e+11;
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[5];
    Kc = exp((g_RT[29] + g_RT[3]) - (g_RT[28] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[115] = q_f - q_r;

    /*reaction 117: C2H + O2 <=> CH + CO2 */
    phi_f = sc[29]*sc[3];
    k_f = 1e-06 * 4.5e+15*exp(-12629.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[22]*sc[11];
    Kc = exp((g_RT[29] + g_RT[3]) - (g_RT[22] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[116] = q_f - q_r;

    /*reaction 118: C2H + O2 <=> HCO + CO */
    phi_f = sc[29]*sc[3];
    k_f = 1e-06 * 2.41e+12;
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[10];
    Kc = exp((g_RT[29] + g_RT[3]) - (g_RT[12] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[117] = q_f - q_r;

    /*reaction 119: CH2OH + H <=> CH2O + H2 */
    phi_f = sc[30]*sc[2];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[6];
    Kc = exp((g_RT[30] + g_RT[2]) - (g_RT[13] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[118] = q_f - q_r;

    /*reaction 120: CH2OH + H <=> CH3 + OH */
    phi_f = sc[30]*sc[2];
    k_f = 1e-06 * 1.75e+14*exp(-1407.33/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[4];
    Kc = exp((g_RT[30] + g_RT[2]) - (g_RT[15] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[119] = q_f - q_r;

    /*reaction 121: CH2OH + OH <=> CH2O + H2O */
    phi_f = sc[30]*sc[4];
    k_f = 1e-06 * 2.4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[7];
    Kc = exp((g_RT[30] + g_RT[4]) - (g_RT[13] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[120] = q_f - q_r;

    /*reaction 122: CH2OH + O2 <=> CH2O + HO2 */
    phi_f = sc[30]*sc[3];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[8];
    Kc = exp((g_RT[30] + g_RT[3]) - (g_RT[13] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[121] = q_f - q_r;

    /*reaction 123: CH2OH + M <=> CH2O + H + M */
    phi_f = sc[30];
    alpha = mixture + 1.4*sc[6] + 14.4*sc[7] + 0.8*sc[10] + 2.6*sc[11];
    k_f = 1e-06 * alpha * 5e+13*exp(-12641.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = refC * exp((g_RT[30]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[122] = q_f - q_r;

    /*reaction 124: CH3O + M <=> CH2OH + M */
    phi_f = sc[19];
    alpha = mixture + 1.5*sc[6] + 11*sc[7] + 0.9*sc[10] + 2.8*sc[11];
    k_f = 1e-06 * alpha * 1e+14*exp(-9622.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30];
    Kc = exp((g_RT[19]) - (g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[123] = q_f - q_r;

    /*reaction 125: CH2CO + OH <=> CH2OH + CO */
    phi_f = sc[27]*sc[4];
    k_f = 1e-06 * 1.02e+13;
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[10];
    Kc = exp((g_RT[27] + g_RT[4]) - (g_RT[30] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[124] = q_f - q_r;

    /*reaction 126: CH3OH + OH <=> CH2OH + H2O */
    phi_f = sc[31]*sc[4];
    k_f = 1e-06 * 1.44e+06*exp(2*tc[0]+422.199/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[7];
    Kc = exp((g_RT[31] + g_RT[4]) - (g_RT[30] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[125] = q_f - q_r;

    /*reaction 127: CH3OH + OH <=> CH3O + H2O */
    phi_f = sc[31]*sc[4];
    k_f = 1e-06 * 6.3e+06*exp(2*tc[0]-757.796/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[7];
    Kc = exp((g_RT[31] + g_RT[4]) - (g_RT[19] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[126] = q_f - q_r;

    /*reaction 128: CH3OH + H <=> CH2OH + H2 */
    phi_f = sc[31]*sc[2];
    k_f = 1e-06 * 1.64e+07*exp(2*tc[0]-2273.38/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[6];
    Kc = exp((g_RT[31] + g_RT[2]) - (g_RT[30] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[127] = q_f - q_r;

    /*reaction 129: CH3OH + H <=> CH3O + H2 */
    phi_f = sc[31]*sc[2];
    k_f = 1e-06 * 3.83e+07*exp(2*tc[0]-2946.98/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[6];
    Kc = exp((g_RT[31] + g_RT[2]) - (g_RT[19] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[128] = q_f - q_r;

    /*reaction 130: CH3OH + O <=> CH2OH + OH */
    phi_f = sc[31]*sc[5];
    k_f = 1e-06 * 1e+13*exp(-2357.58/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[4];
    Kc = exp((g_RT[31] + g_RT[5]) - (g_RT[30] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[129] = q_f - q_r;

    /*reaction 131: CH3OH + HO2 <=> CH2OH + H2O2 */
    phi_f = sc[31]*sc[8];
    k_f = 1e-06 * 6.2e+12*exp(-9755.09/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[9];
    Kc = exp((g_RT[31] + g_RT[8]) - (g_RT[30] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[130] = q_f - q_r;

    /*reaction 132: CH3OH + O2 <=> CH2OH + HO2 */
    phi_f = sc[31]*sc[3];
    k_f = 1e-06 * 2e+13*exp(-22613.5/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[8];
    Kc = exp((g_RT[31] + g_RT[3]) - (g_RT[30] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[131] = q_f - q_r;

    /*reaction 133: C3H4 + O <=> C2H4 + CO */
    phi_f = sc[32]*sc[5];
    k_f = 1e-06 * 2e+07*exp(1.8*tc[0]-503.271/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[10];
    Kc = exp((g_RT[32] + g_RT[5]) - (g_RT[18] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[132] = q_f - q_r;

    /*reaction 134: CH3 + C2H2 <=> C3H4 + H */
    phi_f = sc[15]*sc[23];
    k_f = 1e-06 * 2.56e+09*exp(1.1*tc[0]-6866.57/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[2];
    Kc = exp((g_RT[15] + g_RT[23]) - (g_RT[32] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[133] = q_f - q_r;

    /*reaction 135: C3H4 + O <=> HCCO + CH3 */
    phi_f = sc[32]*sc[5];
    k_f = 1e-06 * 7.3e+12*exp(-1132.36/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[15];
    Kc = exp((g_RT[32] + g_RT[5]) - (g_RT[28] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[134] = q_f - q_r;

    /*reaction 136: C3H3 + H (+M) <=> C3H4 (+M) */
    phi_f = sc[33]*sc[2];
    alpha = mixture;
    k_f = 1e-06 * 3e+13;
    redP = 1e-12 * alpha / k_f * 9e+15*exp(1*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-1e+30))+ (0.5*exp(T/-1e-30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[32];
    Kc = 1.0 / (refC) * exp((g_RT[33] + g_RT[2]) - (g_RT[32]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[135] = q_f - q_r;

    /*reaction 137: C3H3 + HO2 <=> C3H4 + O2 */
    phi_f = sc[33]*sc[8];
    k_f = 1e-06 * 2.5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[3];
    Kc = exp((g_RT[33] + g_RT[8]) - (g_RT[32] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[136] = q_f - q_r;

    /*reaction 138: C3H4 + OH <=> C3H3 + H2O */
    phi_f = sc[32]*sc[4];
    k_f = 1e-06 * 5.3e+06*exp(2*tc[0]-1006.54/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[33]*sc[7];
    Kc = exp((g_RT[32] + g_RT[4]) - (g_RT[33] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[137] = q_f - q_r;

    /*reaction 139: C3H3 + O2 <=> CH2CO + HCO */
    phi_f = sc[33]*sc[3];
    k_f = 1e-06 * 3e+10*exp(-1443.42/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[12];
    Kc = exp((g_RT[33] + g_RT[3]) - (g_RT[27] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[138] = q_f - q_r;

    /*reaction 140: C3H4 + H (+M) <=> C3H5 (+M) */
    phi_f = sc[32]*sc[2];
    alpha = mixture;
    k_f = 1e-06 * 4e+13;
    redP = 1e-12 * alpha / k_f * 3e+24*exp(-2*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2*exp(T/-1e+30))+ (0.8*exp(T/-1e-30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[34];
    Kc = 1.0 / (refC) * exp((g_RT[32] + g_RT[2]) - (g_RT[34]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[139] = q_f - q_r;

    /*reaction 141: C3H5 + H <=> C3H4 + H2 */
    phi_f = sc[34]*sc[2];
    k_f = 1e-06 * 1.8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[6];
    Kc = exp((g_RT[34] + g_RT[2]) - (g_RT[32] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[140] = q_f - q_r;

    /*reaction 142: C3H5 + O2 <=> C3H4 + HO2 */
    phi_f = sc[34]*sc[3];
    k_f = 1e-06 * 4.99e+15*exp(-1.4*tc[0]-11287.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[8];
    Kc = exp((g_RT[34] + g_RT[3]) - (g_RT[32] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[141] = q_f - q_r;

    /*reaction 143: C3H5 + CH3 <=> C3H4 + CH4 */
    phi_f = sc[34]*sc[15];
    k_f = 1e-06 * 3e+12*exp(-0.32*tc[0]+65.9185/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[14];
    Kc = exp((g_RT[34] + g_RT[15]) - (g_RT[32] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[142] = q_f - q_r;

    /*reaction 144: C2H2 + CH3 (+M) <=> C3H5 (+M) */
    phi_f = sc[23]*sc[15];
    alpha = mixture;
    k_f = 1e-06 * 6e+08;
    redP = 1e-12 * alpha / k_f * 2e+09*exp(1*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-1e+30))+ (0.5*exp(T/-1e-30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[34];
    Kc = 1.0 / (refC) * exp((g_RT[23] + g_RT[15]) - (g_RT[34]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[143] = q_f - q_r;

    /*reaction 145: C3H5 + OH <=> C3H4 + H2O */
    phi_f = sc[34]*sc[4];
    k_f = 1e-06 * 6e+12;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[7];
    Kc = exp((g_RT[34] + g_RT[4]) - (g_RT[32] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[144] = q_f - q_r;

    /*reaction 146: C3H3 + HCO <=> C3H4 + CO */
    phi_f = sc[33]*sc[12];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[10];
    Kc = exp((g_RT[33] + g_RT[12]) - (g_RT[32] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[145] = q_f - q_r;

    /*reaction 147: C3H3 + HO2 <=> OH + CO + C2H3 */
    phi_f = sc[33]*sc[8];
    k_f = 1e-06 * 8e+11;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[10]*sc[24];
    Kc = refC * exp((g_RT[33] + g_RT[8]) - (g_RT[4] + g_RT[10] + g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[146] = q_f - q_r;

    /*reaction 148: C3H4 + O2 <=> CH3 + HCO + CO */
    phi_f = sc[32]*sc[3];
    k_f = 1e-06 * 4e+14*exp(-21049.8/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[12]*sc[10];
    Kc = refC * exp((g_RT[32] + g_RT[3]) - (g_RT[15] + g_RT[12] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[147] = q_f - q_r;

    /*reaction 149: C3H6 + O <=> C2H5 + HCO */
    phi_f = sc[35]*sc[5];
    k_f = 1e-06 * 3.5e+07*exp(1.65*tc[0]+489.557/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[12];
    Kc = exp((g_RT[35] + g_RT[5]) - (g_RT[20] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[148] = q_f - q_r;

    /*reaction 150: C3H6 + OH <=> C3H5 + H2O */
    phi_f = sc[35]*sc[4];
    k_f = 1e-06 * 3.1e+06*exp(2*tc[0]+150.116/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[7];
    Kc = exp((g_RT[35] + g_RT[4]) - (g_RT[34] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[149] = q_f - q_r;

    /*reaction 151: C3H6 + O <=> CH2CO + CH3 + H */
    phi_f = sc[35]*sc[5];
    k_f = 1e-06 * 1.2e+08*exp(1.65*tc[0]-164.791/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[15]*sc[2];
    Kc = refC * exp((g_RT[35] + g_RT[5]) - (g_RT[27] + g_RT[15] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[150] = q_f - q_r;

    /*reaction 152: C3H6 + H <=> C3H5 + H2 */
    phi_f = sc[35]*sc[2];
    k_f = 1e-06 * 170000*exp(2.5*tc[0]-1254.57/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[6];
    Kc = exp((g_RT[35] + g_RT[2]) - (g_RT[34] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[151] = q_f - q_r;

    /*reaction 153: C3H5 + H (+M) <=> C3H6 (+M) */
    phi_f = sc[34]*sc[2];
    alpha = mixture + -0.3*sc[1] + sc[6] + 5*sc[7] + 0.5*sc[10] + sc[11] + sc[14] + 2*sc[21];
    k_f = 1e-06 * 2e+14;
    redP = 1e-12 * alpha / k_f * 1.33e+60*exp(-12*tc[0]-3003.51/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.98*exp(T/-1097))+ (0.02*exp(T/-1097))+ (exp(-6860/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[35];
    Kc = 1.0 / (refC) * exp((g_RT[34] + g_RT[2]) - (g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[152] = q_f - q_r;

    /*reaction 154: C3H5 + HO2 <=> C3H6 + O2 */
    phi_f = sc[34]*sc[8];
    k_f = 1e-06 * 2.66e+12;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[3];
    Kc = exp((g_RT[34] + g_RT[8]) - (g_RT[35] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[153] = q_f - q_r;

    /*reaction 155: C3H5 + HO2 <=> OH + C2H3 + CH2O */
    phi_f = sc[34]*sc[8];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[24]*sc[13];
    Kc = refC * exp((g_RT[34] + g_RT[8]) - (g_RT[4] + g_RT[24] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[154] = q_f - q_r;

    /*reaction 156: C2H3 + CH3 (+M) <=> C3H6 (+M) */
    phi_f = sc[24]*sc[15];
    alpha = mixture + -0.3*sc[1] + sc[6] + 5*sc[7] + 0.5*sc[10] + sc[11] + sc[14] + 2*sc[21];
    k_f = 1e-06 * 2.5e+13;
    redP = 1e-12 * alpha / k_f * 4.27e+58*exp(-11.94*tc[0]-4917.24/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.825*exp(T/-1341))+ (0.175*exp(T/-60000))+ (exp(-10140/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[35];
    Kc = 1.0 / (refC) * exp((g_RT[24] + g_RT[15]) - (g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[155] = q_f - q_r;

    /*reaction 157: C3H6 + H <=> C2H4 + CH3 */
    phi_f = sc[35]*sc[2];
    k_f = 1e-06 * 1.6e+22*exp(-2.39*tc[0]-5629.33/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[15];
    Kc = exp((g_RT[35] + g_RT[2]) - (g_RT[18] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[156] = q_f - q_r;

    /*reaction 158: CH3 + C2H3 <=> C3H5 + H */
    phi_f = sc[15]*sc[24];
    k_f = 1e-06 * 1.5e+24*exp(-2.83*tc[0]-9370.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[2];
    Kc = exp((g_RT[15] + g_RT[24]) - (g_RT[34] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[157] = q_f - q_r;

    /*reaction 159: C3H8 (+M) <=> CH3 + C2H5 (+M) */
    phi_f = sc[36];
    alpha = mixture;
    k_f = 1 * 1.1e+17*exp(-42472.5/tc[1]);
    redP = 1e-12 * alpha / k_f * 7.83e+18*exp(-32701.8/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.24*exp(T/-1946))+ (0.76*exp(T/-38)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[20];
    Kc = refC * exp((g_RT[36]) - (g_RT[15] + g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[158] = q_f - q_r;

    /*reaction 160: C3H8 + O2 <=> I-C3H7 + HO2 */
    phi_f = sc[36]*sc[3];
    k_f = 1e-06 * 4e+13*exp(-23905.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[8];
    Kc = exp((g_RT[36] + g_RT[3]) - (g_RT[37] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[159] = q_f - q_r;

    /*reaction 161: C3H8 + O2 <=> N-C3H7 + HO2 */
    phi_f = sc[36]*sc[3];
    k_f = 1e-06 * 4e+13*exp(-25632.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[8];
    Kc = exp((g_RT[36] + g_RT[3]) - (g_RT[38] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[160] = q_f - q_r;

    /*reaction 162: C3H8 + H <=> I-C3H7 + H2 */
    phi_f = sc[36]*sc[2];
    k_f = 1e-06 * 1.3e+06*exp(2.4*tc[0]-2250.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[6];
    Kc = exp((g_RT[36] + g_RT[2]) - (g_RT[37] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[161] = q_f - q_r;

    /*reaction 163: C3H8 + H <=> N-C3H7 + H2 */
    phi_f = sc[36]*sc[2];
    k_f = 1e-06 * 1.33e+06*exp(2.54*tc[0]-3402.85/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[6];
    Kc = exp((g_RT[36] + g_RT[2]) - (g_RT[38] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[162] = q_f - q_r;

    /*reaction 164: C3H8 + O <=> I-C3H7 + OH */
    phi_f = sc[36]*sc[5];
    k_f = 1e-06 * 47600*exp(2.71*tc[0]-1060.55/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[4];
    Kc = exp((g_RT[36] + g_RT[5]) - (g_RT[37] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[163] = q_f - q_r;

    /*reaction 165: C3H8 + O <=> N-C3H7 + OH */
    phi_f = sc[36]*sc[5];
    k_f = 1e-06 * 190000*exp(2.68*tc[0]-1871.39/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[4];
    Kc = exp((g_RT[36] + g_RT[5]) - (g_RT[38] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[164] = q_f - q_r;

    /*reaction 166: C3H8 + OH <=> N-C3H7 + H2O */
    phi_f = sc[36]*sc[4];
    k_f = 1e-06 * 1400*exp(2.66*tc[0]-265.35/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[7];
    Kc = exp((g_RT[36] + g_RT[4]) - (g_RT[38] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[165] = q_f - q_r;

    /*reaction 167: C3H8 + OH <=> I-C3H7 + H2O */
    phi_f = sc[36]*sc[4];
    k_f = 1e-06 * 27000*exp(2.39*tc[0]-197.866/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[7];
    Kc = exp((g_RT[36] + g_RT[4]) - (g_RT[37] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[166] = q_f - q_r;

    /*reaction 168: C3H8 + HO2 <=> I-C3H7 + H2O2 */
    phi_f = sc[36]*sc[8];
    k_f = 1e-06 * 9640*exp(2.6*tc[0]-7004.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[9];
    Kc = exp((g_RT[36] + g_RT[8]) - (g_RT[37] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[167] = q_f - q_r;

    /*reaction 169: C3H8 + HO2 <=> N-C3H7 + H2O2 */
    phi_f = sc[36]*sc[8];
    k_f = 1e-06 * 47600*exp(2.55*tc[0]-8299.65/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[9];
    Kc = exp((g_RT[36] + g_RT[8]) - (g_RT[38] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[168] = q_f - q_r;

    /*reaction 170: I-C3H7 + C3H8 <=> N-C3H7 + C3H8 */
    phi_f = sc[37]*sc[36];
    k_f = 1e-06 * 0.0084*exp(4.2*tc[0]-4366.34/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[36];
    Kc = exp((g_RT[37] + g_RT[36]) - (g_RT[38] + g_RT[36]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[169] = q_f - q_r;

    /*reaction 171: C3H6 + H (+M) <=> I-C3H7 (+M) */
    phi_f = sc[35]*sc[2];
    alpha = mixture + -0.3*sc[1] + sc[6] + 5*sc[7] + 0.5*sc[10] + sc[11] + sc[14] + 2*sc[21];
    k_f = 1e-06 * 1.33e+13*exp(-785.46/tc[1]);
    redP = 1e-12 * alpha / k_f * 8.7e+42*exp(-7.5*tc[0]-2381.64/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0*exp(T/-1000))+ (1*exp(T/-645.4))+ (exp(-6844.3/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[37];
    Kc = 1.0 / (refC) * exp((g_RT[35] + g_RT[2]) - (g_RT[37]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[170] = q_f - q_r;

    /*reaction 172: I-C3H7 + O2 <=> C3H6 + HO2 */
    phi_f = sc[37]*sc[3];
    k_f = 1e-06 * 1.3e+11;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[8];
    Kc = exp((g_RT[37] + g_RT[3]) - (g_RT[35] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[171] = q_f - q_r;

    /*reaction 173: N-C3H7 (+M) <=> CH3 + C2H4 (+M) */
    phi_f = sc[38];
    alpha = mixture;
    k_f = 1 * 1.23e+13*exp(-0.1*tc[0]-15204/tc[1]);
    redP = 1e-12 * alpha / k_f * 5.49e+49*exp(-10*tc[0]-18006.6/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((2.17*exp(T/-251))+ /*(-1.17*exp(T/-0))+*/ (exp(-1185/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[18];
    Kc = refC * exp((g_RT[38]) - (g_RT[15] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[172] = q_f - q_r;

    /*reaction 174: H + C3H6 (+M) <=> N-C3H7 (+M) */
    phi_f = sc[2]*sc[35];
    alpha = mixture + -0.3*sc[1] + sc[6] + 5*sc[7] + 0.5*sc[10] + sc[11] + sc[14] + 2*sc[21];
    k_f = 1e-06 * 1.33e+13*exp(-1640.68/tc[1]);
    redP = 1e-12 * alpha / k_f * 6.26e+38*exp(-6.66*tc[0]-3523.14/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0*exp(T/-1000))+ (1*exp(T/-1310))+ (exp(-48097/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[38];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[35]) - (g_RT[38]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[173] = q_f - q_r;

    /*reaction 175: N-C3H7 + O2 <=> C3H6 + HO2 */
    phi_f = sc[38]*sc[3];
    k_f = 1e-06 * 9e+10;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[8];
    Kc = exp((g_RT[38] + g_RT[3]) - (g_RT[35] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[174] = q_f - q_r;

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *kc, double * g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.314 / T;

    /*reaction 1: H + O2 <=> OH + O */
    kc[0] = exp((g_RT[2] + g_RT[3]) - (g_RT[4] + g_RT[5]));

    /*reaction 2: H2 + O <=> OH + H */
    kc[1] = exp((g_RT[6] + g_RT[5]) - (g_RT[4] + g_RT[2]));

    /*reaction 3: H2 + OH <=> H2O + H */
    kc[2] = exp((g_RT[6] + g_RT[4]) - (g_RT[7] + g_RT[2]));

    /*reaction 4: H2O + O <=> 2 OH */
    kc[3] = exp((g_RT[7] + g_RT[5]) - (2 * g_RT[4]));

    /*reaction 5: H + O + M <=> OH + M */
    kc[4] = 1.0 / (refC) * exp((g_RT[2] + g_RT[5]) - (g_RT[4]));

    /*reaction 6: H2 + O2 <=> 2 OH */
    kc[5] = exp((g_RT[6] + g_RT[3]) - (2 * g_RT[4]));

    /*reaction 7: 2 H + M <=> H2 + M */
    kc[6] = 1.0 / (refC) * exp((2 * g_RT[2]) - (g_RT[6]));

    /*reaction 8: H + OH + M <=> H2O + M */
    kc[7] = 1.0 / (refC) * exp((g_RT[2] + g_RT[4]) - (g_RT[7]));

    /*reaction 9: 2 O + M <=> O2 + M */
    kc[8] = 1.0 / (refC) * exp((2 * g_RT[5]) - (g_RT[3]));

    /*reaction 10: H + O2 (+M) <=> HO2 (+M) */
    kc[9] = 1.0 / (refC) * exp((g_RT[2] + g_RT[3]) - (g_RT[8]));

    /*reaction 11: O + OH + M <=> HO2 + M */
    kc[10] = 1.0 / (refC) * exp((g_RT[5] + g_RT[4]) - (g_RT[8]));

    /*reaction 12: HO2 + H <=> 2 OH */
    kc[11] = exp((g_RT[8] + g_RT[2]) - (2 * g_RT[4]));

    /*reaction 13: HO2 + H <=> H2 + O2 */
    kc[12] = exp((g_RT[8] + g_RT[2]) - (g_RT[6] + g_RT[3]));

    /*reaction 14: HO2 + H <=> H2O + O */
    kc[13] = exp((g_RT[8] + g_RT[2]) - (g_RT[7] + g_RT[5]));

    /*reaction 15: HO2 + O <=> OH + O2 */
    kc[14] = exp((g_RT[8] + g_RT[5]) - (g_RT[4] + g_RT[3]));

    /*reaction 16: HO2 + OH <=> H2O + O2 */
    kc[15] = exp((g_RT[8] + g_RT[4]) - (g_RT[7] + g_RT[3]));

    /*reaction 17: 2 OH (+M) <=> H2O2 (+M) */
    kc[16] = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[9]));

    /*reaction 18: 2 HO2 <=> H2O2 + O2 */
    kc[17] = exp((2 * g_RT[8]) - (g_RT[9] + g_RT[3]));

    /*reaction 19: H2O2 + H <=> HO2 + H2 */
    kc[18] = exp((g_RT[9] + g_RT[2]) - (g_RT[8] + g_RT[6]));

    /*reaction 20: H2O2 + H <=> H2O + OH */
    kc[19] = exp((g_RT[9] + g_RT[2]) - (g_RT[7] + g_RT[4]));

    /*reaction 21: H2O2 + OH <=> H2O + HO2 */
    kc[20] = exp((g_RT[9] + g_RT[4]) - (g_RT[7] + g_RT[8]));

    /*reaction 22: H2O2 + O <=> HO2 + OH */
    kc[21] = exp((g_RT[9] + g_RT[5]) - (g_RT[8] + g_RT[4]));

    /*reaction 23: CO + OH <=> CO2 + H */
    kc[22] = exp((g_RT[10] + g_RT[4]) - (g_RT[11] + g_RT[2]));

    /*reaction 24: CO + HO2 <=> CO2 + OH */
    kc[23] = exp((g_RT[10] + g_RT[8]) - (g_RT[11] + g_RT[4]));

    /*reaction 25: HCO + M <=> CO + H + M */
    kc[24] = refC * exp((g_RT[12]) - (g_RT[10] + g_RT[2]));

    /*reaction 26: HCO + H <=> CO + H2 */
    kc[25] = exp((g_RT[12] + g_RT[2]) - (g_RT[10] + g_RT[6]));

    /*reaction 27: HCO + O <=> CO + OH */
    kc[26] = exp((g_RT[12] + g_RT[5]) - (g_RT[10] + g_RT[4]));

    /*reaction 28: HCO + O <=> CO2 + H */
    kc[27] = exp((g_RT[12] + g_RT[5]) - (g_RT[11] + g_RT[2]));

    /*reaction 29: HCO + OH <=> CO + H2O */
    kc[28] = exp((g_RT[12] + g_RT[4]) - (g_RT[10] + g_RT[7]));

    /*reaction 30: HCO + O2 <=> CO + HO2 */
    kc[29] = exp((g_RT[12] + g_RT[3]) - (g_RT[10] + g_RT[8]));

    /*reaction 31: CH2O + M <=> HCO + H + M */
    kc[30] = refC * exp((g_RT[13]) - (g_RT[12] + g_RT[2]));

    /*reaction 32: CH2O + H <=> HCO + H2 */
    kc[31] = exp((g_RT[13] + g_RT[2]) - (g_RT[12] + g_RT[6]));

    /*reaction 33: CH2O + O <=> HCO + OH */
    kc[32] = exp((g_RT[13] + g_RT[5]) - (g_RT[12] + g_RT[4]));

    /*reaction 34: CH2O + OH <=> HCO + H2O */
    kc[33] = exp((g_RT[13] + g_RT[4]) - (g_RT[12] + g_RT[7]));

    /*reaction 35: CH4 + H <=> H2 + CH3 */
    kc[34] = exp((g_RT[14] + g_RT[2]) - (g_RT[6] + g_RT[15]));

    /*reaction 36: CH4 + OH <=> H2O + CH3 */
    kc[35] = exp((g_RT[14] + g_RT[4]) - (g_RT[7] + g_RT[15]));

    /*reaction 37: CH4 + O <=> CH3 + OH */
    kc[36] = exp((g_RT[14] + g_RT[5]) - (g_RT[15] + g_RT[4]));

    /*reaction 38: CH4 + O2 <=> CH3 + HO2 */
    kc[37] = exp((g_RT[14] + g_RT[3]) - (g_RT[15] + g_RT[8]));

    /*reaction 39: CH4 + HO2 <=> CH3 + H2O2 */
    kc[38] = exp((g_RT[14] + g_RT[8]) - (g_RT[15] + g_RT[9]));

    /*reaction 40: CH3 + H <=> T-CH2 + H2 */
    kc[39] = exp((g_RT[15] + g_RT[2]) - (g_RT[16] + g_RT[6]));

    /*reaction 41: CH3 + H <=> S-CH2 + H2 */
    kc[40] = exp((g_RT[15] + g_RT[2]) - (g_RT[17] + g_RT[6]));

    /*reaction 42: CH3 + OH <=> S-CH2 + H2O */
    kc[41] = exp((g_RT[15] + g_RT[4]) - (g_RT[17] + g_RT[7]));

    /*reaction 43: CH3 + O <=> CH2O + H */
    kc[42] = exp((g_RT[15] + g_RT[5]) - (g_RT[13] + g_RT[2]));

    /*reaction 44: CH3 + T-CH2 <=> C2H4 + H */
    kc[43] = exp((g_RT[15] + g_RT[16]) - (g_RT[18] + g_RT[2]));

    /*reaction 45: CH3 + HO2 <=> CH3O + OH */
    kc[44] = exp((g_RT[15] + g_RT[8]) - (g_RT[19] + g_RT[4]));

    /*reaction 46: CH3 + O2 <=> CH2O + OH */
    kc[45] = exp((g_RT[15] + g_RT[3]) - (g_RT[13] + g_RT[4]));

    /*reaction 47: CH3 + O2 <=> CH3O + O */
    kc[46] = exp((g_RT[15] + g_RT[3]) - (g_RT[19] + g_RT[5]));

    /*reaction 48: 2 CH3 <=> C2H4 + H2 */
    kc[47] = exp((2 * g_RT[15]) - (g_RT[18] + g_RT[6]));

    /*reaction 49: 2 CH3 <=> C2H5 + H */
    kc[48] = exp((2 * g_RT[15]) - (g_RT[20] + g_RT[2]));

    /*reaction 50: CH3 + H (+M) <=> CH4 (+M) */
    kc[49] = 1.0 / (refC) * exp((g_RT[15] + g_RT[2]) - (g_RT[14]));

    /*reaction 51: 2 CH3 (+M) <=> C2H6 (+M) */
    kc[50] = 1.0 / (refC) * exp((2 * g_RT[15]) - (g_RT[21]));

    /*reaction 52: S-CH2 + OH <=> CH2O + H */
    kc[51] = exp((g_RT[17] + g_RT[4]) - (g_RT[13] + g_RT[2]));

    /*reaction 53: S-CH2 + O2 <=> CO + OH + H */
    kc[52] = refC * exp((g_RT[17] + g_RT[3]) - (g_RT[10] + g_RT[4] + g_RT[2]));

    /*reaction 54: S-CH2 + CO2 <=> CO + CH2O */
    kc[53] = exp((g_RT[17] + g_RT[11]) - (g_RT[10] + g_RT[13]));

    /*reaction 55: S-CH2 + M <=> T-CH2 + M */
    kc[54] = exp((g_RT[17]) - (g_RT[16]));

    /*reaction 56: T-CH2 + H <=> CH + H2 */
    kc[55] = exp((g_RT[16] + g_RT[2]) - (g_RT[22] + g_RT[6]));

    /*reaction 57: T-CH2 + OH <=> CH2O + H */
    kc[56] = exp((g_RT[16] + g_RT[4]) - (g_RT[13] + g_RT[2]));

    /*reaction 58: T-CH2 + OH <=> CH + H2O */
    kc[57] = exp((g_RT[16] + g_RT[4]) - (g_RT[22] + g_RT[7]));

    /*reaction 59: T-CH2 + O <=> CO + 2 H */
    kc[58] = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[10] + 2 * g_RT[2]));

    /*reaction 60: T-CH2 + O <=> CO + H2 */
    kc[59] = exp((g_RT[16] + g_RT[5]) - (g_RT[10] + g_RT[6]));

    /*reaction 61: T-CH2 + O2 <=> CO2 + H2 */
    kc[60] = exp((g_RT[16] + g_RT[3]) - (g_RT[11] + g_RT[6]));

    /*reaction 62: T-CH2 + O2 <=> CO + OH + H */
    kc[61] = refC * exp((g_RT[16] + g_RT[3]) - (g_RT[10] + g_RT[4] + g_RT[2]));

    /*reaction 63: 2 T-CH2 <=> C2H2 + 2 H */
    kc[62] = refC * exp((2 * g_RT[16]) - (g_RT[23] + 2 * g_RT[2]));

    /*reaction 64: CH + O <=> CO + H */
    kc[63] = exp((g_RT[22] + g_RT[5]) - (g_RT[10] + g_RT[2]));

    /*reaction 65: CH + O2 <=> HCO + O */
    kc[64] = exp((g_RT[22] + g_RT[3]) - (g_RT[12] + g_RT[5]));

    /*reaction 66: CH + H2O <=> CH2O + H */
    kc[65] = exp((g_RT[22] + g_RT[7]) - (g_RT[13] + g_RT[2]));

    /*reaction 67: CH + CO2 <=> HCO + CO */
    kc[66] = exp((g_RT[22] + g_RT[11]) - (g_RT[12] + g_RT[10]));

    /*reaction 68: CH3O + H <=> CH2O + H2 */
    kc[67] = exp((g_RT[19] + g_RT[2]) - (g_RT[13] + g_RT[6]));

    /*reaction 69: CH3O + H <=> S-CH2 + H2O */
    kc[68] = exp((g_RT[19] + g_RT[2]) - (g_RT[17] + g_RT[7]));

    /*reaction 70: CH3O + OH <=> CH2O + H2O */
    kc[69] = exp((g_RT[19] + g_RT[4]) - (g_RT[13] + g_RT[7]));

    /*reaction 71: CH3O + O <=> OH + CH2O */
    kc[70] = exp((g_RT[19] + g_RT[5]) - (g_RT[4] + g_RT[13]));

    /*reaction 72: CH3O + O2 <=> CH2O + HO2 */
    kc[71] = exp((g_RT[19] + g_RT[3]) - (g_RT[13] + g_RT[8]));

    /*reaction 73: CH3O + M <=> CH2O + H + M */
    kc[72] = refC * exp((g_RT[19]) - (g_RT[13] + g_RT[2]));

    /*reaction 74: C2H6 + H <=> C2H5 + H2 */
    kc[73] = exp((g_RT[21] + g_RT[2]) - (g_RT[20] + g_RT[6]));

    /*reaction 75: C2H6 + O <=> C2H5 + OH */
    kc[74] = exp((g_RT[21] + g_RT[5]) - (g_RT[20] + g_RT[4]));

    /*reaction 76: C2H6 + OH <=> C2H5 + H2O */
    kc[75] = exp((g_RT[21] + g_RT[4]) - (g_RT[20] + g_RT[7]));

    /*reaction 77: C2H6 + CH3 <=> C2H5 + CH4 */
    kc[76] = exp((g_RT[21] + g_RT[15]) - (g_RT[20] + g_RT[14]));

    /*reaction 78: C2H6 (+M) <=> C2H5 + H (+M) */
    kc[77] = refC * exp((g_RT[21]) - (g_RT[20] + g_RT[2]));

    /*reaction 79: C2H5 + H <=> C2H4 + H2 */
    kc[78] = exp((g_RT[20] + g_RT[2]) - (g_RT[18] + g_RT[6]));

    /*reaction 80: C2H5 + O <=> C2H4 + OH */
    kc[79] = exp((g_RT[20] + g_RT[5]) - (g_RT[18] + g_RT[4]));

    /*reaction 81: C2H5 + O <=> CH3 + CH2O */
    kc[80] = exp((g_RT[20] + g_RT[5]) - (g_RT[15] + g_RT[13]));

    /*reaction 82: C2H5 + O2 <=> C2H4 + HO2 */
    kc[81] = exp((g_RT[20] + g_RT[3]) - (g_RT[18] + g_RT[8]));

    /*reaction 83: C2H5 (+M) <=> C2H4 + H (+M) */
    kc[82] = refC * exp((g_RT[20]) - (g_RT[18] + g_RT[2]));

    /*reaction 84: C2H4 + H <=> C2H3 + H2 */
    kc[83] = exp((g_RT[18] + g_RT[2]) - (g_RT[24] + g_RT[6]));

    /*reaction 85: C2H4 + OH <=> C2H3 + H2O */
    kc[84] = exp((g_RT[18] + g_RT[4]) - (g_RT[24] + g_RT[7]));

    /*reaction 86: C2H4 + O <=> CH3 + HCO */
    kc[85] = exp((g_RT[18] + g_RT[5]) - (g_RT[15] + g_RT[12]));

    /*reaction 87: C2H4 + O <=> CH2CHO + H */
    kc[86] = exp((g_RT[18] + g_RT[5]) - (g_RT[25] + g_RT[2]));

    /*reaction 88: 2 C2H4 <=> C2H3 + C2H5 */
    kc[87] = exp((2 * g_RT[18]) - (g_RT[24] + g_RT[20]));

    /*reaction 89: C2H4 + O2 <=> C2H3 + HO2 */
    kc[88] = exp((g_RT[18] + g_RT[3]) - (g_RT[24] + g_RT[8]));

    /*reaction 90: C2H4 + HO2 <=> C2H4O + OH */
    kc[89] = exp((g_RT[18] + g_RT[8]) - (g_RT[26] + g_RT[4]));

    /*reaction 91: C2H4O + HO2 <=> CH3 + CO + H2O2 */
    kc[90] = refC * exp((g_RT[26] + g_RT[8]) - (g_RT[15] + g_RT[10] + g_RT[9]));

    /*reaction 92: C2H4 + M <=> C2H3 + H + M */
    kc[91] = refC * exp((g_RT[18]) - (g_RT[24] + g_RT[2]));

    /*reaction 93: C2H4 + M <=> C2H2 + H2 + M */
    kc[92] = refC * exp((g_RT[18]) - (g_RT[23] + g_RT[6]));

    /*reaction 94: C2H3 + H <=> C2H2 + H2 */
    kc[93] = exp((g_RT[24] + g_RT[2]) - (g_RT[23] + g_RT[6]));

    /*reaction 95: C2H3 (+M) <=> C2H2 + H (+M) */
    kc[94] = refC * exp((g_RT[24]) - (g_RT[23] + g_RT[2]));

    /*reaction 96: C2H3 + O2 <=> CH2O + HCO */
    kc[95] = exp((g_RT[24] + g_RT[3]) - (g_RT[13] + g_RT[12]));

    /*reaction 97: C2H3 + O2 <=> CH2CHO + O */
    kc[96] = exp((g_RT[24] + g_RT[3]) - (g_RT[25] + g_RT[5]));

    /*reaction 98: C2H3 + O2 <=> C2H2 + HO2 */
    kc[97] = exp((g_RT[24] + g_RT[3]) - (g_RT[23] + g_RT[8]));

    /*reaction 99: CH2CHO <=> CH2CO + H */
    kc[98] = refC * exp((g_RT[25]) - (g_RT[27] + g_RT[2]));

    /*reaction 100: C2H2 + O <=> HCCO + H */
    kc[99] = exp((g_RT[23] + g_RT[5]) - (g_RT[28] + g_RT[2]));

    /*reaction 101: C2H2 + O <=> T-CH2 + CO */
    kc[100] = exp((g_RT[23] + g_RT[5]) - (g_RT[16] + g_RT[10]));

    /*reaction 102: C2H2 + O2 <=> CH2O + CO */
    kc[101] = exp((g_RT[23] + g_RT[3]) - (g_RT[13] + g_RT[10]));

    /*reaction 103: C2H2 + OH <=> CH2CO + H */
    kc[102] = exp((g_RT[23] + g_RT[4]) - (g_RT[27] + g_RT[2]));

    /*reaction 104: C2H2 + OH <=> C2H + H2O */
    kc[103] = exp((g_RT[23] + g_RT[4]) - (g_RT[29] + g_RT[7]));

    /*reaction 105: CH2CO + H <=> CH3 + CO */
    kc[104] = exp((g_RT[27] + g_RT[2]) - (g_RT[15] + g_RT[10]));

    /*reaction 106: CH2CO + O <=> T-CH2 + CO2 */
    kc[105] = exp((g_RT[27] + g_RT[5]) - (g_RT[16] + g_RT[11]));

    /*reaction 107: CH2CO + O <=> HCCO + OH */
    kc[106] = exp((g_RT[27] + g_RT[5]) - (g_RT[28] + g_RT[4]));

    /*reaction 108: CH2CO + CH3 <=> C2H5 + CO */
    kc[107] = exp((g_RT[27] + g_RT[15]) - (g_RT[20] + g_RT[10]));

    /*reaction 109: HCCO + H <=> S-CH2 + CO */
    kc[108] = exp((g_RT[28] + g_RT[2]) - (g_RT[17] + g_RT[10]));

    /*reaction 110: HCCO + OH <=> HCO + CO + H */
    kc[109] = refC * exp((g_RT[28] + g_RT[4]) - (g_RT[12] + g_RT[10] + g_RT[2]));

    /*reaction 111: HCCO + O <=> 2 CO + H */
    kc[110] = refC * exp((g_RT[28] + g_RT[5]) - (2 * g_RT[10] + g_RT[2]));

    /*reaction 112: HCCO + O2 <=> 2 CO + OH */
    kc[111] = refC * exp((g_RT[28] + g_RT[3]) - (2 * g_RT[10] + g_RT[4]));

    /*reaction 113: HCCO + O2 <=> CO2 + CO + H */
    kc[112] = refC * exp((g_RT[28] + g_RT[3]) - (g_RT[11] + g_RT[10] + g_RT[2]));

    /*reaction 114: C2H + OH <=> HCCO + H */
    kc[113] = exp((g_RT[29] + g_RT[4]) - (g_RT[28] + g_RT[2]));

    /*reaction 115: C2H + O <=> CO + CH */
    kc[114] = exp((g_RT[29] + g_RT[5]) - (g_RT[10] + g_RT[22]));

    /*reaction 116: C2H + O2 <=> HCCO + O */
    kc[115] = exp((g_RT[29] + g_RT[3]) - (g_RT[28] + g_RT[5]));

    /*reaction 117: C2H + O2 <=> CH + CO2 */
    kc[116] = exp((g_RT[29] + g_RT[3]) - (g_RT[22] + g_RT[11]));

    /*reaction 118: C2H + O2 <=> HCO + CO */
    kc[117] = exp((g_RT[29] + g_RT[3]) - (g_RT[12] + g_RT[10]));

    /*reaction 119: CH2OH + H <=> CH2O + H2 */
    kc[118] = exp((g_RT[30] + g_RT[2]) - (g_RT[13] + g_RT[6]));

    /*reaction 120: CH2OH + H <=> CH3 + OH */
    kc[119] = exp((g_RT[30] + g_RT[2]) - (g_RT[15] + g_RT[4]));

    /*reaction 121: CH2OH + OH <=> CH2O + H2O */
    kc[120] = exp((g_RT[30] + g_RT[4]) - (g_RT[13] + g_RT[7]));

    /*reaction 122: CH2OH + O2 <=> CH2O + HO2 */
    kc[121] = exp((g_RT[30] + g_RT[3]) - (g_RT[13] + g_RT[8]));

    /*reaction 123: CH2OH + M <=> CH2O + H + M */
    kc[122] = refC * exp((g_RT[30]) - (g_RT[13] + g_RT[2]));

    /*reaction 124: CH3O + M <=> CH2OH + M */
    kc[123] = exp((g_RT[19]) - (g_RT[30]));

    /*reaction 125: CH2CO + OH <=> CH2OH + CO */
    kc[124] = exp((g_RT[27] + g_RT[4]) - (g_RT[30] + g_RT[10]));

    /*reaction 126: CH3OH + OH <=> CH2OH + H2O */
    kc[125] = exp((g_RT[31] + g_RT[4]) - (g_RT[30] + g_RT[7]));

    /*reaction 127: CH3OH + OH <=> CH3O + H2O */
    kc[126] = exp((g_RT[31] + g_RT[4]) - (g_RT[19] + g_RT[7]));

    /*reaction 128: CH3OH + H <=> CH2OH + H2 */
    kc[127] = exp((g_RT[31] + g_RT[2]) - (g_RT[30] + g_RT[6]));

    /*reaction 129: CH3OH + H <=> CH3O + H2 */
    kc[128] = exp((g_RT[31] + g_RT[2]) - (g_RT[19] + g_RT[6]));

    /*reaction 130: CH3OH + O <=> CH2OH + OH */
    kc[129] = exp((g_RT[31] + g_RT[5]) - (g_RT[30] + g_RT[4]));

    /*reaction 131: CH3OH + HO2 <=> CH2OH + H2O2 */
    kc[130] = exp((g_RT[31] + g_RT[8]) - (g_RT[30] + g_RT[9]));

    /*reaction 132: CH3OH + O2 <=> CH2OH + HO2 */
    kc[131] = exp((g_RT[31] + g_RT[3]) - (g_RT[30] + g_RT[8]));

    /*reaction 133: C3H4 + O <=> C2H4 + CO */
    kc[132] = exp((g_RT[32] + g_RT[5]) - (g_RT[18] + g_RT[10]));

    /*reaction 134: CH3 + C2H2 <=> C3H4 + H */
    kc[133] = exp((g_RT[15] + g_RT[23]) - (g_RT[32] + g_RT[2]));

    /*reaction 135: C3H4 + O <=> HCCO + CH3 */
    kc[134] = exp((g_RT[32] + g_RT[5]) - (g_RT[28] + g_RT[15]));

    /*reaction 136: C3H3 + H (+M) <=> C3H4 (+M) */
    kc[135] = 1.0 / (refC) * exp((g_RT[33] + g_RT[2]) - (g_RT[32]));

    /*reaction 137: C3H3 + HO2 <=> C3H4 + O2 */
    kc[136] = exp((g_RT[33] + g_RT[8]) - (g_RT[32] + g_RT[3]));

    /*reaction 138: C3H4 + OH <=> C3H3 + H2O */
    kc[137] = exp((g_RT[32] + g_RT[4]) - (g_RT[33] + g_RT[7]));

    /*reaction 139: C3H3 + O2 <=> CH2CO + HCO */
    kc[138] = exp((g_RT[33] + g_RT[3]) - (g_RT[27] + g_RT[12]));

    /*reaction 140: C3H4 + H (+M) <=> C3H5 (+M) */
    kc[139] = 1.0 / (refC) * exp((g_RT[32] + g_RT[2]) - (g_RT[34]));

    /*reaction 141: C3H5 + H <=> C3H4 + H2 */
    kc[140] = exp((g_RT[34] + g_RT[2]) - (g_RT[32] + g_RT[6]));

    /*reaction 142: C3H5 + O2 <=> C3H4 + HO2 */
    kc[141] = exp((g_RT[34] + g_RT[3]) - (g_RT[32] + g_RT[8]));

    /*reaction 143: C3H5 + CH3 <=> C3H4 + CH4 */
    kc[142] = exp((g_RT[34] + g_RT[15]) - (g_RT[32] + g_RT[14]));

    /*reaction 144: C2H2 + CH3 (+M) <=> C3H5 (+M) */
    kc[143] = 1.0 / (refC) * exp((g_RT[23] + g_RT[15]) - (g_RT[34]));

    /*reaction 145: C3H5 + OH <=> C3H4 + H2O */
    kc[144] = exp((g_RT[34] + g_RT[4]) - (g_RT[32] + g_RT[7]));

    /*reaction 146: C3H3 + HCO <=> C3H4 + CO */
    kc[145] = exp((g_RT[33] + g_RT[12]) - (g_RT[32] + g_RT[10]));

    /*reaction 147: C3H3 + HO2 <=> OH + CO + C2H3 */
    kc[146] = refC * exp((g_RT[33] + g_RT[8]) - (g_RT[4] + g_RT[10] + g_RT[24]));

    /*reaction 148: C3H4 + O2 <=> CH3 + HCO + CO */
    kc[147] = refC * exp((g_RT[32] + g_RT[3]) - (g_RT[15] + g_RT[12] + g_RT[10]));

    /*reaction 149: C3H6 + O <=> C2H5 + HCO */
    kc[148] = exp((g_RT[35] + g_RT[5]) - (g_RT[20] + g_RT[12]));

    /*reaction 150: C3H6 + OH <=> C3H5 + H2O */
    kc[149] = exp((g_RT[35] + g_RT[4]) - (g_RT[34] + g_RT[7]));

    /*reaction 151: C3H6 + O <=> CH2CO + CH3 + H */
    kc[150] = refC * exp((g_RT[35] + g_RT[5]) - (g_RT[27] + g_RT[15] + g_RT[2]));

    /*reaction 152: C3H6 + H <=> C3H5 + H2 */
    kc[151] = exp((g_RT[35] + g_RT[2]) - (g_RT[34] + g_RT[6]));

    /*reaction 153: C3H5 + H (+M) <=> C3H6 (+M) */
    kc[152] = 1.0 / (refC) * exp((g_RT[34] + g_RT[2]) - (g_RT[35]));

    /*reaction 154: C3H5 + HO2 <=> C3H6 + O2 */
    kc[153] = exp((g_RT[34] + g_RT[8]) - (g_RT[35] + g_RT[3]));

    /*reaction 155: C3H5 + HO2 <=> OH + C2H3 + CH2O */
    kc[154] = refC * exp((g_RT[34] + g_RT[8]) - (g_RT[4] + g_RT[24] + g_RT[13]));

    /*reaction 156: C2H3 + CH3 (+M) <=> C3H6 (+M) */
    kc[155] = 1.0 / (refC) * exp((g_RT[24] + g_RT[15]) - (g_RT[35]));

    /*reaction 157: C3H6 + H <=> C2H4 + CH3 */
    kc[156] = exp((g_RT[35] + g_RT[2]) - (g_RT[18] + g_RT[15]));

    /*reaction 158: CH3 + C2H3 <=> C3H5 + H */
    kc[157] = exp((g_RT[15] + g_RT[24]) - (g_RT[34] + g_RT[2]));

    /*reaction 159: C3H8 (+M) <=> CH3 + C2H5 (+M) */
    kc[158] = refC * exp((g_RT[36]) - (g_RT[15] + g_RT[20]));

    /*reaction 160: C3H8 + O2 <=> I-C3H7 + HO2 */
    kc[159] = exp((g_RT[36] + g_RT[3]) - (g_RT[37] + g_RT[8]));

    /*reaction 161: C3H8 + O2 <=> N-C3H7 + HO2 */
    kc[160] = exp((g_RT[36] + g_RT[3]) - (g_RT[38] + g_RT[8]));

    /*reaction 162: C3H8 + H <=> I-C3H7 + H2 */
    kc[161] = exp((g_RT[36] + g_RT[2]) - (g_RT[37] + g_RT[6]));

    /*reaction 163: C3H8 + H <=> N-C3H7 + H2 */
    kc[162] = exp((g_RT[36] + g_RT[2]) - (g_RT[38] + g_RT[6]));

    /*reaction 164: C3H8 + O <=> I-C3H7 + OH */
    kc[163] = exp((g_RT[36] + g_RT[5]) - (g_RT[37] + g_RT[4]));

    /*reaction 165: C3H8 + O <=> N-C3H7 + OH */
    kc[164] = exp((g_RT[36] + g_RT[5]) - (g_RT[38] + g_RT[4]));

    /*reaction 166: C3H8 + OH <=> N-C3H7 + H2O */
    kc[165] = exp((g_RT[36] + g_RT[4]) - (g_RT[38] + g_RT[7]));

    /*reaction 167: C3H8 + OH <=> I-C3H7 + H2O */
    kc[166] = exp((g_RT[36] + g_RT[4]) - (g_RT[37] + g_RT[7]));

    /*reaction 168: C3H8 + HO2 <=> I-C3H7 + H2O2 */
    kc[167] = exp((g_RT[36] + g_RT[8]) - (g_RT[37] + g_RT[9]));

    /*reaction 169: C3H8 + HO2 <=> N-C3H7 + H2O2 */
    kc[168] = exp((g_RT[36] + g_RT[8]) - (g_RT[38] + g_RT[9]));

    /*reaction 170: I-C3H7 + C3H8 <=> N-C3H7 + C3H8 */
    kc[169] = exp((g_RT[37] + g_RT[36]) - (g_RT[38] + g_RT[36]));

    /*reaction 171: C3H6 + H (+M) <=> I-C3H7 (+M) */
    kc[170] = 1.0 / (refC) * exp((g_RT[35] + g_RT[2]) - (g_RT[37]));

    /*reaction 172: I-C3H7 + O2 <=> C3H6 + HO2 */
    kc[171] = exp((g_RT[37] + g_RT[3]) - (g_RT[35] + g_RT[8]));

    /*reaction 173: N-C3H7 (+M) <=> CH3 + C2H4 (+M) */
    kc[172] = refC * exp((g_RT[38]) - (g_RT[15] + g_RT[18]));

    /*reaction 174: H + C3H6 (+M) <=> N-C3H7 (+M) */
    kc[173] = 1.0 / (refC) * exp((g_RT[2] + g_RT[35]) - (g_RT[38]));

    /*reaction 175: N-C3H7 + O2 <=> C3H6 + HO2 */
    kc[174] = exp((g_RT[38] + g_RT[3]) - (g_RT[35] + g_RT[8]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: N2 */
        species[0] =
            -1.02089990e+03 / tc[1]
            -6.51695000e-01
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 1: AR */
        species[1] =
            -7.45375000e+02 / tc[1]
            -1.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 2: H */
        species[2] =
            +2.54736599e+04 / tc[1]
            +2.94668285e+00
            -2.50000000e+00 * tc[0]
            -3.52666409e-13 * tc[1]
            +3.32653273e-16 * tc[2]
            -1.91734693e-19 * tc[3]
            +4.63866166e-23 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.06394356e+03 / tc[1]
            +1.24780630e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.38153812e+03 / tc[1]
            +4.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
        /*species 5: O */
        species[5] =
            +2.91222592e+04 / tc[1]
            +1.11633364e+00
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 6: H2 */
        species[6] =
            -9.17935173e+02 / tc[1]
            +1.66132088e+00
            -2.34433112e+00 * tc[0]
            -3.99026037e-03 * tc[1]
            +3.24635850e-06 * tc[2]
            -1.67976745e-09 * tc[3]
            +3.68805881e-13 * tc[4];
        /*species 7: H2O */
        species[7] =
            -3.02937267e+04 / tc[1]
            +5.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +2.94808040e+02 / tc[1]
            +5.85135560e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 9: H2O2 */
        species[9] =
            -1.77025821e+04 / tc[1]
            +8.41061950e-01
            -4.27611269e+00 * tc[0]
            +2.71411208e-04 * tc[1]
            -2.78892835e-06 * tc[2]
            +1.79809011e-09 * tc[3]
            -4.31227182e-13 * tc[4];
        /*species 10: CO */
        species[10] =
            -1.43440860e+04 / tc[1]
            +7.11241900e-02
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 11: CO2 */
        species[11] =
            -4.83719697e+04 / tc[1]
            -7.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 12: HCO */
        species[12] =
            +3.83956496e+03 / tc[1]
            +8.26813410e-01
            -4.22118584e+00 * tc[0]
            +1.62196266e-03 * tc[1]
            -2.29665743e-06 * tc[2]
            +1.10953411e-09 * tc[3]
            -2.16884432e-13 * tc[4];
        /*species 13: CH2O */
        species[13] =
            -1.43089567e+04 / tc[1]
            +4.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 14: CH4 */
        species[14] =
            -1.02466476e+04 / tc[1]
            +9.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +1.64449988e+04 / tc[1]
            +2.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 16: T-CH2 */
        species[16] =
            +4.60040401e+04 / tc[1]
            +2.20014682e+00
            -3.76267867e+00 * tc[0]
            -4.84436072e-04 * tc[1]
            -4.65816402e-07 * tc[2]
            +3.20909294e-10 * tc[3]
            -8.43708595e-14 * tc[4];
        /*species 17: S-CH2 */
        species[17] =
            +5.04968163e+04 / tc[1]
            +4.96772308e+00
            -4.19860411e+00 * tc[0]
            +1.18330710e-03 * tc[1]
            -1.37216037e-06 * tc[2]
            +5.57346651e-10 * tc[3]
            -9.71573685e-14 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +5.08977593e+03 / tc[1]
            -1.38129480e-01
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +1.30772484e+03 / tc[1]
            -2.86060362e+00
            -3.71180502e+00 * tc[0]
            +1.40231653e-03 * tc[1]
            -6.27584952e-06 * tc[2]
            +3.94226741e-09 * tc[3]
            -9.32942100e-13 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +1.28416265e+04 / tc[1]
            -4.00743560e-01
            -4.30646568e+00 * tc[0]
            +2.09329446e-03 * tc[1]
            -8.28571345e-06 * tc[2]
            +4.99272172e-09 * tc[3]
            -1.15254502e-12 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            -1.15222055e+04 / tc[1]
            +1.62460176e+00
            -4.29142492e+00 * tc[0]
            +2.75077135e-03 * tc[1]
            -9.99063813e-06 * tc[2]
            +5.90388571e-09 * tc[3]
            -1.34342886e-12 * tc[4];
        /*species 22: CH */
        species[22] =
            +7.07972934e+04 / tc[1]
            +1.40580557e+00
            -3.48981665e+00 * tc[0]
            -1.61917771e-04 * tc[1]
            +2.81498442e-07 * tc[2]
            -2.63514439e-10 * tc[3]
            +7.03045335e-14 * tc[4];
        /*species 23: C2H2 */
        species[23] =
            +2.64289807e+04 / tc[1]
            -1.31310240e+01
            -8.08681094e-01 * tc[0]
            -1.16807815e-02 * tc[1]
            +5.91953025e-06 * tc[2]
            -2.33460364e-09 * tc[3]
            +4.25036487e-13 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +3.48598468e+04 / tc[1]
            -5.29807380e+00
            -3.21246645e+00 * tc[0]
            -7.57395810e-04 * tc[1]
            -4.32015687e-06 * tc[2]
            +2.98048206e-09 * tc[3]
            -7.35754365e-13 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +3.80428530e+02 / tc[1]
            -1.83431519e+01
            -1.01340010e+00 * tc[0]
            -1.13407335e-02 * tc[1]
            +2.62232400e-06 * tc[2]
            -3.37429192e-10 * tc[3]
            -1.47995060e-14 * tc[4];
        /*species 26: C2H4O */
        species[26] =
            -2.15728780e+04 / tc[1]
            +6.26443600e-01
            -4.72945950e+00 * tc[0]
            +1.59664290e-03 * tc[1]
            -7.92248683e-06 * tc[2]
            +4.78821758e-09 * tc[3]
            -1.09655560e-12 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            -7.04291804e+03 / tc[1]
            -1.00798117e+01
            -2.13583630e+00 * tc[0]
            -9.05943605e-03 * tc[1]
            +2.89912457e-06 * tc[2]
            -7.78664640e-10 * tc[3]
            +1.00728807e-13 * tc[4];
        /*species 28: HCCO */
        species[28] =
            +2.00594490e+04 / tc[1]
            -1.02386956e+01
            -2.25172140e+00 * tc[0]
            -8.82751050e-03 * tc[1]
            +3.95485017e-06 * tc[2]
            -1.43964658e-09 * tc[3]
            +2.53324055e-13 * tc[4];
        /*species 29: C2H */
        species[29] =
            +6.68393932e+04 / tc[1]
            -3.33330705e+00
            -2.88965733e+00 * tc[0]
            -6.70498055e-03 * tc[1]
            +4.74615835e-06 * tc[2]
            -2.45659204e-09 * tc[3]
            +5.46657555e-13 * tc[4];
        /*species 30: CH2OH */
        species[30] =
            -3.52476728e+03 / tc[1]
            +1.16920333e+00
            -4.47832317e+00 * tc[0]
            +6.75348435e-04 * tc[1]
            -4.64139512e-06 * tc[2]
            +3.04056164e-09 * tc[3]
            -7.39533875e-13 * tc[4];
        /*species 31: CH3OH */
        species[31] =
            -2.56427656e+04 / tc[1]
            +7.21949405e+00
            -5.71539582e+00 * tc[0]
            +7.61545645e-03 * tc[1]
            -1.08740193e-05 * tc[2]
            +5.92339074e-09 * tc[3]
            -1.30676349e-12 * tc[4];
        /*species 32: C3H4 */
        species[32] =
            +2.15415670e+04 / tc[1]
            -7.61309450e+00
            -2.61304450e+00 * tc[0]
            -6.06128750e-03 * tc[1]
            -3.08998000e-06 * tc[2]
            +2.87709575e-09 * tc[3]
            -7.66753950e-13 * tc[4];
        /*species 33: C3H3 */
        species[33] =
            +4.01057783e+04 / tc[1]
            -1.38547831e+01
            -1.35110927e+00 * tc[0]
            -1.63705612e-02 * tc[1]
            +7.89711892e-06 * tc[2]
            -3.13591507e-09 * tc[3]
            +5.92704615e-13 * tc[4];
        /*species 34: C3H5 */
        species[34] =
            +1.92456290e+04 / tc[1]
            -1.58100305e+01
            -1.36318350e+00 * tc[0]
            -9.90691050e-03 * tc[1]
            -2.08284333e-06 * tc[2]
            +2.77962958e-09 * tc[3]
            -7.92328550e-13 * tc[4];
        /*species 35: C3H6 */
        species[35] =
            +1.07482600e+03 / tc[1]
            -1.46520330e+01
            -1.49330700e+00 * tc[0]
            -1.04625900e-02 * tc[1]
            -7.47799000e-07 * tc[2]
            +1.39076000e-09 * tc[3]
            -3.57907300e-13 * tc[4];
        /*species 36: C3H8 */
        species[36] =
            -1.40579070e+04 / tc[1]
            -1.82970271e+01
            -9.28510930e-01 * tc[0]
            -1.32302830e-02 * tc[1]
            -1.00554077e-06 * tc[2]
            +1.82624608e-09 * tc[3]
            -4.74807720e-13 * tc[4];
        /*species 37: I-C3H7 */
        species[37] =
            +9.42237240e+03 / tc[1]
            -1.86713971e+01
            -1.44491990e+00 * tc[0]
            -1.04995560e-02 * tc[1]
            -1.28393703e-06 * tc[2]
            +1.53968775e-09 * tc[3]
            -3.56414810e-13 * tc[4];
        /*species 38: N-C3H7 */
        species[38] =
            +1.03123460e+04 / tc[1]
            -2.00869167e+01
            -1.04911730e+00 * tc[0]
            -1.30044865e-02 * tc[1]
            -3.92375267e-07 * tc[2]
            +1.63292767e-09 * tc[3]
            -4.68601035e-13 * tc[4];
    } else {
        /*species 0: N2 */
        species[0] =
            -9.22797700e+02 / tc[1]
            -3.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 1: AR */
        species[1] =
            -7.45375000e+02 / tc[1]
            -1.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 2: H */
        species[2] =
            +2.54736599e+04 / tc[1]
            +2.94668292e+00
            -2.50000001e+00 * tc[0]
            +1.15421486e-11 * tc[1]
            -2.69269913e-15 * tc[2]
            +3.94596029e-19 * tc[3]
            -2.49098679e-23 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.08845772e+03 / tc[1]
            -2.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.71885774e+03 / tc[1]
            -2.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
        /*species 5: O */
        species[5] =
            +2.92175791e+04 / tc[1]
            -2.21491786e+00
            -2.56942078e+00 * tc[0]
            +4.29870569e-05 * tc[1]
            -6.99140982e-09 * tc[2]
            +8.34814992e-13 * tc[3]
            -6.14168455e-17 * tc[4];
        /*species 6: H2 */
        species[6] =
            -9.50158922e+02 / tc[1]
            +6.54230251e+00
            -3.33727920e+00 * tc[0]
            +2.47012365e-05 * tc[1]
            -8.32427963e-08 * tc[2]
            +1.49638662e-11 * tc[3]
            -1.00127688e-15 * tc[4];
        /*species 7: H2O */
        species[7] =
            -3.00042971e+04 / tc[1]
            -1.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +1.11856713e+02 / tc[1]
            +2.32108750e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 9: H2O2 */
        species[9] =
            -1.78617877e+04 / tc[1]
            +1.24884623e+00
            -4.16500285e+00 * tc[0]
            -2.45415847e-03 * tc[1]
            +3.16898708e-07 * tc[2]
            -3.09321655e-11 * tc[3]
            +1.43954153e-15 * tc[4];
        /*species 10: CO */
        species[10] =
            -1.41518724e+04 / tc[1]
            -5.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 11: CO2 */
        species[11] =
            -4.87591660e+04 / tc[1]
            +1.58582223e+00
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 12: HCO */
        species[12] =
            +4.01191815e+03 / tc[1]
            -7.02617054e+00
            -2.77217438e+00 * tc[0]
            -2.47847763e-03 * tc[1]
            +4.14076022e-07 * tc[2]
            -4.90968148e-11 * tc[3]
            +2.66754356e-15 * tc[4];
        /*species 13: CH2O */
        species[13] =
            -1.39958323e+04 / tc[1]
            -1.18956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 14: CH4 */
        species[14] =
            -9.46834459e+03 / tc[1]
            -1.83624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +1.67755843e+04 / tc[1]
            -6.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 16: T-CH2 */
        species[16] =
            +4.62636040e+04 / tc[1]
            -3.29709211e+00
            -2.87410113e+00 * tc[0]
            -1.82819646e-03 * tc[1]
            +2.34824328e-07 * tc[2]
            -2.16816291e-11 * tc[3]
            +9.38637835e-16 * tc[4];
        /*species 17: S-CH2 */
        species[17] =
            +5.09259997e+04 / tc[1]
            -6.33446327e+00
            -2.29203842e+00 * tc[0]
            -2.32794318e-03 * tc[1]
            +3.35319912e-07 * tc[2]
            -3.48255000e-11 * tc[3]
            +1.69858182e-15 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +4.93988614e+03 / tc[1]
            -8.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +3.90139164e+02 / tc[1]
            +6.72459266e+00
            -4.75779238e+00 * tc[0]
            -3.72071237e-03 * tc[1]
            +4.49508627e-07 * tc[2]
            -3.65075420e-11 * tc[3]
            +1.31768549e-15 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +1.28575200e+04 / tc[1]
            -1.15077779e+01
            -1.95465642e+00 * tc[0]
            -8.69863610e-03 * tc[1]
            +1.33034445e-06 * tc[2]
            -1.46014741e-10 * tc[3]
            +7.48207880e-15 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            -1.14263932e+04 / tc[1]
            -1.40437292e+01
            -1.07188150e+00 * tc[0]
            -1.08426339e-02 * tc[1]
            +1.67093445e-06 * tc[2]
            -1.84510001e-10 * tc[3]
            +9.50014450e-15 * tc[4];
        /*species 22: CH */
        species[22] =
            +7.10124364e+04 / tc[1]
            -2.60651526e+00
            -2.87846473e+00 * tc[0]
            -4.85456840e-04 * tc[1]
            -2.40742758e-08 * tc[2]
            +1.08906541e-11 * tc[3]
            -8.80396915e-16 * tc[4];
        /*species 23: C2H2 */
        species[23] =
            +2.59359992e+04 / tc[1]
            +5.37785085e+00
            -4.14756964e+00 * tc[0]
            -2.98083332e-03 * tc[1]
            +3.95491420e-07 * tc[2]
            -3.89510143e-11 * tc[3]
            +1.80617607e-15 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +3.46128739e+04 / tc[1]
            -4.77059978e+00
            -3.01672400e+00 * tc[0]
            -5.16511460e-03 * tc[1]
            +7.80137248e-07 * tc[2]
            -8.48027400e-11 * tc[3]
            +4.31303520e-15 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            -7.31993470e+02 / tc[1]
            +7.12953670e+00
            -5.16620060e+00 * tc[0]
            -5.42391300e-03 * tc[1]
            +7.44306133e-07 * tc[2]
            -6.71904567e-11 * tc[3]
            +2.42050965e-15 * tc[4];
        /*species 26: C2H4O */
        species[26] =
            -2.25931220e+04 / tc[1]
            +8.88490250e+00
            -5.40411080e+00 * tc[0]
            -5.86152950e-03 * tc[1]
            +7.04385617e-07 * tc[2]
            -5.69770425e-11 * tc[3]
            +2.04924315e-15 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            -7.55105311e+03 / tc[1]
            +3.87905011e+00
            -4.51129732e+00 * tc[0]
            -4.50179872e-03 * tc[1]
            +6.94899392e-07 * tc[2]
            -7.69454902e-11 * tc[3]
            +3.97419100e-15 * tc[4];
        /*species 28: HCCO */
        species[28] =
            +1.93272150e+04 / tc[1]
            +9.55846530e+00
            -5.62820580e+00 * tc[0]
            -2.04267005e-03 * tc[1]
            +2.65575783e-07 * tc[2]
            -2.38550433e-11 * tc[3]
            +9.70391600e-16 * tc[4];
        /*species 29: C2H */
        species[29] =
            +6.71210650e+04 / tc[1]
            -3.46808823e+00
            -3.16780652e+00 * tc[0]
            -2.37610951e-03 * tc[1]
            +3.06311795e-07 * tc[2]
            -2.53491877e-11 * tc[3]
            +8.86163850e-16 * tc[4];
        /*species 30: CH2OH */
        species[30] =
            -4.05813228e+03 / tc[1]
            +6.94002650e+00
            -5.09312037e+00 * tc[0]
            -2.97379275e-03 * tc[1]
            +3.44160873e-07 * tc[2]
            -2.69172253e-11 * tc[3]
            +9.40625260e-16 * tc[4];
        /*species 31: CH3OH */
        species[31] =
            -2.53748747e+04 / tc[1]
            -1.27126544e+01
            -1.78970791e+00 * tc[0]
            -7.04691460e-03 * tc[1]
            +1.06083472e-06 * tc[2]
            -1.15142571e-10 * tc[3]
            +5.85301100e-15 * tc[4];
        /*species 32: C3H4 */
        species[32] =
            +2.01174950e+04 / tc[1]
            +1.73126382e+01
            -6.31687220e+00 * tc[0]
            -5.56686400e-03 * tc[1]
            +6.60489633e-07 * tc[2]
            -5.29701983e-11 * tc[3]
            +1.89377700e-15 * tc[4];
        /*species 33: C3H3 */
        species[33] =
            +3.89087427e+04 / tc[1]
            +1.97270624e+01
            -7.14221880e+00 * tc[0]
            -3.80951002e-03 * tc[1]
            +4.45766583e-07 * tc[2]
            -3.54095668e-11 * tc[3]
            +1.25737708e-15 * tc[4];
        /*species 34: C3H5 */
        species[34] =
            +1.74824490e+04 / tc[1]
            +1.77438377e+01
            -6.50078770e+00 * tc[0]
            -7.16236550e-03 * tc[1]
            +9.46360533e-07 * tc[2]
            -9.23400083e-11 * tc[3]
            +4.51819435e-15 * tc[4];
        /*species 35: C3H6 */
        species[35] =
            -9.23570300e+02 / tc[1]
            +2.00456070e+01
            -6.73225700e+00 * tc[0]
            -7.45417000e-03 * tc[1]
            +8.24983167e-07 * tc[2]
            -6.01001833e-11 * tc[3]
            +1.88310200e-15 * tc[4];
        /*species 36: C3H8 */
        species[36] =
            -1.65643940e+04 / tc[1]
            +2.53627902e+01
            -7.52441520e+00 * tc[0]
            -9.44914100e-03 * tc[1]
            +1.04868402e-06 * tc[2]
            -7.68012142e-11 * tc[3]
            +2.43422390e-15 * tc[4];
        /*species 37: I-C3H7 */
        species[37] =
            +7.32271930e+03 / tc[1]
            +1.56022956e+01
            -6.51927410e+00 * tc[0]
            -8.61005200e-03 * tc[1]
            +9.56070283e-07 * tc[2]
            -7.01089433e-11 * tc[3]
            +2.22829565e-15 * tc[4];
        /*species 38: N-C3H7 */
        species[38] =
            +7.97622360e+03 / tc[1]
            +2.32250449e+01
            -7.70974790e+00 * tc[0]
            -8.01574250e-03 * tc[1]
            +8.78670633e-07 * tc[2]
            -6.32402933e-11 * tc[3]
            +1.94313595e-15 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void helmholtz(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: N2 */
        species[0] =
            -1.02089990e+03 / tc[1]
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 1: AR */
        species[1] =
            -7.45375000e+02 / tc[1]
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 2: H */
        species[2] =
            +2.54736599e+04 / tc[1]
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -3.52666409e-13 * tc[1]
            +3.32653273e-16 * tc[2]
            -1.91734693e-19 * tc[3]
            +4.63866166e-23 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.06394356e+03 / tc[1]
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.38153812e+03 / tc[1]
            +3.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
        /*species 5: O */
        species[5] =
            +2.91222592e+04 / tc[1]
            +1.16333640e-01
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 6: H2 */
        species[6] =
            -9.17935173e+02 / tc[1]
            +6.61320882e-01
            -2.34433112e+00 * tc[0]
            -3.99026037e-03 * tc[1]
            +3.24635850e-06 * tc[2]
            -1.67976745e-09 * tc[3]
            +3.68805881e-13 * tc[4];
        /*species 7: H2O */
        species[7] =
            -3.02937267e+04 / tc[1]
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +2.94808040e+02 / tc[1]
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 9: H2O2 */
        species[9] =
            -1.77025821e+04 / tc[1]
            -1.58938050e-01
            -4.27611269e+00 * tc[0]
            +2.71411208e-04 * tc[1]
            -2.78892835e-06 * tc[2]
            +1.79809011e-09 * tc[3]
            -4.31227182e-13 * tc[4];
        /*species 10: CO */
        species[10] =
            -1.43440860e+04 / tc[1]
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 11: CO2 */
        species[11] =
            -4.83719697e+04 / tc[1]
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 12: HCO */
        species[12] =
            +3.83956496e+03 / tc[1]
            -1.73186590e-01
            -4.22118584e+00 * tc[0]
            +1.62196266e-03 * tc[1]
            -2.29665743e-06 * tc[2]
            +1.10953411e-09 * tc[3]
            -2.16884432e-13 * tc[4];
        /*species 13: CH2O */
        species[13] =
            -1.43089567e+04 / tc[1]
            +3.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 14: CH4 */
        species[14] =
            -1.02466476e+04 / tc[1]
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +1.64449988e+04 / tc[1]
            +1.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 16: T-CH2 */
        species[16] =
            +4.60040401e+04 / tc[1]
            +1.20014682e+00
            -3.76267867e+00 * tc[0]
            -4.84436072e-04 * tc[1]
            -4.65816402e-07 * tc[2]
            +3.20909294e-10 * tc[3]
            -8.43708595e-14 * tc[4];
        /*species 17: S-CH2 */
        species[17] =
            +5.04968163e+04 / tc[1]
            +3.96772308e+00
            -4.19860411e+00 * tc[0]
            +1.18330710e-03 * tc[1]
            -1.37216037e-06 * tc[2]
            +5.57346651e-10 * tc[3]
            -9.71573685e-14 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +5.08977593e+03 / tc[1]
            -1.13812948e+00
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +1.30772484e+03 / tc[1]
            -3.86060362e+00
            -3.71180502e+00 * tc[0]
            +1.40231653e-03 * tc[1]
            -6.27584952e-06 * tc[2]
            +3.94226741e-09 * tc[3]
            -9.32942100e-13 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +1.28416265e+04 / tc[1]
            -1.40074356e+00
            -4.30646568e+00 * tc[0]
            +2.09329446e-03 * tc[1]
            -8.28571345e-06 * tc[2]
            +4.99272172e-09 * tc[3]
            -1.15254502e-12 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            -1.15222055e+04 / tc[1]
            +6.24601760e-01
            -4.29142492e+00 * tc[0]
            +2.75077135e-03 * tc[1]
            -9.99063813e-06 * tc[2]
            +5.90388571e-09 * tc[3]
            -1.34342886e-12 * tc[4];
        /*species 22: CH */
        species[22] =
            +7.07972934e+04 / tc[1]
            +4.05805570e-01
            -3.48981665e+00 * tc[0]
            -1.61917771e-04 * tc[1]
            +2.81498442e-07 * tc[2]
            -2.63514439e-10 * tc[3]
            +7.03045335e-14 * tc[4];
        /*species 23: C2H2 */
        species[23] =
            +2.64289807e+04 / tc[1]
            -1.41310240e+01
            -8.08681094e-01 * tc[0]
            -1.16807815e-02 * tc[1]
            +5.91953025e-06 * tc[2]
            -2.33460364e-09 * tc[3]
            +4.25036487e-13 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +3.48598468e+04 / tc[1]
            -6.29807380e+00
            -3.21246645e+00 * tc[0]
            -7.57395810e-04 * tc[1]
            -4.32015687e-06 * tc[2]
            +2.98048206e-09 * tc[3]
            -7.35754365e-13 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +3.80428530e+02 / tc[1]
            -1.93431519e+01
            -1.01340010e+00 * tc[0]
            -1.13407335e-02 * tc[1]
            +2.62232400e-06 * tc[2]
            -3.37429192e-10 * tc[3]
            -1.47995060e-14 * tc[4];
        /*species 26: C2H4O */
        species[26] =
            -2.15728780e+04 / tc[1]
            -3.73556400e-01
            -4.72945950e+00 * tc[0]
            +1.59664290e-03 * tc[1]
            -7.92248683e-06 * tc[2]
            +4.78821758e-09 * tc[3]
            -1.09655560e-12 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            -7.04291804e+03 / tc[1]
            -1.10798117e+01
            -2.13583630e+00 * tc[0]
            -9.05943605e-03 * tc[1]
            +2.89912457e-06 * tc[2]
            -7.78664640e-10 * tc[3]
            +1.00728807e-13 * tc[4];
        /*species 28: HCCO */
        species[28] =
            +2.00594490e+04 / tc[1]
            -1.12386956e+01
            -2.25172140e+00 * tc[0]
            -8.82751050e-03 * tc[1]
            +3.95485017e-06 * tc[2]
            -1.43964658e-09 * tc[3]
            +2.53324055e-13 * tc[4];
        /*species 29: C2H */
        species[29] =
            +6.68393932e+04 / tc[1]
            -4.33330705e+00
            -2.88965733e+00 * tc[0]
            -6.70498055e-03 * tc[1]
            +4.74615835e-06 * tc[2]
            -2.45659204e-09 * tc[3]
            +5.46657555e-13 * tc[4];
        /*species 30: CH2OH */
        species[30] =
            -3.52476728e+03 / tc[1]
            +1.69203330e-01
            -4.47832317e+00 * tc[0]
            +6.75348435e-04 * tc[1]
            -4.64139512e-06 * tc[2]
            +3.04056164e-09 * tc[3]
            -7.39533875e-13 * tc[4];
        /*species 31: CH3OH */
        species[31] =
            -2.56427656e+04 / tc[1]
            +6.21949405e+00
            -5.71539582e+00 * tc[0]
            +7.61545645e-03 * tc[1]
            -1.08740193e-05 * tc[2]
            +5.92339074e-09 * tc[3]
            -1.30676349e-12 * tc[4];
        /*species 32: C3H4 */
        species[32] =
            +2.15415670e+04 / tc[1]
            -8.61309450e+00
            -2.61304450e+00 * tc[0]
            -6.06128750e-03 * tc[1]
            -3.08998000e-06 * tc[2]
            +2.87709575e-09 * tc[3]
            -7.66753950e-13 * tc[4];
        /*species 33: C3H3 */
        species[33] =
            +4.01057783e+04 / tc[1]
            -1.48547831e+01
            -1.35110927e+00 * tc[0]
            -1.63705612e-02 * tc[1]
            +7.89711892e-06 * tc[2]
            -3.13591507e-09 * tc[3]
            +5.92704615e-13 * tc[4];
        /*species 34: C3H5 */
        species[34] =
            +1.92456290e+04 / tc[1]
            -1.68100305e+01
            -1.36318350e+00 * tc[0]
            -9.90691050e-03 * tc[1]
            -2.08284333e-06 * tc[2]
            +2.77962958e-09 * tc[3]
            -7.92328550e-13 * tc[4];
        /*species 35: C3H6 */
        species[35] =
            +1.07482600e+03 / tc[1]
            -1.56520330e+01
            -1.49330700e+00 * tc[0]
            -1.04625900e-02 * tc[1]
            -7.47799000e-07 * tc[2]
            +1.39076000e-09 * tc[3]
            -3.57907300e-13 * tc[4];
        /*species 36: C3H8 */
        species[36] =
            -1.40579070e+04 / tc[1]
            -1.92970271e+01
            -9.28510930e-01 * tc[0]
            -1.32302830e-02 * tc[1]
            -1.00554077e-06 * tc[2]
            +1.82624608e-09 * tc[3]
            -4.74807720e-13 * tc[4];
        /*species 37: I-C3H7 */
        species[37] =
            +9.42237240e+03 / tc[1]
            -1.96713971e+01
            -1.44491990e+00 * tc[0]
            -1.04995560e-02 * tc[1]
            -1.28393703e-06 * tc[2]
            +1.53968775e-09 * tc[3]
            -3.56414810e-13 * tc[4];
        /*species 38: N-C3H7 */
        species[38] =
            +1.03123460e+04 / tc[1]
            -2.10869167e+01
            -1.04911730e+00 * tc[0]
            -1.30044865e-02 * tc[1]
            -3.92375267e-07 * tc[2]
            +1.63292767e-09 * tc[3]
            -4.68601035e-13 * tc[4];
    } else {
        /*species 0: N2 */
        species[0] =
            -9.22797700e+02 / tc[1]
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 1: AR */
        species[1] =
            -7.45375000e+02 / tc[1]
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 2: H */
        species[2] =
            +2.54736599e+04 / tc[1]
            +1.94668292e+00
            -2.50000001e+00 * tc[0]
            +1.15421486e-11 * tc[1]
            -2.69269913e-15 * tc[2]
            +3.94596029e-19 * tc[3]
            -2.49098679e-23 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.08845772e+03 / tc[1]
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.71885774e+03 / tc[1]
            -3.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
        /*species 5: O */
        species[5] =
            +2.92175791e+04 / tc[1]
            -3.21491786e+00
            -2.56942078e+00 * tc[0]
            +4.29870569e-05 * tc[1]
            -6.99140982e-09 * tc[2]
            +8.34814992e-13 * tc[3]
            -6.14168455e-17 * tc[4];
        /*species 6: H2 */
        species[6] =
            -9.50158922e+02 / tc[1]
            +5.54230251e+00
            -3.33727920e+00 * tc[0]
            +2.47012365e-05 * tc[1]
            -8.32427963e-08 * tc[2]
            +1.49638662e-11 * tc[3]
            -1.00127688e-15 * tc[4];
        /*species 7: H2O */
        species[7] =
            -3.00042971e+04 / tc[1]
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +1.11856713e+02 / tc[1]
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 9: H2O2 */
        species[9] =
            -1.78617877e+04 / tc[1]
            +2.48846230e-01
            -4.16500285e+00 * tc[0]
            -2.45415847e-03 * tc[1]
            +3.16898708e-07 * tc[2]
            -3.09321655e-11 * tc[3]
            +1.43954153e-15 * tc[4];
        /*species 10: CO */
        species[10] =
            -1.41518724e+04 / tc[1]
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 11: CO2 */
        species[11] =
            -4.87591660e+04 / tc[1]
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 12: HCO */
        species[12] =
            +4.01191815e+03 / tc[1]
            -8.02617054e+00
            -2.77217438e+00 * tc[0]
            -2.47847763e-03 * tc[1]
            +4.14076022e-07 * tc[2]
            -4.90968148e-11 * tc[3]
            +2.66754356e-15 * tc[4];
        /*species 13: CH2O */
        species[13] =
            -1.39958323e+04 / tc[1]
            -1.28956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 14: CH4 */
        species[14] =
            -9.46834459e+03 / tc[1]
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +1.67755843e+04 / tc[1]
            -7.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 16: T-CH2 */
        species[16] =
            +4.62636040e+04 / tc[1]
            -4.29709211e+00
            -2.87410113e+00 * tc[0]
            -1.82819646e-03 * tc[1]
            +2.34824328e-07 * tc[2]
            -2.16816291e-11 * tc[3]
            +9.38637835e-16 * tc[4];
        /*species 17: S-CH2 */
        species[17] =
            +5.09259997e+04 / tc[1]
            -7.33446327e+00
            -2.29203842e+00 * tc[0]
            -2.32794318e-03 * tc[1]
            +3.35319912e-07 * tc[2]
            -3.48255000e-11 * tc[3]
            +1.69858182e-15 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +4.93988614e+03 / tc[1]
            -9.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +3.90139164e+02 / tc[1]
            +5.72459266e+00
            -4.75779238e+00 * tc[0]
            -3.72071237e-03 * tc[1]
            +4.49508627e-07 * tc[2]
            -3.65075420e-11 * tc[3]
            +1.31768549e-15 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +1.28575200e+04 / tc[1]
            -1.25077779e+01
            -1.95465642e+00 * tc[0]
            -8.69863610e-03 * tc[1]
            +1.33034445e-06 * tc[2]
            -1.46014741e-10 * tc[3]
            +7.48207880e-15 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            -1.14263932e+04 / tc[1]
            -1.50437292e+01
            -1.07188150e+00 * tc[0]
            -1.08426339e-02 * tc[1]
            +1.67093445e-06 * tc[2]
            -1.84510001e-10 * tc[3]
            +9.50014450e-15 * tc[4];
        /*species 22: CH */
        species[22] =
            +7.10124364e+04 / tc[1]
            -3.60651526e+00
            -2.87846473e+00 * tc[0]
            -4.85456840e-04 * tc[1]
            -2.40742758e-08 * tc[2]
            +1.08906541e-11 * tc[3]
            -8.80396915e-16 * tc[4];
        /*species 23: C2H2 */
        species[23] =
            +2.59359992e+04 / tc[1]
            +4.37785085e+00
            -4.14756964e+00 * tc[0]
            -2.98083332e-03 * tc[1]
            +3.95491420e-07 * tc[2]
            -3.89510143e-11 * tc[3]
            +1.80617607e-15 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +3.46128739e+04 / tc[1]
            -5.77059978e+00
            -3.01672400e+00 * tc[0]
            -5.16511460e-03 * tc[1]
            +7.80137248e-07 * tc[2]
            -8.48027400e-11 * tc[3]
            +4.31303520e-15 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            -7.31993470e+02 / tc[1]
            +6.12953670e+00
            -5.16620060e+00 * tc[0]
            -5.42391300e-03 * tc[1]
            +7.44306133e-07 * tc[2]
            -6.71904567e-11 * tc[3]
            +2.42050965e-15 * tc[4];
        /*species 26: C2H4O */
        species[26] =
            -2.25931220e+04 / tc[1]
            +7.88490250e+00
            -5.40411080e+00 * tc[0]
            -5.86152950e-03 * tc[1]
            +7.04385617e-07 * tc[2]
            -5.69770425e-11 * tc[3]
            +2.04924315e-15 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            -7.55105311e+03 / tc[1]
            +2.87905011e+00
            -4.51129732e+00 * tc[0]
            -4.50179872e-03 * tc[1]
            +6.94899392e-07 * tc[2]
            -7.69454902e-11 * tc[3]
            +3.97419100e-15 * tc[4];
        /*species 28: HCCO */
        species[28] =
            +1.93272150e+04 / tc[1]
            +8.55846530e+00
            -5.62820580e+00 * tc[0]
            -2.04267005e-03 * tc[1]
            +2.65575783e-07 * tc[2]
            -2.38550433e-11 * tc[3]
            +9.70391600e-16 * tc[4];
        /*species 29: C2H */
        species[29] =
            +6.71210650e+04 / tc[1]
            -4.46808823e+00
            -3.16780652e+00 * tc[0]
            -2.37610951e-03 * tc[1]
            +3.06311795e-07 * tc[2]
            -2.53491877e-11 * tc[3]
            +8.86163850e-16 * tc[4];
        /*species 30: CH2OH */
        species[30] =
            -4.05813228e+03 / tc[1]
            +5.94002650e+00
            -5.09312037e+00 * tc[0]
            -2.97379275e-03 * tc[1]
            +3.44160873e-07 * tc[2]
            -2.69172253e-11 * tc[3]
            +9.40625260e-16 * tc[4];
        /*species 31: CH3OH */
        species[31] =
            -2.53748747e+04 / tc[1]
            -1.37126544e+01
            -1.78970791e+00 * tc[0]
            -7.04691460e-03 * tc[1]
            +1.06083472e-06 * tc[2]
            -1.15142571e-10 * tc[3]
            +5.85301100e-15 * tc[4];
        /*species 32: C3H4 */
        species[32] =
            +2.01174950e+04 / tc[1]
            +1.63126382e+01
            -6.31687220e+00 * tc[0]
            -5.56686400e-03 * tc[1]
            +6.60489633e-07 * tc[2]
            -5.29701983e-11 * tc[3]
            +1.89377700e-15 * tc[4];
        /*species 33: C3H3 */
        species[33] =
            +3.89087427e+04 / tc[1]
            +1.87270624e+01
            -7.14221880e+00 * tc[0]
            -3.80951002e-03 * tc[1]
            +4.45766583e-07 * tc[2]
            -3.54095668e-11 * tc[3]
            +1.25737708e-15 * tc[4];
        /*species 34: C3H5 */
        species[34] =
            +1.74824490e+04 / tc[1]
            +1.67438377e+01
            -6.50078770e+00 * tc[0]
            -7.16236550e-03 * tc[1]
            +9.46360533e-07 * tc[2]
            -9.23400083e-11 * tc[3]
            +4.51819435e-15 * tc[4];
        /*species 35: C3H6 */
        species[35] =
            -9.23570300e+02 / tc[1]
            +1.90456070e+01
            -6.73225700e+00 * tc[0]
            -7.45417000e-03 * tc[1]
            +8.24983167e-07 * tc[2]
            -6.01001833e-11 * tc[3]
            +1.88310200e-15 * tc[4];
        /*species 36: C3H8 */
        species[36] =
            -1.65643940e+04 / tc[1]
            +2.43627902e+01
            -7.52441520e+00 * tc[0]
            -9.44914100e-03 * tc[1]
            +1.04868402e-06 * tc[2]
            -7.68012142e-11 * tc[3]
            +2.43422390e-15 * tc[4];
        /*species 37: I-C3H7 */
        species[37] =
            +7.32271930e+03 / tc[1]
            +1.46022956e+01
            -6.51927410e+00 * tc[0]
            -8.61005200e-03 * tc[1]
            +9.56070283e-07 * tc[2]
            -7.01089433e-11 * tc[3]
            +2.22829565e-15 * tc[4];
        /*species 38: N-C3H7 */
        species[38] =
            +7.97622360e+03 / tc[1]
            +2.22250449e+01
            -7.70974790e+00 * tc[0]
            -8.01574250e-03 * tc[1]
            +8.78670633e-07 * tc[2]
            -6.32402933e-11 * tc[3]
            +1.94313595e-15 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cv_R(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: N2 */
        species[0] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 1: AR */
        species[1] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: H */
        species[2] =
            +1.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 5: O */
        species[5] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 6: H2 */
        species[6] =
            +1.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 7: H2O */
        species[7] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 9: H2O2 */
        species[9] =
            +3.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 10: CO */
        species[10] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 11: CO2 */
        species[11] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 12: HCO */
        species[12] =
            +3.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 13: CH2O */
        species[13] =
            +3.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 14: CH4 */
        species[14] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +2.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 16: T-CH2 */
        species[16] =
            +2.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 17: S-CH2 */
        species[17] =
            +3.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +2.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +2.71180502e+00
            -2.80463306e-03 * tc[1]
            +3.76550971e-05 * tc[2]
            -4.73072089e-08 * tc[3]
            +1.86588420e-11 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +3.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +3.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 22: CH */
        species[22] =
            +2.48981665e+00
            +3.23835541e-04 * tc[1]
            -1.68899065e-06 * tc[2]
            +3.16217327e-09 * tc[3]
            -1.40609067e-12 * tc[4];
        /*species 23: C2H2 */
        species[23] =
            -1.91318906e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +2.21246645e+00
            +1.51479162e-03 * tc[1]
            +2.59209412e-05 * tc[2]
            -3.57657847e-08 * tc[3]
            +1.47150873e-11 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +1.34001000e-02
            +2.26814670e-02 * tc[1]
            -1.57339440e-05 * tc[2]
            +4.04915030e-09 * tc[3]
            +2.95990120e-13 * tc[4];
        /*species 26: C2H4O */
        species[26] =
            +3.72945950e+00
            -3.19328580e-03 * tc[1]
            +4.75349210e-05 * tc[2]
            -5.74586110e-08 * tc[3]
            +2.19311120e-11 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            +1.13583630e+00
            +1.81188721e-02 * tc[1]
            -1.73947474e-05 * tc[2]
            +9.34397568e-09 * tc[3]
            -2.01457615e-12 * tc[4];
        /*species 28: HCCO */
        species[28] =
            +1.25172140e+00
            +1.76550210e-02 * tc[1]
            -2.37291010e-05 * tc[2]
            +1.72757590e-08 * tc[3]
            -5.06648110e-12 * tc[4];
        /*species 29: C2H */
        species[29] =
            +1.88965733e+00
            +1.34099611e-02 * tc[1]
            -2.84769501e-05 * tc[2]
            +2.94791045e-08 * tc[3]
            -1.09331511e-11 * tc[4];
        /*species 30: CH2OH */
        species[30] =
            +3.47832317e+00
            -1.35069687e-03 * tc[1]
            +2.78483707e-05 * tc[2]
            -3.64867397e-08 * tc[3]
            +1.47906775e-11 * tc[4];
        /*species 31: CH3OH */
        species[31] =
            +4.71539582e+00
            -1.52309129e-02 * tc[1]
            +6.52441155e-05 * tc[2]
            -7.10806889e-08 * tc[3]
            +2.61352698e-11 * tc[4];
        /*species 32: C3H4 */
        species[32] =
            +1.61304450e+00
            +1.21225750e-02 * tc[1]
            +1.85398800e-05 * tc[2]
            -3.45251490e-08 * tc[3]
            +1.53350790e-11 * tc[4];
        /*species 33: C3H3 */
        species[33] =
            +3.51109270e-01
            +3.27411223e-02 * tc[1]
            -4.73827135e-05 * tc[2]
            +3.76309808e-08 * tc[3]
            -1.18540923e-11 * tc[4];
        /*species 34: C3H5 */
        species[34] =
            +3.63183500e-01
            +1.98138210e-02 * tc[1]
            +1.24970600e-05 * tc[2]
            -3.33555550e-08 * tc[3]
            +1.58465710e-11 * tc[4];
        /*species 35: C3H6 */
        species[35] =
            +4.93307000e-01
            +2.09251800e-02 * tc[1]
            +4.48679400e-06 * tc[2]
            -1.66891200e-08 * tc[3]
            +7.15814600e-12 * tc[4];
        /*species 36: C3H8 */
        species[36] =
            -7.14890700e-02
            +2.64605660e-02 * tc[1]
            +6.03324460e-06 * tc[2]
            -2.19149530e-08 * tc[3]
            +9.49615440e-12 * tc[4];
        /*species 37: I-C3H7 */
        species[37] =
            +4.44919900e-01
            +2.09991120e-02 * tc[1]
            +7.70362220e-06 * tc[2]
            -1.84762530e-08 * tc[3]
            +7.12829620e-12 * tc[4];
        /*species 38: N-C3H7 */
        species[38] =
            +4.91173000e-02
            +2.60089730e-02 * tc[1]
            +2.35425160e-06 * tc[2]
            -1.95951320e-08 * tc[3]
            +9.37202070e-12 * tc[4];
    } else {
        /*species 0: N2 */
        species[0] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 1: AR */
        species[1] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: H */
        species[2] =
            +1.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 5: O */
        species[5] =
            +1.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 6: H2 */
        species[6] =
            +2.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 7: H2O */
        species[7] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 9: H2O2 */
        species[9] =
            +3.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 10: CO */
        species[10] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 11: CO2 */
        species[11] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 12: HCO */
        species[12] =
            +1.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 13: CH2O */
        species[13] =
            +7.60690080e-01
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 14: CH4 */
        species[14] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +1.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 16: T-CH2 */
        species[16] =
            +1.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 17: S-CH2 */
        species[17] =
            +1.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +1.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +3.75779238e+00
            +7.44142474e-03 * tc[1]
            -2.69705176e-06 * tc[2]
            +4.38090504e-10 * tc[3]
            -2.63537098e-14 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +9.54656420e-01
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +7.18815000e-02
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 22: CH */
        species[22] =
            +1.87846473e+00
            +9.70913681e-04 * tc[1]
            +1.44445655e-07 * tc[2]
            -1.30687849e-10 * tc[3]
            +1.76079383e-14 * tc[4];
        /*species 23: C2H2 */
        species[23] =
            +3.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +2.01672400e+00
            +1.03302292e-02 * tc[1]
            -4.68082349e-06 * tc[2]
            +1.01763288e-09 * tc[3]
            -8.62607041e-14 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +4.16620060e+00
            +1.08478260e-02 * tc[1]
            -4.46583680e-06 * tc[2]
            +8.06285480e-10 * tc[3]
            -4.84101930e-14 * tc[4];
        /*species 26: C2H4O */
        species[26] =
            +4.40411080e+00
            +1.17230590e-02 * tc[1]
            -4.22631370e-06 * tc[2]
            +6.83724510e-10 * tc[3]
            -4.09848630e-14 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            +3.51129732e+00
            +9.00359745e-03 * tc[1]
            -4.16939635e-06 * tc[2]
            +9.23345882e-10 * tc[3]
            -7.94838201e-14 * tc[4];
        /*species 28: HCCO */
        species[28] =
            +4.62820580e+00
            +4.08534010e-03 * tc[1]
            -1.59345470e-06 * tc[2]
            +2.86260520e-10 * tc[3]
            -1.94078320e-14 * tc[4];
        /*species 29: C2H */
        species[29] =
            +2.16780652e+00
            +4.75221902e-03 * tc[1]
            -1.83787077e-06 * tc[2]
            +3.04190252e-10 * tc[3]
            -1.77232770e-14 * tc[4];
        /*species 30: CH2OH */
        species[30] =
            +4.09312037e+00
            +5.94758550e-03 * tc[1]
            -2.06496524e-06 * tc[2]
            +3.23006703e-10 * tc[3]
            -1.88125052e-14 * tc[4];
        /*species 31: CH3OH */
        species[31] =
            +7.89707910e-01
            +1.40938292e-02 * tc[1]
            -6.36500835e-06 * tc[2]
            +1.38171085e-09 * tc[3]
            -1.17060220e-13 * tc[4];
        /*species 32: C3H4 */
        species[32] =
            +5.31687220e+00
            +1.11337280e-02 * tc[1]
            -3.96293780e-06 * tc[2]
            +6.35642380e-10 * tc[3]
            -3.78755400e-14 * tc[4];
        /*species 33: C3H3 */
        species[33] =
            +6.14221880e+00
            +7.61902005e-03 * tc[1]
            -2.67459950e-06 * tc[2]
            +4.24914801e-10 * tc[3]
            -2.51475415e-14 * tc[4];
        /*species 34: C3H5 */
        species[34] =
            +5.50078770e+00
            +1.43247310e-02 * tc[1]
            -5.67816320e-06 * tc[2]
            +1.10808010e-09 * tc[3]
            -9.03638870e-14 * tc[4];
        /*species 35: C3H6 */
        species[35] =
            +5.73225700e+00
            +1.49083400e-02 * tc[1]
            -4.94989900e-06 * tc[2]
            +7.21202200e-10 * tc[3]
            -3.76620400e-14 * tc[4];
        /*species 36: C3H8 */
        species[36] =
            +6.52441520e+00
            +1.88982820e-02 * tc[1]
            -6.29210410e-06 * tc[2]
            +9.21614570e-10 * tc[3]
            -4.86844780e-14 * tc[4];
        /*species 37: I-C3H7 */
        species[37] =
            +5.51927410e+00
            +1.72201040e-02 * tc[1]
            -5.73642170e-06 * tc[2]
            +8.41307320e-10 * tc[3]
            -4.45659130e-14 * tc[4];
        /*species 38: N-C3H7 */
        species[38] =
            +6.70974790e+00
            +1.60314850e-02 * tc[1]
            -5.27202380e-06 * tc[2]
            +7.58883520e-10 * tc[3]
            -3.88627190e-14 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cp_R(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: N2 */
        species[0] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 1: AR */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: H */
        species[2] =
            +2.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 5: O */
        species[5] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 6: H2 */
        species[6] =
            +2.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 7: H2O */
        species[7] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 9: H2O2 */
        species[9] =
            +4.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 10: CO */
        species[10] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 11: CO2 */
        species[11] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 12: HCO */
        species[12] =
            +4.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 13: CH2O */
        species[13] =
            +4.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 14: CH4 */
        species[14] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +3.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 16: T-CH2 */
        species[16] =
            +3.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 17: S-CH2 */
        species[17] =
            +4.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +3.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +3.71180502e+00
            -2.80463306e-03 * tc[1]
            +3.76550971e-05 * tc[2]
            -4.73072089e-08 * tc[3]
            +1.86588420e-11 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +4.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +4.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 22: CH */
        species[22] =
            +3.48981665e+00
            +3.23835541e-04 * tc[1]
            -1.68899065e-06 * tc[2]
            +3.16217327e-09 * tc[3]
            -1.40609067e-12 * tc[4];
        /*species 23: C2H2 */
        species[23] =
            +8.08681094e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +3.21246645e+00
            +1.51479162e-03 * tc[1]
            +2.59209412e-05 * tc[2]
            -3.57657847e-08 * tc[3]
            +1.47150873e-11 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +1.01340010e+00
            +2.26814670e-02 * tc[1]
            -1.57339440e-05 * tc[2]
            +4.04915030e-09 * tc[3]
            +2.95990120e-13 * tc[4];
        /*species 26: C2H4O */
        species[26] =
            +4.72945950e+00
            -3.19328580e-03 * tc[1]
            +4.75349210e-05 * tc[2]
            -5.74586110e-08 * tc[3]
            +2.19311120e-11 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            +2.13583630e+00
            +1.81188721e-02 * tc[1]
            -1.73947474e-05 * tc[2]
            +9.34397568e-09 * tc[3]
            -2.01457615e-12 * tc[4];
        /*species 28: HCCO */
        species[28] =
            +2.25172140e+00
            +1.76550210e-02 * tc[1]
            -2.37291010e-05 * tc[2]
            +1.72757590e-08 * tc[3]
            -5.06648110e-12 * tc[4];
        /*species 29: C2H */
        species[29] =
            +2.88965733e+00
            +1.34099611e-02 * tc[1]
            -2.84769501e-05 * tc[2]
            +2.94791045e-08 * tc[3]
            -1.09331511e-11 * tc[4];
        /*species 30: CH2OH */
        species[30] =
            +4.47832317e+00
            -1.35069687e-03 * tc[1]
            +2.78483707e-05 * tc[2]
            -3.64867397e-08 * tc[3]
            +1.47906775e-11 * tc[4];
        /*species 31: CH3OH */
        species[31] =
            +5.71539582e+00
            -1.52309129e-02 * tc[1]
            +6.52441155e-05 * tc[2]
            -7.10806889e-08 * tc[3]
            +2.61352698e-11 * tc[4];
        /*species 32: C3H4 */
        species[32] =
            +2.61304450e+00
            +1.21225750e-02 * tc[1]
            +1.85398800e-05 * tc[2]
            -3.45251490e-08 * tc[3]
            +1.53350790e-11 * tc[4];
        /*species 33: C3H3 */
        species[33] =
            +1.35110927e+00
            +3.27411223e-02 * tc[1]
            -4.73827135e-05 * tc[2]
            +3.76309808e-08 * tc[3]
            -1.18540923e-11 * tc[4];
        /*species 34: C3H5 */
        species[34] =
            +1.36318350e+00
            +1.98138210e-02 * tc[1]
            +1.24970600e-05 * tc[2]
            -3.33555550e-08 * tc[3]
            +1.58465710e-11 * tc[4];
        /*species 35: C3H6 */
        species[35] =
            +1.49330700e+00
            +2.09251800e-02 * tc[1]
            +4.48679400e-06 * tc[2]
            -1.66891200e-08 * tc[3]
            +7.15814600e-12 * tc[4];
        /*species 36: C3H8 */
        species[36] =
            +9.28510930e-01
            +2.64605660e-02 * tc[1]
            +6.03324460e-06 * tc[2]
            -2.19149530e-08 * tc[3]
            +9.49615440e-12 * tc[4];
        /*species 37: I-C3H7 */
        species[37] =
            +1.44491990e+00
            +2.09991120e-02 * tc[1]
            +7.70362220e-06 * tc[2]
            -1.84762530e-08 * tc[3]
            +7.12829620e-12 * tc[4];
        /*species 38: N-C3H7 */
        species[38] =
            +1.04911730e+00
            +2.60089730e-02 * tc[1]
            +2.35425160e-06 * tc[2]
            -1.95951320e-08 * tc[3]
            +9.37202070e-12 * tc[4];
    } else {
        /*species 0: N2 */
        species[0] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 1: AR */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: H */
        species[2] =
            +2.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 5: O */
        species[5] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 6: H2 */
        species[6] =
            +3.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 7: H2O */
        species[7] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 9: H2O2 */
        species[9] =
            +4.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 10: CO */
        species[10] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 11: CO2 */
        species[11] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 12: HCO */
        species[12] =
            +2.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 13: CH2O */
        species[13] =
            +1.76069008e+00
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 14: CH4 */
        species[14] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +2.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 16: T-CH2 */
        species[16] =
            +2.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 17: S-CH2 */
        species[17] =
            +2.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +2.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +4.75779238e+00
            +7.44142474e-03 * tc[1]
            -2.69705176e-06 * tc[2]
            +4.38090504e-10 * tc[3]
            -2.63537098e-14 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +1.95465642e+00
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +1.07188150e+00
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 22: CH */
        species[22] =
            +2.87846473e+00
            +9.70913681e-04 * tc[1]
            +1.44445655e-07 * tc[2]
            -1.30687849e-10 * tc[3]
            +1.76079383e-14 * tc[4];
        /*species 23: C2H2 */
        species[23] =
            +4.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +3.01672400e+00
            +1.03302292e-02 * tc[1]
            -4.68082349e-06 * tc[2]
            +1.01763288e-09 * tc[3]
            -8.62607041e-14 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +5.16620060e+00
            +1.08478260e-02 * tc[1]
            -4.46583680e-06 * tc[2]
            +8.06285480e-10 * tc[3]
            -4.84101930e-14 * tc[4];
        /*species 26: C2H4O */
        species[26] =
            +5.40411080e+00
            +1.17230590e-02 * tc[1]
            -4.22631370e-06 * tc[2]
            +6.83724510e-10 * tc[3]
            -4.09848630e-14 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            +4.51129732e+00
            +9.00359745e-03 * tc[1]
            -4.16939635e-06 * tc[2]
            +9.23345882e-10 * tc[3]
            -7.94838201e-14 * tc[4];
        /*species 28: HCCO */
        species[28] =
            +5.62820580e+00
            +4.08534010e-03 * tc[1]
            -1.59345470e-06 * tc[2]
            +2.86260520e-10 * tc[3]
            -1.94078320e-14 * tc[4];
        /*species 29: C2H */
        species[29] =
            +3.16780652e+00
            +4.75221902e-03 * tc[1]
            -1.83787077e-06 * tc[2]
            +3.04190252e-10 * tc[3]
            -1.77232770e-14 * tc[4];
        /*species 30: CH2OH */
        species[30] =
            +5.09312037e+00
            +5.94758550e-03 * tc[1]
            -2.06496524e-06 * tc[2]
            +3.23006703e-10 * tc[3]
            -1.88125052e-14 * tc[4];
        /*species 31: CH3OH */
        species[31] =
            +1.78970791e+00
            +1.40938292e-02 * tc[1]
            -6.36500835e-06 * tc[2]
            +1.38171085e-09 * tc[3]
            -1.17060220e-13 * tc[4];
        /*species 32: C3H4 */
        species[32] =
            +6.31687220e+00
            +1.11337280e-02 * tc[1]
            -3.96293780e-06 * tc[2]
            +6.35642380e-10 * tc[3]
            -3.78755400e-14 * tc[4];
        /*species 33: C3H3 */
        species[33] =
            +7.14221880e+00
            +7.61902005e-03 * tc[1]
            -2.67459950e-06 * tc[2]
            +4.24914801e-10 * tc[3]
            -2.51475415e-14 * tc[4];
        /*species 34: C3H5 */
        species[34] =
            +6.50078770e+00
            +1.43247310e-02 * tc[1]
            -5.67816320e-06 * tc[2]
            +1.10808010e-09 * tc[3]
            -9.03638870e-14 * tc[4];
        /*species 35: C3H6 */
        species[35] =
            +6.73225700e+00
            +1.49083400e-02 * tc[1]
            -4.94989900e-06 * tc[2]
            +7.21202200e-10 * tc[3]
            -3.76620400e-14 * tc[4];
        /*species 36: C3H8 */
        species[36] =
            +7.52441520e+00
            +1.88982820e-02 * tc[1]
            -6.29210410e-06 * tc[2]
            +9.21614570e-10 * tc[3]
            -4.86844780e-14 * tc[4];
        /*species 37: I-C3H7 */
        species[37] =
            +6.51927410e+00
            +1.72201040e-02 * tc[1]
            -5.73642170e-06 * tc[2]
            +8.41307320e-10 * tc[3]
            -4.45659130e-14 * tc[4];
        /*species 38: N-C3H7 */
        species[38] =
            +7.70974790e+00
            +1.60314850e-02 * tc[1]
            -5.27202380e-06 * tc[2]
            +7.58883520e-10 * tc[3]
            -3.88627190e-14 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesInternalEnergy(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: N2 */
        species[0] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
        /*species 1: AR */
        species[1] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
        /*species 2: H */
        species[2] =
            +1.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.38153812e+03 / tc[1];
        /*species 5: O */
        species[5] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 / tc[1];
        /*species 6: H2 */
        species[6] =
            +1.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 / tc[1];
        /*species 7: H2O */
        species[7] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 / tc[1];
        /*species 8: HO2 */
        species[8] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 / tc[1];
        /*species 9: H2O2 */
        species[9] =
            +3.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 / tc[1];
        /*species 10: CO */
        species[10] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 / tc[1];
        /*species 11: CO2 */
        species[11] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 / tc[1];
        /*species 12: HCO */
        species[12] =
            +3.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 / tc[1];
        /*species 13: CH2O */
        species[13] =
            +3.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 / tc[1];
        /*species 14: CH4 */
        species[14] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 / tc[1];
        /*species 15: CH3 */
        species[15] =
            +2.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 / tc[1];
        /*species 16: T-CH2 */
        species[16] =
            +2.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 / tc[1];
        /*species 17: S-CH2 */
        species[17] =
            +3.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 / tc[1];
        /*species 18: C2H4 */
        species[18] =
            +2.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 / tc[1];
        /*species 19: CH3O */
        species[19] =
            +2.71180502e+00
            -1.40231653e-03 * tc[1]
            +1.25516990e-05 * tc[2]
            -1.18268022e-08 * tc[3]
            +3.73176840e-12 * tc[4]
            +1.30772484e+03 / tc[1];
        /*species 20: C2H5 */
        species[20] =
            +3.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 / tc[1];
        /*species 21: C2H6 */
        species[21] =
            +3.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 / tc[1];
        /*species 22: CH */
        species[22] =
            +2.48981665e+00
            +1.61917771e-04 * tc[1]
            -5.62996883e-07 * tc[2]
            +7.90543317e-10 * tc[3]
            -2.81218134e-13 * tc[4]
            +7.07972934e+04 / tc[1];
        /*species 23: C2H2 */
        species[23] =
            -1.91318906e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 / tc[1];
        /*species 24: C2H3 */
        species[24] =
            +2.21246645e+00
            +7.57395810e-04 * tc[1]
            +8.64031373e-06 * tc[2]
            -8.94144617e-09 * tc[3]
            +2.94301746e-12 * tc[4]
            +3.48598468e+04 / tc[1];
        /*species 25: CH2CHO */
        species[25] =
            +1.34001000e-02
            +1.13407335e-02 * tc[1]
            -5.24464800e-06 * tc[2]
            +1.01228758e-09 * tc[3]
            +5.91980240e-14 * tc[4]
            +3.80428530e+02 / tc[1];
        /*species 26: C2H4O */
        species[26] =
            +3.72945950e+00
            -1.59664290e-03 * tc[1]
            +1.58449737e-05 * tc[2]
            -1.43646527e-08 * tc[3]
            +4.38622240e-12 * tc[4]
            -2.15728780e+04 / tc[1];
        /*species 27: CH2CO */
        species[27] =
            +1.13583630e+00
            +9.05943605e-03 * tc[1]
            -5.79824913e-06 * tc[2]
            +2.33599392e-09 * tc[3]
            -4.02915230e-13 * tc[4]
            -7.04291804e+03 / tc[1];
        /*species 28: HCCO */
        species[28] =
            +1.25172140e+00
            +8.82751050e-03 * tc[1]
            -7.90970033e-06 * tc[2]
            +4.31893975e-09 * tc[3]
            -1.01329622e-12 * tc[4]
            +2.00594490e+04 / tc[1];
        /*species 29: C2H */
        species[29] =
            +1.88965733e+00
            +6.70498055e-03 * tc[1]
            -9.49231670e-06 * tc[2]
            +7.36977613e-09 * tc[3]
            -2.18663022e-12 * tc[4]
            +6.68393932e+04 / tc[1];
        /*species 30: CH2OH */
        species[30] =
            +3.47832317e+00
            -6.75348435e-04 * tc[1]
            +9.28279023e-06 * tc[2]
            -9.12168492e-09 * tc[3]
            +2.95813550e-12 * tc[4]
            -3.52476728e+03 / tc[1];
        /*species 31: CH3OH */
        species[31] =
            +4.71539582e+00
            -7.61545645e-03 * tc[1]
            +2.17480385e-05 * tc[2]
            -1.77701722e-08 * tc[3]
            +5.22705396e-12 * tc[4]
            -2.56427656e+04 / tc[1];
        /*species 32: C3H4 */
        species[32] =
            +1.61304450e+00
            +6.06128750e-03 * tc[1]
            +6.17996000e-06 * tc[2]
            -8.63128725e-09 * tc[3]
            +3.06701580e-12 * tc[4]
            +2.15415670e+04 / tc[1];
        /*species 33: C3H3 */
        species[33] =
            +3.51109270e-01
            +1.63705612e-02 * tc[1]
            -1.57942378e-05 * tc[2]
            +9.40774520e-09 * tc[3]
            -2.37081846e-12 * tc[4]
            +4.01057783e+04 / tc[1];
        /*species 34: C3H5 */
        species[34] =
            +3.63183500e-01
            +9.90691050e-03 * tc[1]
            +4.16568667e-06 * tc[2]
            -8.33888875e-09 * tc[3]
            +3.16931420e-12 * tc[4]
            +1.92456290e+04 / tc[1];
        /*species 35: C3H6 */
        species[35] =
            +4.93307000e-01
            +1.04625900e-02 * tc[1]
            +1.49559800e-06 * tc[2]
            -4.17228000e-09 * tc[3]
            +1.43162920e-12 * tc[4]
            +1.07482600e+03 / tc[1];
        /*species 36: C3H8 */
        species[36] =
            -7.14890700e-02
            +1.32302830e-02 * tc[1]
            +2.01108153e-06 * tc[2]
            -5.47873825e-09 * tc[3]
            +1.89923088e-12 * tc[4]
            -1.40579070e+04 / tc[1];
        /*species 37: I-C3H7 */
        species[37] =
            +4.44919900e-01
            +1.04995560e-02 * tc[1]
            +2.56787407e-06 * tc[2]
            -4.61906325e-09 * tc[3]
            +1.42565924e-12 * tc[4]
            +9.42237240e+03 / tc[1];
        /*species 38: N-C3H7 */
        species[38] =
            +4.91173000e-02
            +1.30044865e-02 * tc[1]
            +7.84750533e-07 * tc[2]
            -4.89878300e-09 * tc[3]
            +1.87440414e-12 * tc[4]
            +1.03123460e+04 / tc[1];
    } else {
        /*species 0: N2 */
        species[0] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
        /*species 1: AR */
        species[1] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
        /*species 2: H */
        species[2] =
            +1.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.71885774e+03 / tc[1];
        /*species 5: O */
        species[5] =
            +1.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 / tc[1];
        /*species 6: H2 */
        species[6] =
            +2.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 / tc[1];
        /*species 7: H2O */
        species[7] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 / tc[1];
        /*species 8: HO2 */
        species[8] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 / tc[1];
        /*species 9: H2O2 */
        species[9] =
            +3.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 / tc[1];
        /*species 10: CO */
        species[10] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 / tc[1];
        /*species 11: CO2 */
        species[11] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 / tc[1];
        /*species 12: HCO */
        species[12] =
            +1.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 / tc[1];
        /*species 13: CH2O */
        species[13] =
            +7.60690080e-01
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 / tc[1];
        /*species 14: CH4 */
        species[14] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 / tc[1];
        /*species 15: CH3 */
        species[15] =
            +1.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 / tc[1];
        /*species 16: T-CH2 */
        species[16] =
            +1.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 / tc[1];
        /*species 17: S-CH2 */
        species[17] =
            +1.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 / tc[1];
        /*species 18: C2H4 */
        species[18] =
            +1.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 / tc[1];
        /*species 19: CH3O */
        species[19] =
            +3.75779238e+00
            +3.72071237e-03 * tc[1]
            -8.99017253e-07 * tc[2]
            +1.09522626e-10 * tc[3]
            -5.27074196e-15 * tc[4]
            +3.90139164e+02 / tc[1];
        /*species 20: C2H5 */
        species[20] =
            +9.54656420e-01
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 / tc[1];
        /*species 21: C2H6 */
        species[21] =
            +7.18815000e-02
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 / tc[1];
        /*species 22: CH */
        species[22] =
            +1.87846473e+00
            +4.85456840e-04 * tc[1]
            +4.81485517e-08 * tc[2]
            -3.26719623e-11 * tc[3]
            +3.52158766e-15 * tc[4]
            +7.10124364e+04 / tc[1];
        /*species 23: C2H2 */
        species[23] =
            +3.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 / tc[1];
        /*species 24: C2H3 */
        species[24] =
            +2.01672400e+00
            +5.16511460e-03 * tc[1]
            -1.56027450e-06 * tc[2]
            +2.54408220e-10 * tc[3]
            -1.72521408e-14 * tc[4]
            +3.46128739e+04 / tc[1];
        /*species 25: CH2CHO */
        species[25] =
            +4.16620060e+00
            +5.42391300e-03 * tc[1]
            -1.48861227e-06 * tc[2]
            +2.01571370e-10 * tc[3]
            -9.68203860e-15 * tc[4]
            -7.31993470e+02 / tc[1];
        /*species 26: C2H4O */
        species[26] =
            +4.40411080e+00
            +5.86152950e-03 * tc[1]
            -1.40877123e-06 * tc[2]
            +1.70931128e-10 * tc[3]
            -8.19697260e-15 * tc[4]
            -2.25931220e+04 / tc[1];
        /*species 27: CH2CO */
        species[27] =
            +3.51129732e+00
            +4.50179872e-03 * tc[1]
            -1.38979878e-06 * tc[2]
            +2.30836470e-10 * tc[3]
            -1.58967640e-14 * tc[4]
            -7.55105311e+03 / tc[1];
        /*species 28: HCCO */
        species[28] =
            +4.62820580e+00
            +2.04267005e-03 * tc[1]
            -5.31151567e-07 * tc[2]
            +7.15651300e-11 * tc[3]
            -3.88156640e-15 * tc[4]
            +1.93272150e+04 / tc[1];
        /*species 29: C2H */
        species[29] =
            +2.16780652e+00
            +2.37610951e-03 * tc[1]
            -6.12623590e-07 * tc[2]
            +7.60475630e-11 * tc[3]
            -3.54465540e-15 * tc[4]
            +6.71210650e+04 / tc[1];
        /*species 30: CH2OH */
        species[30] =
            +4.09312037e+00
            +2.97379275e-03 * tc[1]
            -6.88321747e-07 * tc[2]
            +8.07516758e-11 * tc[3]
            -3.76250104e-15 * tc[4]
            -4.05813228e+03 / tc[1];
        /*species 31: CH3OH */
        species[31] =
            +7.89707910e-01
            +7.04691460e-03 * tc[1]
            -2.12166945e-06 * tc[2]
            +3.45427713e-10 * tc[3]
            -2.34120440e-14 * tc[4]
            -2.53748747e+04 / tc[1];
        /*species 32: C3H4 */
        species[32] =
            +5.31687220e+00
            +5.56686400e-03 * tc[1]
            -1.32097927e-06 * tc[2]
            +1.58910595e-10 * tc[3]
            -7.57510800e-15 * tc[4]
            +2.01174950e+04 / tc[1];
        /*species 33: C3H3 */
        species[33] =
            +6.14221880e+00
            +3.80951002e-03 * tc[1]
            -8.91533167e-07 * tc[2]
            +1.06228700e-10 * tc[3]
            -5.02950830e-15 * tc[4]
            +3.89087427e+04 / tc[1];
        /*species 34: C3H5 */
        species[34] =
            +5.50078770e+00
            +7.16236550e-03 * tc[1]
            -1.89272107e-06 * tc[2]
            +2.77020025e-10 * tc[3]
            -1.80727774e-14 * tc[4]
            +1.74824490e+04 / tc[1];
        /*species 35: C3H6 */
        species[35] =
            +5.73225700e+00
            +7.45417000e-03 * tc[1]
            -1.64996633e-06 * tc[2]
            +1.80300550e-10 * tc[3]
            -7.53240800e-15 * tc[4]
            -9.23570300e+02 / tc[1];
        /*species 36: C3H8 */
        species[36] =
            +6.52441520e+00
            +9.44914100e-03 * tc[1]
            -2.09736803e-06 * tc[2]
            +2.30403642e-10 * tc[3]
            -9.73689560e-15 * tc[4]
            -1.65643940e+04 / tc[1];
        /*species 37: I-C3H7 */
        species[37] =
            +5.51927410e+00
            +8.61005200e-03 * tc[1]
            -1.91214057e-06 * tc[2]
            +2.10326830e-10 * tc[3]
            -8.91318260e-15 * tc[4]
            +7.32271930e+03 / tc[1];
        /*species 38: N-C3H7 */
        species[38] =
            +6.70974790e+00
            +8.01574250e-03 * tc[1]
            -1.75734127e-06 * tc[2]
            +1.89720880e-10 * tc[3]
            -7.77254380e-15 * tc[4]
            +7.97622360e+03 / tc[1];
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: N2 */
        species[0] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
        /*species 1: AR */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
        /*species 2: H */
        species[2] =
            +2.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.38153812e+03 / tc[1];
        /*species 5: O */
        species[5] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 / tc[1];
        /*species 6: H2 */
        species[6] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 / tc[1];
        /*species 7: H2O */
        species[7] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 / tc[1];
        /*species 8: HO2 */
        species[8] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 / tc[1];
        /*species 9: H2O2 */
        species[9] =
            +4.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 / tc[1];
        /*species 10: CO */
        species[10] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 / tc[1];
        /*species 11: CO2 */
        species[11] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 / tc[1];
        /*species 12: HCO */
        species[12] =
            +4.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 / tc[1];
        /*species 13: CH2O */
        species[13] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 / tc[1];
        /*species 14: CH4 */
        species[14] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 / tc[1];
        /*species 15: CH3 */
        species[15] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 / tc[1];
        /*species 16: T-CH2 */
        species[16] =
            +3.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 / tc[1];
        /*species 17: S-CH2 */
        species[17] =
            +4.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 / tc[1];
        /*species 18: C2H4 */
        species[18] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 / tc[1];
        /*species 19: CH3O */
        species[19] =
            +3.71180502e+00
            -1.40231653e-03 * tc[1]
            +1.25516990e-05 * tc[2]
            -1.18268022e-08 * tc[3]
            +3.73176840e-12 * tc[4]
            +1.30772484e+03 / tc[1];
        /*species 20: C2H5 */
        species[20] =
            +4.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 / tc[1];
        /*species 21: C2H6 */
        species[21] =
            +4.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 / tc[1];
        /*species 22: CH */
        species[22] =
            +3.48981665e+00
            +1.61917771e-04 * tc[1]
            -5.62996883e-07 * tc[2]
            +7.90543317e-10 * tc[3]
            -2.81218134e-13 * tc[4]
            +7.07972934e+04 / tc[1];
        /*species 23: C2H2 */
        species[23] =
            +8.08681094e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 / tc[1];
        /*species 24: C2H3 */
        species[24] =
            +3.21246645e+00
            +7.57395810e-04 * tc[1]
            +8.64031373e-06 * tc[2]
            -8.94144617e-09 * tc[3]
            +2.94301746e-12 * tc[4]
            +3.48598468e+04 / tc[1];
        /*species 25: CH2CHO */
        species[25] =
            +1.01340010e+00
            +1.13407335e-02 * tc[1]
            -5.24464800e-06 * tc[2]
            +1.01228758e-09 * tc[3]
            +5.91980240e-14 * tc[4]
            +3.80428530e+02 / tc[1];
        /*species 26: C2H4O */
        species[26] =
            +4.72945950e+00
            -1.59664290e-03 * tc[1]
            +1.58449737e-05 * tc[2]
            -1.43646527e-08 * tc[3]
            +4.38622240e-12 * tc[4]
            -2.15728780e+04 / tc[1];
        /*species 27: CH2CO */
        species[27] =
            +2.13583630e+00
            +9.05943605e-03 * tc[1]
            -5.79824913e-06 * tc[2]
            +2.33599392e-09 * tc[3]
            -4.02915230e-13 * tc[4]
            -7.04291804e+03 / tc[1];
        /*species 28: HCCO */
        species[28] =
            +2.25172140e+00
            +8.82751050e-03 * tc[1]
            -7.90970033e-06 * tc[2]
            +4.31893975e-09 * tc[3]
            -1.01329622e-12 * tc[4]
            +2.00594490e+04 / tc[1];
        /*species 29: C2H */
        species[29] =
            +2.88965733e+00
            +6.70498055e-03 * tc[1]
            -9.49231670e-06 * tc[2]
            +7.36977613e-09 * tc[3]
            -2.18663022e-12 * tc[4]
            +6.68393932e+04 / tc[1];
        /*species 30: CH2OH */
        species[30] =
            +4.47832317e+00
            -6.75348435e-04 * tc[1]
            +9.28279023e-06 * tc[2]
            -9.12168492e-09 * tc[3]
            +2.95813550e-12 * tc[4]
            -3.52476728e+03 / tc[1];
        /*species 31: CH3OH */
        species[31] =
            +5.71539582e+00
            -7.61545645e-03 * tc[1]
            +2.17480385e-05 * tc[2]
            -1.77701722e-08 * tc[3]
            +5.22705396e-12 * tc[4]
            -2.56427656e+04 / tc[1];
        /*species 32: C3H4 */
        species[32] =
            +2.61304450e+00
            +6.06128750e-03 * tc[1]
            +6.17996000e-06 * tc[2]
            -8.63128725e-09 * tc[3]
            +3.06701580e-12 * tc[4]
            +2.15415670e+04 / tc[1];
        /*species 33: C3H3 */
        species[33] =
            +1.35110927e+00
            +1.63705612e-02 * tc[1]
            -1.57942378e-05 * tc[2]
            +9.40774520e-09 * tc[3]
            -2.37081846e-12 * tc[4]
            +4.01057783e+04 / tc[1];
        /*species 34: C3H5 */
        species[34] =
            +1.36318350e+00
            +9.90691050e-03 * tc[1]
            +4.16568667e-06 * tc[2]
            -8.33888875e-09 * tc[3]
            +3.16931420e-12 * tc[4]
            +1.92456290e+04 / tc[1];
        /*species 35: C3H6 */
        species[35] =
            +1.49330700e+00
            +1.04625900e-02 * tc[1]
            +1.49559800e-06 * tc[2]
            -4.17228000e-09 * tc[3]
            +1.43162920e-12 * tc[4]
            +1.07482600e+03 / tc[1];
        /*species 36: C3H8 */
        species[36] =
            +9.28510930e-01
            +1.32302830e-02 * tc[1]
            +2.01108153e-06 * tc[2]
            -5.47873825e-09 * tc[3]
            +1.89923088e-12 * tc[4]
            -1.40579070e+04 / tc[1];
        /*species 37: I-C3H7 */
        species[37] =
            +1.44491990e+00
            +1.04995560e-02 * tc[1]
            +2.56787407e-06 * tc[2]
            -4.61906325e-09 * tc[3]
            +1.42565924e-12 * tc[4]
            +9.42237240e+03 / tc[1];
        /*species 38: N-C3H7 */
        species[38] =
            +1.04911730e+00
            +1.30044865e-02 * tc[1]
            +7.84750533e-07 * tc[2]
            -4.89878300e-09 * tc[3]
            +1.87440414e-12 * tc[4]
            +1.03123460e+04 / tc[1];
    } else {
        /*species 0: N2 */
        species[0] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
        /*species 1: AR */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
        /*species 2: H */
        species[2] =
            +2.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.71885774e+03 / tc[1];
        /*species 5: O */
        species[5] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 / tc[1];
        /*species 6: H2 */
        species[6] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 / tc[1];
        /*species 7: H2O */
        species[7] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 / tc[1];
        /*species 8: HO2 */
        species[8] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 / tc[1];
        /*species 9: H2O2 */
        species[9] =
            +4.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 / tc[1];
        /*species 10: CO */
        species[10] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 / tc[1];
        /*species 11: CO2 */
        species[11] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 / tc[1];
        /*species 12: HCO */
        species[12] =
            +2.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 / tc[1];
        /*species 13: CH2O */
        species[13] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 / tc[1];
        /*species 14: CH4 */
        species[14] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 / tc[1];
        /*species 15: CH3 */
        species[15] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 / tc[1];
        /*species 16: T-CH2 */
        species[16] =
            +2.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 / tc[1];
        /*species 17: S-CH2 */
        species[17] =
            +2.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 / tc[1];
        /*species 18: C2H4 */
        species[18] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 / tc[1];
        /*species 19: CH3O */
        species[19] =
            +4.75779238e+00
            +3.72071237e-03 * tc[1]
            -8.99017253e-07 * tc[2]
            +1.09522626e-10 * tc[3]
            -5.27074196e-15 * tc[4]
            +3.90139164e+02 / tc[1];
        /*species 20: C2H5 */
        species[20] =
            +1.95465642e+00
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 / tc[1];
        /*species 21: C2H6 */
        species[21] =
            +1.07188150e+00
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 / tc[1];
        /*species 22: CH */
        species[22] =
            +2.87846473e+00
            +4.85456840e-04 * tc[1]
            +4.81485517e-08 * tc[2]
            -3.26719623e-11 * tc[3]
            +3.52158766e-15 * tc[4]
            +7.10124364e+04 / tc[1];
        /*species 23: C2H2 */
        species[23] =
            +4.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 / tc[1];
        /*species 24: C2H3 */
        species[24] =
            +3.01672400e+00
            +5.16511460e-03 * tc[1]
            -1.56027450e-06 * tc[2]
            +2.54408220e-10 * tc[3]
            -1.72521408e-14 * tc[4]
            +3.46128739e+04 / tc[1];
        /*species 25: CH2CHO */
        species[25] =
            +5.16620060e+00
            +5.42391300e-03 * tc[1]
            -1.48861227e-06 * tc[2]
            +2.01571370e-10 * tc[3]
            -9.68203860e-15 * tc[4]
            -7.31993470e+02 / tc[1];
        /*species 26: C2H4O */
        species[26] =
            +5.40411080e+00
            +5.86152950e-03 * tc[1]
            -1.40877123e-06 * tc[2]
            +1.70931128e-10 * tc[3]
            -8.19697260e-15 * tc[4]
            -2.25931220e+04 / tc[1];
        /*species 27: CH2CO */
        species[27] =
            +4.51129732e+00
            +4.50179872e-03 * tc[1]
            -1.38979878e-06 * tc[2]
            +2.30836470e-10 * tc[3]
            -1.58967640e-14 * tc[4]
            -7.55105311e+03 / tc[1];
        /*species 28: HCCO */
        species[28] =
            +5.62820580e+00
            +2.04267005e-03 * tc[1]
            -5.31151567e-07 * tc[2]
            +7.15651300e-11 * tc[3]
            -3.88156640e-15 * tc[4]
            +1.93272150e+04 / tc[1];
        /*species 29: C2H */
        species[29] =
            +3.16780652e+00
            +2.37610951e-03 * tc[1]
            -6.12623590e-07 * tc[2]
            +7.60475630e-11 * tc[3]
            -3.54465540e-15 * tc[4]
            +6.71210650e+04 / tc[1];
        /*species 30: CH2OH */
        species[30] =
            +5.09312037e+00
            +2.97379275e-03 * tc[1]
            -6.88321747e-07 * tc[2]
            +8.07516758e-11 * tc[3]
            -3.76250104e-15 * tc[4]
            -4.05813228e+03 / tc[1];
        /*species 31: CH3OH */
        species[31] =
            +1.78970791e+00
            +7.04691460e-03 * tc[1]
            -2.12166945e-06 * tc[2]
            +3.45427713e-10 * tc[3]
            -2.34120440e-14 * tc[4]
            -2.53748747e+04 / tc[1];
        /*species 32: C3H4 */
        species[32] =
            +6.31687220e+00
            +5.56686400e-03 * tc[1]
            -1.32097927e-06 * tc[2]
            +1.58910595e-10 * tc[3]
            -7.57510800e-15 * tc[4]
            +2.01174950e+04 / tc[1];
        /*species 33: C3H3 */
        species[33] =
            +7.14221880e+00
            +3.80951002e-03 * tc[1]
            -8.91533167e-07 * tc[2]
            +1.06228700e-10 * tc[3]
            -5.02950830e-15 * tc[4]
            +3.89087427e+04 / tc[1];
        /*species 34: C3H5 */
        species[34] =
            +6.50078770e+00
            +7.16236550e-03 * tc[1]
            -1.89272107e-06 * tc[2]
            +2.77020025e-10 * tc[3]
            -1.80727774e-14 * tc[4]
            +1.74824490e+04 / tc[1];
        /*species 35: C3H6 */
        species[35] =
            +6.73225700e+00
            +7.45417000e-03 * tc[1]
            -1.64996633e-06 * tc[2]
            +1.80300550e-10 * tc[3]
            -7.53240800e-15 * tc[4]
            -9.23570300e+02 / tc[1];
        /*species 36: C3H8 */
        species[36] =
            +7.52441520e+00
            +9.44914100e-03 * tc[1]
            -2.09736803e-06 * tc[2]
            +2.30403642e-10 * tc[3]
            -9.73689560e-15 * tc[4]
            -1.65643940e+04 / tc[1];
        /*species 37: I-C3H7 */
        species[37] =
            +6.51927410e+00
            +8.61005200e-03 * tc[1]
            -1.91214057e-06 * tc[2]
            +2.10326830e-10 * tc[3]
            -8.91318260e-15 * tc[4]
            +7.32271930e+03 / tc[1];
        /*species 38: N-C3H7 */
        species[38] =
            +7.70974790e+00
            +8.01574250e-03 * tc[1]
            -1.75734127e-06 * tc[2]
            +1.89720880e-10 * tc[3]
            -7.77254380e-15 * tc[4]
            +7.97622360e+03 / tc[1];
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEntropy(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: N2 */
        species[0] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 1: AR */
        species[1] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600000e+00 ;
        /*species 2: H */
        species[2] =
            +2.50000000e+00 * tc[0]
            +7.05332819e-13 * tc[1]
            -9.97959820e-16 * tc[2]
            +7.66938773e-19 * tc[3]
            -2.31933083e-22 * tc[4]
            -4.46682853e-01 ;
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
            +4.12530561e+00 * tc[0]
            -3.22544939e-03 * tc[1]
            +3.26382346e-06 * tc[2]
            -1.93284548e-09 * tc[3]
            +5.15593447e-13 * tc[4]
            -6.90432960e-01 ;
        /*species 5: O */
        species[5] =
            +3.16826710e+00 * tc[0]
            -3.27931884e-03 * tc[1]
            +3.32153198e-06 * tc[2]
            -2.04268875e-09 * tc[3]
            +5.28164927e-13 * tc[4]
            +2.05193346e+00 ;
        /*species 6: H2 */
        species[6] =
            +2.34433112e+00 * tc[0]
            +7.98052075e-03 * tc[1]
            -9.73907550e-06 * tc[2]
            +6.71906980e-09 * tc[3]
            -1.84402940e-12 * tc[4]
            +6.83010238e-01 ;
        /*species 7: H2O */
        species[7] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
        /*species 8: HO2 */
        species[8] =
            +4.30179801e+00 * tc[0]
            -4.74912051e-03 * tc[1]
            +1.05791445e-05 * tc[2]
            -8.09212980e-09 * tc[3]
            +2.32306281e-12 * tc[4]
            +3.71666245e+00 ;
        /*species 9: H2O2 */
        species[9] =
            +4.27611269e+00 * tc[0]
            -5.42822417e-04 * tc[1]
            +8.36678505e-06 * tc[2]
            -7.19236043e-09 * tc[3]
            +2.15613591e-12 * tc[4]
            +3.43505074e+00 ;
        /*species 10: CO */
        species[10] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 11: CO2 */
        species[11] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
        /*species 12: HCO */
        species[12] =
            +4.22118584e+00 * tc[0]
            -3.24392532e-03 * tc[1]
            +6.88997230e-06 * tc[2]
            -4.43813643e-09 * tc[3]
            +1.08442216e-12 * tc[4]
            +3.39437243e+00 ;
        /*species 13: CH2O */
        species[13] =
            +4.79372315e+00 * tc[0]
            -9.90833369e-03 * tc[1]
            +1.86610004e-05 * tc[2]
            -1.26428420e-08 * tc[3]
            +3.29431630e-12 * tc[4]
            +6.02812900e-01 ;
        /*species 14: CH4 */
        species[14] =
            +5.14987613e+00 * tc[0]
            -1.36709788e-02 * tc[1]
            +2.45900299e-05 * tc[2]
            -1.61581009e-08 * tc[3]
            +4.16734890e-12 * tc[4]
            -4.64130376e+00 ;
        /*species 15: CH3 */
        species[15] =
            +3.67359040e+00 * tc[0]
            +2.01095175e-03 * tc[1]
            +2.86510928e-06 * tc[2]
            -2.29039142e-09 * tc[3]
            +6.35964335e-13 * tc[4]
            +1.60456433e+00 ;
        /*species 16: T-CH2 */
        species[16] =
            +3.76267867e+00 * tc[0]
            +9.68872143e-04 * tc[1]
            +1.39744921e-06 * tc[2]
            -1.28363718e-09 * tc[3]
            +4.21854298e-13 * tc[4]
            +1.56253185e+00 ;
        /*species 17: S-CH2 */
        species[17] =
            +4.19860411e+00 * tc[0]
            -2.36661419e-03 * tc[1]
            +4.11648110e-06 * tc[2]
            -2.22938660e-09 * tc[3]
            +4.85786843e-13 * tc[4]
            -7.69118967e-01 ;
        /*species 18: C2H4 */
        species[18] =
            +3.95920148e+00 * tc[0]
            -7.57052247e-03 * tc[1]
            +2.85495146e-05 * tc[2]
            -2.30529584e-08 * tc[3]
            +6.74710933e-12 * tc[4]
            +4.09733096e+00 ;
        /*species 19: CH3O */
        species[19] =
            +3.71180502e+00 * tc[0]
            -2.80463306e-03 * tc[1]
            +1.88275486e-05 * tc[2]
            -1.57690696e-08 * tc[3]
            +4.66471050e-12 * tc[4]
            +6.57240864e+00 ;
        /*species 20: C2H5 */
        species[20] =
            +4.30646568e+00 * tc[0]
            -4.18658892e-03 * tc[1]
            +2.48571403e-05 * tc[2]
            -1.99708869e-08 * tc[3]
            +5.76272510e-12 * tc[4]
            +4.70720924e+00 ;
        /*species 21: C2H6 */
        species[21] =
            +4.29142492e+00 * tc[0]
            -5.50154270e-03 * tc[1]
            +2.99719144e-05 * tc[2]
            -2.36155428e-08 * tc[3]
            +6.71714427e-12 * tc[4]
            +2.66682316e+00 ;
        /*species 22: CH */
        species[22] =
            +3.48981665e+00 * tc[0]
            +3.23835541e-04 * tc[1]
            -8.44495325e-07 * tc[2]
            +1.05405776e-09 * tc[3]
            -3.51522668e-13 * tc[4]
            +2.08401108e+00 ;
        /*species 23: C2H2 */
        species[23] =
            +8.08681094e-01 * tc[0]
            +2.33615629e-02 * tc[1]
            -1.77585907e-05 * tc[2]
            +9.33841457e-09 * tc[3]
            -2.12518243e-12 * tc[4]
            +1.39397051e+01 ;
        /*species 24: C2H3 */
        species[24] =
            +3.21246645e+00 * tc[0]
            +1.51479162e-03 * tc[1]
            +1.29604706e-05 * tc[2]
            -1.19219282e-08 * tc[3]
            +3.67877182e-12 * tc[4]
            +8.51054025e+00 ;
        /*species 25: CH2CHO */
        species[25] =
            +1.01340010e+00 * tc[0]
            +2.26814670e-02 * tc[1]
            -7.86697200e-06 * tc[2]
            +1.34971677e-09 * tc[3]
            +7.39975300e-14 * tc[4]
            +1.93565520e+01 ;
        /*species 26: C2H4O */
        species[26] =
            +4.72945950e+00 * tc[0]
            -3.19328580e-03 * tc[1]
            +2.37674605e-05 * tc[2]
            -1.91528703e-08 * tc[3]
            +5.48277800e-12 * tc[4]
            +4.10301590e+00 ;
        /*species 27: CH2CO */
        species[27] =
            +2.13583630e+00 * tc[0]
            +1.81188721e-02 * tc[1]
            -8.69737370e-06 * tc[2]
            +3.11465856e-09 * tc[3]
            -5.03644037e-13 * tc[4]
            +1.22156480e+01 ;
        /*species 28: HCCO */
        species[28] =
            +2.25172140e+00 * tc[0]
            +1.76550210e-02 * tc[1]
            -1.18645505e-05 * tc[2]
            +5.75858633e-09 * tc[3]
            -1.26662028e-12 * tc[4]
            +1.24904170e+01 ;
        /*species 29: C2H */
        species[29] =
            +2.88965733e+00 * tc[0]
            +1.34099611e-02 * tc[1]
            -1.42384751e-05 * tc[2]
            +9.82636817e-09 * tc[3]
            -2.73328777e-12 * tc[4]
            +6.22296438e+00 ;
        /*species 30: CH2OH */
        species[30] =
            +4.47832317e+00 * tc[0]
            -1.35069687e-03 * tc[1]
            +1.39241853e-05 * tc[2]
            -1.21622466e-08 * tc[3]
            +3.69766938e-12 * tc[4]
            +3.30911984e+00 ;
        /*species 31: CH3OH */
        species[31] =
            +5.71539582e+00 * tc[0]
            -1.52309129e-02 * tc[1]
            +3.26220578e-05 * tc[2]
            -2.36935630e-08 * tc[3]
            +6.53381745e-12 * tc[4]
            -1.50409823e+00 ;
        /*species 32: C3H4 */
        species[32] =
            +2.61304450e+00 * tc[0]
            +1.21225750e-02 * tc[1]
            +9.26994000e-06 * tc[2]
            -1.15083830e-08 * tc[3]
            +3.83376975e-12 * tc[4]
            +1.02261390e+01 ;
        /*species 33: C3H3 */
        species[33] =
            +1.35110927e+00 * tc[0]
            +3.27411223e-02 * tc[1]
            -2.36913568e-05 * tc[2]
            +1.25436603e-08 * tc[3]
            -2.96352308e-12 * tc[4]
            +1.52058924e+01 ;
        /*species 34: C3H5 */
        species[34] =
            +1.36318350e+00 * tc[0]
            +1.98138210e-02 * tc[1]
            +6.24853000e-06 * tc[2]
            -1.11185183e-08 * tc[3]
            +3.96164275e-12 * tc[4]
            +1.71732140e+01 ;
        /*species 35: C3H6 */
        species[35] =
            +1.49330700e+00 * tc[0]
            +2.09251800e-02 * tc[1]
            +2.24339700e-06 * tc[2]
            -5.56304000e-09 * tc[3]
            +1.78953650e-12 * tc[4]
            +1.61453400e+01 ;
        /*species 36: C3H8 */
        species[36] =
            +9.28510930e-01 * tc[0]
            +2.64605660e-02 * tc[1]
            +3.01662230e-06 * tc[2]
            -7.30498433e-09 * tc[3]
            +2.37403860e-12 * tc[4]
            +1.92255380e+01 ;
        /*species 37: I-C3H7 */
        species[37] =
            +1.44491990e+00 * tc[0]
            +2.09991120e-02 * tc[1]
            +3.85181110e-06 * tc[2]
            -6.15875100e-09 * tc[3]
            +1.78207405e-12 * tc[4]
            +2.01163170e+01 ;
        /*species 38: N-C3H7 */
        species[38] =
            +1.04911730e+00 * tc[0]
            +2.60089730e-02 * tc[1]
            +1.17712580e-06 * tc[2]
            -6.53171067e-09 * tc[3]
            +2.34300518e-12 * tc[4]
            +2.11360340e+01 ;
    } else {
        /*species 0: N2 */
        species[0] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 1: AR */
        species[1] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600000e+00 ;
        /*species 2: H */
        species[2] =
            +2.50000001e+00 * tc[0]
            -2.30842973e-11 * tc[1]
            +8.07809740e-15 * tc[2]
            -1.57838412e-18 * tc[3]
            +1.24549339e-22 * tc[4]
            -4.46682914e-01 ;
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
            +2.86472886e+00 * tc[0]
            +1.05650448e-03 * tc[1]
            -1.29541379e-07 * tc[2]
            +1.01739558e-11 * tc[3]
            -3.32989690e-16 * tc[4]
            +5.70164073e+00 ;
        /*species 5: O */
        species[5] =
            +2.56942078e+00 * tc[0]
            -8.59741137e-05 * tc[1]
            +2.09742295e-08 * tc[2]
            -3.33925997e-12 * tc[3]
            +3.07084227e-16 * tc[4]
            +4.78433864e+00 ;
        /*species 6: H2 */
        species[6] =
            +3.33727920e+00 * tc[0]
            -4.94024731e-05 * tc[1]
            +2.49728389e-07 * tc[2]
            -5.98554647e-11 * tc[3]
            +5.00638440e-15 * tc[4]
            -3.20502331e+00 ;
        /*species 7: H2O */
        species[7] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
        /*species 8: HO2 */
        species[8] =
            +4.01721090e+00 * tc[0]
            +2.23982013e-03 * tc[1]
            -3.16829075e-07 * tc[2]
            +3.80821233e-11 * tc[3]
            -2.69771337e-15 * tc[4]
            +3.78510215e+00 ;
        /*species 9: H2O2 */
        species[9] =
            +4.16500285e+00 * tc[0]
            +4.90831694e-03 * tc[1]
            -9.50696125e-07 * tc[2]
            +1.23728662e-10 * tc[3]
            -7.19770763e-15 * tc[4]
            +2.91615662e+00 ;
        /*species 10: CO */
        species[10] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 11: CO2 */
        species[11] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
        /*species 12: HCO */
        species[12] =
            +2.77217438e+00 * tc[0]
            +4.95695526e-03 * tc[1]
            -1.24222806e-06 * tc[2]
            +1.96387259e-10 * tc[3]
            -1.33377178e-14 * tc[4]
            +9.79834492e+00 ;
        /*species 13: CH2O */
        species[13] =
            +1.76069008e+00 * tc[0]
            +9.20000082e-03 * tc[1]
            -2.21129406e-06 * tc[2]
            +3.35470707e-10 * tc[3]
            -2.20963910e-14 * tc[4]
            +1.36563230e+01 ;
        /*species 14: CH4 */
        species[14] =
            +7.48514950e-02 * tc[0]
            +1.33909467e-02 * tc[1]
            -2.86642905e-06 * tc[2]
            +4.07641783e-10 * tc[3]
            -2.54538075e-14 * tc[4]
            +1.84373180e+01 ;
        /*species 15: CH3 */
        species[15] =
            +2.28571772e+00 * tc[0]
            +7.23990037e-03 * tc[1]
            -1.49357174e-06 * tc[2]
            +1.98561548e-10 * tc[3]
            -1.16788599e-14 * tc[4]
            +8.48007179e+00 ;
        /*species 16: T-CH2 */
        species[16] =
            +2.87410113e+00 * tc[0]
            +3.65639292e-03 * tc[1]
            -7.04472985e-07 * tc[2]
            +8.67265163e-11 * tc[3]
            -4.69318918e-15 * tc[4]
            +6.17119324e+00 ;
        /*species 17: S-CH2 */
        species[17] =
            +2.29203842e+00 * tc[0]
            +4.65588637e-03 * tc[1]
            -1.00595973e-06 * tc[2]
            +1.39302000e-10 * tc[3]
            -8.49290912e-15 * tc[4]
            +8.62650169e+00 ;
        /*species 18: C2H4 */
        species[18] =
            +2.03611116e+00 * tc[0]
            +1.46454151e-02 * tc[1]
            -3.35538958e-06 * tc[2]
            +4.90743077e-10 * tc[3]
            -3.14265152e-14 * tc[4]
            +1.03053693e+01 ;
        /*species 19: CH3O */
        species[19] =
            +4.75779238e+00 * tc[0]
            +7.44142474e-03 * tc[1]
            -1.34852588e-06 * tc[2]
            +1.46030168e-10 * tc[3]
            -6.58842745e-15 * tc[4]
            -1.96680028e+00 ;
        /*species 20: C2H5 */
        species[20] =
            +1.95465642e+00 * tc[0]
            +1.73972722e-02 * tc[1]
            -3.99103334e-06 * tc[2]
            +5.84058963e-10 * tc[3]
            -3.74103940e-14 * tc[4]
            +1.34624343e+01 ;
        /*species 21: C2H6 */
        species[21] =
            +1.07188150e+00 * tc[0]
            +2.16852677e-02 * tc[1]
            -5.01280335e-06 * tc[2]
            +7.38040003e-10 * tc[3]
            -4.75007225e-14 * tc[4]
            +1.51156107e+01 ;
        /*species 22: CH */
        species[22] =
            +2.87846473e+00 * tc[0]
            +9.70913681e-04 * tc[1]
            +7.22228275e-08 * tc[2]
            -4.35626163e-11 * tc[3]
            +4.40198457e-15 * tc[4]
            +5.48497999e+00 ;
        /*species 23: C2H2 */
        species[23] =
            +4.14756964e+00 * tc[0]
            +5.96166664e-03 * tc[1]
            -1.18647426e-06 * tc[2]
            +1.55804057e-10 * tc[3]
            -9.03088033e-15 * tc[4]
            -1.23028121e+00 ;
        /*species 24: C2H3 */
        species[24] =
            +3.01672400e+00 * tc[0]
            +1.03302292e-02 * tc[1]
            -2.34041174e-06 * tc[2]
            +3.39210960e-10 * tc[3]
            -2.15651760e-14 * tc[4]
            +7.78732378e+00 ;
        /*species 25: CH2CHO */
        species[25] =
            +5.16620060e+00 * tc[0]
            +1.08478260e-02 * tc[1]
            -2.23291840e-06 * tc[2]
            +2.68761827e-10 * tc[3]
            -1.21025483e-14 * tc[4]
            -1.96333610e+00 ;
        /*species 26: C2H4O */
        species[26] =
            +5.40411080e+00 * tc[0]
            +1.17230590e-02 * tc[1]
            -2.11315685e-06 * tc[2]
            +2.27908170e-10 * tc[3]
            -1.02462158e-14 * tc[4]
            -3.48079170e+00 ;
        /*species 27: CH2CO */
        species[27] =
            +4.51129732e+00 * tc[0]
            +9.00359745e-03 * tc[1]
            -2.08469817e-06 * tc[2]
            +3.07781961e-10 * tc[3]
            -1.98709550e-14 * tc[4]
            +6.32247205e-01 ;
        /*species 28: HCCO */
        species[28] =
            +5.62820580e+00 * tc[0]
            +4.08534010e-03 * tc[1]
            -7.96727350e-07 * tc[2]
            +9.54201733e-11 * tc[3]
            -4.85195800e-15 * tc[4]
            -3.93025950e+00 ;
        /*species 29: C2H */
        species[29] =
            +3.16780652e+00 * tc[0]
            +4.75221902e-03 * tc[1]
            -9.18935385e-07 * tc[2]
            +1.01396751e-10 * tc[3]
            -4.43081925e-15 * tc[4]
            +6.63589475e+00 ;
        /*species 30: CH2OH */
        species[30] =
            +5.09312037e+00 * tc[0]
            +5.94758550e-03 * tc[1]
            -1.03248262e-06 * tc[2]
            +1.07668901e-10 * tc[3]
            -4.70312630e-15 * tc[4]
            -1.84690613e+00 ;
        /*species 31: CH3OH */
        species[31] =
            +1.78970791e+00 * tc[0]
            +1.40938292e-02 * tc[1]
            -3.18250418e-06 * tc[2]
            +4.60570283e-10 * tc[3]
            -2.92650550e-14 * tc[4]
            +1.45023623e+01 ;
        /*species 32: C3H4 */
        species[32] =
            +6.31687220e+00 * tc[0]
            +1.11337280e-02 * tc[1]
            -1.98146890e-06 * tc[2]
            +2.11880793e-10 * tc[3]
            -9.46888500e-15 * tc[4]
            -1.09957660e+01 ;
        /*species 33: C3H3 */
        species[33] =
            +7.14221880e+00 * tc[0]
            +7.61902005e-03 * tc[1]
            -1.33729975e-06 * tc[2]
            +1.41638267e-10 * tc[3]
            -6.28688537e-15 * tc[4]
            -1.25848436e+01 ;
        /*species 34: C3H5 */
        species[34] =
            +6.50078770e+00 * tc[0]
            +1.43247310e-02 * tc[1]
            -2.83908160e-06 * tc[2]
            +3.69360033e-10 * tc[3]
            -2.25909717e-14 * tc[4]
            -1.12430500e+01 ;
        /*species 35: C3H6 */
        species[35] =
            +6.73225700e+00 * tc[0]
            +1.49083400e-02 * tc[1]
            -2.47494950e-06 * tc[2]
            +2.40400733e-10 * tc[3]
            -9.41551000e-15 * tc[4]
            -1.33133500e+01 ;
        /*species 36: C3H8 */
        species[36] =
            +7.52441520e+00 * tc[0]
            +1.88982820e-02 * tc[1]
            -3.14605205e-06 * tc[2]
            +3.07204857e-10 * tc[3]
            -1.21711195e-14 * tc[4]
            -1.78383750e+01 ;
        /*species 37: I-C3H7 */
        species[37] =
            +6.51927410e+00 * tc[0]
            +1.72201040e-02 * tc[1]
            -2.86821085e-06 * tc[2]
            +2.80435773e-10 * tc[3]
            -1.11414782e-14 * tc[4]
            -9.08302150e+00 ;
        /*species 38: N-C3H7 */
        species[38] =
            +7.70974790e+00 * tc[0]
            +1.60314850e-02 * tc[1]
            -2.63601190e-06 * tc[2]
            +2.52961173e-10 * tc[3]
            -9.71567975e-15 * tc[4]
            -1.55152970e+01 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double * wt)
{
    wt[0] = 28.013400; /*N2 */
    wt[1] = 39.948000; /*AR */
    wt[2] = 1.007970; /*H */
    wt[3] = 31.998800; /*O2 */
    wt[4] = 17.007370; /*OH */
    wt[5] = 15.999400; /*O */
    wt[6] = 2.015940; /*H2 */
    wt[7] = 18.015340; /*H2O */
    wt[8] = 33.006770; /*HO2 */
    wt[9] = 34.014740; /*H2O2 */
    wt[10] = 28.010400; /*CO */
    wt[11] = 44.009800; /*CO2 */
    wt[12] = 29.018370; /*HCO */
    wt[13] = 30.026340; /*CH2O */
    wt[14] = 16.042880; /*CH4 */
    wt[15] = 15.034910; /*CH3 */
    wt[16] = 14.026940; /*T-CH2 */
    wt[17] = 14.026940; /*S-CH2 */
    wt[18] = 28.053880; /*C2H4 */
    wt[19] = 31.034310; /*CH3O */
    wt[20] = 29.061850; /*C2H5 */
    wt[21] = 30.069820; /*C2H6 */
    wt[22] = 13.018970; /*CH */
    wt[23] = 26.037940; /*C2H2 */
    wt[24] = 27.045910; /*C2H3 */
    wt[25] = 43.045310; /*CH2CHO */
    wt[26] = 44.053280; /*C2H4O */
    wt[27] = 42.037340; /*CH2CO */
    wt[28] = 41.029370; /*HCCO */
    wt[29] = 25.029970; /*C2H */
    wt[30] = 31.034310; /*CH2OH */
    wt[31] = 32.042280; /*CH3OH */
    wt[32] = 40.064880; /*C3H4 */
    wt[33] = 39.056910; /*C3H3 */
    wt[34] = 41.072850; /*C3H5 */
    wt[35] = 42.080820; /*C3H6 */
    wt[36] = 44.096760; /*C3H8 */
    wt[37] = 43.088790; /*I-C3H7 */
    wt[38] = 43.088790; /*N-C3H7 */

    return;
}


/*get temperature given internal energy in mass units and mass fracs */
int feeytt_(double * e, double * y, int * iwrk, double * rwrk, double * t)
{
    const int maxiter = 50;
    const double tol  = 0.001;
    double ein  = *e;
    double tmin = 300; // max lower bound for thermo def
    double tmax = 3000; // min upper bound for thermo def
    double e1,emin,emax,cv,t1,dt;
    int i; // loop counter
    CKUBMS(&tmin, y, iwrk, rwrk, &emin);
    CKUBMS(&tmax, y, iwrk, rwrk, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, iwrk, rwrk, &cv);
        *t = tmin - (emin-ein)/cv;
        return 1;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, iwrk, rwrk, &cv);
        *t = tmax - (emax-ein)/cv;
        return 1;
    }
    t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,iwrk,rwrk,&e1);
        CKCVBS(&t1,y,iwrk,rwrk,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100) { dt = 100; }
        else if (dt < -100) { dt = -100; }
        else if (fabs(dt) < tol) break;
        t1 += dt;
    }
    *t = t1;
    return 0;
}


/*convert phi[species] (specific mole nums) to y[species] (mass fracs) */
void fephity_(double * phi, int * iwrk, double * rwrk, double * y)
{
    double XW  = 0; 
    int id; /*loop counter */
    /*Compute mean molecular wt first */
    y[0] = phi[0]*28.013400;   XW += y[0]; /*N2 */
    y[1] = phi[1]*39.948000;   XW += y[1]; /*AR */
    y[2] = phi[2]*1.007970;   XW += y[2]; /*H */
    y[3] = phi[3]*31.998800;   XW += y[3]; /*O2 */
    y[4] = phi[4]*17.007370;   XW += y[4]; /*OH */
    y[5] = phi[5]*15.999400;   XW += y[5]; /*O */
    y[6] = phi[6]*2.015940;   XW += y[6]; /*H2 */
    y[7] = phi[7]*18.015340;   XW += y[7]; /*H2O */
    y[8] = phi[8]*33.006770;   XW += y[8]; /*HO2 */
    y[9] = phi[9]*34.014740;   XW += y[9]; /*H2O2 */
    y[10] = phi[10]*28.010400;   XW += y[10]; /*CO */
    y[11] = phi[11]*44.009800;   XW += y[11]; /*CO2 */
    y[12] = phi[12]*29.018370;   XW += y[12]; /*HCO */
    y[13] = phi[13]*30.026340;   XW += y[13]; /*CH2O */
    y[14] = phi[14]*16.042880;   XW += y[14]; /*CH4 */
    y[15] = phi[15]*15.034910;   XW += y[15]; /*CH3 */
    y[16] = phi[16]*14.026940;   XW += y[16]; /*T-CH2 */
    y[17] = phi[17]*14.026940;   XW += y[17]; /*S-CH2 */
    y[18] = phi[18]*28.053880;   XW += y[18]; /*C2H4 */
    y[19] = phi[19]*31.034310;   XW += y[19]; /*CH3O */
    y[20] = phi[20]*29.061850;   XW += y[20]; /*C2H5 */
    y[21] = phi[21]*30.069820;   XW += y[21]; /*C2H6 */
    y[22] = phi[22]*13.018970;   XW += y[22]; /*CH */
    y[23] = phi[23]*26.037940;   XW += y[23]; /*C2H2 */
    y[24] = phi[24]*27.045910;   XW += y[24]; /*C2H3 */
    y[25] = phi[25]*43.045310;   XW += y[25]; /*CH2CHO */
    y[26] = phi[26]*44.053280;   XW += y[26]; /*C2H4O */
    y[27] = phi[27]*42.037340;   XW += y[27]; /*CH2CO */
    y[28] = phi[28]*41.029370;   XW += y[28]; /*HCCO */
    y[29] = phi[29]*25.029970;   XW += y[29]; /*C2H */
    y[30] = phi[30]*31.034310;   XW += y[30]; /*CH2OH */
    y[31] = phi[31]*32.042280;   XW += y[31]; /*CH3OH */
    y[32] = phi[32]*40.064880;   XW += y[32]; /*C3H4 */
    y[33] = phi[33]*39.056910;   XW += y[33]; /*C3H3 */
    y[34] = phi[34]*41.072850;   XW += y[34]; /*C3H5 */
    y[35] = phi[35]*42.080820;   XW += y[35]; /*C3H6 */
    y[36] = phi[36]*44.096760;   XW += y[36]; /*C3H8 */
    y[37] = phi[37]*43.088790;   XW += y[37]; /*I-C3H7 */
    y[38] = phi[38]*43.088790;   XW += y[38]; /*N-C3H7 */
    for (id = 0; id < 39; ++id) {
        y[id] = y[id]/XW;
    }

    return;
}


/*convert y[species] (mass fracs) to phi[species] (specific mole num) */
void feytphi_(double * y, int * iwrk, double * rwrk, double * phi)
{
    phi[0] = y[0]/ 2.80134000e-02; /*N2 (wt in kg) */
    phi[1] = y[1]/ 3.99480000e-02; /*AR (wt in kg) */
    phi[2] = y[2]/ 1.00797000e-03; /*H (wt in kg) */
    phi[3] = y[3]/ 3.19988000e-02; /*O2 (wt in kg) */
    phi[4] = y[4]/ 1.70073700e-02; /*OH (wt in kg) */
    phi[5] = y[5]/ 1.59994000e-02; /*O (wt in kg) */
    phi[6] = y[6]/ 2.01594000e-03; /*H2 (wt in kg) */
    phi[7] = y[7]/ 1.80153400e-02; /*H2O (wt in kg) */
    phi[8] = y[8]/ 3.30067700e-02; /*HO2 (wt in kg) */
    phi[9] = y[9]/ 3.40147400e-02; /*H2O2 (wt in kg) */
    phi[10] = y[10]/ 2.80104000e-02; /*CO (wt in kg) */
    phi[11] = y[11]/ 4.40098000e-02; /*CO2 (wt in kg) */
    phi[12] = y[12]/ 2.90183700e-02; /*HCO (wt in kg) */
    phi[13] = y[13]/ 3.00263400e-02; /*CH2O (wt in kg) */
    phi[14] = y[14]/ 1.60428800e-02; /*CH4 (wt in kg) */
    phi[15] = y[15]/ 1.50349100e-02; /*CH3 (wt in kg) */
    phi[16] = y[16]/ 1.40269400e-02; /*T-CH2 (wt in kg) */
    phi[17] = y[17]/ 1.40269400e-02; /*S-CH2 (wt in kg) */
    phi[18] = y[18]/ 2.80538800e-02; /*C2H4 (wt in kg) */
    phi[19] = y[19]/ 3.10343100e-02; /*CH3O (wt in kg) */
    phi[20] = y[20]/ 2.90618500e-02; /*C2H5 (wt in kg) */
    phi[21] = y[21]/ 3.00698200e-02; /*C2H6 (wt in kg) */
    phi[22] = y[22]/ 1.30189700e-02; /*CH (wt in kg) */
    phi[23] = y[23]/ 2.60379400e-02; /*C2H2 (wt in kg) */
    phi[24] = y[24]/ 2.70459100e-02; /*C2H3 (wt in kg) */
    phi[25] = y[25]/ 4.30453100e-02; /*CH2CHO (wt in kg) */
    phi[26] = y[26]/ 4.40532800e-02; /*C2H4O (wt in kg) */
    phi[27] = y[27]/ 4.20373400e-02; /*CH2CO (wt in kg) */
    phi[28] = y[28]/ 4.10293700e-02; /*HCCO (wt in kg) */
    phi[29] = y[29]/ 2.50299700e-02; /*C2H (wt in kg) */
    phi[30] = y[30]/ 3.10343100e-02; /*CH2OH (wt in kg) */
    phi[31] = y[31]/ 3.20422800e-02; /*CH3OH (wt in kg) */
    phi[32] = y[32]/ 4.00648800e-02; /*C3H4 (wt in kg) */
    phi[33] = y[33]/ 3.90569100e-02; /*C3H3 (wt in kg) */
    phi[34] = y[34]/ 4.10728500e-02; /*C3H5 (wt in kg) */
    phi[35] = y[35]/ 4.20808200e-02; /*C3H6 (wt in kg) */
    phi[36] = y[36]/ 4.40967600e-02; /*C3H8 (wt in kg) */
    phi[37] = y[37]/ 4.30887900e-02; /*I-C3H7 (wt in kg) */
    phi[38] = y[38]/ 4.30887900e-02; /*N-C3H7 (wt in kg) */

    return;
}


/*reverse of ytcr, useful for rate computations */
void fectyr_(double * c, double * rho, int * iwrk, double * rwrk, double * y)
{
    y[0] = c[0] * 28.013400 / (*rho); 
    y[1] = c[1] * 39.948000 / (*rho); 
    y[2] = c[2] * 1.007970 / (*rho); 
    y[3] = c[3] * 31.998800 / (*rho); 
    y[4] = c[4] * 17.007370 / (*rho); 
    y[5] = c[5] * 15.999400 / (*rho); 
    y[6] = c[6] * 2.015940 / (*rho); 
    y[7] = c[7] * 18.015340 / (*rho); 
    y[8] = c[8] * 33.006770 / (*rho); 
    y[9] = c[9] * 34.014740 / (*rho); 
    y[10] = c[10] * 28.010400 / (*rho); 
    y[11] = c[11] * 44.009800 / (*rho); 
    y[12] = c[12] * 29.018370 / (*rho); 
    y[13] = c[13] * 30.026340 / (*rho); 
    y[14] = c[14] * 16.042880 / (*rho); 
    y[15] = c[15] * 15.034910 / (*rho); 
    y[16] = c[16] * 14.026940 / (*rho); 
    y[17] = c[17] * 14.026940 / (*rho); 
    y[18] = c[18] * 28.053880 / (*rho); 
    y[19] = c[19] * 31.034310 / (*rho); 
    y[20] = c[20] * 29.061850 / (*rho); 
    y[21] = c[21] * 30.069820 / (*rho); 
    y[22] = c[22] * 13.018970 / (*rho); 
    y[23] = c[23] * 26.037940 / (*rho); 
    y[24] = c[24] * 27.045910 / (*rho); 
    y[25] = c[25] * 43.045310 / (*rho); 
    y[26] = c[26] * 44.053280 / (*rho); 
    y[27] = c[27] * 42.037340 / (*rho); 
    y[28] = c[28] * 41.029370 / (*rho); 
    y[29] = c[29] * 25.029970 / (*rho); 
    y[30] = c[30] * 31.034310 / (*rho); 
    y[31] = c[31] * 32.042280 / (*rho); 
    y[32] = c[32] * 40.064880 / (*rho); 
    y[33] = c[33] * 39.056910 / (*rho); 
    y[34] = c[34] * 41.072850 / (*rho); 
    y[35] = c[35] * 42.080820 / (*rho); 
    y[36] = c[36] * 44.096760 / (*rho); 
    y[37] = c[37] * 43.088790 / (*rho); 
    y[38] = c[38] * 43.088790 / (*rho); 

    return;
}


/*ddebdf compatible right hand side of CV burner */
/*rwrk[0] and rwrk[1] should contain rho and ene respectively */
/*working variable phi contains specific mole numbers */
void fecvrhs_(double * time, double * phi, double * phidot, double * rwrk, int * iwrk)
{
    double rho,ene; /*CV Parameters */
    double y[39], wdot[39]; /*temporary storage */
    int i; /*Loop counter */
    double temperature,pressure; /*temporary var */
    rho = rwrk[0];
    ene = rwrk[1];
    fephity_(phi, iwrk, rwrk, y);
    feeytt_(&ene, y, iwrk, rwrk, &temperature);
    CKPY(&rho, &temperature,  y, iwrk, rwrk, &pressure);
    CKWYP(&pressure, &temperature,  y, iwrk, rwrk, wdot);
    for (i=0; i<39; ++i) phidot[i] = wdot[i] / (rho/1000.0); 

    return;
}


/*returns the dimensionality of the cv burner (number of species) */
int fecvdim_()
{
    return 39;
}


/*ddebdf compatible right hand side of ZND solver */
/*rwrk[0] : scaling factor for pressure */
/*rwrk[1] : preshock density (g/cc)  */
/*rwrk[2] : detonation velocity (cm/s)  */
/*solution vector: [P; rho; y0 ... ylast]  */
void fezndrhs_(double * time, double * z, double * zdot, double * rwrk, int * iwrk)
{
    double psc,rho1,udet; /*ZND Parameters */
    double wt[39], hms[39], wdot[39]; /*temporary storage */
    int i; /*Loop counter */
    /*temporary variables */
    double ru, T, uvel, wtm, p, rho, gam, son, xm, sum, drdy, eta, cp, cv ;
    double *y; /*mass frac pointer */

    ru = 8.314e+07;

    psc = rwrk[0];
    rho1 = rwrk[1];
    udet = rwrk[2];

    p = z[0] * psc;
    rho = z[1];

    y = &z[3];

    CKMMWY(y, 0, 0, &wtm);

    T = p * wtm / rho / ru;

    uvel = (rho1 * udet)/ rho;

    CKCPBS(&T, y, 0, 0, &cp);
    CKCVBS(&T, y, 0, 0, &cv);
    gam = cp/cv;

    son = sqrt(fabs(gam*ru*T/wtm));
    xm = uvel/son;

    CKHMS(&T, 0, 0, hms);
    CKWT(0, 0, wt);
    CKWYP(&p, &T, y, 0, 0, wdot);

    sum = 0.0;
    for (i=0; i<39; ++i) {
        zdot[i+3] = wdot[i] * wt[i] / rho;
        drdy = -rho * wtm / wt[i];
        sum += -( drdy + rho * hms[i]/ (cp*T) ) * zdot[i+3];
    }

    eta = 1.0 - xm*xm;
    zdot[0] = -(uvel*uvel/eta/psc)*sum;
    zdot[1] = -sum/eta;
    zdot[2] = uvel;

    return;
}


/*returns the dimensionality of the ZND solver (3+number of species) */
int feznddim_()
{
    return 42;
}


/*returns the name of the source mechanism file  */
char* femechfile_()
{
    return "";
}


/*returns the species number */
int fesymnum_(const char* s1)
{
    if (strcmp(s1, "N2")==0) return 0; 
    if (strcmp(s1, "AR")==0) return 1; 
    if (strcmp(s1, "H")==0) return 2; 
    if (strcmp(s1, "O2")==0) return 3; 
    if (strcmp(s1, "OH")==0) return 4; 
    if (strcmp(s1, "O")==0) return 5; 
    if (strcmp(s1, "H2")==0) return 6; 
    if (strcmp(s1, "H2O")==0) return 7; 
    if (strcmp(s1, "HO2")==0) return 8; 
    if (strcmp(s1, "H2O2")==0) return 9; 
    if (strcmp(s1, "CO")==0) return 10; 
    if (strcmp(s1, "CO2")==0) return 11; 
    if (strcmp(s1, "HCO")==0) return 12; 
    if (strcmp(s1, "CH2O")==0) return 13; 
    if (strcmp(s1, "CH4")==0) return 14; 
    if (strcmp(s1, "CH3")==0) return 15; 
    if (strcmp(s1, "T-CH2")==0) return 16; 
    if (strcmp(s1, "S-CH2")==0) return 17; 
    if (strcmp(s1, "C2H4")==0) return 18; 
    if (strcmp(s1, "CH3O")==0) return 19; 
    if (strcmp(s1, "C2H5")==0) return 20; 
    if (strcmp(s1, "C2H6")==0) return 21; 
    if (strcmp(s1, "CH")==0) return 22; 
    if (strcmp(s1, "C2H2")==0) return 23; 
    if (strcmp(s1, "C2H3")==0) return 24; 
    if (strcmp(s1, "CH2CHO")==0) return 25; 
    if (strcmp(s1, "C2H4O")==0) return 26; 
    if (strcmp(s1, "CH2CO")==0) return 27; 
    if (strcmp(s1, "HCCO")==0) return 28; 
    if (strcmp(s1, "C2H")==0) return 29; 
    if (strcmp(s1, "CH2OH")==0) return 30; 
    if (strcmp(s1, "CH3OH")==0) return 31; 
    if (strcmp(s1, "C3H4")==0) return 32; 
    if (strcmp(s1, "C3H3")==0) return 33; 
    if (strcmp(s1, "C3H5")==0) return 34; 
    if (strcmp(s1, "C3H6")==0) return 35; 
    if (strcmp(s1, "C3H8")==0) return 36; 
    if (strcmp(s1, "I-C3H7")==0) return 37; 
    if (strcmp(s1, "N-C3H7")==0) return 38; 
    /*species name not found */
    return -1;
}


/*returns the species name */
char* fesymname_(int sn)
{
    if (sn==0) return "N2"; 
    if (sn==1) return "AR"; 
    if (sn==2) return "H"; 
    if (sn==3) return "O2"; 
    if (sn==4) return "OH"; 
    if (sn==5) return "O"; 
    if (sn==6) return "H2"; 
    if (sn==7) return "H2O"; 
    if (sn==8) return "HO2"; 
    if (sn==9) return "H2O2"; 
    if (sn==10) return "CO"; 
    if (sn==11) return "CO2"; 
    if (sn==12) return "HCO"; 
    if (sn==13) return "CH2O"; 
    if (sn==14) return "CH4"; 
    if (sn==15) return "CH3"; 
    if (sn==16) return "T-CH2"; 
    if (sn==17) return "S-CH2"; 
    if (sn==18) return "C2H4"; 
    if (sn==19) return "CH3O"; 
    if (sn==20) return "C2H5"; 
    if (sn==21) return "C2H6"; 
    if (sn==22) return "CH"; 
    if (sn==23) return "C2H2"; 
    if (sn==24) return "C2H3"; 
    if (sn==25) return "CH2CHO"; 
    if (sn==26) return "C2H4O"; 
    if (sn==27) return "CH2CO"; 
    if (sn==28) return "HCCO"; 
    if (sn==29) return "C2H"; 
    if (sn==30) return "CH2OH"; 
    if (sn==31) return "CH3OH"; 
    if (sn==32) return "C3H4"; 
    if (sn==33) return "C3H3"; 
    if (sn==34) return "C3H5"; 
    if (sn==35) return "C3H6"; 
    if (sn==36) return "C3H8"; 
    if (sn==37) return "I-C3H7"; 
    if (sn==38) return "N-C3H7"; 
    /*species name not found */
    return "NOTFOUND";
}

/* End of file  */
