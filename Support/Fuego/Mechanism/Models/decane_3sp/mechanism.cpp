#include "chemistry_file.H"

#ifndef AMREX_USE_CUDA
namespace thermo
{
    double fwd_A[0], fwd_beta[0], fwd_Ea[0];
    double low_A[0], low_beta[0], low_Ea[0];
    double rev_A[0], rev_beta[0], rev_Ea[0];
    double troe_a[0],troe_Ts[0], troe_Tss[0], troe_Tsss[0];
    double sri_a[0], sri_b[0], sri_c[0], sri_d[0], sri_e[0];
    double activation_units[0], prefactor_units[0], phase_units[0];
    int is_PD[0], troe_len[0], sri_len[0], nTB[0], *TBid[0];
    double *TB[0];
    std::vector<std::vector<int>> kiv(0); 
    std::vector<std::vector<int>> nuv(0); 

    double fwd_A_DEF[0], fwd_beta_DEF[0], fwd_Ea_DEF[0];
    double low_A_DEF[0], low_beta_DEF[0], low_Ea_DEF[0];
    double rev_A_DEF[0], rev_beta_DEF[0], rev_Ea_DEF[0];
    double troe_a_DEF[0],troe_Ts_DEF[0], troe_Tss_DEF[0], troe_Tsss_DEF[0];
    double sri_a_DEF[0], sri_b_DEF[0], sri_c_DEF[0], sri_d_DEF[0], sri_e_DEF[0];
    double activation_units_DEF[0], prefactor_units_DEF[0], phase_units_DEF[0];
    int is_PD_DEF[0], troe_len_DEF[0], sri_len_DEF[0], nTB_DEF[0], *TBid_DEF[0];
    double *TB_DEF[0];
    std::vector<int> rxn_map;
};

using namespace thermo;
#endif

/* Inverse molecular weights */
/* TODO: check necessity on CPU */
static AMREX_GPU_DEVICE_MANAGED double imw[3] = {
    1.0 / 142.286840,  /*NC10H22 */
    1.0 / 31.998800,  /*O2 */
    1.0 / 28.013400};  /*N2 */

/* Inverse molecular weights */
/* TODO: check necessity because redundant with molecularWeight */
static AMREX_GPU_DEVICE_MANAGED double molecular_weights[3] = {
    142.286840,  /*NC10H22 */
    31.998800,  /*O2 */
    28.013400};  /*N2 */

AMREX_GPU_HOST_DEVICE
void get_imw(double imw_new[]){
    for(int i = 0; i<3; ++i) imw_new[i] = imw[i];
}

/* TODO: check necessity because redundant with CKWT */
AMREX_GPU_HOST_DEVICE
void get_mw(double mw_new[]){
    for(int i = 0; i<3; ++i) mw_new[i] = molecular_weights[i];
}


#ifndef AMREX_USE_CUDA
/* Initializes parameter database */
void CKINIT()
{

    rxn_map = {};

    SetAllDefaults();
}

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<0; ++i) {
        rmap[i] = rxn_map[i] + 1;
    }
}

#include <ReactionData.H>
double* GetParamPtr(int                reaction_id,
                    REACTION_PARAMETER param_id,
                    int                species_id,
                    int                get_default)
{
  printf("No reactions in this model");
  abort();
  return 0;
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

#else
/* TODO: Remove on GPU, right now needed by chemistry_module on FORTRAN */
AMREX_GPU_HOST_DEVICE void CKINIT()
{
}

AMREX_GPU_HOST_DEVICE void CKFINALIZE()
{
}

#endif


/*A few mechanism parameters */
void CKINDX(int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 4;
    *kk = 3;
    *ii = 0;
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
    for (i=0; i<lenkname*3; i++) {
        kname[i] = ' ';
    }

    /* NC10H22  */
    kname[ 0*lenkname + 0 ] = 'N';
    kname[ 0*lenkname + 1 ] = 'C';
    kname[ 0*lenkname + 2 ] = '1';
    kname[ 0*lenkname + 3 ] = '0';
    kname[ 0*lenkname + 4 ] = 'H';
    kname[ 0*lenkname + 5 ] = '2';
    kname[ 0*lenkname + 6 ] = '2';
    kname[ 0*lenkname + 7 ] = ' ';

    /* O2  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = '2';
    kname[ 1*lenkname + 2 ] = ' ';

    /* N2  */
    kname[ 2*lenkname + 0 ] = 'N';
    kname[ 2*lenkname + 1 ] = '2';
    kname[ 2*lenkname + 2 ] = ' ';

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
    XW += x[0]*142.286840; /*NC10H22 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
AMREX_GPU_HOST_DEVICE void CKPY(double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*NC10H22 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


#ifndef AMREX_USE_CUDA
/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<3; n++) {
        for (int i=0; i<(*np); i++) {
            YOW[i] += y[n*(*np)+i] * imw[n];
        }
    }

    for (int i=0; i<(*np); i++) {
        P[i] = rho[i] * 8.31451e+07 * T[i] * YOW[i]; /*P = rho*R*T/W */
    }

    return;
}
#endif


/*Compute P = rhoRT/W(c) */
void CKPC(double *  rho, double *  T, double *  c,  double *  P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*142.286840; /*NC10H22 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*28.013400; /*N2 */

    for (id = 0; id < 3; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double *  P, double *  T, double *  x,  double *  rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*142.286840; /*NC10H22 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*28.013400; /*N2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
AMREX_GPU_HOST_DEVICE void CKRHOY(double *  P, double *  T, double *  y,  double *  rho)
{
    double YOW = 0;
    double tmp[3];

    for (int i = 0; i < 3; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 3; i++)
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
    W += c[0]*142.286840; /*NC10H22 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*28.013400; /*N2 */

    for (id = 0; id < 3; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31451e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT( double *  wt)
{
    get_mw(wt);
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
    double tmp[3];

    for (int i = 0; i < 3; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 3; i++)
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
    XW += x[0]*142.286840; /*NC10H22 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*28.013400; /*N2 */
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
    W += c[0]*142.286840; /*NC10H22 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*28.013400; /*N2 */

    for (id = 0; id < 3; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
AMREX_GPU_HOST_DEVICE void CKYTX(double *  y,  double *  x)
{
    double YOW = 0;
    double tmp[3];

    for (int i = 0; i < 3; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 3; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 3; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


#ifndef AMREX_USE_CUDA
/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int *  np, double *  y,  double *  x)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<3; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<3; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}
#else
/*TODO: remove this on GPU */
void VCKYTX(int *  np, double *  y,  double *  x)
{
}
#endif


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double *  P, double *  T, double *  y,  double *  c)
{
    double YOW = 0;
    double PWORT;

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 3; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 3; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 3; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE void CKYTCR(double *  rho, double *  T, double *  y,  double *  c)
{
    for (int i = 0; i < 3; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE void CKXTY(double *  x,  double *  y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*142.286840; /*NC10H22 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*28.013400; /*N2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*142.286840*XWinv; 
    y[1] = x[1]*31.998800*XWinv; 
    y[2] = x[2]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double *  P, double *  T, double *  x,  double *  c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 3; ++id) {
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
    XW += x[0]*142.286840; /*NC10H22 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 3; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double *  c, double *  y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*142.286840; /*NC10H22 */
    CW += c[1]*31.998800; /*O2 */
    CW += c[2]*28.013400; /*N2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*142.286840*CWinv; 
    y[1] = c[1]*31.998800*CWinv; 
    y[2] = c[2]*28.013400*CWinv; 

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
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
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
    for (id = 0; id < 3; ++id) {
        sml[id] *= 8.31451e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
AMREX_GPU_HOST_DEVICE void CKCVMS(double *  T,  double *  cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 5.843484892910686e+05; /*NC10H22 */
    cvms[1] *= 2.598381814318037e+06; /*O2 */
    cvms[2] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_HOST_DEVICE void CKCPMS(double *  T,  double *  cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 5.843484892910686e+05; /*NC10H22 */
    cpms[1] *= 2.598381814318037e+06; /*O2 */
    cpms[2] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_HOST_DEVICE void CKUMS(double *  T,  double *  ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 3; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
AMREX_GPU_HOST_DEVICE void CKHMS(double *  T,  double *  hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 3; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


#ifndef AMREX_USE_CUDA
/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, double *  T,  double *  hms)
{
    double tc[5], h[3];

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
    }

    for (int n=0; n<3; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31451e+07 * T[i] * imw[n];
        }
    }
}
#else
/*TODO: remove this on GPU */
void VCKHMS(int *  np, double *  T,  double *  hms)
{
}
#endif


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double *  T,  double *  gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gms, tc);
    for (int i = 0; i < 3; i++)
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
    for (int i = 0; i < 3; i++)
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
    sms[0] *= 5.843484892910686e+05; /*NC10H22 */
    sms[1] *= 2.598381814318037e+06; /*O2 */
    sms[2] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *  T, double *  x,  double *  cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[3]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 3; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
AMREX_GPU_HOST_DEVICE void CKCPBS(double *  T, double *  y,  double *  cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[3], tresult[3]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 3; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 3; i++)
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
    double cvor[3]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 3; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
AMREX_GPU_HOST_DEVICE void CKCVBS(double *  T, double *  y,  double *  cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[3]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*NC10H22 */
    result += cvor[1]*y[1]*imw[1]; /*O2 */
    result += cvor[2]*y[2]*imw[2]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *  T, double *  x,  double *  hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[3]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 3; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
AMREX_GPU_HOST_DEVICE void CKHBMS(double *  T, double *  y,  double *  hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[3], tmp[3]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 3; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 3; ++id) {
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
    double uml[3]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 3; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
AMREX_GPU_HOST_DEVICE void CKUBMS(double *  T, double *  y,  double *  ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[3]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*NC10H22 */
    result += y[1]*ums[1]*imw[1]; /*O2 */
    result += y[2]*ums[2]*imw[2]; /*N2 */

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
    double sor[3]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 3; ++id) {
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
    double sor[3]; /* temporary storage */
    double x[3]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*NC10H22 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(142.286840*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(28.013400*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
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
    double gort[3]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 3; ++id) {
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
    double gort[3]; /* temporary storage */
    double x[3]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*NC10H22 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(142.286840*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(28.013400*YOW); 
    gibbs(gort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(gort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(gort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(gort[2]+log((x[2]+1e-100))+logPratio);
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
    double aort[3]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 3; ++id) {
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
    double aort[3]; /* temporary storage */
    double x[3]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*NC10H22 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(142.286840*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(28.013400*YOW); 
    helmholtz(aort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(aort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(aort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(aort[2]+log((x[2]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE void CKWC(double *  T, double *  C,  double *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 3; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double *  P, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*NC10H22 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double *  P, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 3; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
AMREX_GPU_HOST_DEVICE void CKWYR(double *  rho, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int *  np, double *  rho, double *  T,
	    double *  y,
	    double *  wdot)
{
#ifndef AMREX_USE_CUDA
    double c[3*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<3; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<3*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
#endif
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double *  rho, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*142.286840; /*NC10H22 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 3; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double *  T, double *  C, double *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 3; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 3; ++id) {
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
void CKQYP(double *  P, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*NC10H22 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double *  P, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 3; ++id) {
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
void CKQYR(double *  rho, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double *  rho, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[3]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*142.286840; /*NC10H22 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 3; ++id) {
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
void CKNU(int * kdim,  int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 3 * kd; ++ id) {
         nuki[id] = 0; 
    }
}


#ifndef AMREX_USE_CUDA
/*Returns a count of species in a reaction, and their indices */
/*and stoichiometric coefficients. (Eq 50) */
void CKINU(int * i, int * nspec, int * ki, int * nu)
{
    if (*i < 1) {
        /*Return max num species per reaction */
        *nspec = 0;
    } else {
        if (*i > 0) {
            *nspec = -1;
        } else {
            *nspec = kiv[*i-1].size();
            for (int j=0; j<*nspec; ++j) {
                ki[j] = kiv[*i-1][j] + 1;
                nu[j] = nuv[*i-1][j];
            }
        }
    }
}
#endif


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim,  int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 3; ++ id) {
         ncf[id] = 0; 
    }

    /*NC10H22 */
    ncf[ 0 * kd + 0 ] = 10; /*C */
    ncf[ 0 * kd + 1 ] = 22; /*H */

    /*O2 */
    ncf[ 1 * kd + 2 ] = 2; /*O */

    /*N2 */
    ncf[ 2 * kd + 3 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE( double *  a, double *  b, double *  e)
{

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double *  T, double *  C, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[3]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double *  P, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[3]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double *  P, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[3]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double *  rho, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[3]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double *  rho, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[3]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}

#ifdef AMREX_USE_CUDA
/*GPU version of productionRate: no more use of thermo namespace vectors */
/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE inline void  productionRate(double * wdot, double * sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    for (int i = 0; i < 3; ++i) {
        wdot[i] = 0.0;
    }

    return;
}

AMREX_GPU_HOST_DEVICE inline void comp_qfqr(double *  qf, double * qr, double * sc, double * tc, double invT)
{

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 3; ++i) {
        mixture += sc[i];
    }

    /*compute the Gibbs free energy */
    double g_RT[3];
    gibbs(g_RT, tc);

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;

    /* Evaluate the kfs */
    double k_f, Corr;



    return;
}
#endif


#ifndef AMREX_USE_CUDA
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


/*compute the production rate for each species pointwise on CPU */
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

    double qdot, q_f[0], q_r[0];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 3; ++i) {
        wdot[i] = 0.0;
    }

    return;
}

void comp_k_f(double *  tc, double invT, double *  k_f)
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

void comp_Kc(double *  tc, double invT, double *  Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[3];
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

void comp_qfqr(double *  qf, double *  qr, double *  sc, double *  tc, double invT)
{

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 3; ++i) {
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
#endif


#ifndef AMREX_USE_CUDA
/*compute the production rate for each species */
void vproductionRate(int npt, double *  wdot, double *  sc, double *  T)
{
    double k_f_s[0*npt], Kc_s[0*npt], mixture[npt], g_RT[3*npt];
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

    for (int n=0; n<3; n++) {
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
    }
}

void vcomp_gibbs(int npt, double *  g_RT, double *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[3];
        tg[0] = tc[0*npt+i];
        tg[1] = tc[1*npt+i];
        tg[2] = tc[2*npt+i];
        tg[3] = tc[3*npt+i];
        tg[4] = tc[4*npt+i];

        gibbs(g, tg);

        g_RT[0*npt+i] = g[0];
        g_RT[1*npt+i] = g[1];
        g_RT[2*npt+i] = g[2];
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
    }
}
#endif

/*compute an approx to the reaction Jacobian */
AMREX_GPU_HOST_DEVICE void DWDOT_PRECOND(double *  J, double *  sc, double *  Tp, int * HP)
{
    double c[3];

    for (int k=0; k<3; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<3; k++) {
        J[12+k] *= 1.e-6;
        J[k*4+3] *= 1.e6;
    }

    return;
}

/*compute an approx to the SPS Jacobian */
AMREX_GPU_HOST_DEVICE void SLJ_PRECOND_CSC(double *  Jsps, int * indx, int * len, double * sc, double * Tp, int * HP, double * gamma)
{
    double c[3];
    double J[16];
    double mwt[3];

    get_mw(mwt);

    for (int k=0; k<3; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* Change of coord */
    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<3; k++) {
        J[12+k] = 1.e-6 * J[12+k] * mwt[k];
        J[k*4+3] = 1.e6 * J[k*4+3] / mwt[k];
    }
    /* dTdot/dT */
    /* dwdot[l]/[k] */
    for (int k=0; k<3; k++) {
        for (int l=0; l<3; l++) {
            /* DIAG elem */
            if (k == l){
                J[ 4 * k + l] =  J[ 4 * k + l] * mwt[l] / mwt[k];
            /* NOT DIAG and not last column nor last row */
            } else {
                J[ 4 * k + l] =  J[ 4 * k + l] * mwt[l] / mwt[k];
            }
        }
    }

    for (int k=0; k<(*len); k++) {
        Jsps[k] = J[indx[k]];
    }

    return;
}

/*compute the reaction Jacobian */
AMREX_GPU_HOST_DEVICE void DWDOT(double *  J, double *  sc, double *  Tp, int * consP)
{
    double c[3];

    for (int k=0; k<3; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<3; k++) {
        J[12+k] *= 1.e-6;
        J[k*4+3] *= 1.e6;
    }

    return;
}

/*compute the sparsity pattern Jacobian */
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    double c[3];
    double J[16];

    for (int k=0; k<3; k++) {
        c[k] = 1.0/ 3.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<4; k++) {
        for (int l=0; l<4; l++) {
            if(J[ 4 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}



/*compute the sparsity pattern of simplified Jacobian */
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO_PRECOND( int * nJdata, int * consP)
{
    double c[3];
    double J[16];

    for (int k=0; k<3; k++) {
        c[k] = 1.0/ 3.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<4; k++) {
        for (int l=0; l<4; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J[ 4 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    nJdata[0] = nJdata_tmp;

    return;
}


#ifndef AMREX_USE_CUDA
/*compute the sparsity pattern of the simplified precond Jacobian on CPU */
void SPARSITY_PREPROC_PRECOND(int * rowVals, int * colPtrs, int * indx, int * consP)
{
    double c[3];
    double J[16];

    for (int k=0; k<3; k++) {
        c[k] = 1.0/ 3.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<4; k++) {
        for (int l=0; l<4; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 4*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J[4*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 4*k + l;
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        colPtrs[k+1] = nJdata_tmp;
    }

    return;
}
#else

/*compute the sparsity pattern of the simplified precond Jacobian on GPU */
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_PRECOND(int * rowPtr, int * colIndx, int * consP)
{
    double c[3];
    double J[16];

    for (int k=0; k<3; k++) {
        c[k] = 1.0/ 3.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l=0; l<4; l++) {
        for (int k=0; k<4; k++) {
            if (k == l) {
                colIndx[nJdata_tmp-1] = l+1; 
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J[4*k + l] != 0.0) {
                    colIndx[nJdata_tmp-1] = k+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        rowPtr[l+1] = nJdata_tmp;
    }

    return;
}
#endif

/*compute the sparsity pattern of the Jacobian */
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC(int *  rowVals, int *  colPtrs, int * consP, int NCELLS)
{
    double c[3];
    double J[16];
    int offset_row;
    int offset_col;

    for (int k=0; k<3; k++) {
        c[k] = 1.0/ 3.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        offset_row = nc * 4;
        offset_col = nc * 4;
        for (int k=0; k<4; k++) {
            for (int l=0; l<4; l++) {
                if(J[4*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l + offset_row; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
            colPtrs[offset_col + (k + 1)] = nJdata_tmp;
        }
    }

    return;
}


#ifdef AMREX_USE_CUDA
/*compute the reaction Jacobian on GPU */
AMREX_GPU_HOST_DEVICE
void aJacobian(double * J, double * sc, double T, int consP)
{


    for (int i=0; i<16; i++) {
        J[i] = 0.0;
    }

    double wdot[3];
    for (int k=0; k<3; k++) {
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
    for (int k = 0; k < 3; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[3];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[3];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[3];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    double c_R[3], dcRdT[3], e_RT[3];
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
    for (int k = 0; k < 3; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[12+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 3; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 3; ++m) {
            dehmixdc += eh_RT[m]*J[k*4+m];
        }
        J[k*4+3] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[15] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;

return;
}
#endif


#ifndef AMREX_USE_CUDA
/*compute the reaction Jacobian on CPU */
void aJacobian(double *  J, double *  sc, double T, int consP)
{
    for (int i=0; i<16; i++) {
        J[i] = 0.0;
    }

    double wdot[3];
    for (int k=0; k<3; k++) {
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
    for (int k = 0; k < 3; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[3];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[3];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[3];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    double c_R[3], dcRdT[3], e_RT[3];
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
    for (int k = 0; k < 3; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[12+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 3; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 3; ++m) {
            dehmixdc += eh_RT[m]*J[k*4+m];
        }
        J[k*4+3] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[15] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}
#endif


/*compute an approx to the reaction Jacobian */
AMREX_GPU_HOST_DEVICE void aJacobian_precond(double *  J, double *  sc, double T, int HP)
{
    for (int i=0; i<16; i++) {
        J[i] = 0.0;
    }

    double wdot[3];
    for (int k=0; k<3; k++) {
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
    for (int k = 0; k < 3; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[3];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[3];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[3];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    double c_R[3], dcRdT[3], e_RT[3];
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
    for (int k = 0; k < 3; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[12+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 3; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 3; ++m) {
            dehmixdc += eh_RT[m]*J[k*4+m];
        }
        J[k*4+3] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[15] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void dcvpRdT(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 2: N2 */
        species[2] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
    } else {
        /*species 1: O2 */
        species[1] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 2: N2 */
        species[2] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC10H22 */
        species[0] =
            +1.22535012e-01
            -1.55363148e-04 * tc[1]
            +7.49504631e-08 * tc[2]
            -1.29419215e-11 * tc[3];
    } else {
        /*species 0: NC10H22 */
        species[0] =
            +4.77244922e-02
            -3.24552782e-05 * tc[1]
            +7.52889777e-09 * tc[2]
            -5.80863088e-13 * tc[3];
    }
    return;
}


/*compute the progress rate for each reaction */
AMREX_GPU_HOST_DEVICE void progressRate(double *  qdot, double *  sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

#ifndef AMREX_USE_CUDA
    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }
#endif


    return;
}


/*compute the progress rate for each reaction */
AMREX_GPU_HOST_DEVICE void progressRateFR(double *  q_f, double *  q_r, double *  sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
#ifndef AMREX_USE_CUDA

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }
#endif

    comp_qfqr(q_f, q_r, sc, tc, invT);

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *  kc, double *  g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void gibbs(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 2: N2 */
        species[2] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
    } else {
        /*species 1: O2 */
        species[1] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 2: N2 */
        species[2] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC10H22 */
        species[0] =
            -3.430218630000000e+04 * invT
            -4.631018369000000e+01
            +2.084169690000000e+00 * tc[0]
            -6.126750600000000e-02 * tc[1]
            +1.294692898333333e-05 * tc[2]
            -2.081957308333334e-09 * tc[3]
            +1.617740190000000e-13 * tc[4];
    } else {
        /*species 0: NC10H22 */
        species[0] =
            -4.663928400000000e+04 * invT
            +1.724923449000000e+02
            -3.198822390000000e+01 * tc[0]
            -2.386224610000000e-02 * tc[1]
            +2.704606516666666e-06 * tc[2]
            -2.091360491666667e-10 * tc[3]
            +7.260788599999999e-15 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void helmholtz(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 2: N2 */
        species[2] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
    } else {
        /*species 1: O2 */
        species[1] =
            -1.08845772e+03 * invT
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 2: N2 */
        species[2] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC10H22 */
        species[0] =
            -3.43021863e+04 * invT
            -4.73101837e+01
            +2.08416969e+00 * tc[0]
            -6.12675060e-02 * tc[1]
            +1.29469290e-05 * tc[2]
            -2.08195731e-09 * tc[3]
            +1.61774019e-13 * tc[4];
    } else {
        /*species 0: NC10H22 */
        species[0] =
            -4.66392840e+04 * invT
            +1.71492345e+02
            -3.19882239e+01 * tc[0]
            -2.38622461e-02 * tc[1]
            +2.70460652e-06 * tc[2]
            -2.09136049e-10 * tc[3]
            +7.26078860e-15 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void cv_R(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 2: N2 */
        species[2] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
    } else {
        /*species 1: O2 */
        species[1] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 2: N2 */
        species[2] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC10H22 */
        species[0] =
            -3.08416969e+00
            +1.22535012e-01 * tc[1]
            -7.76815739e-05 * tc[2]
            +2.49834877e-08 * tc[3]
            -3.23548038e-12 * tc[4];
    } else {
        /*species 0: NC10H22 */
        species[0] =
            +3.09882239e+01
            +4.77244922e-02 * tc[1]
            -1.62276391e-05 * tc[2]
            +2.50963259e-09 * tc[3]
            -1.45215772e-13 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void cp_R(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 2: N2 */
        species[2] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
    } else {
        /*species 1: O2 */
        species[1] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 2: N2 */
        species[2] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC10H22 */
        species[0] =
            -2.08416969e+00
            +1.22535012e-01 * tc[1]
            -7.76815739e-05 * tc[2]
            +2.49834877e-08 * tc[3]
            -3.23548038e-12 * tc[4];
    } else {
        /*species 0: NC10H22 */
        species[0] =
            +3.19882239e+01
            +4.77244922e-02 * tc[1]
            -1.62276391e-05 * tc[2]
            +2.50963259e-09 * tc[3]
            -1.45215772e-13 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void speciesInternalEnergy(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 2: N2 */
        species[2] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
    } else {
        /*species 1: O2 */
        species[1] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 2: N2 */
        species[2] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC10H22 */
        species[0] =
            -3.08416969e+00
            +6.12675060e-02 * tc[1]
            -2.58938580e-05 * tc[2]
            +6.24587193e-09 * tc[3]
            -6.47096076e-13 * tc[4]
            -3.43021863e+04 * invT;
    } else {
        /*species 0: NC10H22 */
        species[0] =
            +3.09882239e+01
            +2.38622461e-02 * tc[1]
            -5.40921303e-06 * tc[2]
            +6.27408147e-10 * tc[3]
            -2.90431544e-14 * tc[4]
            -4.66392840e+04 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void speciesEnthalpy(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 2: N2 */
        species[2] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
    } else {
        /*species 1: O2 */
        species[1] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 2: N2 */
        species[2] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC10H22 */
        species[0] =
            -2.08416969e+00
            +6.12675060e-02 * tc[1]
            -2.58938580e-05 * tc[2]
            +6.24587193e-09 * tc[3]
            -6.47096076e-13 * tc[4]
            -3.43021863e+04 * invT;
    } else {
        /*species 0: NC10H22 */
        species[0] =
            +3.19882239e+01
            +2.38622461e-02 * tc[1]
            -5.40921303e-06 * tc[2]
            +6.27408147e-10 * tc[3]
            -2.90431544e-14 * tc[4]
            -4.66392840e+04 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void speciesEntropy(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            +3.78245636e+00 * tc[0]
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 2: N2 */
        species[2] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 1: O2 */
        species[1] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 2: N2 */
        species[2] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC10H22 */
        species[0] =
            -2.08416969e+00 * tc[0]
            +1.22535012e-01 * tc[1]
            -3.88407869e-05 * tc[2]
            +8.32782923e-09 * tc[3]
            -8.08870095e-13 * tc[4]
            +4.42260140e+01 ;
    } else {
        /*species 0: NC10H22 */
        species[0] =
            +3.19882239e+01 * tc[0]
            +4.77244922e-02 * tc[1]
            -8.11381955e-06 * tc[2]
            +8.36544197e-10 * tc[3]
            -3.63039430e-14 * tc[4]
            -1.40504121e+02 ;
    }
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
AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_EY(double *  e, double *  y, double *  t, int * ierr)
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
AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_HY(double *  h, double *  y, double *  t, int * ierr)
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

    double   EPS[3];
    double   SIG[3];
    double    wt[3];
    double avogadro = 6.02214199e23;
    double boltzmann = 1.3806503e-16; //we work in CGS
    double Rcst = 83.144598; //in bar [CGS] !

    egtransetEPS(EPS);
    egtransetSIG(SIG);
    get_mw(wt);

    /*species 0: NC10H22 */
    Tci[0] = 1.316 * EPS[0] ; 
    ai[0] = (5.55 * pow(avogadro,2.0) * EPS[0]*boltzmann * pow(1e-8*SIG[0],3.0) ) / (pow(wt[0],2.0)); 
    bi[0] = 0.855 * avogadro * pow(1e-8*SIG[0],3.0) / (wt[0]); 
    acentric_i[0] = 0.0 ;

    /*species 1: O2 */
    /*Imported from NIST */
    Tci[1] = 154.581000 ; 
    ai[1] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[1],2.0) / (pow(31.998800,2.0) * 50.430466); 
    bi[1] = 0.08664 * Rcst * Tci[1] / (31.998800 * 50.430466); 
    acentric_i[1] = 0.022200 ;

    /*species 2: N2 */
    /*Imported from NIST */
    Tci[2] = 126.192000 ; 
    ai[2] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[2],2.0) / (pow(28.013400,2.0) * 33.958000); 
    bi[2] = 0.08664 * Rcst * Tci[2] / (28.013400 * 33.958000); 
    acentric_i[2] = 0.037200 ;

    return;
}


void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 12;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 252;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 3;}


void egtransetNLITE(int* NLITE ) {
    *NLITE = 0;}


/*Patm in ergs/cm3 */
void egtransetPATM(double* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(double* WT ) {
    WT[0] = 1.42286840E+02;
    WT[1] = 3.19988000E+01;
    WT[2] = 2.80134000E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
void egtransetEPS(double* EPS ) {
    EPS[0] = 7.04917000E+02;
    EPS[1] = 1.07400000E+02;
    EPS[2] = 9.75300000E+01;
}


/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(double* SIG ) {
    SIG[0] = 6.67500000E+00;
    SIG[1] = 3.45800000E+00;
    SIG[2] = 3.62100000E+00;
}


/*the dipole moment in Debye */
void egtransetDIP(double* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[2] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
void egtransetPOL(double* POL ) {
    POL[0] = 0.00000000E+00;
    POL[1] = 1.60000000E+00;
    POL[2] = 1.76000000E+00;
}


/*the rotational relaxation collision number at 298 K */
void egtransetZROT(double* ZROT ) {
    ZROT[0] = 1.00000000E+00;
    ZROT[1] = 3.80000000E+00;
    ZROT[2] = 4.00000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 2;
    NLIN[1] = 1;
    NLIN[2] = 1;
}


/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(double* COFETA) {
    COFETA[0] = -1.29489842E+01;
    COFETA[1] = -6.04647944E-01;
    COFETA[2] = 3.07320970E-01;
    COFETA[3] = -1.89227139E-02;
    COFETA[4] = -1.68118868E+01;
    COFETA[5] = 2.52362554E+00;
    COFETA[6] = -2.49309128E-01;
    COFETA[7] = 1.10211025E-02;
    COFETA[8] = -1.62526779E+01;
    COFETA[9] = 2.24839597E+00;
    COFETA[10] = -2.13428438E-01;
    COFETA[11] = 9.46192413E-03;
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(double* COFLAM) {
    COFLAM[0] = -1.38317203E+01;
    COFLAM[1] = 5.32397257E+00;
    COFLAM[2] = -3.02248916E-01;
    COFLAM[3] = 1.51838934E-03;
    COFLAM[4] = -3.01284291E+00;
    COFLAM[5] = 3.37554994E+00;
    COFLAM[6] = -3.43353119E-01;
    COFLAM[7] = 1.51043444E-02;
    COFLAM[8] = 1.15507063E+01;
    COFLAM[9] = -2.91452379E+00;
    COFLAM[10] = 5.55043580E-01;
    COFLAM[11] = -2.75172461E-02;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(double* COFD) {
    COFD[0] = -1.65247199E+01;
    COFD[1] = 1.86376261E+00;
    COFD[2] = 9.20701178E-02;
    COFD[3] = -8.51470516E-03;
    COFD[4] = -2.24760016E+01;
    COFD[5] = 5.50822096E+00;
    COFD[6] = -4.78234659E-01;
    COFD[7] = 1.98512785E-02;
    COFD[8] = -2.21643696E+01;
    COFD[9] = 5.41977231E+00;
    COFD[10] = -4.69168499E-01;
    COFD[11] = 1.95518107E-02;
    COFD[12] = -2.24760016E+01;
    COFD[13] = 5.50822096E+00;
    COFD[14] = -4.78234659E-01;
    COFD[15] = 1.98512785E-02;
    COFD[16] = -1.53110708E+01;
    COFD[17] = 3.37317428E+00;
    COFD[18] = -2.24900439E-01;
    COFD[19] = 9.81228151E-03;
    COFD[20] = -1.50096240E+01;
    COFD[21] = 3.25515933E+00;
    COFD[22] = -2.09710110E-01;
    COFD[23] = 9.15941830E-03;
    COFD[24] = -2.21643696E+01;
    COFD[25] = 5.41977231E+00;
    COFD[26] = -4.69168499E-01;
    COFD[27] = 1.95518107E-02;
    COFD[28] = -1.50096240E+01;
    COFD[29] = 3.25515933E+00;
    COFD[30] = -2.09710110E-01;
    COFD[31] = 9.15941830E-03;
    COFD[32] = -1.47639290E+01;
    COFD[33] = 3.15955654E+00;
    COFD[34] = -1.97590757E-01;
    COFD[35] = 8.64692156E-03;
}


/*List of specs with small weight, dim NLITE */
void egtransetKTDIF(int* KTDIF) {
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(double* COFTD) {
}





#if 0




\\
\\
\\  This is the mechanism file
\\
\\
!*****************************************
! Dummy "mechanism" to represent 
! a simple 3species system for 
! n-decane non-reacting spray
! data culled from dodecane_lu mech 
!*****************************************


ELEMENTS
C H O N 
END

SPECIES
NC10H22 O2     N2
END

TRANS ALL
NC10H22            2   704.917     6.675     0.000     0.000     1.000          
O2                 1   107.400     3.458     0.000     1.600     3.800          
N2                 1    97.530     3.621     0.000     1.760     4.000          
END

REACTIONS
END

\\
\\
\\  This is the therm file
\\
\\
THERMO
   298.000  1000.000  5000.000
N2                121286N   2               G  0300.00   5000.00  1000.00      1
 0.02926640E+02 0.14879768E-02-0.05684760E-05 0.10097038E-09-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.14082404E-02-0.03963222E-04    3
 0.05641515E-07-0.02444854E-10-0.10208999E+04 0.03950372E+02                   4
AR                120186AR  1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.04366000E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366000E+02                   4
HE                L10/90HE  1    0    0    0G   200.000  6000.000 1000.        1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 9.28723974E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 9.28723974E-01 0.00000000E+00    4
O                 L 1/90O   1   00   00   00G   200.000  3500.000  1000.000    1
 2.56942078E+00-8.59741137E-05 4.19484589E-08-1.00177799E-11 1.22833691E-15    2
 2.92175791E+04 4.78433864E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00                   4
O2                TPIS89O   2   00   00   00G   200.000  3500.000  1000.000    1
 3.28253784E+00 1.48308754E-03-7.57966669E-07 2.09470555E-10-2.16717794E-14    2
-1.08845772E+03 5.45323129E+00 3.78245636E+00-2.99673416E-03 9.84730201E-06    3
-9.68129509E-09 3.24372837E-12-1.06394356E+03 3.65767573E+00                   4
H                 L 7/88H   1   00   00   00G   200.000  3500.000   1000.00    1
 2.50000001E+00-2.30842973E-11 1.61561948E-14-4.73515235E-18 4.98197357E-22    2
 2.54736599E+04-4.46682914E-01 2.50000000E+00 7.05332819E-13-1.99591964E-15    3
 2.30081632E-18-9.27732332E-22 2.54736599E+04-4.46682853E-01                   4
H2                TPIS78H   2   00   00   00G   200.000  3500.000   1000.00    1
 3.33727920E+00-4.94024731E-05 4.99456778E-07-1.79566394E-10 2.00255376E-14    2
-9.50158922E+02-3.20502331E+00 2.34433112E+00 7.98052075E-03-1.94781510E-05    3
 2.01572094E-08-7.37611761E-12-9.17935173E+02 6.83010238E-01                   4
OH                S 9/01O   1H   1    0    0G   200.000  6000.000   1000.00    1
 2.86472886E+00 1.05650448E-03-2.59082758E-07 3.05218674E-11-1.33195876E-15    2
 3.71885774E+03 5.70164073E+00 4.12530561E+00-3.22544939E-03 6.52764691E-06    3
-5.79853643E-09 2.06237379E-12 3.38153812E+03-6.90432960E-01                   4
H2O               L 8/89H   2O   1   00   00G   200.000  3500.000  1000.000    1
 3.03399249E+00 2.17691804E-03-1.64072518E-07-9.70419870E-11 1.68200992E-14    2
-3.00042971E+04 4.96677010E+00 4.19864056E+00-2.03643410E-03 6.52040211E-06    3
-5.48797062E-09 1.77197817E-12-3.02937267E+04-8.49032208E-01                   4
HO2               L 5/89H   1O   2   00   00G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00                   4
H2O2              L 7/88H   2O   2   00   00G   200.000  3500.000  1000.000    1
 4.16500285E+00 4.90831694E-03-1.90139225E-06 3.71185986E-10-2.87908305E-14    2
-1.78617877E+04 2.91615662E+00 4.27611269E+00-5.42822417E-04 1.67335701E-05    3
-2.15770813E-08 8.62454363E-12-1.77025821E+04 3.43505074E+00                   4
CH2               L S/93C   1H   2   00   00G   200.000  3500.000  1000.000    1
 2.87410113E+00 3.65639292E-03-1.40894597E-06 2.60179549E-10-1.87727567E-14    2
 4.62636040E+04 6.17119324E+00 3.76267867E+00 9.68872143E-04 2.79489841E-06    3
-3.85091153E-09 1.68741719E-12 4.60040401E+04 1.56253185E+00                   4
CH2*              L S/93C   1H   2   00   00G   200.000  3500.000  1000.000    1
 2.29203842E+00 4.65588637E-03-2.01191947E-06 4.17906000E-10-3.39716365E-14    2
 5.09259997E+04 8.62650169E+00 4.19860411E+00-2.36661419E-03 8.23296220E-06    3
-6.68815981E-09 1.94314737E-12 5.04968163E+04-7.69118967E-01                   4
CH3               L11/89C   1H   3   00   00G   200.000  3500.000  1000.000    1
 2.28571772E+00 7.23990037E-03-2.98714348E-06 5.95684644E-10-4.67154394E-14    2
 1.67755843E+04 8.48007179E+00 3.67359040E+00 2.01095175E-03 5.73021856E-06    3
-6.87117425E-09 2.54385734E-12 1.64449988E+04 1.60456433E+00                   4
CH4               L 8/88C   1H   4   00   00G   200.000  3500.000  1000.000    1
 7.48514950E-02 1.33909467E-02-5.73285809E-06 1.22292535E-09-1.01815230E-13    2
-9.46834459E+03 1.84373180E+01 5.14987613E+00-1.36709788E-02 4.91800599E-05    3
-4.84743026E-08 1.66693956E-11-1.02466476E+04-4.64130376E+00                   4
CO                TPIS79C   1O   1   00   00G   200.000  3500.000  1000.000    1
 2.71518561E+00 2.06252743E-03-9.98825771E-07 2.30053008E-10-2.03647716E-14    2
-1.41518724E+04 7.81868772E+00 3.57953347E+00-6.10353680E-04 1.01681433E-06    3
 9.07005884E-10-9.04424499E-13-1.43440860E+04 3.50840928E+00                   4
CO2               L 7/88C   1O   2   00   00G   200.000  3500.000  1000.000    1
 3.85746029E+00 4.41437026E-03-2.21481404E-06 5.23490188E-10-4.72084164E-14    2
-4.87591660E+04 2.27163806E+00 2.35677352E+00 8.98459677E-03-7.12356269E-06    3
 2.45919022E-09-1.43699548E-13-4.83719697E+04 9.90105222E+00                   4
HCO               L12/89H   1C   1O   1   00G   200.000  3500.000  1000.000    1
 2.77217438E+00 4.95695526E-03-2.48445613E-06 5.89161778E-10-5.33508711E-14    2
 4.01191815E+03 9.79834492E+00 4.22118584E+00-3.24392532E-03 1.37799446E-05    3
-1.33144093E-08 4.33768865E-12 3.83956496E+03 3.39437243E+00                   4
CH2O              L 8/88H   2C   1O   1   00G   200.000  3500.000  1000.000    1
 1.76069008E+00 9.20000082E-03-4.42258813E-06 1.00641212E-09-8.83855640E-14    2
-1.39958323E+04 1.36563230E+01 4.79372315E+00-9.90833369E-03 3.73220008E-05    3
-3.79285261E-08 1.31772652E-11-1.43089567E+04 6.02812900E-01                   4
CH2OH             IU2/03C   1H   3O   1   00G   200.000  6000.00               1
 5.09314370E+00 5.94761260E-03-2.06497460E-06 3.23008173E-10-1.88125902E-14    2
-4.03409640E+03-1.84691493E+00 4.47834367E+00-1.35070310E-03 2.78484980E-05    3
-3.64869060E-08 1.47907450E-11-3.50072890E+03 3.30913500E+00                   4
CH3O              IU1/03C   1H   3O   1     G   200.000  6000.00               1
 4.75779238E+00 7.44142474E-03-2.69705176E-06 4.38090504E-10-2.63537098E-14    2
 3.78111940E+02-1.96680028E+00 3.71180502E+00-2.80463306E-03 3.76550971E-05    3
-4.73072089E-08 1.86588420E-11 1.29569760E+03 6.57240864E+00                   4
CH3OH             L 8/88C   1H   4O   1   00G   200.000  3500.000  1000.000    1
 1.78970791E+00 1.40938292E-02-6.36500835E-06 1.38171085E-09-1.17060220E-13    2
-2.53748747E+04 1.45023623E+01 5.71539582E+00-1.52309129E-02 6.52441155E-05    3
-7.10806889E-08 2.61352698E-11-2.56427656E+04-1.50409823E+00                   4
C2H2              L 1/91C   2H   2   00   00G   200.000  3500.000  1000.000    1
 4.14756964E+00 5.96166664E-03-2.37294852E-06 4.67412171E-10-3.61235213E-14    2
 2.59359992E+04-1.23028121E+00 8.08681094E-01 2.33615629E-02-3.55171815E-05    3
 2.80152437E-08-8.50072974E-12 2.64289807E+04 1.39397051E+01                   4
C2H3              L 2/92C   2H   3   00   00G   200.000  3500.000  1000.000    1
 3.01672400E+00 1.03302292E-02-4.68082349E-06 1.01763288E-09-8.62607041E-14    2
 3.46128739E+04 7.78732378E+00 3.21246645E+00 1.51479162E-03 2.59209412E-05    3
-3.57657847E-08 1.47150873E-11 3.48598468E+04 8.51054025E+00                   4
C2H4              L 1/91C   2H   4   00   00G   200.000  3500.000  1000.000    1
 2.03611116E+00 1.46454151E-02-6.71077915E-06 1.47222923E-09-1.25706061E-13    2
 4.93988614E+03 1.03053693E+01 3.95920148E+00-7.57052247E-03 5.70990292E-05    3
-6.91588753E-08 2.69884373E-11 5.08977593E+03 4.09733096E+00                   4
C2H5              L12/92C   2H   5   00   00G   200.000  3500.000  1000.000    1
 1.95465642E+00 1.73972722E-02-7.98206668E-06 1.75217689E-09-1.49641576E-13    2
 1.28575200E+04 1.34624343E+01 4.30646568E+00-4.18658892E-03 4.97142807E-05    3
-5.99126606E-08 2.30509004E-11 1.28416265E+04 4.70720924E+00                   4
C2H6              L 8/88C   2H   6   00   00G   200.000  3500.000  1000.000    1
 1.07188150E+00 2.16852677E-02-1.00256067E-05 2.21412001E-09-1.90002890E-13    2
-1.14263932E+04 1.51156107E+01 4.29142492E+00-5.50154270E-03 5.99438288E-05    3
-7.08466285E-08 2.68685771E-11-1.15222055E+04 2.66682316E+00                   4
CH2CO             D05/90C   2H   2O   1   00G   200.000  3500.000  1000.000    1
 4.51129732E+00 9.00359745E-03-4.16939635E-06 9.23345882E-10-7.94838201E-14    2
-7.77850000E+03 6.32247205E-01 2.13583630E+00 1.81188721E-02-1.73947474E-05    3
 9.34397568E-09-2.01457615E-12-7.27000000E+03 1.22156480E+01                   4
CH2CHO            D05/83O   1H   3C   2    0G   300.000  5000.000              1
 0.59756699E+01 0.81305914E-02-0.27436245E-05 0.40703041E-09-0.21760171E-13    2
-0.96950000E+03-0.50320879E+01 0.34090624E+01 0.10738574E-01 0.18914925E-05    3
-0.71585831E-08 0.28673851E-11 0.62000000E+02 0.95714535E+01                   4
CH3CO             T 9/92C   2H   3O   1    0G   200.000  6000.0    1000.0      1
 0.59447731E+01 0.78667205E-02-0.28865882E-05 0.47270875E-09-0.28599861E-13    2
-0.37873075E+04-0.50136751E+01 0.41634257E+01-0.23261610E-03 0.34267820E-04    3
-0.44105227E-07 0.17275612E-10-0.26574529E+04 0.73468280E+01                   4
CH3CHO            L 8/88C   2H   4O   1    0G   200.000  6000.0    1000.0      1
 0.54041108E+01 0.11723059E-01-0.42263137E-05 0.68372451E-09-0.40984863E-13    2
-0.22593122E+05-0.34807917E+01 0.47294595E+01-0.31932858E-02 0.47534921E-04    3
-0.57458611E-07 0.21931112E-10-0.21572878E+05 0.41030159E+01                   4
HCCO              SRIC91H   1C   2O   1     G  0300.00   4000.00  1000.00      1
 0.56282058E+01 0.40853401E-02-0.15934547E-05 0.28626052E-09-0.19407832E-13    2
 0.19327215E+05-0.39302595E+01 0.22517214E+01 0.17655021E-01-0.23729101E-04    3
 0.17275759E-07-0.50664811E-11 0.20059449E+05 0.12490417E+02                   4
C3H3              T 5/97C   3H   3    0    0G   200.000  6000.000              1
 7.14221880E+00 7.61902005E-03-2.67459950E-06 4.24914801E-10-2.51475415E-14    2
 3.89087427E+04-1.25848436E+01 1.35110927E+00 3.27411223E-02-4.73827135E-05    3
 3.76309808E-08-1.18540923E-11 4.01057783E+04 1.52058924E+01 4.16139977E+04    4
aC3H4             L 8/89C   3H   4    0    0G   200.000  6000.000              1
 0.63168722E+01 0.11133728E-01-0.39629378E-05 0.63564238E-09-0.37875540E-13    2
 0.20117495E+05-0.10995766E+02 0.26130445E+01 0.12122575E-01 0.18539880E-04    3
-0.34525149E-07 0.15335079E-10 0.21541567E+05 0.10226139E+02 0.22962267E+05    4
nC3H7             P11/94C   3H   7    0    0G   300.000  3000.000              1
 0.77097479E+01 0.16031485E-01-0.52720238E-05 0.75888352E-09-0.38862719E-13    2
 0.79762236E+04-0.15515297E+02 0.10491173E+01 0.26008973E-01 0.23542516E-05    3
-0.19595132E-07 0.93720207E-11 0.10312346E+05 0.21136034E+02                   4
C3H6              120186C   3H   6          G  0300.00   5000.00  1000.00      1
 0.06732257E+02 0.01490834E+00-0.04949899E-04 0.07212022E-08-0.03766204E-12    2
-0.09235703E+04-0.01331335E+03 0.01493307E+02 0.02092518E+00 0.04486794E-04    3
-0.01668912E-06 0.07158146E-10 0.01074826E+05 0.01614534E+03                   4
C2H3CHO           USC/07C   3H   4O   1    0G   300.000  5000.000              1
 0.58111868E+01 0.17114256E-01-0.74834161E-05 0.14252249E-08-0.91746841E-13    2
-0.10784054E+05-0.48588004E+01 0.12713498E+01 0.26231054E-01-0.92912305E-05    3
-0.47837272E-08 0.33480543E-11-0.93357344E+04 0.19498077E+02                   4
aC3H5             PD5/98C   3H   5    0    0G   300.000  3000.000              1
 0.65007877E+01 0.14324731E-01-0.56781632E-05 0.11080801E-08-0.90363887E-13    2
 0.17482449E+05-0.11243050E+02 0.13631835E+01 0.19813821E-01 0.12497060E-04    3
-0.33355555E-07 0.15846571E-10 0.19245629E+05 0.17173214E+02                   4
iC4H5             USC/07C   4H   5O   0    0G   300.000  5000.000              1
 0.69646029E+01 0.18274333E-01-0.78133735E-05 0.15292154E-08-0.10920493E-12    2
 0.34725098E+05-0.10649321E+02 0.11308105E+00 0.40950615E-01-0.35413581E-04    3
 0.15530969E-07-0.23355122E-11 0.36383371E+05 0.23692457E+02                   4
C4H7              USC/07C   4H   7O   0    0G   300.000  5000.000              1
 0.70134835E+01 0.22634558E-01-0.92545470E-05 0.16807927E-08-0.10408617E-12    2
 0.20955008E+05-0.88893080E+01 0.74449432E+00 0.39678857E-01-0.22898086E-04    3
 0.21352973E-08 0.23096375E-11 0.22653328E+05 0.23437878E+02                   4
C4H6              H6W/94C   4H   6    0    0G   300.000  3000.000              1
 0.88673134E+01 0.14918670E-01-0.31548716E-05-0.41841330E-09 0.15761258E-12    2
 0.91338516E+04-0.23328171E+02 0.11284465E+00 0.34369022E-01-0.11107392E-04    3
-0.92106660E-08 0.62065179E-11 0.11802270E+05 0.23089996E+02                   4
pC4H9             USC/07C   4H   9O   0    0G   300.000  5000.000              1
 0.86822395E+01 0.23691071E-01-0.75948865E-05 0.66427136E-09 0.54845136E-13    2
 0.49644058E+04-0.17891747E+02 0.12087042E+01 0.38297497E-01-0.72660509E-05    3
-0.15428547E-07 0.86859435E-11 0.73221040E+04 0.22169268E+02                   4
C4H81             T 6/83C   4H   8    0    0G   300.000  5000.000              1
 0.20535841E+01 0.34350507E-01-0.15883197E-04 0.33089662E-08-0.25361045E-12    2
-0.21397231E+04 0.15543201E+02 0.11811380E+01 0.30853380E-01 0.50865247E-05    3
-0.24654888E-07 0.11110193E-10-0.17904004E+04 0.21062469E+02                   4
PXC5H11    1/ 2/ 7 THERMC   5H  11    0    0G   300.000  5000.000 1390.000    41
 1.52977446E+01 2.39735310E-02-8.18392948E-06 1.26883076E-09-7.35409055E-14    2
-9.80712307E+02-5.44829293E+01 5.24384081E-02 5.60796958E-02-3.31545803E-05    3
 9.77533781E-09-1.14009660E-12 4.71611460E+03 2.87238666E+01                   4
C5H9              000000N   0H   9O   0C   5G       300      5000    1000      1
 1.01386400E+01 2.27141380E-02-7.79104630E-06 1.18765220E-09-6.59324480E-14    2
-1.72183590E+03-3.31258850E+01-2.41901110E+00 4.04303890E-02 6.78023390E-06    3
-3.37247420E-08 1.51167130E-11 2.81218870E+03 3.64592440E+01                   4
C5H10      1/ 2/ 7 THERMC   5H  10    0    0G   300.000  5000.000 1392.000    31
 1.45851539E+01 2.24072471E-02-7.63348025E-06 1.18188966E-09-6.84385139E-14    2
-1.00898205E+04-5.23683936E+01-1.06223481E+00 5.74218294E-02-3.74486890E-05    3
 1.27364989E-08-1.79609789E-12-4.46546666E+03 3.22739790E+01                   4
PXC6H13    1/ 2/ 7 THERMC   6H  13    0    0G   300.000  5000.000 1390.000    51
 1.85385470E+01 2.83107962E-02-9.65307246E-06 1.49547585E-09-8.66336064E-14    2
-5.09299041E+03-7.04490943E+01-2.04871465E-01 6.83801272E-02-4.14447912E-05    3
 1.26155802E-08-1.53120058E-12 1.83280393E+03 3.16075093E+01                   4
C6H12      1/22/ 7 THERMC   6H  12    0    0G   300.000  5000.000 1392.000    41
 1.78337529E+01 2.67377658E-02-9.10036773E-06 1.40819768E-09-8.15124244E-14    2
-1.42062860E+04-6.83818851E+01-1.35275205E+00 6.98655426E-02-4.59408022E-05    3
 1.56967343E-08-2.21296175E-12-7.34368617E+03 3.53120691E+01                   4
PXC7H15    1/ 2/ 7 THERMC   7H  15    0    0G   300.000  5000.000 1390.000    61
 2.17940709E+01 3.26280243E-02-1.11138244E-05 1.72067148E-09-9.96366999E-14    2
-9.20938221E+03-8.64954311E+01-4.99570406E-01 8.08826467E-02-5.00532754E-05    3
 1.56549308E-08-1.96616227E-12-1.04590223E+03 3.46564011E+01                   4
C7H14      1/ 2/ 7 THERMC   7H  14    0    0G   300.000  5000.000 1392.000    51
 2.10898039E+01 3.10607878E-02-1.05644793E-05 1.63405780E-09-9.45598219E-14    2
-1.83260065E+04-8.44391108E+01-1.67720549E+00 8.24611601E-02-5.46504108E-05    3
 1.87862303E-08-2.65737983E-12-1.02168601E+04 3.85068032E+01                   4
NC8H18     1/ 2/ 7 THERMC   8H  18    0    0G   300.000  5000.000 1391.000    71
 2.54710194E+01 3.90887037E-02-1.33038777E-05 2.05867527E-09-1.19167174E-13    2
-3.83962755E+04-1.08361094E+02-1.54218406E+00 9.78112063E-02-6.09318358E-05    3
 1.92005591E-08-2.42996250E-12-2.85395641E+04 3.83327978E+01                   4
PXC8H17    1/ 2/ 7 THERMC   8H  17    0    0G   300.000  5000.000 1390.000    71
 2.50510356E+01 3.69480162E-02-1.25765264E-05 1.94628409E-09-1.12668898E-13    2
-1.33300535E+04-1.02557384E+02-7.72759438E-01 9.32549705E-02-5.84447245E-05    3
 1.85570214E-08-2.37127483E-12-3.92689511E+03 3.76130631E+01                   4
C8H16      1/ 2/ 7 THERMC   8H  16    0    0G   300.000  5000.000 1392.000    61
 2.43540125E+01 3.53666462E-02-1.20208388E-05 1.85855053E-09-1.07522262E-13    2
-2.24485674E+04-1.00537716E+02-1.89226915E+00 9.46066357E-02-6.27385521E-05    3
 2.15158309E-08-3.02718683E-12-1.31074559E+04 4.11878981E+01                   4
C8H15      1/ 2/ 7 THERMC   8H  15    0    0G   300.000  5000.000 1390.000    61
 2.41485380E+01 3.36466429E-02-1.15692934E-05 1.80261649E-09-1.04847473E-13    2
-5.51631151E+03-1.00894875E+02-1.90098561E+00 9.15067740E-02-6.01588113E-05    3
 2.02337556E-08-2.78235289E-12 3.87213558E+03 4.01405031E+01                   4
NC9H20     1/ 2/ 7 THERMC   9H  20    0    0G   300.000  5000.000 1391.000    81
 2.87289600E+01 4.34074576E-02-1.47660985E-05 2.28420987E-09-1.32194796E-13    2
-4.25174479E+04-1.24428751E+02-1.81390458E+00 1.10176644E-01-6.93124463E-05    3
 2.20957601E-08-2.83355715E-12-3.14207716E+04 4.12827220E+01                   4
PXC9H19    1/ 2/ 7 THERMC   9H  19    0    0G   300.000  5000.000 1390.000    81
 2.83097514E+01 4.12657344E-02-1.40383289E-05 2.17174871E-09-1.25692307E-13    2
-1.74516030E+04-1.16837897E+02-1.04387292E+00 1.05617283E-01-6.68199971E-05    3
 2.14486166E-08-2.77404275E-12-6.80818512E+03 4.23518992E+01                   4
C9H18      1/ 2/ 7 THERMC   9H  18    0    0G   300.000  5000.000 1392.000    71
 2.76142176E+01 3.96825287E-02-1.34819446E-05 2.08390452E-09-1.20539294E-13    2
-2.65709061E+04-1.16618623E+02-2.16108263E+00 1.06958297E-01-7.10973244E-05    3
 2.43971077E-08-3.42771547E-12-1.59890847E+04 4.41245128E+01                   4
NC10H22    1/ 2/ 7 THERMC  10H  22    0    0G   300.000  5000.000 1391.000    91
 3.19882239E+01 4.77244922E-02-1.62276391E-05 2.50963259E-09-1.45215772E-13    2
-4.66392840E+04-1.40504121E+02-2.08416969E+00 1.22535012E-01-7.76815739E-05    3
 2.49834877E-08-3.23548038E-12-3.43021863E+04 4.42260140E+01                   4
PXC10H21   1/ 2/ 7 THERMC  10H  21    0    0G   300.000  5000.000 1390.000    91
 3.15697160E+01 4.55818403E-02-1.54994965E-05 2.39710933E-09-1.38709559E-13    2
-2.15737832E+04-1.34708986E+02-1.31358348E+00 1.17972813E-01-7.51843079E-05    3
 2.43331106E-08-3.17522852E-12-9.68967550E+03 4.35010452E+01                   4
C10H20     1/22/ 7 THERMC  10H  20    0    0G   300.000  5000.000 1392.000    81
 3.08753903E+01 4.39971526E-02-1.49425530E-05 2.30917678E-09-1.33551477E-13    2
-3.06937307E+04-1.32705172E+02-2.42901688E+00 1.19305598E-01-7.94489025E-05    3
 2.72736596E-08-3.82718373E-12-1.88708365E+04 4.70571383E+01                   4
NC11H24    1/22/ 7 THERMC  11H  24    0    0G   300.000  5000.000 1391.000    11
 3.52484813E+01 5.20402416E-02-1.76886732E-05 2.73497226E-09-1.58231832E-13    2
-5.07616214E+04-1.56585288E+02-2.35338447E+00 1.34888270E-01-8.60424000E-05    3
 2.78658195E-08-3.63619953E-12-3.71837502E+04 4.71645217E+01                   4
PXC11H23   1/22/ 7 THERMC  11H  23    0    0G   300.000  5000.000 1391.000    11
 3.48306023E+01 4.98967610E-02-1.69601993E-05 2.62239406E-09-1.51722334E-13    2
-2.56964317E+04-1.50793815E+02-1.58230042E+00 1.30323530E-01-8.35408417E-05    3
 2.72125730E-08-3.57529580E-12-1.25713075E+04 4.64373071E+01                   4
C11H22     1/ 2/ 7 THERMC  11H  22    0    0G   300.000  5000.000 1392.000    91
 3.41376800E+01 4.83102359E-02-1.64025319E-05 2.53434308E-09-1.46557255E-13    2
-3.48170932E+04-1.48798197E+02-2.69653994E+00 1.31650370E-01-8.77957172E-05    3
 3.01467965E-08-4.22584486E-12-2.17526343E+04 4.99879917E+01                   4
NC12H26    1/ 2/ 7 THERMC  12H  26    0    0G   300.000  5000.000 1391.000    11
 3.85095037E+01 5.63550048E-02-1.91493200E-05 2.96024862E-09-1.71244150E-13    2
-5.48843465E+04-1.72670922E+02-2.62181594E+00 1.47237711E-01-9.43970271E-05    3
 3.07441268E-08-4.03602230E-12-4.00654253E+04 5.00994626E+01                   4
PXC12H25   1/ 2/ 7 THERMC  12H  25    0    0G   300.000  5000.000 1390.000    11
 3.80921885E+01 5.42107848E-02-1.84205517E-05 2.84762173E-09-1.64731748E-13    2
-2.98194375E+04-1.66882734E+02-1.85028741E+00 1.42670708E-01-9.18916555E-05    3
 3.00883392E-08-3.97454300E-12-1.54530435E+04 4.93702421E+01                   4
SXC12H25   1/ 2/ 7 THERMC  12H  25    0    0G   300.000  5000.000 1385.000    11
 3.79688268E+01 5.38719464E-02-1.82171263E-05 2.80774503E-09-1.62108420E-13    2
-3.12144988E+04-1.65805933E+02-1.36787089E+00 1.37355348E-01-8.24076158E-05    3
 2.36421562E-08-2.47435932E-12-1.67660539E+04 4.83521895E+01                   4
S3XC12H25  1/ 2/ 7 THERMC  12H  25    0    0G   300.000  5000.000 1385.000    11
 3.79688268E+01 5.38719464E-02-1.82171263E-05 2.80774503E-09-1.62108420E-13    2
-3.12144988E+04-1.65805933E+02-1.36787089E+00 1.37355348E-01-8.24076158E-05    3
 2.36421562E-08-2.47435932E-12-1.67660539E+04 4.83521895E+01                   4
C12H24     1/22/ 7 THERMC  12H  24    0    0G   300.000  5000.000 1391.000    11
 3.74002111E+01 5.26230753E-02-1.78624319E-05 2.75949863E-09-1.59562499E-13    2
-3.89405962E+04-1.64892663E+02-2.96342681E+00 1.43992360E-01-9.61384015E-05    3
 3.30174473E-08-4.62398190E-12-2.46345299E+04 5.29158870E+01                   4
C12H26O2   1/24/ 7 THERMC  12H  26O   2    0G   300.000  5000.000 1390.000    31
 4.35257442E+01 5.78016040E-02-1.99632594E-05 3.11946519E-09-1.81796787E-13    2
-6.68122451E+04-1.92119235E+02-5.33975911E-01 1.56648083E-01-1.05028971E-04    3
 3.67984049E-08-5.38662910E-12-5.09380146E+04 4.62325087E+01                   4
C12H25O2          000000N   0H  25O   2C  12G       300      5000    1000      1
 2.84782000E+01 5.37539000E-02-1.68186000E-05 2.51367000E-09-1.47208000E-13    2
-3.74118000E+04-1.09121000E+02 5.31404000E+00 8.93873000E-02 1.45351000E-05    3
-7.49250000E-08 3.35325000E-11-2.98918000E+04 1.69741000E+01                   4

!C12H25O2   1/24/ 7 THERMC  12H  25O   2    0G   300.000  5000.000 1389.000    21
! 4.10720559E+01 5.76223879E-02-1.98682057E-05 3.10098353E-09-1.80567383E-13    2
!-4.84961160E+04-1.77176648E+02 4.13301762E-01 1.47893311E-01-9.69755376E-05    3
! 3.35724323E-08-4.91454646E-12-3.36844450E+04 4.32258761E+01                   4
C12OOH            000000N   0H  25O   2C  12G       300      5000    1000      1
 2.92019000E+01 5.15917000E-02-1.57327000E-05 2.30306000E-09-1.32640000E-13    2
-3.11192000E+04-1.08855000E+02 5.15231000E+00 9.97913000E-02-1.80635000E-05    3
-4.18435000E-08 2.22786000E-11-2.38380000E+04 1.93526000E+01                   4

!C12OOH     1/24/ 7 THERMC  12H  25O   2    0G   300.000  5000.000 1387.000    31
! 4.29994696E+01 5.50833554E-02-1.88955870E-05 2.93978363E-09-1.70823564E-13    2
!-4.30658136E+04-1.85856878E+02 1.60721186E-01 1.49325879E-01-9.67266789E-05    3
! 3.17154092E-08-4.20449417E-12-2.75569298E+04 4.63656051E+01                   4
C12H26O4   1/24/ 7 THERMC  12H  26O   4    0G   300.000  5000.000 1390.000    51
 5.03368769E+01 5.67494383E-02-1.97115026E-05 3.09203723E-09-1.80689139E-13    2
-8.14466175E+04-2.23085533E+02-2.94753190E-01 1.75805891E-01-1.28643781E-04    3
 4.93995330E-08-7.88291511E-12-6.37896715E+04 4.88174790E+01                   4
O2C12H24OOH       000000N   0H  25O   4C  12G       300      5000    1000      1
 3.50907000E+01 5.10590000E-02-1.54345000E-05 2.24627000E-09-1.28901000E-13    2
-5.12675000E+04-1.37750000E+02 4.81972000E-01 1.45020000E-01-9.99308000E-05    3
 2.60422000E-08 1.19358000E-12-4.16875000E+04 4.13429000E+01                   4

!O2C12H24OOH        THERMC  12H  25O   4    0G   300.000  5000.000 1389.000    41
! 4.78237378E+01 5.66012764E-02-1.96226767E-05 3.07404604E-09-1.79468346E-13    2
!-6.30929400E+04-2.07768890E+02 6.91692611E-01 1.67000994E-01-1.20729787E-04    3
! 4.63412670E-08-7.45842745E-12-4.65446291E+04 4.56074372E+01                   4
OC12H23OOH        000000N   0H  24O   3C  12G       300      5000    1000      1
 2.36731000E+01 6.16392000E-02-2.09836000E-05 3.33166000E-09-2.03590000E-13    2
-7.18258000E+04-7.77662000E+01 8.80733000E+00 6.50623000E-02 6.95058000E-05    3
-1.26905000E-07 5.10991000E-11-6.65361000E+04 6.84155000E+00                   4

!OC12H23OOH         THERMC  12H  24O   3    0G   300.000  5000.000 1392.000    21
! 4.49664435E+01 5.45029400E-02-1.88364401E-05 2.94475639E-09-1.71672676E-13    2
!-7.80784935E+04-1.96879294E+02-9.04780358E-01 1.62600770E-01-1.18036891E-04    3
! 4.52507973E-08-7.22772167E-12-6.21001449E+04 4.93825152E+01                   4
ENDOFDATA

#endif
