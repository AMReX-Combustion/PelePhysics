#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/* PURE CPU stuff  */
#ifndef AMREX_USE_GPU
namespace thermo
{
    amrex::Real fwd_A[3], fwd_beta[3], fwd_Ea[3];
    amrex::Real low_A[3], low_beta[3], low_Ea[3];
    amrex::Real rev_A[3], rev_beta[3], rev_Ea[3];
    amrex::Real troe_a[3],troe_Ts[3], troe_Tss[3], troe_Tsss[3];
    amrex::Real sri_a[3], sri_b[3], sri_c[3], sri_d[3], sri_e[3];
    amrex::Real activation_units[3], prefactor_units[3], phase_units[3];
    int is_PD[3], troe_len[3], sri_len[3], nTB[3], *TBid[3];
    amrex::Real *TB[3];
    std::vector<std::vector<amrex::Real>> kiv(3); 
    std::vector<std::vector<amrex::Real>> nuv(3); 
};

using namespace thermo;


/*Returns a count of species in a reaction, and their indices */
/*and stoichiometric coefficients. (Eq 50) */
void CKINU(int * i, int * nspec, int * ki, int * nu)
{
    if (*i < 1) {
        /*Return max num species per reaction */
        *nspec = 4;
    } else {
        if (*i > 3) {
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

/* Vectorized stuff  */

/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int *  np, amrex::Real *  y,  amrex::Real *  x)
{
    amrex::Real YOW[*np];
    amrex::Real imw[6];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<6; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<6; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}

/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tc[5], h[6];
    amrex::Real imw[6];

    get_imw(imw);

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
    }

    for (int n=0; n<6; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31446261815324e+07 * T[i] * imw[n];
        }
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int *  np, amrex::Real *  rho, amrex::Real *  T,
	    amrex::Real *  y,
	    amrex::Real *  wdot)
{
    amrex::Real c[6*(*np)]; /*temporary storage */
    amrex::Real imw[6];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<6; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<6*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW[*np];
    amrex::Real imw[6];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<6; n++) {
        for (int i=0; i<(*np); i++) {
            YOW[i] += y[n*(*np)+i] * imw[n];
        }
    }

    for (int i=0; i<(*np); i++) {
        P[i] = rho[i] * 8.31446261815324e+07 * T[i] * YOW[i]; /*P = rho*R*T/W */
    }

    return;
}

/*compute the production rate for each species */
void vproductionRate(int npt, amrex::Real *  wdot, amrex::Real *  sc, amrex::Real *  T)
{
    amrex::Real k_f_s[3*npt], Kc_s[3*npt], mixture[npt], g_RT[6*npt];
    amrex::Real tc[5*npt], invT[npt];

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

    for (int n=0; n<6; n++) {
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

void vcomp_k_f(int npt, amrex::Real *  k_f_s, amrex::Real *  tc, amrex::Real *  invT)
{
    for (int i=0; i<npt; i++) {
        k_f_s[0*npt+i] = prefactor_units[0] * fwd_A[0] * exp(fwd_beta[0] * tc[i] - activation_units[0] * fwd_Ea[0] * invT[i]);
        k_f_s[1*npt+i] = prefactor_units[1] * fwd_A[1] * exp(fwd_beta[1] * tc[i] - activation_units[1] * fwd_Ea[1] * invT[i]);
        k_f_s[2*npt+i] = prefactor_units[2] * fwd_A[2] * exp(fwd_beta[2] * tc[i] - activation_units[2] * fwd_Ea[2] * invT[i]);
    }
}

void vcomp_gibbs(int npt, amrex::Real *  g_RT, amrex::Real *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        amrex::Real tg[5], g[6];
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
    }
}

void vcomp_Kc(int npt, amrex::Real *  Kc_s, amrex::Real *  g_RT, amrex::Real *  invT)
{
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        amrex::Real refC = (101325. / 8.31451) * invT[i];
        amrex::Real refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refC * exp((3.000000 * g_RT[0*npt+i] + 2.000000 * g_RT[2*npt+i]) - (4.000000 * g_RT[1*npt+i] + 2.000000 * g_RT[3*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((g_RT[0*npt+i] + 2.000000 * g_RT[3*npt+i]) - (2.000000 * g_RT[4*npt+i]));
        Kc_s[2*npt+i] = refC * exp((2.000000 * g_RT[4*npt+i]) - (g_RT[0*npt+i] + 2.000000 * g_RT[3*npt+i]));
    }
}

void vcomp_wdot(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 1: 2.000000 CH4 + 3.000000 O2 => 2.000000 CO + 4.000000 H2O */
        phi_f = pow(sc[0*npt+i], 3.000000)*pow(sc[2*npt+i], 2.000000);
        k_f = k_f_s[0*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= 3.000000 * qdot;
        wdot[1*npt+i] += 4.000000 * qdot;
        wdot[2*npt+i] -= 2.000000 * qdot;
        wdot[3*npt+i] += 2.000000 * qdot;

        /*reaction 2: 2.000000 CO + O2 => 2.000000 CO2 */
        phi_f = sc[0*npt+i]*pow(sc[3*npt+i], 2.000000);
        k_f = k_f_s[1*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] -= 2.000000 * qdot;
        wdot[4*npt+i] += 2.000000 * qdot;

        /*reaction 3: 2.000000 CO2 => 2.000000 CO + O2 */
        phi_f = pow(sc[4*npt+i], 2.000000);
        k_f = k_f_s[2*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[3*npt+i] += 2.000000 * qdot;
        wdot[4*npt+i] -= 2.000000 * qdot;
    }
}

static amrex::Real T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static amrex::Real k_f_save[3];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static amrex::Real Kc_save[3];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif


/*compute the production rate for each species pointwise on CPU */
void productionRate(amrex::Real *  wdot, amrex::Real *  sc, amrex::Real T)
{
    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    amrex::Real qdot, q_f[3], q_r[3];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 6; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[0] -= 3.000000 * qdot;
    wdot[1] += 4.000000 * qdot;
    wdot[2] -= 2.000000 * qdot;
    wdot[3] += 2.000000 * qdot;

    qdot = q_f[1]-q_r[1];
    wdot[0] -= qdot;
    wdot[3] -= 2.000000 * qdot;
    wdot[4] += 2.000000 * qdot;

    qdot = q_f[2]-q_r[2];
    wdot[0] += qdot;
    wdot[3] += 2.000000 * qdot;
    wdot[4] -= 2.000000 * qdot;

    return;
}

void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f)
{
    for (int i=0; i<3; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(amrex::Real *  tc, amrex::Real invT, amrex::Real *  Kc)
{
    /*compute the Gibbs free energy */
    amrex::Real g_RT[6];
    gibbs(g_RT, tc);

    Kc[0] = 3.000000*g_RT[0] - 4.000000*g_RT[1] + 2.000000*g_RT[2] - 2.000000*g_RT[3];
    Kc[1] = g_RT[0] + 2.000000*g_RT[3] - 2.000000*g_RT[4];
    Kc[2] = -g_RT[0] - 2.000000*g_RT[3] + 2.000000*g_RT[4];

    for (int i=0; i<3; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 * invT;
    amrex::Real refCinv = 1 / refC;

    Kc[0] *= refC;
    Kc[1] *= refCinv;
    Kc[2] *= refC;

    return;
}

void comp_qfqr(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT)
{

    /*reaction 1: 2.000000 CH4 + 3.000000 O2 => 2.000000 CO + 4.000000 H2O */
    qf[0] = sc[0]*sc[2];
    qr[0] = 0.0;

    /*reaction 2: 2.000000 CO + O2 => 2.000000 CO2 */
    qf[1] = pow(sc[0], 0.250000)*pow(sc[1], 0.500000)*sc[3];
    qr[1] = 0.0;

    /*reaction 3: 2.000000 CO2 => 2.000000 CO + O2 */
    qf[2] = sc[4];
    qr[2] = 0.0;

    amrex::Real T = tc[1];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 6; ++i) {
        mixture += sc[i];
    }

    amrex::Real Corr[3];
    for (int i = 0; i < 3; ++i) {
        Corr[i] = 1.0;
    }

    for (int i=0; i<3; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the progress rate for each reaction */
void progressRate(amrex::Real *  qdot, amrex::Real *  sc, amrex::Real T)
{
    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    amrex::Real q_f[3], q_r[3];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 3; ++i) {
        qdot[i] = q_f[i] - q_r[i];
    }

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(amrex::Real *  q_f, amrex::Real *  q_r, amrex::Real *  sc, amrex::Real T)
{
    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    return;
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  q_f, amrex::Real *  q_r)
{
    int id; /*loop counter */
    amrex::Real c[6]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 6; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 6; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[6]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[6];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*O2 */
    YOW += y[1]*imw[1]; /*H2O */
    YOW += y[2]*imw[2]; /*CH4 */
    YOW += y[3]*imw[3]; /*CO */
    YOW += y[4]*imw[4]; /*CO2 */
    YOW += y[5]*imw[5]; /*N2 */
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

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[6]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 6; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[6]; /*temporary storage */
    amrex::Real imw[6];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[6]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*28.010550; /*CO */
    XW += x[4]*44.009950; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 6; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}

/*compute the reaction Jacobian on CPU */
void aJacobian(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP)
{
    for (int i=0; i<49; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[6];
    for (int k=0; k<6; k++) {
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
    for (int k = 0; k < 6; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[6];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[6];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[6];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: 2.000000 CH4 + 3.000000 O2 => 2.000000 CO + 4.000000 H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[0], 3.000000)*pow(sc[2], 2.000000);
    k_f = prefactor_units[0] * fwd_A[0]
                * exp(fwd_beta[0] * tc[0] - activation_units[0] * fwd_Ea[0] * invT);
    dlnkfdT = fwd_beta[0] * invT + activation_units[0] * fwd_Ea[0] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[0] -= 3 * q; /* O2 */
    wdot[1] += 4 * q; /* H2O */
    wdot[2] -= 2 * q; /* CH4 */
    wdot[3] += 2 * q; /* CO */
    /* d()/d[O2] */
    dqdci =  + k_f*3.000000*pow(sc[0],2.000000)*pow(sc[2], 2.000000);
    J[0] += -3 * dqdci;           /* dwdot[O2]/d[O2] */
    J[1] += 4 * dqdci;            /* dwdot[H2O]/d[O2] */
    J[2] += -2 * dqdci;           /* dwdot[CH4]/d[O2] */
    J[3] += 2 * dqdci;            /* dwdot[CO]/d[O2] */
    /* d()/d[CH4] */
    dqdci =  + k_f*pow(sc[0], 3.000000)*2.000000*sc[2];
    J[14] += -3 * dqdci;          /* dwdot[O2]/d[CH4] */
    J[15] += 4 * dqdci;           /* dwdot[H2O]/d[CH4] */
    J[16] += -2 * dqdci;          /* dwdot[CH4]/d[CH4] */
    J[17] += 2 * dqdci;           /* dwdot[CO]/d[CH4] */
    /* d()/dT */
    J[42] += -3 * dqdT;           /* dwdot[O2]/dT */
    J[43] += 4 * dqdT;            /* dwdot[H2O]/dT */
    J[44] += -2 * dqdT;           /* dwdot[CH4]/dT */
    J[45] += 2 * dqdT;            /* dwdot[CO]/dT */

    /*reaction 2: 2.000000 CO + O2 => 2.000000 CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*pow(sc[3], 2.000000);
    k_f = prefactor_units[1] * fwd_A[1]
                * exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
    dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[0] -= q; /* O2 */
    wdot[3] -= 2 * q; /* CO */
    wdot[4] += 2 * q; /* CO2 */
    /* d()/d[O2] */
    dqdci =  + k_f*pow(sc[3], 2.000000);
    J[0] -= dqdci;                /* dwdot[O2]/d[O2] */
    J[3] += -2 * dqdci;           /* dwdot[CO]/d[O2] */
    J[4] += 2 * dqdci;            /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[0]*2.000000*sc[3];
    J[21] -= dqdci;               /* dwdot[O2]/d[CO] */
    J[24] += -2 * dqdci;          /* dwdot[CO]/d[CO] */
    J[25] += 2 * dqdci;           /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[42] -= dqdT;                /* dwdot[O2]/dT */
    J[45] += -2 * dqdT;           /* dwdot[CO]/dT */
    J[46] += 2 * dqdT;            /* dwdot[CO2]/dT */

    /*reaction 3: 2.000000 CO2 => 2.000000 CO + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[0] += q; /* O2 */
    wdot[3] += 2 * q; /* CO */
    wdot[4] -= 2 * q; /* CO2 */
    /* d()/d[CO2] */
    dqdci =  + k_f*2.000000*sc[4];
    J[28] += dqdci;               /* dwdot[O2]/d[CO2] */
    J[31] += 2 * dqdci;           /* dwdot[CO]/d[CO2] */
    J[32] += -2 * dqdci;          /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[42] += dqdT;                /* dwdot[O2]/dT */
    J[45] += 2 * dqdT;            /* dwdot[CO]/dT */
    J[46] += -2 * dqdT;           /* dwdot[CO2]/dT */

    amrex::Real c_R[6], dcRdT[6], e_RT[6];
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
    for (int k = 0; k < 6; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[42+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 6; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 6; ++m) {
            dehmixdc += eh_RT[m]*J[k*7+m];
        }
        J[k*7+6] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[48] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

#endif
/* Transport function declarations  */


void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 24;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 846;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 6;}


void egtransetNLITE(int* NLITE ) {
    *NLITE = 0;}


/*Patm in ergs/cm3 */
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 3.19988000E+01;
    WT[1] = 1.80153400E+01;
    WT[2] = 1.60430300E+01;
    WT[3] = 2.80105500E+01;
    WT[4] = 4.40099500E+01;
    WT[5] = 2.80134000E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 1.07400000E+02;
    EPS[1] = 5.72400000E+02;
    EPS[2] = 1.41400000E+02;
    EPS[3] = 9.81000000E+01;
    EPS[4] = 2.44000000E+02;
    EPS[5] = 9.75300000E+01;
}


/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 3.45800000E+00;
    SIG[1] = 2.60500000E+00;
    SIG[2] = 3.74600000E+00;
    SIG[3] = 3.65000000E+00;
    SIG[4] = 3.76300000E+00;
    SIG[5] = 3.62100000E+00;
}


/*the dipole moment in Debye */
void egtransetDIP(amrex::Real* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 1.84400000E+00;
    DIP[2] = 0.00000000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 1.60000000E+00;
    POL[1] = 0.00000000E+00;
    POL[2] = 2.60000000E+00;
    POL[3] = 1.95000000E+00;
    POL[4] = 2.65000000E+00;
    POL[5] = 1.76000000E+00;
}


/*the rotational relaxation collision number at 298 K */
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 3.80000000E+00;
    ZROT[1] = 4.00000000E+00;
    ZROT[2] = 1.30000000E+01;
    ZROT[3] = 1.80000000E+00;
    ZROT[4] = 2.10000000E+00;
    ZROT[5] = 4.00000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 1;
    NLIN[1] = 2;
    NLIN[2] = 2;
    NLIN[3] = 1;
    NLIN[4] = 1;
    NLIN[5] = 1;
}


/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -1.68118868E+01;
    COFETA[1] = 2.52362554E+00;
    COFETA[2] = -2.49309128E-01;
    COFETA[3] = 1.10211025E-02;
    COFETA[4] = -1.17770937E+01;
    COFETA[5] = -8.26742721E-01;
    COFETA[6] = 3.39009079E-01;
    COFETA[7] = -2.00674327E-02;
    COFETA[8] = -1.95453421E+01;
    COFETA[9] = 3.36385478E+00;
    COFETA[10] = -3.56948469E-01;
    COFETA[11] = 1.56210922E-02;
    COFETA[12] = -1.63031240E+01;
    COFETA[13] = 2.26143219E+00;
    COFETA[14] = -2.15114671E-01;
    COFETA[15] = 9.53461976E-03;
    COFETA[16] = -2.36749526E+01;
    COFETA[17] = 4.99775518E+00;
    COFETA[18] = -5.52687718E-01;
    COFETA[19] = 2.34353338E-02;
    COFETA[20] = -1.62526779E+01;
    COFETA[21] = 2.24839597E+00;
    COFETA[22] = -2.13428438E-01;
    COFETA[23] = 9.46192413E-03;
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = -3.01284291E+00;
    COFLAM[1] = 3.37554994E+00;
    COFLAM[2] = -3.43353119E-01;
    COFLAM[3] = 1.51043444E-02;
    COFLAM[4] = 2.28195645E+01;
    COFLAM[5] = -8.72278946E+00;
    COFLAM[6] = 1.49300487E+00;
    COFLAM[7] = -7.41524047E-02;
    COFLAM[8] = 8.43903811E+00;
    COFLAM[9] = -2.78020298E+00;
    COFLAM[10] = 7.09313623E-01;
    COFLAM[11] = -4.04300978E-02;
    COFLAM[12] = 9.92459822E+00;
    COFLAM[13] = -2.28318157E+00;
    COFLAM[14] = 4.73113746E-01;
    COFLAM[15] = -2.40056652E-02;
    COFLAM[16] = -1.24047589E+01;
    COFLAM[17] = 6.34783131E+00;
    COFLAM[18] = -6.37857884E-01;
    COFLAM[19] = 2.37613820E-02;
    COFLAM[20] = 1.15507063E+01;
    COFLAM[21] = -2.91452379E+00;
    COFLAM[22] = 5.55043580E-01;
    COFLAM[23] = -2.75172461E-02;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.53110708E+01;
    COFD[1] = 3.37317428E+00;
    COFD[2] = -2.24900439E-01;
    COFD[3] = 9.81228151E-03;
    COFD[4] = -2.10640014E+01;
    COFD[5] = 5.50980695E+00;
    COFD[6] = -4.78335488E-01;
    COFD[7] = 1.98515434E-02;
    COFD[8] = -1.59878769E+01;
    COFD[9] = 3.66478157E+00;
    COFD[10] = -2.61506432E-01;
    COFD[11] = 1.13466231E-02;
    COFD[12] = -1.50371784E+01;
    COFD[13] = 3.26249588E+00;
    COFD[14] = -2.10658287E-01;
    COFD[15] = 9.20032462E-03;
    COFD[16] = -1.81197354E+01;
    COFD[17] = 4.33684042E+00;
    COFD[18] = -3.44981265E-01;
    COFD[19] = 1.48142449E-02;
    COFD[20] = -1.50096240E+01;
    COFD[21] = 3.25515933E+00;
    COFD[22] = -2.09710110E-01;
    COFD[23] = 9.15941830E-03;
    COFD[24] = -2.10640014E+01;
    COFD[25] = 5.50980695E+00;
    COFD[26] = -4.78335488E-01;
    COFD[27] = 1.98515434E-02;
    COFD[28] = -1.31492641E+01;
    COFD[29] = 1.48004311E+00;
    COFD[30] = 1.60499553E-01;
    COFD[31] = -1.19765679E-02;
    COFD[32] = -2.12953660E+01;
    COFD[33] = 5.52385268E+00;
    COFD[34] = -4.69683833E-01;
    COFD[35] = 1.90677265E-02;
    COFD[36] = -2.08943798E+01;
    COFD[37] = 5.44718652E+00;
    COFD[38] = -4.72082953E-01;
    COFD[39] = 1.96531321E-02;
    COFD[40] = -2.12021508E+01;
    COFD[41] = 5.20775052E+00;
    COFD[42] = -4.07348327E-01;
    COFD[43] = 1.55473283E-02;
    COFD[44] = -2.08123325E+01;
    COFD[45] = 5.42470154E+00;
    COFD[46] = -4.69700416E-01;
    COFD[47] = 1.95706904E-02;
    COFD[48] = -1.59878769E+01;
    COFD[49] = 3.66478157E+00;
    COFD[50] = -2.61506432E-01;
    COFD[51] = 1.13466231E-02;
    COFD[52] = -2.12953660E+01;
    COFD[53] = 5.52385268E+00;
    COFD[54] = -4.69683833E-01;
    COFD[55] = 1.90677265E-02;
    COFD[56] = -1.68532868E+01;
    COFD[57] = 4.00572564E+00;
    COFD[58] = -3.04255586E-01;
    COFD[59] = 1.31384083E-02;
    COFD[60] = -1.56953813E+01;
    COFD[61] = 3.54443207E+00;
    COFD[62] = -2.46133003E-01;
    COFD[63] = 1.06905018E-02;
    COFD[64] = -1.89594784E+01;
    COFD[65] = 4.68735220E+00;
    COFD[66] = -3.87449514E-01;
    COFD[67] = 1.65352261E-02;
    COFD[68] = -1.56756076E+01;
    COFD[69] = 3.54035770E+00;
    COFD[70] = -2.45653442E-01;
    COFD[71] = 1.06717969E-02;
    COFD[72] = -1.50371784E+01;
    COFD[73] = 3.26249588E+00;
    COFD[74] = -2.10658287E-01;
    COFD[75] = 9.20032462E-03;
    COFD[76] = -2.08943798E+01;
    COFD[77] = 5.44718652E+00;
    COFD[78] = -4.72082953E-01;
    COFD[79] = 1.96531321E-02;
    COFD[80] = -1.56953813E+01;
    COFD[81] = 3.54443207E+00;
    COFD[82] = -2.46133003E-01;
    COFD[83] = 1.06905018E-02;
    COFD[84] = -1.48061490E+01;
    COFD[85] = 3.16912473E+00;
    COFD[86] = -1.98792456E-01;
    COFD[87] = 8.69726395E-03;
    COFD[88] = -1.77673000E+01;
    COFD[89] = 4.20234040E+00;
    COFD[90] = -3.28057658E-01;
    COFD[91] = 1.41006192E-02;
    COFD[92] = -1.47850486E+01;
    COFD[93] = 3.16433919E+00;
    COFD[94] = -1.98191564E-01;
    COFD[95] = 8.67209742E-03;
    COFD[96] = -1.81197354E+01;
    COFD[97] = 4.33684042E+00;
    COFD[98] = -3.44981265E-01;
    COFD[99] = 1.48142449E-02;
    COFD[100] = -2.12021508E+01;
    COFD[101] = 5.20775052E+00;
    COFD[102] = -4.07348327E-01;
    COFD[103] = 1.55473283E-02;
    COFD[104] = -1.89594784E+01;
    COFD[105] = 4.68735220E+00;
    COFD[106] = -3.87449514E-01;
    COFD[107] = 1.65352261E-02;
    COFD[108] = -1.77673000E+01;
    COFD[109] = 4.20234040E+00;
    COFD[110] = -3.28057658E-01;
    COFD[111] = 1.41006192E-02;
    COFD[112] = -2.10907727E+01;
    COFD[113] = 5.29211327E+00;
    COFD[114] = -4.56068366E-01;
    COFD[115] = 1.91195062E-02;
    COFD[116] = -1.77350592E+01;
    COFD[117] = 4.19328271E+00;
    COFD[118] = -3.26911461E-01;
    COFD[119] = 1.40520357E-02;
    COFD[120] = -1.50096240E+01;
    COFD[121] = 3.25515933E+00;
    COFD[122] = -2.09710110E-01;
    COFD[123] = 9.15941830E-03;
    COFD[124] = -2.08123325E+01;
    COFD[125] = 5.42470154E+00;
    COFD[126] = -4.69700416E-01;
    COFD[127] = 1.95706904E-02;
    COFD[128] = -1.56756076E+01;
    COFD[129] = 3.54035770E+00;
    COFD[130] = -2.45653442E-01;
    COFD[131] = 1.06717969E-02;
    COFD[132] = -1.47850486E+01;
    COFD[133] = 3.16433919E+00;
    COFD[134] = -1.98191564E-01;
    COFD[135] = 8.67209742E-03;
    COFD[136] = -1.77350592E+01;
    COFD[137] = 4.19328271E+00;
    COFD[138] = -3.26911461E-01;
    COFD[139] = 1.40520357E-02;
    COFD[140] = -1.47639290E+01;
    COFD[141] = 3.15955654E+00;
    COFD[142] = -1.97590757E-01;
    COFD[143] = 8.64692156E-03;
}


/*List of specs with small weight, dim NLITE */
void egtransetKTDIF(int* KTDIF) {
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(amrex::Real* COFTD) {
}

#ifndef AMREX_USE_GPU
/* Initializes parameter database */
void CKINIT()
{

    // (0):  2.000000 CH4 + 3.000000 O2 => 2.000000 CO + 4.000000 H2O
    kiv[0] = {2,0,3,1};
    nuv[0] = {-2.0,-3.0,2.0,4.0};
    // (0):  2.000000 CH4 + 3.000000 O2 => 2.000000 CO + 4.000000 H2O
    fwd_A[0]     = 154500000000000;
    fwd_beta[0]  = 0.5;
    fwd_Ea[0]    = 39895;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = pow(10,-12.000000);
    is_PD[0] = 0;
    nTB[0] = 0;

    // (1):  2.000000 CO + O2 => 2.000000 CO2
    kiv[1] = {3,0,4};
    nuv[1] = {-2.0,-1,2.0};
    // (1):  2.000000 CO + O2 => 2.000000 CO2
    fwd_A[1]     = 199100000000000;
    fwd_beta[1]  = 0;
    fwd_Ea[1]    = 40000;
    prefactor_units[1]  = 3.1622776601683795e-05;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = pow(10,-10.500000);
    is_PD[1] = 0;
    nTB[1] = 0;

    // (2):  2.000000 CO2 => 2.000000 CO + O2
    kiv[2] = {4,3,0};
    nuv[2] = {-2.0,2.0,1};
    // (2):  2.000000 CO2 => 2.000000 CO + O2
    fwd_A[2]     = 250000000;
    fwd_beta[2]  = 0;
    fwd_Ea[2]    = 40000;
    prefactor_units[2]  = 1;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-6.000000);
    is_PD[2] = 0;
    nTB[2] = 0;

}


/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<3; ++i) {
    free(TB[i]); TB[i] = 0; 
    free(TBid[i]); TBid[i] = 0;
    nTB[i] = 0;
  }
}

#else
/* TODO: Remove on GPU, right now needed by chemistry_module on FORTRAN */
void CKINIT()
{
}

void CKFINALIZE()
{
}

#endif


/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 15.999400; /*O */
    awt[1] = 1.007970; /*H */
    awt[2] = 12.011150; /*C */
    awt[3] = 14.006700; /*N */

    return;
}



/*get atomic weight for all elements */
void CKAWT( amrex::Real *  awt)
{
    atomicWeight(awt);
}



/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * ncf)
{
    int id; /*loop counter */
    int kd = 4; 
    /*Zero ncf */
    for (id = 0; id < kd * 6; ++ id) {
         ncf[id] = 0; 
    }

    /*O2 */
    ncf[ 0 * kd + 0 ] = 2; /*O */

    /*H2O */
    ncf[ 1 * kd + 1 ] = 2; /*H */
    ncf[ 1 * kd + 0 ] = 1; /*O */

    /*CH4 */
    ncf[ 2 * kd + 2 ] = 1; /*C */
    ncf[ 2 * kd + 1 ] = 4; /*H */

    /*CO */
    ncf[ 3 * kd + 2 ] = 1; /*C */
    ncf[ 3 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 4 * kd + 2 ] = 1; /*C */
    ncf[ 4 * kd + 0 ] = 2; /*O */

    /*N2 */
    ncf[ 5 * kd + 3 ] = 2; /*N */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "O";
    ename[1] = "H";
    ename[2] = "C";
    ename[3] = "N";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(6);
    kname[0] = "O2";
    kname[1] = "H2O";
    kname[2] = "CH4";
    kname[3] = "CO";
    kname[4] = "CO2";
    kname[5] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(49);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(6);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[49];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<6; k++) {
                c_d[k] = 1.0/ 6.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<7; k++) {
        for (int l=0; l<7; l++) {
            if(J_h[ 7 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(49);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(6);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[49];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<6; k++) {
                c_d[k] = 1.0/ 6.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<7; k++) {
        for (int l=0; l<7; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 7 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}



/*compute the sparsity pattern of the simplified (for preconditioning) system Jacobian */
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, int * consP)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(49);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(6);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[49];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<6; k++) {
                c_d[k] = 1.0/ 6.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<7; k++) {
        for (int l=0; l<7; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 7 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    nJdata[0] = nJdata_tmp;

    return;
}


/*compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0 */
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, int * consP, int NCELLS)
{
    int offset_row;
    int offset_col;

    amrex::Gpu::DeviceVector<amrex::Real> J_v(49);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(6);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[49];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<6; k++) {
                c_d[k] = 1.0/ 6.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        offset_row = nc * 7;
        offset_col = nc * 7;
        for (int k=0; k<7; k++) {
            for (int l=0; l<7; l++) {
                if(J_h[7*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l + offset_row; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
            colPtrs[offset_col + (k + 1)] = nJdata_tmp;
        }
    }

    return;
}

/*compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0 */
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base)
{
    int offset;
    amrex::Gpu::DeviceVector<amrex::Real> J_v(49);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(6);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[49];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<6; k++) {
                c_d[k] = 1.0/ 6.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    if (base == 1) {
        rowPtrs[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 7;
            for (int l=0; l<7; l++) {
                for (int k=0; k<7; k++) {
                    if(J_h[7*k + l] != 0.0) {
                        colVals[nJdata_tmp-1] = k+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
                rowPtrs[offset + (l + 1)] = nJdata_tmp;
            }
        }
    } else {
        rowPtrs[0] = 0;
        int nJdata_tmp = 0;
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 7;
            for (int l=0; l<7; l++) {
                for (int k=0; k<7; k++) {
                    if(J_h[7*k + l] != 0.0) {
                        colVals[nJdata_tmp] = k + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
                rowPtrs[offset + (l + 1)] = nJdata_tmp;
            }
        }
    }

    return;
}

/*compute the sparsity pattern of the system Jacobian */
/*CSR format BASE is user choice */
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, int * consP, int NCELLS, int base)
{
    int offset;
    amrex::Gpu::DeviceVector<amrex::Real> J_v(49);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(6);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[49];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<6; k++) {
                c_d[k] = 1.0/ 6.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif
    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 7;
            for (int l=0; l<7; l++) {
                for (int k=0; k<7; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[7*k + l] != 0.0) {
                            colVals[nJdata_tmp-1] = k+1 + offset; 
                            nJdata_tmp = nJdata_tmp + 1; 
                        }
                    }
                }
                rowPtr[offset + (l + 1)] = nJdata_tmp;
            }
        }
    } else {
        rowPtr[0] = 0;
        int nJdata_tmp = 0;
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 7;
            for (int l=0; l<7; l++) {
                for (int k=0; k<7; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[7*k + l] != 0.0) {
                            colVals[nJdata_tmp] = k + offset; 
                            nJdata_tmp = nJdata_tmp + 1; 
                        }
                    }
                }
                rowPtr[offset + (l + 1)] = nJdata_tmp;
            }
        }
    }

    return;
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU */
/*BASE 0 */
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(49);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(6);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[49];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<6; k++) {
                c_d[k] = 1.0/ 6.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<7; k++) {
        for (int l=0; l<7; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 7*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[7*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 7*k + l;
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        colPtrs[k+1] = nJdata_tmp;
    }

    return;
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian */
/*CSR format BASE is under choice */
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(49);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(6);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[49];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<6; k++) {
                c_d[k] = 1.0/ 6.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<7; l++) {
            for (int k=0; k<7; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[7*k + l] != 0.0) {
                        colVals[nJdata_tmp-1] = k+1; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
            }
            rowPtr[l+1] = nJdata_tmp;
        }
    } else {
        rowPtr[0] = 0;
        int nJdata_tmp = 0;
        for (int l=0; l<7; l++) {
            for (int k=0; k<7; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[7*k + l] != 0.0) {
                        colVals[nJdata_tmp] = k; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
            }
            rowPtr[l+1] = nJdata_tmp;
        }
    }

    return;
}

#endif

