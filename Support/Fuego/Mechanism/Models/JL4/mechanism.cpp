#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/* PURE CPU stuff  */
#ifndef AMREX_USE_GPU
namespace thermo
{
    amrex::Real fwd_A[4], fwd_beta[4], fwd_Ea[4];
    amrex::Real low_A[4], low_beta[4], low_Ea[4];
    amrex::Real rev_A[4], rev_beta[4], rev_Ea[4];
    amrex::Real troe_a[4],troe_Ts[4], troe_Tss[4], troe_Tsss[4];
    amrex::Real sri_a[4], sri_b[4], sri_c[4], sri_d[4], sri_e[4];
    amrex::Real activation_units[4], prefactor_units[4], phase_units[4];
    int is_PD[4], troe_len[4], sri_len[4], nTB[4], *TBid[4];
    amrex::Real *TB[4];
};

using namespace thermo;


/* Vectorized stuff  */

/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int *  np, amrex::Real *  y,  amrex::Real *  x)
{
    amrex::Real YOW[*np];
    amrex::Real imw[7];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<7; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<7; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}

/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tc[5], h[7];
    amrex::Real imw[7];

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
        hms[6*(*np)+i] = h[6];
    }

    for (int n=0; n<7; n++) {
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
    amrex::Real c[7*(*np)]; /*temporary storage */
    amrex::Real imw[7];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<7; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<7*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW[*np];
    amrex::Real imw[7];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<7; n++) {
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
    amrex::Real k_f_s[4*npt], Kc_s[4*npt], mixture[npt], g_RT[7*npt];
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

    for (int n=0; n<7; n++) {
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
        k_f_s[3*npt+i] = prefactor_units[3] * fwd_A[3] * exp(fwd_beta[3] * tc[i] - activation_units[3] * fwd_Ea[3] * invT[i]);
    }
}

void vcomp_gibbs(int npt, amrex::Real *  g_RT, amrex::Real *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        amrex::Real tg[5], g[7];
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
    }
}

void vcomp_Kc(int npt, amrex::Real *  Kc_s, amrex::Real *  g_RT, amrex::Real *  invT)
{
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        amrex::Real refC = (101325. / 8.31451) * invT[i];
        amrex::Real refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = pow(refC,3.000000) * exp((2.000000 * g_RT[0*npt+i] + g_RT[1*npt+i]) - (2.000000 * g_RT[4*npt+i] + 4.000000 * g_RT[6*npt+i]));
        Kc_s[1*npt+i] = pow(refC,2.000000) * exp((g_RT[0*npt+i] + g_RT[2*npt+i]) - (g_RT[4*npt+i] + 3.000000 * g_RT[6*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[1*npt+i] + 2.000000 * g_RT[6*npt+i]) - (2.000000 * g_RT[2*npt+i]));
        Kc_s[3*npt+i] = exp((g_RT[2*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
    }
}

void vcomp_wdot(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 1: 2.000000 CH4 + O2 => 2.000000 CO + 4.000000 H2 */
        phi_f = pow(sc[0*npt+i], 2.000000)*sc[1*npt+i];
        k_f = k_f_s[0*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= 2.000000 * qdot;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += 2.000000 * qdot;
        wdot[6*npt+i] += 4.000000 * qdot;

        /*reaction 2: CH4 + H2O <=> CO + 3.000000 H2 */
        phi_f = sc[0*npt+i]*sc[2*npt+i];
        k_f = k_f_s[1*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*pow(sc[6*npt+i], 3.000000);
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] += 3.000000 * qdot;

        /*reaction 3: 2.000000 H2 + O2 => 2.000000 H2O */
        phi_f = sc[1*npt+i]*pow(sc[6*npt+i], 2.000000);
        k_f = k_f_s[2*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += 2.000000 * qdot;
        wdot[6*npt+i] -= 2.000000 * qdot;

        /*reaction 4: CO + H2O <=> CO2 + H2 */
        phi_f = sc[2*npt+i]*sc[4*npt+i];
        k_f = k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
    }
}

static amrex::Real T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static amrex::Real k_f_save[4];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static amrex::Real Kc_save[4];
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

    amrex::Real qdot, q_f[4], q_r[4];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 7; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[0] -= 2.000000 * qdot;
    wdot[1] -= qdot;
    wdot[4] += 2.000000 * qdot;
    wdot[6] += 4.000000 * qdot;

    qdot = q_f[1]-q_r[1];
    wdot[0] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[6] += 3.000000 * qdot;

    qdot = q_f[2]-q_r[2];
    wdot[1] -= qdot;
    wdot[2] += 2.000000 * qdot;
    wdot[6] -= 2.000000 * qdot;

    qdot = q_f[3]-q_r[3];
    wdot[2] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;

    return;
}

void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f)
{
    for (int i=0; i<4; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(amrex::Real *  tc, amrex::Real invT, amrex::Real *  Kc)
{
    /*compute the Gibbs free energy */
    amrex::Real g_RT[7];
    gibbs(g_RT, tc);

    Kc[0] = 2.000000*g_RT[0] + g_RT[1] - 2.000000*g_RT[4] - 4.000000*g_RT[6];
    Kc[1] = g_RT[0] + g_RT[2] - g_RT[4] - 3.000000*g_RT[6];
    Kc[2] = g_RT[1] - 2.000000*g_RT[2] + 2.000000*g_RT[6];
    Kc[3] = g_RT[2] + g_RT[4] - g_RT[5] - g_RT[6];

    for (int i=0; i<4; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 * invT;
    amrex::Real refCinv = 1 / refC;

    Kc[0] *= pow(refC,3.000000);
    Kc[1] *= pow(refC,2.000000);
    Kc[2] *= refCinv;

    return;
}

void comp_qfqr(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT)
{

    /*reaction 1: 2.000000 CH4 + O2 => 2.000000 CO + 4.000000 H2 */
    qf[0] = pow(sc[0], 0.500000)*pow(sc[1], 1.250000);
    qr[0] = 0.0;

    /*reaction 2: CH4 + H2O <=> CO + 3.000000 H2 */
    qf[1] = sc[0]*sc[2];
    qr[1] = sc[4]*pow(sc[6], 3.000000);

    /*reaction 3: 2.000000 H2 + O2 => 2.000000 H2O */
    qf[2] = pow(sc[1], 1.500000)*pow(sc[6], 0.250000);
    qr[2] = 0.0;

    /*reaction 4: CO + H2O <=> CO2 + H2 */
    qf[3] = sc[2]*sc[4];
    qr[3] = sc[5]*sc[6];

    amrex::Real T = tc[1];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 7; ++i) {
        mixture += sc[i];
    }

    amrex::Real Corr[4];
    for (int i = 0; i < 4; ++i) {
        Corr[i] = 1.0;
    }

    for (int i=0; i<4; i++)
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

    amrex::Real q_f[4], q_r[4];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 4; ++i) {
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
    amrex::Real c[7]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 7; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 4; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 7; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 4; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[7]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[7];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*CH4 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*CO2 */
    YOW += y[6]*imw[6]; /*H2 */
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

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 4; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[7]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 7; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 4; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[7]; /*temporary storage */
    amrex::Real imw[7];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 4; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[7]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*16.043030; /*CH4 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*44.009950; /*CO2 */
    XW += x[6]*2.015940; /*H2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 7; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 4; ++id) {
        qdot[id] *= 1.0e-6;
    }
}

/*compute the reaction Jacobian on CPU */
void aJacobian(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP)
{
    for (int i=0; i<64; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[7];
    for (int k=0; k<7; k++) {
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
    for (int k = 0; k < 7; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[7];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[7];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[7];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: 2.000000 CH4 + O2 => 2.000000 CO + 4.000000 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[0], 2.000000)*sc[1];
    k_f = prefactor_units[0] * fwd_A[0]
                * exp(fwd_beta[0] * tc[0] - activation_units[0] * fwd_Ea[0] * invT);
    dlnkfdT = fwd_beta[0] * invT + activation_units[0] * fwd_Ea[0] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[0] -= 2 * q; /* CH4 */
    wdot[1] -= q; /* O2 */
    wdot[4] += 2 * q; /* CO */
    wdot[6] += 4 * q; /* H2 */
    /* d()/d[CH4] */
    dqdci =  + k_f*2.000000*sc[0]*sc[1];
    J[0] += -2 * dqdci;           /* dwdot[CH4]/d[CH4] */
    J[1] -= dqdci;                /* dwdot[O2]/d[CH4] */
    J[4] += 2 * dqdci;            /* dwdot[CO]/d[CH4] */
    J[6] += 4 * dqdci;            /* dwdot[H2]/d[CH4] */
    /* d()/d[O2] */
    dqdci =  + k_f*pow(sc[0], 2.000000);
    J[8] += -2 * dqdci;           /* dwdot[CH4]/d[O2] */
    J[9] -= dqdci;                /* dwdot[O2]/d[O2] */
    J[12] += 2 * dqdci;           /* dwdot[CO]/d[O2] */
    J[14] += 4 * dqdci;           /* dwdot[H2]/d[O2] */
    /* d()/dT */
    J[56] += -2 * dqdT;           /* dwdot[CH4]/dT */
    J[57] -= dqdT;                /* dwdot[O2]/dT */
    J[60] += 2 * dqdT;            /* dwdot[CO]/dT */
    J[62] += 4 * dqdT;            /* dwdot[H2]/dT */

    /*reaction 2: CH4 + H2O <=> CO + 3.000000 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[1] * fwd_A[1]
                * exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
    dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
    /* reverse */
    phi_r = sc[4]*pow(sc[6], 3.000000);
    Kc = pow(refC,2.000000) * exp(g_RT[0] + g_RT[2] - g_RT[4] - 3.000000*g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[4] + 3.000000*h_RT[6]) - 2.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* CH4 */
    wdot[2] -= q; /* H2O */
    wdot[4] += q; /* CO */
    wdot[6] += 3 * q; /* H2 */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[0] -= dqdci;                /* dwdot[CH4]/d[CH4] */
    J[2] -= dqdci;                /* dwdot[H2O]/d[CH4] */
    J[4] += dqdci;                /* dwdot[CO]/d[CH4] */
    J[6] += 3 * dqdci;            /* dwdot[H2]/d[CH4] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[0];
    J[16] -= dqdci;               /* dwdot[CH4]/d[H2O] */
    J[18] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[20] += dqdci;               /* dwdot[CO]/d[H2O] */
    J[22] += 3 * dqdci;           /* dwdot[H2]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*pow(sc[6], 3.000000);
    J[32] -= dqdci;               /* dwdot[CH4]/d[CO] */
    J[34] -= dqdci;               /* dwdot[H2O]/d[CO] */
    J[36] += dqdci;               /* dwdot[CO]/d[CO] */
    J[38] += 3 * dqdci;           /* dwdot[H2]/d[CO] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[4]*3.000000*pow(sc[6],2.000000);
    J[48] -= dqdci;               /* dwdot[CH4]/d[H2] */
    J[50] -= dqdci;               /* dwdot[H2O]/d[H2] */
    J[52] += dqdci;               /* dwdot[CO]/d[H2] */
    J[54] += 3 * dqdci;           /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[56] -= dqdT;                /* dwdot[CH4]/dT */
    J[58] -= dqdT;                /* dwdot[H2O]/dT */
    J[60] += dqdT;                /* dwdot[CO]/dT */
    J[62] += 3 * dqdT;            /* dwdot[H2]/dT */

    /*reaction 3: 2.000000 H2 + O2 => 2.000000 H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*pow(sc[6], 2.000000);
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[2] += 2 * q; /* H2O */
    wdot[6] -= 2 * q; /* H2 */
    /* d()/d[O2] */
    dqdci =  + k_f*pow(sc[6], 2.000000);
    J[9] -= dqdci;                /* dwdot[O2]/d[O2] */
    J[10] += 2 * dqdci;           /* dwdot[H2O]/d[O2] */
    J[14] += -2 * dqdci;          /* dwdot[H2]/d[O2] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[1]*2.000000*sc[6];
    J[49] -= dqdci;               /* dwdot[O2]/d[H2] */
    J[50] += 2 * dqdci;           /* dwdot[H2O]/d[H2] */
    J[54] += -2 * dqdci;          /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[57] -= dqdT;                /* dwdot[O2]/dT */
    J[58] += 2 * dqdT;            /* dwdot[H2O]/dT */
    J[62] += -2 * dqdT;           /* dwdot[H2]/dT */

    /*reaction 4: CO + H2O <=> CO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H2O */
    wdot[4] -= q; /* CO */
    wdot[5] += q; /* CO2 */
    wdot[6] += q; /* H2 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[4];
    J[18] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[20] -= dqdci;               /* dwdot[CO]/d[H2O] */
    J[21] += dqdci;               /* dwdot[CO2]/d[H2O] */
    J[22] += dqdci;               /* dwdot[H2]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[2];
    J[34] -= dqdci;               /* dwdot[H2O]/d[CO] */
    J[36] -= dqdci;               /* dwdot[CO]/d[CO] */
    J[37] += dqdci;               /* dwdot[CO2]/d[CO] */
    J[38] += dqdci;               /* dwdot[H2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[6];
    J[42] -= dqdci;               /* dwdot[H2O]/d[CO2] */
    J[44] -= dqdci;               /* dwdot[CO]/d[CO2] */
    J[45] += dqdci;               /* dwdot[CO2]/d[CO2] */
    J[46] += dqdci;               /* dwdot[H2]/d[CO2] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[50] -= dqdci;               /* dwdot[H2O]/d[H2] */
    J[52] -= dqdci;               /* dwdot[CO]/d[H2] */
    J[53] += dqdci;               /* dwdot[CO2]/d[H2] */
    J[54] += dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[58] -= dqdT;                /* dwdot[H2O]/dT */
    J[60] -= dqdT;                /* dwdot[CO]/dT */
    J[61] += dqdT;                /* dwdot[CO2]/dT */
    J[62] += dqdT;                /* dwdot[H2]/dT */

    amrex::Real c_R[7], dcRdT[7], e_RT[7];
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
    for (int k = 0; k < 7; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[56+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 7; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 7; ++m) {
            dehmixdc += eh_RT[m]*J[k*8+m];
        }
        J[k*8+7] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[63] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

#endif
/* Transport function declarations  */


void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 29;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 1148;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 7;}


void egtransetNLITE(int* NLITE ) {
    *NLITE = 1;}


/*Patm in ergs/cm3 */
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 1.60430300E+01;
    WT[1] = 3.19988000E+01;
    WT[2] = 1.80153400E+01;
    WT[3] = 2.80134000E+01;
    WT[4] = 2.80105500E+01;
    WT[5] = 4.40099500E+01;
    WT[6] = 2.01594000E+00;
}


/*the lennard-jones potential well depth eps/kb in K */
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 1.41400000E+02;
    EPS[1] = 1.07400000E+02;
    EPS[2] = 5.72400000E+02;
    EPS[3] = 9.75300000E+01;
    EPS[4] = 9.81000000E+01;
    EPS[5] = 2.44000000E+02;
    EPS[6] = 3.80000000E+01;
}


/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 3.74600000E+00;
    SIG[1] = 3.45800000E+00;
    SIG[2] = 2.60500000E+00;
    SIG[3] = 3.62100000E+00;
    SIG[4] = 3.65000000E+00;
    SIG[5] = 3.76300000E+00;
    SIG[6] = 2.92000000E+00;
}


/*the dipole moment in Debye */
void egtransetDIP(amrex::Real* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[2] = 1.84400000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 0.00000000E+00;
    DIP[6] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 2.60000000E+00;
    POL[1] = 1.60000000E+00;
    POL[2] = 0.00000000E+00;
    POL[3] = 1.76000000E+00;
    POL[4] = 1.95000000E+00;
    POL[5] = 2.65000000E+00;
    POL[6] = 7.90000000E-01;
}


/*the rotational relaxation collision number at 298 K */
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 1.30000000E+01;
    ZROT[1] = 3.80000000E+00;
    ZROT[2] = 4.00000000E+00;
    ZROT[3] = 4.00000000E+00;
    ZROT[4] = 1.80000000E+00;
    ZROT[5] = 2.10000000E+00;
    ZROT[6] = 2.80000000E+02;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 2;
    NLIN[1] = 1;
    NLIN[2] = 2;
    NLIN[3] = 1;
    NLIN[4] = 1;
    NLIN[5] = 1;
    NLIN[6] = 1;
}


/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -2.04642214E+01;
    COFETA[1] = 3.75363357E+00;
    COFETA[2] = -4.11756834E-01;
    COFETA[3] = 1.81766117E-02;
    COFETA[4] = -1.78915826E+01;
    COFETA[5] = 2.98311502E+00;
    COFETA[6] = -3.14105508E-01;
    COFETA[7] = 1.40500162E-02;
    COFETA[8] = -1.14441613E+01;
    COFETA[9] = -9.67162014E-01;
    COFETA[10] = 3.58651197E-01;
    COFETA[11] = -2.09789135E-02;
    COFETA[12] = -1.73976942E+01;
    COFETA[13] = 2.73482764E+00;
    COFETA[14] = -2.81916034E-01;
    COFETA[15] = 1.26588405E-02;
    COFETA[16] = -1.74469975E+01;
    COFETA[17] = 2.74728386E+00;
    COFETA[18] = -2.83509015E-01;
    COFETA[19] = 1.27267083E-02;
    COFETA[20] = -2.28110345E+01;
    COFETA[21] = 4.62954710E+00;
    COFETA[22] = -5.00689001E-01;
    COFETA[23] = 2.10012969E-02;
    COFETA[24] = -1.40419527E+01;
    COFETA[25] = 1.08789225E+00;
    COFETA[26] = -6.18592115E-02;
    COFETA[27] = 2.86838304E-03;
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = 1.76957851E+01;
    COFLAM[1] = -6.71274475E+00;
    COFLAM[2] = 1.26299606E+00;
    COFLAM[3] = -6.62748072E-02;
    COFLAM[4] = 5.31578403E-01;
    COFLAM[5] = 1.87067453E+00;
    COFLAM[6] = -1.31586198E-01;
    COFLAM[7] = 5.22416151E-03;
    COFLAM[8] = 1.81612053E+01;
    COFLAM[9] = -6.74137053E+00;
    COFLAM[10] = 1.21372119E+00;
    COFLAM[11] = -6.11027962E-02;
    COFLAM[12] = 7.77703429E+00;
    COFLAM[13] = -1.30957605E+00;
    COFLAM[14] = 3.28842203E-01;
    COFLAM[15] = -1.69485484E-02;
    COFLAM[16] = 8.17515440E+00;
    COFLAM[17] = -1.53836440E+00;
    COFLAM[18] = 3.68036945E-01;
    COFLAM[19] = -1.90917513E-02;
    COFLAM[20] = -8.74831432E+00;
    COFLAM[21] = 4.79275291E+00;
    COFLAM[22] = -4.18685061E-01;
    COFLAM[23] = 1.35210242E-02;
    COFLAM[24] = 4.34729192E+00;
    COFLAM[25] = 1.55347646E+00;
    COFLAM[26] = -1.60615552E-01;
    COFLAM[27] = 9.89934485E-03;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.75618457E+01;
    COFD[1] = 4.30617914E+00;
    COFD[2] = -3.46490389E-01;
    COFD[3] = 1.51071405E-02;
    COFD[4] = -1.68102713E+01;
    COFD[5] = 4.01337907E+00;
    COFD[6] = -3.10488902E-01;
    COFD[7] = 1.36288975E-02;
    COFD[8] = -1.95657633E+01;
    COFD[9] = 4.78813636E+00;
    COFD[10] = -3.65976055E-01;
    COFD[11] = 1.42215137E-02;
    COFD[12] = -1.65706884E+01;
    COFD[13] = 3.92005093E+00;
    COFD[14] = -2.99040611E-01;
    COFD[15] = 1.31607610E-02;
    COFD[16] = -1.65940140E+01;
    COFD[17] = 3.92553905E+00;
    COFD[18] = -2.99706984E-01;
    COFD[19] = 1.31876655E-02;
    COFD[20] = -1.93937264E+01;
    COFD[21] = 4.87146645E+00;
    COFD[22] = -4.13323360E-01;
    COFD[23] = 1.77408400E-02;
    COFD[24] = -1.29544834E+01;
    COFD[25] = 2.96758239E+00;
    COFD[26] = -1.76586224E-01;
    COFD[27] = 7.90559536E-03;
    COFD[28] = -1.68102713E+01;
    COFD[29] = 4.01337907E+00;
    COFD[30] = -3.10488902E-01;
    COFD[31] = 1.36288975E-02;
    COFD[32] = -1.60936570E+01;
    COFD[33] = 3.70633871E+00;
    COFD[34] = -2.71897253E-01;
    COFD[35] = 1.20097588E-02;
    COFD[36] = -1.99035583E+01;
    COFD[37] = 5.01694644E+00;
    COFD[38] = -4.08963011E-01;
    COFD[39] = 1.66143416E-02;
    COFD[40] = -1.58214936E+01;
    COFD[41] = 3.60000113E+00;
    COFD[42] = -2.58255120E-01;
    COFD[43] = 1.14251480E-02;
    COFD[44] = -1.58458281E+01;
    COFD[45] = 3.60600362E+00;
    COFD[46] = -2.59019961E-01;
    COFD[47] = 1.14576923E-02;
    COFD[48] = -1.87634092E+01;
    COFD[49] = 4.61060397E+00;
    COFD[50] = -3.83564503E-01;
    COFD[51] = 1.66168246E-02;
    COFD[52] = -1.22181183E+01;
    COFD[53] = 2.70415313E+00;
    COFD[54] = -1.41236971E-01;
    COFD[55] = 6.32236816E-03;
    COFD[56] = -1.95657633E+01;
    COFD[57] = 4.78813636E+00;
    COFD[58] = -3.65976055E-01;
    COFD[59] = 1.42215137E-02;
    COFD[60] = -1.99035583E+01;
    COFD[61] = 5.01694644E+00;
    COFD[62] = -4.08963011E-01;
    COFD[63] = 1.66143416E-02;
    COFD[64] = -1.16123849E+01;
    COFD[65] = 8.27754782E-01;
    COFD[66] = 2.52262233E-01;
    COFD[67] = -1.62567414E-02;
    COFD[68] = -1.99472346E+01;
    COFD[69] = 5.05636623E+00;
    COFD[70] = -4.17733686E-01;
    COFD[71] = 1.71403506E-02;
    COFD[72] = -1.99647405E+01;
    COFD[73] = 5.05179386E+00;
    COFD[74] = -4.16351103E-01;
    COFD[75] = 1.70488551E-02;
    COFD[76] = -1.82187624E+01;
    COFD[77] = 3.93854160E+00;
    COFD[78] = -2.28424632E-01;
    COFD[79] = 7.18603342E-03;
    COFD[80] = -1.73864044E+01;
    COFD[81] = 4.71143088E+00;
    COFD[82] = -3.95288626E-01;
    COFD[83] = 1.70702272E-02;
    COFD[84] = -1.65706884E+01;
    COFD[85] = 3.92005093E+00;
    COFD[86] = -2.99040611E-01;
    COFD[87] = 1.31607610E-02;
    COFD[88] = -1.58214936E+01;
    COFD[89] = 3.60000113E+00;
    COFD[90] = -2.58255120E-01;
    COFD[91] = 1.14251480E-02;
    COFD[92] = -1.99472346E+01;
    COFD[93] = 5.05636623E+00;
    COFD[94] = -4.17733686E-01;
    COFD[95] = 1.71403506E-02;
    COFD[96] = -1.56019563E+01;
    COFD[97] = 3.51542686E+00;
    COFD[98] = -2.47677471E-01;
    COFD[99] = 1.09841319E-02;
    COFD[100] = -1.56221627E+01;
    COFD[101] = 3.51977302E+00;
    COFD[102] = -2.48210923E-01;
    COFD[103] = 1.10059241E-02;
    COFD[104] = -1.85068873E+01;
    COFD[105] = 4.52122572E+00;
    COFD[106] = -3.73088946E-01;
    COFD[107] = 1.62076520E-02;
    COFD[108] = -1.20381391E+01;
    COFD[109] = 2.61421687E+00;
    COFD[110] = -1.28887086E-01;
    COFD[111] = 5.75609167E-03;
    COFD[112] = -1.65940140E+01;
    COFD[113] = 3.92553905E+00;
    COFD[114] = -2.99706984E-01;
    COFD[115] = 1.31876655E-02;
    COFD[116] = -1.58458281E+01;
    COFD[117] = 3.60600362E+00;
    COFD[118] = -2.59019961E-01;
    COFD[119] = 1.14576923E-02;
    COFD[120] = -1.99647405E+01;
    COFD[121] = 5.05179386E+00;
    COFD[122] = -4.16351103E-01;
    COFD[123] = 1.70488551E-02;
    COFD[124] = -1.56221627E+01;
    COFD[125] = 3.51977302E+00;
    COFD[126] = -2.48210923E-01;
    COFD[127] = 1.10059241E-02;
    COFD[128] = -1.56423580E+01;
    COFD[129] = 3.52412711E+00;
    COFD[130] = -2.48745351E-01;
    COFD[131] = 1.10277551E-02;
    COFD[132] = -1.85324360E+01;
    COFD[133] = 4.52748688E+00;
    COFD[134] = -3.73847542E-01;
    COFD[135] = 1.62384117E-02;
    COFD[136] = -1.20607690E+01;
    COFD[137] = 2.61969379E+00;
    COFD[138] = -1.29638429E-01;
    COFD[139] = 5.79050588E-03;
    COFD[140] = -1.93937264E+01;
    COFD[141] = 4.87146645E+00;
    COFD[142] = -4.13323360E-01;
    COFD[143] = 1.77408400E-02;
    COFD[144] = -1.87634092E+01;
    COFD[145] = 4.61060397E+00;
    COFD[146] = -3.83564503E-01;
    COFD[147] = 1.66168246E-02;
    COFD[148] = -1.82187624E+01;
    COFD[149] = 3.93854160E+00;
    COFD[150] = -2.28424632E-01;
    COFD[151] = 7.18603342E-03;
    COFD[152] = -1.85068873E+01;
    COFD[153] = 4.52122572E+00;
    COFD[154] = -3.73088946E-01;
    COFD[155] = 1.62076520E-02;
    COFD[156] = -1.85324360E+01;
    COFD[157] = 4.52748688E+00;
    COFD[158] = -3.73847542E-01;
    COFD[159] = 1.62384117E-02;
    COFD[160] = -2.05810669E+01;
    COFD[161] = 5.07469434E+00;
    COFD[162] = -4.25340301E-01;
    COFD[163] = 1.76800795E-02;
    COFD[164] = -1.43978662E+01;
    COFD[165] = 3.49721576E+00;
    COFD[166] = -2.45465191E-01;
    COFD[167] = 1.08948372E-02;
    COFD[168] = -1.29544834E+01;
    COFD[169] = 2.96758239E+00;
    COFD[170] = -1.76586224E-01;
    COFD[171] = 7.90559536E-03;
    COFD[172] = -1.22181183E+01;
    COFD[173] = 2.70415313E+00;
    COFD[174] = -1.41236971E-01;
    COFD[175] = 6.32236816E-03;
    COFD[176] = -1.73864044E+01;
    COFD[177] = 4.71143088E+00;
    COFD[178] = -3.95288626E-01;
    COFD[179] = 1.70702272E-02;
    COFD[180] = -1.20381391E+01;
    COFD[181] = 2.61421687E+00;
    COFD[182] = -1.28887086E-01;
    COFD[183] = 5.75609167E-03;
    COFD[184] = -1.20607690E+01;
    COFD[185] = 2.61969379E+00;
    COFD[186] = -1.29638429E-01;
    COFD[187] = 5.79050588E-03;
    COFD[188] = -1.43978662E+01;
    COFD[189] = 3.49721576E+00;
    COFD[190] = -2.45465191E-01;
    COFD[191] = 1.08948372E-02;
    COFD[192] = -1.04285080E+01;
    COFD[193] = 2.23477534E+00;
    COFD[194] = -8.11809423E-02;
    COFD[195] = 3.77342041E-03;
}


/*List of specs with small weight, dim NLITE */
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 6;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(amrex::Real* COFTD) {
    COFTD[0] = 2.98973958E-01;
    COFTD[1] = 2.32230992E-04;
    COFTD[2] = -1.23675023E-07;
    COFTD[3] = 2.01029713E-11;
    COFTD[4] = 3.81864172E-01;
    COFTD[5] = 1.84117353E-04;
    COFTD[6] = -9.79617476E-08;
    COFTD[7] = 1.62542227E-11;
    COFTD[8] = 1.95123509E-03;
    COFTD[9] = 6.69470998E-04;
    COFTD[10] = -3.12148757E-07;
    COFTD[11] = 4.52949938E-11;
    COFTD[12] = 3.88181125E-01;
    COFTD[13] = 1.55380218E-04;
    COFTD[14] = -8.20880914E-08;
    COFTD[15] = 1.37104636E-11;
    COFTD[16] = 3.87405318E-01;
    COFTD[17] = 1.56883797E-04;
    COFTD[18] = -8.29309791E-08;
    COFTD[19] = 1.38460299E-11;
    COFTD[20] = 2.47129011E-01;
    COFTD[21] = 4.49395677E-04;
    COFTD[22] = -2.32030740E-07;
    COFTD[23] = 3.62578797E-11;
    COFTD[24] = 0.00000000E+00;
    COFTD[25] = 0.00000000E+00;
    COFTD[26] = 0.00000000E+00;
    COFTD[27] = 0.00000000E+00;
}

#ifndef AMREX_USE_GPU
/* Initializes parameter database */
void CKINIT()
{

    // (0):  2.000000 CH4 + O2 => 2.000000 CO + 4.000000 H2
    // (0):  2.000000 CH4 + O2 => 2.000000 CO + 4.000000 H2
    fwd_A[0]     = 39100000000000;
    fwd_beta[0]  = 0;
    fwd_Ea[0]    = 30000;
    prefactor_units[0]  = 3.1622776601683795e-05;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = pow(10,-10.500000);
    is_PD[0] = 0;
    nTB[0] = 0;

    // (1):  CH4 + H2O <=> CO + 3.000000 H2
    // (1):  CH4 + H2O <=> CO + 3.000000 H2
    fwd_A[1]     = 300000000000;
    fwd_beta[1]  = 0;
    fwd_Ea[1]    = 30000;
    prefactor_units[1]  = 1.0000000000000002e-06;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = pow(10,-12.000000);
    is_PD[1] = 0;
    nTB[1] = 0;

    // (2):  2.000000 H2 + O2 => 2.000000 H2O
    // (2):  2.000000 H2 + O2 => 2.000000 H2O
    fwd_A[2]     = 6.045e+18;
    fwd_beta[2]  = -1;
    fwd_Ea[2]    = 40000;
    prefactor_units[2]  = 3.1622776601683795e-05;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-10.500000);
    is_PD[2] = 0;
    nTB[2] = 0;

    // (3):  CO + H2O <=> CO2 + H2
    // (3):  CO + H2O <=> CO2 + H2
    fwd_A[3]     = 2750000000000;
    fwd_beta[3]  = 0;
    fwd_Ea[3]    = 20000;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-12.000000);
    is_PD[3] = 0;
    nTB[3] = 0;

}


/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<4; ++i) {
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
    awt[0] = 12.011150; /*C */
    awt[1] = 15.999400; /*O */
    awt[2] = 1.007970; /*H */
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
    for (id = 0; id < kd * 7; ++ id) {
         ncf[id] = 0; 
    }

    /*CH4 */
    ncf[ 0 * kd + 0 ] = 1; /*C */
    ncf[ 0 * kd + 2 ] = 4; /*H */

    /*O2 */
    ncf[ 1 * kd + 1 ] = 2; /*O */

    /*H2O */
    ncf[ 2 * kd + 2 ] = 2; /*H */
    ncf[ 2 * kd + 1 ] = 1; /*O */

    /*N2 */
    ncf[ 3 * kd + 3 ] = 2; /*N */

    /*CO */
    ncf[ 4 * kd + 0 ] = 1; /*C */
    ncf[ 4 * kd + 1 ] = 1; /*O */

    /*CO2 */
    ncf[ 5 * kd + 0 ] = 1; /*C */
    ncf[ 5 * kd + 1 ] = 2; /*O */

    /*H2 */
    ncf[ 6 * kd + 2 ] = 2; /*H */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "C";
    ename[1] = "O";
    ename[2] = "H";
    ename[3] = "N";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(7);
    kname[0] = "CH4";
    kname[1] = "O2";
    kname[2] = "H2O";
    kname[3] = "N2";
    kname[4] = "CO";
    kname[5] = "CO2";
    kname[6] = "H2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(64);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(7);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[64];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<7; l++) {
                c_d[l] = 1.0/ 7.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<8; k++) {
        for (int l=0; l<8; l++) {
            if(J_h[ 8 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(64);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(7);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[64];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<7; k++) {
                c_d[k] = 1.0/ 7.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<8; k++) {
        for (int l=0; l<8; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 8 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(64);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(7);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[64];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<7; k++) {
                c_d[k] = 1.0/ 7.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<8; k++) {
        for (int l=0; l<8; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 8 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(64);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(7);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[64];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<7; k++) {
                c_d[k] = 1.0/ 7.000000 ;
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
        offset_row = nc * 8;
        offset_col = nc * 8;
        for (int k=0; k<8; k++) {
            for (int l=0; l<8; l++) {
                if(J_h[8*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(64);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(7);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[64];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<7; k++) {
                c_d[k] = 1.0/ 7.000000 ;
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
            offset = nc * 8;
            for (int l=0; l<8; l++) {
                for (int k=0; k<8; k++) {
                    if(J_h[8*k + l] != 0.0) {
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
            offset = nc * 8;
            for (int l=0; l<8; l++) {
                for (int k=0; k<8; k++) {
                    if(J_h[8*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(64);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(7);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[64];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<7; k++) {
                c_d[k] = 1.0/ 7.000000 ;
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
            offset = nc * 8;
            for (int l=0; l<8; l++) {
                for (int k=0; k<8; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[8*k + l] != 0.0) {
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
            offset = nc * 8;
            for (int l=0; l<8; l++) {
                for (int k=0; k<8; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[8*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(64);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(7);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[64];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<7; k++) {
                c_d[k] = 1.0/ 7.000000 ;
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
    for (int k=0; k<8; k++) {
        for (int l=0; l<8; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 8*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[8*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 8*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(64);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(7);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[64];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<7; k++) {
                c_d[k] = 1.0/ 7.000000 ;
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
        for (int l=0; l<8; l++) {
            for (int k=0; k<8; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[8*k + l] != 0.0) {
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
        for (int l=0; l<8; l++) {
            for (int k=0; k<8; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[8*k + l] != 0.0) {
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

