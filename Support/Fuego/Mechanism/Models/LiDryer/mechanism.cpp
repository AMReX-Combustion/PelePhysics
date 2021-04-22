#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"

/* PURE CPU stuff  */
#ifndef AMREX_USE_GPU
namespace thermo
{
    amrex::Real fwd_A[21], fwd_beta[21], fwd_Ea[21];
    amrex::Real low_A[21], low_beta[21], low_Ea[21];
    amrex::Real rev_A[21], rev_beta[21], rev_Ea[21];
    amrex::Real troe_a[21],troe_Ts[21], troe_Tss[21], troe_Tsss[21];
    amrex::Real sri_a[21], sri_b[21], sri_c[21], sri_d[21], sri_e[21];
    amrex::Real activation_units[21], prefactor_units[21], phase_units[21];
    int is_PD[21], troe_len[21], sri_len[21], nTB[21], *TBid[21];
    amrex::Real *TB[21];
};

using namespace thermo;


static amrex::Real T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static amrex::Real k_f_save[21];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static amrex::Real Kc_save[21];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif


/*compute the production rate for each species pointwise on CPU */
void productionRate_cpu(amrex::Real *  wdot, amrex::Real *  sc, amrex::Real T)
{
    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    amrex::Real qdot, q_f[21], q_r[21];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 9; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[5] += qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[2]-q_r[2];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[3] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[4] -= qdot;

    qdot = q_f[4]-q_r[4];
    wdot[3] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[9]-q_r[9];
    wdot[2] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[0] += qdot;
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[6] -= qdot;

    qdot = q_f[11]-q_r[11];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[12]-q_r[12];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[13]-q_r[13];
    wdot[1] += qdot;
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[6] -= qdot;

    qdot = q_f[14]-q_r[14];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[0] += qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    return;
}

void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f)
{
    for (int i=0; i<21; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(amrex::Real *  tc, amrex::Real invT, amrex::Real *  Kc)
{
    /*compute the Gibbs free energy */
    amrex::Real g_RT[9];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[1] + g_RT[3] - g_RT[6];
    Kc[1] = -g_RT[5] - g_RT[5] + g_RT[7];
    Kc[2] = g_RT[0] - g_RT[3] - g_RT[3];
    Kc[3] = -g_RT[1] + g_RT[4] + g_RT[4];
    Kc[4] = g_RT[3] + g_RT[4] - g_RT[5];
    Kc[5] = -g_RT[2] + g_RT[3] + g_RT[5];
    Kc[6] = g_RT[1] + g_RT[3] - g_RT[4] - g_RT[5];
    Kc[7] = g_RT[0] - g_RT[3] + g_RT[4] - g_RT[5];
    Kc[8] = g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5];
    Kc[9] = g_RT[2] + g_RT[4] - g_RT[5] - g_RT[5];
    Kc[10] = -g_RT[0] - g_RT[1] + g_RT[3] + g_RT[6];
    Kc[11] = g_RT[3] - g_RT[5] - g_RT[5] + g_RT[6];
    Kc[12] = -g_RT[1] + g_RT[4] - g_RT[5] + g_RT[6];
    Kc[13] = -g_RT[1] - g_RT[2] + g_RT[5] + g_RT[6];
    Kc[14] = -g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7];
    Kc[15] = -g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7];
    Kc[16] = -g_RT[2] + g_RT[3] - g_RT[5] + g_RT[7];
    Kc[17] = -g_RT[0] + g_RT[3] - g_RT[6] + g_RT[7];
    Kc[18] = g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7];
    Kc[19] = -g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7];
    Kc[20] = -g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7];

    for (int i=0; i<21; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 * invT;
    amrex::Real refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refC;
    Kc[2] *= refC;
    Kc[3] *= refCinv;
    Kc[4] *= refCinv;
    Kc[5] *= refCinv;

    return;
}

void comp_qfqr_cpu(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT)
{

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[1]*sc[3];
    qr[0] = sc[6];

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    qf[1] = sc[7];
    qr[1] = sc[5]*sc[5];

    /*reaction 3: H2 + M <=> H + H + M */
    qf[2] = sc[0];
    qr[2] = sc[3]*sc[3];

    /*reaction 4: O + O + M <=> O2 + M */
    qf[3] = sc[4]*sc[4];
    qr[3] = sc[1];

    /*reaction 5: O + H + M <=> OH + M */
    qf[4] = sc[3]*sc[4];
    qr[4] = sc[5];

    /*reaction 6: H + OH + M <=> H2O + M */
    qf[5] = sc[3]*sc[5];
    qr[5] = sc[2];

    /*reaction 7: H + O2 <=> O + OH */
    qf[6] = sc[1]*sc[3];
    qr[6] = sc[4]*sc[5];

    /*reaction 8: O + H2 <=> H + OH */
    qf[7] = sc[0]*sc[4];
    qr[7] = sc[3]*sc[5];

    /*reaction 9: H2 + OH <=> H2O + H */
    qf[8] = sc[0]*sc[5];
    qr[8] = sc[2]*sc[3];

    /*reaction 10: O + H2O <=> OH + OH */
    qf[9] = sc[2]*sc[4];
    qr[9] = sc[5]*sc[5];

    /*reaction 11: HO2 + H <=> H2 + O2 */
    qf[10] = sc[3]*sc[6];
    qr[10] = sc[0]*sc[1];

    /*reaction 12: HO2 + H <=> OH + OH */
    qf[11] = sc[3]*sc[6];
    qr[11] = sc[5]*sc[5];

    /*reaction 13: HO2 + O <=> O2 + OH */
    qf[12] = sc[4]*sc[6];
    qr[12] = sc[1]*sc[5];

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    qf[13] = sc[5]*sc[6];
    qr[13] = sc[1]*sc[2];

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    qf[14] = sc[6]*sc[6];
    qr[14] = sc[1]*sc[7];

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    qf[15] = sc[6]*sc[6];
    qr[15] = sc[1]*sc[7];

    /*reaction 17: H2O2 + H <=> H2O + OH */
    qf[16] = sc[3]*sc[7];
    qr[16] = sc[2]*sc[5];

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    qf[17] = sc[3]*sc[7];
    qr[17] = sc[0]*sc[6];

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    qf[18] = sc[4]*sc[7];
    qr[18] = sc[5]*sc[6];

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    qf[19] = sc[5]*sc[7];
    qr[19] = sc[2]*sc[6];

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    qf[20] = sc[5]*sc[7];
    qr[20] = sc[2]*sc[6];

    amrex::Real T = tc[1];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 9; ++i) {
        mixture += sc[i];
    }

    amrex::Real Corr[21];
    for (int i = 0; i < 21; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        amrex::Real alpha[2];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[2] + (TB[0][2] - 1)*sc[1];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[2];
        for (int i=0; i<2; i++)
        {
            amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;
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

    /* simple three-body correction */
    {
        amrex::Real alpha;
        alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[2];
        Corr[2] = alpha;
        alpha = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[2];
        Corr[3] = alpha;
        alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[2];
        Corr[4] = alpha;
        alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[2];
        Corr[5] = alpha;
    }

    for (int i=0; i<21; i++)
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

    amrex::Real q_f[21], q_r[21];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 21; ++i) {
        qdot[i] = q_f[i] - q_r[i];
    }

    return;
}


/*Returns the rate of progress for each reaction */
void CKQC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[9]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[9];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
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

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[9]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[9]; /*temporary storage */
    amrex::Real imw[9];

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

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[9]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}

/*compute the reaction Jacobian on CPU */
void aJacobian_cpu(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP)
{
    for (int i=0; i<100; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[9];
    for (int k=0; k<9; k++) {
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
    for (int k = 0; k < 9; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[9];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[9];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[9];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[2] + (TB[0][2] - 1)*sc[1];
    /* forward */
    phi_f = sc[1]*sc[3];
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
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[O2]/d[H2] */
        J[3] -= dqdci;                /* dwdot[H]/d[H2] */
        J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
        /* d()/d[O2] */
        dqdci = (TB[0][2] - 1)*dcdc_fac + k_f*sc[3];
        J[11] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[13] -= dqdci;               /* dwdot[H]/d[O2] */
        J[16] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[21] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[26] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[1];
        J[31] -= dqdci;               /* dwdot[O2]/d[H] */
        J[33] -= dqdci;               /* dwdot[H]/d[H] */
        J[36] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[61] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    }
    else {
        dqdc[0] = TB[0][0]*dcdc_fac;
        dqdc[1] = TB[0][2]*dcdc_fac + k_f*sc[3];
        dqdc[2] = TB[0][1]*dcdc_fac;
        dqdc[3] = dcdc_fac + k_f*sc[1];
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+3] -= dqdc[k];
            J[10*k+6] += dqdc[k];
        }
    }
    J[91] -= dqdT; /* dwdot[O2]/dT */
    J[93] -= dqdT; /* dwdot[H]/dT */
    J[96] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[2];
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
    phi_r = pow(sc[5], 2.000000);
    Kc = refC * exp(-g_RT[5] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2.000000*h_RT[5]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[5] += 2 * q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[5] += 2 * dqdci;            /* dwdot[OH]/d[H2] */
        J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[1][1] - 1)*dcdc_fac;
        J[25] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
        J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[OH] */
        dqdci =  - k_r*2.000000*sc[5];
        J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
        J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O2] */
        dqdci =  + k_f;
        J[75] += 2 * dqdci;           /* dwdot[OH]/d[H2O2] */
        J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    }
    else {
        dqdc[0] = TB[1][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[1][1]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac - k_r*2.000000*sc[5];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac + k_f;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+5] += 2 * dqdc[k];
            J[10*k+7] -= dqdc[k];
        }
    }
    J[95] += 2 * dqdT; /* dwdot[OH]/dT */
    J[97] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[2];
    /* forward */
    phi_f = sc[0];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* reverse */
    phi_r = pow(sc[3], 2.000000);
    Kc = refC * exp(g_RT[0] - g_RT[3] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0]) + (2.000000*h_RT[3]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[3] += 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[2][0] - 1)*q_nocor + k_f;
        J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
        J[3] += 2 * dqdci;            /* dwdot[H]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*q_nocor;
        J[20] -= dqdci;               /* dwdot[H2]/d[H2O] */
        J[23] += 2 * dqdci;           /* dwdot[H]/d[H2O] */
        /* d()/d[H] */
        dqdci =  - k_r*2.000000*sc[3];
        J[30] -= dqdci;               /* dwdot[H2]/d[H] */
        J[33] += 2 * dqdci;           /* dwdot[H]/d[H] */
    }
    else {
        dqdc[0] = TB[2][0]*q_nocor + k_f;
        dqdc[1] = q_nocor;
        dqdc[2] = TB[2][1]*q_nocor;
        dqdc[3] = q_nocor - k_r*2.000000*sc[3];
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+0] -= dqdc[k];
            J[10*k+3] += 2 * dqdc[k];
        }
    }
    J[90] -= dqdT; /* dwdot[H2]/dT */
    J[93] += 2 * dqdT; /* dwdot[H]/dT */

    /*reaction 4: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[2];
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[1];
    Kc = refCinv * exp(-g_RT[1] + g_RT[4] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[1]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[4] -= 2 * q; /* O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*q_nocor;
        J[1] += dqdci;                /* dwdot[O2]/d[H2] */
        J[4] += -2 * dqdci;           /* dwdot[O]/d[H2] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[11] += dqdci;               /* dwdot[O2]/d[O2] */
        J[14] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*q_nocor;
        J[21] += dqdci;               /* dwdot[O2]/d[H2O] */
        J[24] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[4];
        J[41] += dqdci;               /* dwdot[O2]/d[O] */
        J[44] += -2 * dqdci;          /* dwdot[O]/d[O] */
    }
    else {
        dqdc[0] = TB[3][0]*q_nocor;
        dqdc[1] = q_nocor - k_r;
        dqdc[2] = TB[3][1]*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f*2.000000*sc[4];
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+1] += dqdc[k];
            J[10*k+4] += -2 * dqdc[k];
        }
    }
    J[91] += dqdT; /* dwdot[O2]/dT */
    J[94] += -2 * dqdT; /* dwdot[O]/dT */

    /*reaction 5: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[2];
    /* forward */
    phi_f = sc[3]*sc[4];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[4]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* H */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[4][0] - 1)*q_nocor;
        J[3] -= dqdci;                /* dwdot[H]/d[H2] */
        J[4] -= dqdci;                /* dwdot[O]/d[H2] */
        J[5] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*q_nocor;
        J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[24] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[25] += dqdci;               /* dwdot[OH]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[33] -= dqdci;               /* dwdot[H]/d[H] */
        J[34] -= dqdci;               /* dwdot[O]/d[H] */
        J[35] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[3];
        J[43] -= dqdci;               /* dwdot[H]/d[O] */
        J[44] -= dqdci;               /* dwdot[O]/d[O] */
        J[45] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[53] -= dqdci;               /* dwdot[H]/d[OH] */
        J[54] -= dqdci;               /* dwdot[O]/d[OH] */
        J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    }
    else {
        dqdc[0] = TB[4][0]*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = TB[4][1]*q_nocor;
        dqdc[3] = q_nocor + k_f*sc[4];
        dqdc[4] = q_nocor + k_f*sc[3];
        dqdc[5] = q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+3] -= dqdc[k];
            J[10*k+4] -= dqdc[k];
            J[10*k+5] += dqdc[k];
        }
    }
    J[93] -= dqdT; /* dwdot[H]/dT */
    J[94] -= dqdT; /* dwdot[O]/dT */
    J[95] += dqdT; /* dwdot[OH]/dT */

    /*reaction 6: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[2];
    /* forward */
    phi_f = sc[3]*sc[5];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[2];
    Kc = refCinv * exp(-g_RT[2] + g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[2]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[3] -= q; /* H */
    wdot[5] -= q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[5][0] - 1)*q_nocor;
        J[2] += dqdci;                /* dwdot[H2O]/d[H2] */
        J[3] -= dqdci;                /* dwdot[H]/d[H2] */
        J[5] -= dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*q_nocor - k_r;
        J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
        J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[32] += dqdci;               /* dwdot[H2O]/d[H] */
        J[33] -= dqdci;               /* dwdot[H]/d[H] */
        J[35] -= dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[3];
        J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
        J[53] -= dqdci;               /* dwdot[H]/d[OH] */
        J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    }
    else {
        dqdc[0] = TB[5][0]*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = TB[5][1]*q_nocor - k_r;
        dqdc[3] = q_nocor + k_f*sc[5];
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor + k_f*sc[3];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+2] += dqdc[k];
            J[10*k+3] -= dqdc[k];
            J[10*k+5] -= dqdc[k];
        }
    }
    J[92] += dqdT; /* dwdot[H2O]/dT */
    J[93] -= dqdT; /* dwdot[H]/dT */
    J[95] -= dqdT; /* dwdot[OH]/dT */

    /*reaction 7: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[1] + g_RT[3] - g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[4] += q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[3];
    J[11] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[13] -= dqdci;               /* dwdot[H]/d[O2] */
    J[14] += dqdci;               /* dwdot[O]/d[O2] */
    J[15] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[1];
    J[31] -= dqdci;               /* dwdot[O2]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[34] += dqdci;               /* dwdot[O]/d[H] */
    J[35] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[41] -= dqdci;               /* dwdot[O2]/d[O] */
    J[43] -= dqdci;               /* dwdot[H]/d[O] */
    J[44] += dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[51] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[53] -= dqdci;               /* dwdot[H]/d[OH] */
    J[54] += dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[O2]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[94] += dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 8: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[0] - g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[3] += q; /* H */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[3] += dqdci;                /* dwdot[H]/d[H2] */
    J[4] -= dqdci;                /* dwdot[O]/d[H2] */
    J[5] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[30] -= dqdci;               /* dwdot[H2]/d[H] */
    J[33] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[O]/d[H] */
    J[35] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[40] -= dqdci;               /* dwdot[H2]/d[O] */
    J[43] += dqdci;               /* dwdot[H]/d[O] */
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[50] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[53] += dqdci;               /* dwdot[H]/d[OH] */
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[93] += dqdT;                /* dwdot[H]/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 9: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[3];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[2] += q; /* H2O */
    wdot[3] += q; /* H */
    wdot[5] -= q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[2] += dqdci;                /* dwdot[H2O]/d[H2] */
    J[3] += dqdci;                /* dwdot[H]/d[H2] */
    J[5] -= dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[20] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[23] += dqdci;               /* dwdot[H]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[2];
    J[30] -= dqdci;               /* dwdot[H2]/d[H] */
    J[32] += dqdci;               /* dwdot[H2O]/d[H] */
    J[33] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[50] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[53] += dqdci;               /* dwdot[H]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[93] += dqdT;                /* dwdot[H]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */

    /*reaction 10: O + H2O <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = pow(sc[5], 2.000000);
    Kc = exp(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2.000000*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H2O */
    wdot[4] -= q; /* O */
    wdot[5] += 2 * q; /* OH */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[4];
    J[22] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[24] -= dqdci;               /* dwdot[O]/d[H2O] */
    J[25] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[2];
    J[42] -= dqdci;               /* dwdot[H2O]/d[O] */
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += 2 * dqdci;           /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[5];
    J[52] -= dqdci;               /* dwdot[H2O]/d[OH] */
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[92] -= dqdT;                /* dwdot[H2O]/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += 2 * dqdT;            /* dwdot[OH]/dT */

    /*reaction 11: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[1];
    Kc = exp(-g_RT[0] - g_RT[1] + g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[0] + h_RT[1]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] += q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[1];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[O2]/d[H2] */
    J[3] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[10] += dqdci;               /* dwdot[H2]/d[O2] */
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[13] -= dqdci;               /* dwdot[H]/d[O2] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[30] += dqdci;               /* dwdot[H2]/d[H] */
    J[31] += dqdci;               /* dwdot[O2]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 12: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = pow(sc[5], 2.000000);
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (2.000000*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* H */
    wdot[5] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[5];
    J[53] -= dqdci;               /* dwdot[H]/d[OH] */
    J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[65] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[95] += 2 * dqdT;            /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 13: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(-g_RT[1] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[5];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[14] -= dqdci;               /* dwdot[O]/d[O2] */
    J[15] += dqdci;               /* dwdot[OH]/d[O2] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[41] += dqdci;               /* dwdot[O2]/d[O] */
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[51] += dqdci;               /* dwdot[O2]/d[OH] */
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[64] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[65] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[2];
    Kc = exp(-g_RT[1] - g_RT[2] + g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (h_RT[1] + h_RT[2]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[2];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[12] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[15] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[21] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[26] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[51] += dqdci;               /* dwdot[O2]/d[OH] */
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[62] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[65] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[1] + h_RT[7]));
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
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[16] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[17] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[71] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[1] + h_RT[7]));
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
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[16] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[17] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[71] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + g_RT[3] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[3] -= q; /* H */
    wdot[5] += q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[25] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[32] += dqdci;               /* dwdot[H2O]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += dqdci;               /* dwdot[OH]/d[H] */
    J[37] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[53] -= dqdci;               /* dwdot[H]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[72] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[73] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[3] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[3] -= q; /* H */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[3] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[30] += dqdci;               /* dwdot[H2]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] += dqdci;               /* dwdot[HO2]/d[H] */
    J[37] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[70] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[73] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    J[46] += dqdci;               /* dwdot[HO2]/d[O] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    J[56] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[65] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               /* dwdot[O]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[6];
    Kc = exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[26] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[56] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[2];
    J[62] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[65] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[72] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[75] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[6];
    Kc = exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[26] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[56] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[2];
    J[62] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[65] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[72] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[75] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    amrex::Real c_R[9], dcRdT[9], e_RT[9];
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
    for (int k = 0; k < 9; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[90+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 9; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 9; ++m) {
            dehmixdc += eh_RT[m]*J[k*10+m];
        }
        J[k*10+9] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[99] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

#endif


/* Initializes parameter database */
void CKINIT()
{
#ifndef AMREX_USE_GPU
    // (0):  H + O2 <=> O + OH
    fwd_A[6]     = 3547000000000000;
    fwd_beta[6]  = -0.40600000000000003;
    fwd_Ea[6]    = 16599;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = pow(10,-12.000000);
    is_PD[6] = 0;
    nTB[6] = 0;
    // (1):  O + H2 <=> H + OH
    fwd_A[7]     = 50800;
    fwd_beta[7]  = 2.6699999999999999;
    fwd_Ea[7]    = 6290;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = pow(10,-12.000000);
    is_PD[7] = 0;
    nTB[7] = 0;
    // (2):  H2 + OH <=> H2O + H
    fwd_A[8]     = 216000000;
    fwd_beta[8]  = 1.51;
    fwd_Ea[8]    = 3430;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = pow(10,-12.000000);
    is_PD[8] = 0;
    nTB[8] = 0;
    // (3):  O + H2O <=> OH + OH
    fwd_A[9]     = 2970000;
    fwd_beta[9]  = 2.02;
    fwd_Ea[9]    = 13400;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = pow(10,-12.000000);
    is_PD[9] = 0;
    nTB[9] = 0;
    // (4):  H2 + M <=> H + H + M
    fwd_A[2]     = 4.577e+19;
    fwd_beta[2]  = -1.3999999999999999;
    fwd_Ea[2]    = 104380;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-6.000000);
    is_PD[2] = 0;
    nTB[2] = 2;
    TB[2] = (amrex::Real *) malloc(2 * sizeof(amrex::Real));
    TBid[2] = (int *) malloc(2 * sizeof(int));
    TBid[2][0] = 0; TB[2][0] = 2.5; // H2
    TBid[2][1] = 2; TB[2][1] = 12; // H2O
    // (5):  O + O + M <=> O2 + M
    fwd_A[3]     = 6165000000000000;
    fwd_beta[3]  = -0.5;
    fwd_Ea[3]    = 0;
    prefactor_units[3]  = 1.0000000000000002e-12;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-12.000000);
    is_PD[3] = 0;
    nTB[3] = 2;
    TB[3] = (amrex::Real *) malloc(2 * sizeof(amrex::Real));
    TBid[3] = (int *) malloc(2 * sizeof(int));
    TBid[3][0] = 0; TB[3][0] = 2.5; // H2
    TBid[3][1] = 2; TB[3][1] = 12; // H2O
    // (6):  O + H + M <=> OH + M
    fwd_A[4]     = 4.714e+18;
    fwd_beta[4]  = -1;
    fwd_Ea[4]    = 0;
    prefactor_units[4]  = 1.0000000000000002e-12;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 0;
    nTB[4] = 2;
    TB[4] = (amrex::Real *) malloc(2 * sizeof(amrex::Real));
    TBid[4] = (int *) malloc(2 * sizeof(int));
    TBid[4][0] = 0; TB[4][0] = 2.5; // H2
    TBid[4][1] = 2; TB[4][1] = 12; // H2O
    // (7):  H + OH + M <=> H2O + M
    fwd_A[5]     = 3.8000000000000004e+22;
    fwd_beta[5]  = -2;
    fwd_Ea[5]    = 0;
    prefactor_units[5]  = 1.0000000000000002e-12;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-12.000000);
    is_PD[5] = 0;
    nTB[5] = 2;
    TB[5] = (amrex::Real *) malloc(2 * sizeof(amrex::Real));
    TBid[5] = (int *) malloc(2 * sizeof(int));
    TBid[5][0] = 0; TB[5][0] = 2.5; // H2
    TBid[5][1] = 2; TB[5][1] = 12; // H2O
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
    phase_units[0]      = pow(10,-12.000000);
    is_PD[0] = 1;
    nTB[0] = 3;
    TB[0] = (amrex::Real *) malloc(3 * sizeof(amrex::Real));
    TBid[0] = (int *) malloc(3 * sizeof(int));
    TBid[0][0] = 0; TB[0][0] = 2; // H2
    TBid[0][1] = 2; TB[0][1] = 11; // H2O
    TBid[0][2] = 1; TB[0][2] = 0.78000000000000003; // O2
    // (9):  HO2 + H <=> H2 + O2
    fwd_A[10]     = 16600000000000;
    fwd_beta[10]  = 0;
    fwd_Ea[10]    = 823;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = pow(10,-12.000000);
    is_PD[10] = 0;
    nTB[10] = 0;
    // (10):  HO2 + H <=> OH + OH
    fwd_A[11]     = 70790000000000;
    fwd_beta[11]  = 0;
    fwd_Ea[11]    = 295;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = pow(10,-12.000000);
    is_PD[11] = 0;
    nTB[11] = 0;
    // (11):  HO2 + O <=> O2 + OH
    fwd_A[12]     = 32500000000000;
    fwd_beta[12]  = 0;
    fwd_Ea[12]    = 0;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = pow(10,-12.000000);
    is_PD[12] = 0;
    nTB[12] = 0;
    // (12):  HO2 + OH <=> H2O + O2
    fwd_A[13]     = 28900000000000;
    fwd_beta[13]  = 0;
    fwd_Ea[13]    = -497;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = pow(10,-12.000000);
    is_PD[13] = 0;
    nTB[13] = 0;
    // (13):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[14]     = 420000000000000;
    fwd_beta[14]  = 0;
    fwd_Ea[14]    = 11982;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = pow(10,-12.000000);
    is_PD[14] = 0;
    nTB[14] = 0;
    // (14):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[15]     = 130000000000;
    fwd_beta[15]  = 0;
    fwd_Ea[15]    = -1629.3;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = pow(10,-12.000000);
    is_PD[15] = 0;
    nTB[15] = 0;
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
    phase_units[1]      = pow(10,-6.000000);
    is_PD[1] = 1;
    nTB[1] = 2;
    TB[1] = (amrex::Real *) malloc(2 * sizeof(amrex::Real));
    TBid[1] = (int *) malloc(2 * sizeof(int));
    TBid[1][0] = 0; TB[1][0] = 2.5; // H2
    TBid[1][1] = 2; TB[1][1] = 12; // H2O
    // (16):  H2O2 + H <=> H2O + OH
    fwd_A[16]     = 24100000000000;
    fwd_beta[16]  = 0;
    fwd_Ea[16]    = 3970;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = pow(10,-12.000000);
    is_PD[16] = 0;
    nTB[16] = 0;
    // (17):  H2O2 + H <=> HO2 + H2
    fwd_A[17]     = 48200000000000;
    fwd_beta[17]  = 0;
    fwd_Ea[17]    = 7950;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = pow(10,-12.000000);
    is_PD[17] = 0;
    nTB[17] = 0;
    // (18):  H2O2 + O <=> OH + HO2
    fwd_A[18]     = 9550000;
    fwd_beta[18]  = 2;
    fwd_Ea[18]    = 3970;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = pow(10,-12.000000);
    is_PD[18] = 0;
    nTB[18] = 0;
    // (19):  H2O2 + OH <=> HO2 + H2O
    fwd_A[19]     = 1000000000000;
    fwd_beta[19]  = 0;
    fwd_Ea[19]    = 0;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = pow(10,-12.000000);
    is_PD[19] = 0;
    nTB[19] = 0;
    // (20):  H2O2 + OH <=> HO2 + H2O
    fwd_A[20]     = 580000000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 9557;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-12.000000);
    is_PD[20] = 0;
    nTB[20] = 0;
#endif
}


/* Finalizes parameter database */
void CKFINALIZE()
{
#ifndef AMREX_USE_GPU
    for (int i=0; i<21; ++i) {
        free(TB[i]); TB[i] = 0; 
        free(TBid[i]); TBid[i] = 0;
        nTB[i] = 0;
    }
#endif
}


/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 1.007970; /*H */
    awt[1] = 15.999400; /*O */
    awt[2] = 14.006700; /*N */

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
    int kd = 3; 
    /*Zero ncf */
    for (id = 0; id < kd * 9; ++ id) {
         ncf[id] = 0; 
    }

    /*H2 */
    ncf[ 0 * kd + 0 ] = 2; /*H */

    /*O2 */
    ncf[ 1 * kd + 1 ] = 2; /*O */

    /*H2O */
    ncf[ 2 * kd + 0 ] = 2; /*H */
    ncf[ 2 * kd + 1 ] = 1; /*O */

    /*H */
    ncf[ 3 * kd + 0 ] = 1; /*H */

    /*O */
    ncf[ 4 * kd + 1 ] = 1; /*O */

    /*OH */
    ncf[ 5 * kd + 1 ] = 1; /*O */
    ncf[ 5 * kd + 0 ] = 1; /*H */

    /*HO2 */
    ncf[ 6 * kd + 0 ] = 1; /*H */
    ncf[ 6 * kd + 1 ] = 2; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 0 ] = 2; /*H */
    ncf[ 7 * kd + 1 ] = 2; /*O */

    /*N2 */
    ncf[ 8 * kd + 2 ] = 2; /*N */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(3);
    ename[0] = "H";
    ename[1] = "O";
    ename[2] = "N";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(9);
    kname[0] = "H2";
    kname[1] = "O2";
    kname[2] = "H2O";
    kname[3] = "H";
    kname[4] = "O";
    kname[5] = "OH";
    kname[6] = "HO2";
    kname[7] = "H2O2";
    kname[8] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(100);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(9);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[100];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<9; l++) {
                c_d[l] = 1.0/ 9.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<10; k++) {
        for (int l=0; l<10; l++) {
            if(J_h[ 10 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(100);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(9);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[100];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<9; k++) {
                c_d[k] = 1.0/ 9.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<10; k++) {
        for (int l=0; l<10; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 10 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(100);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(9);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[100];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<9; k++) {
                c_d[k] = 1.0/ 9.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<10; k++) {
        for (int l=0; l<10; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 10 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(100);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(9);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[100];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<9; k++) {
                c_d[k] = 1.0/ 9.000000 ;
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
        offset_row = nc * 10;
        offset_col = nc * 10;
        for (int k=0; k<10; k++) {
            for (int l=0; l<10; l++) {
                if(J_h[10*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(100);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(9);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[100];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<9; k++) {
                c_d[k] = 1.0/ 9.000000 ;
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
            offset = nc * 10;
            for (int l=0; l<10; l++) {
                for (int k=0; k<10; k++) {
                    if(J_h[10*k + l] != 0.0) {
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
            offset = nc * 10;
            for (int l=0; l<10; l++) {
                for (int k=0; k<10; k++) {
                    if(J_h[10*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(100);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(9);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[100];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<9; k++) {
                c_d[k] = 1.0/ 9.000000 ;
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
            offset = nc * 10;
            for (int l=0; l<10; l++) {
                for (int k=0; k<10; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[10*k + l] != 0.0) {
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
            offset = nc * 10;
            for (int l=0; l<10; l++) {
                for (int k=0; k<10; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[10*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(100);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(9);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[100];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<9; k++) {
                c_d[k] = 1.0/ 9.000000 ;
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
    for (int k=0; k<10; k++) {
        for (int l=0; l<10; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 10*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[10*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 10*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(100);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(9);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[100];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<9; k++) {
                c_d[k] = 1.0/ 9.000000 ;
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
        for (int l=0; l<10; l++) {
            for (int k=0; k<10; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[10*k + l] != 0.0) {
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
        for (int l=0; l<10; l++) {
            for (int k=0; k<10; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[10*k + l] != 0.0) {
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
