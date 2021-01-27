#ifndef MECHANISM_H
#define MECHANISM_H

#include "chemistry_file.H"
#include "chemistry_K.H"

/* PURE CPU stuff */
#if !defined(AMREX_USE_CUDA) && !defined(AMREX_USE_HIP)
namespace thermo
{
    double fwd_A[21], fwd_beta[21], fwd_Ea[21];
    double low_A[21], low_beta[21], low_Ea[21];
    double rev_A[21], rev_beta[21], rev_Ea[21];
    double troe_a[21],troe_Ts[21], troe_Tss[21], troe_Tsss[21];
    double sri_a[21], sri_b[21], sri_c[21], sri_d[21], sri_e[21];
    double activation_units[21], prefactor_units[21], phase_units[21];
    int is_PD[21], troe_len[21], sri_len[21], nTB[21], *TBid[21];
    double *TB[21];
    std::vector<std::vector<double>> kiv(21); 
    std::vector<std::vector<double>> nuv(21); 
};

using namespace thermo;

/* Vectorized stuff */
/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int *  np, double *  y,  double *  x)
{
    double YOW[*np];
    double imw[9];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}

/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, double *  T,  double *  hms)
{
    double tc[5], h[9];
    double imw[9];

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
        hms[7*(*np)+i] = h[7];
        hms[8*(*np)+i] = h[8];
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31446261815324e+07 * T[i] * imw[n];
        }
    }
}

/*Compute P = rhoRT/W(y) */
/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int *  np, double *  rho, double *  T,
	    double *  y,
	    double *  wdot)
{
    double c[9*(*np)]; /*temporary storage */
    double imw[9];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<9*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}

void VCKPY(int *  np, double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    double imw[9];

    get_imw(imw);

    for (int n=0; n<9; n++) {
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
void vproductionRate(int npt, double *  wdot, double *  sc, double *  T)
{
    double k_f_s[21*npt], Kc_s[21*npt], mixture[npt], g_RT[9*npt];
    double tc[5*npt], invT[npt];

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

    for (int n=0; n<9; n++) {
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
    }
}

void vcomp_gibbs(int npt, double *  g_RT, double *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[9];
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
    }
}

void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT)
{
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = (101325. / 8.31451) * invT[i];
        double refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[1*npt+i] = refC * exp((g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[2*npt+i] = refC * exp((g_RT[0*npt+i]) - (g_RT[3*npt+i] + g_RT[3*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[4*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[3*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[5*npt+i] = refCinv * exp((g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i]));
        Kc_s[6*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[7*npt+i] = exp((g_RT[0*npt+i] + g_RT[4*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[8*npt+i] = exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[3*npt+i]));
        Kc_s[9*npt+i] = exp((g_RT[2*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[10*npt+i] = exp((g_RT[3*npt+i] + g_RT[6*npt+i]) - (g_RT[0*npt+i] + g_RT[1*npt+i]));
        Kc_s[11*npt+i] = exp((g_RT[3*npt+i] + g_RT[6*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[12*npt+i] = exp((g_RT[4*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i]));
        Kc_s[13*npt+i] = exp((g_RT[5*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[2*npt+i]));
        Kc_s[14*npt+i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[7*npt+i]));
        Kc_s[15*npt+i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[7*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[0*npt+i] + g_RT[6*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[4*npt+i] + g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[2*npt+i] + g_RT[6*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[2*npt+i] + g_RT[6*npt+i]));
    }
}

void vcomp_wdot(int npt, double *  wdot, double *  mixture, double *  sc,
		double *  k_f_s, double *  Kc_s,
		double *  tc, double *  invT, double *  T)
{
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        double alpha;
        double redP, F;
        double logPred;
        double logFcent, troe_c, troe_n, troe, F_troe;

        /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[0*npt+i] + (TB[0][1] - 1)*sc[2*npt+i] + (TB[0][2] - 1)*sc[1*npt+i];
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
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
        phi_f = sc[7*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[0*npt+i] + (TB[1][1] - 1)*sc[2*npt+i];
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
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 3: H2 + M <=> H + H + M */
        phi_f = sc[0*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[0*npt+i] + (TB[2][1] - 1)*sc[2*npt+i];
        k_f = alpha * k_f_s[2*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[3*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 4: O + O + M <=> O2 + M */
        phi_f = sc[4*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[0*npt+i] + (TB[3][1] - 1)*sc[2*npt+i];
        k_f = alpha * k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;

        /*reaction 5: O + H + M <=> OH + M */
        phi_f = sc[3*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[0*npt+i] + (TB[4][1] - 1)*sc[2*npt+i];
        k_f = alpha * k_f_s[4*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 6: H + OH + M <=> H2O + M */
        phi_f = sc[3*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[0*npt+i] + (TB[5][1] - 1)*sc[2*npt+i];
        k_f = alpha * k_f_s[5*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 7: H + O2 <=> O + OH */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[6*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[5*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 8: O + H2 <=> H + OH */
        phi_f = sc[0*npt+i]*sc[4*npt+i];
        k_f = k_f_s[7*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 9: H2 + OH <=> H2O + H */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        k_f = k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[3*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 10: O + H2O <=> OH + OH */
        phi_f = sc[2*npt+i]*sc[4*npt+i];
        k_f = k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 11: HO2 + H <=> H2 + O2 */
        phi_f = sc[3*npt+i]*sc[6*npt+i];
        k_f = k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[1*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 12: HO2 + H <=> OH + OH */
        phi_f = sc[3*npt+i]*sc[6*npt+i];
        k_f = k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 13: HO2 + O <=> O2 + OH */
        phi_f = sc[4*npt+i]*sc[6*npt+i];
        k_f = k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 14: HO2 + OH <=> H2O + O2 */
        phi_f = sc[5*npt+i]*sc[6*npt+i];
        k_f = k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[2*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[7*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[7*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 17: H2O2 + H <=> H2O + OH */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 18: H2O2 + H <=> HO2 + H2 */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[6*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 19: H2O2 + O <=> OH + HO2 */
        phi_f = sc[4*npt+i]*sc[7*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 20: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[6*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 21: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[6*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
    }
}
#endif


/* INIT and FINALIZE stuff */
#if !defined(AMREX_USE_CUDA) && !defined(AMREX_USE_HIP)
/* Initializes parameter database */
void CKINIT()
{

    // (0):  H + O2 <=> O + OH
    kiv[6] = {3,1,4,5};
    nuv[6] = {-1,-1,1,1};
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
    kiv[7] = {4,0,3,5};
    nuv[7] = {-1,-1,1,1};
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
    kiv[8] = {0,5,2,3};
    nuv[8] = {-1,-1,1,1};
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
    kiv[9] = {4,2,5,5};
    nuv[9] = {-1,-1,1,1};
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
    kiv[2] = {0,3,3};
    nuv[2] = {-1,1,1};
    // (4):  H2 + M <=> H + H + M
    fwd_A[2]     = 4.577e+19;
    fwd_beta[2]  = -1.3999999999999999;
    fwd_Ea[2]    = 104380;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-6.000000);
    is_PD[2] = 0;
    nTB[2] = 2;
    TB[2] = (double *) malloc(2 * sizeof(double));
    TBid[2] = (int *) malloc(2 * sizeof(int));
    TBid[2][0] = 0; TB[2][0] = 2.5; // H2
    TBid[2][1] = 2; TB[2][1] = 12; // H2O

    // (5):  O + O + M <=> O2 + M
    kiv[3] = {4,4,1};
    nuv[3] = {-1,-1,1};
    // (5):  O + O + M <=> O2 + M
    fwd_A[3]     = 6165000000000000;
    fwd_beta[3]  = -0.5;
    fwd_Ea[3]    = 0;
    prefactor_units[3]  = 1.0000000000000002e-12;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-12.000000);
    is_PD[3] = 0;
    nTB[3] = 2;
    TB[3] = (double *) malloc(2 * sizeof(double));
    TBid[3] = (int *) malloc(2 * sizeof(int));
    TBid[3][0] = 0; TB[3][0] = 2.5; // H2
    TBid[3][1] = 2; TB[3][1] = 12; // H2O

    // (6):  O + H + M <=> OH + M
    kiv[4] = {4,3,5};
    nuv[4] = {-1,-1,1};
    // (6):  O + H + M <=> OH + M
    fwd_A[4]     = 4.714e+18;
    fwd_beta[4]  = -1;
    fwd_Ea[4]    = 0;
    prefactor_units[4]  = 1.0000000000000002e-12;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 0;
    nTB[4] = 2;
    TB[4] = (double *) malloc(2 * sizeof(double));
    TBid[4] = (int *) malloc(2 * sizeof(int));
    TBid[4][0] = 0; TB[4][0] = 2.5; // H2
    TBid[4][1] = 2; TB[4][1] = 12; // H2O

    // (7):  H + OH + M <=> H2O + M
    kiv[5] = {3,5,2};
    nuv[5] = {-1,-1,1};
    // (7):  H + OH + M <=> H2O + M
    fwd_A[5]     = 3.8000000000000004e+22;
    fwd_beta[5]  = -2;
    fwd_Ea[5]    = 0;
    prefactor_units[5]  = 1.0000000000000002e-12;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-12.000000);
    is_PD[5] = 0;
    nTB[5] = 2;
    TB[5] = (double *) malloc(2 * sizeof(double));
    TBid[5] = (int *) malloc(2 * sizeof(int));
    TBid[5][0] = 0; TB[5][0] = 2.5; // H2
    TBid[5][1] = 2; TB[5][1] = 12; // H2O

    // (8):  H + O2 (+M) <=> HO2 (+M)
    kiv[0] = {3,1,6};
    nuv[0] = {-1,-1,1};
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
    TB[0] = (double *) malloc(3 * sizeof(double));
    TBid[0] = (int *) malloc(3 * sizeof(int));
    TBid[0][0] = 0; TB[0][0] = 2; // H2
    TBid[0][1] = 2; TB[0][1] = 11; // H2O
    TBid[0][2] = 1; TB[0][2] = 0.78000000000000003; // O2

    // (9):  HO2 + H <=> H2 + O2
    kiv[10] = {6,3,0,1};
    nuv[10] = {-1,-1,1,1};
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
    kiv[11] = {6,3,5,5};
    nuv[11] = {-1,-1,1,1};
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
    kiv[12] = {6,4,1,5};
    nuv[12] = {-1,-1,1,1};
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
    kiv[13] = {6,5,2,1};
    nuv[13] = {-1,-1,1,1};
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
    kiv[14] = {6,6,7,1};
    nuv[14] = {-1,-1,1,1};
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
    kiv[15] = {6,6,7,1};
    nuv[15] = {-1,-1,1,1};
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
    kiv[1] = {7,5,5};
    nuv[1] = {-1,1,1};
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
    TB[1] = (double *) malloc(2 * sizeof(double));
    TBid[1] = (int *) malloc(2 * sizeof(int));
    TBid[1][0] = 0; TB[1][0] = 2.5; // H2
    TBid[1][1] = 2; TB[1][1] = 12; // H2O

    // (16):  H2O2 + H <=> H2O + OH
    kiv[16] = {7,3,2,5};
    nuv[16] = {-1,-1,1,1};
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
    kiv[17] = {7,3,6,0};
    nuv[17] = {-1,-1,1,1};
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
    kiv[18] = {7,4,5,6};
    nuv[18] = {-1,-1,1,1};
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
    kiv[19] = {7,5,6,2};
    nuv[19] = {-1,-1,1,1};
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
    kiv[20] = {7,5,6,2};
    nuv[20] = {-1,-1,1,1};
    // (20):  H2O2 + OH <=> HO2 + H2O
    fwd_A[20]     = 580000000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 9557;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-12.000000);
    is_PD[20] = 0;
    nTB[20] = 0;
}

/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<21; ++i) {
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
void atomicWeight(double *  awt)
{
    awt[0] = 1.007970; /*H */
    awt[1] = 15.999400; /*O */
    awt[2] = 14.006700; /*N */

    return;
}

/* MISC */ 
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


#if !defined(AMREX_USE_CUDA) && !defined(AMREX_USE_HIP)
/*Returns a count of species in a reaction, and their indices */
/*and stoichiometric coefficients. (Eq 50) */
void CKINU(int * i, int * nspec, int * ki, int * nu)
{
    if (*i < 1) {
        /*Return max num species per reaction */
        *nspec = 4;
    } else {
        if (*i > 21) {
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

/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*3; i++) {
        kname[i] = ' ';
    }

    /* H  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 2*lenkname + 0 ] = 'N';
    kname[ 2*lenkname + 1 ] = ' ';

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


/*compute the critical parameters for each species */
//void GET_CRITPARAMS(double *  Tci, double *  ai, double *  bi, double *  acentric_i)
//{
//
//    double   EPS[9];
//    double   SIG[9];
//    double    wt[9];
//    double avogadro = 6.02214199e23;
//    double boltzmann = 1.3806503e-16; //we work in CGS
//    double Rcst = 83.144598; //in bar [CGS] !
//
//    egtransetEPS(EPS);
//    egtransetSIG(SIG);
//    get_mw(wt);
//
//    /*species 0: H2 */
//    /*Imported from NIST */
//    Tci[0] = 33.145000 ; 
//    ai[0] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[0],2.0) / (pow(2.015880,2.0) * 12.964000); 
//    bi[0] = 0.08664 * Rcst * Tci[0] / (2.015880 * 12.964000); 
//    acentric_i[0] = -0.219000 ;
//
//    /*species 1: O2 */
//    /*Imported from NIST */
//    Tci[1] = 154.581000 ; 
//    ai[1] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[1],2.0) / (pow(31.998800,2.0) * 50.430466); 
//    bi[1] = 0.08664 * Rcst * Tci[1] / (31.998800 * 50.430466); 
//    acentric_i[1] = 0.022200 ;
//
//    /*species 2: H2O */
//    /*Imported from NIST */
//    Tci[2] = 647.096000 ; 
//    ai[2] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[2],2.0) / (pow(18.015340,2.0) * 220.640000); 
//    bi[2] = 0.08664 * Rcst * Tci[2] / (18.015340 * 220.640000); 
//    acentric_i[2] = 0.344300 ;
//
//    /*species 3: H */
//    Tci[3] = 1.316 * EPS[3] ; 
//    ai[3] = (5.55 * pow(avogadro,2.0) * EPS[3]*boltzmann * pow(1e-8*SIG[3],3.0) ) / (pow(wt[3],2.0)); 
//    bi[3] = 0.855 * avogadro * pow(1e-8*SIG[3],3.0) / (wt[3]); 
//    acentric_i[3] = 0.0 ;
//
//    /*species 4: O */
//    Tci[4] = 1.316 * EPS[4] ; 
//    ai[4] = (5.55 * pow(avogadro,2.0) * EPS[4]*boltzmann * pow(1e-8*SIG[4],3.0) ) / (pow(wt[4],2.0)); 
//    bi[4] = 0.855 * avogadro * pow(1e-8*SIG[4],3.0) / (wt[4]); 
//    acentric_i[4] = 0.0 ;
//
//    /*species 5: OH */
//    Tci[5] = 1.316 * EPS[5] ; 
//    ai[5] = (5.55 * pow(avogadro,2.0) * EPS[5]*boltzmann * pow(1e-8*SIG[5],3.0) ) / (pow(wt[5],2.0)); 
//    bi[5] = 0.855 * avogadro * pow(1e-8*SIG[5],3.0) / (wt[5]); 
//    acentric_i[5] = 0.0 ;
//
//    /*species 6: HO2 */
//    Tci[6] = 1.316 * EPS[6] ; 
//    ai[6] = (5.55 * pow(avogadro,2.0) * EPS[6]*boltzmann * pow(1e-8*SIG[6],3.0) ) / (pow(wt[6],2.0)); 
//    bi[6] = 0.855 * avogadro * pow(1e-8*SIG[6],3.0) / (wt[6]); 
//    acentric_i[6] = 0.0 ;
//
//    /*species 7: H2O2 */
//    Tci[7] = 1.316 * EPS[7] ; 
//    ai[7] = (5.55 * pow(avogadro,2.0) * EPS[7]*boltzmann * pow(1e-8*SIG[7],3.0) ) / (pow(wt[7],2.0)); 
//    bi[7] = 0.855 * avogadro * pow(1e-8*SIG[7],3.0) / (wt[7]); 
//    acentric_i[7] = 0.0 ;
//
//    /*species 8: N2 */
//    /*Imported from NIST */
//    Tci[8] = 126.192000 ; 
//    ai[8] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[8],2.0) / (pow(28.013400,2.0) * 33.958000); 
//    bi[8] = 0.08664 * Rcst * Tci[8] / (28.013400 * 33.958000); 
//    acentric_i[8] = 0.037200 ;
//
//    return;
//}

/*get atomic weight for all elements */
void CKAWT( double *  awt)
{
    atomicWeight(awt);
}

/*Transport function declarations */
void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 38;}

void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 1854;}

void egtransetNO(int* NO ) {
    *NO = 4;}

void egtransetKK(int* KK ) {
    *KK = 9;}

void egtransetNLITE(int* NLITE ) {
    *NLITE = 2;}

/*Patm in ergs/cm3 */
void egtransetPATM(double* PATM) {
    *PATM =   0.1013250000000000E+07;}

/*the molecular weights in g/mol */
void egtransetWT(double* WT ) {
    WT[0] = 2.01594000E+00;
    WT[1] = 3.19988000E+01;
    WT[2] = 1.80153400E+01;
    WT[3] = 1.00797000E+00;
    WT[4] = 1.59994000E+01;
    WT[5] = 1.70073700E+01;
    WT[6] = 3.30067700E+01;
    WT[7] = 3.40147400E+01;
    WT[8] = 2.80134000E+01;
}

/*the lennard-jones potential well depth eps/kb in K */
void egtransetEPS(double* EPS ) {
    EPS[0] = 3.80000000E+01;
    EPS[1] = 1.07400000E+02;
    EPS[2] = 5.72400000E+02;
    EPS[3] = 1.45000000E+02;
    EPS[4] = 8.00000000E+01;
    EPS[5] = 8.00000000E+01;
    EPS[6] = 1.07400000E+02;
    EPS[7] = 1.07400000E+02;
    EPS[8] = 9.75300000E+01;
}

/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(double* SIG ) {
    SIG[0] = 2.92000000E+00;
    SIG[1] = 3.45800000E+00;
    SIG[2] = 2.60500000E+00;
    SIG[3] = 2.05000000E+00;
    SIG[4] = 2.75000000E+00;
    SIG[5] = 2.75000000E+00;
    SIG[6] = 3.45800000E+00;
    SIG[7] = 3.45800000E+00;
    SIG[8] = 3.62100000E+00;
}

/*the dipole moment in Debye */
void egtransetDIP(double* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[2] = 1.84400000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 0.00000000E+00;
    DIP[6] = 0.00000000E+00;
    DIP[7] = 0.00000000E+00;
    DIP[8] = 0.00000000E+00;
}

/*the polarizability in cubic Angstroms */
void egtransetPOL(double* POL ) {
    POL[0] = 7.90000000E-01;
    POL[1] = 1.60000000E+00;
    POL[2] = 0.00000000E+00;
    POL[3] = 0.00000000E+00;
    POL[4] = 0.00000000E+00;
    POL[5] = 0.00000000E+00;
    POL[6] = 0.00000000E+00;
    POL[7] = 0.00000000E+00;
    POL[8] = 1.76000000E+00;
}

/*the rotational relaxation collision number at 298 K */
void egtransetZROT(double* ZROT ) {
    ZROT[0] = 2.80000000E+02;
    ZROT[1] = 3.80000000E+00;
    ZROT[2] = 4.00000000E+00;
    ZROT[3] = 0.00000000E+00;
    ZROT[4] = 0.00000000E+00;
    ZROT[5] = 0.00000000E+00;
    ZROT[6] = 1.00000000E+00;
    ZROT[7] = 3.80000000E+00;
    ZROT[8] = 4.00000000E+00;
}

/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 1;
    NLIN[1] = 1;
    NLIN[2] = 2;
    NLIN[3] = 0;
    NLIN[4] = 0;
    NLIN[5] = 1;
    NLIN[6] = 2;
    NLIN[7] = 2;
    NLIN[8] = 1;
}

/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(double* COFETA) {
    COFETA[0] = -1.37549435E+01;
    COFETA[1] = 9.65530587E-01;
    COFETA[2] = -4.45720114E-02;
    COFETA[3] = 2.05871810E-03;
    COFETA[4] = -1.68118868E+01;
    COFETA[5] = 2.52362554E+00;
    COFETA[6] = -2.49309128E-01;
    COFETA[7] = 1.10211025E-02;
    COFETA[8] = -1.17770937E+01;
    COFETA[9] = -8.26742721E-01;
    COFETA[10] = 3.39009079E-01;
    COFETA[11] = -2.00674327E-02;
    COFETA[12] = -1.98744496E+01;
    COFETA[13] = 3.41660514E+00;
    COFETA[14] = -3.63206306E-01;
    COFETA[15] = 1.58671021E-02;
    COFETA[16] = -1.48001581E+01;
    COFETA[17] = 1.79491990E+00;
    COFETA[18] = -1.54008440E-01;
    COFETA[19] = 6.86719439E-03;
    COFETA[20] = -1.47696103E+01;
    COFETA[21] = 1.79491990E+00;
    COFETA[22] = -1.54008440E-01;
    COFETA[23] = 6.86719439E-03;
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
}

/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(double* COFLAM) {
    COFLAM[0] = 1.11035551E+01;
    COFLAM[1] = -1.31884089E+00;
    COFLAM[2] = 2.44042735E-01;
    COFLAM[3] = -8.99837633E-03;
    COFLAM[4] = -2.51299170E+00;
    COFLAM[5] = 3.15166792E+00;
    COFLAM[6] = -3.10009273E-01;
    COFLAM[7] = 1.34523096E-02;
    COFLAM[8] = 2.21730076E+01;
    COFLAM[9] = -8.46933644E+00;
    COFLAM[10] = 1.46153515E+00;
    COFLAM[11] = -7.29500936E-02;
    COFLAM[12] = -3.24539191E-01;
    COFLAM[13] = 3.41660514E+00;
    COFLAM[14] = -3.63206306E-01;
    COFLAM[15] = 1.58671021E-02;
    COFLAM[16] = 1.98513952E+00;
    COFLAM[17] = 1.79491990E+00;
    COFLAM[18] = -1.54008440E-01;
    COFLAM[19] = 6.86719439E-03;
    COFLAM[20] = 1.60618734E+01;
    COFLAM[21] = -4.10626869E+00;
    COFLAM[22] = 6.63571339E-01;
    COFLAM[23] = -2.97906324E-02;
    COFLAM[24] = 5.56023763E-01;
    COFLAM[25] = 1.59073590E+00;
    COFLAM[26] = -5.28053839E-02;
    COFLAM[27] = 4.07601571E-04;
    COFLAM[28] = 1.48801857E+00;
    COFLAM[29] = 1.06175905E+00;
    COFLAM[30] = 5.72200266E-02;
    COFLAM[31] = -6.38393531E-03;
    COFLAM[32] = 1.15507419E+01;
    COFLAM[33] = -2.91453917E+00;
    COFLAM[34] = 5.55045765E-01;
    COFLAM[35] = -2.75173485E-02;
}

/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(double* COFD) {
    COFD[0] = -1.02395222E+01;
    COFD[1] = 2.15403244E+00;
    COFD[2] = -6.97480266E-02;
    COFD[3] = 3.23666871E-03;
    COFD[4] = -1.15797750E+01;
    COFD[5] = 2.43235504E+00;
    COFD[6] = -1.02890179E-01;
    COFD[7] = 4.52903603E-03;
    COFD[8] = -1.68758926E+01;
    COFD[9] = 4.49460303E+00;
    COFD[10] = -3.64766132E-01;
    COFD[11] = 1.56457153E-02;
    COFD[12] = -1.11808682E+01;
    COFD[13] = 2.66936727E+00;
    COFD[14] = -1.34411514E-01;
    COFD[15] = 5.92957488E-03;
    COFD[16] = -1.06250182E+01;
    COFD[17] = 2.15849701E+00;
    COFD[18] = -6.53886401E-02;
    COFD[19] = 2.81453370E-03;
    COFD[20] = -1.06283453E+01;
    COFD[21] = 2.15849701E+00;
    COFD[22] = -6.53886401E-02;
    COFD[23] = 2.81453370E-03;
    COFD[24] = -1.15806808E+01;
    COFD[25] = 2.43235504E+00;
    COFD[26] = -1.02890179E-01;
    COFD[27] = 4.52903603E-03;
    COFD[28] = -1.15815344E+01;
    COFD[29] = 2.43235504E+00;
    COFD[30] = -1.02890179E-01;
    COFD[31] = 4.52903603E-03;
    COFD[32] = -1.13253458E+01;
    COFD[33] = 2.31195095E+00;
    COFD[34] = -8.63988037E-02;
    COFD[35] = 3.77573452E-03;
    COFD[36] = -1.15797750E+01;
    COFD[37] = 2.43235504E+00;
    COFD[38] = -1.02890179E-01;
    COFD[39] = 4.52903603E-03;
    COFD[40] = -1.53110708E+01;
    COFD[41] = 3.37317428E+00;
    COFD[42] = -2.24900439E-01;
    COFD[43] = 9.81228151E-03;
    COFD[44] = -2.10640014E+01;
    COFD[45] = 5.50980695E+00;
    COFD[46] = -4.78335488E-01;
    COFD[47] = 1.98515434E-02;
    COFD[48] = -1.43712864E+01;
    COFD[49] = 3.70920439E+00;
    COFD[50] = -2.67274113E-01;
    COFD[51] = 1.15967481E-02;
    COFD[52] = -1.40864894E+01;
    COFD[53] = 3.07458927E+00;
    COFD[54] = -1.86899591E-01;
    COFD[55] = 8.19829781E-03;
    COFD[56] = -1.41066459E+01;
    COFD[57] = 3.07458927E+00;
    COFD[58] = -1.86899591E-01;
    COFD[59] = 8.19829781E-03;
    COFD[60] = -1.53187643E+01;
    COFD[61] = 3.37317428E+00;
    COFD[62] = -2.24900439E-01;
    COFD[63] = 9.81228151E-03;
    COFD[64] = -1.53261114E+01;
    COFD[65] = 3.37317428E+00;
    COFD[66] = -2.24900439E-01;
    COFD[67] = 9.81228151E-03;
    COFD[68] = -1.50096240E+01;
    COFD[69] = 3.25515933E+00;
    COFD[70] = -2.09710110E-01;
    COFD[71] = 9.15941830E-03;
    COFD[72] = -1.68758926E+01;
    COFD[73] = 4.49460303E+00;
    COFD[74] = -3.64766132E-01;
    COFD[75] = 1.56457153E-02;
    COFD[76] = -2.10640014E+01;
    COFD[77] = 5.50980695E+00;
    COFD[78] = -4.78335488E-01;
    COFD[79] = 1.98515434E-02;
    COFD[80] = -1.31492641E+01;
    COFD[81] = 1.48004311E+00;
    COFD[82] = 1.60499553E-01;
    COFD[83] = -1.19765679E-02;
    COFD[84] = -1.93611051E+01;
    COFD[85] = 5.51579726E+00;
    COFD[86] = -4.76061961E-01;
    COFD[87] = 1.96329391E-02;
    COFD[88] = -1.91096797E+01;
    COFD[89] = 5.02608697E+00;
    COFD[90] = -4.26959993E-01;
    COFD[91] = 1.80709910E-02;
    COFD[92] = -1.91256261E+01;
    COFD[93] = 5.02608697E+00;
    COFD[94] = -4.26959993E-01;
    COFD[95] = 1.80709910E-02;
    COFD[96] = -2.04177482E+01;
    COFD[97] = 5.31457079E+00;
    COFD[98] = -4.58216496E-01;
    COFD[99] = 1.91825910E-02;
    COFD[100] = -2.04230073E+01;
    COFD[101] = 5.31457079E+00;
    COFD[102] = -4.58216496E-01;
    COFD[103] = 1.91825910E-02;
    COFD[104] = -2.08123325E+01;
    COFD[105] = 5.42470154E+00;
    COFD[106] = -4.69700416E-01;
    COFD[107] = 1.95706904E-02;
    COFD[108] = -1.11808682E+01;
    COFD[109] = 2.66936727E+00;
    COFD[110] = -1.34411514E-01;
    COFD[111] = 5.92957488E-03;
    COFD[112] = -1.43712864E+01;
    COFD[113] = 3.70920439E+00;
    COFD[114] = -2.67274113E-01;
    COFD[115] = 1.15967481E-02;
    COFD[116] = -1.93611051E+01;
    COFD[117] = 5.51579726E+00;
    COFD[118] = -4.76061961E-01;
    COFD[119] = 1.96329391E-02;
    COFD[120] = -1.43693056E+01;
    COFD[121] = 4.03992999E+00;
    COFD[122] = -3.08044800E-01;
    COFD[123] = 1.32757775E-02;
    COFD[124] = -1.31860117E+01;
    COFD[125] = 3.38003453E+00;
    COFD[126] = -2.25783856E-01;
    COFD[127] = 9.85028660E-03;
    COFD[128] = -1.31877711E+01;
    COFD[129] = 3.38003453E+00;
    COFD[130] = -2.25783856E-01;
    COFD[131] = 9.85028660E-03;
    COFD[132] = -1.43717529E+01;
    COFD[133] = 3.70920439E+00;
    COFD[134] = -2.67274113E-01;
    COFD[135] = 1.15967481E-02;
    COFD[136] = -1.43721922E+01;
    COFD[137] = 3.70920439E+00;
    COFD[138] = -2.67274113E-01;
    COFD[139] = 1.15967481E-02;
    COFD[140] = -1.40298830E+01;
    COFD[141] = 3.55837688E+00;
    COFD[142] = -2.47785790E-01;
    COFD[143] = 1.07555332E-02;
    COFD[144] = -1.06250182E+01;
    COFD[145] = 2.15849701E+00;
    COFD[146] = -6.53886401E-02;
    COFD[147] = 2.81453370E-03;
    COFD[148] = -1.40864894E+01;
    COFD[149] = 3.07458927E+00;
    COFD[150] = -1.86899591E-01;
    COFD[151] = 8.19829781E-03;
    COFD[152] = -1.91096797E+01;
    COFD[153] = 5.02608697E+00;
    COFD[154] = -4.26959993E-01;
    COFD[155] = 1.80709910E-02;
    COFD[156] = -1.31860117E+01;
    COFD[157] = 3.38003453E+00;
    COFD[158] = -2.25783856E-01;
    COFD[159] = 9.85028660E-03;
    COFD[160] = -1.29877365E+01;
    COFD[161] = 2.80841511E+00;
    COFD[162] = -1.52629888E-01;
    COFD[163] = 6.72604927E-03;
    COFD[164] = -1.30027772E+01;
    COFD[165] = 2.80841511E+00;
    COFD[166] = -1.52629888E-01;
    COFD[167] = 6.72604927E-03;
    COFD[168] = -1.40916052E+01;
    COFD[169] = 3.07458927E+00;
    COFD[170] = -1.86899591E-01;
    COFD[171] = 8.19829781E-03;
    COFD[172] = -1.40964661E+01;
    COFD[173] = 3.07458927E+00;
    COFD[174] = -1.86899591E-01;
    COFD[175] = 8.19829781E-03;
    COFD[176] = -1.38756407E+01;
    COFD[177] = 2.98558426E+00;
    COFD[178] = -1.75507216E-01;
    COFD[179] = 7.71173691E-03;
    COFD[180] = -1.06283453E+01;
    COFD[181] = 2.15849701E+00;
    COFD[182] = -6.53886401E-02;
    COFD[183] = 2.81453370E-03;
    COFD[184] = -1.41066459E+01;
    COFD[185] = 3.07458927E+00;
    COFD[186] = -1.86899591E-01;
    COFD[187] = 8.19829781E-03;
    COFD[188] = -1.91256261E+01;
    COFD[189] = 5.02608697E+00;
    COFD[190] = -4.26959993E-01;
    COFD[191] = 1.80709910E-02;
    COFD[192] = -1.31877711E+01;
    COFD[193] = 3.38003453E+00;
    COFD[194] = -2.25783856E-01;
    COFD[195] = 9.85028660E-03;
    COFD[196] = -1.30027772E+01;
    COFD[197] = 2.80841511E+00;
    COFD[198] = -1.52629888E-01;
    COFD[199] = 6.72604927E-03;
    COFD[200] = -1.30182843E+01;
    COFD[201] = 2.80841511E+00;
    COFD[202] = -1.52629888E-01;
    COFD[203] = 6.72604927E-03;
    COFD[204] = -1.41119732E+01;
    COFD[205] = 3.07458927E+00;
    COFD[206] = -1.86899591E-01;
    COFD[207] = 8.19829781E-03;
    COFD[208] = -1.41170372E+01;
    COFD[209] = 3.07458927E+00;
    COFD[210] = -1.86899591E-01;
    COFD[211] = 8.19829781E-03;
    COFD[212] = -1.38948667E+01;
    COFD[213] = 2.98558426E+00;
    COFD[214] = -1.75507216E-01;
    COFD[215] = 7.71173691E-03;
    COFD[216] = -1.15806808E+01;
    COFD[217] = 2.43235504E+00;
    COFD[218] = -1.02890179E-01;
    COFD[219] = 4.52903603E-03;
    COFD[220] = -1.53187643E+01;
    COFD[221] = 3.37317428E+00;
    COFD[222] = -2.24900439E-01;
    COFD[223] = 9.81228151E-03;
    COFD[224] = -2.04177482E+01;
    COFD[225] = 5.31457079E+00;
    COFD[226] = -4.58216496E-01;
    COFD[227] = 1.91825910E-02;
    COFD[228] = -1.43717529E+01;
    COFD[229] = 3.70920439E+00;
    COFD[230] = -2.67274113E-01;
    COFD[231] = 1.15967481E-02;
    COFD[232] = -1.40916052E+01;
    COFD[233] = 3.07458927E+00;
    COFD[234] = -1.86899591E-01;
    COFD[235] = 8.19829781E-03;
    COFD[236] = -1.41119732E+01;
    COFD[237] = 3.07458927E+00;
    COFD[238] = -1.86899591E-01;
    COFD[239] = 8.19829781E-03;
    COFD[240] = -1.53265780E+01;
    COFD[241] = 3.37317428E+00;
    COFD[242] = -2.24900439E-01;
    COFD[243] = 9.81228151E-03;
    COFD[244] = -1.53340417E+01;
    COFD[245] = 3.37317428E+00;
    COFD[246] = -2.24900439E-01;
    COFD[247] = 9.81228151E-03;
    COFD[248] = -1.50168028E+01;
    COFD[249] = 3.25515933E+00;
    COFD[250] = -2.09710110E-01;
    COFD[251] = 9.15941830E-03;
    COFD[252] = -1.15815344E+01;
    COFD[253] = 2.43235504E+00;
    COFD[254] = -1.02890179E-01;
    COFD[255] = 4.52903603E-03;
    COFD[256] = -1.53261114E+01;
    COFD[257] = 3.37317428E+00;
    COFD[258] = -2.24900439E-01;
    COFD[259] = 9.81228151E-03;
    COFD[260] = -2.04230073E+01;
    COFD[261] = 5.31457079E+00;
    COFD[262] = -4.58216496E-01;
    COFD[263] = 1.91825910E-02;
    COFD[264] = -1.43721922E+01;
    COFD[265] = 3.70920439E+00;
    COFD[266] = -2.67274113E-01;
    COFD[267] = 1.15967481E-02;
    COFD[268] = -1.40964661E+01;
    COFD[269] = 3.07458927E+00;
    COFD[270] = -1.86899591E-01;
    COFD[271] = 8.19829781E-03;
    COFD[272] = -1.41170372E+01;
    COFD[273] = 3.07458927E+00;
    COFD[274] = -1.86899591E-01;
    COFD[275] = 8.19829781E-03;
    COFD[276] = -1.53340417E+01;
    COFD[277] = 3.37317428E+00;
    COFD[278] = -2.24900439E-01;
    COFD[279] = 9.81228151E-03;
    COFD[280] = -1.53416186E+01;
    COFD[281] = 3.37317428E+00;
    COFD[282] = -2.24900439E-01;
    COFD[283] = 9.81228151E-03;
    COFD[284] = -1.50236516E+01;
    COFD[285] = 3.25515933E+00;
    COFD[286] = -2.09710110E-01;
    COFD[287] = 9.15941830E-03;
    COFD[288] = -1.13253458E+01;
    COFD[289] = 2.31195095E+00;
    COFD[290] = -8.63988037E-02;
    COFD[291] = 3.77573452E-03;
    COFD[292] = -1.50096240E+01;
    COFD[293] = 3.25515933E+00;
    COFD[294] = -2.09710110E-01;
    COFD[295] = 9.15941830E-03;
    COFD[296] = -2.08123325E+01;
    COFD[297] = 5.42470154E+00;
    COFD[298] = -4.69700416E-01;
    COFD[299] = 1.95706904E-02;
    COFD[300] = -1.40298830E+01;
    COFD[301] = 3.55837688E+00;
    COFD[302] = -2.47785790E-01;
    COFD[303] = 1.07555332E-02;
    COFD[304] = -1.38756407E+01;
    COFD[305] = 2.98558426E+00;
    COFD[306] = -1.75507216E-01;
    COFD[307] = 7.71173691E-03;
    COFD[308] = -1.38948667E+01;
    COFD[309] = 2.98558426E+00;
    COFD[310] = -1.75507216E-01;
    COFD[311] = 7.71173691E-03;
    COFD[312] = -1.50168028E+01;
    COFD[313] = 3.25515933E+00;
    COFD[314] = -2.09710110E-01;
    COFD[315] = 9.15941830E-03;
    COFD[316] = -1.50236516E+01;
    COFD[317] = 3.25515933E+00;
    COFD[318] = -2.09710110E-01;
    COFD[319] = 9.15941830E-03;
    COFD[320] = -1.47639290E+01;
    COFD[321] = 3.15955654E+00;
    COFD[322] = -1.97590757E-01;
    COFD[323] = 8.64692156E-03;
}

/*List of specs with small weight, dim NLITE */
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 1;
    KTDIF[1] = 4;
}

/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(double* COFTD) {
    COFTD[0] = 0.00000000E+00;
    COFTD[1] = 0.00000000E+00;
    COFTD[2] = 0.00000000E+00;
    COFTD[3] = 0.00000000E+00;
    COFTD[4] = 4.42739084E-01;
    COFTD[5] = 7.11770818E-05;
    COFTD[6] = -3.84768062E-08;
    COFTD[7] = 6.86323437E-12;
    COFTD[8] = 6.02028221E-02;
    COFTD[9] = 5.61561867E-04;
    COFTD[10] = -2.55372862E-07;
    COFTD[11] = 3.63389913E-11;
    COFTD[12] = -1.52534742E-01;
    COFTD[13] = -5.46404022E-05;
    COFTD[14] = 2.93412470E-08;
    COFTD[15] = -4.87091914E-12;
    COFTD[16] = 4.15583337E-01;
    COFTD[17] = 1.09738399E-05;
    COFTD[18] = -3.96021963E-09;
    COFTD[19] = 1.14414443E-12;
    COFTD[20] = 4.21932443E-01;
    COFTD[21] = 1.11414935E-05;
    COFTD[22] = -4.02072219E-09;
    COFTD[23] = 1.16162418E-12;
    COFTD[24] = 4.44452569E-01;
    COFTD[25] = 7.14525507E-05;
    COFTD[26] = -3.86257187E-08;
    COFTD[27] = 6.88979640E-12;
    COFTD[28] = 4.46070183E-01;
    COFTD[29] = 7.17126069E-05;
    COFTD[30] = -3.87662996E-08;
    COFTD[31] = 6.91487226E-12;
    COFTD[32] = 4.45261966E-01;
    COFTD[33] = 4.94697174E-05;
    COFTD[34] = -2.63023442E-08;
    COFTD[35] = 4.90306217E-12;
    COFTD[36] = 1.52534742E-01;
    COFTD[37] = 5.46404022E-05;
    COFTD[38] = -2.93412470E-08;
    COFTD[39] = 4.87091914E-12;
    COFTD[40] = 2.20482843E-01;
    COFTD[41] = 4.80164288E-04;
    COFTD[42] = -2.32927944E-07;
    COFTD[43] = 3.46470436E-11;
    COFTD[44] = -1.41883744E-01;
    COFTD[45] = 7.66558810E-04;
    COFTD[46] = -3.06550003E-07;
    COFTD[47] = 4.02959502E-11;
    COFTD[48] = 0.00000000E+00;
    COFTD[49] = 0.00000000E+00;
    COFTD[50] = 0.00000000E+00;
    COFTD[51] = 0.00000000E+00;
    COFTD[52] = 2.70010150E-01;
    COFTD[53] = 3.61555093E-04;
    COFTD[54] = -1.80744752E-07;
    COFTD[55] = 2.75321248E-11;
    COFTD[56] = 2.72041664E-01;
    COFTD[57] = 3.64275376E-04;
    COFTD[58] = -1.82104647E-07;
    COFTD[59] = 2.77392722E-11;
    COFTD[60] = 2.20907853E-01;
    COFTD[61] = 4.81089870E-04;
    COFTD[62] = -2.33376944E-07;
    COFTD[63] = 3.47138305E-11;
    COFTD[64] = 2.21308399E-01;
    COFTD[65] = 4.81962174E-04;
    COFTD[66] = -2.33800100E-07;
    COFTD[67] = 3.47767730E-11;
    COFTD[68] = 2.40744421E-01;
    COFTD[69] = 4.45343451E-04;
    COFTD[70] = -2.18173874E-07;
    COFTD[71] = 3.26958506E-11;
}

/* SPARSE INFORMATION */ 
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
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<9; k++) {
                c_d[k] = 1.0/ 9.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
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
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<9; k++) {
                c_d[k] = 1.0/ 9.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
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
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<9; k++) {
                c_d[k] = 1.0/ 9.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
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
    double c[9];
    double J[100];
    int offset_row;
    int offset_col;

    //for (int k=0; k<9; k++) {
    //    c[k] = 1.0/ 9.000000 ;
    //}

    //amrex::single_task( [=] 
    //AMREX_GPU_DEVICE () noexcept 
    //{
    //    aJacobian(J, c, 1500.0, *consP);
    //});

    //colPtrs[0] = 0;
    //int nJdata_tmp = 0;
    //for (int nc=0; nc<NCELLS; nc++) {
    //    offset_row = nc * 10;
    //    offset_col = nc * 10;
    //    for (int k=0; k<10; k++) {
    //        for (int l=0; l<10; l++) {
    //            if(J[10*k + l] != 0.0) {
    //                rowVals[nJdata_tmp] = l + offset_row; 
    //                nJdata_tmp = nJdata_tmp + 1; 
    //            }
    //        }
    //        colPtrs[offset_col + (k + 1)] = nJdata_tmp;
    //    }
    //}

    return;
}

/*compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0 */
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base)
{
    double c[9];
    double J[100];
    int offset;

    for (int k=0; k<9; k++) {
        c[k] = 1.0/ 9.000000 ;
    }

    //amrex::single_task( [=] 
    //AMREX_GPU_DEVICE () noexcept 
    //{
    //    aJacobian(J, c, 1500.0, *consP);
    //});

    //if (base == 1) {
    //    rowPtrs[0] = 1;
    //    int nJdata_tmp = 1;
    //    for (int nc=0; nc<NCELLS; nc++) {
    //        offset = nc * 10;
    //        for (int l=0; l<10; l++) {
    //            for (int k=0; k<10; k++) {
    //                if(J[10*k + l] != 0.0) {
    //                    colVals[nJdata_tmp-1] = k+1 + offset; 
    //                    nJdata_tmp = nJdata_tmp + 1; 
    //                }
    //            }
    //            rowPtrs[offset + (l + 1)] = nJdata_tmp;
    //        }
    //    }
    //} else {
    //    rowPtrs[0] = 0;
    //    int nJdata_tmp = 0;
    //    for (int nc=0; nc<NCELLS; nc++) {
    //        offset = nc * 10;
    //        for (int l=0; l<10; l++) {
    //            for (int k=0; k<10; k++) {
    //                if(J[10*k + l] != 0.0) {
    //                    colVals[nJdata_tmp] = k + offset; 
    //                    nJdata_tmp = nJdata_tmp + 1; 
    //                }
    //            }
    //            rowPtrs[offset + (l + 1)] = nJdata_tmp;
    //        }
    //    }
    //}

    return;
}

/*compute the sparsity pattern of the system Jacobian */
/*CSR format BASE is user choice */
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, int * consP, int NCELLS, int base)
{
    double c[9];
    double J[100];
    int offset;

    for (int k=0; k<9; k++) {
        c[k] = 1.0/ 9.000000 ;
    }

    //amrex::single_task( [=] 
    //AMREX_GPU_DEVICE () noexcept 
    //{
    //    aJacobian(J, c, 1500.0, *consP);
    //});

    //if (base == 1) {
    //    rowPtr[0] = 1;
    //    int nJdata_tmp = 1;
    //    for (int nc=0; nc<NCELLS; nc++) {
    //        offset = nc * 10;
    //        for (int l=0; l<10; l++) {
    //            for (int k=0; k<10; k++) {
    //                if (k == l) {
    //                    colVals[nJdata_tmp-1] = l+1 + offset; 
    //                    nJdata_tmp = nJdata_tmp + 1; 
    //                } else {
    //                    if(J[10*k + l] != 0.0) {
    //                        colVals[nJdata_tmp-1] = k+1 + offset; 
    //                        nJdata_tmp = nJdata_tmp + 1; 
    //                    }
    //                }
    //            }
    //            rowPtr[offset + (l + 1)] = nJdata_tmp;
    //        }
    //    }
    //} else {
    //    rowPtr[0] = 0;
    //    int nJdata_tmp = 0;
    //    for (int nc=0; nc<NCELLS; nc++) {
    //        offset = nc * 10;
    //        for (int l=0; l<10; l++) {
    //            for (int k=0; k<10; k++) {
    //                if (k == l) {
    //                    colVals[nJdata_tmp] = l + offset; 
    //                    nJdata_tmp = nJdata_tmp + 1; 
    //                } else {
    //                    if(J[10*k + l] != 0.0) {
    //                        colVals[nJdata_tmp] = k + offset; 
    //                        nJdata_tmp = nJdata_tmp + 1; 
    //                    }
    //                }
    //            }
    //            rowPtr[offset + (l + 1)] = nJdata_tmp;
    //        }
    //    }
    //}

    return;
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU */
/*BASE 0 */
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP)
{
    double c[9];
    double J[100];

    for (int k=0; k<9; k++) {
        c[k] = 1.0/ 9.000000 ;
    }

    //amrex::single_task( [=] 
    //AMREX_GPU_DEVICE () noexcept 
    //{
    //    aJacobian_precond(J, c, 1500.0, *consP);
    //});

    //colPtrs[0] = 0;
    //int nJdata_tmp = 0;
    //for (int k=0; k<10; k++) {
    //    for (int l=0; l<10; l++) {
    //        if (k == l) {
    //            rowVals[nJdata_tmp] = l; 
    //            indx[nJdata_tmp] = 10*k + l;
    //            nJdata_tmp = nJdata_tmp + 1; 
    //        } else {
    //            if(J[10*k + l] != 0.0) {
    //                rowVals[nJdata_tmp] = l; 
    //                indx[nJdata_tmp] = 10*k + l;
    //                nJdata_tmp = nJdata_tmp + 1; 
    //            }
    //        }
    //    }
    //    colPtrs[k+1] = nJdata_tmp;
    //}

    return;
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian */
/*CSR format BASE is under choice */
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base)
{
    double c[9];
    double J[100];

    for (int k=0; k<9; k++) {
        c[k] = 1.0/ 9.000000 ;
    }

    //amrex::single_task( [=] 
    //AMREX_GPU_DEVICE () noexcept 
    //{
    //    aJacobian_precond(J, c, 1500.0, *consP);
    //});

    //if (base == 1) {
    //    rowPtr[0] = 1;
    //    int nJdata_tmp = 1;
    //    for (int l=0; l<10; l++) {
    //        for (int k=0; k<10; k++) {
    //            if (k == l) {
    //                colVals[nJdata_tmp-1] = l+1; 
    //                nJdata_tmp = nJdata_tmp + 1; 
    //            } else {
    //                if(J[10*k + l] != 0.0) {
    //                    colVals[nJdata_tmp-1] = k+1; 
    //                    nJdata_tmp = nJdata_tmp + 1; 
    //                }
    //            }
    //        }
    //        rowPtr[l+1] = nJdata_tmp;
    //    }
    //} else {
    //    rowPtr[0] = 0;
    //    int nJdata_tmp = 0;
    //    for (int l=0; l<10; l++) {
    //        for (int k=0; k<10; k++) {
    //            if (k == l) {
    //                colVals[nJdata_tmp] = l; 
    //                nJdata_tmp = nJdata_tmp + 1; 
    //            } else {
    //                if(J[10*k + l] != 0.0) {
    //                    colVals[nJdata_tmp] = k; 
    //                    nJdata_tmp = nJdata_tmp + 1; 
    //                }
    //            }
    //        }
    //        rowPtr[l+1] = nJdata_tmp;
    //    }
    //}

    return;
}

/* PROD RATE STUFF */
#if !defined(AMREX_USE_CUDA) && !defined(AMREX_USE_HIP)
static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[21];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[21];
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

    double qdot, q_f[21], q_r[21];
    comp_qfqr(q_f, q_r, sc, tc, invT);

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

void comp_k_f(double *  tc, double invT, double *  k_f)
{
    for (int i=0; i<21; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double *  tc, double invT, double *  Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[9];
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
    double refC = 101325 / 8.31446 * invT;
    double refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refC;
    Kc[2] *= refC;
    Kc[3] *= refCinv;
    Kc[4] *= refCinv;
    Kc[5] *= refCinv;

    return;
}

void comp_qfqr(double *  qf, double *  qr, double *  sc, double *  tc, double invT)
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

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 9; ++i) {
        mixture += sc[i];
    }

    double Corr[21];
    for (int i = 0; i < 21; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        double alpha[2];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[2] + (TB[0][2] - 1)*sc[1];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[2];
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

    /* simple three-body correction */
    {
        double alpha;
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

    double q_f[21], q_r[21];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 21; ++i) {
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

/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}

/*Returns the rate of progress for each reaction */
void CKQC(double *  T, double *  C, double *  qdot)
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
void CKQYP(double *  P, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    double imw[9];

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
void CKQXP(double *  P, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

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
void CKQYR(double *  rho, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double imw[9];

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
void CKQXR(double *  rho, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
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
void aJacobian(double *  J, double *  sc, double T, int consP)
{
    for (int i=0; i<100; i++) {
        J[i] = 0.0;
    }

    double wdot[9];
    for (int k=0; k<9; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31446 / T;
    double refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 9; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[9];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[9];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[9];
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

    double c_R[9], dcRdT[9], e_RT[9];
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
    for (int k = 0; k < 9; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[90+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
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
#endif
