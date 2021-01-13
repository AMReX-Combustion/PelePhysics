#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/* PURE CPU stuff  */
#ifndef AMREX_USE_GPU
namespace thermo
{
    amrex::Real fwd_A[29], fwd_beta[29], fwd_Ea[29];
    amrex::Real low_A[29], low_beta[29], low_Ea[29];
    amrex::Real rev_A[29], rev_beta[29], rev_Ea[29];
    amrex::Real troe_a[29],troe_Ts[29], troe_Tss[29], troe_Tsss[29];
    amrex::Real sri_a[29], sri_b[29], sri_c[29], sri_d[29], sri_e[29];
    amrex::Real activation_units[29], prefactor_units[29], phase_units[29];
    int is_PD[29], troe_len[29], sri_len[29], nTB[29], *TBid[29];
    amrex::Real *TB[29];
};

using namespace thermo;


/* Vectorized stuff  */

/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int *  np, amrex::Real *  y,  amrex::Real *  x)
{
    amrex::Real YOW[*np];
    amrex::Real imw[12];

    get_imw(imw);

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

/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tc[5], h[12];
    amrex::Real imw[12];

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
        hms[9*(*np)+i] = h[9];
        hms[10*(*np)+i] = h[10];
        hms[11*(*np)+i] = h[11];
    }

    for (int n=0; n<12; n++) {
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
    amrex::Real c[12*(*np)]; /*temporary storage */
    amrex::Real imw[12];

    get_imw(imw);

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


/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW[*np];
    amrex::Real imw[12];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<12; n++) {
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
    amrex::Real k_f_s[29*npt], Kc_s[29*npt], mixture[npt], g_RT[12*npt];
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

void vcomp_k_f(int npt, amrex::Real *  k_f_s, amrex::Real *  tc, amrex::Real *  invT)
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

void vcomp_gibbs(int npt, amrex::Real *  g_RT, amrex::Real *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        amrex::Real tg[5], g[12];
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

void vcomp_Kc(int npt, amrex::Real *  Kc_s, amrex::Real *  g_RT, amrex::Real *  invT)
{
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        amrex::Real refC = (101325. / 8.31451) * invT[i];
        amrex::Real refCinv = 1.0 / refC;

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

void vcomp_wdot(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;
        amrex::Real redP, F;
        amrex::Real logPred;
        amrex::Real logFcent, troe_c, troe_n, troe, F_troe;

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

static amrex::Real T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static amrex::Real k_f_save[29];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static amrex::Real Kc_save[29];
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

    amrex::Real qdot, q_f[29], q_r[29];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

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

void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f)
{
    for (int i=0; i<29; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(amrex::Real *  tc, amrex::Real invT, amrex::Real *  Kc)
{
    /*compute the Gibbs free energy */
    amrex::Real g_RT[12];
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

    for (int i=0; i<29; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 * invT;
    amrex::Real refCinv = 1 / refC;

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

void comp_qfqr(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT)
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

    amrex::Real T = tc[1];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 12; ++i) {
        mixture += sc[i];
    }

    amrex::Real Corr[29];
    for (int i = 0; i < 29; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        amrex::Real alpha[2];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[4] + (TB[0][2] - 1)*sc[1] + (TB[0][3] - 1)*sc[8] + (TB[0][4] - 1)*sc[9];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[4] + (TB[1][2] - 1)*sc[8] + (TB[1][3] - 1)*sc[9];
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

    /* Lindemann */
    {
        amrex::Real alpha;
        alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[4] + (TB[2][2] - 1)*sc[8] + (TB[2][3] - 1)*sc[9];
        amrex::Real redP = alpha / k_f_save[2] * phase_units[2] * low_A[2] * exp(low_beta[2] * tc[0] - activation_units[2] * low_Ea[2] * invT);
        Corr[2] = redP / (1. + redP);
    }

    /* simple three-body correction */
    {
        amrex::Real alpha;
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

    amrex::Real q_f[29], q_r[29];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 29; ++i) {
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
    amrex::Real c[12]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

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


/*Returns the rate of progress for each reaction */
void CKQC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  qdot)
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
/*Given P, T, and mass fractions */
void CKQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[12]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[12];

    get_imw(imw);

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
void CKQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[12]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

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
void CKQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[12]; /*temporary storage */
    amrex::Real imw[12];

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
void CKQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[12]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
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

/*compute the reaction Jacobian on CPU */
void aJacobian(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP)
{
    for (int i=0; i<169; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[12];
    for (int k=0; k<12; k++) {
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
    for (int k = 0; k < 12; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[12];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[12];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[12];
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
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[6]) + 1.000000);
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
    phi_r = pow(sc[3], 2.000000);
    Kc = refC * exp(-g_RT[3] - g_RT[3] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2.000000*h_RT[3]) - 1.000000);
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
        dqdci =  - k_r*2.000000*sc[3];
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
        dqdc[3] = dcdc_fac - k_r*2.000000*sc[3];
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
    dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[9]) + 1.000000);
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
    phi_r = pow(sc[5], 2.000000);
    Kc = refC * exp(g_RT[0] - g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0]) + (2.000000*h_RT[5]) - 1.000000);
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
        dqdci =  - k_r*2.000000*sc[5];
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
        dqdc[5] = q_nocor - k_r*2.000000*sc[5];
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
    phi_f = pow(sc[2], 2.000000);
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[1];
    Kc = refCinv * exp(-g_RT[1] + g_RT[2] + g_RT[2]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2]) + (h_RT[1]) + 1.000000);
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
        dqdci =  + k_f*2.000000*sc[2];
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
        dqdc[2] = q_nocor + k_f*2.000000*sc[2];
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
    dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[3]) + 1.000000);
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
    dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[4]) + 1.000000);
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
    dlnKcdT = invT * (-(h_RT[10]) + (h_RT[5] + h_RT[8]) - 1.000000);
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
    phi_r = pow(sc[3], 2.000000);
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[3] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2.000000*h_RT[3]));
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
    dqdci =  - k_r*2.000000*sc[3];
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
    phi_r = pow(sc[3], 2.000000);
    Kc = exp(-g_RT[3] - g_RT[3] + g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (2.000000*h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += 2 * q; /* OH */
    wdot[5] -= q; /* H */
    wdot[6] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[3];
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
    phi_f = pow(sc[6], 2.000000);
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
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
    J[14] += dqdci;               /* dwdot[O2]/d[O2] */
    J[19] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[20] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
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
    phi_f = pow(sc[6], 2.000000);
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
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
    J[14] += dqdci;               /* dwdot[O2]/d[O2] */
    J[19] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[20] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
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

    amrex::Real c_R[12], dcRdT[12], e_RT[12];
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
    for (int k = 0; k < 12; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[156+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
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

#endif
/* Transport function declarations  */


void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 50;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 3156;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 12;}


void egtransetNLITE(int* NLITE ) {
    *NLITE = 2;}


/*Patm in ergs/cm3 */
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(amrex::Real* WT ) {
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
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 3.80000000E+01;
    EPS[1] = 1.07400000E+02;
    EPS[2] = 8.00000000E+01;
    EPS[3] = 8.00000000E+01;
    EPS[4] = 5.72400000E+02;
    EPS[5] = 1.45000000E+02;
    EPS[6] = 1.07400000E+02;
    EPS[7] = 1.07400000E+02;
    EPS[8] = 9.81000000E+01;
    EPS[9] = 2.44000000E+02;
    EPS[10] = 4.98000000E+02;
    EPS[11] = 9.75300000E+01;
}


/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 2.92000000E+00;
    SIG[1] = 3.45800000E+00;
    SIG[2] = 2.75000000E+00;
    SIG[3] = 2.75000000E+00;
    SIG[4] = 2.60500000E+00;
    SIG[5] = 2.05000000E+00;
    SIG[6] = 3.45800000E+00;
    SIG[7] = 3.45800000E+00;
    SIG[8] = 3.65000000E+00;
    SIG[9] = 3.76300000E+00;
    SIG[10] = 3.59000000E+00;
    SIG[11] = 3.62100000E+00;
}


/*the dipole moment in Debye */
void egtransetDIP(amrex::Real* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[2] = 0.00000000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[4] = 1.84400000E+00;
    DIP[5] = 0.00000000E+00;
    DIP[6] = 0.00000000E+00;
    DIP[7] = 0.00000000E+00;
    DIP[8] = 0.00000000E+00;
    DIP[9] = 0.00000000E+00;
    DIP[10] = 0.00000000E+00;
    DIP[11] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 7.90000000E-01;
    POL[1] = 1.60000000E+00;
    POL[2] = 0.00000000E+00;
    POL[3] = 0.00000000E+00;
    POL[4] = 0.00000000E+00;
    POL[5] = 0.00000000E+00;
    POL[6] = 0.00000000E+00;
    POL[7] = 0.00000000E+00;
    POL[8] = 1.95000000E+00;
    POL[9] = 2.65000000E+00;
    POL[10] = 0.00000000E+00;
    POL[11] = 1.76000000E+00;
}


/*the rotational relaxation collision number at 298 K */
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 2.80000000E+02;
    ZROT[1] = 3.80000000E+00;
    ZROT[2] = 0.00000000E+00;
    ZROT[3] = 0.00000000E+00;
    ZROT[4] = 4.00000000E+00;
    ZROT[5] = 0.00000000E+00;
    ZROT[6] = 1.00000000E+00;
    ZROT[7] = 3.80000000E+00;
    ZROT[8] = 1.80000000E+00;
    ZROT[9] = 2.10000000E+00;
    ZROT[10] = 0.00000000E+00;
    ZROT[11] = 4.00000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 1;
    NLIN[1] = 1;
    NLIN[2] = 0;
    NLIN[3] = 1;
    NLIN[4] = 2;
    NLIN[5] = 0;
    NLIN[6] = 2;
    NLIN[7] = 2;
    NLIN[8] = 1;
    NLIN[9] = 1;
    NLIN[10] = 2;
    NLIN[11] = 1;
}


/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(amrex::Real* COFETA) {
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
void egtransetCOFLAM(amrex::Real* COFLAM) {
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
void egtransetCOFD(amrex::Real* COFD) {
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
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 0;
    KTDIF[1] = 5;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(amrex::Real* COFTD) {
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

#ifndef AMREX_USE_GPU
/* Initializes parameter database */
void CKINIT()
{

    // (0):  H + O2 <=> O + OH
    // (0):  H + O2 <=> O + OH
    fwd_A[8]     = 3547000000000000;
    fwd_beta[8]  = -0.40600000000000003;
    fwd_Ea[8]    = 16599;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = pow(10,-12.000000);
    is_PD[8] = 0;
    nTB[8] = 0;

    // (1):  O + H2 <=> H + OH
    // (1):  O + H2 <=> H + OH
    fwd_A[9]     = 50800;
    fwd_beta[9]  = 2.6699999999999999;
    fwd_Ea[9]    = 6290;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = pow(10,-12.000000);
    is_PD[9] = 0;
    nTB[9] = 0;

    // (2):  H2 + OH <=> H2O + H
    // (2):  H2 + OH <=> H2O + H
    fwd_A[10]     = 216000000;
    fwd_beta[10]  = 1.51;
    fwd_Ea[10]    = 3430;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = pow(10,-12.000000);
    is_PD[10] = 0;
    nTB[10] = 0;

    // (3):  O + H2O <=> OH + OH
    // (3):  O + H2O <=> OH + OH
    fwd_A[11]     = 2970000;
    fwd_beta[11]  = 2.02;
    fwd_Ea[11]    = 13400;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = pow(10,-12.000000);
    is_PD[11] = 0;
    nTB[11] = 0;

    // (4):  H2 + M <=> H + H + M
    // (4):  H2 + M <=> H + H + M
    fwd_A[3]     = 4.577e+19;
    fwd_beta[3]  = -1.3999999999999999;
    fwd_Ea[3]    = 104380;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-6.000000);
    is_PD[3] = 0;
    nTB[3] = 4;
    TB[3] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[3] = (int *) malloc(4 * sizeof(int));
    TBid[3][0] = 0; TB[3][0] = 2.5; // H2
    TBid[3][1] = 4; TB[3][1] = 12; // H2O
    TBid[3][2] = 8; TB[3][2] = 1.8999999999999999; // CO
    TBid[3][3] = 9; TB[3][3] = 3.7999999999999998; // CO2

    // (5):  O + O + M <=> O2 + M
    // (5):  O + O + M <=> O2 + M
    fwd_A[4]     = 6165000000000000;
    fwd_beta[4]  = -0.5;
    fwd_Ea[4]    = 0;
    prefactor_units[4]  = 1.0000000000000002e-12;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 0;
    nTB[4] = 4;
    TB[4] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[4] = (int *) malloc(4 * sizeof(int));
    TBid[4][0] = 0; TB[4][0] = 2.5; // H2
    TBid[4][1] = 4; TB[4][1] = 12; // H2O
    TBid[4][2] = 8; TB[4][2] = 1.8999999999999999; // CO
    TBid[4][3] = 9; TB[4][3] = 3.7999999999999998; // CO2

    // (6):  O + H + M <=> OH + M
    // (6):  O + H + M <=> OH + M
    fwd_A[5]     = 4.714e+18;
    fwd_beta[5]  = -1;
    fwd_Ea[5]    = 0;
    prefactor_units[5]  = 1.0000000000000002e-12;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-12.000000);
    is_PD[5] = 0;
    nTB[5] = 4;
    TB[5] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[5] = (int *) malloc(4 * sizeof(int));
    TBid[5][0] = 0; TB[5][0] = 2.5; // H2
    TBid[5][1] = 4; TB[5][1] = 12; // H2O
    TBid[5][2] = 8; TB[5][2] = 1.8999999999999999; // CO
    TBid[5][3] = 9; TB[5][3] = 3.7999999999999998; // CO2

    // (7):  H + OH + M <=> H2O + M
    // (7):  H + OH + M <=> H2O + M
    fwd_A[6]     = 3.8000000000000004e+22;
    fwd_beta[6]  = -2;
    fwd_Ea[6]    = 0;
    prefactor_units[6]  = 1.0000000000000002e-12;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = pow(10,-12.000000);
    is_PD[6] = 0;
    nTB[6] = 4;
    TB[6] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[6] = (int *) malloc(4 * sizeof(int));
    TBid[6][0] = 0; TB[6][0] = 2.5; // H2
    TBid[6][1] = 4; TB[6][1] = 12; // H2O
    TBid[6][2] = 8; TB[6][2] = 1.8999999999999999; // CO
    TBid[6][3] = 9; TB[6][3] = 3.7999999999999998; // CO2

    // (8):  H + O2 (+M) <=> HO2 (+M)
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
    nTB[0] = 5;
    TB[0] = (amrex::Real *) malloc(5 * sizeof(amrex::Real));
    TBid[0] = (int *) malloc(5 * sizeof(int));
    TBid[0][0] = 0; TB[0][0] = 2; // H2
    TBid[0][1] = 4; TB[0][1] = 11; // H2O
    TBid[0][2] = 1; TB[0][2] = 0.78000000000000003; // O2
    TBid[0][3] = 8; TB[0][3] = 1.8999999999999999; // CO
    TBid[0][4] = 9; TB[0][4] = 3.7999999999999998; // CO2

    // (9):  HO2 + H <=> H2 + O2
    // (9):  HO2 + H <=> H2 + O2
    fwd_A[12]     = 16600000000000;
    fwd_beta[12]  = 0;
    fwd_Ea[12]    = 823;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = pow(10,-12.000000);
    is_PD[12] = 0;
    nTB[12] = 0;

    // (10):  HO2 + H <=> OH + OH
    // (10):  HO2 + H <=> OH + OH
    fwd_A[13]     = 70790000000000;
    fwd_beta[13]  = 0;
    fwd_Ea[13]    = 295;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = pow(10,-12.000000);
    is_PD[13] = 0;
    nTB[13] = 0;

    // (11):  HO2 + O <=> O2 + OH
    // (11):  HO2 + O <=> O2 + OH
    fwd_A[14]     = 32500000000000;
    fwd_beta[14]  = 0;
    fwd_Ea[14]    = 0;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = pow(10,-12.000000);
    is_PD[14] = 0;
    nTB[14] = 0;

    // (12):  HO2 + OH <=> H2O + O2
    // (12):  HO2 + OH <=> H2O + O2
    fwd_A[15]     = 28900000000000;
    fwd_beta[15]  = 0;
    fwd_Ea[15]    = -497;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = pow(10,-12.000000);
    is_PD[15] = 0;
    nTB[15] = 0;

    // (13):  HO2 + HO2 <=> H2O2 + O2
    // (13):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[16]     = 420000000000000;
    fwd_beta[16]  = 0;
    fwd_Ea[16]    = 11982;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = pow(10,-12.000000);
    is_PD[16] = 0;
    nTB[16] = 0;

    // (14):  HO2 + HO2 <=> H2O2 + O2
    // (14):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[17]     = 130000000000;
    fwd_beta[17]  = 0;
    fwd_Ea[17]    = -1629.3;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = pow(10,-12.000000);
    is_PD[17] = 0;
    nTB[17] = 0;

    // (15):  H2O2 (+M) <=> OH + OH (+M)
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
    nTB[1] = 4;
    TB[1] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[1] = (int *) malloc(4 * sizeof(int));
    TBid[1][0] = 0; TB[1][0] = 2.5; // H2
    TBid[1][1] = 4; TB[1][1] = 12; // H2O
    TBid[1][2] = 8; TB[1][2] = 1.8999999999999999; // CO
    TBid[1][3] = 9; TB[1][3] = 3.7999999999999998; // CO2

    // (16):  H2O2 + H <=> H2O + OH
    // (16):  H2O2 + H <=> H2O + OH
    fwd_A[18]     = 24100000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 3970;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = pow(10,-12.000000);
    is_PD[18] = 0;
    nTB[18] = 0;

    // (17):  H2O2 + H <=> HO2 + H2
    // (17):  H2O2 + H <=> HO2 + H2
    fwd_A[19]     = 48200000000000;
    fwd_beta[19]  = 0;
    fwd_Ea[19]    = 7950;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = pow(10,-12.000000);
    is_PD[19] = 0;
    nTB[19] = 0;

    // (18):  H2O2 + O <=> OH + HO2
    // (18):  H2O2 + O <=> OH + HO2
    fwd_A[20]     = 9550000;
    fwd_beta[20]  = 2;
    fwd_Ea[20]    = 3970;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-12.000000);
    is_PD[20] = 0;
    nTB[20] = 0;

    // (19):  H2O2 + OH <=> HO2 + H2O
    // (19):  H2O2 + OH <=> HO2 + H2O
    fwd_A[21]     = 1000000000000;
    fwd_beta[21]  = 0;
    fwd_Ea[21]    = 0;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = pow(10,-12.000000);
    is_PD[21] = 0;
    nTB[21] = 0;

    // (20):  H2O2 + OH <=> HO2 + H2O
    // (20):  H2O2 + OH <=> HO2 + H2O
    fwd_A[22]     = 580000000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 9557;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = pow(10,-12.000000);
    is_PD[22] = 0;
    nTB[22] = 0;

    // (21):  CO + O (+M) <=> CO2 (+M)
    // (21):  CO + O (+M) <=> CO2 (+M)
    fwd_A[2]     = 18000000000;
    fwd_beta[2]  = 0;
    fwd_Ea[2]    = 2384;
    low_A[2]     = 1.5500000000000001e+24;
    low_beta[2]  = -2.79;
    low_Ea[2]    = 4191;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-12.000000);
    is_PD[2] = 1;
    nTB[2] = 4;
    TB[2] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[2] = (int *) malloc(4 * sizeof(int));
    TBid[2][0] = 0; TB[2][0] = 2.5; // H2
    TBid[2][1] = 4; TB[2][1] = 12; // H2O
    TBid[2][2] = 8; TB[2][2] = 1.8999999999999999; // CO
    TBid[2][3] = 9; TB[2][3] = 3.7999999999999998; // CO2

    // (22):  CO + O2 <=> CO2 + O
    // (22):  CO + O2 <=> CO2 + O
    fwd_A[23]     = 2530000000000;
    fwd_beta[23]  = 0;
    fwd_Ea[23]    = 47700;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = pow(10,-12.000000);
    is_PD[23] = 0;
    nTB[23] = 0;

    // (23):  CO + HO2 <=> CO2 + OH
    // (23):  CO + HO2 <=> CO2 + OH
    fwd_A[24]     = 30100000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 23000;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = pow(10,-12.000000);
    is_PD[24] = 0;
    nTB[24] = 0;

    // (24):  CO + OH <=> CO2 + H
    // (24):  CO + OH <=> CO2 + H
    fwd_A[25]     = 222900;
    fwd_beta[25]  = 1.8899999999999999;
    fwd_Ea[25]    = -1158.7;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = pow(10,-12.000000);
    is_PD[25] = 0;
    nTB[25] = 0;

    // (25):  HCO + M <=> H + CO + M
    // (25):  HCO + M <=> H + CO + M
    fwd_A[7]     = 474850000000;
    fwd_beta[7]  = 0.65900000000000003;
    fwd_Ea[7]    = 14874;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = pow(10,-6.000000);
    is_PD[7] = 0;
    nTB[7] = 4;
    TB[7] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[7] = (int *) malloc(4 * sizeof(int));
    TBid[7][0] = 0; TB[7][0] = 2.5; // H2
    TBid[7][1] = 4; TB[7][1] = 6; // H2O
    TBid[7][2] = 8; TB[7][2] = 1.8999999999999999; // CO
    TBid[7][3] = 9; TB[7][3] = 3.7999999999999998; // CO2

    // (26):  HCO + O2 <=> CO + HO2
    // (26):  HCO + O2 <=> CO + HO2
    fwd_A[26]     = 7580000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 410;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = pow(10,-12.000000);
    is_PD[26] = 0;
    nTB[26] = 0;

    // (27):  HCO + H <=> CO + H2
    // (27):  HCO + H <=> CO + H2
    fwd_A[27]     = 72300000000000;
    fwd_beta[27]  = 0;
    fwd_Ea[27]    = 0;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = pow(10,-12.000000);
    is_PD[27] = 0;
    nTB[27] = 0;

    // (28):  HCO + O <=> CO2 + H
    // (28):  HCO + O <=> CO2 + H
    fwd_A[28]     = 30000000000000;
    fwd_beta[28]  = 0;
    fwd_Ea[28]    = 0;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = pow(10,-12.000000);
    is_PD[28] = 0;
    nTB[28] = 0;

}


/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<29; ++i) {
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
    awt[1] = 1.007970; /*H */
    awt[2] = 15.999400; /*O */
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


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "C";
    ename[1] = "H";
    ename[2] = "O";
    ename[3] = "N";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(12);
    kname[0] = "H2";
    kname[1] = "O2";
    kname[2] = "O";
    kname[3] = "OH";
    kname[4] = "H2O";
    kname[5] = "H";
    kname[6] = "HO2";
    kname[7] = "H2O2";
    kname[8] = "CO";
    kname[9] = "CO2";
    kname[10] = "HCO";
    kname[11] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(169);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(12);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[169];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<12; l++) {
                c_d[l] = 1.0/ 12.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<13; k++) {
        for (int l=0; l<13; l++) {
            if(J_h[ 13 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(169);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(12);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[169];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<12; k++) {
                c_d[k] = 1.0/ 12.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<13; k++) {
        for (int l=0; l<13; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 13 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(169);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(12);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[169];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<12; k++) {
                c_d[k] = 1.0/ 12.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<13; k++) {
        for (int l=0; l<13; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 13 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(169);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(12);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[169];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<12; k++) {
                c_d[k] = 1.0/ 12.000000 ;
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
        offset_row = nc * 13;
        offset_col = nc * 13;
        for (int k=0; k<13; k++) {
            for (int l=0; l<13; l++) {
                if(J_h[13*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(169);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(12);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[169];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<12; k++) {
                c_d[k] = 1.0/ 12.000000 ;
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
            offset = nc * 13;
            for (int l=0; l<13; l++) {
                for (int k=0; k<13; k++) {
                    if(J_h[13*k + l] != 0.0) {
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
            offset = nc * 13;
            for (int l=0; l<13; l++) {
                for (int k=0; k<13; k++) {
                    if(J_h[13*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(169);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(12);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[169];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<12; k++) {
                c_d[k] = 1.0/ 12.000000 ;
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
            offset = nc * 13;
            for (int l=0; l<13; l++) {
                for (int k=0; k<13; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[13*k + l] != 0.0) {
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
            offset = nc * 13;
            for (int l=0; l<13; l++) {
                for (int k=0; k<13; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[13*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(169);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(12);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[169];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<12; k++) {
                c_d[k] = 1.0/ 12.000000 ;
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
    for (int k=0; k<13; k++) {
        for (int l=0; l<13; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 13*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[13*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 13*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(169);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(12);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[169];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<12; k++) {
                c_d[k] = 1.0/ 12.000000 ;
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
        for (int l=0; l<13; l++) {
            for (int k=0; k<13; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[13*k + l] != 0.0) {
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
        for (int l=0; l<13; l++) {
            for (int k=0; k<13; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[13*k + l] != 0.0) {
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

