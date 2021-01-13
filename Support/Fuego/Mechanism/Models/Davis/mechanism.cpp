#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/* PURE CPU stuff  */
#ifndef AMREX_USE_GPU
namespace thermo
{
    amrex::Real fwd_A[38], fwd_beta[38], fwd_Ea[38];
    amrex::Real low_A[38], low_beta[38], low_Ea[38];
    amrex::Real rev_A[38], rev_beta[38], rev_Ea[38];
    amrex::Real troe_a[38],troe_Ts[38], troe_Tss[38], troe_Tsss[38];
    amrex::Real sri_a[38], sri_b[38], sri_c[38], sri_d[38], sri_e[38];
    amrex::Real activation_units[38], prefactor_units[38], phase_units[38];
    int is_PD[38], troe_len[38], sri_len[38], nTB[38], *TBid[38];
    amrex::Real *TB[38];
};

using namespace thermo;


/* Vectorized stuff  */

/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int *  np, amrex::Real *  y,  amrex::Real *  x)
{
    amrex::Real YOW[*np];
    amrex::Real imw[14];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<14; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<14; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}

/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tc[5], h[14];
    amrex::Real imw[14];

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
        hms[12*(*np)+i] = h[12];
        hms[13*(*np)+i] = h[13];
    }

    for (int n=0; n<14; n++) {
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
    amrex::Real c[14*(*np)]; /*temporary storage */
    amrex::Real imw[14];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<14; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<14*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW[*np];
    amrex::Real imw[14];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<14; n++) {
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
    amrex::Real k_f_s[38*npt], Kc_s[38*npt], mixture[npt], g_RT[14*npt];
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

    for (int n=0; n<14; n++) {
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
        k_f_s[29*npt+i] = prefactor_units[29] * fwd_A[29] * exp(fwd_beta[29] * tc[i] - activation_units[29] * fwd_Ea[29] * invT[i]);
        k_f_s[30*npt+i] = prefactor_units[30] * fwd_A[30] * exp(fwd_beta[30] * tc[i] - activation_units[30] * fwd_Ea[30] * invT[i]);
        k_f_s[31*npt+i] = prefactor_units[31] * fwd_A[31] * exp(fwd_beta[31] * tc[i] - activation_units[31] * fwd_Ea[31] * invT[i]);
        k_f_s[32*npt+i] = prefactor_units[32] * fwd_A[32] * exp(fwd_beta[32] * tc[i] - activation_units[32] * fwd_Ea[32] * invT[i]);
        k_f_s[33*npt+i] = prefactor_units[33] * fwd_A[33] * exp(fwd_beta[33] * tc[i] - activation_units[33] * fwd_Ea[33] * invT[i]);
        k_f_s[34*npt+i] = prefactor_units[34] * fwd_A[34] * exp(fwd_beta[34] * tc[i] - activation_units[34] * fwd_Ea[34] * invT[i]);
        k_f_s[35*npt+i] = prefactor_units[35] * fwd_A[35] * exp(fwd_beta[35] * tc[i] - activation_units[35] * fwd_Ea[35] * invT[i]);
        k_f_s[36*npt+i] = prefactor_units[36] * fwd_A[36] * exp(fwd_beta[36] * tc[i] - activation_units[36] * fwd_Ea[36] * invT[i]);
        k_f_s[37*npt+i] = prefactor_units[37] * fwd_A[37] * exp(fwd_beta[37] * tc[i] - activation_units[37] * fwd_Ea[37] * invT[i]);
    }
}

void vcomp_gibbs(int npt, amrex::Real *  g_RT, amrex::Real *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        amrex::Real tg[5], g[14];
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
    }
}

void vcomp_Kc(int npt, amrex::Real *  Kc_s, amrex::Real *  g_RT, amrex::Real *  invT)
{
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        amrex::Real refC = (101325. / 8.31451) * invT[i];
        amrex::Real refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[11*npt+i]) - (g_RT[8*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[12*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[10*npt+i]) - (g_RT[13*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[1*npt+i]) - (g_RT[0*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (g_RT[9*npt+i]));
        Kc_s[5*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[6*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[5*npt+i]) - (g_RT[11*npt+i]));
        Kc_s[7*npt+i] = refC * exp((g_RT[7*npt+i]) - (g_RT[1*npt+i] + g_RT[10*npt+i]));
        Kc_s[8*npt+i] = exp((g_RT[1*npt+i] + g_RT[11*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[9*npt+i] = exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[1*npt+i] + g_RT[6*npt+i]));
        Kc_s[10*npt+i] = exp((g_RT[0*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i]));
        Kc_s[11*npt+i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
        Kc_s[12*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[1*npt+i] + g_RT[1*npt+i]) - (g_RT[0*npt+i] + g_RT[0*npt+i]));
        Kc_s[13*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[1*npt+i] + g_RT[9*npt+i]) - (g_RT[0*npt+i] + g_RT[9*npt+i]));
        Kc_s[14*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[1*npt+i] + g_RT[13*npt+i]) - (g_RT[0*npt+i] + g_RT[13*npt+i]));
        Kc_s[15*npt+i] = exp((g_RT[0*npt+i] + g_RT[11*npt+i]) - (g_RT[1*npt+i] + g_RT[8*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[1*npt+i] + g_RT[8*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[1*npt+i] + g_RT[8*npt+i]) - (g_RT[6*npt+i] + g_RT[6*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[5*npt+i] + g_RT[8*npt+i]) - (g_RT[6*npt+i] + g_RT[11*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[6*npt+i] + g_RT[8*npt+i]) - (g_RT[9*npt+i] + g_RT[11*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[6*npt+i] + g_RT[8*npt+i]) - (g_RT[9*npt+i] + g_RT[11*npt+i]));
        Kc_s[21*npt+i] = exp((g_RT[8*npt+i] + g_RT[8*npt+i]) - (g_RT[11*npt+i] + g_RT[12*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[8*npt+i] + g_RT[8*npt+i]) - (g_RT[11*npt+i] + g_RT[12*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[1*npt+i] + g_RT[12*npt+i]) - (g_RT[0*npt+i] + g_RT[8*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[1*npt+i] + g_RT[12*npt+i]) - (g_RT[6*npt+i] + g_RT[9*npt+i]));
        Kc_s[25*npt+i] = exp((g_RT[5*npt+i] + g_RT[12*npt+i]) - (g_RT[6*npt+i] + g_RT[8*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[6*npt+i] + g_RT[12*npt+i]) - (g_RT[8*npt+i] + g_RT[9*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[6*npt+i] + g_RT[12*npt+i]) - (g_RT[8*npt+i] + g_RT[9*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[6*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[29*npt+i] = exp((g_RT[6*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[30*npt+i] = exp((g_RT[10*npt+i] + g_RT[11*npt+i]) - (g_RT[5*npt+i] + g_RT[13*npt+i]));
        Kc_s[31*npt+i] = exp((g_RT[8*npt+i] + g_RT[10*npt+i]) - (g_RT[6*npt+i] + g_RT[13*npt+i]));
        Kc_s[32*npt+i] = exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[0*npt+i] + g_RT[10*npt+i]));
        Kc_s[33*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[6*npt+i] + g_RT[10*npt+i]));
        Kc_s[34*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[35*npt+i] = exp((g_RT[6*npt+i] + g_RT[7*npt+i]) - (g_RT[9*npt+i] + g_RT[10*npt+i]));
        Kc_s[36*npt+i] = refC * exp((g_RT[7*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i] + g_RT[10*npt+i]));
        Kc_s[37*npt+i] = exp((g_RT[7*npt+i] + g_RT[11*npt+i]) - (g_RT[8*npt+i] + g_RT[10*npt+i]));
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
        phi_f = sc[1*npt+i]*sc[11*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[11*npt+i] + (TB[0][1] - 1)*sc[9*npt+i] + (TB[0][2] - 1)*sc[10*npt+i] + (TB[0][3] - 1)*sc[13*npt+i] + (TB[0][4] - 1)*sc[2*npt+i] + (TB[0][5] - 1)*sc[4*npt+i] + (TB[0][6] - 1)*sc[0*npt+i];
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
        phi_r = sc[8*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[0*npt+i] + (TB[1][1] - 1)*sc[9*npt+i] + (TB[1][2] - 1)*sc[10*npt+i] + (TB[1][3] - 1)*sc[13*npt+i] + (TB[1][4] - 1)*sc[2*npt+i] + (TB[1][5] - 1)*sc[4*npt+i];
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
        phi_r = sc[12*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 3: CO + O (+M) <=> CO2 (+M) */
        phi_f = sc[5*npt+i]*sc[10*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[0*npt+i] + (TB[2][1] - 1)*sc[9*npt+i] + (TB[2][2] - 1)*sc[10*npt+i] + (TB[2][3] - 1)*sc[13*npt+i] + (TB[2][4] - 1)*sc[2*npt+i] + (TB[2][5] - 1)*sc[4*npt+i];
        k_f = k_f_s[2*npt+i];
        redP = alpha / k_f * phase_units[2] * low_A[2] * exp(low_beta[2] * tc[i] - activation_units[2] * low_Ea[2] * invT[i]);
        F = redP / (1 + redP);
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 4: H + H + M <=> H2 + M */
        phi_f = sc[1*npt+i]*sc[1*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[0*npt+i] + (TB[3][1] - 1)*sc[9*npt+i] + (TB[3][2] - 1)*sc[13*npt+i] + (TB[3][3] - 1)*sc[2*npt+i] + (TB[3][4] - 1)*sc[4*npt+i];
        k_f = alpha * k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[1*npt+i] -= qdot;

        /*reaction 5: H + OH + M <=> H2O + M */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[0*npt+i] + (TB[4][1] - 1)*sc[9*npt+i] + (TB[4][2] - 1)*sc[10*npt+i] + (TB[4][3] - 1)*sc[13*npt+i] + (TB[4][4] - 1)*sc[2*npt+i] + (TB[4][5] - 1)*sc[4*npt+i];
        k_f = alpha * k_f_s[4*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 6: O + H + M <=> OH + M */
        phi_f = sc[1*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[0*npt+i] + (TB[5][1] - 1)*sc[9*npt+i] + (TB[5][2] - 1)*sc[10*npt+i] + (TB[5][3] - 1)*sc[13*npt+i] + (TB[5][4] - 1)*sc[2*npt+i] + (TB[5][5] - 1)*sc[4*npt+i];
        k_f = alpha * k_f_s[5*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 7: O + O + M <=> O2 + M */
        phi_f = sc[5*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[6][0] - 1)*sc[0*npt+i] + (TB[6][1] - 1)*sc[9*npt+i] + (TB[6][2] - 1)*sc[10*npt+i] + (TB[6][3] - 1)*sc[13*npt+i] + (TB[6][4] - 1)*sc[2*npt+i] + (TB[6][5] - 1)*sc[4*npt+i];
        k_f = alpha * k_f_s[6*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 8: HCO + M <=> CO + H + M */
        phi_f = sc[7*npt+i];
        alpha = mixture[i] + (TB[7][0] - 1)*sc[0*npt+i] + (TB[7][1] - 1)*sc[9*npt+i] + (TB[7][2] - 1)*sc[10*npt+i] + (TB[7][3] - 1)*sc[13*npt+i];
        k_f = alpha * k_f_s[7*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[10*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 9: H + O2 <=> O + OH */
        phi_f = sc[1*npt+i]*sc[11*npt+i];
        k_f = k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 10: O + H2 <=> H + OH */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        k_f = k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[6*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 11: OH + H2 <=> H + H2O */
        phi_f = sc[0*npt+i]*sc[6*npt+i];
        k_f = k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 12: OH + OH <=> O + H2O */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[9*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 13: H + H + H2 <=> H2 + H2 */
        phi_f = sc[0*npt+i]*sc[1*npt+i]*sc[1*npt+i];
        k_f = k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[0*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] += qdot;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[1*npt+i] -= qdot;

        /*reaction 14: H + H + H2O <=> H2 + H2O */
        phi_f = sc[1*npt+i]*sc[1*npt+i]*sc[9*npt+i];
        k_f = k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[9*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[1*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 15: H + H + CO2 <=> H2 + CO2 */
        phi_f = sc[1*npt+i]*sc[1*npt+i]*sc[13*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[13*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[1*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 16: H2 + O2 <=> HO2 + H */
        phi_f = sc[0*npt+i]*sc[11*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[8*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 17: HO2 + H <=> O + H2O */
        phi_f = sc[1*npt+i]*sc[8*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[9*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 18: HO2 + H <=> OH + OH */
        phi_f = sc[1*npt+i]*sc[8*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[6*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 19: HO2 + O <=> OH + O2 */
        phi_f = sc[5*npt+i]*sc[8*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[11*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 20: HO2 + OH <=> O2 + H2O */
        phi_f = sc[6*npt+i]*sc[8*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[11*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 21: HO2 + OH <=> O2 + H2O */
        phi_f = sc[6*npt+i]*sc[8*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[11*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
        phi_f = sc[8*npt+i]*sc[8*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[12*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
        phi_f = sc[8*npt+i]*sc[8*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[12*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 24: H2O2 + H <=> HO2 + H2 */
        phi_f = sc[1*npt+i]*sc[12*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[8*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 25: H2O2 + H <=> OH + H2O */
        phi_f = sc[1*npt+i]*sc[12*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[9*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 26: H2O2 + O <=> OH + HO2 */
        phi_f = sc[5*npt+i]*sc[12*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[8*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 27: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[6*npt+i]*sc[12*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[9*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 28: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[6*npt+i]*sc[12*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[9*npt+i];
        Kc = Kc_s[27*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 29: CO + OH <=> CO2 + H */
        phi_f = sc[6*npt+i]*sc[10*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[13*npt+i];
        Kc = Kc_s[28*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 30: CO + OH <=> CO2 + H */
        phi_f = sc[6*npt+i]*sc[10*npt+i];
        k_f = k_f_s[29*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[13*npt+i];
        Kc = Kc_s[29*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 31: CO + O2 <=> CO2 + O */
        phi_f = sc[10*npt+i]*sc[11*npt+i];
        k_f = k_f_s[30*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[13*npt+i];
        Kc = Kc_s[30*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 32: CO + HO2 <=> CO2 + OH */
        phi_f = sc[8*npt+i]*sc[10*npt+i];
        k_f = k_f_s[31*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[13*npt+i];
        Kc = Kc_s[31*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 33: HCO + H <=> CO + H2 */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        k_f = k_f_s[32*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[10*npt+i];
        Kc = Kc_s[32*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 34: HCO + O <=> CO + OH */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[33*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[10*npt+i];
        Kc = Kc_s[33*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 35: HCO + O <=> CO2 + H */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[34*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[13*npt+i];
        Kc = Kc_s[34*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 36: HCO + OH <=> CO + H2O */
        phi_f = sc[6*npt+i]*sc[7*npt+i];
        k_f = k_f_s[35*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[10*npt+i];
        Kc = Kc_s[35*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 37: HCO + H2O <=> CO + H + H2O */
        phi_f = sc[7*npt+i]*sc[9*npt+i];
        k_f = k_f_s[36*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i]*sc[10*npt+i];
        Kc = Kc_s[36*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 38: HCO + O2 <=> CO + HO2 */
        phi_f = sc[7*npt+i]*sc[11*npt+i];
        k_f = k_f_s[37*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[10*npt+i];
        Kc = Kc_s[37*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
    }
}

static amrex::Real T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static amrex::Real k_f_save[38];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static amrex::Real Kc_save[38];
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

    amrex::Real qdot, q_f[38], q_r[38];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 14; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[1]-q_r[1];
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[5] -= qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[1] -= qdot;

    qdot = q_f[4]-q_r[4];
    wdot[1] -= qdot;
    wdot[6] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[1] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[5] -= qdot;
    wdot[5] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[1] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[9]-q_r[9];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[0] -= qdot;
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[1] -= qdot;

    qdot = q_f[13]-q_r[13];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[1] -= qdot;
    wdot[9] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[1] -= qdot;
    wdot[13] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[16]-q_r[16];
    wdot[1] -= qdot;
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[6] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[19]-q_r[19];
    wdot[6] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;
    wdot[11] += qdot;

    qdot = q_f[20]-q_r[20];
    wdot[6] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;
    wdot[11] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[8] -= qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;
    wdot[12] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[8] -= qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;
    wdot[12] += qdot;

    qdot = q_f[23]-q_r[23];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[9] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[8] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[6] -= qdot;
    wdot[8] += qdot;
    wdot[9] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[27]-q_r[27];
    wdot[6] -= qdot;
    wdot[8] += qdot;
    wdot[9] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[29]-q_r[29];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[30]-q_r[30];
    wdot[5] += qdot;
    wdot[10] -= qdot;
    wdot[11] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[31]-q_r[31];
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[32]-q_r[32];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[7] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[34]-q_r[34];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[7] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[6] -= qdot;
    wdot[7] -= qdot;
    wdot[9] += qdot;
    wdot[10] += qdot;

    qdot = q_f[36]-q_r[36];
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[9] -= qdot;
    wdot[9] += qdot;
    wdot[10] += qdot;

    qdot = q_f[37]-q_r[37];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    return;
}

void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f)
{
    for (int i=0; i<38; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(amrex::Real *  tc, amrex::Real invT, amrex::Real *  Kc)
{
    /*compute the Gibbs free energy */
    amrex::Real g_RT[14];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[1] - g_RT[8] + g_RT[11];
    Kc[1] = g_RT[6] + g_RT[6] - g_RT[12];
    Kc[2] = g_RT[5] + g_RT[10] - g_RT[13];
    Kc[3] = -g_RT[0] + g_RT[1] + g_RT[1];
    Kc[4] = g_RT[1] + g_RT[6] - g_RT[9];
    Kc[5] = g_RT[1] + g_RT[5] - g_RT[6];
    Kc[6] = g_RT[5] + g_RT[5] - g_RT[11];
    Kc[7] = -g_RT[1] + g_RT[7] - g_RT[10];
    Kc[8] = g_RT[1] - g_RT[5] - g_RT[6] + g_RT[11];
    Kc[9] = g_RT[0] - g_RT[1] + g_RT[5] - g_RT[6];
    Kc[10] = g_RT[0] - g_RT[1] + g_RT[6] - g_RT[9];
    Kc[11] = -g_RT[5] + g_RT[6] + g_RT[6] - g_RT[9];
    Kc[12] = g_RT[0] - g_RT[0] - g_RT[0] + g_RT[1] + g_RT[1];
    Kc[13] = -g_RT[0] + g_RT[1] + g_RT[1] + g_RT[9] - g_RT[9];
    Kc[14] = -g_RT[0] + g_RT[1] + g_RT[1] + g_RT[13] - g_RT[13];
    Kc[15] = g_RT[0] - g_RT[1] - g_RT[8] + g_RT[11];
    Kc[16] = g_RT[1] - g_RT[5] + g_RT[8] - g_RT[9];
    Kc[17] = g_RT[1] - g_RT[6] - g_RT[6] + g_RT[8];
    Kc[18] = g_RT[5] - g_RT[6] + g_RT[8] - g_RT[11];
    Kc[19] = g_RT[6] + g_RT[8] - g_RT[9] - g_RT[11];
    Kc[20] = g_RT[6] + g_RT[8] - g_RT[9] - g_RT[11];
    Kc[21] = g_RT[8] + g_RT[8] - g_RT[11] - g_RT[12];
    Kc[22] = g_RT[8] + g_RT[8] - g_RT[11] - g_RT[12];
    Kc[23] = -g_RT[0] + g_RT[1] - g_RT[8] + g_RT[12];
    Kc[24] = g_RT[1] - g_RT[6] - g_RT[9] + g_RT[12];
    Kc[25] = g_RT[5] - g_RT[6] - g_RT[8] + g_RT[12];
    Kc[26] = g_RT[6] - g_RT[8] - g_RT[9] + g_RT[12];
    Kc[27] = g_RT[6] - g_RT[8] - g_RT[9] + g_RT[12];
    Kc[28] = -g_RT[1] + g_RT[6] + g_RT[10] - g_RT[13];
    Kc[29] = -g_RT[1] + g_RT[6] + g_RT[10] - g_RT[13];
    Kc[30] = -g_RT[5] + g_RT[10] + g_RT[11] - g_RT[13];
    Kc[31] = -g_RT[6] + g_RT[8] + g_RT[10] - g_RT[13];
    Kc[32] = -g_RT[0] + g_RT[1] + g_RT[7] - g_RT[10];
    Kc[33] = g_RT[5] - g_RT[6] + g_RT[7] - g_RT[10];
    Kc[34] = -g_RT[1] + g_RT[5] + g_RT[7] - g_RT[13];
    Kc[35] = g_RT[6] + g_RT[7] - g_RT[9] - g_RT[10];
    Kc[36] = -g_RT[1] + g_RT[7] + g_RT[9] - g_RT[9] - g_RT[10];
    Kc[37] = g_RT[7] - g_RT[8] - g_RT[10] + g_RT[11];

    for (int i=0; i<38; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 * invT;
    amrex::Real refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refCinv;
    Kc[2] *= refCinv;
    Kc[3] *= refCinv;
    Kc[4] *= refCinv;
    Kc[5] *= refCinv;
    Kc[6] *= refCinv;
    Kc[7] *= refC;
    Kc[12] *= refCinv;
    Kc[13] *= refCinv;
    Kc[14] *= refCinv;
    Kc[36] *= refC;

    return;
}

void comp_qfqr(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT)
{

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[1]*sc[11];
    qr[0] = sc[8];

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    qf[1] = sc[6]*sc[6];
    qr[1] = sc[12];

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    qf[2] = sc[5]*sc[10];
    qr[2] = sc[13];

    /*reaction 4: H + H + M <=> H2 + M */
    qf[3] = sc[1]*sc[1];
    qr[3] = sc[0];

    /*reaction 5: H + OH + M <=> H2O + M */
    qf[4] = sc[1]*sc[6];
    qr[4] = sc[9];

    /*reaction 6: O + H + M <=> OH + M */
    qf[5] = sc[1]*sc[5];
    qr[5] = sc[6];

    /*reaction 7: O + O + M <=> O2 + M */
    qf[6] = sc[5]*sc[5];
    qr[6] = sc[11];

    /*reaction 8: HCO + M <=> CO + H + M */
    qf[7] = sc[7];
    qr[7] = sc[1]*sc[10];

    /*reaction 9: H + O2 <=> O + OH */
    qf[8] = sc[1]*sc[11];
    qr[8] = sc[5]*sc[6];

    /*reaction 10: O + H2 <=> H + OH */
    qf[9] = sc[0]*sc[5];
    qr[9] = sc[1]*sc[6];

    /*reaction 11: OH + H2 <=> H + H2O */
    qf[10] = sc[0]*sc[6];
    qr[10] = sc[1]*sc[9];

    /*reaction 12: OH + OH <=> O + H2O */
    qf[11] = sc[6]*sc[6];
    qr[11] = sc[5]*sc[9];

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    qf[12] = sc[0]*sc[1]*sc[1];
    qr[12] = sc[0]*sc[0];

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    qf[13] = sc[1]*sc[1]*sc[9];
    qr[13] = sc[0]*sc[9];

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    qf[14] = sc[1]*sc[1]*sc[13];
    qr[14] = sc[0]*sc[13];

    /*reaction 16: H2 + O2 <=> HO2 + H */
    qf[15] = sc[0]*sc[11];
    qr[15] = sc[1]*sc[8];

    /*reaction 17: HO2 + H <=> O + H2O */
    qf[16] = sc[1]*sc[8];
    qr[16] = sc[5]*sc[9];

    /*reaction 18: HO2 + H <=> OH + OH */
    qf[17] = sc[1]*sc[8];
    qr[17] = sc[6]*sc[6];

    /*reaction 19: HO2 + O <=> OH + O2 */
    qf[18] = sc[5]*sc[8];
    qr[18] = sc[6]*sc[11];

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    qf[19] = sc[6]*sc[8];
    qr[19] = sc[9]*sc[11];

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    qf[20] = sc[6]*sc[8];
    qr[20] = sc[9]*sc[11];

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    qf[21] = sc[8]*sc[8];
    qr[21] = sc[11]*sc[12];

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    qf[22] = sc[8]*sc[8];
    qr[22] = sc[11]*sc[12];

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    qf[23] = sc[1]*sc[12];
    qr[23] = sc[0]*sc[8];

    /*reaction 25: H2O2 + H <=> OH + H2O */
    qf[24] = sc[1]*sc[12];
    qr[24] = sc[6]*sc[9];

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    qf[25] = sc[5]*sc[12];
    qr[25] = sc[6]*sc[8];

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    qf[26] = sc[6]*sc[12];
    qr[26] = sc[8]*sc[9];

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    qf[27] = sc[6]*sc[12];
    qr[27] = sc[8]*sc[9];

    /*reaction 29: CO + OH <=> CO2 + H */
    qf[28] = sc[6]*sc[10];
    qr[28] = sc[1]*sc[13];

    /*reaction 30: CO + OH <=> CO2 + H */
    qf[29] = sc[6]*sc[10];
    qr[29] = sc[1]*sc[13];

    /*reaction 31: CO + O2 <=> CO2 + O */
    qf[30] = sc[10]*sc[11];
    qr[30] = sc[5]*sc[13];

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    qf[31] = sc[8]*sc[10];
    qr[31] = sc[6]*sc[13];

    /*reaction 33: HCO + H <=> CO + H2 */
    qf[32] = sc[1]*sc[7];
    qr[32] = sc[0]*sc[10];

    /*reaction 34: HCO + O <=> CO + OH */
    qf[33] = sc[5]*sc[7];
    qr[33] = sc[6]*sc[10];

    /*reaction 35: HCO + O <=> CO2 + H */
    qf[34] = sc[5]*sc[7];
    qr[34] = sc[1]*sc[13];

    /*reaction 36: HCO + OH <=> CO + H2O */
    qf[35] = sc[6]*sc[7];
    qr[35] = sc[9]*sc[10];

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    qf[36] = sc[7]*sc[9];
    qr[36] = sc[1]*sc[9]*sc[10];

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    qf[37] = sc[7]*sc[11];
    qr[37] = sc[8]*sc[10];

    amrex::Real T = tc[1];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 14; ++i) {
        mixture += sc[i];
    }

    amrex::Real Corr[38];
    for (int i = 0; i < 38; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        amrex::Real alpha[2];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[11] + (TB[0][1] - 1)*sc[9] + (TB[0][2] - 1)*sc[10] + (TB[0][3] - 1)*sc[13] + (TB[0][4] - 1)*sc[2] + (TB[0][5] - 1)*sc[4] + (TB[0][6] - 1)*sc[0];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[9] + (TB[1][2] - 1)*sc[10] + (TB[1][3] - 1)*sc[13] + (TB[1][4] - 1)*sc[2] + (TB[1][5] - 1)*sc[4];
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
        alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[9] + (TB[2][2] - 1)*sc[10] + (TB[2][3] - 1)*sc[13] + (TB[2][4] - 1)*sc[2] + (TB[2][5] - 1)*sc[4];
        amrex::Real redP = alpha / k_f_save[2] * phase_units[2] * low_A[2] * exp(low_beta[2] * tc[0] - activation_units[2] * low_Ea[2] * invT);
        Corr[2] = redP / (1. + redP);
    }

    /* simple three-body correction */
    {
        amrex::Real alpha;
        alpha = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[9] + (TB[3][2] - 1)*sc[13] + (TB[3][3] - 1)*sc[2] + (TB[3][4] - 1)*sc[4];
        Corr[3] = alpha;
        alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[9] + (TB[4][2] - 1)*sc[10] + (TB[4][3] - 1)*sc[13] + (TB[4][4] - 1)*sc[2] + (TB[4][5] - 1)*sc[4];
        Corr[4] = alpha;
        alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[9] + (TB[5][2] - 1)*sc[10] + (TB[5][3] - 1)*sc[13] + (TB[5][4] - 1)*sc[2] + (TB[5][5] - 1)*sc[4];
        Corr[5] = alpha;
        alpha = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[9] + (TB[6][2] - 1)*sc[10] + (TB[6][3] - 1)*sc[13] + (TB[6][4] - 1)*sc[2] + (TB[6][5] - 1)*sc[4];
        Corr[6] = alpha;
        alpha = mixture + (TB[7][0] - 1)*sc[0] + (TB[7][1] - 1)*sc[9] + (TB[7][2] - 1)*sc[10] + (TB[7][3] - 1)*sc[13];
        Corr[7] = alpha;
    }

    for (int i=0; i<38; i++)
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

    amrex::Real q_f[38], q_r[38];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 38; ++i) {
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
    amrex::Real c[14]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 14; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 14; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 38; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[14]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[14];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*AR */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*HE */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*HCO */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H2O */
    YOW += y[10]*imw[10]; /*CO */
    YOW += y[11]*imw[11]; /*O2 */
    YOW += y[12]*imw[12]; /*H2O2 */
    YOW += y[13]*imw[13]; /*CO2 */
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
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[14]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 14; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[14]; /*temporary storage */
    amrex::Real imw[14];

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
    c[12] = 1e6 * (*rho) * y[12]*imw[12]; 
    c[13] = 1e6 * (*rho) * y[13]*imw[13]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[14]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*39.948000; /*AR */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*4.002600; /*HE */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*29.018520; /*HCO */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*18.015340; /*H2O */
    XW += x[10]*28.010550; /*CO */
    XW += x[11]*31.998800; /*O2 */
    XW += x[12]*34.014740; /*H2O2 */
    XW += x[13]*44.009950; /*CO2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 14; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        qdot[id] *= 1.0e-6;
    }
}

/*compute the reaction Jacobian on CPU */
void aJacobian(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP)
{
    for (int i=0; i<225; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[14];
    for (int k=0; k<14; k++) {
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
    for (int k = 0; k < 14; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[14];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[14];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[14];
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
    alpha = mixture + (TB[0][0] - 1)*sc[11] + (TB[0][1] - 1)*sc[9] + (TB[0][2] - 1)*sc[10] + (TB[0][3] - 1)*sc[13] + (TB[0][4] - 1)*sc[2] + (TB[0][5] - 1)*sc[4] + (TB[0][6] - 1)*sc[0];
    /* forward */
    phi_f = sc[1]*sc[11];
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
    phi_r = sc[8];
    Kc = refCinv * exp(g_RT[1] - g_RT[8] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[11]) + (h_RT[8]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[8] += q; /* HO2 */
    wdot[11] -= q; /* O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[0][6] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[8] += dqdci;                /* dwdot[HO2]/d[H2] */
        J[11] -= dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[11];
        J[16] -= dqdci;               /* dwdot[H]/d[H] */
        J[23] += dqdci;               /* dwdot[HO2]/d[H] */
        J[26] -= dqdci;               /* dwdot[O2]/d[H] */
        /* d()/d[AR] */
        dqdci = (TB[0][4] - 1)*dcdc_fac;
        J[31] -= dqdci;               /* dwdot[H]/d[AR] */
        J[38] += dqdci;               /* dwdot[HO2]/d[AR] */
        J[41] -= dqdci;               /* dwdot[O2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[0][5] - 1)*dcdc_fac;
        J[61] -= dqdci;               /* dwdot[H]/d[HE] */
        J[68] += dqdci;               /* dwdot[HO2]/d[HE] */
        J[71] -= dqdci;               /* dwdot[O2]/d[HE] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[121] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
        J[131] -= dqdci;              /* dwdot[O2]/d[HO2] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[143] += dqdci;              /* dwdot[HO2]/d[H2O] */
        J[146] -= dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[0][2] - 1)*dcdc_fac;
        J[151] -= dqdci;              /* dwdot[H]/d[CO] */
        J[158] += dqdci;              /* dwdot[HO2]/d[CO] */
        J[161] -= dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[O2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac + k_f*sc[1];
        J[166] -= dqdci;              /* dwdot[H]/d[O2] */
        J[173] += dqdci;              /* dwdot[HO2]/d[O2] */
        J[176] -= dqdci;              /* dwdot[O2]/d[O2] */
        /* d()/d[CO2] */
        dqdci = (TB[0][3] - 1)*dcdc_fac;
        J[196] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[203] += dqdci;              /* dwdot[HO2]/d[CO2] */
        J[206] -= dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] = TB[0][6]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[11];
        dqdc[2] = TB[0][4]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[0][5]*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac - k_r;
        dqdc[9] = TB[0][1]*dcdc_fac;
        dqdc[10] = TB[0][2]*dcdc_fac;
        dqdc[11] = TB[0][0]*dcdc_fac + k_f*sc[1];
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[0][3]*dcdc_fac;
        for (int k=0; k<14; k++) {
            J[15*k+1] -= dqdc[k];
            J[15*k+8] += dqdc[k];
            J[15*k+11] -= dqdc[k];
        }
    }
    J[211] -= dqdT; /* dwdot[H]/dT */
    J[218] += dqdT; /* dwdot[HO2]/dT */
    J[221] -= dqdT; /* dwdot[O2]/dT */

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[9] + (TB[1][2] - 1)*sc[10] + (TB[1][3] - 1)*sc[13] + (TB[1][4] - 1)*sc[2] + (TB[1][5] - 1)*sc[4];
    /* forward */
    phi_f = pow(sc[6], 2.000000);
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
    phi_r = sc[12];
    Kc = refCinv * exp(g_RT[6] + g_RT[6] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[12]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[6] -= 2 * q; /* OH */
    wdot[12] += q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[6] += -2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[12] += dqdci;               /* dwdot[H2O2]/d[H2] */
        /* d()/d[AR] */
        dqdci = (TB[1][4] - 1)*dcdc_fac;
        J[36] += -2 * dqdci;          /* dwdot[OH]/d[AR] */
        J[42] += dqdci;               /* dwdot[H2O2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[1][5] - 1)*dcdc_fac;
        J[66] += -2 * dqdci;          /* dwdot[OH]/d[HE] */
        J[72] += dqdci;               /* dwdot[H2O2]/d[HE] */
        /* d()/d[OH] */
        dqdci =  + k_f*2.000000*sc[6];
        J[96] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
        J[102] += dqdci;              /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[1][1] - 1)*dcdc_fac;
        J[141] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
        J[147] += dqdci;              /* dwdot[H2O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[1][2] - 1)*dcdc_fac;
        J[156] += -2 * dqdci;         /* dwdot[OH]/d[CO] */
        J[162] += dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[H2O2] */
        dqdci =  - k_r;
        J[186] += -2 * dqdci;         /* dwdot[OH]/d[H2O2] */
        J[192] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[CO2] */
        dqdci = (TB[1][3] - 1)*dcdc_fac;
        J[201] += -2 * dqdci;         /* dwdot[OH]/d[CO2] */
        J[207] += dqdci;              /* dwdot[H2O2]/d[CO2] */
    }
    else {
        dqdc[0] = TB[1][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[1][4]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[1][5]*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac + k_f*2.000000*sc[6];
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[1][1]*dcdc_fac;
        dqdc[10] = TB[1][2]*dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac - k_r;
        dqdc[13] = TB[1][3]*dcdc_fac;
        for (int k=0; k<14; k++) {
            J[15*k+6] += -2 * dqdc[k];
            J[15*k+12] += dqdc[k];
        }
    }
    J[216] += -2 * dqdT; /* dwdot[OH]/dT */
    J[222] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[9] + (TB[2][2] - 1)*sc[10] + (TB[2][3] - 1)*sc[13] + (TB[2][4] - 1)*sc[2] + (TB[2][5] - 1)*sc[4];
    /* forward */
    phi_f = sc[5]*sc[10];
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
    phi_r = sc[13];
    Kc = refCinv * exp(g_RT[5] + g_RT[10] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[10]) + (h_RT[13]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[10] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[2][0] - 1)*dcdc_fac;
        J[5] -= dqdci;                /* dwdot[O]/d[H2] */
        J[10] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[13] += dqdci;               /* dwdot[CO2]/d[H2] */
        /* d()/d[AR] */
        dqdci = (TB[2][4] - 1)*dcdc_fac;
        J[35] -= dqdci;               /* dwdot[O]/d[AR] */
        J[40] -= dqdci;               /* dwdot[CO]/d[AR] */
        J[43] += dqdci;               /* dwdot[CO2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[2][5] - 1)*dcdc_fac;
        J[65] -= dqdci;               /* dwdot[O]/d[HE] */
        J[70] -= dqdci;               /* dwdot[CO]/d[HE] */
        J[73] += dqdci;               /* dwdot[CO2]/d[HE] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[10];
        J[80] -= dqdci;               /* dwdot[O]/d[O] */
        J[85] -= dqdci;               /* dwdot[CO]/d[O] */
        J[88] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*dcdc_fac;
        J[140] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[145] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[148] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[2][2] - 1)*dcdc_fac + k_f*sc[5];
        J[155] -= dqdci;              /* dwdot[O]/d[CO] */
        J[160] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[163] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[2][3] - 1)*dcdc_fac - k_r;
        J[200] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[205] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    }
    else {
        dqdc[0] = TB[2][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[2][4]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[2][5]*dcdc_fac;
        dqdc[5] = dcdc_fac + k_f*sc[10];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[2][1]*dcdc_fac;
        dqdc[10] = TB[2][2]*dcdc_fac + k_f*sc[5];
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[2][3]*dcdc_fac - k_r;
        for (int k=0; k<14; k++) {
            J[15*k+5] -= dqdc[k];
            J[15*k+10] -= dqdc[k];
            J[15*k+13] += dqdc[k];
        }
    }
    J[215] -= dqdT; /* dwdot[O]/dT */
    J[220] -= dqdT; /* dwdot[CO]/dT */
    J[223] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 4: H + H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[9] + (TB[3][2] - 1)*sc[13] + (TB[3][3] - 1)*sc[2] + (TB[3][4] - 1)*sc[4];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[0];
    Kc = refCinv * exp(-g_RT[0] + g_RT[1] + g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1]) + (h_RT[0]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*q_nocor - k_r;
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*2.000000*sc[1];
        J[15] += dqdci;               /* dwdot[H2]/d[H] */
        J[16] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[AR] */
        dqdci = (TB[3][3] - 1)*q_nocor;
        J[30] += dqdci;               /* dwdot[H2]/d[AR] */
        J[31] += -2 * dqdci;          /* dwdot[H]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[3][4] - 1)*q_nocor;
        J[60] += dqdci;               /* dwdot[H2]/d[HE] */
        J[61] += -2 * dqdci;          /* dwdot[H]/d[HE] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*q_nocor;
        J[135] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[136] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
        /* d()/d[CO2] */
        dqdci = (TB[3][2] - 1)*q_nocor;
        J[195] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[196] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
    }
    else {
        dqdc[0] = TB[3][0]*q_nocor - k_r;
        dqdc[1] = q_nocor + k_f*2.000000*sc[1];
        dqdc[2] = TB[3][3]*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = TB[3][4]*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = TB[3][1]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = TB[3][2]*q_nocor;
        for (int k=0; k<14; k++) {
            J[15*k+0] += dqdc[k];
            J[15*k+1] += -2 * dqdc[k];
        }
    }
    J[210] += dqdT; /* dwdot[H2]/dT */
    J[211] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 5: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[9] + (TB[4][2] - 1)*sc[10] + (TB[4][3] - 1)*sc[13] + (TB[4][4] - 1)*sc[2] + (TB[4][5] - 1)*sc[4];
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[9];
    Kc = refCinv * exp(g_RT[1] + g_RT[6] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[9]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[6] -= q; /* OH */
    wdot[9] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[4][0] - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[6] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[9] += dqdci;                /* dwdot[H2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[6];
        J[16] -= dqdci;               /* dwdot[H]/d[H] */
        J[21] -= dqdci;               /* dwdot[OH]/d[H] */
        J[24] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[AR] */
        dqdci = (TB[4][4] - 1)*q_nocor;
        J[31] -= dqdci;               /* dwdot[H]/d[AR] */
        J[36] -= dqdci;               /* dwdot[OH]/d[AR] */
        J[39] += dqdci;               /* dwdot[H2O]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[4][5] - 1)*q_nocor;
        J[61] -= dqdci;               /* dwdot[H]/d[HE] */
        J[66] -= dqdci;               /* dwdot[OH]/d[HE] */
        J[69] += dqdci;               /* dwdot[H2O]/d[HE] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[1];
        J[91] -= dqdci;               /* dwdot[H]/d[OH] */
        J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*q_nocor - k_r;
        J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[4][2] - 1)*q_nocor;
        J[151] -= dqdci;              /* dwdot[H]/d[CO] */
        J[156] -= dqdci;              /* dwdot[OH]/d[CO] */
        J[159] += dqdci;              /* dwdot[H2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[4][3] - 1)*q_nocor;
        J[196] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[201] -= dqdci;              /* dwdot[OH]/d[CO2] */
        J[204] += dqdci;              /* dwdot[H2O]/d[CO2] */
    }
    else {
        dqdc[0] = TB[4][0]*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[6];
        dqdc[2] = TB[4][4]*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = TB[4][5]*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor + k_f*sc[1];
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = TB[4][1]*q_nocor - k_r;
        dqdc[10] = TB[4][2]*q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = TB[4][3]*q_nocor;
        for (int k=0; k<14; k++) {
            J[15*k+1] -= dqdc[k];
            J[15*k+6] -= dqdc[k];
            J[15*k+9] += dqdc[k];
        }
    }
    J[211] -= dqdT; /* dwdot[H]/dT */
    J[216] -= dqdT; /* dwdot[OH]/dT */
    J[219] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 6: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[9] + (TB[5][2] - 1)*sc[10] + (TB[5][3] - 1)*sc[13] + (TB[5][4] - 1)*sc[2] + (TB[5][5] - 1)*sc[4];
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[5] -= q; /* O */
    wdot[6] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[5][0] - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[5] -= dqdci;                /* dwdot[O]/d[H2] */
        J[6] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[16] -= dqdci;               /* dwdot[H]/d[H] */
        J[20] -= dqdci;               /* dwdot[O]/d[H] */
        J[21] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[AR] */
        dqdci = (TB[5][4] - 1)*q_nocor;
        J[31] -= dqdci;               /* dwdot[H]/d[AR] */
        J[35] -= dqdci;               /* dwdot[O]/d[AR] */
        J[36] += dqdci;               /* dwdot[OH]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[5][5] - 1)*q_nocor;
        J[61] -= dqdci;               /* dwdot[H]/d[HE] */
        J[65] -= dqdci;               /* dwdot[O]/d[HE] */
        J[66] += dqdci;               /* dwdot[OH]/d[HE] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[1];
        J[76] -= dqdci;               /* dwdot[H]/d[O] */
        J[80] -= dqdci;               /* dwdot[O]/d[O] */
        J[81] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[91] -= dqdci;               /* dwdot[H]/d[OH] */
        J[95] -= dqdci;               /* dwdot[O]/d[OH] */
        J[96] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*q_nocor;
        J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[140] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[141] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[5][2] - 1)*q_nocor;
        J[151] -= dqdci;              /* dwdot[H]/d[CO] */
        J[155] -= dqdci;              /* dwdot[O]/d[CO] */
        J[156] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[5][3] - 1)*q_nocor;
        J[196] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[200] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[201] += dqdci;              /* dwdot[OH]/d[CO2] */
    }
    else {
        dqdc[0] = TB[5][0]*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[5];
        dqdc[2] = TB[5][4]*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = TB[5][5]*q_nocor;
        dqdc[5] = q_nocor + k_f*sc[1];
        dqdc[6] = q_nocor - k_r;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = TB[5][1]*q_nocor;
        dqdc[10] = TB[5][2]*q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = TB[5][3]*q_nocor;
        for (int k=0; k<14; k++) {
            J[15*k+1] -= dqdc[k];
            J[15*k+5] -= dqdc[k];
            J[15*k+6] += dqdc[k];
        }
    }
    J[211] -= dqdT; /* dwdot[H]/dT */
    J[215] -= dqdT; /* dwdot[O]/dT */
    J[216] += dqdT; /* dwdot[OH]/dT */

    /*reaction 7: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[9] + (TB[6][2] - 1)*sc[10] + (TB[6][3] - 1)*sc[13] + (TB[6][4] - 1)*sc[2] + (TB[6][5] - 1)*sc[4];
    /* forward */
    phi_f = pow(sc[5], 2.000000);
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[11];
    Kc = refCinv * exp(g_RT[5] + g_RT[5] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[5]) + (h_RT[11]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= 2 * q; /* O */
    wdot[11] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[6][0] - 1)*q_nocor;
        J[5] += -2 * dqdci;           /* dwdot[O]/d[H2] */
        J[11] += dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[AR] */
        dqdci = (TB[6][4] - 1)*q_nocor;
        J[35] += -2 * dqdci;          /* dwdot[O]/d[AR] */
        J[41] += dqdci;               /* dwdot[O2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[6][5] - 1)*q_nocor;
        J[65] += -2 * dqdci;          /* dwdot[O]/d[HE] */
        J[71] += dqdci;               /* dwdot[O2]/d[HE] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[5];
        J[80] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[86] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (TB[6][1] - 1)*q_nocor;
        J[140] += -2 * dqdci;         /* dwdot[O]/d[H2O] */
        J[146] += dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[6][2] - 1)*q_nocor;
        J[155] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[161] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[170] += -2 * dqdci;         /* dwdot[O]/d[O2] */
        J[176] += dqdci;              /* dwdot[O2]/d[O2] */
        /* d()/d[CO2] */
        dqdci = (TB[6][3] - 1)*q_nocor;
        J[200] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[206] += dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] = TB[6][0]*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = TB[6][4]*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = TB[6][5]*q_nocor;
        dqdc[5] = q_nocor + k_f*2.000000*sc[5];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = TB[6][1]*q_nocor;
        dqdc[10] = TB[6][2]*q_nocor;
        dqdc[11] = q_nocor - k_r;
        dqdc[12] = q_nocor;
        dqdc[13] = TB[6][3]*q_nocor;
        for (int k=0; k<14; k++) {
            J[15*k+5] += -2 * dqdc[k];
            J[15*k+11] += dqdc[k];
        }
    }
    J[215] += -2 * dqdT; /* dwdot[O]/dT */
    J[221] += dqdT; /* dwdot[O2]/dT */

    /*reaction 8: HCO + M <=> CO + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[7][0] - 1)*sc[0] + (TB[7][1] - 1)*sc[9] + (TB[7][2] - 1)*sc[10] + (TB[7][3] - 1)*sc[13];
    /* forward */
    phi_f = sc[7];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[10];
    Kc = refC * exp(-g_RT[1] + g_RT[7] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (h_RT[1] + h_RT[10]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* HCO */
    wdot[10] += q; /* CO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[7][0] - 1)*q_nocor;
        J[1] += dqdci;                /* dwdot[H]/d[H2] */
        J[7] -= dqdci;                /* dwdot[HCO]/d[H2] */
        J[10] += dqdci;               /* dwdot[CO]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[10];
        J[16] += dqdci;               /* dwdot[H]/d[H] */
        J[22] -= dqdci;               /* dwdot[HCO]/d[H] */
        J[25] += dqdci;               /* dwdot[CO]/d[H] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[106] += dqdci;              /* dwdot[H]/d[HCO] */
        J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        J[115] += dqdci;              /* dwdot[CO]/d[HCO] */
        /* d()/d[H2O] */
        dqdci = (TB[7][1] - 1)*q_nocor;
        J[136] += dqdci;              /* dwdot[H]/d[H2O] */
        J[142] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        J[145] += dqdci;              /* dwdot[CO]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[7][2] - 1)*q_nocor - k_r*sc[1];
        J[151] += dqdci;              /* dwdot[H]/d[CO] */
        J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
        J[160] += dqdci;              /* dwdot[CO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[7][3] - 1)*q_nocor;
        J[196] += dqdci;              /* dwdot[H]/d[CO2] */
        J[202] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        J[205] += dqdci;              /* dwdot[CO]/d[CO2] */
    }
    else {
        dqdc[0] = TB[7][0]*q_nocor;
        dqdc[1] = q_nocor - k_r*sc[10];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor + k_f;
        dqdc[8] = q_nocor;
        dqdc[9] = TB[7][1]*q_nocor;
        dqdc[10] = TB[7][2]*q_nocor - k_r*sc[1];
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = TB[7][3]*q_nocor;
        for (int k=0; k<14; k++) {
            J[15*k+1] += dqdc[k];
            J[15*k+7] -= dqdc[k];
            J[15*k+10] += dqdc[k];
        }
    }
    J[211] += dqdT; /* dwdot[H]/dT */
    J[217] -= dqdT; /* dwdot[HCO]/dT */
    J[220] += dqdT; /* dwdot[CO]/dT */

    /*reaction 9: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[1] - g_RT[5] - g_RT[6] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[11]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[5] += q; /* O */
    wdot[6] += q; /* OH */
    wdot[11] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[11];
    J[16] -= dqdci;               /* dwdot[H]/d[H] */
    J[20] += dqdci;               /* dwdot[O]/d[H] */
    J[21] += dqdci;               /* dwdot[OH]/d[H] */
    J[26] -= dqdci;               /* dwdot[O2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[6];
    J[76] -= dqdci;               /* dwdot[H]/d[O] */
    J[80] += dqdci;               /* dwdot[O]/d[O] */
    J[81] += dqdci;               /* dwdot[OH]/d[O] */
    J[86] -= dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[91] -= dqdci;               /* dwdot[H]/d[OH] */
    J[95] += dqdci;               /* dwdot[O]/d[OH] */
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    J[101] -= dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[166] -= dqdci;              /* dwdot[H]/d[O2] */
    J[170] += dqdci;              /* dwdot[O]/d[O2] */
    J[171] += dqdci;              /* dwdot[OH]/d[O2] */
    J[176] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[211] -= dqdT;               /* dwdot[H]/dT */
    J[215] += dqdT;               /* dwdot[O]/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */
    J[221] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 10: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[6];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[1] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[5] -= q; /* O */
    wdot[6] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[5] -= dqdci;                /* dwdot[O]/d[H2] */
    J[6] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[6];
    J[15] -= dqdci;               /* dwdot[H2]/d[H] */
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[20] -= dqdci;               /* dwdot[O]/d[H] */
    J[21] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[75] -= dqdci;               /* dwdot[H2]/d[O] */
    J[76] += dqdci;               /* dwdot[H]/d[O] */
    J[80] -= dqdci;               /* dwdot[O]/d[O] */
    J[81] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[90] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[91] += dqdci;               /* dwdot[H]/d[OH] */
    J[95] -= dqdci;               /* dwdot[O]/d[OH] */
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[210] -= dqdT;               /* dwdot[H2]/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[215] -= dqdT;               /* dwdot[O]/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 11: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[6] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[6] -= q; /* OH */
    wdot[9] += q; /* H2O */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[6] -= dqdci;                /* dwdot[OH]/d[H2] */
    J[9] += dqdci;                /* dwdot[H2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[15] -= dqdci;               /* dwdot[H2]/d[H] */
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[21] -= dqdci;               /* dwdot[OH]/d[H] */
    J[24] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[90] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[91] += dqdci;               /* dwdot[H]/d[OH] */
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[135] -= dqdci;              /* dwdot[H2]/d[H2O] */
    J[136] += dqdci;              /* dwdot[H]/d[H2O] */
    J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[210] -= dqdT;               /* dwdot[H2]/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 12: OH + OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O */
    wdot[6] -= 2 * q; /* OH */
    wdot[9] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[9];
    J[80] += dqdci;               /* dwdot[O]/d[O] */
    J[81] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[84] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[6];
    J[95] += dqdci;               /* dwdot[O]/d[OH] */
    J[96] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[140] += dqdci;              /* dwdot[O]/d[H2O] */
    J[141] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[215] += dqdT;               /* dwdot[O]/dT */
    J[216] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*pow(sc[1], 2.000000);
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000);
    Kc = refCinv * exp(g_RT[0] - g_RT[0] - g_RT[0] + g_RT[1] + g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + 2.000000*h_RT[1]) + (2.000000*h_RT[0]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*2.000000*sc[0];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[0]*2.000000*sc[1];
    J[15] += dqdci;               /* dwdot[H2]/d[H] */
    J[16] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/dT */
    J[210] += dqdT;               /* dwdot[H2]/dT */
    J[211] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[9];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[9];
    Kc = refCinv * exp(-g_RT[0] + g_RT[1] + g_RT[1] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[9]) + (h_RT[0] + h_RT[9]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[9];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[9];
    J[15] += dqdci;               /* dwdot[H2]/d[H] */
    J[16] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[135] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[136] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
    /* d()/dT */
    J[210] += dqdT;               /* dwdot[H2]/dT */
    J[211] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[13];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[13];
    Kc = refCinv * exp(-g_RT[0] + g_RT[1] + g_RT[1] + g_RT[13] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[13]) + (h_RT[0] + h_RT[13]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[13];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[13];
    J[15] += dqdci;               /* dwdot[H2]/d[H] */
    J[16] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[CO2] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[195] += dqdci;              /* dwdot[H2]/d[CO2] */
    J[196] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
    /* d()/dT */
    J[210] += dqdT;               /* dwdot[H2]/dT */
    J[211] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 16: H2 + O2 <=> HO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[11];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[8];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[8] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[11]) + (h_RT[1] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[8] += q; /* HO2 */
    wdot[11] -= q; /* O2 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[11];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[8] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[11] -= dqdci;               /* dwdot[O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[8];
    J[15] -= dqdci;               /* dwdot[H2]/d[H] */
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[23] += dqdci;               /* dwdot[HO2]/d[H] */
    J[26] -= dqdci;               /* dwdot[O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[1];
    J[120] -= dqdci;              /* dwdot[H2]/d[HO2] */
    J[121] += dqdci;              /* dwdot[H]/d[HO2] */
    J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[131] -= dqdci;              /* dwdot[O2]/d[HO2] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[165] -= dqdci;              /* dwdot[H2]/d[O2] */
    J[166] += dqdci;              /* dwdot[H]/d[O2] */
    J[173] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[176] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[210] -= dqdT;               /* dwdot[H2]/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[218] += dqdT;               /* dwdot[HO2]/dT */
    J[221] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 17: HO2 + H <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[8];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(g_RT[1] - g_RT[5] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[8]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[5] += q; /* O */
    wdot[8] -= q; /* HO2 */
    wdot[9] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[16] -= dqdci;               /* dwdot[H]/d[H] */
    J[20] += dqdci;               /* dwdot[O]/d[H] */
    J[23] -= dqdci;               /* dwdot[HO2]/d[H] */
    J[24] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[9];
    J[76] -= dqdci;               /* dwdot[H]/d[O] */
    J[80] += dqdci;               /* dwdot[O]/d[O] */
    J[83] -= dqdci;               /* dwdot[HO2]/d[O] */
    J[84] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[121] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[125] += dqdci;              /* dwdot[O]/d[HO2] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[129] += dqdci;              /* dwdot[H2O]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[140] += dqdci;              /* dwdot[O]/d[H2O] */
    J[143] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[211] -= dqdT;               /* dwdot[H]/dT */
    J[215] += dqdT;               /* dwdot[O]/dT */
    J[218] -= dqdT;               /* dwdot[HO2]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[8];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = pow(sc[6], 2.000000);
    Kc = exp(g_RT[1] - g_RT[6] - g_RT[6] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[8]) + (2.000000*h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[6] += 2 * q; /* OH */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[16] -= dqdci;               /* dwdot[H]/d[H] */
    J[21] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[23] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[6];
    J[91] -= dqdci;               /* dwdot[H]/d[OH] */
    J[96] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[121] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[126] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[211] -= dqdT;               /* dwdot[H]/dT */
    J[216] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[218] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 19: HO2 + O <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[8];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[11];
    Kc = exp(g_RT[5] - g_RT[6] + g_RT[8] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[8]) + (h_RT[6] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[6] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[11] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[8];
    J[80] -= dqdci;               /* dwdot[O]/d[O] */
    J[81] += dqdci;               /* dwdot[OH]/d[O] */
    J[83] -= dqdci;               /* dwdot[HO2]/d[O] */
    J[86] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[11];
    J[95] -= dqdci;               /* dwdot[O]/d[OH] */
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[101] += dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[125] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[126] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[131] += dqdci;              /* dwdot[O2]/d[HO2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[6];
    J[170] -= dqdci;              /* dwdot[O]/d[O2] */
    J[171] += dqdci;              /* dwdot[OH]/d[O2] */
    J[173] -= dqdci;              /* dwdot[HO2]/d[O2] */
    J[176] += dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[215] -= dqdT;               /* dwdot[O]/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */
    J[218] -= dqdT;               /* dwdot[HO2]/dT */
    J[221] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[8];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[11];
    Kc = exp(g_RT[6] + g_RT[8] - g_RT[9] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[9] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[9] += q; /* H2O */
    wdot[11] += q; /* O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[101] += dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[6];
    J[126] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[129] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[131] += dqdci;              /* dwdot[O2]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[143] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[146] += dqdci;              /* dwdot[O2]/d[H2O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[9];
    J[171] -= dqdci;              /* dwdot[OH]/d[O2] */
    J[173] -= dqdci;              /* dwdot[HO2]/d[O2] */
    J[174] += dqdci;              /* dwdot[H2O]/d[O2] */
    J[176] += dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[218] -= dqdT;               /* dwdot[HO2]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */
    J[221] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[8];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[11];
    Kc = exp(g_RT[6] + g_RT[8] - g_RT[9] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[9] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[9] += q; /* H2O */
    wdot[11] += q; /* O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[101] += dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[6];
    J[126] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[129] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[131] += dqdci;              /* dwdot[O2]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[143] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[146] += dqdci;              /* dwdot[O2]/d[H2O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[9];
    J[171] -= dqdci;              /* dwdot[OH]/d[O2] */
    J[173] -= dqdci;              /* dwdot[HO2]/d[O2] */
    J[174] += dqdci;              /* dwdot[H2O]/d[O2] */
    J[176] += dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[218] -= dqdT;               /* dwdot[HO2]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */
    J[221] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[8], 2.000000);
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[12];
    Kc = exp(g_RT[8] + g_RT[8] - g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[8]) + (h_RT[11] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= 2 * q; /* HO2 */
    wdot[11] += q; /* O2 */
    wdot[12] += q; /* H2O2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[8];
    J[128] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    J[131] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[132] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[12];
    J[173] += -2 * dqdci;         /* dwdot[HO2]/d[O2] */
    J[176] += dqdci;              /* dwdot[O2]/d[O2] */
    J[177] += dqdci;              /* dwdot[H2O2]/d[O2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[11];
    J[188] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[191] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[192] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[218] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[221] += dqdT;               /* dwdot[O2]/dT */
    J[222] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[8], 2.000000);
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[12];
    Kc = exp(g_RT[8] + g_RT[8] - g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[8]) + (h_RT[11] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= 2 * q; /* HO2 */
    wdot[11] += q; /* O2 */
    wdot[12] += q; /* H2O2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[8];
    J[128] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    J[131] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[132] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[12];
    J[173] += -2 * dqdci;         /* dwdot[HO2]/d[O2] */
    J[176] += dqdci;              /* dwdot[O2]/d[O2] */
    J[177] += dqdci;              /* dwdot[H2O2]/d[O2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[11];
    J[188] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[191] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[192] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[218] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[221] += dqdT;               /* dwdot[O2]/dT */
    J[222] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[12];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[8];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[8] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[12]) + (h_RT[0] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[8] += q; /* HO2 */
    wdot[12] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[8];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[8] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[12] -= dqdci;               /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[12];
    J[15] += dqdci;               /* dwdot[H2]/d[H] */
    J[16] -= dqdci;               /* dwdot[H]/d[H] */
    J[23] += dqdci;               /* dwdot[HO2]/d[H] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[120] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[121] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[132] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[180] += dqdci;              /* dwdot[H2]/d[H2O2] */
    J[181] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[192] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[210] += dqdT;               /* dwdot[H2]/dT */
    J[211] -= dqdT;               /* dwdot[H]/dT */
    J[218] += dqdT;               /* dwdot[HO2]/dT */
    J[222] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 25: H2O2 + H <=> OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[12];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[9];
    Kc = exp(g_RT[1] - g_RT[6] - g_RT[9] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[12]) + (h_RT[6] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[6] += q; /* OH */
    wdot[9] += q; /* H2O */
    wdot[12] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[12];
    J[16] -= dqdci;               /* dwdot[H]/d[H] */
    J[21] += dqdci;               /* dwdot[OH]/d[H] */
    J[24] += dqdci;               /* dwdot[H2O]/d[H] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[91] -= dqdci;               /* dwdot[H]/d[OH] */
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[102] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[141] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[147] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[181] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[186] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[189] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[192] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[211] -= dqdT;               /* dwdot[H]/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */
    J[222] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[12];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[8];
    Kc = exp(g_RT[5] - g_RT[6] - g_RT[8] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[12]) + (h_RT[6] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[6] += q; /* OH */
    wdot[8] += q; /* HO2 */
    wdot[12] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[12];
    J[80] -= dqdci;               /* dwdot[O]/d[O] */
    J[81] += dqdci;               /* dwdot[OH]/d[O] */
    J[83] += dqdci;               /* dwdot[HO2]/d[O] */
    J[87] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[8];
    J[95] -= dqdci;               /* dwdot[O]/d[OH] */
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    J[98] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[102] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[6];
    J[125] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[126] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[132] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[185] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[186] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[192] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[215] -= dqdT;               /* dwdot[O]/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */
    J[218] += dqdT;               /* dwdot[HO2]/dT */
    J[222] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[12];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[9];
    Kc = exp(g_RT[6] - g_RT[8] - g_RT[9] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[8] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[8] += q; /* HO2 */
    wdot[9] += q; /* H2O */
    wdot[12] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[98] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[102] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[9];
    J[126] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[129] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[132] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[143] += dqdci;              /* dwdot[HO2]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[147] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[6];
    J[186] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[189] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[192] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[218] += dqdT;               /* dwdot[HO2]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */
    J[222] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[12];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[9];
    Kc = exp(g_RT[6] - g_RT[8] - g_RT[9] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[8] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[8] += q; /* HO2 */
    wdot[9] += q; /* H2O */
    wdot[12] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[98] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[102] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[9];
    J[126] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[129] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[132] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[143] += dqdci;              /* dwdot[HO2]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[147] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[6];
    J[186] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[189] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[192] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[218] += dqdT;               /* dwdot[HO2]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */
    J[222] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[10];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[10] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[10]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[6] -= q; /* OH */
    wdot[10] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[21] -= dqdci;               /* dwdot[OH]/d[H] */
    J[25] -= dqdci;               /* dwdot[CO]/d[H] */
    J[28] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[91] += dqdci;               /* dwdot[H]/d[OH] */
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[100] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[103] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[151] += dqdci;              /* dwdot[H]/d[CO] */
    J[156] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[160] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[163] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[196] += dqdci;              /* dwdot[H]/d[CO2] */
    J[201] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[205] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[220] -= dqdT;               /* dwdot[CO]/dT */
    J[223] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 30: CO + OH <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[10];
    k_f = prefactor_units[29] * fwd_A[29]
                * exp(fwd_beta[29] * tc[0] - activation_units[29] * fwd_Ea[29] * invT);
    dlnkfdT = fwd_beta[29] * invT + activation_units[29] * fwd_Ea[29] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[10] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[10]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[6] -= q; /* OH */
    wdot[10] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[21] -= dqdci;               /* dwdot[OH]/d[H] */
    J[25] -= dqdci;               /* dwdot[CO]/d[H] */
    J[28] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[91] += dqdci;               /* dwdot[H]/d[OH] */
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[100] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[103] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[151] += dqdci;              /* dwdot[H]/d[CO] */
    J[156] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[160] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[163] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[196] += dqdci;              /* dwdot[H]/d[CO2] */
    J[201] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[205] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[220] -= dqdT;               /* dwdot[CO]/dT */
    J[223] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 31: CO + O2 <=> CO2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[11];
    k_f = prefactor_units[30] * fwd_A[30]
                * exp(fwd_beta[30] * tc[0] - activation_units[30] * fwd_Ea[30] * invT);
    dlnkfdT = fwd_beta[30] * invT + activation_units[30] * fwd_Ea[30] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[13];
    Kc = exp(-g_RT[5] + g_RT[10] + g_RT[11] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[11]) + (h_RT[5] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O */
    wdot[10] -= q; /* CO */
    wdot[11] -= q; /* O2 */
    wdot[13] += q; /* CO2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[13];
    J[80] += dqdci;               /* dwdot[O]/d[O] */
    J[85] -= dqdci;               /* dwdot[CO]/d[O] */
    J[86] -= dqdci;               /* dwdot[O2]/d[O] */
    J[88] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[11];
    J[155] += dqdci;              /* dwdot[O]/d[CO] */
    J[160] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[161] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[163] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[170] += dqdci;              /* dwdot[O]/d[O2] */
    J[175] -= dqdci;              /* dwdot[CO]/d[O2] */
    J[176] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[178] += dqdci;              /* dwdot[CO2]/d[O2] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[5];
    J[200] += dqdci;              /* dwdot[O]/d[CO2] */
    J[205] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[206] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[215] += dqdT;               /* dwdot[O]/dT */
    J[220] -= dqdT;               /* dwdot[CO]/dT */
    J[221] -= dqdT;               /* dwdot[O2]/dT */
    J[223] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = prefactor_units[31] * fwd_A[31]
                * exp(fwd_beta[31] * tc[0] - activation_units[31] * fwd_Ea[31] * invT);
    dlnkfdT = fwd_beta[31] * invT + activation_units[31] * fwd_Ea[31] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[13];
    Kc = exp(-g_RT[6] + g_RT[8] + g_RT[10] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[10]) + (h_RT[6] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[10] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[100] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[103] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[10];
    J[126] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[130] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[133] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[8];
    J[156] += dqdci;              /* dwdot[OH]/d[CO] */
    J[158] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[160] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[163] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[6];
    J[201] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[203] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[205] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */
    J[218] -= dqdT;               /* dwdot[HO2]/dT */
    J[220] -= dqdT;               /* dwdot[CO]/dT */
    J[223] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 33: HCO + H <=> CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = prefactor_units[32] * fwd_A[32]
                * exp(fwd_beta[32] * tc[0] - activation_units[32] * fwd_Ea[32] * invT);
    dlnkfdT = fwd_beta[32] * invT + activation_units[32] * fwd_Ea[32] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[10];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[7] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[0] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[7] -= q; /* HCO */
    wdot[10] += q; /* CO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[10];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[7] -= dqdci;                /* dwdot[HCO]/d[H2] */
    J[10] += dqdci;               /* dwdot[CO]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[15] += dqdci;               /* dwdot[H2]/d[H] */
    J[16] -= dqdci;               /* dwdot[H]/d[H] */
    J[22] -= dqdci;               /* dwdot[HCO]/d[H] */
    J[25] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[105] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[106] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[115] += dqdci;              /* dwdot[CO]/d[HCO] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[0];
    J[150] += dqdci;              /* dwdot[H2]/d[CO] */
    J[151] -= dqdci;              /* dwdot[H]/d[CO] */
    J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[160] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[210] += dqdT;               /* dwdot[H2]/dT */
    J[211] -= dqdT;               /* dwdot[H]/dT */
    J[217] -= dqdT;               /* dwdot[HCO]/dT */
    J[220] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 34: HCO + O <=> CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[33] * fwd_A[33]
                * exp(fwd_beta[33] * tc[0] - activation_units[33] * fwd_Ea[33] * invT);
    dlnkfdT = fwd_beta[33] * invT + activation_units[33] * fwd_Ea[33] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[10];
    Kc = exp(g_RT[5] - g_RT[6] + g_RT[7] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[6] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[6] += q; /* OH */
    wdot[7] -= q; /* HCO */
    wdot[10] += q; /* CO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[80] -= dqdci;               /* dwdot[O]/d[O] */
    J[81] += dqdci;               /* dwdot[OH]/d[O] */
    J[82] -= dqdci;               /* dwdot[HCO]/d[O] */
    J[85] += dqdci;               /* dwdot[CO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[10];
    J[95] -= dqdci;               /* dwdot[O]/d[OH] */
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] -= dqdci;               /* dwdot[HCO]/d[OH] */
    J[100] += dqdci;              /* dwdot[CO]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[110] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[111] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[115] += dqdci;              /* dwdot[CO]/d[HCO] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[155] -= dqdci;              /* dwdot[O]/d[CO] */
    J[156] += dqdci;              /* dwdot[OH]/d[CO] */
    J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[160] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[215] -= dqdT;               /* dwdot[O]/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */
    J[217] -= dqdT;               /* dwdot[HCO]/dT */
    J[220] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 35: HCO + O <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[34] * fwd_A[34]
                * exp(fwd_beta[34] * tc[0] - activation_units[34] * fwd_Ea[34] * invT);
    dlnkfdT = fwd_beta[34] * invT + activation_units[34] * fwd_Ea[34] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[5] + g_RT[7] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[5] -= q; /* O */
    wdot[7] -= q; /* HCO */
    wdot[13] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[20] -= dqdci;               /* dwdot[O]/d[H] */
    J[22] -= dqdci;               /* dwdot[HCO]/d[H] */
    J[28] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[76] += dqdci;               /* dwdot[H]/d[O] */
    J[80] -= dqdci;               /* dwdot[O]/d[O] */
    J[82] -= dqdci;               /* dwdot[HCO]/d[O] */
    J[88] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[106] += dqdci;              /* dwdot[H]/d[HCO] */
    J[110] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[118] += dqdci;              /* dwdot[CO2]/d[HCO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[196] += dqdci;              /* dwdot[H]/d[CO2] */
    J[200] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[202] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[215] -= dqdT;               /* dwdot[O]/dT */
    J[217] -= dqdT;               /* dwdot[HCO]/dT */
    J[223] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 36: HCO + OH <=> CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[7];
    k_f = prefactor_units[35] * fwd_A[35]
                * exp(fwd_beta[35] * tc[0] - activation_units[35] * fwd_Ea[35] * invT);
    dlnkfdT = fwd_beta[35] * invT + activation_units[35] * fwd_Ea[35] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[10];
    Kc = exp(g_RT[6] + g_RT[7] - g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[7]) + (h_RT[9] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[7] -= q; /* HCO */
    wdot[9] += q; /* H2O */
    wdot[10] += q; /* CO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[97] -= dqdci;               /* dwdot[HCO]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[100] += dqdci;              /* dwdot[CO]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[6];
    J[111] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[114] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[115] += dqdci;              /* dwdot[CO]/d[HCO] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[10];
    J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[142] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[145] += dqdci;              /* dwdot[CO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[9];
    J[156] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[159] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[160] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[217] -= dqdT;               /* dwdot[HCO]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */
    J[220] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[9];
    k_f = prefactor_units[36] * fwd_A[36]
                * exp(fwd_beta[36] * tc[0] - activation_units[36] * fwd_Ea[36] * invT);
    dlnkfdT = fwd_beta[36] * invT + activation_units[36] * fwd_Ea[36] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9]*sc[10];
    Kc = refC * exp(-g_RT[1] + g_RT[7] + g_RT[9] - g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[9]) + (h_RT[1] + h_RT[9] + h_RT[10]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* HCO */
    wdot[10] += q; /* CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9]*sc[10];
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[22] -= dqdci;               /* dwdot[HCO]/d[H] */
    J[25] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[9];
    J[106] += dqdci;              /* dwdot[H]/d[HCO] */
    J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[115] += dqdci;              /* dwdot[CO]/d[HCO] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[7] - k_r*sc[1]*sc[10];
    J[136] += dqdci;              /* dwdot[H]/d[H2O] */
    J[142] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    J[145] += dqdci;              /* dwdot[CO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[9];
    J[151] += dqdci;              /* dwdot[H]/d[CO] */
    J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[160] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[217] -= dqdT;               /* dwdot[HCO]/dT */
    J[220] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[11];
    k_f = prefactor_units[37] * fwd_A[37]
                * exp(fwd_beta[37] * tc[0] - activation_units[37] * fwd_Ea[37] * invT);
    dlnkfdT = fwd_beta[37] * invT + activation_units[37] * fwd_Ea[37] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[10];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[10] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[11]) + (h_RT[8] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* HCO */
    wdot[8] += q; /* HO2 */
    wdot[10] += q; /* CO */
    wdot[11] -= q; /* O2 */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[11];
    J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[113] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[115] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[116] -= dqdci;              /* dwdot[O2]/d[HCO] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[10];
    J[127] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[130] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[131] -= dqdci;              /* dwdot[O2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[8];
    J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[158] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[160] += dqdci;              /* dwdot[CO]/d[CO] */
    J[161] -= dqdci;              /* dwdot[O2]/d[CO] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[7];
    J[172] -= dqdci;              /* dwdot[HCO]/d[O2] */
    J[173] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[175] += dqdci;              /* dwdot[CO]/d[O2] */
    J[176] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[217] -= dqdT;               /* dwdot[HCO]/dT */
    J[218] += dqdT;               /* dwdot[HO2]/dT */
    J[220] += dqdT;               /* dwdot[CO]/dT */
    J[221] -= dqdT;               /* dwdot[O2]/dT */

    amrex::Real c_R[14], dcRdT[14], e_RT[14];
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
    for (int k = 0; k < 14; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[210+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 14; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 14; ++m) {
            dehmixdc += eh_RT[m]*J[k*15+m];
        }
        J[k*15+14] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[224] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

#endif
/* Transport function declarations  */


void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 59;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 4270;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 14;}


void egtransetNLITE(int* NLITE ) {
    *NLITE = 3;}


/*Patm in ergs/cm3 */
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 2.01594000E+00;
    WT[1] = 1.00797000E+00;
    WT[2] = 3.99480000E+01;
    WT[3] = 2.80134000E+01;
    WT[4] = 4.00260000E+00;
    WT[5] = 1.59994000E+01;
    WT[6] = 1.70073700E+01;
    WT[7] = 2.90185200E+01;
    WT[8] = 3.30067700E+01;
    WT[9] = 1.80153400E+01;
    WT[10] = 2.80105500E+01;
    WT[11] = 3.19988000E+01;
    WT[12] = 3.40147400E+01;
    WT[13] = 4.40099500E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 3.80000000E+01;
    EPS[1] = 1.45000000E+02;
    EPS[2] = 1.36500000E+02;
    EPS[3] = 9.75300000E+01;
    EPS[4] = 1.02000000E+01;
    EPS[5] = 8.00000000E+01;
    EPS[6] = 8.00000000E+01;
    EPS[7] = 4.98000000E+02;
    EPS[8] = 1.07400000E+02;
    EPS[9] = 5.72400000E+02;
    EPS[10] = 9.81000000E+01;
    EPS[11] = 1.07400000E+02;
    EPS[12] = 1.07400000E+02;
    EPS[13] = 2.44000000E+02;
}


/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 2.92000000E+00;
    SIG[1] = 2.05000000E+00;
    SIG[2] = 3.33000000E+00;
    SIG[3] = 3.62100000E+00;
    SIG[4] = 2.57600000E+00;
    SIG[5] = 2.75000000E+00;
    SIG[6] = 2.75000000E+00;
    SIG[7] = 3.59000000E+00;
    SIG[8] = 3.45800000E+00;
    SIG[9] = 2.60500000E+00;
    SIG[10] = 3.65000000E+00;
    SIG[11] = 3.45800000E+00;
    SIG[12] = 3.45800000E+00;
    SIG[13] = 3.76300000E+00;
}


/*the dipole moment in Debye */
void egtransetDIP(amrex::Real* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[2] = 0.00000000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 0.00000000E+00;
    DIP[6] = 0.00000000E+00;
    DIP[7] = 0.00000000E+00;
    DIP[8] = 0.00000000E+00;
    DIP[9] = 1.84400000E+00;
    DIP[10] = 0.00000000E+00;
    DIP[11] = 0.00000000E+00;
    DIP[12] = 0.00000000E+00;
    DIP[13] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 7.90000000E-01;
    POL[1] = 0.00000000E+00;
    POL[2] = 0.00000000E+00;
    POL[3] = 1.76000000E+00;
    POL[4] = 0.00000000E+00;
    POL[5] = 0.00000000E+00;
    POL[6] = 0.00000000E+00;
    POL[7] = 0.00000000E+00;
    POL[8] = 0.00000000E+00;
    POL[9] = 0.00000000E+00;
    POL[10] = 1.95000000E+00;
    POL[11] = 1.60000000E+00;
    POL[12] = 0.00000000E+00;
    POL[13] = 2.65000000E+00;
}


/*the rotational relaxation collision number at 298 K */
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 2.80000000E+02;
    ZROT[1] = 0.00000000E+00;
    ZROT[2] = 0.00000000E+00;
    ZROT[3] = 4.00000000E+00;
    ZROT[4] = 0.00000000E+00;
    ZROT[5] = 0.00000000E+00;
    ZROT[6] = 0.00000000E+00;
    ZROT[7] = 0.00000000E+00;
    ZROT[8] = 1.00000000E+00;
    ZROT[9] = 4.00000000E+00;
    ZROT[10] = 1.80000000E+00;
    ZROT[11] = 3.80000000E+00;
    ZROT[12] = 3.80000000E+00;
    ZROT[13] = 2.10000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 1;
    NLIN[1] = 0;
    NLIN[2] = 0;
    NLIN[3] = 1;
    NLIN[4] = 0;
    NLIN[5] = 0;
    NLIN[6] = 1;
    NLIN[7] = 2;
    NLIN[8] = 2;
    NLIN[9] = 2;
    NLIN[10] = 1;
    NLIN[11] = 1;
    NLIN[12] = 2;
    NLIN[13] = 1;
}


/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -1.37549435E+01;
    COFETA[1] = 9.65530587E-01;
    COFETA[2] = -4.45720114E-02;
    COFETA[3] = 2.05871810E-03;
    COFETA[4] = -1.98744496E+01;
    COFETA[5] = 3.41660514E+00;
    COFETA[6] = -3.63206306E-01;
    COFETA[7] = 1.58671021E-02;
    COFETA[8] = -1.86067598E+01;
    COFETA[9] = 3.27402596E+00;
    COFETA[10] = -3.45827972E-01;
    COFETA[11] = 1.51622680E-02;
    COFETA[12] = -1.62526779E+01;
    COFETA[13] = 2.24839597E+00;
    COFETA[14] = -2.13428438E-01;
    COFETA[15] = 9.46192413E-03;
    COFETA[16] = -1.11555213E+01;
    COFETA[17] = 2.18772782E-01;
    COFETA[18] = 5.60263799E-02;
    COFETA[19] = -2.36018246E-03;
    COFETA[20] = -1.48001581E+01;
    COFETA[21] = 1.79491990E+00;
    COFETA[22] = -1.54008440E-01;
    COFETA[23] = 6.86719439E-03;
    COFETA[24] = -1.47696103E+01;
    COFETA[25] = 1.79491990E+00;
    COFETA[26] = -1.54008440E-01;
    COFETA[27] = 6.86719439E-03;
    COFETA[28] = -2.11306792E+01;
    COFETA[29] = 3.26961843E+00;
    COFETA[30] = -2.51355092E-01;
    COFETA[31] = 7.35605058E-03;
    COFETA[32] = -1.67963797E+01;
    COFETA[33] = 2.52362554E+00;
    COFETA[34] = -2.49309128E-01;
    COFETA[35] = 1.10211025E-02;
    COFETA[36] = -1.17770937E+01;
    COFETA[37] = -8.26742721E-01;
    COFETA[38] = 3.39009079E-01;
    COFETA[39] = -2.00674327E-02;
    COFETA[40] = -1.63031240E+01;
    COFETA[41] = 2.26143219E+00;
    COFETA[42] = -2.15114671E-01;
    COFETA[43] = 9.53461976E-03;
    COFETA[44] = -1.68118868E+01;
    COFETA[45] = 2.52362554E+00;
    COFETA[46] = -2.49309128E-01;
    COFETA[47] = 1.10211025E-02;
    COFETA[48] = -1.67813391E+01;
    COFETA[49] = 2.52362554E+00;
    COFETA[50] = -2.49309128E-01;
    COFETA[51] = 1.10211025E-02;
    COFETA[52] = -2.36749526E+01;
    COFETA[53] = 4.99775518E+00;
    COFETA[54] = -5.52687718E-01;
    COFETA[55] = 2.34353338E-02;
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = 1.15899058E+01;
    COFLAM[1] = -1.52427727E+00;
    COFLAM[2] = 2.72840752E-01;
    COFLAM[3] = -1.03392618E-02;
    COFLAM[4] = -3.24539191E-01;
    COFLAM[5] = 3.41660514E+00;
    COFLAM[6] = -3.63206306E-01;
    COFLAM[7] = 1.58671021E-02;
    COFLAM[8] = -2.73648952E+00;
    COFLAM[9] = 3.27402596E+00;
    COFLAM[10] = -3.45827972E-01;
    COFLAM[11] = 1.51622680E-02;
    COFLAM[12] = 1.15507063E+01;
    COFLAM[13] = -2.91452379E+00;
    COFLAM[14] = 5.55043580E-01;
    COFLAM[15] = -2.75172461E-02;
    COFLAM[16] = 7.01538340E+00;
    COFLAM[17] = 2.18772782E-01;
    COFLAM[18] = 5.60263799E-02;
    COFLAM[19] = -2.36018246E-03;
    COFLAM[20] = 1.98513952E+00;
    COFLAM[21] = 1.79491990E+00;
    COFLAM[22] = -1.54008440E-01;
    COFLAM[23] = 6.86719439E-03;
    COFLAM[24] = 1.60618734E+01;
    COFLAM[25] = -4.10626869E+00;
    COFLAM[26] = 6.63571339E-01;
    COFLAM[27] = -2.97906324E-02;
    COFLAM[28] = 1.71462501E+00;
    COFLAM[29] = -1.73273532E-01;
    COFLAM[30] = 3.32419185E-01;
    COFLAM[31] = -2.31119057E-02;
    COFLAM[32] = 5.56023763E-01;
    COFLAM[33] = 1.59073590E+00;
    COFLAM[34] = -5.28053839E-02;
    COFLAM[35] = 4.07601571E-04;
    COFLAM[36] = 2.28195645E+01;
    COFLAM[37] = -8.72278946E+00;
    COFLAM[38] = 1.49300487E+00;
    COFLAM[39] = -7.41524047E-02;
    COFLAM[40] = 9.92459822E+00;
    COFLAM[41] = -2.28318157E+00;
    COFLAM[42] = 4.73113746E-01;
    COFLAM[43] = -2.40056652E-02;
    COFLAM[44] = -3.01284291E+00;
    COFLAM[45] = 3.37554994E+00;
    COFLAM[46] = -3.43353119E-01;
    COFLAM[47] = 1.51043444E-02;
    COFLAM[48] = 6.27051982E-01;
    COFLAM[49] = 1.43139617E+00;
    COFLAM[50] = 1.80509282E-03;
    COFLAM[51] = -3.55624900E-03;
    COFLAM[52] = -1.24047589E+01;
    COFLAM[53] = 6.34783131E+00;
    COFLAM[54] = -6.37857884E-01;
    COFLAM[55] = 2.37613820E-02;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.02395222E+01;
    COFD[1] = 2.15403244E+00;
    COFD[2] = -6.97480266E-02;
    COFD[3] = 3.23666871E-03;
    COFD[4] = -1.11808682E+01;
    COFD[5] = 2.66936727E+00;
    COFD[6] = -1.34411514E-01;
    COFD[7] = 5.92957488E-03;
    COFD[8] = -1.20638601E+01;
    COFD[9] = 2.63303536E+00;
    COFD[10] = -1.29792632E-01;
    COFD[11] = 5.73363738E-03;
    COFD[12] = -1.13253458E+01;
    COFD[13] = 2.31195095E+00;
    COFD[14] = -8.63988037E-02;
    COFD[15] = 3.77573452E-03;
    COFD[16] = -9.86429034E+00;
    COFD[17] = 2.05348746E+00;
    COFD[18] = -5.90289007E-02;
    COFD[19] = 2.89596157E-03;
    COFD[20] = -1.06250182E+01;
    COFD[21] = 2.15849701E+00;
    COFD[22] = -6.53886401E-02;
    COFD[23] = 2.81453370E-03;
    COFD[24] = -1.06283453E+01;
    COFD[25] = 2.15849701E+00;
    COFD[26] = -6.53886401E-02;
    COFD[27] = 2.81453370E-03;
    COFD[28] = -1.57161204E+01;
    COFD[29] = 3.96062263E+00;
    COFD[30] = -2.98964970E-01;
    COFD[31] = 1.29322565E-02;
    COFD[32] = -1.15806808E+01;
    COFD[33] = 2.43235504E+00;
    COFD[34] = -1.02890179E-01;
    COFD[35] = 4.52903603E-03;
    COFD[36] = -1.68758926E+01;
    COFD[37] = 4.49460303E+00;
    COFD[38] = -3.64766132E-01;
    COFD[39] = 1.56457153E-02;
    COFD[40] = -1.13541075E+01;
    COFD[41] = 2.31999438E+00;
    COFD[42] = -8.75064804E-02;
    COFD[43] = 3.82656365E-03;
    COFD[44] = -1.15797750E+01;
    COFD[45] = 2.43235504E+00;
    COFD[46] = -1.02890179E-01;
    COFD[47] = 4.52903603E-03;
    COFD[48] = -1.15815344E+01;
    COFD[49] = 2.43235504E+00;
    COFD[50] = -1.02890179E-01;
    COFD[51] = 4.52903603E-03;
    COFD[52] = -1.35545239E+01;
    COFD[53] = 3.13878730E+00;
    COFD[54] = -1.94980335E-01;
    COFD[55] = 8.53744486E-03;
    COFD[56] = -1.11808682E+01;
    COFD[57] = 2.66936727E+00;
    COFD[58] = -1.34411514E-01;
    COFD[59] = 5.92957488E-03;
    COFD[60] = -1.43693056E+01;
    COFD[61] = 4.03992999E+00;
    COFD[62] = -3.08044800E-01;
    COFD[63] = 1.32757775E-02;
    COFD[64] = -1.51208119E+01;
    COFD[65] = 3.99904647E+00;
    COFD[66] = -3.03517220E-01;
    COFD[67] = 1.31117363E-02;
    COFD[68] = -1.40298830E+01;
    COFD[69] = 3.55837688E+00;
    COFD[70] = -2.47785790E-01;
    COFD[71] = 1.07555332E-02;
    COFD[72] = -9.71338331E+00;
    COFD[73] = 2.17561180E+00;
    COFD[74] = -7.28270090E-02;
    COFD[75] = 3.38302182E-03;
    COFD[76] = -1.31860117E+01;
    COFD[77] = 3.38003453E+00;
    COFD[78] = -2.25783856E-01;
    COFD[79] = 9.85028660E-03;
    COFD[80] = -1.31877711E+01;
    COFD[81] = 3.38003453E+00;
    COFD[82] = -2.25783856E-01;
    COFD[83] = 9.85028660E-03;
    COFD[84] = -1.95312724E+01;
    COFD[85] = 5.47046983E+00;
    COFD[86] = -4.74577605E-01;
    COFD[87] = 1.97408822E-02;
    COFD[88] = -1.43717529E+01;
    COFD[89] = 3.70920439E+00;
    COFD[90] = -2.67274113E-01;
    COFD[91] = 1.15967481E-02;
    COFD[92] = -1.93611051E+01;
    COFD[93] = 5.51579726E+00;
    COFD[94] = -4.76061961E-01;
    COFD[95] = 1.96329391E-02;
    COFD[96] = -1.40524065E+01;
    COFD[97] = 3.56261348E+00;
    COFD[98] = -2.48287981E-01;
    COFD[99] = 1.07752947E-02;
    COFD[100] = -1.43712864E+01;
    COFD[101] = 3.70920439E+00;
    COFD[102] = -2.67274113E-01;
    COFD[103] = 1.15967481E-02;
    COFD[104] = -1.43721922E+01;
    COFD[105] = 3.70920439E+00;
    COFD[106] = -2.67274113E-01;
    COFD[107] = 1.15967481E-02;
    COFD[108] = -1.72993972E+01;
    COFD[109] = 4.71931868E+00;
    COFD[110] = -3.91258152E-01;
    COFD[111] = 1.66866639E-02;
    COFD[112] = -1.20638601E+01;
    COFD[113] = 2.63303536E+00;
    COFD[114] = -1.29792632E-01;
    COFD[115] = 5.73363738E-03;
    COFD[116] = -1.51208119E+01;
    COFD[117] = 3.99904647E+00;
    COFD[118] = -3.03517220E-01;
    COFD[119] = 1.31117363E-02;
    COFD[120] = -1.68944722E+01;
    COFD[121] = 3.94346012E+00;
    COFD[122] = -2.96835271E-01;
    COFD[123] = 1.28438696E-02;
    COFD[124] = -1.57236706E+01;
    COFD[125] = 3.51447210E+00;
    COFD[126] = -2.42579007E-01;
    COFD[127] = 1.05506318E-02;
    COFD[128] = -1.08140177E+01;
    COFD[129] = 2.11737538E+00;
    COFD[130] = -6.46167749E-02;
    COFD[131] = 2.99827695E-03;
    COFD[132] = -1.47082523E+01;
    COFD[133] = 3.30683499E+00;
    COFD[134] = -2.16378602E-01;
    COFD[135] = 9.44670561E-03;
    COFD[136] = -1.47298720E+01;
    COFD[137] = 3.30683499E+00;
    COFD[138] = -2.16378602E-01;
    COFD[139] = 9.44670561E-03;
    COFD[140] = -2.11661448E+01;
    COFD[141] = 5.40762714E+00;
    COFD[142] = -4.67856822E-01;
    COFD[143] = 1.95051950E-02;
    COFD[144] = -1.59677692E+01;
    COFD[145] = 3.60186887E+00;
    COFD[146] = -2.53302622E-01;
    COFD[147] = 1.09893496E-02;
    COFD[148] = -2.10785324E+01;
    COFD[149] = 5.51573149E+00;
    COFD[150] = -4.78177665E-01;
    COFD[151] = 1.98082796E-02;
    COFD[152] = -1.57440433E+01;
    COFD[153] = 3.51861272E+00;
    COFD[154] = -2.43068621E-01;
    COFD[155] = 1.05698368E-02;
    COFD[156] = -1.59592184E+01;
    COFD[157] = 3.60186887E+00;
    COFD[158] = -2.53302622E-01;
    COFD[159] = 1.09893496E-02;
    COFD[160] = -1.59759490E+01;
    COFD[161] = 3.60186887E+00;
    COFD[162] = -2.53302622E-01;
    COFD[163] = 1.09893496E-02;
    COFD[164] = -1.90183510E+01;
    COFD[165] = 4.64763677E+00;
    COFD[166] = -3.82799418E-01;
    COFD[167] = 1.63539171E-02;
    COFD[168] = -1.13253458E+01;
    COFD[169] = 2.31195095E+00;
    COFD[170] = -8.63988037E-02;
    COFD[171] = 3.77573452E-03;
    COFD[172] = -1.40298830E+01;
    COFD[173] = 3.55837688E+00;
    COFD[174] = -2.47785790E-01;
    COFD[175] = 1.07555332E-02;
    COFD[176] = -1.57236706E+01;
    COFD[177] = 3.51447210E+00;
    COFD[178] = -2.42579007E-01;
    COFD[179] = 1.05506318E-02;
    COFD[180] = -1.47639290E+01;
    COFD[181] = 3.15955654E+00;
    COFD[182] = -1.97590757E-01;
    COFD[183] = 8.64692156E-03;
    COFD[184] = -1.01976409E+01;
    COFD[185] = 1.83188320E+00;
    COFD[186] = -2.40547456E-02;
    COFD[187] = 1.08399898E-03;
    COFD[188] = -1.38756407E+01;
    COFD[189] = 2.98558426E+00;
    COFD[190] = -1.75507216E-01;
    COFD[191] = 7.71173691E-03;
    COFD[192] = -1.38948667E+01;
    COFD[193] = 2.98558426E+00;
    COFD[194] = -1.75507216E-01;
    COFD[195] = 7.71173691E-03;
    COFD[196] = -2.02052895E+01;
    COFD[197] = 5.10993120E+00;
    COFD[198] = -4.36931630E-01;
    COFD[199] = 1.84677592E-02;
    COFD[200] = -1.50168028E+01;
    COFD[201] = 3.25515933E+00;
    COFD[202] = -2.09710110E-01;
    COFD[203] = 9.15941830E-03;
    COFD[204] = -2.08123325E+01;
    COFD[205] = 5.42470154E+00;
    COFD[206] = -4.69700416E-01;
    COFD[207] = 1.95706904E-02;
    COFD[208] = -1.47850486E+01;
    COFD[209] = 3.16433919E+00;
    COFD[210] = -1.98191564E-01;
    COFD[211] = 8.67209742E-03;
    COFD[212] = -1.50096240E+01;
    COFD[213] = 3.25515933E+00;
    COFD[214] = -2.09710110E-01;
    COFD[215] = 9.15941830E-03;
    COFD[216] = -1.50236516E+01;
    COFD[217] = 3.25515933E+00;
    COFD[218] = -2.09710110E-01;
    COFD[219] = 9.15941830E-03;
    COFD[220] = -1.77350592E+01;
    COFD[221] = 4.19328271E+00;
    COFD[222] = -3.26911461E-01;
    COFD[223] = 1.40520357E-02;
    COFD[224] = -9.86429034E+00;
    COFD[225] = 2.05348746E+00;
    COFD[226] = -5.90289007E-02;
    COFD[227] = 2.89596157E-03;
    COFD[228] = -9.71338331E+00;
    COFD[229] = 2.17561180E+00;
    COFD[230] = -7.28270090E-02;
    COFD[231] = 3.38302182E-03;
    COFD[232] = -1.08140177E+01;
    COFD[233] = 2.11737538E+00;
    COFD[234] = -6.46167749E-02;
    COFD[235] = 2.99827695E-03;
    COFD[236] = -1.01976409E+01;
    COFD[237] = 1.83188320E+00;
    COFD[238] = -2.40547456E-02;
    COFD[239] = 1.08399898E-03;
    COFD[240] = -7.72963289E+00;
    COFD[241] = 1.13864728E+00;
    COFD[242] = 7.22991035E-02;
    COFD[243] = -3.32491895E-03;
    COFD[244] = -9.70779324E+00;
    COFD[245] = 1.77912272E+00;
    COFD[246] = -1.67349571E-02;
    COFD[247] = 7.45446845E-04;
    COFD[248] = -9.71375861E+00;
    COFD[249] = 1.77912272E+00;
    COFD[250] = -1.67349571E-02;
    COFD[251] = 7.45446845E-04;
    COFD[252] = -1.23116531E+01;
    COFD[253] = 2.62312681E+00;
    COFD[254] = -1.28556314E-01;
    COFD[255] = 5.68221707E-03;
    COFD[256] = -1.03327323E+01;
    COFD[257] = 1.90522472E+00;
    COFD[258] = -3.44812795E-02;
    COFD[259] = 1.57640018E-03;
    COFD[260] = -1.21950642E+01;
    COFD[261] = 2.72222246E+00;
    COFD[262] = -1.41335602E-01;
    COFD[263] = 6.23222872E-03;
    COFD[264] = -1.02057322E+01;
    COFD[265] = 1.83104667E+00;
    COFD[266] = -2.39235907E-02;
    COFD[267] = 1.07741763E-03;
    COFD[268] = -1.03310318E+01;
    COFD[269] = 1.90522472E+00;
    COFD[270] = -3.44812795E-02;
    COFD[271] = 1.57640018E-03;
    COFD[272] = -1.03343373E+01;
    COFD[273] = 1.90522472E+00;
    COFD[274] = -3.44812795E-02;
    COFD[275] = 1.57640018E-03;
    COFD[276] = -1.09328506E+01;
    COFD[277] = 2.05651569E+00;
    COFD[278] = -5.19591463E-02;
    COFD[279] = 2.22384771E-03;
    COFD[280] = -1.06250182E+01;
    COFD[281] = 2.15849701E+00;
    COFD[282] = -6.53886401E-02;
    COFD[283] = 2.81453370E-03;
    COFD[284] = -1.31860117E+01;
    COFD[285] = 3.38003453E+00;
    COFD[286] = -2.25783856E-01;
    COFD[287] = 9.85028660E-03;
    COFD[288] = -1.47082523E+01;
    COFD[289] = 3.30683499E+00;
    COFD[290] = -2.16378602E-01;
    COFD[291] = 9.44670561E-03;
    COFD[292] = -1.38756407E+01;
    COFD[293] = 2.98558426E+00;
    COFD[294] = -1.75507216E-01;
    COFD[295] = 7.71173691E-03;
    COFD[296] = -9.70779324E+00;
    COFD[297] = 1.77912272E+00;
    COFD[298] = -1.67349571E-02;
    COFD[299] = 7.45446845E-04;
    COFD[300] = -1.29877365E+01;
    COFD[301] = 2.80841511E+00;
    COFD[302] = -1.52629888E-01;
    COFD[303] = 6.72604927E-03;
    COFD[304] = -1.30027772E+01;
    COFD[305] = 2.80841511E+00;
    COFD[306] = -1.52629888E-01;
    COFD[307] = 6.72604927E-03;
    COFD[308] = -1.91045619E+01;
    COFD[309] = 4.87977047E+00;
    COFD[310] = -4.10448693E-01;
    COFD[311] = 1.74535827E-02;
    COFD[312] = -1.40916052E+01;
    COFD[313] = 3.07458927E+00;
    COFD[314] = -1.86899591E-01;
    COFD[315] = 8.19829781E-03;
    COFD[316] = -1.91096797E+01;
    COFD[317] = 5.02608697E+00;
    COFD[318] = -4.26959993E-01;
    COFD[319] = 1.80709910E-02;
    COFD[320] = -1.39007410E+01;
    COFD[321] = 2.99164244E+00;
    COFD[322] = -1.76293106E-01;
    COFD[323] = 7.74575100E-03;
    COFD[324] = -1.40864894E+01;
    COFD[325] = 3.07458927E+00;
    COFD[326] = -1.86899591E-01;
    COFD[327] = 8.19829781E-03;
    COFD[328] = -1.40964661E+01;
    COFD[329] = 3.07458927E+00;
    COFD[330] = -1.86899591E-01;
    COFD[331] = 8.19829781E-03;
    COFD[332] = -1.67115577E+01;
    COFD[333] = 3.98859394E+00;
    COFD[334] = -3.02316219E-01;
    COFD[335] = 1.30661099E-02;
    COFD[336] = -1.06283453E+01;
    COFD[337] = 2.15849701E+00;
    COFD[338] = -6.53886401E-02;
    COFD[339] = 2.81453370E-03;
    COFD[340] = -1.31877711E+01;
    COFD[341] = 3.38003453E+00;
    COFD[342] = -2.25783856E-01;
    COFD[343] = 9.85028660E-03;
    COFD[344] = -1.47298720E+01;
    COFD[345] = 3.30683499E+00;
    COFD[346] = -2.16378602E-01;
    COFD[347] = 9.44670561E-03;
    COFD[348] = -1.38948667E+01;
    COFD[349] = 2.98558426E+00;
    COFD[350] = -1.75507216E-01;
    COFD[351] = 7.71173691E-03;
    COFD[352] = -9.71375861E+00;
    COFD[353] = 1.77912272E+00;
    COFD[354] = -1.67349571E-02;
    COFD[355] = 7.45446845E-04;
    COFD[356] = -1.30027772E+01;
    COFD[357] = 2.80841511E+00;
    COFD[358] = -1.52629888E-01;
    COFD[359] = 6.72604927E-03;
    COFD[360] = -1.30182843E+01;
    COFD[361] = 2.80841511E+00;
    COFD[362] = -1.52629888E-01;
    COFD[363] = 6.72604927E-03;
    COFD[364] = -1.91240379E+01;
    COFD[365] = 4.87977047E+00;
    COFD[366] = -4.10448693E-01;
    COFD[367] = 1.74535827E-02;
    COFD[368] = -1.41119732E+01;
    COFD[369] = 3.07458927E+00;
    COFD[370] = -1.86899591E-01;
    COFD[371] = 8.19829781E-03;
    COFD[372] = -1.91256261E+01;
    COFD[373] = 5.02608697E+00;
    COFD[374] = -4.26959993E-01;
    COFD[375] = 1.80709910E-02;
    COFD[376] = -1.39199663E+01;
    COFD[377] = 2.99164244E+00;
    COFD[378] = -1.76293106E-01;
    COFD[379] = 7.74575100E-03;
    COFD[380] = -1.41066459E+01;
    COFD[381] = 3.07458927E+00;
    COFD[382] = -1.86899591E-01;
    COFD[383] = 8.19829781E-03;
    COFD[384] = -1.41170372E+01;
    COFD[385] = 3.07458927E+00;
    COFD[386] = -1.86899591E-01;
    COFD[387] = 8.19829781E-03;
    COFD[388] = -1.67337768E+01;
    COFD[389] = 3.98859394E+00;
    COFD[390] = -3.02316219E-01;
    COFD[391] = 1.30661099E-02;
    COFD[392] = -1.57161204E+01;
    COFD[393] = 3.96062263E+00;
    COFD[394] = -2.98964970E-01;
    COFD[395] = 1.29322565E-02;
    COFD[396] = -1.95312724E+01;
    COFD[397] = 5.47046983E+00;
    COFD[398] = -4.74577605E-01;
    COFD[399] = 1.97408822E-02;
    COFD[400] = -2.11661448E+01;
    COFD[401] = 5.40762714E+00;
    COFD[402] = -4.67856822E-01;
    COFD[403] = 1.95051950E-02;
    COFD[404] = -2.02052895E+01;
    COFD[405] = 5.10993120E+00;
    COFD[406] = -4.36931630E-01;
    COFD[407] = 1.84677592E-02;
    COFD[408] = -1.23116531E+01;
    COFD[409] = 2.62312681E+00;
    COFD[410] = -1.28556314E-01;
    COFD[411] = 5.68221707E-03;
    COFD[412] = -1.91045619E+01;
    COFD[413] = 4.87977047E+00;
    COFD[414] = -4.10448693E-01;
    COFD[415] = 1.74535827E-02;
    COFD[416] = -1.91240379E+01;
    COFD[417] = 4.87977047E+00;
    COFD[418] = -4.10448693E-01;
    COFD[419] = 1.74535827E-02;
    COFD[420] = -1.98983761E+01;
    COFD[421] = 4.38041133E+00;
    COFD[422] = -2.77538214E-01;
    COFD[423] = 9.06748822E-03;
    COFD[424] = -2.05180636E+01;
    COFD[425] = 5.21473296E+00;
    COFD[426] = -4.48646311E-01;
    COFD[427] = 1.89013813E-02;
    COFD[428] = -1.87171417E+01;
    COFD[429] = 4.00967621E+00;
    COFD[430] = -2.21153539E-01;
    COFD[431] = 6.31528745E-03;
    COFD[432] = -2.02361862E+01;
    COFD[433] = 5.11785645E+00;
    COFD[434] = -4.37867828E-01;
    COFD[435] = 1.85047543E-02;
    COFD[436] = -2.05107487E+01;
    COFD[437] = 5.21473296E+00;
    COFD[438] = -4.48646311E-01;
    COFD[439] = 1.89013813E-02;
    COFD[440] = -2.05250441E+01;
    COFD[441] = 5.21473296E+00;
    COFD[442] = -4.48646311E-01;
    COFD[443] = 1.89013813E-02;
    COFD[444] = -2.20758752E+01;
    COFD[445] = 5.52171573E+00;
    COFD[446] = -4.63284984E-01;
    COFD[447] = 1.85570924E-02;
    COFD[448] = -1.15806808E+01;
    COFD[449] = 2.43235504E+00;
    COFD[450] = -1.02890179E-01;
    COFD[451] = 4.52903603E-03;
    COFD[452] = -1.43717529E+01;
    COFD[453] = 3.70920439E+00;
    COFD[454] = -2.67274113E-01;
    COFD[455] = 1.15967481E-02;
    COFD[456] = -1.59677692E+01;
    COFD[457] = 3.60186887E+00;
    COFD[458] = -2.53302622E-01;
    COFD[459] = 1.09893496E-02;
    COFD[460] = -1.50168028E+01;
    COFD[461] = 3.25515933E+00;
    COFD[462] = -2.09710110E-01;
    COFD[463] = 9.15941830E-03;
    COFD[464] = -1.03327323E+01;
    COFD[465] = 1.90522472E+00;
    COFD[466] = -3.44812795E-02;
    COFD[467] = 1.57640018E-03;
    COFD[468] = -1.40916052E+01;
    COFD[469] = 3.07458927E+00;
    COFD[470] = -1.86899591E-01;
    COFD[471] = 8.19829781E-03;
    COFD[472] = -1.41119732E+01;
    COFD[473] = 3.07458927E+00;
    COFD[474] = -1.86899591E-01;
    COFD[475] = 8.19829781E-03;
    COFD[476] = -2.05180636E+01;
    COFD[477] = 5.21473296E+00;
    COFD[478] = -4.48646311E-01;
    COFD[479] = 1.89013813E-02;
    COFD[480] = -1.53265780E+01;
    COFD[481] = 3.37317428E+00;
    COFD[482] = -2.24900439E-01;
    COFD[483] = 9.81228151E-03;
    COFD[484] = -2.04177482E+01;
    COFD[485] = 5.31457079E+00;
    COFD[486] = -4.58216496E-01;
    COFD[487] = 1.91825910E-02;
    COFD[488] = -1.50443569E+01;
    COFD[489] = 3.26249588E+00;
    COFD[490] = -2.10658287E-01;
    COFD[491] = 9.20032462E-03;
    COFD[492] = -1.53187643E+01;
    COFD[493] = 3.37317428E+00;
    COFD[494] = -2.24900439E-01;
    COFD[495] = 9.81228151E-03;
    COFD[496] = -1.53340417E+01;
    COFD[497] = 3.37317428E+00;
    COFD[498] = -2.24900439E-01;
    COFD[499] = 9.81228151E-03;
    COFD[500] = -1.81286555E+01;
    COFD[501] = 4.33684042E+00;
    COFD[502] = -3.44981265E-01;
    COFD[503] = 1.48142449E-02;
    COFD[504] = -1.68758926E+01;
    COFD[505] = 4.49460303E+00;
    COFD[506] = -3.64766132E-01;
    COFD[507] = 1.56457153E-02;
    COFD[508] = -1.93611051E+01;
    COFD[509] = 5.51579726E+00;
    COFD[510] = -4.76061961E-01;
    COFD[511] = 1.96329391E-02;
    COFD[512] = -2.10785324E+01;
    COFD[513] = 5.51573149E+00;
    COFD[514] = -4.78177665E-01;
    COFD[515] = 1.98082796E-02;
    COFD[516] = -2.08123325E+01;
    COFD[517] = 5.42470154E+00;
    COFD[518] = -4.69700416E-01;
    COFD[519] = 1.95706904E-02;
    COFD[520] = -1.21950642E+01;
    COFD[521] = 2.72222246E+00;
    COFD[522] = -1.41335602E-01;
    COFD[523] = 6.23222872E-03;
    COFD[524] = -1.91096797E+01;
    COFD[525] = 5.02608697E+00;
    COFD[526] = -4.26959993E-01;
    COFD[527] = 1.80709910E-02;
    COFD[528] = -1.91256261E+01;
    COFD[529] = 5.02608697E+00;
    COFD[530] = -4.26959993E-01;
    COFD[531] = 1.80709910E-02;
    COFD[532] = -1.87171417E+01;
    COFD[533] = 4.00967621E+00;
    COFD[534] = -2.21153539E-01;
    COFD[535] = 6.31528745E-03;
    COFD[536] = -2.04177482E+01;
    COFD[537] = 5.31457079E+00;
    COFD[538] = -4.58216496E-01;
    COFD[539] = 1.91825910E-02;
    COFD[540] = -1.31492641E+01;
    COFD[541] = 1.48004311E+00;
    COFD[542] = 1.60499553E-01;
    COFD[543] = -1.19765679E-02;
    COFD[544] = -2.08943798E+01;
    COFD[545] = 5.44718652E+00;
    COFD[546] = -4.72082953E-01;
    COFD[547] = 1.96531321E-02;
    COFD[548] = -2.10640014E+01;
    COFD[549] = 5.50980695E+00;
    COFD[550] = -4.78335488E-01;
    COFD[551] = 1.98515434E-02;
    COFD[552] = -2.04230073E+01;
    COFD[553] = 5.31457079E+00;
    COFD[554] = -4.58216496E-01;
    COFD[555] = 1.91825910E-02;
    COFD[556] = -2.12021508E+01;
    COFD[557] = 5.20775052E+00;
    COFD[558] = -4.07348327E-01;
    COFD[559] = 1.55473283E-02;
    COFD[560] = -1.13541075E+01;
    COFD[561] = 2.31999438E+00;
    COFD[562] = -8.75064804E-02;
    COFD[563] = 3.82656365E-03;
    COFD[564] = -1.40524065E+01;
    COFD[565] = 3.56261348E+00;
    COFD[566] = -2.48287981E-01;
    COFD[567] = 1.07752947E-02;
    COFD[568] = -1.57440433E+01;
    COFD[569] = 3.51861272E+00;
    COFD[570] = -2.43068621E-01;
    COFD[571] = 1.05698368E-02;
    COFD[572] = -1.47850486E+01;
    COFD[573] = 3.16433919E+00;
    COFD[574] = -1.98191564E-01;
    COFD[575] = 8.67209742E-03;
    COFD[576] = -1.02057322E+01;
    COFD[577] = 1.83104667E+00;
    COFD[578] = -2.39235907E-02;
    COFD[579] = 1.07741763E-03;
    COFD[580] = -1.39007410E+01;
    COFD[581] = 2.99164244E+00;
    COFD[582] = -1.76293106E-01;
    COFD[583] = 7.74575100E-03;
    COFD[584] = -1.39199663E+01;
    COFD[585] = 2.99164244E+00;
    COFD[586] = -1.76293106E-01;
    COFD[587] = 7.74575100E-03;
    COFD[588] = -2.02361862E+01;
    COFD[589] = 5.11785645E+00;
    COFD[590] = -4.37867828E-01;
    COFD[591] = 1.85047543E-02;
    COFD[592] = -1.50443569E+01;
    COFD[593] = 3.26249588E+00;
    COFD[594] = -2.10658287E-01;
    COFD[595] = 9.20032462E-03;
    COFD[596] = -2.08943798E+01;
    COFD[597] = 5.44718652E+00;
    COFD[598] = -4.72082953E-01;
    COFD[599] = 1.96531321E-02;
    COFD[600] = -1.48061490E+01;
    COFD[601] = 3.16912473E+00;
    COFD[602] = -1.98792456E-01;
    COFD[603] = 8.69726395E-03;
    COFD[604] = -1.50371784E+01;
    COFD[605] = 3.26249588E+00;
    COFD[606] = -2.10658287E-01;
    COFD[607] = 9.20032462E-03;
    COFD[608] = -1.50512053E+01;
    COFD[609] = 3.26249588E+00;
    COFD[610] = -2.10658287E-01;
    COFD[611] = 9.20032462E-03;
    COFD[612] = -1.77673000E+01;
    COFD[613] = 4.20234040E+00;
    COFD[614] = -3.28057658E-01;
    COFD[615] = 1.41006192E-02;
    COFD[616] = -1.15797750E+01;
    COFD[617] = 2.43235504E+00;
    COFD[618] = -1.02890179E-01;
    COFD[619] = 4.52903603E-03;
    COFD[620] = -1.43712864E+01;
    COFD[621] = 3.70920439E+00;
    COFD[622] = -2.67274113E-01;
    COFD[623] = 1.15967481E-02;
    COFD[624] = -1.59592184E+01;
    COFD[625] = 3.60186887E+00;
    COFD[626] = -2.53302622E-01;
    COFD[627] = 1.09893496E-02;
    COFD[628] = -1.50096240E+01;
    COFD[629] = 3.25515933E+00;
    COFD[630] = -2.09710110E-01;
    COFD[631] = 9.15941830E-03;
    COFD[632] = -1.03310318E+01;
    COFD[633] = 1.90522472E+00;
    COFD[634] = -3.44812795E-02;
    COFD[635] = 1.57640018E-03;
    COFD[636] = -1.40864894E+01;
    COFD[637] = 3.07458927E+00;
    COFD[638] = -1.86899591E-01;
    COFD[639] = 8.19829781E-03;
    COFD[640] = -1.41066459E+01;
    COFD[641] = 3.07458927E+00;
    COFD[642] = -1.86899591E-01;
    COFD[643] = 8.19829781E-03;
    COFD[644] = -2.05107487E+01;
    COFD[645] = 5.21473296E+00;
    COFD[646] = -4.48646311E-01;
    COFD[647] = 1.89013813E-02;
    COFD[648] = -1.53187643E+01;
    COFD[649] = 3.37317428E+00;
    COFD[650] = -2.24900439E-01;
    COFD[651] = 9.81228151E-03;
    COFD[652] = -2.10640014E+01;
    COFD[653] = 5.50980695E+00;
    COFD[654] = -4.78335488E-01;
    COFD[655] = 1.98515434E-02;
    COFD[656] = -1.50371784E+01;
    COFD[657] = 3.26249588E+00;
    COFD[658] = -2.10658287E-01;
    COFD[659] = 9.20032462E-03;
    COFD[660] = -1.53110708E+01;
    COFD[661] = 3.37317428E+00;
    COFD[662] = -2.24900439E-01;
    COFD[663] = 9.81228151E-03;
    COFD[664] = -1.53261114E+01;
    COFD[665] = 3.37317428E+00;
    COFD[666] = -2.24900439E-01;
    COFD[667] = 9.81228151E-03;
    COFD[668] = -1.81197354E+01;
    COFD[669] = 4.33684042E+00;
    COFD[670] = -3.44981265E-01;
    COFD[671] = 1.48142449E-02;
    COFD[672] = -1.15815344E+01;
    COFD[673] = 2.43235504E+00;
    COFD[674] = -1.02890179E-01;
    COFD[675] = 4.52903603E-03;
    COFD[676] = -1.43721922E+01;
    COFD[677] = 3.70920439E+00;
    COFD[678] = -2.67274113E-01;
    COFD[679] = 1.15967481E-02;
    COFD[680] = -1.59759490E+01;
    COFD[681] = 3.60186887E+00;
    COFD[682] = -2.53302622E-01;
    COFD[683] = 1.09893496E-02;
    COFD[684] = -1.50236516E+01;
    COFD[685] = 3.25515933E+00;
    COFD[686] = -2.09710110E-01;
    COFD[687] = 9.15941830E-03;
    COFD[688] = -1.03343373E+01;
    COFD[689] = 1.90522472E+00;
    COFD[690] = -3.44812795E-02;
    COFD[691] = 1.57640018E-03;
    COFD[692] = -1.40964661E+01;
    COFD[693] = 3.07458927E+00;
    COFD[694] = -1.86899591E-01;
    COFD[695] = 8.19829781E-03;
    COFD[696] = -1.41170372E+01;
    COFD[697] = 3.07458927E+00;
    COFD[698] = -1.86899591E-01;
    COFD[699] = 8.19829781E-03;
    COFD[700] = -2.05250441E+01;
    COFD[701] = 5.21473296E+00;
    COFD[702] = -4.48646311E-01;
    COFD[703] = 1.89013813E-02;
    COFD[704] = -1.53340417E+01;
    COFD[705] = 3.37317428E+00;
    COFD[706] = -2.24900439E-01;
    COFD[707] = 9.81228151E-03;
    COFD[708] = -2.04230073E+01;
    COFD[709] = 5.31457079E+00;
    COFD[710] = -4.58216496E-01;
    COFD[711] = 1.91825910E-02;
    COFD[712] = -1.50512053E+01;
    COFD[713] = 3.26249588E+00;
    COFD[714] = -2.10658287E-01;
    COFD[715] = 9.20032462E-03;
    COFD[716] = -1.53261114E+01;
    COFD[717] = 3.37317428E+00;
    COFD[718] = -2.24900439E-01;
    COFD[719] = 9.81228151E-03;
    COFD[720] = -1.53416186E+01;
    COFD[721] = 3.37317428E+00;
    COFD[722] = -2.24900439E-01;
    COFD[723] = 9.81228151E-03;
    COFD[724] = -1.81371948E+01;
    COFD[725] = 4.33684042E+00;
    COFD[726] = -3.44981265E-01;
    COFD[727] = 1.48142449E-02;
    COFD[728] = -1.35545239E+01;
    COFD[729] = 3.13878730E+00;
    COFD[730] = -1.94980335E-01;
    COFD[731] = 8.53744486E-03;
    COFD[732] = -1.72993972E+01;
    COFD[733] = 4.71931868E+00;
    COFD[734] = -3.91258152E-01;
    COFD[735] = 1.66866639E-02;
    COFD[736] = -1.90183510E+01;
    COFD[737] = 4.64763677E+00;
    COFD[738] = -3.82799418E-01;
    COFD[739] = 1.63539171E-02;
    COFD[740] = -1.77350592E+01;
    COFD[741] = 4.19328271E+00;
    COFD[742] = -3.26911461E-01;
    COFD[743] = 1.40520357E-02;
    COFD[744] = -1.09328506E+01;
    COFD[745] = 2.05651569E+00;
    COFD[746] = -5.19591463E-02;
    COFD[747] = 2.22384771E-03;
    COFD[748] = -1.67115577E+01;
    COFD[749] = 3.98859394E+00;
    COFD[750] = -3.02316219E-01;
    COFD[751] = 1.30661099E-02;
    COFD[752] = -1.67337768E+01;
    COFD[753] = 3.98859394E+00;
    COFD[754] = -3.02316219E-01;
    COFD[755] = 1.30661099E-02;
    COFD[756] = -2.20758752E+01;
    COFD[757] = 5.52171573E+00;
    COFD[758] = -4.63284984E-01;
    COFD[759] = 1.85570924E-02;
    COFD[760] = -1.81286555E+01;
    COFD[761] = 4.33684042E+00;
    COFD[762] = -3.44981265E-01;
    COFD[763] = 1.48142449E-02;
    COFD[764] = -2.12021508E+01;
    COFD[765] = 5.20775052E+00;
    COFD[766] = -4.07348327E-01;
    COFD[767] = 1.55473283E-02;
    COFD[768] = -1.77673000E+01;
    COFD[769] = 4.20234040E+00;
    COFD[770] = -3.28057658E-01;
    COFD[771] = 1.41006192E-02;
    COFD[772] = -1.81197354E+01;
    COFD[773] = 4.33684042E+00;
    COFD[774] = -3.44981265E-01;
    COFD[775] = 1.48142449E-02;
    COFD[776] = -1.81371948E+01;
    COFD[777] = 4.33684042E+00;
    COFD[778] = -3.44981265E-01;
    COFD[779] = 1.48142449E-02;
    COFD[780] = -2.10907727E+01;
    COFD[781] = 5.29211327E+00;
    COFD[782] = -4.56068366E-01;
    COFD[783] = 1.91195062E-02;
}


/*List of specs with small weight, dim NLITE */
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 0;
    KTDIF[1] = 1;
    KTDIF[2] = 4;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(amrex::Real* COFTD) {
    COFTD[0] = 0.00000000E+00;
    COFTD[1] = 0.00000000E+00;
    COFTD[2] = 0.00000000E+00;
    COFTD[3] = 0.00000000E+00;
    COFTD[4] = -1.52534742E-01;
    COFTD[5] = -5.46404022E-05;
    COFTD[6] = 2.93412470E-08;
    COFTD[7] = -4.87091914E-12;
    COFTD[8] = 4.22530228E-01;
    COFTD[9] = 1.32084268E-04;
    COFTD[10] = -7.12222323E-08;
    COFTD[11] = 1.19516090E-11;
    COFTD[12] = 4.45261966E-01;
    COFTD[13] = 4.94697174E-05;
    COFTD[14] = -2.63023442E-08;
    COFTD[15] = 4.90306217E-12;
    COFTD[16] = 1.61613664E-01;
    COFTD[17] = 4.74155340E-05;
    COFTD[18] = -1.67115247E-08;
    COFTD[19] = -1.88982125E-12;
    COFTD[20] = 4.15583337E-01;
    COFTD[21] = 1.09738399E-05;
    COFTD[22] = -3.96021963E-09;
    COFTD[23] = 1.14414443E-12;
    COFTD[24] = 4.21932443E-01;
    COFTD[25] = 1.11414935E-05;
    COFTD[26] = -4.02072219E-09;
    COFTD[27] = 1.16162418E-12;
    COFTD[28] = 1.61392287E-01;
    COFTD[29] = 5.01084129E-04;
    COFTD[30] = -2.38273963E-07;
    COFTD[31] = 3.49344424E-11;
    COFTD[32] = 4.44452569E-01;
    COFTD[33] = 7.14525507E-05;
    COFTD[34] = -3.86257187E-08;
    COFTD[35] = 6.88979640E-12;
    COFTD[36] = 6.02028221E-02;
    COFTD[37] = 5.61561867E-04;
    COFTD[38] = -2.55372862E-07;
    COFTD[39] = 3.63389913E-11;
    COFTD[40] = 4.44653617E-01;
    COFTD[41] = 5.06631704E-05;
    COFTD[42] = -2.69820900E-08;
    COFTD[43] = 5.01289759E-12;
    COFTD[44] = 4.42739084E-01;
    COFTD[45] = 7.11770818E-05;
    COFTD[46] = -3.84768062E-08;
    COFTD[47] = 6.86323437E-12;
    COFTD[48] = 4.46070183E-01;
    COFTD[49] = 7.17126069E-05;
    COFTD[50] = -3.87662996E-08;
    COFTD[51] = 6.91487226E-12;
    COFTD[52] = 3.25742450E-01;
    COFTD[53] = 3.03633411E-04;
    COFTD[54] = -1.55290330E-07;
    COFTD[55] = 2.41466436E-11;
    COFTD[56] = 1.52534742E-01;
    COFTD[57] = 5.46404022E-05;
    COFTD[58] = -2.93412470E-08;
    COFTD[59] = 4.87091914E-12;
    COFTD[60] = 0.00000000E+00;
    COFTD[61] = 0.00000000E+00;
    COFTD[62] = 0.00000000E+00;
    COFTD[63] = 0.00000000E+00;
    COFTD[64] = 1.65429221E-01;
    COFTD[65] = 5.61238922E-04;
    COFTD[66] = -2.65650544E-07;
    COFTD[67] = 3.88229592E-11;
    COFTD[68] = 2.40744421E-01;
    COFTD[69] = 4.45343451E-04;
    COFTD[70] = -2.18173874E-07;
    COFTD[71] = 3.26958506E-11;
    COFTD[72] = 3.40762433E-01;
    COFTD[73] = -4.04057756E-05;
    COFTD[74] = 3.27879533E-08;
    COFTD[75] = -6.27093812E-12;
    COFTD[76] = 2.70010150E-01;
    COFTD[77] = 3.61555093E-04;
    COFTD[78] = -1.80744752E-07;
    COFTD[79] = 2.75321248E-11;
    COFTD[80] = 2.72041664E-01;
    COFTD[81] = 3.64275376E-04;
    COFTD[82] = -1.82104647E-07;
    COFTD[83] = 2.77392722E-11;
    COFTD[84] = -1.24647991E-01;
    COFTD[85] = 7.96525614E-04;
    COFTD[86] = -3.24998782E-07;
    COFTD[87] = 4.32516736E-11;
    COFTD[88] = 2.20907853E-01;
    COFTD[89] = 4.81089870E-04;
    COFTD[90] = -2.33376944E-07;
    COFTD[91] = 3.47138305E-11;
    COFTD[92] = -1.41883744E-01;
    COFTD[93] = 7.66558810E-04;
    COFTD[94] = -3.06550003E-07;
    COFTD[95] = 4.02959502E-11;
    COFTD[96] = 2.39409939E-01;
    COFTD[97] = 4.47197179E-04;
    COFTD[98] = -2.18951702E-07;
    COFTD[99] = 3.27973510E-11;
    COFTD[100] = 2.20482843E-01;
    COFTD[101] = 4.80164288E-04;
    COFTD[102] = -2.32927944E-07;
    COFTD[103] = 3.46470436E-11;
    COFTD[104] = 2.21308399E-01;
    COFTD[105] = 4.81962174E-04;
    COFTD[106] = -2.33800100E-07;
    COFTD[107] = 3.47767730E-11;
    COFTD[108] = 2.44369385E-02;
    COFTD[109] = 7.18242498E-04;
    COFTD[110] = -3.19718504E-07;
    COFTD[111] = 4.48828685E-11;
    COFTD[112] = -1.61613664E-01;
    COFTD[113] = -4.74155340E-05;
    COFTD[114] = 1.67115247E-08;
    COFTD[115] = 1.88982125E-12;
    COFTD[116] = -3.40762433E-01;
    COFTD[117] = 4.04057756E-05;
    COFTD[118] = -3.27879533E-08;
    COFTD[119] = 6.27093812E-12;
    COFTD[120] = 4.66315159E-01;
    COFTD[121] = -5.60150425E-05;
    COFTD[122] = 4.65987669E-08;
    COFTD[123] = -9.13646318E-12;
    COFTD[124] = 4.22009934E-01;
    COFTD[125] = -4.14042334E-05;
    COFTD[126] = 4.38751613E-08;
    COFTD[127] = -1.02860246E-11;
    COFTD[128] = 0.00000000E+00;
    COFTD[129] = 0.00000000E+00;
    COFTD[130] = 0.00000000E+00;
    COFTD[131] = 0.00000000E+00;
    COFTD[132] = 3.31587939E-01;
    COFTD[133] = -1.96388078E-05;
    COFTD[134] = 3.02388828E-08;
    COFTD[135] = -8.44998018E-12;
    COFTD[136] = 3.42203127E-01;
    COFTD[137] = -2.02675087E-05;
    COFTD[138] = 3.12069259E-08;
    COFTD[139] = -8.72049099E-12;
    COFTD[140] = 3.56629042E-01;
    COFTD[141] = 1.06116875E-04;
    COFTD[142] = -5.72917265E-08;
    COFTD[143] = 9.65264597E-12;
    COFTD[144] = 4.43649137E-01;
    COFTD[145] = -4.87484458E-05;
    COFTD[146] = 4.69718656E-08;
    COFTD[147] = -1.03568760E-11;
    COFTD[148] = 2.84983505E-01;
    COFTD[149] = 1.15460005E-04;
    COFTD[150] = -6.17197869E-08;
    COFTD[151] = 1.01504212E-11;
    COFTD[152] = 4.22171414E-01;
    COFTD[153] = -4.17749918E-05;
    COFTD[154] = 4.39726219E-08;
    COFTD[155] = -1.02672932E-11;
    COFTD[156] = 4.40220831E-01;
    COFTD[157] = -4.83717413E-05;
    COFTD[158] = 4.66088897E-08;
    COFTD[159] = -1.02768430E-11;
    COFTD[160] = 4.46895651E-01;
    COFTD[161] = -4.91051748E-05;
    COFTD[162] = 4.73155940E-08;
    COFTD[163] = -1.04326650E-11;
    COFTD[164] = 4.59663274E-01;
    COFTD[165] = -1.74770868E-05;
    COFTD[166] = 1.42888118E-08;
    COFTD[167] = -2.03610705E-12;
}

#ifndef AMREX_USE_GPU
/* Initializes parameter database */
void CKINIT()
{

    // (0):  H + O2 <=> O + OH
    // (0):  H + O2 <=> O + OH
    fwd_A[8]     = 26440000000000000;
    fwd_beta[8]  = -0.67069999999999996;
    fwd_Ea[8]    = 17041;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = pow(10,-12.000000);
    is_PD[8] = 0;
    nTB[8] = 0;

    // (1):  O + H2 <=> H + OH
    // (1):  O + H2 <=> H + OH
    fwd_A[9]     = 45890;
    fwd_beta[9]  = 2.7000000000000002;
    fwd_Ea[9]    = 6260;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = pow(10,-12.000000);
    is_PD[9] = 0;
    nTB[9] = 0;

    // (2):  OH + H2 <=> H + H2O
    // (2):  OH + H2 <=> H + H2O
    fwd_A[10]     = 173400000;
    fwd_beta[10]  = 1.51;
    fwd_Ea[10]    = 3430;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = pow(10,-12.000000);
    is_PD[10] = 0;
    nTB[10] = 0;

    // (3):  OH + OH <=> O + H2O
    // (3):  OH + OH <=> O + H2O
    fwd_A[11]     = 39730;
    fwd_beta[11]  = 2.3999999999999999;
    fwd_Ea[11]    = -2110;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = pow(10,-12.000000);
    is_PD[11] = 0;
    nTB[11] = 0;

    // (4):  H + H + M <=> H2 + M
    // (4):  H + H + M <=> H2 + M
    fwd_A[3]     = 1.78e+18;
    fwd_beta[3]  = -1;
    fwd_Ea[3]    = 0;
    prefactor_units[3]  = 1.0000000000000002e-12;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-12.000000);
    is_PD[3] = 0;
    nTB[3] = 5;
    TB[3] = (amrex::Real *) malloc(5 * sizeof(amrex::Real));
    TBid[3] = (int *) malloc(5 * sizeof(int));
    TBid[3][0] = 0; TB[3][0] = 0; // H2
    TBid[3][1] = 9; TB[3][1] = 0; // H2O
    TBid[3][2] = 13; TB[3][2] = 0; // CO2
    TBid[3][3] = 2; TB[3][3] = 0.63; // AR
    TBid[3][4] = 4; TB[3][4] = 0.63; // HE

    // (5):  H + H + H2 <=> H2 + H2
    // (5):  H + H + H2 <=> H2 + H2
    fwd_A[12]     = 90000000000000000;
    fwd_beta[12]  = -0.59999999999999998;
    fwd_Ea[12]    = 0;
    prefactor_units[12]  = 1.0000000000000002e-12;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = pow(10,-18.000000);
    is_PD[12] = 0;
    nTB[12] = 0;

    // (6):  H + H + H2O <=> H2 + H2O
    // (6):  H + H + H2O <=> H2 + H2O
    fwd_A[13]     = 5.624e+19;
    fwd_beta[13]  = -1.25;
    fwd_Ea[13]    = 0;
    prefactor_units[13]  = 1.0000000000000002e-12;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = pow(10,-18.000000);
    is_PD[13] = 0;
    nTB[13] = 0;

    // (7):  H + H + CO2 <=> H2 + CO2
    // (7):  H + H + CO2 <=> H2 + CO2
    fwd_A[14]     = 5.5e+20;
    fwd_beta[14]  = -2;
    fwd_Ea[14]    = 0;
    prefactor_units[14]  = 1.0000000000000002e-12;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = pow(10,-18.000000);
    is_PD[14] = 0;
    nTB[14] = 0;

    // (8):  H + OH + M <=> H2O + M
    // (8):  H + OH + M <=> H2O + M
    fwd_A[4]     = 4.4e+22;
    fwd_beta[4]  = -2;
    fwd_Ea[4]    = 0;
    prefactor_units[4]  = 1.0000000000000002e-12;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 0;
    nTB[4] = 6;
    TB[4] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[4] = (int *) malloc(6 * sizeof(int));
    TBid[4][0] = 0; TB[4][0] = 2; // H2
    TBid[4][1] = 9; TB[4][1] = 6.2999999999999998; // H2O
    TBid[4][2] = 10; TB[4][2] = 1.75; // CO
    TBid[4][3] = 13; TB[4][3] = 3.6000000000000001; // CO2
    TBid[4][4] = 2; TB[4][4] = 0.38; // AR
    TBid[4][5] = 4; TB[4][5] = 0.38; // HE

    // (9):  O + H + M <=> OH + M
    // (9):  O + H + M <=> OH + M
    fwd_A[5]     = 9.428e+18;
    fwd_beta[5]  = -1;
    fwd_Ea[5]    = 0;
    prefactor_units[5]  = 1.0000000000000002e-12;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-12.000000);
    is_PD[5] = 0;
    nTB[5] = 6;
    TB[5] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[5] = (int *) malloc(6 * sizeof(int));
    TBid[5][0] = 0; TB[5][0] = 2; // H2
    TBid[5][1] = 9; TB[5][1] = 12; // H2O
    TBid[5][2] = 10; TB[5][2] = 1.75; // CO
    TBid[5][3] = 13; TB[5][3] = 3.6000000000000001; // CO2
    TBid[5][4] = 2; TB[5][4] = 0.69999999999999996; // AR
    TBid[5][5] = 4; TB[5][5] = 0.69999999999999996; // HE

    // (10):  O + O + M <=> O2 + M
    // (10):  O + O + M <=> O2 + M
    fwd_A[6]     = 1.2e+17;
    fwd_beta[6]  = -1;
    fwd_Ea[6]    = 0;
    prefactor_units[6]  = 1.0000000000000002e-12;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = pow(10,-12.000000);
    is_PD[6] = 0;
    nTB[6] = 6;
    TB[6] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[6] = (int *) malloc(6 * sizeof(int));
    TBid[6][0] = 0; TB[6][0] = 2.3999999999999999; // H2
    TBid[6][1] = 9; TB[6][1] = 15.4; // H2O
    TBid[6][2] = 10; TB[6][2] = 1.75; // CO
    TBid[6][3] = 13; TB[6][3] = 3.6000000000000001; // CO2
    TBid[6][4] = 2; TB[6][4] = 0.82999999999999996; // AR
    TBid[6][5] = 4; TB[6][5] = 0.82999999999999996; // HE

    // (11):  H + O2 (+M) <=> HO2 (+M)
    // (11):  H + O2 (+M) <=> HO2 (+M)
    fwd_A[0]     = 5116000000000;
    fwd_beta[0]  = 0.44;
    fwd_Ea[0]    = 0;
    low_A[0]     = 6.328e+19;
    low_beta[0]  = -1.3999999999999999;
    low_Ea[0]    = 0;
    troe_a[0]    = 0.5;
    troe_Tsss[0] = 1.0000000000000001e-30;
    troe_Ts[0]   = 1e+30;
    troe_len[0]  = 3;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = pow(10,-12.000000);
    is_PD[0] = 1;
    nTB[0] = 7;
    TB[0] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[0] = (int *) malloc(7 * sizeof(int));
    TBid[0][0] = 11; TB[0][0] = 0.84999999999999998; // O2
    TBid[0][1] = 9; TB[0][1] = 11.890000000000001; // H2O
    TBid[0][2] = 10; TB[0][2] = 1.0900000000000001; // CO
    TBid[0][3] = 13; TB[0][3] = 2.1800000000000002; // CO2
    TBid[0][4] = 2; TB[0][4] = 0.40000000000000002; // AR
    TBid[0][5] = 4; TB[0][5] = 0.46000000000000002; // HE
    TBid[0][6] = 0; TB[0][6] = 0.75; // H2

    // (12):  H2 + O2 <=> HO2 + H
    // (12):  H2 + O2 <=> HO2 + H
    fwd_A[15]     = 591600;
    fwd_beta[15]  = 2.4329999999999998;
    fwd_Ea[15]    = 53502;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = pow(10,-12.000000);
    is_PD[15] = 0;
    nTB[15] = 0;

    // (13):  OH + OH (+M) <=> H2O2 (+M)
    // (13):  OH + OH (+M) <=> H2O2 (+M)
    fwd_A[1]     = 111000000000000;
    fwd_beta[1]  = -0.37;
    fwd_Ea[1]    = 0;
    low_A[1]     = 2.01e+17;
    low_beta[1]  = -0.58399999999999996;
    low_Ea[1]    = -2293;
    troe_a[1]    = 0.73460000000000003;
    troe_Tsss[1] = 94;
    troe_Ts[1]   = 1756;
    troe_Tss[1]  = 5182;
    troe_len[1]  = 4;
    prefactor_units[1]  = 1.0000000000000002e-06;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = pow(10,-12.000000);
    is_PD[1] = 1;
    nTB[1] = 6;
    TB[1] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[1] = (int *) malloc(6 * sizeof(int));
    TBid[1][0] = 0; TB[1][0] = 2; // H2
    TBid[1][1] = 9; TB[1][1] = 6; // H2O
    TBid[1][2] = 10; TB[1][2] = 1.75; // CO
    TBid[1][3] = 13; TB[1][3] = 3.6000000000000001; // CO2
    TBid[1][4] = 2; TB[1][4] = 0.69999999999999996; // AR
    TBid[1][5] = 4; TB[1][5] = 0.69999999999999996; // HE

    // (14):  HO2 + H <=> O + H2O
    // (14):  HO2 + H <=> O + H2O
    fwd_A[16]     = 3970000000000;
    fwd_beta[16]  = 0;
    fwd_Ea[16]    = 671;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = pow(10,-12.000000);
    is_PD[16] = 0;
    nTB[16] = 0;

    // (15):  HO2 + H <=> OH + OH
    // (15):  HO2 + H <=> OH + OH
    fwd_A[17]     = 74850000000000;
    fwd_beta[17]  = 0;
    fwd_Ea[17]    = 295;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = pow(10,-12.000000);
    is_PD[17] = 0;
    nTB[17] = 0;

    // (16):  HO2 + O <=> OH + O2
    // (16):  HO2 + O <=> OH + O2
    fwd_A[18]     = 40000000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 0;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = pow(10,-12.000000);
    is_PD[18] = 0;
    nTB[18] = 0;

    // (17):  HO2 + OH <=> O2 + H2O
    // (17):  HO2 + OH <=> O2 + H2O
    fwd_A[19]     = 23750000000000;
    fwd_beta[19]  = 0;
    fwd_Ea[19]    = -500;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = pow(10,-12.000000);
    is_PD[19] = 0;
    nTB[19] = 0;

    // (18):  HO2 + OH <=> O2 + H2O
    // (18):  HO2 + OH <=> O2 + H2O
    fwd_A[20]     = 10000000000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 17330;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-12.000000);
    is_PD[20] = 0;
    nTB[20] = 0;

    // (19):  HO2 + HO2 <=> O2 + H2O2
    // (19):  HO2 + HO2 <=> O2 + H2O2
    fwd_A[21]     = 130000000000;
    fwd_beta[21]  = 0;
    fwd_Ea[21]    = -1630;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = pow(10,-12.000000);
    is_PD[21] = 0;
    nTB[21] = 0;

    // (20):  HO2 + HO2 <=> O2 + H2O2
    // (20):  HO2 + HO2 <=> O2 + H2O2
    fwd_A[22]     = 365800000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 12000;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = pow(10,-12.000000);
    is_PD[22] = 0;
    nTB[22] = 0;

    // (21):  H2O2 + H <=> HO2 + H2
    // (21):  H2O2 + H <=> HO2 + H2
    fwd_A[23]     = 6050000;
    fwd_beta[23]  = 2;
    fwd_Ea[23]    = 5200;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = pow(10,-12.000000);
    is_PD[23] = 0;
    nTB[23] = 0;

    // (22):  H2O2 + H <=> OH + H2O
    // (22):  H2O2 + H <=> OH + H2O
    fwd_A[24]     = 24100000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 3970;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = pow(10,-12.000000);
    is_PD[24] = 0;
    nTB[24] = 0;

    // (23):  H2O2 + O <=> OH + HO2
    // (23):  H2O2 + O <=> OH + HO2
    fwd_A[25]     = 9630000;
    fwd_beta[25]  = 2;
    fwd_Ea[25]    = 3970;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = pow(10,-12.000000);
    is_PD[25] = 0;
    nTB[25] = 0;

    // (24):  H2O2 + OH <=> HO2 + H2O
    // (24):  H2O2 + OH <=> HO2 + H2O
    fwd_A[26]     = 2000000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 427;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = pow(10,-12.000000);
    is_PD[26] = 0;
    nTB[26] = 0;

    // (25):  H2O2 + OH <=> HO2 + H2O
    // (25):  H2O2 + OH <=> HO2 + H2O
    fwd_A[27]     = 2.6700000000000001e+41;
    fwd_beta[27]  = -7;
    fwd_Ea[27]    = 37600;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = pow(10,-12.000000);
    is_PD[27] = 0;
    nTB[27] = 0;

    // (26):  CO + O (+M) <=> CO2 (+M)
    // (26):  CO + O (+M) <=> CO2 (+M)
    fwd_A[2]     = 13620000000;
    fwd_beta[2]  = 0;
    fwd_Ea[2]    = 2384;
    low_A[2]     = 1.1729999999999999e+24;
    low_beta[2]  = -2.79;
    low_Ea[2]    = 4191;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-12.000000);
    is_PD[2] = 1;
    nTB[2] = 6;
    TB[2] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[2] = (int *) malloc(6 * sizeof(int));
    TBid[2][0] = 0; TB[2][0] = 2; // H2
    TBid[2][1] = 9; TB[2][1] = 12; // H2O
    TBid[2][2] = 10; TB[2][2] = 1.75; // CO
    TBid[2][3] = 13; TB[2][3] = 3.6000000000000001; // CO2
    TBid[2][4] = 2; TB[2][4] = 0.69999999999999996; // AR
    TBid[2][5] = 4; TB[2][5] = 0.69999999999999996; // HE

    // (27):  CO + OH <=> CO2 + H
    // (27):  CO + OH <=> CO2 + H
    fwd_A[28]     = 800000000000;
    fwd_beta[28]  = 0.14000000000000001;
    fwd_Ea[28]    = 7352;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = pow(10,-12.000000);
    is_PD[28] = 0;
    nTB[28] = 0;

    // (28):  CO + OH <=> CO2 + H
    // (28):  CO + OH <=> CO2 + H
    fwd_A[29]     = 87840000000;
    fwd_beta[29]  = 0.029999999999999999;
    fwd_Ea[29]    = -16;
    prefactor_units[29]  = 1.0000000000000002e-06;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = pow(10,-12.000000);
    is_PD[29] = 0;
    nTB[29] = 0;

    // (29):  CO + O2 <=> CO2 + O
    // (29):  CO + O2 <=> CO2 + O
    fwd_A[30]     = 1119000000000;
    fwd_beta[30]  = 0;
    fwd_Ea[30]    = 47700;
    prefactor_units[30]  = 1.0000000000000002e-06;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = pow(10,-12.000000);
    is_PD[30] = 0;
    nTB[30] = 0;

    // (30):  CO + HO2 <=> CO2 + OH
    // (30):  CO + HO2 <=> CO2 + OH
    fwd_A[31]     = 30100000000000;
    fwd_beta[31]  = 0;
    fwd_Ea[31]    = 23000;
    prefactor_units[31]  = 1.0000000000000002e-06;
    activation_units[31] = 0.50321666580471969;
    phase_units[31]      = pow(10,-12.000000);
    is_PD[31] = 0;
    nTB[31] = 0;

    // (31):  HCO + H <=> CO + H2
    // (31):  HCO + H <=> CO + H2
    fwd_A[32]     = 120000000000000;
    fwd_beta[32]  = 0;
    fwd_Ea[32]    = 0;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = pow(10,-12.000000);
    is_PD[32] = 0;
    nTB[32] = 0;

    // (32):  HCO + O <=> CO + OH
    // (32):  HCO + O <=> CO + OH
    fwd_A[33]     = 30000000000000;
    fwd_beta[33]  = 0;
    fwd_Ea[33]    = 0;
    prefactor_units[33]  = 1.0000000000000002e-06;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = pow(10,-12.000000);
    is_PD[33] = 0;
    nTB[33] = 0;

    // (33):  HCO + O <=> CO2 + H
    // (33):  HCO + O <=> CO2 + H
    fwd_A[34]     = 30000000000000;
    fwd_beta[34]  = 0;
    fwd_Ea[34]    = 0;
    prefactor_units[34]  = 1.0000000000000002e-06;
    activation_units[34] = 0.50321666580471969;
    phase_units[34]      = pow(10,-12.000000);
    is_PD[34] = 0;
    nTB[34] = 0;

    // (34):  HCO + OH <=> CO + H2O
    // (34):  HCO + OH <=> CO + H2O
    fwd_A[35]     = 30200000000000;
    fwd_beta[35]  = 0;
    fwd_Ea[35]    = 0;
    prefactor_units[35]  = 1.0000000000000002e-06;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = pow(10,-12.000000);
    is_PD[35] = 0;
    nTB[35] = 0;

    // (35):  HCO + M <=> CO + H + M
    // (35):  HCO + M <=> CO + H + M
    fwd_A[7]     = 1.87e+17;
    fwd_beta[7]  = -1;
    fwd_Ea[7]    = 17000;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = pow(10,-6.000000);
    is_PD[7] = 0;
    nTB[7] = 4;
    TB[7] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[7] = (int *) malloc(4 * sizeof(int));
    TBid[7][0] = 0; TB[7][0] = 2; // H2
    TBid[7][1] = 9; TB[7][1] = 0; // H2O
    TBid[7][2] = 10; TB[7][2] = 1.75; // CO
    TBid[7][3] = 13; TB[7][3] = 3.6000000000000001; // CO2

    // (36):  HCO + H2O <=> CO + H + H2O
    // (36):  HCO + H2O <=> CO + H + H2O
    fwd_A[36]     = 2.244e+18;
    fwd_beta[36]  = -1;
    fwd_Ea[36]    = 17000;
    prefactor_units[36]  = 1.0000000000000002e-06;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = pow(10,-12.000000);
    is_PD[36] = 0;
    nTB[36] = 0;

    // (37):  HCO + O2 <=> CO + HO2
    // (37):  HCO + O2 <=> CO + HO2
    fwd_A[37]     = 12040000000;
    fwd_beta[37]  = 0.80700000000000005;
    fwd_Ea[37]    = -727;
    prefactor_units[37]  = 1.0000000000000002e-06;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = pow(10,-12.000000);
    is_PD[37] = 0;
    nTB[37] = 0;

}


/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<38; ++i) {
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
    awt[4] = 39.948000; /*AR */
    awt[5] = 4.002600; /*HE */

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
    int kd = 6; 
    /*Zero ncf */
    for (id = 0; id < kd * 14; ++ id) {
         ncf[id] = 0; 
    }

    /*H2 */
    ncf[ 0 * kd + 1 ] = 2; /*H */

    /*H */
    ncf[ 1 * kd + 1 ] = 1; /*H */

    /*AR */
    ncf[ 2 * kd + 4 ] = 1; /*AR */

    /*N2 */
    ncf[ 3 * kd + 3 ] = 2; /*N */

    /*HE */
    ncf[ 4 * kd + 5 ] = 1; /*HE */

    /*O */
    ncf[ 5 * kd + 0 ] = 1; /*O */

    /*OH */
    ncf[ 6 * kd + 0 ] = 1; /*O */
    ncf[ 6 * kd + 1 ] = 1; /*H */

    /*HCO */
    ncf[ 7 * kd + 1 ] = 1; /*H */
    ncf[ 7 * kd + 2 ] = 1; /*C */
    ncf[ 7 * kd + 0 ] = 1; /*O */

    /*HO2 */
    ncf[ 8 * kd + 1 ] = 1; /*H */
    ncf[ 8 * kd + 0 ] = 2; /*O */

    /*H2O */
    ncf[ 9 * kd + 1 ] = 2; /*H */
    ncf[ 9 * kd + 0 ] = 1; /*O */

    /*CO */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 0 ] = 1; /*O */

    /*O2 */
    ncf[ 11 * kd + 0 ] = 2; /*O */

    /*H2O2 */
    ncf[ 12 * kd + 1 ] = 2; /*H */
    ncf[ 12 * kd + 0 ] = 2; /*O */

    /*CO2 */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 0 ] = 2; /*O */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(6);
    ename[0] = "O";
    ename[1] = "H";
    ename[2] = "C";
    ename[3] = "N";
    ename[4] = "AR";
    ename[5] = "HE";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(14);
    kname[0] = "H2";
    kname[1] = "H";
    kname[2] = "AR";
    kname[3] = "N2";
    kname[4] = "HE";
    kname[5] = "O";
    kname[6] = "OH";
    kname[7] = "HCO";
    kname[8] = "HO2";
    kname[9] = "H2O";
    kname[10] = "CO";
    kname[11] = "O2";
    kname[12] = "H2O2";
    kname[13] = "CO2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(14);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<14; l++) {
                c_d[l] = 1.0/ 14.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<15; k++) {
        for (int l=0; l<15; l++) {
            if(J_h[ 15 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(14);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<14; k++) {
                c_d[k] = 1.0/ 14.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<15; k++) {
        for (int l=0; l<15; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 15 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(14);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<14; k++) {
                c_d[k] = 1.0/ 14.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<15; k++) {
        for (int l=0; l<15; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 15 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(14);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<14; k++) {
                c_d[k] = 1.0/ 14.000000 ;
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
        offset_row = nc * 15;
        offset_col = nc * 15;
        for (int k=0; k<15; k++) {
            for (int l=0; l<15; l++) {
                if(J_h[15*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(14);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<14; k++) {
                c_d[k] = 1.0/ 14.000000 ;
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
            offset = nc * 15;
            for (int l=0; l<15; l++) {
                for (int k=0; k<15; k++) {
                    if(J_h[15*k + l] != 0.0) {
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
            offset = nc * 15;
            for (int l=0; l<15; l++) {
                for (int k=0; k<15; k++) {
                    if(J_h[15*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(14);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<14; k++) {
                c_d[k] = 1.0/ 14.000000 ;
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
            offset = nc * 15;
            for (int l=0; l<15; l++) {
                for (int k=0; k<15; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[15*k + l] != 0.0) {
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
            offset = nc * 15;
            for (int l=0; l<15; l++) {
                for (int k=0; k<15; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[15*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(14);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<14; k++) {
                c_d[k] = 1.0/ 14.000000 ;
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
    for (int k=0; k<15; k++) {
        for (int l=0; l<15; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 15*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[15*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 15*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(14);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<14; k++) {
                c_d[k] = 1.0/ 14.000000 ;
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
        for (int l=0; l<15; l++) {
            for (int k=0; k<15; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[15*k + l] != 0.0) {
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
        for (int l=0; l<15; l++) {
            for (int k=0; k<15; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[15*k + l] != 0.0) {
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

