#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/* PURE CPU stuff  */
#ifndef AMREX_USE_GPU
namespace thermo
{
    amrex::Real fwd_A[27], fwd_beta[27], fwd_Ea[27];
    amrex::Real low_A[27], low_beta[27], low_Ea[27];
    amrex::Real rev_A[27], rev_beta[27], rev_Ea[27];
    amrex::Real troe_a[27],troe_Ts[27], troe_Tss[27], troe_Tsss[27];
    amrex::Real sri_a[27], sri_b[27], sri_c[27], sri_d[27], sri_e[27];
    amrex::Real activation_units[27], prefactor_units[27], phase_units[27];
    int is_PD[27], troe_len[27], sri_len[27], nTB[27], *TBid[27];
    amrex::Real *TB[27];
    std::vector<std::vector<amrex::Real>> kiv(27); 
    std::vector<std::vector<amrex::Real>> nuv(27); 
};

using namespace thermo;


/*Returns a count of species in a reaction, and their indices */
/*and stoichiometric coefficients. (Eq 50) */
void CKINU(int * i, int * nspec, int * ki, int * nu)
{
    if (*i < 1) {
        /*Return max num species per reaction */
        *nspec = 5;
    } else {
        if (*i > 27) {
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
    amrex::Real imw[9];

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
void VCKHMS(int *  np, amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tc[5], h[9];
    amrex::Real imw[9];

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


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int *  np, amrex::Real *  rho, amrex::Real *  T,
	    amrex::Real *  y,
	    amrex::Real *  wdot)
{
    amrex::Real c[9*(*np)]; /*temporary storage */
    amrex::Real imw[9];

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


/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW[*np];
    amrex::Real imw[9];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

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
void vproductionRate(int npt, amrex::Real *  wdot, amrex::Real *  sc, amrex::Real *  T)
{
    amrex::Real k_f_s[27*npt], Kc_s[27*npt], mixture[npt], g_RT[9*npt];
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
    }
}

void vcomp_gibbs(int npt, amrex::Real *  g_RT, amrex::Real *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        amrex::Real tg[5], g[9];
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

void vcomp_Kc(int npt, amrex::Real *  Kc_s, amrex::Real *  g_RT, amrex::Real *  invT)
{
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        amrex::Real refC = (101325. / 8.31451) * invT[i];
        amrex::Real refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((2.000000 * g_RT[4*npt+i]) - (g_RT[7*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((2.000000 * g_RT[2*npt+i]) - (g_RT[3*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[2*npt+i]) - (g_RT[4*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i]) - (g_RT[0*npt+i]));
        Kc_s[5*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[6*npt+i] = exp((g_RT[0*npt+i] + g_RT[2*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i]));
        Kc_s[7*npt+i] = exp((g_RT[2*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[8*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[9*npt+i] = refCinv * exp((g_RT[1*npt+i] + 2.000000 * g_RT[3*npt+i]) - (g_RT[3*npt+i] + g_RT[6*npt+i]));
        Kc_s[10*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[11*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[8*npt+i]) - (g_RT[6*npt+i] + g_RT[8*npt+i]));
        Kc_s[12*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[2*npt+i] + g_RT[4*npt+i]));
        Kc_s[13*npt+i] = refCinv * exp((g_RT[0*npt+i] + 2.000000 * g_RT[1*npt+i]) - (2.000000 * g_RT[0*npt+i]));
        Kc_s[14*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[0*npt+i] + g_RT[5*npt+i]));
        Kc_s[15*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (2.000000 * g_RT[4*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[0*npt+i] + g_RT[6*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[0*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i]));
        Kc_s[21*npt+i] = exp((2.000000 * g_RT[4*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[4*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[4*npt+i] + g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[4*npt+i] + g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[25*npt+i] = exp((2.000000 * g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[7*npt+i]));
        Kc_s[26*npt+i] = exp((2.000000 * g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[7*npt+i]));
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

        /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
        phi_f = pow(sc[4*npt+i], 2.000000);
        alpha = mixture[i] + (TB[0][0] - 1)*sc[0*npt+i] + (TB[0][1] - 1)*sc[5*npt+i];
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
        phi_r = sc[7*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= 2.000000 * qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 2: 2.000000 O + M <=> O2 + M */
        phi_f = pow(sc[2*npt+i], 2.000000);
        alpha = mixture[i] + (TB[1][0] - 1)*sc[0*npt+i] + (TB[1][1] - 1)*sc[5*npt+i];
        k_f = alpha * k_f_s[1*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= 2.000000 * qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 3: O + H + M <=> OH + M */
        phi_f = sc[1*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[0*npt+i] + (TB[2][1] - 1)*sc[5*npt+i];
        k_f = alpha * k_f_s[2*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 4: H + O2 + M <=> HO2 + M */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[3*npt+i] + (TB[3][1] - 1)*sc[5*npt+i] + (TB[3][2] - 1)*sc[8*npt+i];
        k_f = alpha * k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 5: 2.000000 H + M <=> H2 + M */
        phi_f = pow(sc[1*npt+i], 2.000000);
        alpha = mixture[i] + (TB[4][0] - 1)*sc[0*npt+i] + (TB[4][1] - 1)*sc[5*npt+i];
        k_f = alpha * k_f_s[4*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2.000000 * qdot;

        /*reaction 6: H + OH + M <=> H2O + M */
        phi_f = sc[1*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[0*npt+i] + (TB[5][1] - 1)*sc[5*npt+i];
        k_f = alpha * k_f_s[5*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 7: O + H2 <=> H + OH */
        phi_f = sc[0*npt+i]*sc[2*npt+i];
        k_f = k_f_s[6*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 8: O + HO2 <=> OH + O2 */
        phi_f = sc[2*npt+i]*sc[6*npt+i];
        k_f = k_f_s[7*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[4*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 9: O + H2O2 <=> OH + HO2 */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 10: H + 2.000000 O2 <=> HO2 + O2 */
        phi_f = sc[1*npt+i]*pow(sc[3*npt+i], 2.000000);
        k_f = k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[6*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= 2.000000 * qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[5*npt+i];
        k_f = k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[8*npt+i];
        k_f = k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[8*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 13: H + O2 <=> O + OH */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[4*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 14: 2.000000 H + H2 <=> 2.000000 H2 */
        phi_f = sc[0*npt+i]*pow(sc[1*npt+i], 2.000000);
        k_f = k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[0*npt+i], 2.000000);
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] += 2.000000 * qdot;
        wdot[1*npt+i] -= 2.000000 * qdot;

        /*reaction 15: 2.000000 H + H2O <=> H2 + H2O */
        phi_f = pow(sc[1*npt+i], 2.000000)*sc[5*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[5*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2.000000 * qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 16: H + HO2 <=> O + H2O */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 17: H + HO2 <=> O2 + H2 */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 18: H + HO2 <=> 2.000000 OH */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[4*npt+i], 2.000000);
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += 2.000000 * qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 19: H + H2O2 <=> HO2 + H2 */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[6*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 20: H + H2O2 <=> OH + H2O */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[5*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 21: OH + H2 <=> H + H2O */
        phi_f = sc[0*npt+i]*sc[4*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 22: 2.000000 OH <=> O + H2O */
        phi_f = pow(sc[4*npt+i], 2.000000);
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= 2.000000 * qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 23: OH + HO2 <=> O2 + H2O */
        phi_f = sc[4*npt+i]*sc[6*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 24: OH + H2O2 <=> HO2 + H2O */
        phi_f = sc[4*npt+i]*sc[7*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 25: OH + H2O2 <=> HO2 + H2O */
        phi_f = sc[4*npt+i]*sc[7*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
        phi_f = pow(sc[6*npt+i], 2.000000);
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[7*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= 2.000000 * qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
        phi_f = pow(sc[6*npt+i], 2.000000);
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[7*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= 2.000000 * qdot;
        wdot[7*npt+i] += qdot;
    }
}

static amrex::Real T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static amrex::Real k_f_save[27];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static amrex::Real Kc_save[27];
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

    amrex::Real qdot, q_f[27], q_r[27];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 9; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[4] -= 2.000000 * qdot;
    wdot[7] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[2] -= 2.000000 * qdot;
    wdot[3] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;

    qdot = q_f[5]-q_r[5];
    wdot[1] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[8]-q_r[8];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[9]-q_r[9];
    wdot[1] -= qdot;
    wdot[3] -= 2.000000 * qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[0] -= qdot;
    wdot[0] += 2.000000 * qdot;
    wdot[1] -= 2.000000 * qdot;

    qdot = q_f[14]-q_r[14];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[16]-q_r[16];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] -= qdot;
    wdot[4] += 2.000000 * qdot;
    wdot[6] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[2] += qdot;
    wdot[4] -= 2.000000 * qdot;
    wdot[5] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[3] += qdot;
    wdot[6] -= 2.000000 * qdot;
    wdot[7] += qdot;

    qdot = q_f[26]-q_r[26];
    wdot[3] += qdot;
    wdot[6] -= 2.000000 * qdot;
    wdot[7] += qdot;

    return;
}

void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f)
{
    for (int i=0; i<27; ++i) {
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

    Kc[0] = 2.000000*g_RT[4] - g_RT[7];
    Kc[1] = 2.000000*g_RT[2] - g_RT[3];
    Kc[2] = g_RT[1] + g_RT[2] - g_RT[4];
    Kc[3] = g_RT[1] + g_RT[3] - g_RT[6];
    Kc[4] = -g_RT[0] + 2.000000*g_RT[1];
    Kc[5] = g_RT[1] + g_RT[4] - g_RT[5];
    Kc[6] = g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4];
    Kc[7] = g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6];
    Kc[8] = g_RT[2] - g_RT[4] - g_RT[6] + g_RT[7];
    Kc[9] = g_RT[1] + 2.000000*g_RT[3] - g_RT[3] - g_RT[6];
    Kc[10] = g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6];
    Kc[11] = g_RT[1] + g_RT[3] - g_RT[6] + g_RT[8] - g_RT[8];
    Kc[12] = g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4];
    Kc[13] = g_RT[0] - 2.000000*g_RT[0] + 2.000000*g_RT[1];
    Kc[14] = -g_RT[0] + 2.000000*g_RT[1] + g_RT[5] - g_RT[5];
    Kc[15] = g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6];
    Kc[16] = -g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6];
    Kc[17] = g_RT[1] - 2.000000*g_RT[4] + g_RT[6];
    Kc[18] = -g_RT[0] + g_RT[1] - g_RT[6] + g_RT[7];
    Kc[19] = g_RT[1] - g_RT[4] - g_RT[5] + g_RT[7];
    Kc[20] = g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5];
    Kc[21] = -g_RT[2] + 2.000000*g_RT[4] - g_RT[5];
    Kc[22] = -g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6];
    Kc[23] = g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7];
    Kc[24] = g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7];
    Kc[25] = -g_RT[3] + 2.000000*g_RT[6] - g_RT[7];
    Kc[26] = -g_RT[3] + 2.000000*g_RT[6] - g_RT[7];

    for (int i=0; i<27; ++i) {
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
    Kc[9] *= refCinv;
    Kc[10] *= refCinv;
    Kc[11] *= refCinv;
    Kc[13] *= refCinv;
    Kc[14] *= refCinv;

    return;
}

void comp_qfqr(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT)
{

    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    qf[0] = pow(sc[4], 2.000000);
    qr[0] = sc[7];

    /*reaction 2: 2.000000 O + M <=> O2 + M */
    qf[1] = pow(sc[2], 2.000000);
    qr[1] = sc[3];

    /*reaction 3: O + H + M <=> OH + M */
    qf[2] = sc[1]*sc[2];
    qr[2] = sc[4];

    /*reaction 4: H + O2 + M <=> HO2 + M */
    qf[3] = sc[1]*sc[3];
    qr[3] = sc[6];

    /*reaction 5: 2.000000 H + M <=> H2 + M */
    qf[4] = pow(sc[1], 2.000000);
    qr[4] = sc[0];

    /*reaction 6: H + OH + M <=> H2O + M */
    qf[5] = sc[1]*sc[4];
    qr[5] = sc[5];

    /*reaction 7: O + H2 <=> H + OH */
    qf[6] = sc[0]*sc[2];
    qr[6] = sc[1]*sc[4];

    /*reaction 8: O + HO2 <=> OH + O2 */
    qf[7] = sc[2]*sc[6];
    qr[7] = sc[3]*sc[4];

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    qf[8] = sc[2]*sc[7];
    qr[8] = sc[4]*sc[6];

    /*reaction 10: H + 2.000000 O2 <=> HO2 + O2 */
    qf[9] = sc[1]*pow(sc[3], 2.000000);
    qr[9] = sc[3]*sc[6];

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    qf[10] = sc[1]*sc[3]*sc[5];
    qr[10] = sc[5]*sc[6];

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    qf[11] = sc[1]*sc[3]*sc[8];
    qr[11] = sc[6]*sc[8];

    /*reaction 13: H + O2 <=> O + OH */
    qf[12] = sc[1]*sc[3];
    qr[12] = sc[2]*sc[4];

    /*reaction 14: 2.000000 H + H2 <=> 2.000000 H2 */
    qf[13] = sc[0]*pow(sc[1], 2.000000);
    qr[13] = pow(sc[0], 2.000000);

    /*reaction 15: 2.000000 H + H2O <=> H2 + H2O */
    qf[14] = pow(sc[1], 2.000000)*sc[5];
    qr[14] = sc[0]*sc[5];

    /*reaction 16: H + HO2 <=> O + H2O */
    qf[15] = sc[1]*sc[6];
    qr[15] = sc[2]*sc[5];

    /*reaction 17: H + HO2 <=> O2 + H2 */
    qf[16] = sc[1]*sc[6];
    qr[16] = sc[0]*sc[3];

    /*reaction 18: H + HO2 <=> 2.000000 OH */
    qf[17] = sc[1]*sc[6];
    qr[17] = pow(sc[4], 2.000000);

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    qf[18] = sc[1]*sc[7];
    qr[18] = sc[0]*sc[6];

    /*reaction 20: H + H2O2 <=> OH + H2O */
    qf[19] = sc[1]*sc[7];
    qr[19] = sc[4]*sc[5];

    /*reaction 21: OH + H2 <=> H + H2O */
    qf[20] = sc[0]*sc[4];
    qr[20] = sc[1]*sc[5];

    /*reaction 22: 2.000000 OH <=> O + H2O */
    qf[21] = pow(sc[4], 2.000000);
    qr[21] = sc[2]*sc[5];

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    qf[22] = sc[4]*sc[6];
    qr[22] = sc[3]*sc[5];

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    qf[23] = sc[4]*sc[7];
    qr[23] = sc[5]*sc[6];

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    qf[24] = sc[4]*sc[7];
    qr[24] = sc[5]*sc[6];

    /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
    qf[25] = pow(sc[6], 2.000000);
    qr[25] = sc[3]*sc[7];

    /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
    qf[26] = pow(sc[6], 2.000000);
    qr[26] = sc[3]*sc[7];

    amrex::Real T = tc[1];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 9; ++i) {
        mixture += sc[i];
    }

    amrex::Real Corr[27];
    for (int i = 0; i < 27; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        amrex::Real alpha[1];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[5];
        for (int i=0; i<1; i++)
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
        alpha = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[5];
        Corr[1] = alpha;
        alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[5];
        Corr[2] = alpha;
        alpha = mixture + (TB[3][0] - 1)*sc[3] + (TB[3][1] - 1)*sc[5] + (TB[3][2] - 1)*sc[8];
        Corr[3] = alpha;
        alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[5];
        Corr[4] = alpha;
        alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[5];
        Corr[5] = alpha;
    }

    for (int i=0; i<27; i++)
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

    amrex::Real q_f[27], q_r[27];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 27; ++i) {
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
    amrex::Real c[9]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 27; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
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

    for (id = 0; id < 27; ++id) {
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
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
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
    for (id = 0; id < 27; ++id) {
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
    for (id = 0; id < 27; ++id) {
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
    for (id = 0; id < 27; ++id) {
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
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
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
    for (id = 0; id < 27; ++id) {
        qdot[id] *= 1.0e-6;
    }
}

/*compute the reaction Jacobian on CPU */
void aJacobian(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP)
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
    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[4], 2.000000);
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
    phi_r = sc[7];
    Kc = refCinv * exp(2.000000*g_RT[4] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[7]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[4] -= 2 * q; /* OH */
    wdot[7] += q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[4] += -2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[7] += dqdci;                /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  + k_f*2.000000*sc[4];
        J[44] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
        J[47] += dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[54] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
        J[57] += dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[H2O2] */
        dqdci =  - k_r;
        J[74] += -2 * dqdci;          /* dwdot[OH]/d[H2O2] */
        J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    }
    else {
        dqdc[0] = TB[0][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac + k_f*2.000000*sc[4];
        dqdc[5] = TB[0][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac - k_r;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+4] += -2 * dqdc[k];
            J[10*k+7] += dqdc[k];
        }
    }
    J[94] += -2 * dqdT; /* dwdot[OH]/dT */
    J[97] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 2: 2.000000 O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[2], 2.000000);
    k_f = prefactor_units[1] * fwd_A[1]
                * exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
    dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(2.000000*g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2]) + (h_RT[3]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[3] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[1][0] - 1)*q_nocor;
        J[2] += -2 * dqdci;           /* dwdot[O]/d[H2] */
        J[3] += dqdci;                /* dwdot[O2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[2];
        J[22] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[23] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[32] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        J[33] += dqdci;               /* dwdot[O2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[1][1] - 1)*q_nocor;
        J[52] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
    }
    else {
        dqdc[0] = TB[1][0]*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*2.000000*sc[2];
        dqdc[3] = q_nocor - k_r;
        dqdc[4] = q_nocor;
        dqdc[5] = TB[1][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+2] += -2 * dqdc[k];
            J[10*k+3] += dqdc[k];
        }
    }
    J[92] += -2 * dqdT; /* dwdot[O]/dT */
    J[93] += dqdT; /* dwdot[O2]/dT */

    /*reaction 3: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[5];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[4]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[2][0] - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[4] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[11] -= dqdci;               /* dwdot[H]/d[H] */
        J[12] -= dqdci;               /* dwdot[O]/d[H] */
        J[14] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[1];
        J[21] -= dqdci;               /* dwdot[H]/d[O] */
        J[22] -= dqdci;               /* dwdot[O]/d[O] */
        J[24] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[41] -= dqdci;               /* dwdot[H]/d[OH] */
        J[42] -= dqdci;               /* dwdot[O]/d[OH] */
        J[44] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*q_nocor;
        J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[52] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[54] += dqdci;               /* dwdot[OH]/d[H2O] */
    }
    else {
        dqdc[0] = TB[2][0]*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[2];
        dqdc[2] = q_nocor + k_f*sc[1];
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor - k_r;
        dqdc[5] = TB[2][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+2] -= dqdc[k];
            J[10*k+4] += dqdc[k];
        }
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[92] -= dqdT; /* dwdot[O]/dT */
    J[94] += dqdT; /* dwdot[OH]/dT */

    /*reaction 4: H + O2 + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[3] + (TB[3][1] - 1)*sc[5] + (TB[3][2] - 1)*sc[8];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[3];
        J[11] -= dqdci;               /* dwdot[H]/d[H] */
        J[13] -= dqdci;               /* dwdot[O2]/d[H] */
        J[16] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[O2] */
        dqdci = (TB[3][0] - 1)*q_nocor + k_f*sc[1];
        J[31] -= dqdci;               /* dwdot[H]/d[O2] */
        J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*q_nocor;
        J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[53] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
        /* d()/d[N2] */
        dqdci = (TB[3][2] - 1)*q_nocor;
        J[81] -= dqdci;               /* dwdot[H]/d[N2] */
        J[83] -= dqdci;               /* dwdot[O2]/d[N2] */
        J[86] += dqdci;               /* dwdot[HO2]/d[N2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*sc[3];
        dqdc[2] = q_nocor;
        dqdc[3] = TB[3][0]*q_nocor + k_f*sc[1];
        dqdc[4] = q_nocor;
        dqdc[5] = TB[3][1]*q_nocor;
        dqdc[6] = q_nocor - k_r;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[3][2]*q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+3] -= dqdc[k];
            J[10*k+6] += dqdc[k];
        }
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[93] -= dqdT; /* dwdot[O2]/dT */
    J[96] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 5: 2.000000 H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[0];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1]);
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
        dqdci = (TB[4][0] - 1)*q_nocor - k_r;
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*2.000000*sc[1];
        J[10] += dqdci;               /* dwdot[H2]/d[H] */
        J[11] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*q_nocor;
        J[50] += dqdci;               /* dwdot[H2]/d[H2O] */
        J[51] += -2 * dqdci;          /* dwdot[H]/d[H2O] */
    }
    else {
        dqdc[0] = TB[4][0]*q_nocor - k_r;
        dqdc[1] = q_nocor + k_f*2.000000*sc[1];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = TB[4][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+0] += dqdc[k];
            J[10*k+1] += -2 * dqdc[k];
        }
    }
    J[90] += dqdT; /* dwdot[H2]/dT */
    J[91] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 6: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[5];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[4]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[5][0] - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[11] -= dqdci;               /* dwdot[H]/d[H] */
        J[14] -= dqdci;               /* dwdot[OH]/d[H] */
        J[15] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[1];
        J[41] -= dqdci;               /* dwdot[H]/d[OH] */
        J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*q_nocor - k_r;
        J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
        J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    }
    else {
        dqdc[0] = TB[5][0]*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[4];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f*sc[1];
        dqdc[5] = TB[5][1]*q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+4] -= dqdc[k];
            J[10*k+5] += dqdc[k];
        }
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[94] -= dqdT; /* dwdot[OH]/dT */
    J[95] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 7: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[1] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[2] -= dqdci;                /* dwdot[O]/d[H2] */
    J[4] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[10] -= dqdci;               /* dwdot[H2]/d[H] */
    J[11] += dqdci;               /* dwdot[H]/d[H] */
    J[12] -= dqdci;               /* dwdot[O]/d[H] */
    J[14] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[20] -= dqdci;               /* dwdot[H2]/d[O] */
    J[21] += dqdci;               /* dwdot[H]/d[O] */
    J[22] -= dqdci;               /* dwdot[O]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[40] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[41] += dqdci;               /* dwdot[H]/d[OH] */
    J[42] -= dqdci;               /* dwdot[O]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[91] += dqdT;                /* dwdot[H]/dT */
    J[92] -= dqdT;                /* dwdot[O]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[22] -= dqdci;               /* dwdot[O]/d[O] */
    J[23] += dqdci;               /* dwdot[O2]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    J[26] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[32] -= dqdci;               /* dwdot[O]/d[O2] */
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[34] += dqdci;               /* dwdot[OH]/d[O2] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[42] -= dqdci;               /* dwdot[O]/d[OH] */
    J[43] += dqdci;               /* dwdot[O2]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[62] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[64] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[92] -= dqdT;                /* dwdot[O]/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[22] -= dqdci;               /* dwdot[O]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    J[26] += dqdci;               /* dwdot[HO2]/d[O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[42] -= dqdci;               /* dwdot[O]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[62] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[64] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[72] -= dqdci;               /* dwdot[O]/d[H2O2] */
    J[74] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] -= dqdT;                /* dwdot[O]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 10: H + 2.000000 O2 <=> HO2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*pow(sc[3], 2.000000);
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[6];
    Kc = refCinv * exp(g_RT[1] + 2.000000*g_RT[3] - g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + 2.000000*h_RT[3]) + (h_RT[3] + h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*pow(sc[3], 2.000000);
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*2.000000*sc[3] - k_r*sc[6];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[5]) + (h_RT[5] + h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[5];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[5];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[53] -= dqdci;               /* dwdot[O2]/d[H2O] */
    J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[8];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[8];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[8] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[8]) + (h_RT[6] + h_RT[8]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[8];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[8];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[8];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[81] -= dqdci;               /* dwdot[H]/d[N2] */
    J[83] -= dqdci;               /* dwdot[O2]/d[N2] */
    J[86] += dqdci;               /* dwdot[HO2]/d[N2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */

    /*reaction 13: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[12] += dqdci;               /* dwdot[O]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[14] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[21] -= dqdci;               /* dwdot[H]/d[O] */
    J[22] += dqdci;               /* dwdot[O]/d[O] */
    J[23] -= dqdci;               /* dwdot[O2]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[32] += dqdci;               /* dwdot[O]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[34] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[41] -= dqdci;               /* dwdot[H]/d[OH] */
    J[42] += dqdci;               /* dwdot[O]/d[OH] */
    J[43] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[92] += dqdT;                /* dwdot[O]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 14: 2.000000 H + H2 <=> 2.000000 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*pow(sc[1], 2.000000);
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000);
    Kc = refCinv * exp(g_RT[0] - 2.000000*g_RT[0] + 2.000000*g_RT[1]);
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
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] += -2 * dqdT;           /* dwdot[H]/dT */

    /*reaction 15: 2.000000 H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[5];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[5];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[5]) + (h_RT[0] + h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[5];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[50] += dqdci;               /* dwdot[H2]/d[H2O] */
    J[51] += -2 * dqdci;          /* dwdot[H]/d[H2O] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] += -2 * dqdT;           /* dwdot[H]/dT */

    /*reaction 16: H + HO2 <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[5] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[12] += dqdci;               /* dwdot[O]/d[H] */
    J[15] += dqdci;               /* dwdot[H2O]/d[H] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[21] -= dqdci;               /* dwdot[H]/d[O] */
    J[22] += dqdci;               /* dwdot[O]/d[O] */
    J[25] += dqdci;               /* dwdot[H2O]/d[O] */
    J[26] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[52] += dqdci;               /* dwdot[O]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[62] += dqdci;               /* dwdot[O]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[92] += dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[3] += dqdci;                /* dwdot[O2]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] += dqdci;               /* dwdot[O2]/d[H] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[30] += dqdci;               /* dwdot[H2]/d[O2] */
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 18: H + HO2 <=> 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = pow(sc[4], 2.000000);
    Kc = exp(g_RT[1] - 2.000000*g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (2.000000*h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[14] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[4];
    J[41] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[64] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[94] += 2 * dqdT;            /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    J[17] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[70] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[71] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[14] += dqdci;               /* dwdot[OH]/d[H] */
    J[15] += dqdci;               /* dwdot[H2O]/d[H] */
    J[17] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[41] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[54] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[71] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[74] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
    J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[10] -= dqdci;               /* dwdot[H2]/d[H] */
    J[11] += dqdci;               /* dwdot[H]/d[H] */
    J[14] -= dqdci;               /* dwdot[OH]/d[H] */
    J[15] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[40] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[41] += dqdci;               /* dwdot[H]/d[OH] */
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[50] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[51] += dqdci;               /* dwdot[H]/d[H2O] */
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[91] += dqdT;                /* dwdot[H]/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */

    /*reaction 22: 2.000000 OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + 2.000000*g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[4] -= 2 * q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[22] += dqdci;               /* dwdot[O]/d[O] */
    J[24] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[25] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[4];
    J[42] += dqdci;               /* dwdot[O]/d[OH] */
    J[44] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[52] += dqdci;               /* dwdot[O]/d[H2O] */
    J[54] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[O]/dT */
    J[94] += -2 * dqdT;           /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[5];
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[34] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[35] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[43] += dqdci;               /* dwdot[O2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[64] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
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
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[46] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
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
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[46] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[7];
    Kc = exp(-g_RT[3] + 2.000000*g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[3] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[37] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[3];
    J[73] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[7];
    Kc = exp(-g_RT[3] + 2.000000*g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[3] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[37] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[3];
    J[73] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

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
/* Transport function declarations  */


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
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 2.01594000E+00;
    WT[1] = 1.00797000E+00;
    WT[2] = 1.59994000E+01;
    WT[3] = 3.19988000E+01;
    WT[4] = 1.70073700E+01;
    WT[5] = 1.80153400E+01;
    WT[6] = 3.30067700E+01;
    WT[7] = 3.40147400E+01;
    WT[8] = 2.80134000E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 3.80000000E+01;
    EPS[1] = 1.45000000E+02;
    EPS[2] = 8.00000000E+01;
    EPS[3] = 1.07400000E+02;
    EPS[4] = 8.00000000E+01;
    EPS[5] = 5.72400000E+02;
    EPS[6] = 1.07400000E+02;
    EPS[7] = 1.07400000E+02;
    EPS[8] = 9.75300000E+01;
}


/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 2.92000000E+00;
    SIG[1] = 2.05000000E+00;
    SIG[2] = 2.75000000E+00;
    SIG[3] = 3.45800000E+00;
    SIG[4] = 2.75000000E+00;
    SIG[5] = 2.60500000E+00;
    SIG[6] = 3.45800000E+00;
    SIG[7] = 3.45800000E+00;
    SIG[8] = 3.62100000E+00;
}


/*the dipole moment in Debye */
void egtransetDIP(amrex::Real* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[2] = 0.00000000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 1.84400000E+00;
    DIP[6] = 0.00000000E+00;
    DIP[7] = 0.00000000E+00;
    DIP[8] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 7.90000000E-01;
    POL[1] = 0.00000000E+00;
    POL[2] = 0.00000000E+00;
    POL[3] = 1.60000000E+00;
    POL[4] = 0.00000000E+00;
    POL[5] = 0.00000000E+00;
    POL[6] = 0.00000000E+00;
    POL[7] = 0.00000000E+00;
    POL[8] = 1.76000000E+00;
}


/*the rotational relaxation collision number at 298 K */
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 2.80000000E+02;
    ZROT[1] = 0.00000000E+00;
    ZROT[2] = 0.00000000E+00;
    ZROT[3] = 3.80000000E+00;
    ZROT[4] = 0.00000000E+00;
    ZROT[5] = 4.00000000E+00;
    ZROT[6] = 1.00000000E+00;
    ZROT[7] = 3.80000000E+00;
    ZROT[8] = 4.00000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 1;
    NLIN[1] = 0;
    NLIN[2] = 0;
    NLIN[3] = 1;
    NLIN[4] = 1;
    NLIN[5] = 2;
    NLIN[6] = 2;
    NLIN[7] = 2;
    NLIN[8] = 1;
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
    COFETA[8] = -1.48001581E+01;
    COFETA[9] = 1.79491990E+00;
    COFETA[10] = -1.54008440E-01;
    COFETA[11] = 6.86719439E-03;
    COFETA[12] = -1.68118868E+01;
    COFETA[13] = 2.52362554E+00;
    COFETA[14] = -2.49309128E-01;
    COFETA[15] = 1.10211025E-02;
    COFETA[16] = -1.47696103E+01;
    COFETA[17] = 1.79491990E+00;
    COFETA[18] = -1.54008440E-01;
    COFETA[19] = 6.86719439E-03;
    COFETA[20] = -1.17770937E+01;
    COFETA[21] = -8.26742721E-01;
    COFETA[22] = 3.39009079E-01;
    COFETA[23] = -2.00674327E-02;
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
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = 1.15899058E+01;
    COFLAM[1] = -1.52427727E+00;
    COFLAM[2] = 2.72840752E-01;
    COFLAM[3] = -1.03392618E-02;
    COFLAM[4] = -3.24539191E-01;
    COFLAM[5] = 3.41660514E+00;
    COFLAM[6] = -3.63206306E-01;
    COFLAM[7] = 1.58671021E-02;
    COFLAM[8] = 1.98513952E+00;
    COFLAM[9] = 1.79491990E+00;
    COFLAM[10] = -1.54008440E-01;
    COFLAM[11] = 6.86719439E-03;
    COFLAM[12] = -3.01284291E+00;
    COFLAM[13] = 3.37554994E+00;
    COFLAM[14] = -3.43353119E-01;
    COFLAM[15] = 1.51043444E-02;
    COFLAM[16] = 1.53490799E+01;
    COFLAM[17] = -3.77958145E+00;
    COFLAM[18] = 6.13516524E-01;
    COFLAM[19] = -2.72295753E-02;
    COFLAM[20] = 2.28195645E+01;
    COFLAM[21] = -8.72278946E+00;
    COFLAM[22] = 1.49300487E+00;
    COFLAM[23] = -7.41524047E-02;
    COFLAM[24] = 5.56023763E-01;
    COFLAM[25] = 1.59073590E+00;
    COFLAM[26] = -5.28053839E-02;
    COFLAM[27] = 4.07601571E-04;
    COFLAM[28] = 6.27051982E-01;
    COFLAM[29] = 1.43139617E+00;
    COFLAM[30] = 1.80509282E-03;
    COFLAM[31] = -3.55624900E-03;
    COFLAM[32] = 1.15507063E+01;
    COFLAM[33] = -2.91452379E+00;
    COFLAM[34] = 5.55043580E-01;
    COFLAM[35] = -2.75172461E-02;
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
    COFD[8] = -1.06250182E+01;
    COFD[9] = 2.15849701E+00;
    COFD[10] = -6.53886401E-02;
    COFD[11] = 2.81453370E-03;
    COFD[12] = -1.15797750E+01;
    COFD[13] = 2.43235504E+00;
    COFD[14] = -1.02890179E-01;
    COFD[15] = 4.52903603E-03;
    COFD[16] = -1.06283453E+01;
    COFD[17] = 2.15849701E+00;
    COFD[18] = -6.53886401E-02;
    COFD[19] = 2.81453370E-03;
    COFD[20] = -1.68758926E+01;
    COFD[21] = 4.49460303E+00;
    COFD[22] = -3.64766132E-01;
    COFD[23] = 1.56457153E-02;
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
    COFD[36] = -1.11808682E+01;
    COFD[37] = 2.66936727E+00;
    COFD[38] = -1.34411514E-01;
    COFD[39] = 5.92957488E-03;
    COFD[40] = -1.43693056E+01;
    COFD[41] = 4.03992999E+00;
    COFD[42] = -3.08044800E-01;
    COFD[43] = 1.32757775E-02;
    COFD[44] = -1.31860117E+01;
    COFD[45] = 3.38003453E+00;
    COFD[46] = -2.25783856E-01;
    COFD[47] = 9.85028660E-03;
    COFD[48] = -1.43712864E+01;
    COFD[49] = 3.70920439E+00;
    COFD[50] = -2.67274113E-01;
    COFD[51] = 1.15967481E-02;
    COFD[52] = -1.31877711E+01;
    COFD[53] = 3.38003453E+00;
    COFD[54] = -2.25783856E-01;
    COFD[55] = 9.85028660E-03;
    COFD[56] = -1.93611051E+01;
    COFD[57] = 5.51579726E+00;
    COFD[58] = -4.76061961E-01;
    COFD[59] = 1.96329391E-02;
    COFD[60] = -1.43717529E+01;
    COFD[61] = 3.70920439E+00;
    COFD[62] = -2.67274113E-01;
    COFD[63] = 1.15967481E-02;
    COFD[64] = -1.43721922E+01;
    COFD[65] = 3.70920439E+00;
    COFD[66] = -2.67274113E-01;
    COFD[67] = 1.15967481E-02;
    COFD[68] = -1.40298830E+01;
    COFD[69] = 3.55837688E+00;
    COFD[70] = -2.47785790E-01;
    COFD[71] = 1.07555332E-02;
    COFD[72] = -1.06250182E+01;
    COFD[73] = 2.15849701E+00;
    COFD[74] = -6.53886401E-02;
    COFD[75] = 2.81453370E-03;
    COFD[76] = -1.31860117E+01;
    COFD[77] = 3.38003453E+00;
    COFD[78] = -2.25783856E-01;
    COFD[79] = 9.85028660E-03;
    COFD[80] = -1.29877365E+01;
    COFD[81] = 2.80841511E+00;
    COFD[82] = -1.52629888E-01;
    COFD[83] = 6.72604927E-03;
    COFD[84] = -1.40864894E+01;
    COFD[85] = 3.07458927E+00;
    COFD[86] = -1.86899591E-01;
    COFD[87] = 8.19829781E-03;
    COFD[88] = -1.30027772E+01;
    COFD[89] = 2.80841511E+00;
    COFD[90] = -1.52629888E-01;
    COFD[91] = 6.72604927E-03;
    COFD[92] = -1.91096797E+01;
    COFD[93] = 5.02608697E+00;
    COFD[94] = -4.26959993E-01;
    COFD[95] = 1.80709910E-02;
    COFD[96] = -1.40916052E+01;
    COFD[97] = 3.07458927E+00;
    COFD[98] = -1.86899591E-01;
    COFD[99] = 8.19829781E-03;
    COFD[100] = -1.40964661E+01;
    COFD[101] = 3.07458927E+00;
    COFD[102] = -1.86899591E-01;
    COFD[103] = 8.19829781E-03;
    COFD[104] = -1.38756407E+01;
    COFD[105] = 2.98558426E+00;
    COFD[106] = -1.75507216E-01;
    COFD[107] = 7.71173691E-03;
    COFD[108] = -1.15797750E+01;
    COFD[109] = 2.43235504E+00;
    COFD[110] = -1.02890179E-01;
    COFD[111] = 4.52903603E-03;
    COFD[112] = -1.43712864E+01;
    COFD[113] = 3.70920439E+00;
    COFD[114] = -2.67274113E-01;
    COFD[115] = 1.15967481E-02;
    COFD[116] = -1.40864894E+01;
    COFD[117] = 3.07458927E+00;
    COFD[118] = -1.86899591E-01;
    COFD[119] = 8.19829781E-03;
    COFD[120] = -1.53110708E+01;
    COFD[121] = 3.37317428E+00;
    COFD[122] = -2.24900439E-01;
    COFD[123] = 9.81228151E-03;
    COFD[124] = -1.41066459E+01;
    COFD[125] = 3.07458927E+00;
    COFD[126] = -1.86899591E-01;
    COFD[127] = 8.19829781E-03;
    COFD[128] = -2.10640014E+01;
    COFD[129] = 5.50980695E+00;
    COFD[130] = -4.78335488E-01;
    COFD[131] = 1.98515434E-02;
    COFD[132] = -1.53187643E+01;
    COFD[133] = 3.37317428E+00;
    COFD[134] = -2.24900439E-01;
    COFD[135] = 9.81228151E-03;
    COFD[136] = -1.53261114E+01;
    COFD[137] = 3.37317428E+00;
    COFD[138] = -2.24900439E-01;
    COFD[139] = 9.81228151E-03;
    COFD[140] = -1.50096240E+01;
    COFD[141] = 3.25515933E+00;
    COFD[142] = -2.09710110E-01;
    COFD[143] = 9.15941830E-03;
    COFD[144] = -1.06283453E+01;
    COFD[145] = 2.15849701E+00;
    COFD[146] = -6.53886401E-02;
    COFD[147] = 2.81453370E-03;
    COFD[148] = -1.31877711E+01;
    COFD[149] = 3.38003453E+00;
    COFD[150] = -2.25783856E-01;
    COFD[151] = 9.85028660E-03;
    COFD[152] = -1.30027772E+01;
    COFD[153] = 2.80841511E+00;
    COFD[154] = -1.52629888E-01;
    COFD[155] = 6.72604927E-03;
    COFD[156] = -1.41066459E+01;
    COFD[157] = 3.07458927E+00;
    COFD[158] = -1.86899591E-01;
    COFD[159] = 8.19829781E-03;
    COFD[160] = -1.30182843E+01;
    COFD[161] = 2.80841511E+00;
    COFD[162] = -1.52629888E-01;
    COFD[163] = 6.72604927E-03;
    COFD[164] = -1.91256261E+01;
    COFD[165] = 5.02608697E+00;
    COFD[166] = -4.26959993E-01;
    COFD[167] = 1.80709910E-02;
    COFD[168] = -1.41119732E+01;
    COFD[169] = 3.07458927E+00;
    COFD[170] = -1.86899591E-01;
    COFD[171] = 8.19829781E-03;
    COFD[172] = -1.41170372E+01;
    COFD[173] = 3.07458927E+00;
    COFD[174] = -1.86899591E-01;
    COFD[175] = 8.19829781E-03;
    COFD[176] = -1.38948667E+01;
    COFD[177] = 2.98558426E+00;
    COFD[178] = -1.75507216E-01;
    COFD[179] = 7.71173691E-03;
    COFD[180] = -1.68758926E+01;
    COFD[181] = 4.49460303E+00;
    COFD[182] = -3.64766132E-01;
    COFD[183] = 1.56457153E-02;
    COFD[184] = -1.93611051E+01;
    COFD[185] = 5.51579726E+00;
    COFD[186] = -4.76061961E-01;
    COFD[187] = 1.96329391E-02;
    COFD[188] = -1.91096797E+01;
    COFD[189] = 5.02608697E+00;
    COFD[190] = -4.26959993E-01;
    COFD[191] = 1.80709910E-02;
    COFD[192] = -2.10640014E+01;
    COFD[193] = 5.50980695E+00;
    COFD[194] = -4.78335488E-01;
    COFD[195] = 1.98515434E-02;
    COFD[196] = -1.91256261E+01;
    COFD[197] = 5.02608697E+00;
    COFD[198] = -4.26959993E-01;
    COFD[199] = 1.80709910E-02;
    COFD[200] = -1.31492641E+01;
    COFD[201] = 1.48004311E+00;
    COFD[202] = 1.60499553E-01;
    COFD[203] = -1.19765679E-02;
    COFD[204] = -2.04177482E+01;
    COFD[205] = 5.31457079E+00;
    COFD[206] = -4.58216496E-01;
    COFD[207] = 1.91825910E-02;
    COFD[208] = -2.04230073E+01;
    COFD[209] = 5.31457079E+00;
    COFD[210] = -4.58216496E-01;
    COFD[211] = 1.91825910E-02;
    COFD[212] = -2.08123325E+01;
    COFD[213] = 5.42470154E+00;
    COFD[214] = -4.69700416E-01;
    COFD[215] = 1.95706904E-02;
    COFD[216] = -1.15806808E+01;
    COFD[217] = 2.43235504E+00;
    COFD[218] = -1.02890179E-01;
    COFD[219] = 4.52903603E-03;
    COFD[220] = -1.43717529E+01;
    COFD[221] = 3.70920439E+00;
    COFD[222] = -2.67274113E-01;
    COFD[223] = 1.15967481E-02;
    COFD[224] = -1.40916052E+01;
    COFD[225] = 3.07458927E+00;
    COFD[226] = -1.86899591E-01;
    COFD[227] = 8.19829781E-03;
    COFD[228] = -1.53187643E+01;
    COFD[229] = 3.37317428E+00;
    COFD[230] = -2.24900439E-01;
    COFD[231] = 9.81228151E-03;
    COFD[232] = -1.41119732E+01;
    COFD[233] = 3.07458927E+00;
    COFD[234] = -1.86899591E-01;
    COFD[235] = 8.19829781E-03;
    COFD[236] = -2.04177482E+01;
    COFD[237] = 5.31457079E+00;
    COFD[238] = -4.58216496E-01;
    COFD[239] = 1.91825910E-02;
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
    COFD[256] = -1.43721922E+01;
    COFD[257] = 3.70920439E+00;
    COFD[258] = -2.67274113E-01;
    COFD[259] = 1.15967481E-02;
    COFD[260] = -1.40964661E+01;
    COFD[261] = 3.07458927E+00;
    COFD[262] = -1.86899591E-01;
    COFD[263] = 8.19829781E-03;
    COFD[264] = -1.53261114E+01;
    COFD[265] = 3.37317428E+00;
    COFD[266] = -2.24900439E-01;
    COFD[267] = 9.81228151E-03;
    COFD[268] = -1.41170372E+01;
    COFD[269] = 3.07458927E+00;
    COFD[270] = -1.86899591E-01;
    COFD[271] = 8.19829781E-03;
    COFD[272] = -2.04230073E+01;
    COFD[273] = 5.31457079E+00;
    COFD[274] = -4.58216496E-01;
    COFD[275] = 1.91825910E-02;
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
    COFD[292] = -1.40298830E+01;
    COFD[293] = 3.55837688E+00;
    COFD[294] = -2.47785790E-01;
    COFD[295] = 1.07555332E-02;
    COFD[296] = -1.38756407E+01;
    COFD[297] = 2.98558426E+00;
    COFD[298] = -1.75507216E-01;
    COFD[299] = 7.71173691E-03;
    COFD[300] = -1.50096240E+01;
    COFD[301] = 3.25515933E+00;
    COFD[302] = -2.09710110E-01;
    COFD[303] = 9.15941830E-03;
    COFD[304] = -1.38948667E+01;
    COFD[305] = 2.98558426E+00;
    COFD[306] = -1.75507216E-01;
    COFD[307] = 7.71173691E-03;
    COFD[308] = -2.08123325E+01;
    COFD[309] = 5.42470154E+00;
    COFD[310] = -4.69700416E-01;
    COFD[311] = 1.95706904E-02;
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
    KTDIF[0] = 0;
    KTDIF[1] = 1;
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
    COFTD[8] = 4.15583337E-01;
    COFTD[9] = 1.09738399E-05;
    COFTD[10] = -3.96021963E-09;
    COFTD[11] = 1.14414443E-12;
    COFTD[12] = 4.42739084E-01;
    COFTD[13] = 7.11770818E-05;
    COFTD[14] = -3.84768062E-08;
    COFTD[15] = 6.86323437E-12;
    COFTD[16] = 4.21932443E-01;
    COFTD[17] = 1.11414935E-05;
    COFTD[18] = -4.02072219E-09;
    COFTD[19] = 1.16162418E-12;
    COFTD[20] = 6.02028221E-02;
    COFTD[21] = 5.61561867E-04;
    COFTD[22] = -2.55372862E-07;
    COFTD[23] = 3.63389913E-11;
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
    COFTD[40] = 0.00000000E+00;
    COFTD[41] = 0.00000000E+00;
    COFTD[42] = 0.00000000E+00;
    COFTD[43] = 0.00000000E+00;
    COFTD[44] = 2.70010150E-01;
    COFTD[45] = 3.61555093E-04;
    COFTD[46] = -1.80744752E-07;
    COFTD[47] = 2.75321248E-11;
    COFTD[48] = 2.20482843E-01;
    COFTD[49] = 4.80164288E-04;
    COFTD[50] = -2.32927944E-07;
    COFTD[51] = 3.46470436E-11;
    COFTD[52] = 2.72041664E-01;
    COFTD[53] = 3.64275376E-04;
    COFTD[54] = -1.82104647E-07;
    COFTD[55] = 2.77392722E-11;
    COFTD[56] = -1.41883744E-01;
    COFTD[57] = 7.66558810E-04;
    COFTD[58] = -3.06550003E-07;
    COFTD[59] = 4.02959502E-11;
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

#ifndef AMREX_USE_GPU
/* Initializes parameter database */
void CKINIT()
{

    // (0):  2.000000 O + M <=> O2 + M
    kiv[1] = {2,3};
    nuv[1] = {-2.0,1};
    // (0):  2.000000 O + M <=> O2 + M
    fwd_A[1]     = 1.2e+17;
    fwd_beta[1]  = -1;
    fwd_Ea[1]    = 0;
    prefactor_units[1]  = 1.0000000000000002e-12;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = pow(10,-12.000000);
    is_PD[1] = 0;
    nTB[1] = 2;
    TB[1] = (amrex::Real *) malloc(2 * sizeof(amrex::Real));
    TBid[1] = (int *) malloc(2 * sizeof(int));
    TBid[1][0] = 0; TB[1][0] = 2.3999999999999999; // H2
    TBid[1][1] = 5; TB[1][1] = 15.4; // H2O

    // (1):  O + H + M <=> OH + M
    kiv[2] = {2,1,4};
    nuv[2] = {-1,-1,1};
    // (1):  O + H + M <=> OH + M
    fwd_A[2]     = 5e+17;
    fwd_beta[2]  = -1;
    fwd_Ea[2]    = 0;
    prefactor_units[2]  = 1.0000000000000002e-12;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-12.000000);
    is_PD[2] = 0;
    nTB[2] = 2;
    TB[2] = (amrex::Real *) malloc(2 * sizeof(amrex::Real));
    TBid[2] = (int *) malloc(2 * sizeof(int));
    TBid[2][0] = 0; TB[2][0] = 2; // H2
    TBid[2][1] = 5; TB[2][1] = 6; // H2O

    // (2):  O + H2 <=> H + OH
    kiv[6] = {2,0,1,4};
    nuv[6] = {-1,-1,1,1};
    // (2):  O + H2 <=> H + OH
    fwd_A[6]     = 50000;
    fwd_beta[6]  = 2.6699999999999999;
    fwd_Ea[6]    = 6290;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = pow(10,-12.000000);
    is_PD[6] = 0;
    nTB[6] = 0;

    // (3):  O + HO2 <=> OH + O2
    kiv[7] = {2,6,4,3};
    nuv[7] = {-1,-1,1,1};
    // (3):  O + HO2 <=> OH + O2
    fwd_A[7]     = 20000000000000;
    fwd_beta[7]  = 0;
    fwd_Ea[7]    = 0;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = pow(10,-12.000000);
    is_PD[7] = 0;
    nTB[7] = 0;

    // (4):  O + H2O2 <=> OH + HO2
    kiv[8] = {2,7,4,6};
    nuv[8] = {-1,-1,1,1};
    // (4):  O + H2O2 <=> OH + HO2
    fwd_A[8]     = 9630000;
    fwd_beta[8]  = 2;
    fwd_Ea[8]    = 4000;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = pow(10,-12.000000);
    is_PD[8] = 0;
    nTB[8] = 0;

    // (5):  H + O2 + M <=> HO2 + M
    kiv[3] = {1,3,6};
    nuv[3] = {-1,-1,1};
    // (5):  H + O2 + M <=> HO2 + M
    fwd_A[3]     = 2.8e+18;
    fwd_beta[3]  = -0.85999999999999999;
    fwd_Ea[3]    = 0;
    prefactor_units[3]  = 1.0000000000000002e-12;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-12.000000);
    is_PD[3] = 0;
    nTB[3] = 3;
    TB[3] = (amrex::Real *) malloc(3 * sizeof(amrex::Real));
    TBid[3] = (int *) malloc(3 * sizeof(int));
    TBid[3][0] = 3; TB[3][0] = 0; // O2
    TBid[3][1] = 5; TB[3][1] = 0; // H2O
    TBid[3][2] = 8; TB[3][2] = 0; // N2

    // (6):  H + 2.000000 O2 <=> HO2 + O2
    kiv[9] = {1,3,6,3};
    nuv[9] = {-1,-2.0,1,1};
    // (6):  H + 2.000000 O2 <=> HO2 + O2
    fwd_A[9]     = 3e+20;
    fwd_beta[9]  = -1.72;
    fwd_Ea[9]    = 0;
    prefactor_units[9]  = 1.0000000000000002e-12;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = pow(10,-18.000000);
    is_PD[9] = 0;
    nTB[9] = 0;

    // (7):  H + O2 + H2O <=> HO2 + H2O
    kiv[10] = {1,3,5,6,5};
    nuv[10] = {-1,-1,-1,1,1};
    // (7):  H + O2 + H2O <=> HO2 + H2O
    fwd_A[10]     = 9.38e+18;
    fwd_beta[10]  = -0.76000000000000001;
    fwd_Ea[10]    = 0;
    prefactor_units[10]  = 1.0000000000000002e-12;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = pow(10,-18.000000);
    is_PD[10] = 0;
    nTB[10] = 0;

    // (8):  H + O2 + N2 <=> HO2 + N2
    kiv[11] = {1,3,8,6,8};
    nuv[11] = {-1,-1,-1,1,1};
    // (8):  H + O2 + N2 <=> HO2 + N2
    fwd_A[11]     = 3.75e+20;
    fwd_beta[11]  = -1.72;
    fwd_Ea[11]    = 0;
    prefactor_units[11]  = 1.0000000000000002e-12;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = pow(10,-18.000000);
    is_PD[11] = 0;
    nTB[11] = 0;

    // (9):  H + O2 <=> O + OH
    kiv[12] = {1,3,2,4};
    nuv[12] = {-1,-1,1,1};
    // (9):  H + O2 <=> O + OH
    fwd_A[12]     = 83000000000000;
    fwd_beta[12]  = 0;
    fwd_Ea[12]    = 14413;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = pow(10,-12.000000);
    is_PD[12] = 0;
    nTB[12] = 0;

    // (10):  2.000000 H + M <=> H2 + M
    kiv[4] = {1,0};
    nuv[4] = {-2.0,1};
    // (10):  2.000000 H + M <=> H2 + M
    fwd_A[4]     = 1e+18;
    fwd_beta[4]  = -1;
    fwd_Ea[4]    = 0;
    prefactor_units[4]  = 1.0000000000000002e-12;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 0;
    nTB[4] = 2;
    TB[4] = (amrex::Real *) malloc(2 * sizeof(amrex::Real));
    TBid[4] = (int *) malloc(2 * sizeof(int));
    TBid[4][0] = 0; TB[4][0] = 0; // H2
    TBid[4][1] = 5; TB[4][1] = 0; // H2O

    // (11):  2.000000 H + H2 <=> 2.000000 H2
    kiv[13] = {1,0,0};
    nuv[13] = {-2.0,-1,2.0};
    // (11):  2.000000 H + H2 <=> 2.000000 H2
    fwd_A[13]     = 90000000000000000;
    fwd_beta[13]  = -0.59999999999999998;
    fwd_Ea[13]    = 0;
    prefactor_units[13]  = 1.0000000000000002e-12;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = pow(10,-18.000000);
    is_PD[13] = 0;
    nTB[13] = 0;

    // (12):  2.000000 H + H2O <=> H2 + H2O
    kiv[14] = {1,5,0,5};
    nuv[14] = {-2.0,-1,1,1};
    // (12):  2.000000 H + H2O <=> H2 + H2O
    fwd_A[14]     = 6e+19;
    fwd_beta[14]  = -1.25;
    fwd_Ea[14]    = 0;
    prefactor_units[14]  = 1.0000000000000002e-12;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = pow(10,-18.000000);
    is_PD[14] = 0;
    nTB[14] = 0;

    // (13):  H + OH + M <=> H2O + M
    kiv[5] = {1,4,5};
    nuv[5] = {-1,-1,1};
    // (13):  H + OH + M <=> H2O + M
    fwd_A[5]     = 2.2e+22;
    fwd_beta[5]  = -2;
    fwd_Ea[5]    = 0;
    prefactor_units[5]  = 1.0000000000000002e-12;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-12.000000);
    is_PD[5] = 0;
    nTB[5] = 2;
    TB[5] = (amrex::Real *) malloc(2 * sizeof(amrex::Real));
    TBid[5] = (int *) malloc(2 * sizeof(int));
    TBid[5][0] = 0; TB[5][0] = 0.72999999999999998; // H2
    TBid[5][1] = 5; TB[5][1] = 3.6499999999999999; // H2O

    // (14):  H + HO2 <=> O + H2O
    kiv[15] = {1,6,2,5};
    nuv[15] = {-1,-1,1,1};
    // (14):  H + HO2 <=> O + H2O
    fwd_A[15]     = 3970000000000;
    fwd_beta[15]  = 0;
    fwd_Ea[15]    = 671;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = pow(10,-12.000000);
    is_PD[15] = 0;
    nTB[15] = 0;

    // (15):  H + HO2 <=> O2 + H2
    kiv[16] = {1,6,3,0};
    nuv[16] = {-1,-1,1,1};
    // (15):  H + HO2 <=> O2 + H2
    fwd_A[16]     = 28000000000000;
    fwd_beta[16]  = 0;
    fwd_Ea[16]    = 1068;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = pow(10,-12.000000);
    is_PD[16] = 0;
    nTB[16] = 0;

    // (16):  H + HO2 <=> 2.000000 OH
    kiv[17] = {1,6,4};
    nuv[17] = {-1,-1,2.0};
    // (16):  H + HO2 <=> 2.000000 OH
    fwd_A[17]     = 134000000000000;
    fwd_beta[17]  = 0;
    fwd_Ea[17]    = 635;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = pow(10,-12.000000);
    is_PD[17] = 0;
    nTB[17] = 0;

    // (17):  H + H2O2 <=> HO2 + H2
    kiv[18] = {1,7,6,0};
    nuv[18] = {-1,-1,1,1};
    // (17):  H + H2O2 <=> HO2 + H2
    fwd_A[18]     = 12100000;
    fwd_beta[18]  = 2;
    fwd_Ea[18]    = 5200;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = pow(10,-12.000000);
    is_PD[18] = 0;
    nTB[18] = 0;

    // (18):  H + H2O2 <=> OH + H2O
    kiv[19] = {1,7,4,5};
    nuv[19] = {-1,-1,1,1};
    // (18):  H + H2O2 <=> OH + H2O
    fwd_A[19]     = 10000000000000;
    fwd_beta[19]  = 0;
    fwd_Ea[19]    = 3600;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = pow(10,-12.000000);
    is_PD[19] = 0;
    nTB[19] = 0;

    // (19):  OH + H2 <=> H + H2O
    kiv[20] = {4,0,1,5};
    nuv[20] = {-1,-1,1,1};
    // (19):  OH + H2 <=> H + H2O
    fwd_A[20]     = 216000000;
    fwd_beta[20]  = 1.51;
    fwd_Ea[20]    = 3430;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-12.000000);
    is_PD[20] = 0;
    nTB[20] = 0;

    // (20):  2.000000 OH (+M) <=> H2O2 (+M)
    kiv[0] = {4,7};
    nuv[0] = {-2.0,1};
    // (20):  2.000000 OH (+M) <=> H2O2 (+M)
    fwd_A[0]     = 74000000000000;
    fwd_beta[0]  = -0.37;
    fwd_Ea[0]    = 0;
    low_A[0]     = 2.3e+18;
    low_beta[0]  = -0.90000000000000002;
    low_Ea[0]    = -1700;
    troe_a[0]    = 0.73460000000000003;
    troe_Tsss[0] = 94;
    troe_Ts[0]   = 1756;
    troe_Tss[0]  = 5182;
    troe_len[0]  = 4;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = pow(10,-12.000000);
    is_PD[0] = 1;
    nTB[0] = 2;
    TB[0] = (amrex::Real *) malloc(2 * sizeof(amrex::Real));
    TBid[0] = (int *) malloc(2 * sizeof(int));
    TBid[0][0] = 0; TB[0][0] = 2; // H2
    TBid[0][1] = 5; TB[0][1] = 6; // H2O

    // (21):  2.000000 OH <=> O + H2O
    kiv[21] = {4,2,5};
    nuv[21] = {-2.0,1,1};
    // (21):  2.000000 OH <=> O + H2O
    fwd_A[21]     = 35700;
    fwd_beta[21]  = 2.3999999999999999;
    fwd_Ea[21]    = -2110;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = pow(10,-12.000000);
    is_PD[21] = 0;
    nTB[21] = 0;

    // (22):  OH + HO2 <=> O2 + H2O
    kiv[22] = {4,6,3,5};
    nuv[22] = {-1,-1,1,1};
    // (22):  OH + HO2 <=> O2 + H2O
    fwd_A[22]     = 29000000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = -500;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = pow(10,-12.000000);
    is_PD[22] = 0;
    nTB[22] = 0;

    // (23):  OH + H2O2 <=> HO2 + H2O
    kiv[23] = {4,7,6,5};
    nuv[23] = {-1,-1,1,1};
    // (23):  OH + H2O2 <=> HO2 + H2O
    fwd_A[23]     = 1750000000000;
    fwd_beta[23]  = 0;
    fwd_Ea[23]    = 320;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = pow(10,-12.000000);
    is_PD[23] = 0;
    nTB[23] = 0;

    // (24):  OH + H2O2 <=> HO2 + H2O
    kiv[24] = {4,7,6,5};
    nuv[24] = {-1,-1,1,1};
    // (24):  OH + H2O2 <=> HO2 + H2O
    fwd_A[24]     = 580000000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 9560;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = pow(10,-12.000000);
    is_PD[24] = 0;
    nTB[24] = 0;

    // (25):  2.000000 HO2 <=> O2 + H2O2
    kiv[25] = {6,3,7};
    nuv[25] = {-2.0,1,1};
    // (25):  2.000000 HO2 <=> O2 + H2O2
    fwd_A[25]     = 130000000000;
    fwd_beta[25]  = 0;
    fwd_Ea[25]    = -1630;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = pow(10,-12.000000);
    is_PD[25] = 0;
    nTB[25] = 0;

    // (26):  2.000000 HO2 <=> O2 + H2O2
    kiv[26] = {6,3,7};
    nuv[26] = {-2.0,1,1};
    // (26):  2.000000 HO2 <=> O2 + H2O2
    fwd_A[26]     = 420000000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 12000;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = pow(10,-12.000000);
    is_PD[26] = 0;
    nTB[26] = 0;

}


/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<27; ++i) {
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
    ncf[ 0 * kd + 1 ] = 2; /*H */

    /*H */
    ncf[ 1 * kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 2 * kd + 0 ] = 1; /*O */

    /*O2 */
    ncf[ 3 * kd + 0 ] = 2; /*O */

    /*OH */
    ncf[ 4 * kd + 0 ] = 1; /*O */
    ncf[ 4 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 5 * kd + 1 ] = 2; /*H */
    ncf[ 5 * kd + 0 ] = 1; /*O */

    /*HO2 */
    ncf[ 6 * kd + 1 ] = 1; /*H */
    ncf[ 6 * kd + 0 ] = 2; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 1 ] = 2; /*H */
    ncf[ 7 * kd + 0 ] = 2; /*O */

    /*N2 */
    ncf[ 8 * kd + 2 ] = 2; /*N */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(3);
    ename[0] = "O";
    ename[1] = "H";
    ename[2] = "N";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(9);
    kname[0] = "H2";
    kname[1] = "H";
    kname[2] = "O";
    kname[3] = "O2";
    kname[4] = "OH";
    kname[5] = "H2O";
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
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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

