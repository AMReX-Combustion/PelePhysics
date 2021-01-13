#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/* PURE CPU stuff  */
#ifndef AMREX_USE_GPU
namespace thermo
{
    amrex::Real fwd_A[124], fwd_beta[124], fwd_Ea[124];
    amrex::Real low_A[124], low_beta[124], low_Ea[124];
    amrex::Real rev_A[124], rev_beta[124], rev_Ea[124];
    amrex::Real troe_a[124],troe_Ts[124], troe_Tss[124], troe_Tsss[124];
    amrex::Real sri_a[124], sri_b[124], sri_c[124], sri_d[124], sri_e[124];
    amrex::Real activation_units[124], prefactor_units[124], phase_units[124];
    int is_PD[124], troe_len[124], sri_len[124], nTB[124], *TBid[124];
    amrex::Real *TB[124];
};

using namespace thermo;


/* Vectorized stuff  */

/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int *  np, amrex::Real *  y,  amrex::Real *  x)
{
    amrex::Real YOW[*np];
    amrex::Real imw[21];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<21; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<21; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}

/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tc[5], h[21];
    amrex::Real imw[21];

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
        hms[14*(*np)+i] = h[14];
        hms[15*(*np)+i] = h[15];
        hms[16*(*np)+i] = h[16];
        hms[17*(*np)+i] = h[17];
        hms[18*(*np)+i] = h[18];
        hms[19*(*np)+i] = h[19];
        hms[20*(*np)+i] = h[20];
    }

    for (int n=0; n<21; n++) {
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
    amrex::Real c[21*(*np)]; /*temporary storage */
    amrex::Real imw[21];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<21; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<21*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW[*np];
    amrex::Real imw[21];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<21; n++) {
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
    amrex::Real k_f_s[124*npt], Kc_s[124*npt], mixture[npt], g_RT[21*npt];
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

    for (int n=0; n<21; n++) {
        for (int i=0; i<npt; i++) {
            mixture[i] += sc[n*npt+i];
            wdot[n*npt+i] = 0.0;
        }
    }

    vcomp_k_f(npt, k_f_s, tc, invT);

    vcomp_gibbs(npt, g_RT, tc);

    vcomp_Kc(npt, Kc_s, g_RT, invT);

    vcomp_wdot_1_50(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
    vcomp_wdot_51_100(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
    vcomp_wdot_101_124(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
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
        k_f_s[38*npt+i] = prefactor_units[38] * fwd_A[38] * exp(fwd_beta[38] * tc[i] - activation_units[38] * fwd_Ea[38] * invT[i]);
        k_f_s[39*npt+i] = prefactor_units[39] * fwd_A[39] * exp(fwd_beta[39] * tc[i] - activation_units[39] * fwd_Ea[39] * invT[i]);
        k_f_s[40*npt+i] = prefactor_units[40] * fwd_A[40] * exp(fwd_beta[40] * tc[i] - activation_units[40] * fwd_Ea[40] * invT[i]);
        k_f_s[41*npt+i] = prefactor_units[41] * fwd_A[41] * exp(fwd_beta[41] * tc[i] - activation_units[41] * fwd_Ea[41] * invT[i]);
        k_f_s[42*npt+i] = prefactor_units[42] * fwd_A[42] * exp(fwd_beta[42] * tc[i] - activation_units[42] * fwd_Ea[42] * invT[i]);
        k_f_s[43*npt+i] = prefactor_units[43] * fwd_A[43] * exp(fwd_beta[43] * tc[i] - activation_units[43] * fwd_Ea[43] * invT[i]);
        k_f_s[44*npt+i] = prefactor_units[44] * fwd_A[44] * exp(fwd_beta[44] * tc[i] - activation_units[44] * fwd_Ea[44] * invT[i]);
        k_f_s[45*npt+i] = prefactor_units[45] * fwd_A[45] * exp(fwd_beta[45] * tc[i] - activation_units[45] * fwd_Ea[45] * invT[i]);
        k_f_s[46*npt+i] = prefactor_units[46] * fwd_A[46] * exp(fwd_beta[46] * tc[i] - activation_units[46] * fwd_Ea[46] * invT[i]);
        k_f_s[47*npt+i] = prefactor_units[47] * fwd_A[47] * exp(fwd_beta[47] * tc[i] - activation_units[47] * fwd_Ea[47] * invT[i]);
        k_f_s[48*npt+i] = prefactor_units[48] * fwd_A[48] * exp(fwd_beta[48] * tc[i] - activation_units[48] * fwd_Ea[48] * invT[i]);
        k_f_s[49*npt+i] = prefactor_units[49] * fwd_A[49] * exp(fwd_beta[49] * tc[i] - activation_units[49] * fwd_Ea[49] * invT[i]);
        k_f_s[50*npt+i] = prefactor_units[50] * fwd_A[50] * exp(fwd_beta[50] * tc[i] - activation_units[50] * fwd_Ea[50] * invT[i]);
        k_f_s[51*npt+i] = prefactor_units[51] * fwd_A[51] * exp(fwd_beta[51] * tc[i] - activation_units[51] * fwd_Ea[51] * invT[i]);
        k_f_s[52*npt+i] = prefactor_units[52] * fwd_A[52] * exp(fwd_beta[52] * tc[i] - activation_units[52] * fwd_Ea[52] * invT[i]);
        k_f_s[53*npt+i] = prefactor_units[53] * fwd_A[53] * exp(fwd_beta[53] * tc[i] - activation_units[53] * fwd_Ea[53] * invT[i]);
        k_f_s[54*npt+i] = prefactor_units[54] * fwd_A[54] * exp(fwd_beta[54] * tc[i] - activation_units[54] * fwd_Ea[54] * invT[i]);
        k_f_s[55*npt+i] = prefactor_units[55] * fwd_A[55] * exp(fwd_beta[55] * tc[i] - activation_units[55] * fwd_Ea[55] * invT[i]);
        k_f_s[56*npt+i] = prefactor_units[56] * fwd_A[56] * exp(fwd_beta[56] * tc[i] - activation_units[56] * fwd_Ea[56] * invT[i]);
        k_f_s[57*npt+i] = prefactor_units[57] * fwd_A[57] * exp(fwd_beta[57] * tc[i] - activation_units[57] * fwd_Ea[57] * invT[i]);
        k_f_s[58*npt+i] = prefactor_units[58] * fwd_A[58] * exp(fwd_beta[58] * tc[i] - activation_units[58] * fwd_Ea[58] * invT[i]);
        k_f_s[59*npt+i] = prefactor_units[59] * fwd_A[59] * exp(fwd_beta[59] * tc[i] - activation_units[59] * fwd_Ea[59] * invT[i]);
        k_f_s[60*npt+i] = prefactor_units[60] * fwd_A[60] * exp(fwd_beta[60] * tc[i] - activation_units[60] * fwd_Ea[60] * invT[i]);
        k_f_s[61*npt+i] = prefactor_units[61] * fwd_A[61] * exp(fwd_beta[61] * tc[i] - activation_units[61] * fwd_Ea[61] * invT[i]);
        k_f_s[62*npt+i] = prefactor_units[62] * fwd_A[62] * exp(fwd_beta[62] * tc[i] - activation_units[62] * fwd_Ea[62] * invT[i]);
        k_f_s[63*npt+i] = prefactor_units[63] * fwd_A[63] * exp(fwd_beta[63] * tc[i] - activation_units[63] * fwd_Ea[63] * invT[i]);
        k_f_s[64*npt+i] = prefactor_units[64] * fwd_A[64] * exp(fwd_beta[64] * tc[i] - activation_units[64] * fwd_Ea[64] * invT[i]);
        k_f_s[65*npt+i] = prefactor_units[65] * fwd_A[65] * exp(fwd_beta[65] * tc[i] - activation_units[65] * fwd_Ea[65] * invT[i]);
        k_f_s[66*npt+i] = prefactor_units[66] * fwd_A[66] * exp(fwd_beta[66] * tc[i] - activation_units[66] * fwd_Ea[66] * invT[i]);
        k_f_s[67*npt+i] = prefactor_units[67] * fwd_A[67] * exp(fwd_beta[67] * tc[i] - activation_units[67] * fwd_Ea[67] * invT[i]);
        k_f_s[68*npt+i] = prefactor_units[68] * fwd_A[68] * exp(fwd_beta[68] * tc[i] - activation_units[68] * fwd_Ea[68] * invT[i]);
        k_f_s[69*npt+i] = prefactor_units[69] * fwd_A[69] * exp(fwd_beta[69] * tc[i] - activation_units[69] * fwd_Ea[69] * invT[i]);
        k_f_s[70*npt+i] = prefactor_units[70] * fwd_A[70] * exp(fwd_beta[70] * tc[i] - activation_units[70] * fwd_Ea[70] * invT[i]);
        k_f_s[71*npt+i] = prefactor_units[71] * fwd_A[71] * exp(fwd_beta[71] * tc[i] - activation_units[71] * fwd_Ea[71] * invT[i]);
        k_f_s[72*npt+i] = prefactor_units[72] * fwd_A[72] * exp(fwd_beta[72] * tc[i] - activation_units[72] * fwd_Ea[72] * invT[i]);
        k_f_s[73*npt+i] = prefactor_units[73] * fwd_A[73] * exp(fwd_beta[73] * tc[i] - activation_units[73] * fwd_Ea[73] * invT[i]);
        k_f_s[74*npt+i] = prefactor_units[74] * fwd_A[74] * exp(fwd_beta[74] * tc[i] - activation_units[74] * fwd_Ea[74] * invT[i]);
        k_f_s[75*npt+i] = prefactor_units[75] * fwd_A[75] * exp(fwd_beta[75] * tc[i] - activation_units[75] * fwd_Ea[75] * invT[i]);
        k_f_s[76*npt+i] = prefactor_units[76] * fwd_A[76] * exp(fwd_beta[76] * tc[i] - activation_units[76] * fwd_Ea[76] * invT[i]);
        k_f_s[77*npt+i] = prefactor_units[77] * fwd_A[77] * exp(fwd_beta[77] * tc[i] - activation_units[77] * fwd_Ea[77] * invT[i]);
        k_f_s[78*npt+i] = prefactor_units[78] * fwd_A[78] * exp(fwd_beta[78] * tc[i] - activation_units[78] * fwd_Ea[78] * invT[i]);
        k_f_s[79*npt+i] = prefactor_units[79] * fwd_A[79] * exp(fwd_beta[79] * tc[i] - activation_units[79] * fwd_Ea[79] * invT[i]);
        k_f_s[80*npt+i] = prefactor_units[80] * fwd_A[80] * exp(fwd_beta[80] * tc[i] - activation_units[80] * fwd_Ea[80] * invT[i]);
        k_f_s[81*npt+i] = prefactor_units[81] * fwd_A[81] * exp(fwd_beta[81] * tc[i] - activation_units[81] * fwd_Ea[81] * invT[i]);
        k_f_s[82*npt+i] = prefactor_units[82] * fwd_A[82] * exp(fwd_beta[82] * tc[i] - activation_units[82] * fwd_Ea[82] * invT[i]);
        k_f_s[83*npt+i] = prefactor_units[83] * fwd_A[83] * exp(fwd_beta[83] * tc[i] - activation_units[83] * fwd_Ea[83] * invT[i]);
        k_f_s[84*npt+i] = prefactor_units[84] * fwd_A[84] * exp(fwd_beta[84] * tc[i] - activation_units[84] * fwd_Ea[84] * invT[i]);
        k_f_s[85*npt+i] = prefactor_units[85] * fwd_A[85] * exp(fwd_beta[85] * tc[i] - activation_units[85] * fwd_Ea[85] * invT[i]);
        k_f_s[86*npt+i] = prefactor_units[86] * fwd_A[86] * exp(fwd_beta[86] * tc[i] - activation_units[86] * fwd_Ea[86] * invT[i]);
        k_f_s[87*npt+i] = prefactor_units[87] * fwd_A[87] * exp(fwd_beta[87] * tc[i] - activation_units[87] * fwd_Ea[87] * invT[i]);
        k_f_s[88*npt+i] = prefactor_units[88] * fwd_A[88] * exp(fwd_beta[88] * tc[i] - activation_units[88] * fwd_Ea[88] * invT[i]);
        k_f_s[89*npt+i] = prefactor_units[89] * fwd_A[89] * exp(fwd_beta[89] * tc[i] - activation_units[89] * fwd_Ea[89] * invT[i]);
        k_f_s[90*npt+i] = prefactor_units[90] * fwd_A[90] * exp(fwd_beta[90] * tc[i] - activation_units[90] * fwd_Ea[90] * invT[i]);
        k_f_s[91*npt+i] = prefactor_units[91] * fwd_A[91] * exp(fwd_beta[91] * tc[i] - activation_units[91] * fwd_Ea[91] * invT[i]);
        k_f_s[92*npt+i] = prefactor_units[92] * fwd_A[92] * exp(fwd_beta[92] * tc[i] - activation_units[92] * fwd_Ea[92] * invT[i]);
        k_f_s[93*npt+i] = prefactor_units[93] * fwd_A[93] * exp(fwd_beta[93] * tc[i] - activation_units[93] * fwd_Ea[93] * invT[i]);
        k_f_s[94*npt+i] = prefactor_units[94] * fwd_A[94] * exp(fwd_beta[94] * tc[i] - activation_units[94] * fwd_Ea[94] * invT[i]);
        k_f_s[95*npt+i] = prefactor_units[95] * fwd_A[95] * exp(fwd_beta[95] * tc[i] - activation_units[95] * fwd_Ea[95] * invT[i]);
        k_f_s[96*npt+i] = prefactor_units[96] * fwd_A[96] * exp(fwd_beta[96] * tc[i] - activation_units[96] * fwd_Ea[96] * invT[i]);
        k_f_s[97*npt+i] = prefactor_units[97] * fwd_A[97] * exp(fwd_beta[97] * tc[i] - activation_units[97] * fwd_Ea[97] * invT[i]);
        k_f_s[98*npt+i] = prefactor_units[98] * fwd_A[98] * exp(fwd_beta[98] * tc[i] - activation_units[98] * fwd_Ea[98] * invT[i]);
        k_f_s[99*npt+i] = prefactor_units[99] * fwd_A[99] * exp(fwd_beta[99] * tc[i] - activation_units[99] * fwd_Ea[99] * invT[i]);
        k_f_s[100*npt+i] = prefactor_units[100] * fwd_A[100] * exp(fwd_beta[100] * tc[i] - activation_units[100] * fwd_Ea[100] * invT[i]);
        k_f_s[101*npt+i] = prefactor_units[101] * fwd_A[101] * exp(fwd_beta[101] * tc[i] - activation_units[101] * fwd_Ea[101] * invT[i]);
        k_f_s[102*npt+i] = prefactor_units[102] * fwd_A[102] * exp(fwd_beta[102] * tc[i] - activation_units[102] * fwd_Ea[102] * invT[i]);
        k_f_s[103*npt+i] = prefactor_units[103] * fwd_A[103] * exp(fwd_beta[103] * tc[i] - activation_units[103] * fwd_Ea[103] * invT[i]);
        k_f_s[104*npt+i] = prefactor_units[104] * fwd_A[104] * exp(fwd_beta[104] * tc[i] - activation_units[104] * fwd_Ea[104] * invT[i]);
        k_f_s[105*npt+i] = prefactor_units[105] * fwd_A[105] * exp(fwd_beta[105] * tc[i] - activation_units[105] * fwd_Ea[105] * invT[i]);
        k_f_s[106*npt+i] = prefactor_units[106] * fwd_A[106] * exp(fwd_beta[106] * tc[i] - activation_units[106] * fwd_Ea[106] * invT[i]);
        k_f_s[107*npt+i] = prefactor_units[107] * fwd_A[107] * exp(fwd_beta[107] * tc[i] - activation_units[107] * fwd_Ea[107] * invT[i]);
        k_f_s[108*npt+i] = prefactor_units[108] * fwd_A[108] * exp(fwd_beta[108] * tc[i] - activation_units[108] * fwd_Ea[108] * invT[i]);
        k_f_s[109*npt+i] = prefactor_units[109] * fwd_A[109] * exp(fwd_beta[109] * tc[i] - activation_units[109] * fwd_Ea[109] * invT[i]);
        k_f_s[110*npt+i] = prefactor_units[110] * fwd_A[110] * exp(fwd_beta[110] * tc[i] - activation_units[110] * fwd_Ea[110] * invT[i]);
        k_f_s[111*npt+i] = prefactor_units[111] * fwd_A[111] * exp(fwd_beta[111] * tc[i] - activation_units[111] * fwd_Ea[111] * invT[i]);
        k_f_s[112*npt+i] = prefactor_units[112] * fwd_A[112] * exp(fwd_beta[112] * tc[i] - activation_units[112] * fwd_Ea[112] * invT[i]);
        k_f_s[113*npt+i] = prefactor_units[113] * fwd_A[113] * exp(fwd_beta[113] * tc[i] - activation_units[113] * fwd_Ea[113] * invT[i]);
        k_f_s[114*npt+i] = prefactor_units[114] * fwd_A[114] * exp(fwd_beta[114] * tc[i] - activation_units[114] * fwd_Ea[114] * invT[i]);
        k_f_s[115*npt+i] = prefactor_units[115] * fwd_A[115] * exp(fwd_beta[115] * tc[i] - activation_units[115] * fwd_Ea[115] * invT[i]);
        k_f_s[116*npt+i] = prefactor_units[116] * fwd_A[116] * exp(fwd_beta[116] * tc[i] - activation_units[116] * fwd_Ea[116] * invT[i]);
        k_f_s[117*npt+i] = prefactor_units[117] * fwd_A[117] * exp(fwd_beta[117] * tc[i] - activation_units[117] * fwd_Ea[117] * invT[i]);
        k_f_s[118*npt+i] = prefactor_units[118] * fwd_A[118] * exp(fwd_beta[118] * tc[i] - activation_units[118] * fwd_Ea[118] * invT[i]);
        k_f_s[119*npt+i] = prefactor_units[119] * fwd_A[119] * exp(fwd_beta[119] * tc[i] - activation_units[119] * fwd_Ea[119] * invT[i]);
        k_f_s[120*npt+i] = prefactor_units[120] * fwd_A[120] * exp(fwd_beta[120] * tc[i] - activation_units[120] * fwd_Ea[120] * invT[i]);
        k_f_s[121*npt+i] = prefactor_units[121] * fwd_A[121] * exp(fwd_beta[121] * tc[i] - activation_units[121] * fwd_Ea[121] * invT[i]);
        k_f_s[122*npt+i] = prefactor_units[122] * fwd_A[122] * exp(fwd_beta[122] * tc[i] - activation_units[122] * fwd_Ea[122] * invT[i]);
        k_f_s[123*npt+i] = prefactor_units[123] * fwd_A[123] * exp(fwd_beta[123] * tc[i] - activation_units[123] * fwd_Ea[123] * invT[i]);
    }
}

void vcomp_gibbs(int npt, amrex::Real *  g_RT, amrex::Real *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        amrex::Real tg[5], g[21];
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
        g_RT[19*npt+i] = g[19];
        g_RT[20*npt+i] = g[20];
    }
}

void vcomp_Kc(int npt, amrex::Real *  Kc_s, amrex::Real *  g_RT, amrex::Real *  invT)
{
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        amrex::Real refC = (101325. / 8.31451) * invT[i];
        amrex::Real refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[4*npt+i]) - (g_RT[7*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[10*npt+i]) - (g_RT[13*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[11*npt+i]) - (g_RT[13*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[6*npt+i] + g_RT[12*npt+i]) - (g_RT[19*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[15*npt+i]) - (g_RT[16*npt+i]));
        Kc_s[5*npt+i] = refC * exp((g_RT[16*npt+i]) - (g_RT[1*npt+i] + g_RT[8*npt+i]));
        Kc_s[6*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[13*npt+i]) - (g_RT[14*npt+i]));
        Kc_s[7*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[13*npt+i]) - (g_RT[19*npt+i]));
        Kc_s[8*npt+i] = refC * exp((g_RT[18*npt+i]) - (g_RT[2*npt+i] + g_RT[16*npt+i]));
        Kc_s[9*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[18*npt+i]) - (g_RT[19*npt+i]));
        Kc_s[10*npt+i] = refC * exp((g_RT[17*npt+i]) - (g_RT[2*npt+i] + g_RT[16*npt+i]));
        Kc_s[11*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[17*npt+i]) - (g_RT[19*npt+i]));
        Kc_s[12*npt+i] = refCinv * exp((g_RT[3*npt+i] + g_RT[8*npt+i]) - (g_RT[9*npt+i]));
        Kc_s[13*npt+i] = refC * exp((g_RT[1*npt+i]) - (2.000000 * g_RT[2*npt+i]));
        Kc_s[14*npt+i] = refCinv * exp((2.000000 * g_RT[3*npt+i]) - (g_RT[4*npt+i]));
        Kc_s[15*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[3*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[16*npt+i] = refC * exp((g_RT[6*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[17*npt+i] = refC * exp((g_RT[15*npt+i]) - (g_RT[2*npt+i] + g_RT[8*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[2*npt+i] + g_RT[4*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[21*npt+i] = exp((g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[6*npt+i]));
        Kc_s[22*npt+i] = exp((2.000000 * g_RT[5*npt+i]) - (g_RT[3*npt+i] + g_RT[6*npt+i]));
        Kc_s[23*npt+i] = refC * exp((g_RT[1*npt+i] + g_RT[20*npt+i]) - (2.000000 * g_RT[2*npt+i] + g_RT[20*npt+i]));
        Kc_s[24*npt+i] = refCinv * exp((2.000000 * g_RT[3*npt+i] + g_RT[20*npt+i]) - (g_RT[4*npt+i] + g_RT[20*npt+i]));
        Kc_s[25*npt+i] = refC * exp((2.000000 * g_RT[6*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (2.000000 * g_RT[5*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[6*npt+i]));
        Kc_s[29*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[30*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[31*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[32*npt+i] = exp((g_RT[4*npt+i] + g_RT[8*npt+i]) - (g_RT[3*npt+i] + g_RT[9*npt+i]));
        Kc_s[33*npt+i] = exp((g_RT[5*npt+i] + g_RT[8*npt+i]) - (g_RT[2*npt+i] + g_RT[9*npt+i]));
        Kc_s[34*npt+i] = exp((g_RT[5*npt+i] + g_RT[8*npt+i]) - (g_RT[2*npt+i] + g_RT[9*npt+i]));
        Kc_s[35*npt+i] = exp((g_RT[7*npt+i] + g_RT[8*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
        Kc_s[36*npt+i] = exp((g_RT[2*npt+i] + g_RT[15*npt+i]) - (g_RT[1*npt+i] + g_RT[8*npt+i]));
        Kc_s[37*npt+i] = exp((g_RT[3*npt+i] + g_RT[15*npt+i]) - (g_RT[5*npt+i] + g_RT[8*npt+i]));
        Kc_s[38*npt+i] = exp((g_RT[3*npt+i] + g_RT[15*npt+i]) - (g_RT[2*npt+i] + g_RT[9*npt+i]));
        Kc_s[39*npt+i] = exp((g_RT[5*npt+i] + g_RT[15*npt+i]) - (g_RT[6*npt+i] + g_RT[8*npt+i]));
        Kc_s[40*npt+i] = exp((g_RT[4*npt+i] + g_RT[15*npt+i]) - (g_RT[7*npt+i] + g_RT[8*npt+i]));
        Kc_s[41*npt+i] = exp((g_RT[3*npt+i] + g_RT[10*npt+i]) - (g_RT[2*npt+i] + g_RT[8*npt+i]));
        Kc_s[42*npt+i] = exp((g_RT[5*npt+i] + g_RT[10*npt+i]) - (g_RT[2*npt+i] + g_RT[15*npt+i]));
        Kc_s[43*npt+i] = exp((g_RT[1*npt+i] + g_RT[10*npt+i]) - (g_RT[2*npt+i] + g_RT[11*npt+i]));
        Kc_s[44*npt+i] = exp((g_RT[6*npt+i] + g_RT[10*npt+i]) - (g_RT[2*npt+i] + g_RT[16*npt+i]));
        Kc_s[45*npt+i] = exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[3*npt+i] + g_RT[15*npt+i]));
        Kc_s[46*npt+i] = exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[2*npt+i] + g_RT[9*npt+i]));
        Kc_s[47*npt+i] = exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[5*npt+i] + g_RT[8*npt+i]));
        Kc_s[48*npt+i] = refC * exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[2*npt+i] + g_RT[3*npt+i] + g_RT[8*npt+i]));
        Kc_s[49*npt+i] = exp((g_RT[9*npt+i] + g_RT[10*npt+i]) - (g_RT[8*npt+i] + g_RT[15*npt+i]));
        Kc_s[50*npt+i] = refC * exp((g_RT[3*npt+i] + g_RT[11*npt+i]) - (2.000000 * g_RT[2*npt+i] + g_RT[8*npt+i]));
        Kc_s[51*npt+i] = exp((g_RT[5*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[16*npt+i]));
        Kc_s[52*npt+i] = exp((g_RT[5*npt+i] + g_RT[11*npt+i]) - (g_RT[6*npt+i] + g_RT[10*npt+i]));
        Kc_s[53*npt+i] = exp((g_RT[7*npt+i] + g_RT[11*npt+i]) - (g_RT[5*npt+i] + g_RT[16*npt+i]));
        Kc_s[54*npt+i] = exp((g_RT[1*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[13*npt+i]));
        Kc_s[55*npt+i] = refC * exp((g_RT[4*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i] + g_RT[8*npt+i]));
        Kc_s[56*npt+i] = refC * exp((g_RT[4*npt+i] + g_RT[11*npt+i]) - (2.000000 * g_RT[2*npt+i] + g_RT[9*npt+i]));
        Kc_s[57*npt+i] = exp((g_RT[4*npt+i] + g_RT[11*npt+i]) - (g_RT[3*npt+i] + g_RT[16*npt+i]));
        Kc_s[58*npt+i] = exp((g_RT[4*npt+i] + g_RT[11*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i]));
        Kc_s[59*npt+i] = exp((g_RT[4*npt+i] + g_RT[11*npt+i]) - (g_RT[6*npt+i] + g_RT[8*npt+i]));
        Kc_s[60*npt+i] = exp((g_RT[0*npt+i] + g_RT[12*npt+i]) - (g_RT[0*npt+i] + g_RT[11*npt+i]));
        Kc_s[61*npt+i] = exp((g_RT[12*npt+i] + g_RT[20*npt+i]) - (g_RT[11*npt+i] + g_RT[20*npt+i]));
        Kc_s[62*npt+i] = exp((g_RT[2*npt+i] + g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[10*npt+i]));
        Kc_s[63*npt+i] = refC * exp((g_RT[3*npt+i] + g_RT[12*npt+i]) - (2.000000 * g_RT[2*npt+i] + g_RT[8*npt+i]));
        Kc_s[64*npt+i] = exp((g_RT[5*npt+i] + g_RT[12*npt+i]) - (g_RT[2*npt+i] + g_RT[16*npt+i]));
        Kc_s[65*npt+i] = exp((g_RT[1*npt+i] + g_RT[12*npt+i]) - (g_RT[2*npt+i] + g_RT[13*npt+i]));
        Kc_s[66*npt+i] = exp((g_RT[4*npt+i] + g_RT[12*npt+i]) - (g_RT[4*npt+i] + g_RT[11*npt+i]));
        Kc_s[67*npt+i] = exp((g_RT[6*npt+i] + g_RT[12*npt+i]) - (g_RT[6*npt+i] + g_RT[11*npt+i]));
        Kc_s[68*npt+i] = exp((g_RT[6*npt+i] + g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[16*npt+i]));
        Kc_s[69*npt+i] = exp((g_RT[8*npt+i] + g_RT[12*npt+i]) - (g_RT[8*npt+i] + g_RT[11*npt+i]));
        Kc_s[70*npt+i] = exp((g_RT[9*npt+i] + g_RT[12*npt+i]) - (g_RT[9*npt+i] + g_RT[11*npt+i]));
        Kc_s[71*npt+i] = exp((g_RT[9*npt+i] + g_RT[12*npt+i]) - (g_RT[8*npt+i] + g_RT[16*npt+i]));
        Kc_s[72*npt+i] = exp((g_RT[2*npt+i] + g_RT[16*npt+i]) - (g_RT[1*npt+i] + g_RT[15*npt+i]));
        Kc_s[73*npt+i] = exp((g_RT[3*npt+i] + g_RT[16*npt+i]) - (g_RT[5*npt+i] + g_RT[15*npt+i]));
        Kc_s[74*npt+i] = exp((g_RT[5*npt+i] + g_RT[16*npt+i]) - (g_RT[6*npt+i] + g_RT[15*npt+i]));
        Kc_s[75*npt+i] = exp((g_RT[4*npt+i] + g_RT[16*npt+i]) - (g_RT[7*npt+i] + g_RT[15*npt+i]));
        Kc_s[76*npt+i] = exp((g_RT[11*npt+i] + g_RT[16*npt+i]) - (g_RT[13*npt+i] + g_RT[15*npt+i]));
        Kc_s[77*npt+i] = exp((g_RT[12*npt+i] + g_RT[16*npt+i]) - (g_RT[13*npt+i] + g_RT[15*npt+i]));
        Kc_s[78*npt+i] = exp((g_RT[3*npt+i] + g_RT[13*npt+i]) - (g_RT[2*npt+i] + g_RT[16*npt+i]));
        Kc_s[79*npt+i] = refC * exp((g_RT[3*npt+i] + g_RT[13*npt+i]) - (g_RT[1*npt+i] + g_RT[2*npt+i] + g_RT[8*npt+i]));
        Kc_s[80*npt+i] = exp((g_RT[5*npt+i] + g_RT[13*npt+i]) - (g_RT[6*npt+i] + g_RT[11*npt+i]));
        Kc_s[81*npt+i] = exp((g_RT[5*npt+i] + g_RT[13*npt+i]) - (g_RT[6*npt+i] + g_RT[12*npt+i]));
        Kc_s[82*npt+i] = exp((g_RT[5*npt+i] + g_RT[13*npt+i]) - (g_RT[1*npt+i] + g_RT[16*npt+i]));
        Kc_s[83*npt+i] = exp((g_RT[7*npt+i] + g_RT[13*npt+i]) - (g_RT[4*npt+i] + g_RT[14*npt+i]));
        Kc_s[84*npt+i] = exp((g_RT[7*npt+i] + g_RT[13*npt+i]) - (g_RT[5*npt+i] + g_RT[18*npt+i]));
        Kc_s[85*npt+i] = exp((g_RT[4*npt+i] + g_RT[13*npt+i]) - (g_RT[3*npt+i] + g_RT[18*npt+i]));
        Kc_s[86*npt+i] = exp((g_RT[4*npt+i] + g_RT[13*npt+i]) - (g_RT[5*npt+i] + g_RT[16*npt+i]));
        Kc_s[87*npt+i] = exp((g_RT[13*npt+i] + g_RT[15*npt+i]) - (g_RT[8*npt+i] + g_RT[14*npt+i]));
        Kc_s[88*npt+i] = exp((g_RT[13*npt+i] + g_RT[16*npt+i]) - (g_RT[14*npt+i] + g_RT[15*npt+i]));
        Kc_s[89*npt+i] = exp((g_RT[2*npt+i] + g_RT[18*npt+i]) - (g_RT[2*npt+i] + g_RT[17*npt+i]));
        Kc_s[90*npt+i] = exp((g_RT[2*npt+i] + g_RT[18*npt+i]) - (g_RT[1*npt+i] + g_RT[16*npt+i]));
        Kc_s[91*npt+i] = exp((g_RT[2*npt+i] + g_RT[18*npt+i]) - (g_RT[5*npt+i] + g_RT[13*npt+i]));
        Kc_s[92*npt+i] = exp((g_RT[2*npt+i] + g_RT[18*npt+i]) - (g_RT[6*npt+i] + g_RT[12*npt+i]));
        Kc_s[93*npt+i] = exp((g_RT[3*npt+i] + g_RT[18*npt+i]) - (g_RT[5*npt+i] + g_RT[16*npt+i]));
        Kc_s[94*npt+i] = exp((g_RT[5*npt+i] + g_RT[18*npt+i]) - (g_RT[6*npt+i] + g_RT[16*npt+i]));
        Kc_s[95*npt+i] = exp((g_RT[4*npt+i] + g_RT[18*npt+i]) - (g_RT[7*npt+i] + g_RT[16*npt+i]));
        Kc_s[96*npt+i] = exp((g_RT[13*npt+i] + g_RT[18*npt+i]) - (g_RT[14*npt+i] + g_RT[16*npt+i]));
        Kc_s[97*npt+i] = exp((g_RT[8*npt+i] + g_RT[18*npt+i]) - (g_RT[9*npt+i] + g_RT[13*npt+i]));
        Kc_s[98*npt+i] = exp((g_RT[2*npt+i] + g_RT[17*npt+i]) - (g_RT[1*npt+i] + g_RT[16*npt+i]));
        Kc_s[99*npt+i] = exp((g_RT[2*npt+i] + g_RT[17*npt+i]) - (g_RT[5*npt+i] + g_RT[13*npt+i]));
        Kc_s[100*npt+i] = exp((g_RT[2*npt+i] + g_RT[17*npt+i]) - (g_RT[6*npt+i] + g_RT[12*npt+i]));
        Kc_s[101*npt+i] = exp((g_RT[3*npt+i] + g_RT[17*npt+i]) - (g_RT[5*npt+i] + g_RT[16*npt+i]));
        Kc_s[102*npt+i] = exp((g_RT[5*npt+i] + g_RT[17*npt+i]) - (g_RT[6*npt+i] + g_RT[16*npt+i]));
        Kc_s[103*npt+i] = exp((g_RT[4*npt+i] + g_RT[17*npt+i]) - (g_RT[7*npt+i] + g_RT[16*npt+i]));
        Kc_s[104*npt+i] = exp((g_RT[13*npt+i] + g_RT[17*npt+i]) - (g_RT[14*npt+i] + g_RT[16*npt+i]));
        Kc_s[105*npt+i] = exp((g_RT[2*npt+i] + g_RT[14*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[106*npt+i] = exp((g_RT[3*npt+i] + g_RT[14*npt+i]) - (g_RT[5*npt+i] + g_RT[13*npt+i]));
        Kc_s[107*npt+i] = exp((g_RT[5*npt+i] + g_RT[14*npt+i]) - (g_RT[6*npt+i] + g_RT[13*npt+i]));
        Kc_s[108*npt+i] = exp((g_RT[11*npt+i] + g_RT[14*npt+i]) - (2.000000 * g_RT[13*npt+i]));
        Kc_s[109*npt+i] = exp((g_RT[12*npt+i] + g_RT[14*npt+i]) - (2.000000 * g_RT[13*npt+i]));
        Kc_s[110*npt+i] = exp((g_RT[2*npt+i] + g_RT[19*npt+i]) - (g_RT[1*npt+i] + g_RT[17*npt+i]));
        Kc_s[111*npt+i] = exp((g_RT[2*npt+i] + g_RT[19*npt+i]) - (g_RT[1*npt+i] + g_RT[18*npt+i]));
        Kc_s[112*npt+i] = exp((g_RT[3*npt+i] + g_RT[19*npt+i]) - (g_RT[5*npt+i] + g_RT[17*npt+i]));
        Kc_s[113*npt+i] = exp((g_RT[3*npt+i] + g_RT[19*npt+i]) - (g_RT[5*npt+i] + g_RT[18*npt+i]));
        Kc_s[114*npt+i] = exp((g_RT[5*npt+i] + g_RT[19*npt+i]) - (g_RT[6*npt+i] + g_RT[17*npt+i]));
        Kc_s[115*npt+i] = exp((g_RT[5*npt+i] + g_RT[19*npt+i]) - (g_RT[6*npt+i] + g_RT[18*npt+i]));
        Kc_s[116*npt+i] = exp((g_RT[4*npt+i] + g_RT[19*npt+i]) - (g_RT[7*npt+i] + g_RT[17*npt+i]));
        Kc_s[117*npt+i] = exp((g_RT[10*npt+i] + g_RT[19*npt+i]) - (g_RT[13*npt+i] + g_RT[16*npt+i]));
        Kc_s[118*npt+i] = exp((g_RT[11*npt+i] + g_RT[19*npt+i]) - (g_RT[13*npt+i] + g_RT[17*npt+i]));
        Kc_s[119*npt+i] = exp((g_RT[11*npt+i] + g_RT[19*npt+i]) - (g_RT[13*npt+i] + g_RT[18*npt+i]));
        Kc_s[120*npt+i] = exp((g_RT[12*npt+i] + g_RT[19*npt+i]) - (g_RT[13*npt+i] + g_RT[18*npt+i]));
        Kc_s[121*npt+i] = exp((g_RT[12*npt+i] + g_RT[19*npt+i]) - (g_RT[13*npt+i] + g_RT[17*npt+i]));
        Kc_s[122*npt+i] = exp((g_RT[13*npt+i] + g_RT[19*npt+i]) - (g_RT[14*npt+i] + g_RT[17*npt+i]));
        Kc_s[123*npt+i] = exp((g_RT[13*npt+i] + g_RT[19*npt+i]) - (g_RT[14*npt+i] + g_RT[18*npt+i]));
    }
}

void vcomp_wdot_1_50(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
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
        phi_f = sc[2*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[0*npt+i] + (TB[0][1] - 1)*sc[1*npt+i] + (TB[0][2] - 1)*sc[4*npt+i] + (TB[0][3] - 1)*sc[6*npt+i] + (TB[0][4] - 1)*sc[8*npt+i] + (TB[0][5] - 1)*sc[9*npt+i] + (TB[0][6] - 1)*sc[14*npt+i] + (TB[0][7] - 1)*sc[16*npt+i] + (TB[0][8] - 1)*sc[19*npt+i] + (TB[0][9] - 1)*sc[20*npt+i];
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
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 2: CH + H2 (+M) <=> CH3 (+M) */
        phi_f = sc[1*npt+i]*sc[10*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[6*npt+i] + (TB[1][1] - 1)*sc[8*npt+i] + (TB[1][2] - 1)*sc[9*npt+i] + (TB[1][3] - 1)*sc[14*npt+i] + (TB[1][4] - 1)*sc[16*npt+i] + (TB[1][5] - 1)*sc[19*npt+i] + (TB[1][6] - 1)*sc[20*npt+i];
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
        phi_r = sc[13*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 3: CH2 + H (+M) <=> CH3 (+M) */
        phi_f = sc[2*npt+i]*sc[11*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[6*npt+i] + (TB[2][1] - 1)*sc[8*npt+i] + (TB[2][2] - 1)*sc[9*npt+i] + (TB[2][3] - 1)*sc[14*npt+i] + (TB[2][4] - 1)*sc[16*npt+i] + (TB[2][5] - 1)*sc[19*npt+i] + (TB[2][6] - 1)*sc[20*npt+i];
        k_f = k_f_s[2*npt+i];
        redP = alpha / k_f * phase_units[2] * low_A[2] * exp(low_beta[2] * tc[i] - activation_units[2] * low_Ea[2] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[2]) > 1.e-100 ? (1.-troe_a[2])*exp(-T[i]/troe_Tsss[2]) : 0.) 
            + (fabs(troe_Ts[2]) > 1.e-100 ? troe_a[2] * exp(-T[i]/troe_Ts[2]) : 0.) 
            + (troe_len[2] == 4 ? exp(-troe_Tss[2] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 4: CH2(S) + H2O (+M) <=> CH3OH (+M) */
        phi_f = sc[6*npt+i]*sc[12*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[6*npt+i] + (TB[3][1] - 1)*sc[8*npt+i] + (TB[3][2] - 1)*sc[9*npt+i] + (TB[3][3] - 1)*sc[14*npt+i] + (TB[3][4] - 1)*sc[16*npt+i] + (TB[3][5] - 1)*sc[19*npt+i];
        k_f = k_f_s[3*npt+i];
        redP = alpha / k_f * phase_units[3] * low_A[3] * exp(low_beta[3] * tc[i] - activation_units[3] * low_Ea[3] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[3]) > 1.e-100 ? (1.-troe_a[3])*exp(-T[i]/troe_Tsss[3]) : 0.) 
            + (fabs(troe_Ts[3]) > 1.e-100 ? troe_a[3] * exp(-T[i]/troe_Ts[3]) : 0.) 
            + (troe_len[3] == 4 ? exp(-troe_Tss[3] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[19*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 5: HCO + H (+M) <=> CH2O (+M) */
        phi_f = sc[2*npt+i]*sc[15*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[6*npt+i] + (TB[4][1] - 1)*sc[8*npt+i] + (TB[4][2] - 1)*sc[9*npt+i] + (TB[4][3] - 1)*sc[14*npt+i] + (TB[4][4] - 1)*sc[16*npt+i] + (TB[4][5] - 1)*sc[19*npt+i] + (TB[4][6] - 1)*sc[20*npt+i];
        k_f = k_f_s[4*npt+i];
        redP = alpha / k_f * phase_units[4] * low_A[4] * exp(low_beta[4] * tc[i] - activation_units[4] * low_Ea[4] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[4]) > 1.e-100 ? (1.-troe_a[4])*exp(-T[i]/troe_Tsss[4]) : 0.) 
            + (fabs(troe_Ts[4]) > 1.e-100 ? troe_a[4] * exp(-T[i]/troe_Ts[4]) : 0.) 
            + (troe_len[4] == 4 ? exp(-troe_Tss[4] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[16*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[15*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 6: CH2O (+M) <=> H2 + CO (+M) */
        phi_f = sc[16*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[6*npt+i] + (TB[5][1] - 1)*sc[8*npt+i] + (TB[5][2] - 1)*sc[9*npt+i] + (TB[5][3] - 1)*sc[14*npt+i] + (TB[5][4] - 1)*sc[16*npt+i] + (TB[5][5] - 1)*sc[19*npt+i] + (TB[5][6] - 1)*sc[20*npt+i];
        k_f = k_f_s[5*npt+i];
        redP = alpha / k_f * phase_units[5] * low_A[5] * exp(low_beta[5] * tc[i] - activation_units[5] * low_Ea[5] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[5]) > 1.e-100 ? (1.-troe_a[5])*exp(-T[i]/troe_Tsss[5]) : 0.) 
            + (fabs(troe_Ts[5]) > 1.e-100 ? troe_a[5] * exp(-T[i]/troe_Ts[5]) : 0.) 
            + (troe_len[5] == 4 ? exp(-troe_Tss[5] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[8*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 7: CH3 + H (+M) <=> CH4 (+M) */
        phi_f = sc[2*npt+i]*sc[13*npt+i];
        alpha = mixture[i] + (TB[6][0] - 1)*sc[0*npt+i] + (TB[6][1] - 1)*sc[4*npt+i] + (TB[6][2] - 1)*sc[6*npt+i] + (TB[6][3] - 1)*sc[8*npt+i] + (TB[6][4] - 1)*sc[9*npt+i] + (TB[6][5] - 1)*sc[14*npt+i] + (TB[6][6] - 1)*sc[16*npt+i] + (TB[6][7] - 1)*sc[19*npt+i] + (TB[6][8] - 1)*sc[20*npt+i];
        k_f = k_f_s[6*npt+i];
        redP = alpha / k_f * phase_units[6] * low_A[6] * exp(low_beta[6] * tc[i] - activation_units[6] * low_Ea[6] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[6]) > 1.e-100 ? (1.-troe_a[6])*exp(-T[i]/troe_Tsss[6]) : 0.) 
            + (fabs(troe_Ts[6]) > 1.e-100 ? troe_a[6] * exp(-T[i]/troe_Ts[6]) : 0.) 
            + (troe_len[6] == 4 ? exp(-troe_Tss[6] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 8: CH3 + OH (+M) <=> CH3OH (+M) */
        phi_f = sc[5*npt+i]*sc[13*npt+i];
        alpha = mixture[i] + (TB[7][0] - 1)*sc[6*npt+i] + (TB[7][1] - 1)*sc[8*npt+i] + (TB[7][2] - 1)*sc[9*npt+i] + (TB[7][3] - 1)*sc[14*npt+i] + (TB[7][4] - 1)*sc[16*npt+i] + (TB[7][5] - 1)*sc[19*npt+i];
        k_f = k_f_s[7*npt+i];
        redP = alpha / k_f * phase_units[7] * low_A[7] * exp(low_beta[7] * tc[i] - activation_units[7] * low_Ea[7] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[7]) > 1.e-100 ? (1.-troe_a[7])*exp(-T[i]/troe_Tsss[7]) : 0.) 
            + (fabs(troe_Ts[7]) > 1.e-100 ? troe_a[7] * exp(-T[i]/troe_Ts[7]) : 0.) 
            + (troe_len[7] == 4 ? exp(-troe_Tss[7] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[19*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 9: CH3O (+M) <=> H + CH2O (+M) */
        phi_f = sc[18*npt+i];
        alpha = mixture[i] + (TB[8][0] - 1)*sc[1*npt+i] + (TB[8][1] - 1)*sc[6*npt+i] + (TB[8][2] - 1)*sc[8*npt+i] + (TB[8][3] - 1)*sc[9*npt+i] + (TB[8][4] - 1)*sc[14*npt+i] + (TB[8][5] - 1)*sc[16*npt+i] + (TB[8][6] - 1)*sc[19*npt+i] + (TB[8][7] - 1)*sc[20*npt+i];
        k_f = k_f_s[8*npt+i];
        redP = alpha / k_f * phase_units[8] * low_A[8] * exp(low_beta[8] * tc[i] - activation_units[8] * low_Ea[8] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[8]) > 1.e-100 ? (1.-troe_a[8])*exp(-T[i]/troe_Tsss[8]) : 0.) 
            + (fabs(troe_Ts[8]) > 1.e-100 ? troe_a[8] * exp(-T[i]/troe_Ts[8]) : 0.) 
            + (troe_len[8] == 4 ? exp(-troe_Tss[8] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[16*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 10: CH3O + H (+M) <=> CH3OH (+M) */
        phi_f = sc[2*npt+i]*sc[18*npt+i];
        alpha = mixture[i] + (TB[9][0] - 1)*sc[6*npt+i] + (TB[9][1] - 1)*sc[8*npt+i] + (TB[9][2] - 1)*sc[9*npt+i] + (TB[9][3] - 1)*sc[14*npt+i] + (TB[9][4] - 1)*sc[16*npt+i] + (TB[9][5] - 1)*sc[19*npt+i];
        k_f = k_f_s[9*npt+i];
        redP = alpha / k_f * phase_units[9] * low_A[9] * exp(low_beta[9] * tc[i] - activation_units[9] * low_Ea[9] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[9]) > 1.e-100 ? (1.-troe_a[9])*exp(-T[i]/troe_Tsss[9]) : 0.) 
            + (fabs(troe_Ts[9]) > 1.e-100 ? troe_a[9] * exp(-T[i]/troe_Ts[9]) : 0.) 
            + (troe_len[9] == 4 ? exp(-troe_Tss[9] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[19*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[18*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 11: CH2OH (+M) <=> H + CH2O (+M) */
        phi_f = sc[17*npt+i];
        alpha = mixture[i] + (TB[10][0] - 1)*sc[1*npt+i] + (TB[10][1] - 1)*sc[6*npt+i] + (TB[10][2] - 1)*sc[8*npt+i] + (TB[10][3] - 1)*sc[9*npt+i] + (TB[10][4] - 1)*sc[14*npt+i] + (TB[10][5] - 1)*sc[16*npt+i] + (TB[10][6] - 1)*sc[19*npt+i] + (TB[10][7] - 1)*sc[20*npt+i];
        k_f = k_f_s[10*npt+i];
        redP = alpha / k_f * phase_units[10] * low_A[10] * exp(low_beta[10] * tc[i] - activation_units[10] * low_Ea[10] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[10]) > 1.e-100 ? (1.-troe_a[10])*exp(-T[i]/troe_Tsss[10]) : 0.) 
            + (fabs(troe_Ts[10]) > 1.e-100 ? troe_a[10] * exp(-T[i]/troe_Ts[10]) : 0.) 
            + (troe_len[10] == 4 ? exp(-troe_Tss[10] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[16*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 12: CH2OH + H (+M) <=> CH3OH (+M) */
        phi_f = sc[2*npt+i]*sc[17*npt+i];
        alpha = mixture[i] + (TB[11][0] - 1)*sc[6*npt+i] + (TB[11][1] - 1)*sc[8*npt+i] + (TB[11][2] - 1)*sc[9*npt+i] + (TB[11][3] - 1)*sc[14*npt+i] + (TB[11][4] - 1)*sc[16*npt+i] + (TB[11][5] - 1)*sc[19*npt+i] + (TB[11][6] - 1)*sc[20*npt+i];
        k_f = k_f_s[11*npt+i];
        redP = alpha / k_f * phase_units[11] * low_A[11] * exp(low_beta[11] * tc[i] - activation_units[11] * low_Ea[11] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[11]) > 1.e-100 ? (1.-troe_a[11])*exp(-T[i]/troe_Tsss[11]) : 0.) 
            + (fabs(troe_Ts[11]) > 1.e-100 ? troe_a[11] * exp(-T[i]/troe_Ts[11]) : 0.) 
            + (troe_len[11] == 4 ? exp(-troe_Tss[11] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[19*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[17*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 13: CO + O (+M) <=> CO2 (+M) */
        phi_f = sc[3*npt+i]*sc[8*npt+i];
        alpha = mixture[i] + (TB[12][0] - 1)*sc[1*npt+i] + (TB[12][1] - 1)*sc[6*npt+i] + (TB[12][2] - 1)*sc[8*npt+i] + (TB[12][3] - 1)*sc[9*npt+i] + (TB[12][4] - 1)*sc[14*npt+i] + (TB[12][5] - 1)*sc[16*npt+i] + (TB[12][6] - 1)*sc[19*npt+i];
        k_f = k_f_s[12*npt+i];
        redP = alpha / k_f * phase_units[12] * low_A[12] * exp(low_beta[12] * tc[i] - activation_units[12] * low_Ea[12] * invT[i]);
        F = redP / (1 + redP);
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 14: H2 + M <=> 2.000000 H + M */
        phi_f = sc[1*npt+i];
        alpha = mixture[i] + (TB[13][0] - 1)*sc[0*npt+i] + (TB[13][1] - 1)*sc[1*npt+i] + (TB[13][2] - 1)*sc[6*npt+i] + (TB[13][3] - 1)*sc[8*npt+i] + (TB[13][4] - 1)*sc[9*npt+i] + (TB[13][5] - 1)*sc[14*npt+i] + (TB[13][6] - 1)*sc[16*npt+i] + (TB[13][7] - 1)*sc[19*npt+i] + (TB[13][8] - 1)*sc[20*npt+i];
        k_f = alpha * k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[2*npt+i], 2.000000);
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += 2.000000 * qdot;

        /*reaction 15: 2.000000 O + M <=> O2 + M */
        phi_f = pow(sc[3*npt+i], 2.000000);
        alpha = mixture[i] + (TB[14][0] - 1)*sc[1*npt+i] + (TB[14][1] - 1)*sc[6*npt+i] + (TB[14][2] - 1)*sc[8*npt+i] + (TB[14][3] - 1)*sc[9*npt+i] + (TB[14][4] - 1)*sc[14*npt+i] + (TB[14][5] - 1)*sc[16*npt+i] + (TB[14][6] - 1)*sc[19*npt+i] + (TB[14][7] - 1)*sc[20*npt+i];
        k_f = alpha * k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= 2.000000 * qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 16: O + H + M <=> OH + M */
        phi_f = sc[2*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + (TB[15][0] - 1)*sc[0*npt+i] + (TB[15][1] - 1)*sc[1*npt+i] + (TB[15][2] - 1)*sc[6*npt+i] + (TB[15][3] - 1)*sc[8*npt+i] + (TB[15][4] - 1)*sc[9*npt+i] + (TB[15][5] - 1)*sc[14*npt+i] + (TB[15][6] - 1)*sc[16*npt+i] + (TB[15][7] - 1)*sc[19*npt+i] + (TB[15][8] - 1)*sc[20*npt+i];
        k_f = alpha * k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 17: H2O + M <=> H + OH + M */
        phi_f = sc[6*npt+i];
        alpha = mixture[i] + (TB[16][0] - 1)*sc[0*npt+i] + (TB[16][1] - 1)*sc[1*npt+i] + (TB[16][2] - 1)*sc[4*npt+i] + (TB[16][3] - 1)*sc[6*npt+i] + (TB[16][4] - 1)*sc[8*npt+i] + (TB[16][5] - 1)*sc[9*npt+i] + (TB[16][6] - 1)*sc[14*npt+i] + (TB[16][7] - 1)*sc[16*npt+i] + (TB[16][8] - 1)*sc[19*npt+i] + (TB[16][9] - 1)*sc[20*npt+i];
        k_f = alpha * k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 18: HCO + M <=> H + CO + M */
        phi_f = sc[15*npt+i];
        alpha = mixture[i] + (TB[17][0] - 1)*sc[0*npt+i] + (TB[17][1] - 1)*sc[1*npt+i] + (TB[17][2] - 1)*sc[4*npt+i] + (TB[17][3] - 1)*sc[6*npt+i] + (TB[17][4] - 1)*sc[8*npt+i] + (TB[17][5] - 1)*sc[9*npt+i] + (TB[17][6] - 1)*sc[14*npt+i] + (TB[17][7] - 1)*sc[16*npt+i] + (TB[17][8] - 1)*sc[19*npt+i] + (TB[17][9] - 1)*sc[20*npt+i];
        k_f = alpha * k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[8*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 19: H + O2 <=> O + OH */
        phi_f = sc[2*npt+i]*sc[4*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 20: O + H2 <=> H + OH */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 21: O + H2 <=> H + OH */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 22: OH + H2 <=> H + H2O */
        phi_f = sc[1*npt+i]*sc[5*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[6*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 23: 2.000000 OH <=> O + H2O */
        phi_f = pow(sc[5*npt+i], 2.000000);
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[6*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= 2.000000 * qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 24: H2 + HE <=> 2.000000 H + HE */
        phi_f = sc[1*npt+i]*sc[20*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[2*npt+i], 2.000000)*sc[20*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += 2.000000 * qdot;
        wdot[20*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 25: 2.000000 O + HE <=> O2 + HE */
        phi_f = pow(sc[3*npt+i], 2.000000)*sc[20*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[20*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= 2.000000 * qdot;
        wdot[4*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 26: 2.000000 H2O <=> H + OH + H2O */
        phi_f = pow(sc[6*npt+i], 2.000000);
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[6*npt+i] -= 2.000000 * qdot;

        /*reaction 27: HO2 + H <=> H2 + O2 */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 28: HO2 + H <=> 2.000000 OH */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[5*npt+i], 2.000000);
        Kc = Kc_s[27*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] += 2.000000 * qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 29: HO2 + H <=> O + H2O */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[6*npt+i];
        Kc = Kc_s[28*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 30: HO2 + O <=> OH + O2 */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[29*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[5*npt+i];
        Kc = Kc_s[29*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 31: HO2 + OH <=> H2O + O2 */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[30*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[30*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 32: HO2 + OH <=> H2O + O2 */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[31*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[31*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 33: CO + O2 <=> O + CO2 */
        phi_f = sc[4*npt+i]*sc[8*npt+i];
        k_f = k_f_s[32*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[9*npt+i];
        Kc = Kc_s[32*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 34: CO + OH <=> H + CO2 */
        phi_f = sc[5*npt+i]*sc[8*npt+i];
        k_f = k_f_s[33*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[9*npt+i];
        Kc = Kc_s[33*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 35: CO + OH <=> H + CO2 */
        phi_f = sc[5*npt+i]*sc[8*npt+i];
        k_f = k_f_s[34*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[9*npt+i];
        Kc = Kc_s[34*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 36: CO + HO2 <=> OH + CO2 */
        phi_f = sc[7*npt+i]*sc[8*npt+i];
        k_f = k_f_s[35*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[9*npt+i];
        Kc = Kc_s[35*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 37: HCO + H <=> H2 + CO */
        phi_f = sc[2*npt+i]*sc[15*npt+i];
        k_f = k_f_s[36*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[8*npt+i];
        Kc = Kc_s[36*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 38: HCO + O <=> OH + CO */
        phi_f = sc[3*npt+i]*sc[15*npt+i];
        k_f = k_f_s[37*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[8*npt+i];
        Kc = Kc_s[37*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 39: HCO + O <=> H + CO2 */
        phi_f = sc[3*npt+i]*sc[15*npt+i];
        k_f = k_f_s[38*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[9*npt+i];
        Kc = Kc_s[38*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 40: HCO + OH <=> H2O + CO */
        phi_f = sc[5*npt+i]*sc[15*npt+i];
        k_f = k_f_s[39*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[8*npt+i];
        Kc = Kc_s[39*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 41: HCO + O2 <=> HO2 + CO */
        phi_f = sc[4*npt+i]*sc[15*npt+i];
        k_f = k_f_s[40*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[8*npt+i];
        Kc = Kc_s[40*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 42: CH + O <=> H + CO */
        phi_f = sc[3*npt+i]*sc[10*npt+i];
        k_f = k_f_s[41*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[8*npt+i];
        Kc = Kc_s[41*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 43: CH + OH <=> H + HCO */
        phi_f = sc[5*npt+i]*sc[10*npt+i];
        k_f = k_f_s[42*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[15*npt+i];
        Kc = Kc_s[42*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 44: CH + H2 <=> H + CH2 */
        phi_f = sc[1*npt+i]*sc[10*npt+i];
        k_f = k_f_s[43*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[11*npt+i];
        Kc = Kc_s[43*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 45: CH + H2O <=> H + CH2O */
        phi_f = sc[6*npt+i]*sc[10*npt+i];
        k_f = k_f_s[44*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[16*npt+i];
        Kc = Kc_s[44*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 46: CH + O2 <=> O + HCO */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[45*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[15*npt+i];
        Kc = Kc_s[45*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 47: CH + O2 <=> CO2 + H */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[46*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[9*npt+i];
        Kc = Kc_s[46*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 48: CH + O2 <=> CO + OH */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[47*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[8*npt+i];
        Kc = Kc_s[47*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 49: CH + O2 => O + H + CO */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[48*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 50: CH + CO2 <=> HCO + CO */
        phi_f = sc[9*npt+i]*sc[10*npt+i];
        k_f = k_f_s[49*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[15*npt+i];
        Kc = Kc_s[49*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
    }
}

void vcomp_wdot_51_100(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 51: CH2 + O => 2.000000 H + CO */
        phi_f = sc[3*npt+i]*sc[11*npt+i];
        k_f = k_f_s[50*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += 2.000000 * qdot;
        wdot[3*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 52: CH2 + OH <=> H + CH2O */
        phi_f = sc[5*npt+i]*sc[11*npt+i];
        k_f = k_f_s[51*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[16*npt+i];
        Kc = Kc_s[51*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 53: CH2 + OH <=> CH + H2O */
        phi_f = sc[5*npt+i]*sc[11*npt+i];
        k_f = k_f_s[52*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[10*npt+i];
        Kc = Kc_s[52*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 54: CH2 + HO2 <=> OH + CH2O */
        phi_f = sc[7*npt+i]*sc[11*npt+i];
        k_f = k_f_s[53*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[16*npt+i];
        Kc = Kc_s[53*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 55: CH2 + H2 <=> H + CH3 */
        phi_f = sc[1*npt+i]*sc[11*npt+i];
        k_f = k_f_s[54*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[13*npt+i];
        Kc = Kc_s[54*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 56: CH2 + O2 => OH + H + CO */
        phi_f = sc[4*npt+i]*sc[11*npt+i];
        k_f = k_f_s[55*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 57: CH2 + O2 => 2.000000 H + CO2 */
        phi_f = sc[4*npt+i]*sc[11*npt+i];
        k_f = k_f_s[56*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += 2.000000 * qdot;
        wdot[4*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 58: CH2 + O2 <=> O + CH2O */
        phi_f = sc[4*npt+i]*sc[11*npt+i];
        k_f = k_f_s[57*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[16*npt+i];
        Kc = Kc_s[57*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 59: CH2 + O2 <=> H2 + CO2 */
        phi_f = sc[4*npt+i]*sc[11*npt+i];
        k_f = k_f_s[58*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i];
        Kc = Kc_s[58*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 60: CH2 + O2 <=> H2O + CO */
        phi_f = sc[4*npt+i]*sc[11*npt+i];
        k_f = k_f_s[59*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[8*npt+i];
        Kc = Kc_s[59*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 61: CH2(S) + N2 <=> CH2 + N2 */
        phi_f = sc[0*npt+i]*sc[12*npt+i];
        k_f = k_f_s[60*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[11*npt+i];
        Kc = Kc_s[60*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 62: CH2(S) + HE <=> CH2 + HE */
        phi_f = sc[12*npt+i]*sc[20*npt+i];
        k_f = k_f_s[61*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[20*npt+i];
        Kc = Kc_s[61*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 63: CH2(S) + H <=> CH + H2 */
        phi_f = sc[2*npt+i]*sc[12*npt+i];
        k_f = k_f_s[62*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[10*npt+i];
        Kc = Kc_s[62*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 64: CH2(S) + O => 2.000000 H + CO */
        phi_f = sc[3*npt+i]*sc[12*npt+i];
        k_f = k_f_s[63*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += 2.000000 * qdot;
        wdot[3*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 65: CH2(S) + OH <=> H + CH2O */
        phi_f = sc[5*npt+i]*sc[12*npt+i];
        k_f = k_f_s[64*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[16*npt+i];
        Kc = Kc_s[64*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 66: CH2(S) + H2 <=> CH3 + H */
        phi_f = sc[1*npt+i]*sc[12*npt+i];
        k_f = k_f_s[65*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[13*npt+i];
        Kc = Kc_s[65*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 67: CH2(S) + O2 <=> CH2 + O2 */
        phi_f = sc[4*npt+i]*sc[12*npt+i];
        k_f = k_f_s[66*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[11*npt+i];
        Kc = Kc_s[66*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 68: CH2(S) + H2O <=> CH2 + H2O */
        phi_f = sc[6*npt+i]*sc[12*npt+i];
        k_f = k_f_s[67*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[11*npt+i];
        Kc = Kc_s[67*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 69: CH2(S) + H2O <=> H2 + CH2O */
        phi_f = sc[6*npt+i]*sc[12*npt+i];
        k_f = k_f_s[68*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[16*npt+i];
        Kc = Kc_s[68*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 70: CH2(S) + CO <=> CH2 + CO */
        phi_f = sc[8*npt+i]*sc[12*npt+i];
        k_f = k_f_s[69*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[11*npt+i];
        Kc = Kc_s[69*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
        phi_f = sc[9*npt+i]*sc[12*npt+i];
        k_f = k_f_s[70*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[11*npt+i];
        Kc = Kc_s[70*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
        phi_f = sc[9*npt+i]*sc[12*npt+i];
        k_f = k_f_s[71*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[16*npt+i];
        Kc = Kc_s[71*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 73: CH2O + H <=> HCO + H2 */
        phi_f = sc[2*npt+i]*sc[16*npt+i];
        k_f = k_f_s[72*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[15*npt+i];
        Kc = Kc_s[72*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 74: CH2O + O <=> OH + HCO */
        phi_f = sc[3*npt+i]*sc[16*npt+i];
        k_f = k_f_s[73*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[15*npt+i];
        Kc = Kc_s[73*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 75: CH2O + OH <=> HCO + H2O */
        phi_f = sc[5*npt+i]*sc[16*npt+i];
        k_f = k_f_s[74*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[15*npt+i];
        Kc = Kc_s[74*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 76: CH2O + O2 <=> HO2 + HCO */
        phi_f = sc[4*npt+i]*sc[16*npt+i];
        k_f = k_f_s[75*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[15*npt+i];
        Kc = Kc_s[75*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 77: CH2O + CH2 <=> CH3 + HCO */
        phi_f = sc[11*npt+i]*sc[16*npt+i];
        k_f = k_f_s[76*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[15*npt+i];
        Kc = Kc_s[76*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 78: CH2O + CH2(S) <=> CH3 + HCO */
        phi_f = sc[12*npt+i]*sc[16*npt+i];
        k_f = k_f_s[77*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[15*npt+i];
        Kc = Kc_s[77*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 79: CH3 + O <=> H + CH2O */
        phi_f = sc[3*npt+i]*sc[13*npt+i];
        k_f = k_f_s[78*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[16*npt+i];
        Kc = Kc_s[78*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 80: CH3 + O => H + H2 + CO */
        phi_f = sc[3*npt+i]*sc[13*npt+i];
        k_f = k_f_s[79*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 81: CH3 + OH <=> CH2 + H2O */
        phi_f = sc[5*npt+i]*sc[13*npt+i];
        k_f = k_f_s[80*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[11*npt+i];
        Kc = Kc_s[80*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 82: CH3 + OH <=> CH2(S) + H2O */
        phi_f = sc[5*npt+i]*sc[13*npt+i];
        k_f = k_f_s[81*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[12*npt+i];
        Kc = Kc_s[81*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 83: CH3 + OH <=> H2 + CH2O */
        phi_f = sc[5*npt+i]*sc[13*npt+i];
        k_f = k_f_s[82*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[16*npt+i];
        Kc = Kc_s[82*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 84: CH3 + HO2 <=> O2 + CH4 */
        phi_f = sc[7*npt+i]*sc[13*npt+i];
        k_f = k_f_s[83*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[14*npt+i];
        Kc = Kc_s[83*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 85: CH3 + HO2 <=> OH + CH3O */
        phi_f = sc[7*npt+i]*sc[13*npt+i];
        k_f = k_f_s[84*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[18*npt+i];
        Kc = Kc_s[84*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;

        /*reaction 86: CH3 + O2 <=> O + CH3O */
        phi_f = sc[4*npt+i]*sc[13*npt+i];
        k_f = k_f_s[85*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[18*npt+i];
        Kc = Kc_s[85*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;

        /*reaction 87: CH3 + O2 <=> OH + CH2O */
        phi_f = sc[4*npt+i]*sc[13*npt+i];
        k_f = k_f_s[86*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[16*npt+i];
        Kc = Kc_s[86*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 88: CH3 + HCO <=> CH4 + CO */
        phi_f = sc[13*npt+i]*sc[15*npt+i];
        k_f = k_f_s[87*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[14*npt+i];
        Kc = Kc_s[87*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 89: CH3 + CH2O <=> HCO + CH4 */
        phi_f = sc[13*npt+i]*sc[16*npt+i];
        k_f = k_f_s[88*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[15*npt+i];
        Kc = Kc_s[88*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 90: CH3O + H <=> H + CH2OH */
        phi_f = sc[2*npt+i]*sc[18*npt+i];
        k_f = k_f_s[89*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[17*npt+i];
        Kc = Kc_s[89*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 91: CH3O + H <=> H2 + CH2O */
        phi_f = sc[2*npt+i]*sc[18*npt+i];
        k_f = k_f_s[90*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[16*npt+i];
        Kc = Kc_s[90*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 92: CH3O + H <=> OH + CH3 */
        phi_f = sc[2*npt+i]*sc[18*npt+i];
        k_f = k_f_s[91*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[13*npt+i];
        Kc = Kc_s[91*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 93: CH3O + H <=> CH2(S) + H2O */
        phi_f = sc[2*npt+i]*sc[18*npt+i];
        k_f = k_f_s[92*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[12*npt+i];
        Kc = Kc_s[92*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 94: CH3O + O <=> OH + CH2O */
        phi_f = sc[3*npt+i]*sc[18*npt+i];
        k_f = k_f_s[93*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[16*npt+i];
        Kc = Kc_s[93*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 95: CH3O + OH <=> H2O + CH2O */
        phi_f = sc[5*npt+i]*sc[18*npt+i];
        k_f = k_f_s[94*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[16*npt+i];
        Kc = Kc_s[94*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 96: CH3O + O2 <=> HO2 + CH2O */
        phi_f = sc[4*npt+i]*sc[18*npt+i];
        k_f = k_f_s[95*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[16*npt+i];
        Kc = Kc_s[95*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 97: CH3O + CH3 <=> CH4 + CH2O */
        phi_f = sc[13*npt+i]*sc[18*npt+i];
        k_f = k_f_s[96*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[16*npt+i];
        Kc = Kc_s[96*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 98: CH3O + CO <=> CH3 + CO2 */
        phi_f = sc[8*npt+i]*sc[18*npt+i];
        k_f = k_f_s[97*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[13*npt+i];
        Kc = Kc_s[97*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 99: CH2OH + H <=> H2 + CH2O */
        phi_f = sc[2*npt+i]*sc[17*npt+i];
        k_f = k_f_s[98*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[16*npt+i];
        Kc = Kc_s[98*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 100: CH2OH + H <=> OH + CH3 */
        phi_f = sc[2*npt+i]*sc[17*npt+i];
        k_f = k_f_s[99*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[13*npt+i];
        Kc = Kc_s[99*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;
    }
}

void vcomp_wdot_101_124(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 101: CH2OH + H <=> CH2(S) + H2O */
        phi_f = sc[2*npt+i]*sc[17*npt+i];
        k_f = k_f_s[100*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[12*npt+i];
        Kc = Kc_s[100*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 102: CH2OH + O <=> OH + CH2O */
        phi_f = sc[3*npt+i]*sc[17*npt+i];
        k_f = k_f_s[101*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[16*npt+i];
        Kc = Kc_s[101*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 103: CH2OH + OH <=> H2O + CH2O */
        phi_f = sc[5*npt+i]*sc[17*npt+i];
        k_f = k_f_s[102*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[16*npt+i];
        Kc = Kc_s[102*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 104: CH2OH + O2 <=> HO2 + CH2O */
        phi_f = sc[4*npt+i]*sc[17*npt+i];
        k_f = k_f_s[103*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[16*npt+i];
        Kc = Kc_s[103*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 105: CH2OH + CH3 <=> CH4 + CH2O */
        phi_f = sc[13*npt+i]*sc[17*npt+i];
        k_f = k_f_s[104*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[16*npt+i];
        Kc = Kc_s[104*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 106: CH4 + H <=> CH3 + H2 */
        phi_f = sc[2*npt+i]*sc[14*npt+i];
        k_f = k_f_s[105*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[13*npt+i];
        Kc = Kc_s[105*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 107: CH4 + O <=> OH + CH3 */
        phi_f = sc[3*npt+i]*sc[14*npt+i];
        k_f = k_f_s[106*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[13*npt+i];
        Kc = Kc_s[106*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 108: CH4 + OH <=> CH3 + H2O */
        phi_f = sc[5*npt+i]*sc[14*npt+i];
        k_f = k_f_s[107*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[13*npt+i];
        Kc = Kc_s[107*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 109: CH4 + CH2 <=> 2.000000 CH3 */
        phi_f = sc[11*npt+i]*sc[14*npt+i];
        k_f = k_f_s[108*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[13*npt+i], 2.000000);
        Kc = Kc_s[108*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[13*npt+i] += 2.000000 * qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 110: CH4 + CH2(S) <=> 2.000000 CH3 */
        phi_f = sc[12*npt+i]*sc[14*npt+i];
        k_f = k_f_s[109*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[13*npt+i], 2.000000);
        Kc = Kc_s[109*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += 2.000000 * qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 111: CH3OH + H <=> CH2OH + H2 */
        phi_f = sc[2*npt+i]*sc[19*npt+i];
        k_f = k_f_s[110*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[17*npt+i];
        Kc = Kc_s[110*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 112: CH3OH + H <=> CH3O + H2 */
        phi_f = sc[2*npt+i]*sc[19*npt+i];
        k_f = k_f_s[111*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[18*npt+i];
        Kc = Kc_s[111*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 113: CH3OH + O <=> OH + CH2OH */
        phi_f = sc[3*npt+i]*sc[19*npt+i];
        k_f = k_f_s[112*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[17*npt+i];
        Kc = Kc_s[112*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 114: CH3OH + O <=> OH + CH3O */
        phi_f = sc[3*npt+i]*sc[19*npt+i];
        k_f = k_f_s[113*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[18*npt+i];
        Kc = Kc_s[113*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 115: CH3OH + OH <=> CH2OH + H2O */
        phi_f = sc[5*npt+i]*sc[19*npt+i];
        k_f = k_f_s[114*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[17*npt+i];
        Kc = Kc_s[114*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 116: CH3OH + OH <=> CH3O + H2O */
        phi_f = sc[5*npt+i]*sc[19*npt+i];
        k_f = k_f_s[115*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[18*npt+i];
        Kc = Kc_s[115*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 117: CH3OH + O2 <=> CH2OH + HO2 */
        phi_f = sc[4*npt+i]*sc[19*npt+i];
        k_f = k_f_s[116*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[17*npt+i];
        Kc = Kc_s[116*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 118: CH3OH + CH <=> CH3 + CH2O */
        phi_f = sc[10*npt+i]*sc[19*npt+i];
        k_f = k_f_s[117*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[16*npt+i];
        Kc = Kc_s[117*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 119: CH3OH + CH2 <=> CH3 + CH2OH */
        phi_f = sc[11*npt+i]*sc[19*npt+i];
        k_f = k_f_s[118*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[17*npt+i];
        Kc = Kc_s[118*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 120: CH3OH + CH2 <=> CH3 + CH3O */
        phi_f = sc[11*npt+i]*sc[19*npt+i];
        k_f = k_f_s[119*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[18*npt+i];
        Kc = Kc_s[119*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 121: CH3OH + CH2(S) <=> CH3 + CH3O */
        phi_f = sc[12*npt+i]*sc[19*npt+i];
        k_f = k_f_s[120*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[18*npt+i];
        Kc = Kc_s[120*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 122: CH3OH + CH2(S) <=> CH3 + CH2OH */
        phi_f = sc[12*npt+i]*sc[19*npt+i];
        k_f = k_f_s[121*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[17*npt+i];
        Kc = Kc_s[121*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 123: CH3OH + CH3 <=> CH2OH + CH4 */
        phi_f = sc[13*npt+i]*sc[19*npt+i];
        k_f = k_f_s[122*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[17*npt+i];
        Kc = Kc_s[122*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 124: CH3OH + CH3 <=> CH3O + CH4 */
        phi_f = sc[13*npt+i]*sc[19*npt+i];
        k_f = k_f_s[123*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[18*npt+i];
        Kc = Kc_s[123*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;
    }
}


static amrex::Real T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static amrex::Real k_f_save[124];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static amrex::Real Kc_save[124];
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

    amrex::Real qdot, q_f[124], q_r[124];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 21; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[2] -= qdot;
    wdot[4] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[1] -= qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[2] -= qdot;
    wdot[11] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[6] -= qdot;
    wdot[12] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[2] -= qdot;
    wdot[15] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[1] += qdot;
    wdot[8] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[2] -= qdot;
    wdot[13] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[5] -= qdot;
    wdot[13] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[2] += qdot;
    wdot[16] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[9]-q_r[9];
    wdot[2] -= qdot;
    wdot[18] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[2] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[11]-q_r[11];
    wdot[2] -= qdot;
    wdot[17] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[3] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[1] -= qdot;
    wdot[2] += 2.000000 * qdot;

    qdot = q_f[14]-q_r[14];
    wdot[3] -= 2.000000 * qdot;
    wdot[4] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[2] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[2] += qdot;
    wdot[8] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[19]-q_r[19];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[20]-q_r[20];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[3] += qdot;
    wdot[5] -= 2.000000 * qdot;
    wdot[6] += qdot;

    qdot = q_f[23]-q_r[23];
    wdot[1] -= qdot;
    wdot[2] += 2.000000 * qdot;
    wdot[20] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[24]-q_r[24];
    wdot[3] -= 2.000000 * qdot;
    wdot[4] += qdot;
    wdot[20] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[25]-q_r[25];
    wdot[2] += qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[6] -= 2.000000 * qdot;

    qdot = q_f[26]-q_r[26];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[27]-q_r[27];
    wdot[2] -= qdot;
    wdot[5] += 2.000000 * qdot;
    wdot[7] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[29]-q_r[29];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[30]-q_r[30];
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[31]-q_r[31];
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[32]-q_r[32];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[34]-q_r[34];
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[5] += qdot;
    wdot[7] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[36]-q_r[36];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[8] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[37]-q_r[37];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[8] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[38]-q_r[38];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[9] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[39]-q_r[39];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[8] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[40]-q_r[40];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[8] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[41]-q_r[41];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[42]-q_r[42];
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[10] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[43]-q_r[43];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[10] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[44]-q_r[44];
    wdot[2] += qdot;
    wdot[6] -= qdot;
    wdot[10] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[45]-q_r[45];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[10] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[46]-q_r[46];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[47]-q_r[47];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[48]-q_r[48];
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[49]-q_r[49];
    wdot[8] += qdot;
    wdot[9] -= qdot;
    wdot[10] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[50]-q_r[50];
    wdot[2] += 2.000000 * qdot;
    wdot[3] -= qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[51]-q_r[51];
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[11] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[52]-q_r[52];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[53]-q_r[53];
    wdot[5] += qdot;
    wdot[7] -= qdot;
    wdot[11] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[54]-q_r[54];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[11] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[55]-q_r[55];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[56]-q_r[56];
    wdot[2] += 2.000000 * qdot;
    wdot[4] -= qdot;
    wdot[9] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[57]-q_r[57];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[11] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[58]-q_r[58];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[9] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[59]-q_r[59];
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[60]-q_r[60];
    wdot[0] -= qdot;
    wdot[0] += qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[61]-q_r[61];
    wdot[11] += qdot;
    wdot[12] -= qdot;
    wdot[20] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[62]-q_r[62];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[10] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[63]-q_r[63];
    wdot[2] += 2.000000 * qdot;
    wdot[3] -= qdot;
    wdot[8] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[64]-q_r[64];
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[12] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[65]-q_r[65];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[66]-q_r[66];
    wdot[4] -= qdot;
    wdot[4] += qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[67]-q_r[67];
    wdot[6] -= qdot;
    wdot[6] += qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[68]-q_r[68];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[12] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[69]-q_r[69];
    wdot[8] -= qdot;
    wdot[8] += qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[70]-q_r[70];
    wdot[9] -= qdot;
    wdot[9] += qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[71]-q_r[71];
    wdot[8] += qdot;
    wdot[9] -= qdot;
    wdot[12] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[72]-q_r[72];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[73]-q_r[73];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[74]-q_r[74];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[75]-q_r[75];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[76]-q_r[76];
    wdot[11] -= qdot;
    wdot[13] += qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[77]-q_r[77];
    wdot[12] -= qdot;
    wdot[13] += qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[78]-q_r[78];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[13] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[79]-q_r[79];
    wdot[1] += qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[8] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[80]-q_r[80];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[81]-q_r[81];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[12] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[82]-q_r[82];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[13] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[83]-q_r[83];
    wdot[4] += qdot;
    wdot[7] -= qdot;
    wdot[13] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[84]-q_r[84];
    wdot[5] += qdot;
    wdot[7] -= qdot;
    wdot[13] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[85]-q_r[85];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[13] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[86]-q_r[86];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[13] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[87]-q_r[87];
    wdot[8] += qdot;
    wdot[13] -= qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[88]-q_r[88];
    wdot[13] -= qdot;
    wdot[14] += qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[89]-q_r[89];
    wdot[2] -= qdot;
    wdot[2] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[90]-q_r[90];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[16] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[91]-q_r[91];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[13] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[92]-q_r[92];
    wdot[2] -= qdot;
    wdot[6] += qdot;
    wdot[12] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[93]-q_r[93];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[16] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[94]-q_r[94];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[16] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[95]-q_r[95];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[16] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[96]-q_r[96];
    wdot[13] -= qdot;
    wdot[14] += qdot;
    wdot[16] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[97]-q_r[97];
    wdot[8] -= qdot;
    wdot[9] += qdot;
    wdot[13] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[98]-q_r[98];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[99]-q_r[99];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[13] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[100]-q_r[100];
    wdot[2] -= qdot;
    wdot[6] += qdot;
    wdot[12] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[101]-q_r[101];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[102]-q_r[102];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[103]-q_r[103];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[104]-q_r[104];
    wdot[13] -= qdot;
    wdot[14] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[105]-q_r[105];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[106]-q_r[106];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[107]-q_r[107];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[108]-q_r[108];
    wdot[11] -= qdot;
    wdot[13] += 2.000000 * qdot;
    wdot[14] -= qdot;

    qdot = q_f[109]-q_r[109];
    wdot[12] -= qdot;
    wdot[13] += 2.000000 * qdot;
    wdot[14] -= qdot;

    qdot = q_f[110]-q_r[110];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[17] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[111]-q_r[111];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[112]-q_r[112];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[17] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[113]-q_r[113];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[114]-q_r[114];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[17] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[115]-q_r[115];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[116]-q_r[116];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[17] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[117]-q_r[117];
    wdot[10] -= qdot;
    wdot[13] += qdot;
    wdot[16] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[118]-q_r[118];
    wdot[11] -= qdot;
    wdot[13] += qdot;
    wdot[17] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[119]-q_r[119];
    wdot[11] -= qdot;
    wdot[13] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[120]-q_r[120];
    wdot[12] -= qdot;
    wdot[13] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[121]-q_r[121];
    wdot[12] -= qdot;
    wdot[13] += qdot;
    wdot[17] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[122]-q_r[122];
    wdot[13] -= qdot;
    wdot[14] += qdot;
    wdot[17] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[123]-q_r[123];
    wdot[13] -= qdot;
    wdot[14] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    return;
}

void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f)
{
    for (int i=0; i<124; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(amrex::Real *  tc, amrex::Real invT, amrex::Real *  Kc)
{
    /*compute the Gibbs free energy */
    amrex::Real g_RT[21];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[2] + g_RT[4] - g_RT[7];
    Kc[1] = g_RT[1] + g_RT[10] - g_RT[13];
    Kc[2] = g_RT[2] + g_RT[11] - g_RT[13];
    Kc[3] = g_RT[6] + g_RT[12] - g_RT[19];
    Kc[4] = g_RT[2] + g_RT[15] - g_RT[16];
    Kc[5] = -g_RT[1] - g_RT[8] + g_RT[16];
    Kc[6] = g_RT[2] + g_RT[13] - g_RT[14];
    Kc[7] = g_RT[5] + g_RT[13] - g_RT[19];
    Kc[8] = -g_RT[2] - g_RT[16] + g_RT[18];
    Kc[9] = g_RT[2] + g_RT[18] - g_RT[19];
    Kc[10] = -g_RT[2] - g_RT[16] + g_RT[17];
    Kc[11] = g_RT[2] + g_RT[17] - g_RT[19];
    Kc[12] = g_RT[3] + g_RT[8] - g_RT[9];
    Kc[13] = g_RT[1] - 2.000000*g_RT[2];
    Kc[14] = 2.000000*g_RT[3] - g_RT[4];
    Kc[15] = g_RT[2] + g_RT[3] - g_RT[5];
    Kc[16] = -g_RT[2] - g_RT[5] + g_RT[6];
    Kc[17] = -g_RT[2] - g_RT[8] + g_RT[15];
    Kc[18] = g_RT[2] - g_RT[3] + g_RT[4] - g_RT[5];
    Kc[19] = g_RT[1] - g_RT[2] + g_RT[3] - g_RT[5];
    Kc[20] = g_RT[1] - g_RT[2] + g_RT[3] - g_RT[5];
    Kc[21] = g_RT[1] - g_RT[2] + g_RT[5] - g_RT[6];
    Kc[22] = -g_RT[3] + 2.000000*g_RT[5] - g_RT[6];
    Kc[23] = g_RT[1] - 2.000000*g_RT[2] + g_RT[20] - g_RT[20];
    Kc[24] = 2.000000*g_RT[3] - g_RT[4] + g_RT[20] - g_RT[20];
    Kc[25] = -g_RT[2] - g_RT[5] + 2.000000*g_RT[6] - g_RT[6];
    Kc[26] = -g_RT[1] + g_RT[2] - g_RT[4] + g_RT[7];
    Kc[27] = g_RT[2] - 2.000000*g_RT[5] + g_RT[7];
    Kc[28] = g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7];
    Kc[29] = g_RT[3] - g_RT[4] - g_RT[5] + g_RT[7];
    Kc[30] = -g_RT[4] + g_RT[5] - g_RT[6] + g_RT[7];
    Kc[31] = -g_RT[4] + g_RT[5] - g_RT[6] + g_RT[7];
    Kc[32] = -g_RT[3] + g_RT[4] + g_RT[8] - g_RT[9];
    Kc[33] = -g_RT[2] + g_RT[5] + g_RT[8] - g_RT[9];
    Kc[34] = -g_RT[2] + g_RT[5] + g_RT[8] - g_RT[9];
    Kc[35] = -g_RT[5] + g_RT[7] + g_RT[8] - g_RT[9];
    Kc[36] = -g_RT[1] + g_RT[2] - g_RT[8] + g_RT[15];
    Kc[37] = g_RT[3] - g_RT[5] - g_RT[8] + g_RT[15];
    Kc[38] = -g_RT[2] + g_RT[3] - g_RT[9] + g_RT[15];
    Kc[39] = g_RT[5] - g_RT[6] - g_RT[8] + g_RT[15];
    Kc[40] = g_RT[4] - g_RT[7] - g_RT[8] + g_RT[15];
    Kc[41] = -g_RT[2] + g_RT[3] - g_RT[8] + g_RT[10];
    Kc[42] = -g_RT[2] + g_RT[5] + g_RT[10] - g_RT[15];
    Kc[43] = g_RT[1] - g_RT[2] + g_RT[10] - g_RT[11];
    Kc[44] = -g_RT[2] + g_RT[6] + g_RT[10] - g_RT[16];
    Kc[45] = -g_RT[3] + g_RT[4] + g_RT[10] - g_RT[15];
    Kc[46] = -g_RT[2] + g_RT[4] - g_RT[9] + g_RT[10];
    Kc[47] = g_RT[4] - g_RT[5] - g_RT[8] + g_RT[10];
    Kc[48] = -g_RT[2] - g_RT[3] + g_RT[4] - g_RT[8] + g_RT[10];
    Kc[49] = -g_RT[8] + g_RT[9] + g_RT[10] - g_RT[15];
    Kc[50] = -2.000000*g_RT[2] + g_RT[3] - g_RT[8] + g_RT[11];
    Kc[51] = -g_RT[2] + g_RT[5] + g_RT[11] - g_RT[16];
    Kc[52] = g_RT[5] - g_RT[6] - g_RT[10] + g_RT[11];
    Kc[53] = -g_RT[5] + g_RT[7] + g_RT[11] - g_RT[16];
    Kc[54] = g_RT[1] - g_RT[2] + g_RT[11] - g_RT[13];
    Kc[55] = -g_RT[2] + g_RT[4] - g_RT[5] - g_RT[8] + g_RT[11];
    Kc[56] = -2.000000*g_RT[2] + g_RT[4] - g_RT[9] + g_RT[11];
    Kc[57] = -g_RT[3] + g_RT[4] + g_RT[11] - g_RT[16];
    Kc[58] = -g_RT[1] + g_RT[4] - g_RT[9] + g_RT[11];
    Kc[59] = g_RT[4] - g_RT[6] - g_RT[8] + g_RT[11];
    Kc[60] = g_RT[0] - g_RT[0] - g_RT[11] + g_RT[12];
    Kc[61] = -g_RT[11] + g_RT[12] + g_RT[20] - g_RT[20];
    Kc[62] = -g_RT[1] + g_RT[2] - g_RT[10] + g_RT[12];
    Kc[63] = -2.000000*g_RT[2] + g_RT[3] - g_RT[8] + g_RT[12];
    Kc[64] = -g_RT[2] + g_RT[5] + g_RT[12] - g_RT[16];
    Kc[65] = g_RT[1] - g_RT[2] + g_RT[12] - g_RT[13];
    Kc[66] = g_RT[4] - g_RT[4] - g_RT[11] + g_RT[12];
    Kc[67] = g_RT[6] - g_RT[6] - g_RT[11] + g_RT[12];
    Kc[68] = -g_RT[1] + g_RT[6] + g_RT[12] - g_RT[16];
    Kc[69] = g_RT[8] - g_RT[8] - g_RT[11] + g_RT[12];
    Kc[70] = g_RT[9] - g_RT[9] - g_RT[11] + g_RT[12];
    Kc[71] = -g_RT[8] + g_RT[9] + g_RT[12] - g_RT[16];
    Kc[72] = -g_RT[1] + g_RT[2] - g_RT[15] + g_RT[16];
    Kc[73] = g_RT[3] - g_RT[5] - g_RT[15] + g_RT[16];
    Kc[74] = g_RT[5] - g_RT[6] - g_RT[15] + g_RT[16];
    Kc[75] = g_RT[4] - g_RT[7] - g_RT[15] + g_RT[16];
    Kc[76] = g_RT[11] - g_RT[13] - g_RT[15] + g_RT[16];
    Kc[77] = g_RT[12] - g_RT[13] - g_RT[15] + g_RT[16];
    Kc[78] = -g_RT[2] + g_RT[3] + g_RT[13] - g_RT[16];
    Kc[79] = -g_RT[1] - g_RT[2] + g_RT[3] - g_RT[8] + g_RT[13];
    Kc[80] = g_RT[5] - g_RT[6] - g_RT[11] + g_RT[13];
    Kc[81] = g_RT[5] - g_RT[6] - g_RT[12] + g_RT[13];
    Kc[82] = -g_RT[1] + g_RT[5] + g_RT[13] - g_RT[16];
    Kc[83] = -g_RT[4] + g_RT[7] + g_RT[13] - g_RT[14];
    Kc[84] = -g_RT[5] + g_RT[7] + g_RT[13] - g_RT[18];
    Kc[85] = -g_RT[3] + g_RT[4] + g_RT[13] - g_RT[18];
    Kc[86] = g_RT[4] - g_RT[5] + g_RT[13] - g_RT[16];
    Kc[87] = -g_RT[8] + g_RT[13] - g_RT[14] + g_RT[15];
    Kc[88] = g_RT[13] - g_RT[14] - g_RT[15] + g_RT[16];
    Kc[89] = g_RT[2] - g_RT[2] - g_RT[17] + g_RT[18];
    Kc[90] = -g_RT[1] + g_RT[2] - g_RT[16] + g_RT[18];
    Kc[91] = g_RT[2] - g_RT[5] - g_RT[13] + g_RT[18];
    Kc[92] = g_RT[2] - g_RT[6] - g_RT[12] + g_RT[18];
    Kc[93] = g_RT[3] - g_RT[5] - g_RT[16] + g_RT[18];
    Kc[94] = g_RT[5] - g_RT[6] - g_RT[16] + g_RT[18];
    Kc[95] = g_RT[4] - g_RT[7] - g_RT[16] + g_RT[18];
    Kc[96] = g_RT[13] - g_RT[14] - g_RT[16] + g_RT[18];
    Kc[97] = g_RT[8] - g_RT[9] - g_RT[13] + g_RT[18];
    Kc[98] = -g_RT[1] + g_RT[2] - g_RT[16] + g_RT[17];
    Kc[99] = g_RT[2] - g_RT[5] - g_RT[13] + g_RT[17];
    Kc[100] = g_RT[2] - g_RT[6] - g_RT[12] + g_RT[17];
    Kc[101] = g_RT[3] - g_RT[5] - g_RT[16] + g_RT[17];
    Kc[102] = g_RT[5] - g_RT[6] - g_RT[16] + g_RT[17];
    Kc[103] = g_RT[4] - g_RT[7] - g_RT[16] + g_RT[17];
    Kc[104] = g_RT[13] - g_RT[14] - g_RT[16] + g_RT[17];
    Kc[105] = -g_RT[1] + g_RT[2] - g_RT[13] + g_RT[14];
    Kc[106] = g_RT[3] - g_RT[5] - g_RT[13] + g_RT[14];
    Kc[107] = g_RT[5] - g_RT[6] - g_RT[13] + g_RT[14];
    Kc[108] = g_RT[11] - 2.000000*g_RT[13] + g_RT[14];
    Kc[109] = g_RT[12] - 2.000000*g_RT[13] + g_RT[14];
    Kc[110] = -g_RT[1] + g_RT[2] - g_RT[17] + g_RT[19];
    Kc[111] = -g_RT[1] + g_RT[2] - g_RT[18] + g_RT[19];
    Kc[112] = g_RT[3] - g_RT[5] - g_RT[17] + g_RT[19];
    Kc[113] = g_RT[3] - g_RT[5] - g_RT[18] + g_RT[19];
    Kc[114] = g_RT[5] - g_RT[6] - g_RT[17] + g_RT[19];
    Kc[115] = g_RT[5] - g_RT[6] - g_RT[18] + g_RT[19];
    Kc[116] = g_RT[4] - g_RT[7] - g_RT[17] + g_RT[19];
    Kc[117] = g_RT[10] - g_RT[13] - g_RT[16] + g_RT[19];
    Kc[118] = g_RT[11] - g_RT[13] - g_RT[17] + g_RT[19];
    Kc[119] = g_RT[11] - g_RT[13] - g_RT[18] + g_RT[19];
    Kc[120] = g_RT[12] - g_RT[13] - g_RT[18] + g_RT[19];
    Kc[121] = g_RT[12] - g_RT[13] - g_RT[17] + g_RT[19];
    Kc[122] = g_RT[13] - g_RT[14] - g_RT[17] + g_RT[19];
    Kc[123] = g_RT[13] - g_RT[14] - g_RT[18] + g_RT[19];

    for (int i=0; i<124; ++i) {
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
    Kc[5] *= refC;
    Kc[6] *= refCinv;
    Kc[7] *= refCinv;
    Kc[8] *= refC;
    Kc[9] *= refCinv;
    Kc[10] *= refC;
    Kc[11] *= refCinv;
    Kc[12] *= refCinv;
    Kc[13] *= refC;
    Kc[14] *= refCinv;
    Kc[15] *= refCinv;
    Kc[16] *= refC;
    Kc[17] *= refC;
    Kc[23] *= refC;
    Kc[24] *= refCinv;
    Kc[25] *= refC;
    Kc[48] *= refC;
    Kc[50] *= refC;
    Kc[55] *= refC;
    Kc[56] *= refC;
    Kc[63] *= refC;
    Kc[79] *= refC;

    return;
}

void comp_qfqr(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT)
{

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[2]*sc[4];
    qr[0] = sc[7];

    /*reaction 2: CH + H2 (+M) <=> CH3 (+M) */
    qf[1] = sc[1]*sc[10];
    qr[1] = sc[13];

    /*reaction 3: CH2 + H (+M) <=> CH3 (+M) */
    qf[2] = sc[2]*sc[11];
    qr[2] = sc[13];

    /*reaction 4: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    qf[3] = sc[6]*sc[12];
    qr[3] = sc[19];

    /*reaction 5: HCO + H (+M) <=> CH2O (+M) */
    qf[4] = sc[2]*sc[15];
    qr[4] = sc[16];

    /*reaction 6: CH2O (+M) <=> H2 + CO (+M) */
    qf[5] = sc[16];
    qr[5] = sc[1]*sc[8];

    /*reaction 7: CH3 + H (+M) <=> CH4 (+M) */
    qf[6] = sc[2]*sc[13];
    qr[6] = sc[14];

    /*reaction 8: CH3 + OH (+M) <=> CH3OH (+M) */
    qf[7] = sc[5]*sc[13];
    qr[7] = sc[19];

    /*reaction 9: CH3O (+M) <=> H + CH2O (+M) */
    qf[8] = sc[18];
    qr[8] = sc[2]*sc[16];

    /*reaction 10: CH3O + H (+M) <=> CH3OH (+M) */
    qf[9] = sc[2]*sc[18];
    qr[9] = sc[19];

    /*reaction 11: CH2OH (+M) <=> H + CH2O (+M) */
    qf[10] = sc[17];
    qr[10] = sc[2]*sc[16];

    /*reaction 12: CH2OH + H (+M) <=> CH3OH (+M) */
    qf[11] = sc[2]*sc[17];
    qr[11] = sc[19];

    /*reaction 13: CO + O (+M) <=> CO2 (+M) */
    qf[12] = sc[3]*sc[8];
    qr[12] = sc[9];

    /*reaction 14: H2 + M <=> 2.000000 H + M */
    qf[13] = sc[1];
    qr[13] = pow(sc[2], 2.000000);

    /*reaction 15: 2.000000 O + M <=> O2 + M */
    qf[14] = pow(sc[3], 2.000000);
    qr[14] = sc[4];

    /*reaction 16: O + H + M <=> OH + M */
    qf[15] = sc[2]*sc[3];
    qr[15] = sc[5];

    /*reaction 17: H2O + M <=> H + OH + M */
    qf[16] = sc[6];
    qr[16] = sc[2]*sc[5];

    /*reaction 18: HCO + M <=> H + CO + M */
    qf[17] = sc[15];
    qr[17] = sc[2]*sc[8];

    /*reaction 19: H + O2 <=> O + OH */
    qf[18] = sc[2]*sc[4];
    qr[18] = sc[3]*sc[5];

    /*reaction 20: O + H2 <=> H + OH */
    qf[19] = sc[1]*sc[3];
    qr[19] = sc[2]*sc[5];

    /*reaction 21: O + H2 <=> H + OH */
    qf[20] = sc[1]*sc[3];
    qr[20] = sc[2]*sc[5];

    /*reaction 22: OH + H2 <=> H + H2O */
    qf[21] = sc[1]*sc[5];
    qr[21] = sc[2]*sc[6];

    /*reaction 23: 2.000000 OH <=> O + H2O */
    qf[22] = pow(sc[5], 2.000000);
    qr[22] = sc[3]*sc[6];

    /*reaction 24: H2 + HE <=> 2.000000 H + HE */
    qf[23] = sc[1]*sc[20];
    qr[23] = pow(sc[2], 2.000000)*sc[20];

    /*reaction 25: 2.000000 O + HE <=> O2 + HE */
    qf[24] = pow(sc[3], 2.000000)*sc[20];
    qr[24] = sc[4]*sc[20];

    /*reaction 26: 2.000000 H2O <=> H + OH + H2O */
    qf[25] = pow(sc[6], 2.000000);
    qr[25] = sc[2]*sc[5]*sc[6];

    /*reaction 27: HO2 + H <=> H2 + O2 */
    qf[26] = sc[2]*sc[7];
    qr[26] = sc[1]*sc[4];

    /*reaction 28: HO2 + H <=> 2.000000 OH */
    qf[27] = sc[2]*sc[7];
    qr[27] = pow(sc[5], 2.000000);

    /*reaction 29: HO2 + H <=> O + H2O */
    qf[28] = sc[2]*sc[7];
    qr[28] = sc[3]*sc[6];

    /*reaction 30: HO2 + O <=> OH + O2 */
    qf[29] = sc[3]*sc[7];
    qr[29] = sc[4]*sc[5];

    /*reaction 31: HO2 + OH <=> H2O + O2 */
    qf[30] = sc[5]*sc[7];
    qr[30] = sc[4]*sc[6];

    /*reaction 32: HO2 + OH <=> H2O + O2 */
    qf[31] = sc[5]*sc[7];
    qr[31] = sc[4]*sc[6];

    /*reaction 33: CO + O2 <=> O + CO2 */
    qf[32] = sc[4]*sc[8];
    qr[32] = sc[3]*sc[9];

    /*reaction 34: CO + OH <=> H + CO2 */
    qf[33] = sc[5]*sc[8];
    qr[33] = sc[2]*sc[9];

    /*reaction 35: CO + OH <=> H + CO2 */
    qf[34] = sc[5]*sc[8];
    qr[34] = sc[2]*sc[9];

    /*reaction 36: CO + HO2 <=> OH + CO2 */
    qf[35] = sc[7]*sc[8];
    qr[35] = sc[5]*sc[9];

    /*reaction 37: HCO + H <=> H2 + CO */
    qf[36] = sc[2]*sc[15];
    qr[36] = sc[1]*sc[8];

    /*reaction 38: HCO + O <=> OH + CO */
    qf[37] = sc[3]*sc[15];
    qr[37] = sc[5]*sc[8];

    /*reaction 39: HCO + O <=> H + CO2 */
    qf[38] = sc[3]*sc[15];
    qr[38] = sc[2]*sc[9];

    /*reaction 40: HCO + OH <=> H2O + CO */
    qf[39] = sc[5]*sc[15];
    qr[39] = sc[6]*sc[8];

    /*reaction 41: HCO + O2 <=> HO2 + CO */
    qf[40] = sc[4]*sc[15];
    qr[40] = sc[7]*sc[8];

    /*reaction 42: CH + O <=> H + CO */
    qf[41] = sc[3]*sc[10];
    qr[41] = sc[2]*sc[8];

    /*reaction 43: CH + OH <=> H + HCO */
    qf[42] = sc[5]*sc[10];
    qr[42] = sc[2]*sc[15];

    /*reaction 44: CH + H2 <=> H + CH2 */
    qf[43] = sc[1]*sc[10];
    qr[43] = sc[2]*sc[11];

    /*reaction 45: CH + H2O <=> H + CH2O */
    qf[44] = sc[6]*sc[10];
    qr[44] = sc[2]*sc[16];

    /*reaction 46: CH + O2 <=> O + HCO */
    qf[45] = sc[4]*sc[10];
    qr[45] = sc[3]*sc[15];

    /*reaction 47: CH + O2 <=> CO2 + H */
    qf[46] = sc[4]*sc[10];
    qr[46] = sc[2]*sc[9];

    /*reaction 48: CH + O2 <=> CO + OH */
    qf[47] = sc[4]*sc[10];
    qr[47] = sc[5]*sc[8];

    /*reaction 49: CH + O2 => O + H + CO */
    qf[48] = sc[4]*sc[10];
    qr[48] = 0.0;

    /*reaction 50: CH + CO2 <=> HCO + CO */
    qf[49] = sc[9]*sc[10];
    qr[49] = sc[8]*sc[15];

    /*reaction 51: CH2 + O => 2.000000 H + CO */
    qf[50] = sc[3]*sc[11];
    qr[50] = 0.0;

    /*reaction 52: CH2 + OH <=> H + CH2O */
    qf[51] = sc[5]*sc[11];
    qr[51] = sc[2]*sc[16];

    /*reaction 53: CH2 + OH <=> CH + H2O */
    qf[52] = sc[5]*sc[11];
    qr[52] = sc[6]*sc[10];

    /*reaction 54: CH2 + HO2 <=> OH + CH2O */
    qf[53] = sc[7]*sc[11];
    qr[53] = sc[5]*sc[16];

    /*reaction 55: CH2 + H2 <=> H + CH3 */
    qf[54] = sc[1]*sc[11];
    qr[54] = sc[2]*sc[13];

    /*reaction 56: CH2 + O2 => OH + H + CO */
    qf[55] = sc[4]*sc[11];
    qr[55] = 0.0;

    /*reaction 57: CH2 + O2 => 2.000000 H + CO2 */
    qf[56] = sc[4]*sc[11];
    qr[56] = 0.0;

    /*reaction 58: CH2 + O2 <=> O + CH2O */
    qf[57] = sc[4]*sc[11];
    qr[57] = sc[3]*sc[16];

    /*reaction 59: CH2 + O2 <=> H2 + CO2 */
    qf[58] = sc[4]*sc[11];
    qr[58] = sc[1]*sc[9];

    /*reaction 60: CH2 + O2 <=> H2O + CO */
    qf[59] = sc[4]*sc[11];
    qr[59] = sc[6]*sc[8];

    /*reaction 61: CH2(S) + N2 <=> CH2 + N2 */
    qf[60] = sc[0]*sc[12];
    qr[60] = sc[0]*sc[11];

    /*reaction 62: CH2(S) + HE <=> CH2 + HE */
    qf[61] = sc[12]*sc[20];
    qr[61] = sc[11]*sc[20];

    /*reaction 63: CH2(S) + H <=> CH + H2 */
    qf[62] = sc[2]*sc[12];
    qr[62] = sc[1]*sc[10];

    /*reaction 64: CH2(S) + O => 2.000000 H + CO */
    qf[63] = sc[3]*sc[12];
    qr[63] = 0.0;

    /*reaction 65: CH2(S) + OH <=> H + CH2O */
    qf[64] = sc[5]*sc[12];
    qr[64] = sc[2]*sc[16];

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    qf[65] = sc[1]*sc[12];
    qr[65] = sc[2]*sc[13];

    /*reaction 67: CH2(S) + O2 <=> CH2 + O2 */
    qf[66] = sc[4]*sc[12];
    qr[66] = sc[4]*sc[11];

    /*reaction 68: CH2(S) + H2O <=> CH2 + H2O */
    qf[67] = sc[6]*sc[12];
    qr[67] = sc[6]*sc[11];

    /*reaction 69: CH2(S) + H2O <=> H2 + CH2O */
    qf[68] = sc[6]*sc[12];
    qr[68] = sc[1]*sc[16];

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    qf[69] = sc[8]*sc[12];
    qr[69] = sc[8]*sc[11];

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    qf[70] = sc[9]*sc[12];
    qr[70] = sc[9]*sc[11];

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    qf[71] = sc[9]*sc[12];
    qr[71] = sc[8]*sc[16];

    /*reaction 73: CH2O + H <=> HCO + H2 */
    qf[72] = sc[2]*sc[16];
    qr[72] = sc[1]*sc[15];

    /*reaction 74: CH2O + O <=> OH + HCO */
    qf[73] = sc[3]*sc[16];
    qr[73] = sc[5]*sc[15];

    /*reaction 75: CH2O + OH <=> HCO + H2O */
    qf[74] = sc[5]*sc[16];
    qr[74] = sc[6]*sc[15];

    /*reaction 76: CH2O + O2 <=> HO2 + HCO */
    qf[75] = sc[4]*sc[16];
    qr[75] = sc[7]*sc[15];

    /*reaction 77: CH2O + CH2 <=> CH3 + HCO */
    qf[76] = sc[11]*sc[16];
    qr[76] = sc[13]*sc[15];

    /*reaction 78: CH2O + CH2(S) <=> CH3 + HCO */
    qf[77] = sc[12]*sc[16];
    qr[77] = sc[13]*sc[15];

    /*reaction 79: CH3 + O <=> H + CH2O */
    qf[78] = sc[3]*sc[13];
    qr[78] = sc[2]*sc[16];

    /*reaction 80: CH3 + O => H + H2 + CO */
    qf[79] = sc[3]*sc[13];
    qr[79] = 0.0;

    /*reaction 81: CH3 + OH <=> CH2 + H2O */
    qf[80] = sc[5]*sc[13];
    qr[80] = sc[6]*sc[11];

    /*reaction 82: CH3 + OH <=> CH2(S) + H2O */
    qf[81] = sc[5]*sc[13];
    qr[81] = sc[6]*sc[12];

    /*reaction 83: CH3 + OH <=> H2 + CH2O */
    qf[82] = sc[5]*sc[13];
    qr[82] = sc[1]*sc[16];

    /*reaction 84: CH3 + HO2 <=> O2 + CH4 */
    qf[83] = sc[7]*sc[13];
    qr[83] = sc[4]*sc[14];

    /*reaction 85: CH3 + HO2 <=> OH + CH3O */
    qf[84] = sc[7]*sc[13];
    qr[84] = sc[5]*sc[18];

    /*reaction 86: CH3 + O2 <=> O + CH3O */
    qf[85] = sc[4]*sc[13];
    qr[85] = sc[3]*sc[18];

    /*reaction 87: CH3 + O2 <=> OH + CH2O */
    qf[86] = sc[4]*sc[13];
    qr[86] = sc[5]*sc[16];

    /*reaction 88: CH3 + HCO <=> CH4 + CO */
    qf[87] = sc[13]*sc[15];
    qr[87] = sc[8]*sc[14];

    /*reaction 89: CH3 + CH2O <=> HCO + CH4 */
    qf[88] = sc[13]*sc[16];
    qr[88] = sc[14]*sc[15];

    /*reaction 90: CH3O + H <=> H + CH2OH */
    qf[89] = sc[2]*sc[18];
    qr[89] = sc[2]*sc[17];

    /*reaction 91: CH3O + H <=> H2 + CH2O */
    qf[90] = sc[2]*sc[18];
    qr[90] = sc[1]*sc[16];

    /*reaction 92: CH3O + H <=> OH + CH3 */
    qf[91] = sc[2]*sc[18];
    qr[91] = sc[5]*sc[13];

    /*reaction 93: CH3O + H <=> CH2(S) + H2O */
    qf[92] = sc[2]*sc[18];
    qr[92] = sc[6]*sc[12];

    /*reaction 94: CH3O + O <=> OH + CH2O */
    qf[93] = sc[3]*sc[18];
    qr[93] = sc[5]*sc[16];

    /*reaction 95: CH3O + OH <=> H2O + CH2O */
    qf[94] = sc[5]*sc[18];
    qr[94] = sc[6]*sc[16];

    /*reaction 96: CH3O + O2 <=> HO2 + CH2O */
    qf[95] = sc[4]*sc[18];
    qr[95] = sc[7]*sc[16];

    /*reaction 97: CH3O + CH3 <=> CH4 + CH2O */
    qf[96] = sc[13]*sc[18];
    qr[96] = sc[14]*sc[16];

    /*reaction 98: CH3O + CO <=> CH3 + CO2 */
    qf[97] = sc[8]*sc[18];
    qr[97] = sc[9]*sc[13];

    /*reaction 99: CH2OH + H <=> H2 + CH2O */
    qf[98] = sc[2]*sc[17];
    qr[98] = sc[1]*sc[16];

    /*reaction 100: CH2OH + H <=> OH + CH3 */
    qf[99] = sc[2]*sc[17];
    qr[99] = sc[5]*sc[13];

    /*reaction 101: CH2OH + H <=> CH2(S) + H2O */
    qf[100] = sc[2]*sc[17];
    qr[100] = sc[6]*sc[12];

    /*reaction 102: CH2OH + O <=> OH + CH2O */
    qf[101] = sc[3]*sc[17];
    qr[101] = sc[5]*sc[16];

    /*reaction 103: CH2OH + OH <=> H2O + CH2O */
    qf[102] = sc[5]*sc[17];
    qr[102] = sc[6]*sc[16];

    /*reaction 104: CH2OH + O2 <=> HO2 + CH2O */
    qf[103] = sc[4]*sc[17];
    qr[103] = sc[7]*sc[16];

    /*reaction 105: CH2OH + CH3 <=> CH4 + CH2O */
    qf[104] = sc[13]*sc[17];
    qr[104] = sc[14]*sc[16];

    /*reaction 106: CH4 + H <=> CH3 + H2 */
    qf[105] = sc[2]*sc[14];
    qr[105] = sc[1]*sc[13];

    /*reaction 107: CH4 + O <=> OH + CH3 */
    qf[106] = sc[3]*sc[14];
    qr[106] = sc[5]*sc[13];

    /*reaction 108: CH4 + OH <=> CH3 + H2O */
    qf[107] = sc[5]*sc[14];
    qr[107] = sc[6]*sc[13];

    /*reaction 109: CH4 + CH2 <=> 2.000000 CH3 */
    qf[108] = sc[11]*sc[14];
    qr[108] = pow(sc[13], 2.000000);

    /*reaction 110: CH4 + CH2(S) <=> 2.000000 CH3 */
    qf[109] = sc[12]*sc[14];
    qr[109] = pow(sc[13], 2.000000);

    /*reaction 111: CH3OH + H <=> CH2OH + H2 */
    qf[110] = sc[2]*sc[19];
    qr[110] = sc[1]*sc[17];

    /*reaction 112: CH3OH + H <=> CH3O + H2 */
    qf[111] = sc[2]*sc[19];
    qr[111] = sc[1]*sc[18];

    /*reaction 113: CH3OH + O <=> OH + CH2OH */
    qf[112] = sc[3]*sc[19];
    qr[112] = sc[5]*sc[17];

    /*reaction 114: CH3OH + O <=> OH + CH3O */
    qf[113] = sc[3]*sc[19];
    qr[113] = sc[5]*sc[18];

    /*reaction 115: CH3OH + OH <=> CH2OH + H2O */
    qf[114] = sc[5]*sc[19];
    qr[114] = sc[6]*sc[17];

    /*reaction 116: CH3OH + OH <=> CH3O + H2O */
    qf[115] = sc[5]*sc[19];
    qr[115] = sc[6]*sc[18];

    /*reaction 117: CH3OH + O2 <=> CH2OH + HO2 */
    qf[116] = sc[4]*sc[19];
    qr[116] = sc[7]*sc[17];

    /*reaction 118: CH3OH + CH <=> CH3 + CH2O */
    qf[117] = sc[10]*sc[19];
    qr[117] = sc[13]*sc[16];

    /*reaction 119: CH3OH + CH2 <=> CH3 + CH2OH */
    qf[118] = sc[11]*sc[19];
    qr[118] = sc[13]*sc[17];

    /*reaction 120: CH3OH + CH2 <=> CH3 + CH3O */
    qf[119] = sc[11]*sc[19];
    qr[119] = sc[13]*sc[18];

    /*reaction 121: CH3OH + CH2(S) <=> CH3 + CH3O */
    qf[120] = sc[12]*sc[19];
    qr[120] = sc[13]*sc[18];

    /*reaction 122: CH3OH + CH2(S) <=> CH3 + CH2OH */
    qf[121] = sc[12]*sc[19];
    qr[121] = sc[13]*sc[17];

    /*reaction 123: CH3OH + CH3 <=> CH2OH + CH4 */
    qf[122] = sc[13]*sc[19];
    qr[122] = sc[14]*sc[17];

    /*reaction 124: CH3OH + CH3 <=> CH3O + CH4 */
    qf[123] = sc[13]*sc[19];
    qr[123] = sc[14]*sc[18];

    amrex::Real T = tc[1];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 21; ++i) {
        mixture += sc[i];
    }

    amrex::Real Corr[124];
    for (int i = 0; i < 124; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        amrex::Real alpha[12];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[1] + (TB[0][2] - 1)*sc[4] + (TB[0][3] - 1)*sc[6] + (TB[0][4] - 1)*sc[8] + (TB[0][5] - 1)*sc[9] + (TB[0][6] - 1)*sc[14] + (TB[0][7] - 1)*sc[16] + (TB[0][8] - 1)*sc[19] + (TB[0][9] - 1)*sc[20];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[6] + (TB[1][1] - 1)*sc[8] + (TB[1][2] - 1)*sc[9] + (TB[1][3] - 1)*sc[14] + (TB[1][4] - 1)*sc[16] + (TB[1][5] - 1)*sc[19] + (TB[1][6] - 1)*sc[20];
        alpha[2] = mixture + (TB[2][0] - 1)*sc[6] + (TB[2][1] - 1)*sc[8] + (TB[2][2] - 1)*sc[9] + (TB[2][3] - 1)*sc[14] + (TB[2][4] - 1)*sc[16] + (TB[2][5] - 1)*sc[19] + (TB[2][6] - 1)*sc[20];
        alpha[3] = mixture + (TB[3][0] - 1)*sc[6] + (TB[3][1] - 1)*sc[8] + (TB[3][2] - 1)*sc[9] + (TB[3][3] - 1)*sc[14] + (TB[3][4] - 1)*sc[16] + (TB[3][5] - 1)*sc[19];
        alpha[4] = mixture + (TB[4][0] - 1)*sc[6] + (TB[4][1] - 1)*sc[8] + (TB[4][2] - 1)*sc[9] + (TB[4][3] - 1)*sc[14] + (TB[4][4] - 1)*sc[16] + (TB[4][5] - 1)*sc[19] + (TB[4][6] - 1)*sc[20];
        alpha[5] = mixture + (TB[5][0] - 1)*sc[6] + (TB[5][1] - 1)*sc[8] + (TB[5][2] - 1)*sc[9] + (TB[5][3] - 1)*sc[14] + (TB[5][4] - 1)*sc[16] + (TB[5][5] - 1)*sc[19] + (TB[5][6] - 1)*sc[20];
        alpha[6] = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[4] + (TB[6][2] - 1)*sc[6] + (TB[6][3] - 1)*sc[8] + (TB[6][4] - 1)*sc[9] + (TB[6][5] - 1)*sc[14] + (TB[6][6] - 1)*sc[16] + (TB[6][7] - 1)*sc[19] + (TB[6][8] - 1)*sc[20];
        alpha[7] = mixture + (TB[7][0] - 1)*sc[6] + (TB[7][1] - 1)*sc[8] + (TB[7][2] - 1)*sc[9] + (TB[7][3] - 1)*sc[14] + (TB[7][4] - 1)*sc[16] + (TB[7][5] - 1)*sc[19];
        alpha[8] = mixture + (TB[8][0] - 1)*sc[1] + (TB[8][1] - 1)*sc[6] + (TB[8][2] - 1)*sc[8] + (TB[8][3] - 1)*sc[9] + (TB[8][4] - 1)*sc[14] + (TB[8][5] - 1)*sc[16] + (TB[8][6] - 1)*sc[19] + (TB[8][7] - 1)*sc[20];
        alpha[9] = mixture + (TB[9][0] - 1)*sc[6] + (TB[9][1] - 1)*sc[8] + (TB[9][2] - 1)*sc[9] + (TB[9][3] - 1)*sc[14] + (TB[9][4] - 1)*sc[16] + (TB[9][5] - 1)*sc[19];
        alpha[10] = mixture + (TB[10][0] - 1)*sc[1] + (TB[10][1] - 1)*sc[6] + (TB[10][2] - 1)*sc[8] + (TB[10][3] - 1)*sc[9] + (TB[10][4] - 1)*sc[14] + (TB[10][5] - 1)*sc[16] + (TB[10][6] - 1)*sc[19] + (TB[10][7] - 1)*sc[20];
        alpha[11] = mixture + (TB[11][0] - 1)*sc[6] + (TB[11][1] - 1)*sc[8] + (TB[11][2] - 1)*sc[9] + (TB[11][3] - 1)*sc[14] + (TB[11][4] - 1)*sc[16] + (TB[11][5] - 1)*sc[19] + (TB[11][6] - 1)*sc[20];
        for (int i=0; i<12; i++)
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
        alpha = mixture + (TB[12][0] - 1)*sc[1] + (TB[12][1] - 1)*sc[6] + (TB[12][2] - 1)*sc[8] + (TB[12][3] - 1)*sc[9] + (TB[12][4] - 1)*sc[14] + (TB[12][5] - 1)*sc[16] + (TB[12][6] - 1)*sc[19];
        amrex::Real redP = alpha / k_f_save[12] * phase_units[12] * low_A[12] * exp(low_beta[12] * tc[0] - activation_units[12] * low_Ea[12] * invT);
        Corr[12] = redP / (1. + redP);
    }

    /* simple three-body correction */
    {
        amrex::Real alpha;
        alpha = mixture + (TB[13][0] - 1)*sc[0] + (TB[13][1] - 1)*sc[1] + (TB[13][2] - 1)*sc[6] + (TB[13][3] - 1)*sc[8] + (TB[13][4] - 1)*sc[9] + (TB[13][5] - 1)*sc[14] + (TB[13][6] - 1)*sc[16] + (TB[13][7] - 1)*sc[19] + (TB[13][8] - 1)*sc[20];
        Corr[13] = alpha;
        alpha = mixture + (TB[14][0] - 1)*sc[1] + (TB[14][1] - 1)*sc[6] + (TB[14][2] - 1)*sc[8] + (TB[14][3] - 1)*sc[9] + (TB[14][4] - 1)*sc[14] + (TB[14][5] - 1)*sc[16] + (TB[14][6] - 1)*sc[19] + (TB[14][7] - 1)*sc[20];
        Corr[14] = alpha;
        alpha = mixture + (TB[15][0] - 1)*sc[0] + (TB[15][1] - 1)*sc[1] + (TB[15][2] - 1)*sc[6] + (TB[15][3] - 1)*sc[8] + (TB[15][4] - 1)*sc[9] + (TB[15][5] - 1)*sc[14] + (TB[15][6] - 1)*sc[16] + (TB[15][7] - 1)*sc[19] + (TB[15][8] - 1)*sc[20];
        Corr[15] = alpha;
        alpha = mixture + (TB[16][0] - 1)*sc[0] + (TB[16][1] - 1)*sc[1] + (TB[16][2] - 1)*sc[4] + (TB[16][3] - 1)*sc[6] + (TB[16][4] - 1)*sc[8] + (TB[16][5] - 1)*sc[9] + (TB[16][6] - 1)*sc[14] + (TB[16][7] - 1)*sc[16] + (TB[16][8] - 1)*sc[19] + (TB[16][9] - 1)*sc[20];
        Corr[16] = alpha;
        alpha = mixture + (TB[17][0] - 1)*sc[0] + (TB[17][1] - 1)*sc[1] + (TB[17][2] - 1)*sc[4] + (TB[17][3] - 1)*sc[6] + (TB[17][4] - 1)*sc[8] + (TB[17][5] - 1)*sc[9] + (TB[17][6] - 1)*sc[14] + (TB[17][7] - 1)*sc[16] + (TB[17][8] - 1)*sc[19] + (TB[17][9] - 1)*sc[20];
        Corr[17] = alpha;
    }

    for (int i=0; i<124; i++)
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

    amrex::Real q_f[124], q_r[124];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 124; ++i) {
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
    amrex::Real c[21]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 124; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 21; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 124; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[21]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[21];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*N2 */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*H */
    YOW += y[3]*imw[3]; /*O */
    YOW += y[4]*imw[4]; /*O2 */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*HO2 */
    YOW += y[8]*imw[8]; /*CO */
    YOW += y[9]*imw[9]; /*CO2 */
    YOW += y[10]*imw[10]; /*CH */
    YOW += y[11]*imw[11]; /*CH2 */
    YOW += y[12]*imw[12]; /*CH2(S) */
    YOW += y[13]*imw[13]; /*CH3 */
    YOW += y[14]*imw[14]; /*CH4 */
    YOW += y[15]*imw[15]; /*HCO */
    YOW += y[16]*imw[16]; /*CH2O */
    YOW += y[17]*imw[17]; /*CH2OH */
    YOW += y[18]*imw[18]; /*CH3O */
    YOW += y[19]*imw[19]; /*CH3OH */
    YOW += y[20]*imw[20]; /*HE */
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
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 
    c[19] = PWORT * y[19]*imw[19]; 
    c[20] = PWORT * y[20]*imw[20]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 124; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[21]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 124; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[21]; /*temporary storage */
    amrex::Real imw[21];

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
    c[14] = 1e6 * (*rho) * y[14]*imw[14]; 
    c[15] = 1e6 * (*rho) * y[15]*imw[15]; 
    c[16] = 1e6 * (*rho) * y[16]*imw[16]; 
    c[17] = 1e6 * (*rho) * y[17]*imw[17]; 
    c[18] = 1e6 * (*rho) * y[18]*imw[18]; 
    c[19] = 1e6 * (*rho) * y[19]*imw[19]; 
    c[20] = 1e6 * (*rho) * y[20]*imw[20]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 124; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[21]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*15.999400; /*O */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*33.006770; /*HO2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*44.009950; /*CO2 */
    XW += x[10]*13.019120; /*CH */
    XW += x[11]*14.027090; /*CH2 */
    XW += x[12]*14.027090; /*CH2(S) */
    XW += x[13]*15.035060; /*CH3 */
    XW += x[14]*16.043030; /*CH4 */
    XW += x[15]*29.018520; /*HCO */
    XW += x[16]*30.026490; /*CH2O */
    XW += x[17]*31.034460; /*CH2OH */
    XW += x[18]*31.034460; /*CH3O */
    XW += x[19]*32.042430; /*CH3OH */
    XW += x[20]*4.002600; /*HE */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 124; ++id) {
        qdot[id] *= 1.0e-6;
    }
}

/*compute the reaction Jacobian on CPU */
void aJacobian(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP)
{
    for (int i=0; i<484; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[21];
    for (int k=0; k<21; k++) {
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
    for (int k = 0; k < 21; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[21];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[21];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[21];
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
    alpha = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[1] + (TB[0][2] - 1)*sc[4] + (TB[0][3] - 1)*sc[6] + (TB[0][4] - 1)*sc[8] + (TB[0][5] - 1)*sc[9] + (TB[0][6] - 1)*sc[14] + (TB[0][7] - 1)*sc[16] + (TB[0][8] - 1)*sc[19] + (TB[0][9] - 1)*sc[20];
    /* forward */
    phi_f = sc[2]*sc[4];
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
    Kc = refCinv * exp(g_RT[2] + g_RT[4] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (h_RT[7]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] -= q; /* O2 */
    wdot[7] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[N2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[2] -= dqdci;                /* dwdot[H]/d[N2] */
        J[4] -= dqdci;                /* dwdot[O2]/d[N2] */
        J[7] += dqdci;                /* dwdot[HO2]/d[N2] */
        /* d()/d[H2] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[24] -= dqdci;               /* dwdot[H]/d[H2] */
        J[26] -= dqdci;               /* dwdot[O2]/d[H2] */
        J[29] += dqdci;               /* dwdot[HO2]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[46] -= dqdci;               /* dwdot[H]/d[H] */
        J[48] -= dqdci;               /* dwdot[O2]/d[H] */
        J[51] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[O2] */
        dqdci = (TB[0][2] - 1)*dcdc_fac + k_f*sc[2];
        J[90] -= dqdci;               /* dwdot[H]/d[O2] */
        J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[95] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[0][3] - 1)*dcdc_fac;
        J[134] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[136] -= dqdci;              /* dwdot[O2]/d[H2O] */
        J[139] += dqdci;              /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[156] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[158] -= dqdci;              /* dwdot[O2]/d[HO2] */
        J[161] += dqdci;              /* dwdot[HO2]/d[HO2] */
        /* d()/d[CO] */
        dqdci = (TB[0][4] - 1)*dcdc_fac;
        J[178] -= dqdci;              /* dwdot[H]/d[CO] */
        J[180] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[183] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[0][5] - 1)*dcdc_fac;
        J[200] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[202] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[205] += dqdci;              /* dwdot[HO2]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[0][6] - 1)*dcdc_fac;
        J[310] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[312] -= dqdci;              /* dwdot[O2]/d[CH4] */
        J[315] += dqdci;              /* dwdot[HO2]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[0][7] - 1)*dcdc_fac;
        J[354] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[356] -= dqdci;              /* dwdot[O2]/d[CH2O] */
        J[359] += dqdci;              /* dwdot[HO2]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[0][8] - 1)*dcdc_fac;
        J[420] -= dqdci;              /* dwdot[H]/d[CH3OH] */
        J[422] -= dqdci;              /* dwdot[O2]/d[CH3OH] */
        J[425] += dqdci;              /* dwdot[HO2]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[0][9] - 1)*dcdc_fac;
        J[442] -= dqdci;              /* dwdot[H]/d[HE] */
        J[444] -= dqdci;              /* dwdot[O2]/d[HE] */
        J[447] += dqdci;              /* dwdot[HO2]/d[HE] */
    }
    else {
        dqdc[0] = TB[0][0]*dcdc_fac;
        dqdc[1] = TB[0][1]*dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[4];
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[0][2]*dcdc_fac + k_f*sc[2];
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[0][3]*dcdc_fac;
        dqdc[7] = dcdc_fac - k_r;
        dqdc[8] = TB[0][4]*dcdc_fac;
        dqdc[9] = TB[0][5]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[0][6]*dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[0][7]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = TB[0][8]*dcdc_fac;
        dqdc[20] = TB[0][9]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+2] -= dqdc[k];
            J[22*k+4] -= dqdc[k];
            J[22*k+7] += dqdc[k];
        }
    }
    J[464] -= dqdT; /* dwdot[H]/dT */
    J[466] -= dqdT; /* dwdot[O2]/dT */
    J[469] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: CH + H2 (+M) <=> CH3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[6] + (TB[1][1] - 1)*sc[8] + (TB[1][2] - 1)*sc[9] + (TB[1][3] - 1)*sc[14] + (TB[1][4] - 1)*sc[16] + (TB[1][5] - 1)*sc[19] + (TB[1][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[10];
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
    phi_r = sc[13];
    Kc = refCinv * exp(g_RT[1] + g_RT[10] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[13]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H2 */
    wdot[10] -= q; /* CH */
    wdot[13] += q; /* CH3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci =  + k_f*sc[10];
        J[23] -= dqdci;               /* dwdot[H2]/d[H2] */
        J[32] -= dqdci;               /* dwdot[CH]/d[H2] */
        J[35] += dqdci;               /* dwdot[CH3]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[133] -= dqdci;              /* dwdot[H2]/d[H2O] */
        J[142] -= dqdci;              /* dwdot[CH]/d[H2O] */
        J[145] += dqdci;              /* dwdot[CH3]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[1][1] - 1)*dcdc_fac;
        J[177] -= dqdci;              /* dwdot[H2]/d[CO] */
        J[186] -= dqdci;              /* dwdot[CH]/d[CO] */
        J[189] += dqdci;              /* dwdot[CH3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[1][2] - 1)*dcdc_fac;
        J[199] -= dqdci;              /* dwdot[H2]/d[CO2] */
        J[208] -= dqdci;              /* dwdot[CH]/d[CO2] */
        J[211] += dqdci;              /* dwdot[CH3]/d[CO2] */
        /* d()/d[CH] */
        dqdci =  + k_f*sc[1];
        J[221] -= dqdci;              /* dwdot[H2]/d[CH] */
        J[230] -= dqdci;              /* dwdot[CH]/d[CH] */
        J[233] += dqdci;              /* dwdot[CH3]/d[CH] */
        /* d()/d[CH3] */
        dqdci =  - k_r;
        J[287] -= dqdci;              /* dwdot[H2]/d[CH3] */
        J[296] -= dqdci;              /* dwdot[CH]/d[CH3] */
        J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[1][3] - 1)*dcdc_fac;
        J[309] -= dqdci;              /* dwdot[H2]/d[CH4] */
        J[318] -= dqdci;              /* dwdot[CH]/d[CH4] */
        J[321] += dqdci;              /* dwdot[CH3]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[1][4] - 1)*dcdc_fac;
        J[353] -= dqdci;              /* dwdot[H2]/d[CH2O] */
        J[362] -= dqdci;              /* dwdot[CH]/d[CH2O] */
        J[365] += dqdci;              /* dwdot[CH3]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[1][5] - 1)*dcdc_fac;
        J[419] -= dqdci;              /* dwdot[H2]/d[CH3OH] */
        J[428] -= dqdci;              /* dwdot[CH]/d[CH3OH] */
        J[431] += dqdci;              /* dwdot[CH3]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[1][6] - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H2]/d[HE] */
        J[450] -= dqdci;              /* dwdot[CH]/d[HE] */
        J[453] += dqdci;              /* dwdot[CH3]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[10];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[1][0]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[1][1]*dcdc_fac;
        dqdc[9] = TB[1][2]*dcdc_fac;
        dqdc[10] = dcdc_fac + k_f*sc[1];
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac - k_r;
        dqdc[14] = TB[1][3]*dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[1][4]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = TB[1][5]*dcdc_fac;
        dqdc[20] = TB[1][6]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+10] -= dqdc[k];
            J[22*k+13] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H2]/dT */
    J[472] -= dqdT; /* dwdot[CH]/dT */
    J[475] += dqdT; /* dwdot[CH3]/dT */

    /*reaction 3: CH2 + H (+M) <=> CH3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[6] + (TB[2][1] - 1)*sc[8] + (TB[2][2] - 1)*sc[9] + (TB[2][3] - 1)*sc[14] + (TB[2][4] - 1)*sc[16] + (TB[2][5] - 1)*sc[19] + (TB[2][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[2]*sc[11];
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
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[2]) > 1.e-100 ? (1.-troe_a[2])*exp(-T/troe_Tsss[2]) : 0.);
    Fcent2 = (fabs(troe_Ts[2]) > 1.e-100 ? troe_a[2] * exp(-T/troe_Ts[2]) : 0.);
    Fcent3 = (troe_len[2] == 4 ? exp(-troe_Tss[2] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[2]) > 1.e-100 ? -Fcent1/troe_Tsss[2] : 0.)
      + (fabs(troe_Ts[2]) > 1.e-100 ? -Fcent2/troe_Ts[2] : 0.)
      + (troe_len[2] == 4 ? Fcent3*troe_Tss[2]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[13];
    Kc = refCinv * exp(g_RT[2] + g_RT[11] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[11]) + (h_RT[13]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[11] -= q; /* CH2 */
    wdot[13] += q; /* CH3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[11];
        J[46] -= dqdci;               /* dwdot[H]/d[H] */
        J[55] -= dqdci;               /* dwdot[CH2]/d[H] */
        J[57] += dqdci;               /* dwdot[CH3]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[2][0] - 1)*dcdc_fac;
        J[134] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[143] -= dqdci;              /* dwdot[CH2]/d[H2O] */
        J[145] += dqdci;              /* dwdot[CH3]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[2][1] - 1)*dcdc_fac;
        J[178] -= dqdci;              /* dwdot[H]/d[CO] */
        J[187] -= dqdci;              /* dwdot[CH2]/d[CO] */
        J[189] += dqdci;              /* dwdot[CH3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[2][2] - 1)*dcdc_fac;
        J[200] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[209] -= dqdci;              /* dwdot[CH2]/d[CO2] */
        J[211] += dqdci;              /* dwdot[CH3]/d[CO2] */
        /* d()/d[CH2] */
        dqdci =  + k_f*sc[2];
        J[244] -= dqdci;              /* dwdot[H]/d[CH2] */
        J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
        J[255] += dqdci;              /* dwdot[CH3]/d[CH2] */
        /* d()/d[CH3] */
        dqdci =  - k_r;
        J[288] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[297] -= dqdci;              /* dwdot[CH2]/d[CH3] */
        J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[2][3] - 1)*dcdc_fac;
        J[310] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[319] -= dqdci;              /* dwdot[CH2]/d[CH4] */
        J[321] += dqdci;              /* dwdot[CH3]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[2][4] - 1)*dcdc_fac;
        J[354] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[363] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
        J[365] += dqdci;              /* dwdot[CH3]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[2][5] - 1)*dcdc_fac;
        J[420] -= dqdci;              /* dwdot[H]/d[CH3OH] */
        J[429] -= dqdci;              /* dwdot[CH2]/d[CH3OH] */
        J[431] += dqdci;              /* dwdot[CH3]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[2][6] - 1)*dcdc_fac;
        J[442] -= dqdci;              /* dwdot[H]/d[HE] */
        J[451] -= dqdci;              /* dwdot[CH2]/d[HE] */
        J[453] += dqdci;              /* dwdot[CH3]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[11];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[2][0]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[2][1]*dcdc_fac;
        dqdc[9] = TB[2][2]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac + k_f*sc[2];
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac - k_r;
        dqdc[14] = TB[2][3]*dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[2][4]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = TB[2][5]*dcdc_fac;
        dqdc[20] = TB[2][6]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+2] -= dqdc[k];
            J[22*k+11] -= dqdc[k];
            J[22*k+13] += dqdc[k];
        }
    }
    J[464] -= dqdT; /* dwdot[H]/dT */
    J[473] -= dqdT; /* dwdot[CH2]/dT */
    J[475] += dqdT; /* dwdot[CH3]/dT */

    /*reaction 4: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[6] + (TB[3][1] - 1)*sc[8] + (TB[3][2] - 1)*sc[9] + (TB[3][3] - 1)*sc[14] + (TB[3][4] - 1)*sc[16] + (TB[3][5] - 1)*sc[19];
    /* forward */
    phi_f = sc[6]*sc[12];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[3] * exp(low_beta[3] * tc[0] - activation_units[3] * low_Ea[3] * invT);
    Pr = phase_units[3] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[3] * invT + activation_units[3] * low_Ea[3] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[3]) > 1.e-100 ? (1.-troe_a[3])*exp(-T/troe_Tsss[3]) : 0.);
    Fcent2 = (fabs(troe_Ts[3]) > 1.e-100 ? troe_a[3] * exp(-T/troe_Ts[3]) : 0.);
    Fcent3 = (troe_len[3] == 4 ? exp(-troe_Tss[3] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[3]) > 1.e-100 ? -Fcent1/troe_Tsss[3] : 0.)
      + (fabs(troe_Ts[3]) > 1.e-100 ? -Fcent2/troe_Ts[3] : 0.)
      + (troe_len[3] == 4 ? Fcent3*troe_Tss[3]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[19];
    Kc = refCinv * exp(g_RT[6] + g_RT[12] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[19]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[6] -= q; /* H2O */
    wdot[12] -= q; /* CH2(S) */
    wdot[19] += q; /* CH3OH */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2O] */
        dqdci = (TB[3][0] - 1)*dcdc_fac + k_f*sc[12];
        J[138] -= dqdci;              /* dwdot[H2O]/d[H2O] */
        J[144] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
        J[151] += dqdci;              /* dwdot[CH3OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[3][1] - 1)*dcdc_fac;
        J[182] -= dqdci;              /* dwdot[H2O]/d[CO] */
        J[188] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
        J[195] += dqdci;              /* dwdot[CH3OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[3][2] - 1)*dcdc_fac;
        J[204] -= dqdci;              /* dwdot[H2O]/d[CO2] */
        J[210] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
        J[217] += dqdci;              /* dwdot[CH3OH]/d[CO2] */
        /* d()/d[CH2(S)] */
        dqdci =  + k_f*sc[6];
        J[270] -= dqdci;              /* dwdot[H2O]/d[CH2(S)] */
        J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
        J[283] += dqdci;              /* dwdot[CH3OH]/d[CH2(S)] */
        /* d()/d[CH4] */
        dqdci = (TB[3][3] - 1)*dcdc_fac;
        J[314] -= dqdci;              /* dwdot[H2O]/d[CH4] */
        J[320] -= dqdci;              /* dwdot[CH2(S)]/d[CH4] */
        J[327] += dqdci;              /* dwdot[CH3OH]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[3][4] - 1)*dcdc_fac;
        J[358] -= dqdci;              /* dwdot[H2O]/d[CH2O] */
        J[364] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
        J[371] += dqdci;              /* dwdot[CH3OH]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[3][5] - 1)*dcdc_fac - k_r;
        J[424] -= dqdci;              /* dwdot[H2O]/d[CH3OH] */
        J[430] -= dqdci;              /* dwdot[CH2(S)]/d[CH3OH] */
        J[437] += dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[3][0]*dcdc_fac + k_f*sc[12];
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[3][1]*dcdc_fac;
        dqdc[9] = TB[3][2]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac + k_f*sc[6];
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[3][3]*dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[3][4]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = TB[3][5]*dcdc_fac - k_r;
        dqdc[20] = dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+6] -= dqdc[k];
            J[22*k+12] -= dqdc[k];
            J[22*k+19] += dqdc[k];
        }
    }
    J[468] -= dqdT; /* dwdot[H2O]/dT */
    J[474] -= dqdT; /* dwdot[CH2(S)]/dT */
    J[481] += dqdT; /* dwdot[CH3OH]/dT */

    /*reaction 5: HCO + H (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[6] + (TB[4][1] - 1)*sc[8] + (TB[4][2] - 1)*sc[9] + (TB[4][3] - 1)*sc[14] + (TB[4][4] - 1)*sc[16] + (TB[4][5] - 1)*sc[19] + (TB[4][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[2]*sc[15];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[4] * exp(low_beta[4] * tc[0] - activation_units[4] * low_Ea[4] * invT);
    Pr = phase_units[4] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[4] * invT + activation_units[4] * low_Ea[4] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[4]) > 1.e-100 ? (1.-troe_a[4])*exp(-T/troe_Tsss[4]) : 0.);
    Fcent2 = (fabs(troe_Ts[4]) > 1.e-100 ? troe_a[4] * exp(-T/troe_Ts[4]) : 0.);
    Fcent3 = (troe_len[4] == 4 ? exp(-troe_Tss[4] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[4]) > 1.e-100 ? -Fcent1/troe_Tsss[4] : 0.)
      + (fabs(troe_Ts[4]) > 1.e-100 ? -Fcent2/troe_Ts[4] : 0.)
      + (troe_len[4] == 4 ? Fcent3*troe_Tss[4]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[16];
    Kc = refCinv * exp(g_RT[2] + g_RT[15] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[15]) + (h_RT[16]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[15] -= q; /* HCO */
    wdot[16] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[15];
        J[46] -= dqdci;               /* dwdot[H]/d[H] */
        J[59] -= dqdci;               /* dwdot[HCO]/d[H] */
        J[60] += dqdci;               /* dwdot[CH2O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[4][0] - 1)*dcdc_fac;
        J[134] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[147] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        J[148] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[4][1] - 1)*dcdc_fac;
        J[178] -= dqdci;              /* dwdot[H]/d[CO] */
        J[191] -= dqdci;              /* dwdot[HCO]/d[CO] */
        J[192] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[4][2] - 1)*dcdc_fac;
        J[200] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[213] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        J[214] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[4][3] - 1)*dcdc_fac;
        J[310] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[323] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        J[324] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[HCO] */
        dqdci =  + k_f*sc[2];
        J[332] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[345] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        J[346] += dqdci;              /* dwdot[CH2O]/d[HCO] */
        /* d()/d[CH2O] */
        dqdci = (TB[4][4] - 1)*dcdc_fac - k_r;
        J[354] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[367] -= dqdci;              /* dwdot[HCO]/d[CH2O] */
        J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[4][5] - 1)*dcdc_fac;
        J[420] -= dqdci;              /* dwdot[H]/d[CH3OH] */
        J[433] -= dqdci;              /* dwdot[HCO]/d[CH3OH] */
        J[434] += dqdci;              /* dwdot[CH2O]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[4][6] - 1)*dcdc_fac;
        J[442] -= dqdci;              /* dwdot[H]/d[HE] */
        J[455] -= dqdci;              /* dwdot[HCO]/d[HE] */
        J[456] += dqdci;              /* dwdot[CH2O]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[15];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[4][0]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[4][1]*dcdc_fac;
        dqdc[9] = TB[4][2]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[4][3]*dcdc_fac;
        dqdc[15] = dcdc_fac + k_f*sc[2];
        dqdc[16] = TB[4][4]*dcdc_fac - k_r;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = TB[4][5]*dcdc_fac;
        dqdc[20] = TB[4][6]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+2] -= dqdc[k];
            J[22*k+15] -= dqdc[k];
            J[22*k+16] += dqdc[k];
        }
    }
    J[464] -= dqdT; /* dwdot[H]/dT */
    J[477] -= dqdT; /* dwdot[HCO]/dT */
    J[478] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 6: CH2O (+M) <=> H2 + CO (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[6] + (TB[5][1] - 1)*sc[8] + (TB[5][2] - 1)*sc[9] + (TB[5][3] - 1)*sc[14] + (TB[5][4] - 1)*sc[16] + (TB[5][5] - 1)*sc[19] + (TB[5][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[16];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[5] * exp(low_beta[5] * tc[0] - activation_units[5] * low_Ea[5] * invT);
    Pr = phase_units[5] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[5] * invT + activation_units[5] * low_Ea[5] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[5]) > 1.e-100 ? (1.-troe_a[5])*exp(-T/troe_Tsss[5]) : 0.);
    Fcent2 = (fabs(troe_Ts[5]) > 1.e-100 ? troe_a[5] * exp(-T/troe_Ts[5]) : 0.);
    Fcent3 = (troe_len[5] == 4 ? exp(-troe_Tss[5] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[5]) > 1.e-100 ? -Fcent1/troe_Tsss[5] : 0.)
      + (fabs(troe_Ts[5]) > 1.e-100 ? -Fcent2/troe_Ts[5] : 0.)
      + (troe_len[5] == 4 ? Fcent3*troe_Tss[5]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[1]*sc[8];
    Kc = refC * exp(-g_RT[1] - g_RT[8] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[16]) + (h_RT[1] + h_RT[8]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[8] += q; /* CO */
    wdot[16] -= q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci =  - k_r*sc[8];
        J[23] += dqdci;               /* dwdot[H2]/d[H2] */
        J[30] += dqdci;               /* dwdot[CO]/d[H2] */
        J[38] -= dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[5][0] - 1)*dcdc_fac;
        J[133] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[140] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[148] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[5][1] - 1)*dcdc_fac - k_r*sc[1];
        J[177] += dqdci;              /* dwdot[H2]/d[CO] */
        J[184] += dqdci;              /* dwdot[CO]/d[CO] */
        J[192] -= dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[5][2] - 1)*dcdc_fac;
        J[199] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[206] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[214] -= dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[5][3] - 1)*dcdc_fac;
        J[309] += dqdci;              /* dwdot[H2]/d[CH4] */
        J[316] += dqdci;              /* dwdot[CO]/d[CH4] */
        J[324] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[5][4] - 1)*dcdc_fac + k_f;
        J[353] += dqdci;              /* dwdot[H2]/d[CH2O] */
        J[360] += dqdci;              /* dwdot[CO]/d[CH2O] */
        J[368] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[5][5] - 1)*dcdc_fac;
        J[419] += dqdci;              /* dwdot[H2]/d[CH3OH] */
        J[426] += dqdci;              /* dwdot[CO]/d[CH3OH] */
        J[434] -= dqdci;              /* dwdot[CH2O]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[5][6] - 1)*dcdc_fac;
        J[441] += dqdci;              /* dwdot[H2]/d[HE] */
        J[448] += dqdci;              /* dwdot[CO]/d[HE] */
        J[456] -= dqdci;              /* dwdot[CH2O]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac - k_r*sc[8];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[5][0]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[5][1]*dcdc_fac - k_r*sc[1];
        dqdc[9] = TB[5][2]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[5][3]*dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[5][4]*dcdc_fac + k_f;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = TB[5][5]*dcdc_fac;
        dqdc[20] = TB[5][6]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] += dqdc[k];
            J[22*k+8] += dqdc[k];
            J[22*k+16] -= dqdc[k];
        }
    }
    J[463] += dqdT; /* dwdot[H2]/dT */
    J[470] += dqdT; /* dwdot[CO]/dT */
    J[478] -= dqdT; /* dwdot[CH2O]/dT */

    /*reaction 7: CH3 + H (+M) <=> CH4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[4] + (TB[6][2] - 1)*sc[6] + (TB[6][3] - 1)*sc[8] + (TB[6][4] - 1)*sc[9] + (TB[6][5] - 1)*sc[14] + (TB[6][6] - 1)*sc[16] + (TB[6][7] - 1)*sc[19] + (TB[6][8] - 1)*sc[20];
    /* forward */
    phi_f = sc[2]*sc[13];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[6] * exp(low_beta[6] * tc[0] - activation_units[6] * low_Ea[6] * invT);
    Pr = phase_units[6] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[6] * invT + activation_units[6] * low_Ea[6] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[6]) > 1.e-100 ? (1.-troe_a[6])*exp(-T/troe_Tsss[6]) : 0.);
    Fcent2 = (fabs(troe_Ts[6]) > 1.e-100 ? troe_a[6] * exp(-T/troe_Ts[6]) : 0.);
    Fcent3 = (troe_len[6] == 4 ? exp(-troe_Tss[6] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[6]) > 1.e-100 ? -Fcent1/troe_Tsss[6] : 0.)
      + (fabs(troe_Ts[6]) > 1.e-100 ? -Fcent2/troe_Ts[6] : 0.)
      + (troe_len[6] == 4 ? Fcent3*troe_Tss[6]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[14];
    Kc = refCinv * exp(g_RT[2] + g_RT[13] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[14]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[13] -= q; /* CH3 */
    wdot[14] += q; /* CH4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[N2] */
        dqdci = (TB[6][0] - 1)*dcdc_fac;
        J[2] -= dqdci;                /* dwdot[H]/d[N2] */
        J[13] -= dqdci;               /* dwdot[CH3]/d[N2] */
        J[14] += dqdci;               /* dwdot[CH4]/d[N2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[13];
        J[46] -= dqdci;               /* dwdot[H]/d[H] */
        J[57] -= dqdci;               /* dwdot[CH3]/d[H] */
        J[58] += dqdci;               /* dwdot[CH4]/d[H] */
        /* d()/d[O2] */
        dqdci = (TB[6][1] - 1)*dcdc_fac;
        J[90] -= dqdci;               /* dwdot[H]/d[O2] */
        J[101] -= dqdci;              /* dwdot[CH3]/d[O2] */
        J[102] += dqdci;              /* dwdot[CH4]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[6][2] - 1)*dcdc_fac;
        J[134] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[145] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[146] += dqdci;              /* dwdot[CH4]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[6][3] - 1)*dcdc_fac;
        J[178] -= dqdci;              /* dwdot[H]/d[CO] */
        J[189] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[190] += dqdci;              /* dwdot[CH4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[6][4] - 1)*dcdc_fac;
        J[200] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[211] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[212] += dqdci;              /* dwdot[CH4]/d[CO2] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[2];
        J[288] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[300] += dqdci;              /* dwdot[CH4]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[6][5] - 1)*dcdc_fac - k_r;
        J[310] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[321] -= dqdci;              /* dwdot[CH3]/d[CH4] */
        J[322] += dqdci;              /* dwdot[CH4]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[6][6] - 1)*dcdc_fac;
        J[354] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[365] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
        J[366] += dqdci;              /* dwdot[CH4]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[6][7] - 1)*dcdc_fac;
        J[420] -= dqdci;              /* dwdot[H]/d[CH3OH] */
        J[431] -= dqdci;              /* dwdot[CH3]/d[CH3OH] */
        J[432] += dqdci;              /* dwdot[CH4]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[6][8] - 1)*dcdc_fac;
        J[442] -= dqdci;              /* dwdot[H]/d[HE] */
        J[453] -= dqdci;              /* dwdot[CH3]/d[HE] */
        J[454] += dqdci;              /* dwdot[CH4]/d[HE] */
    }
    else {
        dqdc[0] = TB[6][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[13];
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[6][1]*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[6][2]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[6][3]*dcdc_fac;
        dqdc[9] = TB[6][4]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac + k_f*sc[2];
        dqdc[14] = TB[6][5]*dcdc_fac - k_r;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[6][6]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = TB[6][7]*dcdc_fac;
        dqdc[20] = TB[6][8]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+2] -= dqdc[k];
            J[22*k+13] -= dqdc[k];
            J[22*k+14] += dqdc[k];
        }
    }
    J[464] -= dqdT; /* dwdot[H]/dT */
    J[475] -= dqdT; /* dwdot[CH3]/dT */
    J[476] += dqdT; /* dwdot[CH4]/dT */

    /*reaction 8: CH3 + OH (+M) <=> CH3OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[7][0] - 1)*sc[6] + (TB[7][1] - 1)*sc[8] + (TB[7][2] - 1)*sc[9] + (TB[7][3] - 1)*sc[14] + (TB[7][4] - 1)*sc[16] + (TB[7][5] - 1)*sc[19];
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[7] * exp(low_beta[7] * tc[0] - activation_units[7] * low_Ea[7] * invT);
    Pr = phase_units[7] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[7] * invT + activation_units[7] * low_Ea[7] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[7]) > 1.e-100 ? (1.-troe_a[7])*exp(-T/troe_Tsss[7]) : 0.);
    Fcent2 = (fabs(troe_Ts[7]) > 1.e-100 ? troe_a[7] * exp(-T/troe_Ts[7]) : 0.);
    Fcent3 = (troe_len[7] == 4 ? exp(-troe_Tss[7] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[7]) > 1.e-100 ? -Fcent1/troe_Tsss[7] : 0.)
      + (fabs(troe_Ts[7]) > 1.e-100 ? -Fcent2/troe_Ts[7] : 0.)
      + (troe_len[7] == 4 ? Fcent3*troe_Tss[7]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[19];
    Kc = refCinv * exp(g_RT[5] + g_RT[13] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[13]) + (h_RT[19]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[5] -= q; /* OH */
    wdot[13] -= q; /* CH3 */
    wdot[19] += q; /* CH3OH */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[OH] */
        dqdci =  + k_f*sc[13];
        J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
        J[123] -= dqdci;              /* dwdot[CH3]/d[OH] */
        J[129] += dqdci;              /* dwdot[CH3OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[7][0] - 1)*dcdc_fac;
        J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[145] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[151] += dqdci;              /* dwdot[CH3OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[7][1] - 1)*dcdc_fac;
        J[181] -= dqdci;              /* dwdot[OH]/d[CO] */
        J[189] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[195] += dqdci;              /* dwdot[CH3OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[7][2] - 1)*dcdc_fac;
        J[203] -= dqdci;              /* dwdot[OH]/d[CO2] */
        J[211] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[217] += dqdci;              /* dwdot[CH3OH]/d[CO2] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[5];
        J[291] -= dqdci;              /* dwdot[OH]/d[CH3] */
        J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[305] += dqdci;              /* dwdot[CH3OH]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[7][3] - 1)*dcdc_fac;
        J[313] -= dqdci;              /* dwdot[OH]/d[CH4] */
        J[321] -= dqdci;              /* dwdot[CH3]/d[CH4] */
        J[327] += dqdci;              /* dwdot[CH3OH]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[7][4] - 1)*dcdc_fac;
        J[357] -= dqdci;              /* dwdot[OH]/d[CH2O] */
        J[365] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
        J[371] += dqdci;              /* dwdot[CH3OH]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[7][5] - 1)*dcdc_fac - k_r;
        J[423] -= dqdci;              /* dwdot[OH]/d[CH3OH] */
        J[431] -= dqdci;              /* dwdot[CH3]/d[CH3OH] */
        J[437] += dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac + k_f*sc[13];
        dqdc[6] = TB[7][0]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[7][1]*dcdc_fac;
        dqdc[9] = TB[7][2]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac + k_f*sc[5];
        dqdc[14] = TB[7][3]*dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[7][4]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = TB[7][5]*dcdc_fac - k_r;
        dqdc[20] = dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+5] -= dqdc[k];
            J[22*k+13] -= dqdc[k];
            J[22*k+19] += dqdc[k];
        }
    }
    J[467] -= dqdT; /* dwdot[OH]/dT */
    J[475] -= dqdT; /* dwdot[CH3]/dT */
    J[481] += dqdT; /* dwdot[CH3OH]/dT */

    /*reaction 9: CH3O (+M) <=> H + CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[8][0] - 1)*sc[1] + (TB[8][1] - 1)*sc[6] + (TB[8][2] - 1)*sc[8] + (TB[8][3] - 1)*sc[9] + (TB[8][4] - 1)*sc[14] + (TB[8][5] - 1)*sc[16] + (TB[8][6] - 1)*sc[19] + (TB[8][7] - 1)*sc[20];
    /* forward */
    phi_f = sc[18];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[8] * exp(low_beta[8] * tc[0] - activation_units[8] * low_Ea[8] * invT);
    Pr = phase_units[8] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[8] * invT + activation_units[8] * low_Ea[8] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[8]) > 1.e-100 ? (1.-troe_a[8])*exp(-T/troe_Tsss[8]) : 0.);
    Fcent2 = (fabs(troe_Ts[8]) > 1.e-100 ? troe_a[8] * exp(-T/troe_Ts[8]) : 0.);
    Fcent3 = (troe_len[8] == 4 ? exp(-troe_Tss[8] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[8]) > 1.e-100 ? -Fcent1/troe_Tsss[8] : 0.)
      + (fabs(troe_Ts[8]) > 1.e-100 ? -Fcent2/troe_Ts[8] : 0.)
      + (troe_len[8] == 4 ? Fcent3*troe_Tss[8]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[2]*sc[16];
    Kc = refC * exp(-g_RT[2] - g_RT[16] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[18]) + (h_RT[2] + h_RT[16]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[16] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[8][0] - 1)*dcdc_fac;
        J[24] += dqdci;               /* dwdot[H]/d[H2] */
        J[38] += dqdci;               /* dwdot[CH2O]/d[H2] */
        J[40] -= dqdci;               /* dwdot[CH3O]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[16];
        J[46] += dqdci;               /* dwdot[H]/d[H] */
        J[60] += dqdci;               /* dwdot[CH2O]/d[H] */
        J[62] -= dqdci;               /* dwdot[CH3O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[8][1] - 1)*dcdc_fac;
        J[134] += dqdci;              /* dwdot[H]/d[H2O] */
        J[148] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        J[150] -= dqdci;              /* dwdot[CH3O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[8][2] - 1)*dcdc_fac;
        J[178] += dqdci;              /* dwdot[H]/d[CO] */
        J[192] += dqdci;              /* dwdot[CH2O]/d[CO] */
        J[194] -= dqdci;              /* dwdot[CH3O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[8][3] - 1)*dcdc_fac;
        J[200] += dqdci;              /* dwdot[H]/d[CO2] */
        J[214] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        J[216] -= dqdci;              /* dwdot[CH3O]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[8][4] - 1)*dcdc_fac;
        J[310] += dqdci;              /* dwdot[H]/d[CH4] */
        J[324] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        J[326] -= dqdci;              /* dwdot[CH3O]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[8][5] - 1)*dcdc_fac - k_r*sc[2];
        J[354] += dqdci;              /* dwdot[H]/d[CH2O] */
        J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        J[370] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
        /* d()/d[CH3O] */
        dqdci =  + k_f;
        J[398] += dqdci;              /* dwdot[H]/d[CH3O] */
        J[412] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
        J[414] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[8][6] - 1)*dcdc_fac;
        J[420] += dqdci;              /* dwdot[H]/d[CH3OH] */
        J[434] += dqdci;              /* dwdot[CH2O]/d[CH3OH] */
        J[436] -= dqdci;              /* dwdot[CH3O]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[8][7] - 1)*dcdc_fac;
        J[442] += dqdci;              /* dwdot[H]/d[HE] */
        J[456] += dqdci;              /* dwdot[CH2O]/d[HE] */
        J[458] -= dqdci;              /* dwdot[CH3O]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[8][0]*dcdc_fac;
        dqdc[2] = dcdc_fac - k_r*sc[16];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[8][1]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[8][2]*dcdc_fac;
        dqdc[9] = TB[8][3]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[8][4]*dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[8][5]*dcdc_fac - k_r*sc[2];
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac + k_f;
        dqdc[19] = TB[8][6]*dcdc_fac;
        dqdc[20] = TB[8][7]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+2] += dqdc[k];
            J[22*k+16] += dqdc[k];
            J[22*k+18] -= dqdc[k];
        }
    }
    J[464] += dqdT; /* dwdot[H]/dT */
    J[478] += dqdT; /* dwdot[CH2O]/dT */
    J[480] -= dqdT; /* dwdot[CH3O]/dT */

    /*reaction 10: CH3O + H (+M) <=> CH3OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[9][0] - 1)*sc[6] + (TB[9][1] - 1)*sc[8] + (TB[9][2] - 1)*sc[9] + (TB[9][3] - 1)*sc[14] + (TB[9][4] - 1)*sc[16] + (TB[9][5] - 1)*sc[19];
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[9] * exp(low_beta[9] * tc[0] - activation_units[9] * low_Ea[9] * invT);
    Pr = phase_units[9] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[9] * invT + activation_units[9] * low_Ea[9] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[9]) > 1.e-100 ? (1.-troe_a[9])*exp(-T/troe_Tsss[9]) : 0.);
    Fcent2 = (fabs(troe_Ts[9]) > 1.e-100 ? troe_a[9] * exp(-T/troe_Ts[9]) : 0.);
    Fcent3 = (troe_len[9] == 4 ? exp(-troe_Tss[9] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[9]) > 1.e-100 ? -Fcent1/troe_Tsss[9] : 0.)
      + (fabs(troe_Ts[9]) > 1.e-100 ? -Fcent2/troe_Ts[9] : 0.)
      + (troe_len[9] == 4 ? Fcent3*troe_Tss[9]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[19];
    Kc = refCinv * exp(g_RT[2] + g_RT[18] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[18]) + (h_RT[19]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[18] -= q; /* CH3O */
    wdot[19] += q; /* CH3OH */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[18];
        J[46] -= dqdci;               /* dwdot[H]/d[H] */
        J[62] -= dqdci;               /* dwdot[CH3O]/d[H] */
        J[63] += dqdci;               /* dwdot[CH3OH]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[9][0] - 1)*dcdc_fac;
        J[134] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[150] -= dqdci;              /* dwdot[CH3O]/d[H2O] */
        J[151] += dqdci;              /* dwdot[CH3OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[9][1] - 1)*dcdc_fac;
        J[178] -= dqdci;              /* dwdot[H]/d[CO] */
        J[194] -= dqdci;              /* dwdot[CH3O]/d[CO] */
        J[195] += dqdci;              /* dwdot[CH3OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[9][2] - 1)*dcdc_fac;
        J[200] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[216] -= dqdci;              /* dwdot[CH3O]/d[CO2] */
        J[217] += dqdci;              /* dwdot[CH3OH]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[9][3] - 1)*dcdc_fac;
        J[310] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[326] -= dqdci;              /* dwdot[CH3O]/d[CH4] */
        J[327] += dqdci;              /* dwdot[CH3OH]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[9][4] - 1)*dcdc_fac;
        J[354] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[370] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
        J[371] += dqdci;              /* dwdot[CH3OH]/d[CH2O] */
        /* d()/d[CH3O] */
        dqdci =  + k_f*sc[2];
        J[398] -= dqdci;              /* dwdot[H]/d[CH3O] */
        J[414] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
        J[415] += dqdci;              /* dwdot[CH3OH]/d[CH3O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[9][5] - 1)*dcdc_fac - k_r;
        J[420] -= dqdci;              /* dwdot[H]/d[CH3OH] */
        J[436] -= dqdci;              /* dwdot[CH3O]/d[CH3OH] */
        J[437] += dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[18];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[9][0]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[9][1]*dcdc_fac;
        dqdc[9] = TB[9][2]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[9][3]*dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[9][4]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac + k_f*sc[2];
        dqdc[19] = TB[9][5]*dcdc_fac - k_r;
        dqdc[20] = dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+2] -= dqdc[k];
            J[22*k+18] -= dqdc[k];
            J[22*k+19] += dqdc[k];
        }
    }
    J[464] -= dqdT; /* dwdot[H]/dT */
    J[480] -= dqdT; /* dwdot[CH3O]/dT */
    J[481] += dqdT; /* dwdot[CH3OH]/dT */

    /*reaction 11: CH2OH (+M) <=> H + CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[10][0] - 1)*sc[1] + (TB[10][1] - 1)*sc[6] + (TB[10][2] - 1)*sc[8] + (TB[10][3] - 1)*sc[9] + (TB[10][4] - 1)*sc[14] + (TB[10][5] - 1)*sc[16] + (TB[10][6] - 1)*sc[19] + (TB[10][7] - 1)*sc[20];
    /* forward */
    phi_f = sc[17];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[10] * exp(low_beta[10] * tc[0] - activation_units[10] * low_Ea[10] * invT);
    Pr = phase_units[10] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[10] * invT + activation_units[10] * low_Ea[10] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[10]) > 1.e-100 ? (1.-troe_a[10])*exp(-T/troe_Tsss[10]) : 0.);
    Fcent2 = (fabs(troe_Ts[10]) > 1.e-100 ? troe_a[10] * exp(-T/troe_Ts[10]) : 0.);
    Fcent3 = (troe_len[10] == 4 ? exp(-troe_Tss[10] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[10]) > 1.e-100 ? -Fcent1/troe_Tsss[10] : 0.)
      + (fabs(troe_Ts[10]) > 1.e-100 ? -Fcent2/troe_Ts[10] : 0.)
      + (troe_len[10] == 4 ? Fcent3*troe_Tss[10]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[2]*sc[16];
    Kc = refC * exp(-g_RT[2] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[17]) + (h_RT[2] + h_RT[16]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[16] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[10][0] - 1)*dcdc_fac;
        J[24] += dqdci;               /* dwdot[H]/d[H2] */
        J[38] += dqdci;               /* dwdot[CH2O]/d[H2] */
        J[39] -= dqdci;               /* dwdot[CH2OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[16];
        J[46] += dqdci;               /* dwdot[H]/d[H] */
        J[60] += dqdci;               /* dwdot[CH2O]/d[H] */
        J[61] -= dqdci;               /* dwdot[CH2OH]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[10][1] - 1)*dcdc_fac;
        J[134] += dqdci;              /* dwdot[H]/d[H2O] */
        J[148] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        J[149] -= dqdci;              /* dwdot[CH2OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[10][2] - 1)*dcdc_fac;
        J[178] += dqdci;              /* dwdot[H]/d[CO] */
        J[192] += dqdci;              /* dwdot[CH2O]/d[CO] */
        J[193] -= dqdci;              /* dwdot[CH2OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[10][3] - 1)*dcdc_fac;
        J[200] += dqdci;              /* dwdot[H]/d[CO2] */
        J[214] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        J[215] -= dqdci;              /* dwdot[CH2OH]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[10][4] - 1)*dcdc_fac;
        J[310] += dqdci;              /* dwdot[H]/d[CH4] */
        J[324] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        J[325] -= dqdci;              /* dwdot[CH2OH]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[10][5] - 1)*dcdc_fac - k_r*sc[2];
        J[354] += dqdci;              /* dwdot[H]/d[CH2O] */
        J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        J[369] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
        /* d()/d[CH2OH] */
        dqdci =  + k_f;
        J[376] += dqdci;              /* dwdot[H]/d[CH2OH] */
        J[390] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
        J[391] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
        /* d()/d[CH3OH] */
        dqdci = (TB[10][6] - 1)*dcdc_fac;
        J[420] += dqdci;              /* dwdot[H]/d[CH3OH] */
        J[434] += dqdci;              /* dwdot[CH2O]/d[CH3OH] */
        J[435] -= dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[10][7] - 1)*dcdc_fac;
        J[442] += dqdci;              /* dwdot[H]/d[HE] */
        J[456] += dqdci;              /* dwdot[CH2O]/d[HE] */
        J[457] -= dqdci;              /* dwdot[CH2OH]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[10][0]*dcdc_fac;
        dqdc[2] = dcdc_fac - k_r*sc[16];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[10][1]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[10][2]*dcdc_fac;
        dqdc[9] = TB[10][3]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[10][4]*dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[10][5]*dcdc_fac - k_r*sc[2];
        dqdc[17] = dcdc_fac + k_f;
        dqdc[18] = dcdc_fac;
        dqdc[19] = TB[10][6]*dcdc_fac;
        dqdc[20] = TB[10][7]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+2] += dqdc[k];
            J[22*k+16] += dqdc[k];
            J[22*k+17] -= dqdc[k];
        }
    }
    J[464] += dqdT; /* dwdot[H]/dT */
    J[478] += dqdT; /* dwdot[CH2O]/dT */
    J[479] -= dqdT; /* dwdot[CH2OH]/dT */

    /*reaction 12: CH2OH + H (+M) <=> CH3OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[11][0] - 1)*sc[6] + (TB[11][1] - 1)*sc[8] + (TB[11][2] - 1)*sc[9] + (TB[11][3] - 1)*sc[14] + (TB[11][4] - 1)*sc[16] + (TB[11][5] - 1)*sc[19] + (TB[11][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[11] * exp(low_beta[11] * tc[0] - activation_units[11] * low_Ea[11] * invT);
    Pr = phase_units[11] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[11] * invT + activation_units[11] * low_Ea[11] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[11]) > 1.e-100 ? (1.-troe_a[11])*exp(-T/troe_Tsss[11]) : 0.);
    Fcent2 = (fabs(troe_Ts[11]) > 1.e-100 ? troe_a[11] * exp(-T/troe_Ts[11]) : 0.);
    Fcent3 = (troe_len[11] == 4 ? exp(-troe_Tss[11] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[11]) > 1.e-100 ? -Fcent1/troe_Tsss[11] : 0.)
      + (fabs(troe_Ts[11]) > 1.e-100 ? -Fcent2/troe_Ts[11] : 0.)
      + (troe_len[11] == 4 ? Fcent3*troe_Tss[11]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[19];
    Kc = refCinv * exp(g_RT[2] + g_RT[17] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[19]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[17] -= q; /* CH2OH */
    wdot[19] += q; /* CH3OH */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[17];
        J[46] -= dqdci;               /* dwdot[H]/d[H] */
        J[61] -= dqdci;               /* dwdot[CH2OH]/d[H] */
        J[63] += dqdci;               /* dwdot[CH3OH]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[11][0] - 1)*dcdc_fac;
        J[134] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[149] -= dqdci;              /* dwdot[CH2OH]/d[H2O] */
        J[151] += dqdci;              /* dwdot[CH3OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[11][1] - 1)*dcdc_fac;
        J[178] -= dqdci;              /* dwdot[H]/d[CO] */
        J[193] -= dqdci;              /* dwdot[CH2OH]/d[CO] */
        J[195] += dqdci;              /* dwdot[CH3OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[11][2] - 1)*dcdc_fac;
        J[200] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[215] -= dqdci;              /* dwdot[CH2OH]/d[CO2] */
        J[217] += dqdci;              /* dwdot[CH3OH]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[11][3] - 1)*dcdc_fac;
        J[310] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[325] -= dqdci;              /* dwdot[CH2OH]/d[CH4] */
        J[327] += dqdci;              /* dwdot[CH3OH]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[11][4] - 1)*dcdc_fac;
        J[354] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[369] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
        J[371] += dqdci;              /* dwdot[CH3OH]/d[CH2O] */
        /* d()/d[CH2OH] */
        dqdci =  + k_f*sc[2];
        J[376] -= dqdci;              /* dwdot[H]/d[CH2OH] */
        J[391] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
        J[393] += dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
        /* d()/d[CH3OH] */
        dqdci = (TB[11][5] - 1)*dcdc_fac - k_r;
        J[420] -= dqdci;              /* dwdot[H]/d[CH3OH] */
        J[435] -= dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
        J[437] += dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[11][6] - 1)*dcdc_fac;
        J[442] -= dqdci;              /* dwdot[H]/d[HE] */
        J[457] -= dqdci;              /* dwdot[CH2OH]/d[HE] */
        J[459] += dqdci;              /* dwdot[CH3OH]/d[HE] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[17];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[11][0]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[11][1]*dcdc_fac;
        dqdc[9] = TB[11][2]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[11][3]*dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[11][4]*dcdc_fac;
        dqdc[17] = dcdc_fac + k_f*sc[2];
        dqdc[18] = dcdc_fac;
        dqdc[19] = TB[11][5]*dcdc_fac - k_r;
        dqdc[20] = TB[11][6]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+2] -= dqdc[k];
            J[22*k+17] -= dqdc[k];
            J[22*k+19] += dqdc[k];
        }
    }
    J[464] -= dqdT; /* dwdot[H]/dT */
    J[479] -= dqdT; /* dwdot[CH2OH]/dT */
    J[481] += dqdT; /* dwdot[CH3OH]/dT */

    /*reaction 13: CO + O (+M) <=> CO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[12][0] - 1)*sc[1] + (TB[12][1] - 1)*sc[6] + (TB[12][2] - 1)*sc[8] + (TB[12][3] - 1)*sc[9] + (TB[12][4] - 1)*sc[14] + (TB[12][5] - 1)*sc[16] + (TB[12][6] - 1)*sc[19];
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[12] * exp(low_beta[12] * tc[0] - activation_units[12] * low_Ea[12] * invT);
    Pr = phase_units[12] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[12] * invT + activation_units[12] * low_Ea[12] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Lindemann form */
    F = 1.0;
    dlogFdlogPr = 0.0;
    dlogFdT = 0.0;
    /* reverse */
    phi_r = sc[9];
    Kc = refCinv * exp(g_RT[3] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[9]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[12][0] - 1)*dcdc_fac;
        J[25] -= dqdci;               /* dwdot[O]/d[H2] */
        J[30] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[31] += dqdci;               /* dwdot[CO2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[8];
        J[69] -= dqdci;               /* dwdot[O]/d[O] */
        J[74] -= dqdci;               /* dwdot[CO]/d[O] */
        J[75] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (TB[12][1] - 1)*dcdc_fac;
        J[135] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[140] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[141] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[12][2] - 1)*dcdc_fac + k_f*sc[3];
        J[179] -= dqdci;              /* dwdot[O]/d[CO] */
        J[184] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[185] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[12][3] - 1)*dcdc_fac - k_r;
        J[201] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[206] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[207] += dqdci;              /* dwdot[CO2]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[12][4] - 1)*dcdc_fac;
        J[311] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[316] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[317] += dqdci;              /* dwdot[CO2]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[12][5] - 1)*dcdc_fac;
        J[355] -= dqdci;              /* dwdot[O]/d[CH2O] */
        J[360] -= dqdci;              /* dwdot[CO]/d[CH2O] */
        J[361] += dqdci;              /* dwdot[CO2]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[12][6] - 1)*dcdc_fac;
        J[421] -= dqdci;              /* dwdot[O]/d[CH3OH] */
        J[426] -= dqdci;              /* dwdot[CO]/d[CH3OH] */
        J[427] += dqdci;              /* dwdot[CO2]/d[CH3OH] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[12][0]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac + k_f*sc[8];
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[12][1]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[12][2]*dcdc_fac + k_f*sc[3];
        dqdc[9] = TB[12][3]*dcdc_fac - k_r;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = TB[12][4]*dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[12][5]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = TB[12][6]*dcdc_fac;
        dqdc[20] = dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+3] -= dqdc[k];
            J[22*k+8] -= dqdc[k];
            J[22*k+9] += dqdc[k];
        }
    }
    J[465] -= dqdT; /* dwdot[O]/dT */
    J[470] -= dqdT; /* dwdot[CO]/dT */
    J[471] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 14: H2 + M <=> 2.000000 H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[13][0] - 1)*sc[0] + (TB[13][1] - 1)*sc[1] + (TB[13][2] - 1)*sc[6] + (TB[13][3] - 1)*sc[8] + (TB[13][4] - 1)*sc[9] + (TB[13][5] - 1)*sc[14] + (TB[13][6] - 1)*sc[16] + (TB[13][7] - 1)*sc[19] + (TB[13][8] - 1)*sc[20];
    /* forward */
    phi_f = sc[1];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = pow(sc[2], 2.000000);
    Kc = refC * exp(g_RT[1] - 2.000000*g_RT[2]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1]) + (2.000000*h_RT[2]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H2 */
    wdot[2] += 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[N2] */
        dqdci = (TB[13][0] - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H2]/d[N2] */
        J[2] += 2 * dqdci;            /* dwdot[H]/d[N2] */
        /* d()/d[H2] */
        dqdci = (TB[13][1] - 1)*q_nocor + k_f;
        J[23] -= dqdci;               /* dwdot[H2]/d[H2] */
        J[24] += 2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*2.000000*sc[2];
        J[45] -= dqdci;               /* dwdot[H2]/d[H] */
        J[46] += 2 * dqdci;           /* dwdot[H]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[13][2] - 1)*q_nocor;
        J[133] -= dqdci;              /* dwdot[H2]/d[H2O] */
        J[134] += 2 * dqdci;          /* dwdot[H]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[13][3] - 1)*q_nocor;
        J[177] -= dqdci;              /* dwdot[H2]/d[CO] */
        J[178] += 2 * dqdci;          /* dwdot[H]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[13][4] - 1)*q_nocor;
        J[199] -= dqdci;              /* dwdot[H2]/d[CO2] */
        J[200] += 2 * dqdci;          /* dwdot[H]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[13][5] - 1)*q_nocor;
        J[309] -= dqdci;              /* dwdot[H2]/d[CH4] */
        J[310] += 2 * dqdci;          /* dwdot[H]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[13][6] - 1)*q_nocor;
        J[353] -= dqdci;              /* dwdot[H2]/d[CH2O] */
        J[354] += 2 * dqdci;          /* dwdot[H]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[13][7] - 1)*q_nocor;
        J[419] -= dqdci;              /* dwdot[H2]/d[CH3OH] */
        J[420] += 2 * dqdci;          /* dwdot[H]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[13][8] - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H2]/d[HE] */
        J[442] += 2 * dqdci;          /* dwdot[H]/d[HE] */
    }
    else {
        dqdc[0] = TB[13][0]*q_nocor;
        dqdc[1] = TB[13][1]*q_nocor + k_f;
        dqdc[2] = q_nocor - k_r*2.000000*sc[2];
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = TB[13][2]*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[13][3]*q_nocor;
        dqdc[9] = TB[13][4]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = TB[13][5]*q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = TB[13][6]*q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = TB[13][7]*q_nocor;
        dqdc[20] = TB[13][8]*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+2] += 2 * dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H2]/dT */
    J[464] += 2 * dqdT; /* dwdot[H]/dT */

    /*reaction 15: 2.000000 O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[14][0] - 1)*sc[1] + (TB[14][1] - 1)*sc[6] + (TB[14][2] - 1)*sc[8] + (TB[14][3] - 1)*sc[9] + (TB[14][4] - 1)*sc[14] + (TB[14][5] - 1)*sc[16] + (TB[14][6] - 1)*sc[19] + (TB[14][7] - 1)*sc[20];
    /* forward */
    phi_f = pow(sc[3], 2.000000);
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(2.000000*g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[3]) + (h_RT[4]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= 2 * q; /* O */
    wdot[4] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[14][0] - 1)*q_nocor;
        J[25] += -2 * dqdci;          /* dwdot[O]/d[H2] */
        J[26] += dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[3];
        J[69] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[70] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[91] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        J[92] += dqdci;               /* dwdot[O2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[14][1] - 1)*q_nocor;
        J[135] += -2 * dqdci;         /* dwdot[O]/d[H2O] */
        J[136] += dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[14][2] - 1)*q_nocor;
        J[179] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[180] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[14][3] - 1)*q_nocor;
        J[201] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[202] += dqdci;              /* dwdot[O2]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[14][4] - 1)*q_nocor;
        J[311] += -2 * dqdci;         /* dwdot[O]/d[CH4] */
        J[312] += dqdci;              /* dwdot[O2]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[14][5] - 1)*q_nocor;
        J[355] += -2 * dqdci;         /* dwdot[O]/d[CH2O] */
        J[356] += dqdci;              /* dwdot[O2]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[14][6] - 1)*q_nocor;
        J[421] += -2 * dqdci;         /* dwdot[O]/d[CH3OH] */
        J[422] += dqdci;              /* dwdot[O2]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[14][7] - 1)*q_nocor;
        J[443] += -2 * dqdci;         /* dwdot[O]/d[HE] */
        J[444] += dqdci;              /* dwdot[O2]/d[HE] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = TB[14][0]*q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor + k_f*2.000000*sc[3];
        dqdc[4] = q_nocor - k_r;
        dqdc[5] = q_nocor;
        dqdc[6] = TB[14][1]*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[14][2]*q_nocor;
        dqdc[9] = TB[14][3]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = TB[14][4]*q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = TB[14][5]*q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = TB[14][6]*q_nocor;
        dqdc[20] = TB[14][7]*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+3] += -2 * dqdc[k];
            J[22*k+4] += dqdc[k];
        }
    }
    J[465] += -2 * dqdT; /* dwdot[O]/dT */
    J[466] += dqdT; /* dwdot[O2]/dT */

    /*reaction 16: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[15][0] - 1)*sc[0] + (TB[15][1] - 1)*sc[1] + (TB[15][2] - 1)*sc[6] + (TB[15][3] - 1)*sc[8] + (TB[15][4] - 1)*sc[9] + (TB[15][5] - 1)*sc[14] + (TB[15][6] - 1)*sc[16] + (TB[15][7] - 1)*sc[19] + (TB[15][8] - 1)*sc[20];
    /* forward */
    phi_f = sc[2]*sc[3];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[2] + g_RT[3] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] -= q; /* O */
    wdot[5] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[N2] */
        dqdci = (TB[15][0] - 1)*q_nocor;
        J[2] -= dqdci;                /* dwdot[H]/d[N2] */
        J[3] -= dqdci;                /* dwdot[O]/d[N2] */
        J[5] += dqdci;                /* dwdot[OH]/d[N2] */
        /* d()/d[H2] */
        dqdci = (TB[15][1] - 1)*q_nocor;
        J[24] -= dqdci;               /* dwdot[H]/d[H2] */
        J[25] -= dqdci;               /* dwdot[O]/d[H2] */
        J[27] += dqdci;               /* dwdot[OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[3];
        J[46] -= dqdci;               /* dwdot[H]/d[H] */
        J[47] -= dqdci;               /* dwdot[O]/d[H] */
        J[49] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[2];
        J[68] -= dqdci;               /* dwdot[H]/d[O] */
        J[69] -= dqdci;               /* dwdot[O]/d[O] */
        J[71] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[112] -= dqdci;              /* dwdot[H]/d[OH] */
        J[113] -= dqdci;              /* dwdot[O]/d[OH] */
        J[115] += dqdci;              /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[15][2] - 1)*q_nocor;
        J[134] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[135] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[137] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[15][3] - 1)*q_nocor;
        J[178] -= dqdci;              /* dwdot[H]/d[CO] */
        J[179] -= dqdci;              /* dwdot[O]/d[CO] */
        J[181] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[15][4] - 1)*q_nocor;
        J[200] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[201] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[203] += dqdci;              /* dwdot[OH]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[15][5] - 1)*q_nocor;
        J[310] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[311] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[313] += dqdci;              /* dwdot[OH]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[15][6] - 1)*q_nocor;
        J[354] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[355] -= dqdci;              /* dwdot[O]/d[CH2O] */
        J[357] += dqdci;              /* dwdot[OH]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[15][7] - 1)*q_nocor;
        J[420] -= dqdci;              /* dwdot[H]/d[CH3OH] */
        J[421] -= dqdci;              /* dwdot[O]/d[CH3OH] */
        J[423] += dqdci;              /* dwdot[OH]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[15][8] - 1)*q_nocor;
        J[442] -= dqdci;              /* dwdot[H]/d[HE] */
        J[443] -= dqdci;              /* dwdot[O]/d[HE] */
        J[445] += dqdci;              /* dwdot[OH]/d[HE] */
    }
    else {
        dqdc[0] = TB[15][0]*q_nocor;
        dqdc[1] = TB[15][1]*q_nocor;
        dqdc[2] = q_nocor + k_f*sc[3];
        dqdc[3] = q_nocor + k_f*sc[2];
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor - k_r;
        dqdc[6] = TB[15][2]*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[15][3]*q_nocor;
        dqdc[9] = TB[15][4]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = TB[15][5]*q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = TB[15][6]*q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = TB[15][7]*q_nocor;
        dqdc[20] = TB[15][8]*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+2] -= dqdc[k];
            J[22*k+3] -= dqdc[k];
            J[22*k+5] += dqdc[k];
        }
    }
    J[464] -= dqdT; /* dwdot[H]/dT */
    J[465] -= dqdT; /* dwdot[O]/dT */
    J[467] += dqdT; /* dwdot[OH]/dT */

    /*reaction 17: H2O + M <=> H + OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[16][0] - 1)*sc[0] + (TB[16][1] - 1)*sc[1] + (TB[16][2] - 1)*sc[4] + (TB[16][3] - 1)*sc[6] + (TB[16][4] - 1)*sc[8] + (TB[16][5] - 1)*sc[9] + (TB[16][6] - 1)*sc[14] + (TB[16][7] - 1)*sc[16] + (TB[16][8] - 1)*sc[19] + (TB[16][9] - 1)*sc[20];
    /* forward */
    phi_f = sc[6];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = refC * exp(-g_RT[2] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6]) + (h_RT[2] + h_RT[5]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[5] += q; /* OH */
    wdot[6] -= q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[N2] */
        dqdci = (TB[16][0] - 1)*q_nocor;
        J[2] += dqdci;                /* dwdot[H]/d[N2] */
        J[5] += dqdci;                /* dwdot[OH]/d[N2] */
        J[6] -= dqdci;                /* dwdot[H2O]/d[N2] */
        /* d()/d[H2] */
        dqdci = (TB[16][1] - 1)*q_nocor;
        J[24] += dqdci;               /* dwdot[H]/d[H2] */
        J[27] += dqdci;               /* dwdot[OH]/d[H2] */
        J[28] -= dqdci;               /* dwdot[H2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[5];
        J[46] += dqdci;               /* dwdot[H]/d[H] */
        J[49] += dqdci;               /* dwdot[OH]/d[H] */
        J[50] -= dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[O2] */
        dqdci = (TB[16][2] - 1)*q_nocor;
        J[90] += dqdci;               /* dwdot[H]/d[O2] */
        J[93] += dqdci;               /* dwdot[OH]/d[O2] */
        J[94] -= dqdci;               /* dwdot[H2O]/d[O2] */
        /* d()/d[OH] */
        dqdci =  - k_r*sc[2];
        J[112] += dqdci;              /* dwdot[H]/d[OH] */
        J[115] += dqdci;              /* dwdot[OH]/d[OH] */
        J[116] -= dqdci;              /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[16][3] - 1)*q_nocor + k_f;
        J[134] += dqdci;              /* dwdot[H]/d[H2O] */
        J[137] += dqdci;              /* dwdot[OH]/d[H2O] */
        J[138] -= dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[16][4] - 1)*q_nocor;
        J[178] += dqdci;              /* dwdot[H]/d[CO] */
        J[181] += dqdci;              /* dwdot[OH]/d[CO] */
        J[182] -= dqdci;              /* dwdot[H2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[16][5] - 1)*q_nocor;
        J[200] += dqdci;              /* dwdot[H]/d[CO2] */
        J[203] += dqdci;              /* dwdot[OH]/d[CO2] */
        J[204] -= dqdci;              /* dwdot[H2O]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[16][6] - 1)*q_nocor;
        J[310] += dqdci;              /* dwdot[H]/d[CH4] */
        J[313] += dqdci;              /* dwdot[OH]/d[CH4] */
        J[314] -= dqdci;              /* dwdot[H2O]/d[CH4] */
        /* d()/d[CH2O] */
        dqdci = (TB[16][7] - 1)*q_nocor;
        J[354] += dqdci;              /* dwdot[H]/d[CH2O] */
        J[357] += dqdci;              /* dwdot[OH]/d[CH2O] */
        J[358] -= dqdci;              /* dwdot[H2O]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[16][8] - 1)*q_nocor;
        J[420] += dqdci;              /* dwdot[H]/d[CH3OH] */
        J[423] += dqdci;              /* dwdot[OH]/d[CH3OH] */
        J[424] -= dqdci;              /* dwdot[H2O]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[16][9] - 1)*q_nocor;
        J[442] += dqdci;              /* dwdot[H]/d[HE] */
        J[445] += dqdci;              /* dwdot[OH]/d[HE] */
        J[446] -= dqdci;              /* dwdot[H2O]/d[HE] */
    }
    else {
        dqdc[0] = TB[16][0]*q_nocor;
        dqdc[1] = TB[16][1]*q_nocor;
        dqdc[2] = q_nocor - k_r*sc[5];
        dqdc[3] = q_nocor;
        dqdc[4] = TB[16][2]*q_nocor;
        dqdc[5] = q_nocor - k_r*sc[2];
        dqdc[6] = TB[16][3]*q_nocor + k_f;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[16][4]*q_nocor;
        dqdc[9] = TB[16][5]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = TB[16][6]*q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = TB[16][7]*q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = TB[16][8]*q_nocor;
        dqdc[20] = TB[16][9]*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+2] += dqdc[k];
            J[22*k+5] += dqdc[k];
            J[22*k+6] -= dqdc[k];
        }
    }
    J[464] += dqdT; /* dwdot[H]/dT */
    J[467] += dqdT; /* dwdot[OH]/dT */
    J[468] -= dqdT; /* dwdot[H2O]/dT */

    /*reaction 18: HCO + M <=> H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[17][0] - 1)*sc[0] + (TB[17][1] - 1)*sc[1] + (TB[17][2] - 1)*sc[4] + (TB[17][3] - 1)*sc[6] + (TB[17][4] - 1)*sc[8] + (TB[17][5] - 1)*sc[9] + (TB[17][6] - 1)*sc[14] + (TB[17][7] - 1)*sc[16] + (TB[17][8] - 1)*sc[19] + (TB[17][9] - 1)*sc[20];
    /* forward */
    phi_f = sc[15];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[8];
    Kc = refC * exp(-g_RT[2] - g_RT[8] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[15]) + (h_RT[2] + h_RT[8]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[8] += q; /* CO */
    wdot[15] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[N2] */
        dqdci = (TB[17][0] - 1)*q_nocor;
        J[2] += dqdci;                /* dwdot[H]/d[N2] */
        J[8] += dqdci;                /* dwdot[CO]/d[N2] */
        J[15] -= dqdci;               /* dwdot[HCO]/d[N2] */
        /* d()/d[H2] */
        dqdci = (TB[17][1] - 1)*q_nocor;
        J[24] += dqdci;               /* dwdot[H]/d[H2] */
        J[30] += dqdci;               /* dwdot[CO]/d[H2] */
        J[37] -= dqdci;               /* dwdot[HCO]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[8];
        J[46] += dqdci;               /* dwdot[H]/d[H] */
        J[52] += dqdci;               /* dwdot[CO]/d[H] */
        J[59] -= dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[O2] */
        dqdci = (TB[17][2] - 1)*q_nocor;
        J[90] += dqdci;               /* dwdot[H]/d[O2] */
        J[96] += dqdci;               /* dwdot[CO]/d[O2] */
        J[103] -= dqdci;              /* dwdot[HCO]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[17][3] - 1)*q_nocor;
        J[134] += dqdci;              /* dwdot[H]/d[H2O] */
        J[140] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[147] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[17][4] - 1)*q_nocor - k_r*sc[2];
        J[178] += dqdci;              /* dwdot[H]/d[CO] */
        J[184] += dqdci;              /* dwdot[CO]/d[CO] */
        J[191] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[17][5] - 1)*q_nocor;
        J[200] += dqdci;              /* dwdot[H]/d[CO2] */
        J[206] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[213] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[17][6] - 1)*q_nocor;
        J[310] += dqdci;              /* dwdot[H]/d[CH4] */
        J[316] += dqdci;              /* dwdot[CO]/d[CH4] */
        J[323] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[332] += dqdci;              /* dwdot[H]/d[HCO] */
        J[338] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[345] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        /* d()/d[CH2O] */
        dqdci = (TB[17][7] - 1)*q_nocor;
        J[354] += dqdci;              /* dwdot[H]/d[CH2O] */
        J[360] += dqdci;              /* dwdot[CO]/d[CH2O] */
        J[367] -= dqdci;              /* dwdot[HCO]/d[CH2O] */
        /* d()/d[CH3OH] */
        dqdci = (TB[17][8] - 1)*q_nocor;
        J[420] += dqdci;              /* dwdot[H]/d[CH3OH] */
        J[426] += dqdci;              /* dwdot[CO]/d[CH3OH] */
        J[433] -= dqdci;              /* dwdot[HCO]/d[CH3OH] */
        /* d()/d[HE] */
        dqdci = (TB[17][9] - 1)*q_nocor;
        J[442] += dqdci;              /* dwdot[H]/d[HE] */
        J[448] += dqdci;              /* dwdot[CO]/d[HE] */
        J[455] -= dqdci;              /* dwdot[HCO]/d[HE] */
    }
    else {
        dqdc[0] = TB[17][0]*q_nocor;
        dqdc[1] = TB[17][1]*q_nocor;
        dqdc[2] = q_nocor - k_r*sc[8];
        dqdc[3] = q_nocor;
        dqdc[4] = TB[17][2]*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = TB[17][3]*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[17][4]*q_nocor - k_r*sc[2];
        dqdc[9] = TB[17][5]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = TB[17][6]*q_nocor;
        dqdc[15] = q_nocor + k_f;
        dqdc[16] = TB[17][7]*q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = TB[17][8]*q_nocor;
        dqdc[20] = TB[17][9]*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+2] += dqdc[k];
            J[22*k+8] += dqdc[k];
            J[22*k+15] -= dqdc[k];
        }
    }
    J[464] += dqdT; /* dwdot[H]/dT */
    J[470] += dqdT; /* dwdot[CO]/dT */
    J[477] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 19: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[2] - g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* O */
    wdot[4] -= q; /* O2 */
    wdot[5] += q; /* OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[4];
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[47] += dqdci;               /* dwdot[O]/d[H] */
    J[48] -= dqdci;               /* dwdot[O2]/d[H] */
    J[49] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[68] -= dqdci;               /* dwdot[H]/d[O] */
    J[69] += dqdci;               /* dwdot[O]/d[O] */
    J[70] -= dqdci;               /* dwdot[O2]/d[O] */
    J[71] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[2];
    J[90] -= dqdci;               /* dwdot[H]/d[O2] */
    J[91] += dqdci;               /* dwdot[O]/d[O2] */
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[93] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[112] -= dqdci;              /* dwdot[H]/d[OH] */
    J[113] += dqdci;              /* dwdot[O]/d[OH] */
    J[114] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[465] += dqdT;               /* dwdot[O]/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 20: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H2 */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[23] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[24] += dqdci;               /* dwdot[H]/d[H2] */
    J[25] -= dqdci;               /* dwdot[O]/d[H2] */
    J[27] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[45] -= dqdci;               /* dwdot[H2]/d[H] */
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[47] -= dqdci;               /* dwdot[O]/d[H] */
    J[49] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[67] -= dqdci;               /* dwdot[H2]/d[O] */
    J[68] += dqdci;               /* dwdot[H]/d[O] */
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[71] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[111] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[112] += dqdci;              /* dwdot[H]/d[OH] */
    J[113] -= dqdci;              /* dwdot[O]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H2]/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 21: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H2 */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[23] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[24] += dqdci;               /* dwdot[H]/d[H2] */
    J[25] -= dqdci;               /* dwdot[O]/d[H2] */
    J[27] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[45] -= dqdci;               /* dwdot[H2]/d[H] */
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[47] -= dqdci;               /* dwdot[O]/d[H] */
    J[49] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[67] -= dqdci;               /* dwdot[H2]/d[O] */
    J[68] += dqdci;               /* dwdot[H]/d[O] */
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[71] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[111] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[112] += dqdci;              /* dwdot[H]/d[OH] */
    J[113] -= dqdci;              /* dwdot[O]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H2]/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 22: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[6];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H2 */
    wdot[2] += q; /* H */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[5];
    J[23] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[24] += dqdci;               /* dwdot[H]/d[H2] */
    J[27] -= dqdci;               /* dwdot[OH]/d[H2] */
    J[28] += dqdci;               /* dwdot[H2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[6];
    J[45] -= dqdci;               /* dwdot[H2]/d[H] */
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[49] -= dqdci;               /* dwdot[OH]/d[H] */
    J[50] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[111] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[112] += dqdci;              /* dwdot[H]/d[OH] */
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[133] -= dqdci;              /* dwdot[H2]/d[H2O] */
    J[134] += dqdci;              /* dwdot[H]/d[H2O] */
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H2]/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 23: 2.000000 OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[5], 2.000000);
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[6];
    Kc = exp(-g_RT[3] + 2.000000*g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[5]) + (h_RT[3] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O */
    wdot[5] -= 2 * q; /* OH */
    wdot[6] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[6];
    J[69] += dqdci;               /* dwdot[O]/d[O] */
    J[71] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[72] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[5];
    J[113] += dqdci;              /* dwdot[O]/d[OH] */
    J[115] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[135] += dqdci;              /* dwdot[O]/d[H2O] */
    J[137] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O]/dT */
    J[467] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 24: H2 + HE <=> 2.000000 H + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[20];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = pow(sc[2], 2.000000)*sc[20];
    Kc = refC * exp(g_RT[1] - 2.000000*g_RT[2] + g_RT[20] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[20]) + (2.000000*h_RT[2] + h_RT[20]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H2 */
    wdot[2] += 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[20];
    J[23] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[24] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*2.000000*sc[2]*sc[20];
    J[45] -= dqdci;               /* dwdot[H2]/d[H] */
    J[46] += 2 * dqdci;           /* dwdot[H]/d[H] */
    /* d()/d[HE] */
    dqdci =  + k_f*sc[1] - k_r*pow(sc[2], 2.000000);
    J[441] -= dqdci;              /* dwdot[H2]/d[HE] */
    J[442] += 2 * dqdci;          /* dwdot[H]/d[HE] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H2]/dT */
    J[464] += 2 * dqdT;           /* dwdot[H]/dT */

    /*reaction 25: 2.000000 O + HE <=> O2 + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[3], 2.000000)*sc[20];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[20];
    Kc = refCinv * exp(2.000000*g_RT[3] - g_RT[4] + g_RT[20] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[3] + h_RT[20]) + (h_RT[4] + h_RT[20]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= 2 * q; /* O */
    wdot[4] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2.000000*sc[3]*sc[20];
    J[69] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[20];
    J[91] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[92] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[HE] */
    dqdci =  + k_f*pow(sc[3], 2.000000) - k_r*sc[4];
    J[443] += -2 * dqdci;         /* dwdot[O]/d[HE] */
    J[444] += dqdci;              /* dwdot[O2]/d[HE] */
    /* d()/dT */
    J[465] += -2 * dqdT;          /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 26: 2.000000 H2O <=> H + OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5]*sc[6];
    Kc = refC * exp(-g_RT[2] - g_RT[5] + 2.000000*g_RT[6] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[2] + h_RT[5] + h_RT[6]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[5] += q; /* OH */
    wdot[6] -= q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5]*sc[6];
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[49] += dqdci;               /* dwdot[OH]/d[H] */
    J[50] -= dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2]*sc[6];
    J[112] += dqdci;              /* dwdot[H]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[116] -= dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  + k_f*2.000000*sc[6] - k_r*sc[2]*sc[5];
    J[134] += dqdci;              /* dwdot[H]/d[H2O] */
    J[137] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[H2O]/dT */

    /*reaction 27: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[4] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[1] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* O2 */
    wdot[7] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[4];
    J[23] += dqdci;               /* dwdot[H2]/d[H2] */
    J[24] -= dqdci;               /* dwdot[H]/d[H2] */
    J[26] += dqdci;               /* dwdot[O2]/d[H2] */
    J[29] -= dqdci;               /* dwdot[HO2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[45] += dqdci;               /* dwdot[H2]/d[H] */
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[48] += dqdci;               /* dwdot[O2]/d[H] */
    J[51] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[1];
    J[89] += dqdci;               /* dwdot[H2]/d[O2] */
    J[90] -= dqdci;               /* dwdot[H]/d[O2] */
    J[92] += dqdci;               /* dwdot[O2]/d[O2] */
    J[95] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[155] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[156] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[158] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[161] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[466] += dqdT;               /* dwdot[O2]/dT */
    J[469] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 28: HO2 + H <=> 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* reverse */
    phi_r = pow(sc[5], 2.000000);
    Kc = exp(g_RT[2] - 2.000000*g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (2.000000*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[5] += 2 * q; /* OH */
    wdot[7] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[49] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[51] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[5];
    J[112] -= dqdci;              /* dwdot[H]/d[OH] */
    J[115] += 2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[117] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[156] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[159] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[161] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[467] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 29: HO2 + H <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
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
    wdot[2] -= q; /* H */
    wdot[3] += q; /* O */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[47] += dqdci;               /* dwdot[O]/d[H] */
    J[50] += dqdci;               /* dwdot[H2O]/d[H] */
    J[51] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[6];
    J[68] -= dqdci;               /* dwdot[H]/d[O] */
    J[69] += dqdci;               /* dwdot[O]/d[O] */
    J[72] += dqdci;               /* dwdot[H2O]/d[O] */
    J[73] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[134] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[135] += dqdci;              /* dwdot[O]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[139] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[156] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[157] += dqdci;              /* dwdot[O]/d[HO2] */
    J[160] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[161] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[465] += dqdT;               /* dwdot[O]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[469] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 30: HO2 + O <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[29] * fwd_A[29]
                * exp(fwd_beta[29] * tc[0] - activation_units[29] * fwd_Ea[29] * invT);
    dlnkfdT = fwd_beta[29] * invT + activation_units[29] * fwd_Ea[29] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[4] += q; /* O2 */
    wdot[5] += q; /* OH */
    wdot[7] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[O2]/d[O] */
    J[71] += dqdci;               /* dwdot[OH]/d[O] */
    J[73] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[5];
    J[91] -= dqdci;               /* dwdot[O]/d[O2] */
    J[92] += dqdci;               /* dwdot[O2]/d[O2] */
    J[93] += dqdci;               /* dwdot[OH]/d[O2] */
    J[95] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[113] -= dqdci;              /* dwdot[O]/d[OH] */
    J[114] += dqdci;              /* dwdot[O2]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[117] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[157] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[158] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[159] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[161] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[O2]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 31: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[30] * fwd_A[30]
                * exp(fwd_beta[30] * tc[0] - activation_units[30] * fwd_Ea[30] * invT);
    dlnkfdT = fwd_beta[30] * invT + activation_units[30] * fwd_Ea[30] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(-g_RT[4] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* O2 */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[6];
    J[92] += dqdci;               /* dwdot[O2]/d[O2] */
    J[93] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[94] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[95] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[114] += dqdci;              /* dwdot[O2]/d[OH] */
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[117] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[136] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[139] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[158] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[159] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[160] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[161] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[O2]/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[469] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 32: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[31] * fwd_A[31]
                * exp(fwd_beta[31] * tc[0] - activation_units[31] * fwd_Ea[31] * invT);
    dlnkfdT = fwd_beta[31] * invT + activation_units[31] * fwd_Ea[31] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(-g_RT[4] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* O2 */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[6];
    J[92] += dqdci;               /* dwdot[O2]/d[O2] */
    J[93] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[94] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[95] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[114] += dqdci;              /* dwdot[O2]/d[OH] */
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[117] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[136] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[139] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[158] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[159] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[160] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[161] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[O2]/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[469] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 33: CO + O2 <=> O + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[8];
    k_f = prefactor_units[32] * fwd_A[32]
                * exp(fwd_beta[32] * tc[0] - activation_units[32] * fwd_Ea[32] * invT);
    dlnkfdT = fwd_beta[32] * invT + activation_units[32] * fwd_Ea[32] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[9];
    Kc = exp(-g_RT[3] + g_RT[4] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[8]) + (h_RT[3] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O */
    wdot[4] -= q; /* O2 */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[9];
    J[69] += dqdci;               /* dwdot[O]/d[O] */
    J[70] -= dqdci;               /* dwdot[O2]/d[O] */
    J[74] -= dqdci;               /* dwdot[CO]/d[O] */
    J[75] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[91] += dqdci;               /* dwdot[O]/d[O2] */
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[96] -= dqdci;               /* dwdot[CO]/d[O2] */
    J[97] += dqdci;               /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[4];
    J[179] += dqdci;              /* dwdot[O]/d[CO] */
    J[180] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[184] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[185] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[3];
    J[201] += dqdci;              /* dwdot[O]/d[CO2] */
    J[202] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[206] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[207] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O]/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[470] -= dqdT;               /* dwdot[CO]/dT */
    J[471] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 34: CO + OH <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[8];
    k_f = prefactor_units[33] * fwd_A[33]
                * exp(fwd_beta[33] * tc[0] - activation_units[33] * fwd_Ea[33] * invT);
    dlnkfdT = fwd_beta[33] * invT + activation_units[33] * fwd_Ea[33] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[9];
    Kc = exp(-g_RT[2] + g_RT[5] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[8]) + (h_RT[2] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[5] -= q; /* OH */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[49] -= dqdci;               /* dwdot[OH]/d[H] */
    J[52] -= dqdci;               /* dwdot[CO]/d[H] */
    J[53] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[112] += dqdci;              /* dwdot[H]/d[OH] */
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[118] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[119] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[5];
    J[178] += dqdci;              /* dwdot[H]/d[CO] */
    J[181] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[184] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[185] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[2];
    J[200] += dqdci;              /* dwdot[H]/d[CO2] */
    J[203] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[206] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[207] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[470] -= dqdT;               /* dwdot[CO]/dT */
    J[471] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 35: CO + OH <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[8];
    k_f = prefactor_units[34] * fwd_A[34]
                * exp(fwd_beta[34] * tc[0] - activation_units[34] * fwd_Ea[34] * invT);
    dlnkfdT = fwd_beta[34] * invT + activation_units[34] * fwd_Ea[34] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[9];
    Kc = exp(-g_RT[2] + g_RT[5] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[8]) + (h_RT[2] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[5] -= q; /* OH */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[49] -= dqdci;               /* dwdot[OH]/d[H] */
    J[52] -= dqdci;               /* dwdot[CO]/d[H] */
    J[53] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[112] += dqdci;              /* dwdot[H]/d[OH] */
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[118] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[119] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[5];
    J[178] += dqdci;              /* dwdot[H]/d[CO] */
    J[181] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[184] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[185] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[2];
    J[200] += dqdci;              /* dwdot[H]/d[CO2] */
    J[203] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[206] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[207] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[470] -= dqdT;               /* dwdot[CO]/dT */
    J[471] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 36: CO + HO2 <=> OH + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[8];
    k_f = prefactor_units[35] * fwd_A[35]
                * exp(fwd_beta[35] * tc[0] - activation_units[35] * fwd_Ea[35] * invT);
    dlnkfdT = fwd_beta[35] * invT + activation_units[35] * fwd_Ea[35] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(-g_RT[5] + g_RT[7] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[8]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* OH */
    wdot[7] -= q; /* HO2 */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[117] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[118] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[119] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[8];
    J[159] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[161] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[162] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[163] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[7];
    J[181] += dqdci;              /* dwdot[OH]/d[CO] */
    J[183] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[184] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[185] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[5];
    J[203] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[205] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[206] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[207] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[HO2]/dT */
    J[470] -= dqdT;               /* dwdot[CO]/dT */
    J[471] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 37: HCO + H <=> H2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[15];
    k_f = prefactor_units[36] * fwd_A[36]
                * exp(fwd_beta[36] * tc[0] - activation_units[36] * fwd_Ea[36] * invT);
    dlnkfdT = fwd_beta[36] * invT + activation_units[36] * fwd_Ea[36] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[8];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[8] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[15]) + (h_RT[1] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[2] -= q; /* H */
    wdot[8] += q; /* CO */
    wdot[15] -= q; /* HCO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[8];
    J[23] += dqdci;               /* dwdot[H2]/d[H2] */
    J[24] -= dqdci;               /* dwdot[H]/d[H2] */
    J[30] += dqdci;               /* dwdot[CO]/d[H2] */
    J[37] -= dqdci;               /* dwdot[HCO]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[15];
    J[45] += dqdci;               /* dwdot[H2]/d[H] */
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[52] += dqdci;               /* dwdot[CO]/d[H] */
    J[59] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1];
    J[177] += dqdci;              /* dwdot[H2]/d[CO] */
    J[178] -= dqdci;              /* dwdot[H]/d[CO] */
    J[184] += dqdci;              /* dwdot[CO]/d[CO] */
    J[191] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[331] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[332] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[338] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[345] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[477] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 38: HCO + O <=> OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = prefactor_units[37] * fwd_A[37]
                * exp(fwd_beta[37] * tc[0] - activation_units[37] * fwd_Ea[37] * invT);
    dlnkfdT = fwd_beta[37] * invT + activation_units[37] * fwd_Ea[37] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[8];
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[8] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[15]) + (h_RT[5] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[8] += q; /* CO */
    wdot[15] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[15];
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[71] += dqdci;               /* dwdot[OH]/d[O] */
    J[74] += dqdci;               /* dwdot[CO]/d[O] */
    J[81] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[8];
    J[113] -= dqdci;              /* dwdot[O]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[118] += dqdci;              /* dwdot[CO]/d[OH] */
    J[125] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[179] -= dqdci;              /* dwdot[O]/d[CO] */
    J[181] += dqdci;              /* dwdot[OH]/d[CO] */
    J[184] += dqdci;              /* dwdot[CO]/d[CO] */
    J[191] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[333] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[335] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[338] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[345] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[477] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 39: HCO + O <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = prefactor_units[38] * fwd_A[38]
                * exp(fwd_beta[38] * tc[0] - activation_units[38] * fwd_Ea[38] * invT);
    dlnkfdT = fwd_beta[38] * invT + activation_units[38] * fwd_Ea[38] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[9];
    Kc = exp(-g_RT[2] + g_RT[3] - g_RT[9] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[15]) + (h_RT[2] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[9] += q; /* CO2 */
    wdot[15] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[47] -= dqdci;               /* dwdot[O]/d[H] */
    J[53] += dqdci;               /* dwdot[CO2]/d[H] */
    J[59] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[15];
    J[68] += dqdci;               /* dwdot[H]/d[O] */
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[75] += dqdci;               /* dwdot[CO2]/d[O] */
    J[81] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[2];
    J[200] += dqdci;              /* dwdot[H]/d[CO2] */
    J[201] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[207] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[213] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[332] += dqdci;              /* dwdot[H]/d[HCO] */
    J[333] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[339] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[345] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[471] += dqdT;               /* dwdot[CO2]/dT */
    J[477] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 40: HCO + OH <=> H2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[15];
    k_f = prefactor_units[39] * fwd_A[39]
                * exp(fwd_beta[39] * tc[0] - activation_units[39] * fwd_Ea[39] * invT);
    dlnkfdT = fwd_beta[39] * invT + activation_units[39] * fwd_Ea[39] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[8];
    Kc = exp(g_RT[5] - g_RT[6] - g_RT[8] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[15]) + (h_RT[6] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[8] += q; /* CO */
    wdot[15] -= q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[15];
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[118] += dqdci;              /* dwdot[CO]/d[OH] */
    J[125] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[140] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[147] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[181] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[182] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[184] += dqdci;              /* dwdot[CO]/d[CO] */
    J[191] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[335] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[336] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[338] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[345] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[477] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 41: HCO + O2 <=> HO2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[15];
    k_f = prefactor_units[40] * fwd_A[40]
                * exp(fwd_beta[40] * tc[0] - activation_units[40] * fwd_Ea[40] * invT);
    dlnkfdT = fwd_beta[40] * invT + activation_units[40] * fwd_Ea[40] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[8];
    Kc = exp(g_RT[4] - g_RT[7] - g_RT[8] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[15]) + (h_RT[7] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O2 */
    wdot[7] += q; /* HO2 */
    wdot[8] += q; /* CO */
    wdot[15] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[15];
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[95] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[96] += dqdci;               /* dwdot[CO]/d[O2] */
    J[103] -= dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[8];
    J[158] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[161] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[162] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[169] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[7];
    J[180] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[183] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[184] += dqdci;              /* dwdot[CO]/d[CO] */
    J[191] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[4];
    J[334] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[337] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[338] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[345] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[469] += dqdT;               /* dwdot[HO2]/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[477] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 42: CH + O <=> H + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[10];
    k_f = prefactor_units[41] * fwd_A[41]
                * exp(fwd_beta[41] * tc[0] - activation_units[41] * fwd_Ea[41] * invT);
    dlnkfdT = fwd_beta[41] * invT + activation_units[41] * fwd_Ea[41] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[8];
    Kc = exp(-g_RT[2] + g_RT[3] - g_RT[8] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[10]) + (h_RT[2] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[8] += q; /* CO */
    wdot[10] -= q; /* CH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[8];
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[47] -= dqdci;               /* dwdot[O]/d[H] */
    J[52] += dqdci;               /* dwdot[CO]/d[H] */
    J[54] -= dqdci;               /* dwdot[CH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[68] += dqdci;               /* dwdot[H]/d[O] */
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[74] += dqdci;               /* dwdot[CO]/d[O] */
    J[76] -= dqdci;               /* dwdot[CH]/d[O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[2];
    J[178] += dqdci;              /* dwdot[H]/d[CO] */
    J[179] -= dqdci;              /* dwdot[O]/d[CO] */
    J[184] += dqdci;              /* dwdot[CO]/d[CO] */
    J[186] -= dqdci;              /* dwdot[CH]/d[CO] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[3];
    J[222] += dqdci;              /* dwdot[H]/d[CH] */
    J[223] -= dqdci;              /* dwdot[O]/d[CH] */
    J[228] += dqdci;              /* dwdot[CO]/d[CH] */
    J[230] -= dqdci;              /* dwdot[CH]/d[CH] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[472] -= dqdT;               /* dwdot[CH]/dT */

    /*reaction 43: CH + OH <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[10];
    k_f = prefactor_units[42] * fwd_A[42]
                * exp(fwd_beta[42] * tc[0] - activation_units[42] * fwd_Ea[42] * invT);
    dlnkfdT = fwd_beta[42] * invT + activation_units[42] * fwd_Ea[42] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[15];
    Kc = exp(-g_RT[2] + g_RT[5] + g_RT[10] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[10]) + (h_RT[2] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[5] -= q; /* OH */
    wdot[10] -= q; /* CH */
    wdot[15] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[15];
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[49] -= dqdci;               /* dwdot[OH]/d[H] */
    J[54] -= dqdci;               /* dwdot[CH]/d[H] */
    J[59] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[112] += dqdci;              /* dwdot[H]/d[OH] */
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[120] -= dqdci;              /* dwdot[CH]/d[OH] */
    J[125] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[5];
    J[222] += dqdci;              /* dwdot[H]/d[CH] */
    J[225] -= dqdci;              /* dwdot[OH]/d[CH] */
    J[230] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[235] += dqdci;              /* dwdot[HCO]/d[CH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[2];
    J[332] += dqdci;              /* dwdot[H]/d[HCO] */
    J[335] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[340] -= dqdci;              /* dwdot[CH]/d[HCO] */
    J[345] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[472] -= dqdT;               /* dwdot[CH]/dT */
    J[477] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 44: CH + H2 <=> H + CH2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = prefactor_units[43] * fwd_A[43]
                * exp(fwd_beta[43] * tc[0] - activation_units[43] * fwd_Ea[43] * invT);
    dlnkfdT = fwd_beta[43] * invT + activation_units[43] * fwd_Ea[43] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[11];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[10] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[2] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H2 */
    wdot[2] += q; /* H */
    wdot[10] -= q; /* CH */
    wdot[11] += q; /* CH2 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[10];
    J[23] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[24] += dqdci;               /* dwdot[H]/d[H2] */
    J[32] -= dqdci;               /* dwdot[CH]/d[H2] */
    J[33] += dqdci;               /* dwdot[CH2]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[11];
    J[45] -= dqdci;               /* dwdot[H2]/d[H] */
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[54] -= dqdci;               /* dwdot[CH]/d[H] */
    J[55] += dqdci;               /* dwdot[CH2]/d[H] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[1];
    J[221] -= dqdci;              /* dwdot[H2]/d[CH] */
    J[222] += dqdci;              /* dwdot[H]/d[CH] */
    J[230] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[231] += dqdci;              /* dwdot[CH2]/d[CH] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[2];
    J[243] -= dqdci;              /* dwdot[H2]/d[CH2] */
    J[244] += dqdci;              /* dwdot[H]/d[CH2] */
    J[252] -= dqdci;              /* dwdot[CH]/d[CH2] */
    J[253] += dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H2]/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[472] -= dqdT;               /* dwdot[CH]/dT */
    J[473] += dqdT;               /* dwdot[CH2]/dT */

    /*reaction 45: CH + H2O <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[10];
    k_f = prefactor_units[44] * fwd_A[44]
                * exp(fwd_beta[44] * tc[0] - activation_units[44] * fwd_Ea[44] * invT);
    dlnkfdT = fwd_beta[44] * invT + activation_units[44] * fwd_Ea[44] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[16];
    Kc = exp(-g_RT[2] + g_RT[6] + g_RT[10] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[10]) + (h_RT[2] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[6] -= q; /* H2O */
    wdot[10] -= q; /* CH */
    wdot[16] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[50] -= dqdci;               /* dwdot[H2O]/d[H] */
    J[54] -= dqdci;               /* dwdot[CH]/d[H] */
    J[60] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[10];
    J[134] += dqdci;              /* dwdot[H]/d[H2O] */
    J[138] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[142] -= dqdci;              /* dwdot[CH]/d[H2O] */
    J[148] += dqdci;              /* dwdot[CH2O]/d[H2O] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[6];
    J[222] += dqdci;              /* dwdot[H]/d[CH] */
    J[226] -= dqdci;              /* dwdot[H2O]/d[CH] */
    J[230] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[236] += dqdci;              /* dwdot[CH2O]/d[CH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[2];
    J[354] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[358] -= dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[362] -= dqdci;              /* dwdot[CH]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[468] -= dqdT;               /* dwdot[H2O]/dT */
    J[472] -= dqdT;               /* dwdot[CH]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 46: CH + O2 <=> O + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[45] * fwd_A[45]
                * exp(fwd_beta[45] * tc[0] - activation_units[45] * fwd_Ea[45] * invT);
    dlnkfdT = fwd_beta[45] * invT + activation_units[45] * fwd_Ea[45] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[15];
    Kc = exp(-g_RT[3] + g_RT[4] + g_RT[10] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[3] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O */
    wdot[4] -= q; /* O2 */
    wdot[10] -= q; /* CH */
    wdot[15] += q; /* HCO */
    /* d()/d[O] */
    dqdci =  - k_r*sc[15];
    J[69] += dqdci;               /* dwdot[O]/d[O] */
    J[70] -= dqdci;               /* dwdot[O2]/d[O] */
    J[76] -= dqdci;               /* dwdot[CH]/d[O] */
    J[81] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[91] += dqdci;               /* dwdot[O]/d[O2] */
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[98] -= dqdci;               /* dwdot[CH]/d[O2] */
    J[103] += dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[4];
    J[223] += dqdci;              /* dwdot[O]/d[CH] */
    J[224] -= dqdci;              /* dwdot[O2]/d[CH] */
    J[230] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[235] += dqdci;              /* dwdot[HCO]/d[CH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[3];
    J[333] += dqdci;              /* dwdot[O]/d[HCO] */
    J[334] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[340] -= dqdci;              /* dwdot[CH]/d[HCO] */
    J[345] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O]/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[472] -= dqdT;               /* dwdot[CH]/dT */
    J[477] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 47: CH + O2 <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[46] * fwd_A[46]
                * exp(fwd_beta[46] * tc[0] - activation_units[46] * fwd_Ea[46] * invT);
    dlnkfdT = fwd_beta[46] * invT + activation_units[46] * fwd_Ea[46] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[9];
    Kc = exp(-g_RT[2] + g_RT[4] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[2] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* O2 */
    wdot[9] += q; /* CO2 */
    wdot[10] -= q; /* CH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[48] -= dqdci;               /* dwdot[O2]/d[H] */
    J[53] += dqdci;               /* dwdot[CO2]/d[H] */
    J[54] -= dqdci;               /* dwdot[CH]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[90] += dqdci;               /* dwdot[H]/d[O2] */
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[97] += dqdci;               /* dwdot[CO2]/d[O2] */
    J[98] -= dqdci;               /* dwdot[CH]/d[O2] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[2];
    J[200] += dqdci;              /* dwdot[H]/d[CO2] */
    J[202] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[207] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[208] -= dqdci;              /* dwdot[CH]/d[CO2] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[4];
    J[222] += dqdci;              /* dwdot[H]/d[CH] */
    J[224] -= dqdci;              /* dwdot[O2]/d[CH] */
    J[229] += dqdci;              /* dwdot[CO2]/d[CH] */
    J[230] -= dqdci;              /* dwdot[CH]/d[CH] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[471] += dqdT;               /* dwdot[CO2]/dT */
    J[472] -= dqdT;               /* dwdot[CH]/dT */

    /*reaction 48: CH + O2 <=> CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[47] * fwd_A[47]
                * exp(fwd_beta[47] * tc[0] - activation_units[47] * fwd_Ea[47] * invT);
    dlnkfdT = fwd_beta[47] * invT + activation_units[47] * fwd_Ea[47] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[8];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[8] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[5] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O2 */
    wdot[5] += q; /* OH */
    wdot[8] += q; /* CO */
    wdot[10] -= q; /* CH */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[93] += dqdci;               /* dwdot[OH]/d[O2] */
    J[96] += dqdci;               /* dwdot[CO]/d[O2] */
    J[98] -= dqdci;               /* dwdot[CH]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[8];
    J[114] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[118] += dqdci;              /* dwdot[CO]/d[OH] */
    J[120] -= dqdci;              /* dwdot[CH]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[180] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[181] += dqdci;              /* dwdot[OH]/d[CO] */
    J[184] += dqdci;              /* dwdot[CO]/d[CO] */
    J[186] -= dqdci;              /* dwdot[CH]/d[CO] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[4];
    J[224] -= dqdci;              /* dwdot[O2]/d[CH] */
    J[225] += dqdci;              /* dwdot[OH]/d[CH] */
    J[228] += dqdci;              /* dwdot[CO]/d[CH] */
    J[230] -= dqdci;              /* dwdot[CH]/d[CH] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[472] -= dqdT;               /* dwdot[CH]/dT */

    /*reaction 49: CH + O2 => O + H + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[48] * fwd_A[48]
                * exp(fwd_beta[48] * tc[0] - activation_units[48] * fwd_Ea[48] * invT);
    dlnkfdT = fwd_beta[48] * invT + activation_units[48] * fwd_Ea[48] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] += q; /* O */
    wdot[4] -= q; /* O2 */
    wdot[8] += q; /* CO */
    wdot[10] -= q; /* CH */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[90] += dqdci;               /* dwdot[H]/d[O2] */
    J[91] += dqdci;               /* dwdot[O]/d[O2] */
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[96] += dqdci;               /* dwdot[CO]/d[O2] */
    J[98] -= dqdci;               /* dwdot[CH]/d[O2] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[4];
    J[222] += dqdci;              /* dwdot[H]/d[CH] */
    J[223] += dqdci;              /* dwdot[O]/d[CH] */
    J[224] -= dqdci;              /* dwdot[O2]/d[CH] */
    J[228] += dqdci;              /* dwdot[CO]/d[CH] */
    J[230] -= dqdci;              /* dwdot[CH]/d[CH] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[465] += dqdT;               /* dwdot[O]/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[472] -= dqdT;               /* dwdot[CH]/dT */

    /*reaction 50: CH + CO2 <=> HCO + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[10];
    k_f = prefactor_units[49] * fwd_A[49]
                * exp(fwd_beta[49] * tc[0] - activation_units[49] * fwd_Ea[49] * invT);
    dlnkfdT = fwd_beta[49] * invT + activation_units[49] * fwd_Ea[49] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[15];
    Kc = exp(-g_RT[8] + g_RT[9] + g_RT[10] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[10]) + (h_RT[8] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] += q; /* CO */
    wdot[9] -= q; /* CO2 */
    wdot[10] -= q; /* CH */
    wdot[15] += q; /* HCO */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[15];
    J[184] += dqdci;              /* dwdot[CO]/d[CO] */
    J[185] -= dqdci;              /* dwdot[CO2]/d[CO] */
    J[186] -= dqdci;              /* dwdot[CH]/d[CO] */
    J[191] += dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[10];
    J[206] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[207] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    J[208] -= dqdci;              /* dwdot[CH]/d[CO2] */
    J[213] += dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[9];
    J[228] += dqdci;              /* dwdot[CO]/d[CH] */
    J[229] -= dqdci;              /* dwdot[CO2]/d[CH] */
    J[230] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[235] += dqdci;              /* dwdot[HCO]/d[CH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[8];
    J[338] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[339] -= dqdci;              /* dwdot[CO2]/d[HCO] */
    J[340] -= dqdci;              /* dwdot[CH]/d[HCO] */
    J[345] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[471] -= dqdT;               /* dwdot[CO2]/dT */
    J[472] -= dqdT;               /* dwdot[CH]/dT */
    J[477] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 51: CH2 + O => 2.000000 H + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = prefactor_units[50] * fwd_A[50]
                * exp(fwd_beta[50] * tc[0] - activation_units[50] * fwd_Ea[50] * invT);
    dlnkfdT = fwd_beta[50] * invT + activation_units[50] * fwd_Ea[50] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += 2 * q; /* H */
    wdot[3] -= q; /* O */
    wdot[8] += q; /* CO */
    wdot[11] -= q; /* CH2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[11];
    J[68] += 2 * dqdci;           /* dwdot[H]/d[O] */
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[74] += dqdci;               /* dwdot[CO]/d[O] */
    J[77] -= dqdci;               /* dwdot[CH2]/d[O] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[3];
    J[244] += 2 * dqdci;          /* dwdot[H]/d[CH2] */
    J[245] -= dqdci;              /* dwdot[O]/d[CH2] */
    J[250] += dqdci;              /* dwdot[CO]/d[CH2] */
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[464] += 2 * dqdT;           /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 52: CH2 + OH <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = prefactor_units[51] * fwd_A[51]
                * exp(fwd_beta[51] * tc[0] - activation_units[51] * fwd_Ea[51] * invT);
    dlnkfdT = fwd_beta[51] * invT + activation_units[51] * fwd_Ea[51] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[16];
    Kc = exp(-g_RT[2] + g_RT[5] + g_RT[11] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[11]) + (h_RT[2] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[5] -= q; /* OH */
    wdot[11] -= q; /* CH2 */
    wdot[16] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[49] -= dqdci;               /* dwdot[OH]/d[H] */
    J[55] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[60] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[112] += dqdci;              /* dwdot[H]/d[OH] */
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[121] -= dqdci;              /* dwdot[CH2]/d[OH] */
    J[126] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[5];
    J[244] += dqdci;              /* dwdot[H]/d[CH2] */
    J[247] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[258] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[2];
    J[354] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[357] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[363] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 53: CH2 + OH <=> CH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = prefactor_units[52] * fwd_A[52]
                * exp(fwd_beta[52] * tc[0] - activation_units[52] * fwd_Ea[52] * invT);
    dlnkfdT = fwd_beta[52] * invT + activation_units[52] * fwd_Ea[52] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[10];
    Kc = exp(g_RT[5] - g_RT[6] - g_RT[10] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[11]) + (h_RT[6] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[10] += q; /* CH */
    wdot[11] -= q; /* CH2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[120] += dqdci;              /* dwdot[CH]/d[OH] */
    J[121] -= dqdci;              /* dwdot[CH2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[10];
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[142] += dqdci;              /* dwdot[CH]/d[H2O] */
    J[143] -= dqdci;              /* dwdot[CH2]/d[H2O] */
    /* d()/d[CH] */
    dqdci =  - k_r*sc[6];
    J[225] -= dqdci;              /* dwdot[OH]/d[CH] */
    J[226] += dqdci;              /* dwdot[H2O]/d[CH] */
    J[230] += dqdci;              /* dwdot[CH]/d[CH] */
    J[231] -= dqdci;              /* dwdot[CH2]/d[CH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[5];
    J[247] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[248] += dqdci;              /* dwdot[H2O]/d[CH2] */
    J[252] += dqdci;              /* dwdot[CH]/d[CH2] */
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[472] += dqdT;               /* dwdot[CH]/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 54: CH2 + HO2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[11];
    k_f = prefactor_units[53] * fwd_A[53]
                * exp(fwd_beta[53] * tc[0] - activation_units[53] * fwd_Ea[53] * invT);
    dlnkfdT = fwd_beta[53] * invT + activation_units[53] * fwd_Ea[53] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[16];
    Kc = exp(-g_RT[5] + g_RT[7] + g_RT[11] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[11]) + (h_RT[5] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* OH */
    wdot[7] -= q; /* HO2 */
    wdot[11] -= q; /* CH2 */
    wdot[16] += q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[16];
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[117] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[121] -= dqdci;              /* dwdot[CH2]/d[OH] */
    J[126] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[159] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[161] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[165] -= dqdci;              /* dwdot[CH2]/d[HO2] */
    J[170] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[7];
    J[247] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[249] -= dqdci;              /* dwdot[HO2]/d[CH2] */
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[258] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[5];
    J[357] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[359] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[363] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[HO2]/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 55: CH2 + H2 <=> H + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = prefactor_units[54] * fwd_A[54]
                * exp(fwd_beta[54] * tc[0] - activation_units[54] * fwd_Ea[54] * invT);
    dlnkfdT = fwd_beta[54] * invT + activation_units[54] * fwd_Ea[54] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[13];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[11] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[11]) + (h_RT[2] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H2 */
    wdot[2] += q; /* H */
    wdot[11] -= q; /* CH2 */
    wdot[13] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[11];
    J[23] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[24] += dqdci;               /* dwdot[H]/d[H2] */
    J[33] -= dqdci;               /* dwdot[CH2]/d[H2] */
    J[35] += dqdci;               /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[45] -= dqdci;               /* dwdot[H2]/d[H] */
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[55] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[57] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[1];
    J[243] -= dqdci;              /* dwdot[H2]/d[CH2] */
    J[244] += dqdci;              /* dwdot[H]/d[CH2] */
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[255] += dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[2];
    J[287] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[288] += dqdci;              /* dwdot[H]/d[CH3] */
    J[297] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H2]/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 56: CH2 + O2 => OH + H + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = prefactor_units[55] * fwd_A[55]
                * exp(fwd_beta[55] * tc[0] - activation_units[55] * fwd_Ea[55] * invT);
    dlnkfdT = fwd_beta[55] * invT + activation_units[55] * fwd_Ea[55] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* O2 */
    wdot[5] += q; /* OH */
    wdot[8] += q; /* CO */
    wdot[11] -= q; /* CH2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[90] += dqdci;               /* dwdot[H]/d[O2] */
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[93] += dqdci;               /* dwdot[OH]/d[O2] */
    J[96] += dqdci;               /* dwdot[CO]/d[O2] */
    J[99] -= dqdci;               /* dwdot[CH2]/d[O2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[244] += dqdci;              /* dwdot[H]/d[CH2] */
    J[246] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[247] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[250] += dqdci;              /* dwdot[CO]/d[CH2] */
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 57: CH2 + O2 => 2.000000 H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = prefactor_units[56] * fwd_A[56]
                * exp(fwd_beta[56] * tc[0] - activation_units[56] * fwd_Ea[56] * invT);
    dlnkfdT = fwd_beta[56] * invT + activation_units[56] * fwd_Ea[56] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += 2 * q; /* H */
    wdot[4] -= q; /* O2 */
    wdot[9] += q; /* CO2 */
    wdot[11] -= q; /* CH2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[90] += 2 * dqdci;           /* dwdot[H]/d[O2] */
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[97] += dqdci;               /* dwdot[CO2]/d[O2] */
    J[99] -= dqdci;               /* dwdot[CH2]/d[O2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[244] += 2 * dqdci;          /* dwdot[H]/d[CH2] */
    J[246] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[251] += dqdci;              /* dwdot[CO2]/d[CH2] */
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[464] += 2 * dqdT;           /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[471] += dqdT;               /* dwdot[CO2]/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 58: CH2 + O2 <=> O + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = prefactor_units[57] * fwd_A[57]
                * exp(fwd_beta[57] * tc[0] - activation_units[57] * fwd_Ea[57] * invT);
    dlnkfdT = fwd_beta[57] * invT + activation_units[57] * fwd_Ea[57] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[16];
    Kc = exp(-g_RT[3] + g_RT[4] + g_RT[11] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[3] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O */
    wdot[4] -= q; /* O2 */
    wdot[11] -= q; /* CH2 */
    wdot[16] += q; /* CH2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[16];
    J[69] += dqdci;               /* dwdot[O]/d[O] */
    J[70] -= dqdci;               /* dwdot[O2]/d[O] */
    J[77] -= dqdci;               /* dwdot[CH2]/d[O] */
    J[82] += dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[91] += dqdci;               /* dwdot[O]/d[O2] */
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[99] -= dqdci;               /* dwdot[CH2]/d[O2] */
    J[104] += dqdci;              /* dwdot[CH2O]/d[O2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[245] += dqdci;              /* dwdot[O]/d[CH2] */
    J[246] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[258] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[3];
    J[355] += dqdci;              /* dwdot[O]/d[CH2O] */
    J[356] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[363] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O]/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 59: CH2 + O2 <=> H2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = prefactor_units[58] * fwd_A[58]
                * exp(fwd_beta[58] * tc[0] - activation_units[58] * fwd_Ea[58] * invT);
    dlnkfdT = fwd_beta[58] * invT + activation_units[58] * fwd_Ea[58] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(-g_RT[1] + g_RT[4] - g_RT[9] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[4] -= q; /* O2 */
    wdot[9] += q; /* CO2 */
    wdot[11] -= q; /* CH2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[9];
    J[23] += dqdci;               /* dwdot[H2]/d[H2] */
    J[26] -= dqdci;               /* dwdot[O2]/d[H2] */
    J[31] += dqdci;               /* dwdot[CO2]/d[H2] */
    J[33] -= dqdci;               /* dwdot[CH2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[89] += dqdci;               /* dwdot[H2]/d[O2] */
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[97] += dqdci;               /* dwdot[CO2]/d[O2] */
    J[99] -= dqdci;               /* dwdot[CH2]/d[O2] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[199] += dqdci;              /* dwdot[H2]/d[CO2] */
    J[202] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[207] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[209] -= dqdci;              /* dwdot[CH2]/d[CO2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[243] += dqdci;              /* dwdot[H2]/d[CH2] */
    J[246] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[251] += dqdci;              /* dwdot[CO2]/d[CH2] */
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[471] += dqdT;               /* dwdot[CO2]/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 60: CH2 + O2 <=> H2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = prefactor_units[59] * fwd_A[59]
                * exp(fwd_beta[59] * tc[0] - activation_units[59] * fwd_Ea[59] * invT);
    dlnkfdT = fwd_beta[59] * invT + activation_units[59] * fwd_Ea[59] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[8];
    Kc = exp(g_RT[4] - g_RT[6] - g_RT[8] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[6] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O2 */
    wdot[6] += q; /* H2O */
    wdot[8] += q; /* CO */
    wdot[11] -= q; /* CH2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[94] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[96] += dqdci;               /* dwdot[CO]/d[O2] */
    J[99] -= dqdci;               /* dwdot[CH2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[136] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[140] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[143] -= dqdci;              /* dwdot[CH2]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[180] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[182] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[184] += dqdci;              /* dwdot[CO]/d[CO] */
    J[187] -= dqdci;              /* dwdot[CH2]/d[CO] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[246] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[248] += dqdci;              /* dwdot[H2O]/d[CH2] */
    J[250] += dqdci;              /* dwdot[CO]/d[CH2] */
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 61: CH2(S) + N2 <=> CH2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[12];
    k_f = prefactor_units[60] * fwd_A[60]
                * exp(fwd_beta[60] * tc[0] - activation_units[60] * fwd_Ea[60] * invT);
    dlnkfdT = fwd_beta[60] * invT + activation_units[60] * fwd_Ea[60] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[11];
    Kc = exp(g_RT[0] - g_RT[0] - g_RT[11] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[12]) + (h_RT[0] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] += q; /* CH2 */
    wdot[12] -= q; /* CH2(S) */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[12] - k_r*sc[11];
    J[11] += dqdci;               /* dwdot[CH2]/d[N2] */
    J[12] -= dqdci;               /* dwdot[CH2(S)]/d[N2] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[0];
    J[253] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[254] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[0];
    J[275] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[473] += dqdT;               /* dwdot[CH2]/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 62: CH2(S) + HE <=> CH2 + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[20];
    k_f = prefactor_units[61] * fwd_A[61]
                * exp(fwd_beta[61] * tc[0] - activation_units[61] * fwd_Ea[61] * invT);
    dlnkfdT = fwd_beta[61] * invT + activation_units[61] * fwd_Ea[61] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[20];
    Kc = exp(-g_RT[11] + g_RT[12] + g_RT[20] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[20]) + (h_RT[11] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] += q; /* CH2 */
    wdot[12] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[20];
    J[253] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[254] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[20];
    J[275] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[HE] */
    dqdci =  + k_f*sc[12] - k_r*sc[11];
    J[451] += dqdci;              /* dwdot[CH2]/d[HE] */
    J[452] -= dqdci;              /* dwdot[CH2(S)]/d[HE] */
    /* d()/dT */
    J[473] += dqdT;               /* dwdot[CH2]/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 63: CH2(S) + H <=> CH + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[12];
    k_f = prefactor_units[62] * fwd_A[62]
                * exp(fwd_beta[62] * tc[0] - activation_units[62] * fwd_Ea[62] * invT);
    dlnkfdT = fwd_beta[62] * invT + activation_units[62] * fwd_Ea[62] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[10];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[10] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[12]) + (h_RT[1] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[2] -= q; /* H */
    wdot[10] += q; /* CH */
    wdot[12] -= q; /* CH2(S) */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[10];
    J[23] += dqdci;               /* dwdot[H2]/d[H2] */
    J[24] -= dqdci;               /* dwdot[H]/d[H2] */
    J[32] += dqdci;               /* dwdot[CH]/d[H2] */
    J[34] -= dqdci;               /* dwdot[CH2(S)]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[12];
    J[45] += dqdci;               /* dwdot[H2]/d[H] */
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[54] += dqdci;               /* dwdot[CH]/d[H] */
    J[56] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    /* d()/d[CH] */
    dqdci =  - k_r*sc[1];
    J[221] += dqdci;              /* dwdot[H2]/d[CH] */
    J[222] -= dqdci;              /* dwdot[H]/d[CH] */
    J[230] += dqdci;              /* dwdot[CH]/d[CH] */
    J[232] -= dqdci;              /* dwdot[CH2(S)]/d[CH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[2];
    J[265] += dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[266] -= dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[274] += dqdci;              /* dwdot[CH]/d[CH2(S)] */
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[472] += dqdT;               /* dwdot[CH]/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 64: CH2(S) + O => 2.000000 H + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[12];
    k_f = prefactor_units[63] * fwd_A[63]
                * exp(fwd_beta[63] * tc[0] - activation_units[63] * fwd_Ea[63] * invT);
    dlnkfdT = fwd_beta[63] * invT + activation_units[63] * fwd_Ea[63] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += 2 * q; /* H */
    wdot[3] -= q; /* O */
    wdot[8] += q; /* CO */
    wdot[12] -= q; /* CH2(S) */
    /* d()/d[O] */
    dqdci =  + k_f*sc[12];
    J[68] += 2 * dqdci;           /* dwdot[H]/d[O] */
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[74] += dqdci;               /* dwdot[CO]/d[O] */
    J[78] -= dqdci;               /* dwdot[CH2(S)]/d[O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[266] += 2 * dqdci;          /* dwdot[H]/d[CH2(S)] */
    J[267] -= dqdci;              /* dwdot[O]/d[CH2(S)] */
    J[272] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[464] += 2 * dqdT;           /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 65: CH2(S) + OH <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[12];
    k_f = prefactor_units[64] * fwd_A[64]
                * exp(fwd_beta[64] * tc[0] - activation_units[64] * fwd_Ea[64] * invT);
    dlnkfdT = fwd_beta[64] * invT + activation_units[64] * fwd_Ea[64] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[16];
    Kc = exp(-g_RT[2] + g_RT[5] + g_RT[12] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[12]) + (h_RT[2] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[5] -= q; /* OH */
    wdot[12] -= q; /* CH2(S) */
    wdot[16] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[49] -= dqdci;               /* dwdot[OH]/d[H] */
    J[56] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[60] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[112] += dqdci;              /* dwdot[H]/d[OH] */
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[122] -= dqdci;              /* dwdot[CH2(S)]/d[OH] */
    J[126] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[5];
    J[266] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[269] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[280] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[2];
    J[354] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[357] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[364] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[12];
    k_f = prefactor_units[65] * fwd_A[65]
                * exp(fwd_beta[65] * tc[0] - activation_units[65] * fwd_Ea[65] * invT);
    dlnkfdT = fwd_beta[65] * invT + activation_units[65] * fwd_Ea[65] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[13];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[12] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[12]) + (h_RT[2] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H2 */
    wdot[2] += q; /* H */
    wdot[12] -= q; /* CH2(S) */
    wdot[13] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[12];
    J[23] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[24] += dqdci;               /* dwdot[H]/d[H2] */
    J[34] -= dqdci;               /* dwdot[CH2(S)]/d[H2] */
    J[35] += dqdci;               /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[45] -= dqdci;               /* dwdot[H2]/d[H] */
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[56] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[57] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[1];
    J[265] -= dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[266] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[277] += dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[2];
    J[287] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[288] += dqdci;              /* dwdot[H]/d[CH3] */
    J[298] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H2]/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 67: CH2(S) + O2 <=> CH2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[12];
    k_f = prefactor_units[66] * fwd_A[66]
                * exp(fwd_beta[66] * tc[0] - activation_units[66] * fwd_Ea[66] * invT);
    dlnkfdT = fwd_beta[66] * invT + activation_units[66] * fwd_Ea[66] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[11];
    Kc = exp(g_RT[4] - g_RT[4] - g_RT[11] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[12]) + (h_RT[4] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] += q; /* CH2 */
    wdot[12] -= q; /* CH2(S) */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[12] - k_r*sc[11];
    J[99] += dqdci;               /* dwdot[CH2]/d[O2] */
    J[100] -= dqdci;              /* dwdot[CH2(S)]/d[O2] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[4];
    J[253] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[254] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[4];
    J[275] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[473] += dqdT;               /* dwdot[CH2]/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 68: CH2(S) + H2O <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[12];
    k_f = prefactor_units[67] * fwd_A[67]
                * exp(fwd_beta[67] * tc[0] - activation_units[67] * fwd_Ea[67] * invT);
    dlnkfdT = fwd_beta[67] * invT + activation_units[67] * fwd_Ea[67] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[11];
    Kc = exp(g_RT[6] - g_RT[6] - g_RT[11] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[6] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] += q; /* CH2 */
    wdot[12] -= q; /* CH2(S) */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[12] - k_r*sc[11];
    J[143] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[144] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[6];
    J[253] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[254] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[6];
    J[275] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[473] += dqdT;               /* dwdot[CH2]/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 69: CH2(S) + H2O <=> H2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[12];
    k_f = prefactor_units[68] * fwd_A[68]
                * exp(fwd_beta[68] * tc[0] - activation_units[68] * fwd_Ea[68] * invT);
    dlnkfdT = fwd_beta[68] * invT + activation_units[68] * fwd_Ea[68] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[12] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[6] -= q; /* H2O */
    wdot[12] -= q; /* CH2(S) */
    wdot[16] += q; /* CH2O */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[16];
    J[23] += dqdci;               /* dwdot[H2]/d[H2] */
    J[28] -= dqdci;               /* dwdot[H2O]/d[H2] */
    J[34] -= dqdci;               /* dwdot[CH2(S)]/d[H2] */
    J[38] += dqdci;               /* dwdot[CH2O]/d[H2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[12];
    J[133] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[138] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[144] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[148] += dqdci;              /* dwdot[CH2O]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[6];
    J[265] += dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[270] -= dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[280] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[353] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[358] -= dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[364] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[468] -= dqdT;               /* dwdot[H2O]/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = prefactor_units[69] * fwd_A[69]
                * exp(fwd_beta[69] * tc[0] - activation_units[69] * fwd_Ea[69] * invT);
    dlnkfdT = fwd_beta[69] * invT + activation_units[69] * fwd_Ea[69] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[11];
    Kc = exp(g_RT[8] - g_RT[8] - g_RT[11] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[8] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] += q; /* CH2 */
    wdot[12] -= q; /* CH2(S) */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[12] - k_r*sc[11];
    J[187] += dqdci;              /* dwdot[CH2]/d[CO] */
    J[188] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[8];
    J[253] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[254] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[8];
    J[275] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[473] += dqdT;               /* dwdot[CH2]/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[12];
    k_f = prefactor_units[70] * fwd_A[70]
                * exp(fwd_beta[70] * tc[0] - activation_units[70] * fwd_Ea[70] * invT);
    dlnkfdT = fwd_beta[70] * invT + activation_units[70] * fwd_Ea[70] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[11];
    Kc = exp(g_RT[9] - g_RT[9] - g_RT[11] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[12]) + (h_RT[9] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] += q; /* CH2 */
    wdot[12] -= q; /* CH2(S) */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[12] - k_r*sc[11];
    J[209] += dqdci;              /* dwdot[CH2]/d[CO2] */
    J[210] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[9];
    J[253] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[254] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[9];
    J[275] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[473] += dqdT;               /* dwdot[CH2]/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[12];
    k_f = prefactor_units[71] * fwd_A[71]
                * exp(fwd_beta[71] * tc[0] - activation_units[71] * fwd_Ea[71] * invT);
    dlnkfdT = fwd_beta[71] * invT + activation_units[71] * fwd_Ea[71] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[16];
    Kc = exp(-g_RT[8] + g_RT[9] + g_RT[12] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[12]) + (h_RT[8] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] += q; /* CO */
    wdot[9] -= q; /* CO2 */
    wdot[12] -= q; /* CH2(S) */
    wdot[16] += q; /* CH2O */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[16];
    J[184] += dqdci;              /* dwdot[CO]/d[CO] */
    J[185] -= dqdci;              /* dwdot[CO2]/d[CO] */
    J[188] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[192] += dqdci;              /* dwdot[CH2O]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[12];
    J[206] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[207] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    J[210] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    J[214] += dqdci;              /* dwdot[CH2O]/d[CO2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[9];
    J[272] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[273] -= dqdci;              /* dwdot[CO2]/d[CH2(S)] */
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[280] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[8];
    J[360] += dqdci;              /* dwdot[CO]/d[CH2O] */
    J[361] -= dqdci;              /* dwdot[CO2]/d[CH2O] */
    J[364] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[471] -= dqdT;               /* dwdot[CO2]/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 73: CH2O + H <=> HCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[16];
    k_f = prefactor_units[72] * fwd_A[72]
                * exp(fwd_beta[72] * tc[0] - activation_units[72] * fwd_Ea[72] * invT);
    dlnkfdT = fwd_beta[72] * invT + activation_units[72] * fwd_Ea[72] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[15];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[16]) + (h_RT[1] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[2] -= q; /* H */
    wdot[15] += q; /* HCO */
    wdot[16] -= q; /* CH2O */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[15];
    J[23] += dqdci;               /* dwdot[H2]/d[H2] */
    J[24] -= dqdci;               /* dwdot[H]/d[H2] */
    J[37] += dqdci;               /* dwdot[HCO]/d[H2] */
    J[38] -= dqdci;               /* dwdot[CH2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[16];
    J[45] += dqdci;               /* dwdot[H2]/d[H] */
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[59] += dqdci;               /* dwdot[HCO]/d[H] */
    J[60] -= dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[331] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[332] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[345] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[346] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[2];
    J[353] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[354] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[367] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[368] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[477] += dqdT;               /* dwdot[HCO]/dT */
    J[478] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 74: CH2O + O <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[16];
    k_f = prefactor_units[73] * fwd_A[73]
                * exp(fwd_beta[73] * tc[0] - activation_units[73] * fwd_Ea[73] * invT);
    dlnkfdT = fwd_beta[73] * invT + activation_units[73] * fwd_Ea[73] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[15];
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[16]) + (h_RT[5] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[15] += q; /* HCO */
    wdot[16] -= q; /* CH2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[71] += dqdci;               /* dwdot[OH]/d[O] */
    J[81] += dqdci;               /* dwdot[HCO]/d[O] */
    J[82] -= dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[113] -= dqdci;              /* dwdot[O]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[125] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[126] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[5];
    J[333] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[335] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[345] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[346] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[3];
    J[355] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[357] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[367] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[368] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[477] += dqdT;               /* dwdot[HCO]/dT */
    J[478] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 75: CH2O + OH <=> HCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[16];
    k_f = prefactor_units[74] * fwd_A[74]
                * exp(fwd_beta[74] * tc[0] - activation_units[74] * fwd_Ea[74] * invT);
    dlnkfdT = fwd_beta[74] * invT + activation_units[74] * fwd_Ea[74] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[15];
    Kc = exp(g_RT[5] - g_RT[6] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[16]) + (h_RT[6] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[15] += q; /* HCO */
    wdot[16] -= q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[16];
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[125] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[126] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[15];
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[147] += dqdci;              /* dwdot[HCO]/d[H2O] */
    J[148] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[6];
    J[335] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[336] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[345] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[346] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[5];
    J[357] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[358] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[367] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[368] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[477] += dqdT;               /* dwdot[HCO]/dT */
    J[478] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 76: CH2O + O2 <=> HO2 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[16];
    k_f = prefactor_units[75] * fwd_A[75]
                * exp(fwd_beta[75] * tc[0] - activation_units[75] * fwd_Ea[75] * invT);
    dlnkfdT = fwd_beta[75] * invT + activation_units[75] * fwd_Ea[75] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[15];
    Kc = exp(g_RT[4] - g_RT[7] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[16]) + (h_RT[7] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O2 */
    wdot[7] += q; /* HO2 */
    wdot[15] += q; /* HCO */
    wdot[16] -= q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[16];
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[95] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[103] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[104] -= dqdci;              /* dwdot[CH2O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[15];
    J[158] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[161] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[169] += dqdci;              /* dwdot[HCO]/d[HO2] */
    J[170] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[7];
    J[334] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[337] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[345] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[346] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[4];
    J[356] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[359] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[367] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[368] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[469] += dqdT;               /* dwdot[HO2]/dT */
    J[477] += dqdT;               /* dwdot[HCO]/dT */
    J[478] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 77: CH2O + CH2 <=> CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[16];
    k_f = prefactor_units[76] * fwd_A[76]
                * exp(fwd_beta[76] * tc[0] - activation_units[76] * fwd_Ea[76] * invT);
    dlnkfdT = fwd_beta[76] * invT + activation_units[76] * fwd_Ea[76] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[15];
    Kc = exp(g_RT[11] - g_RT[13] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[16]) + (h_RT[13] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] -= q; /* CH2 */
    wdot[13] += q; /* CH3 */
    wdot[15] += q; /* HCO */
    wdot[16] -= q; /* CH2O */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[16];
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[255] += dqdci;              /* dwdot[CH3]/d[CH2] */
    J[257] += dqdci;              /* dwdot[HCO]/d[CH2] */
    J[258] -= dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[15];
    J[297] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[301] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[302] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[13];
    J[341] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[343] += dqdci;              /* dwdot[CH3]/d[HCO] */
    J[345] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[346] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[11];
    J[363] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[365] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[367] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[368] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[477] += dqdT;               /* dwdot[HCO]/dT */
    J[478] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 78: CH2O + CH2(S) <=> CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[16];
    k_f = prefactor_units[77] * fwd_A[77]
                * exp(fwd_beta[77] * tc[0] - activation_units[77] * fwd_Ea[77] * invT);
    dlnkfdT = fwd_beta[77] * invT + activation_units[77] * fwd_Ea[77] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[15];
    Kc = exp(g_RT[12] - g_RT[13] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[16]) + (h_RT[13] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[12] -= q; /* CH2(S) */
    wdot[13] += q; /* CH3 */
    wdot[15] += q; /* HCO */
    wdot[16] -= q; /* CH2O */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[16];
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[277] += dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    J[279] += dqdci;              /* dwdot[HCO]/d[CH2(S)] */
    J[280] -= dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[15];
    J[298] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[301] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[302] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[13];
    J[342] -= dqdci;              /* dwdot[CH2(S)]/d[HCO] */
    J[343] += dqdci;              /* dwdot[CH3]/d[HCO] */
    J[345] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[346] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[12];
    J[364] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[365] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[367] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[368] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[477] += dqdT;               /* dwdot[HCO]/dT */
    J[478] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 79: CH3 + O <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[13];
    k_f = prefactor_units[78] * fwd_A[78]
                * exp(fwd_beta[78] * tc[0] - activation_units[78] * fwd_Ea[78] * invT);
    dlnkfdT = fwd_beta[78] * invT + activation_units[78] * fwd_Ea[78] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[16];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[13] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[13]) + (h_RT[2] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[13] -= q; /* CH3 */
    wdot[16] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[46] += dqdci;               /* dwdot[H]/d[H] */
    J[47] -= dqdci;               /* dwdot[O]/d[H] */
    J[57] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[60] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[68] += dqdci;               /* dwdot[H]/d[O] */
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[79] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[82] += dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[288] += dqdci;              /* dwdot[H]/d[CH3] */
    J[289] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[302] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[2];
    J[354] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[355] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[365] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 80: CH3 + O => H + H2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[13];
    k_f = prefactor_units[79] * fwd_A[79]
                * exp(fwd_beta[79] * tc[0] - activation_units[79] * fwd_Ea[79] * invT);
    dlnkfdT = fwd_beta[79] * invT + activation_units[79] * fwd_Ea[79] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[8] += q; /* CO */
    wdot[13] -= q; /* CH3 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[67] += dqdci;               /* dwdot[H2]/d[O] */
    J[68] += dqdci;               /* dwdot[H]/d[O] */
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[74] += dqdci;               /* dwdot[CO]/d[O] */
    J[79] -= dqdci;               /* dwdot[CH3]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[287] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[288] += dqdci;              /* dwdot[H]/d[CH3] */
    J[289] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[294] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[464] += dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 81: CH3 + OH <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = prefactor_units[80] * fwd_A[80]
                * exp(fwd_beta[80] * tc[0] - activation_units[80] * fwd_Ea[80] * invT);
    dlnkfdT = fwd_beta[80] * invT + activation_units[80] * fwd_Ea[80] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[11];
    Kc = exp(g_RT[5] - g_RT[6] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[13]) + (h_RT[6] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[11] += q; /* CH2 */
    wdot[13] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[121] += dqdci;              /* dwdot[CH2]/d[OH] */
    J[123] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[143] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[145] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[6];
    J[247] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[248] += dqdci;              /* dwdot[H2O]/d[CH2] */
    J[253] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[255] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[5];
    J[291] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[292] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[297] += dqdci;              /* dwdot[CH2]/d[CH3] */
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[473] += dqdT;               /* dwdot[CH2]/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 82: CH3 + OH <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = prefactor_units[81] * fwd_A[81]
                * exp(fwd_beta[81] * tc[0] - activation_units[81] * fwd_Ea[81] * invT);
    dlnkfdT = fwd_beta[81] * invT + activation_units[81] * fwd_Ea[81] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[12];
    Kc = exp(g_RT[5] - g_RT[6] - g_RT[12] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[13]) + (h_RT[6] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[12] += q; /* CH2(S) */
    wdot[13] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[122] += dqdci;              /* dwdot[CH2(S)]/d[OH] */
    J[123] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[12];
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[144] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[145] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[6];
    J[269] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[270] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[276] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[277] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[5];
    J[291] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[292] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[298] += dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[474] += dqdT;               /* dwdot[CH2(S)]/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 83: CH3 + OH <=> H2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = prefactor_units[82] * fwd_A[82]
                * exp(fwd_beta[82] * tc[0] - activation_units[82] * fwd_Ea[82] * invT);
    dlnkfdT = fwd_beta[82] * invT + activation_units[82] * fwd_Ea[82] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[5] + g_RT[13] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[13]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[5] -= q; /* OH */
    wdot[13] -= q; /* CH3 */
    wdot[16] += q; /* CH2O */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[16];
    J[23] += dqdci;               /* dwdot[H2]/d[H2] */
    J[27] -= dqdci;               /* dwdot[OH]/d[H2] */
    J[35] -= dqdci;               /* dwdot[CH3]/d[H2] */
    J[38] += dqdci;               /* dwdot[CH2O]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[111] += dqdci;              /* dwdot[H2]/d[OH] */
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[123] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[126] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[5];
    J[287] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[291] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[302] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[353] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[357] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[365] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 84: CH3 + HO2 <=> O2 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[13];
    k_f = prefactor_units[83] * fwd_A[83]
                * exp(fwd_beta[83] * tc[0] - activation_units[83] * fwd_Ea[83] * invT);
    dlnkfdT = fwd_beta[83] * invT + activation_units[83] * fwd_Ea[83] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(-g_RT[4] + g_RT[7] + g_RT[13] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[13]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* O2 */
    wdot[7] -= q; /* HO2 */
    wdot[13] -= q; /* CH3 */
    wdot[14] += q; /* CH4 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[14];
    J[92] += dqdci;               /* dwdot[O2]/d[O2] */
    J[95] -= dqdci;               /* dwdot[HO2]/d[O2] */
    J[101] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[102] += dqdci;              /* dwdot[CH4]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[13];
    J[158] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[161] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[167] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[168] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[290] += dqdci;              /* dwdot[O2]/d[CH3] */
    J[293] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[300] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[4];
    J[312] += dqdci;              /* dwdot[O2]/d[CH4] */
    J[315] -= dqdci;              /* dwdot[HO2]/d[CH4] */
    J[321] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[322] += dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[O2]/dT */
    J[469] -= dqdT;               /* dwdot[HO2]/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH4]/dT */

    /*reaction 85: CH3 + HO2 <=> OH + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[13];
    k_f = prefactor_units[84] * fwd_A[84]
                * exp(fwd_beta[84] * tc[0] - activation_units[84] * fwd_Ea[84] * invT);
    dlnkfdT = fwd_beta[84] * invT + activation_units[84] * fwd_Ea[84] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[18];
    Kc = exp(-g_RT[5] + g_RT[7] + g_RT[13] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[13]) + (h_RT[5] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* OH */
    wdot[7] -= q; /* HO2 */
    wdot[13] -= q; /* CH3 */
    wdot[18] += q; /* CH3O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[18];
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[117] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[123] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[128] += dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[13];
    J[159] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[161] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[167] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[172] += dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[291] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[293] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[304] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[5];
    J[401] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[403] -= dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[409] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[414] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[HO2]/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */
    J[480] += dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 86: CH3 + O2 <=> O + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[13];
    k_f = prefactor_units[85] * fwd_A[85]
                * exp(fwd_beta[85] * tc[0] - activation_units[85] * fwd_Ea[85] * invT);
    dlnkfdT = fwd_beta[85] * invT + activation_units[85] * fwd_Ea[85] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[18];
    Kc = exp(-g_RT[3] + g_RT[4] + g_RT[13] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[3] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O */
    wdot[4] -= q; /* O2 */
    wdot[13] -= q; /* CH3 */
    wdot[18] += q; /* CH3O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[18];
    J[69] += dqdci;               /* dwdot[O]/d[O] */
    J[70] -= dqdci;               /* dwdot[O2]/d[O] */
    J[79] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[84] += dqdci;               /* dwdot[CH3O]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13];
    J[91] += dqdci;               /* dwdot[O]/d[O2] */
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[101] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[106] += dqdci;              /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[289] += dqdci;              /* dwdot[O]/d[CH3] */
    J[290] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[304] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[3];
    J[399] += dqdci;              /* dwdot[O]/d[CH3O] */
    J[400] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[409] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[414] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O]/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */
    J[480] += dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 87: CH3 + O2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[13];
    k_f = prefactor_units[86] * fwd_A[86]
                * exp(fwd_beta[86] * tc[0] - activation_units[86] * fwd_Ea[86] * invT);
    dlnkfdT = fwd_beta[86] * invT + activation_units[86] * fwd_Ea[86] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[16];
    Kc = exp(g_RT[4] - g_RT[5] + g_RT[13] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[5] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O2 */
    wdot[5] += q; /* OH */
    wdot[13] -= q; /* CH3 */
    wdot[16] += q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13];
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[93] += dqdci;               /* dwdot[OH]/d[O2] */
    J[101] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[104] += dqdci;              /* dwdot[CH2O]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[16];
    J[114] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[123] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[126] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[290] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[291] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[302] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[5];
    J[356] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[357] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[365] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 88: CH3 + HCO <=> CH4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[15];
    k_f = prefactor_units[87] * fwd_A[87]
                * exp(fwd_beta[87] * tc[0] - activation_units[87] * fwd_Ea[87] * invT);
    dlnkfdT = fwd_beta[87] * invT + activation_units[87] * fwd_Ea[87] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[14];
    Kc = exp(-g_RT[8] + g_RT[13] - g_RT[14] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[15]) + (h_RT[8] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] += q; /* CO */
    wdot[13] -= q; /* CH3 */
    wdot[14] += q; /* CH4 */
    wdot[15] -= q; /* HCO */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[14];
    J[184] += dqdci;              /* dwdot[CO]/d[CO] */
    J[189] -= dqdci;              /* dwdot[CH3]/d[CO] */
    J[190] += dqdci;              /* dwdot[CH4]/d[CO] */
    J[191] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[15];
    J[294] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[300] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[301] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[8];
    J[316] += dqdci;              /* dwdot[CO]/d[CH4] */
    J[321] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[322] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[323] -= dqdci;              /* dwdot[HCO]/d[CH4] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[13];
    J[338] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[343] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[344] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[345] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[470] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH4]/dT */
    J[477] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 89: CH3 + CH2O <=> HCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[16];
    k_f = prefactor_units[88] * fwd_A[88]
                * exp(fwd_beta[88] * tc[0] - activation_units[88] * fwd_Ea[88] * invT);
    dlnkfdT = fwd_beta[88] * invT + activation_units[88] * fwd_Ea[88] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[15];
    Kc = exp(g_RT[13] - g_RT[14] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[16]) + (h_RT[14] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[13] -= q; /* CH3 */
    wdot[14] += q; /* CH4 */
    wdot[15] += q; /* HCO */
    wdot[16] -= q; /* CH2O */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[300] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[301] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[302] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[15];
    J[321] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[322] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[323] += dqdci;              /* dwdot[HCO]/d[CH4] */
    J[324] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[14];
    J[343] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[344] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[345] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[346] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[13];
    J[365] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[366] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[367] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[368] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH4]/dT */
    J[477] += dqdT;               /* dwdot[HCO]/dT */
    J[478] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 90: CH3O + H <=> H + CH2OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = prefactor_units[89] * fwd_A[89]
                * exp(fwd_beta[89] * tc[0] - activation_units[89] * fwd_Ea[89] * invT);
    dlnkfdT = fwd_beta[89] * invT + activation_units[89] * fwd_Ea[89] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[17];
    Kc = exp(g_RT[2] - g_RT[2] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[18]) + (h_RT[2] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[17] += q; /* CH2OH */
    wdot[18] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18] - k_r*sc[17];
    J[61] += dqdci;               /* dwdot[CH2OH]/d[H] */
    J[62] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[2];
    J[391] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[392] -= dqdci;              /* dwdot[CH3O]/d[CH2OH] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[2];
    J[413] += dqdci;              /* dwdot[CH2OH]/d[CH3O] */
    J[414] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[479] += dqdT;               /* dwdot[CH2OH]/dT */
    J[480] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 91: CH3O + H <=> H2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = prefactor_units[90] * fwd_A[90]
                * exp(fwd_beta[90] * tc[0] - activation_units[90] * fwd_Ea[90] * invT);
    dlnkfdT = fwd_beta[90] * invT + activation_units[90] * fwd_Ea[90] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[16] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[18]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[2] -= q; /* H */
    wdot[16] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[16];
    J[23] += dqdci;               /* dwdot[H2]/d[H2] */
    J[24] -= dqdci;               /* dwdot[H]/d[H2] */
    J[38] += dqdci;               /* dwdot[CH2O]/d[H2] */
    J[40] -= dqdci;               /* dwdot[CH3O]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[45] += dqdci;               /* dwdot[H2]/d[H] */
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[60] += dqdci;               /* dwdot[CH2O]/d[H] */
    J[62] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[353] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[354] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[370] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[2];
    J[397] += dqdci;              /* dwdot[H2]/d[CH3O] */
    J[398] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[412] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[414] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */
    J[480] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 92: CH3O + H <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = prefactor_units[91] * fwd_A[91]
                * exp(fwd_beta[91] * tc[0] - activation_units[91] * fwd_Ea[91] * invT);
    dlnkfdT = fwd_beta[91] * invT + activation_units[91] * fwd_Ea[91] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[13];
    Kc = exp(g_RT[2] - g_RT[5] - g_RT[13] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[18]) + (h_RT[5] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[5] += q; /* OH */
    wdot[13] += q; /* CH3 */
    wdot[18] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[49] += dqdci;               /* dwdot[OH]/d[H] */
    J[57] += dqdci;               /* dwdot[CH3]/d[H] */
    J[62] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[112] -= dqdci;              /* dwdot[H]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[123] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[128] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[5];
    J[288] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[291] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[304] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[2];
    J[398] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[401] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[409] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[414] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[480] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 93: CH3O + H <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = prefactor_units[92] * fwd_A[92]
                * exp(fwd_beta[92] * tc[0] - activation_units[92] * fwd_Ea[92] * invT);
    dlnkfdT = fwd_beta[92] * invT + activation_units[92] * fwd_Ea[92] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[12];
    Kc = exp(g_RT[2] - g_RT[6] - g_RT[12] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[18]) + (h_RT[6] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[6] += q; /* H2O */
    wdot[12] += q; /* CH2(S) */
    wdot[18] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[50] += dqdci;               /* dwdot[H2O]/d[H] */
    J[56] += dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[62] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[12];
    J[134] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[144] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[150] -= dqdci;              /* dwdot[CH3O]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[6];
    J[266] -= dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[270] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[276] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[282] -= dqdci;              /* dwdot[CH3O]/d[CH2(S)] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[2];
    J[398] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[402] += dqdci;              /* dwdot[H2O]/d[CH3O] */
    J[408] += dqdci;              /* dwdot[CH2(S)]/d[CH3O] */
    J[414] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[474] += dqdT;               /* dwdot[CH2(S)]/dT */
    J[480] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 94: CH3O + O <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[18];
    k_f = prefactor_units[93] * fwd_A[93]
                * exp(fwd_beta[93] * tc[0] - activation_units[93] * fwd_Ea[93] * invT);
    dlnkfdT = fwd_beta[93] * invT + activation_units[93] * fwd_Ea[93] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[16];
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[16] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[18]) + (h_RT[5] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[16] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[18];
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[71] += dqdci;               /* dwdot[OH]/d[O] */
    J[82] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[84] -= dqdci;               /* dwdot[CH3O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[16];
    J[113] -= dqdci;              /* dwdot[O]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[126] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[128] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[5];
    J[355] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[357] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[370] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[3];
    J[399] -= dqdci;              /* dwdot[O]/d[CH3O] */
    J[401] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[412] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[414] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */
    J[480] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 95: CH3O + OH <=> H2O + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[18];
    k_f = prefactor_units[94] * fwd_A[94]
                * exp(fwd_beta[94] * tc[0] - activation_units[94] * fwd_Ea[94] * invT);
    dlnkfdT = fwd_beta[94] * invT + activation_units[94] * fwd_Ea[94] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[16];
    Kc = exp(g_RT[5] - g_RT[6] - g_RT[16] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[18]) + (h_RT[6] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[16] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[18];
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[126] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[128] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[16];
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[148] += dqdci;              /* dwdot[CH2O]/d[H2O] */
    J[150] -= dqdci;              /* dwdot[CH3O]/d[H2O] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[6];
    J[357] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[358] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[370] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[5];
    J[401] -= dqdci;              /* dwdot[OH]/d[CH3O] */
    J[402] += dqdci;              /* dwdot[H2O]/d[CH3O] */
    J[412] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[414] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */
    J[480] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 96: CH3O + O2 <=> HO2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[18];
    k_f = prefactor_units[95] * fwd_A[95]
                * exp(fwd_beta[95] * tc[0] - activation_units[95] * fwd_Ea[95] * invT);
    dlnkfdT = fwd_beta[95] * invT + activation_units[95] * fwd_Ea[95] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[16];
    Kc = exp(g_RT[4] - g_RT[7] - g_RT[16] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[18]) + (h_RT[7] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O2 */
    wdot[7] += q; /* HO2 */
    wdot[16] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[18];
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[95] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[104] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[106] -= dqdci;              /* dwdot[CH3O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[16];
    J[158] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[161] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[170] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[172] -= dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[7];
    J[356] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[359] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[370] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[4];
    J[400] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[403] += dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[412] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[414] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[469] += dqdT;               /* dwdot[HO2]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */
    J[480] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 97: CH3O + CH3 <=> CH4 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[18];
    k_f = prefactor_units[96] * fwd_A[96]
                * exp(fwd_beta[96] * tc[0] - activation_units[96] * fwd_Ea[96] * invT);
    dlnkfdT = fwd_beta[96] * invT + activation_units[96] * fwd_Ea[96] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[16];
    Kc = exp(g_RT[13] - g_RT[14] - g_RT[16] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[18]) + (h_RT[14] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[13] -= q; /* CH3 */
    wdot[14] += q; /* CH4 */
    wdot[16] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[18];
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[300] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[302] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[304] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[16];
    J[321] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[322] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[324] += dqdci;              /* dwdot[CH2O]/d[CH4] */
    J[326] -= dqdci;              /* dwdot[CH3O]/d[CH4] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[14];
    J[365] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[366] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[370] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[13];
    J[409] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[410] += dqdci;              /* dwdot[CH4]/d[CH3O] */
    J[412] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[414] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH4]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */
    J[480] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 98: CH3O + CO <=> CH3 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[18];
    k_f = prefactor_units[97] * fwd_A[97]
                * exp(fwd_beta[97] * tc[0] - activation_units[97] * fwd_Ea[97] * invT);
    dlnkfdT = fwd_beta[97] * invT + activation_units[97] * fwd_Ea[97] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[13];
    Kc = exp(g_RT[8] - g_RT[9] - g_RT[13] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[18]) + (h_RT[9] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= q; /* CO */
    wdot[9] += q; /* CO2 */
    wdot[13] += q; /* CH3 */
    wdot[18] -= q; /* CH3O */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[18];
    J[184] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[185] += dqdci;              /* dwdot[CO2]/d[CO] */
    J[189] += dqdci;              /* dwdot[CH3]/d[CO] */
    J[194] -= dqdci;              /* dwdot[CH3O]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[13];
    J[206] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[207] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[211] += dqdci;              /* dwdot[CH3]/d[CO2] */
    J[216] -= dqdci;              /* dwdot[CH3O]/d[CO2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[9];
    J[294] -= dqdci;              /* dwdot[CO]/d[CH3] */
    J[295] += dqdci;              /* dwdot[CO2]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[304] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[8];
    J[404] -= dqdci;              /* dwdot[CO]/d[CH3O] */
    J[405] += dqdci;              /* dwdot[CO2]/d[CH3O] */
    J[409] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[414] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[470] -= dqdT;               /* dwdot[CO]/dT */
    J[471] += dqdT;               /* dwdot[CO2]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[480] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 99: CH2OH + H <=> H2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = prefactor_units[98] * fwd_A[98]
                * exp(fwd_beta[98] * tc[0] - activation_units[98] * fwd_Ea[98] * invT);
    dlnkfdT = fwd_beta[98] * invT + activation_units[98] * fwd_Ea[98] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[2] -= q; /* H */
    wdot[16] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[16];
    J[23] += dqdci;               /* dwdot[H2]/d[H2] */
    J[24] -= dqdci;               /* dwdot[H]/d[H2] */
    J[38] += dqdci;               /* dwdot[CH2O]/d[H2] */
    J[39] -= dqdci;               /* dwdot[CH2OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[17];
    J[45] += dqdci;               /* dwdot[H2]/d[H] */
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[60] += dqdci;               /* dwdot[CH2O]/d[H] */
    J[61] -= dqdci;               /* dwdot[CH2OH]/d[H] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[353] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[354] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[2];
    J[375] += dqdci;              /* dwdot[H2]/d[CH2OH] */
    J[376] -= dqdci;              /* dwdot[H]/d[CH2OH] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[391] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */
    J[479] -= dqdT;               /* dwdot[CH2OH]/dT */

    /*reaction 100: CH2OH + H <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = prefactor_units[99] * fwd_A[99]
                * exp(fwd_beta[99] * tc[0] - activation_units[99] * fwd_Ea[99] * invT);
    dlnkfdT = fwd_beta[99] * invT + activation_units[99] * fwd_Ea[99] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[13];
    Kc = exp(g_RT[2] - g_RT[5] - g_RT[13] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[5] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[5] += q; /* OH */
    wdot[13] += q; /* CH3 */
    wdot[17] -= q; /* CH2OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[17];
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[49] += dqdci;               /* dwdot[OH]/d[H] */
    J[57] += dqdci;               /* dwdot[CH3]/d[H] */
    J[61] -= dqdci;               /* dwdot[CH2OH]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[112] -= dqdci;              /* dwdot[H]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[123] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[127] -= dqdci;              /* dwdot[CH2OH]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[5];
    J[288] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[291] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[303] -= dqdci;              /* dwdot[CH2OH]/d[CH3] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[2];
    J[376] -= dqdci;              /* dwdot[H]/d[CH2OH] */
    J[379] += dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[387] += dqdci;              /* dwdot[CH3]/d[CH2OH] */
    J[391] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[479] -= dqdT;               /* dwdot[CH2OH]/dT */

    /*reaction 101: CH2OH + H <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = prefactor_units[100] * fwd_A[100]
                * exp(fwd_beta[100] * tc[0] - activation_units[100] * fwd_Ea[100] * invT);
    dlnkfdT = fwd_beta[100] * invT + activation_units[100] * fwd_Ea[100] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[12];
    Kc = exp(g_RT[2] - g_RT[6] - g_RT[12] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[6] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[6] += q; /* H2O */
    wdot[12] += q; /* CH2(S) */
    wdot[17] -= q; /* CH2OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[17];
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[50] += dqdci;               /* dwdot[H2O]/d[H] */
    J[56] += dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[61] -= dqdci;               /* dwdot[CH2OH]/d[H] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[12];
    J[134] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[144] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[149] -= dqdci;              /* dwdot[CH2OH]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[6];
    J[266] -= dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[270] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[276] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[281] -= dqdci;              /* dwdot[CH2OH]/d[CH2(S)] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[2];
    J[376] -= dqdci;              /* dwdot[H]/d[CH2OH] */
    J[380] += dqdci;              /* dwdot[H2O]/d[CH2OH] */
    J[386] += dqdci;              /* dwdot[CH2(S)]/d[CH2OH] */
    J[391] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[474] += dqdT;               /* dwdot[CH2(S)]/dT */
    J[479] -= dqdT;               /* dwdot[CH2OH]/dT */

    /*reaction 102: CH2OH + O <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[17];
    k_f = prefactor_units[101] * fwd_A[101]
                * exp(fwd_beta[101] * tc[0] - activation_units[101] * fwd_Ea[101] * invT);
    dlnkfdT = fwd_beta[101] * invT + activation_units[101] * fwd_Ea[101] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[16];
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[17]) + (h_RT[5] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[16] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[71] += dqdci;               /* dwdot[OH]/d[O] */
    J[82] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[83] -= dqdci;               /* dwdot[CH2OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[16];
    J[113] -= dqdci;              /* dwdot[O]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[126] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[127] -= dqdci;              /* dwdot[CH2OH]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[5];
    J[355] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[357] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[3];
    J[377] -= dqdci;              /* dwdot[O]/d[CH2OH] */
    J[379] += dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[391] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */
    J[479] -= dqdT;               /* dwdot[CH2OH]/dT */

    /*reaction 103: CH2OH + OH <=> H2O + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[17];
    k_f = prefactor_units[102] * fwd_A[102]
                * exp(fwd_beta[102] * tc[0] - activation_units[102] * fwd_Ea[102] * invT);
    dlnkfdT = fwd_beta[102] * invT + activation_units[102] * fwd_Ea[102] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[16];
    Kc = exp(g_RT[5] - g_RT[6] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[17]) + (h_RT[6] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[16] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[17];
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[126] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[127] -= dqdci;              /* dwdot[CH2OH]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[16];
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[148] += dqdci;              /* dwdot[CH2O]/d[H2O] */
    J[149] -= dqdci;              /* dwdot[CH2OH]/d[H2O] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[6];
    J[357] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[358] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[5];
    J[379] -= dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[380] += dqdci;              /* dwdot[H2O]/d[CH2OH] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[391] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */
    J[479] -= dqdT;               /* dwdot[CH2OH]/dT */

    /*reaction 104: CH2OH + O2 <=> HO2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[17];
    k_f = prefactor_units[103] * fwd_A[103]
                * exp(fwd_beta[103] * tc[0] - activation_units[103] * fwd_Ea[103] * invT);
    dlnkfdT = fwd_beta[103] * invT + activation_units[103] * fwd_Ea[103] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[16];
    Kc = exp(g_RT[4] - g_RT[7] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[17]) + (h_RT[7] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O2 */
    wdot[7] += q; /* HO2 */
    wdot[16] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[17];
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[95] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[104] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[105] -= dqdci;              /* dwdot[CH2OH]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[16];
    J[158] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[161] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[170] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[171] -= dqdci;              /* dwdot[CH2OH]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[7];
    J[356] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[359] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[4];
    J[378] -= dqdci;              /* dwdot[O2]/d[CH2OH] */
    J[381] += dqdci;              /* dwdot[HO2]/d[CH2OH] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[391] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[469] += dqdT;               /* dwdot[HO2]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */
    J[479] -= dqdT;               /* dwdot[CH2OH]/dT */

    /*reaction 105: CH2OH + CH3 <=> CH4 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[17];
    k_f = prefactor_units[104] * fwd_A[104]
                * exp(fwd_beta[104] * tc[0] - activation_units[104] * fwd_Ea[104] * invT);
    dlnkfdT = fwd_beta[104] * invT + activation_units[104] * fwd_Ea[104] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[16];
    Kc = exp(g_RT[13] - g_RT[14] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[17]) + (h_RT[14] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[13] -= q; /* CH3 */
    wdot[14] += q; /* CH4 */
    wdot[16] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[17];
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[300] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[302] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[303] -= dqdci;              /* dwdot[CH2OH]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[16];
    J[321] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[322] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[324] += dqdci;              /* dwdot[CH2O]/d[CH4] */
    J[325] -= dqdci;              /* dwdot[CH2OH]/d[CH4] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[14];
    J[365] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[366] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[13];
    J[387] -= dqdci;              /* dwdot[CH3]/d[CH2OH] */
    J[388] += dqdci;              /* dwdot[CH4]/d[CH2OH] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[391] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH4]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */
    J[479] -= dqdT;               /* dwdot[CH2OH]/dT */

    /*reaction 106: CH4 + H <=> CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[14];
    k_f = prefactor_units[105] * fwd_A[105]
                * exp(fwd_beta[105] * tc[0] - activation_units[105] * fwd_Ea[105] * invT);
    dlnkfdT = fwd_beta[105] * invT + activation_units[105] * fwd_Ea[105] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[14]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[2] -= q; /* H */
    wdot[13] += q; /* CH3 */
    wdot[14] -= q; /* CH4 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[13];
    J[23] += dqdci;               /* dwdot[H2]/d[H2] */
    J[24] -= dqdci;               /* dwdot[H]/d[H2] */
    J[35] += dqdci;               /* dwdot[CH3]/d[H2] */
    J[36] -= dqdci;               /* dwdot[CH4]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[14];
    J[45] += dqdci;               /* dwdot[H2]/d[H] */
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[57] += dqdci;               /* dwdot[CH3]/d[H] */
    J[58] -= dqdci;               /* dwdot[CH4]/d[H] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[287] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[288] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[300] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[309] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[310] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[321] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[322] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[476] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 107: CH4 + O <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[14];
    k_f = prefactor_units[106] * fwd_A[106]
                * exp(fwd_beta[106] * tc[0] - activation_units[106] * fwd_Ea[106] * invT);
    dlnkfdT = fwd_beta[106] * invT + activation_units[106] * fwd_Ea[106] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[13];
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[14]) + (h_RT[5] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[13] += q; /* CH3 */
    wdot[14] -= q; /* CH4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[14];
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[71] += dqdci;               /* dwdot[OH]/d[O] */
    J[79] += dqdci;               /* dwdot[CH3]/d[O] */
    J[80] -= dqdci;               /* dwdot[CH4]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[113] -= dqdci;              /* dwdot[O]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[123] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[124] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[5];
    J[289] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[291] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[300] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[3];
    J[311] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[313] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[321] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[322] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[476] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 108: CH4 + OH <=> CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[14];
    k_f = prefactor_units[107] * fwd_A[107]
                * exp(fwd_beta[107] * tc[0] - activation_units[107] * fwd_Ea[107] * invT);
    dlnkfdT = fwd_beta[107] * invT + activation_units[107] * fwd_Ea[107] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[13];
    Kc = exp(g_RT[5] - g_RT[6] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[14]) + (h_RT[6] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[13] += q; /* CH3 */
    wdot[14] -= q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[14];
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[123] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[124] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[13];
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[145] += dqdci;              /* dwdot[CH3]/d[H2O] */
    J[146] -= dqdci;              /* dwdot[CH4]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[6];
    J[291] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[292] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[300] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[5];
    J[313] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[314] += dqdci;              /* dwdot[H2O]/d[CH4] */
    J[321] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[322] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[476] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 109: CH4 + CH2 <=> 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[14];
    k_f = prefactor_units[108] * fwd_A[108]
                * exp(fwd_beta[108] * tc[0] - activation_units[108] * fwd_Ea[108] * invT);
    dlnkfdT = fwd_beta[108] * invT + activation_units[108] * fwd_Ea[108] * invT2;
    /* reverse */
    phi_r = pow(sc[13], 2.000000);
    Kc = exp(g_RT[11] - 2.000000*g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[14]) + (2.000000*h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] -= q; /* CH2 */
    wdot[13] += 2 * q; /* CH3 */
    wdot[14] -= q; /* CH4 */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[14];
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[255] += 2 * dqdci;          /* dwdot[CH3]/d[CH2] */
    J[256] -= dqdci;              /* dwdot[CH4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[13];
    J[297] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[299] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[300] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[11];
    J[319] -= dqdci;              /* dwdot[CH2]/d[CH4] */
    J[321] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[322] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[476] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 110: CH4 + CH2(S) <=> 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[14];
    k_f = prefactor_units[109] * fwd_A[109]
                * exp(fwd_beta[109] * tc[0] - activation_units[109] * fwd_Ea[109] * invT);
    dlnkfdT = fwd_beta[109] * invT + activation_units[109] * fwd_Ea[109] * invT2;
    /* reverse */
    phi_r = pow(sc[13], 2.000000);
    Kc = exp(g_RT[12] - 2.000000*g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[14]) + (2.000000*h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[12] -= q; /* CH2(S) */
    wdot[13] += 2 * q; /* CH3 */
    wdot[14] -= q; /* CH4 */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[14];
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[277] += 2 * dqdci;          /* dwdot[CH3]/d[CH2(S)] */
    J[278] -= dqdci;              /* dwdot[CH4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[13];
    J[298] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[299] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[300] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[12];
    J[320] -= dqdci;              /* dwdot[CH2(S)]/d[CH4] */
    J[321] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[322] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[475] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[476] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 111: CH3OH + H <=> CH2OH + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[19];
    k_f = prefactor_units[110] * fwd_A[110]
                * exp(fwd_beta[110] * tc[0] - activation_units[110] * fwd_Ea[110] * invT);
    dlnkfdT = fwd_beta[110] * invT + activation_units[110] * fwd_Ea[110] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[17];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[17] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[19]) + (h_RT[1] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[2] -= q; /* H */
    wdot[17] += q; /* CH2OH */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[17];
    J[23] += dqdci;               /* dwdot[H2]/d[H2] */
    J[24] -= dqdci;               /* dwdot[H]/d[H2] */
    J[39] += dqdci;               /* dwdot[CH2OH]/d[H2] */
    J[41] -= dqdci;               /* dwdot[CH3OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[19];
    J[45] += dqdci;               /* dwdot[H2]/d[H] */
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[61] += dqdci;               /* dwdot[CH2OH]/d[H] */
    J[63] -= dqdci;               /* dwdot[CH3OH]/d[H] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[1];
    J[375] += dqdci;              /* dwdot[H2]/d[CH2OH] */
    J[376] -= dqdci;              /* dwdot[H]/d[CH2OH] */
    J[391] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[393] -= dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[2];
    J[419] += dqdci;              /* dwdot[H2]/d[CH3OH] */
    J[420] -= dqdci;              /* dwdot[H]/d[CH3OH] */
    J[435] += dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[479] += dqdT;               /* dwdot[CH2OH]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 112: CH3OH + H <=> CH3O + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[19];
    k_f = prefactor_units[111] * fwd_A[111]
                * exp(fwd_beta[111] * tc[0] - activation_units[111] * fwd_Ea[111] * invT);
    dlnkfdT = fwd_beta[111] * invT + activation_units[111] * fwd_Ea[111] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[18];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[19]) + (h_RT[1] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[2] -= q; /* H */
    wdot[18] += q; /* CH3O */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[18];
    J[23] += dqdci;               /* dwdot[H2]/d[H2] */
    J[24] -= dqdci;               /* dwdot[H]/d[H2] */
    J[40] += dqdci;               /* dwdot[CH3O]/d[H2] */
    J[41] -= dqdci;               /* dwdot[CH3OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[19];
    J[45] += dqdci;               /* dwdot[H2]/d[H] */
    J[46] -= dqdci;               /* dwdot[H]/d[H] */
    J[62] += dqdci;               /* dwdot[CH3O]/d[H] */
    J[63] -= dqdci;               /* dwdot[CH3OH]/d[H] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[1];
    J[397] += dqdci;              /* dwdot[H2]/d[CH3O] */
    J[398] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[414] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[415] -= dqdci;              /* dwdot[CH3OH]/d[CH3O] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[2];
    J[419] += dqdci;              /* dwdot[H2]/d[CH3OH] */
    J[420] -= dqdci;              /* dwdot[H]/d[CH3OH] */
    J[436] += dqdci;              /* dwdot[CH3O]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H2]/dT */
    J[464] -= dqdT;               /* dwdot[H]/dT */
    J[480] += dqdT;               /* dwdot[CH3O]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 113: CH3OH + O <=> OH + CH2OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[19];
    k_f = prefactor_units[112] * fwd_A[112]
                * exp(fwd_beta[112] * tc[0] - activation_units[112] * fwd_Ea[112] * invT);
    dlnkfdT = fwd_beta[112] * invT + activation_units[112] * fwd_Ea[112] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[17];
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[17] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[19]) + (h_RT[5] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[17] += q; /* CH2OH */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[O] */
    dqdci =  + k_f*sc[19];
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[71] += dqdci;               /* dwdot[OH]/d[O] */
    J[83] += dqdci;               /* dwdot[CH2OH]/d[O] */
    J[85] -= dqdci;               /* dwdot[CH3OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[17];
    J[113] -= dqdci;              /* dwdot[O]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[CH2OH]/d[OH] */
    J[129] -= dqdci;              /* dwdot[CH3OH]/d[OH] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[5];
    J[377] -= dqdci;              /* dwdot[O]/d[CH2OH] */
    J[379] += dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[391] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[393] -= dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[3];
    J[421] -= dqdci;              /* dwdot[O]/d[CH3OH] */
    J[423] += dqdci;              /* dwdot[OH]/d[CH3OH] */
    J[435] += dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[479] += dqdT;               /* dwdot[CH2OH]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 114: CH3OH + O <=> OH + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[19];
    k_f = prefactor_units[113] * fwd_A[113]
                * exp(fwd_beta[113] * tc[0] - activation_units[113] * fwd_Ea[113] * invT);
    dlnkfdT = fwd_beta[113] * invT + activation_units[113] * fwd_Ea[113] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[18];
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[19]) + (h_RT[5] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[18] += q; /* CH3O */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[O] */
    dqdci =  + k_f*sc[19];
    J[69] -= dqdci;               /* dwdot[O]/d[O] */
    J[71] += dqdci;               /* dwdot[OH]/d[O] */
    J[84] += dqdci;               /* dwdot[CH3O]/d[O] */
    J[85] -= dqdci;               /* dwdot[CH3OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[18];
    J[113] -= dqdci;              /* dwdot[O]/d[OH] */
    J[115] += dqdci;              /* dwdot[OH]/d[OH] */
    J[128] += dqdci;              /* dwdot[CH3O]/d[OH] */
    J[129] -= dqdci;              /* dwdot[CH3OH]/d[OH] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[5];
    J[399] -= dqdci;              /* dwdot[O]/d[CH3O] */
    J[401] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[414] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[415] -= dqdci;              /* dwdot[CH3OH]/d[CH3O] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[3];
    J[421] -= dqdci;              /* dwdot[O]/d[CH3OH] */
    J[423] += dqdci;              /* dwdot[OH]/d[CH3OH] */
    J[436] += dqdci;              /* dwdot[CH3O]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O]/dT */
    J[467] += dqdT;               /* dwdot[OH]/dT */
    J[480] += dqdT;               /* dwdot[CH3O]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 115: CH3OH + OH <=> CH2OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[19];
    k_f = prefactor_units[114] * fwd_A[114]
                * exp(fwd_beta[114] * tc[0] - activation_units[114] * fwd_Ea[114] * invT);
    dlnkfdT = fwd_beta[114] * invT + activation_units[114] * fwd_Ea[114] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[17];
    Kc = exp(g_RT[5] - g_RT[6] - g_RT[17] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[19]) + (h_RT[6] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[17] += q; /* CH2OH */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[19];
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[127] += dqdci;              /* dwdot[CH2OH]/d[OH] */
    J[129] -= dqdci;              /* dwdot[CH3OH]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[17];
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[149] += dqdci;              /* dwdot[CH2OH]/d[H2O] */
    J[151] -= dqdci;              /* dwdot[CH3OH]/d[H2O] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[6];
    J[379] -= dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[380] += dqdci;              /* dwdot[H2O]/d[CH2OH] */
    J[391] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[393] -= dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[5];
    J[423] -= dqdci;              /* dwdot[OH]/d[CH3OH] */
    J[424] += dqdci;              /* dwdot[H2O]/d[CH3OH] */
    J[435] += dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[479] += dqdT;               /* dwdot[CH2OH]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 116: CH3OH + OH <=> CH3O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[19];
    k_f = prefactor_units[115] * fwd_A[115]
                * exp(fwd_beta[115] * tc[0] - activation_units[115] * fwd_Ea[115] * invT);
    dlnkfdT = fwd_beta[115] * invT + activation_units[115] * fwd_Ea[115] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[18];
    Kc = exp(g_RT[5] - g_RT[6] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[19]) + (h_RT[6] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[18] += q; /* CH3O */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[19];
    J[115] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[128] += dqdci;              /* dwdot[CH3O]/d[OH] */
    J[129] -= dqdci;              /* dwdot[CH3OH]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[18];
    J[137] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[138] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[150] += dqdci;              /* dwdot[CH3O]/d[H2O] */
    J[151] -= dqdci;              /* dwdot[CH3OH]/d[H2O] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[6];
    J[401] -= dqdci;              /* dwdot[OH]/d[CH3O] */
    J[402] += dqdci;              /* dwdot[H2O]/d[CH3O] */
    J[414] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[415] -= dqdci;              /* dwdot[CH3OH]/d[CH3O] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[5];
    J[423] -= dqdci;              /* dwdot[OH]/d[CH3OH] */
    J[424] += dqdci;              /* dwdot[H2O]/d[CH3OH] */
    J[436] += dqdci;              /* dwdot[CH3O]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[467] -= dqdT;               /* dwdot[OH]/dT */
    J[468] += dqdT;               /* dwdot[H2O]/dT */
    J[480] += dqdT;               /* dwdot[CH3O]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 117: CH3OH + O2 <=> CH2OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[19];
    k_f = prefactor_units[116] * fwd_A[116]
                * exp(fwd_beta[116] * tc[0] - activation_units[116] * fwd_Ea[116] * invT);
    dlnkfdT = fwd_beta[116] * invT + activation_units[116] * fwd_Ea[116] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[17];
    Kc = exp(g_RT[4] - g_RT[7] - g_RT[17] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[19]) + (h_RT[7] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O2 */
    wdot[7] += q; /* HO2 */
    wdot[17] += q; /* CH2OH */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[19];
    J[92] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[95] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[105] += dqdci;              /* dwdot[CH2OH]/d[O2] */
    J[107] -= dqdci;              /* dwdot[CH3OH]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[17];
    J[158] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[161] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[171] += dqdci;              /* dwdot[CH2OH]/d[HO2] */
    J[173] -= dqdci;              /* dwdot[CH3OH]/d[HO2] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[7];
    J[378] -= dqdci;              /* dwdot[O2]/d[CH2OH] */
    J[381] += dqdci;              /* dwdot[HO2]/d[CH2OH] */
    J[391] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[393] -= dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[4];
    J[422] -= dqdci;              /* dwdot[O2]/d[CH3OH] */
    J[425] += dqdci;              /* dwdot[HO2]/d[CH3OH] */
    J[435] += dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[O2]/dT */
    J[469] += dqdT;               /* dwdot[HO2]/dT */
    J[479] += dqdT;               /* dwdot[CH2OH]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 118: CH3OH + CH <=> CH3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[19];
    k_f = prefactor_units[117] * fwd_A[117]
                * exp(fwd_beta[117] * tc[0] - activation_units[117] * fwd_Ea[117] * invT);
    dlnkfdT = fwd_beta[117] * invT + activation_units[117] * fwd_Ea[117] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[16];
    Kc = exp(g_RT[10] - g_RT[13] - g_RT[16] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[19]) + (h_RT[13] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH */
    wdot[13] += q; /* CH3 */
    wdot[16] += q; /* CH2O */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[19];
    J[230] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[233] += dqdci;              /* dwdot[CH3]/d[CH] */
    J[236] += dqdci;              /* dwdot[CH2O]/d[CH] */
    J[239] -= dqdci;              /* dwdot[CH3OH]/d[CH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[16];
    J[296] -= dqdci;              /* dwdot[CH]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[302] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[305] -= dqdci;              /* dwdot[CH3OH]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[13];
    J[362] -= dqdci;              /* dwdot[CH]/d[CH2O] */
    J[365] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[368] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[371] -= dqdci;              /* dwdot[CH3OH]/d[CH2O] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[10];
    J[428] -= dqdci;              /* dwdot[CH]/d[CH3OH] */
    J[431] += dqdci;              /* dwdot[CH3]/d[CH3OH] */
    J[434] += dqdci;              /* dwdot[CH2O]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[472] -= dqdT;               /* dwdot[CH]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[478] += dqdT;               /* dwdot[CH2O]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 119: CH3OH + CH2 <=> CH3 + CH2OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[19];
    k_f = prefactor_units[118] * fwd_A[118]
                * exp(fwd_beta[118] * tc[0] - activation_units[118] * fwd_Ea[118] * invT);
    dlnkfdT = fwd_beta[118] * invT + activation_units[118] * fwd_Ea[118] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[17];
    Kc = exp(g_RT[11] - g_RT[13] - g_RT[17] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[19]) + (h_RT[13] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] -= q; /* CH2 */
    wdot[13] += q; /* CH3 */
    wdot[17] += q; /* CH2OH */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[19];
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[255] += dqdci;              /* dwdot[CH3]/d[CH2] */
    J[259] += dqdci;              /* dwdot[CH2OH]/d[CH2] */
    J[261] -= dqdci;              /* dwdot[CH3OH]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[17];
    J[297] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[303] += dqdci;              /* dwdot[CH2OH]/d[CH3] */
    J[305] -= dqdci;              /* dwdot[CH3OH]/d[CH3] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[13];
    J[385] -= dqdci;              /* dwdot[CH2]/d[CH2OH] */
    J[387] += dqdci;              /* dwdot[CH3]/d[CH2OH] */
    J[391] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[393] -= dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[11];
    J[429] -= dqdci;              /* dwdot[CH2]/d[CH3OH] */
    J[431] += dqdci;              /* dwdot[CH3]/d[CH3OH] */
    J[435] += dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[479] += dqdT;               /* dwdot[CH2OH]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 120: CH3OH + CH2 <=> CH3 + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[19];
    k_f = prefactor_units[119] * fwd_A[119]
                * exp(fwd_beta[119] * tc[0] - activation_units[119] * fwd_Ea[119] * invT);
    dlnkfdT = fwd_beta[119] * invT + activation_units[119] * fwd_Ea[119] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[18];
    Kc = exp(g_RT[11] - g_RT[13] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[19]) + (h_RT[13] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] -= q; /* CH2 */
    wdot[13] += q; /* CH3 */
    wdot[18] += q; /* CH3O */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[19];
    J[253] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[255] += dqdci;              /* dwdot[CH3]/d[CH2] */
    J[260] += dqdci;              /* dwdot[CH3O]/d[CH2] */
    J[261] -= dqdci;              /* dwdot[CH3OH]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[18];
    J[297] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[304] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    J[305] -= dqdci;              /* dwdot[CH3OH]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[13];
    J[407] -= dqdci;              /* dwdot[CH2]/d[CH3O] */
    J[409] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[414] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[415] -= dqdci;              /* dwdot[CH3OH]/d[CH3O] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[11];
    J[429] -= dqdci;              /* dwdot[CH2]/d[CH3OH] */
    J[431] += dqdci;              /* dwdot[CH3]/d[CH3OH] */
    J[436] += dqdci;              /* dwdot[CH3O]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[473] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[480] += dqdT;               /* dwdot[CH3O]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 121: CH3OH + CH2(S) <=> CH3 + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[19];
    k_f = prefactor_units[120] * fwd_A[120]
                * exp(fwd_beta[120] * tc[0] - activation_units[120] * fwd_Ea[120] * invT);
    dlnkfdT = fwd_beta[120] * invT + activation_units[120] * fwd_Ea[120] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[18];
    Kc = exp(g_RT[12] - g_RT[13] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[19]) + (h_RT[13] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[12] -= q; /* CH2(S) */
    wdot[13] += q; /* CH3 */
    wdot[18] += q; /* CH3O */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[19];
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[277] += dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    J[282] += dqdci;              /* dwdot[CH3O]/d[CH2(S)] */
    J[283] -= dqdci;              /* dwdot[CH3OH]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[18];
    J[298] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[304] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    J[305] -= dqdci;              /* dwdot[CH3OH]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[13];
    J[408] -= dqdci;              /* dwdot[CH2(S)]/d[CH3O] */
    J[409] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[414] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[415] -= dqdci;              /* dwdot[CH3OH]/d[CH3O] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[12];
    J[430] -= dqdci;              /* dwdot[CH2(S)]/d[CH3OH] */
    J[431] += dqdci;              /* dwdot[CH3]/d[CH3OH] */
    J[436] += dqdci;              /* dwdot[CH3O]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[480] += dqdT;               /* dwdot[CH3O]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 122: CH3OH + CH2(S) <=> CH3 + CH2OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[19];
    k_f = prefactor_units[121] * fwd_A[121]
                * exp(fwd_beta[121] * tc[0] - activation_units[121] * fwd_Ea[121] * invT);
    dlnkfdT = fwd_beta[121] * invT + activation_units[121] * fwd_Ea[121] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[17];
    Kc = exp(g_RT[12] - g_RT[13] - g_RT[17] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[19]) + (h_RT[13] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[12] -= q; /* CH2(S) */
    wdot[13] += q; /* CH3 */
    wdot[17] += q; /* CH2OH */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[19];
    J[276] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[277] += dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    J[281] += dqdci;              /* dwdot[CH2OH]/d[CH2(S)] */
    J[283] -= dqdci;              /* dwdot[CH3OH]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[17];
    J[298] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[299] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[303] += dqdci;              /* dwdot[CH2OH]/d[CH3] */
    J[305] -= dqdci;              /* dwdot[CH3OH]/d[CH3] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[13];
    J[386] -= dqdci;              /* dwdot[CH2(S)]/d[CH2OH] */
    J[387] += dqdci;              /* dwdot[CH3]/d[CH2OH] */
    J[391] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[393] -= dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[12];
    J[430] -= dqdci;              /* dwdot[CH2(S)]/d[CH3OH] */
    J[431] += dqdci;              /* dwdot[CH3]/d[CH3OH] */
    J[435] += dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[474] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[475] += dqdT;               /* dwdot[CH3]/dT */
    J[479] += dqdT;               /* dwdot[CH2OH]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 123: CH3OH + CH3 <=> CH2OH + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[19];
    k_f = prefactor_units[122] * fwd_A[122]
                * exp(fwd_beta[122] * tc[0] - activation_units[122] * fwd_Ea[122] * invT);
    dlnkfdT = fwd_beta[122] * invT + activation_units[122] * fwd_Ea[122] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[17];
    Kc = exp(g_RT[13] - g_RT[14] - g_RT[17] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[19]) + (h_RT[14] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[13] -= q; /* CH3 */
    wdot[14] += q; /* CH4 */
    wdot[17] += q; /* CH2OH */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[19];
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[300] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[303] += dqdci;              /* dwdot[CH2OH]/d[CH3] */
    J[305] -= dqdci;              /* dwdot[CH3OH]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[17];
    J[321] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[322] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[325] += dqdci;              /* dwdot[CH2OH]/d[CH4] */
    J[327] -= dqdci;              /* dwdot[CH3OH]/d[CH4] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[14];
    J[387] -= dqdci;              /* dwdot[CH3]/d[CH2OH] */
    J[388] += dqdci;              /* dwdot[CH4]/d[CH2OH] */
    J[391] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[393] -= dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[13];
    J[431] -= dqdci;              /* dwdot[CH3]/d[CH3OH] */
    J[432] += dqdci;              /* dwdot[CH4]/d[CH3OH] */
    J[435] += dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH4]/dT */
    J[479] += dqdT;               /* dwdot[CH2OH]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    /*reaction 124: CH3OH + CH3 <=> CH3O + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[19];
    k_f = prefactor_units[123] * fwd_A[123]
                * exp(fwd_beta[123] * tc[0] - activation_units[123] * fwd_Ea[123] * invT);
    dlnkfdT = fwd_beta[123] * invT + activation_units[123] * fwd_Ea[123] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[18];
    Kc = exp(g_RT[13] - g_RT[14] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[19]) + (h_RT[14] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[13] -= q; /* CH3 */
    wdot[14] += q; /* CH4 */
    wdot[18] += q; /* CH3O */
    wdot[19] -= q; /* CH3OH */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[19];
    J[299] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[300] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[304] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    J[305] -= dqdci;              /* dwdot[CH3OH]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[18];
    J[321] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[322] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[326] += dqdci;              /* dwdot[CH3O]/d[CH4] */
    J[327] -= dqdci;              /* dwdot[CH3OH]/d[CH4] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[14];
    J[409] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[410] += dqdci;              /* dwdot[CH4]/d[CH3O] */
    J[414] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[415] -= dqdci;              /* dwdot[CH3OH]/d[CH3O] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[13];
    J[431] -= dqdci;              /* dwdot[CH3]/d[CH3OH] */
    J[432] += dqdci;              /* dwdot[CH4]/d[CH3OH] */
    J[436] += dqdci;              /* dwdot[CH3O]/d[CH3OH] */
    J[437] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[475] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH4]/dT */
    J[480] += dqdT;               /* dwdot[CH3O]/dT */
    J[481] -= dqdT;               /* dwdot[CH3OH]/dT */

    amrex::Real c_R[21], dcRdT[21], e_RT[21];
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
    for (int k = 0; k < 21; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[462+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 21; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 21; ++m) {
            dehmixdc += eh_RT[m]*J[k*22+m];
        }
        J[k*22+21] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[483] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

#endif
/* Transport function declarations  */


void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 87;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 9198;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 21;}


void egtransetNLITE(int* NLITE ) {
    *NLITE = 3;}


/*Patm in ergs/cm3 */
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 2.80134000E+01;
    WT[1] = 2.01594000E+00;
    WT[2] = 1.00797000E+00;
    WT[3] = 1.59994000E+01;
    WT[4] = 3.19988000E+01;
    WT[5] = 1.70073700E+01;
    WT[6] = 1.80153400E+01;
    WT[7] = 3.30067700E+01;
    WT[8] = 2.80105500E+01;
    WT[9] = 4.40099500E+01;
    WT[10] = 1.30191200E+01;
    WT[11] = 1.40270900E+01;
    WT[12] = 1.40270900E+01;
    WT[13] = 1.50350600E+01;
    WT[14] = 1.60430300E+01;
    WT[15] = 2.90185200E+01;
    WT[16] = 3.00264900E+01;
    WT[17] = 3.10344600E+01;
    WT[18] = 3.10344600E+01;
    WT[19] = 3.20424300E+01;
    WT[20] = 4.00260000E+00;
}


/*the lennard-jones potential well depth eps/kb in K */
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 9.75300000E+01;
    EPS[1] = 3.80000000E+01;
    EPS[2] = 1.45000000E+02;
    EPS[3] = 8.00000000E+01;
    EPS[4] = 1.07400000E+02;
    EPS[5] = 8.00000000E+01;
    EPS[6] = 5.72400000E+02;
    EPS[7] = 1.07400000E+02;
    EPS[8] = 9.81000000E+01;
    EPS[9] = 2.44000000E+02;
    EPS[10] = 8.00000000E+01;
    EPS[11] = 1.44000000E+02;
    EPS[12] = 1.44000000E+02;
    EPS[13] = 1.44000000E+02;
    EPS[14] = 1.41400000E+02;
    EPS[15] = 4.98000000E+02;
    EPS[16] = 4.98000000E+02;
    EPS[17] = 4.17000000E+02;
    EPS[18] = 4.17000000E+02;
    EPS[19] = 4.81800000E+02;
    EPS[20] = 1.02000000E+01;
}


/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 3.62100000E+00;
    SIG[1] = 2.92000000E+00;
    SIG[2] = 2.05000000E+00;
    SIG[3] = 2.75000000E+00;
    SIG[4] = 3.45800000E+00;
    SIG[5] = 2.75000000E+00;
    SIG[6] = 2.60500000E+00;
    SIG[7] = 3.45800000E+00;
    SIG[8] = 3.65000000E+00;
    SIG[9] = 3.76300000E+00;
    SIG[10] = 2.75000000E+00;
    SIG[11] = 3.80000000E+00;
    SIG[12] = 3.80000000E+00;
    SIG[13] = 3.80000000E+00;
    SIG[14] = 3.74600000E+00;
    SIG[15] = 3.59000000E+00;
    SIG[16] = 3.59000000E+00;
    SIG[17] = 3.69000000E+00;
    SIG[18] = 3.69000000E+00;
    SIG[19] = 3.62600000E+00;
    SIG[20] = 2.57600000E+00;
}


/*the dipole moment in Debye */
void egtransetDIP(amrex::Real* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[2] = 0.00000000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 0.00000000E+00;
    DIP[6] = 1.84400000E+00;
    DIP[7] = 0.00000000E+00;
    DIP[8] = 0.00000000E+00;
    DIP[9] = 0.00000000E+00;
    DIP[10] = 0.00000000E+00;
    DIP[11] = 0.00000000E+00;
    DIP[12] = 0.00000000E+00;
    DIP[13] = 0.00000000E+00;
    DIP[14] = 0.00000000E+00;
    DIP[15] = 0.00000000E+00;
    DIP[16] = 0.00000000E+00;
    DIP[17] = 1.70000000E+00;
    DIP[18] = 1.70000000E+00;
    DIP[19] = 0.00000000E+00;
    DIP[20] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 1.76000000E+00;
    POL[1] = 7.90000000E-01;
    POL[2] = 0.00000000E+00;
    POL[3] = 0.00000000E+00;
    POL[4] = 1.60000000E+00;
    POL[5] = 0.00000000E+00;
    POL[6] = 0.00000000E+00;
    POL[7] = 0.00000000E+00;
    POL[8] = 1.95000000E+00;
    POL[9] = 2.65000000E+00;
    POL[10] = 0.00000000E+00;
    POL[11] = 0.00000000E+00;
    POL[12] = 0.00000000E+00;
    POL[13] = 0.00000000E+00;
    POL[14] = 2.60000000E+00;
    POL[15] = 0.00000000E+00;
    POL[16] = 0.00000000E+00;
    POL[17] = 0.00000000E+00;
    POL[18] = 0.00000000E+00;
    POL[19] = 0.00000000E+00;
    POL[20] = 0.00000000E+00;
}


/*the rotational relaxation collision number at 298 K */
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 4.00000000E+00;
    ZROT[1] = 2.80000000E+02;
    ZROT[2] = 0.00000000E+00;
    ZROT[3] = 0.00000000E+00;
    ZROT[4] = 3.80000000E+00;
    ZROT[5] = 0.00000000E+00;
    ZROT[6] = 4.00000000E+00;
    ZROT[7] = 1.00000000E+00;
    ZROT[8] = 1.80000000E+00;
    ZROT[9] = 2.10000000E+00;
    ZROT[10] = 0.00000000E+00;
    ZROT[11] = 0.00000000E+00;
    ZROT[12] = 0.00000000E+00;
    ZROT[13] = 0.00000000E+00;
    ZROT[14] = 1.30000000E+01;
    ZROT[15] = 0.00000000E+00;
    ZROT[16] = 2.00000000E+00;
    ZROT[17] = 2.00000000E+00;
    ZROT[18] = 2.00000000E+00;
    ZROT[19] = 1.00000000E+00;
    ZROT[20] = 0.00000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 1;
    NLIN[1] = 1;
    NLIN[2] = 0;
    NLIN[3] = 0;
    NLIN[4] = 1;
    NLIN[5] = 1;
    NLIN[6] = 2;
    NLIN[7] = 2;
    NLIN[8] = 1;
    NLIN[9] = 1;
    NLIN[10] = 1;
    NLIN[11] = 1;
    NLIN[12] = 1;
    NLIN[13] = 1;
    NLIN[14] = 2;
    NLIN[15] = 2;
    NLIN[16] = 2;
    NLIN[17] = 2;
    NLIN[18] = 2;
    NLIN[19] = 2;
    NLIN[20] = 0;
}


/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -1.66504204E+01;
    COFETA[1] = 2.39101243E+00;
    COFETA[2] = -2.29774271E-01;
    COFETA[3] = 1.00508497E-02;
    COFETA[4] = -1.38569412E+01;
    COFETA[5] = 1.00155451E+00;
    COFETA[6] = -4.85884053E-02;
    COFETA[7] = 2.19590453E-03;
    COFETA[8] = -1.99976653E+01;
    COFETA[9] = 3.43887860E+00;
    COFETA[10] = -3.62118101E-01;
    COFETA[11] = 1.56207330E-02;
    COFETA[12] = -1.50421338E+01;
    COFETA[13] = 1.87267379E+00;
    COFETA[14] = -1.61417742E-01;
    COFETA[15] = 7.04848858E-03;
    COFETA[16] = -1.70973603E+01;
    COFETA[17] = 2.61741456E+00;
    COFETA[18] = -2.58645279E-01;
    COFETA[19] = 1.12774641E-02;
    COFETA[20] = -1.50115860E+01;
    COFETA[21] = 1.87267379E+00;
    COFETA[22] = -1.61417742E-01;
    COFETA[23] = 7.04848858E-03;
    COFETA[24] = -1.34815701E+01;
    COFETA[25] = -3.70881686E-02;
    COFETA[26] = 2.18633673E-01;
    COFETA[27] = -1.40233978E-02;
    COFETA[28] = -1.70818532E+01;
    COFETA[29] = 2.61741456E+00;
    COFETA[30] = -2.58645279E-01;
    COFETA[31] = 1.12774641E-02;
    COFETA[32] = -1.67033549E+01;
    COFETA[33] = 2.40514074E+00;
    COFETA[34] = -2.31616642E-01;
    COFETA[35] = 1.01308571E-02;
    COFETA[36] = -2.24069604E+01;
    COFETA[37] = 4.44991134E+00;
    COFETA[38] = -4.74276243E-01;
    COFETA[39] = 1.97165951E-02;
    COFETA[40] = -1.51451999E+01;
    COFETA[41] = 1.87267379E+00;
    COFETA[42] = -1.61417742E-01;
    COFETA[43] = 7.04848858E-03;
    COFETA[44] = -1.98575662E+01;
    COFETA[45] = 3.41753905E+00;
    COFETA[46] = -3.59449710E-01;
    COFETA[47] = 1.55095188E-02;
    COFETA[48] = -1.98575662E+01;
    COFETA[49] = 3.41753905E+00;
    COFETA[50] = -3.59449710E-01;
    COFETA[51] = 1.55095188E-02;
    COFETA[52] = -1.98228691E+01;
    COFETA[53] = 3.41753905E+00;
    COFETA[54] = -3.59449710E-01;
    COFETA[55] = 1.55095188E-02;
    COFETA[56] = -1.96274003E+01;
    COFETA[57] = 3.36875452E+00;
    COFETA[58] = -3.53428625E-01;
    COFETA[59] = 1.52618463E-02;
    COFETA[60] = -1.83953044E+01;
    COFETA[61] = 2.20428985E+00;
    COFETA[62] = -1.15019749E-01;
    COFETA[63] = 1.62760172E-03;
    COFETA[64] = -1.83782315E+01;
    COFETA[65] = 2.20428985E+00;
    COFETA[66] = -1.15019749E-01;
    COFETA[67] = 1.62760172E-03;
    COFETA[68] = -1.90417559E+01;
    COFETA[69] = 2.55787164E+00;
    COFETA[70] = -1.72818191E-01;
    COFETA[71] = 4.59391648E-03;
    COFETA[72] = -1.90417559E+01;
    COFETA[73] = 2.55787164E+00;
    COFETA[74] = -1.72818191E-01;
    COFETA[75] = 4.59391648E-03;
    COFETA[76] = -1.88670420E+01;
    COFETA[77] = 2.43712963E+00;
    COFETA[78] = -1.49140507E-01;
    COFETA[79] = 3.24451900E-03;
    COFETA[80] = -1.33872771E+01;
    COFETA[81] = 1.20385837E+00;
    COFETA[82] = -8.76638818E-02;
    COFETA[83] = 4.56865188E-03;
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = 6.92874302E+00;
    COFLAM[1] = -9.07396722E-01;
    COFLAM[2] = 2.66298040E-01;
    COFLAM[3] = -1.37535710E-02;
    COFLAM[4] = 7.26451078E+00;
    COFLAM[5] = 1.98529241E-01;
    COFLAM[6] = 4.66900697E-02;
    COFLAM[7] = -5.55926170E-04;
    COFLAM[8] = -4.47754886E-01;
    COFLAM[9] = 3.43887860E+00;
    COFLAM[10] = -3.62118101E-01;
    COFLAM[11] = 1.56207330E-02;
    COFLAM[12] = 1.74316386E+00;
    COFLAM[13] = 1.87267379E+00;
    COFLAM[14] = -1.61417742E-01;
    COFLAM[15] = 7.04848858E-03;
    COFLAM[16] = -2.34239854E-02;
    COFLAM[17] = 2.13473638E+00;
    COFLAM[18] = -1.72920629E-01;
    COFLAM[19] = 7.35398083E-03;
    COFLAM[20] = 1.04938033E+01;
    COFLAM[21] = -1.77795622E+00;
    COFLAM[22] = 3.40457181E-01;
    COFLAM[23] = -1.48973838E-02;
    COFLAM[24] = 1.59286176E+01;
    COFLAM[25] = -5.73806043E+00;
    COFLAM[26] = 1.06489376E+00;
    COFLAM[27] = -5.38121257E-02;
    COFLAM[28] = 5.54281185E+00;
    COFLAM[29] = -5.91706022E-01;
    COFLAM[30] = 2.63415906E-01;
    COFLAM[31] = -1.47653948E-02;
    COFLAM[32] = 6.23503597E+00;
    COFLAM[33] = -6.55695747E-01;
    COFLAM[34] = 2.35473686E-01;
    COFLAM[35] = -1.25171406E-02;
    COFLAM[36] = -1.09126200E+01;
    COFLAM[37] = 5.79604027E+00;
    COFLAM[38] = -5.71862549E-01;
    COFLAM[39] = 2.12287076E-02;
    COFLAM[40] = 1.66650049E+01;
    COFLAM[41] = -4.63874469E+00;
    COFLAM[42] = 7.79286990E-01;
    COFLAM[43] = -3.66242103E-02;
    COFLAM[44] = 7.21416618E+00;
    COFLAM[45] = -1.26815519E+00;
    COFLAM[46] = 3.61347900E-01;
    COFLAM[47] = -1.95206181E-02;
    COFLAM[48] = 1.05393643E+01;
    COFLAM[49] = -2.70277569E+00;
    COFLAM[50] = 5.59436717E-01;
    COFLAM[51] = -2.82541893E-02;
    COFLAM[52] = 7.71961055E+00;
    COFLAM[53] = -1.84168413E+00;
    COFLAM[54] = 4.93963837E-01;
    COFLAM[55] = -2.74966690E-02;
    COFLAM[56] = 6.54583551E+00;
    COFLAM[57] = -1.55972790E+00;
    COFLAM[58] = 4.76952301E-01;
    COFLAM[59] = -2.66704082E-02;
    COFLAM[60] = 4.30159965E+00;
    COFLAM[61] = -9.07738309E-01;
    COFLAM[62] = 3.77594510E-01;
    COFLAM[63] = -2.20998684E-02;
    COFLAM[64] = 5.07284441E+00;
    COFLAM[65] = -1.79906331E+00;
    COFLAM[66] = 5.90301854E-01;
    COFLAM[67] = -3.57739705E-02;
    COFLAM[68] = -2.36815065E+00;
    COFLAM[69] = 1.69340902E+00;
    COFLAM[70] = 7.26678084E-02;
    COFLAM[71] = -1.08598206E-02;
    COFLAM[72] = -6.40836595E+00;
    COFLAM[73] = 3.07055865E+00;
    COFLAM[74] = -7.78044853E-02;
    COFLAM[75] = -5.51842532E-03;
    COFLAM[76] = 1.19117343E+00;
    COFLAM[77] = -3.59012104E-01;
    COFLAM[78] = 4.25848680E-01;
    COFLAM[79] = -2.94281098E-02;
    COFLAM[80] = 4.78362759E+00;
    COFLAM[81] = 1.20385837E+00;
    COFLAM[82] = -8.76638818E-02;
    COFLAM[83] = 4.56865188E-03;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.50401107E+01;
    COFD[1] = 3.25712901E+00;
    COFD[2] = -2.08526776E-01;
    COFD[3] = 9.02671623E-03;
    COFD[4] = -1.19620646E+01;
    COFD[5] = 2.57673469E+00;
    COFD[6] = -1.22832347E-01;
    COFD[7] = 5.43498404E-03;
    COFD[8] = -1.43324132E+01;
    COFD[9] = 3.66582060E+00;
    COFD[10] = -2.59939797E-01;
    COFD[11] = 1.11847100E-02;
    COFD[12] = -1.40892707E+01;
    COFD[13] = 3.05464723E+00;
    COFD[14] = -1.82188775E-01;
    COFD[15] = 7.88301876E-03;
    COFD[16] = -1.52493783E+01;
    COFD[17] = 3.33715446E+00;
    COFD[18] = -2.18459311E-01;
    COFD[19] = 9.43818257E-03;
    COFD[20] = -1.41084968E+01;
    COFD[21] = 3.05464723E+00;
    COFD[22] = -1.82188775E-01;
    COFD[23] = 7.88301876E-03;
    COFD[24] = -1.97528430E+01;
    COFD[25] = 4.97033166E+00;
    COFD[26] = -4.05149611E-01;
    COFD[27] = 1.65318083E-02;
    COFD[28] = -1.52565571E+01;
    COFD[29] = 3.33715446E+00;
    COFD[30] = -2.18459311E-01;
    COFD[31] = 9.43818257E-03;
    COFD[32] = -1.50631056E+01;
    COFD[33] = 3.26275332E+00;
    COFD[34] = -2.09251071E-01;
    COFD[35] = 9.05784179E-03;
    COFD[36] = -1.78479793E+01;
    COFD[37] = 4.21958822E+00;
    COFD[38] = -3.27546464E-01;
    COFD[39] = 1.39386153E-02;
    COFD[40] = -1.40212624E+01;
    COFD[41] = 3.05464723E+00;
    COFD[42] = -1.82188775E-01;
    COFD[43] = 7.88301876E-03;
    COFD[44] = -1.59830680E+01;
    COFD[45] = 3.65919925E+00;
    COFD[46] = -2.59123935E-01;
    COFD[47] = 1.11511685E-02;
    COFD[48] = -1.59830680E+01;
    COFD[49] = 3.65919925E+00;
    COFD[50] = -2.59123935E-01;
    COFD[51] = 1.11511685E-02;
    COFD[52] = -1.60059185E+01;
    COFD[53] = 3.65919925E+00;
    COFD[54] = -2.59123935E-01;
    COFD[55] = 1.11511685E-02;
    COFD[56] = -1.59593356E+01;
    COFD[57] = 3.63991451E+00;
    COFD[58] = -2.56710765E-01;
    COFD[59] = 1.10504045E-02;
    COFD[60] = -1.95147678E+01;
    COFD[61] = 4.80230565E+00;
    COFD[62] = -3.91636105E-01;
    COFD[63] = 1.62625554E-02;
    COFD[64] = -1.95230810E+01;
    COFD[65] = 4.80230565E+00;
    COFD[66] = -3.91636105E-01;
    COFD[67] = 1.62625554E-02;
    COFD[68] = -1.93454560E+01;
    COFD[69] = 4.73801681E+00;
    COFD[70] = -3.85376377E-01;
    COFD[71] = 1.60717327E-02;
    COFD[72] = -1.93454560E+01;
    COFD[73] = 4.73801681E+00;
    COFD[74] = -3.85376377E-01;
    COFD[75] = 1.60717327E-02;
    COFD[76] = -1.94759515E+01;
    COFD[77] = 4.78102478E+00;
    COFD[78] = -3.89558560E-01;
    COFD[79] = 1.61988526E-02;
    COFD[80] = -1.06208300E+01;
    COFD[81] = 2.01222690E+00;
    COFD[82] = -4.94878760E-02;
    COFD[83] = 2.27123470E-03;
    COFD[84] = -1.19620646E+01;
    COFD[85] = 2.57673469E+00;
    COFD[86] = -1.22832347E-01;
    COFD[87] = 5.43498404E-03;
    COFD[88] = -1.00024695E+01;
    COFD[89] = 2.03860814E+00;
    COFD[90] = -5.13985906E-02;
    COFD[91] = 2.28151109E-03;
    COFD[92] = -1.13432713E+01;
    COFD[93] = 2.71750671E+00;
    COFD[94] = -1.38242753E-01;
    COFD[95] = 5.97130430E-03;
    COFD[96] = -1.13075366E+01;
    COFD[97] = 2.44841139E+00;
    COFD[98] = -1.06192871E-01;
    COFD[99] = 4.71814628E-03;
    COFD[100] = -1.20172709E+01;
    COFD[101] = 2.60884075E+00;
    COFD[102] = -1.26364058E-01;
    COFD[103] = 5.55838423E-03;
    COFD[104] = -1.13108636E+01;
    COFD[105] = 2.44841139E+00;
    COFD[106] = -1.06192871E-01;
    COFD[107] = 4.71814628E-03;
    COFD[108] = -1.67143725E+01;
    COFD[109] = 4.40335710E+00;
    COFD[110] = -3.48740776E-01;
    COFD[111] = 1.47504969E-02;
    COFD[112] = -1.20181766E+01;
    COFD[113] = 2.60884075E+00;
    COFD[114] = -1.26364058E-01;
    COFD[115] = 5.55838423E-03;
    COFD[116] = -1.19838506E+01;
    COFD[117] = 2.58175403E+00;
    COFD[118] = -1.23504762E-01;
    COFD[119] = 5.46501266E-03;
    COFD[120] = -1.38214342E+01;
    COFD[121] = 3.23223484E+00;
    COFD[122] = -2.05316141E-01;
    COFD[123] = 8.88855284E-03;
    COFD[124] = -1.12948899E+01;
    COFD[125] = 2.44841139E+00;
    COFD[126] = -1.06192871E-01;
    COFD[127] = 4.71814628E-03;
    COFD[128] = -1.24161872E+01;
    COFD[129] = 2.71289279E+00;
    COFD[130] = -1.37653284E-01;
    COFD[131] = 5.94616973E-03;
    COFD[132] = -1.24161872E+01;
    COFD[133] = 2.71289279E+00;
    COFD[134] = -1.37653284E-01;
    COFD[135] = 5.94616973E-03;
    COFD[136] = -1.24204172E+01;
    COFD[137] = 2.71289279E+00;
    COFD[138] = -1.37653284E-01;
    COFD[139] = 5.94616973E-03;
    COFD[140] = -1.23721975E+01;
    COFD[141] = 2.69959900E+00;
    COFD[142] = -1.35941255E-01;
    COFD[143] = 5.87267208E-03;
    COFD[144] = -1.58056616E+01;
    COFD[145] = 3.97677090E+00;
    COFD[146] = -2.98155122E-01;
    COFD[147] = 1.27513621E-02;
    COFD[148] = -1.58067531E+01;
    COFD[149] = 3.97677090E+00;
    COFD[150] = -2.98155122E-01;
    COFD[151] = 1.27513621E-02;
    COFD[152] = -1.56358883E+01;
    COFD[153] = 3.90212104E+00;
    COFD[154] = -2.89102817E-01;
    COFD[155] = 1.23853193E-02;
    COFD[156] = -1.56358883E+01;
    COFD[157] = 3.90212104E+00;
    COFD[158] = -2.89102817E-01;
    COFD[159] = 1.23853193E-02;
    COFD[160] = -1.57259721E+01;
    COFD[161] = 3.94363285E+00;
    COFD[162] = -2.94141037E-01;
    COFD[163] = 1.25892207E-02;
    COFD[164] = -9.63592260E+00;
    COFD[165] = 1.95552778E+00;
    COFD[166] = -4.51254149E-02;
    COFD[167] = 2.24292521E-03;
    COFD[168] = -1.43324132E+01;
    COFD[169] = 3.66582060E+00;
    COFD[170] = -2.59939797E-01;
    COFD[171] = 1.11847100E-02;
    COFD[172] = -1.13432713E+01;
    COFD[173] = 2.71750671E+00;
    COFD[174] = -1.38242753E-01;
    COFD[175] = 5.97130430E-03;
    COFD[176] = -1.45192862E+01;
    COFD[177] = 4.08130047E+00;
    COFD[178] = -3.10731457E-01;
    COFD[179] = 1.32559167E-02;
    COFD[180] = -1.33768653E+01;
    COFD[181] = 3.44083886E+00;
    COFD[182] = -2.31494284E-01;
    COFD[183] = 9.98477470E-03;
    COFD[184] = -1.45390481E+01;
    COFD[185] = 3.75925647E+00;
    COFD[186] = -2.71338319E-01;
    COFD[187] = 1.16482485E-02;
    COFD[188] = -1.33786247E+01;
    COFD[189] = 3.44083886E+00;
    COFD[190] = -2.31494284E-01;
    COFD[191] = 9.98477470E-03;
    COFD[192] = -1.80359243E+01;
    COFD[193] = 4.95606503E+00;
    COFD[194] = -3.97734369E-01;
    COFD[195] = 1.60000521E-02;
    COFD[196] = -1.45395146E+01;
    COFD[197] = 3.75925647E+00;
    COFD[198] = -2.71338319E-01;
    COFD[199] = 1.16482485E-02;
    COFD[200] = -1.43582417E+01;
    COFD[201] = 3.67141023E+00;
    COFD[202] = -2.60626493E-01;
    COFD[203] = 1.12128494E-02;
    COFD[204] = -1.70566802E+01;
    COFD[205] = 4.59542834E+00;
    COFD[206] = -3.70900044E-01;
    COFD[207] = 1.56015844E-02;
    COFD[208] = -1.33701274E+01;
    COFD[209] = 3.44083886E+00;
    COFD[210] = -2.31494284E-01;
    COFD[211] = 9.98477470E-03;
    COFD[212] = -1.55188438E+01;
    COFD[213] = 4.07293745E+00;
    COFD[214] = -3.09702472E-01;
    COFD[215] = 1.32136563E-02;
    COFD[216] = -1.55188438E+01;
    COFD[217] = 4.07293745E+00;
    COFD[218] = -3.09702472E-01;
    COFD[219] = 1.32136563E-02;
    COFD[220] = -1.55210961E+01;
    COFD[221] = 4.07293745E+00;
    COFD[222] = -3.09702472E-01;
    COFD[223] = 1.32136563E-02;
    COFD[224] = -1.54427378E+01;
    COFD[225] = 4.05065497E+00;
    COFD[226] = -3.06956719E-01;
    COFD[227] = 1.31007098E-02;
    COFD[228] = -1.83686542E+01;
    COFD[229] = 4.97400075E+00;
    COFD[230] = -4.04343341E-01;
    COFD[231] = 1.64481075E-02;
    COFD[232] = -1.83692180E+01;
    COFD[233] = 4.97400075E+00;
    COFD[234] = -4.04343341E-01;
    COFD[235] = 1.64481075E-02;
    COFD[236] = -1.81727317E+01;
    COFD[237] = 4.93176922E+00;
    COFD[238] = -4.03773085E-01;
    COFD[239] = 1.66109709E-02;
    COFD[240] = -1.81727317E+01;
    COFD[241] = 4.93176922E+00;
    COFD[242] = -4.03773085E-01;
    COFD[243] = 1.66109709E-02;
    COFD[244] = -1.83526869E+01;
    COFD[245] = 4.97155391E+00;
    COFD[246] = -4.04980248E-01;
    COFD[247] = 1.65119966E-02;
    COFD[248] = -9.42527913E+00;
    COFD[249] = 2.03733405E+00;
    COFD[250] = -5.10946137E-02;
    COFD[251] = 2.26223729E-03;
    COFD[252] = -1.40892707E+01;
    COFD[253] = 3.05464723E+00;
    COFD[254] = -1.82188775E-01;
    COFD[255] = 7.88301876E-03;
    COFD[256] = -1.13075366E+01;
    COFD[257] = 2.44841139E+00;
    COFD[258] = -1.06192871E-01;
    COFD[259] = 4.71814628E-03;
    COFD[260] = -1.33768653E+01;
    COFD[261] = 3.44083886E+00;
    COFD[262] = -2.31494284E-01;
    COFD[263] = 9.98477470E-03;
    COFD[264] = -1.30922593E+01;
    COFD[265] = 2.83211758E+00;
    COFD[266] = -1.53043887E-01;
    COFD[267] = 6.60949529E-03;
    COFD[268] = -1.43234231E+01;
    COFD[269] = 3.15427455E+00;
    COFD[270] = -1.95171921E-01;
    COFD[271] = 8.44797295E-03;
    COFD[272] = -1.31072999E+01;
    COFD[273] = 2.83211758E+00;
    COFD[274] = -1.53043887E-01;
    COFD[275] = 6.60949529E-03;
    COFD[276] = -1.85304660E+01;
    COFD[277] = 4.76438895E+00;
    COFD[278] = -3.87934679E-01;
    COFD[279] = 1.61491697E-02;
    COFD[280] = -1.43285389E+01;
    COFD[281] = 3.15427455E+00;
    COFD[282] = -1.95171921E-01;
    COFD[283] = 8.44797295E-03;
    COFD[284] = -1.41144343E+01;
    COFD[285] = 3.06075189E+00;
    COFD[286] = -1.82984149E-01;
    COFD[287] = 7.91761706E-03;
    COFD[288] = -1.68029591E+01;
    COFD[289] = 4.00529119E+00;
    COFD[290] = -3.01554873E-01;
    COFD[291] = 1.28863482E-02;
    COFD[292] = -1.30380753E+01;
    COFD[293] = 2.83211758E+00;
    COFD[294] = -1.53043887E-01;
    COFD[295] = 6.60949529E-03;
    COFD[296] = -1.50104785E+01;
    COFD[297] = 3.43330286E+00;
    COFD[298] = -2.30540152E-01;
    COFD[299] = 9.94447745E-03;
    COFD[300] = -1.50104785E+01;
    COFD[301] = 3.43330286E+00;
    COFD[302] = -2.30540152E-01;
    COFD[303] = 9.94447745E-03;
    COFD[304] = -1.50286666E+01;
    COFD[305] = 3.43330286E+00;
    COFD[306] = -2.30540152E-01;
    COFD[307] = 9.94447745E-03;
    COFD[308] = -1.49755336E+01;
    COFD[309] = 3.41371902E+00;
    COFD[310] = -2.28061415E-01;
    COFD[311] = 9.83981192E-03;
    COFD[312] = -1.86934473E+01;
    COFD[313] = 4.68745060E+00;
    COFD[314] = -3.80871893E-01;
    COFD[315] = 1.59567102E-02;
    COFD[316] = -1.86994485E+01;
    COFD[317] = 4.68745060E+00;
    COFD[318] = -3.80871893E-01;
    COFD[319] = 1.59567102E-02;
    COFD[320] = -1.82856918E+01;
    COFD[321] = 4.54152847E+00;
    COFD[322] = -3.64741565E-01;
    COFD[323] = 1.53673781E-02;
    COFD[324] = -1.82856918E+01;
    COFD[325] = 4.54152847E+00;
    COFD[326] = -3.64741565E-01;
    COFD[327] = 1.53673781E-02;
    COFD[328] = -1.86643786E+01;
    COFD[329] = 4.67185976E+00;
    COFD[330] = -3.79518917E-01;
    COFD[331] = 1.59240527E-02;
    COFD[332] = -1.02286594E+01;
    COFD[333] = 2.00678618E+00;
    COFD[334] = -4.96791457E-02;
    COFD[335] = 2.32424547E-03;
    COFD[336] = -1.52493783E+01;
    COFD[337] = 3.33715446E+00;
    COFD[338] = -2.18459311E-01;
    COFD[339] = 9.43818257E-03;
    COFD[340] = -1.20172709E+01;
    COFD[341] = 2.60884075E+00;
    COFD[342] = -1.26364058E-01;
    COFD[343] = 5.55838423E-03;
    COFD[344] = -1.45390481E+01;
    COFD[345] = 3.75925647E+00;
    COFD[346] = -2.71338319E-01;
    COFD[347] = 1.16482485E-02;
    COFD[348] = -1.43234231E+01;
    COFD[349] = 3.15427455E+00;
    COFD[350] = -1.95171921E-01;
    COFD[351] = 8.44797295E-03;
    COFD[352] = -1.55035302E+01;
    COFD[353] = 3.43469264E+00;
    COFD[354] = -2.30716101E-01;
    COFD[355] = 9.95190819E-03;
    COFD[356] = -1.43435795E+01;
    COFD[357] = 3.15427455E+00;
    COFD[358] = -1.95171921E-01;
    COFD[359] = 8.44797295E-03;
    COFD[360] = -1.98015548E+01;
    COFD[361] = 4.97304508E+00;
    COFD[362] = -4.02719974E-01;
    COFD[363] = 1.63209016E-02;
    COFD[364] = -1.55112236E+01;
    COFD[365] = 3.43469264E+00;
    COFD[366] = -2.30716101E-01;
    COFD[367] = 9.95190819E-03;
    COFD[368] = -1.52726647E+01;
    COFD[369] = 3.34263738E+00;
    COFD[370] = -2.19141928E-01;
    COFD[371] = 9.46652191E-03;
    COFD[372] = -1.81040759E+01;
    COFD[373] = 4.30855747E+00;
    COFD[374] = -3.37933955E-01;
    COFD[375] = 1.43423089E-02;
    COFD[376] = -1.42524084E+01;
    COFD[377] = 3.15427455E+00;
    COFD[378] = -1.95171921E-01;
    COFD[379] = 8.44797295E-03;
    COFD[380] = -1.62226984E+01;
    COFD[381] = 3.75275046E+00;
    COFD[382] = -2.70550264E-01;
    COFD[383] = 1.16164402E-02;
    COFD[384] = -1.62226984E+01;
    COFD[385] = 3.75275046E+00;
    COFD[386] = -2.70550264E-01;
    COFD[387] = 1.16164402E-02;
    COFD[388] = -1.62465637E+01;
    COFD[389] = 3.75275046E+00;
    COFD[390] = -2.70550264E-01;
    COFD[391] = 1.16164402E-02;
    COFD[392] = -1.62049408E+01;
    COFD[393] = 3.73560081E+00;
    COFD[394] = -2.68469477E-01;
    COFD[395] = 1.15322831E-02;
    COFD[396] = -1.97385723E+01;
    COFD[397] = 4.87482560E+00;
    COFD[398] = -3.99138155E-01;
    COFD[399] = 1.65147118E-02;
    COFD[400] = -1.97474530E+01;
    COFD[401] = 4.87482560E+00;
    COFD[402] = -3.99138155E-01;
    COFD[403] = 1.65147118E-02;
    COFD[404] = -1.95410678E+01;
    COFD[405] = 4.79863422E+00;
    COFD[406] = -3.91277602E-01;
    COFD[407] = 1.62515604E-02;
    COFD[408] = -1.95410678E+01;
    COFD[409] = 4.79863422E+00;
    COFD[410] = -3.91277602E-01;
    COFD[411] = 1.62515604E-02;
    COFD[412] = -1.96964406E+01;
    COFD[413] = 4.85139382E+00;
    COFD[414] = -3.96758178E-01;
    COFD[415] = 1.64367242E-02;
    COFD[416] = -1.05688069E+01;
    COFD[417] = 2.00249781E+00;
    COFD[418] = -4.76135431E-02;
    COFD[419] = 2.16149275E-03;
    COFD[420] = -1.41084968E+01;
    COFD[421] = 3.05464723E+00;
    COFD[422] = -1.82188775E-01;
    COFD[423] = 7.88301876E-03;
    COFD[424] = -1.13108636E+01;
    COFD[425] = 2.44841139E+00;
    COFD[426] = -1.06192871E-01;
    COFD[427] = 4.71814628E-03;
    COFD[428] = -1.33786247E+01;
    COFD[429] = 3.44083886E+00;
    COFD[430] = -2.31494284E-01;
    COFD[431] = 9.98477470E-03;
    COFD[432] = -1.31072999E+01;
    COFD[433] = 2.83211758E+00;
    COFD[434] = -1.53043887E-01;
    COFD[435] = 6.60949529E-03;
    COFD[436] = -1.43435795E+01;
    COFD[437] = 3.15427455E+00;
    COFD[438] = -1.95171921E-01;
    COFD[439] = 8.44797295E-03;
    COFD[440] = -1.31228070E+01;
    COFD[441] = 2.83211758E+00;
    COFD[442] = -1.53043887E-01;
    COFD[443] = 6.60949529E-03;
    COFD[444] = -1.85464124E+01;
    COFD[445] = 4.76438895E+00;
    COFD[446] = -3.87934679E-01;
    COFD[447] = 1.61491697E-02;
    COFD[448] = -1.43489069E+01;
    COFD[449] = 3.15427455E+00;
    COFD[450] = -1.95171921E-01;
    COFD[451] = 8.44797295E-03;
    COFD[452] = -1.41336596E+01;
    COFD[453] = 3.06075189E+00;
    COFD[454] = -1.82984149E-01;
    COFD[455] = 7.91761706E-03;
    COFD[456] = -1.68251782E+01;
    COFD[457] = 4.00529119E+00;
    COFD[458] = -3.01554873E-01;
    COFD[459] = 1.28863482E-02;
    COFD[460] = -1.30515502E+01;
    COFD[461] = 2.83211758E+00;
    COFD[462] = -1.53043887E-01;
    COFD[463] = 6.60949529E-03;
    COFD[464] = -1.50245172E+01;
    COFD[465] = 3.43330286E+00;
    COFD[466] = -2.30540152E-01;
    COFD[467] = 9.94447745E-03;
    COFD[468] = -1.50245172E+01;
    COFD[469] = 3.43330286E+00;
    COFD[470] = -2.30540152E-01;
    COFD[471] = 9.94447745E-03;
    COFD[472] = -1.50432330E+01;
    COFD[473] = 3.43330286E+00;
    COFD[474] = -2.30540152E-01;
    COFD[475] = 9.94447745E-03;
    COFD[476] = -1.49905950E+01;
    COFD[477] = 3.41371902E+00;
    COFD[478] = -2.28061415E-01;
    COFD[479] = 9.83981192E-03;
    COFD[480] = -1.87129234E+01;
    COFD[481] = 4.68745060E+00;
    COFD[482] = -3.80871893E-01;
    COFD[483] = 1.59567102E-02;
    COFD[484] = -1.87191644E+01;
    COFD[485] = 4.68745060E+00;
    COFD[486] = -3.80871893E-01;
    COFD[487] = 1.59567102E-02;
    COFD[488] = -1.83056374E+01;
    COFD[489] = 4.54152847E+00;
    COFD[490] = -3.64741565E-01;
    COFD[491] = 1.53673781E-02;
    COFD[492] = -1.83056374E+01;
    COFD[493] = 4.54152847E+00;
    COFD[494] = -3.64741565E-01;
    COFD[495] = 1.53673781E-02;
    COFD[496] = -1.86845444E+01;
    COFD[497] = 4.67185976E+00;
    COFD[498] = -3.79518917E-01;
    COFD[499] = 1.59240527E-02;
    COFD[500] = -1.02346248E+01;
    COFD[501] = 2.00678618E+00;
    COFD[502] = -4.96791457E-02;
    COFD[503] = 2.32424547E-03;
    COFD[504] = -1.97528430E+01;
    COFD[505] = 4.97033166E+00;
    COFD[506] = -4.05149611E-01;
    COFD[507] = 1.65318083E-02;
    COFD[508] = -1.67143725E+01;
    COFD[509] = 4.40335710E+00;
    COFD[510] = -3.48740776E-01;
    COFD[511] = 1.47504969E-02;
    COFD[512] = -1.80359243E+01;
    COFD[513] = 4.95606503E+00;
    COFD[514] = -3.97734369E-01;
    COFD[515] = 1.60000521E-02;
    COFD[516] = -1.85304660E+01;
    COFD[517] = 4.76438895E+00;
    COFD[518] = -3.87934679E-01;
    COFD[519] = 1.61491697E-02;
    COFD[520] = -1.98015548E+01;
    COFD[521] = 4.97304508E+00;
    COFD[522] = -4.02719974E-01;
    COFD[523] = 1.63209016E-02;
    COFD[524] = -1.85464124E+01;
    COFD[525] = 4.76438895E+00;
    COFD[526] = -3.87934679E-01;
    COFD[527] = 1.61491697E-02;
    COFD[528] = -1.37094056E+01;
    COFD[529] = 1.78861280E+00;
    COFD[530] = 1.07124793E-01;
    COFD[531] = -9.02498266E-03;
    COFD[532] = -1.95448949E+01;
    COFD[533] = 4.93655367E+00;
    COFD[534] = -4.03984062E-01;
    COFD[535] = 1.66045658E-02;
    COFD[536] = -1.97849289E+01;
    COFD[537] = 4.97246743E+00;
    COFD[538] = -4.04791570E-01;
    COFD[539] = 1.64921307E-02;
    COFD[540] = -1.92176773E+01;
    COFD[541] = 4.40341550E+00;
    COFD[542] = -2.99590977E-01;
    COFD[543] = 1.07737062E-02;
    COFD[544] = -1.84732478E+01;
    COFD[545] = 4.76438895E+00;
    COFD[546] = -3.87934679E-01;
    COFD[547] = 1.61491697E-02;
    COFD[548] = -1.97231481E+01;
    COFD[549] = 4.95541446E+00;
    COFD[550] = -3.97855837E-01;
    COFD[551] = 1.60130658E-02;
    COFD[552] = -1.97231481E+01;
    COFD[553] = 4.95541446E+00;
    COFD[554] = -3.97855837E-01;
    COFD[555] = 1.60130658E-02;
    COFD[556] = -1.97423589E+01;
    COFD[557] = 4.95541446E+00;
    COFD[558] = -3.97855837E-01;
    COFD[559] = 1.60130658E-02;
    COFD[560] = -1.97971065E+01;
    COFD[561] = 4.89976265E+00;
    COFD[562] = -3.83572886E-01;
    COFD[563] = 1.51303290E-02;
    COFD[564] = -1.68463676E+01;
    COFD[565] = 3.28774743E+00;
    COFD[566] = -1.29823485E-01;
    COFD[567] = 2.53317474E-03;
    COFD[568] = -1.68528383E+01;
    COFD[569] = 3.28774743E+00;
    COFD[570] = -1.29823485E-01;
    COFD[571] = 2.53317474E-03;
    COFD[572] = -1.65649309E+01;
    COFD[573] = 3.08026112E+00;
    COFD[574] = -9.50383509E-02;
    COFD[575] = 7.95088167E-04;
    COFD[576] = -1.65649309E+01;
    COFD[577] = 3.08026112E+00;
    COFD[578] = -9.50383509E-02;
    COFD[579] = 7.95088167E-04;
    COFD[580] = -1.70409137E+01;
    COFD[581] = 3.36747127E+00;
    COFD[582] = -1.41722400E-01;
    COFD[583] = 3.10239275E-03;
    COFD[584] = -1.23328795E+01;
    COFD[585] = 2.75990215E+00;
    COFD[586] = -1.43692930E-01;
    COFD[587] = 6.20516179E-03;
    COFD[588] = -1.52565571E+01;
    COFD[589] = 3.33715446E+00;
    COFD[590] = -2.18459311E-01;
    COFD[591] = 9.43818257E-03;
    COFD[592] = -1.20181766E+01;
    COFD[593] = 2.60884075E+00;
    COFD[594] = -1.26364058E-01;
    COFD[595] = 5.55838423E-03;
    COFD[596] = -1.45395146E+01;
    COFD[597] = 3.75925647E+00;
    COFD[598] = -2.71338319E-01;
    COFD[599] = 1.16482485E-02;
    COFD[600] = -1.43285389E+01;
    COFD[601] = 3.15427455E+00;
    COFD[602] = -1.95171921E-01;
    COFD[603] = 8.44797295E-03;
    COFD[604] = -1.55112236E+01;
    COFD[605] = 3.43469264E+00;
    COFD[606] = -2.30716101E-01;
    COFD[607] = 9.95190819E-03;
    COFD[608] = -1.43489069E+01;
    COFD[609] = 3.15427455E+00;
    COFD[610] = -1.95171921E-01;
    COFD[611] = 8.44797295E-03;
    COFD[612] = -1.95448949E+01;
    COFD[613] = 4.93655367E+00;
    COFD[614] = -4.03984062E-01;
    COFD[615] = 1.66045658E-02;
    COFD[616] = -1.55190373E+01;
    COFD[617] = 3.43469264E+00;
    COFD[618] = -2.30716101E-01;
    COFD[619] = 9.95190819E-03;
    COFD[620] = -1.52798431E+01;
    COFD[621] = 3.34263738E+00;
    COFD[622] = -2.19141928E-01;
    COFD[623] = 9.46652191E-03;
    COFD[624] = -1.81129960E+01;
    COFD[625] = 4.30855747E+00;
    COFD[626] = -3.37933955E-01;
    COFD[627] = 1.43423089E-02;
    COFD[628] = -1.42568438E+01;
    COFD[629] = 3.15427455E+00;
    COFD[630] = -1.95171921E-01;
    COFD[631] = 8.44797295E-03;
    COFD[632] = -1.62273737E+01;
    COFD[633] = 3.75275046E+00;
    COFD[634] = -2.70550264E-01;
    COFD[635] = 1.16164402E-02;
    COFD[636] = -1.62273737E+01;
    COFD[637] = 3.75275046E+00;
    COFD[638] = -2.70550264E-01;
    COFD[639] = 1.16164402E-02;
    COFD[640] = -1.62514687E+01;
    COFD[641] = 3.75275046E+00;
    COFD[642] = -2.70550264E-01;
    COFD[643] = 1.16164402E-02;
    COFD[644] = -1.62100659E+01;
    COFD[645] = 3.73560081E+00;
    COFD[646] = -2.68469477E-01;
    COFD[647] = 1.15322831E-02;
    COFD[648] = -1.97458872E+01;
    COFD[649] = 4.87482560E+00;
    COFD[650] = -3.99138155E-01;
    COFD[651] = 1.65147118E-02;
    COFD[652] = -1.97549000E+01;
    COFD[653] = 4.87482560E+00;
    COFD[654] = -3.99138155E-01;
    COFD[655] = 1.65147118E-02;
    COFD[656] = -1.93781049E+01;
    COFD[657] = 4.75043284E+00;
    COFD[658] = -3.86575520E-01;
    COFD[659] = 1.61077519E-02;
    COFD[660] = -1.93781049E+01;
    COFD[661] = 4.75043284E+00;
    COFD[662] = -3.86575520E-01;
    COFD[663] = 1.61077519E-02;
    COFD[664] = -1.97041394E+01;
    COFD[665] = 4.85139382E+00;
    COFD[666] = -3.96758178E-01;
    COFD[667] = 1.64367242E-02;
    COFD[668] = -1.05705074E+01;
    COFD[669] = 2.00249781E+00;
    COFD[670] = -4.76135431E-02;
    COFD[671] = 2.16149275E-03;
    COFD[672] = -1.50631056E+01;
    COFD[673] = 3.26275332E+00;
    COFD[674] = -2.09251071E-01;
    COFD[675] = 9.05784179E-03;
    COFD[676] = -1.19838506E+01;
    COFD[677] = 2.58175403E+00;
    COFD[678] = -1.23504762E-01;
    COFD[679] = 5.46501266E-03;
    COFD[680] = -1.43582417E+01;
    COFD[681] = 3.67141023E+00;
    COFD[682] = -2.60626493E-01;
    COFD[683] = 1.12128494E-02;
    COFD[684] = -1.41144343E+01;
    COFD[685] = 3.06075189E+00;
    COFD[686] = -1.82984149E-01;
    COFD[687] = 7.91761706E-03;
    COFD[688] = -1.52726647E+01;
    COFD[689] = 3.34263738E+00;
    COFD[690] = -2.19141928E-01;
    COFD[691] = 9.46652191E-03;
    COFD[692] = -1.41336596E+01;
    COFD[693] = 3.06075189E+00;
    COFD[694] = -1.82984149E-01;
    COFD[695] = 7.91761706E-03;
    COFD[696] = -1.97849289E+01;
    COFD[697] = 4.97246743E+00;
    COFD[698] = -4.04791570E-01;
    COFD[699] = 1.64921307E-02;
    COFD[700] = -1.52798431E+01;
    COFD[701] = 3.34263738E+00;
    COFD[702] = -2.19141928E-01;
    COFD[703] = 9.46652191E-03;
    COFD[704] = -1.50859060E+01;
    COFD[705] = 3.26830336E+00;
    COFD[706] = -2.09964261E-01;
    COFD[707] = 9.08842243E-03;
    COFD[708] = -1.78733657E+01;
    COFD[709] = 4.22576698E+00;
    COFD[710] = -3.28290672E-01;
    COFD[711] = 1.39685390E-02;
    COFD[712] = -1.40464283E+01;
    COFD[713] = 3.06075189E+00;
    COFD[714] = -1.82984149E-01;
    COFD[715] = 7.91761706E-03;
    COFD[716] = -1.60064035E+01;
    COFD[717] = 3.66477550E+00;
    COFD[718] = -2.59811053E-01;
    COFD[719] = 1.11794185E-02;
    COFD[720] = -1.60064035E+01;
    COFD[721] = 3.66477550E+00;
    COFD[722] = -2.59811053E-01;
    COFD[723] = 1.11794185E-02;
    COFD[724] = -1.60292533E+01;
    COFD[725] = 3.66477550E+00;
    COFD[726] = -2.59811053E-01;
    COFD[727] = 1.11794185E-02;
    COFD[728] = -1.59840634E+01;
    COFD[729] = 3.64607380E+00;
    COFD[730] = -2.57480919E-01;
    COFD[731] = 1.10825348E-02;
    COFD[732] = -1.95359070E+01;
    COFD[733] = 4.80619606E+00;
    COFD[734] = -3.92020389E-01;
    COFD[735] = 1.62745646E-02;
    COFD[736] = -1.95442197E+01;
    COFD[737] = 4.80619606E+00;
    COFD[738] = -3.92020389E-01;
    COFD[739] = 1.62745646E-02;
    COFD[740] = -1.93780646E+01;
    COFD[741] = 4.74505267E+00;
    COFD[742] = -3.86053165E-01;
    COFD[743] = 1.60919128E-02;
    COFD[744] = -1.93780646E+01;
    COFD[745] = 4.74505267E+00;
    COFD[746] = -3.86053165E-01;
    COFD[747] = 1.60919128E-02;
    COFD[748] = -1.94963963E+01;
    COFD[749] = 4.78465594E+00;
    COFD[750] = -3.89908546E-01;
    COFD[751] = 1.62093507E-02;
    COFD[752] = -1.06272812E+01;
    COFD[753] = 2.01052763E+00;
    COFD[754] = -4.92110545E-02;
    COFD[755] = 2.25668718E-03;
    COFD[756] = -1.78479793E+01;
    COFD[757] = 4.21958822E+00;
    COFD[758] = -3.27546464E-01;
    COFD[759] = 1.39386153E-02;
    COFD[760] = -1.38214342E+01;
    COFD[761] = 3.23223484E+00;
    COFD[762] = -2.05316141E-01;
    COFD[763] = 8.88855284E-03;
    COFD[764] = -1.70566802E+01;
    COFD[765] = 4.59542834E+00;
    COFD[766] = -3.70900044E-01;
    COFD[767] = 1.56015844E-02;
    COFD[768] = -1.68029591E+01;
    COFD[769] = 4.00529119E+00;
    COFD[770] = -3.01554873E-01;
    COFD[771] = 1.28863482E-02;
    COFD[772] = -1.81040759E+01;
    COFD[773] = 4.30855747E+00;
    COFD[774] = -3.37933955E-01;
    COFD[775] = 1.43423089E-02;
    COFD[776] = -1.68251782E+01;
    COFD[777] = 4.00529119E+00;
    COFD[778] = -3.01554873E-01;
    COFD[779] = 1.28863482E-02;
    COFD[780] = -1.92176773E+01;
    COFD[781] = 4.40341550E+00;
    COFD[782] = -2.99590977E-01;
    COFD[783] = 1.07737062E-02;
    COFD[784] = -1.81129960E+01;
    COFD[785] = 4.30855747E+00;
    COFD[786] = -3.37933955E-01;
    COFD[787] = 1.43423089E-02;
    COFD[788] = -1.78733657E+01;
    COFD[789] = 4.22576698E+00;
    COFD[790] = -3.28290672E-01;
    COFD[791] = 1.39685390E-02;
    COFD[792] = -2.02474634E+01;
    COFD[793] = 4.92580507E+00;
    COFD[794] = -4.03363183E-01;
    COFD[795] = 1.66070493E-02;
    COFD[796] = -1.67253627E+01;
    COFD[797] = 4.00529119E+00;
    COFD[798] = -3.01554873E-01;
    COFD[799] = 1.28863482E-02;
    COFD[800] = -1.87534057E+01;
    COFD[801] = 4.58907621E+00;
    COFD[802] = -3.70178644E-01;
    COFD[803] = 1.55743550E-02;
    COFD[804] = -1.87534057E+01;
    COFD[805] = 4.58907621E+00;
    COFD[806] = -3.70178644E-01;
    COFD[807] = 1.55743550E-02;
    COFD[808] = -1.87794936E+01;
    COFD[809] = 4.58907621E+00;
    COFD[810] = -3.70178644E-01;
    COFD[811] = 1.55743550E-02;
    COFD[812] = -1.87388847E+01;
    COFD[813] = 4.57236051E+00;
    COFD[814] = -3.68275904E-01;
    COFD[815] = 1.55023319E-02;
    COFD[816] = -2.03084913E+01;
    COFD[817] = 4.78985962E+00;
    COFD[818] = -3.62954786E-01;
    COFD[819] = 1.40013180E-02;
    COFD[820] = -2.03187101E+01;
    COFD[821] = 4.78985962E+00;
    COFD[822] = -3.62954786E-01;
    COFD[823] = 1.40013180E-02;
    COFD[824] = -2.04385690E+01;
    COFD[825] = 4.86152696E+00;
    COFD[826] = -3.76518225E-01;
    COFD[827] = 1.47470844E-02;
    COFD[828] = -2.04385690E+01;
    COFD[829] = 4.86152696E+00;
    COFD[830] = -3.76518225E-01;
    COFD[831] = 1.47470844E-02;
    COFD[832] = -2.03726197E+01;
    COFD[833] = 4.81168488E+00;
    COFD[834] = -3.67023353E-01;
    COFD[835] = 1.42232604E-02;
    COFD[836] = -1.15133379E+01;
    COFD[837] = 2.30296368E+00;
    COFD[838] = -8.66719714E-02;
    COFD[839] = 3.84671258E-03;
    COFD[840] = -1.40212624E+01;
    COFD[841] = 3.05464723E+00;
    COFD[842] = -1.82188775E-01;
    COFD[843] = 7.88301876E-03;
    COFD[844] = -1.12948899E+01;
    COFD[845] = 2.44841139E+00;
    COFD[846] = -1.06192871E-01;
    COFD[847] = 4.71814628E-03;
    COFD[848] = -1.33701274E+01;
    COFD[849] = 3.44083886E+00;
    COFD[850] = -2.31494284E-01;
    COFD[851] = 9.98477470E-03;
    COFD[852] = -1.30380753E+01;
    COFD[853] = 2.83211758E+00;
    COFD[854] = -1.53043887E-01;
    COFD[855] = 6.60949529E-03;
    COFD[856] = -1.42524084E+01;
    COFD[857] = 3.15427455E+00;
    COFD[858] = -1.95171921E-01;
    COFD[859] = 8.44797295E-03;
    COFD[860] = -1.30515502E+01;
    COFD[861] = 2.83211758E+00;
    COFD[862] = -1.53043887E-01;
    COFD[863] = 6.60949529E-03;
    COFD[864] = -1.84732478E+01;
    COFD[865] = 4.76438895E+00;
    COFD[866] = -3.87934679E-01;
    COFD[867] = 1.61491697E-02;
    COFD[868] = -1.42568438E+01;
    COFD[869] = 3.15427455E+00;
    COFD[870] = -1.95171921E-01;
    COFD[871] = 8.44797295E-03;
    COFD[872] = -1.40464283E+01;
    COFD[873] = 3.06075189E+00;
    COFD[874] = -1.82984149E-01;
    COFD[875] = 7.91761706E-03;
    COFD[876] = -1.67253627E+01;
    COFD[877] = 4.00529119E+00;
    COFD[878] = -3.01554873E-01;
    COFD[879] = 1.28863482E-02;
    COFD[880] = -1.29891932E+01;
    COFD[881] = 2.83211758E+00;
    COFD[882] = -1.53043887E-01;
    COFD[883] = 6.60949529E-03;
    COFD[884] = -1.49596790E+01;
    COFD[885] = 3.43330286E+00;
    COFD[886] = -2.30540152E-01;
    COFD[887] = 9.94447745E-03;
    COFD[888] = -1.49596790E+01;
    COFD[889] = 3.43330286E+00;
    COFD[890] = -2.30540152E-01;
    COFD[891] = 9.94447745E-03;
    COFD[892] = -1.49760808E+01;
    COFD[893] = 3.43330286E+00;
    COFD[894] = -2.30540152E-01;
    COFD[895] = 9.94447745E-03;
    COFD[896] = -1.49212797E+01;
    COFD[897] = 3.41371902E+00;
    COFD[898] = -2.28061415E-01;
    COFD[899] = 9.83981192E-03;
    COFD[900] = -1.86246288E+01;
    COFD[901] = 4.68745060E+00;
    COFD[902] = -3.80871893E-01;
    COFD[903] = 1.59567102E-02;
    COFD[904] = -1.86298543E+01;
    COFD[905] = 4.68745060E+00;
    COFD[906] = -3.80871893E-01;
    COFD[907] = 1.59567102E-02;
    COFD[908] = -1.82153563E+01;
    COFD[909] = 4.54152847E+00;
    COFD[910] = -3.64741565E-01;
    COFD[911] = 1.53673781E-02;
    COFD[912] = -1.82153563E+01;
    COFD[913] = 4.54152847E+00;
    COFD[914] = -3.64741565E-01;
    COFD[915] = 1.53673781E-02;
    COFD[916] = -1.85933339E+01;
    COFD[917] = 4.67185976E+00;
    COFD[918] = -3.79518917E-01;
    COFD[919] = 1.59240527E-02;
    COFD[920] = -1.02062643E+01;
    COFD[921] = 2.00678618E+00;
    COFD[922] = -4.96791457E-02;
    COFD[923] = 2.32424547E-03;
    COFD[924] = -1.59830680E+01;
    COFD[925] = 3.65919925E+00;
    COFD[926] = -2.59123935E-01;
    COFD[927] = 1.11511685E-02;
    COFD[928] = -1.24161872E+01;
    COFD[929] = 2.71289279E+00;
    COFD[930] = -1.37653284E-01;
    COFD[931] = 5.94616973E-03;
    COFD[932] = -1.55188438E+01;
    COFD[933] = 4.07293745E+00;
    COFD[934] = -3.09702472E-01;
    COFD[935] = 1.32136563E-02;
    COFD[936] = -1.50104785E+01;
    COFD[937] = 3.43330286E+00;
    COFD[938] = -2.30540152E-01;
    COFD[939] = 9.94447745E-03;
    COFD[940] = -1.62226984E+01;
    COFD[941] = 3.75275046E+00;
    COFD[942] = -2.70550264E-01;
    COFD[943] = 1.16164402E-02;
    COFD[944] = -1.50245172E+01;
    COFD[945] = 3.43330286E+00;
    COFD[946] = -2.30540152E-01;
    COFD[947] = 9.94447745E-03;
    COFD[948] = -1.97231481E+01;
    COFD[949] = 4.95541446E+00;
    COFD[950] = -3.97855837E-01;
    COFD[951] = 1.60130658E-02;
    COFD[952] = -1.62273737E+01;
    COFD[953] = 3.75275046E+00;
    COFD[954] = -2.70550264E-01;
    COFD[955] = 1.16164402E-02;
    COFD[956] = -1.60064035E+01;
    COFD[957] = 3.66477550E+00;
    COFD[958] = -2.59811053E-01;
    COFD[959] = 1.11794185E-02;
    COFD[960] = -1.87534057E+01;
    COFD[961] = 4.58907621E+00;
    COFD[962] = -3.70178644E-01;
    COFD[963] = 1.55743550E-02;
    COFD[964] = -1.49596790E+01;
    COFD[965] = 3.43330286E+00;
    COFD[966] = -2.30540152E-01;
    COFD[967] = 9.94447745E-03;
    COFD[968] = -1.70235907E+01;
    COFD[969] = 4.06453982E+00;
    COFD[970] = -3.08668976E-01;
    COFD[971] = 1.31712027E-02;
    COFD[972] = -1.70235907E+01;
    COFD[973] = 4.06453982E+00;
    COFD[974] = -3.08668976E-01;
    COFD[975] = 1.31712027E-02;
    COFD[976] = -1.70406384E+01;
    COFD[977] = 4.06453982E+00;
    COFD[978] = -3.08668976E-01;
    COFD[979] = 1.31712027E-02;
    COFD[980] = -1.69814236E+01;
    COFD[981] = 4.04284341E+00;
    COFD[982] = -3.06002129E-01;
    COFD[983] = 1.30617700E-02;
    COFD[984] = -2.00395681E+01;
    COFD[985] = 4.97362234E+00;
    COFD[986] = -4.04494602E-01;
    COFD[987] = 1.64622755E-02;
    COFD[988] = -2.00450678E+01;
    COFD[989] = 4.97362234E+00;
    COFD[990] = -4.04494602E-01;
    COFD[991] = 1.64622755E-02;
    COFD[992] = -1.98397067E+01;
    COFD[993] = 4.92900010E+00;
    COFD[994] = -4.03574955E-01;
    COFD[995] = 1.66085129E-02;
    COFD[996] = -1.98397067E+01;
    COFD[997] = 4.92900010E+00;
    COFD[998] = -4.03574955E-01;
    COFD[999] = 1.66085129E-02;
    COFD[1000] = -2.00338611E+01;
    COFD[1001] = 4.97086492E+00;
    COFD[1002] = -4.05089917E-01;
    COFD[1003] = 1.65243238E-02;
    COFD[1004] = -1.07429948E+01;
    COFD[1005] = 2.03773313E+00;
    COFD[1006] = -5.11868371E-02;
    COFD[1007] = 2.26800960E-03;
    COFD[1008] = -1.59830680E+01;
    COFD[1009] = 3.65919925E+00;
    COFD[1010] = -2.59123935E-01;
    COFD[1011] = 1.11511685E-02;
    COFD[1012] = -1.24161872E+01;
    COFD[1013] = 2.71289279E+00;
    COFD[1014] = -1.37653284E-01;
    COFD[1015] = 5.94616973E-03;
    COFD[1016] = -1.55188438E+01;
    COFD[1017] = 4.07293745E+00;
    COFD[1018] = -3.09702472E-01;
    COFD[1019] = 1.32136563E-02;
    COFD[1020] = -1.50104785E+01;
    COFD[1021] = 3.43330286E+00;
    COFD[1022] = -2.30540152E-01;
    COFD[1023] = 9.94447745E-03;
    COFD[1024] = -1.62226984E+01;
    COFD[1025] = 3.75275046E+00;
    COFD[1026] = -2.70550264E-01;
    COFD[1027] = 1.16164402E-02;
    COFD[1028] = -1.50245172E+01;
    COFD[1029] = 3.43330286E+00;
    COFD[1030] = -2.30540152E-01;
    COFD[1031] = 9.94447745E-03;
    COFD[1032] = -1.97231481E+01;
    COFD[1033] = 4.95541446E+00;
    COFD[1034] = -3.97855837E-01;
    COFD[1035] = 1.60130658E-02;
    COFD[1036] = -1.62273737E+01;
    COFD[1037] = 3.75275046E+00;
    COFD[1038] = -2.70550264E-01;
    COFD[1039] = 1.16164402E-02;
    COFD[1040] = -1.60064035E+01;
    COFD[1041] = 3.66477550E+00;
    COFD[1042] = -2.59811053E-01;
    COFD[1043] = 1.11794185E-02;
    COFD[1044] = -1.87534057E+01;
    COFD[1045] = 4.58907621E+00;
    COFD[1046] = -3.70178644E-01;
    COFD[1047] = 1.55743550E-02;
    COFD[1048] = -1.49596790E+01;
    COFD[1049] = 3.43330286E+00;
    COFD[1050] = -2.30540152E-01;
    COFD[1051] = 9.94447745E-03;
    COFD[1052] = -1.70235907E+01;
    COFD[1053] = 4.06453982E+00;
    COFD[1054] = -3.08668976E-01;
    COFD[1055] = 1.31712027E-02;
    COFD[1056] = -1.70235907E+01;
    COFD[1057] = 4.06453982E+00;
    COFD[1058] = -3.08668976E-01;
    COFD[1059] = 1.31712027E-02;
    COFD[1060] = -1.70406384E+01;
    COFD[1061] = 4.06453982E+00;
    COFD[1062] = -3.08668976E-01;
    COFD[1063] = 1.31712027E-02;
    COFD[1064] = -1.69814236E+01;
    COFD[1065] = 4.04284341E+00;
    COFD[1066] = -3.06002129E-01;
    COFD[1067] = 1.30617700E-02;
    COFD[1068] = -2.00395681E+01;
    COFD[1069] = 4.97362234E+00;
    COFD[1070] = -4.04494602E-01;
    COFD[1071] = 1.64622755E-02;
    COFD[1072] = -2.00450678E+01;
    COFD[1073] = 4.97362234E+00;
    COFD[1074] = -4.04494602E-01;
    COFD[1075] = 1.64622755E-02;
    COFD[1076] = -1.98397067E+01;
    COFD[1077] = 4.92900010E+00;
    COFD[1078] = -4.03574955E-01;
    COFD[1079] = 1.66085129E-02;
    COFD[1080] = -1.98397067E+01;
    COFD[1081] = 4.92900010E+00;
    COFD[1082] = -4.03574955E-01;
    COFD[1083] = 1.66085129E-02;
    COFD[1084] = -2.00338611E+01;
    COFD[1085] = 4.97086492E+00;
    COFD[1086] = -4.05089917E-01;
    COFD[1087] = 1.65243238E-02;
    COFD[1088] = -1.07429948E+01;
    COFD[1089] = 2.03773313E+00;
    COFD[1090] = -5.11868371E-02;
    COFD[1091] = 2.26800960E-03;
    COFD[1092] = -1.60059185E+01;
    COFD[1093] = 3.65919925E+00;
    COFD[1094] = -2.59123935E-01;
    COFD[1095] = 1.11511685E-02;
    COFD[1096] = -1.24204172E+01;
    COFD[1097] = 2.71289279E+00;
    COFD[1098] = -1.37653284E-01;
    COFD[1099] = 5.94616973E-03;
    COFD[1100] = -1.55210961E+01;
    COFD[1101] = 4.07293745E+00;
    COFD[1102] = -3.09702472E-01;
    COFD[1103] = 1.32136563E-02;
    COFD[1104] = -1.50286666E+01;
    COFD[1105] = 3.43330286E+00;
    COFD[1106] = -2.30540152E-01;
    COFD[1107] = 9.94447745E-03;
    COFD[1108] = -1.62465637E+01;
    COFD[1109] = 3.75275046E+00;
    COFD[1110] = -2.70550264E-01;
    COFD[1111] = 1.16164402E-02;
    COFD[1112] = -1.50432330E+01;
    COFD[1113] = 3.43330286E+00;
    COFD[1114] = -2.30540152E-01;
    COFD[1115] = 9.94447745E-03;
    COFD[1116] = -1.97423589E+01;
    COFD[1117] = 4.95541446E+00;
    COFD[1118] = -3.97855837E-01;
    COFD[1119] = 1.60130658E-02;
    COFD[1120] = -1.62514687E+01;
    COFD[1121] = 3.75275046E+00;
    COFD[1122] = -2.70550264E-01;
    COFD[1123] = 1.16164402E-02;
    COFD[1124] = -1.60292533E+01;
    COFD[1125] = 3.66477550E+00;
    COFD[1126] = -2.59811053E-01;
    COFD[1127] = 1.11794185E-02;
    COFD[1128] = -1.87794936E+01;
    COFD[1129] = 4.58907621E+00;
    COFD[1130] = -3.70178644E-01;
    COFD[1131] = 1.55743550E-02;
    COFD[1132] = -1.49760808E+01;
    COFD[1133] = 3.43330286E+00;
    COFD[1134] = -2.30540152E-01;
    COFD[1135] = 9.94447745E-03;
    COFD[1136] = -1.70406384E+01;
    COFD[1137] = 4.06453982E+00;
    COFD[1138] = -3.08668976E-01;
    COFD[1139] = 1.31712027E-02;
    COFD[1140] = -1.70406384E+01;
    COFD[1141] = 4.06453982E+00;
    COFD[1142] = -3.08668976E-01;
    COFD[1143] = 1.31712027E-02;
    COFD[1144] = -1.70582879E+01;
    COFD[1145] = 4.06453982E+00;
    COFD[1146] = -3.08668976E-01;
    COFD[1147] = 1.31712027E-02;
    COFD[1148] = -1.69996352E+01;
    COFD[1149] = 4.04284341E+00;
    COFD[1150] = -3.06002129E-01;
    COFD[1151] = 1.30617700E-02;
    COFD[1152] = -2.00626921E+01;
    COFD[1153] = 4.97362234E+00;
    COFD[1154] = -4.04494602E-01;
    COFD[1155] = 1.64622755E-02;
    COFD[1156] = -2.00684536E+01;
    COFD[1157] = 4.97362234E+00;
    COFD[1158] = -4.04494602E-01;
    COFD[1159] = 1.64622755E-02;
    COFD[1160] = -1.98633428E+01;
    COFD[1161] = 4.92900010E+00;
    COFD[1162] = -4.03574955E-01;
    COFD[1163] = 1.66085129E-02;
    COFD[1164] = -1.98633428E+01;
    COFD[1165] = 4.92900010E+00;
    COFD[1166] = -4.03574955E-01;
    COFD[1167] = 1.66085129E-02;
    COFD[1168] = -2.00577365E+01;
    COFD[1169] = 4.97086492E+00;
    COFD[1170] = -4.05089917E-01;
    COFD[1171] = 1.65243238E-02;
    COFD[1172] = -1.07504924E+01;
    COFD[1173] = 2.03773313E+00;
    COFD[1174] = -5.11868371E-02;
    COFD[1175] = 2.26800960E-03;
    COFD[1176] = -1.59593356E+01;
    COFD[1177] = 3.63991451E+00;
    COFD[1178] = -2.56710765E-01;
    COFD[1179] = 1.10504045E-02;
    COFD[1180] = -1.23721975E+01;
    COFD[1181] = 2.69959900E+00;
    COFD[1182] = -1.35941255E-01;
    COFD[1183] = 5.87267208E-03;
    COFD[1184] = -1.54427378E+01;
    COFD[1185] = 4.05065497E+00;
    COFD[1186] = -3.06956719E-01;
    COFD[1187] = 1.31007098E-02;
    COFD[1188] = -1.49755336E+01;
    COFD[1189] = 3.41371902E+00;
    COFD[1190] = -2.28061415E-01;
    COFD[1191] = 9.83981192E-03;
    COFD[1192] = -1.62049408E+01;
    COFD[1193] = 3.73560081E+00;
    COFD[1194] = -2.68469477E-01;
    COFD[1195] = 1.15322831E-02;
    COFD[1196] = -1.49905950E+01;
    COFD[1197] = 3.41371902E+00;
    COFD[1198] = -2.28061415E-01;
    COFD[1199] = 9.83981192E-03;
    COFD[1200] = -1.97971065E+01;
    COFD[1201] = 4.89976265E+00;
    COFD[1202] = -3.83572886E-01;
    COFD[1203] = 1.51303290E-02;
    COFD[1204] = -1.62100659E+01;
    COFD[1205] = 3.73560081E+00;
    COFD[1206] = -2.68469477E-01;
    COFD[1207] = 1.15322831E-02;
    COFD[1208] = -1.59840634E+01;
    COFD[1209] = 3.64607380E+00;
    COFD[1210] = -2.57480919E-01;
    COFD[1211] = 1.10825348E-02;
    COFD[1212] = -1.87388847E+01;
    COFD[1213] = 4.57236051E+00;
    COFD[1214] = -3.68275904E-01;
    COFD[1215] = 1.55023319E-02;
    COFD[1216] = -1.49212797E+01;
    COFD[1217] = 3.41371902E+00;
    COFD[1218] = -2.28061415E-01;
    COFD[1219] = 9.83981192E-03;
    COFD[1220] = -1.69814236E+01;
    COFD[1221] = 4.04284341E+00;
    COFD[1222] = -3.06002129E-01;
    COFD[1223] = 1.30617700E-02;
    COFD[1224] = -1.69814236E+01;
    COFD[1225] = 4.04284341E+00;
    COFD[1226] = -3.06002129E-01;
    COFD[1227] = 1.30617700E-02;
    COFD[1228] = -1.69996352E+01;
    COFD[1229] = 4.04284341E+00;
    COFD[1230] = -3.06002129E-01;
    COFD[1231] = 1.30617700E-02;
    COFD[1232] = -1.69549435E+01;
    COFD[1233] = 4.02672888E+00;
    COFD[1234] = -3.04096858E-01;
    COFD[1235] = 1.29867418E-02;
    COFD[1236] = -2.00524142E+01;
    COFD[1237] = 4.97224558E+00;
    COFD[1238] = -4.04841147E-01;
    COFD[1239] = 1.64972603E-02;
    COFD[1240] = -2.00584260E+01;
    COFD[1241] = 4.97224558E+00;
    COFD[1242] = -4.04841147E-01;
    COFD[1243] = 1.64972603E-02;
    COFD[1244] = -1.99742822E+01;
    COFD[1245] = 4.94945626E+00;
    COFD[1246] = -4.04477742E-01;
    COFD[1247] = 1.65805441E-02;
    COFD[1248] = -1.99742822E+01;
    COFD[1249] = 4.94945626E+00;
    COFD[1250] = -4.04477742E-01;
    COFD[1251] = 1.65805441E-02;
    COFD[1252] = -2.00429522E+01;
    COFD[1253] = 4.96726032E+00;
    COFD[1254] = -4.05128456E-01;
    COFD[1255] = 1.65452719E-02;
    COFD[1256] = -1.07391337E+01;
    COFD[1257] = 2.03866253E+00;
    COFD[1258] = -5.14124378E-02;
    COFD[1259] = 2.28241088E-03;
    COFD[1260] = -1.95147678E+01;
    COFD[1261] = 4.80230565E+00;
    COFD[1262] = -3.91636105E-01;
    COFD[1263] = 1.62625554E-02;
    COFD[1264] = -1.58056616E+01;
    COFD[1265] = 3.97677090E+00;
    COFD[1266] = -2.98155122E-01;
    COFD[1267] = 1.27513621E-02;
    COFD[1268] = -1.83686542E+01;
    COFD[1269] = 4.97400075E+00;
    COFD[1270] = -4.04343341E-01;
    COFD[1271] = 1.64481075E-02;
    COFD[1272] = -1.86934473E+01;
    COFD[1273] = 4.68745060E+00;
    COFD[1274] = -3.80871893E-01;
    COFD[1275] = 1.59567102E-02;
    COFD[1276] = -1.97385723E+01;
    COFD[1277] = 4.87482560E+00;
    COFD[1278] = -3.99138155E-01;
    COFD[1279] = 1.65147118E-02;
    COFD[1280] = -1.87129234E+01;
    COFD[1281] = 4.68745060E+00;
    COFD[1282] = -3.80871893E-01;
    COFD[1283] = 1.59567102E-02;
    COFD[1284] = -1.68463676E+01;
    COFD[1285] = 3.28774743E+00;
    COFD[1286] = -1.29823485E-01;
    COFD[1287] = 2.53317474E-03;
    COFD[1288] = -1.97458872E+01;
    COFD[1289] = 4.87482560E+00;
    COFD[1290] = -3.99138155E-01;
    COFD[1291] = 1.65147118E-02;
    COFD[1292] = -1.95359070E+01;
    COFD[1293] = 4.80619606E+00;
    COFD[1294] = -3.92020389E-01;
    COFD[1295] = 1.62745646E-02;
    COFD[1296] = -2.03084913E+01;
    COFD[1297] = 4.78985962E+00;
    COFD[1298] = -3.62954786E-01;
    COFD[1299] = 1.40013180E-02;
    COFD[1300] = -1.86246288E+01;
    COFD[1301] = 4.68745060E+00;
    COFD[1302] = -3.80871893E-01;
    COFD[1303] = 1.59567102E-02;
    COFD[1304] = -2.00395681E+01;
    COFD[1305] = 4.97362234E+00;
    COFD[1306] = -4.04494602E-01;
    COFD[1307] = 1.64622755E-02;
    COFD[1308] = -2.00395681E+01;
    COFD[1309] = 4.97362234E+00;
    COFD[1310] = -4.04494602E-01;
    COFD[1311] = 1.64622755E-02;
    COFD[1312] = -2.00626921E+01;
    COFD[1313] = 4.97362234E+00;
    COFD[1314] = -4.04494602E-01;
    COFD[1315] = 1.64622755E-02;
    COFD[1316] = -2.00524142E+01;
    COFD[1317] = 4.97224558E+00;
    COFD[1318] = -4.04841147E-01;
    COFD[1319] = 1.64972603E-02;
    COFD[1320] = -1.79689576E+01;
    COFD[1321] = 3.62619353E+00;
    COFD[1322] = -1.80590353E-01;
    COFD[1323] = 4.97223074E-03;
    COFD[1324] = -1.79774212E+01;
    COFD[1325] = 3.62619353E+00;
    COFD[1326] = -1.80590353E-01;
    COFD[1327] = 4.97223074E-03;
    COFD[1328] = -1.88020587E+01;
    COFD[1329] = 4.01481658E+00;
    COFD[1330] = -2.39205318E-01;
    COFD[1331] = 7.79933641E-03;
    COFD[1332] = -1.88020587E+01;
    COFD[1333] = 4.01481658E+00;
    COFD[1334] = -2.39205318E-01;
    COFD[1335] = 7.79933641E-03;
    COFD[1336] = -1.81839522E+01;
    COFD[1337] = 3.71302120E+00;
    COFD[1338] = -1.93519877E-01;
    COFD[1339] = 5.59085807E-03;
    COFD[1340] = -1.24465615E+01;
    COFD[1341] = 2.66010328E+00;
    COFD[1342] = -1.30890295E-01;
    COFD[1343] = 5.65747671E-03;
    COFD[1344] = -1.95230810E+01;
    COFD[1345] = 4.80230565E+00;
    COFD[1346] = -3.91636105E-01;
    COFD[1347] = 1.62625554E-02;
    COFD[1348] = -1.58067531E+01;
    COFD[1349] = 3.97677090E+00;
    COFD[1350] = -2.98155122E-01;
    COFD[1351] = 1.27513621E-02;
    COFD[1352] = -1.83692180E+01;
    COFD[1353] = 4.97400075E+00;
    COFD[1354] = -4.04343341E-01;
    COFD[1355] = 1.64481075E-02;
    COFD[1356] = -1.86994485E+01;
    COFD[1357] = 4.68745060E+00;
    COFD[1358] = -3.80871893E-01;
    COFD[1359] = 1.59567102E-02;
    COFD[1360] = -1.97474530E+01;
    COFD[1361] = 4.87482560E+00;
    COFD[1362] = -3.99138155E-01;
    COFD[1363] = 1.65147118E-02;
    COFD[1364] = -1.87191644E+01;
    COFD[1365] = 4.68745060E+00;
    COFD[1366] = -3.80871893E-01;
    COFD[1367] = 1.59567102E-02;
    COFD[1368] = -1.68528383E+01;
    COFD[1369] = 3.28774743E+00;
    COFD[1370] = -1.29823485E-01;
    COFD[1371] = 2.53317474E-03;
    COFD[1372] = -1.97549000E+01;
    COFD[1373] = 4.87482560E+00;
    COFD[1374] = -3.99138155E-01;
    COFD[1375] = 1.65147118E-02;
    COFD[1376] = -1.95442197E+01;
    COFD[1377] = 4.80619606E+00;
    COFD[1378] = -3.92020389E-01;
    COFD[1379] = 1.62745646E-02;
    COFD[1380] = -2.03187101E+01;
    COFD[1381] = 4.78985962E+00;
    COFD[1382] = -3.62954786E-01;
    COFD[1383] = 1.40013180E-02;
    COFD[1384] = -1.86298543E+01;
    COFD[1385] = 4.68745060E+00;
    COFD[1386] = -3.80871893E-01;
    COFD[1387] = 1.59567102E-02;
    COFD[1388] = -2.00450678E+01;
    COFD[1389] = 4.97362234E+00;
    COFD[1390] = -4.04494602E-01;
    COFD[1391] = 1.64622755E-02;
    COFD[1392] = -2.00450678E+01;
    COFD[1393] = 4.97362234E+00;
    COFD[1394] = -4.04494602E-01;
    COFD[1395] = 1.64622755E-02;
    COFD[1396] = -2.00684536E+01;
    COFD[1397] = 4.97362234E+00;
    COFD[1398] = -4.04494602E-01;
    COFD[1399] = 1.64622755E-02;
    COFD[1400] = -2.00584260E+01;
    COFD[1401] = 4.97224558E+00;
    COFD[1402] = -4.04841147E-01;
    COFD[1403] = 1.64972603E-02;
    COFD[1404] = -1.79774212E+01;
    COFD[1405] = 3.62619353E+00;
    COFD[1406] = -1.80590353E-01;
    COFD[1407] = 4.97223074E-03;
    COFD[1408] = -1.79860305E+01;
    COFD[1409] = 3.62619353E+00;
    COFD[1410] = -1.80590353E-01;
    COFD[1411] = 4.97223074E-03;
    COFD[1412] = -1.88108089E+01;
    COFD[1413] = 4.01481658E+00;
    COFD[1414] = -2.39205318E-01;
    COFD[1415] = 7.79933641E-03;
    COFD[1416] = -1.88108089E+01;
    COFD[1417] = 4.01481658E+00;
    COFD[1418] = -2.39205318E-01;
    COFD[1419] = 7.79933641E-03;
    COFD[1420] = -1.81928387E+01;
    COFD[1421] = 3.71302120E+00;
    COFD[1422] = -1.93519877E-01;
    COFD[1423] = 5.59085807E-03;
    COFD[1424] = -1.24486001E+01;
    COFD[1425] = 2.66010328E+00;
    COFD[1426] = -1.30890295E-01;
    COFD[1427] = 5.65747671E-03;
    COFD[1428] = -1.93454560E+01;
    COFD[1429] = 4.73801681E+00;
    COFD[1430] = -3.85376377E-01;
    COFD[1431] = 1.60717327E-02;
    COFD[1432] = -1.56358883E+01;
    COFD[1433] = 3.90212104E+00;
    COFD[1434] = -2.89102817E-01;
    COFD[1435] = 1.23853193E-02;
    COFD[1436] = -1.81727317E+01;
    COFD[1437] = 4.93176922E+00;
    COFD[1438] = -4.03773085E-01;
    COFD[1439] = 1.66109709E-02;
    COFD[1440] = -1.82856918E+01;
    COFD[1441] = 4.54152847E+00;
    COFD[1442] = -3.64741565E-01;
    COFD[1443] = 1.53673781E-02;
    COFD[1444] = -1.95410678E+01;
    COFD[1445] = 4.79863422E+00;
    COFD[1446] = -3.91277602E-01;
    COFD[1447] = 1.62515604E-02;
    COFD[1448] = -1.83056374E+01;
    COFD[1449] = 4.54152847E+00;
    COFD[1450] = -3.64741565E-01;
    COFD[1451] = 1.53673781E-02;
    COFD[1452] = -1.65649309E+01;
    COFD[1453] = 3.08026112E+00;
    COFD[1454] = -9.50383509E-02;
    COFD[1455] = 7.95088167E-04;
    COFD[1456] = -1.93781049E+01;
    COFD[1457] = 4.75043284E+00;
    COFD[1458] = -3.86575520E-01;
    COFD[1459] = 1.61077519E-02;
    COFD[1460] = -1.93780646E+01;
    COFD[1461] = 4.74505267E+00;
    COFD[1462] = -3.86053165E-01;
    COFD[1463] = 1.60919128E-02;
    COFD[1464] = -2.04385690E+01;
    COFD[1465] = 4.86152696E+00;
    COFD[1466] = -3.76518225E-01;
    COFD[1467] = 1.47470844E-02;
    COFD[1468] = -1.82153563E+01;
    COFD[1469] = 4.54152847E+00;
    COFD[1470] = -3.64741565E-01;
    COFD[1471] = 1.53673781E-02;
    COFD[1472] = -1.98397067E+01;
    COFD[1473] = 4.92900010E+00;
    COFD[1474] = -4.03574955E-01;
    COFD[1475] = 1.66085129E-02;
    COFD[1476] = -1.98397067E+01;
    COFD[1477] = 4.92900010E+00;
    COFD[1478] = -4.03574955E-01;
    COFD[1479] = 1.66085129E-02;
    COFD[1480] = -1.98633428E+01;
    COFD[1481] = 4.92900010E+00;
    COFD[1482] = -4.03574955E-01;
    COFD[1483] = 1.66085129E-02;
    COFD[1484] = -1.99742822E+01;
    COFD[1485] = 4.94945626E+00;
    COFD[1486] = -4.04477742E-01;
    COFD[1487] = 1.65805441E-02;
    COFD[1488] = -1.88020587E+01;
    COFD[1489] = 4.01481658E+00;
    COFD[1490] = -2.39205318E-01;
    COFD[1491] = 7.79933641E-03;
    COFD[1492] = -1.88108089E+01;
    COFD[1493] = 4.01481658E+00;
    COFD[1494] = -2.39205318E-01;
    COFD[1495] = 7.79933641E-03;
    COFD[1496] = -1.86584683E+01;
    COFD[1497] = 3.94494873E+00;
    COFD[1498] = -2.30392769E-01;
    COFD[1499] = 7.44357948E-03;
    COFD[1500] = -1.86584683E+01;
    COFD[1501] = 3.94494873E+00;
    COFD[1502] = -2.30392769E-01;
    COFD[1503] = 7.44357948E-03;
    COFD[1504] = -1.89720566E+01;
    COFD[1505] = 4.08273898E+00;
    COFD[1506] = -2.49576601E-01;
    COFD[1507] = 8.30414156E-03;
    COFD[1508] = -1.23085439E+01;
    COFD[1509] = 2.60934266E+00;
    COFD[1510] = -1.25922166E-01;
    COFD[1511] = 5.51617212E-03;
    COFD[1512] = -1.93454560E+01;
    COFD[1513] = 4.73801681E+00;
    COFD[1514] = -3.85376377E-01;
    COFD[1515] = 1.60717327E-02;
    COFD[1516] = -1.56358883E+01;
    COFD[1517] = 3.90212104E+00;
    COFD[1518] = -2.89102817E-01;
    COFD[1519] = 1.23853193E-02;
    COFD[1520] = -1.81727317E+01;
    COFD[1521] = 4.93176922E+00;
    COFD[1522] = -4.03773085E-01;
    COFD[1523] = 1.66109709E-02;
    COFD[1524] = -1.82856918E+01;
    COFD[1525] = 4.54152847E+00;
    COFD[1526] = -3.64741565E-01;
    COFD[1527] = 1.53673781E-02;
    COFD[1528] = -1.95410678E+01;
    COFD[1529] = 4.79863422E+00;
    COFD[1530] = -3.91277602E-01;
    COFD[1531] = 1.62515604E-02;
    COFD[1532] = -1.83056374E+01;
    COFD[1533] = 4.54152847E+00;
    COFD[1534] = -3.64741565E-01;
    COFD[1535] = 1.53673781E-02;
    COFD[1536] = -1.65649309E+01;
    COFD[1537] = 3.08026112E+00;
    COFD[1538] = -9.50383509E-02;
    COFD[1539] = 7.95088167E-04;
    COFD[1540] = -1.93781049E+01;
    COFD[1541] = 4.75043284E+00;
    COFD[1542] = -3.86575520E-01;
    COFD[1543] = 1.61077519E-02;
    COFD[1544] = -1.93780646E+01;
    COFD[1545] = 4.74505267E+00;
    COFD[1546] = -3.86053165E-01;
    COFD[1547] = 1.60919128E-02;
    COFD[1548] = -2.04385690E+01;
    COFD[1549] = 4.86152696E+00;
    COFD[1550] = -3.76518225E-01;
    COFD[1551] = 1.47470844E-02;
    COFD[1552] = -1.82153563E+01;
    COFD[1553] = 4.54152847E+00;
    COFD[1554] = -3.64741565E-01;
    COFD[1555] = 1.53673781E-02;
    COFD[1556] = -1.98397067E+01;
    COFD[1557] = 4.92900010E+00;
    COFD[1558] = -4.03574955E-01;
    COFD[1559] = 1.66085129E-02;
    COFD[1560] = -1.98397067E+01;
    COFD[1561] = 4.92900010E+00;
    COFD[1562] = -4.03574955E-01;
    COFD[1563] = 1.66085129E-02;
    COFD[1564] = -1.98633428E+01;
    COFD[1565] = 4.92900010E+00;
    COFD[1566] = -4.03574955E-01;
    COFD[1567] = 1.66085129E-02;
    COFD[1568] = -1.99742822E+01;
    COFD[1569] = 4.94945626E+00;
    COFD[1570] = -4.04477742E-01;
    COFD[1571] = 1.65805441E-02;
    COFD[1572] = -1.88020587E+01;
    COFD[1573] = 4.01481658E+00;
    COFD[1574] = -2.39205318E-01;
    COFD[1575] = 7.79933641E-03;
    COFD[1576] = -1.88108089E+01;
    COFD[1577] = 4.01481658E+00;
    COFD[1578] = -2.39205318E-01;
    COFD[1579] = 7.79933641E-03;
    COFD[1580] = -1.86584683E+01;
    COFD[1581] = 3.94494873E+00;
    COFD[1582] = -2.30392769E-01;
    COFD[1583] = 7.44357948E-03;
    COFD[1584] = -1.86584683E+01;
    COFD[1585] = 3.94494873E+00;
    COFD[1586] = -2.30392769E-01;
    COFD[1587] = 7.44357948E-03;
    COFD[1588] = -1.89720566E+01;
    COFD[1589] = 4.08273898E+00;
    COFD[1590] = -2.49576601E-01;
    COFD[1591] = 8.30414156E-03;
    COFD[1592] = -1.23085439E+01;
    COFD[1593] = 2.60934266E+00;
    COFD[1594] = -1.25922166E-01;
    COFD[1595] = 5.51617212E-03;
    COFD[1596] = -1.94759515E+01;
    COFD[1597] = 4.78102478E+00;
    COFD[1598] = -3.89558560E-01;
    COFD[1599] = 1.61988526E-02;
    COFD[1600] = -1.57259721E+01;
    COFD[1601] = 3.94363285E+00;
    COFD[1602] = -2.94141037E-01;
    COFD[1603] = 1.25892207E-02;
    COFD[1604] = -1.83526869E+01;
    COFD[1605] = 4.97155391E+00;
    COFD[1606] = -4.04980248E-01;
    COFD[1607] = 1.65119966E-02;
    COFD[1608] = -1.86643786E+01;
    COFD[1609] = 4.67185976E+00;
    COFD[1610] = -3.79518917E-01;
    COFD[1611] = 1.59240527E-02;
    COFD[1612] = -1.96964406E+01;
    COFD[1613] = 4.85139382E+00;
    COFD[1614] = -3.96758178E-01;
    COFD[1615] = 1.64367242E-02;
    COFD[1616] = -1.86845444E+01;
    COFD[1617] = 4.67185976E+00;
    COFD[1618] = -3.79518917E-01;
    COFD[1619] = 1.59240527E-02;
    COFD[1620] = -1.70409137E+01;
    COFD[1621] = 3.36747127E+00;
    COFD[1622] = -1.41722400E-01;
    COFD[1623] = 3.10239275E-03;
    COFD[1624] = -1.97041394E+01;
    COFD[1625] = 4.85139382E+00;
    COFD[1626] = -3.96758178E-01;
    COFD[1627] = 1.64367242E-02;
    COFD[1628] = -1.94963963E+01;
    COFD[1629] = 4.78465594E+00;
    COFD[1630] = -3.89908546E-01;
    COFD[1631] = 1.62093507E-02;
    COFD[1632] = -2.03726197E+01;
    COFD[1633] = 4.81168488E+00;
    COFD[1634] = -3.67023353E-01;
    COFD[1635] = 1.42232604E-02;
    COFD[1636] = -1.85933339E+01;
    COFD[1637] = 4.67185976E+00;
    COFD[1638] = -3.79518917E-01;
    COFD[1639] = 1.59240527E-02;
    COFD[1640] = -2.00338611E+01;
    COFD[1641] = 4.97086492E+00;
    COFD[1642] = -4.05089917E-01;
    COFD[1643] = 1.65243238E-02;
    COFD[1644] = -2.00338611E+01;
    COFD[1645] = 4.97086492E+00;
    COFD[1646] = -4.05089917E-01;
    COFD[1647] = 1.65243238E-02;
    COFD[1648] = -2.00577365E+01;
    COFD[1649] = 4.97086492E+00;
    COFD[1650] = -4.05089917E-01;
    COFD[1651] = 1.65243238E-02;
    COFD[1652] = -2.00429522E+01;
    COFD[1653] = 4.96726032E+00;
    COFD[1654] = -4.05128456E-01;
    COFD[1655] = 1.65452719E-02;
    COFD[1656] = -1.81839522E+01;
    COFD[1657] = 3.71302120E+00;
    COFD[1658] = -1.93519877E-01;
    COFD[1659] = 5.59085807E-03;
    COFD[1660] = -1.81928387E+01;
    COFD[1661] = 3.71302120E+00;
    COFD[1662] = -1.93519877E-01;
    COFD[1663] = 5.59085807E-03;
    COFD[1664] = -1.89720566E+01;
    COFD[1665] = 4.08273898E+00;
    COFD[1666] = -2.49576601E-01;
    COFD[1667] = 8.30414156E-03;
    COFD[1668] = -1.89720566E+01;
    COFD[1669] = 4.08273898E+00;
    COFD[1670] = -2.49576601E-01;
    COFD[1671] = 8.30414156E-03;
    COFD[1672] = -1.83672327E+01;
    COFD[1673] = 3.78564550E+00;
    COFD[1674] = -2.04448538E-01;
    COFD[1675] = 6.11703434E-03;
    COFD[1676] = -1.24034567E+01;
    COFD[1677] = 2.63805043E+00;
    COFD[1678] = -1.28097564E-01;
    COFD[1679] = 5.53976452E-03;
    COFD[1680] = -1.06208300E+01;
    COFD[1681] = 2.01222690E+00;
    COFD[1682] = -4.94878760E-02;
    COFD[1683] = 2.27123470E-03;
    COFD[1684] = -9.63592260E+00;
    COFD[1685] = 1.95552778E+00;
    COFD[1686] = -4.51254149E-02;
    COFD[1687] = 2.24292521E-03;
    COFD[1688] = -9.42527913E+00;
    COFD[1689] = 2.03733405E+00;
    COFD[1690] = -5.10946137E-02;
    COFD[1691] = 2.26223729E-03;
    COFD[1692] = -1.02286594E+01;
    COFD[1693] = 2.00678618E+00;
    COFD[1694] = -4.96791457E-02;
    COFD[1695] = 2.32424547E-03;
    COFD[1696] = -1.05688069E+01;
    COFD[1697] = 2.00249781E+00;
    COFD[1698] = -4.76135431E-02;
    COFD[1699] = 2.16149275E-03;
    COFD[1700] = -1.02346248E+01;
    COFD[1701] = 2.00678618E+00;
    COFD[1702] = -4.96791457E-02;
    COFD[1703] = 2.32424547E-03;
    COFD[1704] = -1.23328795E+01;
    COFD[1705] = 2.75990215E+00;
    COFD[1706] = -1.43692930E-01;
    COFD[1707] = 6.20516179E-03;
    COFD[1708] = -1.05705074E+01;
    COFD[1709] = 2.00249781E+00;
    COFD[1710] = -4.76135431E-02;
    COFD[1711] = 2.16149275E-03;
    COFD[1712] = -1.06272812E+01;
    COFD[1713] = 2.01052763E+00;
    COFD[1714] = -4.92110545E-02;
    COFD[1715] = 2.25668718E-03;
    COFD[1716] = -1.15133379E+01;
    COFD[1717] = 2.30296368E+00;
    COFD[1718] = -8.66719714E-02;
    COFD[1719] = 3.84671258E-03;
    COFD[1720] = -1.02062643E+01;
    COFD[1721] = 2.00678618E+00;
    COFD[1722] = -4.96791457E-02;
    COFD[1723] = 2.32424547E-03;
    COFD[1724] = -1.07429948E+01;
    COFD[1725] = 2.03773313E+00;
    COFD[1726] = -5.11868371E-02;
    COFD[1727] = 2.26800960E-03;
    COFD[1728] = -1.07429948E+01;
    COFD[1729] = 2.03773313E+00;
    COFD[1730] = -5.11868371E-02;
    COFD[1731] = 2.26800960E-03;
    COFD[1732] = -1.07504924E+01;
    COFD[1733] = 2.03773313E+00;
    COFD[1734] = -5.11868371E-02;
    COFD[1735] = 2.26800960E-03;
    COFD[1736] = -1.07391337E+01;
    COFD[1737] = 2.03866253E+00;
    COFD[1738] = -5.14124378E-02;
    COFD[1739] = 2.28241088E-03;
    COFD[1740] = -1.24465615E+01;
    COFD[1741] = 2.66010328E+00;
    COFD[1742] = -1.30890295E-01;
    COFD[1743] = 5.65747671E-03;
    COFD[1744] = -1.24486001E+01;
    COFD[1745] = 2.66010328E+00;
    COFD[1746] = -1.30890295E-01;
    COFD[1747] = 5.65747671E-03;
    COFD[1748] = -1.23085439E+01;
    COFD[1749] = 2.60934266E+00;
    COFD[1750] = -1.25922166E-01;
    COFD[1751] = 5.51617212E-03;
    COFD[1752] = -1.23085439E+01;
    COFD[1753] = 2.60934266E+00;
    COFD[1754] = -1.25922166E-01;
    COFD[1755] = 5.51617212E-03;
    COFD[1756] = -1.24034567E+01;
    COFD[1757] = 2.63805043E+00;
    COFD[1758] = -1.28097564E-01;
    COFD[1759] = 5.53976452E-03;
    COFD[1760] = -1.07347722E+01;
    COFD[1761] = 2.47768290E+00;
    COFD[1762] = -1.24748035E-01;
    COFD[1763] = 6.25546572E-03;
}


/*List of specs with small weight, dim NLITE */
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 1;
    KTDIF[1] = 2;
    KTDIF[2] = 20;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(amrex::Real* COFTD) {
    COFTD[0] = 4.17774559E-01;
    COFTD[1] = 6.24092547E-05;
    COFTD[2] = -1.77112335E-08;
    COFTD[3] = 1.75074941E-12;
    COFTD[4] = 0.00000000E+00;
    COFTD[5] = 0.00000000E+00;
    COFTD[6] = 0.00000000E+00;
    COFTD[7] = 0.00000000E+00;
    COFTD[8] = -1.44876961E-01;
    COFTD[9] = -4.73735161E-05;
    COFTD[10] = 1.69778255E-08;
    COFTD[11] = -1.94522371E-12;
    COFTD[12] = 3.88365021E-01;
    COFTD[13] = 3.55874668E-05;
    COFTD[14] = -6.27608890E-09;
    COFTD[15] = 2.45512494E-13;
    COFTD[16] = 4.16408613E-01;
    COFTD[17] = 7.67787145E-05;
    COFTD[18] = -2.40172383E-08;
    COFTD[19] = 2.56412236E-12;
    COFTD[20] = 3.94298297E-01;
    COFTD[21] = 3.61311571E-05;
    COFTD[22] = -6.37197233E-09;
    COFTD[23] = 2.49263330E-13;
    COFTD[24] = 8.60708754E-02;
    COFTD[25] = 4.26638760E-04;
    COFTD[26] = -1.50918120E-07;
    COFTD[27] = 1.61991441E-11;
    COFTD[28] = 4.18020194E-01;
    COFTD[29] = 7.70758628E-05;
    COFTD[30] = -2.41101896E-08;
    COFTD[31] = 2.57404600E-12;
    COFTD[32] = 4.17260755E-01;
    COFTD[33] = 6.31593768E-05;
    COFTD[34] = -1.80547664E-08;
    COFTD[35] = 1.79630904E-12;
    COFTD[36] = 3.18635144E-01;
    COFTD[37] = 2.36974242E-04;
    COFTD[38] = -8.88743088E-08;
    COFTD[39] = 1.01490459E-11;
    COFTD[40] = 3.66168497E-01;
    COFTD[41] = 3.35535090E-05;
    COFTD[42] = -5.91738678E-09;
    COFTD[43] = 2.31480531E-13;
    COFTD[44] = 3.26128667E-01;
    COFTD[45] = 1.05358772E-04;
    COFTD[46] = -3.76973458E-08;
    COFTD[47] = 4.31702138E-12;
    COFTD[48] = 3.26128667E-01;
    COFTD[49] = 1.05358772E-04;
    COFTD[50] = -3.76973458E-08;
    COFTD[51] = 4.31702138E-12;
    COFTD[52] = 3.32600239E-01;
    COFTD[53] = 1.07449471E-04;
    COFTD[54] = -3.84453973E-08;
    COFTD[55] = 4.40268668E-12;
    COFTD[56] = 3.40319986E-01;
    COFTD[57] = 1.06474046E-04;
    COFTD[58] = -3.79251714E-08;
    COFTD[59] = 4.33688014E-12;
    COFTD[60] = 1.76275701E-01;
    COFTD[61] = 3.80611815E-04;
    COFTD[62] = -1.38887198E-07;
    COFTD[63] = 1.52788486E-11;
    COFTD[64] = 1.77103675E-01;
    COFTD[65] = 3.82399565E-04;
    COFTD[66] = -1.39539557E-07;
    COFTD[67] = 1.53506140E-11;
    COFTD[68] = 1.92067131E-01;
    COFTD[69] = 3.68272498E-04;
    COFTD[70] = -1.35002986E-07;
    COFTD[71] = 1.49134243E-11;
    COFTD[72] = 1.92067131E-01;
    COFTD[73] = 3.68272498E-04;
    COFTD[74] = -1.35002986E-07;
    COFTD[75] = 1.49134243E-11;
    COFTD[76] = 1.85005265E-01;
    COFTD[77] = 3.78580173E-04;
    COFTD[78] = -1.38434695E-07;
    COFTD[79] = 1.52574969E-11;
    COFTD[80] = 1.56168192E-01;
    COFTD[81] = 7.04125423E-05;
    COFTD[82] = -3.38825072E-08;
    COFTD[83] = 1.38888977E-12;
    COFTD[84] = 2.46916034E-01;
    COFTD[85] = 3.39895237E-04;
    COFTD[86] = -1.26075780E-07;
    COFTD[87] = 1.40960289E-11;
    COFTD[88] = 1.44876961E-01;
    COFTD[89] = 4.73735161E-05;
    COFTD[90] = -1.69778255E-08;
    COFTD[91] = 1.94522371E-12;
    COFTD[92] = 0.00000000E+00;
    COFTD[93] = 0.00000000E+00;
    COFTD[94] = 0.00000000E+00;
    COFTD[95] = 0.00000000E+00;
    COFTD[96] = 2.69744936E-01;
    COFTD[97] = 2.77973299E-04;
    COFTD[98] = -1.03904812E-07;
    COFTD[99] = 1.17399211E-11;
    COFTD[100] = 2.29942818E-01;
    COFTD[101] = 3.65660170E-04;
    COFTD[102] = -1.34974583E-07;
    COFTD[103] = 1.50112544E-11;
    COFTD[104] = 2.71774454E-01;
    COFTD[105] = 2.80064726E-04;
    COFTD[106] = -1.04686575E-07;
    COFTD[107] = 1.18282504E-11;
    COFTD[108] = -8.60291246E-02;
    COFTD[109] = 5.92289115E-04;
    COFTD[110] = -1.88355038E-07;
    COFTD[111] = 1.88537586E-11;
    COFTD[112] = 2.30386064E-01;
    COFTD[113] = 3.66365029E-04;
    COFTD[114] = -1.35234765E-07;
    COFTD[115] = 1.50401906E-11;
    COFTD[116] = 2.45774219E-01;
    COFTD[117] = 3.41257025E-04;
    COFTD[118] = -1.26545982E-07;
    COFTD[119] = 1.41441072E-11;
    COFTD[120] = 6.17789746E-02;
    COFTD[121] = 5.46672074E-04;
    COFTD[122] = -1.90206485E-07;
    COFTD[123] = 2.01787924E-11;
    COFTD[124] = 2.62038079E-01;
    COFTD[125] = 2.70031350E-04;
    COFTD[126] = -1.00936157E-07;
    COFTD[127] = 1.14045009E-11;
    COFTD[128] = 1.56650704E-01;
    COFTD[129] = 3.99141229E-04;
    COFTD[130] = -1.44704494E-07;
    COFTD[131] = 1.58303065E-11;
    COFTD[132] = 1.56650704E-01;
    COFTD[133] = 3.99141229E-04;
    COFTD[134] = -1.44704494E-07;
    COFTD[135] = 1.58303065E-11;
    COFTD[136] = 1.58174720E-01;
    COFTD[137] = 4.03024375E-04;
    COFTD[138] = -1.46112288E-07;
    COFTD[139] = 1.59843156E-11;
    COFTD[140] = 1.63071412E-01;
    COFTD[141] = 4.02650019E-04;
    COFTD[142] = -1.46157144E-07;
    COFTD[143] = 1.60057919E-11;
    COFTD[144] = -6.91749304E-02;
    COFTD[145] = 6.13534385E-04;
    COFTD[146] = -1.98522630E-07;
    COFTD[147] = 2.00897552E-11;
    COFTD[148] = -6.93366297E-02;
    COFTD[149] = 6.14968548E-04;
    COFTD[150] = -1.98986685E-07;
    COFTD[151] = 2.01367158E-11;
    COFTD[152] = -4.04572060E-02;
    COFTD[153] = 6.04847437E-04;
    COFTD[154] = -1.99727759E-07;
    COFTD[155] = 2.04733310E-11;
    COFTD[156] = -4.04572060E-02;
    COFTD[157] = 6.04847437E-04;
    COFTD[158] = -1.99727759E-07;
    COFTD[159] = 2.04733310E-11;
    COFTD[160] = -6.44258988E-02;
    COFTD[161] = 6.15910473E-04;
    COFTD[162] = -2.00075276E-07;
    COFTD[163] = 2.02970982E-11;
    COFTD[164] = 3.16553217E-01;
    COFTD[165] = 6.80989933E-06;
    COFTD[166] = 7.49370144E-09;
    COFTD[167] = -2.29633870E-12;
    COFTD[168] = 3.92961041E-01;
    COFTD[169] = 2.82274387E-05;
    COFTD[170] = 1.59434652E-09;
    COFTD[171] = -3.06092665E-12;
    COFTD[172] = -1.56168192E-01;
    COFTD[173] = -7.04125423E-05;
    COFTD[174] = 3.38825072E-08;
    COFTD[175] = -1.38888977E-12;
    COFTD[176] = -3.16553217E-01;
    COFTD[177] = -6.80989933E-06;
    COFTD[178] = -7.49370144E-09;
    COFTD[179] = 2.29633870E-12;
    COFTD[180] = 3.09744733E-01;
    COFTD[181] = 3.76587299E-05;
    COFTD[182] = -6.10579809E-09;
    COFTD[183] = -2.07762098E-12;
    COFTD[184] = 4.09502860E-01;
    COFTD[185] = 2.19674771E-05;
    COFTD[186] = 4.92952585E-09;
    COFTD[187] = -3.26119041E-12;
    COFTD[188] = 3.19660650E-01;
    COFTD[189] = 3.88643060E-05;
    COFTD[190] = -6.30126415E-09;
    COFTD[191] = -2.14413225E-12;
    COFTD[192] = 2.71318472E-01;
    COFTD[193] = 9.80113451E-05;
    COFTD[194] = -3.55483703E-08;
    COFTD[195] = 4.08574079E-12;
    COFTD[196] = 4.12691944E-01;
    COFTD[197] = 2.21385532E-05;
    COFTD[198] = 4.96791550E-09;
    COFTD[199] = -3.28658757E-12;
    COFTD[200] = 3.93083748E-01;
    COFTD[201] = 2.77604915E-05;
    COFTD[202] = 1.80949924E-09;
    COFTD[203] = -3.06844983E-12;
    COFTD[204] = 4.28299006E-01;
    COFTD[205] = 2.18551562E-05;
    COFTD[206] = 1.55872276E-09;
    COFTD[207] = -1.00496866E-12;
    COFTD[208] = 2.73556754E-01;
    COFTD[209] = 3.32589995E-05;
    COFTD[210] = -5.39244781E-09;
    COFTD[211] = -1.83488915E-12;
    COFTD[212] = 2.94498859E-01;
    COFTD[213] = 6.43294281E-06;
    COFTD[214] = 6.95380548E-09;
    COFTD[215] = -2.14685415E-12;
    COFTD[216] = 2.94498859E-01;
    COFTD[217] = 6.43294281E-06;
    COFTD[218] = 6.95380548E-09;
    COFTD[219] = -2.14685415E-12;
    COFTD[220] = 3.06950528E-01;
    COFTD[221] = 6.70493325E-06;
    COFTD[222] = 7.24781844E-09;
    COFTD[223] = -2.23762501E-12;
    COFTD[224] = 3.18173191E-01;
    COFTD[225] = 7.25646834E-06;
    COFTD[226] = 7.44802399E-09;
    COFTD[227] = -2.34801317E-12;
    COFTD[228] = 3.37695457E-01;
    COFTD[229] = 9.55129821E-05;
    COFTD[230] = -3.34759386E-08;
    COFTD[231] = 3.80562349E-12;
    COFTD[232] = 3.40896402E-01;
    COFTD[233] = 9.64183297E-05;
    COFTD[234] = -3.37932500E-08;
    COFTD[235] = 3.84169620E-12;
    COFTD[236] = 3.60878733E-01;
    COFTD[237] = 7.24930124E-05;
    COFTD[238] = -2.33719013E-08;
    COFTD[239] = 2.54482330E-12;
    COFTD[240] = 3.60878733E-01;
    COFTD[241] = 7.24930124E-05;
    COFTD[242] = -2.33719013E-08;
    COFTD[243] = 2.54482330E-12;
    COFTD[244] = 3.50136058E-01;
    COFTD[245] = 9.31803648E-05;
    COFTD[246] = -3.22899670E-08;
    COFTD[247] = 3.65351327E-12;
    COFTD[248] = 0.00000000E+00;
    COFTD[249] = 0.00000000E+00;
    COFTD[250] = 0.00000000E+00;
    COFTD[251] = 0.00000000E+00;
}

#ifndef AMREX_USE_GPU
/* Initializes parameter database */
void CKINIT()
{

    // (0):  H + O2 <=> O + OH
    // (0):  H + O2 <=> O + OH
    fwd_A[18]     = 98410000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 15310;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = pow(10,-12.000000);
    is_PD[18] = 0;
    nTB[18] = 0;

    // (1):  O + H2 <=> H + OH
    // (1):  O + H2 <=> H + OH
    fwd_A[19]     = 3848000000000;
    fwd_beta[19]  = 0;
    fwd_Ea[19]    = 7950;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = pow(10,-12.000000);
    is_PD[19] = 0;
    nTB[19] = 0;

    // (2):  O + H2 <=> H + OH
    // (2):  O + H2 <=> H + OH
    fwd_A[20]     = 668700000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 19180;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-12.000000);
    is_PD[20] = 0;
    nTB[20] = 0;

    // (3):  OH + H2 <=> H + H2O
    // (3):  OH + H2 <=> H + H2O
    fwd_A[21]     = 225600000;
    fwd_beta[21]  = 1.51;
    fwd_Ea[21]    = 3437;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = pow(10,-12.000000);
    is_PD[21] = 0;
    nTB[21] = 0;

    // (4):  2.000000 OH <=> O + H2O
    // (4):  2.000000 OH <=> O + H2O
    fwd_A[22]     = 31610;
    fwd_beta[22]  = 2.4199999999999999;
    fwd_Ea[22]    = -1928;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = pow(10,-12.000000);
    is_PD[22] = 0;
    nTB[22] = 0;

    // (5):  H2 + M <=> 2.000000 H + M
    // (5):  H2 + M <=> 2.000000 H + M
    fwd_A[13]     = 4.58e+19;
    fwd_beta[13]  = -1.3999999999999999;
    fwd_Ea[13]    = 104390;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = pow(10,-6.000000);
    is_PD[13] = 0;
    nTB[13] = 9;
    TB[13] = (amrex::Real *) malloc(9 * sizeof(amrex::Real));
    TBid[13] = (int *) malloc(9 * sizeof(int));
    TBid[13][0] = 0; TB[13][0] = 1.01; // N2
    TBid[13][1] = 1; TB[13][1] = 2.5499999999999998; // H2
    TBid[13][2] = 6; TB[13][2] = 12.02; // H2O
    TBid[13][3] = 8; TB[13][3] = 1.95; // CO
    TBid[13][4] = 9; TB[13][4] = 3.8300000000000001; // CO2
    TBid[13][5] = 14; TB[13][5] = 2; // CH4
    TBid[13][6] = 16; TB[13][6] = 2.5; // CH2O
    TBid[13][7] = 19; TB[13][7] = 3; // CH3OH
    TBid[13][8] = 20; TB[13][8] = 0; // HE

    // (6):  H2 + HE <=> 2.000000 H + HE
    // (6):  H2 + HE <=> 2.000000 H + HE
    fwd_A[23]     = 5.84e+18;
    fwd_beta[23]  = -1.1000000000000001;
    fwd_Ea[23]    = 104390;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = pow(10,-12.000000);
    is_PD[23] = 0;
    nTB[23] = 0;

    // (7):  2.000000 O + M <=> O2 + M
    // (7):  2.000000 O + M <=> O2 + M
    fwd_A[14]     = 6160000000000000;
    fwd_beta[14]  = -0.5;
    fwd_Ea[14]    = 0;
    prefactor_units[14]  = 1.0000000000000002e-12;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = pow(10,-12.000000);
    is_PD[14] = 0;
    nTB[14] = 8;
    TB[14] = (amrex::Real *) malloc(8 * sizeof(amrex::Real));
    TBid[14] = (int *) malloc(8 * sizeof(int));
    TBid[14][0] = 1; TB[14][0] = 2.5; // H2
    TBid[14][1] = 6; TB[14][1] = 12; // H2O
    TBid[14][2] = 8; TB[14][2] = 1.8999999999999999; // CO
    TBid[14][3] = 9; TB[14][3] = 3.7999999999999998; // CO2
    TBid[14][4] = 14; TB[14][4] = 2; // CH4
    TBid[14][5] = 16; TB[14][5] = 2.5; // CH2O
    TBid[14][6] = 19; TB[14][6] = 3; // CH3OH
    TBid[14][7] = 20; TB[14][7] = 0; // HE

    // (8):  2.000000 O + HE <=> O2 + HE
    // (8):  2.000000 O + HE <=> O2 + HE
    fwd_A[24]     = 18900000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = -1788;
    prefactor_units[24]  = 1.0000000000000002e-12;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = pow(10,-18.000000);
    is_PD[24] = 0;
    nTB[24] = 0;

    // (9):  O + H + M <=> OH + M
    // (9):  O + H + M <=> OH + M
    fwd_A[15]     = 4.71e+18;
    fwd_beta[15]  = -1;
    fwd_Ea[15]    = 0;
    prefactor_units[15]  = 1.0000000000000002e-12;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = pow(10,-12.000000);
    is_PD[15] = 0;
    nTB[15] = 9;
    TB[15] = (amrex::Real *) malloc(9 * sizeof(amrex::Real));
    TBid[15] = (int *) malloc(9 * sizeof(int));
    TBid[15][0] = 0; TB[15][0] = 1.3200000000000001; // N2
    TBid[15][1] = 1; TB[15][1] = 2.5; // H2
    TBid[15][2] = 6; TB[15][2] = 15.800000000000001; // H2O
    TBid[15][3] = 8; TB[15][3] = 2.52; // CO
    TBid[15][4] = 9; TB[15][4] = 5.0099999999999998; // CO2
    TBid[15][5] = 14; TB[15][5] = 2; // CH4
    TBid[15][6] = 16; TB[15][6] = 2.5; // CH2O
    TBid[15][7] = 19; TB[15][7] = 3; // CH3OH
    TBid[15][8] = 20; TB[15][8] = 0.75; // HE

    // (10):  H2O + M <=> H + OH + M
    // (10):  H2O + M <=> H + OH + M
    fwd_A[16]     = 6.0599999999999999e+27;
    fwd_beta[16]  = -3.3220000000000001;
    fwd_Ea[16]    = 120800;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = pow(10,-6.000000);
    is_PD[16] = 0;
    nTB[16] = 10;
    TB[16] = (amrex::Real *) malloc(10 * sizeof(amrex::Real));
    TBid[16] = (int *) malloc(10 * sizeof(int));
    TBid[16][0] = 0; TB[16][0] = 2.46; // N2
    TBid[16][1] = 1; TB[16][1] = 3.77; // H2
    TBid[16][2] = 4; TB[16][2] = 1.5; // O2
    TBid[16][3] = 6; TB[16][3] = 0; // H2O
    TBid[16][4] = 8; TB[16][4] = 2.3999999999999999; // CO
    TBid[16][5] = 9; TB[16][5] = 4.6699999999999999; // CO2
    TBid[16][6] = 14; TB[16][6] = 2; // CH4
    TBid[16][7] = 16; TB[16][7] = 2.5; // CH2O
    TBid[16][8] = 19; TB[16][8] = 3; // CH3OH
    TBid[16][9] = 20; TB[16][9] = 1.3300000000000001; // HE

    // (11):  2.000000 H2O <=> H + OH + H2O
    // (11):  2.000000 H2O <=> H + OH + H2O
    fwd_A[25]     = 7.5280000000000001e+25;
    fwd_beta[25]  = -2.4399999999999999;
    fwd_Ea[25]    = 120200;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = pow(10,-12.000000);
    is_PD[25] = 0;
    nTB[25] = 0;

    // (12):  H + O2 (+M) <=> HO2 (+M)
    // (12):  H + O2 (+M) <=> HO2 (+M)
    fwd_A[0]     = 4565000000000;
    fwd_beta[0]  = 0.44;
    fwd_Ea[0]    = 0;
    low_A[0]     = 6.37e+20;
    low_beta[0]  = -1.72;
    low_Ea[0]    = 525;
    troe_a[0]    = 0.5;
    troe_Tsss[0] = 30;
    troe_Ts[0]   = 90000;
    troe_Tss[0]  = 90000;
    troe_len[0]  = 4;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = pow(10,-12.000000);
    is_PD[0] = 1;
    nTB[0] = 10;
    TB[0] = (amrex::Real *) malloc(10 * sizeof(amrex::Real));
    TBid[0] = (int *) malloc(10 * sizeof(int));
    TBid[0][0] = 0; TB[0][0] = 0.95999999999999996; // N2
    TBid[0][1] = 1; TB[0][1] = 1.8700000000000001; // H2
    TBid[0][2] = 4; TB[0][2] = 0.75; // O2
    TBid[0][3] = 6; TB[0][3] = 15.81; // H2O
    TBid[0][4] = 8; TB[0][4] = 1.8999999999999999; // CO
    TBid[0][5] = 9; TB[0][5] = 3.4500000000000002; // CO2
    TBid[0][6] = 14; TB[0][6] = 2; // CH4
    TBid[0][7] = 16; TB[0][7] = 2.5; // CH2O
    TBid[0][8] = 19; TB[0][8] = 3; // CH3OH
    TBid[0][9] = 20; TB[0][9] = 0.70999999999999996; // HE

    // (13):  HO2 + H <=> H2 + O2
    // (13):  HO2 + H <=> H2 + O2
    fwd_A[26]     = 2945000;
    fwd_beta[26]  = 2.0870000000000002;
    fwd_Ea[26]    = -1455;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = pow(10,-12.000000);
    is_PD[26] = 0;
    nTB[26] = 0;

    // (14):  HO2 + H <=> 2.000000 OH
    // (14):  HO2 + H <=> 2.000000 OH
    fwd_A[27]     = 58880000000000;
    fwd_beta[27]  = 0;
    fwd_Ea[27]    = 300;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = pow(10,-12.000000);
    is_PD[27] = 0;
    nTB[27] = 0;

    // (15):  HO2 + H <=> O + H2O
    // (15):  HO2 + H <=> O + H2O
    fwd_A[28]     = 1632000000000;
    fwd_beta[28]  = 0;
    fwd_Ea[28]    = 0;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = pow(10,-12.000000);
    is_PD[28] = 0;
    nTB[28] = 0;

    // (16):  HO2 + O <=> OH + O2
    // (16):  HO2 + O <=> OH + O2
    fwd_A[29]     = 16090000000000;
    fwd_beta[29]  = 0;
    fwd_Ea[29]    = -445;
    prefactor_units[29]  = 1.0000000000000002e-06;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = pow(10,-12.000000);
    is_PD[29] = 0;
    nTB[29] = 0;

    // (17):  HO2 + OH <=> H2O + O2
    // (17):  HO2 + OH <=> H2O + O2
    fwd_A[30]     = 7347000000000;
    fwd_beta[30]  = 0;
    fwd_Ea[30]    = -1093;
    prefactor_units[30]  = 1.0000000000000002e-06;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = pow(10,-12.000000);
    is_PD[30] = 0;
    nTB[30] = 0;

    // (18):  HO2 + OH <=> H2O + O2
    // (18):  HO2 + OH <=> H2O + O2
    fwd_A[31]     = 453400000000000;
    fwd_beta[31]  = 0;
    fwd_Ea[31]    = 10930;
    prefactor_units[31]  = 1.0000000000000002e-06;
    activation_units[31] = 0.50321666580471969;
    phase_units[31]      = pow(10,-12.000000);
    is_PD[31] = 0;
    nTB[31] = 0;

    // (19):  CO + O (+M) <=> CO2 (+M)
    // (19):  CO + O (+M) <=> CO2 (+M)
    fwd_A[12]     = 188000000000;
    fwd_beta[12]  = 0;
    fwd_Ea[12]    = 2430;
    low_A[12]     = 1.4e+21;
    low_beta[12]  = -2.1000000000000001;
    low_Ea[12]    = 5500;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = pow(10,-12.000000);
    is_PD[12] = 1;
    nTB[12] = 7;
    TB[12] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[12] = (int *) malloc(7 * sizeof(int));
    TBid[12][0] = 1; TB[12][0] = 2.5; // H2
    TBid[12][1] = 6; TB[12][1] = 12; // H2O
    TBid[12][2] = 8; TB[12][2] = 1.8999999999999999; // CO
    TBid[12][3] = 9; TB[12][3] = 3.7999999999999998; // CO2
    TBid[12][4] = 14; TB[12][4] = 2; // CH4
    TBid[12][5] = 16; TB[12][5] = 2.5; // CH2O
    TBid[12][6] = 19; TB[12][6] = 3; // CH3OH

    // (20):  CO + O2 <=> O + CO2
    // (20):  CO + O2 <=> O + CO2
    fwd_A[32]     = 1533000000000;
    fwd_beta[32]  = 0;
    fwd_Ea[32]    = 47700;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = pow(10,-12.000000);
    is_PD[32] = 0;
    nTB[32] = 0;

    // (21):  CO + OH <=> H + CO2
    // (21):  CO + OH <=> H + CO2
    fwd_A[33]     = 61870;
    fwd_beta[33]  = 2.0529999999999999;
    fwd_Ea[33]    = -356;
    prefactor_units[33]  = 1.0000000000000002e-06;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = pow(10,-12.000000);
    is_PD[33] = 0;
    nTB[33] = 0;

    // (22):  CO + OH <=> H + CO2
    // (22):  CO + OH <=> H + CO2
    fwd_A[34]     = 5017000000000;
    fwd_beta[34]  = -0.66400000000000003;
    fwd_Ea[34]    = 332;
    prefactor_units[34]  = 1.0000000000000002e-06;
    activation_units[34] = 0.50321666580471969;
    phase_units[34]      = pow(10,-12.000000);
    is_PD[34] = 0;
    nTB[34] = 0;

    // (23):  CO + HO2 <=> OH + CO2
    // (23):  CO + HO2 <=> OH + CO2
    fwd_A[35]     = 157000;
    fwd_beta[35]  = 2.1800000000000002;
    fwd_Ea[35]    = 17944;
    prefactor_units[35]  = 1.0000000000000002e-06;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = pow(10,-12.000000);
    is_PD[35] = 0;
    nTB[35] = 0;

    // (24):  HCO + M <=> H + CO + M
    // (24):  HCO + M <=> H + CO + M
    fwd_A[17]     = 4.8e+17;
    fwd_beta[17]  = -1.2;
    fwd_Ea[17]    = 17734;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = pow(10,-6.000000);
    is_PD[17] = 0;
    nTB[17] = 10;
    TB[17] = (amrex::Real *) malloc(10 * sizeof(amrex::Real));
    TBid[17] = (int *) malloc(10 * sizeof(int));
    TBid[17][0] = 0; TB[17][0] = 1.3100000000000001; // N2
    TBid[17][1] = 1; TB[17][1] = 1.3100000000000001; // H2
    TBid[17][2] = 4; TB[17][2] = 1.3200000000000001; // O2
    TBid[17][3] = 6; TB[17][3] = 15.31; // H2O
    TBid[17][4] = 8; TB[17][4] = 2.04; // CO
    TBid[17][5] = 9; TB[17][5] = 2; // CO2
    TBid[17][6] = 14; TB[17][6] = 2.6000000000000001; // CH4
    TBid[17][7] = 16; TB[17][7] = 3.29; // CH2O
    TBid[17][8] = 19; TB[17][8] = 3; // CH3OH
    TBid[17][9] = 20; TB[17][9] = 1.3100000000000001; // HE

    // (25):  HCO + H <=> H2 + CO
    // (25):  HCO + H <=> H2 + CO
    fwd_A[36]     = 84820000000000;
    fwd_beta[36]  = 0;
    fwd_Ea[36]    = 0;
    prefactor_units[36]  = 1.0000000000000002e-06;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = pow(10,-12.000000);
    is_PD[36] = 0;
    nTB[36] = 0;

    // (26):  HCO + O <=> OH + CO
    // (26):  HCO + O <=> OH + CO
    fwd_A[37]     = 30100000000000;
    fwd_beta[37]  = 0;
    fwd_Ea[37]    = 0;
    prefactor_units[37]  = 1.0000000000000002e-06;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = pow(10,-12.000000);
    is_PD[37] = 0;
    nTB[37] = 0;

    // (27):  HCO + O <=> H + CO2
    // (27):  HCO + O <=> H + CO2
    fwd_A[38]     = 30010000000000;
    fwd_beta[38]  = 0;
    fwd_Ea[38]    = 0;
    prefactor_units[38]  = 1.0000000000000002e-06;
    activation_units[38] = 0.50321666580471969;
    phase_units[38]      = pow(10,-12.000000);
    is_PD[38] = 0;
    nTB[38] = 0;

    // (28):  HCO + OH <=> H2O + CO
    // (28):  HCO + OH <=> H2O + CO
    fwd_A[39]     = 119900000000000;
    fwd_beta[39]  = 0;
    fwd_Ea[39]    = 0;
    prefactor_units[39]  = 1.0000000000000002e-06;
    activation_units[39] = 0.50321666580471969;
    phase_units[39]      = pow(10,-12.000000);
    is_PD[39] = 0;
    nTB[39] = 0;

    // (29):  HCO + O2 <=> HO2 + CO
    // (29):  HCO + O2 <=> HO2 + CO
    fwd_A[40]     = 75620000000;
    fwd_beta[40]  = 0.52100000000000002;
    fwd_Ea[40]    = -521;
    prefactor_units[40]  = 1.0000000000000002e-06;
    activation_units[40] = 0.50321666580471969;
    phase_units[40]      = pow(10,-12.000000);
    is_PD[40] = 0;
    nTB[40] = 0;

    // (30):  CH + O <=> H + CO
    // (30):  CH + O <=> H + CO
    fwd_A[41]     = 57000000000000;
    fwd_beta[41]  = 0;
    fwd_Ea[41]    = 0;
    prefactor_units[41]  = 1.0000000000000002e-06;
    activation_units[41] = 0.50321666580471969;
    phase_units[41]      = pow(10,-12.000000);
    is_PD[41] = 0;
    nTB[41] = 0;

    // (31):  CH + OH <=> H + HCO
    // (31):  CH + OH <=> H + HCO
    fwd_A[42]     = 30000000000000;
    fwd_beta[42]  = 0;
    fwd_Ea[42]    = 0;
    prefactor_units[42]  = 1.0000000000000002e-06;
    activation_units[42] = 0.50321666580471969;
    phase_units[42]      = pow(10,-12.000000);
    is_PD[42] = 0;
    nTB[42] = 0;

    // (32):  CH + H2 <=> H + CH2
    // (32):  CH + H2 <=> H + CH2
    fwd_A[43]     = 161200000000000;
    fwd_beta[43]  = 0;
    fwd_Ea[43]    = 3320;
    prefactor_units[43]  = 1.0000000000000002e-06;
    activation_units[43] = 0.50321666580471969;
    phase_units[43]      = pow(10,-12.000000);
    is_PD[43] = 0;
    nTB[43] = 0;

    // (33):  CH + H2 (+M) <=> CH3 (+M)
    // (33):  CH + H2 (+M) <=> CH3 (+M)
    fwd_A[1]     = 51300000000000;
    fwd_beta[1]  = 0.14999999999999999;
    fwd_Ea[1]    = 0;
    low_A[1]     = 2.4300000000000001e+22;
    low_beta[1]  = -1.6000000000000001;
    low_Ea[1]    = 0;
    troe_a[1]    = 0.51400000000000001;
    troe_Tsss[1] = 152;
    troe_Ts[1]   = 22850;
    troe_Tss[1]  = 10350;
    troe_len[1]  = 4;
    prefactor_units[1]  = 1.0000000000000002e-06;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = pow(10,-12.000000);
    is_PD[1] = 1;
    nTB[1] = 7;
    TB[1] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[1] = (int *) malloc(7 * sizeof(int));
    TBid[1][0] = 6; TB[1][0] = 6; // H2O
    TBid[1][1] = 8; TB[1][1] = 1.5; // CO
    TBid[1][2] = 9; TB[1][2] = 2; // CO2
    TBid[1][3] = 14; TB[1][3] = 2; // CH4
    TBid[1][4] = 16; TB[1][4] = 2.5; // CH2O
    TBid[1][5] = 19; TB[1][5] = 3; // CH3OH
    TBid[1][6] = 20; TB[1][6] = 0.69999999999999996; // HE

    // (34):  CH + H2O <=> H + CH2O
    // (34):  CH + H2O <=> H + CH2O
    fwd_A[44]     = 3430000000000;
    fwd_beta[44]  = 0;
    fwd_Ea[44]    = -884;
    prefactor_units[44]  = 1.0000000000000002e-06;
    activation_units[44] = 0.50321666580471969;
    phase_units[44]      = pow(10,-12.000000);
    is_PD[44] = 0;
    nTB[44] = 0;

    // (35):  CH + O2 <=> O + HCO
    // (35):  CH + O2 <=> O + HCO
    fwd_A[45]     = 184000000;
    fwd_beta[45]  = 1.4299999999999999;
    fwd_Ea[45]    = 1200;
    prefactor_units[45]  = 1.0000000000000002e-06;
    activation_units[45] = 0.50321666580471969;
    phase_units[45]      = pow(10,-12.000000);
    is_PD[45] = 0;
    nTB[45] = 0;

    // (36):  CH + O2 <=> CO2 + H
    // (36):  CH + O2 <=> CO2 + H
    fwd_A[46]     = 278100000;
    fwd_beta[46]  = 1.4299999999999999;
    fwd_Ea[46]    = 1200;
    prefactor_units[46]  = 1.0000000000000002e-06;
    activation_units[46] = 0.50321666580471969;
    phase_units[46]      = pow(10,-12.000000);
    is_PD[46] = 0;
    nTB[46] = 0;

    // (37):  CH + O2 <=> CO + OH
    // (37):  CH + O2 <=> CO + OH
    fwd_A[47]     = 184000000;
    fwd_beta[47]  = 1.4299999999999999;
    fwd_Ea[47]    = 1200;
    prefactor_units[47]  = 1.0000000000000002e-06;
    activation_units[47] = 0.50321666580471969;
    phase_units[47]      = pow(10,-12.000000);
    is_PD[47] = 0;
    nTB[47] = 0;

    // (38):  CH + O2 => O + H + CO
    // (38):  CH + O2 => O + H + CO
    fwd_A[48]     = 278900000;
    fwd_beta[48]  = 1.4299999999999999;
    fwd_Ea[48]    = 1200;
    prefactor_units[48]  = 1.0000000000000002e-06;
    activation_units[48] = 0.50321666580471969;
    phase_units[48]      = pow(10,-12.000000);
    is_PD[48] = 0;
    nTB[48] = 0;

    // (39):  CH + CO2 <=> HCO + CO
    // (39):  CH + CO2 <=> HCO + CO
    fwd_A[49]     = 63800000;
    fwd_beta[49]  = 1.51;
    fwd_Ea[49]    = -715;
    prefactor_units[49]  = 1.0000000000000002e-06;
    activation_units[49] = 0.50321666580471969;
    phase_units[49]      = pow(10,-12.000000);
    is_PD[49] = 0;
    nTB[49] = 0;

    // (40):  CH2 + H (+M) <=> CH3 (+M)
    // (40):  CH2 + H (+M) <=> CH3 (+M)
    fwd_A[2]     = 21300000000000;
    fwd_beta[2]  = 0.32000000000000001;
    fwd_Ea[2]    = 0;
    low_A[2]     = 1.3900000000000001e+34;
    low_beta[2]  = -5.04;
    low_Ea[2]    = 7400;
    troe_a[2]    = 0.40500000000000003;
    troe_Tsss[2] = 258;
    troe_Ts[2]   = 2811;
    troe_Tss[2]  = 9908;
    troe_len[2]  = 4;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-12.000000);
    is_PD[2] = 1;
    nTB[2] = 7;
    TB[2] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[2] = (int *) malloc(7 * sizeof(int));
    TBid[2][0] = 6; TB[2][0] = 6; // H2O
    TBid[2][1] = 8; TB[2][1] = 1.5; // CO
    TBid[2][2] = 9; TB[2][2] = 2; // CO2
    TBid[2][3] = 14; TB[2][3] = 2; // CH4
    TBid[2][4] = 16; TB[2][4] = 2.5; // CH2O
    TBid[2][5] = 19; TB[2][5] = 3; // CH3OH
    TBid[2][6] = 20; TB[2][6] = 0.69999999999999996; // HE

    // (41):  CH2 + O => 2.000000 H + CO
    // (41):  CH2 + O => 2.000000 H + CO
    fwd_A[50]     = 80000000000000;
    fwd_beta[50]  = 0;
    fwd_Ea[50]    = 0;
    prefactor_units[50]  = 1.0000000000000002e-06;
    activation_units[50] = 0.50321666580471969;
    phase_units[50]      = pow(10,-12.000000);
    is_PD[50] = 0;
    nTB[50] = 0;

    // (42):  CH2 + OH <=> H + CH2O
    // (42):  CH2 + OH <=> H + CH2O
    fwd_A[51]     = 28990000000000;
    fwd_beta[51]  = 0.12;
    fwd_Ea[51]    = -162;
    prefactor_units[51]  = 1.0000000000000002e-06;
    activation_units[51] = 0.50321666580471969;
    phase_units[51]      = pow(10,-12.000000);
    is_PD[51] = 0;
    nTB[51] = 0;

    // (43):  CH2 + OH <=> CH + H2O
    // (43):  CH2 + OH <=> CH + H2O
    fwd_A[52]     = 863000;
    fwd_beta[52]  = 2.02;
    fwd_Ea[52]    = 6776;
    prefactor_units[52]  = 1.0000000000000002e-06;
    activation_units[52] = 0.50321666580471969;
    phase_units[52]      = pow(10,-12.000000);
    is_PD[52] = 0;
    nTB[52] = 0;

    // (44):  CH2 + HO2 <=> OH + CH2O
    // (44):  CH2 + HO2 <=> OH + CH2O
    fwd_A[53]     = 20000000000000;
    fwd_beta[53]  = 0;
    fwd_Ea[53]    = 0;
    prefactor_units[53]  = 1.0000000000000002e-06;
    activation_units[53] = 0.50321666580471969;
    phase_units[53]      = pow(10,-12.000000);
    is_PD[53] = 0;
    nTB[53] = 0;

    // (45):  CH2 + H2 <=> H + CH3
    // (45):  CH2 + H2 <=> H + CH3
    fwd_A[54]     = 1265000;
    fwd_beta[54]  = 2;
    fwd_Ea[54]    = 7230;
    prefactor_units[54]  = 1.0000000000000002e-06;
    activation_units[54] = 0.50321666580471969;
    phase_units[54]      = pow(10,-12.000000);
    is_PD[54] = 0;
    nTB[54] = 0;

    // (46):  CH2 + O2 => OH + H + CO
    // (46):  CH2 + O2 => OH + H + CO
    fwd_A[55]     = 2643000000000;
    fwd_beta[55]  = 0;
    fwd_Ea[55]    = 1000;
    prefactor_units[55]  = 1.0000000000000002e-06;
    activation_units[55] = 0.50321666580471969;
    phase_units[55]      = pow(10,-12.000000);
    is_PD[55] = 0;
    nTB[55] = 0;

    // (47):  CH2 + O2 => 2.000000 H + CO2
    // (47):  CH2 + O2 => 2.000000 H + CO2
    fwd_A[56]     = 1844000000000;
    fwd_beta[56]  = 0;
    fwd_Ea[56]    = 1000;
    prefactor_units[56]  = 1.0000000000000002e-06;
    activation_units[56] = 0.50321666580471969;
    phase_units[56]      = pow(10,-12.000000);
    is_PD[56] = 0;
    nTB[56] = 0;

    // (48):  CH2 + O2 <=> O + CH2O
    // (48):  CH2 + O2 <=> O + CH2O
    fwd_A[57]     = 1600000000000;
    fwd_beta[57]  = 0;
    fwd_Ea[57]    = 1000;
    prefactor_units[57]  = 1.0000000000000002e-06;
    activation_units[57] = 0.50321666580471969;
    phase_units[57]      = pow(10,-12.000000);
    is_PD[57] = 0;
    nTB[57] = 0;

    // (49):  CH2 + O2 <=> H2 + CO2
    // (49):  CH2 + O2 <=> H2 + CO2
    fwd_A[58]     = 1836000000000;
    fwd_beta[58]  = 0;
    fwd_Ea[58]    = 1000;
    prefactor_units[58]  = 1.0000000000000002e-06;
    activation_units[58] = 0.50321666580471969;
    phase_units[58]      = pow(10,-12.000000);
    is_PD[58] = 0;
    nTB[58] = 0;

    // (50):  CH2 + O2 <=> H2O + CO
    // (50):  CH2 + O2 <=> H2O + CO
    fwd_A[59]     = 520000000000;
    fwd_beta[59]  = 0;
    fwd_Ea[59]    = 1000;
    prefactor_units[59]  = 1.0000000000000002e-06;
    activation_units[59] = 0.50321666580471969;
    phase_units[59]      = pow(10,-12.000000);
    is_PD[59] = 0;
    nTB[59] = 0;

    // (51):  CH2(S) + N2 <=> CH2 + N2
    // (51):  CH2(S) + N2 <=> CH2 + N2
    fwd_A[60]     = 12000000000000;
    fwd_beta[60]  = 0;
    fwd_Ea[60]    = 471;
    prefactor_units[60]  = 1.0000000000000002e-06;
    activation_units[60] = 0.50321666580471969;
    phase_units[60]      = pow(10,-12.000000);
    is_PD[60] = 0;
    nTB[60] = 0;

    // (52):  CH2(S) + HE <=> CH2 + HE
    // (52):  CH2(S) + HE <=> CH2 + HE
    fwd_A[61]     = 6620000000000;
    fwd_beta[61]  = 0;
    fwd_Ea[61]    = 755;
    prefactor_units[61]  = 1.0000000000000002e-06;
    activation_units[61] = 0.50321666580471969;
    phase_units[61]      = pow(10,-12.000000);
    is_PD[61] = 0;
    nTB[61] = 0;

    // (53):  CH2(S) + H <=> CH + H2
    // (53):  CH2(S) + H <=> CH + H2
    fwd_A[62]     = 30000000000000;
    fwd_beta[62]  = 0;
    fwd_Ea[62]    = 0;
    prefactor_units[62]  = 1.0000000000000002e-06;
    activation_units[62] = 0.50321666580471969;
    phase_units[62]      = pow(10,-12.000000);
    is_PD[62] = 0;
    nTB[62] = 0;

    // (54):  CH2(S) + O => 2.000000 H + CO
    // (54):  CH2(S) + O => 2.000000 H + CO
    fwd_A[63]     = 30000000000000;
    fwd_beta[63]  = 0;
    fwd_Ea[63]    = 0;
    prefactor_units[63]  = 1.0000000000000002e-06;
    activation_units[63] = 0.50321666580471969;
    phase_units[63]      = pow(10,-12.000000);
    is_PD[63] = 0;
    nTB[63] = 0;

    // (55):  CH2(S) + OH <=> H + CH2O
    // (55):  CH2(S) + OH <=> H + CH2O
    fwd_A[64]     = 30000000000000;
    fwd_beta[64]  = 0;
    fwd_Ea[64]    = 0;
    prefactor_units[64]  = 1.0000000000000002e-06;
    activation_units[64] = 0.50321666580471969;
    phase_units[64]      = pow(10,-12.000000);
    is_PD[64] = 0;
    nTB[64] = 0;

    // (56):  CH2(S) + H2 <=> CH3 + H
    // (56):  CH2(S) + H2 <=> CH3 + H
    fwd_A[65]     = 82910000000000;
    fwd_beta[65]  = 0;
    fwd_Ea[65]    = 0;
    prefactor_units[65]  = 1.0000000000000002e-06;
    activation_units[65] = 0.50321666580471969;
    phase_units[65]      = pow(10,-12.000000);
    is_PD[65] = 0;
    nTB[65] = 0;

    // (57):  CH2(S) + O2 <=> CH2 + O2
    // (57):  CH2(S) + O2 <=> CH2 + O2
    fwd_A[66]     = 31300000000000;
    fwd_beta[66]  = 0;
    fwd_Ea[66]    = 0;
    prefactor_units[66]  = 1.0000000000000002e-06;
    activation_units[66] = 0.50321666580471969;
    phase_units[66]      = pow(10,-12.000000);
    is_PD[66] = 0;
    nTB[66] = 0;

    // (58):  CH2(S) + H2O (+M) <=> CH3OH (+M)
    // (58):  CH2(S) + H2O (+M) <=> CH3OH (+M)
    fwd_A[3]     = 2940000000000;
    fwd_beta[3]  = 0.052999999999999999;
    fwd_Ea[3]    = -1897;
    low_A[3]     = 1.6800000000000001e+41;
    low_beta[3]  = -7.1920000000000002;
    low_Ea[3]    = 5777;
    troe_a[3]    = 0.99199999999999999;
    troe_Tsss[3] = 943;
    troe_Ts[3]   = 47310;
    troe_Tss[3]  = 47110;
    troe_len[3]  = 4;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-12.000000);
    is_PD[3] = 1;
    nTB[3] = 6;
    TB[3] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[3] = (int *) malloc(6 * sizeof(int));
    TBid[3][0] = 6; TB[3][0] = 6; // H2O
    TBid[3][1] = 8; TB[3][1] = 1.5; // CO
    TBid[3][2] = 9; TB[3][2] = 2; // CO2
    TBid[3][3] = 14; TB[3][3] = 2; // CH4
    TBid[3][4] = 16; TB[3][4] = 2.5; // CH2O
    TBid[3][5] = 19; TB[3][5] = 3; // CH3OH

    // (59):  CH2(S) + H2O <=> CH2 + H2O
    // (59):  CH2(S) + H2O <=> CH2 + H2O
    fwd_A[67]     = 15100000000000;
    fwd_beta[67]  = 0;
    fwd_Ea[67]    = -431;
    prefactor_units[67]  = 1.0000000000000002e-06;
    activation_units[67] = 0.50321666580471969;
    phase_units[67]      = pow(10,-12.000000);
    is_PD[67] = 0;
    nTB[67] = 0;

    // (60):  CH2(S) + H2O <=> H2 + CH2O
    // (60):  CH2(S) + H2O <=> H2 + CH2O
    fwd_A[68]     = 6.6700000000000005e+21;
    fwd_beta[68]  = -3.1339999999999999;
    fwd_Ea[68]    = 3300;
    prefactor_units[68]  = 1.0000000000000002e-06;
    activation_units[68] = 0.50321666580471969;
    phase_units[68]      = pow(10,-12.000000);
    is_PD[68] = 0;
    nTB[68] = 0;

    // (61):  CH2(S) + CO <=> CH2 + CO
    // (61):  CH2(S) + CO <=> CH2 + CO
    fwd_A[69]     = 9000000000000;
    fwd_beta[69]  = 0;
    fwd_Ea[69]    = 0;
    prefactor_units[69]  = 1.0000000000000002e-06;
    activation_units[69] = 0.50321666580471969;
    phase_units[69]      = pow(10,-12.000000);
    is_PD[69] = 0;
    nTB[69] = 0;

    // (62):  CH2(S) + CO2 <=> CH2 + CO2
    // (62):  CH2(S) + CO2 <=> CH2 + CO2
    fwd_A[70]     = 13300000000000;
    fwd_beta[70]  = 0;
    fwd_Ea[70]    = 0;
    prefactor_units[70]  = 1.0000000000000002e-06;
    activation_units[70] = 0.50321666580471969;
    phase_units[70]      = pow(10,-12.000000);
    is_PD[70] = 0;
    nTB[70] = 0;

    // (63):  CH2(S) + CO2 <=> CO + CH2O
    // (63):  CH2(S) + CO2 <=> CO + CH2O
    fwd_A[71]     = 6620000000000;
    fwd_beta[71]  = 0;
    fwd_Ea[71]    = 0;
    prefactor_units[71]  = 1.0000000000000002e-06;
    activation_units[71] = 0.50321666580471969;
    phase_units[71]      = pow(10,-12.000000);
    is_PD[71] = 0;
    nTB[71] = 0;

    // (64):  HCO + H (+M) <=> CH2O (+M)
    // (64):  HCO + H (+M) <=> CH2O (+M)
    fwd_A[4]     = 191300000000000;
    fwd_beta[4]  = -0.033000000000000002;
    fwd_Ea[4]    = -142;
    low_A[4]     = 4.1900000000000004e+34;
    low_beta[4]  = -5.5330000000000004;
    low_Ea[4]    = 6128;
    troe_a[4]    = 0.78200000000000003;
    troe_Tsss[4] = 271;
    troe_Ts[4]   = 2755;
    troe_Tss[4]  = 6570;
    troe_len[4]  = 4;
    prefactor_units[4]  = 1.0000000000000002e-06;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 1;
    nTB[4] = 7;
    TB[4] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[4] = (int *) malloc(7 * sizeof(int));
    TBid[4][0] = 6; TB[4][0] = 6; // H2O
    TBid[4][1] = 8; TB[4][1] = 1.5; // CO
    TBid[4][2] = 9; TB[4][2] = 2; // CO2
    TBid[4][3] = 14; TB[4][3] = 2; // CH4
    TBid[4][4] = 16; TB[4][4] = 2.8399999999999999; // CH2O
    TBid[4][5] = 19; TB[4][5] = 3; // CH3OH
    TBid[4][6] = 20; TB[4][6] = 0.69999999999999996; // HE

    // (65):  CH2O (+M) <=> H2 + CO (+M)
    // (65):  CH2O (+M) <=> H2 + CO (+M)
    fwd_A[5]     = 37000000000000;
    fwd_beta[5]  = 0;
    fwd_Ea[5]    = 71976;
    low_A[5]     = 4.4000000000000001e+38;
    low_beta[5]  = -6.0999999999999996;
    low_Ea[5]    = 94000;
    troe_a[5]    = 0.93200000000000005;
    troe_Tsss[5] = 197;
    troe_Ts[5]   = 1540;
    troe_Tss[5]  = 10300;
    troe_len[5]  = 4;
    prefactor_units[5]  = 1;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-6.000000);
    is_PD[5] = 1;
    nTB[5] = 7;
    TB[5] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[5] = (int *) malloc(7 * sizeof(int));
    TBid[5][0] = 6; TB[5][0] = 6; // H2O
    TBid[5][1] = 8; TB[5][1] = 1.5; // CO
    TBid[5][2] = 9; TB[5][2] = 2; // CO2
    TBid[5][3] = 14; TB[5][3] = 2; // CH4
    TBid[5][4] = 16; TB[5][4] = 2.5; // CH2O
    TBid[5][5] = 19; TB[5][5] = 3; // CH3OH
    TBid[5][6] = 20; TB[5][6] = 0.69999999999999996; // HE

    // (66):  CH2O + H <=> HCO + H2
    // (66):  CH2O + H <=> HCO + H2
    fwd_A[72]     = 71490000;
    fwd_beta[72]  = 1.8999999999999999;
    fwd_Ea[72]    = 2742;
    prefactor_units[72]  = 1.0000000000000002e-06;
    activation_units[72] = 0.50321666580471969;
    phase_units[72]      = pow(10,-12.000000);
    is_PD[72] = 0;
    nTB[72] = 0;

    // (67):  CH2O + O <=> OH + HCO
    // (67):  CH2O + O <=> OH + HCO
    fwd_A[73]     = 424400000000;
    fwd_beta[73]  = 0.56999999999999995;
    fwd_Ea[73]    = 2762;
    prefactor_units[73]  = 1.0000000000000002e-06;
    activation_units[73] = 0.50321666580471969;
    phase_units[73]      = pow(10,-12.000000);
    is_PD[73] = 0;
    nTB[73] = 0;

    // (68):  CH2O + OH <=> HCO + H2O
    // (68):  CH2O + OH <=> HCO + H2O
    fwd_A[74]     = 83380000;
    fwd_beta[74]  = 1.6299999999999999;
    fwd_Ea[74]    = -1055;
    prefactor_units[74]  = 1.0000000000000002e-06;
    activation_units[74] = 0.50321666580471969;
    phase_units[74]      = pow(10,-12.000000);
    is_PD[74] = 0;
    nTB[74] = 0;

    // (69):  CH2O + O2 <=> HO2 + HCO
    // (69):  CH2O + O2 <=> HO2 + HCO
    fwd_A[75]     = 329700;
    fwd_beta[75]  = 2.5;
    fwd_Ea[75]    = 36460;
    prefactor_units[75]  = 1.0000000000000002e-06;
    activation_units[75] = 0.50321666580471969;
    phase_units[75]      = pow(10,-12.000000);
    is_PD[75] = 0;
    nTB[75] = 0;

    // (70):  CH2O + CH2 <=> CH3 + HCO
    // (70):  CH2O + CH2 <=> CH3 + HCO
    fwd_A[76]     = 0.073999999999999996;
    fwd_beta[76]  = 4.21;
    fwd_Ea[76]    = 1120;
    prefactor_units[76]  = 1.0000000000000002e-06;
    activation_units[76] = 0.50321666580471969;
    phase_units[76]      = pow(10,-12.000000);
    is_PD[76] = 0;
    nTB[76] = 0;

    // (71):  CH2O + CH2(S) <=> CH3 + HCO
    // (71):  CH2O + CH2(S) <=> CH3 + HCO
    fwd_A[77]     = 13300000000000;
    fwd_beta[77]  = 0;
    fwd_Ea[77]    = -550;
    prefactor_units[77]  = 1.0000000000000002e-06;
    activation_units[77] = 0.50321666580471969;
    phase_units[77]      = pow(10,-12.000000);
    is_PD[77] = 0;
    nTB[77] = 0;

    // (72):  CH3 + H (+M) <=> CH4 (+M)
    // (72):  CH3 + H (+M) <=> CH4 (+M)
    fwd_A[6]     = 180100000000000;
    fwd_beta[6]  = 0;
    fwd_Ea[6]    = 0;
    low_A[6]     = 7.9300000000000004e+24;
    low_beta[6]  = -2.1699999999999999;
    low_Ea[6]    = 0;
    troe_a[6]    = 0.124;
    troe_Tsss[6] = 1801;
    troe_Ts[6]   = 33.100000000000001;
    troe_Tss[6]  = 90000;
    troe_len[6]  = 4;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = pow(10,-12.000000);
    is_PD[6] = 1;
    nTB[6] = 9;
    TB[6] = (amrex::Real *) malloc(9 * sizeof(amrex::Real));
    TBid[6] = (int *) malloc(9 * sizeof(int));
    TBid[6][0] = 0; TB[6][0] = 0.58999999999999997; // N2
    TBid[6][1] = 4; TB[6][1] = 0.58999999999999997; // O2
    TBid[6][2] = 6; TB[6][2] = 3.4199999999999999; // H2O
    TBid[6][3] = 8; TB[6][3] = 0.89000000000000001; // CO
    TBid[6][4] = 9; TB[6][4] = 2; // CO2
    TBid[6][5] = 14; TB[6][5] = 2; // CH4
    TBid[6][6] = 16; TB[6][6] = 2.5; // CH2O
    TBid[6][7] = 19; TB[6][7] = 3; // CH3OH
    TBid[6][8] = 20; TB[6][8] = 0.41999999999999998; // HE

    // (73):  CH3 + O <=> H + CH2O
    // (73):  CH3 + O <=> H + CH2O
    fwd_A[78]     = 57220000000000;
    fwd_beta[78]  = 0;
    fwd_Ea[78]    = 0;
    prefactor_units[78]  = 1.0000000000000002e-06;
    activation_units[78] = 0.50321666580471969;
    phase_units[78]      = pow(10,-12.000000);
    is_PD[78] = 0;
    nTB[78] = 0;

    // (74):  CH3 + O => H + H2 + CO
    // (74):  CH3 + O => H + H2 + CO
    fwd_A[79]     = 23840000000000;
    fwd_beta[79]  = 0;
    fwd_Ea[79]    = 0;
    prefactor_units[79]  = 1.0000000000000002e-06;
    activation_units[79] = 0.50321666580471969;
    phase_units[79]      = pow(10,-12.000000);
    is_PD[79] = 0;
    nTB[79] = 0;

    // (75):  CH3 + OH (+M) <=> CH3OH (+M)
    // (75):  CH3 + OH (+M) <=> CH3OH (+M)
    fwd_A[7]     = 62100000000000;
    fwd_beta[7]  = -0.017999999999999999;
    fwd_Ea[7]    = -33;
    low_A[7]     = 7.2399999999999995e+36;
    low_beta[7]  = -6;
    low_Ea[7]    = 3226;
    troe_a[7]    = 0.1855;
    troe_Tsss[7] = 156;
    troe_Ts[7]   = 1675;
    troe_Tss[7]  = 4530;
    troe_len[7]  = 4;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = pow(10,-12.000000);
    is_PD[7] = 1;
    nTB[7] = 6;
    TB[7] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[7] = (int *) malloc(6 * sizeof(int));
    TBid[7][0] = 6; TB[7][0] = 6; // H2O
    TBid[7][1] = 8; TB[7][1] = 1.5; // CO
    TBid[7][2] = 9; TB[7][2] = 2; // CO2
    TBid[7][3] = 14; TB[7][3] = 2; // CH4
    TBid[7][4] = 16; TB[7][4] = 2.5; // CH2O
    TBid[7][5] = 19; TB[7][5] = 3; // CH3OH

    // (76):  CH3 + OH <=> CH2 + H2O
    // (76):  CH3 + OH <=> CH2 + H2O
    fwd_A[80]     = 44640;
    fwd_beta[80]  = 2.5699999999999998;
    fwd_Ea[80]    = 3998;
    prefactor_units[80]  = 1.0000000000000002e-06;
    activation_units[80] = 0.50321666580471969;
    phase_units[80]      = pow(10,-12.000000);
    is_PD[80] = 0;
    nTB[80] = 0;

    // (77):  CH3 + OH <=> CH2(S) + H2O
    // (77):  CH3 + OH <=> CH2(S) + H2O
    fwd_A[81]     = 7810000000000000;
    fwd_beta[81]  = -0.91000000000000003;
    fwd_Ea[81]    = 546;
    prefactor_units[81]  = 1.0000000000000002e-06;
    activation_units[81] = 0.50321666580471969;
    phase_units[81]      = pow(10,-12.000000);
    is_PD[81] = 0;
    nTB[81] = 0;

    // (78):  CH3 + OH <=> H2 + CH2O
    // (78):  CH3 + OH <=> H2 + CH2O
    fwd_A[82]     = 2735000000;
    fwd_beta[82]  = 0.73399999999999999;
    fwd_Ea[82]    = -2177;
    prefactor_units[82]  = 1.0000000000000002e-06;
    activation_units[82] = 0.50321666580471969;
    phase_units[82]      = pow(10,-12.000000);
    is_PD[82] = 0;
    nTB[82] = 0;

    // (79):  CH3 + HO2 <=> O2 + CH4
    // (79):  CH3 + HO2 <=> O2 + CH4
    fwd_A[83]     = 126900;
    fwd_beta[83]  = 2.2280000000000002;
    fwd_Ea[83]    = -3022;
    prefactor_units[83]  = 1.0000000000000002e-06;
    activation_units[83] = 0.50321666580471969;
    phase_units[83]      = pow(10,-12.000000);
    is_PD[83] = 0;
    nTB[83] = 0;

    // (80):  CH3 + HO2 <=> OH + CH3O
    // (80):  CH3 + HO2 <=> OH + CH3O
    fwd_A[84]     = 8821000000000;
    fwd_beta[84]  = 0;
    fwd_Ea[84]    = -590;
    prefactor_units[84]  = 1.0000000000000002e-06;
    activation_units[84] = 0.50321666580471969;
    phase_units[84]      = pow(10,-12.000000);
    is_PD[84] = 0;
    nTB[84] = 0;

    // (81):  CH3 + O2 <=> O + CH3O
    // (81):  CH3 + O2 <=> O + CH3O
    fwd_A[85]     = 8104000000000;
    fwd_beta[85]  = 0;
    fwd_Ea[85]    = 28297;
    prefactor_units[85]  = 1.0000000000000002e-06;
    activation_units[85] = 0.50321666580471969;
    phase_units[85]      = pow(10,-12.000000);
    is_PD[85] = 0;
    nTB[85] = 0;

    // (82):  CH3 + O2 <=> OH + CH2O
    // (82):  CH3 + O2 <=> OH + CH2O
    fwd_A[86]     = 99.769999999999996;
    fwd_beta[86]  = 2.8599999999999999;
    fwd_Ea[86]    = 9768;
    prefactor_units[86]  = 1.0000000000000002e-06;
    activation_units[86] = 0.50321666580471969;
    phase_units[86]      = pow(10,-12.000000);
    is_PD[86] = 0;
    nTB[86] = 0;

    // (83):  CH3 + HCO <=> CH4 + CO
    // (83):  CH3 + HCO <=> CH4 + CO
    fwd_A[87]     = 5300000000000;
    fwd_beta[87]  = 0;
    fwd_Ea[87]    = 0;
    prefactor_units[87]  = 1.0000000000000002e-06;
    activation_units[87] = 0.50321666580471969;
    phase_units[87]      = pow(10,-12.000000);
    is_PD[87] = 0;
    nTB[87] = 0;

    // (84):  CH3 + CH2O <=> HCO + CH4
    // (84):  CH3 + CH2O <=> HCO + CH4
    fwd_A[88]     = 32.130000000000003;
    fwd_beta[88]  = 3.3599999999999999;
    fwd_Ea[88]    = 4310;
    prefactor_units[88]  = 1.0000000000000002e-06;
    activation_units[88] = 0.50321666580471969;
    phase_units[88]      = pow(10,-12.000000);
    is_PD[88] = 0;
    nTB[88] = 0;

    // (85):  CH3O (+M) <=> H + CH2O (+M)
    // (85):  CH3O (+M) <=> H + CH2O (+M)
    fwd_A[8]     = 11300000000;
    fwd_beta[8]  = 1.21;
    fwd_Ea[8]    = 24075;
    low_A[8]     = 60200000000000000;
    low_beta[8]  = -0.54700000000000004;
    low_Ea[8]    = 18024;
    troe_a[8]    = 0.34100000000000003;
    troe_Tsss[8] = 28;
    troe_Ts[8]   = 1000;
    troe_Tss[8]  = 2339;
    troe_len[8]  = 4;
    prefactor_units[8]  = 1;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = pow(10,-6.000000);
    is_PD[8] = 1;
    nTB[8] = 8;
    TB[8] = (amrex::Real *) malloc(8 * sizeof(amrex::Real));
    TBid[8] = (int *) malloc(8 * sizeof(int));
    TBid[8][0] = 1; TB[8][0] = 2; // H2
    TBid[8][1] = 6; TB[8][1] = 6; // H2O
    TBid[8][2] = 8; TB[8][2] = 1.5; // CO
    TBid[8][3] = 9; TB[8][3] = 2; // CO2
    TBid[8][4] = 14; TB[8][4] = 2; // CH4
    TBid[8][5] = 16; TB[8][5] = 2.5; // CH2O
    TBid[8][6] = 19; TB[8][6] = 3; // CH3OH
    TBid[8][7] = 20; TB[8][7] = 0.67000000000000004; // HE

    // (86):  CH3O + H (+M) <=> CH3OH (+M)
    // (86):  CH3O + H (+M) <=> CH3OH (+M)
    fwd_A[9]     = 244000000000;
    fwd_beta[9]  = 0.76000000000000001;
    fwd_Ea[9]    = 0;
    low_A[9]     = 6.7000000000000002e+40;
    low_beta[9]  = -7.3799999999999999;
    low_Ea[9]    = 9177;
    troe_a[9]    = 0.68400000000000005;
    troe_Tsss[9] = 37050;
    troe_Ts[9]   = 41490;
    troe_Tss[9]  = 3980;
    troe_len[9]  = 4;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = pow(10,-12.000000);
    is_PD[9] = 1;
    nTB[9] = 6;
    TB[9] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[9] = (int *) malloc(6 * sizeof(int));
    TBid[9][0] = 6; TB[9][0] = 6; // H2O
    TBid[9][1] = 8; TB[9][1] = 1.5; // CO
    TBid[9][2] = 9; TB[9][2] = 2; // CO2
    TBid[9][3] = 14; TB[9][3] = 2; // CH4
    TBid[9][4] = 16; TB[9][4] = 2.5; // CH2O
    TBid[9][5] = 19; TB[9][5] = 3; // CH3OH

    // (87):  CH3O + H <=> H + CH2OH
    // (87):  CH3O + H <=> H + CH2OH
    fwd_A[89]     = 12900000;
    fwd_beta[89]  = 1.8200000000000001;
    fwd_Ea[89]    = -703;
    prefactor_units[89]  = 1.0000000000000002e-06;
    activation_units[89] = 0.50321666580471969;
    phase_units[89]      = pow(10,-12.000000);
    is_PD[89] = 0;
    nTB[89] = 0;

    // (88):  CH3O + H <=> H2 + CH2O
    // (88):  CH3O + H <=> H2 + CH2O
    fwd_A[90]     = 37900000000000;
    fwd_beta[90]  = 0;
    fwd_Ea[90]    = 596;
    prefactor_units[90]  = 1.0000000000000002e-06;
    activation_units[90] = 0.50321666580471969;
    phase_units[90]      = pow(10,-12.000000);
    is_PD[90] = 0;
    nTB[90] = 0;

    // (89):  CH3O + H <=> OH + CH3
    // (89):  CH3O + H <=> OH + CH3
    fwd_A[91]     = 388000000000000;
    fwd_beta[91]  = -0.26400000000000001;
    fwd_Ea[91]    = -26;
    prefactor_units[91]  = 1.0000000000000002e-06;
    activation_units[91] = 0.50321666580471969;
    phase_units[91]      = pow(10,-12.000000);
    is_PD[91] = 0;
    nTB[91] = 0;

    // (90):  CH3O + H <=> CH2(S) + H2O
    // (90):  CH3O + H <=> CH2(S) + H2O
    fwd_A[92]     = 197000000000;
    fwd_beta[92]  = 0.41399999999999998;
    fwd_Ea[92]    = 243;
    prefactor_units[92]  = 1.0000000000000002e-06;
    activation_units[92] = 0.50321666580471969;
    phase_units[92]      = pow(10,-12.000000);
    is_PD[92] = 0;
    nTB[92] = 0;

    // (91):  CH3O + O <=> OH + CH2O
    // (91):  CH3O + O <=> OH + CH2O
    fwd_A[93]     = 3780000000000;
    fwd_beta[93]  = 0;
    fwd_Ea[93]    = 0;
    prefactor_units[93]  = 1.0000000000000002e-06;
    activation_units[93] = 0.50321666580471969;
    phase_units[93]      = pow(10,-12.000000);
    is_PD[93] = 0;
    nTB[93] = 0;

    // (92):  CH3O + OH <=> H2O + CH2O
    // (92):  CH3O + OH <=> H2O + CH2O
    fwd_A[94]     = 18100000000000;
    fwd_beta[94]  = 0;
    fwd_Ea[94]    = 0;
    prefactor_units[94]  = 1.0000000000000002e-06;
    activation_units[94] = 0.50321666580471969;
    phase_units[94]      = pow(10,-12.000000);
    is_PD[94] = 0;
    nTB[94] = 0;

    // (93):  CH3O + O2 <=> HO2 + CH2O
    // (93):  CH3O + O2 <=> HO2 + CH2O
    fwd_A[95]     = 63200000000;
    fwd_beta[95]  = 0;
    fwd_Ea[95]    = 2603;
    prefactor_units[95]  = 1.0000000000000002e-06;
    activation_units[95] = 0.50321666580471969;
    phase_units[95]      = pow(10,-12.000000);
    is_PD[95] = 0;
    nTB[95] = 0;

    // (94):  CH3O + CH3 <=> CH4 + CH2O
    // (94):  CH3O + CH3 <=> CH4 + CH2O
    fwd_A[96]     = 24000000000000;
    fwd_beta[96]  = 0;
    fwd_Ea[96]    = 0;
    prefactor_units[96]  = 1.0000000000000002e-06;
    activation_units[96] = 0.50321666580471969;
    phase_units[96]      = pow(10,-12.000000);
    is_PD[96] = 0;
    nTB[96] = 0;

    // (95):  CH3O + CO <=> CH3 + CO2
    // (95):  CH3O + CO <=> CH3 + CO2
    fwd_A[97]     = 6000000000000;
    fwd_beta[97]  = 0;
    fwd_Ea[97]    = 11000;
    prefactor_units[97]  = 1.0000000000000002e-06;
    activation_units[97] = 0.50321666580471969;
    phase_units[97]      = pow(10,-12.000000);
    is_PD[97] = 0;
    nTB[97] = 0;

    // (96):  CH2OH (+M) <=> H + CH2O (+M)
    // (96):  CH2OH (+M) <=> H + CH2O (+M)
    fwd_A[10]     = 73700000000;
    fwd_beta[10]  = 0.81100000000000005;
    fwd_Ea[10]    = 39580;
    low_A[10]     = 30100000000000;
    low_beta[10]  = 0.184;
    low_Ea[10]    = 17230;
    troe_a[10]    = 0.001;
    troe_Tsss[10] = 50;
    troe_Ts[10]   = 600;
    troe_Tss[10]  = 2780;
    troe_len[10]  = 4;
    prefactor_units[10]  = 1;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = pow(10,-6.000000);
    is_PD[10] = 1;
    nTB[10] = 8;
    TB[10] = (amrex::Real *) malloc(8 * sizeof(amrex::Real));
    TBid[10] = (int *) malloc(8 * sizeof(int));
    TBid[10][0] = 1; TB[10][0] = 2; // H2
    TBid[10][1] = 6; TB[10][1] = 5.9699999999999998; // H2O
    TBid[10][2] = 8; TB[10][2] = 1.5; // CO
    TBid[10][3] = 9; TB[10][3] = 2; // CO2
    TBid[10][4] = 14; TB[10][4] = 2; // CH4
    TBid[10][5] = 16; TB[10][5] = 2.5; // CH2O
    TBid[10][6] = 19; TB[10][6] = 3; // CH3OH
    TBid[10][7] = 20; TB[10][7] = 0.67000000000000004; // HE

    // (97):  CH2OH + H (+M) <=> CH3OH (+M)
    // (97):  CH2OH + H (+M) <=> CH3OH (+M)
    fwd_A[11]     = 66700000000;
    fwd_beta[11]  = 0.95999999999999996;
    fwd_Ea[11]    = 0;
    low_A[11]     = 1.34e+41;
    low_beta[11]  = -7.3799999999999999;
    low_Ea[11]    = 9177;
    troe_a[11]    = 0.68400000000000005;
    troe_Tsss[11] = 37050;
    troe_Ts[11]   = 41490;
    troe_Tss[11]  = 3980;
    troe_len[11]  = 4;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = pow(10,-12.000000);
    is_PD[11] = 1;
    nTB[11] = 7;
    TB[11] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[11] = (int *) malloc(7 * sizeof(int));
    TBid[11][0] = 6; TB[11][0] = 6; // H2O
    TBid[11][1] = 8; TB[11][1] = 1.5; // CO
    TBid[11][2] = 9; TB[11][2] = 2; // CO2
    TBid[11][3] = 14; TB[11][3] = 2; // CH4
    TBid[11][4] = 16; TB[11][4] = 2.5; // CH2O
    TBid[11][5] = 19; TB[11][5] = 3; // CH3OH
    TBid[11][6] = 20; TB[11][6] = 0.69999999999999996; // HE

    // (98):  CH2OH + H <=> H2 + CH2O
    // (98):  CH2OH + H <=> H2 + CH2O
    fwd_A[98]     = 24400000000000;
    fwd_beta[98]  = 0;
    fwd_Ea[98]    = 0;
    prefactor_units[98]  = 1.0000000000000002e-06;
    activation_units[98] = 0.50321666580471969;
    phase_units[98]      = pow(10,-12.000000);
    is_PD[98] = 0;
    nTB[98] = 0;

    // (99):  CH2OH + H <=> OH + CH3
    // (99):  CH2OH + H <=> OH + CH3
    fwd_A[99]     = 20060000000000;
    fwd_beta[99]  = 0.19800000000000001;
    fwd_Ea[99]    = -241;
    prefactor_units[99]  = 1.0000000000000002e-06;
    activation_units[99] = 0.50321666580471969;
    phase_units[99]      = pow(10,-12.000000);
    is_PD[99] = 0;
    nTB[99] = 0;

    // (100):  CH2OH + H <=> CH2(S) + H2O
    // (100):  CH2OH + H <=> CH2(S) + H2O
    fwd_A[100]     = 128000000000;
    fwd_beta[100]  = 0.51600000000000001;
    fwd_Ea[100]    = 215;
    prefactor_units[100]  = 1.0000000000000002e-06;
    activation_units[100] = 0.50321666580471969;
    phase_units[100]      = pow(10,-12.000000);
    is_PD[100] = 0;
    nTB[100] = 0;

    // (101):  CH2OH + O <=> OH + CH2O
    // (101):  CH2OH + O <=> OH + CH2O
    fwd_A[101]     = 90300000000000;
    fwd_beta[101]  = 0;
    fwd_Ea[101]    = 0;
    prefactor_units[101]  = 1.0000000000000002e-06;
    activation_units[101] = 0.50321666580471969;
    phase_units[101]      = pow(10,-12.000000);
    is_PD[101] = 0;
    nTB[101] = 0;

    // (102):  CH2OH + OH <=> H2O + CH2O
    // (102):  CH2OH + OH <=> H2O + CH2O
    fwd_A[102]     = 24100000000000;
    fwd_beta[102]  = 0;
    fwd_Ea[102]    = 0;
    prefactor_units[102]  = 1.0000000000000002e-06;
    activation_units[102] = 0.50321666580471969;
    phase_units[102]      = pow(10,-12.000000);
    is_PD[102] = 0;
    nTB[102] = 0;

    // (103):  CH2OH + O2 <=> HO2 + CH2O
    // (103):  CH2OH + O2 <=> HO2 + CH2O
    fwd_A[103]     = 72980000000000;
    fwd_beta[103]  = 0;
    fwd_Ea[103]    = 3736;
    prefactor_units[103]  = 1.0000000000000002e-06;
    activation_units[103] = 0.50321666580471969;
    phase_units[103]      = pow(10,-12.000000);
    is_PD[103] = 0;
    nTB[103] = 0;

    // (104):  CH2OH + CH3 <=> CH4 + CH2O
    // (104):  CH2OH + CH3 <=> CH4 + CH2O
    fwd_A[104]     = 24000000000000;
    fwd_beta[104]  = 0;
    fwd_Ea[104]    = 0;
    prefactor_units[104]  = 1.0000000000000002e-06;
    activation_units[104] = 0.50321666580471969;
    phase_units[104]      = pow(10,-12.000000);
    is_PD[104] = 0;
    nTB[104] = 0;

    // (105):  CH4 + H <=> CH3 + H2
    // (105):  CH4 + H <=> CH3 + H2
    fwd_A[105]     = 478100;
    fwd_beta[105]  = 2.5;
    fwd_Ea[105]    = 9588;
    prefactor_units[105]  = 1.0000000000000002e-06;
    activation_units[105] = 0.50321666580471969;
    phase_units[105]      = pow(10,-12.000000);
    is_PD[105] = 0;
    nTB[105] = 0;

    // (106):  CH4 + O <=> OH + CH3
    // (106):  CH4 + O <=> OH + CH3
    fwd_A[106]     = 678600000;
    fwd_beta[106]  = 1.5600000000000001;
    fwd_Ea[106]    = 8485;
    prefactor_units[106]  = 1.0000000000000002e-06;
    activation_units[106] = 0.50321666580471969;
    phase_units[106]      = pow(10,-12.000000);
    is_PD[106] = 0;
    nTB[106] = 0;

    // (107):  CH4 + OH <=> CH3 + H2O
    // (107):  CH4 + OH <=> CH3 + H2O
    fwd_A[107]     = 983900;
    fwd_beta[107]  = 2.1819999999999999;
    fwd_Ea[107]    = 2446;
    prefactor_units[107]  = 1.0000000000000002e-06;
    activation_units[107] = 0.50321666580471969;
    phase_units[107]      = pow(10,-12.000000);
    is_PD[107] = 0;
    nTB[107] = 0;

    // (108):  CH4 + CH2 <=> 2.000000 CH3
    // (108):  CH4 + CH2 <=> 2.000000 CH3
    fwd_A[108]     = 2483000;
    fwd_beta[108]  = 2;
    fwd_Ea[108]    = 8270;
    prefactor_units[108]  = 1.0000000000000002e-06;
    activation_units[108] = 0.50321666580471969;
    phase_units[108]      = pow(10,-12.000000);
    is_PD[108] = 0;
    nTB[108] = 0;

    // (109):  CH4 + CH2(S) <=> 2.000000 CH3
    // (109):  CH4 + CH2(S) <=> 2.000000 CH3
    fwd_A[109]     = 18670000000000;
    fwd_beta[109]  = 0;
    fwd_Ea[109]    = -497;
    prefactor_units[109]  = 1.0000000000000002e-06;
    activation_units[109] = 0.50321666580471969;
    phase_units[109]      = pow(10,-12.000000);
    is_PD[109] = 0;
    nTB[109] = 0;

    // (110):  CH3OH + H <=> CH2OH + H2
    // (110):  CH3OH + H <=> CH2OH + H2
    fwd_A[110]     = 1550000;
    fwd_beta[110]  = 2.351;
    fwd_Ea[110]    = 5912;
    prefactor_units[110]  = 1.0000000000000002e-06;
    activation_units[110] = 0.50321666580471969;
    phase_units[110]      = pow(10,-12.000000);
    is_PD[110] = 0;
    nTB[110] = 0;

    // (111):  CH3OH + H <=> CH3O + H2
    // (111):  CH3OH + H <=> CH3O + H2
    fwd_A[111]     = 5490000;
    fwd_beta[111]  = 2.1469999999999998;
    fwd_Ea[111]    = 11134;
    prefactor_units[111]  = 1.0000000000000002e-06;
    activation_units[111] = 0.50321666580471969;
    phase_units[111]      = pow(10,-12.000000);
    is_PD[111] = 0;
    nTB[111] = 0;

    // (112):  CH3OH + O <=> OH + CH2OH
    // (112):  CH3OH + O <=> OH + CH2OH
    fwd_A[112]     = 24700000000000;
    fwd_beta[112]  = 0;
    fwd_Ea[112]    = 5306;
    prefactor_units[112]  = 1.0000000000000002e-06;
    activation_units[112] = 0.50321666580471969;
    phase_units[112]      = pow(10,-12.000000);
    is_PD[112] = 0;
    nTB[112] = 0;

    // (113):  CH3OH + O <=> OH + CH3O
    // (113):  CH3OH + O <=> OH + CH3O
    fwd_A[113]     = 8200000000000;
    fwd_beta[113]  = 0;
    fwd_Ea[113]    = 9040;
    prefactor_units[113]  = 1.0000000000000002e-06;
    activation_units[113] = 0.50321666580471969;
    phase_units[113]      = pow(10,-12.000000);
    is_PD[113] = 0;
    nTB[113] = 0;

    // (114):  CH3OH + OH <=> CH2OH + H2O
    // (114):  CH3OH + OH <=> CH2OH + H2O
    fwd_A[114]     = 142000;
    fwd_beta[114]  = 2.3700000000000001;
    fwd_Ea[114]    = -965.20000000000005;
    prefactor_units[114]  = 1.0000000000000002e-06;
    activation_units[114] = 0.50321666580471969;
    phase_units[114]      = pow(10,-12.000000);
    is_PD[114] = 0;
    nTB[114] = 0;

    // (115):  CH3OH + OH <=> CH3O + H2O
    // (115):  CH3OH + OH <=> CH3O + H2O
    fwd_A[115]     = 16000;
    fwd_beta[115]  = 2.7000000000000002;
    fwd_Ea[115]    = 53.299999999999997;
    prefactor_units[115]  = 1.0000000000000002e-06;
    activation_units[115] = 0.50321666580471969;
    phase_units[115]      = pow(10,-12.000000);
    is_PD[115] = 0;
    nTB[115] = 0;

    // (116):  CH3OH + O2 <=> CH2OH + HO2
    // (116):  CH3OH + O2 <=> CH2OH + HO2
    fwd_A[116]     = 358000;
    fwd_beta[116]  = 2.27;
    fwd_Ea[116]    = 42760;
    prefactor_units[116]  = 1.0000000000000002e-06;
    activation_units[116] = 0.50321666580471969;
    phase_units[116]      = pow(10,-12.000000);
    is_PD[116] = 0;
    nTB[116] = 0;

    // (117):  CH3OH + CH <=> CH3 + CH2O
    // (117):  CH3OH + CH <=> CH3 + CH2O
    fwd_A[117]     = 9.04e+18;
    fwd_beta[117]  = -1.9299999999999999;
    fwd_Ea[117]    = 0;
    prefactor_units[117]  = 1.0000000000000002e-06;
    activation_units[117] = 0.50321666580471969;
    phase_units[117]      = pow(10,-12.000000);
    is_PD[117] = 0;
    nTB[117] = 0;

    // (118):  CH3OH + CH2 <=> CH3 + CH2OH
    // (118):  CH3OH + CH2 <=> CH3 + CH2OH
    fwd_A[118]     = 32;
    fwd_beta[118]  = 3.2000000000000002;
    fwd_Ea[118]    = 7175;
    prefactor_units[118]  = 1.0000000000000002e-06;
    activation_units[118] = 0.50321666580471969;
    phase_units[118]      = pow(10,-12.000000);
    is_PD[118] = 0;
    nTB[118] = 0;

    // (119):  CH3OH + CH2 <=> CH3 + CH3O
    // (119):  CH3OH + CH2 <=> CH3 + CH3O
    fwd_A[119]     = 14.5;
    fwd_beta[119]  = 3.1000000000000001;
    fwd_Ea[119]    = 6940;
    prefactor_units[119]  = 1.0000000000000002e-06;
    activation_units[119] = 0.50321666580471969;
    phase_units[119]      = pow(10,-12.000000);
    is_PD[119] = 0;
    nTB[119] = 0;

    // (120):  CH3OH + CH2(S) <=> CH3 + CH3O
    // (120):  CH3OH + CH2(S) <=> CH3 + CH3O
    fwd_A[120]     = 7000000000000;
    fwd_beta[120]  = 0;
    fwd_Ea[120]    = -550;
    prefactor_units[120]  = 1.0000000000000002e-06;
    activation_units[120] = 0.50321666580471969;
    phase_units[120]      = pow(10,-12.000000);
    is_PD[120] = 0;
    nTB[120] = 0;

    // (121):  CH3OH + CH2(S) <=> CH3 + CH2OH
    // (121):  CH3OH + CH2(S) <=> CH3 + CH2OH
    fwd_A[121]     = 20000000000000;
    fwd_beta[121]  = 0;
    fwd_Ea[121]    = -550;
    prefactor_units[121]  = 1.0000000000000002e-06;
    activation_units[121] = 0.50321666580471969;
    phase_units[121]      = pow(10,-12.000000);
    is_PD[121] = 0;
    nTB[121] = 0;

    // (122):  CH3OH + CH3 <=> CH2OH + CH4
    // (122):  CH3OH + CH3 <=> CH2OH + CH4
    fwd_A[122]     = 665;
    fwd_beta[122]  = 3.0299999999999998;
    fwd_Ea[122]    = 8720;
    prefactor_units[122]  = 1.0000000000000002e-06;
    activation_units[122] = 0.50321666580471969;
    phase_units[122]      = pow(10,-12.000000);
    is_PD[122] = 0;
    nTB[122] = 0;

    // (123):  CH3OH + CH3 <=> CH3O + CH4
    // (123):  CH3OH + CH3 <=> CH3O + CH4
    fwd_A[123]     = 21500;
    fwd_beta[123]  = 2.27;
    fwd_Ea[123]    = 8710;
    prefactor_units[123]  = 1.0000000000000002e-06;
    activation_units[123] = 0.50321666580471969;
    phase_units[123]      = pow(10,-12.000000);
    is_PD[123] = 0;
    nTB[123] = 0;

}


/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<124; ++i) {
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
    awt[4] = 4.002600; /*HE */

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
    int kd = 5; 
    /*Zero ncf */
    for (id = 0; id < kd * 21; ++ id) {
         ncf[id] = 0; 
    }

    /*N2 */
    ncf[ 0 * kd + 3 ] = 2; /*N */

    /*H2 */
    ncf[ 1 * kd + 1 ] = 2; /*H */

    /*H */
    ncf[ 2 * kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 3 * kd + 0 ] = 1; /*O */

    /*O2 */
    ncf[ 4 * kd + 0 ] = 2; /*O */

    /*OH */
    ncf[ 5 * kd + 0 ] = 1; /*O */
    ncf[ 5 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 6 * kd + 1 ] = 2; /*H */
    ncf[ 6 * kd + 0 ] = 1; /*O */

    /*HO2 */
    ncf[ 7 * kd + 1 ] = 1; /*H */
    ncf[ 7 * kd + 0 ] = 2; /*O */

    /*CO */
    ncf[ 8 * kd + 2 ] = 1; /*C */
    ncf[ 8 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 0 ] = 2; /*O */

    /*CH */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 1; /*H */

    /*CH2 */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 1 ] = 2; /*H */

    /*CH2(S) */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd + 1 ] = 4; /*H */

    /*HCO */
    ncf[ 15 * kd + 2 ] = 1; /*C */
    ncf[ 15 * kd + 1 ] = 1; /*H */
    ncf[ 15 * kd + 0 ] = 1; /*O */

    /*CH2O */
    ncf[ 16 * kd + 1 ] = 2; /*H */
    ncf[ 16 * kd + 2 ] = 1; /*C */
    ncf[ 16 * kd + 0 ] = 1; /*O */

    /*CH2OH */
    ncf[ 17 * kd + 2 ] = 1; /*C */
    ncf[ 17 * kd + 1 ] = 3; /*H */
    ncf[ 17 * kd + 0 ] = 1; /*O */

    /*CH3O */
    ncf[ 18 * kd + 2 ] = 1; /*C */
    ncf[ 18 * kd + 1 ] = 3; /*H */
    ncf[ 18 * kd + 0 ] = 1; /*O */

    /*CH3OH */
    ncf[ 19 * kd + 2 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 4; /*H */
    ncf[ 19 * kd + 0 ] = 1; /*O */

    /*HE */
    ncf[ 20 * kd + 4 ] = 1; /*HE */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(5);
    ename[0] = "O";
    ename[1] = "H";
    ename[2] = "C";
    ename[3] = "N";
    ename[4] = "HE";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(21);
    kname[0] = "N2";
    kname[1] = "H2";
    kname[2] = "H";
    kname[3] = "O";
    kname[4] = "O2";
    kname[5] = "OH";
    kname[6] = "H2O";
    kname[7] = "HO2";
    kname[8] = "CO";
    kname[9] = "CO2";
    kname[10] = "CH";
    kname[11] = "CH2";
    kname[12] = "CH2(S)";
    kname[13] = "CH3";
    kname[14] = "CH4";
    kname[15] = "HCO";
    kname[16] = "CH2O";
    kname[17] = "CH2OH";
    kname[18] = "CH3O";
    kname[19] = "CH3OH";
    kname[20] = "HE";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<21; l++) {
                c_d[l] = 1.0/ 21.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<22; k++) {
        for (int l=0; l<22; l++) {
            if(J_h[ 22 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<22; k++) {
        for (int l=0; l<22; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 22 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<22; k++) {
        for (int l=0; l<22; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 22 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
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
        offset_row = nc * 22;
        offset_col = nc * 22;
        for (int k=0; k<22; k++) {
            for (int l=0; l<22; l++) {
                if(J_h[22*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
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
            offset = nc * 22;
            for (int l=0; l<22; l++) {
                for (int k=0; k<22; k++) {
                    if(J_h[22*k + l] != 0.0) {
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
            offset = nc * 22;
            for (int l=0; l<22; l++) {
                for (int k=0; k<22; k++) {
                    if(J_h[22*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
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
            offset = nc * 22;
            for (int l=0; l<22; l++) {
                for (int k=0; k<22; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[22*k + l] != 0.0) {
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
            offset = nc * 22;
            for (int l=0; l<22; l++) {
                for (int k=0; k<22; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[22*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
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
    for (int k=0; k<22; k++) {
        for (int l=0; l<22; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 22*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[22*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 22*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
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
        for (int l=0; l<22; l++) {
            for (int k=0; k<22; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[22*k + l] != 0.0) {
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
        for (int l=0; l<22; l++) {
            for (int k=0; k<22; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[22*k + l] != 0.0) {
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

